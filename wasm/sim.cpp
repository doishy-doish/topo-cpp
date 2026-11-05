// sim.cpp -- unified SIMP core (2D corpus + 3D cantilever) with a C API for wasm.
// Same numerics as main.cpp / bridge.cpp / mbb.cpp / lbracket.cpp / topo3d.cpp.
// Native-tested by wasm/test_sim.cpp, browser-built with emcc (SINGLE_FILE).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <array>
#include <set>
#include <vector>
#include "../memtrash.hpp"

namespace core2d {
static const double VOLFRAC = 0.5;
static double PENAL = 3.0;
static const double EMIN = 1e-9, E0 = 1.0, NU = 0.3;
// ---- higher-language port (zfergus/topopt, MIT (c) 2019 Zachary Ferguson,
// re-implemented for zero-dep C++): penal continuation + von Mises stress.
// Defaults (PENAL_CONT=0) reproduce the legacy p=3 path bit-for-bit.
static int PENAL_CONT = 0, penalSteps = 0;

// ---- pro upgrades (rodin-port: Heaviside projection + density-filter mode) ----
// Defaults (PROJ_BETA<=0, FILTER_MODE=0) reproduce the legacy sensitivity path
// bit-for-bit: heav() early-returns x, density branch never taken.
static double PROJ_BETA = 0.0, PROJ_ETA = 0.5;
static int FILTER_MODE = 0; // 0 = sensitivity filter (legacy), 1 = density filter
static inline double heav(double x) {
  if (PROJ_BETA <= 0.0) return x;
  double b = PROJ_BETA, e = PROJ_ETA;
  return (tanh(b*e) + tanh(b*(x-e))) / (tanh(b*e) + tanh(b*(1.0-e)));
}
static inline double dheav(double x) {
  if (PROJ_BETA <= 0.0) return 1.0;
  double b = PROJ_BETA, e = PROJ_ETA, c = cosh(b*(x-e));
  return b / (c*c*(tanh(b*e) + tanh(b*(1.0-e))));
}

static double KE[8][8];
static double KEF[64];
static double EBm[3][8]; // E(NU) @ B(1): plane-stress constitutive times strain-displacement

static void build_ke() {
  double A11[4][4] = {{12,3,-6,-3},{3,12,3,0},{-6,3,12,-3},{-3,0,-3,12}};
  double A12[4][4] = {{-6,-3,0,3},{-3,-6,-3,-6},{0,-3,-6,3},{3,-6,3,-6}};
  double B11[4][4] = {{-4,3,-2,9},{3,-4,-9,4},{-2,-9,-4,-3},{9,4,-3,-4}};
  double B12[4][4] = {{2,-3,4,-9},{-3,2,9,-2},{4,9,2,3},{-9,-2,3,2}};
  double s = 1.0 / ((1 - NU * NU) * 24);
  for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) KE[i][j] = 0;
  auto blk = [&](double M[4][4], double w, int r0, int c0) {
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) KE[r0+i][c0+j] += w * M[i][j];
  };
  blk(A11, s, 0, 0); blk(A12, s, 0, 4); blk(A11, s, 4, 4);
  for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) KE[4+i][j] += s * A12[j][i];
  blk(B11, s * NU, 0, 0); blk(B12, s * NU, 0, 4); blk(B11, s * NU, 4, 4);
  for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) KE[4+i][j] += s * NU * B12[j][i];
  for (int i = 0; i < 8; i++) for (int j = 0; j < 8; j++) KEF[i*8+j] = KE[i][j];
  // EB = E(NU) @ B(side=1), after zfergus/topopt von_mises_stress.py (MIT):
  // unit-quad center B rows + plane-stress E, so per-element stress is EB*ue.
  {
    double Bm[3][8] = {{0.5,0,-0.5,0,-0.5,0,0.5,0},
                       {0,0.5,0,0.5,0,-0.5,0,-0.5},
                       {0.5,0.5,0.5,-0.5,-0.5,-0.5,-0.5,0.5}};
    double inv = 1.0 / (1 - NU*NU), q = (1 - NU) / 2.0;
    for (int j = 0; j < 8; j++) {
      EBm[0][j] = (Bm[0][j] + NU*Bm[1][j]) * inv;
      EBm[1][j] = (NU*Bm[0][j] + Bm[1][j]) * inv;
      EBm[2][j] = (q*Bm[2][j]) * inv;
    }
  }
}

struct Sim2 {
  int nx, ny;
  double rmin;
  std::vector<int> holeEls;
  int nel, ndof, nn;
  std::vector<std::array<int,8>> edof;
  std::vector<int> rowPtr, cols, diagIdx, scatP;
  std::vector<unsigned char> scatK;
  std::vector<double> vals, F, x, Hs;
  std::vector<double> xphys, xt_; // physical densities (projection/density-filter scratch)
  std::vector<double> sbuf; // per-element von Mises (display-only, never feeds OC)
  std::vector<char> isFixed;
  std::vector<int> freeDofs;
  std::vector<std::vector<std::pair<int,double>>> Hn;
  int loadNode = 0;

  int nid(int ix, int iy) { return ix * (ny + 1) + iy; }

  Sim2(int nx_, int ny_, double rmin_, int preset) : nx(nx_), ny(ny_), rmin(rmin_) {
    nel = nx * ny; nn = (nx + 1) * (ny + 1); ndof = 2 * nn;
    for (int ex = 0; ex < nx; ex++) for (int ey = 0; ey < ny; ey++) {
      int n1 = nid(ex,ey), n2 = nid(ex+1,ey), n3 = nid(ex+1,ey+1), n4 = nid(ex,ey+1);
      edof.push_back({2*n1,2*n1+1,2*n2,2*n2+1,2*n3,2*n3+1,2*n4,2*n4+1});
    }
    std::vector<std::set<int>> rows(ndof);
    for (auto &e : edof) for (int a : e) for (int b : e) rows[a].insert(b);
    rowPtr.assign(ndof + 1, 0); diagIdx.assign(ndof, 0);
    std::vector<int> c; int nnz = 0;
    for (int i = 0; i < ndof; i++) {
      rowPtr[i] = nnz;
      int k = 0;
      for (int col : rows[i]) { if (col == i) diagIdx[i] = nnz + k; c.push_back(col); k++; }
      nnz += (int)rows[i].size();
    }
    rowPtr[ndof] = nnz;
    cols = c; vals.assign(nnz, 0.0);
    scatP.assign(nel * 64, -1); scatK.assign(nel * 64, 0);
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      for (int a = 0; a < 8; a++) for (int b = 0; b < 8; b++) {
        int r = ed[a], cc = ed[b], lo = rowPtr[r], hi = rowPtr[r+1]-1, pos = -1;
        while (lo <= hi) { int m = (lo+hi)>>1; if (cols[m]==cc){pos=m;break;} else if (cols[m]<cc) lo=m+1; else hi=m-1; }
        scatP[e*64+a*8+b] = pos; scatK[e*64+a*8+b] = (unsigned char)(a*8+b);
      }
    }
    isFixed.assign(ndof, 0);
    F.assign(ndof, 0.0);
    if (preset == 0) { // cantilever: left edge fixed, right-middle downward load
      for (int iy = 0; iy <= ny; iy++) { isFixed[2*nid(0,iy)] = 1; isFixed[2*nid(0,iy)+1] = 1; }
      loadNode = nid(nx, ny/2);
    } else if (preset == 1) { // bridge: both bottom corners fixed, center-bottom load
      isFixed[2*nid(0,0)] = 1; isFixed[2*nid(0,0)+1] = 1;
      isFixed[2*nid(nx,0)] = 1; isFixed[2*nid(nx,0)+1] = 1;
      loadNode = nid(nx/2, 0);
    } else if (preset == 2) { // mbb: pin left-bottom, roller right-bottom, top-center load
      isFixed[2*nid(0,0)] = 1; isFixed[2*nid(0,0)+1] = 1;
      isFixed[2*nid(nx,0)+1] = 1;
      loadNode = nid(nx/2, ny);
    } else { // lbracket: top edge of stem fixed, tip of bottom flange loaded
      for (int ix = 0; ix <= 30; ix++) { isFixed[2*nid(ix,ny)] = 1; isFixed[2*nid(ix,ny)+1] = 1; }
      loadNode = nid(nx, 15);
      for (int ex = 0; ex < nx; ex++) for (int ey = 0; ey < ny; ey++)
        if (ex >= 30 && ey >= 30) holeEls.push_back(ex*ny+ey);
    }
    F[2*loadNode+1] = -1.0;
    for (int i = 0; i < ndof; i++) if (!isFixed[i]) freeDofs.push_back(i);
    // filter
    Hs.assign(nel, 0.0); Hn.resize(nel);
    for (int e1 = 0; e1 < nel; e1++) {
      int i1 = e1 / ny, j1 = e1 % ny, R = (int)ceil(rmin);
      for (int i2 = std::max(0,i1-R); i2 <= std::min(nx-1,i1+R); i2++)
        for (int j2 = std::max(0,j1-R); j2 <= std::min(ny-1,j1+R); j2++) {
          double fac = rmin - hypot((double)(i1-i2), (double)(j1-j2));
          if (fac > 0) { int e2 = i2*ny+j2; Hn[e1].push_back({e2,fac}); Hs[e1] += fac; }
        }
    }
    x.assign(nel, VOLFRAC);
    for (int e : holeEls) x[e] = 0.001;
  }

  void compute_phys(void) {
    if ((int)xphys.size() != nel) { xphys.assign(nel, 0.0); xt_.assign(nel, 0.0); }
    if (FILTER_MODE == 1) {
      for (int e1 = 0; e1 < nel; e1++) {
        double s = 0;
        for (auto &pr : Hn[e1]) s += pr.second * x[pr.first];
        xt_[e1] = s / Hs[e1];
        xphys[e1] = heav(xt_[e1]);
      }
    } else {
      for (int e = 0; e < nel; e++) xphys[e] = heav(x[e]);
    }
  }

  void assemble() {
    compute_phys();
    std::fill(vals.begin(), vals.end(), 0.0);
    for (int e = 0; e < nel; e++) {
      double s = EMIN + pow(xphys[e], PENAL) * (E0 - EMIN);
      int base = e * 64;
      for (int k = 0; k < 64; k++) vals[scatP[base+k]] += s * KEF[scatK[base+k]];
    }
  }

  void matvec(const std::vector<double> &v, std::vector<double> &out) {
    for (int i = 0; i < ndof; i++) {
      double s = 0;
      for (int k = rowPtr[i]; k < rowPtr[i+1]; k++) s += vals[k] * v[cols[k]];
      out[i] = s;
    }
  }

  // returns {displacement, cg_iters}
  std::pair<std::vector<double>,int> solve(double tol = 1e-6, int maxIt = 2000) {
    int n = (int)freeDofs.size();
    std::vector<double> u(ndof, 0.0), full(ndof, 0.0), Ap(ndof, 0.0);
    std::vector<double> r(n), z(n), p(n), D(n);
    for (int i = 0; i < n; i++) r[i] = F[freeDofs[i]];
    for (int i = 0; i < n; i++) { double d = vals[diagIdx[freeDofs[i]]]; D[i] = d != 0 ? 1/d : 1; z[i] = r[i]*D[i]; p[i] = z[i]; }
    double rz = 0; for (int i = 0; i < n; i++) rz += r[i]*z[i];
    double r0 = sqrt(rz); if (r0 == 0) r0 = 1;
    int it = 0;
    for (; it < maxIt; it++) {
      for (int i = 0; i < n; i++) full[freeDofs[i]] = p[i];
      matvec(full, Ap);
      double pAp = 0; for (int i = 0; i < n; i++) pAp += p[i]*Ap[freeDofs[i]];
      if (pAp == 0) break;
      double a = rz / pAp;
      for (int i = 0; i < n; i++) { u[freeDofs[i]] += a*p[i]; r[i] -= a*Ap[freeDofs[i]]; }
      for (int i = 0; i < n; i++) z[i] = r[i]*D[i];
      double rz2 = 0; for (int i = 0; i < n; i++) rz2 += r[i]*z[i];
      if (sqrt(rz2)/r0 < tol) { it++; break; }
      double b = rz2/rz; rz = rz2;
      for (int i = 0; i < n; i++) p[i] = z[i] + b*p[i];
    }
    return {u, it};
  }

  double step(double &change) {
    if (PENAL_CONT) PENAL = std::min(3.0, 1.0 + 0.5 * (penalSteps / 10));
    penalSteps++;
    assemble(); // refreshes xphys
    auto [u, it] = solve();
    (void)it;
    std::vector<double> dc(nel);
    double c = 0;
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      double ue[8]; for (int a = 0; a < 8; a++) ue[a] = u[ed[a]];
      double ee = 0;
      for (int a = 0; a < 8; a++) { double s = 0; for (int b = 0; b < 8; b++) s += KE[a][b]*ue[b]; ee += ue[a]*s; }
      double phys = (FILTER_MODE == 1) ? xt_[e] : x[e];
      c += pow(xphys[e], PENAL) * (E0-EMIN) * ee;
      dc[e] = -PENAL * pow(xphys[e], PENAL-1) * (E0-EMIN) * ee * dheav(phys);
    }
    std::vector<double> dcf(nel, 0.0);
    if (FILTER_MODE == 1) {
      // chain rule through the density filter: scatter dc via the same weights
      for (int e1 = 0; e1 < nel; e1++)
        for (auto &pr : Hn[e1]) dcf[pr.first] += pr.second * dc[e1] / Hs[e1];
    } else {
      for (int e1 = 0; e1 < nel; e1++) {
        double s = 0;
        for (auto &pr : Hn[e1]) s += pr.second * x[pr.first] * dc[pr.first];
        dcf[e1] = s / (Hs[e1] * std::max(1e-3, x[e1]));
      }
    }
    const double move = 0.2;
    double l1 = 0, l2 = 1e9;
    std::vector<double> xn(nel);
    while (l2 - l1 > 1e-4) {
      double lm = 0.5*(l1+l2), vol = 0;
      for (int e = 0; e < nel; e++) {
        double b = x[e] * sqrt(std::max(0.0, -dcf[e]/lm));
        double lo = std::max(0.001, x[e]-move), hi = std::min(1.0, x[e]+move);
        xn[e] = std::min(hi, std::max(lo, std::min(1.0, b)));
        vol += xn[e];
      }
      if (vol - VOLFRAC*nel > 0) l1 = lm; else l2 = lm;
    }
    change = 0;
    for (int e = 0; e < nel; e++) { change = std::max(change, fabs(xn[e]-x[e])); x[e] = xn[e]; }
    for (int e : holeEls) x[e] = 0.001;
    return c;
  }

  void compute_stress(void) { // per-element von Mises, display-only (never feeds the OC update)
    if ((int)sbuf.size() != nel) sbuf.assign(nel, 0.0);
    assemble();
    auto r = solve();
    const std::vector<double> &u = r.first;
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      double ue[8]; for (int a = 0; a < 8; a++) ue[a] = u[ed[a]];
      double s0 = 0, s1 = 0, s2 = 0;
      for (int b = 0; b < 8; b++) { s0 += EBm[0][b]*ue[b]; s1 += EBm[1][b]*ue[b]; s2 += EBm[2][b]*ue[b]; }
      double rho = EMIN + pow(xphys[e], PENAL) * (E0 - EMIN);
      s0 *= rho; s1 *= rho; s2 *= rho;
      sbuf[e] = sqrt(s0*s0 - s0*s1 + s1*s1 + 3*s2*s2);
    }
  }
};

} // namespace core2d

// topo3d.cpp — 3D hex SIMP strangler port. C++17 + OpenMP, zero other deps.
// Differential-tested against bench/reference3d.json (JS reference), contract topo-bench-v2.
// Chunking: elements assemble in parallel (atomic scatter), rows solve in parallel.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace core3d {
static const int NELX = 40, NELY = 12, NELZ = 6;
static const double VOLFRAC = 0.5, PENAL = 3.0, RMIN = 1.5;
static const double EMIN = 1e-9, E0 = 1.0, NU = 0.3;
static const int BENCH_ITERS = 10;

static double KE[24][24];
static double KEF[576];

static void build_ke() {
  double c = 1.0 / ((1 + NU) * (1 - 2 * NU));
  double D[6][6] = {{0}};
  D[0][0] = D[1][1] = D[2][2] = c * (1 - NU);
  D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = c * NU;
  D[3][3] = D[4][4] = D[5][5] = c * (1 - 2 * NU) / 2;
  for (int i = 0; i < 24; i++) for (int j = 0; j < 24; j++) KE[i][j] = 0;
  double g = 1.0 / sqrt(3.0);
  double xs[2] = {-g, g};
  for (double xi : xs) for (double et : xs) for (double ze : xs) {
    double B[6][24] = {{0}};
    for (int a = 0; a < 8; a++) {
      int sx = (a & 4) ? 1 : -1, sy = (a & 2) ? 1 : -1, sz = (a & 1) ? 1 : -1;
      double dx = 2 * sx * (1 + sy * et) * (1 + sz * ze) / 8;
      double dy = 2 * sy * (1 + sx * xi) * (1 + sz * ze) / 8;
      double dz = 2 * sz * (1 + sx * xi) * (1 + sy * et) / 8;
      int b = a * 3;
      B[0][b] = dx; B[1][b+1] = dy; B[2][b+2] = dz;
      B[3][b] = dy; B[3][b+1] = dx;
      B[4][b+1] = dz; B[4][b+2] = dy;
      B[5][b] = dz; B[5][b+2] = dx;
    }
    double M[6][24] = {{0}};
    for (int k = 0; k < 6; k++) for (int j = 0; j < 24; j++) {
      double s = 0; for (int l = 0; l < 6; l++) s += D[k][l] * B[l][j]; M[k][j] = s;
    }
    for (int i = 0; i < 24; i++) for (int j = 0; j < 24; j++) {
      double s = 0; for (int k = 0; k < 6; k++) s += B[k][i] * M[k][j];
      KE[i][j] += s * 0.125;
    }
  }
  for (int i = 0; i < 24; i++) for (int j = 0; j < 24; j++) KEF[i*24+j] = KE[i][j];
}

struct Sim3 {
  int nel, ndof;
  std::vector<std::array<int,24>> edof;
  std::vector<int> rowPtr, cols, diagIdx, scatP;
  std::vector<unsigned short> scatK;
  std::vector<double> vals, F, x, Hs;
  std::vector<char> isFixed;
  std::vector<int> freeDofs;
  std::vector<std::vector<std::pair<int,double>>> Hn;
  int loadNode = 0;
  int nid(int ix, int iy, int iz) { return (ix * (NELY + 1) + iy) * (NELZ + 1) + iz; }

  Sim3() {
    nel = NELX * NELY * NELZ;
    int nn = (NELX + 1) * (NELY + 1) * (NELZ + 1);
    ndof = 3 * nn;
    for (int ex = 0; ex < NELX; ex++) for (int ey = 0; ey < NELY; ey++) for (int ez = 0; ez < NELZ; ez++) {
      std::array<int,24> e; int k = 0;
      for (int dx = 0; dx < 2; dx++) for (int dy = 0; dy < 2; dy++) for (int dz = 0; dz < 2; dz++) {
        int n = nid(ex+dx, ey+dy, ez+dz);
        e[k++] = 3*n; e[k++] = 3*n+1; e[k++] = 3*n+2;
      }
      edof.push_back(e);
    }
    std::vector<std::set<int>> rows(ndof);
    for (auto &e : edof) for (int a : e) for (int b : e) rows[a].insert(b);
    rowPtr.assign(ndof + 1, 0); diagIdx.assign(ndof, 0);
    std::vector<int> c; int nnz = 0;
    for (int i = 0; i < ndof; i++) {
      rowPtr[i] = nnz; int k = 0;
      for (int col : rows[i]) { if (col == i) diagIdx[i] = nnz + k; c.push_back(col); k++; }
      nnz += (int)rows[i].size();
    }
    rowPtr[ndof] = nnz;
    cols = c; vals.assign(nnz, 0.0);
    scatP.assign(nel * 576, -1); scatK.assign(nel * 576, 0);
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      for (int a = 0; a < 24; a++) for (int b = 0; b < 24; b++) {
        int r = ed[a], cc = ed[b], lo = rowPtr[r], hi = rowPtr[r+1]-1, pos = -1;
        while (lo <= hi) { int mm = (lo+hi)>>1; if (cols[mm]==cc){pos=mm;break;} else if (cols[mm]<cc) lo=mm+1; else hi=mm-1; }
        scatP[e*576+a*24+b] = pos; scatK[e*576+a*24+b] = (unsigned short)(a*24+b);
      }
    }
    isFixed.assign(ndof, 0);
    for (int iy = 0; iy <= NELY; iy++) for (int iz = 0; iz <= NELZ; iz++) {
      int n = nid(0, iy, iz);
      isFixed[3*n] = isFixed[3*n+1] = isFixed[3*n+2] = 1;
    }
    F.assign(ndof, 0.0);
    int my = NELY / 2;
    for (int iz = 0; iz <= NELZ; iz++) F[3 * nid(NELX, my, iz) + 1] = -1.0 / (NELZ + 1);
    loadNode = nid(NELX, my, (NELZ + 1) / 2);
    for (int i = 0; i < ndof; i++) if (!isFixed[i]) freeDofs.push_back(i);
    Hs.assign(nel, 0.0); Hn.resize(nel);
    int R = (int)ceil(RMIN);
    for (int e1 = 0; e1 < nel; e1++) {
      int i1 = e1 / (NELY * NELZ), j1 = (e1 / NELZ) % NELY, k1 = e1 % NELZ;
      for (int i2 = std::max(0,i1-R); i2 <= std::min(NELX-1,i1+R); i2++)
        for (int j2 = std::max(0,j1-R); j2 <= std::min(NELY-1,j1+R); j2++)
          for (int k2 = std::max(0,k1-R); k2 <= std::min(NELZ-1,k1+R); k2++) {
            double fac = RMIN - sqrt((double)((i1-i2)*(i1-i2) + (j1-j2)*(j1-j2) + (k1-k2)*(k1-k2)));
            if (fac > 0) { int e2 = (i2*NELY+j2)*NELZ+k2; Hn[e1].push_back({e2,fac}); Hs[e1] += fac; }
          }
    }
    x.assign(nel, VOLFRAC);
  }

  void assemble() {
    std::fill(vals.begin(), vals.end(), 0.0);
#pragma omp parallel for schedule(static)
    for (int e = 0; e < nel; e++) {
      double s = EMIN + pow(x[e], PENAL) * (E0 - EMIN);
      int base = e * 576;
      for (int k = 0; k < 576; k++) {
#pragma omp atomic update
        vals[scatP[base+k]] += s * KEF[scatK[base+k]];
      }
    }
  }

  void matvec(const std::vector<double> &v, std::vector<double> &out) {
#pragma omp parallel for schedule(static)
    for (int i = 0; i < ndof; i++) {
      double s = 0;
      for (int k = rowPtr[i]; k < rowPtr[i+1]; k++) s += vals[k] * v[cols[k]];
      out[i] = s;
    }
  }

  std::pair<std::vector<double>,int> solve(double tol = 1e-5, int maxIt = 2000) {
    int n = (int)freeDofs.size();
    std::vector<double> u(ndof, 0.0), full(ndof, 0.0), Ap(ndof, 0.0);
    std::vector<double> r(n), z(n), p(n), D(n);
    for (int i = 0; i < n; i++) r[i] = F[freeDofs[i]];
    for (int i = 0; i < n; i++) { double d = vals[diagIdx[freeDofs[i]]]; D[i] = d != 0 ? 1/d : 1; z[i] = r[i]*D[i]; p[i] = z[i]; }
    double rz = 0; for (int i = 0; i < n; i++) rz += r[i]*z[i];
    double r0 = sqrt(rz); if (r0 == 0) r0 = 1;
    int it = 0;
    for (; it < maxIt; it++) {
      for (int i = 0; i < n; i++) full[freeDofs[i]] = p[i];
      matvec(full, Ap);
      double pAp = 0; for (int i = 0; i < n; i++) pAp += p[i]*Ap[freeDofs[i]];
      if (pAp == 0) break;
      double a = rz / pAp;
      for (int i = 0; i < n; i++) { u[freeDofs[i]] += a*p[i]; r[i] -= a*Ap[freeDofs[i]]; }
      for (int i = 0; i < n; i++) z[i] = r[i]*D[i];
      double rz2 = 0; for (int i = 0; i < n; i++) rz2 += r[i]*z[i];
      if (sqrt(rz2)/r0 < tol) { it++; break; }
      double b = rz2/rz; rz = rz2;
      for (int i = 0; i < n; i++) p[i] = z[i] + b*p[i];
    }
    return {u, it};
  }

  double step(double &change) {
    assemble();
    auto [u, it] = solve();
    (void)it;
    std::vector<double> dc(nel);
    double c = 0;
#pragma omp parallel for schedule(static) reduction(+:c)
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      double ue[24]; for (int a = 0; a < 24; a++) ue[a] = u[ed[a]];
      double ee = 0;
      for (int a = 0; a < 24; a++) { double s = 0; for (int b = 0; b < 24; b++) s += KE[a][b]*ue[b]; ee += ue[a]*s; }
      c += pow(x[e], PENAL) * (E0-EMIN) * ee;
      dc[e] = -PENAL * pow(x[e], PENAL-1) * (E0-EMIN) * ee;
    }
    std::vector<double> dcf(nel, 0.0);
#pragma omp parallel for schedule(static)
    for (int e1 = 0; e1 < nel; e1++) {
      double s = 0;
      for (auto &pr : Hn[e1]) s += pr.second * x[pr.first] * dc[pr.first];
      dcf[e1] = s / (Hs[e1] * std::max(1e-3, x[e1]));
    }
    const double move = 0.2;
    double l1 = 0, l2 = 1e9;
    std::vector<double> xn(nel);
    while (l2 - l1 > 1e-4) {
      double lm = 0.5*(l1+l2), vol = 0;
      for (int e = 0; e < nel; e++) {
        double b = x[e] * sqrt(std::max(0.0, -dcf[e]/lm));
        double lo = std::max(0.001, x[e]-move), hi = std::min(1.0, x[e]+move);
        xn[e] = std::min(hi, std::max(lo, std::min(1.0, b)));
        vol += xn[e];
      }
      if (vol - VOLFRAC*nel > 0) l1 = lm; else l2 = lm;
    }
    change = 0;
    for (int e = 0; e < nel; e++) { change = std::max(change, fabs(xn[e]-x[e])); x[e] = xn[e]; }
    return c;
  }
};

} // namespace core3d

extern "C" {
static core2d::Sim2* g2 = nullptr;
static core3d::Sim3* g3 = nullptr;
static bool ke_done = false;
static void ensure_ke(void) {
  if (!ke_done) { core2d::build_ke(); core3d::build_ke(); ke_done = true; }
}
static void dims2d(int preset, int* nx, int* ny, double* rmin) {
  if (preset == 3) { *nx = 60; *ny = 60; *rmin = 2.0; }
  else { *nx = 120; *ny = 20; *rmin = 2.4; }
}
int topo2d_build(int preset) {
  if (preset < 0 || preset > 3) return 0;
  ensure_ke();
  // pro opts reset to legacy defaults so CI/native tests stay deterministic
  core2d::PROJ_BETA = 0.0; core2d::PROJ_ETA = 0.5; core2d::FILTER_MODE = 0;
  core2d::PENAL = 3.0; core2d::PENAL_CONT = 0; core2d::penalSteps = 0;
  int nx, ny; double rmin; dims2d(preset, &nx, &ny, &rmin);
  // memtrash: dead sim leaves the registry as it dies; new sim tracked permanent
  // (audited, never swept). Vectors never reallocate post-construction, pointers stable.
  if (g2) {
    MT_UNTRACK(g2->x.data()); MT_UNTRACK(g2->vals.data()); MT_UNTRACK(g2->cols.data());
    delete g2; g2 = nullptr;
  }
  g2 = new core2d::Sim2(nx, ny, rmin, preset);
  MT_TRACK_PERM(g2, sizeof(*g2), "sim2-obj");
  MT_TRACK_PERM(g2->x.data(), g2->x.size() * sizeof(double), "sim2-x");
  MT_TRACK_PERM(g2->vals.data(), g2->vals.size() * sizeof(double), "sim2-vals");
  MT_TRACK_PERM(g2->cols.data(), g2->cols.size() * sizeof(int), "sim2-cols");
  return g2->nel;
}
int topo2d_nel(void) { return g2 ? g2->nel : 0; }
double* topo2d_x(void) { return g2 ? g2->x.data() : nullptr; }
void topo2d_ones(void) {
  if (!g2) return;
  std::fill(g2->x.begin(), g2->x.end(), 1.0);
  for (int e : g2->holeEls) g2->x[e] = 0.001;
}
void topo2d_reset(void) {
  if (!g2) return;
  std::fill(g2->x.begin(), g2->x.end(), core2d::VOLFRAC);
  for (int e : g2->holeEls) g2->x[e] = 0.001;
}
void topo2d_verify(double* out) {
  if (!g2 || !out) return;
  g2->assemble();
  auto r = g2->solve();
  out[0] = fabs(r.first[2*g2->loadNode+1]);
  out[1] = (double)r.second;
}
double topo2d_step(void) { double ch = 0; memtrash::global().tick(); return g2 ? g2->step(ch) : 0.0; }
// pro controls (opt-in; defaults = legacy path, see topo2d_build reset)
void topo2d_set_projection(double beta, double eta) { core2d::PROJ_BETA = beta; core2d::PROJ_ETA = eta; }
void topo2d_set_filter(int mode) { core2d::FILTER_MODE = mode ? 1 : 0; }
void topo2d_get_projection(double* out) { if (out) { out[0] = core2d::PROJ_BETA; out[1] = core2d::PROJ_ETA; } }
int topo2d_get_filter(void) { return core2d::FILTER_MODE; }
// higher-language port controls (opt-in; defaults = legacy path)
void topo2d_set_penal_cont(int on) {
  if (on) { core2d::PENAL_CONT = 1; core2d::PENAL = 1.0; core2d::penalSteps = 0; }
  else { core2d::PENAL_CONT = 0; core2d::PENAL = 3.0; }
}
double topo2d_get_penal(void) { return core2d::PENAL; }
void topo2d_compute_stress(void) { if (g2) g2->compute_stress(); }
double* topo2d_stress(void) { return (g2 && !g2->sbuf.empty()) ? g2->sbuf.data() : nullptr; }
int topo3d_build(void) {
  ensure_ke();
  if (g3) {
    MT_UNTRACK(g3->x.data()); MT_UNTRACK(g3->vals.data()); MT_UNTRACK(g3->cols.data());
    delete g3; g3 = nullptr;
  }
  g3 = new core3d::Sim3();
  MT_TRACK_PERM(g3, sizeof(*g3), "sim3-obj");
  MT_TRACK_PERM(g3->x.data(), g3->x.size() * sizeof(double), "sim3-x");
  MT_TRACK_PERM(g3->vals.data(), g3->vals.size() * sizeof(double), "sim3-vals");
  MT_TRACK_PERM(g3->cols.data(), g3->cols.size() * sizeof(int), "sim3-cols");
  return g3->nel;
}
int topo3d_nel(void) { return g3 ? g3->nel : 0; }
double* topo3d_x(void) { return g3 ? g3->x.data() : nullptr; }
void topo3d_ones(void) { if (g3) std::fill(g3->x.begin(), g3->x.end(), 1.0); }
void topo3d_reset(void) { if (g3) std::fill(g3->x.begin(), g3->x.end(), core3d::VOLFRAC); }
void topo3d_verify(double* out) {
  if (!g3 || !out) return;
  g3->assemble();
  auto r = g3->solve();
  out[0] = fabs(r.first[3*g3->loadNode+1]);
  out[1] = (double)r.second;
}
double topo3d_step(void) { double ch = 0; memtrash::global().tick(); return g3 ? g3->step(ch) : 0.0; }
int topo_mt_live(void) { return (int)memtrash::global().live(); }
int topo_mt_bytes(void) { return (int)memtrash::global().bytes(); }
int topo_mt_peak(void) { return (int)memtrash::global().peak(); }
void topo_mt_tick(void) { memtrash::global().tick(); }
int topo_mt_sweep(int max_idle) { return memtrash::global().sweep((uint64_t)max_idle).released; }
int topo_mt_log(char* out, int cap) {
  if (!out || cap <= 0) return 0;
  return memtrash::global().log_json(out, cap);
}
} // extern "C"
