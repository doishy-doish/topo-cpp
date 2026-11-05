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

int main() {
  build_ke();
  auto t0 = std::chrono::steady_clock::now();
  Sim3 s;
  std::fill(s.x.begin(), s.x.end(), 1.0);
  s.assemble();
  auto [u0, cg0] = s.solve();
  double tip = fabs(u0[3*s.loadNode+1]);
  double I = (double)NELZ * pow((double)NELY, 3) / 12.0;
  double theory = pow((double)NELX, 3) / (3.0 * E0 * I);
  std::fill(s.x.begin(), s.x.end(), VOLFRAC);
  printf("{\"contract\":\"topo-bench-v2\",\"preset\":\"cantilever3d\",\"nelx\":%d,\"nely\":%d,\"nelz\":%d,"
         "\"volfrac\":%.1f,\"penal\":%.1f,\"rmin\":%.1f,\"iters\":%d,",
         NELX, NELY, NELZ, VOLFRAC, PENAL, RMIN, BENCH_ITERS);
  printf("\"tip_full_density\":%.4f,\"beam_theory\":%.4f,\"cg_verify_iters\":%d,", tip, theory, cg0);
  printf("\"compliance_trajectory\":[");
  double ch = 0;
  for (int i = 0; i < BENCH_ITERS; i++) {
    double c = s.step(ch);
    printf("%s%.4f", (i==0?"":","), c);
  }
  auto t1 = std::chrono::steady_clock::now();
  long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
  int nth = 1;
#ifdef _OPENMP
  nth = omp_get_max_threads();
#endif
  printf("],\"cpp_time_ms\":%ld,\"threads\":%d,\"final_change\":%.4f}\n", ms, nth, ch);
  return 0;
}
