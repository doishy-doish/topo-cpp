// lbracket.cpp — L-bracket preset port (60x60 with top-right void). Same SIMP core as main.cpp, own bench contract topo-bench-lbracket-v1.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <set>
#include <vector>

static const int NELX = 60, NELY = 60;
static const double VOLFRAC = 0.5, PENAL = 3.0, RMIN = 2.0;
static const double EMIN = 1e-9, E0 = 1.0, NU = 0.3;
static const int BENCH_ITERS = 15;

static double KE[8][8];
static double KEF[64];

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
}

struct Sim {
  int nel, ndof, nn;
  std::vector<std::array<int,8>> edof;
  std::vector<int> rowPtr, cols, diagIdx, scatP;
  std::vector<unsigned char> scatK;
  std::vector<double> vals, F, x, Hs;
  std::vector<char> isFixed;
  std::vector<int> holeEls;
  std::vector<int> freeDofs;
  std::vector<std::vector<std::pair<int,double>>> Hn;
  int loadNode = 0;

  int nid(int ix, int iy) { return ix * (NELY + 1) + iy; }

  Sim() {
    nel = NELX * NELY; nn = (NELX + 1) * (NELY + 1); ndof = 2 * nn;
    for (int ex = 0; ex < NELX; ex++) for (int ey = 0; ey < NELY; ey++) {
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
    // L-bracket: top edge of stem fixed, tip of bottom flange loaded downward
    isFixed.assign(ndof, 0);
    for (int ix = 0; ix <= 30; ix++) { isFixed[2*nid(ix,NELY)] = 1; isFixed[2*nid(ix,NELY)+1] = 1; }
    F.assign(ndof, 0.0);
    loadNode = nid(NELX, 15);
    F[2*loadNode+1] = -1.0;
    for (int i = 0; i < ndof; i++) if (!isFixed[i]) freeDofs.push_back(i);
    for (int ex = 0; ex < NELX; ex++) for (int ey = 0; ey < NELY; ey++)
      if (ex >= 30 && ey >= 30) holeEls.push_back(ex*NELY+ey);
    for (int e : holeEls) x[e] = 0.001;
    // filter
    Hs.assign(nel, 0.0); Hn.resize(nel);
    for (int e1 = 0; e1 < nel; e1++) {
      int i1 = e1 / NELY, j1 = e1 % NELY, R = (int)ceil(RMIN);
      for (int i2 = std::max(0,i1-R); i2 <= std::min(NELX-1,i1+R); i2++)
        for (int j2 = std::max(0,j1-R); j2 <= std::min(NELY-1,j1+R); j2++) {
          double fac = RMIN - hypot((double)(i1-i2), (double)(j1-j2));
          if (fac > 0) { int e2 = i2*NELY+j2; Hn[e1].push_back({e2,fac}); Hs[e1] += fac; }
        }
    }
    x.assign(nel, VOLFRAC);
  }

  void assemble() {
    std::fill(vals.begin(), vals.end(), 0.0);
    for (int e = 0; e < nel; e++) {
      double s = EMIN + pow(x[e], PENAL) * (E0 - EMIN);
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
    assemble();
    auto [u, it] = solve();
    (void)it;
    std::vector<double> dc(nel);
    double c = 0;
    for (int e = 0; e < nel; e++) {
      auto &ed = edof[e];
      double ue[8]; for (int a = 0; a < 8; a++) ue[a] = u[ed[a]];
      double ee = 0;
      for (int a = 0; a < 8; a++) { double s = 0; for (int b = 0; b < 8; b++) s += KE[a][b]*ue[b]; ee += ue[a]*s; }
      c += pow(x[e], PENAL) * (E0-EMIN) * ee;
      dc[e] = -PENAL * pow(x[e], PENAL-1) * (E0-EMIN) * ee;
    }
    std::vector<double> dcf(nel, 0.0);
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
    for (int e : holeEls) x[e] = 0.001;
    return c;
  }
};

int main() {
  build_ke();
  auto t0 = std::chrono::steady_clock::now();
  Sim s;
  // verify: full density tip at load node
  std::fill(s.x.begin(), s.x.end(), 1.0);
  for (int e : s.holeEls) s.x[e] = 0.001;
  s.assemble();
  auto [u0, cg0] = s.solve();
  double tip = fabs(u0[2*s.loadNode+1]);
  // optimize from uniform 50%
  std::fill(s.x.begin(), s.x.end(), VOLFRAC);
  for (int e : s.holeEls) s.x[e] = 0.001;
  printf("{\"contract\":\"topo-bench-lbracket-v1\",\"preset\":\"lbracket\",\"nelx\":%d,\"nely\":%d,"
         "\"volfrac\":%.1f,\"penal\":%.1f,\"rmin\":%.1f,\"iters\":%d,",
         NELX, NELY, VOLFRAC, PENAL, RMIN, BENCH_ITERS);
  printf("\"tip_full_density\":%.4f,\"cg_verify_iters\":%d,", tip, cg0);
  printf("\"compliance_trajectory\":[");
  double ch = 0;
  for (int i = 0; i < BENCH_ITERS; i++) {
    double c = s.step(ch);
    printf("%s%.4f", (i==0?"":","), c);
  }
  auto t1 = std::chrono::steady_clock::now();
  long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count();
  printf("],\"cpp_time_ms\":%ld,\"final_change\":%.4f}\n", ms, ch);
  return 0;
}
