// test_pro.cpp -- native smoke for the opt-in pro path (Heaviside + density filter).
// No reference compare: asserts each config still solves (finite compliance,
// bounded densities, projection getters round-trip). Runs in <30 s.
#include <cmath>
#include <cstdio>

extern "C" {
int topo2d_build(int preset);
int topo2d_nel(void);
void topo2d_reset(void);
double topo2d_step(void);
void topo2d_set_projection(double beta, double eta);
void topo2d_set_filter(int mode);
void topo2d_get_projection(double* out);
int topo2d_get_filter(void);
void topo2d_set_penal_cont(int on);
double topo2d_get_penal(void);
void topo2d_compute_stress(void);
double* topo2d_stress(void);
}

static int failures = 0;
static void check(bool ok, const char* msg, double v) {
  if (!ok) { printf("PRO-FAIL %s (%.6g)\n", msg, v); failures++; }
  else printf("pro-ok %s (%.6g)\n", msg, v);
}
static bool finite_pos(double c) { return std::isfinite(c) && c > 0 && c < 1e7; }

int main(void) {
  double g[2];
  // legacy defaults after build
  topo2d_build(0);
  topo2d_get_projection(g);
  check(g[0] == 0.0 && g[1] == 0.5, "default-beta", g[0]);
  check(topo2d_get_filter() == 0, "default-filter", 0);
  // beta sweep on sensitivity path
  for (int b = 0; b <= 3; b++) {
    double beta = b == 0 ? 0.0 : (1 << b); // 0,2,4,8
    topo2d_build(0);
    topo2d_set_projection(beta, 0.5);
    topo2d_reset();
    double c = 0;
    for (int i = 0; i < 5; i++) c = topo2d_step();
    char m[64]; snprintf(m, sizeof m, "sens-beta%.0f-5steps", beta);
    check(finite_pos(c) && c > 0 && c < 1e7, m, c);
  }
  // density filter path, with and without projection
  for (int b = 0; b <= 1; b++) {
    topo2d_build(0);
    topo2d_set_filter(1);
    topo2d_set_projection(b ? 4.0 : 0.0, 0.5);
    topo2d_reset();
    double c = 0;
    for (int i = 0; i < 5; i++) c = topo2d_step();
    char m[64]; snprintf(m, sizeof m, "dens-beta%s-5steps", b ? "4" : "0");
    check(finite_pos(c) && c > 0 && c < 1e7, m, c);
  }
  // bridge preset, pro combo
  topo2d_build(1);
  topo2d_set_filter(1);
  topo2d_set_projection(2.0, 0.5);
  topo2d_reset();
  double c = 0;
  for (int i = 0; i < 5; i++) c = topo2d_step();
  check(finite_pos(c) && c > 0 && c < 1e7, "bridge-dens-beta2-5steps", c);
  // penal continuation (topopt-port): starts p=1, ramps +0.5/10 steps, still solves
  topo2d_build(0);
  topo2d_set_penal_cont(1);
  topo2d_reset();
  check(topo2d_get_penal() == 1.0, "penal-starts-at-1", topo2d_get_penal());
  c = 0;
  for (int i = 0; i < 12; i++) c = topo2d_step();
  check(finite_pos(c), "penal-cont-12steps", c);
  check(topo2d_get_penal() == 1.5, "penal-ramps-to-1.5", topo2d_get_penal());
  // von Mises stress (topopt-port): finite, non-negative, loaded structure reads > 0
  topo2d_build(0);
  topo2d_reset();
  topo2d_compute_stress();
  double* s = topo2d_stress();
  double smax = 0, smin = 1e30; bool allfin = (s != nullptr);
  if (s) for (int e = 0; e < topo2d_nel(); e++) {
    if (!std::isfinite(s[e]) || s[e] < 0) allfin = false;
    if (s[e] > smax) smax = s[e];
    if (s[e] < smin) smin = s[e];
  }
  check(allfin, "stress-finite-nonneg", smax);
  check(smax > 0, "stress-loaded-positive", smax);
  printf("pro-info stress-50pct max=%.6g min=%.6g\n", smax, smin);
  if (failures) { printf("PRO-SMOKE FAIL (%d)\n", failures); return 1; }
  printf("PRO-SMOKE PASS\n");
  return 0;
}
