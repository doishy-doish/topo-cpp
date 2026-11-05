#include <cmath>
#include <cstdio>
#include <cstring>

extern "C" {
int topo2d_build(int preset);
int topo2d_nel(void);
void topo2d_ones(void);
void topo2d_reset(void);
void topo2d_verify(double* out);
double topo2d_step(void);
int topo3d_build(void);
void topo3d_ones(void);
void topo3d_reset(void);
void topo3d_verify(double* out);
double topo3d_step(void);
}

// usage:
//   test_sim 2d <preset 0..3>   -> full candidate JSON for that 2D preset
//   test_sim 3d                 -> full candidate JSON for 3D cantilever
int main(int argc, char** argv) {
  if (argc >= 2 && strcmp(argv[1], "3d") == 0) {
    topo3d_build();
    topo3d_ones();
    double v[2]; topo3d_verify(v);
    topo3d_reset();
    printf("{\"contract\":\"topo-bench-v2\",\"preset\":\"cantilever3d\",\"nelx\":40,\"nely\":12,\"nelz\":6,"
           "\"volfrac\":0.5,\"penal\":3.0,\"rmin\":1.5,\"iters\":10,");
    printf("\"tip_full_density\":%.4f,\"cg_verify_iters\":%d,\"compliance_trajectory\":[", v[0], (int)v[1]);
    for (int i = 0; i < 10; i++) printf("%s%.4f", i ? "," : "", topo3d_step());
    printf("]}\n");
    return 0;
  }
  if (argc < 3 || strcmp(argv[1], "2d") != 0) { printf("usage: test_sim 2d <0..3> | test_sim 3d\n"); return 2; }
  int p = atoi(argv[2]);
  const char* contracts[4] = {"topo-bench-v1", "topo-bench-bridge-v1", "topo-bench-mbb-v1", "topo-bench-lbracket-v1"};
  const char* names[4] = {"cantilever", "bridge", "mbb", "lbracket"};
  int dims[4][3] = {{120, 20, 0}, {120, 20, 0}, {120, 20, 0}, {60, 60, 0}}; // nx, ny, (rmin scaled x10: 24 or 20)
  double rmins[4] = {2.4, 2.4, 2.4, 2.0};
  int iters[4] = {40, 20, 20, 15};
  int sampled[4] = {1, 0, 0, 0}; // v1 prints iters 1,5,10,... like main.cpp
  if (p < 0 || p > 3) return 2;
  topo2d_build(p);
  topo2d_ones();
  double v[2]; topo2d_verify(v);
  topo2d_reset();
  printf("{\"contract\":\"%s\",\"preset\":\"%s\",\"nelx\":%d,\"nely\":%d,"
         "\"volfrac\":0.5,\"penal\":3.0,\"rmin\":%.1f,\"iters\":%d,",
         contracts[p], names[p], dims[p][0], dims[p][1], rmins[p], iters[p]);
  printf("\"tip_full_density\":%.4f,\"cg_verify_iters\":%d,\"compliance_trajectory\":[",
         v[0], (int)v[1]);
  int n = 0;
  for (int i = 0; i < iters[p]; i++) {
    double c = topo2d_step();
    if (sampled[p] && !(i == 0 || (i + 1) % 5 == 0)) continue;
    printf("%s%.4f", n++ ? "," : "", c);
  }
  printf("]}\n");
  (void)dims;
  return 0;
}
