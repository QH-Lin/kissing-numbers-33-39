/* All-pairs check of the 605,322-point configuration in R^38 written by d38_export.py.
   Every pair x != y is tested for  <x,y> <= |x||y|/2  in exact integer arithmetic (int32 for the int8 rows,
   __int128 for pairs involving an axis row). A repeated point, or two points on the same ray, fails the test.
   The program also checks that each row is nonzero and that coordinates 38..51 repeat 24..37, so all points lie
   in a 38-dimensional subspace of R^52.
   build: cc -O3 -pthread -o d38_full verify/d38_full.c        (Apple silicon: add -mcpu=native)
   usage: ./d38_full OUTDIR/d38_rows_i8.bin OUTDIR/d38_axes_i64.bin */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#if defined(__ARM_FEATURE_DOTPROD)
#include <arm_neon.h>
#endif
#define W 64   /* bytes per int8 row */
#define A 52   /* coordinates per axis row */
typedef __int128 i128;
static int8_t (*V)[W]; static int32_t *NRM; static long N;
static int64_t (*X)[A]; static i128 *XN; static long NX;
static int nth; static long gbad[256], gbad2[256]; static double gmax[256], gmax2[256];

static inline int32_t dot8(const int8_t *a, const int8_t *b) {
#if defined(__ARM_FEATURE_DOTPROD)
  int32x4_t acc = vdupq_n_s32(0);
  for (int k = 0; k < W; k += 16) acc = vdotq_s32(acc, vld1q_s8(a + k), vld1q_s8(b + k));
  return vaddvq_s32(acc);
#else
  int32_t s = 0;
  for (int k = 0; k < W; k++) s += a[k] * b[k];
  return s;
#endif
}

static void *work8(void *arg) {           /* all pairs of int8 rows */
  long t = (long)arg, bad = 0; double mx = -2;
  for (long i = t; i < N; i += nth) {
    int64_t ni = NRM[i];
    for (long j = i + 1; j < N; j++) {
      int64_t s = dot8(V[i], V[j]);
      if (s > 0) {
        int64_t nj = NRM[j];
        if (4 * s * s > ni * nj) bad++;
        double c = (double)s * s / ((double)ni * nj); if (c > mx) mx = c;
      }
    }
  }
  gbad[t] = bad; gmax[t] = mx; return 0;
}

static void *workx(void *arg) {           /* axis rows against all int8 rows and against later axis rows */
  long t = (long)arg, bad = 0; double mx = -2;
  for (long i = t; i < NX; i += nth) {
    i128 ni = XN[i];
    for (long j = 0; j < N; j++) {
      i128 s = 0;
      for (int k = 0; k < A; k++) s += (i128)X[i][k] * V[j][k];
      if (s > 0) {
        i128 nn = ni * NRM[j];             /* < 2^125 by the bound checked in main */
        if (4 * s * s > nn) bad++;
        double c = (double)s * (double)s / (double)nn; if (c > mx) mx = c;
      }
    }
    for (long j = i + 1; j < NX; j++) {    /* equal norms (checked in main): cos <= 1/2  iff  2s <= |x|^2 */
      i128 s = 0;
      for (int k = 0; k < A; k++) s += (i128)X[i][k] * X[j][k];
      if (s > 0) {
        if (2 * s > ni) bad++;
        double c = (double)s / (double)ni; c *= c; if (c > mx) mx = c;
      }
    }
  }
  gbad2[t] = bad; gmax2[t] = mx; return 0;
}

static void *slurp(const char *fn, long *bytes) {
  FILE *f = fopen(fn, "rb"); if (!f) { perror(fn); exit(1); }
  fseek(f, 0, SEEK_END); *bytes = ftell(f); fseek(f, 0, SEEK_SET);
  void *p = aligned_alloc(64, ((*bytes + 63) / 64) * 64);
  if (!p || fread(p, 1, *bytes, f) != (size_t)*bytes) { fprintf(stderr, "read error %s\n", fn); exit(1); }
  fclose(f); return p;
}

int main(int argc, char **argv) {
  if (argc != 3) { fprintf(stderr, "usage: d38_full d38_rows_i8.bin d38_axes_i64.bin\n"); return 2; }
  long b1, b2;
  V = slurp(argv[1], &b1); X = slurp(argv[2], &b2);
  if (b1 % W || b2 % (A * 8)) { fprintf(stderr, "bad file sizes\n"); return 1; }
  N = b1 / W; NX = b2 / (A * 8);
  NRM = malloc(N * sizeof *NRM); XN = malloc(NX * sizeof *XN);
  long fmt = 0;
  for (long i = 0; i < N; i++) {
    int32_t s = 0;
    for (int k = 0; k < W; k++) s += V[i][k] * V[i][k];
    for (int k = 0; k < 14; k++) if (V[i][24 + k] != V[i][38 + k]) fmt++;
    for (int k = 52; k < W; k++) if (V[i][k]) fmt++;
    if (s == 0) fmt++;
    NRM[i] = s;
  }
  int64_t big = 0;
  for (long i = 0; i < NX; i++) {
    for (int k = 0; k < A; k++) { int64_t v = X[i][k] < 0 ? -X[i][k] : X[i][k]; if (v > big) big = v; }
    for (int k = 0; k < 14; k++) if (X[i][24 + k] != X[i][38 + k]) fmt++;
  }
  if (big >= ((int64_t)1 << 60)) { fprintf(stderr, "axis coordinates too large for the exact bound\n"); return 1; }
  int32_t nmax = 0; for (long i = 0; i < N; i++) if (NRM[i] > nmax) nmax = NRM[i];
  for (long i = 0; i < NX; i++) {
    i128 s = 0; for (int k = 0; k < A; k++) s += (i128)X[i][k] * X[i][k];
    XN[i] = s;
    if (s == 0 || s != XN[0]) fmt++;                          /* nonzero, all axis rows of equal norm */
    if (s > ((i128)1 << 125) / (nmax ? nmax : 1)) fmt++;      /* keeps ni*nj and 4s^2 below 2^127 */
  }
  printf("%ld int8 rows (caps and equator), %ld axis rows, total %ld points; format errors: %ld\n", N, NX, N + NX, fmt);
  if (fmt) { printf("INVALID input\n"); return 1; }

  long cpus = sysconf(_SC_NPROCESSORS_ONLN); nth = cpus < 1 ? 1 : cpus > 256 ? 256 : (int)cpus;
  pthread_t th[256];
  for (long t = 0; t < nth; t++) pthread_create(&th[t], 0, work8, (void *)t);
  long bad = 0; double mx = -2;
  for (int t = 0; t < nth; t++) { pthread_join(th[t], 0); bad += gbad[t]; if (gmax[t] > mx) mx = gmax[t]; }
  printf("int8 pairs: %ld, max cos^2 = %.6f (limit 0.25), violations = %ld\n", N * (N - 1) / 2, mx, bad);
  for (long t = 0; t < nth; t++) pthread_create(&th[t], 0, workx, (void *)t);
  long bad2 = 0; double mx2 = -2;
  for (int t = 0; t < nth; t++) { pthread_join(th[t], 0); bad2 += gbad2[t]; if (gmax2[t] > mx2) mx2 = gmax2[t]; }
  printf("axis pairs: %ld, max cos^2 = %.6f (limit 0.25), violations = %ld\n", NX * N + NX * (NX - 1) / 2, mx2, bad2);
  long T = N + NX;
  printf("pairs checked in total: %ld\n%s kissing configuration of %ld points in R^38\n", T * (T - 1) / 2,
         bad || bad2 ? "INVALID" : "VALID", T);
  return bad || bad2;
}
