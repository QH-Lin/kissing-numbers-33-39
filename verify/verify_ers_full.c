/* Brute-force check of the Edel-Rains-Sloane three-layer kissing configuration in R^n (n <= 48).
   Builds every vector explicitly and tests every pair x != y for  <x,y> <= |x||y|/2  in exact integer arithmetic
   (<x,y> <= 0, or 4<x,y>^2 <= |x|^2 |y|^2). A repeated vector would fail this test.
     dense layer   (+-1) on coordinates 0..m-1 from a binary code (m = 32 for "pad32", m = n for "coset")
     middle layer  +-2 with an even number of minus signs on each weight-8 support
     root layer    +-4 e_i +-4 e_j
   usage:  verify_ers_full n pad32 gen32_17.txt supports.txt                 (0/1 rows, leftmost = coordinate 0)
           verify_ers_full n coset kernel.txt representatives.txt supports.txt   (hex rows, bit i = coordinate i)
   build:  cc -O3 -pthread -o verify_ers_full verify/verify_ers_full.c   (Apple silicon: add -mcpu=native) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#if defined(__ARM_FEATURE_DOTPROD)
#include <arm_neon.h>
#endif
#define L 48
static int8_t (*V)[L];
static int32_t *NRM;
static long N;
static int nth;
static long gbad[256];
static double gmax[256];

static inline int32_t dot(const int8_t *a, const int8_t *b) {
#if defined(__ARM_FEATURE_DOTPROD)
  int32x4_t acc = vdupq_n_s32(0);
  acc = vdotq_s32(acc, vld1q_s8(a), vld1q_s8(b));
  acc = vdotq_s32(acc, vld1q_s8(a + 16), vld1q_s8(b + 16));
  acc = vdotq_s32(acc, vld1q_s8(a + 32), vld1q_s8(b + 32));
  return vaddvq_s32(acc);
#else
  int32_t s = 0;
  for (int k = 0; k < L; k++) s += a[k] * b[k];
  return s;
#endif
}

static void *work(void *arg) {
  long t = (long)arg, bad = 0;
  double mx = -2;
  for (long i = t; i < N; i += nth) {
    int64_t ni = NRM[i];
    for (long j = i + 1; j < N; j++) {
      int64_t s = dot(V[i], V[j]);
      if (s > 0) {
        int64_t nj = NRM[j];
        if (4 * s * s > ni * nj) bad++;
        double c = (double)s * s / ((double)ni * nj);
        if (c > mx) mx = c;
      }
    }
  }
  gbad[t] = bad; gmax[t] = mx;
  return 0;
}

static int read_rows(const char *fn, uint64_t *out, int max, int hexmode, int len) {
  FILE *f = fopen(fn, "r");
  if (!f) { perror(fn); exit(1); }
  char line[512];
  int n = 0;
  while (fgets(line, sizeof line, f)) {
    if (line[0] == '#' || line[0] == '$' || line[0] == '\n' || line[0] == '\r') continue;
    if (n == max) { fprintf(stderr, "too many rows in %s\n", fn); exit(1); }
    if (hexmode) { out[n++] = strtoull(line, 0, 16); continue; }
    uint64_t m = 0; int k = 0;
    for (char *p = line; *p; p++)
      if (*p == '0' || *p == '1') { if (*p == '1') m |= 1ull << k; k++; }
    if (k == 0) continue;
    if (len && k != len) { fprintf(stderr, "bad row length %d in %s\n", k, fn); exit(1); }
    out[n++] = m;
  }
  fclose(f);
  return n;
}

int main(int argc, char **argv) {
  if (argc < 5) { fprintf(stderr, "usage: see the comment at the top of verify_ers_full.c\n"); return 2; }
  int n = atoi(argv[1]);
  static uint64_t G[64], Rp[4096], S[20000];
  int k, nr = 1, ns, m;
  const char *supf;
  if (n < 32 || n > L) { fprintf(stderr, "n must be in 32..%d\n", L); return 2; }
  if (!strcmp(argv[2], "pad32")) { k = read_rows(argv[3], G, 64, 0, 32); Rp[0] = 0; m = 32; supf = argv[4]; }
  else if (!strcmp(argv[2], "coset") && argc >= 6) { k = read_rows(argv[3], G, 64, 1, 0); nr = read_rows(argv[4], Rp, 4096, 1, 0); m = n; supf = argv[5]; }
  else { fprintf(stderr, "mode must be pad32 or coset\n"); return 2; }
  ns = read_rows(supf, S, 20000, 0, n);
  for (int i = 0; i < k; i++) if (G[i] >> m) { fprintf(stderr, "generator row wider than %d\n", m); return 1; }
  for (int i = 0; i < nr; i++) if (Rp[i] >> m) { fprintf(stderr, "representative wider than %d\n", m); return 1; }

  long nd = (long)nr << k, nmid = 128L * ns, nroot = 2L * n * (n - 1);
  N = nd + nmid + nroot;
  printf("n=%d  dense: %d coset(s) x 2^%d = %ld on %d coordinates; middle: 128 x %d = %ld; roots: %ld; total N = %ld\n",
         n, nr, k, nd, m, ns, nmid, nroot, N);
  V = aligned_alloc(64, ((N * L + 63) / 64) * 64);
  NRM = malloc(N * sizeof *NRM);
  if (!V || !NRM) { fprintf(stderr, "out of memory\n"); return 1; }
  memset(V, 0, N * L);
  long c = 0;
  for (int r = 0; r < nr; r++)
    for (long w = 0; w < (1L << k); w++) {
      uint64_t cw = Rp[r];
      for (int b = 0; b < k; b++) if (w >> b & 1) cw ^= G[b];
      for (int j = 0; j < m; j++) V[c][j] = (cw >> j & 1) ? -1 : 1;
      NRM[c++] = m;
    }
  for (int s = 0; s < ns; s++) {
    int idx[8], q = 0;
    if (__builtin_popcountll(S[s]) != 8) { fprintf(stderr, "support %d has weight != 8\n", s); return 1; }
    for (int j = 0; j < n; j++) if (S[s] >> j & 1) idx[q++] = j;
    for (int p = 0; p < 256; p++) {
      if (__builtin_popcount(p) & 1) continue;
      for (int t = 0; t < 8; t++) V[c][idx[t]] = (p >> t & 1) ? -2 : 2;
      NRM[c++] = 32;
    }
  }
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      for (int a = 0; a < 4; a++) {
        V[c][i] = (a & 1) ? -4 : 4; V[c][j] = (a & 2) ? -4 : 4;
        NRM[c++] = 32;
      }
  for (long i = 0; i < N; i++) {
    int s = 0;
    for (int j = 0; j < L; j++) s += V[i][j] * V[i][j];
    if (s != NRM[i]) { fprintf(stderr, "norm error at vector %ld\n", i); return 1; }
  }

  long cpus = sysconf(_SC_NPROCESSORS_ONLN);
  nth = cpus < 1 ? 1 : cpus > 256 ? 256 : (int)cpus;
  pthread_t th[256];
  for (long t = 0; t < nth; t++) pthread_create(&th[t], 0, work, (void *)t);
  long bad = 0;
  double mx = -2;
  for (int t = 0; t < nth; t++) { pthread_join(th[t], 0); bad += gbad[t]; if (gmax[t] > mx) mx = gmax[t]; }
  printf("pairs checked: %ld, max cos^2 between distinct vectors = %.6f (limit 0.25), violations = %ld\n",
         N * (N - 1) / 2, mx, bad);
  printf("%s kissing configuration of %ld points in R^%d\n", bad ? "INVALID" : "VALID", N, n);
  return bad != 0;
}
