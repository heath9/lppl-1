/* Copyright (c) 2011
 * Vincenzo Liberatore, Case Western Reserve University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 * - Neither the name of Case Western Reserve University nor the names of
 *   its contributors may be used to endorse or promote products derived from
 *   this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 */

/* Test cases */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "lm.h"
#include "shared.h"
#include "lppl.h"

int signal_code = 0;

/* Signal handler */
void handler(int signum) {
  signal_code = 1;
}

#define errhndlr(cmdpred, errmsg) \
  if (cmdpred) { \
    fprintf(stderr, errmsg); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
  }

/* Weight functions */

/* Identity function */
double identity(int i, int j) {
  return 1.;
}

#define W 20.
double quadratic(int i, int n) {
  assert(i < n);

  return W / (double) (n - i + 1);
}


void test_lppl() {
  double p[7];

  p[0] = 1.;
  p[1] = 2.;
  p[2] = 9.;
  p[3] = .5;
  p[4] = .1;
  p[5] = M_PI / 4.;
  p[6] = M_PI / 8.;

  CU_ASSERT_DOUBLE_EQUAL(lppl(p, 0.),  -4.6876, 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(lppl(p, 1.),  -4.4082, 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(lppl(p, 4.),  -3.4337, 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(lppl(p, 5.5), -2.8139, 1e-3);
}


void test_levmar_dif() {
  int iterations;
  unsigned int i, j, n;
  double p[7], x[200], covar[7][7], info[LM_INFO_SZ], lb[7], ub[7];
  FILE *pricefile;

  pricefile = safe_fopen("gld.plt", "cannot open test file gld.plt");
  CU_ASSERT_PTR_NOT_NULL(pricefile);
  n = read_price(200, x, pricefile);
  CU_ASSERT_EQUAL(n, 159);
  fclose(pricefile);

  p[0] = 1031.28;
  p[1] = 1026.35;
  p[2] = 175.5471;
  p[3] = 9.19171278272457e-05;
  p[4] = -1.63403e-05;
  p[5] = 11.0322267827269;
  p[6] = 5. * M_PI / 8.;

  for (i = 0; i < 7; i++) {
    lb[i] = -DBL_MAX;
    ub[i] =  DBL_MAX;
  }
  lb[1] = 0.;
  lb[2] = n + 1.;
  lb[3] = 0.;
  ub[3] = 1.;
  lb[5] = 0.;
  lb[6] = 0.;
  ub[6] = M_PI;
  

  iterations = dlevmar_bc_dif(lppl_vector, p, x, 7, n, lb, ub, 10000, NULL, 
                              info, NULL, &(covar[0][0]), NULL);

  printf("\nIterations %i\n", iterations);
  for (i = 0; i < 7; i++) 
    printf("p[%u] = " LG "\n", i, p[i]);
  for (i = 0; i < LM_INFO_SZ; i++) 
    printf("info[%u]: " LG "\n", i, info[i]);
  printf("Covariance\n");
  for (i = 0; i < 7; i++) {
    for (j = 0; j < 7; j++)
      printf(LG "  ", covar[i][j]);
    printf("\n");
  }
}

void test_guess() {
  double x[3], p[7];
  lppl_env_t env;

  env.start = 0;
  env.verbose = 1;
  env.window = 0.;
  env.expfit = 0;

  p[3] = 2e-4;
  p[4] = 1e-3;
  p[6] = 0.;

  x[0] = 4.5;
  x[1] = 4.6;
  x[2] = 5;

  lppl_guess(p, x, 1.1, 1., 2., 3, &env);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 821.3960050306066, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 815.7654577211625, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[2], 11.99999999999999, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[3], 2e-4, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[4], 1e-3, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[5], 65.9235489858395, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[6], 0, 1e-5);
}


void test_parms() {
  unsigned int state;
  double p[7];
  FILE *parmsfile;

  parmsfile = safe_fopen("gld.par", "cannot open test parameters gld.par");
  state = read_lppl_parms(parmsfile, p, 1);
  CU_ASSERT_EQUAL(state, 1);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 1031.28, 1e-1);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 1026.35, 1e-1);
  CU_ASSERT_DOUBLE_EQUAL(p[2], 175.5471, 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(p[3], 9.19171278272457e-05, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[4], -1.63403e-05, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[5], 11.0322267827269, 1e-5);
  CU_ASSERT_DOUBLE_EQUAL(p[6], 1.9635, 1e-3);

  state = read_lppl_parms(parmsfile, p, 1);
  CU_ASSERT_EQUAL(state, 1);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 12.1, 1e-1);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 1026.35, 1e-1);

  state = read_lppl_parms(parmsfile, p, 1);
  CU_ASSERT_EQUAL(state, 0);

  fclose(parmsfile);
}

void test_mmap() {
  int *x, status;
  pid_t kidpid;

  printf("\n");

  x = (int *) scmmap("cannot mmap", 3, sizeof(int));
  x[0] = x[1] = 0;

  kidpid = fork();

  if (kidpid) {
    /* Parent */
    while (x[1] != 2) sleep(1);
    CU_ASSERT_EQUAL(x[0], 1);
  } else {
    /* child */
    sleep(1);
    printf("x0\n");
    x[0] = 1;
    sleep(2);
    printf("x1\n");
    x[1] = 2;
    scmunmap("cannot unmap", (void *) x, 3, sizeof(int));
    exit(EXIT_SUCCESS);
  }

  CU_ASSERT_EQUAL(swait("can't wait", &status), kidpid);
  CU_ASSERT_EQUAL(status, 0);
  scmunmap("cannot unmap", (void *) x, 3, sizeof(int));
}


void test_signal() {
  safe_signal("cannot install handler", SIGINT, handler);

  sleep(5);

  safe_kill("cannot kill me", 0, SIGINT);

  while(signal_code == 0) {
    sleep(1);
  }

  printf("\nReceived SIGINT\n");
}


void test_vector() {
  unsigned int i;
  double x[5], y[5];

  for (i = 0; i < 5; i++) {
    x[i] = cos(i);
    y[i] = sin(i);
  }

  CU_ASSERT_DOUBLE_EQUAL(vector_error(5, x ,y), 4.13756232, 1e-5);
}


extern int dposvx_(char *fact, char *uplo, int *n, int *nrhs, 
                   double *a, int *lda, double *af, int *ldaf, 
                   char *equed, double *s, double *b, int *ldb, 
                   double *x, int *ldx, double *rcond, double *ferr, 
                   double *berr, double *work, int *iwork, int *info);

/* Linear solution: testing */
void test_linear() {
  char char_E = 'E';
  char char_L = 'L';
  int num_1 = 1;
  int num_3 = 3;
  char equed;
  int info, iwork[3];
  unsigned int i, j, k, n;
  double d[3][3], 
         a[3][3], b[3], af[3][3], s[3], x[3], rcond, ferr, berr, work[9],
         value,
         p[7], j0[3], y[200], jv[200][3], z[200], err;
  /* int state; */
  FILE *pricefile, *parmsfile;
  lppl_env_t env;

  debug_print("\n");

  for (i = 0; i < 3; i++)
    b[i] = 1.;
  d[0][0] = 1.;
  d[0][1] = 2.;
  d[0][2] = 4.;
  d[1][0] = 2.;
  d[1][1] = -3.;
  d[1][2] = 4.;
  d[2][0] = -5.;
  d[2][1] = 7.;
  d[2][2] = -3.;
  for (i = 0; i < 3; i++) {
    /* for (j = i; j <= 3; j++) { */
    for (j = 0; j <= 3; j++) {
      a[i][j] = 0.;
      for (k = 0; k < 3; k++) 
        a[i][j] += d[i][k] * d[j][k];
      debug_print("a%u,%u=" LG "\n", i, j, a[i][j]);
    }
  }

  dposvx_(&char_E, &char_L, &num_3, &num_1, &(a[0][0]), &num_3, 
          &(af[0][0]), &num_3, &equed, s, b, &num_3, x, &num_3, 
          &rcond, &ferr, &berr, work, iwork, &info);

  debug_print("info: %i\n", info);
  debug_print("equed: %c\n", equed);
  debug_print("rcond: " LG "\n", rcond);
  debug_print("ferr: " LG "\n", ferr);
  debug_print("berr: " LG "\n", berr);
  for (i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      debug_print("a[%u][%u]= " LG "\n", i, j, a[i][j]);
  for (i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      debug_print("af[%u][%u]= " LG "\n", i, j, af[i][j]);
  for (i = 0; i < 3; i++) 
    debug_print("b[%u] = %g\n", i, b[i]);
  for (i = 0; i < 3; i++) 
    debug_print("x[%u] = %g\n", i, x[i]);
  for (i = 0; i < 3; i++) 
    debug_print("s[%u] = %g\n", i, s[i]);

  for (i = 0; i < 3; i++) {
    value = 0.;
    for (j = 0; j < 3; j++)
      value += a[i][j] * x[j];
    debug_print("value %u=" LG "\n", i, value);
    CU_ASSERT_DOUBLE_EQUAL(value, 1., 1e-7);
  }

  p[0] = p[1] = 1.;
  p[2] = 1000.;
  p[3] = 1e-05;
  p[4] = -1e-05;
  p[5] = 6.28;
  p[6] = 0.;

  pricefile = safe_fopen("gld.plt", "cannot open test file gld.plt");
  CU_ASSERT_PTR_NOT_NULL(pricefile);
  n = read_price(200, y, pricefile);
  CU_ASSERT_EQUAL(n, 159);
  fclose(pricefile);

  jac_linear_lppl(p, j0, 1.);
  CU_ASSERT_DOUBLE_EQUAL(j0[0], 1., 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j0[1], -1.00006906993300, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j0[2], -0.820918984169654, 1e-7);

  env.start = 0;
  env.verbose = 1;
  env.window = 0.;
  env.expfit = 0;
  jac_linear_lppl_vector(p, jv, n, &env);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][0], 1., 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][1], -1.00006906993300, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][2], -0.820918984169654, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][0], 1., 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][1], -1.00006902980982, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][2], -0.806269048302848, 1e-7);

  debug_print("\n");
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  debug_print("Error " LG "\n", vector_error(n, y, z));
  linear_lsq(p, y, n, &env);
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 1388144.777867017, 1e-7);

  parmsfile = safe_fopen("gld.par", "cannot open test parameters gld.par");
  /* state = */ read_lppl_parms(parmsfile, p, 1);
  fclose(parmsfile);
  debug_print("\n");
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 1257353.700889399, 1e-7);
  linear_lsq(p, y, n, &env);
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 727.6028524813674, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 70416.08024416743, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 70294.16260559887, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[2], 175.5471, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[3], 9.19171278272457e-05, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[4], 4.211726262681084e-06 , 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[5], 11.0322267827269, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[6], 1.9635, 1e-7);
}


void test_weight() {
  /* int state; */
  unsigned int n;
  double p[7], y[200], j[1400], jv[200][3], err, z[200], weight;
  FILE *pricefile, *parmsfile;
  lppl_env_t env;

  printf("\n"); 

  env.start = 0;
  env.verbose = 1;
  env.window = 20.;
  env.expfit = 0;

  p[0] = p[1] = 1.;
  p[2] = 1000.;
  p[3] = 1e-05;
  p[4] = -1e-05;
  p[5] = 6.28;
  p[6] = 0.;

  pricefile = safe_fopen("gld.plt", "cannot open test file gld.plt");
  CU_ASSERT_PTR_NOT_NULL(pricefile);
  n = read_price(200, y, pricefile);
  CU_ASSERT_EQUAL(n, 159);
  fclose(pricefile);

  jac_lppl_vector(p, j, 7, n, (void *) &env);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7], 0.1197604790419162, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+1], -0.119767804392718, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+2], 5.983770816431365e-09, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+3], -0.8258807764981168, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+4], -0.09332365698992041, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+5], 5.176469888939391e-06, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(j[12*7+6], 7.50682726545702e-07, 1e-7);

  jac_linear_lppl_vector(p, jv, n, &env);
  weight = levmar_weight(1, n, &env);
  CU_ASSERT_DOUBLE_EQUAL(weight, .11235955056179775280, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][0], weight, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][1], -1.00006906993300*weight, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[1][2], -0.820918984169654*weight, 1e-7);
  weight = levmar_weight(5, n, &env);
  CU_ASSERT_DOUBLE_EQUAL(weight, .11494252873563218390, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][0], 1.*weight, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][1], -1.00006902980982*weight, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(jv[5][2], -0.806269048302848*weight, 1e-7);

  debug_print("\n");
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  debug_print("Error " LG "\n", vector_error(n, y, z));
  linear_lsq(p, y, n, &env);
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 1388143.417988638, 1e-7);

  env.window=5000.;
  parmsfile = safe_fopen("gld.par", "cannot open test parameters gld.par");
  /* state = */ read_lppl_parms(parmsfile, p, 1);
  fclose(parmsfile);
  debug_print("\n");
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 1259294.83881059, 1e-7);
  linear_lsq(p, y, n, &env);
  print_lppl_parms(p);
  lppl_vector(p, z, 7, n, (void *) &env);
  err = vector_error(n, y, z);
  debug_print("Error " LG "\n", err);
  CU_ASSERT_DOUBLE_EQUAL(err, 777.1276490803839, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 56313.68096396694, 1e-1);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 56196.037889812, 1e-1);
  CU_ASSERT_DOUBLE_EQUAL(p[2], 175.5471, 1e-2);
  CU_ASSERT_DOUBLE_EQUAL(p[3], 9.19171278272457e-05, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[4], 4.992633607783836e-06 , 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[5], 11.0322267827269, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[6], 1.9635, 1e-2);
}


void test_iter() {
  int levmar_stat_iter[2];
  unsigned int i;
  double levmar_stat_time[2];
  levmar_iter_t lmiter;

  levmar_stat_time[1] = -1.;
  levmar_stat_iter[1] = -1;
  printf("\n");

  iter_init(4,  &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 0);
  CU_ASSERT_EQUAL(lmiter.iterations, 4);

  iter_init(0, &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 0);
  CU_ASSERT_EQUAL(lmiter.iterations, 1);

  iter_init(-4, &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 4);
  CU_ASSERT_EQUAL(lmiter.start_up, 1);

  iter_update(.1, 1., 4, 1., .01, levmar_stat_time, levmar_stat_iter, &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 8);
  CU_ASSERT_EQUAL(lmiter.start_up, 1);
  
  iter_update(.15, .1, 8, .1, 1., levmar_stat_time, levmar_stat_iter, &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 7);
  CU_ASSERT_EQUAL(lmiter.start_up, 0);
  
  iter_update(.14, .1, 7, 1., 1.43, levmar_stat_time, levmar_stat_iter, 
              &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 6);
  CU_ASSERT_EQUAL(lmiter.start_up, 0);
  
  iter_update(.12, .1, 6, 1., .8, levmar_stat_time, levmar_stat_iter, &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 7);
  CU_ASSERT_EQUAL(lmiter.start_up, 0);

  for (i = 0; i < 10; i++) 
    iter_update(.14, .1, 7, 1., 1.43, levmar_stat_time, levmar_stat_iter,
                &lmiter);
  CU_ASSERT_EQUAL(lmiter.adaptive, 1);
  CU_ASSERT_EQUAL(lmiter.iterations, 1);
  CU_ASSERT_EQUAL(lmiter.start_up, 0);
}


void test_exp_lsq() {
  int state;
  unsigned int n;
  double p[7], y[200];
  FILE *pricefile;
  lppl_env_t env;

  printf("\n");

  env.start = 1;
  env.verbose = 1;
  env.window = 0.;
  env.expfit = 1;

  p[0] = p[1] = 0.;
  p[2] = 400.;

  y[0] = 2;
  y[1] = 3;
  y[2] = 4;
  n = 3;

  state = exp_lsq(p, y, n, &env);
  CU_ASSERT_EQUAL(state, 0);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 401., 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 1., 1e-3);

  p[0] = p[1] = 0.;
  p[2] = 400.;
  p[4] = .5;
  p[5] = .1;
  safe_exp_lsq(p, y, n, &env);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 401., 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 1., 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(p[3], 1., 1e-3);
  CU_ASSERT_DOUBLE_EQUAL(p[4], 0., 1e-3);

  pricefile = safe_fopen("gld.plt", "cannot open test file gld.plt");
  CU_ASSERT_PTR_NOT_NULL(pricefile);
  n = read_price(200, y, pricefile);
  CU_ASSERT_EQUAL(n, 159);
  fclose(pricefile);
  vector_f(n, y, log);

  state = exp_lsq(p, y, n, &env);
  CU_ASSERT_DOUBLE_EQUAL(p[0], 4.795776987795519, 1e-7);
  CU_ASSERT_DOUBLE_EQUAL(p[1], 0.0008147204362796475, 1e-7);
}


void test_rnd() {
  unsigned int i;
  double sum, sum_square, x[2], b[10];
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
  lppl_env_t env;
  
  printf("\n");

  ssrand48(&mutex, 123);

  sum = sum_square = 0.;
  for (i = 0; i < 10000; i += 2) {
    generate_normal(x, &mutex);
    debug_print(LG "\t" LG "\n", x[0], x[1]);
    sum += x[0] + x[1];
    sum_square += x[0] * x[0] + x[1] * x[1];
  }

  sum /= 10000.;
  sum_square /= 10000.;

  CU_ASSERT_DOUBLE_EQUAL(sum, 0., 1e-2);
  CU_ASSERT_DOUBLE_EQUAL(sum_square, 1., 1e-2);

  printf(LG "\t" LG "\n", sum, sum_square);

  env.diffuse = 2.;
  env.relax = 1.;
  ssrand48(&mutex, 123);

  ou_noise(10, b, &mutex, &env);
  for (i = 0; i < 10; i++)
    printf("OU[%u] = " LG "\n", i, b[i]);
  CU_ASSERT_DOUBLE_EQUAL(b[0], 1.183129757096921, 1e-9);
  CU_ASSERT_DOUBLE_EQUAL(b[4], -1.039502098091391, 1e-9);
  CU_ASSERT_DOUBLE_EQUAL(b[8], -1.747013517001366, 1e-9);
  CU_ASSERT_DOUBLE_EQUAL(b[9], 0.5720424643802948, 1e-9);
}


int main(int argc, const char **argv) {
  /* CU_ErrorCode state; */
  CU_pSuite psuite, sys_suite; 
  CU_pTest ptest;

  argv0_init(argv[0]);

  errhndlr(CU_initialize_registry() != CUE_SUCCESS,
           "CUnit cannot initialize registry");

  errhndlr((psuite = CU_add_suite("lppl", NULL, NULL)) == NULL,
           "CUnit cannot initialize test suite");

  errhndlr((sys_suite = CU_add_suite("system", NULL, NULL)) == NULL,
           "CUnit cannot initialize system test suite");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_lppl)) == NULL,
              "CUnit cannot add test 1: lppl function\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_guess)) == NULL,
              "CUnit cannot add test 2: guess parameters\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_parms)) == NULL,
              "CUnit cannot add test 3: read parameters\n");

  errhndlr((ptest = CU_ADD_TEST(sys_suite, test_mmap)) == NULL,
              "CUnit cannot add test 4: mmap\n");

  errhndlr((ptest = CU_ADD_TEST(sys_suite, test_signal)) == NULL,
              "CUnit cannot add test 5: signal\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_vector)) == NULL,
              "CUnit cannot add test 6: vector operations\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_linear)) == NULL,
              "CUnit cannot add test 7: linear algebra (lapack)\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_weight)) == NULL,
              "CUnit cannot add test 8: weighted systems\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_iter)) == NULL,
              "CUnit cannot add test 9: adjust iterations\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_exp_lsq)) == NULL,
              "CUnit cannot add test 10: exp lsq\n");

  errhndlr((ptest = CU_ADD_TEST(psuite, test_rnd)) == NULL,
              "CUnit cannot add test 11: random numbers\n");

  CU_basic_set_mode(CU_BRM_VERBOSE);
  /* state = */ CU_basic_run_tests();

  CU_cleanup_registry();

  return 0;
}
