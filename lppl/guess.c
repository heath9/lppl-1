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

/** \brief Guess the value of LPPL parameters */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "shared.h"
#include "lppl.h"


/** Compute the partial derivative \f$\part f / \part B\f$ 
    @pre b is not NULL */
static inline void dfdB(unsigned int n, double start, double p[], double *b) {
  int i;

  assert(b != NULL);

  #pragma omp parallel for
  for (i = 0; i < n; i++) 
    b[i] = lppl_B(p, start + i);
}


/** Estimate B 
    @pre b is not NULL */
static inline double B_estimate(unsigned int n, double *b, double b_avg, 
                                double x[]) {
  int i;
  double gap, *xgap, *bgap, B_numerator, B_denominator, B;

  assert(b != NULL);

  xgap = (double *) scalloc((size_t) n, sizeof(double),
                            error_msg(HERE, "cannot allocate xgap vector"));
  bgap = (double *) scalloc((size_t) n, sizeof(double), 
                            error_msg(HERE, "cannot allocate bgap vector"));

  #pragma omp parallel for private(gap)
  for (i = 0; i < n; i++) {
    gap = b[i] - b_avg;
    xgap[i] = x[i] * gap;
    bgap[i] = b[i] * gap;
  }

  B_numerator   = stable_sum(n, xgap);
  B_denominator = stable_sum(n, bgap);
  free(xgap);
  free(bgap);

  debug_print("B = " LG " / " LG "\n", B_numerator, B_denominator);
  B = B_denominator != 0. ? B_numerator / B_denominator : 1;
  if (B < 0.)
    B = 0.;

  return B;
}


/** Calculate the linear coefficients A and B of LPPL as a function of the 
    other ones via the extreme point conditions
    @pre env is not NULL */
void linear_lppl_guess(double p[], double x[], unsigned int n, 
                       lppl_env_t *env) {
  double *b, b_avg, x_avg;

  assert(env != NULL);

  /* Compute the partial derivative df / dB */
  b = (double *) scalloc((size_t) n, sizeof(double),
                         error_msg(HERE, "cannot allocate df/dB"));
  dfdB(n, (double) env->start, p, b);
  b_avg = stable_sum(n, b) / (double) n;
  debug_print("b_avg = " LG "\n", b_avg);

  /* Estimate B */
  p[1] = B_estimate(n, b, b_avg, x);

  /* Estimate A */
  x_avg = stable_sum(n , x) / (double) n;
  p[0] = x_avg - p[1] * b_avg;

  free(b);
}


/** Guess values of some LPPL parameters 
    - p[0]=\f$A\f$ (calculated)
    - p[1]=\f$B\f$ (calculated)
    - p[2]=\f$T\f$ (calculated)
    - p[3]=\f$m\f$ (given)
    - p[4]=\f$C\f$ (given)
    - p[5]=\f$\omega\f$ (calculated)
    - p[6]=\f$\phi\f$ (given) 
    @pre env is not NULL */
void lppl_guess(double p[], double x[], double rho, double xmin, double xmax, 
                unsigned int n, lppl_env_t *env) {
  assert(env != NULL);

  p[5] = (M_PI + M_PI) / log(rho);
  p[2] = (rho * xmax - xmin) / (rho - 1.);
  if (p[2] <= n)
    return;

  linear_lppl_guess(p, x, n, env);

  if (env->verbose) {
    printf("\nInitial guess\n");
    print_lppl_parms(p);
  }
}


/* Declaration of lapack dposvx(1) function */
extern int dposvx_(char *fact, char *uplo, int *n, int *nrhs, 
                   double *a, int *lda, double *af, int *ldaf, 
                   char *equed, double *s, double *b, int *ldb, 
                   double *x, int *ldx, double *rcond, double *ferr, 
                   double *berr, double *work, int *iwork, int *info);
extern int dgels_(char *trans, int *m ,int *n, int *nhrs, double *a, int *lda,
                  double *b, int *ldb, double *work, int *lwork, int *info);


/** Compute the left- and right-hand side of the linear least-square problem */
static void linear_lsq_terms(int n, double x[], double j[][3], 
                             double a[][3], double b[]) { 
  int i, l, t;

  #pragma omp parallel for private(i,l,t)
  for (i = 0; i < 3; i++) {
    b[i] = 0.;
    for (t = 0; t < n; t++)
      b[i] += x[t] * j[t][i];
    for (l = i; l < 3; l++) {
      a[i][l] = 0.;
      for (t = 0; t < n; t++)
        a[i][l] += j[t][i] * j[t][l];
    }
  }
}


/** Convert linear LPPL parameters (if available) into regular parameters */
static void linear2lppl(int info, double p[], double q[]) {

  if (info != 0. ||       /* The linear system cannot be solved */ 
      q[1] <= 0.)         /* Constraint violated: B<=0 */ 
    return;

  p[0] = q[0];
  p[1] = q[1];
  p[4] = q[2] / q[1];
  normalize_lppl(p);
}


/** Print the linear parameters */
void print_qparms(double q[]) {
  unsigned int i;

  for (i = 0; i < 3; i++)  {
    debug_print("q[%u] = " LG "\n", i, q[i]);
  }
}


/** Solve the least-square problem assuming that the non-linear components
    are fixed 
    @pre env is not NULL
    @return info, as defined in dposvx(1) */
int linear_lsq(double p[], double x[], int n, lppl_env_t *env) {
  double (*j)[3]; /* Linear jacobian */

  /* lapack variables */
  char char_E, char_L, equed; 
  int num_1, num_3, info, iwork[3]; 
  double a[3][3], b[3], af[3][3], s[3], q[3], rcond, ferr, berr, work[9];

  /* Compute the Jacobian assuming that the non-linear terms are constant */
  j = (double (*)[3]) 
      scalloc((size_t) n, sizeof(double[3]),
              error_msg(HERE, "cannot allocate linear jacobian"));
  jac_linear_lppl_vector(p, j, n, env);

  /* Compute matrix A and vector b for solving the linear least-square */
  linear_lsq_terms(n, x, j, a, b);
  free(j);

  /* Lapack codes */
  char_E = 'E';
  char_L = 'L';
  num_1  = 1;
  num_3  = 3;

  dposvx_(&char_E, &char_L, &num_3, &num_1, &(a[0][0]), &num_3,
          &(af[0][0]), &num_3, &equed, s, b, &num_3, q, &num_3,
          &rcond, &ferr, &berr, work, iwork, &info);

  
  print_qparms(q);

  linear2lppl(info, p, q);

  return info;
}


/** Free variables allocated in exp_lsq */
void free_exp_lsq(double (*a)[2], double *b, double *work) {
  free(a);
  free(b);
  free(work);
}


/** Allocate and initialize the LSQ matrix A 
    @return A */
double (*(exp_lsq_a(unsigned int n, unsigned int start)))[2] {
  int i, n_int;
  double start_double, (*a)[2];
  
  /* Change argument type */
  n_int = (int) n;
  start_double = (double) start;

  a = (double (*)[2]) 
      scalloc((size_t) n, sizeof(double[2]),
              error_msg(HERE, "cannot allocate LSQ matrix"));

  #pragma omp parallel for 
  for (i = 0; i < n_int; i++) {
    a[i][0] = i + start_double; 
    a[i][1] = 1.;
  }

  for (i = 0; i < n_int; i++)
    debug_print("a[%i] = " LG "\t" LG "\n", i, a[i][0], a[i][1]);

  return a;
}


/** Allocate and initialize the right hand side vector 
    @return b */
double *exp_lsq_b(unsigned int n, double x[]) {
  int i, n_int;  /* n, as an int */
  double *b;

  n_int = (int) n;

  b = (double *) 
      scalloc((size_t) n, sizeof(double),
              error_msg(HERE, "cannot allocate LSQ RHS"));

  #pragma omp parallel for 
  for (i = 0; i < n_int; i++)
    b[i] = x[i];

  for (i = 0; i < n_int; i++)
    debug_print("b[%i] = " LG "\n", i, b[i]);

  return b;
}


/* Print exp LSQ data */
static void print_exp_lsq(unsigned int n, double a[][2], double b[]) {
  unsigned int i;

  debug_print("a_exp = " LG " \n", b[0]);
  debug_print("b_exp = " LG " \n", b[1]);
  for (i = 0; i < n; i++)
    debug_print("a[%i] = " LG "\t" LG "\n", i, a[i][0], a[i][1]);
  for (i = 2; i < n; i++)
    debug_print("b[%i]= " LG "\n", i, b[i]);
}


/** Fits prices to an exponential function via linear least squares 
    @pre env is not NULL
    @return info, as defined in dgels(1), or -12 if B<0 */
int exp_lsq(double p[], double x[], int n, lppl_env_t *env) {
  /* lapack variables */
  char char_T;
  int n_int, num_1, num_2, lwork, info;
  double (*a)[2],    /* LSQ left hand side */
         *work, *b, lwork_opt;

  assert(env != NULL);

  #pragma omp parallel 
  { 
  #pragma omp sections nowait 
  {
    #pragma omp section
    a = exp_lsq_a(n, env->start);

    #pragma omp section
    b = exp_lsq_b(n, x);
  } }

  /* Lapack codes */
  char_T = 'T'; /* Must transponse in C-to-Fortran */
  num_1 = 1;
  num_2 = 2;

  /* Find the optimal value of the work area */
  n_int = (int) n;
  lwork = -1;
  work  = NULL;
  dgels_(&char_T, &num_2, &n_int, &num_1, &(a[0][0]), &num_2, b, &n_int, 
         &lwork_opt, &lwork, &info);
  if (info < 0) {
    free_exp_lsq(a, b, work);
    return info;
  }

  /* Allocate the work area */
  lwork = (int) ceil(lwork_opt);
  debug_print("LSQ work area: %i\n", lwork);
  work = (double *) 
         scalloc((size_t) lwork, sizeof(double),
                 error_msg(HERE, "cannot allocate LSQ work area"));

  /* Solve the LSQ problem */
  dgels_(&char_T, &num_2, &n_int, &num_1, &(a[0][0]), &num_2, b, &n_int, 
         &(work[0]), &lwork, &info);
  if (b[0] < 0.) { /* B is out of bound */
    print_exp_lsq(n, a, b);
    info = -12;
  }
  if (info < 0) { /* error in solving the LSQ problem  or B < 0*/
    free_exp_lsq(a, b, work);
    return info;
  }

  print_exp_lsq(n, a, b);

  /* Convert the LSQ solution into LPPL parameters */
  p[1] = b[0];
  p[0] = b[1] + p[1] * (p[2] - 1. + (double) env->start);

  free_exp_lsq(a, b, work);

  return info;
}


/** Exponential least-square, safely
    @rpe env is not NULL  */
void safe_exp_lsq(double p[], double x[], int n, lppl_env_t *env) {
  int info;

  if (!env->expfit) 
    return;

  /* Enforce ancillary parameters to define an exponential solution */
  p[3] = 1.; /* m = 1 in exponential solution */
  p[4] = 0;  /* C = 0 in exponential solution */

  info = exp_lsq(p, x, n, env);

  if (env->verbose) {
    printf("\nLeast-square fit to exponential");
    if (info != 0)
      printf(" failed: %i\n", info);
    else {
      printf(":\n");
      print_lppl_parms(p);
      printf("Average error: " LG "\n", 
             solution_error(n, x, p, env) / (double) (n - 7.));
    }
  }
}
