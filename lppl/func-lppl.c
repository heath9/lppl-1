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

/**
  \brief Define the LPPL function
*/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "shared.h"
#include "lppl.h"

/** Change LPPL parameter values so that:
    - \phi \in [0, \pi] 
    @pre p is not NULL */
void normalize_lppl(double *p) {
  double phi, twopi;

  assert(p != NULL);

  phi = p[6];
  twopi = M_PI + M_PI;

  if (phi < 0.)
    phi += twopi * ceil(-phi / twopi);
  /* Rounding error: phi<0 but phi/twopi=0 */
  if (phi < 0.)
    phi += twopi;
  assert(phi >= 0);

  phi -= twopi * floor(phi / twopi);
  assert(phi >= 0 && phi < 2 * M_PI);

  if (phi >= M_PI) {
    phi  -= M_PI;
    p[4]  = -p[4];
  }
  assert(phi >= 0 && phi < M_PI);

  p[6] = phi;
}

/** LPPL derivative with respect to B 
      \f[ - (T-x)^m (1 + C \cos(\omega \log(T-x) + \phi)) \f]
    @pre p is not NULL, T > x */
inline double lppl_B(double *p, double x) {
  double t, y;

  assert(p != NULL);
  assert(p[2] > x);

  t = p[2] - x;
  y = - pow(t, p[3]);
  if (p[4] != 0.) /* Avoid the computation of log(t) if C=0 */
    y *= 1. + p[4] * cos(p[5]*log(t) + p[6]);
  debug_print("b(" LG ")=" LG "\n", x, y);

  return y;
}


/** LLPL function 
      \f[ A - B (T-x)^m (1 + C \cos(\omega \log(T-x) + \phi)) \f]
    Parameters:
    - p[0]=\f$A\f$
    - p[1]=\f$B\f$
    - p[2]=\f$T\f$
    - p[3]=\f$m\f$
    - p[4]=\f$C\f$
    - p[5]=\f$\omega\f$
    - p[6]=\f$\phi\f$
    @pre p is not NULL
*/
inline double lppl(double *p, double x) {
  double y;

  assert(p != NULL);

  y = p[0] + p[1] * lppl_B(p, x);
  debug_print("lppl(" LG ")=" LG "\n", x, y);

  return y;
}


/** Step LPPL
    Parameters:
    - p[0]=\f$A\f$
    - p[1]=\f$D\f$
    - p[2]=\f$T\f$
    - p[3]=\f$m\f$
    - p[4]=\f$\alpha\f$
    - p[5]=\f$\rho\f$
    @pre p is not NULL, T > x */
double step_lppl(double *p, double x) {

  assert(p != NULL);
  assert(p[2] > x);
  
  return p[0] - p[1] * exp(p[4] * ceil(p[3] * log(p[5] * (p[2] - x)) / p[4]));
}


/** LPPL function in dlevmar vector format. 
    @pre p, hx, adata are not NULL; T > n
*/
void lppl_vector(double *p, double *hx, int m, int n, void *adata) {
  int i;  /* Loop iterator, signed as per OpenMP */
  unsigned int start;
  lppl_env_t *env;

  assert(p   != NULL);
  assert(hx  != NULL);
  assert(adata != NULL);
  assert(p[2] > (double) n);

  env = (lppl_env_t *) adata;
  start = env->start;

  /* Repeat the calculation of lppl */
  #pragma omp parallel for 
  for (i = 0; i < n; i++) 
    hx[i] = levmar_weight((int) i, (int) n, env) * 
            lppl(p, (double) (i + start));
}


/** LPPL Jacobian
    Parameters:
    - p[0]=\f$A\f$
    - p[1]=\f$B\f$
    - p[2]=\f$T\f$
    - p[3]=\f$m\f$
    - p[4]=\f$C\f$
    - p[5]=\f$\omega\f$
    - p[6]=\f$\phi\f$
    @pre T > x */
void jac_lppl(double p[], double j[], double x) {
  double gap, gap2m, log_gap, phase, cos_phase, Bj1;

  assert(p[2] > x);

  j[0] = 1;

  gap = p[2] - x;
  gap2m = pow(gap, p[3]);
  log_gap = log(gap);
  phase = p[5] * log_gap + p[6];
  cos_phase = cos(phase);

  j[1] = - gap2m * (1. + p[4] * cos_phase);
  j[6] = p[1] * p[4] * gap2m * sin(phase);
  Bj1  = p[1] * j[1];

  j[2] = (j[6] * p[5] - Bj1 * p[3]) / gap;
  j[3] = Bj1 * log_gap;
  j[4] = -p[1] * gap2m * cos_phase;
  j[5] = j[6] * log_gap;
}


/** Jacobian in vector format
    @pre p, hx, adata are not NULL; T > n */
void jac_lppl_vector(double *p, double *j, int m, int n, void *adata) {
  int i, h,  /* Loop iterator, signed as per OpenMP */
      n_int; /* n as an int */
  unsigned int start;
  double weight, *j_i;
  lppl_env_t *env;

  assert(p   != NULL);
  assert(j   != NULL);
  assert(adata != NULL);
  assert(p[2] > (double) n);

  env = (lppl_env_t *) adata;
  start = env->start;
  n_int = (int) n;

  /* Repeat the calculation of the jacobian */
  #pragma omp parallel for private(i, h, j_i, weight)
  for (i = 0; i < n_int; i++) {
    j_i = j + 7*i;
    jac_lppl(p, j_i, (double) (i + start));
    weight = levmar_weight((unsigned int) i, n, env);
    for (h = 0; h < 7; h++) 
      j_i[h] *= weight;
  }
}


/** Jacobian of the LPPL function assuming that T, m, omega, and phi
    are fixed, and that the parameters are A, B, and gamma=BC *
    @pre T > x */
void jac_linear_lppl(double p[], double j[], double x) {
  double gap;

  assert(p[2] > x);

  j[0] = 1.;
  gap = p[2] - x;
  j[1] = -pow(gap, p[3]);
  j[2] = j[1] * cos(p[5] * log(gap) + p[6]);
}


/** Jacobian of the linear part of LPPL in vector form
    @pre p, hx, adata are not NULL; T > n */
void jac_linear_lppl_vector(double p[], double j[][3], int n, lppl_env_t *env) {
  int i, h;  /* Loop iterator, signed as per OpenMP */
  unsigned int start;
  double weight, *j_i;

  assert(env != NULL);
  assert(p[2] > (double) n);

  start = env->start;

  /* Repeat the calculation of the jacobian */
  #pragma omp parallel for private(i, h, j_i, weight)
  for (i = 0; i < n; i++) {
    j_i = &(j[i][0]);
    jac_linear_lppl(p, j_i, (double) (i + start));
    weight = levmar_weight((int) i, (int) n, env);
    for (h = 0; h < 3; h++) 
      j[i][h] *= weight;
  }
}


/** Find the error implied by the current set of parameters
    @pre env is not NULL */
double solution_error(unsigned int n, double x[], double p[], lppl_env_t *env) {
  double error, *hx;

  assert(env != NULL);

  hx = (double *) scalloc((size_t) n, sizeof(double),
                          error_msg(HERE, "cannot allocate f(x) vector"));
  lppl_vector(p, hx, 7, n, (void *) env);
  error = vector_error(n, hx, x);
  free(hx);

  return error;
}

