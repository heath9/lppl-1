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

/** \brief Main for gen-lppl */

#include <assert.h>
#include <math.h>
#include <omp.h>
#include <limits.h>
#include "shared.h"
#include "lppl.h"

/** Round a price x down to penny value, but not less than 1c */
static double penny_round(double x) {
  double result;

  result = floor(x * 100.) / 100.;
  if (result <= 0.)
    result = .01;

  return result;
}


/** Normalize prices
    @pre env is not NULL  */
static void normalize_price(unsigned int n, double x[], lppl_env_t *env) {

  assert(env != NULL);

  /* Exponentiate */
  if (!env->linprice)  
    vector_f(n, x, exp);

  /* Round down to the penny */
  if (!env->nopenny) 
    vector_f(n, x, penny_round);
}


int main(int argc, const char **argv) {
  int i,                   /* Loop counter */
      n_int;               /* n, as an int */
  unsigned int n;          /* Price vector length */
  double p[7],             /* LPPL parameters */
         x[N],             /* LPPL time series vector */
         y[N],             /* Temporary price vector */
         z[N],             /* Price vector */
         b[N],             /* Brownian motion time series */
         ou[N];            /* Ornstein-Uhlbeck time series */
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; /* drand mutex */
  lppl_env_t env;

  /* Initialize error handling */
  argv0_init(argv[0]);

  /* Set up default environment */
  env.verbose = env.nopenny = 0;
  env.start = 1;
  env.pricefile = env.parmsfile = env.rangefile = NULL;
  env.maxiter = INT_MAX;
  env.lmiter = 100000;
  env.printiter = 1;
  env.epsilon = 1e-120;
  env.threads = omp_get_max_threads();
  env.window = 0.;
  env.expfit = 0;
  env.linprice = 0;
  env.noise = 1.;
  env.diffuse = 1.;
  env.relax = 20.;
  env.seed = 635987; /* A seed based on atmospheric weather */

  parse_lppl_args(argc, argv, &env);

  /* Set up the number of threads */
  omp_set_num_threads(env.threads);

  /* Set up default parameters */
  p[0] = p[1] = 1.;
  p[2] = 1000.;
  p[3] = 1e-05;
  p[4] = -1e-05;
  p[5] = 6.28;
  p[6] = 0.;

  /* Read starting parameters */
  safe_read_parms(p, &env);
  n = umin(N, (unsigned int) floor(p[2]));
  n_int = (int) n;

  ssrand48(&mutex, env.seed);

  /* Calculate the LPPL time series */
  lppl_vector(p, x, 7, n, (void *) &env);

  /* Calculate the Brownian motion time series */
  brownian_motion(n, b, &mutex, &env);

  /* Calculate the Ornstein-Uhlenbeck time series */
  ou_noise(n, ou, &mutex, &env);

  /* Add the LPPL and noise time series */
  vector_add(n, y, x, b);
  vector_add(n, z, y, ou);

  normalize_price(n, z, &env);

  for (i = 0; i < n_int; i++)
    printf("%u\t" LG "\t" LG "\t" LG "\t" LG "\n", 
           i+1, z[i], exp(x[i]), exp(b[i]), exp(ou[i]));

  return 0;
}
