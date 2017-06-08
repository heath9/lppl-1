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

/** Generate noise as the superposition of Brownian motion with an
    Ornstein-Uhlenbeck process.

    The program generates pure noise:
    - LPPL is not underlying
    - Prices are not exponentiated
    - Prices are not rounded
*/

#include <stdio.h>
#include <limits.h>
#include <omp.h>
#include "shared.h"
#include "lppl.h"

int main(int argc, const char **argv) {
  unsigned int n, /* Sequence length */
               i; /* Loop counter */
  double p[7], /* LPPL parameters */
         *b,   /* Brownian motion */
         *ou,  /* Ornstein-Uhlenbeck process */
         *y;   /* Noise process */
  pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; /* drand mutex */
  lppl_env_t env;

  /* Initialize error handling */
  argv0_init(argv[0]);

  /* Set up default environment */
  env.verbose = 0;
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

  /* Set up default parameters */
  p[2] = 1000.;

  parse_lppl_args(argc, argv, &env);
    
  /* Read starting parameters */
  safe_read_parms(p, &env);
  n = p[2];

  ssrand48(&mutex, env.seed);

  b  = (double *) scalloc(error_msg(HERE, 
                                    "cannot allocate Brownian motion vector"), 
                          (size_t) n, sizeof(double));
  ou = (double *) scalloc(error_msg(HERE, "cannot allocate OU vector"),
                          (size_t) n, sizeof(double));
  y  = (double *) scalloc(error_msg(HERE, "cannot allocate noise vector"),
                          (size_t) n, sizeof(double));

  /* Calculate the Brownian motion time series */
  brownian_motion(n, b, &mutex, &env);
  
  /* Calculate the Ornstein-Uhlenbeck time series */
  ou_noise(n, ou, &mutex, &env);
  
  /* Add the LPPL and noise time series */
  vector_add(n, y, ou, b);

  for (i = 0; i < n; i++) 
    printf("%u\t" LG "\t" LG "\t " LG "\n", i+1, y[i], b[i], ou[i]);

  free(b);
  free(ou);
  free(y);

  return 0;
}
