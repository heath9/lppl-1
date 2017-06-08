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

/** \brief Generate noise processes

    Generate noise processes:
    - Brownian motion
    - Ornstein-Uhlbeck
**/
    
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "lppl.h"
#include "shared.h"

/** Compute a Brownian motion time series 
    @pre mutex, env are not NULL */
void brownian_motion(
              unsigned int n,         /**< Size of b */
              double b[],             /**< Brownian motion time series */
              pthread_mutex_t *mutex, /**< Mutex for drand */
              lppl_env_t *env) {
  int i,     /* Loop counter */
      n_int; /* n, as an int */
  double noise_power,
         rnd_tmp[2];

  assert(mutex != NULL);
  assert(env != NULL);

  generate_normal(rnd_tmp, mutex);

  noise_power = env->noise;

  if (n % 2) {
    /* Since random numbers are generated in pairs, if n is odd,
       the last random number requires special handling */
    b[n-1] = noise_power * rnd_tmp[0];
    n--;
  }
  assert(n % 2 == 0);
  n_int = (int) n;

  /* Generate normal increments */
  #pragma omp parallel for
  for (i = 0; i < n_int; i += 2) {
    generate_normal(b+i, mutex);
    b[i]   *= noise_power;
    b[i+1] *= noise_power;
    debug_print(LG "\n" LG "\n", b[i], b[i+1]);
  }

  /* Warm-up the Brownian motion process: commented out 
  b[0] = sqrt((double) (n + n)) * noise_power * rnd_tmp[1]; */

  /* Generate Brownian motion */
  for (i = 1; i < n_int; i++) {
    b[i] += b[i-1];
    debug_print("B[%u]=" LG "\n", i, b[i]);
  }
}


/** Generate Ornstein-Uhlenbeck noise.
    If env->relax=0, no noise is generated.
    For the algorithm, see Gillespie - Physical review E, 1996
    @pre mutex, env are not NULL */
void ou_noise(unsigned int n,         /**< Size of ou */
              double ou[],            /**< OU time series */
              pthread_mutex_t *mutex, /**< Mutex for drand */
              lppl_env_t *env) {
  unsigned int i,       /* Loop counter */
               n_bound; /* Bound on loop counter */
  double alpha, beta,   /* OU calculation constants */
         rnd_tmp[2];    /* Temporary storage for random variables */

  assert(mutex != NULL);
  assert(env != NULL);

  if (env->relax == 0.)
    return;

  /* Calculate the value of the constant in the OU process */
  alpha = exp(-1. / env->relax);
  beta  = sqrt(.5 * env->diffuse * env->relax * (1. - alpha * alpha));

  /* Initialize the process with the long term probability distribution */
  generate_normal(rnd_tmp, mutex);
  ou[0] = sqrt(.5 * env->diffuse * env->relax) * rnd_tmp[0];

  /* Simulate the process for pairs of consecutive indexes 
     (normal random variables are generated in pairs) */
  n_bound = n - 1;
  for (i = 1; i < n_bound; i += 2) {
    generate_normal(rnd_tmp, mutex);
    ou[i]   = alpha * ou[i-1] + beta * rnd_tmp[0];
    ou[i+1] = alpha * ou[i]   + beta * rnd_tmp[1];
    debug_print(LG "\n" LG "\n", ou[i], ou[i+1]);
  }

  /* If n is even, then the last random value needs to be generated */
  if (i < n) {
    assert(!(n % 2));
    assert(i == n - 1);
    generate_normal(rnd_tmp, mutex);
    ou[i]   = alpha * ou[i-1] + beta * rnd_tmp[0];
    debug_print(LG "\n", ou[i]);
  }
}
