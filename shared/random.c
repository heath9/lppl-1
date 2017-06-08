/* Copyright (c) 2004, 2005, 2009
 * Vincenzo Liberatore
 * Case Western Reserve University
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

/** \brief Random process simulation */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <errno.h>

/** Initialize the rand48 random number generator safely:
    - Use a mutex to protect from multiple concurrent invocations
    - Extract and discard the first random number 
    @pre mutex is not NULL */
void ssrand48(pthread_mutex_t *mutex, unsigned long seed) {
  assert(mutex != NULL);

  pthread_mutex_lock(mutex);

  srand48(seed);

  /* Some random number generators have something funny on their first
     values, such as returning 0, so it is better to clean the pipe */
  (void) drand48();

  pthread_mutex_unlock(mutex);
}


/** Thread-safe drand48 
    @pre mutex is not NULL */
double sdrand48(pthread_mutex_t *mutex) {
  double u;

  assert(mutex != NULL);

  pthread_mutex_lock(mutex);
  u = drand48();
  pthread_mutex_unlock(mutex);

  return u;
}


/* Generate the log of a random number in (0,1) */
static double log_rand(pthread_mutex_t *mutex) {
  double u;

  /* Keep computing the log of a random number until no error is detected */
  do {
    errno = 0;
    u     = log(sdrand48(mutex));
  } while (errno != 0);

  return u;
}


/** Generate a pair of Gaussian random variable with Box-Muller 
    @pre mutex is not NULL */
void generate_normal(double g[], pthread_mutex_t *mutex) {
  double u1, u2, two_pi;

  assert(mutex != NULL);

  u1     = -log_rand(mutex);
  u1    +=  u1;
  u1     = sqrt(u1);
  u2     = sdrand48(mutex);
  two_pi = M_PI + M_PI;
  g[0]   = u1 * cos(two_pi * u2);
  g[1]   = u1 * sin(two_pi * u2);

  assert(finite(g[0]));
  assert(finite(g[1]));
}
