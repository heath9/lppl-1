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
    \brief Weight functions 
    Management of weights used in levmar fitting
*/

#include <stdlib.h>
#include <assert.h>
#include "lppl.h"

/** Levmar weights 
    @pre env is not null, i<n */
double levmar_weight(unsigned int i, unsigned int n, lppl_env_t *env) {
  double window;

  assert(env != NULL);
  assert(i < n);

  window = env->window;

  return window == 0. ? 
         1. :                       /* Unweighted */
         window / (n - i + window); /* Quadratic */
}


/** Multiply a price vector by weights
    @pre x and env are not NULL */
void weigh_price(unsigned int n, double *x, lppl_env_t *env) {
  int i;     /* Loop iterator, signed as per OpenMP */
  int n_int; /* n as an int */

  assert(x != NULL);
  assert(env != NULL);

  n_int = (int) n;
  #pragma omp parallel for 
  for (i = 0; i < n_int; i++)
    x[i] *= levmar_weight(i, n_int, env);
}

