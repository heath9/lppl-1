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

/** \brief Main for lppl */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <omp.h>
#include "lm.h"
#include "lppl.h"
#include "shared.h"


/** Convert prices to logarithmic form 
    @pre env is not NULL */
static void price2log(unsigned int n, double x[], lppl_env_t *env) {
  assert(env != NULL);

  if (!env->linprice) {
    debug_print("Converting prices to log\n");
    vector_f(n, x, log);
  }
}


int main(int argc, const char **argv) {
  unsigned int n,          /* Price vector length */
               ny;         /* Effective price vector length */
  double p[7],             /* LPPL parameters */
         covar[7][7],      /* Covariance matrix */
         x[N],             /* Price vector */
         *y,               /* Effective price vector */
         info[LM_INFO_SZ]; /* Fit information */
  FILE *pricefile;
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

  parse_lppl_args(argc, argv, &env);

  /* Read price vector */
  pricefile = safe_fopen(env.pricefile, "cannot open price file");
  n = read_price(N, x, pricefile);
  fclose(pricefile);
  debug_print("n = %u\n", n);
  assert(n >= env.start);
  price2log(n, x, &env);

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
  y  = x + env.start - 1;
  ny = n - env.start + 1;

  /* Start from an exponential fit */
  safe_exp_lsq(p, y, ny, &env);

  return levmar(p, y, ny, info, covar, &env) >= 0 ?
         0 : -1;
}
