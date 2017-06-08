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

/** \brief Write a set of candidate parameters. (Experimental.) */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include "shared.h"
#include "lppl.h"

#define EXPLOREMAX 2048

int main(int argc, const char **argv) {
  unsigned int n,          /* Price vector length */
               ny;         /* Effective price vector length */
  double p[7],             /* LPPL parameters */
         x[N],             /* Price vector */
         *y;               /* Effective price vector */
         
  FILE *pricefile;
  levmar_t lm[EXPLOREMAX];
  lppl_env_t env;

  /* Initialize error handling */
  argv0_init(argv[0]);

  /* Set up default environment */
  env.verbose = 0;
  env.start = 1;
  env.pricefile = env.parmsfile = NULL;
  env.maxiter = INT_MAX;
  env.printiter = 100000;
  env.epsilon = 1e-120;
  env.threads = omp_get_max_threads();

  parse_lppl_args(argc, argv, &env);

  /* Set up the number of threads */
  omp_set_num_threads(env.threads);

  /* Read price vector */
  pricefile = safe_fopen(env.pricefile, "cannot open price file");
  n = read_price(N, x, pricefile);
  fclose(pricefile);
  debug_print("n = %u\n", n);
  assert(n >= env.start);

  /* Effective price vector */
  y  = x + env.start - 1;
  ny = n - env.start + 1;
  vector_f(ny, y, log);

  /* Set up default parameters */
  p[0] = p[1] = 1.;
  p[2] = 1000.;
  p[3] = 1.;
  p[4] = 0.;
  p[5] = 6.28;
  p[6] = 0.;

  /* Read starting parameters */
  /* k = */
  read_all_parms(p, lm, EXPLOREMAX, &env);

  /*
  pi_over_8 = M_PI / 8.;
  for (j = 0; j < k; j++) {
    print_lppl_parms(lm[j].p);
    printf("E=999.999\n");
    state = linear_lsq(lm[j].p, y, ny, &env);
    if (!state && lm[j].p[1] > 0 && lm[j].p[2] > n) {
      print_lppl_parms(lm[j].p);
      printf("E=999.999\n");
    }
    for (lm[j].p[6] = 0.; lm[j].p[6] < M_PI; lm[j].p[6] += pi_over_8) {
      state = linear_lsq(lm[j].p, y, ny, &env);
      if (!state && lm[j].p[1] > 0 && lm[j].p[2] > n) {
        print_lppl_parms(lm[j].p);
        printf("E=999.999\n");
      }
    }
  }
  */

  /*
  for (p[2] = n + .5; p[2] < n + n; p[2]++) {
    vector_copy(7, q, p);
    state = linear_lsq(q, y, ny, &env);
    if (!state && q[1] > 0) {
      print_lppl_parms(q);
      printf("E=999.999\n");
    }
  }
  */

  lppl_guess(lm[0].p, y, 1.76, 220, 270, ny, &env);
  /* print_lppl_parms(lm[0].p); */
  for (lm[0].p[6] = 0.; lm[0].p[6] < M_PI; lm[0].p[6] += M_PI / 8.) {
    print_lppl_parms(lm[0].p);
    printf("E=999.999\n");
  }
  /*
  xmin = 147;
  xmax = 167;
  for (i = 0; i < 7; i++) {
    p[3] = m[i];
    for (rho = 1.15; rho < 2.; rho += .01) {
      lppl_guess(p, y, rho, xmin, xmax, ny, &env);
      if (p[2] > n) {
        print_lppl_parms(p);
        printf("E=999.999\n");
      }
    }
  }

  xmin = 125;
  xmax = 147;
  for (i = 0; i < 7; i++) {
    p[3] = m[i];
    for (rho = 1.15; rho < 2.; rho += .01) {
      lppl_guess(p, y, rho, xmin, xmax, ny, &env);
      if (p[2] > n) {
        print_lppl_parms(p);
        printf("E=999.999\n");
      }
    }
  }
  */

  return 0;
}
