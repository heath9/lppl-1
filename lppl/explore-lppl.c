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
#include <signal.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "lm.h"
#include "lppl.h"
#include "shared.h"

#define EXPLOREMAX 1024


/** Sign of x-y */
static inline int dsign(double x, double y) {
  return (x > y) - (x < y);
}


/** Compare levmar fit information by residual error */
static int levmar_cmp(const void *x, const void *y) {
  unsigned int lm_error[2];
  double error[2];
  levmar_t *lm[2];

  lm[0] = (levmar_t *) x;
  lm[1] = (levmar_t *) y;

  lm_error[0] = lm[0]->iterations < 0;
  lm_error[1] = lm[1]->iterations < 0;
  if (lm_error[0] || lm_error[1]) 
    return dsign(lm_error[0], lm_error[1]);

  error[0] = lm[0]->info[1];
  error[1] = lm[1]->info[1];
  return dsign(error[0], error[1]);
}


/** Print all fit data
    @pre lm is not NULL */
static void print_all_fit(unsigned int k, unsigned int df, levmar_t *lm) {
  unsigned int i;

  assert(lm != NULL);

  printf("\n\n===================\n");
  for (i = 0; i < k; i++)
    lppl_print(lm[i].iterations, df, lm[i].p, lm[i].info, lm[i].covar);
}

/** Whether the group has been previously killed */
static unsigned int killed;

/* Send a signal to the children */
static void stop_children() {
  if (!get_sigrecv() || killed) 
    return;
  /* Kill the chidren */
  printf("\nSignal received, dismissing children ...\n");
  fflush(stdout);
  safe_kill(error_msg(HERE, "cannot stop explore children"), 0, SIGTERM);
  /* Avoid multiple signals */
  killed = 1;
}


/** Wait for children to complete */
static void wait_all(unsigned int k) {
  unsigned int i, wait_error;

  for (i = killed = 0; i < k; i++) {
    debug_print("wait %i ...\n", i);
    do {
      wait_error = (wait(NULL) == (pid_t) -1);
      if (wait_error) {
        perror_fail(errno != EINTR, error_msg(HERE, "cannot wait"));
        /* If control reaches here is because a (handled) signal was received 
           (a non-handled signal would presumably cause termination) */
        (void) stop_children();
      }
    } while (wait_error);
  }
}


/** Spawn analysis for each parameter set 
    @pre env is not NULL */
void spawn_levmar(unsigned int k, unsigned int n, double x[], levmar_t lm[], 
                  lppl_env_t *env) {
  unsigned int i;

  assert(env != NULL);

  for (i = 0; i < k; i++) {
    debug_print("explore %i ...\n", i);
    if (!sfork(error_msg(HERE, "cannot fork"))) {
      lm[i].iterations = levmar(lm[i].p, x, n, lm[i].info, lm[i].covar, env);
      exit(EXIT_SUCCESS);
    }
  }
}


int main(int argc, const char **argv) {
  unsigned int k;          /* Parameter set count */
  unsigned int n,          /* Price vector length */
               ny;         /* Effective price vector length */
  double p[7],             /* LPPL parameters */
         x[N],             /* Price vector */
         *y;               /* Effective price vector */
  FILE *pricefile;
  levmar_t *lm;            /* Levmar fit data */
  lppl_env_t env;

  /* Initialize error handling */
  argv0_init(argv[0]);

  /* Initialize signal action */
  set_sigrecv(0);
  safe_signal(error_msg(HERE, "cannot initialize signal handler"), 
              SIGINT, handler_sigrecv);
  safe_signal(error_msg(HERE, "cannot initialize signal handler"), 
              SIGTERM, handler_sigrecv);

  /* Allocate shared data */
  lm = (levmar_t *) 
       scmmap(error_msg(HERE, "cannot allocate shared levmar data"),
              EXPLOREMAX, sizeof(levmar_t));

  /* Set up default environment */
  env.verbose = 0;
  env.start = 1;
  env.pricefile = env.parmsfile = env.rangefile = NULL;
  env.maxiter = INT_MAX;
  env.lmiter = 1;
  env.printiter = 100000;
  env.epsilon = 1e-120;
  env.threads = 1;

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
  p[3] = 1e-05;
  p[4] = -1e-05;
  p[5] = 6.28;
  p[6] = 0.;

  /* Read starting parameters */
  k = (int) read_all_parms(p, lm, EXPLOREMAX, &env);

  /* Spawn analysis for each parameter set */
  spawn_levmar(k, ny, y, lm, &env);

  /* Wait for children to complete */
  wait_all(k);

  /* Sort fitness data by goodness of fit */
  qsort(lm, k, sizeof(levmar_t), levmar_cmp);

  /* Print fit data */
  print_all_fit(k, ny - 7, lm);

  /* Unmap shared data */
  scmunmap(error_msg(HERE, "cannot de-allocate shared levmar data"),
           (void *) lm, EXPLOREMAX, sizeof(levmar_t));

  return 0;
}
