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
    \brief Iterative solution of least squares with levmar 
*/

#include <assert.h>
#include <float.h>
#include <sys/time.h>
#include "lm.h"
#include "shared.h"
#include "lppl.h"


/** Redefine the initial value of \mu */
#undef LM_INIT_MU
#define LM_INIT_MU 100.
#define LM_MAX_MU 1000000.


/** Returns the error code after an execution of levmar, or 0 if there
    was no error */
static int levmar_error(int state, double info[]) {
  /* levmar failed */
  if (state < 0) 
    return -1;

  /* levmar worsened the error */
  if (info[1] > info[0])
    return -2;

  return 0;
}


/** Determine whether an increase in mu was requested and accomplished */
static inline unsigned int increase_mu(double opts[], double info[]) {
  return (info[6] == 4 || info[6] == 5) && opts[0] > 0.;
}


/** Determine whether levmar has made progress and progress is still possible 
    in future iterations */
static inline unsigned int levmar_progress(double opts[], double info[]) {
  return (info[1] < info[0]) && 
         (info[6] == 3 || increase_mu(opts, info));
}


/** Determine if either levmar or linear least square are making progress 
    or future progress is still possible */
static inline unsigned int linlev_progress(double opts[], double info[]) {
  return info[6] == 3.5 || levmar_progress(opts, info);
}


/** Establish whether the iterative algorithm should continue */
static inline unsigned int levmar_continue(
                       double opts[], double info[],
                       int iterations, int maxiter) {
  debug_print("Iterations %i\n", iterations);
  return linlev_progress(opts, info) &&
         iterations < maxiter     && /* Iteration count within bound */
         !get_sigrecv();             /* No signal has been received */
}


/** @return initial value of mu upon restarting, or -1 in case of error
    (an increase of mu was required but not carried out)
    @pre init_mu is not NULL */
static inline double restart_mu(
                     double termination_reason, 
                     double stop_mu, /**< Value of mu at termination */
                     double *init_mu /**< Initial value of mu */) {
  double new_mu;

  assert(init_mu != NULL);

  if (termination_reason == 4. || /* Singular matrix. 
                                     Restart from current p with increased mu */
      termination_reason == 5. || /* No further error reduction is possible. 
                                     Restart with increased mu */
      stop_mu == 0.               /* No gradient component */ ) {
    *init_mu += *init_mu;
    if (*init_mu <= LM_MAX_MU) {
      new_mu = *init_mu;
    } else {
      *init_mu = LM_MAX_MU;
      new_mu = -1.;
    }
  } else
    new_mu = stop_mu;

  return new_mu;
}


/** Initialize parameter upper and lower bounds */
static inline void init_bounds(double lb[], double ub[], unsigned int n) {
  unsigned int i;

  for (i = 0; i < 7; i++) {
    lb[i] = -DBL_MAX;
    ub[i] =  DBL_MAX;
  }
  lb[1] = 0.;
  lb[2] = n;
  lb[3] = 0.;
  ub[3] = 1.;
  lb[5] = 0.;
  /* Experimental: no bounds on phi
  lb[6] = 0.;
  ub[6] = M_PI; */
}


/** Periodically prints the current parameter vector 
    @pre printcnt is not NULL 
    @return 1 if something was printed, 0 otherwise */
static unsigned int periodic_lppl_print(int printiter, int *printcnt, 
                                        int iterations,  int df, double p[], 
                                        double info[], double covar[][7],
                                        int itmax, double levmar_time) {
  unsigned int state;

  assert(printcnt != NULL);

  state = 0;

  if (!printiter) /* Intermediate printouts have been disabled */
    return state;

  if (printiter == ++(*printcnt)) {
    lppl_print(iterations, df, p, info, covar);
    printf("Levmar iterations: %i\n", itmax);
    printf("Levmar time: " LG "\n", levmar_time);
    *printcnt = 0;
    state = 1;
  }

  return state;
}

/** Periodically print the timing of the linear system solution */
static void periodic_linear_print(unsigned int print_time, 
                                  double linear_time) {
  if (print_time)
    printf("Linear time: " LG "\n", linear_time);
}



/** Attempt to improve the current parameters by solving a least-square 
    problem in which the non-linear components are assumed to be constant
    @pre p, x, env are not NULL 
    @return the error improvement due to the linear subsystem solution */
static double linear_lppl(double *p, double *x, unsigned int n, double info[],
                          lppl_env_t *env) {
  int lsq_state;
  double lp[7], new_error;

  assert(p != NULL);
  assert(x != NULL);
  assert(env != NULL);

  /* Calculate a new vector of parameters lp by solving the linear subsystem */
  vector_copy(7, lp, p);
  lsq_state = linear_lsq(lp, x, n, env);
  if (lsq_state != 0) return 0.; 

  /* Compute the error of the new solution */
  new_error = solution_error(n, x, lp, env);

  if (new_error < info[1]) {  
    /* There is an improvement */
    debug_print("Linear LPPL improves error: " LG " to " LG "\n",
                info[1], new_error);
    vector_copy(7, p, lp);
    /* Assign termination code.
       Note: info[6] = 3.5 is a new code: fit terminated but progress can
             still be made if non-linear components are kept constant */
    info[6] = new_error <= env->epsilon ? 6 : 3.5; 
  }

  return info[1] - new_error;
}


/** Levmar: see http://www.ics.forth.gr/~lourakis/levmar/ for an explanation 
    of parameters. Here is a list of changes:
    - info[6]=3.5 stopped by using a non-linear method
    @pre p, x, env are not NULL 
    @return the number of levmar iterations if all went well, otherwise:
            - -1 if levmar failed
            - -2 if a levmar iteration increased the solution error */
int levmar(double *p, double *x, unsigned int n, double info[], 
           double covar[][7], lppl_env_t *env) {
  int state, iterations, maxiter, printiter, printcnt, 
      err,                 /* Error code */
      levmar_stat_iter[2]; /* Statistics on levmar iterations */
  unsigned int print_time, /* Whether print should occur at this iteration */
               df;         /* Degrees of freedom */
  double lb[7], ub[7],     /* LPPL parameter lower and upper bounds */
         opts[LM_INFO_SZ], /* Fit options */
         *work,            /* Work area */
         init_mu,          /* Starting value of mu */
         linear_error,     /* Error improvement due to the linear solution */
         levmar_time, linear_time, /* Time taken by levmar, linear */
         levmar_stat_time[2];       /* Statistics on levmar time */
  struct timeval current_time;
  levmar_iter_t lmiter;

  assert(p != NULL);
  assert(x != NULL);
  assert(env != NULL);

  weigh_price(n, x, env);

  /* Set up temporary work area */
  work = (double *) scalloc(LM_DER_WORKSZ(7, n), sizeof(double), 
                            error_msg(HERE, "cannot allocate work area"));

  /* Parameter boundaries */
  init_bounds(lb, ub, n + env->start);
  
  /* Default optimization options */
  opts[0] = init_mu = LM_INIT_MU;
  opts[1] = opts[2] = opts[3] = env->epsilon;
  opts[4] = -LM_DIFF_DELTA; /* Useful only for the old dlevmar_dif */

  /* Local alias of environment parameters */
  maxiter   = env->maxiter;
  printiter = env->printiter;

  iter_init(env->lmiter, &lmiter);

  iterations = printcnt = 0; 
  df = n - 7;
  perror_fail(gettimeofday(&current_time, NULL) == -1, 
              error_msg(HERE, "cannot find initial time")); 
  printf("\nIterations 0\n");
  printf("Average error: " LG "\n", 
         solution_error(n, x, p, env) / (double) df);

  do {
    /* Approximate the Jacobian by finite differences: old code, no work vector
    iterations += dlevmar_bc_dif(lppl_vector, p, x+env.start-1, 7, 
                                 n-env.start+1, 
                                 lb, ub, env.printiter, opts, info, NULL, 
                                 &(covar[0][0]), (void *) &env);

    lppl_print(iterations, p, info, covar); */

    /* Execute lmiter iterations of the Levenberg-Marquardt algorithm */
    state = dlevmar_bc_der(lppl_vector, jac_lppl_vector, p, x, 7, n, 
                           lb, ub, lmiter.iterations, opts, info, work, 
                           &(covar[0][0]), (void *) env);

    err = levmar_error(state, info);
    if (err < 0) 
      return err;
    iterations += state;

    levmar_time = time_gap(error_msg(HERE, "cannot find levmar time"), 
                           &current_time);

    print_time = periodic_lppl_print(printiter, &printcnt, iterations, df, 
                                     p, info, covar, 
                                     lmiter.iterations, levmar_time);

    /* Attempt to improve solution by fixing the non-linear parameters */
    linear_error = linear_lppl(p, x, n, info, env);

    linear_time = time_gap(error_msg(HERE, "cannot find linear time"),
                           &current_time);

    periodic_linear_print(print_time, linear_time);
    fflush(stdout);

    iter_update(levmar_time, info[0] - info[1], state, 
                linear_time, linear_error, 
                levmar_stat_time, levmar_stat_iter, &lmiter);

    /* Reset the value of mu
       Note: force the old value of mu even after jumping to the
       solution of the linear LPPL */
    opts[0] = restart_mu(info[6], info[4], &init_mu);

  } while (levmar_continue(opts, info, iterations, maxiter));

  free(work);

  normalize_lppl(p);
  lppl_print(iterations, n-7, p, info, covar);

  return iterations;
}
