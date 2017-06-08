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

/** \brief Iterations management in levmar */

#include <assert.h>
#include "shared.h"
#include "lppl.h"

/** Initialize a levmar_iter structure 
    @pre lmiter is not NULL */
void iter_init(int init_iterations, levmar_iter_t *lmiter) {
  assert(lmiter != NULL);

  lmiter->adaptive = init_iterations < 0;
  lmiter->start_up = 1;

  if (init_iterations > 0)
    lmiter->iterations = init_iterations;
  else if (init_iterations < 0)
    lmiter->iterations = -init_iterations;
  else 
    lmiter->iterations = 1;
}


/** Adjust the number of levmar iterations during the start-up phase
    @pre lmiter is not NULL */
static void iter_update_startup(unsigned int levmar_better,
                                levmar_iter_t *lmiter) {
  assert(lmiter != NULL);

  if (levmar_better) {
    lmiter->iterations += lmiter->iterations;
  } else {
    lmiter->start_up = 0;
    (lmiter->iterations)--;
  }
}


/** Adjust the number of levmar iterations at regime
    @pre lmiter is not NULL */
static void iter_update_regime(unsigned int levmar_better,
                               levmar_iter_t *lmiter) {
  assert(lmiter != NULL);

  if (levmar_better)
    (lmiter->iterations)++;
  else
    (lmiter->iterations)--;
}


/** Determine whether there is enough information to compute the marginal
    error reduction per unit of time */
static unsigned int marginal_info(
              double levmar_stat_time[],/**< Last 2 levmar times */
              int   levmar_stat_iter[] /**< Last 2 levmar iterations */) {
  return levmar_stat_iter[0] > 0 && 
         levmar_stat_iter[1] != levmar_stat_iter[0] && 
         levmar_stat_time[1] != levmar_stat_time[0];
}


/** Update the statistics on levmar iteration time 
    @return levmar error reduction per unit of time (marginal) */
static double update_iter_stat(
              double levmar_stat_time[], /**< Last 2 levmar times */
              int    levmar_stat_iter[], /**< Last 2 levmar iterations */
              double levmar_time,        /**< New levmar time */
              int levmar_iter,           /**< New levmar iterations */
              double levmar_err          /**< Error reduction */) {

  /* Update rolling statistics */
  if (levmar_iter != levmar_stat_iter[1]) {
    levmar_stat_time[0] = levmar_stat_time[1];
    levmar_stat_iter[0] = levmar_stat_iter[1];
    levmar_stat_iter[1] = levmar_iter;
  }
  levmar_stat_time[1] = levmar_time;

  debug_print("Levmar time (current)  = " LG "\n", levmar_stat_time[1]);
  debug_print("Levmar time (previous) = " LG "\n", levmar_stat_time[0]);
  debug_print("Levmar iter (current)  = " LG "\n", levmar_stat_iter[1]);
  debug_print("Levmar iter (previous) = " LG "\n", levmar_stat_iter[0]);

  /* Calculate marginal unit error reduction */
  return marginal_info(levmar_stat_time, levmar_stat_iter) ?
         levmar_err * (levmar_stat_iter[1] - levmar_stat_iter[0]) / 
           (levmar_stat_iter[1] * (levmar_stat_time[1] - levmar_stat_time[0])) 
         :
         levmar_err / levmar_time;
}


/** Adjust the number of levmar iterations 
    @pre lmiter is not NULL */
void iter_update(double levmar_time, double levmar_err, int levmar_iter,
                 double lin_time,    double lin_err,
                 double levmar_stat_time[], /**< Last 2 levmar times */
                 int    levmar_stat_iter[], /**< Last 2 levmar iterations */
                 levmar_iter_t *lmiter) {
  unsigned int levmar_better;
  double levmar_gain, lin_gain;

  assert(lmiter != NULL);

  if (!lmiter->adaptive) return;

  levmar_gain = update_iter_stat(levmar_stat_time, levmar_stat_iter, 
                                 levmar_time, levmar_iter, levmar_err);
  debug_print("Levmar marginal gain: " LG "\n", levmar_gain);
  lin_gain    = lin_err / lin_time;
  debug_print("Linear marginal gain: " LG "\n", levmar_gain);
  levmar_better = levmar_gain >= lin_gain;

  if (lmiter->start_up) 
    iter_update_startup(levmar_better, lmiter);
  else
    iter_update_regime(levmar_better, lmiter);

  if (!lmiter->iterations)
    lmiter->iterations = 1;

  debug_print("Iterations: %i\n", lmiter->iterations);
}
