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
  \brief Parse command line parameters 
*/

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <popt.h>
#include "shared.h"
#include "lppl.h"


/** If verbose output is desired, print command line arguments 
    @pre env is not NULL */
static void print_lppl_args(lppl_env_t *env) {
  assert(env != NULL);

  if (env->verbose) {
    printf("Price file:        %s\n",     env->pricefile);
    printf("Parameter file:    %s\n",     env->parmsfile);
    printf("Range file:        %s\n",     env->rangefile);
    printf("Exp fit:           %u\n",     env->expfit);
    printf("Linear prices:     %u\n",     env->linprice);
    printf("Weight window:     " LG "\n", env->window);
    printf("Max iterations:    %i\n",     env->maxiter);
    printf("LM iterations:     %i\n",     env->lmiter);
    printf("Print iterations:  %i\n",     env->printiter);
    printf("Epsilon:           " LG "\n", env->epsilon);
    printf("Start fit:         %u\n",     env->start);
    printf("OpenMP threads:    %u\n",     env->threads);
    printf("Interpolation:     %u\n",     env->f == step_lppl);
    printf("Noise power:       " LG "\n", env->noise);
    printf("Diffusion (var):   " LG "\n", env->diffuse);
    printf("Relaxation (tau):  " LG "\n", env->relax);
    printf("Seed:              %u\n",     env->seed);
    printf("Penny round:       %u\n",     !env->nopenny);
  }
}


void parse_lppl_args(int argc, const char **argv, lppl_env_t *env) {
  char c;
  int func;
  poptContext opt_con;
  struct poptOption options_table[] =  {
    { "pricefile", 'f', POPT_ARG_STRING, &(env->pricefile), 0,
      "Price file name", NULL},
    { "parmsfile", 'p', POPT_ARG_STRING, &(env->parmsfile), 0,
      "Parameter file name", NULL},
    { "rangefile", 'r', POPT_ARG_STRING, &(env->rangefile), 0,
      "Range file name (explore only)", NULL},
    { "maxiter", 'i', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->maxiter), 0, "Maximum number of iterations", NULL},
    { "lmiter", 'L', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->lmiter), 0, 
      "Maximum number of iterations in a LM computation block", 
      NULL},
    { "printiter", 'P', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->printiter), 0, 
      "Number of LM iteration blocks between print-outs (0: no printouts)", 
      NULL},
    { "linprice", 'l', POPT_ARG_NONE, &(env->linprice), 0, 
      "Do not take the logarithm of prices", NULL},
    { "expfit", 'E', POPT_ARG_NONE, &(env->expfit), 0, 
      "Change A, B to start from exponential fit", NULL},
    { "window", 'w', POPT_ARG_DOUBLE|POPT_ARGFLAG_SHOW_DEFAULT,
      &(env->window), 0, "Weight window (0 is unweighted)", NULL},
    { "epsilon", 'e', POPT_ARG_DOUBLE|POPT_ARGFLAG_SHOW_DEFAULT,
      &(env->epsilon), 0, "Error at termination", NULL},
    { "start", 'S', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->start), 0, "Initial time step", NULL},
    { "threads", 'T', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->threads), 0, "Threads", NULL},
    { "func", 'F', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &func, 0, "Interpolation function (0: lppl, 1: step) disabled", NULL},
    { "noise", 'N', POPT_ARG_DOUBLE|POPT_ARGFLAG_SHOW_DEFAULT,
      &(env->noise), 0., 
      "Noise power, gen only", NULL},
    { "diffuse", '\0', POPT_ARG_DOUBLE|POPT_ARGFLAG_SHOW_DEFAULT,
      &(env->diffuse), 0., 
      "Diffusion constant (variance) in Ornstein-Uhlenbeck, gen only", NULL},
    { "relax", '\0', POPT_ARG_DOUBLE|POPT_ARGFLAG_SHOW_DEFAULT,
      &(env->relax), 0., 
      "Relaxation constant in Ornstein-Uhlenbeck, gen only", NULL},
    { "seed", '\0', POPT_ARG_INT|POPT_ARGFLAG_SHOW_DEFAULT, 
      &(env->seed), 0, "Random seed (repeatable sequence require -T 1)", NULL},
    { "nopenny", '\0', POPT_ARG_NONE, &(env->nopenny), 0, 
      "Do not round prices to the closest penny (gen only)", NULL},
    { "verbose", 'v', POPT_ARG_NONE, &(env->verbose), 0, 
      "Set verbose output", NULL},
    POPT_AUTOHELP
    POPT_TABLEEND
  };

  /* Initialize defaults */
  func = 0;

  opt_con = poptGetContext(NULL, argc, argv, options_table, 0);

  while ((c  = poptGetNextOpt(opt_con)) >= 0) {
    switch (c) {
    default:
      fprintf(stderr, "%s: strange error %c in command line parsing\n",
              argv[0], c);
      exit(EXIT_FAILURE);
    }
  }

  /* Post-process arguments */
  env->f = func ? step_lppl : lppl;
  if (env->linprice)
    env->expfit = 0;

  print_lppl_args(env);
}
