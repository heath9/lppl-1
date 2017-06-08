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

/** \brief Read and print parameters */

#include <stdio.h>
#include <assert.h>
#include "shared.h"
#include "lppl.h"


/** Print parameters in gnuplot .par format */
void print_lppl_parms(double p[]) {
  printf("A=" LG "\n", p[0]);
  printf("B=" LG "\n", p[1]);
  printf("D=" LG "\n", p[1]);
  printf("T=" LG "\n", p[2]);
  printf("m=" LG "\n", p[3]);
  printf("C=" LG "\n", p[4]);
  printf("w=" LG "\n", p[5]);
  printf("p=" LG "\n", p[6]);
  printf("a=" LG "\n", p[4]);
  printf("r=" LG "\n", p[5]);
}


/** Read parameters in gnuplot .par format: %c=<double>  until end-of-file
    or untile a line with E=<double> 
    @return the number of parameter sets read (1 or 0)
*/
unsigned int read_lppl_parms(FILE *parmsfile, double p[], 
                             unsigned int verbose) {
  char c;
  unsigned int mid_of_list, parms_set;
  double tmpvar;

  parms_set = 0;

  if (parmsfile == NULL) return parms_set;

  mid_of_list = 1;
  while (mid_of_list && fscanf(parmsfile, "%c=%lg\n", &c, &tmpvar) != EOF) {
    switch (c) {
    case 'A':
      p[0] = tmpvar;
      parms_set = 1;
      break;
    case 'T':
      p[2] = tmpvar;
      parms_set = 1;
      break;
    case 'm':
      p[3] = tmpvar;
      parms_set = 1;
      break;
    /* LPPL parameters */
    case 'B':
      p[1] = tmpvar;
      parms_set = 1;
      break;
    case 'C':
      p[4] = tmpvar;
      parms_set = 1;
      break;
    case 'w':
      p[5] = tmpvar;
      parms_set = 1;
      break;
    case 'p':
      p[6] = tmpvar;
      parms_set = 1;
      break;
    /* Step LPPL parameters */
    case 'D':
      p[1] = tmpvar;
      parms_set = 1;
      break;
    case 'a':
      p[4] = tmpvar;
      parms_set = 1;
      break;
    case 'r':
      p[5] = tmpvar;
      parms_set = 1;
      break;
    /* End of paramter list */
    case 'E':
      mid_of_list = 0;
      break;
    default:
      fprintf(stderr, "lppl: unrecognized parameter %c\n", c);
    }
  }

  if (verbose)
    print_lppl_parms(p);

  return parms_set;
}

/** Read starting parameters, and output a warning if no parameter set 
    could be detected.
    @pre env is not NULL */
void safe_read_parms(double p[], lppl_env_t *env) {
  FILE *parmsfile;

  assert(env != NULL);

  if (env->parmsfile == NULL) return;

  parmsfile = safe_fopen(env->parmsfile, "cannot open parameter file");
  if (!read_lppl_parms(parmsfile, p, env->verbose)) 
    fprintf(stderr, "%s\n", 
            error_msg(HERE, "no parameter set in parameter file"));
  fclose(parmsfile);
}


/** Read a set of starting parameters 
    @pre env is not NULL */
unsigned int read_all_parms(
             double p[],              /**< Default parameter value */
             levmar_t lm[],           /**< Parameter set */
             unsigned int exploremax, /**< Maximum number of parameter sets */
             lppl_env_t *env) {
  unsigned int i, state;
  FILE *parmsfile;

  assert(env != NULL);

  i = 0;
  /* No parameter file */
  if (env->parmsfile == NULL) return i;

  /* Read parameter sets from parameter file */
  parmsfile = safe_fopen(env->parmsfile, "cannot open parameter file");
  state = 1;
  do {
    /* Set parameter vector to default values */
    vector_copy(7, lm[i].p, p);
    /* Read parameters */
    state = read_lppl_parms(parmsfile, lm[i].p, env->verbose);
    i += state;
  } while (i < exploremax && state);
  fclose(parmsfile);

  debug_print("Parameter sets: %u\n", i);
  if (i == exploremax)
    fprintf(stderr, "%s\n", 
            error_msg(HERE, "warning: too many parameter sets"));

  return i;
}

