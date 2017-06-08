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

/** \brief Fit printouts */

#include <stdio.h>
#include <math.h>
#include "lm.h"
#include "shared.h"
#include "lppl.h"

/** Print the resulting value after fit  */
void lppl_print(int iterations, unsigned int df, /**< Degrees of freedom */
                double p[], double info[], double covar[][7]) {
  unsigned int i, j;
  double var_prod;

  printf("\nIterations %i\n", iterations);

  print_lppl_parms(p);

  printf("Average error: " LG "\n", info[1] / (double) df);

  for (i = 0; i < LM_INFO_SZ; i++) 
    printf("info[%u]: " LG "\n", i, info[i]);

  printf("Variance\n");
  for (i = 0; i < 7; i++) 
    printf("p[%u]: " LG "\n", i, covar[i][i]);

  printf("Correlation\n");
  for (i = 0; i < 7; i++) {
    for (j = 0; j <= i; j++) {
      var_prod = covar[i][i] * covar[j][j];
      if (var_prod > 0.)
        printf("%.3g\t", covar[i][j] / sqrt(var_prod));
      else if (var_prod < 0.)
        printf("im\t");
      else
        printf("inf\t");
    }
    printf("\n");
  }
  fflush(stdout);
}
