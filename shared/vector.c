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

/** \brief Vector operations */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "shared.h"

typedef struct {
  double sum;     /**< Current sum value */
  unsigned int i; /**< Next array element to be added */
} partial_cursor_t;

/** Computing sums of many relatively small numbers: recursion 
    @pre cursor is not NULL*/
static void partial_sum(unsigned int n, const double y[], 
                        partial_cursor_t *cursor) {
  partial_cursor_t local_cursor;

  assert(cursor != NULL);

  if (cursor->i >= n) return;

  local_cursor.sum = y[cursor->i];
  local_cursor.i   = cursor->i + 1;

  while (local_cursor.i < n && local_cursor.sum <= cursor->sum) 
    partial_sum(n, y, &local_cursor);
  
  cursor->sum += local_cursor.sum;
  cursor->i    = local_cursor.i;
}


/** Compute the sum of many relatively small numbers: front end */
double stable_sum(unsigned int n, const double y[]) {
  partial_cursor_t cursor;

  cursor.sum = 0.;
  cursor.i   = 0;
  while (cursor.i < n) 
    partial_sum(n, y, &cursor);

  return cursor.sum;
}


/** Dot product of a and b (vectors of size n)
    @pre a, b are not NULL */
double dot_product(unsigned int n, double *a, double *b) {
  int i, bound;
  double result;

  assert(a != NULL);
  assert(b != NULL);

  bound = (int) n;
  result = 0.;
  for (i = 0; i < bound; i++)
    result += a[i] * b[i];

  return result;
}


/** Downshift the contents of a vector of size n 
    @pre a is not NULL */
void down_shift(unsigned int n, double *a) {
  unsigned int i;

  assert(a != NULL);

  n--;
  for (i = 0; i < n; i++)
    a[i] = a[i+1];
}


/** Copy a vector into another vector 
    @pre a, b are not NULL */
void vector_copy(unsigned int n, double *dst, const double *src) {
  int i, bound;

  assert(src != NULL);
  assert(dst != NULL);

  bound = (int) n;
  #pragma omp parallel for
  for (i = 0; i < bound; i++)
    dst[i] = src[i];
}


/** Replace each vector elements x_i with f(x_i) */
void vector_f(unsigned int n, double x[], double (*f)(double)) {
  int i, bound;

  bound = (int) n;
  #pragma omp parallel for
  for (i = 0; i < bound; i++)
    x[i] = f(x[i]);
}


/** Compute the L_2 norm of x-y */
double vector_error(unsigned int n, double x[], double y[]) {
  int i, nx;
  double *d, err;

  nx = (int) n;
  d = (double *) scalloc((size_t ) n, sizeof(double),
                         error_msg(HERE, "cannot allocate difference vector"));

  #pragma omp parallel for
  for (i = 0; i < nx; i++) {
    d[i] = x[i] - y[i];
    d[i] *= d[i];
  }
    
  err = stable_sum(n, d);
  free(d);

  return err;
}


/** Add two vectors */
void vector_add(unsigned int n, /**< Size of a, b */
                double c[],     /**< Result c=a+b */
                double a[],     
                double b[]) {
  int i,     /* Loop counter */
      n_int; /* n, as an int */

  n_int = (int) n;

  /* Additive noise */
  #pragma omp parallel for
  for (i = 0; i < n_int; i++)
    c[i] = a[i] + b[i];
}

