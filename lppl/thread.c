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
  \brief Function to start, join threads and parallelize loops
*/

#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "lppl.h"

/** Safe thread creation 
    @pre thread, arg are not NULL*/
void safe_pthread_create(pthread_t *thread, void *(*start_routine)(void *), 
                         void *arg, char *msg) {
  assert(thread != NULL);
  assert(arg != NULL);

  perror_fail(pthread_create(thread, NULL, start_routine, arg), msg);
}


/** Safe thread join */
void safe_pthread_join(pthread_t thread, char *msg) {
  void *value_ptr; /* Dummy pointer for use by pthread_join */

  perror_fail(pthread_join(thread, &value_ptr), msg);
}


/** Descriptor of a sequential for loop */
typedef struct {
  unsigned int start; /**< for loop initial value */
  unsigned int bound; /**< for loop (strict) upper bound */
  void (*body_routine)(unsigned int, void*); 
  void *arg;
} for_sequential_t;


/** Call start_routine in a loop with a function signature as requested by 
    pthread_create(P) 
    @pre arg is not NULL */
static void *for_sequential(void *arg) {
  unsigned int i;
  for_sequential_t *descriptor;

  assert(arg != NULL);

  descriptor = (for_sequential_t *) arg;
  for (i = descriptor->start; i < descriptor->bound; i++) 
    descriptor->body_routine(i, descriptor->arg);

  return NULL;
}


/** Parallelize for loops with independent iterations */
void for_parallel(unsigned int start,      /**< for loop initial value */
                  unsigned int bound,      /**< for loop (strict) upper bound */
                  unsigned int thread_cnt, /**< Number of threads */
                  void (*body_routine)(unsigned int, void*), 
                  void *arg) {
  unsigned int i, step;
  pthread_t *thread;
  for_sequential_t *descriptor;

  /* Iterations per thread */
  step = (bound - start) / thread_cnt;

  /* Maximum possible number of threads, since step might have been rounded
     down */
  thread_cnt++;

  /* Allocate per-thread vectors */
  thread     = (pthread_t *) 
               scalloc(error_msg(__FILE__, __LINE__, __FUNCTION__,
                                 "cannot allocate thread vector"),
                       thread_cnt, sizeof(pthread_t));
  descriptor = (for_sequential_t *) 
               scalloc(error_msg(__FILE__, __LINE__, __FUNCTION__,
                                 "cannot allocate for-loop descriptors"),
                       thread_cnt, sizeof(for_sequential_t));

  for (i = 0; start < bound; i++) {
    /* Write arguments for the ith thread */
    assert(i < thread_cnt);
    descriptor[i].start = start;
    start += step;
    descriptor[i].bound = start;
    descriptor[i].body_routine = body_routine;
    descriptor[i].arg = arg;

    safe_pthread_create(thread+i, for_sequential, (void *) (descriptor+i));
  }

  /* Record the actual number of threads */
  thread_cnt = i;

  /* Join all threads */
  for (i = 0; i < thread_cnt; i++)
    safe_pthread_join(thread[i]);

  free(thread);
  free(arg);
}


#ifdef NODEF
  void (*start_routine)(unsigned int, void*); 
void jac_lppl_marshal(unsigned int i, void *arg) {
  jac_lppl_t *jac_arg;

  assert(arg != NULL);

  
}
/** Jacobian in vector format
    @pre p, hx, adata are not NULL; T > n */
void jac_lppl_vector(double *p, double *j, int m, int n, void *adata) {
  unsigned int i, start;
  lppl_env_t *env;
  jac_lppl_t jac_arg;

  assert(p   != NULL);
  assert(j   != NULL);
  assert(adata != NULL);
  assert(p[2] > (double) n);

  env = (lppl_env_t *) adata;
  start = env->start;

  for_parallel(0, n, THREAD_CNT, jac_lppl_marshal, jac_arg);
}
#endif
