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

/** 

\file smalloc.c
\brief Safe malloc

smalloc and scalloc are like malloc and calloc, but exit with failure
if they cannot allocate memory.
Take an extra const char * for use by perror.
The term "safe" has nothing to do here with thread-safety.

A general convention is that all memory allocation must use safe malloc or 
safe calloc.
*/

#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "shared.h"

/** Global mutex around smalloc functions */
pthread_mutex_t smalloc_mutex = PTHREAD_MUTEX_INITIALIZER;

/** Safe malloc.
    No need to assert(s != NULL) as per perror(3C) */
void *smalloc(size_t size, const char *s) {
  void *ptr;

  slock(&smalloc_mutex, s);
  perror_fail((ptr = malloc(size)) == NULL, s);
  sunlock(&smalloc_mutex, s);

  debug_print("smalloc: %u\n", size);
  return ptr;
}


/** Safe calloc
    No need to assert(s != NULL) as per perror(3C) */
void *scalloc(size_t nelem, size_t elsize, const char *s) { 
  void *ptr;

  slock(&smalloc_mutex, s);
  perror_fail((ptr = calloc(nelem,elsize)) == NULL, s);
  sunlock(&smalloc_mutex, s);

  debug_print("scalloc: %u x %u\n", nelem, elsize);
  return ptr;
}


/** Safe free. If ptr is null, nothing is freed (but the alloc mutex is quickly
    locked and unlocked 
    No need to assert(s != NULL) as per perror(3C),
    no need to assert(ptr != NULL) as per malloc(3) */
void sfree(void *ptr, const char *s) {

  slock(&smalloc_mutex, s);
  free(ptr);
  sunlock(&smalloc_mutex, s);
}
