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

\file slock.c
\brief Safe locks

slock_init, slock, and sunlock are like pthread_muitex_init,
pthread_mutex_lock and pthread_mutex_unlock
but exit with failure if they cannot lock and unlock the mutex.
Take an extra const char * for use by perror.

A general convention is that all locks and unlock must use slock and sunlock.
*/

#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "shared.h"

/** Safe lock initializaton.
    @pre s, mutex are not NULL */
void slock_init(pthread_mutex_t *mutex, const char *s) {
  assert(s     != NULL);
  assert(mutex != NULL);

  perror_fail(pthread_mutex_init(mutex, NULL) != 0, s);
}


/** Safe lock 
    @pre s, mutex are not NULL */
void slock(pthread_mutex_t *mutex, const char *s) {
  assert(s     != NULL);
  assert(mutex != NULL);
  debug_print("LOCK %p\n", (void *) mutex);

  perror_fail(pthread_mutex_lock(mutex) != 0, s);
  debug_print("LOCK %p done\n", (void *) mutex);
}


/** Safe unlock 
    @pre s, mutex are not NULL */
void sunlock(pthread_mutex_t *mutex, const char *s) {
  assert(s     != NULL);
  assert(mutex != NULL);
  debug_print("UNLOCK %p\n", (void *) mutex);

  perror_fail(pthread_mutex_unlock(mutex) != 0, s);
}
