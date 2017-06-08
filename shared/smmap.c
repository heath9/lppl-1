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

/** \brief smalloc emulation in shared memory */

#include <pthread.h>
#include <sys/mman.h>
#include "shared.h"

/** Global mutex around smmap functions */
pthread_mutex_t smmap_mutex = PTHREAD_MUTEX_INITIALIZER;

/** Memory mapping, safely, with an interface purposefully similar to that
    of smalloc */
void *smmap(const char *s, size_t size) {
  void *ptr;

  slock(&smmap_mutex, s);
  perror_fail(
    (ptr =  mmap(NULL, size, PROT_READ | PROT_WRITE, 
                 MAP_SHARED | MAP_ANONYMOUS, -1, 0)) == NULL, s);
  sunlock(&smmap_mutex, s);

  debug_print("smmap: %u\n", size);
  return ptr;
}


/** Memory mapping, safely, but with an interface purposefully similar to that
    of scalloc */
void *scmmap(const char *s, size_t nelem, size_t elsize) {
  return smmap(s, nelem * elsize);
}


/** Memory unmapping, safely */
void smunmap(const char *s, void *ptr, size_t length) {
  slock(&smmap_mutex, s);
  perror_fail(munmap(ptr, length), s);
  sunlock(&smmap_mutex, s);
}


/** Memory unmapping, safely, with a scalloc interface */
void scmunmap(const char *s, void *ptr, size_t nelem, size_t elsize) {
  smunmap(s, ptr, nelem * elsize);
}
