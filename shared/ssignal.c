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

/** \brief Safe signal(2) and kill(2) */

#include <signal.h>
#include "shared.h"

/** Safe emulation of signal(2) with sigaction */
sighandler_t safe_signal(const char *s, int signum, sighandler_t handler) {
  struct sigaction act, oldact;

  /* Prepare signal handler */
  act.sa_handler = handler;
  perror_fail(sigemptyset(&act.sa_mask) == -1, s);
  act.sa_flags   = 0;

  /* Install signal handler */
  perror_fail(sigaction(signum, &act, &oldact) == -1, s);

  return oldact.sa_handler;
}


/** Safe kill */
void safe_kill(const char *s, pid_t pid, int sig) {
  perror_fail(kill(pid, sig) == -1, s);
}


int sigrecv = 0; /**< Whether a signal has been received */

int get_sigrecv() {
  return sigrecv;
}

void set_sigrecv(int init_sigrecv) {
  sigrecv = init_sigrecv;
}

/** Default signal handler: increment an internal counter */
void handler_sigrecv(int signum) {
  sigrecv++;
}

