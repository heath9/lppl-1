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

/** \brief Error handling */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <stdarg.h>

/** Length of the error string */
#define ERROR_LEN 256
/** Error string */
char error_str[ERROR_LEN];

/** Name with which the program was invoked */
const char *argv0 = "Error handling uninitialized";

/** Error indicator */
unsigned int error_state = 0;

/** Initialize the program name */
void argv0_init(const char *init_argv0) {
  argv0 = init_argv0;
}


/** Format an error message of format: 'argv0: msg'
    @pre str, msg are not NULL
    @return The error message */
char *error_msg(const char *file,  /**< File where error occurred */
                int line,          /**< Line where error occurred */
                const char *func,  /**< Name of function where error occurred */
                const char *msg = "" /**< Additional error information */) {
  assert(argv0 != NULL);
  assert(msg != NULL);

  return 
  #ifdef NDEBUG
    snprintf(error_str, ERROR_LEN, "%s: %s", argv0, msg) 
  #else
    snprintf(error_str, ERROR_LEN, "%s (%s:%i, %s): %s", 
             argv0, file, line, func, msg) 
  #endif
  < 0 ? NULL : error_str;
}

/** Same as error_msg but with empty error information.
    It should only be used if there is another way of generating the additional
    error message. */
char *error_msg(const char *file, int line, const char *func) {
  return error_msg(file, line, func, "");
}


/** Conditional (state!=0) fail after calling perror.
    No need to assert s not NULL, as per perror(3C) */
void perror_fail(int state, const char *s) {
  if (state) {
    perror(s);
    exit(EXIT_FAILURE);
  }
}


/** Conditional (state != 0) fail after printing an error message on 
    standard error. The arguments after state follow the same format
    as in printf(3) */
void die(int state, const char *fmt, ...) {
  va_list args;

  if (state) {
    va_start(args, fmt);
    (void) vfprintf(stderr, fmt, args);
    va_end(args);
    exit(EXIT_FAILURE);
  }
}


/** Conditional (state != 0) croak an error message on standard error, 
    and set the error_state variable. 
    The arguments after state follow the same format as in printf(3) */
void croak(int state, const char *fmt, ...) {
  va_list args;

  if (state) {
    va_start(args, fmt);
    (void) vfprintf(stderr, fmt, args);
    va_end(args);
    error_state = 1;
  }
}


/** Terminates the program if the error_state variable is set */
void croak_fail() {
  if (error_state)
    exit(EXIT_FAILURE);
}


/** Clear the error_state variable */
void croak_clear() {
  error_state = 0;
}
