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

/** \brief Shared error handling and system routines */

#include <stdlib.h>
#include <pthread.h>
#include <sys/types.h>

/* Error handling */
#define HERE __FILE__, __LINE__, __FUNCTION__
void argv0_init(const char *);
char *error_msg(const char *file, int line, const char *func, const char *msg);
char *error_msg(const char *file, int line, const char *func);
void perror_fail(int, const char *);
void die(int state, const char *fmt, ...);
void croak(int state, const char *fmt, ...);
void croak_fail();
void croak_clear();

/* Debug print-outs */
void debug_print_internal(const char *fmt, ...);
#ifdef DEBUG
  #define debug_print(...) debug_print_internal(__VA_ARGS__);
#else
  #define debug_print(...) (void)(0)
#endif

/* Locks */
void slock(pthread_mutex_t *mutex, const char *s);
void sunlock(pthread_mutex_t *mutex, const char *s);

/* Memory allocation */
void *smalloc(size_t size, const char *s);
void *scalloc(size_t nelem, size_t elsize, const char *s);
void sfree(void *ptr, const char *s);

/* Interprocess communication */
typedef void (*sighandler_t)(int);
sighandler_t safe_signal(const char *, int, sighandler_t);
void handler_sigrecv(int);
int get_sigrecv();
void set_sigrecv(int);
void safe_kill(const char *, pid_t, int);
void *smmap(const char *, size_t);
void *scmmap(const char *, size_t nelem, size_t elsize);
void smunmap(const char *, void *, size_t);
void scmunmap(const char *, void *, size_t, size_t);
pid_t sfork(const char *);
pid_t swait(const char *, int *);

#define LG "%.16g" /**< Format for printing %g */

/* Time */
double time_gap(const char *, struct timeval *);

/* Vector operations */
void vector_copy(unsigned int, double *dst, const double *src);
double dot_product(unsigned int, double *, double *);
void down_shift(unsigned int, double *);
void vector_f(unsigned int, double x[], double (*f)(double));
double stable_sum(unsigned int, const double y[]);
double vector_error(unsigned int, double x[], double y[]);
void vector_add(unsigned int, double c[], double a[], double b[]);

/* Minimum and maximum */
unsigned int umin(unsigned int, unsigned int);
unsigned int umax(unsigned int, unsigned int);
double fmin(double, double);
double fmax(double, double);

/* Random numbers */
void ssrand48(pthread_mutex_t *, unsigned long seed);
double sdrand48(pthread_mutex_t *);
void generate_normal(double g[], pthread_mutex_t *);
