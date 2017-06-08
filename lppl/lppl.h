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

/** \brief Header file for the LPPL interpolator */

/** Indexing (price vector): 
    Time	Price
    1		x[0]
    2		x[1]
    ...
    n		x[n-1]
    Note that some functions may normalize it so that 
    x = x + env->start - 1
    n = n - env->start + 1
*/

#include <stdio.h>
#include <sys/types.h>
#include "lm.h"

#define N 35000 /**< Maximum price file length */

typedef struct {
  unsigned int verbose;  /**< Whether to use verbose output */
  char *pricefile;       /**< Price file name */
  char *parmsfile;       /**< Parameter file name */
  char *rangefile;       /**< Range file name */
  int maxiter;           /**< Maximum number of iterations */
  int lmiter;            /**< Maximum number of iterations per LM computation */
  int printiter;         /**< Number of LM iteration blocks between 
                              print-outs */
  double epsilon;        /**< Termination error (epsilon) */
  unsigned int start;    /**< Fit starting point */
  unsigned int threads;  /**< OpenMP threads */
  double (*f)(double *, double); /**< Interpolation function */
  unsigned int expfit;   /**< Whether to start from an exponential fit */
  unsigned int linprice; /**< Whether to use price (no log) */
  double window;         /**< Window parameter for weight calculation */
  double noise;          /**< Noise power */
  double diffuse;        /**< Diffusion constant in Ornstein-Uhlenbeck */
  double relax;          /**< Relaxation constant in Ornstein-Uhlenbeck */
  unsigned int seed;     /**< Random generator seed */
  unsigned int nopenny;  /**< Do not round prices to the closest penny */
} lppl_env_t;


/** Output parameters of levmar fit */
typedef struct {
  int iterations;          /**< Fit iterations */
  double p[7];             /**< LPPL parameters */
  double covar[7][7];      /**< Covariance matrix */
  double info[LM_INFO_SZ]; /**< Fit information */
} levmar_t;


/** Levmar iterations data structure */
typedef struct {
  unsigned int adaptive;  /**< Whether the iteration number should dynamically
                               adjust to running time and error reduction */
  unsigned int start_up;  /**< Whether the number of iterations is being 
                               ramped up exponentially in the inital part
                               of the algorithm */
  int iterations;         /**< Current number of iterations */
} levmar_iter_t;


double lppl(double *p, double x);
double lppl_B(double *p, double x);
void normalize_lppl(double *p);
double step_lppl(double *p, double x);
void lppl_vector(double *p, double *hx, int m, int n, void *adata);
void parse_lppl_args(int, const char **, lppl_env_t *);
void lppl_print(int iterations, unsigned int df, double p[], double info[], 
                double covar[][7]);
void print_lppl_parms(double p[]);
unsigned int read_lppl_parms(FILE *, double p[], unsigned int verbose);
void safe_read_parms(double p[], lppl_env_t *);
unsigned int read_all_parms(double p[], levmar_t lm[], unsigned int, 
                            lppl_env_t *);
void jac_lppl(double p[], double j[], double x);
void jac_lppl_vector(double *p, double *j, int m, int n, void *adata);
int levmar(double *p, double *x, unsigned int n, double info[], 
           double covar[][7], lppl_env_t *env);
void linear_lppl_guess(double p[], double x[], unsigned int, lppl_env_t *);
void lppl_guess(double p[], double x[], double rho, double xmin, double xmax, 
                unsigned int n, lppl_env_t *); 
void jac_linear_lppl(double p[], double j[], double x);
void jac_linear_lppl_vector(double p[], double j[][3], int n, lppl_env_t *);
int linear_lsq(double p[], double x[], int n, lppl_env_t *);
double levmar_weight(unsigned int i, unsigned int n, lppl_env_t *);
void weigh_price(unsigned int, double *, lppl_env_t *);
void iter_init(int init_iterations, levmar_iter_t *lmiter);
void iter_update(double levmar_time, double levmar_err, int levmar_iter,
                 double lin_time,    double lin_err,
                 double levmar_stat_time[], int levmar_stat_iter[],
                 levmar_iter_t *);
double solution_error(unsigned int n, double x[], double p[], lppl_env_t *);
int exp_lsq(double p[], double x[], int n, lppl_env_t *);
void safe_exp_lsq(double p[], double x[], int n, lppl_env_t *);
void brownian_motion(unsigned int, double b[], pthread_mutex_t *, lppl_env_t *);
void ou_noise(unsigned int, double ou[], pthread_mutex_t *, lppl_env_t *);

FILE *safe_fopen(const char *filename, char *error_msg);
FILE *safe_popen(char *command, char *type, char *error_msg);
void safe_pclose(FILE *fileptr, char *error_msg);

unsigned int read_price(unsigned int, double price[], FILE *);

double rsquared(unsigned int, const double x[], const double y[]);
double correlation(unsigned int n, const double x[], const double y[]);
double autocorrelation(unsigned int n, double x[], unsigned int step);
