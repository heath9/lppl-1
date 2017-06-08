/** \brief Computation of correlation coefficients */

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "shared.h"


/** Computation of R-squared (coefficient of determination) */
double correlation(unsigned int n, const double x[], const double y[]) {
  unsigned int i;
  double x_mean, y_mean, *x_gap, *y_gap, *xy_gap, x_term, y_term, rho, 
         covariance;

  /* Calculate the mean of the y's */
  x_mean = stable_sum(n, x) / n;
  y_mean = stable_sum(n, y) / n;

  /* Allocate sum-of-square arrays */
  x_gap  = (double *) scalloc(n, sizeof(double), 
                              error_msg(HERE, "cannot allocate x(gap)"));
  y_gap  = (double *) scalloc(n, sizeof(double), 
                              error_msg(HERE, "cannot allocate y(gap)"));
  xy_gap = (double *) scalloc(n, sizeof(double), 
                              error_msg(HERE, "cannot allocate xy(gap)"));

  /* Calculate all the squared errors */
  for (i = 0; i < n; i++) {
    x_term = x[i] - x_mean;
    y_term = y[i] - y_mean;
    x_gap[i]  = x_term * x_term;
    y_gap[i]  = y_term * y_term;
    xy_gap[i] = x_term * y_term;
  }

  /* Calculate correlation */
  covariance = stable_sum(n, xy_gap);
  debug_print("Covariance " LG "\n", covariance);
  rho = covariance / sqrt(stable_sum(n, x_gap) * stable_sum(n, y_gap));
  debug_print("correlation " LG "\n", rho);

  /* Free sum-of-squares arrays */
  free(x_gap);
  free(y_gap);
  free(xy_gap);

  return rho;
}


/** Correlation */
double rsquared(unsigned int n, const double x[], const double y[]) {
  double rho;

  rho = correlation(n, x, y);
  return rho * rho;
}

 
/** Auto-correlation at the given step */
double autocorrelation(unsigned int n, double x[], unsigned int step) {
  return correlation(n - step, x, x + step);
}
