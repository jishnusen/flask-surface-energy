/*
 * File: numerical_integral.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

#ifndef NUMERICAL_INTEGRAL_H
#define NUMERICAL_INTEGRAL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "numerical_integral_types.h"

/* Function Declarations */
extern void numerical_integral(double theta_chemi, double theta_physi, double c,
  double D, double d, double E, double e, double A, double BET, double entropy,
  double heat, double temperature, double b, double temp_gamma, double
  temp_diff_gamma, double f, double res_p_po[10000], double res_theta[10000],
  double res_gamma[10000], double res_uads[10000]);

#endif

/*
 * File trailer for numerical_integral.h
 *
 * [EOF]
 */
