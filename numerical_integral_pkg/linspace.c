/*
 * File: linspace.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "numerical_integral.h"
#include "linspace.h"

/* Function Definitions */

/*
 * Arguments    : double y[10001]
 * Return Type  : void
 */
void linspace(double y[10001])
{
  int k;
  y[10000] = 0.0;
  y[0] = 0.999999999999;
  for (k = 0; k < 9999; k++) {
    y[k + 1] = 0.999999999999 + (1.0 + (double)k) * -9.99999999999E-5;
  }
}

/*
 * File trailer for linspace.c
 *
 * [EOF]
 */
