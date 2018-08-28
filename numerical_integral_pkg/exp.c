/*
 * File: exp.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "numerical_integral.h"
#include "exp.h"

/* Function Definitions */

/*
 * Arguments    : emxArray_real_T *x
 * Return Type  : void
 */
void b_exp(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0];
  for (k = 0; k < nx; k++) {
    x->data[k] = exp(x->data[k]);
  }
}

/*
 * File trailer for exp.c
 *
 * [EOF]
 */
