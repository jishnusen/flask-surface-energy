/*
 * File: rdivide.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "numerical_integral.h"
#include "rdivide.h"
#include "numerical_integral_emxutil.h"

/* Function Definitions */

/*
 * Arguments    : const emxArray_real_T *x
 *                const emxArray_real_T *y
 *                emxArray_real_T *z
 * Return Type  : void
 */
void rdivide(const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T
             *z)
{
  int i0;
  int loop_ub;
  i0 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity_real_T1(z, i0);
  loop_ub = x->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    z->data[i0] = x->data[i0] / y->data[i0];
  }
}

/*
 * File trailer for rdivide.c
 *
 * [EOF]
 */
