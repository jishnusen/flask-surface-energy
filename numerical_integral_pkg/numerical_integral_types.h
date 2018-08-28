/*
 * File: numerical_integral_types.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

#ifndef NUMERICAL_INTEGRAL_TYPES_H
#define NUMERICAL_INTEGRAL_TYPES_H

/* Include Files */
#include "rtwtypes.h"

/* Type Definitions */
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/
#endif

/*
 * File trailer for numerical_integral_types.h
 *
 * [EOF]
 */
