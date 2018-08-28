/*
 * File: numerical_integral_emxutil.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

#ifndef NUMERICAL_INTEGRAL_EMXUTIL_H
#define NUMERICAL_INTEGRAL_EMXUTIL_H

/* Include Files */
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "numerical_integral_types.h"

/* Function Declarations */
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);

#endif

/*
 * File trailer for numerical_integral_emxutil.h
 *
 * [EOF]
 */
