/*
 * File: numerical_integral.c
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 28-Jul-2018 19:18:05
 */

/* Include Files */
#include <math.h>
#include "rt_nonfinite.h"
#include "numerical_integral.h"
#include "numerical_integral_emxutil.h"
#include "exp.h"
#include "power.h"
#include "rdivide.h"
#include "linspace.h"
#include "numerical_integral_rtwutil.h"

/* Function Definitions */

/*
 * Arguments    : double theta_chemi
 *                double theta_physi
 *                double c
 *                double D
 *                double d
 *                double E
 *                double e
 *                double A
 *                double BET
 *                double entropy
 *                double heat
 *                double temperature
 *                double b
 *                double temp_gamma
 *                double temp_diff_gamma
 *                double f
 *                double res_p_po[10000]
 *                double res_theta[10000]
 *                double res_gamma[10000]
 *                double res_uads[10000]
 * Return Type  : void
 */
void numerical_integral(double theta_chemi, double theta_physi, double c,
  double D, double d, double E, double e, double A, double BET, double entropy,
  double heat, double temperature, double b, double temp_gamma, double
  temp_diff_gamma, double f, double res_p_po[10000], double res_theta[10000],
  double res_gamma[10000], double res_uads[10000])
{
  /* setvbuf (stdout, NULL, _IONBF, 0); */
  double x[10001];
  double temp_diff_diff_gamma;
  double temp_uads;
  double temp_diff_uads;
  emxArray_real_T *p_po;
  emxArray_real_T *diff_ugas;
  emxArray_real_T *theta;
  emxArray_real_T *diff_theta;
  emxArray_real_T *diff_diff_theta;
  emxArray_real_T *diff_diff_H;
  emxArray_real_T *SA;
  emxArray_real_T *diff_diff_gamma;
  emxArray_real_T *diff_gamma;
  emxArray_real_T *b_gamma;
  emxArray_real_T *diff_uads;
  emxArray_real_T *uads;
  emxArray_real_T *y;
  emxArray_real_T *b_diff_gamma;
  emxArray_real_T *b_x;
  emxArray_real_T *a;
  double b_a;
  double c_a;
  double c_x;
  double b_theta_physi;
  double b_c;
  double d_a;
  double e_a;
  double f_a;
  double b_y;
  double c_y;
  int n;
  double ndbl;
  int k;
  double apnd;
  double cdiff;
  double absa;
  double absb;
  int nm1d2;
  int b_n;
  unsigned int p_po_idx_0;
  unsigned int b_p_po_idx_0;

  /*  step size */
  /*  K is the number of times we will loop through the numerical integration */
  /*  and will also be the number of data points to save */
  /*  number of points to store */
  /*  linearly space out the p_po data from 1-1E-12 to 0. The choice */
  /*  of 1-1E-12 occurred because we were able to use the data from the first  */
  /*  integration with dx = 1E- 19. This will allow us to use the boundary  */
  /*  conditions without worrying about NAN from inf*0. The number of values in */
  /*  x is k+1 */
  linspace(x);
  temp_diff_diff_gamma = 0.0;
  temp_uads = 0.0;
  temp_diff_uads = 0.0;

  /*  run the integration over k loops */
  emxInit_real_T(&p_po, 1);
  emxInit_real_T(&diff_ugas, 1);
  emxInit_real_T(&theta, 1);
  emxInit_real_T(&diff_theta, 1);
  emxInit_real_T(&diff_diff_theta, 1);
  emxInit_real_T(&diff_diff_H, 1);
  emxInit_real_T(&SA, 1);
  emxInit_real_T(&diff_diff_gamma, 1);
  emxInit_real_T(&diff_gamma, 1);
  emxInit_real_T(&b_gamma, 1);
  emxInit_real_T(&diff_uads, 1);
  emxInit_real_T(&uads, 1);
  emxInit_real_T1(&y, 2);
  emxInit_real_T(&b_diff_gamma, 1);
  emxInit_real_T(&b_x, 1);
  emxInit_real_T(&a, 1);
  b_a = theta_chemi * b;
  c_a = theta_physi * c;
  c_x = theta_chemi * b;
  b_theta_physi = theta_physi * c;
  b_c = (c - 1.0) * (c - 1.0);
  d_a = 3.0 * b;
  e_a = -theta_chemi * b;
  f_a = 2.0 * theta_physi * c;
  b_y = -D / d;
  c_y = E / e;
  for (n = 0; n < 10000; n++) {
    /*  set up the data needed to run the integration for the given p_po */
    /*  range */
    if (x[n] < x[n + 1]) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity_real_T(y, k);
    } else {
      ndbl = floor((x[n + 1] - x[n]) / -1.0E-9 + 0.5);
      apnd = x[n] + ndbl * -1.0E-9;
      cdiff = x[n + 1] - apnd;
      absa = fabs(x[n]);
      absb = fabs(x[n + 1]);
      if ((absa > absb) || rtIsNaN(absb)) {
        absb = absa;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = x[n + 1];
      } else if (cdiff > 0.0) {
        apnd = x[n] + (ndbl - 1.0) * -1.0E-9;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        b_n = (int)ndbl;
      } else {
        b_n = 0;
      }

      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = b_n;
      emxEnsureCapacity_real_T(y, k);
      if (b_n > 0) {
        y->data[0] = x[n];
        if (b_n > 1) {
          y->data[b_n - 1] = apnd;
          nm1d2 = (b_n - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            ndbl = (double)k * -1.0E-9;
            y->data[k] = x[n] + ndbl;
            y->data[(b_n - k) - 1] = apnd - ndbl;
          }

          if (nm1d2 << 1 == b_n - 1) {
            y->data[nm1d2] = (x[n] + apnd) / 2.0;
          } else {
            ndbl = (double)nm1d2 * -1.0E-9;
            y->data[nm1d2] = x[n] + ndbl;
            y->data[nm1d2 + 1] = apnd - ndbl;
          }
        }
      }
    }

    k = p_po->size[0];
    p_po->size[0] = y->size[1];
    emxEnsureCapacity_real_T1(p_po, k);
    nm1d2 = y->size[1];
    for (k = 0; k < nm1d2; k++) {
      p_po->data[k] = y->data[y->size[0] * k];
    }

    /*  dugas/dp_po */
    k = diff_ugas->size[0];
    diff_ugas->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(diff_ugas, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_ugas->data[k] = 2477.572 / p_po->data[k];
    }

    /*  theta */
    power(p_po, diff_diff_gamma);
    power(p_po, diff_gamma);
    k = a->size[0];
    a->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(a, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      a->data[k] = b_a * diff_diff_gamma->data[k];
    }

    k = b_x->size[0];
    b_x->size[0] = diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 1.0 + b * diff_gamma->data[k];
    }

    rdivide(a, b_x, theta);
    k = a->size[0];
    a->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(a, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      a->data[k] = c_a * p_po->data[k];
    }

    k = b_x->size[0];
    b_x->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = (1.0 - p_po->data[k]) * (p_po->data[k] * (c - 1.0) + 1.0);
    }

    rdivide(a, b_x, diff_diff_gamma);
    k = theta->size[0];
    emxEnsureCapacity_real_T1(theta, k);
    nm1d2 = theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      theta->data[k] += diff_diff_gamma->data[k];
    }

    /*  dtheta/dp_po */
    power(p_po, diff_diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 1.0 + b * diff_diff_gamma->data[k];
    }

    b_power(b_x, diff_diff_gamma);
    power(p_po, diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = c_x / diff_diff_gamma->data[k] / 2.0;
    }

    rdivide(b_x, diff_gamma, diff_theta);
    b_power(p_po, diff_diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 1.0 - p_po->data[k];
    }

    b_power(b_x, diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = p_po->data[k] * (c - 1.0) + 1.0;
    }

    b_power(b_x, SA);
    k = b_x->size[0];
    b_x->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = b_theta_physi * (diff_diff_gamma->data[k] * (c - 1.0) + 1.0);
    }

    k = b_diff_gamma->size[0];
    b_diff_gamma->size[0] = diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_diff_gamma, k);
    nm1d2 = diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_diff_gamma->data[k] = diff_gamma->data[k] * SA->data[k];
    }

    rdivide(b_x, b_diff_gamma, diff_diff_gamma);
    k = diff_theta->size[0];
    emxEnsureCapacity_real_T1(diff_theta, k);
    nm1d2 = diff_theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_theta->data[k] += diff_diff_gamma->data[k];
    }

    /*  d^2theta/d(p_po)^2 */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    b_p_po_idx_0 = (unsigned int)p_po->size[0];
    k = SA->size[0];
    SA->size[0] = (int)b_p_po_idx_0;
    emxEnsureCapacity_real_T1(SA, k);
    for (k = 0; k < (int)p_po_idx_0; k++) {
      SA->data[k] = rt_powd_snf(p_po->data[k], 1.5);
    }

    power(p_po, diff_diff_gamma);
    power(p_po, diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 1.0 + b * diff_gamma->data[k];
    }

    c_power(b_x, diff_gamma);
    k = a->size[0];
    a->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(a, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      a->data[k] = e_a * (d_a * diff_diff_gamma->data[k] + 1.0);
    }

    k = b_x->size[0];
    b_x->size[0] = SA->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = SA->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 4.0 * SA->data[k] * diff_gamma->data[k];
    }

    rdivide(a, b_x, diff_diff_theta);
    c_power(p_po, diff_diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = 1.0 - p_po->data[k];
    }

    c_power(b_x, diff_gamma);
    k = b_x->size[0];
    b_x->size[0] = p_po->size[0];
    emxEnsureCapacity_real_T1(b_x, k);
    nm1d2 = p_po->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_x->data[k] = p_po->data[k] * (c - 1.0) + 1.0;
    }

    c_power(b_x, SA);
    k = a->size[0];
    a->size[0] = diff_diff_gamma->size[0];
    emxEnsureCapacity_real_T1(a, k);
    nm1d2 = diff_diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      a->data[k] = f_a * (((diff_diff_gamma->data[k] * b_c + 3.0 * p_po->data[k]
                            * (c - 1.0)) + 2.0) - c);
    }

    k = b_diff_gamma->size[0];
    b_diff_gamma->size[0] = diff_gamma->size[0];
    emxEnsureCapacity_real_T1(b_diff_gamma, k);
    nm1d2 = diff_gamma->size[0];
    for (k = 0; k < nm1d2; k++) {
      b_diff_gamma->data[k] = diff_gamma->data[k] * SA->data[k];
    }

    rdivide(a, b_diff_gamma, diff_diff_gamma);
    k = diff_diff_theta->size[0];
    emxEnsureCapacity_real_T1(diff_diff_theta, k);
    nm1d2 = diff_diff_theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_diff_theta->data[k] += diff_diff_gamma->data[k];
    }

    /*  d^2(H)/d(theta)^2 */
    k = diff_diff_H->size[0];
    diff_diff_H->size[0] = theta->size[0];
    emxEnsureCapacity_real_T1(diff_diff_H, k);
    nm1d2 = theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_diff_H->data[k] = -theta->data[k] / d;
    }

    b_exp(diff_diff_H);
    k = diff_diff_gamma->size[0];
    diff_diff_gamma->size[0] = theta->size[0];
    emxEnsureCapacity_real_T1(diff_diff_gamma, k);
    nm1d2 = theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_diff_gamma->data[k] = -theta->data[k] / e;
    }

    b_exp(diff_diff_gamma);
    b_power(theta, diff_gamma);
    k = SA->size[0];
    SA->size[0] = theta->size[0];
    emxEnsureCapacity_real_T1(SA, k);
    nm1d2 = theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      SA->data[k] = -theta->data[k] / e;
    }

    b_exp(SA);
    k = diff_diff_H->size[0];
    emxEnsureCapacity_real_T1(diff_diff_H, k);
    nm1d2 = diff_diff_H->size[0];
    for (k = 0; k < nm1d2; k++) {
      diff_diff_H->data[k] = (b_y * diff_diff_H->data[k] + E * (f - 2.0 *
        theta->data[k]) * diff_diff_gamma->data[k]) - c_y * (f * theta->data[k]
        - diff_gamma->data[k]) * SA->data[k];
    }

    /*  surface AREA */
    k = SA->size[0];
    SA->size[0] = theta->size[0];
    emxEnsureCapacity_real_T1(SA, k);
    nm1d2 = theta->size[0];
    for (k = 0; k < nm1d2; k++) {
      SA->data[k] = BET + A * theta->data[k];
    }

    /*  dSA/dtheta */
    /*  d^2(gamma)/d(p_po)^2 */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    k = diff_diff_gamma->size[0];
    diff_diff_gamma->size[0] = (int)p_po_idx_0;
    emxEnsureCapacity_real_T1(diff_diff_gamma, k);
    nm1d2 = (int)p_po_idx_0;
    for (k = 0; k < nm1d2; k++) {
      diff_diff_gamma->data[k] = 0.0;
    }

    /*  dgamma/dp_po */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    k = diff_gamma->size[0];
    diff_gamma->size[0] = (int)p_po_idx_0;
    emxEnsureCapacity_real_T1(diff_gamma, k);
    nm1d2 = (int)p_po_idx_0;
    for (k = 0; k < nm1d2; k++) {
      diff_gamma->data[k] = 0.0;
    }

    /*  gamma */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    k = b_gamma->size[0];
    b_gamma->size[0] = (int)p_po_idx_0;
    emxEnsureCapacity_real_T1(b_gamma, k);
    nm1d2 = (int)p_po_idx_0;
    for (k = 0; k < nm1d2; k++) {
      b_gamma->data[k] = 0.0;
    }

    /*  d(uads)/dp_po */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    k = diff_uads->size[0];
    diff_uads->size[0] = (int)p_po_idx_0;
    emxEnsureCapacity_real_T1(diff_uads, k);
    nm1d2 = (int)p_po_idx_0;
    for (k = 0; k < nm1d2; k++) {
      diff_uads->data[k] = 0.0;
    }

    /*  uads */
    p_po_idx_0 = (unsigned int)p_po->size[0];
    k = uads->size[0];
    uads->size[0] = (int)p_po_idx_0;
    emxEnsureCapacity_real_T1(uads, k);
    nm1d2 = (int)p_po_idx_0;
    for (k = 0; k < nm1d2; k++) {
      uads->data[k] = 0.0;
    }

    /*  check if it is the first loop */
    if (1 + n == 1) {
      /*  integrate the data for all values in p_po array */
      for (nm1d2 = 0; nm1d2 < p_po->size[0]; nm1d2++) {
        if (1 + nm1d2 == 1) {
          /*  if n== 1 and i==1 we are at the boundary and thus must */
          /*  apply the boundary conditions from previous integration. */
          b_gamma->data[0] = temp_gamma;
          uads->data[0] = -heat + temperature * entropy;
          diff_gamma->data[0] = temp_diff_gamma;
        } else {
          /*  Else numerically integrate the function. Note the */
          /*  negative sign is due to dx being negative */
          diff_gamma->data[nm1d2] = diff_gamma->data[nm1d2 - 1] - 1.0E-9 *
            diff_diff_gamma->data[nm1d2 - 1];
          b_gamma->data[nm1d2] = b_gamma->data[nm1d2 - 1] - 1.0E-9 *
            diff_gamma->data[nm1d2 - 1];
          diff_uads->data[nm1d2] = -SA->data[nm1d2] / theta->data[nm1d2] *
            diff_gamma->data[nm1d2];
          uads->data[nm1d2] = uads->data[nm1d2 - 1] - 1.0E-9 * diff_uads->
            data[nm1d2 - 1];
        }

        /*  Apply the boundary equation for gamma to use in the next */
        /*  iteration for the data. */
        diff_diff_gamma->data[nm1d2] = diff_theta->data[nm1d2] / SA->data[nm1d2]
          * ((diff_diff_H->data[nm1d2] * diff_theta->data[nm1d2] +
              diff_ugas->data[nm1d2]) - diff_gamma->data[nm1d2] * ((2.0 * A -
               SA->data[nm1d2] / theta->data[nm1d2]) - SA->data[nm1d2] *
              diff_diff_theta->data[nm1d2] / (diff_theta->data[nm1d2] *
               diff_theta->data[nm1d2])));
      }
    } else {
      /*  integrate the data for all values in p_po array expect the first */
      /*  which has already been counted */
      for (nm1d2 = 1; nm1d2 - 1 <= p_po->size[0] - 2; nm1d2++) {
        if (1 + nm1d2 == 2) {
          /*  Else numerically integrate the function. Note the */
          /*  negative sign is due to dx being negative */
          diff_gamma->data[1] = temp_diff_gamma - 1.0E-9 * temp_diff_diff_gamma;
          b_gamma->data[1] = temp_gamma - 1.0E-9 * temp_diff_gamma;
          diff_uads->data[1] = -SA->data[1] / theta->data[1] * diff_gamma->data
            [1];
          uads->data[1] = temp_uads - 1.0E-9 * temp_diff_uads;
        } else {
          /*  Else numerically integrate the function. Note the */
          /*  negative sign is due to dx being negative */
          diff_gamma->data[nm1d2] = diff_gamma->data[nm1d2 - 1] - 1.0E-9 *
            diff_diff_gamma->data[nm1d2 - 1];
          b_gamma->data[nm1d2] = b_gamma->data[nm1d2 - 1] - 1.0E-9 *
            diff_gamma->data[nm1d2 - 1];
          diff_uads->data[nm1d2] = -SA->data[nm1d2] / theta->data[nm1d2] *
            diff_gamma->data[nm1d2];
          uads->data[nm1d2] = uads->data[nm1d2 - 1] - 1.0E-9 * diff_uads->
            data[nm1d2 - 1];
        }

        /*  Apply the boundary equation for gamma to use in the next */
        /*  iteration for the data. */
        diff_diff_gamma->data[nm1d2] = diff_theta->data[nm1d2] / SA->data[nm1d2]
          * ((diff_diff_H->data[nm1d2] * diff_theta->data[nm1d2] +
              diff_ugas->data[nm1d2]) - diff_gamma->data[nm1d2] * ((2.0 * A -
               SA->data[nm1d2] / theta->data[nm1d2]) - SA->data[nm1d2] *
              diff_diff_theta->data[nm1d2] / (diff_theta->data[nm1d2] *
               diff_theta->data[nm1d2])));
      }
    }

    /*  save some temporary variables to use in the next loop so that we can */
    /*  overwrite the old data */
    temp_gamma = b_gamma->data[b_gamma->size[0] - 1];
    temp_diff_gamma = diff_gamma->data[diff_gamma->size[0] - 1];
    temp_diff_diff_gamma = diff_diff_gamma->data[diff_diff_gamma->size[0] - 1];
    temp_diff_uads = diff_uads->data[diff_uads->size[0] - 1];
    temp_uads = uads->data[uads->size[0] - 1];
    res_p_po[n] = p_po->data[p_po->size[0] - 1];
    res_theta[n] = theta->data[theta->size[0] - 1] * 602214.13 / BET;
    res_gamma[n] = b_gamma->data[b_gamma->size[0] - 1];
    res_uads[n] = uads->data[uads->size[0] - 1] / 1000.0;
  }

  emxFree_real_T(&a);
  emxFree_real_T(&b_x);
  emxFree_real_T(&b_diff_gamma);
  emxFree_real_T(&y);
  emxFree_real_T(&uads);
  emxFree_real_T(&diff_uads);
  emxFree_real_T(&b_gamma);
  emxFree_real_T(&diff_gamma);
  emxFree_real_T(&diff_diff_gamma);
  emxFree_real_T(&SA);
  emxFree_real_T(&diff_diff_H);
  emxFree_real_T(&diff_diff_theta);
  emxFree_real_T(&diff_theta);
  emxFree_real_T(&theta);
  emxFree_real_T(&diff_ugas);
  emxFree_real_T(&p_po);
}

/*
 * File trailer for numerical_integral.c
 *
 * [EOF]
 */
