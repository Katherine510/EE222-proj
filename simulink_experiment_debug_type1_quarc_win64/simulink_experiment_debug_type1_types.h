/*
 * simulink_experiment_debug_type1_types.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.6
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr 30 14:09:34 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_simulink_experiment_debug_type1_types_h_
#define RTW_HEADER_simulink_experiment_debug_type1_types_h_
#include "rtwtypes.h"
#ifndef struct_tag_cRiPtJOtbnGUmsbmlLYZUC
#define struct_tag_cRiPtJOtbnGUmsbmlLYZUC

struct tag_cRiPtJOtbnGUmsbmlLYZUC
{
  real_T t_prev;
  real_T p_prev;
  real_T theta_prev;
  real_T V_servo;
  real_T B[4];
  real_T C[8];
  real_T x_hat[4];
  real_T P[16];
};

#endif                                 /* struct_tag_cRiPtJOtbnGUmsbmlLYZUC */

#ifndef typedef_studentControllerInterface_si_T
#define typedef_studentControllerInterface_si_T

typedef struct tag_cRiPtJOtbnGUmsbmlLYZUC studentControllerInterface_si_T;

#endif                             /* typedef_studentControllerInterface_si_T */

/* Parameters (default storage) */
typedef struct P_simulink_experiment_debug_t_T_ P_simulink_experiment_debug_t_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_simulink_experiment_d_T RT_MODEL_simulink_experiment__T;

#endif                 /* RTW_HEADER_simulink_experiment_debug_type1_types_h_ */
