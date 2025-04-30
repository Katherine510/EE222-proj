/*
 * simulink_experiment_debug_type1.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_debug_type1".
 *
 * Model version              : 13.6
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr 30 13:33:16 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "simulink_experiment_debug_type1.h"
#include "rtwtypes.h"
#include "simulink_experiment_debug_type1_types.h"
#include "simulink_experiment_debug_type1_private.h"
#include <math.h>
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <string.h>
#include "simulink_experiment_debug_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_debug_t_T simulink_experiment_debug_typ_B;

/* Block states (default storage) */
DW_simulink_experiment_debug__T simulink_experiment_debug_ty_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_debug_ty_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_debug_ty_M =
  &simulink_experiment_debug_ty_M_;
static void rate_monotonic_scheduler(void);
time_T rt_SimUpdateDiscreteEvents(
  int_T rtmNumSampTimes, void *rtmTimingData, int_T *rtmSampleHitPtr, int_T
  *rtmPerTaskSampleHits )
{
  rtmSampleHitPtr[1] = rtmStepTask(simulink_experiment_debug_ty_M, 1);
  rtmSampleHitPtr[2] = rtmStepTask(simulink_experiment_debug_ty_M, 2);
  UNUSED_PARAMETER(rtmNumSampTimes);
  UNUSED_PARAMETER(rtmTimingData);
  UNUSED_PARAMETER(rtmPerTaskSampleHits);
  return(-1);
}

/*
 *         This function updates active task flag for each subrate
 *         and rate transition flags for tasks that exchange data.
 *         The function assumes rate-monotonic multitasking scheduler.
 *         The function must be called at model base rate so that
 *         the generated code self-manages all its subrates and rate
 *         transition flags.
 */
static void rate_monotonic_scheduler(void)
{
  /* To ensure a deterministic data transfer between two rates,
   * data is transferred at the priority of a fast task and the frequency
   * of the slow task.  The following flags indicate when the data transfer
   * happens.  That is, a rate interaction flag is set true when both rates
   * will run, and false otherwise.
   */

  /* tid 1 shares data with slower tid rate: 2 */
  if (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[1] == 0) {
    simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2 =
      (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] == 0);

    /* update PerTaskSampleHits matrix for non-inline sfcn */
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits[5] =
      simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2;
  }

  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2])++;
  if ((simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2]) > 4) {/* Sample time: [0.01s, 0.0s] */
    simulink_experiment_debug_ty_M->Timing.TaskCounters.TID[2] = 0;
  }
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T tmp;
  real_T tmp_0;
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  __m128d tmp;
  __m128d tmp_0;
  studentControllerInterface_si_T *obj;
  studentControllerInterface_si_T *obj_0;
  real_T A_lin[16];
  real_T P_p[16];
  real_T b[16];
  real_T b_y[16];
  real_T K[8];
  real_T a[8];
  real_T b_0[8];
  real_T S[4];
  real_T x_p[4];
  real_T y[2];
  real_T amp_max;
  real_T br;
  real_T dt;
  real_T phase_sine_end;
  real_T phase_zero2_end;
  real_T phase_zero_end;
  real_T u0;
  real_T u1;
  real_T u2;
  real_T x;
  real_T x3;
  real_T x4;
  int32_T TWO;
  int32_T r1;
  static const int8_T tmp_1[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1 };

  {                                    /* Sample time: [0.0s, 0.0s] */
    rate_monotonic_scheduler();
  }

  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    } else {
      simulink_experiment_debug_typ_B.HILReadEncoderTimebase =
        simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_debug_ty_DW.HILInitialize_Card,
       &simulink_experiment_debug_typ_P.HILReadAnalog_channels, 1,
       &simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }

    simulink_experiment_debug_typ_B.HILReadAnalog =
      simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_debug_typ_B.BB01SensorGainmV =
    simulink_experiment_debug_typ_P.BB01SensorGainmV_Gain *
    simulink_experiment_debug_typ_B.HILReadAnalog;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_debug_typ_B.EncoderCalibrationradcount =
    simulink_experiment_debug_typ_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_debug_typ_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_debug_typ_B.Bias =
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount +
    simulink_experiment_debug_typ_P.Bias_Bias;

  /* Clock: '<Root>/Clock' */
  simulink_experiment_debug_typ_B.Clock =
    simulink_experiment_debug_ty_M->Timing.t[0];

  /* MATLABSystem: '<Root>/MATLAB System' */
  u0 = simulink_experiment_debug_typ_B.Clock;
  u1 = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  u2 = simulink_experiment_debug_typ_B.Bias;
  obj = &simulink_experiment_debug_ty_DW.obj;

  /*         %% Main Controller Interface */
  /*  State Estimation */
  /*         %% State Estimation: Generic Time Stepping */
  obj_0 = obj;

  /*         %% State Estimation: Extended Kalman Filter */
  /*  Get some data */
  y[0] = u1;
  y[1] = u2;
  dt = u0 - obj_0->t_prev;
  x_p[0] = obj_0->x_hat[0];
  x_p[1] = obj_0->x_hat[1];
  x_p[2] = obj_0->x_hat[2];
  x_p[3] = obj_0->x_hat[3];
  if (dt > 0.0) {
    /*  Calculate kalman states */
    x3 = x_p[2];
    x4 = x_p[3];
    phase_zero_end = rt_powd_snf(x4, 4.0);
    x4 = rt_powd_snf(x4, 3.0);
    phase_sine_end = x3;
    phase_sine_end = cos(phase_sine_end);
    br = phase_sine_end * phase_sine_end;
    phase_sine_end = x3;
    phase_sine_end = cos(phase_sine_end);
    phase_zero2_end = x3;
    phase_zero2_end = sin(phase_zero2_end);
    x3 = cos(x3);
    A_lin[1] = 0.0;
    A_lin[5] = 0.0;
    A_lin[9] = 0.0051 * phase_sine_end * phase_zero2_end * phase_zero_end +
      0.4183 * x3;
    A_lin[13] = -0.0102 * x4 * br;
    A_lin[0] = 0.0;
    A_lin[2] = 0.0;
    A_lin[3] = 0.0;
    A_lin[4] = 1.0;
    A_lin[6] = 0.0;
    A_lin[7] = 0.0;
    A_lin[8] = 0.0;
    A_lin[10] = 0.0;
    A_lin[11] = 0.0;
    A_lin[12] = 0.0;
    A_lin[14] = 1.0;
    A_lin[15] = -40.0;

    /*  Predict */
    S[0] = obj_0->x_hat[0];
    S[1] = obj_0->x_hat[1];
    S[2] = obj_0->x_hat[2];
    S[3] = obj_0->x_hat[3];
    for (TWO = 0; TWO <= 2; TWO += 2) {
      tmp = _mm_loadu_pd(&A_lin[TWO]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(S[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&A_lin[TWO + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(S[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&A_lin[TWO + 8]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(S[2]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&A_lin[TWO + 12]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(S[3]));
      tmp = _mm_add_pd(tmp_0, tmp);
      _mm_storeu_pd(&x_p[TWO], tmp);
    }

    S[0] = obj_0->B[0];
    S[1] = obj_0->B[1];
    S[2] = obj_0->B[2];
    S[3] = obj_0->B[3];
    x3 = obj_0->V_servo;
    x4 = x_p[0];
    br = S[0];
    br *= x3;
    x4 += br;
    x4 *= dt;
    x_p[0] = x4;
    x4 = x_p[1];
    br = S[1];
    br *= x3;
    x4 += br;
    x4 *= dt;
    x_p[1] = x4;
    x4 = x_p[2];
    br = S[2];
    br *= x3;
    x4 += br;
    x4 *= dt;
    x_p[2] = x4;
    x4 = x_p[3];
    br = S[3];
    br *= x3;
    x4 += br;
    x4 *= dt;
    x_p[3] = x4;
    x4 = x_p[0];
    x4 += obj_0->x_hat[0];
    x_p[0] = x4;
    x4 = x_p[1];
    x4 += obj_0->x_hat[1];
    x_p[1] = x4;
    x4 = x_p[2];
    x4 += obj_0->x_hat[2];
    x_p[2] = x4;
    x4 = x_p[3];
    x4 += obj_0->x_hat[3];
    x_p[3] = x4;
    for (TWO = 0; TWO < 16; TWO++) {
      b[TWO] = obj_0->P[TWO];
    }

    for (TWO = 0; TWO < 4; TWO++) {
      for (r1 = 0; r1 <= 2; r1 += 2) {
        _mm_storeu_pd(&P_p[r1 + (TWO << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A_lin[r1]);
        tmp = _mm_mul_pd(_mm_set1_pd(b[TWO << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(b[(TWO << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[r1 + (TWO << 2)], tmp);
      }
    }

    for (TWO = 0; TWO < 16; TWO++) {
      b_y[TWO] = obj_0->P[TWO];
    }

    for (TWO = 0; TWO < 4; TWO++) {
      for (r1 = 0; r1 < 4; r1++) {
        b[TWO + (r1 << 2)] = 0.0;
        br = b[(r1 << 2) + TWO];
        br += b_y[TWO] * A_lin[r1];
        b[TWO + (r1 << 2)] = br;
        br = b[(r1 << 2) + TWO];
        br += b_y[TWO + 4] * A_lin[r1 + 4];
        b[TWO + (r1 << 2)] = br;
        br = b[(r1 << 2) + TWO];
        br += b_y[TWO + 8] * A_lin[r1 + 8];
        b[TWO + (r1 << 2)] = br;
        br = b[(r1 << 2) + TWO];
        br += b_y[TWO + 12] * A_lin[r1 + 12];
        b[TWO + (r1 << 2)] = br;
      }
    }

    for (TWO = 0; TWO < 16; TWO++) {
      x3 = P_p[TWO];
      r1 = tmp_1[TWO];
      x3 = (x3 + b[TWO]) + (real_T)r1;
      x3 *= dt;
      P_p[TWO] = x3;
    }

    for (TWO = 0; TWO < 16; TWO++) {
      x3 = P_p[TWO];
      x3 += obj_0->P[TWO];
      P_p[TWO] = x3;
    }

    for (TWO = 0; TWO < 8; TWO++) {
      a[TWO] = obj_0->C[TWO];
    }

    for (TWO = 0; TWO <= 0; TWO += 2) {
      tmp = _mm_loadu_pd(&a[TWO]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(x_p[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&a[TWO + 2]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&a[TWO + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[2]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&a[TWO + 6]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[3]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&y[TWO]);
      tmp = _mm_sub_pd(tmp_0, tmp);
      _mm_storeu_pd(&y[TWO], tmp);
    }

    for (TWO = 0; TWO < 8; TWO++) {
      b_0[TWO] = obj_0->C[TWO];
    }

    for (TWO = 0; TWO < 8; TWO++) {
      a[TWO] = obj_0->C[TWO];
    }

    for (TWO = 0; TWO < 2; TWO++) {
      for (r1 = 0; r1 < 4; r1++) {
        K[TWO + (r1 << 1)] = 0.0;
        K[TWO + (r1 << 1)] += P_p[r1 << 2] * a[TWO];
        K[TWO + (r1 << 1)] += P_p[(r1 << 2) + 1] * a[TWO + 2];
        K[TWO + (r1 << 1)] += P_p[(r1 << 2) + 2] * a[TWO + 4];
        K[TWO + (r1 << 1)] += P_p[(r1 << 2) + 3] * a[TWO + 6];
      }

      for (r1 = 0; r1 < 2; r1++) {
        S[TWO + (r1 << 1)] = 0.0;
        dt = S[(r1 << 1) + TWO];
        dt += K[TWO] * b_0[r1];
        S[TWO + (r1 << 1)] = dt;
        dt = S[(r1 << 1) + TWO];
        dt += K[TWO + 2] * b_0[r1 + 2];
        S[TWO + (r1 << 1)] = dt;
        dt = S[(r1 << 1) + TWO];
        dt += K[TWO + 4] * b_0[r1 + 4];
        S[TWO + (r1 << 1)] = dt;
        dt = S[(r1 << 1) + TWO];
        dt += K[TWO + 6] * b_0[r1 + 6];
        S[TWO + (r1 << 1)] = dt;
      }
    }

    S[0]++;
    S[3]++;
    for (TWO = 0; TWO < 8; TWO++) {
      b_0[TWO] = obj_0->C[TWO];
    }

    for (TWO = 0; TWO < 16; TWO++) {
      b_y[TWO] = obj_0->P[TWO];
    }

    for (TWO = 0; TWO < 4; TWO++) {
      for (r1 = 0; r1 < 2; r1++) {
        a[TWO + (r1 << 2)] = 0.0;
        a[TWO + (r1 << 2)] += b_y[TWO] * b_0[r1];
        a[TWO + (r1 << 2)] += b_y[TWO + 4] * b_0[r1 + 2];
        a[TWO + (r1 << 2)] += b_y[TWO + 8] * b_0[r1 + 4];
        a[TWO + (r1 << 2)] += b_y[TWO + 12] * b_0[r1 + 6];
      }
    }

    TWO = 1;
    phase_zero2_end = S[1];
    dt = fabs(phase_zero2_end);
    x3 = dt;
    phase_zero2_end = S[0];
    dt = fabs(phase_zero2_end);
    if (x3 > dt) {
      r1 = 1;
      TWO = 0;
    } else {
      r1 = 0;
    }

    dt = S[TWO] / S[r1];
    x3 = S[TWO + 2] - S[r1 + 2] * dt;
    K[r1 << 2] = a[0] / S[r1];
    K[TWO << 2] = (a[4] - K[r1 << 2] * S[r1 + 2]) / x3;
    K[r1 << 2] -= K[TWO << 2] * dt;
    K[(r1 << 2) + 1] = a[1] / S[r1];
    K[(TWO << 2) + 1] = (a[5] - K[(r1 << 2) + 1] * S[r1 + 2]) / x3;
    K[(r1 << 2) + 1] -= K[(TWO << 2) + 1] * dt;
    K[(r1 << 2) + 2] = a[2] / S[r1];
    K[(TWO << 2) + 2] = (a[6] - K[(r1 << 2) + 2] * S[r1 + 2]) / x3;
    K[(r1 << 2) + 2] -= K[(TWO << 2) + 2] * dt;
    K[(r1 << 2) + 3] = a[3] / S[r1];
    K[(TWO << 2) + 3] = (a[7] - K[(r1 << 2) + 3] * S[r1 + 2]) / x3;
    K[(r1 << 2) + 3] -= K[(TWO << 2) + 3] * dt;
    for (TWO = 0; TWO <= 2; TWO += 2) {
      tmp = _mm_loadu_pd(&K[TWO]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(y[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&K[TWO + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&x_p[TWO]);
      tmp = _mm_add_pd(tmp_0, tmp);
      _mm_storeu_pd(&x_p[TWO], tmp);
    }

    obj_0->x_hat[0] = x_p[0];
    obj_0->x_hat[1] = x_p[1];
    obj_0->x_hat[2] = x_p[2];
    obj_0->x_hat[3] = x_p[3];
    for (TWO = 0; TWO < 8; TWO++) {
      b_0[TWO] = obj_0->C[TWO];
    }

    for (TWO = 0; TWO < 4; TWO++) {
      for (r1 = 0; r1 <= 2; r1 += 2) {
        _mm_storeu_pd(&b[r1 + (TWO << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&K[r1]);
        tmp = _mm_mul_pd(_mm_set1_pd(b_0[TWO << 1]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&K[r1 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(b_0[(TWO << 1) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
      }
    }

    for (TWO = 0; TWO < 16; TWO++) {
      A_lin[TWO] = 0.0;
    }

    A_lin[0] = 1.0;
    A_lin[5] = 1.0;
    A_lin[10] = 1.0;
    A_lin[15] = 1.0;
    for (TWO = 0; TWO <= 14; TWO += 2) {
      tmp = _mm_loadu_pd(&A_lin[TWO]);
      tmp_0 = _mm_loadu_pd(&b[TWO]);
      tmp = _mm_sub_pd(tmp, tmp_0);
      _mm_storeu_pd(&A_lin[TWO], tmp);
    }

    for (TWO = 0; TWO < 4; TWO++) {
      for (r1 = 0; r1 <= 2; r1 += 2) {
        _mm_storeu_pd(&b[r1 + (TWO << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A_lin[r1]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[TWO << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(TWO << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(TWO << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[r1 + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(TWO << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&b[(TWO << 2) + r1]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&b[r1 + (TWO << 2)], tmp);
      }
    }

    for (TWO = 0; TWO < 16; TWO++) {
      obj_0->P[TWO] = b[TWO];
    }
  }

  /*  Feedback Controller */
  S[0] = obj->x_hat[0];
  S[1] = obj->x_hat[1];
  S[2] = obj->x_hat[2];
  S[3] = obj->x_hat[3];
  for (TWO = 0; TWO < 16; TWO++) {
    P_p[TWO] = obj->P[TWO];
  }

  obj_0 = obj;

  /*         %% Feedback Controller: P Controller */
  br = x_p[0];
  x4 = x_p[1];
  dt = x_p[2];
  x3 = x_p[3];
  if (u0 < 5.0) {
    phase_sine_end = 0.0;
    amp_max = 0.0;
  } else if (u0 < 61.85) {
    phase_sine_end = u0 - 5.0;
    phase_zero2_end = phase_sine_end / 56.85;
    if (phase_zero2_end < 0.5) {
      phase_zero_end = phase_zero2_end / 0.5 * 0.090000000000000011 + 0.05;
      phase_zero2_end = 0.11423973285781065 * phase_sine_end;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 * phase_sine_end -
        0.2094395102393195 * phase_zero2_end / 0.11423973285781065;
      phase_zero2_end = sin(phase_zero2_end);
      amp_max = 0.11423973285781065 * phase_sine_end;
      amp_max = sin(amp_max);
      amp_max = 0.83775804095727813 * phase_sine_end - 0.2094395102393195 *
        amp_max / 0.11423973285781065;
      amp_max = cos(amp_max);
      x = 0.11423973285781065 * phase_sine_end;
      x = cos(x);
      amp_max = (0.83775804095727813 - 0.2094395102393195 * x) * (phase_zero_end
        * amp_max) + 0.00316622691292876 * phase_zero2_end;
    } else {
      phase_zero_end = 0.14;
      phase_zero2_end = 0.11423973285781065 * phase_sine_end;
      phase_zero2_end = sin(phase_zero2_end);
      phase_zero2_end = 0.83775804095727813 * phase_sine_end -
        0.2094395102393195 * phase_zero2_end / 0.11423973285781065;
      phase_zero2_end = cos(phase_zero2_end);
      amp_max = 0.11423973285781065 * phase_sine_end;
      amp_max = cos(amp_max);
      amp_max = (0.83775804095727813 - 0.2094395102393195 * amp_max) * (0.14 *
        phase_zero2_end);
    }

    phase_zero2_end = 0.11423973285781065 * phase_sine_end;
    phase_zero2_end = sin(phase_zero2_end);
    phase_zero2_end = 0.83775804095727813 * phase_sine_end - 0.2094395102393195 *
      phase_zero2_end / 0.11423973285781065;
    phase_zero2_end = sin(phase_zero2_end);
    phase_sine_end = phase_zero_end * phase_zero2_end;
  } else if (u0 < 65.0) {
    phase_sine_end = 0.0;
    amp_max = 0.0;
  } else if (u0 < 85.0) {
    phase_zero_end = u0 - 65.0;
    phase_sine_end = phase_zero_end / 20.0;
    if (phase_sine_end < 0.5) {
      phase_sine_end = 0.05;
    } else {
      phase_sine_end = 0.1;
    }

    phase_zero2_end = 0.62831853071795862 * phase_zero_end;
    phase_zero2_end = sin(phase_zero2_end);
    if (phase_zero2_end < 0.0) {
      phase_zero2_end = -1.0;
    } else {
      phase_zero2_end = (phase_zero2_end > 0.0);
    }

    phase_sine_end *= phase_zero2_end;
    amp_max = 0.0;
  } else {
    phase_sine_end = 0.0;
    amp_max = 0.0;
  }

  br = (br - phase_sine_end) * -6.0;
  x4 = (x4 - amp_max) * -10.0;
  if (!(br <= 0.97738438111682457)) {
    br = 0.97738438111682457;
  }

  if (!(br >= -0.97738438111682457)) {
    br = -0.97738438111682457;
  }

  x3 = (br - dt) * 7.0 + (x4 - x3) * 8.0;
  obj_0->t_prev = u0;

  /*  V_servo = stepImplLQR(obj, t, xg); */
  /*  V_servo = stepImplLQG(obj, t, xg); */
  x4 = x3 - obj->V_servo;
  dt = fabs(x4);
  if (dt > 0.04) {
    if (rtIsNaN(x4)) {
      x4 = (rtNaN);
    } else if (x4 < 0.0) {
      x4 = -1.0;
    } else {
      x4 = (x4 > 0.0);
    }

    x3 = 0.04 * x4 + obj->V_servo;
  } else {
    dt = fabs(x4);
    if (dt > 0.0025) {
      x3 = obj->V_servo;
    }
  }

  if (x3 > 5.0) {
    x3 = 5.0;
  } else if (x3 < -5.0) {
    x3 = -5.0;
  }

  obj->V_servo = x3;
  obj->t_prev = u0;
  obj->p_prev = u1;
  obj->theta_prev = u2;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o1 = x3;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o2[0] = S[0];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[1] = S[1];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[2] = S[2];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[3] = S[3];

  /* MATLABSystem: '<Root>/MATLAB System' */
  memcpy(&simulink_experiment_debug_typ_B.MATLABSystem_o3[0], &P_p[0], sizeof
         (real_T) << 4U);

  /* Saturate: '<Root>/+//-10V' */
  u0 = simulink_experiment_debug_typ_B.MATLABSystem_o1;
  u1 = simulink_experiment_debug_typ_P.u0V_LowerSat;
  u2 = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (u0 > u2) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u2;
  } else if (u0 < u1) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u1;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u0;
  }

  /* End of Saturate: '<Root>/+//-10V' */

  /* Gain: '<S1>/Motor  Gain (V//V)' */
  simulink_experiment_debug_typ_B.MotorGainVV =
    simulink_experiment_debug_typ_P.MotorGainVV_Gain *
    simulink_experiment_debug_typ_B.u0V;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_debug_ty_DW.HILInitialize_Card,
      &simulink_experiment_debug_typ_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_debug_typ_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
    }
  }

  /* MATLAB Function: '<Root>/MATLAB Function' */
  /* MATLAB Function 'MATLAB Function': '<S2>:1' */
  /* '<S2>:1:3' */
  if (simulink_experiment_debug_typ_B.Clock < 5.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 61.85) {
    phase_zero2_end = (simulink_experiment_debug_typ_B.Clock - 5.0) / 56.85;
    if (phase_zero2_end < 0.5) {
      phase_zero_end = phase_zero2_end / 0.5 * 0.090000000000000011 + 0.05;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * phase_zero_end *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195) + sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.00316622691292876;
      u0 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = (cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 12.0 * (0.83775804095727813 - cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 3.1415926535897931 / 15.0) / 1895.0 + sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * ((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.0 / 1895.0 + 0.05) * (u0 * u0)) + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 19.739208802178716) *
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.0 / 1895.0 + 0.05) /
        825.0;
    } else {
      phase_zero_end = 0.14;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.14 *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195);
      phase_sine_end = 0.83775804095727813 - cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 7.0 * (phase_sine_end * phase_sine_end) /
        50.0 + cos(sin((simulink_experiment_debug_typ_B.Clock - 5.0) *
                       6.2831853071795862 / 55.0) * 11.0 / 6.0 -
                   (simulink_experiment_debug_typ_B.Clock - 5.0) *
                   12.566370614359172 / 15.0) * (sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 69.0872308076255) / 20625.0;
    }

    simulink_experiment_debug_typ_B.p_ref = sin
      ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 - sin
       ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065) *
       0.2094395102393195 / 0.11423973285781065) * phase_zero_end;
  } else if (simulink_experiment_debug_typ_B.Clock < 65.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 85.0) {
    if ((simulink_experiment_debug_typ_B.Clock - 65.0) / 20.0 < 0.5) {
      phase_sine_end = 0.05;
    } else {
      phase_sine_end = 0.1;
    }

    u0 = sin((simulink_experiment_debug_typ_B.Clock - 65.0) *
             0.62831853071795862);
    if (rtIsNaN(u0)) {
      u0 = (rtNaN);
    } else if (u0 < 0.0) {
      u0 = -1.0;
    } else {
      u0 = (u0 > 0.0);
    }

    simulink_experiment_debug_typ_B.p_ref = phase_sine_end * u0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* End of MATLAB Function: '<Root>/MATLAB Function' */

  /* Gain: '<Root>/m to cm' */
  /* '<S2>:1:3' */
  simulink_experiment_debug_typ_B.mtocm[0] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.p_ref;
  simulink_experiment_debug_typ_B.mtocm[1] =
    simulink_experiment_debug_typ_P.mtocm_Gain *
    simulink_experiment_debug_typ_B.BB01SensorGainmV;

  /* Gain: '<S3>/Gain' */
  simulink_experiment_debug_typ_B.Gain =
    simulink_experiment_debug_typ_P.Gain_Gain *
    simulink_experiment_debug_typ_B.Bias;

  /* RateTransition: '<Root>/Rate Transition' */
  if (simulink_experiment_debug_ty_M->Timing.RateInteraction.TID1_2) {
    simulink_experiment_debug_ty_DW.RateTransition_Buffer =
      simulink_experiment_debug_typ_B.Clock;

    /* RateTransition: '<Root>/Rate Transition1' */
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer =
      simulink_experiment_debug_typ_B.p_ref;

    /* RateTransition: '<Root>/Rate Transition2' */
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer =
      simulink_experiment_debug_typ_B.MATLABSystem_o1;

    /* RateTransition: '<Root>/Rate Transition3' */
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer =
      simulink_experiment_debug_typ_B.BB01SensorGainmV;

    /* RateTransition: '<Root>/Rate Transition4' */
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer =
      simulink_experiment_debug_typ_B.Bias;
  }

  /* End of RateTransition: '<Root>/Rate Transition' */
}

/* Model update function for TID0 */
void simulink_experiment_debug_type1_update0(void) /* Sample time: [0.0s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick0)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH0;
  }

  simulink_experiment_debug_ty_M->Timing.t[0] =
    simulink_experiment_debug_ty_M->Timing.clockTick0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 +
    simulink_experiment_debug_ty_M->Timing.clockTickH0 *
    simulink_experiment_debug_ty_M->Timing.stepSize0 * 4294967296.0;

  /* Update absolute time */
  /* The "clockTick1" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick1"
   * and "Timing.stepSize1". Size of "clockTick1" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick1 and the high bits
   * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick1)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH1;
  }

  simulink_experiment_debug_ty_M->Timing.t[1] =
    simulink_experiment_debug_ty_M->Timing.clockTick1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 +
    simulink_experiment_debug_ty_M->Timing.clockTickH1 *
    simulink_experiment_debug_ty_M->Timing.stepSize1 * 4294967296.0;
}

/* Model output function for TID2 */
void simulink_experiment_debug_type1_output2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* RateTransition: '<Root>/Rate Transition2' */
  simulink_experiment_debug_typ_B.RateTransition2 =
    simulink_experiment_debug_ty_DW.RateTransition2_Buffer;

  /* RateTransition: '<Root>/Rate Transition1' */
  simulink_experiment_debug_typ_B.RateTransition1 =
    simulink_experiment_debug_ty_DW.RateTransition1_Buffer;

  /* RateTransition: '<Root>/Rate Transition3' */
  simulink_experiment_debug_typ_B.RateTransition3 =
    simulink_experiment_debug_ty_DW.RateTransition3_Buffer;

  /* RateTransition: '<Root>/Rate Transition4' */
  simulink_experiment_debug_typ_B.RateTransition4 =
    simulink_experiment_debug_ty_DW.RateTransition4_Buffer;

  /* RateTransition: '<Root>/Rate Transition' */
  simulink_experiment_debug_typ_B.RateTransition =
    simulink_experiment_debug_ty_DW.RateTransition_Buffer;
}

/* Model update function for TID2 */
void simulink_experiment_debug_type1_update2(void) /* Sample time: [0.01s, 0.0s] */
{
  /* Update absolute time */
  /* The "clockTick2" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick2"
   * and "Timing.stepSize2". Size of "clockTick2" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick2 and the high bits
   * Timing.clockTickH2. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_debug_ty_M->Timing.clockTick2)) {
    ++simulink_experiment_debug_ty_M->Timing.clockTickH2;
  }

  simulink_experiment_debug_ty_M->Timing.t[2] =
    simulink_experiment_debug_ty_M->Timing.clockTick2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 +
    simulink_experiment_debug_ty_M->Timing.clockTickH2 *
    simulink_experiment_debug_ty_M->Timing.stepSize2 * 4294967296.0;
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_output(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_output0();
    break;

   case 2 :
    simulink_experiment_debug_type1_output2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Use this function only if you need to maintain compatibility with an existing static main program. */
void simulink_experiment_debug_type1_update(int_T tid)
{
  switch (tid) {
   case 0 :
    simulink_experiment_debug_type1_update0();
    break;

   case 2 :
    simulink_experiment_debug_type1_update2();
    break;

   default :
    /* do nothing */
    break;
  }
}

/* Model initialize function */
void simulink_experiment_debug_type1_initialize(void)
{
  {
    studentControllerInterface_si_T *b_obj;
    int32_T i;
    static const int8_T tmp[8] = { 1, 0, 0, 0, 0, 1, 0, 0 };

    static const real_T tmp_0[16] = { 0.00467317845216159, 0.00447778693210613,
      1.79600044138968E-6, -3.3479774974499E-9, 0.00447778693304056,
      0.104399968275514, 8.18335148621331E-5, -3.66963222880958E-8,
      1.79600514031296E-6, 8.18335789697825E-5, 0.00447750476234807,
      2.00500225883704E-5, -3.34363679122678E-9, -3.66928669353362E-8,
      2.00500225891637E-5, 0.00124997492495667 };

    /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
    {
      t_int result;
      t_boolean is_switching;
      result = hil_open("q2_usb", "0",
                        &simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      is_switching = false;
      result = hil_set_card_specific_options
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      result = hil_watchdog_clear
        (simulink_experiment_debug_ty_DW.HILInitialize_Card);
      if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AILow);
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AIHigh;
        result = hil_set_analog_input_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] =
          (simulink_experiment_debug_typ_P.HILInitialize_AOLow);
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOHigh;
        result = hil_set_analog_output_ranges
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0],
           &simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_AOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_AOEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOInitial;
        result = hil_write_analog
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_AOReset) {
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
          simulink_experiment_debug_typ_P.HILInitialize_AOWatchdog;
        result = hil_watchdog_set_analog_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_AOChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      result = hil_set_digital_directions
        (simulink_experiment_debug_ty_DW.HILInitialize_Card, NULL, 0U,
         simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        return;
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_DOStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_DOEnter && is_switching))
      {
        {
          int_T i1;
          boolean_T *dw_DOBits =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOBits[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOInitial;
          }
        }

        result = hil_write_digital
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U,
           (t_boolean *) &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if (simulink_experiment_debug_typ_P.HILInitialize_DOReset) {
        {
          int_T i1;
          int32_T *dw_DOStates =
            &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0];
          for (i1=0; i1 < 8; i1++) {
            dw_DOStates[i1] =
              simulink_experiment_debug_typ_P.HILInitialize_DOWatchdog;
          }
        }

        result = hil_watchdog_set_digital_expiration_state
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_DOChannels, 8U, (const
            t_digital_state *)
           &simulink_experiment_debug_ty_DW.HILInitialize_DOStates[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIPStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIPEnter &&
           is_switching)) {
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIQuadrature;
        result = hil_set_encoder_quadrature_mode
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           (t_encoder_quadrature_mode *)
           &simulink_experiment_debug_ty_DW.HILInitialize_QuadratureModes[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }

      if ((simulink_experiment_debug_typ_P.HILInitialize_EIStart &&
           !is_switching) ||
          (simulink_experiment_debug_typ_P.HILInitialize_EIEnter && is_switching))
      {
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[1] =
          simulink_experiment_debug_typ_P.HILInitialize_EIInitial;
        result = hil_set_encoder_counts
          (simulink_experiment_debug_ty_DW.HILInitialize_Card,
           simulink_experiment_debug_typ_P.HILInitialize_EIChannels, 2U,
           &simulink_experiment_debug_ty_DW.HILInitialize_InitialEICounts[0]);
        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
          return;
        }
      }
    }

    /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

    /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
    {
      t_error result;
      result = hil_task_create_encoder_reader
        (simulink_experiment_debug_ty_DW.HILInitialize_Card,
         simulink_experiment_debug_typ_P.HILReadEncoderTimebase_SamplesI,
         &simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Channels, 1,
         &simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task);
      if (result >= 0) {
        result = hil_task_set_buffer_overflow_mode
          (simulink_experiment_debug_ty_DW.HILReadEncoderTimebase_Task,
           (t_buffer_overflow_mode)
           (simulink_experiment_debug_typ_P.HILReadEncoderTimebase_Overflow - 1));
      }

      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
      }
    }

    /* Start for MATLABSystem: '<Root>/MATLAB System' */
    b_obj = &simulink_experiment_debug_ty_DW.obj;
    b_obj->t_prev = -1.0;
    b_obj->p_prev = 0.0;
    b_obj->theta_prev = 0.0;
    b_obj->V_servo = 0.0;
    b_obj->B[0] = 0.0;
    b_obj->B[1] = 0.0;
    b_obj->B[2] = 0.0;
    b_obj->B[3] = 60.0;
    for (i = 0; i < 8; i++) {
      b_obj->C[i] = tmp[i];
    }

    b_obj->x_hat[0] = -0.05;
    b_obj->x_hat[1] = 0.0;
    b_obj->x_hat[2] = 0.0;
    b_obj->x_hat[3] = 0.0;
    for (i = 0; i < 16; i++) {
      b_obj->P[i] = tmp_0[i];
    }

    simulink_experiment_debug_ty_DW.objisempty = true;

    /* End of Start for MATLABSystem: '<Root>/MATLAB System' */
  }
}

/* Model terminate function */
void simulink_experiment_debug_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_debug_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_debug_typ_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_debug_typ_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_debug_typ_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_debug_typ_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_debug_typ_P.HILInitialize_DOFinal;
        }
      }

      num_final_digital_outputs = 8U;
    } else {
      num_final_digital_outputs = 0;
    }

    if (0
        || num_final_analog_outputs > 0
        || num_final_digital_outputs > 0
        ) {
      /* Attempt to write the final outputs atomically (due to firmware issue in old Q2-USB). Otherwise write channels individually */
      result = hil_write(simulink_experiment_debug_ty_DW.HILInitialize_Card
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_debug_ty_DW.HILInitialize_Card,
             simulink_experiment_debug_typ_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_debug_ty_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_debug_ty_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    hil_close(simulink_experiment_debug_ty_DW.HILInitialize_Card);
    simulink_experiment_debug_ty_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_output(tid);
}

void MdlUpdate(int_T tid)
{
  if (tid == 1)
    tid = 0;
  simulink_experiment_debug_type1_update(tid);
}

void MdlInitializeSizes(void)
{
}

void MdlInitializeSampleTimes(void)
{
}

void MdlInitialize(void)
{
}

void MdlStart(void)
{
  simulink_experiment_debug_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_debug_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_debug_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_debug_ty_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          &simulink_experiment_debug_ty_M->Timing.simTimeStep);
    rtsiSetTPtr(&simulink_experiment_debug_ty_M->solverInfo, &rtmGetTPtr
                (simulink_experiment_debug_ty_M));
    rtsiSetStepSizePtr(&simulink_experiment_debug_ty_M->solverInfo,
                       &simulink_experiment_debug_ty_M->Timing.stepSize0);
    rtsiSetErrorStatusPtr(&simulink_experiment_debug_ty_M->solverInfo,
                          (&rtmGetErrorStatus(simulink_experiment_debug_ty_M)));
    rtsiSetRTModelPtr(&simulink_experiment_debug_ty_M->solverInfo,
                      simulink_experiment_debug_ty_M);
  }

  rtsiSetSimTimeStep(&simulink_experiment_debug_ty_M->solverInfo,
                     MAJOR_TIME_STEP);
  rtsiSetIsMinorTimeStepWithModeChange
    (&simulink_experiment_debug_ty_M->solverInfo, false);
  rtsiSetSolverName(&simulink_experiment_debug_ty_M->solverInfo,
                    "FixedStepDiscrete");

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;
    mdlTsMap[1] = 1;
    mdlTsMap[2] = 2;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_debug_ty_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_debug_ty_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_debug_ty_M->Timing.sampleTimes =
      (&simulink_experiment_debug_ty_M->Timing.sampleTimesArray[0]);
    simulink_experiment_debug_ty_M->Timing.offsetTimes =
      (&simulink_experiment_debug_ty_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_debug_ty_M->Timing.sampleTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[1] = (0.002);
    simulink_experiment_debug_ty_M->Timing.sampleTimes[2] = (0.01);

    /* task offsets */
    simulink_experiment_debug_ty_M->Timing.offsetTimes[0] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[1] = (0.0);
    simulink_experiment_debug_ty_M->Timing.offsetTimes[2] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_debug_ty_M,
             &simulink_experiment_debug_ty_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_debug_ty_M->Timing.sampleHitArray;
    int_T *mdlPerTaskSampleHits =
      simulink_experiment_debug_ty_M->Timing.perTaskSampleHitsArray;
    simulink_experiment_debug_ty_M->Timing.perTaskSampleHits =
      (&mdlPerTaskSampleHits[0]);
    mdlSampleHits[0] = 1;
    simulink_experiment_debug_ty_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_debug_ty_M, 20.0);
  simulink_experiment_debug_ty_M->Timing.stepSize0 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize1 = 0.002;
  simulink_experiment_debug_ty_M->Timing.stepSize2 = 0.01;

  /* External mode info */
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (1816959090U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (1372178900U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (4014209278U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (1382840546U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[3];
    simulink_experiment_debug_ty_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    systemRan[1] = &rtAlwaysEnabled;
    systemRan[2] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_debug_ty_M->extModeInfo,
      &simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_debug_ty_M->extModeInfo,
                        simulink_experiment_debug_ty_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_debug_ty_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_debug_ty_M));
  }

  simulink_experiment_debug_ty_M->solverInfoPtr =
    (&simulink_experiment_debug_ty_M->solverInfo);
  simulink_experiment_debug_ty_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_debug_ty_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_debug_ty_M->solverInfo,
                    SOLVER_MODE_MULTITASKING);

  /* block I/O */
  simulink_experiment_debug_ty_M->blockIO = ((void *)
    &simulink_experiment_debug_typ_B);

  {
    int32_T i;
    for (i = 0; i < 16; i++) {
      simulink_experiment_debug_typ_B.MATLABSystem_o3[i] = 0.0;
    }

    simulink_experiment_debug_typ_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_debug_typ_B.HILReadAnalog = 0.0;
    simulink_experiment_debug_typ_B.BB01SensorGainmV = 0.0;
    simulink_experiment_debug_typ_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_debug_typ_B.Bias = 0.0;
    simulink_experiment_debug_typ_B.Clock = 0.0;
    simulink_experiment_debug_typ_B.u0V = 0.0;
    simulink_experiment_debug_typ_B.MotorGainVV = 0.0;
    simulink_experiment_debug_typ_B.mtocm[0] = 0.0;
    simulink_experiment_debug_typ_B.mtocm[1] = 0.0;
    simulink_experiment_debug_typ_B.Gain = 0.0;
    simulink_experiment_debug_typ_B.RateTransition2 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition1 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition3 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition4 = 0.0;
    simulink_experiment_debug_typ_B.RateTransition = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o1 = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o2[0] = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o2[1] = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o2[2] = 0.0;
    simulink_experiment_debug_typ_B.MATLABSystem_o2[3] = 0.0;
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  }

  /* parameters */
  simulink_experiment_debug_ty_M->defaultParam = ((real_T *)
    &simulink_experiment_debug_typ_P);

  /* states (dwork) */
  simulink_experiment_debug_ty_M->dwork = ((void *)
    &simulink_experiment_debug_ty_DW);
  (void) memset((void *)&simulink_experiment_debug_ty_DW, 0,
                sizeof(DW_simulink_experiment_debug__T));
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_debug_ty_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_debug_ty_DW.HILReadAnalog_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition1_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition2_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition3_Buffer = 0.0;
  simulink_experiment_debug_ty_DW.RateTransition4_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_debug_ty_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 22;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_debug_ty_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_debug_ty_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_debug_ty_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_debug_ty_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_debug_ty_M->Sizes.numSampTimes = (3);/* Number of sample times */
  simulink_experiment_debug_ty_M->Sizes.numBlocks = (31);/* Number of blocks */
  simulink_experiment_debug_ty_M->Sizes.numBlockIO = (21);/* Number of block outputs */
  simulink_experiment_debug_ty_M->Sizes.numBlockPrms = (88);/* Sum of parameter "widths" */
  return simulink_experiment_debug_ty_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
