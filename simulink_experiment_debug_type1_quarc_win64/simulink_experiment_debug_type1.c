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
 * C source code generated on : Wed Apr 30 12:53:37 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "simulink_experiment_debug_type1.h"
#include "simulink_experiment_debug_type1_types.h"
#include "rtwtypes.h"
#include <string.h>
#include <math.h>
#include "rt_nonfinite.h"
#include "simulink_experiment_debug_type1_private.h"
#include <emmintrin.h>
#include "simulink_experiment_debug_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_debug_t_T simulink_experiment_debug_typ_B;

/* Block states (default storage) */
DW_simulink_experiment_debug__T simulink_experiment_debug_ty_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_debug_ty_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_debug_ty_M =
  &simulink_experiment_debug_ty_M_;

/* Forward declaration for local functions */
static void simulink_experiment_debug_mrdiv(const real_T A[8], const real_T B[4],
  real_T Y[8]);
static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[32],
  int32_T ix0);
static void simulink_experiment_debu_xgeqp3(const real_T A[32], real_T b_A[32],
  real_T tau[4], int32_T jpvt[4]);
static void simulink_experiment_de_mldivide(const real_T A[32], const real_T B
  [32], real_T Y[16]);
static void studentControllerInterface_step(studentControllerInterface_si_T *obj,
  real_T t, real_T p_ball, real_T theta, real_T *V_servo, real_T x_p[4], real_T
  P_p[16]);
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

static void simulink_experiment_debug_mrdiv(const real_T A[8], const real_T B[4],
  real_T Y[8])
{
  real_T B_0[8];
  real_T A_0[4];
  real_T a21;
  real_T a22;
  int32_T TWO;
  int32_T r1;
  A_0[0] = B[0];
  A_0[1] = B[1];
  A_0[2] = B[2];
  A_0[3] = B[3];
  memcpy(&B_0[0], &A[0], sizeof(real_T) << 3U);
  TWO = 1;
  a22 = A_0[1];
  a22 = fabs(a22);
  a21 = a22;
  a22 = A_0[0];
  a22 = fabs(a22);
  if (a21 > a22) {
    r1 = 1;
    TWO = 0;
  } else {
    r1 = 0;
  }

  a21 = A_0[TWO] / A_0[r1];
  a22 = A_0[TWO + 2] - A_0[r1 + 2] * a21;
  Y[r1 << 2] = B_0[0] / A_0[r1];
  Y[TWO << 2] = (B_0[4] - Y[r1 << 2] * A_0[r1 + 2]) / a22;
  Y[r1 << 2] -= Y[TWO << 2] * a21;
  Y[(r1 << 2) + 1] = B_0[1] / A_0[r1];
  Y[(TWO << 2) + 1] = (B_0[5] - Y[(r1 << 2) + 1] * A_0[r1 + 2]) / a22;
  Y[(r1 << 2) + 1] -= Y[(TWO << 2) + 1] * a21;
  Y[(r1 << 2) + 2] = B_0[2] / A_0[r1];
  Y[(TWO << 2) + 2] = (B_0[6] - Y[(r1 << 2) + 2] * A_0[r1 + 2]) / a22;
  Y[(r1 << 2) + 2] -= Y[(TWO << 2) + 2] * a21;
  Y[(r1 << 2) + 3] = B_0[3] / A_0[r1];
  Y[(TWO << 2) + 3] = (B_0[7] - Y[(r1 << 2) + 3] * A_0[r1 + 2]) / a22;
  Y[(r1 << 2) + 3] -= Y[(TWO << 2) + 3] * a21;
}

static real_T simulink_experiment_debug_xnrm2(int32_T n, const real_T x[32],
  int32_T ix0)
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  int32_T b;
  int32_T kend;
  y = 0.0;
  if (n < 1) {
  } else if (n == 1) {
    absxk = x[ix0 - 1];
    y = fabs(absxk);
  } else {
    scale = 3.3121686421112381E-170;
    kend = n - 1;
    kend += ix0;
    for (b = ix0; b <= kend; b++) {
      absxk = x[b - 1];
      absxk = fabs(absxk);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }
    }

    y = scale * sqrt(y);
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T a;
  real_T b;
  real_T y;
  a = fabs(u0);
  b = fabs(u1);
  if (a < b) {
    a /= b;
    y = sqrt(a * a + 1.0) * b;
  } else if (a > b) {
    b /= a;
    y = sqrt(b * b + 1.0) * a;
  } else if (rtIsNaN(b)) {
    y = (rtNaN);
  } else {
    y = a * 1.4142135623730951;
  }

  return y;
}

static void simulink_experiment_debu_xgeqp3(const real_T A[32], real_T b_A[32],
  real_T tau[4], int32_T jpvt[4])
{
  __m128d tmp;
  real_T A_0[32];
  real_T x[32];
  real_T vn1[4];
  real_T vn2[4];
  real_T work[4];
  real_T RSAFMIN;
  real_T absxk;
  real_T scale;
  real_T t;
  real_T zero;
  int32_T ONE;
  int32_T exitg1;
  int32_T i;
  int32_T ia;
  int32_T itemp;
  int32_T ix;
  int32_T jA;
  int32_T kend;
  int32_T lda;
  int32_T mmi;
  int32_T mmip1;
  int32_T nfxd;
  int32_T nm1;
  int32_T nmi;
  boolean_T exitg2;
  work[0] = 0.0;
  work[1] = 0.0;
  work[2] = 0.0;
  work[3] = 0.0;
  for (i = 0; i < 32; i++) {
    b_A[i] = A[i];
    x[i] = b_A[i];
  }

  for (lda = 0; lda < 4; lda++) {
    i = lda;
    jpvt[i] = i + 1;
    ix = (lda << 3) + 1;
    zero = 0.0;
    scale = 3.3121686421112381E-170;
    kend = ix + 7;
    for (i = ix; i <= kend; i++) {
      absxk = x[i - 1];
      absxk = fabs(absxk);
      if (absxk > scale) {
        t = scale / absxk;
        zero = zero * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        zero += t * t;
      }
    }

    zero = scale * sqrt(zero);
    vn1[lda] = zero;
    vn2[lda] = vn1[lda];
  }

  for (nfxd = 0; nfxd < 4; nfxd++) {
    kend = ((nfxd << 3) + nfxd) + 1;
    nmi = 3 - nfxd;
    mmi = 7 - nfxd;
    mmip1 = mmi + 1;
    lda = nmi + 1;
    ONE = 1;
    if (lda > 1) {
      ix = nfxd;
      absxk = vn1[nfxd];
      scale = absxk;
      zero = scale;
      for (i = 2; i <= lda; i++) {
        ix++;
        absxk = vn1[ix];
        if (absxk > zero) {
          ONE = i;
          zero = absxk;
        }
      }
    }

    lda = (nfxd + ONE) - 1;
    if (lda != nfxd) {
      ix = lda << 3;
      ONE = nfxd << 3;
      for (i = 0; i < 8; i++) {
        absxk = b_A[ix];
        b_A[ix] = b_A[ONE];
        b_A[ONE] = absxk;
        ix++;
        ONE++;
      }

      itemp = jpvt[lda];
      jpvt[lda] = jpvt[nfxd];
      jpvt[nfxd] = itemp;
      vn1[lda] = vn1[nfxd];
      vn2[lda] = vn2[nfxd];
    }

    zero = b_A[kend - 1];
    ix = kend + 1;
    absxk = 0.0;
    nm1 = mmip1 - 1;
    scale = simulink_experiment_debug_xnrm2(nm1, b_A, ix);
    if (scale != 0.0) {
      t = rt_hypotd_snf(zero, scale);
      if (zero >= 0.0) {
        t = -t;
      }

      scale = fabs(t);
      if (scale < 1.0020841800044864E-292) {
        lda = -1;
        do {
          lda++;
          i = nm1 - 1;
          itemp = ix + i;
          ONE = ((((itemp - ix) + 1) / 2) << 1) + ix;
          jA = ONE - 2;
          for (i = ix; i <= jA; i += 2) {
            tmp = _mm_loadu_pd(&b_A[i - 1]);
            tmp = _mm_mul_pd(tmp, _mm_set1_pd(9.9792015476736E+291));
            _mm_storeu_pd(&b_A[i - 1], tmp);
          }

          for (i = ONE; i <= itemp; i++) {
            b_A[i - 1] *= 9.9792015476736E+291;
          }

          t *= 9.9792015476736E+291;
          zero *= 9.9792015476736E+291;
          scale = fabs(t);
        } while ((scale < 1.0020841800044864E-292) && (lda + 1 < 20));

        scale = simulink_experiment_debug_xnrm2(nm1, b_A, ix);
        t = rt_hypotd_snf(zero, scale);
        if (zero >= 0.0) {
          t = -t;
        }

        absxk = t - zero;
        absxk /= t;
        scale = zero - t;
        zero = 1.0 / scale;
        i = nm1 - 1;
        itemp = ix + i;
        ONE = ((((itemp - ix) + 1) / 2) << 1) + ix;
        jA = ONE - 2;
        for (i = ix; i <= jA; i += 2) {
          tmp = _mm_loadu_pd(&b_A[i - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(zero));
          _mm_storeu_pd(&b_A[i - 1], tmp);
        }

        for (i = ONE; i <= itemp; i++) {
          b_A[i - 1] *= zero;
        }

        for (i = 0; i <= lda; i++) {
          t *= 1.0020841800044864E-292;
        }

        zero = t;
      } else {
        absxk = t - zero;
        absxk /= t;
        scale = zero - t;
        zero = 1.0 / scale;
        i = nm1 - 1;
        itemp = ix + i;
        ONE = ((((itemp - ix) + 1) / 2) << 1) + ix;
        jA = ONE - 2;
        for (i = ix; i <= jA; i += 2) {
          tmp = _mm_loadu_pd(&b_A[i - 1]);
          tmp = _mm_mul_pd(tmp, _mm_set1_pd(zero));
          _mm_storeu_pd(&b_A[i - 1], tmp);
        }

        for (i = ONE; i <= itemp; i++) {
          b_A[i - 1] *= zero;
        }

        zero = t;
      }
    }

    tau[nfxd] = absxk;
    b_A[kend - 1] = zero;
    if (nfxd + 1 < 4) {
      zero = b_A[kend - 1];
      b_A[kend - 1] = 1.0;
      t = tau[nfxd];
      lda = kend + 8;
      if (t != 0.0) {
        i = mmip1 - 1;
        i += kend;
        while ((mmip1 > 0) && (b_A[i - 1] == 0.0)) {
          mmip1--;
          i--;
        }

        exitg2 = false;
        while ((!exitg2) && (nmi > 0)) {
          i = nmi - 1;
          i <<= 3;
          ONE = lda + i;
          i = mmip1 - 1;
          i += ONE;
          do {
            exitg1 = 0;
            if (ONE <= i) {
              if (b_A[ONE - 1] != 0.0) {
                exitg1 = 1;
              } else {
                ONE++;
              }
            } else {
              nmi--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        mmip1 = 0;
        nmi = 0;
      }

      if (mmip1 > 0) {
        for (i = 0; i < 32; i++) {
          absxk = b_A[i];
          x[i] = absxk;
          A_0[i] = absxk;
        }

        if (nmi != 0) {
          ONE = mmip1 - 1;
          nm1 = nmi - 1;
          i = nm1;
          itemp = i;
          for (i = 0; i <= itemp; i++) {
            jA = i;
            work[jA] = 0.0;
          }

          jA = 0;
          i = nm1 << 3;
          itemp = lda + i;
          for (i = lda; i <= itemp; i += 8) {
            ix = kend - 1;
            RSAFMIN = 0.0;
            nm1 = i + ONE;
            for (ia = i; ia <= nm1; ia++) {
              scale = A_0[ia - 1];
              absxk = x[ix];
              absxk *= scale;
              RSAFMIN += absxk;
              ix++;
            }

            work[jA] += RSAFMIN;
            jA++;
          }
        }

        t = -t;
        ONE = 0;
        if (!(t == 0.0)) {
          jA = lda - 1;
          for (lda = 0; lda < nmi; lda++) {
            absxk = work[ONE];
            if (absxk != 0.0) {
              absxk *= t;
              ix = kend - 1;
              itemp = jA + 1;
              i = mmip1 + jA;
              for (nm1 = itemp; nm1 <= i; nm1++) {
                b_A[nm1 - 1] += b_A[ix] * absxk;
                ix++;
              }
            }

            ONE++;
            jA += 8;
          }
        }
      }

      b_A[kend - 1] = zero;
    }

    for (nmi = nfxd + 2; nmi < 5; nmi++) {
      kend = ((nmi - 1) << 3) + nfxd;
      if (vn1[nmi - 1] != 0.0) {
        absxk = b_A[kend];
        scale = fabs(absxk);
        absxk = scale / vn1[nmi - 1];
        absxk = 1.0 - absxk * absxk;
        if (absxk < 0.0) {
          absxk = 0.0;
        }

        zero = vn1[nmi - 1] / vn2[nmi - 1];
        zero = zero * zero * absxk;
        if (zero <= 1.4901161193847656E-8) {
          ix = kend + 2;
          memcpy(&x[0], &b_A[0], sizeof(real_T) << 5U);
          zero = 0.0;
          scale = 3.3121686421112381E-170;
          i = mmi - 1;
          kend = ix + i;
          for (i = ix; i <= kend; i++) {
            absxk = x[i - 1];
            absxk = fabs(absxk);
            if (absxk > scale) {
              t = scale / absxk;
              zero = zero * t * t + 1.0;
              scale = absxk;
            } else {
              t = absxk / scale;
              zero += t * t;
            }
          }

          zero = scale * sqrt(zero);
          vn1[nmi - 1] = zero;
          vn2[nmi - 1] = vn1[nmi - 1];
        } else {
          absxk = sqrt(absxk);
          vn1[nmi - 1] *= absxk;
        }
      }
    }
  }
}

static void simulink_experiment_de_mldivide(const real_T A[32], const real_T B
  [32], real_T Y[16])
{
  __m128d tmp;
  __m128d tmp_0;
  real_T A_0[32];
  real_T B_0[32];
  real_T Q[32];
  real_T tau[4];
  real_T b_0;
  real_T scale;
  real_T tauj;
  real_T wj;
  int32_T jpvt[4];
  int32_T b;
  int32_T b_k;
  int32_T c_i;
  int32_T minmn;
  int32_T pj;
  int32_T rankA;
  int32_T scalarLB;
  int32_T vectorUB;
  boolean_T exitg1;
  memcpy(&A_0[0], &A[0], sizeof(real_T) << 5U);
  memcpy(&B_0[0], &B[0], sizeof(real_T) << 5U);
  simulink_experiment_debu_xgeqp3(A_0, Q, tau, jpvt);
  rankA = 0;
  tauj = Q[0];
  wj = fabs(tauj);
  scale = 1.7763568394002505E-14 * wj;
  exitg1 = false;
  while ((!exitg1) && (rankA < 4)) {
    tauj = Q[(rankA << 3) + rankA];
    wj = fabs(tauj);
    if (!(wj <= scale)) {
      rankA++;
    } else {
      exitg1 = true;
    }
  }

  for (pj = 0; pj < 16; pj++) {
    Y[pj] = 0.0;
  }

  memcpy(&A_0[0], &Q[0], sizeof(real_T) << 5U);
  for (pj = 0; pj < 4; pj++) {
    minmn = pj;
    tauj = tau[minmn];
    if (tauj != 0.0) {
      for (b_k = 0; b_k < 4; b_k++) {
        wj = B_0[(b_k << 3) + minmn];
        b = minmn;
        for (c_i = b + 2; c_i < 9; c_i++) {
          scale = A_0[((minmn << 3) + c_i) - 1];
          b_0 = B_0[((b_k << 3) + c_i) - 1];
          scale *= b_0;
          wj += scale;
        }

        wj *= tauj;
        if (wj != 0.0) {
          B_0[minmn + (b_k << 3)] -= wj;
          c_i = minmn;
          scalarLB = ((((7 - c_i) / 2) << 1) + c_i) + 2;
          vectorUB = scalarLB - 2;
          for (b = c_i + 2; b <= vectorUB; b += 2) {
            tmp = _mm_loadu_pd(&A_0[((minmn << 3) + b) - 1]);
            tmp = _mm_mul_pd(tmp, _mm_set1_pd(wj));
            tmp_0 = _mm_loadu_pd(&B_0[((b_k << 3) + b) - 1]);
            tmp = _mm_sub_pd(tmp_0, tmp);
            _mm_storeu_pd(&B_0[(b + (b_k << 3)) - 1], tmp);
          }

          for (b = scalarLB; b < 9; b++) {
            B_0[(b + (b_k << 3)) - 1] -= A_0[((minmn << 3) + b) - 1] * wj;
          }
        }
      }
    }
  }

  for (b_k = 0; b_k < 4; b_k++) {
    for (b = 0; b < rankA; b++) {
      c_i = b;
      Y[(jpvt[c_i] + (b_k << 2)) - 1] = B_0[(b_k << 3) + c_i];
    }

    for (minmn = rankA; minmn >= 1; minmn--) {
      pj = jpvt[minmn - 1] - 1;
      tauj = Y[(b_k << 2) + pj];
      wj = Q[(((minmn - 1) << 3) + minmn) - 1];
      scale = tauj / wj;
      Y[pj + (b_k << 2)] = scale;
      b = minmn - 2;
      for (c_i = 0; c_i <= b; c_i++) {
        Y[(jpvt[c_i] + (b_k << 2)) - 1] -= Q[((minmn - 1) << 3) + c_i] * Y[(b_k <<
          2) + pj];
      }
    }
  }
}

static void studentControllerInterface_step(studentControllerInterface_si_T *obj,
  real_T t, real_T p_ball, real_T theta, real_T *V_servo, real_T x_p[4], real_T
  P_p[16])
{
  __m128d tmp;
  __m128d tmp_0;
  studentControllerInterface_si_T *obj_0;
  real_T A[64];
  real_T B[64];
  real_T Z[64];
  real_T M[32];
  real_T N[32];
  real_T N_0[32];
  real_T A_lin[16];
  real_T G[16];
  real_T W11[16];
  real_T W21[16];
  real_T W22[16];
  real_T K[8];
  real_T a[8];
  real_T b_0[8];
  real_T S[4];
  real_T b[4];
  real_T y[2];
  real_T amp;
  real_T b_x;
  real_T br;
  real_T dt;
  real_T omega_min;
  real_T p_vel;
  real_T theta_vel;
  real_T x3;
  real_T xg_idx_0;
  real_T xg_idx_1;
  real_T xg_idx_3;
  int32_T A_lin_0;
  int32_T ONE;
  int32_T ib0;
  int32_T ipk;
  int32_T ix;
  int32_T jj;
  int32_T jm1;
  int32_T jp1j;
  int32_T jpiv;
  int32_T jpiv_offset;
  int32_T jprow;
  int32_T mmj;
  int8_T ipiv[8];
  int8_T p[8];
  static const int8_T tmp_1[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1 };

  /*         %% Main Controller Interface */
  /*  State Estimation */
  obj_0 = obj;

  /*         %% State Estimation: Generic Time Stepping */
  dt = t - obj_0->t_prev;
  p_vel = (p_ball - obj_0->p_prev) / dt;
  theta_vel = (theta - obj_0->theta_prev) / dt;
  if (dt == 0.0) {
    p_vel = 0.0;
    theta_vel = 0.0;
  }

  xg_idx_0 = p_ball;
  xg_idx_1 = p_vel;
  xg_idx_3 = theta_vel;
  obj_0 = obj;

  /*         %% State Estimation: Extended Kalman Filter */
  /*  Get some data */
  y[0] = p_ball;
  y[1] = theta;
  dt = t - obj_0->t_prev;
  x_p[2] = obj_0->x_hat[2];
  x_p[3] = obj_0->x_hat[3];
  if (dt > 0.0) {
    /*  Calculate kalman states */
    x3 = x_p[2];
    p_vel = x_p[3];
    theta_vel = rt_powd_snf(p_vel, 4.0);
    br = rt_powd_snf(p_vel, 3.0);
    p_vel = x3;
    p_vel = cos(p_vel);
    p_vel *= p_vel;
    amp = x3;
    amp = cos(amp);
    b_x = x3;
    b_x = sin(b_x);
    x3 = cos(x3);
    A_lin[1] = 0.0;
    A_lin[5] = 0.0;
    A_lin[9] = 0.0051 * amp * b_x * theta_vel + 0.4183 * x3;
    A_lin[13] = -0.0102 * br * p_vel;
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
    b[0] = obj_0->x_hat[0];
    b[1] = obj_0->x_hat[1];
    b[2] = obj_0->x_hat[2];
    b[3] = obj_0->x_hat[3];
    for (ix = 0; ix <= 2; ix += 2) {
      tmp = _mm_loadu_pd(&A_lin[ix]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(b[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&A_lin[ix + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(b[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&A_lin[ix + 8]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(b[2]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&A_lin[ix + 12]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(b[3]));
      tmp = _mm_add_pd(tmp_0, tmp);
      _mm_storeu_pd(&x_p[ix], tmp);
    }

    b[0] = obj_0->B[0];
    b[1] = obj_0->B[1];
    b[2] = obj_0->B[2];
    b[3] = obj_0->B[3];
    omega_min = obj_0->V_servo;
    theta_vel = x_p[0];
    p_vel = b[0];
    p_vel *= omega_min;
    theta_vel += p_vel;
    theta_vel *= dt;
    x_p[0] = theta_vel;
    theta_vel = x_p[1];
    p_vel = b[1];
    p_vel *= omega_min;
    theta_vel += p_vel;
    theta_vel *= dt;
    x_p[1] = theta_vel;
    theta_vel = x_p[2];
    p_vel = b[2];
    p_vel *= omega_min;
    theta_vel += p_vel;
    theta_vel *= dt;
    x_p[2] = theta_vel;
    theta_vel = x_p[3];
    p_vel = b[3];
    p_vel *= omega_min;
    theta_vel += p_vel;
    theta_vel *= dt;
    x_p[3] = theta_vel;
    theta_vel = x_p[0];
    theta_vel += obj_0->x_hat[0];
    x_p[0] = theta_vel;
    theta_vel = x_p[1];
    theta_vel += obj_0->x_hat[1];
    x_p[1] = theta_vel;
    theta_vel = x_p[2];
    theta_vel += obj_0->x_hat[2];
    x_p[2] = theta_vel;
    theta_vel = x_p[3];
    theta_vel += obj_0->x_hat[3];
    x_p[3] = theta_vel;
    for (ix = 0; ix < 16; ix++) {
      W11[ix] = obj_0->P[ix];
    }

    for (ix = 0; ix < 4; ix++) {
      for (A_lin_0 = 0; A_lin_0 <= 2; A_lin_0 += 2) {
        _mm_storeu_pd(&P_p[A_lin_0 + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A_lin[A_lin_0]);
        tmp = _mm_mul_pd(_mm_set1_pd(W11[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(W11[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(W11[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(W11[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&P_p[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&P_p[A_lin_0 + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix < 16; ix++) {
      G[ix] = obj_0->P[ix];
    }

    for (ix = 0; ix < 4; ix++) {
      for (A_lin_0 = 0; A_lin_0 < 4; A_lin_0++) {
        W11[ix + (A_lin_0 << 2)] = 0.0;
        theta_vel = W11[(A_lin_0 << 2) + ix];
        theta_vel += G[ix] * A_lin[A_lin_0];
        W11[ix + (A_lin_0 << 2)] = theta_vel;
        theta_vel = W11[(A_lin_0 << 2) + ix];
        theta_vel += G[ix + 4] * A_lin[A_lin_0 + 4];
        W11[ix + (A_lin_0 << 2)] = theta_vel;
        theta_vel = W11[(A_lin_0 << 2) + ix];
        theta_vel += G[ix + 8] * A_lin[A_lin_0 + 8];
        W11[ix + (A_lin_0 << 2)] = theta_vel;
        theta_vel = W11[(A_lin_0 << 2) + ix];
        theta_vel += G[ix + 12] * A_lin[A_lin_0 + 12];
        W11[ix + (A_lin_0 << 2)] = theta_vel;
      }
    }

    for (ix = 0; ix < 16; ix++) {
      theta_vel = P_p[ix];
      A_lin_0 = tmp_1[ix];
      theta_vel = (theta_vel + W11[ix]) + (real_T)A_lin_0;
      theta_vel *= dt;
      P_p[ix] = theta_vel;
    }

    for (ix = 0; ix < 16; ix++) {
      theta_vel = P_p[ix];
      theta_vel += obj_0->P[ix];
      P_p[ix] = theta_vel;
    }

    for (ix = 0; ix < 8; ix++) {
      a[ix] = obj_0->C[ix];
    }

    for (ix = 0; ix <= 0; ix += 2) {
      tmp = _mm_loadu_pd(&a[ix]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(x_p[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&a[ix + 2]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&a[ix + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[2]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&a[ix + 6]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(x_p[3]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&y[ix]);
      tmp = _mm_sub_pd(tmp_0, tmp);
      _mm_storeu_pd(&y[ix], tmp);
    }

    for (ix = 0; ix < 8; ix++) {
      b_0[ix] = obj_0->C[ix];
    }

    for (ix = 0; ix < 8; ix++) {
      a[ix] = obj_0->C[ix];
    }

    for (ix = 0; ix < 2; ix++) {
      for (A_lin_0 = 0; A_lin_0 < 4; A_lin_0++) {
        K[ix + (A_lin_0 << 1)] = 0.0;
        K[ix + (A_lin_0 << 1)] += P_p[A_lin_0 << 2] * a[ix];
        K[ix + (A_lin_0 << 1)] += P_p[(A_lin_0 << 2) + 1] * a[ix + 2];
        K[ix + (A_lin_0 << 1)] += P_p[(A_lin_0 << 2) + 2] * a[ix + 4];
        K[ix + (A_lin_0 << 1)] += P_p[(A_lin_0 << 2) + 3] * a[ix + 6];
      }

      for (A_lin_0 = 0; A_lin_0 < 2; A_lin_0++) {
        S[ix + (A_lin_0 << 1)] = 0.0;
        dt = S[(A_lin_0 << 1) + ix];
        dt += K[ix] * b_0[A_lin_0];
        S[ix + (A_lin_0 << 1)] = dt;
        dt = S[(A_lin_0 << 1) + ix];
        dt += K[ix + 2] * b_0[A_lin_0 + 2];
        S[ix + (A_lin_0 << 1)] = dt;
        dt = S[(A_lin_0 << 1) + ix];
        dt += K[ix + 4] * b_0[A_lin_0 + 4];
        S[ix + (A_lin_0 << 1)] = dt;
        dt = S[(A_lin_0 << 1) + ix];
        dt += K[ix + 6] * b_0[A_lin_0 + 6];
        S[ix + (A_lin_0 << 1)] = dt;
      }
    }

    S[0]++;
    S[3]++;
    for (ix = 0; ix < 8; ix++) {
      b_0[ix] = obj_0->C[ix];
    }

    for (ix = 0; ix < 16; ix++) {
      G[ix] = obj_0->P[ix];
    }

    for (ix = 0; ix < 4; ix++) {
      for (A_lin_0 = 0; A_lin_0 < 2; A_lin_0++) {
        a[ix + (A_lin_0 << 2)] = 0.0;
        a[ix + (A_lin_0 << 2)] += G[ix] * b_0[A_lin_0];
        a[ix + (A_lin_0 << 2)] += G[ix + 4] * b_0[A_lin_0 + 2];
        a[ix + (A_lin_0 << 2)] += G[ix + 8] * b_0[A_lin_0 + 4];
        a[ix + (A_lin_0 << 2)] += G[ix + 12] * b_0[A_lin_0 + 6];
      }
    }

    simulink_experiment_debug_mrdiv(a, S, K);
    for (ix = 0; ix <= 2; ix += 2) {
      tmp = _mm_loadu_pd(&K[ix]);
      tmp = _mm_mul_pd(tmp, _mm_set1_pd(y[0]));
      tmp = _mm_add_pd(tmp, _mm_set1_pd(0.0));
      tmp_0 = _mm_loadu_pd(&K[ix + 4]);
      tmp_0 = _mm_mul_pd(tmp_0, _mm_set1_pd(y[1]));
      tmp = _mm_add_pd(tmp_0, tmp);
      tmp_0 = _mm_loadu_pd(&x_p[ix]);
      tmp = _mm_add_pd(tmp_0, tmp);
      _mm_storeu_pd(&x_p[ix], tmp);
    }

    obj_0->x_hat[0] = x_p[0];
    obj_0->x_hat[1] = x_p[1];
    obj_0->x_hat[2] = x_p[2];
    obj_0->x_hat[3] = x_p[3];
    for (ix = 0; ix < 8; ix++) {
      b_0[ix] = obj_0->C[ix];
    }

    for (ix = 0; ix < 4; ix++) {
      for (A_lin_0 = 0; A_lin_0 <= 2; A_lin_0 += 2) {
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&K[A_lin_0]);
        tmp = _mm_mul_pd(_mm_set1_pd(b_0[ix << 1]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&K[A_lin_0 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(b_0[(ix << 1) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix < 16; ix++) {
      A_lin[ix] = 0.0;
    }

    A_lin[0] = 1.0;
    A_lin[5] = 1.0;
    A_lin[10] = 1.0;
    A_lin[15] = 1.0;
    for (ix = 0; ix <= 14; ix += 2) {
      tmp = _mm_loadu_pd(&A_lin[ix]);
      tmp_0 = _mm_loadu_pd(&W11[ix]);
      tmp = _mm_sub_pd(tmp, tmp_0);
      _mm_storeu_pd(&A_lin[ix], tmp);
    }

    for (ix = 0; ix < 4; ix++) {
      for (A_lin_0 = 0; A_lin_0 <= 2; A_lin_0 += 2) {
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], _mm_set1_pd(0.0));
        tmp = _mm_loadu_pd(&A_lin[A_lin_0]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[ix << 2]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 4]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(ix << 2) + 1]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 8]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(ix << 2) + 2]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
        tmp = _mm_loadu_pd(&A_lin[A_lin_0 + 12]);
        tmp = _mm_mul_pd(_mm_set1_pd(P_p[(ix << 2) + 3]), tmp);
        tmp_0 = _mm_loadu_pd(&W11[(ix << 2) + A_lin_0]);
        tmp = _mm_add_pd(tmp, tmp_0);
        _mm_storeu_pd(&W11[A_lin_0 + (ix << 2)], tmp);
      }
    }

    for (ix = 0; ix < 16; ix++) {
      obj_0->P[ix] = W11[ix];
    }
  }

  /*  Feedback Controller */
  x_p[0] = obj->x_hat[0];
  x_p[1] = obj->x_hat[1];
  x_p[2] = obj->x_hat[2];
  x_p[3] = obj->x_hat[3];
  for (ix = 0; ix < 16; ix++) {
    P_p[ix] = obj->P[ix];
  }

  /*  V_servo = stepImplP(obj, t, xk); */
  obj_0 = obj;

  /*         %% Feedback Controller: LQR */
  /*  fetch the previous values */
  if (t < 5.0) {
    omega_min = 0.0;
    dt = 0.0;
  } else if (t < 61.85) {
    br = t - 5.0;
    x3 = br / 56.85;
    if (x3 < 0.5) {
      amp = x3 / 0.5 * 0.090000000000000011 + 0.05;
      x3 = 0.11423973285781065 * br;
      x3 = sin(x3);
      x3 = 0.83775804095727813 * br - 0.2094395102393195 * x3 /
        0.11423973285781065;
      x3 = sin(x3);
      b_x = 0.11423973285781065 * br;
      b_x = sin(b_x);
      b_x = 0.83775804095727813 * br - 0.2094395102393195 * b_x /
        0.11423973285781065;
      b_x = cos(b_x);
      dt = 0.11423973285781065 * br;
      dt = cos(dt);
      dt = (0.83775804095727813 - 0.2094395102393195 * dt) * (amp * b_x) +
        0.00316622691292876 * x3;
    } else {
      amp = 0.14;
      x3 = 0.11423973285781065 * br;
      x3 = sin(x3);
      x3 = 0.83775804095727813 * br - 0.2094395102393195 * x3 /
        0.11423973285781065;
      x3 = cos(x3);
      b_x = 0.11423973285781065 * br;
      b_x = cos(b_x);
      dt = (0.83775804095727813 - 0.2094395102393195 * b_x) * (0.14 * x3);
    }

    x3 = 0.11423973285781065 * br;
    x3 = sin(x3);
    x3 = 0.83775804095727813 * br - 0.2094395102393195 * x3 /
      0.11423973285781065;
    x3 = sin(x3);
    omega_min = amp * x3;
  } else if (t < 65.0) {
    omega_min = 0.0;
    dt = 0.0;
  } else if (t < 85.0) {
    p_vel = t - 65.0;
    theta_vel = p_vel / 20.0;
    if (theta_vel < 0.5) {
      theta_vel = 0.05;
    } else {
      theta_vel = 0.1;
    }

    x3 = 0.62831853071795862 * p_vel;
    x3 = sin(x3);
    if (x3 < 0.0) {
      x3 = -1.0;
    } else {
      x3 = (x3 > 0.0);
    }

    omega_min = theta_vel * x3;
    dt = 0.0;
  } else {
    omega_min = 0.0;
    dt = 0.0;
  }

  x3 = theta;
  p_vel = xg_idx_3;
  theta_vel = rt_powd_snf(p_vel, 4.0);
  br = rt_powd_snf(p_vel, 3.0);
  p_vel = x3;
  p_vel = cos(p_vel);
  p_vel *= p_vel;
  amp = x3;
  amp = cos(amp);
  b_x = x3;
  b_x = sin(b_x);
  x3 = cos(x3);
  A_lin[1] = 0.0;
  A_lin[5] = 0.0;
  A_lin[9] = 0.0051 * amp * b_x * theta_vel + 0.4183 * x3;
  A_lin[13] = -0.0102 * br * p_vel;
  omega_min += 0.008;
  xg_idx_0 -= omega_min;
  xg_idx_1 -= dt;
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
  b[0] = obj_0->B[0];
  b[1] = obj_0->B[1];
  b[2] = obj_0->B[2];
  b[3] = obj_0->B[3];
  for (ix = 0; ix < 16; ix++) {
    W11[ix] = obj_0->Q[ix];
  }

  p_vel = obj_0->R;

  /* LQR_VIA_MATRIX_SIGN_FUNCTION Computes the LQR gain using matrix sign function method */
  /*  */
  /*    K = LQR_VIA_MATRIX_SIGN_FUNCTION(A, B, Q, R) returns the optimal gain matrix K */
  /*    for the continuous-time LQR problem: */
  /*  */
  /*        minimize J = &#x222B; (x'Qx + u'Ru) dt */
  /*        subject to dx/dt = Ax + Bu */
  /*  */
  /*    Inputs: */
  /*        A - System dynamics matrix (n x n) */
  /*        B - Input matrix (n x m) */
  /*        Q - State cost matrix (n x n), symmetric positive semi-definite */
  /*        R - Input cost matrix (m x m), symmetric positive definite */
  /*  */
  /*    Output: */
  /*        K - Optimal state feedback gain matrix (m x n) */
  /*  Invert R (assumes R is positive definite) */
  theta_vel = 1.0 / p_vel;
  S[0] = b[0] * theta_vel;
  S[1] = b[1] * theta_vel;
  S[2] = b[2] * theta_vel;
  S[3] = b[3] * theta_vel;

  /*  Construct the Hamiltonian matrix Z */
  for (ix = 0; ix < 4; ix++) {
    p_vel = b[ix];
    dt = S[0] * p_vel;
    Z[ix << 3] = A_lin[ix << 2];
    Z[(ix + 4) << 3] = -dt;
    Z[(ix << 3) + 4] = -W11[ix << 2];
    Z[((ix + 4) << 3) + 4] = -A_lin[ix];
    dt = S[1] * p_vel;
    Z[(ix << 3) + 1] = A_lin[(ix << 2) + 1];
    Z[((ix + 4) << 3) + 1] = -dt;
    Z[(ix << 3) + 5] = -W11[(ix << 2) + 1];
    Z[((ix + 4) << 3) + 5] = -A_lin[ix + 4];
    dt = S[2] * p_vel;
    Z[(ix << 3) + 2] = A_lin[(ix << 2) + 2];
    Z[((ix + 4) << 3) + 2] = -dt;
    Z[(ix << 3) + 6] = -W11[(ix << 2) + 2];
    Z[((ix + 4) << 3) + 6] = -A_lin[ix + 8];
    dt = S[3] * p_vel;
    Z[(ix << 3) + 3] = A_lin[(ix << 2) + 3];
    Z[((ix + 4) << 3) + 3] = -dt;
    Z[(ix << 3) + 7] = -W11[(ix << 2) + 3];
    Z[((ix + 4) << 3) + 7] = -A_lin[ix + 12];
  }

  /*  Initialize matrix W for iteration */
  /*  W = Z; */
  /*  Newton iteration to compute the matrix sign function */
  for (A_lin_0 = 0; A_lin_0 < 1000; A_lin_0++) {
    for (ix = 0; ix < 64; ix++) {
      A[ix] = Z[ix];
      B[ix] = 0.0;
    }

    for (ix = 0; ix < 8; ix++) {
      ipiv[ix] = (int8_T)(ix + 1);
    }

    for (ib0 = 0; ib0 < 7; ib0++) {
      ipk = ib0 + 1;
      jm1 = ipk - 1;
      mmj = 8 - ipk;
      jprow = jm1 * 9;
      jj = jprow + 1;
      jp1j = jj + 1;
      jprow = mmj + 1;
      ONE = 1;
      ix = jj - 1;
      x3 = A[jj - 1];
      dt = fabs(x3);
      omega_min = dt;
      for (jpiv = 2; jpiv <= jprow; jpiv++) {
        ix++;
        x3 = A[ix];
        dt = fabs(x3);
        if (dt > omega_min) {
          ONE = jpiv;
          omega_min = dt;
        }
      }

      jpiv_offset = ONE - 1;
      jpiv = (jj + jpiv_offset) - 1;
      if (A[jpiv] != 0.0) {
        if (jpiv_offset != 0) {
          jprow = ipk + jpiv_offset;
          ipiv[ipk - 1] = (int8_T)jprow;
          ONE = jm1;
          jprow = ONE + jpiv_offset;
          for (jpiv = 0; jpiv < 8; jpiv++) {
            x3 = A[ONE];
            A[ONE] = A[jprow];
            A[jprow] = x3;
            ONE += 8;
            jprow += 8;
          }
        }

        jprow = mmj - 1;
        jpiv_offset = jp1j + jprow;
        for (ix = jp1j; ix <= jpiv_offset; ix++) {
          x3 = A[ix - 1];
          dt = A[jj - 1];
          dt = x3 / dt;
          A[ix - 1] = dt;
        }
      }

      jprow = 7 - ipk;
      ONE = jj + 7;
      jj += 9;
      jm1 = jj - 1;
      for (ipk = 0; ipk <= jprow; ipk++) {
        x3 = A[ONE];
        if (x3 != 0.0) {
          x3 = -x3;
          ix = jp1j - 1;
          jpiv_offset = jm1 + 1;
          jj = mmj + jm1;
          for (jpiv = jpiv_offset; jpiv <= jj; jpiv++) {
            A[jpiv - 1] += A[ix] * x3;
            ix++;
          }
        }

        ONE += 8;
        jm1 += 8;
      }
    }

    for (ix = 0; ix < 8; ix++) {
      p[ix] = (int8_T)(ix + 1);
    }

    for (ib0 = 0; ib0 < 7; ib0++) {
      dt = (real_T)ib0 + 1.0;
      ipk = ipiv[(int32_T)dt - 1] - 1;
      if (ipk + 1 > (int32_T)dt) {
        jpiv_offset = p[ipk];
        p[ipk] = p[(int32_T)dt - 1];
        p[(int32_T)dt - 1] = (int8_T)jpiv_offset;
      }
    }

    for (ib0 = 0; ib0 < 8; ib0++) {
      jpiv = ib0;
      jprow = p[jpiv] - 1;
      B[jpiv + (jprow << 3)] = 1.0;
      for (ipk = jpiv + 1; ipk < 9; ipk++) {
        if (B[((jprow << 3) + ipk) - 1] != 0.0) {
          jpiv_offset = ipk + 1;
          for (ix = jpiv_offset; ix < 9; ix++) {
            B[(ix + (jprow << 3)) - 1] -= A[(((ipk - 1) << 3) + ix) - 1] * B
              [((jprow << 3) + ipk) - 1];
          }
        }
      }
    }

    for (ib0 = 0; ib0 < 8; ib0++) {
      ipk = ib0;
      jm1 = ipk << 3;
      for (jpiv = 7; jpiv >= 0; jpiv--) {
        mmj = jpiv << 3;
        if (B[jpiv + jm1] != 0.0) {
          B[jpiv + jm1] /= A[jpiv + mmj];
          jpiv_offset = jpiv - 1;
          for (ipk = 0; ipk <= jpiv_offset; ipk++) {
            ix = ipk;
            B[ix + jm1] -= B[jpiv + jm1] * A[ix + mmj];
          }
        }
      }
    }

    for (ix = 0; ix <= 62; ix += 2) {
      tmp = _mm_loadu_pd(&Z[ix]);
      tmp_0 = _mm_loadu_pd(&B[ix]);
      tmp_0 = _mm_sub_pd(tmp, tmp_0);
      tmp_0 = _mm_mul_pd(_mm_set1_pd(0.5), tmp_0);
      tmp = _mm_sub_pd(tmp, tmp_0);
      _mm_storeu_pd(&Z[ix], tmp);
      _mm_storeu_pd(&B[ix], tmp_0);
    }
  }

  /*  Determine the size of the system */
  /*  Partition W into 4 submatrices */
  for (ix = 0; ix < 4; ix++) {
    W11[ix << 2] = Z[ix << 3];
    G[ix << 2] = Z[(ix + 4) << 3];
    W21[ix << 2] = Z[(ix << 3) + 4];
    W22[ix << 2] = Z[((ix + 4) << 3) + 4];
    W11[(ix << 2) + 1] = Z[(ix << 3) + 1];
    G[(ix << 2) + 1] = Z[((ix + 4) << 3) + 1];
    W21[(ix << 2) + 1] = Z[(ix << 3) + 5];
    W22[(ix << 2) + 1] = Z[((ix + 4) << 3) + 5];
    W11[(ix << 2) + 2] = Z[(ix << 3) + 2];
    G[(ix << 2) + 2] = Z[((ix + 4) << 3) + 2];
    W21[(ix << 2) + 2] = Z[(ix << 3) + 6];
    W22[(ix << 2) + 2] = Z[((ix + 4) << 3) + 6];
    W11[(ix << 2) + 3] = Z[(ix << 3) + 3];
    G[(ix << 2) + 3] = Z[((ix + 4) << 3) + 3];
    W21[(ix << 2) + 3] = Z[(ix << 3) + 7];
    W22[(ix << 2) + 3] = Z[((ix + 4) << 3) + 7];
  }

  /*  Solve for the unique positive semidefinite solution P to the Riccati equation */
  for (ix = 0; ix < 16; ix++) {
    A_lin[ix] = 0.0;
  }

  for (ib0 = 0; ib0 < 4; ib0++) {
    jpiv = ib0;
    A_lin[jpiv + (jpiv << 2)] = 1.0;
    M[ib0 << 3] = G[ib0 << 2];
    M[(ib0 << 3) + 1] = G[(ib0 << 2) + 1];
    M[(ib0 << 3) + 2] = G[(ib0 << 2) + 2];
    M[(ib0 << 3) + 3] = G[(ib0 << 2) + 3];
    M[(ib0 << 3) + 4] = W22[ib0 << 2] + A_lin[ib0 << 2];
    M[(ib0 << 3) + 5] = W22[(ib0 << 2) + 1] + A_lin[(ib0 << 2) + 1];
    M[(ib0 << 3) + 6] = W22[(ib0 << 2) + 2] + A_lin[(ib0 << 2) + 2];
    M[(ib0 << 3) + 7] = W22[(ib0 << 2) + 3] + A_lin[(ib0 << 2) + 3];
  }

  for (ix = 0; ix < 16; ix++) {
    A_lin[ix] = 0.0;
  }

  A_lin[0] = 1.0;
  A_lin[5] = 1.0;
  A_lin[10] = 1.0;
  A_lin[15] = 1.0;
  for (ix = 0; ix < 4; ix++) {
    N[ix << 3] = W11[ix << 2] + A_lin[ix << 2];
    N[(ix << 3) + 4] = W21[ix << 2];
    N[(ix << 3) + 1] = W11[(ix << 2) + 1] + A_lin[(ix << 2) + 1];
    N[(ix << 3) + 5] = W21[(ix << 2) + 1];
    N[(ix << 3) + 2] = W11[(ix << 2) + 2] + A_lin[(ix << 2) + 2];
    N[(ix << 3) + 6] = W21[(ix << 2) + 2];
    N[(ix << 3) + 3] = W11[(ix << 2) + 3] + A_lin[(ix << 2) + 3];
    N[(ix << 3) + 7] = W21[(ix << 2) + 3];
  }

  for (ix = 0; ix <= 30; ix += 2) {
    tmp = _mm_loadu_pd(&N[ix]);
    tmp = _mm_mul_pd(tmp, _mm_set1_pd(-1.0));
    _mm_storeu_pd(&N_0[ix], tmp);
  }

  simulink_experiment_de_mldivide(M, N_0, A_lin);

  /*  Compute the optimal LQR gain */
  p_vel = b[0];
  p_vel *= theta_vel;
  b[0] = p_vel;
  p_vel = b[1];
  p_vel *= theta_vel;
  b[1] = p_vel;
  p_vel = b[2];
  p_vel *= theta_vel;
  b[2] = p_vel;
  p_vel = b[3];
  p_vel *= theta_vel;
  b[3] = p_vel;
  for (ix = 0; ix < 4; ix++) {
    dt = A_lin[ix << 2] * b[0];
    dt += A_lin[(ix << 2) + 1] * b[1];
    dt += A_lin[(ix << 2) + 2] * b[2];
    dt += A_lin[(ix << 2) + 3] * b[3];
    S[ix] = dt;
  }

  /*  coder.extrinsic('lqr') */
  /*  K = [8.6603 10.0662 2.4465 0.0586]; */
  /*  K = lqr(A_lin, obj.B, obj.Q, obj.R); */
  /*  disp(K) */
  S[0] = -S[0];
  S[1] = -S[1];
  S[2] = -S[2];
  S[3] = -S[3];
  *V_servo = S[0] * xg_idx_0;
  *V_servo += S[1] * xg_idx_1;
  *V_servo += S[2] * theta;
  *V_servo += S[3] * xg_idx_3;

  /*  V_servo = stepImplLQG(obj, t, xg); */
  p_vel = *V_servo - obj->V_servo;
  dt = fabs(p_vel);
  if (dt > 0.04) {
    if (rtIsNaN(p_vel)) {
      p_vel = (rtNaN);
    } else if (p_vel < 0.0) {
      p_vel = -1.0;
    } else {
      p_vel = (p_vel > 0.0);
    }

    *V_servo = 0.04 * p_vel + obj->V_servo;
  } else {
    dt = fabs(p_vel);
    if (dt > 0.0025) {
      *V_servo = obj->V_servo;
    }
  }

  if (*V_servo > 5.0) {
    *V_servo = 5.0;
  } else if (*V_servo < -5.0) {
    *V_servo = -5.0;
  }

  obj->V_servo = *V_servo;
  obj->t_prev = t;
  obj->p_prev = p_ball;
  obj->theta_prev = theta;
}

/* Model output function for TID0 */
void simulink_experiment_debug_type1_output0(void) /* Sample time: [0.0s, 0.0s] */
{
  studentControllerInterface_si_T *obj;
  real_T tmp_0[16];
  real_T tmp[4];
  real_T amplitude;
  real_T b_varargout_1;
  real_T u1;
  real_T u2;

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
  amplitude = simulink_experiment_debug_typ_B.Clock;
  u1 = simulink_experiment_debug_typ_B.BB01SensorGainmV;
  u2 = simulink_experiment_debug_typ_B.Bias;
  obj = &simulink_experiment_debug_ty_DW.obj;
  studentControllerInterface_step(obj, amplitude, u1, u2, &b_varargout_1, tmp,
    tmp_0);

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o1 = b_varargout_1;

  /* MATLABSystem: '<Root>/MATLAB System' */
  simulink_experiment_debug_typ_B.MATLABSystem_o2[0] = tmp[0];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[1] = tmp[1];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[2] = tmp[2];
  simulink_experiment_debug_typ_B.MATLABSystem_o2[3] = tmp[3];

  /* MATLABSystem: '<Root>/MATLAB System' */
  memcpy(&simulink_experiment_debug_typ_B.MATLABSystem_o3[0], &tmp_0[0], sizeof
         (real_T) << 4U);

  /* Saturate: '<Root>/+//-10V' */
  amplitude = simulink_experiment_debug_typ_B.MATLABSystem_o1;
  u1 = simulink_experiment_debug_typ_P.u0V_LowerSat;
  u2 = simulink_experiment_debug_typ_P.u0V_UpperSat;
  if (amplitude > u2) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u2;
  } else if (amplitude < u1) {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = u1;
  } else {
    /* Saturate: '<Root>/+//-10V' */
    simulink_experiment_debug_typ_B.u0V = amplitude;
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
    amplitude = (simulink_experiment_debug_typ_B.Clock - 5.0) / 56.85;
    if (amplitude < 0.5) {
      amplitude = amplitude / 0.5 * 0.090000000000000011 + 0.05;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * amplitude *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195) + sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.00316622691292876;
      u1 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
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
        5.0) * 6.0 / 1895.0 + 0.05) * (u1 * u1)) + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 19.739208802178716) *
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.0 / 1895.0 + 0.05) /
        825.0;
    } else {
      amplitude = 0.14;
      simulink_experiment_debug_typ_B.v_ref = cos
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 -
         sin((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065)
         * 0.2094395102393195 / 0.11423973285781065) * 0.14 *
        (0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock - 5.0)
          * 0.11423973285781065) * 0.2094395102393195);
      u1 = 0.83775804095727813 - cos((simulink_experiment_debug_typ_B.Clock -
        5.0) * 6.2831853071795862 / 55.0) * 3.1415926535897931 / 15.0;
      simulink_experiment_debug_typ_B.a_ref = sin(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * 7.0 * (u1 * u1) / 50.0 + cos(sin
        ((simulink_experiment_debug_typ_B.Clock - 5.0) * 6.2831853071795862 /
         55.0) * 11.0 / 6.0 - (simulink_experiment_debug_typ_B.Clock - 5.0) *
        12.566370614359172 / 15.0) * (sin((simulink_experiment_debug_typ_B.Clock
        - 5.0) * 6.2831853071795862 / 55.0) * 69.0872308076255) / 20625.0;
    }

    simulink_experiment_debug_typ_B.p_ref = sin
      ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.83775804095727813 - sin
       ((simulink_experiment_debug_typ_B.Clock - 5.0) * 0.11423973285781065) *
       0.2094395102393195 / 0.11423973285781065) * amplitude;
  } else if (simulink_experiment_debug_typ_B.Clock < 65.0) {
    simulink_experiment_debug_typ_B.p_ref = 0.0;
    simulink_experiment_debug_typ_B.v_ref = 0.0;
    simulink_experiment_debug_typ_B.a_ref = 0.0;
  } else if (simulink_experiment_debug_typ_B.Clock < 85.0) {
    if ((simulink_experiment_debug_typ_B.Clock - 65.0) / 20.0 < 0.5) {
      amplitude = 0.05;
    } else {
      amplitude = 0.1;
    }

    u1 = sin((simulink_experiment_debug_typ_B.Clock - 65.0) *
             0.62831853071795862);
    if (rtIsNaN(u1)) {
      u1 = (rtNaN);
    } else if (u1 < 0.0) {
      u1 = -1.0;
    } else {
      u1 = (u1 > 0.0);
    }

    simulink_experiment_debug_typ_B.p_ref = amplitude * u1;
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

    static const int8_T tmp_0[16] = { 125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0 };

    static const real_T tmp_1[16] = { 0.00467317845216159, 0.00447778693210613,
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

    for (i = 0; i < 16; i++) {
      b_obj->Q[i] = tmp_0[i];
    }

    b_obj->R = 1.0;
    b_obj->x_hat[0] = -0.05;
    b_obj->x_hat[1] = 0.0;
    b_obj->x_hat[2] = 0.0;
    b_obj->x_hat[3] = 0.0;
    for (i = 0; i < 16; i++) {
      b_obj->P[i] = tmp_1[i];
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
  simulink_experiment_debug_ty_M->Sizes.checksums[0] = (2095139769U);
  simulink_experiment_debug_ty_M->Sizes.checksums[1] = (2315151738U);
  simulink_experiment_debug_ty_M->Sizes.checksums[2] = (1886773115U);
  simulink_experiment_debug_ty_M->Sizes.checksums[3] = (3728185577U);

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
