/*
 * simulink_experiment_test_type1.c
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "simulink_experiment_test_type1".
 *
 * Model version              : 1.3
 * Simulink Coder version : 9.8 (R2022b) 13-May-2022
 * C source code generated on : Wed Apr  9 12:45:32 2025
 *
 * Target selection: quarc_win64.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "simulink_experiment_test_type1.h"
#include "simulink_experiment_test_type1_private.h"
#include <string.h>
#include "rt_nonfinite.h"
#include "simulink_experiment_test_type1_dt.h"

/* Block signals (default storage) */
B_simulink_experiment_test_ty_T simulink_experiment_test_type_B;

/* Block states (default storage) */
DW_simulink_experiment_test_t_T simulink_experiment_test_typ_DW;

/* Real-time model */
static RT_MODEL_simulink_experiment__T simulink_experiment_test_typ_M_;
RT_MODEL_simulink_experiment__T *const simulink_experiment_test_typ_M =
  &simulink_experiment_test_typ_M_;

/* Model output function */
void simulink_experiment_test_type1_output(void)
{
  /* S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_read_encoder
      (simulink_experiment_test_typ_DW.HILReadEncoderTimebase_Task, 1,
       &simulink_experiment_test_typ_DW.HILReadEncoderTimebase_Buffer);
    if (result < 0) {
      simulink_experiment_test_type_B.HILReadEncoderTimebase = 0;
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
    } else {
      simulink_experiment_test_type_B.HILReadEncoderTimebase =
        simulink_experiment_test_typ_DW.HILReadEncoderTimebase_Buffer;
    }
  }

  /* S-Function (hil_read_analog_block): '<S1>/HIL Read Analog' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Read Analog (hil_read_analog_block) */
  {
    t_error result = hil_read_analog
      (simulink_experiment_test_typ_DW.HILInitialize_Card,
       &simulink_experiment_test_type_P.HILReadAnalog_channels, 1,
       &simulink_experiment_test_typ_DW.HILReadAnalog_Buffer);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
    }

    simulink_experiment_test_type_B.HILReadAnalog =
      simulink_experiment_test_typ_DW.HILReadAnalog_Buffer;
  }

  /* Gain: '<S1>/BB01 Sensor  Gain (m//V)' */
  simulink_experiment_test_type_B.BB01SensorGainmV =
    simulink_experiment_test_type_P.BB01SensorGainmV_Gain *
    simulink_experiment_test_type_B.HILReadAnalog;

  /* Bias: '<S1>/Position bias' */
  simulink_experiment_test_type_B.Positionbias =
    simulink_experiment_test_type_B.BB01SensorGainmV +
    simulink_experiment_test_type_P.Positionbias_Bias;

  /* Gain: '<S1>/Encoder Calibration  (rad//count)' */
  simulink_experiment_test_type_B.EncoderCalibrationradcount =
    simulink_experiment_test_type_P.EncoderCalibrationradcount_Gain *
    simulink_experiment_test_type_B.HILReadEncoderTimebase;

  /* Bias: '<S1>/Bias' */
  simulink_experiment_test_type_B.Bias =
    simulink_experiment_test_type_B.EncoderCalibrationradcount +
    simulink_experiment_test_type_P.Bias_Bias;

  /* Gain: '<S1>/Motor  Gain (V//V)' incorporates:
   *  Constant: '<Root>/Constant'
   */
  simulink_experiment_test_type_B.MotorGainVV =
    simulink_experiment_test_type_P.MotorGainVV_Gain *
    simulink_experiment_test_type_P.Constant_Value;

  /* S-Function (hil_write_analog_block): '<S1>/HIL Write Analog' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Write Analog (hil_write_analog_block) */
  {
    t_error result;
    result = hil_write_analog(simulink_experiment_test_typ_DW.HILInitialize_Card,
      &simulink_experiment_test_type_P.HILWriteAnalog_channels, 1,
      &simulink_experiment_test_type_B.MotorGainVV);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
    }
  }

  /* Gain: '<Root>/m to cm' */
  simulink_experiment_test_type_B.mtocm =
    simulink_experiment_test_type_P.mtocm_Gain *
    simulink_experiment_test_type_B.Positionbias;

  /* Gain: '<S2>/Gain' */
  simulink_experiment_test_type_B.Gain =
    simulink_experiment_test_type_P.Gain_Gain *
    simulink_experiment_test_type_B.Bias;
}

/* Model update function */
void simulink_experiment_test_type1_update(void)
{
  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++simulink_experiment_test_typ_M->Timing.clockTick0)) {
    ++simulink_experiment_test_typ_M->Timing.clockTickH0;
  }

  simulink_experiment_test_typ_M->Timing.t[0] =
    simulink_experiment_test_typ_M->Timing.clockTick0 *
    simulink_experiment_test_typ_M->Timing.stepSize0 +
    simulink_experiment_test_typ_M->Timing.clockTickH0 *
    simulink_experiment_test_typ_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void simulink_experiment_test_type1_initialize(void)
{
  /* Start for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_int result;
    t_boolean is_switching;
    result = hil_open("q2_usb", "0",
                      &simulink_experiment_test_typ_DW.HILInitialize_Card);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
      return;
    }

    is_switching = false;
    result = hil_set_card_specific_options
      (simulink_experiment_test_typ_DW.HILInitialize_Card,
       "d0=digital;d1=digital;led=auto;update_rate=normal", 50);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
      return;
    }

    result = hil_watchdog_clear
      (simulink_experiment_test_typ_DW.HILInitialize_Card);
    if (result < 0 && result != -QERR_HIL_WATCHDOG_CLEAR) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
      return;
    }

    if ((simulink_experiment_test_type_P.HILInitialize_AIPStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_AIPEnter &&
            is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_AIMinimums[0] =
        (simulink_experiment_test_type_P.HILInitialize_AILow);
      simulink_experiment_test_typ_DW.HILInitialize_AIMinimums[1] =
        (simulink_experiment_test_type_P.HILInitialize_AILow);
      simulink_experiment_test_typ_DW.HILInitialize_AIMaximums[0] =
        simulink_experiment_test_type_P.HILInitialize_AIHigh;
      simulink_experiment_test_typ_DW.HILInitialize_AIMaximums[1] =
        simulink_experiment_test_type_P.HILInitialize_AIHigh;
      result = hil_set_analog_input_ranges
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_AIChannels, 2U,
         &simulink_experiment_test_typ_DW.HILInitialize_AIMinimums[0],
         &simulink_experiment_test_typ_DW.HILInitialize_AIMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if ((simulink_experiment_test_type_P.HILInitialize_AOPStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_AOPEnter &&
            is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_AOMinimums[0] =
        (simulink_experiment_test_type_P.HILInitialize_AOLow);
      simulink_experiment_test_typ_DW.HILInitialize_AOMinimums[1] =
        (simulink_experiment_test_type_P.HILInitialize_AOLow);
      simulink_experiment_test_typ_DW.HILInitialize_AOMaximums[0] =
        simulink_experiment_test_type_P.HILInitialize_AOHigh;
      simulink_experiment_test_typ_DW.HILInitialize_AOMaximums[1] =
        simulink_experiment_test_type_P.HILInitialize_AOHigh;
      result = hil_set_analog_output_ranges
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_AOChannels, 2U,
         &simulink_experiment_test_typ_DW.HILInitialize_AOMinimums[0],
         &simulink_experiment_test_typ_DW.HILInitialize_AOMaximums[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if ((simulink_experiment_test_type_P.HILInitialize_AOStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_AOEnter &&
            is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_test_type_P.HILInitialize_AOInitial;
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_test_type_P.HILInitialize_AOInitial;
      result = hil_write_analog
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_AOChannels, 2U,
         &simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if (simulink_experiment_test_type_P.HILInitialize_AOReset) {
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_test_type_P.HILInitialize_AOWatchdog;
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_test_type_P.HILInitialize_AOWatchdog;
      result = hil_watchdog_set_analog_expiration_state
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_AOChannels, 2U,
         &simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    result = hil_set_digital_directions
      (simulink_experiment_test_typ_DW.HILInitialize_Card, NULL, 0U,
       simulink_experiment_test_type_P.HILInitialize_DOChannels, 8U);
    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
      return;
    }

    if ((simulink_experiment_test_type_P.HILInitialize_DOStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_DOEnter &&
            is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_test_typ_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] =
            simulink_experiment_test_type_P.HILInitialize_DOInitial;
        }
      }

      result = hil_write_digital
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_DOChannels, 8U,
         (t_boolean *) &simulink_experiment_test_typ_DW.HILInitialize_DOBits[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if (simulink_experiment_test_type_P.HILInitialize_DOReset) {
      {
        int_T i1;
        int32_T *dw_DOStates =
          &simulink_experiment_test_typ_DW.HILInitialize_DOStates[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOStates[i1] =
            simulink_experiment_test_type_P.HILInitialize_DOWatchdog;
        }
      }

      result = hil_watchdog_set_digital_expiration_state
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_DOChannels, 8U, (const
          t_digital_state *)
         &simulink_experiment_test_typ_DW.HILInitialize_DOStates[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if ((simulink_experiment_test_type_P.HILInitialize_EIPStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_EIPEnter &&
            is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_QuadratureModes[0] =
        simulink_experiment_test_type_P.HILInitialize_EIQuadrature;
      simulink_experiment_test_typ_DW.HILInitialize_QuadratureModes[1] =
        simulink_experiment_test_type_P.HILInitialize_EIQuadrature;
      result = hil_set_encoder_quadrature_mode
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_EIChannels, 2U,
         (t_encoder_quadrature_mode *)
         &simulink_experiment_test_typ_DW.HILInitialize_QuadratureModes[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }

    if ((simulink_experiment_test_type_P.HILInitialize_EIStart && !is_switching)
        || (simulink_experiment_test_type_P.HILInitialize_EIEnter &&
            is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_InitialEICounts[0] =
        simulink_experiment_test_type_P.HILInitialize_EIInitial;
      simulink_experiment_test_typ_DW.HILInitialize_InitialEICounts[1] =
        simulink_experiment_test_type_P.HILInitialize_EIInitial;
      result = hil_set_encoder_counts
        (simulink_experiment_test_typ_DW.HILInitialize_Card,
         simulink_experiment_test_type_P.HILInitialize_EIChannels, 2U,
         &simulink_experiment_test_typ_DW.HILInitialize_InitialEICounts[0]);
      if (result < 0) {
        msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
          (_rt_error_message));
        rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        return;
      }
    }
  }

  /* Start for S-Function (hil_read_encoder_timebase_block): '<S1>/HIL Read Encoder Timebase' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Read Encoder Timebase (hil_read_encoder_timebase_block) */
  {
    t_error result;
    result = hil_task_create_encoder_reader
      (simulink_experiment_test_typ_DW.HILInitialize_Card,
       simulink_experiment_test_type_P.HILReadEncoderTimebase_SamplesI,
       &simulink_experiment_test_type_P.HILReadEncoderTimebase_Channels, 1,
       &simulink_experiment_test_typ_DW.HILReadEncoderTimebase_Task);
    if (result >= 0) {
      result = hil_task_set_buffer_overflow_mode
        (simulink_experiment_test_typ_DW.HILReadEncoderTimebase_Task,
         (t_buffer_overflow_mode)
         (simulink_experiment_test_type_P.HILReadEncoderTimebase_Overflow - 1));
    }

    if (result < 0) {
      msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
        (_rt_error_message));
      rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
    }
  }
}

/* Model terminate function */
void simulink_experiment_test_type1_terminate(void)
{
  /* Terminate for S-Function (hil_initialize_block): '<S1>/HIL Initialize' */

  /* S-Function Block: simulink_experiment_test_type1/Ball and Beam Hardware Interface/HIL Initialize (hil_initialize_block) */
  {
    t_boolean is_switching;
    t_int result;
    t_uint32 num_final_analog_outputs = 0;
    t_uint32 num_final_digital_outputs = 0;
    hil_task_stop_all(simulink_experiment_test_typ_DW.HILInitialize_Card);
    hil_monitor_stop_all(simulink_experiment_test_typ_DW.HILInitialize_Card);
    is_switching = false;
    if ((simulink_experiment_test_type_P.HILInitialize_AOTerminate &&
         !is_switching) || (simulink_experiment_test_type_P.HILInitialize_AOExit
         && is_switching)) {
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0] =
        simulink_experiment_test_type_P.HILInitialize_AOFinal;
      simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[1] =
        simulink_experiment_test_type_P.HILInitialize_AOFinal;
      num_final_analog_outputs = 2U;
    } else {
      num_final_analog_outputs = 0;
    }

    if ((simulink_experiment_test_type_P.HILInitialize_DOTerminate &&
         !is_switching) || (simulink_experiment_test_type_P.HILInitialize_DOExit
         && is_switching)) {
      {
        int_T i1;
        boolean_T *dw_DOBits =
          &simulink_experiment_test_typ_DW.HILInitialize_DOBits[0];
        for (i1=0; i1 < 8; i1++) {
          dw_DOBits[i1] = simulink_experiment_test_type_P.HILInitialize_DOFinal;
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
      result = hil_write(simulink_experiment_test_typ_DW.HILInitialize_Card
                         ,
                         simulink_experiment_test_type_P.HILInitialize_AOChannels,
                         num_final_analog_outputs
                         , NULL, 0
                         ,
                         simulink_experiment_test_type_P.HILInitialize_DOChannels,
                         num_final_digital_outputs
                         , NULL, 0
                         ,
                         &simulink_experiment_test_typ_DW.HILInitialize_AOVoltages
                         [0]
                         , NULL
                         , (t_boolean *)
                         &simulink_experiment_test_typ_DW.HILInitialize_DOBits[0]
                         , NULL
                         );
      if (result == -QERR_HIL_WRITE_NOT_SUPPORTED) {
        t_error local_result;
        result = 0;

        /* The hil_write operation is not supported by this card. Write final outputs for each channel type */
        if (num_final_analog_outputs > 0) {
          local_result = hil_write_analog
            (simulink_experiment_test_typ_DW.HILInitialize_Card,
             simulink_experiment_test_type_P.HILInitialize_AOChannels,
             num_final_analog_outputs,
             &simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (num_final_digital_outputs > 0) {
          local_result = hil_write_digital
            (simulink_experiment_test_typ_DW.HILInitialize_Card,
             simulink_experiment_test_type_P.HILInitialize_DOChannels,
             num_final_digital_outputs, (t_boolean *)
             &simulink_experiment_test_typ_DW.HILInitialize_DOBits[0]);
          if (local_result < 0) {
            result = local_result;
          }
        }

        if (result < 0) {
          msg_get_error_messageA(NULL, result, _rt_error_message, sizeof
            (_rt_error_message));
          rtmSetErrorStatus(simulink_experiment_test_typ_M, _rt_error_message);
        }
      }
    }

    hil_task_delete_all(simulink_experiment_test_typ_DW.HILInitialize_Card);
    hil_monitor_delete_all(simulink_experiment_test_typ_DW.HILInitialize_Card);
    hil_close(simulink_experiment_test_typ_DW.HILInitialize_Card);
    simulink_experiment_test_typ_DW.HILInitialize_Card = NULL;
  }
}

/*========================================================================*
 * Start of Classic call interface                                        *
 *========================================================================*/
void MdlOutputs(int_T tid)
{
  simulink_experiment_test_type1_output();
  UNUSED_PARAMETER(tid);
}

void MdlUpdate(int_T tid)
{
  simulink_experiment_test_type1_update();
  UNUSED_PARAMETER(tid);
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
  simulink_experiment_test_type1_initialize();
}

void MdlTerminate(void)
{
  simulink_experiment_test_type1_terminate();
}

/* Registration function */
RT_MODEL_simulink_experiment__T *simulink_experiment_test_type1(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)simulink_experiment_test_typ_M, 0,
                sizeof(RT_MODEL_simulink_experiment__T));

  /* Initialize timing info */
  {
    int_T *mdlTsMap =
      simulink_experiment_test_typ_M->Timing.sampleTimeTaskIDArray;
    mdlTsMap[0] = 0;

    /* polyspace +2 MISRA2012:D4.1 [Justified:Low] "simulink_experiment_test_typ_M points to
       static memory which is guaranteed to be non-NULL" */
    simulink_experiment_test_typ_M->Timing.sampleTimeTaskIDPtr = (&mdlTsMap[0]);
    simulink_experiment_test_typ_M->Timing.sampleTimes =
      (&simulink_experiment_test_typ_M->Timing.sampleTimesArray[0]);
    simulink_experiment_test_typ_M->Timing.offsetTimes =
      (&simulink_experiment_test_typ_M->Timing.offsetTimesArray[0]);

    /* task periods */
    simulink_experiment_test_typ_M->Timing.sampleTimes[0] = (0.002);

    /* task offsets */
    simulink_experiment_test_typ_M->Timing.offsetTimes[0] = (0.0);
  }

  rtmSetTPtr(simulink_experiment_test_typ_M,
             &simulink_experiment_test_typ_M->Timing.tArray[0]);

  {
    int_T *mdlSampleHits = simulink_experiment_test_typ_M->Timing.sampleHitArray;
    mdlSampleHits[0] = 1;
    simulink_experiment_test_typ_M->Timing.sampleHits = (&mdlSampleHits[0]);
  }

  rtmSetTFinal(simulink_experiment_test_typ_M, 90.0);
  simulink_experiment_test_typ_M->Timing.stepSize0 = 0.002;

  /* External mode info */
  simulink_experiment_test_typ_M->Sizes.checksums[0] = (1203433501U);
  simulink_experiment_test_typ_M->Sizes.checksums[1] = (1553460907U);
  simulink_experiment_test_typ_M->Sizes.checksums[2] = (3072448241U);
  simulink_experiment_test_typ_M->Sizes.checksums[3] = (1478629948U);

  {
    static const sysRanDType rtAlwaysEnabled = SUBSYS_RAN_BC_ENABLE;
    static RTWExtModeInfo rt_ExtModeInfo;
    static const sysRanDType *systemRan[1];
    simulink_experiment_test_typ_M->extModeInfo = (&rt_ExtModeInfo);
    rteiSetSubSystemActiveVectorAddresses(&rt_ExtModeInfo, systemRan);
    systemRan[0] = &rtAlwaysEnabled;
    rteiSetModelMappingInfoPtr(simulink_experiment_test_typ_M->extModeInfo,
      &simulink_experiment_test_typ_M->SpecialInfo.mappingInfo);
    rteiSetChecksumsPtr(simulink_experiment_test_typ_M->extModeInfo,
                        simulink_experiment_test_typ_M->Sizes.checksums);
    rteiSetTPtr(simulink_experiment_test_typ_M->extModeInfo, rtmGetTPtr
                (simulink_experiment_test_typ_M));
  }

  simulink_experiment_test_typ_M->solverInfoPtr =
    (&simulink_experiment_test_typ_M->solverInfo);
  simulink_experiment_test_typ_M->Timing.stepSize = (0.002);
  rtsiSetFixedStepSize(&simulink_experiment_test_typ_M->solverInfo, 0.002);
  rtsiSetSolverMode(&simulink_experiment_test_typ_M->solverInfo,
                    SOLVER_MODE_SINGLETASKING);

  /* block I/O */
  simulink_experiment_test_typ_M->blockIO = ((void *)
    &simulink_experiment_test_type_B);

  {
    simulink_experiment_test_type_B.HILReadEncoderTimebase = 0.0;
    simulink_experiment_test_type_B.HILReadAnalog = 0.0;
    simulink_experiment_test_type_B.BB01SensorGainmV = 0.0;
    simulink_experiment_test_type_B.Positionbias = 0.0;
    simulink_experiment_test_type_B.EncoderCalibrationradcount = 0.0;
    simulink_experiment_test_type_B.Bias = 0.0;
    simulink_experiment_test_type_B.MotorGainVV = 0.0;
    simulink_experiment_test_type_B.mtocm = 0.0;
    simulink_experiment_test_type_B.Gain = 0.0;
  }

  /* parameters */
  simulink_experiment_test_typ_M->defaultParam = ((real_T *)
    &simulink_experiment_test_type_P);

  /* states (dwork) */
  simulink_experiment_test_typ_M->dwork = ((void *)
    &simulink_experiment_test_typ_DW);
  (void) memset((void *)&simulink_experiment_test_typ_DW, 0,
                sizeof(DW_simulink_experiment_test_t_T));
  simulink_experiment_test_typ_DW.HILInitialize_AIMinimums[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AIMinimums[1] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AIMaximums[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AIMaximums[1] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOMinimums[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOMinimums[1] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOMaximums[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOMaximums[1] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_AOVoltages[1] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_FilterFrequency[0] = 0.0;
  simulink_experiment_test_typ_DW.HILInitialize_FilterFrequency[1] = 0.0;
  simulink_experiment_test_typ_DW.HILReadAnalog_Buffer = 0.0;

  /* data type transition information */
  {
    static DataTypeTransInfo dtInfo;
    (void) memset((char_T *) &dtInfo, 0,
                  sizeof(dtInfo));
    simulink_experiment_test_typ_M->SpecialInfo.mappingInfo = (&dtInfo);
    dtInfo.numDataTypes = 21;
    dtInfo.dataTypeSizes = &rtDataTypeSizes[0];
    dtInfo.dataTypeNames = &rtDataTypeNames[0];

    /* Block I/O transition table */
    dtInfo.BTransTable = &rtBTransTable;

    /* Parameters transition table */
    dtInfo.PTransTable = &rtPTransTable;
  }

  /* Initialize Sizes */
  simulink_experiment_test_typ_M->Sizes.numContStates = (0);/* Number of continuous states */
  simulink_experiment_test_typ_M->Sizes.numY = (0);/* Number of model outputs */
  simulink_experiment_test_typ_M->Sizes.numU = (0);/* Number of model inputs */
  simulink_experiment_test_typ_M->Sizes.sysDirFeedThru = (0);/* The model is not direct feedthrough */
  simulink_experiment_test_typ_M->Sizes.numSampTimes = (1);/* Number of sample times */
  simulink_experiment_test_typ_M->Sizes.numBlocks = (14);/* Number of blocks */
  simulink_experiment_test_typ_M->Sizes.numBlockIO = (9);/* Number of block outputs */
  simulink_experiment_test_typ_M->Sizes.numBlockPrms = (88);/* Sum of parameter "widths" */
  return simulink_experiment_test_typ_M;
}

/*========================================================================*
 * End of Classic call interface                                          *
 *========================================================================*/
