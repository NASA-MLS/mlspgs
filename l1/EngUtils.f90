! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE EngUtils   ! Engineering utilities
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Utils, ONLY: QNan, BigEndianStr, Finite

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

  FUNCTION Cal_freq (meas_freq, hi_freq, lo_freq, hi_cal, lo_cal) RESULT (freq)

    !! Calibrate the measured frequency against the calibration frequencies
    !! and their corresponding parameter values

    !--------Arguments--------!
    INTEGER, INTENT(IN) :: meas_freq  ! measured frequency
    INTEGER, INTENT(IN) :: hi_freq    ! "high" calibration frequency
    INTEGER, INTENT(IN) :: lo_freq    ! "low" calibration frequency
    REAL, INTENT(IN) :: hi_cal        ! "high" calibration parameter
    REAL, INTENT(IN) :: lo_cal        ! "low" calibration parameter

    REAL :: freq                      ! calculated frequency

    IF ((lo_freq == 0) .OR. (hi_freq == 0) .OR. (hi_freq == lo_freq)) THEN
       freq = QNan()
    ELSE
       freq = (REAL(meas_freq - lo_freq) / (hi_freq - lo_freq)) * &
            & (hi_cal - lo_cal) + lo_cal
    ENDIF

  END FUNCTION Cal_freq

  FUNCTION Therm_temp (rin) RESULT (T)

    !! Convert thermistor (ohms) to temperature (C) */

    !! The thermistor is a YSI 44906 in parallel with a 4.99 K ohm resistor.

    !! input:
    !!         rin: measured thermistor resistance (ohms)
    !! output:
    !!         T:   temperature (C)

    !--------Arguments--------!
    REAL, INTENT(IN) :: rin
    REAL :: T

    !----------Local vars----------!
    REAL :: rth_log

    !--------Parameters------------!
    REAL(r8), PARAMETER :: a = 1.286212e-3
    REAL(r8), PARAMETER :: b = 2.355213e-4
    REAL(r8), PARAMETER :: c = 9.826046e-8
    REAL(r8), PARAMETER :: d = 8.835732e-8

    REAL, PARAMETER :: abs_zero = -273.16
    REAL, PARAMETER :: rmax = 4990.0

    IF (.NOT. Finite (rin)) THEN  ! Nothing to be done
       T = Qnan()
    ELSE IF (rin >= rmax) THEN  ! over the maximum ohms
       T = -100.0               ! minimum possible temperature
    ELSE IF (rin <= 0.0) THEN   ! under the minimum ohms
       T = 500.0                ! maximum possible temperature
    ELSE
       rth_log = log ((rmax * rin) / (rmax - rin))
       T = 1.0 / (a + rth_log * (b + rth_log * (c + rth_log * d))) + abs_zero
    ENDIF

  END FUNCTION Therm_temp

  FUNCTION PRD_temp (rin, r0) RESULT (T)

    !! Convert PRD (ohms) to temperature (C)
    !! inputs:
    !!         rin: measured PRD resistance (ohms)
    !!         r0:  PRD resistance at 0 C (ohms)
    !! output:
    !!         T:   temperature (C)

    !--------Arguments--------!
    REAL, INTENT(IN) :: rin
    REAL, INTENT(IN) :: r0
    REAL :: T

    !----------Local vars----------!
    REAL(r8), PARAMETER :: a = 0.48945548411
    REAL(r8), PARAMETER :: b = 7.20107099888e-5
    REAL :: r

    r = rin * 500.0 / r0;
    T = a * (r - 500.0) / (1.0 - b * r)

  END FUNCTION PRD_temp

  FUNCTION Use_equation (mnemonic, tval) RESULT (val)

    ! Use equation to scale analog telemetry

    ! inputs:
    !          mnemonic: telemetry mnemonic
    !          tval: input telemetry value

    !--------Arguments--------!
    CHARACTER (LEN=*), INTENT(IN) :: mnemonic
    REAL, INTENT(IN)  :: tval

    REAL :: val

    val = tval

    ! check mnemonic for equation to use

    IF (INDEX (mnemonic, "_Bus") /= 0) THEN
       IF (INDEX (mnemonic, "_V") /= 0) THEN
          IF (INDEX (mnemonic, "Quiet") /= 0) THEN
             val = (tval - 0.3113) / 0.1183;   ! Quiet Bus Voltage
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = (tval - 0.2638) / 0.1192;   ! Device Bus Voltage
          ENDIF
       ELSE IF (INDEX (mnemonic, "_I") /= 0) THEN
          IF (INDEX (mnemonic, "Quiet") /= 0) THEN
             val = (tval - 0.317) / 0.162;   ! Quiet Bus Current
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = (tval - 0.228) / 1.691;   ! Device Bus Current
          ENDIF
       ENDIF
    ELSE IF (INDEX (mnemonic, "R4_RFE_Tripler_I") /= 0) THEN
       val = (2.5 - tval) / 49.9;   ! Tripler Current
    ELSE IF (INDEX (mnemonic, "R4_HarmMixer_V") /= 0) THEN
       IF (tval > 0.0) THEN
          val = 8.7225 * log (tval) - 29.277;   ! Harmonic Mixer Voltage
       ELSE
          val = QNan()
       ENDIF
    ENDIF

    ! will need to add equation for R4_RFE_Gunn_I later:
    !  24231 * V3out - 16822 * V2out + 4597 * Vout - 291.58

  END FUNCTION Use_equation


  SUBROUTINE ConvertEngCounts

    !! Convert the raw counts into engineering units

    USE EngTbls, ONLY: Eng_tbl, Riu_tbl, EngPkt, last_tlm_pt, temp_cal, &
     & vin_cal1, vin_cal2, prt1_cal1, prt1_cal2, prt2_cal1, prt2_cal2, &
     & ysi_cal1, ysi_cal2, Cal_const

    INTEGER :: i, n, pkt_no, riu_no, byte1
    REAL :: scale

    REAL, PARAMETER :: abs_zero = -273.16

! Extract the raw counts

    DO i = 1, SIZE(Riu_tbl)

       pkt_no = Riu_tbl(i)%pkt_no
       byte1 = Riu_tbl(i)%start_byte
       Riu_tbl(i)%id_word = BigEndianStr (EngPkt(pkt_no)(byte1:byte1+1))

       DO n = Riu_tbl(i)%first_pt, Riu_tbl(i)%last_pt

          byte1 = byte1 + 2    ! next word
          Eng_tbl(n)%counts = BigEndianStr (EngPkt(pkt_no)(byte1:byte1+1))

! Save the calibration counts

         IF (Eng_tbl(n)%cal_indx > 0) THEN
            Riu_tbl(i)%Cal_cnts(Eng_tbl(n)%cal_indx) = Eng_tbl(n)%counts
         ENDIF

       ENDDO

    ENDDO

! Convert the counts

    DO i = 1, last_tlm_pt

       riu_no = Eng_tbl(i)%riu_no
       scale = Eng_tbl(i)%scale

       Eng_tbl(i)%value = QNan()
       IF (Riu_tbl(riu_no)%id_word /= 0 .AND. &
            Eng_tbl(i)%counts > 0) THEN   ! Data available to convert

          SELECT CASE (Eng_tbl(i)%type)

          CASE ("Cal")

             IF (Eng_tbl(i)%cal_indx /= temp_cal) THEN   ! not a temperature
                Eng_tbl(i)%value = Eng_tbl(i)%counts
             ELSE
                Eng_tbl(i)%value = Cal_freq (Eng_tbl(i)%counts, &
                     & Riu_tbl(riu_no)%Cal_cnts(vin_cal1), &
                     & Riu_tbl(riu_no)%Cal_cnts(vin_cal2), &
                     & Cal_const(riu_no)%volts, 0.0) * 100.0 + abs_zero
             ENDIF

          CASE ("Analog")

             Eng_tbl(i)%value = Cal_freq (Eng_tbl(i)%counts, &
                  & Riu_tbl(riu_no)%Cal_cnts(vin_cal1), &
                  & Riu_tbl(riu_no)%Cal_cnts(vin_cal2), &
                  & Cal_const(riu_no)%volts, 0.0)

             IF (scale /= 0.0) THEN
                Eng_tbl(i)%value = scale * Eng_tbl(i)%value
             ELSE
                Eng_tbl(i)%value = Use_equation (Eng_tbl(i)%mnemonic, &
                     & Eng_tbl(i)%value)
             ENDIF

          CASE ("PRT-1")

             Eng_tbl(i)%value = PRD_temp (Cal_freq(Eng_tbl(i)%counts, &
                  & Riu_tbl(riu_no)%Cal_cnts(prt1_cal1), &
                  & Riu_tbl(riu_no)%Cal_cnts(prt1_cal2), &
                  & Cal_const(riu_no)%prt1_hi, Cal_const(riu_no)%prt1_low), &
                  & scale)

          CASE ("PRT-2")

             Eng_tbl(i)%value = PRD_temp (Cal_freq(Eng_tbl(i)%counts, &
                  & Riu_tbl(riu_no)%Cal_cnts(prt2_cal1), &
                  & Riu_tbl(riu_no)%Cal_cnts(prt2_cal2), &
                  & Cal_const(riu_no)%prt2_hi, Cal_const(riu_no)%prt2_low), &
                  & scale)

          CASE ("Therm")

             Eng_tbl(i)%value = Therm_temp (Cal_freq(Eng_tbl(i)%counts, &
                  & Riu_tbl(riu_no)%Cal_cnts(ysi_cal1), &
                  & Riu_tbl(riu_no)%Cal_cnts(ysi_cal2), &
                  & Cal_const(riu_no)%therm_hi, Cal_const(riu_no)%therm_low))

          CASE DEFAULT

             print *, "Unknown tlm type!"  !! use standard routine here

          END SELECT

       ENDIF

    ENDDO

  END SUBROUTINE ConvertEngCounts

!=============================================================================
  FUNCTION GetEngPkt () RESULT (OK)
!=============================================================================

    USE EngTbls, ONLY: ConvEngPkt
    USE OpenInit, ONLY: eng_unit, eng_recno

    LOGICAL OK
    INTEGER :: i, ios

    READ (eng_unit, rec=eng_recno, iostat=ios) ConvEngPkt

    IF (ios /= 0) THEN
       OK = .FALSE.
       RETURN         ! Nothing more to get
    ENDIF

    eng_recno = eng_recno + 1            ! position for next record

    OK = .TRUE.

  END FUNCTION GetEngPkt

!=============================================================================
  SUBROUTINE NextEngMAF (more_data)
!=============================================================================

    USE EngTbls, ONLY: ConvEngPkt, EngPkt, EngMAF, Eng_tbl
    USE L0Utils, ONLY: ReadL0Eng

    !! Get the next MAF's engineering data

    LOGICAL, INTENT (OUT) :: more_data

    LOGICAL :: OK

    more_data = .TRUE.


    CALL ReadL0Eng (EngPkt, EngMAF%MAFno, OK)

    IF (OK) THEN

       CALL ConvertEngCounts

       !! Save the required data for later use:

       EngMAF%Eng%value = Eng_tbl%value
       EngMAF%Eng%mnemonic = Eng_tbl%mnemonic

       !! Will Get the MIFs per MAF from 2.5 and later!!!

       EngMAF%MIFsPerMAF = 148

    ELSE

       more_data = .FALSE.

    ENDIF


  END SUBROUTINE NextEngMAF

!=============================================================================
END MODULE EngUtils
!=============================================================================

! $Log$
! Revision 2.2  2001/02/23 18:55:17  perun
! Version 0.5 commit
!
