! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE EngUtils   ! Engineering utilities
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Utils, ONLY: QNan, BigEndianStr, Finite

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: NextEngMAF

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

    r = rin * 500.0 / r0
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
             val = (tval - 0.3113) / 0.1183   ! Quiet Bus Voltage
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = (tval - 0.2638) / 0.1192   ! Device Bus Voltage
          ENDIF
       ELSE IF (INDEX (mnemonic, "_I") /= 0) THEN
          IF (INDEX (mnemonic, "Quiet") /= 0) THEN
             val = (-0.156313 + SQRT (0.156313 * 0.156313 - 4.0 * 0.00026955 * &
                  (0.361124 - tval))) / (2 * 0.00026955)  ! Quiet Bus Current
          ELSE IF (INDEX (mnemonic, "Device") /= 0) THEN
             val = 1.36 * (tval - 0.4852) / (1.3635 - 0.4852) ! Device Bus I
           ENDIF
       ENDIF
    ELSE IF (INDEX (mnemonic, "R4_IF_voltage") /= 0) THEN
       IF (tval > 0.0) THEN
          val = 8.4459 * log (tval) - 47.755  ! R4 IF Voltage
       ELSE
          val = QNan()
       ENDIF
    ENDIF

  END FUNCTION Use_equation


  SUBROUTINE ConvertEngCounts (GMAB_ON)

    !! Convert the raw counts into engineering units

    USE EngTbls, ONLY: Eng_tbl, Riu_tbl, EngPkt, last_tlm_pt, temp_cal, &
     & vin_cal1, vin_cal2, prt1_cal1, prt1_cal2, prt2_cal1, prt2_cal2, &
     & ysi_cal1, ysi_cal2, Cal_const

    LOGICAL, DIMENSION(:) :: GMAB_ON

    INTEGER :: i, n, dn, pkt_no, riu_no, byte1
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
            IF (Eng_tbl(n)%cal_indx == prt1_cal1) THEN  ! save also as PRT-2's
               Riu_tbl(i)%Cal_cnts(prt2_cal1) = Eng_tbl(n)%counts
            ELSE IF (Eng_tbl(n)%cal_indx == prt1_cal2) THEN
               Riu_tbl(i)%Cal_cnts(prt2_cal2) = Eng_tbl(n)%counts
            ENDIF
         ENDIF

       ENDDO

    ENDDO

! Convert the counts

    DO i = 1, last_tlm_pt

       riu_no = Eng_tbl(i)%riu_no
       scale = Eng_tbl(i)%scale

       Eng_tbl(i)%value = QNan()

       ! Skip the always bad THz sensor

       IF ((INDEX (Eng_tbl(i)%mnemonic, "THzAmbCalTgt_T2") /= 0) .OR. &
            (INDEX (Eng_tbl(i)%mnemonic, "THzAmbCalTgtRT2") /= 0)) CYCLE

       IF (riu_no >= 1 .AND. riu_no <= 4) THEN
          IF (.NOT. GMAB_ON(riu_no)) CYCLE     ! GMriu_no is NOT ON
          IF (MOD (riu_no, 2) == 0) THEN
             IF (GMAB_ON(riu_no-1)) CYCLE      ! Other related GMriu_no is ON
          ELSE
             IF (GMAB_ON(riu_no+1)) CYCLE      ! Other related GMriu_no is ON
          ENDIF
       ENDIF

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

          CASE ("VIN")

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

             dn = Eng_tbl(i)%counts
             IF (dn > 63000) dn = dn - 65536  ! adjust if DN is too high
             Eng_tbl(i)%value = PRD_temp (Cal_freq(dn, &
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

          CASE ("YSI")

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
  SUBROUTINE NextEngMAF (more_data)
!=============================================================================

    USE EngTbls, ONLY: EngPkt, EngMAF, Eng_tbl
    USE L0Utils, ONLY: ReadL0Eng

    !! Get the next MAF's engineering data

    LOGICAL, INTENT (OUT) :: more_data

    LOGICAL :: OK, GMAB_ON(4)
    INTEGER :: maskbit3, maskbit7
    DATA maskbit3 / z'8' /
    DATA maskbit7 / z'80' /

    more_data = .TRUE.

    CALL ReadL0Eng (EngPkt, EngMAF%MAFno, EngMAF%TotalMAF, EngMAF%MIFsPerMAF, &
         EngMAF%secTAI, OK)

    IF (OK) THEN

       ! Determine GM01 to GM04 A/B ON/OFF states:

       GMAB_ON(1) =(IAND (ICHAR (EngPkt(6)(103:103)), maskbit7) == maskbit7)
       GMAB_ON(2) =(IAND (ICHAR (EngPkt(6)(103:103)), maskbit3) == maskbit3)
       GMAB_ON(3) =(IAND (ICHAR (EngPkt(6)(109:109)), maskbit3) == maskbit3)
       GMAB_ON(4) =(IAND (ICHAR (EngPkt(6)(110:110)), maskbit7) == maskbit7)

       ! Determine which side ASE (GM01/GM02) is on:

       IF (GMAB_ON(1) .AND. .NOT. GMAB_ON(2)) THEN
          EngMAF%ASE_Side = "A"
       ELSE IF (.NOT. GMAB_ON(1) .AND. GMAB_ON(2)) THEN
          EngMAF%ASE_Side = "B"
       ELSE
          EngMAF%ASE_Side = "U"
       ENDIF

       ! Determine which side GSM (GM03/GM04) is on:

       IF (GMAB_ON(3) .AND. .NOT. GMAB_ON(4)) THEN
          EngMAF%GSM_Side = "A"
       ELSE IF (.NOT. GMAB_ON(3) .AND. GMAB_ON(4)) THEN
          EngMAF%GSM_Side = "B"
       ELSE
          EngMAF%GSM_Side = "U"
       ENDIF

       ! Convert engineering counts

       CALL ConvertEngCounts (GMAB_ON)

       !! Write eng data to file:

!!$       WRITE (L1BFileInfo%EngId, iostat=ios) EngPkt
!!$       WRITE (L1BFileInfo%EngId, iostat=ios) eng_tbl%value

       !! Save the required data for later use:

       EngMAF%Eng%value = Eng_tbl%value
       EngMAF%Eng%mnemonic = Eng_tbl%mnemonic

       !! Already have the current MIFsPerMAF

       !! EngMAF%MIFsPerMAF = L1Config%Calib%MIFsPerMAF

    ELSE

       more_data = .FALSE.

    ENDIF


  END SUBROUTINE NextEngMAF

!=============================================================================
END MODULE EngUtils
!=============================================================================

! $Log$
! Revision 2.7  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.6  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2001/02/23 18:55:17  perun
! Version 0.5 commit
!
