! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE EngTbls   ! Level 1 engineering tables
!=============================================================================

  USE MLSStrings, ONLY: LinearSearchStringArray
  USE MLSL1Common, ONLY: R8

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: nrius = 32

  CHARACTER (LEN=4), PARAMETER, DIMENSION(nrius) :: Riu_str = (/ &
       "GM01", "GM02", "GM03", "GM04", "GM05", "GM06", "GM07", "GM08", &
       "GM09", "GM10", "GM11", "GM12", "GM13", "GM14", "GM15", "SM01", &
       "SM02", "SM03", "SM04", "SM05", "SM06", "SM07", "SM08", "SM09", &
       "SM10", "SM11", "SM12", "SM13", "SM14", "TM01", "TM02", "TM03" /)

  CHARACTER (LEN=9), PARAMETER, DIMENSION(9) :: Cal_type_str = (/ &
       "Vin_Cal1 ", "Vin_Cal2 ", "PRT1_Cal1", "PRT1_Cal2", "PRT2_Cal1", &
       "PRT2_Cal2", "YSI_Cal1 ", "YSI_Cal2 ", "Temp     " /)

! Indexes for the calibration buffer

  INTEGER, PARAMETER :: vin_cal1 = 1
  INTEGER, PARAMETER :: vin_cal2 = 2
  INTEGER, PARAMETER :: prt1_cal1 = 3
  INTEGER, PARAMETER :: prt1_cal2 = 4
  INTEGER, PARAMETER :: prt2_cal1 = 5
  INTEGER, PARAMETER :: prt2_cal2 = 6
  INTEGER, PARAMETER :: ysi_cal1 = 7
  INTEGER, PARAMETER :: ysi_cal2 = 8
  INTEGER, PARAMETER :: temp_cal = 9

! Calibration constants for each RIU

  TYPE Cal_Const_T
     CHARACTER (LEN=4) :: riu
     REAL :: volts
     REAL :: therm_hi
     REAL :: therm_low
     REAL :: prt1_hi
     REAL :: prt1_low
     REAL :: prt2_hi
     REAL :: prt2_low
  END TYPE Cal_Const_T

  TYPE (Cal_Const_T), PARAMETER :: Cal_const(nrius) = (/ &
       cal_const_t ("GM01", &
        5.0000, 4525.0, 996.0, 620.267, 480.358, 620.267, 480.358 ), &
       cal_const_t ("GM02", &
        4.9978, 4525.0, 996.0, 620.0, 480.0, 620.252, 480.387), &
       cal_const_t ("GM03", &
        4.9989, 4525.0, 996.0, 620.315, 480.404, 620.315, 480.404), &
       cal_const_t ("GM04", &
        4.9981, 4525.0, 996.0, 620.0, 480.0, 620.29, 480.391), &
       cal_const_t ("GM05", &
        4.9964, 4525.0, 996.0, 620.32, 480.438, 620.0, 480.0), &
       cal_const_t ("GM06", &
        5.0005, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("GM07", &
        4.9999, 4525.0, 996.0, 620.3, 480.36, 620.0, 480.0), &
       cal_const_t ("GM08", &
        4.9961, 4525.0, 996.0, 620.0, 480.0, 620.306, 480.555), &
       cal_const_t ("GM09", &
        4.9964, 4525.0, 996.0, 620.0, 480.0, 620.27, 480.361), &
       cal_const_t ("GM10", &
        4.9969, 4525.0, 996.0, 620.0, 480.0, 620.302, 480.384), &
       cal_const_t ("GM11", &
        4.9990, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("GM12", &
        4.9968, 4525.0, 996.0, 620.0, 480.0, 620.302, 480.393), &
       cal_const_t ("GM13", &
        4.9963, 4525.0, 996.0, 620.0, 480.0, 620.308, 480.455), &
       cal_const_t ("GM14", &
        4.9988, 4525.0, 996.0, 620.0, 480.0, 620.271, 480.36), &
       cal_const_t ("GM15", &
        4.9982, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM01", &
        4.9988, 4525.0, 996.0, 620.0, 480.0, 620.322, 480.331), &
       cal_const_t ("SM02", &
        4.9979, 4525.0, 996.0, 620.0, 480.0, 620.235, 480.357), &
       cal_const_t ("SM03", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM04", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM05", &
        4.9982, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM06", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM07", &
        4.9960, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM08", &
        4.9970, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM09", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM10", &
        4.9971, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM11", &
        4.9969, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM12", &
        4.9966, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM13", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("SM14", &
        5.0000, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("TM01", &
        4.9992, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0), &
       cal_const_t ("TM02", &
        4.9998, 4525.0, 996.0, 620.368, 480.388, 620.0, 480.0), &
       cal_const_t ("TM03", &
        4.9960, 4525.0, 996.0, 620.0, 480.0, 620.0, 480.0) /)

! Engineering telemetry table

  INTEGER, PARAMETER :: maxtlm = 516
  INTEGER :: last_tlm_pt

  TYPE Eng_Tbl_T
     INTEGER :: cal_indx
     REAL :: scale
     INTEGER :: riu_no
     INTEGER :: counts
     REAL :: value
     CHARACTER (LEN=16) :: mnemonic
     CHARACTER (LEN=7) :: type
  END TYPE Eng_Tbl_T

  TYPE (Eng_Tbl_T) :: Eng_tbl(maxtlm)

  TYPE Eng_T
     REAL :: value
     CHARACTER (LEN=16) :: mnemonic
  END TYPE Eng_T

  TYPE Eng_MAF_T
     REAL(r8) :: secTAI
     INTEGER :: MAFno
     INTEGER :: MIFsPerMAF
     INTEGER :: TotalMAF
     CHARACTER(LEN=1) :: GSM_Side   ! GHz Switching Mirror
     CHARACTER(LEN=1) :: ASE_Side   ! Antenna Scan Electronics
     TYPE (Eng_T) :: Eng(maxtlm)
  END TYPE Eng_MAF_T

  TYPE (Eng_MAF_T) :: EngMAF

  ! Calibration target indexes

  TYPE CalTgtIndx_T
     INTEGER :: GHzAmb(20)
     INTEGER :: GHzCntl(20)
     INTEGER :: THzAmb(14)
  END TYPE CalTgtIndx_T

  TYPE (CalTgtIndx_T) :: CalTgtIndx

  ! Reflector temperature indexes

  TYPE ReflecIndx_T
     INTEGER :: Pri(9)
     INTEGER :: Sec(5)
     INTEGER :: Ter(3)
  END TYPE ReflecIndx_T

  TYPE (ReflecIndx_T) :: ReflecIndx

  TYPE Reflec_T
     REAL :: Pri, Sec, Ter
  END TYPE Reflec_T
  TYPE (Reflec_T) :: Reflec

! RIU table

  TYPE Riu_Tbl_T
     INTEGER :: pkt_no
     INTEGER :: start_byte
     INTEGER :: first_pt
     INTEGER :: last_pt
     INTEGER :: id_word
     INTEGER :: Cal_cnts(9)
  END TYPE Riu_Tbl_T

  ! Initialize table with packet numbers and start bytes within the packet:

  TYPE (Riu_Tbl_T) :: Riu_tbl(nrius) = (/ &
       Riu_tbl_t (1,  73, 0, 0, 0, 0), &
       Riu_tbl_t (1,  73, 0, 0, 0, 0), &
       Riu_tbl_t (3,  21, 0, 0, 0, 0), &
       Riu_tbl_t (3,  21, 0, 0, 0, 0), &
       Riu_tbl_t (2,  47, 0, 0, 0, 0), &
       Riu_tbl_t (5, 141, 0, 0, 0, 0), &
       Riu_tbl_t (2, 103, 0, 0, 0, 0), &
       Riu_tbl_t (2, 161, 0, 0, 0, 0), &
       Riu_tbl_t (2, 205, 0, 0, 0, 0), &
       Riu_tbl_t (4, 131, 0, 0, 0, 0), &
       Riu_tbl_t (4, 179, 0, 0, 0, 0), &
       Riu_tbl_t (3,  79, 0, 0, 0, 0), &
       Riu_tbl_t (3, 127, 0, 0, 0, 0), &
       Riu_tbl_t (3, 175, 0, 0, 0, 0), &
       Riu_tbl_t (2,  17, 0, 0, 0, 0), &
       Riu_tbl_t (4,  19, 0, 0, 0, 0), &
       Riu_tbl_t (5,  21, 0, 0, 0, 0), &
       Riu_tbl_t (4,  65, 0, 0, 0, 0), &
       Riu_tbl_t (4,  75, 0, 0, 0, 0), &
       Riu_tbl_t (4,  85, 0, 0, 0, 0), &
       Riu_tbl_t (5,  67, 0, 0, 0, 0), &
       Riu_tbl_t (5,  77, 0, 0, 0, 0), &
       Riu_tbl_t (5,  95, 0, 0, 0, 0), &
       Riu_tbl_t (5, 113, 0, 0, 0, 0), &
       Riu_tbl_t (5, 123, 0, 0, 0, 0), &
       Riu_tbl_t (6,  15, 0, 0, 0, 0), &
       Riu_tbl_t (6,  37, 0, 0, 0, 0), &
       Riu_tbl_t (6,  55, 0, 0, 0, 0), &
       Riu_tbl_t (6,  65, 0, 0, 0, 0), &
       Riu_tbl_t (1, 129, 0, 0, 0, 0), &
       Riu_tbl_t (1, 159, 0, 0, 0, 0), &
       Riu_tbl_t (4, 103, 0, 0, 0, 0) /)

! Survival telemetry table

  TYPE Srvl_Tbl_T
     REAL :: r0
     CHARACTER (LEN=7) :: type
     CHARACTER (LEN=16) :: mnemonic
  END TYPE Srvl_Tbl_T

  TYPE (Srvl_Tbl_T) :: Srvl_tbl(16)

! Engineering packets buffer to hold data for 1 MAF

  CHARACTER(LEN=256) :: EngPkt(6)

! Converted engineering packet (From the EM):

  CHARACTER(LEN=678) :: ConvEngPkt

CONTAINS

  SUBROUTINE Load_Eng_tbls (eng_unit, ios)

    INTEGER, INTENT (IN) :: eng_unit
    INTEGER, INTENT (OUT) :: ios

    CHARACTER (LEN=16) :: mnem
    CHARACTER (LEN=7) :: type
    CHARACTER (LEN=4) :: riu
    CHARACTER (LEN=80) :: line
    CHARACTER (LEN=2) :: dummy

    INTEGER :: i, order
    INTEGER :: GHzAmbIndx = 0
    INTEGER :: GHzCntlIndx = 0
    INTEGER :: THzAmbIndx = 0
    INTEGER :: PriReflecIndx = 0
    INTEGER :: SecReflecIndx = 0
    INTEGER :: TerReflecIndx = 0
    INTEGER :: riu_chan, riu_no, old_riu

    REAL :: scale

!! Read Survival table

    i = 0
    DO
       read (eng_unit, '(a)', iostat=ios) line
       IF (ios /= 0) EXIT
       read (line, *, iostat=ios) dummy, mnem, type, dummy, riu_chan, scale

       IF (ios == 0) THEN
          IF (INDEX (TRIM(type), "Digital") == 0) THEN  ! NOT Digital
             i = i + 1
             srvl_tbl(i)%mnemonic = mnem
             srvl_tbl(i)%type = type
             srvl_tbl(i)%r0 = scale
             IF (i == SIZE(srvl_tbl)) EXIT
          ENDIF
       ENDIF
    ENDDO

!! Read ENG table

    old_riu = -1
    i = 0
    DO
       read (eng_unit, '(a)', iostat=ios) line
       IF (ios /= 0) EXIT
       read (line, *, iostat=ios) order, mnem, type, riu, riu_chan, scale

       IF (ios == 0) THEN
          i = i + 1
          riu_no = LinearSearchStringArray (Riu_str, riu)
          Eng_tbl(i)%type = type
          Eng_tbl(i)%mnemonic = mnem
          Eng_tbl(i)%scale = scale
          Eng_tbl(i)%riu_no = riu_no
          IF (Eng_tbl(i)%type == "Cal") THEN
             Eng_tbl(i)%cal_indx = LinearSearchStringArray (cal_type_str, &
                  mnem, testSubstring=.TRUE., listInString=.TRUE.)
          ELSE
             Eng_tbl(i)%cal_indx = 0
          END IF
          IF (riu_no /= old_riu) THEN
             Riu_tbl(riu_no)%first_pt = i
          END IF
          Riu_tbl(riu_no)%last_pt = i
          old_riu = riu_no
          last_tlm_pt = i

          !! Save Cal Target indexes:

          IF (INDEX (mnem, "GHzAmbCalTgt") /= 0) THEN
             GHzAmbIndx = GHzAmbIndx + 1
             CalTgtIndx%GHzAmb(GHzAmbIndx) = i
          ELSE IF (INDEX (mnem, "GHzCntlCalTgt") /= 0) THEN
             GHzCntlIndx = GHzCntlIndx + 1
             CalTgtIndx%GHzCntl(GHzCntlIndx) = i
          ELSE IF (INDEX (mnem, "THzAmbCalTgt") /= 0) THEN
             THzAmbIndx = THzAmbIndx + 1
             CalTgtIndx%THzAmb(THzAmbIndx) = i
          ENDIF

          !! Save Reflector indexes:

          IF (INDEX (mnem, "Pri_Reflec") /= 0) THEN
             PriReflecIndx = PriReflecIndx + 1
             ReflecIndx%Pri(PriReflecIndx) = i
          ELSE IF (INDEX (mnem, "Sec_Reflec") /= 0) THEN
             SecReflecIndx = SecReflecIndx + 1
             ReflecIndx%Sec(SecReflecIndx) = i
          ELSE IF (INDEX (mnem, "Ter_Reflec") /= 0) THEN
             TerReflecIndx = TerReflecIndx + 1
             ReflecIndx%Ter(TerReflecIndx) = i
          ENDIF

       ENDIF

    ENDDO

    IF (i > 1 .AND. i <= SIZE (eng_tbl)) THEN
       ios = 0
    ELSE
       ios = -1
    ENDIF

  END Subroutine Load_Eng_tbls

END MODULE EngTbls

! $Log$
! Revision 2.7  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.6  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.5  2002/08/27 20:34:43  pwagner
! Back to using MLSStrings, this time from lib
!
! Revision 2.4  2002/08/07 18:54:22  jdone
! Changes due to l1/MLSL1Strings.f90
!
! Revision 2.3  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.2  2001/02/23 18:56:17  perun
! Version 0.5 commit
!
