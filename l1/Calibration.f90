! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Calibration ! Calibration data and routines
!=============================================================================

  USE MLSL1Common, ONLY: Chan_R_T, Chan_R8_T, FBchans, FBnum, MBchans, MBnum, &
       WFchans, WFnum, DACSchans, DACSnum, MaxMIFs, Bandwidth, deflt_gain, &
       deflt_zero, R8, tau
  USE L0_sci_tbls, ONLY: Sci_pkt_T
  USE EngTbls, ONLY : Eng_MAF_T
  USE Interpolation, ONLY : QuadInterpW

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(FBchans,FBnum)          ! standard filter banks
     CHARACTER(len=1) :: MB(MBchans,MBnum)          ! mid-band filter banks
     CHARACTER(len=1) :: WF(WFchans,WFnum)          ! wide filters
     CHARACTER(len=1) :: DACS(DACSchans,DACSnum)    ! DACS filters
  END TYPE Chan_type_T

  !! Bank Logical type

  TYPE BankLogical_T
     LOGICAL :: FB(FBnum)          ! standard filter banks
     LOGICAL :: MB(MBnum)          ! mid-band filter banks
     LOGICAL :: WF(WFnum)          ! wide filters
     LOGICAL :: DACS(DACSnum)      ! DACS filters
  END TYPE BankLogical_T

  !! Bank Integer type

  TYPE BankInt_T
     INTEGER :: FB(FBnum)          ! standard filter banks
     INTEGER :: MB(MBnum)          ! mid-band filter banks
     INTEGER :: WF(WFnum)          ! wide filters
     INTEGER :: DACS(DACSnum)      ! DACS filters
  END TYPE BankInt_T

  !! Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (Sci_pkt_T) :: SciPkt(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     TYPE (BankLogical_T) :: BankWall
     TYPE (BankLogical_T) :: Nominal   ! nominal switching flag
     TYPE (BankInt_T) :: BankCalInd(2) ! start & end indexes for calib
     TYPE (BankInt_T) :: WallMIF       ! MIF for start of wall
     INTEGER :: start_index, end_index  ! start & end indexes within cal vectors
     INTEGER :: last_MIF
     INTEGER :: BandSwitch(5)           ! band switch positions
  END TYPE MAFdata_T

  !! Calibration window:

  TYPE CalWin_T
     INTEGER :: size        ! size in MAFs
     INTEGER :: current     ! current index for new data
     INTEGER :: central     ! central index to calibrate
     TYPE (MAFdata_T), ALLOCATABLE, DIMENSION (:) :: MAFdata
  END TYPE CalWin_T

  TYPE (CalWin_T), TARGET :: CalWin
  TYPE (MAFdata_T), POINTER :: CurMAFdata

  !! Space and Target calibration vectors:

  INTEGER, PARAMETER :: max_cal_index = 2047

  TYPE (Chan_R8_T), TARGET :: space_time(0:max_cal_index)   ! Space times
  TYPE (Chan_R8_T), TARGET :: space_counts(0:max_cal_index) ! Space counts
  TYPE (Chan_R8_T), TARGET :: space_interp(0:MaxMIFS-1)     ! Space interpolate
  TYPE (Chan_R8_T), TARGET :: space_err(0:MaxMIFS-1)        ! Space error
  TYPE (Chan_R8_T), TARGET :: target_time(0:max_cal_index)  ! Target times
  TYPE (Chan_R8_T), TARGET :: target_counts(0:max_cal_index)! Target counts
  TYPE (Chan_R8_T), TARGET :: target_interp(0:MaxMIFS-1)    ! Target interpolate
  TYPE (Chan_R8_T), TARGET :: target_err(0:MaxMIFS-1)       ! Target error

  TYPE Chan_Int_T
     INTEGER :: FB(FBchans,FBnum)          ! standard filter banks
     INTEGER :: MB(MBchans,MBnum)          ! mid-band filter banks
     INTEGER :: WF(WFchans,WFnum)          ! wide filters
     INTEGER :: DACS(DACSchans,DACSnum)    ! DACS filters
  END TYPE Chan_Int_T

  !! Quality vectors:

  TYPE (Chan_Int_T), TARGET :: space_qual(0:max_cal_index)
  TYPE (Chan_Int_T), TARGET :: target_qual(0:max_cal_index)

  !! Weighting vectors:

  TYPE (Chan_R8_T) :: space_weight(0:max_cal_index)
  TYPE (Chan_R8_T) :: target_weight(0:max_cal_index)

  !! Limb time indexes:

  TYPE (Chan_R8_T), TARGET :: limb_time(0:max_cal_index)    ! Limb times
  TYPE (Chan_R8_T), TARGET :: limb_counts(0:max_cal_index)  ! Limb counts

  TYPE (Chan_Int_T), TARGET :: dum_qual(0:max_cal_index)    ! Dummy quality

  !! Default comVec vectors:

  TYPE ComVec_T
     REAL(r8), DIMENSION(:), ALLOCATABLE :: Space
     REAL(r8), DIMENSION(:), ALLOCATABLE :: Target
  END TYPE ComVec_T

  TYPE (ComVec_T), DIMENSION(:), ALLOCATABLE, TARGET :: GHz_comVec, THz_comVec

  !! Default errmul vectors:

  TYPE ErrMul_T
     REAL(r8) :: Space
     REAL(r8) :: Target
  END TYPE ErrMul_T

  TYPE (ErrMul_T), DIMENSION(:), ALLOCATABLE, TARGET :: GHz_errmul, THz_errmul

! chi square:

  TYPE (Chan_R_T) :: chi2

CONTAINS

  SUBROUTINE InitCalibWindow

    USE MLSL1Config, ONLY: L1Config, GHz_seq, THz_seq

    !! Initialize Calibration window data structures

    INTEGER :: i, last_MIF, last_cal_index, nVec, status
    INTEGER, POINTER :: window_MAFs, MIFsPerMAF
    REAL(r8) :: MIFno
    REAL(r8), DIMENSION(:,:), ALLOCATABLE :: cal_weight, cal_time
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: cal_qual
    REAL(r8), DIMENSION(:,:), ALLOCATABLE :: comVec
    REAL(r8), DIMENSION(:), ALLOCATABLE :: errmul

    window_MAFs => L1Config%Calib%CalWindow
    MIFsPerMAF => L1Config%Calib%MIFsPerMAF

    !! Allocate space for a calibration window's worth of science & eng data:

    ALLOCATE (CalWin%MAFdata(window_MAFs))

    !! Initialize indexes:

    CalWin%size = window_MAFs
    CalWin%current = 0        ! Indicates nothing in the window
    CalWin%central = window_MAFs / 2 + 1

    !! Initialize nominal data flags

    DO i = 1, window_MAFs
       CalWin%MAFdata(i)%nominal%FB = .TRUE.
       CalWin%MAFdata(i)%nominal%MB = .TRUE.
       CalWin%MAFdata(i)%nominal%WF = .TRUE.
       CalWin%MAFdata(i)%nominal%DACS = .TRUE.
    ENDDO

    !! Initialize space and target weighting vectors:

    DO i = 0, max_cal_index
       space_weight(i)%FB = 1.0d0
       space_weight(i)%MB = 1.0d0
       space_weight(i)%WF = 1.0d0
       space_weight(i)%DACS = 1.0d0
       target_weight(i)%FB = 1.0d0
       target_weight(i)%MB = 1.0d0
       target_weight(i)%WF = 1.0d0
       target_weight(i)%DACS = 1.0d0
    ENDDO

    !! Allocate default comVec and supporting vectors:

    last_MIF = MIFsPerMAF - 1
    last_cal_index = MIFsPerMAF * window_MAFs - 1

    ALLOCATE (GHz_comVec(0:last_MIF))
    ALLOCATE (THz_comVec(0:last_MIF))
    DO i = 0, last_MIF
       ALLOCATE (GHz_comVec(i)%Space(0:last_cal_index))
       ALLOCATE (GHz_comVec(i)%Target(0:last_cal_index))
       ALLOCATE (THz_comVec(i)%Space(0:last_cal_index))
       ALLOCATE (THz_comVec(i)%Target(0:last_cal_index))
    ENDDO
    ALLOCATE (cal_weight(0:last_MIF,window_MAFs))
    ALLOCATE (cal_qual(0:last_MIF,window_MAFs))
    ALLOCATE (cal_time(0:last_MIF,window_MAFs))
    ALLOCATE (comVec(0:last_MIF,0:last_cal_index))
    ALLOCATE (errmul(0:last_MIF))
    ALLOCATE (GHz_errmul(0:last_MIF))
    ALLOCATE (THz_errmul(0:last_MIF))

    MIFno = (CalWin%Central - 1) * MIFsPerMAF

    !! GHz vectors first:

    !! Space:

    cal_weight = RESHAPE (space_weight%FB(1,1), (/ MIFsPerMAF, window_MAFs /))

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, GHz_Seq, "S", &
       last_cal_index, last_MIF, window_MAFs, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       GHz_comVec(i)%Space = comVec(i,:)
       GHz_errmul(i)%Space = errmul(i)
    ENDDO

    !! Target:

    cal_weight = RESHAPE (target_weight%FB(1,1), (/ MIFsPerMAF, window_MAFs /))

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, GHz_Seq, "T", &
       last_cal_index, last_MIF, window_MAFs, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       GHz_comVec(i)%Target = comVec(i,:)
       GHz_errmul(i)%Target = errmul(i)
    ENDDO

    !! THz vectors next:

    !! Space:

    cal_weight = RESHAPE (space_weight%FB(1,1), (/ MIFsPerMAF, window_MAFs /))

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, THz_Seq, "S", &
       last_cal_index, last_MIF, window_MAFs, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       THz_comVec(i)%Space = comVec(i,:)
       THz_errmul(i)%Space = errmul(i)
    ENDDO

    !! Target:

    cal_weight = RESHAPE (target_weight%FB(1,1), (/ MIFsPerMAF, window_MAFs /))

    CALL CalcComVecs (cal_qual, cal_time, cal_weight, THz_Seq, "T", &
       last_cal_index, last_MIF, window_MAFs, MIFno, comVec, errmul)
    DO i = 0, last_MIF
       THz_comVec(i)%Target = comVec(i,:)
       THz_errmul(i)%Target = errmul(i)
    ENDDO

  END SUBROUTINE InitCalibWindow

  SUBROUTINE CalcComVecs (cal_qual, cal_time, cal_weight, Seq, cal_type, &
       last_cal_index, last_MIF, window_MAFs, MIFtime, comVec, errmul)

    INTEGER :: last_cal_index, last_MIF, window_MAFs
    INTEGER :: cal_qual(0:last_MIF,window_MAFs)
    REAL(r8) :: cal_time(0:last_MIF,window_MAFs)
    REAL(r8) :: cal_weight(0:last_MIF,window_MAFs)
    REAL(r8) :: comVec(0:last_MIF,0:last_cal_index)
    REAL(r8) :: errmul(0:last_MIF)
    CHARACTER(len=1) :: Seq(0:last_MIF), cal_type
    REAL(r8) :: MIFtime

    INTEGER :: i, MIFsPerMAF, nVec, status
    REAL(r8) :: MIFno

    MIFsPerMAF = last_MIF + 1
    cal_time = -1
    cal_qual = 0
    DO i = 0, last_MIF
       IF (Seq(i) == cal_type) THEN
          cal_time(i,1) = i
          cal_qual(i,1) = 1
       ENDIF
    ENDDO

    !! Fill in remainder of Calibration window:

    DO i = 2, window_MAFs
       cal_qual(:,i) = cal_qual(:,1)
       WHERE (cal_qual(:,i) == 1)
          cal_time(:,i) = cal_time(:,1) + MIFsPerMAF * (i-1)
       ENDWHERE
    ENDDO

    nVec = MIFsPerMAF * window_MAFs
    MIFno = MIFtime

    DO i = 0, last_MIF

       CALL QuadInterpW (RESHAPE (cal_time, (/ nVec /)), &
          RESHAPE (cal_weight, (/ nVec /)), RESHAPE (cal_qual, (/ nVec /)), &
          MIFno, nVec, comVec(i,:), errmul(i), status)

       !! Next MIF in central MAF:

       MIFno = MIFno + 1

    ENDDO

  END SUBROUTINE CalcComVecs

  SUBROUTINE SetCalVectors (cal_type, last_MIF, MIF_offset, windx, &
       time_index, cnts_index, qual_index)

!! Set the calibration vectors for a particular data type

    CHARACTER(LEN=1), INTENT (IN) :: cal_type 
    INTEGER, INTENT (IN) :: last_MIF
    INTEGER, INTENT (IN) :: MIF_offset
    INTEGER, INTENT (IN) :: windx
    TYPE (Chan_R8_T), DIMENSION(0:), INTENT(OUT) :: time_index
    TYPE (Chan_R8_T), DIMENSION(0:), INTENT(OUT) :: cnts_index
    TYPE (Chan_Int_T), DIMENSION(0:), INTENT (OUT) :: qual_index

    INTEGER :: MIF

    CurMAFdata => CalWin%MAFdata(windx)  ! point to current data

    DO MIF = 0, last_MIF

       WHERE (CurMAFdata%ChanType(MIF)%FB == cal_type)
          time_index(MIF_offset+MIF)%FB = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cnts_index(MIF_offset+MIF)%FB = &
               CurMAFdata%SciPkt(MIF)%FB
          qual_index(MIF_offset+MIF)%FB = 1    ! use for interpolation
       ELSEWHERE
          time_index(MIF_offset+MIF)%FB = -1   ! not available
          cnts_index(MIF_offset+MIF)%FB = 0    ! not available
          qual_index(MIF_offset+MIF)%FB = 0    ! don't use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%MB == cal_type)
          time_index(MIF_offset+MIF)%MB = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cnts_index(MIF_offset+MIF)%MB = &
               CurMAFdata%SciPkt(MIF)%MB
          qual_index(MIF_offset+MIF)%MB = 1    ! use for interpolation
       ELSEWHERE
          time_index(MIF_offset+MIF)%MB = -1   ! not available
          cnts_index(MIF_offset+MIF)%MB = 0    ! not available
          qual_index(MIF_offset+MIF)%MB = 0    ! don't use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%WF == cal_type)
          time_index(MIF_offset+MIF)%WF = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cnts_index(MIF_offset+MIF)%WF = &
               CurMAFdata%SciPkt(MIF)%WF
          qual_index(MIF_offset+MIF)%WF = 1    ! use for interpolation
       ELSEWHERE
          time_index(MIF_offset+MIF)%WF = -1   ! not available
          cnts_index(MIF_offset+MIF)%WF = 0    ! not available
          qual_index(MIF_offset+MIF)%WF = 0    ! don't use for interpolation
       END WHERE

       WHERE (CurMAFdata%ChanType(MIF)%DACS == cal_type)
          time_index(MIF_offset+MIF)%DACS = &
               CurMAFdata%SciPkt(MIF)%MIFno + MIF_offset
          cnts_index(MIF_offset+MIF)%DACS = &
               CurMAFdata%SciPkt(MIF)%DACS
          qual_index(MIF_offset+MIF)%DACS = 1    ! use for interpolation
       ELSEWHERE
          time_index(MIF_offset+MIF)%DACS = -1   ! not available
          cnts_index(MIF_offset+MIF)%DACS = 0    ! not available
          qual_index(MIF_offset+MIF)%DACS = 0    ! don't use for interpolation
       END WHERE

    ENDDO

  END SUBROUTINE SetCalVectors

  SUBROUTINE UpdateCalVectors

!! Update the Calibration Vectors used by the interpolator

    INTEGER :: windx, MIF_offset, last_MIF
    TYPE (Chan_R8_T), DIMENSION(:), POINTER :: time_index
    TYPE (Chan_R8_T), DIMENSION(:), POINTER :: cnts_index
    TYPE (Chan_Int_T), DIMENSION(:), POINTER :: qual_index

    DO windx = 1, CalWin%current

       last_MIF = CalWin%MAFdata(windx)%last_MIF
       MIF_offset = SUM(CalWin%MAFdata(1:windx)%EMAF%MIFsPerMAF) - &
            CalWin%MAFdata(1)%EMAF%MIFsPerMAF  ! MIF offset from beginning of
                                               ! Calibration Window

       !! Update Space vectors:

       time_index => space_time
       cnts_index => space_counts
       qual_index => space_qual
       CALL SetCalVectors ("S", last_MIF, MIF_offset, windx, &
            time_index, cnts_index, qual_index)

       !! Update Target vectors:

       time_index => target_time
       cnts_index => target_counts
       qual_index => target_qual
       CALL SetCalVectors ("T", last_MIF, MIF_offset, windx, &
            time_index, cnts_index, qual_index)

       !! Update Limb vectors:

       time_index => limb_time
       cnts_index => limb_counts
       qual_index => dum_qual
       CALL SetCalVectors ("L", last_MIF, MIF_offset, windx, &
            time_index, cnts_index, qual_index)

    ENDDO

  END SUBROUTINE UpdateCalVectors

  SUBROUTINE InterpCals (nVec, time, cal_time, cal_weight, cal_qual, &
       cal_counts, cal_interp, cal_err, GHz_comVec, THz_comVec, GHz_errmul, &
       THz_errmul, BankCalInd, cal_index)

    USE MLSL1Config, ONLY: L1Config

    INTEGER :: nVec, cal_index(2)
    REAL(r8) :: time, GHz_errmul, THz_errmul
    TYPE (Chan_R8_T) :: cal_time(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_counts(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_weight(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_interp
    TYPE (Chan_R8_T) :: cal_err
    TYPE (Chan_Int_T) :: cal_qual(0:nVec-1)
    REAL(r8) :: GHz_comVec(0:nVec-1), THz_comVec(0:nVec-1)
    TYPE (BankInt_T) :: BankCalInd(2) ! start & end indexes for calib

    INTEGER :: i, j, Istat, cal1, cal2, calen, calMAFs
    REAL(r8) :: errmul

    !! Previous arguments
    
    REAL(r8), SAVE :: tVecP(0:Max_cal_index), comVecP(0:Max_cal_index), &
         comVec(0:Max_cal_index)
    REAL(r8), SAVE :: timeP, errmulP
    INTEGER, SAVE :: qualVecP(0:Max_cal_index)
    INTEGER, SAVE :: statusP

    LOGICAL :: is_same
    INTEGER, PARAMETER :: MinCalMAFs = 3  ! minimum calibration MAFs needed

comVecP(0:nVec-1) = GHz_comVec   !! TEST!!!
errmulP = GHz_errmul   !! TEST!!!
is_same = .true.       !! TEST!!!
! Interpolate calibration values

    DO j = 1, FBnum
       cal1 = BankCalInd(1)%FB(j)
       cal2 = BankCalInd(2)%FB(j)
       calen = cal2 - cal1 + 1
       calMAFs = calen / L1Config%Calib%MIFsPerMAF

       IF (calMAFs < MinCalMAFs .OR. cal1 > cal_index(1) .OR. &
            cal2 < cal_index(2)) THEN
          cal_interp%FB(:,j) = 0.0
          cal_err%FB(:,j) = 0.0
       ELSE
          DO i = 1, FBchans

             !! Check for same as previous inputs:

             !! is_same = (time == timeP)
             !!is_same = (i .ne. 1) .or. (j .ne. 1)
!!$          IF (is_same) is_same = (ALL (cal_time%FB(i,j) == tVecP))
!!$          IF (is_same) is_same = (ALL (cal_qual%FB(i,j) == qualVecP))

             is_same = ALL (CalWin%MAFdata%nominal%FB(j))
             IF (is_same .AND. calen == nVec) THEN
                comVec = comVecP
                errmul = errmulP
                Istat = statusP
             ELSE
                CALL QuadInterpW (cal_time(cal1:cal2)%FB(i,j), &
                     cal_weight (cal1:cal2)%FB(i,j), &
                     cal_qual(cal1:cal2)%FB(i,j), time, calen, &
                     comVec(0:calen-1), errmul, Istat)

                !! Save for next call:

                tVecP = cal_time%FB(i,j)
                qualVecP = cal_qual%FB(i,j)
                timeP = time
                statusP = Istat
             ENDIF
             cal_interp%FB(i,j) = &
                  SUM (comVec(0:calen-1) * cal_counts(cal1:cal2)%FB(i,j))
             cal_err%FB(i,j) = errmul

          ENDDO
       ENDIF
    ENDDO

    DO j = 1, MBnum
       cal1 = BankCalInd(1)%MB(j)
       cal2 = BankCalInd(2)%MB(j)
       calen = cal2 - cal1 + 1
       DO i = 1, MBchans

          !! Check for same as previous inputs:
          
          !! is_same = (time == timeP)
!! is_same = (i .ne. 1) .or. (j .ne. 1)
!!$          IF (is_same) is_same = (ALL (cal_time%MB(i,j) == tVecP))
!!$          IF (is_same) is_same = (ALL (cal_qual%MB(i,j) == qualVecP))

          is_same = ALL (CalWin%MAFdata%nominal%MB(j))
          IF (is_same .AND. calen == nVec) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time(cal1:cal2)%MB(i,j), &
                  cal_weight(cal1:cal2)%MB(i,j), cal_qual(cal1:cal2)%MB(i,j), &
                  time, calen, comVec(0:calen-1), errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%MB(i,j)
             qualVecP = cal_qual%MB(i,j)
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%MB(i,j) = &
               SUM (comVec(0:calen-1) * cal_counts(0:calen-1)%MB(i,j))
          cal_err%MB(i,j) = errmul
       ENDDO
    ENDDO

    DO j = 1, WFnum
       cal1 = BankCalInd(1)%WF(j)
       cal2 = BankCalInd(2)%WF(j)
       calen = cal2 - cal1 + 1
       DO i = 1, WFchans

          !! Check for same as previous inputs:
          
          !! is_same = (time == timeP)
!! is_same = (i .ne. 1) .or. (j .ne. 1)
!!$          IF (is_same) is_same = (ALL (cal_time%WF(i,j) == tVecP))
!!$          IF (is_same) is_same = (ALL (cal_qual%WF(i,j) == qualVecP))

          is_same = ALL (CalWin%MAFdata%nominal%WF(j))
          IF (is_same .AND. calen == nVec) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time(cal1:cal2)%WF(i,j), &
                  cal_weight(cal1:cal2)%WF(i,j), cal_qual(cal1:cal2)%WF(i,j), &
                  time, calen, comVec(0:calen-1), errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%WF(i,j)
             qualVecP = cal_qual%WF(i,j)
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%WF(i,j) = &
               SUM (comVec(0:calen-1) * cal_counts(0:calen-1)%WF(i,j))
          cal_err%WF(i,j) = errmul
       ENDDO
    ENDDO
!!$
    IF (L1Config%Calib%CalibDACS) THEN

       DO j = 1, DACSnum
          cal1 = BankCalInd(1)%DACS(j)
          cal2 = BankCalInd(2)%DACS(j)
          calen = cal2 - cal1 + 1
          DO i = 1, DACSchans

             !! Check for same as previous inputs:

             !! is_same = (time == timeP)
             !! is_same = (i .ne. 1) .or. (j .ne. 1)
!!$          IF (is_same) is_same = (ALL (cal_time%DACS(i,j) == tVecP))
!!$          IF (is_same) is_same = (ALL (cal_qual%DACS(i,j) == qualVecP))

             is_same = ALL (CalWin%MAFdata%nominal%DACS(j))
             IF (is_same .AND. calen == nVec) THEN
                comVec = comVecP
                errmul = errmulP
                Istat = statusP
             ELSE
                CALL QuadInterpW (cal_time(cal1:cal2)%DACS(i,j), &
                     cal_weight(cal1:cal2)%DACS(i,j), &
                     cal_qual(cal1:cal2)%DACS(i,j), time, calen, &
                     comVec(0:calen-1), errmul, Istat)

                !! Save for next call:

                tVecP = cal_time%DACS(i,j)
                qualVecP = cal_qual%DACS(i,j)
                timeP = time
                statusP = Istat
             ENDIF
             cal_interp%DACS(i,j) = &
                  SUM (comVec(0:calen-1) * cal_counts(0:calen-1)%DACS(i,j))
             cal_err%DACS(i,j) = errmul
          ENDDO
       ENDDO

    ENDIF

  END SUBROUTINE InterpCals

  SUBROUTINE ChiSquare (space_counts, space_interp, nlast)

    INTEGER :: nlast
    TYPE (Chan_R8_T) :: space_counts(0:nlast), space_interp(0:nlast)
    INTEGER :: i, j, nspace, astat
    INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: nvec
    REAL(r8), DIMENSION(:), ALLOCATABLE, SAVE :: difspace, difzero
    REAL(r8) :: SumDifS2, SumDifZ2
    INTEGER, SAVE :: lastv = 0
    REAL :: space_tot
    INTEGER :: minmafs = 6

    IF (lastv /= nlast) THEN
       DEALLOCATE (nvec, STAT=astat)
       ALLOCATE (nvec(0:nlast))
       DEALLOCATE (difspace, STAT=astat)
       ALLOCATE (difspace(0:nlast))
       DEALLOCATE (difzero, STAT=astat)
       ALLOCATE (difzero(0:nlast))
       lastv = nlast
    ENDIF

    chi2%FB = 0.0      ! initial value
    DO j = 1, FBnum
       DO i = 1, FBchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%FB(i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%FB(i,j) - space_interp%FB(i,j)
             difzero = space_counts%FB(i,j) - deflt_zero%FB(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%FB(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
                     / (bandwidth%FB(i,j) * tau))
                chi2%FB(i,j) = chi2%FB(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    chi2%MB = 0.0      ! initial value
    DO j = 1, MBnum
       DO i = 1, MBchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%MB(i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%MB(i,j) - space_interp%MB(i,j)
             difzero = space_counts%MB(i,j) - deflt_zero%MB(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%MB(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) / ((SumDifZ2 / nspace) &
                     / (bandwidth%MB(i,j) * tau))
                chi2%MB(i,j) = chi2%MB(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

    chi2%WF = 0.0      ! initial value
    DO j = 1, WFnum
       DO i = 1, WFchans
          nvec = 0
          difspace = 0.0
          difzero = 0.0
          WHERE (space_counts%WF(i,j) /= 0.0)
             nvec = 1
             difspace = space_counts%WF(i,j) - space_interp%WF(i,j)
             difzero = space_counts%WF(i,j) - deflt_zero%WF(i,j)
          END WHERE
          nspace = SUM (nvec)
          IF (nspace >= minmafs) THEN
             SumDifS2 = SUM(difspace**2)
             SumDifZ2 = SUM(difzero**2)
             IF (SumDifS2 > 0.0 .AND. SumDifZ2 > 0.0) THEN
                chi2%WF(i,j) = (SumDifS2 / nspace - &
                     (SUM (difspace) / nspace)**2) /((SumDifZ2 / nspace) &
                     / (bandwidth%WF(i,j) * tau))
                chi2%WF(i,j) = chi2%WF(i,j) * nspace / (nspace - 1)
             ENDIF
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE ChiSquare

  SUBROUTINE Calibrate

    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Rad, ONLY: UpdateRadSignals

!! Calibrate the science data

    INTEGER :: time_index, start_index, end_index, windex
    INTEGER :: nVec, cal_index(2)
    INTEGER :: i, j, MIF_index
    REAL(r8) :: errmul, time, secs, oldsecs

    CHARACTER(len=8) :: date
    CHARACTER (len=10) :: dtime
    CHARACTER (len=5) :: zone
    INTEGER :: values(8)

PRINT *, 'calibrating...'

CALL DATE_AND_TIME (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
oldsecs = secs

    nVec = CalWin%MAFdata(CalWin%size)%end_index + 1
    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
    cal_index(1) = start_index - L1Config%Calib%MIFsPerMAF
    cal_index(2) = start_index + L1Config%Calib%MIFsPerMAF

    CALL UpdateRadSignals (CalWin%MAFdata(windex)%BandSwitch)

    DO time_index = start_index, end_index  ! for every MIF in the MAF

CALL DATE_AND_TIME (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
! if (time_index == end_index) then
! print *, time_index, dtime, (secs-oldsecs)
! else if (time_index == start_index) then
! oldsecs = secs
! endif
       time = time_index                     ! "Time" from start of cal window
       MIF_index = time_index - start_index  ! MIF # within the central MAF

       ! Space cals:

       CALL InterpCals (nVec, time, space_time, space_weight, space_qual, &
            space_counts, space_interp(MIF_index), space_err(MIF_index), &
            GHz_comVec(MIF_index)%Space, Thz_comVec(MIF_index)%Space, &
            GHz_errmul(MIF_index)%Space, Thz_errmul(MIF_index)%Space, &
            CalWin%MAFdata(windex)%BankCalInd, cal_index)

       ! Target cals:

       CALL InterpCals (nVec, time, target_time, target_weight, target_qual, &
            target_counts, target_interp(MIF_index), target_err(MIF_index), &
            GHz_comVec(MIF_index)%Target, THz_comVec(MIF_index)%Target, &
            GHz_errmul(MIF_index)%Target, THz_errmul(MIF_index)%Target, &
            CalWin%MAFdata(windex)%BankCalInd, cal_index)

    ENDDO

    CALL ChiSquare (space_counts(start_index:end_index), &
         space_interp(0:end_index-start_index), (end_index-start_index))

PRINT *, 'end calibrating...'

  END SUBROUTINE Calibrate

!=============================================================================
END MODULE Calibration
!=============================================================================

! $Log$
! Revision 2.3  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.2  2001/09/10 16:16:08  perun
! Changed ALLOCATABLE component to POINTER
!
! Revision 2.1  2001/02/23 18:50:29  perun
! Version 0.5 commit
!
