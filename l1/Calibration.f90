! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Calibration ! Calibration data and routines
!=============================================================================

  USE MLSL1Common
  USE L0_sci_tbls
  USE EngTbls, ONLY : Eng_MAF_T
  USE Interpolation, ONLY : QuadInterpW

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER :: window_MAFs = 6       ! window size in MAFs (user input!)
  INTEGER :: max_window_MAFs = 10  ! window size in MAFs

  !! Channel type (D, L, S, T, Z):

  TYPE Chan_type_T
     CHARACTER(len=1) :: FB(FBchans,FBnum)    ! standard filter banks
     CHARACTER(len=1) :: MB(MBchans,MBnum)    ! mid-band filter banks
     CHARACTER(len=1) :: WF(WFchans,WFnum)    ! wide filters
  END TYPE Chan_type_T

  !! Science and Engineering data for 1 MAF:

  TYPE MAFdata_T
     TYPE (Sci_pkt_T) :: SciPkt(0:(MaxMIFs-1))
     TYPE (Eng_MAF_T) :: EMAF
     TYPE (Chan_type_T) :: ChanType(0:(MaxMIFs-1))
     INTEGER :: start_index, end_index  ! start & end indexes within cal vectors
     INTEGER :: last_mifno
  END TYPE MAFdata_T

  !! Calibration window:

  TYPE CalWin_T
     INTEGER :: size        ! size in MAFs
     INTEGER :: current     ! current index for new data
     INTEGER :: central     ! central index to calibrate
     TYPE (MAFdata_T), ALLOCATABLE, DIMENSION (:) :: MAFdata
  END TYPE CalWin_T

  TYPE (CalWin_T) :: CalWin

  !! Space and Target calibration vectors:

  INTEGER, PARAMETER :: max_cal_index = 2047

  TYPE Chan_R8_T
     REAL(r8) :: FB(FBchans,FBnum)    ! standard filter banks
     REAL(r8) :: MB(MBchans,MBnum)    ! mid-band filter banks
     REAL(r8) :: WF(WFchans,WFnum)    ! wide filters
  END TYPE Chan_R8_T

  TYPE (Chan_R8_T), TARGET :: space_time(0:max_cal_index)   ! Space times
  TYPE (Chan_R8_T), TARGET :: space_counts(0:max_cal_index) ! Space counts
  TYPE (Chan_R8_T), TARGET :: space_interp(0:MaxMIFS-1)     ! Space interpolate
  TYPE (Chan_R8_T), TARGET :: space_err(0:MaxMIFS-1)        ! Space error
  TYPE (Chan_R8_T), TARGET :: target_time(0:max_cal_index)  ! Target times
  TYPE (Chan_R8_T), TARGET :: target_counts(0:max_cal_index)! Target counts
  TYPE (Chan_R8_T), TARGET :: target_interp(0:MaxMIFS-1)    ! Target interpolate
  TYPE (Chan_R8_T), TARGET :: target_err(0:MaxMIFS-1)       ! Target error

  TYPE Chan_Int_T
     INTEGER :: FB(FBchans,FBnum)    ! standard filter banks
     INTEGER :: MB(MBchans,MBnum)    ! mid-band filter banks
     INTEGER :: WF(WFchans,WFnum)    ! wide filters
  END TYPE Chan_Int_T

  !! Quality vectors:

  TYPE (Chan_Int_T), TARGET :: space_qual(0:max_cal_index)
  TYPE (Chan_Int_T), TARGET :: target_qual(0:max_cal_index)

  !! Weighting vectors and sum**2:

  TYPE (Chan_R8_T) :: space_weight(0:max_cal_index)
  TYPE (Chan_R8_T) :: target_weight(0:max_cal_index)

  !! Limb and Zero time indexes:

  TYPE (Chan_R8_T), TARGET :: limb_time(0:max_cal_index)    ! Limb times
  TYPE (Chan_R8_T), TARGET :: limb_counts(0:max_cal_index)  ! Limb counts
  TYPE (Chan_R8_T), TARGET :: zero_time(0:max_cal_index)    ! Zero times
  TYPE (Chan_R8_T), TARGET :: zero_counts(0:max_cal_index)  ! Zero counts

  TYPE (Chan_Int_T), TARGET :: dum_qual(0:max_cal_index)    ! Dummy quality

CONTAINS

  SUBROUTINE InitCalibWindow

    !! Initialize Calibration window data structures

    INTEGER :: i

    !! Allocate space for a calibration window's worth of science & eng data:

    ALLOCATE (CalWin%MAFdata(window_MAFs))

    !! Initialize indexes:

    CalWin%size = window_MAFs
    CalWin%current = 0        ! Indicates nothing in the window
    CalWin%central = window_MAFs / 2 + 1

    !! Initialize space and target weighting vectors:

    DO i = 0, (SIZE (space_weight) - 1)
       space_weight(i)%FB = 1.0d0
       space_weight(i)%MB = 1.0d0
       space_weight(i)%WF = 1.0d0
       target_weight(i)%FB = 1.0d0
       target_weight(i)%MB = 1.0d0
       target_weight(i)%WF = 1.0d0
    ENDDO

  END SUBROUTINE InitCalibWindow

  SUBROUTINE SetCalVectors (cal_type, last_mifno, mif_offset, windx, &
       time_index, cnts_index, qual_index)

!! Set the calibration vectors for a particular data type

    CHARACTER(LEN=1), INTENT (IN) :: cal_type 
    INTEGER, INTENT (IN) :: last_mifno
    INTEGER, INTENT (IN) :: mif_offset
    INTEGER, INTENT (IN) :: windx
    TYPE (Chan_R8_T), DIMENSION(0:), INTENT(OUT) :: time_index
    TYPE (Chan_R8_T), DIMENSION(0:), INTENT(OUT) :: cnts_index
    TYPE (Chan_Int_T), DIMENSION(0:), INTENT (OUT) :: qual_index

    INTEGER :: MIF

    DO MIF = 0, last_mifno

       WHERE (CalWin%MAFdata(windx)%ChanType(MIF)%FB == cal_type)
          time_index(mif_offset+MIF)%FB = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%MIFno + mif_offset
          cnts_index(mif_offset+MIF)%FB = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%FB
          qual_index(mif_offset+MIF)%FB = 1    ! use for interpolation
       ELSEWHERE
          time_index(mif_offset+MIF)%FB = -1   ! not available
          cnts_index(mif_offset+MIF)%FB = 0    ! not available
          qual_index(mif_offset+MIF)%FB = 0    ! don't use for interpolation
       END WHERE

       WHERE (CalWin%MAFdata(windx)%ChanType(MIF)%MB == cal_type)
          time_index(mif_offset+MIF)%MB = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%MIFno + mif_offset
          cnts_index(mif_offset+MIF)%MB = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%MB
          qual_index(mif_offset+MIF)%MB = 1    ! use for interpolation
       ELSEWHERE
          time_index(mif_offset+MIF)%MB = -1   ! not available
          cnts_index(mif_offset+MIF)%MB = 0    ! not available
          qual_index(mif_offset+MIF)%MB = 0    ! don't use for interpolation
       END WHERE

       WHERE (CalWin%MAFdata(windx)%ChanType(MIF)%WF == cal_type)
          time_index(mif_offset+MIF)%WF = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%MIFno + mif_offset
          cnts_index(mif_offset+MIF)%WF = &
               CalWin%MAFdata(windx)%SciPkt(MIF)%WF
          qual_index(mif_offset+MIF)%WF = 1    ! use for interpolation
       ELSEWHERE
          time_index(mif_offset+MIF)%WF = -1   ! not available
          cnts_index(mif_offset+MIF)%WF = 0    ! not available
          qual_index(mif_offset+MIF)%WF = 0    ! don't use for interpolation
       END WHERE

    ENDDO

  END SUBROUTINE SetCalVectors

  SUBROUTINE UpdateCalVectors

!! Update the Calibration Vectors used by the interpolator

    INTEGER :: windx, mif_offset, last_mifno
    TYPE (Chan_R8_T), DIMENSION(:), POINTER :: time_index
    TYPE (Chan_R8_T), DIMENSION(:), POINTER :: cnts_index
    TYPE (Chan_Int_T), DIMENSION(:), POINTER :: qual_index

    DO windx = 1, CalWin%current

       last_mifno = CalWin%MAFdata(windx)%EMAF%MIFsPerMAF - 1
       mif_offset = SUM(CalWin%MAFdata(1:windx)%EMAF%MIFsPerMAF) - &
            CalWin%MAFdata(1)%EMAF%MIFsPerMAF  ! MIF offset from beginning of
                                               ! Calibration Window

       CalWin%MAFdata(windx)%last_mifno = last_mifno
       CalWin%MAFdata(windx)%start_index = mif_offset
       CalWin%MAFdata(windx)%end_index = &
            CalWin%MAFdata(windx)%start_index + last_mifno

       !! Update Space vectors:

       time_index => space_time
       cnts_index => space_counts
       qual_index => space_qual
       CALL SetCalVectors ("S", last_mifno, mif_offset, windx, &
            time_index, cnts_index, qual_index)

       !! Update Target vectors:

       time_index => target_time
       cnts_index => target_counts
       qual_index => target_qual
       CALL SetCalVectors ("T", last_mifno, mif_offset, windx, &
            time_index, cnts_index, qual_index)

       !! Update Limb vectors:

       time_index => limb_time
       cnts_index => limb_counts
       qual_index => dum_qual
       CALL SetCalVectors ("L", last_mifno, mif_offset, windx, &
            time_index, cnts_index, qual_index)

       !! Update Zero vectors:

       time_index => zero_time
       cnts_index => zero_counts
       qual_index => dum_qual
       CALL SetCalVectors ("Z", last_mifno, mif_offset, windx, &
            time_index, cnts_index, qual_index)

    ENDDO
    !! print *, INT(space_counts(0:mif_offset-1)%FB(1,1))

  END SUBROUTINE UpdateCalVectors

  SUBROUTINE InterpCals (nVec, time, cal_time, cal_weight, cal_qual, &
       cal_counts, cal_interp, cal_err)

    INTEGER :: nVec
    REAL(r8) :: time
    TYPE (Chan_R8_T) :: cal_time(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_counts(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_weight(0:nVec-1)
    TYPE (Chan_R8_T) :: cal_interp
    TYPE (Chan_R8_T) :: cal_err
    TYPE (Chan_Int_T) :: cal_qual(0:nVec-1)

    INTEGER :: i, j, Istat
    REAL(r8) :: errmul

    !! Previous arguments
    
    REAL(r8), ALLOCATABLE, DIMENSION(:), SAVE :: tVecP, comVecP, comVec
    REAL(r8), SAVE :: timeP, errmulP
    INTEGER, ALLOCATABLE, DIMENSION(:), SAVE :: qualVecP
    INTEGER, SAVE :: statusP
    INTEGER :: astat
    INTEGER, SAVE :: nsize = 0

    LOGICAL :: is_same

    IF (nsize /= nVec) THEN

       DEALLOCATE (comVec, STAT=astat)
       ALLOCATE (comVec(0:nVec-1))

       !! For next call:

       DEALLOCATE (tVecP, STAT=astat)
       ALLOCATE (tVecP(nVec))
       DEALLOCATE (comVecP, STAT=astat)
       ALLOCATE (comVecP(nVec))
       DEALLOCATE (qualVecP, STAT=astat)
       ALLOCATE (qualVecP(nVec))

       nsize = nVec
    ENDIF

! Interpolate calibration values

    DO j = 1, FBnum
       DO i = 1, FBchans

          !! Check for same as previous inputs:

          is_same = (time == timeP)
!!$          IF (is_same) is_same = (ALL (cal_time%FB(i,j) == tVecP))
          IF (is_same) is_same = (ALL (cal_qual%FB(i,j) == qualVecP))

          IF (is_same) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE

             CALL QuadInterpW (cal_time%FB(i,j), cal_weight%FB(i,j), &
                  cal_qual%FB(i,j), time, nVec, comVec, errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%FB(i,j)
             qualVecP = cal_qual%FB(i,j)
             comVecP = comVec
             errmulP = errmul
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%FB(i,j) = &
               SUM (comVec(0:nVec-1) * cal_counts(0:nVec-1)%FB(i,j))
          cal_err%FB(i,j) = errmul
       ENDDO
    ENDDO

    DO j = 1, MBnum
       DO i = 1, MBchans

          !! Check for same as previous inputs:
          
          is_same = (time == timeP)
!!$          IF (is_same) is_same = (ALL (cal_time%MB(i,j) == tVecP))
          IF (is_same) is_same = (ALL (cal_qual%MB(i,j) == qualVecP))

          IF (is_same) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time%MB(i,j), cal_weight%MB(i,j), &
                  cal_qual%MB(i,j), time, nVec, comVec, errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%MB(i,j)
             qualVecP = cal_qual%MB(i,j)
             comVecP = comVec
             errmulP = errmul
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%MB(i,j) = &
               SUM (comVec(0:nVec-1) * cal_counts(0:nVec-1)%MB(i,j))
          cal_err%MB(i,j) = errmul
       ENDDO
    ENDDO

    DO j = 1, WFnum
       DO i = 1, WFchans

          !! Check for same as previous inputs:
          
          is_same = (time == timeP)
!!$          IF (is_same) is_same = (ALL (cal_time%WF(i,j) == tVecP))
          IF (is_same) is_same = (ALL (cal_qual%WF(i,j) == qualVecP))

          IF (is_same) THEN
             comVec = comVecP
             errmul = errmulP
             Istat = statusP
          ELSE
             CALL QuadInterpW (cal_time%WF(i,j), cal_weight%WF(i,j), &
                  cal_qual%WF(i,j), time, nVec, comVec, errmul, Istat)

             !! Save for next call:
             
             tVecP = cal_time%WF(i,j)
             qualVecP = cal_qual%WF(i,j)
             comVecP = comVec
             errmulP = errmul
             timeP = time
             statusP = Istat
          ENDIF
          cal_interp%WF(i,j) = &
               SUM (comVec(0:nVec-1) * cal_counts(0:nVec-1)%WF(i,j))
          cal_err%WF(i,j) = errmul
       ENDDO
    ENDDO

  END SUBROUTINE InterpCals

  SUBROUTINE Calibrate

!! Calibrate the science data

    INTEGER :: time_index, start_index, end_index, windex
    INTEGER :: nVec, Istat
    INTEGER :: i, j, MIF_index
    REAL(r8) :: errmul, time, secs, oldsecs

    CHARACTER(len=8) :: date
    CHARACTER (len=10) :: dtime
    CHARACTER (len=5) :: zone
    INTEGER :: values(8)

print *, 'calibrating...'

CALL Date_and_time (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
oldsecs = secs

    nVec = CalWin%MAFdata(CalWin%size)%end_index + 1
    windex = CalWin%central
    start_index = CalWin%MAFdata(windex)%start_index  ! MIF 0
    end_index = CalWin%MAFdata(windex)%end_index      ! MIF max
!!$    print *, start_index, end_index, nVec

    DO time_index = start_index, end_index  ! for every MIF in the MAF

CALL Date_and_time (date, dtime, zone, values)
secs = values(5)*3600.0 + values(6)*60.0 + values(7) + values(8)*0.001
!!$print *, time_index, dtime, (secs-oldsecs)
oldsecs = secs
       time = time_index                     ! "Time" from start of cal window
       MIF_index = time_index - start_index  ! MIF # within the central MAF

       ! Space cals:

       CALL InterpCals (nVec, time, space_time, space_weight, space_qual, &
            space_counts, space_interp(MIF_index), space_err(MIF_index))

       ! Target cals:

       CALL InterpCals (nVec, time, target_time, target_weight, target_qual, &
            target_counts, target_interp(MIF_index), target_err(MIF_index))

    ENDDO

!!$print *, 'end calibrating...'
  END SUBROUTINE Calibrate

!=============================================================================
END MODULE Calibration
!=============================================================================

! $Log$
! Revision 2.1  2001/02/23 18:50:29  perun
! Version 0.5 commit
!
