! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SortQualify ! Sort and qualify the L0 data
!=============================================================================

  USE MLSCommon, ONLY : r8
  USE MLSL1Common, ONLY : MAFinfo_T
  USE L0_sci_tbls, ONLY : SciMAF
  USE EngUtils, ONLY : NextEngMAF
  USE EngTbls, ONLY : EngMAF
  USE SciUtils, ONLY : NextSciMAF
  USE Calibration, ONLY : CalWin, MAFdata_T, UpdateCalVectors

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (MAFinfo_T) :: MAFinfo  ! Needed for L1BOA output

CONTAINS

!=============================================================================
  SUBROUTINE UpdateCalWindow (more_data, CalWinFull)
!=============================================================================

    USE MLSL1Config, ONLY : L1Config

    LOGICAL, INTENT (OUT) :: more_data
    LOGICAL, INTENT (OUT) :: CalWinFull

    INTEGER :: sci_MAFno, dif_MAFno, i
    INTEGER, SAVE :: prev_MAFno
    TYPE (MAFdata_T) :: EmptyMAFdata
    REAL(r8), SAVE :: prev_secTAI
    REAL :: MAF_dur, MIF_dur
    INTEGER :: nom_MIFs = 148   ! Nominal number of MIFs per MAF (TEST!)

!! Initialize empty MAF to "D"iscard:

    EmptyMAFdata%SciPkt%GHz_sw_pos = "D"
    EmptyMAFdata%SciPkt%THz_sw_pos = "D"
    DO i = 0, (SIZE (EmptyMAFdata%ChanType) - 1)
       EmptyMAFdata%ChanType(i)%FB = "D"
       EmptyMAFdata%ChanType(i)%MB = "D"
       EmptyMAFdata%ChanType(i)%WF = "D"
    ENDDO
    EmptyMAFdata%EMAF%MIFsPerMAF = nom_MIFs

    CALL NextEngMAF (more_data)

    CALL NextSciMAF (more_data)

    IF (.NOT. more_data) RETURN    !! Nothing more to do

    sci_MAFno = SciMAF(0)%MAFno

    IF (EngMAF%MAFno < sci_MAFno) THEN   ! Catch up to the Science
       DO
          CALL NextEngMAF (more_data)
          IF (EngMAF%MAFno == sci_MAFno .OR. .NOT. more_data) EXIT
       ENDDO
    ELSE IF (EngMAF%MAFno > sci_MAFno) THEN   ! Catch up to the Engineering
       DO
          CALL NextSciMAF (more_data)
          sci_MAFno = SciMAF(0)%MAFno
          IF (EngMAF%MAFno == sci_MAFno .OR. .NOT. more_data) EXIT
       ENDDO
    ENDIF

!! Use first 11 MIFs to determine integration time
!! NOTE: this requires all the 11 MIFs to be available
    
    MIF_dur = (SUM (SciMAF(1:11)%secTAI) - SUM (SciMAF(0:10)%secTAI)) / 10
    MAF_dur = MIF_dur * nom_MIFs   !Nominal duration of MAF

!! Eventually, hope to have a unique MAF counter to use instead since this
!! method is kind of kludgey

    IF (CalWin%current > 0) THEN
       dif_MAFno = sci_MAFno - prev_MAFno
       IF (dif_MAFno < 0) THEN      ! Rolled over
          dif_MAFno = NINT (((SciMAF(0)%secTAI - prev_secTAI - &
               4.0 * MIF_dur) / MAF_dur) + 0.5)
       ENDIF
    ELSE
       dif_MAFno = 1
    ENDIF

    prev_MAFno = sci_MAFno
    prev_secTAI = SciMAF(0)%secTAI
    IF (CalWin%current /= CalWin%size) THEN
       CalWin%current = CalWin%current + dif_MAFno
       IF (CalWin%current > CalWin%size) THEN  ! Beyond the end of the window
          CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, &
               (CalWin%current-CalWin%size), EmptyMAFdata)
          CalWin%current = CalWin%size
       ENDIF
    ELSE
       CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, dif_MAFno, EmptyMAFdata)
    ENDIF

    CalWin%MAFdata(CalWin%current)%SciPkt = SciMAF
    CalWin%MAFdata(CalWin%current)%EMAF = EngMAF

!! Update MAFinfo from central MAF

    MAFinfo%startTAI = CalWin%MAFdata(CalWin%central)%SciPkt(0)%secTAI
    MAFinfo%MIFsPerMAF = nom_MIFS    ! TEST!!!
    MAFinfo%integTime = MIF_dur      ! TEST!!!

!! Determine if CalWin is full

    IF (CalWin%current == CalWin%size .AND. &
         MAFinfo%startTAI >= L1Config%Input_TAI%startTime) THEN
       CalWinFull = .TRUE.
    ELSE
       CalWinFull = .FALSE.
    ENDIF

    more_data =  MAFinfo%startTAI <= L1Config%Input_TAI%endTime

  END SUBROUTINE UpdateCalWindow

  SUBROUTINE QualifyCurrentMAF

!! Qualify the Current MAF data in the cal window

    CHARACTER(len=2) :: GHz_sw_pos, THz_sw_pos, temp_type
    INTEGER :: current, MIF
    CHARACTER(len=2), PARAMETER :: TargetType = "T1" !! Primary target type

    current = CalWin%current

    DO MIF = 0, (SIZE (CalWin%MAFdata(current)%SciPkt) - 1) !! Check each packet

!! Initialize to "N"ot checked:

       CalWin%MAFdata(current)%ChanType(MIF)%FB = "N"
       CalWin%MAFdata(current)%ChanType(MIF)%MB = "N"
       CalWin%MAFdata(current)%ChanType(MIF)%WF = "N"

!! Rule #1: Data quality information:

       IF (.NOT. CalWin%MAFdata(current)%SciPkt(MIF)%CRC_good) THEN

          !! Discard all

          CalWin%MAFdata(current)%ChanType(MIF)%FB = "D"
          CalWin%MAFdata(current)%ChanType(MIF)%MB = "D"
          CalWin%MAFdata(current)%ChanType(MIF)%WF = "D"
          CYCLE                              !! All done for this packet
       ENDIF

!! Rule #2: User Input qualifications:

       !! Set the appropriate user input channels to "D"iscard

!! Rule #3: Disqualify based on engineering ("OFF" state, out-of-lock, etc.)

       !! Set the appropriate channels to "D"iscard

!! Rule #4: Check for "Z"ero data

       !! Set the appropriate channels to "Z"ero

!! Rule #5: Set Switching Mirror position

       !! Save current sw pos:
       
       GHz_sw_pos = CalWin%MAFdata(current)%SciPkt(MIF)%GHz_sw_pos
       THz_sw_pos = CalWin%MAFdata(current)%SciPkt(MIF)%THz_sw_pos

       !! GHz Module:

       !! Filterbanks 1 through 14

       WHERE (CalWin%MAFdata(current)%ChanType(MIF)%FB(:,1:14) == "N")
          CalWin%MAFdata(current)%ChanType(MIF)%FB(:,1:14) = Ghz_sw_pos
          !! NOTE: The THz module could be using FB 8!
       END WHERE
       WHERE (CalWin%MAFdata(current)%ChanType(MIF)%MB == "N")
          CalWin%MAFdata(current)%ChanType(MIF)%MB = Ghz_sw_pos
       END WHERE
       WHERE (CalWin%MAFdata(current)%ChanType(MIF)%WF == "N")
          CalWin%MAFdata(current)%ChanType(MIF)%WF = Ghz_sw_pos
       END WHERE

       !! "D"iscard incorrect Target Type:

       IF (Ghz_sw_pos /= TargetType) THEN
          WHERE (CalWin%MAFdata(current)%ChanType(MIF)%FB(:,1:14) == "T")
             CalWin%MAFdata(current)%ChanType(MIF)%FB(:,1:14) = "D"
          END WHERE
          WHERE (CalWin%MAFdata(current)%ChanType(MIF)%MB == "T")
             CalWin%MAFdata(current)%ChanType(MIF)%MB = "D"
          END WHERE
          WHERE (CalWin%MAFdata(current)%ChanType(MIF)%WF == "T")
             CalWin%MAFdata(current)%ChanType(MIF)%WF = "D"
          END WHERE
       ENDIF

       !! THz Module:

       !! Filterbanks 15 through 19

       WHERE (CalWin%MAFdata(current)%ChanType(MIF)%FB(:,15:19) == "N")
          CalWin%MAFdata(current)%ChanType(MIF)%FB(:,15:19) = Thz_sw_pos
          !! NOTE: The THz module could be using FB 8!
       END WHERE

SELECT CASE (MIF)   !!! TEST!!!

CASE (0:119)
   temp_type = "L"
CASE (123:128)
   temp_type = "T"
CASE (132:143)
   temp_type = "S"
CASE DEFAULT
   temp_type = "D"
END SELECT
CalWin%MAFdata(current)%ChanType(MIF)%FB = temp_type

    ENDDO

!!$print *, CalWin%MAFdata(current)%ChanType(0:149)%FB(1,1)
!!$print *, SciMAF(0)%MAFno, "Sw pos:", CalWin%MAFdata(current)%SciPkt%GHz_sw_pos
!! print *, CalWin%MAFdata(current)%EMAF%Eng(100)%mnemonic, CalWin%MAFdata(current)%EMAF%Eng(100)%value

  END SUBROUTINE QualifyCurrentMAF

  SUBROUTINE QualifyWindow

!! Qualify the calibration window for spikes, walls, etc.

  END SUBROUTINE QualifyWindow

!=============================================================================
  SUBROUTINE SortAndQualify (more_data, do_calib)
!=============================================================================

    LOGICAL, INTENT (OUT) :: more_data
    LOGICAL, INTENT (OUT) :: do_calib

    LOGICAL :: CalWinFull = .FALSE.

    do_calib = .FALSE.

    CALL UpdateCalWindow (more_data, CalWinFull)

    IF (.NOT. more_data) RETURN

    CALL QualifyCurrentMAF

    IF (CalWinFull) CALL QualifyWindow

    CALL UpdateCalVectors

    do_calib = CalWinFull  !! Calibrate if calib window is full

  END SUBROUTINE SortAndQualify

!=============================================================================
END MODULE SortQualify
!=============================================================================

! $Log$
! Revision 2.2  2001/03/05 19:54:41  perun
! Check TAI against input TAI
!
! Revision 2.1  2001/02/23 20:57:10  perun
! Version 0.5 commit
!
