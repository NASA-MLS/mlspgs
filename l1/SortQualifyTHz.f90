! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SortQualifyTHz ! Sort and qualify the L0 data for the THz module
!=============================================================================

  USE L0_sci_tbls, ONLY: THzSciMAF
  USE THzCalibration, ONLY: CalBuf, MAFdata_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SortAndQualifyTHz

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()

CONTAINS

!=============================================================================
  SUBROUTINE FillCalData
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Common, ONLY: MAFinfo, absZero_C
    USE EngTbls, ONLY: EngMAF, CalTgtIndx 
    USE EngUtils, ONLY: NextEngMAF
    USE SciUtils, ONLY: NextSciMAF
    USE MLSL1Utils, ONLY : GetIndexedAvg

    LOGICAL :: more_data

    INTEGER :: sci_MAFno, MIFsPerMAF, CalMAFs, CalMAFno
    REAL :: MAF_dur, MIF_dur

    MIFsPerMAF = L1Config%Calib%MIFsPerMAF
    MIF_dur = L1Config%Calib%MIF_duration
    MAF_dur = MIF_dur * MIFsPerMAF   !Nominal duration of MAF
    calMAFs = NINT ((L1Config%Expanded_TAI%endTime - &
         L1Config%Expanded_TAI%startTime) / MAF_dur) + 10

    ALLOCATE (CalBuf%MAFdata(calMAFs))

    CalMAFno = 0

    DO

       CALL NextEngMAF (more_data)

       IF (.NOT. more_data) EXIT    !! Nothing more to do

       CALL NextSciMAF (more_data)

       IF (.NOT. more_data) EXIT    !! Nothing more to do

       sci_MAFno = THzSciMAF(0)%MAFno

       DO

          IF (EngMAF%MAFno == sci_MAFno) EXIT   ! Nothing more to read

          IF (EngMAF%secTAI <= THzSciMAF(2)%secTAI) THEN       ! Catch the Sci

             CALL NextEngMAF (more_data)
             IF (.NOT. more_data) RETURN    !! Nothing more to do

          ELSE IF (EngMAF%secTAI > THzSciMAF(2)%secTAI) THEN   ! Catch the Eng

             CALL NextSciMAF (more_data)
             IF (.NOT. more_data) RETURN    !! Nothing more to do
             sci_MAFno = THzSciMAF(0)%MAFno

          ENDIF

       ENDDO

print *, "SCI/ENG MAF: ", sci_MAFno, EngMAF%MAFno
       CalMAFno = CalMAFno + 1
       CurMAFdata => CalBuf%MAFdata(CalMAFno)
       CurMAFdata%SciMIF = THzSciMAF
       CurMAFdata%EMAF = EngMAF
       CurMAFdata%last_MIF = MIFsPerMAF - 1
       CurMAFdata%BandSwitch(4:5) =  &
            CurMAFdata%SciMIF(CurMAFdata%last_MIF)%BandSwitch
       CurMAFdata%CalTgtTemp = &
            GetIndexedAvg (EngMAF%eng%value, CalTgtIndx%THzAmb) - absZero_C

       more_data =  MAFinfo%startTAI <= L1Config%Input_TAI%endTime

    ENDDO
    CalBuf%MAFs = CalMAFno

  END SUBROUTINE FillCalData

!=============================================================================
  SUBROUTINE QualifyEachMAF
!=============================================================================

    USE MLSL1Config, ONLY: THz_seq, THz_seq_use

    INTEGER :: MAF, MIF

    CHARACTER(len=1) :: SwMirPos
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"

    print *, 'qualifying all the MAFs'

    DO MAF = 1, CalBuf%MAFs

       CurMAFdata => CalBuf%MAFdata(MAF)

       DO MIF = 0, CurMAFdata%last_MIF

!! Initialize to "U"ndefined:

          CurMAFdata%ChanType(MIF)%FB = undefined

!! Rule #1: Data quality information:

          IF (.NOT. CurMAFdata%SciMIF(MIF)%CRC_good) THEN

             !! Discard all

             CurMAFdata%ChanType(MIF)%FB = discard

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
       
          SwMirPos = CurMAFdata%SciMIF(MIF)%SwMirPos

          !! Possibly use sequence from configuration:

          IF (THz_seq_use == match) THEN           ! Type Match
             IF (SwMirPos /= THz_seq(MIF)) THEN  ! "D"iscard if not a match
                SwMirPos = discard
             ENDIF
          ELSE IF (THz_seq_use == override) THEN   ! Type Override
             IF (SwMirPos /= discard) THEN       ! Override if not a discard
                SwMirPos = THz_seq(MIF)
             ENDIF
          ENDIF

          WHERE (CurMAFdata%ChanType(MIF)%FB == undefined)
             CurMAFdata%ChanType(MIF)%FB = SwMirPos
             !! NOTE: The GHz module could be using FB 12 (FB(:,6)!
          END WHERE

! Reset Sw Pos

          CurMAFdata%SciMIF(MIF)%SwMirPos = SwMirPos

       ENDDO

    ENDDO

  END SUBROUTINE QualifyEachMAF

!=============================================================================
  SUBROUTINE SortAndQualifyTHz
!=============================================================================

    CALL FillCalData

    CALL QualifyEachMAF

  END SUBROUTINE SortAndQualifyTHz

!=============================================================================
END MODULE SortQualifyTHz
!=============================================================================

! $Log$
! Revision 2.3  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.2  2003/02/10 20:32:45  perun
! Increase calbuf size.
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
