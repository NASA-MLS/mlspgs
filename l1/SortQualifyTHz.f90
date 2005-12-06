! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE SortQualifyTHz ! Sort and qualify the L0 data for the THz module
!=============================================================================

  USE L0_sci_tbls, ONLY: THzSciMAF
  USE THzCalibration, ONLY: CalBuf, MAFdata_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SortAndQualifyTHz

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()

CONTAINS

!=============================================================================
  SUBROUTINE FillCalData
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Common, ONLY: absZero_C, OA_counterMAF
    USE EngTbls, ONLY: EngMAF, CalTgtIndx 
    USE EngUtils, ONLY: NextEngMAF
    USE SciUtils, ONLY: NextSciMAF
    USE MLSL1Utils, ONLY : GetIndexedAvg
    USE BrightObjects_m, ONLY: THz_BO_stat

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

       DO    ! Read sci data until after start time

          CALL NextSciMAF (more_data)

          IF (.NOT. more_data) EXIT    !! Nothing more to do

          IF (THzSciMAF(2)%secTAI >= L1Config%Expanded_TAI%startTime) EXIT

       ENDDO
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

PRINT *, "SCI/ENG MAF: ", sci_MAFno, EngMAF%MAFno
       CalMAFno = CalMAFno + 1
       CurMAFdata => CalBuf%MAFdata(CalMAFno)
       CurMAFdata%SciMIF = THzSciMAF
       CurMAFdata%EMAF = EngMAF
       CurMAFdata%last_MIF = MIFsPerMAF - 1
       CurMAFdata%BandSwitch(4:5) =  &
            CurMAFdata%SciMIF(CurMAFdata%last_MIF)%BandSwitch
       CurMAFdata%CalTgtTemp = &
            GetIndexedAvg (EngMAF%eng%value, CalTgtIndx%THzAmb) - absZero_C
       CurMAFdata%LimbCalTgtTemp = &
            GetIndexedAvg (EngMAF%eng%value, CalTgtIndx%THzLimb) - absZero_C
       CurMAFdata%BO_stat = THz_BO_stat(:,CalMAFno)

       more_data =  THzSciMAF(0)%secTAI <= L1Config%Input_TAI%endTime

    ENDDO

    IF (SIZE (OA_counterMAF) > 1) &
         CalMAFno = MIN (SIZE (OA_counterMAF), CalMAFno)
    CalBuf%MAFs = CalMAFno

  END SUBROUTINE FillCalData

!=============================================================================
  SUBROUTINE QualifyEachMAF
!=============================================================================

    USE MLSL1Config, ONLY: THz_seq, THz_seq_use
    USE MLSL1Common, ONLY: L1BFileInfo
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC
    USE BrightObjects_m, ONLY:Test_BO_stat, BO_Match

    INTEGER :: MAF, MIF, n, stat
    CHARACTER(len=80) :: msg
    CHARACTER(len=27) :: asciiUTC
    REAL :: encoder(2)

    CHARACTER(len=1) :: SwMirPos
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: limb = "L"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"
    REAL, PARAMETER :: BadLimbRange(2) = (/ 6.0, 354.0 /)

    PRINT *, 'qualifying all the MAFs'

    DO MAF = 1, CalBuf%MAFs

       CurMAFdata => CalBuf%MAFdata(MAF)

       DO MIF = 0, CurMAFdata%last_MIF

!! Initialize to "U"ndefined:

          CurMAFdata%ChanType(MIF)%FB = undefined

          !! Save current sw pos:
       
          SwMirPos = CurMAFdata%SciMIF(MIF)%SwMirPos

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

!! Discard all bands (for now)

          IF (ANY (CurMAFdata%SciMIF(MIF)%MaxAtten)) SwMirPos = discard

!! Rule #5: Set Switching Mirror position

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

!! Discard if out of "good" limb angle range:

          IF (SwMirPos == limb) THEN
             encoder =  CurMAFdata%SciMIF(MIF)%TSSM_pos
             IF (ANY (encoder > BadLimbRange(1) .AND. &
                  encoder < BadLimbRange(2))) SwMirPos = discard
          ENDIF

! Reset Sw Pos

          CurMAFdata%SciMIF(MIF)%SwMirPos = SwMirPos

       ENDDO

!! Make sure enough calibration data is in MAF:

       IF ((COUNT (CurMAFdata%SciMIF%SwMirPos == 'T') == 0) .OR. &
            (COUNT (CurMAFdata%SciMIF%SwMirPos == 'S') == 0)) &
            CurMAFdata%SciMIF%SwMirPos = discard

!! Check for bright objects in Limb FOV and mark as "D"iscards

       CALL Test_BO_stat (CurMAFdata%BO_stat)

       DO n = 1, BO_Match%Num
          msg = TRIM(BO_Match%Name(n)) // ' in Limb View'
          stat = PGS_TD_TAItoUTC (CurMAFdata%SciMIF(0)%secTAI, asciiUTC)
          CALL MLSMessage (MLSMSG_Warning, ModuleName, &
               TRIM(msg)//' at '//asciiUTC)
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
          WHERE (BO_Match%InFOV(:,n))
             CurMAFdata%SciMIF%SwMirPos = discard
          ENDWHERE
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
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE SortQualifyTHz
!=============================================================================

! $Log$
! Revision 2.11  2005/12/06 19:29:51  perun
! Removed call to Flag_Bright_Objects and added testing BO_stat
!
! Revision 2.10  2005/10/14 15:55:44  perun
! Restrict maximum size of CalBuf to size of OA_counterMAF
!
! Revision 2.9  2005/10/12 17:12:27  perun
! Check for enough cal views in MAF
!
! Revision 2.8  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2005/05/02 16:07:41  perun
! Read sci data until at or after the requested start time
!
! Revision 2.6  2005/01/28 17:07:14  perun
! Get THz FOV bright object flags
!
! Revision 2.5  2004/11/10 15:36:11  perun
! Check for "Bright Objects" in FOV; check encoder for good range (-/+ 6.0)
!
! Revision 2.4  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
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
