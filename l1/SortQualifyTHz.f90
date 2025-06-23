! Copyright 2006, by the California Institute of Technology. ALL
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
  USE MLSFillValues, ONLY: isNaN

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
    USE dates_module, ONLY: tai93s2hid

    LOGICAL :: more_data

    INTEGER :: sci_MAFno, CalMAFs, CalMAFno
    REAL :: MAF_dur, MIF_dur
    logical, parameter :: verbose = .false. ! Thz is off more days than on

    MIF_dur = L1Config%Calib%MIF_duration
    MAF_dur = MIF_dur * L1Config%Calib%MIFsPerMAF   !Nominal duration of MAF
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

       PRINT *, "THz SCI/ENG MAF, time: ", sci_MAFno, EngMAF%MAFno, &
            & tai93s2hid( EngMAF%sectai, leapsec=.true. )
       CalMAFno = CalMAFno + 1
       CurMAFdata => CalBuf%MAFdata(CalMAFno)
       CurMAFdata%SciMIF = THzSciMAF
       CurMAFdata%EMAF = EngMAF
       CurMAFdata%last_MIF = MAXVAL (THzSciMAF%MIFno)
       CurMAFdata%BandSwitch(4:5) =  &
            CurMAFdata%SciMIF(CurMAFdata%last_MIF)%BandSwitch
       CurMAFdata%CalTgtTemp = &
            GetIndexedAvg (EngMAF%eng%value, CalTgtIndx%THzAmb) - absZero_C
       CurMAFdata%LimbCalTgtTemp = &
            GetIndexedAvg (EngMAF%eng%value, CalTgtIndx%THzLimb) - absZero_C
       CurMAFdata%BO_stat = THz_BO_stat(:,CalMAFno)

       more_data =  THzSciMAF(0)%secTAI <= L1Config%Input_TAI%endTime .AND. &
            (CalMAFno < SIZE (OA_counterMAF))
       if ( isNaN( CurMAFdata%CalTgtTemp) .and. verbose ) then
         print *, 'CalTgtTemp is a NaN'
         print *, 'absZero_C ', absZero_C
         CurMAFdata%CalTgtTemp = &
            GetIndexedAvg ( EngMAF%eng%value, CalTgtIndx%THzAmb, debug=.true.) - absZero_C
       endif
       if ( isNaN( CurMAFdata%LimbCalTgtTemp) .and. verbose ) then
         print *, 'LimbCalTgtTemp is a NaN'
         CurMAFdata%LimbCalTgtTemp = &
            GetIndexedAvg ( EngMAF%eng%value, CalTgtIndx%THzLimb, debug=.true.) - absZero_C
       endif
       ! print *, minval(CalTgtIndx%THzAmb), maxval(CalTgtIndx%THzAmb)
       ! print *, minval(CalTgtIndx%THzAmb), maxval(CalTgtIndx%THzLimb)
       ! print *, size(EngMAF%eng%value)
       ! if ( any(isNaN(EngMAF%eng%value)) ) print *, 'At least one value is a NaN'
       IF (.NOT. more_data) EXIT    !! Nothing more to do

    ENDDO

    IF (SIZE (OA_counterMAF) > 1) &
         CalMAFno = MIN (SIZE (OA_counterMAF), CalMAFno)
    CalBuf%MAFs = CalMAFno

  END SUBROUTINE FillCalData

!=============================================================================
  SUBROUTINE QualifyEachMAF
!=============================================================================

    USE BrightObjects_m, ONLY:Test_BO_stat, BO_Match
    USE Constants, ONLY: Rad2Deg
    USE MLSL1Common, ONLY: L1BFileInfo, SC_YPR, THz_GeodAlt
    USE MLSL1Config, ONLY: THz_seq, THz_seq_use
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

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
    REAL, PARAMETER :: SpaceAltRange(2) = (/ 120.0e03, 600.0e03 /)
    REAL, PARAMETER :: good_yaw = 1.0   ! good YAW value in degrees
    REAL, PARAMETER :: bad_bias = 9.0   ! special "bad" bias value

    PRINT *, 'qualifying all the MAFs'

    DO MAF = 1, CalBuf%MAFs

       CurMAFdata => CalBuf%MAFdata(MAF)

!! Initialize MIF Rad Precision signs:

       CurMAFdata%MIFprecSign = 1.0

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
             !IF (SwMirPos /= discard) THEN       ! Override if not a discard
                SwMirPos = THz_seq(MIF)
             !ENDIF
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


!! Discard if "S"pace altitude out of "good" alt range (pitch maneuver):

          ELSE IF (SwMirPos == 'S') THEN
             IF (THz_GeodAlt(MIF, MAF) < SpaceAltRange(1) .OR. &
                  THz_GeodAlt(MIF, MAF) > SpaceAltRange(2)) SwMirPos = discard
          ENDIF

!! Mark bias if out of "good" YAW position:

          IF (ABS(SC_YPR(1,MIF,MAF) * Rad2Deg) >= good_yaw) THEN
             CurMAFdata%SciMIF(MIF)%LLO_Bias = bad_bias
          ENDIF

! Reset Sw Pos

          CurMAFdata%SciMIF(MIF)%SwMirPos = SwMirPos

       ENDDO

!! Make sure enough calibration data is in MAF:

       IF ((COUNT (CurMAFdata%SciMIF%SwMirPos == 'T') == 0) .OR. &
            (COUNT (CurMAFdata%SciMIF%SwMirPos == 'S') == 0)) &
            CurMAFdata%SciMIF%SwMirPos = discard

!! Check for bright objects in Limb FOV and mark "S"pace views as "D"iscards

       CALL Test_BO_stat (CurMAFdata%BO_stat)

       DO n = 1, BO_Match%Num
          msg = TRIM(BO_Match%Name(n)) // ' in Limb View'
          stat = PGS_TD_TAItoUTC (CurMAFdata%SciMIF(0)%secTAI, asciiUTC)
          CALL MLSMessage (MLSMSG_Warning, ModuleName, &
               TRIM(msg)//' at '//asciiUTC)
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
          WHERE (BO_Match%InFOV(:,n) .AND. CurMAFdata%SciMIF%SwMirPos == 'S')
             CurMAFdata%SciMIF%SwMirPos = discard
          ENDWHERE
          WHERE (BO_Match%InFOV(:,n))
             CurMAFdata%MIFprecSign = BO_Match%PrecScale(n)  ! sign of precision
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
! Revision 2.18  2024/10/09 22:37:38  pwagner
! Stop excess printing on dates THz module is off
!
! Revision 2.17  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.16.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.16  2015/01/13 18:45:19  pwagner
! Prints time with MAF nums; warns if Cal Temps are NaNs
!
! Revision 2.15  2009/05/13 20:33:05  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.14  2006/08/02 18:58:45  perun
! Change maximum SpaceAltRange value to 600 km from 250 km
!
! Revision 2.13  2006/06/14 13:49:29  perun
! Protect reading beyond OA data size
!
! Revision 2.12  2006/03/24 15:18:51  perun
! Add sorting based on "good" altitude range and YAW positions
!
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
