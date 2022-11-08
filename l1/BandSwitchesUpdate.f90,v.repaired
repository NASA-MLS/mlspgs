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
MODULE BandSwitchesUpdate   ! Update Band switches database using L0 data
!=============================================================================

  USE SDPToolkit, ONLY: PGS_IO_Gen_closeF, PGS_TD_TAItoUTC
  USE MLSCommon, ONLY: TAI93_Range_T, r8
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning, MLSMSG_Info
  USE BandSwitches, ONLY: bsw_unit, BandSwFmt, OpenBandSwitchesFile

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OpenInitBsw, ExamineSciData, CloseFiles

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (TAI93_Range_T) :: TAI_range

  INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI

CONTAINS

!=============================================================================
  SUBROUTINE OpenInitBsw
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, GetL1Config
    USE InitPCFs, ONLY: L1PCF, GetPCFParameters
    USE OpenInit, ONLY: OpenL0Files
    USE L0_sci_tbls, ONLY: InitSciPointers

    INTEGER :: returnStatus

    PRINT *, 'Updating Band Switches database'

!! Get user parameters from the PCF file

    CALL GetPCFParameters

!! Get the Level 1 configuration from the L1CF file

    CALL GetL1Config

!! TAI Processing range

    returnStatus = PGS_TD_UTCtoTAI (L1PCF%startUTC, TAI_range%startTime)
    returnStatus = PGS_TD_UTCtoTAI (L1PCF%endUTC, TAI_range%endTime)
    L1Config%Input_TAI = TAI_range
    L1Config%Expanded_TAI = L1Config%Input_TAI

!! Open L0 files:

    CALL OpenL0Files

!! Initialize the science data pointers:

    CALL InitSciPointers

!! Open Band Switches file

    CALL OpenBandSwitchesFile ("READWRITE")   ! Read/Write mode

  END SUBROUTINE OpenInitBsw

!=============================================================================
  SUBROUTINE CloseFiles   ! Close production Level 0 files
!=============================================================================

    USE MLSL1Common, ONLY: L0FileInfo
  
    INTEGER :: i, returnStatus

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    ! Close L0 Science Files:

    DO i = 1, 2
       
       returnStatus = PGS_IO_L0_Close (L0FileInfo%sci_unit(i))

       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L0 Science file: '//L0FileInfo%SciFileName(i))

    ENDDO

    returnStatus = PGS_IO_Gen_CloseF (bsw_unit)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed Band Switches file")

  END SUBROUTINE CloseFiles

!=============================================================================
  SUBROUTINE ExamineSciData

!=============================================================================

!! Despite the name, all this routine really does is check and/or
!! update the Band Switches file, which stores the state of the GSN
!! switching network.


    USE MLSL1Common, ONLY: BandSwitch
    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: NextSciMAF
    use MLSFillValues, only: Isnan

    INTEGER :: MIF, n
    INTEGER :: bandno(5)
    CHARACTER(len=27) :: asciiUTC
    LOGICAL :: more_data = .TRUE., check_DAY = .TRUE.
    REAL(r8) :: TAI

! Get initial time and band switches from last record in file:

    BACKSPACE bsw_unit
    READ (bsw_unit, fmt=BandSwFmt) asciiUTC, bandno
    n = PGS_TD_UTCtoTAI (asciiUTC, TAI)
    PRINT *, 'BandSwitchesUpdate::ExamineSciData :', asciiUTC, bandno

    BandSwitch = bandno     ! Initial band values
    DO
       DO
          CALL NextSciMAF (more_data)
          if ( isNaN(SciMAF(0)%secTAI) ) then
            print *, 'NaN in secTAI'
            print *, 'more data? ', more_data
            cycle
          endif
          IF (.NOT. more_data .OR. SciMAF(0)%secTAI >= TAI_range%startTime) EXIT
       ENDDO
       IF (more_data) more_data = SciMAF(0)%secTAI <= TAI_range%endTime
       IF (.NOT. more_data) EXIT

       IF (check_DAY) THEN
          IF (SciMAF(0)%secTAI <= TAI) THEN
             PRINT *, 'BandSwitches file already up to date'
             CALL MLSMessage (MLSMSG_Info, ModuleName, &
                  & 'BandSwitches file already up to date')
             RETURN
          ELSE IF (SciMAF(0)%secTAI  > (TAI+86400.0)) THEN
             PRINT *, 'Last entry in file is more than 24 hours before ' // &
                  'current data'
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & 'Last entry in file is more than 24 hours before ' // &
                  'current data')
          ENDIF
          check_DAY = .FALSE.
          BACKSPACE bsw_unit       ! Prepare to replace last line
       ENDIF

       DO MIF = 0, 146
          IF (ALL (SciMAF(MIF)%BandSwitch > 0)) THEN
             IF (ANY (SciMAF(MIF)%BandSwitch /= bandno)) THEN
                n = PGS_TD_TAItoUTC (SciMAF(MIF)%secTAI, asciiUTC)
                bandno = SciMAF(MIF)%BandSwitch
                WRITE (bsw_unit, fmt=BandSwFmt) asciiUTC, bandno
                PRINT *, 'Change in switch network! ,',asciiUTC, bandno
             ENDIF
          ENDIF
       ENDDO
       n = PGS_TD_TAItoUTC (SciMAF(0)%secTAI, asciiUTC)

    ENDDO

! Write last record:

    WRITE (bsw_unit, fmt=BandSwFmt) asciiUTC, bandno

  END SUBROUTINE ExamineSciData

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE BandSwitchesUpdate
!=============================================================================
! $Log$
! Revision 2.5  2022/11/08 23:48:42  pwagner
! Workaround for NAG signaling NaN in secTAI
!
! Revision 2.4  2018/04/09 22:12:19  whdaffer
! Documentation
!
! Revision 2.3  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.2.6.1  2015/10/09 10:21:37  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.2  2006/04/05 18:10:24  perun
! Remove unused variables
!
! Revision 2.1  2006/03/24 15:05:50  perun
! Initial release for reading L0 data and updating band switches in the Switch Network database
!
