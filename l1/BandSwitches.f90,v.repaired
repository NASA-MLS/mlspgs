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
MODULE BandSwitches   ! Deal with Band switches in L0 data
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, PGS_IO_Gen_closeF
  USE MLSCommon, ONLY: r8
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
       MLSMSG_Info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: bsw_unit, BandSwFmt, GetBandSwitches, OpenBandSwitchesFile

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  INTEGER :: bsw_unit
  CHARACTER (LEN=*), PARAMETER :: BandSwFmt = "(A27, 5x, 5I3)"

  INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI

CONTAINS

!=============================================================================
  SUBROUTINE OpenBandSwitchesFile (mode)
!=============================================================================

    USE MLSPCF1, ONLY: mlspcf_bandsw_start
    USE MLSL1Common, ONLY: FileNameLen


    CHARACTER (len=*) :: mode

    CHARACTER (LEN=FileNameLen) :: PhysicalFilename

    INTEGER :: ios, returnStatus, version

    INTEGER, EXTERNAL :: PGS_IO_Gen_Track_LUN

!! Open band switches table:

    WRITE (PhysicalFilename, "(I5.5)") mlspcf_bandsw_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_bandsw_start, version, &
     PhysicalFilename)

    IF (returnStatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find Band Switches file entry")
    ENDIF

    returnStatus = PGS_IO_Gen_Track_LUN (bsw_unit, 0)

! Open with file positioned at EOF:

    OPEN (unit=bsw_unit, file=PhysicalFilename, ACTION=mode, status="OLD", &
         FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=ios, POSITION="APPEND")

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open Band Switches file: " // PhysicalFilename)
    ENDIF
    
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened Band Switches file: " // PhysicalFilename)

  END SUBROUTINE OpenBandSwitchesFile

!=============================================================================
  SUBROUTINE GetBandSwitches (TAI_in, BandNo) ! Get latest bands from switches
!=============================================================================

    REAL(r8), INTENT (IN) :: TAI_in
    INTEGER, INTENT (OUT) :: BandNo(5)

    CHARACTER(len=27) :: asciiUTC
    INTEGER :: n
    REAL(r8) :: TAI
    LOGICAL :: first_read = .TRUE.

!! Open Band Switches file:

    CALL OpenBandSwitchesFile ("READ")   ! Read mode

!! Read data until first TAI in file is less than requested TAI:

    BACKSPACE bsw_unit   ! Beginning of last record
    DO
       READ (bsw_unit, fmt=BandSwFmt) asciiUTC, bandno
       n = PGS_TD_UTCtoTAI (asciiUTC, TAI)
       IF (TAI <= TAI_in) THEN
          IF (first_read .AND. ((TAI_in - TAI) > 86400.0)) THEN
             PRINT *, 'Switch positions not guaranteed for this time request!'
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & 'Time requested is more than 24 hours beyond ' // &
                  'guaranteed switch positions')
          ENDIF
          EXIT
       ENDIF
       BACKSPACE bsw_unit    ! Record just read
       BACKSPACE bsw_unit    ! Previous record to read
       first_read = .FALSE.  ! First record already read
    ENDDO

!! Close Band Switches file:

    n = PGS_IO_Gen_CloseF (bsw_unit)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed Band Switches file")

  END SUBROUTINE GetBandSwitches

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here

END MODULE BandSwitches
!=============================================================================
! $Log$
! Revision 2.3  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.2.6.1  2015/10/09 10:21:37  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.2  2006/04/05 18:11:12  perun
! Remove unused variables
!
! Revision 2.1  2006/03/24 15:04:27  perun
! Initial release for opening and getting band switches from Switch Network database file
!
