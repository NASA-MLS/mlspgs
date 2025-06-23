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
MODULE OpenInitLog ! Opens input L0 files and output log file
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, &
       PGSd_IO_Gen_RSeqFrm, PGS_IO_Gen_openF, PGS_IO_Gen_closeF, &
       PGSd_IO_Gen_WSeqFrm, PGSd_IO_Gen_USeqFrm
  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OpenAndInitializeLog

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OpenAndInitializeLog
!=============================================================================

    USE OpenInit, ONLY: OpenL0Files, LoadChanDefaults
    USE MLSL1Config, ONLY: L1Config, GetL1Config
    USE InitPCFs, ONLY: L1PCF, GetPCFParameters
    USE MLSPCF1, ONLY: mlspcf_engtbl_start, mlspcf_l1b_log_start, &
         mlspcf_pcf_start, mlspcf_l1cf_start, mlspcf_sciMAF_start, &
         mlspcf_engMAF_start, mlspcf_MAF_data_start, mlspcf_l1b_eng_start, &
         mlspcf_dacsconst_start, mlspcf_defltzeros_start, mlspcf_l1b_radd_start
    USE L0_sci_tbls, ONLY: InitSciPointers
    USE EngTbls, ONLY: Load_Eng_tbls, Eng_tbl, maxtlm
    USE MLSL1Common, ONLY: L1BFileInfo, deflt_zero, BandSwitch, HDFversion,FileNameLen
    USE MLSFiles, ONLY: MLS_openFile, MLS_closeFile
    USE MLSHDF5, ONLY: MLS_h5open
    USE THzUtils, ONLY: LLO_Label
    USE DACSUtils, ONLY: InitDACS_FFT
    USE BandSwitches, ONLY: GetBandSwitches
    USE SDPTOOLKIT, ONLY : PGSd_PC_FILE_PATH_MAX


    CHARACTER (LEN=FileNameLen) :: PhysicalFilename
    CHARACTER (LEN=28) :: asciiUTC_A

    INTEGER :: eng_tbl_unit, error, ios, log_unit, tbl_unit, sd_id
    INTEGER :: returnStatus, version

    REAL :: MAF_dur

    INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI,  PGS_TD_TAItoUTC, &
         PGS_TD_ASCIItime_AtoB, PGS_IO_Gen_Track_LUN

    TYPE (TAI93_Range_T) :: procRange

!! Get PCF and CF filenames:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_pcf_start, version, &
          & PhysicalFilename)
    L1PCF%PCF_filename = physicalfilename

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1cf_start, version, &
          & PhysicalFilename)
    L1PCF%L1CF_filename = physicalfilename

!! Get user parameters from the PCF file

    CALL GetPCFParameters

!! Get the Level 1 configuration from the L1CF file

    CALL GetL1Config

!! Check output versions from CF and PCF
    PRINT *,"About to check CF/L1PCF versions!"
    IF (L1Config%Globals%OutputVersionString /= L1PCF%OutputVersion) THEN
      PRINT *,"CF and PCF OutputVersions do not match!"
      PRINT *,"CF version: "//TRIM(L1Config%Globals%OutputVersionString)
      PRINT *,"PCF version: "//TRIM(L1PCF%OutputVersion)
      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "CF and PCF OutputVersions do not match!")
      CALL MLSMessage (MLSMSG_Info, ModuleName, &
           & "CF version: "//TRIM(L1Config%Globals%OutputVersionString))
      CALL MLSMessage (MLSMSG_Error, ModuleName, &
           & "PCF version: "//TRIM(L1PCF%OutputVersion))

    ENDIF

!! TAI Processing range

    MAF_dur = L1Config%Calib%MIF_duration * L1Config%Calib%MIFsPerMAF

    returnStatus = PGS_TD_UTCtoTAI (L1PCF%startUTC, procRange%startTime)
    returnStatus = PGS_TD_UTCtoTAI (L1PCF%endUTC, procRange%endTime)

    procRange%startTime = procRange%startTime - MAF_dur * (0.5 + &
         L1Config%Calib%MAFexpandNum)
    procRange%endTime = procRange%endTime +  MAF_dur * (0.5 + &
         L1Config%Calib%MAFexpandNum)

    returnStatus = PGS_TD_TAItoUTC (procRange%startTime, asciiUTC_A)
    returnStatus = PGS_TD_ASCIItime_AtoB (asciiUTC_A, L1PCF%startUTC)
    returnStatus = PGS_TD_TAItoUTC (procRange%endTime, asciiUTC_A)
    returnStatus = PGS_TD_ASCIItime_AtoB (asciiUTC_A, L1PCF%endUTC)

    L1Config%Input_TAI = procRange

    L1Config%Expanded_TAI = procRange

!! Expand time range!!!

    L1Config%Expanded_TAI%startTime = L1Config%Expanded_TAI%startTime - 120.0
    L1Config%Expanded_TAI%endTime = L1Config%Expanded_TAI%endTime + 120.0

!! Open and initialize eng table: (<whd>engtlm.tbl, PCF ID: 903.</whd>)

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_engtbl_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_engtbl_start, PGSd_IO_Gen_RSeqFrm, &
         0, eng_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open engineering table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened engineering table file: " // PhysicalFilename)

    CALL Load_Eng_tbls (eng_tbl_unit, ios)

    returnStatus = PGS_IO_Gen_CloseF (eng_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading engineering table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close engineering table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed engineering table file")

!! Open L1BENG File <whd> As far as I can tell, no other software uses this
!! file.</whd>

    WRITE (PhysicalFilename, "(I5.5)") mlspcf_l1b_eng_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_eng_start, version, &
     PhysicalFilename)

    IF (returnStatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BENG file entry")
    ENDIF

    returnStatus = PGS_IO_Gen_Track_LUN (L1BFileInfo%EngId, 0)

    OPEN (unit=L1BFileInfo%EngId, file=PhysicalFilename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open L1B engineering file: " // PhysicalFilename)
    ENDIF
    
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened L1B engineering file: " // PhysicalFilename)

    L1BFileInfo%EngFileName = PhysicalFilename

! write header info

    WRITE (L1BFileInfo%EngId) maxtlm
    WRITE (L1BFileInfo%EngId) eng_tbl%mnemonic

    WRITE (L1BFileInfo%EngId) LLO_Label   ! Laser LO Labels

!! Open L0 files:

    CALL OpenL0Files

!! Initialize the science data pointers:

    CALL InitSciPointers

!! Open log file:

    WRITE (PhysicalFilename, "(I5.5)") mlspcf_l1b_log_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_log_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_log_start: " // &
            PhysicalFilename)
    ENDIF

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_log_start,&
         PGSd_IO_Gen_WSeqFrm, 0, log_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_log_start, &
            PGSd_IO_Gen_USeqFrm, 0, log_unit, version)
       IF (returnstatus /= PGS_S_SUCCESS) THEN
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & "Could not open log file: " // PhysicalFilename)
       ENDIF
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened L1B log file: " // PhysicalFilename)

    L1BFileInfo%LogId = log_unit
    L1BFileInfo%LogFilename = PhysicalFilename

    WRITE (log_unit, *) ''
    WRITE (log_unit, *) '################## Begin MLSL1log ####################'
    WRITE (log_unit, *) ''
    WRITE (log_unit, *) 'PCF filename: '//TRIM (L1PCF%PCF_filename)
    WRITE (log_unit, *) 'L1CF filename: '//TRIM (L1PCF%L1CF_filename)
    WRITE (log_unit, *) 'Input Start/End UTC: '// &
         TRIM (L1PCF%StartUTC)//' to '//TRIM (L1PCF%EndUTC)

!! Init BandSwitches:

    CALL GetBandSwitches (L1Config%Expanded_TAI%startTime, BandSwitch)
    WRITE (log_unit, *) 'Initial BandSwitches: ', BandSwitch

!! Open and initialize default zeros table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_defltzeros_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_defltzeros_start, &
         PGSd_IO_Gen_RSeqFrm, 0, tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open default zeros table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened default zeros table file: " // PhysicalFilename)

    CALL LoadChanDefaults (tbl_unit, deflt_zero, ios)

    returnStatus = PGS_IO_Gen_CloseF (tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading default zeros table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close default zeros table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed default zeros table file")

!! Open and initialize DACS constants table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_dacsconst_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_dacsconst_start, &
         PGSd_IO_Gen_RSeqFrm, 0, tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open DACS constants table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened DACS constants table file: " // PhysicalFilename)

    CALL LoadDACSconsts (tbl_unit, ios)

    returnStatus = PGS_IO_Gen_CloseF (tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading DACS constants table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close DACS constants table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed DACS constants table file")

!! Initialize DACS FFT

    CALL InitDACS_FFT

!! Open Eng/Sci MAF files:

    ! <whd> Temporary file nominally named engMAF_tmp.dat, PCF ID: 921. Eng !
    ! data will be written here, then cobbled together with the matching sci
    ! data later.</whd>
    WRITE (PhysicalFilename, "(I3.3)") mlspcf_engMAF_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_engMAF_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_engMAF_start: " // &
            PhysicalFilename)
    ENDIF

    ! the file read by EngMAF_unit is a temporary file created by
    ! MLSL1log, typically named engMAF_tmp.dat (PCF ID 921). There's
    ! also one for the science data named sciMAF_tmp.dat (PCF ID 920)
    returnStatus = PGS_IO_Gen_Track_LUN (L1BFileInfo%EngMAF_unit, 0)

    OPEN (unit=L1BFileInfo%EngMAF_unit, file=PhysicalFilename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open engMAF file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened engMAF file: " // PhysicalFilename)

    ! <whd> Temporary file nominally named sciMAF_tmp.dat, PCF ID: 920. Science
    ! data will be written here, then cobbled together with the matching eng
    ! data later.</whd>
    WRITE (PhysicalFilename, "(I3.3)") mlspcf_sciMAF_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_sciMAF_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_sciMAF_start: " // &
            PhysicalFilename)
    ENDIF

    returnStatus = PGS_IO_Gen_Track_LUN (L1BFileInfo%SciMAF_unit, 0)

    OPEN (unit=L1BFileInfo%sciMAF_unit, file=PhysicalFilename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open sciMAF file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened sciMAF file: " // PhysicalFilename)

    ! 
    WRITE (PhysicalFilename, "(I3.3)") mlspcf_MAF_data_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_MAF_data_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_MAF_data_start: " // &
            PhysicalFilename)
    ENDIF

    !<whd> quasi-Temporary file nominally named MAF_data_tmp.dat which will
    !contain the colocated eng/sci data as well as the calibration weights
    !flags. PCF ID: 922 </whd>
    returnStatus = PGS_IO_Gen_Track_LUN (L1BFileInfo%MAF_data_unit, 0)

    OPEN (unit=L1BFileInfo%MAF_data_unit, file=PhysicalFilename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", &
         ACTION="READWRITE", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open MAF_data file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened MAF_data file: " // PhysicalFilename)

!! Write some header info for comparisons:

    WRITE (L1BFileInfo%MAF_data_unit) L1PCF%PCF_filename
    WRITE (L1BFileInfo%MAF_data_unit) L1PCF%L1CF_filename

!! Open the HDF 5 Fortran Interface based on CF file

    error = 0
    CALL MLS_h5open (error)
    IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         "Fortran HDF 5 API error on opening.")

    ! Open L1BRADD File

    IF (L1Config%Calib%CalibDACS) THEN

       version = 1
       returnStatus = PGS_PC_getReference (mlspcf_l1b_radd_start, version, &
            PhysicalFilename)

       IF (returnStatus == PGS_S_SUCCESS) THEN

          ! Open the HDF file and initialize the SD interface

          CALL MLS_openFile (PhysicalFilename, 'create', sd_id, hdfVersion)
          CALL MLSMessage (MLSMSG_Info, &
               & ModuleName, "Opened L1BRADD file: "//PhysicalFilename)
          L1BFileInfo%RADDID = sd_id
          L1BFileInfo%RADDFileName = PhysicalFilename

       ELSE
          
          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & "Could not find L1BRADD file entry")

       ENDIF
    ENDIF

  END SUBROUTINE OpenAndInitializeLog

!=============================================================================
  SUBROUTINE LoadDACSconsts (unit, ios)
!=============================================================================

    USE MLSL1Common, ONLY: DACS_const

    INTEGER :: unit, ios
    CHARACTER (len=80) :: line

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:6) == "#DATA") EXIT
    ENDDO

! Read data constants

    DO
       READ (unit, '(A)', IOSTAT=ios, ERR=999) line
       IF (line(1:1) == "#" .OR. ios /= 0) EXIT
    ENDDO

    READ (unit, *, IOSTAT=ios, ERR=999) DACS_const%L

    DO
       READ (unit, '(A)', IOSTAT=ios, ERR=999) line
       IF (line(1:1) == "#" .OR. ios /= 0) EXIT
    ENDDO

    READ (unit, *, IOSTAT=ios, ERR=999) DACS_const%A

999 CONTINUE

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, &
            & ModuleName, "Error reading DACS constants file!")
    ENDIF

   END SUBROUTINE LoadDACSconsts

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE OpenInitLog
!=============================================================================

! $Log$
! Revision 2.11  2018/04/10 16:39:59  whdaffer
! Too many parens!
!
! Revision 2.10  2018/04/09 22:18:34  whdaffer
! Print out disagreement between l1cf versions, as that failure has
! bitten me several times and I've had to spend *way* to much time
! re-figuring it out.
!
! Revision 2.9  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.8.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.8  2007/06/21 21:04:05  perun
! Only create RADD file if DACS calibration is enabled
!
! Revision 2.7  2006/08/02 18:56:35  perun
! Added creation of RADD file in anticipation of writing TPz attribute
!
! Revision 2.6  2006/03/24 15:15:22  perun
! Expand processing times, iniut Band Switches and write startup message to log file
!
! Revision 2.5  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.3  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
