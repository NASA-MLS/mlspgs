! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OpenAndInitializeLog
!=============================================================================

    USE OpenInit, ONLY: OpenL0Files
    USE MLSL1Config, ONLY: L1Config, GetL1Config
    USE InitPCFs, ONLY: L1PCF, GetPCFParameters
    USE MLSPCF1, ONLY: mlspcf_engtbl_start, mlspcf_l1b_log_start, &
         mlspcf_pcf_start, mlspcf_l1cf_start, mlspcf_sciMAF_start, &
         mlspcf_engMAF_start, mlspcf_MAF_data_start, mlspcf_l1b_eng_start
    USE L0_sci_tbls, ONLY: InitSciPointers
    USE EngTbls, ONLY: Load_Eng_tbls, Eng_tbl, maxtlm
    USE MLSL1Common, ONLY: L1BFileInfo

    CHARACTER (LEN=132) :: PhysicalFilename

    INTEGER :: eng_tbl_unit, ios, log_unit
    INTEGER :: returnStatus, version

    INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI, PGS_IO_Gen_Track_LUN

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

    IF (L1Config%Globals%OutputVersionString /= L1PCF%OutputVersion) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "CF and PCF OutputVersions do not match!")
    ENDIF

!! TAI Processing range

    returnStatus = PGS_TD_UTCtoTAI (L1PCF%startUTC, procRange%startTime)
    returnStatus = PGS_TD_UTCtoTAI (L1PCF%endUTC, procRange%endTime)

    L1Config%Input_TAI = procRange

    L1Config%Expanded_TAI = L1Config%Input_TAI

!! Expand time range!!!

    L1Config%Expanded_TAI%startTime = L1Config%Expanded_TAI%startTime - 120.0
    L1Config%Expanded_TAI%endTime = L1Config%Expanded_TAI%endTime + 120.0

!! Open and initialize eng table:

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

!! Open L1BENG File

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

!! Open Eng/Sci MAF files:

    WRITE (PhysicalFilename, "(I3.3)") mlspcf_engMAF_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_engMAF_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_engMAF_start: " // &
            PhysicalFilename)
    ENDIF

    returnStatus = PGS_IO_Gen_Track_LUN (L1BFileInfo%EngMAF_unit, 0)

    OPEN (unit=L1BFileInfo%EngMAF_unit, file=PhysicalFilename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open engMAF file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened engMAF file: " // PhysicalFilename)

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

    WRITE (PhysicalFilename, "(I3.3)") mlspcf_MAF_data_start
    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_MAF_data_start, version, &
          & PhysicalFilename)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find file entry for pcf_MAF_data_start: " // &
            PhysicalFilename)
    ENDIF

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

  END SUBROUTINE OpenAndInitializeLog

!=============================================================================
END MODULE OpenInitLog
!=============================================================================

! $Log$
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
