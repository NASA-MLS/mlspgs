! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens input L0 files and output L1 files
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, PGSd_EOS_PM_GIRD, &
       PGSd_IO_Gen_RSeqFrm, PGS_IO_Gen_openF, PGS_IO_Gen_closeF
  USE MLSCommon, ONLY: r8, TAI93_Range_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

  IMPLICIT NONE

  INTEGER :: eng_unit = 102
  INTEGER :: eng_recno = 1   ! start with first record

  CHARACTER (LEN=1), POINTER :: anText(:)

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OpenAndInitialize
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, GetL1Config
    USE InitPCFs, ONLY: L1PCF, GetPCFParameters
    USE MLSPCF1, ONLY: mlspcf_engtbl_start, mlspcf_bwtbl_start, &
         mlspcf_nomen_start, mlspcf_pcf_start
    USE PCFHdr, ONLY: CreatePCFAnnotation
    USE L0_sci_tbls, ONLY: InitSciPointers
    USE MLSL1Common, ONLY: L1BFileInfo
    USE Orbit, ONLY: Orbit_init, altG, altT, ascTAI, dscTAI, numOrb, &
         orbIncline, orbitNumber, scanRate, scanRateT
    USE Calibration, ONLY: InitCalibWindow
    USE EngTbls, ONLY: Load_Eng_tbls
    USE MLSL1Rad, ONLY: InitRad
    USE MLSSignalNomenclature, ONLY: ReadSignalsDatabase
    USE OutputL1B, ONLY: OutputL1B_create

    CHARACTER (LEN=80) :: eng_file, sci_file, eng_tbl_file, bw_tbl_file

    CHARACTER (LEN=132) :: PhysicalFilename

    INTEGER :: bw_tbl_unit
    INTEGER :: eng_tbl_unit
    INTEGER :: nomen_unit
    INTEGER :: errno, ios, lenarg, lenval
    INTEGER :: returnStatus, version

    INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI

    TYPE (TAI93_Range_T) :: procRange

!! Get user parameters from the PCF file

    CALL GetPCFParameters

!! Get annotation from PCF file

    CALL CreatePCFAnnotation (mlspcf_pcf_start, anText)

!! Get the Level 1 configuration from the L1CF file

    CALL GetL1Config

!! TAI Processing range

    returnStatus = PGS_TD_UTCtoTAI (L1PCF%startUTC, procRange%startTime)
    returnStatus = PGS_TD_UTCtoTAI (L1PCF%endUTC, procRange%endTime)

    L1Config%Input_TAI = procRange

!! Will determine the expanded time later!!!

    L1Config%Expanded_TAI = L1Config%Input_TAI

!! TEST expansion time!!!

    L1Config%Expanded_TAI%startTime = L1Config%Expanded_TAI%startTime - 120.0
    L1Config%Expanded_TAI%endTime = L1Config%Expanded_TAI%endTime + 120.0

!! Init orbit data

    CALL Orbit_init (procRange, L1PCF%startUTC, altG, altT, ascTAI, dscTAI, &
         numOrb, orbIncline, orbitNumber, scanRate, scanRateT)

!! Open and initialize eng table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_engtbl_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_engtbl_start, PGSd_IO_Gen_RSeqFrm, &
         0, eng_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open engineering table file" // PhysicalFilename)
    END IF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened engineering table file" // PhysicalFilename)

    CALL Load_Eng_tbls (eng_tbl_unit, ios)

    returnStatus = PGS_IO_Gen_CloseF (eng_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading engineering table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close engineering table file")
    END IF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed engineering table file")

!! Open and initialize bandwidths table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_bwtbl_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_bwtbl_start, PGSd_IO_Gen_RSeqFrm, &
         0, bw_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open bandwidths table file" // PhysicalFilename)
    END IF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened bandwidths table file" // PhysicalFilename)

    CALL GetBandwidths (bw_tbl_unit, ios)

    returnStatus = PGS_IO_Gen_CloseF (bw_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading bandwidth table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close bandwidths table file")
    END IF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed bandwidths table file")

!! Open and Read nomenclature file

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_nomen_start, PGSd_IO_Gen_RSeqFrm, &
         0, nomen_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open nomenclature file")
    END IF

    CALL ReadSignalsDatabase (nomen_unit)

    returnStatus = PGS_IO_Gen_closeF (nomen_unit)  ! close unit
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close nomenclature file")
    END IF

!! Open L0 files:

    CALL OpenL0Files

!! Open L1B output files

    CALL OpenL1BFiles

!! Initialize the science data pointers:

    CALL InitSciPointers

!! Initialize the calibration window:

    CALL InitCalibWindow

!! Initialize Radiances:

    CALL InitRad
             
!! Define the SD structures in the output files
          
    CALL OutputL1B_create (L1BFileInfo)

  END SUBROUTINE OpenAndInitialize

!=============================================================================
  SUBROUTINE OpenL0File (pcf, unit, filename, filetype, TAI_range)
!=============================================================================

    INTEGER :: pcf, unit
    CHARACTER(LEN=*) :: filename, filetype
    TYPE (TAI93_Range_T), OPTIONAL :: TAI_range

    INTEGER :: returnStatus, version
    CHARACTER(LEN=480) :: msg
    TYPE (TAI93_Range_T) :: TAI

    ! Externals

    INTEGER, EXTERNAL :: PGS_IO_L0_Open

    version = 1
    returnStatus = PGS_PC_getReference (pcf, version, filename)

    IF (returnStatus /= PGS_S_SUCCESS) THEN

       write (msg, "('Error opening ', a, ' file for PCF: ', i5)") filetype, pcf
       CALL MLSMessage (MLSMSG_Error, ModuleName, msg)

    END IF

    returnStatus = PGS_IO_L0_Open (pcf, PGSd_EOS_PM_GIRD, unit, &
         TAI%startTime, TAI%endTime)
    IF (PRESENT (TAI_range)) THEN
       TAI_range%startTime = TAI%startTime
       TAI_range%endTime = TAI%endTime
    ENDIF

    IF (returnStatus /= PGS_S_SUCCESS) THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Error opening L0 "//filetype//" file: "//filename)

    END IF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened L0 "//filetype//" file: "//filename)

  END SUBROUTINE OpenL0File

!=============================================================================
  SUBROUTINE OpenL0Files
!=============================================================================

    USE MLSPCF1, ONLY: mlspcf_APID1732_start, mlspcf_APID1734_start, &
         mlspcf_APID1736_start, mlspcf_APID1738_start, mlspcf_APID1740_start, &
         mlspcf_APID1742_start, mlspcf_APID1744_start, mlspcf_APID1746_start
    USE MLSL1Common, ONLY: L0FileInfo
    USE MLSL1Config, ONLY: L1Config

    INTEGER :: i, returnStatus

    TYPE (TAI93_Range_T) :: TAI_range
    INTEGER, EXTERNAL :: PGS_IO_L0_SetStart
    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    ! Init starting PCFs:

    L0FileInfo%eng_pcf(1) = mlspcf_APID1732_start
    L0FileInfo%eng_pcf(2) = mlspcf_APID1734_start
    L0FileInfo%eng_pcf(3) = mlspcf_APID1736_start
    L0FileInfo%eng_pcf(4) = mlspcf_APID1738_start
    L0FileInfo%eng_pcf(5) = mlspcf_APID1740_start
    L0FileInfo%eng_pcf(6) = mlspcf_APID1742_start

    L0FileInfo%sci_pcf(1) = mlspcf_APID1744_start
    L0FileInfo%sci_pcf(2) = mlspcf_APID1746_start

    ! Open L0 Science Files:

    DO i = 1, 2

       DO
          CALL OpenL0File (L0FileInfo%sci_pcf(i), L0FileInfo%sci_unit(i), &
               L0FileInfo%SciFilename(i), "Science", TAI_range)

          !! Check against start time

          IF (L1Config%Expanded_TAI%startTime > TAI_range%endTime) THEN

             returnStatus = PGS_IO_L0_Close (L0FileInfo%sci_unit(i))           
             CALL MLSMessage (MLSMSG_Info, ModuleName, &
                  & 'Closed L0 Science file: '//L0FileInfo%SciFileName(i))
             L0FileInfo%sci_pcf(i) = L0FileInfo%sci_pcf(i) + 1

          ELSE IF (L1Config%Expanded_TAI%startTime < TAI_range%startTime) THEN

             L1Config%Expanded_TAI%startTime = TAI_range%startTime
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & "Adjusting start time to beginning of first data file")

          ELSE IF (L1Config%Expanded_TAI%startTime < TAI_range%startTime) THEN

             CALL MLSMessage (MLSMSG_Error, ModuleName, &
                  "Requested tart time is before first data file start time")

          ELSE
             EXIT
          ENDIF

       ENDDO

       returnStatus = PGS_IO_L0_SetStart (L0FileInfo%sci_unit(i), &
            L1Config%Expanded_TAI%startTime)

       IF (returnStatus /= PGS_S_SUCCESS) THEN

          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & "Error positioning L0 file: "//L0FileInfo%SciFilename(i))

       END IF

    ENDDO

    ! Open L0 Engineering  Files:

    DO i = 1, 6

       DO
          CALL OpenL0File (L0FileInfo%eng_pcf(i), L0FileInfo%eng_unit(i), &
               L0FileInfo%EngFilename(i), "Engineering", TAI_range)

          !! Check against start time

          IF (L1Config%Expanded_TAI%startTime > TAI_range%endTime) THEN

             returnStatus = PGS_IO_L0_Close (L0FileInfo%eng_unit(i))
             CALL MLSMessage (MLSMSG_Info, ModuleName, &
                  & 'Closed L0 Engineering file: '//L0FileInfo%EngFileName(i))
             L0FileInfo%eng_pcf(i) = L0FileInfo%eng_pcf(i) + 1

          ELSE IF (L1Config%Expanded_TAI%startTime < TAI_range%startTime) THEN

             L1Config%Expanded_TAI%startTime = TAI_range%startTime
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & "Adjusting start time to beginning of first data file")

          ELSE IF (L1Config%Expanded_TAI%startTime < TAI_range%startTime) THEN

             CALL MLSMessage (MLSMSG_Error, ModuleName, &
                  "Requested tart time is before first data file start time")

          ELSE
             EXIT
          ENDIF

       ENDDO

       returnStatus = PGS_IO_L0_SetStart (L0FileInfo%eng_unit(i), &
            L1Config%Expanded_TAI%startTime)

       IF (returnStatus /= PGS_S_SUCCESS) THEN

          CALL MLSMessage (MLSMSG_Error, ModuleName, &
               & "Error positioning L0 file: "//L0FileInfo%EngFilename(i))

       END IF

    ENDDO

  END SUBROUTINE OpenL0Files

!=============================================================================
  SUBROUTINE OpenL1BFiles
!=============================================================================

    USE MLSPCF1, ONLY: mlspcf_l1b_radf_start, mlspcf_l1b_radd_start, &
         mlspcf_l1b_oa_start, mlspcf_l1b_eng_start
    USE MLSL1Common, ONLY: L1BFileInfo
    USE Hdf, ONLY: DFACC_CREATE, sfstart

    CHARACTER (LEN=132) :: PhysicalFilename
    INTEGER :: returnStatus, sd_id, version

    ! Open L1BRADF File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_radd_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN

       ! Open the HDF file and initialize the SD interface

       sd_id = sfstart (PhysicalFilename, DFACC_CREATE)
       IF (sd_id == -1) THEN
          CALL MLSMessage (MLSMSG_Error, &
               & ModuleName, "Error creating L1BRADD file: "//PhysicalFilename)
       ELSE
          CALL MLSMessage (MLSMSG_Info, &
               & ModuleName, "Opened L1BRADD file: "//PhysicalFilename)
          L1BFileInfo%RADDID = sd_id
          L1BFileInfo%RADDFileName = PhysicalFilename
       END IF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BRADD file entry")

    END IF

    ! Open L1BRADF File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_radf_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN

       ! Open the HDF file and initialize the SD interface

       sd_id = sfstart (PhysicalFilename, DFACC_CREATE)
       IF (sd_id == -1) THEN
          CALL MLSMessage (MLSMSG_Error, &
               & ModuleName, "Error creating L1BRADF file: "//PhysicalFilename)
       ELSE
          CALL MLSMessage (MLSMSG_Info, &
               & ModuleName, "Opened L1BRADF file: "//PhysicalFilename)
          L1BFileInfo%RADFID = sd_id
          L1BFileInfo%RADFFileName = PhysicalFilename
       END IF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BRADF file entry")

    END IF

    ! Open L1BOA File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_oa_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN

       ! Open the HDF file and initialize the SD interface

       sd_id = sfstart (PhysicalFilename, DFACC_CREATE)
       IF (sd_id == -1) THEN
          CALL MLSMessage (MLSMSG_Error, &
               & ModuleName, "Error creating L1BOA file: "//PhysicalFilename)
       ELSE
          CALL MLSMessage (MLSMSG_Info, &
               & ModuleName, "Opened L1BOA file: "//PhysicalFilename)
          L1BFileInfo%OAID = sd_id
          L1BFileInfo%OAFileName = PhysicalFilename
       END IF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BOA file entry")

    END IF

  END SUBROUTINE OpenL1BFiles

!=============================================================================
  SUBROUTINE GetBandwidths (bw_tbl_unit, ios)
!=============================================================================

    USE MLSL1Common, ONLY: FBnum, MBnum, WFnum, BandWidth

    INTEGER :: bw_tbl_unit, ios
    CHARACTER (len=40) :: line
    INTEGER ::  i
    REAL, PARAMETER :: scale = 1.0e6    ! MHz scale factor

    READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line      ! Title line
    DO i = 1, FBnum
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Blank line
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Band
       READ (bw_tbl_unit, *, IOSTAT=ios, ERR=999) BandWidth%FB(:,i)    ! Banks
    ENDDO
    BandWidth%FB = BandWidth%FB * scale

    DO i = 1, MBnum
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Blank line
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Band
       READ (bw_tbl_unit, *, IOSTAT=ios, ERR=999) BandWidth%MB(:,i)    ! Banks
    ENDDO
    BandWidth%MB = BandWidth%MB * scale

    DO i = 1, WFnum
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Blank line
       READ (bw_tbl_unit, '(a)', IOSTAT=ios, ERR=999) line   ! Band
       READ (bw_tbl_unit, *, IOSTAT=ios, ERR=999) BandWidth%WF(:,i)    ! Banks
    ENDDO
    BandWidth%WF = BandWidth%WF * scale

999 CONTINUE

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, &
            & ModuleName, "Error reading bandwidths file!")
    ENDIF

  END SUBROUTINE GetBandwidths

!=============================================================================
END MODULE OpenInit
!=============================================================================

! $Log$
! Revision 2.2  2001/03/12 15:59:13  perun
! Add reading PCF file
!
! Revision 2.1  2001/02/23 20:54:11  perun
! Version 0.5 commit
!
