! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens input L0 files and output L1 files
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, PGSd_EOS_AURA, &
       PGSd_IO_Gen_RSeqFrm, PGS_IO_Gen_openF, PGS_IO_Gen_closeF, &
       PGSd_IO_Gen_WSeqUnf, PGSd_IO_Gen_USeqUnf
  USE MLSCommon, ONLY: r8, TAI93_Range_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

  IMPLICIT NONE

  CHARACTER (LEN=1), POINTER :: anTextPCF(:), anTextCF(:)

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
         mlspcf_nomen_start, mlspcf_pcf_start, mlspcf_l1cf_start, &
         mlspcf_defltgains_start, mlspcf_defltzeros_start
    USE PCFHdr, ONLY: CreatePCFAnnotation
    USE L0_sci_tbls, ONLY: InitSciPointers
    USE MLSL1Common, ONLY: L1BFileInfo, deflt_gain, deflt_zero
    USE Orbit, ONLY: Orbit_init, altG, altT, ascTAI, dscTAI, numOrb, &
         orbIncline, orbitNumber, scanRate, scanRateT
    USE Calibration, ONLY: InitCalibWindow
    USE EngTbls, ONLY: Load_Eng_tbls
    USE MLSL1Rad, ONLY: InitRad
    USE MLSSignalNomenclature, ONLY: ReadSignalsDatabase
    USE OutputL1B, ONLY: OutputL1B_create
    USE DACSUtils, ONLY: InitDACS_FFT

    CHARACTER (LEN=132) :: PhysicalFilename

    INTEGER :: bw_tbl_unit
    INTEGER :: eng_tbl_unit
    INTEGER :: gains_tbl_unit
    INTEGER :: nomen_unit
    INTEGER :: zeros_tbl_unit
    INTEGER :: errno, ios, lenarg, lenval
    INTEGER :: returnStatus, version

    INTEGER, EXTERNAL :: PGS_TD_UTCtoTAI

    TYPE (TAI93_Range_T) :: procRange

!! Get user parameters from the PCF file

    CALL GetPCFParameters

!! Get annotation from PCF and CF files

    CALL CreatePCFAnnotation (mlspcf_pcf_start, anTextPCF)

    CALL CreatePCFAnnotation (mlspcf_l1cf_start, anTextCF)

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

!! Will determine the expanded time later!!!

    L1Config%Expanded_TAI = L1Config%Input_TAI

!! TEST expansion time!!!

    L1Config%Expanded_TAI%startTime = L1Config%Expanded_TAI%startTime - 120.0
    L1Config%Expanded_TAI%endTime = L1Config%Expanded_TAI%endTime + 120.0

!! Init orbit data

    numOrb = 1

    IF (L1Config%Globals%ProduceL1BOA) THEN

       CALL Orbit_init (procRange, L1PCF%startUTC, altG, altT, ascTAI, dscTAI, &
            numOrb, orbIncline, orbitNumber, scanRate, scanRateT)

    ENDIF

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

!! Open and initialize bandwidths table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_bwtbl_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_bwtbl_start, PGSd_IO_Gen_RSeqFrm, &
         0, bw_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open bandwidths table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened bandwidths table file: " // PhysicalFilename)

    CALL GetBandwidths (bw_tbl_unit, ios)

    returnStatus = PGS_IO_Gen_CloseF (bw_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading bandwidth table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close bandwidths table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed bandwidths table file")

!! Open and initialize default gains table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_defltgains_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_defltgains_start, &
         PGSd_IO_Gen_RSeqFrm, 0, gains_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open default gains table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened default gains table file: " // PhysicalFilename)

    CALL GetChanDefaults (gains_tbl_unit, deflt_gain, ios)

    returnStatus = PGS_IO_Gen_CloseF (gains_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading default gains table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close default gains table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed default gains table file")

!! Open and initialize default zeros table:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_defltzeros_start, version, &
          & PhysicalFilename)

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_defltzeros_start, &
         PGSd_IO_Gen_RSeqFrm, 0, zeros_tbl_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open default zeros table file: " // PhysicalFilename)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened default zeros table file: " // PhysicalFilename)

    CALL GetChanDefaults (zeros_tbl_unit, deflt_zero, ios)

    returnStatus = PGS_IO_Gen_CloseF (zeros_tbl_unit)

    IF (ios /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         & "Error reading default zeros table file")

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close default zeros table file")
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed default zeros table file")

!! Open and Read nomenclature file

    version = 1
    returnStatus = PGS_IO_Gen_openF (mlspcf_nomen_start, PGSd_IO_Gen_RSeqFrm, &
         0, nomen_unit, version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open nomenclature file")
    ENDIF

    CALL ReadSignalsDatabase (nomen_unit)

    returnStatus = PGS_IO_Gen_closeF (nomen_unit)  ! close unit
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close nomenclature file")
    ENDIF

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

!! Initialize DACS FFT

    CALL InitDACS_FFT
          
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

    ENDIF

    returnStatus = PGS_IO_L0_Open (pcf, PGSd_EOS_AURA, unit, &
         TAI%startTime, TAI%endTime)
    IF (PRESENT (TAI_range)) THEN
       TAI_range%startTime = TAI%startTime
       TAI_range%endTime = TAI%endTime
    ENDIF

    IF (returnStatus /= PGS_S_SUCCESS) THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Error opening L0 "//filetype//" file: "//filename)

    ENDIF

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

       ENDIF

    ENDDO

    ! Open L0 Engineering Files:

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

       ENDIF

    ENDDO

  END SUBROUTINE OpenL0Files

!=============================================================================
  SUBROUTINE OpenL1BFiles
!=============================================================================

    USE MLSPCF1, ONLY: mlspcf_l1b_radf_start, mlspcf_l1b_radd_start, &
         mlspcf_l1b_oa_start, mlspcf_l1b_eng_start, mlspcf_l1b_diag_start
    USE MLSL1Common, ONLY: L1BFileInfo
    USE EngTbls, ONLY: Eng_tbl, maxtlm
    USE Hdf, ONLY: DFACC_CREATE, sfstart

    CHARACTER (LEN=132) :: PhysicalFilename
    INTEGER :: returnStatus, sd_id, version

    ! Open L1BRADD File

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
       ENDIF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BRADD file entry")

    ENDIF

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
       ENDIF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BRADF file entry")

    ENDIF

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
       ENDIF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BOA file entry")

    ENDIF

    ! Open L1BENG File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_eng_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN
       version = 1

!! Will replace with HDF in a future version!

       returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_eng_start, &
            PGSd_IO_Gen_WSeqUnf, 0, sd_id, version)
       IF (returnstatus /= PGS_S_SUCCESS) THEN
          returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_eng_start, &
               PGSd_IO_Gen_USeqUnf, 0, sd_id, version)
          IF (returnstatus /= PGS_S_SUCCESS) THEN
             CALL MLSMessage (MLSMSG_Error, ModuleName, &
                  & "Could not open L1B engineering file: " // PhysicalFilename)
          ENDIF
       ENDIF
    
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & "Opened L1B engineering file: " // PhysicalFilename)
       L1BFileInfo%EngId = sd_id
       L1BFileInfo%EngFileName = PhysicalFilename

! write header info

       WRITE (sd_id) maxtlm
       WRITE (sd_id) eng_tbl%mnemonic

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BENG file entry")

    ENDIF

    ! Open L1BDIAG File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_diag_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN
       version = 1

!! This file is NOT an HDF file.

       returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_diag_start, &
            PGSd_IO_Gen_WSeqUnf, 0, sd_id, version)
       IF (returnstatus /= PGS_S_SUCCESS) THEN
          returnStatus = PGS_IO_Gen_openF (mlspcf_l1b_diag_start, &
               PGSd_IO_Gen_USeqUnf, 0, sd_id, version)
          IF (returnstatus /= PGS_S_SUCCESS) THEN
             CALL MLSMessage (MLSMSG_Error, ModuleName, &
                  & "Could not open L1B diagnostics file: " // PhysicalFilename)
          ENDIF
       ENDIF
    
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & "Opened L1B diagnostics file: " // PhysicalFilename)
       L1BFileInfo%DiagId = sd_id
       L1BFileInfo%DiagFileName = PhysicalFilename

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BDIAG file entry")

    ENDIF

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

  SUBROUTINE GetChanDefaults (unit, deflt_dat, stat)

    USE MLSL1Common, ONLY: Chan_R_T, WFNum

    INTEGER :: unit, stat
    TYPE (Chan_R_T), TARGET :: deflt_dat

    CHARACTER(LEN=2) :: chan_type
    CHARACTER(LEN=80) :: line
    INTEGER :: ios, i
    REAL, POINTER, DIMENSION(:,:) :: chan_dat

    stat = 0

! Read comments until start of data

    DO
       READ (unit, '(A)') line
       IF (line(1:6) == "#DATA") EXIT
    ENDDO

! read data

    i = 1
    chan_type = "FB"
    chan_dat => deflt_dat%fb
    DO

       DO
          READ (unit, '(A)', IOSTAT=ios) line
          IF (line(1:1) == "#" .or. ios /= 0) EXIT
       ENDDO

       IF (ios /= 0) EXIT

       READ (unit, *) chan_dat(:,i)

       i = i + 1

       SELECT CASE (chan_type)

       CASE ("FB")
          IF (I > SIZE (deflt_dat%fb(1,:))) THEN
             i = 1
             chan_dat => deflt_dat%mb
             chan_type = "MB"
          ENDIF
       CASE ("MB")
          IF (I > SIZE (deflt_dat%mb(1,:))) THEN
             i = 1
             chan_dat => deflt_dat%wf
             chan_type = "WF"
          ENDIF
       END SELECT

    ENDDO

    IF (chan_type /= "WF" .AND. i <= WFNum) stat = -1  ! Error!!!

  END SUBROUTINE GetChanDefaults

!=============================================================================
END MODULE OpenInit
!=============================================================================

! $Log$
! Revision 2.5  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.4  2001/03/22 16:46:02  perun
! Check CF and PCF outputVersions
!
! Revision 2.3  2001/03/12 19:36:00  perun
! Read and save CF file as annotation
!
! Revision 2.1  2001/02/23 20:54:11  perun
! Version 0.5 commit
!
