! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInitBOA ! Opens L1 BOA file
!=============================================================================

  USE SDPToolkit, ONLY: PGS_PC_GetReference, PGS_S_SUCCESS, PGSd_EOS_AURA
  USE MLSCommon, ONLY: r8, TAI93_Range_T
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: OpenAndInitBOA, anTextPCF, anTextCF

  CHARACTER (LEN=1), POINTER :: anTextPCF(:), anTextCF(:)

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OpenAndInitBOA
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, GetL1Config
    USE InitPCFs, ONLY: L1PCF, GetPCFParameters
    USE MLSPCF1, ONLY: mlspcf_pcf_start, mlspcf_l1cf_start
    USE PCFHdr, ONLY: CreatePCFAnnotation
    USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo
    USE Orbit, ONLY: Orbit_init, altG, altT, ascTAI, dscTAI, numOrb, &
         orbIncline, orbitNumber, scanRate, scanRateT
    USE OutputL1B, ONLY: OutputL1BOA_create
    USE FOV, ONLY: InitFOVconsts

    CHARACTER (LEN=132) :: PhysicalFilename

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

!! No need to expand time

    L1Config%Expanded_TAI = L1Config%Input_TAI

!! Init MAF info:

    MAFinfo%MIFsPerMAF = L1Config%Calib%MIFsPerMAF
    MAFinfo%startTAI = L1Config%Expanded_TAI%startTime
    MAFinfo%MIF_dur =L1Config%Calib%MIF_duration

!! Make sure this is a simulation!

    L1Config%Globals%SimOA = .TRUE.

!! Init orbit data:

    numOrb = 1

    CALL Orbit_init (procRange, L1PCF%startUTC, altG, altT, ascTAI, dscTAI, &
         numOrb, orbIncline, orbitNumber, scanRate, scanRateT)

!! Open L1B output file

    CALL OpenL1BOAFile

!! Initialize FOV constants

    CALL InitFOVconsts

!! Define the SD structures in the output BOA file

    CALL OutputL1BOA_create (L1BFileInfo, .FALSE.)

  END SUBROUTINE OpenAndInitBOA

!=============================================================================
  SUBROUTINE OpenL1BOAFile
!=============================================================================

    USE MLSPCF1, ONLY: mlspcf_l1b_oa_start
    USE MLSL1Common, ONLY: L1BFileInfo, HDFversion
    USE MLSFiles, ONLY: MLS_openFile
    use MLSHDF5, only: mls_h5open
    USE MLSL1Config, ONLY: L1Config
!    USE H5LIB

    CHARACTER (LEN=132) :: PhysicalFilename
    INTEGER :: error, returnStatus, sd_id, version

!! Open the HDF 5 Fortran Interface

    error = 0
    CALL mls_h5open (error)
    IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         "Fortran HDF 5 API error on opening.")

    ! Initialize IDs

    L1BFileInfo%RADGID = 0
    L1BFileInfo%RADDID = 0
    L1BFileInfo%RADTID = 0
    L1BFileInfo%OAID = 0
    L1BFileInfo%ENGID = 0
    L1BFileInfo%DIAGID = 0

    ! Open L1BOA File

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_l1b_oa_start, version, &
     PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN

       ! Open the HDF file and initialize the SD interface

       CALL MLS_openFile (PhysicalFilename, 'create', sd_id, hdfVersion)
       CALL MLSMessage (MLSMSG_Info, &
            & ModuleName, "Opened L1BOA file: "//PhysicalFilename)
       L1BFileInfo%OAID = sd_id
       L1BFileInfo%OAFileName = PhysicalFilename

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BOA file entry")

    ENDIF

  END SUBROUTINE OpenL1BOAFile

!=============================================================================
END MODULE OpenInitBOA
!=============================================================================

! $Log$
! Revision 2.4  2004/05/06 21:59:23  pwagner
! Uses mls_h5open/close
!
! Revision 2.3  2004/01/09 20:02:57  perun
! Update BOA to HDF 5
!
! Revision 2.2  2004/01/09 18:10:23  perun
! Version 1.4 commit
!
!
