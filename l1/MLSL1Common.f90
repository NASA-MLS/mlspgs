! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Common              ! Common data types for the MLSL1 program
!=============================================================================

  USE MLSCommon, ONLY: r4, r8, FileNameLen

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$$"
  !---------------------------------------------------------------------------

  ! This module contains data types that are common to the MLSL1 program.

  !---------------------------------------------------------------------------
  
  ! The science and engineering Level 0 file info

  TYPE L0FileInfo_T
     INTEGER :: sci_unit(2), sci_pcf(2)
     INTEGER :: eng_unit(6), eng_pcf(6)
     CHARACTER (LEN=FileNameLen) :: SciFileName(2)
     CHARACTER (LEN=FileNameLen) :: EngFileName(6)
  END TYPE L0FileInfo_T

  TYPE (L0FileInfo_T) :: L0FileInfo

  ! Information on the L1B data files in use

  TYPE L1BFileInfo_T
    INTEGER :: OAId            ! The HDF ID (handle) for the L1BOA file
    INTEGER :: RADDId, RADFId  ! Ids for the L1BRAD files
    CHARACTER (LEN=FileNameLen) :: OAFileName  ! L1BOA file name
    CHARACTER (LEN=FileNameLen) :: RADDFileName, RADFFilename
  END TYPE L1BFileInfo_T
  TYPE (L1BFileInfo_T) :: L1BFileInfo

  ! Spectrometer info

  INTEGER, PARAMETER :: FBnum = 19
  INTEGER, PARAMETER :: FBchans = 25
  INTEGER, PARAMETER :: MBnum = 5
  INTEGER, PARAMETER :: MBchans = 11
  INTEGER, PARAMETER :: DACnum = 4
  INTEGER, PARAMETER :: DACchans = 129
  INTEGER, PARAMETER :: WFnum = 3
  INTEGER, PARAMETER :: WFchans = 4

  ! Maximum number of MIFs per MAF

  INTEGER, PARAMETER :: MaxMIFs = 150

  ! MAF information data type (start time, integration time, etc.)

  TYPE MAFinfo_T
     REAL(r8) :: startTAI    ! TAI93 format
     REAL(r4) :: integTime   ! integration time
     INTEGER :: MIFsPerMAF   ! Number of MIFs per MAF
  END TYPE MAFinfo_T

!! Factor to convert 23 bit encoder angle data (shifted by 1 bit) into degrees:

  REAL, PARAMETER :: DEG24 = 360.0 * 2.0**(-24)

!! Switching mirror position ranges (L, S, T1, T2)
!! Note: these will be initialized as part of the user input parameters.

  TYPE SwMir_Range_T
     CHARACTER(len=2) :: pos
     REAL :: low_angle, high_angle
  END TYPE SwMir_Range_T

!! Will probably move this and use assignments

  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range(4) = (/ &
       SwMir_Range_T ("L", 89.0, 91.0), &
       SwMir_Range_T ("S", 0.0, 2.0), &
       SwMir_Range_T ("T1", 269.0, 271.0), &
       SwMir_Range_T ("T2", 179.0, 181.0) /)
  TYPE (SwMir_Range_T), TARGET :: THz_SwMir_Range(3) = (/ &
       SwMir_Range_T ("L", 89.0, 91.0), &
       SwMir_Range_T ("S", 0.0, 2.0), &
       SwMir_Range_T ("T1", 269.0, 271.0) /)

!! Bandwidths for all channels (currently read in during initialization)

  TYPE BandWidth_T
     REAL :: FB(FBchans,FBnum)
     REAL :: MB(MBchans,MBnum)
     REAL :: WF(WFchans,WFnum)
  END TYPE BandWidth_T

  TYPE (BandWidth_T) :: BandWidth

  REAL :: tau = 1.0 / 6.0  ! default data rate (secs)

!! Useful constants:

  REAL, PARAMETER :: absZero_C = -273.16

  ! --------------------------------------------------------------------------
  

!=============================================================================
END MODULE MLSL1Common
!=============================================================================

! $Log$
