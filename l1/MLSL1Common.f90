! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
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
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
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
    INTEGER :: EngId           ! The HDF ID (handle) for the L1BENG file
    INTEGER :: DiagId          ! The ID (non-HDF) (handle) for the L1BDIAG file
    CHARACTER (LEN=FileNameLen) :: OAFileName  ! L1BOA file name
    CHARACTER (LEN=FileNameLen) :: RADDFileName, RADFFilename
    CHARACTER (LEN=FileNameLen) :: EngFileName   ! L1BENG file name
    CHARACTER (LEN=FileNameLen) :: DiagFileName  ! L1BDIAG file name
  END TYPE L1BFileInfo_T
  TYPE (L1BFileInfo_T) :: L1BFileInfo

  ! Spectrometer info

  INTEGER, PARAMETER :: FBnum = 19
  INTEGER, PARAMETER :: FBchans = 25
  INTEGER, PARAMETER :: MBnum = 5
  INTEGER, PARAMETER :: MBchans = 11
  INTEGER, PARAMETER :: DACSnum = 4
  INTEGER, PARAMETER :: DACSchans = 129
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

  ! Real types for all channels:

  TYPE Chan_R_T
     REAL :: FB(FBchans,FBnum)         ! standard filter banks
     REAL :: MB(MBchans,MBnum)         ! mid-band filter banks
     REAL :: WF(WFchans,WFnum)         ! wide filters
     REAL :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE Chan_R_T

  TYPE Chan_R8_T
     REAL(r8) :: FB(FBchans,FBnum)         ! standard filter banks
     REAL(r8) :: MB(MBchans,MBnum)         ! mid-band filter banks
     REAL(r8) :: WF(WFchans,WFnum)         ! wide filters
     REAL(r8) :: DACS(DACSchans,DACSnum)   ! DACS filters
  END TYPE Chan_R8_T

!! Factor to convert 23 bit encoder angle data (shifted by 1 bit) into degrees:

  REAL, PARAMETER :: DEG24 = 360.0 * 2.0**(-24)

!! Switching mirror position ranges (L, S, T1, T2)
!! Note: these will be initialized as part of the user input parameters.

  TYPE SwMir_Range_T
     CHARACTER(len=1) :: pos
     REAL :: low_angle, high_angle
  END TYPE SwMir_Range_T

!! Switching mirror ranges ( currently +/- 0.1 degrees!)

  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range_A(4) = (/ &  ! "A" side
       SwMir_Range_T ("L", 149.4, 149.6), &
       SwMir_Range_T ("S", 329.4, 329.6), &
       SwMir_Range_T ("T", 239.4, 239.6), &          ! Primary target
       SwMir_Range_T ("t", 59.4, 59.6) /)            ! Secondary target
  TYPE (SwMir_Range_T), TARGET :: GHz_SwMir_Range_B(4) = (/ &  ! "B" side
       SwMir_Range_T ("L", 329.499, 329.699), &
       SwMir_Range_T ("S", 149.499, 149.699), &
       SwMir_Range_T ("T", 59.499, 59.699), &        ! Primary target
       SwMir_Range_T ("t", 239.499, 239.699) /)      ! Secondary target
  TYPE (SwMir_Range_T), TARGET :: THz_SwMir_Range(3) = (/ &
       SwMir_Range_T ("L", 89.0, 91.0), &
       SwMir_Range_T ("S", 0.0, 2.0), &
       SwMir_Range_T ("T", 269.0, 271.0) /)
  TYPE (SwMir_Range_T), TARGET :: Discard_SwMir_Range(1) = (/ &
       SwMir_Range_T ("D", 0.0, 360.0) /)            ! no good pointing

!! Band switches and associated filter bank numbers

  INTEGER :: BandSwitch(5) = (/ 0, 0, 0, 0, 0 /)
  INTEGER, PARAMETER :: SwitchBank(5) = (/ 1, 3, 8, 12, 15 /)

!! Bandwidths for all channels (currently read in during initialization)

  TYPE (Chan_R_T) :: BandWidth

!! Useful channel defaults:

  TYPE (Chan_R_T), TARGET :: deflt_gain, deflt_zero

  REAL :: tau = 1.0 / 6.0  ! default data rate (secs)

!! Useful constants:

  REAL, PARAMETER :: absZero_C = -273.16
  REAL(r8), PARAMETER :: planck = 6.62606891d-34
  REAL(r8), PARAMETER :: boltz = 1.380658d-23

!! First LO frequency (Hz) for each radiometer (except for R1, which is center
!!  freq)

  REAL, PARAMETER :: LO1(5) = (/ 118.7530e9, 191.9000e9, 239.6600e9, &
       642.8700e9, 2.5227816e12 /)

  ! --------------------------------------------------------------------------
  

!=============================================================================
END MODULE MLSL1Common
!=============================================================================

! $Log$
! Revision 2.3  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
