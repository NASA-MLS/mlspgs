! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL1Common              ! Common data types for the MLSL1 program
!=============================================================================

  USE MLSCommon

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module simply contains data types that are common to the MLSL1
  ! program.

  !---------------------------------------------------------------------------
  
  ! The science and engineering Level 0 file info
  ! Note: versions 0.5 and above will be dimensioned 3!

  TYPE L0Info_T
     INTEGER :: sci
     INTEGER :: eng
     CHARACTER (LEN=FileNameLen) :: SciFileName
     CHARACTER (LEN=FileNameLen) :: EngFileName
  END TYPE L0Info_T

  ! Level 0 Science Packet (sized for 1 MIF)

  CHARACTER (LEN=2048) :: l0_sci

  ! Spectrometer info

  INTEGER, PARAMETER :: FBnum = 19
  INTEGER, PARAMETER :: FBchans = 25
  INTEGER, PARAMETER :: MBnum = 5
  INTEGER, PARAMETER :: MBchans = 11
  INTEGER, PARAMETER :: DACnum = 4
  INTEGER, PARAMETER :: DACchans = 129
  INTEGER, PARAMETER :: WFnum = 3
  INTEGER, PARAMETER :: WFchans = 4

  ! Level 0 spectrometer science data type (sized for 1 MIF)

  TYPE L0Sci_T
     INTEGER FB(FBchans,FBnum)    ! standard filter banks
     INTEGER MB(MBchans,MBnum)    ! mid-band filter banks
     INTEGER DAC(DACchans,DACnum) ! digital autocorrelation spectrometers
     INTEGER WF(WFchans,WFnum)    ! wide filters
  END TYPE L0Sci_T

  ! Information on the L1B data files in use

  TYPE L1BInfo_T
    INTEGER :: L1BOAId      ! The HDF ID (handle) for the L1BOA file
    INTEGER :: L1BRADIds(2) ! Ids for the L1BRAD files (1=FB, 2=DAC)
    CHARACTER (LEN=FileNameLen) :: L1BOAFileName  ! L1BOA file name
    CHARACTER (LEN=FileNameLen) :: L1BRADFileNames(2)
  END TYPE L1BInfo_T

  ! Maximum number of MIFs per MAF

  INTEGER, PARAMETER :: MaxMIFs = 150

  ! MAF information data type (start time, integration time, etc.)

  TYPE MAFinfo_T
     REAL(r8) :: startTAI    ! TAI93 format
     REAL(r4) :: integTime   ! integration time
     INTEGER :: MIFsPerMAF   ! Number of MIFs per MAF
  END TYPE MAFinfo_T

  ! -------------------------------------------------------------------------- 

!=============================================================================
END MODULE MLSL1Common
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:55:14  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/02/08 20:31:07  perun
! Initial version
!
