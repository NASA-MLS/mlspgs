! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module contains simple definitions which are common to all the MLS PGS
  ! f90 software.

  ! Firstly, these are standard numerical types, copied from HCP
  ! (again with my change in case, sorry Hugh!)

  INTEGER, PARAMETER:: i1=SELECTED_INT_KIND(2)
  INTEGER, PARAMETER:: i2=SELECTED_INT_KIND(4)
  INTEGER, PARAMETER:: i4=SELECTED_INT_KIND(7)
  INTEGER, PARAMETER:: r4=SELECTED_REAL_KIND(5)
  INTEGER, PARAMETER:: r8=SELECTED_REAL_KIND(13)

  ! Now we have the lengths for various strings

  INTEGER, PARAMETER :: NameLen=32
  INTEGER, PARAMETER :: LineLen=132
  INTEGER, PARAMETER :: FileNameLen=132

  ! --------------------------------------------------------------------------
  
  ! The next datatype describes the information on the L1B data files in use

  TYPE L1BInfo_T
    INTEGER :: L1BOAId     ! The HDF ID (handle) for the L1BOA file
    INTEGER, DIMENSION(:), POINTER :: L1BRADIds ! Id(s) for the L1BRAD file(s)
    CHARACTER (LEN=FileNameLen) :: L1BOAFileName  ! L1BOA file name
    CHARACTER (LEN=FileNameLen), DIMENSION(:), POINTER :: L1BRADFileNames
  END TYPE L1BInfo_T

  ! --------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  TYPE MLSChunk_T
    INTEGER :: firstMAFIndex   ! Index of first MAF in the chunk
    INTEGER :: lastMAFIndex    ! Index of last MAF in the chunk
    INTEGER :: noMAFsLowerOverlap ! Number of MAFs in the lower overlap region
    INTEGER :: noMAFsUpperOverlap ! Number of MAFs in the upper overlap region
    INTEGER :: accumulatedMAFs ! Number of non overlapped MAFs before this.
  END TYPE MLSChunk_T

  ! --------------------------------------------------------------------------

  ! The TAI93 time range

  TYPE TAI93_Range_T
    REAL(r8) :: startTime ! TAI93 format
    REAL(r8) :: endTime   ! TAI93 format
  END TYPE TAI93_Range_T
  ! --------------------------------------------------------------------------

!=============================================================================
END MODULE MLSCommon
!=============================================================================

!
! $Log$
! Revision 2.3  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.2  2001/01/26 23:46:35  pwagner
! Restored L1BInfo from l1/MLSL1Common back to lib/MLSCommon
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
