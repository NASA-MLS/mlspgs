! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
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

  ! This enumerated type defines the `modules' in the MLS instrument

  INTEGER, PARAMETER :: MLSInstrumentModule_Invalid=0
  INTEGER, PARAMETER :: MLSInstrumentModule_GHz=1
  INTEGER, PARAMETER :: MLSInstrumentModule_THz=2
  INTEGER, PARAMETER :: MLSInstrumentNoModules=2
  CHARACTER(LEN=3), PARAMETER, DIMENSION(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNames = (/ &
       & "GHz", &
       & "THz" /)
  CHARACTER (LEN=3), PARAMETER, DIMENSION(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNamesUC=(/"GHZ","THZ"/)

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
! Revision 2.1  2001/01/25 19:52:51  perun
! Moved L1BInfo to MLSL1Common
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.14  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
