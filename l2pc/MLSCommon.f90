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
  INTEGER,PARAMETER:: i1=SELECTED_INT_KIND(2)
  INTEGER,PARAMETER:: i2=SELECTED_INT_KIND(4)
  INTEGER,PARAMETER:: i4=SELECTED_INT_KIND(7)
  INTEGER,PARAMETER:: r4=SELECTED_REAL_KIND(5)
  INTEGER,PARAMETER:: r8=SELECTED_REAL_KIND(13)
  ! Now we have the lengths for various strings
  INTEGER,PARAMETER :: NameLen=32
  INTEGER,PARAMETER :: LineLen=132
  INTEGER,PARAMETER :: FileNameLen=132
  ! This enumerated type defines the `modules' in the MLS instrument
  INTEGER, PARAMETER :: MLSInstrumentModule_Invalid=0
  INTEGER, PARAMETER :: MLSInstrumentModule_GHz=1
  INTEGER, PARAMETER :: MLSInstrumentModule_THz=2
  INTEGER, PARAMETER :: MLSInstrumentNoModules=2
  CHARACTER (LEN=3), PARAMETER, DIMENSION(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNames= (/ &
       & "GHz", &
       & "THz"/)
  CHARACTER (LEN=3), PARAMETER, DIMENSION(MLSInstrumentNoModules) :: &
       & MLSInstrumentModuleNamesUC=(/"GHZ","THZ"/)
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
! Revision 1.12  2000/02/08 19:54:57  nakamura
! Moved TAI93_Range_T here from MLSL1Common.
!
! Revision 1.11  2000/02/02 14:59:53  perun
! Add filenames to L1BInfo Type
!
! Revision 1.10  1999/12/21 18:43:56  livesey
! Reinstated accumulatedMAFs
!
! Revision 1.9  1999/12/21 00:15:57  livesey
! Removed accumulatedMAFs from MLSChunk_T.  This functionality is now served by
! accumulatedProfiles in L2GPData_T
!
! Revision 1.8  1999/12/18 01:04:11  livesey
! Changed name to MLSInstrumentNoModules
!
! Revision 1.7  1999/12/17 21:33:45  livesey
! Still have a problem with MLSInstrumentModuleNames here.  It just doesn't
! seem to behave correctly.
!
! Revision 1.6  1999/12/17 00:59:31  livesey
! Nightly checkin
!
! Revision 1.5  1999/12/16 17:46:31  livesey
! Made the InstrumentModuleNames fixed dimensions.
!
! Revision 1.4  1999/12/16 01:27:28  livesey
! Added instrument module stuff
!
! Revision 1.3  1999/12/15 23:59:54  livesey
! Moved L1BInfo_T and MLSChunk_T in from MLSL2Common
!
! Revision 1.2  1999/12/14 01:01:14  livesey
! Whoops, didn't compile!
!
! Revision 1.1  1999/12/14 00:45:37  livesey
! First version, basically just the old `kinds' module.
!
