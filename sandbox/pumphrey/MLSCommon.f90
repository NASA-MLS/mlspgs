! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! This file is Hugh Pumphrey's private version, changed slightly so it 
! compiles in the F language

!=============================================================================
module MLSCommon                ! Common definitions for the MLS software
!=============================================================================

  implicit none
  !private

  !private :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  character(len=256),parameter,private :: Id = &
       "$Id$"
  character (len=*), parameter,private:: ModuleName=&
       "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module contains simple definitions which are common to all the MLS PGS
  ! f90 software.

  ! Firstly, these are standard numerical types, copied from HCP
  ! (again with my change in case, sorry Hugh!)

  integer,parameter,public:: i1=selected_int_kind(2)
  integer,parameter,public:: i2=selected_int_kind(4)
  integer,parameter,public:: i4=selected_int_kind(7)
  integer,parameter,public:: r4=selected_real_kind(5)
  integer,parameter,public:: r8=selected_real_kind(13)

  ! Now we have the lengths for various strings

  integer,parameter,public :: NameLen=32
  integer,parameter,public :: LineLen=132

  ! This enumerated type defines the `modules' in the MLS instrument

  integer, parameter,public :: MLSInstrumentModule_Invalid=0
  integer, parameter,public :: MLSInstrumentModule_GHz=1
  integer, parameter,public :: MLSInstrumentModule_THz=2
  integer, parameter,public :: MLSInstrumentNoModules=2
  character (len=3), parameter,public, dimension(MLSInstrumentNoModules) :: &
        MLSInstrumentModuleNames= (/ &
        "GHz", &
        "THz"/)
  character (len=3), parameter,public, dimension(MLSInstrumentNoModules) :: &
        MLSInstrumentModuleNamesUC=(/"GHZ","THZ"/)

  ! --------------------------------------------------------------------------
  
  ! The next datatype describes the information on the L1B data files in use

  type,public:: L1BInfo_T
     integer :: L1BOAId     ! The HDF ID (handle) for the L1BOA file
     integer, dimension(:), pointer :: L1BRADIDs ! ID(s) for the L1BRAD file(s)
  end type L1BInfo_T

  ! --------------------------------------------------------------------------

  ! This datatype defines the `chunks' into which the input dataset is split

  type,public:: MLSChunk_T
     integer :: firstMAFIndex   ! Index of first MAF in the chunk
     integer :: lastMAFIndex    ! Index of last MAF in the chunk
     integer :: noMAFsLowerOverlap ! Number of MAFs in the lower overlap region
     integer :: noMAFsUpperOverlap ! Number of MAFs in the upper overlap region
     integer :: accumulatedMAFs ! Number of non overlapped MAFs before this.
  end type MLSChunk_T

  ! --------------------------------------------------------------------------

!=============================================================================
end module MLSCommon
!=============================================================================

!
! $Log$
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
!
