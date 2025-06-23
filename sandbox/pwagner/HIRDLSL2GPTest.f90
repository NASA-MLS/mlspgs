! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=================================
PROGRAM HIRDLSL2GPtest ! reads HIRDLS L2GPData file
!=================================

   use L2GPData, only: Dump, L2GPData_T, ReadL2GPData
   use MLSFiles, only: MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF, &
    & HDFVERSION_4, HDFVERSION_5
   use MLSMessageModule, only: MLSMessageConfig
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the L2GPData subroutines.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! rm -f l2gp.h4 l2gptest.out ; echo "1,4" | LF95.Linux.atlas/test > l2gptest.out

! Variables
   type (L2GPData_T)            :: l2gp
   integer                      :: hdfVersion
   integer                      :: returnStatus, swfid, record_length, swid
   character(len=*), parameter  :: l2gpFilename = &
     & '/users/pwagner/docs/mls/l2gpsamples/HIRDLS2-MOSS3_b035_2000d275.he5'
   !  & '/users/pwagner/docs/mls/l2gpsamples/HIRDLS2-Aura12h_b027_2000d275.he5'
   character(len=*), parameter  :: swathName = 'H2O'
   ! Executable
   MLSMessageConfig%logFileUnit = -1
   MLSMessageConfig%useToolkit = .false.
   call h5open_f(returnStatus)

   hdfVersion = HDFVERSION_5
   call ReadL2GPData ( l2gpFilename, swathName, l2gp, &
    & hdfVersion=hdfVersion, HMOT='H' )
   call dump(l2gp)
   call h5close_f(returnStatus)
!==================
END PROGRAM HIRDLSL2GPtest
!==================

! $Log$
! Revision 1.1  2003/04/17 23:02:43  pwagner
! First commit
!
