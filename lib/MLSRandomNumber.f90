! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSRandomNumber              ! Some random number-generating things
  !=============================================================================

  use MLSCommon, only : R8
  use MLSMessageModule, only: MLSMessage,MLSMSG_Error

!------------------------------------------------------------------
!   Random number routines from MATH77 libraries
! ../l2:69% ls *.f
!  drang.f  ranpk1.f  ranpk2.f  srang.f
!  ../l2:71% cat *.f > stuff.sed
!  ../l2:72% sed 's/^[Cc]/\!/' stuff.sed > sed.out
! plus a small amount of subsequent editing
!------------------------------------------------------------------
  implicit none

  private
  public :: srang, drang

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

!     c o n t e n t s
!     - - - - - - - -

! drang             gauss. distribution: 0 mean, 1 s.d. (double)
! srang             gauss. distribution: 0 mean, 1 s.d. (single)

   ! (Just a "stub" until I clean up the f77 a little more)

contains

      double precision function  DRANG ()
         DRANG=0.0
      end function  DRANG
      double precision function  SRANG ()
         SRANG=0.0
      end function  SRANG

!=============================================================================
end module MLSRandomNumber
!=============================================================================

!
! $Log$
! Revision 2.2  2001/09/24 17:27:07  pwagner
! Fixed blunder
!
! Revision 2.1  2001/09/24 17:22:08  pwagner
! First commit
!
