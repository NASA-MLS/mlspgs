! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_HUNT_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: D_HUNT, HUNT
  interface HUNT; module procedure D_HUNT; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  integer, private, parameter :: RK = r8
contains
! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.
  subroutine D_HUNT ( ELEMENT, ARRAY, N, JLO, JHI )
    include 'hunt.f9h'
  end subroutine D_HUNT
end module D_HUNT_M
! $Log$
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
