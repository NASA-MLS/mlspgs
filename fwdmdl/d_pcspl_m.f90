! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_PCSPL_M
  implicit NONE
  private
  public :: D_PCSPL, PCSPL
  interface PCSPL; module procedure D_PCSPL; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  integer, parameter :: RK = kind(0.0d0)
contains
  subroutine D_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    include 'pcspl.f9h'
  end subroutine D_PCSPL
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module D_PCSPL_M
! $Log$
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.5  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
