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
!---------------------------------------------------------------------------
  integer, parameter :: RK = kind(0.0d0)
contains
  subroutine D_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    include 'pcspl.f9h'
  end subroutine D_PCSPL
end module D_PCSPL_M
! $Log$
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
