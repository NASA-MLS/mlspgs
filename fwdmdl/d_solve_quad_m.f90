! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_SOLVE_QUAD_M
  implicit NONE
  private
  public :: D_SOLVE_QUAD, SOLVE_QUAD
  interface SOLVE_QUAD; module procedure D_SOLVE_QUAD; end interface
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
  integer, parameter :: RK = kind(0.0d0)
contains
! Solve for the roots of a quadratic equation: a*x*x+b*x+c=0
!
  subroutine D_SOLVE_QUAD ( A, B, C, X1, X2 )
    include 'solve_quad.f9h'
  end subroutine D_SOLVE_QUAD
end module D_SOLVE_QUAD_M
! $Log$
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:04  vsnyder
! Initial conversion to Fortran 90
!
