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

