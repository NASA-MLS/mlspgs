! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module D_SOLVE_QUAD_M

  implicit NONE

  private
  public :: D_SOLVE_QUAD, SOLVE_QUAD

  interface SOLVE_QUAD; module procedure D_SOLVE_QUAD; end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains
! Solve for the roots of a quadratic equation: a*x*x+b*x+c=0

  subroutine D_SOLVE_QUAD ( A, B, C, X1, X2 )
    integer, parameter :: RK = kind(0.0d0)
    include 'solve_quad.f9h'
  end subroutine D_SOLVE_QUAD

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module D_SOLVE_QUAD_M

! $Log$
! Revision 2.1  2002/10/08 17:08:02  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
