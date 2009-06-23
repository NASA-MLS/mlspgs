! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Hermite_Coeffs_m

! Given the function values and derivatives at the end points of intervals,
! compute coefficients for Hermite cubic interpolation within those intervals.

  implicit NONE
  private

  public :: Hermite_Coeffs, Hermite_Coeffs_D, Hermite_Coeffs_S

  interface Hermite_Coeffs
    module procedure Hermite_Coeffs_D, Hermite_Coeffs_S
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------  Hermite_Coeffs_D  -----
  subroutine Hermite_Coeffs_D ( Funcs, Derivs, Coeffs )
  ! Given Funcs(0:n) and Derivs(0:n), compute Coeffs(0:3,1:n) of monomial
  ! basis polonomials that interpolate the given functions and derivatives
  ! in each of the N intervals.  Each interval is scaled on 0:1, so to
  ! compute f(xx) for x(i-1) <= xx <= x(i) use R = (xx-x(i-1))/(x(i)-x(i-1),
  ! then F = coeffs(0,i) + r*(coeffs(1,i) + r*(coeffs(2,i) + r*coeffs(3,i))).
    integer, parameter :: RK = kind(0.0d0)
    include "Hermite_Coeffs.f9h"
  end subroutine Hermite_Coeffs_D

  ! -------------------------------------------  Hermite_Coeffs_S  -----
  subroutine Hermite_Coeffs_S ( Funcs, Derivs, Coeffs )
  ! Given Funcs(0:n) and Derivs(0:n), compute Coeffs(0:3,1:n) of monomial
  ! basis polonomials that interpolate the given functions and derivatives
  ! in each of the N intervals.  Each interval is scaled on 0:1, so to
  ! compute f(xx) for x(i-1) <= xx <= x(i) use R = (xx-x(i-1))/(x(i)-x(i-1),
  ! then F = coeffs(0,i) + r*(coeffs(1,i) + r*(coeffs(2,i) + r*coeffs(3,i))).
    integer, parameter :: RK = kind(0.0e0)
    include "Hermite_Coeffs.f9h"
  end subroutine Hermite_Coeffs_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Hermite_Coeffs_m

! $Log$
! Revision 2.2  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2006/07/19 22:31:14  vsnyder
! Initial commit
!
