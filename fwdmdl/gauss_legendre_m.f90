module GAUSS_LEGENDRE_M
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: NG, QUAD, QUAD_START

! -----     Private data     -------------------------------------------

  real(r4), save :: ABLEN              ! Half-Length of A-B interval
  real(r4), save :: ABMID              ! Midpoint of A-B interval

 integer(i4), parameter :: NG = 6      ! Degree of quadrature formula

! These are the 6-point-Gauss-Legendre abscissa (X-axis) values in [-1,1]:

  real(r8), parameter :: GX(ng) = (/ &
 &  -9.32469514203152027812d-1, -6.61209386466264513661d-1,   &
 &  -2.38619186083196908631d-1,  2.38619186083196908631d-1,   &
 &   6.61209386466264513661d-1,  9.32469514203152027812d-1 /)

! These are the corresponding 6-point-Gauss-Legendre Weights values:

  Real(r8), parameter :: Gw(ng) = (/ &
 &  -9.32469514203152027812d-1, -6.61209386466264513661d-1,   &
 &  -2.38619186083196908631d-1,  2.38619186083196908631d-1,   &
 &   6.61209386466264513661d-1,  9.32469514203152027812d-1 /)

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains
! =============================================     QUAD_START     =====
  subroutine QUAD_START ( A, B, WHAT, X, INTEGRAL )

  ! Start the process of estimating the integral of F using Gauss-Legendre
  ! quadrature.

    real(r4), intent(in) :: A, B       ! Limits of interval
    integer(i4), intent(out) :: WHAT   ! What to do?  See usage below
    real(r4), intent(out) :: X         ! Abscissa at which to evaluate F
    real(r4), intent(out) :: INTEGRAL  ! Value of integral -- set to zero

  ! Usage:
  !  call quad_start ( a, b, what, x, integral )
  !  do while ( what > 0 )
  !    Evaluate F at X
  !    call quad ( what, x, f, integral )
  !  end do
  !  INTEGRAL is the integral of F here

    ablen = 0.5 * ( b-a )
    abmid = 0.5 * ( a+b )
    integral = 0.0
    if ( ablen == 0.0 ) then
      what = 0
      return
    end if
    what = 1
    x = abmid + ablen * gx(1)
    return
  end subroutine QUAD_START

! ===================================================     QUAD     =====
  subroutine QUAD ( WHAT, X, F, INTEGRAL )

  ! Continue the process of estimating the integral of F using
  ! Gauss-Legendre quadrature.

    real(r4), intent(in) :: F           ! Value of F at X
    integer(i4), intent(inout) :: WHAT  ! What to do?  See usage below
    real(r4), intent(out) :: X          ! Abscissa at which to evaluate F
    real(r4), intent(inout) :: INTEGRAL ! Value of the integral

    if ( what <= 0 .or. what > ng ) then
      print *, '-E-GAUSS_LEGENDRE_M%QUAD entered with WHAT <= 0 or > ', ng
      stop
    end if

    integral = integral + ablen * gw(what) * f
    if ( what == ng ) then              ! Finished
      what = -ng
      return
    end if

    what = what + 1
    x = abmid + ablen * gx(what)
    return

  end subroutine QUAD
end module GAUSS_LEGENDRE_M

! $Log$
