! Copyright 2009, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Polyfit_m

! Least-squares fit a set of data to a polynomial of specified degree.

  use MLSKinds, only: RK => R8

  implicit NONE
  private
  public :: Polyfit, RK

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine Polyfit ( X, Y, MaxOrd, Order, Stdev, Fit, Diff )

    real(rk), intent(in) :: X(:)       ! Independent variable
    real(rk), intent(in) :: Y(:)       ! size(y) <= 2 voids warranty
    integer, intent(in) :: MaxOrd      ! maximum order of the fit
    integer, intent(out) :: Order      ! actual order of the fit
    real(rk), intent(out) :: Stdev     ! norm2(residual)/(m-n)
    real(rk), intent(out) :: Fit(:)    ! The fitted polynomial evaluated at X
    real(rk), intent(out) :: Diff(:)   ! Y - Fit

    real(rk) :: P(maxOrd+3), SD(1), W((maxOrd+3)*(maxOrd+3))

    integer :: I, M

    interface PFIT
      subroutine DPFIT ( M, X, Y, SD, NMAX, SEEKN, COMTRN, CHBBAS, P, NDEG, &
        &                SIGFAC, W )
        integer, intent(in) :: M
        double precision, intent(in) :: X(m), Y(m), SD(*)
        integer, intent(in) :: NMAX
        logical, intent(in) :: SEEKN, COMTRN, CHBBAS
        double precision, intent(inout) :: P(nmax+3)
        integer, intent(out) :: NDEG
        double precision, intent(out) :: SIGFAC
        double precision, intent(inout) :: W((nmax+3)*(nmax+3))
      end subroutine DPFIT
      subroutine SPFIT ( M, X, Y, SD, NMAX, SEEKN, COMTRN, CHBBAS, P, NDEG, &
        &                SIGFAC, W )
        integer, intent(in) :: M
        real, intent(in) :: X(m), Y(m), SD(*)
        integer, intent(in) :: NMAX
        logical, intent(in) :: SEEKN, COMTRN, CHBBAS
        real, intent(inout) :: P(nmax+3)
        integer, intent(out) :: NDEG
        real, intent(out) :: SIGFAC
        real, intent(inout) :: W((nmax+3)*(nmax+3))
      end subroutine SPFIT
    end interface

    interface CPVAL
      double precision function DCPVAL ( P, NDEG, X )
        integer, intent(in) :: NDEG
        double precision, intent(in) :: P(ndeg+3)
        double precision, intent(in) :: X
      end function DCPVAL
      real function SCPVAL ( P, NDEG, X )
        integer, intent(in) :: NDEG
        real, intent(in) :: P(ndeg+3)
        real, intent(in) :: X
      end function SCPVAL
    end interface

    m = size(y)

    ! Compute the polynomial fit
    sd(1) = -1.0
    call pfit ( m, x, y, sd, maxord, .true., .true., .true., p, order, stdev, w )

    ! Evaluate the fitted polynomial
    do i = 1, m
      fit(i) = cpval(p,order,x(i))
    end do

    diff = y - fit

  end subroutine Polyfit

! ------------------------------------------------  not_used_here  -----
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Polyfit_m

! $Log$
