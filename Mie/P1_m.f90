! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module P1_m

  implicit NONE

  private

  public Coeffs, Compute_P1, Compute_P1_Derivs

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  pure subroutine Coeffs ( C1, C2, W )

  ! Compute the coefficients C1(j) = (2*j-1)/(j-1),
  !                          C2(j) = j/(j-1)
  !                          W(j) = (2*j+1)/(j*(j+1)) = C1(j+1)/(j+1)

  ! C1 and C2 are used to compute Legendre functions of the first
  ! kind and degree 1, divided by sin(theta), using
  ! P1(0) = 0, P1(1) = 1, P1(j) = C1(j)*P1(j-1)*cos(theta) - C2(j)*P1(j-2)

  ! W is used to compute the scattering phase functions.

    use MLSKinds, only: RK => R8

    real(rk), intent(out) :: C1(2:), C2(2:), W(:) ! Assumed all the same ubound

    integer :: J

    w(1) = 1.5
    do j = 2, ubound(c1,1)
      c1(j) = (2*j-1.0_rk)/(j-1.0_rk)
      c2(j) = j/(j-1.0_rk)
      w(j) = (2*j+1.0_rk)/(j*(j+1.0_rk))
    end do

  end subroutine Coeffs

  pure subroutine Compute_P1 ( C1, C2, theta, P1 )

  !{ Evaluate Legendre functions of the first kind and degree 1,
  ! $P_j^1/\sin\theta$, using
  ! P1(0) = 0, P1(1) = 1, P1(j) = C1(j)*P1(j-1)*cos(theta) - C2(j)*P1(j-2)
  ! where P1(j) is $P_j^1(\cos\theta)$.

    use MLSKinds, only: RK => R8

    real(rk), intent(in) :: C1(2:), C2(2:) ! From Coeffs
    real(rk), intent(in) :: Theta
    real(rk), intent(out) :: P1(0:) ! Assumed same upper bound as C1 and C2

    integer :: J
    real(rk) :: CT

    ct = cos(theta)
    P1(0) = 0.0
    P1(1) = 1.0
    do j = 2, ubound(c1,1)
      P1(j) = c1(j) * P1(j-1) * ct - c2(j) * P1(j-2)
    end do

  end subroutine Compute_P1

  pure subroutine Compute_P1_Derivs ( C1, C2, theta, P1, dP1_dTheta )

  !{ Evaluate Legendre functions of the first kind and degree 1,
  ! $P_j^1/\sin\theta$, using
  ! $P1(0) = 0, P1(1) = 1, P1(j) = C1(j) P1(j-1) \cos(\theta) - C2(j) P1(j-2)$
  ! where $P1(j)$ is $P_j^1(\cos\theta)$.  Evaluate their derivatives using
  ! $dP_j^1 (\cos(\theta))/d\theta = j \cos(\theta) P(j) - (j+1) P(j-1)$

    use MLSKinds, only: RK => R8

    real(rk), intent(in) :: C1(2:), C2(2:) ! From Coeffs
    real(rk), intent(in) :: Theta
    real(rk), intent(out) :: P1(0:) ! Assumed same upper bound as C1 and C2
    real(rk), intent(out) :: dP1_dTheta(:)

    integer :: J
    real(rk) :: CT

    ct = cos(theta)
    P1(0) = 0.0
    P1(1) = 1.0
    dP1_dTheta(1) = ct
    do j = 2, ubound(c1,1)
      P1(j) = c1(j) * P1(j-1) * ct - c2(j) * P1(j-2)
      dP1_dTheta(j) = j * ct * P1(j) - (j+1) * P1(j-1)
    end do

  end subroutine Compute_P1_Derivs

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module P1_m

! $Log$
! Revision 1.2  2009/07/06 22:08:38  vsnyder
! Added print to not_used_here to preserve Id
!
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
