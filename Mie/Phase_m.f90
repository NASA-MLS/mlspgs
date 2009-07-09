! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Phase_m

  implicit NONE

  private

  public Phase, Phase_Deriv

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  pure function Phase ( Theta, A, B, C1, C2, W, P1, dP1_dTheta )

  !{ Evaluate the phase function
  ! $p(\theta,r) = C \left( | S_1(\theta,r) |^2 + | S_2(\theta,r) |^2 \right)$
  ! where
  ! \begin{equation}\begin{split}
  !  S_1(\theta,r) =\,& \sum_{j=1}^\infty \frac{2j+1}{j(j+1)} \left(
  !     a_j \frac{\text{d} P_j^1(\cos\theta)}{\text{d}\theta} +
  !     b_j \frac{P_j^1(\cos\theta)}{\sin\theta}\right) \text{ and} \\
  !  S_2(\theta,r) =\,& \sum_{j=1}^\infty \frac{2j+1}{j(j+1)} \left(
  !     a_j \frac{P_j^1(\cos\theta)}{\sin\theta} +
  !     b_j \frac{\text{d} P_j^1(\cos\theta)}{\text{d}\theta}\right) \\
  ! \end{split}\end{equation}
  ! and $P_j^1$ is the Legendre function of the first kind and order 1.
  ! $a_j$ and $b_j$ are coefficients used to compute the Mie efficiency
  ! (see wvs-066).  It is assumed that A and B are the same size.

    use MLSKinds, only: RK => R8
    use P1_m, only: Compute_P1_Derivs

    real(rk) :: Phase
    real(rk), intent(in) :: Theta
    complex(rk), intent(in) :: A(:), B(:)      ! Depend on R but not theta
    real(rk), intent(in) :: C1(2:), C2(2:), W(:) ! Coefficients that do not
                                               ! depend upon theta or R,
                                               ! From P1_m % Coeffs
    real(rk), intent(in), optional :: P1(0:)   ! Legendre function and
    real(rk), intent(in), optional :: dP1_dTheta(:) !   its derivative

    real(rk) :: My_P1(0:ubound(a,1)), My_dP1_dTheta(ubound(a,1))
    complex(rk) :: S1, S2

    integer :: J

    s1 = 0.0
    s2 = 0.0
    if ( present(P1) ) then
      do j = 1, ubound(a,1)
        s1 = s1 + w(j) * ( a(j) * dP1_dTheta(j) + b(j) * P1(j) )
        s2 = s2 + w(j) * ( a(j) * P1(j) + b(j) * dP1_dTheta(j) )
      end do
    else
      call compute_P1_derivs ( c1(2:ubound(a,1)), c2, theta, my_P1, my_dP1_dTheta )
      do j = 1, ubound(a,1)
        s1 = s1 + w(j) * ( a(j) * my_dP1_dTheta(j) + b(j) * my_P1(j) )
        s2 = s2 + w(j) * ( a(j) * my_P1(j) + b(j) * my_dP1_dTheta(j) )
      end do
    end if


    phase = real(s1)**2 + aimag(s1)**2 + real(s2)**2 + aimag(s2)**2

  end function Phase

  pure subroutine Phase_Deriv ( Theta, A, B, dA_dT, dB_dT, C1, C2, W, &
    &                           Phase, dPhase_dT, P1, dP1_dTheta )

  !{ Evaluate the phase function and its derivative w.r.t. T.

    use MLSKinds, only: RK => R8
    use P1_m, only: Compute_P1_Derivs

    real(rk), intent(in) :: Theta
    complex(rk), intent(in) :: A(:), B(:), dA_dT(:), dB_dT(:) ! Depend on R
    real(rk), intent(in) :: C1(2:), C2(2:), W(:) ! Coefficients that do not
                                               ! depend upon theta or R,
                                               ! From P1_m % Coeffs
    real(rk), intent(out) :: Phase, dPhase_dT
    real(rk), intent(in), optional :: P1(0:)   ! Legendre function and
    real(rk), intent(in), optional :: dP1_dTheta(:) !   its derivative

    real(rk) :: My_P1(0:ubound(a,1)), My_dP1_dTheta(ubound(a,1))
    complex(rk) :: S1, S2, dS1_dT, dS2_dT

    integer :: J

    s1 = 0.0
    s2 = 0.0
    dS1_dT = 0.0
    dS2_dT = 0.0
    if ( present(P1) ) then
      do j = 1, ubound(a,1)
        s1 = s1 + w(j) * ( a(j) * dP1_dTheta(j) + b(j) * P1(j) )
        s2 = s2 + w(j) * ( a(j) * P1(j) + b(j) * dP1_dTheta(j) )
        dS1_dT = dS1_dT + w(j) * ( da_dT(j) * dP1_dTheta(j) + db_dT(j) * P1(j) )
        dS2_dT = dS2_dT + w(j) * ( da_dT(j) * P1(j) + db_dT(j) * dP1_dTheta(j) )
      end do
    else
      call compute_P1_derivs ( c1(2:ubound(a,1)), c2, theta, my_P1, my_dP1_dTheta )
      do j = 1, ubound(a,1)
        s1 = s1 + w(j) * ( a(j) * my_dP1_dTheta(j) + b(j) * my_P1(j) )
        s2 = s2 + w(j) * ( a(j) * my_P1(j) + b(j) * my_dP1_dTheta(j) )
        dS1_dT = dS1_dT + w(j) * ( da_dT(j) * my_dP1_dTheta(j) + db_dT(j) * my_P1(j) )
        dS2_dT = dS2_dT + w(j) * ( da_dT(j) * my_P1(j) + db_dT(j) * my_dP1_dTheta(j) )
      end do
    end if

    phase = real(s1)**2 + aimag(s1)**2 + real(s2)**2 + aimag(s2)**2
    dPhase_dT = 2.0 * ( real(s1)*real(dS1_dT) + aimag(s1)*aimag(dS1_dT) + &
                        real(s2)*real(dS2_dT) + aimag(s2)*aimag(dS2_dT) )

  end subroutine Phase_Deriv

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Phase_m

! $Log$
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
