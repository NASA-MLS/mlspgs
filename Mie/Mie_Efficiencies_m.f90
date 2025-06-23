! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Mie_Efficiencies_m

!{ Compute Mie efficiencies 
! \begin{equation*}
! \xi_s = \frac2{\chi^2} \sum_{n=1}^{n_\text{cut}}
!       (2 n + 1 ) ( |a_n|^2 + |b_n|^2 ) \text{ and }
! \xi_e = \frac2{\chi^2} \sum_{n=1}^{n_\text{cut}}
!       (2 n + 1 ) \Re( a_n + b_n )
! \end{equation*}
!%
! where $\chi = \frac{2 \pi r}{\lambda}$, $r$ is particle radius,
!%
! $\lambda$ is wavelength,
! \begin{equation*}
! a_n = \frac{(A_n/m + n/\chi) \Re W_n - \Re W_{n-1}}
!            {(A_n/m + n/\chi) W_n - W_{n-1}}\,, \,
! b_n = \frac{(m A_n + n/\chi) \Re W_n - \Re W_{n-1}}
!            {(m A_n + n/\chi) W_n - W_{n-1}}\,,\\
! \end{equation*}
!%
! $W_n = \chi h^{(2)}_n (\chi)$, $h^{(2)}_n (\chi)$ is the spherical
! Hankel function of the second kind and order $n$ = $j_n(\chi) -
! i y_n(\chi)$, $j_n(\chi)$ and $y_n(\chi)$ are spherical Bessel functions
! of the first and second kind, respectively,
!%
! \begin{equation*}
! A_n = -\frac{n}{m \chi} + \frac{j_{n-1}(m \chi)}
!                                {j_n(m \chi)}\,,
! \end{equation*}
!%
! $m = \sqrt{\epsilon}$ is the complex refractive index, and $\epsilon$ is
! the complex dielectric constant.
! 
! $W_n$ satisfies the recurrence
! $W_n = \frac{2 n - 1}{\chi} W_{n-1} - W_{n-2}$
! with initial conditions $W_{-1} = e^{-i \chi}$ and
! $W_0 = i e^{-i \chi}$, which is unstable for the real part $j_n(\chi)$ in
! the forward direction, and unstable for the imaginary $y_n(\chi)$ part in
! the bacward direction.
!
! $A_n$ satisfies the recurrence
! $A_n = -\frac{n}{m \chi} + \left( \frac{n}{m \chi} - A_{n-1} \right)^{-1}$
! with initial condition $A_0 = \cot m \chi$, which is unstable in the forward
! direction.
!
! The subroutine {\tt SPHBESS\_Z} calculates $A_n$ and $W_n$ using recurrences
! carried in the correct directions.

  implicit NONE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Mie_Efficiencies ( M, Chi, Xi_e, Xi_s, A, B, Ord )

    use MLSKinds, only: R8
    use SphBess, only: SphBess_Z

    complex(r8), intent(in) :: M  ! Index of refraction
    real(r8), intent(in) :: Chi   ! Particle size parameter
    real(r8), intent(out) :: Xi_e ! Extinction efficiency
    real(r8), intent(out) :: Xi_s ! Scattering efficiency
    complex(r8), intent(out) :: A(:) ! Coefficient a_n above, a_j in ATBD
    complex(r8), intent(out) :: B(:) ! Coefficient b_n above, b_j in ATBD
    integer, intent(out), optional :: Ord ! Order of Bessel functions actually used

    complex(r8) :: Acap(0:size(a)) ! Coefficient A_n above, A_j in ATBD
    integer :: N                  ! Subscript, loop inductor
    real(r8) :: JR(0:size(a))     ! Spherical Bessel function first kind
    real(r8) :: R_Chi             ! 1 / chi
    complex(r8) :: T              ! A temporary variable
    complex(r8) :: U              ! 1 - 1/a_n or 1 - 1/b_n =
                   ! ( t*y_n(chi)) - y_{n-1}(chi)) ) /
                   ! ( t*j_n(chi)) - j_{n-1}(chi)) )
                   ! Remember: h^{(2)}_n(chi) = j_n(chi) - y_n(chi)
    real(r8) :: XA, XR            ! To be accumulated
    real(r8) :: Y(0:size(a))      ! Spherical Bessel function second kind =
                                  ! -aimag(h^{(2)}(chi))

    r_chi = 1.0 / chi
    call sphbess_z ( chi, m, acap, jr, y )
    xi_e = 0.0
    xi_s = 0.0
    do n = 1, size(a)
      t = Acap(n)/m + n * r_chi
      u = (t * y(n) - y(n-1)) / (t * jr(n) - jr(n-1))
      a(n) = 1.0 / ( 1.0 - cmplx(-aimag(u),real(u),kind=r8) )
      t = Acap(n)*m + n * r_chi
      u = (t * y(n) - y(n-1)) / (t * jr(n) - jr(n-1))
      b(n) = 1.0 / ( 1.0 - cmplx(-aimag(u),real(u),kind=r8) )
      xa = (2*n+1) * (abs(a(n))**2 + abs(b(n))**2)
      xr = (2*n+1) * (real(a(n)) + real(b(n)))
      if ( abs(xa) < abs(xi_s)*epsilon(xi_s) .and. &
        &  abs(xr) < abs(xi_e)*epsilon(xi_e) ) then
        a(n+1:) = 0.0
        b(n+1:) = 0.0
        exit
      end if
      xi_s = xi_s + xa
      xi_e = xi_e + xr
    end do
    xi_e = xi_e * 2 / chi**2
    xi_s = xi_s * 2 / chi**2
    if ( present(ord) ) ord = n - 1

  end subroutine Mie_Efficiencies

  subroutine Mie_Efficiencies_Derivs ( M, Chi, dM_dT, Xi_e, Xi_s, &
    & dXi_e_dT, dXi_s_dT, A, B, dA_dT, dB_dT, Ord, Acap_out, dAcap_out )

    use Constants, only: Pi
    use MLSKinds, only: R8
    use SphBess, only: SphBess_Z

    complex(r8), intent(in) :: M  ! Index of refraction
    real(r8), intent(in) :: Chi   ! Particle size parameter
    complex(r8), intent(in) :: dM_dT ! Derivative of M w.r.t temperature
    real(r8), intent(out) :: Xi_e ! Extinction efficiency
    real(r8), intent(out) :: Xi_s ! Scattering efficiency
    real(r8), intent(out) :: dXi_e_dT, dXi_s_dT ! Temperature derivatives
    complex(r8), intent(out) :: A(:) ! Coefficient a_n above, a_j in ATBD
    complex(r8), intent(out) :: B(:) ! Coefficient b_n above, b_j in ATBD
    complex(r8), intent(out) :: dA_dT(:), dB_dT(:) ! Temperature derivatives
    integer, intent(out), optional :: Ord ! Order of Bessel functions actually used
    complex(r8), intent(out), optional :: Acap_out(0:) ! See Acap
    complex(r8), intent(out), optional :: dAcap_out(:) ! dAcap_dM * dM_dT

    real(r8), parameter :: TwoDPi = 2.0_r8 / pi

    complex(r8) :: Acap(0:ubound(a,1)) ! Coefficient A_n above, A_j in ATBD
    complex(r8) :: AHat, BHat     ! \hat a_n, \hat b_n
    complex(r8) :: dAcap_dM       ! Derivative w.r.t M
    complex(r8) :: dAhat_dT, dBhat_dT ! Derivatives w.r.t. M
    real(r8) :: JR(0:size(a))     ! Spherical Bessel function first kind
    integer :: N                  ! Subscript, loop inductor
    real(r8) :: R_Chi             ! 1 / chi
    real(r8) :: R_Chi_Sq          ! 1 / chi**2
    complex(r8) :: R_M            ! 1/m
    complex(r8) :: R_M_Chi        ! 1/(m*chi)
    complex(r8) :: U              ! 1 - 1/a_n or 1 - 1/b_n =
                   ! ( aHat*y_n(chi)) - y_{n-1}(chi)) ) /
                   ! ( aHat*j_n(chi)) - j_{n-1}(chi)) ) or
                   ! ( bHat*y_n(chi)) - y_{n-1}(chi)) ) /
                   ! ( bHat*j_n(chi)) - j_{n-1}(chi)) )
                   ! Remember: h^{(2)}_n(chi) = j_n(chi) - y_n(chi)
    real(r8) :: XA, XR            ! To be accumulated
    real(r8) :: Y(0:size(a))      ! Spherical Bessel function second kind =
                                  ! -aimag(h^{(2)}(chi))

    r_chi = 1.0 / chi
    r_chi_sq = r_chi**2
    r_m_chi = 1.0 / (m * chi)
    r_m = 1.0 / m
    call sphbess_z ( chi, m, aCap, jr, y )
    if ( present(aCap_out) ) aCap_out = aCap
    xi_e = 0.0
    xi_s = 0.0
    dXi_e_dT = 0.0
    dXi_s_dT = 0.0
    do n = 1, size(a)
      ! State
      aHat = Acap(n) * r_m + n * r_chi                           ! a hat
      u = (aHat * y(n) - y(n-1)) / (aHat * jr(n) - jr(n-1))
      a(n) = 1.0 / ( 1.0 - cmplx(-aimag(u),real(u),kind=r8) )
      bHat = Acap(n) * m + n * r_chi                             ! b hat
      u = (bHat * y(n) - y(n-1)) / (bHat * jr(n) - jr(n-1))
      b(n) = 1.0 / ( 1.0 - cmplx(-aimag(u),real(u),kind=r8) )
      xa = (2*n+1) * (abs(a(n))**2 + abs(b(n))**2)
      xr = (2*n+1) * (real(a(n)) + real(b(n)))
      if ( abs(xa) < abs(xi_s)*epsilon(xi_s) .and. &
        &  abs(xr) < abs(xi_e)*epsilon(xi_e) ) then
        a(n+1:) = 0.0
        b(n+1:) = 0.0
        exit
      end if
      xi_e = xi_e + xr
      xi_s = xi_s + xa
      ! Derivatives
      dAcap_dM = r_m**2 * ( n * r_chi + (2*n - chi * bHat) * bHat ) - chi
      dAhat_dT = dM_dT * ( dAcap_dM - Acap(n) * r_m ) * r_m
      dA_dT(n) = cmplx( -r_chi_sq * aimag(dAhat_dT), &
        &                r_chi_sq * real(dAhat_dT), kind=r8 ) / &
        & ( aHat * cmplx(jr(n),-y(n),kind=r8) - cmplx(jr(n-1),-y(n-1),kind=r8) )**2
      dBhat_dT = dM_dT * ( m * dAcap_dM + Acap(n) )
      dB_dT(n) = cmplx( -r_chi_sq * aimag(dBhat_dT), &
        &                r_chi_sq * real(dBhat_dT), kind=r8 ) / &
        & ( bHat * cmplx(jr(n),-y(n),kind=r8) - cmplx(jr(n-1),-y(n-1),kind=r8) )**2
      dXi_e_dT = dXi_e_dT + (2*n+1) * ( real(da_dt(n)) + real(db_dt(n)) )
      dXi_s_dT = dXi_s_dT + (2*n+1) * &
        & ( real(a(n))*real(da_dt(n)) + aimag(a(n))*aimag(da_dt(n)) + &
        &   real(b(n))*real(db_dt(n)) + aimag(b(n))*aimag(db_dt(n)) )
      if ( present(dAcap_out) ) dAcap_out(n) = dAcap_dM * dM_dT
    end do
    xi_e = xi_e * 2 * r_chi_sq
    xi_s = xi_s * 2 * r_chi_sq
    dXi_e_dT = dXi_e_dT * 2 * r_chi_sq
    dXi_s_dT = dXi_s_dT * 4 * r_chi_sq
    if ( present(ord) ) ord = n - 1

  end subroutine Mie_Efficiencies_Derivs

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Mie_Efficiencies_m

! $Log$
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
