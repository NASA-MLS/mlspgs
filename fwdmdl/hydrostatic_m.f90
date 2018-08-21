! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Hydrostatic_m

  implicit none

  private
  public :: Hydrostatic

  interface Hydrostatic
    module procedure Hydrostatic_All, Hydrostatic_All_ZZ, Hydrostatic_No_Der
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  contains
!---------------------------------------------------------------------------

  ! --------------------------------------------  Hydrostatic_All  -----
  subroutine Hydrostatic_All ( lat, t_basis, t_coeffs, z_grid, z_ref, h_ref, &
                      &    t_grid, h_grid, dhidzi, dhidtq, ddhdhdtq, z_surface )

! Compute a hydrostatic function per L2PC method and return
! geometric heights. Reference height is now an input

    use MLSKinds, only: RP, IP
    use Geometry, only: EarthRadA, EarthRadB, GM, J2, J4, W
    use Piq_Int_m, only: Piq_Int
    use Physics, only: BoltzMeters => Boltz ! Avogadro * k * ln10 / mmm m^2/(K s^2)
    use Sparse_eta_m, only: Sparse_Eta_t
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

! Inputs

    real(rp), intent(in) :: lat         ! geocentric latitude in radians
    real(rp), intent(in) :: t_basis(:)  ! vertical temperature basis, zeta
    real(rp), intent(in) :: t_coeffs(:) ! temperature values
    real(rp), intent(in) :: z_grid(:)   ! zeta = -log10(P) pressures for which
!                                         heights are needed
    real(rp), intent(in) :: z_ref       ! reference pressure in zeta = -log10(P)
    real(rp), intent(in) :: h_ref       ! reference geopotential height in km

! Outputs

    real(rp), intent(out) :: t_grid(:)  ! temperatures on z_grid
    real(rp), intent(out) :: h_grid(:)  ! heights on z_grid (km)
    real(rp), intent(out) :: dhidzi(:)  ! dh/dz on z_grid

    real(rp), optional, intent(out) :: dhidtq(:,:) ! dh/dt on z_grid assuming
      ! dh/dt = 0.0 at the reference ellipse surface  equivalent to h_ref = 0.0
    real(rp), optional, intent(out) :: ddhdhdtq(:,:) ! ddh/dhdt on z_grid,
      ! needs dhidtq
    real(rp), optional, intent(out) :: z_surface

! Internal stuff

    real(rp), parameter :: ERadAsq = EarthRadA**2, ERadBsq = EarthRadB**2
    ! Boltzmann constant in kilometers:
    real(rp), parameter :: Boltz = boltzMeters/1.0e6_rp ! = kln10/m km^2/(K sec^2)

    integer :: Me = -1          ! String index for trace
    integer(ip) :: n_coeffs,iter

!   real(rp) :: cl, sl ! for derivatives of Legendre polynomials dp2 and dp4
    real(rp) :: Clsq, G_ref, GHB
    real(rp) :: R_e, R_eisq, R_eff, Slsq, Z_surf
    real(rp) :: P2, P4    ! Legendre polynomials, for oblateness model
!   real(rp) :: dp2, dp4  ! Derivatives of P2 and P4
    real(rp) :: dh_dz_S, H_calc, Z_old
    real(rp), dimension(size(h_grid),size(t_basis)) :: Piq
    real(rp), dimension(1,size(t_coeffs)) :: Piqa, Piqb
    real(rp), dimension(size(z_grid)) :: Mass_corr
    type(sparse_eta_t) :: Eta

! begin the code

    call trace_begin ( me, 'Hydrostatic_All', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

    n_coeffs = size(t_basis)

    where ( z_grid > 2.5_rp )
!     mass_corr = 1.0_rp / (0.875_rp + 0.1_rp*z_grid - 0.02_rp*z_grid**2)

!{ A series expansion about z\_grid = 5/2 of
!  $\frac1{\frac78 + \frac1{10}z - \frac1{50}z^2}$ is
!  $\sum_{k=0}^{\infty} \left( \frac{(z-\frac52)^2}{50}\right)^k$.  Use the
!  first two terms.
      mass_corr = 1.0_rp + 0.02_rp*(z_grid - 2.5_rp)**2
    elsewhere
      mass_corr = 1.0_rp
    end where

! compute Eta (interpolation coefficients) from T_Basis to Z_Grid, and
! use them to compute t_grid
    call eta%eta_1d ( t_basis, z_grid, create=.true., sorted=.false. )
    call eta%sparse_dot_vec ( t_coeffs, t_grid )

!{ Compute surface acceleration and effective earth radius.  First,
!  evaluate the Legendre polynomials $P_2(\sin\lambda) = \frac32 \sin^2\lambda
!  - \frac12$ and $P_4(\sin\lambda) = \frac{35}8 \sin^4\lambda
!  - \frac{15}4 \sin^2\lambda + \frac38$.

    slsq = sin(lat) ** 2
    clsq = 1.0_rp - slsq
    p2 = (3.0_rp * slsq - 1.0_rp) * 0.5_rp
    p4 = (35.0_rp * slsq**2 - 30.0_rp * slsq + 3.0_rp) * 0.125_rp
!   dp2 = 3.0_rp * sl * cl
!   dp4 = 2.5_rp * sl * cl * ( 7.0_rp * slsq - 3.0_rp )

!{ Compute radius at geoid having potential = $W_0$ (see Geometry module):
!  $r_e^2 = \frac{a^2 b^2}{a^2 \sin^2 \lambda + b^2 \cos^2 \lambda}$ in meters.
!  First compute $r_e^{-2}$ in $m^{-2}$.  The parenthesization used here
!  avoids potential overflow; $a^2 b^2 \approx 1.64 \times 10^{27}$, but
!  the inner factor is $\approx 1$ and $a^2 \approx 4 \times 10^{13}$.

    r_eisq = ( (eRadAsq*slsq + eRadBsq*clsq) / eRadBsq ) / eRadAsq ! (r_e)^{-2}

    r_e = sqrt(1.0_rp / r_eisq) ! in meters

!{ Radial surface acceleration at latitude $\lambda$, k$m/s^2$:
!  \begin{equation*}
!   g_{\text{ref}} = 0.001 \left ( G m
!                     \frac{1 - 3 J_2 P_2(\sin \lambda)
!                       (a/r_e)^2 - 5 J_4 P_4(\sin \lambda) (a/r_e)^4 }{r_e^2}
!                       - \omega^2 \cos^2 \lambda \, r_e \right )
!  \end{equation*}

    g_ref = 0.001_rp * (gm * (1.0_rp - 3.0_rp*j2*p2*eRadAsq*r_eisq &
        & - 5.0_rp*j4*p4*(eRadAsq*r_eisq)**2)*r_eisq - w**2 * clsq * r_e)

!{ Better effective Earth radius:
!  compute $-2\, g_{\text{ref}} / ( d g_{\text{ref}} / d r)$, kilometers:
!  \begin{equation*}
!   r_{\text{eff}} =
!    \frac{2\, g_{\text{ref}}}
!         { 2\, G m \frac{1 - 6 J_2 P_2 (\sin \lambda) (a/r_e)^2
!                         - 15 J_4 P_4 (\sin \lambda) (a/r_e)^4}{r_e^3}
!                          + \omega^2 \cos^2 \lambda }
!  \end{equation*}
!  $d g_{\text{ref}} / d r$ has units $s^{-2}$, because we compute
!  $d g_{\text{ref}}(\text{meters}) / d r$(meters), thereby canceling
!  the factor of 0.001 you might have been expecting in the denominator.
!  The entire expression has units (k$ms^{-2})/s^{-2} = $ k$m$.

    r_eff = 2.0_rp * g_ref / (2.0_rp * gm * (1.0_rp-6.0_rp*j2*p2 &
        & * eRadAsq*r_eisq - 15.0_rp*j4*p4*(eRadAsq*r_eisq)**2) &
        & / r_e**3 + w**2 * clsq)

! find the surface pressure

    ghb = g_ref * h_ref / boltz
    z_old = z_grid(1) ! This is a guess

    iter = 0
    do
!{ Newton iteration for $z_{\text{surf}}$ at $h_{\text{calc}} = -h_{\text{ref}}$:
!  $z_{n+1} = z_n - (h_{\text{ref}} + h_{\text{calc}}) /
!             \left . \frac{\text{d} h}{\text{d} z} \right |_{z=z_n}$

      call piq_int ( (/z_old/), t_basis, z_ref, piqa )
      h_calc = dot_product(piqa(1,:), t_coeffs)

      call piq_int ( (/z_old+0.01_rp/), t_basis, z_ref, piqa )
      call piq_int ( (/z_old-0.01_rp/), t_basis, z_ref, piqb )
      dh_dz_s = dot_product((piqa(1,:) - piqb(1,:)), t_coeffs) * 50.0_rp
      z_surf = z_old - (ghb + h_calc) / dh_dz_s

      iter = iter + 1
      if ( abs(z_surf - z_old) < 0.0001_rp .or. iter == 10 ) exit
      z_old = z_surf

    end do

    if (present(z_surface)) z_surface = z_surf

!{ Compute the {\tt piq} integrals with mass reduction compensation relative to the
!  surface.  Here, {\tt piq} = $P_l(\zeta) = P(\zeta_l,\phi_m)$.

    call piq_int ( z_grid, t_basis, z_surf, piq, Z_MASS=2.5_rp, C_MASS=0.02_rp )

! compute the height vector

! geopotential height * g_ref

!{ Equation (5.27) in the 19 August 2004 ATBD is
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{g_0 \stackrel{\star}{R_0^2}}
!         {g_0 \stackrel{\star}{R_0} - k\, \ln 10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!           ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!    -\stackrel{\star}{R_0} + R_0 - R^\oplus
!  \end{equation*}
!  The last two terms, i.e., $R_0-R^\oplus$, are handled in {\tt metrics}.
!  Putting what remains over a common denominator, we have
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{\stackrel{\star}{R_0} k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!         {g_0 \stackrel{\star}{R_0} - k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!  \end{equation*}
!  Notice that $g_0$ here is not the average or equatorial surface
!  acceleration.  Rather, it is $g_\text{ref}$, which is corrected for
!  latitude, the Earth's figure (up to fourth order zonal harmonics), and
!  centripetal acceleration.
    h_grid = boltz * matmul(piq,t_coeffs)
    h_grid = r_eff * h_grid / (r_eff * g_ref - h_grid)

!{ \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\text{d}h(\zeta,\phi)}{\text{d}\zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!   \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} f^T_{lm} \eta^T_m(\phi)
!    \frac{\partial P_l(\zeta)}{\partial \zeta}\,.
!  \end{equation*}
!
!  Because
!  \begin{equation*}
!  P_l(\zeta) = \int_{\zeta_0}^\zeta \frac{\eta^T_l(z)}{\mathcal{M}(z)}
!                 \text{d}z\,\text{, then }
!  \frac{\partial P_l(\zeta)}{\partial \zeta} =
!  \frac{\eta^T_l(\zeta)}{\mathcal{M}(\zeta)}\, \text{. Therefore,}
!  \end{equation*}
!  \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\partial h(\zeta,\phi)}{\partial \zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!  \frac1{\mathcal{M}(\zeta)} \left(
!   \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} f^T_{lm} \eta^T_m(\phi)
!    \eta^T_l(\zeta) \right)\,.
!  \end{equation*}
!  The factor in parentheses is simply $T(\zeta,\phi)$.  Therefore,
!  \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\partial h(\zeta,\phi)}{\partial \zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!  \frac{T(\zeta,\phi)}{\mathcal{M}(\zeta)}\,.
!  \end{equation*}
!  The variable {\tt Boltz} is $k \ln 10$.  Therefore
!  \begin{equation*}
!  \frac{\partial h(\zeta,\phi)}{\partial \zeta} = \left(
!  h(\zeta,\phi)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \right) \,\frac{T(\zeta,\phi)}{\mathcal{M}(\zeta)}\,.
!  \end{equation*}
!  At first, {\tt dhidzi} is only the factor in parentheses, which is
!  used to compute {\tt dhidtq}.

    dhidzi = (h_grid+r_eff)**2 * boltz / (g_ref * r_eff**2)

    if ( present(dhidtq) ) then

!{ \begin{equation*}
!  \frac{\partial h(\zeta,\phi)}{\partial f^T_{lm}} =
!  h(\zeta,\phi)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} P_l(\zeta) \eta^T_m(\phi)
!  \end{equation*}

      dhidtq = spread(dhidzi,2,n_coeffs) * piq
! this derivative is useful for antenna derivatives
      if ( present(ddhdhdtq) ) then
        ddhdhdtq = (2.0_rp/(spread(h_grid,2,n_coeffs)+r_eff)) * dhidtq
        block
          integer :: i, j, k
          do k = 1, eta%ne
            i = eta%e(k)%r
            j = eta%e(k)%c
            ddhdhdtq(i,j) = ddhdhdtq(i,j) + eta%e(k)%v / t_grid(i)
          end do
        end block
      end if
    end if

!{ \begin{equation*}
!  \frac{\partial h_{lm}}{\partial \zeta_{\,l}} =
!  h(\zeta_{\,l},\phi_m)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \frac{f^T_{lm}}{\mathcal{M}(\zeta_{\,l})}\,.
!  \end{equation*}

    dhidzi = dhidzi * t_grid * mass_corr

!{ If we ever need it:
!  \begin{equation*}\begin{split}
!  \frac{\partial^2 h_{lm}}{\partial \zeta_l \partial f^T_{lm}} = \,&
!  2 h(\zeta_{\,l},\phi_m) \frac{\partial h_{lm}}{\partial \zeta_{\,l}}
!  \text{\tt Boltz} \frac{T_{lm}}{\mathcal{M}(\zeta_{\,l})} +
!  h(\zeta_{\,l},\phi_m)^2\, \frac{\partial T(\zeta,\phi)}{\partial f^T_{lm}}
!  \text{\tt boltz}\,\frac1{\mathcal{M}(\zeta_{\,l})} \\
!  = \,&
!  2 h(\zeta_{\,l},\phi_m) \frac{\partial h_{lm}}{\partial \zeta_{\,l}}
!  \text{\tt Boltz} \frac{T_{lm}}{\mathcal{M}(\zeta_{\,l})} +
!  h(\zeta_{\,l},\phi_m)^2\, \text{\tt boltz}\,
!   \frac{\eta^T_{lm}}{\mathcal{M}(\zeta_{\,l})} \\
!  \end{split}\end{equation*}

    call trace_end ( 'Hydrostatic_All', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

  end subroutine Hydrostatic_All

  ! -----------------------------------------  Hydrostatic_All_ZZ  -----
  subroutine Hydrostatic_All_ZZ ( Lat, T_Basis, T_Coeffs, Z_Grid, Z_Ref, H_Ref, &
    &              Eta_ZZ, T_Grid, H_Grid, dHidZi, dHidTq, ddHdHdTq, Z_Surface )

! Compute a hydrostatic function per L2PC method and return
! geometric heights. Reference height is now an input

    use Geometry, only: EarthRadA, EarthRadB, GM, J2, J4, W
    use MLSKinds, only: RP, IP
    use Physics, only: BoltzMeters => Boltz ! Avogadro * k * ln10 / mmm m^2/(K s^2)
    use Piq_int_m, only: piq_int
    use Sparse_m, only: Sparse_t
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

! Inputs

    real(rp), intent(in) :: Lat         ! geocentric latitude in radians
    real(rp), intent(in) :: T_Basis(:)  ! vertical temperature basis, zeta
    real(rp), intent(in) :: T_Coeffs(:) ! temperature values
    real(rp), intent(in) :: Z_Grid(:)   ! zeta = -log10(P) pressures for which
!                                         heights are needed
    real(rp), intent(in) :: Z_Ref       ! reference pressure in zeta = -log10(P)
    real(rp), intent(in) :: H_Ref       ! reference geopotential height in km
    class(sparse_t), intent(in) :: Eta_ZZ ! Interpolation coefficients
                                        ! from T_Basis to Z_Grid

! Outputs

    real(rp), intent(out) :: T_Grid(:)  ! temperatures on z_grid
    real(rp), intent(out) :: H_Grid(:)  ! heights on z_grid (km)
    real(rp), intent(out) :: dHidZi(:)  ! dh/dz on z_grid

    real(rp), optional, intent(out) :: dHidTq(:,:) ! dH/dT on z_grid assuming
      ! dH/dT = 0.0 at the reference ellipse surface equivalent to H_Ref = 0.0
    real(rp), optional, intent(out) :: ddHdHdTq(:,:) ! ddH/dHdT on Z_Grid,
      ! needs dHidTq
    real(rp), optional, intent(out) :: Z_Surface

! Internal stuff

    real(rp), parameter :: ERadAsq = EarthRadA**2, ERadBsq = EarthRadB**2
    ! Boltzmann constant in kilometers:
    real(rp), parameter :: Boltz = boltzMeters/1.0e6_rp ! = kln10/m km^2/(K sec^2)

    integer :: I, J, K
    integer :: Me = -1          ! String index for trace
    integer(ip) :: N_Coeffs, Iter

!   real(rp) :: cl, sl ! for derivatives of Legendre polynomials dp2 and dp4
    real(rp) :: Clsq, G_ref, GHB
    real(rp) :: R_e, R_eisq, R_eff, Slsq, Z_surf
    real(rp) :: P2, P4    ! Legendre polynomials, for oblateness model
!   real(rp) :: dp2, dp4  ! Derivatives of P2 and P4
    real(rp) :: dH_dz_S, H_calc, Z_old
    real(rp), dimension(size(h_grid),size(t_basis)) :: Piq
    real(rp), dimension(1,size(t_coeffs)) :: Piqa, Piqb
    real(rp), dimension(size(z_grid)) :: Mass_corr

! begin the code

    call trace_begin ( me, 'Hydrostatic_All_ZZ', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

    n_coeffs = size(t_basis)

    where ( z_grid > 2.5_rp )
!     mass_corr = 1.0_rp / (0.875_rp + 0.1_rp*z_grid - 0.02_rp*z_grid**2)

!{ A series expansion about z\_grid = 5/2 of
!  $\frac1{\frac78 + \frac1{10}z - \frac1{50}z^2}$ is
!  $\sum_{k=0}^{\infty} \left( \frac{(z-\frac52)^2}{50}\right)^k$.  Use the
!  first two terms.
      mass_corr = 1.0_rp + 0.02_rp*(z_grid - 2.5_rp)**2
    elsewhere
      mass_corr = 1.0_rp
    end where

! compute t_grid using a sparse matrix-vector multiply
    call eta_zz%sparse_dot_vec ( t_coeffs, t_grid )

!{ Compute surface acceleration and effective earth radius.  First,
!  evaluate the Legendre polynomials $P_2(\sin\lambda) = \frac32 \sin^2\lambda
!  - \frac12$ and $P_4(\sin\lambda) = \frac{35}8 \sin^4\lambda
!  - \frac{15}4 \sin^2\lambda + \frac38$.

    slsq = sin(lat) ** 2
    clsq = 1.0_rp - slsq
    p2 = (3.0_rp * slsq - 1.0_rp) * 0.5_rp
    p4 = (35.0_rp * slsq**2 - 30.0_rp * slsq + 3.0_rp) * 0.125_rp
!   dp2 = 3.0_rp * sl * cl
!   dp4 = 2.5_rp * sl * cl * ( 7.0_rp * slsq - 3.0_rp )

!{ Compute radius at geoid having potential = $W_0$ (see Geometry module):
!  $r_e^2 = \frac{a^2 b^2}{a^2 \sin^2 \lambda + b^2 \cos^2 \lambda}$ in meters.
!  First compute $r_e^{-2}$ in $m^{-2}$.  The parenthesization used here
!  avoids potential overflow; $a^2 b^2 \approx 1.64 \times 10^{27}$, but
!  the inner factor is $\approx 1$ and $a^2 \approx 4 \times 10^{13}$.

    r_eisq = ( (eRadAsq*slsq + eRadBsq*clsq) / eRadBsq ) / eRadAsq ! (r_e)^{-2}

    r_e = sqrt(1.0_rp / r_eisq) ! in meters

!{ Radial surface acceleration at latitude $\lambda$, k$m/s^2$:
!  \begin{equation*}
!   g_{\text{ref}} = 0.001 \left ( G m
!                     \frac{1 - 3 J_2 P_2(\sin \lambda)
!                       (a/r_e)^2 - 5 J_4 P_4(\sin \lambda) (a/r_e)^4 }{r_e^2}
!                       - \omega^2 \cos^2 \lambda \, r_e \right )
!  \end{equation*}

    g_ref = 0.001_rp * (gm * (1.0_rp - 3.0_rp*j2*p2*eRadAsq*r_eisq &
        & - 5.0_rp*j4*p4*(eRadAsq*r_eisq)**2)*r_eisq - w**2 * clsq * r_e)

!{ Better effective Earth radius:
!  compute $-2\, g_{\text{ref}} / ( d g_{\text{ref}} / d r)$, kilometers:
!  \begin{equation*}
!   r_{\text{eff}} =
!    \frac{2\, g_{\text{ref}}}
!         { 2\, G m \frac{1 - 6 J_2 P_2 (\sin \lambda) (a/r_e)^2
!                         - 15 J_4 P_4 (\sin \lambda) (a/r_e)^4}{r_e^3}
!                          + \omega^2 \cos^2 \lambda }
!  \end{equation*}
!  $d g_{\text{ref}} / d r$ has units $s^{-2}$, because we compute
!  $d g_{\text{ref}}(\text{meters}) / d r$(meters), thereby canceling
!  the factor of 0.001 you might have been expecting in the denominator.
!  The entire expression has units (k$ms^{-2})/s^{-2} = $ k$m$.

    r_eff = 2.0_rp * g_ref / (2.0_rp * gm * (1.0_rp-6.0_rp*j2*p2 &
        & * eRadAsq*r_eisq - 15.0_rp*j4*p4*(eRadAsq*r_eisq)**2) &
        & / r_e**3 + w**2 * clsq)

! find the surface pressure

    ghb = g_ref * h_ref / boltz
    z_old = z_grid(1) ! This is a guess

    iter = 0
    do
!{ Newton iteration for $z_{\text{surf}}$ at $h_{\text{calc}} = -h_{\text{ref}}$:
!  $z_{n+1} = z_n - (h_{\text{ref}} + h_{\text{calc}}) /
!             \left . \frac{\text{d} h}{\text{d} z} \right |_{z=z_n}$

      call piq_int ( (/z_old/), t_basis, z_ref, piqa )
      h_calc = dot_product(piqa(1,:), t_coeffs)

      call piq_int ( (/z_old+0.01_rp/), t_basis, z_ref, piqa )
      call piq_int ( (/z_old-0.01_rp/), t_basis, z_ref, piqb )
      dh_dz_s = dot_product((piqa(1,:) - piqb(1,:)), t_coeffs) * 50.0_rp
      z_surf = z_old - (ghb + h_calc) / dh_dz_s

      iter = iter + 1
      if ( abs(z_surf - z_old) < 0.0001_rp .or. iter == 10 ) exit
      z_old = z_surf

    end do

    if (present(z_surface)) z_surface = z_surf

!{ Compute the {\tt piq} integrals with mass reduction compensation relative to the
!  surface.  Here, {\tt piq} = $P_l(\zeta) = P(\zeta_l,\phi_m)$.

    call piq_int ( z_grid, t_basis, z_surf, piq, Z_MASS=2.5_rp, C_MASS=0.02_rp )

! compute the height vector

! geopotential height * g_ref

!{ Equation (5.27) in the 19 August 2004 ATBD is
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{g_0 \stackrel{\star}{R_0^2}}
!         {g_0 \stackrel{\star}{R_0} - k\, \ln 10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!           ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!    -\stackrel{\star}{R_0} + R_0 - R^\oplus
!  \end{equation*}
!  The last two terms, i.e., $R_0-R^\oplus$, are handled in {\tt metrics}.
!  Putting what remains over a common denominator, we have
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{\stackrel{\star}{R_0} k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!         {g_0 \stackrel{\star}{R_0} - k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!  \end{equation*}
!  Notice that $g_0$ here is not the average or equatorial surface
!  acceleration.  Rather, it is $g_\text{ref}$, which is corrected for
!  latitude, the Earth's figure (up to fourth order zonal harmonics), and
!  centripetal acceleration.
    h_grid = boltz * matmul(piq,t_coeffs)
    h_grid = r_eff * h_grid / (r_eff * g_ref - h_grid)

!{ \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\text{d}h(\zeta,\phi)}{\text{d}\zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!   \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} f^T_{lm} \eta^T_m(\phi)
!    \frac{\partial P_l(\zeta)}{\partial \zeta}\,.
!  \end{equation*}
!
!  Because
!  \begin{equation*}
!  P_l(\zeta) = \int_{\zeta_0}^\zeta \frac{\eta^T_l(z)}{\mathcal{M}(z)}
!                 \text{d}z\,\text{, then }
!  \frac{\partial P_l(\zeta)}{\partial \zeta} =
!  \frac{\eta^T_l(\zeta)}{\mathcal{M}(\zeta)}\, \text{. Therefore,}
!  \end{equation*}
!  \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\partial h(\zeta,\phi)}{\partial \zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!  \frac1{\mathcal{M}(\zeta)} \left(
!   \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} f^T_{lm} \eta^T_m(\phi)
!    \eta^T_l(\zeta) \right)\,.
!  \end{equation*}
!  The factor in parentheses is simply $T(\zeta,\phi)$.  Therefore,
!  \begin{equation*}
!  \frac1{h(\zeta,\phi)^2} \frac{\partial h(\zeta,\phi)}{\partial \zeta} =
!  \frac{k \ln 10}{g_0 \stackrel{\star}{R_0^2}}
!  \frac{T(\zeta,\phi)}{\mathcal{M}(\zeta)}\,.
!  \end{equation*}
!  The variable {\tt Boltz} is $k \ln 10$.  Therefore
!  \begin{equation*}
!  \frac{\partial h(\zeta,\phi)}{\partial \zeta} = \left(
!  h(\zeta,\phi)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \right) \,\frac{T(\zeta,\phi)}{\mathcal{M}(\zeta)}\,.
!  \end{equation*}
!  At first, {\tt dhidzi} is only the factor in parentheses, which is
!  used to compute {\tt dhidtq}.

    dhidzi = (h_grid+r_eff)**2 * boltz / (g_ref * r_eff**2)

    if ( present(dhidtq) ) then

!{ \begin{equation*}
!  \frac{\partial h(\zeta,\phi)}{\partial f^T_{lm}} =
!  h(\zeta,\phi)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T} P_l(\zeta) \eta^T_m(\phi)
!  \end{equation*}

      dhidtq = spread(dhidzi,2,n_coeffs) * piq
! this derivative is useful for antenna derivatives
      if ( present(ddhdhdtq) ) then
        ddHdHdTq = (2.0_rp/(spread(h_grid,2,n_coeffs)+r_eff)) * dhidtq
        do k = 1, eta_zz%ne ! All the nonzero elements
            i = eta_zz%e(k)%r
            j = eta_zz%e(k)%c
            ddHdHdTq(i,j) = ddHdHdTq(i,j) + eta_zz%e(k)%v / t_grid(i)
        end do
      end if
    end if

!{ \begin{equation*}
!  \frac{\partial h_{lm}}{\partial \zeta_{\,l}} =
!  h(\zeta_{\,l},\phi_m)^2\, \frac{\text{\tt Boltz}}{g_0 \stackrel{\star}{R_0^2}}
!  \frac{f^T_{lm}}{\mathcal{M}(\zeta_{\,l})}\,.
!  \end{equation*}

    dhidzi = dhidzi * t_grid * mass_corr

!{ If we ever need it:
!  \begin{equation*}\begin{split}
!  \frac{\partial^2 h_{lm}}{\partial \zeta_l \partial f^T_{lm}} = \,&
!  2 h(\zeta_{\,l},\phi_m) \frac{\partial h_{lm}}{\partial \zeta_{\,l}}
!  \text{\tt Boltz} \frac{T_{lm}}{\mathcal{M}(\zeta_{\,l})} +
!  h(\zeta_{\,l},\phi_m)^2\, \frac{\partial T(\zeta,\phi)}{\partial f^T_{lm}}
!  \text{\tt boltz}\,\frac1{\mathcal{M}(\zeta_{\,l})} \\
!  = \,&
!  2 h(\zeta_{\,l},\phi_m) \frac{\partial h_{lm}}{\partial \zeta_{\,l}}
!  \text{\tt Boltz} \frac{T_{lm}}{\mathcal{M}(\zeta_{\,l})} +
!  h(\zeta_{\,l},\phi_m)^2\, \text{\tt boltz}\,
!   \frac{\eta^T_{lm}}{\mathcal{M}(\zeta_{\,l})} \\
!  \end{split}\end{equation*}

    call trace_end ( 'Hydrostatic_All_ZZ', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

  end subroutine Hydrostatic_All_ZZ

  ! -----------------------------------------  Hydrostatic_No_Der  -----
  subroutine Hydrostatic_No_Der ( lat, t_basis, t_coeffs, z_grid, z_ref, h_ref, &
                               &  h_grid )

! Compute a hydrostatic function per L2PC method and return
! geometric heights. Reference height is now an input

    use MLSCommon, only: RP, IP
    use Geometry, only: EarthRadA, EarthRadB, GM, J2, J4, W
    use Piq_int_m, only: piq_int
    use Physics, only: BoltzMeters => Boltz

! Inputs

    real(rp), intent(in) :: lat         ! geocentric latitude in radians
    real(rp), intent(in) :: t_basis(:)  ! vertical temperature basis, zeta
    real(rp), intent(in) :: t_coeffs(:) ! temperature values
    real(rp), intent(in) :: z_grid(:)   ! zeta = -log10(P) pressures for which
!                                         heights are needed
    real(rp), intent(in) :: z_ref       ! reference pressure in zeta = -log10(P)
    real(rp), intent(in) :: h_ref       ! reference geopotential height in km

! Outputs

    real(rp), intent(out) :: h_grid(:)  ! heights on z_grid (km)

! Internal stuff

    real(rp), parameter :: ERadAsq = EarthRadA**2, ERadBsq = EarthRadB**2
    ! Boltzmann constant in kilometers:
    real(rp), parameter :: Boltz = boltzMeters/1.0e6_rp ! = kln10/m km^2/(K sec^2)

    integer(ip) :: iter

!   real(rp) :: cl, sl ! for derivatives of Legendre polynomials dp2 and dp4
    real(rp) :: Clsq, G_ref, GHB
    real(rp) :: R_e, R_eisq, R_eff, Slsq, Z_surf
    real(rp) :: P2, P4    ! Legendre polynomials, for oblateness model
!   real(rp) :: dp2, dp4  ! Derivatives of P2 and P4
    real(rp) :: dh_dz_S, H_calc, Z_old
    real(rp), dimension(size(z_grid),size(t_basis)) :: Piq
    real(rp), dimension(1,size(t_basis)) :: Piqa, Piqb
!   real(rp), dimension(size(z_grid)) :: Mass_corr

! begin the code


!     where ( z_grid > 2.5_rp )
! !     mass_corr = 1.0_rp / (0.875_rp + 0.1_rp*z_grid - 0.02_rp*z_grid**2)

!{ A series expansion about z\_grid = 5/2 of
!  $\frac1{\frac78 + \frac1{10}z - \frac1{50}z^2}$ is
!  $\sum_{k=0}^{\infty} \left( \frac{(z-\frac52)^2}{50}\right)^k$.  Use the
!  first two terms.

!       mass_corr = 1.0_rp + 0.02_rp*(z_grid - 2.5_rp)**2
!     elsewhere
!       mass_corr = 1.0_rp
!     end where

!{ Compute surface acceleration and effective earth radius.  First,
!  evaluate the Legendre polynomials $P_2(\sin\lambda) = \frac32 \sin^2\lambda
!  - \frac12$ and $P_4(\sin\lambda) = \frac{35}8 \sin^4\lambda
!  - \frac{15}4 \sin^2\lambda + \frac38$.

    slsq = sin(lat) ** 2
    clsq = 1.0_rp - slsq
    p2 = (3.0_rp * slsq - 1.0_rp) * 0.5_rp
    p4 = (35.0_rp * slsq**2 - 30.0_rp * slsq + 3.0_rp) * 0.125_rp
!   dp2 = 3.0_rp * sl * cl
!   dp4 = 2.5_rp * sl * cl * ( 7.0_rp * slsq - 3.0_rp )

!{ Compute radius at geoid having potential = $W_0$ (see Geometry module):
!  $r_e^2 = \frac{a^2 b^2}{a^2 \sin^2 \lambda + b^2 \cos^2 \lambda}$ in meters.
!  First compute $r_e^{-2}$ in $m^{-2}$.  The parenthesization used here
!  avoids potential overflow; $a^2 b^2 \approx 1.64 \times 10^{27}$, but
!  the inner factor is $\approx 1$ and $a^2 \approx 4 \times 10^{13}$.

    r_eisq = ( (eRadAsq*slsq + eRadBsq*clsq) / eRadBsq ) / eRadAsq ! (r_e)^{-2}

    r_e = sqrt(1.0_rp / r_eisq) ! in meters

!{ Radial surface acceleration at latitude $\lambda$, k$m/s^2$:
!  \begin{equation*}
!   g_{\text{ref}} = 0.001 \left ( G m
!                     \frac{1 - 3 J_2 P_2(\sin \lambda)
!                       (a/r_e)^2 - 5 J_4 P_4(\sin \lambda) (a/r_e)^4 }{r_e^2}
!                       - \omega^2 \cos^2 \lambda \, r_e \right )
!  \end{equation*}

    g_ref = 0.001_rp * (gm * (1.0_rp - 3.0_rp*j2*p2*eRadAsq*r_eisq &
        & - 5.0_rp*j4*p4*(eRadAsq*r_eisq)**2)*r_eisq - w**2 * clsq * r_e)

!{ Better effective Earth radius:
!  compute $-2\, g_{\text{ref}} / ( d g_{\text{ref}} / d r)$, kilometers:
!  \begin{equation*}
!   r_{\text{eff}} =
!    \frac{2\, g_{\text{ref}}}
!         { 2\, G m \frac{1 - 6 J_2 P_2 (\sin \lambda) (a/r_e)^2
!                         - 15 J_4 P_4 (\sin \lambda) (a/r_e)^4}{r_e^3}
!                          + \omega^2 \cos^2 \lambda }
!  \end{equation*}
!  $d g_{\text{ref}} / d r$ has units $s^{-2}$, because we compute
!  $d g_{\text{ref}}(\text{meters}) / d r$(meters), thereby canceling
!  the factor of 0.001 you might have been expecting in the denominator.
!  The entire expression has units (k$ms^{-2})/s^{-2} = $ k$m$.

    r_eff = 2.0_rp * g_ref / (2.0_rp * gm * (1.0_rp-6.0_rp*j2*p2 &
        & * eRadAsq*r_eisq - 15.0_rp*j4*p4*(eRadAsq*r_eisq)**2) &
        & / r_e**3 + w**2 * clsq)

! find the surface pressure

    ghb = g_ref * h_ref / boltz
    z_old = z_grid(1) ! This is a guess

    iter = 0
    do
!{ Newton iteration for $z_{\text{surf}}$ at $h_{\text{calc}} = -h_{\text{ref}}$:
!  $z_{n+1} = z_n - (h_{\text{ref}} + h_{\text{calc}}) /
!             \left . \frac{\text{d} h}{\text{d} z} \right |_{z=z_n}$

      call piq_int ( (/z_old/), t_basis, z_ref, piqa )
      h_calc = dot_product(piqa(1,:), t_coeffs)

      call piq_int ( (/z_old+0.01_rp/), t_basis, z_ref, piqa )
      call piq_int ( (/z_old-0.01_rp/), t_basis, z_ref, piqb )
      dh_dz_s = dot_product((piqa(1,:) - piqb(1,:)), t_coeffs) * 50.0_rp
      z_surf = z_old - (ghb + h_calc) / dh_dz_s

      iter = iter + 1
      if ( abs(z_surf - z_old) < 0.0001_rp .or. iter == 10 ) exit
      z_old = z_surf

    end do

!{ Compute the {\tt piq} integrals with mass reduction compensation relative to the
!  surface.  Here, {\tt piq} = $P_l(\zeta) = P(\zeta_l,\phi_m)$.

    call piq_int ( z_grid, t_basis, z_surf, piq, Z_MASS=2.5_rp, C_MASS=0.02_rp )

! compute the height vector

! geopotential height * g_ref


!{ Equation (5.27) in the 19 August 2004 ATBD is
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{g_0 \stackrel{\star}{R_0^2}}
!         {g_0 \stackrel{\star}{R_0} - k\, \ln 10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!           ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!    -\stackrel{\star}{R_0} + R_0 - R^\oplus
!  \end{equation*}
!  The last two terms, i.e., $R_0-R^\oplus$, are handled in {\tt metrics}.
!  Putting what remains over a common denominator, we have
!  \begin{equation*}
!   h(\zeta,\phi) =
!    \frac{\stackrel{\star}{R_0} k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!         {g_0 \stackrel{\star}{R_0} - k \, \ln10
!          \sum_{l=1}^{\text{NH}^T} \sum_{m=1}^{\text{NP}^T}
!          ( f^T_{lm} \eta^T_m(\phi) P_l(\zeta))}
!  \end{equation*}
!  Notice that $g_0$ here is not the average or equatorial surface
!  acceleration.  Rather, it is $g_\text{ref}$, which is corrected for
!  latitude, the Earth's figure (up to fourth order zonal harmonics), and
!  centripetal acceleration.
    h_grid = boltz * matmul(piq,t_coeffs)
    h_grid = r_eff * h_grid / (r_eff * g_ref - h_grid)

  end subroutine Hydrostatic_No_Der

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Hydrostatic_m
!---------------------------------------------------
! $Log$
! Revision 2.32  2018/05/14 23:36:36  vsnyder
! Change to sparse eta representation
!
! Revision 2.31  2017/08/28 20:28:08  livesey
! Changed the n,nf,np,nz elements to j,jf,...
!
! Revision 2.30  2017/01/14 02:57:11  vsnyder
! Make Eta_ZZ polymorphic
!
! Revision 2.29  2016/12/02 02:04:50  vsnyder
! Use 'P' Eta list for Eta_ZZ
!
! Revision 2.28  2016/11/29 00:29:18  vsnyder
! Use interpolator in Indexed_Values_m
!
! Revision 2.27  2016/11/23 20:10:59  vsnyder
! Add Hydrostatic_All_ZZ, which uses an Eta list instead of computing it anew
!
! Revision 2.26  2016/04/28 18:04:12  vsnyder
! More TeXnicalities
!
! Revision 2.25  2015/04/11 00:45:03  vsnyder
! Add units (km) in h_grid comment
!
! Revision 2.24  2014/09/05 21:12:31  vsnyder
! Get kinds from MLSKinds instead of MLSCommon.  Add more tracing.
!
! Revision 2.23  2013/06/12 02:26:27  vsnyder
! Cruft removal
!
! Revision 2.22  2013/05/02 19:45:33  vsnyder
! Simplify r_eisq calculation, add a lot of LaTeX
!
! Revision 2.21  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.20  2007/01/17 23:50:15  vsnyder
! Exchange dhidzi, dhidtq, make dhidtq optional
!
! Revision 2.19  2006/11/30 23:28:29  vsnyder
! Correct some local array sizes
!
! Revision 2.18  2005/12/29 01:12:47  vsnyder
! Move some misplaced comments to correct places
!
! Revision 2.17  2005/12/22 20:57:10  vsnyder
! Added Hydrostatic_No_Der, some cannonball polishing
!
! Revision 2.16  2005/12/10 01:53:10  vsnyder
! Use get_eta_matrix_m instead of get_eta_matrix
!
! Revision 2.15  2005/12/07 00:32:21  vsnyder
! Cannonball polishing
!
! Revision 2.14  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.13  2004/03/20 04:06:13  vsnyder
! Moved Boltz from units to physics
!
! Revision 2.12  2003/09/16 18:43:47  vsnyder
! Might as well TeX-ify the comment about the approximation
!
! Revision 2.11  2003/09/16 18:31:02  vsnyder
! Put in a better comment about an approximation
!
! Revision 2.10  2003/09/16 00:21:32  vsnyder
! Remove unused arguments to get_eta
!
! Revision 2.9  2003/02/12 20:52:38  bill
! fixed serious bug in mass correction polynomial
!
! Revision 2.8  2003/02/10 23:41:27  bill
! got rid of a spread statement
!
! Revision 2.7  2003/02/08 00:59:16  bill
! uses latest piq int
!
! Revision 2.6  2003/02/07 00:04:54  bill
! mass function bug fix
!
! Revision 2.5  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/09/27 01:48:36  vsnyder
! Simplify iteration for z_surf
!
! Revision 2.3  2002/09/25 22:52:54  vsnyder
! Move USE from module scope to procedure scope.  Convert allocatable arrays
! to automatic arrays.  Replace sum(reshape(a...)*b) by dot_product.  Make
! constants consistently kind(rp) -- which caused a 1.5 epsilon relative
! change in output radiances.  Simplify Newton iteration.  Do some comments
! with LaTeX.
!
! Revision 2.2  2002/06/25 17:01:08  bill
! added more digits to J2 and added pressure dependent mass--wgr
!
! Revision 2.1  2002/02/02 11:20:08  zvi
! Some cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1.2.3  2001/09/13 22:51:22  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:50  zvi
! Added CVS stuff
!
