! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

module Hydrostatic_m

  implicit none

  private
  public :: Hydrostatic

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  contains
!---------------------------------------------------------------------------
  subroutine Hydrostatic ( lat, t_basis, t_coeffs, z_grid, z_ref, h_ref, &
                      &    t_grid, h_grid, dhidtq, dhidzi, ddhdhdtq, z_surface )

! Compute a hydrostatic function per L2PC method and return
! geometric heights. Reference height is now an input
! This is for EOS prototyping

    use MLSCommon, only: RP, IP
    use Geometry, only: EarthRadA, EarthRadB, GM, J2, J4, W
    use Get_eta_m, only: get_eta
    use Piq_int_m, only: piq_int
    use Units, only: BoltzMeters => Boltz

! Inputs

    real(rp), intent(in) :: lat ! geocentric latitude in radians
    real(rp), intent(in) :: z_ref ! reference pressure in -log10(P)
    real(rp), intent(in) :: h_ref ! reference geopotential height in km
    real(rp), intent(in) :: t_basis(:) ! vertical temperature basis
    real(rp), intent(in) :: t_coeffs(:) ! temperature values
    real(rp), intent(in) :: z_grid(:) ! -log10(P) pressures for which heights
!                                  are needed
! Outputs

    real(rp), intent(out) :: h_grid(:) ! heights on z_grid
    real(rp), intent(out) :: t_grid(:) ! temperatures on z_grid
    real(rp), intent(out) :: dhidzi(:) ! dh/dz on z_grid
    real(rp), intent(out) :: dhidtq(:,:) ! dh/dt on z_grid assuming dh/dt = 0.0
!                 at the reference ellipse surface  equivalent to h_ref = 0.0

    real(rp), optional, intent(out) :: ddhdhdtq(:,:) !ddh/dhdt on z_grid
    real(rp), optional :: z_surface

! Internal stuff

    real(rp), parameter :: ERadAsq = EarthRadA**2, ERadBsq = EarthRadB**2
    real(rp), parameter :: Boltz = boltzMeters/1.0e6_rp ! = kln10/m km^2/(K sec^2)

    integer(ip) :: n_lvls,n_coeffs,iter

!   real(rp) :: cl, sl ! for derivatives of Legendre polynomials dp2 and dp4
    real(rp) :: Clsq, G_ref, GHB
    real(rp) :: R_e, R_eisq, R_eisq_A, R_eff, Slsq, Z_surf
    real(rp) :: P2, P4    ! Legendre polynomials, for oblateness model
!   real(rp) :: dp2, dp4  ! Derivatives of P2 and P4
    real(rp) :: dh_dz_S, H_calc, Z_old
    real(rp), dimension(size(z_grid),size(t_basis)) :: Eta, Piq
    real(rp), dimension(1,size(t_basis)) :: Piqa, Piqb
    real(rp), dimension(size(t_basis)) :: Mass_corr
!
! begin the code
!
    n_lvls = size(z_grid)
    n_coeffs = size(t_basis)

    where ( t_basis > 2.5_rp )
      mass_corr = 1.0_rp / (0.875_rp + 0.1_rp*t_basis - 0.02_rp*t_basis**2)
    elsewhere
      mass_corr = 1.0_rp
    end where

! compute t_grid

    call get_eta ( z_grid, t_basis, n_lvls, n_coeffs, eta )
    t_grid = matmul(eta,t_coeffs)
!
! compute surface acceleration and effective earth radius
!
    slsq = sin(lat) ** 2
    clsq = 1.0_rp - slsq
    p2 = (3.0_rp * slsq - 1.0_rp) * 0.5_rp
    p4 = (35.0_rp * slsq**2 - 30.0_rp * slsq + 3.0_rp) * 0.125_rp
!   dp2 = 3.0_rp * sl * cl
!   dp4 = 2.5_rp * sl * cl * ( 7.0_rp * slsq - 3.0_rp )

!{ compute earth radius having potential = 62636858.0 $m/s^2$:
!  $r_e^2 = \frac{a^2 b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta}$

    r_eisq_a = (eRadAsq*slsq + eRadBsq*clsq) / eRadBsq ! m^{-1}

    r_eisq = r_eisq_a / eRadAsq ! (r_e)^{-2} in m^{-2}

    r_e = sqrt(1.0_rp / r_eisq) ! in meters

!{ radial surface acceleration, k$m/s^2$:\\
!  $g_{\text{ref}} = 0.001 \left ( G m ( 1 - 3 J_2 P_2(\sin \beta)
!                       (a/r_e)^2 - 5 J_4 P_4(\sin \beta) (a/r_e)^4 ) /r_e^2
!                       - \omega^2 \cos^2 \beta \, r_e \right )$

    g_ref = 0.001_rp * (gm * (1.0_rp - 3.0_rp*j2*p2*eRadAsq*r_eisq &
        & - 5.0_rp*j4*p4*(eRadAsq*r_eisq)**2)*r_eisq - w**2 * clsq * r_e)

!{ better effective Earth radius: compute $-d g_{\text{ref}} / d r$, kilometers:\\
!  $r_{\text{eff}} = 2 g_{\text{ref}} / \left ( 2 G m (
!                      ( 1 - 6 J_2 P_2 (\sin \beta) (a/r_e)^2
!                          - 15 J_4 P_4 (\sin \beta) (a/r_e)^4 ) / r_e^3
!                          + \omega^2 \cos^2 \beta \right )$

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

! compute the piq integrals relative to the surface

    call piq_int ( z_grid, t_basis, z_surf, piq )

! compensate piq for mass reduction. Note this is not the
! most rigourous code because it assumes the mass is constant
! across a coefficient.

    piq = piq * SPREAD(mass_corr,1,n_lvls)

! compute the height vector

! geopotential height * g_ref
    h_grid = boltz * matmul(piq,t_coeffs)
    h_grid = r_eff * h_grid / (r_eff * g_ref - h_grid)
    dhidzi = (h_grid+r_eff)**2 * boltz / (g_ref * r_eff**2)
    dhidtq = spread(dhidzi,2,n_coeffs) * piq
    dhidzi = dhidzi * t_grid

! this derivative is useful for antenna derivatives

    if ( present(ddhdhdtq) ) &
      & ddhdhdtq = (2.0_rp/(spread(h_grid,2,n_coeffs)+r_eff)) * dhidtq &
               & + eta / spread(t_grid,2,n_coeffs)

 end subroutine Hydrostatic

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Hydrostatic_m
!---------------------------------------------------
! $Log$
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
