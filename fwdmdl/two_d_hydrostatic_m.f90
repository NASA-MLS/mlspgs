! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

module Two_D_Hydrostatic_m

  implicit none

  private
  public Two_D_Hydrostatic

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
  contains
!---------------------------------------------------------------------------

  subroutine Two_D_Hydrostatic ( Grids_tmp, z_refs, h_refs, z_grid, beta, &
                              &  t_grid, h_grid, dhidzij, dhidtlm, ddhdhdtl0 )

! Compute the 2 dimensional hydrostatic stuff

  use Geometry, only: earthRadA, earthRadB
  use Hydrostatic_m, only: hydrostatic
  use Load_sps_data_m, ONLY: Grids_T
  use MLSCommon, only: RP, IP
  use Units, only: PI

! Inputs:

  type (Grids_T), intent(in) :: Grids_tmp   ! All Temperature's coordinates
  real(rp), intent(in) :: z_refs(:)  !reference pressures
  real(rp), intent(in) :: h_refs(:)  !reference geopotential heights at z_refs
!                                 the horizontal bases for these is aligned
!                                 with p_basis.
  real(rp), intent(in) :: z_grid(:)!pressures for which heights/temps are
!                                   needed
  real(rp), intent(in) :: beta   ! spacecraft beta angle (Radians)

! Outputs:

  real(rp), intent(out):: t_grid(:,:)!computed temperatures
  real(rp), intent(out):: h_grid(:,:)!computed heights
  real(rp), intent(out):: dhidzij(:,:)!derivative of height wrt zeta
  real(rp), intent(out):: dhidtlm(:,:,:) !derivative of height wrt temps
!                                     on outputted phi grid
  real(rp), optional, intent(out):: ddhdhdtl0(:,:,:)!second order derivative
!                             at the tangent only---used for antenna affects
! Internal stuff

  integer(ip) :: Z_coeffs, P_coeffs, I, J1, J2

  real(rp) :: CSQ ! C**2
  real(rp) :: Lat, Red_phi_t, SinBeta, SinBetaSQ, SinPhi, SinPhiSQ

  real(rp), parameter :: PId2=0.5_rp*pi, PI2=2.0_rp*pi, PI3d2=1.5_rp*pi

! Begin execution

  z_coeffs = Grids_tmp%l_z(1) ! - Grids_tmp%l_z(0), which is always zero
  p_coeffs = Grids_tmp%l_p(1) ! - Grids_tmp%l_p(0), which is always zero

!{ Compute the orbit-plane projected minor axis $c$, where
!  $c^2 = \frac{a^2\,b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta}$

  sinBeta = sin(beta)
  sinBetaSQ = sinBeta**2
  csq = (earthrada*earthradb)**2 / (earthrada**2 * sinBetaSQ + &
      &                             earthradb**2 * (1.0_rp - sinBetaSQ) ) ! in meters

! compute the 2 d hydrostatic

  j2 = 0
  do i = 1, p_coeffs

! rephase the phi

    red_phi_t = modulo(Grids_tmp%phi_basis(i),PI2)
    if ( PiD2 < red_phi_t .and. red_phi_t <= Pi3D2 ) then
      red_phi_t = Pi - red_phi_t
    else if ( red_phi_t > Pi3D2 ) then
      red_phi_t = red_phi_t - Pi2
    end if

!{ Compute the geocentric latitude $\lambda\, = \, 
!  \sin^{-1} ( \sin \gamma \, \sin \beta )$, where
!  $\sin^2 \gamma = \frac{c^4 \sin^2 \phi}{a^4 \cos^2 \phi + c^4 \sin^2 \phi}$
!  and $\gamma$ is the geocentric angle in the orbit plane ellipse between
!  $\mathbf{R}^{\oplus}$ and the $x$ axis.

    sinPhi = sin(red_phi_t)
    sinPhiSQ = sinPhi**2
    lat = asin(csq * sinPhi * sinBeta &
      & / sqrt(earthrada**4*(1.0_rp-sinPhiSQ) + csq**2*sinPhiSQ))

    j1 = j2
    j2 = j1 + z_coeffs
    if ( present(ddhdhdtl0) ) then
      call hydrostatic ( lat, Grids_tmp%zet_basis, Grids_tmp%values(j1+1:j2), &
         & z_grid, z_refs(i), h_refs(i), t_grid(:,i), h_grid(:,i), &
         & dhidtlm(:,:,i), dhidzij(:,i), ddhdhdtl0(:,:,i) )
    else
      call hydrostatic ( lat, Grids_tmp%zet_basis, Grids_tmp%values(j1+1:j2), &
         & z_grid, z_refs(i), h_refs(i), t_grid(:,i), h_grid(:,i), &
         & dhidtlm(:,:,i), dhidzij(:,i) )
    end if
  end do

  end subroutine Two_D_Hydrostatic

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module Two_D_Hydrostatic_m
!---------------------------------------------------
! $Log$
! Revision 2.10.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.10  2002/10/10 19:52:45  vsnyder
! Get rid of several array temps.  Cosmetic changes.
!
! Revision 2.9  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/09/26 20:14:01  vsnyder
! Get PI from Units module
!
! Revision 2.7  2002/09/25 22:55:12  vsnyder
! Move USE statements from module scope to procedure scope.  Convert
! allocatable arrays to automatic arrays.  Cosmetic changes.
!
! Revision 2.6  2002/07/05 07:52:53  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.5  2002/06/24 21:11:25  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 2.2  2002/06/07 14:59:00  bill
! fixed latitude calculation--wgr
!
! Revision 2.1  2002/02/02 11:20:08  zvi
! Some cosmetic changes
!
! Revision 2.0  2001/09/17 20:26:28  livesey
! New forward model
!
! Revision 1.1.2.3  2001/09/13 22:51:25  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:55  zvi
! Added CVS stuff
