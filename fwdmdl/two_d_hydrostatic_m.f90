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
!---------------------------------------------------------------------------
  contains
!---------------------------------------------------------------------------

  subroutine Two_D_Hydrostatic ( Grids_tmp, z_refs, h_refs, z_grid, beta, &
                              &  t_grid, h_grid, dhidzij, dhidtlm, ddhdhdtl0 )

! Compute the 2 dimensional hydrostatic stuff

  use Geometry, only: earthRadA, earthRadB
  use Get_eta_m, only: get_eta
  use Hydrostatic_m, only: hydrostatic
  use Load_sps_data_m, ONLY: Grids_T
  use MLSCommon, only: RP, IP
  use Units, only: PI

! inputs:

  type (Grids_T), intent(in) :: Grids_tmp   ! All Temperature's coordinates
  real(rp), intent(in) :: z_refs(:)  !reference pressures
  real(rp), intent(in) :: h_refs(:)  !reference geopotential heights at z_refs
!                                 the horizontal bases for these is aligned
!                                 with p_basis.
  real(rp), intent(in) :: z_grid(:)!pressures for which heights/temps are
!                                   needed
  real(rp), intent(in) :: beta   ! spacecraft beta angle (Radians)

! outputs:

  real(rp), intent(out):: t_grid(:,:)!computed temperatures
  real(rp), intent(out):: h_grid(:,:)!computed heights
  real(rp), intent(out):: dhidzij(:,:)!derivative of height wrt zeta
  real(rp), intent(out):: dhidtlm(:,:,:) !derivative of height wrt temps
!                                     on outputted phi grid
  real(rp), optional, intent(out):: ddhdhdtl0(:,:,:)!second order derivative
!                             at the tangent only---used for antenna affects
! internal stuff

  integer(ip) :: z_coeffs, p_coeffs, i, j1, j2

  real(rp) :: CSQ ! C**2
  real(rp), dimension(size(z_grid),Grids_tmp%no_z(1)) :: ddhdhdtq, dhidtq
  real(rp), dimension(Grids_tmp%no_p(1)) :: lats, red_phi_t
  real(rp), dimension(size(z_grid)) :: t_prfl, h_prfl, dhidzi
  real(rp) :: SinBeta, SinBetaSQ, SinPhi, SinPhiSQ

  real(rp), parameter :: PId2=0.5_rp*pi, PI2=2.0_rp*pi, PI3d2=1.5_rp*pi

! NOTES
! allocate arrays

  z_coeffs = Grids_tmp%no_z(1)
  p_coeffs = Grids_tmp%no_p(1)

  sinBeta = sin(beta)
  sinBetaSQ = sinBeta**2
  csq = (earthrada*earthradb)**2 / (earthrada**2 * sinBetaSQ + &
      &                             earthradb**2 * (1 - sinBetaSQ) ) ! in meters

! rephase the phi

  red_phi_t = modulo(Grids_tmp%phi_basis,PI2)
  where ( PiD2 < red_phi_t .AND. red_phi_t <= Pi3D2 )
    red_phi_t = Pi - red_phi_t
  elsewhere ( red_phi_t > Pi3D2 )
    red_phi_t = red_phi_t - Pi2
  end where

  do i = 1, size(lats) ! so we don't calculate both sin(red_phi_t) and cos(")
    sinPhi = sin(red_phi_t(i))
    sinPhiSQ = sinPhi**2
    lats(i) = asin(csq * sinPhi * sinBeta &
       / sqrt(earthrada**4*(1-sinPhiSQ) + csq**2*sinPhiSQ))
  end do

! compute the 2 d hydrostatic

  if ( present(ddhdhdtl0) ) then
    j2 = 0
    do i = 1, p_coeffs
      j1 = j2 + 1
      j2 = j1 + z_coeffs - 1
      call hydrostatic ( lats(i), Grids_tmp%zet_basis, Grids_tmp%values(j1:j2), &
         & z_grid, z_refs(i), h_refs(i), t_prfl, h_prfl, dhidtq, dhidzi, ddhdhdtq )
      t_grid(:,i) = t_prfl
      h_grid(:,i) = h_prfl
      dhidzij(:,i) = dhidzi
      dhidtlm(:,:,i) = dhidtq
      ddhdhdtl0(:,:,i) = ddhdhdtq
    end do
  else
    j2 = 0
    do i = 1,p_coeffs
      j1 = j2 + 1
      j2 = j1 + z_coeffs - 1
      call hydrostatic ( lats(i), Grids_tmp%zet_basis, Grids_tmp%values(j1:j2), &
         & z_grid, z_refs(i), h_refs(i), t_prfl, h_prfl, dhidtq, dhidzi )
      t_grid(:,i) = t_prfl
      h_grid(:,i) = h_prfl
      dhidzij(:,i) = dhidzi
      dhidtlm(:,:,i) = dhidtq
    end do
  end if

 end subroutine Two_D_Hydrostatic

end module Two_D_Hydrostatic_m
!---------------------------------------------------
! $Log$
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
