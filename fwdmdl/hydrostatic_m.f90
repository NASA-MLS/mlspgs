MODULE hydrostatic_m
!
! This computes a hydrostatic function per L2PC method
! and returns geometric heights
! reference height is now an input
! This is for EOS prototyping
!
  use MLSCommon, only: RP, IP
  USE Geometry, ONLY: earthrada,earthradb
  USE get_eta_m, only: get_eta
  USE piq_int_m, only: piq_int
  IMPLICIT NONE
!
  Private
  Public :: hydrostatic
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
   "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
   "$RCSfile$"
!---------------------------------------------------------------------------
  CONTAINS
!---------------------------------------------------------------------------
  SUBROUTINE hydrostatic(lat,t_basis,t_coeffs,z_grid,z_ref,h_ref,t_grid, &
                  h_grid,dhidtq,dhidzi,ddhdhdtq,z_surface)
!
! Inputs
!
  REAL(rp), INTENT(in) :: lat ! geocentric latitude in radians
  REAL(rp), INTENT(in) :: z_ref ! reference pressure in -log10(P)
  REAL(rp), INTENT(in) :: h_ref ! reference geopotential height in km
  REAL(rp), INTENT(in) :: t_basis(:) ! vertical temperature basis
  REAL(rp), INTENT(in) :: t_coeffs(:) ! temperature values
  REAL(rp), INTENT(in) :: z_grid(:) ! -log10(P) pressures for which heights
!                                  are needed
! Outputs
!
  REAL(rp), INTENT(out) :: h_grid(:) ! heights on z_grid
  REAL(rp), INTENT(out) :: t_grid(:) ! temperatures on z_grid
  REAL(rp), INTENT(out) :: dhidzi(:) ! dh/dz on z_grid
  REAL(rp), INTENT(out) :: dhidtq(:,:) ! dh/dt on z_grid assuming dh/dt = 0.0
!                 at the reference ellipse surface  equivalent to h_ref = 0.0
  REAL(rp), OPTIONAL, INTENT(out) :: ddhdhdtq(:,:) !ddh/dhdt on z_grid
  REAL(rp), OPTIONAL :: z_surface
!
! Internal stuff
! These are the 1980 reference geoid values these should be later
! moved to a separate module (preferably Geometry.f90)
!
  REAL(rp), PARAMETER :: gm = 3.986005e14_rp ! m^3/sec^2
  REAL(rp), PARAMETER :: j2 = 0.00108263_rp
  REAL(rp), PARAMETER :: j4 = -.0000023709122_rp
!
! earth rotational velocity
!
  REAL(rp), PARAMETER :: w = 7.292115e-05_rp ! rad/sec
!
  INTEGER(ip) :: n_lvls,n_coeffs,iter
!
  REAL(rp) :: sl,cl,p2,p4,dp2,dp4,g_ref,r_eff,boltz,r_e,z_surf
  REAL(rp) :: z_new,h_calc,dh_dz_s
  REAL(rp), ALLOCATABLE :: eta(:,:),piq(:,:),piqa(:,:),piqb(:,:)
!
! begin the code
!
  n_lvls = SIZE(z_grid)
  n_coeffs = SIZE(t_basis)
  ALLOCATE(eta(1:n_lvls,1:n_coeffs))
  ALLOCATE(piq(1:n_lvls,1:n_coeffs))
!
! compute t_grid
!
  CALL get_eta(z_grid,t_basis,n_lvls,n_coeffs,eta)
  t_grid = MATMUL(eta,t_coeffs)
!
! compute surface acceleration and effective earth radius
!
  sl = SIN(lat)
  cl = COS(lat)
  p2 = (3.0_rp * sl**2 - 1.0_rp) / 2.0_rp
  p4 = (35.0_rp * sl**4 - 30.0_rp * sl**2 + 3.0_rp) / 8.0_rp
  dp2 = 3.0_rp * sl * cl
  dp4 = 5.0_rp * (7.0_rp * sl**3 * cl - 3.0_rp * sl * cl) / 2.0_rp
!
! compute earth radius having potential = 62636858.0 m/sec^2
!
  r_e = earthrada*earthradb &
      / SQRT((earthrada*sl)**2+(earthradb*cl)**2) !in meters
!
! radial surface acceleration, km/sec^2
!
  g_ref = 0.001 * (gm * (1.0_rp - 3.0_rp*j2*p2*(earthrada/r_e)**2 &
        - 5.0_rp*j4*p4*(earthrada/r_e)**4)/r_e**2 - (w*cl)**2 * r_e)
!
! better effective Earth radius compute -d g_ref / d r, kilometers
!
  r_eff = 2.0_rp * g_ref / (2.0_rp * gm * (1.0_rp-6.0_rp*j2*p2 &
        * (earthrada/r_e)**2 - 15.0_rp*j4*p4*(earthrada/r_e)**4) &
        / r_e**3 + (w*cl)**2)
!
! find the surface pressure
!
  boltz = 0.000660988_rp ! = kln10/m km^2/(K sec^2)
  z_surf = z_grid(1) ! This is a guess

  ALLOCATE(piqa(1:1,1:n_coeffs))
  ALLOCATE(piqb(1:1,1:n_coeffs))

  CALL piq_int((/z_surf/),t_basis,z_ref,piqa)
  h_calc = boltz*SUM(RESHAPE(piqa,(/n_coeffs/)) * t_coeffs) / g_ref
!
! correct
!
  CALL piq_int((/z_surf+0.01_rp/),t_basis,z_ref,piqa)
  CALL piq_int((/z_surf-0.01_rp/),t_basis,z_ref,piqb)
  dh_dz_s = boltz*SUM((RESHAPE(piqa,(/n_coeffs/))-RESHAPE(piqb,(/n_coeffs/))) &
          * t_coeffs) / (0.02_rp*g_ref)
  z_new = z_surf - (h_ref + h_calc) / dh_dz_s
  iter = 1
  DO
    IF(ABS(z_new - z_surf) < 0.0001_rp .OR. iter == 10) EXIT
    z_surf = z_new
    CALL piq_int((/z_surf/),t_basis,z_ref,piqa)
    h_calc = boltz*SUM(RESHAPE(piqa,(/n_coeffs/)) * t_coeffs) / g_ref
!
! correct
!
    CALL piq_int((/z_surf+0.01_rp/),t_basis,z_ref,piqa)
    CALL piq_int((/z_surf-0.01_rp/),t_basis,z_ref,piqb)
    dh_dz_s = boltz*SUM((RESHAPE(piqa,(/n_coeffs/)) &
            - RESHAPE(piqb,(/n_coeffs/))) * t_coeffs) / (0.02_rp*g_ref)
    z_new = z_surf - (h_ref + h_calc) / dh_dz_s
    iter = iter + 1
  END DO
!
  DEALLOCATE(piqb)
  DEALLOCATE(piqa)

  z_surf = z_new
  IF (PRESENT(z_surface)) z_surface = z_surf
!
! compute the piq integrals relative to the surface
!
  CALL piq_int(z_grid,t_basis,z_surf,piq)
!
! compute the height vector
!
  h_grid = boltz * MATMUL(piq,t_coeffs) ! geopotential height * g_ref
  h_grid = r_eff * h_grid / (r_eff * g_ref - h_grid)
  dhidzi = (h_grid+r_eff)**2 * boltz / (g_ref * r_eff**2)
  dhidtq = SPREAD(dhidzi,2,n_coeffs) * piq
  dhidzi = dhidzi * t_grid
!
! this derivative is useful for antenna derivatives
!
  IF(PRESENT(ddhdhdtq)) &
  ddhdhdtq = (2.0_rp/(SPREAD(h_grid,2,n_coeffs)+r_eff)) * dhidtq &
           + eta / SPREAD(t_grid,2,n_coeffs)

  DEALLOCATE(piq)
  DEALLOCATE(eta)
!
 END SUBROUTINE hydrostatic
!
END MODULE hydrostatic_m
!---------------------------------------------------
! $Log$
! Revision 1.1.2.3  2001/09/13 22:51:22  zvi
! Separating allocation stmts
!
! Revision 1.1.2.2  2001/09/12 21:38:50  zvi
! Added CVS stuff
!
