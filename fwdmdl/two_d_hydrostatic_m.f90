! This subroutine computes the 2 dimensional hydrostatic stuff
MODULE two_d_hydrostatic_m
  use MLSCommon, only: RP, IP
  USE Geometry, ONLY: earthRadA, earthRadB, PI
  USE hydrostatic_m, only: hydrostatic
  USE get_eta_m, only: get_eta
  USE Load_sps_data_m, ONLY: Grids_T
  IMPLICIT none
  Private
  Public two_d_hydrostatic
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
  CONTAINS
!---------------------------------------------------------------------------

  SUBROUTINE two_d_hydrostatic(Grids_tmp,z_refs,h_refs,z_grid,beta,t_grid, &
                            &  h_grid,dhidzij,dhidtlm,ddhdhdtl0)
!
! inputs:
!
  type (Grids_T), INTENT(in) :: Grids_tmp   ! All Temperature's coordinates
  REAL(rp), INTENT(in) :: z_refs(:)  !reference pressures
  REAL(rp), INTENT(in) :: h_refs(:)  !reference geopotential heights at z_refs
!                                 the horizontal bases for these is aligned
!                                 with p_basis.
  REAL(rp), INTENT(in) :: z_grid(:)!pressures for which heights/temps are
!                                   needed
  REAL(rp), INTENT(in) :: beta   ! spacecraft beta angle (Radians)
!
! outputs:
!
  REAL(rp), INTENT(out):: t_grid(:,:)!computed temperatures
  REAL(rp), INTENT(out):: h_grid(:,:)!computed heights
  REAL(rp), INTENT(out):: dhidzij(:,:)!derivative of height wrt zeta
  REAL(rp), INTENT(out):: dhidtlm(:,:,:) !derivative of height wrt temps
!                                     on outputted phi grid
  REAL(rp), OPTIONAL, INTENT(out):: ddhdhdtl0(:,:,:)!second order derivative
!                             at the tangent only---used for antenna affects
! internal stuff
!
  INTEGER(ip) :: n_vert,z_coeffs,p_coeffs,i,j1,j2

  REAL(rp) :: c
  REAL(rp), ALLOCATABLE, DIMENSION(:) :: lats,t_prfl,h_prfl,dhidzi,red_phi_t
  REAL(rp), ALLOCATABLE, DIMENSION(:,:) :: ddhdhdtq,dhidtq
  REAL(rp), PARAMETER :: Pi = 3.1415926535897932384626434_rp
!
! NOTES
! allocate arrays
!
  n_vert = SIZE(z_grid)
  z_coeffs = Grids_tmp%no_z(1)
  p_coeffs = Grids_tmp%no_p(1)
!
! allocate arrays
!
  ALLOCATE(lats(1:p_coeffs))
  ALLOCATE(t_prfl(1:n_vert))
  ALLOCATE(h_prfl(1:n_vert))
  ALLOCATE(dhidzi(1:n_vert))
  ALLOCATE(red_phi_t(1:p_coeffs))
  ALLOCATE(dhidtq(1:n_vert,1:z_coeffs))

  c = earthrada*earthradb / SQRT(earthrada**2 &
    * SIN(beta)**2 + earthradb**2*COS(beta)**2) ! in meters
!
! rephase the phi
!
  red_phi_t = MODULO(Grids_tmp%phi_basis,2.0_rp*Pi)
  WHERE(0.5_rp*Pi < red_phi_t .AND. red_phi_t <= 1.5_rp*Pi) &
              &  red_phi_t = Pi - red_phi_t
  WHERE(red_phi_t > 1.5_rp*Pi) red_phi_t = red_phi_t - 2.0_rp*Pi
  lats = ASIN(c**2 * SIN(red_phi_t) * SIN(beta) &
       / SQRT(earthrada**4*COS(red_phi_t)**2 + &
              c**4*SIN(red_phi_t)**2))
!
! compute the 2 d hydrostatic
!
  IF(PRESENT(ddhdhdtl0)) THEN
    ALLOCATE(ddhdhdtq(1:n_vert,1:z_coeffs))
    j2 = 0
    DO i = 1,p_coeffs
      j1 = j2 + 1
      j2 = j1 + z_coeffs - 1
      CALL hydrostatic(lats(i),Grids_tmp%zet_basis,Grids_tmp%values(j1:j2),&
         & z_grid,z_refs(i),h_refs(i),t_prfl,h_prfl,dhidtq,dhidzi,ddhdhdtq)
      t_grid(:,i) = t_prfl
      h_grid(:,i) = h_prfl
      dhidzij(:,i) = dhidzi
      dhidtlm(:,i,:) = dhidtq
      ddhdhdtl0(:,i,:) = ddhdhdtq
    END DO
    DEALLOCATE(ddhdhdtq)
  ELSE
    j2 = 0
    DO i = 1,p_coeffs
      j1 = j2 + 1
      j2 = j1 + z_coeffs - 1
      CALL hydrostatic(lats(i),Grids_tmp%zet_basis,Grids_tmp%values(j1:j2),&
         & z_grid,z_refs(i),h_refs(i),t_prfl,h_prfl,dhidtq,dhidzi)
      t_grid(:,i) = t_prfl
      h_grid(:,i) = h_prfl
      dhidzij(:,i) = dhidzi
      dhidtlm(:,i,:) = dhidtq
    END DO
  ENDIF
!
  DEALLOCATE(red_phi_t)
  DEALLOCATE(dhidzi)
  DEALLOCATE(dhidtq)
  DEALLOCATE(h_prfl)
  DEALLOCATE(t_prfl)
  DEALLOCATE(lats)
!
 END SUBROUTINE two_d_hydrostatic

END MODULE two_d_hydrostatic_m
!---------------------------------------------------
! $Log$
! Revision 2.3  2002/06/24 21:01:28  zvi
! Adding Grids_tmp stracture and modify calling sequences
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
