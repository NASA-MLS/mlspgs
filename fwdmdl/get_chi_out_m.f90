! This module computes the output angles to interpolate to
MODULE get_chi_out_m
  USE MLSCommon, only: I4, RP, IP
  USE Allocate_Deallocate, only: allocate_test, deallocate_test
  USE Load_sps_data_m, ONLY: Grids_T
  USE get_eta_matrix_m, only: get_eta_sparse
  USE get_chi_angles_m, only: get_chi_angles
  USE two_d_hydrostatic_m, only: two_d_hydrostatic
  USE refraction_m, only: refractive_index
  IMPLICIT none
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------------
 CONTAINS
!
  SUBROUTINE get_chi_out(zetatan,phitan,scgeocalt,Grids_tmp,ref_zeta, &
           & ref_gph,orb_inc,elev_offset,req,tan_chi_out,h2o_zeta_basis,&
           & h2o_phi_basis,h2o_coeffs,lin_log,dxdt_tan,d2xdxdt_tan)
!
! inputs
!
  REAL(rp), INTENT(in) :: zetatan(:)   ! tangent zeta profle for input mmaf
  REAL(rp), INTENT(in) :: phitan(:)    ! tangent phi profle for input mmaf
!                                        (in radians)
  REAL(rp), INTENT(in) :: scgeocalt(:) ! spacecraft geocentric altitude profile in km
  type (Grids_T), INTENT(in) :: Grids_tmp  ! All Temperature's coordinates

  REAL(rp), INTENT(in) :: ref_zeta(:) ! zetas for inputted gph's
  REAL(rp), INTENT(in) :: ref_gph(:) ! reference geopotential heights along the
!                          temperature horizontal basis in km
  REAL(rp), INTENT(in) :: orb_inc    ! orbital incline angle in radians
  REAL(rp), INTENT(in) :: elev_offset ! elevation offset in radians
  REAL(rp), INTENT(in) :: req  ! Earth radius in equivalent circle
!                                representation (km)
! Keywords
  REAL(rp), OPTIONAL, INTENT(in) :: h2o_zeta_basis(:) ! h2o zeta basis
  REAL(rp), OPTIONAL, INTENT(in) :: h2o_phi_basis(:) ! h2o phi basis
  REAL(rp), OPTIONAL, INTENT(in) :: h2o_coeffs(:) ! h2o coefficients in long
!                          vector format zetaXphi.
  LOGICAL,  OPTIONAL, INTENT(in) :: lin_log ! representation type
!
! outputs
!
  REAL(rp), INTENT(out) :: tan_chi_out(:) ! computed tangent pointing angles
!                               corrected for refraction in radians
  REAL(rp), OPTIONAL, INTENT(out) :: dxdt_tan(:,:,:) ! computed dchi dt.
  REAL(rp), OPTIONAL, INTENT(out) :: d2xdxdt_tan(:,:,:) ! computed d2chi dxdt.
!
! the matrix hierarchy is n_out,phi,zeta
! internal variables
!
  INTEGER(i4) :: n_out, n_t_phi, n_t_zeta, n_h2o_zeta, n_h2o_phi
  INTEGER(i4) :: sv_p, sv_z, ht_i
!
  REAL(rp), POINTER :: h_tan_out(:)
  REAL(rp), POINTER :: t_tan_out(:)
  REAL(rp), POINTER :: n_tan_out(:)
  REAL(rp), POINTER :: h2o_tan_out(:)
  REAL(rp), POINTER :: temp_tan(:,:)
  REAL(rp), POINTER :: height_tan(:,:)
  REAL(rp), POINTER :: dhdz_tan(:,:)
  REAL(rp), POINTER :: dhdt_tan(:,:,:)
  REAL(rp), POINTER :: d2hdhdt_tan(:,:,:)
  REAL(rp), POINTER :: eta_t(:,:)
  REAL(rp), POINTER :: eta_z(:,:)
  REAL(rp), POINTER :: eta_p(:,:)
!
! dimensions we use
!
  n_out = SIZE(zetatan)
  n_t_phi = Grids_tmp%no_p(1)
  n_t_zeta = Grids_tmp%no_z(1)
!
  NULLIFY(temp_tan,height_tan,dhdz_tan,dhdt_tan,d2hdhdt_tan,eta_t,h_tan_out, &
  & t_tan_out,n_tan_out,eta_p,eta_z,h2o_tan_out)
!
  CALL allocate_test(temp_tan,n_out,n_t_phi,'temp_tan',ModuleName)
  CALL allocate_test(height_tan,n_out,n_t_phi,'height_tan',ModuleName)
  CALL allocate_test(dhdz_tan,n_out,n_t_phi,'dhdz_tan',ModuleName)
  CALL allocate_test(dhdt_tan,n_out,n_t_phi,n_t_zeta,'dhdt_tan',ModuleName)
  CALL allocate_test(d2hdhdt_tan,n_out,n_t_phi,n_t_zeta,'d2hdhdt_tan', &
                   & ModuleName)
  CALL allocate_test(eta_t,n_out,n_t_phi,'eta_t',ModuleName)
  CALL allocate_test(h_tan_out,n_out,'h_tan_out',ModuleName)
  CALL allocate_test(t_tan_out,n_out,'t_tan_out',ModuleName)
  CALL allocate_test(n_tan_out,n_out,'n_tan_out',ModuleName)
!
  CALL two_d_hydrostatic(Grids_tmp,ref_zeta,ref_gph,zetatan,orb_inc, &
           & temp_tan, height_tan, dhdz_tan, dhdt_tan,d2hdhdt_tan)
!
! tangent heights for inputted pressures along phi
!
  CALL get_eta_sparse(Grids_tmp%phi_basis, phitan, eta_t)
  h_tan_out = SUM(height_tan * eta_t, dim=2)
  t_tan_out = SUM(temp_tan * eta_t, dim=2)
!
! compute the tangent water vapor profile along inputted phi_tan
!
  IF (PRESENT(h2o_coeffs)) THEN
    n_h2o_zeta = SIZE(h2o_zeta_basis)
    n_h2o_phi  = SIZE(h2o_phi_basis)
    CALL allocate_test(eta_p, n_out, n_h2o_phi, 'eta_p',ModuleName)
    CALL allocate_test(eta_z, n_out, n_h2o_zeta, 'eta_z',ModuleName)
    CALL allocate_test(h2o_tan_out, n_out, 'h2o_tan_out', &
    & ModuleName)
    CALL get_eta_sparse(h2o_phi_basis, phitan, eta_p)
    CALL get_eta_sparse(h2o_zeta_basis, zetatan, eta_z)
!
! This is actually one less than the true start index
!
    h2o_tan_out(:) = 0.0_rp
    DO sv_p = 1, n_h2o_phi
      DO sv_z = 1, n_h2o_zeta
        h2o_tan_out = h2o_tan_out + h2o_coeffs(sv_z - n_h2o_zeta &
        & + n_h2o_zeta * sv_p) * eta_z(:,sv_z) * eta_p(:,sv_p)
      ENDDO
    ENDDO
    IF (lin_log) h2o_tan_out = EXP(h2o_tan_out)
!
! compute refractive index
!
    CALL refractive_index(10.0**(-zetatan), t_tan_out, n_tan_out, &
    & h2o_path = h2o_tan_out)
    CALL deallocate_test(eta_p,'eta_p',ModuleName)
    CALL deallocate_test(eta_z,'eta_z',ModuleName)
    CALL deallocate_test(h2o_tan_out,'h2o_tan_out',ModuleName)
!
  ELSE
!
    CALL refractive_index(10.0**(-zetatan), t_tan_out, n_tan_out)
!
  ENDIF
!
! compute output angles to interpolate to
!
  IF(PRESENT(dxdt_tan)) THEN
!
    dhdt_tan = dhdt_tan  * SPREAD(eta_t,3,n_t_zeta)
    d2hdhdt_tan = d2hdhdt_tan * SPREAD(eta_t,3,n_t_zeta)
    DO ht_i = 1, n_out
      CALL get_chi_angles(scgeocalt(ht_i),n_tan_out(ht_i), &
          &    h_tan_out(ht_i),phitan(ht_i),Req,elev_offset, &
          &    tan_chi_out(ht_i), RESHAPE(dhdt_tan(ht_i,:,:), &
          &    (/n_t_phi,n_t_zeta/)), RESHAPE(d2hdhdt_tan(ht_i,:,:), &
          &    (/n_t_phi,n_t_zeta/)), dxdt_tan(ht_i,:,:), &
          &    d2xdxdt_tan(ht_i,:,:))
     ENDDO
!
   ELSE
!
    DO ht_i = 1, n_out
      CALL get_chi_angles(scgeocalt(ht_i),n_tan_out(ht_i), &
          &    h_tan_out(ht_i),phitan(ht_i),Req,elev_offset, &
          &    tan_chi_out(ht_i))
     ENDDO
!
  ENDIF
!
  CALL deallocate_test (temp_tan,'temp_tan',ModuleName)
  CALL deallocate_test (height_tan,'height_tan',ModuleName)
  CALL deallocate_test (dhdz_tan,'dhdz_tan',ModuleName)
  CALL deallocate_test (dhdt_tan,'dhdt_tan',ModuleName)
  CALL deallocate_test (d2hdhdt_tan,'dhdt_tan',ModuleName)
  CALL deallocate_test (eta_t,'eta_t',ModuleName)
  CALL deallocate_test (h_tan_out,'h_tan_out',ModuleName)
  CALL deallocate_test (t_tan_out,'t_tan_out',ModuleName)
  CALL deallocate_test (n_tan_out,'n_tan_out',ModuleName)
!
 END SUBROUTINE get_chi_out
!
END MODULE get_chi_out_m
! $Log$
! Revision 2.5  2002/06/24 21:07:13  zvi
! Bug fixing and correcting Log entries
!
! Revision 2.4  2002/06/24 21:01:28  zvi
! Adding Grids_tmp stracture and modify calling sequences
!
! Revision 1.0  2002/06/24 14:59:00  bill
! First release ...
