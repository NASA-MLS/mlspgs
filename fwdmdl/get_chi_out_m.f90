! This module computes the output angles to interpolate to
module Get_Chi_Out_m

  implicit NONE
  private
  public :: Get_Chi_Out

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter :: ModuleName = &
    & "$RCSfile$"
!---------------------------------------------------------------------------
contains

  subroutine Get_Chi_Out ( zetatan, phitan, scgeocalt, Grids_tmp, ref_zeta, &
           & ref_gph, orb_inc, elev_offset, req, tan_chi_out, dhdz_out, &
           & dx_dh_out, h2o_zeta_basis, h2o_phi_basis, h2o_coeffs, lin_log, &
           & dxdt_tan, d2xdxdt_tan )

    use Allocate_Deallocate, only: Allocate_test, Deallocate_test
    use Get_chi_angles_m, only: Get_chi_angles
    use Get_eta_matrix_m, only: Get_eta_sparse
    use Load_sps_data_m, only: Grids_T
    use MLSCommon, only: I4, RP
    use Refraction_m, only: Refractive_index
    use Two_d_hydrostatic_m, only: Two_d_hydrostatic

! inputs

    real(rp), intent(in) :: zetatan(:)   ! tangent zeta profle for input mmaf
    real(rp), intent(in) :: phitan(:)    ! tangent phi profle for input mmaf
!                                          (in radians)
    real(rp), intent(in) :: scgeocalt(:) ! spacecraft geocentric altitude profile in km
    type (grids_t), intent(in) :: Grids_tmp  ! All Temperature's coordinates

    real(rp), intent(in) :: ref_zeta(:) ! zetas for inputted gph's
    real(rp), intent(in) :: ref_gph(:) ! reference geopotential heights along the
!                            temperature horizontal basis in km
    real(rp), intent(in) :: orb_inc    ! orbital incline angle in radians
    real(rp), intent(in) :: elev_offset ! elevation offset in radians
    real(rp), intent(in) :: req(:)  ! Earth radius in equivalent circle
!                                  ! at each phi value representation (km)
! Keywords
    real(rp), optional, intent(in) :: h2o_zeta_basis(:) ! h2o zeta basis
    real(rp), optional, intent(in) :: h2o_phi_basis(:) ! h2o phi basis
    real(rp), optional, intent(in) :: h2o_coeffs(:) ! h2o coefficients in long
!                            vector format zetaXphi.
    logical,  optional, intent(in) :: lin_log ! representation type

! outputs

    real(rp), intent(out) :: tan_chi_out(:) ! computed tangent pointing angles
!                                 corrected for refraction in radians
    real(rp), intent(out) :: dx_dh_out(:)   ! d(chi)/dh on the output grid
    real(rp), intent(out) :: dhdz_out(:)    ! dh/dz on the output grid

    real(rp), optional, intent(out) :: dxdt_tan(:,:) ! computed dchi dt.
    real(rp), optional, intent(out) :: d2xdxdt_tan(:,:) ! computed d2chi dxdt.

! the matrix hierarchy is n_out,phi,zeta
! internal variables

    integer(i4) :: n_out, n_t_phi, n_t_zeta, n_h2o_zeta, n_h2o_phi
    integer(i4) :: sv_p, sv_z, ht_i

    real(rp) :: d2hdhdt_tan(size(zetatan), grids_tmp%no_z(1), grids_tmp%no_p(1))
    real(rp) :: dhdt_tan(size(zetatan), grids_tmp%no_z(1), grids_tmp%no_p(1))
    real(rp) :: dhdz_tan(size(zetatan), grids_tmp%no_p(1))
    real(rp), pointer :: eta_p(:,:)
    real(rp) :: eta_t(size(zetatan), grids_tmp%no_p(1))
    real(rp), pointer :: eta_z(:,:)
    real(rp), pointer :: h2o_tan_out(:)
    real(rp) :: height_tan(size(zetatan), grids_tmp%no_p(1))
    real(rp) :: h_tan_out(size(zetatan))
    real(rp) :: n_tan_out(size(zetatan))
    real(rp) :: temp_tan(size(zetatan), grids_tmp%no_p(1))
    real(rp) :: t_tan_out(size(zetatan))

! dimensions we use

    n_out = size(zetatan)
    n_t_phi = grids_tmp%no_p(1)
    n_t_zeta = grids_tmp%no_z(1)

    nullify ( eta_p, eta_z, h2o_tan_out )

    call two_d_hydrostatic ( grids_tmp, ref_zeta, ref_gph, zetatan, orb_inc, &
             & temp_tan, height_tan, dhdz_tan, dhdt_tan, d2hdhdt_tan )

! tangent heights for input pressures along phi

    call get_eta_sparse ( grids_tmp%phi_basis, phitan, eta_t )

    h_tan_out = sum(height_tan * eta_t, dim=2)
    t_tan_out = sum(temp_tan * eta_t, dim=2)
    dhdz_out = sum(dhdz_tan * eta_t, dim=2)

! compute the tangent water vapor profile along input phi_tan

    if ( present(h2o_coeffs) ) then

      n_h2o_zeta = size(h2o_zeta_basis)
      n_h2o_phi  = size(h2o_phi_basis)

      call allocate_test ( eta_p, n_out, n_h2o_phi, 'eta_p', modulename )
      call allocate_test ( eta_z, n_out, n_h2o_zeta, 'eta_z', modulename )
      call allocate_test ( h2o_tan_out, n_out, 'h2o_tan_out', modulename )

      call get_eta_sparse ( h2o_phi_basis, phitan, eta_p )
      call get_eta_sparse ( h2o_zeta_basis, zetatan, eta_z )

! This is actually one less than the true start index

      h2o_tan_out(:) = 0.0_rp
      do sv_p = 1, n_h2o_phi
        do sv_z = 1, n_h2o_zeta
          h2o_tan_out = h2o_tan_out + &
            & h2o_coeffs(sv_z + n_h2o_zeta * (sv_p - 1)) * &
            &   eta_z(:,sv_z) * eta_p(:,sv_p)
        end do
      end do

      if ( lin_log ) h2o_tan_out = exp(h2o_tan_out)

! compute refractive index

      call refractive_index ( 10.0**(-zetatan), t_tan_out, n_tan_out, &
                          & h2o_path=h2o_tan_out )

      call deallocate_test ( h2o_tan_out, 'h2o_tan_out', modulename )
      call deallocate_test ( eta_z, 'eta_z', modulename )
      call deallocate_test ( eta_p, 'eta_p', modulename )

    else

      call refractive_index ( 10.0**(-zetatan), t_tan_out, n_tan_out )

    end if

! compute output angles to interpolate to

    if ( present(dxdt_tan) ) then

      dhdt_tan = dhdt_tan  * spread(eta_t,2,n_t_zeta)
      d2hdhdt_tan = d2hdhdt_tan * spread(eta_t,2,n_t_zeta)
      do ht_i = 1, n_out
        call get_chi_angles ( scgeocalt(ht_i), n_tan_out(ht_i), &
           & h_tan_out(ht_i), phitan(ht_i), req(ht_i), elev_offset, &
           & tan_chi_out(ht_i), dx_dh_out(ht_i), dhdz_out(ht_i), &
           & reshape(dhdt_tan(ht_i,:,:),(/n_t_zeta*n_t_phi/)), &
           & reshape(d2hdhdt_tan(ht_i,:,:),(/n_t_zeta*n_t_phi/)), &
           & dxdt_tan(ht_i,:), d2xdxdt_tan(ht_i,:) )
      end do

     else

       do ht_i = 1, n_out
         call get_chi_angles(scgeocalt(ht_i),n_tan_out(ht_i), &
           & h_tan_out(ht_i),phitan(ht_i),req(ht_i),elev_offset, &
           & tan_chi_out(ht_i),dx_dh_out(ht_i),dhdz_out(ht_i))
       end do

    end if

  end subroutine Get_Chi_Out

end module Get_Chi_Out_m
! $Log$
! Revision 2.10  2002/10/03 22:25:13  vsnyder
! Move USE statements from module scope to procedure scope.  Convert some
! allocated arrays to automatic arrays.  Cosmetic changes.
!
! Revision 2.9  2002/07/11 20:51:49  bill
! made req an mmif quantity
!
! Revision 2.8  2002/07/05 07:52:48  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.7  2002/06/28 11:06:51  zvi
! computes dx_dh & dh_dz on output grid as well
!
! Revision 2.6  2002/06/24 21:11:24  zvi
! Adding Grids_tmp stracture and modifying calling sequences
!
! Revision 1.0  2002/06/24 22:20:45  bill
! Intitial release ...
