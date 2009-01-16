! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Eta_Docalc_No_Frq_m

! This is a new module to compute the SPS path

  implicit NONE
  private
  public :: Comp_Eta_Docalc_No_Frq

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ---------------------------------------  Comp_Eta_Docalc_No_Frq  -----

  subroutine Comp_Eta_Docalc_No_Frq ( Grids_x, path_zeta, path_phi, &
                                  &   eta_zp, do_calc_zp, sps, tan_pt, &
                                  &   nz_zp, nnz_zp )

    use MLSCommon, only: RP, IP
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse, Multiply_Eta_Column_Sparse
    use Load_Sps_Data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates

    real(rp), intent(in) :: path_zeta(:) ! zeta values along path for which
!                           species vmr is needed.
    real(rp), intent(in) :: path_phi(:) ! phi values along path for which
!                           species vmr is needed.
    integer, intent(in), optional :: SPS ! Only do this species if present
    integer, intent(in), optional :: Tan_Pt ! Tangent point; path_zeta is sorted
!                           before and after tangent point

! Inout:

    logical, intent(inout), optional :: do_calc_zp(:,:) ! Indicates whether there
!                           is a contribution for this state vector element.
!                           This is the same length as values.
    integer, intent(inout), optional :: NZ_ZP(:,:) ! Nonzeros in eta_zp
    integer, intent(inout), optional :: NNZ_ZP(:)  ! Number of nonzeros in eta_zp

! Output:

    real(rp), intent(out) :: eta_zp(:,:) ! Eta_z * Eta_phi for each state
!                           vector element. This is the same length as values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).

! Internal declarations:

    integer(ip) :: N_p, N_z, N_v
    integer(ip) :: Sps_1, Sps_n, Sps_i
    integer(ip) :: My_Tan, P_inda, V_Inda, V_Indb, Z_inda, P_indb, Z_indb

    real(rp) :: Eta_p(1:size(path_phi), &  ! size(path_phi) == size(path_zeta)
      & maxval(Grids_x%l_p(1:ubound(Grids_x%l_p,1))-Grids_x%l_p(0:ubound(Grids_x%l_p,1)-1)))
    real(rp) :: Eta_z(1:size(path_zeta), & ! == size(path_zeta)
      & maxval(Grids_x%l_z(1:ubound(Grids_x%l_z,1))-Grids_x%l_z(0:ubound(Grids_x%l_z,1)-1)))
    integer :: NZ_P(1:size(Eta_p,1),1:size(Eta_p,2))
    integer :: NZ_Z(1:size(Eta_z,1),1:size(Eta_z,2))
    integer :: NNZ_P(1:size(Eta_p,2))
    integer :: NNZ_Z(1:size(Eta_z,2))
    integer :: MY_NZ_ZP(size(eta_zp,1),size(eta_p,2)*size(eta_z,2))
    integer :: MY_NNZ_ZP(size(eta_p,2)*size(eta_z,2))

! Begin executable code:

    if ( .not. present(nz_zp) ) then
      eta_zp = 0.0
      if ( present(do_calc_zp) ) do_calc_zp = .false.
      my_nnz_zp = 0
    end if

    nnz_z = 0
    nnz_p = 0

    my_tan = size(path_zeta) / 2
    if ( present(tan_pt) ) my_tan = tan_pt

    sps_1 = 1
    sps_n = ubound(Grids_x%l_z,1)
    if ( present(sps) ) then
      sps_1 = sps
      sps_n = sps
    end if

    p_indb = Grids_x%l_p(sps_1-1)
    z_indb = Grids_x%l_z(sps_1-1)
    v_indb = 0

    do sps_i = sps_1, sps_n

      p_inda = p_indb
      z_inda = z_indb
      v_inda = v_indb

      p_indb = Grids_x%l_p(sps_i)
      z_indb = Grids_x%l_z(sps_i)

      n_z = z_indb - z_inda
      n_p = p_indb - p_inda
      n_v = n_z * n_p

      v_indb = v_inda + n_v

      call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
      &    eta_z, my_tan, 1, nz_z, nnz_z, .false. )
      call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
      &    eta_z, my_tan+1, size(path_zeta), nz_z, nnz_z, .true. )
      call get_eta_sparse ( Grids_x%phi_basis(p_inda+1:p_indb), path_phi,  &
      &    eta_p, 1, size(path_phi), nz_p, nnz_p, .false. )
      if ( present(nz_zp) ) then
        if ( present(do_calc_zp) ) then
          call multiply_eta_column_sparse ( &
            & eta_z(:,:n_z), nz_z(:,:n_z), nnz_z(:n_z), &
            & eta_p(:,:n_p), nz_p(:,:n_p), nnz_p(:n_p), &
            & eta_zp(:,v_inda+1:v_indb), nz_zp(:,v_inda+1:v_indb), &
            & nnz_zp(v_inda+1:v_indb), do_calc_zp(:,v_inda+1:v_indb) )
        else
          call multiply_eta_column_sparse ( &
            & eta_z(:,:n_z), nz_z(:,:n_z), nnz_z(:n_z), &
            & eta_p(:,:n_p), nz_p(:,:n_p), nnz_p(:n_p), &
            & eta_zp(:,v_inda+1:v_indb), nz_zp(:,v_inda+1:v_indb), &
            & nnz_zp(v_inda+1:v_indb) )
        end if
      else
        if ( present(do_calc_zp) ) then
          call multiply_eta_column_sparse ( &
            & eta_z(:,:n_z), nz_z(:,:n_z), nnz_z(:n_z), &
            & eta_p(:,:n_p), nz_p(:,:n_p), nnz_p(:n_p), &
            & eta_zp(:,v_inda+1:v_indb), my_nz_zp(:,:n_v), my_nnz_zp(:n_v), &
            & do_calc_zp(:,v_inda+1:v_indb) )
        else
          call multiply_eta_column_sparse ( &
            & eta_z(:,:n_z), nz_z(:,:n_z), nnz_z(:n_z), &
            & eta_p(:,:n_p), nz_p(:,:n_p), nnz_p(:n_p), &
            & eta_zp(:,v_inda+1:v_indb), my_nz_zp(:,:n_v), my_nnz_zp(:n_v) )
        end if
      end if

    end do

  end subroutine Comp_Eta_Docalc_No_Frq

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Comp_Eta_Docalc_No_Frq_m

! $Log$
! Revision 2.12  2009/01/16 23:47:47  vsnyder
! Cannonball polishing
!
! Revision 2.11  2007/07/25 20:21:10  vsnyder
! Delete declarations for unused variables
!
! Revision 2.10  2007/06/26 00:37:01  vsnyder
! Use column-sparse eta
!
! Revision 2.9  2007/06/06 01:12:17  vsnyder
! Add tangent point optional argument
!
! Revision 2.8  2005/12/22 20:48:55  vsnyder
! Add a 'this-species-only' optional argument
!
! Revision 2.7  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.6  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.5.2.3  2003/03/22 02:38:01  vsnyder
! Make do_calc_zp optional, don't compute stuff not needed if it's not present
!
! Revision 2.5.2.2  2003/03/22 02:31:20  vsnyder
! Remove a WHERE that didn't save anything, cosmetic changes
!
! Revision 2.5.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.5  2002/10/08 17:08:01  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/09/06 18:18:03  vsnyder
! Cosmetic changes.  Move USEs from module scope to procedure scope.
! Convert some arrays from pointers to automatics.
!
! Revision 2.3  2002/06/04 10:28:01  zvi
! rename n_sps to: no_mol, more correctly
!
! Revision 2.2  2002/01/27 08:37:46  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.1  2001/11/02 10:47:37  zvi
! Implementing frequecy grid
!
! Revision 1.0 2001/10/30 14:00:00  zvi
