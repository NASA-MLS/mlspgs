! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Comp_Eta_Docalc_No_Frq_m

! This is a new module to compute the SPS path

  implicit NONE
  private
  public :: Comp_Eta_Docalc_No_Frq

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------------

contains

! ---------------------------------------  Comp_Eta_Docalc_No_Frq  -----

  subroutine Comp_Eta_Docalc_No_Frq ( Grids_x, path_zeta, path_phi, &
                                  &   do_calc_zp, eta_zp )

    use MLSCommon, only: RP, IP
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_Sps_Data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates

    real(rp), intent(in) :: path_zeta(:) ! zeta values along path for which
!                           species vmr is needed.
    real(rp), intent(in) :: path_phi(:) ! phi values along path for which
!                           species vmr is needed.
! output:

    logical, intent(out) :: do_calc_zp(:,:) ! Indicates whether there is a
!                           contribution for this state vector element. This
!                           is the same length as values.
    real(rp), intent(out) :: eta_zp(:,:) ! Eta_z * Eta_phi for each state
!                           vector element. This is the same length as values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).

! Internal declarations:

    integer(ip) :: N_p, N_z, Npz
    integer(ip) :: Sps_i, Sv_i, Sv_z, Sv_p, Sv_j
    integer(ip) :: P_inda, Z_inda, V_inda, P_indb, Z_indb, V_indb

    real(rp) :: Eta_p(1:size(path_zeta),1:maxval(Grids_x%no_p))
    real(rp) :: Eta_z(1:size(path_zeta),1:maxval(Grids_x%no_z))
    logical :: Not_zero_p(1:size(path_zeta),1:maxval(Grids_x%no_p))
    logical :: Not_zero_z(1:size(path_zeta),1:maxval(Grids_x%no_z))

! Begin executable code:

    eta_zp = 0.0
    do_calc_zp = .false.

    p_inda = 1
    z_inda = 1
    v_inda = 1

    do sps_i = 1 , size(Grids_x%no_z) ! Number of molecules

      n_z = Grids_x%no_z(sps_i)
      n_p = Grids_x%no_p(sps_i)
      npz = n_z * n_p

      z_indb = z_inda + n_z
      p_indb = p_inda + n_p
      v_indb = v_inda + npz

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)

! Compute etas

      call get_eta_sparse ( Grids_x%zet_basis(z_inda:z_indb-1), path_zeta, &
                       &    eta_z, not_zero_z )
      call get_eta_sparse ( Grids_x%phi_basis(p_inda:p_indb-1), path_phi,  &
                       &    eta_p, not_zero_p )

      do sv_i = 0 , npz - 1
        sv_j = v_inda + sv_i
        sv_z = 1 + modulo(sv_i,n_z)
        sv_p = 1 + sv_i / n_z
        where ( not_zero_z(:,sv_z) .and. not_zero_p(:,sv_p) )
          do_calc_zp(:,sv_j) = .true.
          eta_zp(:,sv_j) = eta_z(:,sv_z) * eta_p(:,sv_p)
        end where
      end do

      z_inda = z_indb
      p_inda = p_indb
      v_inda = v_indb

    end do

  end subroutine Comp_Eta_Docalc_No_Frq

end module Comp_Eta_Docalc_No_Frq_m

! $Log$
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
