!
! This is a new module to compute the sps path
!
MODULE comp_eta_docalc_no_frq_m
!
  use MLSCommon, only: RP, IP
  USE get_eta_matrix_m, ONLY: get_eta_sparse
  use Load_sps_data_m, only: Grids_T
!
  IMPLICIT NONE

  Private
  Public :: comp_eta_docalc_no_frq

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
 "$Id$"
  character (LEN=*), parameter :: ModuleName= &
 "$RCSfile$"
!-----------------------------------------------------------------
 CONTAINS
!-----------------------------------------------------------------
!
  SUBROUTINE comp_eta_docalc_no_frq(Grids_x,path_zeta,path_phi, &
                                &   do_calc_zp,eta_zp)
! Input:
!
  type (Grids_T), INTENT(in) :: Grids_x  ! All the needed coordinates

  REAL(rp), INTENT(in) :: path_zeta(:) ! zeta values along path for which
!                         species vmr is needed.
  REAL(rp), INTENT(in) :: path_phi(:) ! phi values along path for which
!                         species vmr is needed.
! Output:
!
  LOGICAL, INTENT(out) :: do_calc_zp(:,:) !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
  REAL(rp), INTENT(out) :: eta_zp(:,:) ! Eta_z x Eta_phi for each state
!                          vector element. This is the same length as values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
!
! Internal declaritions
!
  INTEGER(ip) :: n_p, n_z, npz
  INTEGER(ip) :: i_sps,n_sps,n_path,sv_i,sv_z,sv_p,sv_j
  INTEGER(ip) :: p_inda,z_inda,v_inda,p_indb,z_indb,v_indb

  REAL(rp), ALLOCATABLE :: eta_z(:,:),eta_p(:,:)
  LOGICAL, ALLOCATABLE :: not_zero_z(:,:), not_zero_p(:,:)
!
! Begin executable code:
!
  n_sps = SIZE(Grids_x%no_z)
  n_path = SIZE(path_zeta)
!
  eta_zp = 0.0
  do_calc_zp = .false.
!
  p_inda = 1
  z_inda = 1
  v_inda = 1
!
  DO i_sps = 1 , n_sps
!
    n_z = Grids_x%no_z(i_sps)
    n_p = Grids_x%no_p(i_sps)
    npz = n_z * n_p

    z_indb = z_inda + n_z
    p_indb = p_inda + n_p
    v_indb = v_inda + npz
!
! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)
!
    ALLOCATE(eta_p(1:n_path,1:n_p))
    ALLOCATE(eta_z(1:n_path,1:n_z))
    ALLOCATE(not_zero_p(1:n_path,1:n_p))
    ALLOCATE(not_zero_z(1:n_path,1:n_z))
!
! Compute etas
!
    CALL get_eta_sparse(Grids_x%zet_basis(z_inda:z_indb-1),path_zeta, &
                     &  eta_z,not_zero_z)
    CALL get_eta_sparse(Grids_x%phi_basis(p_inda:p_indb-1),path_phi,  &
                     &  eta_p,not_zero_p)
!
    DO sv_i = 0 , npz - 1
      sv_j = v_inda + sv_i
      sv_z = 1 + MODULO(sv_i,n_z)
      sv_p = 1 + sv_i / n_z
      WHERE(not_zero_z(:,sv_z) .AND. not_zero_p(:,sv_p))
        do_calc_zp(:,sv_j) = .true.
        eta_zp(:,sv_j) = eta_z(:,sv_z) * eta_p(:,sv_p)
      ENDWHERE
    END DO

    z_inda = z_indb
    p_inda = p_indb
    v_inda = v_indb

    DEALLOCATE(eta_z)
    DEALLOCATE(eta_p)
    DEALLOCATE(not_zero_z)
    DEALLOCATE(not_zero_p)
!
  END DO
!
  END SUBROUTINE comp_eta_docalc_no_frq
!
END MODULE comp_eta_docalc_no_frq_m
!
! $Log$
! Revision 2.1  2001/11/02 10:47:37  zvi
! Implementing frequecy grid
!
! ! Revision 1.0 2001/10/30 14:00:00  zvi
