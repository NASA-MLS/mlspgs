!
! This is a new module to compute some various spectroscopy path arrays
!
MODULE eval_spect_path_m
!
  use MLSCommon, only: RP, IP
  USE get_eta_matrix_m, ONLY: get_eta_sparse
!
  IMPLICIT NONE

  Private
  Public :: eval_spect_path

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
 "$Id$"
  character (LEN=*), parameter :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
 CONTAINS
!-----------------------------------------------------------------
!
 SUBROUTINE eval_spect_path(z_basis,phi_basis,no_z,no_phi,path_zeta, &
                        &   path_phi,do_calc,eta_zxp)
! Input:
!
  REAL(rp), INTENT(in) :: z_basis(:) ! a vector of zeta basis loaded
!                             sequentially ie [basis_1,basis_2,..., basis_n]
  REAL(rp), INTENT(in) :: phi_basis(:) ! a vector of orbit plane projected
!                                      horizontal bases, entered squentially
!                             ie [basis_1,basis_2,..., basis_n]
  INTEGER(ip), INTENT(in) :: no_z(:) ! number of z_bases by sps number
  INTEGER(ip), INTENT(in) :: no_phi(:) ! number of phi_bases by sps number
  REAL(rp), INTENT(in) :: path_zeta(:) ! zeta values along path
  REAL(rp), INTENT(in) :: path_phi(:) ! phi values along path
!
! Output:
!
  LOGICAL, INTENT(out) :: do_calc(:,:) !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
  REAL(rp), INTENT(out) :: eta_zxp(:,:) ! Eta_z x Eta_phi for each state
!                          vector element. This is the same length as values.
! Internal declaritions
!
  INTEGER(ip) :: n_p, n_z, npxz
  INTEGER(ip) :: i_sps,n_sps,n_path,sv_i,sv_z,sv_p,sv_j
  INTEGER(ip) :: p_inda,z_inda,v_inda,p_indb,z_indb,v_indb

  REAL(rp), ALLOCATABLE :: eta_z(:,:),eta_p(:,:)
  LOGICAL, ALLOCATABLE :: not_zero_z(:,:), not_zero_p(:,:)
!
! Begin executable code:
!
  n_sps = SIZE(no_z)
  n_path = SIZE(path_zeta)
!
  eta_zxp = 0.0
  do_calc = .false.
!
  p_inda = 1
  z_inda = 1
  v_inda = 1
!
  DO i_sps = 1, n_sps
!
    n_z = no_z(i_sps)
    n_p = no_phi(i_sps)
    npxz = n_z * n_p
    if(npxz == 0) CYCLE

    z_indb = z_inda + n_z
    p_indb = p_inda + n_p

    v_indb = v_inda + npxz
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
    CALL get_eta_sparse(z_basis(z_inda:z_indb-1),path_zeta,eta_z,not_zero_z)
    CALL get_eta_sparse(phi_basis(p_inda:p_indb-1),path_phi,eta_p,not_zero_p)
!
    DO sv_i = 0 , npxz - 1
      sv_j = v_inda + sv_i
      sv_z = 1 + MODULO(sv_i,n_z)
      sv_p = 1 + sv_i / n_z
      WHERE(not_zero_z(:,sv_z) .AND. not_zero_p(:,sv_p))
        do_calc(:,sv_j) = .true.
        eta_zxp(:,sv_j) = eta_z(:,sv_z) * eta_p(:,sv_p)
      ENDWHERE
    END DO
!
    z_inda = z_indb
    p_inda = p_indb
    v_inda = v_indb

    DEALLOCATE(not_zero_z)
    DEALLOCATE(not_zero_p)
    DEALLOCATE(eta_z)
    DEALLOCATE(eta_p)
!
  END DO
!
 END SUBROUTINE eval_spect_path
!
END MODULE eval_spect_path_m
!
! $Log$
! Revision 1.1.2.5  2001/09/13 22:51:21  zvi
! Separating allocation stmts
!
! Revision 1.1.2.4  2001/09/12 21:48:49  zvi
! Beautifying ..
!
! Revision 1.1.2.3  2001/09/12 21:47:22  zvi
! Adding CVS stuff
!
