!
! This is a new module to compute the sps path
!
MODULE comp_sps_path_frq_m
!
  use MLSCommon, only: RP, IP
  USE get_eta_matrix_m, ONLY: get_eta_sparse
  use Load_sps_data_m, only: Grids_T
!
  IMPLICIT NONE

  Private
  Public :: comp_sps_path_frq

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=256) :: Id = &
 "$Id$"
  character (LEN=*), parameter :: ModuleName= &
 "$RCSfile$"
!-----------------------------------------------------------------
 CONTAINS
!-----------------------------------------------------------------
!
  SUBROUTINE comp_sps_path_frq(Grids_x,Frq,values,eta_zp,do_calc_zp, &
              &   skip_eta_frq,lin_log,sps_values,do_calc_fzp,eta_fzp)
!
! Input:
!
  type (Grids_T), INTENT(in) :: Grids_x  ! All the needed coordinates

  REAL(rp), INTENT(in) :: Frq  ! Frequency at which to compute the values
  REAL(rp), INTENT(in) :: values(:) ! A vector of coefficient break-point
!                         values entered sequentially according to the
!                         heirarch, z_basis, phi_basis, sps_number
  REAL(rp), INTENT(in) :: eta_zp(:,:) ! Eta_z x Eta_phi for each state
!                         vector element. This is the same length as values.
  LOGICAL, INTENT(in) :: do_calc_zp(:,:) !logical indicating whether there
!                        is a contribution for this state vector element
  LOGICAL, INTENT(in) :: skip_eta_frq(:) ! Flags for specie with NO Freq. dim.
  LOGICAL, INTENT(in) :: lin_log(:) ! species representation basis type
!
! Output:
!
  REAL(rp), INTENT(out) :: sps_values(:,:) ! vmr values along the path
!                          by species number
  REAL(rp), INTENT(out) :: eta_fzp(:,:) ! Eta_z x Eta_phi x Eta_f for each
!                          state vector element. This is the same length as
!                          values.
  LOGICAL, INTENT(out) :: do_calc_fzp(:,:) !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_values = values.
!
! Internal declaritions
!
  INTEGER(ip) :: n_zp, n_f, nfzp,iw
  INTEGER(ip) :: sps_i,n_sps,sv_i,sv_zp,sv_f,sv_j
  INTEGER(ip) :: v_inda,v_indb,f_inda,f_indb,w_inda,w_indb

  REAL(rp), ALLOCATABLE :: eta_f(:,:)
  LOGICAL, ALLOCATABLE :: not_zero_f(:,:)
!
! Begin executable code:
!
  n_sps = SIZE(Grids_x%no_z)
!
  IF(Frq < 1.0) THEN
    eta_fzp = 0.0_rp
    sps_values = 0.0_rp
    do_calc_fzp = .false.
  ELSE
    DO sps_i = 1, n_sps
      IF(.NOT. skip_eta_frq(sps_i)) sps_values(:,sps_i) = 0.0_rp
    END DO
  ENDIF

!
  f_inda = 1
  v_inda = 1
  w_inda = 1
!
  DO sps_i = 1, n_sps
!
    n_f = Grids_x%no_f(sps_i)
    n_zp = Grids_x%no_z(sps_i) * Grids_x%no_p(sps_i)
    nfzp = n_f * n_zp

    f_indb = f_inda + n_f
    v_indb = v_inda + nfzp
    w_indb = w_inda + n_zp

    IF((Frq > 1.0) .AND. skip_eta_frq(sps_i)) THEN
      f_inda = f_indb
      v_inda = v_indb
      w_inda = w_indb
      CYCLE
    ENDIF
!
! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)
!
    ALLOCATE(eta_f(1:1,1:n_f))
    ALLOCATE(not_zero_f(1:1,1:n_f))
!
! Compute eta:
!
    CALL get_eta_sparse(Grids_x%frq_basis(f_inda:f_indb-1),(/Frq/), &
                    &   eta_f,not_zero_f)
!
    IF(lin_log(sps_i)) THEN

      DO sv_i = 0 , nfzp - 1
        sv_j = v_inda + sv_i
        sv_f = 1 + MODULO(sv_i,n_f)
        sv_zp = w_inda + sv_i / n_f
        IF(not_zero_f(1,sv_f)) THEN
          WHERE(do_calc_zp(:,sv_zp))
            do_calc_fzp(:,sv_j) = .true.
            eta_fzp(:,sv_j) = eta_f(1,sv_f) * eta_zp(:,sv_zp)
            sps_values(:,sps_i) = sps_values(:,sps_i) +  &
                               &  EXP(values(sv_j) * eta_fzp(:,sv_j))
          ENDWHERE
        ENDIF
      END DO

    ELSE

      DO sv_i = 0 , nfzp - 1
        sv_j = v_inda + sv_i
        sv_f = 1 + MODULO(sv_i,n_f)
        sv_zp = w_inda + sv_i / n_f
        IF(not_zero_f(1,sv_f)) THEN
          WHERE(do_calc_zp(:,sv_zp))
            do_calc_fzp(:,sv_j) = .true.
            eta_fzp(:,sv_j) = eta_f(1,sv_f) * eta_zp(:,sv_zp)
            sps_values(:,sps_i) = sps_values(:,sps_i) +  &
                               &  values(sv_j) * eta_fzp(:,sv_j)
          ENDWHERE
        ENDIF
      END DO

    ENDIF
!
    f_inda = f_indb
    v_inda = v_indb
    w_inda = w_indb

    DEALLOCATE(eta_f)
    DEALLOCATE(not_zero_f)
!
  END DO
!
  END SUBROUTINE comp_sps_path_frq
!
END MODULE comp_sps_path_frq_m
!
! $Log$
! Revision 2.3  2001/11/10 00:45:08  zvi
! Fixing a bug..
!
! Revision 2.2  2001/11/07 09:59:12  zvi
! More effective code for sps_path calculations
!
! Revision 1.0  2001/10/30 14:00:00 zvi Exp $"

