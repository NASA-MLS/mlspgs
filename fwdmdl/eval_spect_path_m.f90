!
! This is a new module to compute some various spectroscopy path arrays
!
MODULE eval_spect_path_m
!
  use MLSCommon, only: RP, IP
  USE get_eta_matrix_m, ONLY: get_eta_sparse
  use Load_sps_data_m, only: Grids_T
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
 SUBROUTINE eval_spect_path(Grids_x,path_zeta,path_phi, &
                        &   do_calc,eta_fzp)
! Input:
  type (Grids_T), INTENT(in) :: Grids_x  ! All the needed coordinates
!
  REAL(rp), INTENT(in) :: path_zeta(:) ! zeta values along path
  REAL(rp), INTENT(in) :: path_phi(:)  ! phi values along path
!
! Output:
!
  LOGICAL, INTENT(out) :: do_calc(:,:) !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
  REAL(rp), INTENT(out) :: eta_fzp(:,:) ! Eta_z x Eta_phi x Eta_frq for 
!          each state vector element. This is the same length as values.
! Internal declaritions
!
  INTEGER(ip) :: n_f, n_p, n_z, nfzp
  INTEGER(ip) :: sps_i,n_sps,n_path,sv_i,sv_j,sv_f,sv_z,sv_p
  INTEGER(ip) :: p_inda,z_inda,v_inda,p_indb,z_indb,v_indb,f_inda,f_indb

  REAL(rp) :: Frq      ! ** ZEBUG ** this will have to change later 
                       ! for the actuall frequecies array coming in as input

  REAL(rp), ALLOCATABLE :: eta_z(:,:),eta_p(:,:),eta_f(:,:)
  LOGICAL, ALLOCATABLE :: not_zero_z(:,:),not_zero_p(:,:),not_zero_f(:,:)
!
! Begin executable code:
!
  n_sps = SIZE(Grids_x%no_z)
  n_path = SIZE(path_zeta)
!
  Frq = 0.0
  eta_fzp = 0.0
  do_calc = .false.
!
  f_inda = 1
  p_inda = 1
  z_inda = 1
  v_inda = 1
!
  DO sps_i = 1, n_sps
!
    n_f = Grids_x%no_f(sps_i)
    n_z = Grids_x%no_z(sps_i)
    n_p = Grids_x%no_p(sps_i)
    nfzp = n_z * n_p * n_f
    if(nfzp == 0) CYCLE

    f_indb = f_inda + n_f
    z_indb = z_inda + n_z
    p_indb = p_inda + n_p

    v_indb = v_inda + nfzp
!
! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)
!
    ALLOCATE(eta_z(1:n_path,1:n_z))
    ALLOCATE(eta_p(1:n_path,1:n_p))

    ALLOCATE(not_zero_z(1:n_path,1:n_z))
    ALLOCATE(not_zero_p(1:n_path,1:n_p))

!   ALLOCATE(eta_f(1:nfrq,1:n_f))
    ALLOCATE(eta_f(1:1,1:n_f))         ! ** ZEBUG, use only one freq.

!   ALLOCATE(not_zero_f(1:nfrq,1:n_f))
    ALLOCATE(not_zero_f(1:1,1:n_f))    ! ** ZEBUG, use only one freq.
!
! Compute etas
!
    CALL get_eta_sparse(Grids_x%frq_basis(f_inda:f_indb-1),(/Frq/), &
                     &  eta_f,not_zero_f)
    CALL get_eta_sparse(Grids_x%zet_basis(z_inda:z_indb-1),path_zeta, &
                     &  eta_z,not_zero_z)
    CALL get_eta_sparse(Grids_x%phi_basis(p_inda:p_indb-1),path_phi, &
                     &  eta_p,not_zero_p)
!
    DO sv_i = 1, nfzp
      sv_j = v_inda + sv_i - 1
      Call BrkMod(sv_i,n_f,n_z,n_p,sv_f,sv_z,sv_p)
      IF(not_zero_f(1,sv_f)) THEN
        WHERE(not_zero_z(:,sv_z) .AND. not_zero_p(:,sv_p))
          do_calc(:,sv_j) = .true.
          eta_fzp(:,sv_j) = eta_f(1,sv_f) * eta_z(:,sv_z) * eta_p(:,sv_p)
        ENDWHERE
      ENDIF
    END DO
!
    f_inda = f_indb
    z_inda = z_indb
    p_inda = p_indb
    v_inda = v_indb

    DEALLOCATE(not_zero_f)
    DEALLOCATE(not_zero_z)
    DEALLOCATE(not_zero_p)

    DEALLOCATE(eta_f)
    DEALLOCATE(eta_z)
    DEALLOCATE(eta_p)
!
  END DO
!
 END SUBROUTINE eval_spect_path

!---------------------------------------------------------------------

 SUBROUTINE BrkMod(m,nf,nz,np,jf,jz,jp)

 Integer, intent(in) :: m, nf, nz, np
 Integer, intent(out) :: jf, jz, jp

 Integer :: sv_p,sv_z,sv_f,k
 
 jf = 0
 jz = 0
 jp = 0

 k = 0
 do sv_p = 1, np
   do sv_z = 1, nz
     do sv_f = 1, nf
       k = k + 1
       if(k == m) then
         jf = sv_f
         jz = sv_z
         jp = sv_p
         Return
       endif
     end do
   end do
 end do

 jf = nf
 jz = nz
 jp = np

 END SUBROUTINE BrkMod
!---------------------------------------------------------------------

!
END MODULE eval_spect_path_m
!
! $Log$
! Revision 2.1  2001/11/02 10:48:39  zvi
! Implementing frequecy grid
!
! Revision 2.0  2001/09/17 20:26:26  livesey
! New forward model
!
! Revision 1.1.2.5  2001/09/13 22:51:21  zvi
! Separating allocation stmts
!
! Revision 1.1.2.4  2001/09/12 21:48:49  zvi
! Beautifying ..
!
! Revision 1.1.2.3  2001/09/12 21:47:22  zvi
! Adding CVS stuff
!
