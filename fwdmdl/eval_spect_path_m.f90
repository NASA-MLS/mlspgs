! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Eval_Spect_Path_m

  implicit NONE

  private
  public :: Eval_Spect_Path

!---------------------------- RCS Ident Info -------------------------------
  character (LEN=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (LEN=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!-----------------------------------------------------------------

  subroutine Eval_Spect_Path ( Grids_x, lo, sideband,  path_zeta, path_phi,  &
                        &      do_calc, eta_fzp )

! Compute some various spectroscopy path arrays

    use MLSCommon, only: RP, IP, R8
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_Sps_Data_m, only: Grids_T

! Input:
    type (Grids_T), intent(in) :: Grids_x  ! All the needed coordinates
    real(r8), intent(in) :: LO            ! Local oscillator
    integer, intent(in) :: Sideband       ! -1 or 1


    real(rp), intent(in) :: Path_zeta(:) ! zeta values along path
    real(rp), intent(in) :: Path_phi(:)  ! phi values along path

! Output:

    logical, intent(out) :: Do_calc(:,:) !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
    real(rp), intent(out) :: Eta_fzp(:,:) ! Eta_z x Eta_phi x Eta_frq for 
!          each state vector element. This is the same length as values.
! Internal declaritions

    integer(ip) :: N_f, N_p, N_z, Nfzp
    integer(ip) :: Sps_i, Sv_i, Sv_j, Sv_f, Sv_z, Sv_p
    integer(ip) :: P_inda, Z_inda, V_inda, P_indb, Z_indb, V_indb, F_inda, F_indb

    real(rp) :: Frq      ! ** ZEBUG ** this will have to change later 
                         ! for the actuall frequecies array coming in as input

!   real(rp) :: Eta_f(nfrq,maxval(grids_x%no_f))
    real(rp) :: Eta_f(1,maxval(grids_x%no_f)) ! ** ZEBUG, use only one freq.
    real(rp) :: Eta_p(size(path_zeta),maxval(grids_x%no_p))
    real(rp) :: Eta_z(size(path_zeta),maxval(grids_x%no_z))
!   logical :: Not_zero_f(nfrq,maxval(grids_x%no_f))
    logical :: Not_zero_f(1,maxval(grids_x%no_f)) ! ** ZEBUG, use only one freq.
    logical :: Not_zero_p(size(path_zeta),maxval(grids_x%no_p))
    logical :: Not_zero_z(size(path_zeta),maxval(grids_x%no_z))

! Begin executable code:

    frq = 0.0
    eta_fzp = 0.0
    do_calc = .false.

    f_inda = 1
    p_inda = 1
    z_inda = 1
    v_inda = 1

    do sps_i = 1, size(grids_x%no_z)

      n_f = grids_x%no_f(sps_i)
      n_z = grids_x%no_z(sps_i)
      n_p = grids_x%no_p(sps_i)
      nfzp = n_z * n_p * n_f
      if ( nfzp == 0 ) cycle

      f_indb = f_inda + n_f
      z_indb = z_inda + n_z
      p_indb = p_inda + n_p

      v_indb = v_inda + nfzp

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)

! Compute etas

      call get_eta_sparse ( lo+sideband*grids_x%frq_basis(f_inda:f_indb-1), &
                         &  (/frq/), eta_f, not_zero_f )
      call get_eta_sparse ( grids_x%zet_basis(z_inda:z_indb-1),  &
                         &  path_zeta, eta_z, not_zero_z )
      call get_eta_sparse ( grids_x%phi_basis(p_inda:p_indb-1),  &
                         &  path_phi, eta_p, not_zero_p )

      do sv_i = 1, nfzp
        sv_j = v_inda + sv_i - 1
        call brkmod ( sv_i, n_f, n_z, n_p, sv_f, sv_z, sv_p )
        if ( not_zero_f(1,sv_f) ) then
          where ( not_zero_z(:,sv_z) .and. not_zero_p(:,sv_p) )
            do_calc(:,sv_j) = .true.
            eta_fzp(:,sv_j) = eta_f(1,sv_f) * eta_z(:,sv_z) * eta_p(:,sv_p)
          end where
        end if
      end do

      f_inda = f_indb
      z_inda = z_indb
      p_inda = p_indb
      v_inda = v_indb

    end do

  contains

    subroutine BrkMod ( m, nf, nz, np, jf, jz, jp )

     Integer, intent(in) :: m, nf, nz, np
     Integer, intent(out) :: jf, jz, jp

     Integer :: k

      jf = 1 + MOD(m-1,nf)
      k = (m-jf+nf-1)/nf
      jz = 1 + MOD(k,nz)
      k = (k-jz+nz-1)/nz
      jp = 1 + MOD(k,np)

    end subroutine BrkMod

  end subroutine Eval_Spect_Path

!---------------------------------------------------------------------

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Eval_Spect_Path_m
!
! $Log$
! Revision 2.5  2002/09/07 02:17:52  vsnyder
! Move USEs from module scope to procedure scope.
! Convert some allocatable arrays to automatics.  Cosmetic changes.
!
! Revision 2.4  2002/08/22 23:13:33  livesey
! New frq_Basis based on intermediate frequency
!
! Revision 2.3  2002/07/08 17:45:39  zvi
! cleaner code to find indecies
!
! Revision 2.2  2002/01/27 08:37:48  zvi
! Adding Users selected coefficients for derivatives
!
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
