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
    real(r8), intent(in) :: LO             ! Local oscillator
    integer, intent(in) :: Sideband        ! -1 or 1


    real(rp), intent(in) :: Path_zeta(:)   ! zeta values along path
    real(rp), intent(in) :: Path_phi(:)    ! phi values along path

! Output:

    logical, intent(out) :: Do_calc(:,:)   !logical indicating whether there
!                         is a contribution for this state vector element
!                         This is the same length as values.
    real(rp), intent(out) :: Eta_fzp(:,:)  ! Eta_z x Eta_phi x Eta_frq for 
!          each state vector element. This is the same length as values.

! Internal declarations

    integer(ip) :: Sps_i, Sv_f, Sv_p, Sv_z, Sv_v
    integer(ip) :: P_inda, Z_inda, V_inda, P_indb, Z_indb, V_indb, F_inda, F_indb

    real(rp) :: Frq      ! ** ZEBUG ** this will have to change later 
                         ! for the actuall frequecies array coming in as input

!   Actually, all we need is the largest of any grids_x%l_X(i)-grids_x%l_X(i-1)
!   real(rp) :: Eta_f(nfrq,grids_x%l_f(ubound(grids_x%l_f,1)))
    real(rp) :: Eta_f(1,grids_x%l_f(ubound(grids_x%l_f,1))) ! ** ZEBUG, use only one freq.
    real(rp) :: Eta_p(size(path_zeta),grids_x%l_p(ubound(grids_x%l_p,1)))
    real(rp) :: Eta_z(size(path_zeta),grids_x%l_z(ubound(grids_x%l_z,1)))
!   logical :: Not_zero_f(nfrq,grids_x%l_f(ubound(grids_x%l_f,1)))
    logical :: Not_zero_f(1,grids_x%l_f(ubound(grids_x%l_f,1))) ! ** ZEBUG, use only one freq.
    logical :: Not_zero_p(size(path_zeta),grids_x%l_p(ubound(grids_x%l_p,1)))
    logical :: Not_zero_z(size(path_zeta),grids_x%l_z(ubound(grids_x%l_z,1)))

! Begin executable code:

    frq = 0.0

    f_inda = 0
    p_inda = 0
    z_inda = 0
    v_inda = 0

    do sps_i = 1, ubound(grids_x%l_z,1)

      v_indb = grids_x%l_v(sps_i)
      if ( v_indb == v_inda ) cycle

      f_indb = grids_x%l_f(sps_i)
      p_indb = grids_x%l_p(sps_i)
      z_indb = grids_x%l_z(sps_i)

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease let's do the slow and easy (and certainly more reliable)

! Compute etas

      call get_eta_sparse ( lo+sideband*grids_x%frq_basis(f_inda+1:f_indb), &
                         &  (/frq/), eta_f, not_zero_f )
      call get_eta_sparse ( grids_x%zet_basis(z_inda+1:z_indb),  &
                         &  path_zeta, eta_z, not_zero_z )
      call get_eta_sparse ( grids_x%phi_basis(p_inda+1:p_indb),  &
                         &  path_phi, eta_p, not_zero_p )

      sv_v = v_inda
      do sv_p = p_inda+1, p_indb
        do sv_z = z_inda+1, z_indb
          do sv_f = f_inda+1, f_indb
            sv_v = sv_v + 1
            if ( not_zero_f(1,sv_f) ) then
              do_calc(:,sv_v) = not_zero_z(:,sv_z) .and. not_zero_p(:,sv_p)
              where ( do_calc(:,sv_v) )
                eta_fzp(:,sv_v) = eta_f(1,sv_f) * eta_z(:,sv_z) * eta_p(:,sv_p)
              elsewhere
                eta_fzp(:,sv_v) = 0.0
              end where
            end if
          end do
        end do
      end do

      f_inda = f_indb
      z_inda = z_indb
      p_inda = p_indb
      v_inda = v_indb

    end do

  end subroutine Eval_Spect_Path

!---------------------------------------------------------------------

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Eval_Spect_Path_m
!
! $Log$
! Revision 2.7  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.6.2.1  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.6  2002/10/08 17:08:03  pwagner
! Added idents to survive zealous Lahey optimizer
!
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
