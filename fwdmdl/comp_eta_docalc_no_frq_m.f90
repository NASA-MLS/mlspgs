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
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ---------------------------------------  Comp_Eta_Docalc_No_Frq  -----

  subroutine Comp_Eta_Docalc_No_Frq ( Grids_x, path_zeta, path_phi, &
                                  &   eta_zp, do_calc_zp )

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

    real(rp), intent(out) :: eta_zp(:,:) ! Eta_z * Eta_phi for each state
!                           vector element. This is the same length as values.
    logical, intent(out), optional :: do_calc_zp(:,:) ! Indicates whether there
!                           is a contribution for this state vector element.
!                           This is the same length as values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).

! Internal declarations:

    integer(ip) :: N_p, N_z
    integer(ip) :: Sps_i, Sv_z, Sv_p
    integer(ip) :: P_inda, Z_inda, V_inda, P_indb, Z_indb

    real(rp) :: Eta_p(1:size(path_zeta),1:Grids_x%l_p(ubound(Grids_x%l_p,1)))
    real(rp) :: Eta_z(1:size(path_zeta),1:Grids_x%l_z(ubound(Grids_x%l_z,1)))
    logical :: Not_zero_p(1:size(path_zeta),1:Grids_x%l_p(ubound(Grids_x%l_p,1)))
    logical :: Not_zero_z(1:size(path_zeta),1:Grids_x%l_z(ubound(Grids_x%l_z,1)))

! Begin executable code:

    eta_zp = 0.0

    p_inda = 0
    v_inda = 0
    z_inda = 0

    do sps_i = 1 , ubound(Grids_x%l_z,1) ! Number of molecules

      p_indb = Grids_x%l_p(sps_i)
      z_indb = Grids_x%l_z(sps_i)

      n_z = z_indb - z_inda
      n_p = p_indb - p_inda

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease let's do the slow and easy (and certainly more reliable)

! Compute etas

      if ( present(do_calc_zp) ) then
        call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
                         &    eta_z, not_zero_z )
        call get_eta_sparse ( Grids_x%phi_basis(p_inda+1:p_indb), path_phi,  &
                         &    eta_p, not_zero_p )
      else
        call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
                         &    eta_z )
        call get_eta_sparse ( Grids_x%phi_basis(p_inda+1:p_indb), path_phi,  &
                         &    eta_p )
      end if

      do sv_p = 0, n_p - 1
        do sv_z = 0, n_z - 1
          v_inda = v_inda + 1
          if ( present(do_calc_zp) ) &
          & do_calc_zp(:,v_inda) = not_zero_z(:,sv_z+1) .and. not_zero_p(:,sv_p+1)
          ! removed "where ( do_calc_zp(:,v_inda) )" because a multiply is cheaper
          eta_zp(:,v_inda) = eta_z(:,sv_z+1) * eta_p(:,sv_p+1)
        end do
      end do

      z_inda = z_indb
      p_inda = p_indb

    end do

  end subroutine Comp_Eta_Docalc_No_Frq

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Comp_Eta_Docalc_No_Frq_m

! $Log$
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
