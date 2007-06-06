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
                                  &   eta_zp, do_calc_zp, sps, tan_pt )

    use MLSCommon, only: RP, IP
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
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
    integer(ip) :: Sps_1, Sps_n, Sps_i, Sv_z, Sv_p
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

    sps_1 = 1
    sps_n = ubound(Grids_x%l_z,1)
    if ( present(sps) ) then
      sps_1 = sps
      sps_n = sps
    end if

    do sps_i = sps_1, sps_n ! Number of molecules

      p_indb = Grids_x%l_p(sps_i)
      z_indb = Grids_x%l_z(sps_i)

      n_z = z_indb - z_inda
      n_p = p_indb - p_inda

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease let's do the slow and easy (and certainly more reliable)

! Compute etas

      if ( present(do_calc_zp) ) then
        if ( present(tan_pt) ) then
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), &
                           &    path_zeta(tan_pt:1:-1), eta_z(tan_pt:1:-1,:), &
                           &    not_zero_z(tan_pt:1:-1,:), sorted=.true. )
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), &
                           &    path_zeta(tan_pt+1:), eta_z(tan_pt+1:,:), &
                           &    not_zero_z(tan_pt+1:,:), sorted=.true. )
        else
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
                           &    eta_z, not_zero_z )
        end if
        call get_eta_sparse ( Grids_x%phi_basis(p_inda+1:p_indb), path_phi,  &
                         &    eta_p, not_zero_p, sorted=.true. )
      else
        if ( present(tan_pt) ) then
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), &
                           &    path_zeta(tan_pt:1:-1), eta_z(tan_pt:1:-1,:), &
                           &    sorted=.true. )
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), &
                           &    path_zeta(tan_pt+1:), eta_z(tan_pt+1:,:), &
                           &    sorted=.true. )
        else
          call get_eta_sparse ( Grids_x%zet_basis(z_inda+1:z_indb), path_zeta, &
                           &    eta_z )
        end if
        call get_eta_sparse ( Grids_x%phi_basis(p_inda+1:p_indb), path_phi,  &
                         &    eta_p, sorted=.true. )
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
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Comp_Eta_Docalc_No_Frq_m

! $Log$
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
