! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Comp_Sps_Path_Frq_m

  implicit NONE

  private
  public :: Comp_Sps_Path_Frq, Comp_Sps_Path

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!-----------------------------------------------------------------
 contains
!-----------------------------------------------------------------

! --------------------------------------------  Comp_Sps_Path_Frq  -----
  subroutine Comp_Sps_Path_Frq ( Grids_x, lo, sideband, Frq, eta_zp, &
    & do_calc_zp, sps_path, do_calc_fzp, eta_fzp, DoPFA )

! Compute the SPS path

    use MLSCommon, only: RP, IP, R8
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
    real(r8), intent(in) :: LO             ! Local oscillator freq
    integer, intent(in) :: Sideband        ! -1 or 1

    real(r8), intent(in) :: Frq  ! Frequency at which to compute the values
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z x Eta_phi)
!                         First dimension is same as sps_values.
    logical, intent(in) :: Do_Calc_Zp(:,:) ! logical indicating whether there
!                         is a contribution for this state vector element

! Output:

    real(rp), intent(inout) :: Sps_Path(:,:) ! Path X Species.  vmr values
!                         along the path by species number
    logical, intent(inout) :: Do_Calc_Fzp(:,:) ! indicates whether there
!                         is a contribution for this state vector element.
!                         Same shape as Eta_Fzp.
    real(rp), intent(inout) :: Eta_Fzp(:,:)  ! Path X (Eta_f x Eta_z x Eta_phi)
!                         First dimension is same as sps_values.

! Optional:

    logical, intent(in), optional :: DoPFA ! Doing PFA, default false.

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    logical :: MyPFA
    integer(ip) :: n_f
    integer(ip) :: sps_i, no_mol, sv_zp, sv_f
    integer(ip) :: v_inda, f_inda, f_indb, w_inda, w_indb

    real(rp) :: eta_f(1:1,1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))
    logical :: not_zero_f(1:1,1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))

! Begin executable code:

    myPFA = .false.
    if ( present(doPFA) ) myPFA = doPFA

    no_mol = grids_x%lastNonPFA

    if ( frq <= 1.0_r8 ) then
      do_calc_fzp = .FALSE.
      eta_f = 1.0_rp
      eta_fzp = 0.0_rp
      not_zero_f = .TRUE.
      sps_path = 0.0_rp
    else
      do sps_i = 1, no_mol
        if ( grids_x%l_f(sps_i) - grids_x%l_f(sps_i-1) > 1 ) sps_path(:,sps_i) = 0.0_rp
      end do
    end if

    f_inda = 0
    w_inda = 0

    do sps_i = 1, no_mol

      f_indb = grids_x%l_f(sps_i)
      n_f = f_indb - f_inda

      w_indb = w_inda + (Grids_x%l_z(sps_i) - Grids_x%l_z(sps_i-1)) * &
                        (Grids_x%l_p(sps_i) - Grids_x%l_p(sps_i-1))

      ! n_f == 1 means qty%template%frequencyCoordinate == l_none, i.e.,
      ! sps_path is not frequency dependent.
      if ( (frq <= 1.0 .or. n_f /= 1) .and. &
        &  ( (sps_i > grids_x%lastNonPFA) .eqv. myPFA) ) then

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)

! Compute eta_f:
        if ( frq > 1.0_rp ) then
          if ( sideband == -1 ) then
            call get_eta_sparse ( lo+sideband*Grids_x%frq_basis(f_indb:f_inda+1:-1), &
              & (/Frq/), eta_f(1:1,n_f:1:-1), not_zero_f(1:1,n_f:1:-1) )
          else
            call get_eta_sparse ( lo+sideband*Grids_x%frq_basis(f_inda+1:f_indb), &
              & (/Frq/), eta_f(1:1,1:n_f), not_zero_f(1:1,1:n_f) )
          end if
        end if

! Compute Sps_Path
        v_inda = grids_x%l_v(sps_i-1)
        ! Grids_X%Values are really 3-d: Frequencies X Zeta X Phi
        do sv_zp = w_inda + 1, w_indb
          do sv_f = 1, n_f
            v_inda = v_inda + 1
            if ( not_zero_f(1,sv_f) ) then
              where ( do_calc_zp(:,sv_zp) )
                eta_fzp(:,v_inda) = eta_f(1,sv_f) * eta_zp(:,sv_zp)
                sps_path(:,sps_i) = sps_path(:,sps_i) +  &
                                 &  grids_x%values(v_inda) * eta_fzp(:,v_inda)
                do_calc_fzp(:,v_inda) = Grids_x%deriv_flags(v_inda)
              elsewhere
                eta_fzp ( :, v_inda ) = 0.0_r8
                do_calc_fzp ( :, v_inda ) = .false.
              end where
            else
              eta_fzp ( :, v_inda ) = 0.0_r8
              do_calc_fzp ( :, v_inda ) = .false.
            end if
          end do ! sv_f
        end do ! sv_zp

        if ( grids_x%lin_log(sps_i)) sps_path(:,sps_i) = EXP(sps_path(:,sps_i))

      end if

      f_inda = f_indb
      w_inda = w_indb

    end do

  end subroutine Comp_Sps_Path_Frq

! ------------------------------------------------  Comp_Sps_Path  -----
  subroutine Comp_Sps_Path ( Grids_x, SPS_I, Eta_zp, Sps_path )

! Compute the SPS path for grids that use the Frequency dimension for
! a vector component, e.g. magnetic field.  There's only one multiply-add
! in the loop (hiding inside of matmul), so there's really no point to
! sending in do_calc_zp.

    use MLSCommon, only: RP
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
    integer, intent(in) :: SPS_I           ! Which thing-o in Grids_X
    real(rp), intent(in) :: Eta_zp(:,:)    ! Eta_z x Eta_phi for each path
!                         element and (ZxP). First dimension is same as sps_path.

! Output:

    real(rp), intent(out) :: Sps_Path(:,:) ! Values, path X component

! Local variables:

    integer :: N_C, SV_C, V_Inda, V_Indb

    n_c = grids_x%l_f(sps_i) - grids_x%l_f(sps_i-1) ! # Components

    v_inda = grids_x%l_v(sps_i-1)              ! One before first value
    v_indb = grids_x%l_v(sps_i)                ! Last value

    ! Grids_X%Values are really 3-d: Components X Zeta X Phi.
    ! For each component, multiply the Zeta X Phi part by Eta (on the left).
    ! We could do this as
    !  matmul( eta_zp, &
    !       &  transpose(reshape(grids_x%values,(/n_c,(v_indb-v_inda)/n_c/) )))
    ! but this is more efficient, not much more writing, and not that hard
    ! to grok.
    do sv_c = 1, n_c
      sps_path(:,sv_c) = matmul(eta_zp,grids_x%values(v_inda+sv_c:v_indb:n_c))
    end do

  end subroutine Comp_Sps_Path

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Comp_Sps_Path_Frq_m
!
! $Log$
! Revision 2.17  2004/07/08 21:00:23  vsnyder
! Inching toward PFA
!
! Revision 2.16  2003/07/08 02:01:31  vsnyder
! Speed up a tad by storing zero in ELSE instead of everywhere
!
! Revision 2.15  2003/05/16 02:46:33  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.14  2003/05/05 23:00:25  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.13.2.3  2003/03/22 02:32:05  vsnyder
! Polish the basic algorithm, add comp_sps_path
!
! Revision 2.13.2.2  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.13  2002/11/14 00:52:24  livesey
! Bug fix, changed arguments to intent(inout) and improved initialization
! steps. (Nathaniel and Bill worked on this one together).
!
! Revision 2.12  2002/10/08 17:08:01  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2002/09/27 20:43:06  livesey
! Bug fix for 'backwards' bases in get_eta_sparse for eta_f
!
! Revision 2.10  2002/09/06 20:58:26  vsnyder
! Cosmetic changes, copyright notice, move USEs to procedure scope
!
! Revision 2.9  2002/08/22 23:13:20  livesey
! New intermediate frequency based frq_bases
!
! Revision 2.8  2002/06/13 22:39:42  bill
! some variable name changes--wgr
!
! Revision 2.7  2002/06/04 10:28:00  zvi
! rename n_sps to: no_mol, more correctly
!
! Revision 2.6  2002/02/16 06:37:34  zvi
! New code for derivative flags..
!
! Revision 2.5  2002/01/09 00:30:48  zvi
! Fix a bug with skip_eta_frq
!
! Revision 2.4  2001/11/15 01:21:58  zvi
! Extiction debug fix
!
! Revision 2.3  2001/11/10 00:45:08  zvi
! Fixing a bug..
!
! Revision 2.2  2001/11/07 09:59:12  zvi
! More effective code for sps_path calculations
!
! Revision 1.0  2001/10/30 14:00:00 zvi Exp $"

