! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Comp_Sps_Path_Frq_m

  implicit NONE

  private
  public :: comp_sps_path_frq

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!-----------------------------------------------------------------
 contains
!-----------------------------------------------------------------

! --------------------------------------------  Comp_Sps_Path_Frq  -----
  subroutine Comp_Sps_Path_Frq ( Grids_x, lo, sideband, Frq, eta_zp, &
    & do_calc_zp, sps_path, do_calc_fzp, eta_fzp )

! Compute the SPS path

    use MLSCommon, only: RP, IP, R8
    use get_eta_matrix_m, only: Get_Eta_Sparse
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
    real(r8), intent(in) :: LO            ! Local oscillator freq
    integer, intent(in) :: Sideband       ! -1 or 1

    real(r8), intent(in) :: Frq  ! Frequency at which to compute the values
    real(rp), intent(in) :: Eta_zp(:,:) ! Eta_z x Eta_phi for each state
!                         vector element. This is the same length as sps_values.
    logical, intent(in) :: Do_Calc_Zp(:,:) !logical indicating whether there
!                        is a contribution for this state vector element

! Output:

    real(rp), intent(out) :: Sps_Path(:,:) ! vmr values along the path
!                          by species number
    real(rp), intent(out) :: Eta_Fzp(:,:) ! Eta_z x Eta_phi x Eta_f for each
!                          state vector element. This is the same length as
!                          sps_values.
    logical, intent(out) :: Do_Calc_Fzp(:,:) ! indicates whether there
!                         is a contribution for this state vector element
!                         This is the same length as sps_values.
! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    integer(ip) :: n_zp, n_f, nfzp, f_len
    integer(ip) :: sps_i, no_mol, sv_i, sv_zp, sv_f, sv_j
    integer(ip) :: v_inda, v_indb, f_inda, f_indb, w_inda, w_indb

    real(rp) :: eta_f(1:1,1:maxval(grids_x%no_f))
    logical :: not_zero_f(1:1,1:maxval(grids_x%no_f))

! Begin executable code:

    no_mol = SIZE(Grids_x%no_z)

    if ( frq < 1.0_r8 ) then
      eta_fzp = 0.0_rp
      sps_path = 0.0_rp
      do_calc_fzp = .FALSE.
    else
      do sps_i = 1, no_mol
        if ( grids_x%no_f(sps_i) > 1 ) sps_path(:,sps_i) = 0.0_rp
      end do
    end if

    f_len = 0
    f_inda = 1
    v_inda = 1
    w_inda = 1

    do sps_i = 1, no_mol

      n_f = Grids_x%no_f(sps_i)
      n_zp = Grids_x%no_z(sps_i) * Grids_x%no_p(sps_i)
      nfzp = n_f * n_zp

      f_indb = f_inda + n_f
      v_indb = v_inda + nfzp
      w_indb = w_inda + n_zp

      if ( (frq > 1.0) .AND. grids_x%no_f(sps_i) == 1 ) then
        f_len = f_len + nfzp
      else

! There are two ways to do this (slow and easy vs quick but difficult)
! For ease lets do the slow and easy (and certainly more reliable)

! Compute eta:

        call get_eta_sparse ( lo+sideband*Grids_x%frq_basis(f_inda:f_indb-1), &
                            & (/Frq/), eta_f(1:1,1:n_f), not_zero_f(1:1,1:n_f) )

        do sv_i = 0 , nfzp - 1
          f_len = f_len + 1
          sv_j = v_inda + sv_i
          sv_f = 1 + MODULO(sv_i,n_f)
          sv_zp = w_inda + sv_i / n_f
          if ( not_zero_f(1,sv_f) ) then
            where ( do_calc_zp(:,sv_zp) )
              eta_fzp(:,sv_j) = eta_f(1,sv_f) * eta_zp(:,sv_zp)
              sps_path(:,sps_i) = sps_path(:,sps_i) +  &
                               &  grids_x%values(sv_j) * eta_fzp(:,sv_j)
            end where
            if ( Grids_x%deriv_flags(f_len) ) then
              where ( do_calc_zp(:,sv_zp) ) do_calc_fzp(:,sv_j) = .TRUE.
            end if
          end if
        end do

        if ( grids_x%lin_log(sps_i)) sps_path(:,sps_i) = EXP(sps_path(:,sps_i))

      end if

      f_inda = f_indb
      v_inda = v_indb
      w_inda = w_indb

    end do

  end subroutine Comp_Sps_Path_Frq

end module Comp_Sps_Path_Frq_m
!
! $Log$
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

