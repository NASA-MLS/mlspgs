! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Comp_Sps_Path_Frq_m

  implicit NONE

  private
  public :: Comp_Sps_Path_Frq, Comp_Sps_Path_Frq_nz, Comp_Sps_Path
  public :: Comp_Sps_Path_No_Frq, Comp_1_Sps_Path_No_Frq, Comp_One_Path_Frq

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

! --------------------------------------------  Comp_Sps_Path_Frq  -----
  subroutine Comp_Sps_Path_Frq ( Grids_x, Frq, eta_zp, &
    & do_calc_zp, sps_path, eta_fzp, do_calc_fzp, lo, sideband, &
    & Already_Spread )

! Compute the SPS path, including for species that are frequency dependent.

    use Comp_Eta_Docalc_No_Frq_m, only: Spread_Eta_FZP_from_Eta_ZP
    use MLSCommon, only: RP, IP, R8
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates

    real(r8), intent(in) :: Frq  ! Frequency at which to compute the values
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z x Eta_phi)
            ! First dimension is same as sps_values.
    logical, intent(in) :: Do_Calc_Zp(:,:) ! logical indicating whether there
            ! is a contribution for this state vector element

! Output:

!   These are intent(inout), not intent(out), so they won't be clobbered
!   in the cases when frq > 1.0 and n_f == 1.  Hopefully, the desired values
!   will have been filled by a call (outside the frequency loop) with frq <= 1.
    real(rp), intent(inout) :: Sps_Path(:,:) ! Path X Species.  vmr values
      !       along the path by species number

    real(rp), intent(inout) :: Eta_Fzp(:,:)  ! Path X (Eta_f x Eta_z x Eta_phi)
      !       First dimension is same as sps_values.
    logical, intent(inout) :: Do_Calc_Fzp(:,:) ! indicates whether
      !       there is a contribution for this state vector element. Same
      !       shape as Eta_Fzp.

! Optional inputs -- LO and Sideband not needed if Frq < 1.0:

    real(r8), intent(in), optional :: LO       ! Local oscillator freq
    integer, intent(in), optional :: Sideband  ! -1, 1 or 0.  Zero means
      !       Grids_x%frq_basis is absolute, not I.F.
    logical, intent(in), optional :: Already_Spread ! Eta_zp and Do_Calc_zp
      !       have already been spread out to Eta_fzp and Do_Calc_fzp

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    integer(ip) :: n_f, no_mol
    integer(ip) :: sps_i, sv_zp, sv_f
    integer(ip) :: f_inda, f_indb, v_inda, w_inda, w_indb

    real(rp) :: eta_f(1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))
    logical :: not_zero_f(1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))
    logical :: Spread

! Begin executable code:

    spread = .false.
    if ( present(already_spread) ) spread = already_spread

    if ( .not. spread ) &
      ! Spread out Eta_zp and Do_Calc_zp to Eta_fzp and Do_Calc_fzp,
      ! copying the same values for all frequencies for frequency-dependent
      ! species
      & call Spread_Eta_FZP_from_Eta_ZP ( grids_x, eta_zp, do_calc_zp, &
                                        &          eta_fzp, do_calc_fzp )

    no_mol = ubound(grids_x%l_z,1)

    if ( frq <= 1.0_r8 ) then
      do_calc_fzp = .FALSE.
      eta_fzp = 0.0_rp
    end if

    f_inda = 0
    w_inda = 0

    do sps_i = 1, no_mol

      f_indb = grids_x%l_f(sps_i)
      n_f = f_indb - f_inda

      w_indb = w_inda + (Grids_x%l_z(sps_i) - Grids_x%l_z(sps_i-1)) * &
                        (Grids_x%l_p(sps_i) - Grids_x%l_p(sps_i-1))
      if ( n_f /= 1 ) then

! Compute eta_f:
        if ( sideband == -1 ) then
          call get_eta_sparse ( lo-Grids_x%frq_basis(f_indb:f_inda+1:-1), &
            & Frq, eta_f(n_f:1:-1), not_zero_f(n_f:1:-1) )
        else if ( sideband == +1 ) then
          call get_eta_sparse ( lo+Grids_x%frq_basis(f_inda+1:f_indb), &
            & Frq, eta_f(1:n_f), not_zero_f(1:n_f) )
        else ! sideband == 0 means Grids_x%frq_basis is absolute, not in I.F.
             ! It doesn't mean folded-sideband calculation.
          call get_eta_sparse ( Grids_x%frq_basis(f_inda+1:f_indb), &
            & Frq, eta_f(1:n_f), not_zero_f(1:n_f) )
        end if

! Compute Sps_Path and fill in the parts of Eta_fzp that depend upon frequency.
        sps_path(:,sps_i) = 0.0_rp
        v_inda = grids_x%l_v(sps_i-1)
        ! Grids_X%Values are really 3-d: Frequencies X Zeta X Phi
        do sv_zp = w_inda + 1, w_indb
          do sv_f = 1, n_f
            v_inda = v_inda + 1
            if ( not_zero_f(sv_f) ) then
              ! Calculate eta_fzp from eta_zp with frequency interpolation
              eta_fzp(:,v_inda) = eta_f(sv_f) * eta_zp(:,sv_zp)
              sps_path(:,sps_i) = sps_path(:,sps_i) +  &
                               &  grids_x%values(v_inda) * eta_fzp(:,v_inda)
              do_calc_fzp(:,v_inda) = do_calc_zp(:,sv_zp) .and. Grids_x%deriv_flags(v_inda)
            else
              eta_fzp ( :, v_inda ) = 0.0_r8
              do_calc_fzp ( :, v_inda ) = .false.
            end if
          end do ! sv_f
        end do ! sv_zp

        if ( grids_x%lin_log(sps_i)) sps_path(:,sps_i) = EXP(sps_path(:,sps_i))

      else ! n_f == 1 here
        continue
        ! Hopefully, Comp_Sps_Path_No_Frq has already been called, and
        ! Sps_Path, Eta_Fzp and Do_Calc_Fzp have not been clobbered since then.
      end if

      f_inda = f_indb
      w_inda = w_indb

    end do

  end subroutine Comp_Sps_Path_Frq

! -----------------------------------------  Comp_Sps_Path_Frq_nz  -----
  subroutine Comp_Sps_Path_Frq_nz ( Grids_x, Frq, eta_zp, nz_zp, nnz_zp, &
    & do_calc_zp, sps_path, do_calc_fzp, eta_fzp, nz_fzp, nnz_fzp,       &
    & lo, sideband, Already_Spread )

! Compute the SPS path

    use Comp_Eta_Docalc_No_Frq_m, only: Spread_Eta_FZP_from_Eta_ZP
    use MLSCommon, only: RP, IP, R8
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates

    real(r8), intent(in) :: Frq  ! Frequency at which to compute the values
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z x Eta_phi)
!                         First dimension is same as sps_values.
    integer, intent(in) :: NZ_ZP(:,:)      ! Nonzeros of Eta_ZP
    integer, intent(in) :: NNZ_ZP(:)       ! Numbers of rows in NZ_ZP
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
    integer, intent(inout) :: NZ_FZP(:,:)    ! Nonzeros of Eta_FZP
    integer, intent(inout) :: NNZ_FZP(:)     ! Numbers of rows in NZ_FZP

! Optional inputs -- LO and Sideband not needed if Frq < 1.0:

    real(r8), intent(in), optional :: LO       ! Local oscillator freq
    integer, intent(in), optional :: Sideband  ! -1, 1 or 0.  Zero means
      !       Grids_x%frq_basis is absolute, not I.F.
    logical, intent(in), optional :: Already_Spread ! Eta_zp and Do_Calc_zp
      !       have already been spread out to Eta_fzp and Do_Calc_fzp

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    integer(ip) :: n_f, no_mol
    integer(ip) :: sps_i, sv_zp, sv_f
    integer(ip) :: v_inda, f_inda, f_indb, w_inda, w_indb

    real(rp) :: eta_f(1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))
    logical :: not_zero_f(1:maxval(grids_x%l_f(1:)-grids_x%l_f(0:ubound(grids_x%l_f,1)-1)))
    logical :: Spread

! Begin executable code:

    spread = .false.
    if ( present(already_spread) ) spread = already_spread

    if ( .not. spread ) &
      ! Spread out Eta_zp and Do_Calc_zp to Eta_fzp and Do_Calc_fzp,
      ! copying the same values for all frequencies for frequency-dependent
      ! species
      & call Spread_Eta_FZP_from_Eta_ZP ( grids_x, eta_zp, do_calc_zp, &
                                        &          eta_fzp, do_calc_fzp )

    ! Clear nonzero entries of eta_fzp
    do v_inda = 1, size(nnz_fzp)
      eta_fzp(nz_fzp(:nnz_fzp(v_inda),v_inda),v_inda) = 0.0
      do_calc_fzp(nz_fzp(:nnz_fzp(v_inda),v_inda),v_inda) = .false.
      nnz_fzp(v_inda) = 0
    end do

    no_mol = ubound(grids_x%l_z,1)

    f_inda = 0
    w_inda = 0

    do sps_i = 1, no_mol

      f_indb = grids_x%l_f(sps_i)
      n_f = f_indb - f_inda

      w_indb = w_inda + (Grids_x%l_z(sps_i) - Grids_x%l_z(sps_i-1)) * &
                        (Grids_x%l_p(sps_i) - Grids_x%l_p(sps_i-1))

      if ( n_f /= 1 ) then

! Compute eta_f:
        if ( sideband == -1 ) then
          call get_eta_sparse ( lo-Grids_x%frq_basis(f_indb:f_inda+1:-1), &
            & Frq, eta_f(n_f:1:-1), not_zero_f(n_f:1:-1) )
        else
          call get_eta_sparse ( lo+Grids_x%frq_basis(f_inda+1:f_indb), &
            & Frq, eta_f(1:n_f), not_zero_f(1:n_f) )
        end if

! Compute Sps_Path and fill in the parts of Eta_fzp that depend upon frequency.

        sps_path(:,sps_i) = 0.0_rp
        v_inda = grids_x%l_v(sps_i-1)
        ! Grids_X%Values are really 3-d: Frequencies X Zeta X Phi
        do sv_zp = w_inda + 1, w_indb
          do sv_f = 1, n_f
            v_inda = v_inda + 1
            if ( not_zero_f(sv_f) ) then
              nnz_fzp(v_inda) = nnz_zp(sv_zp)
              nz_fzp(:nnz_fzp(v_inda),v_inda) = nz_zp(:nnz_zp(sv_zp),sv_zp)
              eta_fzp(nz_fzp(:nnz_fzp(v_inda),v_inda),v_inda) = &
                & eta_f(sv_f) * eta_zp(nz_zp(:nnz_zp(sv_zp),sv_zp),sv_zp)
              sps_path(nz_fzp(:nnz_fzp(v_inda),v_inda),sps_i) = &
                & sps_path(nz_fzp(:nnz_fzp(v_inda),v_inda),sps_i) + &
                & grids_x%values(v_inda) * &
                & eta_fzp(nz_fzp(:nnz_fzp(v_inda),v_inda),v_inda)
              do_calc_fzp(nz_fzp(:nnz_fzp(v_inda),v_inda),v_inda) = &
                & do_calc_zp(nz_zp(:nnz_zp(sv_zp),sv_zp),sv_zp) .and. &
                & Grids_x%deriv_flags(v_inda)
            end if
          end do ! sv_f
        end do ! sv_zp

        if ( grids_x%lin_log(sps_i)) sps_path(:,sps_i) = EXP(sps_path(:,sps_i))

      end if

      f_inda = f_indb
      w_inda = w_indb

    end do

  end subroutine Comp_Sps_Path_Frq_nz

! --------------------------------------------  Comp_One_Path_Frq  -----
  subroutine Comp_One_Path_Frq ( Grids_x, Frq, eta_zp, &
    & do_calc_zp, path_qty, do_calc_fzp, eta_fzp, lo, sideband )

! Compute a one-quantity path value

    use MLSCommon, only: RP, IP, R8
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
!                         Assume it's a one-quantity grid

    real(r8), intent(in) :: Frq  ! Frequency at which to compute the values
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z x Eta_phi)
!                         First dimension is same as sps_values.
    logical, intent(in) :: Do_Calc_Zp(:,:) ! logical indicating whether there
!                         is a contribution for this state vector element

! Output:

    real(rp), intent(inout) :: Path_qty(:) ! Grids_x value interpolated to path.
    logical, intent(inout) :: Do_Calc_Fzp(:,:) ! indicates whether there
!                         is a contribution for this state vector element.
!                         Same shape as Eta_Fzp.
    real(rp), intent(inout) :: Eta_Fzp(:,:)  ! Path X (Eta_f x Eta_z x Eta_phi)
!                         First dimension is same as sps_values.

! Optional inputs -- not needed if Frq < 1.0:

    real(r8), intent(in), optional :: LO       ! Local oscillator freq
    integer, intent(in), optional :: Sideband  ! -1 or 1

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of path_qty = sps_values.

! Internal declarations

    integer(ip) :: n_f
    integer(ip) :: sv_zp, sv_f
    integer(ip) :: v_inda, w_indb

    real(rp) :: eta_f(1:grids_x%l_f(1))
    logical :: not_zero_f(1:grids_x%l_f(1))

! Begin executable code:

    if ( frq <= 1.0_r8 ) then
      do_calc_fzp = .FALSE.
      eta_fzp = 0.0_rp
    end if

    n_f = grids_x%l_f(1)

    w_indb = Grids_x%l_z(1) * Grids_x%l_p(1)

! Compute eta_f:
    if ( sideband == -1 ) then
      call get_eta_sparse ( lo-Grids_x%frq_basis(n_f:1:-1), &
        & Frq, eta_f(n_f:1:-1), not_zero_f(n_f:1:-1) )
    else if ( sideband == +1 ) then
      call get_eta_sparse ( lo+Grids_x%frq_basis(1:n_f), &
        & Frq, eta_f(1:n_f), not_zero_f(1:n_f) )
    else ! sideband == 0 means Grids_x%frq_basis is absolute, not in I.F.
         ! It doesn't mean folded-sideband calculation.
      call get_eta_sparse ( Grids_x%frq_basis(1:n_f), &
        & Frq, eta_f(1:n_f), not_zero_f(1:n_f) )
    end if

! Compute path_qty
    path_qty = 0.0_rp
    v_inda = grids_x%l_v(1-1)
    ! Grids_X%Values are really 3-d: Frequencies X Zeta X Phi
    do sv_zp = 1, w_indb
      do sv_f = 1, n_f
        v_inda = v_inda + 1
        if ( not_zero_f(sv_f) ) then
          eta_fzp(:,v_inda) = eta_f(sv_f) * eta_zp(:,sv_zp)
          path_qty = path_qty +  &
                           &  grids_x%values(v_inda) * eta_fzp(:,v_inda)
          do_calc_fzp(:,v_inda) = do_calc_zp(:,sv_zp) .and. Grids_x%deriv_flags(v_inda)
        else
          eta_fzp ( :, v_inda ) = 0.0_r8
          do_calc_fzp ( :, v_inda ) = .false.
        end if
      end do ! sv_f
    end do ! sv_zp

    if ( grids_x%lin_log(1)) path_qty = EXP(path_qty)

  end subroutine Comp_One_Path_Frq

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
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z X Eta_phi)
!                         First dimension is size(Sps_Path,1).

! Output:

    real(rp), intent(out) :: Sps_Path(:,:) ! Values, path X component

! Local variables:

    integer :: N_C, SV_C, V_Inda, V_Indb

    n_c = grids_x%l_f(sps_i) - grids_x%l_f(sps_i-1) ! # Components

    v_inda = grids_x%l_v(sps_i-1)          ! One before first value
    v_indb = grids_x%l_v(sps_i)            ! Last value

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

! -----------------------------------------  Comp_Sps_Path_No_Frq  -----
  subroutine Comp_Sps_Path_No_Frq ( Grids_x, Eta_zp, Sps_Path )

! Compute the SPS path for species that don't use frequency.
! Skip the ones that do.

    use MLSCommon, only: RP
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z X Eta_phi)
!                         First dimension is size(Sps_Path,1).

! Output:

    real(rp), intent(inout) :: Sps_Path(:,:) ! Path X Species.  vmr values
!                         along the path by species number

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    integer :: f_inda, f_indb
    integer :: no_mol, n_f
    integer :: sps_i
    integer :: v_inda, v_indb
    integer :: w_inda, w_indb
    integer :: sv_f, sv_zp

! Begin executable code:

    no_mol = ubound(grids_x%l_z,1)

    f_indb = 0
    f_indb = 0
    v_indb = 0
    w_indb = 0
    do sps_i = 1, no_mol

      ! We could use Comp_1_Sps_Path_No_Frq here, at a tiny increase in cost.

      f_inda = f_indb
      f_indb = grids_x%l_f(sps_i) ! Index of last one for sps_i
      v_inda = v_indb
      v_indb = grids_x%l_v(sps_i) ! Index of last one for sps_i
      w_inda = w_indb
      w_indb = w_inda + (Grids_x%l_z(sps_i) - Grids_x%l_z(sps_i-1)) * &
                        (Grids_x%l_p(sps_i) - Grids_x%l_p(sps_i-1))

      ! Compute Sps_Path only for sps that don't depend on frequency
      if ( f_indb-f_inda == 1 ) then
        sps_path(:,sps_i) = &
          & matmul(eta_zp(:,w_inda+1:w_indb), grids_x%values(v_inda+1:v_indb))

        if ( grids_x%lin_log(sps_i)) sps_path(:,sps_i) = exp(sps_path(:,sps_i))

      end if

    end do

  end subroutine Comp_Sps_Path_No_Frq

! ---------------------------------------  Comp_1_Sps_Path_No_Frq  -----
  subroutine Comp_1_Sps_Path_No_Frq ( Grids_x, The_Sps, Eta_zp, Sps_Path )

! Compute the SPS path for one species that doesn't use frequency.
! Actually, if you have an Eta_ZP prepared for an arbitrary grid in its
! first dimension (not necessary a line-of-sight path), this will put the
! species on that grid.  It's just a matrix multiply, but it knows how to
! dig up the species values from the Grids_X structure given The_Sps.

    use MLSCommon, only: RP
    use Load_sps_data_m, only: Grids_T

! Input:

    type (grids_t), intent(in) :: Grids_x  ! All the needed coordinates
    integer, intent(in) :: The_Sps ! The molecule for which to compute the
!                         path value.
    real(rp), intent(in) :: Eta_zp(:,:)    ! Path X (Eta_z X Eta_phi)
!                         First dimension is size(Sps_Path).

! Output:

    real(rp), intent(inout) :: Sps_Path(:) ! Path.  vmr values
!                         along the path for species The_Sps

! Notes:
! units of z_basis must be same as zeta_path (usually -log(P)) and units of
! phi_basis must be the same as phi_path (either radians or degrees).
! Units of sps_path = sps_values.

! Internal declarations

    integer :: v_inda, v_indb, w_inda, w_indb

! Begin executable code:

    v_inda = grids_x%l_v(the_sps-1) ! Index of last one for the_sps - 1
    v_indb = grids_x%l_v(the_sps)   ! Index of last one for the_sps
    w_inda = w_inda + (Grids_x%l_z(the_sps-1) - Grids_x%l_z(the_sps-1-1)) * &
                      (Grids_x%l_p(the_sps-1) - Grids_x%l_p(the_sps-1-1))
    w_indb = w_inda + (Grids_x%l_z(the_sps) - Grids_x%l_z(the_sps-1)) * &
                      (Grids_x%l_p(the_sps) - Grids_x%l_p(the_sps-1))

    sps_path = &
      & matmul(eta_zp(:,w_inda+1:w_indb), grids_x%values(v_inda+1:v_indb))

    if ( grids_x%lin_log(the_sps)) sps_path = exp(sps_path)

  end subroutine Comp_1_Sps_Path_No_Frq

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Comp_Sps_Path_Frq_m
!
! $Log$
! Revision 2.35  2016/11/17 01:46:45  vsnyder
! Remove feature to calculate without frequency interpolation in the
! routines designed to do it, if Frq < 1.0.  Use the routines that don't
! do frequency interpolation instead.
!
! Revision 2.34  2013/06/12 02:24:35  vsnyder
! Cruft removal
!
! Revision 2.33  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.32  2011/08/20 00:44:58  vsnyder
! Remove declarations for unused variables
!
! Revision 2.31  2011/08/12 18:57:06  vsnyder
! Change Do_Calc_Fzp back to non-optional.  The only place where it wasn't
! used, it was actually needed.  Undo ill-advised change to evaluate the
! outputs when frq>1 and n_f==1.  The outputs were evaluated when frq<=1,
! when Comp_Sps_Path was called from outside the frequency loop.  Doing it
! again increases run time substantially, for no benefit.  Comments have
! been added to this effect.
!
! Revision 2.30  2011/07/29 01:56:15  vsnyder
! Make do_calc_fzp optional
!
! Revision 2.29  2010/09/25 01:11:20  vsnyder
! Interpret sideband==0 in Comp_Sps_Path_Frq, add Comp_One_Path_Frq
!
! Revision 2.28  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.27  2009/01/16 23:47:00  vsnyder
! Give warning in log-lin case if log(some mixing ratio) is zero
!
! Revision 2.26  2007/06/26 00:38:17  vsnyder
! Separate frq>1 and frq<1 cases
!
! Revision 2.25  2007/06/06 01:13:52  vsnyder
! Use 1-d eta routine
!
! Revision 2.24  2007/01/19 02:37:30  vsnyder
! In Comp_Sps_Path_No_Frq compute only for species that have no freq coordinate
!
! Revision 2.23  2005/12/22 20:51:46  vsnyder
! Added Comp_1_Sps_Path_No_Frq
!
! Revision 2.22  2005/09/16 23:41:19  vsnyder
! Cannonball polishing
!
! Revision 2.21  2005/08/03 18:04:09  vsnyder
! Some spectroscopy derivative stuff
!
! Revision 2.20  2005/06/22 18:08:18  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.19  2004/11/01 20:20:46  vsnyder
! Reorganization of representation for molecules and beta groups; PFA broken for now
!
! Revision 2.18  2004/08/03 21:59:49  vsnyder
! Inching further toward PFA
!
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

