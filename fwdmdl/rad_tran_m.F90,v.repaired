! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Rad_Tran_m

  implicit NONE
  private
  public :: Rad_Tran, Rad_Tran_Pol
  public :: dRad_Tran_dF, dRad_Tran_dT, dRad_Tran_dX
  public :: Get_Do_Calc

  public :: Get_Do_Calc_Indexed, Get_Inds ! Used only by Hessians_m

  private :: Get_d_Delta_df, Get_d_Delta_df_f, Get_d_Delta_df_linlog
  private :: Get_d_Delta_df_linlog_f
  private :: N_Eta_Rows

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!------------------------------------------------------  Rad_Tran  -----
! This is the radiative transfer model, radiances only !

  subroutine Rad_Tran ( tan_pt, gl_inds, more_inds, e_rflty, del_zeta, &
                      & alpha_path, ref_cor, incoptdepth, ds_dz_gw,    &
                      & t_script, tau, inc_rad_path, rad, i_stop )

    use GL_Update_Incoptdepth_m
    use MLSKinds, only: RP
    use SCRT_dN_m, ONLY: SCRT

  ! Inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    integer, intent(in) :: more_inds(:)      ! Places in the coarse path
                                             ! where GL is needed
    real(rp), intent(in) :: e_rflty          ! earth reflectivity value (0--1).
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just the part up
                                             ! to the black-out
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
                                             ! length ratios.
    real(rp), intent(inout) :: incoptdepth(:) ! incremental path opacities
                                             ! from one-sided layer calculation
                                             ! on output. it is the full
                                             ! integrated layer opacity.
    real(rp), intent(in), target :: alpha_path(:) ! absorption coefficient on
                                             ! coarse & fine grid.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative * gw.
    real(rp), intent(in) :: t_script(:)      ! differential temperatures (K)
                                             ! on coarse grid.

  ! Outputs

    real(rp), intent(out) :: tau(:)          ! transmission function.
    real(rp), intent(out) :: inc_rad_path(:) ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(out) :: rad             ! radiance (K)
    integer, intent(out) :: i_stop           ! path stop index

  ! Begin code

  ! See if anything needs to be gl-d

    call GL_update_incoptdepth ( gl_inds, more_inds, 1, del_zeta, alpha_path, &
                               & ds_dz_gw, ref_cor, incoptdepth )

    call scrt ( tan_pt, t_script, e_rflty, incoptdepth, tau, rad, inc_rad_path, &
      &         i_stop )

  end subroutine Rad_Tran

!--------------------------------------------------  Rad_Tran_Pol  -----

  subroutine Rad_Tran_Pol ( tan_pt, gl_inds, more_inds, e_rflty, del_zeta,    &
                          & alpha_path, ref_cor, incoptdepth_pol, deltau_pol, &
                          & ds_dz_gw, ct, stcp, stsp, t_script, do_dumps,     &
                          & prod_pol, tau_pol, rad_pol, p_stop )

    ! Polarized radiative transfer.  Radiances only, no derivatives.

    use CS_Expmat_m, only: CS_Expmat
    use Dump_0, only: Dump, Dump_2x2xN
    use GLNP, ONLY: NG, NGP1
    use MCRT_M, ONLY: MCRT
    use MLSKinds, only: RP
    use Opacity_m, only: Opacity

  ! Inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: Gl_inds(:)        ! Gauss-Legendre grid indices
    integer, intent(in) :: More_inds(:)      ! Places in the coarse path
                                             ! where GL is needed
    real(rp), intent(in) :: E_rflty          ! earth reflectivity value (0--1).
    real(rp), intent(in) :: Del_zeta(:)      ! path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just part up to
                                             ! the black-out
    complex(rp), intent(in), target :: Alpha_path(-1:,:)  ! absorption
                                             ! coefficient on composite coarse &
                                             ! fine grid.
    complex(rp), intent(inout) :: Deltau_pol(:,:,:) ! 2 X 2 X path.  Incremental
                                             ! transmissivity on the coarse path.
                                             ! Called E in some notes.
    real(rp), intent(in) :: Ref_cor(:)       ! refracted to unrefracted path
                                             ! length ratios.
    complex(rp), intent(inout) :: Incoptdepth_pol(:,:,:) ! incremental path
                                             ! opacities from one-sided layer
                                             ! calculation on output. it is the
                                             ! full integrated layer opacity.
                                             ! 2x2xPath
    real(rp), intent(in) :: ds_dz_Gw(:)      ! path length wrt zeta derivative *
                                             ! gw on the entire grid.  Only the
                                             ! gl_inds part is used.
    real(rp), intent(in) :: CT(:)            ! Cos theta          for Mag field
    real(rp), intent(in) :: STCP(:)          ! Sin theta Cos Phi  for Mag field
    real(rp), intent(in) :: STSP(:)          ! Sin theta Sin Phi  for Mag field
    real(rp), intent(in) :: T_script(:)      ! differential temperatures (K)
                                             ! on coarse path.
    integer, intent(in) :: Do_Dumps          ! Dump intermediate results if > 0

  ! Outputs

    complex(rp), intent(out) :: Prod_pol(:,:,:) ! product of E matrices. 2x2xPath
    complex(rp), intent(out) :: Tau_pol(:,:,:)  ! transmission function. 2x2xPath
    complex(rp), intent(out) :: Rad_pol(:,:)    ! radiance (K). 2x2.
    integer, intent(out) :: P_Stop           ! path stop index if >= 0, else
                                             ! -index in incoptdepth_pol where
                                             ! cs_expmat failed.

  ! Internals

    integer :: A, AA
    complex(rp), pointer :: Alpha_Path_c(:,:)
    real(rp), save :: E_Stop  = 1.0_rp ! X for which Exp(X) is too small to worry
    complex(rp) :: gl_Delta_Polarized(-1:1,size(gl_inds)/ng)
    complex(rp) :: Incoptdepth_Pol_gl(2,2,size(gl_inds)/ng)
    integer :: I, J
    integer :: N_Path
    integer :: Status ! from cs_expmat

  ! Begin code

    n_path = size(del_zeta)
    alpha_path_c(-1:,1:) => alpha_path ( -1:1, 1::ngp1 )

    if ( e_stop > 0.0_rp ) e_stop = log(epsilon(0.0_rp)) ! only once

  ! See if anything needs to be gl-d

    if ( size(gl_inds) > 0 ) then

      !{ Apply Gauss-Legendre quadrature to the panels indicated by
      !  {\tt More\_inds}.  We remove a singularity (which actually only
      !  occurs at the tangent point) by writing
      !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
      !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
      !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
      !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
      !   \text{d}\zeta$.  The first integral is easy -- it's just
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  We don't use it here.
      !  In the second integral, $G(\zeta)$ is {\tt funct} -- which has
      !  already been evaluated at the appropriate abscissae -- and
      !  $G(\zeta_i)$ is {\tt singularity}.  The weights are {\tt gw}.

      a = 1
      do i = 1, size(more_inds)
        aa = gl_inds(a)
        do j = -1, 1
          gl_delta_polarized(j,i) = del_zeta(more_inds(i)) *                             &
                 &  sum( ( alpha_path(j,aa:aa+ng-1) - &
                 &         alpha_path_c(j,more_inds(i))) * &
                 &       ds_dz_gw(aa:aa+ng-1) )
        end do 
        a = a + ng
      end do

      ! Turn sigma-, pi, sigma+ GL corrections into 2X2 matrix of
      ! GL corrections to incoptdepth_pol
      call opacity ( ct(more_inds), stcp(more_inds), stsp(more_inds), &
        & gl_delta_polarized, incoptdepth_pol_gl )

      ! add GL corrections to incoptdepth_pol
      incoptdepth_pol(:,:,more_inds) = incoptdepth_pol(:,:,more_inds) - &
        & incoptdepth_pol_gl

      if ( do_dumps > 1 ) then
        call dump ( more_inds, name='More_Inds' )
        call dump ( gl_inds, name='GL_Inds' )
        call dump ( gl_delta_polarized, name='GL_Delta_Polarized', width=3, &
          & options='p' ) ! transpose=.TRUE.
        call Dump_2x2xN ( incoptdepth_pol_gl, name='IncoptDepth_Pol_GL' )
        call Dump_2x2xN ( incoptdepth_pol, name='Incoptdepth_Pol' )
        call dump ( ref_cor, name='Ref_Cor' )
        if ( do_dumps > 2 ) stop
      end if
    end if

    ! At this point, incoptdepth_pol(:,:,1:tan_pt_c) should be nearly
    ! identical to incoptdepth_pol(:,:,1:tan_pt_c+1) (tan_pt_c is the
    ! zero-thickness tangent layer).

    do p_stop = 0, n_path-1
      incoptdepth_pol(:,:,p_stop+1) = incoptdepth_pol(:,:,p_stop+1) * &
        &                             ref_cor(p_stop+1)
      ! exp(A) = exp(s) * ((sinh d)/d (A - s I) + cosh d I) where
      ! s is the sum of A's eigenvalues and d is their difference.
      ! (sinh d)/d and cosh d can be large and positive even when
      ! s is large and negative, so we can't stop just because s
      ! is large and negative.  Well, we could if we knew d was small,
      ! but once we have both eigenvalues we've almost finished the
      ! exponential anyway.
      call cs_expmat ( incoptdepth_pol(:,:,p_stop+1), &
        &              deltau_pol(:,:,p_stop+1), status )  
      if ( status /= 0 ) go to 99 ! because we can't change p_stop in the loop
    end do

    call mcrt ( t_script, sqrt(e_rflty), deltau_pol, &
      & tan_pt, p_stop, prod_pol, tau_pol, rad_pol )

    return

  ! Error exit if cs_expmat detected an overflow
 99 if ( do_dumps > 0 ) then
      call dump ( more_inds, name='More_Inds' )
      call dump ( gl_inds, name='GL_Inds' )
      call dump ( gl_delta_polarized, name='GL_Delta_Polarized', width=3 )
      call Dump_2x2xN ( incoptdepth_pol_gl, name='IncoptDepth_Pol_GL' )
      call Dump_2x2xN ( incoptdepth_pol(:,:,:p_stop), name='Incoptdepth_Pol' )
      call dump ( ref_cor, name='Ref_Cor' )
      if ( do_dumps > 2 ) stop
    end if
    p_stop = - p_stop - 1

  end subroutine Rad_Tran_Pol

!--------------------------------------------------  DRad_Tran_df  -----
! This is the radiative transfer derivative wrt mixing ratio model

  subroutine DRad_Tran_df ( gl_inds, del_zeta, Grids_f, Eta_FZP, do_gl,    &
                          & del_s, ref_cor, ds_dz_gw, inc_rad_path,        &
                          & dAlpha_df, i_start, tan_pt, i_stop, drad_df,   &
                          & dB_df,                                         &
                          ! Optional for polarized
                          & Sparse_d_delta_df,                             &
                          ! Optionals for TScat
                          & Tau, alpha_path_c, Beta_c_e, dBeta_c_a_dIWC,   &
                          & dBeta_c_s_dIWC, dTScat_df, W0 )

    use d_T_Script_dTnp_m, only: dT_Script
    use GLNP, only: NGP1
    use Load_SPS_Data_m, ONLY: Grids_t
    use MLSKinds, only: RP
    use SCRT_dN_m, ONLY: dSCRT_dT, dSCRT_dX
    use Sparse_m, only: Sparse_t
    use TScat_Support_m, only: Get_dB_df

! Inputs

    integer, intent(in) :: GL_Inds(:)        ! Gauss-Legendre grid indices
    real(rp), intent(in) :: Del_Zeta(:)      ! path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just the part up
                                             ! to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    class(sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                             ! from state vector to combined
                                             ! coarse & fine path for each sps
    logical, intent(in) :: Do_GL(:)          ! A logical indicating where to
                                             ! do gl integrations
    real(rp), intent(in) :: Ref_Cor(:)       ! refracted to unrefracted path
                                             ! length ratios.
    real(rp), intent(in) :: Del_s(:)         ! unrefracted coarse path panel
                                             ! length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
                                             ! gw on the entire grid.  Only the
                                             ! gl_inds part is used.
    real(rp), intent(in) :: Inc_Rad_Path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(in) :: dAlpha_df(:,:) ! On the
                                             ! composite coarse & GL path
    integer, intent(in) :: I_start           ! path_start_index + 1 in coarse path
    integer, intent(in) :: Tan_Pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: I_stop            ! path stop index in coarse path

    ! Optional for polarized
    class(sparse_t), intent(inout), optional :: Sparse_d_delta_df(:)

    ! Optionals for TScat
    real(rp), intent(inout) :: dB_df(:) ! scratch, on the path, size=0 for no TScat
    real(rp), intent(in), optional :: Tau(:)
    real(rp), intent(in), optional :: Alpha_path_c(:)
    real(rp), intent(in), optional :: Beta_c_e(:)
    real(rp), intent(in), optional :: dBeta_c_a_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in), optional :: dBeta_c_s_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in), optional :: dTScat_df(:,:) ! Path X SV.  On the path, w.r.t.
                                                 ! f == Molecules on the grid
    real(rp), intent(in), optional :: W0(:)

! Outputs

    real(rp), intent(out) :: drad_df(:)         ! derivative of radiances wrt
                                                ! mixing ratio state vector
                                                ! element. (K)

! Internals

    integer, pointer :: All_inds(:)   ! all_inds => part of all_inds_B;
                                      ! Indices on GL grid for stuff
                                      ! used to make GL corrections
    integer, target :: All_inds_B(1:size(del_s))
    real(rp) :: d_Delta_B_df(size(dB_df))
    real(rp) :: d_Delta_df(1:size(del_s)) ! 1:max_c
    logical :: Do_Calc(1:size(del_s)) ! Flags on coarse path where do_calc_fzp
                                      ! for the coarse path or (do_gl and any
                                      ! corresponding do_calc_fzp on the fine
                                      ! path).
    logical :: Do_Calc_FZP(n_eta_rows(eta_fzp,.true.))
    logical :: Do_TScat                   ! Include dependence upon dB_df
    real(rp) :: Eta_FZP_Col(n_eta_rows(eta_fzp,.true.))
    integer, pointer :: Inds(:)       ! inds => part_of_nz_d_delta_df;
                                      ! Indices on coarse path where do_calc.
    integer :: I, I_Begin, SPS_i, SV_i
    ! NZ_ZXP and NNZ are used to project dB_df from the path to the grid
    integer, target :: More_inds_B(1:size(del_s))
    integer, pointer :: More_inds(:)  ! more_inds => part of more_inds_B;
                                      ! Indices on the coarse path where GL
                                      ! corrections get applied.
    integer :: NNZ_d_Delta_df
    integer :: NNZ_ZXP ! number of elements used in NZ_ZPX
    integer :: No_To_GL
    integer, target :: NZ_d_Delta_df(1:size(d_Delta_df))
    integer :: NZ_ZXP(n_eta_rows(eta_fzp,.true.))
    logical :: Save_Sparse            ! Save d_delta_df in sparse_d_delta_df
    real(rp) :: Singularity(1:size(del_s)) ! integrand on left edge of coarse
                                      ! grid panel -- singular at tangent pt.
    integer :: Sparse_Col             ! Column index in Eta_FZP(sps_i)

! Begin code

    save_sparse = present(sparse_d_delta_df)
    if ( save_sparse ) save_sparse = size(sparse_d_delta_df) > 0

    ! Get pointer to coarse grid subset of dAlpha_df

    d_delta_df = 0
    do_calc_fzp = .false.
    do_TScat = size(dB_df) > 0
    eta_fzp_col = 0
    nnz_d_delta_df = 0
    nz_d_delta_df = 0

    do sps_i = 1, ubound(Grids_f%l_z,1)
      sparse_col = 0
      if ( save_sparse ) call sparse_d_delta_df(sps_i)%empty

      if ( Do_TScat ) then

        !{ $\frac{\partial I}{\partial f^\text{iwc}_{lm}} =
        !   \sum_{i=1}^{N_p}
        !   \frac{\partial \overline{\Delta B_i}}{\partial f^\text{iwc}_{lm}}
        !   \mathcal{T}_i + \overline{\Delta B_i}
        !   \frac{\partial \mathcal{T}_i}{\partial f^\text{iwc}_{lm}}$
        !   where $\overline{\Delta B} = \Delta B^g + \Delta B^s$,
        !   $\Delta B^g = \Delta \left[ ( 1-\omega_0 ) B \right]$,
        !   $\Delta B^s = \Delta \left[ \omega_0 T_\text{scat} \right]$
        ! and
        !   $\frac{\partial \mathcal{T}_i}{\partial f^k_{lm}} =
        !    -\mathcal{T}_i \int_{s_0}^{s_m}
        !      \frac{\partial \alpha(\sigma)}{\partial f^k_{lm}}
        !       \, \text{d} \sigma \approx
        !    -\mathcal{T}_i \sum_{j=i}^{N_p}
        !      \frac{\partial \delta^k_j}{\partial f^k_{lm}}
        !   = -\mathcal{T}_i \sum_{j=i}^{N_p} \beta^k_j
        !       \eta^k_{lm}(s_i) \Delta s_j$.
        !
        !   {\tt dB_df(i)} =
        !   $\frac{\partial \overline{\Delta B_i}}{\partial f^k_{lm}}
        !    \text{ where } \frac{\partial \overline{B_i}}
        !                        {\partial f^k_{lm}(\zeta_i)} =
        !    \frac{\partial \omega_{0_i}}{\partial f^k_{lm}(\zeta_i)}
        !     \left( T_{\text{scat}_i} - B_i \right) +
        !    \omega_{0_i} \frac{\partial T_{\text{scat}_i}}
        !                      {\partial f^k_{lm}(\zeta_i)}
        !   = \left(\frac{\partial \omega_{0_i}}{\partial f^k}
        !     \left( T_{\text{scat}_i} - B_i \right) +
        !     \omega_{0_i} \frac{\partial T_{\text{scat}_i}}{\partial f^k}
        !     \right)
        !     \frac{\partial f^k}{\partial f^k_{lm}(\zeta_i)}
        !   = \left( \frac{\partial \omega_{0_i}}{\partial f^k}
        !      \left( T_{\text{scat}_i} - B_i \right) +
        !     \omega_{0_i} \frac{\partial T_{\text{scat}_i}}{\partial f^k}
        !      \right) \eta^k_{lm}(\zeta_i)$ (see wvs-095 and {\tt Get_dB_df}.
        !
        ! {\tt inc_rad_path(i)} = $\mathcal{T}_i \overline{\Delta B}_i$

        call get_dB_df ( alpha_path_c, beta_c_e, dBeta_c_a_dIWC,       &
                       & dBeta_c_s_dIWC, dAlpha_df(1::ngp1,sps_i), w0, &
                       & grids_f%mol(sps_i), dB_df )

      end if

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)
        sparse_col = sparse_col + 1

        ! We keep track of where we create nonzeros in d_delta_df, and replace
        ! them by zeros on the next iteration.  This is done because the vast
        ! majority of d_delta_df elements are zero, and setting all of
        ! d_Delta_df was found to be a significant expense.

        ! Everything in d_delta_df not indexed by nz_d_delta_df is already zero

        d_delta_df(nz_d_delta_df(:nnz_d_delta_df)) = 0
        nnz_d_delta_df = 0

        ! Skip the masked derivatives, according to the l2cf inputs

        drad_df(sv_i) = 0.0
        if ( .not. Grids_f%deriv_flags(sv_i) ) cycle

        ! Get the interpolating coefficients along the path
        call eta_fzp(sps_i)%get_flags ( sparse_col, do_calc_fzp )
        call get_do_calc_indexed ( size(do_gl), tan_pt, do_calc_fzp, gl_inds, &
          & do_gl, do_calc, nnz_d_delta_df, nz_d_delta_df )
        call eta_fzp(sps_i)%clear_flags ( sparse_col, do_calc_fzp )
        if ( nnz_d_delta_df == 0 ) cycle

        inds => nz_d_delta_df(1:nnz_d_delta_df)

        if ( save_sparse ) then
          do i = 1, nnz_d_delta_df
            call sparse_d_delta_df(sps_i)%add_element ( d_delta_df(inds(i)), &
              & inds(i), sparse_col )
          end do
        end if

        no_to_gl = count(do_gl(inds))

        all_inds => all_inds_B(1:no_to_gl)
        more_inds => more_inds_B(1:no_to_gl)

        ! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

        !{ Get d_Delta_df for one state-vector element.  This is
        !  $\frac{\partial \delta_i}{\partial f^k_{lm}}$ where $i$ is the
        !  index of a path point, $k$ is the species index, and $lm$ are indices
        !  for $(\phi_l,\zeta_m)$.  {\tt sv_i} flattens $(k,l,m)$ to one index.
        if ( eta_fzp(sps_i)%cols(sparse_col) /= 0 ) then ! process only non-empty columns
          call eta_fzp(sps_i)%get_col ( sparse_col, eta_fzp_col )
          call get_d_delta_df ( inds, gl_inds, all_inds, more_inds,     &
            & eta_fzp_col, dAlpha_df(:,sps_i), del_s, del_zeta, ds_dz_gw,        &
            & singularity, d_delta_df, ref_cor, grids_f%lin_log(sps_i), grids_f%values(sv_i) )
          call eta_fzp(sps_i)%clear_col ( sparse_col, eta_fzp_col )
        end if

        i_begin = max(i_start,min(nz_d_delta_df(1),i_stop))

        !{ $I(s_m) = \mathcal{T}(s_0,s_m) \left( I(s_0)-B(s_0) \right) +
        !  B(s_m) - \int_{B(s_0)}^{B(s_m)} \mathcal{T}(s,s_m)
        !            \frac{\partial B(s)}{\partial s} \,\text{d} s$ with
        !  $\mathcal{T}(s,s_m) =
        !    \exp\left( -\int_s^{s_m} \alpha(\sigma) \,\text{d} \sigma \right)$,
        !  $s_0$ is the end of the path away from the instrument,
        !  $s_m$ is the location of the instrument, $\beta^k(\sigma)$ is the
        !  absorption coefficient for the $k^\text{th}$ species, and
        !  $\eta^k_{lm}(\sigma)$ is an interpolation coefficient from
        !  $(\phi^k_l,\zeta^k_m)$ to $\sigma$. $(\phi^k_l,\zeta^k_m)$ is
        !  specified by {\tt sv_i}.
        !  
        !  The integral is approximated by
        !  $\sum_{i=1}^{N_p} \mathcal{T}_i \Delta B_i$ where $\Delta B_i =
        !  \frac12 \left( B_{i+1} - B_{i-1} \right)$.  The end-point terms
        !  are incorporated into special values of $\Delta B_0$ and $\Delta B_{s_m}$.
        !
        !  {\tt inc_rad_path(i)} = $\mathcal{T}_i \Delta B_i$.
        !  $\frac{\partial \mathcal{T}}{\partial f^k_{lm}} =
        !  -\mathcal{T} \int_s^{s_m} \frac{\partial \alpha(s)}{\partial f^k_{lm}}.$
        !  {\tt d_delta_df} =
        !  $ \int_s^{s_m} \frac{\partial \alpha(s)}{\partial f^k_{lm}}
        !   \,\text{d} s$

        if ( Do_TScat ) then

          !{ We got $\frac{\partial \omega_{0_i}}{\partial f^k_i}
          !     \left( T_{\text{scat}_i} - B_i \right)$ from {\tt Get_dB_df}.
          !  Now to finish we need to multiply by
          !  $\frac{\partial f^k_i}{\partial f^k_{lm}} = \eta^k_{lm}(\zeta_i)$
          !  and add $\omega_{0_i}
          !  \frac{\partial T_{\text{scat}_i}}{\partial f^k_{lm}}$ to finish
          !  computing $\frac{\partial \overline{B}_i}{\partial f^k_{lm}}$, and
          !  then compute $\Delta \frac{\partial \overline{B}_i}{\partial f^k_{lm}}$
          !  from that.

          nnz_zxp = 0
          if ( eta_fzp(sps_i)%cols(sparse_col) /= 0 ) then ! process only non-empty columns
            ! Copy nonzero elements from Eta_FZP to Eta_FZP_Col, and fill
            ! the sparsity indicators NZ_ZXP and NNZ_ZXP.  We need to
            ! do this because dT_Script doesn't understand Sparse_t.
            call eta_fzp(sps_i)%get_col ( sparse_col, eta_fzp_col, nnz_zxp, nz_zxp )
          end if
          call dt_script ( dB_df, eta_fzp_col, nz_zxp, nnz_zxp, d_delta_B_df, &
            & w0, dTScat_df(:,sv_i) )
          ! Now put zeroes into Eta_FZP_Col were there were nonzeros
          eta_fzp_col(nz_zxp(:nnz_zxp)) = 0

          !{ Now do the path integration
          !  $\sum \frac{\partial}{\partial f^k_{lm}}
          !     \mathcal{T}_i \overline{\Delta B_i}$

          call dscrt_dt ( tan_pt, d_delta_df, tau, inc_rad_path,&
                        & d_delta_B_df, i_begin, i_stop, drad_df(sv_i) )
        else
          call dscrt_dx ( tan_pt, d_delta_df, inc_rad_path, &
                       &  i_begin, i_stop, drad_df(sv_i))
        end if

      end do ! sv_i

    end do ! sps_i

  end subroutine DRad_Tran_df

!--------------------------------------------------  DRad_Tran_dT  -----
! This is the radiative transfer derivative wrt temperature model

  subroutine DRad_Tran_dT ( gl_inds, del_zeta, h_path_c, dh_dt_path,          &
                         &  alpha_path, dAlpha_dT_path, eta_zxp, del_s,       &
                         &  ref_cor, h_tan, tan_pt_f, do_gl, h_path, t_path,  &
                         &  ds_dh, dh_dz_gw, ds_dz_gw, dt_scr_dt, tau,        &
                         &  inc_rad_path, i_start, tan_pt, i_stop,            &
                         &  deriv_flags, pfa_update, drad_dt )

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP
    use SCRT_dN_m, only: dSCRT_dT, dSCRT_dX
    use Sparse_m, only: Sparse_t

! Inputs

    integer, intent(in) :: gl_inds(:)       ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: h_path_c(:)     ! path heights + req on main grid km.
    type(sparse_t), intent(in) :: dh_dt_path ! derivative of path height wrt
                                            ! temperature(km/K) on composite
                                            ! coarse & fine path.
    real(rp), intent(in), target :: alpha_path(:) ! path absorption(km^-1)
                                            ! on composite coarse & fine path
    real(rp), intent(in) :: dAlpha_dT_path(:) ! path dAlpha/dT on composite
                                            ! coarse & fine grid
    class(sparse_t), intent(in) :: Eta_ZxP  ! Interpolating coefficients from
                                            ! state vector to combined coarse &
                                            ! fine path for temperature only
    real(rp), intent(in) :: del_s(:)        ! unrefracted path length.
    real(rp), intent(in) :: ref_cor(:)      ! refracted to unrefracted path
                                            ! length ratios.
    real(rp), intent(in) :: h_tan           ! tangent height + req (km).
    integer, intent(in) :: tan_pt_f         ! Tangent point in the composite path
    logical, intent(in) :: do_gl(:)         ! Indicates where on the coarse path
                                            ! to do gl integrations.
    real(rp), intent(in) :: h_path(:)       ! path heights + req (km) on
                                            ! composite coarse & fine path.
    real(rp), intent(in) :: t_path(:)       ! path temperature(K) on 
                                            ! composite coarse & fine path.
    real(rp), intent(in) :: ds_dh(:)        ! path length wrt height derivative
                                            ! on complete grid.  Only the
                                            ! gl_inds part is used.
    real(rp), intent(in) :: dh_dz_gw(:)     ! path height wrt zeta derivative * gw
                                            ! on complete grid.  Only the
                                            ! gl_inds part is used.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative * gw
                                            ! on complete grid.  Only the
                                            ! gl_inds part is used.
    real(rp), intent(in) :: dt_scr_dt(:,:)  ! d t_script / d T * d T / d eta.
    real(rp), intent(in) :: tau(:)          ! transmission function.
    real(rp), intent(in) :: inc_rad_path(:) ! incremental radiance along the
                                            ! path.  t_script * tau.
    integer, intent(in) :: i_start          ! path start index + 1
    integer, intent(in) :: Tan_pt           ! Tangent point index in Del_Zeta
    integer, intent(in) :: i_stop           ! path stop index
    logical, intent(in) :: deriv_flags(:)   ! Indicates which temperature
                                            ! derivatives to do
    logical, intent(in) :: PFA_Update       ! Use DSCRT_DX instead of DSCRT_DT.

! Output
    real(rp), intent(out) :: drad_dt(:)     ! derivative of radiances wrt
                                            ! temperature state vector
                                            ! element. (K)

! Internals

    integer :: A, B, GA
    real(rp), target :: dh_dt_path_col(dh_dt_path%nRows)
    real(rp), pointer :: Alpha_Path_c(:), dh_dt_path_c(:)
    real(rp) :: dh_dt_tan
    integer :: i, i_begin, n_inds, n_path, no_to_gl, p_i, sv_i
    integer, target, dimension(1:size(inc_rad_path)) :: All_inds_B
    integer, target, dimension(1:size(inc_rad_path)) :: Inds_B, more_inds_B
    integer, pointer :: All_inds(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, pointer :: Inds(:)      ! inds => part_of_inds_B;  Indices
                                     ! on coarse path where do_calc.
    integer, pointer :: More_inds(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.
    integer :: NPF                   ! Number of points in fine path

    real(rp) :: d_Delta_dt(size(del_s,1)) ! path x sve.
      ! derivative of delta (incremental opacity) wrt temperature. (K)
    real(rp) :: Eta_ZxP_col(eta_zxp%nRows)
    real(rp) :: fa, fb
    real(rp) :: S_DEl_S                  ! Running sum of Del_S
    real(rp) :: Singularity(1:size(del_zeta)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.

    logical, pointer :: Do_Calc_c(:)
    logical, target :: Do_Calc_f(eta_zxp%nRows)
    logical, target :: Do_Calc_Hyd(dh_dt_path%nRows) ! On the whole path
    logical, pointer :: Do_Calc_Hyd_c(:)             ! On the coarse path
    logical :: Do_Calc_t(size(GL_Inds))
    logical :: Do_calc(1:size(del_zeta)) ! do_calc_c .or. ( do_gl .and. any
                                         ! of the corresponding do_calc_f ).
    logical :: NeedFA                    ! Need F(A) for hydrostatic

! Begin code

    n_path = size(del_zeta)

    alpha_path_c => alpha_path ( 1 :: ngp1 )
    dh_dt_path_col = 0
    dh_dt_path_c => dh_dt_path_col ( 1 :: ngp1 )
    do_calc_c => do_calc_f ( 1 :: ngp1 )
    do_calc_f = .false.
    do_calc_hyd_c => do_calc_hyd ( 1 :: ngp1 )
    do_calc_hyd_c = .false.
    eta_zxp_col = 0
    npf = size(alpha_path)

! Compute the opacity derivative singularity value

    do sv_i = 1, size(eta_zxp%cols)
      d_delta_dt = 0.0_rp
      drad_dt(sv_i) = 0.0
      if ( .not. deriv_flags(sv_i)) cycle  ! No derivatives for this column
      i_begin = i_start
      call dh_dt_path%get_col ( sv_i, dh_dt_path_col, do_calc_hyd )
      dh_dt_tan = dh_dt_path_col(tan_pt_f)
      if ( eta_zxp%cols(sv_i) /= 0 ) then  ! Column isn't empty
        call eta_zxp%get_col_vec_and_flags ( sv_i, eta_zxp_col, do_calc_f, &
                                           & last=npf )
        do_calc_t = do_calc_f(gl_inds)

! Do the absorption part
! Combine non zeros flags for both the main and gl parts

        call get_do_calc ( do_calc_c, do_calc_t, do_gl, do_calc, n_inds, inds_B )
        if ( n_inds > 0 ) then ! Column isn't empty up to NPF

          inds => inds_B(1:n_inds)

          no_to_gl = count(do_gl(inds))

          all_inds => all_inds_B(1:no_to_gl)
          more_inds => more_inds_B(1:no_to_gl)

          ! see if anything needs to be gl-d
          if ( no_to_gl > 0 ) &
            & call get_inds ( do_gl, do_calc, more_inds, all_inds )

          ! No ref_cor yet, no Lin_Log for temperature
          call get_d_delta_df ( inds, gl_inds, all_inds, more_inds,   &
            & eta_zxp_col, dAlpha_dT_path, del_s, del_zeta, ds_dz_gw, &
            & singularity, d_delta_dt )

          i_begin = max(inds(1)-1, i_start)

        end if ! n_inds > 0

      end if ! eta_zxp%cols(sv_i) /= 0

! Now do the hydrostatic part
! Combine boundaries flags

!????? Can the parts of the path to be calculated, and the Do_Calc flags,
!????? be handled as in Get_D_Deltau_Pol_dT?

      do_calc = do_calc_hyd_c
      do_calc(2:tan_pt) =          do_calc(tan_pt)   .or. do_calc(2:tan_pt)          .or. do_calc(1:tan_pt-1)
      do_calc(tan_pt+1:n_path-1) = do_calc(tan_pt+1) .or. do_calc(tan_pt+1:n_path-1) .or. do_calc(tan_pt+2:n_path)

! This is a layer calculation.  Before the tangent point, boundary J refers
! to the layer from J-1 to J.  After the tangent point, boundary J refers
! to the layer from J to J+1.  Therefore, we must require:

      do_calc((/1,n_path/)) = .false.

! Find where the non zeros are along the path

      n_inds = count(do_calc)
      needFA = .true.
      fa = 0.0_rp ! in case n_path <= 4
      if ( n_inds > 0 ) then
        inds => inds_B(1:n_inds)
        i = 1
        inds = 0
        s_del_s = sum(del_s(2:tan_pt))
        do p_i = 2 , tan_pt - 1
          if ( do_calc(p_i) ) then
            if ( needFA ) then ! only once in this loop
              fa = (h_path_c(p_i-1) * dh_dt_path_c(p_i-1) - &
                 &  h_tan * dh_dt_tan) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s - del_s(p_i)
            fb = (h_path_c(p_i) * dh_dt_path_c(p_i) - &
                & h_tan * dh_dt_tan) / s_del_s
            inds(i) = p_i
            d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i) * (fa - fb)
            fa = fb
            i = i + 1
          else
            s_del_s = s_del_s - del_s(p_i)
          end if
        end do ! p_i

! Special processing at tangent.  fb is zero

        if ( do_calc(tan_pt) ) then
          d_delta_dt(tan_pt) = d_delta_dt(tan_pt) + alpha_path_c(tan_pt) * fa
          inds(i) = tan_pt
          i = i + 1
        end if

        needFA = .not. do_calc(tan_pt+1)
        s_del_s = del_s(tan_pt+1)
        if ( do_calc(tan_pt+1) ) then
          fa = (h_path_c(tan_pt+2) * dh_dt_path_c(tan_pt+2) - &
              & h_tan * dh_dt_tan) / s_del_s
          d_delta_dt(tan_pt+1) = d_delta_dt(tan_pt+1) + &
              &                  alpha_path_c(tan_pt+1) * fa
          inds(i) = tan_pt + 1
          i = i + 1
        end if

        ! Several subscripts in this loop are offset by 1 from the nearly-
        ! identical loop above, and we use (fb-fa) here instead of (fa-fb).
        do p_i = tan_pt + 2, n_path - 1
          if ( do_calc(p_i) ) then
            if ( needFA ) then ! only once in this loop
              fa = (h_path_c(p_i) * dh_dt_path_c(p_i) - &
                 &  h_tan * dh_dt_tan) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s + del_s(p_i)
            fb = (h_path_c(p_i+1)*dh_dt_path_c(p_i+1) - &
               &  h_tan * dh_dt_tan) / s_del_s
            inds(i) = p_i
            d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i) * (fb - fa)
            fa = fb
            i = i + 1
          else
            s_del_s = s_del_s + del_s(p_i)
          end if
        end do ! p_i

        ! Do GL for hydrostatic for any panels that need it; the singularity
        ! correction is alpha_path_c.
        a = 1
        do p_i = 2, n_path-1           ! along the path
          if ( do_gl(p_i) ) then
            b = a + ng
            ga = gl_inds(a)
            ! Don't test do_calc: There may be GL corrections even if
            ! dh_dt_path_c (from whence came do_calc) is zero.
            d_delta_dt(p_i) = d_delta_dt(p_i) +                              &  
              & del_zeta(p_i) *                                              &  
              &  sum( ( alpha_path(ga:ga+ng-1) - alpha_path_c(p_i) ) *       &  
              &  (((2.0_rp*h_path(ga:ga+ng-1)**2 - 3.0_rp*h_tan**2) *        &  
              &    dh_dt_path_col(ga:ga+ng-1) +                              &  
              &    h_path(ga:ga+ng-1) * h_tan * dh_dt_tan) /                 &     
              &   (sqrt(h_path(ga:ga+ng-1)**2 - h_tan**2))**3                &
              &   + eta_zxp_col(ga:ga+ng-1) * ds_dh(ga:ga+ng-1) /            &
              &   t_path(ga:ga+ng-1)) * dh_dz_gw(ga:ga+ng-1) )
            a = b
          end if
        end do ! p_i

        i_begin = min(i_begin,inds(1))

      end if ! n_inds for hydrostatic > 0

! Correct for path length refraction

      d_delta_dt = ref_cor * d_delta_dt

! Accumulate the incremental opacity derivatives to get drad_dt

      if ( PFA_update ) then
        ! If we're doing a PFA update, we do not want to include
        ! dt_scr_dt again.

        call dscrt_dx ( tan_pt, d_delta_dt, inc_rad_path, i_begin, i_stop, &
                     &  drad_dt(sv_i) )

      else

        call dscrt_dt ( tan_pt, d_delta_dt, tau, inc_rad_path,&
                      & dt_scr_dt(:,sv_i),  i_begin, i_stop, drad_dt(sv_i) )

      end if

      call dh_dt_path%clear_col ( sv_i, dh_dt_path_col, do_calc_hyd )
      call eta_zxp%clear_col ( sv_i, eta_zxp_col, do_calc_f )

    end do ! sv_i

  end subroutine DRad_Tran_dT

!-----------------------------------------_______--  dRad_Tran_dX  -----
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )

  subroutine dRad_Tran_dX ( GL_Inds, Del_Zeta, Grids_f, Eta_fzp, Sps_Path,     &
                          & Sps_Map, dBeta_Path_c, dBeta_Path_f, Do_GL, Del_S, &
                          & Ref_Cor, ds_dz_gw, Inc_Rad_Path, Tan_Pt, I_Stop,   &
                          & dRad_dx )

    use Load_SPS_Data_m, ONLY: Grids_t
    use MLSKinds, only: RP
    use SCRT_dN_m, ONLY: dSCRT_dx
    use Sparse_m, only: Sparse_t

! Inputs

    integer, intent(in) :: GL_Inds(:)        ! Gauss-Legendre grid indicies
    real(rp), intent(in) :: Del_Zeta(:)      ! path -log(P) differences on the
                                             ! main grid.  This is for the whole
                                             ! coarse path, not just the part up
                                             ! to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    class(sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                             ! from state vector to combined
                                             ! coarse & fine path for each sps
    real(rp), intent(in) :: Sps_Path(:,:)    ! Path species function, path X species.
    integer, intent(in) :: Sps_Map(:)        ! second-dimension subscripts for sps_path.
    real(rp), intent(in) :: dBeta_Path_c(:,:) ! derivative of beta wrt dx
                                             ! on main grid.
    real(rp), intent(in) :: dBeta_Path_f(:,:) ! derivative of beta wrt dx
    logical, intent(in) :: Do_GL(:)          ! A logical indicating where to
                                             ! do gl integrations
    real(rp), intent(in) :: Del_S(:)         ! unrefracted path length.
    real(rp), intent(in) :: Ref_Cor(:)       ! refracted to unrefracted path
                                             ! length ratios.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
                                             ! gw on the entire grid.  Only the
                                             ! gl_inds part is used.
    real(rp), intent(in) :: Inc_Rad_Path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    integer, intent(in) :: Tan_pt            ! Tangent point index in inc_rad_path
    integer, intent(in) :: I_Stop            ! path stop index

! Outputs

    real(rp), intent(out) :: dRad_dx(:)      ! derivative of radiances wrt x
!                                              state vector element. (K)
! Internals

    integer :: n_inds, no_to_gl, sps_i, sps_m, sps_n, sv_i
    integer, target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer, target, dimension(1:size(inc_rad_path)) :: inds_B, more_inds_B
    integer, pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    real(rp) :: d_delta_dx(1:size(inc_rad_path))  ! derivative of delta
      !              wrt spectroscopy parameter. (K)
    logical :: do_calc(1:size(inc_rad_path)) ! Flags on coarse path where
                                     ! do_calc_fzp on the coarse path or (do_gl
                                     ! and any corresponding do_calc_fzp on the
                                     ! fine path).
    logical :: Do_Calc_FZP(n_eta_rows(eta_fzp,.true.))
    real(rp) :: Eta_FZP_Col(n_eta_rows(eta_fzp,.true.))
    integer, pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                     ! on coarse path where do_calc.
    integer, pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.
    integer :: NNZ_d_Delta_dx
    integer, target :: NZ_d_Delta_dx(1:size(d_Delta_dx))
    real(rp) :: singularity(1:size(inc_rad_path)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.
    integer :: Sparse_Col

! Begin code

    d_delta_dx = 0.0_rp
    do_calc_fzp = .false.
    eta_fzp_col = 0
    nnz_d_delta_dx = 0
    nz_d_delta_dx = 0
    sps_n = ubound(grids_f%l_z,1)

    do sps_i = 1 , sps_n
      sparse_col = 0
      sps_m = sps_map(sps_i)

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)
        sparse_col = sparse_col + 1

        ! We keep track of where we create nonzeros in d_delta_dx, and replace
        ! them by zeros on the next iteration.  This is done because the vast
        ! majority of d_delta_dx elements are zero, and setting all of
        ! d_Delta_dx was found to be a significant expense.

        ! Everything in d_delta_dx not indexed by nz_d_delta_dx is already zero

        d_delta_dx(nz_d_delta_dx(:nnz_d_delta_dx)) = 0
        nnz_d_delta_dx = 0

        ! Skip the masked derivatives, according to the l2cf inputs

        drad_dx(sv_i) = 0.0

        ! Get the interpolating coefficients along the path
        call eta_fzp(sps_i)%get_flags ( sparse_col, do_calc_fzp )
        call get_do_calc_indexed ( size(do_gl), tan_pt, do_calc_fzp, gl_inds, &
          & do_gl, do_calc, nnz_d_delta_dx, nz_d_delta_dx )
        call eta_fzp(sps_i)%clear_flags ( sparse_col, do_calc_fzp )
        if ( nnz_d_delta_dx == 0 ) cycle

        inds => nz_d_delta_dx(1:nnz_d_delta_dx)
        all_inds => all_inds_B(1:no_to_gl)
        more_inds => more_inds_B(1:no_to_gl)

        ! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

        ! Get d_Delta_dx for one state-vector element.

        if ( eta_fzp(sps_i)%cols(sparse_col) /= 0 ) then
          ! process only non-empty columns
          call eta_fzp(sps_i)%get_col ( sparse_col, eta_fzp_col )
          ! We're not really computing d_delta_df for lin-log mixing
          ! ratio.  It turns out that get_d_delta_df_linlog does the
          ! correct computation if we substitute dbeta_path for beta_path.
          call get_d_delta_df_linlog ( inds, gl_inds, all_inds, more_inds, &
            & eta_fzp_col, sps_path(:,sps_m), dbeta_path_c(:,sps_i),   &
            & dbeta_path_f(:,sps_i), del_s, del_zeta, ds_dz_gw, ref_cor,   &
            & grids_f%values(sv_i), singularity, d_delta_dx )
          call eta_fzp(sps_i)%clear_col ( sparse_col, eta_fzp_col )
        end if

        call dscrt_dx ( tan_pt, d_delta_dx, inc_rad_path, &
                      & 1, i_stop, drad_dx(sv_i))

      end do

    end do

  end subroutine dRad_Tran_dX

  ! ------------------------------------------------  Get_Do_Calc  -----
  subroutine Get_Do_Calc ( Do_Calc_c, Do_Calc_fzp, Do_GL, Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_c or Do_GL and any of the corresponding Do_Calc_fzp
  ! flags are set.

    use GLNP, ONLY: Ng

    logical, intent(in) :: Do_Calc_c(:)   ! On the coarse grid
    logical, intent(in) :: Do_Calc_fzp(:) ! On the GL grid
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out), optional :: N_Inds  ! count(do_calc)
    integer, intent(out), optional :: Inds(:) ! Indices where do_calc is true

    integer :: I, P_I

    i = 1
    do_calc = do_calc_c
    do p_i = 1 , size(do_gl)
      if ( do_gl(p_i) ) then
        do_calc(p_i) = do_calc(p_i) .or. any(do_calc_fzp(i:i+ng-1))
        i = i + Ng
      end if
    end do

    if ( present(n_inds) ) then
      n_inds = 0
      do p_i = 1 , size(do_gl)
        if ( do_calc(p_i) ) then
          n_inds = n_inds + 1
          inds(n_inds) = p_i
        end if
      end do
    end if
  end subroutine Get_Do_Calc

  ! .............................................  Get_d_Delta_df  .....
  subroutine Get_d_Delta_df ( Inds, GL_Inds, All_inds, More_inds, eta_fzp, &
    & dAlpha_df_path, Del_s, Del_Zeta, ds_dz_gw, &
    & Singularity, d_delta_df, Ref_cor, lin_log, grids_v )

    ! Get d_delta_df or d_delta_dT.  For species for which beta does not
    ! depend upon mixing ratio this gets d_delta_df if dAlpha_df_path_* is
    ! beta_path_*.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:)       ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)    ! Indices of GL points within combined
                                         ! coarse & fine path that are GL points
                                         ! for panels needing GL -- subset of
                                         ! f_inds (q.v. in  FullForwardModel).
    integer, intent(in) :: All_inds(:)   ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer, intent(in) :: More_inds(:)  ! Indices on the coarse path where
                                         ! GL corrections get applied.
    real(rp), intent(in) :: Eta_fzp(:)   ! Interpolation coefficients from state
                                         ! vector coordinates for one species to
                                         ! points on the combined coarse & fine
                                         ! path.
    real(rp), intent(in), target :: dAlpha_df_path(:) ! dAlpha_df on GL grid within
                                         ! subset of coarse path that needs GL.
    real(rp), intent(in) :: Del_s(:)     ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:)  ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:)  ! ds/dh * dh/dz * GL weights
    real(rp), intent(out) :: Singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_Delta_df(:) ! Derivative of delta.
                               ! intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: Ref_cor(:) ! refracted to unrefracted
                                         !  path length ratios.
    logical, intent(in), optional :: lin_log  ! logarithmic interpolation was used
    real(rp), intent(in), optional :: Grids_v ! Grids_f%values(sv_i)

    integer :: AA, GA, I, II
    real(rp), pointer :: dAlpha_df_path_c(:)

    ! Get a pointer to dAlpha_df on the coarse path only, to compute the
    ! initial rectangular estimate of d_delta_df
    dAlpha_df_path_c => dAlpha_df_path( 1 :: ngp1 )

    do i = 1, size(inds)
      ii = inds(i)
      singularity(ii) = dAlpha_df_path_c(ii) * eta_fzp(ii*ngp1-ng)
      d_delta_df(ii) = singularity(ii) * del_s(ii)
    end do ! i

    !{ Apply Gauss-Legendre quadrature to compute $\int_{\zeta_i}^{\zeta_{i-1}}
    !  \frac{\partial \alpha(s)}{\partial f^k_{lm}} \frac{\text{d}
    !  s}{\text{d} h} \frac{\text{d}h}{\text{d}\zeta} \, \text{d} s$ to the
    !  panels indicated by {\tt more\_inds}.  Here, $\frac{\partial
    !  \alpha(s)}{\partial f^k_{lm}} = \beta^k(s) \eta^k_{lm}(s)$.  We
    !  remove the singularity introduced at the tangent point by
    !  $\frac{\text{d} s}{\text{d} h}$ by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta = G(\zeta_i)
    !  \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !  \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !  \text{d}\zeta$.  The first integral is easy -- it's just $G(\zeta_i)
    !  (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}. In the
    !  second integral, $G(\zeta)$ is {\tt dAlpha_df_path\_f * eta\_zxp\_f
    !  * sps\_path} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}. The weights  are
    !  {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_df(ii) = d_delta_df(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp(ga:ga+ng-1) * dAlpha_df_path(ga:ga+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    if ( present(ref_cor) ) d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds)

    if ( present(lin_log) ) then
      ! Logarithmic interpolation correction
      if ( lin_log ) d_delta_df(inds) = d_delta_df(inds) * exp(-grids_v)
    end if

  end subroutine Get_d_Delta_df

  ! .............................................  Get_d_delta_df_old  .....
  subroutine Get_d_Delta_df_old ( Inds, GL_Inds, All_inds, More_inds, eta_fzp, &
    & dAlpha_df_path_c, dAlpha_df_path_f, Del_s, Del_Zeta, ds_dz_gw, &
    & Singularity, d_delta_df, Ref_cor, lin_log, grids_v )

    ! Get d_delta_df or d_delta_dT.  For species for which beta does not
    ! depend upon mixing ratio this gets d_delta_df if dAlpha_dx_path_* is
    ! beta_path_*.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:)       ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)    ! Indices of GL points within combined
                                         ! coarse & fine path that are GL points
                                         ! for panels needing GL -- subset of
                                         ! f_inds (q.v. in  FullForwardModel).
    integer, intent(in) :: All_inds(:)   ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer, intent(in) :: More_inds(:)  ! Indices on the coarse path where
                                         ! GL corrections get applied.
    real(rp), intent(in) :: Eta_fzp(*)   ! Interpolation coefficients from state
                                         ! vector coordinates for one species to
                                         ! points on the combined coarse & fine
                                         ! path.
    real(rp), intent(in) :: dAlpha_df_path_c(*) ! dAlpha_df on coarse grid.
    real(rp), intent(in) :: dAlpha_df_path_f(*) ! dAlpha_df on GL grid within
                                         ! subset of coarse path that needs GL.
    real(rp), intent(in) :: Del_s(:)     ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:)  ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:)  ! ds/dh * dh/dz * GL weights
    real(rp), intent(out) :: Singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_Delta_df(:) ! Derivative of delta.
                               ! intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: Ref_cor(:) ! refracted to unrefracted
                                         !  path length ratios.
    logical, intent(in), optional :: lin_log  ! logarithmic interpolation was used
    real(rp), intent(in), optional :: Grids_v ! Grids_f%values(sv_i)

    integer :: AA, GA, I, II

    do i = 1, size(inds)
      ii = inds(i)
      singularity(ii) = dAlpha_df_path_c(ii) * eta_fzp(ii*ngp1-ng)
      d_delta_df(ii) = singularity(ii) * del_s(ii)
    end do ! i

    !{ Apply Gauss-Legendre quadrature to compute $\int_{\zeta_i}^{\zeta_{i-1}}
    !  \frac{\partial \alpha(s)}{\partial f^k_{lm}} \frac{\text{d}
    !  s}{\text{d} h} \frac{\text{d}h}{\text{d}\zeta} \, \text{d} s$ to the
    !  panels indicated by {\tt more\_inds}.  Here, $\frac{\partial
    !  \alpha(s)}{\partial f^k_{lm}} = \beta^k(s) \eta^k_{lm}(s)$.  We
    !  remove the singularity introduced at the tangent point by
    !  $\frac{\text{d} s}{\text{d} h}$ by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta = G(\zeta_i)
    !  \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !  \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !  \text{d}\zeta$.  The first integral is easy -- it's just $G(\zeta_i)
    !  (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}. In the
    !  second integral, $G(\zeta)$ is {\tt dAlpha_df_path\_f * eta\_zxp\_f
    !  * sps\_path} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}. The weights  are
    !  {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_df(ii) = d_delta_df(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp(ga:ga+ng-1) * dAlpha_df_path_f(aa:aa+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    if ( present(ref_cor) ) d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds)

    if ( present(lin_log) ) then
      ! Logarithmic interpolation correction
      if ( lin_log ) d_delta_df(inds) = d_delta_df(inds) * exp(-grids_v)
    end if

  end subroutine Get_d_delta_df_old

  ! ...........................................  Get_d_delta_df_f  .....
  subroutine Get_d_delta_df_f ( Inds, GL_Inds, All_inds, More_inds, eta_fzp, &
    & Sps_path, Beta_path_c, Beta_path_f, dBeta_df_c, dBeta_df_f, Del_s,     &
    & Del_Zeta, ds_dz_gw, Ref_cor, Singularity, d_delta_df )

    ! Get d_delta_df for the case of species for which beta
    ! depends upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:)      ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Indices of GL points within combined
                                        ! coarse & fine path that are GL points
                                        ! for panels needing GL -- subset of
                                        ! f_inds (q.v. in  FullForwardModel).
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: eta_fzp(*)  ! representation basis function.
    real(rp), intent(in) :: Sps_path(:) ! Path mixing ratios
    real(rp), intent(in) :: Beta_path_c(*)  ! cross section on coarse grid.
    real(rp), intent(in) :: Beta_path_f(*)  ! cross section on GL grid.
    real(rp), intent(in) :: dBeta_df_c(*)   ! beta depends on mixing ratio
    real(rp), intent(in) :: dBeta_df_f(*)   ! beta depends on mixing ratio
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:)     ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: ref_cor(:)      ! refracted to unrefracted path
                                            !  length ratios.
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_delta_df(:) ! Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)
      ii = inds(i)
      iii = ii*ngp1 - ng
      singularity(ii) = eta_fzp(iii) &
                & * ( beta_path_c(ii) + &
                &     sps_path(iii) * dBeta_df_c(ii) )
      d_delta_df(ii) = singularity(ii) * del_s(ii)
    end do ! i
    !{ Apply Gauss-Legendre quadrature to compute $\int_{\zeta_i}^{\zeta_{i-1}}
    !  \frac{\partial \alpha(s)}{\partial f^k_{lm}} \frac{\text{d}
    !  s}{\text{d} h} \frac{\text{d}h}{\text{d}\zeta} \, \text{d} s$ to the
    !  panels indicated by {\tt more\_inds}.  Here, $\frac{\partial
    !  \alpha(s)}{\partial f^k_{lm}} = \left( \beta^k(s) + f^k(s)
    !  \frac{\partial \beta(s)}{\partial f^k(s)} \right) \eta^k_{lm}(s)$.  We
    !  remove the singularity introduced at the tangent point by
    !  $\frac{\text{d} s}{\text{d} h}$ by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta = G(\zeta_i)
    !  \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !  \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !  \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !  \text{d}\zeta$.  The first integral is easy -- it's just $G(\zeta_i)
    !  (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}. In the second
    !  integral, $G(\zeta)$ is {\tt eta\_zxp\_f * ( beta\_path\_f + sps\_path
    !  * dBeta\_df )} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}.  The weights are
    !  {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_df(ii) = d_delta_df(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp(ga:ga+ng-1) &
             &  * ( beta_path_f(aa:aa+ng-1) &
             &      + sps_path(aa:aa+ng-1) * dBeta_df_f(aa:aa+ng-1) ) &
             &  - singularity(ii)) &
             & * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds)

  end subroutine Get_d_delta_df_f

  ! ......................................  Get_d_delta_df_linlog  .....
  subroutine Get_d_delta_df_linlog ( Inds, GL_Inds, All_inds, More_inds, &
    & eta_fzp, Sps_path, Beta_path_c, Beta_path_f, Del_s, Del_Zeta,      &
    & ds_dz_gw, Ref_cor, Grids_v, Singularity, d_delta_df )

    ! Get d_delta_df for the case of lin_log species for which beta
    ! does not depend upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:) ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: eta_fzp(*)  ! representation basis function.
    real(rp), intent(in) :: Sps_path(:) ! exp(Path mixing ratios)
    real(rp), intent(in) :: Beta_path_c(*)  ! cross section on coarse grid.
    real(rp), intent(in) :: Beta_path_f(*)  ! cross section on GL grid.
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:) ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: ref_cor(:)  ! refracted to unrefracted path
                                        !  length ratios.
    real(rp), intent(in) :: Grids_v     ! Grids_f%values(sv_i)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_delta_df(:) ! Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)
      ii = inds(i)
      iii = ii*ngp1 - ng
      singularity(ii) = eta_fzp(iii) * sps_path(iii) * &
                & beta_path_c(ii)
      d_delta_df(ii) = singularity(ii) * del_s(ii)
    end do ! i

    !{ Apply Gauss-Legendre quadrature to compute $\int_{\zeta_i}^{\zeta_{i-1}}
    !  \frac{\partial \alpha(s)}{\partial f^k_{lm}} \frac{\text{d} s}{\text{d}
    !  h} \frac{\text{d}h}{\text{d}\zeta} \, \text{d} s$ to the panels
    !  indicated by {\tt more\_inds}.  Here, $\frac{\partial
    !  \alpha(s)}{\partial f^k_{lm}} = \hat{f}^k(s) \beta^k(s)
    !  \frac{\eta^k_{lm}(s)}{f^k_{lm}}$, where $\hat{f}^k(s) = \exp\left(
    !  \sum_{lm} \eta^k_{lm}(s) \ln f^k_{lm} \right)$.  We remove the
    !  singularity introduced at the tangent point by $\frac{\text{d}
    !  s}{\text{d} h}$ by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
    !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !   \text{d}\zeta$.  The first integral is easy -- it's just
    !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}.
    !  In the second integral, $G(\zeta)$ is {\tt beta\_path\_f * eta\_zxp\_f *
    !  sps\_path} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}.  The weights
    !  are {\tt gw}.  {\tt sps\_path} is $\hat{f}^k(s)$, and {\tt grids_v}
    !  is $\ln f^k_{lm}$.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_df(ii) = d_delta_df(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp(ga:ga+ng-1) * sps_path(ga:ga+ng-1) &
             &  * beta_path_f(aa:aa+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds) * exp(-grids_v)

  end subroutine Get_d_delta_df_linlog

  ! ....................................  Get_d_delta_df_linlog_f  .....
  subroutine Get_d_delta_df_linlog_f ( Inds, GL_Inds, All_inds, More_inds, &
    & eta_fzp, Sps_path, Beta_path_c, Beta_path_f, dBeta_df_c, dBeta_df_f, &
    & Del_s, Del_Zeta, ds_dz_gw, Ref_cor, Grids_v,  Singularity, d_delta_df )

    ! Get d_delta_df for the case of lin_log species for which beta
    ! depends upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:) ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: eta_fzp(*)  ! representation basis function.
    real(rp), intent(in) :: Sps_path(:) ! exp(Path mixing ratios)
    real(rp), intent(in) :: Beta_path_c(*)  ! cross section on coarse grid.
    real(rp), intent(in) :: Beta_path_f(*)  ! cross section on GL grid.
    real(rp), intent(in) :: dBeta_df_c(*)   ! beta depends on mixing ratio
    real(rp), intent(in) :: dBeta_df_f(*)   ! beta depends on mixing ratio
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:) ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: ref_cor(:)  ! refracted to unrefracted path
                                        !  length ratios.
    real(rp), intent(in) :: Grids_v     ! Grids_f%values(sv_i)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_delta_df(:) ! Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)
      ii = inds(i)
      iii = ii*ngp1 - ng
      singularity(ii) = eta_fzp(iii) * sps_path(iii) * &
                & ( beta_path_c(ii) + sps_path(iii) * dBeta_df_c(ii) )
      d_delta_df(ii) = singularity(ii) * del_s(ii)
    end do ! i

    !{ Apply Gauss-Legendre quadrature to compute $\int_{\zeta_i}^{\zeta_{i-1}}
    !  \frac{\partial \alpha(s)}{\partial f^k_{lm}} \frac{\text{d} s}{\text{d}
    !  h} \frac{\text{d}h}{\text{d}\zeta} \, \text{d} s$ to the panels
    !  indicated by {\tt more\_inds}.  Here, $\frac{\partial
    !  \alpha(s)}{\partial f^k_{lm}} = \left( \beta^k(s) + \hat{f}^k(s)
    !  \frac{\partial \beta(s)}{\partial \hat{f}^k(s)} \right) \hat{f}^k(s)
    !  \frac{\eta^k_{lm}(s)}{f^k_{lm}}$, where $\hat{f}^k(s) = \exp\left(
    !  \sum_{lm} \eta^k_{lm}(s) \ln f^k_{lm} \right)$.  We remove the
    !  singularity introduced at the tangent point by $\frac{\text{d}
    !  s}{\text{d} h}$ by writing
    !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
    !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
    !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
    !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
    !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
    !   \text{d}\zeta$.  The first integral is easy -- it's just
    !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}.
    !  In the second integral, $G(\zeta)$ is {\tt eta\_zxp\_f * sps\_path
    !  ( beta\_path\_f + dBeta_df )} -- which have already been evaluated at
    !  the appropriate abscissae~-- and $G(\zeta_i)$ is {\tt singularity}. 
    !  The weights are {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_df(ii) = d_delta_df(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp(ga:ga+ng-1) * sps_path(ga:ga+ng-1) &
             &  * ( beta_path_f(aa:aa+ng-1) + &
             &      sps_path(ga:ga+ng-1) * dBeta_df_f(aa:aa+ng-1) ) &
             &  - singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds) * exp(-grids_v)

  end subroutine Get_d_delta_df_linlog_f

! =====     Private Procedures     =====================================

  ! ----------------------------------------  Get_Do_Calc_Indexed  -----
  subroutine Get_Do_Calc_Indexed ( N, Tan_Pt_C, Do_Calc_All, F_Inds, Do_GL, &
    & Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_All(1::ngp1), or Do_GL and any of the corresponding
  ! Do_Calc_All(f_inds) flags are set.
  ! Get_Do_Calc_Indexed determines that Get_d_Delta_df needs to integrate a
  ! panel if Do_Calc_All(1::ngp1) is true, which means the interpolating
  ! coefficient is nonzero, or if Do_GL is true and Do_Calc_All(f_inds) is
  ! true, which means that GL was needed, even if some interpolating
  ! coefficient is zero.

    use GLNP, ONLY: Ng, NGP1

    integer, intent(in) :: N              ! sizes on coarse grid
    integer, intent(in) :: Tan_Pt_C       ! Index of tangent point in coarse grid

#if defined NAG
!   Assumed-shape arguments are slower than assumed size
    logical, intent(in) :: Do_Calc_all(*) ! On the entire path
    integer, intent(in) :: F_Inds(*)      ! Indices in Do_Calc_All for fine grid
                                          ! GL points on panels needing GL.
    logical, intent(in) :: Do_GL(*)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(*)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(*)       ! Indices where do_calc is true
#elif defined IFC
!   Contiguous assumed-shape arguments are no slower than assumed size
    logical, contiguous, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, contiguous, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
                                                      ! GL points on panels needing GL.
    logical, contiguous, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, contiguous, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer,             intent(out) :: N_Inds        ! count(do_calc)
    integer, contiguous, intent(out) :: Inds(:)       ! Indices where do_calc is true
#elif defined LF95
!   Assumed-shape arguments are faster than assumed size because of copying
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
                                          ! GL points on panels needing GL.
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(:)       ! Indices where do_calc is true
#else
!   We don't know whether assumed-shape or assumed size arguments are faster
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
                                          ! GL points on panels needing GL.
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(:)       ! Indices where do_calc is true
#endif

    integer :: I, P_I
!     integer :: J
    integer :: K
!     logical :: T_calc(size(f_inds)/ng)

!     do_calc = do_calc_all(1::ngp1)
    i = 1 - Ng
    n_inds = 0
!     do_calc(1:n) = do_calc_all(ngp1-ng:n*ngp1:ngp1)
!     do_calc(tan_pt_c) = .false. ! Don't do GL between tangent points
    do p_i = 1, n
      do_calc(p_i) = do_calc_all(ngp1*p_i-ng)
      if ( do_gl(p_i) ) then
        i = i + Ng
        k = f_inds(i)
!       This is fastest with NAG -g, NAG -O and lf95 -O but it needs to be
!       changed by hand if NG is changed
        do_calc(p_i) = do_calc(p_i) .or. &
          & do_calc_all(k) .or. do_calc_all(k+1) .or. do_calc_all(k+2)
! !       This isn't quite as fast
!         do j = 0, ng-1
!           do_calc(p_i) = do_calc(p_i) .or. do_calc_all(k+j)
!         end do
! !       These are even slower
!         k = f_inds(i)
!         do_calc(p_i) = do_calc(p_i) .or. any(do_calc_all(k:k+ng-1))
!         if ( any(do_calc_all(f_inds(i:i+ng-1))) ) do_calc(p_i)=.true.
      end if
      if ( do_calc(p_i) ) then
        n_inds = n_inds + 1
        inds(n_inds) = p_i
      end if
    end do
! !   This is much slower
!     do p_i = 1, size(t_calc)
!       k = f_inds((p_i-1)*ng+1)
!       t_calc(p_i) = do_calc_all(k)
!       do j = 1, ng-1
!         t_calc(p_i) = t_calc(p_i) .or. do_calc_all(k+j)
!       end do
!     end do
!     do_calc = do_calc .or. unpack(t_calc,do_gl,.false.)

  end subroutine Get_Do_Calc_Indexed

  ! ---------------------------------------------------  Get_Inds  -----
  subroutine Get_Inds ( Do_GL, Do_Calc, More_Inds, All_Inds )

    ! More_Inds are the places in the coarse path where both Do_Calc and Do_GL.
    ! All_Inds are the corresponding places in the GL-extracted fine path.

    ! Get_Inds is similar to Get_Do_Calc_Indexed, but computes indices 
    ! for the coarse path, or the coarse path points in the composite path.

    use GLNP, ONLY: Ng

    implicit NONE

  ! Inputs
    logical, intent(in) :: Do_GL(:)      ! Coarse path flag indicating where to
      !                                    do gl integrations.
    logical, intent(in) :: Do_Calc(:)    ! Coarse path flag indicating where
                                         ! any interpolating coeffient is nonzero
  ! Outputs
    integer, intent(out) :: More_Inds(:) ! Indices on coarse path where GL is
                                         ! to be used.
    integer, intent(out) :: All_Inds(:)  ! First index in GL_Inds of a fine-path
                                         ! point where a GL ordinate is to be
                                         ! computed.  The remaining NG-1 indices
                                         ! are consecutive.

    integer :: I, J, L, P_I
!   integer :: K

    i = 1
    j = 1
    l = 1
    do p_i = 2 , size(do_gl)-1
      if ( do_gl(p_i) ) then
        if ( do_calc(p_i) ) then
          more_inds(i) = p_i
! !         all_inds(j:j+ng-1) = (/ ( (p_i-2)*ng+k, k = 1, ng ) /)
!           all_inds(j:j+ng-1) = (/ ( l + k, k = 0, ng-1 ) /)
!           forall( k = 0: ng - 1 ) all_inds(j+k) = l + k
          ! OK, here be dragons.
          ! Van originally had this line:
          !         all_inds(j:j+ng-1) = (/ ( (p_i-2)*ng+k, k = 1, ng ) /)
          ! Then he replaced it with this one
          !         all_inds(j:j+ng-1) = (/ ( l + k, k = 0, ng-1 ) /)
          ! It turns out, through a bizarre and worrying series of discoveries about Lahey
          ! that that was really slow under some repeatable, but totally unjustifiable circumstances.
          ! hence the regular do loop below - NJL / WVS
!           do k = 0, ng - 1
!             all_inds(j+k) = l+k
!           end do
          ! We actually only need the first of the GL inds, since the rest
          ! are consecutive.
          all_inds(j) = l
          i = i + 1
          j = j + 1
        end if
        l = l + Ng
      end if
    end do

  end subroutine Get_Inds

  pure integer function N_Eta_Rows ( Eta_FZP, Need ) result ( N )
    use Sparse_m, only: Sparse_t
    class(sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                              ! from state vector to combined
                                              ! coarse & fine path for each sps.
    logical, intent(in) :: Need               ! Need to do the computation
    n = 0
    if ( need ) then ! need space for stuff; Some stuff is needed only for TScat
      n = maxval(eta_fzp%nRows)
    end if
  end function N_Eta_Rows

!----------------------------------------------------------------------
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, not_used_here ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module RAD_TRAN_M

! $Log$
! Revision 2.40  2018/11/19 21:55:37  vsnyder
! Make Get_d_Delta_df private
!
! Revision 2.39  2018/09/12 22:50:44  vsnyder
! Delete dRad_Tran_dX.  Change the name of dRad_Tran_dX_Sparse to dRad_Tran_dX.
!
! Revision 2.38  2018/09/12 22:03:38  vsnyder
! Add dRad_Tran_dX_Sparse.  Inline code from do_delta_m.
!
! Revision 2.37  2018/05/24 03:24:36  vsnyder
! Use sparse representation for dh_dt_path
!
! Revision 2.36  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.35  2017/08/09 20:47:34  vsnyder
! Add Tan_Pt_C argument, but it isn't actually used yet
!
! Revision 2.34  2017/03/31 00:46:30  vsnyder
! Remove Get_d_delta_dx because it's subsumed by Get_d_delta_df
!
! Revision 2.33  2017/03/17 20:19:23  vsnyder
! Cannonball polishing -- changed two variable names
!
! Revision 2.32  2014/09/05 21:26:38  vsnyder
! Use the CONTIGUOUS attribute with ifort
!
! Revision 2.31  2013/07/13 00:03:21  vsnyder
! Remove LD argument from get_all_d2_delta_df2
!
! Revision 2.30  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.29  2011/11/09 00:17:44  vsnyder
! Remove non-standard TAB characters
!
! Revision 2.28  2011/08/20 00:44:02  vsnyder
! Get rid of DOS line ends
!
! Revision 2.27  2011/07/29 01:59:24  vsnyder
! Cannonball polishing
!
! Revision 2.26  2011/07/08 21:25:58  yanovsky
! Use d2Alpha_df2
!
! Revision 2.25  2011/06/02 22:43:12  yanovsky
! In d2Rad_tran_df2 subroutine, add computations of analytical 
! mixing ratio Hessians in logarithmic basis
!
! Revision 2.24  2011/03/25 20:46:59  vsnyder
! Delete declarations of unused objects
!
! Revision 2.23  2011/03/24 00:17:34  vsnyder
! Add TScat derivatives to dRad_tran_df
!
! Revision 2.22  2011/03/23 23:45:32  vsnyder
! This log entry is bogus.  Check in again to get the right one.
! FOV_Convolve_m.f90
!
! Revision 2.21  2011/03/11 03:09:08  vsnyder
! Use Get_dAlpha_df
!
! Revision 2.20  2011/03/04 03:41:25  vsnyder
! Remove declaration for unused variable
!
! Revision 2.19  2011/02/12 03:57:40  vsnyder
! Add mixing-ratio dependence for H2O derivatives
!
! Revision 2.18  2011/02/05 01:18:06  vsnyder
! Correct bugs where dBeta_df is used
!
! Revision 2.17  2011/01/28 19:17:11  vsnyder
! Lots of stuff for TScat derivatives
!
! Revision 2.16  2010/12/07 01:20:57  vsnyder
! dRad_tran_dx needs to call dscrt_dx
!
! Revision 2.15  2010/11/05 20:28:13  vsnyder
! Delete unused declarations
!
! Revision 2.14  2010/08/27 23:17:11  vsnyder
! Remove '(ip)' from all integer declarations
!
! Revision 2.13  2010/08/27 05:51:03  yanovsky
! Changed types of indices_c and gl_inds dummy arguments from integer(ip) to integer.
!
! Revision 2.12  2010/08/19 02:04:03  vsnyder
! Change some variable names, restructure dRad_tran_dT
!
! Revision 2.11  2010/06/12 01:27:59  vsnyder
! Use get_d_delta_df_linlog in drad_tran_dx because it does the right
! calculation if it's given dbeta_df instead of beta.
!
! Revision 2.10  2010/06/11 02:20:46  vsnyder
! Integrated Igor's latest changes with mine
!
! Revision 2.9  2010/06/09 23:13:01  yanovsky
! Removed code implementing d_delta_df computations from drad_tran_df. 
! Drad_tran_df now calls get_d_delta_df
!
! Revision 2.8  2010/06/09 22:04:14  yanovsky
! Added d2Rad_tran_df2 and get_d_delta_df subroutines
!
! Revision 2.7  2010/05/14 02:41:08  vsnyder
! Changed intent for [n]nz_d_delta_df from out to inout.  Replaced
! indices_c(ii) by iii at a place where iii == indices(ii)
!
! Revision 2.6  2010/01/23 01:21:24  vsnyder
! Handle derivatives for betas that depend upon mixing ratio
!
! Revision 2.5  2009/06/16 17:37:09  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.4  2009/06/13 01:11:55  vsnyder
! Specify start and end of path, simplify some index calculations
!
! Revision 2.3  2007/07/11 22:26:45  vsnyder
! Dumps
!
! Revision 2.2  2007/06/26 00:40:01  vsnyder
! Minor improvement in Get_Do_Calc_Indexed
!
! Revision 2.1  2007/06/08 22:05:13  vsnyder
! Replacing rad_tran_m.f90 by rad_tran_m.F90
!
!
! Deleted rad_tran_m.f90 to replace it with rad_tran_m.F90
!
! Revision 2.52  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.51  2006/07/20 01:09:27  vsnyder
! Make sure drad_df gets a value
!
! Revision 2.50  2006/04/11 18:31:58  vsnyder
! Cannonball polishing
!
! Revision 2.49  2006/02/08 01:02:01  vsnyder
! More stuff for spectroscopy derivatives
!
! Revision 2.48  2005/11/21 22:57:27  vsnyder
! PFA derivatives stuff
!
! Revision 2.47  2005/11/01 23:02:21  vsnyder
! PFA Derivatives
!
! Revision 2.46  2005/09/17 00:49:54  vsnyder
! Revise arrays for spectroscopic derivatives, plus some cannonball polishing
!
! Revision 2.45  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.44  2005/04/26 15:35:54  livesey
! Minor changes necessitated by wierdo problems with LF95.  Probably won't
! fix the problem but might as well keep them.
!
! Revision 2.43  2005/03/28 20:24:37  vsnyder
! Taus past the black-out are zero, not one!
!
! Revision 2.42  2005/03/03 02:07:42  vsnyder
! Remove USEs for unreferenced symbols
!
! Revision 2.41  2004/11/01 20:25:44  vsnyder
! Reorganization of representation for molecules and beta groups; PFA may be broken for now
!
! Revision 2.40  2004/10/06 21:18:24  vsnyder
! Add rad_tran_PFA
!
! Revision 2.39  2004/08/03 22:06:46  vsnyder
! Inching further toward PFA
!
! Revision 2.38  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.37  2004/03/20 01:15:39  jonathan
!  remove rad_tran_cld
!
! Revision 2.36  2004/03/08 22:58:03  vsnyder
! Remove D_Delta_DT, which was schlepped from drad_tran_dt to
! Get_D_Delta_Pol_DT, but is no longer needed in the latter place.
!
! Revision 2.35  2004/02/03 02:47:55  vsnyder
! Progress (hopefully) on polarized temperature derivatives
!
! Revision 2.34  2004/01/23 01:16:05  vsnyder
! Repair mistakes in polarized radiance calculation:  CS_EXPMAT needs to be
! applied to incoptdepth_pol even if no GL is done, because the earlier
! calculation (in FullForwardModel) didn't have ref_cor.  We can't stop
! doing CS_EXPMAT when the sum of the eigenvalues gets large and negative,
! because the terms in the matrix are large if the difference of the
! eigenvalues is large.  If we had both eigenvalues, we'd be nearly done
! with CS_EXPMAT, so avoiding it would be as expensive as doing it.
!
! Revision 2.33  2003/12/08 17:52:47  jonathan
! update for 2d cldfwm
!
! Revision 2.32  2003/12/03 00:25:32  vsnyder
! Corrections to hydrostatic calculation
!
! Revision 2.31  2003/11/04 02:49:50  vsnyder
! Use GC_INDS calculated in FullForwardModel for more_inds
!
! Revision 2.30  2003/11/04 01:55:50  vsnyder
! Add 'FA = 0.0' in case n_path <= 4, cosmetic changes
!
! Revision 2.29  2003/11/01 03:04:02  vsnyder
! Use ds_dz_gw instead of ds_dh, dh_dz and gw; use del_zeta from FullForwardModel
!
! Revision 2.28  2003/10/30 20:36:41  vsnyder
! Get del_zeta from FullForwardModel
!
! Revision 2.27  2003/10/16 23:06:09  vsnyder
! Polish up some comments
!
! Revision 2.26  2003/10/15 02:04:08  vsnyder
! Simplifications possible after inlining path_opacity.  Cosmetic changes.
! Make Get_Del_Zeta_All public.  Don't bother checking do_calc(1) and
! do_calc(n_path) because we know it's false there.
!
! Revision 2.25  2003/10/09 21:04:38  vsnyder
! Fix typos generated while inlining path_opacity
!
! Revision 2.23  2003/09/25 20:06:03  vsnyder
! Insert TeXnicalities.  Insert many more comments too.  Inline path_opacity,
! which results in substantial savings in derivative calculations because it
! avoids constructing an array temp within a doubly-nested loop.  It doesn't
! make much differences for radiance calculations.
!
! Revision 2.22  2003/09/24 22:19:55  vsnyder
! Get rid of some array temps
!
! Revision 2.21  2003/09/09 22:34:45  vsnyder
! Don't look at status from cs_expmat if it isn't called
!
! Revision 2.20  2003/09/09 00:02:55  vsnyder
! Make deltau_pol inout, only compute it where needed
!
! Revision 2.19  2003/08/15 18:50:22  vsnyder
! Preparing the way for polarized vmr derivatives
!
! Revision 2.18  2003/06/27 22:05:48  vsnyder
! Check status from cs_expmat
!
! Revision 2.17  2003/06/18 17:24:05  bill
! added temperature derivative subsetting
!
! Revision 2.16  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.15  2003/05/20 00:04:28  vsnyder
! Collect common stuff into subroutines
!
! Revision 2.14  2003/05/15 03:29:00  vsnyder
! Moved some stuff up to FullForwardModel because Get_d_Deltau_pol_dT needs it
!
! Revision 2.13  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.12.2.2  2003/03/20 01:42:26  vsnyder
! Revise Grids_T structure
!
! Revision 2.12.2.1  2003/03/05 03:31:13  vsnyder
! Get rid of unused variables, use 'any' instead of 'count'
!
! Revision 2.12  2003/02/07 02:35:48  vsnyder
! OOPS, forgot to move one down
!
! Revision 2.11  2003/02/07 02:12:00  vsnyder
! Move some USE statements down
!
! Revision 2.10  2003/02/03 19:00:52  bill
! changed interface to rad tran to speed up program
!
! Revision 2.9  2003/01/08 00:15:42  vsnyder
! Moved path_contrib to its own module
!
! Revision 2.8  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.7  2002/10/02 20:09:48  vsnyder
! Use automatic arrays to move allocate/deallocate out of loops.  Numerous
! cosmetic changes.
!
! Revision 2.6  2002/07/05 07:52:52  zvi
! Some cosmetic changes
!
! Revision 2.5  2002/06/13 22:38:40  bill
! some variable name changes--wgr
!
! Revision 2.4  2002/06/04 10:28:04  zvi
! rename n_sps to: no_mol, more correctly
!
! Revision 2.3  2002/02/16 06:38:05  zvi
! Some cosmetic changes..
!
! Revision 2.2  2002/01/30 01:11:22  zvi
! Fix bug in user selectable coeff. code
!
! Revision 2.1  2002/01/27 08:37:51  zvi
! Adding Users selected coefficients for derivatives
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.11.2.3  2001/09/13 22:51:24  zvi
! Separating allocation stmts
!
! Revision 1.11.2.2  2001/09/11 00:00:46  zvi
! Adding dt_script code
!
! Revision 1.11.2.1  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
