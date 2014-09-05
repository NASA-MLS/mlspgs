! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module RAD_TRAN_M

  implicit NONE
  private
  public :: RAD_TRAN, RAD_TRAN_POL
  public :: DRAD_TRAN_DF, DRAD_TRAN_DT, DRAD_TRAN_DX
  public :: D2RAD_TRAN_DF2
  public :: Get_all_d_delta_df
  public :: Get_d_delta_df_f, Get_d_delta_df_linlog
  public :: Get_d_delta_df_linlog_f, Get_d_delta_dx
  public :: Get_d2_delta_df2_linlog
  public :: Get_Do_Calc
  private :: Get_Do_Calc_Indexed, Get_Inds

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!------------------------------------------------------  Rad_tran  -----
! This is the radiative transfer model, radiances only !

  subroutine Rad_tran ( tan_pt, gl_inds, more_inds, e_rflty, del_zeta, &
                      & alpha_path_c, ref_cor, incoptdepth,            &
                      & alpha_path_gl, ds_dz_gw, t_script,             &
                      & tau, inc_rad_path, rad, i_stop )

    use GLNP, only: NG
    use MLSKinds, only: RP
    use SCRT_DN_M, ONLY: SCRT

  ! inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    integer, intent(in) :: more_inds(:)      ! Places in the coarse path
  !                                            where GL is needed
    real(rp), intent(in) :: e_rflty          ! earth reflectivity value (0--1).
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: alpha_path_c(:)  ! absorption coefficient on coarse
  !                                            grid.
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
  !                                            length ratios.
    real(rp), intent(inout) :: incoptdepth(:) ! incremental path opacities
  !                            from one-sided layer calculation on output.
  !                            it is the full integrated layer opacity.
    real(rp), intent(in) :: alpha_path_gl(:) ! absorption coefficient on gl
  !                                            grid.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative * gw.
    real(rp), intent(in) :: t_script(:)      ! differential temperatures (K)
  !                                            on coarse grid.
  ! outputs

    real(rp), intent(out) :: tau(:)          ! transmission function.
    real(rp), intent(out) :: inc_rad_path(:) ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(out) :: rad             ! radiance (K)
    integer, intent(out) :: i_stop           ! path stop index

  ! Internals

    integer :: A, AA, I
  !                                            where GL is needed

  ! Begin code

  ! see if anything needs to be gl-d

    if ( size(gl_inds) > 0 ) then

      !{ Apply Gauss-Legendre quadrature to the panels indicated by
      !  {\tt more\_inds}.  We remove a singularity (which actually only
      !  occurs at the tangent point) by writing
      !  $\int_{\zeta_i}^{\zeta_{i-1}} G(\zeta) \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta =
      !  G(\zeta_i) \int_{\zeta_i}^{\zeta_{i-1}} \frac{\text{d}s}{\text{d}h}
      !   \frac{\text{d}h}{\text{d}\zeta} \text{d}\zeta +
      !  \int_{\zeta_i}^{\zeta_{i-1}} \left[ G(\zeta) - G(\zeta_i) \right]
      !   \frac{\text{d}s}{\text{d}h} \frac{\text{d}h}{\text{d}\zeta}
      !   \text{d}\zeta$.  The first integral is easy -- it's just
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt incoptdepth}.
      !  In the second integral, $G(\zeta)$ is {\tt alpha\_path\_gl} --
      !  which has already been evaluated at the appropriate abscissae -- and
      !  $G(\zeta_i)$ is {\tt alpha\_path\_c}.  The weights are {\tt gw}.

      a = 1
      do i = 1, size(more_inds)
        aa = gl_inds(a)
        incoptdepth(more_inds(i)) = incoptdepth(more_inds(i)) + &
          & del_zeta(more_inds(i)) * &
          & dot_product( (alpha_path_gl(a:a+ng-1) - alpha_path_c(more_inds(i))), &
               & ds_dz_gw(aa:aa+ng-1) )
        a = a + ng
      end do ! i

    end if

    incoptdepth = ref_cor * incoptdepth

    call scrt ( tan_pt, t_script, e_rflty, incoptdepth, tau, rad, inc_rad_path, &
      &         i_stop )

  end subroutine Rad_tran

!--------------------------------------------------  Rad_Tran_Pol  -----

  subroutine Rad_tran_Pol ( tan_pt, gl_inds, more_inds, e_rflty, del_zeta, &
                     &  alpha_path_c, ref_cor, incoptdepth_pol, deltau_pol, &
                     &  alpha_path_gl, ds_dz_gw, ct, stcp, stsp, t_script, &
                     &  do_dumps, prod_pol, tau_pol, rad_pol, p_stop )

    ! Polarized radiative transfer.  Radiances only, no derivatives.

    use CS_Expmat_m, only: CS_Expmat
    use DO_DELTA_M, ONLY: POLARIZED_PATH_OPACITY
    use Dump_0, only: Dump, Dump_2x2xN
    use GLNP, ONLY: Ng
    use MCRT_M, ONLY: MCRT
    use MLSKinds, only: RP
    use Opacity_m, only: Opacity

  ! inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    integer, intent(in) :: more_inds(:)      ! Places in the coarse path
  !                                            where GL is needed
    real(rp), intent(in) :: e_rflty          ! earth reflectivity value (0--1).
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    complex(rp), intent(in) :: alpha_path_c(-1:,:)  ! absorption coefficient
      !              on coarse grid.
    complex(rp), intent(inout) :: deltau_pol(:,:,:) ! 2 X 2 X path.  Incremental
      !              transmissivity on the coarse path. Called E in some notes.
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !              length ratios.
    complex(rp), intent(inout) :: incoptdepth_pol(:,:,:) ! incremental path
      !              opacities from one-sided layer calculation on output. it
      !              is the full integrated layer opacity. 2x2xPath
    complex(rp), intent(in) :: alpha_path_gl(-1:,:) ! absorption coefficient on
      !              gl grid.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: CT(:)            ! Cos theta          for Mag field
    real(rp), intent(in) :: STCP(:)          ! Sin theta Cos Phi  for Mag field
    real(rp), intent(in) :: STSP(:)          ! Sin theta Sin Phi  for Mag field
    real(rp), intent(in) :: T_script(:)      ! differential temperatures (K)
      !              on coarse path.
    integer, intent(in) :: Do_Dumps          ! Dump intermediate results if > 0

  ! outputs

    complex(rp), intent(out) :: prod_pol(:,:,:) ! product of E matrices. 2x2xPath
    complex(rp), intent(out) :: tau_pol(:,:,:)  ! transmission function. 2x2xPath
    complex(rp), intent(out) :: rad_pol(:,:)    ! radiance (K). 2x2.
    integer, intent(out) :: p_stop           ! path stop index if >= 0, else
      !              -index in incoptdepth_pol where cs_expmat failed.

  ! Internals

    real(rp), save :: E_Stop  = 1.0_rp ! X for which Exp(X) is too small to worry
    complex(rp) :: gl_delta_polarized(-1:1,size(gl_inds)/ng)
    complex(rp) :: incoptdepth_pol_gl(2,2,size(gl_inds)/ng)
    integer :: N_PATH
    integer :: Status ! from cs_expmat

  ! Begin code

    n_path = size(del_zeta)

    if ( e_stop > 0.0_rp ) e_stop = log(epsilon(0.0_rp)) ! only once

  ! see if anything needs to be gl-d

    if ( size(gl_inds) > 0 ) then

      call polarized_path_opacity ( del_zeta,    &
                 &  alpha_path_c, alpha_path_gl, &
                 &  ds_dz_gw,                    &
                 &  gl_delta_polarized, more_inds, gl_inds )

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

  end subroutine Rad_tran_Pol

!--------------------------------------------------  DRad_tran_df  -----
! This is the radiative transfer derivative wrt mixing ratio model

  subroutine DRad_tran_df ( max_f, gl_inds, del_zeta, Grids_f, eta_zxp,       & 
                          & do_calc_f, do_gl, del_s, ref_cor, ds_dz_gw,       & 
                          & inc_rad_path, dAlpha_df_c, dAlpha_df_f, i_start,  & 
                          & tan_pt, i_stop, LD, d_delta_df, nz_d_delta_df,    &
                          & nnz_d_delta_df, drad_df, dB_df, Tau, nz_zxp,      &
                          & nnz_zxp, alpha_path_c, Beta_c_e,                  &
                          & dBeta_c_a_dIWC, dBeta_c_s_dIWC, dTScat_df, W0 )

    use d_t_script_dtnp_m, only: dT_script
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP
    use SCRT_DN_M, ONLY: DSCRT_DT, DSCRT_DX
    use TScat_Support_m, only: Get_dB_df

! Inputs

    integer, intent(in) :: Max_f            ! Leading dimension of dAlpha_df_f
    integer, intent(in) :: gl_inds(:)       ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: eta_zxp(max_f,*) ! representation basis function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !                                        representation basis function is
      !                                        not zero.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: inc_rad_path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(in) :: dAlpha_df_c(:,:) ! On the coarse path
    real(rp), intent(in) :: dAlpha_df_f(max_f,*) ! On the GL path
    integer, intent(in) :: I_start           ! path_start_index + 1
    integer, intent(in) :: tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: I_stop            ! path stop index

    integer, intent(in) :: LD                ! Leading dimension of D_Delta_dF

    ! Optionals for TScat
    real(rp), intent(inout) :: dB_df(:) ! scratch, on the path, size=0 for no TScat
    real(rp), intent(in), optional :: Tau(:)
    ! To project dB_df from the path to the grid:
    integer, intent(in), optional :: NZ_ZXP(:,:) ! for eta_zxp: path X 1 SV
    integer, intent(in), optional :: NNZ_ZXP(:)  ! for eta_zxp: SV
    real(rp), intent(in), optional :: Alpha_path_c(:)
    real(rp), intent(in), optional :: Beta_c_e(:)
    real(rp), intent(in), optional :: dBeta_c_a_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in), optional :: dBeta_c_s_dIWC(:)  ! on the path, w.r.t. IWC on the path
    real(rp), intent(in), optional :: dTScat_df(:,:) ! Path X SV.  On the path, w.r.t.
                                                 ! f == Molecules on the grid
    real(rp), intent(in), optional :: W0(:)

! Outputs

    real(rp), intent(inout) :: d_delta_df(ld,*) ! path x sve.  derivative of
      !              delta wrt mixing ratio state vector element. (K)
      !              Initially set to zero by caller.
    integer, intent(inout), target :: nz_d_delta_df(:,:) ! Nonzeros in d_delta_df
    integer, intent(inout) :: nnz_d_delta_df(:) ! Column lengths in nz_delta_df
    real(rp), intent(out) :: drad_df(:)      ! derivative of radiances wrt
      !              mixing ratio state vector element. (K)

! Internals

    integer :: i_begin, sps_i, sv_i
    real(rp) :: d_delta_B_df(size(dB_df),1)
    logical Nothing(Grids_f%l_v(ubound(Grids_f%l_z,1))) ! "Nothing to do here"
    logical :: Do_TScat              ! Include dependence upon dB_df

! Begin code

    call get_all_d_delta_df ( max_f, gl_inds, del_zeta, Grids_f, eta_zxp,    & 
                            & do_calc_f, do_gl, del_s, ref_cor, ds_dz_gw,    & 
                            & dAlpha_df_c, dAlpha_df_f, LD, d_delta_df,      & 
                            & nz_d_delta_df, nnz_d_delta_df, nothing )

    Do_TScat = size(dB_df) > 0

    do sps_i = 1, ubound(Grids_f%l_z,1)

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

        call get_dB_df ( alpha_path_c, beta_c_e, dBeta_c_a_dIWC, &
                       & dBeta_c_s_dIWC, dAlpha_df_c(:,sps_i), w0, &
                       & grids_f%mol(sps_i), dB_df )

      end if

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

        if ( nothing(sv_i) ) then
          drad_df(sv_i) = 0.0
          cycle
        end if

        i_begin = max(i_start,min(nz_d_delta_df(1,sv_i),i_stop))

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

          call dt_script ( dB_df, eta_zxp(:,sv_i:sv_i), nz_zxp(:,sv_i:sv_i), &
            & nnz_zxp(sv_i:sv_i), d_delta_B_df, w0, dTScat_df(:,sv_i:sv_i) )

          !{ Now do the path integration
          !  $\sum \frac{\partial}{\partial f^k_{lm}}
          !     \mathcal{T}_i \overline{\Delta B_i}$

          call dscrt_dt ( tan_pt, d_delta_df(:,sv_i), tau, inc_rad_path,&
                        & d_delta_B_df(:,1), i_begin, i_stop, drad_df(sv_i) )
        else
          call dscrt_dx ( tan_pt, d_delta_df(:,sv_i), inc_rad_path, &
                       &  i_begin, i_stop, drad_df(sv_i))
        end if

      end do ! sv_i

    end do ! sps_i

  end subroutine drad_tran_df

!------------------------------------------------  D2Rad_tran_df2  -----
! This is the radiative transfer second derivative wrt mixing ratio model

  subroutine D2Rad_tran_df2 ( max_f, gl_inds, del_zeta, Grids_f, eta_zxp, &
                            & do_calc_f, do_gl, del_s, ref_cor, ds_dz_gw, &
                            & inc_rad_path, d2Alpha_df2_c, d2Alpha_df2_f, &
                            & i_start, tan_pt, i_stop, LD, d_delta_df,    &
                            & nz_d_delta_df, nnz_d_delta_df,              &
                            & d2_delta_df2, d2rad_df2 )

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP
    use SCRT_DN_M, ONLY: D2SCRT_DX2

! Inputs

    integer, intent(in) :: Max_f             ! Leading dimension of d2Alpha_df2_f
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid. This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: eta_zxp(max_f,*)   ! representation basis function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !                                        representation basis function is
      !                                        not zero.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: del_s(:)        ! unrefracted path length.
    real(rp), intent(in) :: ref_cor(:)      ! refracted to unrefracted path
      !                                       length ratios.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative *
      !              gw on the entire grid. Only the gl_inds part is used.
    real(rp), intent(in) :: inc_rad_path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(in) :: d2Alpha_df2_c(:,:) ! On the coarse path
    real(rp), intent(in) :: d2Alpha_df2_f(max_f,*) ! On the GL path
    integer, intent(in) :: I_start           ! path_start_index + 1
    integer, intent(in) :: tan_pt            ! Tangent point index
    integer, intent(in) :: I_stop            ! path stop index
    integer, intent(in) :: LD                ! Leading dimension of D_Delta_dF
    ! integer, intent(in) :: min_nz_d_delta_df(:) ! First Nonzeros in d_delta_df

    real(rp), intent(inout) :: d_delta_df(ld,*) ! path x sve. derivative of
      !              delta wrt mixing ratio state vector element. (K)
      !              Initially set to zero by caller.




! Outputs

    integer, intent(inout), target :: nz_d_delta_df(:,:) ! Nonzeros in d_delta_df
    integer, intent(inout) :: nnz_d_delta_df(:) ! Column lengths in nz_delta_df
    real(rp), intent(inout) :: d2_delta_df2(:,:,:) ! path x sve x sve.  Second 
      !               derivative of delta wrt mixing ratio state vector element.
    real(rp), intent(out) :: d2rad_df2(:,:)    ! second derivative of radiances wrt
                                               ! mixing ratio state vector element. (K)

! Internals

    logical :: Nothing(Grids_f%l_v(ubound(Grids_f%l_z,1))) ! "Nothing to do here"
    integer :: i_begin
    integer :: sps_i, sps_j          ! species indices
    integer :: q, r                  ! state vector indices

! Begin code

    call get_all_d2_delta_df2( max_f, gl_inds, del_zeta, Grids_f, eta_zxp,   &
                             & do_calc_f, do_gl, del_s, ref_cor, ds_dz_gw,   &
                             & d2Alpha_df2_c, d2Alpha_df2_f, nz_d_delta_df,  &
                             & nnz_d_delta_df, d2_delta_df2, nothing )

    do sps_i = 1, ubound(Grids_f%l_z,1)

      do sps_j = 1, ubound(Grids_f%l_z,1)
    
        ! do q = 1,  Grids_f%l_v(sps_i)
        do q = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

          if ( nothing(q) ) then
            d2rad_df2(q,:) = 0.0
            d2rad_df2(:,q) = 0.0
            cycle
          end if

          ! do r = 1,  Grids_f%l_v(sps_j)
          do r = Grids_f%l_v(sps_j-1)+1, Grids_f%l_v(sps_j)

            i_begin = max(i_start,min(max(nz_d_delta_df(1,q),nz_d_delta_df(1,r)),i_stop))            
            !i_begin = max(i_start,min(max(min_nz_d_delta_df(r),min_nz_d_delta_df(q)),i_stop))
            ! i_begin = max(i_start,min(max(inds_q(1),inds_r(1)),i_stop))
            
            
            call d2scrt_dx2 ( tan_pt, d_delta_df(:,q), d_delta_df(:,r), d2_delta_df2(:,q,r), &
                            & inc_rad_path, i_begin, i_stop, d2rad_df2(q,r) )

          end do ! r

        end do ! q
      
      end do ! sps_j

    end do ! sps_i

  end subroutine d2rad_tran_df2

!-------------------------------------------- Get_all_d2_delta_df2 -----

  subroutine Get_all_d2_delta_df2 ( max_f, gl_inds, del_zeta, Grids_f,       &
                              & eta_zxp, do_calc_f, do_gl, del_s, ref_cor,   &
                              & ds_dz_gw, d2Alpha_df2_c, d2Alpha_df2_f,      &
                              & nz_d_delta_df, nnz_d_delta_df, d2_delta_df2, &
                              & nothing )

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP

! Inputs

    integer, intent(in) :: Max_f             ! Leading dimension of d2Alpha_df2_f
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: eta_zxp(max_f,*) ! representation basis function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !                                        representation basis function is
      !                                        not zero.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: d2Alpha_df2_c(:,:)  ! On the coarse path
    real(rp), intent(in) :: d2Alpha_df2_f(max_f,*)  ! On the GL path

! Outputs

    integer, intent(inout), target :: nz_d_delta_df(:,:)  ! Nonzeros in d_delta_df
    integer, intent(inout) :: nnz_d_delta_df(:)           ! Column lengths in nz_delta_df

    ! IGOR: Compare it with declaration  d_delta_df(ld,*)  in Get_all_d_delta_df .
    real(rp), intent(inout) :: d2_delta_df2(:,:,:) ! path x sve x sve.  Second Derivative
      !              of delta wrt mixing ratio state vector elements. (K)
      !              Initially set to zero by caller.
    logical, intent(out) :: nothing(:) ! "Nothing to do for this s.v. element

! Internals

    integer :: n_inds_q, n_inds_r
    integer :: no_to_gl_q, no_to_gl_r
    integer :: sps_i, sps_j          ! species indices
    integer :: sps_n
    integer :: q, r                  ! state vector indices: sv_i, sv_j
    integer :: diracDelta            !   =1 if q=r;  =0 otherwise
    integer, target, dimension(1:size(del_s)) ::  all_inds_B_q,  all_inds_B_r
    integer, target, dimension(1:size(del_s)) :: more_inds_B_q, more_inds_B_r
    integer, pointer :: all_inds_q(:), all_inds_r(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, pointer :: inds_q(:), inds_r(:)      ! inds => part_of_nz_d_delta_df;
                                     ! Indices on coarse path where do_calc.
    integer, pointer :: more_inds_q(:), more_inds_r(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.

    real(rp) :: singularity(1:size(del_s)) ! integrand on left edge of coarse
                                     ! grid panel -- singular at tangent pt.
    logical :: do_calc_q(1:size(del_s)) ! Flags on coarse path where do_calc_c
                                     ! or (do_gl and any corresponding
                                     ! do_calc_f).
    logical :: do_calc_r(1:size(del_s)) ! Flags on coarse path where do_calc_c
                                     ! or (do_gl and any corresponding
                                     ! do_calc_f).

! Begin code

    ! d2_delta_df2 is set to zero by the caller outside of all its loops.

    sps_n = ubound(Grids_f%l_z,1)

    do sps_i = 1, sps_n

      do sps_j = 1, sps_n

        do q = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

          !d_delta_df(nz_d_delta_df(:nnz_d_delta_df(q),q),q) = 0.0
          nnz_d_delta_df(q) = 0

          do r = Grids_f%l_v(sps_j-1)+1, Grids_f%l_v(sps_j)

            !d_delta_df(nz_d_delta_df(:nnz_d_delta_df(r),r),r) = 0.0
            nnz_d_delta_df(r) = 0

            ! Skip the masked derivatives, according to the l2cf inputs
      
            nothing(q) = .not. Grids_f%deriv_flags(q)
            nothing(r) = .not. Grids_f%deriv_flags(r)
            if ( nothing(q) .or. nothing(r) )   cycle
        
            ! find where the non zeros are along the path (for q)

            call get_do_calc_indexed ( size(do_gl), do_calc_f(:,q), &
              & gl_inds, do_gl, do_calc_q, n_inds_q, nz_d_delta_df(:,q) )
            
            nnz_d_delta_df(q) = n_inds_q
            nothing(q) = n_inds_q == 0
            if ( nothing(q) ) cycle

            inds_q => nz_d_delta_df(1:n_inds_q,q)

            no_to_gl_q = count(do_gl(inds_q))

            all_inds_q =>  all_inds_B_q(1:no_to_gl_q)
            more_inds_q => more_inds_B_q(1:no_to_gl_q)

            ! see if anything needs to be gl-d (for q)
            if ( no_to_gl_q > 0 ) &
              & call get_inds ( do_gl, do_calc_q, more_inds_q, all_inds_q )

            !
            ! find where the non zeros are along the path (for r)
            !
            call get_do_calc_indexed ( size(do_gl), do_calc_f(:,r), &
              & gl_inds, do_gl, do_calc_r, n_inds_r, nz_d_delta_df(:,r) )

            nnz_d_delta_df(r) = n_inds_r
            nothing(r) = n_inds_r == 0
            if ( nothing(r) ) cycle

            inds_r => nz_d_delta_df(1:n_inds_r,r)

            no_to_gl_r = count(do_gl(inds_r))

            all_inds_r =>  all_inds_B_r(1:no_to_gl_r)
            more_inds_r => more_inds_B_r(1:no_to_gl_r)

            if ( no_to_gl_r > 0 ) &
            ! see if anything needs to be gl-d (for r)
              & call get_inds ( do_gl, do_calc_r, more_inds_r, all_inds_r )

     ! IGOR - May not need to enter this subroutine for NON log basis sps, 
     ! since d2_delta_df2 for linear basis is zero.

            ! For molecules in logarithmic basis, calculate d2_delta_df2:

            if ( grids_f%lin_log(sps_i) ) then

              if( sps_i == sps_j ) then    ! otherwise, d2_delta_df2 = 0

                ! For same species, the following quantities should be the same
                ! for different sve:
                !   inds, all_inds, more_inds, sps.  Thus, only q quantities
                ! are passed.

                if ( q == r ) then
                  diracDelta = 1.0
                else
                  diracDelta = 0.0
                end if

                call get_d2_delta_df2( diracDelta, inds_q, gl_inds,        &
                  & all_inds_q, more_inds_q, eta_zxp(:,q), eta_zxp(:,r),   &
                  & d2Alpha_df2_c(:,sps_i), d2Alpha_df2_f(:,sps_i), del_s, &
                  & del_zeta, ds_dz_gw, grids_f%values(q), grids_f%values(r), &
                  & singularity, d2_delta_df2(:,q,r), ref_cor )

              end if   ! sps_i == sps_j

            end if   ! lin_log

          end do   ! r

        end do   ! q

      end do   ! sps_j

    end do   ! sps_i

  end subroutine Get_all_d2_delta_df2

!--------------------------------------------------  drad_tran_dt  -----
! This is the radiative transfer derivative wrt temperature model

  subroutine DRad_tran_dt ( gl_inds, del_zeta, h_path_c, dh_dt_path_c,     &
                         &  alpha_path_c, dAlpha_dT_path_c, eta_zxp,       &
                         &  do_calc_t_c, do_calc_hyd_c, del_s, ref_cor,    &
                         &  h_tan, dh_dt_tan, do_gl, h_path_f, t_path_f,   &
                         &  dh_dt_path_f, alpha_path_f, dAlpha_dT_path_f,  &
                         &  do_calc_t_f, ds_dh, dh_dz_gw, ds_dz_gw,        &
                         &  dt_scr_dt, tau, inc_rad_path, i_start, tan_pt, &
                         &  i_stop, deriv_flags, pfa_update, drad_dt )

    use GLNP, only: NG
    use MLSKinds, only: RP
    use SCRT_DN_M, ONLY: DSCRT_DT, DSCRT_DX

! Inputs

    integer, intent(in) :: gl_inds(:)   ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: h_path_c(:)     ! path heights + req on main grid km.
    real(rp), intent(in) :: dh_dt_path_c(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on main grid.
    real(rp), intent(in) :: alpha_path_c(:) ! path absorption(km^-1)
!                                             on main grid.
    real(rp), intent(in) :: dAlpha_dT_path_c(:) ! path dAlpha/dT on main grid
    real(rp), intent(in) :: eta_zxp(:,:)    ! representation basis function
!                                              combined grid.
    logical, intent(in) :: do_calc_t_c(:,:) ! Indicates where the
!                    representation basis function is not zero on main grid.
    logical, intent(in) :: do_calc_hyd_c(:,:) ! Indicates where dh_dt is not
!                                             zero on main grid.
    real(rp), intent(in) :: del_s(:)        ! unrefracted path length.
    real(rp), intent(in) :: ref_cor(:)      ! refracted to unrefracted path
!                                             length ratios.
    real(rp), intent(in) :: h_tan           ! tangent height + req (km).
    real(rp), intent(in) :: dh_dt_tan(:)    ! derivative of path height wrt
!                                             temperature at the tangent (km/K).
    logical, intent(in) :: do_gl(:)         ! Indicates where on the coarse path
!                                             to do gl integrations.
    real(rp), intent(in) :: h_path_f(:)     ! path heights + req on gl grid km.
    real(rp), intent(in) :: t_path_f(:)     ! path temperature(K) on gl grid.
    real(rp), intent(in) :: dh_dt_path_f(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on gl grid.
    real(rp), intent(in) :: alpha_path_f(:) ! path absorption(km^-1) on gl grid.
    real(rp), intent(in) :: dAlpha_dT_path_f(:) ! path dAlpha/dT on gl grid
    logical, intent(in) :: do_calc_t_f(:,:) ! Indicates where the
!                    representation basis function is not zero on gl grid.
    real(rp), intent(in) :: ds_dh(:)        ! path length wrt height derivative
!                                             on complete grid.  Only the
!                                             gl_inds part is used.
    real(rp), intent(in) :: dh_dz_gw(:)     ! path height wrt zeta derivative * gw
!                                             on complete grid.  Only the
!                                             gl_inds part is used.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative * gw
!                                             on complete grid.  Only the
!                                             gl_inds part is used.
    real(rp), intent(in) :: dt_scr_dt(:,:)  ! d t_script / d T * d T / d eta.
    real(rp), intent(in) :: tau(:)          ! transmission function.
    real(rp), intent(in) :: inc_rad_path(:) ! incremental radiance along the
                                            ! path.  t_script * tau.
    integer, intent(in) :: i_start          ! path start index + 1
    integer, intent(in) :: Tan_pt           ! Tangent point index in Del_Zeta
    integer, intent(in) :: i_stop           ! path stop index
    logical, intent(in) :: deriv_flags(:)   ! Indicates which temperature
!                                             derivatives to do
    logical, intent(in) :: PFA_Update       ! Use DSCRT_DX instead of DSCRT_DT.

! Output
    real(rp), intent(out) :: drad_dt(:)     ! derivative of radiances wrt
!                                             temperature state vector
!                                             element. (K)

! Internals

    integer :: A, B, GA
    integer :: i, i_begin, n_inds, n_path, no_to_gl, p_i, sv_i
    integer, target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer, target, dimension(1:size(inc_rad_path)) :: inds_B, more_inds_B
    integer, pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                     ! on coarse path where do_calc.
    integer, pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.

    real(rp) :: d_delta_dt(size(del_s,1),size(eta_zxp,2)) ! path x sve.
      ! derivative of delta (incremental opacity) wrt temperature. (K)

    real(rp) :: fa, fb
    real(rp) :: S_DEl_S                  ! Running sum of Del_S
    real(rp) :: singularity(1:size(del_zeta)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.

    logical :: do_calc(1:size(del_zeta)) ! do_calc_t_c .or. ( do_gl .and. any
                                         ! of the corresponding do_calc_t_f ).
    logical :: NeedFA                    ! Need F(A) for hydrostatic

! Begin code

    n_path = size(del_zeta)

! compute the opacity derivative singularity value

    d_delta_dt = 0.0_rp

    do sv_i = 1, size(eta_zxp,dim=2)
      drad_dt(sv_i) = 0.0
      if ( .not. deriv_flags(sv_i)) cycle
      i_begin = i_start

! do the absorption part
! combine non zeros flags for both the main and gl parts

      call get_do_calc ( do_calc_t_c(:,sv_i), do_calc_t_f(:,sv_i), do_gl, &
        & do_calc, n_inds, inds_B )

      if ( n_inds > 0 ) then

        inds => inds_B(1:n_inds)

        no_to_gl = count(do_gl(inds))

        all_inds => all_inds_B(1:no_to_gl)
        more_inds => more_inds_B(1:no_to_gl)

        ! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

        call get_d_delta_dx ( inds, gl_inds, all_inds, more_inds, &
          & eta_zxp(:,sv_i), dAlpha_dT_path_c, dAlpha_dT_path_f,  &
          & del_s, del_zeta, ds_dz_gw, singularity, &
          & d_delta_dt(:,sv_i) ) ! No ref_cor yet

        i_begin = max(inds(1)-1, i_start)

      end if ! n_inds > 0

! now do the hydrostatic part
! combine boundaries flags

      do_calc = do_calc_hyd_c(:,sv_i)
      do_calc(2:tan_pt) =          do_calc(tan_pt)   .or. do_calc(2:tan_pt)          .or. do_calc(1:tan_pt-1)
      do_calc(tan_pt+1:n_path-1) = do_calc(tan_pt+1) .or. do_calc(tan_pt+1:n_path-1) .or. do_calc(tan_pt+2:n_path)

! since this is a layer boundary calculation we must require

      do_calc((/1,n_path/)) = .false.

! find where the non zeros are along the path

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
              fa = (h_path_c(p_i-1) * dh_dt_path_c(p_i-1,sv_i) - &
                 &  h_tan * dh_dt_tan(sv_i)) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s - del_s(p_i)
            fb = (h_path_c(p_i) * dh_dt_path_c(p_i,sv_i) - &
                & h_tan * dh_dt_tan(sv_i)) / s_del_s
            inds(i) = p_i
            d_delta_dt(p_i,sv_i) = d_delta_dt(p_i,sv_i) + &
              &                    alpha_path_c(p_i) * (fa - fb)
            fa = fb
            i = i + 1
          else
            s_del_s = s_del_s - del_s(p_i)
          end if
        end do ! p_i

! special processing at tangent.  fb is zero

        if ( do_calc(tan_pt) ) then
          d_delta_dt(tan_pt,sv_i) = d_delta_dt(tan_pt,sv_i) + alpha_path_c(tan_pt) * fa
          inds(i) = tan_pt
          i = i + 1
        end if

        needFA = .not. do_calc(tan_pt+1)
        s_del_s = del_s(tan_pt+1)
        if ( do_calc(tan_pt+1) ) then
          fa = (h_path_c(tan_pt+2) * dh_dt_path_c(tan_pt+2,sv_i) - &
              & h_tan * dh_dt_tan(sv_i)) / s_del_s
          d_delta_dt(tan_pt+1,sv_i) = d_delta_dt(tan_pt+1,sv_i) + &
            &                      alpha_path_c(tan_pt+1) * fa
          inds(i) = tan_pt + 1
          i = i + 1
        end if

        ! Several subscripts in this loop are offset by 1 from the nearly-
        ! identical loop above, and we use (fb-fa) here instead of (fa-fb).
        do p_i = tan_pt + 2, n_path - 1
          if ( do_calc(p_i) ) then
            if ( needFA ) then ! only once in this loop
              fa = (h_path_c(p_i) * dh_dt_path_c(p_i,sv_i) - &
                 &  h_tan * dh_dt_tan(sv_i)) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s + del_s(p_i)
            fb = (h_path_c(p_i+1)*dh_dt_path_c(p_i+1,sv_i) - &
               &  h_tan * dh_dt_tan(sv_i)) / s_del_s
            inds(i) = p_i
            d_delta_dt(p_i,sv_i) = d_delta_dt(p_i,sv_i) + &
              &                    alpha_path_c(p_i) * (fb - fa)
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
            d_delta_dt(p_i,sv_i) = d_delta_dt(p_i,sv_i) +            &
              & del_zeta(p_i) *                                      &
              &  sum( ( alpha_path_f(a:b-1) - alpha_path_c(p_i) ) *  &
              &  (((2.0_rp*h_path_f(a:b-1)**2 - 3.0_rp*h_tan**2)     &     
              &    * dh_dt_path_f(a:b-1,sv_i) +                      &     
              &    h_path_f(a:b-1) * h_tan * dh_dt_tan(sv_i)) /      &     
              &   (sqrt(h_path_f(a:b-1)**2 - h_tan**2))**3           &     
              &   + eta_zxp(ga:ga+ng-1,sv_i) * ds_dh(ga:ga+ng-1) /   &     
              &   t_path_f(a:b-1)) * dh_dz_gw(ga:ga+ng-1) )
            a = b
          end if
        end do ! p_i

        i_begin = min(i_begin,inds(1))

      end if ! n_inds for hydrostatic > 0

! correct for path length refraction

      d_delta_dt(:,sv_i) = ref_cor(:) * d_delta_dt(:,sv_i)

! Accumulate the incremental opacity derivatives to get drad_dt

      if ( PFA_update ) then
        ! If we're doing a PFA update, we do not want to include
        ! dt_scr_dt again.

        call dscrt_dx ( tan_pt, d_delta_dt(:,sv_i), inc_rad_path, i_begin, i_stop, &
                     &  drad_dt(sv_i) )

      else

        call dscrt_dt ( tan_pt, d_delta_dt(:,sv_i), tau, inc_rad_path,&
                      & dt_scr_dt(:,sv_i),  i_begin, i_stop, drad_dt(sv_i) )

      end if

    end do ! sv_i

  end subroutine DRad_tran_dt

!--------------------------------------------------  drad_tran_dx  -----
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )

  subroutine DRad_tran_dx ( gl_inds, del_zeta, Grids_f, eta_zxp, sps_path,  &
                         &  sps_map, do_calc_f, dbeta_path_c, dbeta_path_f, &
                         &  do_gl, del_s, ref_cor, ds_dz_gw, inc_rad_path,  &
                         &  tan_pt, i_stop, drad_dx )

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP
    use SCRT_DN_M, ONLY: DSCRT_DX

! Inputs

    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indicies
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: eta_zxp(:,:)     ! representation basis function,
      !                                        composite path.
    real(rp), intent(in) :: sps_path(:,:)    ! Path species function, path X species.
    integer, intent(in) :: sps_map(:)        ! second-dimension subscripts for sps_path.
    logical, intent(in) :: do_calc_f(:,:)    ! Where the representation basis
      !                                        function is not zero, composite
      !                                        path.
    real(rp), intent(in) :: dbeta_path_c(:,:) ! derivative of beta wrt dx
      !                                        on main grid.
    real(rp), intent(in) :: dbeta_path_f(:,:) ! derivative of beta wrt dx
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: inc_rad_path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    integer, intent(in) :: Tan_pt            ! Tangent point index in inc_rad_path
    integer, intent(in) :: i_stop            ! path stop index

! Outputs

    real(rp), intent(out) :: drad_dx(:)      ! derivative of radiances wrt x
!                                              state vector element. (K)
! Internals

    integer :: n_inds, no_to_gl, sps_i, sps_m, sps_n, sv_i
    integer, target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer, target, dimension(1:size(inc_rad_path)) :: inds_B, more_inds_B
    integer, pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                     ! on coarse path where do_calc.
    integer, pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.

    real(rp) :: d_delta_dx(1:size(inc_rad_path))  ! derivative of delta
      !              wrt spectroscopy parameter. (K)

    real(rp) :: singularity(1:size(inc_rad_path)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.

    logical :: do_calc(1:size(inc_rad_path))      ! Flags on coarse path where
                                         ! do_calc_c or (do_gl and any
                                         ! corresponding do_calc_f).

! Begin code

    sps_n = ubound(grids_f%l_z,1)

    do sps_i = 1 , sps_n
      sps_m = sps_map(sps_i)

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

! find where the non zeros are along the path

        call get_do_calc_indexed ( size(do_gl), do_calc_f(:,sv_i), &
          & gl_inds, do_gl, do_calc, n_inds, inds_B )

        d_delta_dx = 0.0_rp

        if ( n_inds == 0 ) then
           drad_dx(sv_i) = 0.0
           cycle
        end if


        inds => inds_B(1:n_inds)

        no_to_gl = count(do_gl(inds))

        all_inds => all_inds_B(1:no_to_gl)
        more_inds => more_inds_B(1:no_to_gl)

! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

        ! We're not really computing d_delta_df for lin-log mixing
        ! ratio.  It turns out that get_d_delta_df_linlog does the
        ! correct computation if we substitute dbeta_path for beta_path.
        call get_d_delta_df_linlog ( inds, gl_inds, all_inds, more_inds, &
          & eta_zxp(:,sv_i), sps_path(:,sps_m), dbeta_path_c(:,sps_i),   &
          & dbeta_path_f(:,sps_i), del_s, del_zeta, ds_dz_gw, ref_cor,   &
          & grids_f%values(sv_i), singularity, d_delta_dx )

        call dscrt_dx ( tan_pt, d_delta_dx, inc_rad_path, &
                     &  1, i_stop, drad_dx(sv_i))

      end do

    end do

  end subroutine DRad_tran_dx

!--------------------------------------------  Get_all_d_delta_df  -----

  subroutine Get_all_d_delta_df ( max_f, gl_inds, del_zeta, Grids_f, eta_zxp, & 
                                & do_calc_f, do_gl, del_s, ref_cor, ds_dz_gw, & 
                                & dAlpha_df_c, dAlpha_df_f, LD, d_delta_df,   & 
                                & nz_d_delta_df, nnz_d_delta_df, nothing )

    !{ Compute
    !  \begin{equation}
    !  \frac{\partial \delta_{i \rightarrow i-1}}{\partial f^k_{lm}} =
    !  \int_{\zeta_i}^{\zeta_i-1} \frac{\partial \alpha(s)}{\partial f^k(s)}
    !  \eta^k_{lm}(s) \,\text{d}s
    !  \end{equation}
    !  where $k$ is a species index, and $lm$ index $(\phi^k_l,\zeta^k_m)$.
    !  The second dimensions of {\tt dAlpha_df_c, dAlpha_df_f, d_delta_df}, and
    !  {\tt nz_d_delta_df} flatten out $(k,l,m)$ to a single index.  The complication
    !  here arises because $\eta^k_{lm}(s)$ is very sparse.

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP

! Inputs

    integer, intent(in) :: Max_f            ! Leading dimension of dAlpha_df_f
    integer, intent(in) :: gl_inds(:)       ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: eta_zxp(max_f,*) ! representation basis function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !                                        representation basis function is
      !                                        not zero.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: dAlpha_df_c(:,:) ! On the coarse path
    real(rp), intent(in) :: dAlpha_df_f(max_f,*) ! On the GL path

    integer, intent(in) :: LD                ! Leading dimension of D_Delta_dF

! Outputs

    real(rp), intent(inout) :: d_delta_df(ld,*) ! path x sve.  derivative of
      !              delta wrt mixing ratio state vector element. (K)
      !              Initially set to zero by caller.
    integer, intent(inout), target :: nz_d_delta_df(:,:) ! Nonzeros in d_delta_df
    integer, intent(inout) :: nnz_d_delta_df(:) ! Column lengths in nz_delta_df
    logical, intent(out) :: nothing(:)      ! "Nothing to do for this s.v. element

! Internals

    integer :: n_inds, no_to_gl, sps_i, sv_i
    integer, target, dimension(1:size(del_s)) :: all_inds_B
    integer, target, dimension(1:size(del_s)) :: more_inds_B
    integer, pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, pointer :: inds(:)      ! inds => part_of_nz_d_delta_df;
                                     ! Indices on coarse path where do_calc.
    integer, pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.

    real(rp) :: singularity(1:size(del_s)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.
    logical :: do_calc(1:size(del_s))    ! Flags on coarse path where do_calc_c
                                         ! or (do_gl and any corresponding
                                         ! do_calc_f).

! Begin code

    ! d_delta_df is set to zero by the caller outside of all its loops.
    ! We keep track of where we create nonzeros, and replace them by zeros
    ! on the next call.  This is done because the vast majority of
    ! d_delta_df elements are zero.

    do sps_i = 1, ubound(Grids_f%l_z,1)

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

        ! Everything in d_delta_df not indexed by nz_d_delta_df is already zero
        d_delta_df(nz_d_delta_df(:nnz_d_delta_df(sv_i),sv_i),sv_i) = 0.0
        nnz_d_delta_df(sv_i) = 0

        ! Skip the masked derivatives, according to the l2cf inputs

        nothing(sv_i) = .not. Grids_f%deriv_flags(sv_i)
        if ( nothing(sv_i) ) cycle

        ! find where the non zeros are along the path

        call get_do_calc_indexed ( size(do_gl), do_calc_f(:,sv_i), &
          & gl_inds, do_gl, do_calc, n_inds, nz_d_delta_df(:,sv_i) )

        nnz_d_delta_df(sv_i) = n_inds
        nothing(sv_i) = n_inds == 0
        if ( nothing(sv_i) ) cycle

        inds => nz_d_delta_df(1:n_inds,sv_i)

        no_to_gl = count(do_gl(inds))

        all_inds => all_inds_B(1:no_to_gl)
        more_inds => more_inds_B(1:no_to_gl)

        ! see if anything needs to be gl-d
        if ( no_to_gl > 0 ) &
          & call get_inds ( do_gl, do_calc, more_inds, all_inds )

        !{ Get d_delta_df for one state-vector element.  This is
        !  $\frac{\partial \delta_i}{\partial f^k_{lm}}$ where $i$ is the
        !  index of a path point, $k$ is the species index, and $lm$ are indices
        !  for $(\phi_l,\zeta_m)$.  {\tt sv_i} flattens $(k,l,m)$ to one index.
        call get_d_delta_df ( inds, gl_inds, all_inds, more_inds, &
          & eta_zxp(:,sv_i), dAlpha_df_c(:,sps_i), dAlpha_df_f(:,sps_i), &
          & del_s, del_zeta, ds_dz_gw, grids_f%lin_log(sps_i),    &
          & grids_f%values(sv_i), singularity, d_delta_df(:,sv_i), ref_cor )

      end do ! sv_i

    end do ! sps_i

  end subroutine Get_all_d_delta_df

  ! ------------------------------------------------  Get_Do_Calc  -----
  subroutine Get_Do_Calc ( Do_Calc_c, Do_Calc_f, Do_GL, Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_c or Do_GL and any of the corresponding Do-Calc_f
  ! flags are set.

    use GLNP, ONLY: Ng

    logical, intent(in) :: Do_Calc_c(:) ! On the coarse grid
    logical, intent(in) :: Do_Calc_f(:) ! On the GL grid
    logical, intent(in) :: Do_GL(:)     ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)  ! Where on coarse grid to do calc.
    integer, intent(out), optional :: N_Inds  ! count(do_calc)
    integer, intent(out), optional :: Inds(:) ! Indices where do_calc is true

    integer :: I, P_I

    i = 1
    do_calc = do_calc_c
    do p_i = 1 , size(do_gl)
      if ( do_gl(p_i) ) then
        do_calc(p_i) = do_calc(p_i) .or. any(do_calc_f(i:i+ng-1))
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

  ! .............................................  Get_d_delta_dx  .....
  subroutine Get_d_delta_dx ( Inds, GL_Inds, All_inds, More_inds, eta_zxp, &
    & dAlpha_dx_path_c, dAlpha_dx_path_f, Del_s, Del_Zeta, ds_dz_gw,       &
    & Singularity, d_delta_dx, Ref_cor )

    ! Get d_delta_dx.  For species for which beta does not depend upon
    ! mixing ratio this gets d_delta_df if dAlpha_dx_path_* is beta_path_*.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:)       ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)    ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)   ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer, intent(in) :: More_inds(:)  ! Indices on the coarse path where
                                         ! GL corrections get applied.
    real(rp), intent(in) :: eta_zxp(*)   ! representation basis function.
    real(rp), intent(in) :: dAlpha_dx_path_c(*) ! on coarse grid.
    real(rp), intent(in) :: dAlpha_dx_path_f(*) ! cross section on GL grid.
    real(rp), intent(in) :: Del_s(:)     ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:)  ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:)  ! ds/dh * dh/dz * GL weights
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_delta_dx(:) ! Derivative of delta.
                               ! intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: Ref_cor(:) ! refracted to unrefracted
                                           !  path length ratios.

    integer :: AA, GA, I, II

    do i = 1, size(inds)
      ii = inds(i)
      singularity(ii) = dAlpha_dx_path_c(ii) * eta_zxp(ii*ngp1-ng)
      d_delta_dx(ii) = singularity(ii) * del_s(ii)
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
    !  second integral, $G(\zeta)$ is {\tt dAlpha_dx_path\_f * eta\_zxp\_f
    !  * sps\_path} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}. The weights  are
    !  {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_dx(ii) = d_delta_dx(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_zxp(ga:ga+ng-1) * dAlpha_dx_path_f(aa:aa+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    if ( present(ref_cor) ) d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)

  end subroutine Get_d_delta_dx

  ! .............................................  Get_d_delta_df  .....
  subroutine Get_d_delta_df ( Inds, GL_Inds, All_inds, More_inds, eta_zxp, &
    & dAlpha_df_path_c, dAlpha_df_path_f, Del_s, Del_Zeta, ds_dz_gw, lin_log, grids_v, &
    & Singularity, d_delta_dx, Ref_cor )

    ! Get d_delta_dx.  For species for which beta does not depend upon
    ! mixing ratio this gets d_delta_df if dAlpha_dx_path_* is beta_path_*.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:)       ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)    ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)   ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer, intent(in) :: More_inds(:)  ! Indices on the coarse path where
                                         ! GL corrections get applied.
    real(rp), intent(in) :: eta_zxp(*)   ! representation basis function on
                                         ! fine grid.
    real(rp), intent(in) :: dAlpha_df_path_c(*) ! dAlpha_df on coarse grid.
    real(rp), intent(in) :: dAlpha_df_path_f(*) ! dAlpha_df on GL grid.
    real(rp), intent(in) :: Del_s(:)     ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:)  ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:)  ! ds/dh * dh/dz * GL weights
    logical, intent(in) :: lin_log       ! logarithmic interpolation was used
    real(rp), intent(in) :: grids_v      ! Grids_f%values(sv_i)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d_delta_dx(:) ! Derivative of delta.
                               ! intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: Ref_cor(:) ! refracted to unrefracted
                                           !  path length ratios.

    integer :: AA, GA, I, II

    do i = 1, size(inds)
      ii = inds(i)
      singularity(ii) = dAlpha_df_path_c(ii) * eta_zxp(ii*ngp1-ng)
      d_delta_dx(ii) = singularity(ii) * del_s(ii)
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
    !  second integral, $G(\zeta)$ is {\tt dAlpha_dx_path\_f * eta\_zxp\_f
    !  * sps\_path} -- which have already been evaluated at the appropriate
    !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}. The weights  are
    !  {\tt gw}.

    do i = 1, size(all_inds)
      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)
      d_delta_dx(ii) = d_delta_dx(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_zxp(ga:ga+ng-1) * dAlpha_df_path_f(aa:aa+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    if ( present(ref_cor) ) d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)

    ! Logarithmic interpolation correction
    if ( lin_log ) d_delta_dx(inds) = d_delta_dx(inds) * exp(-grids_v)

  end subroutine Get_d_delta_df

  ! ...........................................  Get_d_delta_df_f  .....
  subroutine Get_d_delta_df_f ( Inds, GL_Inds, All_inds, More_inds, Eta_zxp, &
    & Sps_path, Beta_path_c, Beta_path_f, dBeta_df_c, dBeta_df_f, Del_s,     &
    & Del_Zeta, ds_dz_gw, Ref_cor, Singularity, d_delta_df )

    ! Get d_delta_df for the case of species for which beta
    ! depends upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: Inds(:) ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: Eta_zxp(*)  ! representation basis function.
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
      singularity(ii) = eta_zxp(iii) &
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
        & sum( (eta_zxp(ga:ga+ng-1) &
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
    & Eta_zxp, Sps_path, Beta_path_c, Beta_path_f, Del_s, Del_Zeta,      &
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
    real(rp), intent(in) :: Eta_zxp(*)  ! representation basis function.
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
      singularity(ii) = eta_zxp(iii) * sps_path(iii) * &
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
        & sum( (eta_zxp(ga:ga+ng-1) * sps_path(ga:ga+ng-1) &
             &  * beta_path_f(aa:aa+ng-1) - &
             &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds) * exp(-grids_v)

  end subroutine Get_d_delta_df_linlog

  ! ....................................  Get_d_delta_df_linlog_f  .....
  subroutine Get_d_delta_df_linlog_f ( Inds, GL_Inds, All_inds, More_inds, &
    & Eta_zxp, Sps_path, Beta_path_c, Beta_path_f, dBeta_df_c, dBeta_df_f, &
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
    real(rp), intent(in) :: Eta_zxp(*)  ! representation basis function.
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
      singularity(ii) = eta_zxp(iii) * sps_path(iii) * &
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
        & sum( (eta_zxp(ga:ga+ng-1) * sps_path(ga:ga+ng-1) &
             &  * ( beta_path_f(aa:aa+ng-1) + &
             &      sps_path(ga:ga+ng-1) * dBeta_df_f(aa:aa+ng-1) ) &
             &  - singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
    end do

    ! Refraction correction
    d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds) * exp(-grids_v)

  end subroutine Get_d_delta_df_linlog_f


! .............................................  Get_d2_delta_df2  .....
  subroutine Get_d2_delta_df2 ( diracDelta, Inds, GL_Inds, All_inds, &
    & More_inds, Eta_zxp_q, Eta_zxp_r, d2Alpha_df2_path_c, d2Alpha_df2_path_f, &
    & Del_s, Del_Zeta, ds_dz_gw, &
    & Grids_v_q, Grids_v_r, Singularity, d2_delta_df2, Ref_cor )

    ! Get d2_delta_df2 for the case of lin_log species for which beta
    ! does not depend upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: diracDelta   !   =1 if q=r;  =0 otherwise
    integer, intent(in) :: Inds(:)      ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: Eta_zxp_q(*), Eta_zxp_r(*)  ! representation basis function.
    real(rp), intent(in) :: d2Alpha_df2_path_c(*)  ! d2Alpha_df2 on coarse grid.
    real(rp), intent(in) :: d2Alpha_df2_path_f(*)  ! d2Alpha_df2 on GL grid.
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
                       !  main grid.  This is for the whole coarse path, not just
                       !  the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:) ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: Grids_v_q, Grids_v_r     ! Grids_f%values(sv_i),  
                                                     ! Grids_f%values(sv_j)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d2_delta_df2(:) ! Second Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: ref_cor(:)  ! refracted to unrefracted path
                                                  !  length ratios.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)

      ii = inds(i)
      iii = ii*ngp1 - ng

      singularity(ii) = d2Alpha_df2_path_c(ii) * eta_zxp_q(iii) * (eta_zxp_r(iii) - diracDelta)
      d2_delta_df2(ii) = singularity(ii) * del_s(ii)

    end do ! i


    do i = 1, size(all_inds)

      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)

      d2_delta_df2(ii) = d2_delta_df2(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_zxp_q(ga:ga+ng-1) * (eta_zxp_r(ga:ga+ng-1) - diracDelta) * &
             & d2Alpha_df2_path_f(aa:aa+ng-1) - singularity(ii)) &
             & * ds_dz_gw(ga:ga+ng-1) )

    end do

    ! Refraction correction
    if( present(ref_cor) ) d2_delta_df2(inds) = ref_cor(inds) * d2_delta_df2(inds)

    d2_delta_df2(inds) = d2_delta_df2(inds) * exp(-grids_v_q) * exp(-grids_v_r)

  end subroutine Get_d2_delta_df2


! ......................................  Get_d2_delta_df2_linlog  .....
  subroutine Get_d2_delta_df2_linlog ( diracDelta, Inds, GL_Inds, All_inds, &
    & More_inds, Eta_zxp_q, Eta_zxp_r, Sps_path, Beta_path_c, Beta_path_f,  &
    & Del_s, Del_Zeta, ds_dz_gw, Ref_cor, Grids_v_q, Grids_v_r, &
    & Singularity, d2_delta_df2 )

    ! Get d2_delta_df2 for the case of lin_log species for which beta
    ! does not depend upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: diracDelta   !   =1 if q=r;  =0 otherwise
    integer, intent(in) :: Inds(:)   ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: Eta_zxp_q(*), Eta_zxp_r(*)  ! representation basis function.
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
    real(rp), intent(in) :: Grids_v_q, Grids_v_r     ! Grids_f%values(sv_i)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d2_delta_df2(:) ! Second Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)

      ii = inds(i)
      iii = ii*ngp1 - ng

      singularity(ii) = eta_zxp_q(iii) * (eta_zxp_r(iii) - diracDelta) * &
                      & sps_path(iii) * beta_path_c(ii)
      d2_delta_df2(ii) = singularity(ii) * del_s(ii)

    end do ! i


    do i = 1, size(all_inds)

      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)

      d2_delta_df2(ii) = d2_delta_df2(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_zxp_q(ga:ga+ng-1) * (eta_zxp_r(ga:ga+ng-1) - diracDelta) * &
             & sps_path(ga:ga+ng-1) * beta_path_f(aa:aa+ng-1) - singularity(ii)) &
             & * ds_dz_gw(ga:ga+ng-1) )

    end do

    ! Refraction correction
    d2_delta_df2(inds) = ref_cor(inds) * d2_delta_df2(inds) * exp(-grids_v_q) * exp(-grids_v_r)

  end subroutine Get_d2_delta_df2_linlog



! =====     Private Procedures     =====================================

  ! ----------------------------------------  Get_Do_Calc_Indexed  -----
  subroutine Get_Do_Calc_Indexed ( N, Do_Calc_all, F_Inds, Do_GL, &
    & Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_All(1::ngp1) or Do_GL and any of the corresponding
  ! Do_Calc_All(f_inds) flags are set.
  ! Get_Do_Calc_Indexed determines that get_d_delta_df needs to integrate a
  ! panel if do_calc_all(1::ngp1) is true, which means the interpolating
  ! coefficient is nonzero, or if do_gl is true and do_calc_all(f_inds) is
  ! true, which means that GL was needed, even if an interpolating
  ! coefficient is nonzero.

    use GLNP, ONLY: Ng, NGP1

    integer, intent(in) :: N ! sizes on coarse grid

#if defined NAG
!   Assumed-shape arguments are slower than assumed size
    logical, intent(in) :: Do_Calc_all(*) ! On the entire path
    integer, intent(in) :: F_Inds(*)      ! Indices in Do_Calc_All for find grid
    logical, intent(in) :: Do_GL(*)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(*)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(*)       ! Indices where do_calc is true
#elif defined IFC
!   Contiguous assumed-shape arguments are no slower than assumed size
    logical, contiguous, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, contiguous, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
    logical, contiguous, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, contiguous, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer,             intent(out) :: N_Inds        ! count(do_calc)
    integer, contiguous, intent(out) :: Inds(:)       ! Indices where do_calc is true
#elif defined LF95
!   Assumed-shape arguments are faster than assumed size because of copying
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for find grid
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(:)       ! Indices where do_calc is true
#else
!   We don't know whether ssumed-shape or assumed size arguments are faster
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer, intent(in) :: F_Inds(:)      ! Indices in Do_Calc_All for fine grid
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
    logical, intent(in) :: Do_GL(:)          ! path flag indicating where to do
      !                                        gl integrations.
    logical, intent(in) :: Do_Calc(:)

  ! Outputs
    integer, intent(out) :: More_Inds(:)
    integer, intent(out) :: All_Inds(:)

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
