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
  public :: Get_Do_Calc
  private ::  Get_Do_Calc_Indexed, Get_Inds

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
                     &  alpha_path_c, ref_cor, incoptdepth, &
                     &  alpha_path_gl, ds_dz_gw, t_script, &
                     &  tau, inc_rad_path, rad, i_stop )

    use GLNP, only: NG
    use MLSKinds, only: RP, IP
    use SCRT_DN_M, ONLY: SCRT

  ! inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    integer(ip), intent(in) :: more_inds(:)  ! Places in the coarse path
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
    integer(ip), intent(out) :: i_stop       ! path stop index

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
    use MLSKinds, only: RP, IP
    use Opacity_m, only: Opacity

  ! inputs

    integer, intent(in) :: Tan_pt            ! Tangent point index in Del_Zeta
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indices
    integer(ip), intent(in) :: more_inds(:)  ! Places in the coarse path
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
    integer(ip), intent(out) :: p_stop       ! path stop index if >= 0, else
      !              -index in incoptdepth_pol where cs_expmat failed.

  ! Internals

    real(rp), save :: E_Stop  = 1.0_rp ! X for which Exp(X) is too small to worry
    complex(rp) :: gl_delta_polarized(-1:1,size(gl_inds)/ng)
    complex(rp) :: incoptdepth_pol_gl(2,2,size(gl_inds)/ng)
    integer(ip) :: N_PATH
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

  subroutine DRad_tran_df ( max_f, indices_c, gl_inds, del_zeta, Grids_f,  &
                         &  beta_path_c, eta_zxp_f, sps_path, do_calc_f,   &
                         &  beta_path_f, do_gl, del_s, ref_cor,            &
                         &  ds_dz_gw, inc_rad_path, i_start, tan_pt,       &
                         &  i_stop, LD, d_delta_df, nz_d_delta_df,         &
                         &  nnz_d_delta_df, drad_df )

    use GLNP, only: NG
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP, IP
    use SCRT_DN_M, ONLY: DSCRT_DX

! Inputs

    integer, intent(in) :: Max_f            ! Leading dimension of Beta_Path_f
    integer(ip), intent(in) :: indices_c(:) ! coarse grid indicies
    integer(ip), intent(in) :: gl_inds(:)   ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: beta_path_c(:,:) ! cross section for each species
      !                                        on coarse grid.
    real(rp), intent(in) :: eta_zxp_f(max_f,*)   ! representation basis function.
    real(rp), intent(in) :: sps_path(:,:)    ! Path species function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !                                        representation basis function is
      !                                        not zero.
    real(rp), intent(in) :: beta_path_f(max_f,*) ! cross section for each species
      !                                        on gl grid.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: inc_rad_path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    integer, intent(in) :: i_start           ! path_start_index + 1
    integer, intent(in) :: tan_pt            ! Tangent point index in Del_Zeta
    integer, intent(in) :: i_stop            ! path stop index

    integer, intent(in) :: LD                ! Leading dimension of D_Delta_dF

! Outputs

    real(rp), intent(inout) :: d_delta_df(ld,*) ! path x sve.  derivative of
      !              delta wrt mixing ratio state vector element. (K)
      !              Initially set to zero by caller.
    integer, intent(out), target :: nz_d_delta_df(:,:) ! Nonzeros in d_delta_df
    integer, intent(out) :: nnz_d_delta_df(:) ! Column lengths in nz_delta_df
    real(rp), intent(out) :: drad_df(:)      ! derivative of radiances wrt
      !              mixing ratio state vector element. (K)

! Internals

    integer(ip) :: AA, GA, I, II, III
    integer(ip) :: i_begin, n_inds, no_to_gl, npc, sps_i, sps_n, sv_i
    integer(ip), target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer(ip), target, dimension(1:size(inc_rad_path)) :: more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_nz_d_delta_df;
                                         ! Indices on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                         ! Indices on the coarse path where GL
                                         ! corrections get applied.

     real(rp) :: singularity(1:size(inc_rad_path)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.
     logical :: do_calc(1:size(inc_rad_path)) ! Flags on coarse path where do_calc_c
                                         ! or (do_gl and any corresponding
                                         ! do_calc_f).

! Begin code

    ! d_delta_df is set to zero by the caller outside of all its loops.
    ! We keep track of where we create nonzeros, and replace them by zeros
    ! on the next call.  This is done because the vast majority of
    ! d_delta_df elements are zero.

    npc = size(indices_c)

    sps_n = ubound(Grids_f%l_z,1)

    do sps_i = 1, sps_n

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

        d_delta_df(nz_d_delta_df(:nnz_d_delta_df(sv_i),sv_i),sv_i) = 0.0
        nnz_d_delta_df(sv_i) = 0

! Skip the masked derivatives, according to the l2cf inputs

        if ( .not. Grids_f%deriv_flags(sv_i) ) then
          drad_df(sv_i) = 0.0
          cycle
        end if

! find where the non zeros are along the path

        call get_do_calc_indexed ( size(do_gl), do_calc_f(:,sv_i), indices_c, &
          & gl_inds, do_gl, do_calc, n_inds, nz_d_delta_df(:,sv_i) )

        nnz_d_delta_df(sv_i) = n_inds
        if ( n_inds > 0 ) then

          inds => nz_d_delta_df(1:n_inds,sv_i)

          no_to_gl = count(do_gl(inds))

          if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

            all_inds => all_inds_B(1:no_to_gl)
            more_inds => more_inds_B(1:no_to_gl)

            call get_inds ( do_gl, do_calc, more_inds, all_inds )
          end if

          if ( grids_f%lin_log(sps_i) ) then

            do i = 1, n_inds ! Don't trust the compiler to fuse loops
              ii = inds(i)
              iii = indices_c(ii)
              singularity(ii) = beta_path_c(ii,sps_i) &
                        & * eta_zxp_f(iii,sv_i) * sps_path(iii,sps_i)
              d_delta_df(ii,sv_i) = singularity(ii) * del_s(ii)
            end do ! i

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
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}.
      !  In the second integral, $G(\zeta)$ is {\tt beta\_path\_f * eta\_zxp\_f *
      !  sps\_path} -- which have already been evaluated at the appropriate
      !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}.  The weights
      !  are {\tt gw}.

            do i = 1, no_to_gl
              aa = all_inds(i)
              ga = gl_inds(aa)
              ii = more_inds(i)
              d_delta_df(ii,sv_i) = d_delta_df(ii,sv_i) + &
                & del_zeta(ii) * &
                & sum( (beta_path_f(aa:aa+ng-1,sps_i) &
                     &   * eta_zxp_f(ga:ga+ng-1,sv_i) * sps_path(ga:ga+ng-1,sps_i) - &
                     &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
            end do

            ! Refraction correction
            d_delta_df(inds,sv_i) = ref_cor(inds) * d_delta_df(inds,sv_i) * &
                                  & exp(-grids_f%values(sv_i))

          else

            do i = 1, n_inds
              ii = inds(i)
              singularity(ii) = beta_path_c(ii,sps_i) &
                        & * eta_zxp_f(indices_c(ii),sv_i)
              d_delta_df(ii,sv_i) = singularity(ii) * del_s(ii)
            end do ! i

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
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_df}.
      !  In the second integral, $G(\zeta)$ is {\tt beta\_path\_f *
      !  eta\_zxp\_f}~-- which have already been evaluated at the appropriate
      !  abscissae~-- and $G(\zeta_i)$ is {\tt singularity}.  The weights
      !  are {\tt gw}.

            do i = 1, no_to_gl
              aa = all_inds(i)
              ga = gl_inds(aa)
              ii = more_inds(i)
              d_delta_df(ii,sv_i) = d_delta_df(ii,sv_i) + &
                & del_zeta(ii) * &
                & sum( (beta_path_f(aa:aa+ng-1,sps_i) * eta_zxp_f(ga:ga+ng-1,sv_i) - &
                     &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
            end do

            ! Refraction correction
            d_delta_df(inds,sv_i) = ref_cor(inds) * d_delta_df(inds,sv_i)

          end if

          i_begin = max(i_start,min(inds(1),i_stop))

          call dscrt_dx ( tan_pt, d_delta_df(:,sv_i), inc_rad_path, &
                       &  i_begin, i_stop, drad_df(sv_i))

        else
          drad_df(sv_i) = 0.0
        end if

      end do ! sv_i

    end do ! sps_i

  end subroutine drad_tran_df

!--------------------------------------------------  drad_tran_dt  -----
! This is the radiative transfer derivative wrt temperature model

  subroutine DRad_tran_dt ( del_zeta, h_path_c, dh_dt_path_c, &
                         &  alpha_path_c, dAlpha_dT_path_c, &
                         &  eta_zxp_c, do_calc_t_c, &
                         &  do_calc_hyd_c, del_s, ref_cor, h_tan, dh_dt_tan, &
                         &  do_gl, gl_inds, h_path_f, t_path_f, dh_dt_path_f, &
                         &  alpha_path_f, dAlpha_dT_path_f, &
                         &  eta_zxp_f, do_calc_t_f, &
                         &  ds_dh, dh_dz_gw, ds_dz_gw, dt_scr_dt, &
                         &  tau, inc_rad_path, i_start, tan_pt, i_stop, &
                         &  deriv_flags, pfa_update, drad_dt )

    use GLNP, only: NG
    use MLSKinds, only: RP, IP
    use SCRT_DN_M, ONLY: DSCRT_DT, DSCRT_DX

! Inputs

    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: h_path_c(:)     ! path heights + req on main grid km.
    real(rp), intent(in) :: dh_dt_path_c(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on main grid.
    real(rp), intent(in) :: alpha_path_c(:) ! path absorption(km^-1)
!                                             on main grid.
    real(rp), intent(in) :: dAlpha_dT_path_c(:) ! path dAlpha/dT on main grid
    real(rp), intent(in) :: eta_zxp_c(:,:)  ! representation basis function
!                                              main grid.
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
    integer, intent(in) :: GL_Inds(:)       ! Where is do_gl true?
    real(rp), intent(in) :: h_path_f(:)     ! path heights + req on gl grid km.
    real(rp), intent(in) :: t_path_f(:)     ! path temperature(K) on gl grid.
    real(rp), intent(in) :: dh_dt_path_f(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on gl grid.
    real(rp), intent(in) :: alpha_path_f(:) ! path absorption(km^-1) on gl grid.
    real(rp), intent(in) :: dAlpha_dT_path_f(:) ! path dAlpha/dT on gl grid
    real(rp), intent(in) :: eta_zxp_f(:,:)  ! representation basis function
!                                             gl grid.
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

    integer(ip) :: A, AA, B, GA
    integer(ip) :: i, ii, i_begin, n_inds, n_path, no_to_gl, p_i, sv_i
    integer(ip), target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer(ip), target, dimension(1:size(inc_rad_path)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                         ! on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                         ! Indices on the coarse path where GL
                                         ! corrections get applied.

    real(rp) :: d_delta_dt(size(eta_zxp_c,1),size(eta_zxp_c,2)) ! path x sve.
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

    do sv_i = 1 , size(eta_zxp_c,dim=2)
      drad_dt(sv_i) = 0.0
      if ( .not. deriv_flags(sv_i)) cycle
      i_begin = i_start

! do the absorption part
! combine non zeros flags for both the main and gl parts

      call get_do_calc ( do_calc_t_c(:,sv_i), do_calc_t_f(:,sv_i), do_gl, &
        & do_calc, n_inds, inds_B )

      if ( n_inds > 0 ) then
        inds => inds_B(1:n_inds)

        i_begin = max(inds(1)-1, i_start)

        do i = 1, n_inds ! Don't trust the compiler to fuse loops
          ii = inds(i)
          singularity(ii) = dAlpha_dT_path_c(ii) * eta_zxp_c(ii,sv_i)
          d_delta_dt(ii,sv_i) = singularity(ii) * del_s(ii)
        end do ! i

! see if anything needs to be gl-d

        no_to_gl = count(do_gl(inds))
        if ( no_to_gl > 0 ) then

          all_inds => all_inds_B(1:no_to_gl)
          more_inds => more_inds_B(1:no_to_gl)

          call get_inds ( do_gl, do_calc, more_inds, all_inds )

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
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_dt}.
      !  In the second integral, $G(\zeta)$ is {\tt alphaxn\_path\_f *
      !  eta\_zxp\_f / t\_path\_f}~-- which have already been evaluated at
      !  the appropriate abscissae~-- and $G(\zeta_i)$ is {\tt
      !  singularity}.  The weights are {\tt gw}.

          do i = 1, no_to_gl
            aa = all_inds(i)
            ga = gl_inds(aa)
            ii = more_inds(i)
            d_delta_dt(ii,sv_i) = d_delta_dt(ii,sv_i) + &
              & del_zeta(ii) * &
              & sum( (dAlpha_dT_path_f(aa:aa+ng-1) * eta_zxp_f(aa:aa+ng-1,sv_i) - &
                   &  singularity(ii)) * &
                   & ds_dz_gw(ga:ga+ng-1) )
          end do
        end if ! no_to_gl > 0

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
            ! Don't test do_calc: There may be GL corrections even if
            ! dh_dt_path_c (from whence came do_calc) is zero.
            d_delta_dt(p_i,sv_i) = d_delta_dt(p_i,sv_i) + &
              & del_zeta(p_i) * &
              & sum( ( alpha_path_f(a:b-1) - alpha_path_c(p_i) ) *   &
              & (((2.0_rp*h_path_f(a:b-1)**2 - 3.0_rp*h_tan**2)      &     
              &   * dh_dt_path_f(a:b-1,sv_i) +                       &     
              &   h_path_f(a:b-1) * h_tan * dh_dt_tan(sv_i)) /       &     
              &  (sqrt(h_path_f(a:b-1)**2 - h_tan**2))**3            &     
              &  + eta_zxp_f(a:b-1,sv_i) * ds_dh(gl_inds(a:b-1)) /   &     
              &  t_path_f(a:b-1)) * dh_dz_gw(gl_inds(a:b-1)) )
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

        call dscrt_dt ( tan_pt, d_delta_dt(:,sv_i), tau, inc_rad_path, dt_scr_dt(:,sv_i), &
                      & i_begin, i_stop, drad_dt(sv_i) )

      end if

    end do ! sv_i

  end subroutine DRad_tran_dt

!--------------------------------------------------  drad_tran_dx  -----
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )

  subroutine DRad_tran_dx ( indices_c, gl_inds, del_zeta, Grids_f,        &
                         &  eta_zxp, sps_path, sps_map, do_calc_f,        &
                         &  dbeta_path_c, dbeta_path_f, do_gl, del_s,     &
                         &  ref_cor, ds_dz_gw, inc_rad_path, tan_pt,      &
                         &  i_stop, drad_dx )
    use GLNP, only: NG
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP, IP
    use SCRT_DN_M, ONLY: DSCRT_DX

! Inputs

    integer(ip), intent(in) :: indices_c(:)  ! coarse grid indicies
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indicies
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
    integer(ip), intent(in) :: i_stop        ! path stop index

! Outputs

    real(rp), intent(out) :: drad_dx(:)      ! derivative of radiances wrt x
!                                              state vector element. (K)
! Internals

    integer(ip) :: AA, GA, I, II, III
    integer(ip) :: i_start, n_inds, no_to_gl, sps_i, sps_m, sps_n, sv_i
    integer(ip), target, dimension(1:size(inc_rad_path)) :: all_inds_B
    integer(ip), target, dimension(1:size(inc_rad_path)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                         ! on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
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

        d_delta_dx = 0.0_rp

! find where the non zeros are along the path

        call get_do_calc_indexed ( size(do_gl), do_calc_f(:,sv_i), indices_c, &
          & gl_inds, do_gl, do_calc, n_inds, inds_B )

        if ( n_inds > 0 ) then

          inds => inds_B(1:n_inds)

          do i = 1, n_inds ! Don't trust the compiler to fuse loops
            ii = inds(i)
            iii = indices_c(ii)
            singularity(ii) = dbeta_path_c(ii,sps_i) &
                            &  * eta_zxp(iii,sv_i) * sps_path(iii,sps_m)
            d_delta_dx(ii) = singularity(ii) * del_s(ii)
          end do ! i

          no_to_gl = count(do_gl(inds))
          if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

            all_inds => all_inds_B(1:no_to_gl)
            more_inds => more_inds_B(1:no_to_gl)

            call get_inds ( do_gl, do_calc, more_inds, all_inds )

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
      !  $G(\zeta_i) (\zeta_{i-1}-\zeta_i)$.  Here, it is {\tt d\_delta\_dx}.
      !  In the second integral, $G(\zeta)$ is {\tt dbeta\_path\_f *
      !  eta\_zxp\_f\_f * sps\_path\_f}~-- which have already been evaluated at
      !  the appropriate abscissae~-- and $G(\zeta_i)$ is {\tt
      !  singularity}.  The weights are {\tt gw}.

            do i = 1, no_to_gl
              aa = all_inds(i)
              ga = gl_inds(aa)
              ii = more_inds(i)
              d_delta_dx(ii) = d_delta_dx(ii) + &
                & del_zeta(ii) * &
                & sum( (dbeta_path_f(aa:aa+ng-1,sps_i) &
                     &   * eta_zxp(ga:ga+ng-1,sv_i) * sps_path(ga:ga+ng-1,sps_m) - &
                     &  singularity(ii)) * ds_dz_gw(ga:ga+ng-1) )
            end do

          end if

          ! Refraction correction
          d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)

          i_start = min(inds(1),i_stop)

          call dscrt_dx ( tan_pt, d_delta_dx, inc_rad_path, i_start, i_stop, &
                       &  drad_dx(sv_i) )

        else
          drad_dx(sv_i) = 0.0
        end if

      end do

    end do

  end subroutine DRad_tran_dx

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

! =====     Private Procedures     =====================================

  ! ----------------------------------------  Get_Do_Calc_Indexed  -----
  subroutine Get_Do_Calc_Indexed ( N, Do_Calc_all, C_Inds, F_Inds, Do_GL, &
    & Do_Calc, N_Inds, Inds )

  ! Set Do_Calc if Do_Calc_All(c_inds) or Do_GL and any of the corresponding
  ! Do_Calc_All(f_inds) flags are set.

    use GLNP, ONLY: Ng
    use MLSKinds, only: IP

    integer, intent(in) :: N ! sizes on coarse grid

#if (defined NAG) || defined (IFC)
!     Assumed-shape arguments are slower than assumed size
    logical, intent(in) :: Do_Calc_all(*) ! On the entire path
    integer(ip), intent(in) :: C_Inds(*)  ! Indices in Do_Calc_All for coarse grid
    integer(ip), intent(in) :: F_Inds(*)  ! Indices in Do_Calc_All for find grid
    logical, intent(in) :: Do_GL(*)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(*)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(*)       ! Indices where do_calc is true
#elif defined LF95
!     Assumed-shape arguments are faster than assumed size
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer(ip), intent(in) :: C_Inds(:)  ! Indices in Do_Calc_All for coarse grid
    integer(ip), intent(in) :: F_Inds(:)  ! Indices in Do_Calc_All for find grid
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(:)       ! Indices where do_calc is true
#else
!     Assumed-shape arguments are faster than assumed size
    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer(ip), intent(in) :: C_Inds(:)  ! Indices in Do_Calc_All for coarse grid
    integer(ip), intent(in) :: F_Inds(:)  ! Indices in Do_Calc_All for fine grid
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.
    integer, intent(out) :: N_Inds        ! count(do_calc)
    integer, intent(out) :: Inds(:)       ! Indices where do_calc is true
#endif

    integer :: I, P_I
!     integer :: J
    integer :: K
!     logical :: T_calc(size(f_inds)/ng)

!     do_calc = do_calc_all(c_inds)
    i = 1 - Ng
    n_inds = 0
    do p_i = 1, n
      do_calc(p_i) = do_calc_all(c_inds(p_i))
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

    use GLNP, ONLY: Ng
    use MLSKinds, only: IP

    implicit NONE

  ! Inputs
    logical, intent(in) :: Do_GL(:)          ! path flag indicating where to do
      !                                        gl integrations.
    logical, intent(in) :: Do_Calc(:)

  ! Outputs
    integer(ip), intent(out) :: More_Inds(:)
    integer(ip), intent(out) :: All_Inds(:)

    integer :: I, J, K, L, P_I

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
