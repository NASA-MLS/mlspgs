! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Get_D_Deltau_Pol_M

  implicit NONE
  private
  public :: Get_D_Deltau_Pol_DF, Get_D_Deltau_Pol_DT

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

! ------------------------------------------  Get_D_Deltau_Pol_DF  -----
  subroutine Get_D_Deltau_Pol_DF ( CT, STCP, STSP, indices_c, del_zeta, Grids_f, &
               &  beta_path_pol, eta_zxp_f, do_calc_f, sps_path, Del_S, &
               &  incoptdepth, ref_cor, &
               &  d_delta_df, D_Deltau_Pol_DF  )

    use DExdT_m, only: dExdT
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSCommon, only: RP, IP
    use Opacity_m, only: Opacity
    use Where_M, only: Where

    ! SVE == # of state vector elements
    real(rp), intent(in) :: CT(:)           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), intent(in) :: STCP(:)         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), intent(in) :: STSP(:)         ! Sin(Theta) Sin(Phi)
    integer(ip), intent(in) :: indices_c(:) ! coarse grid indicies
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    complex(rp), intent(in) :: beta_path_pol(:,:,:) ! -1:1 x path x species.
!                                              cross section for each species
!                                              on coarse grid.
    real(rp), intent(in) :: eta_zxp_f(:,:)   ! fine path x sve
!                                              representation basis function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where
!                                              eta_zxp_f is not zero.
    real(rp), intent(in) :: sps_path(:,:)    ! fine path x species.
!                                              Path species function.
    real(rp), intent(in) :: Del_S(:)        ! unrefracted path length.  This
      !                                       is for the whole coarse path, not
      !                                       just the part up to the black-out
    complex(rp), intent(in) :: Incoptdepth(:,:,:) ! negative of incremental
      !                                       optical depth.  2 x 2 x path
    real(rp), intent(in) :: Ref_cor(:)      ! refracted to unrefracted path
    !                                         length ratios.
    real(rp), intent(in) :: d_delta_df(:,:)  ! derivative of delta wrt
!                                              mixing ratio state vector
!                                              element.

! Outputs

    complex(rp), intent(out) :: D_Deltau_Pol_DF(:,:,:,:) ! 2 x 2 x path x sve.
!                                              derivative of delta Tau wrt
!                                              mixing ratio state vector
!                                              element.

! Internals

    complex(rp) :: D_Delta_DF_Pol(-1:1,size(indices_c))
    complex(rp) :: D_Incoptdepth_df(2,2,size(indices_c))
    integer :: I_Stop                        ! Length of coarse path
    integer :: Inds(size(indices_c))         ! Where on the path to calc
    integer :: N_Inds                        ! Effective size of Inds
    integer :: N_Sps                         ! Number of species
    integer :: P_I                           ! Path index
    integer :: SPS_I                         ! Species index
    integer :: SV_I                          ! State vector index

    n_sps = ubound(Grids_f%l_z,1)
    i_stop = size(indices_c)

    do sps_i = 1, n_sps

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

! Skip the masked derivatives, according to the l2cf inputs

        if ( .not. Grids_f%deriv_flags(sv_i) ) cycle

        n_inds = count(do_calc_f(indices_c,sv_i))
        if ( n_inds == 0 ) cycle

        d_delta_df_pol = 0.0

        call where ( do_calc_f(indices_c,sv_i), inds(:n_inds) )

        if ( grids_f%lin_log(sps_i) ) then

          do p_i = 1, n_inds

            d_delta_df_pol(:,inds(p_i)) = beta_path_pol(:,inds(p_i),sps_i) &
                      & * sps_path(indices_c(inds(p_i)),sps_i) &
                      & / exp(grids_f%values(sv_i))

          end do ! p_i

        else

          d_delta_df_pol(:,inds(:n_inds)) = beta_path_pol(:,inds(:n_inds),sps_i)

        end if

        ! Finish the integration
        do p_i = 1, n_inds
          d_delta_df_pol(:,inds(p_i)) = d_delta_df_pol(:,inds(p_i)) * &
            & eta_zxp_f(indices_c(inds(p_i)),sv_i) * del_s(inds(p_i)) * &
            & ref_cor(inds(p_i))
        end do ! p_i

        ! Now add in contribution from scalar model, 0.25 for +/- sigma,
        ! 0.5 for pi.
        do p_i = 1, i_stop
          d_delta_df_pol(:,p_i) = d_delta_df_pol(:,p_i) + &
            & 0.25_rp * d_delta_df(p_i,sv_i)
          d_delta_df_pol(0,p_i) = d_delta_df_pol(0,p_i) + &
            & 0.25_rp * d_delta_df(p_i,sv_i)
        end do ! p_i
        ! d_delta_df_pol is now really \int incremental opacity ds.

        call opacity ( ct, stcp, stsp, d_delta_df_pol, d_incoptdepth_df )

        do p_i = 1, i_stop             ! along the path
          if ( eta_zxp_f(indices_c(p_i),sv_i) /= 0.0 &
            & .or. d_delta_df(p_i,sv_i) /= 0.0 &
            & .or. do_calc_f(indices_c(p_i),sv_i) ) then
            call dExdT ( incoptdepth(:,:,p_i), -d_incoptdepth_df(:,:,p_i), &
                       & d_deltau_pol_df(:,:,p_i,sv_i) ) ! d exp(incoptdepth) / df
          else
            d_deltau_pol_df(:,:,p_i,sv_i) = 0.0_rp
          end if
        end do ! p_i

      end do ! sv_i

    end do ! sps_i

  end subroutine Get_D_Deltau_Pol_DF

!{\newpage

! ------------------------------------------  Get_D_Deltau_Pol_DT  -----

!{Compute {\tt D\_Deltau\_Pol\_DT}.
!
! Assume $\beta$ can be approximated by $\hat\beta = \beta_0
! \left(\frac{T}{T_0}\right)^n$, but we don't know $n$.  So evaluate $\beta$
! for $T$, $T+$ some $\delta T$, and $T-$ some $\delta T$ (not necessarily the
! same $\delta T$).  Taking logarithms and differences, we get three estimates
! $n = \ln ( \beta(T+\delta T)/\beta(T-\delta T) ) /
!      \ln ( (T+\delta T) / (T-\delta T) )$,
! $n = \ln ( \beta(T+\delta T)/\beta(T) ) /
!      \ln ( (T+\delta T) / T)$, and
! $n = \ln ( \beta(T)/\beta(T-\delta T) ) /
!      \ln ( T / (T-\delta T))$.
! (The terms involving $\beta_0$ and $n \ln T_0$ cancel out.)
!
! Now observe that $\frac{\partial \hat\beta}{\partial T} = \frac{n}T \hat\beta$.
! Even better, use $\frac{n}T \beta$.  In any case, use a weighted average of
! the three values of $n$ computed above.  Also remember that
! $T = \sum_{i=1}^N \eta_i T_i$, where $N$ is the number of elements in the
! representation basis and $\eta_i$ are the coefficients, and what we really
! want is
! $\frac{\partial\beta}{\partial T_m} = 
!  \frac{n}T \beta \frac{\partial  T}{\partial T_m} =
!  \frac{n}T \beta \eta_m$, where $m$ is a state vector element index.
!
! Remembering that $\alpha = \sum_{i=1}^s f_i \beta_i$, and
! having $\frac{\partial \beta}{\partial T_m}$, compute
! $\frac{\partial \alpha}{\partial T_m} =
!  \sum_{i=1}^s f_i \frac{\partial\beta_i}{\partial T_m}$ ($s$ is the number
! of species).
!
! Then compute the derivative of incremental optical depth,
! $\frac{\partial \Delta \delta_{i \rightarrow i-1}^k}{\partial T_m}$ from
! $\frac{\partial \alpha}{\partial T_m}$ using the {\tt Opacity} routine.
!
! Finally, compute {\tt D\_Deltau\_Pol\_DT} = $\frac{\partial \bf E}{\partial
! T_m} = \frac\partial{\partial T} \exp ( -\Delta \delta_{i \rightarrow i-1}^k)$
! using $\Delta \delta_{i \rightarrow i-1}^k$,
! $\frac{\partial \Delta \delta_{i \rightarrow i-1}^k}{\partial T_m}$
! and the {\tt dExDt} routine.
! 

  subroutine Get_D_Deltau_Pol_DT ( Frq, H, CT, STCP, STSP, My_Catalog, &
                & Beta_group, GL_slabs_M, GL_slabs_P, &
                & T_Path_C, T_Path_M, T_Path_P, T_Path_F, &
                & Beta_Path, Beta_Path_F, SPS_Path, Alpha_Path_C, Alpha_Path_f, &
                & Eta_zxp, Eta_zxp_f, Del_S, Path_inds, GL_Inds, &
                & Del_Zeta, Do_Calc_T_c, Do_Calc_T_f, Do_GL, &
                & ds_dh, dh_dz_gw, ds_dz_gw, Incoptdepth, Ref_cor, &
                & H_path_c, H_path_f, dH_dt_path_c, dH_dt_path_f, H_tan, dH_dt_tan, &
                & Do_calc_hyd_c, Deriv_Flags, D_Delta_dT, D_Deltau_Pol_DT )

    use DExdT_m, only: dExdT
    use Get_Beta_Path_m, only: Get_Beta_Path_Polarized
    use Get_Species_Data_m, only: Beta_Group_T
    use GLNP, only: NG
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: R8, RP, IP
    use Opacity_m, only: Opacity
    use Physics, only: H_OVER_K
    use Rad_Tran_m, only: Get_Do_Calc
    use SpectroscopyCatalog_m, only: CATALOG_T

  ! Arguments
    ! SVE == # of state vector elements
    real(r8), intent(in) :: Frq             ! frequency in MHz
    real(rp), intent(in) :: H(:)            ! Magnetic field component in
                                            ! instrument polarization on the path
    real(rp), intent(in) :: CT(:)           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), intent(in) :: STCP(:)         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), intent(in) :: STSP(:)         ! Sin(Theta) Sin(Phi)
    type(catalog_t), intent(in) :: My_Catalog(:)
    type (beta_group_T), intent(in) :: Beta_group(:)
    type (slabs_struct), intent(in) :: GL_slabs_m(:,:), GL_slabs_p(:,:) ! for T -/+ del_T
    real(rp), intent(in) :: T_Path_c(:)     ! path temperatures on coarse grid
    real(rp), intent(in) :: T_Path_m(:), T_Path_p(:) ! path temperatures -/+ del_temp
    !                                         on fine grid -- index with path_inds
    real(rp), intent(in) :: T_Path_f(:)     ! path temperatures on GL grid
    complex(rp), intent(in) :: Beta_Path(-1:,:,:) ! beta * tanh(h nu / 2 k t) on
                                            ! coarse path.  -1:1 x path x sps
    complex(rp), intent(in) :: Beta_Path_f(-1:,:,:) ! beta * tanh(h nu / 2 k t) on
                                            ! GL path.  -1:1 x path x sps
    real(rp), intent(in) :: SPS_Path(:,:)   ! species on whole path, path x sps
    complex(rp), intent(in) :: Alpha_Path_c(-1:,:) ! -1:1 x path on coarse grid
    complex(rp), intent(in) :: Alpha_Path_f(-1:,:) ! -1:1 x path on fine grid
    real(rp), intent(in) :: Eta_zxp(:,:)    ! representation basis function
      !                                       on coarse grid.  path x sve
    real(rp), intent(in) :: Eta_zxp_f(:,:)  ! representation basis function
      !                                       on fine grid.  path x sve
    real(rp), intent(in) :: Del_S(:)        ! unrefracted path length.  This
      !                                       is for the whole coarse path, not
      !                                       just the part up to the black-out
    integer(ip), intent(in) :: Path_inds(:) ! indicies for reading coarse path
      !                                       elements from gl_slabs and sps_path
    integer(ip), intent(in) :: GL_Inds(:)   ! indices for reading fine path
      !                                       elements from sps_path
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    logical, intent(in) :: do_calc_t_c(:,:) ! Indicates where the
      !              representation basis function is not zero on main grid.
    logical, intent(in) :: do_calc_t_f(:,:) ! Indicates where the
      !              representation basis function is not zero on gl grid.
    logical, intent(in) :: do_gl(:)         ! Indicates where on the coarse path
      !                                       to do gl integrations.
    real(rp), intent(in) :: ds_dh(:)        ! path length wrt height derivative
      !                                       on entire grid.  Only the
      !                                       gl_inds part is used.
    real(rp), intent(in) :: dh_dz_gw(:)     ! path height wrt zeta derivative * gw.
      !                                       on entire grid.  Only the
      !                                       gl_inds part is used.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative * gw.
      !                                       on entire grid.  Only the
      !                                       gl_inds part is used.
    complex(rp), intent(in) :: Incoptdepth(:,:,:) ! negative of incremental
      !                                       optical depth.  2 x 2 x path
    real(rp), intent(in) :: Ref_cor(:)      ! refracted to unrefracted path
    !                                         length ratios.

    ! For hydrostatic
    real(rp), intent(in) :: H_path_c(:)     ! path heights + req on the main
      !              grid km. This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: H_path_f(:)     ! path heights + req on find grid
    real(rp), intent(in) :: dH_dt_path_c(:,:) ! derivative of coarse path height
    !                                         wrt temperature(km/K) on main grid.
    real(rp), intent(in) :: dH_dt_path_f(:,:) ! derivative of fine path height
    !                                         wrt temperature(km/K) on main grid.
    real(rp), intent(in) :: H_tan           ! tangent height + req (km).
    real(rp), intent(in) :: dH_dt_tan(:)    ! derivative of path height wrt
    !                                         temperature at the tangent (km/K).
    logical, intent(in) :: Do_calc_hyd_c(:,:) ! Where is the dh_dt function not
    !                                         zero on main grid?
    logical, intent(in) :: deriv_flags(:)   ! Indicates which temperature
!                                             derivatives to do

    ! From nonpolarized model
    real(rp), intent(in) :: D_Delta_DT(:,:) ! Incremental opacity derivatives
    !                                         schlep from drad_tran_dt.  Path x SVE

! Outputs

    complex(rp), intent(out) :: D_Deltau_Pol_DT(:,:,:,:) ! 2 x 2 x path x sve.
!                                              derivative of delta Tau wrt
!                                              temperature state vector
!                                              element. (K)

  ! Local variables
    integer :: A, B                  ! Indices for GL points
    complex(rp) :: Alpha_Path_N(-1:1)       ! alpha_path_n * N at one path point
    complex(rp) :: Alpha_Path_N_F(-1:1,NG)  ! alpha_path_n * N at one set of GL points
    complex(rp) :: Alpha_Path_N_T(-1:1,size(path_inds)) ! Alpha * n/T on the
                                     ! coarse path
    complex(rp) :: Alpha_Path_N_T_F(-1:1,size(t_path_f)) ! Alpha * n/T on the
                                     ! fine path
    complex(rp) :: Beta_0(-1:1), Beta_M(-1:1), Beta_P(-1:1) ! Single elements of
      ! Beta_Path, Beta_Path_M, Beta_Path_P multiplied by Tanh_Path,  Tanh_M, Tanh_P
    complex(rp), dimension(-1:1,size(path_inds),size(beta_group)) :: &
      & Beta_Path_M, &  ! At T_path_M on coarse path
      & Beta_Path_P     ! At T_path_P on coarse path
    complex(rp):: D_Alpha_DT_eta(-1:1,size(path_inds)) ! Singularity * Del_S
    complex(rp) :: D_Incoptdepth_dT(2,2,size(path_inds))
    logical :: Do_calc(1:size(path_inds)) ! do_calc_t_c .or. ( do_gl .and. any
                                     ! of the corresponding do_calc_t_f ).
    real(rp) :: F(ng)                ! Factor in GL that doesn't depend on
                                     ! sigma +/- or pi.
    real(rp) :: Fa, Fb               ! Hydrostatic integrand at ends of
                                     ! path segment
    real(r8) :: FrqHK                ! 0.5 * Frq * H_Over_K
    integer :: H_Stop                ! Stop point for hydrostatic parts
    integer :: I_start               ! Start point, not necessarily 1.
    integer :: I_stop                ! Stop point, which may be before N_Path
    integer :: Inds(size(path_inds)) ! Where on coarse path is DO_CALC and DO_GL?
    integer :: J, K, L
    real(rp) :: L_TTM, L_TPTM, L_TPT ! Logarithms of temperature ratios
    integer :: Mid                   ! tangent index along the path = N_Path/2
    complex(rp) :: N(-1:1)           ! Exponent of (T/T_0) in
               ! approximation to beta.  One each for Sigma_-, Pi and Sigma_+.
    integer :: N_Inds                ! How much of Inds is used?
    integer :: N_Path                ! Total coarse path length.
    integer :: N_Sps                 ! Number of species
    logical :: NeedFA                ! Need FA in hydrostatic calculation
    integer :: P_i                   ! Index on the path
    complex(rp) :: R0M(-1:1), RPM(-1:1), RP0(-1:1)  ! Beta ratios
    real(rp) :: S_Del_S              ! Sum of Del_S
    complex(rp) :: Singularity(-1:1,size(path_inds)) ! n/T * Alpha * Eta on the
                                     ! coarse path
    integer :: SV_i                  ! Index of state vector element
    real(rp) :: Tanh_M, Tanh_P       ! for T -/+ del_T

    frqhk = 0.5_r8 * Frq * H_Over_K
    i_stop = size(path_inds)
    n_path = size(del_zeta)
    mid = n_path / 2
    n_sps = size(sps_path,2)

    call get_beta_path_polarized ( frq, h, my_Catalog, beta_group, gl_slabs_m, &
      & path_inds, beta_path_m )
    call get_beta_path_polarized ( frq, h, my_Catalog, beta_group, gl_slabs_p, &
      & path_inds, beta_path_p )

    a = 1
    b = 1 + ng
    n_inds = 0
    do p_i = 1, i_stop
      alpha_path_n = (0.0_rp,0.0_rp)
      alpha_path_n_f = (0.0_rp,0.0_rp)
      k = path_inds(p_i) ! K is now index for coarse path

      l_ttm = log(t_path_c(p_i)/t_path_m(k))
      l_tptm = log(t_path_p(k)/t_path_m(k))
      l_tpt = log(t_path_p(k)/t_path_c(p_i))

      tanh_m = tanh( frqhk / t_path_m(k) )
      tanh_p = tanh( frqhk / t_path_p(k) )

      do j = 1, n_sps
        ! Solve for n
        beta_0 = beta_path(:,p_i,j) ! * tanh1(p_i) done by caller
        beta_m = beta_path_m(:,p_i,j) * tanh_m
        beta_p = beta_path_p(:,p_i,j) * tanh_p
        where ( beta_m /= 0.0 .and. beta_p /= 0.0 )
          rpm = log(beta_p/beta_m) / l_tptm
          where ( beta_0 /= 0.0 )
            r0m = log(beta_0/beta_m) / l_ttm
            rp0 = log(beta_p/beta_0) / l_tpt
            n = 0.25 * ( r0m + 2.0 * rpm + rp0 )
          elsewhere
            n = rpm
          end where
        elsewhere ( beta_m /= 0.0 .and. beta_0 /= 0.0 )
          n = log(beta_0/beta_m) / l_ttm
        elsewhere ( beta_0 /= 0.0 .and. beta_p /= 0.0 )
          n = log(beta_p/beta_0) / l_tpt
        elsewhere
          n = 0.0
        end where

        ! not quite D alpha, because we haven't divided by T.  beta_path was
        ! multiplied by tanh in FullForwardModel.
        alpha_path_n = alpha_path_n + n * beta_0 * sps_path(k,j)
        ! Use the same N for the GL points.  Getting 3*NG more N's would cost
        ! 6*NG more Beta's.  N varies slowly, so this is probably OK.
        ! alpha_path_n_f is dimensioned (-1:1,NG).
        if ( do_gl(p_i) ) then
          do l = -1, 1
            alpha_path_n_f(l,:) = alpha_path_n_f(l,:) + n(l) * &
              & beta_path_f(l,a:b-1,j) * sps_path(gl_inds(a:b-1),j)
          end do ! l
        end if
      end do ! j = 1, n_sps

      alpha_path_n_t(:,p_i) = alpha_path_n / t_path_c(p_i)
      if ( do_gl(p_i) ) then
        do l = -1, 1
          alpha_path_n_t_f(l,a:b-1) = alpha_path_n_f(l,:) / t_path_f(a:b-1)
        end do
        a = b
        b = b + ng
      end if
    end do ! p_i = 1, i_stop

    do sv_i = 1, size(eta_zxp,2) ! state vector elements
      if ( .not. deriv_flags(sv_i)) then
        d_deltau_pol_dT(:,:,:,sv_i) = 0.0_rp
        cycle
      end if
      i_start = 1

! do the absorption part
! combine non zeros flags for both the main and gl parts

      call get_do_calc ( do_calc_t_c(:,sv_i), do_calc_t_f(:,sv_i), do_gl, &
        & do_calc )

      a = 1
      do p_i = 1, i_stop
        if ( do_calc(p_i) ) then
          !{ d\_alpha\_dT\_eta$_i = 
          !    \frac{\partial \alpha}{\partial T_i} \text{d}s=
          !    \frac{\partial \alpha}{\partial T}
          !    \frac{\partial T}{\partial T_i} \text{d}s =
          !    \frac{\partial \alpha}{\partial T} \eta_i \text{d}s$
          singularity(:,p_i) = alpha_path_n_t(:,p_i) * eta_zxp(p_i,sv_i)
          d_alpha_dT_eta(:,p_i) = singularity(:,p_i) * del_s(p_i)
        else
          singularity(:,p_i) = 0.0_rp
          d_alpha_dT_eta(:,p_i) = 0.0_rp
        end if
        ! Do GL if needed here
        if ( do_gl(p_i) ) then
          b = a + ng
          if ( do_calc(p_i) ) then
            f = ds_dz_gw(gl_inds(a:b-1))
            do l = -1, 1
              d_alpha_dT_eta(l,p_i) = d_alpha_dT_eta(l,p_i) + &
                 & del_zeta(p_i) * &
                 & sum( ( alpha_path_n_t_f(l,a:b-1) * eta_zxp_f(a:b-1,sv_i) - &
                 &        singularity(l,p_i) ) * f )
            end do ! l
          end if
          a = b
        end if
      end do ! p_i

      ! Now do the hydrostatic part
      !??? This only goes as far as I_Stop.  We may     ???
      !??? want to do it this way in DRad_Tran_dT also. ???

      ! First combine boundary flags
      do_calc = do_calc_hyd_c(:,sv_i)
      if ( i_stop < mid ) then
        do_calc(2:i_stop) =   do_calc(2:i_stop) .or. do_calc(1:i_stop-1)
        h_stop = i_stop
      else
        do_calc(2:i_stop) =   do_calc(2:i_stop) .or. do_calc(1:i_stop-1) .or. do_calc(mid)
        h_stop = mid - 1
      end if
      do_calc(1) = .false.
      s_del_s = sum(del_s(2:mid)) ! Yes, this goes to the midpoint of the coarse path
      needFA = .true.
      fa = 0.0_rp ! In case n_path <= 4
      do p_i = 2 , h_stop
        if ( do_calc(p_i) ) then
          if ( needFA ) then
            fa = (h_path_c(p_i-1) * dh_dt_path_c(p_i-1,sv_i) &
            &   - h_tan * dh_dt_tan(sv_i)) / s_del_s
            needFA = .false.
          end if
          s_del_s = s_del_s - del_s(p_i)
          fb = (h_path_c(p_i) * dh_dt_path_c(p_i,sv_i) &
            & - h_tan * dh_dt_tan(sv_i)) / s_del_s
          d_alpha_dT_eta(:,p_i) = d_alpha_dT_eta(:,p_i) + alpha_path_c(:,p_i) * (fa - fb)
          fa = fb
        else
          s_del_s = s_del_s - del_s(p_i)
        end if
      end do ! p_i

      ! special processing at tangent.  fb is zero

      if ( i_stop >= mid ) then
        if ( do_calc(mid) ) &
          & d_alpha_dT_eta(:,mid) = d_alpha_dT_eta(:,mid) + alpha_path_c(:,mid) * fa
      end if
      if ( i_stop > mid + 1 ) then ! mid+1 instead of mid so that mid+2 will be
                                   ! in bounds if i_stop == 2.
        if ( do_calc(mid+1) ) then
          fa = (h_path_c(mid+2) * dh_dt_path_c(mid+2,sv_i) &
            & - h_tan * dh_dt_tan(sv_i)) / del_s(mid+1)
          d_alpha_dT_eta(:,mid+1) = d_alpha_dT_eta(:,mid+1) + alpha_path_c(:,mid+1) * fa
        else
          needFA = .true.
        end if

        do_calc(mid+2:i_stop-1) = do_calc(mid+2:i_stop-1) .or. do_calc(mid+3:i_stop) .or. do_calc(mid+1)
        if (i_stop == 2*mid) then
          h_stop = i_stop - 1
          do_calc(i_stop) = .false.
        else
          h_stop = i_stop
        end if

        s_del_s = del_s(mid+1)
        do p_i = mid + 2, h_stop
          if ( do_calc(p_i) ) then
            if ( needFA ) then
              fa = (h_path_c(p_i) * dh_dt_path_c(p_i,sv_i) &
                & - h_tan * dh_dt_tan(sv_i)) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s + del_s(p_i)
            fb = (h_path_c(p_i+1) * dh_dt_path_c(p_i+1,sv_i) &
              & - h_tan * dh_dt_tan(sv_i)) / s_del_s
            d_alpha_dT_eta(:,p_i) = d_alpha_dT_eta(:,p_i) + alpha_path_c(:,p_i)*(fb - fa)
            fa = fb
          else
            s_del_s = s_del_s + del_s(p_i)
          end if
        end do ! p_i
      end if

      ! Do GL for hydrostatic for any panels that need it; the singularity
      ! correction is alpha_path_c.
      ! Apply refraction correction.
      ! Add in contribution from scalar model, 0.25 for +/- sigma,
      ! 0.5 for pi.
      a = 1
      do p_i = 1, i_stop             ! along the path
        if ( do_gl(p_i) ) then
          b = a + ng
          if ( do_calc(p_i) ) then
            f = (((2.0_rp*h_path_f(a:b-1)**2 - 3.0_rp*dh_dt_tan(sv_i)**2) &
              &   * dh_dt_path_f(a:b-1,sv_i) +                            &
              &   h_path_f(a:b-1) * dh_dt_tan(sv_i) * dh_dt_tan(sv_i)) /  &
              &  (sqrt(h_path_f(a:b-1)**2 - dh_dt_tan(sv_i)**2))**3       &
              &  + eta_zxp_f(a:b-1,sv_i) * ds_dh(gl_inds(a:b-1)) /        &
              &  t_path_f(a:b-1)) * dh_dz_gw(gl_inds(a:b-1))
            do l = -1, 1
              d_alpha_dT_eta(l,p_i) = d_alpha_dT_eta(l,p_i) + &
                 & del_zeta(p_i) * &
                 & sum( ( alpha_path_f(l,a:b-1) - alpha_path_c(l,p_i) ) * f )
            end do ! l
          end if
          a = b
        end if
        d_alpha_dT_eta(:,p_i) = d_alpha_dT_eta(:,p_i) * ref_cor(p_i) + &
          & 0.25_rp * d_delta_dT(p_i,sv_i)
        d_alpha_dT_eta(0,p_i) = d_alpha_dT_eta(0,p_i) + &
          & 0.25_rp * d_delta_dT(p_i,sv_i)
      end do ! p_i
      ! d_alpha_dT_eta is now really \int incremental opacity ds.

      call opacity ( ct, stcp, stsp, d_alpha_dT_eta, d_incoptdepth_dT )

      do p_i = 1, i_stop             ! along the path
        if ( eta_zxp(p_i,sv_i) /= 0.0 .or. d_delta_dT(p_i,sv_i) /= 0.0 .or. &
          & do_calc(p_i) ) then
          call dExdT ( incoptdepth(:,:,p_i), -d_incoptdepth_dT(:,:,p_i), &
                     & d_deltau_pol_dT(:,:,p_i,sv_i) ) ! d exp(incoptdepth) / dT
        else
          d_deltau_pol_dT(:,:,p_i,sv_i) = 0.0_rp
        end if
      end do ! p_i
    end do ! sv_i

  end subroutine Get_D_Deltau_Pol_DT

!-----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Get_D_Deltau_Pol_M

! $Log$
! Revision 2.19  2003/11/24 22:08:30  vsnyder
! Remove an unnecessary variable
!
! Revision 2.18  2003/11/04 02:01:19  vsnyder
! Add 'FA = 0.0' in case n_path <= 4
!
! Revision 2.17  2003/11/04 01:56:19  vsnyder
! Cosmetic changes
!
! Revision 2.16  2003/11/01 03:03:46  vsnyder
! Use ds_dz_gw instead of ds_dh, dh_dz and gw; use del_zeta from FullForwardModel
!
! Revision 2.15  2003/10/30 20:45:51  vsnyder
! Finish removing z_path_c in favor of del_zeta
!
! Revision 2.14  2003/10/30 20:38:16  vsnyder
! Get del_zeta from FullForwardModel.  Code for GL for derivatives.
! Still something wrong -- zero-field run doesn't agree with nonpolarized model
!
! Revision 2.13  2003/09/10 22:35:09  vsnyder
! Avoid a subscript out-of-bounds in case of path length == 2
!
! Revision 2.12  2003/09/09 00:03:48  vsnyder
! Repair an indexing blunder
!
! Revision 2.11  2003/08/15 20:29:26  vsnyder
! Implement polarized VMR derivatives
!
! Revision 2.10  2003/08/15 18:50:22  vsnyder
! Preparing the way for polarized vmr derivatives
!
! Revision 2.9  2003/06/27 22:04:50  vsnyder
! Simplify calculation of N
!
! Revision 2.8  2003/06/13 23:54:41  vsnyder
! Include nonpolarized d_delta_dT even if there's no polarized stuff
!
! Revision 2.6  2003/06/13 00:00:09  vsnyder
! Move multiplication of beta_path by tanh into FullForwardModel
!
! Revision 2.5  2003/06/10 15:07:08  bill
! fixed polarized t-derivs
!
! Revision 2.4  2003/06/09 20:52:37  vsnyder
! More work on polarized derivatives
!
! Revision 2.3  2003/05/24 02:26:18  vsnyder
! More work on polarized temperature derivatives
!
! Revision 2.2  2003/05/15 20:50:34  vsnyder
! Correct some subscript errors -- coarse vs. fine path
!
! Revision 2.1  2003/05/15 03:27:56  vsnyder
! Initial commit
!
