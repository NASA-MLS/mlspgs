! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_d_Deltau_Pol_M

  implicit NONE
  private
  public :: Get_d_Deltau_Pol_dF, Get_d_Deltau_Pol_dT
    integer, parameter :: max_num_debugs = 10
    integer :: num_debugs

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains

! ------------------------------------------  Get_d_Deltau_Pol_DF  -----
  subroutine Get_d_Deltau_Pol_dF ( CT, STCP, STSP, Grids_f, Tan_pt_c, &
               &  Beta_Path_Pol, Tanh1_c, Eta_zxp, Sps_Path, Del_S,   &
               &  IncOptDepth, Ref_Cor, d_Delta_df, d_Deltau_Pol_dF  )

    use DExdT_m, only: dExdT
    use Dump_0, only: Dump
    use Get_do_Calc_Indexed_m, only: Get_do_Calc_Indexed_Coarse
    use GLNP, only: NG, NGP1
    use Load_Sps_Data_m, Only: Grids_t
    use MLSKinds, only: RP
    use Opacity_m, only: Opacity
    use Sparse_m, only: Sparse_t

    ! SVE == # of state vector elements
    real(rp), intent(in) :: CT(:)             ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), intent(in) :: STCP(:)           ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), intent(in) :: STSP(:)           ! Sin(Theta) Sin(Phi)
    type (Grids_T), intent(in) :: Grids_f     ! All the coordinates
    integer, intent(in) :: Tan_pt_c           ! Tangent point in coarse path
    complex(rp), intent(in) :: Beta_Path_Pol(:,:,:) ! -1:1 x path x species.
                                              ! cross section for each species
                                              ! on coarse grid.
    real(rp), intent(in) :: Tanh1_c(:)        ! tanh(h nu / k T) on coarse path
    class(sparse_t), intent(in) :: Eta_ZXP(:) ! Interpolating coefficients
                                              ! from state vector coordinates to
                                              ! combined coarse & fine path for
                                              ! each sps.
    real(rp), intent(in) :: Sps_Path(:,:)     ! fine path x species.
                                              ! Path species function.
    real(rp), intent(in) :: Del_S(:)          ! unrefracted path length.  This
                                              ! is for the whole coarse path, not
                                              ! just the part up to the black-out
    complex(rp), intent(in) :: IncOptDepth(:,:,:) ! negative of incremental
                                              ! optical depth.  2 x 2 x path
    real(rp), intent(in) :: Ref_Cor(:)        ! refracted to unrefracted path
                                              ! length ratios, only up to the
                                              ! black-out point, on coarse path.
    class (sparse_t), intent(in) :: d_Delta_df(:) ! derivative of delta wrt
                                              ! mixing ratio state vector
                                              ! element.  One per SPS.

! Outputs

    complex(rp), intent(out) :: d_Deltau_Pol_dF(:,:,:,:) ! 2 x 2 x coarse path x sve.
                                              ! derivative of delta Tau wrt
                                              ! mixing ratio state vector
                                              ! element.

! Internals

    integer :: B                             ! Which boundary of the layer, (1
                                             ! or 2 for near, far from tangent)
                                             ! or other boundary of the layer,
                                             ! on fine path
    integer :: BC                            ! B on coarse path when B is a
                                             ! boundary, not a boundary index
    complex(rp) :: Beta(-1:1,size(ref_cor))  ! Either Beta or Beta*f*exp(f)
                                             ! on a boundary
    complex(rp) :: d_Delta_df_Pol(-1:1,size(ref_cor)) ! Layer integral
    complex(rp) :: d_Incoptdepth_df(2,2,size(ref_cor))
!! The (:) in the next declaration is not required by the standard, but
!! ifort 17 produces an error complaining that eta_zxp needs to be an array
!  if it's omitted.
    logical :: Do_Calc_f(maxval(eta_zxp(:)%nRows))    ! Where Eta_ZXP_Col /= 0
    real(rp) :: Eta_ZXP_Col(maxval(eta_zxp(:)%nRows))
    integer :: I_Stop                        ! Length of coarse path
    integer :: IC                            ! inds(p_i) converted to coarse path
    integer :: II                            ! inds(p_i)
    integer :: Inds(size(ref_cor),2)         ! Where on the path to calc
    integer :: N_Inds                        ! Effective size of Inds
    integer :: P_I                           ! Path index
    integer :: Sparse_Col
    integer :: SPS_I                         ! Species index
    integer :: SV_I                          ! State vector index

    do_calc_f = .false.
    eta_zxp_col = 0
    i_stop = size(ref_cor)
    num_debugs = 0

! initialize entire array to zero
    d_deltau_pol_df = 0.0_rp

    do sps_i = 1, ubound(grids_f%l_v,1)

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)
        sparse_col = sv_i - Grids_f%l_v(sps_i-1)

! Skip the masked derivatives, according to the l2cf inputs

        if ( .not. Grids_f%deriv_flags(sv_i) ) cycle

        if ( eta_zxp(sps_i)%cols(sparse_col) == 0 ) cycle ! Zero column
        call eta_zxp(sps_i)%get_col ( sparse_col, eta_zxp_col, do_calc_f )

        ! This call is made only for polarized forward models that
        ! do vmr derivatives
        ! The indices it returns cause bounds errors.
        ! Is it the fault of the input (e.g., do_calc_f)
        ! or some coding error in the subroutine itself?
        ! Get indices for layers for which the interpolating coefficient
        ! is nonzero at either boundary.
        !????? Maybe could use Get_Do_Calc_Sparse, but it hasn't been tested
        call Get_do_calc_indexed_coarse ( i_stop, tan_pt_c, do_calc_f, &
          & n_inds, inds )
        if ( any(inds(1:n_inds,:) < 1 ) ) then
          if ( num_debugs < max_num_debugs ) then
            num_debugs = num_debugs + 1
            print *, 'Warning--some of inds are < 1'
            call Dump ( inds(1:n_inds,:), 'inds' )
            ! We note that the bounds on sps_path and eta_zxp_col
            ! are the same as do_calc_f, namely 533
            ! Based on their use below, this seems to be the intended
            ! range for the values of inds.
            ! On the other hand, the bounds of beta, beta_path_pol and 
            ! d_delta_df_pol are 134, so be careful not to use inds for them.
            print *, 'shape(do_calc_f): ', shape(do_calc_f)
            call Dump ( do_calc_f, 'do_calc_f', options='.' )
            print *, 'beta_path_pol: ', shape(beta_path_pol)
            print *, 'sps_path: ', shape(sps_path)
            print *, 'eta_zxp_col: ', shape(eta_zxp_col)
            print *, 'beta: ', shape(beta)
            print *, 'd_delta_df_pol: ', shape(d_delta_df_pol)

            call Get_do_calc_indexed_coarse ( i_stop, tan_pt_c, do_calc_f, &
            & n_inds, inds, debug=.true. )
          endif
        endif

        if ( n_inds == 0 ) cycle

        d_delta_df_pol = 0.0

        !{ Get $\frac{\partial \alpha}{\partial f^k} = \beta^k$ from
        !  $\alpha = \sum_{k} f^k \beta^k$.  Doesn't cater for the case of
        !  $\beta^k} depending upon $f$.

        if ( grids_f%lin_log(sps_i) ) then

          do p_i = 1, n_inds

            if ( any(inds(p_i,:) < 1 ) ) then
              if ( num_debugs < max_num_debugs ) then
                num_debugs = num_debugs + 1
                print *, 'p_i, inds(p_i,:) ', p_i, inds(p_i,:)
                cycle
              endif
            endif
            do b = 1, 2 ! b is boundary index
              ii = inds(p_i,b)
!             ic was too big by 1 (see remarks below)
              ! ic = 1 + ( ii + ng ) / ngp1 ! On coarse path
              ic = ( ii + ng ) / ngp1 ! On coarse path
              beta(:,ic) = beta_path_pol(:,ic,sps_i) &
                         & * sps_path(ii,sps_i) &
                         & * exp(-grids_f%values(sv_i))

            end do ! b
          end do ! p_i

        else

          do p_i = 1, n_inds
            ! print *, 'inds(p_i,1:2) ', inds(p_i,1:2)
            ! print *, '(inds(p_i,1:2) + ng)/ngp1 ', (inds(p_i,1:2) + ng)/ngp1
!             the 2nd index of beta was too big by 1 (see remarks below)
            beta(:, ( inds(p_i,1:2) + ng ) / ngp1) = &
              &  beta_path_pol(:, ( inds(p_i,1:2) + ng ) / ngp1, sps_i)
!             beta(:, 1 + ( inds(p_i,1:2) + ng ) / ngp1) = &
!               &  beta_path_pol(:, 1 + ( inds(p_i,1:2) + ng ) / ngp1, sps_i)
          enddo
! The NAG compiler complained that
!     Left-hand side of assignment has vector subscript
!     [1+(INDS(:N_INDS,1:2)+NG)/NGP1] with duplicate value 18
! so we replaced the following vector assugnment with the above do loop
!          beta(:,[1 + ( inds(:n_inds,1:2) + ng ) / ngp1]) = &
!            & beta_path_pol(:,[1 + ( inds(:n_inds,1:2) + ng ) / ngp1],sps_i)

        end if

        ! Finish the integration.  Include the factor of tanh(h nu / k T)
        ! that was not included in the beta computation.
!??? Should we (can we) also do GL here?  (We don't have Beta on the fine path)
        do p_i = 1, n_inds
          ii = inds(p_i,1)
          ic = 1 + ( ii + ng ) / ngp1 ! On coarse path
          b = inds(p_i,2)
          bc = 1 + ( b + ng ) / ngp1 ! On coarse path
          if ( b < 1 ) then
            print *, 'n_inds, p_i, ii, b ', n_inds, p_i, ii, b
            cycle
          endif
          ! It looks like ic and bc as defined above are too big by one
          ! Instead of ranging from  1 .. 134
          ! they range from 2 .. 135; so
          ic = ic - 1
          bc = bc - 1
! See, here is where we said to be careful not to use inds for d_delta_df_pol
!          d_delta_df_pol(:,ii) = &
          d_delta_df_pol(:,ic) = &
            & ( beta(:,ic) * ( eta_zxp_col(ii) * tanh1_c(ic) ) + &
            &   beta(:,bc) * ( eta_zxp_col(b) * tanh1_c(bc) ) ) * &
            & 0.5 * del_s(ic) * ref_cor(ic)
        end do ! p_i

        ! Now add in contribution from scalar model, 0.25 for sigma +/-,
        ! 0.5 for pi.  We only need the nonzeros from d_delta_df.
        ii = d_delta_df(sps_i)%cols(sparse_col) ! Last element in the column
        if ( ii /= 0 ) then
          do
            ii = d_delta_df(sps_i)%e(ii)%nc ! Element in next row in the column
            p_i = d_delta_df(sps_i)%e(ii)%r ! Row number of element
            d_delta_df_pol(:,p_i) = d_delta_df_pol(:,p_i) + &
              & 0.25_rp * d_delta_df(sps_i)%e(ii)%v
            d_delta_df_pol(0,p_i) = d_delta_df_pol(0,p_i) + &
              & 0.25_rp * d_delta_df(sps_i)%e(ii)%v
            if ( ii == d_delta_df(sps_i)%cols(sparse_col) ) exit
          end do ! ii
        end if
        ! d_delta_df_pol is now \int (d incremental opacity / df) ds.

        call opacity ( ct, stcp, stsp, d_delta_df_pol, d_incoptdepth_df )

        do p_i = 1, i_stop - 1 ! along the path, P_i is a layer index bounded
                               ! by coarse points P_i and P_i + 1

! I think these are mostly redundant to each other so WGR has
! replaced them with only one (left original uncommented in case I
! am wrong)
!          if ( eta_zxp_col(p_i*ngp1-ng,sv_i) /= 0.0 &
!            & .or. d_delta_df(p_i,sv_i) /= 0.0 &
!            & .or. do_calc_f(p_i*ngp1-ng) ) then
!           if ( do_calc_f(p_i*ngp1-ng) ) then
          ! Compare to D_Delta_Pol_dT, which looks at the input to Opacity:
          if ( any ( d_delta_df_pol(:,p_i) /= 0 ) ) then
            call dExdT ( incoptdepth(:,:,p_i), -d_incoptdepth_df(:,:,p_i), &
                       & d_deltau_pol_df(:,:,p_i,sv_i) ) ! d exp(incoptdepth) / df
          else
            d_deltau_pol_df(:,:,p_i,sv_i) = 0.0_rp
          end if

        end do ! p_i

        call eta_zxp(sps_i)%clear_col ( sparse_col, eta_zxp_col, do_calc_f )

      end do ! sv_i

    end do ! sps_i

  end subroutine Get_d_Deltau_Pol_dF

! ------------------------------------------  Get_d_Deltau_Pol_dT  -----

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

  subroutine Get_d_Deltau_Pol_dT ( CT, STCP, STSP, Tan_PT, T_Path, &
                & Alpha_Path, dAlpha_dT_path, dAlpha_dT_polarized_path, &
                & Eta_zxp, P_Stop, Del_S, GL_Inds, Del_Zeta, Do_GL, &
                & ds_dh, dh_dz_gw, ds_dz_gw, Incoptdepth, Ref_cor, &
                & H_path, dH_dt_path, H_tan, tan_pt_f, &
                & Deriv_Flags, D_Deltau_Pol_DT )

    use DExdT_m, only: dExdT
    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP, IP
    use Opacity_m, only: Opacity
    use Rad_Tran_m, only: Get_do_Calc
    use Sparse_m, only: Sparse_t

  ! Arguments
    ! SVE == # of state vector elements
    real(rp), intent(in) :: CT(:)           ! Cos(Theta), where theta
      ! is the angle between the line of sight and magnetic field vectors.
    real(rp), intent(in) :: STCP(:)         ! Sin(Theta) Cos(Phi) where
      ! theta is as for CT and phi (for this purpose only) is the angle
      ! between the plane defined by the line of sight and the magnetic
      ! field vector, and the "instrument field of view plane polarized"
      ! (IFOVPP) X axis.
    real(rp), intent(in) :: STSP(:)         ! Sin(Theta) Sin(Phi)
    integer, intent(in) :: Tan_PT           ! tangent index along the coarse
                                            ! path, usually N_Path/2
    real(rp), intent(in) :: T_Path(:)       ! path temperatures on composite
                                            ! coarse & fine grid
    complex(rp), intent(in), target :: Alpha_Path(-1:,:) ! -1:1 x path on 
                                            ! composite coarse & fine path
    complex(rp), intent(in), target :: dAlpha_dT_polarized_path(-1:,:) ! -1:1 x path on 
                                            ! composite coarse & fine path path
    real(rp), intent(in), target :: dAlpha_dT_path(:) ! nonpolarized dAlpha_dT
                                            ! on composite coarse & fine path
    class(sparse_t), intent(in) :: Eta_ZXP  ! Interpolating coefficients from
                                            ! state vector to path
    integer, intent(in) :: P_Stop           ! Where to stop on coarse path
    real(rp), intent(in) :: Del_S(:)        ! unrefracted path length.  This
                                            ! is for the whole coarse path, not
                                            ! just the part up to the black-out
    integer(ip), intent(in) :: GL_Inds(:)   ! indices for reading fine path
                                            ! elements from sps_path
    real(rp), intent(in) :: Del_Zeta(:)     ! path -log(P) differences on the
                                            ! main grid.  This is for the whole
                                            ! coarse path, not just the part up
                                            ! to the black-out
    logical, intent(in) :: Do_gl(:)         ! Indicates where on the coarse path
                                            ! to do gl integrations.
    real(rp), intent(in) :: ds_dh(:)        ! path length wrt height derivative
                                            ! on entire grid.  Only the
                                            ! gl_inds part is used.
    real(rp), intent(in) :: dh_dz_gw(:)     ! path height wrt zeta derivative * gw.
                                            ! on entire grid.  Only the
                                            ! gl_inds part is used.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative * gw.
                                            ! on entire grid.  Only the
                                            ! gl_inds part is used.
    complex(rp), intent(in) :: Incoptdepth(:,:,:) ! negative of incremental
                                            ! optical depth.  2 x 2 x path
    real(rp), intent(in) :: Ref_cor(:)      ! refracted to unrefracted path
                                            ! length ratios.

    ! For hydrostatic
    real(rp), intent(in), target :: H_path(:) ! path heights + req on composite
                                            ! coarse & fine grid
    type(sparse_t), intent(in) :: dH_dt_path ! derivative of fine path
                                            ! height wrt temperature(km/K) on
                                            ! composite coarse & fine grid
    real(rp), intent(in) :: H_tan           ! tangent height + req (km).
    integer, intent(in) :: tan_pt_f         ! tangent point in dh_dt_path
    logical, intent(in) :: Deriv_flags(:)   ! Indicates which temperature
                                            ! derivatives to do

! Outputs

    complex(rp), intent(out) :: D_Deltau_Pol_DT(:,:,:,:) ! 2 x 2 x path x sve.
                                            ! derivative of delta Tau wrt
                                            ! temperature state vector
                                            ! element. (K)

  ! Local variables
    integer :: A                     ! Index for GL points
    complex(rp), pointer :: Alpha_Path_c(:,:) ! -1:1 x path on coarse grid
    complex(rp):: D_Alpha_DT_eta(-1:1,p_stop) ! Singularity * Del_S
    real(rp), pointer :: dAlpha_dT_path_c(:) ! nonpolarized dAlpha_dT on
                                     ! coarse path
    complex(rp), pointer :: dAlpha_dT_polarized_path_c(:,:)
    real(rp), target :: dH_dt_path_col(dh_dt_path%nRows)
    real(rp), pointer :: dH_dt_path_c(:) ! derivative of coarse path height
                                     ! wrt temperature(km/K) on coarse path.
    real(rp) :: dh_dt_tan            ! dh_dt_path at the tangent point
    complex(rp) :: D_Incoptdepth_dT(2,2,p_stop)
    logical :: Do_Calc(1:p_stop) ! do_calc_t_c .or. ( do_gl .and.
                                     ! any of the corresponding do_calc_t_f ).
    logical :: Do_Calc_Col(eta_zxp%nRows)  ! Where Eta_ZXP /= 0
    logical :: Do_Calc_Col_f ( size(gl_inds) )
    logical, target :: Do_Calc_Hyd(dh_dt_path%nRows) ! On the whole path
    logical, pointer :: Do_Calc_Hyd_c(:)             ! On the coarse path
    real(rp), target :: Eta_ZXP_Col(eta_zxp%nRows) ! One column of Eta_ZXP
    real(rp), pointer :: Eta_ZXP_Col_C(:)
    real(rp) :: F(ng)                ! Factor in GL that doesn't depend on
                                     ! sigma +/- or pi.
    real(rp) :: Fa, Fb               ! Hydrostatic integrand at ends of
                                     ! path segment
    integer :: GA                    ! GL_inds(a)
    real(rp), pointer :: H_path_c(:) ! path heights + req (km) on the coarse
                                     ! grid. This is for the whole coarse path,
                                     ! not just the part up to the black-out
    integer :: H_Stop                ! Stop point for hydrostatic parts
    integer :: I_stop                ! Stop point, which may be before N_Path
    integer :: L
    integer :: N_Path                ! Total coarse path length.
    logical :: NeedFA                ! Need FA in hydrostatic calculation
    integer :: P_i                   ! Index on the coarse path
!   integer :: P_Stop_f              ! P_Stop on fine path
    real(rp) :: S_Del_S              ! Sum of Del_S
    complex(rp) :: Singularity(-1:1,p_stop) ! n/T * Alpha * Eta
                                     ! on the coarse path
    integer :: SV_i                  ! Index of state vector element

    i_stop = p_stop
    n_path = size(del_zeta)
    h_path_c => h_path ( 1 :: ngp1 )
    alpha_path_c(-1:,1:) => alpha_path ( :, 1 :: ngp1 )
    dAlpha_dT_polarized_path_c(-1:,1:) => dAlpha_dT_polarized_path ( :, 1 :: ngp1 )
    dh_dt_path_col = 0
    dH_dt_path_c => dH_dt_path_col ( 1 :: ngp1 )
    dAlpha_dT_path_c => dAlpha_dT_path ( 1 :: ngp1 )
    do_calc_col = .false.
    do_calc_hyd_c => do_calc_hyd ( 1 :: ngp1 )
    do_calc_hyd_c = .false.
    eta_zxp_col = 0
    eta_zxp_col_c => eta_zxp_col(1::ngp1)
!   p_stop_f = (p_stop - 1) * ngp1 + 1

    do sv_i = 1, size(eta_zxp%cols,1) ! state vector elements
      if ( .not. deriv_flags(sv_i)) then
        d_deltau_pol_dT(:,:,:,sv_i) = 0.0_rp
        cycle
      end if

      call dh_dt_path%get_col ( sv_i, dh_dt_path_col, do_calc_hyd )
      dh_dt_tan = dh_dt_path_col(tan_pt_f)

      call eta_zxp%get_col ( sv_i, eta_zxp_col, do_calc_col ) !, last=p_stop_f )
      do_calc_col_f = do_calc_col(gl_inds) ! only the GL points, no coarse points

      ! do the absorption part
      ! combine non zeros flags for both the main and gl parts
      ! Add in contribution from scalar model, 0.25 for +/- sigma,
      ! 0.5 for pi.

      call Get_do_calc ( do_calc_col(::ngp1), do_calc_col_f, do_gl, &
        & do_calc )

      a = 1
      do p_i = 1, i_stop
        if ( do_calc(p_i) ) then
          !{ d\_alpha\_dT\_eta$_i = 
          !    \frac{\partial \alpha}{\partial T_i} \text{d}s=
          !    \frac{\partial \alpha}{\partial T}
          !    \frac{\partial T}{\partial T_i} \text{d}s =
          !    \frac{\partial \alpha}{\partial T} \eta_i \text{d}s$
          singularity(:,p_i) = eta_zxp_col_c(p_i) * &
            & ( dAlpha_dT_polarized_path_c(:,p_i) + &
            &   (/ 0.25, 0.50, 0.25 /) * dAlpha_dT_path_c(p_i) )
          d_alpha_dT_eta(:,p_i) = singularity(:,p_i) * del_s(p_i)
        else
          singularity(:,p_i) = 0.0_rp
          d_alpha_dT_eta(:,p_i) = 0.0_rp
        end if
        ! Do GL if needed here
        if ( do_gl(p_i) ) then
          ga = gl_inds(a)
          if ( do_calc(p_i) ) then
            f = ds_dz_gw(ga:ga+ng-1)
            do l = -1, 1
              d_alpha_dT_eta(l,p_i) = d_alpha_dT_eta(l,p_i) + &
                 & del_zeta(p_i) * &
                 & sum( ( ( dAlpha_dT_polarized_path(l,ga:ga+ng-1) + &
                 &          (0.5-0.25*abs(l)) * dAlpha_dT_path(ga:ga+ng-1) ) * &
                 &        eta_zxp_col(ga:ga+ng-1) - &
                 &        singularity(l,p_i) ) * f )
            end do ! l
          end if
          a = a + ng
        end if
      end do ! p_i

      ! Now do the hydrostatic part
      !??? This only goes as far as I_Stop.  We may     ???
      !??? want to do it this way in DRad_Tran_dT also. ???

      ! First combine boundary flags
      do_calc = do_calc_hyd_c
      if ( i_stop < tan_pt ) then           
        do_calc(2:i_stop) = do_calc(2:i_stop) .or. do_calc(1:i_stop-1)
        h_stop = i_stop
      else
        do_calc(2:tan_pt) = do_calc(tan_pt) .or. do_calc(2:tan_pt) .or. do_calc(1:tan_pt-1)
        h_stop = tan_pt - 1
      end if
      do_calc(1) = .false.
      s_del_s = sum(del_s(2:tan_pt)) ! Yes, this goes to the midpoint of the coarse path
      needFA = .true.
      fa = 0.0_rp ! In case n_path <= 4
      do p_i = 2 , h_stop
        if ( do_calc(p_i) ) then
          if ( needFA ) then
            fa = (h_path_c(p_i-1) * dh_dt_path_c(p_i-1) &
            &   - h_tan * dh_dt_tan) / s_del_s
            needFA = .false.
          end if
          s_del_s = s_del_s - del_s(p_i)
          fb = (h_path_c(p_i) * dh_dt_path_c(p_i) &
            & - h_tan * dh_dt_tan) / s_del_s
          d_alpha_dT_eta(:,p_i) = d_alpha_dT_eta(:,p_i) + alpha_path_c(:,p_i) * (fa - fb)
          fa = fb
        else
          s_del_s = s_del_s - del_s(p_i)
        end if
      end do ! p_i

      ! special processing at tangent.  fb is zero

      if ( i_stop >= tan_pt ) then
        if ( do_calc(tan_pt) ) &
          & d_alpha_dT_eta(:,tan_pt) = d_alpha_dT_eta(:,tan_pt) + alpha_path_c(:,tan_pt) * fa
      end if
      if ( i_stop > tan_pt + 1 ) then ! tan_pt+1 instead of tan_pt so that tan_pt+2 will be
                                      ! in bounds if i_stop == 2.

        do_calc(tan_pt+1:i_stop-1) = do_calc(tan_pt+1:i_stop-1) .or. do_calc(tan_pt+2:i_stop) .or. do_calc(tan_pt+1)
        if ( i_stop == n_path ) then
          h_stop = i_stop - 1
          do_calc(i_stop) = .false.
        else
          h_stop = i_stop
          do_calc(i_stop) = do_calc(i_stop) .or. do_calc(tan_pt+1)
        end if

        needFA = .not. do_calc(tan_pt+1)
        s_del_s = del_s(tan_pt+1)
        if ( do_calc(tan_pt+1) ) then
          fa = (h_path_c(tan_pt+2) * dh_dt_path_c(tan_pt+2) &
            & - h_tan * dh_dt_tan) / s_del_s
          d_alpha_dT_eta(:,tan_pt+1) = d_alpha_dT_eta(:,tan_pt+1) + alpha_path_c(:,tan_pt+1) * fa
        end if

        s_del_s = del_s(tan_pt+1)
        do p_i = tan_pt + 2, h_stop
          if ( do_calc(p_i) ) then
            if ( needFA ) then
              fa = (h_path_c(p_i) * dh_dt_path_c(p_i) &
                & - h_tan * dh_dt_tan) / s_del_s
              needFA = .false.
            end if
            s_del_s = s_del_s + del_s(p_i)
            fb = (h_path_c(p_i+1) * dh_dt_path_c(p_i+1) &
              & - h_tan * dh_dt_tan) / s_del_s
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
      a = 1
      do p_i = 1, i_stop             ! along the path
        if ( do_gl(p_i) ) then
          ga = gl_inds(a)
          ! Don't test do_calc: There might be GL corrections even if
          ! dh_dt_path_c (from whence came do_calc) is zero.
          f = (((2.0_rp*h_path(ga:ga+ng-1)**2 - 3.0_rp*h_tan**2)           &
            &   * dh_dt_path_col(ga:ga+ng-1) +                             &
            &   h_path(ga:ga+ng-1) * h_tan * dh_dt_tan) /                  &     
            &  (sqrt(h_path(ga:ga+ng-1)**2 - h_tan**2))**3                 &
            &  + eta_zxp_col(ga:ga+ng-1) * ds_dh(ga:ga+ng-1) /             &
            &  t_path(ga:ga+ng-1)) * dh_dz_gw(ga:ga+ng-1)
          do l = -1, 1
            d_alpha_dT_eta(l,p_i) = d_alpha_dT_eta(l,p_i) + &
               & del_zeta(p_i) * &
               & sum( ( alpha_path(l,ga:ga+ng-1) - alpha_path_c(l,p_i) ) * f )
          end do ! l
          a = a + ng
        end if

        d_alpha_dT_eta(:,p_i) = d_alpha_dT_eta(:,p_i) * ref_cor(p_i)

      end do ! p_i

      !{ {\tt d\_alpha\_dT\_eta} is now really
      ! $\int \frac{\partial \Delta \delta}{\partial T} \,\text{d}s$,
      ! where $\Delta \delta$ is the incremental opacity.
      !%
      ! Compute {\tt d\_incoptdepth\_dT} =
      ! $\frac{\partial \int {\bf G} \,\text{d}s}{\partial T}$
      call opacity ( ct, stcp, stsp, d_alpha_dT_eta, d_incoptdepth_dT )

      !{ Compute $\frac{\partial \bf E}{\partial T} =
      !           \frac{\partial \, \exp(-{\int \bf G}\, \text{d}s)}{\partial T}$
      ! where {\bf G} is the incremental optical depth matrix.
      do p_i = 1, i_stop             ! along the path
        if ( any( d_alpha_dT_eta(:,p_i) /= 0.0_rp ) ) then
          call dExdT ( incoptdepth(:,:,p_i), -d_incoptdepth_dT(:,:,p_i), &
                     & d_deltau_pol_dT(:,:,p_i,sv_i) ) ! d exp(incoptdepth) / dT
        else
          d_deltau_pol_dT(:,:,p_i,sv_i) = 0.0_rp
        end if
      end do ! p_i

      call dh_dt_path%clear_col ( sv_i, dh_dt_path_col, do_calc_hyd )
      call eta_zxp%clear_col ( sv_i, eta_zxp_col, do_calc_col ) !, last=p_stop_f )

    end do ! sv_i

  end subroutine Get_d_Deltau_Pol_dT

!-----------------------------------------------------------------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Get_d_Deltau_Pol_M

! $Log$
! Revision 2.49  2023/07/06 18:36:42  pwagner
! correctedd off-by-one errors; limit amount of debug printing
!
! Revision 2.48  2023/06/23 20:44:34  pwagner
! In middle of debugging pol fwdmdl
!
! Revision 2.47  2018/11/19 21:53:33  vsnyder
! Correct confusion between coarse and fine path indexing in
! Get_d_Deltau_Pol_dF, which has apparently so far not been used.
!
! Revision 2.46  2018/05/24 03:24:36  vsnyder
! Use sparse representation for dh_dt_path
!
! Revision 2.45  2018/05/14 23:40:58  vsnyder
! Change to sparse eta representation
!
! Revision 2.44  2017/08/09 20:53:13  vsnyder
! Decide which panels to integrate for the mixing-ration case in a way that
! is similar to how it's done for the temperature case.  The old way was wrong.
!
! Revision 2.43  2013/06/12 02:20:19  vsnyder
! Cruft removal
!
! Revision 2.42  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.41  2011/03/11 03:08:00  vsnyder
! Only use the nonzeros in d_delta_df
!
! Revision 2.40  2010/08/19 02:11:16  vsnyder
! Change some variable names, organize some stuff more like dRad_tran_df
!
! Revision 2.39  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.38  2007/11/08 02:02:40  vsnyder
! Remove unused dummy argument
!
! Revision 2.37  2006/12/13 02:32:02  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.36  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.33  2006/04/21 22:12:31  vsnyder
! Fix final bug in mixing ratio derivatives
!
! Revision 2.32  2006/04/19 23:00:48  bill
! I think I fixed the vmr derivative bug
!
! Revision 2.31  2006/04/11 18:36:21  vsnyder
! Include missing factor of tanh(h nu / k T)
!
! Revision 2.30  2006/03/17 00:41:12  vsnyder
! Use ubound instead of size for grids_f%l_v
!
! Revision 2.29  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.28  2004/11/01 20:24:32  vsnyder
! Reorganization of representation for molecules and beta groups
!
! Revision 2.27  2004/07/29 02:31:17  vsnyder
! Simplify by avoiding PFA-related data structures
!
! Revision 2.26  2004/04/19 21:00:53  vsnyder
! Remove unreferenced USE names
!
! Revision 2.25  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.24  2004/04/02 01:00:20  vsnyder
! Inching toward analytic temperature derivatives
!
! Revision 2.23  2004/03/08 22:56:41  vsnyder
! Remove calculation of complex exponent for T/T0 for beta.  Remove
! D_Delta_DT, which had been gotten from drad_tran_dt but is no longer
! needed (we use alpha_xn_path_? instead).
!
! Revision 2.22  2004/02/04 01:18:45  vsnyder
! Remove accidentally-checked-in test/debug stuff
!
! Revision 2.21  2004/02/03 02:47:55  vsnyder
! Progress (hopefully) on polarized temperature derivatives
!
! Revision 2.20  2003/12/03 00:25:32  vsnyder
! Corrections to hydrostatic calculation
!
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
