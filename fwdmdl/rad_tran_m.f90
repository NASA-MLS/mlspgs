! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RAD_TRAN_M

  implicit NONE
  private
  public :: RAD_TRAN, RAD_TRAN_POL, DRAD_TRAN_DF, DRAD_TRAN_DT, DRAD_TRAN_DX
  public :: Get_Do_Calc
  private ::  Get_Do_Calc_Indexed, Get_Inds

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName = "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!------------------------------------------------------  Rad_tran  -----
! This is the radiative transfer model, radiances only !

  subroutine Rad_tran ( gl_inds, more_inds, e_rflty, del_zeta, &
                     &  alpha_path_c, ref_cor, do_gl, incoptdepth, &
                     &  alpha_path_gl, ds_dz_gw, t_script, &
                     &  tau, rad, i_stop )

    use GLNP, only: NG
    use MLSCommon, only: RP, IP
    use SCRT_DN_M, ONLY: SCRT_DN

  ! inputs

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
    logical, intent(in) :: do_gl(:)          ! path flag indicating where to do
  !                                            gl integrations.
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
          & sum( (alpha_path_gl(a:a+ng-1) -  alpha_path_c(more_inds(i))) * &            
               & ds_dz_gw(aa:aa+ng-1) ) 
        a = a + ng
      end do ! i

    end if

    incoptdepth = ref_cor * incoptdepth

    call scrt_dn ( t_script, e_rflty, incoptdepth, tau, rad, i_stop )

  end subroutine Rad_tran

!--------------------------------------------------  Rad_Tran_Pol  -----

  subroutine Rad_tran_Pol ( gl_inds, more_inds, e_rflty, del_zeta, alpha_path_c, &
                     &  ref_cor, do_gl, incoptdepth_pol, deltau_pol, alpha_path_gl, &
                     &  ds_dz_gw, ct, stcp, stsp, t_script, &
                     &  prod_pol, tau_pol, rad_pol, p_stop )

    ! Polarized radiative transfer.  Radiances only, no derivatives.

    use CS_Expmat_m, only: CS_Expmat
    use DO_DELTA_M, ONLY: POLARIZED_PATH_OPACITY
    use GLNP, ONLY: Ng
    use MCRT_M, ONLY: MCRT
    use MLSCommon, only: RP, IP
    use Opacity_m, only: Opacity

  ! inputs

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
    logical, intent(in) :: do_gl(:)          ! path flag indicating where to do
      !              gl integrations.
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

    end if

    ! At this point, incoptdepth_pol(:,:,1:npc/2) should be nearly
    ! identical to incoptdepth_pol(:,:,1:npc/2+1) (npc/2 is the
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
      & p_stop, prod_pol, tau_pol, rad_pol )

    return

  ! Error exit if cs_expmat detected an overflow
  99 p_stop = - p_stop - 1

  end subroutine Rad_tran_Pol
!--------------------------------------------------  DRad_tran_df  -----
! This is the radiative transfer derivative wrt mixing ratio model

  subroutine DRad_tran_df ( indices_c, gl_inds, del_zeta, Grids_f, &
                         &  beta_path_c, eta_zxp_f, sps_path, do_calc_f, &
                         &  beta_path_f, do_gl, del_s, ref_cor, &
                         &  ds_dz_gw, t_script, tau, &
                         &  i_stop, d_delta_df, drad_df )

    use GLNP, only: NG
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSCommon, only: RP, IP
    use SCRT_DN_M, ONLY: GET_DSCRT_NO_T_DN
    use Where_M, only: Where

! Inputs

    integer(ip), intent(in) :: indices_c(:) ! coarse grid indicies
    integer(ip), intent(in) :: gl_inds(:)   ! Gauss-Legendre grid indicies
    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    real(rp), intent(in) :: beta_path_c(:,:) ! cross section for each species
      !                                        on coarse grid.
    real(rp), intent(in) :: eta_zxp_f(:,:)   ! representation basis function.
    real(rp), intent(in) :: sps_path(:,:)    ! Path species function.
    logical, intent(in) :: do_calc_f(:,:)    ! A logical indicating where the
      !              representation basis function is not zero.
    real(rp), intent(in) :: beta_path_f(:,:) ! cross section for each species
      !              on gl grid.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !              do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !              length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: t_script(:)      ! differential temperatures (K).
    real(rp), intent(in) :: tau(:)           ! transmission function.
    integer(ip), intent(in) :: i_stop        ! path stop index

! Outputs

    real(rp), intent(out) :: d_delta_df(:,:) ! path x sve.  derivative of delta
      !              wrt mixing ratio state vector element. (K)
    real(rp), intent(out) :: drad_df(:)      ! derivative of radiances wrt
      !              mixing ratio state vector element. (K)
! Internals

    integer(ip) :: i_start, n_inds, no_mol, no_to_gl, sps_i, sv_i
    integer(ip) :: a, aa, i
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                         ! on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                         ! Indices on the coarse path where GL
                                         ! corrections get applied.

    real(rp) :: singularity(1:size(tau)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.
    logical :: do_calc(1:size(tau))      ! Flags on coarse path where do_calc_c
                                         ! or (do_gl and any corresponding
                                         ! do_calc_f).

! Begin code

    no_mol = ubound(Grids_f%l_z,1)

    drad_df(:) = 0.0_rp

    do sps_i = 1, no_mol

      do sv_i = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

! Skip the masked derivatives, according to the l2cf inputs

        if ( .not. Grids_f%deriv_flags(sv_i) ) cycle

        d_delta_df(:,sv_i) = 0.0_rp

        call get_do_calc_indexed ( do_calc_f(:,sv_i), indices_c, gl_inds, &
          & do_gl, do_calc )

! find where the non zeros are along the path

        n_inds = count(do_calc)
        if ( n_inds > 0 ) then

          inds => inds_B(1:n_inds)

          call where ( do_calc, inds )

          no_to_gl = count(do_gl(inds))

          if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

            all_inds => all_inds_B(1:ng*no_to_gl)
            more_inds => more_inds_B(1:no_to_gl)

            call get_inds ( do_gl, do_calc, more_inds, all_inds )
          end if

          if ( grids_f%lin_log(sps_i) ) then

            do i = 1, n_inds ! Don't trust the compiler to fuse loops
              singularity(inds(i)) = beta_path_c(inds(i),sps_i) &
                        & * eta_zxp_f(indices_c(inds(i)),sv_i) &
                        & * sps_path(indices_c(inds(i)),sps_i)
              d_delta_df(inds(i),sv_i) = singularity(inds(i)) * del_s(inds(i))
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

            a = 1
            do i = 1, no_to_gl
              aa = all_inds(a)
              d_delta_df(more_inds(i),sv_i) = d_delta_df(more_inds(i),sv_i) + &
                & del_zeta(more_inds(i)) * &
                & sum( (beta_path_f(aa:aa+ng-1,sps_i) * &
                     &  eta_zxp_f(gl_inds(aa:aa+ng-1),sv_i) * &
                     &  sps_path(gl_inds(aa:aa+ng-1),sps_i) - &
                     &  singularity(more_inds(i))) * &
                     & ds_dz_gw(gl_inds(aa:aa+ng-1)) )
              a = a + ng
            end do
            d_delta_df(inds,sv_i) = ref_cor(inds) * d_delta_df(inds,sv_i) * &
                                  & exp(-grids_f%values(sv_i))

          else

            do i = 1, n_inds ! Don't trust the compiler to fuse loops
              singularity(inds(i)) = beta_path_c(inds(i),sps_i) &
                        & * eta_zxp_f(indices_c(inds(i)),sv_i)
              d_delta_df(inds(i),sv_i) = singularity(inds(i)) * del_s(inds(i))
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

            a = 1
            do i = 1, no_to_gl
              aa = all_inds(a)
              d_delta_df(more_inds(i),sv_i) = d_delta_df(more_inds(i),sv_i) + &
                & del_zeta(more_inds(i)) * &
                & sum( (beta_path_f(aa:aa+ng-1,sps_i) * &
                     &  eta_zxp_f(gl_inds(aa:aa+ng-1),sv_i) - &
                     &  singularity(more_inds(i))) * &
                     & ds_dz_gw(gl_inds(aa:aa+ng-1)) )
              a = a + ng
            end do
            d_delta_df(inds,sv_i) = ref_cor(inds) * d_delta_df(inds,sv_i)

          end if

          i_start = MIN(inds(1),i_stop)

          call get_dscrt_no_t_dn ( d_delta_df(:,sv_i), t_script, tau, &
                                &  i_start, i_stop, drad_df(sv_i))

        end if

      end do

    end do

  end subroutine drad_tran_df

!--------------------------------------------------  drad_tran_dx  -----
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )

  subroutine DRad_tran_dx ( del_zeta, Grids_f, dbeta_path_c, eta_zxp_f_c, &
                         &  sps_path_c, do_calc_f_c, dbeta_path_f, eta_zxp_f_f, &
                         &  sps_path_f, do_calc_f_f, do_gl, gl_inds, del_s, &
                         &  ref_cor, ds_dz_gw, t_script, tau, i_stop, drad_dx )

    use GLNP, only: NG
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSCommon, only: RP, IP
    use SCRT_DN_M, ONLY: GET_DSCRT_NO_T_DN
    use Where_M, only: Where

! Inputs

    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T) :: Grids_f                ! All the coordinates
    real(rp), intent(in) :: dbeta_path_c(:,:) ! derivative of beta wrt dx
!                                              on main grid.
    real(rp), intent(in) :: eta_zxp_f_c(:,:) ! representation basis function
!                                              main grid.
    real(rp), intent(in) :: sps_path_c(:,:)  ! species function on  main grid.
    logical, intent(in) :: do_calc_f_c(:,:)  ! A logical indicating where the
!                                              representation basis function is
!                                              not zero on main grid.
    real(rp), intent(in) :: dbeta_path_f(:,:) ! derivative of beta wrt dx
!                                              on gl grid.
    real(rp), intent(in) :: eta_zxp_f_f(:,:) ! representation basis function
!                                              gl grid.
    real(rp), intent(in) :: sps_path_f(:,:)  ! species function on gl grid.
    logical, intent(in) :: do_calc_f_f(:,:)  ! A logical indicating where the
!                                              representation basis function is
!                                              not zero on main grid.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
!                                              do rec_tan_inds gl integrations
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indicies
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
!                                              length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative * gw
!                                              on entire grid.
    real(rp), intent(in) :: t_script(:)      ! differential temperatures (K).
    real(rp), intent(in) :: tau(:)           ! transmission function.
    integer(ip), intent(in) :: i_stop        ! path stop index

! Outputs

    real(rp), intent(out) :: drad_dx(:)      ! derivative of radiances wrt x
!                                              state vector element. (K)
! Internals

    integer(ip) :: A, AA, I
    integer(ip) :: i_start, n_inds, no_mol, no_to_gl, sps_i, sv_i, sv_j
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                         ! on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                         ! Indices on the coarse path where GL
                                         ! corrections get applied.

    real(rp) :: d_delta_dx(1:size(tau))
    real(rp) :: singularity(1:size(tau)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.

    logical :: do_calc(1:size(tau))      ! Flags on coarse path where do_calc_c
                                         ! or (do_gl and any corresponding
                                         ! do_calc_f).

! Begin code

    no_mol = ubound(grids_f%l_z,1)

    sv_i = 0
    drad_dx(:) = 0.0_rp

    do sps_i = 1 , no_mol

      do sv_j = 1 , (Grids_f%l_z(sps_i) - Grids_f%l_z(sps_i-1)) * &
        &           (Grids_f%l_p(sps_i) - Grids_f%l_p(sps_i-1))

        sv_i = sv_i + 1
        d_delta_dx = 0.0_rp

        call get_do_calc ( do_calc_f_c(:,sv_i), do_calc_f_f(:,sv_i), do_gl, &
          & do_calc )

! find where the non zeros are along the path

        n_inds = count(do_calc)
        if ( n_inds > 0 ) then

          inds => inds_B(1:n_inds)

          call where ( do_calc, inds )

          do i = 1, n_inds ! Don't trust the compiler to fuse loops
            singularity(inds(i)) = dbeta_path_c(inds(i),sps_i) * eta_zxp_f_c(inds(i),sv_i)  &
                       &  * sps_path_c(inds(i),sps_i)
            d_delta_dx(inds(i)) = singularity(inds(i)) * del_s(inds(i))
          end do ! i

          no_to_gl = count(do_gl(inds))
          if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

            all_inds => all_inds_B(1:Ng*no_to_gl)
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

            a = 1
            do i = 1, no_to_gl
              aa = all_inds(a)
              d_delta_dx(more_inds(i)) = d_delta_dx(more_inds(i)) + &
                & del_zeta(more_inds(i)) * &
                & sum( (dbeta_path_f(aa:aa+ng-1,sps_i) * &
                     &  eta_zxp_f_f(aa:aa+ng-1,sv_i) * &
                     &  sps_path_f(aa:aa+ng-1,sps_i) - &
                     &  singularity(more_inds(i))) * &
                     & ds_dz_gw(gl_inds(aa:aa+ng-1)) )
              a = a + ng
            end do

          end if

          d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)

          i_start = min(inds(1),i_stop)

          call get_dscrt_no_t_dn ( d_delta_dx, t_script, tau, i_start, i_stop, &
                                &  drad_dx(sv_i) )

        end if

      end do

    end do

  end subroutine DRad_tran_dx

!{\newpage

!--------------------------------------------------  drad_tran_dt  -----
! This is the radiative transfer derivative wrt temperature model

  subroutine DRad_tran_dt ( del_zeta, h_path_c, t_path_c, dh_dt_path_c, &
                         &  alpha_path_c, alphaxn_path_c, eta_zxp_c, do_calc_t_c, &
                         &  do_calc_hyd_c, del_s, ref_cor, h_tan, dh_dt_tan, &
                         &  do_gl, gl_inds, h_path_f, t_path_f, dh_dt_path_f, &
                         &  alpha_path_f, alphaxn_path_f, eta_zxp_f, do_calc_t_f, &
                         &  ds_dh, dh_dz_gw, ds_dz_gw, t_script, dt_scr_dt, &
                         &  tau, i_stop, deriv_flags, drad_dt )

    use GLNP, only: NG
    use MLSCommon, only: RP, IP
    use SCRT_DN_M, ONLY: GET_DSCRT_DN
    use Where_M, only: Where

! Inputs

    real(rp), intent(in) :: del_zeta(:)     ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: h_path_c(:)     ! path heights + req on main grid km.
    real(rp), intent(in) :: t_path_c(:)     ! path temperature(K) on main grid.
    real(rp), intent(in) :: dh_dt_path_c(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on main grid.
    real(rp), intent(in) :: alpha_path_c(:) ! path absorption(km^-1)
!                                             on main grid.
    real(rp), intent(in) :: alphaxn_path_c(:) ! path absorption * temperature
!                                             exponent (km^-1) on main grid.
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
    real(rp), intent(in) :: alphaxn_path_f(:) ! path absorption * temperature
!                                             exponent (km^-1) on gl grid.
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
    real(rp), intent(in) :: t_script(:)     ! differential temperatures (K).
    real(rp), intent(in) :: dt_scr_dt(:,:)  ! d t_script / d T * d T / d eta.
    real(rp), intent(in) :: tau(:)          ! transmission function.
    integer(ip), intent(in) :: i_stop       ! path stop index
    logical, intent(in) :: deriv_flags(:)   ! Indicates which temperature
!                                             derivatives to do

! Outputs
    real(rp), intent(out) :: drad_dt(:)     ! derivative of radiances wrt
!                                             mixing ratio state vector
!                                             element. (K)
! Internals

    integer(ip) :: A, AA, B, BB
    integer(ip) :: i, j, i_start, mid, n_inds, n_path, no_to_gl, p_i, sv_i
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part of all_inds_B;
                                         ! Indices on GL grid for stuff
                                         ! used to make GL corrections
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B;  Indices
                                         ! on coarse path where do_calc.
    integer(ip), pointer :: more_inds(:) ! more_inds => part of more_inds_B;
                                         ! Indices on the coarse path where GL
                                         ! corrections get applied.

    real(rp) :: d_delta_dt(size(eta_zxp_c,1),size(eta_zxp_c,2)) ! incremental
                                         ! opacity derivatives
    real(rp) :: fa, fb
    real(rp) :: S_DEl_S                  ! Running sum of Del_S
    real(rp) :: singularity(1:size(del_zeta)) ! integrand on left edge of coarse
                                         ! grid panel -- singular at tangent pt.

    logical :: do_calc(1:size(del_zeta)) ! do_calc_t_c .or. ( do_gl .and. any
                                         ! of the corresponding do_calc_t_f ).
    logical :: NeedFA                    ! Need F(A) for hydrostatic

! Begin code

    n_path = size(del_zeta)
    mid = n_path / 2

! compute the opacity derivative singularity value

    d_delta_dt = 0.0_rp
    drad_dt(:) = 0.0_rp

    do sv_i = 1 , size(eta_zxp_c,dim=2)
      if ( .not. deriv_flags(sv_i)) cycle
      i_start = 1

! do the absorption part
! combine non zeros flags for both the main and gl parts

      call get_do_calc ( do_calc_t_c(:,sv_i), do_calc_t_f(:,sv_i), do_gl, &
        & do_calc )

! find where the non zeros are along the path

      n_inds = count(do_calc)
      if ( n_inds > 0 ) then
        inds => inds_B(1:n_inds)

        call where ( do_calc, inds )
        i_start = max(inds(1)-1,1)

        do i = 1, n_inds ! Don't trust the compiler to fuse loops
          j = inds(i)
          singularity(j) = alphaxn_path_c(j) * eta_zxp_c(j,sv_i) / t_path_c(j)
          d_delta_dt(j,sv_i) = singularity(j) * del_s(j)
        end do ! i
! see if anything needs to be gl-d

        no_to_gl = count(do_gl(inds))
        if ( no_to_gl > 0 ) then

          all_inds => all_inds_B(1:ng*no_to_gl)
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

          a = 1
          do i = 1, no_to_gl
            aa = all_inds(a)
            bb = aa + ng - 1
            j = more_inds(i)
            d_delta_dt(j,sv_i) = d_delta_dt(j,sv_i) + &
              & del_zeta(j) * &
              & sum( (alphaxn_path_f(aa:bb) * &
                   &  eta_zxp_f(aa:bb,sv_i) / &
                   &  t_path_f(aa:bb) - &
                   &  singularity(j)) * &
                   & ds_dz_gw(gl_inds(aa:bb)) )
            a = a + ng
          end do

        end if

      end if

! now do the hydrostatic part
! combine boundaries flags

      do_calc = do_calc_hyd_c(:,sv_i)
      do_calc(2:mid) =          do_calc(mid)   .or. do_calc(2:mid)          .or. do_calc(1:mid-1)
      do_calc(mid+1:n_path-1) = do_calc(mid+1) .or. do_calc(mid+1:n_path-1) .or. do_calc(mid+2:n_path)

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
        s_del_s = sum(del_s(2:mid))
        do p_i = 2 , mid - 1
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

        if ( do_calc(mid) ) then
          d_delta_dt(mid,sv_i) = d_delta_dt(mid,sv_i) + alpha_path_c(mid) * fa
          inds(i) = mid
          i = i + 1
        end if

        needFA = .not. do_calc(mid+1)
        s_del_s = del_s(mid+1)
        if ( do_calc(mid+1) ) then
          fa = (h_path_c(mid+2) * dh_dt_path_c(mid+2,sv_i) - &
              & h_tan * dh_dt_tan(sv_i)) / s_del_s
          d_delta_dt(mid+1,sv_i) = d_delta_dt(mid+1,sv_i) + &
            &                      alpha_path_c(mid+1) * fa
          inds(i) = mid + 1
          i = i + 1
        end if

        ! Several subscripts in this loop are offset by 1 from the nearly-
        ! identical loop above, and we use (fb-fa) here instead of (fa-fb).
        do p_i = mid + 2, n_path - 1
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

        i_start = min(i_start,inds(1))

      end if

! correct the whole thing for path length refraction

      d_delta_dt(:,sv_i) = ref_cor * d_delta_dt(:,sv_i)

! Accumulate the incremental opacity derivatives to get drad_dt

      call get_dscrt_dn ( d_delta_dt(:,sv_i), t_script, tau, dt_scr_dt(:,sv_i), i_start, &
                       &  i_stop, drad_dt(sv_i) )

    end do ! sv_i

  end subroutine DRad_tran_dt

  ! ------------------------------------------------  Get_Do_Calc  -----
  subroutine Get_Do_Calc ( Do_Calc_c, Do_Calc_f, Do_GL, Do_Calc )

  ! Set Do_Calc if Do_Calc_c or Do_GL and any of the corresponding Do-Calc_f
  ! flags are set.

    use GLNP, ONLY: Ng

    logical, intent(in) :: Do_Calc_c(:) ! On the coarse grid
    logical, intent(in) :: Do_Calc_f(:) ! On the GL grid
    logical, intent(in) :: Do_GL(:)     ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)  ! Where on coarse grid to do calc.

    integer :: I, P_I

    i = 1
    do_calc = do_calc_c
    do p_i = 1 , size(do_gl)
      if ( do_gl(p_i) ) then
        if ( any(do_calc_f(i:i+ng-1)) ) do_calc(p_i)=.true.
        i = i + Ng
      end if
    end do

  end subroutine Get_Do_Calc

! =====     Private Procedures     =====================================

  ! ----------------------------------------  Get_Do_Calc_Indexed  -----
  subroutine Get_Do_Calc_Indexed ( Do_Calc_all, C_Inds, F_Inds, Do_GL, Do_Calc )

  ! Set Do_Calc if Do_Calc_c or Do_GL and any of the corresponding Do-Calc_f
  ! flags are set.

    use GLNP, ONLY: Ng
    use MLSCommon, only: IP

    logical, intent(in) :: Do_Calc_all(:) ! On the entire path
    integer(ip), intent(in) :: C_Inds(:)  ! Indices in Do_Calc_All for coarse grid
    integer(ip), intent(in) :: F_Inds(:)  ! Indices in Do_Calc_All for find grid
    logical, intent(in) :: Do_GL(:)       ! Where on coarse grid to do GL
    logical, intent(out) :: Do_Calc(:)    ! Where on coarse grid to do calc.

    integer :: I, P_I

    i = 1
    do_calc = do_calc_all(c_inds)
    do p_i = 1 , size(do_gl)
      if ( do_gl(p_i) ) then
        if ( any(do_calc_all(f_inds(i:i+ng-1))) ) do_calc(p_i)=.true.
        i = i + Ng
      end if
    end do

  end subroutine Get_Do_Calc_Indexed

  ! ---------------------------------------------------  Get_Inds  -----
  subroutine Get_Inds ( Do_GL, Do_Calc, More_Inds, All_Inds )

    ! More_Inds are the places in the coarse path where both Do_Calc and Do_GL.
    ! All_Inds are the corresponding places in the GL-extracted fine path.

    use GLNP, ONLY: Ng
    use MLSCommon, only: IP

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
!         all_inds(j:j+ng-1) = (/ ( (p_i-2)*ng+k, k = 1, ng ) /)
          all_inds(j:j+ng-1) = (/ ( l + k, k = 0, ng-1 ) /)
          i = i + 1
          j = j + Ng
        end if
        l = l + Ng
      end if
    end do

  end subroutine Get_Inds

!----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module RAD_TRAN_M

! $Log$
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
