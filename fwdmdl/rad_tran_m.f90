! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RAD_TRAN_M

  use MLSCommon, only: R8, RP, IP
  use GLNP, ONLY: Ng

  implicit NONE
  private
  public :: RAD_TRAN, DRAD_TRAN_DF, DRAD_TRAN_DT, DRAD_TRAN_DX

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

  subroutine Rad_tran ( indicies_c, gl_inds, frq, s_temp, e_rflty, z_path, &
                     &  t_path, alpha_path_c, ref_cor, do_gl, incoptdepth, &
                     &  alpha_path_gl, ds_dh_gl, dh_dz_gl, t_script, tau, &
                     &  rad, i_stop )

    use DO_DELTA_M, ONLY: PATH_OPACITY
    use DO_T_SCRIPT_M, ONLY: TWO_D_T_SCRIPT
    use SCRT_DN_M, ONLY: SCRT_DN

  ! inputs

    integer(ip), intent(in) :: indicies_c(:) ! coarse grid indicies
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indicies
    real(r8), intent(in) :: frq ! calculation frequency in MHz.
    real(rp), intent(in) :: s_temp ! farside boundary temperature
  !                                usually cosmic space (2.7K).
    real(rp), intent(in) :: e_rflty ! earth reflectivity value (0--1).
    real(rp), intent(in) :: z_path(:) ! path -log(P) on input grid.
    real(rp), intent(in) :: t_path(:) ! path T(K) on input grid.
    real(rp), intent(in) :: alpha_path_c(:) ! absorption coefficient
  !                         on input grid.
    real(rp), intent(in) :: ref_cor(:) ! refracted to unrefracted path
  !                                      length ratios.
    logical, intent(in) :: do_gl(:) ! path flag indicating where to do
  !                                   gl integrations.
    real(rp), intent(inout) :: incoptdepth(:) ! incremental path opacities
  !                            from one-sided layer calculation on output.
  !                            it is the full integrated layer opacity.
    real(rp), intent(in) :: alpha_path_gl(:) ! absorption coefficient on gl
  !                                        grid.
    real(rp), intent(in) :: ds_dh_gl(:) ! path length wrt height derivative on
  !                                        gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
  !                                       gl grid.
  ! outputs

    real(rp), intent(out) :: t_script(:) ! differential temperatures (K).
    real(rp), intent(out) :: tau(:) ! transmission function.
    real(rp), intent(out) :: rad    ! radiances (K)
    integer(ip), intent(out) :: i_stop ! path stop index

  ! Internals

    integer(ip) :: i, j, k, p_i, n_path
    real(rp) :: gl_delta(size(gl_inds)/ng), del_zeta(size(gl_inds)/ng)
    integer(ip) ::  more_inds(size(gl_inds)/ng)

  ! Begin code

    call two_d_t_script ( t_path(indicies_c), s_temp, frq, t_script )

    n_path = size(indicies_c)
    if ( count(do_gl) > 0 ) then

  ! see if anything needs to be gl-d

      i = 1
      j = 1
      do p_i = 1, n_path
        if ( do_gl(p_i) ) then
          more_inds(i) = p_i
          if ( p_i > n_path/2 ) then
            del_zeta(i) = z_path(indicies_c(p_i+1)) - z_path(indicies_c(p_i))
          else
            del_zeta(i) = z_path(indicies_c(p_i-1)) - z_path(indicies_c(p_i))
          end if
          i = i + 1
          j = j + Ng
        end if
      end do

      call path_opacity ( del_zeta, &
                &  alpha_path_c(more_inds), &
                &  alpha_path_gl, ds_dh_gl(gl_inds),  &
                &  dh_dz_gl(gl_inds), gl_delta )
      incoptdepth(more_inds) = incoptdepth(more_inds) + gl_delta

    end if

    incoptdepth = ref_cor * incoptdepth

    call scrt_dn ( t_script, e_rflty, incoptdepth, tau, rad, i_stop )

  end subroutine Rad_tran

!--------------------------------------------------  drad_tran_df  -----
! This is the radiative transfer derivative wrt mixing ratio model

  subroutine DRad_tran_df ( indicies_c, gl_inds, z_path, Grids_f, &
                         &  beta_path_c, eta_zxp_f, sps_path, do_calc_f, &
                         &  beta_path_f, do_gl, del_s, ref_cor, ds_dh_gl, &
                         &  dh_dz_gl, t_script, tau, i_stop, drad_df, ptg_i, &
                         &  frq_i)

    use DO_DELTA_M, ONLY: PATH_OPACITY
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use SCRT_DN_M, ONLY: GET_DSCRT_NO_T_DN

! Inputs

    integer(ip), intent(in) :: indicies_c(:) ! coarse grid indicies
    integer(ip), intent(in) :: gl_inds(:)    ! Gauss-Legendre grid indicies
    real(rp), intent(in) :: z_path(:) ! -log(P) on main grid.
    type (Grids_T), intent(in) :: Grids_f ! All the coordinates
    real(rp), intent(in) :: beta_path_c(:,:) ! cross section for each species
!                                              on coarse grid.
    real(rp), intent(in) :: eta_zxp_f(:,:) ! representation basis function.
    real(rp), intent(in) :: sps_path(:,:) ! Path species function.
    logical, intent(in) :: do_calc_f(:,:) ! A logical indicating where the
!                                           representation basis function is
!                                           not zero.
    real(rp), intent(in) :: beta_path_f(:,:) ! cross section for each species
!                                              on gl grid.
    logical, intent(in) :: do_gl(:) ! A logical indicating where to do gl
!                                     integrations
    real(rp), intent(in) :: ref_cor(:) ! refracted to unrefracted path
!                                        length ratios.
    real(rp), intent(in) :: del_s(:) ! unrefracted path length.
    real(rp), intent(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                          gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                         gl grid.
    real(rp), intent(in) :: t_script(:) ! differential temperatures (K).
    real(rp), intent(in) :: tau(:) ! transmission function.
    integer(ip), intent(in) :: i_stop ! path stop index
    integer(ip), intent(in) :: ptg_i,frq_i ! debugger statements

! Outputs

    real(rp), intent(out) :: drad_df(:)    ! derivative of radiances wrt
!                                            mixing ratio statevector
!                                            element. (K)
! Internals

    integer(ip) :: sv_i, sv_j, sps_i, n_inds, i, j, k, l, i_start, no_mol
    integer(ip) :: n_path, p_i, no_to_gl, mid, n_tot
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part_of_all_inds_B
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B
    integer(ip), pointer :: more_inds(:) ! more_inds => part_of_more_inds_B

    real(rp) :: d_delta_df(1:size(tau))
    real(rp), target, dimension(1:size(tau)) :: del_zeta_B, gl_delta_B
    real(rp), pointer :: del_zeta(:)     ! del_zeta => part_of_del_zeta_B
    real(rp), pointer :: gl_delta(:)     ! gl_delta => part_of_gl_delta_B
    real(rp), target :: singularity_B(1:size(tau))
    real(rp), pointer :: singularity(:)  ! singularity => part_of_singularity_B
    logical :: do_calc(1:size(tau))

! Begin code

    no_mol = size(Grids_f%no_z)
    n_path = size(tau)
    mid = n_path / 2

    sv_i = 0
    drad_df(:) = 0.0_rp

    do sps_i = 1, no_mol

      n_tot = Grids_f%no_f(sps_i)*Grids_f%no_z(sps_i)*Grids_f%no_p(sps_i)

      do sv_j = 1, n_tot

        sv_i = sv_i + 1

! Skip the masked derivatives, according to the l2cf inputs

        if ( .not. Grids_f%deriv_flags(sv_i) ) cycle

        i = 1
        d_delta_df = 0.0_rp
        do_calc = do_calc_f(indicies_c,sv_i)
        do p_i = 1 , n_path
          if ( do_gl(p_i) ) then
            if ( any(do_calc_f(gl_inds(i:i+ng-1),sv_i))) do_calc(p_i)=.true.
            i = i + Ng
          end if
        end do

! find where the non zeros are along the path

        n_inds = count(do_calc)
        if ( n_inds > 0 ) then

          inds => inds_B(1:n_inds)
          singularity => singularity_B(1:n_inds)

          i = 1
          do p_i = 1 , n_path
            if ( do_calc(p_i) ) then
              inds(i) = p_i
              i = i + 1
            end if
          end do

          if ( grids_f%lin_log(sps_i) ) then

            singularity = beta_path_c(inds,sps_i) &
                      & * eta_zxp_f(indicies_c(inds),sv_i) &
                      & * sps_path(indicies_c(inds),sps_i)
            d_delta_df(inds) = singularity * del_s(inds)

            no_to_gl = count(do_gl(inds))
            if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

              all_inds => all_inds_B(1:ng*no_to_gl)
              gl_delta => gl_delta_B(1:no_to_gl)
              del_zeta => del_zeta_B(1:no_to_gl)
              more_inds => more_inds_B(1:no_to_gl)

              i = 1
              j = 1
              l = 1
              do p_i = 1 , n_path
                if ( do_gl(p_i) ) then
                  if ( do_calc(p_i) ) then
                    more_inds(i) = p_i
                    all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                    if ( p_i > mid ) then
                      del_zeta(i) = z_path(indicies_c(p_i+1)) &
                                & - z_path(indicies_c(p_i))
                    else
                      del_zeta(i) = z_path(indicies_c(p_i-1)) &
                                & - z_path(indicies_c(p_i))
                    end if
                    i = i + 1
                    j = j + Ng
                  end if
                  l = l + Ng
                end if
              end do

              call path_opacity ( del_zeta, &
                   pack(singularity,do_gl(inds)),     &
                   beta_path_f(all_inds,sps_i) &
                 & * eta_zxp_f(gl_inds(all_inds),sv_i) &
                 & * sps_path(gl_inds(all_inds),sps_i), &
                 &   ds_dh_gl(gl_inds(all_inds)), &
                 &   dh_dz_gl(gl_inds(all_inds)), gl_delta )
              d_delta_df(more_inds) = d_delta_df(more_inds) + gl_delta

            end if

            d_delta_df(inds) = ref_cor(inds)*d_delta_df(inds) &
                             / exp(grids_f%values(sv_i))

          else

            singularity = beta_path_c(inds,sps_i) &
                      & * eta_zxp_f(indicies_c(inds),sv_i)
            d_delta_df(inds) = singularity * del_s(inds)

            no_to_gl = count(do_gl(inds))
            if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

              all_inds => all_inds_B(1:ng*no_to_gl)
              gl_delta => gl_delta_B(1:no_to_gl)
              del_zeta => del_zeta_B(1:no_to_gl)
              more_inds => more_inds_B(1:no_to_gl)

              i = 1
              j = 1
              l = 1
              do p_i = 1 , n_path
                if ( do_gl(p_i) ) then
                  if ( do_calc(p_i) ) then
                    more_inds(i) = p_i
                    all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                    if ( p_i > mid ) then
                      del_zeta(i) = z_path(indicies_c(p_i+1)) &
                                & - z_path(indicies_c(p_i))
                    else
                      del_zeta(i) = z_path(indicies_c(p_i-1)) &
                                & - z_path(indicies_c(p_i))
                    end if
                    i = i + 1
                    j = j + Ng
                  end if
                  l = l + Ng
                end if
              end do

              call path_opacity(del_zeta, &
                   pack(singularity,do_gl(inds)), &
                   beta_path_f(all_inds,sps_i) &
                 & * eta_zxp_f(gl_inds(all_inds),sv_i), &
                 & ds_dh_gl(gl_inds(all_inds)),dh_dz_gl(gl_inds(all_inds)), &
                 & gl_delta)
              d_delta_df(more_inds) = d_delta_df(more_inds) + gl_delta

            end if

            d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds)

          end if

          i_start = MIN(inds(1),i_stop)

          call get_dscrt_no_t_dn ( d_delta_df, t_script, tau, i_start, i_stop, &
                                &  drad_df(sv_i))

        end if

      end do

    end do

  end subroutine drad_tran_df

!--------------------------------------------------  drad_tran_dx  -----
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )

  subroutine DRad_tran_dx ( z_path_c, Grids_f, dbeta_path_c, eta_zxp_f_c, &
                         &  sps_path_c, do_calc_f_c, dbeta_path_f, eta_zxp_f_f, &
                         &  sps_path_f, do_calc_f_f, do_gl, del_s, ref_cor,     &
                         &  ds_dh_gl, dh_dz_gl, t_script, tau, i_stop, drad_dx, &
                         &  ptg_i, frq_i )

    use DO_DELTA_M, ONLY: PATH_OPACITY
    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use SCRT_DN_M, ONLY: GET_DSCRT_NO_T_DN

! Inputs

    real(rp), intent(in) :: z_path_c(:) ! -log(P) on main grid.
    type (Grids_T) :: Grids_f           ! All the coordinates
    real(rp), intent(in) :: dbeta_path_c(:,:) ! derivative of beta wrt dx
!                                              on main grid.
    real(rp), intent(in) :: eta_zxp_f_c(:,:) ! representation basis function
!                                              main grid.
    real(rp), intent(in) :: sps_path_c(:,:) ! species function on  main grid.
    logical, intent(in) :: do_calc_f_c(:,:) ! A logical indicating where the
!                                           representation basis function is
!                                           not zero on main grid.
    real(rp), intent(in) :: dbeta_path_f(:,:) ! derivative of beta wrt dx
!                                              on gl grid.
    real(rp), intent(in) :: eta_zxp_f_f(:,:) ! representation basis function
!                                              gl grid.
    real(rp), intent(in) :: sps_path_f(:,:) ! species function on gl grid.
    logical, intent(in) :: do_calc_f_f(:,:) ! A logical indicating where the
!                                           representation basis function is
!                                           not zero on main grid.
    logical, intent(in) :: do_gl(:) ! A logical indicating where to do gl
!                                     integrations
    real(rp), intent(in) :: ref_cor(:) ! refracted to unrefracted path
!                                        length ratios.
    real(rp), intent(in) :: del_s(:) ! unrefracted path length.
    real(rp), intent(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                          gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                         gl grid.
    real(rp), intent(in) :: t_script(:) ! differential temperatures (K).
    real(rp), intent(in) :: tau(:) ! transmission function.
    integer(ip), intent(in) :: i_stop ! path stop index
    integer(ip), intent(in) :: ptg_i,frq_i ! debugger statements

! Outputs

    real(rp), intent(out) :: drad_dx(:)    ! derivative of radiances wrt x
!                                            statevector element. (K)
! Internals

    integer(ip) :: sv_i, sv_j, sps_i, n_inds, i, j, k, l, i_start, no_mol
    integer(ip) :: n_path, p_i, no_to_gl, mid
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part_of_all_inds_B
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B
    integer(ip), pointer :: more_inds(:) ! more_inds => part_of_more_inds_B

    real(rp) :: d_delta_dx(1:size(tau))
    real(rp), target, dimension(1:size(tau)) :: del_zeta_B, gl_delta_B
    real(rp), pointer :: del_zeta(:)     ! del_zeta => part_of_del_zeta_B
    real(rp), pointer :: gl_delta(:)     ! gl_delta => part_of_gl_delta_B
    real(rp), target :: singularity_B(1:size(tau))
    real(rp), pointer :: singularity(:)  ! singularity => part_of_singularity_B

    logical :: do_calc(1:size(tau))

! Begin code

    no_mol = size(grids_f%no_z)
    n_path = size(tau)
    mid = n_path / 2

    sv_i = 0
    drad_dx(:) = 0.0_rp

    do sps_i = 1 , no_mol

      do sv_j = 1 , Grids_f%no_z(sps_i) * Grids_f%no_p(sps_i)

        sv_i = sv_i + 1
        d_delta_dx = 0.0_rp

        do_calc = do_calc_f_c(:,sv_i)

        i = 1
        do p_i = 1 , n_path
          if ( do_gl(p_i) ) then
            if ( any(do_calc_f_f(i:i+Ng-1,sv_i))) do_calc(p_i) = .true.
            i = i + Ng
          end if
        end do

! find where the non zeros are along the path

        n_inds = count(do_calc)
        if ( n_inds > 0 ) then

          inds => inds_B(1:n_inds)
          singularity => singularity_B(1:n_inds)

          i = 1
          do p_i = 1 , n_path
            if ( do_calc(p_i) ) then
              inds(i) = p_i
              i = i + 1
            end if
          end do

          singularity = dbeta_path_c(inds,sps_i) * eta_zxp_f_c(inds,sv_i)  &
                     &  * sps_path_c(inds,sps_i)
          d_delta_dx(inds) = singularity * del_s(inds)

          no_to_gl = count(do_gl(inds))
          if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

            all_inds => all_inds_B(1:Ng*no_to_gl)
            more_inds => more_inds_B(1:no_to_gl)
            gl_delta => gl_delta_B(1:no_to_gl)
            del_zeta => del_zeta_B(1:no_to_gl)

            i = 1
            j = 1
            l = 1
            do p_i = 1 , n_path
              if ( do_gl(p_i) ) then
                if ( do_calc(p_i) ) then
                  more_inds(i) = p_i
                  all_inds(j:j+Ng-1) = l + (/(k-1,k=1,Ng)/)
                  if ( p_i > mid ) then
                    del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                  else
                    del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                  end if
                  i = i + 1
                  j = j + Ng
                end if
                l = l + Ng
              end if
            end do

            call path_opacity ( del_zeta, &
                &  pack(singularity,do_gl(inds)),        &
                &  dbeta_path_f(all_inds,sps_i)*eta_zxp_f_f(all_inds,sv_i) * &
                &    sps_path_f(all_inds,sps_i), ds_dh_gl(all_inds),         &
                &  dh_dz_gl(all_inds), gl_delta )
            d_delta_dx(more_inds) = d_delta_dx(more_inds) + gl_delta

          end if

          d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)

          i_start = min(inds(1),i_stop)

          call get_dscrt_no_t_dn ( d_delta_dx, t_script, tau, i_start, i_stop, &
                                &  drad_dx(sv_i) )

        end if

      end do

    end do

  end subroutine DRad_tran_dx
!--------------------------------------------------  drad_tran_dt  -----
! This is the radiative transfer derivative wrt temperature model

  subroutine DRad_tran_dt ( z_path_c, h_path_c, t_path_c, dh_dt_path_c, &
                         &  alpha_path_c, alphaxn_path_c, eta_zxp_c, do_calc_t_c, &
                         &  do_calc_hyd_c, del_s, ref_cor, h_tan, dh_dt_tan, &
                         &  freq, do_gl, h_path_f, t_path_f, dh_dt_path_f, &
                         &  alpha_path_f, alphaxn_path_f, eta_zxp_f, do_calc_t_f, &
                         &  ds_dh_gl, dh_dz_gl, t_script, tau, i_stop, drad_dt,  &
                         &  ptg_i, frq_i )

    use DO_DELTA_M, ONLY: PATH_OPACITY, HYD_OPACITY
    use D_T_SCRIPT_DTNP_M, only: DT_SCRIPT_DT
    use SCRT_DN_M, ONLY: GET_DSCRT_DN

! Inputs

    real(rp), intent(in) :: z_path_c(:) ! path -log(P) on main grid.
    real(rp), intent(in) :: h_path_c(:) ! path heights + req on main grid km.
    real(rp), intent(in) :: t_path_c(:) ! path temperature(K) on main grid.
    real(rp), intent(in) :: dh_dt_path_c(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on main grid.
    real(rp), intent(in) :: alpha_path_c(:) ! path absorption(km^-1)
!                                             on main grid.
    real(rp), intent(in) :: alphaxn_path_c(:) ! path absorption times
!                            temperature power(km^-1) on main grid.
    real(rp), intent(in) :: eta_zxp_c(:,:) ! representation basis function
!                                              main grid.
    logical, intent(in) :: do_calc_t_c(:,:) ! A logical indicating where the
!                    representation basis function is not zero on main grid.
    logical, intent(in) :: do_calc_hyd_c(:,:) ! A logical indicating where the
!                    dh_dt function is not zero on main grid.
    real(rp), intent(in) :: del_s(:) ! unrefracted path length.
    real(rp), intent(in) :: ref_cor(:) ! refracted to unrefracted path
!                                        length ratios.
    real(rp), intent(in) :: h_tan ! tangent height + req (km).
    real(rp), intent(in) :: dh_dt_tan(:) ! derivative of path height wrt
!                                       temperature at the tangent (km/K).
    real(r8), intent(in) :: freq ! calculation frequency (MHz).
    logical, intent(in) :: do_gl(:) ! A logical indicating where to do gl
!                                     integrations.
    real(rp), intent(in) :: h_path_f(:) ! path heights + req on gl grid km.
    real(rp), intent(in) :: t_path_f(:) ! path temperature(K) on gl grid.
    real(rp), intent(in) :: dh_dt_path_f(:,:) ! derivative of path height wrt
!                                               temperature(km/K) on gl grid.
    real(rp), intent(in) :: alpha_path_f(:) ! path absorption(km^-1) on gl grid.
    real(rp), intent(in) :: alphaxn_path_f(:) ! path absorption times
!                            temperature power(km^-1) on gl grid.
    real(rp), intent(in) :: eta_zxp_f(:,:) ! representation basis function
!                                              gl grid.
    logical, intent(in) :: do_calc_t_f(:,:) ! A logical indicating where the
!                    representation basis function is not zero on gl grid.
    real(rp), intent(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                          gl grid.
    real(rp), intent(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                         gl grid.
    real(rp), intent(in) :: t_script(:) ! differential temperatures (K).
    real(rp), intent(in) :: tau(:) ! transmission function.
    integer(ip), intent(in) :: i_stop ! path stop index
    integer(ip), intent(in) :: ptg_i,frq_i ! debugger statements

! Outputs

    real(rp), intent(out) :: drad_dt(:)    ! derivative of radiances wrt
!                                            mixing ratio statevector
!                                            element. (K)
! Internals

    integer(ip) :: sv_i, sv_t, n_inds, i, j, k, l, i_start
    integer(ip) :: n_path, p_i, no_to_gl, mid
    integer(ip), target, dimension(1:Ng*size(tau)) :: all_inds_B
    integer(ip), target, dimension(1:size(tau)) :: inds_B, more_inds_B
    integer(ip), pointer :: all_inds(:)  ! all_inds => part_of_all_inds_B
    integer(ip), pointer :: inds(:)      ! inds => part_of_inds_B
    integer(ip), pointer :: more_inds(:) ! more_inds => part_of_more_inds_B

    real(rp) :: d_delta_dt(1:size(tau))
    real(rp) :: dt_scr_dt(1:size(tau),1:size(eta_zxp_c,dim=2))
    real(rp), target, dimension(1:size(tau)) :: del_zeta_B, gl_delta_B
    real(rp), pointer :: del_zeta(:)     ! del_zeta => part_of_del_zeta_B
    real(rp), pointer :: gl_delta(:)     ! gl_delta => part_of_gl_delta_B
    real(rp), target :: singularity_B(1:size(tau))
    real(rp), pointer :: singularity(:)  ! singularity => part_of_singularity_B
    real(rp) :: fa, fb

    logical :: do_calc(1:size(tau))

! Begin code

    n_path = size(tau)
    sv_t = size(eta_zxp_c,dim=2)
    mid = n_path / 2

! compute the t_script derivative

    call dt_script_dt ( t_path_c, eta_zxp_c, freq, dt_scr_dt )

! compute the opacity derivative singularity value

    drad_dt(:) = 0.0_rp

    do sv_i = 1 , sv_t

      i_start = 1
      d_delta_dt = 0.0_rp

! do the absorption part
! combine non zeros flags for both the main and gl parts

      do_calc = do_calc_t_c(:,sv_i)

      i = 1
      do p_i = 1 , n_path
        if ( do_gl(p_i) ) then
          if ( any(do_calc_t_f(i:i+ng-1,sv_i))) do_calc(p_i) = .true.
          i = i + Ng
        end if
      end do

! find where the non zeros are along the path

      n_inds = count(do_calc)
      if ( n_inds > 0 ) then
        inds => inds_B(1:n_inds)
        singularity => singularity_B(1:n_inds)
        i = 1
        do p_i = 1 , n_path
          if ( do_calc(p_i) ) then
            inds(i) = p_i
            i = i + 1
          end if
        end do

        singularity = alphaxn_path_c(inds)*eta_zxp_c(inds,sv_i) &
                    / t_path_c(inds)
        d_delta_dt(inds) = singularity * del_s(inds)

        no_to_gl = count(do_gl(inds))
        if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

           all_inds => all_inds_B(1:ng*no_to_gl)
           gl_delta => gl_delta_B(1:no_to_gl)
           del_zeta => del_zeta_B(1:no_to_gl)
           more_inds => more_inds_B(1:no_to_gl)

           i = 1
           j = 1
           l = 1
           do p_i = 1 , n_path
             if ( do_gl(p_i) ) then
               if ( do_calc(p_i) ) then
                 more_inds(i) = p_i
                 all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                 if ( p_i > mid ) then
                   del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                 else
                   del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                 end if
                 i = i + 1
                 j = j + Ng
               end if
               l = l + Ng
             end if
           end do

           call path_opacity ( del_zeta, &
             &  pack(singularity,do_gl(inds)), &
             &  alphaxn_path_f(all_inds)*eta_zxp_f(all_inds,sv_i)  &
             &    / t_path_f(all_inds), ds_dh_gl(all_inds), dh_dz_gl(all_inds), &
             &  gl_delta )
           d_delta_dt(more_inds) = d_delta_dt(more_inds) + gl_delta

        end if

        i_start = max(inds(1)-1,1)

      end if

! now do the hydrostatic part
! combine boundaries flags

      do_calc = do_calc_hyd_c(:,sv_i)
      if ( do_calc_hyd_c(mid,sv_i) ) then
        do_calc(2:mid) = .true.
      else
        do p_i = 2 , mid
          if ( do_calc_hyd_c(p_i-1,sv_i)) do_calc(p_i) = .true.
        end do
      end if
      if ( do_calc_hyd_c(mid+1,sv_i) ) then
        do_calc(mid+1:n_path-1) = .true.
      else
        do p_i = mid + 1 , n_path - 1
          if ( do_calc_hyd_c(p_i+1,sv_i)) do_calc(p_i) = .true.
        end do
      end if

! since this is a layer boundary calculation we must require

      do_calc((/1,n_path/)) = .false.

! find where the non zeros are along the path

      n_inds = count(do_calc)
      if ( n_inds > 0 ) then
        inds => inds_B(1:n_inds)
        i = 1
        inds = 0
        do p_i = 2 , mid - 1
          if ( do_calc(p_i) ) then
            if ( i == 1) fa = (h_path_c(p_i-1)*dh_dt_path_c(p_i-1,sv_i) &
                          - h_tan * dh_dt_tan(sv_i)) / SUM(del_s(p_i:mid))
            fb = (h_path_c(p_i)*dh_dt_path_c(p_i,sv_i) &
               - h_tan * dh_dt_tan(sv_i)) /SUM(del_s(p_i+1:mid))
            inds(i) = p_i
            d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i)*(fa - fb)
            fa = fb
            i = i + 1
          end if
        end do

! special processing at tangent

        if ( do_calc(mid) ) then
          d_delta_dt(mid) = d_delta_dt(mid) + alpha_path_c(mid) * fa

! fb is zero

          inds(i) = mid
          i = i + 1
        end if

        if ( do_calc(mid+1) ) then
          fa = (h_path_c(mid+2)*dh_dt_path_c(mid+2,sv_i) - &
              & h_tan * dh_dt_tan(sv_i))/ del_s(mid+1)
          d_delta_dt(mid+1) = d_delta_dt(mid+1) + alpha_path_c(mid+1) * fa

! fb is 0.0

          inds(i) = mid+1
          i = i + 1
        end if

        do p_i = mid + 2 , n_path - 1
          if ( do_calc(p_i) ) then
            if ( inds(max(i-1,1)) < mid+1)  &
                       & fa = (h_path_c(p_i)*dh_dt_path_c(p_i,sv_i) &
                       &    - h_tan * dh_dt_tan(sv_i)) / SUM(del_s(mid+1:p_i-1))
            fb = (h_path_c(p_i+1)*dh_dt_path_c(p_i+1,sv_i) &
               - h_tan * dh_dt_tan(sv_i)) /sum(del_s(mid+1:p_i))
            inds(i) = p_i
            d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i)*(fb - fa)
            fa = fb
            i = i + 1
          end if
        end do

        no_to_gl = count(do_gl(inds))
        if ( no_to_gl > 0 ) then

! see if anything needs to be gl-d

           all_inds => all_inds_B(1:ng*no_to_gl)
           gl_delta => gl_delta_B(1:no_to_gl)
           del_zeta => del_zeta_B(1:no_to_gl)
           more_inds => more_inds_B(1:no_to_gl)

           i = 1
           j = 1
           l = 1
           do p_i = 1 , n_path
             if ( do_gl(p_i) ) then
               if ( do_calc(p_i) ) then
                 more_inds(i) = p_i
                 all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                 if ( p_i > mid ) then
                   del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                 else
                   del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                 end if
                 i = i + 1
                 j = j + Ng
               end if
               l = l + Ng
             end if
           end do

! add special hydrostatic gl routine here
! the singularity point is alpha_path_c(more_inds)

           call hyd_opacity ( del_zeta, alpha_path_c(more_inds),                &
              & alpha_path_f(all_inds), h_path_f(all_inds),                    &
              & dh_dt_path_f(all_inds,sv_i), t_path_f(all_inds), h_tan,        &
              & dh_dt_tan(sv_i), eta_zxp_f(all_inds,sv_i), ds_dh_gl(all_inds), &
              & dh_dz_gl(all_inds), gl_delta )

           d_delta_dt(more_inds) = d_delta_dt(more_inds) + gl_delta

        end if

        i_start = min(i_start,inds(1))

      end if

! correct the whole thing for path length refraction

      d_delta_dt = ref_cor * d_delta_dt
      call get_dscrt_dn ( d_delta_dt, t_script, tau, dt_scr_dt(:,sv_i), i_start, &
                       &  i_stop, drad_dt(sv_i) )

    end do

  end subroutine DRad_tran_dt

!----------------------------------------------------------------------
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module RAD_TRAN_M
! $Log$
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
