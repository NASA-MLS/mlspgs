module RAD_TRAN_M
!
  use MLSCommon, only: R8, RP, IP
  USE GLNP, ONLY: Ng
  USE DO_T_SCRIPT_M, ONLY: TWO_D_T_SCRIPT
  USE D_T_SCRIPT_DTNP_M, ONLY: DT_SCRIPT_DT
  USE DO_DELTA_M, ONLY: PATH_OPACITY,HYD_OPACITY
  USE SCRT_DN_M, ONLY: SCRT_DN, GET_DSCRT_NO_T_DN, GET_DSCRT_DN
  USE LOAD_SPS_DATA_M, ONLY: GRIDS_T
  implicit NONE
  private
  PUBLIC :: PATH_CONTRIB, RAD_TRAN, DRAD_TRAN_DF, DRAD_TRAN_DT, DRAD_TRAN_DX
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This routine computes the contributions (along the path) of each interval
! of the (coarse) pre-selected integration grid

SUBROUTINE path_contrib(alpha_path,del_s,e_rflty,tol,dtaudn, &
                    &   incoptdepth,do_gl)
!
! inputs
!
  REAL(rp), INTENT(in) :: alpha_path(:) ! absorption coefficient along the
!                                         path
  REAL(rp), INTENT(in) :: del_s(:)       ! path lengths
  REAL(rp), INTENT(in) :: tol            ! accuracy target in K
  REAL(rp), INTENT(in) :: e_rflty        ! earth reflectivity
!
! outputs
!
  REAL(rp), INTENT(out) :: dtaudn(:) ! path derivative of the transmission
!                                    function
  REAL(rp), INTENT(out) :: incoptdepth(:) ! layer optical depth
  LOGICAL(ip), INTENT(out) :: do_gl(:) ! set true for indicies to do
!                                         gl computation
! Internal stuff
!
  INTEGER(ip) :: i, i_tan, n_path
  REAL(rp), PARAMETER :: temp = 250.0_rp
!
! start code
! compute the incoptdepth
!
  n_path = SIZE(alpha_path)
  i_tan = n_path / 2
!
  incoptdepth = alpha_path * del_s
!
  dtaudn(1) = 0.0_rp
  DO i = 2 , i_tan
    dtaudn(i) = dtaudn(i-1) - incoptdepth(i)
  ENDDO
!
  dtaudn(i_tan+1) = dtaudn(i_tan)
!
  DO i = i_tan+2,n_path
    dtaudn(i) = dtaudn(i-1) - incoptdepth(i-1)
  ENDDO
!
! compute the tau path derivative
!
  dtaudn = 0.5_rp*(EOSHIFT(dtaudn,1,dtaudn(n_path)) -             &
                   EOSHIFT(dtaudn,-1,dtaudn(1))) * EXP(dtaudn) *  &
           (/(1.0_rp,i=1,i_tan),(e_rflty,i=i_tan+1,n_path)/)
!
! find the indicies
!
  do_gl = .FALSE.
!
  WHERE(dtaudn < -tol/temp) do_gl = .TRUE.
!
! The first and last index must be false
!
  do_gl((/1,n_path/)) = .FALSE.
!
END SUBROUTINE path_contrib
!
!----------------------------------------------------------------------
! This is the radiative transfer model, radiances only !
!
SUBROUTINE rad_tran(frq,s_temp,e_rflty,z_path_c,t_path_c,alpha_path_c, &
                 &  ref_cor,do_gl,incoptdepth,alpha_path_gl,ds_dh_gl,  &
                 &  dh_dz_gl,t_script,tau,rad,i_stop)
!
! inputs
!
  REAL(r8), INTENT(in) :: frq ! calculation frequency in MHz.
  REAL(rp), INTENT(in) :: s_temp ! farside boundary temperature
!                                usually cosmic space (2.7K).
  REAL(rp), INTENT(in) :: e_rflty ! earth reflectivity value (0--1).
  REAL(rp), INTENT(in) :: z_path_c(:) ! path -log(P) on input grid.
  REAL(rp), INTENT(in) :: t_path_c(:) ! path T(K) on input grid.
  REAL(rp), INTENT(in) :: alpha_path_c(:) ! absorption coefficient
!                         on input grid.
  REAL(rp), INTENT(in) :: ref_cor(:) ! refracted to unrefracted path
!                                      length ratios.
  LOGICAL, INTENT(in) :: do_gl(:) ! path flag indicating where to do
!                                   gl integrations.
  REAL(rp), INTENT(inout) :: incoptdepth(:) ! incremental path opacities
!                            from one-sided layer calculation on output.
!                            it is the full integrated layer opacity.
  REAL(rp), INTENT(in) ::  alpha_path_gl(:) ! absorption coefficient on gl
!                                        grid.
  REAL(rp), INTENT(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                        gl grid.
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                       gl grid.
! out puts
!
  REAL(rp), INTENT(out) :: t_script(:) ! differential temperatures (K).
  REAL(rp), INTENT(out) :: tau(:) ! transmission function.
  REAL(rp), INTENT(out) :: rad    ! radiances (K)
  INTEGER(ip), INTENT(out) :: i_stop ! path stop index
!
! Internals
!
  INTEGER(ip) :: i,j,k,p_i,no_to_gl,n_path
  REAL(rp), ALLOCATABLE :: gl_delta(:),del_zeta(:)
  INTEGER(ip), ALLOCATABLE :: all_inds(:),more_inds(:)
!
! Begin code
!
  CALL two_d_t_script(t_path_c,s_temp,frq,t_script)
  n_path = SIZE(z_path_c)
  no_to_gl = COUNT(do_gl)
  IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
    ALLOCATE(all_inds(1:ng*no_to_gl))
    ALLOCATE(gl_delta(1:no_to_gl))
    ALLOCATE(del_zeta(1:no_to_gl))
    ALLOCATE(more_inds(1:no_to_gl))

    i = 1
    j = 1
    DO p_i = 1 , n_path
      IF(do_gl(p_i)) THEN
        more_inds(i) = p_i
        all_inds(j:j+ng-1) = j + (/(k-1,k=1,ng)/)
        IF(p_i > n_path/2) THEN
          del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
        ELSE
          del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
        ENDIF
        i = i + 1
        j = j + Ng
      ENDIF
    ENDDO
!
    CALL path_opacity(del_zeta, &
              &  alpha_path_c(PACK((/(i,i=1,n_path)/),do_gl)), &
              &  alpha_path_gl(all_inds),ds_dh_gl(all_inds), &
              &  dh_dz_gl(all_inds),gl_delta)
    incoptdepth(more_inds) = incoptdepth(more_inds) + gl_delta

    DEALLOCATE(more_inds)
    DEALLOCATE(del_zeta)
    DEALLOCATE(gl_delta)
    DEALLOCATE(all_inds)
!
  ENDIF
!
  incoptdepth = ref_cor * incoptdepth
!
  CALL scrt_dn(t_script,e_rflty,incoptdepth,tau,rad,i_stop)
!
END SUBROUTINE rad_tran
!
!----------------------------------------------------------------------
! This is the radiative transfer derivative wrt mixing ratio model
!
SUBROUTINE drad_tran_df(z_path_c,Grids_f,lin_log,sps_values,            &
                     &  beta_path_c,eta_zxp_f_c,sps_path_c,do_calc_f_c, &
                     &  beta_path_f,eta_zxp_f_f,sps_path_f,do_calc_f_f, &
                     &  do_gl,del_s,ref_cor,ds_dh_gl,dh_dz_gl,t_script,tau, &
                     &  i_stop,drad_df,ptg_i,frq_i)
!
! Inputs
!
  REAL(rp), INTENT(in) :: z_path_c(:) ! -log(P) on main grid.
  Type (Grids_T) :: Grids_f           ! All the coordinates
  LOGICAL, INTENT(in) :: lin_log(:) ! kind of basis to use for each species.
  REAL(rp), INTENT(in) :: sps_values(:) ! sps basis break-point values.
  REAL(rp), INTENT(in) :: beta_path_c(:,:) ! cross section for each species
!                                            on main grid.
  REAL(rp), INTENT(in) :: eta_zxp_f_c(:,:) ! representation basis function
!                                            main grid.
  REAL(rp), INTENT(in) :: sps_path_c(:,:) ! species function on  main grid.
  LOGICAL, INTENT(in) :: do_calc_f_c(:,:) ! A logical indicating where the
!                                         representation basis function is
!                                         not zero on main grid.
  REAL(rp), INTENT(in) :: beta_path_f(:,:) ! cross section for each species
!                                            on gl grid.
  REAL(rp), INTENT(in) :: eta_zxp_f_f(:,:) ! representation basis function
!                                            gl grid.
  REAL(rp), INTENT(in) :: sps_path_f(:,:) ! species function on gl grid.
  LOGICAL, INTENT(in) :: do_calc_f_f(:,:) ! A logical indicating where the
!                                         representation basis function is
!                                         not zero on main grid.
  LOGICAL, INTENT(in) :: do_gl(:) ! A logical indicating where to do gl
!                                   integrations
  REAL(rp), INTENT(in) :: ref_cor(:) ! refracted to unrefracted path
!                                      length ratios.
  REAL(rp), INTENT(in) :: del_s(:) ! unrefracted path length.
  REAL(rp), INTENT(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                        gl grid.
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                       gl grid.
  REAL(rp), INTENT(in) :: t_script(:) ! differential temperatures (K).
  REAL(rp), INTENT(in) :: tau(:) ! transmission function.
  INTEGER(ip), INTENT(in) :: i_stop ! path stop index
  INTEGER(ip), INTENT(in) :: ptg_i,frq_i ! debugger statements
!
! Outputs
!
  REAL(rp), INTENT(out) :: drad_df(:)    ! derivative of radiances wrt
!                                          mixing ratio statevector
!                                          element. (K)
! Internals
!
  INTEGER(ip) :: sv_i, sv_j, sps_i, n_inds, i, j, k, l, i_start, n_sps
  INTEGER(ip) :: n_path, p_i, no_to_gl, mid, n_tot
  INTEGER(ip), ALLOCATABLE :: inds(:),all_inds(:),more_inds(:)
!
! use automatic allocations here
!
  REAL(rp), ALLOCATABLE :: d_delta_df(:),singularity(:),del_zeta(:),gl_delta(:)
  LOGICAL, ALLOCATABLE :: do_calc(:)
!
! Begin code
!
  n_sps = SIZE(Grids_f%no_z)
  n_path = SIZE(tau)
  mid = n_path / 2
!
  ALLOCATE(do_calc(1:n_path))
  ALLOCATE(d_delta_df(1:n_path))
!
  sv_i = 0
  drad_df(:) = 0.0_rp
!
  DO sps_i = 1 , n_sps
!
    n_tot = Grids_f%no_f(sps_i)*Grids_f%no_z(sps_i)*Grids_f%no_p(sps_i)
    DO sv_j = 1 , n_tot
!
      sv_i = sv_i + 1
      if(.NOT. Grids_f%deriv_flags(sv_i)) CYCLE

      i = 1
      d_delta_df = 0.0_rp
      do_calc = do_calc_f_c(:,sv_i)
      DO p_i = 1 , n_path
        IF(do_gl(p_i)) THEN
          IF(ANY(do_calc_f_f(i:i+ng-1,sv_i))) do_calc(p_i) = .TRUE.
          i = i + Ng
        ENDIF
      ENDDO
!
! find where the non zeros are along the path
!
      n_inds = COUNT(do_calc)
      IF(n_inds > 0) THEN
!
        ALLOCATE(inds(1:n_inds))
        ALLOCATE(singularity(1:n_inds))

        i = 1
        DO p_i = 1 , n_path
          IF(do_calc(p_i)) THEN
            inds(i) = p_i
            i = i + 1
          ENDIF
        ENDDO
!
        IF(lin_log(sps_i)) THEN
!
          singularity = beta_path_c(inds,sps_i)*eta_zxp_f_c(inds,sv_i) &
                     &  * sps_path_c(inds,sps_i)
          d_delta_df(inds) = singularity * del_s(inds)
!
          no_to_gl = COUNT(do_gl(inds))
          IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
            ALLOCATE(all_inds(1:ng*no_to_gl))
            ALLOCATE(gl_delta(1:no_to_gl))
            ALLOCATE(del_zeta(1:no_to_gl))
            ALLOCATE(more_inds(1:no_to_gl))

            i = 1
            j = 1
            l = 1
            DO p_i = 1 , n_path
              IF(do_gl(p_i)) THEN
                IF(do_calc(p_i)) THEN
                  more_inds(i) = p_i
                  all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                  IF(p_i > mid) THEN
                    del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                  ELSE
                    del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                  ENDIF
                  i = i + 1
                  j = j + Ng
                ENDIF
                l = l + Ng
              ENDIF
            ENDDO

            CALL path_opacity(del_zeta, &
                 singularity(PACK((/(i,i=1,n_inds)/),do_gl(inds))), &
                 beta_path_f(all_inds,sps_i)*eta_zxp_f_f(all_inds,sv_i) &
                 * sps_path_f(all_inds,sps_i),ds_dh_gl(all_inds), &
                 dh_dz_gl(all_inds),gl_delta)
            d_delta_df(more_inds) = d_delta_df(more_inds) + gl_delta

            DEALLOCATE(more_inds)
            DEALLOCATE(del_zeta)
            DEALLOCATE(gl_delta)
            DEALLOCATE(all_inds)

          ENDIF
!
          d_delta_df(inds) = ref_cor(inds)*d_delta_df(inds) &
                           / EXP(sps_values(sv_i))
!
        ELSE
!
          singularity = beta_path_c(inds,sps_i) * eta_zxp_f_c(inds,sv_i)
          DO p_i = 1 , n_inds
            d_delta_df(inds(p_i)) = singularity(p_i) * del_s(inds(p_i))
          ENDDO
!
          no_to_gl = COUNT(do_gl(inds))
          IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
            ALLOCATE(all_inds(1:ng*no_to_gl))
            ALLOCATE(gl_delta(1:no_to_gl))
            ALLOCATE(del_zeta(1:no_to_gl))
            ALLOCATE(more_inds(1:no_to_gl))

            i = 1
            j = 1
            l = 1
            DO p_i = 1 , n_path
              IF(do_gl(p_i)) THEN
                IF(do_calc(p_i)) THEN
                  more_inds(i) = p_i
                  all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
                  IF(p_i > mid) THEN
                    del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                  ELSE
                    del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                  ENDIF
                  i = i + 1
                  j = j + Ng
                ENDIF
                l = l + Ng
              ENDIF
            ENDDO

            CALL path_opacity(del_zeta, &
                 singularity(PACK((/(i,i=1,n_inds)/),do_gl(inds))), &
                 beta_path_f(all_inds,sps_i)*eta_zxp_f_f(all_inds,sv_i), &
                 ds_dh_gl(all_inds),dh_dz_gl(all_inds),gl_delta)
            d_delta_df(more_inds) = d_delta_df(more_inds) + gl_delta

            DEALLOCATE(all_inds)
            DEALLOCATE(gl_delta)
            DEALLOCATE(del_zeta)
            DEALLOCATE(more_inds)

          ENDIF
!
          d_delta_df(inds) = ref_cor(inds) * d_delta_df(inds)
!
        ENDIF
!
        i_start = MIN(inds(1),i_stop)

        DEALLOCATE(inds)
        DEALLOCATE(singularity)

        CALL get_dscrt_no_t_dn(d_delta_df,t_script,tau,i_start,i_stop, &
                            &  drad_df(sv_i))
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
  DEALLOCATE(do_calc)
  DEALLOCATE(d_delta_df)
!
END SUBROUTINE drad_tran_df
!
!----------------------------------------------------------------------
! This is the radiative transfer derivative wrt spectroscopy model
!  (Here dx could be: dw, dn or dv (dNu0) )
!
SUBROUTINE drad_tran_dx(z_path_c,Grids_f,dbeta_path_c,eta_zxp_f_c, &
                     &  sps_path_c,do_calc_f_c,dbeta_path_f,eta_zxp_f_f, &
                     &  sps_path_f,do_calc_f_f,do_gl,del_s,ref_cor,      &
                     &  ds_dh_gl,dh_dz_gl,t_script,tau,i_stop,drad_dx,   &
                     &  ptg_i,frq_i)
!
! Inputs
!
  REAL(rp), INTENT(in) :: z_path_c(:) ! -log(P) on main grid.
  Type (Grids_T) :: Grids_f           ! All the coordinates
  REAL(rp), INTENT(in) :: dbeta_path_c(:,:) ! derivative of beta wrt dx
!                                            on main grid.
  REAL(rp), INTENT(in) :: eta_zxp_f_c(:,:) ! representation basis function
!                                            main grid.
  REAL(rp), INTENT(in) :: sps_path_c(:,:) ! species function on  main grid.
  LOGICAL, INTENT(in) :: do_calc_f_c(:,:) ! A logical indicating where the
!                                         representation basis function is
!                                         not zero on main grid.
  REAL(rp), INTENT(in) :: dbeta_path_f(:,:) ! derivative of beta wrt dx
!                                            on gl grid.
  REAL(rp), INTENT(in) :: eta_zxp_f_f(:,:) ! representation basis function
!                                            gl grid.
  REAL(rp), INTENT(in) :: sps_path_f(:,:) ! species function on gl grid.
  LOGICAL, INTENT(in) :: do_calc_f_f(:,:) ! A logical indicating where the
!                                         representation basis function is
!                                         not zero on main grid.
  LOGICAL, INTENT(in) :: do_gl(:) ! A logical indicating where to do gl
!                                   integrations
  REAL(rp), INTENT(in) :: ref_cor(:) ! refracted to unrefracted path
!                                      length ratios.
  REAL(rp), INTENT(in) :: del_s(:) ! unrefracted path length.
  REAL(rp), INTENT(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                        gl grid.
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                       gl grid.
  REAL(rp), INTENT(in) :: t_script(:) ! differential temperatures (K).
  REAL(rp), INTENT(in) :: tau(:) ! transmission function.
  INTEGER(ip), INTENT(in) :: i_stop ! path stop index
  INTEGER(ip), INTENT(in) :: ptg_i,frq_i ! debugger statements
!
! Outputs
!
  REAL(rp), INTENT(out) :: drad_dx(:)    ! derivative of radiances wrt x
!                                          statevector element. (K)
! Internals
!
  INTEGER(ip) :: sv_i, sv_j, sps_i, n_inds, i, j, k, l, i_start, n_sps
  INTEGER(ip) :: n_path, p_i, no_to_gl, mid
  INTEGER(ip), ALLOCATABLE :: inds(:),all_inds(:),more_inds(:)
!
! use automatic allocations here
!
  REAL(rp), ALLOCATABLE :: d_delta_dx(:),singularity(:),del_zeta(:),gl_delta(:)
  LOGICAL, ALLOCATABLE :: do_calc(:)
!
! Begin code
!
  n_sps = SIZE(Grids_f%no_z)
  n_path = SIZE(tau)
  mid = n_path / 2
!
  ALLOCATE(d_delta_dx(1:n_path))
  ALLOCATE(do_calc(1:n_path))
!
  sv_i = 0
  drad_dx(:) = 0.0_rp

  DO sps_i = 1 , n_sps
!
    DO sv_j = 1 , Grids_f%no_z(sps_i) * Grids_f%no_p(sps_i)
!
      sv_i = sv_i + 1
      d_delta_dx = 0.0_rp

      do_calc = do_calc_f_c(:,sv_i)

      i = 1
      DO p_i = 1 , n_path
        IF(do_gl(p_i)) THEN
          IF(ANY(do_calc_f_f(i:i+Ng-1,sv_i))) do_calc(p_i) = .TRUE.
          i = i + Ng
        ENDIF
      ENDDO
!
! find where the non zeros are along the path
!
      n_inds = COUNT(do_calc)
      IF(n_inds > 0) THEN
!
        ALLOCATE(inds(1:n_inds))
        ALLOCATE(singularity(1:n_inds))

        i = 1
        DO p_i = 1 , n_path
          IF(do_calc(p_i)) THEN
            inds(i) = p_i
            i = i + 1
          ENDIF
        ENDDO
!
        singularity = dbeta_path_c(inds,sps_i) * eta_zxp_f_c(inds,sv_i) * &
                   &  sps_path_c(inds,sps_i)
        DO p_i = 1 , n_inds
          d_delta_dx(inds(p_i)) = singularity(p_i) * del_s(inds(p_i))
        ENDDO
!
        no_to_gl = COUNT(do_gl(inds))
        IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
          ALLOCATE(all_inds(1:Ng*no_to_gl))
          ALLOCATE(gl_delta(1:no_to_gl))
          ALLOCATE(del_zeta(1:no_to_gl))
          ALLOCATE(more_inds(1:no_to_gl))

          i = 1
          j = 1
          l = 1
          DO p_i = 1 , n_path
            IF(do_gl(p_i)) THEN
              IF(do_calc(p_i)) THEN
                more_inds(i) = p_i
                all_inds(j:j+Ng-1) = l + (/(k-1,k=1,Ng)/)
                IF(p_i > mid) THEN
                  del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
                ELSE
                  del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
                ENDIF
                i = i + 1
                j = j + Ng
              ENDIF
              l = l + Ng
            ENDIF
          END DO

          CALL path_opacity(del_zeta, &
              &  singularity(PACK((/(i,i=1,n_inds)/),do_gl(inds))),        &
              &  dbeta_path_f(all_inds,sps_i)*eta_zxp_f_f(all_inds,sv_i) * &
              &  sps_path_f(all_inds,sps_i),ds_dh_gl(all_inds),            &
              &  dh_dz_gl(all_inds),gl_delta)
          d_delta_dx(more_inds) = d_delta_dx(more_inds) + gl_delta

          DEALLOCATE(more_inds)
          DEALLOCATE(del_zeta)
          DEALLOCATE(gl_delta)
          DEALLOCATE(all_inds)

        ENDIF
!
        d_delta_dx(inds) = ref_cor(inds) * d_delta_dx(inds)
!
        i_start = MIN(inds(1),i_stop)

        DEALLOCATE(singularity)
        DEALLOCATE(inds)

        CALL get_dscrt_no_t_dn(d_delta_dx,t_script,tau,i_start,i_stop, &
                            &  drad_dx(sv_i))
!
      ENDIF
!
    ENDDO
!
  ENDDO
!
  DEALLOCATE(do_calc)
  DEALLOCATE(d_delta_dx)
!
END SUBROUTINE drad_tran_dx
!----------------------------------------------------------------------
! This is the radiative transfer derivative wrt temperature model
!
SUBROUTINE drad_tran_dt(z_path_c,h_path_c,t_path_c,dh_dt_path_c, &
                     &  alpha_path_c,alphaxn_path_c,eta_zxp_c,do_calc_t_c, &
                     &  do_calc_hyd_c,del_s,ref_cor,h_tan,dh_dt_tan, &
                     &  freq,do_gl,h_path_f,t_path_f,dh_dt_path_f, &
                     &  alpha_path_f,alphaxn_path_f,eta_zxp_f,do_calc_t_f, &
                     &  ds_dh_gl,dh_dz_gl,t_script,tau,i_stop,drad_dt,  &
                     &  ptg_i,frq_i)
!
! Inputs
!
  REAL(rp), INTENT(in) :: z_path_c(:) ! path -log(P) on main grid.
  REAL(rp), INTENT(in) :: h_path_c(:) ! path heights + req on main grid km.
  REAL(rp), INTENT(in) :: t_path_c(:) ! path temperature(K) on main grid.
  REAL(rp), INTENT(in) :: dh_dt_path_c(:,:) ! derivative of path height wrt
!                                             temperature(km/K) on main grid.
  REAL(rp), INTENT(in) :: alpha_path_c(:) ! path absorption(km^-1)
!                                           on main grid.
  REAL(rp), INTENT(in) :: alphaxn_path_c(:) ! path absorption times
!                          temperature power(km^-1) on main grid.
  REAL(rp), INTENT(in) :: eta_zxp_c(:,:) ! representation basis function
!                                            main grid.
  LOGICAL, INTENT(in) :: do_calc_t_c(:,:) ! A logical indicating where the
!                  representation basis function is not zero on main grid.
  LOGICAL, INTENT(in) :: do_calc_hyd_c(:,:) ! A logical indicating where the
!                  dh_dt function is not zero on main grid.
  REAL(rp), INTENT(in) :: del_s(:) ! unrefracted path length.
  REAL(rp), INTENT(in) :: ref_cor(:) ! refracted to unrefracted path
!                                      length ratios.
  REAL(rp), INTENT(in) :: h_tan ! tangent height + req (km).
  REAL(rp), INTENT(in) :: dh_dt_tan(:) ! derivative of path height wrt
!                                     temperature at the tangent (km/K).
  REAL(r8), INTENT(in) :: freq ! calculation frequency (MHz).
  LOGICAL, INTENT(in) :: do_gl(:) ! A logical indicating where to do gl
!                                   integrations.
  REAL(rp), INTENT(in) :: h_path_f(:) ! path heights + req on gl grid km.
  REAL(rp), INTENT(in) :: t_path_f(:) ! path temperature(K) on gl grid.
  REAL(rp), INTENT(in) :: dh_dt_path_f(:,:) ! derivative of path height wrt
!                                             temperature(km/K) on gl grid.
  REAL(rp), INTENT(in) :: alpha_path_f(:) ! path absorption(km^-1) on gl grid.
  REAL(rp), INTENT(in) :: alphaxn_path_f(:) ! path absorption times
!                          temperature power(km^-1) on gl grid.
  REAL(rp), INTENT(in) :: eta_zxp_f(:,:) ! representation basis function
!                                            gl grid.
  LOGICAL, INTENT(in) :: do_calc_t_f(:,:) ! A logical indicating where the
!                  representation basis function is not zero on gl grid.
  REAL(rp), INTENT(in) ::  ds_dh_gl(:) ! path length wrt height derivative on
!                                        gl grid.
  REAL(rp), INTENT(in) :: dh_dz_gl(:) ! path height wrt zeta derivative on
!                                       gl grid.
  REAL(rp), INTENT(in) :: t_script(:) ! differential temperatures (K).
  REAL(rp), INTENT(in) :: tau(:) ! transmission function.
  INTEGER(ip), INTENT(in) :: i_stop ! path stop index
  INTEGER(ip), INTENT(in) :: ptg_i,frq_i ! debugger statements
!
! Outputs
!
  REAL(rp), INTENT(out) :: drad_dt(:)    ! derivative of radiances wrt
!                                          mixing ratio statevector
!                                          element. (K)
! Internals
!
  INTEGER(ip) :: sv_i, sv_t, n_inds, i, j, k, l, i_start
  INTEGER(ip) :: n_path, p_i, no_to_gl, mid
  INTEGER(ip), ALLOCATABLE :: inds(:),all_inds(:),more_inds(:)
!
! use automatic allocations here
!
  REAL(rp), ALLOCATABLE :: d_delta_dt(:),singularity(:)
  REAL(rp), ALLOCATABLE :: del_zeta(:),gl_delta(:),dt_scr_dt(:,:)
  REAL(rp) :: fa,fb

  LOGICAL, ALLOCATABLE :: do_calc(:)
!
! Begin code
!
  n_path = SIZE(tau)
  sv_t = SIZE(eta_zxp_c,DIM=2)
  mid = n_path / 2
!
! allocate memory
!
  ALLOCATE(d_delta_dt(1:n_path))
  ALLOCATE(do_calc(1:n_path))
  ALLOCATE(dt_scr_dt(1:n_path,1:sv_t))
!
! compute the t_script derivative
!
  CALL dt_script_dt(t_path_c,eta_zxp_c,freq,dt_scr_dt)
!
! compute the opacity derivative singularity value
!
  drad_dt(:) = 0.0_rp

  DO sv_i = 1 , sv_t
!
    i_start = 1
    d_delta_dt = 0.0_rp
!
! do the absorption part
! combine non zeros flags for both the main and gl parts
!
    do_calc = do_calc_t_c(:,sv_i)

    i = 1
    DO p_i = 1 , n_path
      IF(do_gl(p_i)) THEN
        IF(ANY(do_calc_t_f(i:i+ng-1,sv_i))) do_calc(p_i) = .TRUE.
        i = i + Ng
      ENDIF
    ENDDO
!
! find where the non zeros are along the path
!
    n_inds = COUNT(do_calc)
    IF(n_inds > 0) THEN
      ALLOCATE(inds(1:n_inds))
      ALLOCATE(singularity(1:n_inds))
      i = 1
      DO p_i = 1 , n_path
        IF(do_calc(p_i)) THEN
          inds(i) = p_i
          i = i + 1
        ENDIF
      ENDDO
!
      singularity = alphaxn_path_c(inds)*eta_zxp_c(inds,sv_i) &
                  / t_path_c(inds)
      d_delta_dt(inds) = singularity * del_s(inds)
!
      no_to_gl = COUNT(do_gl(inds))
      IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
         ALLOCATE(all_inds(1:ng*no_to_gl))
         ALLOCATE(gl_delta(1:no_to_gl))
         ALLOCATE(del_zeta(1:no_to_gl))
         ALLOCATE(more_inds(1:no_to_gl))

         i = 1
         j = 1
         l = 1
         DO p_i = 1 , n_path
           IF(do_gl(p_i)) THEN
             IF(do_calc(p_i)) THEN
               more_inds(i) = p_i
               all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
               IF(p_i > mid) THEN
                 del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
               ELSE
                 del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
               ENDIF
               i = i + 1
               j = j + Ng
             ENDIF
             l = l + Ng
           ENDIF
         ENDDO
!
         CALL path_opacity(del_zeta, &
           &  singularity(PACK((/(i,i=1,n_inds)/),do_gl(inds))), &
           &  alphaxn_path_f(all_inds)*eta_zxp_f(all_inds,sv_i) &
           &  / t_path_f(all_inds),ds_dh_gl(all_inds),dh_dz_gl(all_inds), &
           &  gl_delta)
         d_delta_dt(more_inds) = d_delta_dt(more_inds) + gl_delta
!
         DEALLOCATE(more_inds)
         DEALLOCATE(del_zeta)
         DEALLOCATE(gl_delta)
         DEALLOCATE(all_inds)
!
      ENDIF
!
      i_start = MAX(inds(1)-1,1)
!
      DEALLOCATE(singularity)
      DEALLOCATE(inds)
!
    ENDIF
!
! now do the hydrostatic part
! combine boundaries flags
!
    do_calc = do_calc_hyd_c(:,sv_i)
    IF(do_calc_hyd_c(mid,sv_i)) THEN
      do_calc(2:mid) = .TRUE.
    ELSE
      DO p_i = 2 , mid
        IF(do_calc_hyd_c(p_i-1,sv_i)) do_calc(p_i) = .TRUE.
      ENDDO
    ENDIF
    IF(do_calc_hyd_c(mid+1,sv_i)) THEN
      do_calc(mid+1:n_path-1) = .TRUE.
    ELSE
      DO p_i = mid + 1 , n_path - 1
        IF(do_calc_hyd_c(p_i+1,sv_i)) do_calc(p_i) = .TRUE.
      ENDDO
    ENDIF
!
! since this is a layer boundary calculation we must require
!
    do_calc((/1,n_path/)) = .FALSE.
!
! find where the non zeros are along the path
!
    n_inds = COUNT(do_calc)
    IF(n_inds > 0) THEN
      ALLOCATE(inds(1:n_inds))
      i = 1
      inds = 0
      DO p_i = 2 , mid - 1
        IF(do_calc(p_i)) THEN
          IF(i == 1) fa = (h_path_c(p_i-1)*dh_dt_path_c(p_i-1,sv_i) &
                        - h_tan * dh_dt_tan(sv_i)) / SUM(del_s(p_i:mid))
          fb = (h_path_c(p_i)*dh_dt_path_c(p_i,sv_i) &
             - h_tan * dh_dt_tan(sv_i)) /SUM(del_s(p_i+1:mid))
          inds(i) = p_i
          d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i)*(fa - fb)
          fa = fb
          i = i + 1
        ENDIF
      ENDDO
!
! special processing at tangent
!
      IF(do_calc(mid)) THEN
        d_delta_dt(mid) = d_delta_dt(mid) + alpha_path_c(mid) * fa
!
! fb is zero
!
        inds(i) = mid
        i = i + 1
      ENDIF
!
      IF(do_calc(mid+1)) THEN
        fa = (h_path_c(mid+2)*dh_dt_path_c(mid+2,sv_i) &
           - h_tan * dh_dt_tan(sv_i))/ del_s(mid+1)
        d_delta_dt(mid+1) = d_delta_dt(mid+1) + alpha_path_c(mid+1) * fa
!
! fb is 0.0
!
        inds(i) = mid+1
        i = i + 1
      ENDIF
!
      DO p_i = mid + 2 , n_path - 1
        IF(do_calc(p_i)) THEN
          IF(inds(MAX(i-1,1)) < mid+1)  &
                     & fa = (h_path_c(p_i)*dh_dt_path_c(p_i,sv_i) &
                     &    - h_tan * dh_dt_tan(sv_i)) / SUM(del_s(mid+1:p_i-1))
          fb = (h_path_c(p_i+1)*dh_dt_path_c(p_i+1,sv_i) &
             - h_tan * dh_dt_tan(sv_i)) /SUM(del_s(mid+1:p_i))
          inds(i) = p_i
          d_delta_dt(p_i) = d_delta_dt(p_i) + alpha_path_c(p_i)*(fb - fa)
          fa = fb
          i = i + 1
        ENDIF
      ENDDO
!
      no_to_gl = COUNT(do_gl(inds))
      IF(no_to_gl > 0) THEN
!
! see if anything needs to be gl-d
!
         ALLOCATE(all_inds(1:ng*no_to_gl))
         ALLOCATE(gl_delta(1:no_to_gl))
         ALLOCATE(del_zeta(1:no_to_gl))
         ALLOCATE(more_inds(1:no_to_gl))

         i = 1
         j = 1
         l = 1
         DO p_i = 1 , n_path
           IF(do_gl(p_i)) THEN
             IF(do_calc(p_i)) THEN
               more_inds(i) = p_i
               all_inds(j:j+ng-1) = l + (/(k-1,k=1,ng)/)
               IF(p_i > mid) THEN
                 del_zeta(i) = z_path_c(p_i+1) - z_path_c(p_i)
               ELSE
                 del_zeta(i) = z_path_c(p_i-1) - z_path_c(p_i)
               ENDIF
               i = i + 1
               j = j + Ng
             ENDIF
             l = l + Ng
           ENDIF
         ENDDO
!
! add special hydrostatic gl routine here
! the singularity point is alpha_path_c(more_inds)
!
         CALL hyd_opacity(del_zeta,alpha_path_c(more_inds), &
              alpha_path_f(all_inds),h_path_f(all_inds), &
              dh_dt_path_f(all_inds,sv_i),t_path_f(all_inds),h_tan, &
              dh_dt_tan(sv_i),eta_zxp_f(all_inds,sv_i),ds_dh_gl(all_inds), &
              dh_dz_gl(all_inds),gl_delta)
         d_delta_dt(more_inds) = d_delta_dt(more_inds) + gl_delta
!
         DEALLOCATE(more_inds)
         DEALLOCATE(del_zeta)
         DEALLOCATE(gl_delta)
         DEALLOCATE(all_inds)
!
      ENDIF
!
      i_start = MIN(i_start,inds(1))
!
      DEALLOCATE(inds)
!
    ENDIF
!
! correct the whole thing for path length refraction
!
    d_delta_dt = ref_cor * d_delta_dt
    CALL get_dscrt_dn(d_delta_dt,t_script,tau,dt_scr_dt(:,sv_i),i_start, &
                   &  i_stop,drad_dt(sv_i))
!
  ENDDO
!
  DEALLOCATE(d_delta_dt)
  DEALLOCATE(do_calc)
  DEALLOCATE(dt_scr_dt)
!
END SUBROUTINE drad_tran_dt
!
!----------------------------------------------------------------------
End module RAD_TRAN_M
! $Log$
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
