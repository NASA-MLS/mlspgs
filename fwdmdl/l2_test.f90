Program l2_Test
  use GL6P, only: NG
  use MLSCommon, only: I4, R4, R8
  use L2_TEST_STRUCTURES_M
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
  use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi, &
                     MNM => max_no_mmaf
  use ELLIPSE, only: PHI_TAN, ROC
  use L2_LOAD_M, only: L2_LOAD
  use PTG_FRQ_LOAD_M, only: PTG_FRQ_LOAD
  use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use RAD_TRAN_M, only: RAD_TRAN
  use RAD_TRAN_WD_M, only: RAD_TRAN_WD
  use FREQ_AVG_M, only: FREQ_AVG
  use CONVOLVE_ALL_M, only: CONVOLVE_ALL
  use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
  use D_HUNT_M, only: HUNT          ! ** DEBUG
  implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl

Type(fwd_mdl_info) :: FMI
Type(fwd_mdl_config) :: FMC
Type(temporary_fwd_mdl_info) :: T_FMI

Integer(i4) :: i, j, k, kk, kz, ht_i, mnz, no_tan_hts, ch, Spectag, &
               m, prev_npf, ier, mmaf, si, ptg_i, &
               frq_i, io, klo, jj, l, n, brkpt, no_ele, mid, ilo, ihi, &
               k_info_count, gl_count, ld

Integer(i4) :: ch1, ch2, no_pfa_ch, pfa_ch(2)

Type(path_index)  :: ndx_path(Nptg,mnm)
Type(path_vector) :: z_path(Nptg,mnm),t_path(Nptg,mnm),h_path(Nptg,mnm),  &
                     dhdz_path(Nptg,mnm), spsfunc_path(Nsps,Nptg,mnm),    &
                     n_path(Nptg,mnm),phi_path(Nptg,mnm)

Type(path_derivative) :: dh_dt_path(Nptg,mnm)

Real(r8) :: thbs(10),elev_offset
Real(r8) :: t_script(N2lvl),ref_corr(N2lvl,Nptg),tau(N2lvl), &
            tan_dh_dt(Nlvl,mnm,mxco)

Real(r8) :: dx_dt(Nptg,mxco), d2x_dxdt(Nptg,mxco)

Real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
Real(r8) :: dh_dt_glgrid(ngt,mnm,mxco), dhdz_glgrid(ngt,mnp)

Real(r8) :: ptg_angles(Nptg,mnm), center_angle
Real(r8) :: tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)


! Real(r4) :: K_TEMP(Nch,Nptg,mxco,mnp)
! Real(r4) :: K_ATMOS(Nch,Nptg,mxco,mnp,Nsps)
! Real(r4) :: K_SPECT_DW(Nch,Nptg,mxco,mnp,Nsps),  &
!             K_SPECT_DN(Nch,Nptg,mxco,mnp,Nsps),  &
!             K_SPECT_DNU(Nch,Nptg,mxco,mnp,Nsps)

! ** DEBUG, memory limitations force us to have up to 2 channels
!           only (Replacing Nch by: 2)

Real(r4) :: K_TEMP(02,Nptg,mxco,mnp)
Real(r4) :: K_ATMOS(02,Nptg,mxco,mnp,Nsps)
Real(r4) :: K_SPECT_DW(02,Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DN(02,Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DNU(02,Nptg,mxco,mnp,Nsps)

real(r8) :: I_STAR_ALL(Nch,Nptg)

real(r4) :: K_STAR_ALL(02,20,mxco,mnp,Nptg)      ! 02 should be: Nch
Type(k_matrix_info) :: k_star_info(20)

Type(path_derivative) :: k_temp_frq, k_atmos_frq(Nsps), &
                         k_spect_dw_frq(Nsps), k_spect_dn_frq(Nsps), &
                         k_spect_dnu_frq(Nsps)
!
Type(path_beta), DIMENSION(:,:), POINTER :: beta_path

Real(r8) :: Radiances(Nptg,Nch)
Real(r8) :: e_rad, Zeta, Frq, h_tan, Rad, geoc_lat, r
!
Character (LEN=01) :: CA
Character (LEN=08) :: Name
Character (LEN=16) :: Vname
Character (LEN=80) :: Line
Character (LEN=40) :: Ax, Dtm1, Dtm2

Real(r8), DIMENSION(:), ALLOCATABLE :: RadV, F_grid

!  ----------------------

  Line(1:) = ' '
! Line = '/home/zvi/seez'      ! HOME PC, JPL PC
  Line = '/user5/zvi/seez'     ! MLSGATE, SUN
  FMC%Z = Line
!
! Load all the needed L2 data:
!
  Dtm1(1:)=' '
  Dtm2(1:)=' '
  Call Z_DATETIME(Dtm1)

  Call L2_LOAD(FMC, FMI, T_FMI, Ier)
  if(ier /= 0) goto 99

  ch1 = FMC%Channels_range(1)
  ch2 = FMC%Channels_range(2)
  no_pfa_ch = min(2,ch2-ch1+1)
  do i = 1, no_pfa_ch
    pfa_ch(i) = ch1 + i - 1
  end do

  elev_offset = 0.0                         ! Zero elev_offset in any case

  Call Z_DATETIME(Dtm1)

  Call ptg_frq_load(FMC, FMI, Ier)
  if(ier /= 0) goto 99
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  mmaf = 3                     ! Do only this mmaf (middle phi)
  phi_tan = FMC%phi_tan_mmaf(mmaf)
  Call geoc_geod_conv(T_FMI%beta_inc,phi_tan,geoc_lat,E_rad)
!
! Compute the hydrostatic_model on the GL-Grid for all mmaf(s):
!
  thbs(1:) = 0.0
  si = FMI%Surface_index
  thbs(1:si-1) = FMI%Tan_hts_below_surface(1:si-1)
  Call hydrostatic_model(si,FMC%N_lvls,T_FMI%no_t,FMC%no_mmaf,FMC%t_indx, &
       FMC%no_tan_hts,geoc_lat,T_FMI%Href,T_FMI%Zref,FMI%z_grid,thbs, &
       T_FMI%t_zeta_basis, T_FMI%t_coeff, z_glgrid, h_glgrid, t_glgrid, &
       dhdz_glgrid,dh_dt_glgrid,FMI%tan_press,tan_hts,tan_temp,tan_dh_dt, &
       gl_count, Ier)
  IF(ier /= 0) goto 99
!
  Zeta = -1.666667
  no_tan_hts = FMC%no_tan_hts
  Call Hunt(Zeta,FMI%tan_press,no_tan_hts,jj,i)
  IF(ABS(Zeta-FMI%tan_press(i)) < ABS(Zeta-FMI%tan_press(jj))) jj = i
!
! Compute all path entities for all mmafs and tanget pointings
!
  Call comp_path_entities(FMC%n_lvls,T_FMI%no_t,gl_count,ndx_path, &
       z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,        &
       T_FMI%atmospheric,T_FMI%f_zeta_basis,T_FMI%mr_f,            &
       T_FMI%no_coeffs_f,tan_hts,no_tan_hts,FMI%n_sps,             &
       T_FMI%no_phi_f,T_FMI%f_phi_basis,z_path,h_path,t_path,      &
       phi_path,n_path,dhdz_path,dh_dt_path,T_FMI%no_phi_t,        &
       T_FMI%t_phi_basis,spsfunc_path,T_FMI%is_f_log,FMC%no_mmaf,  &
       FMC%phi_tan_mmaf,Ier)
  IF(ier /= 0) goto 99
!
! **********************  MAIN Mmaf Loop *******************

! DO l = 1, FMC%no_mmaf
  DO l = mmaf, mmaf                 ! ** DEBUG, only one mmaf
!
    phi_tan = FMC%phi_tan_mmaf(l)
!
    T_FMI%t_phi_basis(1:T_FMI%no_phi_t) = &
                    T_FMI%T_PHI_BASIS_COPY(1:T_FMI%no_phi_t) + phi_tan

    DO j = 1, FMI%n_sps
      k = T_FMI%no_phi_f(j)
      T_FMI%f_phi_basis(1:k,j) = T_FMI%F_PHI_BASIS_COPY(1:k,j) + phi_tan
    end do

    k = FMI%mfi + 2
    do j = 1, FMI%no_spectro
      FMI%spectroscopic(j)%PHI_BASIS(1:k) = &
     &                T_FMI%S_PHI_BASIS_COPY(1:k,j) + phi_tan
    end do
!
! Compute the ptg_angles (chi) for Antenna convolution, also the derivatives
! of chi w.r.t to T and other parameters
!
    Call get_chi_angles(ndx_path(1:,l),n_path(1:,l),FMI%tan_press,         &
   &     tan_hts(1:,l),tan_temp(1:,l),phi_tan,RoC,T_FMI%h_obs,elev_offset, &
   &     tan_dh_dt(1:,l,1:),no_tan_hts,T_FMI%no_t,T_FMI%t_zeta_basis,si,   &
   &     center_angle,ptg_angles(1:,l),dx_dt,d2x_dxdt,ier)
    IF(ier /= 0) goto 99
!
! Compute the refraction correction scaling matrix for this mmaf:
!
    Call refraction_correction(no_tan_hts, tan_hts(1:,l), h_path(1:,l), &
   &                n_path(1:,l), ndx_path(1:,l), E_rad, ref_corr)
!
    prev_npf = -1
    Radiances(1:Nptg,1:Nch) = 0.0
!
! **********************  MAIN Pointing Loop *******************
!
    DO ptg_i = 1, no_tan_hts-1
!
      k = ptg_i
      h_tan = tan_hts(k,l)
      kk = FMI%no_ptg_frq(k)
!
      if(kk /= prev_npf) then
!
        prev_npf = kk
        DEALLOCATE(k_temp_frq%values,STAT=i)
!
        DEALLOCATE(RadV,F_grid,STAT=i)
        ALLOCATE(RadV(kk),F_grid(kk),STAT=ier)
        IF(ier /= 0) then
          Print *,'** ALLOCATE Error: RadV or F_grid arrays, STAT =',ier
          goto 99
        endif
!
        do j = 1, FMI%n_sps
          DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
        end do
!
        ALLOCATE(k_temp_frq%values(kk,T_FMI%no_t,T_FMI%no_phi_t),STAT=ier)
        IF(ier /= 0) then
          Print *,'** ALLOCATE Error: k_temp_frq array, STAT =',ier
          goto 99
        endif
!
        do j = 1, FMI%n_sps
          m = max(1,T_FMI%no_phi_f(j))
          i = max(1,T_FMI%no_coeffs_f(j))
          ALLOCATE(k_atmos_frq(j)%values(kk,i,m),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: k_atmos_frq, STAT =',ier
            goto 99
          endif
        end do
!
        do m = 1, FMI%n_sps
          j = FMI%spect_atmos(m)
          if(j < 1) CYCLE
          if(.not. FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
          Vname(1:) = ' '
          Spectag = FMI%spectroscopic(j)%Spectag
          DO
            if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
            n = FMI%spectroscopic(j)%no_phi_values
            i = FMI%spectroscopic(j)%no_zeta_values
            CA = FMI%spectroscopic(j)%type
            select case ( CA )
              case ( 'W' )
                Vname = 'k_spect_dw_frq'
                ALLOCATE(k_spect_dw_frq(m)%values(kk,i,n),STAT=ier)
              case ( 'N' )
                Vname = 'k_spect_dn_frq'
                ALLOCATE(k_spect_dn_frq(m)%values(kk,i,n),STAT=ier)
              case ( 'V' )
                Vname = 'k_spect_dnu_frq'
                ALLOCATE(k_spect_dnu_frq(m)%values(kk,i,n),STAT=ier)
              case default
                Ier = -99
                Print *,'** Unknown Spectroscopic element !'
            end select
            IF(ier /= 0) then
              Print *,'** ALLOCATE Error: ',Vname,', STAT =',ier
              goto 99
            ENDIF
            j = j + 1
            if(j > 3 * FMI%n_sps) EXIT
          END DO
        end do

      endif            ! On DEALLOCATE/ALLOCATE cycle
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
      no_ele = ndx_path(ptg_i,l)%total_number_of_elements
      Call get_beta_path(ptg_i,FMI%pfa_spectrum,no_ele,FMI%no_ptg_frq, &
     &     FMI%ptg_frq_grid,z_path(ptg_i,l),t_path(ptg_i,l),beta_path,ier)
      IF(ier /= 0) goto 99
!
      k_temp_frq%values = 0.0
      do j = 1, FMI%n_sps
        k_atmos_frq(j)%values = 0.0
        k_spect_dw_frq(j)%values = 0.0
        k_spect_dn_frq(j)%values = 0.0
        k_spect_dnu_frq(j)%values = 0.0
      end do
!
      RadV(1:kk) = 0.0
      F_grid(1:kk) = FMI%ptg_frq_grid(k)%values(1:kk)
!
      do frq_i = 1, kk
!
        Frq = F_grid(frq_i)
!
        Call Rad_Tran(Frq, FMC%N_lvls, h_tan, FMI%n_sps, ndx_path(k,l),  &
       &    z_path(k,l), h_path(k,l), t_path(k,l), phi_path(k,l),&
       &    dHdz_path(k,l), T_FMI%earth_ref, beta_path(1:,frq_i),      &
       &    spsfunc_path(1:,k,l), ref_corr(1:,k), T_FMI%s_temp, brkpt, &
       &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
        IF(ier /= 0) goto 99
!
        RadV(frq_i) = Rad
!
! Now, Compute the radiances derivatives:
!
        CALL Rad_Tran_WD(frq_i,FMI%band,Frq,FMC%N_lvls,FMI%n_sps, &
       &     z_path(k,l),h_path(k,l),t_path(k,l),phi_path(k,l),   &
       &     dHdz_path(k,l),T_FMI%atmospheric,beta_path(1:,frq_i),&
       &     spsfunc_path(1:,k,l),T_FMI%t_zeta_basis,  &
       &     T_FMI%f_zeta_basis,T_FMI%no_coeffs_f,   &
       &     T_FMI%mr_f,T_FMI%no_t,ref_corr(1:,k),T_FMI%no_phi_f,       &
       &     T_FMI%f_phi_basis,FMC%temp_der,T_FMI%no_phi_t,             &
       &     T_FMI%t_phi_basis,dh_dt_path(k,l),FMI%spect_atmos,         &
       &     FMI%spectroscopic,k_temp_frq,k_atmos_frq,k_spect_dw_frq,   &
       &     k_spect_dn_frq,k_spect_dnu_frq,T_FMI%is_f_log,brkpt,       &
       &     no_ele,mid,ilo,ihi,t_script,tau,ier)
        IF(ier /= 0) goto 99
!
      end do

! Frequency Average the radiances with the appropriate filter shapes
!
      do i = 1, no_pfa_ch
        ch = pfa_ch(i)
        if(FMC%do_frqavg) then
          Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i),  &
         &     FMI%Filter_func(1:,i),RadV,kk,FMI%no_filt_pts, &
         &     Radiances(ptg_i,ch))
        else
          Radiances(ptg_i,ch) = RadV(1)
        endif
      end do

!     if(i > -3) goto 77           ! ** DEBUG, Bypass derivatives avg.
!
! Frequency Average the temperature derivatives with the appropriate
! filter shapes
!
      if(FMC%temp_der) then
        RadV(1:kk) = 0.0
        do i = 1, no_pfa_ch
!         ch = pfa_ch(i)
          ch = i               ! ** DEBUG, memory limitations on MLSGATE
          do j = 1, T_FMI%no_phi_t
            do k = 1, T_FMI%no_t
              if(FMC%do_frqavg) then
                RadV(1:kk) = k_temp_frq%values(1:kk,k,j)
                Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
               &              FMI%Filter_func(1:,i),          &
               &              RadV,kk,FMI%no_filt_pts,r)
              else
                r = k_temp_frq%values(1,k,j)
              endif
              k_temp(ch,ptg_i,k,j) = r
            end do
          end do
        end do
      endif

!     if(i > -3) goto 77           ! ** DEBUG, Bypass derivatives avg.
!
! Frequency Average the atmospheric derivatives with the appropriate
! filter shapes
!
      do i = 1, no_pfa_ch
!       ch = pfa_ch(i)
        ch = i                 ! ** DEBUG, memory limitations on MLSGATE
        do j = 1, FMI%n_sps
          if(T_FMI%atmospheric(j)%der_calc(FMI%band)) THEN
            RadV(1:kk) = 0.0
            do k = 1, T_FMI%no_phi_f(j)
              do n = 1, T_FMI%no_coeffs_f(j)
                if(FMC%do_frqavg) then
                  RadV(1:kk) = k_atmos_frq(j)%values(1:kk,n,k)
                  Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
                 &              FMI%Filter_func(1:,i),          &
                 &              RadV,kk,FMI%no_filt_pts,r)
                else
                  r = k_atmos_frq(j)%values(1,n,k)
                endif
                k_atmos(ch,ptg_i,n,k,j) = r
              end do
            end do
          endif
        end do
      end do

!     if(i > -3) goto 77           ! ** DEBUG, Bypass derivatives avg.
!
! Frequency Average the spectroscopic derivatives with the appropriate
! filter shapes
!
      do i = 1, no_pfa_ch
!       ch = pfa_ch(i)
        ch = i                 ! ** DEBUG, memory limitations on MLSGATE
        do m = 1, FMI%n_sps
          j = FMI%spect_atmos(m)
          if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
          Spectag = FMI%spectroscopic(j)%Spectag
          DO
            if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
            RadV(1:kk) = 0.0
            CA = FMI%spectroscopic(j)%type
            do k = 1, FMI%spectroscopic(j)%no_phi_values
              do n = 1, FMI%spectroscopic(j)%no_zeta_values
                select case ( CA )
                  case ( 'W' )
                    RadV(1:kk) = k_spect_dw_frq(m)%values(1:kk,n,k)
                  case ( 'N' )
                    RadV(1:kk) = k_spect_dn_frq(m)%values(1:kk,n,k)
                  case ( 'V' )
                    RadV(1:kk) = k_spect_dnu_frq(m)%values(1:kk,n,k)
                end select
                if(FMC%do_frqavg) then
                  Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i), &
                 &              FMI%Filter_func(1:,i),&
                 &              RadV,kk,FMI%no_filt_pts,r)
                else
                  r = RadV(1)
                endif
                select case ( CA )
                  case ( 'W' )
                    k_spect_dw(ch,ptg_i,n,k,j) = r
                  case ( 'N' )
                    k_spect_dn(ch,ptg_i,n,k,j) = r
                  case ( 'V' )
                    k_spect_dnu(ch,ptg_i,n,k,j) = r
                end select
              end do
            end do
            j = j + 1
            if(j > 3 * FMI%n_sps) EXIT
          END DO
        end do
      end do
!
 77   j = 0             ! ** DEBUG
!
    END DO              ! Pointing Loop
!
! Complete the radiances's last location, also  complete k_temp last
! location as well as k_atmos last location and k_spect_d? last location:
!
    kk = no_tan_hts
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      Radiances(kk,ch) = Radiances(kk-1,ch)
      k_temp(i,kk,1:T_FMI%no_t,1:T_FMI%no_phi_t) = &
     &            k_temp(i,kk-1,1:T_FMI%no_t,1:T_FMI%no_phi_t)
      do m = 1, FMI%n_sps
        if(T_FMI%atmospheric(m)%der_calc(FMI%band)) then
          k = T_FMI%no_phi_f(m)
          n = T_FMI%no_coeffs_f(m)
          k_atmos(i,kk,1:n,1:k,m)=k_atmos(i,kk-1,1:n,1:k,m)
        endif
      end do
      do m = 1, FMI%n_sps
        j = FMI%spect_atmos(m)
        if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
        Spectag =  FMI%spectroscopic(j)%Spectag
        DO
          if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
          k = FMI%spectroscopic(j)%no_phi_values
          n = FMI%spectroscopic(j)%no_zeta_values
          k_spect_dw(i,kk,1:n,1:k,j)=k_spect_dw(i,kk-1,1:n,1:k,j)
          k_spect_dn(i,kk,1:n,1:k,j)=k_spect_dn(i,kk-1,1:n,1:k,j)
          k_spect_dnu(i,kk,1:n,1:k,j)=k_spect_dnu(i,kk-1,1:n,1:k,j)
          j = j + 1
          if(j > 3 * FMI%n_sps) EXIT
        END DO
      end do
    end do
!
!  Here comes the Convolution code
!
    DO i = 1, no_pfa_ch
!
      ch = pfa_ch(i)
!
      if(FMC%do_conv) then
!
        Call convolve_all(T_FMI%ptg_press,T_FMI%atmospheric,FMI%n_sps,   &
       &     FMC%temp_der,FMI%tan_press,ptg_angles(1:,l),tan_temp(1:,l), &
       &     dx_dt, d2x_dxdt,FMI%band,center_angle,FMI%fft_pts,          &
       &     Radiances(1:,ch),k_temp(i,1:,1:,1:),k_atmos(i,1:,1:,1:,1:), &
       &     k_spect_dw(i,1:,1:,1:,1:),k_spect_dn(i,1:,1:,1:,1:),    &
       &     k_spect_dnu(i,1:,1:,1:,1:),FMI%spect_atmos,no_tan_hts,  &
       &     k_info_count,i_star_all(i,1:),k_star_all(i,1:,1:,1:,1:), &
       &     k_star_info,T_FMI%no_t,T_FMI%no_phi_t,T_FMI%no_phi_f,   &
       &     FMI%spectroscopic,T_FMI%t_zeta_basis,FMI%Xlamda,FMI%Aaap,&
       &     FMI%D1Aaap,FMI%D2Aaap,FMI%Ias,ier)
        IF(ier /= 0) goto 99
!
      else
!
        Call no_conv_at_all(T_FMI%ptg_press,FMI%n_sps,FMI%tan_press, &
       &     FMI%band,Radiances(1:,ch),k_temp(i,1:,1:,1:),           &
       &     k_atmos(i,1:,1:,1:,1:),k_spect_dw(i,1:,1:,1:,1:),       &
       &     k_spect_dn(i,1:,1:,1:,1:),k_spect_dnu(i,1:,1:,1:,1:),   &
       &     FMI%spect_atmos, no_tan_hts,k_info_count,               &
       &     i_star_all(i,1:), k_star_all(i,1:,1:,1:,1:),            &
       &     k_star_info,FMC%temp_der,T_FMI%no_t,T_FMI%no_phi_t,     &
       &     T_FMI%no_phi_f,T_FMI%t_zeta_basis,T_FMI%atmospheric,    &
       &     FMI%spectroscopic)
!
      endif
!
    END DO

  END DO                ! Mmaf Loop
!
  DEALLOCATE(k_temp_frq%values,STAT=i)
  do j = 1, FMI%n_sps
    DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
  end do
!
! *** DEBUG Print
!
    if(FMC%do_conv) then
      Print *,'Convolution: ON'
    else
      Print *,'Convolution: OFF'
    endif
!
    if(FMC%Zfrq > 0.0) then
      Frq = FMC%Zfrq
      write(*,901) Frq
901   format(' Frequency Averaging: OFF',/,  &
        &    ' (All computations done at Frq =',f12.4,')')
    else
      Print *,'Frequency Averaging: ON'
    endif
    Print *
!
    tau(1:Nptg) = 0.0
    kk = T_FMI%ptg_press%no_lin_values
    tau(1:kk) = dble(T_FMI%ptg_press%lin_val(1:kk))
    Call Hunt(Zeta,tau,kk,klo,j)
    IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(klo))) klo=j
!
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      write(*,903) ch,char(92),kk
      write(*,905) (i_star_all(i,k),k=1,kk)
    end do
903 format('ch',i2.2,'_avg_conv_pfa_rad',a1,i2.2)
905 format(4(2x,1pg15.8))
!
    ch = 1
    tau(1:) = 0.0
    do i = 1, k_info_count
      Print *
      Name(1:) = ' '
      Name = k_star_info(i)%name
      if(Name == 'PTAN') CYCLE
      kz = k_star_info(i)%first_dim_index
      mnz = k_star_info(i)%no_zeta_basis
      ht_i = k_star_info(i)%no_phi_basis
      l = LEN_TRIM(Name)
      if(Name(l-1:l) == '_W' .or.  &
     &   Name(l-1:l) == '_N' .or.  &
     &   Name(l-1:l) == '_V' ) then
        Print *,Name
        r = SUM(k_star_all(ch,kz,1:mnz,1:ht_i,klo))
        Print *,'  Sum over all zeta & phi coeff:',sngl(r)
      else
        if(Name == 'TEMP') then
          write(6,913) 'dI_dT',char(92),ht_i
        else
          write(6,913) 'dI_d'//Name(1:l),char(92),ht_i
        endif
        tau(1:) = 0.0
        tau(1:mnz) = k_star_info(i)%zeta_basis(1:mnz)
        Call Hunt(Zeta,tau,mnz,m,j)
        IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(m))) m=j
        Print *,(k_star_all(ch,kz,m,kk,klo),kk=1,ht_i)
      endif
    end do
913 format(a,a1,i2.2)
!
 99  if(io /= 0) then
       Call ErrMsg(Line,io)
       Stop
     endif

     Call Z_DATETIME(Dtm2)
     Print *
     Print *,'** Actual computations started at: ',Dtm1
     Print *,'**            Program finished at: ',Dtm2

     Ax(1:) = ' '
     Read(Dtm1(14:21),'(i2,1x,i2,1x,i2)') i,j,k       ! hh:mm:ss
     Read(Dtm2(14:21),'(i2,1x,i2,1x,i2)') kk,ch,ld    ! hh:mm:ss
     r = Real(ld-k)+60.0*Real(ch-j)+3600.0*Real(kk-i)
     if(r < 0.0) r = r + 86400.00
     i = Int(r/3600.0)
     r = r - 3600.0*Real(i)
     j = Int(r/60.0)
     r = r - 60.0*Real(j)
     k = Int(r)
     Write(Ax(14:21),'(i2.2,'':'',i2.2,'':'',i2.2)') i,j,k
     Print *,'**         Wall-Clock elapse time: ',Ax

     Stop

!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------

SUBROUTINE Z_DATETIME(dtm)

! Output: Date & Time array in the format:

!    dd-mm-ccyy  hh:mm:ss


CHARACTER (LEN=*), INTENT(OUT) :: dtm

CHARACTER (LEN=10) :: tm
CHARACTER (LEN=8) :: dt

INTEGER :: io,iy,jm,id,ih,im,is

CHARACTER (LEN=3) :: months(12) = (/             &
         &  'JAN','FEB','MAR','APR','MAY','JUN', &
         &  'JUL','AUG','SEP','OCT','NOV','DEC'/)

! Begin code:

  dtm(1:) = ' '
  CALL date_and_time(dt,tm)

!  Dt = ccyymmdd

  READ(dt,'(i4,i2,i2)',IOSTAT=io) iy,jm,id
  IF(io /= 0) RETURN

!  Tm = hhmmss.sss

  READ(tm,'(i2,i2,i2)',IOSTAT=io) ih,im,is

  dtm(1:3) = months(jm)
  dtm(4:4) = '/'
  WRITE(dtm(5:6),'(i2.2)') id
  dtm(7:7) = '/'
  WRITE(dtm(8:11),'(i4.4)') iy

  WRITE(dtm(14:21),'(i2.2,'':'',i2.2,'':'',i2.2)') ih,im,is

  RETURN
END SUBROUTINE Z_DATETIME

end Program l2_Test
