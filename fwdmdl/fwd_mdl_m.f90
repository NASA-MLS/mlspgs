module FWD_MDL_M
  use COMPLETE_I_K_M, only: COMPLETE_I_K
  use D_T_SCRIPT_DTNP_M, only: D_T_SCRIPT_DTNP
  use DO_T_SCRIPT_M, only: DO_T_SCRIPT
  use DUMP_K_RECORDS_M, only: DUMP_K_RECORDS
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE, only: EARTHX, HT, HT2, NPHI_TAN, NPHI_S, PHI_S, PHI_TAN, &
                     ROC, RR
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES
  use GET_DELTA_M, only: GET_DELTA
  use GET_ETA_M, only: GET_ETA
  use GL6P, only: NG
  use I_K_STAR_NO_CONV_M, only: I_K_STAR_NO_CONV
  use I_AND_K_STAR_M, only: I_AND_K_STAR
  use L2PC_FILE_PARAMETERS, only: MAX_NO_POINTINGS,                     &
                                  MSVD => max_no_sv_derivatives,        &
                                  MXCO => max_no_elmnts_per_sv_component
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO,     &
                                  L2PC_HEADER_TRI
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS,  &
                                 MAXGEOPHYS, PFA_SLAB, SPECTRO_PARAM
  use L2PCDim, only: MAXFFT, MNP => max_no_phi, NLVL, NSPS, N2LVL
  use MDBETA, only: MNF => max_no_freq, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use PFA_FWD_MDL_M, only: PFA_FWD_MDL
  use S_LINTRP_M, only: LINTRP
  use SCRT_DN_M, only: SCRT_DN
  use TEMPERATURE_DERIV_M, only: TEMPERATURE_DERIV
  use WRITE_K_RECORDS_M, only: WRITE_K_RECORDS
  use WRITE_X_I_RECORDS_M, only: WRITE_X_I_RECORDS
  use ZATMOS_DERIV_M, only: ZATMOS_DERIV
  use ZOPACITY_M, only: ZOPACITY
  implicit NONE
  private
  public :: FWD_MDL

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains
!----------------------------------------------------------------------
! This is the forward model with derivatives
!
  Subroutine FWD_MDL ( header1, header2, header3, N_lvls, fft_pts,         &
 &           no_conv_hts, n_tan, band, ndx_sps, no_filt_pts, no_geom,      &
 &           no_geophys, no_int_frqs, no_pfa_ch, no_sps_tbl, sps_tbl,      &
 &           pfa_ch, nvr, z_path, h_path, t_path, phi_path, dHdz_path,     &
 &           mdb_pres, mdb_temp, mdb_freq, z_grid, t_grid, h_grid,         &
 &           dh_dt_grid, a_grid, n_grid, v_grid, t_indx, conv_hts,         &
 &           conv_press, conv_temp, spsfunc, cen_angle, earth_ref,         &
 &           t_tan,  Acc, Threshold, atmospheric, cs, d2x_dxdt, dx_dt,     &
 &           f_grid_filter, f_grid, filter_func, freq, no_freqs_f,         &
 &           g_basis, geometric, geophysic, f_basis, no_coeffs_f,          &
 &           mr_f,  no_coeffs_g, mr_g, pfa_spectrum, ptg_angle,            &
 &           ptg_hts,  ptg_press, ref_corr, azim_183, azim_205,            &
 &           azim_ref, c_pitch,  c_roll, c_yaw, elev_183, elev_205,        &
 &           geocsrad, s_pitch, s_roll,  sp_tmp, s_yaw, si, Aaap,          &
 &           InDir, keys, runf, time_i, l2pc_lu,  recn, mxbin_recn,        &
 &           lrun, rec_nos, no_atmos, jkey, n_obs,  atmos_index,           &
 &           geom_index, geophys_index, do_conv, ch1, ch2,                 &
 &           spect_atmos, no_phi_g, no_phi_f, phi_basis_f, no_phi_vec,     &
 &           Npath, path_brkpt, no_phi_t, t_phi_basis, dh_dt_path,         &
             mdb_hdr, mdb_rec, spectroscopic, Ier )
!
    type(l2pc_header_one), intent(in) :: header1
    type(l2pc_header_two), intent(in) :: header2
    type(l2pc_header_tri), intent(in) :: header3
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: FFT_PTS
    integer(i4), intent(in) :: NO_CONV_HTS
    integer(i4), intent(in) :: N_TAN(*)
    integer(i4), intent(in) :: BAND
    integer(i4), intent(in) :: NDX_SPS(Nsps,*)
    integer(i4), intent(in) :: NO_FILT_PTS
    integer(i4), intent(in) :: NO_GEOM
    integer(i4), intent(in) :: NO_GEOPHYS
    integer(i4), intent(in) :: NO_INT_FRQS(*)
    integer(i4), intent(in) :: NO_PFA_CH
    integer(i4), intent(in) :: NO_SPS_TBL(*)
    integer(i4), intent(in) :: SPS_TBL(Nsps,*)
    integer(i4), intent(in) :: PFA_CH(*)
    integer(i4), intent(in) :: NVR(*)
    Real(r4), intent(in) :: Z_PATH(Npath,*), H_PATH(Npath,*)
    Real(r4), intent(in) :: T_PATH(Npath,*), PHI_PATH(Npath,*)
    Real(r4), intent(in) :: DHDZ_PATH(Npath,*)
    Real(r4), intent(in) :: MDB_PRES(*)
    Real(r4), intent(in) :: MDB_TEMP(*)
    Real(r8), intent(in) :: MDB_FREQ(mnf,*)
    Real(r4), intent(in) :: Z_GRID(*), T_GRID(*), H_GRID(*)
    Real(r4), intent(in) :: DH_DT_GRID(Nlvl,*), A_GRID(*), N_GRID(*)
    Real(r4), intent(in) :: V_GRID(*)
    integer(i4), intent(in) :: T_INDX
    Real(r4), intent(in) :: CONV_HTS(*)
    Real(r4), intent(in) :: CONV_PRESS(*)
    Real(r4), intent(in) :: CONV_TEMP(*)
    Real(r4), intent(in) :: SPSFUNC(Nlvl,*)
    Real(r4), intent(out) :: CEN_ANGLE
    Real(r4), intent(in) :: EARTH_REF
    Real(r4), intent(in) :: T_TAN(*)
    Real(r4), intent(in) :: ACC(*)
    Real(r4), intent(in) :: THRESHOLD(*)
    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)
    Real(r4), intent(inout) :: CS(Nlvl,no_t_phi,mnf,*)
    Real(r4), intent(out) :: D2X_DXDT(Nlvl,*)
    Real(r4), intent(out) :: DX_DT(Nlvl,*)
    Real(r8), intent(in) :: F_GRID_FILTER(maxfiltpts,*)
    Real(r8), intent(in) :: F_GRID(maxaitkenpts,*)
    Real(r8), intent(in) :: FILTER_FUNC(maxfiltpts,*)
    Real(r8), intent(in) :: FREQ(*)
    integer(i4), intent(in) :: NO_FREQS_F(*)
    Real(r4), intent(in) :: G_BASIS(mxco,*)
    type(geom_param), intent(in) :: GEOMETRIC(*)
    type(geophys_param), intent(in) :: GEOPHYSIC(*)
    Real(r4), intent(in) :: F_BASIS(mxco,*)
    integer(i4), intent(in) :: NO_COEFFS_F(*)
    Real(r4), intent(in) :: MR_F(mxco,mnp,*)
    integer(i4), intent(in) :: NO_COEFFS_G(*)
    Real(r4), intent(in) :: MR_G(mxco,mnp,*)
    type(pfa_slab) :: PFA_SPECTRUM(6,*)
    Real(r4), intent(in) :: PTG_ANGLE(*)
    Real(r4), intent(in) :: PTG_HTS(*)
    type(limb_press) :: PTG_PRESS
    Real(r4), intent(in) :: REF_CORR(Nlvl,*)
    Real(r4), intent(in) :: AZIM_183, AZIM_205, AZIM_REF
    Real(r4), intent(in) :: C_PITCH, C_ROLL, C_YAW
    Real(r4), intent(in) :: ELEV_183, ELEV_205
    Real(r4), intent(in) :: GEOCSRAD
    Real(r4), intent(in) :: S_PITCH, S_ROLL
    Real(r4), intent(in) :: SP_TMP
    Real(r4), intent(in) :: S_YAW
    integer(i4), intent(in) :: SI
    character(len=*), intent(in) :: AAAP
    character(len=*), intent(in) :: INDIR
    character(len=*), intent(inout) :: KEYS(*)
    character(len=*), intent(in) :: RUNF
    integer(i4), intent(in) :: TIME_I
    integer(i4), intent(in) :: L2PC_LU
    integer(i4), intent(inout) :: RECN
    integer(i4), intent(in) :: MXBIN_RECN
    integer(i4), intent(in) :: LRUN
    integer(i4), intent(inout) :: REC_NOS(*)
    integer(i4), intent(in) :: NO_ATMOS
    integer(i4), intent(in) :: JKEY
    integer(i4), intent(in) :: N_OBS
    integer(i4), intent(in) :: ATMOS_INDEX(*)
    integer(i4), intent(in) :: GEOM_INDEX(*)
    integer(i4), intent(in) :: GEOPHYS_INDEX(*)
    Logical*1, intent(in) :: DO_CONV
    integer(i4), intent(in) :: CH1, CH2
    integer(i4), intent(in) :: SPECT_ATMOS(*)
    integer(i4), intent(in) :: NO_PHI_G(*)
    integer(i4), intent(in) :: NO_PHI_F(*)
    Real(r4), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: NO_PHI_VEC(*)
    integer(i4), intent(in) :: NPATH
    integer(i4), intent(in) :: PATH_BRKPT(3,*)
    integer(i4), intent(in) :: NO_PHI_T
    Real(r4), intent(in) :: T_PHI_BASIS(*)
!   Real(r4), intent(in) :: DH_DT_PATH(Npath,mnp,mxco,*)
    Real(r4), intent(in) :: DH_DT_PATH(Npath,mnp,23,*)
    type(eos_mdb_hdr), intent(in) :: MDB_HDR(*)
    type(eos_mdb_rec), intent(in) :: MDB_REC(max_no_lines,*)
    type(spectro_param), intent(in) :: SPECTROSCOPIC(*)
    integer(i4), intent(out) :: Ier
!
    integer, parameter :: Kgp = 2 * Nlvl * (Ng+1)
!
    real(r4) :: CSE
    integer(i4) :: CYCLE
    real(r4) :: DEL_OPACITY(N2lvl), DELTA(N2lvl,mxco,mnp,Nsps)
    real(r4) :: DH, DT_SCRPT_DNP(N2lvl,mxco,mnp), DTAUDA(Nlvl,mxco)
    real(r4) :: ETA(Nlvl,mxco)
    integer(i4) :: FFT_INDEX(maxfft)
    real(r4) :: FFT_PRESS(maxfft)
    real(r8) :: FRQ
    integer(i4) :: H_I
    real(r4) :: H_TAN
    integer(i4) :: HT_PTR(Nsps)
    integer(i4) :: I, ICH, IC1, IC2, IHI, ILO, INDXL, INDXR
    real(r4) :: I_RAW(Nlvl), I_STAR_ALL(max_no_pointings)
    integer(i4) :: IS
    integer(i4) :: JCH, JTAN
    integer(i4) :: K
    real(r4) :: K_STAR_ALL(msvd,mnp,max_no_pointings)
    real(r4) :: K_STAR_ATMOS(Nlvl,mxco,mnp,Nsps)
    real(r4) :: K_STAR_GEOPHYS(Nlvl,mxco,mnp,maxgeophys)
    integer(i4) :: N_SPS, NO_T
    integer(i4) :: PTG_I
    real(r4) :: R
    real(r8) :: RS ! For ellipse calculation
    real(r4) :: S, SCALE_FACTOR(max_no_pointings)
    integer(i4) :: SPS_COEF_LOOP(2,Nlvl,Nsps), SPS_I, SV_I
    real(r4) :: T, T_SCRIPT(N2lvl), TANX, TAU(N2lvl)
!
! Get central angle for the fourier transform,using 'si' (surface index)
! (the index of the earth surface in: conv_hts  array)
!
    Ier = 0
    cen_angle = ptg_angle(si)
    no_t = no_coeffs_g(t_indx)      ! Number of temperature coeff.
!
! Set up: dx_dt, d2x_dxdt arrays for temperature derivative computations
! (NOTE: These entities has NO PHI dimension, so take the center Phi in dh_dt)
!
!  First: Get table of temperature basis functions
!
    Call get_eta(conv_press,geophysic(t_indx)%basis_peaks,no_conv_hts,  &
   &             no_t,Nlvl,Eta)
!
    dtauda(1:si+1,1:no_t) = 0.0
!
    ptg_i = no_conv_hts - si + 1
    do sv_i = 1, no_t
      Call Lintrp(z_grid,conv_press(si:si+N_lvls-1),dh_dt_grid(1:ptg_i,sv_i),             &
   &                dtauda(si:si+N_lvls-1,sv_i),N_lvls,ptg_i)
    end do
!
    k = 2 * N_lvls
    do sv_i = 1, no_t
!
      is = 1
      r = dh_dt_grid(1,sv_i)
!
      do while (conv_hts(is)  <  0.0)
        t = conv_temp(is)
        dh = conv_hts(is) + RoC
        eta(is,sv_i) = 0.0
        tanx = tan(ptg_angle(is))
        cse = tanx * tanx
        h_tan = cse * cse
        dx_dt(is,sv_i) = r * tanx / dh
        d2x_dxdt(is,sv_i) = (2.0+cse)*r/dh + eta(is,sv_i)/t
        is = is + 1
      end do
!
      do h_i = is, no_conv_hts
        t = conv_temp(h_i)
        dh = conv_hts(h_i) + RoC
        tanx = tan(ptg_angle(h_i))
        cse = tanx * tanx
        h_tan = cse * cse
        r = dtauda(h_i,sv_i)
        dx_dt(h_i,sv_i) = r * tanx / dh
        d2x_dxdt(h_i,sv_i) = (2.0+cse)*r/dh + eta(h_i,sv_i)/t
      end do
!
    end do
!
! Initialize all other arrays used in computations:
!
    t_script(1:N2lvl) = 0.0
!
    do sps_i = 1, Nsps
      do k = 1, mnp
        do is = 1, mxco
          do h_i = 1, Nlvl
            delta(h_i,is,k,sps_i) = 0.0
            delta(h_i+Nlvl,is,k,sps_i) = 0.0
            k_star_atmos(h_i,is,k,sps_i) = 0.0
          end do
        end do
      end do
    end do
!
    k_star_geophys(1:Nlvl,1:mxco,1:mnp,1:maxgeophys) = 0.0
!
! Compute the radiative transfer equation at the state vector linearization
! values and selected channels
!
    jch = 0
    n_sps = no_sps_tbl(band)
!
! Define the begining and ending channels of the inputted band:
!
    is = header1%no_channels_per_band
    ic2 = band * is
    ic1 = ic2  - is + 1
!
! Establish the Sps_Coef_Loop entries (To prevent computing 'Zeros'
! in various integration loops ...).
!
! First, Initialize a species pointer
!
    do sps_i = 1, n_sps
      ht_ptr(sps_i) = no_coeffs_f(sps_tbl(sps_i,band)) + 1
    end do
!
! Second, compute the range as a function of pressure & specie
!
    do h_i = 1, N_lvls
!
      r = z_grid(h_i)
      s = z_grid(h_i+1)
!
      do sps_i = 1, n_sps
!
        is = sps_tbl(sps_i,band)
        sv_i = no_coeffs_f(is)
!
! Use a sequential search method
!
        ptg_i = ht_ptr(sps_i)
        do while (f_basis(ptg_i,is) > r .and. ptg_i > 1)
          ptg_i = ptg_i - 1
        end do
        ht_ptr(sps_i) = ptg_i
!
! Now establish range of non zero coefficient values
!
        do while (s > f_basis(ptg_i,is) .and. ptg_i < sv_i)
          ptg_i = ptg_i + 1
        end do
!
        Sps_Coef_Loop(1,h_i,sps_i) = ht_ptr(sps_i)
        Sps_Coef_Loop(2,h_i,sps_i) = ptg_i
!
      end do
!
    end do
!
! *** Start main loop now:
!
    do ich = ic1, ic2
!
! Clean up the "_star_all"   arrays:
!
      do ptg_i = 1, header1%no_pointings
        i_star_all(ptg_i) = 0.0
        scale_factor(ptg_i) = 1.0           ! Initialize to: 1.0
        k_star_all(1:msvd,1:mnp,ptg_i) = 0.0
      end do
!
      jch = jch + 1
!     Frq = freq(ich)
!
      Frq = 0.0
      do k = 1, 18
        Frq = mdb_freq(jch,k)                ! Temporary code
        if ( frq >= 1.0 ) exit
      end do
!
      if (Frq >= 1.0) then
!
        i_raw(1:no_conv_hts) = 0.0
!
! Now, establish the radiative transfer function and its derivatives:
!
        do ptg_i = 1, no_conv_hts - 1
!
! Find the tangent height in the array and store into preselected grid
!
          k = 2 * N_lvls
          EarthX = .false.
          jtan = n_tan(ptg_i)
          h_tan = conv_hts(ptg_i)
          ht = dble(h_tan)
          Rr = ht + RoC
          ht2 = Rr * Rr
          if (h_tan < -0.01) then
            Rr = Rr / RoC
            cse = earth_ref
            Call Earth_Intersection(Rs)
          else
            cse = 1.0
            Rr = 0.0d0
            Phi_s = Phi_tan
            NPhi_s = NPhi_tan
          end if
!
          IndxR = N_lvls - jtan + 1
          IndxL = N_lvls + jtan
!
!  Compute the appropriate t_script & dt_scrpt_dnp along the integration
!  path of this tanget:
!
          Call do_t_script(N_lvls,Ng,Frq,sp_tmp,t_path(1,ptg_i),IndxR,  &
     &                     IndxL,path_brkpt(1,ptg_i),t_script)
!
          if (geophysic(t_indx)%der_calc(band)) then
!
! Create the dt_scrpt_dnp arrays for all coefficients:
!
            do i = 1, no_phi_t
              do k = 1, no_t
                Call d_t_script_dtnp(Frq,g_basis(1,t_indx),t_phi_basis, &
     &               t_path(1,ptg_i),z_path(1,ptg_i),phi_path(1,ptg_i), &
     &               N_lvls,Ng,path_brkpt(1,ptg_i),IndxR,IndxL,k,i,no_t,&
     &               no_phi_t,dt_scrpt_dnp(1,k,i))
              end do
            end do
!
          end if
!
          k = ptg_i
          Call get_delta(z_path(1,k),t_path(1,k),h_path(1,k),           &
     &         phi_path(1,k),dHdz_path(1,k),N_lvls,cs,Frq,n_sps,        &
     &         no_coeffs_f,sps_tbl(1,band),mxco,Nlvl,no_phi_t,          &
     &         f_basis,Sps_Coef_Loop,ref_corr(1,k),mdb_pres,mdb_temp,   &
     &         mdb_freq,mnp,mnf,no_freqs_f,no_phi_f,phi_basis_f,        &
     &         IndxR,IndxL,path_brkpt(1,k),delta,Ier)
          if (Ier /= 0) Return
!
! Initialize the tau & del_opacity arrays:
!
          do k = 1, N2lvl
            tau(k) = 0.0
            del_opacity(k) = 0.0
          end do
!
          Call zopacity(mr_f,sps_tbl(1,band),no_coeffs_f,n_sps,         &
     &         N_lvls,N2lvl,mxco,mnp,no_phi_f,delta,del_opacity,        &
     &         IndxR,IndxL)
!
          Call Scrt_dn(t_script,N_lvls,cse,del_opacity,tau,i_raw(ptg_i),&
     &                 IndxR,IndxL,ilo,ihi)
!
! Compute atmosperic derivatives for this channel
!
          Call zatmos_deriv(atmospheric,ptg_i,sps_tbl(1,band),n_sps,    &
     &         band,no_coeffs_f,IndxR,IndxL,delta,t_script,tau,ilo,     &
     &         ihi,no_phi_f,k_star_atmos)
!
          k = ich - ic1 + 1
!         if (k <= 1) then         ! DEBUG, T_Deriv. for 1 channel(s) only
          if (k <= -1) then        ! DEBUG, Skip derivatives altogether ..
!
! Compute temperature derivative for this channel (if requested)
!
            if (geophysic(t_indx)%der_calc(band)) then
!
              k = ptg_i
              Call temperature_deriv(geophysic,t_indx,band,cs,Frq,        &
     &             f_basis,z_path(1,k),t_path(1,k),h_path(1,k),           &
     &             phi_path(1,k),dHdz_path(1,k),dh_dt_path(1,1,1,k),Kgp,  &
     &             mnp,mdb_freq,mdb_pres,mdb_temp,mr_f,mr_g,              &
     &             no_coeffs_f,no_freqs_f,no_phi_f,no_phi_t,no_t,N_lvls,  &
     &             n_sps,phi_basis_f,k,ref_corr(1,k),sps_tbl,t_phi_basis, &
     &             tau,t_script,dt_scrpt_dnp,mnf,ilo,ihi,IndxR,IndxL,     &
     &             path_brkpt(1,k),K_star_geophys,Ier)
              if (Ier /= 0) Return
!
            end if
          end if
!
          k = 1
!
        end do                         ! On ptg_i
!
! Complete all arrays for the last location (no_conv_hts):
!
        k = no_conv_hts - 1
        Call complete_i_k(k,no_conv_hts,no_geophys,i_raw,k_star_atmos,  &
     &                    k_star_geophys)
!
        if (ich == 61) then
          Print *
          is = no_conv_hts - si + 1
          write(*,'(''conv_press'',a1,i2,''n'')') char(92),is
          write(*,'(6(1x,f10.6))') (conv_press(k),k=si,no_conv_hts)
          write(*,*)
          write(*,'(''UN-CONVOLVED radiances, channel:'',i3)') ich
          write(*,'(''z_star'',a1,i2,''n'')') char(92),is
          write(*,'(6(1x,1pg12.6))') (i_raw(k),k=si,no_conv_hts)
        end if
!
        if (.not. do_conv) then
!
! Put the derivatives and radiances in a storage array for write-out
! NO ANTENNA CONVOLUTION IS DONE IN HERE ....
!
          Ier = Ich
          Call i_k_star_no_conv(header1,ptg_press,geometric,geophysic,  &
     &         atmospheric,no_geom,no_geophys,n_sps,sps_tbl(1,band),    &
     &         conv_press,conv_hts,ptg_angle,band,i_raw,k_star_geophys, &
     &         k_star_atmos,a_grid,t_grid,z_grid,N_lvls,no_conv_hts,    &
     &         c_yaw,s_yaw,c_roll,s_roll,c_pitch,s_pitch,elev_183,      &
     &         elev_205,azim_183,azim_205,azim_ref,geocsrad,i_star_all, &
     &         k_star_all,no_phi_g,no_phi_f,Ier)
          if (Ier /= 0) Return
!
        else
!
!
! Put the derivatives and radiances in a storage array for write-out
! The antenna convolution is done in here also
!
          Ier = Ich
          Call i_and_k_star(header1,ptg_press,geometric,geophysic,      &
       &       atmospheric,no_geom,no_geophys,n_sps,sps_tbl(1:,band),   &
       &       conv_press,conv_hts,conv_temp,ptg_hts,ptg_angle,dx_dt,   &
       &       d2x_dxdt,band,cen_angle,fft_pts,fft_press,i_raw,         &
       &       InDir,Aaap,k_star_geophys,k_star_atmos,a_grid,t_grid,    &
       &       z_grid,N_lvls,no_conv_hts,c_yaw,s_yaw,c_roll,s_roll,     &
       &       c_pitch,s_pitch,elev_183,elev_205,azim_183,azim_205,     &
       &       azim_ref,geocsrad,i_star_all,k_star_all,fft_index,       &
       &       no_phi_g,no_phi_f,Ier)
          if (Ier /= 0) Return
        end if
!
        if (ich == 61) then
          Print *
          is = ptg_press%no_lin_values
          write(*,'(''pointing_press'',a1,i2,''n'')') char(92),is
          write(*,'(6(1x,f10.6))') (ptg_press%lin_val(k),k=1,is)
          write(*,*)
          write(*,'(''CONVOLVED radiances, channel:'',i3)') ich
          write(*,'(''c_star'',a1,i2,''n'')') char(92),is
          write(*,'(6(1x,1pg12.6))') (i_star_all(k),k=1,is)
        end if
!
        if ( no_pfa_ch >= 1 .and. &
           (ich < pfa_ch(1) .or. ich > pfa_ch(no_pfa_ch)) ) then
!
          is = -1
          do k = 1, no_pfa_ch
            if (ich == pfa_ch(k)) then
              is = k
              exit
            end if
          end do
!
          if (is >= 1) then               ! Not a PFA channel
!
! Run the pfa forward model for this channel:
!
            Call pfa_fwd_mdl(header1,N_lvls,fft_pts,no_conv_hts,n_tan,      &
         &           band,no_filt_pts,no_int_frqs,no_pfa_ch,sps_tbl,pfa_ch, &
         &           z_path,h_path,t_path,phi_path,dHdz_path,mdb_pres,      &
         &           mdb_temp,mdb_freq,z_grid,t_grid,a_grid,conv_hts,       &
         &           cen_angle,earth_ref,Threshold,f_grid_filter,f_grid,    &
         &           filter_func,no_freqs_f,f_basis,no_coeffs_f,mr_f,       &
         &           pfa_spectrum,ptg_angle,ptg_hts,ptg_press,ref_corr,     &
         &           azim_183,azim_205,azim_ref,c_pitch,c_roll,c_yaw,       &
         &           elev_183,elev_205,s_pitch,s_roll,sp_tmp,s_yaw,si,Aaap, &
         &           InDir,do_conv,no_phi_f,phi_basis_f,Npath,path_brkpt,   &
         &           no_phi_t,mdb_hdr,mdb_rec,is,n_sps,Sps_Coef_Loop,       &
         &           spect_atmos,spectroscopic,scale_factor,i_star_all,     &
         &           k_star_all,Ier)
            if (Ier /= 0) Return
!
            if (ich == 61) then
              Print *
              is = ptg_press%no_lin_values
              write(*,'(''CONVOLVED PFA radiances, channel:'',i3)') ich
              write(*,'(''i_star'',a1,i2,''n'')') char(92),is
              write(*,'(6(1x,1pg12.6))') (i_star_all(k),k=1,is)
            end if
          end if
        end if ! no_pfa_ch >= 1 .and. ...
      end if ! Frq >= 1.0
!
! Dump the x & i records:
!
      Call write_x_i_records(jch,time_i,l2pc_lu,recn,mxbin_recn,lrun,   &
   &       rec_nos,band,no_geom,no_atmos,no_geophys,mr_g,mr_f,Keys,     &
   &       runf,ptg_press,geometric,geophysic,atmospheric,header1,      &
   &       header2,jkey,i_star_all,atmos_index,geom_index,              &
   &       geophys_index,no_phi_g,no_phi_f,Ier)
      if (Ier /= 0) Return
!
      Call dump_k_records(jch,time_i,band,scale_factor,header1,         &
   &       header2,header3,no_phi_vec,k_star_all,cycle,Ier)
      if (Ier /= 0) Return
!
    end do                                     ! On ich
!
    Call write_k_records(time_i,l2pc_lu,recn,mxbin_recn,lrun,rec_nos,   &
   &     band,Keys,runf,header1,jkey,cycle,Ier)
!
    Return
  End Subroutine FWD_MDL
end module FWD_MDL_M

! $Log$
