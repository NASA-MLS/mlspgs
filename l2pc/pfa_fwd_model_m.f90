module PFA_FWD_MDL_M
  use AITKEN_INT_M, only: AITKEN_INT
  use CREATE_BETA_COEFF_M, only: CREATE_BETA_COEFF
  use CREATE_PFA_BETA_COEFF_M, only: CREATE_PFA_BETA_COEFF
  use DO_T_SCRIPT_M, only: DO_T_SCRIPT
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE, only: EARTHX, HT, HT2, NPHI_TAN, NPHI_S, PHI_S, PHI_TAN, &
                     ROC, RR
  use EOS_MDB, only: EOS_MDB_HDR,EOS_MDB_REC,MAX_NO_LINES,CS_MNF => max_freq
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  use GET_PFA_DELTA_M, only: GET_PFA_DELTA
  use GL6P, only: NG
  use D_HUNT_M, only: HUNT
  use HYDROSTATIC_INTRP, only: GET_PRESSURES
  use L2PC_FILE_PARAMETERS, only: DEG2RAD, MAX_NO_SV_COMPONENTS, &
                                  MSVD => max_no_sv_derivatives, &
                                  MXCO => max_no_elmnts_per_sv_component
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PC_PFA_STRUCTURES, only: LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS, &
                                 MAXPFACH, PFA_SLAB, SPECTRO_PARAM
  use L2PCDim, only: MAXFFT, MNP => max_no_phi, NLVL, NPTG, NSPS, N2LVL
  use MLSCommon, only: I4, R4, R8
  use D_LINTRP_M, only: LINTRP
  use SCRT_DN_M, only: SCRT_DN
  use DCSPLINE_DER_M, only: CSPLINE_DER
  use SPECTRO_DERIVATIVE_M, only: SPECTRO_DERIVATIVE
  use STRINGS, only: STRLWR
  use ZOPACITY_M, only: ZOPACITY
  implicit NONE
  private
  public :: PFA_FWD_MDL
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This is the PFA forward model, called by channel:
!
  Subroutine PFA_FWD_MDL ( header1, N_lvls, fft_pts, no_conv_hts,        &
 &           n_tan, band, no_filt_pts, no_int_frqs, no_pfa_ch, sps_tbl,  &
 &           pfa_ch, z_path, h_path, t_path, phi_path, dHdz_path,        &
 &           mdb_pres, mdb_temp, mdb_freq, z_grid, t_grid, a_grid,       &
 &           conv_hts, cen_angle, earth_ref, Threshold, f_grid_filter,   &
 &           f_grid, filter_func, no_freqs_f, f_basis, no_coeffs_f,      &
 &           mr_f, pfa_spectrum, ptg_angle, ptg_hts, ptg_press,          &
 &           ref_corr, azim_183, azim_205, azim_ref, c_pitch, c_roll,    &
 &           c_yaw, elev_183, elev_205, s_pitch, s_roll, sp_tmp, s_yaw,  &
 &           si, Aaap, InDir, do_conv, no_phi_f, phi_basis_f, Npath,     &
 &           path_brkpt, no_phi_t, mdb_hdr, mdb_rec, jch, n_sps,         &
 &           Sps_Coef_Loop, spect_atmos, spectroscopic, scale_factor,    &
 &           i_star_all, k_star_all, Ier )
!
    integer(i4), intent(in) :: NPATH
    integer(i4), intent(in) :: N_LVLS
    integer(i4), intent(in) :: FFT_PTS
    integer(i4), intent(in) :: NO_CONV_HTS
    integer(i4), intent(in) :: NO_PHI_T
    integer(i4), intent(in) :: BAND
    integer(i4), intent(in) :: NO_PFA_CH
    integer(i4), intent(in) :: NO_FILT_PTS
    integer(i4), intent(in) :: NO_INT_FRQS(*)
    integer(i4), intent(in) :: SPS_TBL(Nsps,*)
    integer(i4), intent(in) :: PFA_CH(*)
    integer(i4), intent(in) :: N_TAN(*)
!
    type(l2pc_header_one), intent(in) :: HEADER1
    Real(r8), intent(in) :: Z_PATH(Npath,*), H_PATH(Npath,*),           &
   &                        T_PATH(Npath,*), PHI_PATH(Npath,*)
    Real(r4), intent(in) :: DHDZ_PATH(Npath,*)
    Real(r8), intent(in) :: MDB_PRES(*)
    Real(r8), intent(in) :: MDB_TEMP(*)
    Real(r8), intent(in) :: MDB_FREQ(cs_mnf,*)
    Real(r8), intent(in) :: Z_GRID(*), T_GRID(*), A_GRID(*)
    Real(r8), intent(in) :: CONV_HTS(*)
    Real(r8), intent(in) :: CEN_ANGLE
    Real(r8), intent(in) :: EARTH_REF
    Real(r8), intent(in) :: THRESHOLD(*)
    Real(r8), intent(in) :: F_GRID_FILTER(maxfiltpts,*)
    Real(r8), intent(in) :: F_GRID(maxaitkenpts,*)
    Real(r8), intent(in) :: FILTER_FUNC(maxfiltpts,*)
    integer(i4), intent(in) :: NO_FREQS_F(*)
    Real(r8), intent(in) :: F_BASIS(mxco,*)
    integer(i4), intent(in) :: NO_COEFFS_F(*)
    Real(r8), intent(in) :: MR_F(mxco,mnp,*)
    type(pfa_slab), intent(in) :: PFA_SPECTRUM(6,*)
    Real(r8), intent(in) :: PTG_ANGLE(*)
    Real(r8), intent(in) :: PTG_HTS(*)
    type(limb_press), intent(in) :: PTG_PRESS
    Real(r8), intent(in) :: REF_CORR(Nlvl,*)
    Real(r8), intent(in) :: AZIM_183, AZIM_205, AZIM_REF
    Real(r8), intent(in) :: C_PITCH, C_ROLL, C_YAW
    Real(r8), intent(in) :: ELEV_183, ELEV_205
    Real(r8), intent(in) :: S_PITCH, S_ROLL
    Real(r8), intent(in) :: SP_TMP
    Real(r8), intent(in) :: S_YAW
    integer(i4), intent(in) :: SI
    character(len=*), intent(in) :: AAAP
    character(len=*), intent(in) :: INDIR
    Logical*1, intent(in) :: DO_CONV
    integer(i4), intent(in) :: NO_PHI_F(*)
    Real(r8), intent(in) :: PHI_BASIS_F(mnp,*)
    integer(i4), intent(in) :: PATH_BRKPT(3,*)
    type(eos_mdb_hdr), intent(in) :: MDB_HDR(*)
    type(eos_mdb_rec), intent(in) :: MDB_REC(max_no_lines,*)
    integer(i4), intent(in) :: JCH
    integer(i4), intent(in) :: N_SPS
    integer(i4), intent(in) :: SPS_COEF_LOOP(2,Nlvl,*)
    integer(i4), intent(in) :: SPECT_ATMOS(*)
    type(spectro_param), intent(in) :: SPECTROSCOPIC(*)
    Real(r8), intent(out) :: SCALE_FACTOR(*)
    Real(r8), intent(inout) :: I_STAR_ALL(*)
    Real(r4), intent(out) :: K_STAR_ALL(msvd,mnp,*)
    integer(i4), intent(out) :: IER
!
    character(len=24) :: AX
    real(r8) :: A, A_63, A2B2, AZIM_ANGLE, B, B2
!   real(r8) :: beta_t_power(N2lvl,maxaitkenpts,maxpfach,Nsps,Nptg)
    real(r8) :: beta_t_power(N2lvl,maxaitkenpts,2,2,Nptg)
    character(len=1) :: CA
    real(r8) :: C_D_ELEV, C_E, E_Z, ELEV_OFFST, CSE
    integer(i4) :: COMP_NDX, CON_TYPE
    real(r8) :: DEL_OPCTY(N2lvl), DELTA(N2lvl,mxco,mnp,Nsps), DRADZ(Nptg)
    real(r8) :: FFT_ANGLES(maxfft), FFT_PRESS(maxfft)
    integer(i4) :: FFT_INDEX(maxfft)
    real(r8) :: FRQ
    real(r8) :: H_TAN
    integer(i4) :: IHI, ILO, INDXL, INDXR, IP, IZ
    integer(i4) :: J, J4, JF, JTAN, JZ
    integer(i4) :: K
!   real(r4) :: k_star_spect_dw(Nptg,mxco,mnp,Nsps),                   &
!  &       k_star_spect_dn(Nptg,mxco,mnp,Nsps),                        &
!  &       k_star_spect_dnu(Nptg,mxco,mnp,Nsps)
!
    real(r4) :: k_star_spect_dw(Nptg,mxco,mnp,2),                      &
   &       k_star_spect_dn(Nptg,mxco,mnp,2),                           &
   &       k_star_spect_dnu(Nptg,mxco,mnp,2)
    integer(i4) :: KCONV, KTR
    integer(i4) :: M
    integer(i4) :: NF, NTR
!   real(r8) :: pfa_beta_coeff(N2lvl,maxaitkenpts,maxpfach,Nsps,Nptg)
    real(r8) :: pfa_beta_coeff(N2lvl,maxaitkenpts,2,2,Nptg)
!   real(r8) :: pfa_dbeta_dn(N2lvl,maxaitkenpts,maxpfach,Nsps,Nptg)
    real(r8) :: pfa_dbeta_dn(N2lvl,maxaitkenpts,2,2,Nptg)
!   real(r8) :: pfa_dbeta_dnu(N2lvl,maxaitkenpts,maxpfach,Nsps,Nptg)
    real(r8) :: pfa_dbeta_dnu(N2lvl,maxaitkenpts,2,2,Nptg)
!   real(r8) :: pfa_dbeta_dw(N2lvl,maxaitkenpts,maxpfach,Nsps,Nptg)
    real(r8) :: pfa_dbeta_dw(N2lvl,maxaitkenpts,2,2,Nptg)
    real(r8) :: PFA_RAD(Nptg)
    integer(i4) :: PTG_I
    real(r8) :: Q
    real(r8) :: R, RAD(maxaitkenpts)
    real(r8) :: Rad_dn(maxaitkenpts,mxco,mnp,Nsps)
    real(r8) :: Rad_dnu(maxaitkenpts,mxco,mnp,Nsps)
    real(r8) :: Rad_dw(maxaitkenpts,mxco,mnp,Nsps)
    real(r8) :: ROOT
    real(r8) :: RS ! For ellipse calculation
    real(r8) :: S_D_ELEV, S_E
    integer(i4) :: S_NP, S_NZ, SA, SHFT, SPECTAG, SV_C, SV_ELMNT
    real(r8) :: T_SCRIPT(N2lvl), TAU(N2lvl), TMPRAD(Nptg), PtP(Nlvl)
    real(r8) :: VRAD(maxfft)
!
! -----     Statement functions     ------------------------------------
!
! These statement functions compensate for the use of ACOSD, COSD and SIND,
! which are intrinsic in some compilers, but are not part of any standard.
!
    real(r8) :: COSD, SIND, X
    cosd(x) = Cos(x*DEG2RAD)
    sind(x) = Sin(x*DEG2RAD)
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
!
! Clear (zero out) delta, beta_coeff arrays and k_star_spect arrays:
!
    delta = 0.0
    beta_t_power = 0.0
    pfa_dbeta_dw = 0.0
    pfa_dbeta_dn = 0.0
    pfa_dbeta_dnu = 0.0
    pfa_beta_coeff = 0.0
    k_star_spect_dw = 0.0
    k_star_spect_dn = 0.0
    k_star_spect_dnu = 0.0
!
!  Create the beta_coeff for each conv_hts over the Aitken freq. grid:
!
    j4 = 4 * No_Int_Frqs(jch) + 1
    do k = 1, No_Conv_Hts
      jtan = n_tan(k)
      IndxR = N_lvls - jtan + 1
      IndxL = N_lvls + jtan
      Call Create_pfa_beta_coeff(jch,band,IndxL,IndxR,N_lvls,           &
   &       sps_tbl(1:,band),z_path(1:,k),t_path(1:,k),path_brkpt(1:,k), &
   &       j4,f_grid(1:,jch),pfa_spectrum,mdb_hdr,mdb_rec,              &
   &       pfa_beta_coeff(1:,1:,1:,1:,k),pfa_dbeta_dw(1:,1:,1:,1:,k),   &
   &       pfa_dbeta_dn(1:,1:,1:,1:,k),pfa_dbeta_dnu(1:,1:,1:,1:,k),    &
   &       beta_t_power(1:,1:,1:,1:,k))
    end do
!
    shft = si-1
    if (do_conv) shft = 0
    kconv = no_conv_hts - shft
!
    pfa_rad(1:Nptg) = 0.0
    k = ptg_press%no_lin_values
    PtP(1:k) = dble(Ptg_Press%lin_val(1:k))
!
! Now, establish the radiative transfer function and its derivatives:
!
    do ptg_i = 1, No_Conv_Hts - 1
!
! Find the tangent height in the array and store into preselected grid
!
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
      Rad(1:j4) = 0.0
!
! *** Start main frequency loop now:
!
      k = ptg_i
      do jf = 1, j4
!
        Frq = f_grid(jf,jch)
!
!  Compute the appropriate t_script along the integration path of
!  this tanget:
!
        Call do_t_script(N_lvls,Ng,Frq,sp_tmp,t_path(1:,k),IndxR,IndxL, &
   &                     path_brkpt(1:,k),t_script)
!
        Call get_pfa_delta(z_path(1:,k),t_path(1:,k),h_path(1:,k),       &
   &         phi_path(1:,k),dHdz_path(1:,k),N_lvls,jf,jch,n_sps,         &
   &         pfa_beta_coeff(1:,1:,1:,1:,k),no_coeffs_f,sps_tbl(1:,band), &
   &         mxco,Nlvl,N2lvl,maxaitkenpts,maxpfach,Nsps,no_phi_t,        &
   &         f_basis,Sps_Coef_Loop,ref_corr(1:,k),mnp,no_freqs_f,        &
   &         no_phi_f,phi_basis_f,IndxR,IndxL,path_brkpt(1:,k),          &
   &         beta_t_power(1:,1:,1:,1:,k),delta,Ier)
        if (Ier /= 0) Return
!
! Initialize the tau & del_opcty arrays:
!
        tau(1:N2lvl) = 0.0
        del_opcty(1:N2lvl) = 0.0
!
        Call zopacity(mr_f,sps_tbl(1:,band),no_coeffs_f,n_sps,N_lvls, &
   &         N2lvl,mxco,mnp,no_phi_f,delta,del_opcty,IndxR,IndxL)
!
        Call Scrt_dn(t_script,N_lvls,cse,del_opcty,tau,Rad(jf),IndxR, &
   &                 IndxL,ilo,ihi)
!
        del_opcty(1:N2lvl) = 0.0
!
!  Get the dI/d(Spectral parameters) (w, n & nu0) FOR EACH SPECIE,
!  and for each Zeta/Phi coefficient.
!
        do jz = 1, n_sps
!
          nf = sps_tbl(jz,band)
          sa = spect_atmos(nf)
          if (sa >= 1) then
!
            Spectag = spectroscopic(sa)%spectag
!
            do while (spectroscopic(sa)%spectag == Spectag)
!
              CA = spectroscopic(sa)%type
              s_np = spectroscopic(sa)%no_phi_values
              s_nz = spectroscopic(sa)%no_zeta_values
!
              if (CA == 'W') then
                Call spectro_derivative(z_path(1:,k),t_path(1:,k),       &
   &                 h_path(1:,k),phi_path(1:,k),dHdz_path(1:,k),N_lvls, &
   &                 jf,jch,nf,no_coeffs_f,mxco,f_basis,mr_f,            &
   &                 ref_corr(1:,k),mnp,no_phi_f,phi_basis_f,IndxR,      &
   &                 IndxL,path_brkpt(1:,k),beta_t_power(1:,1:,1:,1:,k), &
   &                 pfa_dbeta_dw(1:,1:,1:,1:,k),tau,t_script,del_opcty, &
   &                 ilo,ihi,spectroscopic(sa),Rad_dw(1:,1:,1:,nf),Ier)
!
              else if (CA == 'N') then
                Call spectro_derivative(z_path(1:,k),t_path(1:,k),       &
   &                 h_path(1:,k),phi_path(1:,k),dHdz_path(1:,k),N_lvls, &
   &                 jf,jch,nf,no_coeffs_f,mxco,f_basis,mr_f,            &
   &                 ref_corr(1:,k),mnp,no_phi_f,phi_basis_f,IndxR,      &
   &                 IndxL,path_brkpt(1:,k),beta_t_power(1:,1:,1:,1:,k), &
   &                 pfa_dbeta_dn(1:,1:,1:,1:,k),tau,t_script,del_opcty, &
   &                 ilo,ihi,spectroscopic(sa),Rad_dn(1:,1:,1:,nf),Ier)
!
              else if (CA == 'V') then
                Call spectro_derivative(z_path(1:,k),t_path(1:,k),        &
   &                 h_path(1:,k),phi_path(1:,k),dHdz_path(1:,k),N_lvls,  &
   &                 jf,jch,nf,no_coeffs_f,mxco,f_basis,mr_f,             &
   &                 ref_corr(1:,k),mnp,no_phi_f,phi_basis_f,IndxR,       &
   &                 IndxL,path_brkpt(1:,k),beta_t_power(1:,1:,1:,1:,k),  &
   &                 pfa_dbeta_dnu(1:,1:,1:,1:,k),tau,t_script,del_opcty, &
   &                 ilo,ihi,spectroscopic(sa),Rad_dnu(1:,1:,1:,nf),Ier)
!
              end if
!
              if (Ier /= 0) Return
!
              sa = sa + 1
!
            end do
          end if
        end do                       ! On jz (Specie loop)
      end do                         ! On jf (Frequency loop)
!
! Assemble the radiances:
!
      Pfa_rad(ptg_i) = Aitken_int(f_grid(1:,jch),f_grid_filter(1:,jch), &
                       filter_func(1:,jch),Rad,j4,No_filt_pts,Ier)
      if (Ier /= 0) Return
!
! Assemble the derivatives, per specie:
!
      do jz = 1, n_sps
!
        nf = sps_tbl(jz,band)
        sa = spect_atmos(nf)
        if (sa >= 1) then
!
          Spectag = spectroscopic(sa)%spectag
          do while (spectroscopic(sa)%spectag == Spectag)
!
            CA = spectroscopic(sa)%type
            s_np = spectroscopic(sa)%no_phi_values
            s_nz = spectroscopic(sa)%no_zeta_values
            if ((s_nz > 0) .and. (s_np > 0)) then
!
              if (CA == 'W') then
                do iz = 1, s_nz
                  do ip = 1, s_np
                    a = Aitken_int(f_grid(1:,jch),f_grid_filter(1:,jch), &
                           filter_func(1:,jch),Rad_dw(1:,iz,ip,nf),j4,   &
                           No_filt_pts,Ier)
                    if (Ier /= 0) Return
                    k_star_spect_dw(ptg_i,iz,ip,nf) = a
                  end do
                end do
!
              else if (CA == 'N') then
                do iz = 1, s_nz
                  do ip = 1, s_np
                    a = Aitken_int(f_grid(1:,jch),f_grid_filter(1:,jch), &
                           filter_func(1:,jch),Rad_dn(1:,iz,ip,nf),j4,   &
                           No_filt_pts,Ier)
                    if (Ier /= 0) Return
                    k_star_spect_dn(ptg_i,iz,ip,nf) = a
                  end do
                end do
!
              else if (CA == 'V') then
                do iz = 1, s_nz
                  do ip = 1, s_np
                    a = Aitken_int(f_grid(1:,jch),f_grid_filter(1:,jch), &
                           filter_func(1:,jch),Rad_dnu(1:,iz,ip,nf),j4,  &
                           No_filt_pts,Ier)
                    if (Ier /= 0) Return
                    k_star_spect_dnu(ptg_i,iz,ip,nf) = a
                  end do
                end do
!
              end if
!
            end if
!
            sa = sa + 1
!
          end do
        end if
      end do                       ! On jz (Specie loop)
    end do                         ! On ptg_i
!
! Complete all arrays for the last location (No_Conv_Hts):
!
    k = No_Conv_Hts - 1
    Pfa_Rad(No_Conv_Hts) = Pfa_Rad(k)
!
    do jz = 1, n_sps
!
      nf = sps_tbl(jz,band)
      sa = spect_atmos(nf)
      if (sa >= 1) then
!
        Spectag = spectroscopic(sa)%spectag
        do while (spectroscopic(sa)%spectag == Spectag)
!
          CA = spectroscopic(sa)%type
          s_np = spectroscopic(sa)%no_phi_values
          s_nz = spectroscopic(sa)%no_zeta_values
!
          if ((s_nz > 0) .and. (s_np > 0)) then
!
            if (CA == 'W') then
              do iz = 1, s_nz
                do ip = 1, s_np
                  k_star_spect_dw(No_Conv_Hts,iz,ip,nf) =              &
   &                                      k_star_spect_dw(k,iz,ip,nf)
                end do
              end do
!
            else if (CA == 'N') then
              do iz = 1, s_nz
                do ip = 1, s_np
                  k_star_spect_dn(No_Conv_Hts,iz,ip,nf) =              &
   &                                      k_star_spect_dn(k,iz,ip,nf)
                end do
              end do
!
            else if (CA == 'V') then
              do iz = 1, s_nz
                do ip = 1, s_np
                  k_star_spect_dnu(No_Conv_Hts,iz,ip,nf) =             &
   &                                      k_star_spect_dnu(k,iz,ip,nf)
                end do
              end do
!
            end if
!
          end if
!
          sa = sa + 1
!
        end do
      end if
    end do
!
    con_type = 1
    do ptg_i = 1, no_conv_hts
      VRad(ptg_i) = Pfa_Rad(ptg_i)
      fft_angles(ptg_i) = ptg_angle(ptg_i)
    end do
!
    if (do_conv) then
!
! Now compute the convolution of the combined radiances
!
      Call fov_convolve(fft_angles,VRad,cen_angle,con_type,no_conv_hts, &
     &                  band,fft_pts,InDir,Aaap,Ier)
      if (Ier /= 0) Return
!
    end if
!
! Determine radiometer 1 = 63,2 = 205,3 = 183
!
    if (band  >=  5) then                    ! 183
      azim_angle = azim_ref + azim_183
      elev_offst = elev_183
    else if (band  >=  2) then               ! 205
      azim_angle = azim_ref + azim_205
      elev_offst = elev_205
    else                                     ! 63 / pointing reference
      azim_angle = azim_ref
      elev_offst = 0.0
    end if
!
! At this stage we have I vs X,now we want I vs X^63 our common radiometer
! grid.
!
    e_z = s_pitch * c_yaw
    b = c_pitch * c_roll - s_pitch * s_yaw * s_roll
    r = s_pitch * s_yaw * c_roll + c_pitch * s_roll
    a = -cosd(azim_angle) * e_z + sind(azim_angle) * r
!
    b2 = b * b
    a2b2 = a * a + b2
!
! These are the coresponding reference values
!
    a_63 = -cosd(azim_ref) * e_z + sind(azim_ref) * r
!
    ntr = no_conv_hts
    if (do_conv) ntr = 2 ** fft_pts
!
    do ptg_i = 1, ntr
!
! Compute equivalent epsilon for the band
!
      e_z = cos(fft_angles(ptg_i))
!
      q = a2b2 - e_z * e_z
      if (q < 0.0) then
        Ier = 1
        Print *,'** Error in subroutine: pfa_fwd_mdl'
        Print *,'   Taking Sqrt(Arg) for Arg < 0'
        Return
      end if
!
      root = sqrt(q)
      c_e = (e_z * a + b * root) / a2b2
      s_e = (e_z * b - a * root) / a2b2
!
! Adjust the 'fft_angles' accordingly
!
      c_d_elev = cosd(elev_offst)
      s_d_elev = sind(elev_offst)
      q = c_e * (a_63 * c_d_elev - b * s_d_elev) +  &
   &      s_e * (a_63 * s_d_elev + b * c_d_elev)
      if (abs(q) > 1.0) then
        Ier = 1
        Print *,'** Error in subroutine: pfa_fwd_mdl'
        Print *,'   Taking Acos(Arg) for abs(Arg) > 1'
        Return
      end if
!
      fft_angles(ptg_i) = Acos(q)
!
    end do
!
!  Get 'ntr' pressures associated with the fft_angles:
!
    Call get_pressures('a',a_grid,t_grid,z_grid,N_lvls,fft_angles, &
   &                   fft_press,ntr,Ier)
    if (Ier /= 0) Return
!
!  For the NO CONVOLUTION case, make sure the convolution heights BELOW
!  surface do not go into any interpolation routine. Cut the pressure
!  array artificially by making it monotonically DECREASING for all the
!  indecies BELOW surface. The monotonicity check code below will then cut
!  those out.
!
    if(.not.do_conv) then
      do ptg_i = 1, si-1
        fft_press(ptg_i) = 111.0-real(ptg_i)
      end do
    endif
!
! Make sure the Pressures array is MONOTONICALY increasing:
!
    jz = 1
    do while (jz < ntr .and. fft_press(jz) >= fft_press(jz+1))
      jz = jz + 1
    end do
!
    ktr = 1
    fft_index(ktr) = jz
    VRad(ktr) = VRad(jz)
    fft_press(ktr) = fft_press(jz)
!
    do ptg_i = jz+1, ntr
      r = fft_press(ptg_i)
      if (r > fft_press(ktr)) then
        ktr = ktr + 1
        fft_press(ktr) = r
        fft_index(ktr) = ptg_i
        VRad(ktr) = VRad(ptg_i)
      end if
    end do
!
    if (ktr == ntr) fft_index(1) = -2
!
! Interpolate the output values and store the radiances in: TmpRad
! and its derivative w.r.t. pressure in: dRadz
!
    j = ptg_press%no_lin_values
    Call Cspline_der(fft_press,PtP,VRad,TmpRad,dRadz,ktr,j)
!
! Insert "improved" radiance into i_star and get scaling factor for
! the derivatives:
!
    do ptg_i = 1, j
      r =  i_star_all(ptg_i) / TmpRad(ptg_i)
      scale_factor(ptg_i) = r
      i_star_all(ptg_i) = TmpRad(ptg_i)
    end do
!
! Find out if user wants pointing derivatives
!
    if (ptg_press%der_calc(band)) then
!
! Derivatives wanted,get index location in 'k_star_all' and store in the
! derivative array.
!
! ** NOTE: DO NOT scale the derivative here, it is being done in
!          subroutine: dump_k_records  !!!
!
      sv_c = 1
      do while (ptg_press%name /= header1%sv_components(sv_c)  .and.     &
   &           sv_c <  header1%no_sv_components)
        sv_c = sv_c + 1
      end do
!
      jz = header1%sv_component_first_elmnt_index(sv_c)
!
! Copy the derivatives into the array:
!
      do nf = 1, mnp
        k_star_all(jz,nf,1:j) = dRadz(1:j)
      end do
!
    end if
!
! ***************** spectroscopic derivatives ******************
!
    do jz = 1, n_sps
!
! n_sps is the number of species included for this channel
! check to determine if derivative is desired for that species
!
      nf = sps_tbl(jz,band)
      sa = spect_atmos(nf)
      if (sa >= 1) then
!
! Derivatives needed continue to process
! Find index that matches name in header
!
        comp_ndx = 1
        iz = max_no_sv_components
        do while (spectroscopic(sa)%name /=                          &
     &            header1%sv_components(comp_ndx) .and. comp_ndx < iz)
          comp_ndx = comp_ndx + 1
        end do
!
        m = header1%sv_component_first_elmnt_index(comp_ndx)
        iz = m + spectroscopic(sa)%no_zeta_values - 1
        if (iz > msvd) then
          Ier = 1
          print 941,' variable: sv_elmnt',msvd,iz,msvd
          Return
        end if
!
        Spectag = spectroscopic(sa)%spectag
!
        do while (spectroscopic(sa)%spectag == Spectag)
!
          if (spectroscopic(sa)%der_calc(band)) then
!
            CA = spectroscopic(sa)%type
!
            do ip = 1, spectroscopic(sa)%no_phi_values
!
              do iz = 1, spectroscopic(sa)%no_zeta_values
!
! Run through representation basis coefficients
!
                sv_elmnt = m + iz - 1
!
                do ptg_i = 1, kconv
                  if (CA == 'W') then
                    VRad(ptg_i) = k_star_spect_dw(ptg_i+shft,iz,ip,nf)
                  else if (CA == 'N') then
                    VRad(ptg_i) = k_star_spect_dn(ptg_i+shft,iz,ip,nf)
                  else if (CA == 'V') then
                    VRad(ptg_i) = k_star_spect_dnu(ptg_i+shft,iz,ip,nf)
                  end if
                end do
!
                if (do_conv) then
!
! Now Convolve the derivative
!
                  fft_angles(1:kconv) = ptg_angle(1:kconv+shft)
!
                  Call fov_convolve(fft_angles,VRad,cen_angle,con_type, &
     &                 kconv,band,fft_pts,InDir,Aaap,Ier)
                  if (Ier /= 0) Return
!
                  if (fft_index(1) > 0) then
                    do ptg_i = 1, ktr
                      ilo = fft_index(ptg_i)
                      VRad(ptg_i) = VRad(ilo)
                    end do
                  end if
!
                end if
!
! Interpolate onto the output grid, and store in k_star_all ..
!
                Call Lintrp(fft_press,PtP,VRad,TmpRad,ktr,j)
                k_star_all(sv_elmnt,ip,1:j) = TmpRad(1:j)
!
              end do
!
            end do
!
          end if
!
! *** DEBUG
!
          TmpRad(1:Nptg) = 0.0
          Ax = 'di_d'//spectroscopic(sa)%name
          Call StrLwr(Ax)
          s_np = spectroscopic(sa)%no_phi_values
          s_nz = spectroscopic(sa)%no_zeta_values
          r = -1.0
          TmpRad(1:s_nz) = spectroscopic(sa)%zeta_basis(1:s_nz)
          Call Hunt(r,TmpRad,s_nz,iz,ihi)
          if (abs(spectroscopic(sa)%zeta_basis(iz)-r) >                &
     &        abs(spectroscopic(sa)%zeta_basis(ihi)-r)) iz = ihi
          r = spectroscopic(sa)%zeta_basis(iz)
          sv_elmnt = m + iz - 1
          ilo = len_trim(Ax)
          write(*,909) pfa_ch(jch),Ax(5:ilo),iz,r
          do ip = 1, s_np
            write(*,918) Ax(1:ilo),iz,ip,char(92),j
            write(*,908) (k_star_all(sv_elmnt,ip,ptg_i),ptg_i=1,j)
          end do
!
908   Format(6(1x,1pe12.5))
918   Format(a,i2.2,'_phi_',i1,a1,i4.4,'n')
!
909   Format(/,5x,'CONVOLVED k_star_spect_dump for channel:',i3,/,4x,  &
     &'Derivatrives of Radiance with respect to ',a,i2.2,              &
     &' (Zeta=',f7.4,')',/)
!
! *** END DEBUG
!
          sa = sa + 1
!
        end do
      end if
    end do
!
941 format('** Error in subroutine: pfa_fwd_mdl ..',a,/,               &
   &       '   Exceeds maximum number of derivatives:',i5,/,           &
   &       '   New value should be at least:',i4,' large',/,           &
   &       '   Need to change it in: l2pc_file_paramters.inc',/,       &
   &       '   parameter (max_no_sv_derivatives =',i4,')')
!
      Return
  End Subroutine PFA_FWD_MDL
end module PFA_FWD_MDL_M
! $Log$
! Revision 1.2  2000/07/06 00:11:44  zvi
!  This is the Freeze version of Jun/24/2000
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
