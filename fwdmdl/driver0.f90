program DRIVER0

! This driver doesn't actually do anything.  It only serves to
! access the Forward Model, FWD_MDL, to make sure it can link.

  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES
  use FWD_MDL_M, only: FWD_MDL
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO,     &
                                  L2PC_HEADER_TRI
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS,  &
                                 PFA_SLAB, SPECTRO_PARAM
  use L2PCDim, only: MNP => max_no_phi, NLVL, NSPS
  use MDBETA, only: MNF => max_no_freq, NO_T_PHI
  use MLSCommon, only: I4, R4, R8

  integer, parameter :: NPATH = 1

    type(l2pc_header_one) :: header1
    type(l2pc_header_two) :: header2
    type(l2pc_header_tri) :: header3
    integer(i4) :: N_LVLS
    integer(i4) :: FFT_PTS
    integer(i4) :: NO_CONV_HTS
    integer(i4) :: N_TAN(1)
    integer(i4) :: BAND
    integer(i4) :: NDX_SPS(Nsps,1)
    integer(i4) :: NO_FILT_PTS
    integer(i4) :: NO_GEOM
    integer(i4) :: NO_GEOPHYS
    integer(i4) :: NO_INT_FRQS(1)
    integer(i4) :: NO_PFA_CH
    integer(i4) :: NO_SPS_TBL(1)
    integer(i4) :: SPS_TBL(Nsps,1)
    integer(i4) :: PFA_CH(1)
    integer(i4) :: NVR(1)
    Real(r4) :: Z_PATH(Npath,1), H_PATH(Npath,1)
    Real(r4) :: T_PATH(Npath,1), PHI_PATH(Npath,1)
    Real(r4) :: DHDZ_PATH(Npath,1)
    Real(r4) :: MDB_PRES(1)
    Real(r4) :: MDB_TEMP(1)
    Real(r8) :: MDB_FREQ(mnf,1)
    Real(r4) :: Z_GRID(1), T_GRID(1), H_GRID(1)
    Real(r4) :: DH_DT_GRID(Nlvl,1), A_GRID(1), N_GRID(1)
    Real(r4) :: V_GRID(1)
    integer(i4) :: T_INDX
    Real(r4) :: CONV_HTS(1)
    Real(r4) :: CONV_PRESS(1)
    Real(r4) :: CONV_TEMP(1)
    Real(r4) :: SPSFUNC(Nlvl,1)
    Real(r4) :: CEN_ANGLE
    Real(r4) :: EARTH_REF
    Real(r4) :: T_TAN(1)
    Real(r4) :: ACC(1)
    Real(r4) :: THRESHOLD(1)
    type(atmos_comp) :: ATMOSPHERIC(1)
    Real(r4) :: CS(Nlvl,no_t_phi,mnf,1)
    Real(r4) :: D2X_DXDT(Nlvl,1)
    Real(r4) :: DX_DT(Nlvl,1)
    Real(r8) :: F_GRID_FILTER(maxfiltpts,1)
    Real(r8) :: F_GRID(maxaitkenpts,1)
    Real(r8) :: FILTER_FUNC(maxfiltpts,1)
    Real(r8) :: FREQ(1)
    integer(i4) :: NO_FREQS_F(1)
    Real(r4) :: G_BASIS(mxco,1)
    type(geom_param) :: GEOMETRIC(1)
    type(geophys_param) :: GEOPHYSIC(1)
    Real(r4) :: F_BASIS(mxco,1)
    integer(i4) :: NO_COEFFS_F(1)
    Real(r4) :: MR_F(mxco,mnp,1)
    integer(i4) :: NO_COEFFS_G(1)
    Real(r4) :: MR_G(mxco,mnp,1)
    type(pfa_slab) :: PFA_SPECTRUM(6,1)
    Real(r4) :: PTG_ANGLE(1)
    Real(r4) :: PTG_HTS(1)
    type(limb_press) :: PTG_PRESS
    Real(r4) :: REF_CORR(Nlvl,1)
    Real(r4) :: AZIM_183, AZIM_205, AZIM_REF
    Real(r4) :: C_PITCH, C_ROLL, C_YAW
    Real(r4) :: ELEV_183, ELEV_205
    Real(r4) :: GEOCSRAD
    Real(r4) :: S_PITCH, S_ROLL
    Real(r4) :: SP_TMP
    Real(r4) :: S_YAW
    integer(i4) :: SI
    character(len=1) :: AAAP
    character(len=1) :: INDIR
    character(len=1) :: KEYS(1)
    character(len=1) :: RUNF
    integer(i4) :: TIME_I
    integer(i4) :: L2PC_LU
    integer(i4) :: RECN
    integer(i4) :: MXBIN_RECN
    integer(i4) :: LRUN
    integer(i4) :: REC_NOS(1)
    integer(i4) :: NO_ATMOS
    integer(i4) :: JKEY
    integer(i4) :: N_OBS
    integer(i4) :: ATMOS_INDEX(1)
    integer(i4) :: GEOM_INDEX(1)
    integer(i4) :: GEOPHYS_INDEX(1)
    Logical*1 :: DO_CONV
    integer(i4) :: CH1, CH2
    integer(i4) :: SPECT_ATMOS(1)
    integer(i4) :: NO_PHI_G(1)
    integer(i4) :: NO_PHI_F(1)
    Real(r4) :: PHI_BASIS_F(mnp,1)
    integer(i4) :: NO_PHI_VEC(1)
!   integer(i4) :: NPATH
    integer(i4) :: PATH_BRKPT(3,1)
    integer(i4) :: NO_PHI_T
    Real(r4) :: T_PHI_BASIS(1)
!   Real(r4) :: DH_DT_PATH(Npath,mnp,mxco,1)
    Real(r4) :: DH_DT_PATH(Npath,mnp,23,1)
    type(eos_mdb_hdr) :: MDB_HDR(1)
    type(eos_mdb_rec) :: MDB_REC(max_no_lines,1)
    type(spectro_param) :: SPECTROSCOPIC(1)
    integer(i4) :: Ier
    call fwd_mdl ( header1, header2, header3, N_lvls, fft_pts,         &
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
end program DRIVER0
