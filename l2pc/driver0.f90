program DRIVER0
! This driver doesn't actually do anything.  It only serves to
! access the Forward Model, FWD_MDL, to make sure it can link.
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES, MAX_TEMP, &
                     MAX_ZETA, MAX_FREQ, MNF => max_freq
  use FWD_MDL_M, only: FWD_MDL
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component, &
        MAX_NO_POINTINGS, MAX_NO_KEY_ADDR, MAX_NO_SV_ELMNTS, MAX_TABLE_2D
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO,     &
                                  L2PC_HEADER_TRI
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS,  &
                                 PFA_SLAB, SPECTRO_PARAM, MAXPFACH, MAXGEOM, &
                                 MAXGEOPHYS, MAXLINES, MAXRAT
  use L2PCDim, only: MNP => max_no_phi, NLVL, NSPS, NPTG, NCH
  use MDBETA, only: CS_MNF => max_no_freq, NO_T_PHI
  use ELLIPSE
  use MLSCommon, only: I4, R4, R8
!
  IMPLICIT NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
  integer, parameter :: NPATH = 1400
!
    Character(len=32) :: KEYS(max_no_key_addr)
    Character(len=80) :: AAAP
    Character(len=80) :: INDIR
    Character(len=80) :: RUNF
!
    Integer(i4) :: ATMOS_INDEX(Nsps)
    Integer(i4) :: BAND
    Integer(i4) :: CH1, CH2
    Integer(i4) :: FFT_PTS
    Integer(i4) :: GEOM_INDEX(maxgeom)
    Integer(i4) :: GEOPHYS_INDEX(maxgeophys)
    Integer(i4) :: Ier
    Integer(i4) :: JKEY
    Integer(i4) :: L2PC_LU
    Integer(i4) :: LRUN
    Integer(i4) :: MXBIN_NOR
    Integer(i4) :: NDX_SPS(Nsps,maxpfach)
    Integer(i4) :: NO_ATMOS
    Integer(i4) :: NO_COEFFS_F(Nsps)
    Integer(i4) :: NO_COEFFS_G(Nsps)
    Integer(i4) :: NO_CONV_HTS
    Integer(i4) :: NO_FILT_PTS
    Integer(i4) :: NO_FREQS_F(Nsps)
    Integer(i4) :: NO_GEOM, NO_SPECTRO
    Integer(i4) :: NO_GEOPHYS
    Integer(i4) :: NO_INT_FRQS(maxpfach)
    Integer(i4) :: NO_PFA_CH
    Integer(i4) :: NO_PHI_F(Nsps)
    Integer(i4) :: NO_PHI_G(maxgeophys)
    Integer(i4) :: NO_PHI_T
    Integer(i4) :: NO_PHI_VEC(max_no_sv_elmnts)
    Integer(i4) :: NO_SPS_TBL(6)
    Integer(i4) :: NVR(50)
    Integer(i4) :: N_LVLS
    Integer(i4) :: N_OBS
    Integer(i4) :: N_TAN(Nlvl)
    Integer(i4) :: PATH_BRKPT(3,Nptg)
    Integer(i4) :: PFA_CH(maxpfach)
    Integer(i4) :: NOR
    Integer(i4) :: REC_NOS(max_no_key_addr)
    Integer(i4) :: SI
    Integer(i4) :: SPECT_ATMOS(Nsps)
    Integer(i4) :: SPS_TBL(Nsps,6)
    Integer(i4) :: TIME_I
    Integer(i4) :: T_INDEX
!
    Logical*1 :: DO_CONV
!
    Real(r8) :: ACC(maxpfach)
    Real(r8) :: AZIM_183, AZIM_205, AZIM_REF
    Real(r8) :: CEN_ANGLE
    Real(r8) :: CONV_HTS(Nptg)
    Real(r8) :: CONV_PRESS(Nptg)
    Real(r8) :: CONV_TEMP(Nptg)
    Real(r8) :: CS(Nlvl,no_t_phi,cs_mnf,Nsps)
    Real(r8) :: C_PITCH, C_ROLL, C_YAW
    Real(r8) :: DH_DT_GRID(Nlvl,mxco), A_GRID(Nlvl), N_GRID(Nlvl)
!   Real(r4) :: DH_DT_PATH(Npath,mnp,mxco,Nptg)
    Real(r4) :: DH_DT_PATH(Npath,mnp,23,Nptg)
    Real(r4) :: DHDZ_PATH(Npath,Nptg)
    Real(r8) :: DX_DT(Nlvl,mxco)
    Real(r8) :: D2X_DXDT(Nlvl,mxco)
    Real(r8) :: EARTH_REF
    Real(r8) :: ELEV_183, ELEV_205
    Real(r8) :: F_BASIS(mxco,Nsps)
    Real(r8) :: GEOCSRAD
    Real(r8) :: G_BASIS(mxco,Nsps)
    Real(r8) :: MDB_PRES(Nlvl)
    Real(r8) :: MDB_TEMP(Nlvl)
    Real(r8) :: MR_F(mxco,mnp,Nsps)
    Real(r8) :: MR_G(mxco,mnp,Nsps)
    Real(r8) :: PHI_BASIS_F(mnp,Nsps)
    Real(r8) :: PTG_ANGLE(Nlvl)
    Real(r8) :: PTG_HTS(max_no_pointings)
    Real(r8) :: REF_CORR(Nlvl,Nptg)
    Real(r8) :: SPSFUNC(Nlvl,Nsps)
    Real(r8) :: S_TEMP
    Real(r8) :: S_PITCH, S_ROLL
    Real(r8) :: S_YAW
    Real(r8) :: THRESHOLD(maxpfach)
    Real(r8) :: T_PATH(Npath,Nptg), PHI_PATH(Npath,Nptg)
    Real(r8) :: T_PHI_BASIS(mnp)
    Real(r8) :: T_TAN(Nlvl)
    Real(r8) :: V_GRID(Nlvl)
    Real(r8) :: Z_GRID(Nlvl), T_GRID(Nlvl), H_GRID(Nlvl)
    Real(r8) :: Z_PATH(Npath,Nptg), H_PATH(Npath,Nptg)
!
    Real(r8) :: C_FREQ(Nch)
    Real(r8) :: FILTER_FUNC(maxfiltpts,maxpfach)
    Real(r8) :: F_GRID(maxaitkenpts,maxpfach)
    Real(r8) :: F_GRID_FILTER(maxfiltpts,maxpfach)
    Real(r8) :: MDB_FREQ(mnf,Nsps)
!
    Type(l2pc_header_one) :: header1
    Type(l2pc_header_tri) :: header3
    Type(l2pc_header_two) :: header2
!
    Type(atmos_comp) :: ATMOSPHERIC(2)
    Type(eos_mdb_hdr) :: MDB_HDR(2)
    Type(eos_mdb_rec) :: MDB_REC(max_no_lines,2)
    Type(geom_param) :: GEOMETRIC(1)
    Type(geophys_param) :: GEOPHYSIC(1)
    Type(limb_press) :: PTG_PRESS
    Type(pfa_slab) :: PFA_SPECTRUM(6,2)
    Type(spectro_param) :: SPECTROSCOPIC(3)
!
    Integer(i4) :: i,j,k,m,io,Kpath
    Character(len=76) :: Line
!
    Namelist/In1/N_lvls,fft_pts,no_conv_hts,n_tan,band,ndx_sps,  &
 &      no_filt_pts,no_geom,no_geophys,no_int_frqs, no_pfa_ch,   &
 &      no_sps_tbl,sps_tbl,pfa_ch,nvr,z_path,h_path, t_path,     &
 &      phi_path,dHdz_path,mdb_pres,mdb_temp, mdb_freq,z_grid,   &
 &      t_grid,h_grid,dh_dt_grid,a_grid,n_grid, v_grid,t_index,  &
 &      conv_hts,conv_press,conv_temp,spsfunc, earth_ref,t_tan,  &
 &      Acc,Threshold,cs, f_grid_filter,f_grid,filter_func,      &
 &      c_freq,no_freqs_f,g_basis,f_basis, no_coeffs_f,mr_f,     &
 &      no_coeffs_g,mr_g,ptg_angle,ptg_hts,ref_corr,             &
 &      azim_183,azim_205,azim_ref,c_pitch,c_roll,c_yaw,         &
 &      elev_183,elev_205,geocsrad,s_pitch,s_roll,s_temp,        &
 &      s_yaw,si
!
    Namelist/In2/Aaap,InDir,runf,time_i,l2pc_lu,nor,mxbin_nor,      &
 &      lrun,no_spectro,no_atmos,jkey,n_obs,atmos_index,geom_index, &
 &      geophys_index,do_conv,ch1,ch2,spect_atmos,no_phi_g,         &
 &      no_phi_f,phi_basis_f,no_phi_vec,Kpath,path_brkpt,           &
 &      no_phi_t,t_phi_basis
!
! ***  DEBUG
!
  Line(1:)=' '
  Line = '/home/zvi/fwdmdl_dump.asc'
  close(47,iostat=io)
  open(47,file=Line,status='OLD',iostat=io)
  if(io /= 0) goto 99
!
! read structure: header1
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  Print *,'Reading: ',Trim(Line)
  read(47,*,iostat=io) header1%no_bands,header1%no_channels_per_band, &
 &          header1%no_pointings,header1%no_avail_keys, &
 &          header1%no_sv_components,header1%no_coeff_per_component, &
 &          header1%no_k_records_per_bin,header1%no_mag_fields, &
 &          header1%no_b_theta_lin_val
!
  k = header1%no_sv_components
  m = header1%no_b_theta_lin_val
  read(47,'(A)',iostat=io) header1%line1
  if(io /= 0) goto 99
  read(47,'(A)',iostat=io) header1%line2
  if(io /= 0) goto 99
  read(47,'(A)',iostat=io) header1%line3
  if(io /= 0) goto 99
  read(47,'(A)',iostat=io) (header1%avail_keys(i),i=1,header1%no_avail_keys)
  if(io /= 0) goto 99
  read(47,'(A)',iostat=io) (header1%sv_components(i),i=1,k)
  if(io /= 0) goto 99
  read(47,*,iostat=io) ((header1%sv_rtrvl_by_band(i,j),j=1,6),i=1,k)
  if(io /= 0) goto 99
  read(47,*,iostat=io) (header1%no_elmnts_per_sv_component(i),i=1,k)
  if(io /= 0) goto 99
  read(47,*,iostat=io) (header1%sv_component_first_elmnt_index(i),i=1,k)
  if(io /= 0) goto 99
  read(47,*,iostat=io) (header1%pointings(i),i=1,header1%no_pointings)
  if(io /= 0) goto 99
  read(47,*,iostat=io) (header1%b_fields(i),i=1,header1%no_mag_fields)
  if(io /= 0) goto 99
  read(47,*,iostat=io) header1%b_phi_lin_val
  if(io /= 0) goto 99
  if(m > 0) read(47,*,iostat=io) (header1%b_theta_lin_val(i),i=1,m)
  if(io /= 0) goto 99
!
! read structure: header2
!
  read(47,*,iostat=io) m
  if(io /= 0) goto 99
  header2%no_sv_elmnts = m
  if(m > 0) read(47,*,iostat=io) (header2%tri_basis_vert_grid(i),i=1,m)
  if(io /= 0) goto 99
!
! read structure: header3
!
! read(47,*,iostat=io) header3%matdim
! if(io /= 0) goto 99
! read(47,*,iostat=io) (header3%second_der_matrix_bands(i),i=1,max_table_2d)
! if(io /= 0) goto 99
! read(47,*,iostat=io) (header3%second_der_matrix_namid(i),i=1,max_table_2d)
! if(io /= 0) goto 99
!
! read the ELLIPSE Common block data
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  Print *,'Reading: ',Trim(Line)
  read(47,*,iostat=io) a2,c2,c2oa2,cpt,spt,cps,sps,cpts,spts,ht,ht2,   &
 &           Rr,Phi_tan,NPhi_tan,Phi_s,NPhi_s,ps,RoC,XoC,YoC ! ,EarthX
  if(io /= 0) goto 99
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,NML=In1,iostat=io)
  if(io /= 0) goto 99
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,NML=In2,iostat=io)
  if(io /= 0) goto 99
!
! read structure: atmospheric(m),m=1,2
!
  k =  mxco
  do m = 1, 2
    Line(1:)=' '
    read(47,'(A)',iostat=io) Line
    if(io /= 0) goto 99
    Print *,'Reading: ',Trim(Line)
    read(47,'(A)',iostat=io) atmospheric(m)%name
    if(io /= 0) goto 99
    read(47,*,iostat=io) atmospheric(m)%spectag
    if(io /= 0) goto 99
    read(47,*,iostat=io) atmospheric(m)%no_lin_values
    if(io /= 0) goto 99
    read(47,*,iostat=io)(atmospheric(m)%fwd_calc(i),i=1,6)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(atmospheric(m)%der_calc(i),i=1,6)
    if(io /= 0) goto 99
    atmospheric(m)%der_calc(5) = .true.         ! DEBUG
    read(47,*,iostat=io)(atmospheric(m)%lin_val(i),i=1,k)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(atmospheric(m)%basis_peaks(i),i=1,k+2)
    if(io /= 0) goto 99
  end do
!
! read structure: geometric(1)
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,'(A)',iostat=io) geometric(1)%name
  if(io /= 0) goto 99
  read(47,*,iostat=io)(geometric(1)%der_calc(i),i=1,6)
  if(io /= 0) goto 99
  geometric(1)%der_calc(5) = .true.         ! DEBUG
  read(47,*,iostat=io) geometric(1)%lin_val
  if(io /= 0) goto 99
!
! read structure: geophysic(1)
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,'(A)',iostat=io) geophysic(1)%name
  if(io /= 0) goto 99
  read(47,*,iostat=io) geophysic(1)%no_lin_values
  if(io /= 0) goto 99
  read(47,*,iostat=io)(geophysic(1)%der_calc(i),i=1,6)
  if(io /= 0) goto 99
  geophysic(1)%der_calc(5) = .true.         ! DEBUG
  read(47,*,iostat=io)(geophysic(1)%lin_val(i),i=1,k)
  if(io /= 0) goto 99
  read(47,*,iostat=io)(geophysic(1)%basis_peaks(i),i=1,k)
  if(io /= 0) goto 99
!
! read structure: ptg_press
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,'(A)',iostat=io) ptg_press%name
  if(io /= 0) goto 99
  read(47,*,iostat=io) ptg_press%no_lin_values
  if(io /= 0) goto 99
  read(47,*,iostat=io)(ptg_press%der_calc(i),i=1,6)
  if(io /= 0) goto 99
  read(47,*,iostat=io)(ptg_press%lin_val(i),i=1,max_no_pointings)
  if(io /= 0) goto 99
!
! read structure: pfa_spectrum(5,1) & pfa_spectrum(5,2)
!
  do m = 1, 2
    Line(1:)=' '
    read(47,'(A)',iostat=io) Line
    if(io /= 0) goto 99
    Print *,'Reading: ',Trim(Line)
    read(47,'(A)',iostat=io) pfa_spectrum(5,m)%sps_name
    if(io /= 0) goto 99
    read(47,*,iostat=io) pfa_spectrum(5,m)%no_sps
    if(io /= 0) goto 99
    read(47,*,iostat=io) pfa_spectrum(5,m)%sps_spectag
    if(io /= 0) goto 99
    read(47,*,iostat=io) pfa_spectrum(5,m)%no_lines
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%Nrat(i),i=1,Nlvl)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%varM(i),i=1,Nlvl)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_v0(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_el(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_str(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_w(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_ps(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_n(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_n1(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_n2(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_gamma(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_delta(i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_part(1,i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_part(2,i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(pfa_spectrum(5,m)%sps_part(3,i),i=1,maxlines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)((pfa_spectrum(5,m)%Xx(i,j),j=1,Nlvl),i=1,maxrat)
    if(io /= 0) goto 99
    read(47,*,iostat=io)((pfa_spectrum(5,m)%Yy(i,j),j=1,Nlvl),i=1,maxrat)
    if(io /= 0) goto 99
    read(47,*,iostat=io)((pfa_spectrum(5,m)%Dy(i,j),j=1,Nlvl),i=1,maxrat)
    if(io /= 0) goto 99
  end do
!
! read structure: spectroscopic(1)
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  do m = 1, no_spectro
    read(47,'(A)',iostat=io) spectroscopic(m)%type
    if(io /= 0) goto 99
    read(47,'(A)',iostat=io) spectroscopic(m)%name
    if(io /= 0) goto 99
    read(47,*,iostat=io) spectroscopic(m)%spectag
    if(io /= 0) goto 99
    read(47,*,iostat=io) spectroscopic(m)%no_phi_values
    if(io /= 0) goto 99
    read(47,*,iostat=io) spectroscopic(m)%no_zeta_values
    if(io /= 0) goto 99
    read(47,*,iostat=io)(spectroscopic(m)%der_calc(i),i=1,6)
    if(io /= 0) goto 99
    spectroscopic(m)%der_calc(5) = .true.       ! DEBUG
    read(47,*,iostat=io)(spectroscopic(m)%phi_basis(i),i=1,mnp+2)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(spectroscopic(m)%zeta_basis(i),i=1,k+2)
    if(io /= 0) goto 99
  end do
!
! read structure: mdb_hdr(1) & mdb_hdr(2)
!
  do m = 1, 2
    Line(1:)=' '
    read(47,'(A)',iostat=io) Line
    if(io /= 0) goto 99
    Print *,'Reading: ',Trim(Line)
    read(47,*,iostat=io) mdb_hdr(m)%Spectag
    if(io /= 0) goto 99
    read(47,*,iostat=io) mdb_hdr(m)%no_lines
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%no_f_grid(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%el(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%log_i(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%n(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%w(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%delta(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%n1(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%n2(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%gamma(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%v0(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%ps(i),i=1,max_no_lines)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%q_log(i),i=1,3)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%Zeta(i),i=1,max_zeta)
    if(io /= 0) goto 99
    read(47,*,iostat=io)(mdb_hdr(m)%Log_Temp(i),i=1,max_temp)
    if(io /= 0) goto 99
    read(47,*,iostat=io)((mdb_hdr(m)%x_grid(i,j),j=1,max_no_lines), &
 &                                i=1,max_freq)
    if(io /= 0) goto 99
  end do
!
! read structure: mdb_rec(1,1) & mdb_rec(1,2)
!
  do m = 1, 2
    Line(1:)=' '
    read(47,'(A)',iostat=io) Line
    if(io /= 0) goto 99
    Print *,'Reading: ',Trim(Line)
    read(47,*,iostat=io) (((mdb_rec(1,m)%Log_beta(i,j,k),                  &
 &                           k=1,max_freq),j=1,max_temp),i=1,max_zeta)
    if(io /= 0) goto 99
    read(47,*,iostat=io) (((mdb_rec(1,m)%dLog_beta_dw(i,j,k),              &
 &                           k=1,max_freq),j=1,max_temp),i=1,max_zeta)
    if(io /= 0) goto 99
    read(47,*,iostat=io) (((mdb_rec(1,m)%dLog_beta_dn(i,j,k),              &
 &                           k=1,max_freq),j=1,max_temp),i=1,max_zeta)
    if(io /= 0) goto 99
    read(47,*,iostat=io) (((mdb_rec(1,m)%dLog_beta_dNu0(i,j,k),            &
 &                           k=1,max_freq),j=1,max_temp),i=1,max_zeta)
    if(io /= 0) goto 99
    read(47,'(76A1)',iostat=io) (((mdb_rec(1,m)%Log_beta_intrp(i,j,k),     &
 &                           k=1,max_freq),j=1,max_temp),i=1,max_zeta)
    if(io /= 0) goto 99
  end do
!
  Line(1:)=' '
  read(47,'(A)',iostat=io) Line
  if(io /= 0) goto 99
  Print *,'Reading: ',Trim(Line)
  read(47,*,iostat=io) dh_dt_path
  if(io /= 0) goto 99
!
  close(47,iostat=i)
!
! ***  END DEBUG
!
    call fwd_mdl ( header1, header2, header3, N_lvls, fft_pts,             &
 &           no_conv_hts, n_tan, band, ndx_sps, no_filt_pts, no_geom,      &
 &           no_geophys, no_int_frqs, no_pfa_ch, no_sps_tbl, sps_tbl,      &
 &           pfa_ch, nvr, z_path, h_path, t_path, phi_path, dHdz_path,     &
 &           mdb_pres, mdb_temp, mdb_freq, z_grid, t_grid, h_grid,         &
 &           dh_dt_grid, a_grid, n_grid, v_grid, t_index, conv_hts,        &
 &           conv_press, conv_temp, spsfunc, cen_angle, earth_ref,         &
 &           t_tan,  Acc, Threshold, atmospheric, cs, d2x_dxdt, dx_dt,     &
 &           f_grid_filter, f_grid, filter_func,c_freq, no_freqs_f,        &
 &           g_basis, geometric, geophysic, f_basis, no_coeffs_f,          &
 &           mr_f,  no_coeffs_g, mr_g, pfa_spectrum, ptg_angle,            &
 &           ptg_hts,  ptg_press, ref_corr, azim_183, azim_205,            &
 &           azim_ref, c_pitch,  c_roll, c_yaw, elev_183, elev_205,        &
 &           geocsrad, s_pitch, s_roll,  s_temp, s_yaw, si, Aaap,          &
 &           InDir, keys, runf, time_i, l2pc_lu,  nor, mxbin_nor,          &
 &           lrun, rec_nos, no_atmos, jkey, n_obs,  atmos_index,           &
 &           geom_index, geophys_index, do_conv, ch1, ch2,                 &
 &           spect_atmos, no_phi_g, no_phi_f, phi_basis_f, no_phi_vec,     &
 &           Npath, path_brkpt, no_phi_t, t_phi_basis, dh_dt_path,         &
             mdb_hdr, mdb_rec, spectroscopic, Ier )
!
 99   close(47,iostat=i)
!
      if(io /= 0) then
        Call ErrMsg('** Error in Driveri program',io)
      endif
!
end program DRIVER0
! Log: DRIVER0,v $
