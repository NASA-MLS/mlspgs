Program L2PC_FILE
  use UNITS, only: p_vs_h_unit
  use EOS_MDB, only: EOS_MDB_HDR,EOS_MDB_REC,MAX_NO_LINES,CS_MNF => MAX_FREQ
  use FWD_MDL_M, only: FWD_MDL
  use D_HUNT_M, only: HUNT
  use GL6P, only: NG
  use L2PC_FILE_PARAMETERS, only: MAX_NO_BANDS, MAX_NO_KEY_ADDR, &
                                  MAX_NO_POINTINGS, MAX_NO_SV_ELMNTS, &
                                  MXCB => max_no_channels_per_band, &
                                  MXCO => max_no_elmnts_per_sv_component
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO, &
                                  L2PC_HEADER_TRI
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS,  &
                                 PFA_SLAB, SPECTRO_PARAM,  MAXPFACH, &
                                 MAXGEOPHYS, MAXGEOM, MAXLINES, MAXRAT
  use L2PCdim, only: MNP => max_no_phi, NCH, NLVL, NPTG, NSPS
  use MDBETA, only: MAX_NO_ZETA, MNF => max_no_freq, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use READ_UIF_M, only: READ_UIF
  use L2PC_FILE_MNGMT_SW_M, only: OPEN_L2PC, COMPOSE_HEADER_LINES, &
           WRITE_HDR, POSITION_L2PC, CLOSE_L2PC_XX, RE_WRITE_HEADER1
  use TIME_MOD, only: DIFF_TEXT_TIME, DURATION_TEXT, NOW_CCSDS, TK, &
                      Z_DATETIME
  use STRINGS, only: LEFTJ, SQZSTR, STRUPR, STRLWR
  use VELCOR_M, only: VELCOR
  use GET_ANGLES_M, only: GET_ANGLES
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use FWD_MDL_SET_UP_M, only: FWD_MDL_SET_UP
!
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!  "Crash Proof" L2PC Program
!
  integer(i4), parameter :: NPATH = 2 * Nlvl * (Ng+1)
!
  real(r8) :: A_GRID(Nlvl)
  character(len=80) :: AAAP
  real(r8) :: ACC(maxpfach)                         ! PFA
  character(len=40) :: ACON, AKEYS(17)
  integer(i4) :: ATMOS_INDEX(Nsps)
  type (atmos_comp) :: ATMOSPHERIC(Nsps)
  real(r8) :: AZIM_REF, AZIM_183, AZIM_205          ! PFA
  integer(i4) :: BAND, BAND_1, BAND_2
  character(len=80) :: BUFF
  real(r8) :: C_FREQ(Nch)
  real(r8) :: C_PITCH, C_ROLL, C_YAW                ! PFA
  real(r8) :: CEN_ANGLE
  integer(i4) :: CH1, CH2
  real(r8) :: CONV_HTS(Nptg), CONV_HTS_RAW(Nptg), CONV_HTS_SAVE(Nptg)
  integer(i4) :: CONV_INDX(Nptg)
  real(r8) :: CONV_PRESS(Nptg), CONV_TEMP(Nptg)
  real(r8) :: CS(Nlvl,no_t_phi,cs_mnf,Nsps)
  character(len=80) :: DFILENAME
  integer(i4) :: DCH1, DCH2
  real(r8) :: DH_DT_GRID(Nlvl,mxco)
! Real(r4) :: DH_DT_PATH(Npath,mnp,mxco,Nptg)
  Real(r4) :: DH_DT_PATH(Npath,mnp,23,Nptg)
  real(r4) :: DHDZ_PATH(Npath,Nptg)
  integer(i4) :: DMXBIN_NOR
  Logical*1 :: DO_CONV
  character(len=80) :: DUMPF
  real(r8) :: DX_DT(Nlvl,mxco)
  real(r8) :: D2X_DXDT(Nlvl,mxco)
  real(r8) :: EARTH_RADIUS, EARTH_REF
  real(r8) :: ELEV_183, ELEV_205                    ! PFA
  real(r8) :: F_BASIS(mxco,Nsps)
  real(r8) :: F_GRID(maxaitkenpts,maxpfach)         ! PFA
  real(r8) :: F_GRID_FILTER(maxfiltpts,maxpfach)    ! PFA
  integer(i4) :: FFT_PTS
  character(len=80) :: FILENAME
  real(r8) :: FILTER_FUNC(maxfiltpts,maxpfach)      ! PFA
  character(len=80) :: FND
  real(r8) :: G_BASIS(mxco,Nsps)
  real(r8) :: GEOCSRAD                              ! PFA
  integer(i4) :: GEOM_INDEX(maxgeom)
  type (geom_param) :: GEOMETRIC(maxgeom)
  integer(i4) :: GEOPHYS_INDEX(maxgeophys)
  type (geophys_param) :: GEOPHYSIC(maxgeophys)
  real(r8) :: H, H_GNLV(max_no_zeta), H_GRID(Nlvl), H_OBS
  real(r8) :: H_PATH(Npath,Nptg)
  type (l2pc_header_one) :: HEADER1
  type (l2pc_header_two) :: HEADER2
  type (l2pc_header_tri) :: HEADER3
  real(r8) :: HREF
  integer(i4) :: I, IDT, IER
  character(len=80) :: INDIR
  integer(i4) :: IOS
  integer(i4) :: J, JKEY
  integer(i4) :: K
  Character(len=32) :: KEYS(max_no_key_addr)
  integer(i4) :: LD, LDF, LF
  type (limb_press) :: PTG_PRESS
  character(len=80) :: LINES(3)
  integer(i4) :: LRUN
  integer(i4) :: L2PC_LU, L2PC_LU_KEY, L2PC_REC_LENGTH
  integer(i4) :: M
  real(r8) :: MDB_FREQ(cs_mnf,Nsps)
  type (eos_mdb_hdr) :: MDB_HDR(02)                  ! 2 should be: maxsps
  real(r8) :: MDB_PRES(Nlvl)
  type (eos_mdb_rec) :: MDB_REC(max_no_lines,02)     ! 2 should be: maxsps
  real(r8) ::MDB_TEMP(Nlvl)
  real(r8) :: MR_F(mxco,mnp,Nsps), MR_G(mxco,mnp,Nsps)
  integer(i4) :: MXBIN_NOR
  real(r8) :: N_GRID(Nlvl)
  integer(i4) :: N_LVLS, N_OBS
  integer(i4) :: N_TAN(Nlvl)                        ! PFA
  integer(i4) :: NCUP
  integer(i4) :: NDX_SPS(Nsps,maxpfach)             ! PFA
  integer(i4) :: NEXT_BIN, NKEY
  integer(i4) :: NO_ATMOS
  Logical*1 :: NO_CHECK
  integer(i4) :: NO_COEFFS_F(Nsps), NO_COEFFS_G(Nsps)
  integer(i4) :: NO_CONV_HTS, NO_CONV_HTS_SAVE, NO_CRASH
  integer(i4) :: NO_FILT_PTS                        ! PFA
  integer(i4) :: NO_FREQS_F(Nsps)
  integer(i4) :: NO_GEOM, NO_GEOPHYS
  integer(i4) :: NO_INT_FRQS(maxpfach), NO_PFA_CH   ! PFA
  integer(i4) :: NO_PHI_F(Nsps)
  integer(i4) :: NO_PHI_G(maxgeophys),NO_PHI_VEC(max_no_sv_elmnts)
  integer(i4) :: NO_PHI_T, NO_SPECTRO
  integer(i4) :: NO_SPS_TBL(max_no_bands), NOR
  integer(i4) :: NVR(50)                            ! PFA
  integer(i4) :: P_INDX(Nlvl), PATH_BRKPT(3,Nptg)
  integer(i4) :: PFA, PFA_CH(maxpfach)
  real(r8) :: PHI_BASIS_F(mnp,Nsps), PHI_PATH(Npath,Nptg)
  type (pfa_slab) :: PFA_SPECTRUM(6,Nsps)
  character :: PRIMAG
  real(r8) :: PTG_ANGLE(Nlvl), PTG_HTS(max_no_pointings)
  integer(i4) :: REC_NOS(max_no_key_addr)
  real(r8) :: REF_CORR(Nlvl,Nptg)
  real(r8) :: ZROC
  character(len=80) :: RUNF
  real(r8) :: S_PITCH, S_ROLL                       ! PFA
  real(r8) :: S_TEMP
  real(r8) :: S_YAW                                 ! PFA
  integer(i4) :: SI, SPECT_ATMOS(Nsps), SPECT_INDEX(Nsps)
  type (spectro_param) :: SPECTROSCOPIC(3*Nsps)
  real(r8) :: SPSFUNC(Nlvl,Nsps)
  integer(i4) :: SPS_TBL(Nsps,max_no_bands), SV_I
  real(r8) :: T_GNLV(max_no_zeta), T_GRID(Nlvl)
  integer(i4) :: T_INDEX
  real(r8) :: T_PATH(Npath,Nptg), T_PHI_BASIS(mnp)
  real(r8) :: T_TAN(Nlvl), THRESHOLD(maxpfach)      ! PFA
  integer(i4) :: TIME_I
  character(len=8) :: TIME_STAMP, TS, TT
  character(len=80) :: UIF_NAME
  real(r8) :: V_GRID(Nlvl), VEL_Z, ZREF
  real(r8) :: Z, Z_GNLV(max_no_zeta), Z_GRID(Nlvl), Z_PATH(Npath,Nptg)
!
! For timing:
!
! intrinsic :: CPU_TIME           !  (Fortran 95 only !)
! intrinsic :: DATE_AND_TIME      !  (Fortran 95 only !)
  intrinsic DATE_AND_TIME
  real(tk) :: CPU, CPU_END, CPU_START, ELAPSED, TOTCPU, TOTELP
  character(len=21) :: ELAPSE_END, ELAPSE_START
!
!  Begin code:
!
  uif_name = 'tmp.dat'
! Call Gti(' Enter user input filename',uif_name)
  if(uif_name < '!') Stop
!
  do i = 1, max_no_key_addr
    rec_nos(i) = 0
    Keys(i)(1:) = ' '
  end do
!
!  Start actual program:
!
  Buff = uif_name
  Call Read_Uif ( Buff, header1, header2, header3, ptg_press,       &
 &     geometric, no_geom, geophysic, no_geophys, atmospheric,      &
 &     no_atmos, spectroscopic, no_spectro, spect_atmos,            &
 &     conv_hts_save, no_conv_hts, fft_pts, ch1, ch2,               &
 &     l2pc_lu, pfa, InDir, Fnd, Aaap, no_pfa_ch, no_filt_pts,      &
 &     pfa_ch, no_int_frqs, Acc, Threshold, pfa_spectrum, Primag,   &
 &     P_indx, Conv_indx, si, N_lvls, c_yaw, s_yaw, c_roll, s_roll, &
 &     c_pitch, s_pitch, elev_183, elev_205, azim_183, azim_205,    &
 &     azim_ref, geocsrad, atmos_index, geom_index, geophys_index,  &
 &     do_conv, spect_index, l2pc_rec_length, Ier )
  if(ier /= 0) goto 888
!
  filename = header1%line2
  ld = len_trim(InDir)
  lf = len_trim(Fnd)
!
  i = Index(Aaap,'/')
  do while(i > 0)
    Aaap(1:i) = ' '
    Call Leftj(Aaap)
    i = Index(Aaap,'/')
  end do
!
  Call StrLwr(Aaap)
  Call StrLwr(InDir)
!
  Buff(1:)=' '
  Buff = InDir(1:ld)//'eos_p_vs_h_md.dat'
  Call Get_PvsH(Buff,h_gNlv,t_gNlv,z_gNlv,i,Ier)
  if(ier /= 0) goto 888
!
! Get the selected integration grid pressures and make sure all basis_peaks
! for geophysical and atmospheric entities are falling on that grid:
!
  do i = 1, N_lvls
    j = p_indx(i)
    z_grid(i) = z_gNlv(j)
  end do
!
! Check geophysical entities basis_peaks and give warnings if they do
! not fall on the pre-selected integration grid:
!
  do i = 1, no_geophys
    j = 0
    m = -1
    band = -1
    k = geophysic(i)%no_lin_values
    do while(j < k .and. band < 0)
      j = j + 1
      z = geophysic(i)%basis_peaks(j)
      Call hunt(z,z_grid,N_lvls,m,idt)
      if(abs(z_grid(idt)-z) < abs(z_grid(m)-z)) m = idt
      if(abs(z-z_grid(m)) > 1.0e-3) then
        if(band < 1) then
          band = 2
          write (*,910) geophysic(i)%name
        end if
      end if
    end do
  end do
!
! Check atmospheric entities basis_peaks and give warnings if they do
! not fall on the pre-selected integration grid:
!
  do i = 1, no_atmos
    j = 0
    m = -1
    band = -1
    k = atmospheric(i)%no_lin_values
    do while(j < k .and. band < 0)
      j = j + 1
      z = atmospheric(i)%basis_peaks(j+1)
      Call hunt(z,z_grid,N_lvls,m,idt)
      if(abs(z_grid(idt)-z) < abs(z_grid(m)-z)) m = idt
      if(abs(z-z_grid(m)) > 1.0e-3) then
        if(band < 1) then
          band = 2
          write (*,910) atmospheric(i)%name
        end if
      end if
    end do
  end do
!
!  Set the Reference height & pressure to be the first entry in the EOS
!  "P vs. H"  file, so they are hydrostatically consistent..
!
  href = h_gNlv(1)
  zref = z_gNlv(1)
!
!  Set the "Surface" value of the conv_hts array to be: href. It does not
!  necessarily be zero, but it is consistent with zref.
!
  conv_hts_save(si) = href
!
!  Set conv_hts_save according to conv_indx array:
!
  m = si
  k = no_conv_hts
  do i = 1, k
    j = conv_indx(i)
    h = h_gNlv(j)
    if(h > conv_hts_save(m)) then
      m = m + 1
      conv_hts_save(m) = h
    end if
  end do
!
  no_conv_hts_save = m
!
  j = header1%no_channels_per_band
  band_1 = (ch1 + j - 1) / j
  band_2 = (ch2 + j - 1) / j
!
  ncup = 0
  no_crash = 0
  totelp = 0.0
  totcpu = 0.0
!
  nor = 0
  next_bin = 1
  mxbin_nor = 0
!
! Calculate the number of records in the file
!
  j = 0
  do band = band_1, band_2
    do sv_i = 1, header1%no_sv_components
      if(header1%sv_rtrvl_by_band(sv_i,band)) j = j + 1
    end do
  end do
!
  runf(1:)=' '
  runf = 'l2pc_run.dat'
  lrun = len_trim(runf)
!
  dumpf(1:) = ' '
  dumpf = 'l2pc_dump01.dat'
  ldf = len_trim(dumpf)
!
! Check if the run file exists ...
!
  idt = 0
  Acon(1:)=' '
  write (*,'(/,''** Checking run file ..'')')
!
  close (86,iostat=i)
  open (86,file=runf(1:lrun),status='OLD',iostat=ios)
!
  if(ios /= 0) then                 ! run file does not exists !
    no_crash = 0
    totelp = 0.0
    totcpu = 0.0
    no_check = .true.
!
! Write the l2pc file header
!
    dfilename = filename
    Call open_l2pc(filename,l2pc_lu,l2pc_lu_key,'NEW',              &
 &                 l2pc_rec_length,ier)
    if(ier /= 0) goto 888
    Call Compose_Header_Lines(filename,no_pfa_ch,pfa_ch,Lines,ier)
    if(ier /= 0) goto 888
    header1%line1 = Lines(1)
    header1%line2 = Lines(2)
    header1%line3 = Lines(3)
    Call write_hdr(header1,header2,header3,l2pc_lu,ier)
    if(ier /= 0) goto 888
    Print *,'** Writing dump file: ',dumpf(1:ldf)
    open (88,file=dumpf(1:ldf),form='UNFORMATTED',status='UNKNOWN',  &
 &       iostat=ier)
    if(ier /= 0) goto 888
    write (88) mxbin_nor,ch1,ch2,totelp,totcpu,no_crash,dfilename
    close (88,iostat=i)
!
  else                             ! run file exists - analize it
!
    j = 0
    idt = 1
    nkey = 0
    no_check = .false.
    do while(ios+j == 0)
      Buff(1:)=' '
      read(86,'(A)',iostat=ios) Buff
      Call Leftj(Buff)
      if(Buff(1:2) == 'D_') j = 1
    end do
    if(j > 0) then
      Acon(1:)=' '
      Acon(1:40)=Buff(1:40)
      do while(ios == 0)
        Buff(1:)=' '
        read(86,'(A)',iostat=ios) Buff
        Call Leftj(Buff)
        if(Buff(1:2) == 'D_') then
          if(Buff(1:10) /= Acon(1:10)) then
            nkey = nkey + 1
            Akeys(nkey) = Acon
          end if
          Acon(1:)=' '
          Acon(1:40)=Buff(1:40)
        end if
      end do
      nkey = nkey + 1
      Akeys(nkey) = Acon
      m = 0
      i = -1
      j = nkey + 1
      do while(j+m > 1)
        j = j -1
        Acon(1:)=' '
        Acon(1:40)=Akeys(j)
        read(Acon,'(11x,i2)',iostat=k) i
        if(k /= 0) i = -1
        if(i == header1%no_pointings) m = -200
      end do
      if(i == header1%no_pointings) nkey = j
      Acon(1:)=' '
      Acon(1:40)=Akeys(nkey)
      j = 0
      i = header1%no_avail_keys
      do while(j < i)
        j = j + 1
        if(header1%avail_keys(j)(1:7) == Acon(1:7)) then
          next_bin = j
          j = i + 4
        end if
      end do
    end if
!
    Acon(1:)=' '
    j = max(1,next_bin-1)
    Acon = header1%avail_keys(j)
    Print *,'** Last completed bin was:',Acon(1:7)
!
! Re-position the l2pc file to last record according to next_bin value:
!
    j = (next_bin - 1) * mxbin_nor + 4
    Call Position_l2pc(filename,l2pc_lu,l2pc_lu_key,j,                 &
 &                 l2pc_rec_length,ier)
    if(ier /= 0) goto 888
!
  end if
!
  close (86,iostat=i)
!
  if(.not.no_check) then
!
    close (88,iostat=i)
    Print *,'** Checking dump file: ',dumpf(1:ldf)
    open (88,file=dumpf(1:ldf),form='UNFORMATTED',status='OLD',         &
 &       iostat=ios)
!
    if(ios /= 0) then
!
      dfilename=filename
      close (88,iostat=i)
      Print *,'** Writing dump file: ',dumpf(1:ldf)
      open (88,file=dumpf(1:ldf),form='UNFORMATTED',status='UNKNOWN',   &
 &         iostat=ier)
      if(ier /= 0) goto 888
      write (88) mxbin_nor,ch1,ch2,totelp,totcpu,no_crash,dfilename
      close (88,iostat=i)
!
    else
!
      read(88,iostat=ios) dmxbin_nor,dch1,dch2,totelp,totcpu,           &
 &                       no_crash,dfilename
      close (88,iostat=i)
      if(ios /= 0) goto 777
      ios = 1
      if(dch1 == ch1 .and. dch2 == ch2) then
        ios = 2
        if(filename == dfilename) then
          ios = 3
          if(dmxbin_nor == mxbin_nor) ios = 0
        end if
      end if
!
      if(ios /= 0) then
        Print *,'** Mismatch in dump file: ',dumpf(1:ldf)
        if(ios == 1) then
          Print *,'      ch1=',ch1,'      ch2=',ch2
          Print *,' DUMP ch1=',dch1,'     ch2=',dch2
        else if(ios == 2) then
          Print *,'      filename=',filename
          Print *,' DUMP filename=',dfilename
        else if(ios == 3) then
          Print *,'      mxbin_nor=',mxbin_nor
          Print *,' DUMP mxbin_nor=',dmxbin_nor
        end if
        ier = 1
        goto 888
      end if
!
      if(ncup < 1) then
        ncup = 2
        no_crash = no_crash + 1
      end if
!
    end if
!
  end if
!
!  Open the run file to write:
!
  Acon(1:)=' '
  if(idt == 0) then                ! run file does not exist
    Acon = 'Started at:'
    open (86,file=runf(1:lrun),status='NEW',iostat=ios)
  else                             ! run file exists. Open for APPEND
    Acon = 'Re-Started at:'
    open (86,file=runf(1:lrun),position='APPEND',iostat=ios)
  end if
!
  Ts(1:)=' '
  Buff(1:) = ' '
  Lines(1) = Buff
  j = len_trim(Acon)
  Buff = '"Crash Proof" EOS F90/F95 Version L2PC. '//Acon(1:j)
  i = len_trim(Buff)
!
  Call Z_DATETIME(Acon)
  j = len_trim(Acon)
  Buff = Buff(1:i+1)//Acon(1:j)
  i = len_trim(Buff)
  do j = 1, i
    Lines(1)(j:j)='*'
  end do
!
  write (*,'(a,/,a,/,a)') Lines(1)(1:i),Buff(1:i),Lines(1)(1:i)
  write (86,'(a,/,a,/,a)') Lines(1)(1:i),Buff(1:i),Lines(1)(1:i)
!
! Now go through each available atmospheric bin
!
  jkey = 0
  do time_i = next_bin, header1%no_avail_keys
!
    jkey = jkey + 1
    if(nor >= 2) jkey = nor
!
    nor = 0
    time_stamp(1:)=' '
    write (time_stamp,'(a8)') header1%avail_keys(time_i)
    write (*,'('' Data stream for ...'',a8)') time_stamp
    write (86,'('' Data stream for ...'',a8)') time_stamp
!
    if(Ts /= time_stamp) then
!
      Ts = time_stamp
      Tt = time_stamp
      Call StrLwr(Tt)
      Call Velcor(InDir,Tt,vel_z,Ier)
      if(Ier /= 0) goto 888
!
      do i = 1, N_lvls
        v_grid(i) = vel_z
      end do
!
    end if
!
! Run the forward model set up for each time_stamp
!
    no_conv_hts = no_conv_hts_save
    do k = 1, no_conv_hts_save
      conv_hts_raw(k) = conv_hts_save(k)
    end do
!
    Call fwd_mdl_set_up(time_stamp,Primag,href,zref,ptg_hts,z_gNlv,    &
 &       p_indx,earth_radius,z_grid,t_grid,h_grid,dh_dt_grid,          &
 &       n_grid,v_grid,vel_z,mdb_pres,mdb_temp,mdb_freq,               &
 &       no_freqs_f,N_lvls,ptg_press,geometric,no_geom,                &
 &       geophysic,no_geophys,g_basis,mr_g,no_coeffs_g,                &
 &       atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,c_freq,cs,      &
 &       spsfunc,si,t_index,earth_ref,s_temp,h_obs,n_obs,              &
 &       conv_hts_raw,conv_hts,conv_press,conv_temp,no_conv_hts,       &
 &       no_sps_tbl,sps_tbl,no_pfa_ch,no_filt_pts,pfa_ch,              &
 &       no_int_frqs,pfa_spectrum,t_tan,ndx_sps,n_tan,f_grid_filter,   &
 &       f_grid,filter_func,ch1,ch2,InDir,ld,Fnd,lf,header1,           &
 &       no_phi_g,no_phi_f,phi_basis_f,atmos_index,geophys_index,      &
 &       no_phi_vec,z_path,h_path,t_path,phi_path,dHdz_path,           &
 &       dh_dt_path,Npath,path_brkpt,no_phi_t,t_phi_basis,ZRoC,        &
 &       mdb_hdr,mdb_rec,spectroscopic,no_spectro,spect_index,Ier)
    if(ier /= 0) goto 888
!
! Get 'N_lvls' pointing angles
! Establish the convolution angles with refraction correction
!
    Call get_angles(h_grid,conv_hts,n_grid,a_grid,ptg_angle,           &
 &                  ZRoC,h_obs,N_lvls,n_obs,no_conv_hts,Ier)
    if(Ier /= 0) goto 888
!
! Compute the refraction correction scaling matrix:
!
    Call refraction_correction(N_lvls,no_conv_hts,h_grid,n_grid,       &
 &                  conv_hts,earth_radius,ref_corr)
!
! Run the forward model, by bands
!
    do band = band_1, band_2
!
      nvr(1:no_pfa_ch) = -2
!
      k = no_conv_hts
!     call cpu_time ( cpu_start )           ! F95 only
      call now_ccsds ( elapse_start )
      Call fwd_mdl(header1,header2,header3,N_lvls,fft_pts,k,n_tan,     &
 &         band,ndx_sps,no_filt_pts,no_geom,no_geophys,no_int_frqs,    &
 &         no_pfa_ch,no_sps_tbl,sps_tbl,pfa_ch,nvr,z_path,h_path,      &
 &         t_path,phi_path,dHdz_path,mdb_pres,mdb_temp,                &
 &         mdb_freq,z_grid,t_grid,h_grid,dh_dt_grid,a_grid,n_grid,     &
 &         v_grid,t_index,conv_hts,conv_press,conv_temp,spsfunc,       &
 &         cen_angle,earth_ref,t_tan,Acc,Threshold,atmospheric,cs,     &
 &         d2x_dxdt,dx_dt,f_grid_filter,f_grid,filter_func,            &
 &         c_freq,no_freqs_f,g_basis,geometric,geophysic,f_basis,      &
 &         no_coeffs_f,mr_f,no_coeffs_g,mr_g,pfa_spectrum,ptg_angle,   &
 &         ptg_hts,ptg_press,ref_corr,azim_183,azim_205,azim_ref,      &
 &         c_pitch,c_roll,c_yaw,elev_183,elev_205,geocsrad,s_pitch,    &
 &         s_roll,s_temp,s_yaw,si,Aaap,InDir,Keys,runf,time_i,         &
 &         l2pc_lu,nor,mxbin_nor,lrun,rec_nos,no_atmos,jkey,n_obs,     &
 &         atmos_index,geom_index,geophys_index,do_conv,ch1,ch2,       &
 &         spect_atmos,no_phi_g,no_phi_f,phi_basis_f,no_phi_vec,       &
 &         Npath,path_brkpt,no_phi_t,t_phi_basis,dh_dt_path,           &
 &         mdb_hdr,mdb_rec,spectroscopic,Ier)
      if(ier /= 0) goto 888
!
!     call cpu_time ( cpu_end )           ! F95 only
      call now_ccsds ( elapse_end )
      cpu = cpu_end - cpu_start
      totcpu = totcpu + cpu
      elapsed = diff_text_time ( elapse_start, elapse_end )
      totelp = totelp + elapsed
!
!  Update the control variable dump file
!
      close (88,iostat=ios)
      Print *,'** Updating control variables dump file ..'
      open (88,file=dumpf(1:ldf),form='UNFORMATTED',status='UNKNOWN',   &
 &         iostat=ios)
      write (88) mxbin_nor,ch1,ch2,real(totelp,r4),real(totcpu,r4), &
 &         no_crash,dfilename
      close (88,iostat=ios)
!
    end do                         ! On bands
!
    if(mxbin_nor < 1) then
      mxbin_nor = nor                   ! Set # of rec. per bin
      Call re_write_header1(mxbin_nor,l2pc_lu,ier)
      if(ier /= 0) goto 888
    end if
!
  end do                         ! multiple key basis (time_i)
!
  Ier = 0
  Call close_l2pc_xx(l2pc_lu,l2pc_lu_key,Keys,rec_nos,nor)
  goto 888
!
777 Ier = 2001
  Print *,'** Error in reading dump file: ',dumpf(1:ldf)
  if(ios > 0) Call ErrMsg(' ',ios)
!
888  close (l2pc_lu,iostat=ios)
     close (l2pc_lu_key,iostat=ios)
!
  if(ier == 0) then
    Lines(1) = ' '
    Lines(2) = ' '
    call duration_text ( totelp, lines(1) )
    call duration_text ( totcpu, lines(1) )
    write (*,920) trim(Lines(1)),trim(Lines(2)),no_crash
  else if(ier > 1000) then
    Call Delete_L2pc_Files(filename)
  end if
!
  close (86,status='delete',iostat=ios)
  close (88,status='delete',iostat=ios)
!
910 format(/,' ** Warning: The triangular basis of: ',a,/,4x,          &
 &         'Does not fall on the pre-selected integration grid !')
!
920 format(/,' ** Total ELAPSED time for this run: ',a,/,              &
 &         '    Total   CPU   time for this run: ',a,/,                &
 &         '    Number of times program crashed: ',i2)
!
  Stop
contains
!
!---------------------------------------------------------------------
  SUBROUTINE GTI(MSG,TI)
!
    CHARACTER(LEN=*),INTENT(IN) :: MSG
    CHARACTER(LEN=*),INTENT(OUT) :: TI
!
    INTEGER :: I,J,L
!
    L = LEN(TI)
    TI(1:L) = ' '
    I = LEN_TRIM(MSG)
!
    WRITE(*,900,ADVANCE='NO') MSG(1:I)
900 FORMAT(A,': ')
!
    READ(5,'(A)',IOSTAT=J) TI
    J = LEN_TRIM(TI)
    IF(J > 1) TI = ADJUSTL(TI)
!
    RETURN
!
  END SUBROUTINE GTI
! --------------------------------------     DELETE_L2PC_FILES     -----
  Subroutine DELETE_L2PC_FILES ( FILENAME )
!  Delete the l2pc data file and its key file
!
    Character(len=*), intent(inout) :: FILENAME
!
    Integer :: j,io
!
    close(89,iostat=io)
    open(89,file=filename,status='unknown',iostat=io)
    close(89,status='delete',iostat=io)

    j = Index(filename,'.dat')
    filename(j+1:j+3) = 'key '
    open(89,file=filename,status='unknown',iostat=io)
    close(89,status='delete',iostat=io)
!
    Return
  End Subroutine DELETE_L2PC_FILES
!
!---------------------------------------------------------------------
! Reading any p_vs_h file:

SUBROUTINE get_pvsh(pvsh,h_grid,t_grid,p_grid,no_hts,ier)

CHARACTER (LEN=*), INTENT(IN) :: pvsh

Real(r8), INTENT(OUT) :: h_grid(*), t_grid(*), p_grid(*)

INTEGER(i4), INTENT(OUT) :: no_hts
INTEGER(i4), INTENT(OUT) :: ier

INTEGER(i4) :: io, uph
CHARACTER (LEN=8) :: ax
REAL(r8) :: h,t,p,plog,sgn

  ier = 0
  uph = p_vs_h_unit
  CLOSE(uph,IOSTAT=io)
  OPEN(uph,FILE=pvsh,STATUS='OLD',IOSTAT=io)
  IF(io /= 0) THEN
    ier = 1
    PRINT *,' ** Error in subroutine: Get_PvsH'
    PRINT *,'    Couldn''T OPEN FILE: ',pvsh
    RETURN
  END IF

  no_hts = 0
  sgn = 1.0
  READ(uph,'(A)',IOSTAT=io) ax
  IF(io == 0) READ(uph,'(A)',IOSTAT=io) ax
  DO WHILE(io == 0)
    READ(uph,*,IOSTAT=io,END=10) h,t,p,plog
    IF(io /= 0) THEN
      ier = 1
      PRINT *,' ** Error in subroutine: Get_PvsH'
      PRINT *,'    Error while reading file: ',pvsh
      GO TO 10
    END IF
    IF(no_hts == 0.AND.plog > 0.0) sgn = -1.0
    no_hts = no_hts + 1
    h_grid(no_hts) = h
    t_grid(no_hts) = t
    p_grid(no_hts) = plog * sgn
  END DO

10  CLOSE(uph,IOSTAT=io)

  RETURN
END SUBROUTINE get_pvsh

End Program L2PC_FILE
! $Log$
! Revision 1.2  2000/07/06 00:11:44  zvi
!  This is the Freeze version of Jun/24/2000
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
