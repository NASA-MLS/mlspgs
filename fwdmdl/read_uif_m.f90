module READ_UIF_M
  use L2PC_FILE_PARAMETERS, only: DEG2RAD, MAX_NO_BANDS, &
                                  MAX_NO_CHANNELS_PER_BAND, &
                                  MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  MAX_NO_POINTINGS, MAX_TABLE_2D
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE, L2PC_HEADER_TWO, &
                                  L2PC_HEADER_TRI
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                                 LIMB_PRESS, MAXAITKENPTS, MAXGEOM, &
                                 MAXGEOPHYS, MAXLINES, MAXPFACH, MAXRAT, &
                                 PFA_SLAB
  use L2PCDim, only: NLVL, NSPS
  use MACHINE, only: IO_ERROR
  use MLSCommon, only: I4, R4
  use S_BRKLINE_M, only: BRKLINE
  use STRINGS, only: LEFTJ, STRLWR, STRUPR
  implicit NONE
  private
  public :: READ_UIF

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

!----------------------------------------------------------------------
!
  Subroutine READ_UIF ( FN, HEADER1, HEADER2, HEADER3, PTG_PRESS,          &
 &           GEOMETRIC, NO_GEOM, GEOPHYSIC, NO_GEOPHYS, ATMOSPHERIC,       &
 &           NO_ATMOS, CONV_HTS, NO_CONV_HTS, FFT_PTS, CH1, CH2,           &
 &           L2PC_LU, PFA, INDIR, FND, AAAP, NO_PFA_CH, NO_FILT_PTS,       &
 &           PFA_CH, NO_INT_FRQS, ACC, TH, PFA_SPECTRUM, PRIMAG, P_INDX,   &
 &           CONV_INDX, SI, N_LVLS, C_YAW, S_YAW, C_ROLL, S_ROLL, C_PITCH, &
 &           S_PITCH, ELEV_183, ELEV_205, AZIM_183, AZIM_205, AZIM_REF,    &
 &           GEOCSRAD, ATMOS_INDEX, GEOM_INDEX, GEOPHYS_INDEX, DO_CONV,    &
 &           L2PC_REC_LENGTH, IER )

    character(len=*), intent(in) :: FN
    type (l2pc_header_one), intent(out) :: HEADER1
    type (l2pc_header_two), intent(out) :: HEADER2
    type (l2pc_header_tri), intent(out) :: HEADER3
    type (geom_param), intent(out) :: GEOMETRIC(*)
    type (limb_press), intent(out) :: PTG_PRESS
    integer(i4), intent(out) :: NO_GEOM
    type (geophys_param), intent(out) :: GEOPHYSIC(*)
    integer(i4), intent(out) :: NO_GEOPHYS
    type (atmos_comp), intent(out) :: ATMOSPHERIC(*)
    integer(i4), intent(out) :: NO_ATMOS
    real(r4), intent(out) :: CONV_HTS(*)
    integer(i4), intent(out) ::NO_CONV_HTS
    integer(i4), intent(out) :: FFT_PTS
    integer(i4), intent(out) :: CH1, CH2
    integer(i4), intent(in) :: L2PC_LU
    integer(i4), intent(out) :: PFA
    character(len=*), intent(out) :: INDIR
    character(len=*), intent(out) :: FND
    character(len=*), intent(out) :: AAAP
    integer(i4), intent(out) :: NO_PFA_CH
    integer(i4), intent(out) :: NO_FILT_PTS
    integer(i4), intent(out) :: PFA_CH(*)
    integer(i4), intent(out) ::NO_INT_FRQS(*)
    real(r4), intent(out) :: ACC(*)
    real(r4), intent(out) :: TH(*)
    type (pfa_slab), intent(out) :: PFA_SPECTRUM(6,*)
    character, intent(out) :: PRIMAG
    integer(i4), intent(out) :: P_INDX(*)
    integer(i4), intent(out) :: CONV_INDX(*)
    integer(i4), intent(out) :: SI
    integer(i4), intent(out) :: N_LVLS
    real(r4), intent(out) :: C_YAW, S_YAW
    real(r4), intent(out) :: C_ROLL, S_ROLL
    real(r4), intent(out) :: C_PITCH, S_PITCH
    real(r4), intent(out) :: ELEV_183, ELEV_205
    real(r4), intent(out) :: AZIM_183, AZIM_205, AZIM_REF
    real(r4), intent(out) :: GEOCSRAD
    integer(i4), intent(out) :: ATMOS_INDEX(*)
    integer(i4), intent(out) :: GEOM_INDEX(*)
    integer(i4), intent(out) :: GEOPHYS_INDEX(*)
    Logical*1, intent(out) :: DO_CONV
    integer(i4), intent(out) :: L2PC_REC_LENGTH
    integer(i4), intent(out) :: IER
!
    character(len=8) :: AN1
    integer(i4) :: BAND
    character(len=12) :: BN
    Logical :: DO_DER(6), DO_FWD(6)
    integer(i4) :: I, II
    integer(i4) :: IZ
    integer(i4) :: J, JJ
    integer(i4) :: K
    integer(i4) :: L
    character(len=100) :: LINE
    character(len=80) :: MDBD
    integer(i4) :: NC
    integer(i4) :: NO_SPS
    character(len=8) :: NAME
    integer(i4) :: NV
    integer(i4) :: PBAND
    real(r4) :: R
    integer(i4) :: UU
    real(r4) :: V(5)
!
! -----     Statement functions     ------------------------------------
!
! These statement functions compensate for the use of COSD and SIND,
! which are intrinsic in some compilers, but are not part of any standard.
!
    real(r4) :: COSD, SIND, X
    cosd(x) = cos(x*deg2rad)
    sind(x) = sin(x*deg2rad)
!
! -----     begin code     ---------------------------------------------
!
    Ier = 0
    Uu = 45
    Print *
    Primag = '@'
    Close (Uu,iostat=iz)
    Open (Uu,file=Fn,status='OLD',iostat=iz)
    if (iz /= 0) then
      Print *,'* Error: Could not open File:'
      Print *,'   ',trim(Fn(1:j))
      goto 99
    end if
!
    si = 0
    pfa = 1
    ch1 = 0
    ch2 = 0
    band = 0
    Bn(1:)=' '
    no_geom = 0
    no_atmos = 0
    no_geophys = 0
    do_conv = .true.
    l2pc_rec_length = 0
!
    InDir(1:) = ' '
    Fnd(1:) = ' '
    Aaap(1:) = ' '
    mdbd(1:) = ' '
!
    header2%no_sv_elmnts = 0
    header1%no_avail_keys = 0
    ptg_press%no_lin_values = 0
    header1%b_phi_lin_val = 0.0
    header1%no_b_theta_lin_val = 0
    header1%b_theta_lin_val(1) = -1.0
!
    header1%no_bands = max_no_bands
    header1%no_k_records_per_bin = 0
    header1%no_pointings = max_no_pointings
    header1%no_channels_per_band = max_no_channels_per_band
    header1%no_coeff_per_component = max_no_elmnts_per_sv_component
!
! Initialize the header3 record:
!
    An1(1:)=' '
    header3%matdim = 0
    do j = 1, max_table_2d
      header3%second_der_matrix_namid(j) = An1
      header3%second_der_matrix_bands(j) = An1
    end do
!
    iz = 0
    do while (iz == 0)
!
      Read (Uu,'(A)',iostat=iz) Line
      if (iz < 0) goto 30
      if (iz > 0) goto 99
!
      Call Leftj(Line)
      Call StrUpr(Line)
!
      if (index(Line,'NO_CONV_HTS') > 0) Line(1:20) = 'NO_CONV_INDICES'
!
      if (Line(1:14) == 'HEADER_3_LINES') then
        Read (Uu,'(A)',iostat=iz) header1%line1
        if (iz == 0) Read (Uu,'(A)',iostat=iz) header1%line2
        if (iz == 0) Read (Uu,'(A)',iostat=iz) header1%line3
        if (iz /= 0) goto 99
      else if (Line(1:09) == 'INPUT_DIR') then
        Read (Uu,'(A)',iostat=iz) Line
        if (iz /= 0) goto 99
        Call Leftj(Line)
        j = len_trim(Line)
        if (Line(j:j) /= '/') Line(j+1:j+1)='/'
        InDir = Line(1:j+1)
      else if (Line(1:19) == 'MASTER_DATABASE_DIR') then
        Read (Uu,'(A)',iostat=iz) Line
        if (iz /= 0) goto 99
        Call Leftj(Line)
        j = len_trim(Line)
        if (Line(j:j) /= '/') Line(j+1:j+1)='/'
        mdbd = Line(1:j+1)
      else if (Line(1:16) == 'PRIMARY OR IMAGE') then
        Read (Uu,'(A)',iostat=iz) Line
        if (iz /= 0) goto 99
        Call Leftj(Line)
        Call StrLwr(Line)
        Primag = Line(1:1)
      else if (Line(1:23) == 'NO_OF_CHANNELS_PER_BAND') then
        Read (Uu,*,iostat=iz) j
        if (iz /= 0) goto 99
        if (j > 0) header1%no_channels_per_band = j
      else if (Line(1:23) == 'MAX_NO_OF_COEFF_PER_SPS') then
        Read (Uu,*,iostat=iz) j
        if (iz /= 0) goto 99
        if (j > 0) header1%no_coeff_per_component = j
      else if (Line(1:8) == 'CONVOLVE') then
        Read (Uu,'(A)',iostat=iz) Line
        if (iz /= 0) goto 99
        Call Leftj(Line)
        Call StrUpr(Line)
        do_conv = (Line(1:1) == 'Y')
      else if (Line(1:14) == 'CHANNELS_RANGE') then
        Read (Uu,*,iostat=iz) ch1,ch2
        if (iz /= 0) goto 99
        if (ch1 > ch2) then
          j = ch1
          ch1 = ch2
          ch2 = j
        end if
!
!  *** Special code for the PC (AbSoft F77) version: No more than ONE band
!      could be run at one time !!
!
        k = max_no_channels_per_band
        band = (ch1 + k - 1) / k
        if (ch2-ch1+1 > k) then
          ch1 = 1 + (band - 1) * k
          ch2 = min(ch2,ch1+k-1)
        end if
!
!  *** End of Special code
!
      else if (Line(1:12) == 'LOG2_FFT_PTS') then
        Read (Uu,*,iostat=iz) fft_pts
        if (iz /= 0) goto 99
      else if (Line(1:16) == 'NO_FILTER_POINTS') then
        Read (Uu,*,iostat=iz) no_filt_pts
        if (iz /= 0) goto 99
      else if (Line(1:16) == 'ANTENNA_FILENAME') then
        Read (Uu,'(A)',iostat=iz) Aaap
        if (iz /= 0) goto 99
        Call Leftj(Aaap)
      else if (Line(1:09) == 'PTG_PRESS') then
        Read (Uu,'(A)',iostat=iz) Line
        if (iz /= 0) goto 99
        Call Leftj(Line)
        j = index(Line,'.')
        k = min(8,j-1)
        ptg_press%name = Line(1:k)
        Call StrUpr(ptg_press%name)
        Line(1:j-1) = ' '
        Call Leftj(Line)
        Read (Line,*,iostat=iz) (ptg_press%der_calc(j),j=1,6)
        if (iz /= 0) goto 99
        Read (Uu,*,iostat=iz) k,(ptg_press%lin_val(j),j=1,k)
        if (iz /= 0) goto 99
        ptg_press%no_lin_values = k
      else if (Line(1:11) == 'B_PHI_VALUE') then
        Read (Uu,*,iostat=iz) header1%b_phi_lin_val
        if (iz /= 0) goto 99
      else if (Line(1:20) == 'NO_B_THETA_LIN_VALUE') then
        Read (Uu,*,iostat=iz) header1%no_b_theta_lin_val
        if (iz /= 0) goto 99
      else if (Line(1:17) == 'B_THETA_LIN_VALUE') then
        k = header1%no_b_theta_lin_val
        if (k < 1) then
          iz = -6
          Line = 'NO_B_THETA_LIN_VALUE is missing or misplaced ..'
          goto 99
        end if
        Read (Uu,*,iostat=iz) (header1%b_theta_lin_val(j),j=1,k)
        if (iz /= 0) goto 99
      else if (Line(1:12) == 'NO_GEOMETRIC') then
        Read (Uu,*,iostat=iz) no_geom
        if (iz /= 0) goto 99
      else if (Line(1:09) == 'GEOMETRIC') then
        if (no_geom < 1) then
          iz = -6
          Line = 'NO_GEOMETRIC is missing or misplaced ..'
          goto 99
        else if (no_geom > maxgeom) then
          iz = -6
          Line = 'NO_GEOMETRIC too large, must be < maxgeom'
          goto 99
        end if
        do i = 1, no_geom
          Read (Uu,'(A)',iostat=iz) Line
          if (iz /= 0) goto 99
          Call Leftj(Line)
          k = Index(Line,' ')
          geometric(i)%name = Line(1:k-1)
          Call StrUpr(geometric(i)%name)
          Line(1:k) = ' '
          Call Leftj(Line)
          Read (Line,*,iostat=iz) (geometric(i)%der_calc(j),j=1,6)
          if (iz /= 0) goto 99
          geometric(i)%lin_val = r
        end do
      else if (Line(1:14) == 'NO_GEOPHYSICAL') then
        Read (Uu,*,iostat=iz) no_geophys
        if (iz /= 0) goto 99
      else if (Line(1:11) == 'GEOPHYSICAL') then
        if (no_geophys < 1) then
          iz = -6
          Line = 'NO_GEOPHYSICAL is missing or misplaced ..'
          goto 99
        else if (no_geophys > maxgeophys) then
          iz = -6
          Line = 'NO_GEOPHYSICAL is too large, must be < maxgeophys'
          goto 99
        end if
        do i = 1, no_geophys
          Read (Uu,'(A)',iostat=iz) Line
          if (iz /= 0) goto 99
          Call Leftj(Line)
          k = Index(Line,' ')
          geophysic(i)%name = Line(1:k-1)
          Call StrUpr(geophysic(i)%name)
          Line(1:k) = ' '
          Call Leftj(Line)
          Read (Line,*,iostat=iz) (geophysic(i)%der_calc(j),j=1,6)
          if (iz /= 0) goto 99
          Read (Uu,*,iostat=iz) k,(geophysic(i)%basis_peaks(j),j=1,k)
          if (iz /= 0) goto 99
          geophysic(i)%no_lin_values = k
        end do
      else if (Line(1:14) == 'NO_ATMOSPHERIC') then
        Read (Uu,*,iostat=iz) no_atmos
        if (iz /= 0) goto 99
      else if (Line(1:11) == 'ATMOSPHERIC') then
        if (no_atmos < 1) then
          iz = -6
          Line = 'NO_ATMOSPHERIC is missing or misplaced ..'
          goto 99
        else if (band < 1) then
          iz = -6
          Line = 'CHANNELS_RANGE must be given before ATMOSPHERIC ..'
          goto 99
        end if
        pband = no_atmos
        no_atmos = 0
        do i = 1, pband
          j = -1
          do while (j < 1)
            Read (Uu,'(A)',iostat=iz) Line
            if (iz /= 0) goto 99
            j = Index(Line,'.f.')
            if (j < 1) j = Index(Line,'.t.')
          end do
          Call Leftj(Line)
          k = Index(Line,' ')
          Name = Line(1:k-1)
          Call StrUpr(Name)
          Line(1:k) = ' '
          Call Leftj(Line)
          Read (Line,*,iostat=iz) (do_fwd(j),j=1,6)
          if (iz /= 0) goto 99
          Read (Uu,*,iostat=iz) (do_der(j),j=1,6)
          if (iz /= 0) goto 99
          if (do_fwd(band).or.do_der(band)) then
            no_atmos = no_atmos + 1
            if (no_atmos > Nsps) then
              iz = -6
              Line = 'NO_ATMOSPHERIC too large, must be <= Nsps'
              goto 99
            end if
            atmospheric(no_atmos)%name = Name
            do j = 1, 6
              atmospheric(no_atmos)%fwd_calc(j) = do_fwd(j)
              atmospheric(no_atmos)%der_calc(j) = do_der(j)
            end do
            Read (Uu,*,iostat=iz) k,jj,                                 &
   &                (atmospheric(no_atmos)%basis_peaks(j),j=1,jj)
            if (iz /= 0) goto 99
            atmospheric(no_atmos)%spectag = k
            atmospheric(no_atmos)%no_lin_values = jj
          end if
        end do
      else if (Line(1:08) == 'NO_PRESS') then
        Read (Uu,*,iostat=iz) N_lvls
        if (iz /= 0) goto 99
      else if (Line(1:13) == 'PRESS_INDICES') then
        Read (Uu,*,iostat=iz) (p_indx(j),j=1,N_lvls)
        if (iz /= 0) goto 99
      else if (Line(1:15) == 'NO_CONV_INDICES') then
        Read (Uu,*,iostat=iz) no_conv_hts
        if (iz /= 0) goto 99
      else if (Line(1:12) == 'CONV_INDICES') then
        Read (Uu,*,iostat=iz) (conv_indx(j),j=1,no_conv_hts)
        if (iz /= 0) goto 99
      else if (Line(1:13) == 'SURFACE_INDEX') then
        Read (Uu,*,iostat=iz) si
        if (iz /= 0) goto 99
      else if (Line(1:17) == 'HTS_BELOW_SURFACE') then
        if (si < 1) then
          iz = -6
          Line = 'SURFACE_INDEX is missing or misplaced ..'
          goto 99
        end if
        Read (Uu,*,iostat=iz) (conv_hts(j),j=1,si)
        if (iz /= 0) goto 99
      else if (Line(1:07) == 'NO_BINS') then
        Read (Uu,*,iostat=iz) header1%no_avail_keys
        if (iz /= 0) goto 99
      else if (Line(1:04) == 'BINS') then
        j = 0
        ii = header1%no_avail_keys
        if (ii < 1) then
          iz = -6
          Line = 'NO_BINS is missing or misplaced ..'
          goto 99
        end if
        do while (j < ii)
          l = 1
          Read (Uu,'(A)',iostat=iz) Line
          if (iz /= 0) goto 99
          do while (j < ii.and.l > 0)
            Call Leftj(Line)
            l = len_trim(Line)
            if (l > 0) then
              k = index(Line,' ')
              if (k < 1) k = l + 1
              j = j + 1
              header1%avail_keys(j) = Line(1:k-1)
              Line(1:k) = ' '
            end if
          end do
        end do
      else if (Line(1:09) == 'NO_PFA_CH') then
        Read (Uu,*,iostat=iz) no_pfa_ch
        if (iz /= 0) goto 99
        if (no_pfa_ch > maxpfach) then
          iz = -6
          Line = 'NO_PFA_CH too large, must be < maxpfach'
          goto 99
        end if
      else if (Line(1:06) == 'PFA_CH') then
        if (no_pfa_ch > 0) then
          j = 0
          An1(1:)=' '
          write (An1,'(i4)') maxaitkenpts
          Call Leftj(An1)
          do l = 1, no_pfa_ch
            Read (Uu,'(A)',iostat=iz) Line
            if (iz /= 0) goto 99
            i = index(Line,'!')
            if (i > 0) Line(i:)=' '
            i = len_trim(Line)
            Line(i+1:) = ' 10.0'
            Call Brkline(Line,V,nv,iz)
            if (iz /= 0) goto 99
            k = Int(V(1)+0.2)
            if (k >= ch1.and.k <= ch2) then
              j = j + 1
              pfa_ch(j) = k
              no_int_frqs(j) = Int(V(2)+0.2)
              Th(j) = max(0.1,V(3))
              Acc(j) = V(4)
              i = no_int_frqs(j)
              jj = 4 * i + 1
              if (jj > maxaitkenpts) then
                i = (maxaitkenpts - 1) / 4
                Print *,'** Warning: PFA channel:',k
                Print *,'   Maximum Aitkens points exceeded ',An1
                Print *,'   no_int_frqs()  set to:',i
              end if
              no_int_frqs(j) = i
            else
              Print *,'** Warning: PFA channel:',k
              write (Bn,'(i2,''  to: '',i2)') ch1,ch2
              Print *,'   Outside of channels range: ',Bn(1:10)
              Print *,'   This channel will not be processed !!'
            end if
          end do
          no_pfa_ch = j
        end if
      end if
    end do
!
!  Compose pointing information:
!
30 k = ptg_press%no_lin_values
   if (k < 1) then
     iz = -6
     Line = 'Number of PTAN63 pointings is missing or misplaced ..'
     goto 99
   end if
!
   header1%no_pointings = k
   header1%pointings(1:k) = ptg_press%lin_val(1:k)
!
!  Compute l2pc file max. record length (Computed based on the K-star
!  type record which is the longest. Please see write_one_record.f
!  subroutine for detail):
!
    j = header1%no_pointings * header1%no_channels_per_band
    l2pc_rec_length = 68 + 4 * j     ! In bytes (NOT words !!)
!
    if (Index('pppiii',Primag) < 1) then
      iz = -6
      Line = 'UNKNOWN RUN TYPE ! Must be P (Primary) or I (Image) ..'
      goto 99
    end if
!
    nc = 0
    nv = 0
    if (any(ptg_press%der_calc(1:6))) then
      nv = nv + 1
      nc = header2%no_sv_elmnts + 1
      header2%no_sv_elmnts = nc
      header1%sv_components(nv) = ptg_press%name
      header1%sv_component_first_elmnt_index(nv) = nc
      header1%no_elmnts_per_sv_component(nv) = 1
      do j = 1, 6
        header1%sv_rtrvl_by_band(1,j) = ptg_press%der_calc(j)
      end do
      header2%tri_basis_vert_grid(nv) = 0.0
      header1%no_sv_components = nv
    else
      header1%no_sv_components = 0
      header2%no_sv_elmnts     = 0
    end if
!
    header1%no_mag_fields = 0
    header1%no_b_theta_lin_val = 0
!
    if (no_geom < 1) then
      iz = -6
      Line = 'NO_GEOMETRIC is missing or misplaced ..'
      goto 99
    end if
!
!  Compose geometric information:
!
    do i = 1, no_geom
      geom_index(i) = 0
      if (any(geometric(i)%der_calc(1:6))) then
        nv = nv + 1
        geom_index(i) = nv
        header1%no_sv_components = nv
        header1%sv_components(nv) = geometric(i)%name
        nc = header2%no_sv_elmnts + 1
        header2%no_sv_elmnts = nc
        header1%sv_component_first_elmnt_index(nv) = nc
        header1%no_elmnts_per_sv_component(nv) = 1
        header1%sv_rtrvl_by_band(nv,1:6) = geometric(i)%der_calc(1:6)
        j = header1%sv_component_first_elmnt_index(nv)
        header2%tri_basis_vert_grid(j) = 0.0
      end if
    end do
!
    if (no_geophys < 1) then
      iz = -6
      Line = 'NO_GEOPHYSICAL is missing or misplaced ..'
      goto 99
    end if
!
!  Compose geophysical information:
!
    do i = 1, no_geophys
      geophys_index(i) = 0
      if (any(geophysic(i)%der_calc(1:6))) then
        nv = nv + 1
        geophys_index(i) = nv
        header1%no_sv_components = nv
        header1%sv_components(nv) = geophysic(i)%name
        header1%sv_rtrvl_by_band(nv,1:6) = geophysic(i)%der_calc(1:6)
        nc = header2%no_sv_elmnts + 1
        header1%sv_component_first_elmnt_index(nv) = nc
        jj = geophysic(i)%no_lin_values
        header1%no_elmnts_per_sv_component(nv) = jj
        header2%no_sv_elmnts = nc + jj - 1
        ii = nc - 1
        do j = 1, jj
          ii = ii + 1
          header2%tri_basis_vert_grid(ii)=geophysic(i)%basis_peaks(j)
        end do
      end if
    end do
!
    if (no_atmos < 1) then
      iz = -6
      Line = 'NO_ATMOSPHERIC is missing or misplaced ..'
      goto 99
    end if
!
!  Compose atmospheric information:
!
    do i = 1, no_atmos
      atmos_index(i) = 0
      if (any(atmospheric(i)%der_calc(1:6))) then
        nv = nv + 1
        atmos_index(i) = nv
        header1%no_sv_components = nv
        header1%sv_components(nv) = atmospheric(i)%name
        header1%sv_rtrvl_by_band(nv,1:6) = atmospheric(i)%der_calc(1:6)
        nc = header2%no_sv_elmnts + 1
        header1%sv_component_first_elmnt_index(nv) = nc
        jj = atmospheric(i)%no_lin_values
        header1%no_elmnts_per_sv_component(nv) = jj
        header2%no_sv_elmnts = nc + jj - 1
        ii = nc - 1
        do j = 1, jj
          ii = ii + 1
          header2%tri_basis_vert_grid(ii) = atmospheric(i)%basis_peaks(j+1)
        end do
      end if
    end do
!
    do j = 1, no_geom                ! Set geometric parameters
      Name = geometric(j)%name
      r = geometric(j)%lin_val
      if (Name == 'ROLL') then
        c_roll = cosd(r)
        s_roll = sind(r)
      else if (Name == 'PITCH') then
        c_pitch = cosd(r)
        s_pitch = sind(r)
      else if (Name == 'YAW') then
        c_yaw = cosd(r)
        s_yaw = sind(r)
      else if (Name == 'ELEV_183') then
        elev_183 = r
      else if (Name == 'ELEV_205') then
        elev_205 = r
      else if (Name == 'AZIM_183') then
        azim_183 = r
      else if (Name == 'AZIM_205') then
        azim_205 = r
      else if (Name == 'AZIM_REF') then
        azim_ref = r
      else if (Name == 'GEOCSRAD') then
        geocsrad = r
      end if
    end do
!
    jj = no_conv_hts+si-1
    if (jj > Nlvl) then
      iz = -6
      Line = 'NO_CONV_HTS exceeds maximum (zz), reduce & re-run ..'
      i = index(Line,'zz')
      write (Line(i:i+1),'(i2.2)') Nlvl
      goto 99
    end if
!
! Reformat the key code
!
    iz = 0
    do j = 1, header1%no_avail_keys
      Call StrUpr(header1%avail_keys(j))
      header1%avail_keys(j)(8:8) = '_'
    end do
!
99  Close (Uu,iostat=j)
!
    if (iz == 0) then
      if (len_trim(InDir) < 1) then
        iz = -6
        Line = 'INPUT_DIR is missing or misplaced ..'
      end if
    end if
!
    if (iz /= 0) then
      Print *
      if (iz < -2) then
        Print *,'** Error in subroutine Read_Uif'
        Print *,'   ',trim(Line(1:i))
      else
        Call io_error('From subroutine Read_Uif',iz)
      end if
      Ier = 1
      Return
    end if
!
    if (len_trim(mdbd) < 1) then
      Print *,'** Missing or incorrect MASTER_DATABASE_DIR name ..'
      Print *,'   Set to default: /zvi/data/'
      mdbd = '/zvi/data/'
      iz = len_trim(mdbd)
    end if
!
    Fnd = mdbd(1:iz)//'eos_master_database.dat'
!
!  Clear the pfa_spectrum record:
!
    do band = 1, 6
      do i = 1, no_atmos
        pfa_spectrum(band,i)%no_sps = 0
        pfa_spectrum(band,i)%no_lines = 0
        pfa_spectrum(band,i)%sps_spectag = 0
        do l = 1, Nlvl
          pfa_spectrum(band,i)%varM(l) = 0.0d0
          do k = 1, maxrat
            pfa_spectrum(band,i)%Xx(k,l) = 0.0d0
            pfa_spectrum(band,i)%Yy(k,l) = 0.0d0
            pfa_spectrum(band,i)%Dy(k,l) = 0.0d0
          end do
        end do
        do l = 1, maxlines
          pfa_spectrum(band,i)%sps_w(l) = 0.0d0
          pfa_spectrum(band,i)%sps_n(l) = 0.0d0
          pfa_spectrum(band,i)%sps_ps(l) = 0.0d0
          pfa_spectrum(band,i)%sps_n1(l) = 0.0d0
          pfa_spectrum(band,i)%sps_n2(l) = 0.0d0
          pfa_spectrum(band,i)%sps_v0(l) = 0.0d0
          pfa_spectrum(band,i)%sps_el(l) = 0.0d0
          pfa_spectrum(band,i)%sps_str(l) = 0.0d0
          pfa_spectrum(band,i)%sps_gamma(l) = 0.0d0
          pfa_spectrum(band,i)%sps_delta(l) = 0.0d0
          do ii = 1, 3
            pfa_spectrum(band,i)%sps_part(ii,l) = 0.0d0
          end do
        end do
      end do
    end do
!
    if (no_pfa_ch < 1) Return
!
! Make sure PFA channels are in order, sort them and their attachments
!
    if (no_pfa_ch > 1) then
      do j = 2, no_pfa_ch
        i = -1
        do k = no_pfa_ch, j, -1
          if (pfa_ch(k) < pfa_ch(k-1)) then
            i = pfa_ch(k)
            pfa_ch(k) = pfa_ch(k-1)
            pfa_ch(k-1) = i
            i = no_int_frqs(k)
            no_int_frqs(k) = no_int_frqs(k-1)
            no_int_frqs(k-1) = i
            r = Th(k)
            Th(k) = Th(k-1)
            Th(k-1) = r
            r = Acc(k)
            Acc(k) = Acc(k-1)
            Acc(k-1) = r
          end if
        end do
        if (i < 1) goto 50
      end do
    end if
!
! Load the pfa structures now:
!
50  pband = -1
    do ii = 1, no_pfa_ch
      jj = pfa_ch(ii)
      band = (jj + 14) / 15
      if (band /= pband) then
        k = 0
        pband = band
        do i = 1, no_atmos
          An1 = atmospheric(i)%name
          j = atmospheric(i)%spectag
          if (atmospheric(i)%fwd_calc(band)) then
            k = k + 1
            pfa_spectrum(band,k)%sps_name = An1
            pfa_spectrum(band,k)%sps_spectag = j
          end if
        end do
        no_sps = k
        if (no_sps > 0) pfa_spectrum(band,1:no_sps)%no_sps = no_sps
      end if
    end do
!
    Return
  End Subroutine READ_UIF
end module READ_UIF_M

! $Log$
