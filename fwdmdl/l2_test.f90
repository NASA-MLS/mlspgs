Program l2_Test
  use GL6P, only: NG
  use MLSCommon, only: I4, R4, R8
  use STRINGS, only: STRLWR
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  DEG2RAD
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, MAXFILTPTS, &
                                 PFA_SLAB, SPECTRO_PARAM, MAXPFACH, &
                                 MAXGEOM, LIMB_PRESS, K_MATRIX_INFO
  use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi, &
                     MNM => max_no_mmaf
  use ELLIPSE, only: PHI_TAN, ROC
  use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_FILTERS_M, only: GET_FILTERS
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use RAD_TRAN_M, only: RAD_TRAN
  use RAD_TRAN_WD_M, only: RAD_TRAN_WD
  use FREQ_AVG_M, only: FREQ_AVG
  use CONVOLVE_ALL_M, only: CONVOLVE_ALL
  use D_HUNT_M, only: HUNT          ! ** DEBUG
  implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl

Integer(i4), PARAMETER :: mnf = 75

Logical :: temp_der

Integer(i4) :: p_indx(Nlvl), SPECT_ATMOS(Nsps), no_ptg_frq(Nptg), &
               no_coeffs_f(Nsps), no_phi_f(Nsps), SPECT_INDEX(Nsps), &
               no_spectro, pfa_ch(MAXPFACH), t_indx(Nptg)

Integer(i4) :: i, j, k, kk, ht_i, no_t, mnz, no_geom, no_tan_hts, kz, &
               ld, ch, Spectag, no_freqs, no_pfa_ch, prev_npf, n_obs, &
               no_atmos, ch1, ch2, n_lvls, si, mfi, jp, band, n_sps,  &
               ptg_i, frq_i, io, klo, khi, l, n, brkpt, no_ele, nl,   &
               m, mid, ilo, ihi, ier, no_filt_pts, no_phi_t, fft_pts, &
               k_info_count,  no_mmaf, gl_count

Type(path_index)  :: ndx_path(Nptg,mnm)
Type(path_vector) :: z_path(Nptg,mnm),t_path(Nptg,mnm),h_path(Nptg,mnm),  &
                     dhdz_path(Nptg,mnm), spsfunc_path(Nsps,Nptg,mnm),    &
                     n_path(Nptg,mnm),phi_path(Nptg,mnm)
Type(path_vector) :: ptg_frq_grid(Nptg)

Type(path_derivative) :: dh_dt_path(Nptg,mnm)

Real(r8) :: z_gnlv(400),thbs(10),phi_tan_mmaf(mnm),elev_offset
Real(r8) :: href(Nlvl),zref(Nlvl),t_z_basis(mxco),t_script(N2lvl),  &
            t_coeff(mxco,mnm),ref_corr(N2lvl,Nptg),tau(N2lvl),Qlog(3), &
            tan_dh_dt(Nlvl,mnm,mxco),t_phi_basis(mnp),t_phi_basis_copy(mnp)

Real(r8) :: dx_dt(Nptg,mxco), d2x_dxdt(Nptg,mxco)

Real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
Real(r8) :: z_grid(Nlvl),dh_dt_glgrid(ngt,mnp,mxco), dhdz_glgrid(ngt,mnp)

Real(r8) :: ptg_angles(Nptg,mnm), center_angle
Real(r8) :: tan_press(Nptg), tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)

Logical :: IS_F_LOG(Nsps)

real(r8) :: freq_grid(mnf),freqs(Nch)

real(r8) :: MR_F(mxco,mnp,Nsps), F_BASIS(mxco,Nsps)
real(r8) :: F_PHI_BASIS(mnp,Nsps), F_PHI_BASIS_COPY(mnp,Nsps)

real(r8) :: FILTER_FUNC(maxfiltpts,maxpfach)
real(r8) :: F_GRID_FILTER(maxfiltpts,maxpfach)

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

real(r8) :: I_STAR_ALL(Nch)

real(r4) :: K_STAR_ALL(02,20,mxco,mnp,Nptg)      ! 02 should be: Nch
Type(k_matrix_info) :: k_star_info(20)

Type(path_derivative) :: k_temp_frq, k_atmos_frq(Nsps), &
                         k_spect_dw_frq(Nsps), k_spect_dn_frq(Nsps), &
                         k_spect_dnu_frq(Nsps)
!
Type(path_beta) :: beta_path(Nsps,mnf)

Real(r8) :: Radiances(Nptg,Nch)
Real(r8) :: RadV(mnf), f_grid(mnf)
Real(r8) :: s_temp, h_obs, earth_ref, e_rad, zeta, Frq, h_tan, Rad, &
            beta_inc, geoc_lat, q, r

Real(r4) :: elev_183, elev_205
!
Character (LEN=01) :: CA, Primag
Character (LEN=08) :: Name
Character (LEN=16) :: Vname
Character (LEN=80) :: InDir, Aaap, Fnd, Line
Character (LEN=40) :: Ax, Dtm1, Dtm2

Type(limb_press)       :: PTG_PRESS
Type (atmos_comp)      :: ATMOSPHERIC(Nsps)
Type (geom_param)      :: GEOMETRIC(maxgeom)
Type (pfa_slab)        :: PFA_SPECTRUM(Nsps)
Type (spectro_param)   :: SPECTROSCOPIC(3*Nsps)

Real(r8), DIMENSION(:,:), ALLOCATABLE :: s_phi_basis_copy

!  ----------------------
! Read the inputs from a file...
!
  ier = 0
  Fnd(1:) = ' '
  Fnd = '/user5/zvi/mod_seez'
  CLOSE(11,iostat=io)
  OPEN(11,file=Fnd,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99

  p_indx(1:Nlvl) = 0
  t_indx(1:Nptg) = 0
  no_phi_f(1:Nsps) = 0
  no_ptg_frq(1:Nptg) = 0
  pfa_ch(1:MAXPFACH) = 0
  no_coeffs_f(1:Nsps) = 0

  SPECT_INDEX(1:Nsps) = -1
  SPECT_ATMOS(1:Nsps) = -1

  do
    Ax(1:) = ' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    if (Index(Ax,'No_Mmaf') > 0) EXIT
  end do

  read(11,*,iostat=io) j
  if(io /= 0) goto 99

  Primag = 'p'
  no_mmaf = min(j,mnm)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  tau(1:no_mmaf) = 0.0
  read(11,*,iostat=io) (tau(i),i=1,no_mmaf)
  if(io /= 0) goto 99

  phi_tan_mmaf(1:no_mmaf) = tau(1:no_mmaf) * deg2rad

  do
    Ax(1:) = ' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    if (Index(Ax,'Channels_Range') > 0) EXIT
  end do

  read(11,*,iostat=io) ch1, ch2
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_pfa_ch
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (pfa_ch(i),i=1,no_pfa_ch)
  if(io /= 0) goto 99
!
  primag = 'p'
  freqs(1:Nch) = 0.0D0
  DO i = ch1, ch2
    CALL radiometry(i,q,r,zeta,kk)     ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'p') freqs(i) = q     ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'i') freqs(i) = r     ! DEBUG, Added Jan/23/2000, Z.S
  END DO
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) ptg_press%name, (ptg_press%der_calc(i),i=1,6)
  if(io /= 0) goto 99

  read(11,*,iostat=io) j
  if(io /= 0) goto 99

  ptg_press%no_lin_values = j
  read(11,*,iostat=io) (ptg_press%lin_val(i),i=1,j)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_geom
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  do i = 1, no_geom
    read(11,*,iostat=io) geometric(i)
    if(io /= 0) goto 99
    r = geometric(i)%lin_val
    IF(geometric(i)%name == 'ELEV_183') THEN
      elev_183 = r
    ELSE IF(geometric(i)%name == 'ELEV_205') THEN
      elev_205 = r
    ELSE IF(geometric(i)%name == 'EARTHREF') THEN
      earth_ref = r
    ELSE IF(geometric(i)%name == 'SPACE_T') THEN
      s_temp = r
    ELSE IF(geometric(i)%name == 'GEOCSRAD') THEN
      h_obs = r
    END IF
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_atmos
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

! read(11,*,iostat=io) (atmos_index(i),i=1,no_atmos)
  read(11,*,iostat=io) i                              ! ** DUMMY read
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = mxco
  do j = 1, no_atmos
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    atmospheric(j)%NAME = AdjustL(Ax)
    read(11,*,iostat=io) Spectag, ht_i
    if(io /= 0) goto 99
    atmospheric(j)%SPECTAG = Spectag
    atmospheric(j)%NO_LIN_VALUES = ht_i
    read(11,*,iostat=io) (atmospheric(j)%FWD_CALC(i),i=1,6)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (atmospheric(j)%DER_CALC(i),i=1,6)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (atmospheric(j)%LIN_VAL(i),i=1,kk)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (atmospheric(j)%BASIS_PEAKS(i),i=1,kk+2)
    if(io /= 0) goto 99
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) band, n_sps
  if(io /= 0) goto 99

  if (band >= 5) then                       ! 183
    elev_offset = elev_183 * deg2rad
  else if (band >= 2) then                  ! 205
    elev_offset = elev_205 * deg2rad
  else                                      ! 63 / pointing reference
    elev_offset = 0.0
  end if

  m = 0
  do j = 1, no_atmos
    IF(atmospheric(j)%FWD_CALC(band)) THEN
      m = m + 1
      atmospheric(m) = atmospheric(j)
    ENDIF
  end do
  no_atmos = m

  if(no_atmos /= n_sps) then
    io = 3
    Print *,'** Error: New code: n_sps should be equal to no_atmos !'
    Print *,'          n_sps =',n_sps,' no_atmos:',no_atmos
    goto 99
  endif

  m = 0
  pfa_spectrum(1)%NO_SPS = n_sps     ! Make sure we have this

  read(11,'(A)',iostat=io) Ax    ! pfa_spectrum(s)
  if(io /= 0) goto 99

  DO

    if(m == n_sps) then
      do
        read(11,'(A)',iostat=io) Ax
        if(io /= 0) goto 99
        if(Index(Ax,'END_CAT') > 0) EXIT
      end do
      EXIT
    endif

    Line = ' '
    Name = ' '
    read(11,'(A)',iostat=io) Line
    if(io /= 0) goto 99
    if(Index(Line,'END_CAT').gt.0) EXIT
    read(Line,*,iostat=io) Name, Spectag, nl, (Qlog(i),i=1,3)
    if(io /= 0) goto 99

    j = 0
    DO i = 1, n_sps
      if(Name == atmospheric(i)%NAME) then
        j = i
        EXIT
      endif
    END DO

    if(j < 1) then
      do i = 1, nl
        read(11,'(A)',iostat=io) Ax
        if(io /= 0) goto 99
      end do
    else
      m = m + 1
      pfa_spectrum(j)%SPS_NAME = Name
      pfa_spectrum(j)%NO_SPS = n_sps
      pfa_spectrum(j)%NO_LINES = nl
      pfa_spectrum(j)%SPS_SPECTAG = Spectag
      pfa_spectrum(j)%SPS_QLOG(1:3) = Qlog(1:3)
      do i = 1, nl
        read(11,*,iostat=io) (thbs(k),k=1,10)
        if(io /= 0) goto 99
        pfa_spectrum(j)%SPS_V0(i) = thbs(1)
        pfa_spectrum(j)%SPS_EL(i) = thbs(2)
        pfa_spectrum(j)%SPS_STR(i) = thbs(3)
        pfa_spectrum(j)%SPS_W(i) = thbs(4)
        pfa_spectrum(j)%SPS_PS(i) = thbs(5)
        pfa_spectrum(j)%SPS_N(i) = thbs(6)
        pfa_spectrum(j)%SPS_DELTA(i) = thbs(7)
        pfa_spectrum(j)%SPS_N1(i) = thbs(8)
        pfa_spectrum(j)%SPS_GAMMA(i) = thbs(9)
        pfa_spectrum(j)%SPS_N2(i) = thbs(10)
      end do
    endif
!
  END DO
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_spectro, mfi
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (spect_index(i),i=1,no_spectro)
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = mxco
  k = mfi + 2

  DEALLOCATE(s_phi_basis_copy,STAT=i)
  ALLOCATE(s_phi_basis_copy(k,no_spectro),STAT=io)
  IF(io /= 0) then
    Print *,'** ALLOCATE Error: s_phi_basis_copy, STAT =',io
    goto 99
  endif

  do j = 1, no_spectro
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    spectroscopic(j)%TYPE = AdjustL(Ax)
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    spectroscopic(j)%NAME = AdjustL(Ax)
    read(11,*,iostat=io) Spectag
    if(io /= 0) goto 99
    spectroscopic(j)%SPECTAG = Spectag
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    spectroscopic(j)%NO_PHI_VALUES = ht_i
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    spectroscopic(j)%NO_ZETA_VALUES = ht_i
    read(11,*,iostat=io) (spectroscopic(j)%DER_CALC(i),i=1,6)
    if(io /= 0) goto 99
    tau(1:k) = 0.0
    read(11,*,iostat=io) (tau(i),i=1,k)
    if(io /= 0) goto 99
    spectroscopic(j)%PHI_BASIS(1:k) = tau(1:k) * deg2rad
    s_phi_basis_copy(1:k,j) = spectroscopic(j)%PHI_BASIS(1:k)
    read(11,*,iostat=io) (spectroscopic(j)%ZETA_BASIS(i),i=1,kk+2)
    if(io /= 0) goto 99
  end do
!
! Create spect_atmos array:
!
  do k = 1, no_atmos
    Spectag = atmospheric(k)%Spectag
    do i = 1, no_spectro
      if(Spectag == spectroscopic(i)%Spectag) then
        spect_atmos(k) = i
        EXIT
      endif
    end do
  end do
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  InDir(1:)=' '
  read(11,'(A)',iostat=io) InDir
  if(io /= 0) goto 99
  InDir = AdjustL(InDir)

  Aaap(1:) = ' '
  Aaap = 'aaap.umls'

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  Ax(1:)=' '
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99
  Ax = AdjustL(Ax)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) si,n_lvls,no_tan_hts,mnz,no_filt_pts
  if(io /= 0) goto 99
!
! Set n_obs to be: N_lvls always.
!   (Changed, Aug/6/96 Z.Shippony & W.G.Read)

  n_obs = n_lvls

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) q, r, beta_inc
  if(io /= 0) goto 99

  href(1:n_lvls) = q
  zref(1:n_lvls) = r
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (p_indx(i),i=1,n_lvls)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (t_indx(i),i=1,no_tan_hts)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  thbs(1:) = 0.0
  read(11,*,iostat=io) (thbs(i),i=1,si-1)   ! tan_hts_below_surface
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (z_gnlv(i),i=1,mnz)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_t, no_phi_t
  if(io /= 0) goto 99

  if(no_mmaf < no_phi_t) then
    io = -1
    Print *,'** Error: no_mmaf < no_phi_t ...'
    Print *,'   Please correct input file and re-run !'
    goto 99
  endif

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (t_z_basis(i),i=1,no_t)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (t_phi_basis(i),i=1,no_phi_t)
  if(io /= 0) goto 99

  t_phi_basis(1:no_phi_t) = t_phi_basis(1:no_phi_t) * deg2rad
  t_phi_basis_copy(1:no_phi_t) = t_phi_basis(1:no_phi_t)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) ((t_coeff(i,j),j=1,no_phi_t),i=1,no_t)
  if(io /= 0) goto 99
!
! Complete the t_coeff array to have full NO_MMAF cover
!
  khi = no_phi_t + 1
  do
    klo = khi
    khi = klo + no_phi_t - 1
    if(khi > no_mmaf) EXIT
    t_coeff(1:no_t,klo:khi) = t_coeff(1:no_t,1:no_phi_t)
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (no_phi_f(i),i=1,no_atmos)
  if(io /= 0) goto 99

  do i = 1, no_atmos
    if(no_mmaf < no_phi_f(i)) then
      io = -1
      Print *,'** Error: no_mmaf < no_phi_f(i), i =',i
      Print *,'   Please correct input file and re-run !'
      goto 99
    endif
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (no_coeffs_f(i),i=1,no_atmos)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (is_f_log(i),i=1,no_atmos)     ! ** A a new one !
  if(io /= 0) goto 99
!
  tau(1:mnp) = 0.0
  DO m = 1, no_atmos
    kk = no_phi_f(m)
    ht_i = no_coeffs_f(m)
    read(11,*,iostat=io) (f_basis(i,m),i=1,ht_i)
    if(io /= 0) goto 99
    tau(1:kk) = 0.0
    read(11,*,iostat=io) (tau(i),i=1,kk)
    if(io /= 0) goto 99
    f_phi_basis(1:kk,m) = tau(1:kk) * deg2rad
    f_phi_basis_copy(1:kk,m) = f_phi_basis(1:kk,m)
    read(11,*,iostat=io) ((mr_f(i,j,m),j=1,kk),i=1,ht_i)
    if(io /= 0) goto 99
  END DO

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_freqs
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (freq_grid(i),i=1,no_freqs)
  if(io /= 0) goto 99

  CLOSE(11,iostat=i)
!
  Call Z_DATETIME(Dtm1)

  fft_pts = 10

  temp_der = .true.
! temp_der = .false.
!
! Get all the filter's loaded & define:
!
  Call get_filters(no_pfa_ch,no_filt_pts,pfa_ch,f_grid_filter, &
 &                 freqs,filter_func,InDir,ld,primag,ier)
  if(ier /= 0) goto 99
!
! Get the selected integration grid pressures. Also, define the GL
! pressure grid:

  z_grid(1:) = 0.0
  DO i = 1, n_lvls
    j = p_indx(i)
    z_grid(i) = z_gnlv(j)
  END DO
  z_grid(n_lvls+1) = z_grid(n_lvls)
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  phi_tan = phi_tan_mmaf(1)
  Call geoc_geod_conv(beta_inc,phi_tan,geoc_lat,E_rad)
!
! Compute the hydrostatic_model on the GL-Grid for all mmaf(s):
!
  Call hydrostatic_model(si, n_lvls, no_t, no_mmaf, t_indx,        &
       no_tan_hts, geoc_lat, Href, Zref, z_grid, thbs, t_z_basis,  &
       t_coeff, z_glgrid, h_glgrid, t_glgrid, dhdz_glgrid,         &
       dh_dt_glgrid, tan_press, tan_hts, tan_temp, tan_dh_dt,      &
       gl_count, Ier)
  IF(ier /= 0) goto 99
!
! Load the pointing vs. frequencies database for the given band
! (needed for frequency averaging)
!
  Fnd(1:) = ' '
! Fnd = '/home/zvi/ptg_frq_grid_bxx.dat'     ! HOME PC
  Fnd = '/user5/zvi/ptg_frq_grid_bxx.dat'    ! SUN, MLSGATE
  i = index(Fnd,'_bxx.dat')
  write(Fnd(i+2:i+3),'(i2.2)') band
!
  kk = -1
  no_ptg_frq(1:Nptg) = 0
!
  Close(32,iostat=i)
  Open(32,file=Fnd,action='READ',status='OLD',iostat=io)
  if(io /= 0) goto 44
!
! First entry in the file is the 'Band' frequency. All the rest are
! relative to this (center) frequency for this band
!
  Read(32,*,iostat=io) q
  if(io /= 0) goto 44
!
  DO
    Read(32,*,iostat=io) r, jp
    if(io > 0) goto 44
    if(io /= 0) EXIT
    Call Hunt(r,tan_press,no_tan_hts,k,i)
    IF(ABS(r-tan_press(i)) < ABS(r-tan_press(k))) k = i
    if(ABS(r-tan_press(k)) > 0.001) &
   &       Print *,'** Warning: Zeta:',r,' k:',k
    no_ptg_frq(k) = jp
    DEALLOCATE(ptg_frq_grid(k)%values,STAT=i)
    ALLOCATE(ptg_frq_grid(k)%values(jp),STAT=i)
    IF(i /= 0) THEN
      ier = i
      PRINT *,'** Error: ALLOCATION error for ptg_frq_grid ..'
      PRINT *,'   tan_hts index:',k,' STAT =',ier
      do l = 1, k
        DEALLOCATE(ptg_frq_grid(l)%values,STAT=i)
      end do
      goto 99
    ENDIF
    Read(32,*,iostat=io) (ptg_frq_grid(k)%values(i),i=1,jp)
    if(io /= 0) goto 44
    if(kk < 0) kk = k
!
! Add 'band' frequency to ptg_frq_grid to convert to absolute grid
!
    ptg_frq_grid(k)%values(1:jp) = ptg_frq_grid(k)%values(1:jp) + q
!
  END DO
!
  if(kk > 1) then
    jp = no_ptg_frq(kk)
    do k = 1, kk-1
      DEALLOCATE(ptg_frq_grid(k)%values,STAT=i)
      ALLOCATE(ptg_frq_grid(k)%values(jp),STAT=i)
      IF(i /= 0) THEN
        ier = i
        PRINT *,'** Error: ALLOCATION error for ptg_frq_grid ..'
        PRINT *,'   tan_hts index:',k,' STAT =',ier
        do l = 1, k
          DEALLOCATE(ptg_frq_grid(l)%values,STAT=i)
        end do
        goto 99
      ENDIF
      no_ptg_frq(k) = jp
      ptg_frq_grid(k)%values(1:jp) = ptg_frq_grid(kk)%values(1:jp)
    end do
  endif
!
 44 Close(32,iostat=i)
    if(io > 0) then
      ier = io
      goto 99
    else
      io = 0
    endif
!
! Compute all path entities for all mmafs and tanget pointings
!
  Call comp_path_entities(n_lvls,no_t,gl_count,ndx_path,z_glgrid,  &
       t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,atmospheric,     &
       no_atmos,f_basis,mr_f,no_coeffs_f,tan_hts,no_tan_hts,n_sps, &
       no_phi_f,f_phi_basis,z_path,h_path,t_path,phi_path,n_path,  &
       dhdz_path,dh_dt_path,no_phi_t,t_phi_basis,spsfunc_path,     &
       is_f_log,no_mmaf,phi_tan_mmaf,Ier)
  IF(ier /= 0) goto 99
!
! **********************  MAIN Mmaf Loop *******************

! DO l = 1, no_mmaf
  DO l = 1, 1                 ! ** DEBUG, only one mmaf
!
    phi_tan = phi_tan_mmaf(l)
!
    t_phi_basis(1:no_phi_t) = t_phi_basis_copy(1:no_phi_t) + phi_tan

    DO j = 1, no_atmos
      k = no_phi_f(j)
      f_phi_basis(1:k,j) = f_phi_basis_copy(1:k,j) + phi_tan
    end do

    k = mfi + 2
    do j = 1, no_spectro
      spectroscopic(j)%PHI_BASIS(1:k) = s_phi_basis_copy(1:k,j) + phi_tan
    end do
!
! Compute the ptg_angles (chi) for Antenna convolution, also the derivatives
! of chi w.r.t to T and other parameters
!
    Call get_chi_angles(ndx_path(1:,l),n_path(1:,l),tan_press,tan_hts(1:,l),&
   &     tan_temp(1:,l),phi_tan,RoC,h_obs,elev_offset,tan_dh_dt(1:,l,1:),   &
   &     n_lvls,no_tan_hts,no_t,t_z_basis,z_grid,si,center_angle,           &
   &     ptg_angles(1:,l),dx_dt,d2x_dxdt,ier)
    IF(ier /= 0) goto 99

! Compute the refraction correction scaling matrix for this mmaf:
!
    Call refraction_correction(no_tan_hts, tan_hts(1:,l), h_path(1:,l), &
   &                n_path(1:,l), ndx_path(1:,l), E_rad, ref_corr)
!
! *** DEBUG
    klo = -1
    zeta = -1.666667
    Call Hunt(zeta,tan_press,no_tan_hts,klo,khi)
    IF(ABS(zeta-tan_press(khi)) < ABS(zeta-tan_press(klo))) klo = khi
! *** END DEBUG
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
      kk = no_ptg_frq(k)
!
      if(kk /= prev_npf) then
!
        prev_npf = kk
        f_grid(1:mnf) = 0.0
        DEALLOCATE(k_temp_frq%values,STAT=i)
!
        do j = 1, n_sps
          DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
          DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
        end do
!
        ALLOCATE(k_temp_frq%values(kk,no_t,no_phi_t),STAT=ier)
        IF(ier /= 0) then
          Print *,'** ALLOCATE Error: k_temp_frq array, STAT =',ier
          goto 99
        endif
!
        do j = 1, n_sps
          m = max(1,no_phi_f(j))
          i = max(1,no_coeffs_f(j))
          ALLOCATE(k_atmos_frq(j)%values(kk,i,m),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: k_atmos_frq, STAT =',ier
            goto 99
          endif
        end do
!
        do m = 1, n_sps
          j = spect_atmos(m)
          if(j < 1) CYCLE
          if(.not.spectroscopic(j)%DER_CALC(band)) CYCLE
          Vname(1:) = ' '
          Spectag = spectroscopic(j)%Spectag
          DO
            if(spectroscopic(j)%Spectag /= Spectag) EXIT
            n = spectroscopic(j)%no_phi_values
            i = spectroscopic(j)%no_zeta_values
            CA = spectroscopic(j)%type
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
            if(j > 3 * n_sps) EXIT
          END DO
        end do

      endif            ! On DEALLOCATE/ALLOCATE cycle
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
      khi = ndx_path(ptg_i,l)%total_number_of_elements
      Call get_beta_path(ptg_i,pfa_spectrum,khi,no_ptg_frq, &
     &     ptg_frq_grid,z_path(ptg_i,l),t_path(ptg_i,l),beta_path,ier)
      IF(ier /= 0) goto 99
!
      k_temp_frq%values = 0.0
      do j = 1, n_sps
        k_atmos_frq(j)%values = 0.0
        k_spect_dw_frq(j)%values = 0.0
        k_spect_dn_frq(j)%values = 0.0
        k_spect_dnu_frq(j)%values = 0.0
      end do
!
      RadV(1:mnf) = 0.0
      f_grid(1:kk) = ptg_frq_grid(k)%values(1:kk)

      do frq_i = 1, kk
!
        Frq = f_grid(frq_i)
!
        Call Rad_Tran(Frq, N_lvls, h_tan, n_sps, ndx_path(k,l),  &
       &    z_path(k,l), h_path(k,l), t_path(k,l), phi_path(k,l),&
       &    dHdz_path(k,l), earth_ref, beta_path(1:,frq_i),      &
       &    spsfunc_path(1:,k,l), ref_corr(1:,k), s_temp, brkpt, &
       &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
        IF(ier /= 0) goto 99
!
        RadV(frq_i) = Rad
!
! Now, Compute the radiances derivatives:
!
        CALL Rad_Tran_WD(frq_i,band,Frq,N_lvls,n_sps,z_path(k,l),       &
       &     h_path(k,l),t_path(k,l),phi_path(k,l),dHdz_path(k,l),      &
       &     atmospheric,beta_path(1:,frq_i),spsfunc_path(1:,k,l),      &
       &     t_z_basis,f_basis,no_coeffs_f,mr_f,no_t,ref_corr(1:,k),    &
       &     no_phi_f,f_phi_basis,temp_der,no_phi_t,t_phi_basis,        &
       &     dh_dt_path(k,l),spect_atmos,spectroscopic,k_temp_frq,      &
       &     k_atmos_frq,k_spect_dw_frq,k_spect_dn_frq,k_spect_dnu_frq, &
       &     is_f_log,brkpt,no_ele,mid,ilo,ihi,t_script,tau,ier)
        IF(ier /= 0) goto 99
!
      end do

! Frequency Average the radiances with the appropriate filter shapes
!
      do i = 1, no_pfa_ch
        ch = pfa_ch(i)
        Call Freq_Avg(f_grid,F_grid_filter(1:,i),Filter_func(1:,i), &
       &              RadV,kk,maxfiltpts,Rad,Ier)
        IF(ier /= 0) goto 99
        Radiances(ptg_i,ch) = Rad
      end do

!     if(i > -3) goto 77           ! ** DEBUG, Bypass derivatives avg.
!
! Frequency Average the temperature derivatives with the appropriate
! filter shapes
!
      if(temp_der) then
        RadV(1:mnf) = 0.0
        do i = 1, no_pfa_ch
!         ch = pfa_ch(i)
          ch = i               ! ** DEBUG, memory limitations on MLSGATE
          do j = 1, no_phi_t
            do k = 1, no_t
              RadV(1:kk) = k_temp_frq%values(1:kk,k,j)
              Call Freq_Avg(f_grid,F_grid_filter(1:,i),Filter_func(1:,i), &
       &           RadV,kk,maxfiltpts,r,Ier)
              IF(ier /= 0) goto 99
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
        do j = 1, n_sps
          if(atmospheric(j)%der_calc(band)) THEN
            RadV(1:mnf) = 0.0
            do k = 1, no_phi_f(j)
              do n = 1, no_coeffs_f(j)
                RadV(1:kk) = k_atmos_frq(j)%values(1:kk,n,k)
                Call Freq_Avg(f_grid,F_grid_filter(1:,i),Filter_func(1:,i), &
       &             RadV,kk,maxfiltpts,r,Ier)
                IF(ier /= 0) goto 99
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
        do m = 1, n_sps
          j = spect_atmos(m)
          if(.not.spectroscopic(j)%DER_CALC(band)) CYCLE
          Spectag = spectroscopic(j)%Spectag
          DO
            if(spectroscopic(j)%Spectag /= Spectag) EXIT
            RadV(1:mnf) = 0.0
            CA = spectroscopic(j)%type
            do k = 1, spectroscopic(j)%no_phi_values
              do n = 1, spectroscopic(j)%no_zeta_values
                select case ( CA )
                  case ( 'W' )
                    RadV(1:kk) = k_spect_dw_frq(m)%values(1:kk,n,k)
                  case ( 'N' )
                    RadV(1:kk) = k_spect_dn_frq(m)%values(1:kk,n,k)
                  case ( 'V' )
                    RadV(1:kk) = k_spect_dnu_frq(m)%values(1:kk,n,k)
                end select
                Call Freq_Avg(f_grid,F_grid_filter(1:,i),Filter_func(1:,i),&
               &              RadV,kk,maxfiltpts,r,Ier)
                IF(ier /= 0) goto 99
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
            if(j > 3 * n_sps) EXIT
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
      k_temp(i,kk,1:no_t,1:no_phi_t)=k_temp(i,kk-1,1:no_t,1:no_phi_t)
      do m = 1, n_sps
        if(atmospheric(m)%der_calc(band)) then
          k = no_phi_f(m)
          n = no_coeffs_f(m)
          k_atmos(i,kk,1:n,1:k,m)=k_atmos(i,kk-1,1:n,1:k,m)
        endif
      end do
      do m = 1, n_sps
        j = spect_atmos(m)
        if(.not.spectroscopic(j)%DER_CALC(band)) CYCLE
        Spectag =  spectroscopic(j)%Spectag
        DO
          if(spectroscopic(j)%Spectag /= Spectag) EXIT
          k = spectroscopic(j)%no_phi_values
          n = spectroscopic(j)%no_zeta_values
          k_spect_dw(i,kk,1:n,1:k,j)=k_spect_dw(i,kk-1,1:n,1:k,j)
          k_spect_dn(i,kk,1:n,1:k,j)=k_spect_dn(i,kk-1,1:n,1:k,j)
          k_spect_dnu(i,kk,1:n,1:k,j)=k_spect_dnu(i,kk-1,1:n,1:k,j)
          j = j + 1
          if(j > 3 * n_sps) EXIT
        END DO
      end do
    end do
!
!  Here comes the Convolution code
!
    DO i = 1, no_pfa_ch
      ch = pfa_ch(i)
      Call convolve_all(ch, ptg_press, atmospheric, n_sps, temp_der,   &
     &     tan_press,ptg_angles(1:,l),tan_temp(1:,l), dx_dt, d2x_dxdt, &
     &     band,center_angle,fft_pts,Radiances(1:,ch),k_temp(i,1:,1:,1:), &
     &     k_atmos(i,1:,1:,1:,1:), k_spect_dw(i,1:,1:,1:,1:),          &
     &     k_spect_dn(i,1:,1:,1:,1:),k_spect_dnu(i,1:,1:,1:,1:),       &
     &     spect_atmos,no_tan_hts,k_info_count,i_star_all,k_star_all,  &
     &     k_star_info,no_t,no_phi_t,no_phi_f,InDir,Aaap,spectroscopic,&
     &     Ier)
    END DO

  END DO                ! Mmaf Loop
!
  DEALLOCATE(k_temp_frq%values,STAT=i)
  do j = 1, n_sps
    DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
    DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
  end do
!
! *** DEBUG Print
!
    kk = no_tan_hts - si + 1
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      write(*,903) ch,char(92),kk
      write(*,905) (Radiances(k,ch),k=si,no_tan_hts)
    end do
903 format('ch',i2.2,'_avg_pfa_rad',a1,i2.2)
905 format(4(2x,1pg15.8))
!
!   if(kk > -5) goto 99      ! Bypass print of ALL derivatives
!
     frq_i = 1
     Spectag = 18003
!    ch = pfa_ch(1)
     ch = 1                  ! ** DEBUG, memory limitations on MLSGATE
!
     m = -1
     do j = 1, n_sps
       if(atmospheric(j)%SPECTAG == Spectag) then
         m = j
         EXIT
       endif
     end do

!   if(m < -100) then           ! Bypass print of k_atmos derivatives
    if(m > 0) then
      Ax(1:)=' '
      k = no_coeffs_f(m)
      Ax = 'avg_k_'//atmospheric(m)%NAME
      l = Len_Trim(Ax) + 1
      Ax = Ax(1:l-1)//'_'
      Call StrLwr(Ax)
      j = (no_phi_f(m)+1)/2
      Call Hunt(zeta,f_basis(1:,m),k,kz,i)
      IF(ABS(zeta-f_basis(i,m)) < ABS(zeta-f_basis(kz,m))) kz = i
      Print 900,Ax(1:l),frq_i,frq_i
      Print 901, k_atmos(ch,klo,kz,j,m)
    endif
!
!   if(kk > -5) goto 99            ! Bypass print of k_temp derivatives
!
    j = (no_phi_t+1)/2
    Call Hunt(zeta,t_z_basis,no_t,kz,i)
    IF(ABS(zeta-t_z_basis(i)) < ABS(zeta-t_z_basis(kz))) kz = i
    Print 900,'avg_k_temp_',frq_i,frq_i
    Print 901, k_temp(ch,klo,kz,j)
!
!   if(kk > -5) goto 99            ! Bypass print of k_spect derivatives
!
!   ch = pfa_ch(1)
    ch = 1                  ! ** DEBUG, memory limitations on MLSGATE
!
    do m = 1, n_sps
      j = spect_atmos(m)
      if(.not.spectroscopic(j)%DER_CALC(band)) CYCLE
      Spectag =  spectroscopic(j)%Spectag
      DO
        if(spectroscopic(j)%Spectag /= Spectag) EXIT
        CA = spectroscopic(j)%type
        Ax(1:) = ' '
        RadV(1:mnf) = 0.0
        kk = spectroscopic(j)%no_phi_values
        ht_i = spectroscopic(j)%no_zeta_values
        select case ( CA )
          case ( 'W' )
            Ax = 'avg_k_spect_dw_'
            r = SUM(k_spect_dw(ch,klo,1:ht_i,1:kk,j))
          case ( 'N' )
            Ax = 'avg_k_spect_dn_'
            r = SUM(k_spect_dn(ch,klo,1:ht_i,1:kk,j))
          case ( 'V' )
            Ax = 'avg_k_spect_dnu_'
            r = SUM(k_spect_dnu(ch,klo,1:ht_i,1:kk,j))
        end select
        i = Len_Trim(Ax)
        Print 900,Ax(1:i),frq_i,ht_i*kk
        do k = 1, ht_i
          select case ( CA )
            case ( 'W' )
              Print 901,(k_spect_dw(ch,klo,k,i,j),i=1,kk)
            case ( 'N' )
              Print 901,(k_spect_dn(ch,klo,k,i,j),i=1,kk)
            case ( 'V' )
              Print 901,(k_spect_dnu(ch,klo,k,i,j),i=1,kk)
          end select
        end do
        Print *,'  Sum over all zeta & phi coeff:',sngl(r)
        j = j + 1
        if(j > 3 * n_sps) EXIT
      END DO
    end do
!
900 format(a,i2.2,'\',i3.3)
901 format(5(2x,1pe13.6))
!

 99  CLOSE(11,iostat=i)
     if(io /= 0) Call ErrMsg(' ',io)

     Call Z_DATETIME(Dtm2)
     Print *
     Print *,'** Actual computations started at: ',Dtm1
     Print *,'**            Program finished at: ',Dtm2

     Ax(1:40) = ' '
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
!---------------------------------------------------------------------------

SUBROUTINE radiometry(ch, f_p, f_i, db_fi, lmt)

! This subroutine calculates the center frequency of the primary and image
! sideband by channel. It also Returns the bandwidth limits of integration
! and the gain of the primary sideband relative to the image in db units.

INTEGER(i4), INTENT(IN) :: ch
INTEGER(i4), INTENT(OUT) :: lmt

REAL(r8), INTENT(OUT) :: f_p
REAL(r8), INTENT(OUT) :: f_i
REAL(r8), INTENT(OUT) :: db_fi

LOGICAL, SAVE :: sgn_fp(6) = (/                                    &
                 .false., .true., .false., .true., .true., .false./)

INTEGER(i4) :: band, sub_ch, j

Real(r4), SAVE :: db_fi_data(90) =(/                   &
     &    -0.5218,  0.0000,  0.5218,  0.8264,  0.9654, &
     &     1.0332,  1.0668,  1.0890,  1.1111,  1.1440, &
     &     1.2091,  1.3365,  1.5807,  0.3476, -3.6798, &
     &    -1.7854, -1.6204, -1.5519, -1.5142, -1.4966, &
     &    -1.4884, -1.4840, -1.4811, -1.4783, -1.4741, &
     &    -1.4659, -1.4486, -1.4167, -1.3607, -1.2602, &
     &    -1.0409, -1.1206, -1.1641, -1.1874, -1.1995, &
     &    -1.2052, -1.2084, -1.2104, -1.2124, -1.2157, &
     &    -1.2216, -1.2347, -1.2607, -1.3128, -1.4233, &
     &    -1.0921, -1.2329, -1.2863, -1.3121, -1.3234, &
     &    -1.3284, -1.3309, -1.3327, -1.3343, -1.3367, &
     &    -1.3416, -1.3507, -1.3667, -1.3884, -1.4035, &
     &     0.6661,  0.6480,  0.7759,  0.8029,  0.8741, &
     &     0.8307,  0.9736,  0.9273,  0.5800,  0.8906, &
     &     0.9772,  0.9596,  0.9426,  0.8669,  0.7433, &
     &     0.3108,  1.0632,  0.7029,  0.4426,  0.2910, &
     &     0.2261,  0.2396,  0.2244,  0.2237,  0.1811, &
     &     0.1795,  0.1534,  0.1418,  0.3651,  0.5647 /)
!

REAL(r8), SAVE :: f_prime(6) = (/                          &
    63568.418D0, 204352.161D0, 204574.627D0, 206132.067D0, &
    183310.062D0, 184377.788D0/)

REAL(r8), SAVE :: f_image(6) = (/                          &
    62997.812D0, 202181.555D0, 201959.089D0, 200401.648D0, &
    186245.513D0, 185177.788D0/)

REAL(r8) :: ch_offset

  band = (ch - 1) / 15 + 1
  sub_ch = ch - 15 * (band - 1)
!
  IF(sub_ch == 8) THEN
    j = 0
    lmt = 1
  ELSE                      ! Above and below the spectral center
    lmt = 2**(ABS(sub_ch - 8) - 1)
    j = SIGN(3*lmt - 1, sub_ch - 8)
  END IF

  ch_offset = float(j)
  db_fi = db_fi_data(ch)

  IF(sgn_fp(band)) THEN
    f_p = f_prime(band) + ch_offset
    f_i = f_image(band) - ch_offset
  ELSE
    f_p = f_prime(band) - ch_offset
    f_i = f_image(band) + ch_offset
  END IF

  Return
END SUBROUTINE radiometry

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
