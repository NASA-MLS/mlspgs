module L2_LOAD_M
  use GL6P, only: NG
  use MLSCommon, only: I4, R4, R8
  use L2_TEST_STRUCTURES_M
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT, &
                                  DEG2RAD
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, PFA_SLAB, SPECTRO_PARAM, &
                                 LIMB_PRESS
  use L2PCdim, only: N2lvl, Nptg, NCH, MNP => max_no_phi
  use PATH_ENTITIES_M, only: PATH_VECTOR
  use D_HUNT_M, only: HUNT          ! ** DEBUG
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
! ===========================================     L2_LOAD =====
! This subprogram loads all the needed inputs for the l2_test run
!
  SUBROUTINE L2_LOAD(FMC, FMI, T_FMI, Ier)
!---------------------------------------------------------------------------

Type(fwd_mdl_config), INTENT(IN OUT) :: FMC

Type(fwd_mdl_info), OPTIONAL, INTENT(OUT) :: FMI
Type(temporary_fwd_mdl_info), OPTIONAL, INTENT(OUT) :: T_FMI

Integer(i4), INTENT(OUT) :: ier
!
! ---- Local variables ---------
!
Logical :: geom_deriv(6)
Integer(i4) :: i, j, k, kk, ht_i, no_t, mnz, kz, si, n_sps, io, nl, &
               Spectag, m, no_phi_t, no_mmaf, jj

Integer(i4) :: ch1, ch2, no_pfa_ch, pfa_ch(2)

Real(r8) :: dummy(N2lvl), thbs(10), Qlog(3)

real(r8) :: freqs(Nch)

Real(r8) :: z0, zn, q, r, v

!
Character (LEN=08) :: Name
Character (LEN=40) :: Ax
Character (LEN=80) :: Fnd, Line

!  ----------------------
! Read convolution & freq. averaging data from file: tmp.dat
!
  Line = 'tmp.dat'
  CLOSE(13,iostat=i)
  OPEN(13,file=Line,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99
!
  FMC%Zfrq = -1.0
  read(13,*,iostat=io) FMC%do_conv
  if(io /= 0) goto 99
  read(13,*,iostat=io) FMC%do_frqavg
  if(io /= 0) goto 99
  if(.not. FMC%do_frqavg) then
    read(13,*,iostat=io) FMC%Zfrq
    if(io /= 0) goto 99
  endif
  CLOSE(13,iostat=io)
!  ----------------------
! Read the rest of the inputs from a file...
!
  ier = 0
  Fnd(1:) = ' '
  Fnd = FMC%Z
  Line = Fnd

  CLOSE(11,iostat=io)
  OPEN(11,file=Fnd,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99
!
! *** Begining the FMC Section
!
  do
    Ax(1:) = ' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    if (Index(Ax,'No_Mmaf') > 0) EXIT
  end do

  read(11,*,iostat=io) j
  if(io /= 0) goto 99

  FMC%no_mmaf = min(j,100)   ! Assuming mmaf chunk of no more then 100
  no_mmaf = FMC%no_mmaf

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  dummy(1:no_mmaf) = 0.0
  read(11,*,iostat=io) (dummy(i),i=1,no_mmaf)
  if(io /= 0) goto 99

  DEALLOCATE(FMC%phi_tan_mmaf,STAT=i)
  ALLOCATE(FMC%phi_tan_mmaf(no_mmaf),STAT=i)
  if(i /= 0) goto 99

  FMC%phi_tan_mmaf(1:no_mmaf) = dummy(1:no_mmaf) * deg2rad

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) ch1, ch2
  if(io /= 0) goto 99

  no_pfa_ch = min(2,ch2-ch1+1)
  do i = 1, no_pfa_ch
    ch2 = ch1 + i - 1
    pfa_ch(i) = ch2
  end do

  FMC%Channels_range(1) = ch1
  FMC%Channels_range(2) = ch2

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMC%Sideband
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMC%temp_der, FMC%atmos_der, FMC%spect_der
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) k, j
  if(io /= 0) goto 99

  FMC%n_lvls = k
  FMC%no_tan_hts = j

  DEALLOCATE(FMC%p_Indx,FMC%t_Indx,STAT=i)
  ALLOCATE(FMC%p_Indx(k),FMC%t_Indx(j),STAT=i)
  if(i /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (FMC%p_indx(i),i=1,k)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (FMC%t_indx(i),i=1,FMC%no_tan_hts)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  Fnd(1:)=' '
  read(11,'(A)',iostat=io) Fnd
  if(io /= 0) goto 99

  Fnd = AdjustL(Fnd)
  i = LEN_TRIM(Fnd)
  if(Fnd(i:i) /= '/') then
    i = i + 1
    Fnd(i:i) = '/'
  endif

  FMC%InDir = Fnd

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  Fnd(1:)=' '
  read(11,'(A)',iostat=io) Fnd
  if(io /= 0) goto 99

  FMC%B = AdjustL(Fnd)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  Fnd(1:)=' '
  read(11,'(A)',iostat=io) Fnd
  if(io /= 0) goto 99

  FMC%aaap_file = AdjustL(Fnd)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99
!
! *** Begining the FMI Section
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMI%n_sps
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  n_sps = FMI%n_sps
  DEALLOCATE(FMI%Species,STAT=i)
  ALLOCATE(FMI%Species(n_sps),STAT=i)
  if(i /= 0) goto 99

  do i = 1, FMI%n_sps
    read(11,*,iostat=io) FMI%Species(i)
    if(io /= 0) goto 99
  end do
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMI%band
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMI%Surface_index,FMI%max_no_zeta, &
                       FMI%no_filt_pts,FMI%fft_pts
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) z0, zn
  if(io /= 0) goto 99

  k = FMC%n_lvls
  mnz = FMI%Max_no_zeta
  si = FMI%Surface_index

  r = (zn - z0)/ (mnz - 1)

  j = k + si
  DEALLOCATE(FMI%z_grid,FMI%tan_press,STAT=i)
  ALLOCATE(FMI%z_grid(j),FMI%tan_press(j),STAT=i)
  if(i /= 0) goto 99
!
  FMI%z_grid(j) = 0.0
  DO i = 1, k
    j = FMC%p_indx(i)
    FMI%z_grid(i) = z0 + (j - 1) * r
  END DO
  FMI%z_grid(k+1) = FMI%z_grid(k)
!
! Define tan_press as a TRUE subset of z_grid:
!
  kz = si - 1
  FMI%tan_press(1:kz) = -5.0
  do i = 1, FMC%no_tan_hts
    kz = kz + 1
    j = FMC%t_indx(i)
    FMI%tan_press(kz) = FMI%z_grid(j)
  end do
  FMI%no_ptg_frq(1) = kz
!
  j = 2**FMI%fft_pts
  DEALLOCATE(FMI%AAAP,FMI%D1AAAP,FMI%D2AAAP,STAT=i)
  ALLOCATE(FMI%AAAP(j,3),FMI%D1AAAP(j,3),FMI%D2AAAP(j,3),STAT=i)
  if(i /= 0) goto 99

  Call ANTENNA(FMC%aaap_file, FMI%fft_pts, FMI%Xlamda, FMI%AAAP, &
 &             FMI%D1AAAP, FMI%D2AAAP, FMI%IAS, Ier)
  if(ier /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  DEALLOCATE(FMI%Tan_hts_below_surface,STAT=i)
  ALLOCATE(FMI%Tan_hts_below_surface(si),STAT=i)
  if(i /= 0) goto 99

  FMI%Tan_hts_below_surface(1:si) = 0.0
  read(11,*,iostat=io) (FMI%Tan_hts_below_surface(i),i=1,si-1)
  if(io /= 0) goto 99

  n_sps = FMI%n_sps
  DEALLOCATE(FMI%Pfa_spectrum,STAT=i)
  ALLOCATE(FMI%Pfa_spectrum(n_sps),STAT=i)
  if(i /= 0) goto 99

  m = 0
  FMI%Pfa_spectrum(1)%NO_SPS = n_sps   ! Make sure we have this

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
!
    Read(Line,*,iostat=io) Name, Spectag, nl, (Qlog(i),i=1,3)
    if(io /= 0) goto 99
!
    j = 0
    DO i = 1, n_sps
      if(Name == FMI%Species(i)) then
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
      FMI%Pfa_spectrum(j)%SPS_NAME = Name
      FMI%Pfa_spectrum(j)%NO_SPS = n_sps
      FMI%Pfa_spectrum(j)%NO_LINES = nl
      FMI%Pfa_spectrum(j)%SPS_SPECTAG = Spectag
      FMI%Pfa_spectrum(j)%SPS_QLOG(1:3) = Qlog(1:3)
      do i = 1, nl
        read(11,*,iostat=io) (thbs(k),k=1,10)
        if(io /= 0) goto 99
        FMI%Pfa_spectrum(j)%SPS_V0(i) = thbs(1)
        FMI%Pfa_spectrum(j)%SPS_EL(i) = thbs(2)
        FMI%Pfa_spectrum(j)%SPS_STR(i) = thbs(3)
        FMI%Pfa_spectrum(j)%SPS_W(i) = thbs(4)
        FMI%Pfa_spectrum(j)%SPS_PS(i) = 0.0   !  thbs(5) *** DEBUG
        FMI%Pfa_spectrum(j)%SPS_N(i) = thbs(6)
        FMI%Pfa_spectrum(j)%SPS_DELTA(i) = thbs(7)
        FMI%Pfa_spectrum(j)%SPS_N1(i) = thbs(8)
        FMI%Pfa_spectrum(j)%SPS_GAMMA(i) = thbs(9)
        FMI%Pfa_spectrum(j)%SPS_N2(i) = thbs(10)
      end do
    endif
!
  END DO
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) FMI%no_spectro, FMI%mfi
  if(io /= 0) goto 99

  j = FMI%no_spectro
  DEALLOCATE(FMI%spectroscopic,STAT=i)
  ALLOCATE(FMI%spectroscopic(j),STAT=i)
  if(i /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = mxco
  k = FMI%mfi + 2

  Line(1:) = ' '
  do j = 1, FMI%no_spectro
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%TYPE = AdjustL(Ax)
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%NAME = AdjustL(Ax)
    read(11,*,iostat=io) Spectag
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%SPECTAG = Spectag
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%NO_PHI_VALUES = ht_i
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%NO_ZETA_VALUES = ht_i
    read(11,*,iostat=io) (FMI%spectroscopic(j)%DER_CALC(i),i=1,6)
    if(io /= 0) goto 99
    dummy(1:k) = 0.0
    read(11,*,iostat=io) (dummy(i),i=1,k)
    if(io /= 0) goto 99
    FMI%spectroscopic(j)%PHI_BASIS(1:k) = dummy(1:k) * deg2rad
    read(11,*,iostat=io) (FMI%spectroscopic(j)%ZETA_BASIS(i),i=1,kk+2)
    if(io /= 0) goto 99
  end do
!
! *** Begining the T_FMI Section
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) T_FMI%ptg_press%name, &
                      (T_FMI%ptg_press%der_calc(i),i=1,6)
  if(io /= 0) goto 99

  read(11,*,iostat=io) j
  if(io /= 0) goto 99

  T_FMI%ptg_press%no_lin_values = j
  read(11,*,iostat=io) (T_FMI%ptg_press%lin_val(i),i=1,j)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) T_FMI%No_Geometric
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  do i = 1, T_FMI%No_Geometric
    read(11,*,iostat=io) Name, geom_deriv, r
    if(io /= 0) goto 99
    IF(Name == 'ELEV_183') THEN
      T_FMI%elev_183 = r
    ELSE IF(Name == 'ELEV_205') THEN
      T_FMI%elev_205 = r
    ELSE IF(Name == 'EARTHREF') THEN
      T_FMI%earth_ref = r
    ELSE IF(Name == 'SPACE_T') THEN
      T_FMI%s_temp = r
    ELSE IF(Name == 'GEOCSRAD') THEN
      T_FMI%h_obs = r
    END IF
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  DEALLOCATE(T_FMI%Atmospheric,STAT=i)
  ALLOCATE(T_FMI%Atmospheric(n_sps),STAT=i)
  if(i /= 0) goto 99

  kk = mxco
  do k = 1, n_sps
    j = -1
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    Name = AdjustL(Ax)
    do jj = 1, n_sps
      if(Name == FMI%Species(jj)) then
        j = jj
        EXIT
      endif
    end do
    if(j < 0) then
      io = 5
      Print *,'** Error in: l2_load ..'
      Print *,'   Unknown molecule: ',Name,' in T_FMI%Atmospheric input ..'
      goto 99
    endif
    T_FMI%Atmospheric(j)%NAME = Name
    read(11,*,iostat=io) Spectag, ht_i
    if(io /= 0) goto 99
    T_FMI%Atmospheric(j)%SPECTAG = Spectag
    T_FMI%Atmospheric(j)%NO_LIN_VALUES = ht_i
    read(11,*,iostat=io) (T_FMI%Atmospheric(j)%FWD_CALC(i),i=1,6)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (T_FMI%Atmospheric(j)%DER_CALC(i),i=1,6)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (T_FMI%Atmospheric(j)%LIN_VAL(i),i=1,kk)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (T_FMI%Atmospheric(j)%BASIS_PEAKS(i),i=1,kk+2)
    if(io /= 0) goto 99
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) T_FMI%Zref, T_FMI%beta_inc
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  DEALLOCATE(T_FMI%Href,STAT=i)
  ALLOCATE(T_FMI%Href(no_mmaf),STAT=i)

  read(11,*,iostat=io) (T_FMI%Href(i),i=1,no_mmaf)
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99
!
  read(11,*,iostat=io) T_FMI%no_t, T_FMI%no_phi_t
  if(io /= 0) goto 99

  if(FMC%no_mmaf < T_FMI%no_phi_t) then
    io = -1
    Print *,'** Error: no_mmaf < no_phi_t ...'
    Print *,'   Please correct input file and re-run !'
    goto 99
  endif

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  no_t = T_FMI%no_t
  no_phi_t = T_FMI%no_phi_t
  DEALLOCATE(T_FMI%t_zeta_basis,T_FMI%t_phi_basis,T_FMI%t_coeff, &
 &           T_FMI%t_phi_basis_copy,STAT=i)
  ALLOCATE(T_FMI%t_zeta_basis(no_t),T_FMI%t_phi_basis(no_phi_t),   &
 &    T_FMI%t_phi_basis_copy(no_phi_t),T_FMI%t_coeff(no_t,no_mmaf),&
      STAT=i)
  if(i /= 0) goto 99

  read(11,*,iostat=io) (T_FMI%t_zeta_basis(i),i=1,no_t)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  dummy(1:) = 0.0
  read(11,*,iostat=io) (dummy(i),i=1,no_phi_t)
  if(io /= 0) goto 99

  T_FMI%t_phi_basis(1:no_phi_t) = dummy(1:no_phi_t) * deg2rad
  T_FMI%t_phi_basis_copy(1:no_phi_t) = T_FMI%t_phi_basis(1:no_phi_t)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) ((T_FMI%t_coeff(i,j),j=1,no_mmaf),i=1,no_t)
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  DEALLOCATE(T_FMI%no_phi_f,T_FMI%no_coeffs_f,T_FMI%is_f_log, &
 &           STAT=i)
  ALLOCATE(T_FMI%no_phi_f(n_sps),T_FMI%no_coeffs_f(n_sps), &
 &         T_FMI%is_f_log(n_sps),STAT=i)
  if(i /= 0) goto 99

  read(11,*,iostat=io) (T_FMI%no_phi_f(i),i=1,n_sps)
  if(io /= 0) goto 99

  do i = 1, n_sps
    if(no_mmaf < T_FMI%no_phi_f(i)) then
      io = -1
      Print *,'** Error: no_mmaf < no_phi_f(i), i =',i
      Print *,'   Please correct input file and re-run !'
      goto 99
    endif
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (T_FMI%no_coeffs_f(i),i=1,n_sps)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (T_FMI%is_f_log(i),i=1,n_sps)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = MAXVAL(T_FMI%no_phi_f)
  ht_i = MAXVAL(T_FMI%no_coeffs_f)
  DEALLOCATE(T_FMI%f_zeta_basis,T_FMI%f_phi_basis,T_FMI%mr_f, &
             T_FMI%f_phi_basis_copy,STAT=i)
  ALLOCATE(T_FMI%f_zeta_basis(ht_i,n_sps),   &
 &         T_FMI%f_phi_basis(kk,n_sps),      &
 &         T_FMI%f_phi_basis_copy(kk,n_sps), &
 &         T_FMI%mr_f(ht_i,kk,n_sps), STAT=i)
  if(i /= 0) goto 99
!
  DO m = 1, n_sps
    kk = T_FMI%no_phi_f(m)
    ht_i = T_FMI%no_coeffs_f(m)
    read(11,*,iostat=io) (T_FMI%f_zeta_basis(i,m),i=1,ht_i)
    if(io /= 0) goto 99
    dummy(1:kk) = 0.0
    read(11,*,iostat=io) (dummy(i),i=1,kk)
    if(io /= 0) goto 99
    T_FMI%f_phi_basis(1:kk,m) = dummy(1:kk) * deg2rad
    T_FMI%f_phi_basis_copy(1:kk,m) = T_FMI%f_phi_basis(1:kk,m)
    read(11,*,iostat=io) ((T_FMI%mr_f(i,j,m),j=1,kk),i=1,ht_i)
    if(io /= 0) goto 99
  END DO
!
! Create spect_atmos array:
!
  DEALLOCATE(FMI%spect_atmos,STAT=i)
  ALLOCATE(FMI%spect_atmos(n_sps),STAT=i)
  if(i /= 0) goto 99

  do k = 1, n_sps
    Spectag = T_FMI%atmospheric(k)%Spectag
    do i = 1, FMI%no_spectro
      if(Spectag == FMI%spectroscopic(i)%Spectag) then
        FMI%spect_atmos(k) = i
        EXIT
      endif
    end do
  end do
!
  kk = mxco
  k = FMI%mfi + 2

  DEALLOCATE(T_FMI%S_PHI_BASIS_COPY,STAT=i)
  ALLOCATE(T_FMI%S_PHI_BASIS_COPY(k,FMI%no_spectro),STAT=io)
  IF(io /= 0) then
    Print *,'** ALLOCATE Error: T_FMI%S_PHI_BASIS_COPY, STAT =',io
    goto 99
  endif

  do j = 1, FMI%no_spectro
    T_FMI%S_PHI_BASIS_COPY(1:k,j) = FMI%spectroscopic(j)%PHI_BASIS(1:k)
  end do

  CLOSE(11,iostat=i)
!
! Get all the filter's loaded & define:
!
  j = no_pfa_ch
  k = FMI%no_filt_pts
  DEALLOCATE(FMI%F_grid_filter,FMI%Filter_func,STAT=i)
  ALLOCATE(FMI%F_grid_filter(k,j),FMI%Filter_func(k,j),STAT=i)
  if(i /= 0) goto 99

  freqs(1:Nch) = 0.0D0
  DO i = ch1, ch2
    CALL radiometry(i,q,r,v,kk)
    IF(FMC%Sideband < 0) freqs(i) = q
    IF(FMC%Sideband > 0) freqs(i) = r
  END DO

  Call get_filters(no_pfa_ch,FMI%no_filt_pts,pfa_ch, &
 &                 FMI%F_grid_filter,freqs,FMI%Filter_func,  &
 &                 FMC%InDir,ier)
  if(ier /= 0) goto 99
!
 99  CLOSE(11,iostat=i)
     CLOSE(13,iostat=i)
!
     if(io /= 0) then
       Ier = iabs(io)
       Call ErrMsg(Line,io)
     endif

     Return

  END SUBROUTINE L2_LOAD

! *****     Private procedures     *************************************
! ------------------------------------------------ Radiometry   -----

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

! ------------------------------------------------     ANTENNA     -----
! This subroutine reads an external antenna aperture autocorrelation
! file
!
  Subroutine ANTENNA (Fn, M, XLAMDA, AAAP, D1AAP, D2AAP, IAS, IER)

    use GET_LUN, only: AAAP_UNIT
    use MACHINE, only: IO_ERROR

    Integer(i4), parameter :: MaxV= 2048
!
    Real(r8), intent(out) :: AAAP(*),D1AAP(*),D2AAP(*),XLAMDA
    integer(i4), intent(out) :: IAS, IER
!
    Integer(i4) :: M, k, lf, I, J, NTR, L, N
    REAL(r8) :: V(6), VALL(MaxV, 6), dx2p, Q, Q2
!
    Character*(*) Fn
!
    ier = 0
    lf = len_trim(Fn)
    OPEN(aaap_unit,FILE=Fn(1:lf),action='READ',STATUS='OLD',iostat=k)
    if(k /= 0) then
      ier = 1
      Print *,'** File: ',Fn(1:lf)
      call io_error ('Opening file in ANTENNA', k, Fn(1:lf) )
      goto 20
    endif
!
    READ(aaap_unit,*,iostat=k) XLAMDA
    if(k /= 0) then
      ier = 1
      Print *,'** Reading error in file: ',Fn(1:lf)
      call io_error ( 'Reading file in ANTENNA', k, Fn(1:lf) )
      goto 20
    endif
!
    dx2p = 6.28318530717959_r8 * Xlamda         ! 2 * Pi * Lambda
!
    i = 0
    k = 0
!
    do while(k == 0)
!
! This loop MUST exit on an end of file condition
!
      read(aaap_unit,*,iostat=k) v
      if(k == -1) exit
      if(k == 0) then
        if(I == MaxV) then
          ier = 1
          Print *,'** Error in ANTENNA subroutine !!'
          Print *,'   Too many lines in: ',Fn(1:lf)
          Print *,'   Maximum allowed is:',MaxV
          goto 20
        endif
        i = i + 1
        vall(i,:) = v
      else
        ier = 1
        Print *,'** Reading error in file: ',Fn(1:lf)
        call io_error ( 'Reading file in ANTENNA', k, Fn(1:lf) )
        goto 20
      endif
!
    end do
!
    ias = i
    ntr = 2**m
!
    do i = 1, ias
!
      j = 2 * i - 1
      l = 2 * ias + j
      n = 4 * ias + j
!
      v = vall(i,:)
!
      aaap(j)   = v(1)
      aaap(j+1) = v(2)
      aaap(l)   = v(3)
      aaap(l+1) = v(4)
      aaap(n)   = v(5)
      aaap(n+1) = v(6)
!
! Below is the  first derivative field:     ! i*Q * F(S), i = Sqrt(-1)
!
      q = (i - 1) * dx2p
      d1aap(j)    = -v(2) * q
      d1aap(j+1)  =  v(1) * q
      d1aap(l)    = -v(4) * q
      d1aap(l+1)  =  v(3) * q
      d1aap(n)    = -v(6) * q
      d1aap(n+1)  =  v(5) * q
!
! Below is the second derivative field:     ! (i*Q)**2 * F(S), i = Sqrt(-1)
!
      q2 = -q * q                         ! (i*q)**2 = -q*q
      d2aap(j)    =  v(1) * q2
      d2aap(j+1)  =  v(2) * q2
      d2aap(l)    =  v(3) * q2
      d2aap(l+1)  =  v(4) * q2
      d2aap(n)    =  v(5) * q2
      d2aap(n+1)  =  v(6) * q2
!
    end do
!
20  CLOSE(aaap_unit,iostat=k)
!
    Return

  End subroutine ANTENNA

! ------------------------------------------------ GET_FILTERS   -----
! This subroutine loads the filter shapes into memory
!
SUBROUTINE get_filters(no_pfa_ch,no_filt_pts,pfa_ch,f_grid_filter, &
               &       freqs,filter_func,InDir,ier)

  use MLSCommon, only: I4, R8
  use GET_LUN, only: filter_unit

!  ===============================================================
!  Declaration of variables for sub-program: get_filters
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: pfa_ch(*)
Integer(i4), INTENT(IN) :: no_pfa_ch, no_filt_pts

Integer(i4), INTENT(OUT) :: ier

Real(r8), INTENT(IN) :: freqs(*)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid_filter(:,:)

Character (LEN=*), INTENT(IN) :: InDir

!  ----------------
!  Local variables:
!  ----------------

Real(r8) :: df, frq, q
Character (LEN=80) :: Fn
Integer(i4) :: j, ch_i, ld, mch, pch

Integer(i4) :: Channel,Nfp
Real(r8) :: Xlhs,Xrhs,N_Filter(200)

NAMELIST/IN/Channel,Nfp,Xlhs,Xrhs,N_Filter

! Begin code:

  ier = 0

  pch = -1
  Fn(1:)=' '
  ld = LEN_TRIM(InDir)
  Fn = InDir(1:ld)//'normalized_filter_banks.dat'
!
  Close(filter_unit,iostat=j)
  Open(filter_unit,file=Fn,status='OLD',action='READ',iostat=ier)
  IF(ier /= 0) goto 99

! Find the species index in the l2pc mixing ratio database:

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    frq = freqs(mch)
    IF(frq < 1.0D0) THEN
      ier = 1
      WRITE(6,900) mch
      GOTO 99
    END IF

! Read in filter's response function

    if(mch < pch) then
      ier = 0
      REWIND(filter_unit,IOSTAT=j)
    endif

    DO
      READ(UNIT=filter_unit,nml=IN,IOSTAT=ier)
      IF(ier > 0) goto 99
      if(Channel == mch .or. ier < 0) EXIT
    END DO

    if(Channel /= mch) then
      Ier = 1
      WRITE(6,905) mch,Fn(1:ld)
      GOTO 99
    endif

    IF(nfp /= no_filt_pts) THEN
      ier = 1
      WRITE(6,910) Fn(1:ld),no_filt_pts,nfp
      GOTO 99
    END IF

    pch = mch
    df = (xrhs-xlhs)/(nfp-1)
    DO j = 1, nfp
      q = xlhs + (j - 1) * df
      f_grid_filter(j,ch_i) = frq + q
    END DO

    filter_func(1:nfp,ch_i) = N_Filter(1:nfp)

  END DO                         ! On ch_i

 99 Close(filter_unit,iostat=j)

    IF(ier /= 0) CALL errmsg(Fn,ier)

 900 FORMAT(' ** Error in get_filters subroutine **',/, &
            '    Inconsistant User Input.',/, &
            '    PFA Channel:',i3,' not among the non-PFA channels !')

 905 FORMAT(' ** Error in get_filters subroutine **',/, &
            '    Channel:',i3,' not found in file: ',A)

 910 FORMAT(' ** Error in get_filters subroutine **',/, &
            '    Inconsistant Input for file:',A,/, &
            '    Input no_filt_pts:',i4,' no_filt_pts in file:',i4)

  Return

END SUBROUTINE get_filters

!=====================================================================

end module L2_LOAD_M
! $Log$
! Revision 1.1  2001/02/22 18:12:05  ZShippony
! Initial conversion to Fortran 90
