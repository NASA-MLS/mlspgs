Program l2_Test
  use GL6P
  use EOS_MDB, only: CS_MNF => MAX_FREQ
  use L2PC_FILE_PARAMETERS, only: DEG2RAD,                  &
                                  mnsv => MAX_NO_SV_ELMNTS, &
                                  mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, &
                                 MAXAITKENPTS, MAXFILTPTS,  &
                                 PFA_SLAB, SPECTRO_PARAM,  MAXPFACH, &
                                 EARTH_MAJOR, EARTH_MINOR, &
                                 MAXGEOM, MAXLINES
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi
  use MLSCommon, only: I4, R4, R8
  use STRINGS, only: STRLWR
  use MDBETA, only: NO_T_PHI
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
      HT2, RR, PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  use FWD_MDL_SET_UP_M, only: FWD_MDL_SET_UP
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use RAD_TRAN_M, only: RAD_TRAN
  use RAD_TRAN_WD_M, only: RAD_TRAN_WD
  use D_HUNT_M, only: HUNT          ! ** DEBUG
  implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: ngt = (Ng+1) * Nlvl
Integer(i4), PARAMETER :: Npath = 2 * ngt

Integer(i4), PARAMETER :: mnf = 25

Logical :: temp_der

Integer(i4) :: p_indx(Nlvl), N_TAN(Nlvl),no_filt_pts,  &
  ier, no_phi_t, no_spectro, SPECT_ATMOS(Nsps),        &
  no_coeffs_f(Nsps), no_phi_f(Nsps), n_obs, NO_SPS_TBL(6), no_pfa_ch,    &
  SPECT_INDEX(Nsps), ATMOS_INDEX(Nsps), NDX_SPS(Nsps,MAXPFACH),          &
  SPS_TBL(Nsps,6),pfa_ch(MAXPFACH),NO_INT_FRQS(MAXPFACH),NO_PHI_VEC(mnsv)

Integer(i4) :: i, j, k, kk, ht_i, no_t, mnz, no_geom, no_tan_hts, ld, &
               no_atmos, ch1, ch2, geo_i, n_lvls, si, mfi, band, n_sps, &
               no_freqs, ptg_i, frq_i, io, klo, khi, m, l, brkpt, no_ele, &
               mid, ilo, ihi

Type(path_index)  :: ndx_path(Nptg)
Type(path_vector) :: z_path(Nptg),t_path(Nptg),h_path(Nptg), &
                     n_path(Nptg),phi_path(Nptg),dhdz_path(Nptg), &
                     spsfunc_path(Nsps,Nptg)

Type(path_derivative) :: dh_dt_path(Nptg)

Real(r8) :: href(Nlvl), zref(Nlvl), t_tan(Nptg), tan_hts(Nptg),  &
            tan_hts_raw(Nptg), tan_press(Nptg), tan_temp(Nptg), &
            t_coeff(mxco,mnp), t_z_basis(mxco), t_phi_basis(mnp), &
            tan_dh_dt(Nlvl,mxco),z_grid(Nlvl), Frq, h_tan, Rad, &
            t_grid(Nlvl),h_grid(Nlvl),z_gnlv(400), e_rad, cse, Rs, &
            ref_corr(N2lvl,Nptg), t_script(N2lvl), tau(N2lvl)

Logical :: IS_F_LOG(Nsps)

real(r8) :: C_FREQ(Nch), freq_grid(mnf)
real(r8) :: MR_F(mxco,mnp,Nsps)
real(r8) :: F_BASIS(mxco,Nsps), PHI_BASIS_F(mnp,Nsps)
real(r8) :: F_GRID(maxaitkenpts,maxpfach)
real(r8) :: F_GRID_FILTER(maxfiltpts,maxpfach)
real(r8) :: FILTER_FUNC(maxfiltpts,maxpfach)

Real(r4) :: K_TEMP(Nptg,mxco,mnp)
Real(r4) :: K_ATMOS(Nptg,mxco,mnp,Nsps)
Real(r4) :: K_SPECT_DW(Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DN(Nptg,mxco,mnp,Nsps),  &
            K_SPECT_DNU(Nptg,mxco,mnp,Nsps)
!
Type(path_beta) :: beta_path(Nsps,mnf,Nptg)

Real(r8) :: S_TEMP, H_OBS, EARTH_REF, ZRoC, geoc_lat, q, r
Real(r8) :: RadV(Nptg)
!
Character :: Primag
Character (LEN=40) :: Ax, Fdata
Character (LEN=80) :: InDir, Fnd

Type (l2pc_header_one) :: HEADER1
Type (atmos_comp) :: ATMOSPHERIC(Nsps)
Type (geom_param) :: GEOMETRIC(maxgeom)
Type (spectro_param) :: SPECTROSCOPIC(3*Nsps)
Type (pfa_slab)   :: PFA_SPECTRUM(6,Nsps)

!  ----------------------
! Read the inputs from a file...
!
  ier = 0
  Fdata(1:) = ' '
  Fdata = '/user5/zvi/seez'     ! SUN
! Fdata = '/home/zvi/seez'      ! PC
  CLOSE(11,iostat=io)
  OPEN(11,file=Fdata,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99

  no_int_frqs(1:maxpfach) = 0

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = 0
  read(11,*,iostat=io) kk
  if(io /= 0) goto 99

  HEADER1%NO_POINTINGS = kk
  read(11,*,iostat=io) (HEADER1%sv_components(i),i=1,kk)
  if(io /= 0) goto 99
  read(11,*,iostat=io) (HEADER1%no_elmnts_per_sv_component(i),i=1,kk)
  if(io /= 0) goto 99
  read(11,*,iostat=io) (HEADER1%sv_component_first_elmnt_index(i),i=1,kk)
  if(io /= 0) goto 99
  read(11,*,iostat=io) HEADER1%no_bands
  if(io /= 0) goto 99
  read(11,*,iostat=io) HEADER1%no_channels_per_band
  if(io /= 0) goto 99
  read(11,*,iostat=io) HEADER1%no_sv_components
  if(io /= 0) goto 99
  read(11,*,iostat=io) HEADER1%no_coeff_per_component
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

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

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_geom
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  do i = 1, no_geom
    read(11,*,iostat=io) GEOMETRIC(i)
    if(io /= 0) goto 99
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_atmos
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (atmos_index(i),i=1,no_atmos)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = mxco
  do j = 1, no_atmos
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    atmospheric(j)%NAME = AdjustL(Ax)
    read(11,*,iostat=io) k, ht_i
    if(io /= 0) goto 99
    atmospheric(j)%SPECTAG = k
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

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99
  do i = 1, n_sps
    read(11,*,iostat=io) pfa_spectrum(band,i)%NO_SPS
    if(io /= 0) goto 99
    read(11,*,iostat=io) pfa_spectrum(band,i)%NO_LINES
    if(io /= 0) goto 99
    read(11,*,iostat=io) pfa_spectrum(band,i)%SPS_SPECTAG
    if(io /= 0) goto 99
  end do

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_spectro, mfi
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  SPECT_INDEX(1:Nsps) = -1
  SPECT_ATMOS(1:Nsps) = -1

  read(11,*,iostat=io) (spect_index(i),i=1,no_spectro)
  if(io /= 0) goto 99
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  kk = mxco
  do j = 1, no_spectro
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    spectroscopic(j)%TYPE = AdjustL(Ax)
    Ax(1:)=' '
    read(11,'(A)',iostat=io) Ax
    if(io /= 0) goto 99
    spectroscopic(j)%NAME = AdjustL(Ax)
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    spectroscopic(j)%SPECTAG = ht_i
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    spectroscopic(j)%NO_PHI_VALUES = ht_i
    read(11,*,iostat=io) ht_i
    if(io /= 0) goto 99
    spectroscopic(j)%NO_ZETA_VALUES = ht_i
    read(11,*,iostat=io) (spectroscopic(j)%DER_CALC(i),i=1,6)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (spectroscopic(j)%PHI_BASIS(i),i=1,mfi+2)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (spectroscopic(j)%ZETA_BASIS(i),i=1,kk+2)
    if(io /= 0) goto 99
  end do
!
! Create spect_atmos array:
!
  do k = 1, no_atmos
    kk = atmospheric(k)%Spectag
    do i = 1, no_spectro
      if(kk == spectroscopic(i)%Spectag) then
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

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  Fnd(1:)=' '
  read(11,'(A)',iostat=io) Fnd
  if(io /= 0) goto 99
  Fnd = AdjustL(Fnd)

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) si,n_lvls,no_tan_hts,mnz,no_filt_pts
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) q, r, phi_tan, geoc_lat
  if(io /= 0) goto 99

  href(1:n_lvls) = q
  zref(1:n_lvls) = r
!
! Add phi_tan to the spectroscopic phi's and convert to radiance:
!
  do j = 1, no_spectro
    i = spectroscopic(j)%NO_PHI_VALUES
    spectroscopic(j)%PHI_BASIS(1:i) = (spectroscopic(j)%PHI_BASIS(1:i) + &
   &                                   phi_tan) * deg2rad
  end do
!
  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (p_indx(i),i=1,n_lvls)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (z_gnlv(i),i=1,mnz)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (tan_hts_raw(i),i=1,no_tan_hts)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_t, no_phi_t
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (t_z_basis(i),i=1,no_t)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (t_phi_basis(i),i=1,no_phi_t)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) ((t_coeff(i,j),j=1,no_phi_t),i=1,no_t)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (no_phi_f(i),i=1,no_atmos)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (no_coeffs_f(i),i=1,no_atmos)
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (is_f_log(i),i=1,no_atmos)     ! ** A a new one !
  if(io /= 0) goto 99

  DO geo_i = 1, no_atmos
    kk = no_phi_f(geo_i)
    ht_i = no_coeffs_f(geo_i)
    read(11,*,iostat=io) (f_basis(i,geo_i),i=1,ht_i)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (phi_basis_f(i,geo_i),i=1,kk)
    if(io /= 0) goto 99
    read(11,*,iostat=io) ((mr_f(i,j,geo_i),j=1,kk),i=1,ht_i)
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

! Get the selected integration grid pressures. Also, define the GL
! pressure grid:

  DO i = 1, n_lvls
    j = p_indx(i)
    z_grid(i) = z_gnlv(j)
  END DO
  z_grid(n_lvls+1) = z_grid(n_lvls)
!
  Primag = 'p'
  Call fwd_mdl_set_up(Primag,href,zref,e_rad,z_grid,t_grid,h_grid,    &
 &     ndx_path,n_lvls,geometric,no_geom,t_z_basis,t_coeff,no_t,      &
 &     geoc_lat,atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,c_freq, &
 &     si,earth_ref,s_temp,h_obs,n_obs,tan_hts_raw,tan_hts,tan_press, &
 &     tan_temp,no_tan_hts,no_sps_tbl,sps_tbl,no_pfa_ch,no_filt_pts,  &
 &     pfa_ch,no_int_frqs,pfa_spectrum,t_tan,ndx_sps,n_tan,           &
 &     f_grid_filter,f_grid,filter_func,ch1,ch2,InDir,ld,             &
 &     header1,no_phi_f,phi_basis_f,atmos_index,no_phi_vec,z_path,    &
 &     h_path,t_path,phi_path,n_path,dhdz_path,dh_dt_path,Npath,      &
 &     no_phi_t,t_phi_basis,ZRoC,spectroscopic,no_spectro,            &
 &     spect_index,tan_dh_dt,no_freqs,freq_grid,spsfunc_path,         &
 &     is_f_log,beta_path,ier)
  IF(ier /= 0) goto 99
!
! Compute the refraction correction scaling matrix:
!
  CALL refraction_correction(no_tan_hts, tan_hts, h_path, n_path, &
 &                           ndx_path, E_rad, ref_corr)
!
! Now, Compute the radiances alone:
!
  RadV(:Nptg) = 0.0
  temp_der = .true.
! temp_der = .false.
!
! *** DEBUG
     r = -1.66667
     Call Hunt(r,tan_press,no_tan_hts,klo,khi)
     IF(ABS(r-tan_press(khi)) < ABS(r-tan_press(klo))) klo = khi
! *** END DEBUG

  do ptg_i = 1, no_tan_hts-1
!
    k = ptg_i
    h_tan = tan_hts(k)
!
!   do frq_i = 1, no_freqs
    do frq_i = 1, 1               ! ** DEBUG
!
      Frq = freq_grid(frq_i)
!
      Call Rad_Tran(Frq,N_lvls,h_tan,band,n_sps, sps_tbl, ndx_path(k), &
      &    z_path(k), h_path(k), t_path(k), phi_path(k), dHdz_path(k), &
      &    earth_ref, beta_path(1:,frq_i,k), spsfunc_path(1:,k),       &
      &    ref_corr(1:,k), s_temp, brkpt, no_ele, mid, ilo, ihi, cse,  &
      &    Rs, t_script, tau, Rad, Ier)
      IF(ier /= 0) goto 99
!
      if(frq_i == 1) RadV(ptg_i) = Rad
!
! Now, Compute the radiances derivatives:
!
      CALL Rad_Tran_WD(ptg_i,frq_i,Frq,N_lvls,band,n_sps,sps_tbl,    &
   &       z_path(k), h_path(k), t_path(k), phi_path(k), dHdz_path(k),     &
   &       atmospheric, beta_path(1:,frq_i,k), spsfunc_path(1:,k),  &
   &       t_z_basis, f_basis, no_coeffs_f, mr_f, no_t, ref_corr(1:,k),    &
   &       no_phi_f, phi_basis_f, temp_der, no_phi_t, t_phi_basis,         &
   &       dh_dt_path(k), k_atmos, k_temp, spect_atmos,            &
   &       spectroscopic, k_spect_dw, k_spect_dn, k_spect_dnu,             &
   &       is_f_log, brkpt, no_ele, mid, ilo, ihi, t_script, tau, &
   &       Ier)
      IF(ier /= 0) goto 99
!
    end do
!
  end do
!
  RadV(no_tan_hts) = RadV(no_tan_hts-1)
!
! *** DEBUG Print
!
     k = 61
     j = no_tan_hts
     kk = j - si + 1
     write(*,'(''pfa_rad'',a1,i2.2)') char(92),kk
     write(*,905) (RadV(k),k=si,j)
905  format(4(2x,1pg15.8))
!
     frq_i = 1
     k = 18003
!
     m = -1
     do j = 1, n_sps
       if(atmospheric(j)%SPECTAG == k) then
         m = j
         EXIT
       endif
     end do

!   if(m < -100) then               ! Bypass print
    if(m > 0) then
      Ax(1:)=' '
      k = no_coeffs_f(m)
      Ax = 'k_'//atmospheric(m)%NAME
      l = Len_Trim(Ax)
      Ax = Ax(1:l)//'_'
      Call StrLwr(Ax)
      j = (no_phi_f(m)+1)/2
      Print 900,Ax(1:l+1),frq_i,k
      Print 901,(k_atmos(klo,i,j,m),i=1,k)
    endif
!
    j = (no_phi_t+1)/2
    Print 900,'k_temp_',frq_i,no_t
    Print 901,(k_temp(klo,i,j),i=1,no_t)
!
    kk = spectroscopic(1)%no_phi_values
    ht_i = spectroscopic(1)%no_zeta_values
    Print 900,'k_spect_dw_',frq_i,ht_i*kk
    r = SUM(k_spect_dw(klo,1:ht_i,1:kk,1))
    do k = 1, ht_i
      Print 901,(k_spect_dw(klo,k,j,1),j=1,kk)
    end do
    Print *,'  Sum over all zeta & phi coeff:',sngl(r)
!
    kk = spectroscopic(2)%no_phi_values
    ht_i = spectroscopic(2)%no_zeta_values
    r = SUM(k_spect_dn(klo,1:ht_i,1:kk,1))
    Print 900,'k_spect_dn_',frq_i,ht_i*kk
    do k = 1, ht_i
      Print 901,(k_spect_dn(klo,k,j,1),j=1,kk)
    end do
    Print *,'  Sum over all zeta & phi coeff:',sngl(r)
!
    kk = spectroscopic(3)%no_phi_values
    ht_i = spectroscopic(3)%no_zeta_values
    r = SUM(k_spect_dnu(klo,1:ht_i,1:kk,1))
    Print 900,'k_spect_dnu_',frq_i,ht_i*kk
    do k = 1, ht_i
      Print 901,(k_spect_dnu(klo,k,j,1),j=1,kk)
    end do
    Print *,'  Sum over all zeta & phi coeff:',sngl(r)
!
900 format(a,i2.2,'\',i3.3)
901 format(5(2x,1pe13.6))
!

 99  CLOSE(11,iostat=i)
     if(io /= 0) Call ErrMsg(' ',io)

     Stop

contains

Function I2C(n)

Integer n
Character(LEN=12) I2C,Zx

   Zx(1:)=' '
   write(Zx,*) n
   I2C = AdjustL(Zx)

end Function I2C

end Program l2_Test
