Program l2_numder
  use GL6P
  use EOS_MDB
  use L2PC_FILE_PARAMETERS, mnsv => MAX_NO_SV_ELMNTS, &
                            mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES
  use L2PC_FILE_STRUCTURES
  use L2PCdim
  use MLSCommon
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use STRINGS
  use MDBETA
  use ELLIPSE
  use GEO_GRIDS_M
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use FILTER_SW_M, only: FILTER
  use DSIMPSON_MODULE, only: DSIMPS
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  use REFRACTION_M
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use DO_T_SCRIPT_M, only: DO_T_SCRIPT
  use FAST_DELTA_M, only: FAST_DELTA
  use FAST_ZOPACITY_M, only: FAST_ZOPACITY
  use SCRT_DN_M, only: SCRT_DN

  implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: ngt = (Ng+1) * Nlvl
Integer(i4), PARAMETER :: Npath = 2 * ngt

Integer(i4), PARAMETER :: mnf = 25

Logical :: temp_der

Integer(i4) :: p_indx(Nlvl), N_TAN(Nlvl),no_filt_pts,  &
  jpm, ier, no_phi_t, no_spectro, SPECT_ATMOS(Nsps),        &
  no_coeffs_f(Nsps), no_phi_f(Nsps), n_obs, NO_SPS_TBL(6), no_pfa_ch,    &
  SPECT_INDEX(Nsps), ATMOS_INDEX(Nsps), NDX_SPS(Nsps,MAXPFACH),          &
  SPS_TBL(Nsps,6),pfa_ch(MAXPFACH),NO_INT_FRQS(MAXPFACH),NO_PHI_VEC(mnsv)

Integer(i4) :: i, j, k, kk, ht_i, no_t, mnz, no_geom, no_tan_hts, ld, &
               no_atmos, ch1, ch2, geo_i, n_lvls, si, mfi, band, n_sps, &
               no_freqs, ptg_i, frq_i, io, brkpt, no_ele, mid, ilo, ihi

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
! Fdata = '/hpusr/zvi/new_seez'     ! SUN
  Fdata = '/home/zvi/new_seez'      ! PC
  CLOSE(11,iostat=io)
  OPEN(11,file=Fdata,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99

33  jpm = 0
    Call Gti('Enter plus/minus (p/m)',Ax)
    if(Ax < '!') goto 99
    if(Ax(1:1) == 'p') then
      jpm = 1
    else if(Ax(1:1) == 'm') then
      jpm = -1
    else
      goto 33
    endif

  Print *
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
 &     is_f_log,beta_path,jpm,ier)
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
    end do
!
  end do
!
  RadV(no_tan_hts) = RadV(no_tan_hts-1)
!
     k = 61
     j = no_tan_hts
     kk = j - si + 1
     write(*,'(''pfa_rad'',a1,i2.2)') char(92),kk
     write(*,905) (RadV(k),k=si,j)
905  format(4(2x,1pg15.8))
!
 99  CLOSE(11,iostat=i)
     if(io /= 0) Call ErrMsg(' ',io)

     Stop

contains
!----------------------------------------------------------------

Function I2C(n)

Integer n
Character(LEN=12) I2C,Zx

   Zx(1:)=' '
   write(Zx,*) n
   I2C = AdjustL(Zx)

end Function I2C

!----------------------------------------------------------------

SUBROUTINE fwd_mdl_set_up(primag,href,zref,e_rad,z_grid,t_grid,h_grid,  &
           ndx_path,n_lvls,geometric,no_geom,t_z_basis,t_coeff,no_t,    &
           geoc_lat,atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,freq, &
           si,earth_ref,s_temp,h_obs,n_obs,tan_hts_raw,tan_hts,tan_press, &
           tan_temp,no_tan_hts,no_sps_tbl,sps_tbl,no_pfa_ch,no_filt_pts, &
           pfa_ch,no_int_frqs,pfa_spectrum,t_tan,ndx_sps,n_tan_ptr,     &
           f_grid_fltr,f_grid,fltr_func,ch1,ch2,InDir,ld,hdr1,          &
           no_phi_f,phi_basis_f,atmos_index,no_phi_vec,z_path,h_path,   &
           t_path,phi_path,n_path,dhdz_path,dh_dt_path,npath,           &
           no_phi_t,t_phi_basis,rad_cur,spectroscopic,no_spectro,       &
           spect_index,tan_dh_dt,no_freqs,freq_grid,spsfunc_path,       &
           is_f_log,beta_path,jpm,ier)

!  ===============================================================
!  Declaration of variables for sub-program: fwd_mdl_set_up
!  ===============================================================

Integer(i4), PARAMETER :: ngt = (Ng+1) * Nlvl

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_atmos, no_t, ch1, ch2, si,  &
             atmos_index(:), spect_index(:), no_phi_t, no_freqs, &
             no_coeffs_f(:), no_phi_f(:), jpm, &
             no_spectro, npath, no_pfa_ch, no_filt_pts, pfa_ch(:), &
             no_int_frqs(:)
!
Integer(i4), INTENT(IN OUT) :: no_tan_hts, no_geom, ld

Integer(i4), INTENT(OUT) :: n_obs, no_sps_tbl(:), sps_tbl(:,:),  &
             no_phi_vec(:), ier, ndx_sps(:,:), n_tan_ptr(:)
!
Real(r8), INTENT(IN) :: z_grid(:), href(:), zref(:), freq_grid(:), &
          t_phi_basis(:), t_coeff(:,:), t_z_basis(:), &
          mr_f(:,:,:), phi_basis_f(:,:), f_basis(:,:)
!
Real(r8), INTENT(OUT) :: tan_hts(:), tan_dh_dt(:,:)
Real(r8), INTENT(IN OUT) :: tan_hts_raw(:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:)
!

Type(path_index) , INTENT(OUT) :: ndx_path(:)
Type(path_vector), INTENT(OUT) :: z_path(:),t_path(:),h_path(:), &
                   phi_path(:), n_path(:), dhdz_path(:), spsfunc_path(:,:)

Type(path_derivative), INTENT(OUT) :: dh_dt_path(:)

Real(r8), INTENT(OUT) :: e_rad,earth_ref,s_temp,h_obs,freq(:),geoc_lat, &
          tan_temp(:),t_grid(:),h_grid(:),tan_press(:),t_tan(:),rad_cur, &
          f_grid(:,:),fltr_func(:,:),f_grid_fltr(:,:)

Logical, INTENT(IN) :: is_f_log(*)
!
Character (LEN=*), INTENT(IN) :: InDir
Character (LEN=*), INTENT(IN) :: primag
!  ----------------------
!  Local variables:
!  ----------------

Logical :: wet
Integer(i4) :: i, j, k, jp, ncpb, ich, sps_i, geo_i, band, geo_j, &
               ih2o, kk, no_geom_new, band1, band2, gl_count

Real(r4) :: def_val(maxgeom) = (/                                    &
            0.0, 0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.0, 6972.0, 6372.0, &
            0.05, 2.735, 1.0 /)

Real(r8) :: z_glgrid(ngt)

! Real(r8) :: h2o_path(Npath,Nptg)
Type(path_vector) :: h2o_path(Nptg)

Real(r8) :: q,r,p,g,sw,cw,rp,b2,z1,z2,xm,ym,beta_inc,incl,zeta,phi

Character (LEN=8) :: geomname(maxgeom) = (/                         &
         &  'ELEV_183','ELEV_205','AZIM_183','AZIM_205','AZIM_REF', &
         &  'ROLL    ','PITCH   ','YAW     ','GEOCSRAD','GEOCERAD', &
         &  'EARTHREF','SPACE_T ','LOOK_DIR'/)

!  PFA variables:

type (l2pc_header_one), intent(in)  :: hdr1

type (geom_param),    intent(inout) :: geometric(*)
type (atmos_comp),    intent(inout) :: atmospheric(*)
type (spectro_param), intent(inout) :: spectroscopic(*)
type (pfa_slab), intent(inout)      :: pfa_spectrum(6,*)
!
! Set up the default geometric linearization values of earth radius,earth
! emission,background space temperature and satellite position.
! Compute the default geocentric earth radius

  ier = 0

  q = DBLE(earth_major)
  a2 = q * q

  r = DBLE(earth_minor)
  b2 = r * r

  beta_inc = 98.0D0                ! Orbit inclination angle for EOS-MLS
  incl = (beta_inc - 90.0D0) * deg2rad
  q = TAN(incl)
  r = q * q
  c2 = (1.0D0+r)*a2*b2/(a2+b2*r)   ! This is c*c
  q = SQRT(c2)                     ! This is c, Minor axis for 2D ellipse

  c2oa2 = c2 / a2

! Get the Tangent Phi (Geodetic Lat.) from the Geocentric Lat.
!
! r = geoc_lat * deg2rad
! q = SIN(r)
! sw = q * q
! r = COS(incl)
! cw = r * r
! q = a2*a2*sw*cw/(b2*b2+(a2*a2-b2*b2)*sw*cw)
! spt = SQRT(q)                ! Sin(Phi_tan)
! phi_tan = DASIN(spt)
!
! Get the Geocentric Lat. (geoc_lat) from the Tangent Phi (Geodetic Lat.)
!
  r = Phi_tan
  Phi_tan  = r * deg2rad               ! Convert to Radians
  spt = SIN(Phi_tan)
  cpt = COS(Phi_tan)
  cw = cpt * cpt
  sw = spt * spt
  r = a2*a2*cw +b2*b2*sw
  q = spt*b2/Sqrt(r)/COS(incl)
  geoc_lat = DASIN(q) / deg2rad        ! In Degrees

  nphi_tan = a2 / SQRT(c2-(c2-a2)*cw)
  ps = -1.0D0

!  Compute Radius of Curvature Circle (RoC) and its center coordinates:

  r = c2 * cpt / a2
  q = spt * spt + r * r
  roc = nphi_tan * SQRT(q)
  xoc = (nphi_tan - roc) * cpt
  yoc = (c2oa2 * nphi_tan - roc) * spt
  rad_cur = roc

!  Compure Earth Radius (Elliptical)

  q = ((a2*a2)*cw+(b2*b2)*sw)/(a2*cw+b2*sw)
  rp = SQRT(q)

  e_rad = rp
  def_val(10) = e_rad

! Default earth emission

  earth_ref = def_val(11)

! Default satellite position

  h_obs = def_val(9)

! Default space temperature

  s_temp = def_val(12)

! From the selected integration grid pressures define the GL pressure
! grid:

  gl_count = 0
  z2 = z_grid(1)
  DO i = 2, n_lvls
    z1 = z2
    z2 = z_grid(i)
    xm = 0.5D0 * (z2 + z1)
    ym = 0.5D0 * (z2 - z1)
    gl_count = gl_count + 1
    z_glgrid(gl_count) = z1
    DO j = 1, ng
      gl_count = gl_count + 1
      z_glgrid(gl_count) = xm + ym * Gx(j)
    END DO
  END DO

  gl_count = gl_count + 1
  z_glgrid(gl_count) = z2
!
! Initialize the no_phi_vec array. This array stores the number of Phi's for
! each state_vector type entry. Initial value for all is: 1.
! (Later, we will connect between 'no_phi_vec' and 'no_phi_f' and 'no_phi_g'
! entries.)

  DO i = 1, max_no_sv_elmnts
    no_phi_vec(i) = 1
  END DO

! Scan geometric quantities for specific values of the above if user wants
! to override the defaults

  no_geom_new = no_geom
  DO geo_i = 1, maxgeom

! Search the geometric parameters for inclusion of "vital" quantities
! Insert defaults if these do not presently exist

    geo_j = 1
    DO WHILE(geometric(geo_j)%name /= geomname(geo_i) .AND.  &
          geo_j < no_geom )
      geo_j = geo_j + 1
    END DO

    IF(geometric(geo_j)%name /= geomname(geo_i) ) THEN

! Meaning it did not find it

      no_geom_new = no_geom_new + 1
      geometric(no_geom_new)%name = geomname(geo_i)

      DO i = 1, max_no_bands
        geometric(no_geom_new)%der_calc(i) = .false.
      END DO

      geometric(no_geom_new)%lin_val = def_val(geo_i)

      IF(geomname(geo_i) /= 'LOOK_DIR') THEN
        Print *,' WARNING: ',geometric(no_geom_new)%name, &
     &          ' NOT FOUND IN USER INPUTS'
        Print *,' ASSUMING ',geometric(no_geom_new)%lin_val, &
     &          ' FOR THIS QUANTITY'
      END IF

    ELSE

! It is there and set up internal parameters accordingly
! If applicable place the earth radius in the state vector

      IF(geometric(geo_j)%name == 'GEOCERAD') THEN

! Force geometric earth radius value to be value appropriate for the
! latitude bin

        geometric(geo_j)%lin_val = def_val(10)

!           E_rad = geometric(geo_j)%lin_val

! Look for earth emission

      ELSE IF(geometric(geo_j)%name == 'EARTHREF') THEN

        earth_ref = geometric(geo_j)%lin_val

! Look for the background space temperature

      ELSE IF(geometric(geo_j)%name == 'SPACE_T') THEN

        s_temp = geometric(geo_j)%lin_val

! Look for satellite position

      ELSE IF(geometric(geo_j)%name == 'GEOCSRAD') THEN

        h_obs = geometric(geo_j)%lin_val

      END IF

    END IF

  END DO

  no_geom = no_geom_new

! Call the grids program to get the preselected integration heights and
! pointings. Also transfer the pointings expressed in pressure units to
! heights.  Also, set up geophysical parameters (for derivatives):
! [The following routine replaces the older code using the two routines:
!  grids() and geo_basis() ]

! (z_grid,t_grid,h_grid) are the arrays of the preselected integration grid.

  CALL geo_grids(z_glgrid,gl_count,z_grid,t_grid,h_grid,tan_dh_dt,  &
       n_lvls,t_z_basis,t_coeff,si,tan_hts,tan_hts_raw,tan_press,   &
       tan_temp,no_tan_hts,z_path,h_path,t_path,phi_path,dhdz_path, &
       dh_dt_path,Npath,g,href,zref,phi_tan,geoc_lat,roc,n_tan_ptr, &
       t_tan,no_t,no_phi_t,t_phi_basis,ndx_path,ier)
  IF(ier /= 0) RETURN
! Get the crossections if this is the first Call to the subroutine,
! The subroutine knows that this is the first Call if the time counter
! in the main program is one

  DO ich = 1, nch
    freq(ich) = 0.0D0
  END DO

  DO ich = ch1, ch2
    CALL radiometry(ich,q,r,p,kk)       ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'p') freq(ich) = q     ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'i') freq(ich) = r     ! DEBUG, Added Jan/23/2000, Z.S
  END DO

  ncpb = hdr1%no_channels_per_band
  band1 = (ch1 + ncpb - 1) / ncpb     ! Begining Band
  band2 = (ch2 + ncpb - 1) / ncpb     !  Ending  Band

  DO band = band1, band2
    kk = 0
    DO i = 1, no_atmos
      IF(atmospheric(i)%fwd_calc(band)) THEN
        kk = kk + 1
        sps_tbl(kk,band) = i
      END IF
    END DO
    no_sps_tbl(band) = kk
  END DO
!
! Create the specie function along the path for all species
!
  DO k = 1, no_tan_hts
    gl_count = ndx_path(k)%total_number_of_elements
    DO band = band1, band2
      do sps_i = 1, no_sps_tbl(band)
        j = sps_tbl(sps_i,band)
        jp = no_phi_f(j)
        kk = no_coeffs_f(j)
        ALLOCATE(spsfunc_path(j,k)%values(gl_count),STAT=i)
        IF(i /= 0) THEN
          ier = i
          PRINT *,'** Error: ALLOCATION error for spsfunc_path ..'
          PRINT *,'   STAT =',ier
          RETURN
        ENDIF
        do i = 1, gl_count
          zeta = z_path(k)%values(i)
          phi = phi_path(k)%values(i)
          Call TWO_D_POLATE (f_basis(1:,j), mr_f(1:,1:,j), kk,    &
     &                       phi_basis_f(1:,j), jp, zeta, phi, q)
          if (is_f_log(j)) q = exp(q)
          spsfunc_path(j,k)%values(i) = q
        end do
      end do
    end do
  END DO

! Set n_obs to be: N_lvls always.
!   (Changed, Aug/6/96 Z.Shippony & W.G.Read)

  n_obs = n_lvls

! Connect between 'no_phi_vec' and 'no_phi_f' entries.

  DO i = 1, no_atmos
    j = atmos_index(i)
    kk = hdr1%sv_component_first_elmnt_index(j)
    no_phi_vec(kk) = no_phi_f(i)
  END DO

! Connect between 'no_phi_vec' and number of spectral phi's entries.

  DO i = 1, no_spectro
    j = spect_index(i)
    kk = hdr1%sv_component_first_elmnt_index(j)
    no_phi_vec(kk) = spectroscopic(i)%no_phi_values
  END DO

!  *** Special PC code - check if number of Keys will exceed maximum:

  kk = 3                ! x_star key + i_star key + Ptan63 Key
  DO geo_i = 1, no_geom
    IF(geometric(geo_i)%der_calc(band1)) kk = kk + 1
  END DO

  DO sps_i = 1, no_atmos
    j = no_phi_f(sps_i) * no_coeffs_f(sps_i)
    kk = kk + j
  END DO

  DO i = 1, no_spectro
    j = spectroscopic(i)%no_phi_values
    j = j * spectroscopic(i)%no_zeta_values
    kk = kk + j
  END DO

  IF(kk >= max_no_key_addr) THEN
    ier = 1
    PRINT *,'** Error in fwd_mdl_set_up subroutine **'
    PRINT *,'   Number of records exceeded maximum. kk =',kk
    PRINT *,'   Maximum number of records allowed =',max_no_key_addr
    RETURN
  END IF

!  *** End Special PC code

! Compute the relative refractive index minus one.
! Get the water mixing ratio function

  sps_i = 1
  DO WHILE (atmospheric(sps_i)%name /= 'H2O'  .AND.  &
        sps_i <= no_atmos)
    sps_i = sps_i + 1
  END DO

  IF (atmospheric(sps_i)%name == 'H2O') THEN
    j = sps_i
    ih2o = sps_i
  ELSE             ! Ignore water contribution to refractive index (dry air)
    j = 1
    ih2o = 0
  END IF

  wet = (ih2o > 0)
  jp = no_phi_f(j)
  kk = no_coeffs_f(j)
  CALL refractive_index(mr_f(1:,1:,j),f_basis(1:,j),phi_basis_f(1:,j), &
  &                     kk,jp,ndx_path,z_path,t_path,phi_path,n_path,  &
  &                     h2o_path,no_tan_hts,wet)

! Create filter grids & functions for PFA calculations

  IF(no_pfa_ch > 0) THEN
    CALL pfa_prep(atmospheric,no_atmos,no_pfa_ch,no_filt_pts,pfa_ch,   &
         no_int_frqs,pfa_spectrum,ndx_sps,f_grid_fltr,f_grid,freq,     &
         fltr_func,no_tan_hts,no_freqs,ndx_path,z_path,t_path,   &
         beta_path,freq_grid,InDir,ld,primag,jpm,ier)
    IF(ier /= 0) RETURN
  ENDIF

  RETURN
END SUBROUTINE fwd_mdl_set_up

!-----------------------------------------------------------------------

SUBROUTINE radiometry(ch, f_p, f_i, db_fi, lmt)

! This subroutine calculates the center frequency of the primary and image
! sideband by channel. It also returns the bandwidth limits of integration
! and the gain of the primary sideband relative to the image in db units.

INTEGER(i4), INTENT(IN) :: ch
INTEGER(i4), INTENT(OUT) :: lmt

REAL(r8), INTENT(OUT) :: f_p
REAL(r8), INTENT(OUT) :: f_i
REAL(r8), INTENT(OUT) :: db_fi

LOGICAL, save :: sgn_fp(6) = (/                                    &
                 .false., .true., .false., .true., .true., .false./)

INTEGER(i4) :: band, sub_ch, j

Real(r4), save :: db_fi_data(90) =(/                   &
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

REAL(r8), save :: f_prime(6) = (/                          &
    63568.418D0, 204352.161D0, 204574.627D0, 206132.067D0, &
    183310.062D0, 184377.788D0/)

REAL(r8), save :: f_image(6) = (/                          &
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

  RETURN
END SUBROUTINE radiometry

!-----------------------------------------------------------------

SUBROUTINE pfa_prep(atmospheric,no_atmos,no_pfa_ch,no_filt_pts,pfa_ch, &
           no_int_frqs,pfa_spectrum,ndx_sps,f_grid_filter,f_grid,freqs,   &
           filter_func,no_tan_hts,no_freqs,ndx_path,z_path,t_path,  &
           beta_path,freq_grid,InDir,ld,primag,jpm,ier)

!  ===============================================================
!  Declaration of variables for sub-program: pfa_prep
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_atmos, no_pfa_ch, no_filt_pts, pfa_ch(:), &
                           no_int_frqs(:), no_freqs, no_tan_hts, jpm

Integer(i4), INTENT(OUT) :: ier, ld, ndx_sps(:,:)

Real(r8), INTENT(IN) :: freqs(:),freq_grid(:)

Type(path_index), INTENT(IN) :: ndx_path(*)
Type(path_vector), INTENT(IN) :: z_path(*),t_path(*)

Type (atmos_comp), INTENT(IN) :: atmospheric(*)

Type (pfa_slab), INTENT(INOUT) :: pfa_spectrum(6,*)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid(:,:), f_grid_filter(:,:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:)  ! (sps_i,frq_i,ptg_i)

Character (LEN=*), INTENT(IN) :: InDir, primag

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: Tiny = epsilon(freqs(1))
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, m, ch_i, sps_ind, sps_i, pb, no_sps, spectag, &
               j4, mch, band, ptg_i, frq_i, spectags(MAXSPS)

Real(r8) :: xlhs, xrhs, df, q, area, frq

Character (LEN=8) :: sps_name

! Type (eos_mdb_hdr) :: MDB_HDR(MAXSPS)
! Type (eos_mdb_rec) :: MDB_REC(MAX_NO_LINES)
Type (eos_mdb_hdr) :: MDB_HDR(02)      ! ** DEBUG
Type (eos_mdb_rec) :: MDB_REC(02)      ! ** DEBUG

Real(r8), DIMENSION(:), ALLOCATABLE :: values, t_power, dbeta_dw, &
                                       dbeta_dn, dbeta_dnu

! begin code:

  ier = 0

!  Read in the EOS Beta database for all species involved

  pb = 0
  spectags(1:MAXSPS) = -1

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = (mch + 14) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

sps:DO sps_i = 1, no_sps

      Spectag = pfa_spectrum(band,sps_i)%sps_spectag

!  Build: ndx_sps, the index in the mixing ratio array:

      sps_ind = 1
      DO WHILE(sps_ind < no_atmos .AND. spectag /=  &
            atmospheric(sps_ind)%spectag)
        sps_ind = sps_ind + 1
      END DO

!  Check if sps is in the database

      IF(spectag /= atmospheric(sps_ind)%spectag) THEN
        ier = 1
        WRITE(6,900) spectag
        RETURN
      END IF

!  'sps_ind' locates the index in the mixing ratio array:

      ndx_sps(sps_i,ch_i) = sps_ind

      IF(.NOT.atmospheric(sps_ind)%fwd_calc(band)) CYCLE sps

!  Now check if we loaded this Spectag already ..

      IF(pb > 0) THEN
        DO k = 1, pb
          IF(spectag == spectags(k)) CYCLE sps
        END DO
      END IF

      IF(pb == 9) THEN
        ier = 1
        PRINT *,'** Error in routine: read_EOS_db ..'
        PRINT *,'   This limited version does not support more then &
            &nine species..'
        RETURN
      END IF

      pb = pb + 1
      spectags(pb) = Spectag
      CALL read_eos_db(Spectag,sps_name,mdb_hdr(pb),mdb_rec,jpm,Ier)
      IF(ier /= 0) THEN
        PRINT *,'** Error in routine: read_EOS_db, Spectag:',spectag
        PRINT *,'   Called by routine: pfa_prep ..'
        RETURN
      END IF
!
! Now, build the beta arrays along the path
!
      do ptg_i = 1, no_tan_hts
!
        m = ndx_path(ptg_i)%total_number_of_elements
        ALLOCATE(values(m),t_power(m),dbeta_dw(m),dbeta_dn(m), &
                 dbeta_dnu(m),STAT=ier)
        if(ier /= 0) then
          PRINT *,'** Allocation error in routine: pfa_prep ..'
          PRINT *,'   IER =',ier
          RETURN
        endif
!
        do frq_i = 1, no_freqs
!
          Frq = freq_grid(frq_i)
          CALL get_beta_path (Spectag, Frq, m, z_path(ptg_i), t_path(ptg_i), &
       &                mdb_hdr(pb), mdb_rec, values, t_power, dbeta_dw, &
       &                dbeta_dn,dbeta_dnu)
!
          ALLOCATE(beta_path(sps_ind,frq_i,ptg_i)%values(m),    &
       &           beta_path(sps_ind,frq_i,ptg_i)%t_power(m),   &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dw(m),  &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dn(m),  &
       &           beta_path(sps_ind,frq_i,ptg_i)%dbeta_dnu(m), &
                   STAT = ier)
          if(ier /= 0) then
            PRINT *,'** Allocation error in routine: pfa_prep..'
            PRINT *,'   IER =',ier
            RETURN
          endif
!
          beta_path(sps_ind,frq_i,ptg_i)%values(1:m) = values(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%t_power(1:m) = t_power(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dw(1:m) = dbeta_dw(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dn(1:m) = dbeta_dn(1:m)
          beta_path(sps_ind,frq_i,ptg_i)%dbeta_dnu(1:m) = dbeta_dnu(1:m)
!
        end do

        DEALLOCATE(values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu,STAT=i)

      end do

    END DO sps

  END DO

!  Find the species index in the l2pc mixing ratio database:

  m = pb
  pb = -1
  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)
    band = (mch + 14) / 15
    no_sps = pfa_spectrum(band,1)%no_sps

    j4 = 4 * no_int_frqs(ch_i) + 1     ! Total # of Aitken's points

    IF(j4 > 1) THEN

      DO sps_i = 1, no_sps

        spectag = pfa_spectrum(band,sps_i)%sps_spectag
        spectags(sps_i) = spectag

!  Setup the spectrum record (structure), per band:

        j = 1
        k = spectag
        IF(k == 18999 .OR. k == 28964 .OR. k == 28965) j = -1

        IF(j >= 0 .AND. pb /= band) THEN

          k = -sps_i
          IF(sps_i > 1) k = sps_i

        END IF
!
!  Overwrite the pfa_spectrum data with the appropriate EOS Database
!
        j = 0
        k = -1
        DO WHILE(j < m .AND. k < 1)
          j = j + 1
          IF(mdb_hdr(j)%spectag == spectag) k = j
        END DO

        IF(k > 0) THEN
          j = mdb_hdr(k)%no_lines
          pfa_spectrum(band,sps_i)%no_lines = j
        END IF

      END DO

    END IF

    pb = band
    frq = freqs(mch)
    IF(frq < 1.0D0) THEN
      ier = 1
      WRITE(6,905) mch
      RETURN
    END IF

! Set up filter's response function

    q = 0.0
    df = Filter(q,mch,xlhs,xrhs,area,ier,InDir,primag,ld)
    IF(ier /= 0) GO TO 99

    df = (xrhs-xlhs)/(no_filt_pts-1)
    DO j = 1, no_filt_pts
      q = xlhs + (j - 1) * df
      f_grid_filter(j,ch_i) = frq + q
      filter_func(j,ch_i) = Filter(q)
    END DO

!  Normalize the filter's response array:

    CALL dsimps(filter_func(1:,ch_i),df,no_filt_pts,q)
    DO j = 1, no_filt_pts
      filter_func(j,ch_i) = filter_func(j,ch_i) / q
    END DO

!  If needed, Get Aitken's grid points:

    IF(j4 > 1) THEN
      df = (xrhs - xlhs) / (j4 - 1)
      DO j = 1, j4
        q = xlhs + (j - 1) * df
        f_grid(j,ch_i) = frq + q
      END DO
      f_grid(j4+1,ch_i) = frq                 ! *** DEBUG
    END IF

  END DO                         ! On ch_i

 99 IF(ier /= 0) THEN
      PRINT *,' ** Error in Pfa_Prep subroutine **'
      PRINT *,'    After calling subroutine: Get_Filter_Param'
      CALL errmsg(' ',ier)
    END IF

  900  FORMAT(' ** Error in Pfa_Prep subroutine **',/, &
      '    PFA species: ',i7.7,' not among L2PC species database !')
  905  FORMAT(' ** Error in Pfa_Prep subroutine **',/, &
      '    Inconsistant User Input.',/, &
      '    PFA Channel:',i3,' not among the non-PFA channels !')

  RETURN
END SUBROUTINE pfa_prep

!---------------------------------------------------------------------
!  This routine reads the EOS database and returns the structer holding
!  the required Spectag

SUBROUTINE read_eos_db(spectag,sps_name,mdb_hdr,mdb_rec,jpm,ier)

!  ===============================================================
!  Declaration of variables for sub-program: read_eos_db
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: spectag,jpm
Integer(i4), INTENT(OUT) :: ier

Character (LEN=8), intent (OUT) :: sps_name

type (eos_mdb_hdr), intent (OUT) :: mdb_hdr
type (eos_mdb_rec), intent (OUT) :: mdb_rec(*)
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, io, nl, du, iu, ii, jj, kk, no_lines

Integer(i4), save :: init = 0

Integer(i4) :: no_sps = 14
Integer(i4) :: spectags(14) = (/                          &
        32001, 34001, 18003, 44004, 63001, 27001, 51002,  &
        48004, 28001, 52006, 33001, 36001, 97001, 17001/)

Character (LEN=80) :: fhd, fdt, datdir
Character (LEN=8) :: names(14) = (/                                &
           'O2      ','O-18-O  ','H2O     ','N2O     ','HNO3    ', &
           'HCN     ','CLO     ','O3      ','CO      ','HOCL    ', &
           'HO2     ','HCL     ','BR-81-O ','OH      '/)

type (eos_mdb_hdr), save :: mdb_zero_hdr
type (eos_mdb_rec), save :: mdb_zero_rec

! Begin code:

  IF(init < 1) THEN
!
! Initialize mdb_zero_hdr:
!
    init = 5
    mdb_zero_hdr%spectag = 0
    mdb_zero_hdr%no_lines = 0
    mdb_zero_hdr%q_log(1) = 0.0
    mdb_zero_hdr%q_log(2) = 0.0
    mdb_zero_hdr%q_log(3) = 0.0
    DO ii = 1, max_no_lines
      mdb_zero_hdr%n(ii) = 0.0
      mdb_zero_hdr%w(ii) = 0.0
      mdb_zero_hdr%el(ii) = 0.0
      mdb_zero_hdr%n1(ii) = 0.0
      mdb_zero_hdr%n2(ii) = 0.0
      mdb_zero_hdr%v0(ii) = 0.0
      mdb_zero_hdr%ps(ii) = 0.0
      mdb_zero_hdr%delta(ii) = 0.0
      mdb_zero_hdr%gamma(ii) = 0.0
      mdb_zero_hdr%log_i(ii) = 0.0
      mdb_zero_hdr%no_f_grid(ii) = 0.0
      DO jj = 1, max_freq
        mdb_zero_hdr%x_grid(jj,ii) = 0.0
      END DO
    END DO

    DO ii = 1, max_zeta
      mdb_zero_hdr%zeta(ii) = 0.0
    END DO

    DO ii = 1, max_temp
      mdb_zero_hdr%log_temp(ii) = 0.0
    END DO
!
! Initialize mdb_zero_rec:
!
    DO ii = 1,max_zeta
      DO jj = 1, max_temp
        DO kk = 1, max_freq
          mdb_zero_rec%log_beta(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dw(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dn(ii,jj,kk) = 0.0
          mdb_zero_rec%dlog_beta_dnu0(ii,jj,kk) = 0.0
        END DO
      END DO
    END DO

  END IF

  i = 0
  j = 0
  ier = 0
  sps_name(1:) = ' '
  DO WHILE(j < no_sps .AND. i < 1)
    j = j + 1
    IF(spectags(j) == spectag) THEN
      i = j
      j = 23
    END IF
  END DO

  IF(i < 1) THEN
    io = -1
    GO TO 99
  END IF

  sps_name = names(i)

  datdir(1:) = ' '
  if(jpm < 1) then
    datdir= '/home/zvi/temp/minus/'
  else
    datdir= '/home/zvi/temp/plus/'
  endif

  fdt(1:) = ' '
  fhd(1:) = ' '
  j = LEN_TRIM(datdir)
  i = LEN_TRIM(sps_name)
  fdt = datdir(1:j)//sps_name(1:i)//'_eosmdb.dat'
  CALL strlwr(fdt)
  j = LEN_TRIM(fdt)
  i = INDEX(fdt,'.dat')
  fhd = fdt(1:i-1)//'.hdr'

  du = 43
  iu = du + 1
  CLOSE(iu,IOSTAT=io)
  inquire (iolength=k) mdb_hdr
  OPEN(iu,FILE=fhd,FORM='UNFORMATTED',STATUS='OLD',action='READ', &
       RECL=k,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) GO TO 99

  CLOSE(du,IOSTAT=io)
  inquire (iolength=j) mdb_rec(1)
  OPEN(du,FILE=fdt,FORM='UNFORMATTED',STATUS='OLD',action='READ', &
       RECL=j,ACCESS='DIRECT',IOSTAT=io)
  IF(io /= 0) GO TO 99

  mdb_hdr = mdb_zero_hdr                    ! Initialize mdb_hdr
  READ(iu,REC=1,IOSTAT=io) mdb_hdr
  IF(io /= 0) GO TO 99

  no_lines = mdb_hdr%no_lines

  DO nl = 1, no_lines

    mdb_rec(nl) = mdb_zero_rec              ! Initialize mdb_rec
    READ(du,REC=nl,IOSTAT=io) mdb_rec(nl)
    IF(io /= 0) GO TO 99

  END DO

  99   IF(io > 0) THEN
    ier = 1
    CALL errmsg(' ',io)
  END IF

  CLOSE(du,IOSTAT=i)
  CLOSE(iu,IOSTAT=i)

  RETURN
END SUBROUTINE read_eos_db
!
!---------------------------------------------------------------------
!
      SUBROUTINE GTI(MSG,TI)
!
      IMPLICIT NONE

      CHARACTER(LEN=*),INTENT(IN) :: MSG
      CHARACTER(LEN=*),INTENT(OUT) :: TI
!
      INTEGER :: I,J,L

      L=LEN(TI)
      TI(1:L)=' '
      I=LEN_TRIM(MSG)
      WRITE(*,900,ADVANCE='NO') MSG(1:I)
 900  FORMAT(A,': ')
      READ(5,'(A)',IOSTAT=J) TI
      J=LEN_TRIM(TI)
      IF(J > 1) CALL LEFTJ(TI)
!
      END SUBROUTINE GTI

!-----------------------------------------------------------------------

Subroutine Rad_Tran(Frq,N_lvls,h_tan,band,n_sps, sps_tbl, &
      &    ndx_path, z_path, h_path, t_path, phi_path, dHdz_path,     &
      &    earth_ref,beta_path, spsfunc_path, ref_corr, s_temp,brkpt, &
      &    no_ele, mid, ilo, ihi, cse, Rs, t_script, tau, Rad, Ier)
!
    Integer(i4), intent(in) :: N_LVLS, BAND, N_SPS

    Integer(i4), intent(in) :: SPS_TBL(Nsps,*)

    Integer(i4), intent(out) :: BRKPT, NO_ELE, MID, ILO, IHI, IER

    Real(r8), intent(in) :: FRQ, H_TAN, EARTH_REF, S_TEMP

    Real(r8), intent(in) :: REF_CORR(*)

    Type(path_beta), intent(in) :: BETA_PATH(:)   ! (Nsps)

    Type(path_index), intent(in)  :: NDX_PATH
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)

    Real(r8), intent(out) :: T_SCRIPT(*), TAU(*)
    Real(r8), intent(out) :: CSE, RS, RAD
!
    Integer(i4) :: Ngp1

    Real(r8) :: del_opacity(N2lvl), delta(N2lvl,Nsps)
!
! 'brkpt' is the index of the path break-point (when it change from
!         incoming ray to outgoing ray)
! 'no_ele' is the total number of entries in ?_path%values(1...no_ele)

    Ngp1 = Ng + 1
    brkpt = ndx_path%break_point_index
    no_ele = ndx_path%total_number_of_elements

    EarthX = .false.
    ht = h_tan
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
!  Compute the appropriate t_script & dt_scrpt_dnp along the integration
!  path of this tanget:
!
    CALL do_t_script(Ngp1, Frq, s_temp, brkpt, no_ele, t_path, &
   &                 mid, t_script)
!
    Call FAST_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path,beta_path, &
 &       dHdz_path,spsfunc_path,n_sps,N_lvls,sps_tbl(1:,band),Nlvl,    &
 &       ref_corr,delta,Ier)
    if (Ier /= 0) Return
!
! Initialize the tau & del_opacity arrays:
!
    tau(1:N2lvl) = 0.0
    del_opacity(1:N2lvl) = 0.0
!
    CALL FAST_ZOPACITY(sps_tbl(1:,band), n_sps, Ngp1, N2lvl, brkpt, &
   &                   no_ele, delta, del_opacity)
!
    Call Scrt_dn(t_script, N_lvls, cse, del_opacity, tau, Rad, mid, &
   &             ilo, ihi)
!
    Return
  End Subroutine RAD_TRAN

end Program l2_numder
