Program l2_numder
  use GL6P, only: NG
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
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use FILTER_SW_M, only: FILTER
  use DSIMPSON_MODULE, only: SIMPS
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  use REFRACTION_M
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use DO_T_SCRIPT_M, only: DO_T_SCRIPT
  use FAST_DELTA_M, only: FAST_DELTA
  use FAST_ZOPACITY_M, only: FAST_ZOPACITY
  use SCRT_DN_M, only: SCRT_DN
  use FREQ_AVG_M, only: FREQ_AVG
  use D_HUNT_M, only: HUNT

  implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: Npath = (Ng+1) * N2lvl

Integer(i4), PARAMETER :: mnf = 75

Integer(i4) :: p_indx(Nlvl), SPECT_ATMOS(Nsps), no_ptg_frq(Nptg), &
               no_coeffs_f(Nsps), no_phi_f(Nsps), SPECT_INDEX(Nsps), &
               no_spectro, pfa_ch(MAXPFACH), jpm, klo, kz, ma, frn

Integer(i4) :: i, j, k, kk, ht_i, no_t, mnz, no_geom, no_tan_hts, ld, &
               no_atmos, ch1, ch2, n_lvls, si, mfi, band, n_sps, Spectag, &
               no_freqs, ptg_i, frq_i, io, m, brkpt, no_ele, &
               mid, ilo, ihi, no_pfa_ch, n_obs, ier, no_filt_pts, no_phi_t

Type(path_index)  :: ndx_path(Nptg)
Type(path_vector) :: z_path(Nptg),t_path(Nptg),h_path(Nptg), &
                     n_path(Nptg),phi_path(Nptg),dhdz_path(Nptg), &
                     spsfunc_path(Nsps,Nptg),ptg_frq_grid(Nptg)

Type(path_derivative) :: dh_dt_path(Nptg)

Real(r8) :: href(Nlvl), zref(Nlvl), tan_hts(Nptg), dx, zco, &
            tan_hts_raw(Nptg), tan_press(Nptg), tan_temp(Nptg), &
            t_coeff(mxco,mnp), t_z_basis(mxco), t_phi_basis(mnp), &
            tan_dh_dt(Nlvl,mxco),z_grid(Nlvl), Frq, h_tan, Rad, &
            t_grid(Nlvl),h_grid(Nlvl),z_gnlv(400), e_rad, var, daz, &
            ref_corr(N2lvl,Nptg), t_script(N2lvl), tau(N2lvl)

Logical :: IS_F_LOG(Nsps)

real(r8) :: MR_F(mxco,mnp,Nsps)
real(r8) :: C_FREQ(Nch),freq_grid(mnf)
real(r8) :: F_BASIS(mxco,Nsps), PHI_BASIS_F(mnp,Nsps)
real(r8) :: FILTER_FUNC(maxfiltpts,maxpfach)
real(r8) :: F_GRID_FILTER(maxfiltpts,maxpfach)

!
Type(path_beta) :: beta_path(Nsps,mnf,Nptg)

Real(r8) :: Radiances(Nptg)
Real(r8) :: RadV(mnf),f_grid(mnf)
Real(r8) :: S_TEMP, H_OBS, EARTH_REF, q, r
!
Character :: Primag
Character (LEN=08) :: dName,Vname
Character (LEN=40) :: Ax, Fdata
Character (LEN=80) :: InDir, Fnd, Line

Logical :: do_mol

Type (l2pc_header_one) :: HEADER1
Type (atmos_comp)      :: ATMOSPHERIC(Nsps)
Type (geom_param)      :: GEOMETRIC(maxgeom)
Type (pfa_slab)        :: PFA_SPECTRUM(6,Nsps)
Type (spectro_param)   :: SPECTROSCOPIC(3*Nsps)

!  ----------------------
! Read the inputs from a file...
!
  ier = 0
  Fdata(1:) = ' '
  Fdata = '/user5/zvi/seez'     ! SUN, MLSGATE
! Fdata = '/home/zvi/seez'      ! HOME PC
  CLOSE(11,iostat=io)
  OPEN(11,file=Fdata,status='OLD',action='READ',iostat=io)
  if(io /= 0) goto 99

  ma = -1
  kz = -1
  dx = 0.0
  dName = ' '
  Vname = ' '
  do_mol = .false.
  Print *,'Enter name of variable to differentiate by. Choises are:'
  Call Gti(' dT, dMr, dw, dn & dNu',dName)
  if(dName < '!') goto 99
  Call StrLwr(dName)
  Print *

  Call Gti(' Enter step size',Ax)
  if(Ax < '!') goto 99
  read(Ax,*,iostat=io) dx
  if(io /= 0) goto 99
  if(abs(dx) < 1.0e-12) goto 99
  Print *

  Line(1:)=' '
  daz = -1.666667
  Line = 'Enter Zeta level of the coefficient to differentiate w.r.t.'
  i = LEN_TRIM(Line)
  Line = Line(1:i)//'(Default: -1.6667)'
  i = LEN_TRIM(Line)
  Print *,Line(1:i)
  Call Gti(' (Phi coefficient is taken to be the middle one)',Ax)
  if(Ax < '!') Ax = '-1.666667'
  read(Ax,*,iostat=i) daz
  if(i /= 0) daz = -1.666667
  Print *,Ax

  j = -1
  if(dName == 'dmr') j = 1
  if(dName == 'dw' .or. dName == 'dn' .or. dName == 'dnu') j = 2
  if(j > 0) then
    do_mol = .true.
    if(j == 2) zco = daz
    Call Gti(' Enter molecule name',Vname)
    if(Vname < '!') goto 99
    Call StrLwr(Vname)
    Print *,Vname
  endif

  p_indx(1:Nlvl) = 0
  no_phi_f(1:Nsps) = 0
  no_ptg_frq(1:Nptg) = 0
  pfa_ch(1:MAXPFACH) = 0
  no_coeffs_f(1:Nsps) = 0

  SPECT_INDEX(1:Nsps) = -1
  SPECT_ATMOS(1:Nsps) = -1

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
    if(ma < 1 .and. do_mol) then
      Ax = atmospheric(j)%NAME
      Call StrLwr(Ax)
      if(Ax == Vname) ma = j
    endif
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
    read(11,*,iostat=io) (spectroscopic(j)%PHI_BASIS(i),i=1,mfi+2)
    if(io /= 0) goto 99
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

  read(11,*,iostat=io) q, r, phi_tan   !, geoc_lat
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

  if(dName == 'dt') then
    Call Hunt(daz,t_z_basis,no_t,kz,i)
    IF(ABS(daz-t_z_basis(i)) < ABS(daz-t_z_basis(kz))) kz=i
    zco = t_z_basis(kz)
  endif

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

  DO m = 1, no_atmos
    kk = no_phi_f(m)
    ht_i = no_coeffs_f(m)
    read(11,*,iostat=io) (f_basis(i,m),i=1,ht_i)
    if(io /= 0) goto 99
    read(11,*,iostat=io) (phi_basis_f(i,m),i=1,kk)
    if(io /= 0) goto 99
    read(11,*,iostat=io) ((mr_f(i,j,m),j=1,kk),i=1,ht_i)
    if(io /= 0) goto 99
  END DO

  if(do_mol) then
    ht_i = no_coeffs_f(ma)
    Call Hunt(daz,f_basis(1:,ma),ht_i,kz,i)
    IF(ABS(daz-f_basis(i,ma)) < ABS(daz-f_basis(kz,ma))) kz=i
    zco = f_basis(kz,ma)
    if(dName == 'dmr') then
      k = (no_phi_f(ma)+1)/2
      r = 0.05 * abs(mr_f(kz,k,ma))
      q = sign(1.0_r8,dx) * max(1.0d-8,r)
      dx = q
      Print *,'** Modified Step Size:',Sngl(dx)
    endif
  endif

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) no_freqs
  if(io /= 0) goto 99

  read(11,'(A)',iostat=io) Ax
  if(io /= 0) goto 99

  read(11,*,iostat=io) (freq_grid(i),i=1,no_freqs)
  if(io /= 0) goto 99

  CLOSE(11,iostat=i)

  Print *
  r = abs(dx)
  Print *, 'VarName: ',dName
  Print *, 'Molecule: ',Vname
  Print *, 'Step Size: ',Sngl(r)
  Print *, 'Differentiated w.r.t. Coefficient #',kz
  Print *, 'Zeta of differentiated Coefficient:',Sngl(zco)
  Print *

! Get the selected integration grid pressures. Also, define the GL
! pressure grid:

  DO i = 1, n_lvls
    j = p_indx(i)
    z_grid(i) = z_gnlv(j)
  END DO
  z_grid(n_lvls+1) = z_grid(n_lvls)
!
  jpm = 0
  Radiances(1:Nptg) = 0.0

  if(dName == 'dw' .or. dName == 'dn' .or. dName == 'dnu') then
    jpm = 1
    zco = daz
    if(dx < 0.0) jpm = -1
  endif

  if(dName == 'dt') then
    k = (no_phi_t+1)/2
    var = t_coeff(kz,k)
    t_coeff(kz,k) = var + dx
  else if(dName == 'dmr') then
    k = (no_phi_f(ma)+1)/2
    var = mr_f(kz,k,ma)
    mr_f(kz,k,ma) = var + dx
  endif
!
  Primag = 'p'
  Call fwd_mdl_set_up(Primag,href,zref,e_rad,z_grid,t_grid,h_grid,    &
 &     ndx_path,n_lvls,geometric,no_geom,t_z_basis,t_coeff,no_t,      &
 &     atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,c_freq, &
 &     si,earth_ref,s_temp,h_obs,n_obs,tan_hts_raw,tan_hts,tan_press, &
 &     tan_temp,no_tan_hts,n_sps,no_pfa_ch,no_filt_pts,  &
 &     pfa_ch,pfa_spectrum,f_grid_filter,filter_func,ch1,ch2,InDir,ld,&
 &     band,no_phi_f,phi_basis_f,z_path,h_path,t_path,phi_path,n_path,&
 &     dhdz_path,dh_dt_path,Npath,no_phi_t,t_phi_basis,spectroscopic,&
 &     no_spectro,tan_dh_dt,spsfunc_path,no_ptg_frq,ptg_frq_grid,     &
 &     is_f_log,beta_path,jpm,ier)
  IF(ier /= 0) goto 99

  if(dName == 'dt') then
    k = (no_phi_t+1)/2
    t_coeff(kz,k) = var
  else if(dName == 'dmr') then
    k = (no_phi_f(ma)+1)/2
    mr_f(kz,k,ma) = var
  endif
!
! Compute the refraction correction scaling matrix:
!
  CALL refraction_correction(no_tan_hts, tan_hts, h_path, n_path, &
 &                           ndx_path, E_rad, ref_corr)
!
! Now, Compute the radiances alone:
!
  Call Hunt(daz,tan_press,no_tan_hts,klo,i)
  IF(ABS(daz-tan_press(i)) < ABS(daz-tan_press(klo))) klo=i

! do ptg_i = 1, no_tan_hts-1
  do ptg_i = klo, klo                        ! ** DEBUG
!
    k = ptg_i
    h_tan = tan_hts(k)
    kk = no_ptg_frq(k)
!
    RadV(1:mnf) = 0.0
    f_grid(1:kk) = ptg_frq_grid(k)%values(1:kk)

    Print *
    frn = 10                                   ! ** DEBUG
    Print *,'** Frequency used:',f_grid(frn)   ! ** DEBUG
!
!   do frq_i = 1, kk
    do frq_i = frn, frn                        ! ** DEBUG
!
      Frq = f_grid(frq_i)
!
      Call Rad_Tran(Frq, N_lvls, h_tan, n_sps, ndx_path(k),  &
     &    z_path(k), h_path(k), t_path(k), phi_path(k), dHdz_path(k), &
     &    earth_ref, beta_path(1:,frq_i,k), spsfunc_path(1:,k),       &
     &    ref_corr(1:,k), s_temp, brkpt, no_ele, mid, ilo, ihi,       &
     &    t_script, tau, Rad, Ier)
      IF(ier /= 0) goto 99
!
      RadV(frq_i) = Rad
!
    end do
!
! Frequency Average the radiances with the appropriate filter shapes
!
!   i = 1
!   Call Freq_Avg(f_grid,F_grid_filter(1:,i),Filter_func(1:,i), &
!  &              RadV,kk,maxfiltpts,Rad,Ier)
!   IF(ier /= 0) goto 99
!   Radiances(ptg_i) = Rad
    Radiances(ptg_i) = RadV(frn)          ! ** DEBUG

  end do
!
  j = no_tan_hts
  Radiances(j) = Radiances(j-1)

    kk = j - si + 1
    write(*,903) char(92),kk
    write(*,904) (Radiances(k),k=si,j)
903 format('avg_pfa_rad',a1,i2.2)
904 format(4(2x,1pg17.10))

 99  CLOSE(11,iostat=i)
     if(io /= 0) Call ErrMsg(' ',io)

     Stop

!----------------------------------------------------------------
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
           atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,freq, &
           si,earth_ref,s_temp,h_obs,n_obs,tan_hts_raw,tan_hts,tan_press,&
           tan_temp,no_tan_hts,n_sps,no_pfa_ch,no_filt_pts, &
           pfa_ch,pfa_spectrum,f_grid_fltr,fltr_func,ch1,ch2,InDir,ld,   &
           band,no_phi_f,phi_basis_f,z_path,h_path,t_path,phi_path,n_path,&
           dhdz_path,dh_dt_path,npath,no_phi_t,t_phi_basis,   &
           spectroscopic,no_spectro,tan_dh_dt,spsfunc_path,no_ptg_frq,  &
           ptg_frq_grid,is_f_log,beta_path,jpm,ier)

!  ===============================================================
!  Declaration of variables for sub-program: fwd_mdl_set_up
!  ===============================================================

Integer(i4), PARAMETER :: ngt = (Ng+1) * Nlvl

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_atmos, no_t, ch1, ch2, jpm,  &
             no_phi_t, no_coeffs_f(:), no_phi_f(:), band, n_sps, &
             no_spectro, npath, si, no_pfa_ch, no_filt_pts, pfa_ch(:)
!
Integer(i4), INTENT(IN OUT) :: no_tan_hts, no_geom, ld

Integer(i4), INTENT(OUT) :: n_obs,no_ptg_frq(:),ier

!
Real(r8), INTENT(IN) :: z_grid(:), href(:), zref(:), &
          t_phi_basis(:), t_coeff(:,:), t_z_basis(:), &
          mr_f(:,:,:), phi_basis_f(:,:), f_basis(:,:)
!
Real(r8), INTENT(OUT) :: tan_hts(:), tan_dh_dt(:,:)
Real(r8), INTENT(IN OUT) :: tan_hts_raw(:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:)
!

Type(path_index) , INTENT(OUT) :: ndx_path(:)
Type(path_vector), INTENT(OUT) :: z_path(:),t_path(:),h_path(:), n_path(:), &
             phi_path(:), dhdz_path(:), spsfunc_path(:,:), ptg_frq_grid(:)

Type(path_derivative), INTENT(OUT) :: dh_dt_path(:)

Real(r8), INTENT(OUT) :: e_rad,earth_ref,s_temp,h_obs,freq(:), &
          tan_temp(:),t_grid(:),h_grid(:),tan_press(:), &
          fltr_func(:,:),f_grid_fltr(:,:)

Logical, INTENT(IN) :: is_f_log(*)
!
Character (LEN=*), INTENT(IN) :: InDir
Character (LEN=*), INTENT(IN) :: primag
!  ----------------------
!  Local variables:
!  ----------------

Logical :: wet
Integer(i4) :: i, j, k, jp, io, ich, sps_i, geo_i, geo_j, ih2o, kk, &
               no_geom_new, gl_count

Real(r4) :: def_val(maxgeom) = (/                                    &
            0.0, 0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.0, 6972.0, 6372.0, &
            0.05, 2.735, 1.0 /)

Type(path_vector) :: h2o_path(Nptg)

Real(r8) :: q,r,p,g,rp,beta_inc,zeta,phi,geoc_lat

Character (LEN=80) :: Fnd
Character (LEN=8) :: geomname(maxgeom) = (/                         &
         &  'ELEV_183','ELEV_205','AZIM_183','AZIM_205','AZIM_REF', &
         &  'ROLL    ','PITCH   ','YAW     ','GEOCSRAD','GEOCERAD', &
         &  'EARTHREF','SPACE_T ','LOOK_DIR'/)

!  PFA variables:

type (geom_param),    intent(inout) :: GEOMETRIC(*)
type (atmos_comp),    intent(inout) :: ATMOSPHERIC(*)
type (spectro_param), intent(inout) :: SPECTROSCOPIC(*)
type (pfa_slab), intent(inout)      :: PFA_SPECTRUM(6,*)
!
! Set up the default geometric linearization values of earth radius,earth
! emission,background space temperature and satellite position.
! Compute the default geocentric earth radius

  ier = 0

  beta_inc = 98.0D0                ! Orbit inclination angle for EOS-MLS
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  Call geoc_geod_conv('d2c',beta_inc,phi_tan,geoc_lat,rp)
!
! Convert the t_phi_basis to radians and add phi_tan:
!
  do j = 1, no_phi_t
    t_phi_basis(j) = deg2rad * t_phi_basis(j) + phi_tan
  end do
!
! Convert the phi_basis_f to radians and add phi_tan:
!
    do j = 1, n_sps
      i = no_phi_f(j)
      phi_basis_f(1:i,j) = deg2rad * phi_basis_f(1:i,j) + phi_tan
    end do
!
! Convert the spectroscopic to radians and add phi_tan:
!
  do j = 1, no_spectro
    i = spectroscopic(j)%NO_PHI_VALUES
    spectroscopic(j)%PHI_BASIS(1:i) = &
   &        deg2rad * spectroscopic(j)%PHI_BASIS(1:i) + phi_tan
  end do

! Call the grids program to get the preselected integration heights and
! pointings. Also transfer the pointings expressed in pressure units to
! heights.  Also, set up geophysical parameters (for derivatives):
! [The following routine replaces the older code using the two routines:
!  grids() and geo_basis() ]

! (z_grid,t_grid,h_grid) are the arrays of the preselected integration grid.

  CALL geo_grids(z_grid,t_grid,h_grid,tan_dh_dt,n_lvls,t_z_basis,    &
       t_coeff,si,tan_hts,tan_hts_raw,tan_press,tan_temp,no_tan_hts, &
       z_path,h_path,t_path,phi_path,dhdz_path,dh_dt_path,Npath,g,href, &
       zref,geoc_lat,no_t,no_phi_t,t_phi_basis,ndx_path,ier)
  IF(ier /= 0) Return

  ps = -1.0D0

  e_rad = rp
  def_val(10) = e_rad

! Default earth emission

  earth_ref = def_val(11)

! Default satellite position

  h_obs = def_val(9)

! Default space temperature

  s_temp = def_val(12)
!
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
!
  DO ich = 1, nch
    freq(ich) = 0.0D0
  END DO

  DO ich = ch1, ch2
    CALL radiometry(ich,q,r,p,kk)       ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'p') freq(ich) = q     ! DEBUG, Added Jan/23/2000, Z.S
    IF(primag == 'i') freq(ich) = r     ! DEBUG, Added Jan/23/2000, Z.S
  END DO
!
! Create the specie function along the path for all species
!
  DO k = 1, no_tan_hts
    gl_count = ndx_path(k)%total_number_of_elements
    do j = 1, n_sps
      jp = no_phi_f(j)
      kk = no_coeffs_f(j)
      DEALLOCATE(spsfunc_path(j,k)%values,STAT=i)
      ALLOCATE(spsfunc_path(j,k)%values(gl_count),STAT=i)
      IF(i /= 0) THEN
        ier = i
        PRINT *,'** Error: ALLOCATION error for spsfunc_path ..'
        PRINT *,'   STAT =',ier
        Return
      ENDIF
      do i = 1, gl_count
        zeta = z_path(k)%values(i)
        phi = phi_path(k)%values(i)
        if (is_f_log(j)) then
          Call TWO_D_POLATE(f_basis(1:,j), LOG(mr_f(1:kk,1:jp,j)), &
         &                  kk, phi_basis_f(1:,j), jp, zeta, phi, r)
          q = exp(r)
        else
          Call TWO_D_POLATE(f_basis(1:,j), mr_f(1:kk,1:jp,j), kk,  &
         &                  phi_basis_f(1:,j), jp, zeta, phi, q)
        endif
        spsfunc_path(j,k)%values(i) = q
      end do
    end do
  END DO
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
  Open(32,file=Fnd,action='READ',iostat=io)
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
      Return
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
        Return
      ENDIF
      no_ptg_frq(k) = jp
      ptg_frq_grid(k)%values(1:jp) = ptg_frq_grid(kk)%values(1:jp)
    end do
  endif
!
 44 Close(32,iostat=i)
  if(io > 0) then
    ier = io
    Return
  endif
!
! Set n_obs to be: N_lvls always.
!   (Changed, Aug/6/96 Z.Shippony & W.G.Read)

  n_obs = n_lvls

!  *** Special PC code - check if number of Keys will exceed maximum:

  kk = 3                ! x_star key + i_star key + Ptan63 Key
  DO geo_i = 1, no_geom
    IF(geometric(geo_i)%der_calc(band)) kk = kk + 1
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
    Return
  END IF

!  *** End Special PC code

! Compute the relative refractive index minus one.
! Get the water mixing ratio function

  sps_i = 1
  DO WHILE (atmospheric(sps_i)%name /= 'H2O' .AND. sps_i <= no_atmos)
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
    CALL pfa_prep(atmospheric,band,no_atmos,no_pfa_ch,no_filt_pts,pfa_ch, &
         pfa_spectrum,f_grid_fltr,freq,fltr_func,no_tan_hts,  &
         ndx_path,no_ptg_frq,ptg_frq_grid,z_path,t_path,    &
         beta_path,InDir,ld,primag,jpm,ier)
    IF(ier /= 0) Return
  ENDIF

  Return
END SUBROUTINE fwd_mdl_set_up

!-----------------------------------------------------------------------

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

!----------------------------------------------------------------

SUBROUTINE pfa_prep(atmospheric,band,no_atmos,no_pfa_ch,no_filt_pts,    &
           pfa_ch,pfa_spectrum,f_grid_filter,freqs,filter_func,  &
           no_tan_hts,ndx_path,no_ptg_frq,ptg_frq_grid,z_path, &
           t_path,beta_path,InDir,ld,primag,jpm,ier)

!  ===============================================================
!  Declaration of variables for sub-program: pfa_prep
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_atmos, no_pfa_ch, no_filt_pts, pfa_ch(*), &
                           no_tan_hts, band, jpm, no_ptg_frq(*)

Integer(i4), INTENT(OUT) :: ier, ld

Real(r8), INTENT(IN) :: freqs(*)

Type(path_index), INTENT(IN) :: ndx_path(*)
Type(path_vector), INTENT(IN) :: z_path(*), t_path(*), ptg_frq_grid(*)

Type (atmos_comp), INTENT(IN) :: ATMOSPHERIC(*)

Type (pfa_slab), INTENT(INOUT) :: PFA_SPECTRUM(6,*)

Real(r8), INTENT(OUT) :: filter_func(:,:), f_grid_filter(:,:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:)  ! (sps_i,frq_i,ptg_i)

Character (LEN=*), INTENT(IN) :: InDir, primag

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, m, ch_i, sps_i, pb, no_sps, spectag, &
               mch, ptg_i, frq_i, spectags(MAXSPS)

Real(r8) :: xlhs, xrhs, df, q, area, frq

Character (LEN=8) :: sps_name

! Type (eos_mdb_hdr) :: MDB_HDR(MAXSPS)
! Type (eos_mdb_rec) :: MDB_REC(MAX_NO_LINES)
Type (eos_mdb_hdr) :: MDB_HDR(02)      ! ** DEBUG
Type (eos_mdb_rec) :: MDB_REC(02)      ! ** DEBUG

Real(r8), DIMENSION(:), ALLOCATABLE :: values, t_power, dbeta_dw, &
                                       dbeta_dn, dbeta_dnu
! Begin code:

  ier = 0

!  Read in the EOS Beta database for all species involved

  pb = 0
  spectags(1:MAXSPS) = -1
  no_sps = pfa_spectrum(band,1)%no_sps

sps:DO sps_i = 1, no_sps

    Spectag = pfa_spectrum(band,sps_i)%sps_spectag

    j = 1
    DO WHILE(j < no_atmos .AND. spectag /= atmospheric(j)%spectag)
      j = j + 1
    END DO

!  Check if specie is in the database

    IF(spectag /= atmospheric(j)%spectag) THEN
      ier = 1
      WRITE(6,900) spectag
      Return
    END IF

    IF(.NOT.atmospheric(j)%fwd_calc(band)) CYCLE sps

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
      Return
    END IF

    pb = pb + 1
    spectags(pb) = Spectag
    CALL read_eos_db(Spectag,sps_name,mdb_hdr(pb),mdb_rec,jpm,Ier)
    IF(ier /= 0) THEN
      PRINT *,'** Error in routine: read_EOS_db, Spectag:',spectag
      PRINT *,'   Called by routine: pfa_prep ..'
      Return
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
        Return
      endif
!
      k = ptg_i
      do frq_i = 1, no_ptg_frq(k)
!
        DEALLOCATE(beta_path(j,frq_i,ptg_i)%values,    &
     &             beta_path(j,frq_i,ptg_i)%t_power,   &
     &             beta_path(j,frq_i,ptg_i)%dbeta_dw,  &
     &             beta_path(j,frq_i,ptg_i)%dbeta_dn,  &
     &             beta_path(j,frq_i,ptg_i)%dbeta_dnu, &
     &             STAT=i)
!
        Frq = ptg_frq_grid(k)%values(frq_i)
        CALL get_beta_path (Spectag, Frq, m, z_path(k), t_path(k),     &
     &                mdb_hdr(pb), mdb_rec, values, t_power, dbeta_dw, &
     &                dbeta_dn,dbeta_dnu)
!
        ALLOCATE(beta_path(j,frq_i,ptg_i)%values(m),    &
     &           beta_path(j,frq_i,ptg_i)%t_power(m),   &
     &           beta_path(j,frq_i,ptg_i)%dbeta_dw(m),  &
     &           beta_path(j,frq_i,ptg_i)%dbeta_dn(m),  &
     &           beta_path(j,frq_i,ptg_i)%dbeta_dnu(m), &
     &           STAT = ier)
        if(ier /= 0) then
          PRINT *,'** Allocation error in routine: pfa_prep ..'
          PRINT *,'   IER =',ier
          Return
        endif
!
        beta_path(j,frq_i,ptg_i)%values(1:m) = values(1:m)
        beta_path(j,frq_i,ptg_i)%t_power(1:m) = t_power(1:m)
        beta_path(j,frq_i,ptg_i)%dbeta_dw(1:m) = dbeta_dw(1:m)
        beta_path(j,frq_i,ptg_i)%dbeta_dn(1:m) = dbeta_dn(1:m)
        beta_path(j,frq_i,ptg_i)%dbeta_dnu(1:m) = dbeta_dnu(1:m)
!
      end do

      DEALLOCATE(values,t_power,dbeta_dw,dbeta_dn,dbeta_dnu,STAT=i)

    end do

  END DO sps

! Find the species index in the l2pc mixing ratio database:

  m = pb
  pb = -1
  no_sps = pfa_spectrum(band,1)%no_sps

  DO ch_i = 1, no_pfa_ch

    mch = pfa_ch(ch_i)

    pb = band
    frq = freqs(mch)
    IF(frq < 1.0D0) THEN
      ier = 1
      WRITE(6,905) mch
      Return
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

    CALL Simps(filter_func(1:,ch_i),df,no_filt_pts,q)
    DO j = 1, no_filt_pts
      filter_func(j,ch_i) = filter_func(j,ch_i) / q
    END DO

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

  Return

END SUBROUTINE pfa_prep

!---------------------------------------------------------------------
!  This routine reads the EOS database and Returns the structer holding
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

Integer(i4), SAVE :: init = 0

Integer(i4) :: no_sps = 14
Integer(i4) :: spectags(14) = (/                          &
        32001, 34001, 18003, 44004, 63001, 27001, 51002,  &
        48004, 28001, 52006, 33001, 36001, 97001, 17001/)

Character (LEN=80) :: fhd, fdt, datdir
Character (LEN=8) :: names(14) = (/                                &
           'O2      ','O-18-O  ','H2O     ','N2O     ','HNO3    ', &
           'HCN     ','CLO     ','O3      ','CO      ','HOCL    ', &
           'HO2     ','HCL     ','BR-81-O ','OH      '/)

type (eos_mdb_hdr), SAVE :: mdb_zero_hdr
type (eos_mdb_rec), SAVE :: mdb_zero_rec

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
  Call StrLwr(sps_name)

  fdt(1:) = ' '
! fdt = '/home/zvi/temp/'              ! HOME PC
  fdt = '/user5/zvi/linux/temp/'       ! MLSGATE, VANPC
! fdt = '/user5/zvi/temp/'             ! SUN, SGI

  i = LEN_TRIM(fdt)
  if(jpm < 0) then
    datdir = fdt(1:i)//'/minus/'
!   datdir = '/user5/zvi/linux/temp/minus/zero_ps/'  ! ** DEBUG
  else if(jpm > 0) then
    datdir = fdt(1:i)//'/plus/'
!   datdir = '/user5/zvi/linux/temp/plus/zero_ps/'   ! ** DEBUG
  else
!   datdir = '/home/zvi/data/'              ! HOME PC
!   datdir = '/user5/zvi/zvi/eos/data/'     ! SUN, SGI
    datdir = '/user5/zvi/linux/MLS/data/'   ! MLSGATE, VANPC
!   datdir = '/user5/zvi/linux/MLS/data/zero_ps/'    ! ** DEBUG
  endif

  fdt(1:) = ' '
  fhd(1:) = ' '
  j = LEN_TRIM(datdir)
  i = LEN_TRIM(sps_name)
  fdt = datdir(1:j)//sps_name(1:i)//'_eosmdb.dat'
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

  Return
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

!----------------------------------------------------------------------
! This is the radiative transfer model, radiances only !

    Subroutine Rad_Tran(Frq,N_lvls,h_tan,n_sps, &
      &    ndx_path, z_path, h_path, t_path, phi_path, dHdz_path,     &
      &    earth_ref,beta_path, spsfunc_path, ref_corr, s_temp,brkpt, &
      &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
!
    Integer(i4), intent(in) :: N_LVLS, N_SPS

    Integer(i4), intent(out) :: BRKPT, NO_ELE, MID, ILO, IHI, IER

    Real(r8), intent(in) :: FRQ, H_TAN, EARTH_REF, S_TEMP

    Real(r8), intent(in) :: REF_CORR(*)

    Type(path_beta), intent(in) :: BETA_PATH(:)   ! (Nsps)

    Type(path_index), intent(in)  :: NDX_PATH
    Type(path_vector), intent(in) :: Z_PATH, T_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)

    Real(r8), intent(out) :: T_SCRIPT(*), TAU(*)
    Real(r8), intent(out) :: RAD
!
    Integer(i4) :: Ngp1

    Real(r8) :: CSE, RS
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
 &       dHdz_path,spsfunc_path,n_sps,N_lvls,Nlvl, &
 &       ref_corr,delta,Ier)
    if (Ier /= 0) Return
!
! Initialize the tau & del_opacity arrays:
!
    tau(1:N2lvl) = 0.0
    del_opacity(1:N2lvl) = 0.0
!
    CALL FAST_ZOPACITY(n_sps, Ngp1, N2lvl, brkpt, &
   &                   no_ele, delta, del_opacity)
!
    Call Scrt_dn(t_script, N_lvls, cse, del_opacity, tau, Rad, mid, &
   &             ilo, ihi)
!
    Return
  End Subroutine RAD_TRAN

End Program l2_numder
