Program l2_numder
  use GL6P, only: NG
  use MLSCommon, only: I4, R8
  use L2_TEST_STRUCTURES_M
  use STRINGS, only: STRLWR
  use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: LIMB_PRESS
  use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, MNP => max_no_phi, &
                     MNM => max_no_mmaf
  use ELLIPSE, only: PHI_TAN, ROC
  use L2_LOAD_M, only: L2_LOAD
  use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
  use REFRACTION_M, only: REFRACTION_CORRECTION
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                             PATH_DERIVATIVE
  use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
  use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
  use GET_BETA_PATH_M, only: GET_BETA_PATH
  use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
  use RAD_TRAN_M, only: RAD_TRAN
  use FOV_CONVOLVE_M, only: FOV_CONVOLVE
  use FREQ_AVG_M, only: FREQ_AVG
  use D_HUNT_M, only: HUNT          ! ** DEBUG
  use D_CSPLINE_M, only: CSPLINE

  Implicit NONE
!---------------------------------------------------------------------------
!
Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl

Type(fwd_mdl_info) :: FMI
Type(fwd_mdl_config) :: FMC
Type(temporary_fwd_mdl_info) :: T_FMI

Integer(i4) :: SPECT_ATMOS(Nsps)

Integer(i4) :: ma, line_no, kz, ier, ch1, ch2, no_pfa_ch, pfa_ch(2)
Integer(i4) :: i, j, k, kk, jj, ht_i, no_t, no_tan_hts, ch, Spectag, &
               prev_npf, n_obs, si, band, n_sps,  ptg_i, &
               frq_i, io, l, brkpt, no_ele, mid, m, ilo, ihi, no_phi_t, &
               gl_count, mmaf

Type(path_index)  :: ndx_path(Nptg,mnm)
Type(path_vector) :: z_path(Nptg,mnm),t_path(Nptg,mnm),h_path(Nptg,mnm),  &
                     dhdz_path(Nptg,mnm), spsfunc_path(Nsps,Nptg,mnm),    &
                     n_path(Nptg,mnm),phi_path(Nptg,mnm)

Type(path_derivative) :: dh_dt_path(Nptg,mnm)

Real(r8) :: thbs(10),elev_offset
Real(r8) :: t_script(N2lvl),ref_corr(N2lvl,Nptg),tau(N2lvl), &
            tan_dh_dt(Nlvl,mnm,mxco)

Real(r8) :: dx,var,zco,daz
Real(r8) :: dx_dt(Nptg,mxco), d2x_dxdt(Nptg,mxco)

Real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
Real(r8) :: dh_dt_glgrid(ngt,mnm,mxco), dhdz_glgrid(ngt,mnp)

Real(r8) :: ptg_angles(Nptg,mnm), center_angle
Real(r8) :: tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)

Logical :: do_mol, do_spec

Real(r8) :: T_PHI_BASIS_COPY(mnp)
Real(r8) :: F_PHI_BASIS_COPY(mnp,Nsps)

Real(r8), DIMENSION(:), ALLOCATABLE :: RadV, F_grid
Real(r8), DIMENSION(:,:), ALLOCATABLE :: S_PHI_BASIS_COPY

Real(r8) :: I_STAR_ALL(Nptg)
!
Type(path_beta), DIMENSION(:,:), POINTER :: beta_path    ! (sps_i,frq_i)

Real(r8) :: Radiances(Nptg)
Real(r8) :: e_rad, zeta, Frq, h_tan, Rad, geoc_lat, q, r
!
Character (LEN=16) :: dName, Vname
Character (LEN=80) :: Line
Character (LEN=40) :: Ax

!  ----------------------

  Line(1:) = ' '
! Line = '/home/zvi/seez'      ! HOME PC, JPL PC
  Line = '/user5/zvi/seez'     ! MLSGATE, SUN
  FMC%Z = Line
!
! Load all needed data via l2_load routine:
!
  Call L2_LOAD(FMC, FMI, T_FMI, Ier)
  if(ier /= 0) goto 99

  ch1 = FMC%Channels_range(1)
  ch2 = FMC%Channels_range(2)
  no_pfa_ch = min(2,ch2-ch1+1)
  do i = 1, no_pfa_ch
    pfa_ch(i) = ch1 + i - 1
  end do

  elev_offset = 0.0                         ! Zero elev_offset in any case
  SPECT_ATMOS(1:Nsps) = -1
!
! Set n_obs to be: N_lvls always.
!   (Changed, Aug/6/96 Z.Shippony & W.G.Read)

  band = FMI%band
  n_obs = FMC%n_lvls

  no_t = T_FMI%no_t
  no_phi_t = T_FMI%no_phi_t

  T_PHI_BASIS_COPY(1:no_phi_t) = T_FMI%t_phi_basis(1:no_phi_t)
!
  n_sps = FMI%n_sps
  DO m = 1, n_sps
    kk = T_FMI%no_phi_f(m)
    F_PHI_BASIS_COPY(1:kk,m) = T_FMI%f_phi_basis(1:kk,m)
  END DO
!
! Create spect_atmos array:
!
  do k = 1, n_sps
    Spectag = T_FMI%atmospheric(k)%Spectag
    do i = 1, FMI%no_spectro
      if(Spectag == FMI%spectroscopic(i)%Spectag) then
        spect_atmos(k) = i
        EXIT
      endif
    end do
  end do
!
  kk = mxco
  k = FMI%mfi + 2

  DEALLOCATE(S_PHI_BASIS_COPY,STAT=i)
  ALLOCATE(S_PHI_BASIS_COPY(k,FMI%no_spectro),STAT=io)
  IF(io /= 0) then
    Print *,'** ALLOCATE Error: S_PHI_BASIS_COPY, STAT =',io
    goto 99
  endif

  do j = 1, FMI%no_spectro
    S_PHI_BASIS_COPY(1:k,j) = FMI%spectroscopic(j)%PHI_BASIS(1:k)
  end do
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  mmaf = 3                     ! Do only this mmaf (middle phi)
  phi_tan = FMC%phi_tan_mmaf(mmaf)
  Call geoc_geod_conv(T_FMI%beta_inc,phi_tan,geoc_lat,E_rad)
!
  ma = -1
  kz = -1
  dx = 0.0
  dName = ' '
  Vname = ' '
  daz = -1.666667
  do_mol = .false.
  do_spec = .false.
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
    if(j == 2) then
      zco = daz
      do_spec = .true.
    else
      do_mol = .true.
    endif
    Call Gti(' Enter molecule name',Vname)
    if(Vname < '!') goto 99
    Call StrLwr(Vname)
    Print *,Vname
  endif

  if(dName == 'dt') then
    Call Hunt(daz,T_FMI%t_zeta_basis,no_t,kz,i)
    IF(ABS(daz-T_FMI%t_zeta_basis(i)) < ABS(daz-T_FMI%t_zeta_basis(kz))) kz=i
    zco = T_FMI%t_zeta_basis(kz)
  endif

  if(do_mol) then
    ht_i = T_FMI%no_coeffs_f(ma)
    Call Hunt(daz,T_FMI%f_zeta_basis(1:,ma),ht_i,kz,i)
    IF(ABS(daz-T_FMI%f_zeta_basis(i,ma)) <  &
   &   ABS(daz-T_FMI%f_zeta_basis(kz,ma))) kz=i
    zco = T_FMI%f_zeta_basis(kz,ma)
    k = (T_FMI%no_phi_f(ma)+1)/2
    r = 0.05 * abs(T_FMI%mr_f(kz,k,ma))
    q = sign(1.0_r8,dx) * max(1.0d-8,r)
    dx = q
    Print *,'** Modified Step Size:',Sngl(dx)
  endif

  if(do_spec) then
    line_no = 1
    DO j = 1, n_sps
      Ax = T_FMI%atmospheric(j)%NAME
      Call StrLwr(Ax)
      if(Ax == Vname) then
        if(dName == 'dw') then
          var = FMI%pfa_spectrum(j)%SPS_W(line_no)
          r = 0.1 * abs(var)
          q = sign(1.0_r8,dx) * max(1.0d-8,r)
          dx = q
          Print *,'** Modified Step Size:',Sngl(dx)
          FMI%pfa_spectrum(j)%SPS_W(line_no) = var + dx
        else if(dName == 'dn') then
          var = FMI%pfa_spectrum(j)%SPS_N(line_no)
          r = 0.1 * abs(var)
          q = sign(1.0_r8,dx) * max(1.0d-8,r)
          dx = q
          Print *,'** Modified Step Size:',Sngl(dx)
          FMI%pfa_spectrum(j)%SPS_N(line_no) = var + dx
        else if(dName == 'dnu') then
          var = FMI%pfa_spectrum(j)%SPS_V0(line_no)
          Print *,'** Modified Step Size:',Sngl(dx)
          FMI%pfa_spectrum(j)%SPS_V0(line_no) = var + dx
        endif
        EXIT
      endif
    END DO
  endif
!
  Print *
  Print *, 'VarName: ',dName
  Print *, 'Molecule: ',Vname
  Print *, 'Step Size: ',Sngl(dx)
  Print *, 'Differentiated w.r.t. Coefficient #',kz
  Print *, 'Zeta of differentiated Coefficient:',Sngl(zco)
  Print *
!
! Perturbe the middle Phi coefficient for the 'zco' pressure level
!
  if(dName == 'dt') then
    k = (FMC%no_mmaf+1)/2
    var = T_FMI%t_coeff(kz,k)
    T_FMI%t_coeff(kz,k) = var + dx
  else if(dName == 'dmr') then
    k = (T_FMI%no_phi_f(ma)+1)/2
    var = T_FMI%mr_f(kz,k,ma)
    T_FMI%mr_f(kz,k,ma) = var + dx
  endif
!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
  mmaf = 3                       ! Do only this mmaf (The middle phi)
  phi_tan = FMC%phi_tan_mmaf(mmaf)
  Call geoc_geod_conv(T_FMI%beta_inc,phi_tan,geoc_lat,E_rad)
!
! Compute the hydrostatic_model on the GL-Grid for all mmaf(s):
!
  thbs(1:) = 0.0
  si = FMI%Surface_index
  thbs(1:si-1) = FMI%Tan_hts_below_surface(1:si-1)
  Call hydrostatic_model(si,FMC%N_lvls,T_FMI%no_t,FMC%no_mmaf,FMC%t_indx, &
       FMC%no_tan_hts,geoc_lat,T_FMI%Href,T_FMI%Zref,FMI%z_grid,thbs, &
       T_FMI%t_zeta_basis, T_FMI%t_coeff, z_glgrid, h_glgrid, t_glgrid, &
       dhdz_glgrid,dh_dt_glgrid,FMI%tan_press,tan_hts,tan_temp,tan_dh_dt, &
       gl_count, Ier)
  IF(ier /= 0) goto 99

  Zeta = zco
  no_tan_hts = FMC%no_tan_hts
  Call Hunt(Zeta,FMI%tan_press,no_tan_hts,jj,i)
  IF(ABS(Zeta-FMI%tan_press(i)) < ABS(Zeta-FMI%tan_press(jj))) jj = i
!
! Compute all path entities for all mmafs and tanget pointings
!
  Call comp_path_entities(FMC%n_lvls,T_FMI%no_t,gl_count,ndx_path, &
       z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,        &
       T_FMI%atmospheric,T_FMI%f_zeta_basis,T_FMI%mr_f,            &
       T_FMI%no_coeffs_f,tan_hts,no_tan_hts,FMI%n_sps,             &
       T_FMI%no_phi_f,T_FMI%f_phi_basis,z_path,h_path,t_path,phi_path,&
       n_path,dhdz_path,dh_dt_path,T_FMI%no_phi_t,T_FMI%t_phi_basis,  &
       spsfunc_path,T_FMI%is_f_log,FMC%no_mmaf,FMC%phi_tan_mmaf,Ier)
  IF(ier /= 0) goto 99
!
! **********************  MAIN Mmaf Loop *******************

  l = mmaf                    ! ** Do this mmaf only
!
  phi_tan = FMC%phi_tan_mmaf(l)
!
  T_FMI%t_phi_basis(1:no_phi_t) = t_phi_basis_copy(1:no_phi_t) + phi_tan

  DO j = 1, n_sps
    k = T_FMI%no_phi_f(j)
    T_FMI%f_phi_basis(1:k,j) = f_phi_basis_copy(1:k,j) + phi_tan
  end do

  k = FMI%mfi + 2
  do j = 1, FMI%no_spectro
    FMI%spectroscopic(j)%PHI_BASIS(1:k) = s_phi_basis_copy(1:k,j) + phi_tan
  end do
!
! Compute the ptg_angles (chi) for Antenna convolution, also the derivatives
! of chi w.r.t to T and other parameters
!
  Call get_chi_angles(ndx_path(1:,l),n_path(1:,l),FMI%tan_press,         &
 &     tan_hts(1:,l),tan_temp(1:,l),phi_tan,RoC,T_FMI%h_obs,elev_offset, &
 &     tan_dh_dt(1:,l,1:),no_tan_hts,T_FMI%no_t,T_FMI%t_zeta_basis,si,   &
 &     center_angle,ptg_angles(1:,l),dx_dt,d2x_dxdt,ier)
  IF(ier /= 0) goto 99

! Compute the refraction correction scaling matrix for this mmaf:
!
  Call refraction_correction(no_tan_hts, tan_hts(1:,l), h_path(1:,l), &
 &                n_path(1:,l), ndx_path(1:,l), E_rad, ref_corr)
!
! Now, Compute the radiances alone:
!
  prev_npf = -1
  Radiances(1:Nptg) = 0.0
!
! **********************  MAIN Pointing Loop *******************
!
  DO ptg_i = 1, no_tan_hts-1
!
    k = ptg_i
    h_tan = tan_hts(k,l)
    kk = FMI%no_ptg_frq(k)
!
    if(kk /= prev_npf) then
      prev_npf = kk
      DEALLOCATE(RadV,F_grid,STAT=i)
      ALLOCATE(RadV(kk),F_grid(kk),STAT=ier)
      IF(ier /= 0) then
        Print *,'** ALLOCATE Error: RadV or F_grid arrays, STAT =',ier
        goto 99
      endif
    endif
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
    no_ele = ndx_path(ptg_i,l)%total_number_of_elements
    Call get_beta_path(ptg_i,FMI%pfa_spectrum,no_ele,FMI%no_ptg_frq, &
   &     FMI%ptg_frq_grid,z_path(ptg_i,l),t_path(ptg_i,l),beta_path,ier)
    IF(ier /= 0) goto 99
!
    RadV(1:kk) = 0.0
    f_grid(1:kk) = FMI%ptg_frq_grid(k)%values(1:kk)
!
    do frq_i = 1, kk
!
      Frq = f_grid(frq_i)
!
      Call Rad_Tran(Frq, FMC%N_lvls, h_tan, FMI%n_sps, ndx_path(k,l),  &
     &    z_path(k,l), h_path(k,l), t_path(k,l), phi_path(k,l),&
     &    dHdz_path(k,l), T_FMI%earth_ref, beta_path(1:,frq_i),      &
     &    spsfunc_path(1:,k,l), ref_corr(1:,k), T_FMI%s_temp, brkpt, &
     &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
      IF(ier /= 0) goto 99
!
      RadV(frq_i) = Rad
!
    end do

! Frequency Average the radiances with the appropriate filter shapes
!
    i = 1
    ch = pfa_ch(i)
    if(FMC%do_frqavg) then
      Call Freq_Avg(F_grid,FMI%F_grid_filter(1:,i),  &
     &     FMI%Filter_func(1:,i),RadV,kk,FMI%no_filt_pts, &
     &     Radiances(ptg_i))
    else
      Radiances(ptg_i) = RadV(1)
    endif
!
  END DO              ! Pointing Loop
!
! Complete the radiances's last location
!
  kk = no_tan_hts
  Radiances(kk) = Radiances(kk-1)
!
  if(FMC%do_conv) then
!
! Here comes the Convolution code
!
    Call convolve_rad(T_FMI%ptg_press,FMI%tan_press,ptg_angles(1:,l), &
   &     band,center_angle,FMI%fft_pts,Radiances,no_tan_hts,i_star_all, &
   &     FMI%Xlamda,FMI%Aaap,FMI%D1Aaap,FMI%D2Aaap,FMI%Ias,ier)
    IF(ier /= 0) goto 99
!
  else
!
! Here comes the No_Convolution code
!
    Call no_convolve_rad (T_FMI%ptg_press,FMI%tan_press,Radiances, &
   &                      no_tan_hts,i_star_all)
!
  endif
!
! *** DEBUG Print
!
  if(FMC%do_conv) then
    Print *,'Convolution: ON'
  else
    Print *,'Convolution: OFF'
  endif
!
  Frq = FMC%Zfrq
  if(Frq > 0.0) then
    write(*,901) Frq
901 format(' Frequency Averaging: OFF',/,  &
         & ' (All computations done at Frq =',f12.4,')')
  else
    Print *,'Frequency Averaging: ON'
  endif
  Print *
!
  i = 1
  ch = pfa_ch(i)
  kk = T_FMI%ptg_press%no_lin_values
  write(*,903) ch,char(92),kk
  write(*,905) (i_star_all(k),k=1,kk)

903 format('ch',i2.2,'_avg_conv_pfa_rad',a1,i2.2)
905 format(4(2x,1pg15.8))
!
 99  CLOSE(11,iostat=i)
     CLOSE(13,iostat=i)
     CLOSE(32,iostat=i)
!
     if(io /= 0) Call ErrMsg(Line,io)

  Stop

!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
!
Subroutine Gti(Msg,Ti)
!
  Implicit none

  Character(len=*), intent(IN) :: Msg
  Character(len=*), intent(OUT) :: Ti
!
  Integer :: i

  i = Len(ti)
  Ti(1:i) = ' '
  i = Len_Trim(msg)
  write(*,900,Advance='NO') Msg(1:i)
  Read(5,'(a)',iostat=i) Ti
  i = Len_Trim(Ti)
  if(i > 1) Ti = AdjustL(Ti)

 900  Format(A,': ')
!
End Subroutine Gti

!---------------------------------------------------------------------------
! This subroutine transfers the derivatives over from the internal
! convolution grid to the users specified points. This module uses
! cubic spline interpolation to do the job.
!
Subroutine convolve_rad (ptg_press,tan_press,ptg_angles, &
           band,center_angle, fft_pts, i_raw, no_tan_hts, i_star_all,&
           Xlamda,Aaap,D1Aap,D2Aap,Ias,ier)

!
    integer(i4), intent(IN) :: IAS, no_tan_hts, band, fft_pts
!
    real(r8), intent(IN) :: CENTER_ANGLE, Xlamda
    real(r8), intent(IN) :: TAN_PRESS(*), PTG_ANGLES(*)
    real(r8), intent(IN) :: I_RAW(*)

    Real(r8), intent(in) :: AAAP(:,:),D1AAP(:,:),D2AAP(:,:)

!
    type(limb_press), intent(IN) :: PTG_PRESS
!
! -----     Output Variables   ----------------------------------------
!
    integer(i4), intent(OUT) :: IER
!
    real(r8), intent(OUT) :: I_STAR_ALL(:)
!
! -----     Local Variables     ----------------------------------------
!
    integer(i4) :: IS, J, K, NTR
!
    real(r8) :: Q, FFT_ANGLES(2**fft_pts), RAD(2**fft_pts)
!
! -----  Begin the code  -----------------------------------------
!
! Compute the ratio of the strengths
!
! This subroutine is called by channel
!
    Ier = 0
    ntr = 2**fft_pts
!
    Rad(1:no_tan_hts) = i_raw(1:no_tan_hts)
!
! Compute the convolution of the mixed radiances
!
    fft_angles(1:no_tan_hts) = ptg_angles(1:no_tan_hts)
    Call fov_convolve(fft_angles,Rad,center_angle,1,no_tan_hts,band, &
   &                  fft_pts,XLAMDA,AAAP,D1AAP,D2AAP,IAS,Ier)
    if (Ier /= 0) Return
!
! Interpolate the output values and store the radiances in: i_star_all
!
    is = 1
    k = no_tan_hts
    q = ptg_press%lin_val(1)-0.01
    do while(tan_press(is) < q)
      is = is + 1
    end do
    j = k - is + 1
    Call Cspline(fft_angles,ptg_angles(is:k),Rad,i_star_all,ntr,j)
!
    Return
!
  End Subroutine CONVOLVE_RAD

!---------------------------------------------------------------------------
! This subroutine transfers the derivatives over from the internal
! convolution grid to the users specified points. This module uses
! cubic spline interpolation to do the job.
!
Subroutine no_convolve_rad (ptg_press,tan_press,i_raw, no_tan_hts, i_star_all)
!
    Integer(i4), intent(IN) :: no_tan_hts
!
    Real(r8), intent(IN) :: TAN_PRESS(*), I_RAW(*)
!
    Type(limb_press), intent(IN) :: PTG_PRESS
!
! -----     Output Variables   ----------------------------------------
!
    Real(r8), intent(OUT) :: I_STAR_ALL(:)
!
! -----     Local Variables     ----------------------------------------
!
    Integer(i4) :: J, K
!
    Real(r8) :: PtP(Nlvl)
!
! -----  Begin the code  -----------------------------------------
!
    j = ptg_press%no_lin_values
    PtP(1:j) = dble(ptg_press%lin_val(1:j))
!
! Interpolate the output values and store the radiances in: i_star_all
!
    k = no_tan_hts
    Call Cspline(tan_press,PtP,i_raw,i_star_all,k,j)
!
    Return
!
  End Subroutine NO_CONVOLVE_RAD

end Program l2_numder
