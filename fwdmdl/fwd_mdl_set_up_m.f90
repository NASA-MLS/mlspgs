module FWD_MDL_SET_UP_M
  use L2PC_FILE_PARAMETERS, only: MAX_NO_BANDS, MAX_NO_KEY_ADDR, &
                                  DEG2RAD, MAX_NO_SV_ELMNTS
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, PFA_SLAB,        &
                                 SPECTRO_PARAM, MAXAITKENPTS, MAXFILTPTS, &
                                 MAXGEOM, EARTH_MAJOR, EARTH_MINOR
  use L2PCDIM, only: MNP => max_no_phi, NCH, Nlvl, Nptg, NSPS
  use MDBETA, only: MAX_NO_ZETA, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
      HT2, RR, PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  use GEO_GRIDS_M, only: GEO_GRIDS
  use PATH_ENTITIES_M, only: PATH_BETA, PATH_INDEX, PATH_VECTOR, &
                             PATH_DERIVATIVE
  use GL6P, only: NG, GX
  use REFRACTION_M, only: REFRACTIVE_INDEX
  use PFA_PREP_M, only: PFA_PREP
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!---------------------------------------------------------------------

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
           is_f_log,beta_path,ier)

!  ===============================================================
!  Declaration of variables for sub-program: fwd_mdl_set_up
!  ===============================================================

Integer(i4), PARAMETER :: ngt = (Ng+1) * Nlvl

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_atmos, no_t, ch1, ch2, si,  &
             atmos_index(:), spect_index(:), no_phi_t, no_freqs, &
             no_coeffs_f(:), no_phi_f(:), &
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
          if (is_f_log(j)) then
            Call TWO_D_POLATE(f_basis(1:,j), LOG(mr_f(1:kk,1:jp,j)), &
     &                        kk, phi_basis_f(1:,j), jp, zeta, phi, r)
            q = exp(r)
          else
            Call TWO_D_POLATE(f_basis(1:,j), mr_f(1:kk,1:jp,j), kk,  &
     &                        phi_basis_f(1:,j), jp, zeta, phi, q)
          endif
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
         beta_path,freq_grid,InDir,ld,primag,ier)
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

end module FWD_MDL_SET_UP_M
! $Log$
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
