!
module FWD_MDL_SET_UP_M
  use EOS_MDB, only: EOS_MDB_HDR, EOS_MDB_REC, MAX_NO_LINES
  use L2PC_FILE_PARAMETERS, only: MAX_NO_BANDS, MAX_NO_KEY_ADDR, &
                                  DEG2RAD, MAX_NO_SV_ELMNTS
  use L2PC_FILE_STRUCTURES, only: L2PC_HEADER_ONE
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, GEOPHYS_PARAM, &
                           LIMB_PRESS, MAXAITKENPTS, MAXFILTPTS, MAXGEOM, &
                           PFA_SLAB, SPECTRO_PARAM, EARTH_MAJOR, EARTH_MINOR
  use L2PCdim, only: MNP => max_no_phi, NCH, NLVL, NPTG, NSPS
  use MDBETA, only: MAX_NO_ZETA, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
       HT2, RR, PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  use GEO_GRIDS_M, only: GEO_GRIDS
  use GET_CS_M, only: GET_CS
  use ATMOS_BASIS_M, only: ATMOS_BASIS
  use REFRACTION_M, only: REFRACTIVE_INDEX
  use PFA_PREP_M, only: PFA_PREP
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!---------------------------------------------------------------------

SUBROUTINE fwd_mdl_set_up(time_stamp,primag,href,zref,ptg_hts,     &
     &     z_gnlv,p_indx,e_rad,z_grid,t_grid,h_grid,dh_dt_grid,    &
     &     n_grid,v_grid,vel_z,mdb_pres,mdb_temp,mdb_freq,         &
     &     no_freqs_f,n_lvls,ptg_press,geometric,no_geom,          &
     &     geophys,no_geophys,g_basis,mr_g,no_coeffs_g,            &
     &     atmospheric,no_atmos,f_basis,mr_f,no_coeffs_f,freq,cs,  &
     &     spsfunc,si,t_index,earth_ref,s_temp,h_obs,n_obs,        &
     &     conv_hts_raw,conv_hts,conv_press,conv_temp,no_conv_hts, &
     &     no_sps_tbl,sps_tbl,no_pfa_ch,no_filt_pts,pfa_ch,        &
     &     no_int_frqs,pfa_spectrum,t_tan,ndx_sps,n_tan_ptr,       &
     &     f_grid_fltr,f_grid,fltr_func,ch1,ch2,InDir,ld,fnd,      &
     &     lf,hdr1,no_phi_g,no_phi_f,phi_basis_f,atmos_index,      &
     &     geophys_index,no_phi_vec,z_path,h_path,t_path,          &
     &     phi_path,dhdz_path,dh_dt_path,npath,path_brkpt,         &
     &     no_phi_t,t_phi_basis,rad_cur,mdb_hdr,mdb_rec,           &
     &     spectroscopic,no_spectro,spect_index,ier)

!  ===============================================================
!  Declaration of variables for sub-program: fwd_mdl_set_up
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_geophys, no_atmos, ch1, ch2, si,    &
             p_indx(:), atmos_index(:), geophys_index(:), spect_index(:), &
             no_spectro, npath, no_pfa_ch, no_filt_pts, pfa_ch(:), lf,    &
             no_int_frqs(:)
!
Integer(i4), INTENT(IN OUT) :: no_geom, no_conv_hts, ld

Integer(i4), INTENT(OUT) :: n_obs, no_sps_tbl(:), sps_tbl(:,:), t_index,   &
             no_phi_vec(:), no_coeffs_g(:), no_coeffs_f(:), no_freqs_f(:), &
             ier, no_phi_f(:), no_phi_g(:), no_phi_t, path_brkpt(:,:),     &
             ndx_sps(:,:), n_tan_ptr(:)
!
Real(r4), INTENT(OUT) :: dh_dt_path(:,:,:,:)
!
Real(r8), INTENT(IN) :: vel_z, href, zref, z_gnlv(:)
!
Real(r8), INTENT(OUT) :: conv_hts(:)
Real(r8), INTENT(IN OUT) :: n_grid(:), cs(:,:,:,:), conv_hts_raw(:), &
                            v_grid(:)
!
Real(r4), INTENT(OUT) :: dhdz_path(:,:)
Real(r8), INTENT(OUT) :: e_rad, earth_ref, s_temp, h_obs, freq(:), &
          mr_f(:,:,:), phi_basis_f(:,:), mr_g(:,:,:),g_basis(:,:), &
          f_basis(:,:), t_phi_basis(:),t_path(:,:), z_path(:,:),   &
          h_path(:,:),conv_temp(:), phi_path(:,:), mdb_pres(:),    &
          z_grid(:), t_grid(:),h_grid(:), conv_press(:), t_tan(:), &
          rad_cur, f_grid(:,:), fltr_func(:,:), f_grid_fltr(:,:),  &
          dh_dt_grid(:,:), ptg_hts(:), spsfunc(:,:), mdb_temp(:),  &
          mdb_freq(:,:)
!
Character (LEN=*), INTENT(IN) :: fnd
Character (LEN=*), INTENT(IN) :: InDir
Character (LEN=*), INTENT(IN) :: primag
Character (LEN=*), INTENT(IN OUT) :: time_stamp
!  ----------------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, jp, ncpb, ich, sps_i, geo_i, band, ht_i, geo_j, &
               ih2o, kk, no_geom_new, ncs, no_t, band1, band2

Real(r4) :: def_val(maxgeom) = (/                                    &
            0.0, 0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.0, 6972.0, 6372.0, &
            0.05, 2.735, 1.0 /)

Real(r8) :: q, r, p, geoc_lat, sw, cw, h2o_grid(nlvl), h2o_corr(nlvl), &
            g, rp, b2, incl, beta, PtP(nlvl)

Character (LEN=8) :: geomname(maxgeom) = (/                         &
         &  'ELEV_183','ELEV_205','AZIM_183','AZIM_205','AZIM_REF', &
         &  'ROLL    ','PITCH   ','YAW     ','GEOCSRAD','GEOCERAD', &
         &  'EARTHREF','SPACE_T ','LOOK_DIR'/)

!  PFA variables:

type (l2pc_header_one) :: hdr1

type (limb_press) :: ptg_press
type (geom_param) :: geometric(*)
type (geophys_param) :: geophys(*)
type (atmos_comp) :: atmospheric(*)
type (spectro_param) :: spectroscopic(*)
!
type (pfa_slab) :: pfa_spectrum(6,*)
!
type (eos_mdb_hdr) :: mdb_hdr(*)
type (eos_mdb_rec) :: mdb_rec(max_no_lines,*)

! Set up the default geometric linearization values of earth radius,earth
! emission,background space temperature and satellite position.
! Compute the default geocentric earth radius

  ier = 0
  READ(time_stamp(3:4),'(f2.0)') geoc_lat
  IF(time_stamp(5:5) == 'S') geoc_lat = -geoc_lat

  q = DBLE(earth_major)
  a2 = q * q

  r = DBLE(earth_minor)
  b2 = r * r

  beta = 98.0D0                    ! Orbit inclination angle for EOS-MLS

  incl = (beta - 90.0D0) * deg2rad
  q = TAN(incl)
  r = q * q
  c2 = (1.0D0+r)*a2*b2/(a2+b2*r)   ! This is c*c
  q = SQRT(c2)                     ! This is c, Minor axis for 2D ellipse

  c2oa2 = c2 / a2

! Get the Tangent Phi (Geodetic Lat.) from the Geocentric Lat.

  r = geoc_lat * deg2rad
  q = SIN(r)
  sw = q * q
  r = COS(incl)
  cw = r * r
  q = a2*a2*sw*cw/(b2*b2+(a2*a2-b2*b2)*sw*cw)
  spt = SQRT(q)                ! Sin(Phi_tan)
  phi_tan = DASIN(spt)

  cpt = COS(phi_tan)           ! Cos(Phi_tan)
  cw = cpt * cpt
  sw = q
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

! Note that "geodetic" convolution heights must be saved for future calls
! Set up the radiative function evaluted at users heights
! Find temp index (t_index) into the geophys structure

  no_t = 0
  ih2o = 0
  geo_i = 0
  t_index = 0
  DO WHILE(geo_i < no_geophys)
    geo_i = geo_i + 1
    IF(geophys(geo_i)%name(1:4) == 'TEMP') THEN
      t_index = geo_i
      no_t = geophys(geo_i)%no_lin_values
      no_coeffs_g(geo_i) = no_t
      geo_i = no_geophys + 2
    END IF
  END DO

! Call the grids program to get the preselected integration heights and
! pointings. Also transfer the pointings expressed in pressure units to
! heights.  Also, set up geophysical parameters (for derivatives):
! [The following routine replaces the older code using the two routines:
!  grids() and geo_basis() ]

! (z_grid,t_grid,h_grid) are the arrays of the preselected integration grid.
! (ptg_hts,ptg_press) are the arrays used to establish the radiative transfer
! function after antenna convolution at the users instument tangent heights.

  kk = ptg_press%no_lin_values
  PtP(1:kk) = DBLE(ptg_press%lin_val(1:kk))
  CALL geo_grids(time_stamp,geophys,no_geophys,PtP, &
      z_gnlv,p_indx,ptg_hts,kk,z_grid,t_grid,h_grid, &
      dh_dt_grid,v_grid,n_lvls,mr_g,g_basis,no_coeffs_g, &
      no_phi_g,t_index,si,conv_hts,conv_hts_raw,conv_press, &
      conv_temp,no_conv_hts,z_path,h_path,t_path,phi_path, &
      dhdz_path,dh_dt_path,npath,g,geoc_lat,href,zref,phi_tan, &
      roc,InDir,ld,n_tan_ptr,t_tan,path_brkpt,no_t,no_phi_t,t_phi_basis, &
      ier)
  IF(ier /= 0) RETURN

! Connect between 'no_phi_vec' and 'no_phi_g' entries.

  DO i = 1, no_geophys
    j = geophys_index(i)
    kk = hdr1%sv_component_first_elmnt_index(j)
    no_phi_vec(kk) = no_phi_g(i)
  END DO

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

  n_obs = -n_lvls
  DO band = band1, band2
    j = no_sps_tbl(band)
    DO kk = 1, j
      i = sps_tbl(kk,band)
      no_freqs_f(i) = 0
      IF(atmospheric(i)%fwd_calc(band)) THEN
        CALL get_cs(primag, atmospheric(i)%spectag, n_obs, p_indx, ncs, &
             mdb_pres, mdb_temp, mdb_freq(1:,i), cs(1:,1:,1:,i), fnd, &
             ch1, ier)
        IF(ier /= 0) RETURN
        no_freqs_f(i) = ncs
      END IF
    END DO
  END DO

  CLOSE(38,IOSTAT=i)

! Set n_obs to be: N_lvls always.
!   (Changed, Aug/6/96 Z.Shippony & W.G.Read)

  n_obs = n_lvls

! Now get the mixing ratios

  CALL atmos_basis(time_stamp,atmospheric,no_atmos,z_grid,  &
                   spsfunc,n_lvls,mr_f,f_basis,no_coeffs_f, &
                   phi_tan,no_phi_f,phi_basis_f,InDir,ld,ier)
  IF(ier /= 0) RETURN

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

  DO geo_i = 1, no_geophys
    j = no_phi_g(geo_i) * no_coeffs_g(geo_i)
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

    ih2o = sps_i
    kk = no_coeffs_f(ih2o)
    jp = (no_phi_f(ih2o) + 1) / 2            ! Temporary code ..
    CALL refractive_index(mr_f(1:,jp,ih2o),f_basis(1:,ih2o),kk,    &
     &                    z_grid,t_grid,n_grid,h2o_grid,n_lvls)

  ELSE

! Ignore water contribution to refractive index, i.e:  dry air

    DO ht_i = 1, n_lvls
      h2o_grid(ht_i) = 0.0
      p = 10.0**(-z_grid(ht_i))
      n_grid(ht_i) = 7.76E-5 * p / t_grid(ht_i)
    END DO

  END IF

! Create filter grids & functions for PFA calculations

  CALL pfa_prep(atmospheric,no_atmos,no_pfa_ch,n_lvls,no_filt_pts, &
       pfa_ch,no_int_frqs,pfa_spectrum,z_grid,t_grid,h2o_grid,     &
       ndx_sps,f_grid_fltr,f_grid,fltr_func,freq,vel_z,            &
       h2o_corr,InDir,ld,fnd,lf,primag, mdb_hdr,mdb_rec,ier)
  IF(ier /= 0) RETURN

  IF(ih2o < 1) RETURN

!  Correct H2O (18003) cross_sections for the non-linear part:

  jp = no_phi_f(ih2o)
  DO ht_i = 1, n_lvls
    p = h2o_corr(ht_i)
    IF(ABS(1.0-p) >= 1.0E-4) THEN
      j = 0
      DO ich = ch1, ch2
        j = j + 1
        DO kk = 1, jp
          cs(ht_i,kk,j,ih2o) = p * cs(ht_i,kk,j,ih2o)
        END DO
      END DO
    END IF
  END DO

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

!  dfloat not an intrinsic function
!  ch_offset = dfloat(j)
  ch_offset = real(j, kind(ch_offset))
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
! Revision 1.2  2000/08/16 01:02:22  zvi
! Correcting array bounds erros in calling seq. SUN & SGI f90 pass on these
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
