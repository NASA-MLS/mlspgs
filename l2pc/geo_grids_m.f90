!
module GEO_GRIDS_M
  use GL6P
  use STRINGS, only: STRLWR
  use UNITS, only: MODEL_UNIT
  use L2PC_FILE_PARAMETERS, only: MAX_NO_BANDS, MAX_NO_KEY_ADDR, &
                                  MAX_NO_SV_ELMNTS, DEG2RAD, &
                                  MXCO => max_no_elmnts_per_sv_component
  use L2PC_PFA_STRUCTURES, only: GEOPHYS_PARAM
  use L2PCdim, only: MNP => max_no_phi, NCH, NLVL, NPTG, NSPS
  use MDBETA, only: MAX_NO_ZETA, CS_MNF => max_no_freq, NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use D_HUNT_M, only: HUNT
  use D_LINTRP_M, only: LINTRP
  use GET_ATMOS_M, only: READ_ZPM
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use HYDROSTATIC_INTRP, only: GET_HEIGHTS, GET_PRESSURES
  use CH_RUN_SETUP_M, only: CH_RUN_SETUP
  use VERT_TO_PATH_M, only: VERT_TO_PATH

  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------

SUBROUTINE geo_grids(time_stamp,geophys,no_geophys,ptg_press,    &
           z_gnlv,p_indx,ptg_hts,no_ptg,z_grid,t_grid,h_grid,    &
           dh_dt_grid,v_grid,n_lvls,mr_g,g_basis,no_coeffs_g,    &
           no_phi_g,t_index,si,conv_hts,conv_hts_raw,conv_press, &
           conv_temp,no_conv_hts,z_path,h_path,t_path,phi_path,  &
           dhdz_path,dh_dt_path,npath,g,lat,href,zref,phi_tan,   &
           roc,InDir,ld,n_tan,t_tan,path_brkpt,no_t,no_phi_t,    &
           t_phi_basis,ier)

! This module sets up the geophysical coefficients

!  ===============================================================
!  Declaration of variables for sub-program: geo_grids
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Character (LEN=*), INTENT(IN) :: time_stamp, InDir

Integer(i4), INTENT(IN) :: p_indx(:), ld, no_geophys, t_index, &
                           npath, no_ptg, si, n_lvls

Integer(i4), INTENT(IN OUT) :: no_conv_hts

Integer(i4), INTENT(OUT) :: no_t, no_phi_t, no_coeffs_g(:), no_phi_g(:), &
             path_brkpt(:,:), n_tan(:), Ier

Real(r4), INTENT(OUT) :: dh_dt_path(:,:,:,:), dhdz_path(:,:)

Real(r8), INTENT(IN) :: roc, phi_tan, lat, href, zref, z_gnlv(:),  &
          ptg_press(:)

Real(r8), INTENT(OUT) :: g, z_grid(:), t_grid(:), h_grid(:), mr_g(:,:,:),   &
   g_basis(:,:), conv_temp(:), conv_press(:), conv_hts(:), dh_dt_grid(:,:), &
   t_phi_basis(:), t_path(:,:), phi_path(:,:), z_path(:,:), h_path(:,:),    &
   t_tan(:), ptg_hts(:)

Real(r8), INTENT(IN OUT) :: conv_hts_raw(:), v_grid(:)

type (geophys_param) :: geophys(*)
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: ngt = (ng+1) * nlvl
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: g_i, no_phi, no_pts, io, h_i, jp, i, j, k, klo, khi, gl_count

Real(r8) :: reff, const, q, z1, z2, xm, ym, h, z, t, v, dhdt(mxco,mnp),     &
    mr_g_geo(nlvl,mnp),phi,tri_base(mxco),z_glgrid(ngt),h_glgrid(ngt,nlvl), &
    t_glgrid(ngt,nlvl), dh_dt_glgrid(ngt,mnp,nlvl), dhdz_glgrid(ngt,nlvl)

Character (LEN=80) :: fn
Character (LEN=28) :: buff

! *** DEBUG - These are the V51 Temperature values...

Real(r8) :: v51_t(23) = (/                                        &
     &      300.661, 282.276, 263.885, 243.500, 223.101, 211.353, &
     &      199.592, 205.740, 211.886, 216.286, 220.679, 226.842, &
     &      233.007, 240.569, 248.137, 255.747, 263.363, 265.913, &
     &      268.462, 263.453, 258.436, 249.422, 240.412 /)

! *** END DEBUG

! Begin the code here

  ier = 0

  DO g_i = 1, no_geophys

! Read other geophysical parameters from files. Assemble the File name:

    buff(1:)=' '
    IF(g_i == t_index) THEN
      buff = 'zpmodel_d_10n_tempera.data'
      buff(9:13) = time_stamp(1:5)
    ELSE
      WRITE(buff,'(''ZPMODEL_'',A7,''_'',A5,''.DATA'')')  &
          time_stamp(1:7),geophys(g_i)%name
    END IF

    fn(1:)=' '
    fn = InDir(1:ld)//buff(1:)
    CALL strlwr(fn)
    CALL read_zpm(fn,model_unit,no_pts,no_phi,v_grid,mr_g_geo,io)
    IF(io /= 0) THEN
      ier = 1
      PRINT *,' Error in: geo_grids subroutine ..'
      CALL errmsg(fn,io)
      RETURN
    END IF

    no_phi_g(g_i) = no_phi
    jp = (no_phi + 1) / 2
    k = geophys(g_i)%no_lin_values
    no_coeffs_g(g_i) = k

    DO i = 1, k
      tri_base(i) = geophys(g_i)%basis_peaks(i)
    END DO

! Loop ove r Phi's

    DO j = 1, no_phi

! Get the coefficient values

      CALL lintrp(v_grid,tri_base,mr_g_geo(1:,j),mr_g(1:,j,g_i),no_pts,k)

    END DO

    IF(g_i == t_index) THEN

! *** DEBUG ***  "Cooked up" temperature profile, copy of V51 values

      DO i = 1, k
        mr_g(i,jp,g_i) = v51_t(i)
      END DO

! *** END DEBUG ***

! *** DEBUG ***
! Make all Temperature Phi values equal to center Phi, thus No gradient
! in Temperature in the Phi direction

      DO j = 1, no_phi
        IF(j /= jp) THEN
          DO i = 1, k
            mr_g(i,j,g_i) = mr_g(i,jp,g_i)
          END DO
        END IF
      END DO

! *** END DEBUG

! Define the Phi base coefficients:

      no_t = k
      phi = 1.5 * deg2rad
      no_phi_t = no_phi
      DO j = 1, no_phi
        t_phi_basis(j) = phi_tan + (j - jp) * phi
      END DO

! Copy user's coefficients pressure grid for Temperature onto g_basis:

      DO i = 1, k
        g_basis(i,g_i) = tri_base(i)
      END DO

! Get the selected integration grid pressures. Also, define the GL pressure
! grid:

      DO i = 1, n_lvls
        j = p_indx(i)
        v_grid(i) = 0.0
        z_grid(i) = z_gnlv(j)
      END DO

      gl_count = 0
      z2 = DBLE(z_grid(1))
      DO i = 2, n_lvls
        z1 = z2
        z2 = DBLE(z_grid(i))
        xm = 0.5D0 * (z2 + z1)
        ym = 0.5D0 * (z2 - z1)
        gl_count = gl_count + 1
        z_glgrid(gl_count) = z1
        DO j = 1, ng
          gl_count = gl_count + 1
          z_glgrid(gl_count) = xm + ym * gx(j)
        END DO
      END DO

      gl_count = gl_count + 1
      z_glgrid(gl_count) = z2

! *** Make sure there is NO Extrapolation on Temperature at either end of
!     the range. This is the Temperature on the GL grid:

      DO j = 1, no_phi
        klo = 1
        z = mr_g(1,j,g_i)
        DO WHILE(tri_base(1) > z_glgrid(klo))
          t_glgrid(klo,j) = z
          klo = klo + 1
        END DO
        khi = gl_count
        z = mr_g(k,j,g_i)
        DO WHILE(tri_base(k) < z_glgrid(khi))
          t_glgrid(khi,j) = z
          khi = khi - 1
        END DO
        i = khi - klo + 1
        CALL lintrp(tri_base,z_glgrid(klo:),mr_g(1:,j,g_i), &
                    t_glgrid(klo:,j),k,i)
      END DO

! Get Hydrostatically balanced group ([zth]_glgrid, dhdz_glgrid and
! dh_dt_glgrid on the GL grid)

      DO h_i = 1, gl_count
        z = z_glgrid(h_i)
        DO j = 1, no_phi
          phi = t_phi_basis(j)
          CALL get_h_dhdt(mr_g(1:,1:,g_i),tri_base,t_phi_basis,z,phi, &
              zref,href,lat,reff,g,const,nlvl,mxco,k,no_phi,h,dhdt)
          h_glgrid(h_i,j) = h
          q = h + reff
          dhdz_glgrid(h_i,j) = q * q * const * t_glgrid(h_i,j)
          DO i = 1, k
            dh_dt_glgrid(h_i,j,i) = dhdt(i,j)
          END DO
        END DO
      END DO

! Define the [zth]_grid Hydrostatically balanced group as a subset
! of the [zth]_glgrid group, at the "Center Phi"

      i = -ng
      h_i = 0
      j = ng + 1
      DO WHILE(i < gl_count)
        i = i + j
        h_i = h_i + 1
        h_grid(h_i) = h_glgrid(i,jp)
        t_grid(h_i) = t_glgrid(i,jp)
        DO khi = 1, k
          dh_dt_grid(h_i,khi) = dh_dt_glgrid(i,jp,khi)
        END DO
      END DO

      DO i = h_i+1, n_lvls+1
        z_grid(i) = z_glgrid(gl_count)
        h_grid(i) = h_glgrid(gl_count,jp)
        t_grid(i) = t_glgrid(gl_count,jp)
        DO khi = 1, k
          dh_dt_grid(i,khi) = dh_dt_glgrid(gl_count,jp,khi)
        END DO
      END DO

    ELSE IF(geophys(g_i)%name == 'VEL_Z') THEN  ! Create v_grid

      CALL lintrp(ptg_press,z_grid,mr_g_geo(1:,jp),v_grid,no_pts,n_lvls)

    END IF

  END DO

!  Set conv_hts_raw to be a REAL SUBSET of h_grid, so NO height interpolation
!  for abs_cs will occure, also define conv_press as a TRUE subset of z_grid:

  klo = -1
  DO i = si+1, no_conv_hts
    h = conv_hts_raw(i)
    CALL hunt(h,h_grid,n_lvls,klo,khi)
    IF(ABS(h-h_grid(khi)) < ABS(h-h_grid(klo))) klo = khi
    conv_press(i) = z_grid(klo)
    conv_hts_raw(i) = h_grid(klo)
  END DO

!  Make sure conv_hts array is stricktly monotonically increasing:

  j = no_conv_hts
  no_conv_hts = 1
  v = conv_hts_raw(1)
  DO i = 2, j
    z = conv_press(i)
    h = conv_hts_raw(i)
    IF(h > v) THEN
      v = h
      no_conv_hts = no_conv_hts + 1
      conv_press(no_conv_hts) = z
      conv_hts_raw(no_conv_hts) = h
    END IF
  END DO

! Interpolate the hydrostatic grid for output grid heights (ptg_hts)

  CALL get_heights('h',h_grid,t_grid,z_grid,n_lvls,ptg_press, &
                    ptg_hts,no_ptg,ier)
  IF(ier /= 0) RETURN

! Interpolate the hydrostatic grid for conv. grid pressures (conv_press)
! for the values BELOW Earth surface only:

  CALL get_pressures('h',h_grid,t_grid,z_grid,n_lvls,conv_hts_raw, &
      conv_press,si,ier)
  IF(ier /= 0) RETURN

! Interpolate the hydrostatic grid for convolution grid temperatures

  ier = 0
  h_i = 1
  z = z_grid(h_i)
  t = t_grid(h_i)
  DO WHILE(conv_press(h_i) < z)
    conv_temp(h_i) = t
    h_i = h_i + 1
  END DO

  CALL lintrp(z_grid,conv_press(h_i:),t_grid,conv_temp(h_i:), &
              n_lvls,no_conv_hts-h_i+1)

  IF(n_lvls < nlvl) h_grid(n_lvls+1) = h_grid(n_lvls)

! Define the final convolution heights array.

  DO i = 1, no_conv_hts
    conv_hts(i) = conv_hts_raw(i)
  END DO

! Get n_tan & t_tan arrays:

  CALL ch_run_setup(conv_hts,no_conv_hts,h_grid,t_grid,n_lvls, n_tan,t_tan)

! Compute all the various integration paths according to tanget heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  DO k = 1, no_conv_hts
    h = conv_hts(k)
    CALL vert_to_path(n_lvls,ng,npath,gl_count,no_phi_t, &
        no_t,n_tan(k),h,z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid, &
        dh_dt_glgrid,t_phi_basis,z_path(1:,k),h_path(1:,k), &
        t_path(1:,k),phi_path(1:,k),dhdz_path(1:,k), &
        dh_dt_path(1:,1:,1:,k),path_brkpt(1:,k),roc,phi_tan)
  END DO

  RETURN
END SUBROUTINE geo_grids

!----------------------------------------------------------------------
!  The 2 dimensional Hydrostatic integrator

SUBROUTINE get_h_dhdt(t_profile,t_basis,t_phi_basis,zeta,phi,zref, &
           href,lat,reff,g,const,nlvl,nc,nt,no_t_phi,h,dhdt)

!  ===============================================================
!  Declaration of variables for sub-program: get_h_dhdt
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: nlvl, nc     ! Not Used

Integer(i4), INTENT(IN) :: nt, no_t_phi

Real(r8), INTENT(IN) :: lat, href, zref, zeta, phi, t_profile(:,:), &
             t_basis(:), t_phi_basis(:)

Real(r8), INTENT(OUT) :: const, h, g, reff, dhdt(:,:)
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: boltzxln10 = 19.14486942d0
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: iq, j, nv

Real(r8) :: v, denom, sum, eta_phi

Real(r8), save :: gr2, c, prevlat = -500.0

Real(r8) :: pm(50) = 0.0
Real(r8) :: etapv(15) = 0.0

! Begin the code here

  IF(lat /= prevlat) THEN
    prevlat = lat
    CALL get_g_reff(lat,href,g,reff)
    c = g * reff
    gr2 = c * reff
    const = boltzxln10 / (gr2 * 28.964125D0)      ! 28.964... = Mass
  END IF

  nv = 0
  h = href
  DO j = 1, nt
    Call pq_ana(zref,zeta,t_basis,j,nt,v)
    pm(j) = v
    IF(v /= 0.0) nv = nv + 1
    DO iq = 1, no_t_phi
      dhdt(j,iq) = 0.0
    END DO
  END DO

  IF(nv < 1) RETURN

  nv = 0
  sum = 0.0D0
  DO iq = 1, no_t_phi
    CALL get_one_eta(phi,t_phi_basis,no_t_phi,iq,eta_phi)
    etapv(iq) = eta_phi
    IF(eta_phi > 0.0) THEN
      nv = nv + 1
      DO j = 1, nt
        sum = sum + t_profile(j,iq) * eta_phi * pm(j)
      END DO
    END IF
  END DO

  IF(nv < 1) RETURN

  IF(sum /= 0.0D0) THEN
    denom = c - boltzxln10 * sum
    h = gr2 / denom - reff + href
    v = gr2 * boltzxln10 / (denom * denom)
  ELSE
    h = href
    v = boltzxln10 / g
  END IF

  DO iq = 1, no_t_phi
    eta_phi = etapv(iq)
    IF(eta_phi > 0.0) THEN
      DO j = 1, nt
        dhdt(j,iq) = v * eta_phi * pm(j)
      END DO
    END IF
  END DO

  RETURN
END SUBROUTINE get_h_dhdt

!---------------------------------------------------------------------

SUBROUTINE get_g_reff(lat,href,g,reff)

!  ===============================================================
!  Declaration of variables for sub-program: get_g_reff
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Real(r8), INTENT(IN) :: lat, href

Real(r8), INTENT(OUT) :: g, reff
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: g0 = 9.80616d0
Real(r8), PARAMETER :: g1 = 2.6373d-3
Real(r8), PARAMETER :: g2 = 5.9d-6
Real(r8), PARAMETER :: r1 = 3.085462d-3
Real(r8), PARAMETER :: r2 = 2.27d-6
Real(r8), PARAMETER :: r3 = 2.0d-9
Real(r8), PARAMETER :: d2r = 1.74532925199433d-2
!  ----------------
!  Local variables:
!  ----------------
Real(r8) :: ang, cr, crr, a

! Begin code:

  ang  = d2r * lat                                ! Convert to radians
  cr   = DCOS(2.0D0 * ang)
  g    = g0 * (1.0D0 - cr * (g1 - cr * g2))       ! Modified G
  crr  = 2.0D0 * cr * cr - 1.0D0                  ! Cos(4*Ang)
  reff = 2.0D0 * g / (r1 + r2 * cr - r3 * crr)    ! Reff in Kilometers

! Make approriate correction for Href not being the surface

  a = (1.0D0 + href / reff)
  g = g / (a * a)
  reff = reff + href

  RETURN
END SUBROUTINE get_g_reff

!---------------------------------------------------------------------

SUBROUTINE pq_ana(lower_lim,upper_lim,zetabase,iq,nt,v)

!  ===============================================================
!  Declaration of variables for sub-program: pq_ana
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: iq, nt

Real(r8), INTENT(IN) :: lower_lim, upper_lim, zetabase(:)
Real(r8), INTENT(OUT) :: v
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), PARAMETER :: m0 = 28.964125d0
!  ----------------
!  Local variables:
!  ----------------
Real(r8) :: sgn, z_0, z_i, h, l

! Begin code:

  v = 0.0_r8
  IF(ABS(upper_lim-lower_lim) < 1.0E-4) RETURN

  sgn = 1.0
  z_0 = lower_lim
  z_i = upper_lim
  IF(lower_lim > upper_lim) THEN
    sgn = -1.0
    z_0 = upper_lim
    z_i = lower_lim
  END IF

  IF(iq > 1) THEN

! Standard lower triangular integration

    h = MIN(zetabase(iq),MAX(z_i,zetabase(iq-1)))
    l = MAX(zetabase(iq-1),MIN(z_0,zetabase(iq)))
    v = (h - l) * (0.5 * (h + l) - zetabase(iq-1)) /  &
     &  (m0 * (zetabase(iq) - zetabase(iq-1)))

  ELSE

! Special lower rectangular integration

    h = MIN(z_i,zetabase(1))
    l = MIN(z_0,zetabase(1))
    v = (h - l) / m0

  END IF

  IF(iq < nt) THEN

! Standard upper triangular integration

    h = MIN(zetabase(iq+1),MAX(z_i,zetabase(iq)))
    l = MAX(zetabase(iq),MIN(z_0,zetabase(iq+1)))
    v = sgn * (v + (h - l) * (zetabase(iq+1) - 0.5 * (h + l)) /  &
              (m0 * (zetabase(iq+1) - zetabase(iq))) )

  ELSE

! Special upper rectangular integration

    h = MAX(z_i,zetabase(nt))
    l = MAX(z_0,zetabase(nt))
    v = sgn * (v + (h - l) / m0)

  END IF

  RETURN
END SUBROUTINE pq_ana

end module GEO_GRIDS_M
! $Log$
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
