module GEO_GRIDS_M
  use GL6P, only: NG, GX
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use L2PCDIM, only: N2lvl, Nlvl, Nptg, MNP => max_no_phi
  use MDBETA, only: NO_T_PHI
  use MLSCommon, only: I4, R4, R8
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_DERIVATIVE
  use D_HUNT_M, only: HUNT
  use D_LINTRP_M, only: LINTRP
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use HYDROSTATIC_INTRP, only: GET_HEIGHTS, GET_PRESSURES
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

SUBROUTINE geo_grids(z_grid,t_grid,h_grid,tan_dh_dt,n_lvls,t_z_basis,    &
           t_coeff,si,tan_hts,tan_hts_raw,tan_press,tan_temp,no_tan_hts, &
           z_path,h_path,t_path,phi_path,dhdz_path,dh_dt_path,npath,g,   &
           href,zref,geoc_lat,no_t,no_phi_t,t_phi_basis,ndx_path,Ier)

!  ===============================================================
!  Declaration of variables for sub-program: geo_grids
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------

Integer(i4), INTENT(IN) :: npath, si, n_lvls, no_phi_t, no_t

Integer(i4), INTENT(IN OUT) :: no_tan_hts

Integer(i4), INTENT(OUT) :: Ier

Type(path_index), INTENT(OUT)  :: ndx_path(:)
Type(path_vector), INTENT(OUT) :: z_path(:),t_path(:),h_path(:), &
                                  phi_path(:),dhdz_path(:)

Type(path_derivative), INTENT(OUT) :: dh_dt_path(:)

Real(r8), INTENT(IN) :: geoc_lat, href(:), zref(:), &
            z_grid(:), t_z_basis(:), t_coeff(:,:), t_phi_basis(:)

Real(r8), INTENT(OUT) :: g, t_grid(:), h_grid(:), tan_temp(:), tan_hts(:), &
          tan_press(:),  tan_dh_dt(:,:)

Real(r8), INTENT(IN OUT) :: tan_hts_raw(:)

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: no_phi, h_i, jp, i, j, k, klo, khi, gl_count

Real(r4) :: dhdtp(ngt,mnp,mxco)

Real(r8) :: reff, const, q, h, z, t, v, phi, dummy(Nlvl),z1,z2,xm,ym, &
     dhdt(mxco,mnp), h_glgrid(ngt,mnp), t_glgrid(ngt,mnp), &
     dh_dt_glgrid(ngt,mnp,mxco), dhdz_glgrid(ngt,mnp), zpath(ngt), &
     tpath(ngt),hpath(ngt),ppath(ngt),dhdzp(ngt),z_glgrid(ngt/2)

! Begin the code here

  ier = 0

  no_phi = no_phi_t
  jp = (no_phi + 1) / 2

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
! *** Create the Temperature on the GL grid by linear interpolation via
!     the get_eta procedure (Bill's request)
!
  DO j = 1, no_phi
    do i = 1, gl_count
      t = 0.0
      z = z_glgrid(i)
      do h_i = 1, no_t
        CALL get_one_eta(z,t_z_basis,no_t,h_i,v)
        t = t + t_coeff(h_i,j) * v
      end do
      t_glgrid(i,j) = t
    end do
  END DO

! Get Hydrostatically balanced group ([zth]_glgrid, dhdz_glgrid and
! dh_dt_glgrid on the GL grid)

  DO h_i = 1, gl_count
    z = z_glgrid(h_i)
    DO j = 1, no_phi_t
      phi = t_phi_basis(j)
      CALL get_h_dhdt(t_coeff(1:,1:),t_z_basis,t_phi_basis,z,phi, &
           zref(j),href(j),geoc_lat,reff,g,const,no_t,no_phi,h,dhdt)
      h_glgrid(h_i,j) = h
      q = h + reff
      dhdz_glgrid(h_i,j) = q * q * const * t_glgrid(h_i,j)
      DO i = 1, no_t
        dh_dt_glgrid(h_i,j,i) = dhdt(i,j)
      END DO
    END DO
  END DO

! Define the [th]_grid Hydrostatically balanced group as a subset
! of the [th]_glgrid group, at the "Center Phi"  (the (independent) z
! (zeta) component has been define outside of this routine ..)

  j = Ng + 1
  klo = min(Nlvl,n_lvls+1)
  h_i = (gl_count + Ng) / j
  h_grid(1:h_i) = h_glgrid(1:gl_count:j,jp)
  t_grid(1:h_i) = t_glgrid(1:gl_count:j,jp)
  h_grid(h_i+1:klo) = h_glgrid(gl_count,jp)
  t_grid(h_i+1:klo) = t_glgrid(gl_count,jp)
!
!  Set tan_hts_raw to be a REAL SUBSET of h_grid, so NO height interpolation
!  for abs_cs will occure, also define tan_press as a TRUE subset of z_grid:

  klo = -1
  DO i = si+1, no_tan_hts
    h = tan_hts_raw(i)
    CALL hunt(h,h_grid,n_lvls,klo,khi)
    IF(ABS(h-h_grid(khi)) < ABS(h-h_grid(klo))) klo = khi
    tan_hts_raw(i) = h_grid(klo)
!   tan_press(i) = z_grid(klo)
    z = z_grid(klo)
    h = 12.0 * abs(z)
    khi = int(h+0.1)
    v = Real(khi)/12.0
    if(z < 0.0) v = -v
    if(abs(z-v).ge.1.0e-5) v = z
    tan_press(i) = v
  END DO

!  Make sure tan_hts array is stricktly monotonically increasing:

  j = no_tan_hts
  no_tan_hts = 1
  v = tan_hts_raw(1)
  DO i = 2, j
    z = tan_press(i)
    h = tan_hts_raw(i)
    IF(h > v) THEN
      v = h
      no_tan_hts = no_tan_hts + 1
      tan_press(no_tan_hts) = z
      tan_hts_raw(no_tan_hts) = h
    END IF
  END DO

! Interpolate the hydrostatic grid for conv. grid pressures (tan_press)
! for the values BELOW Earth surface only:

  CALL get_pressures('h',h_grid,t_grid,z_grid,n_lvls,tan_hts_raw, &
                     tan_press,si,ier)
  IF(ier /= 0) RETURN

! Interpolate the hydrostatic grid for convolution grid temperatures

  ier = 0
  h_i = 1
  z = z_grid(h_i)
  t = t_grid(h_i)
  DO WHILE(tan_press(h_i) < z)
    tan_temp(h_i) = t
    h_i = h_i + 1
  END DO

  CALL lintrp(z_grid,tan_press(h_i:),t_grid,tan_temp(h_i:), &
              n_lvls,no_tan_hts-h_i+1)

! Interpolate the dh_dt into the tan_press grid:

  DO khi = 1, no_t
    i = -ng
    h_i = 0
    j = ng + 1
    DO WHILE(i < gl_count)
      i = i + j
      h_i = h_i + 1
      dummy(h_i) = dh_dt_glgrid(i,jp,khi)
    END DO
    CALL lintrp(z_grid,tan_press(si:),dummy,tan_dh_dt(si:,khi), &
                h_i,no_tan_hts-si+1)
    tan_dh_dt(1:si-1,khi) = 0.0
    tan_dh_dt(h_i+1:,khi) = dh_dt_glgrid(gl_count,jp,khi)
  END DO

! Define the final convolution heights array.

  DO i = 1, no_tan_hts
    tan_hts(i) = tan_hts_raw(i)
  END DO

! Compute all the various integration paths according to tanget heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  DO k = 1, no_tan_hts
    h = tan_hts(k)
    CALL vert_to_path(n_lvls,Ng,Npath,gl_count,no_phi_t,no_t, &
         h,z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,            &
         dh_dt_glgrid,t_phi_basis,zpath,hpath,tpath,ppath,    &
         dhdzp,dhdtp,klo,khi,Ier)
    IF(ier /= 0) RETURN
    ALLOCATE(z_path(k)%values(khi), h_path(k)%values(khi),   &
             t_path(k)%values(khi), phi_path(k)%values(khi), &
             dhdz_path(k)%values(khi), STAT=j)
    IF(j == 0) ALLOCATE(dh_dt_path(k)%values(khi,no_phi_t,no_t),STAT=j)
    IF(j /= 0) THEN
      ier = j
      PRINT *,'** Error: ALLOCATION error in routine: geo_grids ..'
      PRINT *,'   STAT =',ier
      RETURN
    ENDIF
    ndx_path(k)%break_point_index = klo
    ndx_path(k)%total_number_of_elements = khi
    z_path(k)%values(1:khi) = zpath(1:khi)
    t_path(k)%values(1:khi) = tpath(1:khi)
    h_path(k)%values(1:khi) = hpath(1:khi)
    phi_path(k)%values(1:khi) = ppath(1:khi)
    dhdz_path(k)%values(1:khi) = dhdzp(1:khi)
    dh_dt_path(k)%values(1:khi,1:no_phi_t,1:no_t) = &
   &                              dhdtp(1:khi,1:no_phi_t,1:no_t)
  END DO

  RETURN
END SUBROUTINE geo_grids

!----------------------------------------------------------------------
!  The 2 dimensional Hydrostatic integrator

SUBROUTINE get_h_dhdt(t_profile,t_basis,t_phi_basis,zeta,phi,zref, &
           href,geoc_lat,reff,g,const,nt,no_t_phi,h,dhdt)

!  ===============================================================
!  Declaration of variables for sub-program: get_h_dhdt
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: nt, no_t_phi

Real(r8), INTENT(IN) :: geoc_lat, href, zref, zeta, phi, t_profile(:,:), &
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

Real(r8), SAVE :: gr2, c, prevlat = -500.0

Real(r8) :: pm(50) = 0.0
Real(r8) :: etapv(15) = 0.0

! Begin the code here

  IF(geoc_lat /= prevlat) THEN
    prevlat = geoc_lat
    CALL get_g_reff(geoc_lat,href,g,reff)
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

SUBROUTINE get_g_reff(geoc_lat,href,g,reff)

!  ===============================================================
!  Declaration of variables for sub-program: get_g_reff
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Real(r8), INTENT(IN) :: geoc_lat, href

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
!  ----------------
!  Local variables:
!  ----------------
Real(r8) :: cr, crr, a

! Begin code:

  cr   = DCOS(2.0D0 * geoc_lat)
  g    = g0 * (1.0D0 - cr * (g1 - cr * g2))       ! Modified G
  crr  = 2.0D0 * cr * cr - 1.0D0                  ! Cos(4*geoc_lat)
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
Real(r8) :: sgn, z_0, z_i, h, l, q

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
    q = m0 * (zetabase(iq) - zetabase(iq-1))
    v = (h - l) * (0.5 * (h + l) - zetabase(iq-1)) / q

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
    q = m0 * (zetabase(iq+1) - zetabase(iq))
    v = sgn * (v + (h - l) * (zetabase(iq+1) - 0.5 * (h + l)) / q)

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
