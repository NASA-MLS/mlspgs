module HYDROSTATIC_MODEL_M
  use MLSCommon, only: I4, R8
  use GL6P, only: NG, GX
  use D_HUNT_M, only: HUNT
  use D_LINTRP_M, only: LINTRP
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use HYDROSTATIC_INTRP, only: GET_HEIGHTS, GET_PRESSURES

  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------

SUBROUTINE hydrostatic_model(si,no_mmaf, geoc_lat, Href, &
           Zref, z_grid, t_z_basis, t_coeff, z_glgrid, h_glgrid,     &
           t_glgrid, dhdz_glgrid, dh_dt_glgrid, tan_press, tan_hts,  &
           tan_temp, tan_dh_dt, gl_count, Ier)

!  ===============================================================
!  Declaration of variables for sub-program: hydrostatic_model
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------

Integer(i4), INTENT(IN) :: si,no_mmaf

Integer(i4), INTENT(OUT) :: Ier, gl_count

Real(r8), INTENT(IN) :: geoc_lat(:), Zref, Href(:), z_grid(:)

Real(r8), INTENT(IN) :: t_z_basis(:)

Real(r8), INTENT(IN) :: t_coeff(:,:)

Real(r8), INTENT(OUT) :: z_glgrid(:), h_glgrid(:,:), t_glgrid(:,:), &
                         dhdz_glgrid(:,:), dh_dt_glgrid(:,:,:)

Real(r8), INTENT(IN) :: tan_press(:)
Real(r8), INTENT(OUT) :: tan_hts(:,:),  tan_temp(:,:)

Real(r8), INTENT(OUT) :: tan_dh_dt(:,:,:)

!  ----------------------
!  Local variables:
!  ----------------
Integer(i4) :: t_index(Size(z_grid))
Integer(i4) :: h_i,i,j,k,l,m,jj,cnt,Ngp1,no_tan_hts,n_lvls,no_t

Real(r8) :: G,Reff,const,q,h,z,t,v,z1,z2,xm,ym
Real(r8) :: h_grid(Size(z_grid)),t_grid(Size(z_grid)),dhdt(Size(t_z_basis))

! Begin the code here

  ier = 0
  Ngp1 = Ng + 1
  N_lvls = Size(z_grid)
  no_t = Size(t_z_basis)
  no_tan_hts = Size(tan_press) - si + 1 
!
  j = -1
  do i = 1, no_tan_hts
    z1 = tan_press(i+si-1)
    Call Hunt(z1,z_grid,n_lvls,j,k)
    if(abs(z1-z_grid(j)) > abs(z1-z_grid(k))) j=k
    if(abs(z1-z_grid(j)) > 1.0e-4) then
      Ier = 1
      Print *,'** Error in hydrostatic_model routine ..'
      Print *,'   Tanget array NOT a subset of Integration grid !'
      Print *,'   Nathaniel is a LIAR !'
      Return
    endif
    t_index(i) = j
  end do

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
    DO j = 1, Ng
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
  DO l = 1, no_mmaf
    do i = 1, gl_count
      t = 0.0
      z = z_glgrid(i)
      do h_i = 1, no_t
        CALL get_one_eta(z,t_z_basis,no_t,h_i,v)
        t = t + t_coeff(h_i,l) * v
      end do
      t_glgrid(i,l) = t
    end do
  END DO

! Get Hydrostatically balanced group ([zth]_glgrid, dhdz_glgrid and
! dh_dt_glgrid on the GL grid:

  DO h_i = 1, gl_count
    z = z_glgrid(h_i)
    j = (h_i + Ng) / Ngp1
    DO l = 1, no_mmaf
      CALL get_h_dhdt(t_coeff(1:,l),t_z_basis,z,Zref,Href(l), &
           geoc_lat(l),Reff,G,const,no_t,h,dhdt)
      h_glgrid(h_i,l) = h
      q = h + Reff
      dhdz_glgrid(h_i,l) = q * q * const * t_glgrid(h_i,l)
      dh_dt_glgrid(h_i,l,1:no_t) = dhdt(1:no_t)
    END DO
  END DO
!
! Define tan_hts as a TRUE subset of h_grid for each mmaf:
!
  h_grid(1:) = 0.0
  jj = Size(z_grid)
  cnt = (gl_count + Ng) / Ngp1
  do l = 1, no_mmaf
    tan_temp(1:si-1,l) = t_glgrid(1,l)
    h_grid(1:cnt) = h_glgrid(1:gl_count:Ngp1,l)
    if(cnt < jj) h_grid(cnt+1:jj) = h_glgrid(gl_count,l)
    t_grid(1:cnt) = t_glgrid(1:gl_count:Ngp1,l)
    if(cnt < jj) t_grid(cnt+1:jj) = t_glgrid(gl_count,l)
    do i = 1, no_tan_hts
      j = t_index(i)
      tan_hts(si+i-1,l) = h_grid(j)
      tan_temp(si+i-1,l) = t_grid(j)
    end do
  end do
!
  no_tan_hts = no_tan_hts + si - 1

! Interpolate the hydrostatic grid for conv. grid heights 
! for the values BELOW Earth surface only:

  h_grid(1:) = 0.0
  t_grid(1:) = 0.0
  do l = 1, no_mmaf
    h_grid(1:cnt) = h_glgrid(1:gl_count:Ngp1,l)
    t_grid(1:cnt) = t_glgrid(1:gl_count:Ngp1,l)
    if(cnt < jj) h_grid(cnt+1:jj) = h_glgrid(gl_count,l)
    if(cnt < jj) t_grid(cnt+1:jj) = t_glgrid(gl_count,l)
    CALL get_heights('h',h_grid,t_grid,z_grid,n_lvls,tan_press, &
   &                  tan_hts(:,l),si-1,ier)
    IF(ier /= 0) RETURN
  end do

  k = no_mmaf / 2
  h_grid(1:) = 0.0
  t_grid(1:) = 0.0
  h_grid(1:cnt) = h_glgrid(1:gl_count:Ngp1,k)
  t_grid(1:cnt) = t_glgrid(1:gl_count:Ngp1,k)
  if(cnt < jj) h_grid(cnt+1:jj) = h_glgrid(gl_count,k)
  if(cnt < jj) t_grid(cnt+1:jj) = t_glgrid(gl_count,k)

! Interpolate the dh_dt into the tan_press grid:

  t_grid(1:cnt) = 0.0
  DO l = 1, no_mmaf
    DO m = 1, no_t
      t_grid(1:cnt) = dh_dt_glgrid(1:gl_count:Ngp1,l,m)
      CALL lintrp(z_grid,tan_press(si:),t_grid,tan_dh_dt(si:,l,m), &
     &            cnt,no_tan_hts-si+1)
      tan_dh_dt(1:si-1,l,m) = 0.0
    END DO
  END DO

  RETURN
END SUBROUTINE hydrostatic_model

!----------------------------------------------------------------------
!  The 2 dimensional Hydrostatic integrator

SUBROUTINE get_h_dhdt(t_profile,t_basis,zeta,Zref,Href,geoc_lat, &
           Reff,G,const,nt,h,dhdt)

!  ===============================================================
!  Declaration of variables for sub-program: get_h_dhdt
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: nt

Real(r8), INTENT(IN) :: t_profile(:)
Real(r8), INTENT(IN) :: geoc_lat, Href, Zref, zeta, t_basis(:)

Real(r8), INTENT(OUT) :: const, h, G, Reff, dhdt(:)
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), Parameter :: m0 = 28.964125d0
Real(r8), Parameter :: boltzxln10 = 19.14486942d0
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: j

Real(r8) :: v, denom, sum

Real(r8), SAVE :: gr2, c, prevlat = -500.0

Real(r8) :: pm(50) = 0.0

! Begin the code here

  IF(geoc_lat /= prevlat) THEN
    prevlat = geoc_lat
    CALL get_g_reff(geoc_lat,Href,G,Reff)
    c = G * Reff
    gr2 = c * Reff
    const = boltzxln10 / (gr2 * m0)
  END IF

  dhdt(1:) = 0.0

  DO j = 1, nt
    Call pq_ana(Zref,zeta,t_basis,j,nt,v)
    pm(j) = v
  END DO

  sum = 0.0D0
  DO j = 1, nt
    sum = sum + t_profile(j) * pm(j)
  END DO

  IF(sum /= 0.0D0) THEN
    denom = c - boltzxln10 * sum
    h = gr2 / denom - Reff + Href
    v = gr2 * boltzxln10 / (denom * denom)
  ELSE
    h = Href
    v = boltzxln10 / G
  END IF

  dhdt(1:nt) = v * pm(1:nt)

  RETURN
END SUBROUTINE get_h_dhdt

!---------------------------------------------------------------------

SUBROUTINE get_g_reff(geoc_lat,Href,G,Reff)

!  ===============================================================
!  Declaration of variables for sub-program: get_g_reff
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Real(r8), INTENT(IN) :: geoc_lat, Href

Real(r8), INTENT(OUT) :: G, Reff
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), Parameter :: g0 = 9.80616d0
Real(r8), Parameter :: g1 = 2.6373d-3
Real(r8), Parameter :: g2 = 5.9d-6
Real(r8), Parameter :: r1 = 3.085462d-3
Real(r8), Parameter :: r2 = 2.27d-6
Real(r8), Parameter :: r3 = 2.0d-9
!  ----------------
!  Local variables:
!  ----------------
Real(r8) :: cr, crr, a

! Begin code:

  cr   = Cos(2.0D0 * geoc_lat)
  G    = g0 * (1.0D0 - cr * (g1 - cr * g2))       ! Modified G
  crr  = 2.0D0 * cr * cr - 1.0D0                  ! Cos(4*geoc_lat)
  Reff = 2.0D0 * g / (r1 + r2 * cr - r3 * crr)    ! Reff in Kilometers

! Make approriate correction for Href not being the surface

  a = (1.0D0 + Href / Reff)
  G = G / (a * a)
  Reff = Reff + Href

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
Real(r8), Parameter :: m0 = 28.964125d0
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

end module HYDROSTATIC_MODEL_M
! $Log$
! Revision 1.12  2001/03/29 23:58:48  livesey
! Turned phi_tan etc. into an array
!
! Revision 1.11  2001/03/29 02:27:04  zvi
! *** empty log message ***
!
! Revision 1.10  2001/03/29 02:09:28  zvi
! Fix an error
!
! Revision 1.9  2001/03/29 01:39:42  zvi
! Fixing an error in tan_hts computations
!
! Revision 1.8  2001/03/29 01:27:15  livesey
! Fixed bug with wrong intent for tan_press
!
! Revision 1.7  2001/03/29 01:21:50  zvi
! Interim version
!
! Revision 1.6  2001/03/28 23:50:11  zvi
! Tanget below surface are now in Zeta units..
!
! Revision 1.5  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
