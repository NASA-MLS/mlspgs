! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module HYDROSTATIC_MODEL_M
  use MLSCommon, only: I4, R8
  use GLNP, only: NG, GX
  use D_HUNT_M, only: HUNT
  use D_LINTRP_M, only: LINTRP
  use GET_ETA_M, only: GET_ETA
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

SUBROUTINE hydrostatic_model(si, no_phi_t, geoc_lat, Href, &
           Zref, z_grid, t_z_basis, t_coeff, z_glgrid, h_glgrid,     &
           t_glgrid, dhdz_glgrid, dh_dt_glgrid, tan_press, tan_hts,  &
           tan_temp, tan_dh_dt, gl_count, Ier)

!  ===============================================================
!  Declaration of variables for sub-program: hydrostatic_model
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------

Integer(i4), INTENT(IN) :: si, no_phi_t

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
Integer(i4) :: cast(Size(z_grid))
Integer(i4) :: h_i,i,j,k,l,m,jj,cnt,Ngp1,no_tan_hts,n_lvls,no_t,NLm1

Real(r8) :: zGx(Ng+1)
Real(r8) :: G,Reff,const,q,h,z,t,v,z1,z2
Real(r8) :: h_grid(Size(z_grid)),t_grid(Size(z_grid)),dhdt(Size(t_z_basis))

Real(r8), DIMENSION(:), ALLOCATABLE :: xm, ym
Real(r8), DIMENSION(:,:), ALLOCATABLE :: Eta

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
    cast(i) = 1 + (j - 1) * Ngp1
  end do
!
  NLm1 = N_Lvls - 1
  ALLOCATE(xm(NLm1),ym(NLm1),STAT=ier)
  if(ier /= 0) then
    Print *,'** Error in hydrostatic_model routine ..'
    Print *,'   Allocation error for vaectors: xm and ym, STAT=',ier
    Return
  endif
!
! From the selected integration grid pressures define the GL pressure
! grid:

  zGx(1) = -1.0d0
  zGx(2:Ngp1) = Gx(1:Ng)
!
  xm(1:NLm1) = 0.5d0 * ( z_grid(2:N_Lvls) + z_grid(1:NLm1) )
  ym(1:NLm1) = 0.5d0 * ( z_grid(2:N_Lvls) - z_grid(1:NLm1) )

  k = NLm1 * Ngp1
  z_glgrid(1:k) = RESHAPE ( (SPREAD(xm,1,Ngp1) +                 &
                   &    SPREAD(ym,1,Ngp1) * SPREAD(zGx,2,NLm1)), &
                   &   (/k/))
  gl_count = k + 1
  z_glgrid(gl_count) = z_grid(N_Lvls)
!
  DEALLOCATE(xm,ym,STAT=i)
!
! *** Create the Temperature on the GL grid by linear interpolation via
!     the get_eta procedure (Bill's request)
!
  ALLOCATE(Eta(gl_count,no_t),STAT=ier)
  if(ier /= 0) then
    Print *,'** Error in hydrostatic_model routine ..'
    Print *,'   Allocation error for matrix: Eta, STAT=',ier
    Return
  endif

  Call GET_ETA(z_glgrid,t_z_basis,gl_count,no_t,Eta)

  t_glgrid(1:gl_count,:) = MATMUL(Eta,t_coeff)

! Get Hydrostatically balanced group ([zth]_glgrid, dhdz_glgrid and
! dh_dt_glgrid on the GL grid:

  DO l = 1, no_phi_t
    Call get_h_dhdt(t_coeff(1:,l),t_z_basis,Zref,         &
      &  Href(l),geoc_lat(l),z_glgrid,t_glgrid(1:,l),Eta, &
      &  h_glgrid(1:,l),dhdz_glgrid(:,l),dh_dt_glgrid(:,l,:))
  END DO

  DEALLOCATE(Eta,STAT=i)
!
! Define tan_hts as a TRUE subset of h_grid for each mmaf:
!
  k = no_tan_hts + si - 1
  DO l = 1, no_phi_t
    tan_hts(si:k,l) = h_glgrid((/(cast(i),i=1,no_tan_hts)/),l)
   tan_temp(si:k,l) = t_glgrid((/(cast(i),i=1,no_tan_hts)/),l)
   if(si > 1) tan_temp(1:si-1,l) = t_glgrid(1,l)
  END DO
  no_tan_hts = k

  if(si > 1) then

! Interpolate the hydrostatic grid for conv. grid heights
! for the values BELOW Earth surface only:

    h_grid(1:) = 0.0
    t_grid(1:) = 0.0
    DO l = 1, no_phi_t
      h_grid(1:) = h_glgrid(1:gl_count:Ngp1,l)
      t_grid(1:) = t_glgrid(1:gl_count:Ngp1,l)
      CALL get_heights('h',h_grid,t_grid,z_grid,N_Lvls, &
                        tan_press,tan_hts(:,l),si-1,ier)
      IF(ier /= 0) RETURN
    end do

  endif
!
! Cast the dh_dt into the tan_press grid:
!
  k = no_tan_hts-si+1
  DO m = 1, no_t
    DO l = 1, no_phi_t
      tan_dh_dt(si:no_tan_hts,l,m)=dh_dt_glgrid((/(cast(i),i=1,k)/),l,m)
      if(si > 1) tan_dh_dt(1:si-1,l,m) = 0.0
    END DO
  END DO

  RETURN
END SUBROUTINE hydrostatic_model

!----------------------------------------------------------------------
!  The 2 dimensional Hydrostatic integrator

SUBROUTINE get_h_dhdt(t_profile,t_basis,Zref,Href,geoc_lat,z_grid, &
                   &  t_grid,Pqi,h_grid,dhdz,dhdt)
!  ===============================================================
!  Declaration of variables for sub-program: get_h_dhdt
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------

Real(r8), INTENT(IN) :: t_profile(:), t_basis(:), z_grid(:), t_grid(:)
Real(r8), INTENT(IN) :: geoc_lat, Zref, Href

Real(r8), INTENT(OUT) :: h_grid(:), dhdz(:)
Real(r8), INTENT(OUT) :: dhdt(:,:), Pqi(:,:)

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Real(r8), Parameter :: m0 = 28.964125d0
Real(r8), Parameter :: boltzxln10 = 19.14486942d0
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: j, iter, nt, nz

Real(r8) :: G, Reff, z_prev, z_surf, h, f_prev, f_surf, v, s, dh_dz_s

Real(r8) :: pm(Size(t_basis))

! Begin the code here

!  Compute Surface zeta, given Href & Zref

  nt = Size(t_basis)
  nz = Size(z_grid)
!
  s = 0.0d0
  Pqi = 0.0
  pm(1:nt) = 0.0
  CALL get_g_reff(geoc_lat,s,G,Reff)

  z_surf = z_grid(1)
  z_prev = z_surf + 1.0
  f_surf = 1.0
  iter = 0

  DO

    iter = iter + 1
    IF(ABS(z_surf-z_prev) < 1.0e-5 .OR. iter > 10) EXIT

    DO j = 1, nt
      Call pq_ana(z_surf,Zref,t_basis,j,nt,pm(j))
    END DO

    s = SUM(t_profile(1:nt)*pm(1:nt))
    h = boltzxln10 * s / G
    f_surf = h - Href
    IF(ABS(f_surf) < 1.0e-4) EXIT

    if(iter == 1) then
      z_prev = z_surf
      f_prev = f_surf
      z_surf = z_prev + 0.05
    else
      dh_dz_s = (f_surf - f_prev) / (z_surf - z_prev)
      z_prev = z_surf
      f_prev = f_surf
      z_surf = z_prev - f_prev / dh_dz_s
    endif

  END DO

! compute the piq integrals relative to the surface

  Call do_pqi(z_surf,z_grid,t_basis,nz,nt,Pqi)

! compute the height vector

  v = G * Reff
  s = boltzxln10 / (v * Reff)
  h_grid(1:nz) = boltzxln10*MATMUL(Pqi,t_profile)    ! GPH * G
  h_grid(1:nz) = Reff * h_grid(1:nz) / (v - h_grid(1:nz))
  dhdz(1:nz) = s * ((h_grid(1:nz) + Reff)**2)
  dhdt(1:nz,:) = SPREAD(dhdz(1:nz),2,nt) * Pqi
  dhdz(1:nz) = dhdz(1:nz) * t_grid(1:nz) / m0

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

  if(Href == 0.0) Return

! Make approriate correction for Href not being the surface

  a = (1.0D0 + Href / Reff)
  G = G / (a * a)
  Reff = Reff + Href

  RETURN
END SUBROUTINE get_g_reff

!---------------------------------------------------------------------

SUBROUTINE do_pqi(Zsurf,z_grid,zetabase,nz,nt,pqi)

!  ===============================================================
!  Declaration of variables for sub-program: do_pqi
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: nz,nt

Real(r8), INTENT(IN) :: z_grid(:),zetabase(:), Zsurf
Real(r8), INTENT(OUT) :: pqi(:,:)

Integer(i4) :: i,j

    DO j = 1, nt
      DO i = 1, nz
        Call pq_ana(Zsurf,z_grid(i),zetabase,j,nt,pqi(i,j))
      END DO
    END DO

  RETURN
END SUBROUTINE do_pqi

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
! Revision 1.18  2001/06/07 23:39:31  pwagner
! Added Copyright statement
!
! Revision 1.17  2001/05/11 22:18:03  livesey
! Changed an old fixed dimension variable.
!
! Revision 1.16  2001/04/23 21:43:28  zvi
! Introducing no_phi_t etc.
!
! Revision 1.15  2001/04/09 20:52:07  zvi
! Debugging Derivatives version
!
! Revision 1.14  2001/04/06 23:58:46  zvi
! Fix a bug LF95 was compalining about (tan_hts ..)
!
! Revision 1.13  2001/04/06 20:54:21  zvi
! Fix a small bug concerning initialization of JLO in HUNT
!
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
