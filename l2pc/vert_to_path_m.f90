!
module VERT_TO_PATH_M
  use L2PCDim, only: NLVL, NSPS, N2LVL
  use MLSCommon, only: I4, R4, R8
  use D_LINTRP_M, only: LINTRP
  use I_HUNT_M, only: HUNT
  use ELLIPSE_SW_M, only: H_TO_S_PHI
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
 &    HT2, RR, PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
!    This subroutine computes all the (z,t,h,phi,dh_dz,dh_dt) along a
! GIVEN ray path, which are tanget dependent.
!    The grid points are the Gauss-Legendre points, i.e. there are: gl_count
!  (2*Ng*N_lvls) points per array.

! *** NOTE: This routine is using The Equivalent Circle concept

SUBROUTINE vert_to_path(n_lvls,ng,npath,gl_count,no_phi_t,no_t,    &
           jtan,htan,z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,       &
           dh_dt_glgrid,t_phi_basis,z_path,h_path,t_path,phi_path, &
           dhdz_path,dh_dt_path,path_brkpt,rad_cur,phi_lat)

!  ===============================================================
!  Declaration of variables for sub-program: vert_to_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls,ng,gl_count,no_phi_t,no_t,jtan,npath

Integer(i4), INTENT(OUT) :: path_brkpt(:)

Real(r4), INTENT(OUT) :: dh_dt_path(:,:,:), dhdz_path(:)

Real(r8), INTENT(IN) :: rad_cur, phi_lat, htan, z_glgrid(:), h_glgrid(:,:), &
      t_glgrid(:,:), t_phi_basis(:), dhdz_glgrid(:,:), dh_dt_glgrid(:,:,:)

Real(r8), INTENT(OUT) :: z_path(:), t_path(:), h_path(:), phi_path(:)

!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: kgp = 1800
!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, l, m, n, jp, n_d, ntr, npp, ntl, brkp, ngp1, &
               istop, istart, no_iter, cndx(kgp)

Real(r8) :: h, s, r, dz, rs, phi, rss, sum, dhdz, prev_h, h_a(kgp,5),    &
          dum_z(kgp), dum_h(kgp), dum_t(kgp), dum_dh(kgp), dum_phi(kgp), &
          phi_eta(kgp,5)

! Store the initial guess heights in dum_h. This estimate is the h_glgrid
! at the "Center Phi" :

  n_d = kgp
  DO i = 1, n_d
    cndx(i) = i
    dum_z(i) = 0.0
    dum_h(i) = 0.0
    dum_t(i) = 0.0
    dum_dh(i) = 0.0
    dum_phi(i) = 0.0
    DO j = 1, 5
      h_a(i,j) = 0.0
      phi_eta(i,j) = 0.0
    END DO
  END DO

  r = -999.99
  dz = r / 57.3
  DO i = 1, npath
    h_path(i) = r
    z_path(i) = r
    t_path(i) = r
    dhdz_path(i) = r
    phi_path(i) = dz
  END DO

! Define the various COMMON variables needed for computations:

  ht = htan
  roc = rad_cur
  phi_tan = phi_lat
  earthx = (htan < -0.01)

  IF(htan < -0.01) THEN
    rr = (roc + ht) / roc
    CALL earth_intersection(rs)
  ELSE
    rr = 0.0D0
    phi_s = phi_tan
    nphi_s = nphi_tan
  END IF

! Define the index points of the tangent locations:

  k = 2 * n_lvls
  ntr = n_lvls - jtan + 1
  ntl = n_lvls + jtan

!  Compute the Path zeta on the GL grid:

  i = 0
  npp = 0
  ngp1 = ng + 1
  l = gl_count + 1
  jp = (no_phi_t+1)/2
  DO WHILE(i < ntr)
    i = i + 1
    DO n = 1, ngp1
      l = l - 1
      npp = npp + 1
      dum_z(npp) = z_glgrid(l)
      dum_h(npp) = h_glgrid(l,jp)
    END DO
    IF(l-ngp1 < 1) i = ntr+3
  END DO

  l = l - 1
  npp = npp + 1
  dum_z(npp) = z_glgrid(l)
  dum_h(npp) = h_glgrid(l,jp)

  l = l - 1
  i = ntl - 1
  brkp = npp + 1
  DO WHILE(i < k)
    i = i + 1
    DO n = 1, ngp1
      l = l + 1
      npp = npp + 1
      dum_z(npp) = z_glgrid(l)
      dum_h(npp) = h_glgrid(l,jp)
    END DO
    IF(l+ngp1 > gl_count) i = k+3
  END DO

  l = l + 1
  npp = npp + 1
  dum_z(npp) = z_glgrid(l)
  dum_h(npp) = h_glgrid(l,jp)

!  Cast the h_glgrid onto the path (using liner interpolation)

  DO m = 1, no_phi_t
    CALL lintrp(z_glgrid,dum_z,h_glgrid(1:,m),h_a(1:,m),gl_count,npp)
  END DO

  no_iter = 0
  10   no_iter = no_iter + 1

!  Compute the Phi that goes with the Heights along the path:
!  Right hand side ray:

  l = -1
  ps = -1.0D0
  DO i = 1, brkp-1
    h = dum_h(i)
    phi = dum_phi(i)
    IF(n_d == kgp) THEN
      CALL H_TO_S_PHI(h,s,phi)
    ELSE
      CALL hunt(i,cndx,n_d,l,j)
      IF(cndx(l) == i.OR.cndx(j) == i) CALL H_TO_S_PHI(h,s,phi)
    END IF
    dum_phi(i) = phi
  END DO

!  Left hand side ray:

  ps = 1.0D0
  DO i = brkp, npp
    h = dum_h(i)
    phi = dum_phi(i)
    IF(n_d == kgp) THEN
      CALL H_TO_S_PHI(h,s,phi)
    ELSE
      CALL hunt(i,cndx,n_d,l,j)
      IF(cndx(l) == i.OR.cndx(j) == i) CALL H_TO_S_PHI(h,s,phi)
    END IF
    dum_phi(i) = phi
  END DO

!  Compute the Phi_Eta matrix based on the "Path" Phi's

  DO i = 1, npp
    r = dum_phi(i)
    DO m = 1, no_phi_t
      CALL get_one_eta(r,t_phi_basis,no_phi_t,m,phi_eta(i,m))
    END DO
  END DO

!  Re-Compute the estimate H:

  n_d = 0
  rss = 0.0D0
  DO i = 1, npp
    sum = 0.0
    prev_h = dum_h(i)
    DO m = 1, no_phi_t
      sum = sum + h_a(i,m) * phi_eta(i,m)
    END DO
    dum_h(i) = sum
    r = ABS(sum - prev_h)
    IF(r > 0.01) THEN
      n_d = n_d + 1
      cndx(n_d) = i
      rss = rss + r * r
    END IF
  END DO

! **** ITERATE

  IF(n_d > 0.AND.no_iter < 10) GO TO 10

! For indecies which did not converged within 10 iterations, use liner
! interpolation on H and then recompute Phi and Phi_Eta for those indecies.

  IF(n_d > 1) THEN
    k = MAX(1,cndx(1)-1)
    l = MIN(npp,cndx(n_d)+1)
    dz = dum_z(l) - dum_z(k)
    dhdz = (dum_h(l) - dum_h(k)) / dz
    DO i = 1, n_d
      j = cndx(i)
      ps = -1.0D0
      dz = dum_z(j) - dum_z(k)
      dum_h(j) = dum_h(k) + dhdz * dz
      h = dum_h(j)
      IF(j >= brkp) ps = 1.0D0
      CALL H_TO_S_PHI(h,s,phi)
      dum_phi(j) = phi
      DO m = 1, no_phi_t
        CALL get_one_eta(dum_phi(j),t_phi_basis,no_phi_t,m, phi_eta(j,m))
      END DO
    END DO
  END IF

! **** CONVERGED !!
! Now - Compute the path Temperature, dh_dz  and dH_dTlm

! First, compute the path Temperature:
! Cast the t_glgrid onto the path (using liner interpolation)

  DO m = 1, no_phi_t
    CALL lintrp(z_glgrid,dum_z,t_glgrid(1:,m),h_a(1:,m),gl_count,npp)
  END DO

  DO i = 1, npp
    sum = 0.0
    DO m = 1, no_phi_t
      sum = sum + h_a(i,m) * phi_eta(i,m)
    END DO
    dum_t(i) = sum
  END DO

! Second, compute the path dh_dz:
! Cast the dhdz_glgrid onto the path (using liner interpolation)

  DO m = 1, no_phi_t
    CALL lintrp(z_glgrid,dum_z,dhdz_glgrid(1:,m),h_a(1:,m),gl_count,npp)
  END DO

  DO i = 1, npp
    sum = 0.0
    DO m = 1, no_phi_t
      sum = sum + h_a(i,m) * phi_eta(i,m)
    END DO
    dum_dh(i) = sum
  END DO

!  Compute the break points for the tangent layer for this tanget height

  i = 0
  j = 0
  istop = -1
  istart = -1
  h = htan + 1.0
  DO WHILE(i < brkp-1.AND.h >= htan)
    i = i + 1
    h = dum_h(i)
    IF(h >= htan) THEN
      j = j + 1
      istop = i
      h_path(j) = h
      z_path(j) = dum_z(i)
      t_path(j) = dum_t(i)
      phi_path(j) = dum_phi(i)
      dhdz_path(j) = dum_dh(i)
    END IF
  END DO

  path_brkpt(1) = j

  m = ntr
  DO WHILE(m < ntl-2)
    m = m + 1
    j = j + ngp1
  END DO

  path_brkpt(2) = j + 1

  i = brkp-1
  DO WHILE(i < npp)
    i = i + 1
    j = j + 1
    h = dum_h(i)
    IF(h < htan) THEN
      istart = i + 1
      path_brkpt(2) = j + 1
    ELSE
      h_path(j) = dum_h(i)
      z_path(j) = dum_z(i)
      t_path(j) = dum_t(i)
      phi_path(j) = dum_phi(i)
      dhdz_path(j) = dum_dh(i)
    END IF
  END DO

  path_brkpt(3) = j

! Fourth: compute the path dH_dTlm:
! Cast the dh_dt_glgrid onto the path (using liner interpolation)

  DO j = 1, no_t
    DO m = 1, no_phi_t
      CALL lintrp(z_glgrid,dum_z,dh_dt_glgrid(1:,m,j),dum_h,gl_count,npp)
      DO i = 1, npp
        dum_t(i) = dum_h(i) * phi_eta(i,m)
      END DO
      k = 0
      DO i = 1, istop
        k = k + 1
        dh_dt_path(k,m,j) = dum_t(i)
      END DO
      DO WHILE(k < path_brkpt(2)-1)
        k = k + 1
        dh_dt_path(k,m,j) = -999.99
      END DO
      i = istart - 1
      DO WHILE(i < npp)
        i = i + 1
        k = k + 1
        dh_dt_path(k,m,j) = dum_t(i)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE vert_to_path
end module VERT_TO_PATH_M
! $Log$
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
