module VERT_TO_PATH_M
  use MLSCommon, only: I4, R4, R8
  use D_LINTRP_M, only: LINTRP
  use I_HUNT_M, only: HUNT
  use ELLIPSE_SW_M, only: H_TO_S_PHI
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE_M, only: ELLIPSE
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
!  (2*(Ng+1)*N_lvls) points per array.

! *** NOTE: This routine is using The Equivalent Circle concept

SUBROUTINE vert_to_path(elvar,n_lvls,Ng,ngt,gl_count,no_phi_t,no_t,htan,     &
           z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,t_phi_basis,z_path,  &
           h_path,t_path,phi_path,dhdz_path,phi_eta,brkpt,totnp,Ier)

!  ===============================================================
!  Declaration of variables for sub-program: vert_to_path
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls,Ng,gl_count,no_phi_t,no_t,ngt

Integer(i4), INTENT(OUT) :: totnp,brkpt,Ier

Type(ELLIPSE), intent(in out) :: elvar

Real(r8), INTENT(OUT) :: h_path(:), z_path(:), t_path(:), phi_path(:), &
                         dhdz_path(:), phi_eta(:,:)

Real(r8), INTENT(IN) :: htan, z_glgrid(:)
Real(r8), INTENT(IN) :: h_glgrid(:,:), t_glgrid(:,:), t_phi_basis(:), &
                        dhdz_glgrid(:,:)

!  ----------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, l, m, n, jp, n_d, npp, ibrk, Ngp1, no_iter

Real(r8) :: h, s, r, dz, rs, phi, rss, dhdz, prev_h

Integer(i4), ALLOCATABLE, DIMENSION(:) :: cndx

Real(r8), ALLOCATABLE, DIMENSION(:) :: dum_z, dum_h, dum_phi

Real(r8), ALLOCATABLE, DIMENSION(:,:) :: h_a

! Allocate enough space for the various (temporary) arrys we are going
! to use...

  Ier = 0
  Ngp1 = Ng + 1

  DEALLOCATE(cndx, dum_z, dum_h, dum_phi, STAT=i)
  ALLOCATE(cndx(ngt), dum_z(ngt), dum_h(ngt), dum_phi(ngt), &
 &         STAT = ier)
  IF(ier /= 0) THEN
    Ier = 1
    Print *,'** Error: ALLOCATION error in VERT_TO_PATH routine ..'
    GOTO 99
  endif

  DEALLOCATE(h_a, STAT=i)
  ALLOCATE(h_a(ngt,no_phi_t),STAT=ier)
  IF(ier /= 0) THEN
    Ier = 1
    Print *,'** Error: ALLOCATION error in VERT_TO_PATH routine ..'
    GOTO 99
  endif

! Initialze all arrays:

  cndx = 0
  dum_z = 0.0
  dum_h = 0.0
  dum_phi = 0.0
  DO j = 1, no_phi_t
    h_a(1:ngt,j) = 0.0
    phi_eta(1:ngt,j) = 0.0
  END DO

  r = -999.99
  dz = r / 57.2958
  z_path(1:ngt) = r
  h_path(1:ngt) = r
  t_path(1:ngt) = r
  phi_path(1:ngt) = dz
  dhdz_path(1:ngt) = 0.0

! Define the various ELLIPSE variables needed for computations:

  elvar%ht = htan
  elvar%earthx = (htan < -0.01)

  IF(htan < -0.01) THEN
    elvar%rr = (elvar%roc + elvar%ht) / elvar%roc
    CALL earth_intersection(elvar,rs)
  ELSE
    elvar%rr = 0.0D0
    elvar%phi_s = elvar%phi_tan
    elvar%nphi_s = elvar%nphi_tan
  END IF

! Define the index points of the tangent locations:

  k = 2 * N_lvls

! Store the initial guess heights in dum_h. This estimate is the h_glgrid
! at the "Center Phi", Also compute the Path zeta on the GL grid:

  npp = 0
  l = gl_count + 1
  jp = (no_phi_t + 1) / 2

  DO
    DO n = 1, Ngp1
      l = l - 1
      npp = npp + 1
      dum_z(npp) = z_glgrid(l)
      dum_h(npp) = h_glgrid(l,jp)
    END DO
    h = h_glgrid(l-1,jp)
    if(h <= htan  .OR. l-Ngp1 < 1) EXIT
  END DO

  l = l - 1
  npp = npp + 1
  dum_z(npp) = z_glgrid(l)
  dum_h(npp) = h_glgrid(l,jp)

  l = l - 1
  ibrk = npp + 1

  DO
    DO n = 1, Ngp1
      l = l + 1
      npp = npp + 1
      dum_z(npp) = z_glgrid(l)
      dum_h(npp) = h_glgrid(l,jp)
    END DO
    IF(l+Ngp1 > gl_count) EXIT
  END DO

  l = l + 1
  npp = npp + 1
  dum_z(npp) = z_glgrid(l)
  dum_h(npp) = h_glgrid(l,jp)

!  Cast the h_glgrid onto the path (using liner interpolation)

  DO m = 1, no_phi_t
    CALL lintrp(z_glgrid,dum_z,h_glgrid(1:,m),h_a(1:,m),gl_count,npp)
  END DO

  n_d = ngt
  no_iter = 0

10 no_iter = no_iter + 1

!  Compute the Phi that goes with the Heights along the path:
!  Right hand side ray:

  l = -1
  elvar%ps = -1.0D0
  DO i = 1, ibrk-1
    h = dum_h(i)
    phi = dum_phi(i)
    IF(n_d == ngt) THEN
      CALL H_TO_S_PHI(elvar,h,s,phi)
    ELSE
      CALL hunt(i,cndx,n_d,l,j)
      IF(cndx(l) == i .OR. cndx(j) == i) &
              CALL H_TO_S_PHI(elvar,h,s,phi)
    END IF
    dum_phi(i) = phi
  END DO

!  Left hand side ray:

  l = -1
  elvar%ps = 1.0D0
  DO i = ibrk, npp
    h = dum_h(i)
    phi = dum_phi(i)
    IF(n_d == ngt) THEN
      CALL H_TO_S_PHI(elvar,h,s,phi)
    ELSE
      CALL hunt(i,cndx,n_d,l,j)
      IF(cndx(l) == i .OR. cndx(j) == i) &
               CALL H_TO_S_PHI(elvar,h,s,phi)
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
    prev_h = dum_h(i)
    dum_h(i) = SUM(h_a(i,:)*phi_eta(i,:))
    r = ABS(dum_h(i) - prev_h)
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
      elvar%ps = -1.0D0
      dz = dum_z(j) - dum_z(k)
      dum_h(j) = dum_h(k) + dhdz * dz
      h = dum_h(j)
      IF(j >= ibrk) elvar%ps = 1.0D0
      CALL H_TO_S_PHI(elvar,h,s,phi)
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

  s = 1.0e10
  brkpt = -1
  totnp = npp

!  Also, compute the break points for the tangent layer for this tanget height

  DO i = 1, npp
    t_path(i) = SUM(h_a(i,:)*phi_eta(i,:))
    z_path(i) = dum_z(i)
    h_path(i) = dum_h(i)
    phi_path(i) = dum_phi(i)
    r = abs(h_path(i)-htan)
    if(r < s) then
      s = r
      brkpt = i
    endif
  END DO

! Second, compute the path dh_dz:
! Cast the dhdz_glgrid onto the path (using liner interpolation)

  DO m = 1, no_phi_t
    CALL lintrp(z_glgrid,dum_z,dhdz_glgrid(1:,m),h_a(1:,m),gl_count,npp)
  END DO

  DO i = 1, npp
    dhdz_path(i) = SUM(h_a(i,:)*phi_eta(i,:))
  END DO

 99 DEALLOCATE(h_a, STAT=i)
    DEALLOCATE(cndx, dum_z, dum_h, dum_phi, STAT=i)

  RETURN
END SUBROUTINE vert_to_path
end module VERT_TO_PATH_M
! $Log$
! Revision 1.6  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.5  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.4  2001/03/20 11:03:16  zvi
! Fixing code for "real" data run, increase dim. etc.
!
! Revision 1.3  2001/03/09 00:40:32  zvi
! Correcting an error in HUNT routine
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:18  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
