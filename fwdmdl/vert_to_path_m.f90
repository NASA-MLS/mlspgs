! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module VERT_TO_PATH_M
  use Allocate_Deallocate, only: Allocate_test, Deallocate_Test
  use MLSCommon, only: I4, R4, R8
  use I_HUNT_M, only: HUNT
  use ELLIPSE_SW_M, only: H_TO_S_PHI
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use EARTH_INTERSECTION_M, only: EARTH_INTERSECTION
  use ELLIPSE_M, only: ELLIPSE

  Implicit NONE
  private
  public :: Vert_To_Path

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
!    This subroutine computes all the (z,t,h,phi,dh_dz,dh_dt) along a
! GIVEN ray path, which are tanget dependent.
!    The grid points are the Gauss-Legendre points, i.e. there are: gl_count
!  (2*(Ng+1)*N_lvls) points per array.

! *** NOTE: This routine is using The Equivalent Circle concept

  subroutine Vert_To_Path ( Elvar, N_lvls, Ng, Ngt, gl_count, WinSize, &
    &  MidWin, No_t, Htan, Z_glgrid, T_glgrid, H_glgrid, Dhdz_glgrid, &
    &  T_phi_basis, Z_path, H_path, T_path, Phi_path, Dhdz_path, Phi_eta, &
    &  Brkpt, Totnp, Ier )

    !  ===============================================================
    !  Declaration of variables for sub-program: vert_to_path
    !  ===============================================================
    !  ---------------------------
    !  calling sequence variables:
    !  ---------------------------
    integer, intent(in) :: N_lvls, Ng, gl_count, WinSize, MidWin, No_t, Ngt

    integer, intent(out) :: Totnp, Brkpt, Ier

    type(ellipse), intent(in out) :: Elvar

    real(r8), intent(out) :: h_path(:), z_path(:), t_path(:), phi_path(:), &
                             dhdz_path(:), phi_eta(:,:)

    real(r8), intent(in) :: htan, z_glgrid(:)
    real(r8), intent(in) :: h_glgrid(:,:), t_glgrid(:,:), t_phi_basis(:), &
                            dhdz_glgrid(:,:)

    !  ----------------
    !  Local variables:
    !  ----------------

    integer, dimension(1) :: brk_pt            ! Result of minloc
    integer :: i, j, k, l, m, n, n_d, npp, ibrk, Ngp1, iter

    real(r8) :: h, s, r, dz, rs, phi, rss, dhdz

    integer, pointer, dimension(:) :: cndx=>NULL()
    integer, pointer, dimension(:) :: vtp_ndx=>NULL()

    real(r8), pointer, dimension(:) :: dum_z=>NULL(), dum_h=>NULL(),&
                                       dum_phi=>NULL(), prev_h=>NULL()

    real(r8), pointer, dimension(:,:) :: h_a=>NULL()

    ! Allocate enough space for the various (temporary) arrys we are going
    ! to use...

    Ier = 0
    Ngp1 = Ng + 1

!?? print*,'Ngt is:',ngt
    call Allocate_test ( cndx, ngt, 'cndx', ModuleName )
    call Allocate_test ( vtp_ndx, ngt, 'vtp_ndx', ModuleName )
    call Allocate_test ( dum_z, ngt, 'dum_z', ModuleName )
    call Allocate_test ( dum_h, ngt, 'dum_h', ModuleName )
    call Allocate_test ( prev_h, ngt, 'prev_h', ModuleName )
    call Allocate_test ( dum_phi, ngt, 'dum_phi', ModuleName )
    call Allocate_test ( h_a, ngt, WinSize, 'h_a', ModuleName )
!?? print*,'WinSize is:',WinSize

!   Initialize all arrays:

    n_d = 0
    cndx = 0
    vtp_ndx = 0

    h_a = 0.0
    dum_z = 0.0
    dum_h = 0.0
    dum_phi = 0.0
    phi_eta = 0.0

    r = -999.99
    dz = r / 57.2958
    z_path = r
    h_path = r
    t_path = r
    phi_path = dz
    dhdz_path = 0.0

!   Define the various ELLIPSE variables needed for computations:

    elvar%ht = htan
    elvar%earthx = (htan < -0.01)

    if ( htan < -0.01 ) then
      elvar%rr = (elvar%roc + elvar%ht) / elvar%roc
      call earth_intersection ( elvar, rs )
    else
      elvar%rr = 0.0D0
      elvar%phi_s = elvar%phi_tan
      elvar%nphi_s = elvar%nphi_tan
    end if

!   Define the index points of the tangent locations:

    k = 2 * N_lvls

!   Store the initial guess heights in dum_h. This estimate is the h_glgrid
!   at the "Center Phi", Also compute the Path zeta on the GL grid:

    npp = 0
    l = gl_count + 1

    do
      do n = 1, Ngp1
        l = l - 1
        npp = npp + 1
        vtp_ndx(npp) = l
        dum_z(npp) = z_glgrid(l)
        dum_h(npp) = h_glgrid(l,MidWin)
      end do
      h = h_glgrid(l-1,MidWin)
      if ( abs(h-htan) < 0.001  .OR. l-Ngp1 < 1) exit
    end do

    l = l - 1
    npp = npp + 1
    vtp_ndx(npp) = l
    dum_z(npp) = z_glgrid(l)
    dum_h(npp) = h_glgrid(l,MidWin)

    l = l - 1
    ibrk = npp + 1

    do
      do n = 1, Ngp1
        l = l + 1
        npp = npp + 1
        vtp_ndx(npp) = l
        dum_z(npp) = z_glgrid(l)
        dum_h(npp) = h_glgrid(l,MidWin)
      end do
      if ( l+Ngp1 > gl_count) EXIT
    end do

    l = l + 1
    npp = npp + 1
    vtp_ndx(npp) = l
    dum_z(npp) = z_glgrid(l)
    dum_h(npp) = h_glgrid(l,MidWin)

!  Cast the h_glgrid onto the path

    DO m = 1, WinSize
      h_a(1:npp,m) = h_glgrid((/(vtp_ndx(k),k=1,npp)/),m)
    END DO

    iter = 0

    DO

      iter = iter + 1

      if(iter == 1) then

        elvar%ps = -1.0D0     ! Right hand side ray
        DO i = 1, ibrk-1
          CALL H_TO_S_PHI(elvar,dum_h(i),s,dum_phi(i))
        END DO
        elvar%ps = 1.0D0     ! Left hand side ray
        DO i = ibrk, npp
          CALL H_TO_S_PHI(elvar,dum_h(i),s,dum_phi(i))
        END DO

      else

        DO j = 1, n_d
          i = cndx(j)
          elvar%ps = -1.0D0                  ! Right hand side ray
          IF(i >= ibrk) elvar%ps = 1.0D0     ! Left hand side ray
          CALL H_TO_S_PHI(elvar,dum_h(i),s,dum_phi(i))
        END DO

      endif

!  Compute the Phi_Eta matrix based on the "Path" Phi's

      DO i = 1, npp
        r = dum_phi(i)
        DO m = 1, WinSize
          CALL get_one_eta(r,t_phi_basis,WinSize,m,phi_eta(i,m))
        END DO
      END DO
!
!  Re-Compute the estimate H:
!
      prev_h = dum_h
      dum_h = SUM(h_a*phi_eta,dim=2)
!
! Define the indecies array for which the process did not converged ..
!
!   cndx(1:) = 0
!   cndx = PACK((/(i,i=1,ngt)/),(abs(dum_h-prev_h) > 0.01))
!   n_d = COUNT(cndx > 0)
!
      n_d = 0
      do i = 1, npp
        if(abs(dum_h(i)-prev_h(i)) > 0.01) then
          n_d = n_d + 1
          cndx(n_d) = i
        endif
      end do

! **** ITERATE as needed

      IF(n_d == 0 .OR. iter == 10) EXIT

    END DO           ! On iter loop

!   For indices that did not converge within 10 iterations, use linear
!   interpolation on H and then recompute Phi and Phi_Eta for those indecies.

    IF(n_d > 1) THEN
      k = MAX(1,cndx(1)-1)
      l = MIN(npp,cndx(n_d)+1)
      dz = dum_z(l) - dum_z(k)
      dhdz = (dum_h(l) - dum_h(k)) / dz
      dO i = 1, n_d
        j = cndx(i)
        elvar%ps = -1.0D0
        IF(j >= ibrk) elvar%ps = 1.0D0
        dz = dum_z(j) - dum_z(k)
        dum_h(j) = dum_h(k) + dhdz * dz
        CALL H_TO_S_PHI(elvar,dum_h(j),s,dum_phi(j))
        DO m = 1, WinSize
          call get_one_eta ( dum_phi(j), t_phi_basis, WinSize, m, &
                          &  phi_eta(j,m) )
        END DO
      END DO
    END IF

!   **** CONVERGED !!
!   Now - Compute the path Temperature and path dh_dz

!   First, compute the path Temperature:
!   Cast the t_glgrid onto the path

    DO m = 1, WinSize
      h_a(1:npp,m) = t_glgrid((/(vtp_ndx(k),k=1,npp)/),m)
    END DO
    t_path = SUM(h_a*phi_eta,dim=2)

    totnp = npp
    z_path(1:npp)   = dum_z(1:npp)
    h_path(1:npp)   = dum_h(1:npp)
    phi_path(1:npp) = dum_phi(1:npp)

! Now compute the break points for the tangent layer for this tanget height

    brk_pt = minloc(abs(h_path(1:)-htan))
    brkpt = brk_pt(1)
!
! Second, compute the path dh_dz:
! Cast the dhdz_glgrid onto the path

    DO m = 1, WinSize
      h_a(1:npp,m) = dhdz_glgrid((/(vtp_ndx(k),k=1,npp)/),m)
    END DO
    dhdz_path = SUM(h_a*phi_eta,dim=2)

    call deallocate_test ( h_a, 'h_a', ModuleName )
    call deallocate_test ( dum_phi, 'dump_phi', ModuleName )
    call deallocate_test ( dum_h, 'dum_h', ModuleName )
    call deallocate_test ( prev_h, 'prev_h', ModuleName )
    call deallocate_test ( dum_z, 'dum_z', ModuleName )
    call deallocate_test ( cndx, 'cndx', ModuleName )
    call deallocate_test ( vtp_ndx, 'vtp_ndx', ModuleName )

    Return

  end subroutine Vert_To_Path
end module Vert_To_Path_M
! $Log$
! Revision 1.16  2001/06/07 23:39:32  pwagner
! Added Copyright statement
!
! Revision 1.15  2001/05/21 22:08:35  zvi
! A small modification in searching for tanget
!
! Revision 1.14  2001/04/23 21:43:28  zvi
! Introducing no_phi_t etc.
!
! Revision 1.13  2001/04/13 03:34:46  zvi
! Correcting minor error in allocation, cleaing up comments.
!
! Revision 1.12  2001/04/13 02:00:25  vsnyder
! Undo all of the lmin:lmax nonsense
!
! Revision 1.11  2001/04/13 01:44:36  vsnyder
! Work on moving window
!
! Revision 1.10  2001/04/13 01:13:59  vsnyder
! Use lmin:lmax for more dimensions
!
! Revision 1.9  2001/04/13 00:27:35  vsnyder
! Comment out some of Nathaniel's debugging print
!
! Revision 1.8  2001/04/12 21:41:24  livesey
! Changed allocatable to automatic then pointer.  Still fails on large chunks.
! Van will try it on LF95
!
! Revision 1.7  2001/03/31 23:40:56  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
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
