module VERT_TO_PATH_M
  use Allocate_Deallocate, only: Allocate_test, Deallocate_Test
  use MLSCommon, only: I4, R4, R8
  use D_LINTRP_M, only: LINTRP
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

  subroutine Vert_To_Path ( Elvar, N_lvls, Ng, Ngt, Gl_count, No_phi_t, No_t, &
    &  Htan, Z_glgrid, T_glgrid, H_glgrid, Dhdz_glgrid, T_phi_basis, Z_path, &
    &  H_path, T_path, Phi_path, Dhdz_path, Phi_eta, Brkpt, Totnp, Ier )

    !  ===============================================================
    !  Declaration of variables for sub-program: vert_to_path
    !  ===============================================================
    !  ---------------------------
    !  calling sequence variables:
    !  ---------------------------
    integer, intent(in) :: n_lvls,Ng,gl_count,no_phi_t,no_t,ngt

    integer, intent(out) :: totnp,brkpt,Ier

    type(ellipse), intent(in out) :: elvar

    real(r8), intent(out) :: h_path(:), z_path(:), t_path(:), phi_path(:), &
                             dhdz_path(:), phi_eta(:,:)

    real(r8), intent(in) :: htan, z_glgrid(:)
    real(r8), intent(in) :: h_glgrid(:,:), t_glgrid(:,:), t_phi_basis(:), &
                            dhdz_glgrid(:,:)

    !  ----------------
    !  Local variables:
    !  ----------------

    integer :: i, j, k, l, m, n, jp, n_d, npp, ibrk, Ngp1, no_iter

    real(r8) :: h, s, r, dz, rs, phi, rss, dhdz, prev_h

    integer, pointer, dimension(:) :: cndx=>NULL()

    real(r8), pointer, dimension(:) :: dum_z=>NULL(), dum_h=>NULL(), dum_phi=>NULL()

    real(r8), pointer, dimension(:,:) :: h_a=>NULL()

    ! Allocate enough space for the various (temporary) arrys we are going
    ! to use...

    Ier = 0
    Ngp1 = Ng + 1

!?? print*,'Ngt is:',ngt
    call Allocate_test ( cndx, ngt, 'cndx', ModuleName )
    call Allocate_test ( dum_z, ngt, 'dum_z', ModuleName )
    call Allocate_test ( dum_h, ngt, 'dum_h', ModuleName )
    call Allocate_test ( dum_phi, ngt, 'dum_phi', ModuleName )
    call Allocate_test ( h_a,ngt,no_phi_t,'h_a',ModuleName )
    print*,'no_phi_t is:',no_phi_t

!     DEALLOCATE(cndx, dum_z, dum_h, dum_phi, STAT=i)
!      ALLOCATE(cndx(ngt), dum_z(ngt), dum_h(ngt), dum_phi(ngt), &
!     &         STAT = ier)
!      if ( ier /= 0 ) then
!        Ier = 1
!        Print *,'** Error: ALLOCATION error in VERT_TO_PATH routine ..'
!        GOTO 99
!      end if
!     DEALLOCATE(h_a, STAT=i)
!      ALLOCATE(h_a(ngt,no_phi_t),STAT=ier)
!      if ( ier /= 0 ) then
!        Ier = 1
!        Print *,'** Error: ALLOCATION error in VERT_TO_PATH routine ..'
!        GOTO 99
!      end if

!   Initialize all arrays:

    cndx = 0
    dum_z = 0.0
    dum_h = 0.0
    dum_phi = 0.0
    do j = 1, no_phi_t
      h_a(1:ngt,j) = 0.0
      phi_eta(1:ngt,j) = 0.0
    end do

    r = -999.99
    dz = r / 57.2958
    z_path(1:ngt) = r
    h_path(1:ngt) = r
    t_path(1:ngt) = r
    phi_path(1:ngt) = dz
    dhdz_path(1:ngt) = 0.0

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
    jp = (no_phi_t + 1) / 2

    do
      do n = 1, Ngp1
        l = l - 1
        npp = npp + 1
        dum_z(npp) = z_glgrid(l)
        dum_h(npp) = h_glgrid(l,jp)
      end do
      h = h_glgrid(l-1,jp)
      if ( h <= htan  .OR. l-Ngp1 < 1) exit
    end do

    l = l - 1
    npp = npp + 1
    dum_z(npp) = z_glgrid(l)
    dum_h(npp) = h_glgrid(l,jp)

    l = l - 1
    ibrk = npp + 1

    do
      do n = 1, Ngp1
        l = l + 1
        npp = npp + 1
        dum_z(npp) = z_glgrid(l)
        dum_h(npp) = h_glgrid(l,jp)
      end do
      if ( l+Ngp1 > gl_count) exit
    end do

    l = l + 1
    npp = npp + 1
    dum_z(npp) = z_glgrid(l)
    dum_h(npp) = h_glgrid(l,jp)

!    Cast the h_glgrid onto the path (using liner interpolation)

    do m = 1, no_phi_t
      call lintrp ( z_glgrid, dum_z, h_glgrid(1:,m), h_a(1:,m), gl_count, npp )
    end do

    n_d = ngt
    no_iter = 0

    do
      no_iter = no_iter + 1

!     Compute the Phi that goes with the Heights along the path:
!     Right hand side ray:

      l = -1
      elvar%ps = -1.0D0
      do i = 1, ibrk-1
        h = dum_h(i)
        phi = dum_phi(i)
        if ( n_d == ngt ) then
          call H_TO_S_PHI ( elvar, h, s, phi )
        else
          call hunt ( i, cndx, n_d, l, j )
          if ( cndx(l) == i .OR. cndx(j) == i) &
                  call H_TO_S_PHI ( elvar, h, s, phi )
        end if
        dum_phi(i) = phi
      end do

!     Left hand side ray:

      l = -1
      elvar%ps = 1.0D0
      do i = ibrk, npp
        h = dum_h(i)
        phi = dum_phi(i)
        if ( n_d == ngt ) then
          call H_TO_S_PHI ( elvar, h, s, phi )
        else
          call hunt ( i, cndx, n_d, l, j )
          if ( cndx(l) == i .OR. cndx(j) == i) &
                   call H_TO_S_PHI ( elvar, h, s, phi )
        end if
        dum_phi(i) = phi
      end do

!     Compute the Phi_Eta matrix based on the "Path" Phi's

      do i = 1, npp
        r = dum_phi(i)
        do m = 1, no_phi_t
          call get_one_eta ( r, t_phi_basis, no_phi_t, m, phi_eta(i,m) )
        end do
      end do

!     Re-Compute the estimate H:

      n_d = 0
      rss = 0.0D0
      do i = 1, npp
        prev_h = dum_h(i)
        dum_h(i) = SUM(h_a(i,:)*phi_eta(i,:))
        r = ABS(dum_h(i) - prev_h)
        if ( r > 0.01 ) then
          n_d = n_d + 1
          cndx(n_d) = i
          rss = rss + r * r
        end if
      end do

!   **** ITERATE

      if ( n_d <= 0 .OR. no_iter >= 10 ) exit
    end do

!   For indices that did not converge within 10 iterations, use linear
!   interpolation on H and then recompute Phi and Phi_Eta for those indecies.

    if ( n_d > 1 ) then
      k = MAX(1,cndx(1)-1)
      l = MIN(npp,cndx(n_d)+1)
      dz = dum_z(l) - dum_z(k)
      dhdz = (dum_h(l) - dum_h(k)) / dz
      do i = 1, n_d
        j = cndx(i)
        elvar%ps = -1.0D0
        dz = dum_z(j) - dum_z(k)
        dum_h(j) = dum_h(k) + dhdz * dz
        h = dum_h(j)
        if ( j >= ibrk) elvar%ps = 1.0D0
        call H_TO_S_PHI ( elvar, h, s, phi )
        dum_phi(j) = phi
        do m = 1, no_phi_t
          call get_one_eta ( dum_phi(j), t_phi_basis, no_phi_t, m, phi_eta(j,m) )
        end do
      end do
    end if

!   **** CONVERGED !!
!   Now - Compute the path Temperature, dh_dz  and dH_dTlm

!   First, compute the path Temperature:
!   Cast the t_glgrid onto the path (using liner interpolation)

    do m = 1, no_phi_t
      call lintrp ( z_glgrid, dum_z, t_glgrid(1:,m), h_a(1:,m), gl_count, npp )
    end do

    s = 1.0e10
    brkpt = -1
    totnp = npp

!    Also, compute the break points for the tangent layer for this tanget height

    do i = 1, npp
      t_path(i) = SUM(h_a(i,:)*phi_eta(i,:))
      z_path(i) = dum_z(i)
      h_path(i) = dum_h(i)
      phi_path(i) = dum_phi(i)
      r = abs(h_path(i)-htan)
      if ( r < s ) then
        s = r
        brkpt = i
      end if
    end do

!   Second, compute the path dh_dz:
!   Cast the dhdz_glgrid onto the path (using liner interpolation)

    do m = 1, no_phi_t
      call lintrp ( z_glgrid, dum_z, dhdz_glgrid(1:,m), h_a(1:,m), gl_count, npp )
    end do

    do i = 1, npp
      dhdz_path(i) = SUM(h_a(i,:)*phi_eta(i,:))
    end do
 
!?? print*,'Hello Im deallocating!'
    call deallocate_test ( h_a,'h_a',ModuleName )
    call deallocate_test ( dum_phi, 'dump_phi', ModuleName )
    call deallocate_test ( dum_h, 'dum_h', ModuleName )
    call deallocate_test ( dum_z, 'dum_z', ModuleName )
    call deallocate_test ( cndx, 'cndx', ModuleName )

!    99 DEALLOCATE(h_a, STAT=i)
!       DEALLOCATE(cndx, dum_z, dum_h, dum_phi, STAT=i)

    return
  end subroutine Vert_To_Path
end module Vert_To_Path_M
! $Log$
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
