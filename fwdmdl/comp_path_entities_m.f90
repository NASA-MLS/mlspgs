module Comp_Path_Entities_M
  use MLSCommon, only: R8
  use GL6P, only: NG
  use ELLIPSE_M, only: ELLIPSE
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_VECTOR_2D
  use VERT_TO_PATH_M, only: VERT_TO_PATH
  use VectorsModule, only: VectorValue_T

  implicit NONE
  private
  public :: Comp_Path_Entities

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------

  subroutine Comp_Path_Entities ( radiance, temperature, closestInstances, &
    &        N_lvls, No_t, Gl_count, Ndx_path, Z_glgrid, &
    &        T_glgrid, H_glgrid, Dhdz_glgrid, Tan_hts, No_tan_hts, Z_path, &
    &        H_path, T_path, Phi_path, Dhdz_path, Eta_phi, No_phi_t, &
    &        T_phi_basis, NoMAFs, PhiWindow, Elvar, Ier )

  !  ===============================================================
  !  Declaration of variables for sub-program: comp_path_entities
  !  ===============================================================

  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------
  type (VectorValue_T), intent(in) :: radiance
  type (VectorValue_T), intent(in) :: temperature
  integer, dimension(:), intent(in) :: closestInstances

  integer, intent(in) :: no_t, no_phi_t, N_lvls, gl_count, PhiWindow, NoMAFs
  !
  integer, intent(in out) :: No_tan_hts

  type(ellipse), intent(in out) :: Elvar(:)

  integer, intent(out) :: Ier
  !
  real(r8), intent(in) :: Z_glgrid(:), H_glgrid(:,:), T_glgrid(:,:)
  real(r8), intent(in) :: Dhdz_glgrid(:,:)

  real(r8), intent(in) :: T_phi_basis(:)
  real(r8), intent(in) :: Tan_hts(:,:)

  type(path_index) , intent(inout) :: Ndx_path(:,:)
  type(path_vector), intent(inout) :: Z_path(:,:), T_path(:,:), H_path(:,:), &
                                      Phi_path(:,:), Dhdz_path(:,:)

  type(path_vector_2d), intent(inout) :: Eta_phi(:,:)
  !
  !  ----------------------
  !  Local variables:
  !  ----------------

  Integer :: i, k, l, jj, kk, lmin, lmax, klo, khi, ngt, maf, &
             WinSize, MidWin

  Real(r8) :: h, q, r, zeta, phi

  Real(r8), dimension(:)  , allocatable :: zpath, tpath, hpath, ppath, dhdzp
  Real(r8), dimension(:,:), allocatable :: phi_eta

  ier = 0
  ngt = 2 * (Ng+1) * (N_lvls+1)

! Compute all the various integration paths according to tangent heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  allocate ( zpath(ngt), hpath(ngt), tpath(ngt), ppath(ngt), dhdzp(ngt), &
             & stat=ier )
  if(ier /= 0) then
    print *,'** Allocation Error in comp_path_entities: ?path ...'
    print *,'** Error: ALLOCATION error in MAIN ..'
    print *,'   STAT =',ier
    Return
  end if


  do maf = 1, NoMAFs
    l = closestInstances(maf)
    lmin = max(1,l-phiWindow/2)
    lmax = min(no_phi_t,l+phiWindow/2)
    WinSize = lmax-lmin+1
    MidWin = l - lmin + 1
    print*,'Hello there',l,lmin,lmax,winSize,midwin
    allocate ( phi_eta(ngt,lmin:lmax), stat=ier )
    if ( ier /= 0 ) then
      print *,'** Error Allocating phi_eta in routine: comp_path_entities..'
      print *,'   STAT =',ier
      return
    endif
    print*,'Elvar:',elvar(l)%phi_tan,elvar(maf)%rr
    do k = 1, no_tan_hts
      h = tan_hts(k,l)
      call vert_to_path ( elvar(l), n_lvls, Ng, ngt, gl_count, WinSize, &
        & MidWin, no_t, h, z_glgrid, t_glgrid(1:,lmin:lmax), &
        & h_glgrid(1:,lmin:lmax), dhdz_glgrid(1:,lmin:lmax), &
        & t_phi_basis(lmin:lmax), zpath, hpath, &
        & tpath, ppath, dhdzp, phi_eta, klo, khi, ier )
      if(ier /= 0) return
      allocate ( z_path(k,maf)%values(khi), h_path(k,maf)%values(khi),   &
        &        t_path(k,maf)%values(khi), phi_path(k,maf)%values(khi), &
        &        dhdz_path(k,maf)%values(khi), stat=ier )
      if ( ier == 0 ) allocate ( eta_phi(k,maf)%values(khi,lmin:lmax),&
                                & stat=ier )
      if ( ier /= 0 ) then
        print *,'** Error: ALLOCATION error in routine: comp_path_entities ..'
        print *,'   STAT =',ier
        return
      endif
      ndx_path(k,maf)%break_point_index = klo
      ndx_path(k,maf)%total_number_of_elements = khi
      z_path(k,maf)%values(1:khi) = zpath(1:khi)
      t_path(k,maf)%values(1:khi) = tpath(1:khi)
      h_path(k,maf)%values(1:khi) = hpath(1:khi)
      phi_path(k,maf)%values(1:khi) = ppath(1:khi)
      dhdz_path(k,maf)%values(1:khi) = dhdzp(1:khi)
      eta_phi(k,maf)%values(1:khi,lmin:lmax) = phi_eta(1:khi,lmin:lmax)
    end do ! k = 1, no_tan_hts
    deallocate ( phi_eta, stat=i )
  end do ! maf = 1, NoMAFs

 99  deallocate ( zpath, hpath, tpath, ppath, dhdzp, stat=i )

  return

end subroutine Comp_Path_Entities

end module Comp_Path_Entities_M
! $Log$
! Revision 1.31  2001/04/25 00:38:48  zvi
! Fix bug in elvar..
!
! Revision 1.30  2001/04/23 22:12:08  livesey
! Whoops, bug fix.
!
! Revision 1.29  2001/04/23 21:56:25  livesey
! Accepts closestInstances as parameter
!
! Revision 1.28  2001/04/23 21:43:28  zvi
! Introducing no_phi_t etc.
!
! Revision 1.27  2001/04/19 22:09:18  livesey
! Modified to deal with no mafs /= T%noInstances
!
! Revision 1.26  2001/04/19 06:48:13  zvi
! Fixing memory leaks..
!
! Revision 1.25  2001/04/13 03:34:45  zvi
! Correcting minor error in allocation, cleaing up comments.
!
! Revision 1.24  2001/04/13 02:00:55  vsnyder
! Finish changing phi_eta back to allocatable
!
! Revision 1.23  2001/04/13 02:00:10  vsnyder
! Change phi_eta back to allocatable
!
! Revision 1.22  2001/04/13 01:44:36  vsnyder
! Work on moving window
!
! Revision 1.21  2001/04/13 00:27:20  vsnyder
! Limit amount of phi_eta assigned
!
! Revision 1.20  2001/04/12 21:43:50  livesey
! Some attempts to fix array bounds in vert_to_path call
!
! Revision 1.19  2001/04/12 01:02:11  zvi
! Adding phiwindow to calling seq.
!
! Revision 1.18  2001/04/07 23:51:17  zvi
! New code - move the spsfunc & refraction along the path to get_path_spsfunc
!
! Revision 1.17  2001/04/07 01:21:10  livesey
! Removed print statement
!
! Revision 1.16  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.15  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.14  2001/03/30 01:40:24  livesey
! Removed some arguments
!
! Revision 1.13  2001/03/30 00:07:57  livesey
! Removed more arguments
!
! Revision 1.12  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.11  2001/03/26 21:06:31  zvi
! *** empty log message ***
!
! Revision 1.10  2001/03/26 17:56:14  zvi
! New codes to deal with dh_dt_path issue.. now being computed on the fly
!
! Revision 1.9  2001/03/21 18:40:09  livesey
! Removed dump statements
!
! Revision 1.7  2001/03/21 06:30:10  livesey
! Minor change, still wrong
!
! Revision 1.6  2001/03/21 02:12:54  livesey
! Interim version, got mr_f in but not working yet.
!
! Revision 1.5  2001/03/21 01:09:42  livesey
! Commented out dh_dt_path to save memory
!
! Revision 1.4  2001/03/20 11:03:15  zvi
! Fixing code for "real" data run, increase dim. etc.
!
! Revision 1.3  2001/03/05 21:37:20  zvi
! New filter format
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
