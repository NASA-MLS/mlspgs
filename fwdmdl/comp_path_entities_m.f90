module COMP_PATH_ENTITIES_M
  use MLSCommon, only: I4, R8
  use GL6P, only: NG
  use ELLIPSE_M, only: ELLIPSE
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_VECTOR_2D
  use VERT_TO_PATH_M, only: VERT_TO_PATH
  implicit NONE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------
contains
!---------------------------------------------------------------------

SUBROUTINE comp_path_entities(n_lvls,no_t,gl_count,ndx_path,z_glgrid,    &
           t_glgrid, h_glgrid, dhdz_glgrid, tan_hts, no_tan_hts, z_path, &
           h_path,t_path,phi_path,dhdz_path,eta_phi,no_phi_t,t_phi_basis,&
           no_mmaf,phiWindow,elvar,Ier)

!  ===============================================================
!  Declaration of variables for sub-program: comp_path_entities
!  ===============================================================

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_t,n_lvls,gl_count,phiWindow,no_mmaf,no_phi_t
!
Integer(i4), INTENT(IN OUT) :: no_tan_hts

Type(ELLIPSE), INTENT(IN OUT) :: elvar(:)

Integer(i4), INTENT(OUT) :: ier
!
Real(r8), INTENT(IN) :: z_glgrid(:), h_glgrid(:,:), t_glgrid(:,:)
Real(r8), INTENT(IN) :: dhdz_glgrid(:,:)

Real(r8), INTENT(IN) :: t_phi_basis(:)
Real(r8), INTENT(IN) :: tan_hts(:,:)

Type(path_index) , INTENT(OUT) :: ndx_path(:,:)
Type(path_vector), INTENT(OUT) :: z_path(:,:),t_path(:,:),h_path(:,:), &
                                  phi_path(:,:),dhdz_path(:,:)

Type(path_vector_2d), INTENT(OUT) :: eta_phi(:,:)
!
!  ----------------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, l, jj, kk, lmin, lmax, klo, khi, ngt

Real(r8) :: h, q, r, zeta, phi

Real(r8), DIMENSION(:)  , ALLOCATABLE :: zpath,tpath,hpath,ppath,dhdzp
Real(r8), DIMENSION(:,:), ALLOCATABLE :: phi_eta

  ier = 0
  ngt = 2 * (Ng+1) * (N_lvls+1)

! Compute all the various integration paths according to tanget heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,phi_eta,STAT=i)
  ALLOCATE(zpath(ngt),hpath(ngt),tpath(ngt),ppath(ngt),dhdzp(ngt), &
 &         phi_eta(ngt,no_phi_t),STAT=ier)
  IF(ier /= 0) THEN
    Print *,'** Allocation Error in comp_path_entities: ?path ...'
    PRINT *,'** Error: ALLOCATION error in MAIN ..'
    PRINT *,'   STAT =',ier
    Return
  ENDIF

  DO l = 1, no_mmaf
    lmin = max(1,l-phiWindow)
    lmax = min(l+phiWindow/2,no_mmaf)
    DO k = 1, no_tan_hts
      h = tan_hts(k,l)
      CALL vert_to_path(elvar(l),n_lvls,Ng,ngt,gl_count,no_phi_t,no_t,h,&
           z_glgrid,t_glgrid(1:,lmin:lmax),h_glgrid(1:,lmin:lmax),  &
           dhdz_glgrid(1:,lmin:lmax),t_phi_basis,zpath,hpath,tpath, &
           ppath,dhdzp,phi_eta,klo,khi,Ier)
      IF(ier /= 0) RETURN
      DEALLOCATE(z_path(k,l)%values, h_path(k,l)%values,   &
                 t_path(k,l)%values, phi_path(k,l)%values, &
                 dhdz_path(k,l)%values, eta_phi(k,l)%values,STAT=i)
      ALLOCATE(z_path(k,l)%values(khi), h_path(k,l)%values(khi),   &
               t_path(k,l)%values(khi), phi_path(k,l)%values(khi), &
               dhdz_path(k,l)%values(khi), STAT=j)
      IF(j == 0) &
     &    ALLOCATE(eta_phi(k,l)%values(khi,no_phi_t),STAT=j)
      IF(j /= 0) THEN
        ier = j
        PRINT *,'** Error: ALLOCATION error in routine: comp_path_entities ..'
        PRINT *,'   STAT =',ier
        RETURN
      ENDIF
      ndx_path(k,l)%break_point_index = klo
      ndx_path(k,l)%total_number_of_elements = khi
      z_path(k,l)%values(1:khi) = zpath(1:khi)
      t_path(k,l)%values(1:khi) = tpath(1:khi)
      h_path(k,l)%values(1:khi) = hpath(1:khi)
      phi_path(k,l)%values(1:khi) = ppath(1:khi)
      dhdz_path(k,l)%values(1:khi) = dhdzp(1:khi)
      eta_phi(k,l)%values(1:khi,1:no_phi_t) = phi_eta(1:khi,1:no_phi_t)
    END DO
  END DO

 99  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,phi_eta,STAT=i)

  Return

END SUBROUTINE comp_path_entities

end module COMP_PATH_ENTITIES_M
! $Log$
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
