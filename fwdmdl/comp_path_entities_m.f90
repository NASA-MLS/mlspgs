module COMP_PATH_ENTITIES_M
  use MLSCommon, only: I4, R4, R8
  use L2PCDIM, only: N2lvl, MNP => max_no_phi
  use GL6P, only: NG
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component, &
                                  DEG2RAD
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
      HT2, RR, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_DERIVATIVE, &
                             PATH_VECTOR_2D
  use REFRACTION_M, only: REFRACTIVE_INDEX
  use VERT_TO_PATH_M, only: VERT_TO_PATH
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  use Molecules, only: l_h2o
  use Intrinsic, only: l_vmr
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  implicit NONE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------

SUBROUTINE comp_path_entities(fwdModelIn, fwdModelExtra, molecules, &
           n_lvls,no_t,gl_count,ndx_path,z_glgrid,t_glgrid,h_glgrid,&
           dhdz_glgrid,tan_hts,&
           no_tan_hts,n_sps,z_path,h_path,     &
           t_path,phi_path,n_path,dhdz_path,eta_phi,no_phi_t,       &
           t_phi_basis,spsfunc_path,no_mmaf,Ier)

!  ===============================================================
!  Declaration of variables for sub-program: comp_path_entities
!  ===============================================================

type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra
integer, dimension(:), intent(in) :: molecules

Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_t, n_sps, n_lvls, gl_count, &
             no_mmaf, no_phi_t
!
Integer(i4), INTENT(IN OUT) :: no_tan_hts

Integer(i4), INTENT(OUT) :: ier
!
Real(r8), INTENT(IN) :: z_glgrid(:), h_glgrid(:,:), t_glgrid(:,:)
Real(r8), INTENT(IN) :: dhdz_glgrid(:,:)

Real(r8), INTENT(IN) :: t_phi_basis(:)
Real(r8), INTENT(IN) :: tan_hts(:,:)

Type(path_index) , INTENT(OUT) :: ndx_path(:,:)
Type(path_vector), INTENT(OUT) :: z_path(:,:),t_path(:,:),h_path(:,:), &
           n_path(:,:),phi_path(:,:),dhdz_path(:,:), spsfunc_path(:,:,:)

Type(path_vector_2d), INTENT(OUT) :: eta_phi(:,:)
!
!  ----------------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, l, jp, sps_i, jj, kk, ih2o, lmin, lmax, &
               klo, khi

Real(r8) :: h, q, r, zeta, phi

Real(r8), DIMENSION(:)  , ALLOCATABLE :: zpath,tpath,hpath,ppath,dhdzp
Real(r8), DIMENSION(:,:), ALLOCATABLE :: phi_eta

type (VectorValue_T), pointer :: f, h2o

!  PFA variables:

  ier = 0

! Compute all the various integration paths according to tanget heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,STAT=i)
  ALLOCATE(zpath(ngt),hpath(ngt),tpath(ngt),ppath(ngt),dhdzp(ngt), &
 &         phi_eta(ngt,no_phi_t),STAT=ier)
  IF(ier /= 0) THEN
    Print *,'** Allocation Error in comp_path_entities: ?path ...'
    PRINT *,'** Error: ALLOCATION error in MAIN ..'
    PRINT *,'   STAT =',ier
    Return
  ENDIF

  jp = no_phi_t/2
  DO l = 1, no_mmaf
    lmin = max(1,l-jp)
    lmax = lmin + 2 * jp
    if(lmax > no_mmaf) then
      lmax = no_mmaf
      lmin = max(1,lmax - 2 * jp)
    endif
    print*,'Hello there:',jp,l,no_mmaf,lmin,lmax,shape(t_glgrid),lbound(t_glgrid,2)
    DO k = 1, no_tan_hts
      h = tan_hts(k,l)
      CALL vert_to_path(n_lvls,Ng,ngt,gl_count,no_phi_t,no_t,h,   &
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

! Create the specie function along the path for all species
!
  do j = 1, size(molecules)
    f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_vmr, molecule=molecules(j))
    jp = f%template%noInstances
    kk = f%template%noSurfs
    DO l = 1, no_mmaf
      DO k = 1, no_tan_hts
        jj = ndx_path(k,l)%total_number_of_elements
        DEALLOCATE(spsfunc_path(j,k,l)%values,STAT=i)
        ALLOCATE(spsfunc_path(j,k,l)%values(jj),STAT=ier)
        IF(ier /= 0) THEN
          PRINT *,'** Error: ALLOCATION error for spsfunc_path ..'
          PRINT *,'   STAT =',ier
          goto 99
        ENDIF
        do i = 1, jj
          zeta = z_path(k,l)%values(i)
          phi = phi_path(k,l)%values(i)
          if (f%template%logBasis) then
            Call TWO_D_POLATE(f%template%surfs(:,1), &
              & log(f%values), &
           &           kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, r)
            q = exp(r)
          else
            Call TWO_D_POLATE(f%template%surfs(:,1), &
              & f%values, kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, q)
          endif
          spsfunc_path(j,k,l)%values(i) = q
        end do
      end do
    END DO
  END DO

  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,phi_eta,STAT=i)
!
! Compute the relative refractive index minus one.
! Get the water mixing ratio function

  h2o => GetVectorQuantityByType( fwdModelIn, fwdModelExtra, &
    & quantityType=l_vmr, molecule=l_h2o, noError=.true.)
  print*,'Wet would have been:',associated(h2o)
  if (associated(h2o)) then
    jp = h2o%template%noInstances
    kk = h2o%template%noSurfs
  else
    jp = 0
    kk = 0
  endif

  DO l = 1, no_mmaf
    CALL refractive_index(h2o%values,h2o%template%surfs(:,1),&
      &             h2o%template%phi(1,:), &
  &                 kk,jp,ndx_path(1:,l),z_path(1:,l),t_path(1:,l),    &
  &                 phi_path(1:,l),n_path(1:,l),associated(h2o),no_tan_hts)
  END DO

 99  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,STAT=i)

  RETURN

END SUBROUTINE comp_path_entities

end module COMP_PATH_ENTITIES_M
! $Log$
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
! Revision 1.8  2001/03/21 18:39:43  livesey
! Fixed deg2rad bug
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
