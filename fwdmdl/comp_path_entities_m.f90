module COMP_PATH_ENTITIES_M
  use MLSCommon, only: I4, R4, R8
  use L2PCDIM, only: N2lvl, MNP => max_no_phi
  use GL6P, only: NG
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP, GEOM_PARAM, PFA_SLAB
  use L2PC_FILE_PARAMETERS, only: MXCO => max_no_elmnts_per_sv_component
  use ELLIPSE, only: A2, C2, C2OA2, CPT, SPT, CPS, SPS, CPTS, SPTS, HT, &
      HT2, RR, PHI_TAN, NPHI_TAN, PHI_S, NPHI_S, PS, ROC, XOC, YOC, EARTHX
  use PATH_ENTITIES_M, only: PATH_BETA, PATH_INDEX, PATH_VECTOR, &
                             PATH_DERIVATIVE
  use REFRACTION_M, only: REFRACTIVE_INDEX
  use PFA_PREP_M, only: PFA_PREP
  use VERT_TO_PATH_M, only: VERT_TO_PATH
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------

SUBROUTINE comp_path_entities(primag,n_lvls,no_t,gl_count,ndx_path,z_glgrid,&
           t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,atmospheric,no_atmos, &
           f_basis,mr_f,no_coeffs_f,freq,tan_hts,no_tan_hts,n_sps,no_pfa_ch,&
           no_filt_pts,pfa_ch,pfa_spectrum,f_grid_fltr,fltr_func,InDir,     &
           ld,band,no_phi_f,f_phi_basis,z_path,h_path,t_path,phi_path,      &
           n_path,dhdz_path,dh_dt_path,no_phi_t,t_phi_basis,spsfunc_path,   &
           no_ptg_frq,ptg_frq_grid,is_f_log,beta_path,no_mmaf,phi_tan_mmaf, &
           ier)

!  ===============================================================
!  Declaration of variables for sub-program: comp_path_entities
!  ===============================================================

Integer(i4), PARAMETER :: ngt = (Ng+1) * N2lvl

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_atmos, no_t, n_sps, n_lvls, gl_count, &
             no_mmaf, no_phi_t, no_coeffs_f(*), no_phi_f(:), band, &
             no_pfa_ch, no_filt_pts, pfa_ch(*)
!
Integer(i4), INTENT(IN OUT) :: no_tan_hts, ld

Integer(i4), INTENT(OUT) :: no_ptg_frq(*),ier
!
Real(r8), INTENT(IN) :: phi_tan_mmaf(*)

Real(r8), INTENT(IN) :: z_glgrid(:), h_glgrid(:,:), t_glgrid(:,:)
Real(r8), INTENT(IN) :: dh_dt_glgrid(:,:,:), dhdz_glgrid(:,:)

Real(r8), INTENT(IN) :: mr_f(:,:,:), f_basis(:,:)

Real(r8), INTENT(IN) :: t_phi_basis(:)
Real(r8), INTENT(IN) :: f_phi_basis(:,:), tan_hts(:,:)

Type(path_beta), INTENT(OUT) :: beta_path(:,:,:,:)

Type(path_index) , INTENT(OUT) :: ndx_path(:,:)
Type(path_vector), INTENT(OUT) :: z_path(:,:),t_path(:,:),h_path(:,:), &
                   n_path(:,:),phi_path(:,:),dhdz_path(:,:),ptg_frq_grid(:),&
                   spsfunc_path(:,:,:)

Type(path_derivative), INTENT(OUT) :: dh_dt_path(:,:)

Real(r8), INTENT(OUT) :: freq(*)

Real(r8), INTENT(OUT) :: fltr_func(:,:),f_grid_fltr(:,:)

Logical, INTENT(IN) :: is_f_log(*)
!
Character (LEN=*), INTENT(IN) :: InDir
Character (LEN=*), INTENT(IN) :: primag
!  ----------------------
!  Local variables:
!  ----------------

Logical :: wet
Integer(i4) :: i, j, k, l, jp, sps_i, jj, kk, ih2o, lmin, lmax, &
               klo, khi

Real(r4) :: dhdtp(ngt,mnp,mxco)

Real(r8) :: h, q, r, zeta, phi

Real(r8), DIMENSION(:)  , ALLOCATABLE :: t_phi_tan
Real(r8), DIMENSION(:)  , ALLOCATABLE :: zpath,tpath,hpath,ppath,dhdzp

Real(r8), DIMENSION(:,:), ALLOCATABLE :: f_phi_tan

!  PFA variables:

type (atmos_comp), intent(inout) :: ATMOSPHERIC(*)
type (pfa_slab)  , intent(inout) :: PFA_SPECTRUM(6,*)
!
  ier = 0

! Compute all the various integration paths according to tanget heights.
! Get the z, t, h, phi, dhdz & dh_dt arrays on these paths.

  DEALLOCATE(t_phi_tan,STAT=i)
  ALLOCATE(t_phi_tan(no_phi_t),STAT=ier)
  if(ier /= 0) then
    Print *,'** Allocation Error in comp_path_entities: t_phi_tan..'
    Return
  endif

  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,STAT=i)
  ALLOCATE(zpath(ngt),hpath(ngt),tpath(ngt),ppath(ngt),dhdzp(ngt), &
  &         STAT=ier)
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
      lmin = lmax - 2 * jp
    endif
    phi_tan = phi_tan_mmaf(l)
    t_phi_tan(1:no_phi_t) = t_phi_basis(1:no_phi_t) + phi_tan
    DO k = 1, no_tan_hts
      h = tan_hts(k,l)
      CALL vert_to_path(n_lvls,Ng,ngt,gl_count,no_phi_t,no_t,h,   &
           z_glgrid,t_glgrid(1:,lmin:lmax),h_glgrid(1:,lmin:lmax),  &
           dhdz_glgrid(1:,lmin:lmax),dh_dt_glgrid(1:,lmin:lmax,1:), &
           t_phi_tan,zpath,hpath,tpath,ppath,dhdzp,dhdtp,klo,khi,Ier)
      IF(ier /= 0) RETURN
      DEALLOCATE(z_path(k,l)%values, h_path(k,l)%values,   &
                 t_path(k,l)%values, phi_path(k,l)%values, &
                 dhdz_path(k,l)%values, dh_dt_path(k,l)%values,STAT=i)
      ALLOCATE(z_path(k,l)%values(khi), h_path(k,l)%values(khi),   &
               t_path(k,l)%values(khi), phi_path(k,l)%values(khi), &
               dhdz_path(k,l)%values(khi), STAT=j)
      IF(j == 0) &
     &     ALLOCATE(dh_dt_path(k,l)%values(khi,no_phi_t,no_t),STAT=j)
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
      dh_dt_path(k,l)%values(1:khi,1:no_phi_t,1:no_t) = &
     &                              dhdtp(1:khi,1:no_phi_t,1:no_t)
    END DO
  END DO

  DEALLOCATE(t_phi_tan,STAT=i)
  DEALLOCATE(f_phi_tan,STAT=i)
!
  l = MAXVAL(no_phi_f,n_sps)
  ALLOCATE(f_phi_tan(l,n_sps),STAT=ier)
  IF(ier /= 0) then
    Print *,'** ALLOCATE Error, comp_path_entities: f_phi_tan, STAT =',ier
    goto 99
  endif
!
! Create the specie function along the path for all species
!
  DO l = 1, no_mmaf
    phi_tan = phi_tan_mmaf(l)
    DO k = 1, no_tan_hts
      jj = ndx_path(k,l)%total_number_of_elements
      do j = 1, n_sps
        jp = no_phi_f(j)
        kk = no_coeffs_f(j)
        f_phi_tan(1:jp,j) = f_phi_basis(1:jp,j) + phi_tan
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
          if (is_f_log(j)) then
            Call TWO_D_POLATE(f_basis(1:,j), LOG(mr_f(1:kk,1:jp,j)), &
           &           kk, f_phi_tan(1:,j), jp, zeta, phi, r)
            q = exp(r)
          else
            Call TWO_D_POLATE(f_basis(1:,j), mr_f(1:kk,1:jp,j), kk,  &
           &                  f_phi_tan(1:,j), jp, zeta, phi, q)
          endif
          spsfunc_path(j,k,l)%values(i) = q
        end do
      end do
    END DO
  END DO

  DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,STAT=i)
!
! Compute the relative refractive index minus one.
! Get the water mixing ratio function

  sps_i = 1
  DO WHILE (atmospheric(sps_i)%name /= 'H2O' .AND. sps_i <= no_atmos)
    sps_i = sps_i + 1
  END DO

  IF (atmospheric(sps_i)%name == 'H2O') THEN
    j = sps_i
    ih2o = sps_i
  ELSE             ! Ignore water contribution to refractive index (dry air)
    j = 1
    ih2o = 0
  END IF

  wet = (ih2o > 0)
  jp = no_phi_f(j)
  kk = no_coeffs_f(j)

  DO l = 1, no_mmaf
    phi_tan = phi_tan_mmaf(l)
    f_phi_tan(1:jp,j) = f_phi_basis(1:jp,j) + phi_tan
    CALL refractive_index(mr_f(1:,1:,j),f_basis(1:,j),f_phi_tan(1:,j), &
  &                 kk,jp,ndx_path(1:,l),z_path(1:,l),t_path(1:,l),    &
  &                 phi_path(1:,l),n_path(1:,l),wet,no_tan_hts)
  END DO

! Create filter grids & functions for PFA calculations

  IF(no_pfa_ch > 0) THEN
    CALL pfa_prep(atmospheric,band,no_atmos,no_pfa_ch,no_filt_pts,pfa_ch, &
         pfa_spectrum,f_grid_fltr,freq,fltr_func,no_tan_hts,ndx_path,     &
         no_ptg_frq,ptg_frq_grid,z_path,t_path,beta_path,InDir,ld,primag, &
         no_mmaf,ier)
    IF(ier /= 0) goto 99
  ENDIF

 99  DEALLOCATE(f_phi_tan,t_phi_tan,STAT=i)
     DEALLOCATE(zpath,hpath,tpath,ppath,dhdzp,STAT=i)

  RETURN

END SUBROUTINE comp_path_entities

end module COMP_PATH_ENTITIES_M
! $Log$
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
