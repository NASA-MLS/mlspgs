module GET_PATH_SPSFUNC_NGRID_M
  use MLSCommon, only: I4, R8
  use L2PC_FILE_PARAMETERS, only: DEG2RAD
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR
  use REFRACTION_M, only: REFRACTIVE_INDEX
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

SUBROUTINE get_path_spsfunc_ngrid(fwdModelIn, fwdModelExtra, molecules, &
           ndx_path, no_tan_hts, z_path, t_path, phi_path, n_path,     &
           spsfunc_path, Ier)

!  =================================================================
!  Declaration of variables for sub-program: get_path_spsfunc_ngrid
!  =================================================================

type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra
integer, dimension(:), intent(in) :: molecules

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_tan_hts
!
Type(path_index) , INTENT(IN) :: ndx_path(:)
Type(path_vector), INTENT(IN) :: z_path(:),t_path(:),phi_path(:)

Integer(i4), INTENT(OUT) :: ier
Type(path_vector), INTENT(INOUT) :: n_path(:), spsfunc_path(:,:)
!
!  ----------------------
!  Local variables:
!  ----------------

Integer(i4) :: i, j, k, jp, jj, kk

Real(r8) :: q, r, zeta, phi

type (VectorValue_T), pointer :: f, h2o

  ier = 0
!
! Create the specie function along the path for all species
!
  do j = 1, size(molecules)
    f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_vmr, molecule=molecules(j))
    jp = f%template%noInstances
    kk = f%template%noSurfs
    DO k = 1, no_tan_hts
      jj = ndx_path(k)%total_number_of_elements
      DEALLOCATE(spsfunc_path(j,k)%values,STAT=i)
      ALLOCATE(spsfunc_path(j,k)%values(jj),STAT=ier)
      IF(ier /= 0) THEN
        PRINT *,'** Error: ALLOCATION error for spsfunc_path ..'
        PRINT *,'   STAT =',ier
        Return
      ENDIF
      do i = 1, jj
        zeta = z_path(k)%values(i)
        phi = phi_path(k)%values(i)
        if (f%template%logBasis) then
          Call TWO_D_POLATE(f%template%surfs(:,1), Log(f%values), &
         &         kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, r)
          q = exp(r)
        else
          Call TWO_D_POLATE(f%template%surfs(:,1), f%values, &
            &      kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, q)
        endif
        spsfunc_path(j,k)%values(i) = q
      end do
    end do
  END DO
!
! Compute the relative refractive index minus one.
! Get the water mixing ratio function

  h2o => GetVectorQuantityByType( fwdModelIn, fwdModelExtra, &
    & quantityType=l_vmr, molecule=l_h2o, noError=.true.)
  if (associated(h2o)) then
    jp = h2o%template%noInstances
    kk = h2o%template%noSurfs
  else
    jp = 0
    kk = 0
  endif

  CALL refractive_index(h2o%values,h2o%template%surfs(:,1),&
  &    h2o%template%phi(1,:),kk,jp,ndx_path,z_path,t_path, &
  &    phi_path,n_path,associated(h2o),no_tan_hts)

  Return

 END SUBROUTINE get_path_spsfunc_ngrid

end module GET_PATH_SPSFUNC_NGRID_M
! $Log$
! Revision 1.3  2001/04/19 06:48:14  zvi
! Fixing memory leaks..
!
! Revision 1.2  2001/04/07 23:51:17  zvi
! New code - move the spsfunc & refraction along the path to get_path_spsfunc
!
! Revision 1.1  2001/04/07 23:45:49  zvi
! New routine to calculate spsfunc & refraction,(not done in comp_path )
!
! Initial version, 2001/04/07 14:00:00  Z. Shippony
