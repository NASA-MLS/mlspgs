! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_PATH_SPSFUNC_NGRID_M
  use MLSCommon, only: I4, R8
  use GLNP, only: NG
  use L2PC_FILE_PARAMETERS, only: DEG2RAD
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_INT_VECTOR_2D
  use REFRACTION_M, only: REFRACTIVE_INDEX
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use Molecules, only: l_h2o
  use Intrinsic, only: l_vmr
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  implicit NONE

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------

  subroutine Get_path_spsfunc_ngrid ( fwdModelIn, fwdModelExtra, molecules, &
    &        ndx_path, no_tan_hts, z_path, t_path, phi_path, n_path,     &
    &        spsfunc_path, sps_zeta_loop, sps_phi_loop, Ier )

  !  =================================================================
  !  Declaration of variables for sub-program: get_path_spsfunc_ngrid
  !  =================================================================

  type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra
  integer, dimension(:), intent(in) :: molecules

  !  ---------------------------
  !  Calling sequence variables:
  !  ---------------------------
  Integer(i4), INTENT(IN) :: no_tan_hts

  Type(path_index) , intent(in) :: ndx_path(:)
  Type(path_vector), intent(in) :: z_path(:),t_path(:),phi_path(:)

  Integer(i4), intent(out) :: ier
  Type(path_vector), intent(inout) :: n_path(:), spsfunc_path(:,:)

  type(path_int_vector_2d), intent(inout), dimension(:,:) :: sps_phi_loop
  type(path_int_vector_2d), intent(inout), dimension(:,:) :: sps_zeta_loop

  !  ----------------------
  !  Local variables:
  !  ----------------

  Integer(i4) :: i, j, k, jp, jj, kk, brkpt, no_ele, Ngp1

  Real(r8) :: q, r, zeta, phi

  type (VectorValue_T), pointer :: f, h2o

  ier = 0

  ! Create the specie function along the path for all species

  do j = 1, size(molecules)
    f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_vmr, molecule=molecules(j))
    jp = f%template%noInstances
    kk = f%template%noSurfs
    do k = 1, no_tan_hts
      no_ele = ndx_path(k)%total_number_of_elements
      deallocate ( spsfunc_path(j,k)%values, STAT=i )
      allocate ( spsfunc_path(j,k)%values(no_ele), STAT=ier )
      if ( ier /= 0 ) then
        print *,'** Error: ALLOCATION error for spsfunc_path ..'
        print *,'   STAT =',ier
        Return
      end if
      do i = 1, no_ele
        zeta = z_path(k)%values(i)
        phi = phi_path(k)%values(i)
        if (f%template%logBasis) then
          Call TWO_D_POLATE(f%template%surfs(:,1), Log(f%values), &
         &         kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, r)
          q = exp(r)
        else
          Call TWO_D_POLATE(f%template%surfs(:,1), f%values, &
            &      kk, Deg2Rad*f%template%phi(1,:), jp, zeta, phi, q)
        end if
        spsfunc_path(j,k)%values(i) = q
      end do
    end do
  end do
!
! Establish the Specie coefficient limits arrays - which will prevent
! using zeta/phi which produce only zeros (This is the Species coeff. loop)
!
  Ngp1 = Ng+1
  do j = 1, size(molecules)
    f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=molecules(j))
    jp = f%template%noInstances
    kk = f%template%noSurfs
    do k = 1, no_tan_hts
      brkpt = ndx_path(k)%break_point_index
      no_ele = ndx_path(k)%total_number_of_elements
      jj = (no_ele + Ng) / Ngp1
      allocate(sps_phi_loop(j,k)%values(jj,2), &
            &  sps_zeta_loop(j,k)%values(jj,2),stat=ier)
      IF(ier /= 0) THEN
        PRINT *,'** Error: ALLOCATION error for Sps_loop ..'
        PRINT *,'   STAT =',ier
        Return
      ENDIF
      call path_sweep(brkpt, no_ele, kk, ngp1, z_path(k), &
             &  f%template%surfs(:,1), sps_zeta_loop(j,k)%values)
      call path_sweep(brkpt, no_ele, jp, ngp1, phi_path(k), &
             &  deg2rad*f%template%phi(1,:), sps_phi_loop(j,k)%values)
    end do
  end do

  ! Compute the relative refractive index minus one.
  ! Get the water mixing ratio function

  h2o => GetVectorQuantityByType( fwdModelIn, fwdModelExtra, &
    & quantityType=l_vmr, molecule=l_h2o, noError=.true.)
  if ( associated(h2o) ) then
    jp = h2o%template%noInstances
    kk = h2o%template%noSurfs
  else
    jp = 0
    kk = 0
  end if

  call refractive_index ( h2o%values, h2o%template%surfs(:,1), &
    &  h2o%template%phi(1,:), kk, jp, ndx_path, z_path, t_path, &
    &  phi_path, n_path, associated(h2o), no_tan_hts )

  Return

!---------------------------------------------------------------------------
contains
  ! *****     Internal procedures     **********************************
!
!-------------------------------------------------------------------
 Subroutine path_sweep(brkpt, no_ele, no_coeff, Ngp1, v_path, &
                    &  basis, values)

  Integer(i4), intent(in) :: brkpt, no_ele, no_coeff, Ngp1
  Integer(i4), intent(out) :: values(:,:)

  Real(r8), intent(in) :: basis(:)
  Type(path_vector), INTENT(IN) :: v_path

  Integer(i4) :: i, mp, klo, khi, kk
  Real(r8) :: q, v
!
   values(:,1) = 0
   values(:,2) = 0

   i = 0
   mp = 1 - Ngp1
   do
     mp = mp + Ngp1
     if(mp > brkpt) EXIT
     klo = 0
     khi = 0
     i = i + 1
     v = v_path%values(mp)
     do kk = 1, no_coeff
       Call get_one_eta(v,basis,no_coeff,kk,q)
       if(q /= 0.0) then
         khi = kk + 1
         if(klo == 0) klo = kk
       endif
     end do
     values(i,1) = max(1,klo)
     values(i,2) = min(no_coeff,khi+1)
   end do
   i = i - 1
   mp = brkpt + 1 - Ngp1
   do
     mp = mp + Ngp1
     if(mp >= no_ele) EXIT
     klo = 0
     khi = 0
     i = i + 1
     v = v_path%values(mp)
     do kk = 1, no_coeff
       Call get_one_eta(v,basis,no_coeff,kk,q)
       if(q /= 0.0) then
         khi = kk + 1
         if(klo == 0) klo = kk
       endif
     end do
     values(i,1) = max(1,klo)
     values(i,2) = min(no_coeff,khi+1)
   end do
!
   Return

 End subroutine path_sweep

 End subroutine Get_path_spsfunc_ngrid

end module GET_PATH_SPSFUNC_NGRID_M
! $Log$
! Revision 1.5  2001/05/03 22:23:29  vsnyder
! Insert copyright notice, some cosmetic changes
!
! Revision 1.4  2001/04/24 00:04:35  livesey
! Removed ancient print statement
!
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
