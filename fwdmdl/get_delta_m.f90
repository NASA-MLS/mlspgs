! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GET_DELTA_M
  use MLSCommon, only: I4, R8
  use ELLIPSE_M, only: ELLIPSE
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA, PATH_INT_VECTOR_2D
  use GENERIC_DELTA_INTEGRAL_M, only: GENERIC_DELTA_INTEGRAL
  use ForwardModelConfig, only: ForwardModelConfig_T
  use VectorsModule, only: Vector_T, VectorValue_T, GetVectorQuantityByType
  use Intrinsic, only: l_vmr
  use Units, only: Deg2Rad
  implicit NONE
  private
  public :: GET_DELTA
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrate in ZETA Space !
!
  Subroutine GET_DELTA(ForwardModelConfig, FwdModelExtra, FwdModelIn,    &
 &           mid,brkpt,no_ele,z_path,h_path,phi_path,beta_path,dHdz_path,&
 &           n_sps,N_lvls,ref_corr,spsfunc_path,elvar,midval_ndx,        &
 &           no_midval_ndx,gl_ndx,no_gl_ndx,Sps_zeta_loop,Sps_phi_loop,  &
 &           midval_delta,delta,Ier)
!
    type(forwardModelConfig_T), intent(in) :: forwardModelConfig
    type (Vector_T), intent(in) :: fwdModelIn, fwdModelExtra
!
    Integer(i4), intent(in) :: n_sps, N_LVLS
    Integer(i4), intent(in) :: mid, brkpt, no_ele
!
    Integer(i4), intent(in) :: gl_ndx(:,:),no_gl_ndx
    Integer(i4), intent(in) :: midval_ndx(:,:),no_midval_ndx

    Real(r8), intent(in) :: REF_CORR(:)
    Real(r8), intent(in) :: midval_delta(:,:)

    Real(r8), intent(inout) :: DELTA(:,:,:,:)

    Integer(i4), intent(out) :: IER

    Type(ELLIPSE), intent(in out) :: elvar

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, DHDZ_PATH

    Type(path_int_vector_2d), intent(in) :: SPS_PHI_LOOP(:), &
                                        &   SPS_ZETA_LOOP(:)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: j, ip, iz, nco, npf

    Real(r8) :: Q
    Real(r8), allocatable, dimension(:) :: Integrand

    type (VectorValue_T), pointer :: f
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
    allocate ( Integrand(no_ele), STAT=ier )
    IF(ier /= 0) THEN
      Ier = 1
      Print *,'** Error: ALLOCATION error in GET_DELTA routine ..'
      goto 99
    end if
!
! Start delta array computations:
!
    delta = 0.0
    do j = 1, n_sps
!
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
     &     quantityType=l_vmr, molecule=forwardModelConfig%molecules(j))

      IF (.not. forwardModelConfig%moleculeDerivatives(j)) CYCLE
!
      nco = f%template%noSurfs
      npf = f%template%noInstances
!
      integrand(1:no_ele) = beta_path(j)%values(1:no_ele)
!
      if (f%template%logBasis) integrand(1:no_ele) = &
     &           integrand(1:no_ele) * spsfunc_path(j)%values(1:no_ele)
!
! Loop over the specie's Phi's
!
      do ip = 1, npf
!
! Loop over the specie's zeta's
!
        do iz = 1, nco
!
          q = 1.0
          if (f%template%logBasis) q = 1.0 / f%values(iz,ip)
!
          Call generic_delta_integral(mid,brkpt,no_ele,z_path,      &
         &     h_path,phi_path,dhdz_path,N_lvls,ref_corr,integrand, &
         &     f%template%surfs(:,1),Deg2Rad*f%template%phi(1,:),   &
         &     nco,npf,iz,ip,q,elvar,midval_delta(1:,j),midval_ndx, &
         &     no_midval_ndx,gl_ndx,no_gl_ndx,Sps_zeta_loop(j),     &
         &     Sps_phi_loop(j),1,delta(1:,iz,ip,j),Ier)
          IF(ier /= 0) goto 99
!
        end do
!
      end do
!
    end do
!
 99  deallocate ( Integrand, STAT=j )
!
    Return
!
  End Subroutine GET_DELTA
!
end module GET_DELTA_M
! $Log$
! Revision 1.9  2001/05/03 22:19:36  vsnyder
! Insert copyright notice, some cosmetic changes
!
! Revision 1.8  2001/04/09 20:52:07  zvi
! Debugging Derivatives version
!
! Revision 1.7  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.6  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.5  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
