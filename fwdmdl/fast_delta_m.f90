module FAST_DELTA_M
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA
  use GENERIC_DELTA_INTEGRAL_M, only: GENERIC_DELTA_INTEGRAL
  implicit NONE
  private
  public :: FAST_DELTA
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrate in ZETA Space !
!
  Subroutine FAST_DELTA(mid,brkpt,no_ele,z_path,h_path,phi_path, &
 &           beta_path,dHdz_path,spsfunc_path,n_sps,N_lvls,      &
 &           Nlvl,ref_corr,delta,Ier)
!
    Integer(i4), intent(in) :: NLVL, N_SPS, N_LVLS, MID, BRKPT, NO_ELE

    Real(r8), intent(in) :: REF_CORR(:)

    Real(r8), intent(inout) :: DELTA(:,:)     ! (N2lvl,Nsps)

    Integer(i4), intent(out) :: IER

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, DHDZ_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
!
! -----     Local variables     ----------------------------------------
!
    Integer(i4) :: J, K
!
    Integer(i4), SAVE :: NZ=1, NP=1, IZ=1, IP=1
!
    Real(r8), SAVE :: FQ = 1.0
    Real(r8), SAVE, DIMENSION(2) :: Z_BASIS = (/ -5.0, 5.0/)
    Real(r8), SAVE, DIMENSION(2) :: PHI_BASIS = (/-50.0, 50.0/)
!
    Real(r8), ALLOCATABLE, DIMENSION(:) :: Integrand
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
    DEALLOCATE(Integrand, STAT=k)
!
    ALLOCATE(Integrand(no_ele), STAT=ier)
    IF(ier /= 0) THEN
      Ier = 1
      Print *,'** Error: ALLOCATION error in FAST_DELTA routine ..'
      goto 99
    endif
!
    do j = 1, n_sps
!
      integrand(1:no_ele) = spsfunc_path(j)%values(1:no_ele) *  &
     &                      beta_path(j)%values(1:no_ele)
!
      Call generic_delta_integral(mid, brkpt, no_ele, z_path, h_path, &
 &         phi_path, dhdz_path, N_lvls, ref_corr, integrand, z_basis, &
 &         phi_basis, nz, np, iz, ip, fq, delta(1:,j), Ier)
      IF(ier /= 0) goto 99

    end do
!
 99  DEALLOCATE(Integrand, STAT=k)

    Return

  End Subroutine FAST_DELTA
!
end module FAST_DELTA_M
! $Log$
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:13  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
