module PFA_DB_DELTA_M
  use MLSCOmmon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR
  use GENERIC_DELTA_INTEGRAL_M, only: GENERIC_DELTA_INTEGRAL
  implicit NONE
  private
  public :: PFA_DB_DELTA
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta_db function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrates in ZETA Space !
!
  Subroutine PFA_DB_DELTA (mid, brkpt, no_ele, z_path, h_path, phi_path, &
 &           dhdz_path, N_lvls, ref_corr, spsfunc_s, pfa_dbeta_s,        &
 &           z_basis, phi_basis, nz, np, iz, ip, delta_s, Ier )
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, NO_ELE, NZ, NP, &
   &                           IZ, IP
!
    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, DHDZ_PATH
!
    Real(r8), intent(in) :: REF_CORR(*)
    Real(r8), intent(in) :: SPSFUNC_S(*)
    Real(r8), intent(in) :: PFA_DBETA_S(*)

    Real(r8), intent(in) :: Z_BASIS(:)
    Real(r8), intent(in) :: PHI_BASIS(:)

    Real(r8), intent(out) :: DELTA_S(*)

    Integer(i4), intent(out) :: IER
!
!  Local variables:
!
    Integer(i4) :: I
!
    Real(r8) :: Q
    Real(r8), ALLOCATABLE, DIMENSION(:) :: Integrand
!
! Start code:
!
    Ier = 0
    DEALLOCATE(Integrand, STAT=i)
!
    ALLOCATE(Integrand(no_ele), STAT=ier)
    IF(ier /= 0) THEN
      Ier = 1
      Print *,'** Error: ALLOCATION error in PFA_DB_DELTA routine ..'
      goto 99
    endif
!
!  Initialize all arrays:
!
    i = 2 * N_lvls
    delta_s(1:i) = 0.0
!
    q = 1.0
    integrand(1:no_ele) = pfa_dbeta_s(1:no_ele) * spsfunc_s(1:no_ele)
!
    Call generic_delta_integral(mid, brkpt, no_ele, z_path, h_path, &
   &     phi_path, dhdz_path, N_lvls, ref_corr, integrand, z_basis, &
   &     phi_basis, nz, np, iz, ip, q, delta_s, Ier )
!
 99  DEALLOCATE(Integrand, STAT=i)
!
    Return
!
  End Subroutine PFA_DB_DELTA
!
end module PFA_DB_DELTA_M
! $Log$
! Revision 1.1  2000/06/21 21:56:16  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
