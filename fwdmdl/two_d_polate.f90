module TWO_D_POLATE_M
  use MLSCommon, only: I4, R4
  use S_GET_ONE_ETA_M, only: GET_ONE_ETA
  implicit NONE
  private
  public :: TWO_D_POLATE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains
!  Interpolating in 2 dimensions - Zeta & Phi
!
  Subroutine TWO_D_POLATE ( f_basis, f_profile, Nc, no_coeffs, phi_basis, &
 &                          no_phi, z, phi, V )
!
    Real(r4), intent(in) :: F_BASIS(*)
    integer(i4), intent(in) :: NC
    Real(r4), intent(in) :: F_PROFILE(Nc,*)
    integer(i4), intent(in) :: NO_COEFFS
    Real(r4), intent(in) :: PHI_BASIS(*)
    integer(i4), intent(in) :: NO_PHI
    Real(r4), intent(in) :: Z
    Real(r4), intent(in) :: PHI
    Real(r4), intent(out) :: V
! -----     Local variables     ----------------------------------------
    integer :: J, K
    Real(r4) :: ETAP, ETAZ
! -----     Executable statements     ----------------------------------
    V = 0.0
    do k = 1, no_coeffs
      Call get_one_eta(z,f_basis,no_coeffs,k,etaz)
      if (etaz > 0.0) then
        do j = 1, no_phi
          Call get_one_eta(phi,phi_basis,no_phi,j,etap)
          V = V + f_profile(k,j) * etaz * etap
        end do
      end if
    end do
!
    Return
  End Subroutine TWO_D_POLATE
end module TWO_D_POLATE_M

! $Log$
