module ZATMOS_DERIV_M
  use GET_DRAD_NOTDER_M, only: GET_DRAD_NOTDER
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_DERIVATIVE
  implicit NONE
  private
  public :: ZATMOS_DERIV
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Subroutine ZATMOS_DERIV(atmospheric, n_sps, band, frq_i,no_coeffs_f, &
 &           mid, delta, t_script, tau, ilo, ihi, no_phi_f, k_atmos,Ier)
!
    type(atmos_comp), intent(in) :: ATMOSPHERIC(:)

    integer(i4), intent(in) :: N_SPS, BAND, mid, ILO, IHI, frq_i
    integer(i4), intent(in) :: NO_COEFFS_F(:), NO_PHI_F(:)

    real(r8), intent(in) :: DELTA(:,:,:,:)
    real(r8), intent(in) :: T_SCRIPT(:), TAU(:)

    integer(i4), intent(out) :: Ier

    Type(path_derivative), intent(in out) :: k_atmos(:)

    Real(r8) :: r
    Integer(i4) :: J, IZ, IP, NO_Z, NO_PHI
!
    Ier = 0
    do j = 1, n_sps
!
      IF(atmospheric(j)%der_calc(band)) THEN
!
        no_phi = no_phi_f(j)
        no_z = no_coeffs_f(j)
!
        do ip = 1, no_phi
!
          do iz = 1, no_z
!
            Call get_drad_notder (delta(1:,iz,ip,j), t_script, tau, &
   &                              mid, ilo, ihi, r)
            k_atmos(j)%values(frq_i,iz,ip) = r
!
          end do
!
        end do
!
      ENDIF
!
    end do
!
    Return
  End Subroutine ZATMOS_DERIV
end module ZATMOS_DERIV_M
