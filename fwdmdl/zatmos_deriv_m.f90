module ZATMOS_DERIV_M
  use GET_DRAD_M, only: GET_DRAD
  use L2PCdim, only: MAX_NO_PHI, N2LVL, Nptg
  use L2PC_FILE_PARAMETERS, only: MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP
  use MLSCommon, only: I4, R4, R8
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
  Subroutine ZATMOS_DERIV(atmospheric, ptg_i, sps_tbl, n_sps, band,      &
 &           no_coeffs_f ,mid, delta, t_script, tau, ilo, ihi, no_phi_f, &
 &           k_atmos)
!
    Integer(i4), Parameter :: mnp=max_no_phi
    Integer(i4), Parameter :: mxco=max_no_elmnts_per_sv_component

    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)

    integer(i4), intent(in) :: PTG_I, N_SPS, BAND, mid, ILO, IHI
    integer(i4), intent(in) :: NO_COEFFS_F(*), SPS_TBL(*), NO_PHI_F(*)

    real(r8), intent(in) :: DELTA(N2lvl,mxco,mnp,*), T_SCRIPT(*), TAU(*)

    real(r4), intent(out) :: K_ATMOS(Nptg,mxco,mnp,*)

    Integer(i4) :: SPS_I, SV_I, NJ, NF, K, J

    Real(r8) :: r
    Real(r8), SAVE :: Ary_zero(N2lvl) = 0.0_r8
!
    do sps_i = 1, n_sps
!
      j = sps_tbl(sps_i)
      IF(atmospheric(j)%der_calc(band)) THEN
!
        nj = no_coeffs_f(j)
        nf = no_phi_f(sps_i)
!
        do k = 1, nf
!
          do sv_i = 1, nj
!
            Call get_drad (delta(1:,sv_i,k,j), t_script, tau, Ary_zero, &
   &                       mid, ilo, ihi, r)
            k_atmos(ptg_i,sv_i,k,j) = r
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
