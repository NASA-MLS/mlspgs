module ZATMOS_DERIV_M
  use GET_DRAD_M, only: GET_DRAD
  use L2PCdim, only: MAX_NO_PHI, N2LVL, NLVL
  use L2PC_FILE_PARAMETERS, only: MAX_NO_ELMNTS_PER_SV_COMPONENT
  use L2PC_PFA_STRUCTURES, only: ATMOS_COMP
  use MLSCommon, only: I4, R4
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

  Subroutine ZATMOS_DERIV ( atmospheric, ptg_i, sps_tbl, n_sps, band,   &
 &           no_coeffs_f, IndxR, IndxL, delta, t_script, tau, ilo, ihi, &
 &           no_phi_f, k_star_atmos )
!
    Integer(i4), Parameter :: mnp=max_no_phi
    Integer(i4), Parameter :: mxco=max_no_elmnts_per_sv_component

    type(atmos_comp), intent(in) :: ATMOSPHERIC(*)
    integer(i4), intent(in) :: PTG_I, SPS_TBL(*), N_SPS, BAND
    integer(i4), intent(in) :: NO_COEFFS_F(*), IndxR, IndxL
    real(r4), intent(in) :: DELTA(N2lvl,mxco,mnp,*), T_SCRIPT(*), TAU(*)
    integer(i4), intent(in) :: ILO, IHI, NO_PHI_F(*)
    real(r4), intent(out) :: K_STAR_ATMOS(Nlvl,mxco,mnp,*)

    Integer(i4) :: SPS_I, SV_I, NJ, NF, K, IS
    Real(r4), save :: Ary_zero(N2lvl) = 0.0
!
    do sps_i = 1, n_sps
!
      is = sps_tbl(sps_i)
      if ( atmospheric(is)%der_calc(band) ) then
!
        nj = no_coeffs_f(is)
        nf = no_phi_f(sps_i)
!
        do k = 1, nf
!
          do sv_i = 1, nj
!
            Call get_drad ( delta(1, sv_i, k, is), t_script, tau, Ary_zero, &
   &                        IndxR, IndxL, ilo, ihi,                         &
   &                        k_star_atmos(ptg_i,sv_i,k,is) )
          end do
!
        end do
!
      endif
!
    end do
!
    Return
  End Subroutine ZATMOS_DERIV
end module ZATMOS_DERIV_M
