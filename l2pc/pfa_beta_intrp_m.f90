module PFA_BETA_INTRP_M
  use MLSCommon, only: I4, R4, R8
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public :: PFA_BETA_INTRP
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!--------------------------------------------------------------------
! Interpolate the beta using the Power interpolation model
!
  Subroutine PFA_BETA_INTRP ( z, t, z_path, t_path, ihmin, ihmax, n1, n2, &
 &           Ngp1, pfa_beta_coeff, beta_t_power, pfa_beta_Ntrp, t_power )
    real(r8), intent(in) :: Z
    real(r8), intent(in) :: T
    real(r8), intent(in) :: Z_PATH(*)
    real(r8), intent(in) :: T_PATH(*)
    integer(i4), intent(in) :: IHMIN, IHMAX
    integer(i4), intent(in) :: N1, N2
    integer(i4), intent(in) :: NGP1
    real(r8), intent(in) :: PFA_BETA_COEFF(*)
    real(r8), intent(in) :: BETA_T_POWER(*)
    real(r8), intent(out) :: PFA_BETA_NTRP
    real(r8), intent(out) :: T_POWER
!
    Integer(i4) :: IH, J, JH, L, SKNG, ZKLO
    Real(r8) :: D, Q, R, VP, VL, W
!
!  Find the closest zeta (Log(pressure)) index along the path:
!
    l = n1-1
    Call Hunt(z,z_path(n1:n2),n2-n1+1,zklo,j)
    if (abs(z-z_path(l+j)) < abs(z-z_path(l+zklo))) zklo = j
!
    zklo = zklo + l
!
    l = (zklo-n1)/Ngp1
    ih = ihmin + l
    skng = n1 + l * Ngp1
!
    t_power = beta_t_power(ih)
    Pfa_beta_Ntrp = pfa_beta_coeff(ih)
    if (Pfa_beta_Ntrp == 0.0) Return
!
    r = Pfa_beta_Ntrp
!
    j = min(ih+1,ihmax)
    if (j == ih) j = ih-1
    jh = skng + (j - ih) * Ngp1
!
!  Interpolation in temperature
!
    vl = t / t_path(skng)
    if (abs(vl-1.0_r8) >= 1.0e-3_r8) &  ! ** Temperature correction needed ..
   &         Pfa_beta_Ntrp = Pfa_beta_Ntrp * (vl**t_power)
!
!  Interpolation in zeta
!
    d = z - z_path(skng)
    if (abs(d) >= 1.0e-3_r8) then       ! ** Pressure correction needed ..
!
      vl = t_path(skng) / t_path(jh)
      vp = pfa_beta_coeff(j) * (vl**t_power)
      q = vp / r
      w = abs(q-1.0_r8)
      if (w > 1.0e-7_r8) then           ! Power rule for pressure
        w = d / (z_path(jh) - z_path(skng))
        Pfa_beta_Ntrp = Pfa_beta_Ntrp * (q**w)
      end if
!
    end if
!
    Return
  End Subroutine PFA_BETA_INTRP
end module PFA_BETA_INTRP_M
! $Log$
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
