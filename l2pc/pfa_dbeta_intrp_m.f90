module PFA_DBETA_INTRP_M
  use MLSCommon, only: I4, R4, R8
  use D_HUNT_M, only: HUNT
  implicit NONE
  private
  public :: PFA_DBETA_INTRP
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!--------------------------------------------------------------------
! Interpolate the beta derivatives w.r.t. to a spectroscopic (w,n or: nu0)
! using the Power interpolation model (for temperature variations) and
! Linear interpolation (for pressure variations)
!
  Subroutine PFA_DBETA_INTRP ( z, t, z_path, t_path, ihmin, ihmax, n1, n2, &
 &           Ngp1, beta_t_power, pfa_dbeta_ds, dbeta_ds_Ntrp )
!
    real(r8), intent(in) :: Z
    real(r8), intent(in) :: T
    real(r8), intent(in) :: Z_PATH(*)
    real(r8), intent(in) :: T_PATH(*)
    integer(i4), intent(in) :: IHMIN, IHMAX
    integer(i4), intent(in) :: N1, N2
    integer(i4), intent(in) :: NGP1
    real(r8), intent(in) :: BETA_T_POWER(*)
    real(r8), intent(in) :: PFA_DBETA_DS(*)
    real(r8), intent(out) :: DBETA_DS_NTRP
    integer(i4) :: IH, IK, J, JZ, L, M
    real(r8) :: Q, R, V, VL
    integer(i4) :: SKNG
    real(r8) :: T_P
!
!  Find the closest zeta (Log(pressure)) index along the path:
!
    l = n1-1
    m = n2-n1+1
    Call Hunt(z,z_path(n1:n2),m,ik,j)
    if (abs(z-z_path(l+j)) < abs(z-z_path(l+ik))) ik = j
!
    ik = ik + l
    l = (ik-n1)/Ngp1
    ih = ihmin + l
    skng = n1 + l * Ngp1
!
    t_p = beta_t_power(ih)
!
!  For the interpolation in temperature (Use Power rule eqn.)
!
    q = 1.0
    vl = t / t_path(skng)
    if (abs(vl-1.0) >= 1.0d-4) q = vl**t_p     ! Temp. Power rule
!
!  Interpolation in zeta (Using Linear interpolation scheme)
!
    v = abs(z_path(skng)-z)
!
    if (v < 1.0d-3) then             ! ** Zeta correction NOT needed ..
!
      dbeta_ds_Ntrp = pfa_dbeta_ds(ih) * q     ! Temp. Power rule
!
    else                             ! ** Zeta correction needed ..
!
      if (ih > ihmin) then
        j = ih - 1
      else
        j = ih + 1
      end if
!
      l = j - ihmin
      jz = n1 + l * Ngp1
      v = (z - z_path(skng)) / (z_path(jz) - z_path(skng))
!
      r = pfa_dbeta_ds(ih) * (1.0 - v) + pfa_dbeta_ds(j) * v
      dbeta_ds_Ntrp = r * q          ! Apply Temp. corrections
!
    end if
!
      Return
  End Subroutine PFA_DBETA_INTRP
end module PFA_DBETA_INTRP_M
! $Log$
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
