module SCRT_DN_M
  use MLSCommon, only: I4, R4, R8
  implicit NONE
  private
  public :: SCRT_DN
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
  Subroutine Scrt_Dn ( t_script, N_lvls, earth_ref, del_opacity, tau, &
 &                     radiance, IndxR, IndxL, ilo, ihi )
!
!     SCRT = Scalar Condensed Radiative Transfer
!
! This subroutine computes the radiative transfer using a
! condensed algorithm. Enter t_script(*) in a quantity that is linear
! in power.
!
    real(r8), intent(in) :: T_SCRIPT(*)
    integer(i4), intent(in) :: N_lvls
    real(r8), intent(in) :: EARTH_REF, DEL_OPACITY(*)
    real(r8), intent(out) :: TAU(*), RADIANCE
    integer(i4), intent(in) :: IndxR, IndxL
    integer(i4), intent(out) :: ILO, IHI
    Real(r8), Parameter :: cut_off = -15.0
    Integer(i4) i, k
    Real(r8) :: r, total_opacity
!
! First segment initialization (From the spacecraft to the tangent point)
!
    ihi = 0
    ilo = 2
    k = 2 * N_lvls
!
    tau(1) = 1.0
    total_opacity = 0.0
    radiance = t_script(1)
!
    i = 2
    do while (total_opacity > cut_off .and. i <= IndxR)
      total_opacity = total_opacity - del_opacity(i-1)
      r = exp(total_opacity)
      tau(i) = r
      radiance = radiance + t_script(i) * r
      ilo = i
      i = i + 1
    end do
!
    if (total_opacity <= cut_off) Return
!
    tau(IndxL) = earth_ref * tau(IndxR)
    radiance = radiance + t_script(IndxL) * tau(IndxL)
!
    i = IndxL + 1
    ihi = i
!
    do while (total_opacity > cut_off .and. i <= k)
      total_opacity = total_opacity - del_opacity(i-2)
      r = earth_ref * exp(total_opacity)
      tau(i) = r
      radiance = radiance + t_script(i) * r
      ihi = i
      i = i + 1
    end do
!
    Return
  End Subroutine Scrt_Dn
end module SCRT_DN_M
! $Log$
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
