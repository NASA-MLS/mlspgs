module SCRT_DN_M
  use MLSCommon, only: I4, R8
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
  Subroutine Scrt_Dn (t_script, N_lvls, earth_ref, del_opacity, tau, &
 &                    radiance, mid, ilo, ihi)
!
!     SCRT = Scalar Condensed Radiative Transfer
!
! This subroutine computes the radiative transfer using a
! condensed algorithm. Enter t_script(*) in a quantity that is linear
! in power.
!
    integer(i4), intent(in) :: N_lvls
    integer(i4), intent(in) :: mid

    real(r8), intent(in) :: T_SCRIPT(*)
    real(r8), intent(in) :: EARTH_REF, DEL_OPACITY(*)

    real(r8), intent(out) :: TAU(*), RADIANCE
    integer(i4), intent(out) :: ILO, IHI

    Real(r8), Parameter :: cut_off = -15.0

    Integer(i4) :: i, j, k
    Real(r8)    :: r, total_opacity
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
    do while (total_opacity > cut_off .and. i <= mid)
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
    i = mid + 1
    tau(i) = earth_ref * tau(mid)
    radiance = radiance + t_script(i) * tau(i)
!
    i = mid + 2
    ihi = i
!
    j = k + 1
    do 
      j = j - 1
      if(t_script(j) /= 0.0) EXIT 
      if(del_opacity(j-2) /= 0.0) EXIT 
      ihi = j - 1
    end do
!
    j = ihi
    ihi = i
    do while (total_opacity > cut_off .and. i <= j)
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
