module GET_DRAD_NOTDER_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: GET_DRAD_NOTDER
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName = &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This subroutine computes the scalarized condensed radiative transfer
! derivatives, without Temperature derivatives
!
  Subroutine GET_DRAD_NOTDER (DDER, T_SCRIPT, TAU, mid, ILO, IHI, RDRAD)
!
    Integer(i4), intent(in) :: mid, ILO, IHI
    Real(r8), intent(in) :: DDER(:), T_SCRIPT(:), TAU(:)
    Real(r8), intent(out) :: RDRAD

    Integer(i4) :: h_i, IndxL
    Real(r8) :: drad, q, w
!
    w = 0.0_r8
    drad = 0.0_r8
!
    h_i = 2
    do while (h_i <= ilo)
      w = w - dder(h_i-1)
      q = t_script(h_i) * w
      drad = drad + q * tau(h_i)
      h_i = h_i + 1
    end do
!
    if (ilo > ihi) then
      rdrad = drad
      Return
    endif
!
    IndxL = mid + 1
    q = t_script(IndxL) * w
    drad = drad + q * tau(IndxL)
!
    h_i = IndxL + 1
    do while (h_i <= ihi)
      w = w - dder(h_i-2)
      q = t_script(h_i) * w
      drad = drad + q * tau(h_i)
      h_i = h_i + 1
    end do
!
    rdrad = drad
!
    Return
  End Subroutine GET_DRAD_NOTDER
end module GET_DRAD_NOTDER_M
! $Log$
! Revision 1.1  2001/03/31 18:12:05  Z.Shippony
! Initial conversion to Fortran 90
