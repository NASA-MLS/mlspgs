!
module TRAPEZ_INT_M
  use MLSCommon, only: I4, R4, R8
  Implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  Trapezoidal integration routine.

SUBROUTINE trapez_int(n,x,f,r)

INTEGER(i4), INTENT(IN) :: n

REAL(r8), INTENT(IN) :: x(*), f(*)

REAL(r8), INTENT(OUT) :: r

INTEGER(i4) :: i
REAL(r8) :: xi, zi, xim1, zim1

r = 0.0D0
IF(n < 2) RETURN

  xi = x(1)
  zi = f(1)
  DO i = 2, n
    xim1 = xi
    xi = x(i)
    zim1 = zi
    zi = f(i)
    r = r + 0.5_r8 * (zi + zim1) * (xi - xim1)
  END DO

  RETURN
END SUBROUTINE trapez_int
end module TRAPEZ_INT_M
! $Log$
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
