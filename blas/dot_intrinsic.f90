! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Dot_M

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  interface DOT
    module procedure DDOT, SDOT, MSDDOT, MDSDOT
  end interface

contains

  double precision function DDOT ( N, X, INCX, Y, INCY )
    integer, intent(in) :: N, INCX, INCY
    double precision, intent(in) :: X(*), Y(*)
    ddot = dot_product(x(1:1+(n-1)*incx:incx), y(1:1+(n-1)*incy:incy))
  end function DDOT

  real function SDOT ( N, X, INCX, Y, INCY )
    integer, intent(in) :: N, INCX, INCY
    real, intent(in) :: X(*), Y(*)
    sdot = dot_product(x(1:1+(n-1)*incx:incx), y(1:1+(n-1)*incy:incy))
  end function SDOT

  double precision function MSDDOT ( N, X, INCX, Y, INCY )
    integer, intent(in) :: N, INCX, INCY
    real, intent(in) :: X(*)
    double precision, intent(in) :: Y(*)
    ! double precision, dimension(:), pointer :: dbleX
    if ( n <= 0 ) then
      sddot = 0.d0
    elseif ( n == 1 ) then
      sddot = x(1) * y(1)
    else
      ! allocate(dbleX(n))
      ! dbleX = x(1:1+(n-1)*incx:incx)
      ! sddot = ddot(n, dbleX, 1, y, incy)
      ! deallocate(dbleX)
      msddot = dot_product(x(1:1+(n-1)*incx:incx), y(1:1+(n-1)*incy:incy))
    endif
  end function MSDDOT

  double precision function MDSDOT ( N, X, INCX, Y, INCY )
    integer, intent(in) :: N, INCX, INCY
    double precision, intent(in) :: X(*)
    real, intent(in) :: Y(*)
    ! mdsdot = msddot( n, y, incy, x, incx)
    mdsdot = dot_product(x(1:1+(n-1)*incx:incx), y(1:1+(n-1)*incy:incy))
  end function MDSDOT

end module Dot_M

! $Log$
! Revision 1.3  2002/09/13 18:04:31  pwagner
! Added mixed type dot products: a_r4 . b_r8 and a_r8 . b_r4
!
! Revision 1.2  2001/11/03 00:56:36  pwagner
! Added Copyright statement
!
