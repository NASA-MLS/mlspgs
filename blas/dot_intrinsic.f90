! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
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
    module procedure DDOT, SDOT
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

end module Dot_M

! $Log$
