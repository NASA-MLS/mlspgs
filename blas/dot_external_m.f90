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
    real function SDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      real, intent(in) :: X, Y
    end function SDOT
    double precision function DDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      double precision, intent(in) :: X, Y
    end function DDOT
    double precision function MSDDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      real, intent(in) :: X
      double precision, intent(in) :: Y
    end function MSDDOT
    double precision function MDSDOT ( N, X, INCX, Y, INCY )
      integer, intent(in) :: N, INCX, INCY
      double precision, intent(in) :: X
      real, intent(in) :: Y
    end function MDSDOT
  end interface
end module Dot_M

! $Log$
! Revision 1.4  2002/09/13 18:04:31  pwagner
! Added mixed type dot products: a_r4 . b_r8 and a_r8 . b_r4
!
! Revision 1.3  2001/11/08 02:37:50  vsnyder
! Make real argumemnts of ?DOT scalar
!
! Revision 1.2  2001/11/03 00:56:36  pwagner
! Added Copyright statement
!
