! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Dot_M

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
! dot         Returns dot product of two vectors
! === (end of toc) ===

! === (start of api) ===
! value dot ( int n, val1 x(:), int incx, val2 y(:), int incy )
!      where val1, val2 can be any of the types:
!      { real, double precision }
! === (end of api) ===
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
contains 
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Dot_M

! $Log$
! Revision 1.7  2005/06/22 19:45:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.6  2002/10/10 23:48:57  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.5  2002/09/13 22:50:51  pwagner
! Change external names to MSDDOT and MDSDOT; removed pointer copies
!
! Revision 1.4  2002/09/13 18:04:31  pwagner
! Added mixed type dot products: a_r4 . b_r8 and a_r8 . b_r4
!
! Revision 1.3  2001/11/08 02:37:50  vsnyder
! Make real argumemnts of ?DOT scalar
!
! Revision 1.2  2001/11/03 00:56:36  pwagner
! Added Copyright statement
!
