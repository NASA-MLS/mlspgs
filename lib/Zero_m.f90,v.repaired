! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ZERO_M

  ! Univariate zero finder from Math 77.

  implicit NONE

  private

  public DZERO, SZERO, ZERO

  interface ZERO
    module procedure DZERO, SZERO
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DZERO ( X1, F1, X2, F2, MODE, TOL )
    use Mess_M, only: Mess
    integer, parameter :: RK = kind(0.0d0)
    include 'zero.f9h'
  end subroutine DZERO

  subroutine SZERO ( X1, F1, X2, F2, MODE, TOL )
    use Mess_M, only: Mess
    integer, parameter :: RK = kind(0.0e0)
    include 'zero.f9h'
  end subroutine SZERO

!------------------------------------------------------------------

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ZERO_M

! $Log$
! Revision 2.2  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.1  2007/09/07 01:35:04  vsnyder
! Initial commit
!
