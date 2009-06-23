! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DREXP_M

! ----------------------------------------------------------------------
!            EVALUATE THE FUNCTION EXP(X) - 1
!-----------------------------------------------------------------------

  implicit NONE
  private
  public :: REXP, DREXP

  interface REXP; module procedure DREXP; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  function DREXP ( X ) result ( REXP )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X
    real(rk) :: REXP
    include 'rexp.f9h'
  end function DREXP

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module DREXP_M

! $Log$
! Revision 2.5  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.4  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.3  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.2.2.1  2003/04/16 20:12:33  vsnyder
! Move argument declarations from rexp.f9h to here
!
! Revision 2.2  2002/10/10 20:36:50  vsnyder
! Darn, I forgot the CVS log comment
!
