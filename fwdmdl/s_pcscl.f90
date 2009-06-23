! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module S_PCSPL_M

  implicit NONE

  private
  public :: S_PCSPL, PCSPL

  interface PCSPL; module procedure S_PCSPL; end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine S_PCSPL ( TAU, C, N, IBCBEG, IBCEND )
    integer, parameter :: RK = kind(0.0e0)
    include 'pcspl.f9h'
  end subroutine S_PCSPL

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module S_PCSPL_M

! $Log$
! Revision 2.3  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2005/06/22 18:08:20  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
