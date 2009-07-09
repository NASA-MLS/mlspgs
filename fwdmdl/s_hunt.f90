! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module S_HUNT_M
  implicit NONE
  private
  public :: S_HUNT, HUNT
  interface HUNT; module procedure S_HUNT; end interface
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
! A binary search routine with a hunt procedure, to start from last known
! location (if 0 < JLO < N) or from the begining otherwise.
  pure subroutine S_HUNT ( ELEMENT, ARRAY, N, JLO, JHI )
    integer, parameter :: RK = kind(0.0e0)
    include 'hunt.f9h'
  end subroutine S_HUNT
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module S_HUNT_M

! $Log$
! Revision 2.5  2009/07/09 23:54:38  vsnyder
! Remove need for MLSCommon or MLSKinds
!
! Revision 2.4  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2005/03/29 01:58:17  vsnyder
! Make stuff pure
!
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
