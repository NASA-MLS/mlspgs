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
    use MLSCommon, only: I4, R4
    integer, parameter :: RK = r4
    include 'hunt.f9h'
  end subroutine S_HUNT
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module S_HUNT_M

! $Log$
! Revision 2.3  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2005/03/29 01:58:17  vsnyder
! Make stuff pure
!
! Revision 2.1  2003/11/08 02:01:01  vsnyder
! Initial commit -- created just for symmetry
!
