! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TOGGLES_LR

  implicit NONE
  private

  integer, parameter, public :: MAX_TOGGLE = 255
  integer, public :: TOGGLE(0:MAX_TOGGLE) = 0 ! Indexed by IACHAR(...)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

end module TOGGLES_LR

! $Log$
! Revision 1.2  2013/10/26 00:11:53  vsnyder
! Add CVS stuff
!
