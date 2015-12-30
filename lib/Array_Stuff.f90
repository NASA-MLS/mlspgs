! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Array_Stuff
!=============================================================================

  ! Compute array element order position from subscripts, or vice versa.

  implicit NONE
  private

  public :: Element_Position, Subscripts

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  pure integer function Element_Position ( Lbounds, Ubounds, Subscripts )
    ! Compute position in array element order from subscripts.
    integer, intent(in), contiguous :: Lbounds(:)
    integer, intent(in) :: Ubounds(size(lbounds))
    integer, intent(in) :: Subscripts(size(lbounds))
    integer :: D
    integer :: I
    element_position = 0
    do i = size(lbounds), 1, -1
      d = ubounds(i) - lbounds(i) + 1
      if ( d <= 0 ) then
        element_position = 0
        return
      end if
      element_position = element_position * d + subscripts(i) - lbounds(i)
    end do
    element_position = element_position + 1
    
  end function Element_Position

  pure integer function Subscripts ( Lbounds, Ubounds, Element_Position )
    ! Compute subscripts from position in array element order.
    integer, intent(in), contiguous :: Lbounds(:)
    integer, intent(in) :: Ubounds(size(lbounds))
    integer, intent(in) :: Element_Position
    dimension :: Subscripts(size(lbounds))
    integer :: D, E, I
    e = element_position - 1
    do i = 1, size(lbounds)
      d = ubounds(i) - lbounds(i) + 1
      if ( d <= 0 ) then
        subscripts = lbounds - 1
        return
      end if
      subscripts(i) = mod(e,d) + lbounds(i)
      e = e / d
    end do
  end function Subscripts

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Array_Stuff

! $Log$
! Revision 2.1  2015/12/30 22:14:37  vsnyder
! Initial commit
!
