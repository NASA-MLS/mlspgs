! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Trace

  ! Print trace beginnings and ends

  implicit NONE
  private

  public :: Trace_Begin, Trace_Begin_Array, Trace_End

  integer, private :: Depth = -1
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Trace_Begin_Array ( Who, Number, Text, For )
    use Output_m, only: NewLine, Output
    character(len=*), intent(in) :: Who
    integer, intent(in), optional :: Number(:)
    character(len=*), intent(in), optional :: Text
    character(len=*), intent(in), optional :: For
    integer :: I
    depth = depth + 1
    call output ( repeat('.',max(depth,0)) )
    call output ( 'Enter ' )
    call output ( trim(who) )
    if ( present(for) ) call output ( for )
    if ( present(text) ) call output ( text )
    if ( present(number) ) then
      do i = 1, size(number)
        call output ( number(i), before=' ' )
      end do
    end if
    call newLine
  end subroutine Trace_Begin_Array

  subroutine Trace_Begin ( Who, Number, Text, For )
    use Output_m, only: NewLine, Output
    character(len=*), intent(in) :: Who
    integer, intent(in), optional :: Number
    character(len=*), intent(in), optional :: Text
    character(len=*), intent(in), optional :: For
    depth = depth + 1
    call output ( repeat('.',max(depth,0)) )
    call output ( 'Enter ' )
    call output ( trim(who) )
    if ( present(for) ) call output ( for )
    if ( present(text) ) call output ( text )
    if ( present(number) ) call output ( number )
    call newLine
  end subroutine Trace_Begin

  subroutine Trace_End ( Who, Number, Text )
    use Output_m, only: NewLine, Output
    character(len=*), intent(in) :: Who
    integer, intent(in), optional :: Number
    character(len=*), intent(in), optional :: Text
    call output ( repeat('.',max(depth,0)) )
    depth = depth - 1
    call output ( 'Exit  ' )
    call output ( trim(who) )
    if ( present(text) ) call output ( text )
    if ( present(number) ) call output ( number )
    call newLine
  end subroutine Trace_End

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Trace

! $Log$
