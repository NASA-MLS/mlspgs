! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Print_List

  implicit NONE
  private
  public :: PNTLST

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine PNTLST (IPTR, LSTART, LCONT)
    use LISTS, only: LIST
    use Output_m, only: Blanks, NewLine
    use String_Table, only: Display_String, String_Length
    use Tables, only: VOCAB
    implicit NONE

  ! Print the context set pointed to by IPTR.
  ! Start printing the list at position LSTART in LINE.
  ! If LINE becomes filled, print it and continue printing in LCONT.

    integer, intent(in) :: IPTR, LSTART, LCONT

  ! *****     Local Variables     ************************************

  ! IP      is a pointer that steps through the list rooted at IPTR.
  ! L       is a subscript for LINE.

    integer IP, L

  ! *****     Procedures     *****************************************

    ip = iptr
    l = lstart
    do while (ip /= 0)
      if (string_length(vocab(list(ip)%item))+l >= 114) then
        call newLine
        call blanks ( lcont )
        l = lcont
      end if
      call display_string (vocab(list(ip)%item) , before=' ')
      l = l + string_length(vocab(list(ip)%item)) + 1
      ip = list(ip)%next
    end do
    call newLine

    return

  end subroutine PNTLST

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Print_List

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
