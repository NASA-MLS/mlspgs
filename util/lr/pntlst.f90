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

  subroutine PNTLST (IPTR, LINE, LSTART, LCONT)
    use IO, only: OUTPUT
    use LISTS, only: ITEM, NEXT
    use S1, only: LENGTH, MOVSTR
    use S3, only: VOCAB
    implicit NONE

  ! Print the context set pointed to by IPTR.
  ! Start printing the list at position LSTART in LINE.
  ! If LINE becomes filled, print it and continue printing in LCONT.

    integer IPTR, LSTART, LCONT
    character(len=*) :: LINE

  ! *****     External References     ********************************

  ! LENGTH  calculates the length of a vocabulary item.
  ! MOVSTR  moves a vocabulary item from the symbol table.

  ! *****     Local Variables     ************************************

  ! IP      is a pointer that steps through the list rooted at IPTR.
  ! L       is a subscript for LINE.

    integer IP, L

  ! *****     Procedures     *****************************************

    ip = iptr
    l = lstart
    do while (ip /= 0)
      if (length(vocab(item(ip)))+l >= 114) then
        call output (line(1:l-1))
        l = lcont
      end if
      call movstr (vocab(item(ip)), line, l, 120)
      l = l + 1
      ip = next(ip)
    end do
    call output (line(1:l-1))

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
