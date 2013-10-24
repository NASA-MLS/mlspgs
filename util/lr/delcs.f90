! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Delete

  implicit NONE
  private
  public :: DELCS

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine DELCS (IPTR)
    use IO, only: OUTPUT
    use LISTS, only: ITEM, NEXT, REL
    use S3, only: HEADCS
    use TOGGLES
    implicit NONE

  ! Delete the context set pointed to by IPTR.  The first cell of the
  ! list is a counter.  If the counter is greater than 1, simply
  ! decrement it.  Otherwise release the list.

    integer IPTR

  ! *****     Local Variables     ************************************

  ! I       Points to the node being examined.
  ! LAST    Points to the last node examined.
  ! LINE    Is used for message assembly.

    integer :: I, LAST
    character(len=120) :: LINE

  !     *****     Procedures     *****************************************

    if (item(iptr) > 1) then
      item(iptr) = item(iptr) - 1
    else
      if (toggle(ichar('3')) /= 0) then
        line=' Destroy context set at'
        write ( line(24:28), '(i5)' ) iptr
        call output (line(1:28))
      end if
      i = headcs(item(next(iptr)))
      last = 0
      do while (i /= 0)
        if (item(i) .eq. iptr) then
          if (last  /=  0) then
            next(last) = next(i)
          else
            headcs(item(next(iptr))) = next(i)
          end if
          next(i) = 0
          call rel (i)
          exit
        end if
        last = i
        i = next(i)
      end do
      call rel (iptr)
    end if

  end subroutine DELCS

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Delete

! $Log$
