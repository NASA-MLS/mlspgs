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
    use Output_m, only: OUTPUT
    use LISTS, only: LIST, REL
    use Tables, only: HEADCS
    use TOGGLES_LR
    implicit NONE

  ! Delete the context set pointed to by IPTR.  The first cell of the
  ! list is a counter.  If the counter is greater than 1, simply
  ! decrement it.  Otherwise release the list.

    integer IPTR

  ! *****     Local Variables     ************************************

  ! I       Points to the node being examined.
  ! LAST    Points to the last node examined.

    integer :: I, LAST

  !     *****     Procedures     *****************************************

    if (list(iptr)%item > 1) then
      list(iptr)%item = list(iptr)%item - 1
    else
      if (toggle(iachar('3')) /= 0) then
        call output ( iptr, before=' Destroy context set at ', advance='yes')
      end if
      i = headcs(list(list(iptr)%next)%item)
      last = 0
      do while (i /= 0)
        if (list(i)%item == iptr) then
          if (last  /=  0) then
            list(last)%next = list(i)%next
          else
            headcs(list(list(iptr)%next)%item) = list(i)%next
          end if
          list(i)%next = 0
          call rel (i)
          exit
        end if
        last = i
        i = list(i)%next
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
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
