! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module New_Context_Set

  implicit NONE
  private
  public :: NEWCS

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine NEWCS (IS, IPTR)

    use LISTS, only: LCOMPR, LIST, NEW, REL
    use Output_m, only: OUTPUT
    use Print_List, only: PNTLST
    use Tables, only: HEADCS
    use TOGGLES_LR

    implicit NONE

  ! Create a new context set having elements in the list pointed to
  ! by IS.  Return the pointer to the context set in IPTR.

    integer, intent(in) :: IS
    integer, intent(out) :: IPTR

  ! *****     External References     ********************************

  ! LCOMPR  compares two lists.
  ! PNTLST  Prints a context set.

  ! *****     Local Variables     ************************************

  ! I       is a list pointer subscript.
  ! IHEAD   is HEADCS(list(IS)%item).
  ! IIS     is list(IS)%item.

    integer I, IHEAD, IIS

  ! *****     Procedures     *****************************************

    iis = list(is)%item
    ihead = headcs(iis)
    i = ihead
    do while (i /= 0)
      iptr=list(i)%item
      if ( lcompr(list(iptr)%next, is) ) then
        ! Use the list at I, increase its reference count.
        list(iptr)%item = list(iptr)%item + 1
        if (toggle(iachar('3')) /= 0) then
          call output ( list(iptr)%next, before=' Use context set at ' )
          call output ( is, before=', destroy set at ', advance='yes' )
        end if
        call rel (is)
        return
      end if
      i = list(i)%next
    end do
    call new (iptr)
    list(iptr)%item = 1
    list(iptr)%next = is
    call new (i)
    list(i)%item = iptr
    list(i)%next = ihead
    headcs(iis) = i
    if (toggle(iachar('3')) /= 0) then
      call output ( is, before=' Create a context set at ' )
      call output ( ' = ' )
      call pntlst ( is, 33, 3 )
    end if

  end subroutine NEWCS

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module New_Context_Set

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
