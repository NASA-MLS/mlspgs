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

    use IO, only: OUTPUT
    use LISTS, only: ITEM, LCOMPR, NEW, NEXT, REL
    use Print_List, only: PNTLST
    use S3, only: HEADCS
    use TOGGLES

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
  ! IHEAD   is HEADCS(ITEM(IS)).
  ! IIS     is ITEM(IS).
  ! LINE    Is used for message assembly.

    integer I, IHEAD, IIS
    character(len=120) :: LINE

  ! *****     Procedures     *****************************************

    iis = item(is)
    ihead = headcs(iis)
    i = ihead
    do while (i /= 0)
      iptr=item(i)
      if ( lcompr(next(iptr), is) ) then
        ! Use the list at I, increase its reference count.
        item(iptr) = item(iptr) + 1
        if (toggle(ichar('3')) /= 0) then
          write ( line(1:44), '(a,i5,a,i4)' ) &
            ' Use context set at', next(iptr), ', destroy set at', is
          call output (line(1:44))
        end if
        call rel (is)
        return
      end if
      i = next(i)
    end do
    call new (iptr)
    item(iptr) = 1
    next(iptr) = is
    call new (i)
    item(i) = iptr
    next(i) = ihead
    headcs(iis) = i
    if (toggle(ichar('3')) /= 0) then
      line=' Create a context set at      ='
      write ( line(25:29), '(i5)' ) is
      call pntlst (is,line,33,3)
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
