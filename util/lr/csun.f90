! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Union

  implicit NONE
  private
  public ::CSUN

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine CSUN ( IPTR1, IPTR2, CHANGE )

    use Delete, only: DELCS
    use LISTS, only: ADDLTL, COPYL, LIST
    use New_Context_Set, only: NEWCS

    implicit NONE

  ! Union the contexts sets pointed to by IPTR1 and IPTR2, leaving the
  ! result in IPTR2.   CHANGE is returned true iff there are any
  ! elements of IPTR1 that are not members of IPTR2.

    integer, intent(in) :: IPTR1
    integer, intent(inout) :: IPTR2
    logical, intent(out) :: CHANGE

  ! *****     External References     ********************************

  ! ADDLTL  adds the elements of one list to another.
  ! COPYL   copies a list.
  ! DELCS   deletes a context set.
  ! NEWCS   creates a new context set.

  ! *****     Local Variables     ************************************

  ! HEAD    is the head of a the list that represents the union of IPTR1
  !         and IPTR2.

    integer head

  !     *****     Procedures     *****************************************

    change = .false.
    if (iptr1 /= iptr2) then
      call copyl ( list(iptr2)%next, head )
      call addltl ( list(iptr1)%next, head, change )
      call delcs ( iptr2 )
      call newcs ( head, iptr2 )
    end if

    return

  end subroutine CSUN

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Union

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
