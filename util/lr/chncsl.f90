! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Chain_Context_Lists

  implicit NONE
  private
  public :: CHNCSL
  integer,public :: LENCSL, LISTCS, LSETS, NCSETS = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine CHNCSL

    use Error_m, only: Error
    use LISTS, only: LIST, NEW, REL
    use Tables, only: HEADCS, NTERMS

    implicit NONE

  ! Chain the context sets together.

  ! *****     External References     ********************************

  ! ERROR   prints error messages.
  ! NEW     fetches a new list node.
  ! REL     releases a list node.

  ! *****     Local Variables     ************************************

  ! I       is a loop inductor and subscript for HEADCS.
  ! IPTR    is the root of the list taken from HEADCS(I).
  ! LAST    is the previous value of IPTR during the list search.
  ! LPTR    points to the last list node fetched.
  ! NPTR    is the pointer to the list item fetched from NEW.

    integer I, IPTR, LAST, LPTR, NPTR

  !     *****     Procedures     *****************************************

  ! Find the first context set chain.

    lsets = 0
    ncsets = 0

    do i = 1, nterms
      if ( headcs(i) /= 0 ) go to 10
    end do
    call error ('CHNCSL failed to find a context set chain',2)
  10 continue

  ! Start the context set length list pointed to by LENCSL.
  ! Save the head of the context set chain and start.

    call new (nptr)
    lencsl = nptr
    listcs = headcs(i)
  o:do
      iptr = headcs(i)
      do while (iptr /= 0)
        ncsets = ncsets + 1
        list(list(iptr)%item)%item = ncsets
        ! Store the list length in list(nptr)%item.
        list(nptr)%item = 0
        lptr = list(list(iptr)%item)%next
        do while (lptr > 0)
          list(nptr)%item = list(nptr)%item + 1
          lptr = list(lptr)%next
        end do
        lsets = lsets + list(nptr)%item
        lptr = nptr
        call new ( list(nptr)%next )
        nptr = list(nptr)%next
        last = iptr
        iptr = list(iptr)%next
      end do

      ! Search for the next context set chain.  Link it into the previous
      ! chain and start on the new chain.

      do while (i < nterms)
        i = i + 1
        if (headcs(i) /= 0) then
          list(last)%next = headcs(i)
          cycle o
        end if
      end do
      exit o
    end do o
    list(lptr)%next = 0
    call rel (nptr)

  end subroutine CHNCSL

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Chain_Context_Lists

! $Log$
! Revision 1.1  2013/10/24 22:41:13  vsnyder
! Initial commit
!
