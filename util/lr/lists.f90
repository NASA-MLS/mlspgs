! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module LISTS
! A package to manage lists.  ITEM is the array of list items.  NEXT is
! the array of pointers-to-next-items.
! ADDLTL adds a list to a list, eliminating duplicates
! LINT determines whether two lists intersect, that is, have at least
!      one common item.
! LISTS_INIT initializes the lists package.
! NEW provides a list node.
! REL releases a list or list node.

  implicit NONE
  public
  integer, parameter :: MAXLST = 5000
  integer, save :: ITEM(MAXLST), NEXT(MAXLST)

  ! GARBAG  points to the free list.  GARBAG IS A SAVE VARIABLE.

  integer, save, private :: GARBAG

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
! ===================================================     ADDLTL     =====
  subroutine ADDLTL (LPTR1, LPTR2, ICHNGE)

  ! Add the items in the list pointed to by LPTR1 to the list pointed to
  ! by LPTR2.  Both lists are initially in order.  Keep the list pointed
  ! to by LPTR2 in order.  Set ICHNGE if any items are added to the list
  ! pointed to by LPTR2.  If an item is added at the head of the list
  ! pointed to by LPTR2, LPTR2 will be changed to point to that item.

    integer, intent(in) :: LPTR1
    integer, intent(inout) :: LPTR2
    integer, intent(out) :: ICHNGE

    !     *****     Local Variables     **********************************

    ! I       is the pointer to an element of the list rooted at LPTR1.
    ! J       is the pointer to an element of the list rooted at LPTR2.
    ! LAST    is the previous list pointer during the search of the list
    !         rooted at LPTR2.
    ! NPTR    is the new list item fetched from NEW.

    integer I, J, LAST, NPTR

    ! *****     Procedures     *******************************************

    ichnge = 0
    i = lptr1
    j = lptr2
    last = 0

    ! Scan both lists as for a merge - they are always in order.

    do while (i  /=  0)
      if (j == 0) then
      ! Add remaining elements of list i to list j.
        ichnge = 1
        do while (i /= 0)
          call new (nptr)
          item(nptr) = item(i)
          if (last  /=  0) then
            next(last) = nptr
          else
            lptr2 = nptr
          end if
          last = nptr
          i = next(i)
        end do
        exit
      end if
      if (item(i) > item(j)) then
        last = j
        j = next(j)
      else if (item(i) == item(j)) then
        i = next(i)
        last = j
        j = next(j)
      else
        ichnge = 1
        call new (nptr)
        item(nptr) = item(i)
        next(nptr) = j
        if (last  /=  0) then
          next(last) = nptr
        else
          lptr2 = nptr
        end if
        last = nptr
        i = next(i)
      end if
    end do
    return
  end subroutine ADDLTL

! ====================================================     COPYL     =====

  subroutine COPYL (LPTR1, LPTR2)
  ! Create a new list rooted at LPTR2 that is a copy of the list
  ! rooted at LPTR1.

    integer, intent(in) :: LPTR1
    integer, intent(out) :: LPTR2

    ! *****     Local Variables     **************************************

    ! I       is the pointer to the node being copied.
    ! NPTR    is the pointer to the newly created node.

    integer I, NPTR

    !     *****     Procedures     ***************************************

    if (lptr1 == 0) then
      lptr2 = 0
    else
      ! Copy the first node.
      call new (lptr2)
      nptr = lptr2
      item(nptr) = item(lptr1)
      i = next(lptr1)
      ! Copy the rest of the list.
      do while (i > 0)
        call new (next(nptr))
        nptr = next(nptr)
        item(nptr) = item(i)
        i = next(i)
      end do
    end if

    return

  end subroutine COPYL

! ===================================================     LCOMPR     =====

  logical function LCOMPR (IPTR1, IPTR2)
    ! Compare the lists pointed to by IPTR1 and IPTR2.  Return .TRUE. if
    ! they are equal, and .FALSE. if they are unequal.

    integer, intent(in) :: IPTR1, IPTR2

    ! *****     Local Variables     **************************************

    ! I1      points to an element of IPTR1.
    ! I2      points to an element of IPTR2.

    integer I1, I2

    ! *****     Procedures     *******************************************

    lcompr = .true.
    if (iptr1 == iptr2) return
    lcompr = .false.
    i1 = iptr1
    i2 = iptr2
    do while (i1 /= 0)
      if (i2 == 0) return
      if (item(i1) /= item(i2)) return
      i1 = next(i1)
      i2 = next(i2)
    end do
    if (i2 .eq. 0) lcompr = .true.

    return

  end function LCOMPR

! =====================================================     LINT     =====

  logical function LINT (IPTR1, IPTR2)
  ! Returns .TRUE. if there is a common item in the lists pointed to by
  ! IPTR1 and IPTR2, .FALSE. otherwise.

    integer, intent(in) :: IPTR1, IPTR2

    ! *****     Local Variables     **************************************

    ! I1      points to an element of the list starting at IPTR1.
    ! I2      points to an element of the list starting at IPTR1.

    integer I1, I2

    ! *****     Procedures     *******************************************

    i1 = iptr1
    i2 = iptr2
    lint = .false.
    do
      if (i1 == 0 .or. i2 == 0) return
      lint = item(i1) == item(i2)
      if ( lint ) return
      if (item(i1) < item(i2)) then
        i1 = next(i1)
      else
        i2 = next(i2)
      end if
    end do
    return
  end function LINT

! ===============================================     LISTS_INIT     =====

  subroutine LISTS_INIT
    integer :: I
    do i = 1, maxlst-1
      next(i) = i + 1
      item(i) = 0
    end do
    next(maxlst) = 0
    item(maxlst) = 0
    garbag = 1
    return
  end subroutine LISTS_INIT

! ======================================================     NEW     =====

  subroutine NEW (IPTR)

  ! Fetch a new list node and return its subscript in IPTR.

    use Error_Handler, only: Error

    integer, intent(out) :: IPTR

    if (garbag == 0) call error ('LIST space overflow',2)
    if (item(garbag) .ne. 0) then
      iptr = item(garbag)
      item(garbag) = next(iptr)
    else
      iptr = garbag
      garbag = next(garbag)
    end if
    item(iptr) = 0
    next(iptr) = 0
    return
  end subroutine NEW

! ======================================================     REL     =====

  subroutine REL (IPTR)

  ! Release the list or list node at IPTR.

    integer, intent(in) :: IPTR

    if ( iptr > 0 ) then
      item(iptr) = next(iptr)
      next(iptr) = garbag
      garbag = iptr
    end if

  end subroutine REL

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module LISTS

! $Log$
