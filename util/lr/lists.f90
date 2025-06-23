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

  type :: Item_T
    integer :: Item = 0
    integer :: Next = 0
  end type Item_T

  type(item_t), allocatable, save :: List(:)

  integer, parameter, private :: Init_List = 1000

  ! GARBAG  points to the free list.  GARBAG IS A SAVE VARIABLE.

  integer, save, private :: GARBAG

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains
! ===================================================     ADDLTL     =====

  subroutine ADDLTL ( LPTR1, LPTR2, CHANGE )

  ! Add the items in the list pointed to by LPTR1 to the list pointed to
  ! by LPTR2.  Both lists are initially in order.  Keep the list pointed
  ! to by LPTR2 in order.  Set change if any items are added to the list
  ! pointed to by LPTR2.  If an item is added at the head of the list
  ! pointed to by LPTR2, LPTR2 will be changed to point to that item.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: LPTR1
    integer, intent(inout) :: LPTR2
    logical, intent(out) :: CHANGE

    !     *****     Local Variables     **********************************

    ! I       is the pointer to an element of the list rooted at LPTR1.
    ! J       is the pointer to an element of the list rooted at LPTR2.
    ! LAST    is the previous list pointer during the search of the list
    !         rooted at LPTR2.
    ! NPTR    is the new list item fetched from NEW.

    integer I, J, LAST, NPTR

    ! *****     Procedures     *******************************************

    if ( levels(gen) > 1 ) call trace_begin ( 'ADDLTL' )

    change = .false.
    i = lptr1
    j = lptr2
    last = 0

    ! Scan both lists as for a merge - they are always in order.

    do while ( i  /=  0 )
      if (j == 0) then
      ! Add remaining elements of list i to list j.
        change = .true.
        do while (i /= 0)
          call new ( nptr )
          list(nptr)%item = list(i)%item
          if ( last  /=  0 ) then
            list(last)%next = nptr
          else
            lptr2 = nptr
          end if
          last = nptr
          i = list(i)%next
        end do
        exit
      end if
      if ( list(i)%item > list(j)%item ) then
        last = j
        j = list(j)%next
      else if ( list(i)%item == list(j)%item ) then
        i = list(i)%next
        last = j
        j = list(j)%next
      else
        change = .true.
        call new ( nptr )
        list(nptr)%item = list(i)%item
        list(nptr)%next = j
        if ( last  /=  0 ) then
          list(last)%next = nptr
        else
          lptr2 = nptr
        end if
        last = nptr
        i = list(i)%next
      end if
    end do

    if ( levels(gen) > 1 ) call trace_end ( 'ADDLTL' )

  end subroutine ADDLTL

! ====================================================     COPYL     =====

  subroutine COPYL ( LPTR1, LPTR2 )
  ! Create a new list rooted at LPTR2 that is a copy of the list
  ! rooted at LPTR1.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: LPTR1
    integer, intent(out) :: LPTR2

    ! *****     Local Variables     **************************************

    ! I       is the pointer to the node being copied.
    ! NPTR    is the pointer to the newly created node.

    integer I, NPTR

    !     *****     Procedures     ***************************************

    if ( levels(gen) > 1 ) call trace_begin ( 'COPYL' )

    if ( lptr1 == 0 ) then
      lptr2 = 0
    else
      ! Copy the first node.
      call new ( lptr2 )
      nptr = lptr2
      list(nptr)%item = list(lptr1)%item
      i = list(lptr1)%next
      ! Copy the rest of the list.
      do while ( i > 0 )
        call new ( list(nptr)%next )
        nptr = list(nptr)%next
        list(nptr)%item = list(i)%item
        i = list(i)%next
      end do
    end if

    if ( levels(gen) > 1 ) call trace_end ( 'COPYL' )

  end subroutine COPYL

! ================================================     Dump_List     =====

  subroutine Dump_List ( IPTR )
    ! Dump the list starting at IPTR
    use Output_m, only: Output
    integer, intent(in) :: IPTR  ! Index of list to dump
    integer :: PTR               ! Index of list item
    ptr = iptr
    call output ( iptr, before='Dumping list starting at ', advance='yes' )
    do while ( ptr > 0 )
      call output ( ptr, 5 )
      call output ( list(ptr)%item, 5, before=': ' )
      call output ( list(ptr)%next, 5, advance='yes' )
      ptr = list(ptr)%next
    end do
    call output ( iptr, before='End of list starting at ', advance='yes' )
  end subroutine Dump_List

! ===================================================     LCOMPR     =====

  logical function LCOMPR ( IPTR1, IPTR2 )
    ! Compare the lists pointed to by IPTR1 and IPTR2.  Return .TRUE. if
    ! they are equal, and .FALSE. if they are unequal.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: IPTR1, IPTR2

    ! *****     Local Variables     **************************************

    ! I1      points to an element of IPTR1.
    ! I2      points to an element of IPTR2.

    integer I1, I2

    ! *****     Procedures     *******************************************

    if ( levels(gen) > 1 ) call trace_begin ( 'LCOMPR' )

    lcompr = .true.
    if ( iptr1 == iptr2 ) go to 9
    lcompr = .false.
    i1 = iptr1
    i2 = iptr2
    do while (i1 /= 0)
      if ( i2 == 0 ) go to 9
      if ( list(i1)%item /= list(i2)%item ) go to 9
      i1 = list(i1)%next
      i2 = list(i2)%next
    end do
    if ( i2 == 0 ) lcompr = .true.

  9 if ( levels(gen) > 1 ) call trace_end ( 'LCOMPR' )

  end function LCOMPR

! =====================================================     LINT     =====

  logical function LINT ( IPTR1, IPTR2 )
  ! Returns .TRUE. if there is a common item in the lists pointed to by
  ! IPTR1 and IPTR2, .FALSE. otherwise.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: IPTR1, IPTR2

    ! *****     Local Variables     **************************************

    ! I1      points to an element of the list starting at IPTR1.
    ! I2      points to an element of the list starting at IPTR1.

    integer I1, I2

    ! *****     Procedures     *******************************************

    if ( levels(gen) > 1 ) call trace_begin ( 'LINT' )

    i1 = iptr1
    i2 = iptr2
    lint = .false.
    do
      if ( i1 == 0 .or. i2 == 0 ) exit
      lint = list(i1)%item == list(i2)%item
      if ( lint ) exit
      if ( list(i1)%item < list(i2)%item ) then
        i1 = list(i1)%next
      else
        i2 = list(i2)%next
      end if
    end do

    if ( levels(gen) > 1 ) call trace_end ( 'LINT' )

  end function LINT

! ===============================================     LISTS_INIT     =====

  subroutine LISTS_INIT ( Start )

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in), optional :: Start ! Where to start chaining GARBAG
    integer :: I, MyStart, N

    if ( levels(gen) > 1 ) call trace_begin ( 'LISTS_INIT' )

    myStart = 1
    if ( present(start) ) myStart = start
    if ( .not. allocated(list) ) allocate ( list(init_list) )
    n = size(list)
    do i = myStart, n-1
      list(i)%next = i + 1
      list(i)%item = 0
    end do
    list(n)%next = 0
    list(n)%item = 0
    garbag = myStart

    if ( levels(gen) > 1 ) call trace_end ( 'LISTS_INIT' )

  end subroutine LISTS_INIT

! ======================================================     NEW     =====

  subroutine NEW (IPTR)

  ! Fetch a new list node and return its subscript in IPTR.

    use Output_m, only: Output
    use Toggles, only: Gen, Levels
    use Toggles_LR, only: Toggle_LR => Toggle
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(out) :: IPTR

    integer :: N
    type(item_t), allocatable :: Temp(:)

    if ( levels(gen) > 1 ) call trace_begin ( 'NEW' )

    if (garbag == 0) then
      if ( allocated(list) ) then
        n = size(list)
        allocate ( temp(2*n) )
        temp(1:n) = list
        call move_alloc ( temp, list )
        call lists_init ( n+1 )
      else
        call lists_init
      end if
    end if
    if ( list(garbag)%item /= 0) then
      iptr = list(garbag)%item
      list(garbag)%item = list(iptr)%next
    else
      iptr = garbag
      garbag = list(garbag)%next
    end if
    list(iptr)%item = 0
    list(iptr)%next = 0
    if ( toggle_LR(iachar('4')) /= 0 ) &
      & call output ( iptr, before=' New list at ', advance='yes' )

    if ( levels(gen) > 1 ) call trace_end ( 'NEW' )

  end subroutine NEW

! ======================================================     REL     =====

  subroutine REL (IPTR)

  ! Release the list or list node at IPTR.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: IPTR

    if ( levels(gen) > 1 ) call trace_begin ( 'REL' )

    if ( iptr > 0 ) then
      list(iptr)%item = list(iptr)%next
      list(iptr)%next = garbag
      garbag = iptr
    end if

    if ( levels(gen) > 1 ) call trace_end ( 'REL' )

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
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
