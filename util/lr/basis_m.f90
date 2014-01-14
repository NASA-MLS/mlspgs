! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Basis_m

  ! Data structures and procedures related to the basis of a configuration
  ! set, sometimes called a state.

  implicit NONE

  private

  ! A BASIS record contains information about the configuration:
  type, public :: Basis_T
    integer :: Item      ! Index of first item in ITEMS for the configuration
                         ! set
    integer :: Next      ! Index is Basis of the next configuration set to be
                         ! processed, or -1 if the configuration set is not
                         ! on the queue.  It is not added to the queue again
                         ! if it is already on the queue.
    integer :: Same      ! Index is Basis of the next configuration having the
                         ! same entrance symbol.
    integer :: Tran      ! Index in TRAN of first transition for the state.
                         ! Transitions are stored in TRAN as positions in the
                         ! BASIS array.
    integer :: Red       ! Index in RED of first reduction from the state.
                         ! Reductions are stored in RED as a pair consisting
                         ! of the number of the production to be reduced and
                         ! a pointer to the context set.
  end type Basis_T

  type(basis_t), public, save, allocatable :: Basis(:)
  integer, parameter :: Init_Basis = 4000
  integer, save :: Basis_Size = 0

  ! Let P be the pointer to the first item = BASIS(INDEX)%Item.  Then
  ! for the I'th item in the configuration set, Items(P) contains:
  type, public :: Item_T
    integer :: Prod      ! Production number
    integer :: Dot       ! Position of the dot in Prod
    integer :: Set       ! Index of context set list
  end type Item_T

  type(item_t), public, save, allocatable :: Items(:)
  integer, parameter :: Init_Items = 1000
  integer, save :: Items_Size = 0 ! Items starts out deallocated

  integer, public :: IFINAL, INDBAS = 1

  public :: ADDBAS, DEQUE, Dump_Basis, Dump_Items, ENQUE, Increase_Items, NEWBAS

  ! QHEAD   points to the head of the queue.
  ! QTAIL   points to the tail of the queue.
  !         QHEAD AND QTAIL ARE SAVE VARIABLES.

  integer, private, save :: QHEAD = 0, QTAIL = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! ===================================================     ADDBAS     =====
  subroutine ADDBAS (IPTR, NPR, NDOT, NSET)

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    ! Add the item consisting of production NPR, having the
    ! dot before the NDOT'th item in its right hand side and context
    ! NSET, to the basis at IPTR.

    ! The items of the basis of the configuration set are stored in
    ! reverse order starting at the end of the BASIS array.

    integer, intent(in) :: IPTR, NPR, NDOT, NSET

    ! *****     Local Variables     ************************************

    ! INDITM  is the index for the next available space to store an
    !         item.  INDITM is available from BASIS(IPTR+5).

    integer INDITM

    ! *****     Procedures     *****************************************

    if ( levels(gen) > 1 ) call trace_begin ( 'ADDBAS' )

    inditm = basis(iptr+1)%item
    ! Leave room for one more basis state in the bottom of the BASIS
    ! array, because the extent of the items in the BASIS array,
    ! the extent of transitions in the TRAN array, and the extent of
    ! reductions in the RED array are indicated by the difference
    ! between the starting pointers in two consecutive state records.

    if ( inditm > items_size ) &
      & call increase_items ( items, items_size, for=' for Basis' )
    items(inditm) = item_t(npr,ndot,nset)
    basis(iptr+1)%item = inditm + 1

    if ( levels(gen) > 1 ) call trace_end ( 'ADDBAS' )

  end subroutine ADDBAS

! ====================================================     DEQUE     =====
  subroutine DEQUE ( IPTR )
    !   Removes the basis at the head of the queue and returns a pointer
    !   to it in IPTR.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(out) :: IPTR

    if ( levels(gen) > 1 ) call trace_begin ( 'DEQUE' )

    iptr = qhead
    if (qhead  /=  0) then
       qhead = basis(qhead)%next
       basis(iptr)%next = -1
    end if

    if ( levels(gen) > 1 ) call trace_end ( 'DEQUE' )

  end subroutine DEQUE

! ===============================================     Dump_Basis     =====
  subroutine Dump_Basis ( First, Last, Before )
    use Output_m, only: Output
    integer, intent(in) :: First, Last ! Which ones to dump
    character(*), intent(in), optional :: Before
    integer :: I
    if ( present(before) ) call output ( before, advance='yes' )
    call output ( '       items next same tran  red', advance='yes' )
    do i = first, last
      call output ( i, 5 )
      call output ( basis(i)%item, 5, before=': ' )
      call output ( basis(i)%next, 5 )
      call output ( basis(i)%same, 5 )
      call output ( basis(i)%tran, 5 )
      call output ( basis(i)%red, 5, advance='yes' )
    end do
  end subroutine Dump_Basis

! ===============================================     Dump_Items     =====
  subroutine Dump_Items ( Items, First, Last, Before )

    use Output_m, only: Output

    type(item_t), intent(in) :: Items(:)
    integer, intent(in) :: First, Last
    character(*), intent(in), optional :: Before
    integer :: I

    if ( present(before) ) call output ( before, advance='yes' )
    call output ( '        Prod Dot  Set', advance='yes' )
    do i = first, last
      call output ( i, 5 )
      call output ( items(i)%prod, 5, before=': ' )
      call output ( items(i)%dot,  5 )
      call output ( items(i)%set,  5, advance='yes' )
    end do

  end subroutine Dump_Items

! ====================================================     ENQUE     =====
  subroutine ENQUE ( IPTR )
    !   Enters the basis pointed at by IPTR at the tail of the queue of
    !   BASIS objects to be processed.

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer, intent(in) :: IPTR

    if ( levels(gen) > 1 ) call trace_begin ( 'ENQUE' )

    if ( basis(iptr)%next  ==  -1) then
      basis(iptr)%next = 0
      if (qhead  ==  0) then
        qhead = iptr
      else
        basis(qtail)%next = iptr
        end if
        qtail = iptr
     end if

    if ( levels(gen) > 1 ) call trace_end ( 'ENQUE' )

  end subroutine ENQUE

! ===========================================     Increase_Basis     =====
  subroutine Increase_Basis

    ! If BASIS is not allocated, allocate it with Init_Basis elements.
    ! Otherwise, double it.

    use Error_m, only: Error
    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    integer :: N                            ! Current size of BASIS
    integer :: Stat                         ! From ALLOCATE
    type(basis_t), allocatable :: Temp(:)   ! Used for doubling BASIS

    if ( levels(gen) > 1 ) &
      & call trace_begin ( 'Increase_Basis', number=basis_size, text=' from ' )

    if ( .not. allocated(basis) ) then
      n = 0
      allocate ( basis(init_basis), stat=stat )
      basis(1) = basis_t(item=1,same=0,next=0,tran=1,red=1)
    else
      n = size(basis)
      allocate ( temp(2*n), stat=stat )
    end if
    if ( stat /= 0 ) &
      & call error ( 'Unable to allocate new BASIS space', 9 ) ! stops
    if ( n > 0 ) then
      temp(:n) = basis(:n)
      call move_alloc ( temp, basis )
    end if
    basis_size = size(basis)

    if ( levels(gen) > 1 ) &
      & call trace_end ( 'Increase_Basis', number=basis_size, text=' to ' )

  end subroutine Increase_Basis


! ===========================================     Increase_Items     =====
  subroutine Increase_Items ( Items, NewSize, For )

    ! If ITEMS is not allocated, allocate it with Init_Items elements.
    ! Otherwise, double it.

    use Error_m, only: Error
    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    type(item_t), allocatable, intent(inout) :: Items(:)
    integer, intent(out) :: NewSize
    character(len=*), optional :: For

    integer :: N                            ! Current size of Items
    integer :: Stat                         ! From ALLOCATE
    type(item_t), allocatable :: Temp(:)    ! Used for doubling Items

    n = 0
    if ( allocated(items) ) n = size(items)

    if ( levels(gen) > 1 ) &
      & call trace_begin ( 'Increase_Items', number=n, text=' from ', for=for )

    if ( .not. allocated(items) ) then
      allocate ( items(init_items), stat=stat )
    else
      allocate ( temp(2*n), stat=stat )
    end if
    if ( stat /= 0 ) &
      & call error ( 'Unable to allocate new scratch space', 9 ) ! stops
    if ( n > 0 ) then
      temp(:n) = items(:n)
      call move_alloc ( temp, items )
    end if
    newSize = size(items)

    if ( levels(gen) > 1 ) &
      & call trace_end ( 'Increase_Items', number=newSize, text=' to ' )

  end subroutine Increase_Items

! ===================================================     NEWBAS     =====
  subroutine NEWBAS ( INDEX )

    use Toggles, only: Gen, Levels
    use Trace, only: Trace_Begin, Trace_End

    implicit NONE

    ! Prepare space for a new basis.  Return the pointer to the space in
    ! INDEX.

    integer, intent(out) :: INDEX

    ! *****     Procedures     *****************************************

    ! Make sure there is room for the current configuration set, the
    ! next configuration set (which if this is the last is used to
    ! calculate the number of items, transitions and reductions),
    ! and one item.

    if ( levels(gen) > 1 ) call trace_begin ( 'NEWBAS' )

    if ( indbas+2 > basis_size ) call increase_basis
    basis(indbas+1) = basis_t ( basis(indbas)%item, 0, 0, huge(0), huge(0) )
    ! Store pointer to next basis to be processed, -1 means this basis
    ! is not on a queue.
    basis(indbas)%next = -1
    ! Store pointer to next basis having same entry symbol.
    basis(indbas)%same = 0
    index = indbas
    indbas = indbas + 1

    if ( levels(gen) > 1 ) call trace_end ( 'NEWBAS' )

  end subroutine NEWBAS

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Basis_m

! $Log$
! Revision 1.1  2013/10/24 22:41:14  vsnyder
! Initial commit
!
