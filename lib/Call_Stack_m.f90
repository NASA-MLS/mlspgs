! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Call_Stack_m

  implicit NONE
  private

  public :: Stack_t
  public :: Dump_Stack, Get_Frame, Pop_Stack, Push_Stack, Stack_Depth, Top_Stack

  type :: Stack_t
    real :: Clock = 0.0                ! Whatever Time_Now returns (CPU or ?)
    integer :: Index = -1              ! Extra info, e.g. MAF number
    double precision :: Memory = 0.0d0 ! In use, in Memory_Units
    character(len=10) :: Now = ''      ! Date and time
    integer :: Text = 0                ! Index in string table
    integer :: Tree = 0                ! Where in l2cf, -1 if stack not allocated,
                                       ! -2 if stack index < 1, -3 if stack index
                                       ! > stack_ptr.
  end type

  interface Push_Stack
    module procedure Push_Stack_b, Push_Stack_c, Push_Stack_i
  end interface

  integer, parameter, private :: StartingStackSize = 100
  type(stack_t), allocatable, save, private :: Stack(:)
  integer, save, private :: Stack_Ptr = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

! ---------------------------------------------------  Dump_Stack  -----
  subroutine Dump_Stack ( Top, Before, Where, Size, CPU, DoDepth, Rev, Advance )
    ! Dump the call stack.  If Top is present and true, dump only the top
    ! frame.  If Before is present, print it after "depth" dots and before
    ! anything else.  If Where is present and true, and the Tree component
    ! is not zero, print the line and column number from the configuration
    ! file.

    use Allocate_Deallocate, only: Memory_Units
    use Lexer_Core, only: Print_Source
    use Output_m, only: DumpSize, Output
    use String_Table, only: Display_String
    use Tree, only: Dump_Tree_Node, Source_Ref

    logical, intent(in), optional :: Top   ! Dump only the top frame
    character(len=*), intent(in), optional :: Before ! first thing output
    logical, intent(in), optional :: Where ! Dump tree location
    logical, intent(in), optional :: Size  ! Dump memory size (default true)
    logical, intent(in), optional :: CPU   ! Print CPU (default false)
    logical, intent(in), optional :: DoDepth ! Print "depth" dots (default true)
    logical, intent(in), optional :: Rev   ! Print in reverse order (default false)
    character(len=*), intent(in), optional :: Advance ! Default 'yes'

    character(len=3) :: MyAdvance
    integer :: Depth, First, I, Inc, Last
    logical :: Error
    logical :: MyCPU, MyDoDepth, MyRev, MySize, MyTop, MyWhere

    myAdvance = 'yes'
    myCPU = .false.; myDoDepth = .true.; myRev = .false.; mySize = .true.
    myTop = .false.; myWhere = .false.

    if ( present(top) ) myTop = top
    if ( present(advance) ) myAdvance = advance
    if ( stack_ptr /= merge(stack_ptr,lbound(stack,1),myTop) ) myAdvance = 'yes'
    if ( present(cpu) ) myCPU = cpu
    if ( present(doDepth) ) myDoDepth = doDepth
    if ( present(rev) ) myRev = rev
    if ( present(size) ) mySize = size
    if ( present(where) ) myWhere = Where

    error = .not. allocated(stack)
    if ( .not. error ) error = stack_ptr < lbound(stack,1)
    if ( error ) then
      call output ( "*****" )
      if ( present(before) ) call output ( " " // before )
      call output ( " There is no stack to dump *****", advance='yes' )
      return
    end if

    if ( myRev ) then ! Dump stack bottop-up
      first = merge(stack_ptr,lbound(stack,1),myTop)
      last = stack_ptr
      inc = 1
    else              ! Dump stack top-down
      first = stack_ptr
      last = merge(stack_ptr,lbound(stack,1),myTop)
      inc = -1
    end if

    do depth = first, last, inc
      do i = lbound(stack,1), merge(depth,-1,myDoDepth)
        call output ( '.' )
      end do
      if ( present(before) ) call output ( before )
      if ( stack(depth)%text /= 0 ) call display_string ( stack(depth)%text )
      if ( stack(depth)%index >= 0 ) call output ( stack(depth)%index )
      if ( myWhere .and. stack(depth)%tree /= 0 ) &
        & call print_source ( source_ref(stack(depth)%tree), before=', ' )
      call output ( ' at ' // stack(depth)%now(1:2) // ':' // &
        & stack(depth)%now(3:4) // ':' // trim(stack(depth)%now(5:)) )
      if ( mySize ) call dumpSize ( memory_units * stack(depth)%memory, &
        & before = ' Memory: ' )
      if ( myCPU ) call output ( stack(depth)%clock, format='(g10.3)', &
        & before=' CPU: ' )
      call output ( '', advance=myAdvance )
    end do

  end subroutine Dump_Stack

! ----------------------------------------------------  Get_Frame  -----
  type(stack_t) function Get_Frame ( Depth )
    ! Get stack(depth) if 1 <= Depth <= Stack_Ptr, otherwise result is
    ! default initialized except if Stack is not allocated, %tree= -1, if
    ! Depth < lbound(stack,1), %tree= -2, and if Depth > Stack_Ptr, %tree= -3.
    integer, intent(in) :: Depth

    if ( .not. allocated(stack) ) then
      get_frame%tree = -1
    else if ( depth < lbound(stack,1) ) then
      get_frame%tree = -2
    else if ( depth > stack_ptr ) then
      get_frame%tree = -3
    else
      get_frame = stack(depth)
    end if

  end function Get_Frame

! ----------------------------------------------------  Pop_Stack  -----
  subroutine Pop_Stack ( Before, Where )
    ! Pop the stack.  If Before or Where are present, dump the top frame first.

    use Allocate_Deallocate, only: Memory_Units, NoBytesAllocated
    use Output_m, only: DumpSize, NewLine, Output
    use Time_m, only: Time_Now

    character(len=*), intent(in), optional :: Before
    logical, intent(in), optional :: Where

    double precision :: Delta
    real :: T
    character(len=10) :: Used

    if ( present(before) .or. present(where) ) then
      call dump_stack ( .true., before, where, size=.false., advance='no' )
      if ( .not. allocated(stack) ) return
      if ( stack_ptr < lbound(stack,1) ) return
      call time_now ( t )
      t = t - stack(stack_ptr)%clock
      write ( used, '(g10.3)' ) t
      call output ( ' used ' // trim(adjustl(used)) //  ' cpu' )
      if ( stack(stack_ptr)%memory /= noBytesAllocated ) then
        delta = memory_units * (noBytesAllocated - stack(stack_ptr)%memory)
        if ( abs(delta) < huge(1) ) then
          call dumpSize ( int(delta), before=', Memory changed by ' )
        else
          call dumpSize ( delta, before=', Memory changed by ' )
        end if
        call dumpSize ( memory_units * noBytesAllocated, before = ' to ' )
      end if
      call newLine
    end if

    stack_ptr = stack_ptr - 1

  end subroutine Pop_Stack

! -------------------------------------------------  Push_Stack_B  -----
  subroutine Push_Stack_B ( Name_I, Name_C, Root, Index, Before, Where )
    ! If Name_I <= 0, use Create_String ( Name_C ) to give it a value.
    ! We assume the actual argument is a SAVE variable.
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use String_Table, only: Create_String

    integer, intent(inout) :: Name_I
    character(len=*), intent(in) :: Name_C
    integer, optional, intent(in) :: Root  ! Where in configuration tree
    integer, optional, intent(in) :: Index ! Whatever caller wants to send
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where

    if ( name_i <= 0 ) name_i = create_string ( name_c )
    call push_stack ( name_i, Root, Index, Before, Where )

  end subroutine Push_Stack_B

! -------------------------------------------------  Push_Stack_C  -----
  subroutine Push_Stack_C ( Name, Root, Index, Before, Where )
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use String_Table, only: Create_String

    character(len=*), intent(in) :: Name
    integer, optional, intent(in) :: Root  ! Where in configuration tree
    integer, optional, intent(in) :: Index ! Whatever caller wants to send
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where

    call push_stack ( create_string ( name ), Root, Index, Before, Where )

  end subroutine Push_Stack_C

! -------------------------------------------------  Push_Stack_I  -----
  subroutine Push_Stack_I ( Name, Root, Index, Before, Where )
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use Allocate_Deallocate, only: Test_Allocate, NoBytesAllocated
    use Time_m, only: Time_Now

    integer, intent(in) :: Name
    integer, optional, intent(in) :: Root  ! Where in configuration tree
    integer, optional, intent(in) :: Index ! Whatever caller wants to send
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where

    integer :: Stat
    type(stack_t), allocatable :: Temp_Stack(:)
    intrinsic :: Date_And_Time, Storage_Size

    if ( .not. allocated(stack) ) then
      ! If you allocate with lbound < 0, other stuff won't work.
      allocate ( stack(startingStackSize), stat=stat )
      call test_allocate ( stat, moduleName, 'Stack', &
        & ubounds=(/startingStackSize/), elementSize=storage_size(stack) )
      stack_ptr = lbound(stack,1) - 1
    end if
    if ( stack_ptr >= ubound(stack,1) ) then
      allocate ( temp_stack(2*stack_ptr), stat=stat )
      call test_allocate ( stat, moduleName, 'Temp_Stack', &
        & ubounds=(/2*stack_ptr/), elementSize=storage_size(stack) )
      temp_stack(:min(stack_ptr,ubound(stack,1))) = stack
      call move_alloc ( temp_stack, stack )
    end if
    stack_ptr = stack_ptr + 1
    stack(stack_ptr) = stack_t ( memory=noBytesAllocated, text=name )
    if ( present(root) ) stack(stack_ptr)%tree = root
    if ( present(index) ) stack(stack_ptr)%index = index
    call time_now ( stack(stack_ptr)%clock )
    call date_and_time ( stack(stack_ptr)%now )
    if ( present(before) .or. present(where) ) &
      & call dump_stack ( .true., before, where, advance='yes' )
  end subroutine Push_Stack_I

! --------------------------------------------------  Stack_Depth  -----
  integer function Stack_Depth ()
    ! Return -1 if the stack is not allocated, else Stack_Ptr
    stack_depth = merge(stack_ptr, -1, allocated(stack) )
  end function Stack_Depth

! ----------------------------------------------------  Top_Stack  -----
  subroutine Top_Stack ( Top )
    type(stack_t), intent(out) :: Top
    if ( .not. allocated(stack) ) then
      top%tree = -1
    else if ( stack_ptr < lbound(stack,1) ) then
      top%tree = -2
    else
      top = stack(stack_ptr)
    end if
  end subroutine Top_Stack

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Call_Stack_m

! $Log$
! Revision 2.5  2013/08/29 19:47:52  vsnyder
! Don't output anything in Pop_Stack if Before and Where absent
!
! Revision 2.4  2013/08/23 03:16:52  vsnyder
! Simplify Push_Stack, some cannonball polishing (I couldn't resist)
!
! Revision 2.3  2013/08/23 02:49:19  vsnyder
! Add Get_Frame, Rev, and Stack_Depth
!
! Revision 2.2  2013/08/17 03:08:12  vsnyder
! Some cannonball polishing (already!)
!
! Revision 2.1  2013/08/17 02:57:41  vsnyder
! Initial commit
!
