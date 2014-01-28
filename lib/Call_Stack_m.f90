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

  public :: STACK_T
  public :: DUMP_STACK, GET_FRAME, POP_STACK, PUSH_STACK, STACK_DEPTH, TOP_STACK
  public :: SAY_WHEN

  type :: Stack_t
    real :: Clock = 0.0                ! Whatever Time_Now returns (CPU or ?)
    integer :: Index = -1              ! Extra info, e.g. MAF number
    double precision :: Memory = 0.0d0 ! In use, in Memory_Units
    integer :: String = 0              ! Index in string table
    integer :: Text = 0                ! Index in string table
    integer :: Tree = 0                ! Where in l2cf, -1 if stack not allocated,
                                       ! -2 if stack index < 1, -3 if stack index
                                       ! > stack_ptr.
  end type

  interface Push_Stack
    module procedure Push_Stack_b, Push_Stack_c, Push_Stack_i
  end interface

  ! Display date and time during push and pop displays
  logical, save :: Say_When = .false.

  logical, save, private :: StaySilent = .false. ! E.g., in case it becomes corrupted
  logical, save, private :: Verbose = .false. ! Print each push, pop
  integer, parameter, private :: StartingStackSize = 100
  integer, parameter, private :: MAXDOUBLINGS = 10 ! To limit max stack size
  type(stack_t), allocatable, save, private :: Stack(:)
  integer, save, private :: Stack_Ptr = 0
  integer, save, private :: Stack_Doublings = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! ====     Public Procedures     ==============================

! ---------------------------------------------------  Dump_Stack  -----
  subroutine Dump_Stack ( Top, Before, Where, Size, CPU, DoDepth, Rev, Index, &
                        & String, StringIndex, Advance )
    ! Dump the call stack.  If Top is present and true, dump only the top
    ! frame.  If Before is present, print it after "depth" dots and before
    ! anything else.  If Where is present and true, and the Tree component
    ! is not zero, print the line and column number from the configuration
    ! file.

    use ALLOCATE_DEALLOCATE, only: MEMORY_UNITS
    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use HIGHOUTPUT, only: DUMPSIZE
    use STRING_TABLE, only: DISPLAY_STRING
    use TREE, only: WHERE_AT=>Where

    logical, intent(in), optional :: Top   ! Dump only the top frame
    character(len=*), intent(in), optional :: Before ! first thing output
    logical, intent(in), optional :: Where ! Dump tree location
    logical, intent(in), optional :: Size  ! Dump memory size (default true)
    logical, intent(in), optional :: CPU   ! Print CPU (default false)
    logical, intent(in), optional :: DoDepth ! Print "depth" dots (default true)
    logical, intent(in), optional :: Rev   ! Print in reverse order (default false)
    integer, intent(in), optional :: Index ! Print this instead of from stack
    character(len=*), optional, intent(in) :: String
    integer, optional, intent(in) :: StringIndex
    character(len=*), intent(in), optional :: Advance ! Default 'yes'

    character(len=3) :: MyAdvance
    integer :: Depth, First, I, Inc, Last
    logical :: Error
    logical :: MyCPU, MyDoDepth, MyRev, MySize, MyTop, MyWhere
    ! Executable
    if ( StaySilent ) return
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
      if ( stack(depth)%text > 0 ) call display_string ( stack(depth)%text )
      if ( stack(depth)%string > 0 ) &
        & call display_string ( stack(depth)%string, before=' ' )
      if ( present(index) ) then
        call output ( index, before=' ' )
      else if ( stack(depth)%index >= 0 ) then
        call output ( stack(depth)%index, before=' ' )
      end if
      if ( present(string) ) call output ( ' ' // string )
      if ( present(stringIndex) ) then
        if ( stringIndex > 0 ) call display_string ( stringIndex, before=' ' )
      end if
      if ( myWhere .and. stack(depth)%tree > 0 ) then
        call output ( stack(depth)%tree, before=', tree at ' )
        call print_source ( where_at(stack(depth)%tree), before=', ' )
      end if
      if ( say_when ) call show_when
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
  subroutine Pop_Stack ( Before, Where, Frame, Index, String, StringIndex, &
    & Silent )
    ! Pop the stack.  If Before or Where are present, dump the top frame first.

    use ALLOCATE_DEALLOCATE, only: MEMORY_UNITS, NOBYTESALLOCATED
    use OUTPUT_M, only: NEWLINE, OUTPUT
    use HIGHOUTPUT, only: DUMPSIZE
    use STRING_TABLE, only: DISPLAY_STRING
    use TIME_M, only: TIME_NOW

    character(len=*), intent(in), optional :: Before
    logical, intent(in), optional :: Where
    type(stack_t), intent(out), optional :: Frame
    integer, intent(in), optional :: Index
    character(len=*), optional, intent(in) :: String
    integer, optional, intent(in) :: StringIndex
    logical, intent(in), optional :: Silent

    double precision :: Delta
    logical :: HaveStack
    logical :: MySilent
    real :: T
    character(len=10) :: Used

    ! Executable
    if ( Verbose ) then
      call output( 'Popping ', advance='no' )
      if ( stack(stack_ptr)%string > 0 ) &
        & call display_string ( stack(stack_ptr)%string, before=' ' )
      if ( stack(stack_ptr)%text > 0 ) &
        & call display_string ( stack(stack_ptr)%text, before=' ' )
      call newLine
    endif
    if ( StaySilent ) return
    mySilent = .false.
    if ( present(silent) ) mySilent = silent
    haveStack = allocated(stack)
    if ( haveStack ) haveStack = stack_ptr >= lbound(stack,1)

    if ( .not. mySilent .and. (present(before) .or. present(where)) ) then
      ! We call dump_stack even without haveStack because it prints
      ! an error message if we don't have a stack.
      call dump_stack ( .true., before, where, size=.false., index=index, &
        & string=string, stringIndex=stringIndex, advance='no' )
      if ( .not. haveStack ) return
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

    if ( haveStack ) then
      if ( present(frame) ) frame = stack(stack_ptr)
      stack_ptr = stack_ptr - 1
    end if

  end subroutine Pop_Stack

! -------------------------------------------------  Push_Stack_B  -----
  subroutine Push_Stack_B ( Name_I, Name_C, Root, Index, String, Before, Where, &
    & Advance )
    ! If Name_I <= 0, use Create_String ( Name_C ) to give it a value.
    ! We assume the actual argument is a SAVE variable.
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use STRING_TABLE, only: CREATE_STRING

    integer, intent(inout) :: Name_I
    character(len=*), intent(in) :: Name_C
    integer, optional, intent(in) :: Root   ! Where in configuration tree
    integer, optional, intent(in) :: Index  ! Whatever caller wants to send
    integer, optional, intent(in) :: String ! To print after Name_C
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where
    character(len=*), intent(in), optional :: Advance

    if ( name_i <= 0 ) name_i = create_string ( trim(name_c) )
    call push_stack ( name_i, Root, Index, String, Before, Where, Advance )

  end subroutine Push_Stack_B

! -------------------------------------------------  Push_Stack_C  -----
  subroutine Push_Stack_C ( Name, Root, Index, String, Before, Where, Advance )
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use STRING_TABLE, only: CREATE_STRING

    character(len=*), intent(in) :: Name
    integer, optional, intent(in) :: Root   ! Where in configuration tree
    integer, optional, intent(in) :: Index  ! Whatever caller wants to send
    integer, optional, intent(in) :: String ! To print after Name
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where
    character(len=*), intent(in), optional :: Advance

    call push_stack ( create_string ( trim(name) ), Root, Index, String, &
      & Before, Where, Advance )

  end subroutine Push_Stack_C

! -------------------------------------------------  Push_Stack_I  -----
  subroutine Push_Stack_I ( Name, Root, Index, String, Before, Where, Advance )
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use ALLOCATE_DEALLOCATE, only: TEST_ALLOCATE, NOBYTESALLOCATED
    use OUTPUT_M, only: NEWLINE, OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    use TIME_M, only: TIME_NOW

    integer, intent(in) :: Name
    integer, optional, intent(in) :: Root  ! Where in configuration tree
    integer, optional, intent(in) :: Index ! Whatever caller wants to send
    integer, optional, intent(in) :: String ! To print after Name
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where
    character(len=*), intent(in), optional :: Advance ! Default yes in Dump_Stack

    integer :: Stat
    type(stack_t), allocatable :: Temp_Stack(:)
    intrinsic :: Storage_Size
    ! Executable
    if ( StaySilent ) return

    if ( .not. allocated(stack) ) then
      ! If you allocate with lbound < 0, other stuff won't work.
      allocate ( stack(startingStackSize), stat=stat )
      if ( stat /= 0 ) then
        call output ( 'Unable to alloxate temp_stack', advance='yes' )
        StaySilent = .true.
        return
      endif
      call test_allocate ( stat, moduleName, 'Stack', &
        & ubounds=(/startingStackSize/), elementSize=storage_size(stack) )
      stack_ptr = lbound(stack,1) - 1
    end if
    if ( stack_ptr >= ubound(stack,1) ) then
      ! Must increase stack size
      ! so we double it
      ! But limit number doublings to MAXDOUBLINGS
      Stack_Doublings = Stack_Doublings + 1
      if ( Stack_Doublings > MAXDOUBLINGS ) then
        call output( 'Maximum number of stack size doublings exceeded', advance='yes' )
        StaySilent = .true.
        return
      endif
      allocate ( temp_stack(2*stack_ptr), stat=stat )
      if ( stat /= 0 ) then
        call output ( 'Unable to alloxate temp_stack', advance='yes' )
        StaySilent = .true.
        return
      endif
      call test_allocate ( stat, moduleName, 'Temp_Stack', &
        & ubounds=(/2*stack_ptr/), elementSize=storage_size(stack) )
      temp_stack(:min(stack_ptr,ubound(stack,1))) = stack
      call move_alloc ( temp_stack, stack )
    end if
    stack_ptr = stack_ptr + 1
    stack(stack_ptr) = stack_t ( memory=noBytesAllocated, text=name )
    if ( present(root) ) stack(stack_ptr)%tree = root
    if ( present(index) ) stack(stack_ptr)%index = index
    if ( present(string) ) stack(stack_ptr)%string = string
    call time_now ( stack(stack_ptr)%clock )
    if ( present(before) .or. present(where) ) &
      & call dump_stack ( .true., before, where, advance=advance )
    if ( Verbose ) then
      call output( 'Pushing ', advance='no' )
      if ( stack(stack_ptr)%string > 0 ) &
        & call display_string ( stack(stack_ptr)%string, before=' ' )
      if ( stack(stack_ptr)%text > 0 ) &
        & call display_string ( stack(stack_ptr)%text, before=' ' )
      call newLine
    endif

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

! =====     Private Procedures     =====================================

  subroutine Show_When
    use OUTPUT_M, only: OUTPUT
    character(8) :: Date
    character(10) :: Time
    call date_and_time ( date=date, time=time )
    call output ( ' at ' // time(1:2) // ':' // time(3:4) // ':' // &
      & time(5:) // ' on ' // date(1:4) // '/' // date(5:6) // '/' // date(7:8) )
  end subroutine Show_When

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
! Revision 2.19  2014/01/28 02:59:18  vsnyder
! Add StringIndex argument to PopStack and DumpStack
!
! Revision 2.18  2014/01/11 01:41:02  vsnyder
! Decruftification
!
! Revision 2.17  2014/01/09 00:25:06  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.16  2013/11/01 00:03:29  pwagner
! Added safeguards limiting stacksize, dumping corrupted stacks, etc.
!
! Revision 2.15  2013/10/26 00:41:30  vsnyder
! Cannonball polishing
!
! Revision 2.14  2013/10/02 01:28:21  vsnyder
! Add 'string' argument to pop_stack, dump_stack
!
! Revision 2.13  2013/09/24 23:27:14  vsnyder
! Use Get_Where or Print_Source to start error messages
!
! Revision 2.12  2013/09/19 23:33:10  vsnyder
! More graceful handling of zero source ref
!
! Revision 2.11  2013/09/12 03:11:18  vsnyder
! Add Index argument to Pop_Stack and Dump_Stack
!
! Revision 2.10  2013/09/04 02:49:00  vsnyder
! Add 'Frame' argument to Pop_Stack to return top frame before popping
!
! Revision 2.9  2013/08/31 02:25:57  vsnyder
! Output a blank before Index component in Dump_Stack
!
! Revision 2.8  2013/08/31 01:23:13  vsnyder
! Don't try to display strings or trees with negative indices
!
! Revision 2.7  2013/08/30 23:14:26  pwagner
! pop_stack may pop silently
!
! Revision 2.6  2013/08/30 03:55:00  vsnyder
! Add "string" to stack.  Delete "now" from stack.  Add switch to control
! printing current date and time in Dump_Stack.  More repairs in Pop_Stack.
!
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
