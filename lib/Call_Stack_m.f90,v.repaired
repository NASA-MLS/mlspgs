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
  use MLSCommon, only: NoBytesAllocated
  ! Datatype and related procedures for a Stack capable of tracing execution
  ! may dump a walkback of callers, memory usage, timings, etc.
  
  ! See also MLSMessage
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! Say_Date                Display date at each push, pop
! Say_When                Display time at each push, pop
! Show_Final_Summary      Display summary of allocates, deallocates, etc.
! Show_Sys_Memory         Display memory usage at each push, pop
! Sys_memory_ch           Memory usage units (e.g., 'kb')
! Sys_memory_convert      Memory usage units converted to kb
!
!     (data type)
! Stack_t                 The stack type

!     (subroutines and functions)
! Deallocate_Stack        Recover storage from stack array
! Dump_Stack              Print walkback
! Get_Frame               Get an item from within the stack without altering it
! Pop_Stack               Pop top frame off stack
! Invert_Stack            Invert the order of frames in the stack stack
! Push_Stack              Push new frame onto stack
! Stack_Depth             Get stack depth
! Top_Stack               Get the top frame from the stack without altering it

! === (end of toc) ===
  implicit none
  private

  public :: Stack_T
  public :: Deallocate_Stack, Dump_Stack, Get_Frame
  public :: Invert_Stack, Pop_Stack, Push_Stack, Stack_Depth, Top_Stack
  public :: Show_Sys_Memory, sys_memory_ch, sys_memory_convert, sys_memory_max
  public :: Show_Final_Summary

  type :: Stack_t
    real :: Clock = 0.0                ! Whatever Time_Now returns (CPU or ?)
    integer :: Index = -1              ! Extra info, e.g. MAF number
    double precision :: Memory = 0.0d0 ! As accounted in Allocate_Deallocate
    integer :: Sys_Memory = 0          ! In use, in kB (1024), as accounted by
                                       ! the system and accessed by Memory_Used
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
  logical, parameter :: Say_Date         = .false.
  logical, save      :: Say_When         = .false.

  ! Display system memory changes during push and pop displays
  logical, save ::    Show_Final_Summary = .false. ! Should be .false.
  logical, save ::    Show_Sys_Memory    = .true.
  character(len=2) :: sys_memory_ch      = 'kB'
  real             :: sys_memory_convert = 1.0
  real, save       :: sys_memory_max     = 0.0  ! max so far
  
  logical, save, private :: StaySilent = .false. ! E.g., in case it becomes corrupted
  logical, save, private :: Verbose    = .false. ! Print each push, pop

  ! NoBytesAllocated at previous stack operation, so we can report memory
  ! size changes that occurred during untraced procedure exits.
  double precision, save, public, protected :: PrevBytes = 0.0d0

  integer, parameter, private :: MaxDoublings = 10 ! To limit max stack size
  integer, parameter, private :: StartingStackSize = 100
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
  subroutine Dump_Stack ( Top, Before, Where, Size, SysSize, CPU, DoDepth, &
    & Rev, Index, String, StringIndex, ShowTime, Used, Advance, &
    & PrintMemoryReport, StackIsEmpty )
    ! Dump the call stack.  If Top is present and true, dump only the top
    ! frame.  If Before is present, print it after "depth" dots and before
    ! anything else.  If Where is present and true, and the Tree component
    ! is not zero, print the line and column number from the configuration
    ! file.

    use Lexer_Core, only: Print_Source
    use Output_m, only: NewLine, Output
    use HighOutput, only: DumpSize, FinalMemoryReport
    use String_Table, only: Display_String
    use Tree, only: Node_in_tree, Where_At=>Where

    logical, intent(in), optional :: Top   ! Dump only the top frame
    character(len=*), intent(in), optional :: Before ! first thing output
    logical, intent(in), optional :: Where ! Dump tree location
    logical, intent(in), optional :: Size  ! Dump memory size (default true)
    logical, intent(in), optional :: SysSize ! Dump memory size, as the system
                                           ! accounts for it, in kB (default
                                           ! Show_Sys_Memory)
    logical, intent(in), optional :: CPU   ! Print CPU (default false)
    logical, intent(in), optional :: DoDepth ! Print "depth" dots (default true)
    logical, intent(in), optional :: Rev   ! Print in reverse order (default false)
    integer, intent(in), optional :: Index ! Print this instead of from stack
    character(len=*), optional, intent(in) :: String
    integer, optional, intent(in) :: StringIndex
    logical, intent(in), optional :: ShowTime   ! Show time when we dumped
    character(len=*), optional :: Used
    character(len=*), intent(in), optional :: Advance ! Default 'yes'
    logical, intent(in), optional :: PrintMemoryReport! Report on mem usage so far
    logical, intent(out), optional :: StackIsEmpty ! Return true if stack is empty

    character(len=3) :: DoAdvance, MyAdvance
    integer :: Depth, First, I, Inc, Last
    logical :: Error
    logical :: MyCPU, MyDoDepth, MyReport, MyRev, MySize, MySysSize, MyTop, MyWhere
    logical :: MyShowTime

    ! Executable

    if ( StaySilent ) return
    error = .not. allocated(stack)
    if ( .not. error ) error = stack_ptr < lbound(stack,1)
    if ( error ) then
      if ( present( StackIsEmpty ) ) then
        StackIsEmpty = .true.
      else
        call output ( "*****" )
        if ( present(before) ) call output ( " " // before )
        call output ( " There is no stack to dump *****", advance='yes' )
      endif
      return
    end if
    if ( present( StackIsEmpty ) ) StackIsEmpty = .false.

    myAdvance = 'yes'
    myCPU = .false.; myDoDepth = .true.; myRev = .false.
    myShowTime = .false.; mySize = .true.; MySysSize = show_sys_memory
    myTop = .false.; myWhere = .false.

    if ( present(advance) ) myAdvance = advance
    if ( present(cpu) ) myCPU = cpu
    if ( present(doDepth) ) myDoDepth = doDepth
    if ( present(rev) ) myRev = rev
    if ( present(showTime) ) myShowTime = showTime
    if ( present(size) ) mySize = size
    if ( present(sysSize) ) mySysSize = sysSize
    if ( present(top) ) myTop = top
    if ( present(where) ) myWhere = Where

    if ( myRev ) then ! Dump stack bottop-up
      first = merge(stack_ptr,lbound(stack,1),myTop)
      last = stack_ptr
      inc = 1
    else              ! Dump stack top-down
      first = stack_ptr
      last = merge(stack_ptr,lbound(stack,1),myTop)
      inc = -1
    end if
    
    myReport = Show_Final_Summary
    if ( present( PrintMemoryReport ) ) myReport = myReport .or. PrintMemoryReport

    doAdvance = 'yes'
    do depth = first, last, inc
      if ( depth == last ) doAdvance = myAdvance
      do i = lbound(stack,1), merge(depth,-1,myDoDepth)
        call output ( '.' )
      end do
      if ( present(before) ) call output ( before )
      if ( stack(depth)%text > 0 ) call display_string ( stack(depth)%text )
      if ( present(index) ) then
        call output ( index, before=' ' )
      else if ( stack(depth)%index >= 0 ) then
        call output ( stack(depth)%index, before=' ' )
      end if
      if ( stack(depth)%string > 0 ) &
        & call display_string ( stack(depth)%string, before=' ' )
      if ( say_when .or. myShowTime ) then
        call show_when
        if ( present(used) ) call output ( ' used ' // trim(adjustl(used)) //  ' cpu' )
      end if
      if ( present(string) ) call output ( ' ' // string )
      if ( present(stringIndex) ) then
        if ( stringIndex > 0 ) call display_string ( stringIndex, before=' ' )
      end if
      if ( myWhere .and. node_in_tree(stack(depth)%tree) ) then
        call output ( stack(depth)%tree, before=', tree at ' )
        call print_source ( where_at(stack(depth)%tree), before=', ' )
      elseif ( myWhere ) then
        call output (' (tree reference outside range) ' )
      end if
      if ( mySysSize ) call output ( sys_memory_convert*stack(depth)%sys_memory, &
        & before = ' Sys_Memory (' // sys_memory_ch // '): ' )
      if ( mySize ) call dumpSize ( stack(depth)%memory, &
        & before = ' Memory: ' )
      if ( myCPU ) call output ( stack(depth)%clock, format='(g10.3)', &
        & before=' CPU: ' )
      call output ( '', advance=doAdvance )
      sys_memory_max = max( sys_memory_max, stack(depth)%sys_memory*1.0 )
      if ( myReport ) then
        call NewLine ( dont_make_blank_line=.true. )
        call FinalMemoryReport ( Show_Final_Summary )
      endif
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

! ----------------------------------------------------  Invert_Stack  -----
  subroutine Invert_Stack
    ! Invert the stack; i.e. topmost now bottommost
    ! 2nd from top now 2nd from bottom, etc.
    ! Method: 
    ! One-by-one pop items from input Stack, each time pushing them onto
    ! inverted stack
    ! Internal variables
    integer                               :: i, n, status
    type(stack_t)                         :: Frame
    type(stack_t), allocatable            :: InvertStack(:)
    ! Executable
    if ( .not. allocated(Stack) ) return
    n = size(Stack)
    if ( n < 2 ) return
    allocate( InvertStack(n), stat = status )
    if ( status /= 0 ) then
      print *, 'Unable to allocate InvertStack'
      return
    endif
    do i=1, n
      call Pop_Stack( Frame=Frame, silent=.true. )
      InvertStack(i) = Frame
    enddo
    call Deallocate_Stack
    allocate( Stack(n), stat = status )
    if ( status /= 0 ) then
      print *, 'Unable to re-allocate Stack in Invert_Stack'
      return
    endif
    Stack = InvertStack
    deallocate( InvertStack, stat = status )
  end subroutine Invert_Stack

! ----------------------------------------------------  Pop_Stack  -----
  subroutine Pop_Stack ( Before, Where, Frame, Index, String, StringIndex, &
    & Silent, ShowTime, SysSize, PrintMemoryReport )
    ! Pop the stack.  If Before or Where are present, dump the top frame first.

    ! use Allocate_Deallocate, only: NoBytesAllocated
    use Output_m, only: FlushOutputLines, Newline, Output
    use HighOutput, only: DumpSize
    use Memory_m, only: Memory_Used
    use String_Table, only: Display_String
    use Time_m, only: Time_Now

    character(len=*), intent(in), optional :: Before
    logical, intent(in), optional :: Where
    type(stack_t), intent(out), optional :: Frame
    integer, intent(in), optional :: Index
    character(len=*), optional, intent(in) :: String
    integer, optional, intent(in) :: StringIndex
    logical, intent(in), optional :: Silent
    logical, intent(in), optional :: ShowTime ! Show time when we dumped
    logical, intent(in), optional :: SysSize  ! Dump memory size, as the system
                                              ! accounts for it, in kB (default
                                              ! Show_Sys_Memory)
    logical, intent(in), optional :: PrintMemoryReport! Report on mem usage so far

    double precision :: Delta ! Memory change, as accounted by Allocate_Deallocate
    logical :: HaveStack
    integer :: IDelta         ! Memory change, as accounted by system in kB
    logical :: MySilent
    logical :: MySysSize
    real :: T
    integer :: Total_used     ! Memory used in kilobytes (1024)
    character(len=10) :: Used

    ! Executable
    if ( stack_ptr < 1 ) then
      print *, 'Attempt to pop an already empty stack failed'
      return
    endif
    if ( Verbose ) then
      call output( 'Popping' )
      if ( stack(stack_ptr)%string > 0 ) &
        & call display_string ( stack(stack_ptr)%string, before=' ' )
      if ( stack(stack_ptr)%text > 0 ) &
        & call display_string ( stack(stack_ptr)%text, before=' ' )
      call output ( stack(stack_ptr)%memory, before=' Memory = ' )
      call output ( noBytesAllocated, before=' NoBytesAllocated = ', &
        & advance='yes' )
    end if
    if ( StaySilent ) return
    mySilent = .false.
    if ( present(silent) ) mySilent = silent
    mySysSize = Show_Sys_Memory
    if ( present(sysSize) ) mySysSize = sysSize
    haveStack = allocated(stack)
    if ( haveStack ) haveStack = stack_ptr >= lbound(stack,1)

    if ( .not. mySilent .and. (present(before) .or. present(where)) ) then
      call time_now ( t )
      t = t - stack(stack_ptr)%clock
      write ( used, '(g10.3)' ) t
      ! We call dump_stack even without haveStack because it prints
      ! an error message if we don't have a stack.
      call dump_stack ( .true., before, where, size=.false., sysSize=.false., &
        & index=index, string=string, stringIndex=stringIndex, &
        & showTime=showTime, used=used, PrintMemoryReport=PrintMemoryReport, &
        & advance='no' )
      if ( haveStack ) then
        if ( mySysSize ) then
          call memory_used ( total=total_used )
          sys_memory_max = max( sys_memory_max, total_used*1.0 )
          iDelta = total_used - stack(stack_ptr)%sys_memory
          if ( iDelta /= 0 ) then
            call output ( sys_memory_convert*iDelta, &
              & before=', Sys_Memory (' // sys_memory_ch // ') ' // &
                & trim(merge('+', ' ', iDelta > 0 )) )
            call output ( sys_memory_convert*total_used, before =' to ' )
          end if
        end if
        delta = noBytesAllocated - stack(stack_ptr)%memory
        idelta = noBytesAllocated - prevBytes
        if ( abs(delta) + abs(idelta) > 0 ) &
          & call output ( ', Memory ' )
        if ( delta /= 0 ) then
          call dumpSize ( abs(delta), signed=.true. )
          if ( abs(idelta) > abs(delta) ) call output ( ' (traced) ' )
        end if
        if ( abs(idelta) > abs(delta) ) &
          & call dumpSize ( idelta, after= ' (untraced)' )
        if ( abs(delta) + abs(idelta) > 0 ) then
          call dumpSize ( noBytesAllocated, before = ' to ' )
        else
          call output ( '(Memory unchanged)', advance='no' )
        end if
      end if
      call newLine
      call flushOutputLines ! unbuffer even if unbuffered not requested
    end if

    if ( haveStack ) then
      if ( present(frame) ) frame = stack(stack_ptr)
      stack_ptr = stack_ptr - 1
    end if

    prevBytes = noBytesAllocated

  end subroutine Pop_Stack

! -------------------------------------------------  Push_Stack_B  -----
  subroutine Push_Stack_B ( Name_I, Name_C, Root, Index, String, Before, Where, &
    & Advance )
    ! If Name_I <= 0, use Create_String ( Name_C ) to give it a value.
    ! We assume the actual argument is a SAVE variable.
    ! Push the stack.  If Before or Where are present, dump the new top frame.
    use String_table, only: Create_String

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
    use String_Table, only: Create_String

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
    ! use Allocate_Deallocate, only: Test_Allocate, NoBytesAllocated
    use HighOutput, only: OutputNamedValue
    use Memory_m, only: Memory_Used
    use Output_m, only: Output
    use String_Table, only: Display_String
    use Time_m, only: Time_Now

    integer, intent(in) :: Name
    integer, optional, intent(in) :: Root  ! Where in configuration tree
    integer, optional, intent(in) :: Index ! Whatever caller wants to send
    integer, optional, intent(in) :: String ! To print after Name
    character(len=*), optional, intent(in) :: Before  ! Dump top stack frame after push
    logical, intent(in), optional :: Where
    character(len=*), intent(in), optional :: Advance ! Default yes in Dump_Stack

    integer :: Stat
    type(stack_t), allocatable :: Temp_Stack(:)
    integer :: Total_Used ! memory, in kB (1024), as accounted by the system
    intrinsic :: Storage_Size
    ! Executable
    if ( StaySilent ) return

    if ( .not. allocated(stack) ) then
      ! If you allocate with lbound < 0, other stuff won't work.
      allocate ( stack(startingStackSize), stat=stat )
      if ( stat /= 0 ) then
        call output ( 'Unable to allocate temp_stack', advance='yes' )
        StaySilent = .true.
        return
      elseif ( verbose ) then
        call outputnamedValue ( 'allocated stack with size', startingStackSize )
      end if
      ! call test_allocate ( stat, moduleName, 'Stack', &
      !   & ubounds=(/startingStackSize/), elementSize=storage_size(stack) / 8 )
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
      end if
      allocate ( temp_stack(2*stack_ptr), stat=stat )
      if ( stat /= 0 ) then
        call output ( 'Unable to allocate temp_stack', advance='yes' )
        StaySilent = .true.
        return
      end if
      ! The increase in size will be 2*stack_ptr - size(stack) because
      ! temp_stack's allocation will be moved to stack, which will deallocate
      ! stack (without accounting for its size).
      ! call test_allocate ( stat, moduleName, 'Temp_Stack', &
      !   & ubounds= 2*stack_ptr - size(stack), elementSize=storage_size(stack) / 8 )
      temp_stack(:min(stack_ptr,ubound(stack,1))) = stack
      call move_alloc ( temp_stack, stack )
    end if
    stack_ptr = stack_ptr + 1
    call memory_used ( total=total_used )
    sys_memory_max = max( sys_memory_max, total_used*1.0 )
    stack(stack_ptr) = stack_t ( memory=noBytesAllocated, sys_memory=total_used, &
                               & text=name )
    if ( present(root) ) stack(stack_ptr)%tree = root
    if ( present(index) ) stack(stack_ptr)%index = index
    if ( present(string) ) stack(stack_ptr)%string = string
    call time_now ( stack(stack_ptr)%clock )
    if ( present(before) .or. present(where) ) &
      & call dump_stack ( .true., before, where, advance=advance )
    if ( Verbose ) then
      call output( 'Pushing', advance='no' )
      if ( stack(stack_ptr)%string > 0 ) &
        & call display_string ( stack(stack_ptr)%string, before=' ' )
      if ( stack(stack_ptr)%text > 0 ) &
        & call display_string ( stack(stack_ptr)%text, before=' ' )
      call output ( stack(stack_ptr)%memory, before=' Memory = ' )
      call output ( noBytesAllocated, before=' NoBytesAllocated = ', &
        & advance='yes' )
    end if

    prevBytes = noBytesAllocated

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

! ----------------------------------------------------  Deallocate_Stack  -----
  subroutine Deallocate_Stack
    ! use Allocate_Deallocate, only: Test_Deallocate
    ! internal variables
    integer :: stat
    ! Executable
    if ( .not. allocated(stack) ) return
    deallocate ( stack, stat=stat )
    ! call test_deallocate ( stat, ModuleName, 'Stack' )
    stack_ptr = 0
    Stack_Doublings = 0
  end subroutine Deallocate_Stack

! =====     Private Procedures     =====================================

  subroutine Show_When
    use Output_m, only: Output
    character(8) :: Date
    character(10) :: Time
    call date_and_time ( date=date, time=time )
    call output ( ' at ' // time(1:2) // ':' // time(3:4) // ':' // &
      & time(5:) )
    if ( say_date ) call output( ' on ' // date(1:4) // '/' // date(5:6) &
      & // '/' // date(7:8) )
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
! Revision 2.40  2021/09/02 22:45:44  pwagner
! Improve memory reporting
!
! Revision 2.39  2021/08/20 15:55:27  pwagner
! Made Show_Final_Summary public
!
! Revision 2.38  2021/02/05 05:09:40  pwagner
! Prevents index error if popping an empty stack
!
! Revision 2.37  2020/04/30 23:15:46  pwagner
! Added optional arg PrintMemoryReport to Dump_Stack and Pop_Stack
!
! Revision 2.36  2019/08/19 22:00:48  pwagner
! Avoid USE-ing Allocate_Deallocate due to circular dependency
!
! Revision 2.35  2018/08/06 20:00:59  vsnyder
! Test the status Deallocate_Stack
!
! Revision 2.34  2018/08/04 00:29:16  vsnyder
! Output Index after Name or Text instead of after When
!
! Revision 2.33  2017/01/04 19:16:53  pwagner
! Optional arg StackIsEmpty to Dump_Stack prevents printing if stackIsEmpty
!
! Revision 2.32  2016/03/25 00:36:34  pwagner
! More diagnostics if verbose
!
! Revision 2.31  2016/02/11 21:16:07  pwagner
! Added Deallocate_Stack
!
! Revision 2.30  2015/09/17 22:48:50  pwagner
! Will now catch error where tree ref outside bounds
!
! Revision 2.29  2015/08/25 18:37:38  vsnyder
! Call FlushOutputLines at the end of PopStack
!
! Revision 2.28  2014/10/06 23:06:49  vsnyder
! Shorten memory change trace line
!
! Revision 2.27  2014/09/29 20:48:57  pwagner
! sys_memory_max tracks max sys memory usage
!
! Revision 2.26  2014/09/29 20:23:35  vsnyder
! Cannonball polishing, change memory usage printing
!
! Revision 2.25  2014/09/11 18:25:13  pwagner
! May show sys memory in units other than kB
!
! Revision 2.24  2014/09/04 23:32:54  vsnyder
! Add system memory tracking (default off).  Add untraced memory size changes.
! Some cannonball polishing.
!
! Revision 2.23  2014/08/05 00:18:30  pwagner
! May print memory_used instead of noBytesAllocated
!
! Revision 2.22  2014/04/21 16:57:26  pwagner
! Tries harder not reference stack if not allocated
!
! Revision 2.21  2014/03/15 00:04:40  vsnyder
! Correct misspallig in two ertor masseges
!
! Revision 2.20  2014/02/13 00:04:01  pwagner
! Revert to older appearance on Exit WALK_TREE..
!
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
