! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program CondenseLeakLog

  ! Condense the leak log output by mpatrol.

  use Machine, only: HP, IO_Error, GETARG

  integer, parameter :: In = 11              ! Input logical unit number
  integer, parameter :: TraceLineLen = 80

  type CharListItem
    integer :: Nalloc = 0, Nkill = 0, Talloc = 0 ! for Kill list only.  Alas ...
    !                                          If only we had Fortran 2000's
    !                                          polymorphic types now....
    type(charListItem), pointer :: Next => NULL()
    character(len=traceLineLen) :: Line
  end type CharListItem

  type StackTrace
    integer :: Allocations = 0               ! Number of allocations
    type(charListItem), pointer :: Buffers => NULL()   ! Points to the last one
    !                                          in a circular list
    character(len=traceLineLen), pointer, dimension(:) :: Lines
    type(stackTrace), pointer :: Next        ! in the circular list
    integer :: TotalSize = 0                 ! of all allocations
  end type StackTrace

  type(charListItem), pointer :: After => NULL()  ! A stack of strings
  !                                            indicating that anything after
  !                                            is to be trimmed.
  character(len=127) :: Arg                  ! From the command line
  logical :: Before = .false.                ! Trim before adding to
  !                                            StackTraces list
  logical :: Done                            ! Done reading input
  integer :: HowManyInputTraces = 0
  integer :: HowManyTraces = 0
  integer :: I
  character(len=traceLineLen), pointer, dimension(:) :: InputTrace => NULL()
  type(charListItem), pointer :: Kill => NULL()   ! A stack of strings
  !                                            indicating that a stack trace
  !                                            is to be killed (not printed)
  !                                            if any element in it has a
  !                                            substring listed here
  integer :: LinesInTrace
  integer :: NumPrint                        ! After successful trim
  type(stackTrace), pointer :: StackTraces => NULL()   ! Points to the last one
    !                                               in a circular list
  integer :: Status                          ! From I/O or Allocate
  logical :: Summarize = .false.             ! Summarize allocations for each
  !                                            stack trace instead of keeping
  !                                            and printing all of them.
  type(charListItem), pointer :: TempList => NULL()    ! For walking lists
  type(stackTrace), pointer :: TempStackTraces    ! For walking list
  type(charListItem), pointer :: TrimText => NULL()    ! A stack of strings
  !                                            indicating a trace is to be
  !                                            trimmed.

!---------------------------- RCS Ident Info -------------------------------
  character(len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! =====     Executable Statements     ================================

  numPrint = huge(numPrint) / 2
  i = hp
  do
    i = i + 1
    call getarg ( i, arg )
    if ( arg == ' ' ) then
      print *, 'Enter leak trace file name: '
      read ( *, '(a)', end=999 ) arg
      exit
    end if
    if ( arg(1:1) /= '-' ) exit
    if ( arg(1:2) == '-a' ) then
      if ( arg(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, arg )
      else
        arg = arg(3:)
      end if
      call add_to_char_list ( after, arg )
    else if ( arg(1:3) == '-b ' ) then
      before = .true.
    else if ( arg(1:2) == '-p' ) then
      if ( arg(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, arg )
      else
        arg(1:2) = ' '
      end if
      read ( arg, *, iostat=status ) numPrint
      if ( status /= 0 ) then
        call io_error ( 'While converting argument of -p', status )
        numPrint = huge(numPrint)
      end if
    else if ( arg(1:3) == '-s ' ) then
      summarize = .true.
    else if ( arg(1:2) == '-t' ) then
      if ( arg(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, arg )
      else
        arg = arg(3:)
      end if
      call add_to_char_list ( trimText, arg )
    else if ( arg(1:2) == '-x' ) then
      if ( arg(3:3) == ' ' ) then
        i = i + 1
        call getarg ( i, arg )
      else
        arg = arg(3:)
      end if
      call add_to_char_list ( kill, arg )
    else
      call getarg ( hp, arg )
      print *, 'Usage: ', trim(arg), ' [options] input_file'
      print *, ' Options: -a[ ]string: Stop printing (storing if -b selected) the'
      print *, '              stack trace if it contains "string".'
      print *, '          -b: Trim stack trace Before putting it in the database.'
      print *, '              This has the additional effect that traces arriving'
      print *, '              at the same trimmed (see -t option) place and identical'
      print *, '              for the next numPrint (see -p option) places are'
      print *, '              considered to be identical.'
      print *, '          -p[ ]numPrint: How many frames (default infinite) of'
      print *, '              the stack trace are to be printed after trimming'
      print *, '              (see -t option).  The last specification of this'
      print *, '              option has effect.'
      print *, '          -s: Summarize the number and size of allocations instead'
      print *, '              of listing them all.'
      print *, '          -t[ ]string: Trim leading stack trace items that have'
      print *, '              "string" within them.  This option can be specified'
      print *, '              any number of times.'
      print *, '          -x[ ]string: Kill (don''t print) a stack trace if it'
      print *, '              has "string" in it anywhere.'
      print *, '          -[anything else]: this output.'
      stop
    end if
  end do

  open ( in, file=arg, status='old', form='formatted', access='sequential', &
    & iostat=status )
  if ( status /= 0 ) then
    call io_error ( 'Opening leak trace file', status, arg )
    go to 999
  end if

  allocate ( inputTrace(100), stat=status )
  if ( status /= 0 ) then
    call io_error ( 'Allocating InputTrace array', status )
    go to 999
  end if

  do
    read ( in, '(a)', end=999 ) inputTrace(1)
    if ( index(inputTrace(1),'unfreed allocations') /= 0 ) exit
    if ( inputtrace(1)(1:13) == 'unfreed alloc' ) then
!     print *, 'Index failed'
      exit
    end if
  end do

  do
    call get_a_trace
    if ( linesInTrace > 2 ) call add_to_trace_list
    if ( done ) exit
  end do

  if ( associated(stackTraces) ) then
    tempStackTraces => stackTraces ! end of the circular list
    do
      tempStackTraces => tempStackTraces%next
      call dump_one_trace
      if ( associated(tempStackTraces, stackTraces) ) exit
    end do
  end if

  tempList => kill
  if ( associated(kill) ) write ( *, * )
  do while ( associated(tempList) )
    write ( *, * ) 'Killed ', tempList%nKill, ' stack traces containing ', &
      & trim(tempList%line), ' accounting for ', tempList%nAlloc, &
      & ' allocations totalling ', tempList%tAlloc, ' bytes'
    tempList => tempList%next
  end do
    
999 stop

contains

  ! -------------------------------------------  Add_To_Char_List  -----
  subroutine Add_To_Char_List ( List, Line )
    type(charListItem), pointer :: List
    character(len=*) :: Line
    allocate ( tempList, stat=status )
    if ( status /= 0 ) then
      call io_error ( 'Allocating TempList', status )
      stop
    end if
    tempList%line = line
    tempList%next => list
    list => tempList
    nullify ( tempList )
  end subroutine Add_To_Char_List

  ! ------------------------------------------  Add_To_Trace_List  -----
  subroutine Add_To_Trace_List
    integer :: I1, I2, N
    type(charListItem), pointer :: TempItem
    tempStackTraces => findTrace()
    if ( associated(tempStackTraces) ) then
      if ( .not. summarize ) then
        ! Trace is the same, add buffer to list
        allocate ( tempItem, stat=status )
        if ( status /= 0 ) then
          call io_error ( 'Allocating new "buffer" for existing trace', &
            & status )
          stop
        end if
        tempItem%line = inputTrace(1)
        tempItem%next => tempStackTraces%buffers%next
        tempStackTraces%buffers%next => tempItem
        tempStackTraces%buffers => tempItem
      end if
    else
      call create_new_trace
      tempStackTraces => stackTraces
    end if
    tempStackTraces%allocations = tempStackTraces%allocations + 1
    i1 = index(inputTrace(1),'(') + 1
    i2 = index(inputTrace(1),' bytes')
    if ( i1 /= 0 .and. i2 /= 0 ) then
      read ( inputTrace(1)(i1:i2), *, iostat=status ) n
      if ( status == 0 ) &
        & tempStackTraces%totalSize = tempStackTraces%totalSize + n
    end if
  end subroutine Add_To_Trace_List

  ! -------------------------------------------  Create_New_Trace  -----
  subroutine Create_New_Trace
    allocate ( tempStackTraces, stat=status )
    if ( status /= 0 ) then
      call io_error ( 'Allocating a new trace', status )
      stop
    end if
    allocate ( tempStackTraces%lines(linesInTrace-2), stat=status )
    if ( status /= 0 ) then
      call io_error ( 'Allocating lines in a new trace', status )
      stop
    end if
    allocate ( tempStackTraces%buffers, stat=status )
    if ( status /= 0 ) then
      call io_error ( 'Allocating "buffers" in a new trace', status )
      stop
    end if
    howManyTraces = howManyTraces + 1
    tempStackTraces%buffers%line = inputTrace(1)
    tempStackTraces%buffers%next => tempStackTraces%buffers
    tempStackTraces%lines = inputTrace(2:linesInTrace-1)
    if ( associated(stackTraces) ) then
      tempStackTraces%next => stackTraces%next
      stackTraces%next => tempStackTraces
      stackTraces => tempStackTraces
    else
      stackTraces => tempStackTraces
      stackTraces%next => tempStackTraces
    end if
  end subroutine Create_New_Trace

  ! ---------------------------------------------  Dump_One_Trace  -----
  subroutine Dump_One_Trace
    integer :: I
    integer :: Printed
    type(charListItem), pointer :: TempItem
    logical :: Trimmed
    printed = 0
    trimmed = .false.
    if ( associated(kill) ) then
      do i = 1, size(tempStackTraces%lines)
        tempItem => test(kill,tempStackTraces%lines(i))
        if ( associated(tempItem) ) then
          tempItem%nKill = tempItem%nKill + 1
          tempItem%nAlloc = tempItem%nAlloc + tempStackTraces%allocations
          tempItem%tAlloc = tempItem%tAlloc + tempStackTraces%totalSize
          return
        end if
      end do
    end if
    write ( *, * )
    if ( summarize ) then
      write ( *, * ) tempStackTraces%allocations, ' allocations totalling ', &
        & tempStackTraces%totalSize, ' bytes'
    else
      tempItem => tempStackTraces%buffers
      do
        tempItem => tempItem%next
        write ( *, '(a)' ) trim(tempItem%line)
        if ( associated(tempItem, tempStackTraces%buffers) ) exit
      end do
    end if
    do i = 1, size(tempStackTraces%lines)
      if ( .not. before .and. printed == 0 .and. &
        & associated(test(trimText,tempStackTraces%lines(i))) ) then
        trimmed = .true.
        cycle
      end if
      if ( trimmed .and. printed >= numPrint ) exit
      printed = printed + 1
      write ( *, '(a)' ) trim(tempStackTraces%lines(i))
      if ( associated(test(after,tempStackTraces%lines(i))) ) exit
    end do
  end subroutine Dump_One_Trace

  ! --------------------------------------------------  FindTrace  -----
  function FindTrace ()
    type(stackTrace), pointer :: findTrace
    integer :: I
    findTrace => stackTraces ! end of the circular list
    if ( .not. associated(findTrace) ) return
    do ! stackTraces
      findTrace => findTrace%next
      if ( size(findTrace%lines) == linesInTrace-2 ) then
        do i = 1, size(findTrace%lines)
          if ( findTrace%lines(i) /= inputTrace(i+1) ) go to 88
        end do ! i
        return
      end if
88    if ( associated(findTrace, stackTraces) ) exit
    end do ! stackTraces
    nullify ( findTrace )
  end function FindTrace

  ! ------------------------------------------------  Get_A_Trace  -----
  subroutine Get_A_Trace
    logical :: CutAfter                 ! Remove the rest of the trace
    character(len=traceLineLen) :: Line ! To read a line of the trace
    character(len=traceLineLen), pointer, dimension(:) :: TempTrace
    logical :: Trimmed                  ! A prefix of the trace was trimmed
    integer :: old_size                 ! Size of trace (before doubling)
    cutAfter = .false.
    linesInTrace = 1
    trimmed = .false.
    old_size = size(inputTrace)
    do
      if ( linesInTrace > size(inputTrace) ) then
        tempTrace => inputTrace
        allocate ( inputTrace(2*old_size), stat=status )
        if ( status /= 0 ) then
          call io_error ( 'Increasing InputTrace', status )
          stop
        end if
        inputTrace(:size(tempTrace)) = tempTrace
        deallocate ( tempTrace, stat=status )
      end if
      read ( in, '(a)', iostat=status ) line
      done = status < 0
      if ( done ) exit
      if ( status > 0 ) then
        call io_error ( 'Reading leak trace file', status, arg )
        stop
      end if
      if ( line == ' ' ) then
        if ( linesInTrace == 1 ) cycle
        exit
      end if
      if ( before .and. linesInTrace == 2 .and. &
        & associated(test(trimText,line)) ) then
        trimmed = .true.
        cycle
      end if
      if ( cutAfter .or. trimmed .and. linesInTrace > numPrint + 1 ) cycle
      inputTrace(linesInTrace) = line
      linesInTrace = linesInTrace + 1
      if ( before ) cutAfter = associated(test(after,line))
    end do
    inputTrace(linesInTrace) = ' '
    howManyInputTraces = howManyInputTraces + 1
  end subroutine Get_A_Trace

  ! -------------------------------------------------------  Test  -----
  function Test ( List, Line )
    ! Returns the element of List that is a substring of Line, or returns
    ! a disassociated pointer otherwise
    type(charListItem), pointer :: Test
    type(charListItem), pointer :: List
    character(len=*), intent(in) :: Line
    test => list
    if ( .not. associated(test) ) return
    do while ( associated(test) )
      if ( index(line, trim(test%line)) /= 0 ) return
      test => test%next
    end do
  end function Test

end program CondenseLeakLog

! $Log$
! Revision 1.1  2001/04/20 16:51:47  vsnyder
! Initial commit
!
