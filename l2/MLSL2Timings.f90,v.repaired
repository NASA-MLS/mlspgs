! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSL2Timings              !  Timings for the MLSL2 program sections
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, FinalMemoryReport
  use Call_Stack_M, only: Sys_Memory_Ch, Sys_Memory_Convert, Sys_Memory_Max
  use HighOutput, only: Banner, BeVerbose, OutputNamedValue
  use Init_Tables_Module, only: F_Additional, F_Debug, F_Options, F_Reset, &
    & F_Silent, F_Skipdirectwrites, F_Skipdirectwritesif, &
    & F_Skipretrieval, F_Skipretrievalif, F_Stamp, F_Verbose, &
    & Field_First, Field_Last
  use Intrinsic, only: L_Hours, L_Minutes, L_Seconds
  use L2Parinfo, only: Parallel
  use MLSCommon, only: MLSDebug, MLSVerbose, &
    & MLSDebugSticky, MLSVerboseSticky
  use MLSL2Options, only: L2Options, OriginalOptions, &
    & DumpOptions, DumpMacros, &
    & ProcessOptions, RestartWarnings, RestoreDefaults, RuntimeValues, &
    & SectionTimingUnits, SkipDirectwrites, SkipDirectwritesOriginal, &
    & StopAfterSection
  use MLSMessageModule, only: MLSMessageConfig, MLSMessage, MLSMessageReset, &
    & MLSMSG_Error
  use MLSStrings, only: Lowercase
  use MLSStringLists, only: BooleanValue, Catlists, GetStringElement, &
    & NumStringElements, StringElementNum, SwitchDetail
  use MoreTree, only: Get_Boolean
  use Output_M, only: StampOptions, Blanks, IsOutputSuspended, NewLine, &
    & Output, RestoreSettings, ResumeOutput, SuspendOutput
  use String_Table, only: Get_String
  use Time_M, only: Time_Now
  use Toggles, only: Switches
  use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree

  implicit none

  public :: Section_times, total_times, &
    &       Add_to_directwrite_timing, &
    &       Add_to_retrieval_timing, Add_to_section_timing, &
    &       Addphasetophasenames, &
    &       Dump_section_timings, Run_start_time
  public :: Finishtimings, Filltimings, Restarttimings
  public :: Showtimingnames
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!
!     (Data)
! run_start_time         The starting point in calculating timings;
!                          can be reset by call to restartTimings
! section_times          Show times in each section
! total_times            Show total times from start
!
!     (subroutines and functions)
! add_to_directwrite_timing  
!                        Contribute latest directwrite times to totals
! add_to_retrieval_timing 
!                        Contribute latest retrieval times to totals
! add_to_section_timing 
!                        Contribute latest section times to totals
! addPhaseToPhaseNames   Process the phase specification statement of the l2cf
! dump_section_timings   dump accumulated elapsed timings with detailed breakdown
! finishTimings          Finish accumulating timings
! fillTimings            Return accumulated timings in array
! restartTimings         Zero out accumulating timings; reinitialize flags,
!                          run_start_time
! showTimingNames        Return list of phase names, section names, or both
! === (end of toc) ===

! === (start of api) ===
! add_to_directwrite_timing ( char* section_name, real t1 )
! add_to_retrieval_timing ( char* section_name, real t1 )
! add_to_section_timing ( char* section_name, real t1, [log now_stop] )
! addPhaseToPhaseNames ( int name, int root )
! dump_section_timings
! finishTimings ( char* which, [int returnStatus] )
! fillTimings ( num timings(:), char* which, char* names, [log othersToo] )
! restartTimings ( char* which )
! char* showTimingNames ( char* which, [log othersToo] )
! === (end of api) ===

  interface filltimings
    module procedure filltimings_single
    module procedure filltimings_double
  end interface
  
  ! This module simply contains initial settings and accumulated values.

  ! The following public settings are stored here; they may be set by MLSL2
  logical          :: SECTION_TIMES = .false.  ! Show times in each section
  logical          :: TOTAL_TIMES = .false.    ! Show total times from start
  logical, save    :: FINISHEDPHASETIMES = .false.  ! times in each phase done
  logical, save    :: FINISHEDSECTIONTIMES = .false.  ! times in each section done

  logical, private   :: COUNTEMPTY = .false.     ! Any sections named ' '?
  ! dimension of the following is big enough to allow adding unknown
  ! section names and unknown retrieval names
  integer, parameter :: MAXNUMSECTIONTIMES = 60
  integer, parameter :: SECTIONNAMELEN = 32
  integer, parameter :: MAXNUMPHASES = 40
  integer, parameter :: PHASENAMESLEN = MAXNUMPHASES*SECTIONNAMELEN
  ! The next two are no longer parameters, instead calculated whenever needed
  ! (someday Fortran may allow me to call procedures during initializing)
  integer :: num_section_times ! = NumStringElements(section_names, countEmpty)
  integer :: num_retrieval_times ! = NumStringElements(retrieval_names, countEmpty)

  character(len=*), parameter        :: section_names = &
    & 'main,open_init,global_settings,signals,spectroscopy,' // &
    & 'read_apriori,merge_grids,chunk_divide,construct,fill,retrieve,join,' // &
    & 'directwrite,algebra,output,master'

  character(len=*), parameter        :: retrieval_names = &
    & 'newton_solver,cholesky_factor,cholesky_solver,cholesky_invert,' // &
    & 'baseline,hybrid,polar_linear,switching_mirror,' // &
    & 'full_fwm,fullcloud_fwm,scan_fwm,twod_scan_fwm,linear_fwm,' // &
    & 'low_cloud,high_cloud,sids,form_normeq,tikh_reg'

  character(len=*), parameter        :: directwrite_names = &
    & 'writing,waiting'
  ! This should be the number of elements in the above ----------------
  integer, parameter                 :: num_directwrite_times = 2  ! <--|
  real, dimension(MAXNUMSECTIONTIMES), &
    & save                           :: section_timings = 0.
  character(len=PHASENAMESLEN), save :: phaseNames = ' '
  integer, save :: num_phases = 0
  real, dimension(MAXNUMPHASES), save :: phase_memory = 0.
  real, dimension(MAXNUMPHASES), save :: phase_timings = 0.
  real, public, save :: ChunkMemoryMax = 0.
  real, save :: run_start_time = 0.
  integer, parameter :: MAXNAMESLENGTH = PHASENAMESLEN+LEN(section_names)

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------  add_to_directwrite_timing  -----
  subroutine add_to_directwrite_timing( section_name, t1 )
  ! Add current elapsed directwrite section time to total so far for section_name

  ! Formal arguments
    character(len=*), intent(in):: section_name   ! One of the dw. sect_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    integer                     :: elem_offset
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
      if ( trim(lowercase(section_name)) == 'restart' ) then
        call time_now ( myLastTime )
        return
      endif
      if ( present(t1) ) myLastTime = t1
      num_section_times = NumStringElements(section_names, countEmpty)
      num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
      elem_offset = num_section_times + num_retrieval_times

      elem = StringElementNum(directwrite_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_directwrite_times ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find directwrite section name ' // section_name // &
        & ' among list ' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(elem_offset+elem) = &
          & section_timings(elem_offset+elem) + t2 - myLastTime
      endif
      myLastTime = t2
      if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_directwrite_timing

  ! -----------------------------------------  add_to_phase_timing  -----
  subroutine add_to_phase_timing( phase_name, t1 )
  ! Add current elapsed phase time to total so far for last phase

  ! Method:
  ! (1) If phase is not one we've seen before (among stored phaseNames)
  !     let it augment phaseNames
  ! (2) Compute elapsed time (using either arg t1 or stored lastTime)
  ! (3) Add elapsed time to amount for last phase name (unless 1st phase ever)
  ! Formal arguments
    character(len=*), intent(in):: phase_name   ! One of the phases, e.g. 'core'
    real, optional, intent(inout)  :: t1        ! Prior time_now, then current

  ! Private
    integer                     :: elem
    integer, save               :: myLastElem = 0
    real                        :: t2
    real, save                  :: myLastTime = 0.

    ! Executable
      if ( trim(lowercase(phase_name)) == 'restart' ) then
        call time_now ( myLastTime )
        myLastElem=0
        return
      endif
      if ( present(t1) ) myLastTime = t1
      elem = StringElementNum(LowerCase(trim(phaseNames)), &
        & LowerCase(phase_name), countEmpty)
      if ( elem < 1 .and. phase_name /= ' ' ) then
        num_phases = num_phases + 1
        if ( num_phases > MAXNUMPHASES ) &
          & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Too many phases--unable add to phase name ' // phase_name // &
          & ' to list ' // trim(phaseNames) )
        phaseNames = catLists(trim(phaseNames), trim(phase_name))
        elem = num_phases
      endif
      call time_now ( t2 )
      if ( myLastElem > 0 ) then
        phase_timings(myLastElem) = &
          & phase_timings(myLastElem) + t2 - myLastTime
        phase_memory(myLastElem) = max( phase_memory(myLastElem), sys_memory_max )
      endif
      ChunkMemoryMax = max( ChunkMemoryMax, sys_memory_max )
      myLastTime = t2
      myLastElem = elem
      if ( present(t1) ) call time_now ( t1 )
      if ( len_trim(phase_name) < 1 ) return
      select case ( switchDetail(switches, 'phase') )
      case ( 0 )
        call announce_phase(trim(phase_name))
      case ( 1: )
        call banner( 'Beginning phase ' // trim(phase_name) ) 
      end select
  end subroutine add_to_phase_timing

  ! -----------------------------------------  add_to_retrieval_timing  -----
  subroutine add_to_retrieval_timing( section_name, t1 )
  ! Add current elapsed retrieval section time to total so far for section_name

  ! Formal arguments
    character(len=*), intent(in):: section_name   ! One of the retr. sect_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current

  ! Private
    integer                     :: elem
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
      if ( trim(lowercase(section_name)) == 'restart' ) then
        call time_now ( myLastTime )
        return
      endif
      if ( present(t1) ) myLastTime = t1

      num_section_times = NumStringElements(section_names, countEmpty)
      num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
      elem = StringElementNum(retrieval_names, LowerCase(section_name), countEmpty)
      if ( elem < 1 .or. elem > num_retrieval_times ) then
          call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Unable to find retrieval section name ' // section_name // &
        & ' among list ' // retrieval_names )
      else
        call time_now ( t2 )
        section_timings(num_section_times+elem) = &
          & section_timings(num_section_times+elem) + t2 - myLastTime
      endif
      myLastTime = t2
      if ( present(t1) ) call time_now ( t1 )
  end subroutine add_to_retrieval_timing

  ! -----------------------------------------------  add_to_section_timing  -----
  subroutine add_to_section_timing( section_name, t1, now_stop )
  ! Add current elapsed section time to total so far for section_name
  ! (or possibly, one of the retrieval sections)

  ! Formal arguments
    character(len=*), intent(in):: section_name   ! One of the section_names
    real, optional, intent(inout)  :: t1          ! Prior time_now, then current
    logical, intent(out), optional :: now_stop

  ! Private
    integer                     :: elem
    real                        :: t2
    real, save                  :: myLastTime = 0.

  ! Executable
      if ( trim(lowercase(section_name)) == 'restart' ) then
        call time_now ( myLastTime )
        return
      endif
    if ( present(t1) ) myLastTime = t1
    num_section_times = NumStringElements(section_names, countEmpty)
    num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
    elem = StringElementNum(section_names, LowerCase(section_name), countEmpty)
    if ( elem < 1 .or. elem > num_section_times ) then
      print *, 'elem ', elem
      print *, 'num_section_times ', num_section_times
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & 'Unable to find section name ' // section_name // &
      & ' among list ' // section_names // ',' // retrieval_names )
    else
      call time_now ( t2 )
      section_timings(elem) = section_timings(elem) + t2 - myLastTime
    endif
    myLastTime = t2
    if ( present(t1) ) call time_now ( t1 )
    if ( present(now_stop) ) then
      now_stop = ( lowerCase(section_name) == lowercase(stopAfterSection) )
      if ( len_trim(stopAfterSection) > 0 .and. MLSverbose ) then
        print *, trim(stopAfterSection), ' is section to stop after'
        print *, trim(section_name), ' is section we just finished'
        print *, ' Stop now? ', now_stop
      endif
      if ( now_stop ) call Banner ( '** Stopping after ' // &
      & trim(section_name) // ' as requested **' )
    endif
  end subroutine add_to_section_timing

  ! -----------------------------------  addPhaseToPhaseNames  -----
  ! Process the phase specification statement of the l2cf
  ! In particular, set up an entry in the phase names database
  ! to hold timing info
  ! set options controllable by runtime flags, like whether
  ! to skip retrievals, skip directwrites, verboseness, etc.
  
  ! This may be called in response to either
  ! name: phase, ..
  ! or
  ! changeSettings, ..
  subroutine addPhaseToPhaseNames ( name, root )
    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of phase subtree
    ! Local variables
    logical :: additional
    logical :: debug
    integer :: detail
    integer :: field
    integer :: field_index
    integer :: fieldValue
    logical, dimension(field_first:field_last) :: GOT
    ! integer :: interval
    integer :: keyNo
    logical, save :: LASTPHASEOVERWROTEOPTS = .false.
    character(len=1) :: null
    character(len=128) :: OPTIONS
    character(len=80) :: PHASESTRING    ! E.g., 'Core'
    character(len=80) :: BOOLEANSTRING  ! E.g., 'BAND13_OK'
    character(len=8)  :: CHUNKSTRING  ! E.g., '257'
    logical :: reset
    logical :: silent
    integer :: son
    logical :: stamp
    logical :: verbose
    ! Executable
    L2Options = OriginalOptions ! Restore Original Options
    additional = .false.
    detail = switchDetail( switches, 'phase' )
    null = achar(0)
    reset = .false.
    silent = .false.
    stamp = detail > 1 ! E.g., -Sphase2; was .false.
    skipDirectwrites = skipDirectWritesoriginal
    ! skipRetrieval = skipRetrievalOriginal
    options = ' '
    Phasestring = ' '
    Chunkstring = ' '
    got= .false.
    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)
      got(field_Index) = .true.
      select case ( field_index )
      case ( f_additional )
        additional = get_boolean ( fieldValue )
      case ( f_debug )
        debug = get_boolean ( fieldValue )
        if ( stamp ) call output( 'Processing debug field', advance='yes' )
      case ( f_options )
        call get_string ( sub_rosa(subtree(2,son)), options, strip=.true. )
      case ( f_reset )
        reset = get_boolean ( fieldValue )
      case ( f_silent )
        silent = get_boolean ( fieldValue )
      case ( f_skipDirectwrites )
        skipDirectwrites = skipDirectWritesOriginal .or. &
          & get_boolean ( fieldValue )
      case ( f_skipDirectwritesif )
        call get_string( sub_rosa(subtree(2,son)), booleanString )
        if ( stamp ) call output( 'told to skipDirectwrites: ', advance='no' )
        if ( stamp ) call output( BooleanValue ( lowercase(booleanString), &
          & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep ), &
          advance='yes' )
        skipDirectwrites = skipDirectWritesOriginal .or. &
          & BooleanValue ( lowercase(booleanString), &
          & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep )
        call output( 'skipDirectwrites: ', advance='no' )
        call output( skipDirectwrites, advance='yes' )
      case ( f_skipRetrieval )
        ! Don't turn retrievals back on if cmdline options said to skip them
        if ( L2Options%SkipRetrievalOriginal ) cycle
        L2Options%SkipRetrieval = get_boolean ( fieldValue )
      case ( f_skipRetrievalif )
        ! Don't turn retrievals back on if cmdline options said to skip them
        if ( L2Options%SkipRetrievalOriginal ) cycle
        call get_string( sub_rosa(subtree(2,son)), booleanString )
        if ( stamp ) call output( 'told to skipRetrieval: ', advance='no' )
        if ( stamp ) call output( BooleanValue ( lowercase(booleanString), &
          & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep ), &
          & advance='yes' )
        L2Options%SkipRetrieval =  BooleanValue ( lowercase(booleanString), &
          & runTimeValues%lkeys, runTimeValues%lvalues, runTimeValues%sep )
        if ( stamp ) call output( 'skipRetrieval: ', advance='no' )
        if ( stamp ) call output( L2Options%SkipRetrieval, advance='yes' )
        if ( stamp ) call output( trim(booleanString), advance='yes' )
      case ( f_stamp )
        stamp = stamp .or. get_boolean ( fieldValue )
      case ( f_verbose )
        verbose = get_boolean ( fieldValue )
        if ( stamp ) call output( 'Processing verbose field', advance='yes' )
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    if ( name > 0 ) then
      call get_string(name, phaseString)
      L2Options%CurrentPhaseName = phaseString
    endif
    ! Restore settings if last one overwrote them (unless additional)
    if ( LASTPHASEOVERWROTEOPTS .and. .not. additional ) then
      call restoredefaults
      booleanstring = processOptions( trim( L2Options%Originalcmds ), null )
      MLSDebugSticky = .false.
      MLSVerboseSticky = .false.
      call output( 'Restoring default command-line args', advance='yes' )
      call outputNamedValue( 'command_line', trim(L2Options%Command_line) )
      call outputNamedValue( 'Switches', trim(Switches) )
    endif
    ! Does this phase overwrite settings?
    LASTPHASEOVERWROTEOPTS = .false. ! If not, we won't need to restore them next time
    if ( options /= ' ' ) then
      booleanstring = processOptions( trim(options ) )
      LASTPHASEOVERWROTEOPTS = .true.
      call output( 'Overwriting default command-line args', advance='yes' )
      call outputNamedValue( 'command_line', trim(L2Options%Command_line) )
      call outputNamedValue( 'Switches', trim(Switches) )
    endif
    if( got(f_debug) ) then
      MLSDebug = debug
      MLSDebugSticky = .true.
      LASTPHASEOVERWROTEOPTS = .true.
    endif
    if( got(f_additional) ) then
      LASTPHASEOVERWROTEOPTS = .false.
    endif
    if( got(f_verbose) ) then
      MLSVerbose = verbose
      MLSVerboseSticky = .true.
      LASTPHASEOVERWROTEOPTS = .true.
    endif
    
    if ( RESTARTWARNINGS .and. name > 0 ) call MLSMessageReset( Warnings=.true. )
    ! This will cause Warnings and Errors to print the phase name
    ! where they occurred
    ! and, if we're a slave, the Chunk number
    if ( parallel%slave .and. .not. parallel%fwmParallel ) then
      write(chunkString, '(i3)') parallel%ChunkNo
      MLSMessageConfig%Info = trim(phaseString) // ' (' // &
        & trim(adjustl(chunkString)) // ') '
    elseif ( name > 0 ) then
      MLSMessageConfig%Info = phaseString
    endif
    if ( name > 0 ) MLSMessageConfig%Warning = MLSMessageConfig%Info
    if ( silent ) then
      call suspendOutput
    else
      call resumeOutput
    endif

    if ( .not. got(f_stamp) ) then
      ! No stamp, no separate headline
    elseif ( detail < 2 ) then
      StampOptions%interval = 25 ! print a separate headline once every 25 lines
    elseif ( detail < 3 ) then
      ! -Sphase2
      StampOptions%interval = 10 ! print a separate headline once every 10 lines
    else
      ! -Sphase3 or more
      StampOptions%interval = 1 ! stamp every line with time, phase name
    endif

    If ( reset ) OriginalOptions = L2Options ! Make these Original Options
    if ( BeVerbose('opt', -1) .and. .not. IsOutputSuspended() ) &
      & call DumpOptions
    
    if ( name < 1 ) return
    if ( got(f_stamp) ) then
      stampOptions%neverStamp = .false.
      stampOptions%textCode = ': : ' // trim(phaseString) // ' : :'
    else
      call restoreSettings( 'stamp' )
    endif
    call add_to_phase_timing( trim(phaseString) )
    call FinalMemoryReport
    ChunkMemoryMax = max( ChunkMemoryMax, sys_memory_max )
    call outputNamedValue( 'Resetting to 0 sys_memory which was', sys_memory_max )
    sys_memory_max     = 0.0
    if ( switchDetail ( switches, 'bool' ) > 0 ) &
      & call dumpMacros
  end subroutine addPhaseToPhaseNames

  ! -----------------------------------------------  dump_section_timings  -----
  subroutine dump_section_timings
  ! (1)dump accumulated elapsed timings for section_names
  ! (2) breakdown timings by subsection in retrieve, directwrite
  ! (3) breakdown timings by phase
  ! (4) print sys memory usage by phase and overall
  ! Private
    character(len=*), parameter     :: TIMEFORMSMALL = '(F10.2)'
    character(len=*), parameter     :: TIMEFORMBIG = '(1PE10.2)'
    character(len=*), parameter     :: PCTFORM='(F10.0)'
    logical, parameter              :: PRINTCHUNKNUMWITHPHASES=.true.
    integer                         :: elem
    integer                         :: retrElem
    character(LEN=16)               :: section_name   ! One of the section_names
    real                            :: final
    real                            :: total
    real                            :: retrFinal
    real                            :: retrTotal
    real                            :: dwFinal
    real                            :: dwTotal
    real                            :: maximumMemory
    real                            :: phaseTotal
    real                            :: percent
    real                            :: elem_memory
    real                            :: elem_time
    character(len=len(timeformbig)) :: memoryForm
    character(len=len(timeformbig)) :: timeForm
    integer                         :: timeDivisor

  ! Executable
    num_section_times = NumStringElements(section_names, countEmpty)
    num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
    select case ( sectionTimingUnits )
    case ( l_seconds )
      timeDivisor = 1
    case ( l_minutes )
      timeDivisor = 60
    case ( l_hours )
      timeDivisor = 3600
    case default
      call output ( 'Unrecognized section timing unit : ', advance='no' )
      call output ( sectionTimingUnits, advance='yes' )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Should be among list l_(seconds, minutes, hours)' )
    end select
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Level 2 section timings : ', advance='yes' )
    call blanks ( 8, advance='no' )
    call printTaskType
    call output ( '==========================================', advance='yes' )
    call output ( 'section name ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 8, advance='no' )
    call output ( 'percent ', advance='yes' )
    call finishTimings('sections')
    total = sum(section_timings(1:num_section_times)) ! + &

    call time_now ( final )
    final = max(final-run_start_time, total)
    if ( final/timeDivisor <= 0.0 ) final = 1.0       ! Just so we don't divide by 0
    if ( final/timeDivisor < 0.5 .or. final/timeDivisor > 99999.99 ) then
      TIMEFORM = TIMEFORMBIG
    else
      TIMEFORM = TIMEFORMSMALL
    endif
    if ( sys_memory_convert*maxval(phase_memory) > 99999.99 ) then
      memoryForm = timeFormBig
    else
      memoryForm = timeFormSmall
    endif
    do elem = 1, num_section_times
      elem_time = section_timings(elem)
      call GetStringElement(section_names, section_name, elem, countEmpty)
      if ( section_name == 'master' .and. .not. parallel%master ) cycle
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    enddo
    call output ( '==========================================', advance='yes' )
    percent = 100 * total / final
    call blanks ( 2, advance='no' )
    call output ( '(subtotal)', advance='no' )
    call blanks ( 5, advance='no' )
    call output ( total/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    percent = 100 * (final-total) / final
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( (final-total)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    call output ( '==========================================', advance='yes' )
    percent = 100
    call blanks ( 4, advance='no' )
    call output ( '(final)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( final/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='no' )
    call blanks ( 4, advance='no' )
    call printTaskType
    ! Subdivision of Retrieval section
    if ( STOPAFTERSECTION /= ' ' ) then
      call output ( '(Some sections skipped) ', advance='yes' )
      return
    endif
    retrElem = StringElementNum(section_names, 'retrieve', countEmpty)
    if ( retrElem == 0 ) then
      call output ( '(Illegal section name--spelling?) ', advance='yes' )
      return
    endif
    retrFinal = section_timings(retrElem) 
    if ( L2Options%SkipRetrieval ) then
      call output ( '(Retrieval section skipped) ', advance='yes' )
    elseif ( retrFinal == 0.0 ) then
      call output ( '(Retrieval section number ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( retrFinal, advance='yes' )
      call output ( '(No Retrieval section timings breakdown) ', advance='yes' )
    else
      call output ( '==========================================', advance='yes' )
      call blanks ( 8, advance='no' )
      call output ( 'Retrieval section timings : ', advance='yes' )
      call output ( '==========================================', advance='yes' )
      call output ( 'subsection name ', advance='no' )
      call blanks ( 8, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'percent of total retrieval time', advance='yes' )
      retrTotal = sum(section_timings(1+num_section_times:num_section_times+num_retrieval_times)) ! - &
      do elem = 1, num_retrieval_times  ! num_elems
          elem_time = section_timings(num_section_times+elem)
          call GetStringElement(retrieval_names, section_name, elem, countEmpty)
        percent = 100 * elem_time / retrFinal
        call output ( section_name, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
      enddo
      call blanks ( 3, advance='no' )
      call output ( '(others)', advance='no' )
      call blanks ( 6, advance='no' )
      call output ( (retrFinal-retrTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 100*(retrFinal-retrTotal)/retrFinal, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    endif

    ! Subdivision of DirectWrite section
    elem = StringElementNum(section_names, 'directwrite', countEmpty)
    if ( elem == 0 ) then
      call output ( '(Illegal section name--spelling?) ', advance='yes' )
      return
    endif
    dwFinal = section_timings(elem) 
    if ( SKIPDIRECTWRITES ) then
      call output ( '(DirectWrite section skipped) ', advance='yes' )
    elseif ( dwFinal == 0.0 ) then
      call output ( '(DirectWrite section number ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( dwFinal, advance='yes' )
      call output ( '(No DirectWrite section timings breakdown) ', &
        & advance='yes' )
    else
      call output ( '==========================================', advance='yes' )
      call blanks ( 8, advance='no' )
      call output ( 'DirectWrite section timings : ', advance='yes' )
      call output ( '==========================================', advance='yes' )
      call output ( 'subsection name ', advance='no' )
      call blanks ( 8, advance='no' )
      call output ( 'time ', advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 'percent of total DirectWrite time', advance='yes' )
      dwTotal = sum(section_timings(1+num_section_times+num_retrieval_times:)) 
      do elem = 1, num_directwrite_times
        elem_time = section_timings(num_section_times+num_retrieval_times+elem)
        call GetStringElement(directwrite_names, section_name, elem, countEmpty)
        percent = 100 * elem_time / dwFinal
        call output ( section_name, advance='no' )
        call blanks ( 2, advance='no' )
        call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, &
          & advance='no' )
        call blanks ( 2, advance='no' )
        call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
      enddo
      call blanks ( 3, advance='no' )
      call output ( '(others)', advance='no' )
      call blanks ( 6, advance='no' )
      call output ( (dwFinal-dwTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, &
        & advance='no' )
      call blanks ( 2, advance='no' )
      call output ( 100*(dwFinal-dwTotal)/dwFinal, FORMAT=PCTFORM, &
        & LOGFORMAT=PCTFORM, advance='yes' )
    endif

    ! Subdivision of Phases
    if ( num_phases == 0 ) then
      call output ( '(No phase timings breakdown) ', advance='yes' )
      return
    endif
    call finishTimings('phases')
    ! call dump( phase_timings, 'phase_timings' )
    ! call dump ( sys_memory_convert*phase_memory, sys_memory_ch )
    call output ( '==========================================', advance='yes' )
    call blanks ( 8, advance='no' )
    call output ( 'Phase timings : ', advance='yes' )
    call output ( '==========================================', advance='yes' )
    call output ( 'phase name ', advance='no' )
    call blanks ( 11, advance='no' )
    call output ( 'time ', advance='no' )
    call blanks ( 5, advance='no' )
    call output ( 'pct of total', advance='no' )
    if ( any(phase_memory > 0. ) ) call output ( '  sys_mem max', advance='no' )
    call newLine
    maximumMemory = sys_memory_convert*ChunkMemoryMax ! 0.
    phaseTotal = 0.
    do elem = 1, num_phases
      elem_time   = phase_timings(elem)
      elem_memory = sys_memory_convert*phase_memory(elem)
      call GetStringElement(trim(phaseNames), section_name, elem, countEmpty)
      percent = 100 * elem_time / final
      call output ( section_name, advance='no' )
      call blanks ( 2, advance='no' )
      call output ( elem_time/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
      call output ( percent, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='no' )
      if ( elem_memory > 0. ) then
        call blanks ( 4, advance='no' )
        call output ( elem_memory, FORMAT=MEMORYFORM, LOGFORMAT=MEMORYFORM, advance='no' )
        call blanks ( 1 )
        call output ( sys_memory_ch, advance='no' )
      endif
      if ( PRINTCHUNKNUMWITHPHASES ) then
        call blanks ( 2, advance='no' )
        call printTaskType
      else
        call output ( ' ', advance='yes' )
      endif
      phaseTotal = phaseTotal + elem_time
      maximumMemory = max ( maximumMemory, elem_memory )
    enddo
    call blanks ( 3, advance='no' )
    call output ( '(others)', advance='no' )
    call blanks ( 6, advance='no' )
    call output ( (final-phaseTotal)/timeDivisor, FORMAT=TIMEFORM, LOGFORMAT=TIMEFORM, advance='no' )
    call blanks ( 2, advance='no' )
    call output ( 100*(final-phaseTotal)/final, FORMAT=PCTFORM, LOGFORMAT=PCTFORM, advance='yes' )
    if ( maximumMemory > 0. ) then
      call output ( 'Maximum memory used ', advance='no' )
      call output ( maximumMemory, advance='no' )
      call blanks ( 1 )
      call output ( sys_memory_ch, advance='no' )
      if ( PRINTCHUNKNUMWITHPHASES ) then
        call blanks ( 2, advance='no' )
        call printTaskType
      else
        call output ( ' ', advance='yes' )
      endif
    endif

    contains
    ! Internal subprograms
    subroutine printTaskType
      if ( parallel%master ) then
        call output ( '(Master Task) ', advance='yes' )
      elseif ( parallel%slave .and. .not. parallel%fwmParallel ) then
        call output ( '(Slave: chunk ', advance='no' )
        call output ( parallel%ChunkNo, advance='no' )
        call output ( ' ) ', advance='yes' )
      else
        call output ( '(Serial Task) ', advance='yes' )
      endif
    end subroutine printTaskType
  end subroutine dump_section_timings

  ! -----------------------------------------  finishTimings  -----
  subroutine finishTimings( which, returnStatus )
  ! Finish accumulating timings
  ! Args
  character(len=*), intent(in) :: which ! phases, sections, or all
  integer, intent(out), optional :: returnStatus ! 0 unless already finished
  ! Internal variables
  logical :: sections, phases
  integer :: dwElem
  integer :: Elem
  integer :: joinElem
  integer :: status
  character(len=32) :: section_name
  ! Executable
  status = 0
  num_section_times = NumStringElements(section_names, countEmpty)
  num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
  ! Do sections, phases, both?
  sections = (StringElementNum('all,both,sections', trim(LowerCase(which)), &
    & countEmpty) > 0)
  phases = (StringElementNum('all,both,phases', trim(LowerCase(which)), &
    & countEmpty) > 0)
  joinElem = StringElementNum(section_names, 'join', countEmpty)
  dwElem = StringElementNum(section_names, 'directwrite', countEmpty)
  if ( sections .and. .not. FINISHEDSECTIONTIMES ) then
    ! A trick:
    ! The DirectWrite section doesn't automatically include the waiting time
    ! (due to lazy coding in Join) so add it in now
    elem = StringElementNum(directwrite_names, 'waiting', countEmpty)
    section_timings(dwElem) = section_timings(dwElem) + &
      & section_timings(num_section_times+num_retrieval_times+elem)
    FINISHEDSECTIONTIMES = .true.
    sections = .false.
    ! Another trick:
    ! Adjust Join section timing to exclude directwrite timings
    ! (otherwise they would be counted twice)
    if ( joinElem > 0 .and. dwElem > 0 ) then
      section_timings(joinElem) = &
      &         section_timings(joinElem) - section_timings(dwElem)
    endif
  endif
  if ( sections ) status = max(status, 1)

  if ( phases .and. .not. FINISHEDPHASETIMES ) then
    ! A trick! Add final elapsed time to last phase
    call add_to_phase_timing(' ')
    FINISHEDPHASETIMES = .true.
    ! print *, 'num_phases: ', num_phases
    do elem = 1, num_phases
      call GetStringElement(trim(phaseNames), section_name, elem, countEmpty)
      ! print *, trim(section_name), phase_timings(elem)
    enddo
    phases = .false.
  endif
  if ( phases ) status = max(status, 1)
  if ( present(returnStatus) ) returnStatus = status
  end  subroutine finishTimings

  ! -----------------------------------------  filltimings_double  -----
  subroutine filltimings_double( timings, which, names, othersToo )
  ! Return accumulated timings in array
  ! Args
  double precision, dimension(:), intent(out) :: timings
  character(len=*), intent(in) :: which      ! 'phase' or 'section'
  character(len=*), intent(in) :: names      ! phase or section names, or 'all'
  logical, intent(in), optional :: othersToo
  ! Internal variables
  real, dimension(:), pointer :: singleTimings
  integer :: timingSize
  ! Executable
  timingSize = size(timings)
  if ( timingSize < 1 ) return
  nullify(singleTimings)
  call Allocate_test ( singleTimings, timingSize, 'singleTimings', ModuleName )
  call filltimings_single( singleTimings, which, names, othersToo )
  timings = singleTimings
  call Deallocate_test ( singleTimings, 'singleTimings', ModuleName )
  end subroutine filltimings_double

  ! -----------------------------------------  filltimings_single  -----
  subroutine filltimings_single( timings, which, names, othersToo )
  ! Return accumulated timings in array
  ! Args
  real, dimension(:), intent(out) :: timings
  character(len=*), intent(in) :: which      ! 'phase' or 'section'
  character(len=*), intent(in) :: names      ! phase or section names, or 'all'
  logical, intent(in), optional :: othersToo
  ! Internal variables
  real :: final, total
  logical :: sections, phases
  integer :: k
  integer :: elem
  integer :: iName
  integer :: numNames
  character(len=MAXNAMESLENGTH) :: myNames
  character(len=SECTIONNAMELEN) :: name      ! will be lower case
  logical :: myOthersToo
  ! Executable
  myOthersToo = .false.
  if ( present( othersToo) ) myOthersToo = othersToo
  num_section_times = NumStringElements(section_names, countEmpty)
  num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
  call time_now ( final )
  final = final - run_start_time
  ! Do sections, phases, both?
  sections = (StringElementNum('all,both,sections', LowerCase(which), &
    & countEmpty) > 0)
  phases = (StringElementNum('all,both,phases', LowerCase(which), &
    & countEmpty) > 0)
  timings = 0.
  k = 0
  if ( sections ) then
    total = 0.
    if ( trim(lowercase(names)) == 'all' ) then
      myNames=section_names
    else
      myNames=names
    endif
    numNames = NumStringElements(trim(myNames), countEmpty)
    do iName=1, numNames
      call GetStringElement(trim(LowerCase(myNames)), name, iName, countEmpty)
      elem = StringElementNum(LowerCase(section_names), trim(name), countEmpty)
      if ( elem > 0 ) then
        k = min(k+1, size(timings))
        timings(k) = section_timings(elem)
        total = total + timings(k)
      endif     
    enddo
    if ( myOthersToo ) then
      k = min(k+1, size(timings))
      timings(k) = max(final, total) - total
    endif
  endif
  if ( phases ) then
    total = 0.
    if ( trim(lowercase(names)) == 'all' ) then
      myNames=trim(phaseNames)
    else
      myNames=names
    endif
    numNames = NumStringElements(trim(myNames), countEmpty)
    do iName=1, numNames
      call GetStringElement(trim(LowerCase(myNames)), name, iName, countEmpty)
      elem = StringElementNum(trim(LowerCase(phaseNames)), &
        & trim(name), countEmpty)
      if ( elem > 0 ) then
        k = min(k+1, size(timings))
        timings(k) = phase_timings(elem)
      endif     
    enddo
    if ( myOthersToo ) then
      k = min(k+1, size(timings))
      timings(k) = max(final, total) - total
    endif
  endif
    
  end subroutine filltimings_single

  ! -----------------------------------------  restartTimings  -----
  subroutine restartTimings( which )
  ! Zero out accumulating timings unless 'flags'-only
  ! re-initialize some specific flags and parameters
  ! Args
  character(len=*), intent(in) :: which ! phases, sections, or all
  ! Internal variables
  logical :: sections, phases
  ! Executable
  num_section_times = NumStringElements(section_names, countEmpty)
  num_retrieval_times = NumStringElements(retrieval_names, countEmpty)
  if ( lowercase(which) == 'flags' ) then
    FINISHEDSECTIONTIMES = .false.
    FINISHEDPHASETIMES = .false.
    return
  endif
  ! Do sections, phases, both?
  sections = (StringElementNum('all,both,sections', LowerCase(which), &
    & countEmpty) > 0)
  phases = (StringElementNum('all,both,phases', LowerCase(which), &
    & countEmpty) > 0)

  if ( sections ) then
    section_timings = 0.
    call add_to_section_timing('restart')
    call add_to_retrieval_timing('restart')
    call add_to_directwrite_timing('restart')
    call time_now(run_start_time)
    FINISHEDSECTIONTIMES = .false.
  endif
  if ( phases ) then
    phase_timings = 0.
    call add_to_phase_timing('restart')
    FINISHEDPHASETIMES = .false.
  endif
  end subroutine restartTimings
  
  function showTimingNames(which, othersToo)  result (names)
  ! Args
  character(len=*), intent(in) :: which ! phases, sections, or all
  character(len=MAXNAMESLENGTH) :: names
  logical, intent(in), optional :: othersToo
  ! Internal variables
  logical :: sections, phases, myOthersToo
  ! Executable
  myOthersToo = .false.
  if ( present(othersToo) ) myOthersToo = othersToo
  ! Do sections, phases, both?
  sections = (StringElementNum('all,both,sections', LowerCase(which), &
    & countEmpty) > 0)
  phases = (StringElementNum('all,both,phases', LowerCase(which), &
    & countEmpty) > 0)
  names = ' '
  if ( sections ) then
    names = trim(section_names)
  else
    names = ' '
  endif
  if ( phases ) then
    names = catLists(trim(names), trim(phaseNames))
  endif
  if ( myOthersToo ) names = catLists(trim(names), '(others)')
  end function showTimingNames

!=============================================================================
  subroutine announce_phase(phase_name)
    character(len=*), intent(in) :: phase_name
    character(len=*), parameter :: BLAZON=' '
    if ( phase_name == ' ' ) return
    if ( BLAZON /= ' ' ) call output ( BLAZON // '  ' , advance='no' )
    call output ( 'Beginning phase ', advance='no' )
    call output ( phase_name , advance='no' )
    if ( BLAZON /= ' ' ) then
      call output ( '  '  // BLAZON  , advance='no' )
    else
      call output ( ' ', advance='yes', DONT_STAMP=.true. )
    endif
  end subroutine announce_phase

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

END MODULE MLSL2Timings
!=============================================================================

!
! $Log$
! Revision 2.78  2020/07/29 23:35:01  pwagner
! Attempt to keep closer track of maximumMemory used
!
! Revision 2.77  2020/04/30 23:31:09  pwagner
! Add call To FinalMemoryReport at end of each phase
!
! Revision 2.76  2019/09/19 16:10:28  pwagner
! restartTimings will optionally reset only flags, not timings
!
! Revision 2.75  2019/09/05 17:52:54  pwagner
! cmdline --skipretrievals no longer reversible
!
! Revision 2.74  2019/07/17 20:21:11  pwagner
! Dont DumpOptions while output is suspended
!
! Revision 2.73  2018/10/27 01:20:08  vsnyder
! Give initial value to Reset, so undefined-checking runs don't bomb
!
! Revision 2.72  2018/09/13 20:27:46  pwagner
! Moved changeable options to new L2Options; added DumpOptions; added /reset flag to phase commands to make current l2Options permanent
!
! Revision 2.71  2018/03/14 22:16:43  pwagner
! Initialize Phasestring, Chunkstring to prevent printing nulls
!
! Revision 2.70  2017/11/30 20:52:08  pwagner
! Added optional /stamp field to phase spec
!
! Revision 2.69  2017/11/15 00:18:11  pwagner
! Stamp with phase name if /stamp; dont reset interval unless /stamp
!
! Revision 2.68  2017/10/18 00:01:08  pwagner
! interval had been a temporary local variable, instead of StampOptions component; fixed
!
! Revision 2.67  2016/07/28 01:44:41  vsnyder
! Remove unused USE
!
! Revision 2.66  2016/05/27 00:06:14  pwagner
! Should now correctly process options containing an embedded space
!
! Revision 2.65  2016/05/19 23:22:58  pwagner
! Corrected some mistakes and misspellings in comments
!
! Revision 2.64  2015/10/13 23:50:05  pwagner
! Delete unneeded debug printing
!
! Revision 2.63  2015/09/17 23:18:28  pwagner
! Repaired bug using currentPhaseName before being defined; addPhaseToPhaseNames can now be called with name=0
!
! Revision 2.62  2015/09/10 17:50:35  pwagner
! Fixed bug counting DirectWrite timing twice
!
! Revision 2.61  2015/07/16 22:11:24  pwagner
! Will print calculation of now_stop if verbose
!
! Revision 2.60  2015/05/06 20:48:14  pwagner
! Slaves prefix Warnings or worse with phase and chunk num
!
! Revision 2.59  2014/09/29 20:50:03  pwagner
! Summarizes sys memory usagee by phase
!
! Revision 2.58  2014/08/12 23:29:40  pwagner
! /additional flag added to phase commands applies to options field
!
! Revision 2.57  2014/08/06 23:33:18  vsnyder
! Remove USE for CurrentChunkNumber, which is not referenced
!
! Revision 2.56  2014/06/25 20:43:25  pwagner
! Fixed error in evaluating runTimeFlags
!
! Revision 2.55  2014/04/10 00:44:21  pwagner
! Moved currentChunkNumber, currentPhaseName from MLSL2Timings to MLSL2Options
!
! Revision 2.54  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.53  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.52  2013/11/20 00:58:08  pwagner
! Reduce printing during non-verbose processing
!
! Revision 2.51  2013/11/18 22:26:07  pwagner
! phase spec takes optional /debug /verbose fields
!
! Revision 2.50  2013/06/14 18:49:46  vsnyder
! Decruftification
!
! Revision 2.49  2013/06/13 00:43:08  pwagner
! Make actual used of interval variable set between printing header
!
! Revision 2.48  2013/06/12 02:37:49  vsnyder
! Cruft removal
!
! Revision 2.47  2013/05/17 00:53:57  pwagner
! Use dumpMacros to dump r/t macros
!
! Revision 2.46  2013/04/05 23:26:04  pwagner
! Made 'master' a 'section' for timings summary
!
! Revision 2.45  2013/02/04 22:02:28  pwagner
! Less verbose; trimmed commented-out stuff
!
! Revision 2.44  2012/07/02 20:33:41  pwagner
! -Sphasen: n > 0 to bannerize phasename; n > 1 to stamp stdout with time,phase
!
! Revision 2.43  2012/06/27 18:00:34  pwagner
! May overwrite command line options with options field to phase spec
!
! Revision 2.42  2012/04/26 23:14:30  pwagner
! Now tracks currentPhaseName and currentChunkNumber (is there a better place?)
!
! Revision 2.41  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.40  2007/09/06 23:32:47  pwagner
! Fixed repeated phase name bug in Info component of MLSMessageConfig
!
! Revision 2.39  2007/08/31 00:04:06  pwagner
! Make chunk number part of Warning string
!
! Revision 2.38  2007/08/13 17:40:55  pwagner
! Warnings and Errors automatically print phase where they occur
!
! Revision 2.37  2007/06/21 00:54:08  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.36  2007/06/04 23:24:46  pwagner
! Global TRUE skipDirectWrite, skipRetrievals not overriden by phase settings
!
! Revision 2.35  2006/07/21 20:11:54  pwagner
! Can select what section to stop after
!
! Revision 2.34  2006/06/24 23:10:17  pwagner
! Remove unneeded thing from output_m
!
! Revision 2.33  2006/06/12 18:44:25  pwagner
! Must always obey if originally told to skip
!
! Revision 2.32  2006/03/04 00:16:38  pwagner
! May skip retrieval, directWrites depending on runtime Booleans
!
! Revision 2.31  2006/02/16 00:12:02  pwagner
! Added stamp boolean field to phase asks for printing phase names, times
!
! Revision 2.30  2006/02/10 21:11:29  pwagner
! May specify skipRetrivel for particular Phases
!
! Revision 2.29  2006/01/06 01:15:32  pwagner
! Added addPhaseToPhaseNames
!
! Revision 2.28  2005/09/22 23:39:38  pwagner
! time_config and retry_config now hold config settings
!
! Revision 2.27  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.26  2004/12/14 21:40:20  pwagner
! Notes skipped sections
!
! Revision 2.25  2004/08/16 17:13:06  pwagner
! Commented out debug printing
!
! Revision 2.24  2004/08/04 23:19:58  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.23  2004/06/29 00:08:41  pwagner
! Now can fill timings
!
! Revision 2.22  2004/05/06 20:42:24  pwagner
! Announces beginning of each phase if phase switch set
!
! Revision 2.21  2004/02/10 19:26:23  pwagner
! Fixed bug in phase names
!
! Revision 2.20  2004/01/14 18:49:58  vsnyder
! Stuff to support the Algebra section
!
! Revision 2.19  2003/12/05 00:50:44  pwagner
! Numerous small bugixes related to wall clock time use
!
! Revision 2.18  2003/10/23 22:20:16  pwagner
! A few bugfixes in dump_section_timings
!
! Revision 2.17  2003/10/22 21:19:14  pwagner
! Added timings breakdown by phases
!
! Revision 2.16  2003/10/20 18:21:45  pwagner
! Timings breakdown added for directWrite
!
! Revision 2.15  2003/08/11 23:24:48  pwagner
! Chunk no. printed if slave task
!
! Revision 2.14  2003/06/09 22:51:36  pwagner
! Renamed scan_divide to chunk_divide in timings table
!
! Revision 2.13  2003/02/27 21:56:07  pwagner
! Passes LOGFORMAT along with FORMAT
!
! Revision 2.12  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.11  2002/09/24 18:15:12  pwagner
! Prepared to allow unknown names; add_to_section_timing now like retrieval
!
! Revision 2.10  2002/09/19 19:07:05  vsnyder
! Only update t1 if it's present!
!
! Revision 2.9  2002/09/18 23:56:01  vsnyder
! Call time_now at end of add_to_retrieval_timing
!
! Revision 2.8  2002/07/23 00:06:05  pwagner
! No upper-case allowed in section names
!
! Revision 2.7  2002/07/22 22:53:10  pwagner
! Added names of 2d scan model, form norm eq, and tikh reg to retrieval
!
! Revision 2.6  2001/11/27 23:34:49  pwagner
! Split forward model timings into four types
!
! Revision 2.5  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.4  2001/10/01 23:30:50  pwagner
! Fixed bug in spelling cholesky_solver
!
! Revision 2.3  2001/10/01 22:54:22  pwagner
! Added subsection timings for Retrieval section
!
! Revision 2.2  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.1  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
