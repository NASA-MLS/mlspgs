! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program MLSL2
  use Allocate_Deallocate, only: SET_GARBAGE_COLLECTION
  use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES
  use L2PARINFO, only: PARALLEL, INITPARALLEL, ACCUMULATESLAVEARGUMENTS
  use LEXER_CORE, only: INIT_LEXER
  use LEXER_M, only: CapIdentifiers
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use MLSCOMMON, only: FILENAMELEN
  use MLSL2Options, only: PCF_FOR_INPUT, PCF, OUTPUT_PRINT_UNIT, &
    & QUIT_ERROR_THRESHOLD, TOOLKIT, CREATEMETADATA, CURRENT_VERSION_ID, &
    & PENALTY_FOR_NO_METADATA, PUNISH_FOR_INVALID_PCF, NORMAL_EXIT_STATUS, &
    & GARBAGE_COLLECTION_BY_CHUNK, &
    & DEFAULT_HDFVERSION_READ, DEFAULT_HDFVERSION_WRITE
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, &
    & ADD_TO_SECTION_TIMING, DUMP_SECTION_TIMINGS
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMSG_Debug, &
    & MLSMSG_Error, MLSMSG_Severity_to_quit, MLSMessageExit
  use MLSPCF2, only: MLSPCF_L2CF_START
  use OBTAIN_MLSCF, only: Close_MLSCF, Open_MLSCF
  use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
  use PARSER, only: CONFIGURATION
  use PVM, only: ClearPVMArgs, FreePVMArgs
  use SDPToolkit, only: UseSDPToolkit
  use SnoopMLSL2, only: SNOOPINGACTIVE
  use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
    & DESTROY_STRING_TABLE, DO_LISTING, INUNIT
  use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
  use Time_M, only: Time_Now, Use_Wall_Clock
  use TOGGLES, only: CON, EMIT, GEN, LEVELS, LEX, PAR, SYN, SWITCHES, TAB, &
    & TOGGLE
  use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE, PRINT_SUBTREE
  use TREE_CHECKER, only: CHECK_TREE
  use TREE_WALKER, only: WALK_TREE_TO_DO_MLS_L2
  use MATRIXMODULE_0, only: CHECKBLOCKS

  ! Main program for level 2 processing
  ! (It is assumed that mlsl1 has already been run successfully)
  ! Usage:
  ! mlsl2 [options] [<] [l2cf]
  ! where l2cf is an ascii file which comes from one of
  ! (i)   a file named by a line in the pcf 
  !        (if and only if pcf_for_input is TRUE)
  ! (ii)  a file named on the command line w/o the '<' redirection
  ! (iii) stdin or a file redirected as stdin using '<'
  ! and where we expand 'mlsl2 [options]' below
  ! mlsl2 -slo1 -slo2 .. --mlo1 --mlo2 ..
  !    where slon are single-letter options and mlon are multiletter options
  !    E.g., mlsl2 -m -p --nmeta -S"glo jac"
  !    For a list of available options enter 'mlsl2 --help'
  !    For a list of available switches enter 'mlsl2 -S"?"'
  ! In case the l2cf is named, say, "file_name" by (ii), we will try: 
  ! First, to find file_name in the current directory
  ! If that fails, we will try file_name.l2cf in the same directory
  ! If that fails, woe is us and we exit with an error

  ! Overview
  ! Level 2 must accomplish some operational tasks (opening and reading files),
  ! prior to any scientific tasks (processing data)
  ! The operational tasks are under the control of any and all of
  ! (a) hard-wired options (see MLSL2Options module)
  ! (b) command-line options (see /mlspgs/notes/options and 
  !         /mlspgs/notes/switches)
  ! (c) the l2cf 
  ! (d) the pcf (unless the variable 'pcf' is FALSE) 
  ! Tasks
  ! (1) Accept hard-wired options 
  ! (2) Initialize parser structures
  ! (3) Read command-line options
  ! (4) Open the l2cf (if not stdin)
  ! (5) Parse the l2cf
  ! (6) Close the l2cf (if not stdin)
  ! (7) Call walk_tree_to_do_MLS_L2 to do the (mostly) scientific tasks
  ! (8) Close down some remaining things and exit with an appropriate status

  ! Parallel processing
  ! In parallel with this, there is a possible cooperation where
  ! a single master process starts up, then tries to pass messages
  ! to slave processes which may be running. Master and slaves, all,
  ! will be running mlsl2

  implicit NONE

  integer, parameter :: L2CF_UNIT = 20  ! Unit # if L2CF is opened by Fortran

  character(len=2048) :: command_line ! All the opts
  logical :: COPYARG               ! Copy this argument to parallel command line
  logical :: COUNTCHUNKS = .false. ! Just count the chunks and quit
  logical :: CHECKL2CF = .false.   ! Just check the l2cf and quit
  logical :: DO_DUMP = .false.     ! Dump declaration table
  logical :: DUMP_TREE = .false.   ! Dump tree after parsing
  integer :: ERROR                 ! Error flag from check_tree
  integer :: FIRST_SECTION         ! Index of son of root of first n_cf node
  logical :: garbage_collection_by_dt = .false. ! Collect garbage after each deallocate_test?
  integer :: I                     ! counter for command line arguments
  integer :: J                     ! index within option
  character(len=2048) :: LINE           ! Into which is read the command args
  integer :: N                     ! Offset for start of --'s text
  integer :: ROOT                  ! of the abstract syntax tree
  integer :: SINGLECHUNK = 0       ! Just run one chunk
  integer :: STATUS                ! From OPEN
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T1, T2               ! For timing
  logical :: Timing = .false.      ! -T option is set
  character(len=FILENAMELEN) :: L2CF_file       ! Some text
  character(len=2048) :: WORD       ! Some text
  character(len=*), parameter :: L2CFNAMEEXTENSION = ".l2cf"

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
     "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  !---------------- Task (1) ------------------

  call time_now ( t0 )

! Where to send output, how severe an error to quit
   prunit = OUTPUT_PRINT_UNIT
   MLSMSG_Severity_to_quit = MAX(QUIT_ERROR_THRESHOLD, MLSMSG_Debug+1)

! Clear the command line arguments we're going to accumulate to pass
! to slave tasks
   call ClearPVMArgs
   
  !---------------- Task (2) ------------------
! Initialize the lexer, symbol table, and tree checker's tables:
!  ( Under some circumstances, you may need to increase these )
  call init_lexer ( n_chars=10000, n_symbols=1000, hash_table_size=4051 )
  call allocate_decl ( ndecls=1000 )
  call allocate_tree ( n_tree=60000 )
  call init_tables

  ! We set up a mirror command line for launching slaves
  call getarg ( hp, parallel%executable )
  i = 1+hp

  !---------------- Task (3) ------------------
  command_line = ' '
  do ! Process the command line options to set toggles
    copyArg = .true.
    call getarg ( i, line )
    command_line = trim(command_line) // ' ' // trim(line)
    if ( line(1:4) == '-Wl,' ) then     ! skip Lahey/Fujitsu run-time options
      i = i + 1
      call AccumulateSlaveArguments(line)
    end if
    if ( line(1:2) == '--' ) then       ! "word" options
      n = 0
      switch = .true.
      if ( line(3:3) == 'n' .or. line(3:3) == 'N' ) then
        switch = .false.
        n = 1
      end if
      if ( line(3+n:8+n) == 'cfpcf ' ) then
        pcf_for_input = switch
      else if ( line(3+n:7+n) == 'check ' ) then
        checkl2cf = switch
      else if ( line(3+n:7+n) == 'chunk' ) then
        copyArg=.false.
        if ( line(8+n:) /= ' ' ) then
          line(:7+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
        end if
        read ( line, *, iostat=status ) singleChunk
        if ( status /= 0 ) then
          call io_error ( "After --chunk option", status, line )
          stop
        end if
      else if ( line(3+n:7+n) == 'ckbk ' ) then
        checkBlocks = switch
      else if ( line(3+n:14+n) == 'countChunks ' ) then
        countChunks = switch
      else if ( line(3+n:7+n) == 'gcch ' ) then
        garbage_collection_by_chunk = switch
      else if ( line(3+n:7+n) == 'gcdt ' ) then
        garbage_collection_by_dt = switch
      else if ( line(3+n:6+n) == 'kit ' ) then
        MLSMessageConfig%useToolkit = switch
      else if ( line(3+n:7+n) == 'meta ' ) then
        createMetadata = switch
      else if ( line(3+n:9+n) == 'master ' ) then
        copyArg = .false.
        parallel%master = .true.
        i = i + 1
        call getarg ( i, line )
        parallel%slaveFilename = trim ( line )
        call InitParallel ( 0 )
        word = '--slave'
        write ( word(len_trim(word)+1:), * ) parallel%myTid
        call AccumulateSlaveArguments(word)
      else if ( line(3+n:6+n) == 'pcf ' ) then
        pcf = switch
      else if ( line(3+n:9+n) == 'submit ' ) then
        copyArg = .false.
        i = i + 1
        call getarg ( i, line )
        parallel%submit = trim ( line )
      else if ( line(3+n:7+n) == 'slave' ) then
        copyArg=.false.
        parallel%slave = .true.
        if ( line(8+n:) /= ' ' ) then
          line(:7+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
        end if
        read ( line, *, iostat=status ) parallel%masterTid
        if ( status /= 0 ) then
          call io_error ( "After --slave option", status, line )
          stop
        end if
      else if ( line(3+n:8+n) == 'snoop ' ) then
        snoopingActive = .true.
      else if ( line(3+n:5+n) == 'tk ' ) then
        toolkit = switch
      else if ( line(3+n:10+n) == 'version ' ) then
        do j=1, size(current_version_id)
          print*, current_version_id(j)
        enddo
        stop
      else if ( line(3+n:7+n) == 'wall ' ) then
        use_wall_clock = switch
      else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
        i = i + 1
        call getarg ( i, line )
        call AccumulateSlaveArguments(line)
        exit
      else
        print *, 'unrecognized option ', trim(line), ' ignored.'
        call option_usage
      end if
    else if ( line(1:1) == '-' ) then   ! "letter" options
      j = 1
      do while ( j < len(line) )
        j = j + 1
        select case ( line(j:j) )
        case ( ' ' )
          exit
        case ( 'A' ); dump_tree = .true.
        case ( 'a' ); toggle(syn) = .true.
        case ( 'c' ); toggle(con) = .true.
        case ( 'd' ); do_dump = .true.
        case ( 'f' )
          toggle(emit) = .true.
          levels(emit) = 0
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              j = j + 1
              levels(emit) = ichar(line(j:j)) - ichar('0')
            end if
          end if
        case ( 'g' )
          toggle(gen) = .true.
          levels(gen) = 0
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              j = j + 1
              levels(gen) = ichar(line(j:j)) - ichar('0')
            end if
          end if
        case ( 'h', 'H', '?' )     ! Describe command line usage
          call option_usage
        case ( 'K' ); capIdentifiers = .true.
        case ( 'k' ); capIdentifiers = .false.
        case ( 'l' ); toggle(lex) = .true.
        case ( 'M' ); prunit = -2
        case ( 'm' ); prunit = -1
        case ( 'p' ); toggle(par) = .true.
        case ( 'S' )
          switches = trim(switches) // line(j+1:)
          exit ! Took the rest of the string, so there can't be more options
        case ( 'T' )
          timing = .true.
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              if( line(j+1:j+1) /= '0' ) &
                & switches = trim(switches) // 'time'
              section_times = &
                & ( index(switches, 'time') /= 0 .and. (line(j+1:j+1) /= '1') )
              total_times = section_times .and. (line(j+1:j+1) /= '2')
              j = j + 1
            end if
          end if
        case ( 't' ); toggle(tab) = .true.
        case ( 'v' ); do_listing = .true.
        case default
          print *, 'Unrecognized option -', line(j:j), ' ignored.'
          call option_usage
        end select
      end do
    else    
      call AccumulateSlaveArguments(line)
      exit ! This must be the l2cf filename
    end if
    i = i + 1
    if ( copyArg ) call AccumulateSlaveArguments(line)
  end do

  if( index(switches, '?') /= 0 .or. index(switches, 'hel') /= 0 ) then
   call switch_usage
  end if
  call Set_garbage_collection(garbage_collection_by_dt)
! Done with command-line parameters; enforce cascading negative options
! (waited til here in case any were (re)set on command line)

   if ( .not. toolkit ) then
      pcf = .false.
      prunit = max(-1, prunit)   ! stdout or Fortran unit
   end if

   if( index(switches, 'log') /= 0 .or. .not. toolkit ) then
      MLSMessageConfig%LogFileUnit = -1
   else
      MLSMessageConfig%LogFileUnit = -2   ! the default in MLSMessageModule
   end if

   UseSDPToolkit = pcf    ! Redundant, but may be needed in lib

   if ( .not. pcf ) then
      pcf_for_input = .false.
      punish_for_invalid_pcf = .false.
      createMetadata = .false.
      penalty_for_no_metadata = 0
   end if

  if( index(switches, 'opt') /= 0 ) then
    call dump_settings
  endif

  ! Setup the parallel stuff.  Register our presence with the master if we're a
  ! slave.
  if ( parallel%slave ) call InitParallel ( singleChunk )
! Parse the L2CF, producing an abstract syntax tree

  !---------------- Task (4) ------------------
  status = 0
  L2CF_file = '<STDIN>'
  if ( line /= ' ' ) then
    L2CF_file = line
    open ( l2cf_unit, file=line, status='old', &
      & form='formatted', access='sequential', iostat=status )
    if ( status /= 0 ) then
      L2CF_file = trim(line) // L2CFNAMEEXTENSION
      open ( l2cf_unit, file=trim(line) // L2CFNAMEEXTENSION, status='old', &
        & form='formatted', access='sequential', iostat=status )
    end if
    if ( status /= 0 ) then
      call io_error ( "While opening L2CF", status, line )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file: " // trim(line) )
    end if
    inunit = l2cf_unit
  else if ( pcf_for_input ) then
    call open_MLSCF ( MLSPCF_L2CF_Start, inunit, L2CF_file, status )
    if(status /= 0) then
      call output( 'Non-zero status returned from open_MLSCF: ', &
      & advance='no')
      call output(status, advance='yes')
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file named in pcf" )
    end if
  end if
  error = status
  call time_now ( t1 )

  !---------------- Task (5) ------------------
  if (error == 0) then
    call configuration ( root )
  else
    root = -1
  endif
  if ( timing ) call sayTime ( 'Parsing the L2CF' )

  !---------------- Task (6) ------------------
  if ( PCF_FOR_INPUT .and. error==0) then
    call close_MLSCF ( inunit, error )
  else
    if ( inunit >= 0 ) close ( inunit )  ! Don't worry about the status
  end if
  if ( error /= 0) then
    call output ( &
      'An io error occurred with the l2cf -- there is no abstract syntax tree', &
      advance='yes' )
  elseif ( root <= 0 ) then
    call output ( &
      'A syntax error occurred -- there is no abstract syntax tree', &
      advance='yes' )
      error = 1
  else
    if ( dump_tree ) then
      call output ( 'Begin un-type-checked abstract syntax tree:', &
        & advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End un-type-checked abstract syntax tree', &
        & advance='yes' )
    end if

    ! Check that supra-syntactic conditions are met, e.g. correct
    ! types for fields of commands, correct command order, etc.
    call time_now ( t1 )
    call check_tree ( root, error, first_section )
    if(error /= 0) then
       call MLSMessage(MLSMSG_Error, ModuleName, &
       & 'error in check_tree: probably need to repair l2cf ' )
    end if
    if ( timing ) call sayTime ( 'Type checking the L2CF' )
    if ( do_dump ) call dump_decl
    if ( toggle(syn) ) then
      call output ( 'Begin type-checked abstract syntax tree:', advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End type-checked abstract syntax tree', advance='yes' )
    end if

    call add_to_section_timing( 'main', t0 )

  !---------------- Task (7) ------------------
    if ( error == 0 .and. first_section /= 0 .and. .not. checkl2cf ) then
      ! Now do the L2 processing.
      call time_now ( t1 )
      if ( timing ) call output ( "-------- Processing Begun ------ ", advance='yes' )
      call walk_tree_to_do_MLS_L2 ( root, error, first_section, countChunks, singleChunk )
      if ( timing ) then
        call output ( "-------- Processing Ended ------ ", advance='yes' )
        call sayTime ( 'Processing' )
      end if
    end if
  end if

  call time_now ( t0 )
  t1 = t0
  !---------------- Task (8) ------------------
  call destroy_char_table
  call destroy_hash_table
  call destroy_string_table
  call destroy_symbol_table
  call deallocate_decl
  call deallocate_tree
  call FreePVMArgs
  if ( timing ) call sayTime ( 'Closing and deallocating' )
  call add_to_section_timing( 'main', t0 )
  if ( index(switches, 'time') /= 0 ) call dump_section_timings
  if(error /= 0) then
     call MLSMessageExit(1)
  elseif(NORMAL_EXIT_STATUS /= 0 .and. .not. parallel%slave) then
     call MLSMessageExit(NORMAL_EXIT_STATUS)
  else                  
     call MLSMessageExit 
  endif                 

contains
  subroutine SayTime ( What )
    character(len=*), intent(in) :: What
    call time_now ( t2 )
    if ( total_times ) then
      call output ( "Total time = " )
      call output ( dble(t2), advance = 'no' )
      call blanks ( 4, advance = 'no' )
    endif
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - t1), advance = 'yes' )
  end subroutine SayTime

  ! The following offer some information on the options
  ! available either on the command-line or via the PCF
  ! Note the unashamed use of 'print' statements which
  ! are officially discouraged in favor of calls to MLSMessage.
  ! Unfortunately, we have not yet decided which method to use
  ! until *after* processing all the options.
  
  subroutine switch_usage
    print *, 'Switch usage: -S"sw1 sw2 .. swn" or -Ssw1 -Ssw2 ..'
    print *, ' where each of the swk may be one of the following'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   ! === (start of automatic usage lines) ===
    print *, '  A => AntennaPatterns'
   ! === (end of automatic usage lines) ===
    stop
  end subroutine switch_usage

  subroutine option_usage
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [L2CF-name]'
    print *, ' Options:'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   ! === (start of automatic option lines) ===
    print *, '  -A: Dump the un-decorated abstract syntax tree.'
   ! === (end of automatic option lines) ===
    stop
  end subroutine option_usage

  subroutine dump_settings
  ! Show current run-time settings resulting from
  ! command-line, MLSL2Options, etc.
    call output(' mlsl2 called with command line options: ', advance='no')
    call output(trim(command_line), advance='yes')
    call output(' l2cf file:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(L2CF_file), advance='yes')                            
    if( index(switches, 'opt1') /= 0 ) then                                 
      call output(' -------------- Summary of run time options'      , advance='no')
      call output(' -------------- ', advance='yes')
      call output(' Use toolkit panoply:                            ', advance='no')
      call blanks(4, advance='no')
      call output(toolkit, advance='yes')
      call output(' Use PCF file:                                   ', advance='no')
      call blanks(4, advance='no')
      call output(pcf, advance='yes')
      call output(' Get l2cf from pcf:                              ', advance='no')
      call blanks(4, advance='no')
      call output(pcf_for_input, advance='yes')
      call output(' Punish for errors in pcf:                       ', advance='no')
      call blanks(4, advance='no')
      call output(punish_for_invalid_pcf, advance='yes')
      call output(' Create metadata for each output file:           ', advance='no')
      call blanks(4, advance='no')
      call output(createMetadata, advance='yes')
      call output(' Punish for metadata creation errors:            ', advance='no')
      call blanks(5, advance='no')
      call output(penalty_for_no_metadata, advance='yes')
      call output(' Error threshold before halting:                 ', advance='no')
      call blanks(5, advance='no')
      call output(quit_error_threshold, advance='yes')
      call output(' Status on normal exit:                          ', advance='no')
      call blanks(5, advance='no')
      call output(normal_exit_status, advance='yes')
      call output(' Default hdf version on reads:                   ', advance='no')
      call blanks(5, advance='no')
      call output(default_hdfversion_read, advance='yes')
      call output(' Default hdf version on writes:                  ', advance='no')
      call blanks(5, advance='no')
      call output(default_hdfversion_write, advance='yes')
      call output(' Manually collect garbage after each chunk:      ', advance='no')
      call blanks(4, advance='no')
      call output(garbage_collection_by_chunk, advance='yes')
      call output(' Manually collect garbage after each deallocate: ', advance='no')
      call blanks(4, advance='no')
      call output(garbage_collection_by_dt, advance='yes')
      call output(' Is this the master task in pvm?:                ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%master, advance='yes')
      call output(' Is this a slave task in pvm?:                   ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%slave, advance='yes')
      call output(' Using wall clock instead of cpu time?:          ', advance='no')
      call blanks(4, advance='no')
      call output(use_wall_clock, advance='yes')
      call output(' Standard output unit:                           ', advance='no')
      call blanks(4, advance='no')
      call output(PrUnit, advance='yes')
      call output(' Log file unit:                                  ', advance='no')
      call blanks(4, advance='no')
      call output(MLSMessageConfig%LogFileUnit, advance='yes')
      call output(' ----------------------------------------------------------', &
        & advance='yes')
    endif
  end subroutine dump_settings

end program MLSL2

! $Log$
! Revision 2.71  2002/04/24 20:21:13  livesey
! Upped some string sizes.
!
! Revision 2.70  2002/04/24 16:54:02  livesey
! Changes for submit option
!
! Revision 2.69  2002/03/20 00:48:29  pwagner
! Option -check just checks l2cf then quits
!
! Revision 2.68  2002/02/20 00:29:04  pwagner
! Retries FN+.l2cf; tracks successful l2cf file name
!
! Revision 2.67  2002/02/12 00:25:00  pwagner
! New switch -opt[n] and new --version option
!
! Revision 2.66  2002/02/05 00:44:03  pwagner
! Added garbage collection stuff
!
! Revision 2.65  2002/01/18 18:55:25  livesey
! Added the --chunk option
!
! Revision 2.64  2002/01/09 22:56:46  livesey
! Tidied up parsing of master option.
!
! Revision 2.63  2002/01/09 00:00:48  pwagner
! Fixed small comment; added others explaining unavoidable use of print
!
! Revision 2.62  2001/12/13 23:21:20  livesey
! Added countChunks option
!
! Revision 2.61  2001/11/09 23:17:22  vsnyder
! Use Time_Now instead of CPU_TIME
!
! Revision 2.60  2001/11/09 18:12:38  livesey
! Added --nckbk option.
!
! Revision 2.59  2001/10/12 23:11:28  pwagner
! Clarified what Processing time means?
!
! Revision 2.58  2001/10/09 22:37:55  livesey
! Increased tree size and hash table size to accomodate new
! spectroscopy database from Bill
!
! Revision 2.57  2001/10/04 23:50:25  livesey
! Added the ckbk option
!
! Revision 2.56  2001/10/04 00:16:45  pwagner
! Increased hash table size; added note that size(s) may need to grow
!
! Revision 2.55  2001/09/28 23:59:20  pwagner
! Fixed various timing problems
!
! Revision 2.54  2001/09/28 17:50:30  pwagner
! MLSL2Timings module keeps timing info
!
! Revision 2.53  2001/09/19 23:43:49  livesey
! Added --snoop option
!
! Revision 2.52  2001/07/19 22:00:09  pwagner
! Better behaved when l2cf is stdin
!
! Revision 2.51  2001/07/18 23:56:24  pwagner
! Gets error from close_mlscf
!
! Revision 2.50  2001/07/18 00:16:54  pwagner
! Better control over when to exit with status
!
! Revision 2.49  2001/07/16 23:43:15  pwagner
! With settable NORMAL_EXIT_STATUS
!
! Revision 2.48  2001/05/25 01:03:47  livesey
! Working parallel version
!
! Revision 2.47  2001/05/23 23:21:51  pwagner
! Did the same for options, switches
!
! Revision 2.46  2001/05/23 22:31:30  pwagner
! Automatically updates based on notes/switches
!
! Revision 2.45  2001/05/23 21:59:43  livesey
! Interim version, almost there
!
! Revision 2.44  2001/05/23 01:44:24  livesey
! Parallel code starting to fit into place
!
! Revision 2.43  2001/05/18 01:14:21  vsnyder
! Add a 'stop' in 'switch_usage', plus cosmetic changes
!
! Revision 2.42  2001/05/17 22:34:55  pwagner
! output and MLSMessage modules cooperate better w/o toolkit; switch_usage
!
! Revision 2.41  2001/05/15 23:46:07  pwagner
! Removed 2 settings from MLSL2Opts; now in switches
!
! Revision 2.40  2001/05/11 23:47:00  pwagner
! (Re)Sets prunit depending on toolkit
!
! Revision 2.39  2001/05/11 17:34:31  vsnyder
! Improve built-in usage display
!
! Revision 2.38  2001/05/09 23:33:00  pwagner
! Sets new MLSL2Options
!
! Revision 2.37  2001/05/08 20:33:41  vsnyder
! Moved usage display into a subroutine
!
! Revision 2.36  2001/05/07 23:30:51  pwagner
! Sets USESDPToolkit
!
! Revision 2.35  2001/05/07 21:53:28  vsnyder
! Improve built-in usage display
!
! Revision 2.34  2001/05/07 21:05:03  vsnyder
! Separated '[n]pcf' and [n]cfpcf'
!
! Revision 2.33  2001/05/07 18:16:11  vsnyder
! Added "do [not] output MLSMessage messages through the toolkit" option.
! Resolved conflicts on merge.
!
! Revision 2.32  2001/05/07 17:17:31  pwagner
! Calls MLSMessage to exit with error if check_tree fails
!
! Revision 2.31  2001/05/04 22:55:36  pwagner
! Added cascading negatives; new command line options
!
! Revision 2.30  2001/05/03 01:58:52  vsnyder
! Add error checking for --slave option
!
! Revision 2.29  2001/05/02 23:22:48  livesey
! Added parallel stuff
!
! Revision 2.28  2001/05/01 17:51:47  vsnyder
! Ignore Lahey/Fujitsu run-time library's command-line arguments
!
! Revision 2.27  2001/04/28 01:44:47  vsnyder
! Provide to set levels(emit)
!
! Revision 2.26  2001/04/27 20:59:16  vsnyder
! Change -G option to -S
!
! Revision 2.25  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.24  2001/04/24 23:04:42  vsnyder
! Add -f[digit] to set toggle(emit) and levels(emit)
!
! Revision 2.23  2001/04/21 01:44:05  vsnyder
! Make the timing message prettier
!
! Revision 2.22  2001/04/21 01:42:21  vsnyder
! Add -T option for timing
!
! Revision 2.21  2001/04/20 20:44:18  pwagner
! Sets MLSMSG_Severity_to_quit
!
! Revision 2.20  2001/04/17 22:08:56  pwagner
! Sets prunit according to OUTPUT_PRINT_UNIT
!
! Revision 2.19  2001/04/16 23:46:58  pwagner
! Gets PCF_FOR_INPUT flag from MLSL2Options
!
! Revision 2.18  2001/04/06 20:11:53  vsnyder
! Add -e option
!
! Revision 2.17  2001/04/05 01:33:46  vsnyder
! Increase initial sizes for sevaral paraer tables
!
! Revision 2.16  2001/03/28 01:29:48  vsnyder
! Add description of -G to -h output
!
! Revision 2.15  2001/03/16 21:47:57  vsnyder
! Add -G option
!
! Revision 2.14  2001/03/14 18:59:03  vsnyder
! Add K and k options to control whether the lexer capitalizes identifiers
!
! Revision 2.13  2001/03/08 00:39:37  vsnyder
! Improve some debugging output
!
! Revision 2.12  2001/03/02 02:38:17  vsnyder
! Expand LINE, alphabetize USEs
!
! Revision 2.11  2001/02/28 03:01:48  vsnyder
! Make presence of L2CF-name on command line take precedence over --[n]pcf
!
! Revision 2.10  2001/02/28 02:52:32  vsnyder
! Improve usage description
!
! Revision 2.9  2001/02/28 02:44:24  vsnyder
! Identify abstract syntax tree dumps, show default --[n]pcf in usage
!
! Revision 2.8  2001/02/28 01:59:29  vsnyder
! Access Open_MLSCF and Close_MLSCF from Obtain_MLSCF instead of Open_Init
!
! Revision 2.7  2001/02/23 02:39:56  vsnyder
! Add description of --pcf option to usage instructions.
!
! Revision 2.6  2001/02/23 02:38:34  vsnyder
! Open L2CF either by PCF or by Fortran OPEN or expect it on stdin
!
! Revision 2.5  2001/02/22 23:51:00  vsnyder
! Improved usage messages
!
! Revision 2.4  2001/02/22 23:05:12  vsnyder
! Display usage if -h, -H or -? option is present.
!
! Revision 2.3  2000/10/12 00:33:47  vsnyder
! Insert CVS variables and copyright notice
!
