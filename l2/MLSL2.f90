! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program MLSL2
  use Allocate_Deallocate, only: SET_GARBAGE_COLLECTION
  use DECLARATION_TABLE, only: ALLOCATE_DECL, DEALLOCATE_DECL, DUMP_DECL
  use INIT_TABLES_MODULE, only: INIT_TABLES
  use L2PARINFO, only: PARALLEL, INITPARALLEL, ACCUMULATESLAVEARGUMENTS
  use LEXER_CORE, only: INIT_LEXER
  use LEXER_M, only: CapIdentifiers
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use MLSCOMMON, only: FILENAMELEN, MLSFile_T
  USE MLSFiles, only: WILDCARDHDFVERSION, HDFVERSION_4, HDFVERSION_5, &
    & MLS_IO_GEN_OPENF, ADDFILETODATABASE, Deallocate_filedatabase
  !    & MLSFile_T, MLS_IO_GEN_OPENF, ADDFILETODATABASE, Deallocate_filedatabase
  use MLSL2Options, only: OUTPUT_PRINT_UNIT, &
    & QUIT_ERROR_THRESHOLD, TOOLKIT, CURRENT_VERSION_ID, &
    & PENALTY_FOR_NO_METADATA, NORMAL_EXIT_STATUS, &
    & GARBAGE_COLLECTION_BY_CHUNK, &
    & DEFAULT_HDFVERSION_READ, &
    & DEFAULT_HDFVERSION_WRITE, DEFAULT_HDFVERSION_READ, LEVEL1_HDFVERSION, &
    & SKIPRETRIEVAL, CHECKPATHS
  use MLSL2Timings, only: SECTION_TIMES, TOTAL_TIMES, &
    & ADD_TO_SECTION_TIMING, DUMP_SECTION_TIMINGS
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMSG_Debug, &
    & MLSMSG_Error, MLSMSG_Severity_to_quit, MLSMessageExit
  use MLSPCF2, only: MLSPCF_L2CF_START
  use MLSStrings, only: GetUniqueList
  use OBTAIN_MLSCF, only: Close_MLSCF, Open_MLSCF
  use OUTPUT_M, only: BLANKS, OUTPUT, PRUNIT
  use PARSER, only: CONFIGURATION
  use PVM, only: ClearPVMArgs, FreePVMArgs
  use SDPToolkit, only: UseSDPToolkit, PGSD_IO_GEN_RSEQFRM
  use SnoopMLSL2, only: SNOOPINGACTIVE, SNOOPNAME
  use STRING_TABLE, only: DESTROY_CHAR_TABLE, DESTROY_HASH_TABLE, &
    & DESTROY_STRING_TABLE, DO_LISTING, INUNIT
  use SYMBOL_TABLE, only: DESTROY_SYMBOL_TABLE
  use Time_M, only: Time_Now, Use_Wall_Clock
  use TOGGLES, only: CON, EMIT, GEN, LEVELS, LEX, PAR, SYN, SWITCHES, TAB, &
    & TOGGLE
  use TREE, only: ALLOCATE_TREE, DEALLOCATE_TREE, PRINT_SUBTREE
  use TREE_CHECKER, only: CHECK_TREE
  use TREE_WALKER, only: WALK_TREE_TO_DO_MLS_L2
  use MATRIXMODULE_0, only: CHECKBLOCKS, SUBBLOCKLENGTH

  ! === (start of toc) ===
  !     c o n t e n t s
  !     - - - - - - - -

  ! Main program for level 2 processing
  ! === (end of toc) ===
  ! (It is assumed that mlsl1 has already been run successfully)
  ! Usage:
  ! mlsl2 [options] [<] [l2cf]
  ! where l2cf is an ascii file which comes from one of
  ! (i)   a file named by a line in the pcf 
  !        (if and only if toolkit is TRUE)
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
  ! (1) Accept hard-wired options (see module MLSL2OPtions.f90)
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
  ! The interprocess communication is handled by pvm

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
  character(len=2048) :: LINE      ! Into which is read the command args
  integer :: N                     ! Offset for start of --'s text
  integer :: NUMFILES
  integer :: NUMSWITCHES
  integer :: RECL = 20000          ! Record length for l2cf (but see --recl opt)
  integer :: RECORD_LENGTH
  integer :: ROOT                  ! of the abstract syntax tree
  integer :: SINGLECHUNK = 0       ! Just run one chunk; unless lastChunk nonzero
  integer :: LastCHUNK = 0         ! Just run range [SINGLECHUNK, LastCHUNK]
  integer :: SLAVEMAF = 0          ! Slave MAF for fwmParallel mode
  integer :: STATUS                ! From OPEN
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T1, T2               ! For timing
  logical :: Timing = .false.      ! -T option is set
  character(len=FILENAMELEN) :: L2CF_file       ! Some text
  character(len=len(switches)) :: tempSwitches
  character(len=2048) :: WORD      ! Some text
  character(len=1) :: arg_rhs      ! 'n' part of 'arg=n'
  character(len=*), parameter :: L2CFNAMEEXTENSION = ".l2cf"

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
     "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
  type (MLSFile_T)                        ::     theFile
  nullify (filedatabase)
  !---------------- Task (1) ------------------

  call time_now ( t0 )
  call h5open_f(error)
  if (error /= 0) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to h5_open_f" )
  endif
  ! Before looking at command-line options, TOOLKIT is set to SIPS_VERSION
  ! So here's a good place to put any SIPS-specific settings overriding defaults
  if ( TOOLKIT ) then
    ! SIPS_VERSION
    parallel%maxFailuresPerMachine = 2
    parallel%maxFailuresPerChunk = 1
    switches=''
    DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
  else
    ! SCF_VERSION
    switches='0sl'
  endif
! Initialize the lexer, symbol table, and tree checker's tables:
!  ( Under some circumstances, you may need to increase these )
  call init_lexer ( n_chars=80000, n_symbols=4000, hash_table_size=611957 )
  call allocate_decl ( ndecls=8000 )
  call allocate_tree ( n_tree=360000 )
  call init_tables

  !---------------- Task (2) ------------------
! Where to send output, how severe an error to quit
   prunit = OUTPUT_PRINT_UNIT
   MLSMSG_Severity_to_quit = MAX(QUIT_ERROR_THRESHOLD, MLSMSG_Debug+1)

! Clear the command line arguments we're going to accumulate to pass
! to slave tasks
   call ClearPVMArgs

  ! We set up a mirror command line for launching slaves
  call getarg ( hp, parallel%executable )

  !---------------- Task (3) ------------------
  i = 1+hp
  command_line = ' '
  do ! Process Lahey/Fujitsu run-time options; they begin with "-Wl,"
    call getarg ( i, line )
    if ( line(1:4) /= '-Wl,' ) exit
    call AccumulateSlaveArguments(line) ! pass them to slave processes
    i = i + 1
  end do
  do ! Process the command line options to set toggles
    copyArg = .true.
    call getarg ( i, line )
    command_line = trim(command_line) // ' ' // trim(line)
    if ( line(1:2) == '--' ) then       ! "word" options
      n = 0
      switch = .true.
      if ( line(3:3) == 'n' .or. line(3:3) == 'N' ) then
        switch = .false.
        n = 1
      end if
      if ( line(3+n:7+n) == 'check ' ) then
        checkl2cf = switch
      else if ( line(3+n:12+n) == 'checkpaths' ) then
        checkPaths = switch
      else if ( line(3+n:7+n) == 'chunk' ) then
        call AccumulateSlaveArguments ( line )
        if ( line(8+n:) /= ' ' ) then
          copyArg = .false.
          line(:7+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
        end if
        if ( singleChunk == 0 ) then
          read ( line, *, iostat=status ) singleChunk
          if ( status /= 0 ) then
            call io_error ( "After --chunk option", status, line )
            stop
          end if
        else
          read ( line, *, iostat=status ) lastChunk
          if ( status /= 0 ) then
            call io_error ( "After --chunk option", status, line )
            stop
          end if
        end if
      else if ( line(3+n:7+n) == 'ckbk ' ) then
        checkBlocks = switch
      else if ( line(3+n:14+n) == 'countChunks ' ) then
        countChunks = switch
      else if ( line(3+n:14+n) == 'fwmParallel ' ) then
        parallel%fwmParallel = .true.
      else if ( line(3+n:7+n) == 'gcch ' ) then
        garbage_collection_by_chunk = switch
      else if ( line(3+n:7+n) == 'gcdt ' ) then
        garbage_collection_by_dt = switch
      else if ( line(3+n:6+n) == 'kit ' ) then
        MLSMessageConfig%useToolkit = switch
      else if ( line(3+n:6+n) == 'l1b=' ) then
        arg_rhs = line(7+n:7+n)
        select case(arg_rhs)
        case('4')
          LEVEL1_HDFVERSION = HDFVERSION_4
        case('5')
          LEVEL1_HDFVERSION = HDFVERSION_5
        case default
          LEVEL1_HDFVERSION = WILDCARDHDFVERSION
        end select
      else if ( line(3+n:8+n) == 'l2gpr=' ) then
        arg_rhs = line(9+n:9+n)
        select case(arg_rhs)
        case('4')
          DEFAULT_HDFVERSION_READ = HDFVERSION_4
        case('5')
          DEFAULT_HDFVERSION_READ = HDFVERSION_5
        case default
          DEFAULT_HDFVERSION_READ = WILDCARDHDFVERSION
        end select
      else if ( line(3+n:8+n) == 'l2gpw=' ) then
        arg_rhs = line(9+n:9+n)
        select case(arg_rhs)
        case('4')
          DEFAULT_HDFVERSION_WRITE = HDFVERSION_4
        case('5')
          DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
        case default
          DEFAULT_HDFVERSION_WRITE = HDFVERSION_5
        end select
      else if ( line(3+n:9+n) == 'master ' ) then
        copyArg = .false.
        parallel%master = .true.
        i = i + 1
        call getarg ( i, line )
        command_line = trim(command_line) // ' ' // trim(line)
        parallel%slaveFilename = trim ( line )
        call InitParallel ( 0, 0 )
        word = '--slave'
        write ( word(len_trim(word)+1:), * ) parallel%myTid
        call AccumulateSlaveArguments(word)
      else if ( line(3+n:5+n) == 'pge ' ) then
        i = i + 1
        call getarg ( i, line )
        parallel%pgeName = trim(line)
        command_line = trim(command_line) // ' ' // trim(adjustl(line))
      else if ( line(3+n:6+n) == 'recl' ) then
        if ( line(7+n:) /= ' ' ) then
          line(:6+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
        end if
        read ( line, *, iostat=status ) recl
        if ( status /= 0 ) then
          call io_error ( "After --recl option", status, line )
          stop
        end if
      else if ( line(3+n:10+n) == 'skipRetr' ) then
        SKIPRETRIEVAL = switch
      else if ( line(3+n:10+n) == 'slaveMAF' ) then
        copyArg=.false.
        if ( line(11+n:) /= ' ' ) then
          line(:10+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(adjustl(line))
        end if
        read ( line, *, iostat=status ) slaveMAF
        if ( status /= 0 ) then
          call io_error ( "After --slaveMAF option", status, line )
          stop
        end if
      else if ( line(3+n:7+n) == 'slave' ) then
        copyArg=.false.
        parallel%slave = .true.
        if ( line(8+n:) /= ' ' ) then
          line(:7+n) = ' '
        else
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(adjustl(line))
        end if
        read ( line, *, iostat=status ) parallel%masterTid
        if ( status /= 0 ) then
          call io_error ( "After --slave option", status, line )
          stop
        end if
      else if ( line(3+n:8+n) == 'snoop ' ) then
        snoopingActive = .true.
      else if ( line(3+n:12+n) == 'snoopname' ) then
        call AccumulateSlaveArguments ( line )
        i = i + 1
        call getarg ( i, line )
        command_line = trim(command_line) // ' ' // trim(line)
        snoopName = line
      else if ( line(3+n:9+n) ==  'stgmem ' ) then
        parallel%stageInMemory = .true.
      else if ( line(3+n:11+n) == 'subblock ' ) then
        call AccumulateSlaveArguments ( line )
        i = i + 1
        call getarg ( i, line )
        command_line = trim(command_line) // ' ' // trim(line)
        read ( line, *, iostat=status ) subBlockLength
        if ( status /= 0 ) then
          call io_error ( "After --subblock option", status, line )
          stop
        end if
      else if ( line(3+n:9+n) == 'submit ' ) then
        copyArg = .false.
        i = i + 1
        call getarg ( i, line )
        command_line = trim(command_line) // ' ' // trim(line)
        parallel%submit = trim ( line )
      else if ( line(3+n:5+n) == 'tk ' ) then
        toolkit = switch
      else if ( line(3+n:10+n) == 'version ' ) then
        do j=1, size(current_version_id)
          print *, current_version_id(j)
        end do
        stop
      else if ( line(3+n:7+n) == 'wall ' ) then
        use_wall_clock = switch
      else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
        i = i + 1
        call getarg ( i, line )
        command_line = trim(command_line) // ' ' // trim(line)
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
          switches = trim(switches) // ',' // line(j+1:)
          exit ! Took the rest of the string, so there can't be more options
        case ( 'T' )
          timing = .true.
          if ( j < len(line) ) then
            if ( line(j+1:j+1) >= '0' .and. line(j+1:j+1) <= '9' ) then
              if( line(j+1:j+1) /= '0' ) &
                & switches = trim(switches) // ',' // 'time'
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
  call GetUniqueList(switches, tempSwitches, numSwitches, countEmpty=.true., &
        & ignoreLeadingSpaces=.true.)
  switches = tempSwitches
  call Set_garbage_collection(garbage_collection_by_dt)
! Done with command-line parameters; enforce cascading negative options
! (waited til here in case any were (re)set on command line)

  if ( .not. toolkit ) then
     prunit = max(-1, prunit)   ! stdout or Fortran unit
  elseif (parallel%slave .or. parallel%master) then
     prunit = -2          ! output both logged and sent to stdout
  end if

  if( index(switches, 'log') /= 0 .or. .not. toolkit ) then
     MLSMessageConfig%LogFileUnit = -1
  else
     MLSMessageConfig%LogFileUnit = -2   ! the default in MLSMessageModule
  end if

  UseSDPToolkit = toolkit    ! Redundant, but may be needed in lib

  if ( .not. toolkit ) then
     penalty_for_no_metadata = 0
  end if
  
  ! If checking paths, run as a single-chunk case in serial mode
  if ( checkPaths ) then
    parallel%master = .false.
    parallel%slave = .false.
    singleChunk = 1
    lastChunk = 0
  endif
  ! Setup the parallel stuff.  Register our presence with the master if we're a
  ! slave.
  if ( parallel%master .and. parallel%myTid <= 0 ) &
    & call MLSMessage ( MLSMSG_Error, ModuleName, &
    & 'master Tid <= 0; probably pvm trouble' )
  if ( parallel%fwmParallel .and. parallel%master .and. singleChunk == 0 ) &
    & call MLSMessage ( MLSMSG_Error, ModuleName, &
    & 'fwmParallel mode can only be run for a single chunk' )
  if ( parallel%slave ) then
    if ( parallel%masterTid <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'masterTid of this slave <= 0' )
    call InitParallel ( singleChunk, slaveMAF )
    if ( parallel%myTid <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'slave Tid <= 0; probably pvm trouble' )
  endif
  !---------------- Task (4) ------------------
  ! Open the L2CF
  status = 0
  L2CF_file = '<STDIN>'
  if ( line /= ' ' ) then
    L2CF_file = line
    open ( l2cf_unit, file=line, status='old', &
     & form='formatted', access='sequential', recl=recl, iostat=status )
    ! inunit = mls_io_gen_openF('op', .true., status, &
     ! & record_length, PGSd_IO_Gen_RSeqFrm, FileName=trim(line), &
     ! & inp_rec_length=recl)
    if ( status /= 0 ) then
      L2CF_file = trim(line) // L2CFNAMEEXTENSION
       open ( l2cf_unit, file=trim(line) // L2CFNAMEEXTENSION, status='old', &
        & form='formatted', access='sequential', recl=recl, iostat=status )
    ! inunit = mls_io_gen_openF('op', .true., status, &
     ! & record_length, PGSd_IO_Gen_RSeqFrm, &
     ! & FileName=trim(line) // L2CFNAMEEXTENSION, &
     ! & inp_rec_length=recl)
    end if
    if ( status /= 0 ) then
      call io_error ( "While opening L2CF", status, line )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file: " // trim(line) )
    elseif(index(switches, 'pro') /= 0) then                            
      call announce_success(L2CF_file, l2cf_unit)               
    end if
    inunit = l2cf_unit
  else if ( TOOLKIT ) then
    call open_MLSCF ( MLSPCF_L2CF_Start, inunit, L2CF_file, status, recl )
    if(status /= 0) then
      call output( 'Non-zero status returned from open_MLSCF: ', &
      & advance='no')
      call output(status, advance='yes')
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open L2CF file named in pcf" )
    elseif(index(switches, 'pro') /= 0) then                            
      call announce_success(L2CF_file, inunit)               
    end if
  end if
  error = status
  numfiles = AddFileToDataBase(filedatabase, theFile)
  
  call time_now ( t1 )

  if( index(switches, 'opt') /= 0 ) then
    do j=1, size(current_version_id)
      call output(trim(current_version_id(j)), advance='yes')
    end do
    call dump_settings
  end if

  !---------------- Task (5) ------------------
  if (error == 0) then
    ! Parse the L2CF, producing an abstract syntax tree
    call configuration ( root )
  else
    root = -1
  endif
  if ( timing ) call sayTime ( 'Parsing the L2CF' )

  !---------------- Task (6) ------------------
  if ( TOOLKIT .and. error==0) then
    call close_MLSCF ( inunit, error )
  else
    if ( inunit >= 0 ) close ( inunit )  ! Don't worry about the status
  end if
  if ( error /= 0) then
    call output ( &
      'An io error occurred with the l2cf -- there is no abstract syntax tree', &
      advance='yes' )
  else if ( root <= 0 ) then
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
    if ( timing ) call sayTime ( 'Type checking the L2CF' )
    if ( do_dump ) call dump_decl
    if ( toggle(syn) ) then
      call output ( 'Begin type-checked abstract syntax tree:', advance='yes' )
      call print_subtree ( root, 0 )
      call output ( 'End type-checked abstract syntax tree', advance='yes' )
    end if

    t2 = t0
    call add_to_section_timing( 'main', t2 )

    if(error /= 0) then
       call MLSMessage(MLSMSG_Error, ModuleName, &
       & 'error in check_tree: probably need to repair l2cf ' )
    end if

  !---------------- Task (7) ------------------
    if ( error == 0 .and. first_section /= 0 .and. .not. checkl2cf ) then
      ! Now do the L2 processing.
      call time_now ( t1 )
      if ( timing ) call output ( "-------- Processing Begun ------ ", advance='yes' )
      call walk_tree_to_do_MLS_L2 ( root, error, first_section, countChunks, &
        & singleChunk, lastChunk, filedatabase )
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
  call Deallocate_filedatabase(filedatabase)
  if ( timing ) call sayTime ( 'Closing and deallocating' )
  call add_to_section_timing( 'main', t0 )
  if ( index(switches, 'time') /= 0 ) call dump_section_timings
  if ( error == 0 ) then
     call h5close_f(error)
     if (error /= 0) then
       call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to h5_close_f" )
     endif    
  endif    
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
  
  subroutine Switch_usage
    print *, 'Switch usage: -S"sw1,sw2,..,swn" or -Ssw1 -Ssw2 ..'
    print *, ' where each of the swk may be one of the following'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   !  See mlspgs/notes/switches
   ! === (start of automatic usage lines) ===
    print *, '  A => AntennaPatterns'
   ! === (end of automatic usage lines) ===
    stop
  end subroutine Switch_usage

  subroutine Option_usage
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [L2CF-name]'
    print *, ' Options:'
   ! (This incorporates automatic source code replacement by
   !  a custom build command in the Makefile --
   !  Please don't edit/remove the following 3 lines)
   !  See mlspgs/notes/options
   ! === (start of automatic option lines) ===
    print *, '  -A: Dump the un-decorated abstract syntax tree.'
   ! === (end of automatic option lines) ===
    stop
  end subroutine Option_usage

  subroutine Dump_settings
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
      call output(' Error threshold before halting:                 ', advance='no')
      call blanks(5, advance='no')
      call output(quit_error_threshold, advance='yes')
      call output(' Status on normal exit:                          ', advance='no')
      call blanks(5, advance='no')
      call output(normal_exit_status, advance='yes')
      call output(' Default hdf version for l1b files:              ', advance='no')
      call blanks(5, advance='no')
      call output(level1_hdfversion, advance='yes')
      call output(' Default hdfeos version on reads:                ', advance='no')
      call blanks(5, advance='no')
      call output(default_hdfversion_read, advance='yes')
      call output(' Default hdfeos version on writes:               ', advance='no')
      call blanks(5, advance='no')
      call output(default_hdfversion_write, advance='yes')
      if ( singleChunk /= 0 .and. lastChunk == 0 ) then
      call output(' Compute only the single chunk:                  ', advance='no') 
      call blanks(5, advance='no')                                                   
      call output(singleChunk, advance='yes')
      elseif ( singleChunk /= 0 .and. lastChunk /= 0 ) then
      call output(' Compute chunks in range:                        ', advance='no') 
      call blanks(5, advance='no')                                                   
      call output(singleChunk, advance='no')
      call blanks(5, advance='no')                                                   
      call output(lastChunk, advance='yes')
      endif                      
      call output(' Is this run in forward model parallel?:         ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%fwmParallel, advance='yes')
      call output(' Is this the master task in pvm?:                ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%master, advance='yes')
      if ( parallel%master ) then
      call output(' Master task number:                             ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(parallel%myTid, advance='yes')
      call output(' Command line sent to slaves:                    ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(trim(parallel%pgeName), advance='yes')
      call output(' Command to queue slave tasks:                   ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(trim(parallel%submit), advance='yes')
      call output(' Maximum failures per chunk:                     ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(parallel%maxFailuresPerChunk, advance='yes')
      call output(' Maximum failures per machine:                   ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(parallel%maxFailuresPerMachine, advance='yes')
      endif                      
      call output(' Is this a slave task in pvm?:                   ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%slave, advance='yes')
      if ( parallel%slave ) then
      call output(' Master task number:                             ', advance='no') 
      call blanks(4, advance='no')                                                   
      call output(parallel%masterTid, advance='yes')
      endif                      
      call output(' Preflight check paths?:                         ', advance='no')
      call blanks(4, advance='no')
      call output(checkPaths, advance='yes')
      call output(' Skip all retrievals?:                           ', advance='no')
      call blanks(4, advance='no')
      call output(SKIPRETRIEVAL, advance='yes')
      call output(' Stage in memory instead of a file?:             ', advance='no')
      call blanks(4, advance='no')
      call output(parallel%stageInMemory, advance='yes')
      call output(' Using wall clock instead of cpu time?:          ', advance='no')
      call blanks(4, advance='no')
      call output(use_wall_clock, advance='yes')
      call output(' Number of switches set:                         ', advance='no')
      call blanks(4, advance='no')
      call output(numSwitches, advance='yes')
      if ( switches /= ' ' ) then
        call output(' (All switches)', advance='no')
        call blanks(4, advance='no')
        call output(trim(switches), advance='yes')
      endif
      call output(' Standard output unit:                           ', advance='no')
      call blanks(4, advance='no')
      call output(PrUnit, advance='yes')
      call output(' Log file unit:                                  ', advance='no')
      call blanks(4, advance='no')
      call output(MLSMessageConfig%LogFileUnit, advance='yes')
      call output(' ----------------------------------------------------------', &
        & advance='yes')
    endif
  end subroutine Dump_settings

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, unit_number )
    character(LEN=*), intent(in)   :: Name
    integer, intent(in)            :: unit_number

    call output ( 'Level 2 configuration file ' )
    call output ( 'name : ' )
    call blanks(4)
    call output ( trim(Name), advance='no')
    call blanks(10)
    call output ( 'unit number : ' )
    call blanks(4)
    call output ( unit_number, advance='yes')
  end subroutine announce_success

end program MLSL2

! $Log$
! Revision 2.105  2003/11/07 00:46:51  pwagner
! New quicker preflight option: --checkPaths
!
! Revision 2.104  2003/11/05 21:27:54  pwagner
! Can enter range of chunks to be processed instead of single
!
! Revision 2.103  2003/10/21 00:02:32  pwagner
! Revert to writing just once if toolkit and parallel
!
! Revision 2.102  2003/10/14 18:17:02  pwagner
! Fixed problem with reducing switches to unique list
!
! Revision 2.101  2003/10/09 23:57:35  pwagner
! A few more SIPS-related tweaks
!
! Revision 2.100  2003/10/09 23:35:21  pwagner
! Treats SIPS version special; switches ,-separated
!
! Revision 2.99  2003/09/05 23:22:08  pwagner
! Takes in new --skipRetrieval option
!
! Revision 2.98  2003/08/01 20:26:01  pwagner
! gets slave pge name from command line
!
! Revision 2.97  2003/06/13 20:02:50  vsnyder
! Put snoopMAF before snoop, so it can be found, futzing
!
! Revision 2.96  2003/06/09 22:49:33  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.95  2003/05/14 00:59:40  livesey
! Increased hash table size.
!
! Revision 2.94  2003/05/13 04:48:20  livesey
! Added stgmem option
!
! Revision 2.93  2003/02/27 18:40:29  pwagner
! recl passed to let NAG open l2cf with long lines
!
! Revision 2.92  2003/02/08 00:30:57  pwagner
! Increased default RECL due to new l2cfs in lib
!
! Revision 2.91  2002/12/19 11:54:45  livesey
! Upped tree size
!
! Revision 2.90  2002/12/10 00:39:06  pwagner
! Announces success a la read_apriori
!
! Revision 2.89  2002/12/06 22:33:41  livesey
! Added the snoop name stuff
!
! Revision 2.88  2002/12/05 19:45:20  pwagner
! Moved MLSFile_T from MLSFiles to MLSCommon
!
! Revision 2.87  2002/12/04 01:18:21  pwagner
! First halting steps toward using filedatabase
!
! Revision 2.86  2002/11/13 01:08:40  pwagner
! Prints more info when called asked to print opts
!
! Revision 2.85  2002/10/30 00:56:51  livesey
! Increased size of hash table
!
! Revision 2.84  2002/10/08 17:41:20  livesey
! Now passes chunk argument to slaves if present.  Needed for
! fwmParallel mode.
!
! Revision 2.83  2002/10/05 00:44:14  livesey
! Included the FMWParallel stuff
!
! Revision 2.82  2002/10/03 23:00:03  pwagner
! You can set l1b, l2gp hdfversions on command line
!
! Revision 2.81  2002/09/24 18:17:20  pwagner
! Consistent with add_to_section_timing now calling time_now at its end
!
! Revision 2.80  2002/08/29 21:47:51  livesey
! Added subblock option
!
! Revision 2.79  2002/08/21 18:22:26  vsnyder
! Revise processing of Lahey/Fujitsy runtime (-Wl) options
!
! Revision 2.78  2002/08/15 21:48:51  pwagner
! h5open(close)_f now called at start (end)
!
! Revision 2.77  2002/07/23 23:15:05  pwagner
! Moved dump_settings call after learning name of l2cf file
!
! Revision 2.76  2002/07/18 21:59:28  vsnyder
! Cosmetic changes, move some stuff around so PRUNIT is stdout in init_tables
!
! Revision 2.75  2002/05/28 22:34:59  livesey
! Increased tree size
!
! Revision 2.74  2002/05/24 17:27:18  pwagner
! Added recl=recl to second-chance open of l2cf
!
! Revision 2.73  2002/05/23 20:57:57  vsnyder
! Add --recl option, some cosmetic changes
!
! Revision 2.72  2002/05/14 00:26:49  livesey
! Larger table sizes (may be unnecessary)
!
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
