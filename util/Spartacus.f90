! Copyright 2008, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Spartacus
  use allocate_deallocate, only: allocate_test, deallocate_test
  use highOutput, only: output_date_and_time, timestamp
  use L2Parinfo, only: parallel, initparallel
  use L2Parinfo, only: machine_t, parallel, &
    & petitiontag, giveuptag, grantedtag, notifytag, &
    & sig_finished, sig_register, &
    & sig_hostdied, sig_releasehost, sig_requesthost, sig_thankshost, &
    & machineNameLen, getMachines, getNiceTIDString, &
    & dump, addMachineToDatabase
  use Machine, only: io_error, USleep
  use MLSCommon, only: fileNameLen
  use MLSL2Options, only: current_version_id
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, &
    & MLSMsg_error, &
    & MLSMsg_warning, PVMErrormessage
  use MLSFinds, only: findFirst
  use MLSStrings, only: lowercase, unwraplines
  use output_m, only: blanks, newline, output, outputoptions
  use PrintIt_m, only: Set_Config
  use PVM, only: clearPVMArgs, myPVMSpawn, nextPVMArg, &
    & PVMdatadefault, PVMFinitsend, PVMF90pack, PVMFkill, PVMFmytid, &
    & PVMF90Unpack, &
    & PVMFsend, PVMFnotify, PVMTaskExit, PVMTaskHost
  use time_m, only: time_now, time_config
  use toggles, only: gen, levels, toggle

  ! === (start of toc) ===
  !     c o n t e n t s
  !     - - - - - - - -

  ! Launch materless tasks on hosts requested from l2q queue manager
  ! each task is a line in a supplied text file
  ! === (end of toc) ===
  ! (It is assumed that pvm is already up and running)
  ! Usage:
  ! Spartacus [options] [--] [<] [list]
  ! where list is an ascii file which comes from one of
  ! (i)  a file named on the command line w/o the '<' redirection
  ! (ii) stdin or a file redirected as stdin using '<'
  
  ! For a list of options 
  ! l2q --help
  
  ! Bugs and limitations:
  ! --------------------
  ! (1) The commands in the file list must sensibly handle pwd, io
  !      and any needed environment settings

  implicit none

  character(len=2048) :: command_line ! All the opts
  logical, parameter :: DEEBUG = .false.
  integer :: INFO
  integer :: INUNIT = -1       ! Input unit, * if < 0
  character(len=2048) :: LINE      ! Into which is read the command args
  integer, parameter :: LIST_UNIT = 20  ! Unit # for hosts file if not stdin
  integer :: STATUS                ! From OPEN
  logical :: success               ! From INQUIRE
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T2, T_CONVERSION ! For timing
  integer :: TID                   ! Our own TID

  character(len=*), parameter :: GROUPNAME = "mlsl2"
  integer, parameter          :: AVOIDSELECTEDHOSTSTAG = GIVEUPTAG - 1
  integer, parameter          :: CHECKREVIVEDHOSTSTAG = AVOIDSELECTEDHOSTSTAG - 1
  integer, parameter          :: CHECKSELECTEDHOSTSTAG = CHECKREVIVEDHOSTSTAG - 1
  integer, parameter          :: CLEANMASTERDBTAG = CHECKSELECTEDHOSTSTAG - 1
  integer, parameter          :: DUMPDBTAG = CLEANMASTERDBTAG - 1
  integer, parameter          :: DUMPMASTERSDBTAG = DUMPDBTAG - 1

  integer, parameter          :: DUMPUNIT = LIST_UNIT + 1
!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]

  ! Our data type for the master tasks we'll be communicating with via pvm
  type master_T
    character(len=MACHINENAMELEN) :: name = ' '
    character(len=16            ) :: date = ' '
    integer                       :: tid = 0
    integer                       :: numCmds = 0
    integer                       :: numHosts = 0
    integer                       :: numFreed = 0
    integer, dimension(:), pointer:: hosts => null() ! hosts assigned to this master
    logical                       :: needs_host = .false.
    logical                       :: owes_thanks = .false.
    logical                       :: finished = .false.
  end type master_T

  type options_T
    logical            :: verbose = .false.
    character(len=16)  :: command = 'run'         ! {run, kill, dumphdb, dumpmdb
    logical            :: debug = DEEBUG          !   dump, checkh, clean}
    logical            :: dryrun = .false.
    integer            :: errorLevel = MLSMSG_Warning
    logical            :: exitOnError = .false.   ! quit if an error occurs
    character(len=FILENAMELEN) &
      &                :: cmds_file = '<STDIN>'   ! name of cmds file
    character(len=FILENAMELEN) &
      &                :: dump_file = '<STDIN>'   ! name of dump file
    character(len=FILENAMELEN) &
      &                :: rslts_file = ''  ! name of results file
    character(len=FILENAMELEN) &
      &                :: prefix     = ''  ! to put in front of every cmd
    character(len=FILENAMELEN) &
      &                :: suffix     = ''  ! to put at end of every cmd
    logical :: bufferedDumpFile = .true.          ! lack of buffering slows sips
    logical :: Rescue = .false.                   ! -R option is set
    logical :: Timing = .false.                   ! -T option is set
    logical :: date_and_times = .false.
    character(len=1)   :: timingUnits = 's'       ! l_s, l_m, l_h
    character(len=2048):: selectedHosts =''       ! E.g., 'c0-1,c0-23,c0-55'
    character(len=FILENAMELEN) &
      &                :: PMFile = ''        ! name of periodic Master dump file
    character(len=FILENAMELEN) &
      &                :: PHFile = ''        ! name of periodic host dump file
  end type options_T
  
  integer, parameter :: FIXDELAYFORSLAVETOFINISH   = 1500000 ! 15000000 ! 15 sec
  integer, parameter :: MAXCMDLEN                  = 51200 ! was 4096
  integer, parameter :: MAXLINELEN                 = 512
  integer, parameter :: MAXNUMCMDS                 = 6000
  integer, parameter :: MAXNUMLINES                = 20000
  integer :: BUFFERIDRCV              ! From PVM
  integer, dimension(:), pointer :: CmdMACHINES ! Machine indices for Cmds
  integer, dimension(:), pointer :: CmdTIDS ! Tids for Cmds
  character(len=16), dimension(:), pointer :: CmdNICETIDS ! Tids for Cmds
  integer :: CmdID
  character(len=MAXCMDLEN) :: Cmd
  character(len=MAXCMDLEN), dimension(MAXNUMCMDS) :: Cmds ! Cmds
  logical, dimension(:), pointer :: CmdSCOMPLETED ! Cmds completed
  logical, dimension(:), pointer :: CmdSSTARTED ! Cmds being processed
  logical, dimension(:), pointer :: CmdSABANDONED ! Cmds kept failing
  integer :: L2QTid
  character(len=MAXLINELEN), dimension(MAXNUMLINES) :: lines
  integer :: machNum
  logical :: machineRequestQueued
  type (Machine_T),dimension(:), pointer :: Machines
  integer :: nextCmd
  integer :: noCmds
  integer :: noMachines
  type ( options_T ) :: options
  character(len=MAXCMDLEN), dimension(MAXNUMCMDS) :: rslts ! Result files
  logical :: SKIPDELAY                ! Don't wait before doing the next go round
  logical :: SKIPDEATHWATCH           ! Don't check for deaths
  integer :: TIDARR(1)                ! One tid
  !
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  MLSMessageConfig%CrashOnAnyError = .true.
  time_config%use_wall_clock = .true.
  parallel%master = .true.  ! not merely master, but master of masters
  parallel%slaveFilename = 'pvm' ! for later cures only
  parallel%delay = 10000 ! musec between beats
  !---------------- Task (0) ------------------
  call get_options
  options%verbose = options%verbose .or. options%debug
  select case (options%timingUnits)
  case ('m')
    t_conversion = 1./60
  case ('h')
    t_conversion = 1./3600
  case default
    t_conversion = 1.
  end select
  ! By default we won't to quit if an error occurs
  ! (unless you want us to)
  if ( options%exitOnError ) options%errorLevel = MLSMSG_Error
  ! Do we use buffered output?
  if ( options%dump_file /= '<STDIN>' ) then
    OutputOptions%prunit = DUMPUNIT
    OutputOptions%buffered = options%bufferedDumpFile ! .false.
    ! OutputOptions%opened = .true.
    OutputOptions%name = options%dump_file
    ! print *, 'Opening ', prunit, ' as ', trim(options%dump_file)
    open(OutputOptions%prunit, file=trim(options%dump_file), &
      & status='replace', form='formatted')
  endif
  !---------------- Task (1) ------------------
  ! Is l2q already alive and running?
  call pvmfgettid(GROUPNAME, 0, tid)
  if ( tid < 1 ) then
    call output( 'Sorry--l2q not running', advance='yes' )
    stop
  endif
  !---------------- Task (2) ------------------

  call time_now ( t0 )

  call InitParallel ( 0, 0 )
  ! Learn our own TID
  call pvmfmyTID( tid )
  ! Clear the command line arguments we're going to accumulate to pass
  ! to slave tasks
  call ClearPVMArgs

  !---------------- Task (3) ------------------
  ! Open the List of commands
  ! Store them temporarily as Lines
  ! In case any span multiple lines
  status = 0
  options%cmds_file = '<STDIN>'
  Lines = '(no more)'
  inunit = 0
  line = adjustl(line)
  if ( line /= ' ' ) then
    options%cmds_file = line
    call output( 'commands file name: ' )
    call output( trim(line), advance='yes' )
    call read_textfile_arr( trim(options%cmds_file), Lines )
    call output( trim(Lines(1)), advance='yes' )
    noCmds = FindFirst( Lines == '(no more)' )
    noCmds = noCmds - 1
  else
    ! Loop over commands
    noCmds = 0
    do
      read( *, fmt='(a512)', iostat=status ) cmd
      call output( status, advance='yes' )
      call output( trim(cmd), advance='yes' )
      if ( status /= 0 ) exit
      noCmds = noCmds + 1
      cmd = adjustl(cmd)
      if ( cmd == ' ' .or. cmd(1:1) == '#' ) cycle
      Lines(noCmds) = cmd
    enddo
  endif
  call dump_settings
  call output( ' Before unwrapping lines', advance='yes' )
  do CmdID=1, noCmds
     call output( trim(lines(CmdID)), advance='yes' )
  enddo
  ! Now unsplit any commands that span multiple lines
  call unWrapLines ( Lines(1:noCmds), Cmds, noCmds, escape="\", comment="#" )
  if ( noCmds < 1 ) then
    call output( 'Sorry-no commands found', advance='yes' )
    stop
  endif
  ! Any prefix or suffix required?
  if ( len_trim(options%prefix) > 0 ) then
    do CmdID=1, noCmds
       Cmds(CmdID) = trim(options%prefix) // " " // trim(Cmds(CmdID))
    enddo
  endif
  if ( len_trim(options%suffix) > 0 ) then
    do CmdID=1, noCmds
       Cmds(CmdID) = trim(Cmds(CmdID)) // " " // trim(options%suffix)
    enddo
  endif
  call output( ' After unwrapping lines, adding pre- and suffixes', advance='yes' )
  do CmdID=1, noCmds
    call output( trim(Cmds(CmdID)), advance='yes' )
  enddo
  ! A corresponding results file? 
  ! Unless cmd crashed, a result should be created
  Rslts = ' '
  if ( len_trim(options%rslts_file) > 1 ) then
    call output( 'results file name: ' )
    call output( trim(options%rslts_file), advance='yes' )
    call read_textfile_arr( trim(options%rslts_file), Rslts )
    call output( trim(Rslts(1)), advance='yes' )
  endif
  nullify( CmdsCompleted, CmdsStarted, CmdsAbandoned, CmdMachines, CmdTids, &
    & CmdNiceTids, machines )
  call Allocate_test ( CmdsAbandoned, noCmds, 'CmdsAbandoned', ModuleName )
  call Allocate_test ( CmdsStarted, noCmds, 'CmdsStarted', ModuleName )
  call Allocate_test ( CmdMachines, noCmds, 'CmdMachines', ModuleName )
  call Allocate_test ( CmdTids, noCmds, 'CmdTids', ModuleName )
  call Allocate_test ( CmdNiceTids, noCmds, 'CmdNiceTids', ModuleName )
  call Allocate_test ( CmdsCompleted, noCmds, 'CmdsCompleted', ModuleName )
  
  call RegisterWithL2Q( noCmds, machines, L2QTID )
  if ( options%verbose ) then
    call output( 'registered with l2q ', advance='no' )
    call output( L2QTID, advance='yes' )
  endif
  noMachines = size(machines)
  if ( noMachines < 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
  & 'No machines available for master to assign to slave tasks' )
  if ( .not. any(machines%OK) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
  & 'No machines OK for master to assign to slave tasks' )
  if ( options%verbose ) call dump ( machines )
  ! Loop until all Cmds are done
  CmdsCompleted = options%dryrun   ! .false.
  CmdsStarted = .false.
  CmdsAbandoned = .false.
  CmdTids = 0
  CmdNICETIDS = ' '
  CmdMachines = 0
  machineRequestQueued = .false. ! Request one machine at a time from L2Q
  masterLoop: do ! --------------------------- Master loop -----------------------
    skipDelay = .false.               ! Assume we're going to delay
    skipDeathWatch = .false.          ! Assume we'll listen for slave deaths
    ! This loop is in two main parts.

    ! In the first part, we look to see if there are any Cmds still to be
    ! started done, and any vacant machines to do them.
    ! --------------------------------------------------------- Start new jobs? --
    do while ( CmdAndMachineReady() ) ! (nextCmd, machNum) )
      if ( nextCmd < 1 ) then
        ! Should have returned false
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Illegal Cmd number' )
      endif
      info = SpawnRebellion ( trim(cmds(nextCmd)), PvmTaskHost, &
        & trim(machines(machNum)%Name), 1, tidarr )
      if ( options%verbose ) then
        call output ( 'Tried to spawn ' )
        call TimeStamp ( trim(cmds(nextCmd)), advance='yes' )
        call output ( 'PvmTaskHost ' )
        call output ( PvmTaskHost, advance='yes' )
        call output ( 'on machine ' )
        call output ( trim(machines(machNum)%Name), advance='yes' )
        call output ( 'result was ' )
        call output ( info, advance='yes' )
      end if

      ! Did this launch work
      if ( info == 1 ) then
        machines(machNum)%free = .false.
        machines(machNum)%tid = tidArr(1)
        machines(machNum)%Chunk = nextCmd
        CmdMachines(nextCmd) = machNum
        CmdTids(nextCmd) = tidArr(1)
        CmdNiceTids(nextCmd) = GetNiceTidString(CmdTids(nextCmd))
        CmdsStarted(nextCmd) = .true.
        if ( options%verbose ) then
          if ( options%debug ) then
            call output ( tidArr(1) )
            call output ( ' ' )
          endif
          call output ( 'Launched Cmd ' )
          call output ( nextCmd )
          call TimeStamp ( ' on slave ' // trim(machines(machNum)%name) // &
            & ' ' // trim(CmdNiceTids(nextCmd)), &
            & advance='yes' )
        end if
        call WelcomeSlave ( nextCmd, CmdTids(nextCmd) )
        call ThankL2Q(machines(machNum), L2Qtid)
        skipDeathWatch = .true.
      else
        ! Couldn't start this job, mark this machine as unreliable
        if ( options%verbose ) then
          call output ( 'Unable to start slave task on ' // &
            & trim(machines(machNum)%Name) // ' info=' )
          if ( info == -6 ) then
            call TimeStamp ( '(pvm demon dispirited)', advance='yes' )
          elseif ( info == -7 ) then
            call TimeStamp ( '(file system problem)', advance='yes' )
          elseif ( info < 0 ) then
            call TimeStamp ( info, advance='yes' )
          else
            call TimeStamp ( tidArr(1), advance='yes' )
          end if
          call TimeStamp ( 'Marking this machine as not usable', advance='yes')
        end if
        ! Mark all instances of this machine as not to be used.
        where ( machines%Name == machines(machNum)%Name )
          machines%OK = .false.
        end where
        call ThankL2Q( machines(machNum), L2Qtid )
        ! Send bad news back to l2 queue manager
        call TellL2QMachineDied( machines(machNum), L2Qtid )
      end if
    end do

    if ( options%debug ) &
      & call output( 'Listening for administrative message ', advance='yes' )
    ! ----------------------------------------------------- Administrative messages?
    ! Listen out for any message telling us to quit now
    call PVMFNRecv ( -1, GiveUpTag, bufferIDRcv )
    if ( bufferIDRcv > 0 ) then
      if ( options%verbose ) then
        call TimeStamp ( 'Received an external message to give up, so finishing now', &
          & advance='yes' )
      end if
      exit masterLoop
    end if

    ! Listen out for any message that a slave task has died
    ! ( But only if we're sure there isn't a finished message from it
    !   queued up and waiting for us to read )
    if ( skipDeathWatch ) then
      bufferIDRcv = 0
    else
      call PVMFNRecv ( -1, NotifyTAG, bufferIDRcv )
    endif
    if ( bufferIDRcv > 0 ) then
      ! So we got a message.  There may be another one following on so don't
      ! delay before going round this loop again.
      skipDelay = .true.
      Tid = 0
      CmdID = 0
      machNum = 0
      ! Get the TID for the dead task
      call PVMF90Unpack ( Tid, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking deadTid' )
      ! Now this may well be a legitimate exit, detectable by one of 2 cases
      ! either
      ! (1) we won't know about this tid any more
      ! (2) the machine status was reset to free after a finished signal
      ! Otherwise we need to tidy up.
      machNum = FindFirst ( machines%tid, Tid )
      if ( machNum > 0 ) then
        ! On the other hand
        if ( machines(machNum)%free ) &
          & machNum = 0
      endif
      if ( machNum > 0 ) then
        ! Now, to get round a memory management bug, we'll ignore this
        ! if, as far as we're concerned, the task was finished anyway.
        ! Cmd = machines(machNum)%Chunk
        CmdID = machines(machNum)%Chunk
        machines(machNum)%free = .true.
        if ( options%verbose ) then
          if ( options%debug ) then
            call output ( TID )
            call output ( ' ' )
          endif
          call output ( 'The run of Cmd number ' )
          call output ( CmdID )
          call output ( ' ' )
          call output ( 'on ' // trim(machines(machNum)%Name) // ' ' )
          call TimeStamp ( trim(GetNiceTidString(Tid)) // &
            & ' finished.', advance='yes' )
        end if

        if ( options%verbose ) then
          call output ( 'Got a finished message from ' )
          call output ( trim(machines(machNum)%Name) // ' ' )
          call output ( trim(GetNiceTidString(Tid)) // &
            & ' processing Cmd ' )
          call TimeStamp ( CmdID, advance='yes')
        endif

        ! Now do we have a Result file we should have been expecting?
        if ( len_trim( Rslts(CmdID) ) > 0 ) then
          inquire( file=trim(Rslts(CmdID)), exist=success )
          if ( .not. success ) then
            ! Uh-oh, we crashed .. must tell l2q and try to relaunch task
            call TellL2QMachineDied( machines(machNum), L2Qtid )
            if ( options%verbose ) then
              call output ( 'Bad news about command ' )
              call TimeStamp ( CmdID, advance='yes' )
            end if
            ! Now update our information
            CmdsCompleted(CmdID) = .false.
            CmdsStarted(CmdID) = .false.
            CmdTids(CmdID) = 0
            machines(machNum)%free = .false.
            machines(machNum)%tid = 0
            machines(machNum)%Chunk = 0
            if ( options%verbose ) then
              call printMasterStatus
            end if
            cycle
          endif
        endif
        ! Send news back to l2 queue manager
        call TellL2QMachineFinished( &
          & trim(machines(machNum)%name), machines(machNum)%tid, L2Qtid, &
          & FIXDELAYFORSLAVETOFINISH )
        ! Now update our information
        CmdsCompleted(CmdID) = .true.
        CmdTids(CmdID) = 0
        machines(machNum)%free = .true.
        machines(machNum)%tid = 0
        machines(machNum)%Chunk = 0
        parallel%numCompletedChunks = parallel%numCompletedChunks + 1
        if ( options%verbose ) then
          call printMasterStatus
        end if
      end if

    else if ( bufferIDRcv < 0 ) then
      call PVMErrorMessage ( info, "checking for Notify message" )
    end if

    ! Have we finished all cmds?
    if ( all(CmdSCOMPLETED .or. CmdsAbandoned) ) exit

    ! Now, rather than chew up cpu time on the master machine, we'll wait a
    ! bit here.
    if ( .not. skipDelay .and. parallel%delay > 0 ) &
      & call usleep ( parallel%delay )
  end do masterLoop ! --------------------- End of master loop --------------

  if ( options%verbose ) then
    call TimeStamp ( '   Finished', advance='yes' )
    call dump( CmdNiceTids(1:noCmds), 'CmdNiceTids', options='t' )
    call dump( CmdTids(1:noCmds), 'CmdTids' )
    call dump( machines%tid, 'machines%Tid', format='(i8)' )
  endif
  ! First kill any children still running (only if we got a give up message).
  do CmdID = 1, noCmds
    if ( CmdTids(CmdID) /= 0 ) then
      call usleep ( 100*parallel%delay )
      machNum = CmdMachines(CmdID)
      call TellL2QMachineFinished( &
        & trim(machines(machNum)%name), machines(machNum)%tid, L2Qtid, 0 )
      if ( options%verbose ) then
        call output( 'tid: ', advance='no' )
        call output( machines(machNum)%tid, advance='no' )
        call output( 'machine name: ', advance='no' )
        call output( trim(machines(machNum)%name), advance='no' )
        call TimeStamp ( '   released', &
        & advance='yes' )
      endif
      call usleep ( parallel%delay )
      call pvmfkill ( CmdTids(CmdID), info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'killing slave' )
    end if
  end do

  call TellL2QMasterFinished(L2Qtid)
  if ( options%verbose ) &
    & call TimeStamp ( 'Telling l2q we are finished', &
    & advance='yes' )

  call DeAllocate_test ( CmdsAbandoned, 'CmdsAbandoned', ModuleName )
  call DeAllocate_test ( CmdsStarted,   'CmdsStarted',   ModuleName )
  call DeAllocate_test ( CmdMachines,   'CmdMachines',   ModuleName )
  call DeAllocate_test ( CmdTids,       'CmdTids',       ModuleName )
  call DeAllocate_test ( CmdNiceTids,   'CmdNiceTids',   ModuleName )
  call DeAllocate_test ( CmdsCompleted, 'CmdsCompleted', ModuleName )

contains

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, unit_number, items )
    character(LEN=*), intent(in)             :: Name
    integer, intent(in)                      :: unit_number
    character(LEN=*), optional, intent(in)   :: items

    call output ( 'List of ' )
    if ( present(items) ) then
      call output ( trim(items) )
    else
    call output ( 'hosts' )
    endif
    call output ( ' file name : ' )
    call blanks(4)
    call output ( trim(Name), advance='no')
    call blanks(10)
    call output ( 'unit number : ' )
    call blanks(4)
    call timestamp ( unit_number, advance='yes')
  end subroutine announce_success

  ! ---------------------------------------------  Dump_settings  -----
  subroutine Dump_settings
  ! Show current run-time settings
    call output(' Spartacus called with command line options: ', advance='no')
    call output(trim(command_line), advance='yes')
    call output(' cmds file:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(options%cmds_file), advance='yes')                            
    call output(' io Unit  :', advance='no')  
    call blanks(4, advance='no')                                     
    call output(inunit, advance='yes')                            
    call output(' results file:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(options%rslts_file), advance='yes')                            
    call output(' prefix:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(options%prefix), advance='yes')                            
    call output(' suffix:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(options%suffix), advance='yes')                            
    call output(' -------------- Summary of run time options'      , advance='no')
    call output(' -------------- ', advance='yes')
    call output(' Debug               :                           ', advance='no')
    call blanks(4, advance='no')
    call output(options%debug, advance='yes')
    call output(' dryrun               :                          ', advance='no')
    call blanks(4, advance='no')
    call output(options%dryrun, advance='yes')
    call output(' Verbose             :                           ', advance='no')
    call blanks(4, advance='no')
    call output(options%verbose, advance='yes')
    call output(' Time stamp each report:                         ', advance='no')
    call blanks(4, advance='no')
    call output(options%Timing, advance='yes')
    call output(' Show date and times:                            ', advance='no')
    call blanks(4, advance='no')
    call output(options%date_and_times, advance='yes')
    call output(' Stamp time in what units?:                       ', advance='no')
    call blanks(4, advance='no')
    call output(options%timingUnits, advance='yes' )
    call output(' Using wall clock instead of cpu time?:          ', advance='no')
    call blanks(4, advance='no')
    call output(time_config%use_wall_clock, advance='yes')
    call output(' ----------------------------------------------------------', &
      & advance='yes')
  end subroutine Dump_settings

  ! ---------------------------------------------  get_options  -----
  subroutine get_options
  use machine ! At least HP for command lines, and maybe GETARG, too
    ! Internal variables
    integer :: i
    integer :: j
    integer :: n
    ! Executable
    i = 1+hp
    command_line = ' '
    do ! Process Lahey/Fujitsu run-time options; they begin with "-Wl,"
      call getarg ( i, line )
      if ( line(1:4) /= '-Wl,' ) exit
      i = i + 1
    end do
    do ! Process the command line options to set toggles
      call getarg ( i, line )
      ! call output( 'Processing ' // trim(line), advance='yes' )
      command_line = trim(command_line) // ' ' // trim(line)
      if ( line(1:2) == '--' ) then       ! "word" options
        n = 0
        switch = .true.
        if ( line(3:3) == 'n' .or. line(3:3) == 'N' ) then
          switch = .false.
          n = 1
        end if
        if ( lowerCase(line(3+n:10+n)) == 'help' ) then
          call option_usage
        elseif ( lowerCase(line(3+n:10+n)) == 'dryrun' ) then
          options%dryrun = .true.
        else if ( line(3+n:10+n) == 'version ' ) then
          do j=1, size(current_version_id)
            print*, current_version_id(j)
          enddo
          stop
        else if ( line(3+n:5+n) == 'pre' ) then
          i = i + 1
          call getarg ( i, options%prefix )
        else if ( line(3+n:5+n) == 'suf' ) then
          i = i + 1
          call getarg ( i, options%suffix )
        else if ( line(3+n:7+n) == 'wall ' ) then
          time_config%use_wall_clock = switch
        else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
          exit
        else
          print *, 'unrecognized option ', trim(line), ' ignored.'
          call option_usage
        end if
      else if ( line(1:1) == '-' ) then   ! "letter" options
        j = 1
        j = j + 1
        select case ( line(j:j) )
        case ( ' ' )
          exit
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
        case ( 'r' )
          i = i + 1
          call getarg ( i, line )
          ! call output( 'Processing ' // trim(line), advance='yes' )
          command_line = trim(command_line) // ' ' // trim(line)
          options%rslts_file = line
        case ( 'T' )
          options%timing = .true.
          do
            if ( j >= len(line) ) exit
            if ( line(j+1:j+1) > '0' .and. line(j+1:j+1) <= '9' ) then
              options%date_and_times = .true.
            elseif ( lowercase(line(j+1:j+1)) == 's' ) then
              options%timingUnits = 's'  ! l_seconds
            elseif ( lowercase(line(j+1:j+1)) == 'm' ) then
              options%timingUnits = 'm'
            elseif ( lowercase(line(j+1:j+1)) == 'h' ) then
              options%timingUnits = 'h'
            end if
            j = j + 1
          enddo
        case ( 'v' ); options%verbose = .true.
        case ( 'd' ); options%debug = .true.
        case default
          print *, 'Unrecognized option -', line(j:j), ' ignored.'
          call option_usage
        end select
      else    
        exit ! This must be the hosts filename
      end if
      i = i + 1
    end do

  end subroutine get_options
   
  ! ---------------------------------------------  Option_usage  -----
  ! The following offer some information on the options
  ! available on the command-line
  
  subroutine Option_usage
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [cmds-file]'
    print *, ' Options:'
    print *, ' --dryrun:        show commands that would be executed'
    print *, ' --pre "string":  prefix each cmd with string'
    print *, ' --suf "string":  suffix each cmd with string'
    print *, ' --help:          show help; stop'
    print *, ' --version:       print version string; stop'
    print *, ' --wall:          show times according to wall clock (if T[0] set)'
    print *, ' -d:              debug'
    print *, ' -v:              verbose'
    print *, ' -r rslts-file:   rerun each cmd until corresponding file from'
    print *, '                     rslts-file appears'
    print *, ' -h:              show help; stop'
    print *, ' -T[0][smh]:      show timing [in s, m, h]'
    print *, ' -T1:             show dates with timing'
    stop
  end subroutine Option_usage

  ! ---------------------------------------------  proclaim  -----
  subroutine proclaim( Event, Name, Signal, advance )
    ! Yet another routine to print output
    ! This one sspecialized for signaling noteworthy events
    character(LEN=*), intent(in)   :: Event
    character(LEN=*), intent(in), optional   :: Name
    character(LEN=*), intent(in), optional   :: ADVANCE
    integer, intent(in), optional   :: Signal
    ! Internal variables
    character(len=3) :: myadvance
    ! Executable
    myadvance = 'yes'
    if ( present(advance) ) myAdvance=advance
    if ( options%timing ) call sayTime
    call output(trim(event), advance='no')
    if ( present(name) ) then
      call output(' associated with ', advance='no')
      call output(trim(name), advance='no')
    endif
    if ( present(signal) ) then
      call output(' signal ', advance='no')
      call output(signal, advance='no')
    endif
    if ( lowercase(myadvance) == 'yes' ) call newline
  end subroutine proclaim

  ! ---------------------------------------------  SayTime  -----
  subroutine SayTime ( What )
    character(len=*), intent(in), optional :: What
    ! Internal variables
    ! real :: t_conversion
    ! Executable
    call time_now ( t2 )
    if ( options%date_and_times ) then
      call output_date_and_time(advance='no')
    else
      call output ( "time (now) " )
      call output ( dble(t_conversion*t2), advance = 'no' )
    endif
    call blanks ( 4, advance = 'no' )
    ! call output ( "time (elapsed) " // what // " = " )
    ! call output ( dble(t2 - t1), advance = 'yes' )
    if ( present(What) ) call output ( trim(what), advance = 'yes' )
  end subroutine SayTime

  logical function cmdAndMachineReady()
    ! Return .true. if next Cmd remaining to be done
    ! and a suitable machine can be found
    ! Also set nextCmd, machineName values (which could be made arguments)
    character(len=MACHINENAMELEN) :: machineName
    type (Machine_T)                       :: thisMachine
    nextCmd = 0
    machNum = 0
    machineName = ' '
    ! 1st check if any Cmds remain
    CmdAndMachineReady = .not. all(CmdsStarted .or. CmdsAbandoned)
    if ( options%debug ) then
      call output( 'Command and machine yet ready?', advance='no' )
      call output( CmdAndMachineReady, advance='yes' )
    endif
    if ( .not. CmdAndMachineReady ) return
    nextCmd = FindFirst ( &
      & (.not. CmdsStarted) .and. (.not. CmdsAbandoned) )
    if ( options%debug ) then
      call output( 'NextCmd ', advance='no' )
      call output( NextCmd, advance='no' )
    endif
    CmdAndMachineReady = .not. machineRequestQueued
    if ( .not. CmdAndMachineReady ) then
      if ( options%debug ) call output( 'Hope to receive', advance='yes' )
      call ReceiveMachineFromL2Q(machineName)
    else
      if ( options%debug ) call output( 'Making request', advance='yes' )
      call RequestMachineFromL2Q(L2Qtid, nextCmd, machineName)
    endif
    if ( len_trim(machineName) < 1 ) then
      ! Sorry--no machine ready yet; but our request is/remains queued
      CmdAndMachineReady = .false.
      machineRequestQueued = .true.
      if ( options%debug ) call output('Sorry, l2q scorns us', advance='yes')
    else
      ! Good news--a machine was/has become ready
      CmdAndMachineReady = .true.
      machineRequestQueued = .false.
      ! What if the machine has been added after this master began?
      ! Or, troubles other prevented machine from being added
      machNum = FindFirst( (trim(machineName) == machines%Name) &
        & .and. machines%free )
      if ( machNum < 1 ) then
        thisMachine%name = machineName
        machNum = AddMachineToDataBase(machines, thisMachine)
        if ( options%verbose ) &
          & call output('Added machine to db', advance='yes')
        ! machine = AddMachineNameToDataBase(machines%Name, machineName)
      endif
      if ( options%debug ) call output('Good news, l2q likes us', advance='yes')
    endif
  end function cmdAndMachineReady

  subroutine RegisterWithL2Q(noCmds, machines, L2Qtid)
    integer, intent(in)                                  :: noCmds
    type(machine_t), dimension(:), pointer :: MACHINES
    integer, intent(out)                                 :: L2Qtid
    !
    integer :: BUFFERID                 ! From PVM
    character(len=16) :: DATESTRING
    integer :: INFO                     ! From PVM
    character(len=16) :: L2QSTRING
    !
    call pvmfgettid(GROUPNAME, 0, L2Qtid)
    if ( L2Qtid < 1 ) then
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'l2q queue manager not running--dead or not started yet?' )
    else
      write ( L2QSTRING, * ) L2QTID
      call TimeStamp('Registering with l2q (l2qtid=' // &
        & trim(L2QSTRING) // ')', advance='yes')
    endif
    call GetMachines ( machines )
    ! Register ourselves with the l2 queue manager
    ! Identify ourselves
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_Register, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing registration' )
    call PVMF90Pack ( noCmds, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing number of Cmds' )
    ! Now send our data date as a string
    dateString = '(unknown)'
    call PVMF90Pack ( dateString, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing date' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    ! call PVMF90Pack ( noCmds, info )
    ! if ( info /= 0 ) &
    !   & call PVMErrorMessage ( info, 'packing L2QTid' )
    ! call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
  end subroutine RegisterWithL2Q

  subroutine RequestMachineFromL2Q(L2Qtid, nextCmd, machineName)
    !
    ! request a free host from the l2 queue manager
    integer, intent(in)                               :: L2Qtid
    integer, intent(in)                               :: nextCmd
    character(len=MACHINENAMELEN), intent(out)        :: MACHINENAME
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Identify ourselves
    if ( options%verbose ) &
      & call TimeStamp('Requesting host from l2q', advance='yes')
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( sig_requestHost, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing registration' )
    call PVMF90Pack ( nextCmd, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing number of Cmd' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
    call ReceiveMachineFromL2Q(machineName)
  end subroutine requestMachineFromL2Q

  subroutine ReceiveMachineFromL2Q(machineName)
    !
    ! request a free host from the l2 queue manager
    character(len=MACHINENAMELEN), intent(out)        :: MACHINENAME
    !
    integer :: BEAT
    integer :: INFO                     ! From PVM
    integer :: BUFFERIDRCV              ! From PVM
    ! How many times to listen for queue manager's granting of request
    ! (Setting NBEATS=0 means merely register the request with queue manager;
    ! subsequent trips through master event loop will pick up message
    ! when and if queue manager grants our request)
    integer, parameter :: NBEATS = 0  ! 16   ! ~2x MAXNUMMASTERS  
    ! Possibly a machine is free and the queue manager will
    ! respond quickly
    ! (This idea is tentative, however)
    ! Now let's give queue manager a chance to receive our request
    ! and deliberate a while (but not too long)
    ! before giving up and going back to checking on our own slaves
    
    ! In fact this was a bad idea, or it turned out bad, anyway:
    ! It slowed processing routine communication by the master to no more
    ! than one message per NBEATS*delaytime, so long as the master was
    ! hoping to hear from the queue manager. This meant between 3 and 4 seconds
    ! between each signal being processed, so signals piled up in a huge
    ! line waiting to be received and acted upon.
    machineName = ' '
    do beat = 1, max(NBEATS, 1)
      call PVMFNRecv( -1, GRANTEDTAG, bufferIDRcv )
      if ( bufferIDRcv < 0 ) then
        call PVMErrorMessage ( info, "checking for Granted message" )
      else if ( bufferIDRcv > 0 ) then
        if ( options%verbose ) &
          & call TimeStamp('Weve been granted a host', advance='yes')
        ! Granted us a machine
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking machine name" )
        endif
        if ( len_trim(machineName) < 1 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'granted blank host name' )
        return
      else
        if ( beat < max(NBEATS, 1)) call usleep ( parallel%delay )
      endif
    enddo
    if ( options%verbose .and. len_trim(machineName) > 0 ) &
      & call TimeStamp('Received from l2q ' // trim(machineName), advance='yes')
  end subroutine ReceiveMachineFromL2Q

  subroutine TellL2QMachineDied(machine, L2Qtid)
    type(machine_T), intent(in)                         :: machine
    integer, intent(in)                                 :: L2Qtid
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    !
    ! Identify ourselves
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_HostDied, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing signal of dead machine' )
    ! call PVMF90Pack ( machineName, info )
    call PVMF90Pack ( machine%tid, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing dead machine tid' )
    if ( len_trim(machine%name) < 1 ) then
      call PVMF90Pack ( ' ', info )
    else 
      call PVMF90Pack ( trim(machine%name), info )
    endif
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing dead machine name' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
    if ( options%verbose ) then
      call output('Telling l2q about death of host ', advance='no')
      call output(trim(machine%name), advance='no')
      call blanks(2)
      call TimeStamp(machine%tid, advance='yes')
    endif
  end subroutine TellL2QMachineDied

  subroutine TellL2QMachineFinished(machineName, tid, L2Qtid, delay)
    character(len=*), intent(in)                        :: MACHINENAME
    integer, intent(in)                                 :: tid
    integer, intent(in)                                 :: L2Qtid
    integer, intent(in)                                 :: delay
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    !
    ! If slave finished normally, give it enough time to deconstruct gracefully
    call usleep( delay )
    if ( options%verbose ) then
      call output( 'Releasing ' )
      call output( tid )
      call output( ' ' // trim(machinename), advance='yes' )
    endif
    ! Identify ourselves
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_RELEASEHOST, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing signal of finished machine' )
    call PVMF90Pack ( tid, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing finished tid' )
    call PVMF90Pack ( machineName, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing finished machine name' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
  end subroutine TellL2QMachineFinished

  subroutine TellL2QMasterFinished(L2Qtid)
    integer, intent(in)                                 :: L2Qtid
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    !
    ! Identify ourselves
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_FINISHED, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing signal we finished' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
  end subroutine TellL2QMasterFinished

  subroutine ThankL2Q(machine, L2Qtid)
    type(machine_T), intent(in)         :: machine
    integer, intent(in)                 :: L2Qtid
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    !
    ! Identify ourselves
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_ThanksHost, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing signal of thanks for machine' )
    call PVMF90Pack ( machine%tid, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing tid' )
    call PVMF90Pack ( machine%Name, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing machineName' )
    call PVMF90Pack ( machine%Chunk, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing machineCmd' )
    ! info = INDEX (L2PCF%startUTC, "T")
    ! if ( info > 0 ) machine%master_Date = L2PCF%startUTC(1:info-1)
    call PVMF90Pack ( machine%master_date, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing machineMasterDate' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
    if ( options%verbose ) then
      call output('Thanking l2q for host', advance='no')
      call TimeStamp(machine%tid, advance='yes')
    endif
  end subroutine ThankL2Q

  subroutine printMasterStatus
    ! Print status: completed, underway, abandoned
     call TimeStamp ( 'Master status:', advance='yes' )
     call output ( count(cmdsCompleted) )
     call output ( ' of ' )
     call output ( noCmds )
     call output ( ' commands completed, ')
     call output ( count(cmdsStarted .and. .not. cmdsCompleted) )
     call output ( ' underway, ' )
     call output ( count(cmdsAbandoned) )
     call output ( ' abandoned, ' )
     call output ( count(.not. &
       & (cmdsStarted .or. cmdsCompleted .or. cmdsAbandoned ) ) )
     call TimeStamp ( ' left. ', advance='yes' )
     call output ( count ( .not. machines%Free ) )
     call output ( ' of ' )
     call output ( noMachines )
     call output ( ' machines busy, with ' )
     call output ( count ( .not. machines%OK ) )
     call TimeStamp ( ' being avoided.', advance='yes' )
  end subroutine printMasterStatus

  function spawnRebellion( act, PvmTaskHost, machineName, ntask, tids ) &
    & result(info)
    ! Prepare for myPVMSpawn by separating command from args
    use MLSStringLists, only: numStringElements, stringElement
    ! Args
    character(len=*), intent(in)               :: act
    integer, intent(in)                        :: PvmTaskHost
    character(len=*), intent(in)               :: machineName
    integer, intent(in)                        :: ntask
    integer, dimension(:), intent(out)         :: tids
    integer                                    :: info
    ! Internal variables
    character(len=MAXCMDLEN) :: arg
    integer :: elem
    integer :: nargs
    ! Executable
    call ClearPVMArgs
    nargs = NumStringElements( act, .false., ' ' )
    if ( nargs < 2 ) then
      info = myPVMSpawn ( trim(act), PvmTaskHost, machineName, ntask, tids )
    else
      do elem=2, nargs
        arg = StringElement( act, elem, .false., ' ' )
        call NextPVMArg ( trim(arg) )
      enddo
      arg = StringElement( act, 1, .false., ' ' )
      info = myPVMSpawn ( trim(arg), PvmTaskHost, machineName, ntask, tids )
    endif
    call FreePVMArgs
  end function spawnRebellion

  subroutine WelcomeSlave ( chunk, tid )
    ! This routine welcomes a slave into the fold and tells it stuff
    integer, intent(in) :: CHUNK      ! Chunk we're asking it to do
    integer, intent(in) :: TID        ! Task ID for slave

    ! Local variables
    integer :: info

    ! call SendChunkInfoToSlave ( chunks, chunk, tid )
    ! Now ask to be notified when this task exits
    call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ tid /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'setting up notify' )
  end subroutine WelcomeSlave

  subroutine READ_TEXTFILE_arr ( File, string )
  ! read a textfile into string
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), dimension(:), intent(inout) :: string    ! its contents
    ! Internal variables
    integer :: lun
    integer :: recrd
    integer :: status
    character(len=512) :: tempStr
    ! Try to read the textfile
    write(*,*) 'Now in read_textfile_arr'
    call my_GET_LUN ( LUN )
    write(*,*) 'About to open ', lun
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE-E- Unable to open textfile'
      return
    endif
    recrd = 0
    write(*,*) 'About to read command lines from ', trim(File)
    do
      read( UNIT=lun, fmt='(a512)', IOSTAT=status ) tempStr
      write(*,*) recrd, status, trim(tempStr)
      if ( status /= 0 ) exit
      tempStr = adjustl(tempStr)
      if ( tempStr == ' ' .or. tempStr(1:1) == '#' ) cycle
      recrd = recrd + 1
      string(recrd) = tempStr
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine READ_TEXTFILE_arr

  subroutine my_GET_LUN ( LUN, MSG )
  ! Find a Fortran logical unit number that's not in use.
    integer, intent(out) :: LUN          ! The logical unit number
    logical, intent(in), optional :: MSG ! Print failure message if absent or .true.
    logical :: EXIST, OPENED             ! Used to inquire about the unit
    write(*,*) 'Trying to find available io unit number'
    do lun = 30, 100
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) return
    end do
    write(*,*) 'None found!'
    lun = -1
    if ( present(msg) ) then
      if ( .not. msg ) return
    end if
    write(*,*) 'IO_STUFF%GET_LUN-E- Unable to get a logical unit number'
    return
  end subroutine my_GET_LUN

end program Spartacus

module BOGUS_MODULE
  implicit none
  private
  character (len=*), parameter :: ModuleName= &
       "BOGUS"
contains
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------
end module BOGUS_MODULE

! $Log$
! Revision 1.12  2016/06/06 21:15:56  pwagner
! Gets USleep from machine module; renamed machine counter to machNum
!
! Revision 1.11  2014/03/04 18:50:28  pwagner
! Must allow more commands and lines; should there even be a limit?
!
! Revision 1.10  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.9  2013/08/23 02:51:47  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.8  2013/08/12 23:50:59  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 1.7  2012/02/13 23:47:20  pwagner
! Needed longer MAXCMDLEN; bogus_module needed so idents shows version id
!
! Revision 1.6  2011/06/24 00:40:31  pwagner
! Add -dryrun option; repair some cosmetic errors
!
! Revision 1.5  2011/06/18 00:35:34  pwagner
! Added new --pre and --suf options
!
! Revision 1.4  2009/06/23 19:43:54  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 1.3  2008/10/24 23:16:04  pwagner
! May check each cmd against corresponding line in rslts file for crash
!
! Revision 1.2  2008/04/22 17:59:26  pwagner
! Hope we fixed the bug freeing the same wrong host over and over
!
! Revision 1.1  2008/04/09 17:13:54  pwagner
! First commit
!
