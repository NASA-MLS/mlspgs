! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program l2q
  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test, &
    & Test_Allocate, Test_Deallocate
  use Dates_Module, only: Dateform, Reformatdate
  use Dump_1, only: Dump
  use HighOutput, only: Output_Date_And_Time, OutputNamedValue, TimeStamp
  use L2ParInfo, only: Parallel, InitParallel
  use L2ParInfo, only: Machine_T, Parallel, &
    & Petitiontag, Giveuptag, Grantedtag, Masterdumptag, Notifytag, &
    & Sig_Finished, Sig_Register, Sig_Swearallegiance, Sig_Switchallegiance, &
    & Sig_Hostdied, Sig_Releasehost, Sig_Requesthost, Sig_Thankshost, &
    & Machinenamelen, Getmachinenames, &
    & Dump, AddmachinetoDatabase
  use Machine ! At Least Hp For Command Lines, And Maybe Getarg, Too
  use MLSCommon, only: FileNameLen
  use MLSL2Options, only: Current_Version_Id
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMessageExit, &
    & MLSMSG_Debug, MLSMSG_Error, &
    & MLSMSG_Info, MLSMSG_Success, MLSMSG_Warning, PVMErrorMessage
  use MLSFinds, only: Findfirst
  use MLSStringLists, only: CatLists, GetStringElement, NumStringElements, &
    & StringElementnum
  use MLSStrings, only: Lowercase, ReadIntsFromChars, Streq
  use Output_M, only: Blanks, Newline, &
    & Output, OutputOptions
  use PrintIt_M, only: Set_Config
  use PVM, only: PVMOK, &
    & ClearPVMArgs, Getmachinenamefromtid, &
    & PvmDatadefault, Pvmfinitsend, Pvmf90pack, Pvmfkill, Pvmfmytid, &
    & Pvmf90unpack, Pvmfpstat, &
    & Pvmfsend, Pvmfnotify, Pvmtaskexit, &
    & Pvmffreebuf
  use Sort_M, only: Sort
  use Time_M, only: Time_Now, Time_Config
  use Toggles, only: Gen, Levels, &
    & Toggle

  ! === (start of toc) ===
  !     c o n t e n t s
  !     - - - - - - - -

  ! Main program for queue manager of level 2 mastertasks
  ! === (end of toc) ===
  ! (It is assumed that pvm is already up and running)
  ! Usage:
  ! l2q [options] [--] [<] [list]
  ! where list is an ascii file which comes from one of
  ! (i)  a file named on the command line w/o the '<' redirection
  ! (ii) stdin or a file redirected as stdin using '<'
  
  ! For a list of options 
  ! l2q --help
  
  ! Bugs and limitations:
  ! --------------------
  ! (1) The dump_file, if used, is buffered by default (otherwise everything
  !       is slowed to an unacceptable pace during production)
  ! (2) The hosts named in the file list must match the names returned
  !       by the pvm subroutine PVMFConfig which normally takes its names
  !       from pvm's own input file;
  !       to be clearer, guard against the case where one file names its nodes
  !       "c0-nnn" and the other file something like "jet-0-nnn" because
  !       l2q won't recognize that they refer to the same machines

  implicit none

  integer :: BUFFERID
  character(len=2048) :: command_line ! All the opts
  logical, parameter :: COUNTEMPTY = .true.
  logical, parameter :: DEEBUG = .false.
  logical, parameter :: DUMPDBSONDEBUG = .false.
  integer :: ERROR
  integer :: I
  integer :: INFO
  integer :: INUNIT = -1       ! Input unit, * if < 0
  character(len=2048) :: LINE      ! Into which is read the command args
  integer, parameter :: LIST_UNIT = 20  ! Unit # for hosts file if not stdin
  character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES => null()
  integer, parameter :: MASTERNAMELEN = 16
  character(len=MasterNameLen) :: MASTERSNAME
  integer, parameter :: MAXNUMMASTERS = 100 ! Mas num running simultaneously
  integer :: RECL = 10000          ! Record length for list
  integer :: STATUS                ! From OPEN
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T1, T2, T_CONVERSION ! For timing
  integer :: TAG
  integer :: TID                   ! Our own TID
  character(len=32) :: TIDSTR

  character(len=*), parameter :: GROUPNAME = "mlsl2"
  character(len=*), parameter :: LISTNAMEEXTENSION = ".txt"
  integer, parameter          :: AVOIDSELECTEDHOSTSTAG = GIVEUPTAG - 1
  integer, parameter          :: CHECKREVIVEDHOSTSTAG = AVOIDSELECTEDHOSTSTAG - 1
  integer, parameter          :: CHECKSELECTEDHOSTSTAG = CHECKREVIVEDHOSTSTAG - 1
  integer, parameter          :: CLEANMASTERDBTAG = CHECKSELECTEDHOSTSTAG - 1
  integer, parameter          :: DUMPDBTAG = CLEANMASTERDBTAG - 1
  integer, parameter          :: DUMPMASTERSDBTAG = DUMPDBTAG - 1
  integer, parameter          :: DUMPHOSTSDBTAG = DUMPMASTERSDBTAG - 1
  integer, parameter          :: FREEANYHOSTSTAG = DUMPHOSTSDBTAG - 1
  integer, parameter          :: FREEHOSTSTAG = FREEANYHOSTSTAG - 1
  integer, parameter          :: KILLMASTERSTAG = FREEHOSTSTAG - 1
  integer, parameter          :: SUICIDETAG = KILLMASTERSTAG - 1
  integer, parameter          :: SWITCHDUMPFILETAG = SUICIDETAG - 1
  integer, parameter          :: TELLMASTERDUMPTAG = SWITCHDUMPFILETAG - 1
  integer, parameter          :: TURNREVIVALSONTAG = TELLMASTERDUMPTAG - 1
  integer, parameter          :: TURNREVIVALSOFFTAG = TURNREVIVALSONTAG - 1

  ! These are special tid values
  integer, parameter          :: UNASSIGNED = -1
  integer, parameter          :: AWAITINGREVIVAL = UNASSIGNED - 1
  integer, parameter          :: AWAITINGTHANKS = AWAITINGREVIVAL - 1

  integer, parameter          :: DUMPUNIT = LIST_UNIT + 1
  integer, parameter          :: TEMPUNIT = DUMPUNIT + 1
  integer, parameter          :: HDBUNIT = TEMPUNIT + 1
  integer, parameter          :: MDBUNIT = HDBUNIT + 1
!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  ! Our data type for the master tasks we'll be communicating with via pvm
  type master_T
    character(len=MACHINENAMELEN) :: name = ' '
    character(len=16            ) :: date = ' '
    integer                       :: tid = 0
    integer                       :: numChunks = 0
    integer                       :: numHosts = 0
    integer, dimension(:), pointer:: hosts => null() ! hosts assigned to this master
    logical                       :: needs_host = .false.
    logical                       :: owes_thanks = .false.
    logical                       :: finished = .false.
  end type master_T

!   ! This datatype logs a request for a host
!   type HostRequest_T
!     integer :: Master=0                 ! Which master made the request
!     integer :: CHUNK=0                  ! Which chunk is it for
!     integer :: TICKET=0                 ! What ticket number is it
!     integer :: STATUS=DW_INVALID        ! One of the DW_Status types
!   end type HostRequest_T

  type options_T
    logical            :: checkList = .false.
    logical            :: cleanMasterDB = .false. ! Regularly clean db of done
    logical            :: deferToElders = .true.  ! New masters defer to old
    logical            :: dumpEachNewMaster = .false.
    logical            :: verbose = .false.
    logical            :: freeHosts = .false.     ! Regularly check for abandoned
    logical            :: reviveHosts = .false.   ! Regularly check for revivals
    character(len=16)  :: command = 'run'         ! {run, kill, dumphdb, dumpmdb
    logical            :: debug = DEEBUG          !   dump, checkh, clean}
    integer            :: errorLevel = MLSMSG_Warning
    logical            :: exitOnError = .false.   ! quit if an error occurs
    character(len=FILENAMELEN) &
      &                :: LIST_file = '<STDIN>'   ! name of hosts file
    character(len=FILENAMELEN) &
      &                :: dump_file = '<STDIN>'   ! name of dump file
    character(len=FILENAMELEN) &
      &                :: HDBfile = ''            ! name of hosts DB file
    character(len=FILENAMELEN) &
      &                :: MDBfile = ''            ! name of masters DB file
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
    real               :: MDBPeriod = 600.   ! Period (s) of Master dumps
    real               :: HDBPeriod = 600.   ! Periodic (s) host dumps
  end type options_T
  
  type ( options_T ) :: options
  type(Machine_T), dimension(:), pointer :: hosts => null()
  type(master_T), dimension(:), pointer :: masters => null()
  !
  call set_config ( useToolkit = .false., logFileUnit = -1 )
  MLSMessageConfig%CrashOnAnyError = .true.
  outputOptions%skipmlsmsglogging = .true.
  time_config%use_wall_clock = .true.
  parallel%master = .true.  ! not merely master, but master of masters
  parallel%slaveFilename = 'pvm' ! for later cures only
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
    OutputOptions%name = options%dump_file
    ! print *, 'Opening ', prunit, ' as ', trim(options%dump_file)
    open(OutputOptions%prunit, file=trim(options%dump_file), &
      & status='replace', form='formatted')
  endif
  !---------------- Task (1) ------------------
  ! Is l2q already alive and running?
  call pvmfgettid(GROUPNAME, 0, tid)
  if ( tid > 0 ) then
    ! Have we been ordered to communicate with it?
    if ( options%command == 'inquire' ) then
      call timestamp( 'l2q queue manager already running' , advance='yes')
      call MLSMessageExit 
    elseif ( options%command /= 'run' ) then
      ! Our target will be listening for something like:
      ! call PVMFNRecv ( -1, GiveUpTag, bufferIDRcv )
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMF90Pack ( GROUPNAME, info )
      select case( trim(options%command) )
      case ( 'avoid' )
        tag = avoidSelectedHostsTag
      case ( 'check' )
        tag = checkSelectedHostsTag
      case ( 'checkall' )
        tag = checkRevivedHostsTag
      case ( 'clean' )
        tag = cleanMasterDBTag
      case ( 'dumpdb' )
        tag = dumpDBTag
      case ( 'dumphdb' )
        tag = dumpHostsDBTag
      case ( 'dumpmdb' )
        tag = dumpMastersDBTag
      case ( 'free' )
        tag = freeHostsTag
      case ( 'freeany' )
        tag = freeAnyHostsTag
      case ( 'kill' )
        tag = killMastersTag
      case ( 'revive=on' )
        tag = turnRevivalsOnTag
      case ( 'revive=off' )
        tag = turnRevivalsOffTag
      case ( 'suicide' )
        tag = suicideTag ! GiveUpTag
      case ( 'switch' )
        tag = switchDumpFileTag
      case ( 'tellmdump' )
        tag = tellMasterDumpTag
      case default
        call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'l2q would not recognize unknown command '//trim(options%command) )
      end select
      if ( any(tag == &
        & (/avoidSelectedHostsTag, checkSelectedHostsTag, freeHostsTag, &
        & killMastersTag, tellMasterDumpTag/) )) &
        & then
        call PVMF90Pack ( trim(options%selectedHosts), info )
      else
        call PVMF90Pack ( trim(options%dump_file), info )
      endif
      call PVMFSend ( tid, Tag, info )
      if ( options%timing ) call sayTime
      call output('Already-running l2q tid: ', advance='no')
      call timestamp(tid, advance='yes')
      call timestamp(' commanded: '//trim(options%command), advance='yes')
      if ( any(tag == &
        & (/avoidSelectedHostsTag, checkSelectedHostsTag, freeHostsTag, &
        & killMastersTag, tellMasterDumpTag/) )) &
        & then
        call timestamp(trim(options%selectedHosts), advance='yes')
      elseif ( any(tag == (/dumpHostsDBTag, dumpMastersDBTag, dumpDBTag, &
        & switchDumpFileTag/) )) then
        call timestamp(trim(options%dump_file), advance='yes')
      endif
      call MLSMessageExit 
    else
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'l2q queue manager already running' )
    endif
  endif

  if ( options%command == 'inquire' ) then                                 
    call timestamp( 'l2q queue manager not running' , advance='yes')   
    call MLSMessageExit                                                    
  elseif ( options%command /= 'run' ) then
    call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'l2q not running--ignoring command '//trim(options%command) )
  endif
  !---------------- Task (2) ------------------

  call time_now ( t0 )

  call InitParallel ( 0, 0 )
  call GetMachineNames ( machineNames )
  if ( .not. associated(machineNames) ) &
      & call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'unable to get machine names' )
  if ( size(machineNames) < 1 ) call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'machine names array of zero size' )
  ! Join (actually start) our group
  call PVMFJoinGroup ( GROUPNAME, status )
  if ( status < 0 ) &
    & call myPVMErrorMessage ( status, 'Joining '//GROUPNAME//' group' )
  ! Learn our own TID
  call pvmfmyTID( tid )
! Clear the command line arguments we're going to accumulate to pass
! to slave tasks
  call ClearPVMArgs

  !---------------- Task (3) ------------------
  ! Open the List of hosts
  status = 0
  options%LIST_file = '<STDIN>'
  if ( line /= ' ' ) then
    options%LIST_file = line
    open ( list_unit, file=line, status='old', &
     & form='formatted', access='sequential', recl=recl, iostat=status )
    if ( status /= 0 ) then
      options%LIST_file = trim(line) // LISTNAMEEXTENSION
       open ( list_unit, file=trim(line) // LISTNAMEEXTENSION, status='old', &
        & form='formatted', access='sequential', recl=recl, iostat=status )
    end if
    if ( status /= 0 ) then
      call io_error ( "While opening LIST", status, line )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open LIST file: " // trim(line) )
    elseif( options%verbose ) then                            
      call announce_success(options%LIST_file, list_unit)               
    end if
    inunit = list_unit
  end if
  error = status
  if ( .not. associated(machineNames) ) call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'unable to get machine names after reading list' )
  
  write( tidStr, * ) tid
  call output_date_and_time( msg='starting l2q  TID=' // trim(tidStr) )
  call time_now ( t1 )

  if( options%verbose  ) call dump_settings

  ! Read list of prospective hosts that will run slave tasks for masters
  call read_list
  call cure_host_database( hosts, silent=.true. )
  ! Are we rescuing an older l2q that died?
  if ( options%Rescue ) then
    ! Read Masters db
    open ( MDBUnit, file=options%MDBFile, status='old', &
     & form='formatted', access='sequential', recl=recl, iostat=status )
    if ( status /= 0 ) then
      call io_error ( "While opening MDBFile", status, options%MDBFile )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open MDB file: " // trim(options%MDBFile) )
    elseif( options%verbose ) then                            
      call announce_success( options%MDBFile, MDBUnit, 'rescued masters' ) 
    end if
    call read_master_database( Masters )
    close ( MDBUnit )
    ! Read Hosts db
    open ( HDBUnit, file=options%HDBFile, status='old', &
     & form='formatted', access='sequential', recl=recl, iostat=status )
    if ( status /= 0 ) then
      call io_error ( "While opening HDBFile", status, options%HDBFile )
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open HDB file: " // trim(options%HDBFile) )
    elseif( options%verbose ) then                            
      call announce_success( options%HDBFile, HDBUnit, 'rescued hosts'  )               
    end if
    call read_host_database( Hosts )
    close ( HDBUnit )

    ! This was just a test
    ! call output_date_and_time(msg='ending l2q')
    ! if ( options%dump_file /= '<STDIN>' ) close(OutputOptions%prunit)
    ! stop

    ! call dump_master_database( Masters )
    do i=1, size(masters)
      call SwitchMastersAllegiance( Masters(i) )
    enddo
  end if
  if ( options%checkList ) then
    call dump(hosts)
  !---------------- Task (4) ------------------
  elseif (error == 0) then
    call dump(hosts, details=0)
    if ( .not. associated(machineNames) ) call MLSMessage( MLSMSG_Error, ModuleName, &
    & 'unable to get machine names before event_loop' )
    if ( size(machineNames) < 1 ) call MLSMessage( MLSMSG_Error, ModuleName, &
    & 'machine names array of zero size before event_loop' )
    if( options%verbose .or. DEEBUG ) call timestamp('Entering event loop', advance='yes')
    call event_loop( machineNames )
  else
    call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'error reading list of hosts' )
  end if

  call time_now ( t0 )
  t1 = t0
  if ( options%timing ) call sayTime ( 'Closing and deallocating' )
  call output_date_and_time(msg='ending l2q')
  if ( options%dump_file /= '<STDIN>' ) close(OutputOptions%prunit)

contains

  ! ---------------------------------------------  addMasterToDatabase  -----
  integer function  addMasterToDatabase( DATABASE, ITEM )
    ! This function adds a master data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where the new item is put.

    ! Dummy arguments
    type (master_t), dimension(:), pointer :: DATABASE
    type (master_t), intent(in) :: ITEM

    ! Local variables
    integer :: S
    type (Master_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    addMasterToDatabase = newSize
  end function  addMasterToDatabase

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

  ! ---------------------------------------------  assignHostToMaster  -----
  subroutine assignHostToMaster(host, master, hostsID)
    type(Machine_T) :: host
    type(master_t) :: master
    integer, intent(in) :: hostsID
    ! Internal variables
    integer :: bufferID
    integer :: hcid
    integer :: info
    integer, dimension(:), pointer :: tempHosts
    ! Executable
    ! 1st--check that host is free and OK
    if ( .not. (host%free .and. host%OK) ) then
      call dump(host)
      call dump_master(master)
      call myMLSMessage( options%errorLevel, ModuleName, &
      & 'Tried to assign an unqualified host to this master' )
    elseif( master%owes_thanks ) then
      call myMLSMessage( options%errorLevel, ModuleName, &
      & 'Tried to assign host to an ungrateful master' )
    endif      
    host%master_tid = master%tid
    host%free = .false.
    host%tid = AWAITINGTHANKS
    if ( options%debug ) then
       call output ('Number of machines free: ', advance='no')
       call timestamp (count(hosts%free .and. hosts%OK), advance='yes')
    endif
    hcid = 1
    if ( .not. associated(master%hosts) ) then
      nullify(master%hosts)
      call allocate_test(master%hosts, 1, 'master%hosts', moduleName)
    elseif( size(master%hosts) < 1 ) then
      call allocate_test(master%hosts, 1, 'master%hosts', moduleName)
    else
      ! We'll need a larger hosts component in master, so ..
      nullify(tempHosts)
      hcid = size(master%hosts) + 1
      call allocate_test(tempHosts, hcid, 'temphosts', moduleName)
      tempHosts(1:hcid-1) = master%hosts
      call deallocate_test(master%hosts, 'master%hosts', moduleName)
      master%hosts => tempHosts
    endif
    master%hosts(hcid) = hostsID
    ! Now use PVM to tell master he's got a new slave host
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( trim(host%name), info )
    mastersname = masterNameFun(master%Tid)
    if ( info /= 0 ) &
      & call myPVMErrorMessage ( info, 'packing host name' )
    if ( options%debug ) &
      & call proclaim('Sending new host ' // trim(host%name), &
      & mastersname)
    call PVMFSend ( master%tid, grantedTag, info )
    master%owes_thanks = .true.
    master%needs_host = .false.
    master%numHosts = master%numHosts + 1
  end subroutine assignHostToMaster
  
  ! ---------------------------------------------  clean_master_database  -----
  subroutine clean_master_database(Masters)
    ! clean master database of finished tasks
    type(Master_T), dimension(:), pointer :: Masters
    ! Internal variables
    integer :: i
    ! Executable
    ! call output('Trying to clean master db of finished tasks')
    do i=1, size(Masters)
      if (Masters(i)%finished) then
        call deallocate_test(Masters(i)%hosts, 'hosts', ModuleName)
        ! We can't delete the db entry because
        ! everybody knows masters by their db index
        ! All we can do is reclaim the little bit of the memory allocated
        ! for their hosts component
        ! call rmMasterFromDatabase(Masters, Masters(i))
      endif
    enddo
  end subroutine clean_master_database

  ! ---------------------------------------------  cure_host_database  -----
  subroutine cure_host_database(hosts, silent, selectedHosts)
    ! cure any recovered hosts in database--declare them useable ("OK")
    ! If optional selectedHosts is present, cure only them
    type(Machine_T), intent(inout), dimension(:) :: hosts
    logical, optional, intent(in) :: silent
    character(len=*), optional, intent(in) :: selectedHosts
    ! Internal variables
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES => null()
    logical :: mySilent
    integer :: i
    integer :: status
    ! Executable
    mySilent = .not. options%verbose
    if ( present(silent) ) mySilent = silent .or. mySilent
    ! No need to cure if none have died
    if ( .not. mySilent ) then
      call outputNamedValue( 'number hosts alive', count(hosts%OK) )
      call outputNamedValue( 'number hosts dead', count(.not. hosts%OK) )
    endif
    if ( all(hosts%OK) ) return    
    if ( .not. mySilent .and. present(selectedHosts) ) &
      & call dump( selectedHosts, 'selectedHosts' )
    call GetMachineNames ( machineNames )
    if ( .not. associated(machineNames) ) &
        & call myMLSMessage( MLSMSG_Error, ModuleName, &
        & 'unable to get machine names' )
    if ( size(machineNames) < 1 ) call myMLSMessage( MLSMSG_Error, ModuleName, &
        & 'machine names array of zero size' )
    if ( .not. mySilent .and. DEEBUG ) &
      & call dump( machineNames, 'machineNames' )
    do i=1, size(hosts)
      if ( hosts(i)%OK ) cycle
      ! Warning--following use of variable 'status' is non-standard
      ! status is good if non-zero, bad if 0
      status = 1
      if ( present(selectedHosts) ) &
        & status = StringElementNum(selectedHosts, trim(hosts(i)%name), .true.)
      if ( status < 1 ) cycle
      ! hosts(i)%OK = any( hosts(i)%name == machineNames )
      hosts(i)%OK = any( streq(hosts(i)%name, machineNames, '-f') )
      if ( hosts(i)%OK .and. .not. mySilent ) &
        & call proclaim('Host Cured', hosts(i)%Name)
      if ( .not. hosts(i)%OK .and. .not. mySilent ) &
        & call proclaim('Sorry--host still not useable', hosts(i)%Name)
      if ( .not. hosts(i)%OK .and. DEEBUG ) then
        call output( hosts(i)%name, advance='yes' )
        call dump( machineNames, 'machinenames' )
      endif
    enddo
    call deAllocate_test(machineNames, 'machineNames', moduleName )
  end subroutine cure_host_database

  ! ---------------------------------------------  dump_his_hosts  -----
  subroutine dump_his_hosts( hostIDs, sortUs )
    integer, dimension(:), intent(in)   :: hostIDs
    logical, intent(in), optional   :: sortUs
    ! Internal variables
    integer :: i
    integer, dimension(size(hostIDs)) :: myHostIDs
    character(len=MACHINENAMELEN), dimension(size(hostIDs)) :: myHostNames
    logical :: mySort
    integer :: maxNameLength
    ! Executable
    mySort = .false.
    if ( present(sortUs) ) mySort=sortUs
    myHostIDs = hostIDs
    if ( mySort ) call sort(myHostIDs, 1, size(hostIDs))
    maxNameLength = 1
    do i=1, size(hostIDs)
      myHostNames(i) = hosts(myHostIDs(i))%name
      maxNameLength = max(maxNameLength, len_trim(myHostNames(i)) + 1)
    enddo
    maxNameLength = min(maxNameLength, MACHINENAMELEN)
    call dump( myHostNames(:)(1:maxNameLength) )
  end subroutine dump_his_hosts

  ! ---------------------------------------------  dump_master_database  -----
  subroutine dump_master_database(Masters, Details)
    ! dump data type
    type(Master_T), dimension(:), pointer :: Masters
    integer, intent(in), optional :: DETAILS  ! 1: dump each master, too
    ! Internal variables
    integer :: i
    integer :: myDetails
    ! Executable
    if ( .not. associated(Masters) ) then
      call output('Master database not associated', advance='yes')
      return
    endif
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ('Size of Master database: ', advance='no')
    call output (size(Masters), advance='yes')
    call output ('Number of Masters active: ', advance='no')
    call output (count(.not. Masters%finished), advance='yes')
    call output ('Number of Masters waiting for hosts: ', advance='no')
    call output (count(Masters%needs_host), advance='yes')
    if ( myDetails < 1 ) return
    do i = 1, size(Masters)
      call dump_Master(Masters(i))
    enddo
  end subroutine dump_Master_database

  ! ---------------------------------------------  dump_Master  -----
  subroutine dump_Master(Master)
    ! dump data type
    type(Master_T), intent(in) :: Master
    call output('machine name: ', advance = 'no')
    call output(trim(Master%name), advance = 'yes')
    call output('date: ', advance = 'no')
    call output(trim(Master%date), advance = 'yes')
    call output('tid: ', advance = 'no')
    call output(Master%tid, advance = 'yes')
    call output('needs host?: ', advance = 'no')
    call output(Master%needs_host, advance = 'yes')
    call output('owes thanks?: ', advance = 'no')
    call output(Master%owes_thanks, advance = 'yes')
    call output('finished?: ', advance = 'no')
    call output(Master%finished, advance = 'yes')
  end subroutine dump_Master
   
  ! ---------------------------------------------  Dump_settings  -----
  subroutine Dump_settings
  ! Show current run-time settings
    call output(' l2q called with command line options: ', advance='no')
    call output(trim(command_line), advance='yes')
    call output(' LIST file:', advance='no')  
    call blanks(4, advance='no')                                     
    call output(trim(options%LIST_file), advance='yes')                            
    call output(' -------------- Summary of run time options'      , advance='no')
    call output(' -------------- ', advance='yes')
    call output(' Just check the host list:                       ', advance='no')
    call blanks(4, advance='no')
    call output(options%checkList, advance='yes')
    call output(' Debug               :                           ', advance='no')
    call blanks(4, advance='no')
    call output(options%debug, advance='yes')
    call output(' Verbose             :                           ', advance='no')
    call blanks(4, advance='no')
    call output(options%verbose, advance='yes')
    call output(' Dump each new master:                           ', advance='no')
    call blanks(4, advance='no')
    call output(options%dumpEachNewMaster, advance='yes')
    call output(' Regularly check for revived hosts:              ', advance='no')
    call blanks(4, advance='no')
    call output(options%reviveHosts, advance='yes')
    call output(' Regularly check and free abandoned hosts:       ', advance='no')
    call blanks(4, advance='no')
    call output(options%freeHosts, advance='yes')
    call output(' Regularly clean db of finished masters:         ', advance='no')
    call blanks(4, advance='no')
    call output(options%cleanMasterDB, advance='yes')
    call output(' New masters defer to old:                       ', advance='no')
    call blanks(4, advance='no')
    call output(options%deferToElders, advance='yes')
    if ( options%PMFile /= '') then
    call output(' Periodically dump Master DB to:                 ', advance='no')
    call blanks(4, advance='no')
    call output(trim(options%PMFile), advance='yes')
    call output(' Period (s):                                     ', advance='no')
    call blanks(4, advance='no')
    call output(options%MDBPeriod, advance='yes')
    endif
    if ( options%PHFile /= '') then
    call output(' Periodically dump host DB to:                   ', advance='no')
    call blanks(4, advance='no')
    call output(trim(options%PHFile), advance='yes')
    call output(' Period (s):                                     ', advance='no')
    call blanks(4, advance='no')
    call output(options%HDBPeriod, advance='yes')
    endif
    call output(' Dump other output to:                           ', advance='no')
    call blanks(4, advance='no')
    call output(trim(options%dump_File), advance='yes')
    call output(' Using unit number:                              ', advance='no')
    call blanks(4, advance='no')
    call output(OutputOptions%prunit, advance='yes')
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
    call output(' Using buffered output?:                         ', advance='no')
    call blanks(4, advance='no')
    call output(Options%bufferedDumpFile, advance='yes')
    call output(' ----------------------------------------------------------', &
      & advance='yes')
  end subroutine Dump_settings

  ! ---------------------------------------------  event_loop  -----
  subroutine event_loop( machineNames )
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    ! The main event loop
    ! l2q spends eternity circling endlessly, responding
    ! to the occasional promptings from master tasks via pvm
    ! or else initiating a contact with one
    ! or responding to commands sent by a fraternal l2q launched with
    ! "-c command" option
    integer :: bufferIDRcv
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    logical :: CHECKREVIVEDHOSTS
    logical :: CLEANMASTERDB
    logical :: DUMPHOSTS
    logical :: DUMPMASTERS
    logical :: FREEABANDONEDHOSTS
    integer :: grandMastersID           ! index into database of an older master
    integer :: host
    integer :: hostsID                  ! index into database of a host
    integer :: i
    integer, dimension(MAXNUMMASTERS) :: IDs
    integer :: INFO                     ! From PVM
    character(len=MachineNameLen)  :: MACHINENAME
    integer :: mastersID                ! index into database of a master
    integer :: masterTid                ! TID of a master
    logical :: mayAssignAHost
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: nextFree
    integer :: numHosts
    integer :: numMasters
    logical :: opened
    integer :: SIGNAL                   ! From a master
    character(len=16) :: SIGNALSTRING
    logical :: SIGNIFICANTEVENT
    logical :: SKIPDELAY            ! Don't wait before doing the next go round
    integer :: STATUS
    character(len=FileNameLen)  :: tempfile
    integer :: TID
    real    :: tLastHDBDump
    real    :: tLastMDBDump
    type(master_t) :: aMaster
    ! Executable
    if ( .not. associated(machineNames) ) call MLSMessage( MLSMSG_Error, ModuleName, &
    & 'unable to get machine names in event_loop' )
    if ( size(machineNames) < 1 ) call MLSMessage( MLSMSG_Error, ModuleName, &
    & 'machine names array of zero size in event_loop' )
    numMasters = 0
    tLastHDBDump = 0.
    tLastMDBDump = 0.
    ! call timeStamp( 'Looping over events', advance='yes' )
    do
      checkrevivedhosts = .false.
      dumpHosts = .false.
      dumpMasters = .false.
      freeabandonedhosts = .false.
      mayAssignAHost = .true.
      significantEvent = .false.
      skipDelay = .false.               ! Assume we're going to delay
      ! Listen out for pvm events
      call PVMFNRecv( -1, PETITIONTAG, bufferIDRcv )
      if ( bufferIDRcv < 0 ) then
        call myPVMErrorMessage ( info, "checking for Petition message" )
      else if ( bufferIDRcv > 0 ) then
        ! So we got a message.  There may be another one following on so don't 
        ! delay before going round this loop again.
        skipDelay = .true.
        ! Who sent this?
        call PVMFBufInfo ( bufferIDRcv, bytes, msgTag, masterTid, info )
        if ( info /= 0 ) &
          & call myPVMErrorMessage ( info, "calling PVMFBufInfo" )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) then
          call myPVMErrorMessage ( info, "unpacking signal" )
        endif
        ! if ( options%debug ) call proclaim('Master sent message', signal=signal)
        select case (signal) 

        case ( sig_register ) ! ----------------- Master registering ------
          call PVMF90Unpack ( aMaster%NumChunks, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking number of chunks" )
          endif
          call PVMF90Unpack ( aMaster%date, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking date string" )
          endif
          aMaster%tid = masterTid
          call GetMachineNameFromTid ( masterTid, aMaster%Name, info )
          if ( info == -1 ) & 
            & call myMLSMessage( options%errorLevel, ModuleName, &
            & 'Unable to get machine name from tid' )
          aMaster%NumHosts = 0
          aMaster%needs_Host = .false.
          aMaster%owes_thanks = .false.
          aMaster%finished = .false.
          numMasters = addMasterToDatabase(masters, aMaster)
          significantEvent = .true.
          ! masterName = catlists('m-', numMasters)
          ! Now ask to be notified when this master exits
          call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ masterTid /), info )
          if ( info /= 0 ) call myPVMErrorMessage ( info, 'setting up notify' )
          if ( options%verbose ) then
            ! call proclaim('New Master', aMaster%Name, signal=MasterTid, advance='no')
            mastersname = masterNameFun(masterTid)
            call proclaim('New Master', mastersname, advance='no')
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
          if ( options%dumpEachNewMaster ) call dump_master(aMaster)
          dumpMasters = .true.
        case ( sig_requestHost ) ! ----------------- Request for host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          grandMastersID = FindFirst(masters%needs_host .and. &
            & .not. masters%owes_thanks)
          if ( options%deferToElders ) &
            grandMastersID = FindFirst(.not. masters%finished .and. &
            &  masters%needs_host )
          if ( grandMastersID < 1 ) grandMastersID = size(masters) + 100
          mastersname = masterNameFun(masterTid)
          if ( options%debug ) &
            & call proclaim('Master requested host', mastersname)
          ! Any free hosts?
          nextFree = FindFirst( hosts%free .and. hosts%OK )
          if ( masters(mastersID)%owes_thanks ) then
            ! Don't assign hosts to ungrateful masters
            mastersname = masterNameFun(masterTid)
            if ( options%verbose ) &
              & call proclaim('Ungrateful master requested another host', &
              & mastersname, advance='yes')
            masters(mastersID)%needs_host = .true.
          elseif ( nextFree == 0 ) then
            ! No hosts free; master put into needs_host state
            masters(mastersID)%needs_host = .true.
          elseif ( grandMastersID < mastersID ) then
            ! An older master has priority; master put into needs_host state
            masters(mastersID)%needs_host = .true.
          else
            ! Assign free host to master
            call assignHostToMaster( hosts(nextFree), masters(mastersID), nextFree )
            mayAssignAHost = .false. ! Already did one this cycle
          endif
        case ( sig_releaseHost ) ! ----------------- Done with this host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          call IDHostFromMaster(sig_releaseHost, tid, hostsID)
          if ( mastersID < 1 ) then
            call myMLSMessage( MLSMSG_Warning, ModuleName, &
              & 'Unknown master attempted to free a host' )
          elseif ( hostsID > 0 ) then
            call releaseHostFromMaster( hosts(hostsID), masters(mastersID), &
              & hostsID )
            if ( options%debug ) then
              mastersname = masterNameFun(masterTid)
              call proclaim('Master ' // trim(mastersname) // &
                & ' freed host', hosts(hostsID)%name, advance='no')
              call output(' (Still has ', advance='no')
              call output(masters(mastersID)%numHosts, advance='no')
              call timestamp(' remaining)', advance='yes')
              call output ('Number of machines free: ', advance='no')
              call timestamp (count(hosts%free .and. hosts%OK), advance='yes')
            endif
            mayAssignAHost = ( mastersID < size(masters) ) ! .false. ! Don't assign to a later master
          else
            call timestamp ( 'Unknown host freed by master', advance='yes')
          endif

          ! In case the master's thanks was somehow lost
          if ( mastersID > 0 ) masters(mastersID)%owes_thanks = .false.
          dumpHosts = .true.
        case ( sig_ThanksHost ) ! -------------- Happy with this host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          masters(mastersID)%owes_thanks = .false.
          ! if ( options%debug ) &
          !   & call proclaim('Master thanks for host', signal=masterTid)
          call PVMF90Unpack ( tid, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking tid" )
          endif
          call PVMF90Unpack ( machineName, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking machineName" )
          endif
          hostsID = FindFirst(  hosts%name==machineName .and. &
            & hosts%tid == AWAITINGTHANKS .and. hosts%master_tid == masterTid )
          if ( hostsID < 1 ) then
            call output(' mastersID ', advance='no')
            call output(mastersID , advance='yes')
            call output(' machineName ', advance='no')
            call output(trim(machineName) , advance='yes')
            call myMLSMessage( options%errorLevel, ModuleName, &
               & 'Master thanked us for unknown host' )
          endif
          hosts(hostsID)%tid = tid
          call PVMF90Unpack ( hosts(hostsID)%chunk, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking machine%chunk" )
          endif
          call PVMF90Unpack ( hosts(hostsID)%master_date, info )
          if ( info /= 0 ) then
            call myPVMErrorMessage ( info, "unpacking machine%master_date" )
          endif
          ! masters(mastersID)%date = hosts(hostsID)%master_date
          hosts(hostsID)%master_name = masterNameFun(masterTID)
          if ( options%verbose ) then
            mastersname = masterNameFun(masterTid)
            call proclaim('Master ' // trim(mastersname) // &
            & ' thanks for host', hosts(hostsID)%Name, advance='no')
            call output(' ; tid ', advance='no')
            call output(tid, advance='no')
            call output(' (Now has ', advance='no')
            call output(masters(mastersID)%numHosts, advance='no')
            call timestamp(' )', advance='yes')
          endif
          dumpHosts = .true.
        case ( sig_HostDied ) ! -------------- Done away with this host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          mastersname = masterNameFun(masterTid)
          if ( options%debug ) &
            & call proclaim('Master ' // trim(mastersname) // &
            & ' reports host died', advance='no')
          call IDHostFromMaster(sig_HostDied, tid, hostsID)
          if ( hostsID > 0 ) then
            call releaseHostFromMaster( hosts(hostsID), masters(mastersID), &
              & hostsID )
            ! Mark this host and any others sharing the same node
            ! hosts(hostsID)%OK = .false.
            call mark_host_database(hosts, trim(hosts(hostsID)%name))
            significantEvent = .true.
            if ( options%debug ) then
              call output( ' Host ' // trim(hosts(hostsID)%Name), advance='no')
              call output(' (Now has ', advance='no')
              call output(masters(mastersID)%numHosts, advance='no')
              call timestamp(' )', advance='yes')
            elseif ( options%verbose ) then
              call proclaim('Host Died', hosts(hostsID)%Name)
            endif
          endif
          dumpHosts = .true.
          ! In case we were waiting for grateful acknowledgement
          masters(mastersID)%owes_thanks = .false.
        case ( sig_Finished ) ! ----------------- Done with this master ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          numHosts = 0
          if ( associated(masters(mastersID)%hosts) ) &
            & numHosts = size(masters(mastersID)%hosts) ! masters(mastersID)%numHosts
          if ( numHosts > 0 ) then
            do host = 1, numHosts
              hostsID = masters(mastersID)%hosts(host)
              if ( hostsID > 0 .and. hostsID <= size(hosts) ) &
                & call releaseHostFromMaster( &
                & hosts(hostsID), masters(mastersID), hostsID )
            enddo
          endif
          masters(mastersID)%finished = .true.
          masters(mastersID)%needs_host = .false.
          masters(mastersID)%owes_thanks = .false.
          ! call deAllocate_test(masters(mastersID)%hosts, 'hosts', moduleName )
          significantEvent = .true.
          if ( options%verbose ) then
            mastersname = masterNameFun(masterTid)
            call proclaim('Master ' // trim(mastersname) // &
              & ' Finished', advance='no')
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
          dumpMasters = .true.
        case default
          write(SIGNALSTRING, *) signal
          mastersname = masterNameFun(masterTid)
          if ( options%debug ) &
            & call proclaim( 'Unrecognized signal from Master ' // &
            & mastersname )
          call myMLSMessage( options%errorLevel, ModuleName, &
             & 'Unrecognized signal from master' )
        end select
        call PVMFFreeBuf ( bufferIDRcv, info )
        if ( options%debug .and. info /= 0 ) &
          & call timestamp('Trouble freeing buf', advance='yes')
        if ( info /= 0 ) &
          & call myPVMErrorMessage ( info, 'freeing receive buffer' )
      endif
      ! ----------------------------------------------------- Master Task Died?
      ! Listen out for any message telling us one of our masters died
      call PVMFNRecv ( -1, NotifyTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        ! Get the TID for the dead task
        call PVMF90Unpack ( masterTid, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking dead masterTid' )
        ! Now this may well be a legitimate exit, 
        ! in which case, master finished normally.  
        mastersID = FindFirst(masters%tid, masterTid)
        if ( mastersID < 1 ) call myMLSMessage( options%errorLevel, ModuleName, &
           & 'Exit signal from unrecognized master' )
        if ( .not. masters(mastersID)%finished ) then
          ! Otherwise we need to tidy up.
          numHosts = 0
          if ( associated(masters(mastersID)%hosts) ) &
            & numHosts = size(masters(mastersID)%hosts) ! masters(mastersID)%numHosts
          if ( numHosts > 0 ) then
            do host = 1, numHosts
              hostsID = masters(mastersID)%hosts(host)
              if ( hostsID > 0 .and. hostsID <= size(hosts) ) &
                & call releaseHostFromMaster( &
                & hosts(hostsID), masters(mastersID), hostsID )
            enddo
          endif
          masters(mastersID)%finished = .true.
          masters(mastersID)%needs_host = .false.
          masters(mastersID)%owes_thanks = .false.
          significantEvent = .true.
          if ( options%verbose ) then
            mastersname = masterNameFun(masterTid)
            call proclaim('Master Died', mastersname)
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
        endif
        dumpMasters = .true.
      endif
      ! ----------------------------------------------------- Administrative messages?
      ! Listen out for any message telling *us* to quit now
      call PVMFNRecv ( -1, suicideTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        if ( options%timing ) call sayTime
        call timestamp ( 'Received an external message to give up, so finishing now', &
          & advance='yes' )
        exit ! event_loop
      end if

      ! Listen out for any message telling us to stop regular revivals
      call PVMFNRecv ( -1, turnRevivalsOffTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call timestamp ( 'Received an external message to stop regular revivals', &
          & advance='yes' )
        options%reviveHosts = .false.
      end if
      ! Listen out for any message telling us to resume regular revivals
      call PVMFNRecv ( -1, turnRevivalsOnTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call timestamp ( 'Received an external message to resume regular revivals', &
          & advance='yes' )
        options%reviveHosts = .true.
      end if

      cleanMasterDB = significantEvent .and. options%cleanMasterDB
      ! Listen out for any message telling us to clean master database
      call PVMFNRecv ( -1, CleanMasterDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) &
        & call timestamp ( 'Received an external message to clean masterDB', &
          & advance='yes' )
      if ( bufferIDRcv > 0 .or. cleanMasterDB ) then
        if ( options%debug .and. bufferIDRcv < 1 ) &
          & call timestamp('Cleaning database of finished masters', advance='yes')
        call clean_master_database(Masters)
      endif

      checkrevivedhosts  = significantEvent .and. options%reviveHosts
      freeabandonedhosts = significantEvent .and. options%freeHosts
      ! Listen out for any message telling us to dump databases
      call PVMFNRecv ( -1, DumpDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking dumpfile' )
        call redumpHosts( hosts, tempFile, tempUnit )
        call redumpMasters( masters, tempFile, tempUnit )
        dumpHosts = .true.
        dumpMasters = .true.
      end if
      call PVMFNRecv ( -1, DumpHostsDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking dumpfile' )
        call redumpHosts( hosts, tempFile, tempUnit )
        dumpHosts = .true.
      end if
      call PVMFNRecv ( -1, DumpMastersDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking dumpfile' )
        call redumpMasters( masters, tempFile, tempUnit )
        dumpMasters = .true.
      end if
      ! Listen out for any message telling us to kill selected masters
      call PVMFNRecv ( -1, killMastersTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to kill masters ' // &
          & trim(options%selectedHosts), advance='yes' )
        call IDMasterFromName( masters, options%selectedHosts, IDs )
        do i = 1, count( IDs > 0 )
          if ( Ids(i) < 1 .or. Ids(i) > size(masters) ) then
            call timestamp ( 'DB lacks any master matching ' // &
              & trim(options%selectedHosts), advance='yes' )
          else
            call timestamp ( 'killing master ' // &
              & trim(masters(Ids(i))%name), advance='yes' )
            if ( options%verbose ) call dump_Master( masters(Ids(i)) )
            call PVMFSend ( masters(Ids(i))%tid, giveUpTag, info )
            masters(Ids(i))%needs_host = .false. ! So we don't assign it hosts
          endif
        enddo
        dumpMasters = .true.
      end if
      ! Listen out for any message telling us to avoid selected hosts
      call PVMFNRecv ( -1, avoidSelectedHostsTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to avoid ' // &
          & trim(options%selectedHosts), advance='yes' )
        call mark_host_database(hosts, options%selectedHosts)
        dumpHosts = .true.
      end if
      ! Listen out for any message telling us to free hosts
      ! belonging to finished masters
      call PVMFNRecv ( -1, freeAnyHostsTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call timestamp ( 'Received an external message to free hosts' &
        &  // ' belonging to finished masters', advance='yes' )
      endif
      if ( bufferIDRcv > 0 .or. freeabandonedhosts ) then
        call free_hosts( hosts, machineNames, masters=masters )
        dumpHosts = .true.
      endif
      ! Listen out for any message telling us to free specified hosts 
      call PVMFNRecv ( -1, freeHostsTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to free ' // &
          & trim(options%selectedHosts), advance='yes' )
        call free_hosts( hosts, machineNames )
        dumpHosts = .true.
      endif
      ! Listen out for any message telling us to check for revived hosts
      call PVMFNRecv ( -1, checkRevivedHostsTag, bufferIDRcv )
      ! Do it if so commanded, or if part of regular drill
      if ( bufferIDRcv > 0 ) call timestamp ( &
        & 'Received an external message to check for revivals', advance='yes' )
      if ( bufferIDRcv > 0 .or. checkrevivedhosts ) then
        call cure_host_database( hosts )
        dumpHosts = .true.
      endif
      ! Listen out for any message telling us to revive selected hosts
      call PVMFNRecv ( -1, checkSelectedHostsTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to check ' // &
          & trim(options%selectedHosts), advance='yes' )
        call cure_host_database(hosts, .false., options%selectedHosts)
        dumpHosts = .true.
      end if
      ! Listen out for any message telling us to flush current dumpfile
      ! and switch future output to new one
      call PVMFNRecv ( -1, switchDumpFileTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        if ( OutputOptions%prunit == DUMPUNIT ) then
          if ( options%timing ) call sayTime
          tempfile = options%dump_file
          call output('Switching dump file from '//trim(options%dump_file))
          call output(' to ')
          call PVMF90Unpack ( options%dump_file, info )
          if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking newdumpfile' )
          call timestamp(trim(options%dump_file), advance='yes')
          close(DUMPUNIT)
          ! print *, 'Opening ', prunit, ' as ', trim(options%dump_file)
          open( OutputOptions%prunit, file=trim(options%dump_file), &
            & status='replace', form='formatted', iostat=status )
          if ( status /= 0 ) then
            call myMLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'No new dumpfile; unable to open ' // &
              & trim(options%dump_file) )
            options%dump_file = tempfile
            open( OutputOptions%prunit, file=trim(options%dump_file), &
              & status='replace', form='formatted', iostat=status )
            if ( status /= 0 ) call myMLSMessage( options%errorLevel, ModuleName, &
              & 'Now we cant even revert to old dumpfile; unable to open ' // &
              & trim(options%dump_file) )
          else
            call timestamp('(new dump file opened)', advance='yes')
          endif
        else
          call timestamp('(switchDumpFileTag ignored--dumping to <STDIN>)', &
            & advance='yes')
        end if
      end if
      ! Listen out for any message telling us to tell selected masters to dump
      call PVMFNRecv ( -1, tellMasterDumpTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call myPVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to tell masters to dump ' // &
          & trim(options%selectedHosts), advance='yes' )
        call IDMasterFromName( masters, options%selectedHosts, IDs )
        do i = 1, count( IDs > 0 )
          if ( Ids(i) < 1 .or. Ids(i) > size(masters) ) then
            call timestamp ( 'DB lacks any master matching ' // &
              & trim(options%selectedHosts), advance='yes' )
          else
            call timestamp ( 'telling master ' // &
              & trim(masters(Ids(i))%name) // ' to dump', advance='yes' )
            if ( options%verbose ) call dump_Master( masters(Ids(i)) )
            call PVMFSend ( masters(Ids(i))%tid, masterDumpTag, info )
          endif
        enddo
        dumpMasters = .true.
      end if
      ! Unless there's a good reason not to do so,
      ! check if any hosts are currently free and any masters need hosts
      if ( mayAssignAHost .and. associated(masters) ) then
        ! Find oldest master needing a host and grateful for the ones he has
        mastersID = FindFirst( masters%needs_host .and. &
          & .not. masters%owes_thanks )
        grandMastersID = mastersID
        if ( options%deferToElders ) &
          grandMastersID = FindFirst(.not. masters%finished .and. &
            &  masters%needs_host )
        if ( grandMastersID < 1 ) grandMastersID = size(masters) + 100
        if ( mastersID > 0 .and. mastersID == grandMastersID ) then
          ! Find first free host
          nextFree = FindFirst( hosts%free .and. hosts%OK )
          if ( nextFree > 0 ) then
            call assignHostToMaster( hosts(nextFree), masters(mastersID), &
              & nextFree )
            skipDelay = .true.
          endif
        endif
      endif
      ! Now perform any other housekeeping, administration, etc. tasks
      ! Perhaps a periodic dump of hosts
      if ( options%PHFile /= ' ' ) then
        call time_now ( t2 )
        if ( (t2-tLastHDBDump > options%HDBPeriod) .or. dumpHosts ) then
          call reDumpHosts ( hosts, options%PHFile, tempUnit )
          tLastHDBDump = t2
          if ( options%debug .and. DUMPDBSONDEBUG ) then
            inquire(unit=OutputOptions%prunit, opened=opened)
            if ( .not. opened ) call myMLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Miscellaneous dump file not opened ' // &
              & trim(options%dump_file) )
            call dump(hosts)
          endif
        endif
      endif
      ! Perhaps a periodic dump of masters
      if ( options%PMFile /= ' ' ) then
        call time_now ( t2 )
        if ( (t2-tLastMDBDump > options%MDBPeriod) .or. dumpMasters ) then
          call reDumpMasters( masters, options%PMFile, tempUnit )
          tLastMDBDump = t2
          if ( options%debug .and. DUMPDBSONDEBUG ) then
            inquire(unit=OutputOptions%prunit, opened=opened)
            if ( .not. opened ) call myMLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Miscellaneous dump file not opened ' // &
              & trim(options%dump_file) )
            call dump_master_database(masters)
          endif
        endif
      endif
      ! Now, rather than chew up cpu time on this machine, we'll wait a
      ! bit here.
      if ( .not. skipDelay .and. parallel%delay > 0 ) then
        call usleep ( parallel%delay )
      endif
    enddo
  end subroutine event_loop
  
  ! ---------------------------------------------  free_hosts  -----
  subroutine free_hosts(hosts, machineNames, masters, silent)
    ! free any hosts in database 
    ! (1) specified in options%selectedHosts; or
    ! (2) belonging to finished masters
    type(Machine_T), intent(inout), dimension(:)       :: hosts
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    type(Master_T), optional, intent(in), dimension(:) :: masters
    logical, optional, intent(in) :: silent
    ! Internal variables
    integer :: i
    integer :: masterID
    logical :: mySilent
    integer :: selected
    ! Executable
    mySilent = .not. options%verbose
    if ( present(silent) ) mySilent = silent .or. mySilent
    if ( size(hosts) < 1 ) return
    ! No need to free if all are
    if ( all(hosts%free) ) return
    if ( .not. present(masters) ) then
      ! usage (1) 
      if ( .not. associated(machineNames) ) call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'unable to get machine names in free_hosts' )
      if ( size(machineNames) < 1 ) call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'machine names array of zero size in free_hosts' )
      if ( DEEBUG ) call timeStamp('selected hosts to free: ' // trim(options%selectedHosts), advance='yes' )
      do i=1, size(hosts)
        if ( hosts(i)%free ) cycle
        if ( DEEBUG ) call timeStamp('checking: '// trim(hosts(i)%name) // ' ' )
        selected = StringElementNum(options%selectedHosts, trim(hosts(i)%name), .true.)
        if ( DEEBUG ) call timeStamp(selected, advance='yes' )
        if ( selected < 1 ) cycle
        ! hosts(i)%free = any( hosts(i)%name == machineNames )
        hosts(i)%free = any( streq(hosts(i)%name, machineNames, '-f') )
        if ( hosts(i)%free .and. .not. mySilent ) &
          & call proclaim('Host freed', hosts(i)%Name)
        if ( .not. hosts(i)%free .and. DEEBUG ) then
          call output( hosts(i)%name, advance='yes' )
          call dump( machineNames, 'machinenames' )
        endif
      enddo
      return
    endif
    ! usage (2) 
    if ( size(masters) < 1 ) return
    ! 
    if ( .not. any(masters%finished) ) return
    do i=1, size(hosts)
      if ( hosts(i)%free ) cycle
      masterID = findFirst( hosts(i)%master_tid == masters%tid )
      if ( masterID > 0 ) &
        & hosts(i)%free = masters(masterID)%finished
    enddo
  end subroutine free_hosts

  ! ---------------------------------------------  get_options  -----
  subroutine get_options
    ! Internal variables
    integer :: i
    integer :: j
    integer :: lineVal
    integer :: n
    character(len=8) :: subcase
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
        if ( line(3+n:5+n) == 'buf' ) then
          options%bufferedDumpFile = switch
        elseif ( line(3+n:7+n) == 'clean ' ) then
          options%cleanMasterDB = switch
        elseif ( line(3+n:7+n) == 'check ' ) then
          options%checklist = switch
        elseif ( line(3+n:7+n) == 'defer' ) then
          options%deferToElders = switch
        elseif ( line(3+n:5+n) == 'dhp' ) then
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
          call readIntsFromChars(line, lineVal)
          options%HDBPeriod = lineVal
        elseif ( line(3+n:5+n) == 'dmp' ) then
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
          call readIntsFromChars(line, lineVal)
          options%MDBPeriod = lineVal
        elseif ( lowerCase(line(3+n:10+n)) == 'dumpnewm' ) then
          options%dumpEachNewMaster = switch
        elseif ( lowerCase(line(3+n:10+n)) == 'help' ) then
          call option_usage
        elseif ( lowerCase(line(3+n:10+n)) == 'inquire' ) then
          options%command = 'inquire'
        else if ( line(3+n:10+n) == 'version ' ) then
          do j=1, size(current_version_id)
            print*, current_version_id(j)
          enddo
          stop
        else if ( line(3+n:6+n) == 'free' ) then
          options%freeHosts = switch
        else if ( line(3+n:7+n) == 'reviv' ) then
          options%reviveHosts = switch
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
        case ( 'k' ); options%command = 'suicide' ! options%killer = .true.
        case ( 'R' ); options%Rescue = .true.
        case ( 'c' )
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
          options%command = lowercase(line)
          ! Special commands take yet another arg
          if ( index('check,avoid,kill,free,tellmdump', trim(options%command)) > 0 ) then
            i = i + 1
            call getarg ( i, line )
            command_line = trim(command_line) // ' ' // trim(line)
            options%selectedHosts = lowercase(line)
          endif
        case ( 'o' )
          subcase = line(j:j+2)
          i = i + 1
          call getarg ( i, line )
          ! print *, 'subcase: ', subcase
          ! print *, 'file: ', trim(line)
          command_line = trim(command_line) // ' ' // trim(line)
          if ( subcase == 'o') then
          options%dump_file = line
          elseif ( subcase == 'oph') then
          options%PHFile = line
          elseif ( subcase == 'opm') then
          options%PMFile = line
          else
          call option_usage
          endif
        case ( 'r' )
          subcase = line(j:j+3)
          i = i + 1
          call getarg ( i, line )
          ! print *, 'subcase: ', subcase
          ! print *, 'file: ', trim(line)
          command_line = trim(command_line) // ' ' // trim(line)
          if ( subcase == 'r') then
            ! We have no provision for a bare "r"
          call option_usage
          elseif ( subcase == 'rhdb') then
          options%HDBFile = line
          call output('Setting rescue HDB to: ' // trim(options%HDBFile), advance='yes' )
          elseif ( subcase == 'rmdb') then
          options%MDBFile = line
          call output('Setting rescue MDB to: ' // trim(options%MDBFile), advance='yes' )
          else
          call option_usage
          endif
        case ( 's' )
          options%command = 'switch'
          i = i + 1
          call getarg ( i, line )
          command_line = trim(command_line) // ' ' // trim(line)
          options%dump_file = line
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
   
  ! ---------------------------------------------  IDHostFromMaster  -----
  subroutine IDHostFromMaster(signal, tid, hostsID)
    ! Figure out which hosts ID a master is talking about
    ! Arguments
    integer, intent(in)  :: signal
    integer, intent(out) :: tid
    integer, intent(out) :: hostsID
    ! Internal variables
    integer :: info
    character(len=MACHINENAMELEN) :: machineName
    ! Executable
    hostsID = 0
    call PVMF90Unpack ( tid, info )
    if ( info /= 0 ) then
      call myPVMErrorMessage ( info, "unpacking tid" )
    endif
    call PVMF90Unpack ( machineName, info )
    if ( info /= 0 ) then
      call myPVMErrorMessage ( info, "unpacking machineName" )
    endif
    if ( options%debug  ) then
       if ( signal == SIG_RELEASEHOST ) then
         call output ('Freed host ', advance='no')
       elseif ( signal == sig_HostDied ) then
         call output ('Host died ', advance='no')
       endif
       call output ('tid: ', advance='no')
       call output (tid, advance='no')
       call output ('   machineName: ', advance='no')
       call output (trim(machineName), advance='yes')
    endif
    if ( tid < 1 ) then
      ! This means master couldn't start slave task on host
      ! So we'll look at machineName
      hostsID = FindFirst(hosts%name, machineName)
    else
      hostsID = FindFirst(hosts%tid, tid)
    endif
    if ( hostsID < 1 ) then
      ! Oops--the tid doesn't match any of the ones we've had so far
      if ( options%verbose ) then
        if ( options%timing ) call sayTime
        if ( signal == sig_HostDied ) then
          call output('Warning--unknown tid died', advance='no')
        else
          call output('Warning--unknown tid freed', advance='no')
        endif
        call output('  tid ', advance='no')
        call output(tid, advance='no')
        call output('  machineName ', advance='no')
        call output(trim(machineName), advance='no')
        call output('  hostsID ', advance='no')
        call timestamp(hostsID, advance='yes')
      endif
      if ( tid > 0 ) then
        ! if ( options%verbose ) call dump (hosts%tid, format='(i10)')
        ! We'll try again--using the machine name
        hostsID = FindFirst(hosts%name, machineName)
        if ( hostsID < 1 ) call myMLSMessage( MLSMSG_Warning, ModuleName, &
          & 'Even the machine name ' // trim(machineName) // &
          & ' yields an ID outside list of hosts' )
      else
        call myMLSMessage( MLSMSG_Warning, ModuleName, &
          & 'machine name yields an ID outside list of hosts' )
      endif
    endif
  end subroutine IDHostFromMaster
  
  ! ---------------------------------------------  IDMasterFromName  -----
  subroutine IDMasterFromName ( masters, names, IDs )
    ! Find which masterIDs match given list of names
    type(master_t), dimension(:)       :: masters
    character(*), intent(in)           :: names
    integer, dimension(:), intent(out) :: IDs
    ! Internal variables
    logical :: byName
    integer :: i, ivalue, n
    character(len=16) :: inputFormat
    character(len=16) :: masterFormat
    character(len=16) :: masDate
    character(len=16) :: masName
    ! Executable
    if ( options%debug ) then
    call output(' len_trim(names) ', advance='no')
    call output( len_trim(names) , advance='no')
    call output(' size(IDs) ', advance='no')
    call output( size(IDs) , advance='no')
    call output('   size(masters) ', advance='no')
    call output( size(masters) , advance='yes')
    endif
    if ( size(IDs) < 1 ) return
    IDs = 0
    if ( size(masters) < 1 .or. len_trim(names) < 1 ) return
    masterFormat = dateForm(masters(1)%date)
    n = NumStringElements( names,  countEmpty )
    n = min( n, size(IDs) )
    if ( options%debug ) then
    call timestamp ( 'n: ', advance='no' )
    call timestamp ( n, advance='yes' )
    endif
    do i=1, n
      call GetStringElement( names, masName, i, countEmpty )
      ! Two possible formats:
      ! master name (possibly prefixed by 'm-')
      ! We assume that we'll never have more than 10^6 masters in our db
      ! (a limitation, true, but unreasonable)
      ! and that date formats require at least 8 characters
      ! (which is captured by the next logical)
      byName = (len_trim(masName) < 8)
      if ( byName ) then
        if ( index(lowerCase(masName), 'm-') > 0 ) then
          ! e.g., 'm-1'
          call readIntsFromChars ( trim(masName(3:)), ivalue )
        else
          !  e.g., '1'
          call readIntsFromChars ( trim(masName), ivalue )
        endif
        if ( ivalue < 1 ) then
          ! Already zero
        elseif ( ivalue > size(masters) ) then
          ! Already zero
        elseif ( masters(ivalue)%finished ) then
          ! Already zero
        else
          IDs(i) = ivalue
        endif
      else
      ! master date; e.g. '2006-171'
        inputFormat = dateForm(masName)
        masDate = reFormatDate(masName, &
          & fromForm=trim(inputFormat), toForm=masterFormat)
        Ids(i) = findFirst( (masters%date == trim(masDate)) .and. &
          & .not. masters%finished )
        call output(' (inputFormat) ', advance='no')
        call output( trim(inputFormat) , advance='no')
        call output('   (masterFormat) ', advance='no')
        call output( trim(masterFormat) , advance='yes')
        call output(' (input) ', advance='no')
        call output( trim(masName) , advance='no')
        call output('   (reformatted) ', advance='no')
        call output( trim(masDate) , advance='yes')
      endif
      call timestamp ( 'ID: ', advance='no' )
      call timestamp ( Ids(i), advance='yes' )
    enddo
  end subroutine IDMasterFromName

  ! ---------------------------------------------  mark_host_database  -----
  subroutine mark_host_database(hosts, selectedHosts)
    ! mark selected hosts in database as not to be used
    type(Machine_T), intent(inout), dimension(:) :: hosts
    character(len=*), intent(in) :: selectedHosts
    ! Internal variables
    integer :: i
    integer :: status
    ! Executable
    ! No need to mark if none can be used
    if ( all(.not. hosts%OK) ) return
    do i=1, size(hosts)
      if ( .not. hosts(i)%OK ) cycle
      ! Warning--following use of variable 'status' is non-standard
      ! status is good if non-zero, bad if 0
      status = StringElementNum(selectedHosts, trim(hosts(i)%name), .true.)
      if ( status > 0 ) hosts(i)%OK = .false.
    enddo
  end subroutine mark_host_database

  ! ---------------------------------------------  masterNameFun  -----
  function masterNameFun ( masterTid ) result(name)
    character(LEN=MachineNameLen)  :: Name
    integer, intent(in)            :: masterTid
    ! Internal variables
    integer :: mastersID
    !
    mastersID = FindFirst(masters%tid, masterTid)
    name = catLists('m', mastersID, inseparator='-')
  end function masterNameFun

  ! ---------------------------------------------  myMLSMessage  -----
  subroutine myMLSMessage ( Severity, ModuleNameIn, Message )
    ! A 'safer' substitute for MLSMessage
    ! (in case we don't trust the toolkit)
    integer, intent(in) :: Severity ! e.g. MLSMSG_Error
    character (len=*), intent(in) :: ModuleNameIn ! Name of module (see below)
    character (len=*), intent(in) :: Message ! Line of text

    select case ( severity )
    case ( MLSMSG_Success )
      call timeStamp('Success: ' // trim(Message), advance='yes')
    case ( MLSMSG_Debug )
      call timeStamp('Debug: ' // trim(Message), advance='yes')
    case ( MLSMSG_Info )
      call timeStamp('Info: ' // trim(Message), advance='yes')
    case ( MLSMSG_Warning )
      call timeStamp('Warning: ' // trim(Message), advance='yes')
    case ( MLSMSG_Error )
      call MLSMessage ( Severity, ModuleNameIn, Message )
    case default
      call timeStamp('Default: ' // trim(Message), advance='yes')
    end select
  end subroutine myMLSMessage

  ! ---------------------------------------------  myPVMErrorMessage  -----
  subroutine myPVMErrorMessage ( info, place )
    ! This routine is called to log a PVM error
    ! possibly because we may not wish to exit on such an error
    integer, intent(in) :: INFO
    character (LEN=*) :: PLACE
    character (LEN=132) :: LINE

    write (line, * ) info
    if ( options%exitOnError ) then
      call PVMErrorMessage ( info, place )
    else
      call myMLSMessage ( options%errorLevel, Place, &
        & 'PVM Error:  Info='//trim(adjustl(line)))
    endif
  end subroutine myPVMErrorMessage

  ! ---------------------------------------------  Option_usage  -----
  ! The following offer some information on the options
  ! available on the command-line
  
  subroutine Option_usage
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [LIST-name]'
    print *, ' Options:'
    print *, ' --[n]buf:    do [not] buffer dumpfile (if used)'
    print *, ' --check:     check LIST of hosts'
    print *, ' --[n]clean:  do [not] regularly clean db of finished masters'
    print *, ' --defer:     new masters defer to oldest (has enough hosts)'
    print *, ' --dhp n :    dump host db every n seconds (need phfile)'
    print *, ' --dmp m :    dump master db every m seconds (need pmfile)'
    print *, ' --dumpnewm:  dump each new master'
    print *, ' --help:      show help; stop'
    print *, ' --inquire:   show whether an instance of l2q already running'
    print *, ' --[n]free:   do [not] regularly check, free abandoned hosts'
    print *, ' --[n]revive: do [not] regularly check for revived hosts'
    print *, ' --version:   print version string; stop'
    print *, ' --wall:      show times according to wall clock (if T[0] set)'
    print *, ' -k:          kill currently-running l2q'
    print *, ' -oph phfile: dump db of hosts periodically to phfile'
    print *, ' -opm pmfile: dump db of masters periodically to pmfile'
    print *, ' -o dumpfile: direct most other output to dumpfile'
    print *, '               (if modifying an already-running l2q, direct only'
    print *, '               the requested output to dumpfile)'
    print *, ' -s newfile:  flush current dumpfile; direct later output'
    print *, '                to newfile ( same as -o newfile -c switch )'
    print *, ' -v:          verbose'
    print *, ' -h:          show help; stop'
    print *, ' -R:          Rescue: take over from a dead l2q'
    print *, ' -rhdb file:  read db of hosts from file (required by a rescue)'
    print *, ' -rmdb file:  read db of masters from file (required by a rescue)'
    print *, ' -T[0][smh]:  show timing [in s, m, h]'
    print *, ' -T1:         show dates with timing'
    print *, ' -c command [arg]:  issue command to currently-running l2q'
    print *, '        command may be one of'
    print *, ' avoid h1[,h2,..]:'
    print *, '              mark selected host[s] as unuseable'
    print *, ' check h1[,h2,..]:'
    print *, '              check for revival of selected host[s]'
    print *, ' checkall:    check for all revived hosts'
    print *, ' clean:       clean the db of finished master tasks'
    print *, ' freeany:     free any hosts held by finished master tasks'
    print *, ' free h1[,h2,..]:'
    print *, '              free specified hosts'
    print *, ' revive=on:   begin regularly checking for revived hosts'
    print *, ' revive=off:  cease regularly checking for revived hosts'
    print *, ' kill m1[,m2,..]:'
    print *, '              kill selected masters cleanly'
    print *, ' tellmdump m1[,m2,..]:'
    print *, '              tell selected masters to dump current status'
    print *, ' switch:      flush current dumpfile'
    print *, '              ( switch to one specified by -o newfile )'
    print *, '  ( if the following commands are accompanied by the option'
    print *, '    -o newfile, the output requested will be dumped to newfile )'
    print *, ' dumphdb:     dump host database'
    print *, ' dumpmdb:     dump master database'
    print *, ' dumpdb:      dump host and master databases'
    print *, ' N o t e :'
    print *, ' only the options:'
    print *, ' -o, -v, -T, --check, --clean, --defer, '
    print *, ' --dumpnewm, --[n]revive, --wall'
    print *, ' are appropriate for l2q when first launched'
    print *, ' all other options will modify an already-running l2q'
    print *, ' (and will generate an error if l2q not already running)'
    print *, ' '
    print *, ' commanding an already-running l2q to dump output to a file'
    print *, ' without specifying a path will select that working directory'
    print *, ' in force when the original l2q was launched '
    print *, ' (and not necessarily the current one)'
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

  ! ---------------------------------------------  read_host_database  -----
  subroutine read_host_database( Hosts )
    type(machine_T), dimension(:), pointer :: Hosts
    ! Internal variables
    integer :: i
    character(len=32) :: irrelevant
    integer :: hostSize
    integer :: status
    ! Executable
    deallocate( hosts, stat=status )
    if ( status /= 0 ) call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'Unable to deallocate Hosts before reading db' )
    hostSize = 0
    ! 1st--we need to go past the 7 lines at the start of the file
    if ( DEEBUG ) call timeStamp('1st--we need to go past the 7 lines at the start of the file', advance='yes')
    call read_oneline ( HDBUnit, status, string=irrelevant )
    if ( DEEBUG ) then
    call timeStamp( 'Read irrelevant: ' // trim(irrelevant), advance='yes' )
    call output( 'status after read_oneline' )
    call blanks(4)
    call timeStamp( status, advance='yes' )
    endif
    if ( status /= 0 ) then
      call output('Hosts database still empty', advance='yes')
      return
    endif
    call read_oneline ( HDBUnit, status, int=hostSize )
    if ( DEEBUG )  then
    call output( 'host size: ' )
    call output( hostsize, advance='yes' )
    endif
    do i=1, 5
      call read_oneline ( HDBUnit, status, string=irrelevant )
    enddo
    do i=1, hostSize
      call read_host(Hosts)
    enddo
    if ( DEEBUG ) call dump( Hosts, Details=1 )
  end subroutine read_host_database
  
  ! ---------------------------------------------  read_host  -----
  subroutine read_host( hosts )
    type(machine_T), dimension(:), pointer :: hosts
    ! Internal variables
    type(machine_t) :: ahost
    integer :: numhosts
    integer :: status
    ! Executable
    call read_oneline ( HDBUnit, status, ahost%Name )
    if ( DEEBUG ) then
    call timeStamp( 'Read aHost%name: ' // trim(aHost%name) )
    call output( 'status after read_oneline' )
    call blanks(4)
    call timeStamp( status, advance='yes' )
    endif
    if ( status == 0 ) call read_oneline ( HDBUnit, status, int=ahost%master_tid )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, string=ahost%master_name )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, string=ahost%master_date )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, int=ahost%tid )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, int=ahost%chunk )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, log=ahost%OK )
    if ( status == 0 ) call read_oneline ( HDBUnit, status, log=ahost%free )
    if ( status == 0 ) numhosts = addMachineToDatabase(hosts, ahost)
  end subroutine read_host

  ! ---------------------------------------------  read_list  -----
  subroutine read_list
    type(Machine_T) :: newhost
    integer :: iostat
    integer :: nhosts
    nhosts = 0
    do
      if ( inunit >= 0 ) then
        read ( inunit, '(a)', advance='no', eor=100, end=200, err=400, &
          iostat=iostat ) newhost%name
      else
        read ( *, '(a)', advance='no', eor=100, end=200, err=400, &
          iostat=iostat ) newhost%name
      end if
100   newhost%name = adjustl( newhost%name )
      ! Ignore commented-out lines and comments
      if ( newhost%name(1:1) == '#' ) cycle
      if ( options%verbose ) &
          & call output(trim(newhost%name) // ' added', advance='yes')
      newhost%master_tid = -1
      ! newhost%tid = -1
      newhost%tid = UNASSIGNED
      newhost%chunk = -1
      newhost%OK = .false.
      newhost%free = .true.
      nhosts = addMachineToDatabase(hosts, newhost)
    enddo
200 if ( options%verbose ) then
      call output(nhosts, advance='no')
      call timestamp(' hosts read', advance='yes')
    endif
    return
400 call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'error in reading list of hosts' )
  end subroutine read_list

  ! ---------------------------------------------  read_master_database  -----
  subroutine read_master_database( Masters )
    type(Master_T), dimension(:), pointer :: Masters
    ! Internal variables
    integer :: i
    character(len=32) :: irrelevant
    integer :: masterSize
    integer :: status
    ! Executable
    masterSize = 0
    ! 1st--we need to go past the 4 lines at the start of the file
    if ( DEEBUG ) call timeStamp('1st--we need to go past the 4 lines at the start of the file', advance='yes')
    call read_oneline ( MDBUnit, status, string=irrelevant )
    if ( DEEBUG ) then
    call timeStamp( 'Read irrelevant: ' // trim(irrelevant), advance='yes' )
    call output( 'status after read_oneline' )
    call blanks(4)
    call timeStamp( status, advance='yes' )
    endif
    if ( status /= 0 ) then
      call output('Master database still empty', advance='yes')
      return
    endif
    call read_oneline ( MDBUnit, status, int=masterSize )
    call read_oneline ( MDBUnit, status, string=irrelevant )
    call read_oneline ( MDBUnit, status, string=irrelevant )
    if ( DEEBUG ) then
    call output( 'size of master db' )
    call blanks(4)
    call timeStamp( masterSize, advance='yes' )
    endif
    do i=1, masterSize
      call read_master(Masters)
    enddo
    if ( DEEBUG ) call dump_master_database( Masters, Details=1 )
  end subroutine read_master_database
  
  ! ---------------------------------------------  read_master  -----
  subroutine read_master( Masters )
    type(Master_T), dimension(:), pointer :: Masters
    ! Internal variables
    type(master_t) :: aMaster
    integer :: numMasters
    integer :: status
    ! Executable
    call read_oneline ( MDBUnit, status, aMaster%Name )
    if ( DEEBUG ) then
    call timeStamp( 'Read aMaster%name: ' // trim(aMaster%name) )
    call output( 'status after read_oneline' )
    call blanks(4)
    call timeStamp( status, advance='yes' )
    endif
    if ( status == 0 ) call read_oneline ( MDBUnit, status, string=amaster%Date )
    if ( status == 0 ) call read_oneline ( MDBUnit, status, int=amaster%tid )
    if ( status == 0 ) call read_oneline ( MDBUnit, status, log=amaster%needs_host )
    if ( status == 0 ) call read_oneline ( MDBUnit, status, log=amaster%owes_thanks )
    if ( status == 0 ) call read_oneline ( MDBUnit, status, log=amaster%finished )
    if ( status == 0 ) numMasters = addMasterToDatabase(masters, aMaster)
  end subroutine read_master

  ! ---------------------------------------------  read_oneline  -----
  subroutine read_oneline ( unit, status, string, int, log )
    ! Read one line from unit assuming it was formatted as
    ! Something: value
    ! returning value either as a string or as an int
    !
    ! status will be 0 unless EOR, ERROR, or somesuch
    integer, intent(in)                     :: unit
    integer, intent(out)                    :: status
    character(len=*), optional, intent(out) :: string
    integer, optional, intent(out)          :: int
    logical, optional, intent(out)          :: log
    ! Internal variables
    integer :: iColon
    character(len=256) :: last
    character(len=256) :: line
    ! logical, parameter :: DEEBUG = .false.
    ! Executable
    last = ' '
    line = ' '
    ! do
    read ( unit, '(a)', advance='no', eor=100, end=200, err=400, &
      & iostat=status ) line
    if ( DEEBUG ) call output( 'oneline: ' // trim(line), advance='yes' )
100 iColon = index ( line, ':' )
    if ( iColon > 0 ) then
      status = 0
      last = line(iColon+1:)
      if ( DEEBUG ) call output( 'value: ', advance='no' )
      if ( present(string) ) string = adjustl(last)
      if ( present(int) ) call readIntsFromChars( last, int )
      if ( present(log) ) log = ( index( lowerCase(last), 't' ) > 0 )
      if ( DEEBUG ) then
      if ( present(string) ) call output( trim(string), advance='yes' )
      if ( present(int) ) call output( int, advance='yes' )
      if ( present(log) ) call output( log, advance='yes' )
      endif
    else
      status = 1
    endif
    if ( DEEBUG ) then
    call output(trim(line), advance='no')
    call blanks(4)
    call timestamp(status, advance='yes')
    endif
    return
    ! enddo
200 status = 99
    if ( .not. DEEBUG ) return
    call output(trim(line), advance='no')
    call blanks(4)
    call timestamp(status, advance='yes')
    return
400 status = 999
    if ( .not. DEEBUG ) return
    call output(trim(line), advance='no')
    call blanks(4)
    call timestamp(status, advance='yes')
  end subroutine read_oneline
  
  ! ---------------------------------------------  reDumpHosts  -----
  subroutine reDumpHosts(hosts, tempFile, tempUnit)
    type(Machine_T), dimension(:), pointer :: hosts
    character(len=*), intent(in)           :: tempFile
    integer, intent(in)                    :: tempUnit
    ! Internal variables
    character(len=FILENAMELEN)             :: oldPrName
    integer                                :: oldPrUnit
    integer :: status
    logical :: switched
    ! Executable
    switched = .false.
    status = 0
    ! call timeStamp( 'dumping Hosts DB', advance='yes' )
    if ( tempfile /= '<STDIN>' ) then
      if ( options%dump_file /= '<STDIN>' ) then
        switched = .true.
        oldPrName = OutputOptions%name
        close(outputOptions%prunit)
        OutputOptions%name = tempFile
        ! print *, 'switching to ', trim(tempFile)
      endif
      oldPrUnit = OutputOptions%prunit
      OutputOptions%prunit = tempUnit
      ! print *, 'Opening ', prunit, ' as ', trim(tempfile)
      open( OutputOptions%prunit, file=trim(tempfile), &
        & status='replace', form='formatted', iostat=status )
    endif
    call timestamp ( 'Received an external message to dump hostDB', &
      & advance='yes' )
    if ( status == 0 ) then
      call dump(hosts)
    else
      call myMLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Ignoring message; unable to open ' // &
        & trim(tempfile) )
    endif
    if ( tempfile /= '<STDIN>' ) then
      close(OutputOptions%prunit)
      if ( switched ) then
        OutputOptions%name = oldPrName
        open( oldPrUnit, file=trim(options%dump_file), &
          & position='append', &
          & status='old', form='formatted', iostat=status )
        ! print *, 'switching back to ', trim(options%dump_file)
      endif
      OutputOptions%prunit = oldPrUnit
    endif
    if ( options%dump_file /= '<STDIN>' .and. DEEBUG ) &
    & call timeStamp( 'Reverting to ' // trim(options%dump_file), advance='yes' )
  end subroutine reDumpHosts

  ! ---------------------------------------------  reDumpMasters  -----
  subroutine reDumpMasters(masters, tempFile, tempUnit)
    type(Master_T), dimension(:), pointer :: masters
    character(len=*), intent(in)           :: tempFile
    integer, intent(in)                    :: tempUnit
    ! Internal variables
    character(len=FILENAMELEN)             :: oldPrName
    integer                                :: oldPrUnit
    integer :: status
    logical :: switched
    ! Executable
    switched = .false.
    ! call timeStamp( 'Dumping masters DB', advance='yes' )
    status = 0
    if ( tempfile /= '<STDIN>' ) then
      if ( options%dump_file /= '<STDIN>' ) then
        switched = .true.
        oldPrName = OutputOptions%name
        close(outputOptions%prunit)
        OutputOptions%name = tempFile
        ! print *, 'switching to ', trim(tempFile)
      endif
      oldPrUnit = OutputOptions%prunit
      OutputOptions%prunit = tempUnit
      ! print *, 'Opening ', prunit, ' as ', trim(tempfile)
      open( OutputOptions%prunit, file=trim(tempfile), &
        & status='replace', form='formatted', iostat=status )
    endif
    call timestamp ( 'Received an external message to dump masterDB', &
      & advance='yes' )
    if ( status == 0 ) then
      call dump_master_database(masters)
    else
      call myMLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Ignoring message; unable to open ' // &
        & trim(tempfile) )
    endif
    if ( tempfile /= '<STDIN>' ) then
      close(OutputOptions%prunit)
      if ( switched ) then
        OutputOptions%name = oldPrName
        open( oldPrUnit, file=trim(options%dump_file), &
          & position='append', &
          & status='old', form='formatted', iostat=status )
        ! print *, 'switching back to ', trim(options%dump_file)
      endif
      OutputOptions%prunit = oldPrUnit
    endif
    if ( options%dump_file /= '<STDIN>' .and. DEEBUG ) &
    & call timeStamp( 'Reverting to ' // trim(options%dump_file), advance='yes' )
  end subroutine reDumpMasters

  ! ---------------------------------------------  releaseHostFromMaster  -----
  subroutine releaseHostFromMaster(host, master, hostsID)
    type(Machine_T) :: host
    type(master_t) :: master
    integer, intent(in) :: hostsID
    ! Internal variables
    integer :: hcid
    integer, dimension(:), pointer :: tempHosts
    integer :: i, status
    ! Executable
    if ( options%debug ) &
      & call output( 'Releasing ' // trim(host%name) // '   ', advance='yes')
    if( size(master%hosts) < 1 .and. .not. master%owes_thanks ) then
      call myMLSMessage( MLSMSG_Warning, ModuleName, &
        & 'master lacks any hosts' )
      return
    elseif ( hostsID < 1 ) then
      call myMLSMessage( options%errorLevel, ModuleName, &
        & 'Programming error--hostsID=0 in releaseHostFromMaster' )
    elseif ( hostsID > size(hosts) ) then
      call myMLSMessage( options%errorLevel, ModuleName, &
        & 'Programming error--hostsID>size(hosts) in releaseHostFromMaster' )
    elseif ( .not. associated(master%hosts) ) then
      ! call myMLSMessage( MLSMSG_Warning, ModuleName, &
      !  & 'master lacks any hosts' )
      call timestamp('master lacks any hosts', advance='yes' )
      return
    elseif ( host%master_tid /= master%tid ) then
      call output('host%master_tid ', advance='no')
      call output(host%master_tid, advance='no')
      call output('   master%tid ', advance='no')
      call output(master%tid, advance='no')
      call timestamp('Host not assigned to that master', advance='yes' )
      return
    elseif ( .not. any(master%hosts == hostsID) ) then
      call output('hostsID ', advance='no')
      call output(hostsID, advance='no')
      call output('   name ', advance='no')
      call output(trim(host%name) // '   ', advance='no')
      call newline
      call dump_his_hosts(master%hosts, sortUs=.true.)
      ! call myMLSMessage( MLSMSG_Warning, ModuleName, &
      !  & 'master lacks that host' )
      call timestamp( '; master lacks that host, or did he forget?', advance='yes' )
      return
    elseif ( master%owes_thanks .and. .false. ) then
      call myMLSMessage( MLSMSG_Warning, ModuleName, &
        & 'master still owes thanks' )
    endif
    if ( host%tid > 0 ) then
      ! In case slave task is still running, try to kill it
      call pvmfpstat(host%tid, status)
      if ( status == PVMOK ) then
        call pvmfkill(host%tid, status)
        if ( status /= 0 ) then
          call proclaim('Failed to kill running task', &
            & trim(host%name), signal=host%tid)
        elseif ( options%debug ) then
          call timeStamp( 'Forcibly killed task ' // trim(host%name), advance='yes' )
        endif
      endif
    endif
    host%tid = UNASSIGNED
    ! host%tid = 0
    host%master_tid = 0
    host%free = .true.
    host%chunk = -1
    host%master_name = ' '
    host%master_date = ' '
    hcid = 0
    if( .not. associated(master%hosts) ) then
      call myMLSMessage( MLSMSG_Warning, ModuleName, &
        & 'this master has no hosts--how could it have released one?' )
    elseif( size(master%hosts) < 2 ) then
      call deallocate_test(master%hosts, 'master%hosts', moduleName)
    elseif( .not. any(master%hosts == hostsID) ) then
      call myMLSMessage( MLSMSG_Warning, ModuleName, &
        & 'this master does not have this host--how could it have released it?' )
    else
      ! We'll need a smaller hosts component in master, so ..
      nullify(tempHosts)
      hcid = size(master%hosts) - 1
      call allocate_test(tempHosts, hcid, 'temphosts', moduleName)
      hcid = 0
      do i=1, size(master%hosts)
        if ( master%hosts(i) == hostsID ) cycle
        hcid = hcid + 1
        tempHosts(hcid) = master%hosts(i)
      enddo
      call deallocate_test(master%hosts, 'master%hosts', moduleName)
      master%hosts => tempHosts
    endif
    master%numHosts = master%numHosts - 1
  end subroutine releaseHostFromMaster
  
  ! ---------------------------------------------  rmMasterFromDatabase  -----
  subroutine  rmMasterFromDatabase( DATABASE, ITEM )
    ! This function removes a master data type from a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the reduced size
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy arguments
    type (master_t), dimension(:), pointer :: DATABASE
    type (master_t), intent(in)            :: ITEM

    ! Local variables
    integer :: S
    type (Master_T), dimension(:), pointer :: tempDatabase
    logical, parameter                     :: okToDeallocEmptyDB = .false.
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "rmItemFromDatabase.f9h" 

  end subroutine  rmMasterFromDatabase

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

  ! ---------------------------------------------  SwitchMastersAllegiance  -----
  subroutine SwitchMastersAllegiance( Master )
    type(Master_T) :: Master
    ! Internal variables
    integer :: bufferID
    integer :: i
    integer :: info
    ! Executable
    ! call dump_master( master )
    if ( master%finished ) then
      call output(master%tid)
      call timeStamp( ' This master already finished', advance='yes' )
      ! Sometimes we can't trust this claim, so try to send msg anyway
      ! return
    endif
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMFSend ( master%tid, SIG_Switchallegiance, info )
    call timeStamp( ' Sending notice to switch to master', advance='yes' )
    do i=1, 1000
      call PVMFNRecv ( master%tid, sig_swearallegiance, info )
      if ( info > 0 ) exit
      call usleep ( parallel%delay )
    enddo
    if ( i > 1000 ) then
      ! Here we face a choice:
      ! (1) We could kill the unresponsive master; e.g.
      !      call PVMFSend ( master%tid, giveUpTag, info )
      ! (2) We could just mark it as finished, relieving us of further worry
      ! (except that we haven't freed any hosts used by this master)
      master%finished = .true.
      call output(master%tid)
      call timeStamp( ' Received no oath from this master; finished? ', advance='yes' )
    else
      ! Got a response; was this master requesting a host?
      call PVMF90Unpack ( master%needs_Host, info )
      if ( info /= 0 ) &
        & call myPVMErrorMessage ( info, "unpacking signal" )
      call output(master%tid)
      call output( ' Received oath of allegiance from master; needs host? ' )
      call timeStamp( master%needs_Host, advance='yes' )
    endif
    ! In case master sent thanks after prior l2q died, we'll assume 
    ! it does not owe thanks anymore
    master%owes_thanks = .false.
  end subroutine SwitchMastersAllegiance

end program L2Q

! $Log$
! Revision 1.39  2016/08/12 16:18:03  pwagner
! Made consistent with our split of Dump_0
!
! Revision 1.38  2015/04/29 16:18:17  pwagner
! Fixed bug due to changes in rmItemFromDatabase.f9h
!
! Revision 1.37  2014/09/22 18:02:10  pwagner
! Corrected to conform with new add/rmItemFromDatabase
!
! Revision 1.36  2014/01/09 00:31:26  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 1.35  2013/08/23 02:51:48  vsnyder
! Move PrintItOut to PrintIt_m
!
! Revision 1.34  2013/08/14 17:27:26  pwagner
! Moved FindSomethings into MLSFinds
!
! Revision 1.33  2013/07/18 22:36:38  pwagner
! Consistent with having deleted OutputOptions%opened
!
! Revision 1.32  2013/02/14 19:05:29  pwagner
! Added way for l2q to tell master to dump status
!
! Revision 1.31  2012/07/12 17:55:16  pwagner
! Added mdline option --free to regularly check, free abandoned hosts
!
! Revision 1.30  2010/06/03 23:36:07  pwagner
! Tried again to fix problem of freezing
!
! Revision 1.29  2009/09/02 22:05:36  pwagner
! Tried to fix bug preventing some hosts from being revived
!
! Revision 1.28  2009/08/17 20:25:46  pwagner
! Should fix the forgotten needs_host when a master is killed
!
! Revision 1.27  2009/07/21 20:37:52  pwagner
! Print tid when logging masters thanks for hostl2q.f90
!
! Revision 1.26  2009/06/09 17:06:19  pwagner
! Prevent using masters before it becomes allocated
!
! Revision 1.25  2009/04/13 20:34:20  pwagner
! Ignore comments in hosts list file
!
! Revision 1.24  2009/02/10 17:52:04  pwagner
! Fixed syntax error
!
! Revision 1.23  2008/06/17 00:03:13  pwagner
! More tweaking, less freezing
!
! Revision 1.22  2007/09/06 23:37:36  pwagner
! Delay before forcibly killing slaves is masters responsibility, not ours
!
! Revision 1.21  2007/08/31 00:05:26  pwagner
! Fixed bugs preventing mesgs after event_loop from going to stdout
!
! Revision 1.20  2007/06/29 21:00:21  pwagner
! Fixed bug causing crashes when checking for revivals
!
! Revision 1.19  2007/05/18 23:50:15  pwagner
! Now can successfully revive dead hosts
!
! Revision 1.18  2007/02/09 21:24:22  pwagner
! Fixed an error only NAG caught
!
! Revision 1.17  2007/01/18 23:25:59  pwagner
! Changed to use unbuffered output with -o dump_file
!
! Revision 1.16  2007/01/12 00:41:16  pwagner
! May start a new l2q to rescue things if old one dies
!
! Revision 1.15  2006/11/22 18:34:24  pwagner
! Tracks hosts freed by killed masters better
!
! Revision 1.14  2006/11/03 21:27:26  pwagner
! Another last-minute freezing bug fixed (or so we hope)
!
! Revision 1.13  2006/11/01 20:45:29  pwagner
! Fixed another freezing bug
!
! Revision 1.12  2006/09/29 00:31:21  pwagner
! Many life-prolonging changes; may kill masters if so commanded
!
! Revision 1.11  2006/08/02 22:48:22  pwagner
! prunit now a component of OutputOptions
!
! Revision 1.10  2005/11/15 22:39:58  pwagner
! More robust against failures while openeing dump files
!
! Revision 1.9  2005/09/23 21:01:13  pwagner
! use_wall_clock now a component of time_config
!
! Revision 1.8  2005/06/22 19:27:33  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.7  2005/04/12 20:32:19  pwagner
! May periodically autodump masters, hosts files
!
! Revision 1.6  2005/03/18 00:59:36  pwagner
! Now gets PVMERRORMESSAGE from MLSMessageModule
!
! Revision 1.5  2005/02/03 19:10:11  pwagner
! Receives master_date, master_time data from masters for each host
!
! Revision 1.4  2005/01/20 00:54:25  pwagner
! Implementing changes suggested at design review Jan 14 2005
!
! Revision 1.3  2005/01/14 21:39:04  pwagner
! Changes bring us into line with pw presentation of 20050114
!
! Revision 1.2  2004/12/23 00:19:30  pwagner
! New options; more convenient dumping to a file
!
! Revision 1.1  2004/12/09 00:42:46  pwagner
! First commit
!
