! Copyright (c) 2005, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contracts NAS7-1407/NAS7-03001 is acknowledged.

program L2Q
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use L2PARINFO, only: PARALLEL, INITPARALLEL, ACCUMULATESLAVEARGUMENTS
  use L2ParInfo, only: MACHINE_T, PARALLEL, &
    & PETITIONTAG, GIVEUPTAG, GRANTEDTAG, &
    & SIG_FINISHED, SIG_REGISTER, NOTIFYTAG, &
    & SIG_HOSTDIED, SIG_RELEASEHOST, SIG_REQUESTHOST, SIG_THANKSHOST, &
    & MACHINENAMELEN, GETMACHINENAMES, &
    & DW_INVALID, DUMP, addMachineToDatabase
  use MACHINE ! At least HP for command lines, and maybe GETARG, too
  use MLSCOMMON, only: FILENAMELEN
  use MLSL2Options, only: CURRENT_VERSION_ID
  use MLSMessageModule, only: MLSMessage, MLSMessageConfig, MLSMessageExit, &
    & MLSMSG_Allocate, MLSMSG_DeAllocate, MLSMSG_Debug, MLSMSG_Error, &
    & MLSMSG_Info, MLSMSG_Warning
  use MLSSETS, only: FINDFIRST, FINDALL
  use MLSSTRINGLISTS, only: CATLISTS, STRINGELEMENTNUM
  use MLSSTRINGS, only: LOWERCASE
  use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT, OUTPUT_DATE_AND_TIME, PRUNIT, &
    & TIMESTAMP
  use PVM, only: PVMOK, &
    & ClearPVMArgs, FreePVMArgs, GETMACHINENAMEFROMTID, &
    & PVMDATADEFAULT, PVMFINITSEND, PVMF90PACK, PVMFKILL, PVMFMYTID, &
    & PVMF90UNPACK, PVMERRORMESSAGE, PVMFPSTAT, &
    & PVMFCATCHOUT, PVMFSEND, PVMFNOTIFY, PVMTASKEXIT, &
    & PVMFFREEBUF
  use Sort_M, only: SORT
  use Time_M, only: Time_Now, Use_Wall_Clock
  use TOGGLES, only: CON, GEN, LEVELS, PAR, SYN, TAB, &
    & TOGGLE

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

  implicit NONE

  character(len=1) :: arg_rhs       ! 'n' part of 'arg=n'
  integer :: BUFFERID
  character(len=2048) :: command_line ! All the opts
  integer :: ERROR
  integer :: INFO
  integer :: INUNIT = -1       ! Input unit, * if < 0
  character(len=2048) :: LINE      ! Into which is read the command args
  integer, parameter :: LIST_UNIT = 20  ! Unit # for hosts file if not stdin
  integer, parameter :: MAXNUMMULTIPROCS = 8 ! For some architectures > 1000 
  integer :: RECL = 10000          ! Record length for list
  integer :: RECORD_LENGTH
  integer :: STATUS                ! From OPEN
  logical :: SWITCH                ! "First letter after -- was not n"
  real :: T0, T1, T2, T_CONVERSION ! For timing
  integer :: TAG
  integer :: TID

  character(len=*), parameter :: GROUPNAME = "mlsl2"
  character(len=*), parameter :: LISTNAMEEXTENSION = ".txt"
  integer, parameter          :: AVOIDSELECTEDHOSTSTAG = GIVEUPTAG - 1
  integer, parameter          :: CHECKREVIVEDHOSTSTAG = AVOIDSELECTEDHOSTSTAG - 1
  integer, parameter          :: CHECKSELECTEDHOSTSTAG = CHECKREVIVEDHOSTSTAG - 1
  integer, parameter          :: CLEANMASTERDBTAG = CHECKSELECTEDHOSTSTAG - 1
  integer, parameter          :: DUMPDBTAG = CLEANMASTERDBTAG - 1
  integer, parameter          :: DUMPMASTERSDBTAG = DUMPDBTAG - 1
  integer, parameter          :: DUMPHOSTSDBTAG = DUMPMASTERSDBTAG - 1
  integer, parameter          :: SWITCHDUMPFILETAG = DUMPHOSTSDBTAG - 1
  integer, parameter          :: TURNREVIVALSONTAG = SWITCHDUMPFILETAG - 1
  integer, parameter          :: TURNREVIVALSOFFTAG = TURNREVIVALSONTAG - 1
  integer, parameter          :: DUMPUNIT = LIST_UNIT + 1
  integer, parameter          :: TEMPUNIT = DUMPUNIT + 1
  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
     "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it
! LF95.Linux/test [options] [input files]

  ! Our data type for the master tasks we'll be communicating with via pvm
  type master_T
    character(len=MACHINENAMELEN) :: name = ' '
    integer                       :: tid = 0
    integer                       :: numChunks = 0
    integer                       :: numHosts = 0
    integer                       :: numFreed = 0
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
    logical            :: reviveHosts = .false.   ! Regularly check for revivals
    character(len=16)  :: command = 'run'         ! {run, kill, dumphdb, dumpmdb
    logical            :: debug = .false.         !   dump, checkh, clean}
    character(len=FILENAMELEN) &
      &                :: LIST_file               ! name of hosts file
    character(len=FILENAMELEN) &
      &                :: dump_file = '<STDIN>'   ! name of dump file
    logical :: Timing = .false.                   ! -T option is set
    logical :: date_and_times = .false.
    character(len=1) :: timingUnits = 's'         ! l_s, l_m, l_h
    character(len=2048) :: selectedHosts =''      ! E.g., 'c0-1,c0-23,c0-55'
  end type options_T
  
  type ( options_T ) :: options
  type(Machine_T), dimension(:), pointer :: hosts => null()
  type(master_T), dimension(:), pointer :: masters => null()
  !
  MLSMessageConfig%useToolkit = .false.
  MLSMessageConfig%logFileUnit = -1
  MLSMessageConfig%CrashOnAnyError = .true.
  use_wall_clock = .true.
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
      case ( 'kill' )
        tag = GiveUpTag
      case ( 'revive=on' )
        tag = turnRevivalsOnTag
      case ( 'revive=off' )
        tag = turnRevivalsOffTag
      case ( 'switch' )
        tag = switchDumpFileTag
      case default
        call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'l2q would not recognize unknown command '//trim(options%command) )
      end select
      if ( any(tag == (/avoidSelectedHostsTag, checkSelectedHostsTag/) )) then
        call PVMF90Pack ( trim(options%selectedHosts), info )
      else
        call PVMF90Pack ( trim(options%dump_file), info )
      endif
      call PVMFSend ( tid, Tag, info )
      if ( options%timing ) call sayTime
      call output('Already-running l2q tid: ', advance='no')
      call timestamp(tid, advance='yes')
      call timestamp(' commanded: '//trim(options%command), advance='yes')
      if ( any(tag == (/avoidSelectedHostsTag, checkSelectedHostsTag/) )) then
        call timestamp(trim(options%selectedHosts), advance='yes')
      elseif ( any(tag == (/dumpHostsDBTag,dumpMastersDBTag, dumpDBTag, &
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
  ! Join (actually start) our group
  call PVMFJoinGroup ( GROUPNAME, status )
  if ( status < 0 ) &
    & call PVMErrorMessage ( status, 'Joining '//GROUPNAME//' group' )
! Clear the command line arguments we're going to accumulate to pass
! to slave tasks
  call ClearPVMArgs

  if ( options%dump_file /= '<STDIN>' ) then
    prunit = DUMPUNIT
    open(prunit, file=trim(options%dump_file), &
      & status='replace', form='formatted')
  endif
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
  
  call output_date_and_time(msg='starting l2q')
  call time_now ( t1 )

  if( options%verbose ) then
    call dump_settings
  end if

  call read_list
  call cure_host_database(hosts, silent=.true.)
  if ( options%checkList ) then
    call dump(hosts)
  !---------------- Task (4) ------------------
  elseif (error == 0) then
    call dump(hosts, details=0)
    if( options%verbose ) call timestamp('Entering event loop', advance='yes')
    call event_loop
  else
    call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'error reading list of hosts' )
  end if

  call time_now ( t0 )
  t1 = t0
  if ( options%timing ) call sayTime ( 'Closing and deallocating' )
  call output_date_and_time(msg='ending l2q')
  if ( options%dump_file /= '<STDIN>' ) close(prunit)

contains
!   integer function  addHostToDatabase( DATABASE, ITEM )
!     ! This function adds a host data type to a database of said types,
!     ! creating a new database if it doesn't exist.  The result value is
!     ! the size -- where the new item is put.
! 
!     ! Dummy arguments
!     type (Machine_T), dimension(:), pointer :: DATABASE
!     type (Machine_T), intent(in) :: ITEM
! 
!     ! Local variables
!     type (Machine_T), dimension(:), pointer :: tempDatabase
!     !This include causes real trouble if you are compiling in a different 
!     !directory.
!     include "addItemToDatabase.f9h" 
! 
!     AddHostToDatabase = newSize
!   end function  addHostToDatabase

  subroutine cure_host_database(hosts, silent, selectedHosts)
    ! cure any recovered hosts in database--declare them fit for action
    ! If optional selectedHosts is present, cure only them
    type(Machine_T), intent(inout), dimension(:) :: hosts
    logical, optional, intent(in) :: silent
    character(len=*), optional, intent(in) :: selectedHosts
    ! Internal variables
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    logical :: mySilent
    integer :: i
    integer :: status
    ! Executable
    mySilent = .not. options%verbose
    if ( present(silent) ) mySilent = silent .or. mySilent
    ! No need to cure if none have died
    if ( all(hosts%OK) ) return
    call GetMachineNames ( machineNames )
    do i=1, size(hosts)
      if ( hosts(i)%OK ) cycle
      ! Warning--following use of variable 'status' is non-standard
      ! status is good if non-zero, bad if 0
      status = 1
      if ( present(selectedHosts) ) &
        & status = StringElementNum(selectedHosts, trim(hosts(i)%name), .true.)
      if ( status < 1 ) cycle
      hosts(i)%OK = any( hosts(i)%name == machineNames )
      if ( hosts(i)%OK .and. .not. mySilent ) &
        & call proclaim('Host Cured', hosts(i)%Name)
    enddo
    call deAllocate_test(machineNames, 'machineNames', moduleName )
  end subroutine cure_host_database

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

  integer function  addMasterToDatabase( DATABASE, ITEM )
    ! This function adds a master data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where the new item is put.

    ! Dummy arguments
    type (master_t), dimension(:), pointer :: DATABASE
    type (master_t), intent(in) :: ITEM

    ! Local variables
    type (Master_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    addMasterToDatabase = newSize
  end function  addMasterToDatabase

  subroutine  rmMasterFromDatabase( DATABASE, ITEM )
    ! This function removess a master data type from a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the reduced size

    ! Dummy arguments
    type (master_t), dimension(:), pointer :: DATABASE
    type (master_t), intent(in)            :: ITEM

    ! Local variables
    type (Master_T), dimension(:), pointer :: tempDatabase
    logical, parameter                     :: okToDeallocEmptyDB = .false.
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "rmItemFromDatabase.f9h" 

  end subroutine  rmMasterFromDatabase

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

  subroutine dump_Master(Master)
    ! dump data type
    integer :: i
    type(Master_T), intent(in) :: Master
    call output('machine name: ', advance = 'no')
    call output(trim(Master%name), advance = 'yes')
    call output('tid: ', advance = 'no')
    call output(Master%tid, advance = 'yes')
    call output('needs host?: ', advance = 'no')
    call output(Master%needs_host, advance = 'yes')
    call output('owes thanks?: ', advance = 'no')
    call output(Master%owes_thanks, advance = 'yes')
    call output('finished?: ', advance = 'no')
    call output(Master%finished, advance = 'yes')
  end subroutine dump_Master
   
  subroutine get_options
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
      command_line = trim(command_line) // ' ' // trim(line)
      if ( line(1:2) == '--' ) then       ! "word" options
        n = 0
        switch = .true.
        if ( line(3:3) == 'n' .or. line(3:3) == 'N' ) then
          switch = .false.
          n = 1
        end if
        if ( line(3+n:7+n) == 'clean ' ) then
          options%cleanMasterDB = switch
        elseif ( line(3+n:7+n) == 'check ' ) then
          options%checklist = switch
        elseif ( line(3+n:6+n) == 'defer' ) then
          options%deferToElders = switch
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
        else if ( line(3+n:7+n) == 'reviv' ) then
          options%reviveHosts = switch
        else if ( line(3+n:7+n) == 'wall ' ) then
          use_wall_clock = switch
        else if ( line(3:) == ' ' ) then  ! "--" means "no more options"
          i = i + 1
          call getarg ( i, line )
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
        case ( 'k' ); options%command = 'kill' ! options%killer = .true.
        case ( 'c' )
          i = i + 1
          call getarg ( i, line )
          options%command = lowercase(line)
          ! Special commands take yet another arg
          if ( index('check,avoid', trim(options%command)) > 0 ) then
            i = i + 1
            call getarg ( i, line )
            options%selectedHosts = lowercase(line)
          endif
        case ( 'o' )
          i = i + 1
          call getarg ( i, line )
          options%dump_file = line
        case ( 's' )
          options%command = 'switch'
          i = i + 1
          call getarg ( i, line )
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
   
  subroutine event_loop
    ! The main event loop
    ! l2q spends eternity circling endlessly, responding
    ! to the occasional promptings from master tasks via pvm
    ! or else initiating a contact with one
    ! or responding to commands sent by a fraternal l2q launched with
    ! "-c command" option
    integer :: bufferIDRcv
    integer :: BUFFERIDSND              ! From PVM
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    logical :: CHECKREVIVEDHOSTS
    logical :: CLEANMASTERDB
    logical :: DUMPDB
    integer :: grandMastersID           ! index into database of an older master
    integer :: host
    integer :: hostsID                  ! index into database of a host
    integer :: INFO                     ! From PVM
    character(len=MachineNameLen)  :: MACHINENAME
    integer :: mastersID                ! index into database of a master
    integer :: masterTid                ! TID of a master
    logical :: mayAssignAHost
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: nextFree
    integer :: NTIDS
    integer :: numHosts
    integer :: numMasters
    integer :: oldPrUnit
    integer :: SIGNAL                   ! From a master
    logical :: SIGNIFICANTEVENT
    logical :: SKIPDELAY            ! Don't wait before doing the next go round
    character(len=FileNameLen)  :: tempfile
    integer :: TID
    integer, dimension(MAXNUMMULTIPROCS) :: TIDS
    type(master_t) :: aMaster
    ! Executable
    numMasters = 0
    do
      checkrevivedhosts = .false.
      dumpdb = .false.
      mayAssignAHost = .true.
      significantEvent = .false.
      skipDelay = .false.               ! Assume we're going to delay
      ! Listen out for pvm events
      call PVMFNRecv( -1, PETITIONTAG, bufferIDRcv )
      if ( bufferIDRcv < 0 ) then
        call PVMErrorMessage ( info, "checking for Petition message" )
      else if ( bufferIDRcv > 0 ) then
        ! So we got a message.  There may be another one following on so don't 
        ! delay before going round this loop again.
        skipDelay = .true.
        ! Who sent this?
        call PVMFBufInfo ( bufferIDRcv, bytes, msgTag, masterTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking signal" )
        endif
        ! if ( options%debug ) call proclaim('Master sent message', signal=signal)
        select case (signal) 

        case ( sig_register ) ! ----------------- Master registering ------
          call PVMF90Unpack ( aMaster%NumChunks, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking number of chunks" )
          endif
          aMaster%tid = masterTid
          call GetMachineNameFromTid ( masterTid, aMaster%Name )
          aMaster%NumHosts = 0
          aMaster%needs_Host = .false.
          aMaster%owes_thanks = .false.
          aMaster%finished = .false.
          numMasters = addMasterToDatabase(masters, aMaster)
          significantEvent = .true.
          ! masterName = catlists('m-', numMasters)
          ! Now ask to be notified when this master exits
          call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ masterTid /), info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'setting up notify' )
          if ( options%verbose ) then
            ! call proclaim('New Master', aMaster%Name, signal=MasterTid, advance='no')
            call proclaim('New Master', masterNameFun(masterTid), advance='no')
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
          if ( options%dumpEachNewMaster ) call dump_master(aMaster)
        case ( sig_requestHost ) ! ----------------- Request for host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          grandMastersID = FindFirst(masters%needs_host .and. &
            & .not. masters%owes_thanks)
          if ( options%deferToElders ) &
            grandMastersID = FindFirst(.not. masters%finished .and. &
              & masters%numChunks > (masters%numHosts+masters%numFreed) )
          if ( grandMastersID < 1 ) grandMastersID = size(masters) + 100
          if ( options%debug ) &
            & call proclaim('Master requested host', masterNameFun(masterTid))
          ! Any free hosts?
          nextFree = FindFirst( hosts%free .and. hosts%OK )
          if ( masters(mastersID)%owes_thanks ) then
            ! Don't assign hosts to ungrateful masters
            if ( options%verbose ) &
              & call proclaim('Ungrateful master requested another host', &
              & masterNameFun(masterTid), advance='yes')
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
          call IDHostFromMaster(sig_releaseHost, tid, mastersID, hostsID)
          if ( hostsID > 0 ) then
            call releaseHostFromMaster( hosts(hostsID), masters(mastersID), &
              & hostsID )
            if ( options%debug ) then
              call proclaim('Master ' // trim(masterNameFun(masterTid)) // &
                & ' freed host', hosts(hostsID)%name, advance='no')
              call output(' (Still has ', advance='no')
              call output(masters(mastersID)%numHosts, advance='no')
              call timestamp(' remaining)', advance='yes')
              call output ('Number of machines free: ', advance='no')
              call timestamp (count(hosts%free), advance='yes')
            endif
            mayAssignAHost = .false. ! Don't assign to a later master
          endif
        case ( sig_ThanksHost ) ! -------------- Happy with this host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          masters(mastersID)%owes_thanks = .false.
          ! if ( options%debug ) &
          !   & call proclaim('Master thanks for host', signal=masterTid)
          call PVMF90Unpack ( tid, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking tid" )
          endif
          call PVMF90Unpack ( machineName, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking machineName" )
          endif
          hostsID = FindFirst(  hosts%name==machineName .and. hosts%tid < 1 )
          if ( hostsID < 1 ) then
            call output(' mastersID ', advance='no')
            call output(mastersID , advance='yes')
            call output(' machineName ', advance='no')
            call output(trim(machineName) , advance='yes')
            call MLSMessage( MLSMSG_Error, ModuleName, &
               & 'Master thanked us for unknown host' )
          endif
          hosts(hostsID)%tid = tid
          call PVMF90Unpack ( hosts(hostsID)%chunk, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking machine%chunk" )
          endif
          call PVMF90Unpack ( hosts(hostsID)%master_date, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking machine%master_date" )
          endif
          ! call PVMF90Unpack ( hosts(hostsID)%master_name, info )
          ! if ( info /= 0 ) then
          !   call PVMErrorMessage ( info, "unpacking machine%master_name" )
          ! endif
          hosts(hostsID)%master_name = masterNameFun(masterTID)
          if ( options%verbose ) then
            call proclaim('Master ' // trim(masterNameFun(masterTID)) // &
            & ' thanks for host', hosts(hostsID)%Name, advance='no')
            call output(' (Now has ', advance='no')
            call output(masters(mastersID)%numHosts, advance='no')
            call timestamp(' )', advance='yes')
          endif
        case ( sig_HostDied ) ! -------------- Done away with this host ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          if ( options%debug ) &
            & call proclaim('Master ' // trim(masterNameFun(masterTID)) // &
            & ' reports host died', advance='no')
          call IDHostFromMaster(sig_HostDied, tid, mastersID, hostsID)
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
            ! The next snippet has been superseded by mark_hosts_database
            ! Now see if any other hosts are associated with this machine name
            ! (we'll assume the processors on amultiprocessor node are either 
            !  all good or all bad)
            ! if ( tid < 1 ) then
            !   call FindAll(hosts%name, machineName, Tids, nTids)
            !   if ( nTids > 1 ) then
            !     do host=2, nTids
            !       hosts(host)%OK = .false.
            !     enddo
            !   endif
            ! endif
          endif
        case ( sig_Finished ) ! ----------------- Done with this master ------
          ! Find master's index into masters array
          mastersID = FindFirst(masters%tid, masterTid)
          numHosts = masters(mastersID)%numHosts
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
            call proclaim('Master ' // trim(masterNameFun(masterTid)) // &
              & ' Finished', advance='no')
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
        case default
          if ( options%debug ) &
            & call proclaim( 'Unrecognized signal from Master ' // &
            & masterNameFun(masterTid) )
          call MLSMessage( MLSMSG_Error, ModuleName, &
             & 'Unrecognized signal from master' )
        end select
        call PVMFFreeBuf ( bufferIDRcv, info )
        if ( options%debug .and. info /= 0 ) &
          & call timestamp('Trouble freeing buf', advance='yes')
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'freeing receive buffer' )
      endif
      ! ----------------------------------------------------- Master Task Died?
      ! Listen out for any message telling us one of our masters died
      call PVMFNRecv ( -1, NotifyTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        ! Get the TID for the dead task
        call PVMF90Unpack ( masterTid, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking dead masterTid' )
        ! Now this may well be a legitimate exit, 
        ! in which case, master finished normally.  
        mastersID = FindFirst(masters%tid, masterTid)
        if ( mastersID < 1 ) call MLSMessage( MLSMSG_Error, ModuleName, &
           & 'Exit signal from unrecognized master' )
        if ( .not. masters(mastersID)%finished ) then
          ! Otherwise we need to tidy up.
          numHosts = masters(mastersID)%numHosts
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
            call proclaim('Master Died', masterNameFun(masterTid))
            call output(' (', advance='no')
            call output(count(.not. masters%finished), advance='no')
            call output(' active ', advance='no')
            call output(numMasters, advance='no')
            call timestamp(' total)', advance='yes')
          endif
        endif
      endif
      ! ----------------------------------------------------- Administrative messages?
      ! Listen out for any message telling *us* to quit now
      call PVMFNRecv ( -1, GiveUpTag, bufferIDRcv )
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

      checkrevivedhosts = significantEvent .and. options%reviveHosts
      ! Listen out for any message telling us to dump databases
      call PVMFNRecv ( -1, DumpDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking dumpfile' )
        if ( tempfile /= '<STDIN>' ) then
          oldPrUnit = prunit
          prUnit = tempUnit
          open(prunit, file=trim(tempfile), &
            & status='replace', form='formatted')
        endif
        call timestamp ( 'Received an external message to dump databases', &
          & advance='yes' )
        call dump_master_database(masters)
        call dump(hosts)
        if ( tempfile /= '<STDIN>' ) then
          close(prunit)
          prunit = oldPrUnit
        endif
      end if
      call PVMFNRecv ( -1, DumpHostsDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking dumpfile' )
        if ( tempfile /= '<STDIN>' ) then
          oldPrUnit = prunit
          prUnit = tempUnit
          open(prunit, file=trim(tempfile), &
            & status='replace', form='formatted')
        endif
        call timestamp ( 'Received an external message to dump hostDB', &
          & advance='yes' )
        call dump(hosts)
        if ( tempfile /= '<STDIN>' ) then
          close(prunit)
          prunit = oldPrUnit
        endif
      end if
      call PVMFNRecv ( -1, DumpMastersDBTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( tempfile, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking dumpfile' )
        if ( tempfile /= '<STDIN>' ) then
          oldPrUnit = prunit
          prUnit = tempUnit
          open(prunit, file=trim(tempfile), &
            & status='replace', form='formatted')
        endif
        call timestamp ( 'Received an external message to dump masterDB', &
          & advance='yes' )
        call dump_master_database(masters)
        if ( tempfile /= '<STDIN>' ) then
          close(prunit)
          prunit = oldPrUnit
        endif
      end if
      ! Listen out for any message telling us to avoid selected hosts
      call PVMFNRecv ( -1, avoidSelectedHostsTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to avoid ' // &
          & trim(options%selectedHosts), advance='yes' )
        call mark_host_database(hosts, options%selectedHosts)
      end if
      ! Listen out for any message telling us to check for revived hosts
      call PVMFNRecv ( -1, checkRevivedHostsTag, bufferIDRcv )
      ! Do it if so commanded, or if part of regular drill
      if ( bufferIDRcv > 0 ) call timestamp ( &
        & 'Received an external message to check for revivals', advance='yes' )
      if ( bufferIDRcv > 0 .or. checkrevivedhosts ) &
        & call cure_host_database(hosts)
      ! Listen out for any message telling us to revive selected hosts
      call PVMFNRecv ( -1, checkSelectedHostsTag, bufferIDRcv )
      ! Do it if so commanded
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        call PVMF90Unpack ( options%selectedHosts, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking selectedHosts' )
        call timestamp ( 'Received an external message to check ' // &
          & trim(options%selectedHosts), advance='yes' )
        call cure_host_database(hosts, .true., options%selectedHosts)
      end if
      ! Listen out for any message telling us to flush current dumpfile
      ! and switch future output to new one
      call PVMFNRecv ( -1, switchDumpFileTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking GROUPNAME' )
        if ( prunit == DUMPUNIT ) then
          if ( options%timing ) call sayTime
          call output('Switching dump file from '//trim(options%dump_file))
          call output(' to ')
          call PVMF90Unpack ( options%dump_file, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking newdumpfile' )
          call timestamp(trim(options%dump_file), advance='yes')
          close(DUMPUNIT)
          open(prunit, file=trim(options%dump_file), &
            & status='replace', form='formatted')
          call timestamp('(new dump file opened)', advance='yes')
        else
          call timestamp('(switchDumpFileTag ignored--dumping to <STDIN>)', &
            & advance='yes')
        end if
      end if
      ! Unless there's a good reason not to do so,
      ! check if any hosts are currently free and any masters need hosts
      if ( mayAssignAHost ) then
        ! Find oldest master needing a host and grateful for the ones he has
        mastersID = FindFirst( masters%needs_host .and. &
          & .not. masters%owes_thanks )
        grandMastersID = mastersID
        if ( options%deferToElders ) &
          grandMastersID = FindFirst(.not. masters%finished .and. &
            & masters%numChunks > (masters%numHosts+masters%numFreed) )
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
      ! Now, rather than chew up cpu time on this machine, we'll wait a
      ! bit here.
      if ( .not. skipDelay .and. parallel%delay > 0 ) then
        call usleep ( parallel%delay )
      endif
    enddo
  end subroutine event_loop
  
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
      call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'Tried to assign an unqualified host to this master' )
    elseif( master%owes_thanks ) then
      call MLSMessage( MLSMSG_Error, ModuleName, &
      & 'Tried to assign host to an ungrateful master' )
    endif      
    host%master_tid = master%tid
    host%free = .false.
    host%tid = 0
    if ( options%debug) then
       call output ('Number of machines free: ', advance='no')
       call timestamp (count(hosts%free), advance='yes')
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
    master%hosts(hcid) = hostsid
    ! Now use PVM to tell master he's got a new slave host
    call PVMFInitSend ( PvmDataDefault, bufferID )
    ! call PVMF90Pack ( sig_requestHost, info )
    ! if ( info /= 0 ) &
    !   & call PVMErrorMessage ( info, 'packing registration' )
    call PVMF90Pack ( trim(host%name), info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing host name' )
    if ( options%debug ) &
      & call proclaim('Sending new host ' // trim(host%name), &
      & masterNameFun(master%Tid))
    call PVMFSend ( master%tid, grantedTag, info )
    master%owes_thanks = .true.
    master%needs_host = .false.
    master%numHosts = master%numHosts + 1
  end subroutine assignHostToMaster
  
  subroutine IDHostFromMaster(signal, tid, mastersID, hostsID)
    ! Figure out which master is talking about
    ! Arguments
    integer, intent(in)  :: signal
    integer, intent(out) :: tid
    integer, intent(in)  :: mastersID
    integer, intent(out) :: hostsID
    ! Internal variables
    integer :: info
    character(len=MACHINENAMELEN) :: machineName
    ! Executable
    hostsID = 0
    call PVMF90Unpack ( tid, info )
    if ( info /= 0 ) then
      call PVMErrorMessage ( info, "unpacking tid" )
    endif
    call PVMF90Unpack ( machineName, info )
    if ( info /= 0 ) then
      call PVMErrorMessage ( info, "unpacking machineName" )
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
          call output('Warning--unknown tid freed', advance='no')
        else
          call output('Warning--unknown tid died', advance='no')
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
        if ( hostsID < 1 ) call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'Even the machine name ' // trim(machineName) // &
          & ' yields an ID outside list of hosts' )
      else
        call MLSMessage( MLSMSG_Warning, ModuleName, &
          & 'machine name yields an ID outside list of hosts' )
      endif
    endif
  end subroutine IDHostFromMaster

  subroutine releaseHostFromMaster(host, master, hostsID)
    type(Machine_T) :: host
    type(master_t) :: master
    integer, intent(in) :: hostsID
    ! Internal variables
    integer :: hcid
    integer, dimension(:), pointer :: tempHosts
    integer :: i, status
    ! Executable
    if( size(master%hosts) < 1 ) then
      call MLSMessage( MLSMSG_Warning, ModuleName, &
        & 'master lacks any hosts' )
      return
    elseif ( hostsID < 1 ) then
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'Programming error--hostsID=0 in releaseHostFromMaster' )
    elseif ( hostsID > size(hosts) ) then
      call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'Programming error--hostsID>size(hosts) in releaseHostFromMaster' )
    elseif ( .not. associated(master%hosts) ) then
      ! call MLSMessage( MLSMSG_Warning, ModuleName, &
      !  & 'master lacks any hosts' )
      call timestamp('master lacks any hosts', advance='yes' )
      return
    elseif ( host%master_tid /= master%tid ) then
      call timestamp('Host not assigned to that master', advance='yes' )
      return
    elseif ( .not. any(master%hosts == hostsID) ) then
      call output('hostsID ', advance='no')
      call output(hostsID, advance='no')
      call output('   name ', advance='no')
      call output(trim(host%name) // '   ', advance='no')
      call newline
      call dump_his_hosts(master%hosts, sortUs=.true.)
      ! call MLSMessage( MLSMSG_Warning, ModuleName, &
      !  & 'master lacks that host' )
      call timestamp( '; master lacks that host', advance='yes' )
      return
    endif
    if ( host%tid > 0 ) then
      ! In case slave task is still running, try to kill it
      call pvmfpstat(host%tid, status)
      if ( status == PVMOK ) then
        call pvmfkill(host%tid, status)
        if ( status /= 0 ) then
          call proclaim('Failed to kill running task', &
            & trim(host%name), signal=host%tid)
        endif
      endif
    endif
    host%tid = 0
    host%master_tid = 0
    host%free = .true.
    host%chunk = -1
    host%master_name = ' '
    host%master_date = ' '
    hcid = 0
    if( size(master%hosts) == 1 ) then
      call deallocate_test(master%hosts, 'master%hosts', moduleName)
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
    master%numFreed = master%numFreed + 1
  end subroutine releaseHostFromMaster
  
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
100     if ( options%verbose ) &
          & call output(trim(newhost%name) // ' added', advance='yes')
      newhost%master_tid = -1
      newhost%tid = -1
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

  ! The following offer some information on the options
  ! available either on the command-line or via the PCF
  ! Note the unashamed use of 'print' statements which
  ! are officially discouraged in favor of calls to MLSMessage.
  ! Unfortunately, we have not yet decided which method to use
  ! until *after* processing all the options.
  
  subroutine Option_usage
    call getarg ( 0+hp, line )
    print *, 'Usage: ', trim(line), ' [options] [--] [LIST-name]'
    print *, ' Options:'
    print *, ' --check:     check LIST of hosts'
    print *, ' --[n]clean:  do [not] regularly clean db of finished masters'
    print *, ' --defer:     new masters defer to oldest (has enough hosts)'
    print *, ' --dumpnewm:  dump each new master'
    print *, ' --help:      show help; stop'
    print *, ' --inquire:   show whether an instance of l2q already running'
    print *, ' --[n]revive: do [not] regularly check for revived hosts'
    print *, ' --version:   print version string; stop'
    print *, ' --wall:      show times according to wall clock (if T[0] set)'
    ! print *, ' -g:          turn tracing on'
    print *, ' -k:          kill currently-running l2q'
    print *, '               ( same as -c kill )'
    print *, ' -o dumpfile: direct most output to dumpfile'
    print *, '               (if modifying an already-running l2q, direct only'
    print *, '               the requested output to dumpfile)'
    print *, ' -s newfile:  flush current dumpfile; direct later output'
    print *, '                to newfile ( same as -o newfile -c switch )'
    print *, ' -v:          verbose'
    print *, ' -h:          show help; stop'
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
    print *, ' revive=on:   begin regularly checking for revived hosts'
    print *, ' revive=off:  cease regularly checking for revived hosts'
    print *, ' kill:        kill currently-running l2q'
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
    call output(' Regularly clean db of finished masters:         ', advance='no')
    call blanks(4, advance='no')
    call output(options%cleanMasterDB, advance='yes')
    call output(' New masters defer to old:                       ', advance='no')
    call blanks(4, advance='no')
    call output(options%deferToElders, advance='yes')
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
    call output(use_wall_clock, advance='yes')
    call output(' ----------------------------------------------------------', &
      & advance='yes')
  end subroutine Dump_settings

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

  ! ---------------------------------------------  proclaim  -----
  subroutine proclaim( Event, Name, Signal, advance )
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

  ! ---------------------------------------------  announce_success  -----
  subroutine announce_success ( Name, unit_number )
    character(LEN=*), intent(in)   :: Name
    integer, intent(in)            :: unit_number

    call output ( 'List of hosts file ' )
    call output ( 'name : ' )
    call blanks(4)
    call output ( trim(Name), advance='no')
    call blanks(10)
    call output ( 'unit number : ' )
    call blanks(4)
    call timestamp ( unit_number, advance='yes')
  end subroutine announce_success

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

end program L2Q

! $Log$
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
