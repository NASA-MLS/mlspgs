! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module L2ParInfo
  ! This module provides definitions needed by L2Parallel and other modules to
  ! manage the parallel aspects of the L2 code.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: Dump
  use HighOutput, only: OutputNamedValue
  use Io_Stuff, only: Get_Lun
  use Machine, only: Usleep
  use MLSFinds, only: Findfirst
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, PVMErrorMessage
  use MLSStats1, only: Allstats
  use MLSStrings, only: Lowercase
  use MLSStringLists, only: Getuniqueints, Getuniquestrings, NumstringElements, &
    & Removenumfromlist, StringElement, Switchdetail
  use MorePVM, only: PVMpackstringindex, PVMunpackstringindex
  use Output_M, only: Beep, Output
  use PVM, only: Infotag, &
    & PVMfmytid, PVMfinitsend, PVMf90pack, PVMfsend, &
    & PVMDatadefault, PVMf90unpack, NextPVMarg, PVMtaskexit, Sig_Abouttodie
  use PVMidl, only: PVMidlpack
  use QuantityPVM, only: PVMsendquantity
  use Toggles, only: Switches
  use Vectorsmodule, only: VectorValue_T

  implicit none
  private

  public :: L2parallelinfo_T, Parallel, Machine_T
  public :: AddmachinenametoDatabase, AddmachinetoDatabase
  public :: Initparallel, Closeparallel
  public :: Slavejoin
  public :: Getnicetidstring, Slavearguments
  public :: Accumulateslavearguments, Logdirectwriterequest
  public :: Sniplastslaveargument, Transmitslavearguments
  public :: Finisheddirectwrite, Machinenamelen, Getmachinenames, Getmachines
  public :: Fwmslavegroup, Directwriterequest_T
  public :: Inflatedirectwriterequestdb, Waitfordirectwritepermission
  public :: Compactdirectwriterequestdb, Dump
  public :: SigToName, TagToName
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Parameters

  ! So requests don't pile up delaying notification of the finishings
  integer, parameter :: DELAYFOREACHSLAVEDWREQUEST   = 0 ! unneeded 200000
  ! To leave time for slaves to flush stdout before master closes pvm
  integer, parameter :: DELAYFOREACHSLAVESTDOUTBUFFER   = 500
  integer, parameter :: FIXDELAYFORSLAVESTDOUTBUFFER   = 500000
  integer, public, parameter :: CHUNKTAG   = InfoTag + 1  ! Master => slave: chunkinfo
  integer, public, parameter :: NOTIFYTAG  = ChunkTag + 1  ! pvm => master: Slave exited
  integer, public, parameter :: PETITIONTAG  = NotifyTag + 1 ! master => l2q
  integer, public, parameter :: GRANTEDTAG  = PetitionTag + 1 ! l2q => master
  integer, public, parameter :: MASTERDUMPTAG  = GrantedTag + 1 ! l2q => master
  integer, public, parameter :: MACHINEFIXEDTAG = 800
  integer, public, parameter :: GIVEUPTAG  = 999

  integer, public, parameter :: SIG_TOJOIN = SIG_AboutToDie + 1
  integer, public, parameter :: SIG_FINISHED = SIG_toJoin + 1
  integer, public, parameter :: SIG_ACKFINISH = SIG_finished + 1
  integer, public, parameter :: SIG_REGISTER = SIG_AckFinish + 1
  integer, public, parameter :: SIG_REQUESTDIRECTWRITE = SIG_Register + 1
  integer, public, parameter :: SIG_DIRECTWRITEGRANTED = SIG_RequestDirectWrite + 1
  integer, public, parameter :: SIG_DIRECTWRITEFINISHED = SIG_DirectWriteGranted + 1
  integer, public, parameter :: SIG_NEWSETUP = SIG_DirectWriteFinished + 1
  integer, public, parameter :: SIG_RUNMAF = SIG_NewSetup + 1
  integer, public, parameter :: SIG_SENDRESULTS = SIG_RunMAF + 1
  integer, public, parameter :: SIG_REQUESTHOST = SIG_SendResults + 1
  integer, public, parameter :: SIG_RELEASEHOST = SIG_RequestHost + 1
  integer, public, parameter :: SIG_HOSTDIED = SIG_ReleaseHost + 1
  integer, public, parameter :: SIG_THANKSHOST = SIG_HostDied + 1
  integer, public, parameter :: SIG_SWEARALLEGIANCE = SIG_ThanksHost + 1
  integer, public, parameter :: SIG_SWITCHALLEGIANCE = SIG_SwearAllegiance + 1

  integer, parameter :: MACHINENAMELEN = 64 ! Max length of name of machine

  ! Name of fwm slave group
  character (len=*), parameter :: FWMSLAVEGROUP = "MLSL2FWMSlaves"

  ! This datatype defines configuration for the parallel code
  type L2ParallelInfo_T
    logical :: fwmParallel = .false.    ! Set if we are in forward model parallel mode
    logical :: master = .false.         ! Set if this is a master task
    logical :: slave = .false.          ! Set if this is a slave task
    integer :: verbosity = 0            ! Set > 0 for extra output 
    integer :: myTid                    ! My task ID in pvm       | than a file
    integer :: masterTid                ! task ID in pvm
    integer :: noFWMSlaves              ! No. slaves in pvm system for fwm cases
    integer :: chunkNo                  ! Chunk No. if this is a slave
    integer :: Delay = 200000           ! For Usleep, no. microsecs
    character(len=132) :: pgeName="mlsl2"   ! command name, if not 'mlsl2'
    character(len=132) :: slaveFilename ! Filename with list of slaves
    character(len=132) :: executable    ! Executable filename
    character(len=132) :: submit=""     ! Submit comand for batch queue system
    ! run range in comma-separated list with possible ranges
    ! e.g. '1,2-6+2,9-11,15' expands to '1,2,4,6,9,10,11,15'
    character(len=4096) :: chunkRange='' ! if blank, runs all chunks
    ! chunks are ','-separated ints; e.g. '2,5,129'
    character(len=4096) :: failedChunks='' ! if blank, no chunks failed
    ! Machs are string list; e.g. 'c0-1,c0-66,c0-66'
    character(len=4096) :: failedMachs='' ! blank if usingSubmit
    ! Msgs are string list; e.g. 'msg 1\msg 2\..\msg n'
    character(len=8192) :: failedMsgs='' ! in same order (hopefully)
    integer :: maxFailuresPerMachine = 1 ! More than this then don't use it | staging
    integer :: maxFailuresPerChunk = 10 ! More than this then give up on getting it
    integer :: numFailedChunks = 0
    integer :: numCompletedChunks = 0
  end type L2ParallelInfo_T

  ! This enumerated type describes the state that directWrites can be in
  integer, public, parameter :: DW_INVALID = 0
  integer, public, parameter :: DW_PENDING = DW_INVALID + 1
  integer, public, parameter :: DW_INPROGRESS = DW_PENDING + 1
  integer, public, parameter :: DW_COMPLETED = DW_INPROGRESS + 1

  ! This datatype logs a directWrite request
  type DirectWriteRequest_T
    integer :: CHUNK=0                  ! Which chunk made the request
    integer :: MACHINE=0                ! What machine is that on
    integer :: NODE=0                   ! What tree node was it
    integer :: FILEINDEX=0              ! Which file was it
    integer :: TICKET=0                 ! What ticket number is it
    integer :: VALUE=0                  ! Workspace for master
    integer :: STATUS=DW_INVALID        ! One of the DW_... above
    real    :: WHENMADE
    real    :: WHENGRANTED
    real    :: WHENFINISHED = -999.99
  end type DirectWriteRequest_T

  ! This datatype describes the machines that will host parallel tasks
  type machine_T
    character(len=MACHINENAMELEN) :: name = ' '
    character(len=MACHINENAMELEN) :: master_name = ' '
    character(len=16            ) :: master_date = ' '
    integer                       :: master_tid = 0
    integer                       :: tid = 0
    integer                       :: chunk = 0
    integer                       :: jobsKilled = 0
    logical                       :: OK = .true.
    logical                       :: free = .true.
  end type machine_T

  ! Shared variables
  type (L2ParallelInfo_T), save :: parallel
  character ( len=2048 ), save :: slaveArguments = ""

  interface dump
    module procedure DumpDirectWriteRequest, DumpAllDirectWriteRequests
    module procedure Dump_Machine, Dump_Machine_DataBase
  end interface

  interface LogDirectWriteRequest
    module procedure LogDirectWriteRequest_int, LogDirectWriteRequest_str
  end interface
  
  logical, parameter :: countEmpty = .false.
  logical, parameter :: DEEBUG = .false.

contains ! ==================================================================

  ! --------------------------------------------- AccumulateSlaveArguments ------
  subroutine AccumulateSlaveArguments ( arg )
    ! This routine accumulates the command line arguments for the slaves
    character (len=*), intent(in) :: arg
    ! Executable code
    ! call NextPVMArg ( trim(arg) )
    if ( len_trim(slaveArguments) + len_trim(arg) +1 > len(slaveArguments) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Argument list for slave too long.' )
    slaveArguments = trim(slaveArguments)//' '//trim(arg)
  end subroutine AccumulateSlaveArguments
    
  ! ---------------------------------- AddMachineNameToDataBase ------------
  integer function AddMachineNameToDataBase ( DATABASE, ITEM )
    ! Add a machine name to a databsse of machine names
    ! an array of machine names

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    character(len=MachineNameLen), dimension(:), pointer :: DATABASE
    character(len=MachineNameLen), intent(in) :: ITEM
    ! Local variables
    character(len=MachineNameLen), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddMachineNameToDatabase = newSize
  end function AddMachineNameToDataBase

  ! ---------------------------------- AddMachineToDataBase ------------
  integer function AddMachineToDataBase ( DATABASE, ITEM )
    ! Add a machine to a database of machines

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    type(machine_t), dimension(:), pointer :: DATABASE
    type(machine_t), intent(in) :: ITEM
    ! Local variables
    type(machine_t), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddMachineToDatabase = newSize
  end function AddMachineToDataBase

  ! ---------------------------------------------- InitParallel -------------
  subroutine InitParallel ( chunkNo, MAFNo )
    ! This routine initialises the parallel code
    integer, intent(in) :: chunkNo      ! Chunk number asked to do.
    integer, intent(in) :: MAFNo        ! MAF number asked to do in fwmParallel mode.
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Executable code
    if ( parallel%master .or. parallel%slave ) then
      call output ( 'About to get my tid', advance='yes' )
      call PVMFMyTid ( parallel%myTid )
      call output ( 'Task ID: ' )
      call output ( trim(GetNiceTidString ( parallel%myTid ) ), advance='yes' )
    end if
    if ( parallel%slave ) then
      ! Register ourselves with the master
      ! Identify ourselves
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMF90Pack ( SIG_Register, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing registration' )
      call PVMF90Pack ( chunkNo, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing chunkNumber' )
      call PVMF90Pack ( MAFNo, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing MAFNumber' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending finish packet' )
      if ( parallel%fwmParallel ) then
        call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ parallel%masterTid /), info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'setting up notify' )
        call PVMFJoinGroup ( FWMSlaveGroup, info )
        if ( info < 0 ) &
          & call PVMErrorMessage ( info, 'Joining '//FWMSlaveGroup//' group' )
      end if
    end if
  end subroutine InitParallel

  ! --------------------------------------------- CloseParallel -------------
  subroutine CloseParallel ( noSlaves, dontWait )
    ! This routine closes down any parallel stuff

    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    integer, intent(in) :: noSlaves     ! How many slaves
    logical, intent(in), optional :: dontWait
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    integer :: Me = -1                  ! String index for trace

    integer :: SIGNAL                   ! From acknowledgement packet
    integer :: slave
    logical :: myDontWait

    ! Executable code

    call trace_begin ( me, "CloseParallel", cond=toggle(gen) )
    myDontWait = .false.
    if ( present(dontWait) ) myDontWait = dontWait
    if ( parallel%slave ) then
      myDontWait = .true.
      ! call MLSMessage ( MLSMSG_Info, ModuleName, &
      !  & 'About to call PVMFInitSend' )
      ! Send a request to finish
      call PVMFInitSend ( PvmDataDefault, bufferID )
      ! call MLSMessage ( MLSMSG_Info, ModuleName, &
      !  & 'About to call PVMF90Pack' )
      call PVMF90Pack ( SIG_Finished, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing finished signal' )
      ! call MLSMessage ( MLSMSG_Info, ModuleName, &
      !  & 'About to call PVMFSend' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending finish packet' )
      ! call MLSMessage ( MLSMSG_Error, ModuleName, &
      !  & 'Hey--we did it' )
      if ( myDontWait ) go to 9
      ! Wait for an acknowledgement from the master
      call PVMFrecv ( parallel%masterTid, InfoTag, bufferID )
      if ( bufferID <= 0 ) &
        & call PVMErrorMessage ( bufferID, &
        & 'receiving finish acknowledgement' )
      call PVMF90Unpack ( signal, info )
      if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking finish acknowledgement')
      if ( signal /= SIG_AckFinish ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Got unrecognised signal from master' )
      if ( DEEBUG ) call output ( 'Master acknowledged our completion', advance='yes' )
      if ( DEEBUG ) call Beep ( 'Master acknowledged our completion' )
    else if ( parallel%master ) then
      call usleep(FIXDELAYFORSLAVESTDOUTBUFFER)
      if ( noSlaves < 1 ) go to 9
      do slave=1, noSlaves+40
        call usleep(DELAYFOREACHSLAVESTDOUTBUFFER)
      end do
    end if
9   continue
    call trace_end ( "CloseParallel", cond=toggle(gen) )
  end subroutine CloseParallel

  ! ----------------------------------- CompactDirectWriteRequestDB --------
  subroutine CompactDirectWriteRequestDB ( database, noRequests )
    type ( DirectWriteRequest_T), dimension(:), pointer :: DATABASE
    integer, intent(out) :: NOREQUESTS
    ! Move all the non completed requests up the list
    noRequests = count ( database%status /= DW_Completed .and. &
      & database%status /= DW_Invalid )
    database ( 1 : noRequests ) = pack ( database, &
      & database%status /= DW_Completed .and. database%status /= DW_Invalid )
    database(noRequests+1:)%status = DW_Invalid
  end subroutine CompactDirectWriteRequestDB

  ! --------------------------------------- dump_machine_database ---------
  subroutine dump_machine_database(machines, Details)
    ! dump data type
    type(Machine_T), dimension(:), pointer :: machines
    integer, intent(in), optional :: DETAILS  ! 1: dump each machine, too
    ! Internal variables
    integer :: i
    integer :: myDetails
    integer, dimension(:), pointer :: master_tids
    character(len=MACHINENAMELEN), dimension(:), pointer :: master_names
    ! Executable
    if ( .not. associated(machines) ) then
      call output('machine database not associated', advance='yes')
      return
    endif
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ('Size of machine database: ', advance='no')
    call output (size(machines), advance='yes')
    call output ('Number of machines alive: ', advance='no')
    call output (count(machines%OK), advance='yes')
    call output ('Number of machines free: ', advance='no')
    call output (count(machines%free .and. machines%OK), advance='yes')
    if ( size(machines) < 1 ) return
    nullify ( master_tids, master_names )
    call Allocate_test ( master_tids, size(machines), 'machine_tids', moduleName )
    call Allocate_test ( master_names, size(machines), 'machine_names', moduleName )
    master_tids = 0
    call GetUniqueInts(machines%master_tid, master_tids, i, minValue=1)
    call output ('Number of unique master_ids: ', advance='no')
    call output (i, advance='yes')
    call GetUniqueStrings(machines%master_name, master_names, i, fillValue=' ')
    call output ('Number of unique master_names: ', advance='no')
    call output (i, advance='yes')
    call GetUniqueStrings(machines%master_date, master_names, i, fillValue=' ')
    call output ('Number of unique master_dates: ', advance='no')
    call output (i, advance='yes')
    call DeAllocate_test ( master_tids,  'machine_tids', moduleName )
    call DeAllocate_test ( master_names, 'machine_names', moduleName )
    if ( myDetails < 1 ) return
    do i = 1, size(machines)
      call dump_machine(machines(i))
    enddo
  end subroutine dump_machine_database

  ! --------------------------------------- dump_machine ---------
  subroutine dump_machine(machine)
    ! dump data type
    type(Machine_T), intent(in) :: machine
    call output('machine name: ', advance = 'no')
    call output(trim(machine%name), advance = 'yes')
    call output('master tid: ', advance = 'no')
    call output(machine%master_tid, advance = 'yes')
    call output('master name: ', advance = 'no')
    call output(trim(machine%master_name), advance = 'yes')
    call output('master date: ', advance = 'no')
    call output(trim(machine%master_date), advance = 'yes')
    call output('tid: ', advance = 'no')
    call output(machine%tid, advance = 'yes')
    call output('chunk: ', advance = 'no')
    call output(machine%chunk, advance = 'yes')
    call output('alive?: ', advance = 'no')
    call output(machine%OK, advance = 'yes')
    call output('free?: ', advance = 'no')
    call output(machine%free, advance = 'yes')
  end subroutine dump_machine
   
  ! --------------------------------------- DumpDirectWriteRequest ---------
  subroutine DumpDirectWriteRequest ( request )
    type(DirectWriteRequest_T), intent(in) :: REQUEST
    ! Executable code
    call output ( 'Chunk=' )
    call output ( request%chunk )
    call output ( ', machine=' )
    call output ( request%machine )
    call output ( ', tree node=' )
    call output ( request%node, advance='yes' )
    call output ( 'File index=' )
    call output ( request%fileIndex )
    call output ( ', ticket=' )
    call output ( request%ticket )
    call output ( ', value=' )
    call output ( request%value )
    call output ( ' Status: ' )
    select case ( request%status )
    case ( dw_invalid )
      call output ( 'invalid', advance='yes' )
    case ( dw_pending )
      call output ( 'pending', advance='yes' )
    case ( dw_inProgress )
      call output ( 'in Progress', advance='yes' )
    case ( dw_completed )
      call output ( 'completed', advance='yes' )
    case default
      call output ( request%status, advance='no' )
      call output ( ' (unrecognized)', advance='yes' )
    end select
    if ( request%whenFinished < 0. ) return
    call output ( 'grant delay: ' )
    call output ( request%whenGranted - request%whenMade, advance='yes' )
    call output ( 'writing time: ' )
    call output ( request%whenFinished - request%whenGranted, advance='yes' )
  end subroutine DumpDirectWriteRequest

  ! --------------------------------------- DumpAllDirectWriteRequests -----
  subroutine DumpAllDirectWriteRequests ( requests, statsOnly )
    type(DirectWriteRequest_T), intent(in), dimension(:) :: REQUESTS
    logical, optional, intent(in) :: statsOnly
    ! Internal variables
    logical :: mystatsOnly
    integer :: I, numRequests
    real :: min, max, mean, stddev, rms
    ! Executable code
    mystatsOnly = .false.
    if ( present(statsOnly) ) mystatsOnly = statsOnly
    numRequests = FindFirst(requests%status, DW_INVALID) - 1
    if ( numRequests < 1 ) then
      call output ( 'No valid directWrite requests among ' )
      call output ( size(requests) )
      call output ( ' direct write requests', advance='yes' )
      return
    endif
    call output ( 'Dumping ' )
    call output ( numRequests )
    call output ( ' direct write requests', advance='yes' )
    call allstats(requests(1:numRequests)%whenGranted &
      & - requests(1:numRequests)%whenMade, &
      & min=min, max=max, mean=mean, stddev=stddev, rms=rms)
    call output ( 'grant delay statistics: ', advance='yes' )
    call output ( 'min: ', advance='no' )
    call output ( min, advance='no' )
    call output ( '   max: ', advance='no' )
    call output ( max, advance='no' )
    call output ( '   mean: ', advance='no' )
    call output ( mean, advance='yes' )
    call allstats(requests(1:numRequests)%whenFinished &
      & - requests(1:numRequests)%whenGranted, &
      & min=min, max=max, mean=mean, stddev=stddev, rms=rms)
    call output ( 'writing time statistics: ', advance='yes' )
    call output ( 'min: ', advance='no' )
    call output ( min, advance='no' )
    call output ( '   max: ', advance='no' )
    call output ( max, advance='no' )
    call output ( '   mean: ', advance='no' )
    call output ( mean, advance='yes' )
    if ( mystatsOnly ) return
    ! Now dump each dw request (this could be   h u g e  ) !    
    do i = 1, numRequests ! size(requests)
      call output ( 'Request #' )
      call output ( i, advance='yes' )
      call dump ( requests(i) )
    end do
  end subroutine DumpAllDirectWriteRequests

  ! --------------------------------------- FinishedDirectWrite ------------
  subroutine FinishedDirectWrite ( ticket )
    integer, intent(in) :: TICKET       ! Ticket number
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Executable code
    if ( parallel%verbosity > 0 ) then
      call output ( "Sending finished on ticket " )
      call output ( ticket, advance='yes' )
    endif
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_DirectWriteFinished, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write finished flag" )
    call PVMF90Pack ( ticket, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write finished ticket" )

    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write finished packet" )
    
  end subroutine FinishedDirectWrite

  ! ---------------------------------------- GetMachineNames ------------
  subroutine GetMachineNames ( machineNames )
    ! This reads the parallel slave name file and returns
    ! an array of machine names
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES

    ! Local variables
    character(len=MachineNameLen) :: ARCH ! A line from the file
    integer :: DTID                     ! From PVMFConfig
    integer :: FIRST                    ! First machine in file to use
    integer :: FIRSTCOLONPOS            ! Position in string
    logical :: GOTFIRSTLAST             ! Got a range
    integer :: I                        ! Loop inductor
    integer :: INFO                     ! From PVMFConfig
    integer :: LAST                     ! Last machine in file to use
    character(len=MachineNameLen) :: LINE ! A line from the file
    integer :: LUN                      ! Logical unit number
    integer :: MACHINE                  ! Counter
    integer :: NARCH                    ! From PVMFConfig
    integer :: NOLINES                  ! Number of lines in file
    integer :: NOMACHINES               ! Array size
    character(len=MachineNameLen) :: ORIGINAL ! A line from the file
    integer :: SECONDCOLONPOS           ! Position in string
    integer :: SPEED                    ! From PVMFConfig
    integer :: STAT                     ! Status flag from read

    ! Executable code
    ! The slave filename may contain starting and ending line numbers

    original = parallel%slaveFilename

    if ( trim(original) == 'pvm' ) then
      ! If the user specifies 'pvm' then just get the names of the machines
      ! in the pvm system
      machine = 1
      hostLoop: do
        call PVMFConfig ( noMachines, narch, dtid, line, arch, speed, info )
        if ( info < 0 ) call PVMErrorMessage ( info, &
          & 'Calling PVMFConfig' )
        if ( machine == 1 ) then
          call Allocate_test ( machineNames, noMachines, 'machineNames', moduleName )
        end if
        machineNames(machine) = trim(line)
        machine = machine + 1
        if ( machine == noMachines + 1 ) exit hostLoop
      end do hostLoop
    else
      gotFirstLast = .false.
      firstColonPos = index ( parallel%slaveFilename, ':' )
      if ( firstColonPos /= 0 ) then
        gotFirstLast = .true.
        parallel%slaveFilename(firstColonPos:firstColonPos) = ' '
        secondColonPos = index ( parallel%slaveFilename,':' )
        if ( secondColonPos < firstColonPos ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Incorrect syntax for slave info:'//trim(original) )
        parallel%slaveFilename(secondColonPos:secondColonPos) = ' '
        read (parallel%slaveFilename(firstColonPos+1:),*,iostat=stat) first, last
        if ( stat /= 0 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Incorrect syntax for slave info:'//trim(original) )
      else
        firstColonPos = len_trim(parallel%slaveFilename)+1
      end if

      ! Find a free logical unit number
      ! lun = get_lun ()
      call get_lun( lun, msg=.false. )
      if ( lun < 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "No logical unit numbers available" )
      open ( unit=lun, file=parallel%slaveFilename(1:firstColonPos-1),&
        & status='old', form='formatted', &
        & access='sequential', iostat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
        & "Unable to open slave file " // parallel%slaveFilename(1:firstColonPos-1) )

      ! Now read the file and count the lines
      noMachines = 0
      noLines = 0
      firstTimeRound: do
        read ( unit=lun, fmt=*, iostat=stat ) line
        if ( stat < 0 ) exit firstTimeRound
        noLines = noLines + 1
        line = adjustl ( line )
        if ( line(1:1) /= '#' ) noMachines = noMachines + 1
      end do firstTimeRound

      ! Now setup the result array
      if ( .not. gotFirstLast ) then
        first = 1
        last = noMachines
      else
        first = min( max ( first, 1 ), noMachines )
        last = max( 1, min ( last, noMachines ) )
      end if
      noMachines = last - first + 1

      call Allocate_test ( machineNames, noMachines, 'machineNames', moduleName )
      ! Now rewind the file and read the names
      rewind ( lun )
      machine = 1
      do i = 1, noLines
        read ( unit=lun, fmt=* ) line
        line = adjustl ( line )
        if ( line(1:1) /= '#' ) then
          if ( (machine >= first) .and. (machine <= last) ) &
            & machineNames(machine-first+1) = line
          machine = machine + 1
        end if
      end do

      close ( unit=lun )
    end if
    if ( switchDetail(switches,'mach') /=0 ) &
      & call dump ( machineNames, options=what_options(trim=.true.) )
  end subroutine GetMachineNames

  ! ---------------------------------------- GetMachines ------------
  subroutine GetMachines ( machines )
    ! This reads the parallel slave name file and returns
    ! an array of machine types
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type(machine_t), dimension(:), pointer :: MACHINES

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: noMachines, stat
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    ! Executable
    nullify(machineNames, machines)
    call GetMachineNames(machineNames)
    noMachines = size(machineNames)
    if ( noMachines < 1 ) return
    allocate ( machines(noMachines), stat=stat )
    addr = 0
    if ( stat == 0 ) addr = transfer(c_loc(machines(1)), addr)
    call test_allocate ( stat, ModuleName, 'machines', ubounds=noMachines, &
      & elementSize = storage_size(machines) / 8, address=addr )
    machines%name = machineNames
  end subroutine GetMachines

  ! ------------------------------------ InflateDirectWriteRequestDB --
  integer function InflateDirectWriteRequestDB ( database, extra )
    ! Make the directWrite database bigger by extra
    ! return index of first new element

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type(DirectWriteRequest_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type(DirectWriteRequest_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateDirectWriteRequestDB = firstNewItem
  end function InflateDirectWriteRequestDB

  ! ---------------------------------------- LogDirectWriteRequest_int --
  subroutine LogDirectWriteRequest_int ( filename, node )
    integer, intent(in) :: FILENAME     ! String index of filename
    integer, intent(in) :: NODE         ! Node for directWrite line

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM

    ! Executable code

    ! Pack and dispatch
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_RequestDirectWrite, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request flag" )
    call PVMPackStringIndex ( filename, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request filename" )
    call PVMF90Pack ( node, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request node" )

    if ( DELAYFOREACHSLAVEDWREQUEST > 0 ) &
      & call usleep ( DELAYFOREACHSLAVEDWREQUEST )
    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write request packet" )

  end subroutine LogDirectWriteRequest_int

  ! ---------------------------------------- LogDirectWriteRequest_str --
  subroutine LogDirectWriteRequest_str ( filename, node )
    character(len=*), intent(in) :: FILENAME     ! String filename
    integer, intent(in) :: NODE         ! Node for directWrite line

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM

    ! Executable code

    ! Pack and dispatch
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_RequestDirectWrite, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request flag" )
    ! call PVMPackStringIndex ( filename, info )
    call PVMIDLPack ( trim(filename), info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request filename" )
    call PVMF90Pack ( node, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request node" )

    if ( DELAYFOREACHSLAVEDWREQUEST > 0 ) &
      & call usleep ( DELAYFOREACHSLAVEDWREQUEST )
    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write request packet" )

  end subroutine LogDirectWriteRequest_str

  ! -------------------------------------------- WaitForDirectWritePermission --
  subroutine WaitForDirectWritePermission ( node, ticket, theFile, createFile )
    integer, intent(out) :: NODE        ! Which line was granted
    integer, intent(out) :: TICKET      ! What is the ticket number
    ! integer, intent(out) :: THEFILE   ! What is the file name
    character(len=*), intent(out), optional :: THEFILE  ! file name
    logical, intent(out) :: CREATEFILE  ! Do we have to create the file?
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    integer :: SIGNAL                   ! The signal from the master
    integer :: I4(4)                    ! Information from master
    ! Executable code
    call PVMFRecv ( parallel%masterTid, InfoTag, bufferID )
    if ( bufferID <= 0 ) call PVMErrorMessage ( bufferID, &
      & 'receiving direct write permission/wait' )

    ! Once we have the reply unpack it
    call PVMF90Unpack ( signal, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking direct write permission signal')

    if ( signal == SIG_DirectWriteGranted ) then
      call PVMF90Unpack ( I4, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking direct write permission information')
      call PVMUnpackStringIndex ( I4(4), info, outWord=theFile )
      if ( info /= 0 )  call PVMErrorMessage ( info, &
        & "unpacking direct write request filename" )
      node = I4(1)
      ticket = I4(2)
      createFile = I4(3) /= 0
      ! theFile = I4(4)

    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Got unrecognised signal from master' )
    end if

  end subroutine WaitForDirectWritePermission

  ! ------------------------------------------- SlaveJoin ---------------
  subroutine SlaveJoin ( quantity, precisionQuantity, hdfName, key )
    ! This simply sends one or more vector quantities down a pvm spigot.
    type (VectorValue_T), pointer :: QUANTITY    ! Quantity to join
    type (VectorValue_T), pointer :: PRECISIONQUANTITY ! Its precision
    character(len=*), intent(in) :: HDFNAME ! Swath / sd name
    integer, intent(in) :: KEY          ! Tree node

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: GOTPRECISION             ! really boolean
    integer :: INFO                     ! Flag from PVM

    ! Executable code
    gotPrecision = 0
    if ( associated ( precisionQuantity) ) gotPrecision = 1

    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_ToJoin, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing kind" )
    call PVMF90Pack ( (/key, gotPrecision /), info )
    if ( info /= 0 ) &
       & call PVMErrorMessage ( info, "packing key, gotPrecision" )
    call PVMIDLPack ( hdfName, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing hdfName" )

    call PVMSendQuantity ( quantity, justPack=.true., noMask=.true. )
    if ( associated ( precisionQuantity ) ) &
      call PVMSendQuantity ( precisionQuantity, justPack=.true., noMask=.true. )
    
    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage( info, 'sending join packet' )

  end subroutine SlaveJoin

  ! --------------------------------------------- SnipLastSlaveArgument ------
  subroutine SnipLastSlaveArgument
    ! This routine snips the last command line arguments for the slaves
    integer :: n
    character ( len=len(slaveArguments) ), save :: tempArguments = ""
    character(len=1), parameter :: SPACE = ' '
    ! Executable code
    if ( DEEBUG ) call output( 'Planning to SnipLastSlaveArgument from ' // trim(slaveArguments), &
      & advance='yes' )
    ! How many elements?
    n = NumStringElements( slaveArguments, countEmpty, inseparator=space )
    if ( DEEBUG ) call outputNamedValue ( 'n', n )
    if ( n < 1 ) return
    tempArguments = slaveArguments
    if ( DEEBUG ) call output( 'Before snip: ' // trim(slaveArguments), advance='yes' )
    call RemoveNumFromList( tempArguments, slaveArguments, n, &
      & inseparator=space )
    if ( DEEBUG ) call output( 'After snip: ' // trim(slaveArguments), advance='yes' )
  end subroutine SnipLastSlaveArgument
    
  ! --------------------------------------------- SigToName ------
  function SigToName( Sig ) result ( name )
    integer, intent(in)           :: Sig
    character(len=32)             :: name
    ! Executable
    select case ( Sig )
    case ( SIG_AboutToDie )
      name = 'AboutToDie'
    case ( SIG_ToJoin )
      name = 'ToJoin'
    case ( SIG_Finished )
      name = 'Finished'
    case ( SIG_AckFinish )
      name = 'AckFinish'
    case ( SIG_Register )
      name = 'Register'
    case ( SIG_RequestDirectWrite )
      name = 'RequestDirectWrite'
    case ( SIG_DirectWriteGranted )
      name = 'DirectWriteGranted'
    case ( SIG_DirectWriteFinished )
      name = 'DirectWriteFinished'
    case ( SIG_NewSetUp )
      name = 'NewSetUp'
    case ( SIG_RunMAF )
      name = 'RunMAF'
    case ( SIG_SendResults )
      name = 'SendResults'
    case ( SIG_RequestHost )
      name = 'RequestHost'
    case ( SIG_ReleaseHost )
      name = 'ReleaseHost'
    case ( SIG_HostDied )
      name = 'HostDied'
    case ( SIG_ThanksHost )
      name = 'ThanksHost'
    case ( SIG_SwearAllegiance )
      name = 'SwearAllegiance'
    case ( SIG_SwitchAllegiance )
      name = 'SwitchAllegiance'
    case default
      name = 'Unknown'
    end select
  end function SigToName
    
  ! --------------------------------------------- TagToName ------
  function TagToName( tag ) result ( name )
    integer, intent(in)           :: tag
    character(len=24)             :: name
    ! Executable
    select case ( tag )
    case ( Infotag )
      name = 'Info'
    case ( Chunktag )
      name = 'Chunk'
    case ( Notifytag )
      name = 'Notify'
    case ( Petitiontag )
      name = 'Petition'
    case ( Grantedtag )
      name = 'Granted'
    case ( MasterDumptag )
      name = 'MasterDump'
    case ( MachineFixedtag )
      name = 'MachineFixed'
    case ( GiveUptag )
      name = 'GiveUp'
    case default
      name = 'Unknown'
    end select
  end function TagToName
    
  ! --------------------------------------------- TransmitSlaveArguments ------
  subroutine TransmitSlaveArguments
    ! This routine transmits command line arguments to the slaves via pvm
    integer :: i, n
    character(len=1), parameter :: SPACE = ' '
    ! Executable code
    if ( parallel%verbosity > 0 ) call output( 'slave arguments: ' // &
      & trim(slaveArguments), advance='yes' )
    ! How many elements?
    n = NumStringElements( slaveArguments, countEmpty, inseparator=space )
    do i=1, n
      call NextPVMArg ( trim(StringElement ( slaveArguments, &
        & i, countEmpty, inseparator=space ) ) )
    enddo
  end subroutine TransmitSlaveArguments
    
  ! ----------------------------------------- GetNiceTidString -----
  character(len=16) function GetNiceTidString ( tid )
    integer, intent(in) :: tid

    ! Executable code
    write ( GetNiceTidString, '(z0)' ) tid
    GetNiceTidString = '[t'//trim(LowerCase ( GetNiceTidString ))//']'
  end function GetNiceTidString

  function what_options( clean, transpose, trim ) result( options )
    use MLSStrings, only: Trim_Safe
    logical, optional, intent(in) :: clean
    logical, optional, intent(in) :: transpose
    logical, optional, intent(in) :: trim
    character(len=8) :: options
    options = ' '
    if ( present(clean) ) then
      if ( clean ) options = trim_safe(options) // 'c'
    endif
    if ( present(transpose) ) then
      if ( transpose ) options = trim_safe(options) // 'p'
    endif
    if ( present(trim) ) then
      if ( trim ) options = trim_safe(options) // 't'
    endif
  end function what_options

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2ParInfo

! $Log$
! Revision 2.70  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.69  2018/02/09 01:07:07  pwagner
! Added functions to convert tags and sigs to corresponding names
!
! Revision 2.68  2017/12/22 00:31:05  pwagner
! Made Use statements CamelCase
!
! Revision 2.67  2016/03/23 20:12:42  pwagner
! Usleep not an external any more
!
! Revision 2.66  2016/02/29 19:50:26  pwagner
! Usleep got from machine module instead of being an external
!
! Revision 2.65  2015/03/28 02:48:02  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.64  2014/09/05 01:05:03  vsnyder
! Add some tracing.  Move some USEs from module to procedure scope.
!
! Revision 2.63  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.62  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.61  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.60  2013/02/14 19:03:32  pwagner
! Added way for l2q to tell master to dump status
!
! Revision 2.59  2012/07/10 15:19:45  pwagner
! Adapted to new api for RemoveNumFromList
!
! Revision 2.58  2012/07/02 20:39:59  pwagner
! Once-routine output now requires verbosity > 0
!
! Revision 2.57  2012/06/27 18:10:41  pwagner
! Added TransmitSlaveArguments, SnipLastSlaveArgument
!
! Revision 2.56  2012/03/28 20:05:10  pwagner
! Removed staging file--slave tasks lost ability to join quantities
!
! Revision 2.55  2011/05/09 18:20:55  pwagner
! Converted to using switchDetail
!
! Revision 2.54  2010/08/06 23:08:48  pwagner
! Pass Hessians, matrices to DumpCommand
!
! Revision 2.53  2009/07/21 20:35:38  pwagner
! Nullify pointers before call ing allocate_test with them
!
! Revision 2.52  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.51  2009/06/16 17:41:37  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.50  2007/10/24 00:16:38  pwagner
! Removed unused declarations
!
! Revision 2.49  2007/04/03 20:55:22  pwagner
! Uses get_lun from lib/io_stuff instead of internal function
!
! Revision 2.48  2007/01/12 00:38:57  pwagner
! Switches allegiance to a replacement l2q
!
! Revision 2.47  2006/06/24 23:10:00  pwagner
! Remove unneeded thing from output_m
!
! Revision 2.46  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.45  2005/03/24 21:20:33  pwagner
! New fileds in directWriteRequests figure grant delay, writing time
!
! Revision 2.44  2005/03/15 23:53:03  pwagner
! PVMERRORMESSAGE now part of MLSMessageModule; INFOTag, SIG_ABOUTTODIE from lib/PVM
!
! Revision 2.43  2005/02/03 19:07:35  pwagner
! master_name and master_date fields added to machine type
!
! Revision 2.42  2005/01/13 00:01:31  pwagner
! Dumping machineDB shows only free machines that are still alive
!
! Revision 2.41  2005/01/07 17:27:21  pwagner
! Details an optional arg to dump_machine_database to bypass dumping each machine
!
! Revision 2.40  2004/12/14 21:53:30  pwagner
! Added machine_T, tags and signals related to l2q
!
! Revision 2.39  2004/09/16 23:57:32  pwagner
! Now tracks machine names of failed chunks
!
! Revision 2.38  2004/09/16 00:18:03  pwagner
! Keeps record of completed, failed chunks
!
! Revision 2.37  2004/08/05 22:47:47  pwagner
! New --chunkRange option to run selected chunks in parallel mode
!
! Revision 2.36  2004/01/22 00:56:35  pwagner
! Fixed many bugs in auto-distribution of DirectWrites
!
! Revision 2.35  2004/01/02 23:36:00  pwagner
! DirectWrites may choose files automatically from db
!
! Revision 2.34  2003/12/11 23:00:58  pwagner
! Make master task wait for slaves stdout buffers to flush
!
! Revision 2.33  2003/11/14 23:37:13  pwagner
! Lets user change masterLoop delay via commandline option
!
! Revision 2.32  2003/08/11 23:23:31  pwagner
! ChunkNo component added to L2ParallelInfo_T
!
! Revision 2.31  2003/08/01 20:26:53  pwagner
! slave command name saved as component pgeName
!
! Revision 2.30  2003/07/07 17:32:10  livesey
! New approach to DirectWrite
!
! Revision 2.29  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.28  2003/05/12 19:04:50  livesey
! Added the staging file to L2ParallelInfo_T
!
! Revision 2.27  2003/01/27 17:18:58  livesey
! Made master task output it's tid as well as slaves.
!
! Revision 2.26  2003/01/17 21:54:21  livesey
! Added the machineFixed stuff
!
! Revision 2.25  2002/12/11 01:58:33  livesey
! Added noFWMSlaves
!
! Revision 2.24  2002/11/23 00:13:18  livesey
! Increased maxFailuresPerChunk to reflect nature of Raytheon cluster (at
! least for the moment).
!
! Revision 2.23  2002/10/17 18:19:13  livesey
! Added GiveupTag
!
! Revision 2.22  2002/10/08 20:33:51  livesey
! Added notify stuff for FWMParallel
!
! Revision 2.21  2002/10/08 17:40:56  livesey
! Minor bug fix in FWMParallel stuff
!
! Revision 2.20  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.19  2002/10/06 22:22:47  livesey
! Removed MAFTAG and MAF communication stuff
!
! Revision 2.18  2002/10/06 01:10:17  livesey
! More progress on the fwmParallel stuff
!
! Revision 2.17  2002/10/05 00:43:20  livesey
! Started work on the fwmParallel stuff
!
! Revision 2.16  2002/07/19 06:07:32  livesey
! Cut down MaxFailuresPerChunk
!
! Revision 2.15  2002/07/17 20:02:17  livesey
! Bug fix
!
! Revision 2.14  2002/07/17 20:01:42  livesey
! Put an advance=yes in there
!
! Revision 2.13  2002/07/17 19:54:52  livesey
! Made slaves at least identify themselves
!
! Revision 2.12  2002/05/22 00:48:28  livesey
! Added direct write stuff
!
! Revision 2.11  2002/04/24 16:53:50  livesey
! Changes to implement submit.
!
! Revision 2.10  2002/03/21 01:23:36  livesey
! Changed thresholds
!
! Revision 2.9  2001/10/30 01:45:21  livesey
! Some modifications/fixes to parallel join
!
! Revision 2.8  2001/09/08 00:21:44  pwagner
! Revised to work for new column Abundance in lone swaths
!
! Revision 2.7  2001/05/30 23:53:54  livesey
! Vectors now sent within info packets
!
! Revision 2.6  2001/05/29 20:42:18  pwagner
! Added Log at bottom of source
!
