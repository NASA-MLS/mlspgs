! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2ParInfo
  ! This module provides definitions needed by L2Parallel and other modules to
  ! manage the parallel aspects of the L2 code.

  use Allocate_Deallocate, only: ALLOCATE_TEST
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_Allocate, &
    & MLSMSG_Deallocate
  use PVM, only: PVMFMYTID, PVMFINITSEND, PVMF90PACK, PVMFSEND, &
    & PVMDATADEFAULT, PVMERRORMESSAGE, PVMF90UNPACK, NEXTPVMARG, PVMTASKEXIT
  use PVMIDL, only: PVMIDLPACK
  use MorePVM, only: PVMPACKSTRINGINDEX, PVMUNPACKSTRINGINDEX
  use VectorsModule, only: VECTORVALUE_T
  use QuantityPVM, only: PVMSENDQUANTITY
  use MLSStrings, only: LowerCase
  use Output_M, only: Output
  use Toggles, only: SWITCHES

  implicit none
  private

  public :: L2ParallelInfo_T, parallel, InitParallel, CloseParallel
  public :: SIG_ToJoin, SIG_Finished, SIG_Register, ChunkTag, InfoTag, SlaveJoin
  public :: SIG_AckFinish, SIG_RequestDirectWrite, SIG_DirectWriteGranted
  public :: SIG_DirectWriteFinished
  public :: SIG_NewSetup, SIG_RunMAF, SIG_SendResults
  public :: NotifyTag, GetNiceTidString, GiveUpTag, SlaveArguments
  public :: AccumulateSlaveArguments, LogDirectWriteRequest
  public :: FinishedDirectWrite, MachineNameLen, GetMachineNames
  public :: FWMSlaveGroup, MachineFixedTag, DirectWriteRequest_T
  public :: DW_Invalid, DW_Pending, DW_InProgress, DW_Completed
  public :: InflateDirectWriteRequestDB, WaitForDirectWritePermission
  public :: CompactDirectWriteRequestDB, Dump
  
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! Parameters

  integer, parameter :: DELAYFOREACHSLAVESTDOUTBUFFER   = 500
  integer, parameter :: FIXDELAYFORSLAVESTDOUTBUFFER   = 500000
  integer, parameter :: CHUNKTAG   = 10
  integer, parameter :: INFOTAG    = ChunkTag + 1
  integer, parameter :: NOTIFYTAG  = InfoTag + 1
  integer, parameter :: MACHINEFIXEDTAG = 800
  integer, parameter :: GIVEUPTAG  = 999

  integer, parameter :: SIG_TOJOIN = 1
  integer, parameter :: SIG_FINISHED = SIG_toJoin + 1
  integer, parameter :: SIG_ACKFINISH = SIG_finished + 1
  integer, parameter :: SIG_REGISTER = SIG_AckFinish + 1
  integer, parameter :: SIG_REQUESTDIRECTWRITE = SIG_Register + 1
  integer, parameter :: SIG_DIRECTWRITEGRANTED = SIG_RequestDirectWrite + 1
  integer, parameter :: SIG_DIRECTWRITEFINISHED = SIG_DirectWriteGranted + 1
  integer, parameter :: SIG_NEWSETUP = SIG_DirectWriteFinished + 1
  integer, parameter :: SIG_RUNMAF = SIG_NewSetup + 1
  integer, parameter :: SIG_SENDRESULTS = SIG_RunMAF + 1

  integer, parameter :: MACHINENAMELEN = 64 ! Max length of name of machine

  ! Name of fwm slave group
  character (len=*), parameter :: FWMSLAVEGROUP = "MLSL2FWMSlaves"
  character (len=*), parameter :: DEFAULTSTAGINGFILE = &
    & "MLS-Aura_L2Staging-Full_v0-0-0_0000d000.h5"

  ! This datatype defines configuration for the parallel code
  type L2ParallelInfo_T
    logical :: fwmParallel = .false.    ! Set if we are in forward model parallel mode
    logical :: master = .false.         ! Set if this is a master task
    logical :: slave = .false.          ! Set if this is a slave task
    logical :: stageInMemory = .false.  ! Set if master stages to memory rather
    integer :: myTid                    ! My task ID in pvm       | than a file
    integer :: masterTid                ! task ID in pvm
    integer :: noFWMSlaves              ! No. slaves in pvm system for fwm cases
    integer :: chunkNo                  ! Chunk No. if this is a slave
    integer :: Delay = 200000           ! For Usleep, no. microsecs
    character(len=132) :: pgeName="mlsl2"   ! command name, if not 'mlsl2'
    character(len=132) :: slaveFilename ! Filename with list of slaves
    character(len=132) :: executable    ! Executable filename
    character(len=132) :: submit=""     ! Submit comand for batch queue system
    character(len=132) :: stagingFile=DefaultStagingFile ! Filename for possible
    integer :: maxFailuresPerMachine = 1 ! More than this then don't use it | staging
    integer :: maxFailuresPerChunk = 10 ! More than this then give up on getting it
  end type L2ParallelInfo_T

  ! This enumerated type describes the state that directWrites can be in
  integer, parameter :: DW_INVALID = 0
  integer, parameter :: DW_PENDING = DW_INVALID + 1
  integer, parameter :: DW_INPROGRESS = DW_PENDING + 1
  integer, parameter :: DW_COMPLETED = DW_INPROGRESS + 1

  ! This datatype logs a directWrite request
  type DirectWriteRequest_T
    integer :: CHUNK=0                  ! Which chunk made the request
    integer :: MACHINE=0                ! What machine is that on
    integer :: NODE=0                   ! What tree node was it
    integer :: FILEINDEX=0              ! Which file was it
    integer :: TICKET=0                 ! What ticket number is it
    integer :: VALUE=0                  ! Workspace for master
    integer :: STATUS=DW_INVALID        ! One of the DW_... above
    real    :: WHENFIRSTPERMITTED
  end type DirectWriteRequest_T

  ! Shared variables
  type (L2ParallelInfo_T), save :: parallel
  character ( len=2048 ), save :: slaveArguments = ""

  interface dump
    module procedure DumpDirectWriteRequest, DumpAllDirectWriteRequests
  end interface

contains ! ==================================================================

  ! --------------------------------------------- AccumulateSlaveArguments ------
  subroutine AccumulateSlaveArguments ( arg )
    ! This routine accumulates the command line arguments for the slaves
    character (len=*), intent(in) :: arg
    ! Executable code
    call NextPVMArg ( arg )
    if ( len_trim(slaveArguments) + len_trim(arg) +1 > len(slaveArguments) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Argument list for slave too long.' )
    slaveArguments = trim(slaveArguments)//' '//trim(arg)
  end subroutine AccumulateSlaveArguments
    
  ! ---------------------------------------------- InitParallel -------------
  subroutine InitParallel ( chunkNo, MAFNo )
    use Output_m, only: PRUNIT
    ! This routine initialises the parallel code
    integer, intent(in) :: chunkNo      ! Chunk number asked to do.
    integer, intent(in) :: MAFNo        ! MAF number asked to do in fwmParallel mode.
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Executable code
    if ( parallel%master .or. parallel%slave ) then
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
  subroutine CloseParallel(noSlaves)
    ! This routine closes down any parallel stuff
    integer, intent(in) :: noSlaves     ! How many slaves
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM

    integer :: SIGNAL                   ! From acknowledgement packet
    integer :: slave

    ! Executable code
    if ( parallel%slave ) then
      ! Send a request to finish
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMF90Pack ( SIG_Finished, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing finished signal' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending finish packet' )

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
    elseif ( parallel%master ) then
      call usleep(FIXDELAYFORSLAVESTDOUTBUFFER)
      if ( noSlaves < 1 ) return
      do slave=1, noSlaves
        call usleep(DELAYFOREACHSLAVESTDOUTBUFFER)
      enddo
    end if
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
    end select
  end subroutine DumpDirectWriteRequest

  ! --------------------------------------- DumpAllDirectWriteRequests -----
  subroutine DumpAllDirectWriteRequests ( requests )
    type(DirectWriteRequest_T), intent(in), dimension(:) :: REQUESTS
    integer :: I
    ! Executable code
    call output ( 'Dumping ' )
    call output ( size(requests) )
    call output ( ' direct write requests', advance='yes' )
    do i = 1, size(requests)
      call output ( 'Request #' )
      call output ( i, advance='yes' )
      call dump ( requests(i) )
    end do
  end subroutine DumpAllDirectWriteRequests

  ! --------------------------------------- FinishedDirectWrite ------------
  subroutine FinishedDirectWrite ( ticket )
    integer, intent(in) :: TICKET       ! Ticket number
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Local variables
    call output ( "Sending finished on ticket " )
    call output ( ticket, advance='yes' )
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
    
    ! Executable code
  end subroutine FinishedDirectWrite

  ! ---------------------------------------- GetMachineNames ------------
  subroutine GetMachineNames ( machineNames )
    ! This reads the parallel slave name file and returns
    ! an array of machine names
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES

    ! Local variables
    character(len=MachineNameLen) :: LINE ! A line from the file
    character(len=MachineNameLen) :: ORIGINAL ! A line from the file
    character(len=MachineNameLen) :: ARCH ! A line from the file
    logical :: EXIST                    ! Flag from inquire
    logical :: GOTFIRSTLAST             ! Got a range
    logical :: OPENED                   ! Flag from inquire

    integer :: DTID                     ! From PVMFConfig
    integer :: FIRST                    ! First machine in file to use
    integer :: FIRSTCOLONPOS            ! Position in string
    integer :: I                        ! Loop inductor
    integer :: INFO                     ! From PVMFConfig
    integer :: LAST                     ! Last machine in file to use
    integer :: LUN                      ! Logical unit number
    integer :: MACHINE                  ! Counter
    integer :: NARCH                    ! From PVMFConfig
    integer :: NOLINES                  ! Number of lines in file
    integer :: NOMACHINES               ! Array size
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
      lun = get_lun ()
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

      allocate ( machineNames(noMachines), stat=stat )
      if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate//' machineNames')

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
        endif
      end do

      close ( unit=lun )
    end if
  end subroutine GetMachineNames

  ! ------------------------------------ InflateDirectWriteRequestDB --
  integer function InflateDirectWriteRequestDB ( database, extra )
    ! Make the directWrite database bigger by extra
    ! return index of first new element

    ! Dummy arguments
    type(DirectWriteRequest_T), dimension(:), pointer :: DATABASE
    integer, intent(in) :: EXTRA

    ! Local variables
    type(DirectWriteRequest_T), dimension(:), pointer :: TEMPDATABASE

    include "inflateDatabase.f9h"
    InflateDirectWriteRequestDB = firstNewItem
  end function InflateDirectWriteRequestDB

  ! ---------------------------------------- LogDirectWriteRequest --
  subroutine LogDirectWriteRequest ( filename, node )
    integer, intent(in) :: FILENAME     ! String index of filename
    integer, intent(in) :: NODE         ! Node for directWrite line

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    integer :: SIGNAL                   ! From Master

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

    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write request packet" )

  end subroutine LogDirectWriteRequest

  ! -------------------------------------------- WaitForDirectWritePermission --
  subroutine WaitForDirectWritePermission ( node, ticket, theFile, createFile )
    integer, intent(out) :: NODE        ! Which line was granted
    integer, intent(out) :: TICKET      ! What is the ticket number
    integer, intent(out) :: THEFILE     ! What is the file name
    logical, intent(out) :: CREATEFILE  ! Do we have to create the file?
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    integer :: SIGNAL                   ! The signal from the master
    integer :: CREATE                   ! Integer version of createFile
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
      node = I4(1)
      ticket = I4(2)
      createFile = I4(3) /= 0
      theFile = I4(4)

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

  ! ----------------------------------------- GetNiceTidString -----
  character(len=16) function GetNiceTidString ( tid )
    integer, intent(in) :: tid

    ! Executable code
    write ( GetNiceTidString, '(z0)' ) tid
    GetNiceTidString = '[t'//trim(LowerCase ( GetNiceTidString ))//']'
  end function GetNiceTidString

  ! --------------------------------------- get_lun -----
  integer function get_lun ()
    ! Local variables
    integer :: LUN
    logical :: EXIST
    logical :: OPENED
    ! Executable code
    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    get_lun = lun
  end function get_lun

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module L2ParInfo

! $Log$
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
