! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module L2Parallel
  ! This module contains low level routines and stuff for dealing with parallel
  ! invocations of the MLSL2 program.

  ! Level 2 programs can be masters or slaves, or neither, not both.  A task
  ! which is neither simply runs through the l2cf as normal.  A master task
  ! does the chunk divide and then launches slave tasks for each chunk, and
  ! awaits their results in the Join section.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use BitStuff, only: BitsToBooleans, BooleansToBits
  use Chunks_M, only: Dump, MLSChunk_T
  use ChunkDivide_M, only: ChunkDivideConfig
  use Dump_0, only: Dump
  use HighOutput, only: BeVerbose, Banner, OutputNamedValue, TimeStamp
  use L2ParInfo, only: Machine_T, Parallel, &
    & Chunktag, Giveuptag, Grantedtag, Notifytag, Masterdumptag, Petitiontag, &
    & Sig_Tojoin, Sig_Finished, Sig_Ackfinish, Sig_Register, &
    & Sig_Requestdirectwrite, Sig_Swearallegiance, Sig_Switchallegiance, &
    & Sig_Directwritegranted, Sig_Directwritefinished, &
    & Sig_Hostdied, Sig_Releasehost, Sig_Requesthost, Sig_Thankshost, &
    & Getnicetidstring, Slavearguments, Machinenamelen, Getmachines, &
    & Machinefixedtag, Directwriterequest_T, &
    & Dw_Pending, Dw_Inprogress, Dw_Completed, &
    & Inflatedirectwriterequestdb, Compactdirectwriterequestdb, Dump, &
    & AddmachinetoDatabase, SigToName
  use Machine, only: Shell_Command, Usleep
  use MLSKinds, only: R8
  use MLSL2options, only: L2Options, OriginalOptions, MLSL2Message
  use MLSMessagemodule, only: MLSMSG_Error, MLSMSG_Warning, PVMErrorMessage
  use MLSFinds, only: Findall, Findfirst
  use MLSStringlists, only: Catlists, Expandstringrange, Removenumfromlist, &
    & Replacesubstring, Switchdetail
  use MLSStrings, only: Lowercase
  use MorePVM, only: PVMunpackstringindex, PVMpackstringindex
  use MLSStrings, only: NAppearances
  use Output_M, only: Blanks, Output
  use PVM, only: Infotag, &
    & PVMDatadefault, PVMFinitSend, PVMf90Pack, PVMFKill, &
    & PVMF90unpack, PVMtaskhost, &
    & MyPVMSpawn, PVMFCatchout, PVMFSend, PVMFNotify, PVMTaskExit, &
    & GetMachineNameFromTID, PVMFFreeBuf, Sig_Abouttodie
  use PVMidl, only: PVMidlunpack
  use String_Table, only: Display_String
  use Time_M, only: Time_Now
  use Toggles, only: Gen, Switches, Toggle
  use Trace_M, only: Trace_Begin, Trace_End
  use WriteMetaData, only: L2pcf

  implicit none
  private

  public :: L2MasterTask
  public :: GetChunkinfoFromMaster
  
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Parameters
  integer, parameter :: HDFNAMELEN = 132 ! Max length of name of swath/sd
  character(len=*), parameter :: &
    &                   GROUPNAME = "mlsl2" ! Set by l2q
  logical, parameter :: NOTFORGOTTEN = .false. ! Note the death of forgottens
  integer, save      :: FIXDELAYFORSLAVETOFINISH   = 1500000 ! 15000000 ! 15 sec
  integer, save      :: KILLINGSLAVESDELAY   = 1000000 ! 1 sec
  ! ***********************************************************************
  ! The following must be longer than the ???? setting in
  ! slavettmplt.sh
  integer, save      :: WaitBeforeKillingSlaves   = 130*1000000 ! 130 sec
  !                      Note: This must be > pgekilldelay in slavetmplt.sh
  ! ***********************************************************************

  ! Private types
  type StoredResult_T
    integer :: key                      ! Tree node for this join
    logical :: gotPrecision             ! If set have precision as well as value
    integer, dimension(:), pointer :: valInds=>NULL() ! Array vec. dtbs. inds (noChunks)
    integer, dimension(:), pointer :: precInds=>NULL() ! Array vec. dtbs. inds (noChunks)
    character(len=HDFNameLen) :: hdfName ! Name of swath/sd
  end type StoredResult_T

contains 
! ================================ Public Procedures ======================
  ! ----------------------------------------------- GetChunkInfoFromMaster ------
  subroutine GetChunkInfoFromMaster ( chunks, chunkNo )
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    ! This function gets a chunk sent by SendChunkToSlave

    ! Dummy arguments and result
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS
    integer, intent(out) :: CHUNKNO

    ! Local parameters
    integer, parameter :: noChunkTerms = 5 ! Number of components in MLSChunk_T

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: BITS
    integer :: BUFFERID                 ! ID for buffer for receive
    integer :: CHUNK                    ! Loop counter
    integer :: INFO                     ! Flag from PVM
    integer :: STATUS                   ! From allocate
    integer :: NOCHUNKS                 ! Size.
    integer :: NOHGRIDS                 ! Size of hGrid information
    integer, dimension(2) :: HEADER     ! No chunks, chunkNo
    logical, dimension(2) :: LOGICALS   ! Allow overlaps outside proc. rnge?
    integer, dimension(noChunkTerms) :: VALUES ! Chunk as integer array

    ! Executable code
    if ( .not. parallel%slave ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Only slave tasks can receive a chunk' )

    call PVMFrecv ( parallel%masterTid, ChunkTAG, bufferID )
    if ( bufferID <= 0 ) &
      & call PVMErrorMessage ( info, 'receiving chunkInfo' )

    call PVMF90Unpack ( header, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking chunkInfo header')
    noChunks = header(1)
    chunkNo = header(2)
    
    allocate ( chunks ( noChunks ), STAT=status)
    addr = 0
    if ( status == 0 .and. noChunks>0 ) addr = transfer(c_loc(chunks(1)), addr)
    call test_allocate ( status, moduleName, 'chunks', uBounds = noChunks, &
      & elementSize = storage_size(chunks) / 8, address=addr )

    do chunk = 1, noChunks
      call PVMF90Unpack ( values, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking one chunk')
      chunks(chunk)%firstMAFIndex = values(1)
      chunks(chunk)%lastMAFIndex = values(2)
      chunks(chunk)%noMAFsLowerOverlap = values(3)
      chunks(chunk)%noMAFsUpperOverlap = values(4)
      chunks(chunk)%chunkNumber = values(5)

      call PVMF90Unpack ( noHGrids, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking noHGrids')
      call Allocate_test ( chunks(chunk)%hGridOffsets, noHGrids, &
        & 'hGridOffsets', ModuleName )
      call PVMF90Unpack ( chunks(chunk)%hGridOffsets, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking hGridOffsets')

      call PVMF90Unpack ( noHGrids, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking noHGrids')
      call Allocate_test ( chunks(chunk)%hGridTotals, noHGrids, &
        & 'hGridTotals', ModuleName )
      call PVMF90Unpack ( chunks(chunk)%hGridTotals, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking hGridTotals')

      ! The latest components of the MLSChunk_T
      call PVMF90UnPack ( chunks(chunk)%StartTime, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking Start time of chunk' )
      call PVMF90UnPack ( chunks(chunk)%EndTime, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'unpacking End time of chunk' )
    end do

    call PVMF90Unpack ( bits, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking overlaps outside proc. rnge')
    call BitsToBooleans( Bits, logicals )
    ChunkDivideConfig%allowPriorOverlaps = logicals(1)
    ChunkDivideConfig%allowPostOverlaps  = logicals(2)
    if( BeVerbose( 'opt', -1 ) .or. parallel%verbosity > 0 ) then
      call output(' (chunk info received from master) ', advance='yes')
      call output(' Chunk number                       :            ', advance='no')
      call blanks(4, advance='no')
      call output(chunkno, advance='yes')
      call output(' Allow overlaps outside proc. range?:            ', advance='no')
      call blanks(4, advance='no')
      call output(ChunkDivideConfig%allowPriorOverlaps, advance='no')
      call blanks(4, advance='no')
      call output(ChunkDivideConfig%allowPostOverlaps, advance='yes')
    endif
    if ( BeVerbose( 'chu', -1 ) ) call dump ( chunks )
     L2Options%currentChunkNumber = chunkno
     OriginalOptions%currentChunkNumber = chunkno

  end subroutine GetChunkInfoFromMaster

  ! --------------------------------------------- L2MasterTask ----------
  subroutine L2MasterTask ( chunks )
    ! This is a `master' task for the l2 software
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS

    ! Local parameter
    integer, parameter :: MAXDIRECTWRITEFILES=200 ! For internal array sizing
    integer, parameter :: DATABASEINFLATION=100

    ! External (C) function

    ! Local variables
    logical :: MACHINEREQUESTQUEUED     ! Set if waiting for a free machine
    logical :: SKIPDELAY                ! Don't wait before doing the next go round
    logical :: SKIPDEATHWATCH           ! Don't check for deaths
    logical :: USINGOLDSUBMIT           ! Set if using the old submit mechanism
    logical :: USINGSUBMIT              ! Set if using the submit or l2q
    logical :: USINGL2Q                 ! Set if using the l2q queue manager
    character(len=MachineNameLen) :: THISNAME
    character(len=8) :: CHUNKNOSTR
    character(len=2048) :: COMMANDLINE
    character(len=16) :: DATESTRING

    integer :: BUFFERIDRCV              ! From PVM
    integer :: BUFFERIDSND              ! From PVM
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    integer :: CHUNK                    ! Loop counter
    integer :: CREATEFILE               ! Flag for direct writes
    integer :: DEADCHUNK                ! A chunk from a dead task
    integer :: DEADMACHINE              ! A machine for a dead task
    integer :: DEADTID                  ! A task that's died
    integer :: DUMMY                    ! From inflate database
    logical :: DUMPFULLDWREQS
    integer :: FILEINDEX                ! Index for a direct write
    integer :: INDX                     ! Where "T" occurs in utc string
    integer :: INFO                     ! From PVM
    integer :: L2QTID                   ! TID of queue manager
    integer :: MACHINE                  ! Index
    integer :: Me = -1                  ! String index for trace
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: NEXTCHUNK                ! A chunk number
    integer :: NEXTTICKET               ! For direct write handling
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: NODE                     ! A tree node
    integer :: NODIRECTWRITEFILES       ! Need to keep track of filenames
    integer :: NODIRECTWRITEREQUESTS    ! Number of (relevantish) directWrite requests
    integer :: NOMACHINES               ! Number of slaves
    integer :: REQUESTEDFILE            ! String index from slave
    integer :: REQUESTINDEX             ! Index of direct write request
    integer :: REQUESTINDEXA(1)         ! Result of minloc
    integer :: RETURNEDTICKET           ! Ticket for completed direct write
    character ( len=256 ) :: sbmtdSlaveArg1, sbmtdSlaveArg2
    character ( len=2048 ) :: sbmtdSlaveArguments ! To remove --chunks if present
    integer :: SIGNAL                   ! From slave
    integer :: SLAVETID                 ! One slave
    ! integer :: STAGEFILEID              ! From HDF5
    integer :: TIDARR(1)                ! One tid

    integer, dimension(size(chunks)) :: CHUNKMACHINES ! Machine indices for chunks
    integer, dimension(size(chunks)) :: CHUNKTIDS ! Tids for chunks
    character(len=16), dimension(size(chunks)) :: CHUNKNICETIDS ! Tids for chunks

    integer, dimension(maxDirectWriteFiles) :: DIRECTWRITEFILENAMES
    logical, dimension(maxDirectWriteFiles) :: DIRECTWRITEFILEBUSY
    real(r8), dimension(maxDirectWriteFiles) :: TIMEDWFILEBEGAN
    real(r8), dimension(maxDirectWriteFiles) :: MAXTIMEDWFILETOOK
    integer, dimension(maxDirectWriteFiles) :: NODIRECTWRITECHUNKS

    logical, dimension(size(chunks)) :: CHUNKSCOMPLETED ! Chunks completed
    logical, dimension(size(chunks)) :: CHUNKSSTARTED ! Chunks being processed
    logical, dimension(size(chunks)) :: CHUNKSABANDONED ! Chunks kept failing
    integer, dimension(size(chunks)) :: CHUNKFAILURES ! Failure count
    logical, dimension(size(chunks)) :: CHUNKSWRITING ! Which chunks are writing
    real(r8), dimension(maxDirectWriteFiles) :: MAXTIMECHUNKSPENTWRITING

    logical, save :: FINISHED = .false. ! This will be called multiple times

    type (Machine_T),dimension(:), pointer :: Machines
    type (Machine_T)                       :: thisMachine
    ! Local vector database
    type (StoredResult_T), dimension(:), pointer :: storedResults
    ! Map into the above arrays
    type (DirectWriteRequest_T), dimension(:), pointer :: DIRECTWRITEREQUESTS
    type (DirectWriteRequest_T), pointer :: REQUEST

    ! Executable code --------------------------------
    KILLINGSLAVESDELAY   = 100*parallel%delay
    ! First, is this the first call. The first call does all the work
    ! so if it's not then quit.
    if ( finished ) return

    call trace_begin ( me, "L2MasterTask", cond=toggle(gen))
    usingSubmit = trim(parallel%submit) /= ''
    usingL2Q = ( index(parallel%submit, 'l2q') > 0 )
    USINGOLDSUBMIT = USINGSUBMIT .and. .not. usingL2Q

    ! Setup some stuff
    noChunks = size(chunks)
    nullify ( machines, storedResults, &
      & directWriteRequests )
    noDirectWriteRequests = 0
    directWriteFileNames = 0
    noDirectWriteChunks = 0
    directWriteFileBusy = .false.
    timeDWFileBegan = 0.d0
    maxTimeDWFileTook = 0.d0
    maxTimeChunkSpentWriting = 0.d0
    noDirectWriteFiles = 0
    nextTicket = 1

    ! Setup the directWrite request database with a default size.
    dummy = InflateDirectWriteRequestDB ( directWriteRequests, DatabaseInflation )

    ! Work out the information on our virtual machine
    if ( usingL2Q ) then
      call RegisterWithL2Q(noChunks, machines, L2QTID)
      noMachines = size(machines)
      if ( noMachines < 1 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'No machines available for master to assign to slave tasks' )
      if ( .not. any(machines%OK) ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'No machines OK for master to assign to slave tasks' )
      if ( BeVerbose( 'mach', -1 ) .or. &
        & parallel%verbosity > 1 ) call dump ( machines )
    elseif ( .not. usingSubmit ) then
      call GetMachines ( machines )
      noMachines = size(machines)
      if ( noMachines < 1 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'No machines available for master to assign to slave tasks' )
      if ( .not. any(machines%OK) ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'No machines OK for master to assign to slave tasks' )
      if ( BeVerbose( 'mach', -1 ).or. &
        & parallel%verbosity > 1 ) call dump ( machines )
    else
      noMachines = 0
    end if

    ! Loop until all chunks are done
    chunksCompleted = .false.
    chunksStarted = .false.
    chunksAbandoned = .false.
    chunksWriting = .false.
    chunkFailures = 0
    chunkTids = 0
    CHUNKNICETIDS = ' '
    chunkMachines = 0
    ! Special switches to control which chunks to process
    ! by pre-abandoning the others right off the bat
    ! Note that this leaves output files possibly as big
    ! as if you ran all the chunks anyway
    
    ! Why isn't '3' actually '2' in the following test?
    if ( size(chunks) < 3 ) then
    elseif ( parallel%chunkRange /= '' ) then
      chunksAbandoned = .true.
      call ExpandStringRange(trim(parallel%chunkRange), chunksAbandoned, &
        & sense=.false.)
      if ( all(chunksAbandoned) ) then
        call KillSlaves ( .true., &
          & 'No chunks to process-was your request within range?' )
      endif
    endif
    if ( BeVerbose( 'l2q', -1 ) ) then
      call outputnamedValue( 'size(chunks)', size(chunks) )
      call outputnamedValue( 'parallel%chunkRange', trim(parallel%chunkRange) )
      call dump( chunksAbandoned, 'chunksAbandoned', options='N' )
    endif

    machineRequestQueued = .false. ! Request one machine at a time from L2Q
    masterLoop: do ! --------------------------- Master loop -----------------------
      skipDelay = .false.               ! Assume we're going to delay
      skipDeathWatch = .false.          ! Assume we'll listen for slave deaths
      ! This loop is in two main parts.

      ! In the first part, we look to see if there are any chunks still to be
      ! started done, and any vacant machines to do them.
      ! --------------------------------------------------------- Start new jobs? --
      do while ( chunkAndMachineReady() ) ! (nextChunk, machine) )
        if ( nextChunk < 1 ) then
          ! Should have returned false
          call KillSlaves (.true., 'Illegal chunk number' )
        elseif ( USINGOLDSUBMIT ) then ! -- Using a batch system
          ! We must remove any --chunk chunkRange from slaveArguments
          info = index(lowerCase(slaveArguments), '--chunk')
          if ( info > 0 ) then
            if ( info == 1 ) then
              sbmtdSlaveArg1 = ' '
            else
              sbmtdSlaveArg1 = slaveArguments(1:info-1)
            endif
            call RemoveNumFromList( slaveArguments(info:), sbmtdSlaveArguments, 1, ' ' )
            call RemoveNumFromList( sbmtdSlaveArguments, sbmtdSlaveArg2, 1, ' ' )
            sbmtdSlaveArguments = trim(sbmtdSlaveArg1) // ' ' // sbmtdSlaveArg2
          else
            sbmtdSlaveArguments = slaveArguments
          endif
          write ( chunkNoStr, '(i0)' ) nextChunk
          commandLine = &
            & trim(parallel%submit) // ' ' // &
            & trim(parallel%executable) // ' ' // &
            & ' --chunk ' // trim(chunkNoStr) // ' ' // &
            & trim(sbmtdSlaveArguments)
          call shell_command ( trim(commandLine) )
          chunksStarted(nextChunk) = .true.
          skipDeathWatch = .true.
          ! We'll have to wait for it to come one line later
          if ( parallel%verbosity > 0 ) then
            call output ( 'Submitted chunk ' )
            call TimeStamp ( nextChunk, advance='yes' )
          end if
        elseif ( machine < 1 ) then
          ! Just go on; should have returned .false. anyway
          exit
        else ! ------------------------------------- Start job using pvmspawn
          commandLine = trim(parallel%pgeName)   ! was 'mlsl2'
          if ( BeVerbose( 'slv', -1 ) ) then
            call PVMFCatchOut ( 1, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, "calling catchout" )
          end if
          info = myPVMSpawn ( trim(commandLine), PvmTaskHost, &
            & trim(machines(machine)%Name), &
            & 1, tidarr )
          call output (tidarr(1), advance = 'no')
          call output (' result is ', advance='no')
          call output (info, advance = 'yes')
          if ( parallel%verbosity > 0 ) then
            call output ( 'Tried to spawn ' )
            call TimeStamp ( trim(commandLine), advance='yes' )
            call output ( 'PvmTaskHost ' )
            call output ( PvmTaskHost, advance='yes' )
            call output ( 'on machine ' )
            call output ( trim(machines(machine)%Name), advance='yes' )
            call output ( 'result was ' )
            call output ( info, advance='yes' )
          end if

          ! Did this launch work
          if ( info == 1 ) then
            machines(machine)%free = .false.
            machines(machine)%tid = tidArr(1)
            machines(machine)%chunk = nextChunk
            indx = INDEX (L2PCF%startUTC, "T")
            if ( indx > 0 ) then
              DateString = L2PCF%startUTC(1:indx-1)
              call ReplaceSubString ( DateString, machines(machine)%master_Date, &
                & '-', 'd' )
            endif
            chunkMachines(nextChunk) = machine
            chunkTids(nextChunk) = tidArr(1)
            chunkNiceTids(nextChunk) = GetNiceTidString(chunkTids(nextChunk))
            chunksStarted(nextChunk) = .true.
            if ( parallel%verbosity > 0 ) then
              if ( BeVerbose( 'l2q', 0 ) ) then
                call output ( tidArr(1) )
                call output ( ' ' )
              endif
              call output ( 'Launched chunk ' )
              call output ( nextChunk )
              call TimeStamp ( ' on slave ' // trim(machines(machine)%name) // &
                & ' ' // trim(chunkNiceTids(nextChunk)), &
                & advance='yes' )
            end if
            if ( BeVerbose( 'l2q', 0 ) ) &
              & call output ( chunkTids(nextChunk), advance='yes' )
            call WelcomeSlave ( nextChunk, chunkTids(nextChunk) )
            if ( usingL2Q ) call ThankL2Q(machines(machine), L2Qtid)
            skipDeathWatch = .true.
          else
            ! Couldn't start this job, mark this machine as unreliable
            if ( parallel%verbosity > 0 ) then
              call output ( 'Unable to start slave task on ' // &
                & trim(machines(machine)%Name) // ' info=' )
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
            where ( machines%Name == machines(machine)%Name )
              machines%OK = .false.
            end where
            ! Send bad news back to l2 queue manager
            if ( usingL2Q ) then
              call TellL2QMachineDied( machines(machine), L2Qtid )
            endif
          end if
        end if
      end do

      ! --------------------------------------------------- Messages from jobs?
      ! In this next part, we listen out for communication from the slaves and
      ! process it accordingly.
      call PVMFNRecv( -1, InfoTag, bufferIDRcv )
      if ( bufferIDRcv < 0 ) then
        call PVMErrorMessage ( info, "checking for Info message" )
      else if ( bufferIDRcv > 0 ) then
        ! So we got a message.  There may be another one following on so don't delay
        ! before going round this loop again.
        skipDelay = .true.
        skipDeathWatch = .true.
        ! Who sent this?
        call PVMFBufInfo ( bufferIDRcv, bytes, msgTag, slaveTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking signal" )
        endif

        ! Who did this come from
        ! chunk = FindFirst ( chunkTids, slaveTid )
        machine = 0
        if ( associated(machines) ) &
          & machine = FindFirst ( machines%tid, slaveTid )
        if ( USINGOLDSUBMIT ) then
          chunk = FindFirst ( chunkTids, slaveTid )
          if ( chunk == 0 .and. signal /= sig_register ) then
            call output ( 'Signal is:' )
            call output ( signal )
            call outputNamedValue ( ' meaning', SigToName(signal) )
            call TimeStamp ( ' Tid: ' // trim ( GetNiceTidString ( slaveTid ) ), &
              & advance='yes' )
            call MLSL2Message ( MLSMSG_Warning, ModuleName, &
              & "Got a message from an unknown slave" )
          endif
        ! if ( chunk == 0 .and. &
        elseif ( machine == 0 ) then
          call output ( 'Signal is:' )
          call output ( signal )
          call outputNamedValue ( ' meaning', SigToName(signal) )
          call TimeStamp ( ' Tid: ' // trim ( GetNiceTidString ( slaveTid ) ), &
            & advance='yes' )
          call MLSL2Message ( MLSMSG_Warning, ModuleName, &
            & "Got a message from an unknown slave" )
          call dump ( chunkNiceTids, 'chunkNiceTids', options='t' )
          call TimeStamp ( slaveTid, advance='yes' )
          call dump ( chunkTids, 'chunkTids', format='places=10' )
          call dump ( machines%tid, 'machines%Tid', format='places=10' )
          cycle masterLoop
        else
          ! Unpack the first integer in the buffer
          if ( .not. USINGOLDSUBMIT .and. signal /= sig_register ) &
            & chunk = machines(machine)%chunk
            ! & machine = chunkMachines(chunk)
        end if

        select case (signal) 

        case ( sig_abouttodie ) ! --------------- Slave is about to die ---------
          call PVMF90Unpack ( commandLine, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking last gasp message" )
          endif
          write ( chunkNoStr, '(i0)' ) chunk
          commandLine = '(' // trim(chunkNoStr) // ')' // commandLine
          parallel%failedMsgs = catLists(parallel%failedMsgs, commandLine, '\')

        case ( sig_register ) ! ----------------- Chunk registering ------
          call PVMF90Unpack ( chunk, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking chunk number" )
          endif
          ! Note, we'll ignore the slave MAFNumber sent for fwmParallel stuff
          if ( USINGOLDSUBMIT ) then
            ! We only really care about this message if we're using submit
            
            ! A crude hack to fix a sometimes mlssubmit/pvm error
            ! in which multiple jobs are submitted for the same chunk
            ! (Why??)
            ! No time to explore, we're in the middle of v3 development
            !
            ! What to do?
            
            ! Kill the excess slave and continue
            ! (You're a cruel one, paw)
            if ( chunkTids(chunk) > 0 ) then
              call pvmfkill ( slaveTid, info )
              if ( info /= 0 ) &
                & call PVMErrorMessage ( info, 'killing slave' )
              call output ( 'Killed off excess slave ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' running chunk ' )
              call TimeStamp ( chunk, advance='yes' )
              cycle
            endif
            chunkTids(chunk) = slaveTid
            call GetMachineNameFromTid ( slaveTid, thisName, info )
            if ( info == -1 ) &
              & call KillSlaves (.true., 'Unable to get machine name from tid' )
            call WelcomeSlave ( chunk, slaveTid )
            if ( parallel%verbosity > 0 ) then
              call output ( 'Welcomed task ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' running chunk ' )
              call TimeStamp ( chunk, advance='yes' )
            end if
          endif

        case ( sig_tojoin ) ! --------------- Got a join request ---------
          call MLSL2Message ( MLSMSG_Error, ModuleName, &
            & 'This version does not support slave Join commands' )

        case ( sig_RequestDirectWrite ) ! ------- Direct write permission --
          ! What file did they ask for?
          call PVMUnpackStringIndex ( requestedFile, info )
          if ( info /= 0 )  call PVMErrorMessage ( info, &
            & "unpacking direct write request filename" )
          call PVMF90Unpack ( node, info )
          if ( info /= 0 )  call PVMErrorMessage ( info, &
            & "unpacking direct write request node" )
          ! Is this a new file?
          fileIndex = FindFirst ( directWriteFilenames(1:noDirectWriteFiles), &
            & requestedFile )
          if ( fileIndex == 0 ) then
            ! Clearly, if we don't know about this file it's new
            noDirectWriteFiles = noDirectWriteFiles + 1
            if ( noDirectWriteFiles > maxDirectWriteFiles ) &
              & call KillSlaves (.true., 'Too many direct write files, ' &
              & // 'increase limit maxDirectWriteFiles in ' // ModuleName )
            fileIndex = noDirectWriteFiles
            directWriteFilenames ( fileIndex ) = requestedFile
          end if
          if ( size(directWriteRequests) == noDirectWriteRequests ) then
            dummy = InflateDirectWriteRequestDB ( directWriteRequests, DatabaseInflation )
          end if
          noDirectWriteRequests = noDirectWriteRequests + 1
          request => directWriteRequests ( noDirectWriteRequests )
          ! Record this request, and have it 'take a ticket'
          request%chunk = chunk
          if ( .not. USINGOLDSUBMIT ) request%machine = machine
          request%node = node
          request%fileIndex = fileIndex
          request%ticket = nextTicket
          request%status = DW_Pending
          call time_now(request%whenMade)
          nextTicket = nextTicket + 1

          if ( parallel%verbosity > 1 ) then
            call output ( 'Direct write request from ' )
            if ( .not. USINGOLDSUBMIT ) &
              & call output ( trim(machines(machine)%Name) // ', ' )
            call output ( trim(GetNiceTidString(slaveTid)) )
            call output ( ' chunk ' )
            call output ( chunk )
            call output ( ' ticket ' )
            call output ( nextTicket - 1, advance='no' )
            call output ( ' index ' )
            call output ( fileIndex, advance='no' )
            call output ( ' file ' )
            call TimeStamp ( directWriteFilenames ( fileIndex ), advance='yes' )
            call display_string ( requestedFile, strip=.true., advance='yes' )
          end if

        case ( sig_DirectWriteFinished ) ! - Finished with direct write -
          ! Unpack the ticket number we got back
          call PVMF90Unpack ( returnedTicket, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, &
            & "unpacking returned ticket" )
          ! Record that the chunk has finished direct write
          requestIndex = FindFirst ( directWriteRequests%ticket, returnedTicket )

          request => directWriteRequests(requestIndex)
          request%status = DW_Completed
          call time_now(request%whenFinished)
          fileIndex = request%fileIndex
          directWriteFileBusy ( fileIndex ) = .false.
          chunksWriting ( chunk ) = .false.
          noDirectWriteChunks ( fileIndex ) = noDirectWriteChunks ( fileIndex ) + 1
          if ( parallel%verbosity > 1 ) then
            call output ( 'Direct write finished from ' )
            if ( .not. USINGOLDSUBMIT ) &
              & call output ( trim(machines(machine)%Name) // ', ' )
            call output ( trim(GetNiceTidString(slaveTid)) )
            call output ( ' chunk ' )
            call output ( chunk )
            call output ( ' ticket ' )
            call TimeStamp ( returnedTicket, advance='yes')
            if ( BeVerbose( 'dwtime', -1 ) ) then
              call output ( ' after ' )
              call output( timeDWHasBeenWriting( fileIndex, chunk ) )
              call output ( ' (s) ' )
            endif
            call display_string ( directWriteFilenames(fileIndex), &
              & strip=.true., advance='yes' )
          end if

        case ( sig_finished ) ! -------------- Got a finish message ----
          if ( parallel%verbosity > 0 ) then
            call output ( 'Got a finished message from ' )
            if ( .not. USINGOLDSUBMIT ) &
              & call output ( trim(machines(machine)%Name) // ' ' )
            call output ( trim(GetNiceTidString(slaveTid)) // &
              & ' processing chunk ' )
            call TimeStamp ( chunk, advance='yes')
          endif

          ! Send an acknowledgement
          call PVMFInitSend ( PVMDataDefault, bufferIDSnd )
          if ( bufferIdSnd < 0 ) &
            & call PVMErrorMessage ( bufferIDSnd, 'setting up finish ack.' )
          call PVMF90Pack ( SIG_AckFinish, info )
          if ( info /= 0 ) &
            & call PVMErrorMessage ( info, 'packing finish ack.' )
          call PVMFSend ( slaveTid, InfoTag, info )
          if ( info /= 0 ) &
            & call PVMErrorMessage ( info, 'sending finish ack.' )
          if ( parallel%verbosity > 1 ) &
            & call TimeStamp ( 'Acknowledgment sent', advance='yes')

          ! Now update our information
          chunksCompleted(chunk) = .true.
          chunkTids(chunk) = 0
          if ( .not. USINGOLDSUBMIT ) then
            machines(machine)%free = .true.
            ! Must wait on updating the following--
            ! pvm may tell us later this tid has quit
            ! machines(machine)%tid = 0
            ! machines(machine)%chunk = 0
            chunkMachines(chunk) = 0
          end if
          parallel%numCompletedChunks = parallel%numCompletedChunks + 1
          if ( parallel%verbosity > 0 ) then
            call printMasterStatus
          end if
          ! Send news back to l2 queue manager
          if ( usingL2Q ) then
            call TellL2QMachineFinished( &
              & trim(machines(machine)%name), machines(machine)%tid, L2Qtid, &
              & FIXDELAYFORSLAVETOFINISH )
          endif

        case default
          call KillSlaves (.true., 'Unkown signal from slave')
        end select

        ! Free the receive buffer
        call PVMFFreeBuf ( bufferIDRcv, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'freeing receive buffer' )
      end if

      ! ----------------------------------------------------- Administrative messages?
      ! Listen out for any message telling us to quit now
      call PVMFNRecv ( -1, GiveUpTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        if ( parallel%verbosity > 0 ) then
          call TimeStamp ( 'Received an external message to give up, so finishing now', &
            & advance='yes' )
        end if
        ! We're going to be quitting anyway, so might as well cut the delays
        FIXDELAYFORSLAVETOFINISH = parallel%delay
        KILLINGSLAVESDELAY = parallel%delay
        exit masterLoop
      end if

      ! Listen out for any message telling us to dump current status
      call PVMFNRecv ( -1, masterDumpTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        if ( parallel%verbosity > 0 ) then
          call TimeStamp ( 'Received an external message to dump our current status', &
            & advance='yes' )
        end if
        call dumpMastersStatus( chunkmachines, chunknicetids, &
          & directwritefilenames, directwritefilebusy, &
          & chunkscompleted, chunksabandoned, chunkswriting, machines )
      end if

      ! Listen out for any message telling us to switch to a new l2q
      call PVMFNRecv ( -1, sig_switchallegiance, bufferIDRcv )
      if ( bufferIDRcv > 0 .and. usingL2Q ) then
        if ( parallel%verbosity > 0 ) then
          call TimeStamp ( 'Received a message to switch to a new l2q', &
            & advance='yes' )
        end if
        call pvmfgettid(GROUPNAME, 0, L2Qtid)
        if ( L2Qtid < 1 ) then
          call KillSlaves (.true., 'switch l2q queue manager not running--dead or not started yet?')
        else
          call TimeStamp ( 'new l2q tid: ' )
          call TimeStamp ( L2Qtid, advance='yes' )
        endif
        ! Send an acknowledgement, pack machineRequestQueued, swear allegiance
        call PVMFInitSend ( PVMDataDefault, bufferIDSnd )
        if ( bufferIdSnd < 0 ) &
          & call PVMErrorMessage ( bufferIDSnd, 'setting up swearallegiance' )
        call PVMF90Pack ( machineRequestQueued, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing machineRequestQueued' )
        call PVMFSend ( L2Qtid, SIG_swearAllegiance, info )
        call TimeStamp ( 'Sending oath of allegiance; do we need host?: ' )
        call TimeStamp ( machineRequestQueued, advance='yes' )
      end if

      ! Listen out for any message telling us that a machine is OK again
      call PVMFNRecv ( -1, MachineFixedTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then 
        ! So we got a message.  There may be another one following on so don't delay
        ! before going round this loop again.
        skipDelay = .true.
        call PVMIDLUnpack ( thisName, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'unpacking machine fixed message' )
        if ( USINGOLDSUBMIT ) then
          call MLSL2Message ( MLSMSG_Warning, ModuleName, &
            & 'Got unexpected MachineFixed message but using submit method' )
        else
          if ( parallel%verbosity > 0 ) then
            call TimeStamp ( 'Received an external message to trust ' // &
              & trim(thisName) // ' again.' , advance='yes' )
          end if
          where ( machines%Name == thisName )
            machines%OK = .true.
            machines%jobsKilled = 0
          end where
        end if
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
        deadTid = 0
        deadChunk = 0
        deadMachine = 0
        ! Get the TID for the dead task
        call PVMF90Unpack ( deadTid, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking deadTid' )
        ! Now this may well be a legitimate exit, detectable by one of 2 cases
        ! either
        ! (1) we won't know about this tid any more
        ! (2) the machine status was reset to free after a finished signal
        ! Otherwise we need to tidy up.
        ! if ( any ( chunkTids == deadTid ) ) then
        !   deadChunk = FindFirst ( chunkTids, deadTid )
        if ( .not. USINGOLDSUBMIT ) then
          deadMachine = FindFirst ( machines%tid, deadTid )
          if ( deadMachine > 0 ) then
            ! On the other hand
            if ( .not. USINGOLDSUBMIT .and. machines(deadMachine)%free ) &
              & deadMachine = 0
          endif
        else
          deadChunk = FindFirst ( chunkTids, deadTid )
          deadMachine = deadChunk ! A trick--only deadChunk matters
        endif

        if ( deadMachine > 0 ) then
          ! Now, to get round a memory management bug, we'll ignore this
          ! if, as far as we're concerned, the task was finished anyway.
          ! if ( deadMachine /= 0 ) then
            if ( .not. USINGOLDSUBMIT ) then
              ! deadMachine = chunkMachines(deadChunk)
              deadChunk = machines(deadMachine)%chunk
              machines(deadMachine)%free = .true.
            end if
            if ( parallel%verbosity > 0 ) then
              if ( BeVerbose( 'l2q', -1 ) ) then
                call output ( deadTID )
                call output ( ' ' )
              endif
              call output ( 'The run of chunk ' )
              call output ( deadChunk )
              call output ( ' ' )
              if ( .not. USINGOLDSUBMIT ) &
                & call output ( 'on ' // trim(machines(deadMachine)%Name) // ' ' )
              call TimeStamp ( trim(GetNiceTidString(deadTid)) // &
                & ' died, try again.', advance='yes' )
            end if
            call CleanUpDeadChunksOutput ( deadChunk, storedResults )
            chunksStarted(deadChunk) = .false.
            chunkFailures(deadChunk) = chunkFailures(deadChunk) + 1
            if ( .not. USINGOLDSUBMIT ) then
              where ( machines(deadMachine)%Name == machines%Name )
                machines%jobsKilled = machines%jobsKilled + 1
              end where
            end if

            ! If the chunk posesses any direct writes, free them up.
            chunksWriting ( deadChunk ) = .false.
            ! First free any files that are in progress
            requestIndex = FindFirst ( directWriteRequests%chunk == deadChunk .and. &
              & directWriteRequests%status == DW_InProgress )
            if ( requestIndex /= 0 ) then
              request => directWriteRequests(requestIndex)
              directWriteFileBusy ( request%fileIndex ) = .false.
            if ( parallel%verbosity > 0 ) then
                call output ( 'Direct write died from ' )
                if ( .not. USINGOLDSUBMIT ) &
                  & call output ( trim(machines(deadMachine)%Name) // ', ' )
                call output ( trim(GetNiceTidString(deadTid)) )
                call output ( ' chunk ' )
                call output ( deadChunk )
                call output ( ' ticket ' )
                call TimeStamp ( request%ticket, advance='yes')
                if ( BeVerbose( 'dwtime', -1 ) ) then
                  call output ( ' after ' )
                  call output( timeDWHasBeenWriting( fileIndex, deadChunk ) )
                  call output ( ' (s) ' )
                endif
                call display_string ( directWriteFilenames(request%fileIndex), &
                  & strip=.true., advance='yes' )
              end if
            end if
            ! Now forget all the requests the dead chunk had pending
            where ( directWriteRequests%chunk == deadChunk )
              directWriteRequests%status = DW_Completed
            end where

            ! Does this chunk keep failing, if so, give up.
            if ( chunkFailures(deadChunk) >= &
              & parallel%maxFailuresPerChunk ) then
              chunksAbandoned(deadChunk) = .true.
            if ( parallel%verbosity > 0 ) then
                call output ( 'Chunk ' )
                call output ( deadChunk )
                call TimeStamp ( ' keeps dying.  Giving up on it.', &
                  & advance='yes' )
                call printMasterStatus
              end if
            end if

            ! Does this machine have a habit of killing jobs.  If so
            ! mark it as not OK.  We can't do much about it though if
            ! we're using submit.
            if ( .not. USINGOLDSUBMIT ) then
              if ( machines(deadMachine)%jobsKilled >= &
                & parallel%maxFailuresPerMachine ) then
              if ( parallel%verbosity > 0 ) &
                  & call TimeStamp ('The machine ' // &
                  & trim(machines(deadMachine)%Name) // &
                  & ' keeps killing things, marking it bad', &
                  & advance='yes' )
                where ( machines(deadMachine)%Name == machines%Name )
                  machines%OK = .false.
                end where
              end if
            end if
            
            ! Save dead chunk number, increment casualty figure
            parallel%failedChunks = catLists(parallel%failedChunks, deadChunk)
            if ( .not. USINGOLDSUBMIT ) &
            & parallel%failedMachs = &
            & catLists(parallel%failedMachs, trim(machines(deadMachine)%Name))
            parallel%numFailedChunks = parallel%numFailedChunks + 1
          else
            ! Otherwise we'd already forgotten about this slave, it told
            ! us it had finished.
            if ( parallel%verbosity > 0 .and. NOTFORGOTTEN ) call TimeStamp ( &
              & "A slave task died after giving results, " // &
              & "we won't worry about it.", &
              & advance='yes' )
          end if
        ! end if
        ! Send bad news back to l2 queue manager
        if ( usingL2Q .and. deadChunk /= 0 ) then
          ! call TellL2QMachineDied( trim(machines(machine)%name), L2Qtid )
          if ( deadTid /= machines(deadMachine)%tid ) then
            call output('deadChunk ', advance='no')
            call output(deadChunk , advance='yes')
            call output('deadTID ', advance='no')
            call output(deadTID , advance='yes')
            call output('machine%tid ', advance='no')
            call output(machines(deadMachine)%tid , advance='yes')
            call KillSlaves (.true., 'Dead slave tid doesnt match machine tid')
          endif
          call TellL2QMachineDied( machines(deadMachine), L2Qtid )
          if ( BeVerbose( 'l2q', -1 ) ) then
            call output ( 'Bad news about chunk ' )
            call output ( deadChunk )
            call TimeStamp ( ' on slave ' // trim(machines(deadMachine)%name) // &
              & ' ' // trim(chunkNiceTids(deadChunk)), &
              & advance='yes' )
          end if
          ! Must zero out its Tid 
          ! (so we won't try to free it when we're finished)
          chunkTids(deadChunk) = 0
        elseif ( usingL2Q ) then
          if ( BeVerbose( 'l2q', -1 ) ) then
            call output ( 'tID ' )
            call output ( deadTID )
            call TimeStamp ( ' finished normally ', advance='yes' )
          end if
        endif
        ! Update info about dead machine
        if ( deadMachine > 0 .and. .not. USINGOLDSUBMIT ) then
          machines(deadMachine)%ok = .false.
          machines(deadMachine)%tid = 0
          machines(deadMachine)%chunk = 0
        endif        
      else if ( bufferIDRcv < 0 ) then
        call PVMErrorMessage ( info, "checking for Notify message" )
      end if

      ! -------------------------------------------- Direct write permission logic ---

      ! Now hand out whatever direct write permissions we are able to give
      ! Give each request a value the same as their ticket number
      directWriteRequests%value = directWriteRequests%ticket
      ! Invalidate all the requests that are not pending
      where ( directWriteRequests%status /= DW_Pending )
        directWriteRequests%value = nextTicket + 1
      end where
      ! Invalidate all requests from a chunk that is busy
      do requestIndex = 1, noDirectWriteRequests
        if ( directWriteRequests(requestIndex)%status == DW_Pending ) then
          if ( chunksWriting(directWriteRequests(requestIndex)%chunk) ) &
            & directWriteRequests(requestIndex)%value = nextTicket + 1
        end if
      end do
      ! Now invalidate all corresponding to busy files
      do fileIndex = 1, noDirectWriteFiles
        if ( directWriteFileBusy ( fileIndex ) ) then
          where ( directWriteRequests%fileIndex == fileIndex )
            directWriteRequests%value = nextTicket + 1
          end where
        end if
      end do
      ! Now find the one with the best cost
      requestIndexA = minloc ( directWriteRequests%value )
      request => directWriteRequests ( requestIndexA(1) )
      ! If we can do this one then do so
      if ( request%value <= nextTicket ) then
        ! So we can do a direct write let's not delay the next time round the loop
        ! in case there are more things we can do immediately
        skipDelay = .true.
        ! OK, we can grant this request, record that in our information
        request%status = DW_InProgress
        call time_now(request%whenGranted)
        fileIndex = request%fileIndex
        directWriteFileBusy ( fileIndex ) = .true.
        call time_now ( timeDWFileBegan(fileIndex) )
        chunksWriting ( request%chunk ) = .true.
        if ( noDirectWriteChunks(fileIndex) == 0 ) then
          createFile = 1
        else
          createFile = 0
        end if
        if ( parallel%verbosity > 1 ) then
          call output ( 'Direct write granted to ' )
          if ( .not. USINGOLDSUBMIT ) &
            & call output ( trim(machines(request%machine)%Name) // ' ' )
          call output ( trim(GetNiceTidString(chunkTids(request%chunk))) )
          call output ( ' chunk ' )
          call output ( request%chunk )
          call output ( ' ticket ' )
          call output ( request%ticket, advance='no' )
          call output ( ' index ' )
          call TimeStamp ( request%fileIndex, advance='yes' )
          call output ( ' file ' )
          call output ( directWriteFilenames(request%fileIndex), advance='yes' )
          call display_string ( directWriteFilenames(request%fileIndex), strip=.true., &
            & advance='yes' )
        end if
        
        ! Tell the slave to go ahead
        call PVMFInitSend ( PvmDataDefault, bufferIdSnd )
        if ( bufferIdSnd < 0 ) &
          & call PVMErrorMessage ( bufferIDSnd, 'setting up direct write granted' )
        call PVMF90Pack ( SIG_DirectWriteGranted, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing direct write granted flag' )
        call PVMF90Pack ( (/ request%node, request%ticket, createFile, &
          & directWriteFilenames(request%fileIndex) /),&
          &  info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing direct write granted info' )
        
        call PVMPackStringIndex ( directWriteFilenames(request%fileIndex), info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "packing direct write request filename" )
        call PVMFSend ( chunkTids(request%chunk), InfoTag, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'sending direct write granted' )
      end if

      ! Perhaps compact the direct write request database
      if ( count ( directWriteRequests%status /= DW_Completed ) < &
        & noDirectWriteRequests / 2 .and. switchDetail(switches,'dwreq') == -1 ) then
        call CompactDirectWriteRequestDB ( directWriteRequests, noDirectWriteRequests )
        skipDelay = .true. ! This may take time so don't hang around later
      end if

      ! Have we abandoned everything?
      if ( all(chunksAbandoned) ) call KillSlaves ( .true., 'All chunks abandoned' )

      ! If we're done then exit
      if (all(chunksCompleted .or. chunksAbandoned)) then
        if ( parallel%verbosity > 0 ) &
          & call TimeStamp ( 'All chunks either processed or abandoned', advance='yes' )
        exit masterLoop
      end if

      ! Now deal with the case where we have no machines left to use.
      ! When all the chunks that have started have finished (surely this is a
      ! necessary condition anyway?) finish.
      ! Actually it is not a necessary condition as a machine could be listed more
      ! than one time (e.g. when it has multiple processors), and could
      ! have killed one too many jobs, but still have another job running.
      
      ! Oops--this was a bug--if using l2q, someone may revive one or more
      ! dead hosts so we can run, run again
      if ( .not. ( USINGOLDSUBMIT .or. USINGL2Q ) ) then
        if ( .not. any(machines%OK) .and. &
          &  all ( chunksStarted .eqv. chunksCompleted ) ) then
          if ( parallel%verbosity > 0 ) &
            & call TimeStamp ( 'No machines left to do the remaining work', &
            & advance='yes' )
          exit masterLoop
        end if
      end if

      ! Now, rather than chew up cpu time on the master machine, we'll wait a
      ! bit here.
      if ( .not. skipDelay .and. parallel%delay > 0 ) &
        & call usleep ( parallel%delay )
    end do masterLoop ! --------------------- End of master loop --------------

    if ( BeVerbose( 'l2q', -1 )) then
      call TimeStamp ( '   Finished', advance='yes' )
      call dump( chunkNiceTids, 'chunkNiceTids', options='t' )
      call dump( chunkTids, 'chunkTids' )
      call dump( machines%tid, 'machines%Tid', format='(i10)' )
    endif
    ! First kill any children still running (only if we got a give up message).
    call KillSlaves(.false., '')  ! use any string here because isKillMaster = false

    if ( usingL2Q ) then
      call TellL2QMasterFinished(L2Qtid)
      if ( parallel%verbosity > 0 ) &
        & call TimeStamp ( 'Telling l2q we are finished', advance='yes' )
    endif
    ! Now, we have to tidy some stuff up here to ensure we can join things
    where ( .not. chunksCompleted )
      chunksAbandoned = .true.
    end where

    if ( count(chunksCompleted) == 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'No chunks were processed successfully.' )

    if ( BeVerbose( 'dwtime', -1 ) ) then
      call dump( maxTimeDWFileTook, 'max time to do direct writes' )
      call dump( maxTimeChunkSpentWriting, 'max time chunks spent writing' )
    endif

    if ( parallel%verbosity > 0 ) then
      call TimeStamp ( 'All chunks processed, starting join task', advance='yes' )
    endif
    ! Now we join up all our results into l2gp and l2aux quantities

    call DestroyStoredResultsDatabase ( storedResults )

    finished = .true.
    if ( parallel%verbosity > 0 ) then
      call TimeStamp ( 'All chunks joined', advance='yes' )
    endif
    dumpfulldwreqs = switchDetail(switches,'dwreq1') > -1
    if ( switchDetail(switches,'dwreq') > -1 ) &
      & call dump(directWriteRequests, statsOnly=.not. dumpfulldwreqs)
    call trace_end ( "L2MasterTask", cond=toggle(gen) )

  contains

    logical function chunkAndMachineReady()  ! (nextChunk, machine, machineName)
      ! Return .true. if next chunk remaining to be done
      ! and a suitable machine can be found
      ! Also set nextChunk, machineName values (which could be made arguments)
      character(len=MACHINENAMELEN) :: machineName
      nextChunk = 0
      machine = 0
      machineName = ' '
      ! 1st check if any chunks remain
      chunkAndMachineReady = .not. all(chunksStarted .or. chunksAbandoned)
      if ( .not. chunkAndMachineReady ) return
      ! Now divide into 3 cases
      if ( .not. (USINGOLDSUBMIT .or. usingL2Q) ) then
      ! Case (1): Not using submit, nor l2q: master controls choice of machine
        chunkAndMachineReady = any(machines%Free .and. machines%OK)
        if ( .not. chunkAndMachineReady ) return
        nextChunk = FindFirst ( &
          & (.not. chunksStarted) .and. (.not. chunksAbandoned) )
        machine = FindFirst(machines%Free .and. machines%OK)
      elseif ( usingL2Q ) then
      ! Case (2): using l2q
        nextChunk = FindFirst ( &
          & (.not. chunksStarted) .and. (.not. chunksAbandoned) )
        chunkAndMachineReady = .not. machineRequestQueued
        if ( .not. chunkAndMachineReady ) then
          call ReceiveMachineFromL2Q(machineName)
        else
          call RequestMachineFromL2Q(L2Qtid, nextChunk, machineName)
        endif
        if ( len_trim(machineName) < 1 ) then
          ! Sorry--no machine ready yet; but our request is/remains queued
          chunkAndMachineReady = .false.
          machineRequestQueued = .true.
          ! if ( switchDetail(switches,'l2q') > -1 ) call output('Sorry, l2q scorns us', advance='yes')
        else
          ! Good news--a machine was/has become ready
          chunkAndMachineReady = .true.
          machineRequestQueued = .false.
          ! What if the machine has been added after this master began?
          ! Or, troubles other prevented machine from being added
          machine = FindFirst( (trim(machineName) == machines%Name) &
            & .and. machines%free )
          if ( machine < 1 ) then
            thisMachine%name = machineName
            machine = AddMachineToDataBase(machines, thisMachine)
            if ( BeVerbose( 'l2q', -1 ) ) &
              & call output('Added machine to db', advance='yes')
            ! machine = AddMachineNameToDataBase(machines%Name, machineName)
          endif
          ! if ( switchDetail(switches,'l2q') > -1 ) call output('Good news, l2q likes us', advance='yes')
          if ( BeVerbose( 'l2q', -1 ) ) then
            call outputNamedValue( 'machine', machine )
          endif
        endif
      else
      ! Case (3): Using submit
        nextChunk = FindFirst ( &
          & (.not. chunksStarted) .and. (.not. chunksAbandoned) )
      endif
    end function chunkAndMachineReady

    subroutine printMasterStatus
      ! Print status: completed, underway, abandoned
       call TimeStamp ( 'Master status:', advance='yes' )
       call output ( count(chunksCompleted) )
       call output ( ' of ' )
       call output ( noChunks )
       call output ( ' chunks completed, ')
       call output ( count(chunksStarted .and. .not. chunksCompleted) )
       call output ( ' underway, ' )
       call output ( count(chunksAbandoned) )
       call output ( ' abandoned, ' )
       call output ( count(.not. &
         & (chunksStarted .or. chunksCompleted .or. chunksAbandoned ) ) )
       call TimeStamp ( ' left. ', advance='yes' )
       if ( .not. ( USINGOLDSUBMIT .or. USINGL2Q ) ) then
         call output ( count ( .not. machines%Free ) )
         call output ( ' of ' )
         call output ( noMachines )
         call output ( ' machines busy, with ' )
         call output ( count ( .not. machines%OK ) )
         call TimeStamp ( ' being avoided.', advance='yes' )
       end if
    end subroutine printMasterStatus

    function timeDWHasBeenWriting ( fileIndex, chunk ) result ( howLong )
      ! How long has this chunk taken to write to this file?
      ! Args
      integer, intent(in) :: fileIndex
      integer, intent(in) :: chunk
      real(r8) :: howLong
      ! Internal variables
      real(r8) :: t2
      ! Executable
      call time_now ( t2 )
      howLong = t2 - timeDWFileBegan( fileIndex )
      maxTimeDWFileTook( fileIndex ) = max( howLong, maxTimeDWFileTook( fileIndex ) )
      maxTimeChunkSpentWriting( chunk ) = max( howLong, maxTimeChunkSpentWriting( chunk ) )
    end function timeDWHasBeenWriting

    subroutine WelcomeSlave ( chunk, tid )
      ! This routine welcomes a slave into the fold and tells it stuff
      integer, intent(in) :: CHUNK      ! Chunk we're asking it to do
      integer, intent(in) :: TID        ! Task ID for slave

      ! Local variables
      integer :: info

      call SendChunkInfoToSlave ( chunks, chunk, tid )
      ! Now ask to be notified when this task exits
      call PVMFNotify ( PVMTaskExit, NotifyTag, 1, (/ tid /), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, 'setting up notify' )
    end subroutine WelcomeSlave

    subroutine KillSlaves (isKillMaster, killMasterMsg)
      ! This routine kills all running slaves, and if isKillMaster is
      ! true, this routine will invoke MLSMessage with MLSError code
      logical, intent(in) :: isKillMaster
      character (len=*), intent(in) :: killMasterMsg
      integer :: Chunk ! Loop index

      call usleep( WaitBeforeKillingSlaves ) ! This is 130 s
      do chunk = 1, noChunks
        if ( chunkTids(chunk) /= 0 ) then
          call usleep ( KILLINGSLAVESDELAY )
          if ( usingL2Q ) then
            machine = chunkMachines(chunk)
            call TellL2QMachineFinished( &
              & trim(machines(machine)%name), machines(machine)%tid, L2Qtid, 0 )
            if ( BeVerbose( 'l2q', -1 ) ) then
              call output( 'tid: ', advance='no' )
              call output( machines(machine)%tid, advance='no' )
              call output( 'machine name: ', advance='no' )
              call output( trim(machines(machine)%name), advance='no' )
              call TimeStamp ( '   released', &
                & advance='yes' )
            endif
            call usleep ( parallel%delay )
          endif
          call output (chunkTids(chunk), advance='yes')
          call pvmfkill ( chunkTids(chunk), info )
          if ( info /= 0 ) &
            & call PVMErrorMessage ( info, 'killing slave' )
        end if
      end do

      if (isKillMaster) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, killMasterMsg )

    end subroutine KillSlaves

  end subroutine L2MasterTask

  ! ---------------------------------------- Private Procedures ----------
  ! -------------------------------------------- AddStoredQuantityToDatabase ----
  integer function AddStoredResultToDatabase ( database, item )
    ! Add a storedQuantity to a database, or create the database if it
    ! doesn't yet exist

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (StoredResult_T), dimension(:), pointer :: database
    type (StoredResult_T), intent(in) :: item

    ! Local variables
    type (StoredResult_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddStoredResultToDatabase = newSize
  end function AddStoredResultToDatabase

  ! --------------------------------- CleanUpDeadChunksOutput -------------
  subroutine CleanUpDeadChunksOutput ( chunk, storedResults )
    ! This destroys the vector info etc. associated with any output
    ! from a dead chunk.
    integer, intent(in) :: CHUNK        ! Index of chunk
    type (StoredResult_T), dimension(:), pointer :: STOREDRESULTS

    ! Local variables
    integer :: PRECIND                  ! Index for precision (if any)
    integer :: VALIND                   ! Index for value
    integer :: RES                      ! Loop counter

    ! Executable code
    if ( .not. associated(storedResults) ) return
    
    do res = 1, size ( storedResults )
      valInd = storedResults(res)%valInds(chunk)
      if ( valInd /= 0 ) then
        storedResults(res)%valInds(chunk) = 0
      end if

      if ( storedResults(res)%gotPrecision ) then
        precInd = storedResults(res)%precInds(chunk)
        if ( precInd /= 0 ) then
          storedResults(res)%precInds(chunk) = 0
        end if
      end if
    end do
  end subroutine CleanUpDeadChunksOutput
  
  ! ----------------------------------------- DestroyStoredResultsDatabase ---
  subroutine DestroyStoredResultsDatabase ( database )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    type (StoredResult_T), dimension(:), pointer :: database

    ! Local variabels
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: i, s, status

    ! Executable code
    if ( .not. associated ( database ) ) return
    do i = 1, size(database)
      call deallocate_test ( database(i)%valInds, 'database(?)%valInds', ModuleName)
      if ( associated ( database(i)%precInds ) ) &
        & call deallocate_test ( database(i)%precInds,&
        & 'database(?)%precInds', ModuleName)
    end do
    s = size(database) * storage_size(database) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
    deallocate ( database, STAT=status )
    call test_deallocate ( status, ModuleName, 'stored quantities database', s, &
      & address=addr )

  end subroutine DestroyStoredResultsDatabase

  ! ----------------------------------------- DumpMastersStatus ---
  subroutine DumpMastersStatus ( CHUNKMACHINES, CHUNKNICETIDS, &
    & DIRECTWRITEFILENAMES, DIRECTWRITEFILEBUSY, &
    & CHUNKSCOMPLETED, CHUNKSABANDONED, CHUNKSWRITING, MACHINES )
    integer, dimension(:) :: CHUNKMACHINES ! Machine indices for chunks
    character(len=16), dimension(:) :: CHUNKNICETIDS ! Tids for chunks
    integer, dimension(:) :: DIRECTWRITEFILENAMES
    logical, dimension(:) :: DIRECTWRITEFILEBUSY
    logical, dimension(:) :: CHUNKSCOMPLETED ! Chunks completed
    logical, dimension(:) :: CHUNKSABANDONED ! Chunks kept failing
    logical, dimension(:) :: CHUNKSWRITING ! Which chunks are writing
    type (Machine_T),dimension(:), pointer :: Machines
    integer, dimension(size(chunkmachines)) :: which, whichNot

    ! Local variabels
    integer :: i, m, n

    ! Executable code
    call banner ( 'Master Status Dumped', (/ 1, 80 /), 'C' )
    call output( 'Chunks', advance='yes' )
    m = 0
    do i = 1, size(chunkmachines)
      if ( chunkmachines(i) == 0 ) cycle
      m = m + 1
      call output( chunkmachines(i), advance='no' )
      call blanks( 3 )
      call output( chunksCompleted(i), advance='no' )
      call blanks( 3 )
      call output( chunksAbandoned(i), advance='no' )
      call blanks( 3 )
      call output( chunksWriting(i), advance='no' )
      call blanks( 3 )
      call output( trim(chunkNiceTids(i)), advance='yes' )
    enddo
    call output( 'Chunks abandoned', advance='yes' )
    call FindAll( chunksAbandoned, which, n, whichNot )
    if ( n > 0 ) &
      & call output( which(1:n), advance='yes' )
    call output( 'Chunks not completed', advance='yes' )
    call FindAll( chunksCompleted, which, n, whichNot )
    if ( n < m ) &
      & call output( whichNot(1:m-n), advance='yes' )
    call output( 'Chunks still writing', advance='yes' )
    call FindAll( chunksWriting, which, n, whichNot )
    if ( n > 0 ) call output( which(1:n), advance='yes' )
    call output( 'Direct Writes', advance='yes' )
    n = 0
    do i = 1, size(directwritefilenames)
      if ( directWriteFilenames(i) > 0 ) then
        n = n + 1
        call display_string ( directWriteFilenames(i), &
          & strip=.true., advance='no' )
        call blanks( 3 )
        call output( directwritefilebusy(i), advance='yes' )
      endif
    enddo
    call FindAll( directwritefilebusy, which, n, whichNot )
    if ( n > 0 ) then
      call output( 'Files still busy writing', advance='yes' )
      call output( directwritefilebusy(1:n), onlyIf = .true., advance='yes' )
    endif
    call dump ( machines )
  end subroutine DumpMastersStatus

  ! -------------------------------------- RecieveMachineFromL2Q --------------
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
        if ( BeVerbose( 'l2q', -1 ) ) &
          & call TimeStamp('Weve been granted a host', advance='yes')
        ! Granted us a machine
        call PVMF90Unpack ( machineName, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking machine name" )
        endif
        if ( len_trim(machineName) < 1 ) &
          & call MLSL2Message ( MLSMSG_Error, ModuleName, &
          & 'granted blank host name' )
        return
      else
        if ( beat < max(NBEATS, 1)) call usleep ( parallel%delay )
      endif
    enddo
    if ( BeVerbose( 'l2q', -1 ) .and. len_trim(machineName) > 0 ) &
      & call TimeStamp('Received from l2q ' // trim(machineName), advance='yes')
  end subroutine ReceiveMachineFromL2Q

  ! -------------------------------------- RegisterWithL2Q --------------
  subroutine RegisterWithL2Q(noChunks, machines, L2Qtid)
    integer, intent(in)                                  :: noChunks
    type(machine_t), dimension(:), pointer :: MACHINES
    ! character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    integer, intent(out)                                 :: L2Qtid
    !
    integer :: BUFFERID                 ! From PVM
    character(len=16) :: DATESTRING
    integer :: INDX
    integer :: INFO                     ! From PVM
    character(len=16) :: L2QSTRING
    !
    call pvmfgettid(GROUPNAME, 0, L2Qtid)
    if ( L2Qtid < 1 ) then
      call MLSL2Message( MLSMSG_Error, ModuleName, &
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
    call PVMF90Pack ( noChunks, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing number of chunks' )
    indx = INDEX (L2PCF%startUTC, "T")
    if ( len_trim(l2PCF%RunID) > 0 ) then
      dateString = l2PCF%RunID
    elseif ( indx > 0 ) then
    ! Now send our data date as a string
    ! We might reasonably expect one of two formats for startUTC:
    ! 1996-051T00:00:00.000000Z (using yyyy-dayofyear)
    ! 2006-05-07T00:00:00.000000Z (using yyyy-month-dayofmonth)
    ! (e.g., '2006d121')
      if ( all( NAppearances(L2PCF%startUTC(1:indx-1), (/'-'/)) > 1 ) ) then
        dateString = L2PCF%startUTC(1:indx-1) ! it was yyyy-month-dayofmonth
      else
        call ReplaceSubString ( L2PCF%startUTC(1:indx-1), dateString, &
        & '-', 'd' )  ! it was yyyy-dayofyear
      endif
    else
      dateString = '(unknown)'
    endif
    call PVMF90Pack ( dateString, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing date' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    ! call PVMF90Pack ( noChunks, info )
    ! if ( info /= 0 ) &
    !   & call PVMErrorMessage ( info, 'packing L2QTid' )
    ! call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
  end subroutine RegisterWithL2Q

  ! -------------------------------------- RequestMachineFromL2Q --------------
  subroutine RequestMachineFromL2Q(L2Qtid, nextChunk, machineName)
    !
    ! request a free host from the l2 queue manager
    integer, intent(in)                               :: L2Qtid
    integer, intent(in)                               :: nextChunk
    character(len=MACHINENAMELEN), intent(out)        :: MACHINENAME
    !
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Identify ourselves
    if ( BeVerbose( 'l2q', -1 ) ) &
      & call TimeStamp('Requesting host from l2q', advance='yes')
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( sig_requestHost, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing registration' )
    call PVMF90Pack ( nextChunk, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing number of chunk' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
    call ReceiveMachineFromL2Q(machineName)
  end subroutine requestMachineFromL2Q

  ! --------------------------------------------- SendChunkInfoToSlave ------
  subroutine SendChunkInfoToSlave ( chunks, chunkNo, slaveTid )
    ! This routine sends a chunk to a slave task

    ! Dummy arguments
    type (MLSChunk_T), intent(in), dimension(:) :: CHUNKS ! The chunk
    integer, intent(in) :: CHUNKNO      ! Which one is it.
    integer, intent(in) :: SLAVETID     ! Who to send it to

    ! Local variables
    integer :: BUFFERID                 ! ID for buffer to send
    integer :: INFO                     ! Flag from PVM
    integer :: CHUNK                    ! Loop counter
    integer :: BITS

    ! Executable code
    call BooleansToBits( (/ ChunkDivideConfig%allowPriorOverlaps, &
      & ChunkDivideConfig%allowPostOverlaps /), Bits )
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( (/ size(chunks), chunkNo /), info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing noChunks and chunkNo' )
    do chunk = 1, size(chunks)
      call PVMF90Pack ( &
        & (/ chunks(chunk)%firstMAFIndex, &
        &    chunks(chunk)%lastMAFIndex, &
        &    chunks(chunk)%noMAFsLowerOverlap, &
        &    chunks(chunk)%noMAFsUpperOverlap, &
        &    chunks(chunk)%chunkNumber /), info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing one chunk header' )

      if ( associated ( chunks(chunk)%hGridOffsets ) ) then
        call PVMF90Pack ( size ( chunks(chunk)%hGridOffsets ), info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing size(hGridOffsets)' )
        call PVMF90Pack ( chunks(chunk)%hGridOffsets, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing hGridOffsets chunk' )
      else
        call PVMF90Pack ( 0, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing 0 for size(hGridOffsets)' )
      end if

      if ( associated ( chunks(chunk)%hGridTotals ) ) then
        call PVMF90Pack ( size ( chunks(chunk)%hGridTotals ), info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing size(hGridTotals)' )
        call PVMF90Pack ( chunks(chunk)%hGridTotals, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing hGridTotals chunk' )
      else
        call PVMF90Pack ( 0, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing 0 for size(hGridTotals)' )
      end if
      
      ! The latest components added to the MLSChunk_T
      call PVMF90Pack ( chunks(chunk)%StartTime, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing Start time of chunk' )
      call PVMF90Pack ( chunks(chunk)%EndTime, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing End time of chunk' )
      

    end do
    call PVMF90Pack ( Bits, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing overlaps outside proc. rnge' )
    call PVMFSend ( slaveTid, ChunkTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending chunk' )
  end subroutine SendChunkInfoToSlave

  ! -------------------------------------- TellL2QMachineDied --------------
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
    if ( BeVerbose( 'l2q', -1 ) ) then
      call output('Telling l2q about death of host ', advance='no')
      call output(trim(machine%name), advance='no')
      call blanks(2)
      call TimeStamp(machine%tid, advance='yes')
    endif
  end subroutine TellL2QMachineDied

  ! -------------------------------------- TellL2QMachineFinished --------------
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
    if ( BeVerbose( 'l2q', -1 ) ) then
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

  ! -------------------------------------- TellL2QMasterFinished --------------
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

  ! -------------------------------------- ThankL2Q --------------
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
    call PVMF90Pack ( machine%chunk, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing machinechunk' )
    ! info = INDEX (L2PCF%startUTC, "T")
    ! if ( info > 0 ) machine%master_Date = L2PCF%startUTC(1:info-1)
    call PVMF90Pack ( machine%master_date, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'packing machineMasterDate' )
    call PVMFSend ( L2Qtid, petitionTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending finish packet' )
    if ( BeVerbose( 'l2q', -1 ) ) then
      call output('Thanking l2q for host ', advance='no')
      call TimeStamp(machine%tid, advance='yes')
    endif
  end subroutine ThankL2Q

  ! -------------------------------------- Never used here --------------
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module L2Parallel

!
! $Log$
! Revision 2.121  2020/07/22 22:57:14  pwagner
! Master now passes chunk Start, End times to slave
!
! Revision 2.120  2019/07/22 22:30:29  pwagner
! chunksAbandoned now will display the indices where true and the indices where false if array is logical-valued
!
! Revision 2.119  2018/09/13 20:26:07  pwagner
! Now properly sets L2Options%currentChunkNumber
!
! Revision 2.118  2018/07/27 23:18:48  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.117  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.116  2018/02/09 01:27:31  pwagner
! Output the meaning along with the signal number
!
! Revision 2.115  2017/12/22 00:34:41  pwagner
! WaitBeforeKillingSlaves; must be bigger than pgekilldelay
!
! Revision 2.114  2017/12/07 01:01:23  vsnyder
! Don't use host-associated variable as a DO index
!
! Revision 2.113  2016/02/29 19:50:46  pwagner
! Usleep got from machine module instead of being an external
!
! Revision 2.112  2015/03/28 02:48:22  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.111  2014/09/05 01:06:08  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.110  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.109  2014/04/22 18:17:22  pwagner
! Uses MLSMessage from MLSL2Options module
!
! Revision 2.108  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.107  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.106  2013/12/05 01:45:32  pwagner
! Pass RunID to l2q if non-blank; started using BeVerbose
!
! Revision 2.105  2013/08/30 02:45:43  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.104  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.103  2013/04/05 23:31:46  pwagner
! Added more output controllable by l2q switch setting
!
! Revision 2.102  2013/02/14 19:02:27  pwagner
! Added way for l2q to tell mster to dump status; removed unused stuff
!
! Revision 2.101  2012/07/02 20:38:02  pwagner
! Fix 'no live hosts so must stop' error with l2q; remove outdated staging stuff
!
! Revision 2.100  2012/03/28 20:05:54  pwagner
! slave tasks lost ability to join quantities
!
! Revision 2.99  2011/08/03 21:57:23  pwagner
! Improved routine printing and dumps by master task
!
! Revision 2.98  2011/05/26 20:39:34  pwagner
! Widened fields for dumping machinesTid
!
! Revision 2.97  2011/05/09 18:21:38  pwagner
! Converted to using switchDetail
!
! Revision 2.96  2010/05/13 23:47:19  pwagner
! Dropped superseded options run[first][last]
!
! Revision 2.95  2010/04/22 23:34:47  pwagner
! Print more useful error mesg if no chunks to process
!
! Revision 2.94  2009/09/18 00:30:25  pwagner
! Dont print misleading lines about busy, avoided machines if using l2q
!
! Revision 2.93  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.92  2009/06/16 17:41:19  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.91  2009/03/18 23:08:27  honghanh
! Kill all running slaves before letting the master dies upon encountering errors.
! Added subroutine KillSlaves
!
! Revision 2.90  2009/01/12 19:22:04  pwagner
! Alphabetized procedures; shorten time to finish when told to give up
!
! Revision 2.89  2008/06/12 22:28:27  pwagner
! Added options to print direct write timings
!
! Revision 2.88  2008/04/03 00:13:39  pwagner
! Not the nicest way to fix mlssubmit/pvm launching multiple instances of same slave
!
! Revision 2.87  2008/01/23 21:26:01  pwagner
! filters --chunks from slaveArguments
!
! Revision 2.86  2007/12/07 01:51:08  pwagner
! Removed unused dummy variables, etc.
!
! Revision 2.85  2007/11/01 23:31:47  pwagner
! Print clearer msgs when unable to launch slave
!
! Revision 2.84  2007/09/06 23:34:06  pwagner
! Delays a bit to prevent l2q from killing finishing slave
!
! Revision 2.83  2007/02/27 00:01:46  pwagner
! Change log message so grep -i died gives number died, not 2x
!
! Revision 2.82  2007/02/14 17:29:51  pwagner
! Lengthened time allowed for slaves to complete final tasks
!
! Revision 2.81  2007/02/13 21:36:19  pwagner
! Lengthen delays before telling l2q we finished so slaves may complete timings summary
!
! Revision 2.80  2007/01/12 00:38:57  pwagner
! Switches allegiance to a replacement l2q
!
! Revision 2.79  2006/11/22 18:11:04  pwagner
! Help l2q track hosts freed by killed masters
!
! Revision 2.78  2006/10/11 00:30:15  pwagner
! Fixed bug in length of arg MACHINENAME
!
! Revision 2.77  2006/09/29 00:29:42  pwagner
! Fixes bug where masters try to free dead hosts when finished
!
! Revision 2.76  2006/08/05 02:12:27  vsnyder
! Add ForWhom argument to ConstructVectorTemplate
!
! Revision 2.75  2006/04/20 23:21:37  pwagner
! Pass more chunk info from master to slave
!
! Revision 2.74  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.73  2005/03/24 21:22:32  pwagner
! Tried to fix slowing of masters relying on l2q
!
! Revision 2.72  2005/03/15 23:55:48  pwagner
! Saves error messages from about-to-die chunks
!
! Revision 2.71  2005/02/03 19:09:15  pwagner
! Receives master_date, master_time data from masters for each host
!
! Revision 2.70  2004/12/27 23:04:21  pwagner
! Fixed bug affecting old submit method
!
! Revision 2.69  2004/12/14 21:54:23  pwagner
! Changes related to l2q
!
! Revision 2.68  2004/09/16 23:57:31  pwagner
! Now tracks machine names of failed chunks
!
! Revision 2.67  2004/09/16 00:18:03  pwagner
! Keeps record of completed, failed chunks
!
! Revision 2.66  2004/08/05 22:47:47  pwagner
! New --chunkRange option to run selected chunks in parallel mode
!
! Revision 2.65  2004/08/03 23:14:56  pwagner
! Unseemly switch to allow running first last, or first and last chunks only
!
! Revision 2.64  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.63  2004/05/19 19:16:11  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.62  2004/04/16 00:47:54  livesey
! Now kills slaves on finishing (probably only happens when told to quit).
!
! Revision 2.61  2004/01/22 00:56:35  pwagner
! Fixed many bugs in auto-distribution of DirectWrites
!
! Revision 2.60  2004/01/02 23:36:00  pwagner
! DirectWrites may choose files automatically from db
!
! Revision 2.59  2003/11/15 00:34:38  pwagner
! Now stops after exactly maxNumFailures
!
! Revision 2.58  2003/11/14 23:37:13  pwagner
! Lets user change masterLoop delay via commandline option
!
! Revision 2.57  2003/09/24 23:41:55  pwagner
! Moved some USE statements around to workaround IFC bugs
!
! Revision 2.56  2003/08/01 20:28:05  pwagner
! Spawns slave tasks with parallel%pgeName as command
!
! Revision 2.55  2003/07/11 21:50:02  livesey
! Bug fix in a diagnostic message
!
! Revision 2.54  2003/07/11 01:23:59  livesey
! Minor bug fix and more informative output
!
! Revision 2.53  2003/07/08 17:30:43  livesey
! Bug fix in ticket issuing
!
! Revision 2.52  2003/07/07 17:32:00  livesey
! New approach to directWrite
!
! Revision 2.51  2003/07/02 00:54:10  livesey
! Various tidy ups and bug fixes
!
! Revision 2.50  2003/07/01 16:30:52  livesey
! Changed error to warning.
!
! Revision 2.49  2003/06/20 19:38:25  pwagner
! Allows direct writing of output products
!
! Revision 2.48  2003/06/05 23:53:34  livesey
! Made the diagnostic output less verbose.
!
! Revision 2.47  2003/05/22 02:23:59  livesey
! More informative error message
!
! Revision 2.46  2003/05/13 04:48:06  livesey
! Can now stage to temporary hdf5 file instead of in memory.
!
! Revision 2.45  2003/05/12 02:06:48  livesey
! Changed to use the inflation of vectors etc for efficiency.
!
! Revision 2.44  2003/05/10 22:30:12  livesey
! Tidy up a message
!
! Revision 2.43  2003/02/20 20:32:56  livesey
! Added the 'take a ticket' approach to direct write conflicts.
!
! Revision 2.42  2003/01/17 21:54:12  livesey
! Added the machineFixed stuff.
!
! Revision 2.41  2002/11/08 21:24:00  livesey
! Minor tidy ups associated with the non-submit mode of doing things.
! Now manages dead macines in a more transparent way.
!
! Revision 2.40  2002/10/17 18:19:24  livesey
! Added the GiveUp signal capability.
!
! Revision 2.39  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.38  2002/10/05 00:43:34  livesey
! Started work on the FwmParallel stuff
!
! Revision 2.37  2002/08/20 04:31:39  livesey
! Obscure limitation fixed.  Used to hang if chunk died while doing direct
! write.
!
! Revision 2.36  2002/07/18 02:15:29  livesey
! Bug fix, uninitialised variable
!
! Revision 2.35  2002/05/29 22:43:22  livesey
! Embarssing bug fix, finished=.true. was wrong side of if statement.
!
! Revision 2.34  2002/05/29 21:55:11  livesey
! Bug fixes, diagnostic message emitted when not -Smas, and associated
! checks around some deallocates.
!
! Revision 2.33  2002/05/22 00:48:36  livesey
! Added direct write stuff
!
! Revision 2.32  2002/05/21 01:12:05  livesey
! Got rid of machine stuff not relevant in submit case
!
! Revision 2.31  2002/05/08 16:19:29  livesey
! Bug fix.
!
! Revision 2.30  2002/04/24 20:20:53  livesey
! Submit now working, also if pvm specified as host filename, takes
! hostname list from pvm configuration
!
! Revision 2.29  2002/04/24 16:53:38  livesey
! Reordered arrays, on the way to having submit working
!
! Revision 2.28  2002/04/08 20:49:17  pwagner
! Swath name optionally passed to JoinL2GPQuantities
!
! Revision 2.27  2002/03/14 00:15:43  livesey
! Minor bug fix
!
! Revision 2.26  2002/01/10 01:07:43  livesey
! Another bug fix, slave wasn't receiving chunks correctly.
!
! Revision 2.25  2002/01/09 23:16:34  livesey
! Whoops, left the SendChunkInfoToSlave call in the wrong place!
!
! Revision 2.24  2002/01/09 22:55:57  livesey
! Now sends slaves all the chunks as regular HGrids need them.
!
! Revision 2.23  2001/11/14 18:03:45  livesey
! Changed for new FindFirst behavior
!
! Revision 2.22  2001/11/05 23:21:31  livesey
! Can now specify subset of slaves using <filename>:a:b.
!
! Revision 2.21  2001/10/30 01:45:21  livesey
! Some modifications/fixes to parallel join
!
! Revision 2.20  2001/09/10 23:39:12  livesey
! Minor change, not sure what it did, sorry.
!
! Revision 2.19  2001/09/09 02:52:43  livesey
! Now gets find first from a different place
!
! Revision 2.18  2001/09/08 00:21:44  pwagner
! Revised to work for new column Abundance in lone swaths
!
! Revision 2.17  2001/09/05 20:34:56  pwagner
! Reverted to pre-columnAbundance state
!
! Revision 2.16  2001/08/02 23:59:22  pwagner
! Compiles, but incomplete teatment of column abundance(s)
!
! Revision 2.15  2001/08/02 00:18:55  pwagner
! Began adding column quantities; incomplete
!
! Revision 2.14  2001/06/22 05:19:46  livesey
! Moved FindFirst into MLSL2Common
!
! Revision 2.13  2001/06/19 22:58:07  pwagner
! Eliminated duplicate declaration of MLSCommon
!
