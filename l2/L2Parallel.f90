! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2Parallel
  ! This module contains low level routines and stuff for dealing with parallel
  ! invocations of the MLSL2 program.

  ! Level 2 programs can be masters or slaves, or neither, not both.  A task
  ! which is neither simply runs through the l2cf as normal.  A master task
  ! does the chunk divide and then launches slave tasks for each chunk, and
  ! awaits their results in the Join section.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use Dumper, only: DUMP
  use Join, only: JOINL2GPQUANTITIES, JOINL2AUXQUANTITIES
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMF90PACK, PVMFMYTID, &
    & PVMF90UNPACK, PVMERRORMESSAGE, PVMTASKHOST, PVMFSPAWN, &
    & MYPVMSPAWN, PVMFCATCHOUT, PVMFSEND, PVMFNOTIFY, PVMTASKEXIT, &
    & GETMACHINENAMEFROMTID, PVMFFREEBUF
  use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
  use QuantityPVM, only: PVMSENDQUANTITY, PVMRECEIVEQUANTITY
  use MLSCommon, only: R8, MLSCHUNK_T, FINDFIRST
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VECTORTEMPLATE_T, &
    & CONSTRUCTVECTORTEMPLATE, &
    & CREATEVECTOR, DESTROYVECTORINFO, DESTROYVECTORTEMPLATEINFO, &
    & INFLATEVECTORTEMPLATEDATABASE, INFLATEVECTORDATABASE
  use Machine, only: SHELL_COMMAND
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE, &
    & MLSMSG_Deallocate, MLSMSG_WARNING
  use L2GPData, only: L2GPDATA_T
  use L2AUXData, only: L2AUXDATA_T
  use L2ParInfo, only: L2PARALLELINFO_T, PARALLEL, INFOTAG, CHUNKTAG, GIVEUPTAG, &
    & SIG_TOJOIN, SIG_FINISHED, SIG_ACKFINISH, SIG_REGISTER, NOTIFYTAG, &
    & SIG_REQUESTDIRECTWRITE, SIG_DIRECTWRITEGRANTED, SIG_DIRECTWRITEFINISHED, &
    & GETNICETIDSTRING, SLAVEARGUMENTS, MACHINENAMELEN, GETMACHINENAMES, &
    & MACHINEFIXEDTAG, DIRECTWRITEREQUEST_T, DW_PENDING, DW_INPROGRESS, DW_COMPLETED, &
    & INFLATEDIRECTWRITEREQUESTDB, COMPACTDIRECTWRITEREQUESTDB, DUMP
  use QuantityTemplates, only: QUANTITYTEMPLATE_T, &
    & DESTROYQUANTITYTEMPLATECONTENTS, INFLATEQUANTITYTEMPLATEDATABASE, &
    & NULLIFYQUANTITYTEMPLATE, DESTROYQUANTITYTEMPLATEDATABASE
  use Toggles, only: Gen, Switches, Toggle
  use Output_m, only: Output
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_STRING
  use String_table, only: Display_String
  use Init_Tables_Module, only: S_L2GP, S_L2AUX
  use MoreTree, only: Get_Spec_ID
  use MorePVM, only: PVMUNPACKSTRINGINDEX, PVMPACKSTRINGINDEX
  use VectorHDF5, only: WRITEVECTORASHDF5, READVECTORFROMHDF5
  use HDF5, only: H5FCREATE_F, H5FCLOSE_F, H5FOPEN_F, H5F_ACC_RDONLY_F, H5F_ACC_TRUNC_F

  use VectorsModule, only: CHECKINTEGRITY

  implicit none
  private

  public :: L2MasterTask
  public :: GetChunkInfoFromMaster
  
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! Parameters
  integer, parameter :: HDFNAMELEN = 132 ! Max length of name of swath/sd

  integer :: counterStart
  parameter ( counterStart = 2 * (huge (0) / 4 ) )

  ! Private types
  type StoredResult_T
    integer :: key                      ! Tree node for this join
    logical :: gotPrecision             ! If set have precision as well as value
    integer, dimension(:), pointer :: valInds=>NULL() ! Array vec. dtbs. inds (noChunks)
    integer, dimension(:), pointer :: precInds=>NULL() ! Array vec. dtbs. inds (noChunks)
    character(len=HDFNameLen) :: hdfName ! Name of swath/sd
  end type StoredResult_T

contains ! ================================ Procedures ======================

  ! --------------------------------------------- SendChunkInfoToSlave ----------
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

    ! Executable code
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

    end do
    call PVMFSend ( slaveTid, ChunkTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'sending chunk' )
  end subroutine SendChunkInfoToSlave

  ! ----------------------------------------------- GetChunkInfoFromMaster ------
  subroutine GetChunkInfoFromMaster ( chunks, chunkNo )
    ! This function gets a chunk sent by SendChunkToSlave

    ! Dummy arguments and result
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS
    integer, intent(out) :: CHUNKNO

    ! Local parameters
    integer, parameter :: noChunkTerms = 5 ! Number of components in MLSChunk_T

    ! Local variables
    integer :: BUFFERID                 ! ID for buffer for receive
    integer :: CHUNK                    ! Loop counter
    integer :: INFO                     ! Flag from PVM
    integer :: STATUS                   ! From allocate
    integer :: NOCHUNKS                 ! Size.
    integer :: NOHGRIDS                 ! Size of hGrid information
    integer, dimension(2) :: HEADER     ! No chunks, chunkNo
    integer, dimension(noChunkTerms) :: VALUES ! Chunk as integer array

    ! Executable code
    if ( .not. parallel%slave ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
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
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to allocate chunks' )

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

    end do

    if ( index(switches,'chu') /=0 ) call dump ( chunks )

  end subroutine GetChunkInfoFromMaster

  ! --------------------------------------------- L2MasterTask ----------
  subroutine L2MasterTask ( chunks, l2gpDatabase, l2auxDatabase )
    ! This is a `master' task for the l2 software
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    type (L2AuxData_T), dimension(:), pointer :: L2AUXDATABASE

    ! Local parameter
    integer, parameter :: DELAY = 200000  ! For Usleep, no. microsecs
    integer, parameter :: MAXDIRECTWRITEFILES=200 ! For internal array sizing
    integer, parameter :: DATABASEINFLATION=100

    ! External (C) function
    external :: Usleep

    ! Local variables
    logical :: USINGSUBMIT              ! Set if using the submit mechanism
    logical :: SKIPDELAY                ! Don't wait before doing the next go round
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    character(len=MachineNameLen) :: THISNAME
    character(len=8) :: CHUNKNOSTR
    character(len=2048) :: COMMANDLINE

    integer :: BUFFERIDRCV              ! From PVM
    integer :: BUFFERIDSND              ! From PVM
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    integer :: CHUNK                    ! Loop counter
    integer :: COMPLETEDFILE            ! String index from slave
    integer :: CREATEFILE               ! Flag for direct writes
    integer :: DEADCHUNK                ! A chunk from a dead task
    integer :: DEADFILE                 ! ID of file writen by dead task
    integer :: DEADMACHINE              ! A machine for a dead task
    integer :: DEADTID                  ! A task that's died
    integer :: DUMMY                    ! From inflate database
    integer :: FILEINDEX                ! Index for a direct write
    integer :: HDFNAMEINDEX             ! String index
    integer :: INFO                     ! From PVM
    integer :: MACHINE                  ! Index
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: NEXTCHUNK                ! A chunk number
    integer :: NEXTTICKET               ! For direct write handling
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: NODE                     ! A tree node
    integer :: NODIRECTWRITEFILES       ! Need to keep track of filenames
    integer :: NODIRECTWRITEREQUESTS    ! Number of (relevantish) directWrite requests
    integer :: NOMACHINES               ! Number of slaves
    integer :: NOQUANTITIESACCUMULATED  ! Running counter / index
    integer :: PRECIND                  ! Array index
    integer :: REQUESTEDFILE            ! String index from slave
    integer :: REQUESTINDEX             ! Index of direct write request
    integer :: REQUESTINDEXA(1)         ! Result of minloc
    integer :: RESIND                   ! Loop counter
    integer :: RETURNEDTICKET           ! Ticket for completed direct write
    integer :: SIGNAL                   ! From slave
    integer :: SLAVETID                 ! One slave
    integer :: STAGEFILEID              ! From HDF5
    integer :: STATUS                   ! From deallocate etc.
    integer :: TIDARR(1)                ! One tid
    integer :: VALIND                   ! Array index

    integer, dimension(size(chunks)) :: CHUNKMACHINES ! Machine indices for chunks
    integer, dimension(size(chunks)) :: CHUNKTIDS ! Tids for chunks

    integer, dimension(maxDirectWriteFiles) :: DIRECTWRITEFILENAMES
    logical, dimension(maxDirectWriteFiles) :: DIRECTWRITEFILEBUSY
    integer, dimension(maxDirectWriteFiles) :: NODIRECTWRITECHUNKS

    logical, dimension(:), pointer :: MACHINEFREE ! Is this machine busy
    logical, dimension(:), pointer :: MACHINEOK ! Is this machine working?
    integer, dimension(:), pointer :: JOBSMACHINEKILLED ! Counts failures per machine.
    logical, dimension(size(chunks)) :: CHUNKSCOMPLETED ! Chunks completed
    logical, dimension(size(chunks)) :: CHUNKSSTARTED ! Chunks being processed
    logical, dimension(size(chunks)) :: CHUNKSABANDONED ! Chunks kept failing
    integer, dimension(size(chunks)) :: CHUNKFAILURES ! Failure count
    logical, dimension(size(chunks)) :: CHUNKSWRITING ! Which chunks are writing

    logical, save :: FINISHED = .false. ! This will be called multiple times
    logical :: INTEGRITY

    type (QuantityTemplate_T),dimension(:), pointer :: joinedQuantities 
    ! Local quantity template database
    type (VectorTemplate_T), dimension(:), pointer  :: joinedVectorTemplates
    ! Local vector template database
    type (Vector_T), dimension(:), pointer :: joinedVectors
    ! Local vector database
    type (StoredResult_T), dimension(:), pointer :: storedResults
    ! Map into the above arrays
    type (VectorValue_T), pointer :: QTY
    type (VectorValue_T), pointer :: PRECQTY
    type (Vector_T), target :: MYVALUES
    type (Vector_T), target :: MYPRECISIONS
    type (DirectWriteRequest_T), dimension(:), pointer :: DIRECTWRITEREQUESTS
    type (DirectWriteRequest_T), pointer :: REQUEST

    ! Executable code --------------------------------

    ! First, is this the first call. The first call does all the work
    ! so if it's not then quit.
    if ( finished ) return

    usingSubmit = trim(parallel%submit) /= ''

    ! Setup some stuff
    noChunks = size(chunks)
    nullify ( joinedQuantities, joinedVectorTemplates, joinedVectors, &
      & machineNames, machineFree, storedResults, machineOK, jobsMachineKilled, &
      & directWriteRequests )
    noDirectWriteRequests = 0
    directWriteFileNames = 0
    noDirectWriteChunks = 0
    directWriteFileBusy = .false.
    noDirectWriteFiles = 0
    nextTicket = 1
    noQuantitiesAccumulated = 0

    ! Setup the directWrite request database with a default size.
    dummy = InflateDirectWriteRequestDB ( directWriteRequests, DatabaseInflation )

    ! Work out the information on our virtual machine
    if ( .not. usingSubmit ) then
      call GetMachineNames ( machineNames )
      noMachines = size(machineNames)
    else
      noMachines = 0
    end if
    call Allocate_test ( machineFree, noMachines, 'machineFree', ModuleName )
    call Allocate_test ( machineOK, noMachines, 'machineOK', ModuleName )
    call Allocate_test ( jobsMachineKilled, noMachines, 'jobsMachineKilled', ModuleName )
    machineFree = .true.
    machineOK = .true.
    jobsMachineKilled = 0

    ! Setup the staging file if we're using one
    if ( .not. parallel%stageInMemory ) then
      call H5FCreate_F ( trim(parallel%stagingFile), H5F_ACC_TRUNC_F, stageFileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to open staging file: '//trim(parallel%stagingFile) )
    end if

    ! Loop until all chunks are done
    chunksCompleted = .false.
    chunksStarted = .false.
    chunksAbandoned = .false.
    chunksWriting = .false.
    chunkFailures = 0
    chunkTids = 0
    chunkMachines = 0

    masterLoop: do ! --------------------------- Master loop -----------------------
      skipDelay = .false.               ! Assume we're going to delay
      ! This loop is in two main parts.

      ! In the first part, we look to see if there are any chunks still to be
      ! started done, and any vacant machines to do them.
      ! --------------------------------------------------------- Start new jobs? --
      do while ( (.not. all(chunksStarted .or. chunksAbandoned)) .and. &
        & ( any(machineFree .and. machineOK) .or. usingSubmit ) )
        nextChunk = FindFirst ( (.not. chunksStarted) .and. (.not. chunksAbandoned) )
        if ( usingSubmit ) then ! --------------------- Using a batch system
          write ( chunkNoStr, '(i0)' ) nextChunk
          commandLine = &
            & trim(parallel%submit) // ' ' // &
            & trim(parallel%executable) // ' ' // &
            & ' --chunk ' // trim(chunkNoStr) // ' ' // &
            & trim(slaveArguments)
          call shell_command ( trim(commandLine) )
          chunksStarted(nextChunk) = .true.
          ! We'll have to wait for it to come one line later
          if ( index(switches,'mas') /= 0 ) then
            call output ( 'Submitted chunk ' )
            call output ( nextChunk, advance='yes' )
          end if
        else ! ----------------------------------------- Start job using pvmspawn
          machine = FindFirst(machineFree .and. machineOK)
          commandLine = 'mlsl2'
          if ( index(switches,'slv') /= 0 ) then
            call PVMFCatchOut ( 1, info )
            if ( info /= 0 ) call PVMErrorMessage ( info, "calling catchout" )
          end if
          info = myPVMSpawn ( trim(commandLine), PvmTaskHost, &
            & trim(machineNames(machine)), &
            & 1, tidarr )

          ! Did this launch work
          if ( info == 1 ) then
            machineFree(machine) = .false.
            chunkMachines(nextChunk) = machine
            chunkTids(nextChunk) = tidArr(1)
            chunksStarted(nextChunk) = .true.
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'Launched chunk ' )
              call output ( nextChunk )
              call output ( ' on slave ' // trim(machineNames(machine)) // &
                & ' ' // trim(GetNiceTidString(chunkTids(nextChunk))), &
                & advance='yes' )
            end if
            call WelcomeSlave ( nextChunk, chunkTids(nextChunk) )
          else
            ! Couldn't start this job, mark this machine as unreliable
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'Unable to start slave task on ' // &
                & trim(machineNames(machine)) // ' info=' )
              if ( info < 0 ) then
                call output ( info, advance='yes' )
              else
                call output ( tidArr(1), advance='yes' )
              end if
              call output ( 'Marking this machine as not usable', advance='yes')
            end if
            ! Mark all instances of this machine as not to be used.
            where ( machineNames == machineNames(machine) )
              machineOK = .false.
            end where
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
        ! Who sent this?
        call PVMFBufInfo ( bufferIDRcv, bytes, msgTag, slaveTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking signal" )
        endif

        ! Who did this come from
        chunk = FindFirst ( chunkTids == slaveTid )
        if ( chunk == 0 .and. &
          &  (.not. usingSubmit .or. signal /= sig_register) ) then
          call output ( 'Signal is:' )
          call output ( signal )
          call output ( ' Tid: ' // trim ( GetNiceTidString ( slaveTid ) ), &
            & advance='yes' )
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & "Got a message from an unknown slave")
          cycle masterLoop
        else
          ! Unpack the first integer in the buffer
          if ( .not. usingSubmit .and. signal /= sig_register ) &
            & machine = chunkMachines(chunk)
        end if

        select case (signal) 

        case ( sig_register ) ! ----------------- Chunk registering ------
          call PVMF90Unpack ( chunk, info )
          if ( info /= 0 ) then
            call PVMErrorMessage ( info, "unpacking chunk number" )
          endif
          ! Note, we'll ignore the slave MAFNumber sent for fwmParallel stuff
          if ( usingSubmit ) then
            ! We only really care about this message if we're using submit
            chunkTids(chunk) = slaveTid
            call GetMachineNameFromTid ( slaveTid, thisName )
            call WelcomeSlave ( chunk, slaveTid )
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'Welcomed task ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' running chunk ' )
              call output ( chunk, advance='yes' )
            end if
          endif

        case ( sig_tojoin ) ! --------------- Got a join request ---------
          call StoreSlaveQuantity( joinedQuantities, &
            & joinedVectorTemplates, joinedVectors, &
            & storedResults, chunk, noChunks, slaveTid, noQuantitiesAccumulated, &
            & stageFileID )

        case ( sig_RequestDirectWrite ) ! ------- Direct write permission --
          ! What file did they ask for?
          call PVMUnpackStringIndex ( requestedFile, info )
          if ( info /= 0 )  call PVMErrorMessage ( info, &
            & "unpacking direct write request filename" )
          call PVMF90Unpack ( node, info )
          if ( info /= 0 )  call PVMErrorMessage ( info, &
            & "unpacking direct write request node" )
          ! Is this a new file?
          fileIndex = FindFirst ( directWriteFilenames(1:noDirectWriteFiles) == &
            & requestedFile )
          if ( fileIndex == 0 ) then
            ! Clearly, if we don't know about this file it's new
            noDirectWriteFiles = noDirectWriteFiles + 1
            if ( noDirectWriteFiles > maxDirectWriteFiles ) &
              & call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Too many direct write files, increase limit maxDirectWriteFiles in ' &
              & // ModuleName )
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
          request%machine = machine
          request%node = node
          request%fileIndex = fileIndex
          request%ticket = nextTicket
          request%status = DW_Pending
          nextTicket = nextTicket + 1

          if ( index ( switches, 'mas' ) /= 0 ) then
            call output ( 'Direct write request from ' )
            if ( .not. usingSubmit ) &
              & call output ( trim(machineNames(machine)) // ', ' )
            call output ( trim(GetNiceTidString(slaveTid)) )
            call output ( ' chunk ' )
            call output ( chunk )
            call output ( ' ticket ' )
            call output ( nextTicket - 1, advance='yes' )
            call display_string ( requestedFile, strip=.true., advance='yes' )
          end if

        case ( sig_DirectWriteFinished ) ! - Finished with direct write -
          ! Unpack the ticket number we got back
          call PVMF90Unpack ( returnedTicket, info )
          if ( info /= 0 ) call PVMErrorMessage ( info, &
            & "unpacking returned ticket" )
          ! Record that the chunk has finished direct write
          requestIndex = FindFirst ( directWriteRequests%ticket == returnedTicket )

          request => directWriteRequests(requestIndex)
          request%status = DW_Completed
          fileIndex = request%fileIndex
          directWriteFileBusy ( fileIndex ) = .false.
          chunksWriting ( chunk ) = .false.
          noDirectWriteChunks ( fileIndex ) = noDirectWriteChunks ( fileIndex ) + 1
          if ( index ( switches, 'mas' ) /= 0 ) then
            call output ( 'Direct write finished from ' )
            if ( .not. usingSubmit ) &
              & call output ( trim(machineNames(machine)) // ', ' )
            call output ( trim(GetNiceTidString(slaveTid)) )
            call output ( ' chunk ' )
            call output ( chunk )
            call output ( ' ticket ' )
            call output ( returnedTicket, advance='yes')
            call display_string ( directWriteFilenames(fileIndex), &
              & strip=.true., advance='yes' )
          end if

        case ( sig_finished ) ! -------------- Got a finish message ----
          if ( index(switches,'mas') /= 0 ) then
            call output ( 'Got a finished message from ' )
            if ( .not. usingSubmit ) &
              & call output ( trim(machineNames(machine)) // ' ' )
            call output ( trim(GetNiceTidString(slaveTid)) // &
              & ' processing chunk ' )
            call output ( chunk, advance='yes')
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

          ! Now update our information
          chunksCompleted(chunk) = .true.
          chunkTids(chunk) = 0
          if ( .not. usingSubmit ) then
            machineFree(machine) = .true.
            chunkMachines(chunk) = 0
          end if
          if ( index(switches,'mas') /= 0 ) then
            call output ( 'Master status:', advance='yes' )
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
            call output ( ' left. ', advance='yes' )
            if ( .not. usingSubmit ) then
              call output ( count ( .not. machineFree ) )
              call output ( ' of ' )
              call output ( noMachines )
              call output ( ' machines busy, with ' )
              call output ( count ( .not. machineOK ) )
              call output ( ' being avoided.', advance='yes' )
            end if
          end if

        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unkown signal from slave' )
        end select

        ! Free the recieve buffer
        call PVMFFreeBuf ( bufferIDRcv, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'freeing receive buffer' )
      end if

      ! ----------------------------------------------------- Administrative messages?
      ! Listen out for any message telling us to quit now
      call PVMFNRecv ( -1, GiveUpTag, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        if ( index(switches,'mas') /= 0 ) then
          call output ( 'Received an external message to give up, so finishing now', &
            & advance='yes' )
        end if
        exit masterLoop
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
        if ( usingSubmit ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Got unexpected MachineFixed message but using submit method' )
        else
          if ( index(switches,'mas') /= 0 ) then
            call output ( 'Received an external message to trust ' // &
              & trim(thisName) // ' again.' , advance='yes' )
          end if
          where ( machineNames == thisName )
            machineOK = .true.
            jobsMachineKilled = 0
          end where
        end if
      end if

      ! Listen out for any message that a slave task has died
      call PVMFNRecv ( -1, NotifyTAG, bufferIDRcv )
      if ( bufferIDRcv > 0 ) then
        ! So we got a message.  There may be another one following on so don't delay
        ! before going round this loop again.
        skipDelay = .true.
        ! Get the TID for the dead task
        call PVMF90Unpack ( deadTid, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, 'unpacking deadTid' )
        ! Now this may well be a legitimate exit, in which case, we won't
        ! know about this tid any more.  Otherwise we need to tidy up.
        if ( any ( chunkTids == deadTid ) ) then
          deadChunk = FindFirst ( chunkTids == deadTid )

          ! Now, to get round a memory management bug, we'll ignore this
          ! if, as far as we're concerned, the task was finished anyway.
          if ( deadChunk /= 0 ) then
            if ( .not. usingSubmit ) then
              deadMachine = chunkMachines(deadChunk)
              machineFree(deadMachine) = .true.
            end if
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'The run of chunk ' )
              call output ( deadChunk )
              call output ( ' ' )
              if ( .not. usingSubmit ) &
                & call output ( 'on ' // trim(machineNames(deadMachine)) // ' ' )
              call output ( trim(GetNiceTidString(deadTid)) // &
                & ' died, try again.', advance='yes' )
            end if
            call CleanUpDeadChunksOutput ( deadChunk, joinedQuantities, &
              & joinedVectorTemplates, joinedVectors, storedResults )
            chunksStarted(deadChunk) = .false.
            chunkFailures(deadChunk) = chunkFailures(deadChunk) + 1
            if ( .not. usingSubmit ) then
              where ( machineNames(deadMachine) == machineNames )
                jobsMachineKilled = jobsMachineKilled + 1
              end where
            end if

            ! If the chunk posesses any direct writes, free them up.
            chunksWriting ( deadChunk ) = .false.
            ! First free any files that are in progress
            requestIndex = FindFirst ( directWriteRequests%chunk == deadChunk .and. &
              & directWriteRequests%status == DW_InProgress )
            if ( requestIndex /= 0 ) directWriteFileBusy ( &
              & directWriteRequests(requestIndex)%fileIndex ) = .false.
            ! Now forget all the requests we had pending
            where ( directWriteRequests%chunk == deadChunk )
              directWriteRequests%status = DW_Completed
            end where

            ! Does this chunk keep failing, if so, give up.
            if ( chunkFailures(deadChunk) > &
              & parallel%maxFailuresPerChunk ) then
              if ( index(switches,'mas') /= 0 ) then
                call output ( 'Chunk ' )
                call output ( deadChunk )
                call output ( ' keeps dying.  Giving up on it.', &
                  & advance='yes' )
              end if
              chunksAbandoned(deadChunk) = .true.
            end if

            ! Does this machine have a habit of killing jobs.  If so
            ! mark it as not OK.  We can't do much about it though if
            ! we're using submit.
            if ( .not. usingSubmit ) then
              if ( jobsMachineKilled(deadMachine) > &
                & parallel%maxFailuresPerMachine ) then
                if ( index(switches,'mas') /= 0 ) &
                  & call output ('The machine ' // &
                  & trim(machineNames(deadMachine)) // &
                  & ' keeps killing things, marking it bad', &
                  & advance='yes' )
                where ( machineNames(deadMachine) == machineNames )
                  machineOK = .false.
                end where
              end if
            end if
          else
            ! Otherwise we'd already forgotten about this slave, it told
            ! us it had finished.
            if ( index(switches,'mas') /= 0 ) call output ( &
              & "A slave task died after giving results, " // &
              & "we won't worry about it.", &
              & advance='yes' )
          end if
        end if

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
        fileIndex = request%fileIndex
        directWriteFileBusy ( fileIndex ) = .true.
        chunksWriting ( request%chunk ) = .true.
        if ( noDirectWriteChunks(fileIndex) == 0 ) then
          createFile = 1
        else
          createFile = 0
        end if
        if ( index(switches,'mas') /= 0 ) then
          call output ( 'Direct write granted to ' )
          if ( .not. usingSubmit ) &
            & call output ( trim(machineNames(request%machine)) // ' ' )
          call output ( trim(GetNiceTidString(chunkTids(request%chunk))) )
          call output ( ' chunk ' )
          call output ( request%chunk )
          call output ( ' ticket ' )
          call output ( request%ticket, advance='yes' )
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
        call PVMF90Pack ( (/ request%node, request%ticket, createFile /), info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'packing direct write granted info' )
        
        call PVMFSend ( chunkTids(request%chunk), InfoTag, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, 'sending direct write granted' )
      end if

      ! Perhaps compact the direct write request database
      if ( count ( directWriteRequests%status /= DW_Completed ) < &
        & noDirectWriteRequests / 2 ) then
        call CompactDirectWriteRequestDB ( directWriteRequests, noDirectWriteRequests )
        skipDelay = .true. ! This may take time so don't hang around later
      end if

      ! Have we abandoned everything?
      if ( all(chunksAbandoned) ) call MLSMessage( MLSMSG_Error, ModuleName, &
        & 'All chunks abandoned' )

      ! If we're done then exit
      if (all(chunksCompleted .or. chunksAbandoned)) then
        if ( index(switches,'mas') /= 0 ) &
          & call output ( 'All chunks either processed or abandoned', advance='yes' )
        exit masterLoop
      end if

      ! Now deal with the case where we have no machines left to use.
      ! When all the chunks that have started have finished (surely this is a
      ! necessary condition anyway?) finish.
      ! Actually it is not a necessary condition as a machine could be listed more
      ! than one time (e.g. when it has multiple processors), and could
      ! have killed one too many jobs, but still have another job running.
      if ( .not. usingSubmit ) then
        if ( .not. any(machineOK) .and. &
          &  all ( chunksStarted .eqv. chunksCompleted ) ) then
          if ( index(switches,'mas') /= 0 ) &
            & call output ( 'No machines left to do the remaining work', &
            & advance='yes' )
          exit masterLoop
        end if
      end if

      ! Now, rather than chew up cpu time on the master machine, we'll wait a
      ! bit here.
      if ( .not. skipDelay ) call usleep ( delay )
    end do masterLoop ! --------------------- End of master loop --------------

    ! Now, we have to tidy some stuff up here to ensure we can join things
    where ( .not. chunksCompleted )
      chunksAbandoned = .true.
    end where

    if ( count(chunksCompleted) == 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'No chunks were processed successfully.' )

    if ( index(switches,'mas') /= 0 ) then
      call output ( 'All chunks processed, starting join task', advance='yes' )
    endif
    ! Now we join up all our results into l2gp and l2aux quantities

    ! If we're staging to a file, close the file and reopen it.
    if ( .not. parallel%stageInMemory ) then
      call H5FClose_F ( stageFileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close staging file:'//parallel%stagingFile )
      call H5FOpen_F ( trim(parallel%stagingFile), H5F_ACC_RDONLY_F, stageFileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to reopen staging file:'//parallel%stagingFile )
    end if

    do resInd = 1, size ( storedResults )
      hdfNameIndex = enter_terminal ( trim(storedResults(resInd)%hdfName), t_string )
      if ( index(switches,'mas') /= 0 ) then
        call output ( 'Joining ' // trim ( storedResults(resInd)%hdfName ), &
          & advance='yes' )
      endif
      
      do chunk = 1, noChunks
        if (.not. chunksAbandoned(chunk) ) then
          ! Setup for this quantity
          if ( parallel%stageInMemory ) then
            valInd = storedResults(resInd)%valInds(chunk)
            qty => joinedVectors(valInd)%quantities(1)

            if ( storedResults(resInd)%gotPrecision ) then
              precInd = storedResults(resInd)%precInds(chunk)
              precQty => joinedVectors(precInd)%quantities(1)
            else
              nullify ( precQty )
            endif
          else
            write ( thisName, '(i0)' ) storedResults(resInd)%valInds(chunk)
            call ReadVectorFromHDF5 ( stageFileID, trim(thisName), &
              & myValues, joinedQuantities )
            qty => myValues%quantities(1)
            if ( storedResults(resInd)%gotPrecision ) then
              write ( thisName, '(i0)' ) storedResults(resInd)%precInds(chunk)
              call ReadVectorFromHDF5 ( stageFileID, trim(thisName), &
                & myPrecisions, joinedQuantities )
              precQty => myPrecisions%quantities(1)
            else
              nullify ( precQty )
            end if
          end if

          select case ( get_spec_id ( storedResults(resInd)%key ) )
          case ( s_l2gp )
            call JoinL2GPQuantities ( storedResults(resInd)%key, hdfNameIndex, &
              & qty, precQty, l2gpDatabase, chunk, &
              & nameString=trim(storedResults(resInd)%hdfName))
          case ( s_l2aux )
            call JoinL2AuxQuantities ( storedResults(resInd)%key, hdfNameIndex, &
              & qty, l2auxDatabase, chunk, chunks )
            ! Ignore timing and direct writes
          end select

          ! Now destroy this vector.  We'll do this as we go along to make
          ! life easier for the computer.
          if ( parallel%stageInMemory ) then
            call DestroyVectorInfo ( joinedVectors(valInd) )
            call DestroyVectorTemplateInfo ( joinedVectorTemplates(valInd) )
            call DestroyQuantityTemplateContents ( joinedQuantities(valInd) )
            if ( storedResults(resInd)%gotPrecision ) then
              call DestroyVectorInfo ( joinedVectors(precInd) )
              call DestroyVectorTemplateInfo ( joinedVectorTemplates(precInd) )
              call DestroyQuantityTemplateContents ( joinedQuantities(precInd) )
            end if
          else
            call DestroyQuantityTemplateDatabase ( joinedQuantities )
            call DestroyVectorTemplateInfo ( myValues%template )
            call DestroyVectorInfo ( myValues )
            if ( storedResults(resInd)%gotPrecision ) then
              call DestroyVectorTemplateInfo ( myPrecisions%template )
              call DestroyVectorInfo ( myPrecisions )
            end if
          end if
        end if                          ! Didn't give up on this chunk
      end do
    end do

    ! Now clean up and quit
    if ( .not. parallel%stageInMemory ) then
      call H5FClose_F ( stageFileID, status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & 'Unable to close staging file:'//parallel%stagingFile )
    end if

    call DestroyStoredResultsDatabase ( storedResults )

    if ( associated ( joinedQuantities ) ) then
      deallocate ( joinedQuantities, STAT=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'joinedQuantities' )
    end if

    if ( associated ( joinedVectorTemplates ) ) then
      deallocate ( joinedVectorTemplates, STAT=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'joinedVectorTemplates' )
    end if

    if ( associated ( joinedVectors ) ) then
      deallocate ( joinedVectors, STAT=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'joinedVectors' )
    end if

    finished = .true.
    if ( index(switches,'mas') /= 0 ) then
      call output ( 'All chunks joined', advance='yes' )
    endif

    call Deallocate_test ( machineFree, 'machineFree', ModuleName )
    call Deallocate_test ( machineOK, 'machineOK', ModuleName )
    call Deallocate_test ( jobsMachineKilled, 'jobsMachineKilled', ModuleName )

  contains

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

  end subroutine L2MasterTask


  ! -------------------------------------------- AddStoredQuantityToDatabase ----
  integer function AddStoredResultToDatabase ( database, item )
    ! Add a storedQuantity to a database, or create the database if it
    ! doesn't yet exist

    ! Dummy arguments
    type (StoredResult_T), dimension(:), pointer :: database
    type (StoredResult_T), intent(in) :: item

    ! Local variables
    type (StoredResult_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddStoredResultToDatabase = newSize
  end function AddStoredResultToDatabase

  ! ----------------------------------------- DestroyStoredQuantitiesDatabase ---
  subroutine DestroyStoredResultsDatabase ( database )
    type (StoredResult_T), dimension(:), pointer :: database

    ! Local variabels
    integer :: i, status

    ! Executable code
    if ( .not. associated ( database ) ) return
    do i = 1, size(database)
      call deallocate_test ( database(i)%valInds, 'database(?)%valInds', ModuleName)
      if ( associated ( database(i)%precInds ) ) &
        & call deallocate_test ( database(i)%precInds,&
        & 'database(?)%precInds', ModuleName)
    end do
    deallocate ( database, STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'stored quantities database' )
  end subroutine DestroyStoredResultsDatabase

  ! -------------------------------------- StoreSlaveQuantity --------------
  subroutine StoreSlaveQuantity ( joinedQuantities, joinedVectorTemplates, &
    & joinedVectors, storedResults, chunk, noChunks, tid, noQuantitiesAccumulated,&
    & stageFileID )
    ! This routine reads a vector from a slave and stores it in
    ! an appropriate place.

    type (QuantityTemplate_T), dimension(:), pointer :: JOINEDQUANTITIES
    type (VectorTemplate_T), dimension(:), pointer :: JOINEDVECTORTEMPLATES
    type (Vector_T), dimension(:), pointer :: JOINEDVECTORS
    type (StoredResult_T), dimension(:), pointer :: STOREDRESULTS
    integer, intent(in) :: CHUNK        ! Index of chunk
    integer, intent(in) :: NOCHUNKS     ! Number of chunks
    integer, intent(in) :: TID          ! Slave tid
    integer, intent(inout) :: NOQUANTITIESACCUMULATED ! Running counter /index
    integer, intent(in) :: STAGEFILEID  ! For using staging file

    ! Local parameters
    integer, parameter :: DATABASEINFLATION = 500
    ! Local saved variables
    integer, parameter :: NINJOINPACKET   = 2 ! Num of ints in join packet

    ! Local variables
    integer, dimension(NINJOINPACKET) :: I2         ! Array to unpack
    integer :: INFO                                 ! Flag from pvm
    character(len=132) :: HDFNAME                   ! Output name for quantity
    character(len=132) :: GNAME                     ! HDF Group name for quantity
    integer :: KEY                                  ! Tree node
    integer :: GOTPRECISION                         ! Flag
    integer :: VIND                                 ! Index for vector
    integer :: PVIND                                ! Index for precision vector
    integer :: I                                    ! Loop inductor
    integer :: NOSTOREDRESULTS                      ! Array size
    logical :: SEENTHISBEFORE                       ! Flag
    logical :: NEEDTOINFLATE                        ! Flag
    real (r8), dimension(:,:), pointer :: VALUES    ! Values for this vector quantity
    character(len=1), dimension(:,:), pointer :: MASK ! Mask for this vector quantity
    type (QuantityTemplate_T) :: qt                 ! A quantity template
    type (VectorTemplate_T) :: vt                   ! A vector template
    type (Vector_T) :: v                            ! A vector
    type (StoredResult_T), pointer :: THISRESULT    ! Pointer to a storedResult
    type (StoredResult_T) :: ONERESULT              ! A new single storedResult_T

    ! Executable code

    ! First we unpack the rest of the information in the packet we were sent
    call PVMF90Unpack ( i2, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking integers from join packet" )
    call PVMIDLUnpack ( hdfName, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "unpacking hdf name from" )
    key = i2(1)
    gotPrecision = i2(2)

    ! Now get the quantity itself, possibly the precision and tropopause
    do i = 1, gotPrecision+1
      ! Increment our counter
      noQuantitiesAccumulated = noQuantitiesAccumulated + 1
      ! Get the quantity
      ! Just going to nullify qt, etc. I know it's intent in PVMRecieve..., but just
      ! to be sure.  Ditto for values and mask
      call NullifyQuantityTemplate ( qt )
      nullify ( values, mask )
      call PVMReceiveQuantity ( qt, values, mask=mask, justUnpack=.true. )

      ! Now we do different things depending on whether we're staging
      ! in memory or to a file
      if ( parallel%stageInMemory ) then
        ! Perhaps inflate our databases, some extra logic to avoig doing size on
        ! unassociated pointer
        needToInflate = noQuantitiesAccumulated == 1
        if ( .not. needToInflate ) &
          & needToInflate = noQuantitiesAccumulated > size ( joinedQuantities )
        if ( needToInflate ) then
          noQuantitiesAccumulated = &
            & InflateQuantityTemplateDatabase ( joinedQuantities, databaseInflation )
          noQuantitiesAccumulated = &
            & InflateVectorTemplateDatabase ( joinedVectorTemplates, databaseInflation )
          noQuantitiesAccumulated = &
            & InflateVectorDatabase ( joinedVectors, databaseInflation )
        end if
        
        ! Now add its template to our template database
        joinedQuantities ( noQuantitiesAccumulated ) = qt
        
        ! Now make a vector template up for this
        call ConstructVectorTemplate ( 0, joinedQuantities, &
          & (/ noQuantitiesAccumulated /), vt )
        joinedVectorTemplates ( noQuantitiesAccumulated ) = vt
        
        ! Now make a vector up for this
        v = CreateVector ( 0, joinedVectorTemplates( noQuantitiesAccumulated ), &
          & joinedQuantities, vectorNameText='joined' )
      else ! ------------------- Staging in a file

        call ConstructVectorTemplate ( 0, (/ qt /), (/ 1 /), vt )
        v = CreateVector ( 0, vt, (/ qt /), vectorNameText='joined', noValues=.true. )
      end if
      
      v%quantities(1)%values => values
      v%quantities(1)%mask => mask

      if ( parallel%stageInMemory ) then
        joinedVectors ( noQuantitiesAccumulated ) = v
      else
        ! Write this in our staging file
        write ( gName, '(i0)' ) noQuantitiesAccumulated
        call WriteVectorAsHDF5 ( stageFileID, gName, v )
        ! Destroy the information we got
        call DestroyQuantityTemplateContents ( qt )
        call DestroyVectorTemplateInfo ( vt )
        call DestroyVectorInfo ( v )
      end if

      select case ( i )
      case ( 1 )
        vInd = noQuantitiesAccumulated
      case ( 2 )
        pvInd = noQuantitiesAccumulated
      end select
    end do

    ! Now update our stored result stuff
    seenThisBefore = associated ( storedResults ) ! Perhaps, anyway
    if ( seenThisBefore ) seenThisBefore = any (storedResults%key == key )

    if ( seenThisBefore ) then
      thisResult => storedResults ( FindFirst ( storedResults%key == key ) )
    else
      ! We haven't seen this one before
      oneResult%key = key
      oneResult%gotPrecision = ( gotPrecision == 1 )
      oneResult%hdfName = hdfName
      call allocate_test ( oneResult%valInds, noChunks, &
        & 'valInds', ModuleName )
      if ( oneResult%gotPrecision ) &
        & call allocate_test ( oneResult%precInds, noChunks, &
        & 'precInds', ModuleName )
      ! Shouldn't need these really, but just in case
      oneResult%valInds = 0
      if ( oneResult%gotPrecision ) oneResult%precInds = 0
      noStoredResults = AddStoredResultToDatabase ( storedResults, oneResult )
      thisResult => storedResults ( noStoredResults )
    end if

    thisResult%valInds(chunk) = vInd
    if ( thisResult%gotPrecision ) thisResult%precInds(chunk) = pvInd

  end subroutine StoreSlaveQuantity

  ! --------------------------------- CleanUpDeadChunksOutput -------------
  subroutine CleanUpDeadChunksOutput ( chunk, joinedQuantities, &
    & joinedVectorTemplates, joinedVectors, storedResults )
    ! This destroys the vector info etc. associated with any output
    ! from a dead chunk.
    integer, intent(in) :: CHUNK        ! Index of chunk
    type (QuantityTemplate_T), dimension(:), pointer :: JOINEDQUANTITIES
    type (VectorTemplate_T), dimension(:), pointer :: JOINEDVECTORTEMPLATES
    type (Vector_T), dimension(:), pointer :: JOINEDVECTORS
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
        if ( parallel%stageInMemory ) then
          call DestroyVectorInfo ( joinedVectors(valInd) )
          call DestroyVectorTemplateInfo ( joinedVectorTemplates(valInd) )
          call DestroyQuantityTemplateContents ( joinedQuantities(valInd) )
        end if
        storedResults(res)%valInds(chunk) = 0
      end if

      if ( storedResults(res)%gotPrecision ) then
        precInd = storedResults(res)%precInds(chunk)
        if ( precInd /= 0 ) then
          if ( parallel%stageInMemory ) then
            call DestroyVectorInfo ( joinedVectors(precInd) )
            call DestroyVectorTemplateInfo ( joinedVectorTemplates(precInd) )
            call DestroyQuantityTemplateContents ( joinedQuantities(precInd) )
          end if
          storedResults(res)%precInds(chunk) = 0
        end if
      end if
    end do
  end subroutine CleanUpDeadChunksOutput

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module L2Parallel

!
! $Log$
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
