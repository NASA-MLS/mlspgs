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
    & GETMACHINENAMEFROMTID
  use PVMIDL, only: PVMIDLPACK, PVMIDLUNPACK
  use QuantityPVM, only: PVMSENDQUANTITY, PVMRECEIVEQUANTITY
  use MLSCommon, only: R8, MLSCHUNK_T, FINDFIRST
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VECTORTEMPLATE_T, &
    & ADDVECTORTEMPLATETODATABASE, CONSTRUCTVECTORTEMPLATE, ADDVECTORTODATABASE, &
    & CREATEVECTOR, DESTROYVECTORINFO, DESTROYVECTORTEMPLATEINFO
  use Machine, only: SHELL_COMMAND
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE, &
    & MLSMSG_Deallocate
  use L2GPData, only: L2GPDATA_T
  use L2AUXData, only: L2AUXDATA_T
  use L2ParInfo, only: L2PARALLELINFO_T, PARALLEL, INFOTAG, CHUNKTAG, &
    & SIG_TOJOIN, SIG_FINISHED, SIG_ACKFINISH, SIG_REGISTER, NOTIFYTAG, &
    & SIG_REQUESTDIRECTWRITE, SIG_DIRECTWRITEGRANTED, SIG_DIRECTWRITEFINISHED, &
    & GETNICETIDSTRING, SLAVEARGUMENTS
  use QuantityTemplates, only: QUANTITYTEMPLATE_T, ADDQUANTITYTEMPLATETODATABASE, &
    & DESTROYQUANTITYTEMPLATECONTENTS
  use Toggles, only: Gen, Switches, Toggle
  use Output_m, only: Output
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_STRING
  use String_table, only: Display_String
  use Init_Tables_Module, only: S_L2GP, S_L2AUX
  use MoreTree, only: Get_Spec_ID

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
  !---------------------------------------------------------------------------

  ! Parameters
  integer, parameter :: MACHINENAMELEN = 64 ! Max length of name of machine
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
        &    chunks(chunk)%accumulatedMAFs /), info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing one chunk' )
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
      chunks(chunk) = MLSChunk_T ( &
        & values(1), values(2), values(3), values(4), values(5) )
    end do

    if ( index(switches,'chu') /=0 ) call dump ( chunks )

  end subroutine GetChunkInfoFromMaster

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
      do lun = 20, 99
        inquire ( unit=lun, exist=exist, opened=opened )
        if ( exist .and. .not. opened ) exit
      end do
      if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
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

  ! --------------------------------------------- L2MasterTask ----------
  subroutine L2MasterTask ( chunks, l2gpDatabase, l2auxDatabase )
    ! This is a `master' task for the l2 software
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    type (L2AuxData_T), dimension(:), pointer :: L2AUXDATABASE

    ! Local parameter
    integer, parameter :: DELAY = 200000  ! For Usleep, no. microsecs
    integer, parameter :: MAXDIRECTWRITEFILES=200 ! For internal array sizing

    ! External (C) function
    external :: Usleep

    ! Local variables
    logical :: USINGSUBMIT              ! Set if using the submit mechanism
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    character(len=MachineNameLen) :: THISNAME
    character(len=8) :: CHUNKNOSTR
    character(len=2048) :: COMMANDLINE

    integer :: BUFFERID                 ! From PVM
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    integer :: CHUNK                    ! Loop counter
    integer :: COMPLETEDFILE            ! String index from slave
    logical :: CREATEFILE               ! Flag for direct writes
    integer :: DEADCHUNK                ! A chunk from a dead task
    integer :: DEADMACHINE              ! A machine for a dead task
    integer :: DEADTID                  ! A task that's died
    integer :: HDFNAMEINDEX             ! String index
    integer :: INFO                     ! From PVM
    integer :: MACHINE                  ! Index
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: NEXTCHUNK                ! A chunk number
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: NODIRECTWRITEFILES       ! Need to keep track of filenames
    integer :: NOMACHINES               ! Number of slaves
    integer :: PENDINGCHUNK             ! Chunk waiting for direct write permission
    integer :: PRECIND                  ! Array index
    integer :: RESIND                   ! Loop counter
    integer :: REQUESTEDFILE            ! String index from slave
    integer :: SIGNAL                   ! From slave
    integer :: SLAVETID                 ! One slave
    integer :: STATUS                   ! From deallocate
    integer :: TIDARR(1)                ! One tid
    integer :: VALIND                   ! Array index

    integer, dimension(size(chunks)) :: CHUNKMACHINES ! Machine indices for chunks
    integer, dimension(size(chunks)) :: CHUNKTIDS ! Tids for chunks

    integer, dimension(size(chunks)) :: DIRECTWRITESTATUS
    integer, dimension(maxDirectWriteFiles) :: DIRECTWRITEFILES

    logical, dimension(:), pointer :: MACHINEFREE ! Is this machine busy
    logical, dimension(:), pointer :: MACHINEOK ! Is this machine working?
    integer, dimension(:), pointer :: JOBSMACHINEKILLED ! Counts failures per machine.
    logical, dimension(size(chunks)) :: CHUNKSCOMPLETED ! Chunks completed
    logical, dimension(size(chunks)) :: CHUNKSSTARTED ! Chunks being processed
    logical, dimension(size(chunks)) :: CHUNKSABANDONED ! Chunks kept failing
    integer, dimension(size(chunks)) :: CHUNKFAILURES ! Failure count

    logical, save :: FINISHED = .false. ! This will be called multiple times
    
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
    
    ! Executable code --------------------------------

    ! First, is this the first call. The first call does all the work
    ! so if it's not then quit.
    if ( finished ) return

    usingSubmit = trim(parallel%submit) /= ''

    ! Setup some stuff
    noChunks = size(chunks)
    nullify ( joinedQuantities, joinedVectorTemplates, joinedVectors, &
      & machineNames, machineFree, storedResults, machineOK, jobsMachineKilled )
    directWriteFiles = 0
    directWriteStatus = 0

    ! Work out the information on our virtual machine
    call GetMachineNames ( machineNames )
    noMachines = size(machineNames)
    call Allocate_test ( machineFree, noMachines, 'machineFree', ModuleName )
    call Allocate_test ( machineOK, noMachines, 'machineOK', ModuleName )
    call Allocate_test ( jobsMachineKilled, noMachines, 'jobsMachineKilled', ModuleName )
    machineFree = .true.
    machineOK = .true.
    jobsMachineKilled = 0

    ! Loop until all chunks are done
    chunksCompleted = .false.
    chunksStarted = .false.
    chunksAbandoned = .false.
    chunkFailures = 0
    chunkTids = 0
    chunkMachines = 0

    masterLoop: do ! --------------------------- Master loop -----------------------
      ! This loop is in two main parts.

      ! In the first part, we look to see if there are any chunks still to be
      ! started done, and any vacant machines to do them.
      ! --------------------------------------------------------- Start new jobs? --
      if ( (.not. all(chunksStarted .or. chunksAbandoned)) .and. &
        & ( any(machineFree .and. machineOK) .or. usingSubmit ) ) then
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
            machineOK(machine) = .false.
          end if
        end if
      end if

      ! --------------------------------------------------- Messages from jobs?
      ! In this next part, we listen out for communication from the slaves and
      ! process it accordingly.
      receiveInfoLoop: do
        call PVMFNRecv( -1, InfoTag, bufferID )
        if ( bufferID == 0 ) exit receiveInfoLoop
        if ( bufferID < 0 ) then
          call PVMErrorMessage ( info, "checking for Info message" )
        else if ( bufferID > 0 ) then
          ! Who sent this?
          call PVMFBufInfo ( bufferID, bytes, msgTag, slaveTid, info )
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
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "Got a message from an unknown slave")
          else
            ! Unpack the first integer in the buffer
            if ( signal /= sig_register ) machine = chunkMachines(chunk)
          end if

          select case (signal) 
            
          case ( sig_register ) ! ----------------- Chunk registering ------
            call PVMF90Unpack ( chunk, info )
            if ( info /= 0 ) then
              call PVMErrorMessage ( info, "unpacking chunk number" )
            endif
            if ( usingSubmit ) then
              ! We only really care about this message if we're using submit
              chunkTids(chunk) = slaveTid
              call GetMachineNameFromTid ( slaveTid, thisName )
              machine = FindFirst ( trim(thisName) == machineNames )
              if ( machine == 0 ) &
                & call MLSMessage ( MLSMSG_Error, ModuleName, &
                & "An unknown machine sent a registration message!" )
              chunkMachines(chunk) = machine
              call WelcomeSlave ( chunk, slaveTid )
              if ( index(switches,'mas') /= 0 ) then
                call output ( 'Welcomed task ' // &
                  & trim(GetNiceTidString(slaveTid)) // &
                  & ' running chunk ' )
                call output ( chunk )
                call output ( ' on machine ' // &
                  & trim(machineNames(machine)), advance='yes' )
              end if
            endif
            
          case ( sig_tojoin ) ! --------------- Got a join request ---------
            call StoreSlaveQuantity( joinedQuantities, &
              & joinedVectorTemplates, joinedVectors, &
              & storedResults, chunk, noChunks, slaveTid )

          case ( sig_RequestDirectWrite ) ! ------- Direct write permission --
            ! What file did they ask for?
            call PVMF90Unpack ( requestedFile, info )
            if ( info /= 0 )  call PVMErrorMessage ( info, &
              & "unpacking direct write request" )
            ! Is this a new file?
            createFile = .not. any ( directWriteFiles == requestedFile )
            if ( createFile ) then
              noDirectWriteFiles = noDirectWriteFiles + 1
              if ( noDirectWriteFiles > maxDirectWriteFiles ) &
                & call MLSMessage ( MLSMSG_Error, ModuleName, &
                & 'Too many direct write files, increase limit' )
              directWriteFiles ( noDirectWriteFiles ) = requestedFile
            end if
            if ( index ( switches, 'mas' ) /= 0 ) then
              call output ( 'Direct write request for file ' )
              call output ( requestedFile )
              call output ( ' from ' // &
                & trim(machineNames(machine)) // ' ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' chunk ' )
              call output ( chunk, advance='yes')
            end if
            ! Is anyone else using this file?
            if ( any ( directWriteStatus == requestedFile ) ) then
              ! If so, log a request for it by setting our status to 
              ! -requestedFile
              if ( index ( switches, 'mas' ) /= 0 ) then
                call output ( 'Request pending', &
                  & advance='yes' )
              end if
              directWriteStatus(chunk) = -requestedFile
            else
              ! Otherwise, go ahead
              if ( index ( switches, 'mas' ) /= 0 ) then
                call output ( 'Request was granted' )
                if ( createFile ) then
                  call output ( ' (new file)', advance='yes' )
                else
                  call output ( ' (old file)', advance='yes' )
                end if
              end if
              call GrantDirectWrite ( slaveTid, createFile )
              directWriteStatus(chunk) = requestedFile
            end if

          case ( sig_DirectWriteFinished ) ! - Finished with direct write -
            ! Record that the chunk has finished direct write
            completedFile = directWriteStatus(chunk)
            directWriteStatus(chunk) = 0
            if ( index ( switches, 'mas' ) /= 0 ) then
              call output ( 'Direct write finished on file ' )
              call output ( completedFile )
              call output ( ' by ' // &
                & trim(machineNames(machine)) // ' ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' chunk ' )
              call output ( chunk, advance='yes')
            end if
            ! Is anyone else waiting for this file?
            pendingChunk = FindFirst ( directWriteStatus == -completedFile )
            if ( pendingChunk /= 0 ) then
              ! We know createFile is false here.
              if ( index ( switches, 'mas' ) /= 0 ) then
                call output ( 'Permission granted to ' // &
                  & trim(GetNiceTidString(chunkTids(pendingChunk))) // &
                  & ' chunk ' )
                call output ( pendingChunk, advance='yes' )
              end if
              call GrantDirectWrite ( chunkTids(pendingChunk), .false. )
              directWriteStatus ( pendingChunk ) = completedFile
            end if

          case ( sig_finished ) ! -------------- Got a finish message ----
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'Got a finished message from ' // &
                & trim(machineNames(machine)) // ' ' // &
                & trim(GetNiceTidString(slaveTid)) // &
                & ' processing chunk ' )
              call output ( chunk, advance='yes')
            endif
            
            ! Send an acknowledgement
            call PVMFInitSend ( PVMDataDefault, bufferID )
            if ( bufferId < 0 ) &
              & call PVMErrorMessage ( bufferID, 'setting up finish ack.' )
            call PVMF90Pack ( SIG_AckFinish, info )
            if ( info /= 0 ) &
              & call PVMErrorMessage ( info, 'packing finish ack.' )
            call PVMFSend ( slaveTid, InfoTag, info )
            if ( info /= 0 ) &
              & call PVMErrorMessage ( info, 'sending finish ack.' )
            
            ! Now update our information
            chunksCompleted(chunk) = .true.
            if ( .not. usingSubmit ) machineFree(machine) = .true.
            chunkTids(chunk) = 0
            chunkMachines(chunk) = 0
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
        end if
      end do receiveInfoLoop

      ! Listen out for any message that a slave task has died
      call PVMFNRecv ( -1, NotifyTAG, bufferID )
      if ( bufferID > 0 ) then
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
            deadMachine = chunkMachines(deadChunk)
            if ( .not. usingSubmit ) machineFree(deadMachine) = .true.
            if ( index(switches,'mas') /= 0 ) then
              call output ( 'The run of chunk ' )
              call output ( deadChunk )
              call output ( ' on ' // trim(machineNames(deadMachine)) // &
                & ' ' // trim(GetNiceTidString(deadTid)) // &
                & ' died, try again.', advance='yes' )
            end if
            call CleanUpDeadChunksOutput ( deadChunk, joinedQuantities, &
              & joinedVectorTemplates, joinedVectors, storedResults )
            chunksStarted(deadChunk) = .false.
            chunkFailures(deadChunk) = chunkFailures(deadChunk) + 1
            directWriteStatus(deadChunk) = 0
            where ( machineNames(deadMachine) == machineNames )
              jobsMachineKilled = jobsMachineKilled + 1
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
            if ( jobsMachineKilled(deadMachine) > &
              & parallel%maxFailuresPerMachine .and. .not. usingSubmit) then
              if ( index(switches,'mas') /= 0 ) &
                & call output ('The machine ' // &
                & trim(machineNames(deadMachine)) // &
                & ' keeps killing things, marking it bad', &
                & advance='yes' )
              where ( machineNames(deadMachine) == machineNames )
                machineOK = .false.
              end where
            end if
          else
            if ( index(switches,'mas') /= 0 ) call output ( &
              & "A slave task died after giving results, " // &
              & "we won't worry about it.", &
              & advance='yes' )
          end if
        end if
      else if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for Notify message" )
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
      if ( .not. any(machineOK) .and. &
        &  all ( chunksStarted .eqv. chunksCompleted ) ) then
        if ( index(switches,'mas') /= 0 ) &
          & call output ( 'No machines left to do the remaining work', &
          & advance='yes' )
        exit masterLoop
      end if

      ! Now, rather than chew up cpu time on the master machine, we'll wait a
      ! bit here.
      call usleep ( delay )
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
    do resInd = 1, size ( storedResults )
      hdfNameIndex = enter_terminal ( trim(storedResults(resInd)%hdfName), t_string )
      do chunk = 1, noChunks

        if (.not. chunksAbandoned(chunk) ) then
          ! Setup for this quantity
          valInd = storedResults(resInd)%valInds(chunk)
          qty => joinedVectors(valInd)%quantities(1)
          
          if ( storedResults(resInd)%gotPrecision ) then
            precInd = storedResults(resInd)%precInds(chunk)
            precQty => joinedVectors(precInd)%quantities(1)
          else
            nullify ( precQty )
          endif
          
          if ( index(switches,'mas') /= 0 .and. chunk==1 ) then
            call output ( 'Joining ' )
            call display_string ( qty%template%name, advance='yes' )
            call output ( 'Minor frame:' )
            call output ( qty%template%minorFrame, advance='yes' )
          endif
          
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
          call DestroyVectorInfo ( joinedVectors(valInd) )
          call DestroyVectorTemplateInfo ( joinedVectorTemplates(valInd) )
          call DestroyQuantityTemplateContents ( joinedQuantities(valInd) )

          if ( storedResults(resInd)%gotPrecision ) then
            call DestroyVectorInfo ( joinedVectors(precInd) )
            call DestroyVectorTemplateInfo ( joinedVectorTemplates(precInd) )
            call DestroyQuantityTemplateContents ( joinedQuantities(precInd) )
          end if

        end if                          ! Didn't give up on this chunk
      end do
    end do
    
    ! Now clean up and quit
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
      finished = .true.
    end if

    if ( index(switches,'mas') /= 0 ) then
      call output ( 'All chunks joined', advance='yes' )
    endif

    call Deallocate_test ( machineFree, 'machineFree', ModuleName )
    call Deallocate_test ( machineOK, 'machineOK', ModuleName )
    call Deallocate_test ( jobsMachineKilled, 'jobsMachineKilled', ModuleName )

  contains

    subroutine GrantDirectWrite ( tid, createFile )
      ! This routine sends a direct write grant message to a slave
      integer, intent(in) :: TID        ! Slave tid
      logical, intent(in) :: CREATEFILE ! Set if new file
      ! Local variables
      integer :: INFO                   ! From PVM
      integer :: BUFFERID               ! From PVM
      integer :: CREATE                 ! Integer version of createFile
      ! Executable code
      create = 0
      if ( createFile ) create = 1
      call PVMFInitSend ( PvmDataDefault, bufferID ) 
      call PVMF90Pack ( (/ SIG_DirectWriteGranted, create /), info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing direct write permission' )
      call PVMFSend ( tid, InfoTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending direct write permission' )
    end subroutine GrantDirectWrite

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

  ! -------------------------------------- StoreSlaveQuantity -------------------
  subroutine StoreSlaveQuantity ( joinedQuantities, joinedVectorTemplates, &
    & joinedVectors, storedResults, chunk, noChunks, tid )
    ! This routine reads a vector from a slave and stores it in
    ! an appropriate place.

    type (QuantityTemplate_T), dimension(:), pointer :: JOINEDQUANTITIES
    type (VectorTemplate_T), dimension(:), pointer :: JOINEDVECTORTEMPLATES
    type (Vector_T), dimension(:), pointer :: JOINEDVECTORS
    type (StoredResult_T), dimension(:), pointer :: STOREDRESULTS
    integer, intent(in) :: CHUNK        ! Index of chunk
    integer, intent(in) :: NOCHUNKS     ! Number of chunks
    integer, intent(in) :: TID          ! Slave tid

    ! Local saved variables
    integer, save      :: JOINEDQTCOUNTER = CounterStart ! To place in qt%id
    integer, save      :: JOINEDVTCOUNTER = CounterStart ! To place in vt%id
    integer, parameter :: NINJOINPACKET   = 2 ! Num of ints in join packet

    ! Local variables
    integer, dimension(NINJOINPACKET) :: I2         ! Array to unpack
    integer :: INFO                                 ! Flag from pvm
    character(len=132) :: HDFNAME                   ! HDF name for quantity
    integer :: KEY                                  ! Tree node
    integer :: GOTPRECISION                         ! Flag
    integer :: QTIND                                ! Index for quantity template
    integer :: VTIND                                ! Index for vector template
    integer :: VIND                                 ! Index for vector
    integer :: PVIND                                ! Index for precision vector
    integer :: I                                    ! Loop inductor
    integer :: NOSTOREDRESULTS                      ! Array size
    logical :: SEENTHISBEFORE                       ! Flag
    real (r8), dimension(:,:), pointer :: VALUES    ! Values for this vector quantity
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
      call PVMReceiveQuantity ( qt, values, justUnpack=.true. )
      
      ! Now add its template to our template database
      qt%id = joinedQTCounter
      joinedQTCounter = joinedQTCounter + 1
      qtInd = AddQuantityTemplateToDatabase ( joinedQuantities, qt )
      
      ! Now make a vector template up for this
      call ConstructVectorTemplate ( 0, joinedQuantities, (/ qtInd /), vt )
      vt%id = joinedVTCounter
      joinedVTCounter = joinedVTCounter + 1
      vtInd = AddVectorTemplateToDatabase ( joinedVectorTemplates, vt )
      
      ! Now make a vector up for this
      v = CreateVector ( 0, joinedVectorTemplates(vtInd), &
        & joinedQuantities, vectorNameText='joined' )
      v%quantities(1)%values => values
      select case ( i )
      case ( 1 )
        vInd = AddVectorToDatabase ( joinedVectors, v )
      case ( 2 )
        pvInd = AddVectorToDatabase ( joinedVectors, v )
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
        call DestroyVectorInfo ( joinedVectors(valInd) )
        call DestroyVectorTemplateInfo ( joinedVectorTemplates(valInd) )
        call DestroyQuantityTemplateContents ( joinedQuantities(valInd) )
        storedResults(res)%valInds(chunk) = 0
      end if

      if ( storedResults(res)%gotPrecision ) then
        precInd = storedResults(res)%precInds(chunk)
        if ( precInd /= 0 ) then
          call DestroyVectorInfo ( joinedVectors(precInd) )
          call DestroyVectorTemplateInfo ( joinedVectorTemplates(precInd) )
          call DestroyQuantityTemplateContents ( joinedQuantities(precInd) )
          storedResults(res)%precInds(chunk) = 0
        end if
      end if
    end do
  end subroutine CleanUpDeadChunksOutput

end module L2Parallel

!
! $Log$
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
