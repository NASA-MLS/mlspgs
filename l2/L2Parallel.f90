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
  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMF90PACK, PVMFMYTID, &
    & PVMF90UNPACK, PVMERRORMESSAGE, PVMTASKHOST, PVMFSPAWN, &
    & MYPVMSPAWN, PVMFCATCHOUT
  use PVMIDL, only: PVMIDLPACK
  use QuantityPVM, only: PVMSENDQUANTITY, PVMRECEIVEQUANTITY
  use MLSCommon, only: R8, MLSCHUNK_T
  use VectorsModule, only: VECTOR_T, VECTORVALUE_T, VECTORTEMPLATE_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE
  use L2GPData, only: L2GPDATA_T
  use L2AUXData, only: L2AUXDATA_T
  use QuantityTemplates, only: QuantityTemplate_T
  use MLSCommon, only: MLSCHUNK_T


  implicit none
  private

  public :: L2ParallelInfo_T, parallel, InitParallel, L2MasterTask
  public :: SlaveJoin, GetChunkFromMaster, CloseParallel
  
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This datatype defines configuration for the parallel code
  type L2ParallelInfo_T
    logical :: master = .false.         ! Set if this is a master task
    logical :: slave = .false.          ! Set if this is a slace task
    integer :: myTid                    ! My task ID in pvm
    integer :: masterTid                ! task ID in pvm
    character(len=132) :: slaveFilename ! Filename with list of slaves
  end type L2ParallelInfo_T

  ! Shared variables

  type (L2ParallelInfo_T), save :: parallel

  ! Private parameters

  integer, parameter :: CHUNKTAG   = 10
  integer, parameter :: INFOTAG    = ChunkTag + 1

  integer, parameter :: SIG_TOJOIN = 1
  integer, parameter :: SIG_FINISHED = SIG_toJoin + 1

  ! Parameters
  integer, parameter :: MACHINENAMELEN = 64 ! Name of machine

  ! Private types
  type StoredResult_T
    logical :: gotPrecision             ! If set have precision as well as value
    integer, dimension(:), pointer :: valInds ! Array vec. dtbs. inds (noChunks)
    integer, dimension(:), pointer :: precInds ! Array vec. dtbs. inds (noChunks)
  end type StoredResult_T

contains ! ================================ Procedures ======================

  ! ---------------------------------------------- InitParallel -------------
  subroutine InitParallel
    ! This routine initialises the parallel code
    ! Executable code
    if ( parallel%master .or. parallel%slave ) then
      call PVMFMyTid ( parallel%myTid )
    end if
  end subroutine InitParallel

  ! --------------------------------------------- CloseParallel -------------
  subroutine CloseParallel
    ! This routine closes down any parallel stuff
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM

    ! Exeuctable code
    if ( parallel%slave ) then
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMF90Pack ( (/ SIG_Finished /), info )
      if ( info /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to pack finished signal' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to send finish packet' )
    end if
  end subroutine CloseParallel

  ! --------------------------------------------- SendChunkToSlave ----------
  subroutine SendChunkToSlave ( chunk, chunkNo, slaveTid )
    ! This routine sends a chunk to a slave task

    ! Dummy arguments
    type (MLSChunk_T), intent(in) :: CHUNK ! The chunk
    integer, intent(in) :: CHUNKNO      ! Which one is it.
    integer, intent(in) :: SLAVETID     ! Who to send it to

    ! Local variables
    integer :: BUFFERID                 ! ID for buffer to send
    integer :: INFO                     ! Flag from PVM

    ! Executable code
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( &
      & (/ chunkNo, &
      &    chunk%firstMAFIndex, &
      &    chunk%lastMAFIndex, &
      &    chunk%noMAFsLowerOverlap, &
      &    chunk%noMAFsUpperOverlap, &
      &    chunk%accumulatedMAFs /), info )
    if ( info /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to pack chunk' )
    call PVMFSend ( slaveTid, ChunkTag, info )
    if ( info /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to send chunk' )
  end subroutine SendChunkToSlave

  ! ----------------------------------------------- GetChunkFromMaster ------
  subroutine GetChunkFromMaster ( chunks )
    ! This function gets a chunk sent by SendChunkToSlave

    ! Dummy arguments and result
    type (MLSChunk_T), dimension(:), pointer :: CHUNKS

    ! Local parameters
    integer, parameter :: noChunkTerms = 5 ! Number of components in MLSChunk_T

    ! Local variables
    integer :: BUFFERID                 ! ID for buffer for receive
    integer :: INFO                     ! Flag from PVM
    integer :: STATUS                   ! From allocate
    integer, dimension(noChunkTerms+1) :: VALUES ! Chunk as integer array

    ! Executable code
    if ( .not. parallel%slave ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Only slave tasks can receive a chunk' )

    call PVMFrecv ( parallel%masterTid, ChunkTAG, bufferID )
    if ( bufferID <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to receive chunk' )
    call PVMF90Unpack ( values, info )
    if ( info /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to unpack chunk')
    
    allocate ( chunks ( values(1):values(1) ), STAT=status)
    if ( status /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Failed to allocate one chunk' )

    chunks(lbound(chunks)) = MLSChunk_T ( &
      & values(2), values(3), values(4), values(5), values(6) )

  end subroutine GetChunkFromMaster

  ! ------------------------------------------- SlaveJoin ---------------
  subroutine SlaveJoin ( quantity, precisionQuantity, hdfName, key )
    ! This simply sends a vector quantity (or two) down a pvm spigot.
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
    call PVMF90Pack ( (/ SIG_ToJoin, key, gotPrecision /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing kind" )
    call PVMIDLPack ( hdfName, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing hdfName" )
    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to send join packet' )

    call PVMSendQuantity ( quantity, parallel%masterTid )
    if ( associated ( precisionQuantity ) ) &
      call PVMSendQuantity ( precisionQuantity, parallel%masterTid )
    
  end subroutine SlaveJoin

  ! ---------------------------------------- GetMachineNames ------------
  subroutine GetMachineNames ( machineNames )
    ! This reads the parallel slave name file and returns
    ! an array of machine names
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES

    ! Local variables
    integer :: LUN                      ! Logical unit number
    character(len=MachineNameLen) :: LINE ! A line from the file
    integer :: MACHINE                  ! Counter
    integer :: NOMACHINES               ! Array size
    logical :: EXIST                    ! Flag from inquire
    logical :: OPENED                   ! Flag from inquire
    integer :: STAT                     ! Status flag from read

    ! Executable code

    ! Find a free logical unit number
    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=parallel%slaveFilename, status='old', form='formatted', &
      & access='sequential', iostat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open slave file " // parallel%slaveFilename )

    ! Now read the file and count the lines
    noMachines = 0
    firstTimeRound: do
      read ( unit=lun, fmt=*, iostat=stat ) line
      if ( stat < 0 ) exit firstTimeRound
      noMachines = noMachines + 1
    end do firstTimeRound

    ! Now setup the result array
    allocate ( machineNames(noMachines), stat=stat )
    if ( stat /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//' machineNames')

    ! Now rewind the file and read the names
    rewind ( lun )
    do machine = 1, noMachines
      read ( unit=lun, fmt=* ) machineNames(machine)
    end do

    close ( unit=lun )
  end subroutine GetMachineNames

  ! -------------------------------------------- FindFirst --------------
  integer function FindFirst ( condition )
    ! Find the first logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION

    ! Local variables
    integer :: I                        ! Loop counter

    ! Executable code
    FindFirst = -1
    do i = 1, size(condition)
      if ( condition(i) ) then
        FindFirst = i
        return
      end if
    end do
  end function FindFirst

  ! --------------------------------------------- L2MasterTask ----------
  subroutine L2MasterTask ( chunks, l2gpDatabase, l2auxDatabase )
    ! This is a `master' task for the l2 software
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    type (L2GPData_T), dimension(:), pointer :: L2GPDATABASE
    type (L2AuxData_T), dimension(:), pointer :: L2AUXDATABASE

    ! Local parameter
    integer, parameter :: DELAY=200000  ! For Usleep, no. microsecs

    ! External (C) function
    external :: Usleep

    ! Local variables
    character(len=MachineNameLen), dimension(:), pointer :: MACHINENAMES
    character(len=132) :: COMMANDLINE

    integer :: BUFFERID                 ! From PVM
    integer :: BYTES                    ! Dummy from PVMFBufInfo
    integer :: INFO                     ! From PVM
    integer :: LUN                      ! File handle
    integer :: MACHINE                  ! Index
    integer :: MSGTAG                   ! Dummy from PVMFBufInfo
    integer :: NEXTCHUNK                ! Which chunk is next to be done
    integer :: NOCHUNKS                 ! Number of chunks
    integer :: NOMACHINES               ! Number of slaves
    integer :: SIGNAL                   ! From slave
    integer :: SLAVETID                 ! One slave
    integer :: TIDARR(1)                ! One tid

    integer, dimension(:), pointer :: SLAVECHUNKS ! Chunks for machines
    integer, dimension(:), pointer :: SLAVETIDS ! Taks ids for machines

    logical, dimension(:), pointer :: MACHINEFREE ! Is this machine busy
    logical, dimension(size(chunks)) :: COMPLETED

    logical, save :: FINISHED = .false. ! This will be called multiple times
    
    type (QuantityTemplate_T),dimension(:), pointer :: joinedQuantities 
    ! Local quantity template database
    type (VectorTemplate_T), dimension(:), pointer  :: joinedVectorTemplates
    ! Local vector template database
    type (Vector_T), dimension(:), pointer :: joinedVectors
    ! Local vector database
    
    ! Executable code --------------------------------

    ! First, is this the first call. The first call does all the work
    ! so if it's not then quit.
    if ( finished ) return

    ! Setup some stuff
    noChunks = size(chunks)
    nullify ( joinedQuantities, joinedVectorTemplates, joinedVectors, &
      & machineNames, slaveTids, slaveChunks, machineFree )

    ! Work out the information on our virtual machine
    call GetMachineNames ( machineNames )
    noMachines = size(machineNames)
    call Allocate_test ( slaveTids, noMachines, 'slaveTids', ModuleName )
    call Allocate_test ( slaveChunks, noMachines, 'slaveChunks', ModuleName )
    call Allocate_test ( machineFree, noMachines, 'machineFree', ModuleName )
    machineFree = .true.

    ! Loop until all chunks are done
    completed = .false.
    nextChunk = 1
    masterLoop: do
      ! This loop is in two main parts.

      ! In the first part, we look to see if there are any chunks still to be
      ! done, and any vacant machines to do them.
      if ( nextChunk < noChunks .and. any(machineFree) ) then 
        machine = FindFirst(machineFree)
        commandLine = 'mlsl2'

        print*,'Calling catchout'
        call PVMFCatchOut ( 1, info )
        if ( info /= 0 ) call PVMErrorMessage ( info, "calling catchout" )
        print*,'Launching '//trim(commandLine)//' on '//trim(machineNames(machine))
        info = myPVMSpawn ( trim(commandLine), PvmTaskHost, &
          & trim(machineNames(machine)), &
          & 1, tidarr )
        slaveTids(machine) = tidArr(1)
        if ( info < 0 ) call PVMErrorMessage ( info, "launching slave" )
        if ( info == 0 ) call PVMErrorMessage ( slaveTids(machine), "launching slave" )
        call SendChunkToSlave ( chunks(nextChunk), nextChunk, &
          & slaveTids(machine) )
        machineFree(machine) = .false.
        slaveChunks(machine) = nextChunk
        nextChunk = nextChunk + 1
      endif

      ! In this next part, we listen out for communication from the slaves and
      ! process it accordingly.
      call PVMFNRecv( -1, InfoTag, bufferID )
      if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for Info message" )
      else if ( bufferID > 0 ) then
        ! Who sent this?
        call PVMFBufInfo ( bufferID, bytes, msgTag, slaveTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        machine = FindFirst ( slaveTids == slaveTid )
        if ( machine == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Got a message from an unknown slave")
        ! Unpack the first integer in the buffer
        call PVMF90Unpack ( signal, info )
        if ( info /= 0 ) then
          call PVMErrorMessage ( info, "unpacking signal" )
        else
          select case (signal) 
          case ( sig_tojoin )
          case ( sig_finished )
            completed(slaveChunks(machine)) = .true.
            machineFree(machine) = .true.
          case default
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Unkown signal from slave' )
          end select
        end if
      end if
      
      ! If we're done then exit
      if (all(completed)) exit masterLoop

      ! Now, rather than chew up cpu time on the master process, we'll wait a
      ! bit here.
      call usleep ( delay )

    end do masterLoop

    ! Now we join up all our results into l2gp and l2aux quantities

    finished = .true.

  end subroutine L2MasterTask

end module L2Parallel




