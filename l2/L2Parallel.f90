! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2Parallel
  ! This module contains low level routines and stuff for dealing with parallel
  ! invocations of the MLSL2 program.

  ! Level 2 programs can be masters or slaves, or neither, not both.  A task
  ! which is neither simply runs through the l2cf as normal.  A master task
  ! does the chunk divide and then launches slave tasks for each chunk, and
  ! awaits their results in the Join section.

  use PVM, only: PVMDATADEFAULT, PVMFINITSEND, PVMF90PACK, PVMFMYTID, &
    & PVMF90UNPACK
  use MLSCommon, only: R8, MLSCHUNK_T
  use VectorsModule, only: VECTOR_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none
  private

  public :: L2ParallelInfo_T, parallel, InitParallel, GetChunkFromMaster
  
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
  end type L2ParallelInfo_T

  ! Private parameters

  integer, parameter :: ChunkTag = 10
  integer, parameter :: VectorTag = ChunkTag + 1

  ! Shared variables

  type (L2ParallelInfo_T), save :: parallel

contains ! ================================ Procedures ======================

  ! ---------------------------------------------- InitParallel -------------
  subroutine InitParallel
    ! This routine initialises the parallel code
    ! Executable code
    if ( parallel%master .or. parallel%slave ) then
      call PVMFMyTid ( parallel%myTid )
    end if
  end subroutine InitParallel

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

end module L2Parallel




