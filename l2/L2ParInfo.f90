! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2ParInfo
  ! This module provides definitions needed by L2Parallel and other modules to
  ! manage the parallel aspects of the L2 code.

  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use PVM, only: PVMFMYTID, PVMFINITSEND, PVMF90PACK, PVMFSEND, &
    & PVMDATADEFAULT, PVMERRORMESSAGE
  use PVMIDL, only: PVMIDLPACK
  use VectorsModule, only: VECTORVALUE_T
  use QuantityPVM, only: PVMSENDQUANTITY

  implicit none
  private

  public :: L2ParallelInfo_T, parallel, InitParallel, CloseParallel
  public :: SIG_ToJoin, SIG_Finished, ChunkTag, InfoTag, SlaveJoin
  
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Parameters

  integer, parameter :: CHUNKTAG   = 10
  integer, parameter :: INFOTAG    = ChunkTag + 1

  integer, parameter :: SIG_TOJOIN = 1
  integer, parameter :: SIG_FINISHED = SIG_toJoin + 1

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

contains ! ==================================================================

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
      call PVMF90Pack ( SIG_Finished, info )
      if ( info /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to pack finished signal' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to send finish packet' )
    end if
  end subroutine CloseParallel

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
    call PVMF90Pack ( SIG_ToJoin, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing kind" )
    call PVMF90Pack ( (/key, gotPrecision /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing key, gotPrecision" )
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

end module L2ParInfo
