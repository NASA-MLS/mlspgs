! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module L2ParInfo
  ! This module provides definitions needed by L2Parallel and other modules to
  ! manage the parallel aspects of the L2 code.

  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use PVM, only: PVMFMYTID, PVMFINITSEND, PVMF90PACK, PVMFSEND, &
    & PVMDATADEFAULT, PVMERRORMESSAGE, PVMF90UNPACK, NEXTPVMARG
  use PVMIDL, only: PVMIDLPACK
  use VectorsModule, only: VECTORVALUE_T
  use QuantityPVM, only: PVMSENDQUANTITY
  use MLSStrings, only: LowerCase

  implicit none
  private

  public :: L2ParallelInfo_T, parallel, InitParallel, CloseParallel
  public :: SIG_ToJoin, SIG_Finished, SIG_Register, ChunkTag, InfoTag, SlaveJoin
  public :: SIG_AckFinish, SIG_RequestDirectWrite, SIG_DirectWriteGranted
  public :: SIG_DirectWriteFinished
  public :: NotifyTag, GetNiceTidString, SlaveArguments
  public :: AccumulateSlaveArguments, RequestDirectWritePermission
  public :: FinishedDirectWrite
  
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
  integer, parameter :: NOTIFYTAG  = InfoTag + 1

  integer, parameter :: SIG_TOJOIN = 1
  integer, parameter :: SIG_FINISHED = SIG_toJoin + 1
  integer, parameter :: SIG_ACKFINISH = SIG_finished + 1
  integer, parameter :: SIG_REGISTER = SIG_AckFinish + 1
  integer, parameter :: SIG_REQUESTDIRECTWRITE = SIG_Register + 1
  integer, parameter :: SIG_DIRECTWRITEGRANTED = SIG_RequestDirectWrite + 1
  integer, parameter :: SIG_DIRECTWRITEFINISHED = SIG_DirectWriteGranted + 1

  ! This datatype defines configuration for the parallel code
  type L2ParallelInfo_T
    logical :: master = .false.         ! Set if this is a master task
    logical :: slave = .false.          ! Set if this is a slace task
    integer :: myTid                    ! My task ID in pvm
    integer :: masterTid                ! task ID in pvm
    character(len=132) :: slaveFilename ! Filename with list of slaves
    character(len=132) :: executable    ! Executable filename
    character(len=132) :: submit=""     ! Submit comand for batch queue system
    integer :: maxFailuresPerMachine = 1 ! More than this then don't use it
    integer :: maxFailuresPerChunk = 4 ! More than this then give up on getting it
  end type L2ParallelInfo_T

  ! Shared variables

  type (L2ParallelInfo_T), save :: parallel
  character ( len=2048 ), save :: slaveArguments = ""

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
  subroutine InitParallel ( chunkNo )
    ! This routine initialises the parallel code
    integer, intent(in) :: chunkNo      ! Chunk number asked to do.
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Executable code
    if ( parallel%master .or. parallel%slave ) then
      call PVMFMyTid ( parallel%myTid )
    end if
    if ( parallel%slave ) then
      ! Register ourselves with the master
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMF90Pack ( SIG_Register, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing registration' )
      call PVMF90Pack ( chunkNo, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'packing chunkNumber' )
      call PVMFSend ( parallel%masterTid, InfoTag, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, 'sending finish packet' )
    end if
  end subroutine InitParallel

  ! --------------------------------------------- CloseParallel -------------
  subroutine CloseParallel
    ! This routine closes down any parallel stuff
    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM

    integer :: SIGNAL                   ! From acknowledgement packet

    ! Exeuctable code
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
    end if
  end subroutine CloseParallel

  ! --------------------------------------- FinishedDirectWrite ------------
  subroutine FinishedDirectWrite
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    ! Local variables
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_DirectWriteFinished, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write finished flag" )

    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write finished packet" )
    
    ! Executable code
  end subroutine FinishedDirectWrite
    

  ! ---------------------------------------- RequestDirectWritePermission --
  subroutine RequestDirectWritePermission (filename, createFile )
    integer, intent(in) :: FILENAME
    logical, intent(out) :: CREATEFILE

    ! Local variables
    integer :: BUFFERID                 ! From PVM
    integer :: INFO                     ! From PVM
    integer :: SIGNAL                   ! From Master
    integer :: CREATEFLAG               ! Returned by PVM

    ! Executable code

    ! Pack and dispatch
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMF90Pack ( SIG_RequestDirectWrite, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write request flag" )
    call PVMF90Pack ( filename, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "packing direct write information" )

    call PVMFSend ( parallel%masterTid, InfoTag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, "sending direct write request packet" )

    ! Now wait for a reply.  This may mean waiting for other chunks to
    ! finish writing their part of the file.
    call PVMFRecv ( parallel%masterTid, InfoTag, bufferID )
    if ( bufferID <= 0 ) call PVMErrorMessage ( bufferID, &
      & 'receiving direct write permission' )

    ! Once we have the reply unpack it
    call PVMF90Unpack ( signal, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking direct write permission')
    if ( signal /= SIG_DirectWriteGranted ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Got unrecognised signal from master' )
    call PVMF90Unpack ( createFlag, info )
    if ( info /= 0 ) &
      & call PVMErrorMessage ( info, 'unpacking create file flag')
    createFile = createFlag == 1
  end subroutine RequestDirectWritePermission

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

end module L2ParInfo

! $Log$
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
