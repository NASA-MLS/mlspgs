! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module SnoopMLSL2               ! Interface between MLSL2 and IDL snooper via pvm3.

  ! This module is the fortran end of an interaction between the MLS level 2
  ! programmer and a user running an interactive IDL program to diagnose the
  ! progress and behavior of the MLS level 2 program.  The idea is that one can
  ! place a call to a Snoop routine, passing it any subset of relevant
  ! information (e.g. vectors and matrix databases). The software will look for
  ! an IDL program ready to accept snoop requests and if one is around will do
  ! as instructed.

  ! This will make heavy use of some low level routines in the PVM and PVMIDL
  ! modules to do much of the communication.

  ! Also note that multiple IDL snoopers can talk to one or many f90 procedures,
  ! the little extra book keeping this involves is worth it.

  use VectorsModule, only: Vector_T, VectorValue_T
  use QuantityTemplates, only: QuantityTemplate_T
  use MatrixModule_1, ONLY: Matrix_T

  use Intrinsic, only: LIT_INDICES
  use MLSCommon, only: R4, R8, I4
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
    MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate
  use PVM, only: PVMDataDefault, PVMFinitsend, PVMFmyTid, PVMFgSize, &
    & PVMErrorMessage, PVMF90Unpack, PVMTaskExit, PVMFNotify, PVMFSend
  use PVMIDL, only:  IDLMsgTag, PVMIDLPack, PVMIDLReceive, PVMIDLSend, PVMIDLUnpack
  use QuantityPVM, only: PVMSENDQUANTITY
  use TREE, only:  DUMP_TREE_NODE, SUB_ROSA, SUBTREE, SOURCE_REF
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use LEXER_CORE, only: PRINT_SOURCE
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_IDENTIFIER
  use MLSSignals_M, only: GETSIGNALNAME
  use MLSCommon, only: FINDFIRST

  implicit none
  private

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! =============================================================================

  ! The first main thing is a data type that describes each of the active
  ! snoopers.

  ! Part of this the following enumerated type, describing the `mode' of the
  ! snooper.  The mode can be `observing', or `controling'. In observing mode
  ! the snooper can only request copies of vectors etc.  In controling mode,
  ! the snooper can write values back for these and to some extent control the
  ! program (give up on iterations, exit etc.)

  ! Note that snoopers communicate amongst themselves to arbitrate requests to
  ! become controling, the succuesful one contacts the l2 code.

  integer, parameter :: SnooperObserving =            0
  integer, parameter :: SnooperControling =           SnooperObserving + 1

  integer, parameter :: SnoopTag = 300
  integer, parameter :: SnooperDiedTag = 301

  character (LEN=*), parameter :: Level2CodeGroupName="MLSL2Executable"

  ! Now the type definitions for snooping

  type SnooperInfo_T
    integer :: tid                      ! Task ID of snooper
    integer :: mode                     ! Mode of the snooper
  end type SnooperInfo_T

  logical, save :: SNOOPINGACTIVE = .false.

  public :: SNOOPINGACTIVE
  public :: SNOOP

contains ! ========  Public Procedures ==========================================

  ! ---------------------------------------  GetSnooperModeString  -----
  subroutine GetSnooperModeString ( MODE, STRING )
    integer, intent(in) :: MODE
    character (len=*), intent(out) :: STRING

    select case (mode)
    case (SnooperObserving)
      string='Observing'
    case (SnooperControling)
      string='Controling'
    case default
    end select
  end subroutine GetSnooperModeString

  ! ------------------------------------- SendVectorsListToSnooper -----
  subroutine SendVectorsListToSnooper ( SNOOPER, VectorDatabase, &
    & AnotherVectorDatabase )

    ! Arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER
    type (Vector_T), dimension(:), optional, pointer :: VectorDatabase
    type (Vector_T), dimension(:), optional, pointer :: AnotherVectorDatabase

    ! Local variables
    logical :: AnyMoreVectors, AnyVectors    ! Flags
    integer :: BUFFERID                      ! For PVM
    integer :: INFO                          ! Flag
    character (len=132) :: LINE              ! A line of text
    integer :: VECTOR, QUANTITY              ! Loop counters
    
    anyVectors = present(vectorDatabase)
    if ( anyVectors ) anyVectors = associated(vectorDatabase)
    if ( anyVectors ) anyVectors = size(vectorDatabase) > 0
    anyMoreVectors = present(anotherVectorDatabase)
    if ( anyMoreVectors ) anyMoreVectors = associated(anotherVectorDatabase)
    if ( anyMoreVectors ) anyMoreVectors = size(anotherVectorDatabase) > 0
    if ( anyVectors .or. anyMoreVectors) then
      call PVMFInitSend ( PvmDataDefault, bufferID )

      call PVMIDLPack ( "Vectors", info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing 'Vectors'" )

      call PVMIDLPack ( size(vectorDatabase), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing size(vectorDatabase)" )

      ! Send vector names and quantity names
      do vector = 1, size(vectorDatabase)
        call get_string ( vectorDatabase(vector)%name, line )
        call PVMIDLPack ( TRIM(line), info )
        if ( info /=0 ) call PVMErrorMessage ( info, "packing vector name:" &
          & // TRIM(line) )
        
        call PVMIDLPack ( size(vectorDatabase(vector)%quantities), info )
        if ( info /=0 ) call PVMErrorMessage ( info, "packing no quantities for " &
          & // TRIM(line) )
        
        do quantity = 1, size(vectorDatabase(vector)%quantities)
          call get_string &
            & ( vectorDatabase(vector)%quantities(quantity)%template%name, line )
          call PVMIDLPack ( TRIM(line), info )
          if ( info /=0 ) call PVMErrorMessage ( info, "packing quantity name:" &
            & // TRIM(line) )
        end do
      end do

      ! Send more vector names and quantity names
      do vector = 1, size(anotherVectorDatabase)
        call get_string ( anotherVectorDatabase(vector)%name, line )
        call PVMIDLPack ( TRIM(line), info )
        if ( info /=0 ) call PVMErrorMessage ( info, "packing vector name:" &
          & // TRIM(line) )
        
        call PVMIDLPack ( size(anotherVectorDatabase(vector)%quantities), info )
        if ( info /=0 ) call PVMErrorMessage ( info, "packing no quantities for " &
          & // TRIM(line) )
        
        do quantity = 1, size(anotherVectorDatabase(vector)%quantities)
          call get_string &
            & ( anotherVectorDatabase(vector)%quantities(quantity)%template%name, line )
          call PVMIDLPack ( TRIM(line), info )
          if ( info /=0 ) call PVMErrorMessage ( info, "packing quantity name:" &
            & // TRIM(line) )
        end do
      end do

      ! Now send this buffer
      call PVMFSend ( snooper%tid, SnoopTag, info )
      if (info /= 0) call PVMErrorMessage ( info, "sending vector information" )

    else                                ! No vectors to send
      call PVMIDLSend ( "No Vectors", snooper%tid, info, msgTag=SnoopTag )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending 'No Vectors'" )
    endif

  end subroutine SendVectorsListToSnooper

  ! ---------------------------------------- AddSnooperToDatabase ------
  integer function AddSnooperToDatabase ( database, item )
    type (SnooperInfo_T), dimension(:), pointer :: DATABASE
    type (SnooperInfo_T), intent(in) :: ITEM
    ! Local variables
    type (SnooperInfo_T), dimension(:), pointer :: TEMPDATABASE
    include "addItemToDatabase.f9h"
    AddSnooperToDatabase=newSize
  end function AddSnooperToDatabase

  ! ------------------------------------------- RegisterNewSnooper -----
  subroutine RegisterNewSnooper ( snoopers, snooperTid ) 
    type (SnooperInfo_T), dimension(:), pointer :: snoopers
    integer, intent(in) :: snooperTid

    ! Local variables
    integer :: INFO
    type (SnooperInfo_T) :: NEWSNOOPER
    integer :: NOSNOOPERS

    ! Executable code

    ! Setup the new snooper information
    newSnooper%tid = snooperTid
    newSnooper%mode = SnooperObserving
    noSnoopers = AddSnooperToDatabase ( snoopers, newSnooper )

    ! Ask to be notified of its death
    call PVMFNotify ( PVMTaskExit, SnooperDiedTag, 1, &
      (/ snooperTid /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, &
      & "calling PVMFNotify" )

    ! Now if this is the first, let's have it controling us
    if ( size(snoopers) == 1 ) then
      snoopers(1)%mode = SnooperControling
      call PVMIDLSend ( 'ControlAccept', snooperTid, info, msgTag=SnoopTag )
      if ( info /= 0 ) call PVMErrorMessage ( info, &
        & "sending control information" )
    endif
    
  end subroutine RegisterNewSnooper

  ! ------------------------------------------------ ForgetSnooper -----
  subroutine ForgetSnooper ( snoopers, snooper )
    type (SnooperInfo_T), dimension(:), pointer :: SNOOPERS
    integer, intent(in) :: SNOOPER

    ! Local variables
    type (SnooperInfo_T), dimension(:), pointer :: NEWSNOOPERS
    integer :: STATUS

    ! Executable code
    nullify ( newSnoopers )
    allocate ( newSnoopers(size(snoopers)-1), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'newSnoopers' )

    if ( size(newSnoopers) > 0 ) then
      newSnoopers(1:snooper-1) = snoopers(1:snooper-1)
      newSnoopers(snooper:) = snoopers(snooper+1:)
    endif
    deallocate ( snoopers, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate//'snoopers' )

    snoopers => newSnoopers

  end subroutine ForgetSnooper

  ! ------------------------------------------ SendStatusToSnooper -----
  subroutine SendStatusToSnooper ( SNOOPER, STATUS, LOCATION, COMMENT, &
    & CONTROLED )
    type (SnooperInfo_T), intent(INOUT) :: SNOOPER
    character (len=*), intent(in) :: STATUS
    character (len=*), intent(in) :: LOCATION
    character (len=*), intent(in) :: COMMENT
    logical, intent(in) :: CONTROLED

    ! Local variables
    integer :: BUFFERID                 ! ID for PVM
    integer :: INFO                     ! Flag from PVM

    ! Executable code
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMIDLPack ( status, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing status" )
    call PVMIDLPack ( location, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing location" )
    call PVMIDLPack ( comment, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing comment" )
    call PVMIDLPack ( controled, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing controled" )
    call PVMFSend ( snooper%tid, SnoopTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "sending status etc." )
  end subroutine SendStatusToSnooper
    
  ! ------------------------------------------------------  SNOOP  -----

  ! This is the main routine in the module.  It can be called from anywhere
  ! within the code, and takes optional arguments which the user can supply
  ! to pass to and from the IDL end of the snooper.

  subroutine Snoop ( KEY, VectorDatabase, AnotherVectorDatabase )

    ! Arguments
    integer, intent(in), optional :: KEY ! Tree node where snoop called
    type (Vector_T), dimension(:), pointer, optional :: VectorDatabase
    type (Vector_T), dimension(:), pointer, optional :: AnotherVectorDatabase
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: MATRIXDATABASE
    
    ! Local parameter
    integer, parameter :: DELAY=100000  ! For Usleep, no. microsecs

    ! External (C) function
    external :: Usleep

    ! Local variables, first the more exciting ones.
    integer, save :: MYTID=0            ! Local task ID under PVM
    type (SnooperInfo_T), dimension(:), pointer, save :: SNOOPERS => NULL()
    type (SnooperInfo_T), dimension(:), pointer, save :: OLDSNOOPERS => NULL()
    ! For add/del ops.

    ! Now the more mundane items
    integer :: BYTES
    integer :: BUFFERID, INFO           ! Flags and ids from PVM
    character (len=132) :: COMMENT      ! Comment field to snoop command
    integer :: CONTROLINGSNOOPER        ! This one is controling
    integer :: INUM                     ! Index in group
    character (len=132) :: LINE         ! Line of text received
    character (len=132) :: NEXTLINE     ! Line of text received
    integer :: MSGTAG                   ! Incomming Message tag
    character (len=132) :: LOCATION     ! Line of text to send off
    integer :: SNOOPER                  ! Loop counter
    integer :: SNOOPERTID               ! Task ID for snooper
    integer :: STATUS                   ! Status from allocate/deallocate
    logical :: KEEPWAITING                ! First time snoop called
    logical :: GOTSOMETHING             ! Set if we got a message

    ! Executable code
    if ( .not. snoopingActive ) return

    if ( present(key) ) then
      write ( location, * ) source_ref(key)/256
      location=adjustl(trim(location))
      call get_string ( sub_rosa(subtree(2,subtree(2,key))), comment, strip=.true. )
    else
      location = '????'
      comment = 'Unknown'
    endif

    ! If this is the very first call enroll in PVM for the first time, allocate
    ! 0 snoopers to start with.
    if ( myTid==0 ) then
      keepWaiting = .true.
      call PVMfmytid ( myTid )
      if ( myTid<=0 ) call PVMErrorMessage ( myTid, "Enroling in PVM" )
      call PVMfjoingroup ( Level2CodeGroupName, inum )
      if ( inum<0 ) call PVMErrorMessage ( inum, "Joining group " &
        & // Level2CodeGroupName )
      allocate ( snoopers(0), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & MLSMSG_Allocate // "snoopers(0)" )
    else
      keepWaiting = .false.
    end if

    ! Tell all the snoopers we're ready to talk to them
    do snooper = 1, size(snoopers)
      call SendStatusToSnooper ( snoopers(snooper), 'Ready', location, comment, &
        & any ( snoopers%mode == SnooperControling ) )
      call SendVectorsListToSnooper ( snoopers(snooper), vectorDatabase, &
        & anotherVectorDatabase )
    end do

    snoopEventLoop: do ! ---------------------------- Snoop event loop -----

      gotSomething = .false.
      ! Now try to receive a message
      call PVMFNRecv ( -1, SnoopTag, bufferID )
      if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for snoop message" )
      else if ( bufferID > 0 ) then
        ! We have something to look at
        ! Get first line and find out who sent it.
        gotSomething = .true.
        call PVMFBufInfo ( bufferID, bytes, msgTag, snooperTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        snooper = FindFirst ( snoopers%tid == snooperTid )
        call PVMIDLUnpack ( line, info )
        select case ( trim(line) )

        case ( 'NewSnooper' )
          keepWaiting = .false.
          call RegisterNewSnooper ( snoopers, snooperTid )
          snooper = FindFirst ( snoopers%tid == snooperTid )
          ! Tell it where we are
          call SendStatusToSnooper ( snoopers(snooper), 'Ready', location, comment, &
            & any(snoopers%mode == SnooperControling ) )
          ! Send it vectors
          call SendVectorsListToSnooper ( snoopers(snooper), vectorDatabase, &
            anotherVectorDatabase )
          
        case ( 'Finishing' )
          call ForgetSnooper ( snoopers, snooper )

        case ( 'Continue' )
          if ( snoopers(snooper)%mode == SnooperControling ) &
            exit SnoopEventLoop

        case ( 'Control' )
          if ( any ( snoopers%mode == SnooperControling ) ) then
            call PVMIDLSend ( 'ControlReject', snooperTid, info, msgTag=SnoopTag )
          else
            call PVMIDLSend ( 'ControlAccept', snooperTid, info, msgTag=SnoopTag )
            snoopers(snooper)%mode = SnooperControling
          end if
          if ( info /= 0 ) call PVMErrorMessage ( info, &
            & "sending control response" )

        case ( 'Relinquish' )
          snoopers(snooper)%mode = SnooperObserving
          call PVMIDLSend ( 'Relinquished', snooperTid, info, msgTag=SnoopTag )
          if ( info /= 0 ) call PVMErrorMessage ( info, &
            & "sending relinquish response" )

        case ( 'Quantity' )
          call PVMIDLUnpack ( nextLine, info )
          if ( info /=0 ) call PVMErrorMessage ( info, &
            & "unpacking response from snooper" )
          call SnooperRequestedQuantity ( snoopers(snooper), &
            & trim(nextLine), vectorDatabase )
          if ( present(anotherVectorDatabase) ) &
            & call SnooperRequestedQuantity ( snoopers(snooper), &
              & trim(nextLine), anotherVectorDatabase )

        case default
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Got unexpected response from snooper:'//trim(line) )

        end select
      end if

      ! Now check for dead task messages
      call PVMFNRecv ( -1, SnooperDiedTag, bufferID )
      if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for snooper died message" )
      else if ( bufferID > 0 ) then     ! Got a message
        gotSomething = .true.
        call PVMF90Unpack ( snooperTid, info )
        if ( info < 0 ) then
          call PVMErrorMessage ( info, "unpacking dead snooper tid" )
        end if
        snooper = FindFirst ( snoopers%tid == snooperTid )
        if ( snooper /= 0 ) then
          call ForgetSnooper ( snoopers, snooper )
        endif
      endif ! Got a message

      ! Shall we quit the loop?
      if ( size(snoopers) == 0 ) then
        if (.not. keepWaiting) exit SnoopEventLoop
      else
        if ( all ( snoopers%mode /= SnooperControling ) ) exit SnoopEventLoop
      end if

      if (.not. gotSomething) call usleep ( delay )
    end do snoopEventLoop ! -------------- End of snoop event loop -----

    ! Tell all the snoopers we're off and running again
    do snooper = 1, size(snoopers)
      call SendStatusToSnooper ( snoopers(snooper), 'Running', location, comment, &
        & any ( snoopers%mode == SnooperControling ) )
    end do
    
  end subroutine Snoop

  ! ------------------------------------- Snooper requested matrix -----
  subroutine SnooperRequestedMatrix ( SNOOPER, LINE, MATRIXDATABASE )
    ! This routine sends a matrix (including the vectors associated
    ! with it) to a snooping task.  It will probably take a while in many
    ! cases.

    ! Dummy arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER ! This snooper
    character (len=*), intent(in) :: LINE ! Name of matrix to send
    type (Matrix_T), dimension(:), pointer :: MATRIXDATABASE

    ! Local variables
    integer :: MATRIXNAME               ! String index of matrix name
    integer :: MATRIX                   ! Index of matrix in database

    ! Executable code

    matrixName = enter_terminal ( line, t_identifier )
    do matrix = 1, size(matrixDatabase)
      if ( matrixDatabase(matrix)%name == matrixName ) exit
    end do
    if ( matrixDatabase(matrix)%name /= matrixName ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to find requested matrix:'//trim(line) )

    !call PVMSendMatrix ( matrixDatabase(matrix), snooper%tid )

  end subroutine SnooperRequestedMatrix

  ! ----------------------------------- Snooper requested quantity -----
  subroutine SnooperRequestedQuantity ( SNOOPER, LINE, VectorDatabase )
    ! This routine sends a quantity (including its template) to 
    ! a snooping task.

    ! Dummy arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER ! This snooper
    character (len=*), intent(in) :: LINE ! Name of quantity to send
    type (Vector_T), dimension(:), pointer :: VectorDatabase

    ! Local variables
    integer :: BUFFERID                 ! ID for PVM
    integer :: INFO                     ! Flag from PVM
    integer :: DOTPOS                   ! For parsing string
    integer :: VECTORNAME               ! String index of vector name
    integer :: QUANTITYNAME             ! String index of quantity name
    integer :: VECTOR                   ! index
    integer :: QUANTITY                 ! index

    type (VectorValue_T), pointer :: q  ! This vector quantity
    type (QuantityTemplate_T), pointer :: qt ! The quantity template

    character (len=132) :: word         ! A line of text to send.

    ! Executable code

    dotPos = index ( line, '.' )
    vectorName = enter_terminal ( line(1:dotPos-1), t_identifier )
    quantityName = enter_terminal ( line(dotPos+1:), t_identifier )

    do vector = 1, size(vectorDatabase)
      if ( vectorDatabase(vector)%name == vectorName) exit
    end do
    if ( vector > size(vectorDatabase) ) return

    do quantity = 1, size(vectorDatabase(vector)%quantities)
      if ( vectorDatabase(vector)%quantities(quantity)%template%name == &
        & quantityName ) exit
    end do
    if ( quantity > size(vectorDatabase(vector)%quantities ) ) return

    q => vectorDatabase(vector)%quantities(quantity)

    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMIDLPack ( 'Quantity', info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'packing "Quantity"' )
    call PVMIDLPack ( trim(line), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'packing quantity name' )
    call PVMSendQuantity( q, snooper%tid, justPack=.true., noMask=.true. )
    call PVMFSend ( snooper%tid, SnoopTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, 'sending quantity' )

  end subroutine SnooperRequestedQuantity
  
end module SnoopMLSL2

! $Log$
! Revision 2.20  2001/10/05 17:32:35  vsnyder
! Add 'AnotherVectorDatabase' argument; cosmetic changes
!
! Revision 2.19  2001/10/02 00:44:46  livesey
! Bug fix
!
! Revision 2.18  2001/09/28 23:39:09  livesey
! Now doesn't give error if unknown quantity/vector requested
!
! Revision 2.17  2001/09/27 23:39:22  livesey
! Bug fix for going from one to zero snoopers
!
! Revision 2.16  2001/09/22 15:49:54  livesey
! Made delay in loop optional, only if nothing received
!
! Revision 2.15  2001/09/21 22:04:16  livesey
! Removed some print statements
!
! Revision 2.14  2001/09/20 23:08:21  livesey
! New version, minor tidy ups etc.
!
! Revision 2.13  2001/09/20 00:27:21  livesey
! Minor changes
!
! Revision 2.12  2001/09/19 23:47:39  livesey
! Compilable version
!
! Revision 2.11  2001/09/19 23:33:16  livesey
! New version of communication protocol
!
