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
  use PVM, only: PVMFbcast, PVMDataDefault, PVMFinitsend, PVMFmyTid, PVMFgSize, &
    & PVMErrorMessage
  use PVMIDL, only:  IDLMsgTag, PVMIDLPack, PVMIDLReceive, PVMIDLSend, PVMIDLUnpack
  use QuantityPVM, only: PVMSENDQUANTITY
  use TREE, only:  DUMP_TREE_NODE, SUB_ROSA, SUBTREE, SOURCE_REF
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use LEXER_CORE, only: PRINT_SOURCE
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_IDENTIFIER
  use MLSSignals_M, only: GETSIGNALNAME

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

  character (LEN=*), parameter :: ReceptiveSnoopersGroupName="MLSL2ReceptiveSnoopers"
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

  ! ---------------------------------------  GETSNOOPERMODESTRING  -----
  subroutine GetSnooperModeString ( MODE, STRING )
    integer, intent(in) :: MODE
    character (len=*), intent(out) :: STRING

    select case (mode)
    case (SnooperObserving)
      string='Observing'
    case (SnooperControling)
      string='Controling'
    case (SnooperFinishing)
      string='Finishing'
    case default
    end select
  end subroutine GetSnooperModeString

  ! ------------------------------------------- SENDVECTORSLISTTOSNOOPER --------
  subroutine SendVectorsListToSnooper ( SNOOPER, VECTORDATABASE )

    ! Arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER
    type (Vector_T), dimension(:), optional, pointer :: VECTORDATABASE

    ! Local variables
    logical :: ANYVECTORS               ! Flag
    integer :: BUFFERID                 ! For PVM
    integer :: INFO                     ! Flag
    character (len=132) :: LINE         ! A line of text
    integer :: VECTOR, QUANTITY         ! Loop counters
    
    anyVectors = present(vectorDatabase)
    if ( anyVectors ) anyVectors = associated(vectorDatabase)
    if ( anyVectors ) then
      call PVMFInitSend ( PvmDataDefault, bufferID )

      call PVMIDLPack ( "Vectors", info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing 'Vectors'" )

      call PVMIDLPack ( size(vectorDatabase), info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing size(vectorDatabase)" )

      ! Send vector names and quantity names
      do vector = 1,size(vectorDatabase)
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

      ! Now send this buffer
      call PVMFSend ( snooper%tid, IDLMsgTag, info )
      if (info /= 0) call PVMErrorMessage ( info, "sending vector information" )

    else                                ! No vectors to send
      call PVMIDLSend ( "No Vectors", snooper%tid, info )
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

  ! --------------------------------------- RegisterNewSnooper ---------
  subroutine RegisterNewSnooper ( snoopers, snooperTid ) 
    type (SnooperInfo_T), dimension(:), pointer :: snoopers
    integer, intent(in) :: snooperTid

    ! Local variables
    type (SnooperInfo_T) :: newSnooper
    integer :: noSnoopers

    ! Executable code

    ! Setup the new snooper information
    newSnooper%tid = snooperTid
    newSnooper%mode = SnooperObserving
    noSnoopers = AddSnooperToDatabase ( snooper, newSnooper )

    ! Now if this is the first, let's have it controling us
    if ( size(snoopers) == 1 ) then
      snoopers(1)%mode = SnooperControling
      call PVMIDLSend ( 'ControlAccept', snooperTid, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, &
        & "sending control information" )
    endif
    
  end subroutine RegisterNewSnooper

  ! ---------------------------------------- ForgetSnooper -------------
  subroutine ForgetSnooper ( snoopers, snooper )
    type (SnooperInfo_T), dimension(:), pointer :: SNOOPERS
    integer, intent(in) :: SNOOPER

    ! Local variables
    type (SnooperInfo_T), dimension(:), pointer :: NEWSNOOPERS
    integer :: STATUS

    ! Executable code
    nullify ( newSnoopers )
    allocate ( newSnoopers, size(snoopers)-1, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'newSnoopers' )

    newSnoopers(1:snooper-1) = snoopers(1:snooper-1)
    newSnoopers(snooper:) = snoopers(snooper+1:)
    deallocate ( snoopers, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_DeAllocate//'snoopers' )

    snoopers => newSnoopers

  end subroutine ForgetSnooper

  ! ---------------------------------------- SendStatusToSnooper -------
  subroutine SendStatusToSnooper ( SNOOPER, LOCATION, COMMENT, STATUS, CONTROLED )
    type (SnooperInfo_T), intent(INOUT) :: SNOOPER
    character (len=*), intent(in) :: LOCATION
    character (len=*), intent(in) :: COMMENT
    character (len=*), intent(in) :: STATUS
    logical, intent(in) :: CONTROLED

    ! Local variables
    integer :: BUFFERID                 ! ID for PVM
    integer :: INFO                     ! Flag from PVM

    ! Executable code
    call PVMFInitSend ( PvmDataDefault, bufferID )
    call PVMIDLPack ( location, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing location" )
    call PVMIDLPack ( comment, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing comment" )
    call PVMIDLPack ( status, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing status" )
    call PVMIDLPack ( controled, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing controled" )
    call PVMFSend ( snooper%tid, IDLMsgTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "sending status etc." )
  end subroutine SendStatusToSnooper
    
  ! ------------------------------------------------------  SNOOP  -----

  ! This is the main routine in the module.  It can be called from anywhere
  ! within the code, and takes optional arguments which the user can supply to
  ! pass to and from the IDL end of the snooper. 

  subroutine Snoop ( KEY, VECTORDATABASE )

    ! Arguments
    integer, intent(in), optional :: KEY ! Tree node where snoop called
    type (Vector_T), dimension(:), pointer, optional :: VECTORDATABASE
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
    integer :: BUFFERID, INFO           ! Flags and ids from PVM
    character (len=132) :: COMMENT      ! Comment field to snoop command
    integer :: CONTROLINGSNOOPER        ! This one is controling
    logical :: DONESNOOPINGFORNOW       ! Flag to end loop
    logical :: FIRSTTIME                ! First time round the polling loop
    integer :: INUM                     ! Index in group
    character (len=132) :: LOCATION     ! Line of text to send off
    integer :: SNOOPER                  ! Loop counter
    integer :: STATUS                   ! Status from allocate/deallocate

    ! Executable code
    if ( .not. snoopingActive ) return

    if ( present(key) ) then
      write ( location, * ) source_ref(key)/256
      location='Line '//adjustl(trim(location))
      call get_string ( sub_rosa(subtree(2,subtree(2,key))), comment, strip=.true. )
    else
      location = 'Unknown'
      comment = 'Unknown'
    endif

    ! If this is the very first call enroll in PVM for the first time, allocate
    ! 0 snoopers to start with.
    if ( myTid==0 ) then
      call PVMfmytid ( myTid )
      if ( myTid<=0 ) call PVMErrorMessage ( myTid, "Enroling in PVM" )
      call PVMfjoingroup ( Level2CodeGroupName, inum )
      if ( inum<0 ) call PVMErrorMessage ( inum, "Joining group " &
        & // Level2CodeGroupName )
      allocate ( snoopers(0), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & MLSMSG_Allocate // "snoopers(0)" )
    end if

    ! Tell all the snoopers we're ready to talk to them
    do snooper = 1, size(snoopers)
      call SendStatusToSnooper ( snoopers(snooper), location, comment, &
        & 'Ready', any ( snoopers%mode == SnooperControling ) )
      call SendVectorsListToSnooper ( snoopers(snooper), vectorDatabase )
    end do

    doneSnoopingForNow = size(snoopers) /= 0 .and. &
      & all ( snoopers%mode /= SnooperControling )

    snoopEventLoop: do ! --------------------------- Snoop event loop ------

      ! Now try to receive a message
      call PVMFNRecv ( -1, SnoopTag, bufferID )
      if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for snoop message" )
      else if ( bufferID > 0 ) then
        ! We have something to look at
        ! Get first line and find out who sent it.
        call PVMFBufInfo ( bufferID, bytes, msgTag, snooperTid, info )
        if ( info /= 0 ) &
          & call PVMErrorMessage ( info, "calling PVMFBufInfo" )
        snooper = FindFirst ( snoopers%tid, snooperTid )
        call PVMIDLUnpack ( line, info )

        select case ( trim(line) )

        case ( 'NewSnooper' )
          call RegisterNewSnooper ( snoopers, snooperTid )
          
        case ( 'Finishing' )
          call ForgetSnooper ( snoopers, snooper )

        case ( 'Continue' )
          if ( snoopers(snooper)%mode == SnooperControling ) &
            & doneSnoopingForNow = .true.

        case ( 'Control' )
          if ( any ( snoopers%mode == SnooperControling ) ) then
            call PVMIDLSend ( 'ControlReject', snooperTid, info )
          else
            call PVMIDLSend ( 'ControlAccept', snooperTid, info )
            snoopers(snooper)%mode = SnooperControling
          end if
          if ( info /= 0 ) call PVMErrorMessage ( info, &
            & "sending control response" )

        case ( 'Relinquish' )
          snoopers(snooper)%mode = SnooperObserving

        case ( 'Quantity' )
          call PVMIDLUnpack ( nextLine, info )
          if ( info /=0 ) call PVMErrorMessage ( info, &
            & "unpacking response from snooper" )
          call SnooperRequestedQuantity(snoopers(snooper), &
            & trim(nextLine), vectorDatabase)

        case default
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Got unexpected response from snooper:'//line )

        end select
      end if

      ! Now check for dead task messages
      call PVMFNRecv ( -1, SnooperDiedTag, bufferID )
      if ( bufferID < 0 ) then
        call PVMErrorMessage ( info, "checking for snooper died message" )
      else if ( bufferID > 0 ) then     ! Got a message
        call PVMF90Unpack ( snooperTid, info )
        if ( info < 0 ) then
          call PVMErrorMessage ( info, "unpacking dead snooper tid" )
        end if
        snooper = FindFirst ( snoopers%tid, snooperTid )
        if ( snooper /= 0 ) then
          call ForgetSnooper ( snoopers, snooper )
        endif
      endif ! Got a message

      ! End of the loop
      if ( doneSnoopingForNow ) exit snoopEventLoop
      call usleep ( delay )
    end do snoopEventLoop ! ------------------ End of snoop event loop -----

    ! Tell all the snoopers we're off and running again
    do snooper = 1, size(snoopers)
      call SendStatusToSnooper ( snoopers(snooper), location, comment, &
        &  'Running', any ( snoopers%mode == SnooperControling ) )
    end do
    
  end subroutine Snoop

  ! ----------------------------------------- Snooper requested matrix -----
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

  ! ----------------------------------------- Snooper requested quantity ---
  subroutine SnooperRequestedQuantity ( SNOOPER, LINE, VECTORDATABASE )
    ! This routine sends a quantity (including its template) to 
    ! a snooping task.

    ! Dummy arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER ! This snooper
    character (len=*), intent(in) :: LINE ! Name of quantity to send
    type (Vector_T), dimension(:), pointer :: VECTORDATABASE

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
    if ( vectorDatabase(vector)%name /= vectorName ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to find requested vector: '//trim(line) )

    do quantity = 1, size(vectorDatabase(vector)%quantities)
      if ( vectorDatabase(vector)%quantities(quantity)%template%name == &
        & quantityName ) exit
    end do
    if ( vectorDatabase(vector)%quantities(quantity)%template%name /= &
      & quantityName ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to find requested vector quantity: '//line )

    q => vectorDatabase(vector)%quantities(quantity)

    call PVMSendQuantity( q, snooper%tid )

  end subroutine SnooperRequestedQuantity
  
end module SnoopMLSL2

! $Log$
! Revision 2.11  2001/09/19 23:33:16  livesey
! New version of communication protocol
!
