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
  !  USE MatrixModule_1, ONLY: Matrix_T

  use Intrinsic, only: LIT_INDICES
  use MLSCommon, only: R4, R8, I4
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
    MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate
  use PVM, only: PVMFbcast, PVMDataDefault, PVMFinitsend, PVMFmyTid, PVMFgSize
  use PVMIDL, only:  IDLMsgTag, PVMIDLPack, PVMIDLReceive, PVMIDLSend, PVMIDLUnpack
  use TREE, only:  DUMP_TREE_NODE, SUB_ROSA, SUBTREE, SOURCE_REF
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use LEXER_CORE, only: PRINT_SOURCE
  use Symbol_Table, only: ENTER_TERMINAL
  use Symbol_Types, only: T_IDENTIFIER
  use MLSSignals_M, only: GETSIGNALNAME

  implicit none

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
  integer, parameter :: SnooperFinishing =            SnooperControling + 1

  character (LEN=*), parameter :: ReceptiveSnoopersGroupName="MLSL2ReceptiveSnoopers"
  character (LEN=*), parameter :: Level2CodeGroupName="MLSL2Executable"

  ! Now the type definitions for snooping

  type SnooperInfo_T
    integer :: tid                      ! Task ID of snooper
    integer :: mode                     ! Mode of the snooper
    logical :: new             ! Used as flag
  end type SnooperInfo_T

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

  ! --------------------------------------------  PVMERRORMESSAGE  -----
  subroutine PVMErrorMessage ( INFO, PLACE )
    ! This routine is called to log a PVM error
    integer, intent(IN) :: INFO
    character (LEN=*) :: PLACE

    character (LEN=132) :: LINE

    write (line, * ) info
    call MLSMessage(MLSMSG_Error,ModuleName,'PVM error '//trim(place)//&
      ' Info='//trim(adjustl(line)))
  end subroutine PVMErrorMessage

  ! --------------------------------------------  LOOKFORSNOOPERS  -----
  subroutine LookForSnoopers ( MYTID, SNOOPERS )
    ! This routine is called by snoop each time to keep an eye on snoopers. The
    ! ins and outs of the protocol are discussed in the routine itself

    ! Arguments
    integer, intent(IN) :: MYTID
    type (SnooperInfo_T), dimension(:), pointer :: SNOOPERS ! Info about each snooper

    ! Local paramaters
    integer, parameter :: SETUPMSGTAG = 50

    ! Local variables
    integer :: BUFFERID, INFO           ! PVM stuff
    character (LEN=132) :: FRAGMENT     ! Fragment of said messages
    integer :: GROUPMEMBER              ! Loop counter
    integer :: GROUPSIZE                ! Number of receptive snoopers
    character (LEN=132) :: MESSAGE      ! Messages that pass by
    integer :: NOSNOOPERS               ! Number of active snoopers
    integer :: STATUS                   ! From ALLOCATE/DEALLOCATE
    integer :: THISSNOOPERTID           ! TaskID

    type (SnooperInfo_T), dimension(:), pointer :: TEMPSNOOPERS

    ! Executable code

    ! First, we're going to see if there are any receptive snoopers
    call PVMFgsize(ReceptiveSnoopersGroupName,groupSize)

    ! If there are we're going to broadcast to them.
    if ( groupSize > 0 ) then

      call PVMFInitSend ( PVMDataDefault, bufferID )
      write ( message,* ) myTid, "SolicitingSnoopers"
      call PVMIDLPack ( message, info )
      if ( info /= 0 ) &
        & call PVMErrorMessage ( info, ' packing SolicitingSnoopers buffer.')
      call PVMFbcast ( ReceptiveSnoopersGroupName, SetupMsgTag, info )


      ! Now wait for a reply from each of them
      if ( info == 0 ) then
        do groupMember = 1,groupSize
          call PVMIDLReceive ( message, -1, info, msgTag=SetupMsgTag )
          if (info /= 0) call PVMErrorMessage ( info, &
            ' getting reply from listening snooper' )
          read ( message,* ) thisSnooperTid, fragment
          if ( trim(fragment)/="ReadyToSnoop" ) &
            call MLSMessage ( MLSMSG_Error, ModuleName, &
            "Unexpected message from potential snooper: "//fragment )

          ! If really a new snooper add it.
          if ( all(thisSnooperTid /= snoopers%tid) ) then
            tempSnoopers => snoopers
            noSnoopers = SIZE(snoopers)+1
            allocate ( snoopers(noSnoopers), STAT=status )
            if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,&
              MLSMSG_Allocate//"snoopers in LookForSnoopers" )
            if ( size(tempSnoopers) /= 0 ) snoopers(1:noSnoopers-1) = tempSnoopers
            snoopers(noSnoopers)%tid = thisSnooperTid
            if ( noSnoopers == 1 ) then
              snoopers(noSnoopers)%mode = SnooperControling
            else
              snoopers(noSnoopers)%mode = SnooperObserving
            endif
            snoopers(noSnoopers)%new = .true.
            if (size(snoopers)/=0) then
              deallocate ( tempSnoopers, STAT=status )
              if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,&
                MLSMSG_Allocate//"tempSnoopers in LookForSnoopers")
            end if

            ! Now send the snooper a message indicating its mode
            call GetSnooperModeString ( snoopers(noSnoopers)%mode, message )
            call PVMIDLSend ( message, snoopers(noSnoopers)%tid, info, SetupMsgTag )
            if ( info /= 0 ) call PVMErrorMessage ( info, &
              & 'sending snooper mode' )
          end if
        end do
      else
        call PVMErrorMessage(info,' broadcasting.')
      end if
    end if
  end subroutine LookForSnoopers

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

  ! -----------------------------------------  DEALWITHONESNOOPER  -----

  subroutine DealWithOneSnooper ( MYTID, FIRSTTIME, NONECONTROLING, LOCATION, &
    & COMMENT, SNOOPER, VECTORDATABASE, DONESNOOPINGFORNOW )

    ! Arguments
    integer, intent(IN) :: MYTID
    logical, intent(IN) :: FIRSTTIME    ! First time for this snooper at this point
    logical, intent(IN) :: NONECONTROLING ! So we can take control if we want
    character (len=*), intent(in) :: LOCATION
    character (len=*), intent(in) :: COMMENT
    type (SnooperInfo_T), intent(INOUT) :: SNOOPER
    type (Vector_T), dimension(:), pointer, optional :: VECTORDATABASE
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: MATRIXDATABASE
    logical, intent(INOUT), optional :: DONESNOOPINGFORNOW

    ! Local variables
    integer :: BUFFERID                 ! ID for PVM
    integer :: INFO                     ! Flag from PVM
    character (len=132) :: Line         ! Temporary string
    character (len=132) :: NextLine     ! Temporary string
    logical :: MyDone                   ! I'm done

    ! Executable code
    
    ! First we send a message to the snooper to indicate our presence and our
    ! location/comment.

    if ( firstTime .or. snooper%new ) then
      snooper%new = .false.
      call PVMFInitSend ( PvmDataDefault, bufferID )
      call PVMIDLPack ( location, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing location" )
      call PVMIDLPack ( comment, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing comment" )
      call PVMFSend ( snooper%tid, IDLMsgTag, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "sending location comment" )

      ! Now we send our list of current vectors
      call SendVectorsListToSnooper ( snooper, vectorDatabase )
      
      ! Now we send our list of current matrices
      ! call SendMatricesListToSnooper(snooper, matrixDatabase)
    endif

    ! Now see if there are some instructions from the snooper
    call PVMFNRecv ( snooper%tid, IDLMsgTag, bufferID )
    if ( bufferID > 0 ) then
      ! Something to hear
      call PVMIDLUnpack ( line, info )
      if ( info /=0 ) call PVMErrorMessage ( info, &
        & "unpacking response from snooper" )
      
      select case ( trim(line) )
      case ( 'Continue' )                ! Finished, program can continue
        if ( present(doneSnoopingForNow) ) doneSnoopingForNow=.true.
        myDone = .true.
      case ( 'Control' )
        if ( noneControling ) then
          call PVMIDLSend ( 'OK', snooper%tid, info )
          snooper%mode = SnooperControling
        else
          call PVMIDLSend ( 'NO', snooper%tid, info )
        endif
        if ( info /= 0 ) call PVMErrorMessage ( info, "sending control response" )
      case ( 'Relinquish' )
        snooper%mode = SnooperObserving
      case ( 'Finishing' )
        snooper%mode = SnooperFinishing
      case ( 'Quantity' )
        call PVMIDLUnpack ( nextLine, info )
        if ( info /=0 ) call PVMErrorMessage ( info, &
          & "unpacking response from snooper" )
        call SnooperRequestedQuantity(snooper, trim(nextLine), vectorDatabase)
      case default
      end select

    endif                               ! Something to hear

  end subroutine DealWithOneSnooper

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
    integer, parameter :: DELAY=200000  ! For Usleep, no. microsecs

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

    ! Now look for snoopers
    call LookForSnoopers ( myTid, snoopers )

    ! Now if we have some snoopers talk to them, set loop to exit by default
    ! for the case where there are now controling snoopers.
    if ( size(snoopers) /= 0 ) then
      doneSnoopingForNow = .false.
      firstTime = .true.
      snoopLoop: do
        controlingSnooper=0

        ! Go through all snoopers, deal with any non-controling ones first
        do snooper = 1, size(snoopers)
          if ( snoopers(snooper)%mode==SnooperControling ) then
            controlingSnooper = snooper
          else
            call DealWithOneSnooper ( myTid, firstTime, &
              & all(snoopers%mode /= SnooperControling), &
              & trim(location), trim(comment), &
              & snoopers(snooper), vectorDatabase )
          end if
        end do

        ! Deal with controling snooper if any
        if ( controlingSnooper /= 0 ) then
           call DealWithOneSnooper ( myTid, firstTime, .false., &
                & TRIM(location), TRIM(comment), &
                & snoopers(controlingSnooper), vectorDatabase, &
                & doneSnoopingForNow=doneSnoopingForNow )
        end if

        ! Deal with any snoopers about to finish
        if ( any(snoopers%mode==SnooperFinishing) ) then 
          allocate ( oldSnoopers(size(snoopers)), STAT=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_Allocate//"oldSnoopers" )
          oldSnoopers = snoopers
          allocate ( snoopers(count(snoopers%mode/=SnooperFinishing)), STAT=status )
          if ( status /= 0 ) call MLSMessage(MLSMSG_Error,ModuleName, &
            MLSMSG_Allocate // "snoopers")
          snoopers=pack(oldSnoopers,oldSnoopers%mode/=SnooperFinishing)
          deallocate(oldSnoopers,STAT=status)
          if (status/=0) call MLSMessage ( MLSMSG_Error, ModuleName, &
            MLSMSG_DeAllocate // "oldSnoopers" )
        endif

        ! This shouldn't be necessary, but just to be sure...
        if ( (size(snoopers)==0) .or. (all(snoopers%mode/=SnooperControling)) ) &
          & doneSnoopingForNow = .true.
        firstTime = .false.
        if ( doneSnoopingForNow ) exit snoopLoop

        ! Now look for and new snoopers
        call LookForSnoopers ( myTid, snoopers )
        call usleep ( delay )
      end do snoopLoop
    end if
  end subroutine Snoop

  ! ----------------------------------------- Snooper requested quantity ---
  subroutine SnooperRequestedQuantity ( SNOOPER, LINE, VECTORDATABASE )
    ! This routine sends a quantity (including its template) to 
    ! a snooping task.

    ! Dummy arguments
    type (SnooperInfo_T), intent(in) :: SNOOPER ! This snooper
    character (len=*), intent(in) :: LINE ! Name of quantity to watch
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
      & 'Unable to find requested vector: '//line )

    do quantity = 1, size(vectorDatabase(vector)%quantities)
      if ( vectorDatabase(vector)%quantities(quantity)%template%name == &
        & quantityName ) exit
    end do
    if ( vectorDatabase(vector)%quantities(quantity)%template%name /= &
      & quantityName ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Unable to find requested vector quantity: '//line )

    q => vectorDatabase(vector)%quantities(quantity)
    qt => q%template

    ! Now we simply pack the quantity up and send it down the pvm spigot
    call PVMFInitSend ( PvmDataDefault, bufferID )

    call PVMIDLPack ( (/ qt%noInstances, qt%noSurfs, qt%noChans /), &
      & info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity dimensions." )

    call PVMIDLPack ( (/ qt%coherent, qt%stacked, qt%regular, qt%minorFrame, &
      & qt%logBasis /), info ) 
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantity flags" )

    call PVMIDLPack ( (/ qt%noInstancesLowerOverlap, &
      & qt%noInstancesUpperOverlap, qt%sideband /), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing misc quantity stuff" )

    ! Now pack some strings

    call Get_String( lit_indices(qt%quantityType), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing quantityType" )

    call Get_String( lit_indices(qt%verticalCoordinate), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing verticalCoordinate" )

    call Get_String( lit_indices(qt%unit), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing unit" )

    call Get_String( lit_indices(qt%frequencyCoordinate), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing frequencyCoordinate" )

    call GetSignalName( lit_indices(qt%signal), word, sideband=qt%sideband )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing signal" )

    call Get_String( lit_indices(qt%instrumentModule), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing instrumentModule" )

    call Get_String( lit_indices(qt%radiometer), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing radiometer" )

    call Get_String( lit_indices(qt%molecule), word, noError=.true. )
    call PVMIDLPack ( trim(word), info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing molecule" )

    ! Now pack the arrays

    call PVMIDLPack ( qt%surfs, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing surfs" )

    call PVMIDLPack ( qt%phi, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing phi" )

    call PVMIDLPack ( qt%geodLat, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing geodLat" )

    call PVMIDLPack ( qt%lon, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing lon" )

    call PVMIDLPack ( qt%time, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing time" )

    call PVMIDLPack ( qt%solarTime, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing solarTime" )

    call PVMIDLPack ( qt%solarZenith, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing solarZenith" )

    call PVMIDLPack ( qt%losAngle, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing losAngle" )

    if ( associated ( qt%mafIndex ) ) then
      call PVMIDLPack ( qt%mafIndex, info )
    else
      call PVMIDLPack ( (/ 0 /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing mafIndex" )

    if ( associated ( qt%mafCounter ) ) then
      call PVMIDLPack ( qt%mafCounter, info )
    else
      call PVMIDLPack ( (/ 0 /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing mafCounter" )

    if ( associated ( qt%frequencies ) ) then
      call PVMIDLPack ( qt%frequencies, info )
    else
      call PVMIDLPack ( (/ 0.0_r8 /), info )
    end if
    if ( info /= 0 ) call PVMErrorMessage ( info, "packing frequencies" )

    ! Finally the two arrays for irregular quantities

    if ( .not. qt%regular ) then
      call PVMIDLPack ( qt%surfIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing surfIndex" )

      call PVMIDLPack ( qt%chanIndex, info )
      if ( info /= 0 ) call PVMErrorMessage ( info, "packing chanIndex" )
    end if

    ! Now we're going to send this to the snooper.

    call PVMFSend ( snooper%tid, IDLMSGTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector template" )

    ! Now we're going to send the values in a separate message
    call PVMFInitSend ( PVMDataDefault, bufferID)

    ! Pack the values
    call PVMIDLPack ( q%values, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "sending values" )

    ! Skip the mask for the moment.

    ! Send this buffer
    call PVMFSend ( snooper%tid, IDLMSGTag, info )
    if ( info /= 0 ) call PVMErrorMessage ( info, "sending vector values" )
      
  end subroutine SnooperRequestedQuantity
  
end module SnoopMLSL2
