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

  use VectorsModule, only: Vector_T
  !  USE MatrixModule_1, ONLY: Matrix_T

  use MLSCommon, only: R4,R8,I4
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
    MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate
  use PVM, only: PVMFbcast, PVMDataDefault, PVMFinitsend, PVMFmyTid, PVMFgSize
  use PVMIDL, only:  IDLMsgTag, PVMIDLPack, PVMIDLReceive, PVMIDLSend
  use TREE, only:  DUMP_TREE_NODE, SUB_ROSA, SUBTREE, SOURCE_REF
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use LEXER_CORE, only: PRINT_SOURCE

  !  USE VectorsModule, ONLY: 
  implicit none

  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=130), private :: Id = &
    "$Id$"
  character(LEN=*), parameter, private :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

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
  integer, parameter :: SnooperRequestingControling = SnooperControling + 1
  integer, parameter :: SnooperFinishing =            SnooperRequestingControling + 1

  character (LEN=*), parameter :: ReceptiveSnoopersGroupName="MLSL2ReceptiveSnoopers"
  character (LEN=*), parameter :: Level2CodeGroupName="MLSL2Executable"

  ! Now the type definitions for snooping

  type SnooperInfo_T
    integer :: tid             ! Task ID of snooper
    integer :: mode            ! Mode of the snooper
  end type SnooperInfo_T

contains ! ========  Public Procedures ==========================================

  ! -------------------------------------------------- PVMERRORMESSAGE ----------
  subroutine PVMErrorMessage(INFO,PLACE)
    ! This routine is called to log a PVM error
    integer, intent(IN) :: INFO
    character (LEN=*) :: PLACE

    character (LEN=132) :: LINE

    write (line, * ) info
    call MLSMessage(MLSMSG_Error,ModuleName,'PVM error '//trim(place)//&
      ' Info='//trim(adjustl(line)))
  end subroutine PVMErrorMessage

  ! -------------------------------------------------- LOOKFORSNOOPERS ----------
  subroutine LookForSnoopers(MYTID,SNOOPERS)
    ! This routine is called by snoop each time to keep an eye on snoopers. The
    ! ins and outs of the protocol are discussed in the routine itself

    ! Arguments
    integer, intent(IN) :: MYTID
    type (SnooperInfo_T), dimension(:), pointer :: SNOOPERS ! Info about each snooper

    ! Local paramaters
    integer, parameter :: SETUPMSGTAG = 50

    ! Local variables
    integer :: BUFFERID, INFO           ! PVM stuff
    character (LEN=132) :: MESSAGE      ! Messages that pass by
    character (LEN=132) :: FRAGMENT     ! Fragment of said messages
    integer :: GROUPSIZE                ! Number of receptive snoopers
    integer :: GROUPMEMBER              ! Loop counter
    integer :: STATUS                   ! From ALLOCATE/DEALLOCATE
    integer :: THISSNOOPERTID           ! TaskID
    integer :: NOSNOOPERS               ! Number of active snoopers

    type (SnooperInfo_T), dimension(:), pointer :: TEMPSNOOPERS

    ! Executable code

    ! First, we're going to see if there are any receptive snoopers
    call PVMFgsize(ReceptiveSnoopersGroupName,groupSize)

    ! If there are we're going to broadcast to them.
    if (groupSize>0) then

      call PVMFInitSend(PVMDataDefault,bufferID)
      write (message,*) myTid,"SolicitingSnoopers"
      call PVMIDLPack(message,info)
      if (info /= 0) call PVMErrorMessage(info,' packing SolicitingSnoopers buffer.')
      call PVMFbcast(ReceptiveSnoopersGroupName,SetupMsgTag,info)


      ! Now wait for a reply from each of them
      if (info == 0) then
        do groupMember=1,groupSize
          call PVMIDLReceive(message,-1,info,msgTag=SetupMsgTag)
          if (info /= 0) call PVMErrorMessage(info,&
            ' getting reply from listening snooper')
          read (message,*) thisSnooperTid,fragment
          if (trim(fragment)/="ReadyToSnoop") call MLSMessage(MLSMSG_Error, &
            ModuleName, &
            "Unexpected message from potential snooper: "//fragment)

          ! If really a new snooper add it.
          if (all(thisSnooperTid /= snoopers%tid)) then
            tempSnoopers=>snoopers
            noSnoopers=SIZE(snoopers)+1
            allocate(snoopers(noSnoopers),STAT=status)
            if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
              MLSMSG_Allocate//"snoopers in LookForSnoopers")
            if (size(tempSnoopers) /= 0) snoopers(1:noSnoopers-1) = tempSnoopers
            snoopers(noSnoopers)%tid = thisSnooperTid
            snoopers(noSnoopers)%mode = SnooperObserving
            if (size(snoopers)/=0) then
              deallocate(tempSnoopers, STAT=status)
              if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,&
                MLSMSG_Allocate//"tempSnoopers in LookForSnoopers")
            end if
          end if
        end do
      else
        call PVMErrorMessage(info,' broadcasting.')
      end if
    end if
  end subroutine LookForSnoopers

  ! ---------------------------------------------- DEALWITHONESNOOPER -----------

  subroutine DealWithOneSnooper(MYTID, LOCATION, COMMENT, &
    & SNOOPER, VECTORDATABASE,DONESNOOPINGFORNOW)

    ! Arguments
    integer, intent(IN) :: MYTID
    character (len=*), intent(in) :: LOCATION
    character (len=*), intent(in) :: COMMENT
    type (SnooperInfo_T), intent(INOUT) :: SNOOPER
    type (Vector_T), dimension(:), pointer, optional :: VECTORDATABASE
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: MATRIXDATABASE
    logical, intent(OUT), optional :: DONESNOOPINGFORNOW

    ! Local variables
    integer :: BUFFERID                 ! ID for PVM
    integer :: INFO                     ! Flag from PVM
    integer :: vector                   ! Loop counter
    integer :: quantity                 ! Loop counter
    character (len=132) :: line         ! Temporary string

    ! Executable code
    
    ! First we send a message to the snooper to indicate our presence and our
    ! location/comment.

    call PVMFInitSend(PvmDataDefault, bufferID)
    call PVMIDLPack(location, info)
    if (info /= 0) call PVMErrorMessage(info, "packing location")
    call PVMIDLPack(comment, info)
    if (info /= 0) call PVMErrorMessage(info, "packing comment")
    call PVMFSend(snooper%tid, IDLMsgTag, info)
    if (info /= 0) call PVMErrorMessage(info, "sending location comment")

    if (present(vectorDatabase)) then
      call PVMFInitSend(PvmDataDefault, bufferID)

      call PVMIDLPack("Vectors", info)
      if (info /= 0) call PVMErrorMessage(info, "packing 'Vectors'")

      call PVMIDLPack(size(vectorDatabase), info)
      if (info /= 0) call PVMErrorMessage(info, "packing size(vectorDatabase)")

      ! Send vector names and quantity names
      do vector=1,size(vectorDatabase)
        call get_string(vectorDatabase(vector)%name,line)
        call PVMIDLPack(TRIM(line), info)
        if (info /=0) call PVMErrorMessage(info, "packing vector name:" &
          & // TRIM(line))
        
        call PVMIDLPack(size(vectorDatabase(vector)%quantities), info)
        if (info /=0) call PVMErrorMessage(info, "packing no quantities for " &
          & // TRIM(line))
        
        do quantity=1,size(vectorDatabase(vector)%quantities)
          call get_string(vectorDatabase(vector)%quantities(quantity)%template%name,line)
          call PVMIDLPack(TRIM(line), info)
          if (info /=0) call PVMErrorMessage(info, "packing quantity name:" &
            & // TRIM(line))
        end do
      end do

      ! Now send this buffer
      call PVMFSend(snooper%tid, IDLMsgTag, info)
      if (info /= 0) call PVMErrorMessage(info, "sending vector information")

    else                                ! No vectors to send
      call PVMIDLSend("No Vectors", snooper%tid, info)
      if (info /= 0) call PVMErrorMessage(info, "sending 'No Vectors'")
    endif

    
    

  end subroutine DealWithOneSnooper

  ! ------------------------------------------------------ SNOOP ----------------

  ! This is the main routine in the module.  It can be called from anywhere
  ! within the code, and takes optional arguments which the user can supply to
  ! pass to and from the IDL end of the snooper. 

  subroutine Snoop(KEY, VECTORDATABASE)

    ! Arguments
    integer, intent(in), optional :: KEY ! Tree node where snoop called
    type (Vector_T), dimension(:), pointer, optional :: VECTORDATABASE
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: MATRIXDATABASE

    ! Local variables, first the more exciting ones.
    integer, save :: MYTID=0            ! Local task ID under PVM
    type (SnooperInfo_T), dimension(:), pointer, save :: SNOOPERS => null()
    type (SnooperInfo_T), dimension(:), pointer, save :: OLDSNOOPERS => null()
    ! For add/del ops.

    ! Now the more mundane items
    integer :: BUFFERID, INFO           ! Flags and ids from PVM
    integer :: INUM                     ! Index in group
    integer :: SNOOPER                  ! Loop counter
    integer :: CONTROLINGSNOOPER        ! This one is controling
    integer :: STATUS                   ! Status from allocate/deallocate
    logical :: DONESNOOPINGFORNOW       ! Flag to end loop

    character (len=132) :: COMMENT      ! Comment field to snoop command
    character (len=132) :: LOCATION     ! Line of text to send off
    ! Executable code

    if ( present(key) ) then
      write (location,*) source_ref(key)/256
      location='Line '//adjustl(trim(location))
      call get_string(sub_rosa(subtree(2,subtree(2,key))),comment,strip=.true.)
    else
      location='Unknown'
      comment='Unknown'
    endif

    ! If this is the very first call enroll in PVM for the first time, allocate
    ! 0 snoopers to start with.
    if ( myTid==0 ) then
      call pvmfmytid(myTid)
      if (myTid<=0) call PVMErrorMessage(myTid,"Enroling in pvm")
      call pvmfjoingroup(Level2CodeGroupName,inum)
      if (inum<0) call PVMErrorMessage(inum,"Joining group "&
        &//Level2CodeGroupName)
      allocate(snoopers(0), stat=status)
      if (status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName,&
        & MLSMSG_Allocate // "snoopers(0)")
    end if

    ! Now look for snoopers
    call LookForSnoopers(myTid,snoopers)

    ! Now if we have some snoopers talk to them, set loop to exit by default
    ! for the case where there are now controling snoopers.
    if (size(snoopers)/=0) then
      doneSnoopingForNow=.true.
      snoopLoop: do
        controlingSnooper=0

        ! Go through all snoopers, deal with any non-controling ones first
        do snooper=1,size(snoopers)
          if (snoopers(snooper)%mode==SnooperControling) then
            controlingSnooper=snooper
          else
            call DealWithOneSnooper(myTid, TRIM(location), TRIM(comment), &
              & snoopers(snooper), vectorDatabase)
          endif
        end do

        ! Deal with controling snooper if any
        if (controlingSnooper/=0) &
          call DealWithOneSnooper(myTid, TRIM(location), TRIM(comment), &
          & snoopers(controlingSnooper), vectorDatabase, &
          & doneSnoopingForNow=doneSnoopingForNow)

        ! Deal with any snoopers requesting controling
        if (any(snoopers%mode==SnooperRequestingControling)) then
          if (controlingSnooper/=0) snoopers(controlingSnooper)%mode=SnooperObserving
          where (snoopers%mode==SnooperRequestingControling) &
            snoopers%mode=SnooperControling
        end if

        ! Deal with any snoopers about to finish
        if (any(snoopers%mode==SnooperFinishing)) then 
          allocate(oldSnoopers(size(snoopers)),STAT=status)
          if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,&
            MLSMSG_Allocate//"oldSnoopers")
          oldSnoopers=snoopers
          allocate(snoopers(count(snoopers%mode/=SnooperFinishing)),STAT=status)
          if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,&
            MLSMSG_Allocate//"snoopers")
          snoopers=pack(oldSnoopers,snoopers%mode/=SnooperFinishing)
          deallocate(oldSnoopers,STAT=status)
          if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,&
            MLSMSG_DeAllocate//"oldSnoopers")
        endif

        ! This shouldn't be necessary, but just to be sure...
        if (size(snoopers)==0) doneSnoopingForNow= .true.

        if (doneSnoopingForNow) exit snoopLoop
      end do snoopLoop
    end if
  end subroutine Snoop

end module SnoopMLSL2
