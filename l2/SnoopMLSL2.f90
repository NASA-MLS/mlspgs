! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

MODULE SnoopMLSL2               ! Interface between MLSL2 and IDL snooper via pvm3.

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

  USE VectorsModule, ONLY: Vector_T
  !  USE MatrixModule_1, ONLY: Matrix_T
  USE L2GPData, ONLY: L2GPData_T
  USE L2AUXData, ONLY: L2AUXData_T

  USE MLSCommon, ONLY: R4,R8,I4
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Warning, &
       MLSMSG_Info, MLSMSG_Allocate, MLSMSG_DeAllocate
  USE PVMIDL

  !  USE VectorsModule, ONLY: 

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130), PRIVATE :: Id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! =============================================================================

  ! The first main thing is a data type that describes each of the active
  ! snoopers.

  ! Part of this the following enumerated type, describing the `mode' of the
  ! snooper.  The mode can be `disinterested', `observing', or `controling'.  In 
  ! disinterested mode, the snooper and the l2 code will exchange no
  ! communication except a request to change to observing or controling on behalf of
  ! the snooper.  In observing mode, the snooper expects updates of status and
  ! possibly other parameters each invocation of Snoop.  In controling mode, the
  ! l2 code is expected to wait for instructions from the snooper before
  ! continuing.  Controling snoopers also have the right to change the contents
  ! of vectors and matrices etc. Only one snooper may be controling at a time.

  ! Note that snoopers communicate amongst themselves to arbitrate requests to
  ! become controling, the succuesful one contacts the l2 code.

  INTEGER, PARAMETER :: SnooperDisinterested=0
  INTEGER, PARAMETER :: SnooperObserving=1
  INTEGER, PARAMETER :: SnooperControling=2
  INTEGER, PARAMETER :: SnooperRequestingControling=3
  INTEGER, PARAMETER :: SnooperFinishing=4

  CHARACTER (LEN=*), PARAMETER :: ReceptiveSnoopersGroupName="MLSL2ReceptiveSnoopers"
  CHARACTER (LEN=*), PARAMETER :: Level2CodeGroupName="MLSL2Executable"

  ! Now the type definitions for snooping

  TYPE SnooperInfo_T
     INTEGER :: tid             ! Task ID of snooper
     INTEGER :: mode            ! Mode of the snooper
  END TYPE SnooperInfo_T

CONTAINS ! ========  Public Procedures ==========================================

  ! This routine is called to log a PVM error

  SUBROUTINE PVMErrorMessage(info,place)
    INTEGER, INTENT(IN) :: info
    CHARACTER (LEN=*) :: place

    CHARACTER (LEN=132) :: line

    WRITE (line,'(I)') info
    CALL MLSMessage(MLSMSG_Error,ModuleName,'PVM error '//TRIM(place)//&
         ' Info='//TRIM(ADJUSTL(line)))
  END SUBROUTINE PVMErrorMessage

  ! This routine is called by snoop each time to keep an eye on snoopers. The
  ! ins and outs of the protocol are discussed in the routine itself

  SUBROUTINE LookForSnoopers(myTid,snoopers,noSnoopers)

    ! Arguments
    INTEGER, INTENT(IN) :: myTid
    INTEGER, INTENT(INOUT) :: noSnoopers ! Number of snoopers we know about
    TYPE (SnooperInfo_T), DIMENSION(:), POINTER :: snoopers ! Info about each

    ! Local variables
    INTEGER :: bufferID, info   ! PVM stuff
    CHARACTER (LEN=132) :: message ! Messages that pass by
    CHARACTER (LEN=132) :: fragment ! Fragment of said messages
    INTEGER :: groupSize        ! Number of receptive snoopers
    INTEGER :: groupMember      ! Loop counter
    INTEGER :: status           ! From ALLOCATE/DEALLOCATE

    TYPE (SnooperInfo_T), DIMENSION(:), POINTER :: tempSnoopers

    ! Executable code

    ! First, we're going to see if there are any receptive snoopers
    CALL PVMFgsize(ReceptiveSnoopersGroupName,groupSize)

    ! If there are we're going to broadcast to them.
    IF (groupSize>0) THEN

       CALL PVMFInitSend(PVMDataDefault,bufferID)
       WRITE (message,'(I,A)') myTid,"SolicitingSnoopers"
       CALL PVMIDLPack(message,info)
       IF (info==0) CALL PVMErrorMessage(info,' packing SolicitingSnoopers buffer.')
       CALL PVMFbcast(ReceptiveSnoopersGroupName,IDLMsgTag,info)

       ! Now wait for a reply from each of them
       IF (info==0) THEN
          DO groupMember=1,groupSize
             CALL PVMIDLReceive(message,-1,info)
             IF (info/=0) CALL PVMErrorMessage(info,&
                  ' getting reply from listening snooper')
             READ (message,'(I,A)') thisSnooperTid,fragment
             IF (TRIM(fragment)/="ReadyToSnoop") CALL MLSMessage(MLSMSG_Error, &
                  ModuleName, &
                  "Unexpected message from potential snooper: "//fragment)

             ! If really a new snooper add it.
             IF (ALL(thisSnooperTid /= snoopers%tid)) THEN
                noSnoopers=noSnoopers+1
                IF (noSnoopers/=1) tempSnoopers=>snoopers
                ALLOCATE(snoopers(noSnoopers),STAT=status)
                IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
                     MLSMSG_Allocate//"snoopers in LookForSnoopers")
                IF (noSnoopers/=1) snoopers(1:noSnoopers-1)=tempSnoopers
                snoopers(noSnoopers)%tid=thisSnooperTid
                snoopers(noSnoopers)%mode=SnooperDisinterested
                IF (noSnoopers/=1) THEN
                   DEALLOCATE(tempSnoopers,STAT=status)
                   IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
                        MLSMSG_Allocate//"tempSnoopers in LookForSnoopers")
                END IF
             END IF
          END DO
       END IF
    END IF
  END SUBROUTINE LookForSnoopers

  ! -----------------------------------------------------------------------------

  SUBROUTINE DealWithOneSnooper(myTid,snooper, vectorDatabase,l2gpDatabase,&
       l2auxDatabase,doneSnoopingForNow)

    ! Arguments
    INTEGER, INTENT(IN) :: myTid
    TYPE (SnooperInfo_T), INTENT(INOUT) :: snooper
    TYPE (Vector_T), DIMENSION(:), POINTER, OPTIONAL :: vectorDatabase
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: matrixDatabase
    TYPE (L2GPData_T), DIMENSION(:), POINTER, OPTIONAL :: l2gpDatabase
    TYPE (L2AUXData_T), DIMENSION(:), POINTER, OPTIONAL :: l2auxDatabase
    LOGICAL, INTENT(OUT), OPTIONAL :: doneSnoopingForNow

    ! Local variables

    ! Executable code
  END SUBROUTINE DealWithOneSnooper

  ! -----------------------------------------------------------------------------

  ! This is the main routine in the module.  It can be called from anywhere
  ! within the code, and takes optional arguments which the user can supply to
  ! pass to and from the IDL end of the snooper. 

  SUBROUTINE Snoop(vectorDatabase, l2gpdatabase, l2auxDatabase)

    ! Arguments
    TYPE (Vector_T), DIMENSION(:), POINTER, OPTIONAL :: vectorDatabase
    !  TYPE (Matrix_T), DIMENSION(:), POINTER, OPTIONAL :: matrixDatabase
    TYPE (L2GPData_T), DIMENSION(:), POINTER, OPTIONAL :: l2gpDatabase
    TYPE (L2AUXData_T), DIMENSION(:), POINTER, OPTIONAL :: l2auxDatabase

    ! Local variables, first the more exciting ones.
    INTEGER, SAVE :: myTID=0                 ! Local task ID under PVM
    INTEGER, SAVE :: noSnoopers=0            ! Number of snoopers
    TYPE (SnooperInfo_T), DIMENSION(:), POINTER, SAVE :: snoopers
    TYPE (SnooperInfo_T), DIMENSION(:), POINTER, SAVE :: oldSnoopers 
    ! For add/del ops.

    ! Now the more mundane items
    INTEGER (i4) :: bufferID, info           ! Flags and ids from PVM
    INTEGER :: snooper                       ! Loop counter
    INTEGER :: controlingSnooper               ! This one is controling
    INTEGER :: status                        ! Status from allocate/deallocate
    LOGICAL :: doneSnoopingForNow            ! Flag to end loop

    ! Executable code

    ! If this is the very first call enroll in PVM for the first time.
    IF ( myTid==0 ) THEN
       CALL pvmfmytid(myTid)
       IF (myTid<=0) CALL PVMErrorMessage(myTid,"Enroling in pvm")
       CALL pvmfjoingroup(Level2CodeGroupName,inum)
       IF (inum<0) CALL PVMErrorMessage(inum,"Joining group "//Level2CodeGroupName)
    END IF

    ! Now look for snoopers
    CALL LookForSnoopers(myTid,snoopers,noSnoopers)

    ! Now if we have some snoopers talk to them, set loop to exit by default
    ! for the case where there are now controling snoopers.
    IF (noSnoopers/=0) THEN
       doneSnoopingForNow=.TRUE.
       snoopLoop: DO
          controlingSnooper=0
          
          ! Go through all snoopers, deal with any non-controling ones first
          DO snooper=1,noSnoopers
             IF (snoopers(snooper)%mode==SnooperControling) THEN
                controlingSnooper=snooper
             ELSE
                CALL DealWithOneSnooper(myTid,snoopers(snooper),vectorDatabase,&
                     l2gpDatabase,l2auxDatabase)
             ENDIF
          END DO
          
          ! Deal with controling snooper if any
          IF (controlingSnooper/=0) &
               CALL DealWithOneSnooper(myTid ,snoopers(controlingSnooper),&
               vectorDatabase, l2gpDatabase, l2auxDatabase, &
               doneSnoopingForNow=doneSnoopingForNow)
          
          ! Deal with any snoopers requesting controling
          IF (ANY(snoopers%mode==RequestingControling)) THEN
             IF (controlingSnooper/=0) snoopers(controlingSnooper)%mode=SnooperObserving
             WHERE (snoopers%mode==SnooperRequestingControling) &
                  snoopers%mode=SnooperControling
          END IF
          
          ! Deal with any snoopers about to finish
          IF (ANY(snoopers%mode==SnooperFinising)) THEN 
             ALLOCATE(oldSnoopers(noSnoopers),STAT=status)
             IF (status/=0) CALL MLSMessage(MLSMSG_Errors,ModuleName,&
                  MLSMSG_Allocate//"oldSnoopers")
             oldSnoopers=snoopers
             noSnoopers=COUNT(snoopers%mode/=SnooperFinishing)
             ALLOCATE(snoopers(noSnoopers),STAT=status)
             IF (status/=0) CALL MLSMessage(MLSMSG_Errors,ModuleName,&
                  MLSMSG_Allocate//"snoopers")
             snoopers=PACK(oldSnoopers,snoopers%mode/=SnooperFinishing)
             DEALLOCATE(oldSnoopers,STAT=status)
             IF (status/=0) CALL MLSMessage(MLSMSG_Errors,ModuleName,&
                  MLSMSG_DeAllocate//"oldSnoopers")
          ENDIF
          
          ! This shouldn't be necessary, but just to be sure...
          IF (noSnoopers==0) doneSnoopingForNow= .TRUE.
          
          IF (doneSnoopingForNow) EXIT snoopLoop
       END DO snoopLoop
    END IF
  END SUBROUTINE Snoop

END MODULE SnoopMLSL2
