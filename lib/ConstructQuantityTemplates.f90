! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ConstructQuantityTemplates ! Construct templates from user supplied info
!=============================================================================

  USE MLSCommon
  USE MLSMessageModule
  USE QuantityTemplates
  USE MLSCF
  USE VGrid
  USE HGrid
  USE VerticalCoordinate
  USE Units
  USE MLSSignalNomenclature

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module has various functionality for constructing quantity templates.

CONTAINS
  
  ! First we have some general support routines

  ! --------------------------------------------------------------------------

  ! This routine copies HGrid information into an already defined quantity

  SUBROUTINE CopyHGridInfoIntoQuantity(hGrid,qty)

    ! Dummy arguments
    TYPE (hGrid_T), INTENT(IN) :: hGrid
    TYPE (QuantityTemplate_T), INTENT(INOUT) :: qty

    ! Executable code

    IF (qty%noSubVectors/=hGrid%noProfs) CALL MLSMessage(MLSMSG_Error,&
         & ModuleName,"Size of HGrid not compatible with size of quantity")

    IF (.NOT. qty%stacked) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Cannot copy hGrids into unstacked quantities")

    qty%phi(1,:)=hGrid%phi
    qty%geodLat(1,:)=hGrid%geodLat
    qty%lon(1,:)=hGrid%lon
    qty%time(1,:)=hGrid%time
    qty%solarTime(1,:)=hGrid%solarTime
    qty%solarZenith(1,:)=hGrid%solarZenith
    qty%losAngle(1,:)=hGrid%losAngle

  END SUBROUTINE CopyHGridInfoIntoQuantity
    
  ! --------------------------------------------------------------------------

  ! This similar routine copies VGrid information into an already 
  ! defined quantity

  SUBROUTINE CopyVGridInfoIntoQuantity(vGrid,qty)

    ! Dummy arguments
    TYPE (VGrid_T), INTENT(IN) :: vGrid
    TYPE (QuantityTemplate_T), INTENT(INOUT) :: qty

    ! Executable code

    IF (vGrid%noSurfs/=qty%noSurfs) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Size of vGrid not compatible with size of quantity")

    IF (.NOT. qty%coherent) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Cannot copy vGrid information into incoherent quantities")

    qty%verticalCoordinate=vGrid%verticalCoordinate
    qty%surfs(:,1)=vGrid%surfs

  END SUBROUTINE CopyVGridInfoIntoQuantity

  ! --------------------------------------------------------------------------

  ! This routine constructs a minor frame based quantity.

  SUBROUTINE ConstructMinorFrameQuantity(l1bInfo,chunk,instrumentModule,qty,&
       & noChans,regular,subVectorLen,storeByChannel,mifGeolocation)

    ! Dummy arguments
    TYPE (L1BInfo_T), INTENT(IN) :: l1bInfo ! File handles for l1bdata
    TYPE (MLSChunk_T), INTENT(IN) :: chunk ! The chunk under consideration
    INTEGER, INTENT(IN) :: instrumentModule ! THz or GHz?
    TYPE (QuantityTemplate_T), INTENT(OUT) :: qty ! Resulting quantity
    INTEGER, INTENT(IN), OPTIONAL :: noChans
    LOGICAL, INTENT(IN), OPTIONAL :: regular
    INTEGER, INTENT(IN), OPTIONAL :: subVectorLen
    LOGICAL, INTENT(IN), OPTIONAL :: storeByChannel
    TYPE (QuantityTemplate_T), INTENT(IN), DIMENSION(:), OPTIONAL :: &
         & mifGeolocation

    ! Local parameters, note the similarity to CreateHGridFromMLSCFInfo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IF MODIFYING THIS SECTION PLEASE TAKE CARE, SEE BELOW!
    INTEGER, PARAMETER :: NoL1BItemsToRead=7
    CHARACTER (LEN=15), DIMENSION(NoL1BItemsToRead), &
         &   PARAMETER :: L1bItemsToRead=&
         & (/"time           ","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpLosAngle     "/)
    INTEGER, PARAMETER :: TransitionToModularItems=2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz.  The layout of the above array is critically
    ! bound to the SELECT CASE(l1bItem) code below.  So TAKE CARE! when
    ! modifing it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables

    TYPE (L1BData_T) :: l1bField
    CHARACTER (LEN=NameLen) :: l1bItemName

    INTEGER :: noMAFs,l1bFlag,l1bItem

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    IF (PRESENT(mifGeolocation)) THEN

       ! We have geolocation information, setup the quantity as a clone of that.

       CALL SetupNewQuantityTemplate(qty,&
            & source=mifGeolocation(instrumentModule), &
            & noChans=noChans, regular=regular,subVectorLen=subVectorLen,&
            & storeByChannel=storeByChannel)

       ! Now we're going to deal with a VGrid for this quantity

       qty%verticalCoordinate=VC_Altitude
       qty%surfs=>mifGeolocation(instrumentModule)%surfs

       ! Now we're going to fill in the hGrid information

       qty%time=>       MIFGeolocation(instrumentModule)%time
       qty%geodLat=>    MIFGeolocation(instrumentModule)%geodLat
       qty%lon=>        MIFGeolocation(instrumentModule)%lon
       qty%phi=>        MIFGeolocation(instrumentModule)%phi
       qty%solarZenith=>MIFGeolocation(instrumentModule)%solarZenith
       qty%solarTime=>  MIFGeolocation(instrumentModule)%solarTime
       qty%losAngle=>   MIFGeolocation(instrumentModule)%losAngle
    ELSE
       ! We have no geolocation information, we have to read it ourselves
       ! from the l1boa file.

       ! First we read tpGeodalt to get the size of the quantity.

       l1bItemName=TRIM(MLSInstrumentModuleNames(instrumentModule))//"tpGeodAlt"
       CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,l1bField,noMAFs,l1bFlag, &
            & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex)
       IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_L1BRead//l1bItemName)

       ! Now noMAFs qty%noSubVectors, l1bField%maxMIFs is no surfs.

       CALL SetupNewQuantityTemplate(qty,noSubVectors=noMAFs, &
            & noSurfs=l1bField%maxMIFs, noChans=noChans, coherent=.FALSE.,&
            & stacked=.FALSE.,regular=regular,subVectorLen=subVectorLen,&
            & storeByChannel=storeByChannel,minorFrame=.TRUE.)

       ! Now we're going to deal with a VGrid for this quantity

       qty%verticalCoordinate=VC_Altitude
       qty%surfs=l1bField%dpField(1,:,:)  ! Vert coord is tpGeodAlt read above.

       ! Now we're going to fill in the hGrid information

       DO l1bItem=1,NoL1BItemsToRead
          ! Get the name of the item to read
          l1bItemName=l1bItemsToRead(l1bItem)
          IF (l1bItem>=TransitionToModularItems) l1bItemName=&
               & MLSInstrumentModuleNames(instrumentModule)//l1bItemName

          ! Read it from the l1boa file
          CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,l1bField,noMAFs, &
               & l1bFlag,firstMAF=chunk%firstMafIndex,&
               & lastMAF=chunk%lastMafIndex)
          IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & MLSMSG_L1BRead//l1bItemName)

          ! Now we have to save this field in the quantity data.
          ! This is rather a
          ! kludgy way of doing it but this worked out the least boring way to
          ! write the code.  See the definition of L1BItemsToRead above for
          ! reference.

          SELECT CASE(l1bItem)
          CASE(1)
             qty%time=l1bField%dpField(1,:,:)
          CASE(2)
             qty%geodLat=l1bField%dpField(1,:,:)
          CASE(3)
             qty%lon=l1bField%dpField(1,:,:)
          CASE(4)
             qty%phi=l1bField%dpField(1,:,:)
          CASE(5)
             qty%solarZenith=l1bField%dpField(1,:,:)
          CASE(6)
             qty%solarTime=l1bField%dpField(1,:,:)
          CASE(7)
             qty%losAngle=l1bField%dpField(1,:,:)
          END SELECT
       END DO                      ! Loop over l1b quantities
    ENDIF

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.

  END SUBROUTINE ConstructMinorFrameQuantity

  ! --------------------------------------------------------------------------

  ! This routine constructs a vector quantity template based on instructions
  ! passed in an mlscf line.
  
  SUBROUTINE CreateQtyTemplateFromMLSCFInfo(qty,cfInfo,hGrids,vGrids, &
       & l1bInfo,chunk,mifGeolocation)

    ! Dummy arguments
    TYPE (QuantityTemplate_T), INTENT(OUT) :: qty
    TYPE (MLSCFEntry_T), INTENT(IN) :: cfInfo
    TYPE (HGrid_T), DIMENSION(:), INTENT(IN) :: hGrids
    TYPE (VGrid_T), DIMENSION(:), INTENT(IN) :: vGrids
    TYPE (l1bInfo_T), INTENT(IN) :: l1bInfo
    TYPE (MLSChunk_T), INTENT(IN) :: chunk
    TYPE (QuantityTemplate_T), DIMENSION(:), INTENT(IN), OPTIONAL :: &
         & mifGeolocation

    ! Local parameters
    CHARACTER (LEN=*), PARAMETER :: BadUnitMessage = &
         & "Unexpected or absent unit for "

    ! Local variables

    INTEGER :: quantityType=QTY_Invalid
    CHARACTER (LEN=NameLen) :: name
    INTEGER :: family=PHYQ_Invalid
    REAL(r8) :: scaleFactor
    INTEGER :: hGridIndex=0,vGridIndex=0
    CHARACTER (LEN=NameLen) :: molecule="",radiometer="", band=""
    TYPE (MLSCFCell_T) :: cell
    INTEGER :: keyNo            ! Loop counter
    LOGICAL :: needVHGrids, badUnit
    TYPE(MLSSignal_T), DIMENSION(:), POINTER :: signals
    INTEGER :: instrumentModule
    INTEGER :: noChans=1
    LOGICAL :: storeByChannel=.FALSE.

    ! Executable code

    ! First we'll loop over the mlscf keys and parse them.

    DO keyNo=1,cfInfo%mlscfEntryNoKeys
       cell=cfInfo%cells(keyNo)
       SELECT CASE(TRIM(cell%keyword))
       CASE ("NAME")
          name=cell%charValue
       CASE ("HGRID")
          hGridIndex=LinearSearchStringArray(hGrids%name,cell%charValue,&
               & caseInsensitive=.FALSE.)
          IF (hGridIndex==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & "Unknown hGrid: "//cell%charValue)
       CASE ("VGRID")
          vGridIndex=LinearSearchStringArray(vGrids%name,cell%charValue,&
               & caseInsensitive=.FALSE.)
          IF (vGridIndex==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & "Unknown vGrid: "//cell%charValue)
       CASE ("TYPE")
          quantityType=LinearSearchStringArray(QTYTypeNames,cell%charValue,&
               & caseInsensitive=.FALSE.)
          IF (quantityType==QTY_Invalid) &
               & CALL MLSMessage(MLSMSG_Error,ModuleName,&
               &   "No such quantity type: "//cell%charValue)
       CASE ("UNIT")
          CALL ParseUnitName(cell%charValue,family,scaleFactor)
          IF (family==PHYQ_Invalid) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & "Unrecognised unit: "//cell%charValue)
       CASE("MOLECULE")
          molecule=cell%charValue
       CASE("RADIOMETER")
          molecule=cell%charValue
       CASE("BAND")
          band=cell%charValue
       CASE("STOREBYCHANNEL")
          storeByChannel=(cell%intValue==1)
       END SELECT
    END DO

    ! Now we know what the user asked for, try to make sense of it.
    ! First see that the `type' has been defined

    IF (LEN_TRIM(name) == 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "No name given for quantity template")
    IF (quantityType==QTY_Invalid) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "No type specified for "//name)

    ! Now, depending on the type, check out stuff and see if's ok to first
    ! order.

    badUnit=.FALSE.
    SELECT CASE(quantityType)
    CASE(QTY_Temperature)
       needVHGrids=.TRUE.
       IF (family/=PHYQ_Temperature) badUnit=.TRUE.
    CASE(QTY_Vmr)
       needVHGrids=.TRUE.
       IF (family/=PHYQ_vmr) badUnit=.TRUE.
    CASE(QTY_Radiance)
       needVHGrids=.FALSE.
       IF (family/=PHYQ_Temperature) badUnit=.TRUE.
    CASE(QTY_Ptan)
       IF (family/=PHYQ_Zeta) badUnit=.TRUE.
       needVHGrids=.FALSE.
    CASE(QTY_Baseline)
       IF (family/=PHYQ_Temperature) badUnit=.TRUE.
       needVHGrids=.FALSE.
    CASE(QTY_Extinction)
       ! Need to think about a family here
       needVHGrids=.TRUE.
    CASE DEFAULT
       CALL MLSMessage(MLSMSG_Error,ModuleName,"Unknown quantityType")
    END SELECT

    IF (badUnit) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Incorrect/absent unit for "//name)

    ! Here the code diverges.  In the case where vGrids and hGrids are required
    ! we set them up.  Where they're not we assume it's a minor frame quantity
    ! and work from there.

    IF (needVHGrids) THEN
       ! This is not a minor frame quantity, set it up from VGrids and HGrids
       IF ((hGridIndex==0).OR.(vGridIndex==0)) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName, "Quantity "//TRIM(name)//&
            & " needs a vGrid and/or an hGrid")

       CALL SetupNewQuantityTemplate(qty, &
            & noSubVectors=hGrids(hGridIndex)%noProfs, &
            & noSurfs=vGrids(vGridIndex)%noSurfs, &
            & coherent=.TRUE.,stacked=.FALSE.,regular=.TRUE.)
       ! Note in later versions we'll need to think about channels here

       CALL CopyHGridInfoIntoQuantity(hGrids(hGridIndex),qty)
       CALL CopyVGridInfoIntoQuantity(vGrids(vGridIndex),qty)

    ELSE
       ! This is a minor frame type quantity.
       IF ((hGridIndex/=0).OR.(vGridIndex/=0)) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName, "Quantity "//TRIM(name)//&
            & " has an unnecessary vGrid and/or hGrid")

       SELECT CASE(quantityType)
       CASE(QTY_Ptan)
          CALL ParseMLSSignalRequest(radiometer,signals)
          instrumentModule=signals(1)%instrumentModule
          CALL DestroyMLSSignalsInfo(signals)
       CASE (QTY_Radiance)
          CALL ParseMLSSignalRequest(band,signals)
          IF (SIZE(signals)>1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & "Only one matching signal allowed: "//band)
          noChans=signals(1)%noChannelsInBand
       CASE DEFAULT
          CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Inappropriate quantity for this version")
       END SELECT
       
       ! Construct an empty quantity

       CALL ConstructMinorFrameQuantity(l1bInfo,chunk,instrumentModule,qty, &
            & noChans=noChans,storeByChannel=storeByChannel, &
            & mifGeolocation=mifGeolocation)

       ! Fill what information we can

       IF (quantityType==QTY_Radiance) THEN
          qty%signal=signals(1)
          CALL DestroyMLSSignalsInfo(signals)
       END IF

    ENDIF

    ! Now fill up the remaining items, e.g. name etc.

    qty%name=name
    qty%quantityType=quantityType
    qty%scaleFactor=scaleFactor
    ! Add a unit name field here probably.

  END SUBROUTINE CreateQtyTemplateFromMLSCFInfo

!=============================================================================
END MODULE ConstructQuantityTemplates
!=============================================================================

!
! $Log$
! Revision 1.3  2000/01/11 22:51:34  livesey
! Dealt with ramifications of change from read_parse_l2cf to MLSCF
!
! Revision 1.2  1999/12/18 01:07:00  livesey
! Change vGrids and hGrids from pointer to intent(in)
!
! Revision 1.1  1999/12/18 00:35:40  livesey
! First version
!
!
