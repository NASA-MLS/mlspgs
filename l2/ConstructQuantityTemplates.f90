! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ConstructQuantityTemplates ! Construct templates from user supplied info
!=============================================================================

  ! This module has various functionality for constructing quantity templates.

  use HGrid, only: hGrid_T
  use INIT_TABLES_MODULE, only: F_BAND, F_HGRID, &
    F_MOLECULE, F_RADIOMETER, F_TYPE, F_UNIT, F_VGRID, FIRST_LIT, LAST_LIT, &
    L_BASELINE, L_EXTINCTION, L_GEODALTITUDE, L_GPH, L_PTAN, L_RADIANCE, &
    L_TEMPERATURE, L_TRUE, L_VMR, PHYQ_LENGTH, PHYQ_TEMPERATURE, PHYQ_VMR, &
    PHYQ_ZETA, F_INSTRUMENTMODULE
  use L1BData, only: L1BData_T, READL1BDATA
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Info, MLSMSG_L1BRead
  use MLSSignalNomenclature, only: DestroyMLSSignalsInfo, MLSSignal_T, &
    & ParseMLSSignalRequest
  use OUTPUT_M, only: OUTPUT
  use QuantityTemplates, only: QuantityTemplate_T,SetupNewQuantityTemplate
  use STRING_TABLE, only: GET_STRING
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_SET_ONE
  use VGrid, only: VGrid_T
  USE Intrinsic, ONLY: L_NONE, L_THz, L_GHz
  use init_tables_module, ONLY: LIT_INDICES
  implicit none
  private
  public :: ConstructMinorFrameQuantity, CreateQtyTemplateFromMLSCFInfo

! -----     Private declarations     -----------------------------------

  integer :: ERROR

! Error codes for "announce_error"
  integer, parameter :: BadUnitMessage = 1
  integer, parameter :: InappropriateQuantity = 2
  integer, parameter :: NeedGrid = 3
  integer, parameter :: NoQuantityType = 4
  integer, parameter :: UnnecessaryGrid = 5

  !---------------------------- RCS Ident Info -------------------------------
  character (len=256) :: Id = &
       "$Id$"
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================
  ! -----------------------------  CreateQtyTemplateFromMLSCFInfo  -----
  type (QuantityTemplate_T) function CreateQtyTemplateFromMLSCFInfo ( &
    & name, root, hGrids, vGrids, l1bInfo, chunk, mifGeolocation ) &
    result ( QTY )

  ! This routine constructs a vector quantity template based on instructions
  ! passed in an mlscf line.

    ! Dummy arguments
    integer, intent(in) :: NAME              ! Sub-rosa index of name
    integer, intent(in) :: ROOT              ! Root of QuantityTemplate subtree
    type (HGrid_T), dimension(:), intent(in) :: hGrids
    type (VGrid_T), dimension(:), intent(in) :: vGrids
    type (l1bInfo_T), intent(in) :: l1bInfo
    type (MLSChunk_T), intent(in) :: chunk
    type (QuantityTemplate_T), dimension(:), intent(in), optional :: &
      & mifGeolocation

    ! Local variables

    logical :: BADUNIT
    integer :: BAND           ! String index of BAND= value
    integer :: FAMILY
    logical :: FIRSTINDEXCHANNEL
    integer :: HGRIDINDEX
    integer :: I              ! Loop counter
    integer :: INSTRUMENTMODULE
    integer :: KEY            ! Field name, F_... from Init_Tables_Module
    character(len=127) :: LIT_TEXT
    integer :: MOLECULE
    integer :: NATURAL_UNITS(first_lit:last_lit)
    integer :: NOINSTANCES
    integer :: NOSURFS
    logical :: minorFrame       ! Is a minor frame quantity
    integer :: NOCHANS
    integer :: QUANTITYTYPE
    integer :: RADIOMETER     ! String index of RADIOMETER= value
    real(r8) :: SCALEFACTOR
    type (MLSSignal_T), dimension(:), pointer :: SIGNALS
    integer :: SON            ! A Son of Root -- an n_assign node
    integer :: TYPE_FIELD     ! Index in subtree of "type"
    integer :: VALUE          ! Node index of value of field of spec
    integer :: VGRIDINDEX

    ! Executable code

    if ( toggle(gen) ) &
      & call trace_begin ( "CreateQtyTemplateFromMLSCFInfo", root )

! ??? Do we need a GOT_FIELD check like in VGrid, e.g. ???

    error = 0
    family = 0
    hGridIndex = 0
    vGridIndex = 0
    molecule = 0
    natural_units = 0
    natural_units( (/     l_baseline,       l_extinction, &
      & l_gph,            l_ptan,           l_radiance, &
      & l_temperature,    l_vmr /) ) = &
                   (/     PHYQ_Temperature, PHYQ_Temperature, &
      & PHYQ_length,      PHYQ_Zeta,        PHYQ_Temperature, &
      & PHYQ_Temperature, PHYQ_Vmr /)
    noChans = 1
    quantitytype = 0
    radiometer = 0
    instrumentModule= 0
    scaleFactor = 1.0

    ! First we'll loop over the MLSCF keys.

    do i = 2, nsons(root)
      son = subtree(i,root)
      key = subtree(1,son)
      if ( node_id(key) == n_set_one ) then
        key = subtree(1,key)
        value = l_true
      else
        value = decoration(subtree(2,son))
      end if

      select case ( decoration(key) )
      case ( f_hgrid )
        hGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_vgrid )
        vGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_type )
        quantityType = value
        type_field = son
      case ( f_unit );              scaleFactor = value
      case ( f_molecule );          molecule = value
      case ( f_radiometer );        radiometer = sub_rosa(subtree(2,son))
      case ( f_instrumentmodule);   instrumentModule = sub_rosa(subtree(2,son))
      case ( f_band );              band = sub_rosa(subtree(2,son))
      end select
    end do

    ! Now we know what the user asked for, try to make sense of it.
    ! First see that the `type' has been defined

    if ( quantityType == 0 ) call announce_error ( root, noQuantityType )

    ! Now, depending on the type, check out stuff and see if it's ok to
    ! first order.

    if ( family == 0 ) family = natural_units(quantityType)
    select case ( quantityType )
    case ( l_Baseline )
      minorFrame=.FALSE.
    case ( l_Extinction )
      ! ??? Need to think about a family here
      minorFrame=.FALSE.
    case ( l_Gph )
      ! ??? Need to think about a family here
      minorFrame=.FALSE.
    case ( l_Ptan )
      minorFrame=.TRUE.
    case ( l_Radiance )
      minorFrame=.TRUE.
    case ( l_Temperature )
      minorFrame=.FALSE.
    case ( l_Vmr )
      minorFrame=.FALSE.
    case default ! Can't get here if tree_walker works correctly
    end select

    ! Here the code splits, for minor frame quantities, we take the information
    ! from the previously constructed MIFGeolocation information.  Otherwise,
    ! we'll probably need to use any supplied vGrid/hGrid information.

    if ( minorFrame ) then

      ! This is a minor frame type quantity.
      if ( (hGridIndex/=0) .OR. (vGridIndex/=0)) &
        &  call announce_error ( root, unnecessaryGrid )

      select case ( quantityType )
      case ( l_Ptan )
        call get_string ( radiometer, lit_text, cap=.true. )
        call ParseMLSSignalRequest ( lit_text(2:len_trim(lit_text)-1), signals )
        instrumentModule = signals(1)%instrumentModule
        call DestroyMLSSignalsInfo ( signals )
      case ( l_Radiance )
        call get_string ( band, lit_text, cap=.true. )
        call ParseMLSSignalRequest ( lit_text(2:len_trim(lit_text)-1), signals )
        if ( SIZE(signals)>1) call MLSMessage(MLSMSG_Error,ModuleName,&
             & "Only one matching signal allowed: "//lit_text )
        noChans = signals(1)%noChannelsInBand
!       call DestroyMLSSignalsInfo ( signals ) ! Done later
      case default
        call announce_error ( type_field, inappropriateQuantity )
      end select

      ! Construct an empty quantity
      call ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
        & qty, noChans=noChans, mifGeolocation=mifGeolocation )

      ! Fill what information we can

      if ( quantityType==l_Radiance ) then
        qty%signal = signals(1)
        call DestroyMLSSignalsInfo ( signals )
      end if

    else

      ! This is not a minor frame quantity, set it up from VGrids and HGrids

      if ( hGridIndex/=0 ) then
        noInstances=hGrids(hGridIndex)%noProfs
      else
        noInstances=1
      endif

      if ( vGridIndex/=0 ) then
        noSurfs=vGrids(vGridIndex)%noSurfs
      else
        noSurfs=1
      endif

      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.TRUE., stacked=.TRUE., regular=.TRUE. )
      ! ??? Note in later versions we'll need to think about channels here

      qty%frequencyCoordinate = L_None

      if (hGridIndex /=0 ) then
        call CopyHGridInfoIntoQuantity ( hGrids(hGridIndex), qty )
      else                      ! Set `empty' values
        qty%phi=0.0
        qty%geodLAt=0.0
        qty%lon=0.0
        qty%time=0.0
        qty%solarTime=0.0
        qty%solarZenith=0.0
        qty%losAngle=0.0
      endif

      if (vGridIndex /=0 ) then
        call CopyVGridInfoIntoQuantity ( vGrids(vGridIndex), qty )
      else
        qty%surfs=0.0
        qty%verticalCoordinate=L_None
      endif

    end if

    ! Now fill up the remaining items, e.g. name etc.

    qty%unit = family
    qty%molecule = molecule
    qty%name = name
    qty%quantityType = quantityType
    qty%instrumentmodule= instrumentmodule
    qty%radiometer = radiometer
    qty%scaleFactor = scaleFactor

    if ( toggle(gen) ) call trace_end ( "CreateQtyTemplateFromMLSCFInfo" )

  end function CreateQtyTemplateFromMLSCFInfo

! =====     Private Procedures     =====================================

! -----------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( badUnitMessage )
      call output ( "Incorrect or absent unit.", advance='yes' )
    case ( InappropriateQuantity )
      call output ( "A quantity inappropriate for this version is specified.", &
                    advance='yes' )
    case ( needGrid )
      call output ( "Quantity needs a vGrid or hGrid.", advance='yes' )
    case ( noQuantityType )
      call output ( "No quantity type specified.", advance='yes' )
    case ( unnecessaryGrid )
      call output ( "Quantity has an unnecessary vGrid or hGrid.", advance='yes' )
    end select
  end subroutine ANNOUNCE_ERROR

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  subroutine ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
    & qty, noChans, regular, instanceLen, mifGeolocation )

  ! This routine constructs a minor frame based quantity.

    ! Dummy arguments
    type (L1BInfo_T), intent(in) :: l1bInfo ! File handles for l1bdata
    type (MLSChunk_T), intent(in) :: chunk ! The chunk under consideration
    integer, intent(in) :: instrumentModule ! L_THz or L_GHz?
    type (QuantityTemplate_T), intent(out) :: qty ! Resulting quantity
    integer, intent(in), optional :: noChans
    logical, intent(in), optional :: regular
    integer, intent(in), optional :: instanceLen
    type (QuantityTemplate_T), intent(in), dimension(:), optional :: &
         & mifGeolocation

    ! Local parameters
    real(r8), parameter :: SIXTH = 1.0_r8 / 6.0_r8
    ! Note the similarity to CreateHGridFromMLSCFInfo:
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IF MODIFYING THIS SECTION PLEASE TAKE CARE, SEE BELOW!
    integer, parameter :: NoL1BItemsToRead=7
    character (len=15), dimension(NoL1BItemsToRead), parameter :: &
         & L1bItemsToRead = &
         & (/"MAFStartTimeTAI","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpLosAngle     "/)
    integer, parameter :: TransitionToModularItems = 2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz.  The layout of the above array is critically
    ! bound to the select case(l1bItem) code below.  So TAKE CARE! when
    ! modifing it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables

    type (L1BData_T) :: l1bField
    character (len=NameLen) :: l1bItemName

    integer :: noMAFs, l1bFlag, l1bItem, mafIndex, mifIndex, instrumentModuleIndex
    integer, DIMENSION(1) :: instrumentModuleIndexArray

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    if ( present(mifGeolocation) ) then

      ! We have geolocation information, setup the quantity as a clone of that.
       instrumentModuleIndexArray=PACK( (/1,2/), &
            instrumentModule==mifGeolocation%instrumentModule)
       instrumentModuleIndex=instrumentModuleIndexArray(1)

      call SetupNewQuantityTemplate ( qty, &
        & source=mifGeolocation(instrumentModuleIndex), &
        & noChans=noChans,  regular=regular, instanceLen=instanceLen )

      ! Now we're going to deal with a VGrid for this quantity

      qty%verticalCoordinate = l_geodAltitude
      qty%surfs=>mifGeolocation(instrumentModuleIndex)%surfs

      ! Now we're going to fill in the hGrid information

      qty%time =>        MIFGeolocation(instrumentModuleIndex)%time
      qty%geodLat =>     MIFGeolocation(instrumentModuleIndex)%geodLat
      qty%lon =>         MIFGeolocation(instrumentModuleIndex)%lon
      qty%phi =>         MIFGeolocation(instrumentModuleIndex)%phi
      qty%solarZenith => MIFGeolocation(instrumentModuleIndex)%solarZenith
      qty%solarTime =>   MIFGeolocation(instrumentModuleIndex)%solarTime
      qty%losAngle =>    MIFGeolocation(instrumentModuleIndex)%losAngle
      qty%MAFCounter =>  MIFGeolocation(instrumentModuleIndex)%mafCounter
      qty%mafIndex =>    MIFGeolocation(instrumentModuleIndex)%mafIndex
      qty%noInstancesLowerOverlap = &
        & MIFGeolocation(instrumentModuleIndex)%noInstancesLowerOverlap
      qty%noInstancesUpperOverlap = &
        & MIFGeolocation(instrumentModuleIndex)%noInstancesUpperOverlap
    else
      ! We have no geolocation information, we have to read it ourselves
      ! from the L1BOA file.

      ! First we read tpGeodalt to get the size of the quantity.

      CALL Get_String(lit_indices(instrumentModule),l1bItemName)
      l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
      
      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
      if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )

      ! Now noMAFs qty%noInstances, l1bField%maxMIFs is no surfs.

      call SetupNewQuantityTemplate ( qty, noInstances=noMAFs, &
        & noSurfs=l1bField%maxMIFs, noChans=noChans, coherent=.FALSE., &
        & stacked=.FALSE., regular=regular, instanceLen=instanceLen, &
        & minorFrame=.TRUE. )

      qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
      qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap

      ! Now we're going to deal with a VGrid for this quantity

      qty%verticalCoordinate = l_geodAltitude
      qty%surfs = l1bField%dpField(1,:,:)  ! Vert coord is tpGeodAlt read above.
      qty%mafCounter = l1bField%counterMAF
      do mafIndex = chunk%firstMAFIndex, chunk%lastMAFIndex
        qty%mafIndex(mafIndex-chunk%firstMAFIndex+1) = mafIndex
      end do

      ! Now we're going to fill in the hGrid information

      do l1bItem = 1, NoL1BItemsToRead
        ! Get the name of the item to read
        l1bItemName=l1bItemsToRead(l1bItem)
        if ( l1bItem>=TransitionToModularItems ) then
           call Get_String(lit_indices(instrumentModule),l1bItemName)
           l1bItemName = trim(l1bItemName)//'.'//l1bItemsToRead(l1bItem)
        else
           l1bItemName=l1bItemsToRead(l1bItem)
        endif

        ! Read it from the l1boa file
        call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
          & l1bFlag, firstMAF=chunk%firstMafIndex, &
          & lastMAF=chunk%lastMafIndex )
        if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )

        ! Now we have to save this field in the quantity data.
        ! This is rather a kludgy way of doing it but this worked out the
        ! least boring way to write the code.  See the definition of
        ! L1BItemsToRead above for reference.

        select case(l1bItem)
        case(1)
          ! For time we have to do something a little more complicated.
          ! ******* This is a real kludge, and we have to find a way
          ! to do it better in 0.5. Probably simply have time as a minor
          ! frame quantity in L1, or MIF duration.
          !
          do mafIndex = 1, noMAFs
            do mifIndex = 1, l1bField%maxMIFs
              qty%time(mifIndex,mafIndex) = &
                & l1bField%dpField(1,1,mafIndex) + &
                & (mifIndex-1) * sixth
            end do
          end do
        case(2)
          qty%geodLat=l1bField%dpField(1,:,:)
        case(3)
          qty%lon=l1bField%dpField(1,:,:)
        case(4)
          qty%phi=l1bField%dpField(1,:,:)
        case(5)
          qty%solarZenith=l1bField%dpField(1,:,:)
        case(6)
          qty%solarTime=l1bField%dpField(1,:,:)
        case(7)
          qty%losAngle=l1bField%dpField(1,:,:)
        end select
      end do                      ! Loop over l1b quantities
    endif

    qty%frequencyCoordinate = L_None

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.

  end subroutine ConstructMinorFrameQuantity

  ! ----------------------------------  CopyHGridInfoIntoQuantity  -----
  subroutine CopyHGridInfoIntoQuantity ( hGrid, qty )

  ! This routine copies HGrid information into an already defined quantity

    ! Dummy arguments
    type (hGrid_T), intent(in) :: hGrid
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code

    if ( qty%noInstances/=hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )

    if ( .NOT. qty%stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy hGrids into unstacked quantities")

    qty%phi(1,:) = hGrid%phi
    qty%geodLat(1,:) = hGrid%geodLat
    qty%lon(1,:) = hGrid%lon
    qty%time(1,:) = hGrid%time
    qty%solarTime(1,:) = hGrid%solarTime
    qty%solarZenith(1,:) = hGrid%solarZenith
    qty%losAngle(1,:) = hGrid%losAngle
    qty%noInstancesLowerOverlap = hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = hGrid%noProfsUpperOverlap

  end subroutine CopyHGridInfoIntoQuantity

  ! ----------------------------------  CopyVGridInfoIntoQuantity  -----
  subroutine CopyVGridInfoIntoQuantity ( vGrid, qty )

  ! This similar routine copies VGrid information into an already 
  ! defined quantity

    ! Dummy arguments
    type (VGrid_T), intent(in) :: vGrid
    type (QuantityTemplate_T), intent(inout) :: qty

    ! Executable code

    if ( vGrid%noSurfs /= qty%noSurfs ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, "Size of vGrid not compatible with size of quantity" )

    if ( .NOT. qty%coherent ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy vGrid information into incoherent quantities")

    qty%verticalCoordinate = vGrid%verticalCoordinate
    qty%surfs(:,1) = vGrid%surfs

  end subroutine CopyVGridInfoIntoQuantity


!=============================================================================
end module ConstructQuantityTemplates
!=============================================================================

!
! $Log$
! Revision 2.4  2001/02/21 01:09:00  livesey
! Allowed for quantities with no h/v grid
!
! Revision 2.3  2001/02/20 18:43:50  livesey
! Removed all references to firstIndexChannel
!
! Revision 2.2  2001/02/14 00:12:45  livesey
! Removed firstIndexChannel
!
! Revision 2.1  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.0  2000/09/05 18:57:02  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:03  vsnyder
! Initial entry
!
! Revision 1.13  2000/05/17 23:19:07  lungu
! Added "." between MLSInstrumentModuleNames and l1bItemName.
! Made hGridIndex=0 and vGridIndex=0 upon entry, so it does not "inherit" attributes from previous call.
! Made caseInsensitive=.TRUE. for all searches.
! Added type for QTY_Gph.
! Made stacked=.TRUE. so that CopyHGridInfoIntoQuantity works.
!
