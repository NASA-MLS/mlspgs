! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ConstructQuantityTemplates ! Construct templates from user supplied info
!=============================================================================

  ! This module has various functionality for constructing quantity templates.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: Dump
  use EXPR_M, only: EXPR
  use HGrid, only: hGrid_T
  use INIT_TABLES_MODULE, only: F_BAND, F_GEODANGLE, F_HGRID, F_INCLINATION, &
    & F_LOGBASIS, F_MODULE, F_MOLECULE, F_NOMIFS, F_RADIOMETER, &
    & F_SIGNAL, F_TYPE, F_UNIT, F_VGRID
  use INIT_TABLES_MODULE, only: FIELD_FIRST, FIELD_LAST, &
    FIRST_LIT, LAST_LIT, L_BASELINE, L_CHANNEL, L_CloudIce,  L_EARTHREFL, &
    L_ELEVOFFSET, L_EXTINCTION, L_GEODALTITUDE, L_GPH, &
    L_HEIGHTOFFSET, L_LOSTRANSFUNC, L_LOSVEL, L_NONE, L_ORBITINCLINATION, &
    L_PTAN, L_RADIANCE, &
    L_REFGPH, L_SCANRESIDUAL, L_SCECI, L_SCGEOCALT, L_SCVEL, L_SIDEBANDRATIO, &
    L_SPACERADIANCE, &
    L_TEMPERATURE, L_TNGTECI, L_TNGTGEOCALT, L_TNGTGEODALT, L_TRUE,&
    L_VMR, L_XYZ, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, PHYQ_EXTINCTION, &
    PHYQ_IceDensity, PHYQ_LENGTH, &
    PHYQ_TEMPERATURE, PHYQ_VELOCITY, PHYQ_VMR, PHYQ_ZETA
  use L1BData, only: L1BData_T, READL1BDATA, DEALLOCATEL1BDATA
  use LEXER_CORE, only: PRINT_SOURCE
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Info, MLSMSG_L1BRead
  use MLSSignals_m, only:  IsModuleSpacecraft, &
    & GetModuleFromRadiometer, GetModuleFromSignal, GetModuleName, &
    & GetRadiometerName, GetRadiometerFromSignal, GetSignal, GetSignalName,&
    & Signal_T
  use OUTPUT_M, only: OUTPUT
  use Parse_Signal_m, only: PARSE_SIGNAL
  use QuantityTemplates, only: QuantityTemplate_T,SetupNewQuantityTemplate, &
    & QuantityTemplateCounter
  use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
  use TOGGLES, only: GEN, TOGGLE
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TREE, only: DECORATION, NODE_ID, NSONS, SOURCE_REF, SUB_ROSA, SUBTREE
  use TREE_TYPES, only: N_SET_ONE
  use Units, only: DEG2RAD, RAD2DEG
  use VGridsDatabase, only: VGrid_T

  implicit none
  private
  public :: ConstructMinorFrameQuantity, CreateQtyTemplateFromMLSCFInfo
  public :: ForgeMinorFrames

! -----     Private declarations     -----------------------------------

  integer :: ERROR

! Error codes for "announce_error"
  integer, parameter :: BadUnitMessage = 1
  integer, parameter :: InappropriateQuantity = 2
  integer, parameter :: NeedGrid = 3
  integer, parameter :: NoQuantityType = 4
  integer, parameter :: UnnecessaryGrid = 5
  integer, parameter :: noModule=6

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================
  ! -----------------------------  CreateQtyTemplateFromMLSCFInfo  -----
  type (QuantityTemplate_T) function CreateQtyTemplateFromMLSCFInfo ( &
    & Name, Root, HGrids, VGrids, L1bInfo, Chunk, MifGeolocation ) &
    result ( QTY )

  ! This routine constructs a vector quantity template based on instructions
  ! passed in an mlscf line.

    ! Dummy arguments
    integer, intent(in) :: NAME              ! Sub-rosa index of name
    integer, intent(in) :: ROOT              ! Root of QuantityTemplate subtree
    type (HGrid_T), dimension(:), pointer :: HGrids
    type (VGrid_T), dimension(:), pointer :: VGrids
    type (l1bInfo_T), intent(in) :: L1bInfo
    type (MLSChunk_T), intent(in) :: Chunk
    type (QuantityTemplate_T), dimension(:), intent(in), optional :: &
      & MifGeolocation

    ! Local variables

    logical :: BadUnit
    integer :: Band                     ! Tree index
    integer :: Family
    integer :: FrequencyCoordinate
    integer :: HGridIndex
    integer :: I                        ! Loop counter
    integer :: InstrumentModule         ! Database index
    integer :: Key            ! Field name, F_... from Init_Tables_Module
    character(len=127) :: SIGNALSTRING
    logical :: LOGBASIS                 ! To place in quantity
    integer :: Molecule
    integer :: Natural_Units(first_lit:last_lit)
    integer :: NoInstances
    integer :: NoSurfs
    logical :: MinorFrame               ! Is a minor frame quantity
    integer :: NoChans
    integer :: QuantityType
    integer :: Radiometer               ! Database index
    real(r8) :: ScaleFactor
    integer :: Sideband
    integer :: Signal                   ! Database index
    integer, dimension(:), pointer :: SignalInds ! From parse signal
    type (signal_T) :: SignalInfo       ! Details of the appropriate signal
    integer :: Son                      ! A Son of Root -- an n_assign node
    character (len=80) :: Str
    integer :: Type_Field               ! Index in subtree of "type"
    integer :: Value                    ! Node index of value of field of spec
    integer :: VGridIndex

    ! Executable code

    if ( toggle(gen) ) &
      & call trace_begin ( "CreateQtyTemplateFromMLSCFInfo", root )

! ??? Do we need a GOT_FIELD check like in VGrid, e.g. ???

    nullify ( signalInds )
    error = 0
    family = 0
    hGridIndex = 0
    instrumentModule = 0
    logBasis = .false.
    molecule = 0

    natural_units = 0
    natural_units(l_baseline) =       PHYQ_Temperature
    natural_units(l_cloudice) =       PHYQ_IceDensity
    natural_units(l_earthRefl) =      PHYQ_Dimensionless
    natural_units(l_elevOffset) =     PHYQ_Angle
    natural_units(l_extinction) =     PHYQ_Extinction
    natural_units(l_gph) =            PHYQ_Length
    natural_units(l_heightOffset ) =  PHYQ_Length
    natural_units(l_losTransFunc) =   PHYQ_Dimensionless
    natural_units(l_losVel) =         PHYQ_Velocity
    natural_units(l_orbitInclination) =   PHYQ_Angle
    natural_units(l_ptan) =           PHYQ_Zeta
    natural_units(l_radiance) =       PHYQ_Temperature
    natural_units(l_refGPH) =         PHYQ_Length
    natural_units(l_scGeocAlt ) =     PHYQ_Length
    natural_units(l_scVel) =          PHYQ_Velocity
    natural_units(l_scanResidual ) =  PHYQ_Length
    natural_units(l_spaceRadiance) =  PHYQ_Temperature
    natural_units(l_temperature) =    PHYQ_Temperature
    natural_units(l_tngtGeocAlt) =    PHYQ_Length
    natural_units(l_tngtGeodAlt) =    PHYQ_Length
    natural_units(l_vmr) =            PHYQ_Vmr

    noChans = 1
    quantitytype = 0
    radiometer = 0
    scaleFactor = 1.0
    sideband = 0
    signal = 0
    vGridIndex = 0

    ! First we'll loop over the MLSCF keys.

    do i = 2, nsons(root)
      son = subtree(i,root)
      key = subtree(1,son)
      if ( node_id(son) == n_set_one ) then
        value = l_true
      else
        value = decoration(subtree(2,son))
      end if

      select case ( decoration(key) )
      case ( f_hgrid )
        hGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_logBasis )
        logBasis = (value == l_true)
      case ( f_module)
        instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_molecule );          molecule = value
      case ( f_radiometer )
        radiometer = decoration(decoration(subtree(2,son)))
        instrumentModule = GetModuleFromRadiometer(radiometer)
      case ( f_signal )
        !??? For the moment it is simple, later we'll be more intelligent here
        !??? for example, letting the user choose either R1A or R1B.
        call get_string( sub_rosa(subtree(2,son)), signalString, strip=.true. )
        !??? Here we would do intelligent stuff to work out which bands
        !??? are present, for the moment choose the first
        call parse_Signal ( signalString, signalInds, &
          & tree_index=son, sideband=sideband )
        if ( .not. associated(signalInds) ) then ! A parse error occurred
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Unable to parse signal string' )
        end if
        signal = signalInds(1)
        call deallocate_test ( signalInds, 'signalInds', ModuleName )
        instrumentModule = GetModuleFromSignal(signal)
        radiometer = GetRadiometerFromSignal(signal)
      case ( f_type )
        quantityType = value
        type_field = son
      case ( f_unit );              scaleFactor = value
      case ( f_vgrid )
        vGridIndex = decoration(value) ! node_id(value) == n_spec_args
      end select
    end do

    ! Now we know what the user asked for, try to make sense of it.
    ! First see that the `type' has been defined

    if ( quantityType == 0 ) call announce_error ( root, noQuantityType )

    ! Now, depending on the type, check out stuff and see if it's ok to
    ! first order.

    if ( family == 0 ) family = natural_units(quantityType)
    minorFrame = any(quantityType == (/ l_Baseline, l_Ptan, l_Radiance, &
      & l_tngtECI, l_tngtGeodAlt, l_tngtGeocAlt, l_scECI, l_scGeocAlt,&
      & l_scVel, l_losTransFunc, l_losVel, l_heightOffset, l_scanResidual/) )

    ! Set defaults for other parameters
    frequencyCoordinate = L_None
    noChans = 1

    ! Here the code splits, for minor frame quantities, we take the information
    ! from the previously constructed MIFGeolocation information.  Otherwise,
    ! we'll probably need to use any supplied vGrid/hGrid information.

    if ( minorFrame ) then

      ! This is a minor frame type quantity.
      if ( ((hGridIndex /= 0) .OR. (vGridIndex /= 0 )) &
         &  .AND. quantityType /= l_losTransFunc ) then
        call announce_error ( root, unnecessaryGrid )
      end if
      if ( instrumentModule == 0 ) then 
        call announce_error ( root, noModule )
      end if

      ! Work out the channel information
      if ( signal /= 0 ) then
        signalInfo = GetSignal(signal)
        noChans = size(signalInfo%frequencies)
        frequencyCoordinate = l_channel
      end if
      
      ! For some cases we know the quantity is an xyz vector
      if ( any(quantityType == (/ l_tngtECI, l_scECI, l_scVel /)) ) then
        noChans = 3
        frequencyCoordinate = l_xyz
      end if

      ! For quantity type l_losTransFunc
      if ( quantityType == l_losTransFunc ) then
        noChans = vGrids(vGridIndex)%noSurfs
        frequencyCoordinate = l_losTransFunc
      end if

      ! Construct an empty quantity
      call ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
        & qty, noChans=noChans, mifGeolocation=mifGeolocation )
    else

      ! This is not a minor frame quantity, set it up from VGrids and HGrids

      if ( hGridIndex/=0 ) then
        noInstances=hGrids(hGridIndex)%noProfs
      else
        noInstances=1
      end if

      if ( vGridIndex/=0 ) then
        noSurfs=vGrids(vGridIndex)%noSurfs
      else
        noSurfs=1
      end if
      
      ! Some special cases for certain quantities
      select case (quantityType)
      case (l_sidebandRatio)
        frequencyCoordinate = l_channel
        signalInfo = GetSignal(signal)
        noChans = size ( signalInfo%frequencies ) 
      case default
      end select

      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.TRUE., stacked=.TRUE., regular=.TRUE.,&
        & noChans=noChans )
      ! ??? Note in later versions we'll need to think about channels here

      if ( hGridIndex /=0 ) then
        call CopyHGridInfoIntoQuantity ( hGrids(hGridIndex), qty )
      else                      ! Set `empty' values
        qty%phi = 0.0
        qty%geodLAt = 0.0
        qty%lon = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
      end if

      if ( vGridIndex /= 0 ) then
        call CopyVGridInfoIntoQuantity ( vGrids(vGridIndex), qty )
      else
        qty%surfs = 0.0
        qty%verticalCoordinate = L_None
      end if

    end if

    ! Now fill up the remaining items, e.g. name etc.

    qty%badValue = -999.99              ! Think more about this later NJL !????
    qty%frequencyCoordinate = frequencyCoordinate
    qty%instrumentmodule = instrumentmodule
    qty%logBasis = logBasis
    qty%molecule = molecule
    qty%name = name
    qty%quantityType = quantityType
    qty%radiometer = radiometer
    qty%scaleFactor = scaleFactor
    qty%sideband = sideband
    qty%signal = signal
    qty%unit = family

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
    case ( noModule )
      call output ( "No module given or deducible.", advance='yes' )
    end select
  end subroutine ANNOUNCE_ERROR

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  subroutine ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
    & qty, noChans, regular, instanceLen, mifGeolocation )

  ! This routine constructs a minor frame based quantity.

    ! Dummy arguments
    type (L1BInfo_T), intent(in) :: l1bInfo  ! File handles for l1bdata
    type (MLSChunk_T), intent(in) :: chunk   ! The chunk under consideration
    integer, intent(in) :: instrumentModule  ! Database index
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
    integer :: noMAFs, l1bFlag, l1bItem, mafIndex, mifIndex

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    if ( present(mifGeolocation) ) then
      ! -------------------------------------- Got mifGeolocation ------------
      ! We have geolocation information, setup the quantity as a clone of that.
      qty = mifGeolocation(instrumentModule)
      qty%noChans = noChans
      qty%instanceLen = qty%noChans*qty%noSurfs
      ! Increment the id counter and set the id field
      quantityTemplateCounter = quantityTemplateCounter + 1
      qty%id = quantityTemplateCounter

    else
      ! -------------------------------------- Not Got mifGeolocation ------------
      ! We have no geolocation information, we have to read it ourselves
      ! from the L1BOA file.

      ! First we read xxGeodalt to get the size of the quantity.
      if ( IsModuleSpacecraft(instrumentModule) ) then
        l1bItemName = 'scGeocAlt' 
      else
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
      endif

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

      if ( .not. IsModuleSpacecraft(instrumentModule) ) then
        call GetModuleName ( instrumentModule, l1bItemName )
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
        call DeallocateL1BData(l1bfield)

        ! Now we're going to fill in the hGrid information
        
        do l1bItem = 1, NoL1BItemsToRead
          ! Get the name of the item to read
          l1bItemName = l1bItemsToRead(l1bItem)
          if ( l1bItem >= TransitionToModularItems ) then
            call GetModuleName ( instrumentModule, l1bItemName )
            l1bItemName = trim(l1bItemName)//'.'//l1bItemsToRead(l1bItem)
          else
            l1bItemName = l1bItemsToRead(l1bItem)
          end if
          
          ! Read it from the l1boa file
          call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
            & l1bFlag, firstMAF=chunk%firstMafIndex, &
            & lastMAF=chunk%lastMafIndex )
          if ( l1bFlag == -1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & MLSMSG_L1BRead//l1bItemName )
          
          ! Now we have to save this field in the quantity data.
          ! This is rather a kludgy way of doing it but this worked out the
          ! least boring way to write the code.  See the definition of
          ! L1BItemsToRead above for reference.
          
          select case(l1bItem)
          case ( 1 )
            !??? For time we have to do something a little more complicated.
            !??? This is a real kludge, and we have to find a way
            !??? to do it better in 0.5. Probably simply have time as a minor
            !??? frame quantity in L1, or MIF duration. !???????
            !??? Also note that it fills in times even for non existant MIFs
            do mafIndex = 1, noMAFs
              do mifIndex = 1, qty%noSurfs 
                qty%time(mifIndex,mafIndex) = &
                  & l1bField%dpField(1,1,mafIndex) + &
                  & (mifIndex-1) * sixth
              end do
            end do
          case ( 2 )
            qty%geodLat = l1bField%dpField(1,:,:)
          case ( 3 )
            qty%lon = l1bField%dpField(1,:,:)
          case ( 4 )
            qty%phi = l1bField%dpField(1,:,:)
          case ( 5 )
            qty%solarZenith = l1bField%dpField(1,:,:)
          case ( 6 )
            qty%solarTime = l1bField%dpField(1,:,:)
          case ( 7 )
            qty%losAngle = l1bField%dpField(1,:,:)
          end select

          call DeallocateL1BData ( l1bField )
        end do                          ! Loop over l1b quantities
      else                              ! Spacecraft module
        ! just zero stuff out.
        qty%surfs = 0.0
        qty%phi = 0.0
        qty%geodLat = 0.0
        qty%lon = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
        
        qty%mafCounter = l1bField%counterMAF
        do mafIndex = chunk%firstMAFIndex, chunk%lastMAFIndex
          qty%mafIndex(mafIndex-chunk%firstMAFIndex+1) = mafIndex
        end do
        
        call DeallocateL1BData(l1bfield)
        
      end if
      
    end if
    qty%frequencyCoordinate = L_None
    qty%instrumentModule = instrumentModule

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.

  end subroutine ConstructMinorFrameQuantity

  ! ------------------------------------ ForgeMinorFrames --------------
  subroutine ForgeMinorFrames ( root, chunk, mifGeolocation )
    ! This routine is used when we're in the mode of creating l2pc files
    ! and we want to invent a set of minor frame quantities with no
    ! reference to the l1 files
    
    ! Dummy arguments
    integer, intent(in) :: ROOT         ! Tree vertex
    type (MLSChunk_T), intent(in) :: CHUNK ! This chunk
    type (QuantityTemplate_T), dimension(:), intent(inout) :: MIFGEOLOCATION

    ! Local variables
    integer :: GEODANGLENODE            ! Tree vertex
    integer :: I                        ! Loop counter
    integer :: INSTRUMENTMODULE         ! Database index
    integer :: KEY                      ! Tree vertex
    integer :: MAF                      ! Loop counter
    integer :: NOMAFS                   ! Dimension
    integer :: NOMIFS                   ! Dimension
    integer :: SON                      ! Tree vertex

    integer, dimension(2) :: EXPR_UNITS ! From tree

    real(r8), dimension(2) :: EXPR_VALUE ! From tree
    real(r8) :: incline                 ! Orbital inclination

    ! Executable code
    do i = 2, nsons( root )
      son = subtree(i,root)
      key = subtree(1,son)

      select case ( decoration(key) )
      case ( f_module )
        instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_geodAngle )
        geodAngleNode = son
      case ( f_noMIFs )
        call expr ( subtree(2,son), expr_units, expr_value )
        noMIFs = expr_value(1)
      case ( f_inclination )
        call expr (subtree ( 2, son), expr_units, expr_value )
        incline = expr_value(1)
      case default ! Can't get here if tree_checker works correctly
      end select
    end do
    ! The tree checker will ensure we get all of these

    noMAFs = nsons(geodAngleNode) - 1

    call SetupNewQuantityTemplate ( mifGeolocation(instrumentModule), &
      & noInstances=noMAFs, noSurfs=noMIFs, noChans=1,&
      & coherent=.false., stacked=.false., regular=.true.,&
      & minorFrame=.true. )

    mifGeolocation(instrumentModule)%instrumentModule = instrumentModule
    mifGeolocation(instrumentModule)%noInstancesLowerOverlap = 0
    mifGeolocation(instrumentModule)%noInstancesUpperOverlap = 0
    mifGeolocation(instrumentModule)%verticalCoordinate = l_none
    mifGeolocation(instrumentModule)%time = 0.0
    mifGeolocation(instrumentModule)%solarTime = 0.0
    mifGeolocation(instrumentModule)%solarZenith = 0.0
    mifGeolocation(instrumentModule)%lon = 0.0
    mifGeolocation(instrumentModule)%losAngle = 0.0
    mifGeolocation(instrumentModule)%mafIndex = 0.0
    mifGeolocation(instrumentModule)%mafCounter = 0.0

    mifGeolocation(instrumentModule)%surfs = 0.0
    do maf = 1, noMAFs
      call expr (subtree ( maf+1, geodAngleNode), expr_units, expr_value )
      mifGeolocation(instrumentModule)%phi(:,maf) = expr_value(1)
      mifGeolocation(instrumentModule)%geodLat(:,maf) = &
        & mifGeolocation(instrumentModule)%phi(:,maf) ! Sort this out later
      mifGeolocation(instrumentModule)%mafIndex(maf) = maf + chunk%firstMAFIndex - 1
      mifGeolocation(instrumentModule)%mafCounter(maf) = maf + chunk%firstMAFIndex - 1
    end do

  end subroutine ForgeMinorFrames

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
! Revision 2.43  2001/07/16 18:24:45  dwu
! add feature for losTransFunc type of quantities
!
! Revision 2.42  2001/07/13 18:41:59  dwu
! fix problem after adding losTransFunc
!
! Revision 2.41  2001/07/13 18:10:03  dwu
! add quantity losTransFunc
!
! Revision 2.40  2001/07/11 21:40:00  livesey
! More bug fixes
!
! Revision 2.39  2001/07/10 23:45:16  jonathan
! added cloudicedensity and template for cloudsfwm, paul/jonathan
!
! Revision 2.38  2001/07/09 22:37:23  livesey
! Embarassing memory leak caught!  It's our old friend
! mifGeolocation again.  I'm going to regret trying
! to be this clever!
!
! Revision 2.37  2001/05/31 19:53:56  livesey
! Whoops, debug stuff left in.
!
! Revision 2.36  2001/05/30 23:59:51  livesey
! Thought I'd made this change already.  I'm confused
!
! Revision 2.35  2001/05/30 23:55:28  livesey
! Previous one was debug version, this is correct one.
!
! Revision 2.34  2001/05/30 23:53:01  livesey
! For new version of L1Bdata
!
! Revision 2.33  2001/05/29 23:21:35  livesey
! Changed l_orbitincline to l_orbitinclination
!
! Revision 2.32  2001/05/26 00:20:20  livesey
! Cosmetic changes
!
! Revision 2.31  2001/05/10 01:08:53  livesey
! Changed hGrids and vGrids to pointers, rather than intent(in)
! to allow them to be empty.
!
! Revision 2.30  2001/05/03 23:08:36  livesey
! Added stuff to support scan model items.
!
! Revision 2.29  2001/05/03 20:30:09  vsnyder
! Add a 'nullify' and some cosmetic changes
!
! Revision 2.28  2001/04/26 02:45:25  vsnyder
! Fix up CVS stuff
!
! Revision 2.27  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.26  2001/04/25 20:33:07  livesey
! Minor bug fix, Forge now also zeros surfs.
!
! Revision 2.25  2001/04/25 19:29:49  livesey
! Fixed bug in forge, now sets mafCounter and mafIndex correctly.
!
! Revision 2.24  2001/04/25 00:01:23  livesey
! Bug fix, no default units for scGeocAlt
!
! Revision 2.23  2001/04/24 22:21:17  livesey
! Gave up on latitude for forge
!
! Revision 2.22  2001/04/23 23:25:10  livesey
! Fixed bug in forge
!
! Revision 2.21  2001/04/20 23:11:39  livesey
! Added forge stuff for minor frames
!
! Revision 2.20  2001/04/19 20:30:06  livesey
! Added specific stuff for sideband ratio
!
! Revision 2.19  2001/04/12 23:25:29  vsnyder
! Give "Sideband" an initial value
!
! Revision 2.18  2001/04/12 21:41:42  livesey
! Signal now a string.
!
! Revision 2.17  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.16  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
!
! Revision 2.15  2001/03/28 23:48:13  livesey
! Bug fixes zero out some stuff.
!
! Revision 2.14  2001/03/28 18:33:19  livesey
! Fixed bug with logBasis (wasn't initialised!)
!
! Revision 2.13  2001/03/21 02:13:30  livesey
! Bug with logBasis, put in a work around. Will need to fix later
!
! Revision 2.12  2001/03/17 02:23:55  livesey
! Added logBasis (and set value for badData)
!
! Revision 2.11  2001/03/15 21:07:47  vsnyder
! Cross-references between databases are by database index, not tree index
!
! Revision 2.10  2001/03/15 18:41:17  livesey
! Tidied up the frequency coordinate stuff.
!
! Revision 2.9  2001/03/08 21:49:26  livesey
! Added elev_offset
!
! Revision 2.8  2001/03/03 00:08:09  livesey
! Lots of changes mostly with minor frame quantities
!
! Revision 2.7  2001/03/02 01:28:23  livesey
! New quantity types etc.
!
! Revision 2.6  2001/02/28 01:17:04  livesey
! Interim version, on the way to using proper signals stuff
!
! Revision 2.5  2001/02/22 23:37:24  livesey
! Really removed all references to firstIndexChannel
!
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
