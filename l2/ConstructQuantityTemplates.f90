! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ConstructQuantityTemplates

  ! This module is responsible for constructing templates for quantities.
  ! This version is a rewrite, aimed at tidying up a lot of the codebase
  use Init_tables_module, only: FIRST_LIT, LAST_LIT

  implicit none
  private

  public :: ConstructMinorFrameQuantity, CreateQtyTemplateFromMLSCFInfo
  public :: ForgeMinorFrames

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  ! The various properties has/can have
  integer, parameter :: NEXT = -1
  integer, parameter :: P_CHUNKED            = 1
  integer, parameter :: P_MAJORFRAME         = P_CHUNKED + 1
  integer, parameter :: P_MINORFRAME         = P_MAJORFRAME + 1
  integer, parameter :: P_MUSTBEZETA         = P_MINORFRAME + 1
  integer, parameter :: P_FGRID              = P_MUSTBEZETA + 1
  integer, parameter :: P_FGRIDOPTIONAL      = P_FGRID + 1
  integer, parameter :: P_HGRID              = P_FGRIDOPTIONAL + 1
  integer, parameter :: P_MODULE             = P_HGRID + 1
  integer, parameter :: P_MOLECULE           = P_MODULE + 1
  integer, parameter :: P_SGRID              = P_MOLECULE + 1
  integer, parameter :: P_VGRID              = P_SGRID + 1
  integer, parameter :: P_RADIOMETER         = P_VGRID + 1
  integer, parameter :: P_RADIOMETEROPTIONAL = P_RADIOMETER + 1
  integer, parameter :: P_REFLECTOR          = P_RADIOMETEROPTIONAL + 1
  integer, parameter :: P_SCMODULE           = P_REFLECTOR + 1
  integer, parameter :: P_SIGNAL             = P_SCMODULE + 1
  integer, parameter :: P_SUPPRESSCHANNELS   = P_SIGNAL + 1
  integer, parameter :: P_XYZ                = P_SUPPRESSCHANNELS + 1
  integer, parameter :: P_MATRIX3X3          = P_XYZ + 1

  integer, parameter :: NOPROPERTIES = P_MATRIX3X3

  ! Local, saved variables (constant tables really)
  logical, save, dimension ( noProperties, first_lit : last_lit ) :: &
    & PROPERTYTABLE
  integer, save, dimension ( first_lit : last_lit ) :: UNITSTABLE
  logical, save :: FIRSTCALL = .true.

contains ! ============= Public procedures ===================================

  ! ------------------------------------------- CreateQtyTemplateFromMLSCFInfo ----
  type (QuantityTemplate_T) function CreateQtyTemplateFromMLSCFInfo ( &
    & Name, Root, FGrids, VGrids, HGrids, L1bInfo, Chunk, MifGeolocation ) &
    & result ( QTY )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use EXPR_M, only: EXPR
    use FGrid, only: fGrid_T
    use HGrid, only: hGrid_T
    use INIT_TABLES_MODULE, only:  F_FGRID, F_GEODANGLE, F_HGRID, F_IRREGULAR, &
      & F_LOGBASIS, F_MINVALUE, F_MODULE, F_MOLECULE, F_RADIOMETER, F_SGRID, &
      & F_SIGNAL, F_TYPE, F_UNIT, F_VGRID, F_REFLECTOR, FIELD_FIRST, FIELD_LAST, &
      & L_TRUE, L_ZETA, L_XYZ, L_MATRIX3X3, L_CHANNEL, L_LOSTRANSFUNC, L_NONE
    use Intrinsic, only: LIT_INDICES
    use MLSCommon, only: L1BInfo_T, MLSChunk_T, RK => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only:GetModuleFromRadiometer, GetModuleFromSignal, &
      & GetRadiometerFromSignal, GetSignal, Signal_T, SIGNALS, MODULES, IsModuleSpacecraft
    use OUTPUT_M, only: OUTPUT
    use Parse_Signal_m, only: PARSE_SIGNAL
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate, &
      & NullifyQuantityTemplate
    use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
    use TOGGLES, only: GEN, TOGGLE, SWITCHES, LEVELS
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NODE_ID, NSONS, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_SET_ONE
    use VGridsDatabase, only: VGrid_T

    ! Dummy arguments
    integer, intent(in) :: NAME              ! Sub-rosa index of name
    integer, intent(in) :: ROOT              ! Root of QuantityTemplate subtree
    type (FGrid_T), dimension(:), pointer :: FGrids
    type (VGrid_T), dimension(:), pointer :: VGrids
    type (HGrid_T), dimension(:), pointer :: HGrids
    type (l1bInfo_T), intent(in) :: L1bInfo
    type (MLSChunk_T), intent(in) :: Chunk
    type (QuantityTemplate_T), dimension(:), intent(in), optional :: &
      & MifGeolocation

    ! Local variables
    logical :: LOGBASIS                 ! To place in quantity
    logical :: MAJORFRAME               ! Is a major frame quantity
    logical :: MINORFRAME               ! Is a minor frame quantity
    logical :: REGULAR                  ! Flag
    logical, dimension(noProperties) :: PROPERTIES ! Properties for this quantity type
    logical, dimension(field_first:field_last) :: GOT ! Fields
    character(len=127) :: SIGNALSTRING

    integer :: FGRIDINDEX               ! Index of frequency grid
    integer :: FREQUENCYCOORDINATE      ! Literal
    integer :: HGRIDINDEX               ! Index of horizontal grid
    integer :: I                        ! Loop counter
    integer :: INSTRUMENTMODULE         ! Database index
    integer :: KEY                      ! Field name, F_...
    integer :: MOLECULE                 ! Literal
    integer :: NOCHANS                  ! Quantity dimension
    integer :: NOINSTANCES              ! Quantity dimension
    integer :: NOSURFS                  ! Quantity dimension
    integer :: QUANTITYTYPE             ! Literal
    integer :: RADIOMETER               ! Database index
    integer :: REFLECTOR                ! Reflector literal
    integer :: SGRIDINDEX               ! Index for 'sGrid'
    integer :: SIDEBAND                 ! -1, 0, 1
    integer :: SIGNAL                   ! Database index
    integer :: SON                      ! A Son of Root -- an n_assign node
    integer :: VGRIDINDEX               ! Index in database
    integer :: VALUE                    ! Node index of value of field of spec
    integer :: S_INDEX                  ! Loop counter

    integer, dimension(2) :: EXPR_UNITS
    integer, dimension(:), pointer :: SignalInds ! From parse signal

    real(rk) :: MINVALUE                ! Minimum value allowed for quantity in fwm
    real(rk) :: SCALEFACTOR             ! Probably not used
    real(rk), dimension(2) :: EXPR_VALUE
    type (signal_T) :: SIGNALINFO       ! Details of the appropriate signal

    ! Executable code
    if ( firstCall ) then
      call InitQuantityTemplates
      firstCall = .true.
    end if

    ! Set appropriate defaults
    call nullifyQuantityTemplate ( qty ) ! for Sun's rubbish compiler
    nullify ( signalInds )
    fGridIndex = 0
    hGridIndex = 0
    instrumentModule = 0
    logBasis = .false.
    minValue = -huge ( 0.0_rk )
    molecule = 0
    noChans = 1
    quantitytype = 0
    radiometer = 0
    reflector = 0
    regular = .true.
    sGridIndex = 0
    scaleFactor = 1.0
    sideband = 0
    signal = 0
    vGridIndex = 0
    signalString = ''

    ! Go through the l2cf command line and parse it.
    got = .false.
    do i = 2, nsons(root)
      son = subtree(i,root)
      key = subtree(1,son)
      if ( node_id(son) == n_set_one ) then
        value = l_true
      else
        value = decoration(subtree(2,son))
      end if
      got ( decoration(key) ) = .true.

      select case ( decoration(key) )
      case ( f_fgrid )
        fGridIndex = decoration(value)
      case ( f_hgrid )
        hGridIndex = decoration(value)
      case ( f_logBasis )
        logBasis = (value == l_true)
      case ( f_irregular )
        regular = (value /= l_true)
      case ( f_minValue )
        call expr ( subtree(2,son), expr_units, expr_value )
        minValue = expr_value(1)
      case ( f_module)
        instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_molecule )
        molecule = value
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
        if ( size(signalInds) == 1 .or. .not. associated(L1bInfo%L1BRADIds) ) then
          signal = signalInds(1)
        else
          ! Seek a signal with any precision values !< 0
          do s_index=1, size(signalInds)
            if ( AnyGoodSignalData ( signalInds(s_index), sideband, &
              & l1bInfo, Chunk) ) exit
          end do
          if ( s_index > size(signalInds) ) then
            signal = signalInds(1)
          else
            signal = signalInds(s_index)
          end if
        end if
        call deallocate_test ( signalInds, 'signalInds', ModuleName )
        instrumentModule = GetModuleFromSignal(signal)
        radiometer = GetRadiometerFromSignal(signal)
      case ( f_reflector )
        reflector = value
      case ( f_sgrid )
        sGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_type )
        quantityType = value
      case ( f_vgrid )
        vGridIndex = decoration(value) ! node_id(value) == n_spec_args
      end select
    end do

    ! Do a very low level sanity check, irregular quantities not supported
    if ( .not. regular ) call announce_error ( root,&
      & 'Inappropriate irregular quantity request' )

    ! Now get the properties for this quantity type
    properties = propertyTable ( :, quantityType )

    ! Now check various things out first check that required fields are present
    ! according to the quantity type
    ! First those that are fairly clear cut to check
    if ( got ( f_hGrid ) .neqv. properties ( p_hGrid ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_hGrid) ) ) // ' hGrid for quantity type ', quantityType )
    if ( got ( f_vGrid ) .neqv. properties ( p_vGrid ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_vGrid) ) ) // ' vGrid for quantity type ', quantityType )
    if ( got ( f_sGrid ) .neqv. properties ( p_sGrid ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_sGrid) ) ) // ' sGrid for quantity type ', quantityType )
    if ( got ( f_Molecule ) .neqv. properties ( p_Molecule ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Molecule) ) ) // ' molecule for quantity type ', quantityType )
    if ( got ( f_Signal ) .neqv. properties ( p_Signal ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Signal) ) ) // ' signal for quantity type ', quantityType )
    if ( got ( f_Reflector ) .neqv. properties ( p_Reflector ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Reflector) ) ) // ' reflector for quantity type ', quantityType )

    ! These ones need a little more thought
    if ( .not. properties ( p_fGridOptional ) ) then
      if ( got ( f_fGrid ) .neqv. properties ( p_fGrid ) ) &
        & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
        & got(f_fGrid) ) ) // ' fGrid for quantity type ', quantityType )
    end if
    if ( .not. properties ( p_radiometerOptional ) ) then
      if ( got ( f_Radiometer ) .neqv. properties ( p_Radiometer ) ) &
        & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
        & got(f_Radiometer) ) ) // ' radiometer for quantity type ', quantityType )
    end if
    if ( got ( f_Module ) .neqv. ( properties ( p_Module ) .or. properties ( p_scModule) ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Module) ) ) // ' module for quantity type ', quantityType )

    ! Now do more derived checking / setting up
    if ( properties ( p_mustBeZeta ) ) then
      if ( vGrids(vGridIndex)%verticalCoordinate /= l_zeta ) &
        & call Announce_error ( root, 'Expecting log pressure coordinates for', &
        & quantityType )
    end if
    if ( properties ( p_scModule ) ) then
      if ( .not. IsModuleSpacecraft ( instrumentModule ) ) &
        & call Announce_error ( root, 'Module must be spacecraft' )
    end if

    ! Now establish the frequency coordinate system
    if ( got(f_fGrid) ) then
      frequencyCoordinate = fGrids(fGridIndex)%frequencyCoordinate
      noChans = fGrids(fGridIndex)%noChans
    else if ( properties(p_xyz) ) then
      ! XYZ quantity (e.g. ECI/ECR stuff)
      frequencyCoordinate = l_xyz
      noChans = 3
    else if ( properties(p_matrix3x3) ) then
      ! XYZ^2 quantity (e.g. rotation matrix)
      frequencyCoordinate = l_matrix3x3
      noChans = 9
    else if ( got(f_signal) .and. .not. properties(p_suppressChannels) ) then
      ! This is a channel based quantity
      signalInfo = GetSignal ( signal )
      frequencyCoordinate = l_channel
      noChans = size ( signalInfo%frequencies )
    else if ( got(f_sGrid) ) then
      ! Uses an sGrid
      noChans = size ( vGrids(sGridIndex)%surfs )
      frequencyCoordinate = l_losTransFunc
    else
      ! No frequency variation
      noChans = 1
      frequencyCoordinate = l_none
    end if

    ! Now do the setup for the different families of quantities
    if ( properties(p_minorFrame) ) then
      call ConstructMinorFrameQuantity ( l1bInfo, chunk, &
        & instrumentModule, qty, noChans=noChans, mifGeolocation=mifGeolocation, &
        & regular=regular )
      ! Setup a minor frame quantity
    else if ( properties(p_majorFrame) ) then
      ! Setup a major frame quantity
      call ConstructMajorFrameQuantity ( chunk, instrumentModule, &
        & qty, noChans, mifGeolocation )
    else
      ! Setup a non major/minor frame quantity
      noInstances = 1
      if ( properties ( p_hGrid ) ) noInstances = hGrids(hGridIndex)%noProfs
      noSurfs =1
      if ( properties ( p_vGrid ) ) noSurfs = vGrids(vGridIndex)%noSurfs

      ! Setup the quantity template
      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.true., stacked=.true., regular=.true.,&
        & noChans=noChans )

      ! Setup the horizontal coordinates or zero if irrelevant
      if ( properties(p_hGrid) ) then
        call CopyHGridInfoIntoQuantity ( hGrids(hGridIndex), qty )
      else                      ! Set `empty' values
        qty%phi = 0.0
        qty%geodLat = 0.0
        qty%lon = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
      end if

      ! Setup the vertical coordiantes or zero if irrelvant
      if ( properties(p_vGrid) ) then 
        call CopyVGridInfoIntoQuantity ( vGrids(vGridIndex), qty )
      else
        qty%surfs = 0.0
        qty%verticalCoordinate = l_none
      end if
    end if

    ! Fill the frequency information if appropriate
    if ( got ( f_fGrid ) ) then
      call Allocate_test ( qty%frequencies, qty%noChans, 'qty%frequencies', &
        & ModuleName )
      qty%frequencies = fGrids(fGridIndex)%values
    end if
    if ( properties ( p_sGrid ) ) then
      call Allocate_test ( qty%frequencies, noChans, 'qty%frequencies', ModuleName )
      qty%frequencies = VGrids(sGridIndex)%surfs
    end if

    ! Set up the remaining stuff
    qty%badValue = -999.99
    qty%frequencyCoordinate = frequencyCoordinate
    qty%instrumentmodule = instrumentmodule
    qty%logBasis = logBasis
    qty%minValue = minValue
    qty%molecule = molecule
    qty%name = name
    qty%quantityType = quantityType
    qty%radiometer = radiometer
    qty%reflector = reflector
    qty%scaleFactor = scaleFactor
    qty%sideband = sideband
    qty%signal = signal
    qty%unit = unitsTable ( quantityType )
    
    if ( index(switches,'qtmp') /= 0 ) then
      call output( 'Template name: ', advance='no' )
      call display_string ( qty%name, advance='no' )
      call output( '   quantityType: ', advance='no' )
      call display_string( lit_indices(qty%quantityType), advance='yes<' )
      call output( '   major frame? ', advance='no' )
      call output( qty%majorFrame, advance='no' )
      call output( '   minor frame? ', advance='no' )
      call output( qty%minorFrame, advance='yes' )
      call output( '    signal name: ', advance='no' )
      call output( trim(signalString), advance='yes' )
      call output( '   Instrument module: ', advance='no' )
      if ( qty%instrumentModule < 1 ) then
        call output( ' (none) ', advance='yes' )
      else
        call display_string ( modules(qty%instrumentModule)%name, advance='yes' )
      end if
      call output( '    Molecule: ', advance='no' )
      if  ( qty%molecule < 1 ) then
        call output ( ' (none)', advance='yes' )
      else
        call display_string ( qty%Molecule, advance='yes' )
      end if
      call output ( '   noChans = ' )
      call output ( qty%noChans, advance='no' )
      call output ( ' noSurfs = ' )
      call output ( qty%noSurfs, advance='no' )
      call output ( ' noInstances = ' )
      call output ( qty%noInstances, advance='no' )
      call output ( ' instanceLen = ' )
      call output ( qty%instanceLen, advance='yes' )
      call output ( ' verticalCoordinate = ' )
      call display_string ( lit_indices(qty%verticalCoordinate) )
      call output ( ' frequencyCoordinate = ' )
      call display_string ( lit_indices(qty%frequencyCoordinate), advance='yes' )
    end if

  end function CreateQtyTemplateFromMLSCFInfo

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  subroutine ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
    & qty, noChans, regular, instanceLen, mifGeolocation )

    use INIT_TABLES_MODULE, only: L_GEODALTITUDE, L_NONE
    use L1BData, only: L1BData_T, READL1BDATA, DEALLOCATEL1BDATA, &
      & AssembleL1BQtyName
    use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, RK => R8
    use MLSFiles, only: MLS_HDF_Version
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_L1BRead
    use MLSSignals_m, only:  IsModuleSpacecraft, GetModuleName
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate, &
      & QuantityTemplateCounter

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
    real(rk), parameter :: SIXTH = 1.0_rk / 6.0_rk
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
    integer :: noMAFs, l1bFlag, l1bItem, mafIndex, mifIndex, hdfVersion

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    hdfVersion = mls_hdf_version(trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION)
    if ( hdfversion <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )

    if ( present(mifGeolocation) ) then
      ! -------------------------------------- Got mifGeolocation ------------
      if ( .not. (present(noChans)) ) &
         call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'You must supply NoChans to reuse geolocation information' )
      ! We have geolocation information, setup the quantity as a clone of that.
      qty = mifGeolocation(instrumentModule)
      if ( present ( regular ) ) qty%regular = regular
      qty%noChans = noChans
      if ( qty%regular ) then
        qty%instanceLen = qty%noChans*qty%noSurfs
      else
        qty%instanceLen = 0
      end if
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
      end if
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )

      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & hdfVersion=hdfVersion )
      if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
      
      call SetupNewQuantityTemplate ( qty, noInstances=noMAFs, &
        & noSurfs=l1bField%maxMIFs, noChans=noChans, coherent=.false., &
        & stacked=.false., regular=regular, instanceLen=instanceLen, &
        & minorFrame=.true. )
      qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
      qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap

      if ( .not. IsModuleSpacecraft(instrumentModule) ) then
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
        
        l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
        call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
          & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & hdfVersion=hdfVersion )
        if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        
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
          l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
          call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
            & l1bFlag, firstMAF=chunk%firstMafIndex, &
            & lastMAF=chunk%lastMafIndex, hdfVersion=hdfVersion )
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

    use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_FGRID, F_GEODANGLE, F_MODULE, F_NOMIFS, &
      & F_SOLARTIME, F_SOLARZENITH
    use INIT_TABLES_MODULE, only: L_NONE, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, &
      & PHYQ_TIME
    use MLSCommon, only: MLSChunk_T, RK => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate
    use TREE, only: DECORATION, NSONS, SUBTREE

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
    integer :: NODE                     ! A tree node
    integer :: NOMAFS                   ! Dimension
    integer :: NOMIFS                   ! Dimension
    integer :: PARAM                    ! Loop counter
    integer :: SOLARTIMENODE            ! Tree vertex
    integer :: SOLARZENITHNODE          ! Tree vertex
    integer :: SON                      ! Tree vertex
    integer :: UNITS                    ! Units for node

    real(rk), dimension(:,:), pointer :: VALUES ! An array to fill
    integer, dimension(2) :: EXPR_UNITS ! From tree
    real(rk), dimension(2) :: EXPR_VALUE ! From tree

    ! Executable code
    solarTimeNode = 0
    solarZenithNode = 0
    geodAngleNode = 0
    do i = 2, nsons( root )
      son = subtree(i,root)
      key = subtree(1,son)

      select case ( decoration(key) )
      case ( f_module )
        instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_geodAngle )
        geodAngleNode = son
      case ( f_solarTime )
        solarTimeNode = son
      case ( f_solarZenith )
        solarZenithNode = son
      case ( f_noMIFs )
        call expr ( subtree(2,son), expr_units, expr_value )
        noMIFs = expr_value(1)
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    ! Work out how many MAFs we're after
    if ( geodAngleNode /= 0 ) noMAFs = nsons ( geodAngleNode ) - 1
    if ( solarTimeNode /= 0 ) then
      if ( noMAFs /= 0 ) then
        if ( nsons ( solarTimeNode ) - 1 /= noMAFs ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Inconsistent explicit Forge definition' )
      else
        noMAFs = nsons ( solarTimeNode ) - 1
      end if
    end if
    if ( solarZenithNode /= 0 ) then
      if ( noMAFs /= 0 ) then
        if ( nsons ( solarZenithNode ) - 1 /= noMAFs ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Inconsistent explicit Forge definition' )
      else
        noMAFs = nsons ( solarZenithNode ) - 1
      end if
    end if

    ! Setup the minor frame quantity template
    call SetupNewQuantityTemplate ( mifGeolocation(instrumentModule), &
      & noInstances=noMAFs, noSurfs=noMIFs, noChans=1,&
      & coherent=.false., stacked=.false., regular=.true.,&
      & minorFrame=.true. )

    ! Put zeros etc. in the appropriate places, may overwrite these later
    mifGeolocation(instrumentModule)%instrumentModule = instrumentModule
    mifGeolocation(instrumentModule)%noInstancesLowerOverlap = 0
    mifGeolocation(instrumentModule)%noInstancesUpperOverlap = 0
    mifGeolocation(instrumentModule)%verticalCoordinate = l_none
    mifGeolocation(instrumentModule)%time = 0.0
    mifGeolocation(instrumentModule)%phi = 0.0
    mifGeolocation(instrumentModule)%geodLat = 0.0
    mifGeolocation(instrumentModule)%solarTime = 0.0
    mifGeolocation(instrumentModule)%solarZenith = 0.0
    mifGeolocation(instrumentModule)%lon = 0.0
    mifGeolocation(instrumentModule)%losAngle = 0.0
    mifGeolocation(instrumentModule)%mafIndex = 0.0
    mifGeolocation(instrumentModule)%mafCounter = 0.0
    mifGeolocation(instrumentModule)%surfs = 0.0

    ! Loop over the parameters we might have
    do param = 1, 3
      select case ( param )
      case ( 1 )
        values => mifGeolocation(instrumentModule)%phi
        node = geodAngleNode
        units = phyq_angle
      case ( 2 )
        values => mifGeolocation(instrumentModule)%solarTime
        node = solarTimeNode
        units = phyq_time
      case ( 3 )
        values => mifGeolocation(instrumentModule)%solarZenith
        node = solarZenithNode
        units = phyq_angle
      end select

      if ( node /= 0 ) then
        do maf = 1, noMAFs
          call expr ( subtree ( maf+1, node), expr_units, expr_value )
          if ( all ( expr_units(1) /= (/ phyq_dimensionless, units /) ) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Invalid units for explicit hGrid' )
          values(:,maf) = expr_value(1)
        end do
      end if
    end do
    ! Make the latitude the same as the geod angle
    ! This is a bit of a hack, but it should be OK in all cases.
    mifGeolocation(instrumentModule)%geodLat = &
      & mifGeolocation(instrumentModule)%phi
    
  end subroutine ForgeMinorFrames

  ! ======================================= Private proceedures ========
  
  ! -----------------------------------------------  Announce_Error  -----
  subroutine Announce_Error ( where, message, extra )

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: SOURCE_REF
    use Intrinsic, only: LIT_INDICES
    use String_Table, only: DISPLAY_STRING
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    call output ( message )
    if ( present ( extra ) ) call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem in Construct' )
  end subroutine Announce_Error

  ! -----------------------------------------------  AnyGoodSignalData  -----
  function AnyGoodSignalData ( signal, sideband, l1bInfo, Chunk )  result (answer)
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  ! Arguments

    use Allocate_Deallocate, only: Deallocate_Test
    use L1BData, only: L1BData_T, READL1BDATA, &
      & FindL1BData, AssembleL1BQtyName, PRECISIONSUFFIX
    use MLSCommon, only: L1BInfo_T, MLSChunk_T, RK => R8
    use MLSFiles, only: MLS_HDF_Version
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSSignals_m, only: GetSignalName

    integer, intent(in) :: signal
    integer, intent(in) :: sideband
    logical             :: answer
    type (MLSChunk_T), intent(in) :: Chunk
    type (l1bInfo_T), intent(in) :: L1bInfo
  ! Private
    integer :: FileID, flag, noMAFs
    character(len=127)  :: namestring
    type (l1bData_T) :: MY_L1BDATA
    integer :: hdfVersion

  ! Executable
    hdfVersion = mls_hdf_version(trim(l1bInfo%L1BOAFileName), LEVEL1_HDFVERSION)
    if ( hdfversion <= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Illegal hdf version for l1boa file (file missing or non-hdf?)' )
    call GetSignalName ( signal, nameString, &                   
    & sideband=sideband, noChannels=.TRUE. )                     
    nameString = AssembleL1BQtyName ( nameString, hdfVersion, .false. )
    nameString = trim(nameString) // PRECISIONSUFFIX
    fileID = FindL1BData (l1bInfo%l1bRadIDs, nameString, hdfVersion )
    if ( fileID <= 0 ) then
      answer = .false.
      return
    end if
    ! print *, 'About to read ', trim(nameString)
    ! print *, 'From Fileid ', fileID
    call ReadL1BData ( fileID , nameString, my_l1bData, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & NeverFail= .true., hdfVersion=hdfVersion )
    if ( flag == 0 ) then
      answer = .not. all (my_l1bData%DpField < 0._rk)
      call deallocate_test(my_l1bData%DpField, trim(nameString), ModuleName)
    else
      answer = .false.
    end if
  end function AnyGoodSignalData

  ! ----------------------------------  CopyHGridInfoIntoQuantity  -----
  subroutine CopyHGridInfoIntoQuantity ( My_hGrid, qty )

  ! This routine copies HGrid information into an already defined quantity

    use HGrid, only: hGrid_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_T

    ! Dummy arguments
    type (hGrid_T), intent(in) :: My_hGrid
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code

    if ( qty%noInstances/=my_hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )

    if ( .NOT. qty%stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy hGrids into unstacked quantities")

    qty%phi(1,:) = my_hGrid%phi
    qty%geodLat(1,:) = my_hGrid%geodLat
    qty%lon(1,:) = my_hGrid%lon
    qty%time(1,:) = my_hGrid%time
    qty%solarTime(1,:) = my_hGrid%solarTime
    qty%solarZenith(1,:) = my_hGrid%solarZenith
    qty%losAngle(1,:) = my_hGrid%losAngle
    qty%noInstancesLowerOverlap = my_hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = my_hGrid%noProfsUpperOverlap

  end subroutine CopyHGridInfoIntoQuantity

  ! ----------------------------------  CopyVGridInfoIntoQuantity  -----
  subroutine CopyVGridInfoIntoQuantity ( vGrid, qty )

  ! This similar routine copies VGrid information into an already 
  ! defined quantity

    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_T
    use VGridsDatabase, only: VGrid_T

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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

  ! --------------------------------  ConstructMajorFrameQuantity  -----
  subroutine ConstructMajorFrameQuantity( chunk, instrumentModule, qty, noChans, &
    & mifGeolocation )
    ! Dummy arguments
    use INIT_TABLES_MODULE, only: L_NONE
    use MLSCommon, only: MLSChunk_T, RK => R8
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate

    type (MLSChunk_T), intent(in) :: CHUNK
    integer, intent(in) :: INSTRUMENTMODULE
    type (QuantityTemplate_T), intent(out) :: QTY
    type (QuantityTemplate_T), dimension(:), target :: MIFGEOLOCATION
    integer, intent(in) :: NOCHANS
    ! Local variables
    type (QuantityTemplate_T), pointer :: source
    
    ! Executable code
    source => mifGeolocation(instrumentModule)
    call SetupNewQuantityTemplate ( qty, noInstances=source%noInstances, &
      & noSurfs=1, coherent=.true., stacked=.true., regular=.true., &
      & noChans=noChans )
    qty%phi => source%phi(1:1,:)
    qty%geodLat => source%geodLat(1:1,:)
    qty%lon => source%lon(1:1,:)
    qty%time => source%time(1:1,:)
    qty%solarTime => source%solarTime(1:1,:)
    qty%solarZenith => source%solarZenith(1:1,:)
    qty%losAngle => source%losAngle(1:1,:)
    qty%mafIndex => source%mafIndex
    qty%mafCounter => source%mafCounter

    qty%verticalCoordinate = l_none
    qty%majorFrame = .true.
    qty%minorFrame = .false.

    qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
    qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap
  end subroutine ConstructMajorFrameQuantity
  
  ! ----------------------------------------------- InitQuantityTemplates ----
  subroutine InitQuantityTemplates
    ! This routine initializes the quantity template properties
    ! This is the routine one needs to update when one introduces a new quantity type.
    use Init_Tables_Module, only:  L_BASELINE, L_BOUNDARYPRESSURE, &
      L_CALSIDEBANDFRACTION, &
      L_CHISQBINNED, L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, L_CLOUDICE, &
      L_CLOUDINDUCEDRADIANCE, L_CLOUDEXTINCTION, L_CLOUDRADSENSITIVITY, &
      L_CLOUDWATER, L_COLUMNABUNDANCE, &
      L_DNWT_AJN, L_DNWT_AXMAX, L_DNWT_CAIT, &
      L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, &
      L_DNWT_DIAG, L_DNWT_DXDX, L_DNWT_DXDXL, &
      L_DNWT_DXN, L_DNWT_DXNL, L_DNWT_FLAG, L_DNWT_FNMIN, &
      L_DNWT_FNORM, L_DNWT_GDX, L_DNWT_GFAC, &
      L_DNWT_GRADN, L_DNWT_SQ, L_DNWT_SQT,&
      L_EARTHREFL, L_ECRTOFOV, L_EFFECTIVEOPTICALDEPTH, &
      L_ELEVOFFSET, L_EXTINCTION, L_GPH, L_HEIGHTOFFSET, &
      L_ISOTOPERATIO, L_JACOBIAN_COLS, L_JACOBIAN_ROWS, &
      L_LIMBSIDEBANDFRACTION, L_LOSTRANSFUNC, L_LOSVEL, &
      L_MASSMEANDIAMETERICE, L_MASSMEANDIAMETERWATER, L_MAGNETICFIELD, &
      L_NOISEBANDWIDTH, L_NORADSPERMIF, L_NORADSBINNED, &
      L_NUMJ, L_OPTICALDEPTH, &
      L_ORBITINCLINATION, L_PHITAN, L_PTAN, L_RADIANCE, L_EARTHRADIUS,&
      L_REFGPH, L_REFLTEMP, L_REFLTRANS, L_REFLREFL, L_REFLSPILL, &
      L_RHI, L_SIZEDISTRIBUTION, &
      L_SCANRESIDUAL, L_SCECI, L_SCVEL, L_SCVELECI, &
      L_SCVELECR, L_SCGEOCALT, &
      L_SPACERADIANCE, L_STRAYRADIANCE, L_SURFACETYPE, L_SYSTEMTEMPERATURE, &
      L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, L_TNGTGEOCALT, &
      L_TOTALEXTINCTION, L_VMR
    use Init_Tables_Module, only: PHYQ_EXTINCTION, PHYQ_FREQUENCY,&
      & PHYQ_GAUSS, PHYQ_IceDensity, PHYQ_LENGTH, &
      & PHYQ_PCTRHI, PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_VELOCITY, &
      & PHYQ_VMR, PHYQ_ZETA, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, PHYQ_DOBSONUNITS
    use MLSMessageModule, only: MLSMSG_Error, MLSMessage
    use Intrinsic, only: LIT_INDICES
    use Output_M, only: OUTPUT
    use String_Table, only: DISPLAY_STRING

    ! Local variables
    integer :: I                        ! Loop counter
    integer, dimension(0), parameter :: NONE = (/ ( 0, i=1,0 ) /)
    logical :: VALID                    ! Flag
    character(len=132) :: MESSAGE       ! An error message

    ! Executable code
    ! Basically here, we're going to go through and populate the various tables

    propertyTable = .false.
    unitsTable = 0

    call DefineQtyTypes ( (/ &
      l_baseline, phyq_temperature, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                  p_mustBeZeta, next, &
      l_boundaryPressure, phyq_pressure, p_hGrid, next, &
      l_calSidebandFraction, phyq_dimensionless, p_signal, next, &
      l_chisqBinned, phyq_dimensionless, p_hGrid, p_vGrid, &
                     p_signal, p_suppressChannels, p_mustBeZeta, next, &
      l_chisqChan, phyq_dimensionless, p_majorFrame, p_signal, next, &
      l_chisqMMAF, phyq_dimensionless, p_majorFrame, p_signal, &
                   p_suppressChannels, next, &
      l_chisqMMIF, phyq_dimensionless, p_minorFrame, p_signal, &
                   p_suppressChannels, next, &
      l_cloudExtinction, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_cloudIce, phyq_IceDensity, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_cloudInducedRadiance, phyq_temperature, p_minorFrame, p_signal, next, &
      l_cloudRadSensitivity, phyq_temperature, p_minorFrame, p_signal, next, &
      l_cloudWater, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_columnAbundance, phyq_dobsonunits, p_hGrid, p_molecule, next, &
      l_dnwt_ajn, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_axmax, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_cait, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_chisqminnorm, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_chisqnorm, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_diag, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_dxdx, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_dxdxl, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_dxn, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_dxnl, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_flag, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_fnmin, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_fnorm, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_gdx, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_gfac, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_gradn, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_sq, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_dnwt_sqt, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_earthRadius, phyq_length, p_hGrid, next, &
      l_earthRefl, phyq_dimensionless, none /) )

    call DefineQtyTypes ( (/ &
      l_ecrToFOV, phyq_dimensionless, p_minorFrame, p_module, p_matrix3x3, next, &
      l_effectiveOpticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_elevOffset, phyq_angle, p_signal, next, &
      l_extinction, phyq_extinction, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                    p_mustBeZeta, next, &
      l_gph, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_heightOffset, phyq_length, p_hGrid, p_vGrid, next, &
      l_isotopeRatio, phyq_dimensionless, p_molecule, next, &
      l_jacobian_cols, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_jacobian_rows, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_limbSidebandFraction, phyq_dimensionless, p_signal, next, &
      l_losTransFunc, phyq_dimensionless, p_minorFrame, p_sGrid, p_module, next, &
      l_losVel, phyq_dimensionless, p_minorFrame, p_module, next, &
      l_magneticField, phyq_gauss, p_vGrid, p_hGrid, p_xyz, p_mustBeZeta, next, &
      l_massMeanDiameterIce, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_massMeanDiameterWater, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_noRadsBinned, phyq_dimensionless, p_vGrid, p_hGrid, &
                      p_signal, p_suppressChannels, p_mustBeZeta, next, &
      l_noRadsPerMIF, phyq_dimensionless, p_minorFrame, p_signal, &
                      p_suppressChannels, next, &
      l_noiseBandwidth, phyq_frequency, p_signal, next, &
      l_numJ, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, &
      l_opticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_orbitInclination, phyq_angle, p_minorFrame, p_scModule, next, &
      l_phiTan, phyq_angle, p_minorFrame, p_module, next, & 
      l_ptan, phyq_zeta, p_minorFrame, p_module, next /) )

    call DefineQtyTypes ( (/ &
      l_radiance, phyq_temperature, p_minorFrame, p_signal, next, & 
      l_refltemp, phyq_temperature, p_majorFrame, p_reflector, p_module, next, &
      l_refltrans, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflrefl, phyq_dimensionless, p_reflector, next, &
      l_reflspill, phyq_temperature, p_signal, p_majorframe, p_reflector, next, &
      l_refGPH, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_rhi, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, p_mustBeZeta, next, &
      l_scECI, phyq_length, p_minorFrame, p_scModule, p_xyz, next, &
      l_scGeocAlt, phyq_length, p_minorFrame, p_scModule, next, &
      l_scVel, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECI, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECR, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scanResidual, phyq_length, p_minorFrame, p_module, next, &
      l_sizeDistribution, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, & 
      l_spaceRadiance, phyq_temperature, none, next, &
      l_strayRadiance, phyq_temperature, p_signal, p_majorFrame, next, &
      l_surfaceType, phyq_dimensionless, p_hGrid, next, & 
      l_systemTemperature, phyq_temperature, p_signal, next, &
      l_temperature, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, next, &
      l_tngtECI, phyq_length, p_minorFrame, p_module, p_xyz, next, &
      l_tngtGeocAlt, phyq_length, p_minorFrame, p_module, next, &
      l_tngtGeodAlt, phyq_length, p_minorFrame, p_module, next, &
      l_vmr, phyq_vmr, p_hGrid, p_vGrid, p_fGridOptional, p_molecule, &
             p_radiometerOptional, p_mustbezeta, next /) )

    ! Do a bit of checking
    do i = first_lit, last_lit
      valid = .true.
      message =  ''
      ! Check it's not both major and minor frame
      if ( count ( propertyTable ( (/ p_minorFrame, p_majorFrame /), i ) ) > 1 ) then
        valid = .false.
        message = "Quantity cannot be both major and minor frame"
      end if
      ! Check that we can identify the module for major/minor frame
      if ( ( propertyTable ( p_minorFrame, i ) .or. &
        &    propertyTable ( p_majorFrame, i ) ) .and. &
        & .not. any ( propertyTable ( &
        & (/ p_radiometer, p_module, p_scModule, p_signal /), i ) ) ) then
        valid = .false.
        message = "Badly defined major/minor frame quantity"
      end if
      ! Check that mustBeZeta quantities have a vGrid
      if ( propertyTable ( p_mustBeZeta, i ) .and. .not. propertyTable ( p_vGrid, i ) ) then
        valid = .false.
        message = "Quantity must have vGrid if it must be on log-pressure"
      end if
      ! Print out any message
      if ( .not. valid ) then
        call output ( "Offending quantity: " )
        call display_string ( lit_indices ( i ), strip=.true., advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, message )
      end if
    end do
  contains
    ! --------------------------- Internal subroutine
    subroutine DefineQtyTypes ( info )
      integer, dimension(:), intent(in) :: INFO
      ! Local variables
      integer :: I                      ! Location
      integer :: QTYTYPE                ! Index
      ! Executable code
      qtyType = 0
      i = 1
      defineLoop: do
        if ( i > size(info) ) exit defineLoop
        if ( qtyType == 0 ) then
          qtyType = info ( i )
          propertyTable ( :, qtyType ) = .false.
          i = i + 1
          if ( i > size(info) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Malformed call to DefineQtyTypes' )
          unitsTable ( qtyType ) = info ( i )
        else
          if ( info(i) /= next ) then
            propertyTable ( info(i), qtyType ) = .true.
          else
            qtyType = 0
          end if
        end if
        i = i + 1
      end do defineLoop
    end subroutine DefineQtyTypes
    
  end subroutine InitQuantityTemplates

end module ConstructQuantityTemplates
