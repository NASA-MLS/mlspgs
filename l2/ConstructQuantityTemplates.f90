! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
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
  integer, parameter :: P_FLEXIBLEVHGRID     = P_FGRIDOPTIONAL + 1
  integer, parameter :: P_HGRID              = P_FLEXIBLEVHGRID + 1
  integer, parameter :: P_MODULE             = P_HGRID + 1
  integer, parameter :: P_MOLECULE           = P_MODULE + 1
  integer, parameter :: P_SGRID              = P_MOLECULE + 1
  integer, parameter :: P_VGRID              = P_SGRID + 1
  integer, parameter :: P_RADIOMETER         = P_VGRID + 1
  integer, parameter :: P_RADIOMETEROPTIONAL = P_RADIOMETER + 1
  integer, parameter :: P_REFLECTOR          = P_RADIOMETEROPTIONAL + 1
  integer, parameter :: P_SCMODULE           = P_REFLECTOR + 1
  integer, parameter :: P_SIGNAL             = P_SCMODULE + 1
  integer, parameter :: P_SIGNALOPTIONAL     = P_SIGNAL + 1
  integer, parameter :: P_SUPPRESSCHANNELS   = P_SIGNALOPTIONAL + 1
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
    use Chunks_m, only: MLSChunk_T
    use EXPR_M, only: EXPR
    use FGrid, only: fGrid_T
    use HGridsDatabase, only: hGrid_T
    use INIT_TABLES_MODULE, only:  F_FGRID, F_HGRID, F_IRREGULAR, &
      & F_LOGBASIS, F_MINVALUE, F_MODULE, F_MOLECULE, F_RADIOMETER, F_SGRID, &
      & F_SIGNAL, F_TYPE, F_VGRID, F_REFLECTOR, FIELD_FIRST, FIELD_LAST, &
      & L_TRUE, L_ZETA, L_XYZ, L_MATRIX3X3, L_CHANNEL, L_LOSTRANSFUNC, L_NONE
    use Intrinsic, only: LIT_INDICES
    use MLSCommon, only: L1BInfo_T, RK => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSSignals_m, only:GetModuleFromRadiometer, GetModuleFromSignal, &
      & GetRadiometerFromSignal, GetSignal, Signal_T, MODULES, IsModuleSpacecraft
    use OUTPUT_M, only: OUTPUT
    use Parse_Signal_m, only: PARSE_SIGNAL
    use QuantityTemplates, only: QuantityTemplate_T, &
      & DUMP, SetupNewQuantityTemplate, &
      & NullifyQuantityTemplate
    use STRING_TABLE, only: GET_STRING, DISPLAY_STRING
    use TOGGLES, only: SWITCHES
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
    logical :: ISMINORFRAME             ! Is a minor frame quantity
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
    real(rk), dimension(2) :: EXPR_VALUE
    type (signal_T) :: SIGNALINFO       ! Details of the appropriate signal

    ! Executable code
    if ( firstCall ) then
      call InitQuantityTemplates
      firstCall = .true.
    end if

    ! Set appropriate defaults
    call NullifyQuantityTemplate ( qty ) ! for Sun's rubbish compiler
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
    if ( .not. properties ( p_flexibleVHGrid ) ) then
      if ( got ( f_hGrid ) .neqv. properties ( p_hGrid ) ) &
        & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
        & got(f_hGrid) ) ) // ' hGrid for quantity type ', quantityType, &
        & severity='nonfatal' )
      if ( got ( f_vGrid ) .neqv. properties ( p_vGrid ) ) &
        & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
        & got(f_vGrid) ) ) // ' vGrid for quantity type ', quantityType )
    endif
    if ( got ( f_sGrid ) .neqv. properties ( p_sGrid ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_sGrid) ) ) // ' sGrid for quantity type ', quantityType )
    if ( got ( f_Molecule ) .neqv. properties ( p_Molecule ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Molecule) ) ) // ' molecule for quantity type ', quantityType )
    if ( got ( f_Reflector ) .neqv. properties ( p_Reflector ) ) &
      & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
      & got(f_Reflector) ) ) // ' reflector for quantity type ', quantityType )

    ! These ones need a little more thought
    if ( .not. properties ( p_signalOptional ) ) then
      if ( got ( f_Signal ) .neqv. properties ( p_Signal ) ) &
        & call Announce_error ( root, trim ( merge ( 'unexpected', 'need      ', &
        & got(f_Signal) ) ) // ' signal for quantity type ', quantityType )
    end if
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
    if ( properties ( p_mustBeZeta ) .and. got(f_vGrid) ) then
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
      qty%frequencies => fGrids(fGridIndex)%values
      noChans = fGrids(fGridIndex)%noChans
      ! print *, '1st try: fGridIndex ', fGridIndex
      ! print *, fGrids(fGridIndex)%values
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

    ! Now deal with the fleixbleVHGrid quantities
    if ( properties(p_flexibleVHGrid) ) then
      if ( got ( f_hGrid ) .neqv. got ( f_vGrid ) ) &
        & call Announce_Error ( root, 'Must supply both or neither vGrid and hGrid' )
      isMinorFrame = .not. got ( f_hGrid )
    else
      isMinorFrame = properties(p_minorFrame)
    endif

    ! Now do the setup for the different families of quantities
    if ( isMinorFrame ) then
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
      if ( got ( f_hGrid ) ) noInstances = hGrids(hGridIndex)%noProfs
      noSurfs = 1
      if ( got ( f_vGrid ) ) noSurfs = vGrids(vGridIndex)%noSurfs

      ! Setup the quantity template
      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.true., stacked=.true., regular=.true.,&
        & noChans=noChans, &
        & sharedVGrid=.true., sharedHGrid=.true., sharedFGrid=.true. )

      ! Setup the horizontal coordinates
      if ( got(f_hGrid) ) then
        qty%hGridIndex = hGridIndex
        call PointQuantityToHGrid ( hGrids(hGridIndex), qty )
      else
        call SetupEmptyHGridForQuantity ( qty )
      end if
      ! Work out the instance offset
      if ( associated ( chunk%hGridOffsets ) ) then
        if ( got ( f_hGrid )  ) then
          qty%instanceOffset = chunk%hGridOffsets(hGridIndex)
          qty%grandTotalInstances = chunk%hGridTotals(hGridIndex)
        else
          ! Must have a single instance per chunk
          qty%instanceOffset = chunk%chunkNumber - 1
          ! -1 because it's an offset remember, not an origin.
          qty%grandTotalInstances = 0
        end if
      end if

      ! Setup the vertical coordiantes
      if ( got (f_vGrid) ) then
        qty%vGridIndex = vGridIndex
        qty%verticalCoordinate = vGrids(vGridIndex)%verticalCoordinate
        qty%surfs => vGrids(vGridIndex)%surfs
      else
        call SetupEmptyVGridForQuantity ( qty )
      end if
    end if

    ! Fill the frequency information if appropriate
    if ( got ( f_fGrid ) ) then
      ! print *, 'About to: qty%GridIndex was ', qty%fGridIndex
      ! print *, 'qty%sharedFGrid ', qty%sharedFGrid
      ! print *, fGrids(fGridIndex)%values
      qty%fGridIndex = fGridIndex
      qty%frequencies => fGrids(fGridIndex)%values
      if ( .not. qty%sharedFGrid ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName,&
          & 'sharedFGrid needed to be reset' )
        qty%sharedFGrid = .true.
      endif
      ! print *, '2nd try: fGridIndex ', fGridIndex
      ! print *, 'qty%sharedFGrid ', qty%sharedFGrid
      ! print *, fGrids(fGridIndex)%values
    end if
    if ( properties ( p_sGrid ) ) then
      qty%sharedFGrid = .false.
      nullify(qty%frequencies) ! Lest we deallocate a database entry
      call Allocate_test ( qty%frequencies, qty%noChans, 'qty%frequencies', &
        & ModuleName )
      qty%fGridIndex = -sGridIndex        ! Use -ve value to denote sGrid
      qty%frequencies = vGrids(sGridIndex)%surfs(:,1)
      ! print *, 'Resetting: fGridIndex ', qty%fGridIndex
      ! print *, qty%frequencies
    end if

    ! Set up the remaining stuff
    qty%frequencyCoordinate = frequencyCoordinate
    qty%instrumentmodule = instrumentmodule
    qty%logBasis = logBasis
    qty%minValue = minValue
    qty%molecule = molecule
    qty%name = name
    qty%quantityType = quantityType
    qty%radiometer = radiometer
    qty%reflector = reflector
    qty%sideband = sideband
    qty%signal = signal
    qty%unit = unitsTable ( quantityType )
    
    if ( index(switches,'qtmp') /= 0 ) then
      call output( 'Template name: ', advance='no' )
      call display_string ( qty%name, advance='no' )
      call output( '   quantityType: ', advance='no' )
      call display_string( lit_indices(qty%quantityType), advance='yes' )
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

    if ( index(switches, 'qtmp') > 0 ) call dump(qty, details=0, noL2CF=.true.)
  end function CreateQtyTemplateFromMLSCFInfo

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  subroutine ConstructMinorFrameQuantity ( l1bInfo, chunk, instrumentModule, &
    & qty, noChans, regular, instanceLen, mifGeolocation )

    use Chunks_m, only: MLSChunk_T
    use INIT_TABLES_MODULE, only: L_GEODALTITUDE, L_NONE
    use L1BData, only: L1BData_T, READL1BDATA, DEALLOCATEL1BDATA, &
      & AssembleL1BQtyName
    use MLSCommon, only: L1BInfo_T, NameLen, RK => R8
    use MLSFiles, only: MLS_HDF_Version
    use MLSL2Options, only: LEVEL1_HDFVERSION
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_L1BRead
    use MLSSignals_m, only:  IsModuleSpacecraft, GetModuleName
    use Output_m, only: OUTPUT
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate
    use TOGGLES, only: SWITCHES

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
      qty%sharedHGrid = .true.
      qty%sharedVGrid = .true.
      if ( present ( regular ) ) qty%regular = regular
      qty%noChans = noChans
      if ( qty%regular ) then
        qty%instanceLen = qty%noChans*qty%noSurfs
      else
        qty%instanceLen = 0
      end if

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

      ! Now we're going to fill in the hGrid information
      qty%instanceOffset = chunk%firstMAFIndex + chunk%noMAFsLowerOverlap
      qty%grandTotalInstances = 0
      if ( index(switches,'qtmp') /= 0 ) then
        call output ( "Instance offset for minor frame quantity is:" )
        call output ( qty%instanceOffset, advance='yes' )
      endif

      if ( .not. IsModuleSpacecraft(instrumentModule) ) then
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // "tpGeodAlt"
        
        l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
        call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
          & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & hdfVersion=hdfVersion )
        if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        
        ! Now we're going to deal with a VGrid for this quantity
        qty%verticalCoordinate = l_geodAltitude
        qty%surfs = l1bField%dpField(1,:,:)  ! Vert coord is tpGeodAlt read above.

        call DeallocateL1BData(l1bfield)

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
        call DeallocateL1BData(l1bfield)
      end if
    end if
    qty%frequencyCoordinate = L_None
    qty%instrumentModule = instrumentModule

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.
  end subroutine ConstructMinorFrameQuantity

  ! ------------------------------------ ForgeMinorFrames --------------
  subroutine ForgeMinorFrames ( root, chunk, mifGeolocation, vGrids )
    ! This routine is used when we're in the mode of creating l2pc files
    ! and we want to invent a set of minor frame quantities with no
    ! reference to the l1 files

    use Chunks_m, only: MLSChunk_T
    use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: F_GEODANGLE, F_MODULE, F_NOMIFS, &
      & F_SOLARTIME, F_SOLARZENITH, F_GEODALT
    use INIT_TABLES_MODULE, only: L_GEODALTITUDE, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, &
      & PHYQ_TIME
    use MLSCommon, only: RK => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_T, SetupNewQuantityTemplate
    use TREE, only: DECORATION, NSONS, SUBTREE
    use VgridsDatabase, only: VGRID_T

    ! Dummy arguments
    integer, intent(in) :: ROOT         ! Tree vertex
    type (MLSChunk_T), intent(in) :: CHUNK ! This chunk
    type (QuantityTemplate_T), dimension(:), intent(inout) :: MIFGEOLOCATION
    type (VGrid_T), dimension(:), pointer :: VGRIDS

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
    type(VGrid_T), pointer :: GEODALT

    ! Executable code
    nullify ( geodAlt )
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
      case ( f_geodAlt )
        geodAlt => vGrids(decoration(decoration(subtree(2,son))))
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
    ! Note this will destroy the old ones contents bit by bit.
    call SetupNewQuantityTemplate ( mifGeolocation(instrumentModule), &
      & noInstances=noMAFs, noSurfs=noMIFs, noChans=1,&
      & coherent=.false., stacked=.false., regular=.true.,&
      & minorFrame=.true. )

    ! Put zeros etc. in the appropriate places, may overwrite these later
    mifGeolocation(instrumentModule)%instrumentModule = instrumentModule
    mifGeolocation(instrumentModule)%noInstancesLowerOverlap = 0
    mifGeolocation(instrumentModule)%noInstancesUpperOverlap = 0
    mifGeolocation(instrumentModule)%verticalCoordinate = l_geodAltitude
    mifGeolocation(instrumentModule)%time = 0.0
    mifGeolocation(instrumentModule)%phi = 0.0
    mifGeolocation(instrumentModule)%geodLat = 0.0
    mifGeolocation(instrumentModule)%solarTime = 0.0
    mifGeolocation(instrumentModule)%solarZenith = 0.0
    mifGeolocation(instrumentModule)%lon = 0.0
    mifGeolocation(instrumentModule)%losAngle = 0.0
    mifGeolocation(instrumentModule)%surfs = 0.0

    ! Check the geodAlt information if supplied
    if ( associated ( geodAlt ) ) then 
      if ( noMIFs /= geodAlt%noSurfs ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Incorrect number of minor frames for geodAlt' )
      mifGeolocation(instrumentModule)%surfs = spread ( geodAlt%surfs(:,1), 2, noMAFs )
    end if

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
  subroutine Announce_Error ( where, message, extra, severity )

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: BLANKS, OUTPUT
    use TREE, only: SOURCE_REF
    use Intrinsic, only: LIT_INDICES
    use String_Table, only: DISPLAY_STRING
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA
    character(len=*), intent(in), optional :: SEVERITY ! 'nonfatal' or 'fatal'
    character(len=8) :: mySeverity

    mySeverity = 'fatal'
    if ( present(severity) ) then
      if (index('WwNn', severity(1:1)) > 0 ) mySeverity='warning'
    endif
    if ( mySeverity /= 'fatal' ) then
      call blanks(5)
      call output ( ' (warning) ' )
    else
      call blanks(5, fillChar='*')
      call output ( ' (fatal) ' )
    endif
    call output ( 'At ' )
    if ( where > 0 ) then
      call print_source ( source_ref(where) )
    else
      call output ( '(no lcf tree available)' )
    end if
    call output ( ': ' )
    call output ( message )
    if ( present ( extra ) ) call display_string ( lit_indices ( extra ), strip=.true. )
    call output ( '', advance='yes' )
    if ( mySeverity == 'fatal') &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Problem in Construct' )
  end subroutine Announce_Error

  ! -----------------------------------------------  AnyGoodSignalData  -----
  function AnyGoodSignalData ( signal, sideband, l1bInfo, Chunk )  result (answer)
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  ! Arguments

    use Chunks_m, only: MLSChunk_T
    use Allocate_Deallocate, only: Deallocate_Test
    use L1BData, only: L1BData_T, READL1BDATA, &
      & FindL1BData, AssembleL1BQtyName, PRECISIONSUFFIX
    use MLSCommon, only: L1BInfo_T, RK => R8
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

  ! ----------------------------------  PointQuantityToHGrid  -----
  subroutine PointQuantityToHGrid ( hGrid, qty )

  ! This routine copies HGrid information into an already defined quantity

    use HGridsDatabase, only: hGrid_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_T

    ! Dummy arguments
    type (hGrid_T), intent(in) :: HGRID
    type (QuantityTemplate_T), intent(inout) :: QTY

    ! Executable code
    if ( qty%noInstances/=hGrid%noProfs ) call MLSMessage ( MLSMSG_Error,&
      & ModuleName, "Size of HGrid not compatible with size of quantity" )
    if ( .not. qty%stacked ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Cannot copy hGrids into unstacked quantities")

    qty%phi => hGrid%phi
    qty%geodLat => hGrid%geodLat
    qty%lon => hGrid%lon
    qty%time => hGrid%time
    qty%solarTime => hGrid%solarTime
    qty%solarZenith => hGrid%solarZenith
    qty%losAngle => hGrid%losAngle
    qty%noInstancesLowerOverlap = hGrid%noProfsLowerOverlap
    qty%noInstancesUpperOverlap = hGrid%noProfsUpperOverlap

  end subroutine PointQuantityToHGrid

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

  ! --------------------------------  ConstructMajorFrameQuantity  -----
  subroutine ConstructMajorFrameQuantity( chunk, instrumentModule, qty, noChans, &
    & mifGeolocation )
    ! Dummy arguments
    use Chunks_m, only: MLSChunk_T
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
      & noChans=noChans, sharedHGrid=.true., sharedVGrid=.true. )
    call SetupEmptyVGridForQuantity ( qty )

    qty%phi => source%phi(1:1,:)
    qty%geodLat => source%geodLat(1:1,:)
    qty%lon => source%lon(1:1,:)
    qty%time => source%time(1:1,:)
    qty%solarTime => source%solarTime(1:1,:)
    qty%solarZenith => source%solarZenith(1:1,:)
    qty%losAngle => source%losAngle(1:1,:)

    qty%majorFrame = .true.
    qty%minorFrame = .false.
    qty%instanceOffset = source%instanceOffset
    qty%instrumentModule = source%instrumentModule
    qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
    qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap
  end subroutine ConstructMajorFrameQuantity
  
  ! ----------------------------------------------- InitQuantityTemplates ----
  subroutine InitQuantityTemplates
    ! This routine initializes the quantity template properties
    ! This is the routine one needs to update when one introduces a new quantity type.
    use Init_Tables_Module, only:  L_ADOPTED, L_BASELINE, L_BOUNDARYPRESSURE, &
      L_CALSIDEBANDFRACTION, &
      L_CHISQBINNED, L_CHISQCHAN, L_CHISQMMAF, L_CHISQMMIF, L_CLOUDICE, &
      L_CLOUDINDUCEDRADIANCE, L_CLOUDEXTINCTION, L_CLOUDMINMAX, L_CLOUDRADSENSITIVITY, &
      L_CLOUDWATER, L_COLUMNABUNDANCE, &
      L_DNWT_AJN, L_DNWT_AXMAX, L_DNWT_CAIT, &
      L_DNWT_CHISQMINNORM, L_DNWT_CHISQNORM, &
      L_DNWT_DIAG, L_DNWT_DXDX, L_DNWT_DXDXL, &
      L_DNWT_DXN, L_DNWT_DXNL, L_DNWT_FLAG, L_DNWT_FNMIN, &
      L_DNWT_FNORM, L_DNWT_GDX, L_DNWT_GFAC, &
      L_DNWT_GRADN, L_DNWT_SQ, L_DNWT_SQT,&
      L_EARTHRADIUS, L_EARTHREFL, L_ECRTOFOV, L_EFFECTIVEOPTICALDEPTH, &
      L_ELEVOFFSET, L_EXTINCTION, &
      L_FIELDAZIMUTH, L_FIELDELEVATION, L_FIELDSTRENGTH, &
      L_GPH, L_HEIGHTOFFSET, &
      L_ISOTOPERATIO, L_JACOBIAN_COLS, L_JACOBIAN_ROWS, &
      L_L1BMAFBASELINE, L_LIMBSIDEBANDFRACTION, L_LOSTRANSFUNC, L_LOSVEL, &
      L_MASSMEANDIAMETERICE, L_MASSMEANDIAMETERWATER, L_MAGNETICFIELD, &
      L_NOISEBANDWIDTH, L_NORADSPERMIF, L_NORADSBINNED, &
      L_NUMJ, L_OPTICALDEPTH, L_ORBITINCLINATION, &
      L_PHASETIMING, L_PHITAN, L_PTAN, L_QUALITY, L_RADIANCE, &
      L_REFGPH, L_REFLTEMP, L_REFLTRANS, L_REFLREFL, L_REFLSPILL, &
      L_RHI, L_SINGLECHANNELRADIANCE, L_SIZEDISTRIBUTION, &
      L_SCANRESIDUAL, L_SCECI, L_SCVEL, L_SCVELECI, &
      L_SCVELECR, L_SCGEOCALT, L_SPACERADIANCE, &
      L_STATUS, L_STRAYRADIANCE, L_SURFACETYPE, L_SYSTEMTEMPERATURE, &
      L_TEMPERATURE, L_TNGTECI, L_TNGTGEODALT, L_TNGTGEOCALT, &
      L_VMR
    use Init_Tables_Module, only: PHYQ_EXTINCTION, PHYQ_FREQUENCY,&
      & PHYQ_GAUSS, PHYQ_IceDensity, PHYQ_LENGTH, &
      & PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_VELOCITY, &
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
      l_adopted, phyq_dimensionless, none, next, &
      l_baseline, phyq_temperature, p_flexibleVHGrid, p_fGrid, p_radiometerOptional, &
                  p_signalOptional, p_suppressChannels, p_mustBeZeta, next, &
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
      l_cloudMinMax, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, &
                     p_signal, p_suppressChannels, next, &
      l_cloudRadSensitivity, phyq_temperature, p_minorFrame, p_signal, next, &
      l_cloudWater, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_columnAbundance, phyq_dobsonunits, p_hGrid, p_molecule, next, &
      l_dnwt_ajn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_axmax, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_cait, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_chisqminnorm, phyq_dimensionless, p_vGrid /) )

    call DefineQtyTypes ( (/ &
      l_dnwt_chisqnorm, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_diag, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdx, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdxl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxnl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_flag, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_fnmin, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_fnorm, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gdx, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gfac, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_gradn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_sq, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_sqt, phyq_dimensionless, p_vGrid, next, &
      l_earthRadius, phyq_length, p_hGrid, next, &
      l_earthRefl, phyq_dimensionless, none /) )

    call DefineQtyTypes ( (/ &
      l_ecrToFOV, phyq_dimensionless, p_minorFrame, p_module, p_matrix3x3, next, &
      l_effectiveOpticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_elevOffset, phyq_angle, p_signal, next, &
      l_extinction, phyq_extinction, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                    p_mustBeZeta, next, &
      l_fieldAzimuth, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldElevation, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldStrength, phyq_gauss, p_hGrid, p_vGrid, next, &
      l_gph, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_heightOffset, phyq_length, p_hGrid, p_vGrid, next, &
      l_isotopeRatio, phyq_dimensionless, p_molecule, next, &
      l_jacobian_cols, phyq_dimensionless, p_vGrid, next, &
      l_jacobian_rows, phyq_dimensionless, p_vGrid, next, &
      l_l1bMAFBaseline, phyq_temperature, p_majorFrame, p_signal, next, &
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
      l_numJ, phyq_dimensionless, p_vGrid, next, &
      l_opticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_orbitInclination, phyq_angle, p_minorFrame, p_scModule, next, &
      l_phiTan, phyq_angle, p_minorFrame, p_module, next, & 
      l_ptan, phyq_zeta, p_minorFrame, p_module, next /) )

    call DefineQtyTypes ( (/ & 
      l_quality, phyq_dimensionless, p_hGrid, next, &
      l_radiance, phyq_temperature, p_minorFrame, p_signal, next, & 
      l_refltemp, phyq_temperature, p_majorFrame, p_reflector, p_module, next, &
      l_refltrans, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflrefl, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflspill, phyq_temperature, p_signal, p_majorframe, p_reflector, next, &
      l_refGPH, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_rhi, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, p_mustBeZeta, next, &
      l_scECI, phyq_length, p_minorFrame, p_scModule, p_xyz, next, &
      l_scGeocAlt, phyq_length, p_minorFrame, p_scModule, next, &
      l_scVel, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECI, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECR, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scanResidual, phyq_length, p_minorFrame, p_module, next, &
      l_singleChannelRadiance, phyq_temperature, p_minorFrame, p_signal, &
                               p_suppressChannels, next, &
      l_sizeDistribution, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, & 
      l_spaceRadiance, phyq_temperature, none, next, &
      l_status, phyq_dimensionless, p_hGrid, next, &
      l_strayRadiance, phyq_temperature, p_signal, p_majorFrame, next, &
      l_surfaceType, phyq_dimensionless, p_hGrid, next, & 
      l_systemTemperature, phyq_temperature, p_signal, next, &
      l_temperature, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, next, &
      l_tngtECI, phyq_length, p_minorFrame, p_module, p_xyz, next, &
      l_tngtGeocAlt, phyq_length, p_minorFrame, p_module, next, &
      l_tngtGeodAlt, phyq_length, p_minorFrame, p_module, next, &
      l_vmr, phyq_vmr, p_hGrid, p_vGrid, p_fGridOptional, p_molecule, &
             p_radiometerOptional, p_mustbezeta, next /) )

    call DefineQtyTypes ( (/ & 
      l_phaseTiming, phyq_dimensionless, p_vGrid, next/) )
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
      if ( propertyTable ( p_mustBeZeta, i ) .and. .not. &
        & ( propertyTable ( p_vGrid, i ) .or. propertyTable ( p_flexibleVHGrid, i ) ) ) then
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

  ! ---------------------------------- SetupEmptyHGridForQuantity
  subroutine SetupEmptyHGridForQuantity ( qty ) 
    use Allocate_Deallocate, only: ALLOCATE_TEST
    use QuantityTemplates, only: QuantityTemplate_T
    ! Dummy arguments
    type ( QuantityTemplate_T ), intent(inout) :: QTY
    ! Executable code
    qty%sharedHGrid = .false.
    qty%hGridIndex = 0
    nullify( qty%frequencies, qty%geodLat, qty%lon, qty%time, qty%solarTime, &
      & qty%solarZenith, qty%losAngle ) ! Lest we deallocate a database entry
    call Allocate_test ( qty%phi, 1, 1, 'qty%phi(1,1)', ModuleName )
    call Allocate_test ( qty%geodLat, 1, 1, 'qty%geodLat(1,1)', ModuleName )
    call Allocate_test ( qty%lon, 1, 1, 'qty%lon(1,1)', ModuleName )
    call Allocate_test ( qty%time, 1, 1, 'qty%time(1,1)', ModuleName )
    call Allocate_test ( qty%solarTime, 1, 1, 'qty%solarTime(1,1)', ModuleName )
    call Allocate_test ( qty%solarZenith, 1, 1, 'qty%solarZenith(1,1)', ModuleName )
    call Allocate_test ( qty%losAngle, 1, 1, 'qty%losAngle(1,1)', ModuleName )
    qty%phi = 0.0
    qty%geodLat = 0.0
    qty%lon = 0.0
    qty%time = 0.0
    qty%solarTime = 0.0
    qty%solarZenith = 0.0
    qty%losAngle = 0.0
  end subroutine SetupEmptyHGridForQuantity

  ! ---------------------------------- SetupEmptyVGridForQuantity
  subroutine SetupEmptyVGridForQuantity ( qty ) 
    use Allocate_Deallocate, only: ALLOCATE_TEST
    use Intrinsic, only: L_NONE
    use QuantityTemplates, only: QuantityTemplate_T
    ! Dummy arguments
    type ( QuantityTemplate_T ), intent(inout) :: QTY
    ! Executable code
    qty%sharedVGrid = .false.
    qty%vGridIndex = 0
    qty%verticalCoordinate = l_none
    nullify(qty%surfs) ! Lest we deallocate a database entry
    call Allocate_test ( qty%surfs, 1, 1, 'qty%surfs(1,1)', ModuleName )
  end subroutine SetupEmptyVGridForQuantity

end module ConstructQuantityTemplates
!
! $Log$
! Revision 2.119  2005/01/07 01:03:19  vsnyder
! Remove unused declarations
!
! Revision 2.118  2004/10/16 17:25:55  livesey
! Added signalOptional property and more flexible handling of baseline
!
! Revision 2.117  2004/10/13 02:24:33  livesey
! Added ability to set altitude in forge
!
! Revision 2.116  2004/09/27 20:11:05  livesey
! Added L1BMAFBaseline
!
! Revision 2.115  2004/08/26 18:51:03  pwagner
! Failsafe feature for shared fgrids
!
! Revision 2.114  2004/08/16 23:42:54  livesey
! Added p_flexibleVHGrid in order to allow minor frame dependent
! baselines.
!
! Revision 2.113  2004/06/29 00:09:25  pwagner
! New phaseTiming type
!
! Revision 2.112  2004/06/21 23:58:40  pwagner
! numJ also mad e diagnostic qty that may be written to dgm file
!
! Revision 2.111  2004/06/17 23:17:00  pwagner
! Retrieval diagnostics freed from need to be l2gp
!
! Revision 2.110  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.109  2004/04/16 00:48:52  livesey
! Added singleChannelRadiance type
!
! Revision 2.108  2004/03/17 17:16:25  livesey
! New cloudMinMax type.
!
! Revision 2.107  2004/02/10 21:17:24  livesey
! Added status and quality as valid quantity types
!
! Revision 2.106  2004/01/24 01:03:46  livesey
! Added the adopted option etc.
!
! Revision 2.105  2003/08/08 23:06:53  livesey
! Added the fieldStrength etc. quantities.
!
! Revision 2.104  2003/07/07 20:22:37  livesey
! Minor bug fix in minor frame quantities
!
! Revision 2.103  2003/07/01 23:18:04  livesey
! Added setting of grandTotalInstances for minor frame quantities.
!
! Revision 2.102  2003/07/01 19:29:32  livesey
! Added filling of grandTotalInstances
!
! Revision 2.101  2003/06/27 00:07:07  pwagner
! Made print offset depend on qtmp switch; restored some logs
!
! (logging disabled from Rev. 2.92-2.100)
! Revision 2.91  2003/05/22 04:04:07  livesey
! elevOffset is now a channel dependent rather than radiometer dependent
! quantity.
!
! Revision 2.90  2003/05/10 23:39:25  livesey
! Fixed problems with chisqs and noRads
!
! Revision 2.89  2003/05/07 01:00:03  livesey
! More stuff got missed in the merge, was I asleep or something?
!
! Revision 2.88  2003/02/13 19:05:39  vsnyder
! Move USEs from module to procedure scope, cosmetic changes
!
! Revision 2.87  2003/02/12 21:53:48  pwagner
! Fix to errant ANY_GOOD_SIGNALDATA in case of hdf5 files
!
! Revision 2.86  2003/02/07 03:42:50  vsnyder
! Fill NATURAL_UNITS on first call only -- it's a SAVE variable now
!
! Revision 2.85  2003/02/07 03:38:05  vsnyder
! Cosmetic change
!
! Revision 2.84  2003/02/07 00:41:41  livesey
! Bug fix/workaround
!
! Revision 2.83  2003/02/06 23:31:00  livesey
! New approach to Forge
!
! Revision 2.82  2003/01/14 00:40:29  pwagner
! Moved GetQuantityAttributes to L2AUXData
!
! Revision 2.81  2003/01/09 00:09:15  pwagner
! routine GetQuantityAttributes added
!
! Revision 2.80  2003/01/08 23:50:44  livesey
! Added irregular argument
!
! Revision 2.79  2003/01/07 23:46:53  livesey
! Added magnetic field stuff
!
! Revision 2.78  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.77  2002/11/26 23:37:50  livesey
! Better handling of major frame quantities
!
! Revision 2.76  2002/11/22 12:16:08  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.75  2002/11/13 01:05:03  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.74  2002/10/08 17:36:19  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.73  2002/09/26 18:03:06  livesey
! Changed extinction to a vmr
!
! Revision 2.72  2002/09/25 20:07:55  livesey
! Made -g less verbose
!
! Revision 2.71  2002/09/24 21:37:44  livesey
! Added minValue stuff
!
! Revision 2.70  2002/09/24 00:27:16  pwagner
! Wont bomb if no l1brads; nor whine if no good signals
!
! Revision 2.69  2002/09/18 22:48:32  pwagner
! Chooses signals first band with any good data
!
! Revision 2.68  2002/08/04 16:01:19  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.67  2002/06/14 16:39:49  livesey
! Orbital inclination now minor frame
!
! Revision 2.66  2002/06/04 22:07:35  livesey
! Added phiTan as a state vector element
!
! Revision 2.65  2002/05/22 19:06:32  jonathan
! added units for cloudextinction(m-1), totalextinction(m-1), and
! massmeandiameterice(micron m) as dimentionless for now, may define
! more clear later
!
! Revision 2.64  2002/05/14 00:27:42  livesey
! New code for system temperatures and noise bandwidths
!
! Revision 2.63  2002/05/07 20:02:54  livesey
! Added noise bandwidth
!
! Revision 2.62  2002/04/10 17:44:22  pwagner
! Added rhi quantity (but is this enough?)
!
! Revision 2.61  2002/03/19 00:51:32  pwagner
! Added new scVel quantity types
!
! Revision 2.60  2002/02/09 21:35:52  livesey
! Added optical depth stuff
!
! Revision 2.59  2001/11/08 00:13:38  livesey
! Sorted out extinction stuff
!
! Revision 2.58  2001/10/31 19:07:25  livesey
! Added fGrid stuff
!
! Revision 2.57  2001/10/12 23:15:05  pwagner
! Fixed biggest erros in diagnostic quantity templates
!
! Revision 2.56  2001/10/02 23:12:50  pwagner
! More chi^2 fixes
!
! Revision 2.55  2001/10/02 20:50:54  livesey
! No longer asserts baseline to be minor frame
!
! Revision 2.54  2001/09/17 23:11:50  pwagner
! Tiny changes for chi^2
!
! Revision 2.53  2001/09/17 21:58:50  livesey
! Added allocate of frequencies if needed
!
! Revision 2.52  2001/09/14 23:30:33  pwagner
! Now constructs major frame quantity templates
!
! Revision 2.51  2001/08/01 00:04:29  dwu
! add qty%frequencies = VGrids(sGridIndex)%surfs for quantity l_losTransFunc
!
! Revision 2.50  2001/07/30 23:28:38  pwagner
! Added columnAbundances scaffolding--needs fleshing out
!
! Revision 2.49  2001/07/26 17:34:25  jonathan
! add DTcir, etc, jonathan
!
! Revision 2.48  2001/07/19 22:17:44  jonathan
! add cloud stuff , jonathan/dwu
!
! Revision 2.47  2001/07/19 17:42:31  dwu
! add sGrid field
!
! Revision 2.46  2001/07/18 23:17:30  dwu
! rename l_radiusofearth as l_earthradius
!
! Revision 2.45  2001/07/18 18:42:19  dwu
! add radiusofearth quantity type
!
! Revision 2.44  2001/07/17 23:23:05  dwu
! make l_losTransFunc as non-minorframe but minorframe-like quantity
!
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
