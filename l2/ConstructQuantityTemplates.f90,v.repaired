! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ConstructQuantityTemplates

  ! This module is responsible for constructing templates for quantities.
  ! This version is a rewrite, aimed at tidying up a lot of the codebase
  use init_tables_module, only: first_lit, last_lit

  implicit none

  private 

  ! The various properties has/can have
  integer, parameter ::         Next                 = -1
  integer, public, parameter :: FirstProperty        = 1
  integer, public, parameter :: P_Chunked            = FirstProperty
  integer, public, parameter :: P_MajorFrame         = P_Chunked + 1
  integer, public, parameter :: P_MinorFrame         = P_MajorFrame + 1
  integer, public, parameter :: P_MustBeZeta         = P_MinorFrame + 1
  integer, public, parameter :: P_MustBeGeocalt      = P_MustBeZeta + 1
  integer, public, parameter :: P_FGrid              = P_MustBeGeocAlt + 1
  integer, public, parameter :: P_FGridOptional      = P_FGrid + 1
  integer, public, parameter :: P_FlexibleVHGrid     = P_FGridOptional + 1
  integer, public, parameter :: P_HGrid              = P_FlexibleVHGrid + 1
  integer, public, parameter :: P_Module             = P_HGrid + 1
  integer, public, parameter :: P_Molecule           = P_Module + 1
  integer, public, parameter :: P_SGrid              = P_Molecule + 1
  integer, public, parameter :: P_VGrid              = P_SGrid + 1
  integer, public, parameter :: P_Radiometer         = P_VGrid + 1
  integer, public, parameter :: P_RadiometerOptional = P_Radiometer + 1
  integer, public, parameter :: P_Reflector          = P_RadiometerOptional + 1
  integer, public, parameter :: P_SCmodule           = P_Reflector + 1
  integer, public, parameter :: P_Signal             = P_SCmodule + 1
  integer, public, parameter :: P_SignalOptional     = P_Signal + 1
  integer, public, parameter :: P_SuppressChannels   = P_SignalOptional + 1
  integer, public, parameter :: P_XYZ                = P_SuppressChannels + 1
  integer, public, parameter :: P_Matrix3x3          = P_XYZ + 1

  integer, public, parameter :: LastProperty         = P_Matrix3x3

  ! Local, saved variables (constant tables really)
  logical, public, save :: &
    & PROPERTYTABLE ( firstProperty: LastProperty, first_lit : last_lit )
  integer, public, save :: UNITSTABLE ( first_lit : last_lit )

  public :: AnyGoodSignalData, ConstructMinorFrameQuantity, ConstructMajorFrameQuantity
  public :: CreateQtyTemplateFromMLSCFInfo
  public :: ForgeMinorFrames
  public :: InitQuantityTemplates, GetQtyTypeIndex

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  logical, save :: FIRSTCALL = .true.

contains ! ============= Public procedures ===================================

  ! ------------------------------------------- CreateQtyTemplateFromMLSCFInfo ----
  type (QuantityTemplate_T) function CreateQtyTemplateFromMLSCFInfo ( &
    & Name, Root, FGrids, HGrids, filedatabase, Chunk, MifGeolocation ) &
    & result ( QTY )
    use allocate_deallocate, only: allocate_test, deallocate_test
    use chunks_m, only: mlschunk_t
!   use chunkdivide_m, only: chunkdivideconfig
    use expr_m, only: expr
    use FGrid, only: FGrid_t
    use HGridsDatabase, only: HGrids_t
    use HighOutput, only: BeVerbose, outputNamedValue
    use init_tables_module, only:  f_badvalue, f_coordinate, f_fgrid, f_hgrid, &
      & f_irregular, f_keepchannels, f_logbasis, f_minvalue, f_module, &
      & f_molecule, f_radiometer, f_reflector, f_sgrid, f_signal, f_stacked, &
      & f_type, f_vgrid, f_xgrid, field_first, field_last, l_channel, &
      & l_explicit, l_geocaltitude, l_lostransfunc, l_matrix3x3, l_none, &
      & l_phitan, l_true, l_xyz, l_zeta
    use MLSCommon, only: MLSFile_t
    use MLSKinds, only: rk => r8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSSignals_m, only:getModuleFromRadiometer, getModuleFromSignal, &
      & getRadiometerFromSignal, getSignal, signal_t, &
      & isModuleSpacecraft
    use moreTree, only: get_boolean
    use parse_signal_m, only: parse_signal
    use quantityTemplates, only: nullifyQuantityTemplate, pointQuantityToHGrid, &
      & quantityTemplate_t, setupNewQuantityTemplate
    use string_table, only: get_string
    use toggles, only: gen, levels, toggle
    use trace_m, only: trace_begin, trace_end
    use tree, only: decoration, node_id, nsons, sub_rosa, subtree
    use tree_types, only: n_set_one
    use VGridsDatabase, only: VGrids

    ! Dummy arguments
    integer, intent(in) :: NAME              ! Sub-rosa index of name
    integer, intent(in) :: ROOT              ! Root of QuantityTemplate subtree
    type (FGrid_T), dimension(:), pointer :: FGrids
    type (HGrids_T), dimension(:), pointer :: HGrids
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (MLSChunk_T), intent(in) :: Chunk
    type (QuantityTemplate_T), dimension(:), intent(in), optional, target :: &
      & MifGeolocation ! TARGET attribute needed so pointers to MifGeolocation
                       ! created at lower levels in the call tree don't become
                       ! undefined when those procedures return.

    ! Local variables
    logical, pointer :: Channels(:)     ! From Parse_Signal
    logical :: LOGBASIS                 ! To place in quantity
    logical :: ISMINORFRAME             ! Is a minor frame quantity
    logical :: KeepChannels             ! From /channels, means keep the channels
                                        ! information from the signal
    logical :: REGULAR                  ! Flag
    logical :: PROPERTIES(firstProperty: LastProperty) ! Properties for this quantity type
    logical :: GOT(field_first:field_last) ! Fields
    character(len=127) :: SIGNALSTRING

    integer :: Coordinate               ! For vertical coordinate, in case we
                                        ! don't want L_None for MIF qty, e.g.
                                        ! tngtGeocAlt
    integer :: FGRIDINDEX               ! Index of frequency grid
    integer :: FREQUENCYCOORDINATE      ! Literal
    integer :: HGRIDINDEX               ! Index of horizontal grid
    integer :: HORIZONTALCOORDINATE     ! Literal
    integer :: I                        ! Loop counter
    integer :: INSTRUMENTMODULE         ! Database index
    integer :: KEY                      ! Field name, F_...
    integer :: Me = -1                  ! String index for trace
    integer :: MOLECULE                 ! Literal
    integer :: NOCHANS                  ! Quantity dimension
    integer :: NoCrossTrack             ! Number of cross-track angles
    integer :: NOINSTANCES              ! Quantity dimension
    integer :: NOSURFS                  ! Quantity dimension
    integer :: QUANTITYTYPE             ! Literal
    integer :: RADIOMETER               ! Database index
    integer :: REFLECTOR                ! Reflector literal
    integer :: SGRIDINDEX               ! Index for 'sGrid'
    integer :: SIDEBAND                 ! -1, 0, 1
    integer :: SIGNAL                   ! Database index
    integer :: SON                      ! A Son of Root -- an n_assign node
    logical :: Stacked                  ! All heights at a particular instance
                                        ! have same Lat and Lon
    integer :: VGRIDINDEX               ! Index in database
    integer :: VALUE                    ! Node index of value of field of spec
    integer :: S_INDEX                  ! Loop counter
    logical :: verbose
    integer :: xGridIndex               ! Index in hGrid database

    integer, dimension(2) :: EXPR_UNITS
    integer, dimension(:), pointer :: SignalInds ! From parse signal

    real(rk) :: BadValue
    real(rk) :: MINVALUE                ! Minimum value allowed for quantity in fwm
    real(rk), dimension(2) :: EXPR_VALUE
    type (signal_T) :: SIGNALINFO       ! Details of the appropriate signal

    ! Executable code
    if ( firstCall ) then
      call InitQuantityTemplates
      firstCall = .false.
    end if
    verbose = BeVerbose( 'qtmp', 0 )

    call trace_begin ( me, "CreateQtyTemplate", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 )

    ! Set appropriate defaults
    call NullifyQuantityTemplate ( qty ) ! for Sun's rubbish compiler
    nullify ( channels, signalInds )
    qty%name = name
    fGridIndex = 0
    hGridIndex = 0
    horizontalCoordinate = l_phiTan
    instrumentModule = 0
    keepChannels = .false.
    logBasis = .false.
    minValue = -huge ( 0.0_rk )
    molecule = 0
    noChans = 1
    noCrossTrack = 1
    quantitytype = 0
    radiometer = 0
    reflector = 0
    regular = .true.
    sGridIndex = 0
    sideband = 0
    signal = 0
    signalString = ''
    stacked = .true.
    vGridIndex = 0
    xGridIndex = 0

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
      case ( f_badValue )
        call expr ( subtree(2,son), expr_units, expr_value )
        badValue = expr_value(1)
      case ( f_coordinate )
        coordinate = value
      case ( f_fgrid )
        fGridIndex = decoration(value)
      case ( f_hgrid )
        hGridIndex = decoration(value)
      case ( f_keepChannels )
        keepChannels = get_boolean(son)
      case ( f_logBasis )
        logBasis = get_boolean(son)
      case ( f_irregular )
        regular = get_boolean(son)
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
      case ( f_reflector )
        reflector = value
      case ( f_sgrid )
        sGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_signal )
        !??? For the moment it is simple, later we'll be more intelligent here
        !??? for example, letting the user choose either R1A or R1B.
        call get_string( sub_rosa(subtree(2,son)), signalString, strip=.true. )
        !??? Here we would do intelligent stuff to work out which bands
        !??? are present, for the moment choose the first
        call parse_Signal ( signalString, signalInds, &
          & tree_index=son, sideband=sideband, channels=channels )
        if ( .not. associated(signalInds) ) then ! A parse error occurred
          call MLSMessage ( MLSMSG_Error, ModuleName,&
            & 'Unable to parse signal string' )
        end if
        if ( size(signalInds) == 1 .or. .not. associated(filedatabase) ) then
          signal = signalInds(1)
        else
          ! Seek a signal with any precision values !< 0
          do s_index=1, size(signalInds)
            if ( AnyGoodSignalData ( signalInds(s_index), sideband, &
              & filedatabase, Chunk) ) exit
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
      case ( f_stacked )
        stacked = get_boolean(son)
      case ( f_type )
        quantityType = value
      case ( f_vgrid )
        vGridIndex = decoration(value) ! node_id(value) == n_spec_args
      case ( f_xgrid ) ! an hGrid
        xGridIndex = decoration(value) ! node_id(value) == n_spec_args
        if ( hGrids(xGridIndex)%the_hGrid%type /= l_explicit ) &
          & call Announce_error ( root, 'XGrid is not explicit', quantityType )
        noCrossTrack = size(hgrids(xGridIndex)%the_hGrid%phi)
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
    end if
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
    if ( properties ( p_mustBeGeocAlt ) .and. got(f_vGrid) ) then
      if ( vGrids(vGridIndex)%verticalCoordinate /= l_geocAltitude ) &
        & call Announce_error ( root, 'Expecting geocentric altitude coordinates for ', &
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
      if ( keepChannels ) then
        qty%channels => channels
        nullify ( channels ) ! so we don't deallocate out from under qty%channels
      end if
    else if ( got(f_sGrid) ) then
      ! Uses an sGrid.  This is faked as frequency, but it's not
      noChans = size ( vGrids(sGridIndex)%surfs )
      frequencyCoordinate = l_losTransFunc
    else
      ! No frequency variation
      noChans = 1
      frequencyCoordinate = l_none
    end if

    ! Now deal with the flexibleVHGrid quantities
    if ( properties(p_flexibleVHGrid) ) then
      if ( got ( f_hGrid ) .neqv. got ( f_vGrid ) ) &
        & call Announce_Error ( root, 'Must supply both or neither vGrid and hGrid' )
      isMinorFrame = .not. got ( f_hGrid )
    else
      isMinorFrame = properties(p_minorFrame)
    end if
    
    ! horizontalCoordinate?
    if ( got(f_hGrid) ) then
      if (  hGridIndex > 0 ) &
        & horizontalCoordinate = HGrids(hGridIndex)%the_hGrid%masterCoordinate
    else if ( got(f_xGrid) ) then
      call Announce_error ( root, 'XGrid specified without HGrid', quantityType )
    end if

    ! Now do the setup for the different families of quantities
    if ( isMinorFrame ) then
      ! Setup a minor frame quantity
      if ( got(f_coordinate) ) qty%verticalCoordinate = coordinate
      qty%quantityType = quantityType
      call ConstructMinorFrameQuantity ( instrumentModule, qty, &
        noChans=noChans, noCrossTrack=noCrossTrack, filedatabase=filedatabase, &
        chunk=chunk, mifGeolocation=mifGeolocation, &
        regular=regular )
    else if ( properties(p_majorFrame) ) then
      ! Setup a major frame quantity
      call ConstructMajorFrameQuantity ( chunk, instrumentModule, &
        & qty, noChans, mifGeolocation, noCrossTrack )
    else
      ! Setup a non major/minor frame quantity
      noInstances = 1
      if ( got ( f_hGrid ) ) noInstances = hGrids(hGridIndex)%the_hGrid%noProfs
      noSurfs = 1
      if ( got ( f_vGrid ) ) noSurfs = vGrids(vGridIndex)%noSurfs

      ! Setup the quantity template
      call SetupNewQuantityTemplate ( qty, noInstances=noInstances, &
        & noSurfs=noSurfs, coherent=.true., stacked=stacked, regular=.true.,&
        & noChans=noChans, noCrossTrack=noCrossTrack, &
        & sharedVGrid=.true., sharedFGrid=.true. )

      ! Setup the horizontal coordinates
      if ( got(f_hGrid) ) then
        qty%the_hGrid => hGrids(hGridIndex)%the_hGrid
        qty%hGridIndex = hGridIndex
        qty%xGridIndex = xGridIndex
        call PointQuantityToHGrid ( qty )
      else
        call SetupEmptyHGridForQuantity ( qty )
      end if
      ! Work out the instance offset
      if ( associated ( chunk%hGridOffsets ) ) then
        if ( got ( f_hGrid )  ) then
          qty%instanceOffset = chunk%hGridOffsets(hGridIndex)
          qty%grandTotalInstances = chunk%hGridTotals(hGridIndex)
          if ( verbose ) call outputnamedvalue ( 'grandTotalInstances from hgrid)', &
            & qty%grandTotalInstances )
          ! Subtract any instances outside processing range
          !   if ( ChunkDivideConfig%allowPriorOverlaps ) &
          !     & qty%grandTotalInstances = qty%grandTotalInstances - &
          !     & hGrids(hGridIndex)%the_hGrid%noProfsLowerOverlap
          !   if ( ChunkDivideConfig%allowPostOverlaps ) &
          !     & qty%grandTotalInstances = qty%grandTotalInstances - &
          !     & hGrids(hGridIndex)%the_hGrid%noProfsUpperOverlap
        else
          ! Must have a single instance per chunk
          qty%instanceOffset = chunk%chunkNumber - 1
          ! -1 because it's an offset remember, not an origin.
          qty%grandTotalInstances = 0
        end if
      end if

      ! Setup the vertical coordinates
      if ( got (f_vGrid) ) then
        qty%vGridIndex = vGridIndex
        qty%verticalCoordinate = vGrids(vGridIndex)%verticalCoordinate
        qty%surfs = vGrids(vGridIndex)%surfs
      else
        call SetupEmptyVGridForQuantity ( qty )
      end if
    end if

    ! Setup the cross-track coordinates
    nullify ( qty%crossAngles )
    if ( xGridIndex /= 0 ) then
      call allocate_test ( qty%crossAngles, &
        & size(hGrids(xGridIndex)%the_hGrid%phi,2), 'qty%crossAngles', moduleName )
      qty%crossAngles = hGrids(xGridIndex)%the_hGrid%phi(1,:)
    else
      call allocate_test ( qty%crossAngles, 1, 'qty%crossAngles', moduleName )
      qty%crossAngles = 0.0
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
      end if
      ! print *, '2nd try: fGridIndex ', fGridIndex
      ! print *, 'qty%sharedFGrid ', qty%sharedFGrid
      ! print *, fGrids(fGridIndex)%values
    end if
    if ( properties ( p_sGrid ) ) then ! Faked as frequency, but it's not
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
    qty%name = name
    qty%frequencyCoordinate = frequencyCoordinate
    qty%instrumentmodule = instrumentmodule
    qty%logBasis = logBasis
    qty%minValue = minValue
    qty%molecule = molecule
    qty%quantityType = quantityType
    qty%radiometer = radiometer
    qty%reflector = reflector
    qty%sideband = sideband
    qty%signal = signal
    qty%unit = unitsTable ( quantityType )
    qty%horizontalCoordinate = horizontalCoordinate
    if ( got(f_badValue) ) qty%badValue = badValue

    call deallocate_test ( channels, 'Channels', moduleName )

    call trace_end ( "CreateQtyTemplate", cond=toggle(gen) .and. levels(gen) > 1 )

  end function CreateQtyTemplateFromMLSCFInfo

  ! --------------------------------  ConstructMinorFrameQuantity  -----
  ! I think we should make filedatabase and chunk optional as well
  ! because we don't need it when mifGeolocation is present -haley
  subroutine ConstructMinorFrameQuantity ( instrumentModule, &
    & qty, noChans, regular, instanceLen, NoCrossTrack, &
    & filedatabase, chunk, mifGeolocation )

    use Chunks_m, only: MLSChunk_t
    use Dump_0, only: Dump
    use Dump_1, only: Dump
    use HighOutput, only: BeVerbose, LetsDebug, OutputNamedValue
    use Init_tables_module, only: l_geodaltitude, l_none
    use L1BData, only: L1BData_t, readL1BData, deallocateL1BData, &
      & AssembleL1bQtyName
    use MLSCommon, only: MLSFile_t, nameLen
    use MLSKinds, only: rk => r8
    use MLSFiles, only: HDFVersion_5, Dump, GetMLSFileByType
    use MLSHDF5, only: GetAllHDF5DSNames, IsHDF5DSPresent
    use MLSMessageModule, only: MLSMessage, &
      & MLSMSG_Error, MLSMSG_L1BRead, MLSMSG_Warning
    use MLSSignals_m, only:  Dump_Modules, GetModuleName, &
      & IsAnyModuleSpacecraft, isModuleSpacecraft
    use Output_m, only: Output
    use QuantityTemplates, only: QuantityTemplate_t, &
      & Dump, setupNewQuantityTemplate
    use String_table, only: Display_string
    use Toggles, only: Gen, Levels, Toggle
    use Trace_m, only: Trace_begin, Trace_end

    ! This routine constructs a minor frame based quantity.

    ! Dummy arguments
    integer, intent(in) :: instrumentModule  ! Database index
    type (QuantityTemplate_T), intent(inout) :: qty ! Resulting quantity
    integer, intent(in), optional :: noChans
    logical, intent(in), optional :: regular
    integer, intent(in), optional :: instanceLen
    integer, intent(in), optional :: NoCrossTrack ! Number of cross-track angles
    type (MLSFile_T), dimension(:), pointer, optional ::     FILEDATABASE
    type (MLSChunk_T), intent(in), optional :: chunk   ! The chunk under consideration
    type (QuantityTemplate_T), intent(in), dimension(:), optional :: &
         & mifGeolocation

    ! Local parameters
    real(rk), parameter :: SIXTH = 1.0_rk / 6.0_rk
    ! Note the similarity to CreateHGridFromMLSCFInfo:
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    type :: Item_T
      Integer :: Type
      character(15) :: Name
      logical :: Modular = .true. ! Tangent qty?
      logical :: Bare = .true.  ! Not under any groupname?
    end type Item_T
    enum, bind(C)
      enumerator :: FirstSCItem = 1, scGeodLat = 1, scLon, &
        & MAFStartTimeTAI, tpGeodLat, tpLon, tpGeodAngle, &
        & tpSolarZenith, tpSolarTime, tpLosAngle
    end enum
    integer, parameter :: LastSCItem = scLon
    integer, parameter :: LastL1BItem = tpLosAngle
    ! First, SC items, then instrument items
    type (item_t), parameter :: L1bItemsToRead(firstSCItem:lastL1BItem) = (/ &
      !        Type              Name              Modular
      & item_t(scGeodLat,       "scGeodLat      ", .false., .false.), &
      & item_t(scLon,           "scLon          ", .false., .false.), &
      & item_t(MAFStartTimeTAI, "MAFStartTimeTAI", .false.,  .true.), &
      & item_t(tpGeodLat,       "tpGeodLat      ", .true.,  .false.), &
      & item_t(tpLon,           "tpLon          ", .true.,  .false.), &
      & item_t(tpGeodAngle,     "tpGeodAngle    ", .true.,  .false.), &
      & item_t(tpSolarZenith,   "tpSolarZenith  ", .true.,  .false.), &
      & item_t(tpSolarTime,     "tpSolarTime    ", .true.,  .false.), &
      & item_t(tpLosAngle,      "tpLosAngle     ", .true.,  .false.) /)

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables
    character(len=1024) :: DSNames
    character (len=NameLen) :: instrumentModuleName
    type (L1BData_T) :: l1bField
    type (MLSFile_T), pointer :: L1BFile
    character (len=NameLen) :: l1bItemName
    integer :: hdfVersion, l1bFlag, l1bItem, noMAFs, mafIndex, mifIndex
    integer :: Me = -1      ! For tracing
    logical :: MissingOK
    integer :: Start, Stop  ! For reading L1B quantities, depending upon whether
                            ! they're for the instrument or the tangent point
    logical :: UseMIFGeolocation
    logical :: verbose
    logical :: verboser

    ! Executable code. There are basically two cases here. If we have a
    ! MIFGeolocation argument this conveys all the geolocation for this
    ! quantity.  Otherwise, we have to read it all from the l1boa file
    ! ourselves.

    call trace_begin ( me, "ConstructMinorFrameQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

    MissingOK = .false. ! .not. AURA_L1BFILES
    verbose = BeVerbose( 'qtmp', -1 )
    verboser = LetsDebug( 'qtmp', 0 )
    MissingOK = BeVerbose( 'missOK', -1 )

    call GetModuleName( instrumentModule, instrumentModuleName )
    ! Now due to some error deper than we wish to delve,
    ! The s/c instrument module name is coming out as
    ! 'SC', but the l1boa file thinks it's just 'sc'.
    if ( instrumentModuleName == 'SC' ) instrumentModuleName = 'sc'
    if ( verbose ) then
      if ( qty%name > 0 ) then
        call output ( 'name: ', advance='no' )
        call display_string ( qty%name, advance='yes' )
      endif
      call outputnamedValue( 'instrumentModule index', instrumentModule )
      call outputnamedValue( 'instrumentModule name', trim(instrumentModuleName) )
    endif
    useMIFGeolocation = present(mifGeolocation)
    if ( useMIFGeolocation ) &
      & useMIFGeolocation = mifGeolocation(instrumentModule)%verticalCoordinate == &
                            qty%verticalCoordinate
    if ( verbose ) call outputnamedValue( 'useMIFGeolocation', useMIFGeolocation )
    if ( useMIFGeolocation ) then
      ! --  Got mifGeolocation and vertical coordinates match --
      if ( .not. (present(noChans)) ) &
         call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'You must supply NoChans to reuse geolocation information' )
      if ( verboser ) then
        call output( 'Dump of mifGeolocation; provides lat, etc.', advance = 'yes' )
        call Dump ( mifGeolocation(instrumentModule) )
      endif
      ! We have geolocation information, setup the quantity as a clone of that.
      qty = mifGeolocation(instrumentModule)
      qty%sharedVGrid = .true.
      if ( present ( regular ) ) qty%regular = regular
      qty%noChans = noChans
      if ( qty%regular ) then
        qty%instanceLen = qty%noChans * qty%noSurfs * qty%noCrossTrack
      else
        qty%instanceLen = 0
      end if

    else if (present(chunk) .and. present(filedatabase)) then
      if ( verbose ) call output( 'We have no geolocation information', advance='yes' )
      ! --  Not Got mifGeolocation or vertical coordinates do not match  --
      ! We have no geolocation information, or we don't want to use it.
      ! We have to read it ourselves from the L1BOA file.

      ! First we read xxGeodAlt to get the number of MAFs, and the surfs
      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      hdfversion = L1BFile%HDFVersion
      if ( IsModuleSpacecraft(instrumentModule) ) then
        if ( .not. isAnyModuleSpacecraft() ) then
          call output( 'You should never have come here like this', advance='yes' )
          stop
        endif
        l1bItemName = merge("scGeodAlt","scGeocAlt", &
                           &qty%verticalCoordinate==l_geodAltitude)
      else
        call GetModuleName ( instrumentModule, l1bItemName )
        l1bItemName = TRIM(l1bItemName) // "." // &
          & merge("tpGeodAlt","tpGeocAlt",qty%verticalCoordinate==l_geodAltitude)
      end if
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
      if ( verboser ) then
        call outputnamedValue ( 'l1bItemName#1', trim(l1bItemName) )
        call outputnamedValue ( 'isModuleSpacecraft', isModuleSpacecraft(instrumentModule) )
        call outputnamedValue ( 'isAnyModuleSpacecraft', isAnyModuleSpacecraft() )
        call outputnamedValue ( 'firstMAFIndex', chunk%firstMAFIndex )
        call outputnamedValue ( 'lastMAFIndex', chunk%lastMAFIndex )
        call Dump_Modules
      endif
      if ( verbose ) call outputnamedValue( 'hdfversion', L1BFile%HDFVersion )
      if ( L1BFile%HDFVersion == HDFVersion_5 ) then
        ! OK, so what if GeodAlt isn't here?
        ! 2nd bet: 
        if ( .not. IsHDF5DSPresent( L1BFile, trim(l1bItemName) ) ) &
          & l1bItemName = AssembleL1BQtyName ( 'tpGeodAlt', HDFVersion, &
          & .true., instrumentModuleName )
        if ( .not. IsHDF5DSPresent( L1BFile, trim(l1bItemName) ) ) then

          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & "Will not know number of mafs; unable to find " &
            & // trim(l1bItemName) )
          call Dump ( L1BFile, details=2 )
          ! call Dump ( filedatabase, details=2 )
          call GetAllHDF5DSNames ( L1BFile, DSNames )
          call Dump ( DSNames, 'DSNames' )
          qty%name = 0
          qty%instanceOffset = chunk%firstMAFIndex + chunk%noMAFsLowerOverlap
          qty%grandTotalInstances = 0
          call SetupNewQuantityTemplate ( qty, &
            & noInstances=chunk%lastMAFIndex-chunk%firstMAFIndex + 1, &
            & noSurfs=0, noChans=0, &
            & noCrossTrack=noCrossTrack, coherent=.false., &
            & stacked=.false., regular=regular, instanceLen=0, &
            & minorFrame=.true., verticalCoordinate=qty%verticalCoordinate )
          call outputNamedValue ( 'Set num of MAFs to', &
            & chunk%lastMAFIndex-chunk%firstMAFIndex + 1 )
          call Dump( qty, details=2 )
          call trace_end ( "ConstructMinorFrameQuantity", &
            & cond=toggle(gen) .and. levels(gen) > 2 )
          return
        endif
      endif
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & neverfail=MissingOK )
      if ( l1bFlag /= 0 .and. .not. MissingOK ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
      else if ( l1bFlag /= 0 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        ! return
      end if
      if ( verboser ) then
        call outputnamedValue ( 'noMAFs', noMAFs )
        call outputnamedValue ( 'l1bFlag', l1bFlag )
      endif

      qty%name = 0
      call SetupNewQuantityTemplate ( qty, noInstances=noMAFs, &
        & noSurfs=l1bField%maxMIFs, noChans=noChans, &
        & noCrossTrack=noCrossTrack, coherent=.false., &
        & stacked=.false., regular=regular, instanceLen=instanceLen, &
        & minorFrame=.true., verticalCoordinate=qty%verticalCoordinate )
      qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
      qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap

      ! Now we're going to fill in the hGrid information
      qty%instanceOffset = chunk%firstMAFIndex + chunk%noMAFsLowerOverlap
      qty%grandTotalInstances = 0
      if ( verbose ) then
        call output ( "Instance offset for minor frame quantity is:" )
        call output ( qty%instanceOffset, advance='yes' )
      end if
      
      ! In case we didn't find the proper item in the l1boa files
      if ( l1bFlag /= 0 ) then
        call trace_end ( "ConstructMinorFrameQuantity", &
          & cond=toggle(gen) .and. levels(gen) > 1 )
        return
      endif

      ! Work out which items to read
      if ( IsModuleSpacecraft(instrumentModule) ) then
        start = firstSCitem
        stop = lastSCitem
        ! Zero out stuff that we don't read.  See L1bItemsToRead above.
        qty%phi = 0.0
        qty%time = 0.0
        qty%solarTime = 0.0
        qty%solarZenith = 0.0
        qty%losAngle = 0.0
      else
        start = lastSCitem + 1
        stop = lastL1Bitem
      end if

      ! Now we're going to deal with a VGrid for this quantity.
      ! This is either the scGeocAlt or tpGeo[cd]Alt read above
      qty%surfs = l1bField%dpField(1,:,:)

      call DeallocateL1BData ( l1bField )

      do l1bItem = start, stop
        if ( verbose ) call outputNamedValue ( 'l1bItem num', l1bItem )
        ! Get the name of the item to read
        l1bItemName = L1bItemsToRead(l1bItem)%name
        if ( verbose ) call outputNamedValue ( 'its name', trim(l1bItemName) )
        if ( L1bItemsToRead(l1bItem)%Bare ) then
          ! This should be the name as stored in the l1boa file
        elseif ( L1bItemsToRead(l1bItem)%modular ) then
          if ( verbose ) call output ( 'it is modular', advance='yes' )
          call GetModuleName ( instrumentModule, l1bItemName )
          if ( IsModuleSpacecraft(instrumentModule) ) &
            & call GetModuleName ( instrumentModule+1, instrumentModuleName )
          if ( instrumentModuleName == 'SC' ) instrumentModuleName = 'sc'
          l1bItemName = &
            & trim(instrumentModuleName)//'.'//L1bItemsToRead(l1bItem)%name
          if ( verboser ) &
            & call outputnamedValue ( 'before assembly', trim(l1bItemName) )
          l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
        elseif ( .not. isAnyModuleSpacecraft() ) then
          cycle
        else
          l1bItemName = L1bItemsToRead(l1bItem)%name
          if ( verboser ) &
            & call outputnamedValue ( 'before assembly', trim(l1bItemName) )
          if ( l1bItemName(1:2) == 'sc' ) l1bItemName = l1bItemName(3:)
          l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false., &
            & InstrumentModuleName )
        end if

        ! Read it from the l1boa file
        if ( verboser ) call outputNamedValue ( 'MissingOK', MissingOK )
        call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
          & l1bFlag, firstMAF=chunk%firstMafIndex, &
          & lastMAF=chunk%lastMafIndex, neverfail=MissingOK )
        if ( l1bFlag /= 0 ) then
          call MLSMessage ( merge(MLSMSG_Warning,MLSMSG_Error,missingOK), &
            & ModuleName, MLSMSG_L1BRead//l1bItemName )
          cycle
        end if

        ! Now we have to save this field in the quantity data.
        ! This is rather a kludgy way of doing it but this worked out the
        ! least boring way to write the code.  See the definition of
        ! L1BItemsToRead above for reference.
        
        select case ( L1bItemsToRead(l1bItem)%type )
        case ( MAFStartTimeTAI )
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
        case ( tpGeodLat, scGeodLat )
          qty%geodLat = l1bField%dpField(1,:,:)
        case ( tpLon, scLon )
          qty%lon = l1bField%dpField(1,:,:)
        case ( tpGeodAngle )
          if ( associated(l1bField%intField) ) then
            qty%phi = l1bField%intField(1,:,:)
          else
            qty%phi = l1bField%dpField(1,:,:)
          endif
        case ( tpSolarZenith )
          qty%solarZenith = l1bField%dpField(1,:,:)
        case ( tpSolarTime )
          qty%solarTime = l1bField%dpField(1,:,:)
        case ( tpLosAngle )
          qty%losAngle = l1bField%dpField(1,:,:)
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "No code to read L1B item " // trim(L1bItemsToRead(l1bItem)%name) )
        end select
        if ( verboser ) then
          if ( associated(l1bField%intField) ) then
            call dump( l1bField%intField(1,:,:) )
          else
            call dump( l1bField%dpField(1,:,:) )
          endif
        endif
        call DeallocateL1BData ( l1bField )
      end do                          ! Loop over l1b quantities
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
                        "Must supply either mifGeolocation " &
                        // "or both filedatabase and chunk")
    end if
    qty%frequencyCoordinate = L_None
    qty%instrumentModule = instrumentModule

    ! In later versions we'll probably need to think about FILL_VALUEs and
    ! setting things to the badData flag.
    
    ! We thought about about them, but not too hard apparently.

    call trace_end ( "ConstructMinorFrameQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

  end subroutine ConstructMinorFrameQuantity

  ! ------------------------------------ ForgeMinorFrames --------------
  subroutine ForgeMinorFrames ( root, mifGeolocation )
    ! This routine is used when we're in the mode of creating l2pc files
    ! and we want to invent a set of minor frame quantities with no
    ! reference to the l1 files

    use expr_m, only: expr
    use init_tables_module, only: f_geodalt, f_geodangle, f_module, f_nomifs, &
      & f_solartime, f_solarzenith, f_truncate
    use init_tables_module, only: l_geodaltitude, phyq_angle, phyq_dimensionless, &
      & phyq_time
    use MLSKinds, only: rk => r8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MoreTree, only: Get_Boolean
    use QuantityTemplates, only: CopyQuantityTemplate, &
      & DestroyQuantityTemplateContents, QuantityTemplate_T, &
      & SetupNewQuantityTemplate
    use tree, only: decoration, nsons, subtree
    use vgridsdatabase, only: vgrid_t, vgrids

    ! Dummy arguments
    integer, intent(in) :: ROOT         ! Tree vertex
    type (QuantityTemplate_T), intent(inout), target :: MIFGEOLOCATION(:)

    ! Local variables
    integer :: GEODANGLENODE            ! Tree vertex
    integer :: I                        ! Loop counter
    integer :: INSTRUMENTMODULE         ! Database index
    integer :: KEY                      ! Tree vertex
    integer :: MAF                      ! Loop counter
    integer :: NODE                     ! A tree node
    integer :: NoInstances              ! In a truncated MIF geolocation
    integer :: NOMAFS                   ! Dimension
    integer :: NOMIFS                   ! Dimension
    integer :: NoSurfs                  ! In a truncated MIF geolocation
    integer :: PARAM                    ! Loop counter
    integer :: SOLARTIMENODE            ! Tree vertex
    integer :: SOLARZENITHNODE          ! Tree vertex
    integer :: SON                      ! Tree vertex
    logical :: Truncate                 ! Truncate existing MIF Geolocation
                                        ! instead of filling with zeros.
    integer :: UNITS                    ! Units for node

    real(rk), dimension(:,:), pointer :: VALUES ! An array to fill
    integer, dimension(2) :: EXPR_UNITS ! From tree
    real(rk), dimension(2) :: EXPR_VALUE ! From tree
    type ( quantityTemplate_t ) :: NewMIFGeolocation
    type(VGrid_T), pointer :: GEODALT

    ! Executable code
    nullify ( geodAlt )
    solarTimeNode = 0
    solarZenithNode = 0
    geodAngleNode = 0
    truncate = .false.

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
      case ( f_truncate )
        truncate = get_boolean ( son )
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
    ! Note this will destroy the old one's contents bit by bit.
    if ( truncate ) then
      noSurfs = min(noMIFS, mifGeolocation(instrumentModule)%noSurfs)
      noInstances = min(noMAFS, mifGeolocation(instrumentModule)%noInstances)
      call SetupNewQuantityTemplate ( NewMIFGeolocation, &
        & noInstances=noInstances, noSurfs=noSurfs, noChans=1, &
        & noCrossTrack=1, &
        & coherent=.false., stacked=.false., regular=.true., &
        & minorFrame=.true. )
      ! First put zeros in the new  arrays, in case the old ones are shorter.
      newMIFGeolocation%time = 0.0
      newMIFGeolocation%phi = 0.0
      newMIFGeolocation%geodLat = 0.0
      newMIFGeolocation%solarTime = 0.0
      newMIFGeolocation%solarZenith = 0.0
      newMIFGeolocation%lon = 0.0
      newMIFGeolocation%losAngle = 0.0
      newMIFGeolocation%surfs = 0.0
      newMIFGeolocation%instrumentModule = mifGeolocation(instrumentModule)%instrumentModule
      newMIFGeolocation%noInstancesLowerOverlap = mifGeolocation(instrumentModule)%noInstancesLowerOverlap
      newMIFGeolocation%noInstancesUpperOverlap = mifGeolocation(instrumentModule)%noInstancesUpperOverlap
      newMIFGeolocation%verticalCoordinate = mifGeolocation(instrumentModule)%verticalCoordinate
      newMIFGeolocation%time = mifGeolocation(instrumentModule)%time(:noSurfs,:noInstances)
      newMIFGeolocation%geodLat = mifGeolocation(instrumentModule)%geodLat(:noSurfs,:noInstances)
      newMIFGeolocation%solarTime = mifGeolocation(instrumentModule)%solarTime(:noSurfs,:noInstances)
      newMIFGeolocation%solarZenith = mifGeolocation(instrumentModule)%solarZenith(:noSurfs,:noInstances)
      newMIFGeolocation%lon = mifGeolocation(instrumentModule)%lon(:noSurfs,:noInstances)
      newMIFGeolocation%losAngle = mifGeolocation(instrumentModule)%losAngle(:noSurfs,:noInstances)
      newMIFGeolocation%surfs = mifGeolocation(instrumentModule)%surfs(:noSurfs,:noInstances)
      call copyQuantityTemplate ( mifGeolocation(instrumentModule), newMIFGeolocation )
      call destroyQuantityTemplateContents ( newMIFGeolocation )
  else
      mifGeolocation(instrumentModule)%name = 0
      call SetupNewQuantityTemplate ( mifGeolocation(instrumentModule), &
        & noInstances=noMAFs, noSurfs=noMIFs, noChans=1, &
        & noCrossTrack=1, &
        & coherent=.false., stacked=.false., regular=.true., &
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
    end if

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

    use lexer_core, only: print_source
    use output_m, only: blanks, output
    use tree, only: where_at=>where
    use Intrinsic, only: lit_indices
    use string_table, only: display_string
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    character (LEN=*), intent(in) :: MESSAGE
    integer, intent(in), optional :: EXTRA
    character(len=*), intent(in), optional :: SEVERITY ! 'nonfatal' or 'fatal'
    character(len=8) :: mySeverity

    mySeverity = 'fatal'
    if ( present(severity) ) then
      if (index('WwNn', severity(1:1)) > 0 ) mySeverity='warning'
    end if
    if ( mySeverity /= 'fatal' ) then
      call blanks(5)
      call output ( ' (warning) ' )
    else
      call blanks(5, fillChar='*')
      call output ( ' (fatal) ' )
    end if
    call output ( 'At ' )
    if ( where > 0 ) then
      call print_source ( where_at(where) )
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
  function AnyGoodSignalData ( SIGNAL, SIDEBAND, FILEDATABASE, CHUNK, MAFWISE ) &
    &  result (ANSWER)
  ! Read precision of signal
  ! if all values < 0.0, return FALSE
  ! if no precision data in file, return FALSE
  ! otherwise return true
  ! Arguments

    use chunks_m, only: MLSChunk_t
    use allocate_deallocate, only: deallocate_test
    use L1BData, only: L1BData_t, readL1BData, getl1bfile, &
      & assemblel1bqtyname, precisionsuffix
    use MLSCommon, only: MLSFile_t
    use MLSKinds, only: rk => r8
    use MLSFiles, only: geTMLSFileByType
    use MLSSignals_m, only: getSignalName

    integer, intent(in)                         :: signal
    integer, intent(in)                         :: sideband
    logical                                     :: answer
    type (MLSChunk_T), intent(in)               :: Chunk
    type (MLSFile_T), dimension(:), pointer     :: FILEDATABASE
    logical, dimension(:), optional, intent(out):: MAFWISE
  ! Private
    ! integer :: FileID
    integer :: flag
    integer :: i
    integer :: noMAFs
    character(len=127)  :: namestring
    type (l1bData_T) :: MY_L1BDATA
    type (MLSFile_T), pointer ::     L1BFile

  ! Executable
  answer = .false.
  if ( present(mafwise) ) mafwise = .false.
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) return
    call GetSignalName ( signal, nameString, &                   
    & sideband=sideband, noChannels=.TRUE. )                     
    nameString = AssembleL1BQtyName ( nameString, L1BFile%hdfVersion, .false. )
    nameString = trim(nameString) // PRECISIONSUFFIX
    L1BFile => GetL1bFile(filedatabase, namestring)
    ! fileID = FindL1BData ( filedatabase, nameString )
    if ( .not. associated(L1BFile) ) then
      answer = .false.
      return
    end if
    ! print *, 'About to read ', trim(nameString)
    ! print *, 'From Fileid ', fileID
    call ReadL1BData ( L1BFile, nameString, my_l1bData, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & NeverFail= .true. )
    if ( flag == 0 ) then
      answer = .not. all (my_l1bData%DpField < 0._rk)
      if ( present( mafwise) ) then
        do i = 1, chunk%lastMAFIndex - chunk%firstMAFIndex + 1
          mafwise(i) = .not. all (my_l1bData%DpField(:,:,i) < 0._rk)
        enddo
      endif
      call deallocate_test(my_l1bData%DpField, trim(nameString), ModuleName)
    else
      answer = .false.
    end if
  end function AnyGoodSignalData

  ! --------------------------------  ConstructMajorFrameQuantity  -----
  subroutine ConstructMajorFrameQuantity( chunk, instrumentModule, qty, noChans, &
    & mifGeolocation, NoCrossTrack )
    ! Dummy arguments
    use chunks_m, only: MLSChunk_t
    use quantityTemplates, only: CreateGeolocationFields, quantityTemplate_t, &
      & setUpNewQuantityTemplate
    use toggles, only: gen, levels, toggle
    use trace_m, only: trace_begin, trace_end

    type (MLSChunk_T), intent(in) :: CHUNK
    integer, intent(in) :: INSTRUMENTMODULE
    type (QuantityTemplate_T), intent(out) :: QTY
    type (QuantityTemplate_T), dimension(:), intent(in), target :: MIFGEOLOCATION
    integer, intent(in) :: NOCHANS
    integer, intent(in), optional :: NoCrossTrack
    ! Local variables

    integer :: Me = -1     ! For tracing
    type (QuantityTemplate_T), pointer :: source
    
    ! Executable code
    call trace_begin ( me, "ConstructMajorFrameQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 2 )

    source => mifGeolocation(instrumentModule)
    qty%name = 0
    call SetupNewQuantityTemplate ( qty, noInstances=source%noInstances, &
      & noSurfs=1, coherent=.true., stacked=.true., regular=.true., &
      & noChans=noChans, noCrossTrack=noCrossTrack, &
      & sharedVGrid=.true. )
    call SetupEmptyVGridForQuantity ( qty )
    ! In some rare cases, e.g. non-Aura satellite data, source might be a dud
    if ( allocated(source%phi) ) then
      call createGeolocationFields ( qty, qty%noSurfs, 'Qty' )
      qty%geodLat = source%geodLat(1:1,:)
      qty%lon = source%lon(1:1,:)
      qty%phi = source%phi(1:1,:)
      qty%time => source%time(1:1,:)
      qty%solarTime => source%solarTime(1:1,:)
      qty%solarZenith => source%solarZenith(1:1,:)
      qty%losAngle => source%losAngle(1:1,:)
    end if

    qty%majorFrame = .true.
    qty%minorFrame = .false.
    qty%instanceOffset = source%instanceOffset
    qty%instrumentModule = source%instrumentModule
    qty%noInstancesLowerOverlap = chunk%noMAFsLowerOverlap
    qty%noInstancesUpperOverlap = chunk%noMAFsUpperOverlap

    call trace_end ( "ConstructMajorFrameQuantity", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

  end subroutine ConstructMajorFrameQuantity

  ! ----------------------------------------------  GetQtyTypeIndex  -----
  subroutine GetQtyTypeIndex(string_text, QtyType)
    use Intrinsic, only: lit_indices
    use MLSStrings, only: lowercase
    use string_table, only: get_string
    ! Returns the lit index,  given QtyType name in mixed case
    ! Returns 0 if QtyType name not found
    ! (inverse function: GetQtyTypeName)
    integer, intent(out) :: QtyType
    character (len=*), intent(in) :: string_text
    ! Local variables
    character (len=len(string_text))             :: string_test
    do QtyType=first_lit, last_lit
      ! Skip if not quantity type
      if ( .not. any(PROPERTYTABLE(:, QtyType) ) ) cycle
      call get_string ( lit_indices(QtyType), string_test, strip=.true. )
      if ( LowerCase(trim(string_text)) == LowerCase(trim(string_test))) &
       & return
    end do
    QtyType = 0
  end subroutine GetQtyTypeIndex

  ! ----------------------------------------------- InitQuantityTemplates ----
  subroutine InitQuantityTemplates
    ! This routine initializes the quantity template properties
    ! This is the routine one needs to update when one introduces a new quantity type.
    use Init_Tables_Module, only:  l_adopted, l_adopted, l_baseline, &
      l_boundarypressure, l_calsidebandfraction, &
      l_chisqbinned, l_chisqchan, l_chisqmmaf, l_chisqmmif, l_cloudice, &
      l_cloudinducedradiance, l_cloudextinction, l_cloudminmax, l_cloudradsensitivity, &
      l_cloudtemperature, l_CloudTopPressure, l_cloudwater, l_columnabundance, &
      l_dnwt_abandoned, l_dnwt_ajn, l_dnwt_axmax, l_dnwt_cait, &
      l_dnwt_chisqminnorm, l_dnwt_chisqnorm, l_dnwt_chisqratio, &
      l_dnwt_count, l_dnwt_diag, l_dnwt_dxdx, l_dnwt_dxdxl, &
      l_dnwt_dxn, l_dnwt_dxnl, l_dnwt_flag, l_dnwt_fnmin, &
      l_dnwt_fnorm, l_dnwt_gdx, l_dnwt_gfac, &
      l_dnwt_gradn, l_dnwt_sq, l_dnwt_sqt,&
      l_earthradius, l_earthrefl, l_ecrtofov, l_effectiveopticaldepth, &
      l_elevoffset, l_extinction, l_extinctionv2, &
      l_fieldazimuth, l_fieldelevation, l_fieldstrength, &
      l_geolocation, l_gph, l_ghzazim, l_heightoffset, l_isotoperatio, l_iwc, &
      l_jacobian_cols, l_jacobian_rows, &
      l_l1bmafbaseline, l_l1bmif_tai, l_limbsidebandfraction, &
      l_linecenter, l_linewidth, l_linewidth_tdep, &
      l_lostransfunc, l_losvel, l_lowestretrievedpressure, &
      l_massmeandiameterice, l_massmeandiameterwater, l_magneticfield, &
      l_mifdeadtime, l_mifextinction, l_mifextinctionextrapolation, &
      l_mifExtinctionform, l_mifExtinctionv2, &
      l_mifRadC, l_mifRHI, &
      l_noisebandwidth, l_noradspermif, l_noradsbinned, &
      l_numgrad, l_numj, l_numnewt, &
      l_opticaldepth, l_orbitinclination, l_ascdescmode, l_geoheight, &
      l_phasetiming, l_phitan, l_ptan, l_quality, l_radiance, &
      l_refgph, l_refltemp, l_refltrans, l_reflrefl, l_reflspill, &
      l_rhi, l_singlechannelradiance, l_sizedistribution, &
      l_scanresidual, l_scatteringangle, l_sceci, l_instecr, l_scecr, l_scveleci, &
      l_scvelecr, l_scgeocalt, l_spaceradiance, l_status, &
      l_strayradiance, l_surfaceheight, l_surfacetype, l_systemtemperature, &
      l_temperature, l_tngteci, l_tngtecr, l_tngtgeodalt, l_tngtgeocalt, &
      l_totalpowerweight, l_tscat, l_vmr
    use Init_Tables_Module, only:  phyq_angle, phyq_colmabundance, &
     & phyq_dimensionless, phyq_extinction, phyq_frequency,&
     & phyq_gauss, phyq_icedensity, phyq_length, &
     & phyq_pressure, phyq_temperature, phyq_time, phyq_velocity, &
     & phyq_vmr, phyq_zeta

    use MLSMessageModule, only: MLSMSG_Error, MLSMessage
    use Intrinsic, only: Lit_Indices
    use Output_M, only: Output
    use String_Table, only: Display_String

    ! Local variables
    integer :: I                        ! Loop counter
    integer, dimension(1:0), parameter :: None = (/ ( 0, i=1,0 ) /)
    logical :: Valid                    ! Flag
    character(len=132) :: Message       ! An error message

    ! Executable code
    ! Basically here, we're going to go through and populate the various tables

    propertyTable = .false.
    unitsTable = 0

    ! DefineQtyTypes creates the arrays PropertyTable and UnitsTable. For
    ! each quantity type, the first thing in the list is the quantity's lit
    ! index.  The next is its physical dimension.  After that there is a
    ! list of properties.  The list of properties for a quantity ends with
    ! either "next,", after which information for another quantity may
    ! begin, or the end of the array.  The property "none" is a zero-extent
    ! array, provided for its documentary value.

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
      l_cloudTemperature, phyq_temperature, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_CloudTopPressure, phyq_pressure, p_hGrid, next, &
      l_cloudWater, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_columnAbundance, phyq_colmabundance, p_hGrid, p_molecule, next, &
      l_dnwt_abandoned, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_ajn, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_axmax, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_cait, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_chisqminnorm, phyq_dimensionless, p_vGrid /) )

    call DefineQtyTypes ( (/ &
      l_dnwt_chisqnorm, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_chisqratio, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_count, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_diag, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdxl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxdx, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxnl, phyq_dimensionless, p_vGrid, next, &
      l_dnwt_dxn, phyq_dimensionless, p_vGrid, next, &
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
      l_extinctionv2, phyq_extinction, p_hGrid, p_vGrid, p_fGrid, p_radiometer, &
                      p_mustBeZeta, next, &
      l_fieldAzimuth, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldElevation, phyq_angle, p_hGrid, p_vGrid, next, &
      l_fieldStrength, phyq_gauss, p_hGrid, p_vGrid, next, &
      l_geolocation, phyq_dimensionless, p_majorframe, p_module, next, &
      l_gph, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_GHzAzim, phyq_angle, p_minorFrame, p_module, next, &
      l_heightOffset, phyq_length, p_hGrid, p_vGrid, next, &
      l_isotopeRatio, phyq_dimensionless, p_molecule, next, &
      l_iwc, phyq_icedensity, p_hGrid, p_vGrid, next, &
      l_jacobian_cols, phyq_dimensionless, p_vGrid, next, &
      l_jacobian_rows, phyq_dimensionless, p_vGrid, next /) )

    call DefineQtyTypes ( (/ & 
      l_l1bMAFBaseline, phyq_temperature, p_majorFrame, p_signal, next, &
      l_l1bMIF_TAI, phyq_time, p_minorFrame, p_scmodule, next, &
      l_limbSidebandFraction, phyq_dimensionless, p_signal, next, &
      l_lineCenter, phyq_frequency, p_hGrid, p_vGrid, p_molecule, next, &
      l_lineWidth,  phyq_frequency, p_hGrid, p_vGrid, p_molecule, next, &
      l_lineWidth_TDep, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, next, &
      l_losTransFunc, phyq_dimensionless, p_minorFrame, p_sGrid, p_module, next, &
      l_losVel, phyq_dimensionless, p_minorFrame, p_module, next, & ! ??? Really phyq_dimensionless
      l_lowestRetrievedPressure, phyq_zeta, none, next, &
      l_magneticField, phyq_gauss, p_vGrid, p_hGrid, p_xyz, next, &
      l_massMeanDiameterIce, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, & ! ??? Really phyq_dimensionless
      l_massMeanDiameterWater, phyq_dimensionless, p_vGrid, p_hGrid, p_mustBeZeta, next, & ! ??? Really phyq_dimensionless
      l_mifDeadTime, phyq_time, next, &
      l_mifExtinction, phyq_extinction, p_flexibleVHGrid, &
        & p_minorFrame, p_radiometer, p_mustBeZeta, next, &
      l_mifExtinctionExtrapolation, phyq_dimensionless, none, next, &
      l_mifExtinctionForm, phyq_dimensionless, none, next, &
      l_mifExtinctionV2, phyq_extinction, p_flexibleVHGrid, &
        & p_minorFrame, p_radiometer, p_mustBeZeta, next, &
      phyq_length, p_minorFrame, p_module, next, &
      l_mifRadC, phyq_length, p_minorFrame, p_module, next, &
      l_mifRHI, phyq_dimensionless, p_flexibleVHGrid, &
        & p_minorFrame, p_radiometer, p_mustBeZeta, next, &
      l_noiseBandwidth, phyq_frequency, p_signal, next, &
      l_noRadsBinned, phyq_dimensionless, p_vGrid, p_hGrid, &
        & p_signal, p_suppressChannels, p_mustBeZeta, next, &
      l_noRadsPerMIF, phyq_dimensionless, p_minorFrame, p_signal, &
        & p_suppressChannels, next, &
      l_numNewt, phyq_dimensionless, p_vGrid, next, &
      l_numGrad, phyq_dimensionless, p_vGrid, next, &
      l_numJ, phyq_dimensionless, p_vGrid, next /) )

    call DefineQtyTypes ( (/ &
      l_opticalDepth, phyq_dimensionless, p_minorFrame, p_signal, next, &
      l_orbitInclination, phyq_angle, p_minorFrame, p_scModule, next, &
      l_AscDescMode, phyq_dimensionless, p_hGrid, next, & 
      l_geoheight, phyq_dimensionless, p_hGrid, p_vGrid, next, &
      l_phaseTiming, phyq_dimensionless, p_vGrid, next, &
      l_phiTan, phyq_angle, p_minorFrame, p_module, next, &
      l_ptan, phyq_zeta, p_minorFrame, p_module, next, &
      l_quality, phyq_dimensionless, p_hGrid, next, &
      l_radiance, phyq_temperature, p_minorFrame, p_signal, next, & 
      l_refGPH, phyq_length, p_hGrid, p_vGrid, p_mustBeZeta, next, &
      l_reflrefl, phyq_dimensionless, p_signal, p_reflector, next, &
      l_reflspill, phyq_temperature, p_signal, p_majorframe, p_reflector, next, &
      l_refltemp, phyq_temperature, p_majorFrame, p_reflector, p_module, next, &
      l_refltrans, phyq_dimensionless, p_signal, p_reflector, next, &
      l_rhi, phyq_dimensionless, p_hGrid, p_vGrid, p_molecule, p_mustBeZeta, next, &
      l_scanResidual, phyq_length, p_minorFrame, p_module, next, &
      l_scatteringAngle, phyq_angle, p_vGrid, next, &
      l_scECI, phyq_length, p_minorFrame, p_scModule, p_xyz, next, &
      l_instECR, phyq_length, p_minorFrame, p_Module, p_xyz, next, &
      l_scECR, phyq_length, p_minorFrame, p_scModule, p_xyz, next, &
      l_scGeocAlt, phyq_length, p_minorFrame, p_scModule, next, &
      l_scVelECI, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_scVelECR, phyq_velocity, p_minorFrame, p_scModule, p_xyz, next, &
      l_singleChannelRadiance, phyq_temperature, p_minorFrame, p_signal, &
                               p_suppressChannels, next, &
      l_sizeDistribution, phyq_dimensionless, p_hGrid, p_vGrid, p_mustBeZeta, next, & 
      l_spaceRadiance, phyq_temperature, none, next, &
      l_status, phyq_dimensionless, p_hGrid, next, &
      l_strayRadiance, phyq_temperature, p_signal, p_majorFrame, next, &
      l_surfaceHeight, phyq_length, p_hGrid, next, &
      l_surfaceType, phyq_dimensionless, p_hGrid, next, & 
      l_systemTemperature, phyq_temperature, p_signal, next, &
      l_temperature, phyq_temperature, p_hGrid, p_vGrid, p_mustbezeta, next, &
      l_tngtECI, phyq_length, p_minorFrame, p_module, p_xyz, next, &
      l_tngtECR, phyq_length, p_minorFrame, p_module, p_xyz, next, &
      l_tngtGeocAlt, phyq_length, p_minorFrame, p_module, next, &
      l_tngtGeodAlt, phyq_length, p_minorFrame, p_module, next, &
      l_totalPowerWeight, phyq_dimensionless, p_signal, next, &
      l_TScat, phyq_temperature, p_hGrid, p_signal, p_vGrid, next, &
      l_vmr, phyq_vmr, p_hGrid, p_vGrid, p_fGridOptional, &
             p_molecule, p_radiometerOptional, p_mustbezeta, next /) )

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
      ! Check that mustBeZeta or mustBeGeocAlt quantities have a vGrid
      if ( ( propertyTable ( p_mustBeZeta, i ) .or. &
           & propertyTable ( p_mustBeGeocAlt, i ) ).and. .not. &
        & ( propertyTable ( p_vGrid, i ) .or. propertyTable ( p_flexibleVHGrid, i ) ) ) then
        valid = .false.
        message = "Quantity must have vGrid if it must be on log-pressure or geocentric altitude"
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
      do while ( i <= size(info) )
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
      end do
    end subroutine DefineQtyTypes
    
  end subroutine InitQuantityTemplates

  ! ---------------------------------- SetupEmptyHGridForQuantity
  subroutine SetupEmptyHGridForQuantity ( qty ) 
    use Allocate_Deallocate, only: allocate_test
    use QuantityTemplates, only: QuantityTemplate_T
    ! Dummy arguments
    type ( QuantityTemplate_T ), intent(inout) :: QTY
    ! Executable code
    qty%hGridIndex = 0
    qty%xGridIndex = 0
    ! Lest we deallocate a database entry:
    nullify( qty%frequencies, qty%time, qty%solarTime, &
      & qty%solarZenith, qty%losAngle, qty%crossAngles, qty%the_hGrid )
    call allocate_test ( qty%phi, 1, 1, 'qty%phi(1,1)', ModuleName )
    call allocate_test ( qty%geodLat, 1, 1, 'qty%geodLat(1,1)', ModuleName )
    call allocate_test ( qty%lon, 1, 1, 'qty%lon1(1,1)', ModuleName )
    call allocate_test ( qty%time, 1, 1, 'qty%time(1,1)', ModuleName )
    call allocate_test ( qty%solarTime, 1, 1, 'qty%solarTime(1,1)', ModuleName )
    call allocate_test ( qty%solarZenith, 1, 1, 'qty%solarZenith(1,1)', ModuleName )
    call allocate_test ( qty%losAngle, 1, 1, 'qty%losAngle(1,1)', ModuleName )
    call allocate_test ( qty%crossAngles, 1, 'qty%crossAngles(1)', ModuleName )
    qty%phi = 0.0
    qty%geodLat = 0.0
    qty%lon = 0.0
    qty%time = 0.0
    qty%solarTime = 0.0
    qty%solarZenith = 0.0
    qty%losAngle = 0.0
    qty%crossAngles = 0.0
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
    call Allocate_test ( qty%surfs, 1, 1, 'qty%surfs(1,1)', ModuleName )
    qty%surfs = 0. ! We used to have impossible values for bnd. prs.
  end subroutine SetupEmptyVGridForQuantity

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ConstructQuantityTemplates
!
! $Log$
! Revision 2.209  2024/02/02 21:31:45  pwagner
! Added CloudTopPressure as new Quantity type
!
! Revision 2.208  2020/01/27 18:02:46  pwagner
! Worked around the error that made instrumentModuleName allcaps
!
! Revision 2.207  2018/04/11 17:47:05  pwagner
! Should work better with ASMLS mif geolocations
!
! Revision 2.206  2018/02/23 22:04:45  mmadatya
! Added l_instECR for ASMLS
!
! Revision 2.205  2017/12/15 18:43:56  mmadatya
! Added geoHeight as new quantity template
!
! Revision 2.204  2017/10/27 01:29:46  vsnyder
! Remove L_MIFLOS
!
! Revision 2.203  2017/06/01 22:49:15  vsnyder
! Remove unused USE names
!
! Revision 2.202  2016/11/04 23:05:48  pwagner
! Define some potentially undefined qty components
!
! Revision 2.201  2016/10/20 23:13:39  pwagner
! When warning of missing GeodAlt, say how big we set numMAFs
!
! Revision 2.200  2016/10/14 00:05:46  pwagner
! Show qty name, too, when verbosely Constructing its template
!
! Revision 2.199  2016/09/21 00:41:32  pwagner
! Default to printing less
!
! Revision 2.198  2016/08/12 00:33:05  pwagner
! Seems to restore tthe gold brick
!
! Revision 2.197  2016/08/09 21:51:57  pwagner
! Survives encounter with non-satellite data
!
! Revision 2.196  2016/07/28 19:56:00  pwagner
! Extra care in CConstructMinorFrameQs
!
! Revision 2.195  2016/07/28 00:42:46  vsnyder
! Remove unused USE
!
! Revision 2.194  2016/07/27 23:03:44  pwagner
! Works better with Aircraft-borne instrument data
!
! Revision 2.193  2016/07/25 23:41:53  pwagner
! Avoid reading l1b s/c fields if no module is s/c
!
! Revision 2.192  2016/07/21 20:28:54  pwagner
! Now remembers to call trace_end before returning from a failed l1b read
!
! Revision 2.191  2016/05/27 20:59:04  vsnyder
! Remove l_mifIncline and l_mifTanHt because they are redundant to
! l_OrbitInclination and l_tngtGeodAlt
!
! Revision 2.190  2016/05/27 01:24:07  vsnyder
! Add L_mifTanHt
!
! Revision 2.189  2016/05/24 01:26:01  vsnyder
! Add mifIncline, mifLOS, mifRadC
!
! Revision 2.188  2016/05/18 01:37:30  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.187  2016/05/04 18:31:39  pwagner
! Requires -Sqtmp1 to dump l1b arrays when constructing minor fram qty
!
! Revision 2.186  2015/09/25 02:16:20  vsnyder
! Preserve the quantity template vertical coordinate.  Read the appropriate
! kind of altitude quantity from the L1BOA file, depending on the quantity
! template's vertical coordinate.
!
! Revision 2.185  2015/09/22 01:57:57  vsnyder
! Add GHzAzim and ScECR quantity types
!
! Revision 2.184  2015/09/17 23:21:08  pwagner
! Turning on verbose(r) prints more now in ConstructMinorFrameQuantity
!
! Revision 2.183  2015/07/29 00:29:54  vsnyder
! Convert Phi from pointer to allocatable
!
! Revision 2.182  2015/07/27 22:28:06  vsnyder
! Make cross angles always associated, with value zero if there's no xGrid
!
! Revision 2.181  2015/06/04 03:14:14  vsnyder
! Make Surfs component of quantity template allocatable
!
! Revision 2.180  2015/06/03 00:01:51  vsnyder
! Allocate rank-2 latitude and longitude, instead of rank-1 ones and then
! remapping to rank-2 and rank-3.
!
! Revision 2.179  2015/05/28 18:25:16  vsnyder
! Eliminate shared HGrid, get PointQuantityToHGrid from QuantityTemplates
!
! Revision 2.178  2015/03/28 02:29:39  vsnyder
! Added P_MustBeGeocalt, Coordinate, NoCrossTrack, Stacked, xGridIndex,
! xGrid, noCrossTrack, tracing, scGeodLat, scLon.  Support cross-track grids.
! Use Get_Boolean from More_Tree.  Added truncate field to Forge command --
! maybe should have been called "save" -- to save any geolocations gotten
! from L1BOA.  Deleted Azimuth.  Deleted mustBeZeta requirement of
! magnetic field (because we now use height for the vertical coordinate).
! Support for 3-d GeodLat and Lon.
!
! Revision 2.177  2014/09/05 00:39:49  vsnyder
! Add some tracing
!
! Revision 2.176  2014/04/24 23:58:21  pwagner
! If hGrid supplied, use it to set horizontalCoordinate
!
! Revision 2.175  2014/04/07 18:01:12  pwagner
! Added new quantity type AscDescMode
!
! Revision 2.174  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.173  2013/08/17 00:23:35  pwagner
! New cmdline arg relaxes some for non-Aura l1b datasets
!
! Revision 2.172  2013/07/25 00:23:00  vsnyder
! Add MIFRHI quantity type
!
! Revision 2.171  2013/07/19 01:21:49  vsnyder
! Sort some stuff, turn off the first time flag so that InitQuantityTemplates
! is not done on every call.
!
! Revision 2.170  2013/07/18 01:12:11  vsnyder
! Remove scVel since it's ambiguous whether it's ECI or ECR, and nobody
! uses it anyway.
!
! Revision 2.169  2013/05/31 00:42:35  vsnyder
! Get geolocation if requested for L1B fill
!
! Revision 2.168  2013/05/21 23:52:47  vsnyder
! Add MIFExtinctionExtrapolation and MIFExtinctionForm
!
! Revision 2.167  2013/02/21 21:37:10  pwagner
! New optional arg mafwise to return anysignaldata maf-by-maf
!
! Revision 2.166  2012/01/05 01:20:47  pwagner
! Capitalized USEd stuff
!
! Revision 2.165  2011/12/21 01:40:57  vsnyder
! Add LowestRetrievedPressure, MIFExtinction[v2], and a description of the
! DefineQtyTypes table.
!
! Revision 2.164  2011/08/20 00:49:05  vsnyder
! Remove unused use names
!
! Revision 2.163  2011/05/09 18:05:32  pwagner
! Converted to using switchDetail
!
! Revision 2.162  2010/04/05 17:40:55  honghanh
! Fixed bug introduced by making filedatabase optional
!
! Revision 2.160  2010/03/31 19:59:55  honghanh
! Removing l_geodAltitude quantityType as there is l_tngtgeodaltitude
! quantityType already
!
! Revision 2.159  2010/03/24 17:35:13  honghanh
! Just a note to fix ConstructMinorFrameQuantity
!
! Revision 2.158  2010/03/09 22:44:08  honghanh
! Remove cfm_CreateQtyTemplate, set InitQuantityTemplate, and previously private constants to public
!
! Revision 2.157  2010/03/09 22:40:50  honghanh
! Change to support cfm creation of quantity template
!
! Revision 2.156  2010/03/05 20:17:51  honghanh
! Put ConstructMajorFrameQuantity to the list of public subroutine
!
! Revision 2.155  2010/03/02 01:09:20  pwagner
! Added geodAltitude as a quantity type
!
! Revision 2.154  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.153  2010/01/23 01:02:37  vsnyder
! Remove LogIWC
!
! Revision 2.152  2009/09/25 02:40:28  vsnyder
! Add badValue to specify badValue field, keepChannels to allocate a
! channels(:) pointer and save the channels output from Parse_Signal
!
! Revision 2.151  2009/09/22 17:02:31  pwagner
! NAG complained of too many continuation statements
!
! Revision 2.150  2009/09/19 00:33:44  vsnyder
! Add LogIWC
!
! Revision 2.149  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.148  2009/03/14 02:44:46  honghanh
! Add dnwt_count and dnwt_abandoned
!
! Revision 2.147  2008/10/03 16:39:42  livesey
! Added EXTINCTIONV2
!
! Revision 2.146  2008/09/30 22:34:19  vsnyder
! Remove AuxGrid and an unnecessary dump routine
!
! Revision 2.145  2008/08/22 01:04:18  vsnyder
! Don't allow multiple AuxGrids fields
!
! Revision 2.144  2008/06/06 01:58:31  vsnyder
! Make Aux grids VGrids, define ScatteringAngle grid type
!
! Revision 2.143  2008/06/05 20:00:45  vsnyder
! auxGrid never required
!
! Revision 2.142  2008/06/05 02:14:29  vsnyder
! Added AuxGrid field to Quantity.  Defined CloudTemperature and TScat types.
!
! Revision 2.141  2008/05/28 21:03:58  pwagner
! New quantity type to hold chunk number[maf]
!
! Revision 2.140  2008/04/26 00:39:46  livesey
! Added total power stuff
!
! Revision 2.139  2007/03/08 01:33:33  pwagner
! We should never use undefined values for surfs even if no vgrid
!
! Revision 2.138  2007/01/24 02:17:29  vsnyder
! Add TARGET attribute for MifGeolocation to prevent dangling pointer
!
! Revision 2.137  2007/01/11 20:44:27  vsnyder
! Add SurfaceHeight
!
! Revision 2.136  2006/10/02 23:05:31  pwagner
! May Fill chi^2 ratio to measure convergence
!
! Revision 2.135  2006/08/11 20:36:24  vsnyder
! Add INTENT(IN) to MIFGeolocation
!
! Revision 2.134  2006/08/04 20:52:20  pwagner
! Restore quantity name (if clobbered by SetUp..)
!
! Revision 2.133  2006/08/03 01:57:42  vsnyder
! Make sure qty%name is defined
!
! Revision 2.132  2006/07/20 23:39:53  vsnyder
! Remove unused declarations and USEs
!
! Revision 2.131  2006/06/01 03:06:46  vsnyder
! Define numGrad and numNewt
!
! Revision 2.130  2006/04/11 23:32:08  pwagner
! Fixed bug which added excess profiles
!
! Revision 2.129  2006/03/04 00:18:02  pwagner
! AnyGoodSignalData made public
!
! Revision 2.128  2006/01/11 18:00:12  pwagner
! Consistent with new abstract phys quant colmAbundance
!
! Revision 2.127  2005/09/14 00:13:30  pwagner
! Uses ChunkDivideConfig%allowPriorOverlaps in calculating qty%noInstancesLowerOverlap
!
! Revision 2.126  2005/09/02 21:57:23  vsnyder
! Add spectroscopy parameter quantities
!
! Revision 2.125  2005/08/09 00:03:04  pwagner
! hdfVersion not left undefined in AnyGoodSignalData
!
! Revision 2.124  2005/08/04 02:59:54  vsnyder
! Correct definitions for L1BMIF_TAI and MIFDeadTime
!
! Revision 2.123  2005/08/03 18:08:35  vsnyder
! Add L1BMIF_TAI and MifDeadTime for scan averaging
!
! Revision 2.122  2005/06/03 02:05:29  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of passing as an argument.
!
! Revision 2.121  2005/06/01 17:39:26  pwagner
! Dont read L1bFile if unassocated
!
! Revision 2.120  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
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
