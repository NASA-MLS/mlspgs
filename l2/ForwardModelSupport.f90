! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelSupport

  ! Set up the forward model stuff.

  implicit none
  private
  public :: ConstructForwardModelConfig, ForwardModelGlobalSetup, &
    & CreateBinSelectorFromMLSCFInfo, printForwardModelTiming, &
    & resetForwardModelTiming

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Error codes

  integer, parameter :: AllocateError        = 1
  integer, parameter :: BadBinSelectors      = AllocateError + 1
  integer, parameter :: BadMolecule          = BadBinSelectors + 1
  integer, parameter :: DefineSignalsFirst   = BadMolecule + 1
  integer, parameter :: DefineMoleculesFirst = DefineSignalsFirst + 1
  integer, parameter :: IncompleteBinSelectors    = DefineMoleculesFirst + 1
  integer, parameter :: IncompleteFullFwm    = IncompleteBinSelectors + 1
  integer, parameter :: IncompleteLinearFwm  = IncompleteFullFwm + 1
  integer, parameter :: IrrelevantFwmParameter = IncompleteLinearFwm + 1
  integer, parameter :: LinearSidebandHasUnits = IrrelevantFwmParameter + 1
  integer, parameter :: TangentNotSubset     = LinearSidebandHasUnits + 1
  integer, parameter :: ToleranceNotK        = TangentNotSubset + 1
  integer, parameter :: TooManyHeights       = ToleranceNotK + 1
  integer, parameter :: TooManyCosts         = TooManyHeights + 1
  integer, parameter :: BadHeightUnit        = TooManyCosts + 1
  integer, parameter :: NoMolecule           = BadHeightUnit + 1
  integer, parameter :: BadQuantityType      = NoMolecule + 1
  integer, parameter :: WrongUnitsForWindow  = BadQuantityType + 1

  integer :: Error            ! Error level -- 0 = OK

contains ! =====     Public Procedures     =============================

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root, any_errors )
    ! Process the forwardModel specification to produce ForwardModelInfo.

    use AntennaPatterns_m, only: OPEN_ANTENNA_PATTERNS_FILE, &
      & READ_ANTENNA_PATTERNS_FILE, CLOSE_ANTENNA_PATTERNS_FILE
    use FilterShapes_m, only: OPEN_FILTER_SHAPES_FILE, &
      & READ_FILTER_SHAPES_FILE, READ_DACS_FILTER_SHAPES_FILE, &
      & CLOSE_FILTER_SHAPES_FILE
    use Init_Tables_Module, only: F_ANTENNAPATTERNS, F_DACSFILTERSHAPES, &
      & F_FILTERSHAPES, F_L2PC, F_POINTINGGRIDS
    use L2ParInfo, only: PARALLEL
    use L2PC_m, only: OPEN_L2PC_FILE, CLOSE_L2PC_FILE, READ_L2PC_FILE, &
      & READCOMPLETEHDF5L2PCFILE
    use MLSFiles, only: GetPCFromRef, split_path_name
    use MLSL2Options, only: TOOLKIT
    use MLSPCF2, only: MLSPCF_antpats_start, MLSPCF_filtshps_start, &
      &          mlspcf_dacsfltsh_start, MLSPCF_ptggrids_start, &
      &          mlspcf_l2pc_start, mlspcf_l2pc_end
    use MoreTree, only: Get_Field_ID
    use PointingGrid_m, only: Close_Pointing_Grid_File, &
      & Open_Pointing_Grid_File, Read_Pointing_Grid_File
    use String_Table, only: Get_String
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use Tree, only: Nsons, Sub_Rosa, Subtree

    integer, intent(in) :: Root         ! of the forwardModel specification.
    !                                     Indexes a "spec_args" vertex.
    integer, intent(out) :: any_errors  ! non-zero means trouble

    logical, parameter :: DEBUG = .false.
    character(len=255) :: FileName      ! Duh
    integer :: I, J                     ! Loop inductor, subscript
    integer :: Lun                      ! Unit number for reading a file
    integer :: Son                      ! Some subtree of root.
    integer :: Version

    ! Error message codes

    ! We skip this stage if we're just a master task
    any_errors = 0                      ! At least clear error if we're master
    if ( parallel%master .and. .not. parallel%fwmParallel ) return

    error = 0
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_begin ( 'ForwardModelGlobalSetup', root )

    ! "Root" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.
    ! Collect data from the fields.

    do i = 2, nsons(root)
      Version = 1
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_antennaPatterns )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_antpats_start, &
            & 'Antenna Patterns File not found in PCF' )
          call open_antenna_patterns_file ( fileName, lun )
          call read_antenna_patterns_file ( lun )
          call close_antenna_patterns_file ( lun )
        end do
      case ( f_filterShapes )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_filtshps_start, &
            & 'Filter Shapes File not found in PCF' )
          call open_filter_shapes_file ( fileName, lun )
          call read_filter_shapes_file ( lun )
          call close_filter_shapes_file ( lun )
        end do
      case ( f_DACSfilterShapes )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_dacsfltsh_start, &
            & 'DACS Filter Shapes File not found in PCF' )
          call open_filter_shapes_file ( fileName, lun )
          call read_DACS_filter_shapes_file ( lun )
          call close_filter_shapes_file ( lun )
        end do
      case ( f_pointingGrids )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_ptggrids_start, &
            & 'Pointing Grids File not found in PCF' )
          call open_pointing_grid_file ( fileName, lun )
          call read_pointing_grid_file ( lun )
          call close_pointing_grid_file ( lun )
        end do
      case ( f_l2pc )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_l2pc_start, 'L2PC File not found in PCF', &
            & mlspcf_l2pc_end )
          if ( index ( fileName, '.txt' ) /= 0 ) then
            call open_l2pc_file ( fileName, lun)
            call read_l2pc_file ( lun )
            call close_l2pc_file ( lun )
          else
            call ReadCompleteHDF5L2PCFile ( fileName )
          end if
        end do
      case default
        ! Can't get here if the type checker worked
      end select
    end do

    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_end ( 'ForwardModelGlobalSetup' )
    any_errors = error

  contains

    ! ............................................  Get_File_Name  .....
    subroutine Get_File_Name ( pcfCode, MSG, pcfEndCode )
      integer, intent(in) :: pcfCode
      character(len=*), intent(in) :: MSG ! in case of error
      integer, intent(in), optional :: pcfEndCode

      character(len=255) :: PCFFileName, path
      integer :: returnStatus             ! non-zero means trouble
      integer :: mypcfEndCode

      call get_string ( sub_rosa(subtree(j,son)), fileName, strip=.true. )
      if ( TOOLKIT ) then
        mypcfEndCode = pcfCode
        if ( present(pcfEndCode) ) mypcfEndCode = pcfEndCode
        PCFFileName = fileName
        call split_path_name(PCFFileName, path, fileName)
        lun = GetPCFromRef(fileName, pcfCode, &
          & mypcfEndCode, &
          & TOOLKIT, returnStatus, Version, DEBUG, &
          & exactName=PCFFileName)
        if ( returnStatus /= 0 ) then
          call AnnounceError ( 0, son, extraMessage=MSG )
        else
          fileName = PCFFileName
        end if
      end if
    end subroutine Get_File_Name

  end subroutine ForwardModelGlobalSetup

  ! -----------------------------  CreateBinSelectorFromMLSCFINFO  -----
  type (BinSelector_T) function CreateBinSelectorFromMLSCFINFO ( root ) &
    & result ( binSelector )

    use Expr_M, only: EXPR
    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: L_NAMEFRAGMENT, L_VMR, L_TEMPERATURE, &
      & L_LATITUDE, L_SZA
    use Init_Tables_Module, only: F_COST, F_HEIGHT, F_MOLECULE, F_TYPE, &
      & F_NAMEFRAGMENT, F_EXACT
    use Intrinsic, only: T_NUMERIC_RANGE, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, &
      & PHYQ_INVALID, PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_VMR
    use L2PC_m, only: BINSELECTOR_T, BINSELECTORS, CREATEDEFAULTBINSELECTORS
    use MLSCommon, only: R8
    use MoreTree, only: Get_Field_ID, GET_BOOLEAN
    use Tree, only: Decoration, Nsons, Sub_Rosa, Subtree

    integer, intent(in) :: ROOT         ! Tree node
    ! Local variables
    integer :: SON                      ! Tree node
    integer :: GSON                     ! Tree node
    integer :: I                        ! Loop counters
    integer :: FIELD                    ! Field identifier
    logical :: GOT(field_first:field_last) ! "Got this field already"
    integer :: TYPE                     ! Type of value returned by expr
    integer :: EXPR_UNITS(2)            ! Units from expr
    real(r8) :: VALUE(2)                ! Value from expr
    integer :: COSTUNIT                 ! Units for cost
    integer :: WANTEDUNIT               ! Units wanted for cost

    ! Executable code

    ! Make sure the default ones are created first
    if ( .not. associated ( binSelectors ) ) &
      & call CreateDefaultBinSelectors

    ! Set up appropriate initial values
    binSelector%molecule = 0
    binSelector%nameFragment = 0
    binSelector%heightRange = 0.0
    binSelector%cost = 0.0
    binSelector%exact = .false.

    do i = 2, nsons(root)               ! Skip binSelector command
      son = subtree ( i, root )
      field = get_field_id ( son )
      if ( nsons(son) == 2 ) gson = subtree(2,son)
      got(field) = .true.
      select case ( field )
      case ( f_type )
        binSelector%selectorType = decoration(gson)
      case ( f_molecule )
        binSelector%molecule = decoration(gson)
      case ( f_nameFragment )
        binSelector%nameFragment = sub_rosa(gson)
      case ( f_exact )
        binSelector%exact = get_boolean(gson)
      case ( f_height )
        if ( nsons(son) > 2 ) call AnnounceError ( TooManyHeights, son, &
          & f_height )
        call expr ( gson, expr_units, value, type )
        if ( type /= t_numeric_range ) call AnnounceError ( 0, son, f_height, &
          & 'Height range expected' )
        if ( any ( expr_units /= phyq_pressure .and. expr_units /= phyq_dimensionless ) .or. &
          & all ( expr_units /= phyq_pressure ) ) &
          & call AnnounceError ( BadHeightUnit, son, f_height )
        binSelector%heightRange = value
      case ( f_cost )
        if ( nsons(son) > 2 ) call AnnounceError ( TooManyCosts, son, &
          & f_cost )
        call expr ( gson, expr_units, value, type )
        if ( type == t_numeric_range ) call AnnounceError ( 0, son, &
          & f_cost, 'Cost must not be a range' ) 
        ! Some units checking should probably go here in the long run !???? NJL
        binSelector%cost = value(1)
        costUnit = expr_units(1)
      end select
    end do

    wantedUnit = phyq_invalid
    select case ( binSelector%selectorType )
    case ( l_vmr )
      wantedUnit = phyq_vmr
    case ( l_temperature )
      wantedUnit = phyq_temperature
    case ( l_latitude, l_sza )
      wantedUnit = phyq_angle
    end select

    ! Just check a few last details
    if ( ( binSelector%selectorType == l_vmr ) &
      & .and. ( .not. got(f_molecule) ) ) call AnnounceError ( &
      & NoMolecule, son )
    if ( ( binSelector%selectorType == l_nameFragment ) .and. &
      & ( .not. got(f_nameFragment) ) ) call AnnounceError ( &
      0, son, extraMessage='No name fragment supplied' )
    if ( ( wantedUnit /= phyq_invalid ) .and. ( wantedUnit /= costUnit ) ) &
      & call AnnounceError ( 0, son, extraMessage='Wrong units for cost' )

  end function CreateBinSelectorFromMLSCFINFO

  ! --------------------------------  ConstructForwardModelConfig  -----
  type (forwardModelConfig_T) function ConstructForwardModelConfig &
    & ( NAME, ROOT, VGRIDS, GLOBAL ) result ( info )
    ! Process the forwardModel specification to produce ForwardModelConfig to add
    ! to the database

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Expr_M, only: EXPR
    use ForwardModelConfig, only: ForwardModelConfig_T, NullifyForwardModelConfig
    use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
    use Init_Tables_Module, only: L_FULL, L_SCAN, L_LINEAR, L_CLOUDFULL, L_HYBRID, &
      & L_POLARLINEAR
    use Init_Tables_Module, only: F_ALLLINESFORRADIOMETER, F_ATMOS_DER, &
      & F_BINSELECTORS, F_CHANNELS, F_CLOUD_DER, &
      & F_DEFAULT_spectroscopy, F_DIFFERENTIALSCAN, F_DO_BASELINE, F_DO_CONV, &
      & F_DO_FREQ_AVG, F_DO_1D, F_FREQUENCY, F_I_SATURATION, F_INCL_CLD, &
      & F_FORCESIDEBANDFRACTION, F_INTEGRATIONGRID, F_LOCKBINS, F_MODULE, F_MOLECULES, &
      & F_MOLECULEDERIVATIVES, F_NABTERMS, &
      & F_NAZIMUTHANGLES, F_NCLOUDSPECIES, F_NMODELSURFS, F_NSCATTERINGANGLES, &
      & F_NSIZEBINS, F_PHIWINDOW, F_POLARIZED, F_SIGNALS, F_SKIPOVERLAPS, &
      & F_SPECIFICQUANTITIES, F_SPECT_DER, F_SWITCHINGMIRROR, F_TANGENTGRID, F_TEMP_DER, &
      & F_TOLERANCE, F_TYPE, F_LINEARSIDEBAND
    use Intrinsic, only: L_NONE, L_CLEAR, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, &
      & PHYQ_PROFILES, PHYQ_TEMPERATURE
    use L2PC_m, only: BINSELECTORS, DEFAULTSELECTOR_LATITUDE, CREATEDEFAULTBINSELECTORS
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: HUNT
    use MLSSignals_M, only: Signals
    use MoreTree, only: Get_Boolean, Get_Field_ID
    use Parse_Signal_m, only: PARSE_SIGNAL
    use String_Table, only: Get_String
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use Tree, only: Decoration, Node_ID, Nsons, Sub_Rosa, Subtree
    use Tree_Types, only: N_Array, N_named
    use VGridsDatabase, only: VGrid_T

    integer, intent(in) :: NAME         ! The name of the config
    integer, intent(in) :: ROOT         ! of the forwardModel specification.
    !                                     Indexes either a "named" or
    !                                     "spec_args" vertex. Local variables
    type (vGrid_T), dimension(:), target :: vGrids ! vGrid database
    logical, intent(in) :: GLOBAL       ! Goes into info%globalConfig

    logical, dimension(:), pointer :: Channels   ! From Parse_Signal
    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: J                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: MoleculeSign             ! in the info%molecules array.
    integer :: Mol                      ! Tree indices of f_molecules, ...Pol
    integer :: NELTS                    ! Number of elements of an array tree
    integer :: SIDEBAND                 ! Returned from Parse_Signal
    integer, dimension(:), pointer :: SIGNALINDS ! From Parse_Signal
    character (len=80) :: SIGNALSTRING  ! E.g. R1A....
    integer :: Son                      ! Some subtree of root.
    integer :: STATUS                   ! From allocates etc.
    integer :: TANGENT                  ! Loop counter
    integer :: THISMOLECULE             ! Tree index.
    integer :: type                     ! Type of value returned by EXPR
    integer :: Expr_Units(2)            ! Units of value returned by EXPR
    real (r8) :: Value(2)               ! Value returned by EXPR
    integer :: WANTED                   ! Which signal do we want?

    ! Error message codes

    ! Nullify some pointers so allocate_test doesn't try to deallocate them.
    ! Don't initialize them with =>NULL() because that makes them SAVEd.

    nullify ( channels, signalInds )
    call NullifyForwardModelConfig ( info ) ! for Sun's rubbish compiler

    error = 0
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_begin ( "ConstructForwardModelConfig", root )

    ! Set sensible defaults
    info%allLinesForRadiometer = .false.
    info%atmos_der = .false.
    info%cloud_der = l_none
    info%DEFAULT_spectroscopy = .false.
    info%differentialScan = .false.
    info%do_1d = .false.
    info%do_baseline = .false.
    info%do_conv = .false.
    info%do_freq_avg = .false.
    info%forceSidebandFraction = .false.
    info%globalConfig = global
    info%incl_cld = .false.
    info%i_saturation = l_clear
    info%lockBins = .false.
    info%linearSideband = 0
    info%name = name
    info%no_cloud_species = 2
    info%no_model_surfs = 640
    info%num_ab_terms = 50
    info%num_azimuth_angles = 8
    info%num_scattering_angles = 16
    info%num_size_bins = 40
    info%phiwindow = 5
    info%polarized = .false.
    info%sideBandStart = -1
    info%sideBandStop = 1
    info%skipOverlaps = .false.
    info%spect_der = .false.
    info%switchingMirror= .false.
    info%temp_der = .false.
    info%windowUnits = phyq_profiles

    got = .false.
    info%tolerance = -1.0 ! Kelvins, in case the tolerance field is absent
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case ( f_allLinesForRadiometer )
        info%allLinesForRadiometer = get_boolean(son)
      case ( f_atmos_der )
        info%atmos_der = get_boolean(son)
      case  ( f_binSelectors )
        call Allocate_test ( info%binSelectors, nsons(son)-1, &
          & 'info%binSelectors', ModuleName )
        do j = 1, nsons(son) - 1
          info%binSelectors(j) = decoration ( decoration ( subtree ( j+1, son ) ) )
        end do
      case ( f_cloud_der )
        info%cloud_der = decoration ( subtree(2,son) )
      case ( f_DEFAULT_spectroscopy )
        info%DEFAULT_spectroscopy = get_boolean(son)
      case ( f_differentialScan )
        info%differentialScan = get_boolean(son)
      case ( f_do_baseline )
        info%do_baseline = get_boolean(son)
      case ( f_do_conv )
        info%do_conv = get_boolean(son)
      case ( f_do_freq_avg )
        info%do_freq_avg = get_boolean(son)
      case ( f_do_1d )
        info%do_1d = get_boolean(son)
      case ( f_forceSidebandFraction )
        info%forceSidebandFraction = get_boolean(son)
      case ( f_i_saturation )
        info%i_saturation = decoration(subtree(2,son))
      case ( f_incl_cld )
        info%incl_cld = get_boolean(son)
      case ( f_integrationGrid )
        info%integrationGrid => vGrids(decoration(decoration(subtree(2,son))))
      case ( f_linearSideband )
        call expr ( subtree(2,son), expr_units, value, type )
        info%linearSideband = nint ( value(1) )
        if ( expr_units(1) /= phyq_dimensionless ) &
          & call AnnounceError ( toleranceNotK, root )
      case ( f_lockBins )
        info%lockBins = get_boolean(son)
      case ( f_module )
        info%instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_molecules )
        mol = son
        nelts = 0
        do j = 2, nsons(son)
          call countElements ( subtree(j,son), nelts )
        end do
        call allocate_test ( info%molecules, nelts, "info%molecules", &
          & ModuleName )
        call allocate_test ( info%moleculeDerivatives, nelts, &
          & "info%moleculeDerivatives", ModuleName )
        info%moleculeDerivatives = .false.
        nelts = 0
        do j = 2, nsons(son)
          moleculeSign = +1 ! Indicate "root of a tree of molecules"
          call fillElements ( subtree(j,son), nelts, 0, info%molecules )
        end do
      case ( f_moleculeDerivatives )
        if ( .not. associated(info%molecules) ) then
          call announceError( DefineMoleculesFirst, root)
        else
          do j = 1, nsons(son)-1
            thisMolecule = decoration( subtree( j+1, son ) )
            if ( .not. any(info%molecules == thisMolecule ) ) &
              & call announceError( BadMolecule, root )
            where ( info%molecules == thisMolecule )
              info%moleculeDerivatives = .true.
            end where
          end do                          ! End loop over listed species
        end if
      case ( f_nabterms )
        call expr ( subtree(2,son), expr_units, value, type )
        info%NUM_AB_TERMS = nint( value(1) )
      case ( f_nazimuthangles )
        call expr ( subtree(2,son), expr_units, value, type )
      case ( f_ncloudspecies )
        call expr ( subtree(2,son), expr_units, value, type )
        info%no_cloud_species = nint( value(1) )
      case ( f_nmodelsurfs )
        call expr ( subtree(2,son), expr_units, value, type )
        info%no_model_surfs = nint( value(1) )
      case ( f_nscatteringangles )
        call expr ( subtree(2,son), expr_units, value, type )
        info%NUM_SCATTERING_ANGLES = nint( value(1) )
        info%NUM_AZIMUTH_ANGLES = nint( value(1) )
      case ( f_nsizebins )
        call expr ( subtree(2,son), expr_units, value, type )
        info%NUM_SIZE_BINS = nint( value(1) )
      case ( f_phiWindow )
        call expr ( subtree(2,son), expr_units, value, type )
        info%phiWindow = value(1)
        if ( all ( expr_units(1) /= (/ PHYQ_Profiles, PHYQ_Angle /) ) ) &
          call AnnounceError ( WrongUnitsForWindow, root )
        info%windowUnits = expr_units(1)
      case ( f_polarized )
        info%polarized = get_boolean(son)
      case ( f_signals )
        allocate ( info%signals (nsons(son)-1), stat = status )
        if ( status /= 0 ) call announceError( AllocateError, root )
        do j = 1, nsons(son)-1
          call get_string ( sub_rosa(subtree(j+1,son)), signalString, &
            & strip=.true.)
          call parse_Signal ( signalString, signalInds, &
            & tree_index=son, sideband=sideband, channels=channels )
          if ( .not. associated(signalInds) ) then ! A parse error occurred
            error = max(error,1)
            exit
          end if
          ! Later on choose the `right' one from the match
          ! For the moment choose first !????
          wanted=1
          info%signals(j) = signals(signalInds(wanted))
          info%signals(j)%sideband = sideband
          ! Don't hose channels in database, though shouldn't be an issue
          nullify ( info%signals(j)%channels ) 

          call allocate_Test ( info%signals(j)%channels, &
            & size(info%signals(j)%frequencies), 'info%signals%channels', &
            & ModuleName )
          if ( associated(channels) ) then
            info%signals(j)%channels(1:lbound(channels,1)-1) = .false.
            info%signals(j)%channels(lbound(channels,1):ubound(channels,1)) = &
              channels
            info%signals(j)%channels(ubound(channels,1)+1:) = .false.
          else
            info%signals(j)%channels = .true.
          end if
          call deallocate_test ( channels, 'channels', ModuleName )
          call deallocate_test ( signalInds, 'signalInds', ModuleName )
        end do                          ! End loop over listed signals
      case ( f_skipOverlaps )
        info%skipOverlaps = get_boolean(son)
      case  ( f_specificQuantities )
        call Allocate_test ( info%specificQuantities, nsons(son)-1, &
          & 'info%specificQuantities', ModuleName )
        do j = 1, nsons(son) - 1
          info%specificQuantities(j) = decoration ( decoration ( subtree ( j+1, son ) ) )
        end do
      case ( f_spect_der )
        info%spect_der = get_boolean(son)
      case ( f_switchingMirror )
        info%switchingMirror = get_boolean(son)
      case ( f_tangentGrid )
        info%tangentGrid => vGrids(decoration(decoration(subtree(2,son))))
      case ( f_temp_der )
        info%temp_der = get_boolean(son)
      case ( f_tolerance )
        call expr ( subtree(2,son), expr_units, value, type )
        info%tolerance = value(1)
        if ( expr_units(1) /= phyq_temperature ) &
          & call AnnounceError ( toleranceNotK, root )
      case ( f_type )
        info%fwmType = decoration(subtree(2,son))
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do ! i = 2, nsons(root)

    ! Now some more error checking
    select case ( info%fwmType )
    case ( l_full, l_hybrid )
      if ( .not. all(got( (/ f_molecules, f_signals, f_integrationGrid, &
        & f_tangentGrid /) )) ) call AnnounceError ( IncompleteFullFwm, root )

      ! Now identify the Earth's surface in the tangent grid
      call Hunt ( info%tangentGrid%surfs(:,1), info%integrationGrid%surfs(1,1), &
        &  info%surfaceTangentIndex )

      ! Ensure that points in tangentGrid at and above the surface are a subset
      ! of integration grid
      do tangent = info%surfaceTangentIndex, info%tangentGrid%noSurfs
        if ( .not. any ( abs( info%tangentGrid%surfs(tangent,1) - &
          & info%integrationGrid%surfs(:,1) ) < 1e-4 ) ) &
          & call AnnounceError ( TangentNotSubset, root )
      end do

      ! Make sure signal specifications make sense; get sideband Start/Stop
      call validateSignals
    case ( l_cloudfull )

      ! full cloud forward model

    case ( l_scan )
      ! Add 1d/2d method later probably !??? NJL
      if ( any(got( (/f_atmos_der,f_channels,f_do_conv,f_do_baseline,f_do_freq_avg,&
        & f_do_1d, f_incl_cld, f_frequency, f_molecules, f_moleculeDerivatives, &
        & f_signals, f_spect_der, f_temp_der /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root )

    case ( l_linear, l_polarLinear )
      if ( .not. all(got( (/f_signals/) )) ) & ! Maybe others later
        & call AnnounceError ( IncompleteLinearFwm, root )
      if ( any(got( (/f_do_conv,f_do_freq_avg,f_do_1d,f_incl_cld,f_frequency /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root )
      ! Make sure we get a default bin selector
      if ( .not. associated ( info%binSelectors ) ) then
        if ( .not. associated ( binSelectors ) ) call CreateDefaultBinSelectors
        call Allocate_test ( info%binSelectors, 1, 'info%binSelectors', ModuleName )
        info%binSelectors = DefaultSelector_Latitude
      end if

    end select

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An error occured' )
    if ( toggle(gen) .and. levels(gen) > 0 ) &
      & call trace_end ( "ConstructForwardModelConfig" )

  contains
    ! ............................................  CountElements  .....
    recursive subroutine CountElements ( root, count )
      integer, intent(in) :: ROOT       ! of a subtree
      integer, intent(inout) :: COUNT   ! Number of array elements
      integer :: I                      ! Subtree index, loop inductor
      if ( node_id(root) == n_array ) then
        do i = 1, nsons(root)
          call countElements ( subtree(i,root), count )
        end do
      else
        count = count + 1
      end if
    end subroutine CountElements

    ! .............................................  FillElements  .....
    recursive subroutine FillElements ( root, count, depth, molecules )
      integer, intent(in) :: ROOT       ! of a subtree
      integer, intent(inout) :: COUNT   ! of array elements processed
      integer, intent(in) :: DEPTH      ! in the array tree
      integer, intent(inout) :: MOLECULES(:)      ! The array to be filled
      integer :: I                      ! Subtree index, loop inductor
      if ( node_id(root) == n_array ) then
        do i = 1, nsons(root)
          call fillElements ( subtree(i,root), count, depth+1, molecules )
        end do
      else
        count = count + 1
        molecules(count) = moleculeSign * decoration(root)
        moleculeSign = -1 ! Indicate "Part of a tree of molecules"
      end if
    end subroutine FillElements

    ! ..........................................  ValidateSignals  .....
    subroutine ValidateSignals
      use MLSSignals_m, only: Signal_t
      use MLSNumerics, only: EssentiallyEqual
      type (Signal_T), pointer :: FirstSignal

      firstSignal => info%signals(1)

      ! Make sure all the signals we're dealing with are same module,
      ! radiometer and sideband.
      if ( any( info%signals%sideband /= firstSignal%sideband ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  "Can't have mixed sidebands in forward model config" )
      if ( .not. all ( EssentiallyEqual ( info%signals%lo, firstSignal%lo ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  "Can't have mixed radiometers in forward model config" )

      ! Think about sidebands
      if ( ( firstSignal%sideband == 0 ) .and.&
        &  ( firstSignal%singleSideband == 0 ) ) then
        ! Do a folded measurement
        info%sidebandStart = -1
        info%sidebandStop = 1
      else
        ! It's either a single sideband radiometer, or the user requested a
        ! specific sideband.
        ! Check sanity, if they are both non zero they should be the same.
        if ( ( firstSignal%singleSideband /= 0 ) .and. &
          &  ( firstSignal%sideband /= 0 ) .and. &
          &  ( firstSignal%singleSideband /= &
          &    firstSignal%sideband ) ) call MLSMessage ( &
          &      MLSMSG_Error, ModuleName, &
          &      "User requested a sideband that doesn't exist" )
        ! OK, use whichever one is given
        if ( firstSignal%singleSideband /= 0 ) then
          info%sidebandStart = firstSignal%singleSideband
        else
          info%sidebandStart = firstSignal%sideband
        end if
        info%sidebandStop = info%sidebandStart
      end if
    end subroutine ValidateSignals

  end function ConstructForwardModelConfig

  ! ------------------------------------  PrintForwardModelTiming  -----
  subroutine PrintForwardModelTiming ( FWModelConfig )
  !  Print mean, std_dev for timing FullforwardModel
                                                                                
    use ForwardModelConfig, only: ForwardModelConfig_T
    use Output_m, only: BLANKS, Output
    use String_Table, only: DISPLAY_STRING, GET_STRING
                                                                                
    ! Dummy argument
    type(ForwardModelConfig_T), intent(inout) :: FWModelConfig
                                                                                
    ! Local variables
    real :: mean, mean_sqDelta, meanTimes, std_dev
    character(len=25) :: thisName

    mean = FWModelConfig%sum_DeltaTime/FWModelConfig%Ntimes
    mean_sqDelta =  FWModelConfig%sum_squareDeltaTime / FWModelConfig%Ntimes 
    if (FWModelConfig%Ntimes <= 1) then
       meanTimes = 1.0 
    else
       meanTimes = FWModelConfig%Ntimes / (FWModelConfig%Ntimes - 1)
    end if
    std_dev = sqrt(abs(meanTimes * (mean_sqDelta - (mean*mean)))) 
    call get_string ( FWModelConfig%name, thisName )
    call output ( thisName, advance='no')
    call output ( FWModelConfig%Ntimes, format = '(i6)',advance = 'no' )
    call blanks (12, advance = 'no')
    if (FWModelconfig%Ntimes == 0) then
	mean = 0.0
	std_dev = 0.0
    end if
    call output ( mean, format='(f8.2)', advance = 'no' )
    call blanks (10, advance = 'no')
    call output ( std_dev, format='(f8.2)', advance = 'yes' )
                                                                                
  end subroutine PrintForwardModelTiming
                                                                                
  ! ------------------------------------  ResetForwardModelTiming  -----
  subroutine ResetForwardModelTiming ( FWModelConfig )
  ! reset the values for timing FullforwardModel
                                                                                
    use ForwardModelConfig, only: ForwardModelConfig_T
                                                                                
    ! Dummy argument
    type(ForwardModelConfig_T), intent(inout) :: FWModelConfig
                                                                                
    FWModelConfig%Ntimes = 0
    FWModelConfig%sum_DeltaTime = 0.0
    FWModelConfig%sum_squareDeltaTime = 0.0
                                                                                
  end subroutine ResetForwardModelTiming

  ! =====     Private Procedures     ===================================
  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( Code, where, FieldIndex, extraMessage )

    use Intrinsic, only: Lit_Indices
    use Lexer_Core, only: PRINT_SOURCE
    use Output_M, only: Output
    use String_Table, only: Display_String
    use Tree, only: Source_Ref

    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: where      ! Where in the tree did the error occur?
    integer, intent(in), optional :: FieldIndex ! f_...
    character (LEN=*), optional :: extraMessage

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref ( where ) )
    call output ( ' ForwardModelSupport complained: ' )
    select case ( code )
    case ( AllocateError )
      call output ( 'allocation error.', advance='yes' )
    case ( BadMolecule )
      call output ( 'asked for derivatives for an unlisted molecule.', &
        & advance='yes' )
    case ( BadBinSelectors )
      call output ('cannot have a fieldAzimuth binSelector for polarlinear model', &
        & advance='yes' )
    case ( DefineMoleculesFirst )
      call output ( 'molecule must be defined before moleules derivatives.', &
        & advance='yes')
    case ( DefineSignalsFirst )
      call output ( 'signals must be defined before channels.',advance='yes')
    case ( IncompleteBinSelectors )
      call output ('must have some binSelectors for the polarlinear model',advance='yes' )
    case ( IncompleteFullFwm )
      call output ('incomplete full foward model specification',advance='yes' )
    case ( IncompleteLinearFwm )
      call output ( 'incomplete linear foward model specification', &
        & advance='yes' )
    case ( IrrelevantFwmParameter )
      call output ( 'irrelevant parameter for this forward model type', &
        & advance='yes' )
    case ( LinearSidebandHasUnits )
      call output ( 'irrelevant units for this linear sideband', &
        & advance='yes' )
    case ( TangentNotSubset )
      call output ('non subsurface tangent grid not a subset of integration&
        & grid', advance='yes' )
    case ( ToleranceNotK )
      call output ( 'tolerance does not have dimensions of temperature/radiance',&
        & advance='yes' )
    case ( TooManyHeights )
      call output ( 'Bin Selectors can only refer to one height range',&
        & advance='yes' )
    case ( BadHeightUnit )
      call output ( 'Inappropriate units for height in binSelector',&
        & advance='yes' )
    case ( TooManyCosts )
      call output ( 'Bin Selectors can only have one cost',&
        & advance='yes' )
    case ( BadQuantityType )
      call output ( 'Bin Selectors cannot apply to this quantity type',&
        & advance='yes' )
    case ( NoMolecule )
      call output ( 'A bin selector of type vmr must have a molecule',&
        & advance='yes' )
    case ( WrongUnitsForWindow )
      call output ( 'phiWindow must be in degrees or profiles', &
        & advance='yes' )
   case default
      call output ( '(no specific description of this error)', advance='yes' )
    end select
    if ( present(extraMessage) ) call output ( extraMessage, advance='yes' )
  end subroutine AnnounceError

  logical function NOT_USED_HERE()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function NOT_USED_HERE

end module ForwardModelSupport

! $Log$
! Revision 2.83  2003/08/27 20:31:07  livesey
! Removed the prevention of tolerance>0.0 on polarized runs
!
! Revision 2.82  2003/08/21 21:15:18  cvuu
! Change the output format for fullForwardModel Timing
!
! Revision 2.81  2003/08/19 05:51:31  livesey
! Extra call to CreateDefaultBinSelectors
!
! Revision 2.80  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.79  2003/08/15 20:28:44  vsnyder
! Remove check for and prohibition against polarized VMR derivatives
!
! Revision 2.78  2003/08/14 20:25:06  livesey
! Added the exact bin selector stuff
!
! Revision 2.77  2003/08/13 00:49:40  livesey
! Added the polarLinear model and now ensures that a default bin selector
! is provided to the linear model if the user doesn't supply one.
!
! Revision 2.76  2003/08/12 18:16:13  livesey
! Forbid non negative tolerance in polarized model.
!
! Revision 2.75  2003/08/12 17:12:08  livesey
! Added check to ensure you don't ask for mixing ratio derivatives from
! the polarized model
!
! Revision 2.74  2003/08/11 22:35:34  livesey
! Now more forgiving of running e.g. R5H and R5V together.
!
! Revision 2.73  2003/07/21 18:13:44  pwagner
! Looks farther than just 1st l2pc file in pcf
!
! Revision 2.72  2003/07/18 22:54:35  pwagner
! Ignores path in file names if TOOLKIT
!
! Revision 2.71  2003/07/16 21:51:29  pwagner
! Uses mlspcf_dacsfltsh_start
!
! Revision 2.70  2003/07/16 01:27:35  vsnyder
! Futzing
!
! Revision 2.69  2003/07/16 01:06:36  vsnyder
! Add DACS filter shapes
!
! Revision 2.68  2003/07/15 22:10:59  livesey
! Added support for hybrid model
!
! Revision 2.67  2003/07/15 18:17:50  livesey
! Better handling of configuration name
!
! Revision 2.66  2003/06/30 22:55:01  cvuu
! Find mean, std dev timing of fullForwardModel calls
!
! Revision 2.65  2003/06/26 23:15:07  vsnyder
! Reinstate revision 2.63, which inscrutably got lost
!
! Revision 2.64  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.63  2003/06/18 01:57:12  vsnyder
! Move checking that all signals in a config are for the same radiometer,
! module and sideband here from FullForwardModel.  Move computation for
! SidebandStart and SidebandStop here from FullForwardModel.
!
! Revision 2.62  2003/06/09 22:50:13  pwagner
! Reduced everything (PCF, PUNISH.., etc.) to TOOLKIT
!
! Revision 2.61  2003/05/29 16:42:19  livesey
! New switchingMirror stuff
!
! Revision 2.60  2003/05/05 23:00:34  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.56.2.2  2003/04/08 23:40:50  jonathan
! remove cloud_fov
!
! Revision 2.56.2.1  2003/02/22 00:47:26  vsnyder
! Delete moleculesPol, moleculeDerivativesPol, add Polarized to ForwardModelConfig
!
! Revision 2.56  2003/02/08 05:28:47  vsnyder
! Squash a bug (looking at type before evaluating the expression).
! Move USE statements from module scope to procedure scope, so as not to
! cause TREE_WALKER to be recompiled every time there's a minor change here.
!
! Revision 2.55  2003/02/06 22:04:48  vsnyder
! Add f_moleculesPol, f_moleculeDerivativesPol, delete f_polarized
!
! Revision 2.54  2003/02/06 00:45:51  livesey
! Added new binSelector stuff
!
! Revision 2.53  2003/02/06 00:15:01  jonathan
! sorry no_cloud_sps=2
!
! Revision 2.52  2003/02/05 23:27:32  jonathan
! change default no_cloud_species=1
!
! Revision 2.51  2003/02/05 21:56:39  livesey
! New binSelectors stuff
!
! Revision 2.50  2003/02/04 22:02:18  jonathan
! set i_saturation==0 as default
!
! Revision 2.49  2003/02/04 19:03:06  livesey
! Default tolerance now -1.0
!
! Revision 2.48  2003/01/30 17:28:21  jonathan
! add logical incl_cld
!
! Revision 2.47  2003/01/29 01:48:29  vsnyder
! Add 'polarized' field to forwardModel
!
! Revision 2.46  2003/01/27 16:51:08  livesey
! Added initialisation for windowUnits
!
! Revision 2.45  2003/01/26 04:42:55  livesey
! Added units for phiWindow
!
! Revision 2.44  2003/01/16 00:55:41  jonathan
! add do_1d
!
! Revision 2.43  2003/01/13 17:17:04  jonathan
!  change cloud_width to i_saturation
!
! Revision 2.42  2003/01/03 21:03:02  pwagner
! l2pc filenames now inputtable via PCF
!
! Revision 2.41  2002/11/22 12:20:13  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.40  2002/11/15 01:33:24  livesey
! Added allLinesForRadiometer stuff
!
! Revision 2.39  2002/10/18 22:44:11  vsnyder
! Remove some unreferenced USE names
!
! Revision 2.38  2002/10/18 18:01:58  livesey
! Ensure l2pc files etc. read in fwmParallel master mode.
!
! Revision 2.37  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.36  2002/09/25 20:08:26  livesey
! Added globalConfig and specificQuantities
!
! Revision 2.35  2002/08/21 23:31:52  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.34  2002/08/04 16:02:23  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.33  2002/07/17 06:02:36  livesey
! New HDF5 l2pc stuff
!
! Revision 2.32  2002/06/12 17:01:54  livesey
! Stuff to support change to real phiWindow from integer
!
! Revision 2.31  2002/05/14 22:31:31  livesey
! Minor bug fix associated with running in parallel mode.
!
! Revision 2.30  2002/03/21 16:42:34  livesey
! Made it skip reading l2pc files etc for parallel master tasks.
!
! Revision 2.29  2002/03/15 21:22:31  livesey
! Dealt with new BinSelector type
!
! Revision 2.28  2002/03/14 23:30:21  pwagner
! Changed id of who announcesError
!
! Revision 2.27  2002/03/12 23:44:12  vsnyder
! Removed private attribute from CVS stuff because private is module's default
!
! Revision 2.26  2002/03/07 17:17:57  livesey
! Removed frqGap
!
! Revision 2.25  2002/02/14 23:02:43  livesey
! Added a safety net to guard against deallocating channels in signals
! database (shouldn't occur in current version, though may be relevant
! later).
!
! Revision 2.24  2002/02/13 00:08:50  livesey
! Added differential scan model
!
! Revision 2.23  2002/02/08 22:52:21  livesey
! Hooked up bin selectors
!
! Revision 2.22  2002/02/04 23:24:49  livesey
! Removed dumps
!
! Revision 2.21  2002/01/21 21:13:28  livesey
! Added binSelector parsing
!
! Revision 2.20  2001/12/17 18:26:37  vsnyder
! Improve method to put '-' sign on 'part of a tree of molecules'
!
! Revision 2.19  2001/11/29 00:27:56  vsnyder
! Fix blunders in arrays-of-arrays, alphabetize USEs
!
! Revision 2.18  2001/11/28 03:50:07  vsnyder
! Allow array of arrays for 'molecules' field
!
! Revision 2.17  2001/11/15 23:49:40  jonathan
! rename DF_spectroscopy to default_spectroscopy
!
! Revision 2.16  2001/11/15 20:55:55  jonathan
! add df_spectroscopy
!
! Revision 2.15  2001/11/08 00:12:52  livesey
! Removed an obsolete variable
!
! Revision 2.14  2001/10/31 15:26:59  livesey
! Allowed filename arguments to forward model global stuff to be arrays
!
! Revision 2.13  2001/10/02 20:34:50  livesey
! Added do_baseline stuff
!
! Revision 2.12  2001/09/04 15:58:02  jonathan
! add cloud_fov, jonathan
!
! Revision 2.11  2001/07/17 22:36:19  jonathan
! add cloud_width, jonathan/paul
!
! Revision 2.10  2001/07/16 22:07:21  jonathan
! change cloud_der to integer-type, jonathan
!
! Revision 2.9  2001/07/12 23:27:48  livesey
! Got rid of the s_cloudForwardModel stuff
!
! Revision 2.8  2001/07/09 22:51:37  pwagner
! Picks up cloud-related config components
!
! Revision 2.7  2001/06/21 20:06:42  vsnyder
! Make the tolerance field of the ForwardModel spec optional
!
! Revision 2.6  2001/06/21 15:05:20  livesey
! Added tolerance field to config
!
! Revision 2.5  2001/06/19 22:48:42  pwagner
! Eliminated duplicate declaration of PointingGrid_m
!
! Revision 2.4  2001/06/04 22:42:24  livesey
! Added belowRef to forwardModelIntermediate_T
!
! Revision 2.3  2001/05/31 22:08:06  livesey
! Added cloud_der flag
!
! Revision 2.2  2001/05/30 23:05:39  pwagner
! Uses PCF for 3 fwdmdl files
!
! Revision 2.1  2001/05/29 23:18:18  livesey
! First version, was ForwardModelInterface
!

