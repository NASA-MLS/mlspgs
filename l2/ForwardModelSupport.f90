! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module ForwardModelSupport

  ! Set up the forward model stuff.

  implicit none
  private
  public :: CONSTRUCTFORWARDMODELCONFIG, FORWARDMODELGLOBALSETUP, &
    & CREATEBINSELECTORFROMMLSCFINFO, PRINTFORWARDMODELTIMING, &
    & RESETFORWARDMODELTIMING, SHOWFWDMODELNAMES, FILLFWDMODELTIMINGS

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Error codes

  integer, parameter :: AllocateError          = 1
  integer, parameter :: BadBinSelectors        = AllocateError + 1
  integer, parameter :: BadHeightUnit          = BadBinSelectors + 1
  integer, parameter :: BadMoleculeGroup       = BadHeightUnit + 1
  integer, parameter :: BadQuantityType        = BadMoleculeGroup + 1
  integer, parameter :: CloudHas               = BadQuantityType + 1
  integer, parameter :: CloudLBL               = CloudHas + 1
  integer, parameter :: CloudNeeds             = CloudLBL + 1
  integer, parameter :: CloudNot               = CloudNeeds + 1
  integer, parameter :: DerivSansMolecules     = CloudNot  + 1
  integer, parameter :: DuplicateMolecule      = DerivSansMolecules + 1
  integer, parameter :: FirstSansFirst1        = DuplicateMolecule + 1
  integer, parameter :: FirstSansFirst2        = FirstSansFirst1 + 1
  integer, parameter :: Hess_notJac            = FirstSansFirst2 + 1
  integer, parameter :: IncompleteBinSelectors = Hess_notJac + 1
  integer, parameter :: IncompleteFullFwm      = IncompleteBinSelectors + 1
  integer, parameter :: IncompleteLinearFwm    = IncompleteFullFwm + 1
  integer, parameter :: IrrelevantFwmParameter = IncompleteLinearFwm + 1
  integer, parameter :: LBLandPFA              = IrrelevantFwmParameter + 1
  integer, parameter :: LinearSidebandHasUnits = LBLandPFA + 1
  integer, parameter :: LineNotMolecule        = LinearSidebandHasUnits + 1
  integer, parameter :: LineParamTwice         = LineNotMolecule + 1
  integer, parameter :: MIFTransformation_signals = LineParamTwice + 1
  integer, parameter :: NeedBothXYStar         = MIFTransformation_signals + 1
  integer, parameter :: Nested                 = NeedBothXYStar + 1
  integer, parameter :: NoArray                = Nested + 1
  integer, parameter :: NoBetaGroup            = NoArray + 1
  integer, parameter :: NoMolecule             = NoBetaGroup + 1
  integer, parameter :: NoPolarizedAndPFA      = NoMolecule + 1
  integer, parameter :: PFANotMolecule         = NoPolarizedAndPFA + 1
  integer, parameter :: PFATwice               = PFANotMolecule + 1
  integer, parameter :: PolarizedAndAllLines   = PFATwice + 1
  integer, parameter :: SecondSansFirst        = PolarizedAndAllLines  + 1
  integer, parameter :: SecondSansSecond1      = SecondSansFirst + 1
  integer, parameter :: SecondSansSecond2      = SecondSansSecond1 + 1
  integer, parameter :: TangentNotSubset       = SecondSansSecond2 + 1
  integer, parameter :: TooManyCosts           = TangentNotSubset + 1
  integer, parameter :: TooManyHeights         = TooManyCosts + 1
  integer, parameter :: WrongUnitsForWindow    = TooManyHeights + 1

  integer :: Error            ! Error level -- 0 = OK

contains ! =====     Public Procedures     =============================

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root, any_errors, fileDataBase )
    ! Process the forwardModel specification to produce ForwardModelInfo.

    use ANTENNAPATTERNS_M, only: OPEN_ANTENNA_PATTERNS_FILE, &
      & READ_ANTENNA_PATTERNS_FILE, CLOSE_ANTENNA_PATTERNS_FILE
    use FILTERSHAPES_M, only: OPEN_FILTER_SHAPES_FILE, &
      & READ_FILTER_SHAPES_FILE, READ_DACS_FILTER_SHAPES_FILE, &
      & CLOSE_FILTER_SHAPES_FILE
    use INIT_TABLES_MODULE, only: F_ANTENNAPATTERNS, F_DACSFILTERSHAPES, &
      & F_FILTERSHAPES, F_L2PC, F_MIETABLES, F_PFAFILES, F_POINTINGGRIDS
    use INTRINSIC, only: L_ASCII, L_HDF
    use L2PARINFO, only: PARALLEL
    use L2PC_M, only: READCOMPLETEHDF5L2PCFILE
    use MLSCOMMON, only: MLSFILE_T
    use MLSPCF2, only: MLSPCF_ANTPATS_START, MLSPCF_FILTSHPS_START, &
      &          MLSPCF_DACSFLTSH_START, MLSPCF_PTGGRIDS_START, &
      &          MLSPCF_L2PC_START, MLSPCF_L2PC_END, &
      &          MLSPCF_MIETABLES_START, &
      &          MLSPCF_PFA_START, MLSPCF_PFA_END
    use MORETREE, only: GET_FIELD_ID
    use PFADATABASE_M, only: PROCESS_PFA_FILE
    use POINTINGGRID_M, only: CLOSE_POINTING_GRID_FILE, &
      & OPEN_POINTING_GRID_FILE, READ_POINTING_GRID_FILE
    use READ_MIE_M, only: READ_MIE
    use TOGGLES, only: GEN, LEVELS, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: NSONS, SUB_ROSA, SUBTREE

    integer, intent(in) :: Root         ! of the forwardModel specification.
    !                                     Indexes a "spec_args" vertex.
    integer, intent(out) :: any_errors  ! non-zero means trouble
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE

    ! Internal variables
    logical, parameter :: DEBUG = .false.
    integer :: FileIndex                ! In the string table
    character(len=255) :: FileName      ! Duh
    integer :: I, J                     ! Loop inductor, subscript
    integer :: Lun                      ! Unit number for reading a file
    integer :: Me = -1                  ! String index for trace
    type (MLSFile_T), pointer   :: MLSFile
    integer :: Son                      ! Some subtree of root.
    integer :: Version
    integer, save :: last_l2pc = mlspcf_l2pc_start - 1
    integer, save :: last_pfa  = mlspcf_pfa_start  - 1

    ! Error message codes

    ! We skip this stage if we're just a master task
    any_errors = 0                      ! At least clear error if we're master
    if ( parallel%master .and. .not. parallel%fwmParallel ) return

    error = 0
    call trace_begin ( me, 'ForwardModelGlobalSetup', root, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

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
            & get_field_id(son), filedatabase, MLSFile, &
            & 'Antenna Patterns File not found in PCF' )
          call open_antenna_patterns_file ( fileName, lun )
          call read_antenna_patterns_file ( lun, subtree(j,son) )
          call close_antenna_patterns_file ( lun )
        end do
      case ( f_DACSfilterShapes )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_dacsfltsh_start, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'DACS Filter Shapes File not found in PCF' )
          call open_filter_shapes_file ( fileName, lun, fileIndex )
          call read_DACS_filter_shapes_file ( lun, fileIndex, subtree(j,son) )
          call close_filter_shapes_file ( lun )
        end do
      case ( f_filterShapes )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_filtshps_start, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'Filter Shapes File not found in PCF' )
          call open_filter_shapes_file ( fileName, lun, fileIndex )
          call read_filter_shapes_file ( lun, fileIndex, subtree(j,son) )
          call close_filter_shapes_file ( lun )
        end do
      case ( f_l2pc )
        last_l2pc = last_l2pc + 1
        do j = 2, nsons(son)
          call get_file_name ( last_l2pc, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'L2PC File not found in PCF', mlspcf_l2pc_end )
          call ReadCompleteHDF5L2PCFile ( MLSFile, subtree(j,son) )
        end do
      case ( f_MieTables )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_MieTables_start, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'Mie tables File not found in PCF' )
          call read_mie ( fileName )
        end do
      case ( f_PFAFiles )
        last_pfa = last_pfa + 1
        do j = 2, nsons(son)
          call get_file_name ( last_pfa, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'PFA File not found in PCF', mlspcf_pfa_end )
          if ( index ( fileName, '.h5' ) /= 0 ) then
            if ( process_PFA_File ( filename, &
            & subtree(j,son) ) /= 0 ) continue
          endif
        end do
      case ( f_pointingGrids )
        do j = 2, nsons(son)
          call get_file_name ( mlspcf_ptggrids_start, &
            & get_field_id(son), filedatabase, MLSFile, &
            & 'Pointing Grids File not found in PCF' )
          call open_pointing_grid_file ( fileName, lun )
          call read_pointing_grid_file ( lun, subtree(j,son) )
          call close_pointing_grid_file ( lun )
        end do
      case default
        ! Can't get here if the type checker worked
      end select
    end do

    call trace_end ( 'ForwardModelGlobalSetup', &
      & cond=toggle(gen) .and. levels(gen) > 0 )
    any_errors = error

  contains

    ! ............................................  Get_File_Name  .....
    subroutine Get_File_Name ( pcfCode, &
      & fileType, fileDataBase, MLSFile, MSG, pcfEndCode )
      use HDF, only: DFACC_RDONLY
      use INIT_TABLES_MODULE, only: FIELD_INDICES
      use MLSCOMMON, only: MLSFILE_T
      use MLSFILES, only: HDFVERSION_5, &
        & ADDINITIALIZEMLSFILE, GETPCFROMREF, SPLIT_PATH_NAME
      use MLSL2OPTIONS, only: TOOLKIT
      use SDPTOOLKIT, only: PGS_PC_GETREFERENCE
      use STRING_TABLE, only: GET_STRING
      ! Dummy args
      integer, intent(in) :: pcfCode
      integer, intent(in) :: fileType ! f_l2pc, f_antennaPatterns, etc.
      type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
      type (MLSFile_T), pointer   :: MLSFile
      character(len=*), intent(in) :: MSG ! in case of error
      integer, intent(in), optional :: pcfEndCode
      ! Internal variables
      character(len=255) :: fileTypeStr, PCFFileName, path, shortName
      integer :: returnStatus             ! non-zero means trouble
      integer :: mypcfEndCode
      ! Executable
      mypcfEndCode = 0
      lun = 0
      call get_string ( sub_rosa(subtree(j,son)), shortName, strip=.true. )
      fileName = shortName
      call get_string ( field_indices(fileType), fileTypeStr, strip=.true. )
      if ( TOOLKIT ) then
        mypcfEndCode = pcfCode
        if ( present(pcfEndCode) ) mypcfEndCode = pcfEndCode
        if ( fileName == ' ' ) then
          returnStatus = Pgs_pc_getReference(pcfCode, version, &
            & fileName)
          lun = pcfCode
        else
          PCFFileName = fileName
          call split_path_name ( PCFFileName, path, fileName )
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
      end if
      if ( fileType == f_l2pc ) then
        MLSFile => AddInitializeMLSFile(filedatabase, &
          & content=fileTypeStr, &
          & name=Filename, shortName=shortName, &
          & type=l_hdf, access=dfacc_rdonly, HDFVersion=HDFVERSION_5)
      else
        MLSFile => AddInitializeMLSFile(filedatabase, &
          & content=fileTypeStr, &
          & name=Filename, shortName=shortName, &
          & type=l_ascii, access=dfacc_rdonly)
      endif      
      MLSFile%PCFId = lun
    end subroutine Get_File_Name

  end subroutine ForwardModelGlobalSetup

  ! -----------------------------  CreateBinSelectorFromMLSCFINFO  -----
  type (BinSelector_T) function CreateBinSelectorFromMLSCFINFO ( root ) &
    & result ( binSelector )

    use EXPR_M, only: EXPR
    use INIT_TABLES_MODULE, only: FIELD_FIRST, FIELD_LAST
    use INIT_TABLES_MODULE, only: L_NAMEFRAGMENT, L_VMR, L_TEMPERATURE, &
      & L_LATITUDE, L_SZA
    use INIT_TABLES_MODULE, only: F_COST, F_HEIGHT, F_MOLECULE, F_TYPE, &
      & F_NAMEFRAGMENT, F_EXACT
    use INTRINSIC, only: T_NUMERIC_RANGE, PHYQ_ANGLE, PHYQ_DIMENSIONLESS, &
      & PHYQ_INVALID, PHYQ_PRESSURE, PHYQ_TEMPERATURE, PHYQ_VMR
    use L2PC_M, only: BINSELECTOR_T, BINSELECTORS, CREATEDEFAULTBINSELECTORS
    use MLSKINDS, only: R8
    use MORETREE, only: GET_FIELD_ID, GET_BOOLEAN
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE

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
        if ( nsons(son) > 2 ) call AnnounceError ( TooManyHeights, son )
        call expr ( gson, expr_units, value, type )
        if ( type /= t_numeric_range ) call AnnounceError ( 0, son, &
          & extraMessage='Height range expected' )
        if ( any ( expr_units /= phyq_pressure .and. expr_units /= phyq_dimensionless ) .or. &
          & all ( expr_units /= phyq_pressure ) ) &
          & call AnnounceError ( BadHeightUnit, son )
        binSelector%heightRange = value
      case ( f_cost )
        if ( nsons(son) > 2 ) call AnnounceError ( TooManyCosts, son )
        call expr ( gson, expr_units, value, type )
        if ( type == t_numeric_range ) call AnnounceError ( 0, son, &
          & extraMessage='Cost must not be a range' )
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
    & ( name, root, global ) result ( info )
    ! Process the forwardModel specification to produce ForwardModelConfig to
    ! add to the database

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use EXPR_M, only: EXPR
    use FORWARDMODELCONFIG, only: DUMP, FORWARDMODELCONFIG_T, &
      & LINECENTER, LINEWIDTH, LINEWIDTH_TDEP, &
      & NULLIFYFORWARDMODELCONFIG, SPECTROPARAM_T
    use INIT_TABLES_MODULE, only: FIELD_FIRST, FIELD_LAST
    use INIT_TABLES_MODULE, only: L_FULL, L_SCAN, L_LINEAR, L_CLOUDFULL, L_HYBRID, &
      & L_POLARLINEAR
    use INIT_TABLES_MODULE, only:  F_ALLLINESFORRADIOMETER, &
      & F_ALLLINESINCATALOG, F_ATMOS_DER, F_ATMOS_SECOND_DER, &
      & F_BINSELECTORS, F_CHANNELS, F_CLOUD_DER, F_DEFAULT_SPECTROSCOPY, &
      & F_DIFFERENTIALSCAN, F_DO_1D, F_DO_BASELINE, F_DO_CONV, &
      & F_DO_FREQ_AVG, F_FORCESIDEBANDFRACTION, F_FREQUENCY, F_FRQTOL, &
      & F_IGNOREHESSIAN, F_INCL_CLD, F_INTEGRATIONGRID, F_I_SATURATION, &
      & F_LINEARSIDEBAND, F_LINECENTER, F_LINEWIDTH, F_LINEWIDTH_TDEP, &
      & F_LOCKBINS, F_LSBLBLMOLECULES, F_LSBPFAMOLECULES, F_MODULE, &
      & F_MOLECULEDERIVATIVES, F_MOLECULES, F_MOLECULESECONDDERIVATIVES, &
      & F_NABTERMS, F_NAZIMUTHANGLES, F_NCLOUDSPECIES, F_NMODELSURFS, &
      & F_NO_DUP_MOL, F_NSCATTERINGANGLES, F_NSIZEBINS, F_PATHNORM, &
      & F_PHIWINDOW, F_POLARIZED, F_ReferenceMIF, F_REFRACT, F_SCANAVERAGE, &
      & F_SIGNALS, F_SKIPOVERLAPS, F_SPECIFICQUANTITIES, F_SPECT_DER, &
      & F_SWITCHINGMIRROR, F_TANGENTGRID, F_TEMP_DER, F_TOLERANCE, &
      & F_TRANSFORMMIFEXTINCTION, F_TRANSFORMMIFRHI, F_TSCATMIF, F_TYPE, &
      & F_USBLBLMOLECULES, F_USBPFAMOLECULES, F_useTSCAT, F_XSTAR, F_YSTAR
    use INTRINSIC, only: L_NONE, L_CLEAR, PHYQ_ANGLE, PHYQ_PROFILES
    use L2PC_M, only: BINSELECTORS, DEFAULTSELECTOR_LATITUDE, CREATEDEFAULTBINSELECTORS
    use MLSKINDS, only: R8
    use MLSL2OPTIONS, only: L2CFNODE, MLSMESSAGE
    use MLSMESSAGEMODULE, only: MLSMSG_ERROR, MLSMSG_WARNING
    use MLSNUMERICS, only: HUNT
    use MLSSIGNALS_M, only: SIGNALS
    use MOLECULES, only: L_CLOUDICE
    use MORETREE, only: GET_BOOLEAN, GET_FIELD_ID
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use PARSE_SIGNAL_M, only: PARSE_SIGNAL
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NODE_ID, NSONS, NULL_TREE, SUB_ROSA, SUBTREE
    use TREE_TYPES, only: N_ARRAY
    use VGRIDSDATABASE, only: VGRIDS

    integer, intent(in) :: NAME         ! The name of the config
    integer, intent(in) :: ROOT         ! of the forwardModel specification.
    !                                     Indexes either a "named" or
    !                                     "spec_args" vertex. Local variables
    logical, intent(in) :: GLOBAL       ! Goes into info%globalConfig

    integer :: B                        ! Index of a beta group
    logical, dimension(:), pointer :: Channels   ! From Parse_Signal
    integer :: DerivTree                ! Tree index of f_MoleculeDerivatives
    integer :: Expr_Units(2)            ! Units of value returned by EXPR
    integer :: Field                    ! Field index -- f_something
    integer :: Found                    ! Where something is found in a list
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: GSON                     ! Son of son
    integer :: I, J, K                  ! Subscript and loop inductor.
    integer :: LBLTrees(2)              ! Tree indices of f_[lu]sbLBLMolecules
    type(spectroParam_t), pointer :: LineStru(:)
    integer :: LineTrees(3)             ! Tree indices of f_line...
    integer :: Me = -1                  ! String index for trace
    integer :: MoleculeTree             ! Tree index of f_Molecules
    integer, dimension(:), pointer :: MyMolecules ! In a LineTree
    integer :: NELTS                    ! Number of elements of an array tree
    logical :: No_Dup_Mol               ! Duplicate molecules => error
    integer :: NumPFA, NumLBL           ! Numbers of such molecules in a beta group
    integer :: PFATrees(2)              ! Tree indices of f_[lu]sbPFAMolecules
    integer :: S                        ! Sideband index, for PFAtrees
    integer :: S1, S2                   ! Sideband Start/Stop (1..2, not -1..1).
    integer :: SecondDerivTree          ! Tree index of f_MoleculeSecondDerivatives
    integer :: SIDEBAND                 ! Returned from Parse_Signal
    integer, dimension(:), pointer :: SIGNALINDS ! From Parse_Signal
    character (len=80) :: SIGNALSTRING  ! E.g. R1A....
    integer :: Son                      ! Some subtree of root.
    integer :: STATUS                   ! From allocates etc.
    integer :: TANGENT                  ! Loop counter
    integer, pointer :: TempLBL(:), TempPFA(:) ! Used to separate LBL and PFA
    integer :: TheTree                  ! Either pfaTrees(s) or lblTrees(s)
    integer :: THISMOLECULE             ! Tree index.
    integer :: type                     ! Type of value returned by EXPR
    real (r8) :: Value(2)               ! Value returned by EXPR
    integer :: WANTED                   ! Which signal do we want?

    ! Nullify some pointers so allocate_test doesn't try to deallocate them.
    ! Don't initialize them with =>NULL() because that makes them SAVEd.

    nullify ( channels, signalInds )
    call NullifyForwardModelConfig ( info ) ! for Sun's rubbish compiler

    error = 0
    call trace_begin ( me, "ConstructForwardModelConfig", root, &
      & cond=toggle(gen) .and. levels(gen) > 0 )

    ! Set sensible defaults
    info%allLinesForRadiometer = .false.
    info%allLinesInCatalog = .false.
    info%anyLBL = .false.
    info%anyPFA = .false.
    info%atmos_der = .false.
    info%atmos_second_der = .false.
    info%cloud_der = l_none
    info%default_spectroscopy = .false.
    info%differentialScan = .false.
    info%do_1d = .false.
    info%do_baseline = .false.
    info%do_conv = .false.
    info%do_freq_avg = .false.
    info%do_path_norm = .false.
    info%forceFoldedOutput = .false.
    info%forceSidebandFraction = .false.
    info%frqTol = 600.0 ! MHz
    info%GenerateTScat = .false.
    info%globalConfig = global
    info%incl_cld = .false.
    info%ignoreHessian = .false.
    info%instrumentModule = 0
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
    info%refract = switchDetail(switches,'norf') < 0 ! Default .true.
    info%scanAverage = .false.
    info%sideBandStart = -1
    info%sideBandStop = 1
    info%skipOverlaps = .false.
    info%spect_der = .false.
    info%switchingMirror= .false.
    info%temp_der = .false.
    info%transformMIFextinction = .false.
    info%transformMIFRHI = .false.
    info%TScatMIF = 1
    info%useTScat = .false.
    info%where = root
    info%windowUnits = phyq_profiles
    info%xStar = 0
    info%yStar = 0

    got = .false.
    info%tolerance = -1.0 ! Kelvins, in case the tolerance field is absent
    lblTrees = null_tree
    lineTrees = null_tree
    no_dup_mol = .false.
    pfaTrees = null_tree
    do i = 2, nsons(root)
      son = subtree(i,root)
      L2CFNODE = son
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case ( f_allLinesForRadiometer )
        info%allLinesForRadiometer = get_boolean(son)
      case ( f_allLinesInCatalog )
        info%allLinesInCatalog = get_boolean(son)
      case ( f_atmos_der )
        info%atmos_der = get_boolean(son)
      case ( f_atmos_second_der )
        info%atmos_second_der = get_boolean(son)
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
      case ( f_frqTol )
        call expr ( subtree(2,son), expr_units, value, type )
        info%frqTol = value(1)
      case ( f_i_saturation )
        info%i_saturation = decoration(subtree(2,son))
      case ( f_ignoreHessian )
        info%ignoreHessian = get_boolean(son)
      case ( f_incl_cld )
        info%incl_cld = get_boolean(son)
      case ( f_integrationGrid )
        info%integrationGrid => vGrids(decoration(decoration(subtree(2,son))))
      case ( f_linearSideband )
        call expr ( subtree(2,son), expr_units, value, type )
        info%linearSideband = nint ( value(1) )
      case ( f_lineCenter )
        lineTrees(lineCenter) = son
      case ( f_lineWidth )
        lineTrees(lineWidth) = son
      case ( f_lineWidth_TDep )
        lineTrees(lineWidth_TDep) = son
      case ( f_lockBins )
        info%lockBins = get_boolean(son)
      case ( f_lsbLBLMolecules )
        LBLTrees(1) = son
      case ( f_lsbPFAMolecules )
        PFATrees(1) = son
      case ( f_module )
        info%instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_moleculeDerivatives )
        derivTree = son
      case ( f_molecules )
        moleculeTree = son
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
      case ( f_no_dup_mol )
        no_dup_mol = get_boolean(son)
      case ( f_nscatteringangles )
        call expr ( subtree(2,son), expr_units, value, type )
        info%NUM_SCATTERING_ANGLES = nint( value(1) )
        info%NUM_AZIMUTH_ANGLES = nint( value(1) )
      case ( f_nsizebins )
        call expr ( subtree(2,son), expr_units, value, type )
        info%NUM_SIZE_BINS = nint( value(1) )
      case ( f_pathNorm )
        info%do_path_norm = get_boolean(son)
      case ( f_phiWindow )
        call expr ( subtree(2,son), expr_units, value, type )
        info%phiWindow = value(1)
        if ( all ( expr_units(1) /= (/ PHYQ_Profiles, PHYQ_Angle /) ) ) &
          & call AnnounceError ( WrongUnitsForWindow, son )
        info%windowUnits = expr_units(1)
      case ( f_polarized )
        info%polarized = get_boolean(son)
      case ( f_referenceMIF )
        call expr ( subtree(2,son), expr_units, value, type )
        info%referenceMIF = nint(value(1))
      case ( f_refract )
        info%refract = get_boolean(son)
      case ( f_scanAverage )
        info%scanAverage = get_boolean(son)
      case ( f_moleculeSecondDerivatives )
        secondDerivTree = son
      case ( f_signals )
        info%noUsedChannels = 0
        allocate ( info%signals (nsons(son)-1), stat = status )
        if ( status /= 0 ) call announceError( AllocateError, root )
        info%signals%index = -1 ! Indicate error, covered up if successful
        call allocate_test ( info%signalIndices, nsons(son)-1, &
          & 'Info%SignalIndices', moduleName )
        do j = 1, nsons(son)-1
          gson = subtree(j+1,son)
          call get_string ( sub_rosa(gson), signalString, strip=.true.)
          call parse_Signal ( signalString, signalInds, &
            & tree_index=gson, sideband=sideband, channels=channels )
          if ( .not. associated(signalInds) ) then ! A parse error occurred
            error = max(error,1)
            cycle
          end if
          ! Later on choose the `right' one from the match
          ! For the moment choose first !????
          wanted=1
          info%signals(j) = signals(signalInds(wanted))
          info%signals(j)%sideband = sideband
          info%signalIndices(j) = signalInds(wanted)
          ! Don't hose channels in database, though shouldn't be an issue
          nullify ( info%signals(j)%channels )

          call allocate_Test ( info%signals(j)%channels, &
            & size(info%signals(j)%frequencies), 'info%signals%channels', &
            & ModuleName ) ! , lowBound=lbound(info%signals(j)%frequencies,1) )
          if ( associated(channels) ) then
            info%signals(j)%channels(1:lbound(channels,1)-1) = .false.
            info%signals(j)%channels(lbound(channels,1):ubound(channels,1)) = &
              channels
            info%noUsedChannels = info%noUsedChannels + count(channels)
            info%signals(j)%channels(ubound(channels,1)+1:) = .false.
          else
            info%signals(j)%channels = .true.
            info%noUsedChannels = info%noUsedChannels + &
              & size(info%signals(j)%channels)
          end if
          call deallocate_test ( channels, 'channels', ModuleName )
          call deallocate_test ( signalInds, 'signalInds', ModuleName )
        end do                          ! End loop over listed signals
        ! Make sure signal specifications make sense; get sideband Start/Stop
        call validateSignals
      case ( f_skipOverlaps )
        info%skipOverlaps = get_boolean(son)
      case ( f_specificQuantities )
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
      case ( f_transformMIFextinction )
        info%transformMIFextinction = get_Boolean(son)
      case ( f_transformMIFRHI )
        info%transformMIFRHI = get_Boolean(son)
      case ( f_TScatMIF )
        call expr ( subtree(2,son), expr_units, value, type )
        info%TScatMIF = nint(value(1))
      case ( f_type )
        info%fwmType = decoration(subtree(2,son))
      case ( f_usbLBLMolecules )
        LBLTrees(2) = son
      case ( f_usbPFAMolecules )
        PFATrees(2) = son
      case ( f_useTScat )
        info%useTScat = get_boolean(son)
      case ( f_xStar )
        info%xStar = decoration ( decoration ( subtree ( 2, son ) ) )
      case ( f_yStar )
        info%yStar = decoration ( decoration ( subtree ( 2, son ) ) )
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do ! i = 2, nsons(root)

    ! call dump(info)

    if ( ( got(f_lsbPFAMolecules) .and. got(f_lsbLBLMolecules) ) ) &
      & call announceError ( LBLandPFA, lblTrees(1) )
    if ( ( got(f_usbPFAMolecules) .and. got(f_usbLBLMolecules) ) ) &
      & call announceError ( LBLandPFA, lblTrees(2) )

    if ( ( .not. info%atmos_der ) .and. info%atmos_second_der ) then
      call announceError( Hess_notJac, root )
    end if

    s1 = (info%sidebandStart+3)/2; s2 = (info%sidebandStop+3)/2

    ! Now the molecule lists.
    nelts = 0    ! Total number of molecules, not counting group names
    if ( got(f_molecules) ) then
      ! Create the Beta_Group structure.  Temporarily assume all the
      ! molecules are LBL.  We'll separate them later, when processing
      ! the LBL or PFA Molecules trees.
      allocate ( info%beta_group(nsons(moleculeTree)-1), stat = status )
      if ( status /= 0 ) call announceError( AllocateError, moleculeTree )
      do b = 1, nsons(moleculeTree) - 1
        son = subtree(b+1,moleculeTree)
        info%beta_group(b)%group = node_id(son) == n_array
        if ( info%beta_group(b)%group ) then
          if ( info%fwmType /= l_full ) &
            & call announceError ( noBetaGroup, son )
          call allocate_test ( info%beta_group(b)%lbl(1)%molecules, nsons(son)-1, &
            & 'info%beta_group(b)%lbl(1)%molecules', moduleName )
          call allocate_test ( info%beta_group(b)%lbl(2)%molecules, nsons(son)-1, &
            & 'info%beta_group(b)%lbl(2)%molecules', moduleName )
          info%beta_group(b)%molecule = decoration(subtree(1,son)) ! group name
          j = nsons(son)
          if ( j < 2 ) call announceError ( badMoleculeGroup, son )
          nelts = nelts + j - 1
          do j = 2, j
            gson = subtree(j,son)
            if ( node_id(gson) == n_array ) then
              call announceError ( nested, gson )
            else
              info%beta_group(b)%lbl(1)%molecules(j-1) = decoration(gson)
            end if
          end do
          info%beta_group(b)%lbl(2)%molecules = info%beta_group(b)%lbl(1)%molecules
        else ! type checker guarantees a molecule name here
          nelts = nelts + 1
          info%beta_group(b)%molecule = decoration(son)
          call allocate_test ( info%beta_group(b)%lbl(1)%molecules, 1, &
            & 'info%beta_group(b)%lbl(1)%molecules', moduleName )
          call allocate_test ( info%beta_group(b)%lbl(2)%molecules, 1, &
            & 'info%beta_group(b)%lbl(2)%molecules', moduleName )
          info%beta_group(b)%lbl(1)%molecules(1) = info%beta_group(b)%molecule
          info%beta_group(b)%lbl(2)%molecules(1) = info%beta_group(b)%molecule
        end if
        do s = s1, s2 ! both sidebands
          if ( max(pfaTrees(s),lblTrees(s)) /= null_tree ) then
            ! Allocate PFA_Molecules temporarily for a "this is a PFA molecule" flag
            ! If it turns out to be an LBL molecules list, we'll swap them later.
            call allocate_test ( info%beta_group(b)%pfa(s)%molecules, &
              & size(info%beta_group(b)%lbl(s)%molecules), 'PFA Molecules', &
              & moduleName, fill=0 )
          else
            call allocate_test ( info%beta_group(b)%pfa(s)%molecules, 0, &
              & 'PFA Molecules', moduleName )
          end if
        end do ! s = s1, s2
      end do ! b = 1, nsons(moleculeTree) - 1
    else
      allocate ( info%beta_group(0), stat = status )
      if ( status /= 0 ) call announceError( AllocateError, moleculeTree )
    end if
    info%molecules => info%beta_group%molecule

    ! Announce duplicate molecules, but don't make it an error.
    ! (unless we have specifically disallowed them with no_dup_mol)
    ! All of the molecules are in ...%lbl(1) at this moment.
    do b = 1, size(info%beta_group)
      do i = 1, size(info%beta_group(b)%lbl(1)%molecules)
        ! First in the same group
        do j = i+1, size(info%beta_group(b)%lbl(1)%molecules)
          if ( info%beta_group(b)%lbl(1)%molecules(i) == &
            &  info%beta_group(b)%lbl(1)%molecules(j) ) then
            call announceError ( duplicateMolecule, moleculeTree, &
              what=info%beta_group(b)%lbl(1)%molecules(i), &
              warn=.not. no_dup_mol)
          end if
        end do
        ! Now in other groups
        do k = b+1, size(info%beta_group)
          do j = 1, size(info%beta_group(k)%lbl(1)%molecules)
            if ( info%beta_group(b)%lbl(1)%molecules(i) == &
              &  info%beta_group(k)%lbl(1)%molecules(j) ) then
              call announceError ( duplicateMolecule, moleculeTree, &
                what=info%beta_group(b)%lbl(1)%molecules(i), &
                warn=.not. no_dup_mol )
            end if
          end do
        end do
      end do
    end do

    ! Now the LBLMolecules or PFAMolecules lists
    info%cat_size = 0
    do s = s1, s2
      theTree = merge(lblTrees(s),pfaTrees(s),lblTrees(s)/=null_tree)
      if ( theTree /= null_tree ) then
        ! Verify that the LBL or PFA molecules are all listed molecules.
        ! Mark the ones that are LBL or PFA.
op:     do j = 2, nsons(theTree)
          son = subtree( j, theTree )
          if ( node_id(son) == n_array ) then
            call announceError ( NoArray, son )
            cycle
          end if
          thisMolecule = decoration( son )
          do b = 1, size(info%beta_group)
            do k = 1, size(info%beta_group(b)%lbl(s)%molecules)
              if ( thisMolecule == info%beta_group(b)%lbl(s)%molecules(k) ) then
                ! Duplicate in list is warning, not error.
!                 if ( info%beta_group(b)%pfa(s)%molecules(k) /= 0 ) &
!                   & call announceError ( PFATwice, son, warn=.true. )
                info%beta_group(b)%pfa(s)%molecules(k) = -1
                cycle op
              end if
            end do ! k = 1, size(info%beta_group(b)%lbl(s)%molecules)
          end do ! b = 1, size(info%beta_group)
          call announceError ( PFANotMolecule, son )
        end do op ! j = 2, nsons(theTree)
        ! Divide the LBL and PFA molecules into separate lists
        do b = 1, size(info%beta_group)
          tempLBL => info%beta_group(b)%lbl(s)%molecules
          tempPFA => info%beta_group(b)%pfa(s)%molecules
          if ( pfaTrees(s) /= null_tree ) then
            numLBL = count(tempPFA == 0)
            numPFA = size(tempPFA) - numLBL
          else
            numPFA = count(tempPFA == 0)
            numLBL = size(tempPFA) - numPFA
          end if
          nullify ( info%beta_group(b)%lbl(s)%molecules, &
            &       info%beta_group(b)%pfa(s)%molecules )
          call allocate_test ( info%beta_group(b)%lbl(s)%molecules, numLBL, &
            & 'LBL Molecules', moduleName )
          call allocate_test ( info%beta_group(b)%pfa(s)%molecules, numPFA, &
            & 'PFA Molecules', moduleName )
          if ( pfaTrees(s) /= null_tree ) then
            info%beta_group(b)%lbl(s)%molecules = pack(tempLBL,tempPFA == 0)
            info%beta_group(b)%pfa(s)%molecules = pack(tempLBL,tempPFA /= 0)
          else
            info%beta_group(b)%lbl(s)%molecules = pack(tempLBL,tempPFA /= 0)
            info%beta_group(b)%pfa(s)%molecules = pack(tempLBL,tempPFA == 0)
          end if
          call deallocate_test ( tempLBL, 'TempLBL', moduleName )
          call deallocate_test ( tempPFA, 'TempPFA', moduleName )
          if ( numLBL /= 0 ) info%anyLBL(s) = .true.
          if ( numPFA /= 0 ) info%anyPFA(s) = .true.
        end do ! b
      else
        info%anyLBL(s) = .true.
      end if

      ! Now the cat_index and isotope ratio fields
      info%cat_size(s) = 0
      do b = 1, size(info%beta_group)
        k = size(info%beta_group(b)%lbl(s)%molecules)
        info%cat_size(s) = info%cat_size(s) + k
        call allocate_test ( info%beta_group(b)%lbl(s)%cat_index, &
          &                  k, 'beta_group(b)%Cat_Index', moduleName )
        info%beta_group(b)%lbl(s)%cat_index = 0 ! in case somebody asks for a dump
        call allocate_test ( info%beta_group(b)%lbl(s)%ratio, &
          &                  k, 'LBL Ratio', moduleName )
        info%beta_group(b)%lbl(s)%ratio = 1.0
        call allocate_test ( info%beta_group(b)%pfa(s)%ratio, &
          &                  size(info%beta_group(b)%pfa(s)%molecules), &
          &                  'PFA Ratio', moduleName )
        info%beta_group(b)%pfa(s)%ratio = 1.0
      end do ! b = 1, size(info%beta_group)

    end do ! s = s1, s2

    ! Now the spectroscopy parameters.  They only make sense for LBL.
    do i = lineCenter, lineWidth_TDep
      ! Make a list of the molecules
      allocate ( lineStru(max(nsons(lineTrees(i))-1,0)), stat=status )
      if ( status /= 0 ) call announceError( AllocateError, lineTrees(i) )
      select case ( i )
      case ( lineCenter )
        info%lineCenter => lineStru
      case ( lineWidth )
        info%lineWidth => lineStru
      case ( lineWidth_TDep )
        info%lineWidth_TDep => lineStru
      end select
      if ( lineTrees(i) == null_tree ) cycle
      myMolecules => lineStru%molecule
      do j = 2, nsons(LineTrees(i))
        son = subtree( j, LineTrees(i) )
        if ( node_id(son) == n_array ) then
          call announceError ( NoArray, son )
          cycle
        end if
        myMolecules(j-1) = decoration( son )
        ! Look for duplicates
        do k = 2, j-1
          if ( myMolecules(k-1) == thisMolecule ) &
            & call announceError ( lineParamTwice, son )
        end do  ! k = 2, j-1
      end do ! j = 2, nsons(LineTrees(i))
      ! Check that molecules are in some LBL's molecule list
      do s = s1, s2 ! Sideband
        do b = 1, size(info%beta_group)
          if ( size(info%beta_group(b)%lbl(s)%molecules) > 0 ) then
            do found = 1, size(myMolecules)
              if ( abs(myMolecules(found)) == info%beta_group(b)%lbl(s)%molecules(1) ) then
                lineStru%beta(s) = b ! Store beta group index
                myMolecules(found) = -abs(myMolecules(found)) ! Mark it as used
                info%beta_group(b)%lbl(s)%spect_der_ix(i) = found
                exit
              end if
            end do ! found
          end if
        end do ! b = 1, size(info%beta_group)
      end do ! s = s1, s2
      ! Check for unused ones, make them positive again
      do j = 1, size(myMolecules)
        if ( myMolecules(j) > 0 ) & ! not used, so not in any LBL list
          & call announceError ( LineNotMolecule, subtree(j+1,lineTrees(i)) )
        myMolecules(j) = abs(myMolecules(j))
      end do
    end do ! i = lineCenter, lineWidth_TDep

    info%moleculeDerivatives => info%beta_group%derivatives
    info%moleculeDerivatives = .false.

    info%moleculeSecondDerivatives => info%beta_group%secondDerivatives
    info%moleculeSecondDerivatives = .false.

    ! Get info%moleculeDerivatives 
    if ( got(f_moleculeDerivatives) ) then
      if ( .not. associated(info%molecules) ) &
        & call announceError ( derivSansMolecules, derivTree )
      
      if ( .not. info%atmos_der ) &
        & call announceError ( firstSansFirst1, derivTree )

      do j = 2, nsons(derivTree)
        thisMolecule = decoration( subtree( j, derivTree ) )
        if ( .not. any(info%molecules == thisMolecule) ) &
          & call announceError ( derivSansMolecules, subtree(j,derivTree) )
        if ( got(f_molecules) ) where ( info%molecules == thisMolecule ) &
          & info%moleculeDerivatives = .true.
      end do                          ! End loop over listed species
    else if ( info%atmos_der ) then
      call announceError ( firstSansFirst2, derivTree )
    end if

    ! Get info%moleculeSecondDerivatives 
    if ( got(f_moleculeSecondDerivatives) ) then
      if ( .not. got(f_moleculeDerivatives) ) &
        & call announceError ( secondSansFirst, secondDerivTree )

      if ( .not. info%atmos_second_der ) &
        & call announceError ( secondSansSecond1, secondDerivTree )

      do j = 2, nsons(secondDerivTree)
        thisMolecule = decoration( subtree( j, secondDerivTree ) )
        if ( .not. any(info%molecules == thisMolecule) ) &
          & call announceError ( derivSansMolecules, subtree(j,secondDerivTree) )
        if ( got(f_molecules) ) where ( info%molecules == thisMolecule ) &
          & info%moleculeSecondDerivatives = .true.
      end do                          ! End loop over listed species
    else if ( info%atmos_second_der ) then
       call announceError ( secondSansSecond2, secondDerivTree )
    end if

    if ( (.not. got(f_moleculeSecondDerivatives)) .and. info%atmos_second_der ) then
      info%moleculeSecondDerivatives = info%moleculeDerivatives
    end if

    ! Now some more error checking
    ! MIFExtinction transformation needs signals
    if ( ( info%transformMIFextinction .or. info%transformMIFRhi ) &
         & .and. .not. associated(info%signals) ) &
      & call announceError ( MIFTransformation_signals, root )

    ! If any PFA, can't do polarized
    if ( any(info%anyPFA(s1:s2)) .and. info%polarized ) &
      & call announceError ( noPolarizedAndPFA, root )

    if ( ( info%xStar == 0 ) .neqv. ( info%yStar == 0 ) ) &
      & call AnnounceError ( NeedBothXYStar, root )

    select case ( info%fwmType )
    case ( l_full, l_hybrid )
      info%isRadianceModel = .true.
      if ( .not. ( got(f_molecules) .and. &
        & all(got( (/ f_signals, f_integrationGrid, f_tangentGrid /) )) ) ) &
        & call AnnounceError ( IncompleteFullFwm, root )

      if ( .not. associated(info%integrationGrid%surfs) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        &   'How can the integration grid not be associated?' )

      ! Now identify the Earth's surface in the tangent grid
      call Hunt ( info%tangentGrid%surfs(:,1), info%integrationGrid%surfs(1,1), &
        &  info%surfaceTangentIndex )

      ! Ensure that points in tangentGrid at and above the surface are a subset
      ! of integration grid
      if ( .not. associated(info%tangentGrid,info%integrationGrid) ) then
        do tangent = info%surfaceTangentIndex, info%tangentGrid%noSurfs
          if ( all ( abs( info%tangentGrid%surfs(tangent,1) - &
            & info%integrationGrid%surfs(:,1) ) > 1.0e-4 ) ) &
            & call AnnounceError ( TangentNotSubset, root )
        end do
      end if

      if ( info%sidebandStart == 1 .and. info%anyPFA(1) .or. &
         & info%sidebandStop == -1 .and. info%anyPFA(2) ) &
         & call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & 'Signal is SSB, but PFA is requested for the other sideband' )

      ! Cannot specify allLinesInCatalog and polarized
      if ( info%allLinesInCatalog .and. info%polarized ) &
        & call AnnounceError ( PolarizedAndAllLines, root )

      if ( info%incl_cld ) call checkCloud

    case ( l_cloudfull ) ! full cloud forward model

      info%isRadianceModel = .true.
      if ( .not. info%default_spectroscopy ) call checkCloud

    case ( l_scan )
      info%isRadianceModel = .false.
      ! Add 1d/2d method later probably !??? NJL
      if ( any(got( (/ f_channels, f_frequency, f_lineCenter, f_lineWidth, &
        &              f_lineWidth_TDep, f_molecules, f_moleculeDerivatives, &
        &              f_signals /) )) .or. &
        & any( (/ info%atmos_der, info%do_conv, info%do_baseline, &
        &         info%do_freq_avg, info%do_1d, info%incl_cld, &
        &         info%temp_der /) ) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root, &
          & "channels, frequency, lineCenter, lineWidth, lineWidth_TDep, "// &
          & "molecules, moleculeDerivatives, signals" )

    case ( l_linear, l_polarLinear )
      info%isRadianceModel = .true.
      if ( .not. all(got( (/f_signals/) )) ) & ! Maybe others later
        & call AnnounceError ( IncompleteLinearFwm, root )
      if ( any(got( (/f_do_conv,f_do_freq_avg,f_do_1d,f_incl_cld,f_frequency /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root, &
        & "do_conv, do_freq_avg, do_1d, incl_cld, frequency" )

    case default
      info%isRadianceModel = .false.
    end select

    if ( any ( info%fwmType == (/ l_linear, l_polarLinear, l_hybrid /) ) ) then
      ! Make sure we get a default bin selector
      if ( .not. associated ( info%binSelectors ) ) then
        if ( .not. associated ( binSelectors ) ) call CreateDefaultBinSelectors
        call Allocate_test ( info%binSelectors, 1, 'info%binSelectors', ModuleName )
        info%binSelectors = DefaultSelector_Latitude
      end if
    end if

    do i = 1, size(info%beta_group)
      do s = s1, s2 ! Sideband
        if ( any(info%beta_group(i)%lbl(s)%molecules == l_cloudIce) ) &
             & call announceError ( CloudLBL, root )
      end do
    end do

    if ( error /= 0 ) then
      call dump ( info, 'ConstructForwardModelConfig' )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'An error occured; see **** above' )
    end if

    call trace_end ( "ConstructForwardModelConfig", &
      & cond=toggle(gen) .and. levels(gen) > 0 )

  contains
    ! ...............................................  CheckCloud  .....
    subroutine CheckCloud
      use Intrinsic, only: Lit_Indices
      use Molecules, only: L_H2O, L_H2O_18, L_N2, L_N2O, L_O2, L_O3, L_O_18_O
      integer :: J
      logical :: GotH2O, GotO3
      gotH2O = .false.; gotO3 = .false.
      do j = 1, size(info%molecules) - 1
        select case ( info%molecules(j) )
        case ( l_h2o )
          gotH2O = .true.
        case ( l_o3 )
          gotO3 = .true.
        case ( l_n2o )
        case default
          select case ( info%molecules(j) )
          case ( l_n2, l_o2, l_h2o_18, l_o_18_o)
            call announceError ( cloudHas, root, &
              & what=lit_indices( info%molecules(j) ), warn=.true. )
          case default
            call announceError ( cloudNot, root, &
              & what=lit_indices( info%molecules(j) ) )
          end select
        end select
      end do
      if ( .not. GotH2O .or. .not. GotO3 ) &
        & call announceError ( cloudNeeds, root )
    end subroutine CheckCloud

    ! ..........................................  ValidateSignals  .....
    subroutine ValidateSignals
      use MLSSignals_m, only: GetSidebandStartStop
      use MLSFillValues, only: ESSENTIALLYEQUAL

      ! Make sure all the signals we're dealing with are same module,
      ! radiometer and sideband.
      if ( any( info%signals%sideband /= info%signals(1)%sideband ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  "Can't have mixed sidebands in forward model config" )
      if ( .not. all ( EssentiallyEqual ( info%signals%lo, info%signals(1)%lo ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        &  "Can't have mixed radiometers in forward model config" )

      ! Think about sidebands
      call getSidebandStartStop ( &
        & info%signals(1), info%sidebandStart, info%sidebandStop )

    end subroutine ValidateSignals

  end function ConstructForwardModelConfig

  ! ------------------------------FillFwdModelTimings  -----
  subroutine FillFwdModelTimings ( timings, FWModelConfig, which )
  !  Fill and return an array of time, mean, std_dev for timing FullforwardModel

    use ForwardModelConfig, only: ForwardModelConfig_T

    ! Dummy argument
    !real, pointer :: timings(:)
    double precision, dimension(:) :: timings
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    character(len=*), intent(in) :: which    ! 'fwdTiming', 'mean, or 'stdDev'

    ! Local variables
    real :: mean_sqDelta, meanTimes, tmp_mean
    integer :: i

    if (which == 'fwdTiming') then
      do i =1, size(FWModelConfig)
        timings(i) = FWModelConfig(i)%sum_DeltaTime
      end do
    end if

    if (which == 'mean') then
      do i =1, size(FWModelConfig)
        if (FWModelConfig(i)%sum_DeltaTime == 0.0 ) then
          timings(i) = 0.0
        else
          timings(i) = FWModelConfig(i)%sum_DeltaTime/FWModelConfig(i)%Ntimes
        end if
      end do
     end if

    if (which == 'stdDev') then
      do i =1, size(FWModelConfig)
        tmp_mean = FWModelConfig(i)%sum_DeltaTime/FWModelConfig(i)%Ntimes
        mean_sqDelta =  FWModelConfig(i)%sum_squareDeltaTime / &
          & FWModelConfig(i)%Ntimes
        if (FWModelConfig(i)%Ntimes <= 1) then
          meanTimes = 1.0
        else
          meanTimes = FWModelConfig(i)%Ntimes / (FWModelConfig(i)%Ntimes - 1)
        end if
        if (FWModelConfig(i)%sum_DeltaTime == 0.0 .AND. &
            & FWModelConfig(i)%sum_squareDeltaTime == 0.0 ) then
          timings(i) = 0.0
        else
          timings(i) = sqrt(abs(meanTimes * (mean_sqDelta - (tmp_mean * tmp_mean))))
        end if
      end do
    end if

  end subroutine FillFwdModelTimings

  ! ------------------------------ShowFwdModelNames  -----
  function ShowFwdModelNames ( FWModelConfig ) result (fwdNames)

  !  Fill and return an array of forward Model Names

    use ForwardModelConfig, only: ForwardModelConfig_T
    use MLSStringLists, only: catLists
    use String_Table, only: GET_STRING

    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    character(len=2000) :: fwdNames

    character(len=30) :: thisNames
    integer :: i

    fwdNames = ' '
    do i =1, size(FWModelConfig)
       call get_string ( FWModelConfig(i)%name, thisNames )
       fwdNames = catLists(trim(fwdNames), trim(thisNames), ' ')
       ! print *,' fwdNames: ', trim(fwdNames)
    end do

  end function ShowFwdModelNames

  ! ------------------------------------  PrintForwardModelTiming  -----
  subroutine PrintForwardModelTiming ( FWModelConfig )
  !  Print mean, std_dev for timing FullforwardModel

    use ForwardModelConfig, only: ForwardModelConfig_T
    use Output_m, only: BLANKS, Output
    use String_Table, only: GET_STRING

    ! Dummy argument
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig

    ! Local variables
    real :: mean, mean_sqDelta, meanTimes, std_dev
    character(len=25) :: thisName
    integer :: i

    call output ( "======= printForwardModelTiming =========", advance='yes')
    call output ( " ", advance = 'yes')
    call output ( "Name", advance='no')
    call blanks (18, advance ='no')
    call output ("| Invocation ", advance = 'no')
    call output ( "| Mean_time / s ", advance='no')
    call output ("| St. dev. / s", advance='yes')
    call output ("-------------------------------------------&
         &-----------------------------", advance= 'yes')
    do i =1, size(FWModelConfig)

       mean = FWModelConfig(i)%sum_DeltaTime/FWModelConfig(i)%Ntimes
       mean_sqDelta =  FWModelConfig(i)%sum_squareDeltaTime / &
                & FWModelConfig(i)%Ntimes
       if (FWModelConfig(i)%Ntimes <= 1) then
          meanTimes = 1.0
       else
          meanTimes = FWModelConfig(i)%Ntimes / (FWModelConfig(i)%Ntimes - 1)
       end if
       std_dev = sqrt(abs(meanTimes * (mean_sqDelta - (mean*mean))))
       call get_string ( FWModelConfig(i)%name, thisName )
       call output ( thisName, advance='no')
       call output ( FWModelConfig(i)%Ntimes, format = '(i6)',advance = 'no' )
       call blanks (8, advance = 'no')
       if (FWModelconfig(i)%Ntimes == 0) then
          mean = 0.0
          std_dev = 0.0
       end if
       call output ( mean, format='(f8.2)', advance = 'no' )
       call blanks (6, advance = 'no')
       call output ( std_dev, format='(f8.2)', advance = 'yes' )
       call output ( " ", advance = 'yes')
       call resetForwardModelTiming ( FWModelConfig(i))
     end do
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
  subroutine AnnounceError ( Code, where, extraMessage, what, warn )

    use Intrinsic, only: Lit_Indices
    use MoreTree, only: StartErrorMessage
    use Output_M, only: Output
    use String_Table, only: Display_String
    use Tree, only: Decoration

    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: where      ! Where in the tree did the error occur?
    character (LEN=*), optional, intent(in) :: extraMessage
    integer, intent(in), optional :: What ! Optional extra, usually string index
    logical, optional, intent(in) :: Warn ! Warning if TRUE
    ! Internal variables
    logical :: onlyWarn
    ! Executable
    onlyWarn = .false.
    if ( present(warn) ) onlyWarn = warn
    if ( .not. onlyWarn ) error = max(error,1)
    ! if ( .not. present(warn) ) error = max(error,1)
    call startErrorMessage ( where )
    call output ( ' ForwardModelSupport complained: ' )
    select case ( code )
    case ( AllocateError )
      call output ( 'Allocation error.', advance='yes' )
    case ( BadBinSelectors )
      call output ('Cannot have a fieldAzimuth binSelector for polarlinear model', &
        & advance='yes' )
    case ( BadHeightUnit )
      call output ( 'Inappropriate units for height in binSelector', &
        & advance='yes' )
    case ( BadMoleculeGroup )
      call output ( 'A molecule group has to have a name and a molecule', &
        & advance='yes' )
    case ( BadQuantityType )
      call output ( 'Bin Selectors cannot apply to this quantity type', &
        & advance='yes' )
    case ( CloudHas )
      call display_string ( what, before='Cloud forward model internally has the ' )
      call output ( ' molecule', advance='yes' )
    case ( CloudLBL )
      call output ( 'Cloud_A and Cloud_S cannot be LBL', advance='yes' )
    case ( CloudNeeds )
      call output ( 'Cloud forward model needs both H2O and O3 molecules', &
        & advance='yes' )
    case ( CloudNot )
      call display_string ( what, before='Cloud forward model cannot accept the ' )
      call output ( ' molecule', advance='yes' )
    case ( DerivSansMolecules )
      call output ( 'Derivative(s) requested for molecule(s) not specified.', &
        & advance='yes')
    case ( DuplicateMolecule )
      call display_string ( lit_indices(what), before='Duplicate molecule ', &
        & advance='yes' )
    case ( FirstSansFirst1 )
      call output ('moleculeDerivatives IS present, but atmos_der is NOT present', &
        & advance='yes' )
    case ( FirstSansFirst2 )
      call output ('atmos_der IS present, but moleculeDerivatives is NOT present', &
        & advance='yes' )
    case ( Hess_notJac )
      call output ('Atmospheric Jacobian not present, while atmospheric Hessian present', &
        & advance='yes' )
    case ( IncompleteBinSelectors )
      call output ('Must have some binSelectors for the polarlinear model',advance='yes' )
    case ( IncompleteFullFwm )
      call output ('Incomplete full foward model specification',advance='yes' )
    case ( IncompleteLinearFwm )
      call output ( 'Incomplete linear foward model specification', &
        & advance='yes' )
    case ( LBLandPFA )
      call output ( 'Cannot have both LBL and PFA molecule lists for one sideband', &
        & advance='yes' )
    case ( IrrelevantFwmParameter )
      call output ( 'Irrelevant parameter for this forward model type', &
        & advance='yes' )
    case ( LinearSidebandHasUnits )
      call output ( 'Irrelevant units for this linear sideband', &
        & advance='yes' )
    case ( LineNotMolecule )
      call display_string ( lit_indices(decoration(where)), &
        before='Spectral parameter requested for ' )
      call output ( ' but it is not an LBL molecule group', advance='yes' )
    case ( LineParamTwice )
      call display_string ( lit_indices(decoration(where)), before='Molecule '  )
      call output ( ' listed twice for spectral parameter', advance='yes' )
    case ( MIFTransformation_signals )
      call output ( 'MIF transformation needs signals', advance='yes' )
    case ( NeedBothXYStar )
      call output ( 'X/YStar must either be both present or both absent', &
        & advance='yes' )
    case ( Nested )
      call output ( 'Group within group is not allowed', advance='yes' )
    case ( NoArray )
      call output ( 'Nested array not allowed here', advance='yes' )
    case ( NoBetaGroup )
      call output ( 'Beta grouping allowed only for full clear-sky model', &
        & advance='yes' )
    case ( NoMolecule )
      call output ( 'A bin selector of type vmr must have a molecule', &
        & advance='yes' )
    case ( NoPolarizedAndPFA )
      call output ( "Sorry, don't yet know how to combine polarized and PFA", &
        & advance='yes' )
    case ( PFANotMolecule )
      call output ( 'LBL or PFA requested for ' )
      call display_string ( lit_indices(decoration(where)) )
      call output ( ' but it is not in the molecule list', advance='yes' )
    case ( PFATwice )
      call display_string ( lit_indices(decoration(where)), &
        & before='Molecule ' )
      call output ( ' listed twice for LBL or PFA', advance='yes' )
    case ( PolarizedAndAllLines )
      call output ( 'Cannot specify both polarized and allLinesInCatalog', &
        & advance='yes' )
    case ( SecondSansFirst )
      call output ('Second requests molecule with no first derivative', &
        & advance='yes' )
    case ( SecondSansSecond1 )
      call output ('moleculeSecondDerivatives IS present, but atmos_second_der is NOT present', &
        & advance='yes' )
    case ( SecondSansSecond2 )
      call output ('atmos_second_der IS present, but moleculeSecondDerivatives is NOT present', &
        & advance='yes' )
    case ( TangentNotSubset )
      call output ('Non subsurface tangent grid not a subset of integration grid', &
        & advance='yes' )
    case ( TooManyHeights )
      call output ( 'Bin Selectors can only refer to one height range', &
        & advance='yes' )
    case ( TooManyCosts )
      call output ( 'Bin Selectors can only have one cost', &
        & advance='yes' )
    case ( WrongUnitsForWindow )
      call output ( 'PhiWindow must be in degrees or profiles', &
        & advance='yes' )
    case default
      call output ( '(no specific description of this error)', advance='yes' )
    end select
    if ( present(extraMessage) ) call output ( extraMessage, advance='yes' )
    if ( .not. onlyWarn ) &
      & call output ( '(Set to stop due to error)', advance='yes' )
  end subroutine AnnounceError

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ForwardModelSupport

! $Log$
! Revision 2.173  2013/08/30 02:45:39  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.172  2013/08/23 00:32:47  pwagner
! Initialize transformMIFRHI
!
! Revision 2.171  2013/08/16 02:34:46  vsnyder
! Remove Model_Plane_MIF
!
! Revision 2.170  2013/08/09 01:03:59  vsnyder
! Add ReferenceMIF component
!
! Revision 2.169  2013/07/25 00:23:41  vsnyder
! Replace TransformRHI with TransformMIFRHI
!
! Revision 2.168  2013/07/19 01:19:46  vsnyder
! Add TransformRHI field
!
! Revision 2.167  2013/07/12 23:44:28  vsnyder
! Move units checking to type checker
!
! Revision 2.166  2013/07/12 23:25:28  vsnyder
! Bugus checkin: Remove unreferenced error messages
!
! Revision 2.165  2013/06/12 02:37:14  vsnyder
! Cruft removal
!
! Revision 2.164  2012/08/16 18:07:43  pwagner
! Exploit level 2-savvy MLSMessage
!
! Revision 2.163  2012/08/01 00:10:02  pwagner
! Prevent references to undefined cat_sizes
!
! Revision 2.162  2012/05/01 22:22:58  vsnyder
! Set IsRadianceModel component
!
! Revision 2.161  2012/03/28 00:56:49  vsnyder
! Move check for signals with MIF extinction from Wrappers to Support
!
! Revision 2.160  2012/03/07 02:13:18  vsnyder
! Add transformMIFextinction switch
!
! Revision 2.159  2011/08/20 02:03:34  vsnyder
! Remove unused use name
!
! Revision 2.158  2011/07/29 02:00:44  vsnyder
! Remove TScatMolecules and TScatMoleculeDerivatives.  Make CloudIce a
! molecule.  Remove Cloud_A and Cloud_S.
!
! Revision 2.157  2011/05/10 17:11:28  pwagner
! Corrected goof that changed info%refract default
!
! Revision 2.156  2011/05/09 18:09:30  pwagner
! Converted to using switchDetail
!
! Revision 2.155  2011/03/31 19:47:54  vsnyder
! Validate signals when the signals field is processed, so that sidebands
! are set.  Only check that cloud_a and cloud_s are PFA for the desired
! sidebands.
!
! Revision 2.154  2011/01/29 00:52:51  vsnyder
! Allow PFA without frequency averaging
!
! Revision 2.153  2010/11/08 19:26:22  pwagner
! Slight improvement of unhelpful error message
!
! Revision 2.152  2010/08/27 06:20:47  yanovsky
! Added atmos_second_der, moleculeSecondDerivatives.
! Added error handling: Hess_notJac, SecondDerivTree, FirstSansFirst1,
! FirstSansFirst2, SecondSansFirst, SecondSansSecond1, SecondSansSecond2.
!
! Revision 2.151  2010/06/07 23:30:12  vsnyder
! Add TScatMolecules, TScatMoleculeDerivatives, Use_Tscat.  Change
! PhaseFrqTol to FrqTol.
!
! Revision 2.150  2010/04/30 01:53:26  vsnyder
! Add description to 'Irrelevant parameter' error message
!
! Revision 2.149  2010/03/26 23:16:12  vsnyder
! Add ignoreHessian field to forward model config
!
! Revision 2.148  2010/02/25 18:17:48  pwagner
! Removed outmoded ascii l2pc file support
!
! Revision 2.147  2010/02/09 16:22:32  pwagner
! Give info%GenerateTScat an initial value
!
! Revision 2.146  2010/01/22 01:00:20  vsnyder
! Require Cloud_A and Cloud_S "molecules" to be handled as PFA, not LBL.
!
! Revision 2.145  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.144  2009/04/16 21:56:43  pwagner
! Added print to not_used_here
!
! Revision 2.143  2008/09/30 22:37:51  vsnyder
! Use nint(value) for TScatMIF
!
! Revision 2.142  2008/08/21 23:42:46  vsnyder
! Remove GenerateTScat from ForwardModel; use TScat on Sids
!
! Revision 2.141  2008/07/30 19:08:00  vsnyder
! Add PhaseFrqTol field to ForwardModel spec
!
! Revision 2.140  2008/06/06 22:52:53  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.139  2008/05/20 00:28:06  vsnyder
! Process MieTables field from ForwardModelGlobal
!
! Revision 2.138  2008/05/02 00:33:30  vsnyder
! Delete unused symbol
!
! Revision 2.137  2008/05/01 01:56:33  vsnyder
! Don't check grids subset if they're associated
!
! Revision 2.136  2007/11/07 03:10:48  vsnyder
! Add pathNorm field to forward model config
!
! Revision 2.135  2007/10/03 23:59:04  vsnyder
! Add 'where' for tracing
!
! Revision 2.134  2006/11/29 01:08:58  vsnyder
! Allocate spectroscopy derivative stuff with zero size if not used
!
! Revision 2.133  2006/07/20 23:39:53  vsnyder
! Remove unused declarations and USEs
!
! Revision 2.132  2006/06/03 01:46:20  vsnyder
! Remove no_dup_mol flag from config structure
!
! Revision 2.131  2006/05/11 19:37:32  pwagner
! Added option to disallow duplicate molecules
!
! Revision 2.130  2006/04/18 00:08:45  pwagner
! Allow abbreviated, pathless PFA files with PCF
!
! Revision 2.129  2006/03/30 18:58:47  vsnyder
! Allow duplicates in [LU]SBpfaMolecules and [LU]SBlblMolecules lists
!
! Revision 2.128  2006/03/22 02:23:46  vsnyder
! Add lsbLBLmolecules, useLBLmolecules
!
! Revision 2.127  2006/02/23 00:52:09  vsnyder
! Make sure instrumentModule has a value, report all errors before quitting
!
! Revision 2.126  2006/02/08 21:33:43  vsnyder
! Announce warning for duplicate molecules
!
! Revision 2.125  2006/02/07 00:19:07  vsnyder
! Prohibit conspiracy of polarized and PFA
!
! Revision 2.124  2006/01/11 01:56:48  vsnyder
! Allow PFA for sidebands not in signals, but ignore it -- of course
!
! Revision 2.123  2006/01/04 21:54:26  vsnyder
! 'norf' switch sets refraction default to false
!
! Revision 2.122  2005/12/29 01:11:08  vsnyder
! Add boolean 'refract' field to ForwardModel spec
!
! Revision 2.121  2005/12/22 21:08:15  vsnyder
! Require frequency averaging if both PFA and LBL
!
! Revision 2.120  2005/11/02 21:37:19  vsnyder
! Hoist some stuff out of FullForwardModel
!
! Revision 2.119  2005/10/14 23:14:28  vsnyder
! Require Frequency Averaging if any PFA molecules
!
! Revision 2.118  2005/09/16 23:39:07  vsnyder
! Add spect_der field to ForwardModel
!
! Revision 2.117  2005/09/03 01:20:41  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.116  2005/08/03 18:07:39  vsnyder
! Scan averaging, some spectroscopy derivative stuff
!
! Revision 2.115  2005/06/29 00:43:45  pwagner
! Utilizes new interface to ReadCompleteHDF5L2PCFile
!
! Revision 2.114  2005/06/14 20:41:55  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.113  2005/06/03 02:07:56  vsnyder
! New copyright notice, move Id to not_used_here to avoid cascades,
! get VGrids from VGridsDatabase instead of an argument, add SignalIndices
! component to config.
!
! Revision 2.112  2005/05/27 17:56:07  vsnyder
! Access Source_Ref in the right place
!
! Revision 2.111  2005/05/26 22:35:48  vsnyder
! Add PFAFiles field to ForwardModelGlobal
!
! Revision 2.110  2005/03/28 20:29:09  vsnyder
! Better error checking and reporting, some PFA stuff
!
! Revision 2.109  2005/02/17 02:34:41  vsnyder
! Fix a blunder in PFA setup
!
! Revision 2.108  2005/02/16 23:16:24  vsnyder
! Revise data structures for split-sideband PFA
!
! Revision 2.107  2005/01/27 21:23:13  vsnyder
! Inching toward separate LSB and USB PFA
!
! Revision 2.106  2004/12/28 00:22:34  vsnyder
! Remove unused declarations
!
! Revision 2.105  2004/12/13 20:15:40  vsnyder
! Add filter file names to string table
!
! Revision 2.104  2004/11/16 02:56:01  vsnyder
! Dump imputed to wrong requestor
!
! Revision 2.103  2004/11/05 19:39:06  vsnyder
! Moved some stuff here from fwdmdl/Get_Species_Data
!
! Revision 2.102  2004/11/04 03:42:33  vsnyder
! Provide for both LBL_Ratio and PFA_Ratio in beta_group
!
! Revision 2.101  2004/11/01 20:27:55  vsnyder
! Reorganization of representation for molecules and beta groups
!
! Revision 2.100  2004/10/13 02:24:56  livesey
! Bug fix, vGrids now pointer in case none defined
!
! Revision 2.99  2004/10/06 21:14:37  vsnyder
! Require Molecules or PFAMolecules
!
! Revision 2.98  2004/08/05 21:01:59  vsnyder
! Add sentinel at end of %molecules
!
! Revision 2.97  2004/08/04 23:19:57  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.96  2004/07/22 20:40:10  cvuu
! Add 2 subroutines FillFwdModelTimings and ShowFwdModelnames
!
! Revision 2.95  2004/07/17 02:27:24  vsnyder
! Better error message for PFA and non-PFA conflict
!
! Revision 2.94  2004/07/08 02:35:46  vsnyder
! Put all line-by-line molecules before PFA ones
!
! Revision 2.93  2004/06/12 00:42:35  vsnyder
! Make sure PFAMolecules is associated -- at least with zero size
!
! Revision 2.92  2004/05/18 01:25:14  vsnyder
! Hopefully finish pfaMolecules field support
!
! Revision 2.91  2004/05/01 04:05:50  vsnyder
! Added pfaMolecules -- but more work needed
!
! Revision 2.90  2004/03/22 18:24:40  livesey
! Added handling of AllLinesInCatalog flag.
!
! Revision 2.89  2004/03/05 18:32:55  livesey
! Bug fix, but commented out for the moment, as I need to get a test
! going.
!
! Revision 2.88  2003/10/29 00:44:33  livesey
! Made sure that default bin selectors also applied to hybrid model.
!
! Revision 2.87  2003/10/28 23:44:43  livesey
! Added initialization of forceFoldedOutput
!
! Revision 2.86  2003/10/15 16:59:25  pwagner
! Should allow null filenames if TOOLKIT
!
! Revision 2.85  2003/09/11 23:15:42  livesey
! Added handling of xStar / yStar arguments to l2pc models.
!
! Revision 2.84  2003/09/03 16:07:52  cvuu
! Add all of the printout stuff into routine PrintForwardModelTiming
!
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

