! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
  !=============================================================================

  ! Set up the forward model.  Interface from the retrieve step to the
  ! forward model.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use AntennaPatterns_m, only: AntennaPatterns
  use AntennaPatterns_m, only: Close_Antenna_Patterns_File, &
    & Open_Antenna_Patterns_File, Read_Antenna_Patterns_File
  use Declaration_Table, only: NUM_VALUE, RANGE
  use Dump_0, only: Dump
  use Expr_M, only: EXPR
  use ForwardModelConfig, only: AddForwardModelConfigToDatabase, Dump, &
    & ForwardModelConfig_T
  use ForwardModelIntermediate, only: ForwardModelIntermediate_T, ForwardModelStatus_T
  use SpectroscopyCatalog_m, only: Catalog_T, Lines, Catalog
  use FilterShapes_m, only: Close_Filter_Shapes_File, &
    & Open_Filter_Shapes_File, Read_Filter_Shapes_File, FilterShapes
  ! We're going to use lots of things from init_tables_module, so let's sort
  ! them into some sort of order
  ! First admin stuff
  use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
  ! Now fields
  use Init_Tables_Module, only: F_ANTENNAPATTERNS, F_ATMOS_DER, F_CHANNELS, &
    & F_DO_CONV, F_DO_FREQ_AVG, F_FILTERSHAPES, F_FREQUENCY, F_FRQGAP,&
    & F_INTEGRATIONGRID, F_L2PC, F_MOLECULES, F_MOLECULEDERIVATIVES, F_PHIWINDOW, &
    & F_POINTINGGRIDS, F_SIGNALS, F_SPECT_DER, F_TANGENTGRID, F_TEMP_DER, F_TYPE,&
    & F_MODULE
  ! Now literals
  use Init_Tables_Module, only: L_CHANNEL, L_EARTHREFL, L_ELEVOFFSET, L_FULL, &
    & L_LINEAR, L_LOSVEL, L_NONE, L_ORBITINCLINE, L_PTAN, L_RADIANCE,&
    & L_REFGPH, L_SCAN, L_SCGEOCALT, L_SIDEBANDRATIO, L_SPACERADIANCE, &
    & L_TEMPERATURE, L_VMR
  ! That's it for Init_Tables_Module
  use Lexer_Core, only: Print_Source
  use L2PC_M, only: OPEN_L2PC_FILE, READ_L2PC_FILE, CLOSE_L2PC_FILE
  use ManipulateVectorQuantities, only: FindClosestInstances
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error
  use MLSNumerics, only: Hunt
  use MLSSignals_m, only: ARESIGNALSSUPERSET, GetSignal, MaxSigLen, Signal_T, GetSignalName,&
    & MATCHSIGNAL, DUMP
  use Molecules, only: spec_tags
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use Parse_Signal_m, only: PARSE_SIGNAL
  use PointingGrid_m, only: Close_Pointing_Grid_File, &
    & Open_Pointing_Grid_File, Read_Pointing_Grid_File, PointingGrids
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Emit, Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named
  use Units, only: Deg2Rad, PHYQ_FREQUENCY
  use VectorsModule, only: GetVectorQuantityByType, ValidateVectorQuantity, &
    & Vector_T, VectorValue_T
  use VGridsDatabase, only: VGrid_T

  implicit none
  private
  public :: ConstructForwardModelConfig, FullForwardModel, &
    ForwardModelGlobalSetup

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Error codes

  integer, parameter :: AllocateError        = 1
  integer, parameter :: BadMolecule          = AllocateError + 1
  integer, parameter :: DefineSignalsFirst   = BadMolecule + 1
  integer, parameter :: DefineMoleculesFirst = DefineSignalsFirst + 1
  integer, parameter :: IncompleteFullFwm    = DefineMoleculesFirst + 1
  integer, parameter :: IncompleteLinearFwm  = IncompleteFullFwm + 1
  integer, parameter :: IrrelevantFwmParameter = IncompleteLinearFwm + 1
  integer, parameter :: TangentNotSubset     =  IrrelevantFwmParameter + 1
  integer, parameter :: PhiWindowMustBeOdd   = TangentNotSubset + 1
  integer, parameter :: FrqGapNotFrq         = PhiWindowMustBeOdd + 1

  integer :: Error            ! Error level -- 0 = OK

contains ! =====     Public Procedures     =============================

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root )
    ! Process the forwardModel specification to produce ForwardModelInfo.

    integer, intent(in) :: Root         ! of the forwardModel specification.
    !                                     Indexes a "spec_args" vertex.

    integer :: I                        ! Loop inductor, subscript
    integer :: Lun                      ! Unit number for reading a file
    character(len=255) :: FileName      ! Duh
    integer :: Son                      ! Some subtree of root.

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( 'ForwardModelGlobalSetup', root )

    ! "Root" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.
    ! Collect data from the fields.

    do i = 2, nsons(root)
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_antennaPatterns )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_antenna_patterns_file ( fileName, lun )
        call read_antenna_patterns_file ( lun )
        call close_antenna_patterns_file ( lun )
      case ( f_filterShapes )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_filter_shapes_file ( fileName, lun )
        call read_filter_shapes_file ( lun )
        call close_filter_shapes_file ( lun )
      case ( f_pointingGrids )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_pointing_grid_file ( fileName, lun )
        call read_pointing_grid_file ( lun )
        call close_pointing_grid_file ( lun )
      case ( f_l2pc )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_l2pc_file (fileName, lun)
        call read_l2pc_file ( lun )
        call close_l2pc_file ( lun )
      case default
        ! Can't get here if the type checker worked
      end select
    end do

    if ( toggle(gen) ) call trace_end ( 'ForwardModelGlobalSetup' )
  end subroutine ForwardModelGlobalSetup

  ! ------------------------------------------  ConstructForwardModelConfig  -----
  type (forwardModelConfig_T) function ConstructForwardModelConfig &
    & ( ROOT, VGRIDS ) result ( info )
    ! Process the forwardModel specification to produce ForwardModelConfig to add
    ! to the database
    use MLSSignals_M, only: Signals

    integer, intent(in) :: ROOT         ! of the forwardModel specification.
    !                                     Indexes either a "named" or
    !                                     "spec_args" vertex. Local variables
    type (vGrid_T), dimension(:), target :: vGrids ! vGrid database

    logical, dimension(:), pointer :: Channels   ! From Parse_Signal
    integer :: COMMONSIZE               ! Dimension
    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: J                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
    ! if any, else zero.
    integer :: NoChannelsSpecs          ! Number of channel specs we've had
    integer :: SIDEBAND                 ! Returned from Parse_Signal
    integer, dimension(:), pointer :: SIGNALINDS ! From Parse_Signal
    character (len=80) :: SIGNALSTRING  ! E.g. R1A....
    integer :: Son                      ! Some subtree of root.
    integer :: STATUS                   ! From allocates etc.
    integer :: TANGENT                  ! Loop counter
    integer :: THISMOLECULE             ! Tree index.
    integer :: type                     ! Type of value returned by EXPR
    integer :: Units(2)                 ! Units of value returned by EXPR
    real (r8) :: Value(2)               ! Value returned by EXPR
    integer :: WANTED                   ! Which signal do we want?

    ! Error message codes

    ! Nullify some pointers so allocate_test doesn't try to deallocate them.
    ! Don't initialize them with =>NULL() because that makes them SAVEd.

    nullify ( channels, signalInds )

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ConstructForwardModelConfig", root )
    if ( node_id(root) == n_named ) then
      name = subtree(1, root)
      key = subtree(2, root)
    else
      name = 0
      key = root
    end if

    ! Set sensible defaults
    info%do_conv = .false.
    info%do_freq_avg = .false.
    info%temp_der = .false.
    info%atmos_der = .false.
    info%spect_der = .false.
    info%phiwindow = 5
    info%frqGap = 0.0                   ! Default to everything

    noChannelsSpecs=0

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    got = .false.
    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      got(field) = .true.
      select case ( field )
      case ( f_type )
        info%fwmType = decoration(subtree(2,son))
      case ( f_atmos_der )
        info%atmos_der = get_boolean(son)
      case ( f_do_conv )
        info%do_conv = get_boolean(son)
      case ( f_do_freq_avg )
        info%do_freq_avg = get_boolean(son)
      case ( f_frqGap )
        call expr ( subtree(2,son), units, value, type )
        info%frqGap = value(1)
        if ( units(1) /= phyq_frequency ) &
          & call AnnounceError ( frqGapNotFrq, key )
      case ( f_module )
        info%instrumentModule = decoration(decoration(subtree(2,son)))
      case ( f_molecules )
        call allocate_test ( info%molecules, nsons(son)-1, "info%molecules", &
          & ModuleName )
        call allocate_test ( info%moleculeDerivatives, nsons(son)-1, &
          & "info%moleculeDerivatives", ModuleName )
        info%moleculeDerivatives = .false.
        do j = 1, nsons(son)-1
          info%molecules(j) = decoration( subtree( j+1, son ) )
        end do                          ! End loop over listed signals
      case ( f_moleculeDerivatives )
        if ( .not. associated(info%molecules) ) then
          call announceError( DefineMoleculesFirst, key)
        else
          do j = 1, nsons(son)-1
            thisMolecule = decoration( subtree( j+1, son ) )
            if ( .not. any(info%molecules == thisMolecule ) ) &
              & call announceError( BadMolecule, key )
            where ( info%molecules == thisMolecule )
              info%moleculeDerivatives = .true.
            end where
          end do                          ! End loop over listed signals
        end if
      case ( f_signals )
        allocate ( info%signals (nsons(son)-1), stat = status )
        if ( status /= 0 ) call announceError( AllocateError, key )
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
      case ( f_phiWindow )
        call expr ( subtree(2,son), units, value, type )
        info%phiWindow = nint( value(1) )
        if ( mod(info%phiWindow,2) /= 1 ) &
          & call AnnounceError ( phiWindowMustBeOdd, key )
      case ( f_spect_der )
        info%spect_der = get_boolean(son)
      case ( f_temp_der )
        info%temp_der = get_boolean(son)
      case ( f_integrationGrid )
        info%integrationGrid => vGrids(decoration(decoration(subtree(2,son))))
      case ( f_tangentGrid )
        info%tangentGrid => vGrids(decoration(decoration(subtree(2,son))))
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do ! i = 2, nsons(key)

    ! Now some more error checking
    select case ( info%fwmType )
    case ( l_full )
      if ( .not. all(got( (/ f_molecules, f_signals, f_integrationGrid, &
        & f_tangentGrid /) )) ) call AnnounceError ( IncompleteFullFwm, root )

      ! Now identify the Earth's surface in the tangent grid
      call Hunt(info%tangentGrid%surfs, info%integrationGrid%surfs(1), &
        &  info%surfaceTangentIndex)

      ! Ensure that points in tangentGrid at and above the surface are a subset
      ! of integration grid
      do tangent = info%surfaceTangentIndex, info%tangentGrid%noSurfs
        if ( .not. any ( abs(info%tangentGrid%surfs(tangent) - &
          & info%integrationGrid%surfs) < 1e-4) ) &
          & call AnnounceError ( TangentNotSubset, root )
      end do

      ! Check parameters needed only for linear/scan are not included
      !????
    case ( l_scan )
      ! Add 1d/2d method later probably !??? NJL
      if ( any(got( (/f_atmos_der, f_channels, f_do_conv, &
        & f_do_freq_avg, f_frequency, f_molecules, f_moleculeDerivatives, &
        & f_signals, f_spect_der, f_temp_der /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root )
    case ( l_linear)
      if ( .not. all(got( (/f_molecules, f_signals/) )) ) & ! Maybe others later
        & call AnnounceError ( IncompleteLinearFwm, root )
      if ( any(got( (/f_do_conv, f_do_freq_avg, f_frequency /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root )
    end select

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An error occured' )
    if ( toggle(gen) ) call trace_end ( "ConstructForwardModelConfig" )

  end function ConstructForwardModelConfig

  ! -----------------------------------------------  ForwardModel  -----
  subroutine FullForwardModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    &                       FwdModelOut, Ifm, FmStat, Jacobian )

    use GL6P, only: NG
    use MLSCommon, only: I4, R4, R8
    use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
    use ELLIPSE_M, only: ELLIPSE
    use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
    use GET_PATH_SPSFUNC_NGRID_M, only: GET_PATH_SPSFUNC_NGRID
    use REFRACTION_M, only: REFRACTION_CORRECTION
    use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
      PATH_DERIVATIVE, PATH_VECTOR_2D
    use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
    use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
    use GET_BETA_PATH_M, only: GET_BETA_PATH
    use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
    use RAD_TRAN_M, only: RAD_TRAN
    use RAD_TRAN_WD_M, only: RAD_TRAN_WD
    use FREQ_AVG_M, only: FREQ_AVG
    use CONVOLVE_ALL_M, only: CONVOLVE_ALL
    use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
    use D_LINTRP_M, only: LINTRP
    use D_HUNT_M, only: hunt_zvi => HUNT

    ! Dummy arguments --------------------------------------------------------

    ! From ForwardModelSetup
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: Ifm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Local parameters ---------------------------------------------------------

    character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

    ! Local variables ----------------------------------------------------------

    ! First the old stuff which we hope to get rid of or redefine
    integer(i4) :: brkpt, ch, frq_i, i, ier, ihi, ilo, j, k, lmax, m, maf, &
                   max_phi_dim, max_zeta_dim, mid, n, no_ele, no_tan_hts, &
                   ptg_i, si, Spectag

!    real(r4) :: K_STAR_ALL(25,20,mxco,mnp,Nptg)
!    type(k_matrix_info) :: K_star_info(20)


    type(path_derivative) :: K_temp_frq
    type(path_derivative), allocatable, dimension(:) :: K_atmos_frq

    type(path_beta), dimension(:,:), pointer :: Beta_path

    real(r8) :: Frq, Geod_lat, H_tan, Phi_tan, Rad

    ! This is the `legit stuff' we hope will stay; they are all pointers to
    ! VectorValue_T's containing vector quantities.
    type (VectorValue_T), pointer :: EARTHREFL     ! Earth reflectivity
    type (VectorValue_T), pointer :: ELEVOFFSET    ! Elevation offset quantity
    type (VectorValue_T), pointer :: F             ! An arbitrary species
    type (VectorValue_T), pointer :: LOSVEL        ! Line of sight velocity
    type (VectorValue_T), pointer :: ORBINCLINE    ! Orbital inclination (beta)
    type (VectorValue_T), pointer :: PTAN          ! PTAN quantity
    type (VectorValue_T), pointer :: FIRSTRADIANCE ! One radiance quantity to be filled
    type (VectorValue_T), pointer :: THISRADIANCE ! One radiance quantity to be filled
    type (VectorValue_T), pointer :: REFGPH        ! Reference GPH, (zRef and hRef)
    type (VectorValue_T), pointer :: SCGEOCALT     ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: SIDEBANDRATIO ! Sideband ratio for radiance
    type (VectorValue_T), pointer :: SPACERADIANCE ! Space radiance
    type (VectorValue_T), pointer :: TEMP          ! Temperature quantity

    integer :: CHANNEL                  ! Loop counter
    integer :: INSTANCE                 ! Loop counter
    integer :: MAFTINSTANCE             ! Temperature instance closest to this MAF
    integer :: MAXNOFREQS               ! Used for sizing arrays
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXSUPERSET              ! Max. value of superset
    integer :: MAXVERT                  ! Number of points in gl grid
    integer :: N2LVL                    ! Twice size of tangent grid
    integer :: NLVL                     ! Size of tangent grid
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOMIFS                   ! Number of minor frames
    integer :: NOSPECIES                ! Number of molecules we're considering
    integer :: NOUSEDCHANNELS           ! Number of channels to output
    integer :: NO_PHI_T                 ! No. of Temp. profiles in the chunk
    integer :: PHIWINDOW                ! Copy of forward model config%phiWindow
    integer :: SHAPEIND                 ! Index into filter shapes
    integer :: SIDEBANDSTART            ! Loop limit
    integer :: SIDEBANDSTEP             ! Loop step
    integer :: SIDEBANDSTOP             ! Loop limit
    integer :: SIGIND                   ! Loop counter
    integer :: SPECIE                   ! Loop counter
    integer :: STATUS                   ! From allocates etc.
    integer :: SURFACE                  ! Loop counter
    integer :: THISSIDEBAND             ! Loop counter
    integer :: TOTALSIGNALS             ! Used when hunting for pointing grids
    integer :: WHICHPATTERN             ! Index of antenna pattern
    integer :: WHICHPOINTINGGRID        ! Index of pointing grid
    integer :: WINDOWFINISH             ! Range of window
    integer :: WINDOWSTART              ! Range of window

    real (r8) :: CENTERFREQ             ! Of band
    real (r8) :: CENTER_ANGLE           ! For angles
    real (r8) :: R                      ! To convert the kind of output from
    !                                     Freq_Avg
    real (r8) :: THISRATIO              ! A sideband ratio

    integer, dimension(:), pointer :: CHANNELINDEX ! E.g. 1..25
    integer, dimension(:), pointer :: GRIDS ! Frq grid for each tan_press
    integer, dimension(:), pointer :: SUPERSET ! Result of AreSignalsSuperset
    integer, dimension(:), pointer :: USEDCHANNELS ! Array of indices used
    integer, dimension(:), pointer :: USEDSIGNALS ! Array of indices used

    logical :: FOUNDINFIRST                     ! Flag to indicate derivatives

    real(r8), dimension(:,:),   pointer :: D2X_DXDT    ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:,:,:), pointer :: DH_DT_PATH  ! (pathSize, Tsurfs, Tinstance)
    real(r8), dimension(:), allocatable :: Dum
    real(r8), dimension(:,:),   pointer :: DX_DT       ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:),     pointer :: FREQUENCIES ! Frequency points
    real(r8), dimension(:,:),   pointer :: I_STAR_ALL    ! (noMIFs,noChans)
    real(r4), dimension(:,:,:,:,:), pointer :: K_ATMOS ! (channel,Nptg,mxco,mnp,Nsps)
    real(r4), dimension(:,:,:,:), pointer :: K_TEMP    ! (channel,Nptg,mxco,mnp)
    real(r8), dimension(:),     pointer :: PTG_ANGLES  ! (no_tan_hts)
    real(r8), dimension(:,:),   pointer :: RADIANCES     ! (Nptg,noChans)
    real(r8), dimension(:),     pointer :: RadV
    real(r8), dimension(:,:),   pointer :: REF_CORR    ! (n2lvl, no_tan_hts)
    real(kind(k_temp_frq%values)), &
      & dimension(:), pointer :: TOAVG ! Stuff to be passed to frq.avg.
    real(r8), dimension(:),     pointer :: T_SCRIPT    ! (n2lvl)
    real(r8), dimension(:),     pointer :: TAU         ! (n2lvl)
    
    integer, dimension(1) :: WHICHPOINTINGGRIDASARRAY ! Result of minloc
    integer, dimension(1) :: WHICHPATTERNASARRAY ! Result of minloc

    type(path_vector), dimension(:), allocatable :: N_PATH    ! (No_tan_hts)

    ! dimensions of SPSFUNC_PATH are: (Nsps,No_tan_hts)
    type(path_vector), allocatable, dimension(:,:) :: SPSFUNC_PATH
    type(signal_t) :: FirstSignal
    type(signal_t) :: ThisSignal
    type(catalog_T), dimension(:), pointer :: My_Catalog

    ! Executable code --------------------------------------------------------

    if ( toggle(emit) ) call trace_begin ( 'ForwardModel' )

    ! Nullify a bunch of pointers so that Allocate_Test doesn't try to
    ! deallocate them.  We don't want them to be initialized NULL()
    ! because that makes them SAVEd.

    nullify (beta_path, channelIndex, d2x_dxdt, dh_dt_path, dx_dt, &
      & frequencies, grids, i_star_all, k_atmos, k_temp, my_Catalog, &
      & ptg_angles, radiances, radV, ref_corr, superset, t_script, &
      & tau, usedChannels, usedSignals )

    ! Identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    firstSignal = fwdModelConf%signals(1)

    ! Now make sure all the signals we're dealing with are same module,
    ! radiometer and sideband.
    if ( any( fwdModelConf%signals%sideband .ne. &
      & firstSignal%sideband ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed sidebands in forward model config")
    if ( any( fwdModelConf%signals%radiometer .ne. &
      & firstSignal%radiometer ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      &  "Can't have mixed radiometers in forward model config")

    ! Now from that we identify the radiance quantity we'll be outputting
    firstRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal=firstSignal%index, sideband=firstSignal%sideband )

    ! Identify the appropriate state vector components, save vmrs for later
    temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature )
    ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=firstSignal%instrumentModule )
    elevOffset => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_elevOffset, radiometer=firstSignal%radiometer )
    orbIncline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_orbitIncline )
    spaceRadiance => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_spaceRadiance )
    earthRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_earthRefl )
    refGPH => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_refGPH )
    losVel => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_losVel, instrumentModule=firstSignal%instrumentModule )
    scGeocAlt => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scGeocAlt )
    sidebandRatio => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_sidebandRatio, signal= firstSignal%index, noError=.true. )

    ! We won't seek for molecules here as we can't have an array of pointers.
    ! When we do want molecule i we would do something like
    ! vmr => GetVectorQuantityBytype (fwdModelIn, fwdModelExtra, &
    !   quantityType=l_vmr, molecule=fwdModelConf.molecules(i))

    ! Now we're going to validate the quantities we've been given, don't forget
    ! we already know what their quantityType's are as that's how we found them
    !, so we don't need to check that.
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )
    if ( .not. ValidateVectorQuantity(elevOffset, verticalCoordinate=(/l_none/), &
      & frequencyCoordinate=(/l_none/), noInstances=(/1/)) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'elevOffset' )
    ! There will be more to come here.

    noSpecies = size(fwdModelConf%molecules)

    !  Create a subset of the catalog composed only of those molecules to be
    !  used for this run

    maxNoFSurfs = 0
    do specie = 1, noSpecies
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=fwdModelConf%molecules(specie) )
      maxNoFSurfs = max(maxNoFSurfs, f%template%noSurfs)
    end do

    allocate ( My_Catalog(noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )

    do j = 1, noSpecies
      Spectag = spec_tags(fwdModelConf%molecules(j))
      do i = 1, Size(Catalog)
        if ( Catalog(i)%Spec_Tag == Spectag ) then
          My_Catalog(j) = Catalog(i)
          exit
        end if
      end do
    end do

    ! Get the max. dimension in zeta coeff. space and phi coeff. space
    ! (To be used later in rad_tran_wd, for automatic arrays asignement)
    max_phi_dim = 1
    max_zeta_dim = 1
    if ( FwdModelConf%temp_der ) then
      max_zeta_dim = temp%template%noSurfs
      max_phi_dim = temp%template%noInstances
    end if

    if ( fwdModelConf%atmos_der ) then
      do k = 1, noSpecies
        if ( fwdModelConf%moleculeDerivatives(k) ) then
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            &     quantityType=l_vmr, molecule=fwdModelConf%molecules(k))
          j = f%template%noInstances
          max_phi_dim = max(max_phi_dim,j)
          j = f%template%noSurfs
          max_zeta_dim = max(max_zeta_dim,j)
        end if
      end do
    end if

    ! Deal with fmStat%rows
    if ( present(Jacobian) .and. ( .not. associated (fmStat%rows) ) ) then
      call Allocate_test ( fmStat%rows, Jacobian%row%nb, 'fmStat%rows', &
        & ModuleName)
      fmStat%rows = .false.
    endif

    ! Get some dimensions that we'll use a lot
    noMAFs = firstRadiance%template%noInstances
    noMIFs = firstRadiance%template%noSurfs
    no_phi_t = temp%template%noInstances
    no_tan_hts = FwdModelConf%TangentGrid%nosurfs
    maxVert = 2 * (NG+1) * size(FwdModelConf%integrationGrid%surfs)
    nlvl=size(FwdModelConf%integrationGrid%surfs)
    n2lvl=2*nlvl
    phiWindow = FwdModelConf%phiWindow

    if ( toggle(emit) ) then
      print*,'Dimensions:'
      print*,'noMAFs:',noMAFs
      print*,'no_phi_t:',no_phi_t
      print*,'no_tan_hts:',no_tan_hts
      print*,'maxVert:',maxVert
      print*,'nlvl:',nlvl
      print*,'n2lvl:',n2lvl
      print*,'phiWindow:',phiWindow
      print*,'noSpecies:',noSpecies
      print*,'maxNoFSurfs:',maxNoFSurfs
      print*,'MAF:',fmStat%maf
    end if

    ! Work out which channels are used, also check we have radiances for them.
    noUsedChannels = 0
    do sigInd = 1, size(fwdModelConf%signals)
      thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
        & signal=fwdModelConf%signals(sigInd)%index, sideband=firstSignal%sideband )
      if ( .not. ValidateVectorQuantity(thisRadiance, minorFrame=.true.,&
        & frequencyCoordinate=(/l_channel/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'radiance' )
      noUsedChannels = noUsedChannels + &
        & count( fwdModelConf%signals(sigInd)%channels )
    end do
    call allocate_test ( usedChannels, noUsedChannels, &
      & 'usedChannels', ModuleName )
    call allocate_test ( usedSignals, noUsedChannels, &
      & 'usedSignals', ModuleName )
    channel = 1
    do sigInd = 1, size(fwdModelConf%signals)
      do i = 1, size(fwdModelConf%signals(sigInd)%frequencies)
        if (fwdModelConf%signals(sigInd)%channels(i)) then
          usedChannels(channel) = i
          usedSignals(channel) = sigInd
          channel = channel + 1
        end if
      end do
    end do

    ! --------------- Hydrostatic stuff ---------------------------------
    if ( fmStat%newHydros ) then

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.hydrostatic' )

      ! Now we're going to create the many temporary arrays we need
      allocate ( ifm%ndx_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'ndx_path' )
      allocate ( ifm%dhdz_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'dhdz_path' )
      allocate ( ifm%h_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'h_path' )
      allocate ( ifm%phi_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'phi_path' )
      allocate ( ifm%t_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'t_path' )
      allocate ( ifm%z_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'z_path' )

      allocate ( ifm%eta_phi(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'eta_phi' )

      allocate ( ifm%elvar(no_phi_t), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'elvar' )

      call allocate_test ( ifm%geoc_lat, no_phi_t, 'geoc_lat', ModuleName )
      call allocate_test ( ifm%e_rad, no_phi_t, 'e_rad', ModuleName )

      call allocate_test ( ifm%z_glgrid, maxVert/2, 'z_glgrid', ModuleName )
      call allocate_test ( ifm%h_glgrid, maxVert, no_phi_t, 'h_glgrid', &
                        &  ModuleName )
      call allocate_test ( ifm%t_glgrid, maxVert, no_phi_t, 't_glgrid', &
                        &  ModuleName )
      call allocate_test ( ifm%dh_dt_glgrid, maxVert, no_phi_t, &
        & temp%template%noSurfs,'dh_dt_glgrid', ModuleName )
      call allocate_test ( ifm%dhdz_glgrid, maxVert, no_phi_t, &
                        &  'dhdz_glgrid', ModuleName )
      call allocate_test ( ifm%tan_hts, no_tan_hts, no_phi_t, 'tan_hts', &
                        &  ModuleName )
      call allocate_test ( ifm%tan_temp, no_tan_hts, no_phi_t, 'tan_hts', &
                        &  ModuleName )
      call allocate_test ( ifm%tan_dh_dt, no_tan_hts, no_phi_t, &
                        & temp%template%noSurfs, 'tan_dh_dt', ModuleName )
      call Allocate_test( ifm%closestInstances, noMAFs, 'closestInstances', ModuleName)

      ! Setup for hydrostatic calculation
      call FindClosestInstances ( temp, firstRadiance, ifm%closestInstances )

      do i = 1, no_phi_t
        phi_tan = Deg2Rad*temp%template%phi(1,i)
        geod_lat= Deg2Rad*temp%template%geodLat(1,i)
        call geoc_geod_conv ( ifm%elvar(i), orbIncline%values(1,1), &
          &  phi_tan, geod_lat, ifm%geoc_lat(i), ifm%E_rad(i) )
      end do

      ! Now compute a hydrostatic grid given the temperature and refGPH
      ! information.
      call hydrostatic_model ( FwdModelConf%SurfaceTangentIndex, &
        &  no_phi_t, ifm%geoc_lat, 0.001*refGPH%values(1,:), &
        &  refGPH%template%surfs(1,1), &
        &  FwdModelConf%integrationGrid%surfs, &
        &  temp%template%surfs(:,1), temp%values, &
        &  ifm%z_glgrid, ifm%h_glgrid, ifm%t_glgrid, &
        &  ifm%dhdz_glgrid, ifm%dh_dt_glgrid, &
        &  FwdModelConf%TangentGrid%surfs, &
        &  ifm%tan_hts, ifm%tan_temp, ifm%tan_dh_dt, &
        &  ifm%gl_count, Ier )
      if ( ier /= 0 ) goto 99

      ! Now compute stuff along the path given this hydrostatic grid.
      call comp_path_entities ( temp, ifm%closestInstances, &
        &  FwdModelConf%integrationGrid%noSurfs, &
        &  temp%template%noSurfs, ifm%gl_count, ifm%ndx_path, ifm%z_glgrid, &
        &  ifm%t_glgrid, ifm%h_glgrid, ifm%dhdz_glgrid, ifm%tan_hts,        &
        &  no_tan_hts, ifm%z_path, ifm%h_path, ifm%t_path, ifm%phi_path,    &
        &  ifm%dhdz_path, ifm%eta_phi, no_phi_t,                            &
        &  temp%template%phi(1,:)*Deg2Rad, noMAFs, phiWindow, ifm%elvar, Ier )
      if ( ier /= 0 ) goto 99

      fmStat%newHydros = .false.

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.hydrostatic' )
    end if

    ! ------ End of hydrostatic setup stuff --------------------------
    ! ------ Begin main MAF Specific stuff ---------------------------

    ! Now allocate other stuff
    call allocate_test ( t_script, n2lvl, 't_srcipt', ModuleName )
    call allocate_test ( ref_corr, n2lvl, no_tan_hts, 'ref_corr', ModuleName )
    call allocate_test ( tau, n2lvl, 'tau', ModuleName )
    call allocate_test ( ptg_angles, no_tan_hts, 'ptg_angles', ModuleName )

    allocate ( k_atmos_frq(noSpecies), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_atmos_frq' )

    allocate ( n_path(No_tan_hts), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'n_path' )
    allocate ( spsfunc_path(noSpecies,No_tan_hts), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'spsfunc_path' )
    call allocate_test ( dx_dt, No_tan_hts, temp%template%noSurfs, &
      & 'dx_dt', ModuleName )
    call allocate_test ( d2x_dxdt, No_tan_hts, temp%template%noSurfs, &
      & 'd2x_dxdt', ModuleName )

    call allocate_test ( radiances, no_tan_hts, noUsedChannels, &
      & 'Radiances', ModuleName )
    call allocate_test ( i_star_all, noUsedChannels, noMIFs, &
      & 'i_star_all', ModuleName )

    ! Now set radiances to zero, forward model just adds in terms
    do i = 1, noUsedChannels
      thisRadiance = GetVectorQuantityByType (fwdModelOut, &
        & quantityType=l_radiance, &
        & signal=fwdModelConf%signals(usedSignals(i))%index, &
        & sideband=firstSignal%sideband )
      ch = usedChannels(i)
      do j = 1, noMIFs
        thisRadiance%values( ch + (j-1)*thisRadiance%template%noChans, fmStat%maf) = 0.0
      end do
    end do

    if ( fwdModelConf%signals(1)%sideband == 0 ) then
      if (.not. associated (sidebandRatio) ) &
        & call MLSMessage(MLSMSG_Error,ModuleName, &
        & "No sideband ratio supplied")
      sidebandStart = -1
      sidebandStop = 1
      sidebandStep = 2
    else
      sidebandStart = fwdModelConf%signals(1)%sideband
      sidebandStop = sideBandStart
      sidebandStep = 1
    endif

    ! ----------------- Begin loop over sidebands -----------------------
    do thisSideband = sidebandStart, sidebandStop, sidebandStep


      if ( toggle(emit) .and. levels(emit) > 0 ) then
        call trace_begin ( 'ForwardModel.sideband' )
        call output ( ' Doing sideband ' )
        call output ( thisSideband )
        call output ( ' (' ); call output ( sidebandStart )
        call output ( ', ' ); call output ( sidebandStop )
        call output ( ')', advance='yes' )
      end if

      ! Now code splits into two sections, one for when we're doing frequency
      ! averaging, and one when we're not.
      if ( fwdModelConf%do_freq_avg ) then ! --- Doing freq. avg. ---
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.FreqAvg' )

        call allocate_test ( superset, size(pointingGrids), &
          & 'superset', ModuleName )
        do i = 1, size(pointingGrids)
         superset(i) = AreSignalsSuperset ( pointingGrids(i)%signals, &
            & fwdModelConf%signals, sideband=thisSideband )
        end do
        if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "No matching pointing frequency grids." )
        
        maxSuperset = maxval ( superset )
        where ( superset < 0 )
          superset = maxSuperset + 1
        end where
        whichPointingGridAsArray = minloc ( superset )
        whichPointingGrid = whichPointingGridAsArray(1)
        call deallocate_test ( superset, 'superset', ModuleName )

        if ( toggle(emit) ) then
          call output ( 'Using pointing frequency grid: ' )
          call output ( whichPointingGrid, advance='yes' )
        end if
        
        ! Now we've identified the pointing grids.  Locate the tangent grid
        ! within it.
        call allocate_test ( grids, FwdModelConf%TangentGrid%nosurfs, &
          "Grids", ModuleName )
        call Hunt ( PointingGrids(whichPointingGrid)%oneGrid%height, &
          & FwdModelConf%TangentGrid%surfs, grids, allowTopValue=.true. )
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.FreqAvg' )

      else ! ------------------------- Not frequency averaging ---------
        
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_begin ( 'ForwardModel.NotFreqAvg' )

        call allocate_test ( frequencies,noFreqs, "frequencies", ModuleName )
        do channel = 1, noUsedchannels
          frequencies(channel) = &
            & fwdModelConf%signals(usedSignals(channel))%centerFrequency + &
            & fwdModelConf%signals(usedSignals(channel))% &
            &  frequencies(usedChannels(channel))
        end do
        select case ( thisSideband )
        case ( -1 )
          frequencies = firstSignal%lo - frequencies
        case ( +1 )
          frequencies = firstSignal%lo + frequencies
        case ( 0 )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Folded signal requested in forward model' )
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Bad value of signal%sideband' )
        end select
        noFreqs = noUsedChannels
        if ( toggle(emit) .and. levels(emit) > 0 ) &
          & call trace_end ( 'ForwardModel.NotFreqAvg' )
      end if

      ! ----------- Done the gnarly frequency stuff ----------

      maf=fmStat%maf

      ! Now work out what `window' we're inside.  This will need to be changed
      ! a bit in later versions to avoid the noMAFS==noTemp/f instances
      ! assertion

      mafTInstance = ifm%closestInstances(maf)

      windowStart  = max(1, mafTInstance - phiWindow/2)
      windowFinish = min(mafTInstance + phiWindow/2, no_phi_t)

      if ( toggle(emit) .and. levels(emit) > 0 ) then
        print *, 'Doing MAF: ', maf
        Print *, 'mafTInstance:',mafTInstance
        print *, 'WindowStart:',WindowStart
        print *, 'WindowFinish:',WindowFinish
      end if

      allocate ( k_temp(noUsedChannels, no_tan_hts, temp%template%noSurfs, &
        & windowStart:windowFinish), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'k_temp' )
      allocate ( k_atmos(noUsedChannels, no_tan_hts, maxNoFSurfs, &
        & windowStart:windowFinish, noSpecies), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'k_atmos' )

      ! Compute the specie function (spsfunc) and the refraction along
      ! all the paths for the current maf

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.get_path_spsfunc_ngrid' )
      Call get_path_spsfunc_ngrid ( fwdModelIn, fwdModelExtra, &
        &  fwdModelConf%molecules, ifm%ndx_path(:,maf), no_tan_hts, &
        &  ifm%z_path(:,maf), ifm%t_path(:,maf), ifm%phi_path(:,maf), n_path, &
        &  spsfunc_path, Ier )
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.get_path_spsfunc_ngrid' )
      if ( ier /= 0 ) goto 99

      !??? Choose better value for phi_tan later
      phi_tan = Deg2Rad * temp%template%phi(1,mafTInstance)

      ! Compute the ptg_angles (chi) for Antenna convolution, also the
      ! derivatives of chi w.r.t to T and other parameters
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.get_chi_angles' )
      call get_chi_angles ( ifm%ndx_path(:,maf), n_path, &
        &  fwdModelConf%tangentGrid%surfs, &
        &  ifm%tan_hts(:,mafTInstance),ifm%tan_temp(:,mafTInstance),&
        &  phi_tan,ifm%elvar(maf)%Roc,&
        &  0.001*scGeocAlt%values(1,1),  &
        &  elevOffset%values(1,1), &
        &  ifm%tan_dh_dt(:,mafTInstance,:), no_tan_hts, &
        &  temp%template%noSurfs, &
        &  temp%template%surfs(:,1), &
        &  fwdModelConf%SurfaceTangentIndex, &
        &  center_angle, ptg_angles, dx_dt, d2x_dxdt, ier )
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.get_chi_angles' )
      if ( ier /= 0 ) goto 99

      ! Compute the refraction correction scaling matrix for this mmaf:
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.refraction_correction' )
      call refraction_correction ( no_tan_hts, ifm%tan_hts(:,mafTInstance), &
        &  ifm%h_path(:,maf), n_path, ifm%ndx_path(:,maf),      &
        &  ifm%E_rad(mafTInstance), ref_corr )
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.refraction_correction' )

      Radiances = 0.0

      ! If we're not doing frequency averaging, instead outputting radiances
      ! corresponding to delta function responses, we can set up the frequency
      ! information here.  In the more common case where we are doing the
      ! averaging, the frequency grid varies from pointing to pointing, and is
      ! allocated inside the pointing loop.

      ! First we have a mini loop over pointings to work out an upper limit
      ! for the number of frequencies we're going to be dealing with
      if ( fwdModelConf%do_freq_avg ) then
        maxNoFreqs = size(PointingGrids(whichPointingGrid)%OneGrid(grids(1))%Frequencies)
        do ptg_i = 2, no_tan_hts - 1
          maxNoFreqs = max ( maxNoFreqs, size(PointingGrids(whichPointingGrid) &
            & %OneGrid(grids(ptg_i))%Frequencies) )
        end do
      else
        maxNoFreqs = noFreqs
      end if

      ! Now allocate arrays this size
      if ( fwdModelConf%temp_der ) then
        allocate ( k_temp_frq%values( maxNoFreqs, temp%template%noSurfs, &
          & windowStart:windowFinish), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
          & MLSMSG_Allocate//'k_temp_frq' )
        k_temp_frq%values = 0.0_r8
      end if

      call allocate_test ( radV,maxNoFreqs, 'radV', ModuleName )

      do specie = 1, noSpecies
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_vmr, molecule=fwdModelConf%molecules(specie) )

        ! Allocate intermediate space for vmr derivatives
        if ( fwdModelConf%moleculeDerivatives(specie) ) then
          allocate ( k_atmos_frq(specie)%values(maxNoFreqs,f%template%noSurfs,&
            & windowStart:windowFinish), stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//'k_atmos_frq' )
        end if

      end do ! End loop over species

      ! Now we can go ahead and loop over pointings
      ! ------------------------------ Begin loop over pointings --------
      do ptg_i = 1, no_tan_hts - 1
        if ( toggle(emit) .and. levels(emit) > 1 ) then
          call trace_begin ( 'ForwardModel.Pointing' )
          call output ( 'Ptg = ' ); call output ( ptg_i, advance='yes' )
        end if
        k = ptg_i
        h_tan = ifm%tan_hts(k,mafTInstance)
        lmax = ubound(ifm%eta_phi(ptg_i,maf)%values,2)

        ! Compute the beta's along the path, for this tanget hight and this mmaf:

        no_ele = ifm%ndx_path(ptg_i,maf)%total_number_of_elements

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.
        if ( FwdModelConf%do_freq_avg ) then
          frequencies => PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
          noFreqs = size(frequencies)
        end if ! If not, we dealt with this outside the loop

        call get_beta_path ( frequencies, my_Catalog, no_ele, &
          &                  ifm%z_path(ptg_i,maf), ifm%t_path(ptg_i,maf), &
          &                  beta_path, 0.001*losVel%values(1,maf), &
          &                  fwdModelConf%frqGap,             &
          &                  fwdModelConf%temp_der,           &
          &                  fwdModelConf%spect_der, Ier)
        if ( ier /= 0 ) goto 99

        ! Define the dh_dt_path for this pointing and this MAF:

        ! Need to allocate this even if no derivatives as we pass it

        call allocate_test ( dh_dt_path, no_ele, phiWindow, &
          & temp%template%noSurfs, "dh_dt_path", ModuleName )

        if ( fwdModelConf%temp_der ) then
          allocate ( dum(no_ele), stat=ier )
          if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & MLSMSG_Allocate // 'dum' )
          do j = 1, temp%template%noSurfs
            do i = 1, phiWindow
              m = min(lmax,i+windowStart-1)
              call Lintrp ( ifm%z_glgrid, ifm%z_path(ptg_i,maf)%values,&
                &           ifm%dh_dt_glgrid(:,m,j), dum, ifm%gl_count,&
                &           no_ele )
              dh_dt_path(:,i,j) = dum(:) * ifm%eta_phi(ptg_i,maf)%values(:,m)
            end do
          end do
          deallocate ( dum, stat=ier )
          if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
            & MLSMSG_DeAllocate // 'dum' )
        end if

        ! ------------------------------- Begin loop over frequencies ------
        do frq_i = 1, noFreqs

          Frq = frequencies(frq_i)
          if ( toggle(emit) .and. levels(emit) > 2 ) then
            call trace_begin ( 'ForwardModel.Frequencies' )
            call output ( 'Frq = ' ); call output ( frq_i, advance='yes' )
          end if

          Call Rad_Tran ( ifm%elvar(maf), Frq, &
            & fwdModelConf%integrationGrid%noSurfs, h_tan, &
            & noSpecies, ifm%ndx_path(k,maf), ifm%z_path(k,maf), &
            & ifm%h_path(k,maf), ifm%t_path(k,maf), ifm%phi_path(k,maf), &
            & ifm%dHdz_path(k,maf), earthRefl%values(1,1), beta_path(:,frq_i), &
            & spsfunc_path(:,k), ref_corr(:,k), spaceRadiance%values(1,1), &
            & brkpt, no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier )
          if ( ier /= 0 ) goto 99

          RadV(frq_i) = Rad

          ! Now, Compute the radiance derivatives:

!??? Do we need to do this if there's no Jacobian or no derivatives requested ???
          Call Rad_Tran_WD ( FwdModelConf, FwdModelExtra, FwdModelIn, &
            &  ifm%elvar(maf), frq_i, Frq, noSpecies, ifm%z_path(k,maf), &
            &  ifm%h_path(k,maf), ifm%t_path(k,maf), ifm%phi_path(k,maf), &
            &  ifm%dHdz_path(k,maf), beta_path(:,frq_i), spsfunc_path(:, &
            &  k), temp%template%surfs(:,1), temp%template%noSurfs, &
            &  ref_corr(:,k), temp%template%noInstances, &
            &  temp%template%phi(1,:)*Deg2Rad, dh_dt_path, k_temp_frq, &
            &  k_atmos_frq, brkpt, no_ele, mid, ilo, ihi, t_script, tau, &
            &  max_zeta_dim, max_phi_dim, ier )
          if ( ier /= 0 ) goto 99

          if ( toggle(emit) .and. levels(emit) > 2 ) &
            & call trace_end ( 'ForwardModel.Frequencies' )
        end do                          ! Frequency loop

        ! ----------------------------- End loop over frequencies ----

        ! Here we either frequency average to get the unconvolved radiances, or
        ! we just store what we have as we're using delta funciton channels

        if ( toggle(emit) .and. levels(emit) > 1 ) &
          & call trace_begin ( 'ForwardModel.FrequencyAvg' )
        if ( fwdModelConf%do_freq_avg ) then
          do i = 1, noUsedChannels
            sigInd = usedSignals(i)
            ch = usedChannels(i)
            if ( toggle(emit) .and. levels(emit) > 2 ) then
              call output ( 'Channel = ' )
              call output ( i )
              call output ( ' ( ' )
              call output ( sigInd )
              call output ( ':' )
              call output ( ch )
              call output ( ' )', advance='yes' )
            end if
            centerFreq = firstSignal%lo + &
              & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
            shapeInd = MatchSignal ( filterShapes%signal, &
              & fwdModelConf%signals(sigInd), sideband = thisSideband )
            if ( shapeInd == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
              & "No matching channel shape information" )
            if ( toggle(emit) .and. levels(emit) > 1 ) then
              call output ( 'Using filter shape:' )
              call output ( shapeInd, advance='yes' )
            endif
                          
            call Freq_Avg ( frequencies, &
              & centerFreq+thisSideband * &
              & FilterShapes(shapeInd)%FilterGrid(ch,:), &
              & FilterShapes(shapeInd)%FilterShape(ch,:), RadV, noFreqs,  &
              & Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), Radiances(ptg_i,i) )
          end do
        else
          Radiances(ptg_i,:) = RadV(1:noFreqs)
        end if

        ! Frequency Average the temperature derivatives with the appropriate
        ! filter shapes
!??? Do we need to do this if there's no Jacobian ???
        if ( fwdModelConf%temp_der ) then
          if ( fwdModelConf%do_freq_avg ) then
            do i = 1, noUsedChannels
              sigInd = usedSignals(i)
              ch = usedChannels(i)
              centerFreq = firstSignal%lo + &
                & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), &
                & sideband = thisSideband, channel=ch )
              do instance = lbound(k_temp_frq%values,3), &
                & ubound(k_temp_frq%values,3)
                do surface = 1, temp%template%noSurfs
                  ToAvg => k_temp_frq%values(1:noFreqs,surface,instance)
                  call Freq_Avg ( frequencies,                        &
                    & centerFreq+thisSideband* &
                    & FilterShapes(shapeInd)%FilterGrid(ch,:), &
                    & FilterShapes(shapeInd)%FilterShape(ch,:),&
                    & real(ToAvg,r8), noFreqs, &
                    & Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), r )
                  k_temp(i,ptg_i,surface,instance) = r
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! Channel loop
          else
            do i = 1, noUsedChannels
              k_temp(i,ptg_i,1:temp%template%noSurfs,:) = &
                &  k_temp_frq%values(i,1:temp%template%noSurfs,:)
            end do
          end if
        end if

        ! Frequency Average the atmospheric derivatives with the appropriate
        ! filter shapes
!??? Do we need to do this if there's no Jacobian ???
        do specie = 1, noSpecies
          if ( fwdModelConf%moleculeDerivatives(specie) ) then
            f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
              &  quantityType=l_vmr, molecule=fwdModelConf%molecules(specie))
            if ( fwdModelConf%do_freq_avg ) then
              do i = 1, noUsedChannels
              sigInd = usedSignals(i)
              ch = usedChannels(i)
              centerFreq = firstSignal%lo + &
                & thisSideband * fwdModelConf%signals(sigInd)%centerFrequency
              shapeInd = MatchSignal ( filterShapes%signal, &
                & fwdModelConf%signals(sigInd), &
                & sideband = thisSideband, channel=ch )
              do instance = lbound(k_atmos_frq(specie)%values,3),&
                  & ubound(k_atmos_frq(specie)%values,3)
                  do surface = 1, f%template%noSurfs
                    ToAvg => k_atmos_frq(specie)%values(1:noFreqs,surface,instance)
                    call Freq_Avg ( frequencies,                      &
                      & centerFreq+thisSideband * &
                      & FilterShapes(shapeInd)%FilterGrid(ch,:), &
                      & FilterShapes(shapeInd)%FilterShape(ch,:), &
                      & real(ToAvg,r8),  &
                      & noFreqs, Size(FilterShapes(shapeInd)%FilterGrid(ch,:)), r )
                      k_atmos(i,ptg_i,surface,instance,specie) = r
                  end do                ! Surface loop
                end do                  ! Instance loop
              end do                    ! Channel loop
            else                        ! Else not frequency averaging
              surface = f%template%noSurfs
              do i = 1, noUsedChannels
                k_atmos(i,ptg_i,1:surface,:,specie) = &
                  &  k_atmos_frq(specie)%values(i,1:surface,:)
              end do
            end if                      ! Frequency averaging or not
          end if                        ! Want derivatives for this
        end do                          ! Loop over species

        if ( toggle(emit) .and. levels(emit) > 1 ) &
          & call trace_end ( 'ForwardModel.FrequencyAvg' )

        call deallocate_test ( dh_dt_path, 'dh_dt_path', ModuleName )

        if ( toggle(emit) .and. levels(emit) > 1 ) &
          & call trace_end ( 'ForwardModel.Pointing' )

      end do                            ! Pointing Loop
      ! ---------------------------------- End of Pointing Loop ---------------

      ! Complete the radiance's last location; also complete k_temp last
      ! location.

      do i = 1, noUsedChannels
        ch = usedChannels(i)
        Radiances(no_tan_hts,i) = Radiances(no_tan_hts-1,i)
        if ( FwdModelConf%temp_der ) then
          n = temp%template%noSurfs
          k_temp(i,no_tan_hts,1:n,:) = k_temp(i,no_tan_hts-1,1:n,:)
        end if
        if ( FwdModelConf%atmos_der ) then
          do m = 1, noSpecies
            f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_vmr, molecule=fwdModelConf%molecules(m),&
              & foundInFirst=foundInFirst )
            if ( foundInFirst ) then
              n = f%template%noSurfs
              k_atmos(i,no_tan_hts,1:n,:,m)= k_atmos(i,no_tan_hts-1,1:n,:,m)
            end if
          end do
        end if
      end do

      !  Here comes the Convolution code
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_begin ( 'ForwardModel.Convolution' )
      call allocate_test ( superset, size(antennaPatterns), &
        & 'superset', ModuleName )
      do i = 1, noUsedChannels

        ch = usedChannels(i)
        sigInd = usedSignals(i)
        thisRadiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
          & signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=firstSignal%sideband )

        if ( sidebandStart /= sidebandStop ) then ! We're folding
          thisRatio = sidebandRatio%values(ch,1)
          if ( thisSideband == 1 ) thisRatio = 1.0 - thisRatio
        else ! Otherwise, want just unfolded signal
          thisRatio = 1.0
        end if

        if ( FwdModelConf%do_conv ) then

          do j = 1, size(antennaPatterns)
            superset(j) = AreSignalsSuperset ( antennaPatterns(j)%signals, &
              & fwdModelConf%signals, sideband=thisSideband, channel=ch )
          end do
          if ( all( superset < 0 ) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "No matching antenna patterns." )
          
          maxSuperset = maxval ( superset )
          where ( superset < 0 )
            superset = maxSuperset + 1
          end where
          whichPatternAsArray = minloc ( superset )
          whichPattern = whichPatternAsArray(1)
          if ( toggle(emit) .and. levels(emit) > 1 ) then
            call output ( 'Using antenna pattern: ' )
            call output ( whichPattern, advance='yes' )
          end if

          call convolve_all ( fwdModelConf, fwdModelIn, maf, ch, &
            &     windowStart, windowFinish, mafTInstance-windowStart+1, &
            &     temp, ptan, thisRadiance, &
            &     FwdModelConf%tangentGrid%surfs, ptg_angles, &
            &     ifm%tan_temp(:,maf), dx_dt, d2x_dxdt, si, center_angle, &
            &     Radiances(:, i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
            &     thisRatio, Jacobian, fmStat%rows,  &
            &     antennaPatterns(whichPattern), ier )
          !??? Need to choose some index other than 1 for AntennaPatterns ???
          if ( ier /= 0 ) goto 99
        else
          call no_conv_at_all ( fwdModelConf, fwdModelIn, maf, ch,  &
            &     windowStart, windowFinish, temp, ptan, thisRadiance, &
            &     FwdModelConf%tangentGrid%surfs, &
            &     Radiances(:,i), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
            &     thisRatio, Jacobian, fmStat%rows )

        end if

      end do                            ! Channel loop
      call deallocate_test ( superset, 'superset', ModuleName )
      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.Convolution' )

      if ( toggle(emit) .and. levels(emit) > 0 ) &
        & call trace_end ( 'ForwardModel.sideband' )
    end do
    ! ---------------------------- End of loop over sideband ------------------

    if ( maf == noMAFs ) fmStat%finished = .true.
    ! ------------------------------ End of Major Frame Specific stuff --------

    if ( associated(beta_path) ) then
      do i = 1, size(beta_path,1)
        do j = 1, size(beta_path,2)
          deallocate ( beta_path(i,j)%values, beta_path(i,j)%t_power, &
            & beta_path(i,j)%dbeta_dw, beta_path(i,j)%dbeta_dn, &
            & beta_path(i,j)%dbeta_dnu, STAT=k )
        end do
      end do
    end if

    deallocate ( beta_path, STAT=k )

    if ( FwdModelConf%temp_der ) call deallocate_test &
      & ( k_temp_frq%values, "k_temp_frq%values", ModuleName )
    if ( FwdModelConf%atmos_der ) then
      do j = 1, noSpecies
        call deallocate_test ( k_atmos_frq(j)%values, "k_atmos_frq(j)%values", &
          & ModuleName )
      end do
    end if

    if ( index(switches,'rad') /= 0 ) then
      ! *** DEBUG Print
      if ( FwdModelConf%do_conv ) then
        print *,'Convolution: ON'
      else
        print *,'Convolution: OFF'
      end if

      if ( FwdModelConf%do_freq_avg ) then
        print *,'Frequency Averaging: ON'
      else
        print *,'Frequency Averaging: OFF'
        print '(A,f12.4,a)', ' (All computations done at Frq =',Frequencies(1),')'
      end if
      print *

      do i = 1, noUsedChannels
        ch = usedChannels(i)
        print 903, ch, char(92), ptan%template%noSurfs
903     format ( 'ch', i2.2, '_pfa_rad', a1, i3.3 )
        print 905, ( firstRadiance%values(ch+(k-1)*firstRadiance%template%noChans,maf),&
          & k = 1, ptan%template%noSurfs )
905     format ( 4(2x, 1pg15.8) )
      end do
    end if

    if ( .not. fwdModelConf%do_freq_avg ) call deallocate_test ( &
      & frequencies, "frequencies", ModuleName )

    call Deallocate_test ( usedChannels, 'usedChannels', ModuleName )

99  continue

    do j = 1, No_tan_hts
      deallocate ( n_path(j)%values, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_Deallocate//'n_path%values' )
      do i = 1, noSpecies
        deallocate ( spsfunc_path(i,j)%values, stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_Deallocate//'spsfunc_path%values' )
      end do
    end do

    deallocate ( n_path, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'n_path' )

    deallocate ( spsfunc_path, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'spsfunc_path' )

    call deallocate_test ( dx_dt, 'dx_dt', ModuleName )
    call deallocate_test ( d2x_dxdt, 'd2x_dxdt', ModuleName )

    call deallocate_test ( t_script, 't_srcipt', ModuleName )
    call deallocate_test ( ref_corr, 'ref_corr', ModuleName )
    call deallocate_test ( tau, 'tau', ModuleName )
    call deallocate_test ( ptg_angles, 'ptg_angles', ModuleName )

    deallocate ( k_atmos_frq, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'k_atmos_frq' )
    deallocate ( k_temp, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'k_temp' )
    deallocate ( k_atmos, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'k_atmos' )
    deallocate ( My_Catalog, stat=status)
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'My_catalog' )

    call deallocate_test ( radiances, 'Radiances', ModuleName )
    call deallocate_test ( i_star_all, 'i_star_all', ModuleName )
    call deallocate_test ( radv, 'rad_v', ModuleName )
    call deallocate_test ( grids, 'grids',  ModuleName )

!    if ( i > -22) Stop      ! DEBUG, Zvi

    if ( toggle(emit) ) call trace_end ( 'ForwardModel' )

    Return

  end subroutine FullForwardModel

  ! =====     Private Procedures     =====================================
  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( Code, where, FieldIndex )
    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: where      ! Where in the tree did the error occur?
    integer, intent(in), optional :: FieldIndex ! f_...

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref ( where ) )
    call output ( ' ForwardModelSetup complained: ' )
    select case ( code )
    case ( AllocateError )
      call output ( 'allocation error.', advance='yes' )
    case ( BadMolecule )
      call output ( 'asked for derivatives for an unlisted molecule.', &
        & advance='yes' )
    case ( DefineMoleculesFirst )
      call output ( 'molecule must be defined before moleules derivatives.', &
        & advance='yes')
    case ( DefineSignalsFirst )
      call output ( 'signals must be defined before channels.',advance='yes')
    case ( IncompleteFullFwm )
      call output ('incomplete full foward model specification',advance='yes' )
    case ( IncompleteLinearFwm )
      call output ( 'incomplete linear foward model specification', &
        & advance='yes' )
    case ( IrrelevantFwmParameter )
      call output ( 'irrelevant parameter for this forward model type', &
        & advance='yes' )
    case ( TangentNotSubset )
      call output ('non subsurface tangent grid not a subset of integration&
        & grid', advance='yes' )
    case ( PhiWindowMustBeOdd )
      call output ( 'phiWindow is not odd', advance='yes' )
    case ( FrqGapNotFrq )
      call output ( 'frqGap does not have dimensions of frequency', advance='yes' )
    end select
  end subroutine AnnounceError

end module ForwardModelInterface

! $Log$
! Revision 2.133  2001/05/17 00:49:54  livesey
! Interim version.  Slight problem somewhere with convolution for some bands.
!
! Revision 2.132  2001/05/16 23:03:59  livesey
! New version, now gets correct antenna pattern too.
!
! Revision 2.131  2001/05/16 01:20:08  livesey
! Interim version.  Does pointing grids and channel shapes right, not
! yet antenna patterns.  Also can do forward model with multiple signals
!
! Revision 2.130  2001/05/15 03:46:31  zvi
! Adding derivative flag to beta calculations
!
! Revision 2.129  2001/05/14 23:18:26  livesey
! Added frqGap parameter
!
! Revision 2.128  2001/05/11 22:18:56  livesey
! Changed first dimension of ifm%tan_dh_dt from nlvl to no_tan_hts.
! Also tidied up allocates of ifm stuff.
!
! Revision 2.127  2001/05/03 23:34:44  livesey
! More stuff to support scan model
!
! Revision 2.126  2001/05/03 20:31:34  vsnyder
! Thought I needed to add a nullify, ended up with only cosmetic changes
!
! Revision 2.125  2001/05/02 20:31:43  livesey
! Removed frequency from config
!
! Revision 2.124  2001/05/02 20:28:35  livesey
! Removed some dead variables.
!
! Revision 2.123  2001/04/28 17:47:42  livesey
! Passes row flags to convolve and no convolve
!
! Revision 2.122  2001/04/28 01:44:09  vsnyder
! Make debugging output conditional on toggle(emit) etc.
!
! Revision 2.121  2001/04/26 22:53:37  zvi
! Fixing some phiwindow bug
!
! Revision 2.120  2001/04/26 20:02:09  livesey
! Made l2pc database a saved array in L2PC_m
!
! Revision 2.119  2001/04/26 19:47:53  livesey
! Renamed main routine to full forward model
!
! Revision 2.118  2001/04/26 02:50:47  vsnyder
! Cosmetic changes
!
! Revision 2.117  2001/04/26 02:49:52  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic (again)
!
! Revision 2.116  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.115  2001/04/26 00:06:58  livesey
! Added l2pc reading to global setup
!
! Revision 2.114  2001/04/25 00:50:33  livesey
! Changed maxPath to maxVert, better name
!
! Revision 2.113  2001/04/25 00:09:48  vsnyder
! Use 'levels(emit)' to control output detail
!
! Revision 2.112  2001/04/24 23:11:04  vsnyder
! Use 'emit' toggle (set by -f option or @E) to control output
!
! Revision 2.111  2001/04/24 19:44:48  vsnyder
! Add 'toggle(gen)' stuff in ForwardModel
!
! Revision 2.110  2001/04/24 00:03:36  livesey
! Bug fix
!
! Revision 2.109  2001/04/23 22:22:50  livesey
! Working version.  Can now have no_phi_t /= noMAFs
!
! Revision 2.108  2001/04/23 22:09:54  zvi
! Re-Introducing no_phi_t etc.
!
! Revision 2.107  2001/04/23 21:56:06  livesey
! Pass closest instances to comp_path_entities
!
! Revision 2.106  2001/04/23 21:42:46  zvi
! Introducing no_phi_t etc.
!
! Revision 2.105  2001/04/21 01:21:29  livesey
! Fixed memory leak properly!
!
! Revision 2.104  2001/04/20 23:34:54  livesey
! Fixed sideband code
!
! Revision 2.103  2001/04/20 23:08:55  livesey
! Cleaned up confusion in multi-channel cases
!
! Revision 2.102  2001/04/20 02:59:40  livesey
! Whoops typo!
!
! Revision 2.101  2001/04/20 02:55:19  livesey
! Now writes back derivatives in Jacobian!!! Not folded yet though,
! but that will be easy!
!
! Revision 2.100  2001/04/19 23:55:04  livesey
! Interim version, new convolve etc. but no derivatives
!
! Revision 2.99  2001/04/19 22:24:09  livesey
! More moving window stuff sorted out, now uses FindClosestInstance
!
! Revision 2.98  2001/04/19 22:09:32  livesey
! New calling sequence for comp_path_entities
!
! Revision 2.97  2001/04/19 20:59:50  livesey
! Reordered a loop
!
! Revision 2.96  2001/04/19 20:42:34  livesey
! Minor bug in calls to freq_avg fixed.  Going to change agains soon
! anyway
!
! Revision 2.95  2001/04/19 20:31:14  livesey
! Removed stuff that destroyed ifm (now upto calling code, e.g. in
! SidsModule) to do it.  Also added sideband folding loop, works
! for radiances, not derivatives yet.
!
! Revision 2.94  2001/04/19 08:13:06  zvi
! Some more leaks..
!
! Revision 2.92  2001/04/17 09:16:12  zvi
! Taking care of a whole buch of deallocation statements ..
!
! Revision 2.91  2001/04/17 01:01:34  vsnyder
! Add ??? Deallocated ??? comments
!
! Revision 2.90  2001/04/13 23:29:36  livesey
! Sorted out selection of appropriate pointing frequency grid.
!
! Revision 2.89  2001/04/13 21:40:22  vsnyder
! Replace pointing-grid stuff by STOP -- Nathaniel will fix it.
!
! Revision 2.88  2001/04/12 21:41:59  livesey
! Interim version.
!
! Revision 2.87  2001/04/12 17:48:31  livesey
! Moved maf increment to calling code, left finished flag here though.
!
! Revision 2.86  2001/04/12 16:55:08  livesey
! Fixed arguments to comp_path_entities
!
! Revision 2.85  2001/04/12 01:50:02  vsnyder
! Explicitly nullify instead of =>NULL()
!
! Revision 2.84  2001/04/11 02:09:46  vsnyder
! Handle 'Parse_Signal' error
!
! Revision 2.83  2001/04/11 01:18:37  vsnyder
! Check channel number range
!
! Revision 2.82  2001/04/11 00:50:06  livesey
! Another interim version, the `moving window' is better implemented
!
! Revision 2.81  2001/04/10 23:51:18  livesey
! Another working version.  More temporary arrays now alloctable/pointer
!
! Revision 2.80  2001/04/10 23:15:55  livesey
! Reverse communication seems to be working. Needs a bit more tidying up though.
!
! Revision 2.79  2001/04/10 22:04:16  livesey
! Intermediate version, slight problem with signals.
!
! Revision 2.78  2001/04/10 18:51:02  vsnyder
! Finish removing sideband stuff
!
! Revision 2.77  2001/04/10 10:15:48  zvi
! fixing bug conneced with convolve
!
! Revision 2.76  2001/04/10 02:46:16  livesey
! Working version, no more FMI/TFMI
!
! Revision 2.75  2001/04/10 02:24:55  livesey
! Stable derivative code, still not sure about vmr.  About to remove FMI, TFMI!!  :-)
!
! Revision 2.74  2001/04/10 01:16:10  livesey
! Another interim version.
!
! Revision 2.73  2001/04/09 22:21:41  livesey
! An interim version, derivatives get right numbers.
!
! Revision 2.72  2001/04/09 21:05:40  vsnyder
! Remove unneeded explicit conversion to double
!
! Revision 2.71  2001/04/09 20:51:03  zvi
! Debugging Derivatives version
!
! Revision 2.70  2001/04/07 23:59:32  zvi
! Ellimination the second dimension from ptg_angles (not MAF dependant)
!
! Revision 2.69  2001/04/07 23:49:54  zvi
! Code modified to do spsfunc & refraction inside the MAF loop
!
! Revision 2.68  2001/04/07 01:50:48  vsnyder
! Move some of VGrid to lib/VGridsDatabase.  Move FwdModelConf_T and
! some related stuff to fwdmdl/FwdModelConf.
!
! Revision 2.67  2001/04/07 01:38:22  livesey
! Another interim working version
!
! Revision 2.66  2001/04/06 21:53:40  vsnyder
! Move duplicate-field checking to init_tables
!
! Revision 2.65  2001/04/05 23:02:31  zvi
! Implementing Anntena, FilterShape & Spectroscopy l2cf inputs instead of FMI
!
! Revision 2.64  2001/04/05 22:53:20  vsnyder
! Use AntennaPatterns_m
!
! Revision 2.63  2001/04/01 00:08:52  zvi
! *** empty log message ***
!
! Revision 2.62  2001/03/31 01:49:45  zvi
! *** empty log message ***
!
! Revision 2.61  2001/03/30 20:55:25  zvi
! Remove the need for COMMON BLOCK..(ELLIPSE)
!
! Revision 2.60  2001/03/30 03:05:49  vsnyder
! Add 'antennaPatterns' field to 'forwardModelGlobal'
!
! Revision 2.59  2001/03/30 02:45:23  livesey
! Numbers agree again, was velocity.
!
! Revision 2.58  2001/03/30 01:45:08  livesey
! Some changes and debug stuff, still no agreement.
!
! Revision 2.57  2001/03/30 00:36:57  livesey
! Interim version, doesn't quite get the same numbers as Zvi, but we
! think we know why.
!
! Revision 2.56  2001/03/30 00:07:36  livesey
! Removed more FMC/TFMI stuff
!
! Revision 2.55  2001/03/29 23:56:49  livesey
! Added phi Window
!
! Revision 2.54  2001/03/29 23:42:55  vsnyder
! Add 'filterShapes' field to forwardModelGlobal
!
! Revision 2.53  2001/03/29 22:07:16  livesey
! Added phiWindow
!
! Revision 2.52  2001/03/29 12:11:16  zvi
! Fixing Bug seeting surface index erroniously
!
! Revision 2.51  2001/03/29 02:37:30  livesey
! Modified convolution to use new grids
!
! Revision 2.50  2001/03/29 01:21:25  zvi
! Interim version
!
! Revision 2.49  2001/03/29 00:53:54  livesey
! Modified error message.
!
! Revision 2.48  2001/03/29 00:33:21  livesey
! Added some error checking for tangentGrid
!
! Revision 2.47  2001/03/28 23:51:15  zvi
! Tanget below surface are in Zeta units
!
! Revision 2.46  2001/03/28 22:41:10  livesey
! Got rid of a print statement that was annoying zvi
!
! Revision 2.45  2001/03/28 22:00:10  livesey
! Interim version, now uses more allocatables etc.
!
! Revision 2.44  2001/03/28 21:22:22  vsnyder
! Use Deg2Rad from Units
!
! Revision 2.43  2001/03/28 01:31:39  livesey
! Got rid of a dump statement
!
! Revision 2.42  2001/03/28 01:31:16  livesey
! Got convolution working.
!
! Revision 2.41  2001/03/28 00:39:15  zvi
! Fixing up the convolution call
!
! Revision 2.40  2001/03/27 00:21:18  livesey
! Got frequency averaging working.
!
! Revision 2.39  2001/03/26 21:13:02  livesey
! Stableish version, frequency averaging still highly suspect.
!
! Revision 2.38  2001/03/26 18:01:20  zvi
! New code to deal with dh_dt_path being computed on the fly
!
! Revision 2.37  2001/03/25 00:50:31  livesey
! Interim version, bug with frequency averaging
!
! Revision 2.36  2001/03/24 00:33:38  livesey
! Bug fix (Typo)
!
! Revision 2.35  2001/03/24 00:32:56  livesey
! Modified use of FMI%f_grid_filter
!
! Revision 2.34  2001/03/23 23:55:54  livesey
! Another interim version.  Frequency averaging doesn't crash but
! produces bad numbers. Note that the handling of k_... when passed to
! the convolution (or noconvolution) routine is highly nefarious.
!
! Revision 2.33  2001/03/23 19:00:14  livesey
! Interim version, tidied some stuff up, still gets same numbers as Zvi
!
! Revision 2.32  2001/03/22 01:01:12  livesey
! Interim version, no rights radiances out to a vector.
!
! Revision 2.31  2001/03/21 02:14:01  livesey
! Interim version mr_f in, but not quite working yet.
!
! Revision 2.30  2001/03/21 01:10:09  livesey
! Interim version before dealing with mr_f
!
! Revision 2.29  2001/03/21 01:07:45  livesey
! Before moving away from mr_f
!
! Revision 2.28  2001/03/20 23:25:54  livesey
! Copied changes from Zvi
!
! Revision 2.27  2001/03/20 11:03:43  zvi
! Fixing code, increase dim. etc.
!
! Revision 2.26  2001/03/20 02:28:58  livesey
! Interim version, gets same numbers as Zvi
!
! Revision 2.25  2001/03/19 17:10:35  livesey
! Added more checks etc.
!
! Revision 2.24  2001/03/18 00:55:50  livesey
! Interim version
!
! Revision 2.23  2001/03/17 21:08:01  livesey
! Added forward model type stuff to FwdModelConf_T and parser thereof
!
! Revision 2.22  2001/03/17 03:24:23  vsnyder
! Work on forwardModelGlobalSetup
!
! Revision 2.21  2001/03/17 01:05:46  livesey
! OK, I've sorted it out, but problems may remain in forwardmodelglobalsetup
!
! Revision 2.20  2001/03/17 00:58:22  livesey
! Fixed bug in previous commit, have had to comment out line 139 to
! let it compile.
!
! Revision 2.19  2001/03/17 00:50:57  livesey
! New fwdModelConf stuff and merge from Van
!
! Revision 2.18  2001/03/16 21:05:22  vsnyder
! Move dumping of pointing grid database to PointingGrid_m
!
! Revision 2.17  2001/03/15 12:18:37  zvi
! Adding the velocity effect on line center frequency
!
! Revision 2.16  2001/03/13 00:43:12  zvi
! *** empty log message ***
!
! Revision 2.15  2001/03/13 00:23:41  zvi
! Correction to no_tan_hts for hydrostatic_model repeating calls
!
! Revision 2.14  2001/03/09 02:46:15  vsnyder
! Make sure "io" has a value
!
! Revision 2.13  2001/03/09 02:27:20  zvi
! *** empty log message ***
!
! Revision 2.12  2001/03/09 01:49:31  zvi
! *** empty log message ***
!
! Revision 2.11  2001/03/09 01:08:07  zvi
! *** empty log message ***
!
! Revision 2.10  2001/03/09 00:54:00  zvi
! *** empty log message ***
!
! Revision 2.9  2001/03/08 21:59:52  vsnyder
! Mostly just cosmetic rearranging
!
! Revision 2.8  2001/03/08 20:11:19  zvi
! *** empty log message ***
!
! Revision 2.7  2001/03/08 19:22:12  zvi
! New ForwardModelInterface with Zvi's code in it ..
!
! Revision 2.6  2001/03/08 03:23:45  vsnyder
! More stuff to work with L2_Load
!
! Revision 2.5  2001/03/08 00:42:09  vsnyder
! Add temporary stuff to use with L2_Load
!
! Revision 2.4  2001/03/07 23:59:52  vsnyder
! Add stuff for SIDS.
!
! Revision 2.3  2001/02/21 00:07:57  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.2  2001/02/08 00:56:11  vsnyder
! Periodic commit.  Still needs a lot of work.
!
! Revision 2.1  2001/02/07 00:52:27  vsnyder
! Initial commit
