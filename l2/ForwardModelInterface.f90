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
  use Init_Tables_Module, only: FIELD_FIRST, FIELD_INDICES, FIELD_LAST, &
    & LIT_INDICES, SPEC_INDICES
  ! Now fields
  use Init_Tables_Module, only: F_ANTENNAPATTERNS, F_ATMOS_DER, F_CHANNELS, &
    & F_DO_CONV, F_DO_FREQ_AVG, F_FILTERSHAPES, F_FREQUENCY, &
    & F_INTEGRATIONGRID, F_MOLECULES, F_MOLECULEDERIVATIVES, F_PHIWINDOW, &
    & F_POINTINGGRIDS, F_SIGNALS, F_SPECT_DER, F_TANGENTGRID, F_TEMP_DER, F_TYPE
  ! Now literals
  use Init_Tables_Module, only: L_CHANNEL, L_EARTHREFL, L_ELEVOFFSET, L_FULL, &
    & L_LINEAR, L_LOSVEL, L_NONE, L_ORBITINCLINE, L_PTAN, L_RADIANCE,&
    & L_REFGPH, L_SCAN, L_SCGEOCALT, L_SPACERADIANCE, L_TEMPERATURE, &
    & L_VMR
  ! That's it for Init_Tables_Module
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error
  use MLSNumerics, only: Hunt
  use MLSSignals_m, only: GetSignal, MaxSigLen, Signal_T, GetSignalName,&
    & MATCHSIGNAL, DUMP
  use Molecules, only: spec_tags
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use Parse_Signal_m, only: PARSE_SIGNAL
  use PointingGrid_m, only: Close_Pointing_Grid_File, &
    & Open_Pointing_Grid_File, Read_Pointing_Grid_File, PointingGrids
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named
  use Units, only: Deg2Rad
  use VectorsModule, only: GetVectorQuantityByType, ValidateVectorQuantity, &
    & Vector_T, VectorValue_T
  use VGridsDatabase, only: VGrid_T
  use dump_0, only: Dump

  implicit none
  private
  public :: ConstructForwardModelConfig, Dump, ForwardModel, &
    ForwardModelGlobalSetup

  interface Dump
    module procedure MyDump_ForwardModelConfigs
  end interface

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
    if ( toggle(gen) ) call trace_begin ( "ForwardModelGlobalSetup", root )

    ! "Root" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.
    ! Collect data from the fields.

    do i = 2, nsons(root)
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_antennaPatterns )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_antenna_patterns_file ( fileName, lun )
        call read_antenna_patterns_file ( lun, spec_indices )
        call close_antenna_patterns_file ( lun )
      case ( f_filterShapes )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_filter_shapes_file ( fileName, lun )
        call read_filter_shapes_file ( lun, spec_indices )
        call close_filter_shapes_file ( lun )
      case ( f_pointingGrids )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        call open_pointing_grid_file ( fileName, lun )
        call read_pointing_grid_file ( lun, spec_indices )
        call close_pointing_grid_file ( lun )
      case default
        ! Can't get here if the type checker worked
      end select
    end do

    if ( toggle(gen) ) call trace_end ( "ForwardModelGlobalSetup" )
  end subroutine ForwardModelGlobalSetup

  ! ------------------------------------------  ConstructForwardModelConfig  -----
  type (forwardModelConfig_T) function ConstructForwardModelConfig &
    & ( ROOT, VGRIDS ) result(info)
    ! Process the forwardModel specification to produce ForwardModelConfig to add
    ! to the database
    use MLSSignals_M, only: Signals

    integer, intent(in) :: ROOT         ! of the forwardModel specification.
    !                                     Indexes either a "named" or
    !                                     "spec_args" vertex. Local variables
    type (vGrid_T), dimension(:), target :: vGrids ! vGrid database

    type (signal_T) :: thisSignal ! A signal
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
    integer :: Son                      ! Some subtree of root.
    integer :: STATUS                   ! From allocates etc.
    integer :: THISMOLECULE             ! Tree index.
    integer :: type                     ! Type of value returned by EXPR
    integer :: TANGENT                  ! Loop counter
    integer :: Units(2)                 ! Units of value returned by EXPR
    integer :: WANTED                   ! Which signal do we want?
    real (r8) :: Value(2)               ! Value returned by EXPR
    character (len=80) :: SIGNALSTRING  ! E.g. R1A....
    logical, dimension(:), pointer :: channels   ! From Parse_Signal
    integer, dimension(:), pointer :: SIGNALINDS ! From Parse_Signal

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
    info%the_freq = 0.0
    info%temp_der = .false.
    info%atmos_der = .false.
    info%spect_der = .false.
    info%phiwindow = 5

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
      case ( f_frequency )
        call expr ( subtree(2,son), units, value, type )
        info%the_freq = value(1)
        !       case ( f_channels )
        !         if ( .not. associated( info%signals ) ) then
        !           call AnnounceError(DefineSignalsFirst, key)
        !         else
        !           noChannelsSpecs=noChannelsSpecs+1
        !           thisSignal => info%signals(noChannelsSpecs)
        !           ! Now default to none included
        !           thisSignal%channels = .false.
        !           do j = 2, nsons(son)
        !             call expr ( subtree(j,son), units, value, type )
        !             select case (type)
        !             case (num_value)
        !               thisFWMSignal%channelIncluded(int(value(1))) = .true.
        !             case (range)
        !               thisFWMSignal%channelIncluded(int(value(1)):int(value(2))) =&
        !                 & .true.
        !             case default
        !               ! Shouldn't get here if parser worked
        !             end select
        !           end do
        !         end if
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
          call parse_Signal ( signalString, signalInds, spec_indices, &
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
      if ( any(got( (/f_atmos_der, f_do_conv, f_do_freq_avg, f_frequency /) )) ) &
        & call AnnounceError ( IrrelevantFwmParameter, root )
    end select

    if ( error /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'An error occured' )
    if ( toggle(gen) ) call trace_end ( "ConstructForwardModelConfig" )

  end function ConstructForwardModelConfig

  ! -----------------------------------------------  ForwardModel  -----
  subroutine ForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
    &                       FwdModelOut, Ifm, FmStat, Jacobian )

    use GL6P, only: NG
    use MLSCommon, only: I4, R4, R8
    use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
    use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
    use L2PCdim, only: NSPS, Nptg, MNP => max_no_phi, &
      MNM => max_no_mmaf
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
    type(forwardModelConfig_T), intent(inout) :: forwardModelConfig
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: Ifm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Local parameters ---------------------------------------------------------

    character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

    ! Local variables ----------------------------------------------------------

    ! First the old stuff which we hope to get rid of or redefine
    integer(i4) :: i, j, k, kz, ht_i, mnz, no_tan_hts, ch, Spectag, &
      m, ier, maf, si, ptg_i, frq_i, klo, n, brkpt, no_ele, &
      mid, ilo, ihi, k_info_count, ld, &
      max_phi_dim, max_zeta_dim

    real(r8) :: I_STAR_ALL(25,Nptg)
    real(r4) :: K_STAR_ALL(25,20,mxco,mnp,Nptg)
    type(k_matrix_info) :: K_star_info(20)


    type(path_derivative) :: K_temp_frq
    type(path_derivative), allocatable, dimension(:) :: K_atmos_frq

    type(path_beta), dimension(:,:), pointer :: Beta_path

    real(r8) :: Zeta, Frq, H_tan, Rad, Geod_lat, Phi_tan, R

    Real(r8), dimension(:), allocatable :: dum

    character (len=01) :: CA
    character (len=08) :: Name
    character (len=16) :: Vname
    character (len=80) :: Line
    character (len=40) :: Ax, Dtm1, Dtm2

    real(r8), dimension(:), pointer :: RadV

    ! This is the `legit stuff' we hope will stay; they are all pointers to
    ! VectorValue_T's containing vector quantities.
    type (VectorValue_T), pointer :: RADIANCE      ! Radiance quantity to be filled
    type (VectorValue_T), pointer :: TEMP          ! Temperature quantity
    type (VectorValue_T), pointer :: PTAN          ! PTAN quantity
    type (VectorValue_T), pointer :: ELEVOFFSET    ! Elevation offset quantity
    type (VectorValue_T), pointer :: ORBINCLINE    ! Orbital inclination (beta)
    type (VectorValue_T), pointer :: SPACERADIANCE ! Space radiance
    type (VectorValue_T), pointer :: EARTHREFL     ! Earth reflectivity
    type (VectorValue_T), pointer :: REFGPH        ! Reference GPH, (zRef and hRef)
    type (VectorValue_T), pointer :: LOSVEL        ! Line of sight velocity
    type (VectorValue_T), pointer :: SCGEOCALT     ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: F             ! An arbitrary species

    integer :: WHICHPOINTINGGRID        ! Index of poiting grid
    integer :: MIF                      ! Loop counter
    integer :: CHANNEL                  ! Loop counter
    integer :: MAXNOFREQS               ! Used for sizing arrays
    integer :: MAXNOFSURFS              ! Max. no. surfaces for any molecule
    integer :: MAXPATH                  ! Number of points on longest path
    integer :: NLVL                     ! Size of tangent grid
    integer :: N2LVL                    ! Twice size of tangent grid
    integer :: NOUSEDCHANNELS           ! Number of channels to output
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOSPECIES                ! Number of molecules we're considering
    integer :: NAMELEN                  ! Length of string
    integer :: PHIWINDOW                ! Copy of forward model config%phiWindow
    integer :: SIG0                     ! Lower part of signal range
    integer :: SIG1                     ! Upper part of signal range
    integer :: SPECIE                   ! Loop counter
    integer :: SURFACE                  ! Loop counter
    integer :: STATUS                   ! From allocates etc.
    integer :: TOTALSIGNALS             ! Used when hunting for pointing grids
    integer :: INSTANCE                 ! Loop counter
    integer :: WINDOWFINISH             ! Range of window
    integer :: WINDOWSTART              ! Range of window

    real (r8) :: CENTERFREQ             ! Of band
    real (r8) :: CENTER_ANGLE           ! For angles
    real (r8) :: SENSE                  ! Multiplier (+/-1)

    integer, dimension(:), pointer :: CHANNELINDEX ! E.g. 1..25
    integer, dimension(:), pointer :: SIGNALSGRID ! Used in ptg grid hunt
    integer, dimension(:), pointer :: USEDCHANNELS ! Array of indices used

    logical, dimension(:), pointer :: ALLMATCH ! Used in pointing grid hunt
    logical, dimension(:), pointer :: THISMATCH ! USed in pointing grid hunt

    real(r4), dimension(:),     pointer :: TOAVG       ! Stuff to be passed to frq.avg.
    real(r8), dimension(:),     pointer :: FREQUENCIES ! Frequency points
    real(r8), dimension(:,:,:), pointer :: DH_DT_PATH  ! (pathSize, Tsurfs, Tinstance)
    real(r8), dimension(:,:),   pointer :: DX_DT       ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:,:),   pointer :: D2X_DXDT    ! (No_tan_hts, Tsurfs)
    real(r8), dimension(:),     pointer :: T_SCRIPT    ! (n2lvl)
    real(r8), dimension(:,:),   pointer :: REF_CORR    ! (n2lvl, no_tan_hts)
    real(r8), dimension(:),     pointer :: TAU         ! (n2lvl)
    real(r8), dimension(:),     pointer :: PTG_ANGLES  ! (no_tan_hts)

    real(r4), dimension(:,:,:,:), pointer :: K_TEMP    ! (channel,Nptg,mxco,mnp)
    real(r4), dimension(:,:,:,:,:), pointer :: K_ATMOS ! (channel,Nptg,mxco,mnp,Nsps)
    real(r8), dimension(:,:), pointer :: Radiances     ! (Nptg,25)

    integer, pointer, dimension(:) :: GRIDS            ! Frq grid for each tan_press

    logical :: FOUNDINFIRST             ! Flag to indicate derivatives

    type(Signal_T), pointer, dimension(:) :: ALLSIGNALS ! Used in ptg grid hunt

    type(path_vector), allocatable, dimension(:) :: N_PATH    ! (No_tan_hts)

    ! dimensions of SPSFUNC_PATH are: (Nsps,No_tan_hts)
    type(path_vector), allocatable, dimension(:,:) :: SPSFUNC_PATH
    type(signal_t) :: Signal
    type(catalog_T), pointer, dimension(:) :: My_Catalog

    ! Executable code --------------------------------------------------------

    ! Nullify a bunch of pointers so that Allocate_Test doesn't try to
    ! deallocate them.  We don't want them to be initialized NULL()
    ! because that makes them SAVEd.

    nullify ( radV, channelIndex, usedChannels, frequencies, dh_dt_path, &
      & dx_dt, d2x_dxdt, t_script, ref_corr, tau, ptg_angles, k_temp, &
      & k_atmos, radiances, grids, my_Catalog, beta_path, allMatch, thisMatch,&
      & signalsGrid )

    ! Identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    ! Deal with multiple signals in later versions !??? NJL
    if ( size(forwardModelConfig%signals) > 1 ) call MLSMessage ( &
      & MLSMSG_Error,ModuleName, &
      & "Can't yet have multiple signals in forward model" )
    signal = forwardModelConfig%signals(1)

    ! Now from that we identify the radiance quantity we'll be outputting
    radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal= signal%index, sideband=signal%sideband )

    ! Identify the appropriate state vector components, save vmrs for later
    temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature )
    ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=signal%instrumentModule )
    elevOffset => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_elevOffset, radiometer=signal%radiometer )
    orbIncline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_orbitIncline )
    spaceRadiance => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_spaceRadiance )
    earthRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_earthRefl )
    refGPH => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_refGPH )
    losVel => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_losVel, instrumentModule=signal%instrumentModule )
    scGeocAlt => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scGeocAlt )

    ! We won't seek for molecules here as we can't have an array of pointers.
    ! When we do want molecule i we would do something like
    ! vmr => GetVectorQuantityBytype (fwdModelIn, fwdModelExtra, &
    !   quantityType=l_vmr, molecule=forwardModelConfig.molecules(i))

    ! Now we're going to validate the quantities we've been given, don't forget
    ! we already know what their quantityType's are as that's how we found them
    !, so we don't need to check that.
    if ( .not. ValidateVectorQuantity(radiance, minorFrame=.true.,&
      & frequencyCoordinate=(/l_channel/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'radiance' )
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

    noSpecies = size(forwardModelConfig%molecules)
    !  Create a subset of the catalog composed only of those molecules to be
    !  used for this run

    maxNoFSurfs = 0
    do specie = 1, noSpecies
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=forwardModelConfig%molecules(specie) )
      maxNoFSurfs = max(maxNoFSurfs, f%template%noSurfs)
    end do

    allocate ( My_Catalog(noSpecies), stat=ier )
    if ( ier /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'my_catalog' )

    do j = 1, noSpecies
      Spectag = spec_tags(forwardModelConfig%molecules(j))
      do i = 1, Size(Catalog)
        if ( Catalog(i)%Spec_Tag == Spectag ) then
          My_Catalog(j) = Catalog(i)
          exit
        end if
      end do
    end do
!
    ! Get the max. dimension in zeta coeff. space and phi coeff. space
    ! (To be used later in rad_tran_wd, for automatic arrays asignement)
    max_phi_dim = 1
    max_zeta_dim = 1
    if ( ForwardModelConfig%temp_der ) then
      max_zeta_dim = temp%template%noSurfs
      max_phi_dim = temp%template%noInstances
    end if
!
    if ( forwardModelConfig%atmos_der ) then
      do k = 1, noSpecies
        if ( forwardModelConfig%moleculeDerivatives(k) ) then
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            &     quantityType=l_vmr, molecule=forwardModelConfig%molecules(k))
          j = f%template%noInstances
          max_phi_dim = max(max_phi_dim,j)
          j = f%template%noSurfs
          max_zeta_dim = max(max_zeta_dim,j)
        end if
      end do
    end if
!
    ! Work out what signal we're after
    signal = forwardModelConfig%signals(1)

    ! Get some dimensions which we'll use a lot
    noMAFs = radiance%template%noInstances
    no_tan_hts = ForwardModelConfig%TangentGrid%nosurfs
    maxPath = 2 * (NG+1) * size(ForwardModelConfig%integrationGrid%surfs)
    nlvl=size(ForwardModelConfig%integrationGrid%surfs)
    n2lvl=2*nlvl
    phiWindow = ForwardModelConfig%phiWindow
    print*,'Dimensions:'
    print*,'noMAFs:',noMAFs
    print*,'no_tan_hts:',no_tan_hts
    print*,'maxPath:',maxPath
    print*,'nlvl:',nlvl
    print*,'n2lvl:',n2lvl
    print*,'phiWindow:',phiWindow
    print*,'noSpecies:',noSpecies
    print*,'maxNoFSurfs:',maxNoFSurfs
!   fmStat%maf = 3                      ! DEBUG, Zvi
    print*,'MAF:',fmStat%maf

    ! Work out which channels are used
    call allocate_test ( channelIndex, size(signal%frequencies), 'channelIndex', &
      & ModuleName )
    noUsedChannels = count ( signal%channels )
    call allocate_test ( usedChannels, noUsedChannels, 'channelIndex', ModuleName )
    do channel = 1, size( signal%frequencies)
      channelIndex(channel) = channel
    end do
    usedChannels = pack ( channelIndex, signal%channels )
    call deallocate_test ( channelIndex,'channelIndex',ModuleName )

    if ( fmStat%newHydros ) then

      print*,'(re)computing hydrostatic stuff.'
      ! Now we're going to create the many temporary arrays we need
      allocate (ifm%ndx_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'ndx_path' )
      allocate (ifm%dhdz_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'dhdz_path' )
      allocate (ifm%h_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'h_path' )
      allocate (ifm%phi_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'phi_path' )
      allocate (ifm%t_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'t_path' )
      allocate (ifm%z_path(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'z_path' )

      allocate (ifm%eta_phi(No_tan_hts,noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'eta_phi' )

      allocate ( ifm%elvar(noMAFs), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'elvar' )

      call allocate_test ( ifm%geoc_lat, noMAFs, 'geoc_lat', ModuleName )
      call allocate_test ( ifm%e_rad, noMAFs, 'e_rad', ModuleName )

      call allocate_test ( ifm%h_glgrid, maxPath, noMAFs, 'h_glgrid', ModuleName )
      call allocate_test ( ifm%t_glgrid, maxPath, noMAFs, 't_glgrid', ModuleName )
      call allocate_test ( ifm%z_glgrid, maxPath/2, 'z_glgrid', ModuleName )
      call allocate_test ( ifm%dh_dt_glgrid, maxPath, noMAFs, &
        & temp%template%noSurfs,'dh_dt_glgrid', ModuleName )
      call allocate_test ( ifm%dhdz_glgrid, maxPath, noMAFs, 'dhdz_glgrid', ModuleName )
      call allocate_test ( ifm%tan_hts, &
        & size(ForwardModelConfig%tangentGrid%surfs), noMAFs, 'tan_hts', ModuleName )
      call allocate_test ( ifm%tan_temp, &
        & size(ForwardModelConfig%tangentGrid%surfs), noMAFs, 'tan_hts', ModuleName )
      call allocate_test ( ifm%tan_dh_dt, nlvl, noMAFs, &
        & temp%template%noSurfs, 'tan_dh_dt', ModuleName )

      ! Setup for hydrostatic calculation
      ! Assert radiance%template%noInstances=temp%template%noInstances

      if ( temp%template%noInstances /= noMAFs ) &
        & call MLSMessage(MLSMSG_Error,ModuleName,'no temperature profiles /= no maf')
      do maf = 1, noMAFs
        phi_tan = Deg2Rad*temp%template%phi(1,maf)
        ! ??? For the moment, change this soon.
        geod_lat= Deg2Rad*temp%template%geodLat(1,maf)
        call geoc_geod_conv ( ifm%elvar(maf), orbIncline%values(1,1), &
          &  phi_tan, geod_lat, ifm%geoc_lat(maf), ifm%E_rad(maf) )
      end do

      ! Now compute a hydrostatic grid given the temperature and refGPH
      ! information.
      call hydrostatic_model ( ForwardModelConfig%SurfaceTangentIndex, &
        &  noMAFs, ifm%geoc_lat, 0.001*refGPH%values(1,:), &
        &  refGPH%template%surfs(1,1), &
        &  ForwardModelConfig%integrationGrid%surfs, &
        &  temp%template%surfs(:,1), temp%values, &
        &  ifm%z_glgrid, ifm%h_glgrid, ifm%t_glgrid, &
        &  ifm%dhdz_glgrid, ifm%dh_dt_glgrid, &
        &  ForwardModelConfig%TangentGrid%surfs, &
        &  ifm%tan_hts, ifm%tan_temp, ifm%tan_dh_dt, &
        &  ifm%gl_count, Ier )
      if ( ier /= 0 ) goto 99

      ! Now compute stuff along the path given this hydrostatic grid.
      call comp_path_entities ( ForwardModelConfig%integrationGrid%noSurfs, &
        &  temp%template%noSurfs, ifm%gl_count, ifm%ndx_path, ifm%z_glgrid, &
        &  ifm%t_glgrid, ifm%h_glgrid, ifm%dhdz_glgrid, ifm%tan_hts,        &
        &  no_tan_hts, ifm%z_path, ifm%h_path, ifm%t_path, ifm%phi_path,    &
        &  ifm%dhdz_path, ifm%eta_phi, temp%template%noInstances,           &
        &  temp%template%phi(1,:)*Deg2Rad, noMAFs, phiWindow, ifm%elvar, Ier )
      if ( ier /= 0 ) goto 99

      fmStat%newHydros = .false.
    end if

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

    ! The first part of the forward model dealt with the chunks as a whole.
    ! This next part is more complex, and is performed within a global outer
    ! loop over major frame (maf)

    ! First we need to identify the pointing grid we're going to be using.
    ! I might as well go the whole way and deal with the multiple signals case
    ! while I'm about it.  The surrounding code can't yet support that though

    totalSignals = 0
    do i = 1, size(pointingGrids)
      totalSignals = totalSignals + size(pointingGrids(i)%signals)
    end do

    allocate (allSignals(totalSignals), STAT=status)
    call allocate_test( signalsGrid,   totalSignals, 'signalsGrid',  &
      & ModuleName)

    totalSignals = 0
    do i = 1, size(pointingGrids)
!     print*,3,i
      sig0 = totalSignals + 1
      sig1 = totalSignals + size(pointingGrids(i)%signals)
      allSignals(sig0:sig1) = pointingGrids(i)%signals
      signalsGrid(sig0:sig1) = i
      totalSignals = sig1
    end do

    call allocate_test( allMatch, totalSignals, 'allMatch', ModuleName )
    call allocate_test( thisMatch, totalSignals, 'thisMatch', ModuleName )

    allMatch = .true.
    do i = 1, size ( forwardModelConfig%signals)
      j = MatchSignal ( allSignals, forwardModelConfig%signals(i), thisMatch )
      allMatch = allMatch .and. thisMatch
    end do

    if ( count (allMatch) == 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
      & 'No matching pointing frequency grids' )

    ! For the moment take the first match, later we'll be cleverer and choose
    ! the smallest.  Or maybe we'll turn this whole section into another routine.
    do i = 1, totalSignals
      if ( allMatch(i) ) exit
    end do
    whichPointingGrid = i

    call deallocate_test( thisMatch, 'thisMatch', ModuleName )
    call deallocate_test( allMatch, 'allMatch', ModuleName )

    deallocate (allSignals, STAT=status)
    if ( status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
      & MLSMSG_DeAllocate//'allSignals')
    call deallocate_test(signalsGrid, 'signalsGrid', ModuleName)

    if ( whichPointingGrid <= 0 ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "There is no pointing grid for the desired signal" )

    ! Now we've identified the pointing grids.  Locate the tangent grid within
    ! it.
    call allocate_test ( grids, ForwardModelConfig%TangentGrid%nosurfs, &
      "Grids", ModuleName )
    call Hunt ( PointingGrids(whichPointingGrid)%oneGrid%height, &
      & ForwardModelConfig%TangentGrid%surfs, grids, allowTopValue=.true. )

    ! ---------------------------- Begin main MAF Specific stuff --------

    maf=fmStat%maf
    print*,'Doing maf:',maf

    ! Now work out what `window' we're inside.  This will need to be changed
    ! a bit in later versions to avoid the noMAFS==noTemp/f instances assertion
    windowStart = max(1,maf-phiWindow/2)
    windowFinish = min(maf+phiWindow/2, temp%template%noInstances)

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

    Call get_path_spsfunc_ngrid ( fwdModelIn, fwdModelExtra, &
      &     forwardModelConfig%molecules, ifm%ndx_path(:,maf), no_tan_hts, &
      &     ifm%z_path(:,maf), ifm%t_path(:,maf), ifm%phi_path(:,maf), n_path, &
      &     spsfunc_path, Ier )
    if ( ier /= 0 ) goto 99

    !??? Choose better value for phi_tan later
    phi_tan = Deg2Rad * temp%template%phi(1,maf)

    ! Compute the ptg_angles (chi) for Antenna convolution, also the
    ! derivatives of chi w.r.t to T and other parameters
    call get_chi_angles ( ifm%ndx_path(:,maf), n_path, &
      &     forwardModelConfig%tangentGrid%surfs,      &
      &     ifm%tan_hts(:,maf), ifm%tan_temp(:,maf), phi_tan, ifm%elvar(maf)%Roc, &
      &     0.001*scGeocAlt%values(1,1),  &
      &     elevOffset%values(1,1), &
      &     ifm%tan_dh_dt(:,maf,:), no_tan_hts, temp%template%noSurfs, &
      &     temp%template%surfs(:,1), &
      &     forwardModelConfig%SurfaceTangentIndex, &
      &     center_angle, ptg_angles, dx_dt, d2x_dxdt, ier )
    if ( ier /= 0 ) goto 99

    ! Compute the refraction correction scaling matrix for this mmaf:
    call refraction_correction(no_tan_hts, ifm%tan_hts(:,maf), &
      &  ifm%h_path(:,maf), n_path, ifm%ndx_path(:,maf),      &
      &  ifm%E_rad(maf), ref_corr)

    Radiances = 0.0

    ! If we're not doing frequency averaging, instead outputting radiances
    ! corresponding to delta function responses, we can setup the frequency
    ! information here.  In the more common case where we are doing the
    ! averaging the frequency grid varies from pointing to pointing, and is
    ! allocated inside the pointing loop.
    if ( .not. forwardModelConfig%do_freq_avg ) then
      ! Think later about multiple signals case!???
      noFreqs = count( forwardModelConfig%signals(1)%channels )
      call allocate_test ( frequencies,noFreqs, "frequencies", ModuleName )
      frequencies = pack( signal%frequencies, &
        & forwardModelConfig%signals(1)%channels )
!%%%%%%%%%%%%% DEBUG NJL
      print*,'Signal sideband:',signal%sideband
      signal%sideband=-1
!%%%%%%%%%%%%%%
      select case ( signal%sideband )
      case ( -1 )
        frequencies = signal%lo - (signal%centerFrequency+frequencies)
      case ( +1 )
        frequencies = signal%lo + (signal%centerFrequency+frequencies)
      case ( 0 )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Folded signal requested in forward model' )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Bad value of signal%sideband' )
      end select
      noFreqs = size(frequencies)
    end if

    ! First we have a mini loop over pointings to work out an upper limit
    ! for the number of frequencies we're going to be dealing with
    maxNoFreqs = size(PointingGrids(whichPointingGrid)%OneGrid(grids(1))%Frequencies)
    do ptg_i = 2, no_tan_hts - 1
      maxNoFreqs = max ( maxNoFreqs, size(PointingGrids(whichPointingGrid) &
        & %OneGrid(grids(ptg_i))%Frequencies) )
    end do

    ! Now allocate arrays this size
    if ( forwardModelConfig%temp_der ) then
      allocate ( k_temp_frq%values( maxNoFreqs, temp%template%noSurfs, &
        & windowStart:windowFinish), stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
        & MLSMSG_Allocate//'k_temp_frq' )
    end if

    call allocate_test ( radV,maxNoFreqs, 'radV', ModuleName )

    do specie = 1, noSpecies
      f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=forwardModelConfig%molecules(specie) )

      ! Allocate intermediate space for vmr derivatives
      if ( forwardModelConfig%moleculeDerivatives(specie) ) then
        allocate ( k_atmos_frq(specie)%values( maxNoFreqs, f%template%noSurfs, &
          & windowStart:windowFinish), stat=status )
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName,&
          & MLSMSG_Allocate//'k_atmos_frq' )
      end if

    end do ! End loop over speices

    ! Now we can go ahead and loop over pointings
    ! ------------------------------ Begin loop over pointings --------
    do ptg_i = 1, no_tan_hts - 1
      k = ptg_i
      h_tan = ifm%tan_hts(k,maf)

      ! Compute the beta's along the path, for this tanget hight and this mmaf:

      no_ele = ifm%ndx_path(ptg_i,maf)%total_number_of_elements

      ! If we're doing frequency averaging, get the frequencies we need for
      ! this pointing.
      if ( ForwardModelConfig%do_freq_avg ) then
        frequencies => PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
        noFreqs = size(frequencies)
      end if ! If not, we dealt with this outside the loop

      call get_beta_path ( frequencies, my_Catalog, no_ele, &
        &                  ifm%z_path(ptg_i,maf), ifm%t_path(ptg_i,maf), &
        &                  beta_path, 0.001*losVel%values(1,maf), ier )
      if ( ier /= 0 ) goto 99

      ! Define the dh_dt_path for this pointing and this MAF:

      ! Need to allocate this even if no derivatives as we pass it

      call allocate_test ( dh_dt_path, no_ele, temp%template%noInstances, &
        & temp%template%noSurfs, "dh_dt_path", ModuleName )

      if ( forwardModelConfig%temp_der ) then
        allocate ( dum(no_ele), stat=ier )
        if ( ier /= 0 ) then
          Print *,'** ALLOCATE Error: dum or dh_dt_path, stat =',ier
          goto 99
        else
          do j = 1,  temp%template%noSurfs
            do i = 1, temp%template%noInstances
              call Lintrp ( ifm%z_glgrid, ifm%z_path(ptg_i,maf)%values,       &
                &           ifm%dh_dt_glgrid(:,i,j), dum, ifm%gl_count, no_ele )
              dh_dt_path(:,i,j) = dum(:) * ifm%eta_phi(ptg_i,maf)%values(:,i)
            end do
          end do
          deallocate ( dum, stat=i )
        end if
      end if

      ! ------------------------------- Begin loop over frequencies ------
      do frq_i = 1, noFreqs

        Frq = frequencies(frq_i)

        Call Rad_Tran ( ifm%elvar(maf), Frq, &
          & forwardModelConfig%integrationGrid%noSurfs, h_tan, &
          & noSpecies, ifm%ndx_path(k,maf), ifm%z_path(k,maf), &
          & ifm%h_path(k,maf), ifm%t_path(k,maf), ifm%phi_path(k,maf), &
          & ifm%dHdz_path(k,maf), earthRefl%values(1,1),beta_path(:,frq_i), &
          & spsfunc_path(:,k), ref_corr(:,k), spaceRadiance%values(1,1), &
          & brkpt, no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier )
        if ( ier /= 0 ) goto 99

        RadV(frq_i) = Rad

        ! Now, Compute the radiance derivatives:

        if ( ForwardModelConfig%temp_der ) k_temp_frq%values = 0.0_r8
        Call Rad_Tran_WD ( ForwardModelConfig, FwdModelExtra, FwdModelIn, &
          &  ifm%elvar(maf), frq_i, Frq, noSpecies, ifm%z_path(k, maf), &
          &  ifm%h_path(k, maf), ifm%t_path(k, maf), ifm%phi_path(k, maf), &
          &  ifm%dHdz_path(k, maf), beta_path(:, frq_i), spsfunc_path(:, &
          &  k), temp%template%surfs(:, 1), temp%template%noSurfs, &
          &  ref_corr(:, k), temp%template%noInstances, &
          &  temp%template%phi(1, :)*Deg2Rad, dh_dt_path, k_temp_frq, &
          &  k_atmos_frq, brkpt, no_ele, mid, ilo, ihi, t_script, tau, &
          &  max_zeta_dim, max_phi_dim, ier )
        if ( ier /= 0 ) goto 99

      end do                          ! Frequency loop

      ! ----------------------------- End loop over frequencies ----

      ! Here we either frequency average to get the unconvolved radiances, or
      ! we just store what we have as we're using delta funciton channels

      if ( forwardModelConfig%do_freq_avg ) then
        if ( signal%sideband == 0 ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Asked for folded in wrong place' )
        else
          sense = signal%sideband
          centerFreq = signal%lo + sense * signal%centerFrequency
        end if
        do i = 1, noUsedChannels
          ch = usedChannels(i)
          call Freq_Avg ( frequencies,                           &
            &       centerFreq+sense*FilterShapes(i)%FilterGrid, &
            &       FilterShapes(i)%FilterShape, RadV, noFreqs,  &
            &       Size(FilterShapes(i)%FilterGrid), Radiances(ptg_i,ch) )
        end do
      else
        Radiances(ptg_i,usedChannels) = RadV(1:noFreqs)
      end if

      ! Frequency Average the temperature derivatives with the appropriate
      ! filter shapes
      if ( forwardModelConfig%temp_der ) then
        if ( forwardModelConfig%do_freq_avg ) then
          do i = 1, noUsedChannels
            ch = usedChannels(i)
            do instance = lbound(k_temp_frq%values,3), &
              & ubound(k_temp_frq%values,3)
              do surface = 1, temp%template%noSurfs
                ToAvg => k_temp_frq%values(1:noFreqs,surface,instance)
                call Freq_Avg ( frequencies,                        &
                  &        centerFreq+sense*FilterShapes(i)%FilterGrid, &
                  &        FilterShapes(i)%FilterShape, real(ToAvg,r8), &
                  &        noFreqs, Size(FilterShapes(i)%FilterGrid), r )
                k_temp(ch,ptg_i,surface,instance) = r
              end do                  ! Surface loop
            end do                    ! Instance loop
          end do                      ! Channel loop
        else
          do i = 1, noUsedChannels
            k_temp(i,ptg_i,1:temp%template%noSurfs,:) = &
              &  k_temp_frq%values(1,1:temp%template%noSurfs,:)
          end do
        end if
      end if

      ! Frequency Average the atmospheric derivatives with the appropriate
      ! filter shapes
      do specie = 1, noSpecies
        if ( forwardModelConfig%moleculeDerivatives(specie) ) then
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            &  quantityType=l_vmr, molecule=forwardModelConfig%molecules(specie))
          if ( forwardModelConfig%do_freq_avg ) then
            do i = 1, noUsedChannels
              ch = usedChannels(i)
              do instance = lbound(k_atmos_frq(specie)%values,3),&
                & ubound(k_atmos_frq(specie)%values,3)
                do surface = 1, f%template%noSurfs
                  ToAvg => k_atmos_frq(specie)%values(1:noFreqs,surface,instance)
                  call Freq_Avg ( frequencies,                      &
                    &          centerFreq+sense*FilterShapes(i)%FilterGrid, &
                    &          FilterShapes(i)%FilterShape, real(ToAvg,r8),  &
                    &          noFreqs, Size(FilterShapes(i)%FilterGrid), r )
                  k_atmos(ch,ptg_i,surface,instance,specie) = r
                end do                ! Surface loop
              end do                  ! Instance loop
            end do                    ! Channel loop
          else                        ! Else not frequency averaging
            surface = f%template%noSurfs
            do i = 1, noUsedChannels
              k_atmos(i,ptg_i,1:surface,:,specie) = &
                &  k_atmos_frq(specie)%values(1,1:surface,:)
            end do
          end if                      ! Frequency averaging or not
        end if                        ! Want derivatives for this
      end do                          ! Loop over species


      call deallocate_test ( dh_dt_path, 'dh_dt_path', ModuleName )

    end do                            ! Pointing Loop
    ! ---------------------------------- End of Pointing Loop ---------------

    ! Complete the radiances's last location, also complete k_temp last
    ! location.

    do i = 1, noUsedChannels
      ch = usedChannels(i)
      Radiances(no_tan_hts,ch) = Radiances(no_tan_hts-1,ch)
      if ( ForwardModelConfig%temp_der ) then
        n = temp%template%noSurfs
        k_temp(i,no_tan_hts,1:n,1:phiWindow) = &
          & k_temp(i,no_tan_hts-1,1:n,1:phiWindow)
      end if
      if ( ForwardModelConfig%atmos_der ) then
        do m = 1, noSpecies
          f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_vmr, molecule=forwardModelConfig%molecules(m),&
            & foundInFirst=foundInFirst )
          if ( foundInFirst ) then
            k = f%template%noInstances
            n = f%template%noSurfs
            k_atmos(i,no_tan_hts,1:n,1:phiWindow,m)= &
              & k_atmos(i,no_tan_hts-1,1:n,1:phiWindow,m)
          end if
        end do
      end if
    end do

    !  Here comes the Convolution code
    do i = 1, noUsedChannels

      ch = usedChannels(i)

      if ( ForwardModelConfig%do_conv ) then

        ! Note I am replacing the i's in the k's with 1's (enclosed in
        ! brackets to make it clear.)  We're not wanting derivatives anyway
        ! so it shouldn't matter
        call convolve_all ( forwardModelConfig, fwdModelIn, &
          &     ptan%values(:,maf), noSpecies, &
          &     ForwardModelConfig%tangentGrid%surfs, ptg_angles, &
          &     ifm%tan_temp(:,maf), dx_dt, d2x_dxdt,si, center_angle, &
          &     Radiances(:,ch), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
          &     no_tan_hts, k_info_count, i_star_all(i,:), &
          &     k_star_all(i,:,:,:,:), k_star_info, &
          &     temp%template%noSurfs, temp%template%noInstances, &
          &     temp%template%surfs(:,1), antennaPatterns(1), ier )
        !??? Need to choose some index other than 1 for AntennaPatterns ???
        if ( ier /= 0 ) goto 99
      else

        call no_conv_at_all ( forwardModelConfig, fwdModelIn, &
          &     ptan%values(:,maf), noSpecies,  &
          &     ForwardModelConfig%tangentGrid%surfs, &
          &     Radiances(:,ch), k_temp(i,:,:,:), k_atmos(i,:,:,:,:), &
          &     no_tan_hts, k_info_count, i_star_all(i,:), &
          &     k_star_all(i,:,:,:,:), k_star_info, &
          &     temp%template%noSurfs, temp%template%noInstances, &
          &     temp%template%surfs(:,1) )

      end if

    end do                            ! Channel loop

    do mif = 1, radiance%template%noSurfs
      do channel = 1, noUsedChannels
        k = usedChannels(channel)+(mif-1)*radiance%template%noChans
        radiance%values(k,maf) = i_star_all(channel,mif)
      end do
    end do

    ! ------------------------------ End of Major Frame Specific stuff --------

    if ( maf == noMAFs ) fmStat%finished = .true.

    if ( associated(beta_path) ) then
      do i = 1, size(beta_path,1)
        do j = 1, size(beta_path,2)
          DEALLOCATE ( beta_path(i,j)%values, beta_path(i,j)%t_power, &
            & beta_path(i,j)%dbeta_dw, beta_path(i,j)%dbeta_dn, &
            & beta_path(i,j)%dbeta_dnu, STAT=k )
        end do
      end do
    end if

    deallocate(beta_path,STAT=k)

    if(ForwardModelConfig%temp_der) call deallocate_test (k_temp_frq%values,&
      & "k_temp_frq%values", ModuleName )
    if ( ForwardModelConfig%atmos_der ) then
      do j = 1, noSpecies
        call deallocate_test (k_atmos_frq(j)%values, "k_atmos_frq(j)%values",&
          & ModuleName )
      end do
    end if

    ! *** DEBUG Print

    if ( ForwardModelConfig%do_conv ) then
      print *,'Convolution: ON'
    else
      print *,'Convolution: OFF'
    end if

    if ( .not. ForwardModelConfig%do_freq_avg ) then
      Frq = Frequencies(1)
      print *,'Frequency Averaging: OFF'
      write(6,'(A,f12.4,a1)') ' (All computations done at Frq =',Frq,')'
    else
      print *,'Frequency Averaging: ON'
    end if
    print *

    if ( .not. forwardModelConfig%do_freq_avg) call deallocate_test ( &
      & frequencies, "frequencies", ModuleName )

    do i = 1, noUsedChannels
      ch = usedChannels(i)
      write(*,903) ch, char(92), ptan%template%noSurfs
      write(*,905) ( i_star_all(i,k), k = 1, ptan%template%noSurfs )
    end do
903 format('ch',i2.2,'_pfa_rad',a1,i3.3)
905 format(4(2x,1pg15.8))

    call Deallocate_test ( usedChannels, 'usedChannels', ModuleName )

    ! ** DEBUG, Zvi
    !   if ( i > -22) Stop
    ! ** END DEBUG

    if ( .not. any((/ForwardModelConfig%temp_der,&
      & ForwardModelConfig%atmos_der,ForwardModelConfig%spect_der/)) ) goto 99

    m = ptan%template%noSurfs
    tau(1:m) = ptan%values(1:m,3)
    tau(m+1:) = 0.0

    klo = -1
    Zeta = -1.666667
    call Hunt_zvi ( Zeta, tau, m, klo, j )
    if ( abs(Zeta-tau(j)) < abs(Zeta-tau(klo))) klo=j

    m = -1
    ch = 1
    tau = 0.0
    do i = 1, k_info_count
      print *
      Name = ' '
      Name = k_star_info(i)%name
      if ( Name == 'PTAN') cycle
      kz = k_star_info(i)%first_dim_index
      mnz = k_star_info(i)%no_zeta_basis
      ht_i = k_star_info(i)%no_phi_basis
      nameLen = len_trim(Name)
      if ( Name(nameLen-1:nameLen) == '_W' .or.  &
        &   Name(nameLen-1:nameLen) == '_N' .or.  &
        &   Name(nameLen-1:nameLen) == '_V' ) then
        print *,Name
        r = sum(k_star_all(ch,kz,1:mnz,1:ht_i,klo))
        print *,'  Sum over all zeta & phi coeff:',sngl(r)
      else
        if ( Name == 'TEMP' ) then
          write(6,913) 'dI_dT',char(92),ht_i
        else
          write(6,913) 'dI_d'//Name(1:nameLen),char(92),ht_i
        end if
        tau = 0.0
        tau(1:mnz) = k_star_info(i)%zeta_basis(1:mnz)
        call Hunt_zvi ( Zeta, tau, mnz, m, j )
        if ( abs(Zeta-tau(j)) < abs(Zeta-tau(m))) m=j
        print *,(k_star_all(ch,kz,m,j,klo),j=1,ht_i)
      end if
    end do

99  continue

913 format(a,a1,i2.2)

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

    if ( fmStat%Finished ) then

      deallocate (ifm%ndx_path, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'ndx_path' )
!
      do j = 1, noMAFs
        do i = 1, No_tan_hts
          deallocate (ifm%z_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'z_path%values' )
          deallocate (ifm%h_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'h_path%values' )
          deallocate (ifm%t_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'t_path%values' )
          deallocate (ifm%phi_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'phi_path%values' )
          deallocate (ifm%dhdz_path(i,j)%values, stat=status )
          if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'dhdz_path%values' )
          deallocate (ifm%eta_phi(i,j)%values, stat=status )
          if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//'eta_phi%values' )
        end do
      end do
!
      deallocate (ifm%z_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'z_path' )
      deallocate (ifm%h_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'h_path' )
      deallocate (ifm%t_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'t_path' )
      deallocate (ifm%phi_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'phi_path' )
      deallocate (ifm%dhdz_path, stat=status )
      if( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'dhdz_path' )
      deallocate (ifm%eta_phi, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'eta_phi' )

      deallocate ( ifm%elvar, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error,ModuleName, &
        & MLSMSG_Allocate//'elvar' )

      call deallocate_test ( ifm%geoc_lat, 'geoc_lat', ModuleName )
      call deallocate_test ( ifm%e_rad, 'e_rad', ModuleName )

      call deallocate_test ( ifm%h_glgrid, 'h_glgrid', ModuleName )
      call deallocate_test ( ifm%t_glgrid, 't_glgrid', ModuleName )
      call deallocate_test ( ifm%z_glgrid, 'z_glgrid', ModuleName )
      call deallocate_test ( ifm%dh_dt_glgrid, 'dh_dt_glgrid', ModuleName )
      call deallocate_test ( ifm%dhdz_glgrid, 'dhdz_glgrid', ModuleName )
      call deallocate_test ( ifm%tan_hts,'tan_hts', ModuleName )
      call deallocate_test ( ifm%tan_dh_dt, 'tan_dh_dt', ModuleName )

    end if

    ! ** DEBUG, Zvi
    !    if ( i > -22) Stop
    ! ** END DEBUG

    Return

  end subroutine ForwardModel

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
      call output ( 'signals must be defined before channels.', advance='yes')
    case ( IncompleteFullFwm )
      call output ( 'incomplete full foward model specification', advance='yes' )
    case ( IncompleteLinearFwm )
      call output ( 'incomplete linear foward model specification', &
        & advance='yes' )
    case ( IrrelevantFwmParameter )
      call output ( 'irrelevant parameter for this forward model type', &
        & advance='yes' )
    case ( TangentNotSubset )
      call output ( 'non subsurface tangent grid not a subset of integration grid', &
        & advance='yes' )
    case ( PhiWindowMustBeOdd )
      call output ( 'phiWindow is not odd', advance='yes' )
    end select
  end subroutine AnnounceError

  ! ---------------------------------  MyDump_ForwardModelConfigs  -----
  subroutine MyDump_ForwardModelConfigs ( ForwardModelConfigs )
    type(forwardModelConfig_T), pointer, dimension(:) :: ForwardModelConfigs
    call dump ( forwardModelConfigs, lit_indices )
  end subroutine MyDump_ForwardModelConfigs

end module ForwardModelInterface

! $Log$
! Revision 2.93  2001/04/19 06:46:35  zvi
! Fixing Memory leaks ..
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
! Move some of VGrid to lib/VGridsDatabase.  Move ForwardModelConfig_T and
! some related stuff to fwdmdl/ForwardModelConfig.
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
! Added forward model type stuff to ForwardModelConfig_T and parser thereof
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
! New forwardModelConfig stuff and merge from Van
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
