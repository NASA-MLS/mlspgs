! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
  !=============================================================================

  ! Set up the forward model.  Interface from the retrieve step to the
  ! forward model.

  !??? Do we want a forward model database ???

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use AntennaPatterns_m, only: Close_Antenna_Patterns_File, &
    & Open_Antenna_Patterns_File, Read_Antenna_Patterns_File
  use Declaration_Table, only: NUM_VALUE, RANGE
  use Dump_0, only: DUMP
  use Expr_M, only: EXPR
  use FilterShapes_m, only: Close_Filter_Shapes_File, &
    & Open_Filter_Shapes_File, Read_Filter_Shapes_File
  ! We're going to use lots of things from init_tables_module, so let's sort
  ! them into some sort of order
  ! First admin stuff
  use Init_Tables_Module, only: FIELD_FIRST, FIELD_INDICES, FIELD_LAST, &
    & LIT_INDICES, SPEC_INDICES
  ! Now fields
  use Init_Tables_Module, only: F_ANTENNAPATTERNS, F_ATMOS_DER, F_CHANNELS, &
    & F_DO_CONV, F_DO_FREQ_AVG, F_FILTERSHAPES, F_FREQUENCY, &
    & F_INTEGRATIONGRID, F_MOLECULES, F_MOLECULEDERIVATIVES, F_PHIWINDOW, &
    & F_POINTINGGRIDS,F_SIGNALS, F_SPECT_DER, F_TANGENTGRID, F_TEMP_DER, F_TYPE
  ! Now literals
  use Init_Tables_Module, only: L_CHANNEL, L_EARTHREFL, L_ELEVOFFSET, L_FULL, L_FOLDED, &
    & L_LINEAR, L_LOSVEL, L_LOWER, L_NONE, L_ORBITINCLINE, L_PTAN, L_RADIANCE,&
    & L_REFGPH, L_SCAN, L_SCGEOCALT, L_SPACERADIANCE, L_TEMPERATURE, L_UPPER,&
    & L_VMR
  ! Now temporary stuff we hope not to have to use in the end
  use Init_Tables_Module, only: F_ZVI
  ! That's it for Init_Tables_Module
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
    & MLSMSG_Error
  use MLSNumerics, only: Hunt
  use MLSSignals_m, only: GetSignal, MaxSigLen, Signal_T, GetSignalName
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use PointingGrid_m, only: Close_Pointing_Grid_File, &
    & Open_Pointing_Grid_File, Read_Pointing_Grid_File, PointingGrids
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named
  use Units, only: DegToRad => Deg2Rad
  use VectorsModule, only: GetVectorQuantityByType, ValidateVectorQuantity, &
    & Vector_T, VectorValue_T
  use VGrid, only: VGRID_T, DUMP

  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO

  implicit none
  private
  public :: AddForwardModelConfigToDatabase, ConstructForwardModelConfig, &
    Dump, DestroyFWMConfigDatabase, ForwardModel, ForwardModelGlobalSetup, &
    ForwardModelConfig_T, ForwardModelSignalInfo_T

  interface DUMP
    module procedure DUMP_FORWARDMODELCONFIGS
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  type ForwardModelSignalInfo_T
    integer :: signal                   ! The signal we're considering
    logical, dimension(:), pointer :: channelIncluded=>null() ! Which channels to use
  end type ForwardModelSignalInfo_T

  type ForwardModelConfig_T
    integer :: fwmType        ! l_linear, l_full or l_scan
    logical :: Atmos_Der      ! Do atmospheric derivatives
    logical :: Do_Conv        ! Do convolution
    logical :: Do_Freq_Avg    ! Do Frequency averaging
    integer, dimension(:), pointer :: molecules=>NULL() ! Which molecules to consider
    logical, dimension(:), pointer :: moleculeDerivatives=>NULL() ! Want jacobians
    real(r8) :: The_Freq      ! Frequency to use if .not. do_freq_avg
    type (ForwardModelSignalInfo_T), dimension(:), pointer :: siginfo=>NULL()
    logical :: Spect_Der      ! Do spectroscopy derivatives
    logical :: Temp_Der       ! Do temperature derivatives
    type(vGrid_T), pointer :: integrationGrid ! Zeta grid for integration
    type(vGrid_T), pointer :: tangentGrid     ! Zeta grid for integration
    integer :: surfaceTangentIndex  ! Index in Tangentgrid of Earth's surface
    integer :: phiWindow            ! Window size for examining stuff
  end type ForwardModelConfig_T

  ! Error codes

  integer, parameter :: AllocateError        = 1
  integer, parameter :: BadMolecule          = AllocateError + 1
  integer, parameter :: DefineSignalsFirst   = BadMolecule + 1
  integer, parameter :: DefineMoleculesFirst = DefineSignalsFirst + 1
  integer, parameter :: DuplicateField       = DefineMoleculesFirst + 1
  integer, parameter :: IncompleteFullFwm    = DuplicateField + 1
  integer, parameter :: IncompleteLinearFwm  = IncompleteFullFwm + 1
  integer, parameter :: IrrelevantFwmParameter = IncompleteLinearFwm + 1
  integer, parameter :: TangentNotSubset     =  IrrelevantFwmParameter + 1
  integer, parameter :: PhiWindowMustBeOdd   = TangentNotSubset + 1

  integer :: Error            ! Error level -- 0 = OK

contains

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
  type (ForwardModelConfig_T) function ConstructForwardModelConfig &
    & ( ROOT, VGRIDS ) result(info)
    ! Process the forwardModel specification to produce ForwardModelConfig to add
    ! to the database

    integer, intent(in) :: ROOT         ! of the forwardModel specification.
    !                                     Indexes either a "named" or
    !                                     "spec_args" vertex. Local variables
    type (vGrid_T), dimension(:), target :: vGrids ! vGrid database

    type (Signal_T) :: thisSignal ! A signal
    type (ForwardModelSignalInfo_T), pointer :: thisFWMSignal=>NULL() ! A signal
    integer :: Field                    ! Field index -- f_something
    logical :: Got(field_first:field_last)   ! "Got this field already"
    integer :: I                        ! Subscript and loop inductor.
    integer :: J                        ! Subscript and loop inductor.
    integer :: Key                      ! Indexes the spec_args vertex.
    integer :: Name                     ! sub_rosa of label of specification,
    ! if any, else zero.
    integer :: NoChannelsSpecs          ! Number of channel specs we've had
    integer :: Son                      ! Some subtree of root.
    integer :: STATUS                   ! From allocates etc.
    integer :: THISMOLECULE             ! Tree index.
    integer :: type                     ! Type of value returned by EXPR
    integer :: TANGENT                  ! Loop counter
    integer :: Units(2)                 ! Units of value returned by EXPR
    real (r8) :: Value(2)               ! Value returned by EXPR

    ! Error message codes

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
    info%atmos_der = .false.
    info%do_conv = .false.
    info%do_freq_avg = .false.
    info%the_freq = 0.0
    info%spect_der = .false.
    info%temp_der = .false.

    noChannelsSpecs=0

    ! "Key" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.

    got = .false.
    do i = 2, nsons(key)
      son = subtree(i,key)
      field = get_field_id(son)
      if ( got(field) .and. (field /= f_channels) ) &
        &  call AnnounceError( DuplicateField, key, field)
      got(field) = .true.
      select case ( field )
      case (f_type)
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
      case ( f_channels )
        if ( .not. associated( info%siginfo ) ) then
          call AnnounceError(DefineSignalsFirst, key)
        else
          noChannelsSpecs=noChannelsSpecs+1
          thisFWMSignal => info%sigInfo(noChannelsSpecs)
          ! Now default to none included
          thisFWMSignal%channelIncluded = .false.
          do j = 2, nsons(son)
            call expr ( subtree(j,son), units, value, type )
            select case (type)
            case (num_value)
              thisFWMSignal%channelIncluded(int(value(1))) = .true.
            case (range)
              thisFWMSignal%channelIncluded(int(value(1)):int(value(2))) =&
                & .true.
            case default
              ! Shouldn't get here if parser worked
            end select
          end do
        end if
      case ( f_molecules )
        allocate ( info%molecules(nsons(son)-1), stat = status)
        if (status /= 0) call AnnounceError( AllocateError, key )
        allocate ( info%moleculeDerivatives (nsons(son)-1), stat = status)
        if (status /= 0) call AnnounceError( AllocateError, key )
        info%moleculeDerivatives = .false.
        do j = 1, nsons(son)-1
          info%molecules(j) = decoration( subtree( j+1, son ) )
        end do                          ! End loop over listed signals
      case ( f_moleculeDerivatives )
        if (.not. associated(info%molecules)) then
          call AnnounceError( DefineMoleculesFirst, key)
        else
          do j = 1, nsons(son)-1
            thisMolecule = decoration( subtree( j+1, son ) )
            if (.not. any(info%molecules == thisMolecule ) ) &
              & call AnnounceError( BadMolecule, key )
            where (info%molecules == thisMolecule)
              info%moleculeDerivatives = .true.
            end where
          end do                          ! End loop over listed signals
        end if
      case ( f_signals )
        allocate ( info%sigInfo (nsons(son)-1), stat = status)
        if (status /= 0) call AnnounceError( AllocateError, key )
        do j = 1, nsons(son)-1
          info%sigInfo(j)%signal = &
            & decoration(decoration( subtree(j+1, son )))
          ! Now allocate the channels information
          thisSignal = GetSignal( info%sigInfo(j)%signal )
          allocate ( info%sigInfo(j)%channelIncluded( &
            & lbound(thisSignal%frequencies,1):ubound(thisSignal%frequencies,1)),&
            & stat=status)
          if (status /= 0) then
            call AnnounceError ( AllocateError, key )
            exit
          end if
          ! Default to include this channel
          info%sigInfo(j)%channelIncluded = .true.
        end do                          ! End loop over listed signals
      case ( f_phiWindow )
        call expr( subtree(2,son), units, value, type )
        info%phiWindow = nint( value(1) )
        if (mod(info%phiWindow,2) /= 1) call AnnounceError (phiWindowMustBeOdd, key)
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
    select case (info%fwmType)
    case (l_full)
      if (.not. all(got( (/f_molecules, f_signals/) ))) &
        & call AnnounceError (IncompleteFullFwm, root)

      ! Now identify the Earth's surface in the tangent grid
      call Hunt(info%tangentGrid%surfs, info%integrationGrid%surfs(1), &
        &  info%surfaceTangentIndex)

      ! Ensure that points in tangentGrid at and above the surface are a subset
      ! of integration grid
      do tangent = info%surfaceTangentIndex, info%tangentGrid%noSurfs
        if (.not. any ( abs(info%tangentGrid%surfs(tangent) - &
          & info%integrationGrid%surfs) < 1e-4) ) &
          & call AnnounceError (TangentNotSubset, root)
      end do

      ! Check parameters needed only for linear/scan are not included
      !????
    case (l_scan)
      ! Add 1d/2d method later probably !??? NJL
      if (any(got( (/f_atmos_der, f_channels, f_do_conv, &
        & f_do_freq_avg, f_frequency, f_molecules, f_moleculeDerivatives, &
        & f_signals, f_spect_der, f_temp_der /) ))) &
        & call AnnounceError (IrrelevantFwmParameter, root)
    case (l_linear)
      if (.not. all(got( (/f_molecules, f_signals/) ))) & ! Maybe others later
        & call AnnounceError (IncompleteLinearFwm, root)
      if (any(got( (/f_atmos_der, f_do_conv, f_do_freq_avg, f_frequency /) ))) &
        & call AnnounceError (IrrelevantFwmParameter, root)
    end select

    if (error /= 0) call MLSMessage(MLSMSG_Error, ModuleName, 'An error occured')
    if ( toggle(gen) ) call trace_end ( "ConstructForwardModelConfig" )

  end function ConstructForwardModelConfig

  ! -----------------------------------------------  ForwardModel  -----
  ! subroutine ForwardModel ( ForwardModelConfig, FwdModelExtra, FwdModelIn, &
  !   &                       Jacobian, RowBlock, FwdModelOut )
  subroutine ForwardModel ( ForwardModelConfig, FwdModelExtra, FwdModelIn, &
    &                       Jacobian, RowBlock, FwdModelOut, FMI, TFMI)

    use GL6P, only: NG
    use MLSCommon, only: I4, R4, R8
    use L2_TEST_STRUCTURES_M
    use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
    use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
    use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi, &
      MNM => max_no_mmaf
    use ELLIPSE_M, only: ELLIPSE
    use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
    use REFRACTION_M, only: REFRACTION_CORRECTION
    use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
      PATH_DERIVATIVE, PATH_VECTOR_2D
    use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
    use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
    use GET_BETA_PATH_M, only: GET_BETA_PATH
    use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
    use RAD_TRAN_M, only: RAD_TRAN
    ! use RAD_TRAN_WD_M, only: RAD_TRAN_WD
    use FREQ_AVG_M, only: FREQ_AVG
    use CONVOLVE_ALL_M, only: CONVOLVE_ALL
    use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
    use D_LINTRP_M, only: LINTRP
    use D_HUNT_M, only: hunt_zvi => HUNT

    !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    ! Dummy arguments --------------------------------------------------------

    type(forwardModelConfig_T), intent(in) :: forwardModelConfig ! From ForwardModelSetup
    type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn ! ???
    type(matrix_T), intent(inout), optional :: Jacobian
    integer, intent(in), optional :: RowBlock          ! With which block of
    ! rows of F and Jacobian are we computing? All of them if absent.
    type(vector_T), intent(inout), optional :: FwdModelOut  ! Radiances, etc.

    !??? Begin temporary stuff to start up the forward model
    !   type(fwd_mdl_info), dimension(:), pointer, optional :: FMI
    type(fwd_mdl_info),                        optional :: FMI
    !   type(temporary_fwd_mdl_info), dimension(:), pointer, optional :: TFMI
    type(temporary_fwd_mdl_info),                        optional :: TFMI
    !??? End of temporary stuff to start up the forward model

    ! Local parameters ---------------------------------------------------------

    character, parameter :: INVALIDQUANTITY="Invalid vector quantity for "

    ! Local variables ----------------------------------------------------------

    ! First the old stuff which we hope to get rid of or redefine
    integer(i4), parameter :: NGT = (Ng+1) * N2lvl

    integer(i4) :: i, j, k, kz, ht_i, mnz, no_tan_hts, ch, Spectag, &
      m, prev_npf, ier, maf, si, ptg_i, &
      frq_i, klo, n, brkpt, no_ele, mid, ilo, ihi, &
      k_info_count, gl_count, ld

    real(r8) :: t_script(N2lvl),ref_corr(N2lvl,Nptg),tau(N2lvl), &
                tan_dh_dt(Nlvl,mnm,mxco)

    real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
    real(r8) :: dh_dt_glgrid(ngt,mnm,mxco), dhdz_glgrid(ngt,mnp)

    real(r8) :: ptg_angles(Nptg,mnm), center_angle
    real(r8) :: tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)


    !   Real(r4) :: K_TEMP(Nch,Nptg,mxco,mnp)
    !   Real(r4) :: K_ATMOS(Nch,Nptg,mxco,mnp,Nsps)
    !   Real(r4) :: K_SPECT_DW(Nch,Nptg,mxco,mnp,Nsps),  &
    !               K_SPECT_DN(Nch,Nptg,mxco,mnp,Nsps),  &
    !               K_SPECT_DNU(Nch,Nptg,mxco,mnp,Nsps)

    ! ** DEBUG, memory limitations force us to have up to 2 channels
    !           only (Replacing Nch by: 2)

    real(r4) :: K_TEMP(02,Nptg,mxco,mnp)
    real(r4) :: K_ATMOS(02,Nptg,mxco,mnp,Nsps)
    real(r4) :: K_SPECT_DW(02,Nptg,mxco,mnp,Nsps),  &
                K_SPECT_DN(02,Nptg,mxco,mnp,Nsps),  &
                K_SPECT_DNU(02,Nptg,mxco,mnp,Nsps)

    real(r8) :: I_STAR_ALL(Nch,Nptg)

    real(r4) :: K_STAR_ALL(02,20,mxco,mnp,Nptg)      ! 02 should be: Nch
    type(k_matrix_info) :: k_star_info(20)

    type(path_derivative) :: k_temp_frq, k_atmos_frq(Nsps), &
      k_spect_dw_frq(Nsps), k_spect_dn_frq(Nsps), &
      k_spect_dnu_frq(Nsps)

    type(path_beta), dimension(:,:), pointer :: beta_path => null()

    real(r8) :: Radiances(Nptg,Nch)
    real(r8) :: Zeta, Frq, h_tan, Rad, geod_lat, phi_tan, r

    Real(r8), DIMENSION(:), ALLOCATABLE :: dum

    character (LEN=01) :: CA
    character (LEN=08) :: Name
    character (LEN=16) :: Vname
    character (LEN=80) :: Line
    character (LEN=40) :: Ax, Dtm1, Dtm2

    real(r8), dimension(:), pointer :: RadV=>NULL()

    ! This is the `legit stuff' we hope will stay they are all pointers to
    ! VectorValue_T's containing vector quantities
    type (VectorValue_T), pointer :: RADIANCE=>NULL() ! Radiance quantity to be filled
    type (VectorValue_T), pointer :: TEMP=>NULL() ! Temperature quantity
    type (VectorValue_T), pointer :: PTAN=>NULL() ! PTAN quantity
    type (VectorValue_T), pointer :: ELEVOFFSET=>NULL() ! Elevation offset quantity
    type (VectorValue_T), pointer :: ORBINCLINE=>NULL() ! Orbital inclination (beta)
    type (VectorValue_T), pointer :: SPACERADIANCE=>NULL() ! Space radiance
    type (VectorValue_T), pointer :: EARTHREFL=>NULL() ! Earth reflectivity
    type (VectorValue_T), pointer :: REFGPH=>NULL() ! Reference GPH, (zRef and hRef)
    type (VectorValue_T), pointer :: LOSVEL=>NULL() ! Line of sight velocity
    type (VectorValue_T), pointer :: SCGEOCALT=>NULL() ! Geocentric spacecraft altitude
    type (VectorValue_T), pointer :: F=>NULL() ! An arbitrary species

    integer :: WHICHPOINTINGGRID        ! Index of poiting grid
    integer :: MIF                      ! Loop counter
    integer :: CHANNEL                  ! Loop counter
    integer :: MAXNOFREQS               ! Used for sizing arrays
    integer :: NOUSEDCHANNELS           ! Number of channels to output
    integer :: NOFREQS                  ! Number of frequencies for a pointing
    integer :: NOMAFS                   ! Number of major frames
    integer :: NOSPECIES                ! Number of molecules we're considering
    integer :: NAMELEN                  ! Length of string
    integer :: SPECIE                   ! Loop counter
    integer :: SURFACE                  ! Loop counter
    integer :: STATUS                   ! From allocates etc.
    integer :: INSTANCE                 ! Loop counter

    real (r8) :: CENTERFREQ             ! Of band
    real (r8) :: SENSE                  ! Multiplier (+/-1)

    integer, dimension(:), pointer :: CHANNELINDEX=>NULL() ! E.g. 1..25
    integer, dimension(:), pointer :: USEDCHANNELS=>NULL() ! Array of indices used

    real(r8), dimension(:), pointer :: GEOC_LAT=>NULL() ! Geocentric latitude of maf (radians)
    real(r8), dimension(:), pointer :: E_RAD=>NULL() ! Effective earth radius at maf (km)
    real(r4), dimension(:), pointer :: TOAVERAGE=>NULL()   ! Stuff to be passed to frq.avg.
    real(r8), dimension(:), pointer :: FREQUENCIES=>NULL() ! Frequency points
    real(r8), dimension(:,:,:), allocatable :: DH_DT_PATH  ! (pathSize, Tsurfs, Tinstance)
    real(r8), dimension(:,:), pointer :: DX_DT=>NULL() ! (Nptg, Tsurfs)
    real(r8), dimension(:,:), pointer :: D2X_DXDT=>NULL() ! (Nptg, Tsurfs)

    integer, pointer, dimension(:) :: GRIDS=>NULL()       ! Frq grid for each tan_press

    type(path_index), allocatable, dimension(:,:) :: NDX_PATH ! (Nptg,mnm)

    type(path_vector), allocatable, dimension(:,:) :: DHDZ_PATH ! (Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:) :: H_PATH    ! (Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:) :: N_PATH    ! (Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:) :: PHI_PATH  ! (Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:) :: T_PATH    ! (Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:) :: Z_PATH    ! (Nptg,mnm)

! dimensions of SPSFUNC_PATH are: (Nsps,Nptg,mnm)
    type(path_vector), allocatable, dimension(:,:,:) :: SPSFUNC_PATH

    Type(path_vector_2d), allocatable, dimension(:,:) :: ETA_PHI ! (Nptg,mnm)

    Type(ELLIPSE), allocatable, dimension(:)  :: ELVAR    ! mnm

    type(signal_t) :: Signal

    ! Executable code --------------------------------------------------------

    ! First we identify the vector quantities we're going to need.
    ! The key is to identify the signal we'll be working with first
    ! Deal with multiple signals in later versions !??? NJL
    if (size(forwardModelConfig%sigInfo) > 1) call MLSMessage ( &
      & MLSMSG_Error,ModuleName, &
      & "Can't yet have multiple signals in forward model")
    signal = GetSignal(forwardModelConfig%sigInfo(1)%signal)

    ! Now from that we identify the radiance quantity we'll be outputting
    radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal= forwardModelConfig%sigInfo(1)%signal )

    ! Identify the appropriate state vector components, save vmrs for later
    temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature)
    ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, instrumentModule=signal%instrumentModule)
    elevOffset => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_elevOffset, radiometer=signal%radiometer)
    orbIncline => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_orbitIncline)
    spaceRadiance => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_spaceRadiance)
    earthRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_earthRefl)
    refGPH => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_refGPH)
    losVel => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_losVel, instrumentModule=signal%instrumentModule)
    scGeocAlt => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scGeocAlt)

    print*,'Just some checks:',temp%template%noInstances, &
      & radiance%template%noInstances, ptan%template%noInstances
    ! We won't seek for molecules here as we can't have an array of pointers.
    ! When we do want molecule i we would do something like
    ! vmr => GetVectorQuantityBytype (fwdModelIn, fwdModelExtra, &
    !   quantityType=l_vmr, molecule=forwardModelConfig.molecules(i))

    ! Now we're going to validate the quantities we've been given, don't forget
    ! we already know what their quantityType's are as that's how we found them
    !, so we don't need to check that.
    if (.not. ValidateVectorQuantity(radiance, minorFrame=.true.,&
      & frequencyCoordinate=(/l_channel/))) call MLSMessage(MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'radiance')
    if (.not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/))) call MLSMessage(MLSMSG_Error, ModuleName,&
      & InvalidQuantity//'temperature')
    if (.not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/))) call MLSMessage(MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'ptan')
    if (.not. ValidateVectorQuantity(elevOffset, verticalCoordinate=(/l_none/), &
      & frequencyCoordinate=(/l_none/), noInstances=(/1/))) &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & InvalidQuantity//'elevOffset')
    ! There will be more to come here.

    ! Work out what signal we're after
    signal = getSignal ( forwardModelConfig%sigInfo(1)%signal )

    ! Get some dimensions which we'll use a lot
    noMAFs = radiance%template%noInstances

    ! Work out which channels are used
    call allocate_test (channelIndex, size(signal%frequencies), 'channelIndex', ModuleName)
    noUsedChannels = count ( forwardModelConfig%sigInfo(1)%channelIncluded)
    call allocate_test (usedChannels, noUsedChannels, 'channelIndex', ModuleName)
    do channel = 1, size( signal%frequencies)
      channelIndex(channel) = channel
    end do
    usedChannels = pack ( channelIndex, forwardModelConfig%sigInfo(1)%channelIncluded)
    call deallocate_test(channelIndex,'channelIndex',ModuleName)

    ! Look at the species
    noSpecies = size(forwardModelConfig%molecules)

!    Nptg = forwardModelConfig%tangentGrid%noSurfs

    ! Now we're going to create the many temporary arrays we need
    allocate (ndx_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'ndx_path')
    allocate (dhdz_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'dhdz_path')
    allocate (h_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'h_path')
    allocate (n_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'n_path')
    allocate (phi_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'phi_path')
    allocate (t_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'t_path')
    allocate (z_path(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'z_path')

    allocate (spsfunc_path(noSpecies,Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'spsfunc_path')
    allocate (eta_phi(Nptg,noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'eta_phi')

    allocate(elvar(noMAFs), STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Allocate//'elvar')

    call allocate_test(dx_dt, Nptg, temp%template%noSurfs, &
      & 'dx_dt', ModuleName)
    call allocate_test(d2x_dxdt, Nptg, temp%template%noSurfs, &
      & 'd2x_dxdt', ModuleName)

    call allocate_test(geoc_lat, noMAFs, 'geoc_lat', ModuleName)
    call allocate_test(e_rad, noMAFs, 'e_rad', ModuleName)

    ! Assert radiance%template%noInstances=temp%template%noInstances
    if (temp%template%noInstances /= noMAFs) &
      & call MLSMessage(MLSMSG_Error,ModuleName,'no temperature profiles /= no maf')
    do maf = 1, noMAFs
      phi_tan = degToRad*temp%template%phi(1,maf)
         ! ??? For the moment, change this soon.
      print*,'MAF ',maf,' phi_tan ',phi_tan
      geod_lat= degToRad*temp%template%geodLat(1,maf)
      call geoc_geod_conv(elvar(maf),orbIncline%values(1,1), &
        &  phi_tan,geod_lat, geoc_lat(maf),E_rad(maf))
    end do

    ! Compute the hydrostatic_model on the GL-Grid for all maf(s):
    ! First extend the grid below the surface.

    ! Now compute a hydrostatic grid given the temperature and refGPH
    ! information.
    call hydrostatic_model(ForwardModelConfig%SurfaceTangentIndex, &
      &  noMAFs, geoc_lat,refGPH%values(1,:)/1.0e3, &
      &  refGPH%template%surfs(1,1), &
      &  ForwardModelConfig%integrationGrid%surfs, &
      &  temp%template%surfs(:,1),temp%values,z_glgrid,h_glgrid,t_glgrid, &
      &  dhdz_glgrid,dh_dt_glgrid, &
      &  ForwardModelConfig%TangentGrid%surfs, &
      &  tan_hts,tan_temp,tan_dh_dt, &
      &  gl_count, Ier)
    if(ier /= 0) goto 99

    ! Now we have the full information about the number of tangent heights,
    ! including the subsrufaces ones.

    whichPointingGrid = signal%pointingGrid
    if ( whichPointingGrid <= 0 ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
      & "There is no pointing grid for the desired signal" )
    call allocate_test ( grids, ForwardModelConfig%TangentGrid%nosurfs, &
    "Grids", moduleName )
    call Hunt(PointingGrids(whichPointingGrid)%oneGrid%height, &
      & ForwardModelConfig%TangentGrid%surfs, grids, allowTopValue=.true.)

    no_tan_hts = ForwardModelConfig%TangentGrid%nosurfs
!    phi_window = ForwardModelConfig%phiWindow
    ! Now compute stuff along the path given this hydrostatic grid.
    call comp_path_entities(fwdModelIn, fwdModelExtra, &
      &  forwardModelConfig%molecules, &
      &  ForwardModelConfig%integrationGrid%noSurfs,temp%template%noSurfs,&
      &  gl_count,ndx_path,z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid, &
      &  tan_hts,no_tan_hts, z_path,h_path,t_path, phi_path,n_path,&
      &  dhdz_path,eta_phi,temp%template%noInstances, &
      &  temp%template%phi(1,:)*degToRad,spsfunc_path,noMAFs,elvar,Ier)
    if(ier /= 0) goto 99

    ! The first part of the forward model dealt with the chunks as a whole.
    ! This next part is more complex, and is performed within a global outer
    ! loop over major frame (maf)

    ! ---------------------------- Begin main Major Frame loop --------

    do maf = 3, 3
      !do maf = 1, noMAFs
      print*,'Doing maf:',maf

      phi_tan = degtorad*temp%template%phi(1,maf) !??? Choose better value later

      ! Compute the ptg_angles (chi) for Antenna convolution, also the
      ! derivatives of chi w.r.t to T and other parameters
      call get_chi_angles(ndx_path(:,maf),n_path(:,maf),&
        &     forwardModelConfig%tangentGrid%surfs,        &
        &     tan_hts(:,maf),tan_temp(:,maf),phi_tan,elvar(maf)%Roc,&
        &     1.0e-3*scGeocAlt%values(1,1),  &
        &     elevOffset%values(1,1), &
        &     tan_dh_dt(:,maf,:),no_tan_hts,temp%template%noSurfs,temp%template%surfs(:,1),&
        &     forwardModelConfig%SurfaceTangentIndex, &
        &     center_angle,ptg_angles(:,maf),dx_dt,d2x_dxdt,ier)
      if(ier /= 0) goto 99

      ! Compute the refraction correction scaling matrix for this mmaf:
      call refraction_correction(no_tan_hts, tan_hts(:,maf), h_path(:,maf), &
        &                n_path(:,maf), ndx_path(:,maf), E_rad(maf), ref_corr)

      prev_npf = -1
      Radiances(1:Nptg,1:Nch) = 0.0

      ! If we're not doing frequency averaging, instead outputting radiances
      ! corresponding to delta function responses, we can setup the frequency
      ! information here.  In the more common case where we are doing the
      ! averaging the frequency grid varies from pointing to pointing, and is
      ! allocated inside the pointing loop.
      if (.not. forwardModelConfig%do_freq_avg) then
        ! Think later about multiple signals case!???
        noFreqs = count ( forwardModelConfig%sigInfo(1)%channelIncluded )
        call allocate_test(frequencies,noFreqs,"frequencies",ModuleName)
        frequencies = pack ( signal%frequencies, &
          & forwardModelConfig%sigInfo(1)%channelIncluded )
        select case (signal%sideband)
        case ( l_lower )
          frequencies = signal%lo - (signal%centerFrequency+frequencies)
        case ( l_upper )
          frequencies = signal%lo + (signal%centerFrequency+frequencies)
        case ( l_folded )
          call MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Folded signal requested in forward model')
        case default
          call MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Bad value of signal%sideband')
        end select
        noFreqs = size(frequencies)
!       print*,'Doing frequencies'
!       print*,frequencies
      endif

      ! First we have a mini loop over pointings to work out an upper limit
      ! for the number of frequencies we're going to be dealing with
      maxNoFreqs = size(PointingGrids(whichPointingGrid)%OneGrid(grids(1))%Frequencies)
      do ptg_i = 2, no_tan_hts - 1
        maxNoFreqs = max ( maxNoFreqs, size(PointingGrids(whichPointingGrid) &
          & %OneGrid(grids(ptg_i))%Frequencies) )
      end do

      ! Now allocate arrays this size

      if ( forwardModelConfig%temp_der ) then
        call Allocate_Test(k_temp_frq%values, maxNoFreqs, &
          &   temp%template%noSurfs, temp%template%noInstances, &
          &   'k_temp_frq', ModuleName)
      end if

      call Allocate_Test(radV,maxNoFreqs, 'radV', ModuleName)

      do specie = 1, noSpecies
        f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_vmr, molecule=forwardModelConfig%molecules(specie))

        ! Allocate intermediate space for vmr derivatives
        if ( forwardModelConfig%moleculeDerivatives(specie) ) then
          call Allocate_Test( k_atmos_frq(specie)%values, &
            & maxNoFreqs, f%template%noSurfs, f%template%noInstances, &
            & 'k_atmos_frq(..)', ModuleName )
        end if

        ! Allocate intermediate space for spectroscopy derivatives
        if ( forwardModelConfig%spect_der ) then
          call Allocate_Test( k_spect_dw_frq(specie)%values, &
            & maxNoFreqs, f%template%noSurfs, f%template%noInstances, &
            & 'k_spect_dw_frq(..)', ModuleName )
          call Allocate_Test( k_spect_dn_frq(specie)%values, &
            & maxNoFreqs, f%template%noSurfs, f%template%noInstances, &
            & 'k_spect_dn_frq(..)', ModuleName )
          call Allocate_Test( k_spect_dnu_frq(specie)%values, &
            & maxNoFreqs, f%template%noSurfs, f%template%noInstances, &
            & 'k_spect_dnu_frq(..)', ModuleName )
        endif
      end do ! End loop over speices

      ! Now we can go ahead and loop over pointings
      ! ------------------------------ Begin loop over pointings --------
      do ptg_i = 1, no_tan_hts - 1
!       print*,'  Doing pointing:',ptg_i
        k = ptg_i
        h_tan = tan_hts(k,maf)

        ! Compute the beta's along the path, for this tanget hight and this mmaf:

        no_ele = ndx_path(ptg_i,maf)%total_number_of_elements

        ! If we're doing frequency averaging, get the frequencies we need for
        ! this pointing.
        if (ForwardModelConfig%do_freq_avg) then
!           call Allocate_Test(frequencies, &
!             & size(PointingGrids(whichPointingGrid)%&
!             &      oneGrid(grids(ptg_i))%frequencies),&
!             &      'frequencies', ModuleName)
!          Note that this deallocates the ones from the previous go round
!          print*,'Center frequency is:',PointingGrids(whichPointingGrid)%CenterFrequency
!          print*,'Offsets:'
!          call dump(PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies )
          frequencies => PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
          noFreqs = size(frequencies)
        endif ! If not, we dealt with this outside the loop

        call get_beta_path(frequencies,&
          & FMI%pfa_spectrum,no_ele, z_path(ptg_i,maf),t_path(ptg_i,maf), &
          & beta_path, 1.0e-3*losVel%values(1,maf),ier)
        if(ier /= 0) goto 99
!
!  Define the dh_dt_path for this pointing and this MAF:
!
        if ( forwardModelConfig%temp_der ) then
          ALLOCATE(dh_dt_path(no_ele,temp%template%noInstances, &
                 & temp%template%noSurfs),STAT=ier)
          if(ier == 0) ALLOCATE(dum(no_ele),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: dum or dh_dt_path, STAT =',ier
            goto 99
          else
            do j = 1,  temp%template%noSurfs
              do i = 1, temp%template%noInstances
                CALL Lintrp(z_glgrid, z_path(ptg_i,maf)%values,       &
               &            dh_dt_glgrid(:,i,j), dum, gl_count, no_ele)
                dh_dt_path(:,i,j) = dum(:) * eta_phi(ptg_i,maf)%values(:,i)
              end do
            end do
            DEALLOCATE(dum,STAT=i)
          endif
        end if
!

        ! ------------------------------- Begin loop over frequencies ------
        do frq_i = 1, noFreqs

          Frq = frequencies(frq_i)

          call Rad_Tran(elvar(maf), Frq, &
            &    forwardModelConfig%integrationGrid%noSurfs, &
            &    h_tan, noSpecies, ndx_path(k,maf),  &
            &    z_path(k,maf), h_path(k,maf), t_path(k,maf), phi_path(k,maf),&
            &    dHdz_path(k,maf), earthRefl%values(1,1), beta_path(:,frq_i),      &
            &    spsfunc_path(:,k,maf), ref_corr(:,k), spaceRadiance%values(1,1), brkpt, &
            &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
          if(ier /= 0) goto 99

          RadV(frq_i) = Rad

          ! Now, Compute the radiance derivatives:
          !         CALL Rad_Tran_WD(elvar(maf),frq_i,FMI%band,Frq, &
          !        &     ForwardModelConfig%integrationGrid%noSrfs,FMI%n_sps, &
          !        &     ForwardModelConfig%temp_der,forwardModelConfig%atmos_der, &
          !        &     ForwardModelConfig%spect_der,            &
          !        &     z_path(k,maf),h_path(k,maf),t_path(k,maf),phi_path(k,maf),   &
          !        &     dHdz_path(k,maf),TFMI%atmospheric,beta_path(:,frq_i),&
          !        &     spsfunc_path(:,k,maf),temp%template%surfs(:,1),  &
          !        &     TFMI%f_zeta_basis,TFMI%no_coeffs_f,   &
          !        &     TFMI%mr_f,TFMI%no_t,ref_corr(:,k),TFMI%no_phi_f,       &
          !        &     TFMI%f_phi_basis,temp%template%noInstances, &
          !        &     temp%template%phi(1,:)*degToRad,  &
          !        &     dh_dt_path,FMI%spect_atmos,         &
          !        &     FMI%spectroscopic,k_temp_frq,k_atmos_frq,k_spect_dw_frq,   &
          !        &     k_spect_dn_frq,k_spect_dnu_frq,TFMI%is_f_log,brkpt,       &
          !        &     no_ele,mid,ilo,ihi,t_script,tau,ier)
          !         IF(ier /= 0) goto 99

        end do                          ! Frequency loop
        ! ----------------------------- End loop over frequencies ----

        ! Here we either frequency average to get the unconvolved radiances, or
        ! we just store what we have as we're using delta funciton channels

        if ( forwardModelConfig%do_freq_avg ) then
          select case (signal%sideband)
          case (l_lower)
            centerFreq = signal%lo - signal%centerFrequency
            sense = -1.0
          case (l_upper)
            centerFreq = signal%lo + signal%centerFrequency
            sense = +1.0
          case default              ! Put error message here later
            call MLSMessage(MLSMSG_Error, ModuleName, &
              & 'Asked for folded in wrong place')
          end select

          do i = 1, noUsedChannels
            ch = usedChannels(i)
            call Freq_Avg(frequencies,centerFreq+sense*FMI%F_grid_filter(:,i),  &
              &     FMI%Filter_func(:,ch),RadV,noFreqs,FMI%no_filt_pts, &
              &     Radiances(ptg_i,ch))
          end do
        else
          Radiances(ptg_i, usedChannels) = RadV
        endif

        ! Frequency Average the temperature derivatives with the appropriate
        ! filter shapes
        if ( forwardModelConfig%temp_der) then
          if ( forwardModelConfig%do_freq_avg) then
            do i = 1, noUsedChannels
              ch = usedChannels(i)
              do instance = 1, temp%template%noInstances
                do surface = 1, temp%template%noSurfs
                  toAverage => k_temp_frq%values( &
                    & 1:noFreqs,surface,instance)
                  call Freq_Avg(frequencies, centerFreq+sense*FMI%F_grid_filter(:,i), &
                    & FMI%Filter_func(:,i), real(toAverage,r8), &
                    & noFreqs,FMI%no_filt_pts,r)
                  k_temp(ch,ptg_i,surface,instance) = r
                end do                  ! Surface loop
              end do                    ! Instance loop
            end do                      ! Channel loop
          else
            k_temp(:,ptg_i,:,:) = k_temp_frq%values
          endif
        endif

        ! Frequency Average the atmospheric derivatives with the appropriate
        ! filter shapes
        do specie = 1, noSpecies
          if ( forwardModelConfig%moleculeDerivatives(specie) ) then
            f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_vmr, molecule=forwardModelConfig%molecules(specie))

            if ( forwardModelConfig%do_freq_avg) then
              do i = 1, noUsedChannels
                ch = usedChannels(i)
                do instance = 1, f%template%noInstances
                  do surface = 1, f%template%noSurfs
                    toAverage => k_atmos_frq(specie)%values(1:noFreqs,surface,instance)
                    call Freq_Avg(frequencies,centerFreq+sense*FMI%F_grid_filter(:,i), &
                      & FMI%Filter_func(:,i), real(toAverage, r8), &
                      & noFreqs,FMI%no_filt_pts,r)
                    k_atmos(ch,ptg_i,surface,instance,specie) = r
                  end do                ! Surface loop
                end do                  ! Instance loop
              end do                    ! Channel loop
            else                        ! Else not frequency averaging
              k_atmos(:,ptg_i,:,:,specie)=  k_atmos_frq(specie)%values
            end if                      ! Frequency averaging or not
          end if                        ! Want derivatives for this
        end do                          ! Loop over species

!        if(ForwardModelConfig%spect_der) then
!           ! Frequency Average the spectroscopic derivatives with the appropriate
!           ! filter shapes
!           do i = 1, noUsedChannels
!             ch = usedChannels(i)
!             do m = 1, FMI%n_sps
!               j = FMI%spect_atmos(m)
!               if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) cycle
!               Spectag = FMI%spectroscopic(j)%Spectag
!               do
!                 if(FMI%spectroscopic(j)%Spectag /= Spectag) exit
!                 RadV(1:noFreqs) = 0.0
!                 CA = FMI%spectroscopic(j)%type
!                 do k = 1, FMI%spectroscopic(j)%no_phi_values
!                   do n = 1, FMI%spectroscopic(j)%no_zeta_values
!                     select case ( CA )
!                     case ( 'W' )
!                       RadV(1:noFreqs) = k_spect_dw_frq(m)%values(1:noFreqs,n,k)
!                     case ( 'N' )
!                       RadV(1:noFreqs) = k_spect_dn_frq(m)%values(1:noFreqs,n,k)
!                     case ( 'V' )
!                       RadV(1:noFreqs) = k_spect_dnu_frq(m)%values(1:noFreqs,n,k)
!                     end select
!                     if(ForwardModelConfig%do_freq_avg) then
!                       call Freq_Avg(frequencies,centerFrequency+sense*FMI%F_grid_filter(:,i), &
!                         &              FMI%Filter_func(:,i),&
!                         &              RadV,noFreqs,FMI%no_filt_pts,r)
!                     else
!                       r = RadV(1)
!                     endif
!                     select case ( CA )
!                     case ( 'W' )
!                       k_spect_dw(ch,ptg_i,n,k,j) = r
!                     case ( 'N' )
!                       k_spect_dn(ch,ptg_i,n,k,j) = r
!                     case ( 'V' )
!                       k_spect_dnu(ch,ptg_i,n,k,j) = r
!                     end select
!                   end do
!                 end do
!                 j = j + 1
!                 if(j > 3 * FMI%n_sps) exit
!               end do
!             end do
!           end do
!        endif

      end do                            ! Pointing Loop
      ! ---------------------------------- End of Pointing Loop ---------------

      ! Complete the radiances's last location, also  complete k_temp last
      ! location as well as k_atmos last location and k_spect_d? last location:

      do i = 1, noUsedChannels
        ch = usedChannels(i)
        Radiances(no_tan_hts,ch) = Radiances(no_tan_hts-1,ch)
        if(ForwardModelConfig%temp_der) then
          k_temp(i,no_tan_hts,1:temp%template%noSurfs,1:temp%template%noInstances) = &
            &              k_temp(i,no_tan_hts-1,1:temp%template%noSurfs,&
            &              1:temp%template%noInstances)
        endif
        if(ForwardModelConfig%atmos_der) then
          do m = 1, noSpecies
            f => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
              & quantityType=l_vmr, molecule=forwardModelConfig%molecules(m))

            if(TFMI%atmospheric(m)%der_calc(FMI%band)) then
              k = f%template%noInstances
              n = f%template%noSurfs
              k_atmos(i,no_tan_hts,1:n,1:k,m)=k_atmos(i,no_tan_hts-1,1:n,1:k,m)
            endif
          end do
        endif
        if(ForwardModelConfig%spect_der) then
          do m = 1, noSpecies
            j = FMI%spect_atmos(m)
            if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) cycle
            Spectag =  FMI%spectroscopic(j)%Spectag
            do
              if(FMI%spectroscopic(j)%Spectag /= Spectag) exit
              k = FMI%spectroscopic(j)%no_phi_values
              n = FMI%spectroscopic(j)%no_zeta_values
              k_spect_dw(i,no_tan_hts,1:n,1:k,j)=k_spect_dw(i,no_tan_hts-1,1:n,1:k,j)
              k_spect_dn(i,no_tan_hts,1:n,1:k,j)=k_spect_dn(i,no_tan_hts-1,1:n,1:k,j)
              k_spect_dnu(i,no_tan_hts,1:n,1:k,j)=k_spect_dnu(i,no_tan_hts-1,1:n,1:k,j)
              j = j + 1
              if(j > 3 * FMI%n_sps) exit
            end do
          end do
        endif
      end do

      !  Here comes the Convolution code
      do i = 1, noUsedChannels

        ch = usedChannels(i)

        if(ForwardModelConfig%do_conv) then

          ! Note I am replacing the i's in the k's with 1's (enclosed in
          ! brackets to make it clear.)  We're not wanting derivatives anyway
          ! so it shouldn't matter
          call convolve_all(ptan%values(:,maf),TFMI%atmospheric, &
            &     noSpecies, &
            &     ForwardModelConfig%temp_der,ForwardModelConfig%atmos_der,&
            &     ForwardModelConfig%spect_der,                   &
            &     ForwardModelConfig%tangentGrid%surfs,&
            &     ptg_angles(:,maf),tan_temp(:,maf), &
            &     dx_dt, d2x_dxdt,FMI%band,si,center_angle,FMI%fft_pts,   &
            &     Radiances(:,ch),k_temp((1),:,:,:),k_atmos((1),:,:,:,:), &
            &     k_spect_dw((1),:,:,:,:),k_spect_dn((1),:,:,:,:),    &
            &     k_spect_dnu((1),:,:,:,:),FMI%spect_atmos,&
            &     ForwardModelConfig%tangentGrid%noSurfs,  &
            &     k_info_count,i_star_all(i,:),k_star_all((1),:,:,:,:), &
            &     k_star_info,temp%template%noSurfs,temp%template%noInstances,&
            &     TFMI%no_phi_f, FMI%spectroscopic,temp%template%surfs(:,1), &
            &     FMI%Xlamda,FMI%Aaap, FMI%D1Aaap,FMI%D2Aaap,FMI%Ias,ier)
          if(ier /= 0) goto 99
        else
          ! Note I am replacing the i's in the k's with 1's (enclosed in
          ! brackets to make it clear.)  We're not wanting derivatives anyway
          ! so it shouldn't matter
          call no_conv_at_all(ptan%values(:,maf),noSpecies, &
            &     ForwardModelConfig%tangentGrid%surfs, &
            &     FMI%band,ForwardModelConfig%temp_der,&
            &     ForwardModelConfig%atmos_der,ForwardModelConfig%spect_der,      &
            &     Radiances(:,ch),k_temp((1),:,:,:),                    &
            &     k_atmos((1),:,:,:,:),k_spect_dw((1),:,:,:,:),       &
            &     k_spect_dn((1),:,:,:,:),k_spect_dnu((1),:,:,:,:),   &
            &     FMI%spect_atmos, ForwardModelConfig%tangentGrid%noSurfs,&
            &     k_info_count,               &
            &     i_star_all(i,:), k_star_all((1),:,:,:,:),            &
            &     k_star_info,temp%template%noSurfs,temp%template%noInstances, &
            &     TFMI%no_phi_f,temp%template%surfs(:,1),TFMI%atmospheric,    &
            &     FMI%spectroscopic)
        endif

      end do                            ! Channel loop

      do mif = 1, radiance%template%noSurfs
        do channel = 1, noUsedChannels
          radiance%values( usedChannels(channel)+ &
            & (mif-1)*radiance%template%noChans, maf ) = i_star_all( channel, mif )
        end do
      end do

    end do                              ! MAF loop
    ! ------------------------------ End of Major Frame Loop -----------

    if(ForwardModelConfig%temp_der) deallocate(k_temp_frq%values,STAT=i)
    do j = 1, noSpecies
      if(ForwardModelConfig%atmos_der) deallocate(k_atmos_frq(j)%values,STAT=i)
      if(ForwardModelConfig%spect_der) then
        deallocate(k_spect_dw_frq(j)%values,STAT=i)
        deallocate(k_spect_dn_frq(j)%values,STAT=i)
        deallocate(k_spect_dnu_frq(j)%values,STAT=i)
      endif
    end do

    ! *** DEBUG Print

    if(ForwardModelConfig%do_conv) then
      print *,'Convolution: ON'
    else
      print *,'Convolution: OFF'
    endif

    if(.not. ForwardModelConfig%do_freq_avg) then
      Frq = Frequencies(1)
      print *,'Frequency Averaging: OFF'
      write(6,'(A,f12.4,a1)') ' (All computations done at Frq =',Frq,')'
    else
      print *,'Frequency Averaging: ON'
    endif
    print *

    if (.not. forwardModelConfig%do_freq_avg) deallocate(frequencies)

    tau(1:Nptg) = 0.0
    tau(1:ptan%template%noSurfs) = dble(ptan%values(:,3))

    klo = -1
    Zeta = -1.666667
    call Hunt_zvi(Zeta,tau,ptan%template%noSurfs,klo,j)
    if(abs(Zeta-tau(j)) < abs(Zeta-tau(klo))) klo=j

    do i = 1, noUsedChannels
      ch = usedChannels(i)
      write(*,903) ch,char(92),ptan%template%noSurfs
      write(*,905) (i_star_all(i,k),k=1,ptan%template%noSurfs)
    end do
903 format('ch',i2.2,'_pfa_rad',a1,i3.3)
905 format(4(2x,1pg15.8))

    call Deallocate_test(usedChannels, 'usedChannels', ModuleName)
!   print*,'At the end radiances are:'
!   call dump(radiance%values)
!  *** Zvi DEBUG   
     if(i > -2) Stop
!  *** END Zvi DEBUG   

    if(.not. any((/ForwardModelConfig%temp_der,&
      & ForwardModelConfig%atmos_der,ForwardModelConfig%spect_der/))) goto 99

    m = -1
    ch = 1
    tau = 0.0
    do i = 1, k_info_count
      print *
      Name = ' '
      Name = k_star_info(i)%name
      if(Name == 'PTAN') cycle
      kz = k_star_info(i)%first_dim_index
      mnz = k_star_info(i)%no_zeta_basis
      ht_i = k_star_info(i)%no_phi_basis
      nameLen = len_trim(Name)
      if(Name(nameLen-1:nameLen) == '_W' .or.  &
        &   Name(nameLen-1:nameLen) == '_N' .or.  &
        &   Name(nameLen-1:nameLen) == '_V' ) then
        print *,Name
        r = sum(k_star_all(ch,kz,1:mnz,1:ht_i,klo))
        print *,'  Sum over all zeta & phi coeff:',sngl(r)
      else
        if(Name == 'TEMP') then
          write(6,913) 'dI_dT',char(92),ht_i
        else
          write(6,913) 'dI_d'//Name(1:nameLen),char(92),ht_i
        endif
        tau = 0.0
        tau(1:mnz) = k_star_info(i)%zeta_basis(1:mnz)
        call Hunt_zvi(Zeta,tau,mnz,m,j)
        if(abs(Zeta-tau(j)) < abs(Zeta-tau(m))) m=j
        print *,(k_star_all(ch,kz,m,ptan%template%noSurfs,klo),k=1,ht_i)
      endif
    end do

99  continue

913 format(a,a1,i2.2)

    deallocate (ndx_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'ndx_path')
    deallocate (dhdz_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'dhdz_path')
    deallocate (h_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'h_path')
    deallocate (n_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'n_path')
    deallocate (phi_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'phi_path')
    deallocate (t_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'t_path')
    deallocate (z_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'z_path')
    deallocate (spsfunc_path, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'spsfunc_path')
    deallocate (eta_phi, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'eta_phi')
    deallocate (elvar, STAT=status)
    if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName, &
      & MLSMSG_Deallocate//'elvar')

    call deallocate_test(dx_dt, 'dx_dt', ModuleName)
    call deallocate_test(d2x_dxdt, 'd2x_dxdt', ModuleName)
    call deallocate_test(geoc_lat, 'geoc_lat', ModuleName)
    call deallocate_test(e_rad, 'e_rad', ModuleName)

    return

  end subroutine ForwardModel

  ! ----------------------------  AddForwardModelConfigToDatabase  -----
  integer function AddForwardModelConfigToDatabase ( database, item )

    ! Add a quantity template to a database, or create the database if it
    ! doesn't yet exist

    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: database
    type (ForwardModelConfig_T), intent(in) :: item

    ! Local variables
    type (ForwardModelConfig_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddForwardModelConfigToDatabase = newSize
  end function AddForwardModelConfigToDatabase

  ! --------------------------  DestroyForwardModelConfigDatabase  -----
  subroutine DestroyFWMConfigDatabase(database)
    ! Dummy arguments
    type (ForwardModelConfig_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer :: CONFIG                   ! Loop counter
    integer :: SIGNAL                   ! Loop counter
    integer :: STATUS                   ! Flag

    if (associated(database)) then
      do config = 1, size(database)
        do signal = 1, size(database(config)%sigInfo)
          deallocate(database(config)%sigInfo(signal)%channelIncluded,stat=status)
          if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Deallocate//"database%signals%channelIncluded")
        end do
        deallocate (database(config)%sigInfo, stat=status)
        if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
          & MLSMSG_Deallocate//"database%signals")
        deallocate (database(config)%molecules, stat=status)
        if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
          & MLSMSG_Deallocate//"database%molecules")
      end do

      deallocate (database, stat=status)
      if (status /= 0) call MLSMessage(MLSMSG_Error,ModuleName,&
        & MLSMSG_Deallocate//"database")
    end if
  end subroutine DestroyFWMConfigDatabase

  ! =====     Private Procedures     =====================================
  ! ------------------------------------  DUMP_FOWARDMODELCONFIGS  -----
  subroutine Dump_ForwardModelConfigs ( database )
    type (ForwardModelConfig_T), dimension(:), pointer :: database

    ! Local variables
    integer :: I, J                          ! Loop counters
    character (len=MaxSigLen) :: SignalName  ! A line of text

    ! executable code
    if ( associated(database) ) then
      do i = 1, size(database)
        call output ( 'FowardModelConfig: ' )
        call output ( i, advance = 'yes' )
        call output ( '  Atmos_der:' )
        call output ( database(i)%atmos_der, advance='yes' )
        call output ( '  Do_conv:' )
        call output ( database(i)%do_conv, advance='yes' )
        call output ( '  Do_freq_avg:' )
        call output ( database(i)%do_freq_avg, advance='yes' )
        call output ( '  Spect_der:' )
        call output ( database(i)%spect_der, advance='yes' )
        call output ( '  Temp_der:' )
        call output ( database(i)%temp_der, advance='yes' )
        call output ( '  The_freq:' )
        call output ( database(i)%the_freq, advance='yes' )
        call output ( '  Molecules: ', advance='yes' )
        do j = 1, size(database(i)%molecules)
          call output ( '    ' )
          call display_string(lit_indices(database(i)%molecules(j)))
          if (database(i)%moleculeDerivatives(j)) then
            call output (' compute derivatives', advance='yes')
          else
            call output (' no derivatives', advance='yes')
          end if
        end do
        call output ( '  Signals:', advance='yes')
        do j = 1, size(database(i)%sigInfo)
          call output ( '    ' )
          call GetSignalName( database(i)%sigInfo(j)%signal, signalName)
          call output ( signalName//' channelIncluded:', advance='yes')
          call dump ( database(i)%sigInfo(j)%channelIncluded )
        end do
      end do
    end if
  end subroutine Dump_ForwardModelConfigs

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
    case (AllocateError)
      call output ( 'allocation error.', advance='yes')
    case (BadMolecule)
      call output ( 'asked for derivatives for an unlisted molecule.')
    case (DefineMoleculesFirst)
      call output ( 'molecule must be defined before moleules derivatives.', advance='yes')
    case (DefineSignalsFirst)
      call output ( 'signals must be defined before channels.', advance='yes')
    case (DuplicateField)
      call output ( 'duplicate field specified:' )
      call display_string(field_indices(FieldIndex), advance='yes')
    case (IncompleteFullFwm)
      call output ( 'incomplete full foward model specification' )
    case (IncompleteLinearFwm)
      call output ( 'incomplete linear foward model specification' )
    case (IrrelevantFwmParameter)
      call output ( 'irrelevant parameter for this forward model type' )
    case (TangentNotSubset)
      call output ( 'non subsurface tangent grid not a subset of integration grid' )
    case (PhiWindowMustBeOdd)
      call output ( 'phiWindow is not odd' )
    end select
  end subroutine AnnounceError

end module ForwardModelInterface

! $Log$
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
!
