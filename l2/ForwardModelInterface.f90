! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
!=============================================================================

! Set up the forward model.  Interface from the retrieve step to the
! forward model.

!??? Do we want a forward model database ???

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Declaration_Table, only: NUM_VALUE, RANGE
  use Dump_0, only: DUMP
  use Expr_M, only: EXPR
  ! We're going to use lots of things from init_tables_module, so lets sort
  ! them into some sort of order
  ! First admin stuff
  use Init_Tables_Module, only: FIELD_FIRST, FIELD_INDICES, FIELD_LAST, &
    & LIT_INDICES, SPEC_INDICES
  ! Now fields
  use Init_Tables_Module, only: F_ATMOS_DER, F_CHANNELS, F_DO_CONV, F_DO_FREQ_AVG, &
    F_FREQUENCY, F_MOLECULES, F_MOLECULEDERIVATIVES, F_POINTINGGRIDS, F_SIGNALS, &
    F_SPECT_DER, F_TEMP_DER, F_TYPE
  ! Now literals
  use Init_Tables_Module, only: L_CHANNEL, L_EARTHREFL, L_ELEVOFFSET, L_FULL, L_LINEAR,&
    & L_LOSVEL, L_NONE, L_ORBITINCLINE, L_PTAN, L_RADIANCE, L_REFGPH, L_SCAN,&
    & L_SCGEOCALT, L_SPACERADIANCE, L_TEMPERATURE
  ! Now temporary stuff we hope not to have to use in the end
  use Init_Tables_Module, only: F_ZVI
  ! That's it for Init_Tables_Module
  use Lexer_Core, only: Print_Source
  use MatrixModule_1, only: Matrix_Database_T, Matrix_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, MLSMSG_Error
  use MLSSignals_m, only: GetSignal, Signal_T, GetSignalName
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use PointingGrid_m, only: Close_Pointing_Grid_File, &
    & Dump_Pointing_Grid_Database, Get_Grids_Near_Tan_Pressures, &
    & Open_Pointing_Grid_File, Read_Pointing_Grid_File, PointingGrids
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Gen, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named
  use VectorsModule, only: GetVectorQuantityByType, ValidateVectorQuantity, &
    & Vector_T, VectorValue_T

  !??? The next USE statement is Temporary for l2load:
  use L2_TEST_STRUCTURES_M, only: FWD_MDL_CONFIG, FWD_MDL_INFO, &
    & TEMPORARY_FWD_MDL_INFO

  implicit NONE
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

  integer :: Error            ! Error level -- 0 = OK

  real (r8), parameter :: DegToRad=1.74532925e-2

contains

  ! ------------------------------------  ForwardModelGlobalSetup  -----
  subroutine ForwardModelGlobalSetup ( Root )
  ! Process the forwardModel specification to produce ForwardModelInfo.

    integer, intent(in) :: Root         ! of the forwardModel specification.
                                        ! Indexes a "spec_args" vertex.

    integer :: Lun                      ! Unit number for pointing grid file
    character(len=80) :: PointingGridsFile   ! Duh
    integer :: Son                      ! Some subtree of root.

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( "ForwardModelGlobalSetup", root )

    ! "Root" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.
    ! Collect data from the fields.

    ! The only field so far is the PointingGrids field, and it is required.
    son = subtree(2,root)
    call get_string ( sub_rosa(subtree(2,son)), pointingGridsFile, strip=.true. )

    ! The ExtraHeights and PointingGrids fields are required, so we don't
    ! need to verify that they were provided.
    call open_pointing_grid_file ( pointingGridsFile, lun )
    call read_pointing_grid_file ( lun, spec_indices )
    call close_pointing_grid_file ( lun )

    if ( toggle(gen) ) call trace_end ( "ForwardModelGlobalSetup" )
  end subroutine ForwardModelGlobalSetup

  ! ------------------------------------------  ConstructForwardModelConfig  -----
  type (ForwardModelConfig_T) function ConstructForwardModelConfig ( ROOT ) result(info)
    ! Process the forwardModel specification to produce ForwardModelConfig to add
    ! to the database

    integer :: ROOT                     ! of the forwardModel specification.
                                        ! Indexes either a "named" or
                                        ! "spec_args" vertex.
    ! Local variables

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
    integer :: Type                     ! Type of value returned by EXPR
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
      if ( got(field) .and. (field /= f_channels) )&
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
      case ( f_spect_der )
        info%spect_der = get_boolean(son)
      case ( f_temp_der )
        info%temp_der = get_boolean(son)
      case default
        ! Shouldn't get here if the type checker worked
      end select

    end do ! i = 2, nsons(key)

    ! Now some more error checking
    select case (info%fwmType)
    case (l_full)
      if (.not. all(got( (/f_molecules, f_signals/) ))) &
        & call AnnounceError (IncompleteFullFwm, root)
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
    &                       Jacobian, RowBlock, FwdModelOut, FMC, FMI, TFMI )

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    use GL6P, only: NG
    use MLSCommon, only: I4, R4, R8
    use L2_TEST_STRUCTURES_M
    use L2PC_FILE_PARAMETERS, only: mxco => MAX_NO_ELMNTS_PER_SV_COMPONENT
    use L2PC_PFA_STRUCTURES, only: K_MATRIX_INFO
    use L2PCdim, only: Nlvl, N2lvl, NSPS, Nptg, NCH, MNP => max_no_phi, &
                       MNM => max_no_mmaf
    use ELLIPSE, only: PHI_TAN, ROC
    use L2_LOAD_M, only: L2_LOAD
    use PTG_FRQ_LOAD_M, only: PTG_FRQ_LOAD
    use COMP_PATH_ENTITIES_M, only: COMP_PATH_ENTITIES
    use REFRACTION_M, only: REFRACTION_CORRECTION
    use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR, PATH_BETA, &
                               PATH_DERIVATIVE
    use HYDROSTATIC_MODEL_M, only: HYDROSTATIC_MODEL
    use GET_CHI_ANGLES_M, only: GET_CHI_ANGLES
    use GET_BETA_PATH_M, only: GET_BETA_PATH
    use GEOC_GEOD_CONV_M, only: GEOC_GEOD_CONV
    use RAD_TRAN_M, only: RAD_TRAN
  ! use RAD_TRAN_WD_M, only: RAD_TRAN_WD
    use FREQ_AVG_M, only: FREQ_AVG
    use CONVOLVE_ALL_M, only: CONVOLVE_ALL
    use NO_CONV_AT_ALL_M, only: NO_CONV_AT_ALL
    use D_HUNT_M, only: HUNT

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    type(forwardModelConfig_T), intent(in) :: forwardModelConfig ! From ForwardModelSetup
    type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn ! ???
    type(matrix_T), intent(inout), optional :: Jacobian
    integer, intent(in), optional :: RowBlock          ! With which block of
    ! rows of F and Jacobian are we computing? All of them if absent.
    type(vector_T), intent(inout), optional :: FwdModelOut  ! Radiances, etc.

!??? Begin temporary stuff to start up the forward model
    type(fwd_mdl_config), optional :: FMC
!   type(fwd_mdl_info), dimension(:), pointer, optional :: FMI
    type(fwd_mdl_info),                        optional :: FMI
!   type(temporary_fwd_mdl_info), dimension(:), pointer, optional :: TFMI
    type(temporary_fwd_mdl_info),                        optional :: TFMI
!??? End of temporary stuff to start up the forward model

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    ! Local parameters
    character, parameter :: InvalidQuantity="Invalid vector quantity for "

    ! This is the `legit stuff' we hope will stay they are all pointers to
    ! VectorValue_T's containing vector quantities
    type (VectorValue_T), pointer :: radiance=>NULL() ! Radiance quantity to be filled
    type (VectorValue_T), pointer :: temp=>NULL() ! Temperature quantity
    type (VectorValue_T), pointer :: ptan=>NULL() ! PTAN quantity
    type (VectorValue_T), pointer :: elevOffset=>NULL() ! Elevation offset quantity
    type (VectorValue_T), pointer :: orbIncline=>NULL() ! Orbital inclination (beta)
    type (VectorValue_T), pointer :: spaceRadiance=>NULL() ! Space radiance
    type (VectorValue_T), pointer :: earthRefl=>NULL() ! Earth reflectivity
    type (VectorValue_T), pointer :: refGPH=>NULL() ! Reference GPH, (zRef and hRef)
    type (VectorValue_T), pointer :: losVel=>NULL() ! Line of sight velocity
    type (VectorValue_T), pointer :: scGeocAlt=>NULL() ! Geocentric spacecraft altitude

    Integer(i4), parameter :: ngt = (Ng+1) * N2lvl

    Integer(i4) :: i, j, k, kk, kz, ht_i, mnz, no_tan_hts, ch, Spectag, &
                   m, prev_npf, ier, mmaf, si, ptg_i, &
                   frq_i, io, klo, l, n, brkpt, no_ele, mid, ilo, ihi, &
                   k_info_count, gl_count, ld

    Integer(i4) :: ch1, ch2, no_pfa_ch, pfa_ch(2)

    Integer(i4) :: WhichPointingGrid

    Type(path_index)  :: ndx_path(Nptg,mnm)
    Type(path_vector) :: z_path(Nptg,mnm),t_path(Nptg,mnm),h_path(Nptg,mnm),  &
                         dhdz_path(Nptg,mnm), spsfunc_path(Nsps,Nptg,mnm),    &
                         n_path(Nptg,mnm),phi_path(Nptg,mnm)

    Type(path_derivative) :: dh_dt_path(Nptg,mnm)

    Real(r8) :: thbs(10),elev_offset
    Real(r8) :: t_script(N2lvl),ref_corr(N2lvl,Nptg),tau(N2lvl), &
                tan_dh_dt(Nlvl,mnm,mxco)

    Real(r8) :: dx_dt(Nptg,mxco), d2x_dxdt(Nptg,mxco)

    Real(r8) :: h_glgrid(ngt,mnm), t_glgrid(ngt,mnm), z_glgrid(ngt/2)
    Real(r8) :: dh_dt_glgrid(ngt,mnm,mxco), dhdz_glgrid(ngt,mnp)

    Real(r8) :: ptg_angles(Nptg,mnm), center_angle
    Real(r8) :: tan_hts(Nptg,mnm), tan_temp(Nptg,mnm)


!   Real(r4) :: K_TEMP(Nch,Nptg,mxco,mnp)
!   Real(r4) :: K_ATMOS(Nch,Nptg,mxco,mnp,Nsps)
!   Real(r4) :: K_SPECT_DW(Nch,Nptg,mxco,mnp,Nsps),  &
!               K_SPECT_DN(Nch,Nptg,mxco,mnp,Nsps),  &
!               K_SPECT_DNU(Nch,Nptg,mxco,mnp,Nsps)

! ** DEBUG, memory limitations force us to have up to 2 channels
!           only (Replacing Nch by: 2)

    Real(r4) :: K_TEMP(02,Nptg,mxco,mnp)
    Real(r4) :: K_ATMOS(02,Nptg,mxco,mnp,Nsps)
    Real(r4) :: K_SPECT_DW(02,Nptg,mxco,mnp,Nsps),  &
                K_SPECT_DN(02,Nptg,mxco,mnp,Nsps),  &
                K_SPECT_DNU(02,Nptg,mxco,mnp,Nsps)

    real(r8) :: I_STAR_ALL(Nch,Nptg)

    real(r4) :: K_STAR_ALL(02,20,mxco,mnp,Nptg)      ! 02 should be: Nch
    Type(k_matrix_info) :: k_star_info(20)

    Type(path_derivative) :: k_temp_frq, k_atmos_frq(Nsps), &
                             k_spect_dw_frq(Nsps), k_spect_dn_frq(Nsps), &
                             k_spect_dnu_frq(Nsps)
    !
    Type(path_beta), DIMENSION(:,:), POINTER :: beta_path => null()

    Real(r8) :: Radiances(Nptg,Nch)
    Real(r8) :: e_rad, Zeta, Frq, h_tan, Rad, geoc_lat, geod_lat, r
    !
    Character (LEN=01) :: CA
    Character (LEN=08) :: Name
    Character (LEN=16) :: Vname
    Character (LEN=80) :: Line
    Character (LEN=40) :: Ax, Dtm1, Dtm2

    Real(r8), dimension(:), allocatable :: RadV
    real(r8), dimension(:), pointer :: frequencies=>NULL()

    real(r8), parameter :: TOL = 0.001            ! How close to a tan_press
    !                                               must a frq grid be in order
    !                                               to correspond?
    integer, pointer, dimension(:) :: Grids=>NULL()       ! Frq grid for each tan_press
    type(signal_t) :: Signal

!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

    print*,'Entering forwardModel'

    ! First we identify the vector quantities we're going to need.  The key is to
    ! identify the signal we'll be working with first
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

    print*,'Got part of the way through'

    fmc%atmos_der = forwardModelConfig%atmos_der
    fmc%do_conv = forwardModelConfig%do_conv
    fmc%do_frqavg = forwardModelConfig%do_Freq_Avg
    fmc%spect_Der = forwardModelConfig%spect_Der
    fmc%temp_Der = forwardModelConfig%temp_Der
    fmc%zfrq = forwardModelConfig%the_Freq

    io = 0
    ch1 = FMC%Channels_range(1)
    ch2 = FMC%Channels_range(2)
    no_pfa_ch = min(2,ch2-ch1+1)
    do i = 1, no_pfa_ch
      pfa_ch(i) = ch1 + i - 1
    end do

    elev_offset = 0.0                         ! Zero elev_offset in any case

    signal = getSignal ( forwardModelConfig%sigInfo(1)%signal )
    whichPointingGrid = signal%pointingGrid
    if ( whichPointingGrid <= 0 ) &
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "There is no pointing grid for the desired signal" )
    call allocate_test ( grids, size(FMI%tan_press), "Grids", moduleName )
    call get_grids_near_tan_pressures ( whichPointingGrid, FMI%tan_press, &
    & tol, grids )

!
! Convert GeoDetic Latitude to GeoCentric Latitude, and convert both to
! Radians (instead of Degrees). Also compute the effective earth radius.
!
    mmaf = 3                     ! Do only this mmaf (middle phi)
    phi_tan = FMC%phi_tan_mmaf(mmaf)
    geod_lat= TFMI%t_geod_lat(mmaf)
    Call geoc_geod_conv(orbIncline%values(1,1),phi_tan,geod_lat, geoc_lat,E_rad)
!
! Compute the hydrostatic_model on the GL-Grid for all mmaf(s):
!
    thbs = 0.0
    si = FMI%Surface_index
    no_tan_hts = FMC%no_tan_hts
    thbs(1:si-1) = FMI%Tan_hts_below_surface(1:si-1)
    Call hydrostatic_model(si,FMC%N_lvls,temp%template%noSurfs, &
      & radiance%template%noInstances,FMC%t_indx, &
      & no_tan_hts,geoc_lat,refGPH%values(1,:)/1e3, &
      & refGPH%template%surfs(1,1),FMI%z_grid,thbs, &
      & temp%template%surfs(:,1), temp%values, z_glgrid, h_glgrid, t_glgrid, &
      dhdz_glgrid,dh_dt_glgrid,FMI%tan_press,tan_hts,tan_temp,tan_dh_dt, &
      gl_count, Ier)
    IF(ier /= 0) goto 99
!
! Compute all path entities for all mmafs and tanget pointings
!
    Call comp_path_entities(FMC%n_lvls,temp%template%noSurfs,&
      & gl_count,ndx_path, &
   &     z_glgrid,t_glgrid,h_glgrid,dhdz_glgrid,dh_dt_glgrid,       &
   &     TFMI%atmospheric,TFMI%f_zeta_basis,TFMI%mr_f,              &
   &     TFMI%no_coeffs_f,tan_hts,no_tan_hts,FMI%n_sps,             &
   &     TFMI%no_phi_f,TFMI%f_phi_basis,z_path,h_path,t_path,       &
   &     phi_path,n_path,dhdz_path,dh_dt_path,temp%template%noInstances,        &
   &     temp%template%phi(1,:)*degToRad,spsfunc_path,TFMI%is_f_log,FMC%no_mmaf,Ier)
    IF(ier /= 0) goto 99

!
! **********************  MAIN Mmaf Loop *******************
!
    l = mmaf
!   DO l = 1, FMC%no_mmaf
!   DO l = 1, radiance%template%noInstances   ! ** DEBUG, only one mmaf
!
      phi_tan = fmc%phi_tan_mmaf(l)
!
! Compute the ptg_angles (chi) for Antenna convolution, also the derivatives
! of chi w.r.t to T and other parameters
!
      Call get_chi_angles(ndx_path(:,l),n_path(:,l),FMI%tan_press,        &
     &     tan_hts(:,l),tan_temp(:,l),phi_tan,RoC,scGeocAlt%values(1,1),  &
     &     elevOffset%values(1,1), &
     &     tan_dh_dt(:,l,:),no_tan_hts,temp%template%noSurfs,temp%template%surfs(:,1),&
     &     si,    &
     &     center_angle,ptg_angles(:,l),dx_dt,d2x_dxdt,ier)
      IF(ier /= 0) goto 99
!
! Compute the refraction correction scaling matrix for this mmaf:
!
      Call refraction_correction(no_tan_hts, tan_hts(:,l), h_path(:,l), &
     &                n_path(:,l), ndx_path(:,l), E_rad, ref_corr)
!
      prev_npf = -1
      Radiances(1:Nptg,1:Nch) = 0.0
!
! **********************  MAIN Pointing Loop *******************
!
      DO ptg_i = 1, no_tan_hts-1

        k = ptg_i
        h_tan = tan_hts(k,l)
        kk = SIZE(PointingGrids(whichPointingGrid)%OneGrid(grids(ptg_i))%Frequencies)
        FMI%no_ptg_frq(k)=kk !??? Remove this later...
!        kk = FMI%no_ptg_frq(k)
!
        if(kk /= prev_npf) then
!
          prev_npf = kk
          DEALLOCATE(k_temp_frq%values,STAT=i)
!
          DEALLOCATE(RadV,STAT=i)
          ALLOCATE(RadV(kk),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: RadV array, STAT =',ier
            goto 99
          endif
!
          do j = 1, FMI%n_sps
            DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
            DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
            DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
            DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
          end do
!
          ALLOCATE(k_temp_frq%values(kk,temp%template%noSurfs,&
            & temp%template%noInstances),STAT=ier)
          IF(ier /= 0) then
            Print *,'** ALLOCATE Error: k_temp_frq array, STAT =',ier
            goto 99
          endif
!
          do j = 1, FMI%n_sps
            m = max(1,TFMI%no_phi_f(j))
            i = max(1,TFMI%no_coeffs_f(j))
            ALLOCATE(k_atmos_frq(j)%values(kk,i,m),STAT=ier)
            IF(ier /= 0) then
              Print *,'** ALLOCATE Error: k_atmos_frq, STAT =',ier
              goto 99
            endif
          end do
!
          do m = 1, FMI%n_sps
            j = FMI%spect_atmos(m)
            if(j < 1) CYCLE
            if(.not. FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
            Vname = ' '
            Spectag = FMI%spectroscopic(j)%Spectag
            DO
              if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
              n = FMI%spectroscopic(j)%no_phi_values
              i = FMI%spectroscopic(j)%no_zeta_values
              CA = FMI%spectroscopic(j)%type
              select case ( CA )
                case ( 'W' )
                  Vname = 'k_spect_dw_frq'
                  ALLOCATE(k_spect_dw_frq(m)%values(kk,i,n),STAT=ier)
                case ( 'N' )
                  Vname = 'k_spect_dn_frq'
                  ALLOCATE(k_spect_dn_frq(m)%values(kk,i,n),STAT=ier)
                case ( 'V' )
                  Vname = 'k_spect_dnu_frq'
                  ALLOCATE(k_spect_dnu_frq(m)%values(kk,i,n),STAT=ier)
                case default
                  Ier = -99
                  Print *,'** Unknown Spectroscopic element !'
              end select
              IF(ier /= 0) then
                Print *,'** ALLOCATE Error: ',Vname,', STAT =',ier
                goto 99
              ENDIF
              j = j + 1
              if(j > 3 * FMI%n_sps) EXIT
            END DO
          end do

        endif            ! On DEALLOCATE/ALLOCATE cycle
!
! Compute the beta's along the path, for this tanget hight and this mmaf:
!
        no_ele = ndx_path(ptg_i,l)%total_number_of_elements

        if (FMC%do_frqavg) then
          frequencies => PointingGrids(whichPointingGrid)%oneGrid(grids(ptg_i))%frequencies
        else
          call allocate_test(frequencies,1,"frequencies",ModuleName)
          frequencies(1)=FMC%zfrq
        endif
        kk = size(frequencies)

        Call get_beta_path(frequencies,&
          & FMI%pfa_spectrum,no_ele, z_path(ptg_i,l),t_path(ptg_i,l), &
          & beta_path, FMC%vel_z_mmaf(l),ier)
        IF(ier /= 0) goto 99
!
        k_temp_frq%values = 0.0
        do j = 1, FMI%n_sps
          k_atmos_frq(j)%values = 0.0
          k_spect_dw_frq(j)%values = 0.0
          k_spect_dn_frq(j)%values = 0.0
          k_spect_dnu_frq(j)%values = 0.0
        end do
!
          
        do frq_i = 1, kk
!
          Frq = frequencies(frq_i)
!
!           call output("For pointing: ")
!           call output(ptg_i)
!           call output(" and frequency ")
!           call output(frq_i,advance='yes')
!           call output("beta_bath",advance='yes')
!           call dump(beta_path(1,frq_i)%values)
          
          Call Rad_Tran(Frq, FMC%N_lvls, h_tan, FMI%n_sps, ndx_path(k,l),  &
         &    z_path(k,l), h_path(k,l), t_path(k,l), phi_path(k,l),&
         &    dHdz_path(k,l), earthRefl%values(1,1), beta_path(:,frq_i),      &
         &    spsfunc_path(:,k,l), ref_corr(:,k), spaceRadiance%values(1,1), brkpt, &
         &    no_ele, mid, ilo, ihi, t_script, tau, Rad, Ier)
          IF(ier /= 0) goto 99
!
          RadV(frq_i) = Rad
!
! Now, Compute the radiances derivatives:
!
!         CALL Rad_Tran_WD(frq_i,FMI%band,Frq,FMC%N_lvls,FMI%n_sps, &
!        &     FMC%temp_der,FMC%atmos_der,FMC%spect_der,            &
!        &     z_path(k,l),h_path(k,l),t_path(k,l),phi_path(k,l),   &
!        &     dHdz_path(k,l),TFMI%atmospheric,beta_path(:,frq_i),&
!        &     spsfunc_path(:,k,l),temp%template%surfs(:,1),  &
!        &     TFMI%f_zeta_basis,TFMI%no_coeffs_f,   &
!        &     TFMI%mr_f,TFMI%no_t,ref_corr(:,k),TFMI%no_phi_f,       &
!        &     TFMI%f_phi_basis,temp%template%noInstances,temp%template%phi(1,:)*degToRad,  &
!        &     dh_dt_path(k,l),FMI%spect_atmos,         &
!        &     FMI%spectroscopic,k_temp_frq,k_atmos_frq,k_spect_dw_frq,   &
!        &     k_spect_dn_frq,k_spect_dnu_frq,TFMI%is_f_log,brkpt,       &
!        &     no_ele,mid,ilo,ihi,t_script,tau,ier)
!         IF(ier /= 0) goto 99
!
        end do

!
! Frequency Average the radiances with the appropriate filter shapes
!
        do i = 1, no_pfa_ch
          ch = pfa_ch(i)
          if(FMC%do_frqavg) then
            print*,'Doing frequency averaging'
            Call Freq_Avg(frequencies,FMI%F_grid_filter(:,i),  &
           &     FMI%Filter_func(:,i),RadV,kk,FMI%no_filt_pts, &
           &     Radiances(ptg_i,ch))
          else
            Radiances(ptg_i,ch) = RadV(1)
          endif
        end do
!
        if(FMC%temp_der) then
!
! Frequency Average the temperature derivatives with the appropriate
! filter shapes
!
          RadV(1:kk) = 0.0
          do i = 1, no_pfa_ch
!           ch = pfa_ch(i)
            ch = i               ! ** DEBUG, memory limitations on MLSGATE
            do j = 1, temp%template%noInstances
              do k = 1, temp%template%noSurfs
                if(FMC%do_frqavg) then
                  RadV(1:kk) = k_temp_frq%values(1:kk,k,j)
                  Call Freq_Avg(frequencies,FMI%F_grid_filter(:,i), &
                 &              FMI%Filter_func(:,i),          &
                 &              RadV,kk,FMI%no_filt_pts,r)
                else
                  r = k_temp_frq%values(1,k,j)
                endif
                k_temp(ch,ptg_i,k,j) = r
              end do
            end do
          end do
!
        endif

        if(FMC%atmos_der) then
!
! Frequency Average the atmospheric derivatives with the appropriate
! filter shapes
!
          do i = 1, no_pfa_ch
!           ch = pfa_ch(i)
            ch = i                 ! ** DEBUG, memory limitations on MLSGATE
            do j = 1, FMI%n_sps
              if(TFMI%atmospheric(j)%der_calc(FMI%band)) THEN
                RadV(1:kk) = 0.0
                do k = 1, TFMI%no_phi_f(j)
                  do n = 1, TFMI%no_coeffs_f(j)
                    if(FMC%do_frqavg) then
                      RadV(1:kk) = k_atmos_frq(j)%values(1:kk,n,k)
                      Call Freq_Avg(frequencies,FMI%F_grid_filter(:,i), &
                     &              FMI%Filter_func(:,i),          &
                     &              RadV,kk,FMI%no_filt_pts,r)
                    else
                      r = k_atmos_frq(j)%values(1,n,k)
                    endif
                    k_atmos(ch,ptg_i,n,k,j) = r
                  end do
                end do
              endif
            end do
          end do
!
        endif

        if(FMC%spect_der) then
!
! Frequency Average the spectroscopic derivatives with the appropriate
! filter shapes
!
          do i = 1, no_pfa_ch
!           ch = pfa_ch(i)
            ch = i                 ! ** DEBUG, memory limitations on MLSGATE
            do m = 1, FMI%n_sps
              j = FMI%spect_atmos(m)
              if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
              Spectag = FMI%spectroscopic(j)%Spectag
              DO
                if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
                RadV(1:kk) = 0.0
                CA = FMI%spectroscopic(j)%type
                do k = 1, FMI%spectroscopic(j)%no_phi_values
                  do n = 1, FMI%spectroscopic(j)%no_zeta_values
                    select case ( CA )
                      case ( 'W' )
                        RadV(1:kk) = k_spect_dw_frq(m)%values(1:kk,n,k)
                      case ( 'N' )
                        RadV(1:kk) = k_spect_dn_frq(m)%values(1:kk,n,k)
                      case ( 'V' )
                        RadV(1:kk) = k_spect_dnu_frq(m)%values(1:kk,n,k)
                    end select
                    if(FMC%do_frqavg) then
                        Call Freq_Avg(frequencies,FMI%F_grid_filter(:,i), &
                     &              FMI%Filter_func(:,i),&
                     &              RadV,kk,FMI%no_filt_pts,r)
                    else
                      r = RadV(1)
                    endif
                    select case ( CA )
                      case ( 'W' )
                        k_spect_dw(ch,ptg_i,n,k,j) = r
                      case ( 'N' )
                        k_spect_dn(ch,ptg_i,n,k,j) = r
                      case ( 'V' )
                        k_spect_dnu(ch,ptg_i,n,k,j) = r
                    end select
                  end do
                end do
                j = j + 1
                if(j > 3 * FMI%n_sps) EXIT
              END DO
            end do
          end do
!
        endif
!
      END DO              ! Pointing Loop

!
! Complete the radiances's last location, also  complete k_temp last
! location as well as k_atmos last location and k_spect_d? last location:
!
      kk = no_tan_hts
      do i = 1, no_pfa_ch
        ch = pfa_ch(i)
        Radiances(kk,ch) = Radiances(kk-1,ch)
        if(FMC%temp_der) then
          k_temp(i,kk,1:temp%template%noSurfs,1:temp%template%noInstances) = &
       &              k_temp(i,kk-1,1:temp%template%noSurfs,&
       &              1:temp%template%noInstances)
        endif
        if(FMC%atmos_der) then
          do m = 1, FMI%n_sps
            if(TFMI%atmospheric(m)%der_calc(FMI%band)) then
              k = TFMI%no_phi_f(m)
              n = TFMI%no_coeffs_f(m)
              k_atmos(i,kk,1:n,1:k,m)=k_atmos(i,kk-1,1:n,1:k,m)
            endif
          end do
        endif
        if(FMC%spect_der) then
          do m = 1, FMI%n_sps
            j = FMI%spect_atmos(m)
            if(.not.  FMI%spectroscopic(j)%DER_CALC(FMI%band)) CYCLE
            Spectag =  FMI%spectroscopic(j)%Spectag
            DO
              if(FMI%spectroscopic(j)%Spectag /= Spectag) EXIT
              k = FMI%spectroscopic(j)%no_phi_values
              n = FMI%spectroscopic(j)%no_zeta_values
              k_spect_dw(i,kk,1:n,1:k,j)=k_spect_dw(i,kk-1,1:n,1:k,j)
              k_spect_dn(i,kk,1:n,1:k,j)=k_spect_dn(i,kk-1,1:n,1:k,j)
              k_spect_dnu(i,kk,1:n,1:k,j)=k_spect_dnu(i,kk-1,1:n,1:k,j)
              j = j + 1
              if(j > 3 * FMI%n_sps) EXIT
            END DO
          end do
        endif
      end do
!
!  Here comes the Convolution code
!
      DO i = 1, no_pfa_ch
!
        ch = pfa_ch(i)
!
        if(FMC%do_conv) then
!
          Call convolve_all(ptan%values(:,l),TFMI%atmospheric,FMI%n_sps,   &
         &     FMC%temp_der,FMC%atmos_der,FMC%spect_der,                   &
         &     FMI%tan_press,ptg_angles(:,l),tan_temp(:,l), &
         &     dx_dt, d2x_dxdt,FMI%band,center_angle,FMI%fft_pts,          &
         &     Radiances(:,ch),k_temp(i,:,:,:),k_atmos(i,:,:,:,:), &
         &     k_spect_dw(i,:,:,:,:),k_spect_dn(i,:,:,:,:),    &
         &     k_spect_dnu(i,:,:,:,:),FMI%spect_atmos,no_tan_hts,  &
         &     k_info_count,i_star_all(i,:),k_star_all(i,:,:,:,:), &
         &     k_star_info,temp%template%noSurfs,temp%template%noInstances,&
         &     TFMI%no_phi_f,   &
         &     FMI%spectroscopic,temp%template%surfs(:,1),FMI%Xlamda,FMI%Aaap,&
         &     FMI%D1Aaap,FMI%D2Aaap,FMI%Ias,ier)
          IF(ier /= 0) goto 99
!
        else
!
          Call no_conv_at_all(ptan%values(:,l),FMI%n_sps,FMI%tan_press, &
         &     FMI%band,FMC%temp_der,FMC%atmos_der,FMC%spect_der,      &
         &     Radiances(:,ch),k_temp(i,:,:,:),                    &
         &     k_atmos(i,:,:,:,:),k_spect_dw(i,:,:,:,:),       &
         &     k_spect_dn(i,:,:,:,:),k_spect_dnu(i,:,:,:,:),   &
         &     FMI%spect_atmos, no_tan_hts,k_info_count,               &
         &     i_star_all(i,:), k_star_all(i,:,:,:,:),            &
         &     k_star_info,temp%template%noSurfs,temp%template%noInstances, &
         &     TFMI%no_phi_f,temp%template%surfs(:,1),TFMI%atmospheric,    &
         &     FMI%spectroscopic)
!
        endif
!
      END DO
!
!   END DO                ! Mmaf Loop
!
    if(FMC%temp_der) DEALLOCATE(k_temp_frq%values,STAT=i)
    do j = 1, FMI%n_sps
      if(FMC%atmos_der) DEALLOCATE(k_atmos_frq(j)%values,STAT=i)
      if(FMC%spect_der) then
        DEALLOCATE(k_spect_dw_frq(j)%values,STAT=i)
        DEALLOCATE(k_spect_dn_frq(j)%values,STAT=i)
        DEALLOCATE(k_spect_dnu_frq(j)%values,STAT=i)
      endif
    end do
!
! *** DEBUG Print
!
    if(FMC%do_conv) then
      Print *,'Convolution: ON'
    else
      Print *,'Convolution: OFF'
    endif
!
    if(FMC%Zfrq > 0.0) then
      Frq = FMC%Zfrq
      write(*,901) Frq
901   format(' Frequency Averaging: OFF',/,  &
        &    ' (All computations done at Frq =',f12.4,')')
    else
      Print *,'Frequency Averaging: ON'
    endif
    Print *
    if (.not. FMC%do_frqavg) deallocate(frequencies)
!
    tau(1:Nptg) = 0.0
    kk = TFMI%ptg_press%no_lin_values
    tau(1:kk) = dble(TFMI%ptg_press%lin_val(1:kk))

    klo = -1
    Zeta = -1.666667
    Call Hunt(Zeta,tau,kk,klo,j)
    IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(klo))) klo=j
!
    do i = 1, no_pfa_ch
      ch = pfa_ch(i)
      write(*,903) ch,char(92),kk
      write(*,905) (i_star_all(i,k),k=1,kk)
    end do
903 format('ch',i2.2,'_pfa_rad',a1,i3.3)
905 format(4(2x,1pg15.8))
!
    if(.not. ANY((/FMC%temp_der,FMC%atmos_der,FMC%spect_der/))) goto 99
!
    m = -1
    ch = 1
    tau = 0.0
    do i = 1, k_info_count
      Print *
      Name = ' '
      Name = k_star_info(i)%name
      if(Name == 'PTAN') CYCLE
      kz = k_star_info(i)%first_dim_index
      mnz = k_star_info(i)%no_zeta_basis
      ht_i = k_star_info(i)%no_phi_basis
      l = LEN_TRIM(Name)
      if(Name(l-1:l) == '_W' .or.  &
     &   Name(l-1:l) == '_N' .or.  &
     &   Name(l-1:l) == '_V' ) then
        Print *,Name
        r = SUM(k_star_all(ch,kz,1:mnz,1:ht_i,klo))
        Print *,'  Sum over all zeta & phi coeff:',sngl(r)
      else
        if(Name == 'TEMP') then
          write(6,913) 'dI_dT',char(92),ht_i
        else
          write(6,913) 'dI_d'//Name(1:l),char(92),ht_i
        endif
        tau = 0.0
        tau(1:mnz) = k_star_info(i)%zeta_basis(1:mnz)
        Call Hunt(Zeta,tau,mnz,m,j)
        IF(ABS(Zeta-tau(j)) < ABS(Zeta-tau(m))) m=j
        Print *,(k_star_all(ch,kz,m,kk,klo),kk=1,ht_i)
      endif
    end do
913 format(a,a1,i2.2)
!
 99  if (io /= 0) then
       Call ErrMsg(Line,io)
    endif

    Return
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz

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
          if (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Deallocate//"database%signals%channelIncluded")
        end do
        deallocate (database(config)%sigInfo, stat=status)
        if (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
          & MLSMSG_Deallocate//"database%signals")
        deallocate (database(config)%molecules, stat=status)
        if (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
          & MLSMSG_Deallocate//"database%molecules")
      end do
      
      deallocate (database, stat=status)
      if (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
        & MLSMSG_Deallocate//"database")
    end if
  end subroutine DestroyFWMConfigDatabase

! =====     Private Procedures     =====================================
  ! ------------------------------------  DUMP_FOWARDMODELCONFIGS  -----
  subroutine Dump_ForwardModelConfigs ( database )
    type (ForwardModelConfig_T), dimension(:), pointer :: database

    ! Local variables
    integer :: I,J                      ! Loop counters
    character (len=80) :: signalName    ! A line of text

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
  subroutine AnnounceError ( Code, Where, FieldIndex )
    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: Where      ! Where in the tree did the error occur?
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
    end select
  end subroutine AnnounceError

end module ForwardModelInterface

! $Log$
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
