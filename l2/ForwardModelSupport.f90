! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelSupport

  ! Set up the forward model stuff.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use ForwardModelConfig, only: AddForwardModelConfigToDatabase, Dump, &
    & ForwardModelConfig_T
  use Init_Tables_Module, only: FIELD_FIRST, FIELD_LAST
  use Init_Tables_Module, only: L_FULL, L_SCAN, L_LINEAR
  use Init_Tables_Module, only: F_ANTENNAPATTERNS, F_ATMOS_DER, F_CHANNELS, &
    & F_DO_CONV, F_DO_FREQ_AVG, F_FILTERSHAPES, F_FREQUENCY, F_FRQGAP,&
    & F_INTEGRATIONGRID, F_L2PC, F_MOLECULES, F_MOLECULEDERIVATIVES, F_PHIWINDOW, &
    & F_POINTINGGRIDS, F_SIGNALS, F_SPECT_DER, F_TANGENTGRID, F_TEMP_DER, F_TYPE,&
    & F_MODULE, F_SKIPOVERLAPS
  use MLSFiles, only: GetPCFromRef, MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF
  use MLSCommon, only: R8
  use MLSL2Options, only: PCF, PCFL2CFSAMECASE, PUNISH_FOR_INVALID_PCF
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
     & MLSMSG_Error
  use MLSPCF2, only: mlspcf_antpats_start, mlspcf_filtshps_start, &
     &          mlspcf_ptggrids_start
  use MoreTree, only: Get_Boolean, Get_Field_ID
  use Output_M, only: Output
  use Parse_Signal_m, only: PARSE_SIGNAL
  use String_Table, only: Display_String, Get_String
  use Toggles, only: Emit, Gen, Levels, Switches, Toggle
  use Trace_M, only: Trace_begin, Trace_end
  use Tree, only: Decoration, Node_ID, Nsons, Source_Ref, Sub_Rosa, Subtree
  use Tree_Types, only: N_named
  use Units, only: Deg2Rad, PHYQ_FREQUENCY
  use VGridsDatabase, only: VGrid_T
  use PointingGrid_m, only: READ_POINTING_GRID_FILE, CLOSE_POINTING_GRID_FILE
  use L2PC_m, only: OPEN_L2PC_FILE, CLOSE_L2PC_FILE, READ_L2PC_FILE
  use AntennaPatterns_m, only: OPEN_ANTENNA_PATTERNS_FILE, READ_ANTENNA_PATTERNS_FILE,&
    & CLOSE_ANTENNA_PATTERNS_FILE
  use FilterShapes_m, only: OPEN_FILTER_SHAPES_FILE, READ_FILTER_SHAPES_FILE,&
    & CLOSE_FILTER_SHAPES_FILE
  use Expr_M, only: EXPR
  use Lexer_Core, only: PRINT_SOURCE
  use MLSNumerics, only: HUNT
  use PointingGrid_m, only: Close_Pointing_Grid_File, &
    & Open_Pointing_Grid_File, Read_Pointing_Grid_File


  implicit none
  private
  public :: ConstructForwardModelConfig, ForwardModelGlobalSetup

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
  subroutine ForwardModelGlobalSetup ( Root, any_errors )
    ! Process the forwardModel specification to produce ForwardModelInfo.

    integer, intent(in) :: Root         ! of the forwardModel specification.
    !                                     Indexes a "spec_args" vertex.
    integer, intent(out) :: any_errors  ! non-zero means trouble

    logical, parameter :: DEBUG = .false.
    character(len=255) :: FileName      ! Duh
    integer :: I                        ! Loop inductor, subscript
    integer :: Lun                      ! Unit number for reading a file
    character(len=255) :: PCFFileName
    integer :: returnStatus             ! non-zero means trouble
    integer :: Son                      ! Some subtree of root.
    integer :: Version

    ! Error message codes

    error = 0
    if ( toggle(gen) ) call trace_begin ( 'ForwardModelGlobalSetup', root )

    ! "Root" now indexes an n_spec_args vertex.  See "Configuration file
    ! parser users' guide" for pictures of the trees being analyzed.
    ! Collect data from the fields.

    do i = 2, nsons(root)
      Version = 1
      son = subtree(i,root)
      select case ( get_field_id(son) )
      case ( f_antennaPatterns )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if ( PCF ) then
            lun = GetPCFromRef(fileName, mlspcf_antpats_start, &
            & mlspcf_antpats_start, &
            & PCFL2CFSAMECASE, returnStatus, Version, DEBUG, &
            & exactName=PCFFileName)
            if ( returnStatus /= 0 .and. PUNISH_FOR_INVALID_PCF) then
               call AnnounceError(0, son, &
               & extraMessage='Antenna Patterns File not found in PCF')
            elseif( returnStatus == 0) then
               fileName = PCFFileName
            endif
        endif
        call open_antenna_patterns_file ( fileName, lun )
        call read_antenna_patterns_file ( lun )
        call close_antenna_patterns_file ( lun )
      case ( f_filterShapes )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if ( PCF ) then
            lun = GetPCFromRef(fileName, mlspcf_filtshps_start, &
            & mlspcf_filtshps_start, &
            & PCFL2CFSAMECASE, returnStatus, Version, DEBUG, &
            & exactName=PCFFileName)
            if ( returnStatus /= 0 .and. PUNISH_FOR_INVALID_PCF) then
               call AnnounceError(0, son, &
               & extraMessage='Filter Shapes File not found in PCF')
            elseif( returnStatus == 0) then
               fileName = PCFFileName
            endif
        endif
        call open_filter_shapes_file ( fileName, lun )
        call read_filter_shapes_file ( lun )
        call close_filter_shapes_file ( lun )
      case ( f_pointingGrids )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        if ( PCF ) then
            lun = GetPCFromRef(fileName, mlspcf_ptggrids_start, &
            & mlspcf_ptggrids_start, &
            & PCFL2CFSAMECASE, returnStatus, Version, DEBUG, &
            & exactName=PCFFileName)
            if ( returnStatus /= 0 .and. PUNISH_FOR_INVALID_PCF) then
               call AnnounceError(0, son, &
               & extraMessage='Pointing Grids File not found in PCF')
            elseif( returnStatus == 0) then
               fileName = PCFFileName
            endif
        endif
        call open_pointing_grid_file ( fileName, lun )
        call read_pointing_grid_file ( lun )
        call close_pointing_grid_file ( lun )
      case ( f_l2pc )
        call get_string ( sub_rosa(subtree(2,son)), fileName, strip=.true. )
        ! l2pc files not in PCF (yet)
        call open_l2pc_file (fileName, lun)
        call read_l2pc_file ( lun )
        call close_l2pc_file ( lun )
      case default
        ! Can't get here if the type checker worked
      end select
    end do

    if ( toggle(gen) ) call trace_end ( 'ForwardModelGlobalSetup' )
    any_errors = error
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
    info%skipOverlaps = .false.
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
      case ( f_skipOverlaps )
        info%skipOverlaps = get_boolean(son)
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

  ! =====     Private Procedures     =====================================
  ! ----------------------------------------------  AnnounceError  -----
  subroutine AnnounceError ( Code, where, FieldIndex, extraMessage )
    integer, intent(in) :: Code       ! Index of error message
    integer, intent(in) :: where      ! Where in the tree did the error occur?
    integer, intent(in), optional :: FieldIndex ! f_...
    character (LEN=*), optional :: extraMessage

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
    case default
      call output ( '(no specific description of this error)', advance='yes' )
    end select
    if ( present(extraMessage) ) call output ( extraMessage, advance='yes' )
  end subroutine AnnounceError

end module ForwardModelSupport

! $Log$
! Revision 2.2  2001/05/30 23:05:39  pwagner
! Uses PCF for 3 fwdmdl files
!
! Revision 2.1  2001/05/29 23:18:18  livesey
! First version, was ForwardModelInterface
!
