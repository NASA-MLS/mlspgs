! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ForwardModelInterface
  !=============================================================================

  ! Set up the forward model.  Interface from the retrieve step to the
  ! forward model.

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
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate,&
     & MLSMSG_Error
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
! Revision 2.137  2001/05/29 23:09:59  livesey
! Last version, before being renamed! :-(
!
! Revision 2.136  2001/05/25 20:25:56  livesey
! Now optionally skips MAFs in overlap regions
!
! Revision 2.135  2001/05/21 23:58:10  livesey
! Moved some of the emits to higher numbers
!
! Revision 2.134  2001/05/17 01:44:16  livesey
! Bug fix for non frequency avg. case
!
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
