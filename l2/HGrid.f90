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
module HGrid                    ! Horizontal grid information
!=============================================================================

  use HGridsDatabase, only: HGrid_T, HGridGeolocations, &
    & L1BGeolocation, L1BSubsample
  use MLSCommon, only: MLSFile_T, NameLen, TAI93_Range_T
  use MLSFiles, only: HDFVersion_5, GetMLSFileByType
  use MLSFinds, only: FindAll
  use MLSKinds, only: Rk => R8, R8
  use MLSSignals_M, only: GetModuleName
  
  implicit none
  private
  public :: CreateHGridfromMLSCFInfo, ComputeNextChunksHGridOffsets, &
   & ComputeAllHGridOffsets, DealWithObstructions, DestroyHGridGeoLocations

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! -----     Private declarations     ---------------------------------

  logical, save, private :: computingOffsets = .false.
  integer, private       :: error                     

  interface PlaceArray
    module procedure PlaceArray_r4
    module procedure PlaceArray_r8
  end interface

  real :: T0, T1, T2               ! For timing

! Error codes for "announce_error"
  integer, private, parameter :: BadTime = 1
  integer, private, parameter :: NoFraction = BadTime + 1
  integer, private, parameter :: NoHeight = NoFraction + 1
  integer, private, parameter :: NoL1Bfiles = NoHeight + 1
  integer, private, parameter :: NoMIF = NoL1Bfiles + 1
  integer, private, parameter :: NoModule = NoMIF + 1
  integer, private, parameter :: NoPolygon = NoModule + 1
  integer, private, parameter :: NoSpacingOrigin = NoPolygon + 1


contains ! =====     Public Procedures     =============================

  ! -----------------------------------  CreateHGridFromMLSCFInfo  -----
  type(hGrid_T) function CreateHGridFromMLSCFInfo &
    & ( name, root, filedatabase, l2gpDatabase, &
    & processingRange, chunk, onlyComputingOffsets, check ) result ( hGrid )

    use Allocate_Deallocate, only: Deallocate_Test
    use Chunks_M, only: MLSChunk_T
    use Constants, only: Ln2
    use Dates_Module, only: Tai93s2hid
    use Expr_M, only: Expr
    use HGridsDatabase, only: HGrid_T, CreateEmptyHGrid, NullifyHGrid
    use HighOutput, only: BeVerbose, LetsDebug, OutputNamedValue
    use Init_Tables_Module, only: F_Coordinate, F_Date, &
      & F_Extendible, F_Forbidoverspill, F_Fraction, F_Geodangle, F_Geodlat, &
      & F_Height, F_Inclination, F_Insetoverlaps, F_Interpolationfactor, &
      & F_Lon, F_Losangle, F_Maxloweroverlap, F_Maxupperoverlap, F_Mif, &
      & F_Module, F_Origin, F_QTMlevel, &
      & F_Single, F_Solartime, F_Solarzenith, F_Sourcel2GP, F_Spacing, &
      & F_Time, F_Type, &
      & Field_First, Field_Last, &
      & L_Explicit, L_Fixed, L_Fractional, L_Height, &
      & L_L2gp, L_QTM, L_Regular
    use Intrinsic, only: PHYQ_Angle, PHYQ_Dimensionless, PHYQ_Length
    use L1BData, only: CheckForCorruptFileDatabase
    use L2GPData, only: L2GPData_T
    use MLSHDF5, only: IsHDF5DSInFile
    use MLSL2options, only: Need_L1BFiles
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: Hunt
    use MLSStringLists, only: SwitchDetail
    use MoreMessage, only: MLSMessage
    use MoreTree, only: Get_Boolean
    use Output_M, only: Output
    use Polygon_M, only: Polygon_Inside, Polygon_Vertices
    use QTM_M, only: QTM_Depth ! Maximum Depth That Will Fit In One Integer
    ! use String_Table, only: Get_String
    use Time_M, only: SayTime, Time_Now
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Decoration, NSons, Sub_Rosa, Subtree, Where

  ! This routine creates an hGrid based on the user requests.

    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of hGrid subtree
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (L2GPData_T), pointer, dimension(:) :: L2GPDATABASE
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    ! logical, intent(in), optional :: SUPPRESSGEOMETRYDUMP
    logical, intent(in), optional :: onlyComputingOffsets
    logical, intent(in), optional :: check
    type (MLSChunk_T), intent(in) :: CHUNK ! The chunk

    ! Local variables
    integer :: DATE                     ! Tree node
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    real(rk) :: fraction, height
    integer :: hGridType
    real(rk) :: interpolationFactor
    integer :: instrumentModule
    type (L2GPData_T), pointer :: L2GP  ! The l2gp to use 
    integer :: Me = -1                  ! String index for trace
    real(rk) :: spacing, origin

    integer :: keyNo                    ! Entry in the mlscf line
    integer :: fieldValue               ! Node in the tree

    real(rk) :: MINTIME, MAXTIME        ! Span for a chunk
    logical :: MYSUPPRESSGEOMETRYDUMP
    integer, dimension(2) :: PROFRANGE  ! Profile range
    integer :: A,B                      ! Elements of profile range
    logical :: DeeBug
    logical :: EXTENDIBLE               ! If set don't lose profiles between chunks
    integer :: FIELD                    ! Subtree index of "field" node
    integer :: FIELD_INDEX              ! F_..., see Init_Tables_Module
    double precision, dimension(:), pointer :: FullArray
    integer :: GEODANGLENODE            ! Tree node
    integer :: GEODLATNODE            ! Tree node
    logical :: GOT_FIELD(field_first:field_last)
    real(rk) :: Incline                 ! Orbit inclination, degrees
    logical :: INSETOVERLAPS            ! Flag
    integer :: LONNODE                  ! Tree node
    integer :: LOSANGLENODE             ! Tree node
    integer :: MAXLOWEROVERLAP          ! For some hGrids
    integer :: MAXUPPEROVERLAP          ! For some hGrids
    integer :: MIF                      ! For fixed hGrids
    logical :: myCheck
    integer :: NOMAFS                   ! Number of MAFs of L1B data read
    integer :: QTM_Level                ! From QTMLevel field
    logical :: SINGLE                   ! Just one profile please
    integer :: SON                      ! Son of Root
    integer :: SOLARTIMENODE            ! Tree node
    integer :: SOLARZENITHNODE          ! Tree node
    double precision, dimension(:), pointer :: TAI
    integer :: TIMENODE                 ! Tree node
    logical :: verbose

    character (len=NameLen) :: InstrumentModuleName

    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code
    call trace_begin ( me, "CreateHGridFromMLSCFInfo", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

    nullify ( fullArray, TAI )

    deebug = LetsDebug ( 'hgrid', 2 )
    verbose = beVerbose ( 'hgrid', 1 )
    computingOffsets       = .false.               
    mySuppressGeometryDump = .false.               
    if ( present ( onlyComputingOffsets ) ) then   
      mySuppressGeometryDump = onlyComputingOffsets
      computingOffsets       = onlyComputingOffsets
    endif
    if ( verbose ) then
      call time_now ( t1 )
    endif
    myCheck = .false.
    if ( present(check) ) myCheck = check

    L1BFile => GetMLSFileByType( filedatabase, content='l1boa' )
    if ( .not. associated(L1BFile) .and. NEED_L1BFILES ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
      & // 'CreateHGridFromMLSCFInfo', &
      & "Didn't I warn you about not having an L1BOA file?" )

    call nullifyHGrid ( hgrid ) ! for Sun's rubbish compiler
    if ( mycheck ) then
      call outputNamedValue( 'Before creating HGrid', chunk%chunkNumber )
      call CheckForCorruptFileDatabase( filedatabase )
    endif
    
    hGrid%name = name

    date = 0
    extendible = .false.
    geodAngleNode = 0
    geodLatNode = 0
    got_field = .false.
    incline = 90 ! Latitude same as phi, not projected onto orbit plane
    insetOverlaps = .false.
    interpolationFactor = 1.0
    lonNode = 0
    LosAngleNode = 0
    maxLowerOverlap = 999 ! was -1 but now -1 is a legal value
    maxUpperOverlap = 999 ! was -1 but ..
    QTM_level = 6             ! QTM level 6 gives 156 km resolution
    single = .false.
    solarTimeNode = 0
    solarZenithNode = 0
    timeNode = 0

    do keyNo = 2, nsons(root)
      son = subtree(keyNo,root)
      field = subtree(1,son)
      if ( nsons(son) > 1 ) then
        fieldValue = decoration(subtree(2,son)) ! The field's value
      else
        fieldValue = son
      end if
      field_index = decoration(field)
      ! The tree checker prevents duplicate fields
      got_field(field_index) = .true.

      select case ( field_index )
      case ( f_type ) ! Required
        hGridType = decoration(subtree(2,son))
      case ( f_coordinate )
        hGrid%masterCoordinate = decoration(subtree(2,son))
      case ( f_date )
        date = sub_rosa(subtree(2,son))
      case  ( f_extendible )
        extendible = get_boolean ( fieldValue )
      case  ( f_forbidOverspill )
        hGrid%forbidOverspill = get_boolean ( fieldValue )
      case ( f_fraction )
        call expr ( subtree(2,son), expr_units, expr_value )
        fraction = expr_value(1)
      case ( f_height )
        call expr ( subtree(2,son), expr_units, expr_value )
        height = expr_value(1)
      case ( f_insetOverlaps )
        insetOverlaps = get_boolean ( fieldValue )
      case ( f_interpolationFactor )
        call expr ( subtree(2,son), expr_units, expr_value )
        interpolationFactor = expr_value(1)
      case ( f_geodAngle )
        geodAngleNode = son
      case ( f_geodLat )
        geodLatNode = son
      case ( f_inclination )
        call expr (subtree ( 2, son), expr_units, expr_value )
        incline = expr_value(1)
      case ( f_lon )
        lonNode = son
      case ( f_losAngle )
        losAngleNode = son
      case ( f_maxLowerOverlap )
        call expr ( subtree(2,son), expr_units, expr_value )
        maxLowerOverlap = nint(expr_value(1))
      case ( f_maxUpperOverlap )
        call expr ( subtree(2,son), expr_units, expr_value )
        maxUpperOverlap = nint(expr_value(1))
      case ( f_mif )
        call expr ( subtree(2,son), expr_units, expr_value )
        mif = nint(expr_value(1))
      case ( f_module )
        ! instrumentModule = sub_rosa(subtree(2,son))
        ! call get_string ( instrumentModule , instrumentModuleName )
        instrumentModule = decoration(decoration(subtree(2,son)))
        call GetModuleName ( instrumentModule , instrumentModuleName )
        hGrid%module = instrumentModule
        if ( deebug ) then
          call output('insrumentModuleName in the case f_module is ' )
          call output(instrumentModuleName, advance='yes')
        endif
      case ( f_origin )
        call expr ( subtree(2,son), expr_units, expr_value )
        origin = expr_value(1)
      case ( f_QTMlevel )
        call expr ( subtree(2,son), expr_units, expr_value )
          QTM_level = nint(expr_value(1))
        select case ( expr_units(1) )
        case ( phyq_angle ) ! expr_value is in degrees
          QTM_level = ceiling ( log ( 180.0d0 / expr_value(1) ) / ln2 )
        case ( phyq_dimensionless )
          QTM_level = nint(expr_value(1))
        case ( phyq_length ) ! expr_value is in meters
          QTM_level = ceiling ( log (20.0d6 / expr_value(1)) / ln2 )
        case default ! Don't really need this because units are checked already
          call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
            & // 'CreateHGridFromMLSCFInfo', &
            & "Units of QTMLevel not dimensionless, length, or angle." )
        end select
        if ( QTM_level > QTM_Depth ) then
          call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
            & // 'CreateHGridFromMLSCFInfo', &
            & "QTM resolution is too fine; maximum depth is %D.", QTM_depth, &
            & where=where(son), advance='no' )
          call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
            & // 'CreateHGridFromMLSCFInfo', &
            & " QTM resolution is 20000/2**level km or 180/2**level degrees." )
        end if
      case ( f_single )
        single = get_boolean ( fieldValue )
      case ( f_solarTime )
        solarTimeNode = son
      case ( f_solarZenith )
        solarZenithNode = son
      case ( f_sourceL2gp )
        l2gp => l2gpDatabase(decoration(decoration(subtree(2,son))))
      case ( f_spacing )
        call expr ( subtree(2,son), expr_units, expr_value )
        spacing = expr_value(1)
      case ( f_time )
        timeNode = son
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    hGrid%type = hGridType
    select case (hGridType)

    case ( l_height, l_fractional, l_fixed ) ! --Fixed, Fractional or Height --
      if (.not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      else if (.not. NEED_L1BFILES ) then
        call announce_error ( root, NoL1BFILES )
      else
        if ( verbose ) call output ( 'Creating MIF-Based HGrid', advance='yes' )
        ! 1st--did we set any explicit geolocations?
        if ( any ( got_field ( (/ f_geodAngle, f_geodLat, f_losAngle, f_lon, &
          & f_solarTime, f_solarZenith, f_time /) ) ) ) then
          call CreateExplicitHGrid ( son, date, geodAngleNode, geodLatNode, &
            & solarTimeNode, solarZenithNode, lonNode, losAngleNode, incline, &
            & processingRange%startTime, timeNode, hGrid )
        end if
        ! call Dump ( filedatabase, details=2 )
        if ( mycheck ) then
          call outputNamedValue( 'Before creating mif-based MGrid', chunk%chunkNumber )
          call CheckForCorruptFileDatabase( filedatabase )
        endif
        call CreateMIFBasedHGrids ( filedatabase, hGridType, chunk, &
          & got_field, root, height, fraction, interpolationFactor, &
          & instrumentModuleName, mif, maxLowerOverlap, maxUpperOverlap, hGrid )
        if ( mycheck ) then
          call outputNamedValue( 'After creating mif-based MGrid', chunk%chunkNumber )
          call CheckForCorruptFileDatabase( filedatabase )
        endif
      end if

    case ( l_explicit ) ! ---------------------- Explicit --------------
      if ( verbose ) call output ( 'Creating Explicit HGrid', advance='yes' )
      call CreateExplicitHGrid ( son, date, geodAngleNode, geodLatNode, &
        & solarTimeNode, solarZenithNode, lonNode, losAngleNode, incline, &
        & processingRange%startTime, timeNode, hGrid )

    case ( l_l2gp) ! --------------------------- L2GP ------------------
      
      if ( verbose ) call output ( 'Creating L2GP-Based HGrid', advance='yes' )
      ! Get the time from the l1b file
      call L1BGeoLocation ( filedatabase, "MAFStartTimeTAI", &
        & instrumentModuleName, fullArray )
      call L1BSubsample ( chunk, FullArray, values=TAI )
      noMAFs = size(TAI)
      minTime = TAI(1) ! l1bField%dpField(1,1,1)
      maxTime = TAI(noMAFs) ! l1bField%dpField(1,1,noMAFs)

      call Hunt ( l2gp%time, (/minTime, maxTime/), profRange, allowTopValue=.true. )
      profRange(1) = min ( profRange(1)+1, l2gp%nTimes )
      
      a = profRange(1)
      b = profRange(2)
      hGrid%noProfs = b - a + 1
      call CreateEmptyHGrid(hGrid)
      hGrid%phi(1,:) =         l2gp%geodAngle(a:b)
      hGrid%geodLat(1,:) =     l2gp%latitude(a:b)
      hGrid%lon(1,:) =         l2gp%longitude(a:b) 
      hGrid%time(1,:) =        l2gp%time(a:b)
      hGrid%solarTime(1,:) =   l2gp%solarTime(a:b)
      hGrid%solarZenith(1,:) = l2gp%solarZenith(a:b)
      hGrid%losAngle(1,:) =    l2gp%losAngle(a:b)

      call Deallocate_test ( fullArray, 'fullArray', ModuleName )
      call Deallocate_test ( TAI, 'TAI', ModuleName )
      ! call deallocateL1BData ( l1bField ) ! Avoid memory leaks

    case ( l_QTM ) ! --------------------------- Regular ---------------
      if ( verbose ) call output ( 'Creating QTM HGrid', advance='yes' )
      if ( .not. allocated(polygon_vertices) ) then
        call announce_error ( root, noPolygon )
      else
        call output('insrumentModuleName is ' )
        call output(instrumentModuleName, advance='yes')
        call create_QTM_hgrid ( filedatabase, polygon_inside, polygon_vertices, QTM_level, &
          & hGrid, trim(instrumentModuleName) )
      end if

    case ( l_regular ) ! ----------------------- Regular ---------------
      if ( verbose ) call output ( 'Creating Regular HGrid', advance='yes' )
      if ( .not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      else if ( .not. all(got_field((/f_spacing, f_origin/)))) then
        call announce_error ( root, NoSpacingOrigin )
      else if (.not. NEED_L1BFILES ) then
        call announce_error ( root, NoL1BFILES )
      else
        call CreateRegularHGrid ( filedatabase, processingRange, chunk, &
          & spacing, origin, trim(instrumentModuleName), extendible, &
          & maxLowerOverlap, maxUpperOverlap, insetOverlaps, single, hGrid, &
          & onlyComputingOffsets )
      end if

    end select
    
    if ( mycheck ) then
      call outputNamedValue( 'Before finding nearest maf', chunk%chunkNumber )
      call CheckForCorruptFileDatabase( filedatabase )
    end if
    if ( hGridType /= l_QTM ) then
      ! Find nearest maf based on time
      if ( L1BFile%hdfVersion /= HDFVERSION_5 ) then
      elseif ( .not. IsHDF5DSInFile ( L1BFile%name, "MAFStartTimeTAI" ) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
          & // 'CreateHGridFromMLSCFInfo', &
          & "MAFStartTimeTAI not found in l1b file" )
        call trace_end ( "CreateHGridFromMLSCFInfo", &
          & cond=toggle(gen) .and. ( levels(gen) > 1 .or. .not. computingOffsets ) )
        return
      endif
      if ( verbose ) call output ( 'Reading MAFStartTimeTAI', advance='yes' )
      call L1BGeoLocation ( filedatabase, "MAFStartTimeTAI", &
          & instrumentModuleName, fullArray )
      if ( verbose ) call output ( 'Subsampling MAFStartTimeTAI', advance='yes' )
      call L1BSubsample ( chunk, fullArray, values=TAI )
      if ( verbose ) call output ( 'Hunting MAFStartTimeTAI', advance='yes' )
      call Hunt ( TAI, hgrid%time(1,:), hgrid%maf, &
        & allowTopValue=.true. )
    end if

    if ( DeeBug ) then
      call outputNamedValue ( 'firstMAF', chunk%firstMAFIndex )
      call outputNamedValue ( 'lastMAF',  chunk%lastMAFIndex )
      call outputNamedValue ( 'shape(MAFStartTimeTAI)', &
        &                      shape(TAI) )
      call outputNamedValue ( 'MAFStartTimeTAI (hid)',  &
        & tai93s2hid(TAI, leapsec=.true.) )
      call outputNamedValue ( 'hgrid time (hid)', &
        & tai93s2hid(hgrid%time(1,:), leapsec=.true.) )
      call outputNamedValue ( 'nearest maf', hgrid%maf + chunk%firstMAFIndex )
      call outputNamedValue ( 'nearest maf time', &
        & tai93s2hid(TAI(hgrid%maf), leapsec=.true.) )
    end if
    ! No--we already use hgrid%maf as a relativee maf numbr elsewhere in l2
    ! hgrid%maf = hgrid%maf + chunk%firstMAFIndex - 1
    ! call deallocateL1BData ( l1bField )
    if ( hGridType /= l_QTM ) then
      call Deallocate_test ( fullArray, 'fullArray', ModuleName )
      call Deallocate_test ( TAI, 'TAI', ModuleName )
      if ( switchDetail(switches, 'geom') >= 0 .and. .not. mySuppressGeometryDump ) &
        & call DumpChunkHGridGeometry ( hGrid, chunk, &
        & trim(instrumentModuleName), filedatabase )
    end if
    if ( verbose ) then
      call sayTime ( 'Constructing this HGrid' )
      call time_now( t1 )
    end if

    if ( error /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
      & // 'CreateHGridFromMLSCFInfo', &
      & "See ***** above for error message" )
    call trace_end ( "CreateHGridFromMLSCFInfo", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )
  end function CreateHGridFromMLSCFInfo

  ! ----------------------------------------  CreateExplicitHGrid  -----
  subroutine CreateExplicitHGrid ( key, date, geodAngleNode, geodLatNode, &
        & solarTimeNode, solarZenithNode, lonNode, losAngleNode, incline, &
        & Time, timeNode, hGrid )

    use Dates_Module, only: UTC2TAI93s
    use Error_Handler, only: Simple, Error_Intro
    use Expr_M, only: Expr
    use HighOutput, only: BeVerbose, OutputNamedValue
    use Geometry, only: Phi_To_Lat_Deg
    use Global_Settings, only: LeapSecFilename
    use HgridsDatabase, only: CreateEmptyHGrid, HGrid_T
    use MLSKinds, only: Rk => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use SDPToolkit, only: MLS_UTCtoTAI
    use String_Table, only: Get_String
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: NSons, Source_Ref, Subtree

    ! dummy arguments
    integer, intent(in) :: KEY           ! Tree node
    integer, intent(in) :: DATE          ! Date if any
    integer, intent(in) :: GEODANGLENODE ! Geod angle if any
    integer, intent(in) :: GEODLATNODE   ! Geod latitude if any
    integer, intent(in) :: LONNODE       ! longitude if any
    integer, intent(in) :: LOSANGLENODE  ! los angle if any
    real(rk), intent(in) :: Incline      ! Orbit inclination, degrees
    integer, intent(in) :: SOLARTIMENODE ! Solar time if any
    integer, intent(in) :: SOLARZENITHNODE ! Solar zenith angle if any
    integer, intent(in) :: TIMENODE ! Explicit d.p. numeric time field if any
    real(rk), intent(in) :: TIME        ! Time for HGrid
    type (HGrid_T), intent(inout) :: HGRID
    ! Local variables
    character(len=80) :: DATESTRING
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    real(rk) :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: NOPROFS                  ! Number of profiles
    integer :: NODE                     ! A tree node
    integer :: Me = -1                  ! String index for tracing
    integer :: PARAM                    ! Loop counter
    integer :: PROF                     ! Loop counter
    integer :: RETURNSTATUS             ! Flag
    real(rk), dimension(:,:), pointer :: VALUES
    logical :: verbose
    ! The following mysterious integers allow us
    ! to engage in some contemptible trickery to fill
    ! the HGrid's geolocation fields--a better
    ! method would call a well-written procedure, passing it appropriate
    ! args for each field
    ! However, we inherited this piece-o-work and just
    ! expanded it include three more fields
    integer, parameter :: PHIPARAM = 1
    integer, parameter :: SOLARTIMEPARAM = PHIPARAM + 1
    integer, parameter :: SOLARZENITHPARAM = SOLARTIMEPARAM + 1
    integer, parameter :: TIMEPARAM = SOLARZENITHPARAM + 1
    integer, parameter :: LONPARAM = TIMEPARAM + 1
    integer, parameter :: LOSANGLEPARAM = LONPARAM + 1
    integer, parameter :: GEODLATPARAM = LOSANGLEPARAM + 1
    integer, parameter :: NUMPARAMS = GEODLATPARAM

    integer :: ParamNode(numParams)

    ! Executable code

    call trace_begin ( me, "CreateExplicitHGrid", key, &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )
    verbose = BeVerbose( 'hgrid', -1 )

    paramNode = [ &
      & geodAngleNode             , & ! PHIPARAM
      & solarTimeNode             , & ! SOLARTIMEPARAM
      & solarZenithNode           , & ! SOLARZENITHPARAM
      & TimeNode                  , & ! TIMEPARAM
      & lonNode                   , & ! LONPARAM
      & losAngleNode              , & ! LOSANGLEPARAM
      & geodLatNode               ]   ! GEODLATPARAM

    ! 1st--try to get the number of instances by rude count
    noProfs = 0
    do param = 1, NUMPARAMS
      node = paramNode(param)
      if ( node /= 0 ) then
        if ( noProfs /= 0 ) then
          if ( nSons(node) - 1 /= noProfs ) then
            call error_intro ( simple, source_ref(node) )
            call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Explicit HGrid fields have different sizes' )
          end if
        end if
        noProfs = nSons ( node ) - 1
      end if
    end do

    ! Create an hGrid with this many profiles and no overlaps
    hGrid%noProfs = noProfs
    hGrid%noProfsLowerOverlap = 0
    hGrid%noProfsUpperOverlap = 0

    call CreateEmptyHGrid(hGrid)
    ! Fill up the obvious stuff
    hGrid%lon = 0.0_rk
    hGrid%losAngle = 0.0_rk
    if ( date /= 0 ) then
      call get_string ( date, dateString, strip=.true. )
      if ( verbose ) call outputNamedValue( 'dateString', dateString )
      if ( verbose ) call outputNamedValue( 'LeapSecFileName', LeapSecFileName )
      if ( len_trim(LeapSecFileName) > 1 ) then
        returnStatus = mls_utctotai &
          & ( trim(LeapSecFileName), trim(dateString), hGrid%time(1,1) )
      else
        ! Without a leapsec file, let's use date_module's built-in feature
        hGrid%time(1,1) = utc2tai93s ( dateString, leapsec=.true. )
      end if
      if ( returnStatus /= 0 ) call announce_error( key, badTime )
      hGrid%time = hGrid%time(1,1)
    else
      hGrid%time = time                   ! Get it from input time range
    end if
    ! Set defaults for the others, our expressions may overwrite them
    hGrid%phi = 0.0_rk
    hGrid%geodLat = 0.0_rk
    hGrid%solarTime = 0.0_rk
    hGrid%solarZenith = 0.0_rk
    hGrid%lon = 0.0_rk
    hGrid%losAngle = 0.0_rk

    ! Loop over the parameters we might have
    ! which may override any set above
    do param = 1, NUMPARAMS
      node = paramNode(param)
      if ( node /= 0 ) then
        select case ( param )
        case ( PHIPARAM )
          values => hGrid%phi
        case ( SOLARTIMEPARAM )
          values => hGrid%solarTime
        case ( SOLARZENITHPARAM )
          values => hGrid%solarZenith
        case ( TIMEPARAM )
          values => hGrid%Time
        case ( LONPARAM )
          values => hGrid%lon
        case ( LOSANGLEPARAM )
          values => hGrid%losAngle
        case ( GEODLATPARAM )
          values => hGrid%GeodLat
        end select
        do prof = 1, hGrid%noProfs
          call expr ( subtree ( prof+1, node), expr_units, expr_value )
          ! Units are checked in the type checker.
          values(1,prof) = expr_value(1)
        end do
      end if
    end do

    ! If the geodetic latitude isn't explicitly specified, compute it from
    ! orbit geodetic angle and orbit inclination.
    if ( geodLatNode == 0 ) &
      & hGrid%geodLat = Phi_To_Lat_Deg ( hGrid%phi, incline )

    call trace_end ( "CreateExplicitHGrid", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

  end subroutine CreateExplicitHGrid

  ! ---------------------------------------  CreateMIFBasedHGrids  -----
  subroutine CreateMIFBasedHGrids ( filedatabase, hGridType, &
    & chunk, got_field, root, height, fraction, interpolationFactor,&
    & instrumentModuleName, mif, maxLowerOverlap, maxUpperOverlap, hGrid )
    ! This is part of ConstructHGridFromMLSCFInfo

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Chunks_M, only: MLSCHunk_T
    use Dump_0, only: Dump
    use HGridsDatabase, only: CreateEmptyHGrid, Dump, HGrid_T, TrimHGrid
    use HighOutput, only: LetsDebug, OutputNamedValue
    use Init_Tables_Module, only: F_Fraction, F_Geodangle, F_Geodlat, F_Height, &
      & F_Lon, F_Losangle, F_Mif, F_Time, &
      & F_Solartime, F_Solarzenith, L_Fixed, L_Fractional, L_Height, L_Mif
    use L1BData, only: DeallocateL1BData, L1BData_T, ReadL1BData, &
      & AssembleL1BQtyName, CheckForCorruptFileDatabase
    use MLSHDF5, only: IsHDF5DSInFile
    use MLSKinds, only: Rk => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_L1BRead, &
      & MLSMSG_Warning
    use MLSNumerics, only: Hunt, InterpolateValues
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer, intent(in)               :: HGRIDTYPE
    type (MLSChunk_T), intent(in)     :: CHUNK
    integer                           :: error
    logical, dimension(:), intent(in) :: GOT_FIELD
    integer, intent(in)               :: ROOT
    real(rk), intent(in)              :: INTERPOLATIONFACTOR
    real(rk), intent(in)              :: HEIGHT
    real(rk), intent(in)              :: FRACTION
    character (len=*), intent(in)     :: INSTRUMENTMODULENAME
    integer, intent(in)               :: MIF
    integer, intent(in)               :: MAXLOWEROVERLAP
    integer, intent(in)               :: MAXUPPEROVERLAP
    type (HGrid_T), intent(inout)     :: HGRID ! Needs inout as name set by caller

    ! Local variables / parameters
    integer, parameter :: L1B_MAFSTARTTIMETAI = 1
    integer, parameter :: L1B_TPGEODLAT       = l1b_MAFStartTimeTAI + 1
    integer, parameter :: L1B_TPLON           = l1b_tpGeodLat + 1
    integer, parameter :: L1B_TPGEODANGLE     = l1b_tpLon + 1
    integer, parameter :: L1B_TPSOLARZENITH   = l1b_tpGeodAngle + 1
    integer, parameter :: L1B_TPSOLARTIME     = l1b_tpSolarZenith + 1
    integer, parameter :: L1B_TPLOSANGLE      = l1b_tpSolarTime + 1
    integer, parameter :: NOL1BITEMSTOREAD=l1b_tpLosAngle
    integer, parameter :: FIRSTMODULARITEM=l1b_tpGeodLat

    real(rk), parameter :: SIXTH = 1.0_rk / 6.0_rk

    character (len=15), dimension(noL1BItemsToRead) :: L1bItemNames
    ! Entries in the above array following FirstModularItem are prefixed
    ! with either GHz or THz. 

    logical :: deeBug
    integer :: L1BFLAG                  ! Flag
    integer :: L1BITEM                  ! Loop counter
    integer :: MAF                      ! Loop counter etc.
    integer :: Me = -1                  ! String index for trace
    real(rk) :: MinAngle, MaxAngle
    logical :: MissingOK
    logical :: MissingThisOK
    logical :: myCheck
    integer :: NOMAFS                   ! Dimension
    real(rk), dimension(:,:,:), pointer :: TpGeodAlt, TpGeodAngle

    ! MIFs it would choose in the non over/undersampled case
    real(rk), dimension(:), pointer :: defaultField, interpolatedField

    type (L1BData_T) :: l1bField ! L1B data
    integer, dimension(:), pointer :: defaultMIFs
    real(rk), dimension(:), pointer :: defaultIndex
    real(rk), dimension(:), pointer :: interpolatedIndex
    integer ::  hdfVersion
    ! integer ::  status
    character (len=NameLen) :: L1BItemName
    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code
    deebug = LetsDebug ( 'hgrid', 1 )
    call trace_begin ( me, "CreateMIFBasedHGrids", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    hdfversion = L1BFile%HDFVersion
    myCheck = .false. ! (chunk%chunkNumber == 1)
    if ( mycheck ) then
      call outputNamedValue( 'On entering mif-based MGrid', chunk%chunkNumber )
      call CheckForCorruptFileDatabase( filedatabase )
    endif
    ! call MLS_closeFile ( L1BFile, status )
    ! call outputNamedValue ( 'Status on close', status )
    ! call Dump ( L1BFile, details=2 )
    ! call MLS_openFile ( L1BFile, status )
    ! call outputNamedValue ( 'Status on open', status )
    ! call Dump ( L1BFile, details=2 )
    nullify ( tpGeodAlt, tpGeodAngle, defaultField, interpolatedField, &
      & defaultMIFs, defaultIndex, interpolatedIndex )

    l1bItemNames(l1b_mafstarttimetai ) = 'MAFStartTimeTAI'
    l1bItemNames(l1b_tpgeodlat       ) = 'tpGeodLat'
    l1bItemNames(l1b_tplon           ) = 'tpLon'
    l1bItemNames(l1b_tpgeodangle     ) = 'tpGeodAngle'
    l1bItemNames(l1b_tpsolarzenith   ) = 'tpSolarZenith'
    l1bItemNames(l1b_tpsolartime     ) = 'tpSolarTime'
    l1bItemNames(l1b_tplosangle      ) = 'tpLosAngle'

    select case ( hGridType )
    case ( l_Fractional )
      if ( .not. got_field(f_Fraction) ) &
        & call announce_error ( root, noFraction )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAngle"
    case ( l_Height )
      if ( .not. got_field(f_height) ) call announce_error ( root, noHeight )
      l1bItemName = trim(instrumentModuleName)//"."//"tpGeodAlt"
    case ( l_Fixed )
      if ( .not. got_field(f_mif) ) call announce_error ( root, noMIF )
    case ( l_mif )
      ! ??? Does something go here?
      ! ??? Probably not, since L_MIF isn't in hGridType in init_tables_module
    end select

    if ( hGridType /= l_Fixed ) then
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex, &
        & dontPad=.true. )
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
    else
      noMAFs = chunk%lastMafIndex - chunk%firstMafIndex + 1
    end if

    ! allocate default MIFs
    call Allocate_test ( defaultMIFs, noMAFs, 'defaultMIFs', ModuleName )
    
    ! Work out which MIF should have the profile for each MAF.
    select case ( hGridType )
    case ( l_Fractional )
      ! A fractional hGrid, we need to read the tangent point phi
      tpGeodAngle=> l1bField%dpField
      
      ! Loop over the MAFs
      do maf = 1, noMAFs
        ! ??? Think about missing data here! ***
        ! Probably need to do a pack on tpGeodAngle and then unpack on
        ! defaultMIFs
        
        minAngle=minval(tpGeodAngle(1,:,maf))
        maxAngle=maxval(tpGeodAngle(1,:,maf))
        
        call Hunt ( tpGeodAngle(1,:,maf), &
          & minAngle+fraction*(maxAngle-minAngle) , defaultMIFs(maf) )
      end do
    case ( l_Height )
      tpGeodAlt => l1bField%dpField
      
      ! Loop over the MAFs
      do maf = 1, noMAFs
        ! ??? Think about missing data here! ***
        ! Probably need to do a pack on tpGeodAngle and then unpack on
        ! defaultMIFs
        
        call Hunt ( tpGeodAlt(1,:,maf), height, defaultMIFs(maf) )
      end do
    case ( l_fixed )
      defaultMIFs = mif
    end select
    
    if ( mycheck ) then
      call outputNamedValue( 'After reading 1st l1bdata', chunk%chunkNumber )
      call CheckForCorruptFileDatabase( filedatabase )
    endif
    ! Done with this piece of l1b data for the moment
    call DeallocateL1BData ( l1bField )
    
    ! Now we have a default MIFs array; this is a list of the `standard'
    ! MIFs we would choose in the interpolationFactor=1 case.
    ! Work out how many profiles this is going to be.
    
    if ( .not. any( got_field ( (/ f_geodAngle, f_geodLat, f_losAngle, f_lon, &
      & f_solarTime, f_solarZenith, f_time /) ) ) ) then
      ! Create an empty hGrid
      hGrid%noProfs = NINT(noMAFs*interpolationFactor)
      call CreateEmptyHGrid( hGrid )
      MissingOK = .false.
    else
      MissingOK = .true.
      ! call dump( hGrid )
    end if
    
    hGrid%noProfsLowerOverlap = 0
    hGrid%noProfsUpperOverlap = 0
    
    ! Setup some arrays
    call allocate_test ( defaultField, noMAFs, 'defaultField', ModuleName )
    call allocate_test ( interpolatedField, hGrid%noProfs, 'interpolatedField', ModuleName )
    call allocate_test ( defaultIndex, noMAFs, 'defaultIndex', ModuleName )
    call allocate_test ( interpolatedIndex, hGrid%noProfs, 'interpolatedIndex', ModuleName )

    ! Now we go through all the important geolocation quantities, read them
    ! in, interpolate them if required and store the result in the hGrid
    error = 0
    do l1bItem = 1, NoL1BItemsToRead
      ! Get the name of the item to read
      l1bItemName = l1bItemNames(l1bItem)
      MissingThisOK = MissingOK .and. l1bItem /= l1b_mafstarttimetai
      if ( l1bItem >= firstModularItem ) l1bItemName = &
        & trim(instrumentModuleName)//"."//l1bItemName
      
      ! Read it from the l1boa file
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
      if ( deebug ) then
        call outputNamedValue( 'item', trim(l1bItemName) )
        call outputNamedValue( 'MissingOK', MissingOK )
      endif
      if ( .not. IsHDF5DSInFile( L1BFile%name, l1bItemName ) &
        & .and. MissingthisOK ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        cycle
      endif
      if ( mycheck ) then
        call outputNamedValue( 'Before next reading l1bdata', trim(l1bItemName) )
        call CheckForCorruptFileDatabase( filedatabase )
      endif
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex, &
        & neverfail=(MissingOK .or. l1bItem == l1b_tplosangle), dontPad=.true. )
      if ( mycheck ) then
        call outputNamedValue( 'After next reading l1bdata', trim(l1bItemName) )
        call CheckForCorruptFileDatabase( filedatabase )
      endif
      if ( deebug ) call outputNamedValue( 'l1bFlag', l1bFlag )
      if ( l1bFlag==-1 .and. .not. MissingthisOK ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & MLSMSG_L1BRead//l1bItemName )
        error = 1
        stop
        cycle
      elseif ( .not. IsHDF5DSInFile ( L1BFile%name, trim(l1bItemName) ) ) then
        call outputNamedValue ( 'item missing from l1b file', trim(l1bItemName) )
        cycle
      endif
      if ( deebug .and. index(l1bItemName, 'MAFStartTimeTAI') > 0 ) then
        call dump(l1bField%DpField(1,1,:), 'MAFStartTimeTAI (before interpolating)')
      end if
      
      if ( l1bItem == 1 ) then       ! do something special for time
        do maf = 1, noMAFs
          defaultField(maf) = l1bField%dpField(1,1,maf) + &
            & (defaultMIFs(maf)-1)*sixth
        end do
      else                         ! Otherwise this is fairly "easy". ???
        if ( associated(l1bField%intField) ) then
          do maf = 1, noMAFs
            defaultField(maf) = l1bField%intField(1,defaultMIFs(maf),maf)
          end do
        else
          do maf = 1, noMAFs
            defaultField(maf) = l1bField%dpField(1,defaultMIFs(maf),maf)
          end do
        endif
      end if

      if ( mycheck ) then
        call outputNamedValue( 'After setting default field', trim(l1bItemName) )
        call CheckForCorruptFileDatabase( filedatabase )
      endif
      if ( deebug .and. index(l1bItemName, 'MAFStartTimeTAI') > 0 ) then
        call dump(defaultField, 'MAFStartTimeTAI, corrected to defaultMIFs')
      end if
      call DeallocateL1BData(l1bField)
      
      if ( interpolationFactor == 1.0 ) then
        interpolatedField = defaultField
      else
        ! Interpolations are performed based on MAF index and periodic variables
        ! like lon, solartime, losAngle are handled specially so that it can cope
        ! with the period jump between two adjacent MAFs. The range of the periods
        ! of these variable must be specified, i.e., (lower, upper) bounds. 

        ! create index of the old and new hGrid indices
        do maf=1,noMAFs
          defaultIndex(maf)=maf*1.0_rk 
        end do
        do maf=1,hGrid%noProfs 
          interpolatedIndex(maf)=maf/interpolationFactor
        end do

        select case ( l1bItem )
        case ( l1b_tpLon )
          call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
                interpolatedField,method='Linear',rangeofPeriod=(/-180.0_rk,180.0_rk/))
        case ( l1b_tpSolarTime )
          call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
                interpolatedField,method='Linear',rangeofPeriod=(/0.0_rk,24.0_rk/))
        case ( l1b_tpLosAngle )
          call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
                interpolatedField,method='Linear',rangeofPeriod=(/0.0_rk,360.0_rk/))
        case default
          call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
                interpolatedField,method='Linear')
        end select

      end if
      
      if ( mycheck ) then
        call outputNamedValue( 'After interpolating', trim(l1bItemName) )
        call CheckForCorruptFileDatabase( filedatabase )
      endif
      if ( deebug .and. index(l1bItemName, 'MAFStartTimeTAI') > 0 ) then
        call dump(InterpolatedField, 'MAFStartTimeTAI (after interpolating)')
      end if
      select case ( l1bItem )
      case ( l1b_MAFStartTimeTAI )
        hGrid%time(1,:) = interpolatedField
      case ( l1b_tpGeodLat )
        hGrid%geodLat(1,:) = interpolatedField
      case ( l1b_tpLon )
        hGrid%lon(1,:) = interpolatedField
      case ( l1b_tpGeodAngle )
        hGrid%phi(1,:) = interpolatedField
      case ( l1b_tpSolarZenith )
        hGrid%solarZenith(1,:) = interpolatedField
      case ( l1b_tpSolarTime )
        hGrid%solarTime(1,:) = interpolatedField
      case ( l1b_tpLosAngle )
        hGrid%losAngle(1,:) = interpolatedField
      end select
      if ( mycheck ) then
        call outputNamedValue( 'After filling HGrid with this l1bdata', chunk%chunkNumber )
        call CheckForCorruptFileDatabase( filedatabase )
      endif
    end do
    
    call Deallocate_test ( defaultMIFs, 'defaultMIFs', ModuleName )
    call Deallocate_test ( defaultField, 'defaultField', ModuleName )
    call Deallocate_test ( defaultIndex, 'defaultIndex', ModuleName )
    call Deallocate_test ( interpolatedField, 'interpolatedField', ModuleName )
    call Deallocate_test ( interpolatedIndex, 'interpolatedIndex', ModuleName )
    
    ! ??? This calculation may need attention! ***
    if ( deeBug ) call outputNamedValue( 'interpolationFactor', interpolationFactor )
    hGrid%noProfsLowerOverlap = &
      & nint(chunk%noMAFsLowerOverlap*interpolationFactor)
    hGrid%noProfsUpperOverlap = &
      & nint(chunk%noMAFsUpperOverlap*interpolationFactor)
    if ( deeBug ) call outputNamedValue( 'overlaps', &
      & (/ hGrid%noProfsLowerOverlap, hGrid%noProfsUpperOverlap/) )

    ! Now deal with the user requests
    if ( maxLowerOverlap >= 0 .and. &
      & ( hGrid%noProfsLowerOverlap > maxLowerOverlap ) .and. error == 0 ) &
      call TrimHGrid ( hGrid, -1, hGrid%noProfsLowerOverlap - maxLowerOverlap )
    if ( maxUpperOverlap >= 0 .and. &
      & ( hGrid%noProfsUpperOverlap > maxUpperOverlap ) .and. error == 0 ) &
      call TrimHGrid ( hGrid, -1, hGrid%noProfsUpperOverlap - maxUpperOverlap )

    if ( mycheck ) then
      call outputNamedValue( 'After trimming hgrid', chunk%chunkNumber )
      call CheckForCorruptFileDatabase( filedatabase )
    endif
    call trace_end ( "CreateMIFBasedHGrids", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

  end subroutine CreateMIFBasedHGrids

  ! -------------------------------------------  Create_QTM_HGrid  -----
  subroutine Create_QTM_HGrid ( filedatabase, Polygon_Inside, Polygon_Vertices, Level, HGrid, instrumentModuleName )

    ! Create a QTM within the specified polygon with the specified level
    ! of refinement.

    ! Presumably we don't get here until it's been checked whether the
    ! polygon exists.

    ! We need to know both the polygon boundary and a point inside the
    ! polygon because the concept "inside a polygon" is ambiguous on the
    ! surface of a sphere.

    use Allocate_Deallocate, only: Test_Allocate
    use Generate_QTM_M, only: Generate_QTM
    use Geolocation_0, only: H_T
    use HGridsDatabase, only: CreateEmptyHGrid
    use MLSStringLists, only: SwitchDetail
    use QTM_Output, only: Write_QTM_Unformatted
    use Toggles, only: Switches
    use Output_M, only: Output
    use L1bData, only: L1bData_T, ReadL1BData, Namelen, &
      & AssembleL1Bqtyname, Dump
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Geometry, only: To_XYZ
    
    type (MLSFile_T), dimension(:), pointer ::     FileDatabase
    type(h_t), intent(in) :: Polygon_Inside ! A point inside the polygon
    type(h_t), intent(in) :: Polygon_Vertices(:)
    integer, intent(in) :: Level
    type(hGrid_t), intent(inout) :: HGrid
    character (len=*), intent(in) :: InstrumentModuleName

    integer :: I, QTMFile, Stat
    type (L1BData_T) :: lats            ! lats data
    type (L1BData_T) :: lons            ! lons data
    type (MLSFile_T), pointer             :: L1BFile
    integer :: L1BFLAG                  ! error Flag
    integer :: NoMAFs
    integer :: MAFsofar
    character (len=NameLen) :: L1BItemName
    real(r8), dimension(:, :), pointer :: xyz
    real(r8), dimension(:, :), pointer :: xyz_p
    integer :: MAF
    integer :: noP
    real(r8), dimension(:),  pointer :: r
    real(r8) :: rsofar
    type (L1BData_T) :: SolarTimes            ! SolarTime data
    type (L1BData_T) ::SolarZenith
    type (L1BData_T) ::MAFStartTimeTAI
     
    
    call output ('Inside Create_QTM_HGrid', advance='yes' )
    ! 1st-- read lats and lons from the l1b file
    
    L1BFile => GetMLSFileByType( filedatabase, content='l1boa' )  
    if ( .not. associated(L1BFile) ) &
          & call MLSMessage  ( MLSMSG_Error, ModuleName, &
          & 'Can not make progress in  Create_QTM_HGrid without L1BOA files' )
    l1bItemName = AssembleL1BQtyName (  instrumentModuleName//".tpGeodLat", L1BFile%HDFVersion, .false. )
    call output (trim(l1bItemName), advance ='yes')
    call ReadL1BData ( L1BFile, l1bItemName, lats, noMAFs, &
        & l1bFlag )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Create_QTM_HGrid Where is ' // l1bItemName )
    call output ('lats is ', advance='yes' )
    call output( lats%dpField(1,1,1::NoMAFS), advance='yes')
    call dump ( lats, 0)
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpLon", L1BFile%HDFVersion, .false. )
    call output (trim(l1bItemName), advance ='yes')
    call ReadL1BData ( L1BFile, l1bItemName, lons, noMAFs, &
        & l1bFlag )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Create_QTM_HGrid Where is ' // l1bItemName )
    call output ('lons is ' )
    call output( lons%dpField(1,1,1::NoMAFS), advance='yes')
    call dump ( lons, 0 ) 
    
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpSolarTime", L1BFile%HDFVersion, .false. )
    call output (trim(l1bItemName), advance ='yes')
    call ReadL1BData ( L1BFile, l1bItemName, SolarTimes, noMAFs, &
        & l1bFlag )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Create_QTM_HGrid Where is ' // l1bItemName )
    call output ('SolarTimes is ' )
    call output( SolarTimes%dpField(1,1,1::NoMAFS), advance='yes')
    call dump ( SolarTimes, 0 ) 
    
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpSolarZenith", L1BFile%HDFVersion, .false. )
    call output (trim(l1bItemName), advance ='yes')
    call ReadL1BData ( L1BFile, l1bItemName, SolarZenith, noMAFs, &
        & l1bFlag )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Create_QTM_HGrid Where is ' // l1bItemName )
    call output ('SolarZenith is ' )
    call output( SolarZenith%dpField(1,1,1::NoMAFS), advance='yes')
    call dump ( SolarZenith, 0 ) 
   
    l1bItemName = AssembleL1BQtyName ( "MAFStartTimeTAI", L1BFile%HDFVersion, .false. )
    call output (trim(l1bItemName), advance ='yes')
    call ReadL1BData ( L1BFile, l1bItemName, MAFStartTimeTAI, noMAFs, &
        & l1bFlag )
    if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Create_QTM_HGrid Where is ' // l1bItemName )
    call output ('MAFStartTimeTAI is ' )
    call output( MAFStartTimeTAI%dpField(1,1,1::NoMAFS), advance='yes')
    call dump ( MAFStartTimeTAI, 0 )
  
    !Converting lats, lons to XYZ
    call output ('Create_QTM_HGrid before do MAF loop', advance='yes')
    allocate (xyz(3,NoMAFs), stat = stat)
    call test_allocate (stat, moduleName, "xyz", &
        & [3], [NoMAFS], storage_size(3*NoMAFs ) )
    do MAF=1, NoMAFs
       call output('lat = ')
       call output( lats%dpField(1,1,MAF) )
       call output ('  lon = ')   
       call output( lons%dpField(1,1,MAF), advance='yes')
       xyz(1:3,MAF) = to_xyz ( lats%dpField(1,1,MAF), lons%dpField(1,1,MAF))
       call output ('xyz is ')
       call output (xyz(1:3, MAF), advance='yes')
    enddo
    
    call output('Testing xyz array after population', advance='yes')
    do MAF = 1, NoMAFs
       call output ('MAF = ')
       call output (MAF)
       call output (' , XYZ is ')
       call output ( xyz(1:3, MAF), advance='yes')
    enddo
               
    allocate ( hGrid%QTM_tree, stat=stat )
    call test_allocate ( stat, moduleName, "hGrid%QTM_tree", &
      & [1], [1], storage_size(hGrid%QTM_tree) / 8 )
    hGrid%QTM_tree%level = level
    hGrid%QTM_tree%in_geo = polygon_inside
    ! Explicit allocation won't be necessary when compilers support
    ! automatic allocation for assignment to allocatable arrays.
    allocate ( hGrid%QTM_tree%polygon_geo(size(polygon_vertices)), stat=stat )
    call test_allocate ( stat, moduleName, "hGrid%Polygon_geo", &
      & [1], [size(polygon_vertices)], storage_size(polygon_vertices) / 8 )
    hGrid%QTM_tree%polygon_geo = polygon_vertices

    call generate_QTM ( hGrid%QTM_tree )

    hGrid%noProfs = hGrid%QTM_tree%n_in
    call CreateEmptyHGrid(hGrid)

    ! Put imaginative values in HGrid fields so they're not undefined.
    ! They might actually be useful.  Some might get more meaningful
    ! values elsewhere.

    hGrid%maf = [ ( i, i = 1, hGrid%noProfs ) ]
    hGrid%phi(1,:) = 0 ! Phi is kind of meaningless for QTM
    hGrid%geodLat(1,:) = hGrid%QTM_tree%geo_in%lat
    hGrid%lon(1,:) = hGrid%QTM_tree%geo_in%lon%d
    hGrid%time(1,:) = 0
    hGrid%solarTime(1,:) = 0
    hGrid%solarZenith(1,:) = 0
    hGrid%losAngle(1,:) = 0

    allocate (xyz_p(3,hGrid%noProfs), stat = stat)
    call test_allocate (stat, moduleName, "xyz_p", &
        & [3], [hGrid%noProfs], storage_size(3*hGrid%noProfs ) )
     
    do noP=1, hGrid%noProfs
       call output('lat for profile = ')
       call output( hGrid%geodLat(1,noP) )
       call output ('  lon for profile  = ')   
       call output( hGrid%lon(1,noP), advance='yes')
       xyz_p(1:3,noP) = to_xyz ( hGrid%geodLat(1,noP), hGrid%lon(1,noP))
       call output ('xyz_p is ')
       call output (xyz_p(1:3,noP), advance='yes')   
    enddo 
     
    allocate (r(NoMAFs), stat = stat)
    call test_allocate (stat, moduleName, "r", &
        & [NoMAFs], [0], storage_size(NoMAFs) )
        
    call output('Ready to compute distance r', advance='yes')
    
    do  noP=1, hGrid%noProfs   
       r = (xyz(1,:) - xyz_p(1,noP))**2 + &
           &(xyz(2,:) - xyz_p(2,noP))**2 +  &
           &(xyz(3,:) - xyz_p(3,noP))**2
       call output('distance for ')
       call output(noP)
       call output ('  profile is ')
       call output (r, advance='yes') 
       rsofar = r(1)
       MAFsofar = 1
       do  MAF=2, noMAFS 
         if(r(MAF) .GE. rsofar) cycle
         rsofar = r(MAF)
         MAFsofar = MAF
       enddo        
       call output ('rsofar is ')
       call output (rsofar, advance = 'yes') 
       call output ('MAFsofar is ')
       call output (MAFsofar, advance = 'yes') 
       hGrid%time(1,noP) =  MAFStartTimeTAI%dpField(1,1,MAFsofar)
       hGrid%solarTime(1,noP) =  SolarTimes%dpField(1,1,MAFsofar)
       hGrid%solarZenith(1,noP) =  SolarZenith%dpField(1,1,MAFsofar)      
    enddo
       
    
    QTMFile = switchDetail ( switches, 'QTMFile' )
    if ( QTMFile > 0 ) then
      ! Write a Fortran unformatted file to make it easier to look at a QTM
      ! using other tools, e.g., IDL
      call write_QTM_unformatted ( hGrid%QTM_tree, QTMFile )
    end if

  end subroutine Create_QTM_HGrid

  ! -----------------------------------------  CreateRegularHGrid  -----
  subroutine CreateRegularHGrid ( filedatabase, processingRange, chunk, &
    & spacing, origin, instrumentModuleName, extendible, &
    & maxLowerOverlap, maxUpperOverlap, insetOverlaps, single, hGrid, &
    & onlyComputingOffsets )

    ! Creates an HGrid with coordinates laid out in
    ! a regular spacing w.r.t. master coordinate, phi, aka Geodetic Angle

    ! Depends on an approximation to the Earth's shape (Empirical Geometry)
    ! and interpolation
    
    ! With older l1b files (pre v2.0), some coordinates
    ! (solar time, solar zenith angle) are mean rather than apparent local
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use ChunkDivide_M, only: ChunkDivideConfig
    use Chunks_M, only: MLSChunk_T, Dump
    use Dump_0, only: Dump
    use Diff_1, only: Selfdiff
    use EmpiricalGeometry, only: EmpiricalLongitude, ChooseOptimumLon0
    use Geometry, only: Phi_To_Lat_Deg
    use HGridsDatabase, only: CreateEmptyHGrid, HGrid_T, TrimHGrid, &
      & FindClosestMatch
    use HighOutput, only: BeVerbose, LetsDebug, OutputNamedValue
    use MLSFillValues, only: IsFillValue, Monotonize
    use MLSKinds, only: Rk => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MLSNumerics, only: Hunt, InterpolateValues
    use MLSStringLists, only: SwitchDetail
    use Monotone, only: IsMonotonic
    use Output_M, only: Output
    use String_Table, only: Display_String
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End

    type (MLSFile_T), dimension(:), pointer ::     FileDatabase
    type (TAI93_Range_T), intent(in) :: ProcessingRange
    type (MLSChunk_T), intent(in) :: Chunk
    real(rk), intent(in) :: Spacing
    real(rk), intent(in) :: Origin
    character (len=*), intent(in) :: InstrumentModuleName
    logical, intent(in) :: Extendible
    integer, intent(in) :: MaxLowerOverlap
    integer, intent(in) :: MaxUpperOverlap
    logical, intent(in) :: InsetOverlaps
    logical, intent(in) :: Single
    type (HGrid_T), intent(inout) :: HGrid ! Needs inout as name set by caller
    logical, intent(in), optional :: onlyComputingOffsets

    ! Local variables/parameters
    real(rk), dimension(:,:), pointer :: AllGeodAngle ! For every MIF
    logical :: deeBug
    real(rk), parameter :: SecondsInDay = 24*60*60
    ! Note this next one is ONLY NEEDED for the case where we have only
    ! one MAF in the chunk
    real(rk), parameter :: OrbitalPeriod = 98.8418*60.0

    logical :: DeebugHere
    real(rk) :: Delta                   ! A change in angle
    integer :: Extra                    ! How many profiles over 1 are we
    real(rk) :: First                   ! First point in of hGrid
    integer :: FirstProfInRun           ! Index of first profile in processing time
    double precision, dimension(:), pointer :: FullArray
    double precision, dimension(:,:), pointer :: FullArray2d
    integer :: I                        ! Loop counter
    real(rk) :: Incline                 ! Mean orbital inclination
    type (MLSFile_T), pointer             :: L1BFile
    integer :: LastProfInRun            ! Index of last profile in processing time
    real(rk) :: Last                    ! Last point in hGrid
    integer :: Left                     ! How many profiles to delete from the LHS in single
    real(rk) :: MaxAngle                ! Largest angle in chunk
    real(rk) :: MaxAngleFirstMAF        ! Gives 'range' of first maf
    integer :: Me = -1                  ! String index for tracing
    real(rk), dimension(:), pointer :: MIF1GEODANGLE ! For first mif
    real(rk) :: MinAngle                ! Smallest angle in chunk
    real(rk) :: MinAngleLastMAF         ! Gives 'range' of last maf
    integer :: N                        ! Guess at number of profiles
    integer :: NoMAFs                   ! How many in chunk
    integer :: NoMIFs
    real(rk) :: NextAngle               ! First non ovl. MAF for next chunk
    integer :: Right                     ! How many profiles to delete from the RHS in single
    double precision, dimension(:), pointer :: GeodAngle
    double precision, dimension(:,:), pointer :: GeodAngle2d
    double precision, dimension(:), pointer :: LosAngle
    double precision, dimension(:), pointer :: scOrbIncl
    double precision, dimension(:), pointer :: SolarZenith
    double precision, dimension(:,:), pointer :: SolarZenith2d
    double precision, dimension(:), pointer :: TAI
    real(rk), dimension(:), pointer :: TmpAngle ! A temporary array for the single case
    logical :: Verbose
    logical :: WarnIfNoProfs

    ! Executable code
    deebug = LetsDebug ( 'hgrid', 1 )
    call trace_begin ( me, "CreateRegularHGrid", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )
    deebughere = deebug .or. ( switchDetail(switches, 'hgrid') > 0 ) ! e.g., 'hgrid1' 
    warnIfNoProfs = single ! .or. ChunkDivideConfig%maxLength < 2 ! .false.
    ! call outputNamedValue ( 'maxLength', ChunkDivideConfig%maxLength, options='--Banner' )
    ! call outputNamedValue ( 'warnIfNoProfs', warnIfNoProfs, options='--Banner' )
    verbose = BeVerbose ( 'hgrid', -1 ) ! deebughere
    if ( present(onlyComputingOffsets) ) then
      ! verbose = .not. onlyComputingOffsets
      warnIfNoProfs = warnIfNoProfs .or. onlyComputingOffsets
      deebughere = deebug .or. ( switchDetail(switches, 'hgrid') > 1 ) ! e.g., 'hgrid2' 
      ! verbose = verbose .or. deebughere
      verbose = BeVerbose ( 'hgrid', 1 ) ! deebughere
    end if
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file not in database' )

    ! Setup the empircal geometry estimate of lon0
    ! (it makes sure it's not done twice
    call ChooseOptimumLon0 ( filedatabase, chunk, instrumentModuleName )

    ! First we're going to work out the geodetic angle range
    ! Read an extra MAF if possible, as we may use it later
    ! when computing the overlaps.
    call L1BGeoLocation ( filedatabase, "tpGeodAngle", instrumentModuleName, &
      & fullArray, FullArray2d )
    call L1BSubsample ( chunk, fullArray, fullArray2d, GeodAngle, GeodAngle2d )
    if ( any(isFillValue(FullArray2d) ) ) then
      call output( 'Fill values among full array', advance='yes' )
    endif
    if ( any(isFillValue(GeodAngle2d) ) ) then
      call output( 'Fill values among reduced array', advance='yes' )
    endif
    noMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
    noMIFs = size(GeodAngle2d(:,1))
    if ( deebughere ) then
      call outputnamedValue ( 'noMAFs', noMAFs )
      call outputnamedValue ( 'noMIFs', noMIFs )
      call outputnamedValue ( 'shape(GeodAngle)', shape(GeodAngle) )
      call outputnamedValue ( 'shape(GeodAngle2d)', shape(GeodAngle2d) )
      call Dump ( chunk )
      call Dump ( GeodAngle2d, 'GeodAngle2d' )
    endif
    minAngle = minval ( GeodAngle2d(:,1) )
    maxAngleFirstMAF = maxval ( GeodAngle2d(:,1) )
    maxAngle = maxval ( GeodAngle2d(:,noMAFs) )
    minAngleLastMAF = minval ( GeodAngle2d(:,noMAFs) )
    nullify ( mif1GeodAngle, allgeodangle )
    call Allocate_test ( mif1GeodAngle, noMAFs, 'mif1Geodangle', ModuleName )
    call Allocate_test ( allGeodAngle, noMIFs, noMAFs, 'allGeodangle', ModuleName )
    mif1GeodAngle = GeodAngle2d(1,1:noMAFs)
    allGeodAngle = GeodAngle2d(1:noMIFs,1:noMAFs)
    if ( .not. isMonotonic(mif1GeodAngle) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'mif1GeodAngle is not monotonic--will try anyway' )
      call dump( mif1GeodAngle, 'mif1GeodAngle (before)' )
      call monotonize(mif1GeodAngle)
      call dump( mif1GeodAngle, 'mif1GeodAngle (after)' )
    end if

    ! Get or guess the start of the next chunk.
    i = noMAFs - chunk%noMAFsUpperOverlap + 1
    if ( deebughere ) then
      call outputNamedValue( 'i', i )
      call outputNamedValue( 'NoMAFs', size(FullArray2d, 2) )
      call outputNamedValue( 'maxAngle + spacing', maxAngle + spacing )
    endif
    if ( i < 1 ) then
      if ( deebughere ) then
        call output ( 'While constructing regular hGrid ', advance='yes' )
        call output ( 'minAngle: ' )
        call output ( minAngle, format='(F7.2)' )
        call output ( ' maxAngle: ' )
        call output ( maxAngle, format='(F7.2)' )
        call output ( 'i: ' )
        call output ( i, advance='yes' )
        call output ( 'noMAFs: ' )
        call output ( noMAFs, advance='yes' )
        call output ( 'chunk%noMAFsUpperOverlap: ' )
        call output ( chunk%noMAFsUpperOverlap, advance='yes' )
        call dump(chunk)
      endif
      if ( .not. associated(L1BFile) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Obviously impossible guess where to start next chunk' )
      
    else if ( i < size(GeodAngle2d, 2) + 1 ) then ! noMAFsInFile + 1 ) then
      nextAngle = GeodAngle2d(1, i)
    else
      nextAngle = maxAngle + spacing
    end if

    ! Now choose the geodetic angles for the hGrid
    if ( .not. single ) then 
      ! First identify the first point - 
      ! the one closest to the start of the first MAF
      first = origin + spacing * int ( (minAngle-origin)/spacing )
      delta = first - minAngle            ! So +ve means first could be smaller
      if ( delta > spacing/2 ) then
        first = first - spacing
      else if ( delta < -spacing/2 ) then
        first = first + spacing
      end if
    
      ! Now work out the last point in a similar manner
      if ( extendible ) then
        last = origin + spacing * (1 + int ( (maxAngle-origin)/spacing ) )
        delta = last - maxAngle            ! So +ve means last could be smaller
        if ( delta > 3*spacing/2 ) then
          last = last - spacing
        else if ( delta < -spacing/2 ) then
          last = last + spacing
        end if
      else
        last = origin + spacing * int ( (maxAngle-origin)/spacing )
        delta = last - maxAngle            ! So +ve means last could be smaller
        if ( delta > spacing/2 ) then
          last = last - spacing
        else if ( delta < -spacing/2 ) then
          last = last + spacing
        end if
      end if
    else
      ! The 'single' option is typically used for running single profile retrievals
      ! in a debug mode.  In order to ensure we choose the same profile for each MAF
      ! that a phiWindow=0 forward model would select, we have to do some extra work.

      ! First check this is a sane request
      if ( noMAFs /= 1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Single hGrid option set but more than one MAF in chunk' )
      ! Now construct a temporary HGrid that spans an orbit either side of this MAF
      ! Overkill I know but not harmfull.
      first = origin + spacing * int ( (minAngle-origin)/spacing ) - 360.0_rk
      last = first + 720.0_rk
      n = ( last - first ) / spacing
      nullify ( tmpAngle )
      call Allocate_test ( tmpAngle, n, 'tmpAngle', ModuleName )
      do i = 1, n
        tmpAngle ( i ) = first + (i-1) * spacing
      end do
      i = FindClosestMatch ( tmpAngle, GeodAngle2d, 1 )
      first = tmpAngle ( i )
      last = first
      call Deallocate_test ( tmpAngle, 'tmpAngle', ModuleName )
    end if

    ! Done with the L1B data
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )
    call Deallocate_test ( fullArray2d, 'fullArray2d', ModuleName )

    ! Now in the case where we have overlaps, let's try and have the
    ! first and last profile inside the MAF range
    if ( .not. single ) then
      if ( chunk%noMAFsLowerOverlap > 0 ) then
        setFirstLoop: do
          if ( first >= last ) exit setFirstLoop
          if ( first > maxAngleFirstMAF ) exit setFirstLoop
          first = first + spacing
        end do setFirstLoop
      else 
        if ( .not. insetOverlaps .and. first > minAngle ) first = first - spacing
      end if
      if ( chunk%noMAFsUpperOverlap > 0 ) then
        setLastLoop: do
          if ( last <= first ) exit setLastLoop
          if ( extendible ) then
            if ( last < max( nextAngle, minAngleLastMAF ) ) exit setLastLoop
          else
            if ( last < minAngleLastMAF ) exit setLastLoop
          end if
          last = last - spacing
        end do setLastLoop
      else
        if ( .not. insetOverlaps .and. last < maxAngle ) last = last + spacing
      end if
    end if

    ! Now work out how many profiles that is and lay them down
    hGrid%noProfs = nint( (last-first) / spacing ) + 1
    call CreateEmptyHGrid ( hGrid )
    do i = 1, hGrid%noProfs
      hGrid%phi(1,i) = first + (i-1)*spacing
    end do

    if ( verbose ) then
      call output ( 'Constructing regular hGrid ' )
      if ( hgrid%name /= 0 ) call display_string ( hgrid%name )
      call output ( '', advance='yes' )
      call outputNamedValue ( 'Number of profiles', hGrid%noProfs )
      call output ( 'minAngle: ' )
      call output ( minAngle, format='(F7.2)' )
      call output ( ' maxAngle: ' )
      call output ( maxAngle, format='(F7.2)' )
      call output ( ' nextAngle: ' )
      call output ( nextAngle, format='(F7.2)', advance='yes' )
      call output ( 'Num profiles: ' )
      call output ( hGrid%noProfs )
      call output ( ' Spacing: ' )
      call output ( spacing )
      call output ( ' first: ' )
      call output ( first, format='(F7.2)' )
      call output ( ' last: ' )
      call output ( last, format='(F7.2)', advance='yes' )
      call output ( ' firstMAFIndex: ' )
      call output ( chunk%firstMAFIndex )
      call output ( ' lastMAFIndex: ' )
      call output ( chunk%lastMAFIndex )
      call output ( ' extendible: ' )
      call output ( extendible, advance='yes' )
      call output ( ' forbidoverspill: ' )
      call output ( hGrid%forbidoverspill, advance='yes' )
      call output ( ' allowPriorOverlaps: ' )
      call output ( ChunkDivideConfig%allowPriorOverlaps, advance='yes' )
      call output ( ' allowPostOverlaps: ' )
      call output ( ChunkDivideConfig%allowPostOverlaps, advance='yes' )
    end if

    ! Now fill the other geolocation information, first latitude
    ! Get orbital inclination
    ! l1bItemName = AssembleL1BQtyName ( "scOrbIncl", hdfVersion, .false. )
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file not in database' )
    call L1BGeoLocation ( filedatabase, "scOrbIncl", &
      & instrumentModuleName, fullArray )
    call L1BSubsample ( chunk, FullArray, values=scOrbIncl )
    if ( deebughere ) then
      call dump(scOrbIncl, "scOrbIncl")
    end if
    ! Use the average of all the first MIFs to get inclination for chunk
    incline = sum ( scOrbIncl ) / size(scOrbIncl) ! noMAFs
    if ( noMAFs /= size(scOrbIncl) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'num of MAFs /= size(orb inclination qty)' )
    if ( deebughere ) then
      call dump( (/incline/), 'Average inclination')
    end if
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )
    call Deallocate_test ( scOrbIncl, 'scOrbIncl', ModuleName )
    hGrid%geodLat = phi_to_lat_deg ( hGrid%phi, incline )

    ! Now longitude
    call EmpiricalLongitude ( hGrid%phi(1,:), hGrid%lon(1,:) )

    call L1BGeoLocation ( filedatabase, "MAFStartTimeTAI", &
      & instrumentModuleName, fullArray )
    call L1BSubsample ( chunk, FullArray, values=TAI )
    if ( deebughere ) then
      call dump( TAI, "MAFStartTimeTAI" // ' (before interpolating)')
      call selfdiff( TAI, "MAFStartTimeTAI"  )
    end if
    if ( chunk%firstMAFIndex /= chunk%lastMAFIndex ) then
      if ( deebughere ) call outputNamedValue( 'shape(hGrid%phi)', size(hGrid%phi) )
      if ( deebughere ) call outputNamedValue( 'noProfs', hGrid%noProfs )
      call InterpolateValues ( mif1GeodAngle, TAI, &
        & hGrid%phi(1,:), hGrid%time(1,:), &
        & method='Linear', extrapolate='Allow' )
    else
      ! Case where only single MAF per chunk, treat it specially
      hGrid%time = TAI(1) + &
        & ( OrbitalPeriod/360.0 ) * ( hGrid%phi - mif1GeodAngle(1) )
    end if
    if ( deebughere ) then
      call dump( hGrid%time(1,:), "MAFStartTimeTAI" // ' (after interpolating)' )
      call output( 'geod angle, before and after interpolating', &
        & advance='yes' )
      call dump( mif1GeodAngle, 'before' )
      call dump( HGrid%phi(1,:), 'after' )
    end if
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )
    call Deallocate_test ( TAI, 'TAI', ModuleName )
      
    ! Solar time
    ! First get fractional day, note this neglects leap seconds.
    ! Perhaps fix this later !???????? NJL. We do have access to the
    ! UTC ascii time field, perhaps we could use that?
    hGrid%solarTime = modulo ( hGrid%time, secondsInDay ) / secondsInDay
    ! Now correct for longitude and convert to hours
    hGrid%solarTime = 24.0 * ( hGrid%solarTime + hGrid%lon/360.0 )
    hGrid%solarTime = modulo ( hGrid%solarTime, 24.0_rk )
    if ( deebughere ) call output ( 'Content with l1boa solartime', advance='yes' )

    ! Solar zenith
    call L1BGeoLocation ( filedatabase, "tpSolarZenith", &
      & instrumentModuleName,  fullArray, fullArray2d )
    call L1BSubsample ( chunk, FullArray, fullArray2d, SolarZenith, SolarZenith2d )
    if ( deebughere ) then
      call dump( SolarZenith, "SolarZenith" // &
        & ' (before interpolating)')
    end if
    ! This we'll have to do with straight interpolation
    call InterpolateValues ( mif1GeodAngle, SolarZenith, &
      & hGrid%phi(1,:), hGrid%solarZenith(1,:), &
      & method='Linear', extrapolate='Allow' )
    if ( deebughere ) then
      call dump( hGrid%solarZenith(1,:), "SolarZenith" // &
        & ' (after interpolating)' )
    end if
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )
    call Deallocate_test ( SolarZenith, 'SolarZenith', ModuleName )
    call Deallocate_test ( fullArray2d, 'fullArray2d', ModuleName )
    call Deallocate_test ( SolarZenith2d, 'SolarZenith2d', ModuleName )

    ! Line of sight angle
    call L1BGeoLocation ( filedatabase, "tpLosAngle", &
      & instrumentModuleName, fullArray )
    call L1BSubsample ( chunk, FullArray, values=LosAngle )
    if ( deebughere ) then
      call dump( LosAngle, 'LosAngle' // ' (before interpolating)')
    end if
    call InterpolateValues ( mif1GeodAngle, LosAngle, &
      & hGrid%phi(1,:), hGrid%losAngle(1,:), &
      & method='Linear', extrapolate='Allow' )
    if ( deebughere ) then
      call dump( hGrid%losAngle(1,:), 'LosAngle' // ' (after interpolating)')
    end if
    call Deallocate_test ( fullArray, 'fullArray', ModuleName )
    call Deallocate_test ( LosAngle, 'LosAngle', ModuleName )
    hGrid%losAngle = modulo ( hGrid%losAngle, 360.0_rk )

    ! Now work out how much of this HGrid is overlap
    ! The deal will be the first legitimate profile is the first one who's phi
    ! is above the first non overlapped MAF.
    call Hunt ( hGrid%phi(1,:), mif1GeodAngle(chunk%noMAFsLowerOverlap+1), &
       & hGrid%noProfsLowerOverlap, allowTopValue=.true., allowBelowValue=.true. )
    if ( verbose ) &
      & call outputNamedValue ( 'hGrid%noProfsLowerOverlap (after Hunt)', &
        & hGrid%noProfsLowerOverlap )
    ! So the hunt returns the index of the last overlapped, which is
    ! the number we want to be in the overlap.
    call Hunt ( hGrid%phi(1,:), nextAngle, &
       & hGrid%noProfsUpperOverlap, allowTopValue=.true., allowBelowValue=.true. )
    ! Here the hunt returns the index of the last non overlapped profile
    ! So we do a subtraction to get the number in the overlap.
    hGrid%noProfsUpperOverlap = hGrid%noProfs - hGrid%noProfsUpperOverlap
    ! call outputNamedValue ( 'hGrid%noProfsUpperOverlap (after Hunt)', &
    !  & hGrid%noProfsUpperOverlap )

    if ( verbose .or. hGrid%noProfs == 0 ) then
      call output ( 'Initial Hgrid size: ' )
      call output ( hGrid%noProfs ) 
      call output ( ', overlaps: ' )
      call output ( hGrid%noProfsLowerOverlap )
      call output ( ', ' )
      call output ( hGrid%noProfsUpperOverlap, advance='yes' )
    end if

    ! Now, we want to ensure we don't spill beyond the processing time range.
    ! First an optional 'hard' limit.  This guarantees that no part of the
    ! chunk (even an overlap) spills beyond the processing range.  This is
    ! important for runs using sids L2GP data which has hard limits (would
    ! crash Fill otherwise).  However, it's not important for other runs, so is
    ! optional.
    if ( hGrid%forbidOverspill ) then
      call Hunt ( hGrid%time(1,:), processingRange%startTime, &
        & firstProfInRun, allowTopValue=.true., allowBelowValue=.true. )
      if ( deebughere ) then
        call dump( hGrid%time(1,:), 'Hgrid times')
        call dump( hGrid%time(1,:)-hGrid%time(1,1), 'Hgrid delta times')
        call output ( 'First profile in Run: ' )
        call output ( firstProfInRun, advance='no' )
        call output ( '    processingRange%startTime: ' )
        call output ( processingRange%startTime, advance='yes' )
      end if
      ! if ( firstProfInRun > 0 .and. .not. ChunkDivideConfig%allowPriorOverlaps ) then
      if ( firstProfInRun > 0 ) then
        if ( verbose ) &
          & call output ( 'hGrid starts before run trimming start.', &
          & advance='yes' )
        call TrimHGrid ( hGrid, -1, firstProfInRun )
      end if

      call Hunt ( hGrid%time(1,:), processingRange%endTime, &
        & lastProfInRun, allowTopValue=.true., allowBelowValue=.true. )
      if ( deebughere ) then
        call output ( 'Last profile in Run: ' )
        call output ( lastProfInRun, advance='no' )
        call output ( '    processingRange%endTime: ' )
        call output ( processingRange%endTime, advance='yes' )
      end if
      ! if ( lastProfInRun < hGrid%noProfs .and. .not. ChunkDivideConfig%allowPostOverlaps ) then
      if ( lastProfInRun < hGrid%noProfs ) then
        if ( verbose ) &
          & call output ( 'hGrid ends after run trimming end.',&
          & advance='yes' )
        call TrimHGrid ( hGrid, 1, hGrid%noProfs-lastProfInRun )
      end if
    end if
    ! Now a 'softer' limit that applies to all cases, this just moves the
    ! overlap regions around if necessary to deal with overspill.
    if ( hGrid%noProfsLowerOverlap+1 <= hGrid%noProfs ) then
      if ( hGrid%time(1,hGrid%noProfsLowerOverlap+1) < processingRange%startTime ) then
        if ( verbose ) &
          & call output ( &
          & 'Non overlapped part of hGrid starts before run, extending overlap.', &
          & advance='yes' )
        call Hunt ( hGrid%time(1,:), processingRange%startTime, &
          &   hGrid%noProfsLowerOverlap, allowTopValue=.true., allowBelowValue=.true. )
      end if
    end if

    if ( hGrid%noProfs-hGrid%noProfsUpperOverlap >= 1 ) then
      if ( hGrid%time(1,hGrid%noProfs-hGrid%noProfsUpperOverlap) > &
        & processingRange%endTime ) then
        if ( verbose ) &
          & call output ( &
          & 'Non overlapped part of hGrid end after run, extending overlap.', &
          & advance='yes' )
        call Hunt ( hGrid%time(1,:), processingRange%endTime, &
          & hGrid%noProfsUpperOverlap, allowTopValue=.true., allowBelowValue=.true. )
        hGrid%noProfsUpperOverlap = hGrid%noProfs - hGrid%noProfsUpperOverlap
      end if
    end if

    ! Now deal with the user requests
    if ( single ) then
       ! Set up for no overlaps
       if ( verbose ) call outputNamedValue( 'Before TrimHGrid', hGrid%noProfs )
       hGrid%noProfsLowerOverlap = 0
       hGrid%noProfsUpperOverlap = 0
       ! Delete all but the 'center' profile
       extra = hGrid%noProfs - 1
       left = extra / 2
       right = extra - left
       if ( verbose ) call outputNamedValue( '/single: left, right', (/ left, right /) )
       if ( left > 0 ) call TrimHGrid ( hGrid, -1, left )
       if ( right > 0 ) call TrimHGrid ( hGrid, 1, right )
    else if (hgrid%noProfs > 1) then
       if ( verbose ) call outputNamedValue( 'Before TrimHGrid', hGrid%noProfs )
       if ( maxLowerOverlap >= 0 .and. &
          & ( hGrid%noProfsLowerOverlap > maxLowerOverlap ) ) then
          if ( verbose ) call outputNamedValue( 'Lower overlap too big', hGrid%noProfs )
          call TrimHGrid ( hGrid, -1, hGrid%noProfsLowerOverlap - maxLowerOverlap )
       end if
       if ( maxUpperOverlap >= 0 .and. &
          & ( hGrid%noProfsUpperOverlap > maxUpperOverlap ) ) then
          if ( verbose ) call outputNamedValue( 'Upper overlap too big', hGrid%noProfs )
          call TrimHGrid ( hGrid, 1, hGrid%noProfsUpperOverlap - maxUpperOverlap )
       end if
       ! Setting maxLowerOverlap to -m means we always trim away 
       ! the first m profiles
       if ( maxLowerOverlap < 0 .and. &
         & (hgrid%noProfs > 2) ) then
         if ( verbose ) call outputNamedValue( 'Trimmed away first m profile; m', -maxLowerOverlap )
         call TrimHGrid ( hGrid, -1, -maxLowerOverlap )
       endif
       ! Setting maxUpperOverlap to -m means we always trim away 
       ! the last m profiles
       if ( maxUpperOverlap < 0 .and. &
         & (hgrid%noProfs > 2) ) then
         if ( verbose ) call outputNamedValue( 'Trimmed away last m profile; m', -maxUpperOverlap )
         call TrimHGrid ( hGrid, 1, -maxUpperOverlap )
       endif
    else ! if there is only one profile, then don't care about overlap
       hGrid%noProfsLowerOverlap = 0
       hGrid%noProfsUpperOverlap = 0
    end if

    if ( hGrid%noProfs == 0 ) then
      call dump(chunk)
      if ( warnIfNoProfs ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, 'No profiles in hGrid' )
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, 'No profiles in hGrid' )
      end if
    end if

    ! Finally we're done.
    if ( verbose .or. hGrid%noProfs == 0 ) then
      call output ( 'Final Hgrid size: ' )
      call output ( hGrid%noProfs )
      call output ( ', overlaps: ' )
      call output ( hGrid%noProfsLowerOverlap )
      call output ( ', ' )
      call output ( hGrid%noProfsUpperOverlap, advance='yes' )
    end if

    ! That's it
    call Deallocate_test ( mif1GeodAngle, 'mif1GeodAngle', ModuleName )
    call Deallocate_test ( allGeodAngle, 'allGeodAngle', ModuleName )

    call trace_end ( "CreateRegularHGrid", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

  end subroutine CreateRegularHGrid

  ! --------------------------------------- DealWithObstructions -----
  subroutine DealWithObstructions ( HGrid, obstructions, DestroyOld )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Chunkdivide_M, only: Obstruction_T
    use HgridsDatabase, only: Hgrid_T, Createemptyhgrid, Destroyhgridcontents
    use HighOutput, only: LetsDebug
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    ! Args
    type (HGRID_T), pointer :: HGRID
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    logical, optional, intent(in) :: DestroyOld
    ! This routine modifies the chunks according to the information
    ! given in the obstructions.

    ! Local variables
    logical :: deeBug
    type (HGRID_T), target :: NEWHGrid ! The one we'll create create
    integer :: FIRSTMAF                 ! Index of first MAF in range
    integer :: LASTMAF                  ! Index of last MAF in range
    integer :: Me = -1                  ! String index for tracing
    integer :: MAF                      ! Index of MAF for wall
    logical :: mayDestroyOld
    integer :: newProfile               ! Counter in newHGrid
    integer :: OBSTRUCTION              ! Loop counter
    integer :: PROFILE                  ! Loop counter
    logical, dimension(:), pointer :: obstructed => null()

    ! Executable code
    deebug = LetsDebug ( 'hgrid', 1 )
    ! 1st--some short-circuits
    if ( hGrid%noProfs < 1 ) return
    if ( .not. associated(obstructions) ) return
    if ( .not. associated(hGrid%phi) ) return
    call trace_begin ( me, "DealWithObstructions", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )
    mayDestroyOld = .false.
    if ( present(DestroyOld) ) mayDestroyOld = DestroyOld
    
    ! Next--fill obstructed array
    nullify (obstructed)
    call Allocate_Test ( obstructed, hGrid%noProfs, &
      & 'obstructed', ModuleName )
    obstructed = .false.
    do obstruction=1, size(obstructions)
      if ( obstructions(obstruction)%range ) then
        ! A range obstruction
        firstMAF = obstructions(obstruction)%mafs(1)
        lastMAF = obstructions(obstruction)%mafs(2)
        if ( firstMAF < 0 .or. firstMAF+1 > hGrid%noProfs ) cycle
        if ( lastMAF < firstMAF .or. lastMAF+1 > hGrid%noProfs ) cycle
        obstructed(firstMAF+1:lastMAF+1) = .true. ! MAF indexes start at 0
      else
        ! A wall obstruction
        maf = obstructions(obstruction)%mafs(1)
        if ( maf < 0 .or. maf+1 > hGrid%noProfs ) cycle
        obstructed(maf+1:maf+1) = .true. ! MAF indexes start at 0
      end if
    end do
    if ( any(obstructed) ) then
      newHGrid%noProfs = count( .not. obstructed )
      call CreateEmptyHGrid(newHGrid)
      newHGrid%name = hgrid%name
      newProfile = 0
      do profile=1, hGrid%noProfs
        if ( obstructed(profile) ) cycle
        newProfile = newProfile + 1
        newHGrid%phi        (1, newProfile) = hgrid%phi        (1, profile)
        newHGrid%geodLat    (1, newProfile) = hgrid%geodLat    (1, profile)
        newHGrid%lon        (1, newProfile) = hgrid%lon        (1, profile)
        newHGrid%time       (1, newProfile) = hgrid%time       (1, profile)
        newHGrid%solarTime  (1, newProfile) = hgrid%solarTime  (1, profile)
        newHGrid%solarZenith(1, newProfile) = hgrid%solarZenith(1, profile)
        newHGrid%losAngle   (1, newProfile) = hgrid%losAngle   (1, profile)
      end do
      if ( hgrid%noProfsLowerOverlap > 0 ) then
        newHGrid%noProfsLowerOverlap = hgrid%noProfsLowerOverlap &
          & - count( obstructed(1:hgrid%noProfsLowerOverlap) )
      end if
      if ( hgrid%noProfsUpperOverlap > 0 ) then
        newHGrid%noProfsUpperOverlap = hgrid%noProfsUpperOverlap &
          & - count( obstructed(hgrid%noProfs-hgrid%noProfsUpperOverlap+1:hgrid%noProfs) )
      end if
      if ( mayDestroyOld ) then
        call DestroyHGridContents(hGrid)
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Destroying old HGrid known to corrupt other HGrid' )
      end if
      hGrid => newHGrid
    end if
    call DeAllocate_Test( obstructed, 'obstructed', ModuleName )

    call trace_end ( "DealWithObstructions", &
      & cond=toggle(gen) .and. levels(gen) > 1 .and. .not. computingOffsets )

  end subroutine DealWithObstructions

  ! ---------------------------------------------  DestroyHGridGeoLocations  -----
  ! Deallocate all components of HGridGeolocations.
  subroutine DestroyHGridGeoLocations
    use Allocate_Deallocate, only: Deallocate_Test
    call Deallocate_test ( HGridGeolocations%MAFStartTimeTAI, 'MAFStartTimeTAI', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzGeodAngle   , 'GHzGeodAngle   ', ModuleName )
    call Deallocate_test ( HGridGeolocations%Orbincl        , 'Orbincl        ', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzGeodAlt     , 'GHzGeodAlt     ', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzGeodLat     , 'GHzGeodLat     ', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzLon         , 'GHzLon         ', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzSolarTime   , 'GHzSolarTime   ', ModuleName )
    call Deallocate_test ( HGridGeolocations%GHzSolarZenith , 'GHzSolarZenith ', ModuleName )
  end subroutine DestroyHGridGeoLocations

  ! -------------------------------------  DumpChunkHGridGeometry  -----
  subroutine DumpChunkHGridGeometry ( hGrid, chunk, &
    & instrumentModuleName, filedatabase )

    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Chunks_M, only: MLSChunk_T
    use HgridsDatabase, only: Hgrid_T
    use L1BData, only: DeallocateL1BData, L1BData_T, ReadL1BData, &
      & AssembleL1Bqtyname
    use MLSKinds, only: R8
    use Output_M, only: Output
    use String_Table, only: Display_String

    type (HGrid_T), intent(in) :: HGRID
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    character (len=*) :: INSTRUMENTMODULENAME
    type (MLSChunk_T), intent(in) :: CHUNK

    ! Local parameters
    real (r8), parameter :: BINSIZE=0.05 ! Width of one character in degrees
    integer, parameter :: NOLINES=8     ! Number of lines that make up a scan
    integer, parameter :: WIDTH=75      ! Number of columns to print at a time

    ! Local variables
    character (len=1), dimension(:,:), pointer :: TEXT
    integer :: CHAR                     ! Loop counter
    integer :: CHARMAX                  ! Last character to fill
    integer :: CHARMIN                  ! First character to fill
    integer :: FIRSTMIF                 ! First mif to consider
    integer :: FLAG                     ! From L1B read
    integer :: LASTMIF                  ! Last mif to consider
    integer :: LINE                     ! Loop counter
    integer :: MAF                      ! Loop counter
    integer :: MIFSPERLINE              ! Number of mifs covered by a line
    integer :: NOEMPTYWINDOWS           ! Number of empty windows
    integer :: NOBINS                   ! Number of characters to print
    integer :: NOMAFS                   ! From L1B read
    integer :: NOMIFS                   ! Deduced
    integer :: NOWINDOWS                ! Number of blocks of text
    integer :: PROF                     ! Loop counter
    integer :: START                    ! Index
    integer :: WINDOW                   ! Loop counter
    integer :: WINDOWSIZE               ! Number of characters in this window
    real (r8) :: PHIMAX                 ! Maximum value to consider
    real (r8) :: PHIMIN                 ! Minimum value to consider
    real (r8) :: THISPHIMAX             ! Max phi in a line for a maf
    real (r8) :: THISPHIMIN             ! Min phi in a line for a maf
    real (r8), dimension(:,:), pointer :: MIFPHI ! Tangent phis
    type (L1BData_T) :: L1BFIELD

    integer ::  hdfVersion
    character(len=NameLen) :: l1bItemName
    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code
    call output ( "Dumping geometry for HGrid: " )
    if ( hGrid%name /= 0 ) then
      call display_string ( hGrid%name, strip=.true., advance='yes' )
    else
      call output ( "<no-name>", advance='yes' )
    end if

      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      hdfversion = L1BFile%HDFVersion

    ! Read the geodetic angle from the L1Bfile
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpGeodAngle", hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex, &
      & dontPad=.true. )
    mifPhi => l1bField%dpField(1,:,:)

    phiMin = minval ( mifPhi )
    phiMax = maxval ( mifPhi )
    if ( hGrid%noProfs > 0 ) then
      phiMin = min ( phiMin, hGrid%phi(1,1) )
      phiMax = max ( phiMax, hGrid%phi(1,hGrid%noProfs) )
    end if

    ! Now setup the text to print
    noBins = ( phiMax - phiMin ) / binSize
    call allocate_test ( text, noBins, noLines, 'text', moduleName, fill='' )

    ! Trivial stuff to setup
    noMIFs = size ( mifPhi, 1 )
    mifsPerLine = noMIFs / noLines
    
    ! Now loop over the mafs and fill them up
    do maf = 1, noMAFs
      do line = 1, noLines
        firstMIF = (line-1) * mifsPerLine + 1
        lastMIF = firstMIF + mifsPerLine - 1
        if ( line == noLines ) lastMIF = noMIFs
        thisPhiMin = minval ( mifPhi ( firstMIF:lastMIF, maf ) )
        thisPhiMax = maxval ( mifPhi ( firstMIF:lastMIF, maf ) )
        charMin = ( thisPhiMin - phiMin ) / binSize + 1
        charMax = ( thisPhiMax - phiMin ) / binSize + 1
        charMin = min ( max ( charMin, 1 ), noBins )
        charMax = min ( max ( charMax, 1 ), noBins )
        if ( ( maf < chunk%noMAFsLowerOverlap+1 ) .or. &
          &  ( maf > chunk%lastMAFIndex - chunk%firstMAFIndex + 1 - &
          & chunk%noMAFsUpperOverlap) ) then
          text ( charMin:charMax , line ) = '*'
        else
          text ( charMin:charMax , line ) = '#'
        end if
      end do
    end do

    ! Now loop over the profiles and place them
    do prof = 1, hGrid%noProfs
      charMin = ( hGrid%phi(1,prof) - phiMin ) / binSize + 1
      charMin = min ( max ( charMin, 1 ), noBins )
      if ( ( prof < hGrid%noProfsLowerOverlap+1 ) .or. &
        &  ( prof > hGrid%noProfs-hGrid%noProfsUpperOverlap) ) then
        text ( charMin, : ) = ':'
      else
        text ( charMin, : ) = '|'
      end if
    end do

    ! Now print them out
    noWindows = noBins / width
    noEmptyWindows = 0
    if ( mod ( noBins, width ) /= 0 ) noWindows = noWindows + 1
    do window = 1, noWindows
      windowSize = width
      start = ( window-1 ) * width + 1
      if ( window == noWindows ) windowSize = mod ( noBins, width )
      if ( all ( text ( start : start+windowSize-1, : ) == ' ' ) ) then
        noEmptyWindows = noEmptyWindows + 1
      else
        if ( noEmptyWindows > 0 ) then
          call output ( noEmptyWindows )
          call output ( ' empty windows suppressed', advance='yes' )
          noEmptyWindows = 0
        end if
        call output ( ' ', advance='yes' )
        do line = noLines, 1, -1
          do char = start, start+windowSize-2
            call output ( text ( char, line ) )
          end do
          call output  ( text ( start+windowSize-1, line ), advance='yes' ) 
        end do
        do char = start, start+windowSize-2
          call output ( '-' )
        end do
        call output ( '-', advance='yes' )
      end if
    end do

    ! Finish off
    call deallocate_test ( text, 'text', moduleName )

    call DeallocateL1BData ( l1bField )

  end subroutine DumpChunkHGridGeometry

  ! -------------------------------------  ComputeAllHGridOffsets  -----
  subroutine ComputeAllHGridOffsets ( root, first_section, chunks, filedatabase, &
    & l2gpDatabase, processingRange )
    ! This routine goes through the L2CF up to an Output section to accumulate
    ! HGrid sizes from the Construct sections, and through the L1 file to work
    ! out how big each HGrid is going to be
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
    use Chunks_M, only: MLSChunk_T
    use ChunkDivide_M, only: ChunkDivideConfig
    use Dump_0, only: Dump
    use HGridsDatabase, only: HGrid_T, CopyHGrid, DestroyHGridContents, Dump
    use HighOutput, only: BeVerbose, LetsDebug, OutputNamedValue
    use Init_Tables_Module, only: Z_Construct, S_Hgrid, Z_Output
    use L1BData, only: CheckForCorruptFileDatabase
    use L2GPData, only: L2GPData_T
    use MLSKinds, only: Rk => R8
    use MLSL2Options, only: SpecialDumpFile
    use Moretree, only: Get_Spec_Id
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Next_Tree_Node_M, only: Init_Next_Tree_Node, Next_Tree_Node, &
      & Next_Tree_Node_State
    use Output_M, only: Output, RevertOutput, SwitchOutput
    use Time_M, only: SayTime, Time_Now
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use Tree, only: Subtree, Node_Id, Decoration
    use Tree_Types, only: N_Named
    ! Dummy arguments
    integer, intent(in) :: Root         ! of the entire tree
    integer, intent(in) :: First_Section ! First subtree after definitions
    type(MLSChunk_T), dimension(:), intent(inout) :: CHUNKS
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type(L2GPData_T), dimension(:), pointer :: L2GPDatabase
    type(TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    ! Local variables
    integer :: CHUNK                    ! Loop counter
    integer :: C                        ! Inner loop counter
    logical :: deeBug
    integer :: HGRID                    ! Loop counter
    integer :: How_many
    integer :: Me = -1                  ! String index for trace
    integer :: NOHGRIDS                 ! Number of hGrids
    integer :: SON                      ! Tree node
    integer :: sum1, sum2, sum3
    integer :: GSON                     ! son of son
    integer :: KEY                      ! Tree node
    integer :: KEYINDEX
    type(HGrid_T) :: dummyHGrid         ! A temporary hGrid
    type(HGrid_T), dimension(size(chunks)) :: FirstHGrid     ! The "std" HGrid
    type(next_tree_node_state) :: State1 ! while hunting for Construct sections
    type(next_tree_node_state) :: State2 ! within Construct sections
    integer, dimension(100) :: keyArray
    integer, dimension(:), pointer :: LowerOverlaps => null()
    integer :: nkeys
    integer :: nsections
    logical :: verbose
    logical :: verboser
    real(rk), dimension(:), pointer :: MAFStartTimeTAI, OrbIncl, &
      & GeodAngle, GeodAlt, GeodLat, Lon, LosAngle, SolarTime, SolarZenith
    integer :: InstrumentModule
    character (len=NameLen) :: InstrumentModuleName
    integer, dimension(size(chunks)) :: ChunksWithoutProfiles
    ! Executable code
    computingOffsets = .true.
    deebug = LetsDebug ( 'hgrid', 1 )
    verbose = beVerbose ( 'hgrid', 1 )
    verboser = beVerbose ( 'hgrid', 2 )
    call OutputnamedValue ( 'deebug', deebug )
    call OutputnamedValue ( 'verbose', verbose )
    call OutputnamedValue ( 'verboser', verboser )
    call trace_begin ( me, "ComputeAllHGridOffsets", root, &
      & cond=toggle(gen) )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ! Slightly clever loop here, the zero run doesn't do anything except
    ! count the number of hGrids.  The other run works out how big each HGrid
    ! is going to be for each chunk.  This is stored in chunks%hGridOffsets.
    ! Finally we accumulate these to get offsets.
    noHGrids = 0
    nkeys = 0
    if ( deebug .and. .false. ) then
      call outputNamedValue ( 'size(chunks)', size(chunks) )
      call output( 'Before loop of chunks', advance='yes' )
      call CheckForCorruptFileDatabase( filedatabase )
    endif
    do chunk = 0, size(chunks)
      if ( chunk > 0 ) chunks(chunk)%hGridOffsets = 0 ! Initialize
      call time_now ( t0 )
      call time_now ( t1 )
      hGrid = 1
      ! if ( chunk > 20) stop
      ! Loop over all the setions in the l2cf, look for construct sections
      if ( verbose ) &
        & call outputNamedValue ( 'Computing offsets for chunk number', chunk, &
        & Before='*   ', After = '  *' )
      call init_next_tree_node ( state1 )
      nsections = 0
      keyindex = 0
      sectionLoop: do
        nsections = nsections + 1
        son = next_tree_node ( root, state1, start=first_section, traceLevel=4 )
        if ( son == 0 ) exit
        call trace_begin ( me, "ComputeAllHGridOffsets.loop", son, &
          & cond=toggle(gen) .and. levels(gen) > 1 .and. verboser )
        select case ( decoration ( subtree ( 1, son ) ) )
        case ( z_construct )
          ! Now loop through the construct section and identify the hGrids
          call init_next_tree_node ( state2 )
          do ! if chunk == 0, loop over every statement; if chunk > 0, just the keys
            if ( chunk == 0 ) then
              ! For the 'zeroth' pass just count up the hgrids
              gson = next_tree_node ( son, state2, traceLevel=5 )
              if ( gson == 0 ) exit ! Past end-of-section?
              if ( node_id(gson) /= n_named ) cycle ! Is spec labeled?
              key = subtree(2,gson)
              if ( get_spec_id(key) /= s_hGrid ) cycle ! Is it an hgrid?
              nkeys = nkeys + 1
              keyArray(nkeys) = key
              noHGrids = noHGrids + 1
            else
              keyindex = keyindex + 1
              if ( keyindex > nkeys ) exit ! did all keys?
              key = keyArray(keyIndex)
              ! nullify ( dummyHGrid )
              if ( verbose ) then
                call output ( 'Creating dummyHGrid', advance='yes' )
                call sayTime( 'Getting set to call CreateHGrid' )
                call outputnamedValue ( 'This section, keyindex #', &
                  & (/ nsections, keyindex /) )
                call time_now( t1 )
              endif
              ! call Dump ( filedatabase, details=2 )
              if ( chunk == 1 .and. deebug .and. .false. ) then
                call outputNamedValue( 'Before making dummyHGrid', chunk )
                call CheckForCorruptFileDatabase( filedatabase )
              endif
              dummyHGrid = CreateHGridFromMLSCFInfo ( 0, key, filedatabase, l2gpDatabase, &
                & processingRange, chunks(chunk), onlyComputingOffsets=.true., &
                & check=.false. )
                ! & check=(chunk == 1) )
              if ( chunk == 1 .and. deebug.and. .false.  ) then
                call outputNamedValue( 'During chunk#', chunk )
                call CheckForCorruptFileDatabase( filedatabase )
              endif
              if ( keyindex == 1 ) then
                if ( verbose ) call output ( 'Creating firstHGrid', advance='yes' )
                  ! firstHGrid(chunk) = CreateHGridFromMLSCFInfo ( 0, key, filedatabase, &
                  !& l2gpDatabase, processingRange, chunks(chunk), &
                  !  & onlyComputingOffsets=.true. )
                call time_now ( t1 )
                call copyHGrid( dummyHGrid, firstHGrid(chunk) )
                if ( verbose ) then
                  call Dump ( firstHGrid(chunk) )
                  call outputNamedValue( 'num Profs copied', dummyHGrid%noProfs )
                  call sayTime( 'Copying HGrid' )
                endif
                InstrumentModule = dummyHGrid%module
                call time_now( t1 )
              endif
              if ( chunk == 1 ) then
                LowerOverlaps(hGrid) = dummyHGrid%noProfsLowerOverlap
                if ( ChunkDivideConfig%allowPriorOverlaps .and. .false. ) then
                  chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                  & dummyHGrid%noProfsUpperOverlap
                else
                  chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                  & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                end if
                if ( deebug .and. .false. ) then
                  call outputNamedValue( 'During chunk#', chunk )
                  call CheckForCorruptFileDatabase( filedatabase )
                endif
              else if ( chunk == size(chunks) ) then
                if ( ChunkDivideConfig%allowPostOverlaps ) then
                  chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                  & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                else
                  chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                  & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                end if
              else
                chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                if ( verbose .and. &
                  & any( dummyHGrid%phi(1,:) /= 0.) .and. key==9927 ) then
                  call outputNamedValue( 'num Profs', dummyHGrid%noProfs )
                  call outputNamedValue( 'key', key )
                  call outputNamedValue( 'phis', dummyHGrid%phi(1,:) )
                endif
              end if
              ! Sometimes we have been computing negative offsets which is .. offsetting
              if ( chunk > 0 ) &
                & chunks(chunk)%hGridOffsets(hGrid) = max( 0, chunks(chunk)%hGridOffsets(hGrid) )
              call time_now ( t1 )
              if ( DEEBUG ) &
                & call dump(dummyHGrid)
              call DestroyHGridContents ( dummyHGrid )
              if ( verbose ) call sayTime ( 'Dumping and destroying this HGrid' )
              call time_now( t1 )
              hGrid = hGrid + 1
            end if
          end do
        case ( z_output )
          ! exit sectionLoop ! No, don't exit the section loop w/o trace_end
        case default
        end select
        call trace_end ( "ComputeAllHGridOffsets.loop", &
          & cond=toggle(gen) .and. levels(gen) > 1 .and. verboser )
      end do sectionLoop
      if ( verbose ) then
        call outputNamedValue( 'number of keys', nKeys )
        call outputNamedValue( 'number of sections', nSections )
        call Dump ( keyArray(:nKeys), 'keys', width=5 )
      endif

      ! If this is the first time round, we now know how many hGrids there
      ! are, set up our arrays
      ! call outputNamedValue ( 'noHGrids', noHGrids )
      ! call outputNamedValue ( 'hGrid', hGrid )
      if ( chunk == 0 ) then
        do c = 1, size ( chunks )
          call Allocate_Test ( chunks(c)%hGridOffsets, noHGrids, &
            & 'chunks(?)%hGridOffsets', ModuleName )
          call Allocate_Test ( chunks(c)%hGridTotals, noHGrids, &
            & 'chunks(?)%hGridTotals', ModuleName )
          call Allocate_Test ( LowerOverlaps, noHGrids, &
            & 'LowerOverlaps', ModuleName )
        end do
        if ( deebug.and. .false.  ) then
          call output( 'After 0th chunk', advance='yes' )
          call CheckForCorruptFileDatabase( filedatabase )
        endif
      elseif ( nkeys /= noHGrids ) then
        call outputnamedValue ( 'nKeys, noHGrids', (/ nkeys, noHGrids /) )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Got a different number of hGrids than expected!' )
      elseif ( chunk < 10 .and. deebug.and. .false.  ) then
        call outputNamedValue( 'After chunk#', chunk )
        call CheckForCorruptFileDatabase( filedatabase )
      end if
    end do                              ! Chunk loop
    ! call output( 'After loop of chunks', advance='yes' )
    ! call CheckForCorruptFileDatabase( filedatabase )
    ! call output( 'Done with chunk Loop', advance='yes' )
    
    if ( DEEBUG ) then
      sum1 = 0
      sum2 = 0
      sum3 = 0
      call output ( "number of profiles in each chunk/grid: " , advance='yes')
      do chunk = 1, size ( chunks )
        call output ( chunk, before='Chunk ', after=': ' )
        call output ( chunks(chunk)%hGridOffsets, advance='yes' )
        sum1 = sum1 + chunks(chunk)%hGridOffsets(1)
        if ( size(chunks(chunk)%hGridOffsets) > 1 ) then
          sum2 = sum2 + chunks(chunk)%hGridOffsets(2)
          if ( size(chunks(chunk)%hGridOffsets) > 2 ) &
            & sum3 = sum3 + chunks(chunk)%hGridOffsets(3)
        end if
      end do
      call output ( "Total number of profiles" , advance='yes')
      call output ( (/ sum1, sum2, sum3 /), advance='yes' )
    end if
    ! Now accumulate hGridOffsets which currently contains the number
    ! of non-overlap profiles in each chunk/hGrid.  After this it will contain
    ! the accumulated number.  This is equivalent to storing the
    ! index of the last profile in each chunk.
    if ( verbose ) &
      & call output( chunks(1)%hGridOffsets(1), before='Chunk 1: ', advance='yes' )
    do chunk = 2, size ( chunks )
      chunks(chunk)%hGridOffsets = chunks(chunk)%hGridOffsets + &
        & chunks(chunk-1)%hGridOffsets
      if ( verbose ) then
        call output ( chunk, before='Chunk ', after=': ' )
        call output ( chunks(chunk)%hGridOffsets(1), advance='yes' )
      end if
    end do
    ! Now fill (i.e. spread) the hGridTotals array
    do chunk = 1, size ( chunks )
      chunks(chunk)%hGridTotals = chunks(size(chunks))%hGridOffsets
    end do
    ! Now move hGridOffsets back one to get the offsets we want.
    ! Each hGridOffsets will now contain the last profile index in the preceding
    ! chunk, which is the offset we desire.
    do chunk = size ( chunks ), 2, -1
      chunks(chunk)%hGridOffsets = chunks(chunk-1)%hGridOffsets
    end do
    chunks(1)%hGridOffsets = 0
    if ( DEEBUG .or. verbose ) then
      call output ( 'chunks(1)%hGridOffsets: ', advance='no' )
      call output ( chunks(1)%hGridOffsets, advance='yes' )
    end if
    if ( ChunkDivideConfig%allowPriorOverlaps .and. &
      & chunks(1)%noMAFsLowerOverlap > 0 .and. firstHGrid(1)%forbidOverspill ) then
      ! This dubious sequence disallows the very Prior overlaps you pretended to
      ! allow
      call dump( LowerOverlaps, 'LowerOverlaps' )
      chunks(1)%hGridOffsets = 0
      do chunk=1, size(chunks)
        chunks(chunk)%hGridOffsets = chunks(chunk)%hGridOffsets + LowerOverlaps
      end do
    else
      chunks(1)%hGridOffsets = 0
    end if
    
    if ( verbose ) then
      call output ( "Dumping offsets, hgridTotals for all chunks: " , &
        & advance='yes')
      do chunk = 1, size ( chunks )
        call output ( chunk, before='Chunk ', after=': ' )
        call output ( chunks(chunk)%hGridOffsets )
        call output ( ' totals ' )
        call output ( chunks(chunk)%hGridTotals, advance='yes' )
      end do
    end if

    if ( verbose ) then
      instrumentModuleName = 'none'
      if ( instrumentModule > 0 ) &
        & call GetModuleName ( instrumentModule, instrumentModuleName )
      ! Read the l1boa items we will need
      nullify( MAFStartTimeTAI, GeodAngle, GeodAlt, GeodLat, &
        & LosAngle, OrbIncl, SolarTime, SolarZenith )
      call L1BGeoLocation ( filedatabase, 'MAFStartTimeTAI  ', instrumentModuleName,  MAFStartTimeTAI, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/GeodAngle    ', instrumentModuleName,  GeodAngle, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/GeodAlt      ', instrumentModuleName,  GeodAlt, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/GeodLat      ', instrumentModuleName,  GeodLat, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/Lon          ', instrumentModuleName,  Lon, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/LosAngle     ', instrumentModuleName,  LosAngle, neverFail=.true. )
      ! call L1BGeoLocation ( filedatabase, 'sc/OrbIncl       ', instrumentModuleName,  OrbIncl, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/SolarTime    ', instrumentModuleName,  SolarTime, neverFail=.true. )
      call L1BGeoLocation ( filedatabase, 'GHz/SolarZenith  ', instrumentModuleName,  SolarZenith, neverFail=.true. )
      call output ( "Comparing start, end geolocations all chunks: " , &
        & advance='yes')
      do chunk = 1, size ( chunks )
        call outputnamedValue( 'Chunk', chunk )
        if ( firstHGrid(chunk)%noProfs < 1 ) then
          call output ( '(No profiles in this chunk', advance='yes' )
        else
          if ( all ( [ associated(MAFStartTimeTAI), associated(GeodAngle), &
                     & associated(GeodAlt), associated(GeodLat), &
                     & associated(SolarTime) ] ) ) &
            & call CompareWithChunk( chunks(chunk), &
                & chunks(min(size ( chunks ), chunk+1)), firstHGrid(chunk), &
                & MAFStartTimeTAI, GeodAngle, GeodAlt, GeodLat, SolarTime )
        end if
      end do
      call deAllocate_Test ( MAFStartTimeTAI, 'MAFStartTimeTAI', ModuleName )
      ! call deAllocate_Test ( OrbIncl        , 'OrbIncl        ', ModuleName )
      call deAllocate_Test ( GeodAngle      , 'GHz/GeodAngle  ', ModuleName )
      call deAllocate_Test ( GeodAlt        , 'GHz/GeodAlt    ', ModuleName )
      call deAllocate_Test ( GeodLat        , 'GHz/GeodLat    ', ModuleName )
      call deAllocate_Test ( GeodLat        , 'GHz/GeodLat    ', ModuleName )
      call deAllocate_Test ( Lon            , 'GHz/Lon        ', ModuleName )
      call deAllocate_Test ( LosAngle       , 'GHz/LosAngle   ', ModuleName )
      call deAllocate_Test ( SolarTime      , 'GHz/SolarTime  ', ModuleName )
      call deAllocate_Test ( SolarZenith    , 'GHz/SolarZenith', ModuleName )
    end if

    ! Before destroying anything, 
    ! let's reveal which any chunks have no profiles.
    ! We might need this info if we later choose to run just a few, select chunks
    if ( verboser ) call Dump( firstHGrid%noProfs, 'Num of Profiles in each chunk' )
    if ( any(firstHGrid%noProfs == 0) ) then
      call FindAll ( firstHGrid%noProfs, 0, ChunksWithoutProfiles, how_many )
      call Dump( ChunksWithoutProfiles(1:How_many), 'Chunks with no profiles' )
    endif
    do chunk = 1, size ( chunks )
      if ( verbose ) call Dump ( firstHGrid(chunk), Details=1 )
      call DestroyHGridContents ( firstHGrid(chunk) )
    enddo
    call deAllocate_Test ( LowerOverlaps, 'LowerOverlaps', ModuleName )
    if ( specialDumpFile /= ' ' ) call revertOutput
    call trace_end ( "ComputeAllHGridOffsets", &
      & cond=toggle(gen) )
    ! stop
  end subroutine ComputeAllHGridOffsets

  ! ------------------------------  ComputeNextChunksHGridOffsets  -----
  ! This routine is CURRENTLY NOT USED!
  subroutine ComputeNextChunksHGridOffsets ( chunks, chunkNo, hGrids )
    use Allocate_Deallocate, only: Allocate_Test
    use Chunks_m, only: MLSChunk_T
    use HGridsDatabase, only: HGrids_T
    ! Dummy arguments
    type(MLSChunk_T), dimension(:), intent(inout) :: CHUNKS
    integer, intent(in) :: CHUNKNO
    type(HGrids_T), dimension(:), intent(in) :: HGrids
    integer :: I
    ! Executable code
    ! Do nothing if this is the last chunk
    if ( chunkNo == size(chunks) ) return
    ! Allocate the array
    call Allocate_Test ( chunks(chunkNo+1)%hGridOffsets, size(hGrids), &
      & 'chunks(?)%hGridOffsets', ModuleName )
    ! Compute our number of output instances
    do i = 1, size(hGrids)
      chunks(chunkNo+1)%hGridOffsets = hGrids(i)%the_hGrid%noProfs - &
        & hGrids(i)%the_hGrid%noProfsLowerOverlap - &
        & hGrids(i)%the_hGrid%noProfsUpperOverlap
    end do
    ! For later chunks, add on the accumulated previous stuff
    if ( chunkNo > 1 ) &
      & chunks(chunkNo+1)%hGridOffsets = chunks(chunkNo+1)%hGridOffsets + &
      & chunks(chunkNo)%hGridOffsets + 1
  end subroutine ComputeNextChunksHGridOffsets
  
  ! ------------------------------------- PlaceHGridContents --
  subroutine PlaceHGridContents ( HGrid1, HGrid2, offset )
    ! Place the contents of one Hgrid1 inside HGrid2, possibly offset
    use HGridsDatabase, only: HGrid_t
    ! Args
    type(HGrid_T), intent(in)     :: HGrid1
    type(HGrid_T), intent(inout)  :: HGrid2
    integer, optional, intent(in) :: offset ! where in HGrid2 to place HGrid1
    ! Internal variables
    integer :: myOffset
    ! Executable
    myOffset = 0
    if ( present(offset) ) myOffset=offset
    call PlaceArray(HGrid1%phi, HGrid2%phi, myOffset)
    call PlaceArray(HGrid1%geodLat, HGrid2%geodLat, myOffset)
    call PlaceArray(HGrid1%lon, HGrid2%lon, myOffset)
    call PlaceArray(HGrid1%time, HGrid2%time, myOffset)
    call PlaceArray(HGrid1%solarTime, HGrid2%solarTime, myOffset)
    call PlaceArray(HGrid1%solarZenith, HGrid2%solarZenith, myOffset)
    call PlaceArray(HGrid1%losAngle, HGrid2%losAngle, myOffset)
  end subroutine PlaceHGridContents
  
! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )

    use Lexer_Core, only: Print_Source
    use Output_M, only: Output
    use Tree, only: Where_At => Where

    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( where_at(where) )
    call output ( ': ' )
    select case ( code )
    case ( noFraction )
      call output ( "TYPE = FRACTIONAL but no fraction is specified", &
        & advance='yes' )
    case ( noHeight )
      call output ( "TYPE = HEIGHT but no height is specified", advance='yes' )
    case ( noL1Bfiles )
      call output ( "This type of hgrid needs an l1boa file", advance='yes' )
    case ( noMIF )
      call output ( "TYPE = FIXED but no MIF is specified", advance='yes' )
    case ( noModule )
      call output ( "Instrument module must be specified", advance='yes' )
    case ( noPolygon )
      call output ( "A polygon has not been defined", advance='yes' )
    case ( noSpacingOrigin )
      call output ( "TYPE = Regular but no spacing and/or origin is specified", &
        & advance='yes' )
    case ( badTime )
      call output ( "Bad information given for date in explicit hGrid", &
        & advance='yes' )
    end select
  end subroutine ANNOUNCE_ERROR
    

  ! ---------------------------------------------  CompareWithChunk  -----
  subroutine CompareWithChunk ( chunk, nextChunk, hgrid, &
    & MAFStartTimeTAI, GeodAngle, GeodAlt, GeodLat, SolarTime )

     use Chunks_M, only: MLSChunk_T ! , Dump
     use Dates_Module, only: Gethid
     use HGridsDatabase, only: HGrid_T
     use HighOutput, only: BlanksToColumn, OutputNamedValue
     ! use Machine, only: Crash_Burn
     use MLSKinds, only: Rk => R8
     use Output_M, only: Blanks, NewLine, Output

    ! Args
    type (MLSChunk_T), intent(in)      :: chunk
    type (MLSChunk_T), intent(in)      :: nextChunk
    type (HGrid_T), intent(in)         :: hgrid
    real(rk), dimension(:), pointer    :: MAFStartTimeTAI, GeodAngle, &
      &                                   GeodAlt, GeodLat, SolarTime
    ! Local variables
    integer :: firstMAF
    integer :: lastMAF
    integer :: lowMAF
    integer :: highMAF
    integer :: firstProfile
    integer :: firstProfIndx
    integer :: lastProfile
    integer :: lastProfIndx
    real(rk) :: firstVal
    real(rk) :: lastVal
    ! Executable
    ! call Dump( chunk )
    ! call Dump( nextChunk )
    call output( 'First, last MAFs (then times)', advance='no' )
    call blanks( 20 )
    call output ( (/chunk%firstMAFIndex, chunk%lastMAFIndex/), advance='no' )
    call gethid( MAFStartTimeTAI ( chunk%firstMAFIndex + 1 ), leapsec=.true. , hid=firstVal )
    call gethid( MAFStartTimeTAI ( chunk%lastMAFIndex + 1 ), leapsec=.true. , hid=lastVal  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call newLine
!     call crash_burn
    firstMAF = chunk%firstMAFIndex + chunk%noMAFsLowerOverlap
    lastMAF = chunk%lastMAFIndex - chunk%noMAFsUpperOverlap
    firstProfile = 1 + chunk%HGridOffsets(1)
    if ( chunk%lastMAFIndex <  nextChunk%lastMAFIndex ) then
      lastProfile = nextChunk%HGridOffsets(1)
    else
      lastProfile = firstProfile + lastMAF - firstMAF
    endif
    call output( 'First, last non-overlap MAFs', advance='no' )
    call blanks( 6 )
    call output ( (/firstMAF, lastMAF/), advance='yes' )
    
    call output( 'First, last profiles', advance='no' )
    call blanks( 6 )
    call output ( (/firstProfile, lastProfile/), advance='yes' )
    if ( lastProfile > chunk%hGridTotals(1) ) then
      call output ( 'Uh-oh, last profile beyond chunk grand total', advance='yes' )
    elseif ( lastProfile < firstProfile ) then
      call output ( '(No profiles in this chunk)', advance='yes' )
      return
    endif
    ! Must reset profile first, last prfile numbers
    ! before using them as indexes into HGrid
    firstProfIndx = firstProfile - chunk%HGridOffsets(1)
    lastProfIndx = lastProfile - chunk%HGridOffsets(1)
    call outputnamedValue ( 'numMAFs (non-overlap)', lastMAF - firstMAF + 1 )
    call outputnamedValue ( 'num Profiles', lastProfile - firstProfile + 1 )
    call outputnamedValue ( 'noProfs', hGrid%noProfs )
    call outputnamedValue ( 'net Profs', hGrid%noProfs &
      & - hGrid%noProfsLowerOverlap - hGrid%noProfsUpperOverlap )
    if ( lastProfIndx > size(hGrid%time, 2) ) then
      call output ( 'Uh-oh, last profile beyond hgrid upper bound', advance='yes' )
      lastProfIndx = size(hGrid%time, 2)
    endif
    lowMAF = hgrid%maf(1) + chunk%firstMAFIndex - 1
    highMAF = hgrid%maf(hGrid%noProfs) + chunk%firstMAFIndex - 1
    highMAF = min(highMAF, size(MAFStartTimeTAI) - 1 )
    call output( 'First, last MAFs matching grid (then times)', advance='no' )
    call blanks( 6 )
    call output ( (/lowMAF, highMAF/), advance='no' )
    call gethid( MAFStartTimeTAI ( lowMAF + 1 ), leapsec=.true. , hid=firstVal )
    call gethid( MAFStartTimeTAI ( highMAF + 1 ), leapsec=.true. , hid=lastVal  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call newLine
    
    call blanksToColumn( 24 )
    call output( 'Chunk', advance='no' )
    call blanksToColumn( 48 )
    call output( 'HGrid', advance='yes' )

    call blanksToColumn( 4 )
    call output( 'First', advance='no' )
    call blanksToColumn( 18 )
    call output( 'Last', advance='no' )
    call blanksToColumn( 30 )
    call output( 'First', advance='no' )
    call blanksToColumn( 46 )
    call output( 'Last', advance='yes' )
    ! Time
    call gethid( MAFStartTimeTAI ( firstMAF + 1 ), leapsec=.true. , hid=firstVal )
    call gethid( MAFStartTimeTAI ( lastMAF + 1  ), leapsec=.true. , hid=lastVal  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 28 )
    call gethid( hGrid%time ( 1,firstProfIndx ), leapsec=.true. , hid=firstVal )
    call gethid( hGrid%time ( 1,lastProfIndx  ), leapsec=.true. , hid=lastVal  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 54 )
    call output( 'hours in day', advance='yes' )

    ! GeodAngle
    firstVal = GeodAngle ( firstMAF + 1 )
    lastVal  = GeodAngle ( lastMAF + 1 )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 28 )
    firstVal = hGrid%phi ( 1,firstProfIndx )
    lastVal  = hGrid%phi ( 1,lastProfIndx  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 54 )
    call output( 'geodangle', advance='yes' )

    ! geodLat
    firstVal = geodLat ( firstMAF + 1 )
    lastVal  = geodLat ( lastMAF + 1 )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 28 )
    firstVal = hGrid%geodLat ( 1,firstProfIndx )
    lastVal  = hGrid%geodLat ( 1,lastProfIndx  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 54 )
    call output( 'geodLat', advance='yes' )

    ! SolarTime
    firstVal = SolarTime ( firstMAF + 1 )
    lastVal  = SolarTime ( lastMAF + 1 )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 28 )
    firstVal = hGrid%SolarTime ( 1,firstProfIndx )
    lastVal  = hGrid%SolarTime ( 1,lastProfIndx  )
    call output ( firstVal, format='(f9.4)', advance='no' )
    call blanks ( 3 )
    call output ( lastVal, format='(f9.4)', advance='no' )
    call blanksToColumn( 54 )
    call output( 'SolarTime', advance='yes' )
  end subroutine CompareWithChunk

  subroutine PlaceArray_r4(array1, array2, offset)
    ! place contents of array1 inside array2, possibly offset
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: R4 = kind(0.0e0)
    ! Args
    real(r4), dimension(:,:), intent(in)     :: array1
    real(r4), dimension(:,:), intent(inout)  :: array2
    integer, intent(in) :: offset ! Where inside array2 to place array
    ! Internal variables
    integer :: N
    ! Executable
    N = size(array1, 2)
    if ( size(array1, 1) /= size(array2, 1) ) &
      & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '/PlaceArray_r4', &
      & 'array sizes mismatched in first index' )
    if ( N+offset > size(array2, 2) ) &
      & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '/PlaceArray_r4', &
      & 'array1 too big to place in array2 with given offset' )
    array2(:, 1+offset:N+offset) = array1
  end subroutine PlaceArray_r4

  subroutine PlaceArray_r8(array1, array2, offset)
    ! place contents of array1 inside array2, possibly offset
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: R8 = kind(0.0d0)
    ! Args
    real(r8), dimension(:,:), intent(in)     :: array1
    real(r8), dimension(:,:), intent(inout)  :: array2
    integer, intent(in) :: offset ! Where inside array2 to place array
    ! Internal variables
    integer :: N
    ! Executable
    N = size(array1, 2)
    if ( size(array1, 1) /= size(array2, 1) ) &
      & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '/PlaceArray_r4', &
      & 'array sizes mismatched in first index' )
    if ( N+offset > size(array2, 2) ) &
      & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '/PlaceArray_r4', &
      & 'array1 too big to place in array2 with given offset' )
    array2(:, 1+offset:N+offset) = array1
  end subroutine PlaceArray_r8

  subroutine SayTimeHere ( What )
    use Output_M, only: Blanks, Output
    use Time_M, only: Time_Now
    character(len=*), intent(in) :: What
    call time_now ( t2 )
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - t1), advance = 'no' )
    if ( .true. ) then
      call blanks ( 4, advance = 'no' )
      call output ( "Total = " )
      call output ( dble(t2-t0), advance = 'yes' )
    end if
  end subroutine SayTimeHere

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HGrid
!=============================================================================

!
! $Log$
! Revision 2.157  2021/11/03 23:47:58  pwagner
! maxLowerOverlap and -Upper many now be negative, meaning to trim away profiles from either end
!
! Revision 2.156  2019/11/22 00:38:43  pwagner
! Be less prone to printing lots; show which chunks have no profiles
!
! Revision 2.155  2019/11/14 22:28:32  pwagner
! Now requires more coaxing to trace sectionLoop in ComputeAllHGridOffsets
!
! Revision 2.154  2019/07/09 22:26:19  pwagner
! Reduce unwanted Checks ForCorruptFileDatabase
!
! Revision 2.153  2018/08/06 20:24:15  pwagner
! Dont exit section loop in ComputeAllHGridOffsets w/o trace_end
!
! Revision 2.152  2018/08/03 23:45:43  vsnyder
! Add tracing to several routines.  Improve some debugging output.  Add
! NeverFail=.true. to L1BGeolcation calls in verbose case to prevent crashing
! if irrelevant stuff is attempted to be found.
!
! Revision 2.151  2018/04/19 00:49:06  vsnyder
! Remove USE statements for unused names
!
! Revision 2.150  2018/04/16 22:20:21  pwagner
! More thorough CamelCase in use statements
!
! Revision 2.149  2018/03/02 00:58:17  pwagner
! Reduce non-debug printing
!
! Revision 2.148  2018/02/23 22:20:51  mmadatya
! Updates for Create_QTM_HGrid for ASMLS
!
! Revision 2.147  2018/01/03 01:17:56  pwagner
! Removed disused solving of a quadratic equation
!
! Revision 2.146  2016/11/04 19:36:24  pwagner
! begin transition to sayTime from time_m
!
! Revision 2.145  2016/10/19 00:31:18  pwagner
! Trying to avoid certain crashes; may be signs of deeper problems
!
! Revision 2.144  2016/10/14 00:04:28  pwagner
! Avoid negative hGridOffsets
!
! Revision 2.143  2016/10/01 01:53:28  vsnyder
! Fill hGrid fields after creating QTM
!
! Revision 2.142  2016/10/01 01:37:36  vsnyder
! Make QTM_Tree component of HGrid_t allocatable
!
! Revision 2.141  2016/09/14 20:11:42  vsnyder
! Move writing QTM to QTM_Output module
!
! Revision 2.140  2016/09/03 00:06:01  vsnyder
! Turn off some debug printing
!
! Revision 2.139  2016/09/02 00:48:08  vsnyder
! Use Inclination field of explicit HGrid to compute latitude from phi.
! Simplify method of filling fields of explicit HGrid.  Verify that fields
! of explicit HGrid have the same sizes.  Use Phi_To_Lat_Deg from Geometry.
!
! Revision 2.138  2016/08/23 00:43:34  vsnyder
! Components within or adjacent to the polygon are now within the QTM_Tree_t
! structure instead of the HGrid_t structure.
!
! Revision 2.137  2016/08/17 00:47:52  vsnyder
! Allow QTM resolution to be defined by level, length along meridian, or
! degrees along meridian.  Don't get the instrument module name in
! ComputeAllHGridOffsets unless verbose is set.
!
! Revision 2.136  2016/08/09 21:12:50  pwagner
! Survives encounter with non-satellite data
!
! Revision 2.135  2016/07/28 02:02:26  vsnyder
! Removed unreferenceed USE
!
! Revision 2.134  2016/07/28 01:45:07  vsnyder
! Refactor dump and diff
!
! Revision 2.133  2016/07/27 23:02:20  pwagner
! Works better with Aircraft-borne instrument data
!
! Revision 2.132  2016/07/22 20:04:46  pwagner
! Permit using hdf4 l1b files
!
! Revision 2.131  2016/07/21 20:55:09  pwagner
! Deal smootly when items missing from l1b file
!
! Revision 2.130  2016/05/18 01:37:30  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.129  2016/02/26 02:07:18  vsnyder
! Add QTM support
!
! Revision 2.128  2016/02/12 20:10:53  pwagner
! Better error checking, more complete DestroyHGridGeoLocations
!
! Revision 2.127  2015/10/14 23:23:01  pwagner
! Warn instead of crashing if no profiles in chunk; housekeeping
!
! Revision 2.126  2015/07/14 23:34:22  pwagner
! May easily debug cases when gaps in counterMAF
!
! Revision 2.125  2015/06/19 21:12:44  pwagner
! Sped up Computing HGrid offsets
!
! Revision 2.124  2015/06/03 23:10:51  pwagner
! Prevent certain crashes
!
! Revision 2.123  2015/05/05 16:45:13  pwagner
! Merged changes in branch v4.21
!
! Revision 2.122  2015/04/29 01:16:34  vsnyder
! Cosmetic changes
!
! Revision 2.121  2015/03/28 02:45:56  vsnyder
! Get IsMonotonic from Monotone instead of MLSFillValues.  Save HGrid type.
!
! Revision 2.120  2015/03/10 23:40:18  pwagner
! Revert maf component of HGrid to being relative, not absolute
!
! Revision 2.119  2015/03/10 00:20:53  pwagner
! Many diagnostics added; HGrid%maf now stores absolute maf
!
! Revision 2.118  2014/09/05 01:02:36  vsnyder
! More complete and accurate allocate/deallocate size tracking.
! Add some tracing.  Correct some names to be reported to Allocate_Test.
! Plug a potential memory leak.
!
! Revision 2.117  2014/09/05 00:49:06  vsnyder
! EmpiricalGeometry.f90 -- Wrong comment
!
! Revision 2.116  2014/08/06 23:30:29  vsnyder
! Remove CoordIndex, which is set but not referenced
!
! Revision 2.115  2014/08/01 01:45:52  vsnyder
! Remove unreferenced USE names
!
! Revision 2.114  2014/04/24 23:56:43  pwagner
! May set master coordinate in hGrid specification
!
! Revision 2.113  2014/03/18 17:15:09  pwagner
! Can get leapseconds from dates module if run sans toolkit
!
! Revision 2.112  2014/03/01 03:10:56  vsnyder
! Move units checking to init_tables_module
!
! Revision 2.111  2014/01/11 01:44:18  vsnyder
! Decruftification
!
! Revision 2.110  2013/12/12 02:11:26  vsnyder
! Use iterator to handle variables, and IF and SELECT constructs
!
! Revision 2.109  2013/10/01 22:17:51  pwagner
! Added maf component to HGrid_T
!
! Revision 2.108  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.107  2013/08/30 02:45:41  vsnyder
! Revise calls to trace_begin and trace_end
!
! Revision 2.106  2013/08/21 00:25:04  pwagner
! Removed a debugging remnant dumped when overriding geolocations
!
! Revision 2.105  2013/08/17 00:20:25  pwagner
! Fixed HGrids may have geolocations overridden
!
! Revision 2.104  2013/08/13 01:27:42  vsnyder
! Get kind type parameters from MLSKinds instead of MLSCommon
!
! Revision 2.103  2013/08/13 00:58:54  vsnyder
! Move SolveQuadratic into MLSNumerics
!
! Revision 2.102  2013/06/12 02:37:14  vsnyder
! Cruft removal
!
! Revision 2.101  2012/04/25 20:32:24  pwagner
! Inserting missing profiles after chunk end now an option controlled by 'extendible' field
!
! Revision 2.100  2012/04/20 01:09:22  pwagner
! Regular HGrids no longer drop profiles between chunks
!
! Revision 2.99  2011/06/29 21:54:26  pwagner
! Some cases may safely omit l1b files
!
! Revision 2.98  2011/03/10 21:39:11  pwagner
! May now specify time in explicit hGrids
!
! Revision 2.97  2010/03/23 23:26:03  honghanh
! Add a case where single option is not set, but hgrid%noProfs is 1
! then set overlap fields to 0.
!
! Revision 2.96  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.95  2009/06/16 17:42:16  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.94  2009/05/13 20:41:55  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.93  2008/12/02 23:29:41  pwagner
! Added print to not_used_here
!
! Revision 2.92  2007/10/02 22:40:51  vsnyder
! Increase trace level for CreateHGridFromMLSCFInfo
!
! Revision 2.91  2007/06/21 00:54:07  vsnyder
! Remove tabs, which are not part of the Fortran standard
!
! Revision 2.90  2006/07/12 20:42:30  pwagner
! 0 profiles in an HGrid will stop unless computing offsets
!
! Revision 2.89  2006/06/29 21:56:25  pwagner
! Some debugging stuff removed from routine processing
!
! Revision 2.88  2006/05/31 17:48:49  pwagner
! Another bug fix relating to extra profiles
!
! Revision 2.87  2006/04/19 20:48:13  pwagner
! Undid most of the changes regarding extra MAFs; perhaps fixed bugs
!
! Revision 2.86  2006/04/11 23:33:51  pwagner
! Fixed bug which added excess profiles
!
! Revision 2.85  2006/03/07 23:23:28  vsnyder
! Crash gently if there's bo L1BOA file
!
! Revision 2.84  2006/02/21 19:11:30  pwagner
! Some tweaks to where, when to dump
!
! Revision 2.83  2006/02/07 00:56:26  pwagner
! Now allows overlaps after data end time
!
! Revision 2.82  2006/01/10 23:52:11  pwagner
! Fixed segment fault when hdf4 l1boa file
!
! Revision 2.81  2005/12/16 00:07:12  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.80  2005/12/14 01:55:33  pwagner
! Inadvertantly omitted some of the statements from r2.78
!
! Revision 2.79  2005/12/14 01:43:33  pwagner
! Now stores local apparent solar time
!
! Revision 2.78  2005/12/13 22:15:30  livesey
! New approach to the single option in regular hGrids.  Now it uses
! exactly the same approach to that chosen by forward models when
! phiWindow=0 to allow for useful truth in truth out 1D tests.
!
! Revision 2.77  2005/12/13 21:26:02  livesey
! Minor buglet fix for single case in regular hGrid construction.
!
! Revision 2.76  2005/11/15 00:20:18  pwagner
! Should catch error arising from chunks with bad data
!
! Revision 2.75  2005/11/04 18:52:36  pwagner
! Added warning, correction if mif1GeodAngle non-monotonic
!
! Revision 2.74  2005/10/22 00:49:26  pwagner
! Should throw error when l1boa file not in filedatabase
!
! Revision 2.73  2005/10/19 00:06:09  pwagner
! Fixed bug causing DealWithObstructions to corrupt HGrid
!
! Revision 2.72  2005/09/21 23:19:47  pwagner
! Added DealWithObstructions
!
! Revision 2.71  2005/09/14 00:12:33  pwagner
! Uses ChunkDivideConfig%allowPriorOverlaps in calculating offsets
!
! Revision 2.70  2005/08/31 19:41:16  livesey
! Added option to suppress the geometry dump in that first run through
! generating the HGrids that is done by tree walker to assess where each
! direct write places stuff in the output files.
!
! Revision 2.69  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.68  2005/05/31 17:51:17  pwagner
! Began switch from passing file handles to passing MLSFiles
!
! Revision 2.67  2004/12/27 23:05:47  vsnyder
! Remove unreferenced use names
!
! Revision 2.66  2004/08/20 17:58:42  livesey
! Made dontpad=true on most calls.
!
! Revision 2.65  2004/08/16 17:10:26  pwagner
! Passes dontPad option to readL1BData
!
! Revision 2.64  2004/08/05 20:04:48  livesey
! Changed some spline interpolations to linear to allow for more
! stability.
!
! Revision 2.63  2004/07/30 00:16:36  livesey
! Minor update in the 'single' option.
!
! Revision 2.62  2004/06/14 20:00:36  livesey
! A bit more helpful information in dumping and geometry dumping
!
! Revision 2.61  2004/05/19 19:16:10  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.60  2004/03/24 17:55:19  livesey
! Just tidying up after myself
!
! Revision 2.59  2004/03/24 01:30:16  livesey
! Bug fix in hGrid%time
!
! Revision 2.58  2004/03/24 01:03:23  livesey
! Added date option to explicit hGrid
!
! Revision 2.57  2003/08/28 23:52:36  livesey
! Bug fix for computing total size of hGrid, and tidied up the geometry
! dumper.
!
! Revision 2.56  2003/08/11 20:55:20  livesey
! Changed 0 profiles error to warning
!
! Revision 2.55  2003/08/11 18:08:27  livesey
! Added the single option for regular hGrids
!
! Revision 2.54  2003/06/25 22:27:02  livesey
! Fixed an undefined variable
!
! Revision 2.53  2003/06/25 22:05:31  vsnyder
! Don't run off the end of the tree if there's no output section
!
! Revision 2.52  2003/06/24 23:30:30  livesey
! Got ComputeAllHGridOffsets working (on the surface at least)
!
! Revision 2.51  2003/06/20 19:37:06  pwagner
! Quanities now share grids stored separately in databses
!
! Revision 2.50  2003/05/11 01:43:17  livesey
! Made the hgrid switch a little less voluable
!
! Revision 2.49  2003/05/05 23:00:34  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.48  2003/03/07 00:41:24  pwagner
! DeeBug is turned on by switch
!
! Revision 2.47  2003/02/13 19:05:39  vsnyder
! Move USEs from module to procedure scope, cosmetic changes
!
! Revision 2.46.2.2  2003/03/06 19:26:17  vsnyder
! Delete unreferenced USEd name, show hGrid's label on debugging output
!
! Revision 2.46.2.1  2003/02/13 20:36:06  livesey
! Changes merged in from HEAD
!
! Revision 2.47  2003/02/13 19:05:39  vsnyder
! Move USEs from module to procedure scope, cosmetic changes
!
! Revision 2.46  2003/02/07 00:41:32  livesey
! Bug fix (or at least workaround)
!
! Revision 2.45  2003/02/06 23:30:50  livesey
! New approach for explicit hGrids
!
! Revision 2.44  2003/01/06 20:13:46  livesey
! New handling of overlaps
!
! Revision 2.43  2002/12/11 22:17:05  pwagner
! Added error checks on hdf version
!
! Revision 2.42  2002/12/06 01:07:54  pwagner
! A lot of extra debugging output possible
!
! Revision 2.41  2002/12/05 02:21:08  livesey
! Cosmetic changes
!
! Revision 2.40  2002/11/22 12:20:42  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.39  2002/11/13 01:07:04  pwagner
! Actually reads hdf5 radiances
!
! Revision 2.38  2002/10/08 17:36:20  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.37  2002/09/11 17:40:38  livesey
! Bug fix
!
! Revision 2.36  2002/08/09 16:56:37  livesey
! Modified the edge handling to avoid having 'orphaned' profiles beyond
! the edges of the scan range.
!
! Revision 2.35  2002/08/04 16:03:33  mjf
! Added some nullify statements for Sun's rubbish compiler.
!
! Revision 2.34  2002/08/01 17:19:56  livesey
! Flipped geometry dump over the right way.
!
! Revision 2.33  2002/07/17 06:02:50  livesey
! More conservative settings
!
! Revision 2.32  2002/07/01 23:57:06  livesey
! Explicit HGrids now inherit the processing start time as their time.
!
! Revision 2.31  2002/07/01 23:42:42  vsnyder
! Plug a memory leak
!
! Revision 2.30  2002/06/29 06:11:36  livesey
! Typo!
!
! Revision 2.29  2002/06/29 05:55:19  livesey
! Added the geom diagnostic
!
! Revision 2.28  2002/06/18 22:41:26  livesey
! More fixes to do with nasty aspects of regular hGrids on day boundaries
!
! Revision 2.27  2002/05/24 20:56:53  livesey
! Some fixes for cases where the chunk is only 1 MAF long
!
! Revision 2.26  2002/05/24 16:47:39  livesey
! Added some diagnostics
!
! Revision 2.25  2002/05/06 22:31:28  livesey
! Fixed nullify stuff
!
! Revision 2.24  2002/05/06 21:59:28  livesey
! Fixed get_Boolean bug
!
! Revision 2.23  2002/05/06 21:37:40  livesey
! Added forbidOverspill option
!
! Revision 2.22  2001/12/16 00:58:24  livesey
! Working version. Deals with first and last chunks in day properly.
!
! Revision 2.21  2001/12/14 01:43:02  livesey
! Various bug fixes
!
! Revision 2.20  2001/12/10 20:21:36  livesey
! Added code for regular HGrids
!
! Revision 2.19  2001/07/09 20:15:07  livesey
! Fixed an embarassing memory leak I thought I caught before. Also
! changed some allocatables to pointers to let us use Allocate_Deallocate
!
! Revision 2.18  2001/07/06 21:33:23  dwu
! forgot to deallocate variables
!
! Revision 2.17  2001/07/06 18:48:16  dwu
! Add codes to make interpolationFactor functioning
!
! Revision 2.16  2001/05/30 23:53:15  livesey
! For new version of L1BData
!
! Revision 2.15  2001/05/12 00:17:24  livesey
! Brief tidy up of constructing from l2gp.  Tree walker currently prevents
! this however.
!
! Revision 2.14  2001/05/03 20:32:19  vsnyder
! Cosmetic changes
!
! Revision 2.13  2001/04/26 02:44:17  vsnyder
! Moved *_indices declarations from init_tables_module to intrinsic
!
! Revision 2.12  2001/04/24 22:21:05  livesey
! Gave up on latitude stuff
!
! Revision 2.11  2001/04/23 23:25:26  livesey
! Changed l2gpDatabase to pointer
!
! Revision 2.10  2001/04/21 01:24:55  livesey
! New version, tidied up, can create hGrid from L2GP now
!
! Revision 2.9  2001/04/20 23:11:48  livesey
! Added explicit functionality
!
! Revision 2.8  2001/03/02 01:27:06  livesey
! For new MLSSignals
!
! Revision 2.7  2001/02/22 23:44:29  livesey
! Typo
!
! Revision 2.6  2001/02/22 23:43:26  livesey
! Nullified pointer elements of HGrid_T
!
! Revision 2.5  2001/02/21 01:09:24  livesey
! Tidied stuff up a bit
!
! Revision 2.4  2001/02/09 19:30:16  vsnyder
! Move checking for required and duplicate fields to init_tables_module
!
! Revision 2.3  2001/02/09 00:38:22  livesey
! Various updates
!
! Revision 2.2  2001/02/08 01:50:11  vsnyder
! Move duplicate field checking to tree_checker, set by init_tables
!
! Revision 2.1  2000/12/04 23:34:38  vsnyder
! Move more of addItemToDatabase into the include.
!
! Revision 2.0  2000/09/11 19:18:01  ahanzel
! Changing revision to 2.0.
!
! Revision 1.1  2000/09/07 17:36:29  vsnyder
! Initial version 2.0
!
! Revision 1.9  2000/05/17 23:33:51  lungu
! Added dots between MLSInstrumentModuleName and l1bItemName so that is consistent with L1BOA file.
! Added check "if ( ASSOCIATED(database))deallocate(database)" so it doesn't chrash trying to dealocate
! an "empty" database.
!
! Revision 1.8  2000/05/17 18:15:23  livesey
! Finished off interaction with l2cf.
!

