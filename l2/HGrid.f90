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

  implicit none
  private
  public :: CREATEHGRIDFROMMLSCFINFO, COMPUTENEXTCHUNKSHGRIDOFFSETS, &
    & COMPUTEALLHGRIDOFFSETS, DEALWITHOBSTRUCTIONS

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR
  logical, parameter     :: DEEBUG = .false.
  logical, parameter :: DONTPAD = .false.

  interface PlaceArray
    module procedure PlaceArray_r4
    module procedure PlaceArray_r8
  end interface

! Error codes for "announce_error"
  integer, private, parameter :: NoFraction = 1
  integer, private, parameter :: NoHeight = NoFraction + 1
  integer, private, parameter :: NoL1Bfiles = NoHeight + 1
  integer, private, parameter :: NoModule = NoL1Bfiles + 1
  integer, private, parameter :: NoMIF = NoModule + 1
  integer, private, parameter :: NoSpacingOrigin = NoMIF + 1
  integer, private, parameter :: BadTime = NoSpacingOrigin + 1

contains ! =====     Public Procedures     =============================

  ! -----------------------------------  CreateHGridFromMLSCFInfo  -----
  type(hGrid_T) function CreateHGridFromMLSCFInfo &
    & ( name, root, filedatabase, l2gpDatabase, processingRange, chunk, &
    & onlyComputingOffsets ) result ( hGrid )

    use CHUNKS_M, only: MLSCHUNK_T
    use EXPR_M, only: EXPR
    use HGRIDSDATABASE, only: HGRID_T, CREATEEMPTYHGRID, NULLIFYHGRID
    use init_tables_module, only: f_coordinate, f_date, &
      & f_extendible, f_forbidOverspill, f_fraction, f_geodAngle, f_geodLat, &
      & f_height, f_inclination, f_insetOverlaps, f_interpolationFactor, &
      & f_lon, f_losAngle, f_maxLowerOverlap, f_maxUpperOverlap, f_mif, &
      & f_module, f_origin, &
      & f_single, f_solarTime, f_solarZenith, f_sourceL2gp, f_spacing, &
      & f_time, f_type, &
      & field_first, field_last, &
      & l_explicit, l_fixed, l_fractional, l_height, &
      & l_l2gp, l_regular
    use L1BDATA, only: DEALLOCATEL1BDATA, L1BDATA_T, READL1BDATA, &
      & ASSEMBLEL1BQTYNAME
    use L2GPDATA, only: L2GPDATA_T
    use MLSCOMMON, only: MLSFILE_T, NAMELEN, TAI93_RANGE_T
    use MLSFILES, only: GETMLSFILEBYTYPE
    use MLSKINDS, only: RK => R8
    use MLSL2OPTIONS, only: NEED_L1BFILES
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_L1BREAD
    use MLSNUMERICS, only: HUNT
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use MORETREE, only: GET_BOOLEAN
    use STRING_TABLE, only: GET_STRING
    use TOGGLES, only: GEN, LEVELS, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: DECORATION, NSONS, SUB_ROSA, SUBTREE

  ! This routine creates an hGrid based on the user requests.

    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of hGrid subtree
    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (L2GPData_T), pointer, dimension(:) :: L2GPDATABASE
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    ! logical, intent(in), optional :: SUPPRESSGEOMETRYDUMP
    logical, intent(in), optional :: onlyComputingOffsets
    type (MLSChunk_T), intent(in) :: CHUNK ! The chunk

    ! Local variables
    integer :: coordIndex               ! Tree node
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

    type (L1BData_T) :: l1bField ! L1B data

    real(rk) :: incline                 ! Orbital inclination / degrees
    real(rk) :: MINTIME, MAXTIME        ! Span for a chunk
    logical :: MYSUPPRESSGEOMETRYDUMP
    integer, dimension(2) :: PROFRANGE  ! Profile range
    integer :: A,B                      ! Elements of profile range

    logical :: EXTENDIBLE               ! If set don't lose profiles between chunks
    integer :: FIELD                    ! Subtree index of "field" node
    integer :: FIELD_INDEX              ! F_..., see Init_Tables_Module
    logical :: FORBIDOVERSPILL          ! If set don't allow overlaps beyond L1B
    integer :: GEODANGLENODE            ! Tree node
    integer :: GEODLATNODE            ! Tree node
    logical :: GOT_FIELD(field_first:field_last)
    logical :: INSETOVERLAPS            ! Flag
    integer :: L1BFLAG
    integer :: LONNODE                  ! Tree node
    integer :: LOSANGLENODE             ! Tree node
    integer :: MAXLOWEROVERLAP          ! For some hGrids
    integer :: MAXUPPEROVERLAP          ! For some hGrids
    integer :: MIF                      ! For fixed hGrids
    integer :: NOMAFS                   ! Number of MAFs of L1B data read
    logical :: SINGLE                   ! Just one profile please
    integer :: SON                      ! Son of Root
    integer :: SOLARTIMENODE            ! Tree node
    integer :: SOLARZENITHNODE          ! Tree node
    integer :: TIMENODE                 ! Tree node

    character (len=NameLen) :: InstrumentModuleName

    integer ::  hdfVersion
    character(len=NameLen) :: l1bItemName
    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code
    mySuppressGeometryDump = .false.
    if ( present ( onlyComputingOffsets ) ) &
      & mySuppressGeometryDump = onlyComputingOffsets

    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) .and. NEED_L1BFILES ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
      & // 'CreateHGridFromMLSCFInfo', &
      & "Didn't I warn you about not having an L1BOA file?" )
    hdfversion = L1BFile%HDFVersion

    call nullifyHGrid ( hgrid ) ! for Sun's rubbish compiler
    
    hGrid%name = name
    call trace_begin ( me, "CreateHGridFromMLSCFInfo", root, &
      & cond=toggle(gen) .and. levels(gen) > 1 )

    got_field = .false.
    interpolationFactor = 1.0
    extendible = .false.
    forbidOverspill = .false.
    maxLowerOverlap = -1
    maxUpperOverlap = -1
    insetOverlaps = .false.
    single = .false.
    solarTimeNode = 0
    date = 0
    solarZenithNode = 0
    timeNode = 0
    geodAngleNode = 0
    geodLatNode = 0
    LosAngleNode = 0
    lonNode = 0

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
      case ( f_type )
        hGridType = decoration(subtree(2,son))
      case ( f_coordinate )
        coordIndex = field
        hGrid%masterCoordinate = decoration(subtree(2,son))
      case ( f_module )
        instrumentModule = sub_rosa(subtree(2,son))
        call get_string ( instrumentModule , instrumentModuleName )
      case  ( f_extendible )
        extendible = get_boolean ( fieldValue )
      case  ( f_forbidOverspill )
        forbidOverspill = get_boolean ( fieldValue )
      case ( f_height )
        call expr ( subtree(2,son), expr_units, expr_value )
        height = expr_value(1)
      case ( f_mif )
        call expr ( subtree(2,son), expr_units, expr_value )
        mif = nint(expr_value(1))
      case ( f_maxLowerOverlap )
        call expr ( subtree(2,son), expr_units, expr_value )
        maxLowerOverlap = nint(expr_value(1))
      case ( f_maxUpperOverlap )
        call expr ( subtree(2,son), expr_units, expr_value )
        maxUpperOverlap = nint(expr_value(1))
      case ( f_fraction )
        call expr ( subtree(2,son), expr_units, expr_value )
        fraction = expr_value(1)
      case ( f_insetOverlaps )
        insetOverlaps = get_boolean ( fieldValue )
      case ( f_interpolationFactor )
        call expr ( subtree(2,son), expr_units, expr_value )
        interpolationFactor = expr_value(1)
      case ( f_spacing )
        call expr ( subtree(2,son), expr_units, expr_value )
        spacing = expr_value(1)
      case ( f_origin )
        call expr ( subtree(2,son), expr_units, expr_value )
        origin = expr_value(1)
      case ( f_geodAngle )
        geodAngleNode = son
      case ( f_geodLat )
        geodLatNode = son
      case ( f_losAngle )
        losAngleNode = son
      case ( f_lon )
        lonNode = son
      case ( f_single )
        single = get_boolean ( fieldValue )
      case ( f_solarTime )
        solarTimeNode = son
      case ( f_solarZenith )
        solarZenithNode = son
      case ( f_time )
        timeNode = son
      case ( f_date )
        date = sub_rosa(subtree(2,son))
      case ( f_sourceL2gp )
        l2gp => l2gpDatabase(decoration(decoration(subtree(2,son))))
      case ( f_inclination )
        call expr (subtree ( 2, son), expr_units, expr_value )
        incline = expr_value(1) !??? Never used
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    select case (hGridType)

    case ( l_height, l_fractional, l_fixed ) ! --Fixed, Fractional or Height --
      if (.not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      elseif (.not. NEED_L1BFILES ) then
        call announce_error ( root, NoL1BFILES )
      else
        ! 1st--did we set any explicit geolocations?
        if ( any ( got_field ( (/ f_geodAngle, f_geodLat, f_losAngle, f_lon, &
          & f_solarTime, f_solarZenith, f_time /) ) ) ) then
          call CreateExplicitHGrid ( son, date, geodAngleNode, geodLatNode, &
            & solarTimeNode, solarZenithNode, lonNode, losAngleNode, &
            & processingRange%startTime, timeNode, hGrid )
        endif
        call CreateMIFBasedHGrids ( filedatabase, hGridType, chunk, &
          & got_field, root, height, fraction, interpolationFactor, &
          & instrumentModuleName, mif, maxLowerOverlap, maxUpperOverlap, hGrid )
      end if

    case ( l_explicit ) ! ----------------- Explicit ------------------
      call CreateExplicitHGrid ( son, date, geodAngleNode, geodLatNode, &
        & solarTimeNode, solarZenithNode, lonNode, losAngleNode, &
        & processingRange%startTime, timeNode, hGrid )

    case ( l_regular ) ! ----------------------- Regular --------------
      if (.not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      else if ( .not. all(got_field((/f_spacing, f_origin/)))) then
        call announce_error ( root, NoSpacingOrigin )
      elseif (.not. NEED_L1BFILES ) then
        call announce_error ( root, NoL1BFILES )
      else
        call CreateRegularHGrid ( filedatabase, processingRange, chunk, &
          & spacing, origin, trim(instrumentModuleName), extendible, &
          & forbidOverspill, &
          & maxLowerOverlap, maxUpperOverlap, insetOverlaps, single, hGrid, &
          & onlyComputingOffsets )
      end if

    case ( l_l2gp) ! -------------------- L2GP ------------------------
      
      ! Get the time from the l1b file
      l1bItemName = AssembleL1BQtyName ( "MAFStartTimeTAI", hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMAFIndex, &
        & lastMAF=chunk%lastMAFIndex, &
        & dontPad=DONTPAD )
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//"MAFStartTimeTAI" )
      
      minTime = l1bField%dpField(1,1,1)
      maxTime = l1bField%dpField(1,1,noMAFs)

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

      call deallocateL1BData ( l1bField ) ! Avoid memory leaks
    end select
    
    ! Find nearest maf based on time
    l1bItemName = AssembleL1BQtyName ( "MAFStartTimeTAI", hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, l1bField, noMAFs, &
      & l1bFlag, firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex, &
      & dontPad=DONTPAD )
    call Hunt ( l1bField%dpField(1,1,:), hgrid%time(1,:), hgrid%maf, allowTopValue=.true. )

    if ( switchDetail(switches, 'geom') >= 0 .and. .not. mySuppressGeometryDump ) &
      & call DumpChunkHGridGeometry ( hGrid, chunk, &
      & trim(instrumentModuleName), filedatabase )

    if ( error /= 0 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName // '/' &
      & // 'CreateHGridFromMLSCFInfo', &
      & "See ***** above for error message" )
    call trace_end ( "CreateHGridFromMLSCFInfo", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

  end function CreateHGridFromMLSCFInfo

  ! ----------------------------------------  CreateExplicitHGrid  -----
  subroutine CreateExplicitHGrid ( key, date, geodAngleNode, geodLatNode, &
        & solarTimeNode, solarZenithNode, lonNode, losAngleNode, &
        & Time, timeNode, hGrid )

    use DATES_MODULE, only: UTC2TAI93S
    use EXPR_M, only: EXPR
    use HIGHOUTPUT, only: BEVERBOSE, OUTPUTNAMEDVALUE
    use GLOBAL_SETTINGS, only: LEAPSECFILENAME
    use HGRIDSDATABASE, only: CREATEEMPTYHGRID, HGRID_T
    use INIT_TABLES_MODULE, only: PHYQ_ANGLE, PHYQ_DIMENSIONLESS, PHYQ_TIME
    use MLSKINDS, only: RK => R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use SDPTOOLKIT, only: MLS_UTCTOTAI
    use STRING_TABLE, only: GET_STRING
    use TREE, only: NSONS, SUBTREE

    ! dummy arguments
    integer, intent(in) :: KEY          ! Tree node
    integer, intent(in) :: DATE          ! Date if any
    integer, intent(in) :: GEODANGLENODE ! Geod angle if any
    integer, intent(in) :: GEODLATNODE   ! Geod latitude if any
    integer, intent(in) :: LONNODE       ! longitude if any
    integer, intent(in) :: LOSANGLENODE  ! los angle if any
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
    integer :: PARAM                    ! Loop counter
    integer :: PROF                     ! Loop counter
    integer :: RETURNSTATUS             ! Flag
    integer :: UNITS                    ! Units
    real(rk), dimension(:,:), pointer :: VALUES
    logical :: verbose
    ! The following mysterious integers allow us
    ! to engage in some contemptible trickery to fill
    ! the HGrid's geolocation fields--a better
    ! method would call a well-written procedure, passing it appropriate
    ! args for each field
    ! However, we inherited this piece-o-work and just
    ! expanded it include three more fields
    integer, parameter :: NUMPARAMS = 7
    integer, parameter :: PHIPARAM = 1
    integer, parameter :: SOLARTIMEPARAM = PHIPARAM + 1
    integer, parameter :: SOLARZENITHPARAM = SOLARTIMEPARAM + 1
    integer, parameter :: TIMEPARAM = SOLARZENITHPARAM + 1
    integer, parameter :: LONPARAM = TIMEPARAM + 1
    integer, parameter :: LOSANGLEPARAM = LONPARAM + 1
    integer, parameter :: GEODLATPARAM = LOSANGLEPARAM + 1

    ! Executable code
    verbose = BeVerbose( 'hgrid', -1 )
    ! 1st--try to get the number of instances by rude count
    noProfs = 0
    if ( geodAngleNode /= 0 ) then
      noProfs = nsons ( geodAngleNode ) - 1
    elseif ( geodLatNode /= 0 ) then
      noProfs = nsons ( geodLatNode ) - 1
    elseif ( solarTimeNode /= 0 ) then
      noProfs = nsons ( solarTimeNode ) - 1
    elseif ( solarZenithNode /= 0 ) then
      noProfs = nsons ( solarZenithNode ) - 1
    elseif ( lonNode /= 0 ) then
      noProfs = nsons ( lonNode ) - 1
    elseif ( losAngleNode /= 0 ) then
      noProfs = nsons ( losAngleNode ) - 1
    endif

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
      endif
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
      select case ( param )
      case ( PHIPARAM )
        values => hGrid%phi
        node = geodAngleNode
        units = phyq_angle
      case ( SOLARTIMEPARAM )
        values => hGrid%solarTime
        node = solarTimeNode
        units = phyq_time
      case ( SOLARZENITHPARAM )
        values => hGrid%solarZenith
        node = solarZenithNode
        units = phyq_angle
      case ( TIMEPARAM )
        values => hGrid%Time
        node = TimeNode
        units = phyq_time
      case ( LONPARAM )
        values => hGrid%lon
        node = lonNode
        units = phyq_angle
      case ( LOSANGLEPARAM )
        values => hGrid%losAngle
        node = losAngleNode
        units = phyq_angle
      case ( GEODLATPARAM )
        values => hGrid%GeodLat
        node = geodLatNode
        units = phyq_angle
      end select
      if ( node /= 0 ) then
        do prof = 1, hGrid%noProfs
          call expr ( subtree ( prof+1, node), expr_units, expr_value )
          if ( all ( expr_units(1) /= (/ phyq_dimensionless, units /) ) ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Invalid units for explicit hGrid' )
          values(1,prof) = expr_value(1)
        end do
      end if
    end do
    ! Make the latitude the same as the geod angle
    ! This is a bit of a hack, but it should be OK in all cases.
    hGrid%geodLat = hGrid%phi

  end subroutine CreateExplicitHGrid

  ! ---------------------------------------  CreateMIFBasedHGrids  -----
  subroutine CreateMIFBasedHGrids ( filedatabase, hGridType, &
    & chunk, got_field, root, height, fraction, interpolationFactor,&
    & instrumentModuleName, mif, maxLowerOverlap, maxUpperOverlap, hGrid )
    ! This is part of ConstructHGridFromMLSCFInfo

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use CHUNKS_M, only: MLSCHUNK_T
    use DUMP_0, only: DUMP
    use HGRIDSDATABASE, only: CREATEEMPTYHGRID, DUMP, HGRID_T, TRIMHGRID
    use INIT_TABLES_MODULE, only: F_FRACTION, F_GEODANGLE, F_GEODLAT, F_HEIGHT, &
      & F_LON, F_LOSANGLE, F_MIF, F_TIME, &
      & F_SOLARTIME, F_SOLARZENITH, L_FIXED, L_FRACTIONAL, L_HEIGHT, L_MIF
    use L1BDATA, only: DEALLOCATEL1BDATA, L1BDATA_T, READL1BDATA, &
      & ASSEMBLEL1BQTYNAME
    use MLSCOMMON, only: MLSFILE_T, NAMELEN
    use MLSFILES, only: GETMLSFILEBYTYPE
    use MLSKINDS, only: RK => R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_L1BREAD
    use MLSNUMERICS, only: HUNT, INTERPOLATEVALUES

    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    integer, intent(in)               :: HGRIDTYPE
    type (MLSChunk_T), intent(in)     :: CHUNK
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

    integer :: L1BFLAG                  ! Flag
    integer :: L1BITEM                  ! Loop counter
    integer :: MAF                      ! Loop counter etc.
    real(rk) :: MinAngle, MaxAngle
    logical :: MissingOK
    integer :: NOMAFS                   ! Dimension
    real(rk), dimension(:,:,:), pointer :: TpGeodAlt, TpGeodAngle

    ! MIFs it would choose in the non over/undersampled case
    real(rk), dimension(:), pointer :: defaultField, interpolatedField

    type (L1BData_T) :: l1bField ! L1B data
    integer, dimension(:), pointer :: defaultMIFs
    real(rk), dimension(:), pointer :: defaultIndex
    real(rk), dimension(:), pointer :: interpolatedIndex
    integer ::  hdfVersion
    character (len=NameLen) :: L1BItemName
    type (MLSFile_T), pointer             :: L1BFile

    ! Executable code

      L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
      hdfversion = L1BFile%HDFVersion
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
    endif
    
    hGrid%noProfsLowerOverlap = 0
    hGrid%noProfsUpperOverlap = 0
    
    ! Setup some arrays
    call allocate_test ( defaultField, noMAFs, 'defaultField', ModuleName )
    call allocate_test ( interpolatedField, hGrid%noProfs, 'defaultField', ModuleName )
    call allocate_test ( defaultIndex, noMAFs, 'defaultField', ModuleName )
    call allocate_test ( interpolatedIndex, hGrid%noProfs, 'defaultField', ModuleName )

    ! Now we go through all the important geolocation quantities, read them
    ! in, interpolate them if required and store the result in the hGrid
    
    do l1bItem = 1, NoL1BItemsToRead
      ! Get the name of the item to read
      l1bItemName = l1bItemNames(l1bItem)
      if ( l1bItem >= firstModularItem ) l1bItemName = &
        & trim(instrumentModuleName)//"."//l1bItemName
      
      ! Read it from the l1boa file
      l1bItemName = AssembleL1BQtyName ( l1bItemName, hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, l1bField,noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex, &
        & neverfail=(MissingOK .or. l1bItem == l1b_tplosangle), dontPad=.true. )
      if ( l1bFlag==-1 .and. .not. MissingOK) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
      if ( deebug .and. index(l1bItemName, 'MAFStartTimeTAI') > 0 ) then
        call dump(l1bField%DpField(1,1,:), 'MAFStartTimeTAI (before interpolating)')
      end if
      
      if ( l1bItem==1 ) then       ! do something special for time
        do maf = 1, noMAFs
          defaultField(maf) = l1bField%dpField(1,1,maf) + &
            & (defaultMIFs(maf)-1)*sixth
        end do
      else                         ! Otherwise this is fairly easy.
        do maf = 1, noMAFs
          defaultField(maf) = l1bField%dpField(1,defaultMIFs(maf),maf)
        end do
      end if

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
    end do
    
    call Deallocate_test ( defaultMIFs, 'defaultMIFs', ModuleName )
    call Deallocate_test ( defaultField, 'defaultField', ModuleName )
    call Deallocate_test ( defaultIndex, 'defaultIndex', ModuleName )
    call Deallocate_test ( interpolatedField, 'interpolatedField', ModuleName )
    call Deallocate_test ( interpolatedIndex, 'interpolatedIndex', ModuleName )
    
    ! ??? This calculation may need attention! ***
    hGrid%noProfsLowerOverlap = &
      & nint(chunk%noMAFsLowerOverlap*interpolationFactor)
    hGrid%noProfsUpperOverlap = &
      & nint(chunk%noMAFsUpperOverlap*interpolationFactor)

    ! Now deal with the user requests
    if ( maxLowerOverlap >= 0 .and. &
      & ( hGrid%noProfsLowerOverlap > maxLowerOverlap ) ) &
      call TrimHGrid ( hGrid, -1, hGrid%noProfsLowerOverlap - maxLowerOverlap )
    if ( maxUpperOverlap >= 0 .and. &
      & ( hGrid%noProfsUpperOverlap > maxUpperOverlap ) ) &
      call TrimHGrid ( hGrid, -1, hGrid%noProfsUpperOverlap - maxUpperOverlap )
    
  end subroutine CreateMIFBasedHGrids

  ! -----------------------------------------  CreateRegularHGrid  -----
  subroutine CreateRegularHGrid ( filedatabase, processingRange, chunk, &
    & spacing, origin, instrumentModuleName, extendible, forbidOverspill, &
    & maxLowerOverlap, maxUpperOverlap, insetOverlaps, single, hGrid, &
    & onlyComputingOffsets )

    ! Creates an HGrid with coordinates laid out in
    ! a regular spacing w.r.t. master coordinate, phi, aka Geodetic Angle

    ! Depends on an approximation to the Earth's shape (Empirical Geometry)
    ! and interpolation
    
    ! With older l1b files (pre v2.0), some coordinates
    ! (solar time, solar zenith angle) are mean rather than apparent local
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use CHUNKDIVIDE_M, only: CHUNKDIVIDECONFIG
    use CHUNKS_M, only: MLSCHUNK_T, DUMP
    use CONSTANTS, only: DEG2RAD, RAD2DEG
    use DATES_MODULE, only: UTC_TO_TIME
    use DUMP_0, only: DIFF, DUMP, selfdiff
    use EMPIRICALGEOMETRY, only: EMPIRICALLONGITUDE, CHOOSEOPTIMUMLON0
    use HGRIDSDATABASE, only: CREATEEMPTYHGRID, HGRID_T, TRIMHGRID, FINDCLOSESTMATCH
    use HIGHOUTPUT, only: OUTPUTNAMEDVALUE
    use L1BDATA, only: DEALLOCATEL1BDATA, L1BDATA_T, READL1BDATA, &
      & ASSEMBLEL1BQTYNAME
    use MLSCOMMON, only: MLSFILE_T, NAMELEN, TAI93_RANGE_T
    use MLSFILES, only: HDFVERSION_5, DUMP, GETMLSFILEBYTYPE
    use MLSFILLVALUES, only: ISMONOTONIC, MONOTONIZE
    use MLSHDF5, only: ISHDF5ATTRIBUTEINFILE
    use MLSKINDS, only: RK => R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSNUMERICS, only: HUNT, INTERPOLATEVALUES, SOLVEQUADRATIC
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use MLSSTRINGS, only: HHMMSS_VALUE
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING
    use TOGGLES, only: SWITCHES

    type (MLSFile_T), dimension(:), pointer ::     FILEDATABASE
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (MLSChunk_T), intent(in) :: CHUNK
    real(rk), intent(in) :: SPACING
    real(rk), intent(in) :: ORIGIN
    character (len=*), intent(in) :: INSTRUMENTMODULENAME
    logical, intent(in) :: EXTENDIBLE
    logical, intent(in) :: FORBIDOVERSPILL
    integer, intent(in) :: MAXLOWEROVERLAP
    integer, intent(in) :: MAXUPPEROVERLAP
    logical, intent(in) :: INSETOVERLAPS
    logical, intent(in) :: SINGLE
    type (HGrid_T), intent(inout) :: HGRID ! Needs inout as name set by caller
    logical, intent(in), optional :: onlyComputingOffsets

    ! Local variables/parameters
    real(rk), dimension(:,:), pointer :: ALLGEODANGLE ! For every mif
    real(rk), parameter :: SECONDSINDAY = 24*60*60
    ! Note this next one is ONLY NEEDED for the case where we have only
    ! one MAF in the chunk
    real(rk), parameter :: ORBITALPERIOD = 98.8418*60.0

    logical :: DEEBUGHERE
    real(rk) :: DELTA                   ! A change in angle
    integer :: EXTRA                    ! How many profiles over 1 are we
    real(rk) :: FIRST                   ! First point in of hGrid
    integer :: FIRSTPROFINRUN           ! Index of first profile in processing time
    integer :: FLAG                     ! From ReadL1B
    integer ::  hdfVersion
    integer :: I                        ! Loop counter
    real(rk) :: INCLINE                 ! Mean orbital inclination
    type (L1BData_T) :: L1BFIELD        ! A field read from L1 file
    type (MLSFile_T), pointer             :: L1BFile
    character(len=NameLen) :: l1bItemName
    integer :: LASTPROFINRUN            ! Index of last profile in processing time
    real(rk) :: LAST                    ! Last point in hGrid
    integer :: LEFT                     ! How many profiles to delete from the LHS in single
    real(rk) :: MAXANGLE                ! Largest angle in chunk
    real(rk) :: MAXANGLEFIRSTMAF        ! Gives 'range' of first maf
    real(rk), dimension(:), pointer :: MIF1GEODANGLE ! For first mif
    real(rk) :: MINANGLE                ! Smallest angle in chunk
    real(rk) :: MINANGLELASTMAF         ! Gives 'range' of last maf
    integer :: N                        ! Guess at number of profiles
    integer :: NOMAFS                   ! How many in chunk
    integer :: NOMAFSINFILE             ! From ReadL1B
    integer :: NOMIFS
    real(rk) :: NEXTANGLE               ! First non ovl. MAF for next chunk
    real(rk), dimension(:,:), pointer :: OldMethodValues=> null() ! For comparing with
    logical :: PreV2Oh                  ! Old way (mean) or new (apparent)?
    integer :: RIGHT                     ! How many profiles to delete from the RHS in single
    real(rk), dimension(:), pointer :: TMPANGLE ! A temporary array for the single case
    logical :: verbose
    logical :: warnIfNoProfs

    real(rk) :: a
    real(rk) :: b
    real(rk) :: c
    real(rk) :: ImPart
    real(rk) :: r1
    real(rk) :: r2

    ! Executable code
    nullify(OldMethodValues)
    deebughere = deebug .or. ( switchDetail(switches, 'hgrid') > 0 ) ! e.g., 'hgrid1' 
    warnIfNoProfs = .false.
    verbose = deebughere
    if ( present(onlyComputingOffsets) ) then
      verbose = .not. onlyComputingOffsets
      warnIfNoProfs = onlyComputingOffsets
      deebughere = deebug .or. ( switchDetail(switches, 'hgrid') > 1 ) ! e.g., 'hgrid2' 
      verbose = verbose .or. deebughere
    endif
    if ( deebughere .and. .false.) then
      print *, 'Checking quadratic solution'
      do i=1, 10
        a = 0.25
        b = i - 5
        c = 3
        print *, 'a, b, c: ', a, b, c
        print *, 'b^2 - 4 a c ', b**2 - 4*a*c
        call SolveQuadratic( a, b, c, &
              & r1, r2, imPart )
        print *, 'r1, r2, ImPart: ', r1, r2, imPart
      enddo
      stop
    endif
    L1BFile => GetMLSFileByType(filedatabase, content='l1boa')
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file not in database' )
    hdfversion = L1BFile%HDFVersion
    if ( ( L1BFile%hdfVersion == HDFVERSION_5 ) .and. .false. ) then
      PreV2Oh = .not. IsHDF5AttributeInFile(L1BFile%name, 'BO_name')
    else
      PreV2Oh= .true.
    endif
    ! Setup the empircal geometry estimate of lon0
    ! (it makes sure it's not done twice
    call ChooseOptimumLon0 ( filedatabase, chunk )

    ! First we're going to work out the geodetic angle range
    ! Read an extra MAF if possible, as we may use it later
    ! when computing the overlaps.
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpGeodAngle", &
      & hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFsInFile, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex+1, dontPad=.true. )
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file nullified during read' )
    noMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
    noMIFs = size(l1bField%dpField(1,:,1))
    minAngle = minval ( l1bField%dpField(1,:,1) )
    maxAngleFirstMAF = maxval ( l1bField%dpField(1,:,1) )
    maxAngle = maxval ( l1bField%dpField(1,:,noMAFs) )
    minAngleLastMAF = minval ( l1bField%dpField(1,:,noMAFs) )
    nullify ( mif1GeodAngle, allgeodangle )
    call Allocate_test ( mif1GeodAngle, noMAFs, 'mif1Geodangle', ModuleName )
    call Allocate_test ( allGeodAngle, noMIFs, noMAFs, 'allGeodangle', ModuleName )
    mif1GeodAngle = l1bField%dpField(1,1,1:noMAFs)
    allGeodAngle = l1bField%dpField(1,1:noMIFs,1:noMAFs)
    if ( .not. isMonotonic(mif1GeodAngle) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'mif1GeodAngle is not monotonic--will try anyway' )
      call dump( mif1GeodAngle, 'mif1GeodAngle (before)' )
      call monotonize(mif1GeodAngle)
      call dump( mif1GeodAngle, 'mif1GeodAngle (after)' )
    endif

    ! Get or guess the start of the next chunk.
    i = noMAFs - chunk%noMAFsUpperOverlap + 1
    if ( i < 1 ) then
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
      if ( .not. associated(L1BFile) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Obviously impossible guess where to start next chunk' )
      
    elseif ( i < noMAFsInFile + 1 ) then
      nextAngle = l1bField%dpField(1,1,i)
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
      i = FindClosestMatch ( tmpAngle, l1bField%dpField(1,:,:), 1 )
      first = tmpAngle ( i )
      last = first
      call Deallocate_test ( tmpAngle, 'tmpAngle', ModuleName )
    endif

    ! Done with the L1B data
    call DeallocateL1BData ( l1bField )

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
          endif
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
      call output ( forbidoverspill, advance='yes' )
      call output ( ' allowPriorOverlaps: ' )
      call output ( ChunkDivideConfig%allowPriorOverlaps, advance='yes' )
      call output ( ' allowPostOverlaps: ' )
      call output ( ChunkDivideConfig%allowPostOverlaps, advance='yes' )
    end if

    ! Now fill the other geolocation information, first latitude
    ! Get orbital inclination
    l1bItemName = AssembleL1BQtyName ( "scOrbIncl", hdfVersion, .false. )
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file not in database' )
    ! print *, 'About to try to read ', trim(l1bItemName)
    ! call dump(L1BFile)
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex, &
      & dontPad=.true. )
    if ( .not. associated(L1BFile) ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'l1boa file nullified in read' )
    if ( deebughere ) then
      call dump(l1bField%DpField(1,1,:), l1bItemName)
    end if
    ! Use the average of all the first MIFs to get inclination for chunk
    incline = sum ( l1bField%dpField(1,1,:) ) / noMAFs
    if ( deebughere ) then
      call dump( (/incline/), 'Average inclination')
    end if
    call DeallocateL1BData ( l1bField )
    hGrid%geodLat = rad2deg * asin ( sin( deg2Rad*hGrid%phi ) * &
      & sin ( deg2Rad*incline ) )

    ! Now longitude
    call EmpiricalLongitude ( hGrid%phi(1,:), hGrid%lon(1,:) )

    ! Now time, because this is important to get right, I'm going to put in
    ! special code for the case where the chunk is of length one.
    l1bItemName = AssembleL1BQtyName ( "MAFStartTimeTAI", hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & dontPad=.true. )
    if ( deebughere ) then
      call dump( l1bField%DpField(1,1,:), trim(l1bItemName) // ' (before interpolating)')
      call selfdiff( l1bField%DpField(1,1,:), trim(l1bItemName)  )
    end if
    if ( chunk%firstMAFIndex /= chunk%lastMAFIndex ) then
      if ( deebughere ) call outputNamedValue( 'shape(hGrid%phi)', size(hGrid%phi) )
      if ( deebughere ) call outputNamedValue( 'noProfs', hGrid%noProfs )
      call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
        & hGrid%phi(1,:), hGrid%time(1,:), &
        & method='Linear', extrapolate='Allow' )
    else
      ! Case where only single MAF per chunk, treat it specially
      hGrid%time = l1bField%dpField(1,1,1) + &
        & ( OrbitalPeriod/360.0 ) * ( hGrid%phi - mif1GeodAngle(1) )
    end if
    if ( deebughere ) then
      call dump(hGrid%time(1,:), trim(l1bItemName) // ' (after interpolating)')
      call output('geod angle, before and after interpolating', &
        & advance='yes')
      call dump(mif1GeodAngle, 'before')
      call dump(HGrid%phi(1,:), 'after')
    end if
    call DeallocateL1BData ( l1bField )
      
    ! Solar time
    if ( Prev2Oh ) then
      ! First get fractional day, note this neglects leap seconds.
      ! Perhaps fix this later !???????? NJL. We do have access to the
      ! UTC ascii time field, perhaps we could use that?
      hGrid%solarTime = modulo ( hGrid%time, secondsInDay ) / secondsInDay
      ! Now correct for longitude and convert to hours
      hGrid%solarTime = 24.0 * ( hGrid%solarTime + hGrid%lon/360.0 )
      hGrid%solarTime = modulo ( hGrid%solarTime, 24.0_rk )
    else
      if ( deebughere ) then
        call output ( 'About to attempt to get solar time', advance='yes')
        call output ( 'L1Bfile ' // trim(L1BFile%name), advance='yes')
        call output ( 'Module ' // trim(instrumentModuleName), advance='yes')
        call output ( 'firstMAF ')
        call output ( chunk%firstMAFIndex, advance='yes')
        call output ( 'lastMAF ')
        call output ( chunk%lastMAFIndex, advance='yes')
        call output ( 'shape(solarTime) ')
        call output ( shape(hGrid%solarTime), advance='yes')
      endif
      call GetApparentLocalSolarTime( L1BFile, instrumentModuleName, chunk, &
        & hGrid )
      ! What the heck, let's compute the old way, too, and compare
      
      noMIFs = size(hGrid%solarTime, 1)
      if ( deebughere ) then
        call output('NoMAFs ')
        call output(NoMAFs, advance='yes' )
        call output('NoMIFs ')
        call output(NoMIFs, advance='yes' )
        call output('About to allocate', advance='yes')
        call allocate_test(OldMethodValues, noMIFs, size(hGrid%solarTime, 2), &
          & 'OldMethodValues', ModuleName)
        call output('Managed to allocate', advance='yes')
        OldMethodValues = modulo ( hGrid%time, secondsInDay ) / secondsInDay
        OldMethodValues = 24.0 * ( OldMethodValues + hGrid%lon/360.0 )
        OldMethodValues = modulo ( OldMethodValues, 24.0_rk )
        call output('About to diff solartimes', advance='yes')
        call output('Shape(v1.51) ')
        call output(Shape(OldMethodValues), advance='yes')
        call output('Shape(v2.0) ')
        call output(Shape(hGrid%solarTime), advance='yes')
        call dump(OldMethodValues, 'v1.51')
        call dump(hGrid%solarTime, 'v2.0')
        call diff( OldMethodValues, 'v1.51', hGrid%solarTime, 'v2.0', &
          & options='rs' )
        call deallocate_test(OldMethodValues, 'OldMethodValues', ModuleName)
      endif
    endif

    ! Solar zenith
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpSolarZenith", &
      & hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & dontPad=.true. )
    if ( deebughere ) then
      call dump(l1bField%DpField(1,1,:), trim(l1bItemName) // &
        & ' (before interpolating)')
    end if
    if ( Prev2Oh .or. .true. ) then ! Haven't been able to get the other to work
      ! This we'll have to do with straight interpolation
      call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
        & hGrid%phi(1,:), hGrid%solarZenith(1,:), &
        & method='Linear', extrapolate='Allow' )
      if ( deebughere ) then
        call dump(hGrid%solarZenith(1,:), trim(l1bItemName) // &
          & ' (after interpolating)')
      end if
    else
!       call GetApparentLocalSolarZenith( L1BFile, hGrid, &
!         & l1bField%dpField(1,1,:), chunk )
      call closestApparentLocalSolarZenith( allGeodAngle, l1bField%dpField(1,:,:), &
        & hGrid%phi(1,:), hGrid%solarZenith(1,:) )
      ! What the heck, let's compute the old way, too, and compare
      call allocate_test(OldMethodValues, noMIFs, size(hGrid%solarTime, 2), &
        & 'OldMethodValues', ModuleName)
      call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
        & hGrid%phi(1,:), OldMethodValues(1,:), &
        & method='Linear', extrapolate='Allow' )
      call output('About to diff solar zeniths', advance='yes')
      call output('Shape(v1.51) ')
      call output(Shape(OldMethodValues), advance='yes')
      call output('Shape(v2.0) ')
      call output(Shape(hGrid%solarZenith), advance='yes')
      call dump(OldMethodValues, 'v1.51')
      call dump(hGrid%solarZenith, 'v2.0')
      call diff( OldMethodValues, 'v1.51', hGrid%solarZenith, 'v2.0', &
        & options='rs' )
      call deallocate_test(OldMethodValues, 'OldMethodValues', ModuleName)
      
    endif
    call DeallocateL1BData ( l1bField )

    ! Line of sight angle
    ! This we'll have to do with straight interpolation
    l1bItemName = AssembleL1BQtyName ( instrumentModuleName//".tpLosAngle", &
      & hdfVersion, .false. )
    call ReadL1BData ( L1BFile, l1bItemName, &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
      & dontPad=.true. )
    if ( deebughere ) then
      call dump(l1bField%DpField(1,1,:), trim(l1bItemName) // ' (before interpolating)')
    end if
    call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
      & hGrid%phi(1,:), hGrid%losAngle(1,:), &
      & method='Linear', extrapolate='Allow' )
    if ( deebughere ) then
      call dump(hGrid%losAngle(1,:), trim(l1bItemName) // ' (after interpolating)')
    end if
    call DeallocateL1BData ( l1bField )
    hGrid%losAngle = modulo ( hGrid%losAngle, 360.0_rk )

    ! Now work out how much of this HGrid is overlap
    ! The deal will be the first legitimate profile is the first one who's phi
    ! is above the first non overlapped MAF.
    call Hunt ( hGrid%phi(1,:), mif1GeodAngle(chunk%noMAFsLowerOverlap+1), &
       & hGrid%noProfsLowerOverlap, allowTopValue=.true., allowBelowValue=.true. )

    ! So the hunt returns the index of the last overlapped, which is
    ! the number we want to be in the overlap.
    call Hunt ( hGrid%phi(1,:), nextAngle, &
       & hGrid%noProfsUpperOverlap, allowTopValue=.true., allowBelowValue=.true. )
    ! Here the hunt returns the index of the last non overlapped profile
    ! So we do a subtraction to get the number in the overlap.
    hGrid%noProfsUpperOverlap = hGrid%noProfs - hGrid%noProfsUpperOverlap

    if ( verbose ) then
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
    if ( forbidOverspill ) then
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
       hGrid%noProfsLowerOverlap = 0
       hGrid%noProfsUpperOverlap = 0
       ! Delete all but the 'center' profile
       extra = hGrid%noProfs - 1
       left = extra / 2
       right = extra - left
       if ( left > 0 ) call TrimHGrid ( hGrid, -1, left )
       if ( right > 0 ) call TrimHGrid ( hGrid, 1, right )
    else if (hgrid%noProfs > 1) then
       if ( maxLowerOverlap >= 0 .and. &
          & ( hGrid%noProfsLowerOverlap > maxLowerOverlap ) ) &
          call TrimHGrid ( hGrid, -1, hGrid%noProfsLowerOverlap - maxLowerOverlap )
       if ( maxUpperOverlap >= 0 .and. &
          & ( hGrid%noProfsUpperOverlap > maxUpperOverlap ) ) &
          call TrimHGrid ( hGrid, 1, hGrid%noProfsUpperOverlap - maxUpperOverlap )
    else ! if there is only one profile, then don't care about overlap
       hGrid%noProfsLowerOverlap = 0
       hGrid%noProfsUpperOverlap = 0
    end if

    if ( hGrid%noProfs == 0 ) then
      if ( warnIfNoProfs ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, 'No profiles in hGrid' )
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, 'No profiles in hGrid' )
      end if
    end if

    ! Finally we're done.
    if ( verbose ) then
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

  contains
    subroutine GetApparentLocalSolarTime( L1BFile, instrumentModuleName, chunk, &
        & HGrid )
      character (len=*), intent(in) :: INSTRUMENTMODULENAME
      type (MLSFile_T), pointer     :: L1BFile
      type (MLSChunk_T), intent(in) :: CHUNK
      type(HGrid_T) :: hGrid
      ! real(rk), dimension(:,:)      :: solarTime
      ! Local variables
      type (L1BData_T) :: L1BFIELD        ! A field read from L1 file
      real (rk), dimension(:), pointer :: GMT => null()
      real (rk), dimension(:,:), pointer :: MAFLongitude => null()
      real (rk), dimension(:,:), pointer :: MAFTime => null()
      integer :: maf
      integer :: noMAFs
      integer :: noMIFs
      character(len=1), parameter :: SEPARATOR = ':'
      real (rk), dimension(:), pointer :: SolarLongitude => null()
      real (rk), dimension(:), pointer :: SolarTimeCalibration => null()
      integer :: status
      logical, parameter :: STRICT = .true.
      character(len=27) :: time
      ! Executable
      nullify( GMT, MAFLongitude, MAFTime, &
        & SolarLongitude, SolarTimeCalibration )
      
      ! The 1st formula:

      ! (LST) = (GMT) + ( (LON) - (LON_S) ) / 15

      ! Where (LST)  Local Solar Time
      ! (GMT)   Greenwich mean time (aka utc)
      ! (LON)   Local longitude
      ! (LON_S) Solar longitude (which we will infer as a 1st step)
      l1bItemName = AssembleL1BQtyName ( "Lon", &
        & hdfVersion, .false., instrumentModuleName )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read Lon from l1boa', MLSFile=L1BFile )
      noMIFs = size(l1bField%dpField, 2)
      call allocate_test(GMT, noMAFs, 'GMT', ModuleName)
      call allocate_test(MAFLongitude, noMIFs, noMAFs, 'MAFLongitude', ModuleName)
      call allocate_test(MAFTime, noMIFs, noMAFs, 'MAFTime', ModuleName)
      call allocate_test(SolarLongitude, noMAFs, 'SolarLongitude', ModuleName)
      call allocate_test(SolarTimeCalibration, noMAFs, 'SolarLongitude', ModuleName)
      MAFLongitude = l1bField%dpField(1,:,:)
      call DeallocateL1BData ( l1bField )

      l1bItemName = AssembleL1BQtyName ( "MAFStartTimeTAI", &
        & hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read MAFStartTimeTAI from l1boa', MLSFile=L1BFile )
      MAFTime = l1bField%dpField(1,:,:)

      l1bItemName = AssembleL1BQtyName ( "MAFStartTimeUTC", &
        & hdfVersion, .false. )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read MAFStatTimeUTC from l1boa', MLSFile=L1BFile )
      ! This converts the utc string into time into the day (in hours)
      do maf=1, noMAFs
          call utc_to_time(l1bField%CharField(1,1,maf), status, time, strict)
          GMT(maf) = hhmmss_value ( time, status, separator, strict ) / 3600
      enddo
      call DeallocateL1BData ( l1bField )

      l1bItemName = AssembleL1BQtyName ( "SolarTime", &
        & hdfVersion, .false., instrumentModuleName )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read SolarTime from l1boa', MLSFile=L1BFile )
      SolarTimeCalibration = l1bField%dpField(1,1,:)
      call DeallocateL1BData ( l1bField )

      ! 1st we will find the calibrated solar longitude (in deg)
      ! assuming our 1st formula above holds for all mifs
      SolarLongitude = MAFLongitude(1,:) - &
        & (SolarTimeCalibration-GMT(:))*360./24

      ! 2nd we need the GMT corrected for the regular HGrid offsets from MAFs
      GMT = GMT + ( HGrid%time(1,:) - MAFTime(1,:) ) / 3600
      do maf=1, noMAFs
        HGrid%SolarTime(:,maf) = GMT(maf) + &
          & ( HGrid%Lon(:,maf) - SolarLongitude(maf) ) / 15
      enddo

      ! Now deallocate all the junk we allocated
      call deallocate_test(GMT, 'GMT', ModuleName)
      call deallocate_test(MAFLongitude, 'MAFLongitude', ModuleName)
      call deallocate_test(MAFTime, 'MAFTime', ModuleName)
      call deallocate_test(SolarLongitude, 'SolarLongitude', ModuleName)
      call deallocate_test(SolarTimeCalibration, 'SolarTimeCalibration', ModuleName)
    end subroutine GetApparentLocalSolarTime

    subroutine closestApparentLocalSolarZenith( allGeodAngle, allSolarZenith, &
      & GeodAngle, solarZenith )
    use MLSNUMERICS, only: CLOSESTELEMENT
      ! Another approach--compare against linear interpolation
      ! Dummy arguments
      real(rk), dimension(:,:), intent(in) :: allGeodAngle
      real(rk), dimension(:,:), intent(in) :: allSolarZenith
      real(rk), dimension(:), intent(in)   :: GeodAngle
      real(rk), dimension(:), intent(out)  :: solarZenith
      ! Local variables
      integer :: indices(2)
      integer :: maf
      ! Executable
      do maf = 1, size(GeodAngle)
        call closestElement( GeodAngle(maf), allGeodAngle, indices )
        if ( any(indices < 1) .or. &
          & ( indices(1) > size(allGeodAngle, 1) ) .or. &
          & ( indices(2) > size(allGeodAngle, 2) ) ) then
          call dump(indices, 'indices')
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Unable to find closest element' )
        endif
        solarZenith(maf) = allSolarZenith( indices(1), indices(2) )
      enddo
    end subroutine closestApparentLocalSolarZenith

    subroutine GetApparentLocalSolarZenith(L1BFile, hGrid, solarZenithAtMIF1, &
      & chunk )
      ! Alas, was unable to get this to work properly
      ! There is a web page at
      ! http://www.pcigeomatics.com/cgi-bin/pcihlp/AVHRRAD%7CDETAILS%7CANGLE+GENERATION
      ! that includes a summary of some of the same equations we're trying to
      ! implement here, but starting from the assumed equations
      ! for solar declination and right ascension
      ! In any case, we're doing no better than linear interpolation
      ! for other quantities, so we'll resort to that for solarzenithangle, too
      type (MLSFile_T), pointer     :: L1BFile
      type(HGrid_T) :: hGrid
      real(rk), dimension(:), intent(in) :: solarZenithAtMIF1
      type (MLSChunk_T), intent(in) :: CHUNK
      ! The key equation to be considered is

      ! cos(SZA) = sin(LAT) sin(LAT_S) + cos(LAT) cos(LAT_S) cos( 15 (LST-12) )

      ! where SZA is the solar zenit angle we are seeking
      ! LAT is the local latitude, LAT_S is the sun's latitude, and
      ! LST is the local solar time
      ! We will have to solve for LAT_S and assume it's constant over a maf
      ! Note that all angles are assumed to be in degrees, not radians

      ! Local variables
      real(rk) :: a ! sin(local latitude)
      real(rk) :: b ! cos(local latitude) cos( 15 (local solar time - 12) )
      real(rk) :: c ! cos(SZA calibrating)
      ! The next variables are used if we forego the quadratic solution
      ! and instead rely on solving two equations in two "independent"
      ! variables, s=sin(LAT_S) and c=cos(LAT_S)
      ! So we solve
      ! y1 = cos(SZA[1]) = s sin(Lat[1]) + c Q[1] cos(Lat[1])
      ! y2 = cos(SZA[2]) = s sin(Lat[2]) + c Q[2] cos(Lat[2])
      ! where
      ! Q[i] = cos( 15 (LST[i]-12) )
      ! The solution is of course
      !
      !  s        sin(Lat[1])  Q[1] cos(Lat[1])            y1
      ! ( )  =  (                               ) ^ (-1)  (  )
      !  c        sin(Lat[2])  Q[2] cos(Lat[2])            y2
      !
      ! Rewrite the matrix above as
      ! a11  a12
      !(        )
      ! a21  a22
      ! and let its determinant delta = ( a11 a22 - a12 a21 )
      ! then the inverse will be
      ! a22  -a12
      !(         ) / delta
      ! a21   a11
      real(rk) :: delta ! determinant
      real(rk) :: a11, a12, a21, a22, y1, y2, s
      type (L1BData_T) :: L1BFIELD        ! A field read from L1 file
      real (rk), dimension(:), pointer :: LocalLatitude => null()
      real (rk), dimension(:), pointer :: LocalSolarTime => null()
      integer :: maf
      integer :: mif
      integer :: noMAFs
      integer :: noMIFs
      real(rk) :: r1, r2, imPart, root
      real (rk), dimension(:,:), pointer :: sinSZA => null()
      real (rk), dimension(:), pointer :: SolarLatitudeCalibration => null()
      integer :: status
      logical, parameter :: TESTFIRST = .true. ! Test whether key eq'n true
      real (rk), dimension(:,:), pointer :: testLAT => null()
      real (rk), dimension(:,:), pointer :: testLST => null()
      real (rk), dimension(:,:), pointer :: testSZA => null()
      real (rk), dimension(:,:), pointer :: testSSL => null() ! sin solar lat
      logical, parameter                 :: USEQUADRATIC = .false.
      ! Executable
      ! (We'll put stuff here when we're ready)
      ! call output('Darn! Noone coded this part yet .. tell paw', advance='yes')
      ! hGrid%solarZenith = -999.99

      ! Try to solve for sin(LAT_S) knowing that
      ! sin^2(LAT_S) + cos^2(LAT_S) = 1
      l1bItemName = AssembleL1BQtyName ( "GeodLat", &
        & hdfVersion, .false., instrumentModuleName )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read Lon from l1boa', MLSFile=L1BFile )
      noMIFs = size(l1bField%dpField, 2)
      call allocate_test(LocalLatitude,  noMAFs, 'LocalLatitude', ModuleName)
      call allocate_test(LocalSolarTime, noMAFs, 'LocalSolarTime', ModuleName)
      call allocate_test(SolarLatitudeCalibration, noMAFs, 'SolarLatitudeCalibration', ModuleName)
      LocalLatitude = l1bField%dpField(1,1,:)

      if ( TESTFIRST .or. .not. USEQUADRATIC ) then
        ! Test whether whatever is stored in l1boa file satisfies
        ! key equation
        call allocate_test(testLAT,  noMIFs, noMAFs, 'testLAT', ModuleName)
        call allocate_test(testLST,  noMIFs, noMAFs, 'testLST', ModuleName)
        call allocate_test(testSZA,  noMIFs, noMAFs, 'testSZA', ModuleName)
        call allocate_test(testSSL,  noMIFs, noMAFs, 'testSSL', ModuleName)
        testLAT = l1bField%dpField(1,:,:)
        call dump(testLAT, 'Latitude')
        call DeallocateL1BData ( l1bField )

        l1bItemName = AssembleL1BQtyName ( "SolarTime", &
          & hdfVersion, .false., instrumentModuleName )
        call ReadL1BData ( L1BFile, l1bItemName, &
          & l1bField, noMAFs, status, &
          & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & dontPad=.true. )
        testLST = l1bField%dpField(1,:,:)
        call dump(testLST, 'Solar Time')
        call DeallocateL1BData ( l1bField )

        l1bItemName = AssembleL1BQtyName ( "SolarZenith", &
          & hdfVersion, .false., instrumentModuleName )
        call ReadL1BData ( L1BFile, l1bItemName, &
          & l1bField, noMAFs, status, &
          & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
          & dontPad=.true. )
        testSZA = l1bField%dpField(1,:,:)
        call dump(testSZA, 'Solar Zenith')
        call DeallocateL1BData ( l1bField )

        if ( .not. USEQUADRATIC ) then
          ! Just pick two different points and solve for
          ! sin(LAT_S) and cos(LAT_S) as if they were independent variables
          noMAFs = HGrid%noProfs
          y1 = cos( deg2rad * testSZA(1, 1) )
          y2 = cos( deg2rad * testSZA(1, noMAFs) )
          a11 = sin( deg2rad * testLAT(1, 1) )
          a21 = sin( deg2rad * testLAT(1, noMAFs) )
          a12 = cos( deg2rad * testLAT(1, 1) ) * &
            &   cos( deg2rad * 15 * (testLST(1, 1) - 12) )
          a22 = cos( deg2rad * testLAT(1, 1) ) * &
            &   cos( deg2rad * 15 * (testLST(1, noMAFs) - 12) )
          delta = a11*a22 - a12*a21
          if ( delta == 0._rk)  call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'determinant vanishes in GetApparentLocalSolarZenith', &
          & MLSFile=L1BFile )
          s = ( a22*y1 - a12*y2 ) / delta
          c = ( a21*y1 + a11*y2 ) / delta
          ! How reasonable are these?
          print *, 's, c: ', s, c
          print *, 's^2 + c^2: ', s**2 + c**2
          print *, 'atan(s/c): ', rad2deg*atan2(s, c)
          ! Now with s and c in hand, we can go ahead and apply our formula
          ! (and get answers we hope to be correct)
          HGrid%solarZenith = acos( &
            & s * sin( deg2rad * HGrid%geodLat ) + &
            & c * cos( deg2rad * HGrid%geodLat ) * &
            &     cos( deg2rad * 15 * (HGrid%solarTime - 12) ) &
            & )
          return
        endif
        do maf=1, noMAFs
          do mif=1, noMIFs
            a = sin( deg2rad * testLAT(mif, maf) )
            b = cos( deg2rad * testLAT(mif, maf) ) * &
              & cos( deg2rad * 15 * (testLST(mif, maf) - 12) )
            c = cos( deg2rad * testSZA(mif, maf) )
            call SolveQuadratic( (a**2 + b**2), -2*a*c, c**2 - b**2, &
              & r1, r2, imPart )
            if ( imPart /= 0._rk ) then
              call output('mif, maf: ', advance='no')
              call output(mif, advance='no')
              call output(maf, advance='yes')
              call output('testLAT ', advance='no')
              call output(testLAT(mif, maf), advance='yes')
              call output('testLST ', advance='no')
              call output(testLST(mif, maf), advance='yes')
              call output('testSZA ', advance='no')
              call output(testSZA(mif, maf), advance='yes')
              call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'quadratic solver returned complex roots during test', &
              & MLSFile=L1BFile )
            endif
            if ( (c - a*r1)/b > 0. ) then
              testSSL(mif, maf) = asin( rad2deg * r1 )
            else
              testSSL(mif, maf) = asin( rad2deg * r2 )
            endif
          enddo
        enddo
        call dump(testSSL, 'test solar latitudes')
        call deallocate_test(testLAT, 'testLAT', ModuleName)
        call deallocate_test(testLST, 'testLST', ModuleName)
        call deallocate_test(testSZA, 'testSZA', ModuleName)
        call deallocate_test(testSSL, 'testSSL', ModuleName)
      else
        call DeallocateL1BData ( l1bField )
      endif

      l1bItemName = AssembleL1BQtyName ( "SolarTime", &
        & hdfVersion, .false., instrumentModuleName )
      call ReadL1BData ( L1BFile, l1bItemName, &
        & l1bField, noMAFs, status, &
        & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex, &
        & dontPad=.true. )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Unable to read SolarTime from l1boa', MLSFile=L1BFile )
      LocalSolarTime = l1bField%dpField(1,1,:)
      call DeallocateL1BData ( l1bField )

      ! The 1st thing to do is to compute calibrating LAT_S
      do maf=1, noMAFs
      
        ! Now solve the key equation above for LAT_S given the l1b data
        a = sin( deg2rad * LocalLatitude(maf) )
        b = cos( deg2rad * LocalLatitude(maf) ) * &
          & cos( deg2rad * 15 * (LocalSolarTime(maf) - 12) )
        c = cos( deg2rad * solarZenithAtMIF1(maf) )
        ! b = 0 is a special case--don't need to solve quadratic
        call output('LocalLatitude: ')
        call output(LocalLatitude(maf), advance='yes')
        call output('LocalSolarTime: ')
        call output(LocalSolarTime(maf), advance='yes')
        call output('solarZenithAtMIF1: ')
        call output(solarZenithAtMIF1(maf), advance='yes')
        call output('a, b, c: ')
        call output( (/a, b, c/), advance='yes')
        if ( b == 0._rk ) then
          root = c / a
        else
          call SolveQuadratic( (a**2 + b**2), -2*a*c, c**2 - b**2, &
            & r1, r2, imPart )
          if ( imPart /= 0._rk ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'quadratic solver returned complex roots', MLSFile=L1BFile )
          ! Only one of the roots satisfies (c - a x) / b > 0
          if ( (c - a*r1)/b > 0._rk ) then
            root = r1
          elseif ( (c - a*r2)/b > 0._rk ) then
            root = r2
          else
            call MLSMessage ( MLSMSG_Error, ModuleName, &
              & 'Neither root satisfies (c - a x)/b > 0', MLSFile=L1BFile )
          endif
        endif
        if ( abs(root) > 1._rk ) call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'abs(root) > 1', MLSFile=L1BFile )
        SolarLatitudeCalibration(maf) = rad2deg * asin(root)
      enddo
      call deallocate_test(LocalLatitude, 'LocalLatitude', ModuleName)
      call deallocate_test(LocalSolarTime, 'LocalSolarTime', ModuleName)

      call allocate_test(sinSZA, noMIFs, noMAFs, 'sinSZA', ModuleName)
      ! Now use key equation to find SZA
      do maf=1, noMAFs
        sinSZA(:,maf) = sin( deg2rad * hGrid%geodLat(:,maf) ) * &
          &             sin( deg2rad * SolarLatitudeCalibration(maf) ) &
          &   + &
          &             cos( deg2rad * hGrid%geodLat(:,maf) ) * &
          &             cos( deg2rad * SolarLatitudeCalibration(maf) ) * &
          &             cos( deg2rad * 15 * (hGrid%SolarTime(:,maf) - 12) )
      enddo
      hGrid%solarZenith = rad2deg * asin(sinSZA)
      call deallocate_test(sinSZA, 'sinSZA', ModuleName)
      call deallocate_test(SolarLatitudeCalibration, 'SolarLatitudeCalibration', ModuleName)
    end subroutine GetApparentLocalSolarZenith

  end subroutine CreateRegularHGrid

  ! --------------------------------------- DealWithObstructions -----
  subroutine DealWithObstructions ( HGrid, obstructions, DestroyOld )
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use CHUNKDIVIDE_M, only: OBSTRUCTION_T
    use HGRIDSDATABASE, only: HGRID_T, CREATEEMPTYHGRID, DESTROYHGRIDCONTENTS
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_WARNING
    ! Args
    type (HGRID_T), pointer :: HGRID
    type (Obstruction_T), dimension(:), pointer :: OBSTRUCTIONS
    logical, optional, intent(in) :: DestroyOld
    ! This routine modifies the chunks according to the information
    ! given in the obstructions.

    ! Local variables
    type (HGRID_T), target :: NEWHGrid ! The one we'll create create
    integer :: FIRSTMAF                 ! Index of first MAF in range
    integer :: LASTMAF                  ! Index of last MAF in range
    integer :: MAF                      ! Index of MAF for wall
    logical :: mayDestroyOld
    integer :: newProfile               ! Counter in newHGrid
    integer :: OBSTRUCTION              ! Loop counter
    integer :: PROFILE                  ! Loop counter
    logical, dimension(:), pointer :: obstructed => null()

    ! Executable code
    ! 1st--some short-circuits
    if ( hGrid%noProfs < 1 ) return
    if ( .not. associated(obstructions) ) return
    if ( .not. associated(hGrid%phi) ) return
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
      endif
    enddo
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
      enddo
      if ( hgrid%noProfsLowerOverlap > 0 ) then
        newHGrid%noProfsLowerOverlap = hgrid%noProfsLowerOverlap &
          & - count( obstructed(1:hgrid%noProfsLowerOverlap) )
      endif
      if ( hgrid%noProfsUpperOverlap > 0 ) then
        newHGrid%noProfsUpperOverlap = hgrid%noProfsUpperOverlap &
          & - count( obstructed(hgrid%noProfs-hgrid%noProfsUpperOverlap+1:hgrid%noProfs) )
      endif
      if ( mayDestroyOld ) then
        call DestroyHGridContents(hGrid)
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Destroying old HGrid known to corrupt other HGrid' )
      endif
      hGrid => newHGrid
    endif
    call DeAllocate_Test( obstructed, 'obstructed', ModuleName )
  end subroutine DealWithObstructions

  ! -------------------------------------  DumpChunkHGridGeometry  -----
  subroutine DumpChunkHGridGeometry ( hGrid, chunk, &
    & instrumentModuleName, filedatabase )

    use CHUNKS_M, only: MLSCHUNK_T
    use HGRIDSDATABASE, only: HGRID_T
    use L1BDATA, only: DEALLOCATEL1BDATA, L1BDATA_T, READL1BDATA, &
      & ASSEMBLEL1BQTYNAME
    use MLSCOMMON, only: MLSFILE_T, NAMELEN
    use MLSFILES, only: GETMLSFILEBYTYPE
    use MLSKINDS, only: R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ALLOCATE, &
      & MLSMSG_DEALLOCATE, MLSMSG_ERROR
    use OUTPUT_M, only: OUTPUT
    use STRING_TABLE, only: DISPLAY_STRING

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
    allocate ( text ( noBins, noLines ), stat=flag )
    if ( flag /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'text' )
    text = ' '

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
    deallocate ( text, stat=flag )
    if ( flag /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'text' )
    call DeallocateL1BData ( l1bField )

  end subroutine DumpChunkHGridGeometry

  ! -------------------------------------  ComputeAllHGridOffsets  -----
  subroutine ComputeAllHGridOffsets ( root, first_section, chunks, filedatabase, &
    & l2gpDatabase, processingRange )
    ! This routine goes through the L2CF up to an Output section to accumulate
    ! HGrid sizes from the Construct sections, and through the L1 file to work
    ! out how big each HGrid is going to be
    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use CHUNKS_M, only: MLSCHUNK_T
    use CHUNKDIVIDE_M, only: CHUNKDIVIDECONFIG
    use HGRIDSDATABASE, only: HGRID_T, DESTROYHGRIDCONTENTS, DUMP
    use INIT_TABLES_MODULE, only: Z_CONSTRUCT, S_HGRID, Z_OUTPUT
    use L2GPDATA, only: L2GPDATA_T
    use MLSCOMMON, only: MLSFILE_T, TAI93_RANGE_T
    use MLSL2OPTIONS, only: SPECIALDUMPFILE
    use MORETREE, only: GET_SPEC_ID
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
    use MLSSTRINGLISTS, only: SWITCHDETAIL
    use Next_Tree_Node_m, only: Init_Next_Tree_Node, Next_Tree_Node, &
      & Next_Tree_Node_State
    use OUTPUT_M, only: BLANKS, OUTPUT, REVERTOUTPUT, SWITCHOUTPUT
    use TOGGLES, only: GEN, SWITCHES, TOGGLE
    use TRACE_M, only: TRACE_BEGIN, TRACE_END
    use TREE, only: SUBTREE, NODE_ID, DECORATION
    use TREE_TYPES, only: N_NAMED
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
    integer :: HGRID                    ! Loop counter
    integer :: Me = -1                  ! String index for trace
    integer :: NOHGRIDS                 ! Number of hGrids
    integer :: SON                      ! Tree node
    integer :: sum1, sum2, sum3
    integer :: GSON                     ! son of son
    integer :: KEY                      ! Tree node
    type(HGrid_T) :: DUMMYHGRID         ! A temporary hGrid
    type(next_tree_node_state) :: State1 ! while hunting for Construct sections
    type(next_tree_node_state) :: State2 ! within Construct sections
    integer, dimension(:), pointer :: LowerOverlaps => null()
    ! Executable code
    call trace_begin ( me, "ComputeAllHGridOffsets", root, cond=toggle(gen) )
    if ( specialDumpFile /= ' ' ) &
      & call switchOutput( specialDumpFile, keepOldUnitOpen=.true. )
    ! Slightly clever loop here, the zero run doesn't do anything except
    ! count the number of hGrids.  The other run works out how big each HGrid
    ! is going to be for each chunk.  This is stored in chunks%hGridOffsets.
    ! Finally we accumulate these to get offsets.
    noHGrids = 0
    do chunk = 0, size(chunks)
      hGrid = 1
      ! Loop over all the setions in the l2cf, look for construct sections
      call init_next_tree_node ( state1 )
      sectionLoop: do
        son = next_tree_node ( root, state1, start=first_section, traceLevel=4 )
        if ( son == 0 ) exit
        select case ( decoration ( subtree ( 1, son ) ) )
        case ( z_construct )
          ! Now loop through the construct section and identify the hGrids
          call init_next_tree_node ( state2 )
          do
            gson = next_tree_node ( son, state2, traceLevel=5 )
            if ( gson == 0 ) exit
            if ( node_id(gson) == n_named ) then ! Is spec labeled?
              key = subtree(2,gson)
              if ( get_spec_id(key) == s_hGrid ) then
                ! This is an hGrid definition
                if ( chunk == 0 ) then
                  ! For the 'zeroth' pass just count up the hgrids
                  noHGrids = noHGrids + 1
                else
                  dummyHGrid = CreateHGridFromMLSCFInfo ( 0, key, filedatabase, l2gpDatabase, &
                    & processingRange, chunks(chunk), onlyComputingOffsets=.true. )
                  if ( chunk == 1 ) then
                    LowerOverlaps(hGrid) = dummyHGrid%noProfsLowerOverlap
                    if ( ChunkDivideConfig%allowPriorOverlaps .and. .false. ) then
                      chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                      & dummyHGrid%noProfsUpperOverlap
                    else
                      chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                      & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                    endif
                  elseif ( chunk == size(chunks) ) then
                    if ( ChunkDivideConfig%allowPostOverlaps ) then
                      chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                      & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                    else
                      chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                      & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                    endif
                  else
                    chunks(chunk)%hGridOffsets(hGrid) = dummyHGrid%noProfs - &
                    & dummyHGrid%noProfsLowerOverlap - dummyHGrid%noProfsUpperOverlap
                  endif
                  if ( switchDetail(switches, 'hgrid') >= 0 .and. DEEBUG ) &
                    & call dump(dummyHGrid)
                  call DestroyHGridContents ( dummyHGrid )
                  hGrid = hGrid + 1
                end if
              end if
            else ! Son is n_spec_args
              key = gson
            end if
          end do
        case ( z_output )
          exit sectionLoop
        case default
        end select
      end do sectionLoop

      ! If this is the first time round, we now know how many hGrids there
      ! are, setup our arrays
      if ( chunk == 0 ) then
        do c = 1, size ( chunks )
          call Allocate_Test ( chunks(c)%hGridOffsets, noHGrids, &
            & 'chunks(?)%hGridOffsets', ModuleName )
          call Allocate_Test ( chunks(c)%hGridTotals, noHGrids, &
            & 'chunks(?)%hGridTotals', ModuleName )
          call Allocate_Test ( LowerOverlaps, noHGrids, &
            & 'LowerOverlaps', ModuleName )
        end do
      else
        ! Otherwise, at least check we got the same number of hGrids each chunk
        if ( hGrid /= noHGrids + 1 ) &
          & call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Got a different number of hGrids for different chunks, that makes no sense!' )
      end if
    end do                              ! Chunk loop
    
    if ( switchDetail(switches, 'hgrid') >= 0 .and. DEEBUG ) then
      sum1 = 0
      sum2 = 0
      sum3 = 0
      call output ( "number of profiles in each chunk/grid: " , advance='yes')
      do chunk = 1, size ( chunks )
        call output ( chunk )
        call blanks ( 3 )
        call output ( chunks(chunk)%hGridOffsets, advance='yes' )
        sum1 = sum1 + chunks(chunk)%hGridOffsets(1)
        sum2 = sum2 + chunks(chunk)%hGridOffsets(2)
        sum3 = sum3 + chunks(chunk)%hGridOffsets(3)
      end do
      call output ( "Total number of profiles" , advance='yes')
      call output ( sum1 , advance='no' )
      call output ( sum2 , advance='no' )
      call output ( sum3 , advance='yes' )
    endif
    ! Now accumulate hGridOffsets which currently contains the number
    ! of non-overlap profiles each chunk/hGrid.  After this it will contain
    ! the accumulated number.  This is equivalent to storing the
    ! index of the last profile in each chunk.
    do chunk = 2, size ( chunks )
      chunks(chunk)%hGridOffsets = chunks(chunk)%hGridOffsets + &
        & chunks(chunk-1)%hGridOffsets
    end do
    ! Now fill (i.e. spread) the hGridTotals array
    do chunk = 1, size ( chunks )
      chunks(chunk)%hGridTotals = chunks(size(chunks))%hGridOffsets
    end do
    ! Now move hGridOffsets back one to get the offsets we want.
    ! Each hGridOffsets will now contain the last profile index in the preceeding
    ! chunk, which is the offset we desire.
    do chunk = size ( chunks ), 2, -1
      chunks(chunk)%hGridOffsets = chunks(chunk-1)%hGridOffsets
    end do
    if ( switchDetail(switches, 'pro') >= 0 .or. DEEBUG ) then
      call output ( 'chunks(1)%hGridOffsets: ', advance='no' )
      call output ( chunks(1)%hGridOffsets, advance='yes' )
    endif
    if ( ChunkDivideConfig%allowPriorOverlaps .and. &
      & chunks(1)%noMAFsLowerOverlap > 0 ) then
      chunks(1)%hGridOffsets = 0
      do chunk=1, size(chunks)
        chunks(chunk)%hGridOffsets = chunks(chunk)%hGridOffsets + LowerOverlaps
      enddo
    else
      chunks(1)%hGridOffsets = 0
    endif
    
    if ( switchDetail(switches, 'pro') >= 0 .or. DEEBUG ) then
      call output ( "Dumping offsets, hgridTotals for all chunks: " , &
        & advance='yes')
      do chunk = 1, size ( chunks )
        call output ( chunk )
        call blanks ( 3 )
        call output ( chunks(chunk)%hGridOffsets )
        call blanks ( 3 )
        call output ( chunks(chunk)%hGridTotals, advance='yes' )
      end do
    endif

    call deAllocate_Test ( LowerOverlaps, 'LowerOverlaps', ModuleName )
    if ( specialDumpFile /= ' ' ) call revertOutput
    call trace_end ( "ComputeAllHGridOffsets", cond=toggle(gen) )
  end subroutine ComputeAllHGridOffsets

  ! ------------------------------  ComputeNextChunksHGridOffsets  -----
  ! This routine is CURRENTLY NOT USED!
  subroutine ComputeNextChunksHGridOffsets ( chunks, chunkNo, hGrids )
    use Allocate_Deallocate, only :ALLOCATE_TEST
    use Chunks_m, only: MLSChunk_T
    use HGridsDatabase, only: HGRID_T
    ! Dummy arguments
    type(MLSChunk_T), dimension(:), intent(inout) :: CHUNKS
    integer, intent(in) :: CHUNKNO
    type(HGrid_T), dimension(:), intent(in) :: HGRIDS
    ! Executable code
    ! Do nothing if this is the last chunk
    if ( chunkNo == size(chunks) ) return
    ! Allocate the array
    call Allocate_Test ( chunks(chunkNo+1)%hGridOffsets, size(hGrids), &
      & 'chunks(?)%hGridOffsets', ModuleName )
    ! Compute our number of output instances
    chunks(chunkNo+1)%hGridOffsets = hGrids%noProfs - &
      & hGrids%noProfsLowerOverlap - &
      & hGrids%noProfsUpperOverlap
    ! For later chunks, add on the accumulated previous stuff
    if ( chunkNo > 1 ) &
      & chunks(chunkNo+1)%hGridOffsets = chunks(chunkNo+1)%hGridOffsets + &
      & chunks(chunkNo)%hGridOffsets + 1
  end subroutine ComputeNextChunksHGridOffsets
  
  ! ------------------------------------- PlaceHGridContents --
  subroutine PlaceHGridContents ( HGrid1, HGrid2, offset )
    ! Place the contents of one Hgrid1 inside HGrid2, possibly offset
    use HGRIDSDATABASE, only: HGRID_T
    ! Args
    type(HGRID_T), intent(in)     :: HGrid1
    type(HGRID_T), intent(inout)  :: HGrid2
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

    use LEXER_CORE, only: PRINT_SOURCE
    use OUTPUT_M, only: OUTPUT
    use TREE, only: WHERE_AT => WHERE

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
    case ( noSpacingOrigin )
      call output ( "TYPE = Regular but no spacing and/or origin is specified", &
        & advance='yes' )
    case ( badTime )
      call output ( "Bad information given for date in explicit hGrid", &
        & advance='yes' )
    end select
    end subroutine ANNOUNCE_ERROR
    
    subroutine PlaceArray_r4(array1, array2, offset)
      ! place contents of array1 inside array2, possibly offset
      use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
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
      use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR
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

