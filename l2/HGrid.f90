! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module HGrid                    ! Horizontal grid information
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use EmpiricalGeometry, only: EmpiricalLongitude, ChooseOptimumLon0
  use EXPR_M, only: EXPR
  use Dump_0, only: DUMP
  use INIT_TABLES_MODULE, only: F_FRACTION, F_HEIGHT, F_INCLINATION, &
    & F_INTERPOLATIONFACTOR, F_MIF, F_MODULE, F_TYPE, F_VALUES, FIELD_FIRST,&
    & F_SPACING, F_ORIGIN, F_FORBIDOVERSPILL, &
    & F_SOURCEL2GP, FIELD_LAST, L_EXPLICIT, L_FIXED, L_FRACTIONAL, L_HEIGHT, L_L2GP,&
    & L_MIF, L_REGULAR, PHYQ_DIMENSIONLESS, PHYQ_LENGTH, PHYQ_ANGLE
  use LEXER_CORE, only: PRINT_SOURCE
  use L1BData, only: DeallocateL1BData, L1BData_T, ReadL1BData
  use L2GPData, only: L2GPDATA_T
  use MLSCommon, only: L1BInfo_T, MLSChunk_T, NameLen, R8, TAI93_RANGE_T
  use MLSMessageModule, only: MLSMessage, MLSMSG_allocate, &
    & MLSMSG_DeAllocate, MLSMSG_Error, MLSMSG_Info, MLSMSG_L1BRead
  use MLSNumerics, only: HUNT, InterpolateValues
  use MoreTree, only: GET_BOOLEAN
  use OUTPUT_M, only: OUTPUT
  use STRING_TABLE, only: GET_STRING
  use TRACE_M, only: TRACE_BEGIN, TRACE_END
  use TOGGLES, only: GEN, TOGGLE, SWITCHES
  use TREE, only: DECORATION, DUMP_TREE_NODE, NSONS, NULL_TREE, SOURCE_REF, &
                  SUB_ROSA, SUBTREE
  use UNITS, only: DEG2RAD, RAD2DEG

  implicit none

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  type HGrid_T
    integer :: NAME                 ! String index of name.
    integer :: noProfs              ! Number of profiles in this grid
    integer :: noProfsLowerOverlap  ! Number of profiles in the lower overlap
    integer :: noProfsUpperOverlap  ! Number of profiles in the upper overlap

    ! Now the various coordinates in the HGrid, all dimensioned (noProfs)
    real(r8), dimension(:), pointer :: phi => NULL()
    real(r8), dimension(:), pointer :: geodLat => NULL()
    real(r8), dimension(:), pointer :: lon => NULL()
    real(r8), dimension(:), pointer :: time => NULL()
    real(r8), dimension(:), pointer :: solarTime => NULL()
    real(r8), dimension(:), pointer :: solarZenith => NULL()
    real(r8), dimension(:), pointer :: losAngle => NULL()
  end type HGrid_T

! -----     Private declarations     ---------------------------------

  integer, private :: ERROR

! Error codes for "announce_error"
  integer, private, parameter :: AngleUnitMessage = 1
  integer, private, parameter :: LengthUnitMessage = AngleUnitMessage + 1
  integer, private, parameter :: NoFraction = LengthUnitMessage + 1
  integer, private, parameter :: NoHeight = NoFraction + 1
  integer, private, parameter :: UnitlessMessage = NoHeight + 1
  integer, private, parameter :: NoModule = UnitlessMessage + 1
  integer, private, parameter :: NoMIF = NoModule + 1
  integer, private, parameter :: NoSpacingOrigin = NoMIF + 1

contains ! =====     Public Procedures     =============================

  ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( database, item )

    ! Dummy arguments
    type (HGrid_T), dimension(:), pointer :: database
    type (HGrid_T), intent(in) :: item

    ! Local variables
    type (HGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddHGridToDatabase = newSize
  end function AddHGridToDatabase

  ! -----------------------------------  CreateHGridFromMLSCFInfo  -----
  type(hGrid_T) function CreateHGridFromMLSCFInfo &
    & ( name, root, l1bInfo, l2gpDatabase, processingRange, chunks, chunkNo ) result ( hGrid )

  ! This routine creates an hGrid based on the user requests.

    ! Dummy arguments
    integer, intent(in) :: NAME               ! String index of name
    integer, intent(in) :: ROOT               ! Root of hGrid subtree
    type (L1BInfo_T), intent(in) :: L1BINFO   ! File handles for l1b data
    type (L2GPData_T), pointer, dimension(:) :: L2GPDATABASE
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (MLSChunk_T), intent(in), dimension(:) :: CHUNKS ! The chunks
    integer, intent(in) :: CHUNKNO

    ! Local variables
    integer :: EXPR_UNITS(2)            ! Output from Expr subroutine
    double precision :: EXPR_VALUE(2)   ! Output from Expr subroutine
    integer :: hGridType
    real(r8) :: interpolationFactor
    integer :: instrumentModule
    real(r8) :: fraction, height
    real(r8) :: spacing, origin
    type (L2GPData_T), pointer :: L2GP  ! The l2gp to use 

    integer :: keyNo                    ! Entry in the mlscf line
    integer :: fieldValue               ! Node in the tree

    type (L1BData_T) :: l1bField ! L1B data

    real(r8) :: incline                 ! Orbital inclination / degrees
    real(r8) :: MINTIME, MAXTIME        ! Span for a chunk
    integer, dimension(2) :: PROFRANGE  ! Profile range
    integer :: A,B                      ! Elements of profile range

    integer :: FIELD                    ! Subtree index of "field" node
    integer :: FIELD_INDEX              ! F_..., see Init_Tables_Module
    logical :: GOT_FIELD(field_first:field_last)
    integer :: L1BFLAG
    integer :: L1BITEM                  ! Loop counter
    integer :: MAF                      ! Loop counters
    integer :: MIF                      ! For fixed hGrids
    integer :: NOMAFS                   ! Number of MAFs of L1B data read
    integer :: PROF                     ! Loop counter
    integer :: SON                      ! Son of Root
    integer :: STATUS                   ! From Allocate, ReadL1B... etc.
    integer :: VALUESNODE               ! Node of tree for explicit
    logical :: FORBIDOVERSPILL          ! If set don't allow overlaps beyond L1B

    character (len=NameLen) :: InstrumentModuleName

    ! Executable code

    hGrid%name = name
    if ( toggle(gen) ) call trace_begin ( "CreateHGridFromMLSCFInfo", root )

    got_field = .false.
    interpolationFactor = 1.0
    forbidOverspill = .false.

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
      case ( f_module )
        instrumentModule = sub_rosa(subtree(2,son))
        call get_string ( instrumentModule , instrumentModuleName )
      case  ( f_forbidOverspill )
        forbidOverspill = get_boolean ( fieldValue )
      case ( f_height )
        call expr ( subtree(2,son), expr_units, expr_value )
        height = expr_value(1)
        if ( expr_units(1) /= PHYQ_Length) &
          & call announce_error ( field, lengthUnitMessage )
      case ( f_mif )
        call expr ( subtree(2,son), expr_units, expr_value )
        mif = expr_value(1)
        if ( expr_units(1) /= PHYQ_Dimensionless) &
          & call announce_error ( field, lengthUnitMessage )
      case ( f_fraction )
        call expr ( subtree(2,son), expr_units, expr_value )
        fraction = expr_value(1)
        if ( expr_units(1) /= PHYQ_Dimensionless) &
          & call announce_error ( field, unitlessMessage )
      case ( f_interpolationFactor )
        call expr ( subtree(2,son), expr_units, expr_value )
        interpolationFactor = expr_value(1)
        if ( expr_units(1) /= PHYQ_Dimensionless) &
          & call announce_error ( field, unitlessMessage )
      case ( f_spacing )
        call expr ( subtree(2,son), expr_units, expr_value )
        spacing = expr_value(1)
        if ( expr_units(1) /= PHYQ_Angle) &
          & call announce_error ( field, angleUnitMessage )
      case ( f_origin )
        call expr ( subtree(2,son), expr_units, expr_value )
        origin = expr_value(1)
        if ( expr_units(1) /= PHYQ_Angle) &
          & call announce_error ( field, angleUnitMessage )
      case ( f_values )
        valuesNode = son
      case ( f_sourceL2gp )
        l2gp => l2gpDatabase(decoration(decoration(subtree(2,son))))
      case ( f_inclination )
        call expr (subtree ( 2, son), expr_units, expr_value )
        incline = expr_value(1)
      case default ! Can't get here if tree_checker works correctly
      end select
    end do

    select case (hGridType)

    case ( l_height, l_fractional, l_fixed ) ! ----- Fractional or Height ------
      if (.not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      else
        call CreateMIFBasedHGrids ( l1bInfo, hGridType, chunks(chunkNo), &
          & got_field, root, height, fraction, interpolationFactor, &
          & instrumentModuleName, mif, hGrid )
      end if

    case ( l_explicit ) ! ----------------- Explicit ------------------
      ! For explicit hGrids, do things differently
      hGrid%noProfs = nsons(valuesNode)-1
      hGrid%noProfsLowerOverlap = 0
      hGrid%noProfsUpperOverlap = 0
      call CreateEmptyHGrid(hGrid)
      hGrid%lon = 0.0_r8
      hGrid%time = 0.0_r8
      hGrid%solarTime = 0.0_r8
      hGrid%losAngle = 0.0_r8
      hGrid%solarZenith = 0.0_r8
      
      do prof = 1, hGrid%noProfs
        call expr ( subtree ( prof+1, valuesNode), expr_units, expr_value )
        hGrid%phi(prof) = expr_value(1)
        hGrid%geodLat(prof) = hGrid%phi(prof) !???? Sort this out later!
      end do

    case ( l_regular ) ! ----------------------- Regular --------------
      if (.not. got_field(f_module) ) then
        call announce_error ( root, NoModule )
      else if ( .not. all(got_field((/f_spacing, f_origin/)))) then
        call announce_error ( root, NoSpacingOrigin )
      else
        call CreateRegularHGrid ( l1bInfo, processingRange, chunks, chunkNo, &
          & spacing, origin, trim(instrumentModuleName), forbidOverspill, hGrid )
      end if

    case ( l_l2gp) ! -------------------- L2GP ------------------------
      
      ! Get the time from the l1b file
      call ReadL1BData ( l1bInfo%l1boaID, "MAFStartTimeTAI", l1bField, noMAFs, &
        & l1bFlag, &
        & firstMAF=chunks(chunkNo)%firstMAFIndex, &
        & lastMAF=chunks(chunkNo)%lastMAFIndex)
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
      hGrid%phi =         l2gp%geodAngle(a:b)
      hGrid%geodLat =     l2gp%latitude(a:b)
      hGrid%lon =         l2gp%longitude(a:b)
      hGrid%time =        l2gp%time(a:b)
      hGrid%solarTime =   l2gp%solarTime(a:b)
      hGrid%solarZenith = l2gp%solarZenith(a:b)
      hGrid%losAngle =    l2gp%losAngle(a:b)

    end select
    
    if ( toggle(gen) ) call trace_end ( "CreateHGridFromMLSCFInfo" )

    if ( index ( switches, 'geom' ) /= 0 ) &
      & call DumpChunkHGridGeometry ( hGrid, chunks(chunkNo), &
      & trim(instrumentModuleName), l1bInfo )

  end function CreateHGridFromMLSCFInfo

  ! ------------------------------------- CreateMIFBasedHGrids -----------
  subroutine CreateMIFBasedHGrids ( l1bInfo, hGridType, &
    & chunk, got_field, root, height, fraction, interpolationFactor,&
    & instrumentModuleName, mif, hGrid )
    ! This is part of ConstructHGridFromMLSCFInfo
    type (L1BInfo_T), intent(in)      :: L1BINFO
    integer, intent(in)               :: HGRIDTYPE
    type (MLSChunk_T), intent(in)     :: CHUNK
    logical, dimension(:), intent(in) :: GOT_FIELD
    integer, intent(in)               :: ROOT
    real(r8), intent(in)              :: INTERPOLATIONFACTOR
    real(r8), intent(in)              :: HEIGHT
    real(r8), intent(in)              :: FRACTION
    character (len=*), intent(in)     :: INSTRUMENTMODULENAME
    integer, intent(in)               :: MIF
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

    real(r8), parameter :: SIXTH = 1.0_r8 / 6.0_r8

    character (len=15), dimension(noL1BItemsToRead) :: L1bItemNames
    ! Entries in the above array following FirstModularItem are prefixed
    ! with either GHz or THz. 

    integer :: L1BFLAG                  ! Flag
    integer :: L1BITEM                  ! Loop counter
    integer :: MAF                      ! Loop counter etc.
    real(r8) :: MinAngle, MaxAngle
    integer :: NOMAFS                   ! Dimension
    integer :: STATUS                   ! Flag
    real(r8), dimension(:,:,:), pointer :: TpGeodAlt, TpGeodAngle

    ! MIFs it would choose in the non over/undersampled case
    real(r8), dimension(:), pointer :: defaultField, interpolatedField

    type (L1BData_T) :: l1bField ! L1B data
    integer, dimension(:), pointer :: defaultMIFs
    real(r8), dimension(:), pointer :: defaultIndex
    real(r8), dimension(:), pointer :: interpolatedIndex
    character (len=NameLen) :: L1BItemName

    ! Executable code

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

    end select

    if ( hGridType /= l_Fixed ) then
      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField, noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex )
      if ( l1bFlag==-1) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
    else
      noMAFs = chunk%lastMafIndex - chunk%firstMafIndex + 1
    endif

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
    
    ! Create an empty hGrid
    hGrid%noProfs = NINT(noMAFs*interpolationFactor)
    call CreateEmptyHGrid(hGrid)
    
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
      call ReadL1BData ( l1bInfo%l1boaid, l1bItemName, l1bField,noMAFs, &
        & l1bFlag, firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex )
      if ( l1bFlag==-1 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_L1BRead//l1bItemName )
      
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
		defaultIndex(maf)=maf*1.0_r8 
	end do
	do maf=1,hGrid%noProfs 
		interpolatedIndex(maf)=maf/interpolationFactor
	end do

      	select case ( l1bItem )
        case ( l1b_tpLon )
	  call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
		interpolatedField,method='Linear',rangeofPeriod=(/-180.0_r8,180.0_r8/))
        case ( l1b_tpSolarTime )
	  call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
		interpolatedField,method='Linear',rangeofPeriod=(/0.0_r8,24.0_r8/))
        case ( l1b_tpLosAngle )
	  call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
		interpolatedField,method='Linear',rangeofPeriod=(/0.0_r8,360.0_r8/))
	case default
	  call InterpolateValues(defaultIndex,defaultField,interpolatedIndex,&
		interpolatedField,method='Linear')
        end select

      end if
      
      select case ( l1bItem )
      case ( l1b_MAFStartTimeTAI )
        hGrid%time = interpolatedField
      case ( l1b_tpGeodLat )
        hGrid%geodLat = interpolatedField
      case ( l1b_tpLon )
        hGrid%lon = interpolatedField
      case ( l1b_tpGeodAngle )
        hGrid%phi = interpolatedField
      case ( l1b_tpSolarZenith )
        hGrid%solarZenith = interpolatedField
      case ( l1b_tpSolarTime )
        hGrid%solarTime = interpolatedField
      case ( l1b_tpLosAngle )
        hGrid%losAngle = interpolatedField
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
    
  end subroutine CreateMIFBasedHGrids

  ! ----------------------------------- CreateRegularHGrid ------------
  subroutine CreateRegularHGrid ( l1bInfo, processingRange, chunks, chunkNo, &
    & spacing, origin, instrumentModuleName, forbidOverspill, hGrid )
    type (L1BInfo_T), intent(in) :: L1BINFO
    type (TAI93_Range_T), intent(in) :: PROCESSINGRANGE
    type (MLSChunk_T), intent(in), dimension(:) :: CHUNKS
    integer, intent(in) :: CHUNKNO
    real(r8), intent(in) :: SPACING
    real(r8), intent(in) :: ORIGIN
    character (len=*), intent(in) :: INSTRUMENTMODULENAME
    logical, intent(in) :: FORBIDOVERSPILL
    type (HGrid_T), intent(inout) :: HGRID ! Needs inout as name set by caller

    ! Local variables/parameters
    real(r8), parameter :: SECONDSINDAY = 24*60*60
    ! Note this next one is ONLY NEEDED for the case where we have only
    ! one MAF in the chunk
    real(r8), parameter :: ORBITALPERIOD = 98.8418*60.0

    integer :: NOMAFS                   ! From ReadL1B
    integer :: FLAG                     ! From ReadL1B
    integer :: I                        ! Loop counter
    integer :: FIRSTPROFINRUN           ! Index of first profile in processing time
    integer :: LASTPROFINRUN            ! Index of last profile in processing time

    real(r8) :: MINANGLE                ! Smallest angle in chunk
    real(r8) :: MAXANGLE                ! Largest angle in chunk
    real(r8) :: FIRST                   ! First point in of hGrid
    real(r8) :: LAST                    ! Last point in hGrid
    real(r8), dimension(:), pointer :: MIF1GEODANGLE ! For first mif
    real(r8) :: DAYSTART                ! Start of day
    real(r8) :: INCLINE                 ! Mean orbital inclination
    real(r8) :: DELTA                   ! A change in angle
    real(r8) :: NEXTANGLE               ! First non ovl. MAF for next chunk

    type (L1BData_T) :: L1BFIELD        ! A field read from L1 file
    type (MLSChunk_T) :: CHUNK

    ! Executable code
    chunk = chunks ( chunkNo )

    ! Setup the empircal geometry estimate of lon0
    ! (it makes sure it's not done twice
    call ChooseOptimumLon0 ( l1bInfo, chunk )

    ! First we're going to work out the geodetic angle range
    ! Read an extra MAF if possible, as we may use it later
    ! when computing the overlaps.
    call ReadL1BData ( l1bInfo%l1bOAID, instrumentModuleName//".tpGeodAngle", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex+1 )
    noMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
    minAngle = minval ( l1bField%dpField(1,:,1) )
    maxAngle = maxval ( l1bField%dpField(1,:,noMAFs) )
    nullify ( mif1GeodAngle )
    call Allocate_test ( mif1GeodAngle, noMAFs, 'mif1Geodangle', ModuleName )
    mif1GeodAngle = l1bField%dpField(1,1,1:noMAFs)

    ! Get or guess the start of the next chunk.
    i = noMAFs - chunk%noMAFsUpperOverlap + 1
    if ( i < noMAFs + 1 ) then
      nextAngle = l1bField%dpField(1,1,i)
    else
      nextAngle = maxAngle + spacing
    endif

    call DeallocateL1BData ( l1bField )

    ! Now choose the geodetic angles for the hGrid
    ! First identify the first point - the one closest to the start of the first MAF
    first = origin + spacing * int ( (minAngle-origin)/spacing )
    delta = first - minAngle            ! So +ve means first could be smaller
    if ( delta > spacing/2 ) then
      first = first - spacing
    else if ( delta < -spacing/2 ) then
      first = first + spacing
    end if

    ! Now work out the last point in a similar manner
    last = origin + spacing * int ( (maxAngle-origin)/spacing )
    delta = last - maxAngle            ! So +ve means last could be smaller
    if ( delta > spacing/2 ) then
      last = last - spacing
    else if ( delta < -spacing/2 ) then
      last = last + spacing
    end if

    ! Now outset by one to be sure
    first = first - spacing
    last = last + spacing

    ! Now work out how many profiles that is and lay them down
    hGrid%noProfs = nint( (last-first) / spacing ) + 1
    call CreateEmptyHGrid ( hGrid )
    do i = 1, hGrid%noProfs
      hGrid%phi(i) = first + (i-1)*spacing
    end do

    if ( index ( switches, 'hgrid' ) /= 0 ) then
      call output ( 'Constructing regular hGrid', advance='yes' )
      call output ( 'minAngle: ' )
      call output ( minAngle, format='(F7.2)' )
      call output ( ' maxAngle: ' )
      call output ( maxAngle, format='(F7.2)' )
      call output ( ' nextAngle: ' )
      call output ( nextAngle, format='(F7.2)', advance='yes' )
      call output ( 'Spacing: ' )
      call output ( spacing )
      call output ( ' first: ' )
      call output ( first, format='(F7.2)' )
      call output ( ' last: ' )
      call output ( last, format='(F7.2)', advance='yes' )
    end if

    ! Now fill the other geolocation information, first latitude
    ! Get orbital inclination
    call ReadL1BData ( l1bInfo%l1bOAID, "scOrbIncl", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex )
    ! Use the average of all the first MIFs to get inclination for chunk
    incline = sum ( l1bField%dpField(1,1,:) ) / noMAFs
    call DeallocateL1BData ( l1bField )
    hGrid%geodLat = rad2deg * asin ( sin( deg2Rad*hGrid%phi ) * &
      & sin ( deg2Rad*incline ) )

    ! Now longitude
    call EmpiricalLongitude ( hGrid%phi, hGrid%lon )

    ! Now time, because this is important to get right, I'm going to put in
    ! special code for the case where the chunk is of length one.
    call ReadL1BData ( l1bInfo%l1bOAID, "MAFStartTimeTAI", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    if ( chunk%firstMAFIndex /= chunk%lastMAFIndex ) then
      call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
        & hGrid%phi, hGrid%time, &
        & method='Spline', extrapolate='Allow' )
    else
      ! Case where only single MAF per chunk, treat it specially
      hGrid%time = l1bField%dpField(1,1,1) + &
        & ( OrbitalPeriod/360.0 ) * ( hGrid%phi - mif1GeodAngle(1) )
    end if
    call DeallocateL1BData ( l1bField )
      
    ! Solar time
    ! First get fractional day, note this neglects leap seconds.
    ! Perhaps fix this later !???????? NJL. We do have access to the
    ! UTC ascii time field, perhaps we could use that?
    hGrid%solarTime = modulo ( hGrid%time, secondsInDay ) / secondsInDay
    ! Now correct for longitude and convert to hours
    hGrid%solarTime = 24.0 * ( hGrid%solarTime + hGrid%lon/360.0 )
    hGrid%solarTime = modulo ( hGrid%solarTime, 24.0_r8 )

    ! Solar zenith
    ! This we'll have to do with straight interpolation
    call ReadL1BData ( l1bInfo%l1bOAID, instrumentModuleName//".tpSolarZenith", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
      & hGrid%phi, hGrid%solarZenith, &
      & method='Spline', extrapolate='Allow' )
    call DeallocateL1BData ( l1bField )

    ! Line of sight angle
    ! This we'll have to do with straight interpolation
    call ReadL1BData ( l1bInfo%l1bOAID, instrumentModuleName//".tpLosAngle", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, lastMAF=chunk%lastMAFIndex )
    call InterpolateValues ( mif1GeodAngle, l1bField%dpField(1,1,:), &
      & hGrid%phi, hGrid%losAngle, &
      & method='Spline', extrapolate='Allow' )
    call DeallocateL1BData ( l1bField )
    hGrid%losAngle = modulo ( hGrid%losAngle, 360.0_r8 )

    ! Now work out how much of this HGrid is overlap
    ! The deal will be the first legitimate profile is the first one who's phi
    ! is above the first non overlapped MAF.
    ! The exceptions are for the first and last chunks, where the later test
    ! for the processing time range is the limiting factor
    if ( chunkNo > 1 ) then
      call Hunt ( hGrid%phi, mif1GeodAngle(chunk%noMAFsLowerOverlap+1), &
        & hGrid%noProfsLowerOverlap, allowTopValue=.true., allowBelowValue=.true. )
      ! So the hunt returns the index of the last overlapped, which is
      ! the number we want to be in the overlap.
    else
      hGrid%noProfsLowerOverlap = 0
    end if

    if ( chunkNo < size(chunks) ) then
      call Hunt ( hGrid%phi, nextAngle, &
        & hGrid%noProfsUpperOverlap, allowTopValue=.true., allowBelowValue=.true. )
      ! Here the hunt returns the index of the last non overlapped profile
      ! So we do a subtraction to get the number in the overlap.
      hGrid%noProfsUpperOverlap = hGrid%noProfs - hGrid%noProfsUpperOverlap
    else
      hGrid%noProfsUpperOverlap = 0
    end if

    if ( index ( switches, 'hgrid' ) /= 0 ) then
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
      call Hunt ( hGrid%time, processingRange%startTime, &
        & firstProfInRun, allowTopValue=.true., allowBelowValue=.true. )
      if ( firstProfInRun > 0 .and. forbidOverspill ) then
        if ( index ( switches, 'hgrid' ) /= 0 ) &
          & call output ( 'hGrid starts before run trimming start.', &
          & advance='yes' )
        call TrimHGrid ( hGrid, -1, firstProfInRun )
      end if
      
      call Hunt ( hGrid%time, processingRange%endTime, &
        & lastProfInRun, allowTopValue=.true., allowBelowValue=.true. )
      if ( lastProfInRun < hGrid%noProfs .and. forbidOverspill ) then
        if ( index ( switches, 'hgrid' ) /= 0 ) &
          & call output ( 'hGrid ends after run trimming end.',&
          & advance='yes' )
        call TrimHGrid ( hGrid, 1, hGrid%noProfs-lastProfInRun )
      end if
    end if

    ! Now a 'softer' limit that applies to all cases, this just moves the
    ! overlap regions around if necessary to deal with overspill.
    if ( hGrid%time(hGrid%noProfsLowerOverlap+1) < processingRange%startTime ) then
      if ( index ( switches, 'hgrid' ) /= 0 ) &
        & call output ( &
        & 'Non overlapped part of hGrid starts before run, extending overlap.', &
        & advance='yes' )
      call Hunt ( hGrid%time, processingRange%startTime, &
        &   hGrid%noProfsLowerOverlap, allowTopValue=.true., allowBelowValue=.true. )
    end if

    if ( hGrid%time(hGrid%noProfs-hGrid%noProfsUpperOverlap) > &
      & processingRange%endTime ) then
      if ( index ( switches, 'hgrid' ) /= 0 ) &
        & call output ( &
        & 'Non overlapped part of hGrid end after run, extending overlap.', &
        & advance='yes' )
      call Hunt ( hGrid%time, processingRange%endTime, &
        & hGrid%noProfsUpperOverlap, allowTopValue=.true., allowBelowValue=.true. )
      hGrid%noProfsUpperOverlap = hGrid%noProfs - hGrid%noProfsUpperOverlap
    end if

    ! Finally we're done.
    if ( index ( switches, 'hgrid' ) /= 0 ) then
      call output ( 'Final Hgrid size: ' )
      call output ( hGrid%noProfs )
      call output ( ', overlaps: ' )
      call output ( hGrid%noProfsLowerOverlap )
      call output ( ', ' )
      call output ( hGrid%noProfsUpperOverlap, advance='yes' )
    end if

    ! That's it
    call Deallocate_test ( mif1GeodAngle, 'mif1GeodAngle', ModuleName )

  end subroutine CreateRegularHGrid

  ! ----------------------------------- CreateEmptyHGrid ---------------
  subroutine CreateEmptyHGrid ( hGrid )
    ! Just does allocates etc.
    type (HGrid_T), intent(inout) :: HGRID

    ! Executable code
    call Allocate_Test ( hGrid%phi, hGrid%noProfs, 'hGrid%phi', ModuleName)
    call Allocate_Test ( hGrid%geodLat, hGrid%noProfs, 'hGrid%geodLat', ModuleName)
    call Allocate_Test ( hGrid%lon, hGrid%noProfs, 'hGrid%lon', ModuleName)
    call Allocate_Test ( hGrid%time, hGrid%noProfs, 'hGrid%time', ModuleName)
    call Allocate_Test ( hGrid%solarTime, hGrid%noProfs, 'hGrid%solarTime', ModuleName)
    call Allocate_Test ( hGrid%solarZenith, hGrid%noProfs, 'hGrid%solarZenith', ModuleName)
    call Allocate_Test ( hGrid%losAngle, hGrid%noProfs, 'hGrid%losAngle', ModuleName)

  end subroutine CreateEmptyHGrid

  ! ---------------------------------------  TrimHGrid ---
  subroutine TrimHGrid ( hGrid, side, NOTODELETE )
    type (HGrid_T), intent(inout) :: HGrid
    integer, intent(in) :: SIDE         ! -1 = lower, 1 = upper
    integer, intent(in) :: NOTODELETE ! How many to delete, default all
    ! Local variables
    integer :: newNoProfs
    real(r8), dimension(:), pointer :: temp
    integer :: first, last

    ! Executable code
    select case ( side )
    case ( -1 )
      newNoProfs = hGrid%noProfs - noToDelete
      first = noToDelete + 1
      last = hGrid%noProfs
      hGrid%noProfsLowerOverlap = max ( hGrid%noProfsLowerOverlap - noToDelete, 0 )
    case ( 1 )
      newNoProfs = hGrid%noProfs - noToDelete
      first = 1
      last = newNoProfs
      hGrid%noProfsUpperOverlap = max ( hGrid%noProfsUpperOverlap - noToDelete, 0 )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Invalid side argument to TrimHGrid' )
    end select
    if ( newNoProfs <= 0 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Too many profiles to delete' )

    nullify ( temp )
    call allocate_test ( temp, hGrid%noProfs, 'temp', ModuleName )
    hGrid%noProfs = newNoProfs

    ! Now allocate each entry and trim it
    temp = hGrid%phi                    ! ------------------------- Phi
    call Allocate_Test ( hGrid%phi, newNoProfs, 'hGrid%phi', ModuleName)
    hGrid%phi = temp ( first : last )
    temp = hGrid%geodLat                ! ------------------------- GeodLat
    call Allocate_Test ( hGrid%geodLat, newNoProfs, 'hGrid%geodLat', ModuleName)
    hGrid%geodLat = temp ( first : last )
    temp = hGrid%lon                    ! ------------------------- Lon
    call Allocate_Test ( hGrid%lon, newNoProfs, 'hGrid%lon', ModuleName)
    hGrid%lon = temp ( first : last )
    temp = hGrid%time                   ! ------------------------- Time
    call Allocate_Test ( hGrid%time, newNoProfs, 'hGrid%time', ModuleName)
    hGrid%time = temp ( first : last )
    temp = hGrid%solarTime              ! ------------------------- SolarTime
    call Allocate_Test ( hGrid%solarTime, newNoProfs, 'hGrid%solarTime', ModuleName)
    hGrid%solarTime = temp ( first : last )
    temp = hGrid%solarZenith            ! ------------------------- SolarZenith
    call Allocate_Test ( hGrid%solarZenith, newNoProfs, 'hGrid%solarZenith', ModuleName)
    hGrid%solarZenith = temp ( first : last )
    temp = hGrid%losAngle               ! ------------------------- LosAngle
    call Allocate_Test ( hGrid%losAngle, newNoProfs, 'hGrid%losAngle', ModuleName)
    hGrid%losAngle = temp ( first : last )

    call Deallocate_test ( temp, 'temp', ModuleName )
  end subroutine TrimHGrid
    
  ! ---------------------------------------  DestroyHGridContents  -----
  subroutine DestroyHGridContents ( hGrid )

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(inout) :: hGrid

    ! Local Variables
    integer :: STATUS    ! From Deallocate

    ! Executable code
    
    call deallocate_test ( hGrid%phi, 'hGrid%phi', ModuleName )
    call deallocate_test ( hGrid%geodLat, 'hGrid%geodLat', ModuleName )
    call deallocate_test ( hGrid%lon, 'hGrid%lon', ModuleName )
    call deallocate_test ( hGrid%time, 'hGrid%time', ModuleName )
    call deallocate_test ( hGrid%solarTime, 'hGrid%solarTime', ModuleName )
    call deallocate_test ( hGrid%solarZenith, 'hGrid%solarZenith', ModuleName )
    call deallocate_test ( hGrid%losAngle, 'hGrid%losAngle', ModuleName )
    hGrid%noProfs = 0
  end subroutine DestroyHGridContents

  ! ---------------------------------------  DestroyHGridDatabase  -----
  subroutine DestroyHGridDatabase ( database )

  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (HGrid_T), dimension(:), pointer :: database

    ! Local variables
    integer :: hGridIndex, status
    if ( associated(database) ) then
      do hGridIndex=1,SIZE(database)
        call DestroyHGridContents ( database(hGridIndex) )
      end do
      deallocate ( database, stat=status )
      if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // "database" )
    end if
  end subroutine DestroyHGridDatabase

  ! ----------------------------------------- DumpChunkHGridGeometry
  subroutine DumpChunkHGridGeometry ( hGrid, chunk, &
    & instrumentModuleName, l1bInfo )
    type (HGrid_T), intent(in) :: HGRID
    type (L1BInfo_T), intent(in) :: L1BINFO
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

    ! Executable code

    ! Read the geodetic angle from the L1Bfile
    call ReadL1BData ( l1bInfo%l1bOAID, instrumentModuleName//".tpGeodAngle", &
      & l1bField, noMAFs, flag, &
      & firstMAF=chunk%firstMAFIndex, &
      & lastMAF=chunk%lastMAFIndex )
    mifPhi => l1bField%dpField(1,:,:)

    phiMin = minval ( mifPhi )
    phiMax = maxval ( mifPhi )
    phiMin = min ( phiMin, hGrid%phi(1) )
    phiMax = max ( phiMax, hGrid%phi(hGrid%noProfs) )

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
      charMin = ( hGrid%phi(prof) - phiMin ) / binSize + 1
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
    if ( mod ( noBins, width ) /= 0 ) noWindows = noWindows + 1
    do window = 1, noWindows
      call output ( ' ', advance='yes' )
      windowSize = width
      start = ( window-1 ) * width + 1
      if ( window == noWindows ) windowSize = mod ( noBins, width )
      do line = 1, noLines
        do char = start, start+windowSize-2
          call output ( text ( char, line ) )
        end do
        call output  ( text ( start+windowSize-1, line ), advance='yes' ) 
      end do
      do char = start, start+windowSize-2
        call output ( '-' )
      end do
      call output ( '-', advance='yes' )
    end do

    ! Finish off
    deallocate ( text, stat=flag )
    if ( flag /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'text' )
    call DeallocateL1BData ( l1bField )

  end subroutine DumpChunkHGridGeometry
 
! =====     Private Procedures     =====================================

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( WHERE, CODE )
    integer, intent(in) :: WHERE   ! Tree node where error was noticed
    integer, intent(in) :: CODE    ! Code for error message

    error = max(error,1)
    call output ( '***** At ' )
    call print_source ( source_ref(where) )
    call output ( ': ' )
    select case ( code )
    case ( angleUnitMessage )
      call output ( "Value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to be an angle, e.g. degrees", &
        advance='yes' )
    case ( lengthUnitMessage )
      call output ( "Value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to be a length, e.g. km", &
        advance='yes' )
    case ( noFraction )
      call output ( "TYPE = FRACTIONAL but no fraction is specified", &
        & advance='yes' )
    case ( noHeight )
      call output ( "TYPE = HEIGHT but no height is specified", advance='yes' )
    case ( noMIF )
      call output ( "TYPE = FIXED but no MIF is specified", advance='yes' )
    case ( noModule )
      call output ( "Instrument module must be specified", advance='yes' )
    case ( noSpacingOrigin )
      call output ( "TYPE = Regular but no spacing and/or origin is specified", &
        & advance='yes' )
    case ( unitlessMessage )
      call output ( "Value for the " )
      call dump_tree_node ( where, 0 )
      call output ( " field is required to dimensionless", advance='yes' )
    end select
    end subroutine ANNOUNCE_ERROR

!=============================================================================
end module HGrid
!=============================================================================

!
! $Log$
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

