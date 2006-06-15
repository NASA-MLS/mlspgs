! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MergeGridsModule

  ! This module contains code for merging operational gridded data with apriori
  ! information.
  ! Secondary operations may be performed directly on the gridded data--
  ! e.g., calculating wmo tropopause pressures from eta-level temperatures

  implicit none
  private

  public :: MergeGrids

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains ! =================================== Public procedures

  ! ----------------------------------------- MergeGrid

  subroutine MergeGrids ( root, griddedDataBase, l2gpDatabase )

    use GriddedData, only: GRIDDEDDATA_T, RGR, &
      & ADDGRIDDEDDATATODATABASE, &
      & CONCATENATEGRIDDEDDATA, CONVERTFROMETALEVELGRIDS, COPYGRID, &
      & NULLIFYGRIDDEDDATA, SETUPNEWGRIDDEDDATA, WRAPGRIDDEDDATA
    use Init_tables_module, only: S_CONCATENATE, S_CONVERTETATOP, &
      & S_DELETE, S_MERGE, S_WMOTROP
    use L2GPData, only: L2GPDATA_T
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
    use MoreTree, only: GET_SPEC_ID
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATE, DECORATION, NODE_ID, SUB_ROSA
    use Tree_Types, only: N_NAMED
    use Toggles, only: GEN, TOGGLE
    use VGridsDatabase, only: AddVGridToDatabase, VGrid_T, VGrids
    use VGrid, only: CREATEVGRIDFROMMLSCFINFO

    integer, intent(in) :: ROOT         ! Tree root
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    type ( l2gpData_T), dimension(:), pointer :: L2GPDATABASE

    ! Local variables
    integer :: I                        ! Loop counter
    integer :: returnStatus
    integer :: SON                      ! Tree node
    integer :: KEY                      ! Another node
    integer :: NAME                     ! Index into string table

    ! excutable code
    if ( toggle(gen) ) call trace_begin ( "MergeGrids", root )
    do i = 2, nsons(root) - 1           ! Skip the begin and end stuff
      son = subtree ( i, root )
      if ( node_id(son) == n_named ) then ! Is spec labed?
        key = subtree ( 2, son )
        name = sub_rosa ( subtree(1,son) )
      else
        key = son
      end if

      select case ( get_spec_id(key) )
      case ( s_merge )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & MergeOneGrid ( key, griddedDataBase ) ) )
      case ( s_concatenate )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & Concatenate ( key, griddedDataBase ) ) )
      case ( s_ConvertEtaToP )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & ConvertEtaToP ( key, griddedDataBase ) ) )
      case ( s_delete )
        call DeleteGriddedData ( key, griddedDatabase )
!       case ( s_vGrid )
!         call decorate ( son, AddVGridToDatabase ( vGrids, &
!           & CreateVGridFromMLSCFInfo ( name, son, l2gpDatabase, returnStatus ) ) )
      case ( s_wmoTrop )
        call decorate ( key, AddgriddedDataToDatabase ( griddedDataBase, &
          & wmoTropFromGrid ( key, griddedDataBase ) ) )
      case default
        ! Shouldn't get here is parser worked?
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Unrecognized command in MergeGrids section' )
      end select
    end do
    if ( toggle(gen) ) call trace_end ( "MergeGrids" )

  end subroutine MergeGrids

  ! ---------------------------------------- ConvertEtaToP --
  type (griddedData_T) function ConvertEtaToP ( root, griddedDataBase ) &
    & result ( newGrid )
    use GriddedData, only: GRIDDEDDATA_T, DUMP, NULLIFYGRIDDEDDATA, &
      & CONVERTFROMETALEVELGRIDS
    use Init_tables_module, only: F_A, F_B, F_GRID
    use output_m, only: output
    use Toggles, only: GEN, TOGGLE
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATION
    ! use VGridsDatabase, only: VGrid_T, VGrids, ConvertVGrid
    
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    ! This routine parses the l2cf instructions that
    ! convert two gridded data on eta surfaces, one of them pressures,
    ! to pressure surfaces

    ! Local variables
    ! Local variables
    integer :: SON                    ! Tree node
    integer :: FIELD                  ! Another tree node
    integer :: FIELD_INDEX            ! Type of tree node
    integer :: VALUE                  ! Tree node
    integer :: I                      ! Loop counter

    type (griddedData_T), pointer :: A ! Temperatures on eta surfaces
    type (griddedData_T), pointer :: B ! Pressures on eta surfaces
    type (griddedData_T), pointer :: V ! Grid with proper pressure surfaces
!    type (VGrid_T), pointer       :: V ! Desired pressure surfaces

    logical, parameter :: DEEBUG = .false.
    ! Executable code
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "ConvertEtaToP", root )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_a ) 
        a => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_b )
        b => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_grid )
        v => griddedDataBase ( decoration ( decoration ( value ) ) )
!       case ( f_VGrid )
!         v => VGrids ( decoration ( decoration ( value ) ) )
      end select
    end do
    if ( DEEBUG ) then
    call output( 'a grid', advance='yes' )
    call dump( a, details=0 )
    call output( 'b grid', advance='yes' )
    call dump( b, details=0 )
    call output( 'v grid', advance='yes' )
    call dump( v, details=0 )
    endif
    call ConvertFromEtaLevelGrids ( a, b, V, newGrid )
    newGrid%sourceFileName      = a%sourceFileName
    newGrid%quantityName        = a%quantityName
    newGrid%description         = a%description
    newGrid%units               = a%units
    newGrid%verticalCoordinate  = v%verticalCoordinate
    newGrid%missingValue        = a%missingValue

  end function ConvertEtaToP

  ! ---------------------------------------- Concatenate --
  function Concatenate ( root, griddedDataBase ) &
    & result ( newGrid )
    use GriddedData, only: GRIDDEDDATA_T, DUMP, &
      & CONCATENATEGRIDDEDDATA, COPYGRID, DestroyGriddedData, NULLIFYGRIDDEDDATA
    use Init_tables_module, only: F_A, F_B, F_GRID
    use Toggles, only: GEN, TOGGLE
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATION
    
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database
    type (griddedData_T), target :: newGrid
    ! This routine parses the l2cf instructions that request
    ! a grid concatenation, and then performs the concatenation

    ! Local variables
    type (griddedData_T), pointer :: A
    type (griddedData_T), pointer :: B
    integer :: db_index
    logical, parameter            :: DEEBUG = .false.
    integer :: FIELD                  ! Another tree node
    integer :: FIELD_INDEX            ! Type of tree node
    integer :: GRIDS_NODE
    integer :: I                      ! Loop counter
    integer :: SON                    ! Tree node
    type (griddedData_T), target :: Intermediate
    integer :: VALUE                  ! Tree node

    ! Executable code
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    call nullifyGriddedData ( Intermediate ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "Concatenate", root )

    ! Get the information from the l2cf    
    grids_node = 0
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_a ) 
        a => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_b )
        b => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_Grid )
        grids_node = son
      end select
    end do

    ! Do the concatenation unless one or other is empty
    if ( grids_node > 0 ) then
      ! Method:
      ! At any step let the result of all prior steps be held in "Intermediate"
      ! Then at each step concatenate the next gridded data with Intermediate
      ! When done, copy Intermediate into result
      do i=2, nsons(grids_node)
        db_index = decoration(decoration(subtree(i, grids_node )))
        b => griddedDataBase ( db_index )
        if ( DEEBUG ) then
          print *, ' '
          print *, 'db_index: ', db_index
          call dump( b, details=-1 )
        endif
        if ( i == 2 ) then
          call CopyGrid ( Intermediate, b )
        else
          call ConcatenateGriddedData ( A, B, Intermediate )
          if ( DEEBUG ) then
            print *, ' '
            print *, 'Result of intermediate concatenate'
            call dump( Intermediate, details=-1 )
          endif
        endif
        call CopyGrid ( newGrid, Intermediate )
        a => newGrid
      enddo
      ! call CopyGrid ( newGrid, Intermediate )
      call DestroyGriddedData ( Intermediate )
      if ( DEEBUG ) call dump( newGrid, details=-1 )
    elseif ( .not. a%empty .and. .not. b%empty ) then
      call ConcatenateGriddedData ( A, B, newGrid )
    else if ( a%empty ) then
      ! Copy B into the result, of course, that may be empty too
      ! in which case the result is empty, no problem!
      call CopyGrid ( newGrid, b )
    else
      ! Otherwise a must be full, b empty
      call CopyGrid ( newGrid, a )
    end if
    newGrid%sourceFileName      = a%sourceFileName
    newGrid%quantityName        = a%quantityName
    newGrid%description         = a%description
    newGrid%units               = a%units
    newGrid%verticalCoordinate  = a%verticalCoordinate
    newGrid%missingValue        = a%missingValue

    if ( toggle(gen) ) call trace_end ( "Concatenate" )
  end function Concatenate

  ! ------------------------------------ DeleteGriddedData ---
  subroutine DeleteGriddedData ( root, griddedDataBase )
    use Tree, only: NSONS, SUBTREE, DECORATION
    use GriddedData, only: DESTROYGRIDDEDDATA, GRIDDEDDATA_T
    use Init_Tables_Module, only: F_GRID
    ! This routine deletes the grid indicated by the l2cf
    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: GRIDDEDDATABASE ! Database
    ! Local variables
    type (griddedData_T), pointer :: GRID
    integer :: FIELD                    ! Tree node
    integer :: FIELD_INDEX              ! Tree node type
    integer :: I                        ! Counter
    integer :: SON                      ! Tree node
    integer :: VALUE                    ! Tree node

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    ! In this case there is only one argument anyway
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_grid ) 
        grid => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do
    call DestroyGriddedData ( grid )
  end subroutine DeleteGriddedData

  ! ----------------------------------------- MergeOneGrid
  type (griddedData_T) function MergeOneGrid ( root, griddedDataBase ) &
    & result ( newGrid )
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use Expr_m, only: EXPR
    use GriddedData, only: GRIDDEDDATA_T, RGR, V_IS_PRESSURE, &
      & COPYGRID, NULLIFYGRIDDEDDATA, &
      & SETUPNEWGRIDDEDDATA, SLICEGRIDDEDDATA, WRAPGRIDDEDDATA
    use Init_tables_module, only: F_CLIMATOLOGY, F_HEIGHT, F_OPERATIONAL, &
      & F_SCALE
    use Intrinsic, only: PHYQ_Length, PHYQ_Pressure
    use L3ASCII, only: L3ASCII_INTERP_FIELD
    use MLSCommon, only: R8
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_ALLOCATE
    use MLSNumerics, only: ESSENTIALLYEQUAL
    use output_m, only: BLANKS, OUTPUT
    use Toggles, only: GEN, TOGGLE
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATION

    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! This routine creates a new grid being a merge of two others.
    ! The operational grid forms the bottom of the dataset
    ! The climatology grid the top.  The result has the horizontal
    ! coordinates to operational and the vertical coordiantes of climatology

    ! Note! This routine is far from efficient. Not least because it
    ! uses l3atascii_interp_field which isn't terribly efficient either.
    ! But hey, this isn't a key part of the software when it comes
    ! to a desire for speed (or at least it shouldn't be)

    ! I'll need to think about missing data at some point.

    ! Local parameters
    real (r8), parameter :: SCALEHEIGHT = 16.0e3_r8 ! Approximate scale height / m

    ! Local variables
    integer :: DAY                      ! Loop counter
    logical, parameter :: DEEBUG = .false.
    integer :: FIELD                    ! Another tree node
    integer :: FIELD_INDEX              ! Type of tree node
    integer :: I                        ! Loop inductor
    integer :: LAT                      ! Loop counter
    integer :: LON                      ! Loop counter
    integer :: LST                      ! Loop counter
    integer :: SON                      ! Tree node
    integer :: STATUS                   ! Flag from allocate
    integer :: SURF                     ! Loop counter
    integer :: SZA                      ! Loop counter
    integer :: EXPRUNITS(2)             ! Units for expr
    integer :: VALUE                    ! Value of tree node

    real (r8) :: CLIWEIGHT              ! Climatological 'weight'
    real (r8) :: HEIGHT                 ! Transition height
    real (r8) :: OPWEIGHT               ! Operational 'weight'
    real (r8) :: SCALE                  ! Transition scale
    real (r8) :: TOTALWEIGHT            ! Total weight
    real (r8) :: EXPRVALUES(2)          ! Value of expr
    real (r8) :: ZTRANS                 ! Transition 'height'
    real (r8) :: Z                      ! One 'height'
    real (r8) :: Z1, Z2                 ! Range of transition region

    real (rgr) :: CLIVAL                ! One interpolated value
    real (rgr) :: OPVAL                 ! One interpolated value
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: CLIMAPPED
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: OPERMAPPED
    real (r8), dimension(:), pointer :: MEANDATES ! Mean dates for new grid

    type (griddedData_T), pointer :: OPERATIONAL
    type (griddedData_T), pointer :: CLIMATOLOGY

    ! Executable code

    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "MergeOneGrid", root )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_operational ) 
        operational => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_climatology )
        climatology => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_height )
        call expr ( value, exprUnits, exprValues )
        if ( exprUnits(1) /= phyq_pressure ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only pressure is allowed for the height field' )
        height = exprValues(1)
      case ( f_scale )
        call expr ( value, exprUnits, exprValues )
        if ( exprUnits(1) /= phyq_length ) call MLSMessage ( &
          & MLSMSG_Error, ModuleName, &
          & 'Only altitude is allowed for the scale field' )
        scale = exprValues(1)
      end select
    end do

    ! Think about cases where one or other grid is empty
    if ( climatology%empty ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'The climatology grid for the merge is empty' )
    if ( operational%empty ) then
      ! If no operational data, then just use climatology
      call CopyGrid ( newGrid, climatology )
      return
    end if
    if ( DEEBUG ) then
    call output( 'operational%verticalCoordinate: ', advance='no' )
    call output( operational%verticalCoordinate, advance='no' )
    call blanks(3)
    call output( v_is_pressure, advance='yes' )
    call output( 'climatology%verticalCoordinate: ', advance='no' )
    call output( climatology%verticalCoordinate, advance='no' )
    call blanks(3)
    call output( v_is_pressure, advance='yes' )
    endif
    ! Do some final sanity checks
    if ( operational%verticalCoordinate /= v_is_pressure ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Operational grid not on pressure surfaces' )
    if ( climatology%verticalCoordinate /= v_is_pressure ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Climatology grid not on pressure surfaces' )
    !     if ( climatology%units /= operational%units ) &
    !       & call MLSMessage ( MLSMSG_Error, ModuleName, &
    !       & 'The climatology and operational data describe different physical quantities' )
    if ( climatology%equivalentLatitude .neqv. operational%equivalentLatitude ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'The climatology and operational data are mixed latitude/equivalent latitude.' )

    ! OK, now we're ready to go.
    ! First we're going to 'wrap' the climatology to be sure that we can
    ! interpolate it in longitude.  The chances are that it has no longitudinal
    ! variation anyway, so this won't actually do anything
    call WrapGriddedData ( climatology )
    ! Create the result.  It has the same vertical coordinates as climatology
    ! But the same horizontal coordinates as operational.
    call SetupNewGriddedData ( newGrid, source=operational, &
      & noHeights=climatology%noHeights )
    ! Setup the rest of the quantity
    newGrid%sourceFileName = 'Result of merge'
    newGrid%quantityName   = 'Result of merge'
    newGrid%description    = 'Result of merge'
    newGrid%units          = climatology%units
    newGrid%verticalCoordinate = v_is_pressure
    newGrid%equivalentLatitude = climatology%equivalentLatitude
    newGrid%heights = climatology%heights
    newGrid%lats = operational%lats
    newGrid%lons = operational%lons
    newGrid%lsts = operational%lsts
    newGrid%szas = operational%szas
    newGrid%dateStarts = operational%dateStarts
    newGrid%dateEnds = operational%dateEnds

    ! Get the 'mean' dates for the result
    nullify ( meanDates )
    call Allocate_test ( meanDates, newGrid%noDates, 'meanDates', ModuleName )
    meanDates = ( newGrid%dateStarts + newGrid%dateEnds ) / 2.0

    ! Now create two fields the same shape as the new field that contain
    ! the operational and climatological data interpolated to our new locations.
    allocate ( operMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped' )
    allocate ( cliMapped ( &
      & newGrid%noHeights, newGrid%noLats, newGrid%noLons, &
      & newGrid%noLsts, newGrid%noSzas, newGrid%noDates ), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped' )

    call SliceGriddedData ( operational, operMapped, &
      & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
      & newGrid%szas, meanDates, missingValue=newGrid%missingValue )
    call SliceGriddedData ( climatology, cliMapped, &
      & newGrid%heights, newGrid%lats, newGrid%lons, newGrid%lsts, &
      & newGrid%szas, meanDates, missingValue=newGrid%missingValue )

    zTrans = scaleHeight * ( 3.0 - log10 ( height ) )
    z1 = zTrans - scale/2.0
    z2 = zTrans + scale/2.0

    ! Now we're going to fill in the rest of the field
    do day = 1, newGrid%noDates
      do sza = 1, newGrid%noSzas
        do lst = 1, newGrid%noLsts
          do lon = 1, newGrid%noLons
            do lat = 1, newGrid%noLats
              do surf = 1, newGrid%noHeights
                ! Get the values
                cliVal = cliMapped ( surf, lat, lon, lst, sza, day )
                opVal = operMapped ( surf, lat, lon, lst, sza, day )
                ! Weight them by height
                z = scaleHeight * ( 3.0 - log10 ( newGrid%heights(surf) ) )
                if ( scale /= 0.0 ) then
                  cliWeight = ( z - z1 ) / ( z2-z1 )
                  opWeight = 1.0 - cliWeight
                end if

                cliWeight = min ( max ( cliWeight, 0.0_r8 ), 1.0_r8 )
                opWeight = min ( max ( opWeight, 0.0_r8 ), 1.0_r8 )

                ! Check for bad data in operational dataset
                if ( EssentiallyEqual ( opVal, newGrid%missingValue ) ) opWeight = 0.0

                ! Check for bad data in the climatology
                if ( EssentiallyEqual ( cliVal, newGrid%missingValue ) ) &
                  & call MLSMessage ( MLSMSG_Error, ModuleName, &
                  & 'There is a bad data point in the climatology field' )

                totalWeight = cliWeight + opWeight
                if ( totalWeight == 0.0 ) then
                  ! Presumably was in the region where operational was supposed
                  ! to dominate, but it's bad, so switch to a priori
                  cliWeight = 1.0
                  totalWeight = 1.0
                end if

                ! OK, store this value
                newGrid%field ( surf, lat, lon, lst, sza, day ) = &
                  & ( cliWeight*cliVal + opWeight*opVal ) / totalWeight
              end do
            end do
          end do
        end do
      end do
    end do

    ! Tidy up
    call Deallocate_test ( meanDates, 'meanDates', ModuleName )
    deallocate ( cliMapped, operMapped, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'operMapped or cliMapped' )

    if ( toggle(gen) ) call trace_end ( "MergeOneGrid" )
  end function MergeOneGrid

  ! ----------------------------------------- wmoTropFromGrid
  type (griddedData_T) function wmoTropFromGrid ( root, griddedDataBase ) &
    & result ( newGrid )
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use GriddedData, only: GRIDDEDDATA_T, DUMP, RGR, V_IS_PRESSURE, V_IS_ETA, &
      & COPYGRID, NULLIFYGRIDDEDDATA, &
      & DOGRIDDEDDATAMATCH, &
      & SETUPNEWGRIDDEDDATA, SLICEGRIDDEDDATA, WRAPGRIDDEDDATA
    use Init_tables_module, only: F_A, F_B, F_GRID
    use MLSCommon, only: DEFAULTUNDEFINEDVALUE
    use MLSFillValues, only: IsFillValue, RemoveFillValues
    use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR, MLSMSG_WARNING
    use MLSStrings, only: LOWERCASE
    use output_m, only: output
    use Toggles, only: GEN, TOGGLE
    use Trace_M, only: TRACE_BEGIN, TRACE_END
    use Tree, only: NSONS, SUBTREE, DECORATION
    use wmoTropopause, only: ExtraTropics, twmo
    ! Implements the algorithm published in GRL

    integer, intent(in) :: ROOT         ! Tree node
    type (griddedData_T), dimension(:), pointer :: griddedDataBase ! Database

    ! This routine creates a new gridded data by finding wmo Tropopause
    ! pressure levels among the temperatures of another gridded data,
    ! possibly using the vertical coords of that gridded data
    ! (if entered via the "grid=" field) 
    ! or else using the corresponding pressures stored
    ! as the field values of a second gridded data
    ! (if entered via "a=first_gridded, b=second_gridded")
    ! The new gridded has only one "level" per horizontal grid point
    ! with tropopause pressures stored in the values field
    
    ! We'll assume some things:
    ! (1) verticalCoordinate is either pressure or eta
    ! (2) If Pressures grid supplied, grids match
    !     and Pressures units either Pa or hPa
    
    ! Local variables
    integer :: field
    integer :: field_index
    real, dimension(:), pointer :: h ! hPa
    integer :: i
    integer :: iDate
    integer :: iLst
    integer :: invert
    integer :: iSza
    integer :: lat
    integer :: lon
    real :: missingValue
    integer :: nLev
    integer :: nValid
    real, dimension(:), pointer :: p ! Pa
    real, parameter :: pliml = 65.*100 ! in Pa
    real, parameter :: plimlex = 65.*100 ! in Pa
    real, parameter :: plimu = 550.*100 ! in Pa
    type (griddedData_T), pointer :: Placeholder    => null()
    type (griddedData_T), pointer :: Pressures    => null()
    integer :: son
    real, dimension(:), pointer :: t
    type (griddedData_T), pointer :: Temperatures => null()
    real :: scale
    real :: trp
    integer :: value
    real, dimension(:), pointer :: xyTemp, xyPress

    ! Executable code
    call nullifyGriddedData ( newGrid ) ! for Sun's still useless compiler
    if ( toggle(gen) ) call trace_begin ( "MergeOneGrid", root )
    MISSINGVALUE = REAL( DEFAULTUNDEFINEDVALUE )

    ! Get the information from the l2cf    
    ! Note that init_tables_module has insisted that we have all
    ! arguments so we don't need a 'got' type arrangement
    do i = 2, nsons(root)
      son = subtree(i,root)
      field = subtree(1,son)
      value = subtree(2,son)
      field_index = decoration(field)
      select case ( field_index )
      case ( f_a ) 
        Temperatures => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_b ) 
        Pressures    => griddedDataBase ( decoration ( decoration ( value ) ) )
      case ( f_grid ) 
        Temperatures => griddedDataBase ( decoration ( decoration ( value ) ) )
      end select
    end do

    if ( associated(Temperatures) .and. associated(Pressures) ) then
    ! What if you reversed the sense of a (Temperatures) and b (Pressures)?
    ! We must switch them
      if ( index( 'hpa,mb', lowercase(trim(Temperatures%units)) ) > 0 .and. &
        &  index( 'hpa,mb', lowercase(trim(Pressures%units)) ) < 1 ) then
        Placeholder  => Pressures
        Pressures    => Temperatures
        Temperatures => Placeholder
      endif
    ! What if Temperatures and Pressures don't match
      if ( .not. doGriddeddataMatch( Temperatures, Pressures ) ) then
        call output( 'Vert. coords match? ', advance='no' )
        call output( Temperatures%verticalCoordinate==Pressures%verticalcoordinate, &
          & advance='yes' )
        call output( 'starting dates match? ', advance='no' )
        call output( any( Temperatures%DateStarts==Pressures%DateStarts ), &
          & advance='yes' )
        call output( 'ending dates match? ', advance='no' )
        call output( any( Temperatures%DateEnds==Pressures%DateEnds ), &
          & advance='yes' )
        call output( 'Lsts match? ', advance='no' )
        call output( any( Temperatures%Lsts==Pressures%Lsts ), &
          & advance='yes' )
        call output( 'Szas match? ', advance='no' )
        call output( any( Temperatures%Szas==Pressures%Szas ), &
          & advance='yes' )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Gridded T,P data must match to calculate wmo Tropopause' )
      endif
    endif
    call output( 'Temperatures grid', advance='yes' )
    call dump( Temperatures, details=0 )
    call output( 'Pressures grid', advance='yes' )
    call dump( Pressures, details=0 )
    if ( .not. associated(Temperatures) ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'No associated Temperatures grid for calculating wmo tropopause' )
      return
    endif
    nlev = Temperatures%noHeights
    if ( nlev < 2 ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'Too few levels on Temperatures grid for calculating wmo tropopause' )
      return
    endif
    ! Right now we can't read eta levels, only pressures
    ! but when we move to GEOS5 GMAO we'll have no choice:
    ! Must read eta-level files
    if ( .not. any( &
      & Temperatures%verticalCoordinate == (/ v_is_pressure, v_is_eta /) ) ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Temperatures must be on eta or pressure calculating wmo tropopause' )
    if (  Temperatures%verticalCoordinate /= v_is_pressure  .and. &
      & .not. associated(Pressures) ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Temperatures illegal verticalcoordinate calculating wmo tropopause' )
    ! For now, just crudely assume the heights is in units of Pa
    ! If not, we'll need to check heightsUnits
    if ( lowercase(Temperatures%heightsUnits) /= 'pa' .and. &
      & .not. associated(Pressures) ) &
      & call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Temperatures illegal heightsUnits calculating wmo tropopause' )

    call Allocate_test (h, nlev, 'h', ModuleName )
    call Allocate_test (p, nlev, 'p', ModuleName )
    call Allocate_test (t, nlev, 't', ModuleName )
    ! check vertical orientation of data
    ! (twmo expects ordered from top downward)

    h = Temperatures%heights
    if (h(1) .gt. h(2)) then
      invert=1
      p = h(nlev:1:-1)*100.  ! hPa > Pa
    else
      invert=0
      p = h(:)*100.         ! hPa > Pa
    endif
    
    call SetupNewGriddedData ( newGrid, source=Temperatures, &
      & noHeights=1, noDates=1 )
    ! Setup the rest of the quantity
    newGrid%sourceFileName     = 'Gridded Temperatures'
    newGrid%quantityName       = 'wmo Tropopause'
    newGrid%description        = 'wmo Tropopause'
    newGrid%units              = 'hPa' ! If we want 'Pa', restore /100 below
    newGrid%verticalCoordinate = v_is_pressure
    newGrid%equivalentLatitude = Temperatures%equivalentLatitude
    newGrid%heights            = Temperatures%missingValue
    newGrid%lats               = Temperatures%lats
    newGrid%lons               = Temperatures%lons
    newGrid%lsts               = Temperatures%lsts
    newGrid%szas               = Temperatures%szas
    newGrid%dateEnds           = Temperatures%dateEnds(1)
    newGrid%dateStarts         = Temperatures%dateStarts(1)
    newGrid%field              = MISSINGVALUE
    ! Now actually calculate the tropopause
    ! for every "horizontal" point
    do idate=1, 1 ! size( Temperatures%field, 6 )
      do iSza=1, size( Temperatures%field, 5 )
        do iLst=1, size( Temperatures%field, 4 )
          do lon=1, size( Temperatures%field, 3 )
            do lat=1, size( Temperatures%field, 2 )
              ! Do we have a second, pressure, gridded data?
              if ( associated(Pressures) ) then
                select case (lowercase(Pressures%units))
                case ('pa', 'b')
                  scale = 100. ! To convert Pa to hPa
                case ('hpa', 'mb')
                  scale = 1.
                case default
                  call output( 'Pressures%units: ', advance='no' )
                  call output( trim(Pressures%units), advance='yes' )
                  call MLSMessage ( MLSMSG_Error, moduleName, &
                    & 'Pressures units must be Pa, hPa, or mb calculating wmo tropopause' )
                end select
                h = Pressures%field(nlev:1:-1,lat,lon,iLst,iSza,idate)
                if (h(1) .gt. h(2)) then
                  invert=1
                  p = Pressures%field(nlev:1:-1,lat,lon,iLst,iSza,idate) * scale
                else
                  invert=0
                  p = Pressures%field(nlev:1:-1,lat,lon,iLst,iSza,idate) * scale
                endif
              endif
              if ( invert == 1 ) then
                t = temperatures%field(nlev:1:-1,lat,lon,iLst,iSza,idate)
              else
                t = temperatures%field(:,lat,lon,iLst,iSza,idate)
              endif
              where ( t < 0. .or. t > 100000. )
                t = MissingValue
              end where
              nvalid = count( .not. isFillValue(t) )
              if ( nvalid < 2 ) cycle
              call Allocate_test (xyTemp, nvalid, 'xyTemp', ModuleName )
              call Allocate_test (xyPress, nvalid, 'xyPress', ModuleName )
              call RemoveFillValues( t, MISSINGVALUE, xyTemp, &
                & p, xyPress )
              call twmo(nvalid, xyTemp, xyPress, plimu, pliml, trp)
              ! Don't let tropopause sink too low in "extra tropics"
              if ( trp < plimlex .and. &
                & extraTropics(temperatures%lats(lat)) )  &
                & trp = MISSINGVALUE
              if ( trp > 0. .and. trp < 100000000. ) &
                & newGrid%field(1, lat,lon,iLst,iSza,idate) = trp ! /100 for 'Pa'
              call Deallocate_test ( xyTemp, 'xyTemp', ModuleName )
              call Deallocate_test ( xyPress, 'xyPress', ModuleName )
            enddo ! Lats
          enddo ! Lons
        enddo ! Lsts
      enddo ! Szas
    enddo ! dates
    call Deallocate_test ( h, 'h', ModuleName )
    call Deallocate_test ( p, 'p', ModuleName )
    call Deallocate_test ( t, 't', ModuleName )
  end function wmoTropFromGrid

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MergeGridsModule

! $Log$
! Revision 2.21  2006/06/15 00:02:33  pwagner
! Should work with geos5: convert then concatenate
!
! Revision 2.20  2006/06/13 22:13:12  pwagner
! changed interface to ConvertFromEtaLevelGrids
!
! Revision 2.19  2006/05/12 21:26:37  pwagner
! Added extra debugging statements
!
! Revision 2.18  2006/05/09 16:42:02  pwagner
! May find wmo p trop with two eta-level grids
!
! Revision 2.17  2006/05/04 23:04:59  pwagner
! May convertEtaToP and create a VGrid in MergeGrids section
!
! Revision 2.16  2006/02/11 00:14:08  pwagner
! May calculate wmoTropopause in this section directly
!
! Revision 2.15  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units
!
! Revision 2.13  2003/06/06 01:06:59  livesey
! Added DeleteGrids stuff
!
! Revision 2.12  2003/05/09 02:13:05  livesey
! Removed a dump
!
! Revision 2.11  2003/05/09 01:55:14  livesey
! Sped up Merge by using new SliceGriddedData routine.
!
! Revision 2.10  2003/04/04 00:08:26  livesey
! Added Concatenate capability, various reorganizations.
!
! Revision 2.9  2003/02/28 02:33:28  livesey
! Bug fix, careless with the old emacs.
!
! Revision 2.8  2003/02/28 02:25:50  livesey
! First working version.
!
! Revision 2.7  2003/02/19 19:15:13  pwagner
! Consistent with new GriddedData_T and rgr
!
! Revision 2.6  2002/11/22 12:21:14  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.5  2002/10/08 17:36:21  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.4  2002/08/22 20:26:09  vsnyder
! Move another USE from module scope to procedure scope
!
! Revision 2.3  2002/08/21 02:23:39  vsnyder
! Move USE statements from module scope to procedure scope
!
! Revision 2.2  2002/01/26 00:10:54  livesey
! It compiles at least
!
! Revision 2.1  2002/01/24 00:58:03  livesey
! First version, not much more than a stub
!
