! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GriddedData ! Contains the derived TYPE GriddedData_T

  use MLSCommon, only: RGR=>R4, R8, LINELEN, NAMELEN, DEFAULTUNDEFINEDVALUE
  ! r4 corresponds to sing. prec. :: same as stored in files
  ! (except for dao dimensions)
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_ERROR, &
    & MLSMSG_DEALLOCATE
  use Output_m, only: OUTPUT
  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST

  implicit NONE
  private

  !------------------------------- RCS Ident Info ------------------------------
  character(LEN=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)), private :: Id = idParm
  character(LEN=*), parameter :: ModuleName="$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------

  public :: AddGriddedDataToDatabase, ConcatenateGriddedData, CopyGrid, &
    & DestroyGriddedData, DestroyGriddedDataDatabase, Dump, GriddedData_T, &
    & NullifyGriddedData, RGR, SetupNewGriddedData, SliceGriddedData, &
    & WrapGriddedData

  logical, private, parameter :: MAYDUMPFIELDVALUES = .true.

  interface DUMP
    module procedure DumpGriddedData
    module procedure DumpGriddedDatabase
  end interface

  ! These are 'enumerated types' consistent with hph's
  ! work in l3ascii_read_field
  public :: V_is_pressure, V_is_altitude, V_is_GPH, V_is_theta

  integer, parameter :: V_is_pressure = 1
  integer, parameter :: V_is_altitude = v_is_pressure+1
  integer, parameter :: V_is_GPH = v_is_altitude+1
  integer, parameter :: V_is_theta = v_is_gph+1
  
  ! If dumping gridded data, always give some details of any matching these
  character(len=*), parameter :: ALWAYSDUMPTHESE = 'dao,ncep' ! -> ' '
  ! and for these automatically dumped ones, this level of detail for multi-dim
  integer, parameter :: AUTOMATICDETAILS = 0 ! 1 means dump, 0 means no

  ! This type reflects the format of the Level 3 ASCII files, though note that
  ! these files can store multiple quantities such as these.

  type GriddedData_T

    ! First the comment line(s) from the relevant input file
    logical :: EMPTY                    ! Set for example when file read failed
    character (LEN=LineLen), pointer, dimension(:) :: fileComments => NULL()

    ! Now the name, description and units information
    character (LEN=LineLen) :: sourceFileName ! Input file name
    character (LEN=NameLen) :: quantityName ! From input file
    character (LEN=LineLen) :: description ! Quantity description
    character (LEN=NameLen) :: units ! Units for quantity

    ! Now define the various coordinate systems, first vertical
    integer :: verticalCoordinate ! An 'enumerated' type
    integer :: noHeights         ! Number of surfaces
    real (rgr), pointer, dimension(:) :: heights  => NULL()
    ! Surfaces (e.g. pressures etc.) [noHeights]

    ! Now the latitudinal coordinate
    logical :: equivalentLatitude       ! If set, coordinate is equivalent latitude
    logical :: noYear                   ! If set, field is for any year
    integer :: noLats                   ! Number of latitudes
    real (rgr), pointer, dimension(:) :: Lats => NULL() ! Latitudes [noLats]
    integer :: noLons                   ! Number of longitudes
    real (rgr), pointer, dimension(:) :: Lons => NULL() ! Longitudes [noLons]
    integer noLsts                      ! Number of local times
    real (rgr), pointer, dimension(:) :: Lsts => NULL() ! Local times [noLsts]
    integer noSzas                      ! Number of solar zenith angles
    real (rgr), pointer, dimension(:) :: Szas => NULL() ! Zenith angles [noSzas]
    integer noDates                     ! Number of dates in data
    real (r8), pointer, dimension(:) :: DateStarts => NULL()
    ! Starting dates in SDP toolkit format
    real (r8), pointer, dimension(:) :: DateEnds => NULL()
    ! Ending dates in SDP toolkit format

    ! The data itself.  This is stored as
    !  [noHeights, noLats, noLons, noLsts, noSzas, noDates]
    real (rgr), pointer, dimension(:,:,:,:,:,:) :: field => NULL()
    real (rgr) :: missingValue

  end type GriddedData_T

! ======================================================================
contains

  ! -----------------------------------  AddGriddedDataToDatabase  -----
  integer function AddGriddedDataToDatabase ( Database, Item )
  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

    ! Dummy arguments
    type (GriddedData_T), dimension(:), pointer :: Database
    type (GriddedData_T), intent(in) :: Item

    ! Local variables
    type (GriddedData_T), dimension(:), pointer :: tempDatabase

    ! Executable code

    include "addItemToDatabase.f9h"
    AddGriddedDataToDatabase = newSize

  end function AddGriddedDataToDatabase

  ! -------------------------------------------- ConcatenateGriddedData
  subroutine ConcatenateGriddedData ( A, B, X )
    ! This routine takes two grids A and B, B dated after A and tries
    ! to produce a third grid which is a combination of A and B
    use MLSNumerics, only: EssentiallyEqual
    type ( GriddedData_T ), intent(in) :: A
    type ( GriddedData_T ), intent(in) :: B
    type ( GriddedData_T ), intent(inout) :: X ! inout to let us deallocate it
    ! Local variables

    ! Executable code
    ! First, check that the grids A and B are conformable.
    if ( a%verticalCoordinate /= b%verticalCoordinate .or. &
      & a%equivalentLatitude .neqv. b%equivalentLatitude .or. &
      & a%noYear .neqv. b%noYear ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate are not compatible' )
    if ( .not. all ( EssentiallyEqual ( a%heights, b%heights ) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate do not share heights' )
    if ( .not. all ( EssentiallyEqual ( a%lats, b%lats ) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate do not share lats' )
    if ( .not. all ( EssentiallyEqual ( a%lons, b%lons ) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate do not share lons' )
    if ( .not. all ( EssentiallyEqual ( a%lsts, b%lsts ) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate do not share lsts' )
    if ( .not. all ( EssentiallyEqual ( a%szas, b%szas ) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate do not share szas' )
    if ( .not. EssentiallyEqual ( a%missingValue, b%missingValue ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for Concatenate have different missing values' )

    ! Check that the dates are going to work out
    if ( maxval ( a%dateEnds ) > minval ( b%dateStarts ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for concatenation are not in correct time order' )

    ! OK, now we're ready
    call DestroyGriddedData ( X )
    call SetupNewGriddedData ( X, source=A, noDates= a%noDates + b%noDates )
    ! Copy over the unchanged position data
    x%heights = a%heights
    x%lats = a%lats
    x%lons = a%lons
    x%lsts = a%lsts
    x%szas = a%szas
    x%dateStarts = (/ a%dateStarts, b%dateStarts /)
    x%dateEnds = (/ a%dateEnds, b%dateEnds /) 

    x%field ( :, :, :, :, :, 1:a%noDates ) = a%field
    x%field ( :, :, :, :, :, a%noDates+1:x%noDates ) = b%field
  end subroutine ConcatenateGriddedData

  ! ---------------------------------------- CopyGrid ------------------
  subroutine CopyGrid ( Z, X )
    ! Does a deep copy of X into Z
    type ( GriddedData_T ), intent(inout) :: Z
    type ( GriddedData_T ), intent(in) :: X
    ! Executable code
    call DestroyGriddedData ( Z )
    call SetupNewGriddedData ( Z, source=X )
    ! Copy the information over
    z%heights = x%heights
    z%lats = x%lats
    z%lons = x%lons
    z%lsts = x%lsts
    z%szas = x%szas
    z%dateStarts = x%dateStarts
    z%dateEnds = x%dateEnds
    z%field = x%field
  end subroutine CopyGrid

  ! -----------------------------------------  DestroyGriddedData  -----
  subroutine DestroyGriddedData ( Qty )
    ! This subroutine destroys a quantity template

    ! Dummy argument
    type (GriddedData_T), intent(INOUT) :: Qty

    ! Local variables
    integer :: STATUS

    ! Executable code
    call Deallocate_test ( qty%heights, "qty%heights", ModuleName )
    call Deallocate_test ( qty%lats, "qty%lats", ModuleName )
    call Deallocate_test ( qty%lons, "qty%lons", ModuleName )
    call Deallocate_test ( qty%lsts, "qty%lsts", ModuleName )
    call Deallocate_test ( qty%szas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Deallocate_test ( qty%dateStarts, "qty%dateStarts", ModuleName )
    call Deallocate_test ( qty%dateEnds, "qty%dateEnds", ModuleName )

    ! Now the data itself
    if ( associated(qty%field) ) then
      deallocate(qty%field, STAT=status)

      if (status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//"qty%field")
    end if

  end subroutine DestroyGriddedData

  ! ----------------------------------  DestroyGriddedDataDatabase -----
  subroutine DestroyGriddedDataDatabase ( Database )
  ! This subroutine destroys a quantity template database

    use Toggles, only: TOGGLE, GEN
    use Trace_M, only: TRACE_BEGIN, TRACE_END

    ! Dummy argument
    type (GriddedData_T), dimension(:), pointer :: Database

    ! Local variables
    integer :: qtyIndex, status

    if ( toggle(gen) ) call trace_begin ( "DestroyGriddedDataDatabase" )

    if (associated(database)) then
      do qtyIndex=1,size(database)
        call DestroyGriddedData(database(qtyIndex))
      enddo
      deallocate(database, stat=status)
      if (status /= 0) call MLSMessage ( MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate//"database")
    endif
    if ( toggle(gen) ) then
      call trace_end ( "DestroyGriddedDataDatabase" )
    end if
  end subroutine DestroyGriddedDataDatabase

  ! ----------------------------------------  DumpGriddedDatabase  -----
  subroutine DumpGriddedDatabase ( GriddedData, Details )

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao

    type (GriddedData_T), dimension(:), pointer :: GriddedData 

    integer, intent(in), optional :: DETAILS
    
    ! Local Variables
    integer            :: i

    if ( .not. associated(GriddedData)) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Gridded database still null')

    call output ( '============ Gridded Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    call output ( 'database: a priori grids: SIZE = ' )
    call output ( size(GriddedData), advance='yes' )
    if ( size(GriddedData) < 1 ) return
    do i = 1, size(GriddedData)

      call output ( 'item number ' )
      call output ( i, advance='yes' )

      call DumpGriddedData(GriddedData(i), Details)
    end do ! i
  end subroutine DumpGriddedDatabase

  ! --------------------------------------------  DumpGriddedData  -----
  subroutine DumpGriddedData ( GriddedData, Details )
    use Dump_0, only: Dump

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao
    type (GriddedData_T) :: GriddedData 
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0

    ! Local Variables
    integer :: MYDETAILS
    integer :: FIELDVALUESDETAILS

    ! Executable code
    myDetails = 0
    if ( present(details) ) myDetails = details
    fieldvaluesdetails = myDetails
    if ( index(ALWAYSDUMPTHESE, trim(GriddedData%description)) > 0 ) then
      myDetails = 1
      fieldvaluesdetails = max(fieldvaluesdetails, AUTOMATICDETAILS)
    endif
    if ( GriddedData%empty ) then
      call output('This Gridded quantity was empty (perhaps the file name' &
        & // ' was wrong)', advance='yes')
      return
    endif
    call output('Gridded quantity name ' // GriddedData%quantityName, advance='yes')
      if ( myDetails < -1 ) return
    call output('description ' // GriddedData%description, advance='yes')
    call output('units ' // GriddedData%units, advance='yes')
    call output('missing value ', advance='no')
    call output(GriddedData%missingValue, advance='yes')

    call output ( ' ************ Geometry ********** ' ,advance='yes')

    call output ( ' Vertical coordinate = ' )
    call output ( GriddedData%verticalCoordinate, advance='yes' )
    call output ( ' No. of heights = ' )
    call output ( GriddedData%noHeights, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%heights, &
      & '    Heights =' )

    call output ( ' Equivalent latitude = ' )
    call output ( GriddedData%equivalentLatitude, advance='yes' )
    call output ( ' No. of latitudes = ' )
    call output ( GriddedData%noLats, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lats, &
      & '    latitudes =' )

    call output ( ' No. of longitudes = ' )
    call output ( GriddedData%noLons, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lons, &
      & '    longitudes =' )

    call output ( ' No. of local times = ' )
    call output ( GriddedData%noLsts, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%lsts, &
      & '    local times =' )

    call output ( ' No. of solar zenith angles = ' )
    call output ( GriddedData%noSzas, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%szas, &
      & '    solar zenith angles =' )

    call output ( ' No. of dates = ' )
    call output ( GriddedData%noDates, advance='yes' )
    if ( myDetails >= 0 ) call dump ( GriddedData%dateStarts, &
      & '    starting dates =' )
    if ( myDetails >= 0 ) call dump ( GriddedData%dateEnds, &
      & '    ending dates =' )

    if ( MAYDUMPFIELDVALUES .and. fieldvaluesdetails > 0 ) then
      call output ( ' ************ tabulated field values ********** ' ,advance='yes')
     ! May dump a 3-d slice of 6-d array
      if ( GriddedData%noDates == 1 .and. GriddedData%noSzas == 1 &
        & .and. GriddedData%noLsts == 1 ) then
        call dump ( GriddedData%field(:,:,:,1,1,1), &
          & '    gridded field values =' )
      elseif ( GriddedData%noDates == 1 .and. GriddedData%noSzas == 1 ) then
        call dump ( GriddedData%field(:,:,:,1,1,1), &
          & '    gridded field values (1st solar time) =' )
        call dump ( GriddedData%field(:,:,1,:,1,1), &
          & '    gridded field values (1st longitude) =' )
        call dump ( GriddedData%field(:,1,:,:,1,1), &
          & '    gridded field values (1st latitude) =' )
        call dump ( GriddedData%field(1,:,:,:,1,1), &
          & '    gridded field values (1st height) =' )
      else
        ! No dump for 6-dimensional double arrays yet, anyway
        !     call dump ( GriddedData%field, &
        !      & '    gridded field values =' )
        call output ( ' *(Sorry, dump_6d_... not yet coded)* ' , &
          & advance='yes')
      endif
    endif

  end subroutine DumpGriddedData

  ! ----------------------------------------  SetupNewGriddedData  -----
  subroutine SetupNewGriddedData ( Qty, Source, NoHeights, NoLats, &
    & NoLons, NoLsts, NoSzas, NoDates, missingValue, &
    & Empty )
  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.
    ! Dummy arguments
    type (GriddedData_T) :: QTY ! Result
    type (GriddedData_T), optional, intent(in) :: SOURCE ! Template
    integer, optional, intent(in) :: NOHEIGHTS, NOLATS, NOLONS, NOLSTS, NOSZAS, NODATES
    logical, optional, intent(in) :: EMPTY
    real(rgr), optional, intent(in) :: missingValue
    ! Local parameters
    real(rgr), parameter :: DefaultMissingValue = DEFAULTUNDEFINEDVALUE ! -999.99
    ! Local variables
    integer :: status           ! Status from allocates etc.
    logical :: myEmpty                  ! Copy of empty possibly

    ! Executable code

    myEmpty = .false.
    if ( present ( empty ) ) myEmpty = empty
    if ( present ( source ) .and. .not. present ( empty )  ) myEmpty = source%empty
    if ( myEmpty ) then
      qty%empty = .true.
      qty%noHeights = 0
      qty%noLats = 0
      qty%noLons = 0
      qty%noLsts = 0
      qty%noSzas = 0
      qty%noDates = 0
    else
      qty%empty = .false.
      
      ! Set sensible values for some parameters
      ! First, if we have a template setup according to that
      if (present(source)) then
        qty%noHeights = source%noHeights
        qty%noLats = source%noLats
        qty%noLons = source%noLons
        qty%noLsts = source%noLsts
        qty%noSzas = source%noSzas
        qty%noDates = source%noDates
        qty%equivalentLatitude = source%equivalentLatitude
        qty%noYear = source%noYear
        qty%missingValue = source%missingValue
        qty%verticalCoordinate = source%verticalCoordinate
      else ! We have no template, setup a very bare quantity
        qty%noHeights = 1
        qty%noLats = 1
        qty%noLons = 1
        qty%noLsts = 1
        qty%noSzas = 1
        qty%noDates = 1
        qty%equivalentLatitude = .false.
        qty%noYear = .false.
        qty%missingValue = defaultMissingValue
      endif
      
      ! Now, see if the user asked for modifications to this
      if (present(noHeights)) qty%noHeights = noHeights
      if (present(noLats)) qty%noLats = noLats
      if (present(noLons)) qty%noLons = noLons
      if (present(noLsts)) qty%noLsts = noLsts
      if (present(noSzas)) qty%noSzas = noSzas
      if (present(noDates)) qty%noDates = noDates
      if ( present ( missingValue ) ) qty%missingValue = missingValue
    end if

    ! First the vertical/horizontal coordinates
    call Allocate_test ( qty%heights, qty%noHeights, "qty%heights", ModuleName )
    call Allocate_test ( qty%lats, qty%noLats, "qty%lats", ModuleName )
    call Allocate_test ( qty%lons, qty%noLons, "qty%lons", ModuleName )
    call Allocate_test ( qty%lsts, qty%noLsts, "qty%lsts", ModuleName )
    call Allocate_test ( qty%szas, qty%noSzas, "qty%szas", ModuleName )

    ! Now the temporal coordinates
    call Allocate_test ( qty%dateStarts, qty%noDates, "qty%dateStarts", ModuleName )
    call Allocate_test ( qty%dateEnds, qty%noDates, "qty%dateEnds", ModuleName )

    ! Now the data itself
    allocate(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
      qty%noLsts, qty%noSzas, qty%noDates), STAT=status)
    if ( status /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'qty%field')
  
  end subroutine SetupNewGriddedData

  ! ----------------------------------------NullifyGriddedData -----
  subroutine NullifyGriddedData ( G )
    ! Given a GriddedData, nullify all the pointers associated with it
    type ( GriddedData_T ) :: G

    ! Executable code
    nullify ( g%fileComments )
    nullify ( g%heights )
    nullify ( g%lats )
    nullify ( g%lons )
    nullify ( g%lsts )
    nullify ( g%szas )
    nullify ( g%dateStarts )
    nullify ( g%dateEnds )
    nullify ( g%field )
  end subroutine NullifyGriddedData

  ! ---------------------------------------- SliceGriddedData ------
  subroutine SliceGriddedData ( grid, slice, &
    & heights, lats, lons, lsts, szas, dates, missingValue )

    use MLSNumerics, only: HUNT, ESSENTIALLYEQUAL

    type ( GriddedData_T ), intent(inout) :: GRID ! Input grid
    real(rgr), dimension(:,:,:,:,:,:), intent(out) :: SLICE ! Result
    real(rgr), dimension(:), intent(in) :: HEIGHTS
    real(rgr), dimension(:), intent(in) :: LATS
    real(rgr), dimension(:), intent(in) :: LONS
    real(rgr), dimension(:), intent(in) :: LSTS
    real(rgr), dimension(:), intent(in) :: SZAS
    real(r8),  dimension(:), intent(in) :: DATES
    real(rgr), optional, intent(in) :: MISSINGVALUE

    ! Local parameters
    real(r8), parameter :: NOSECONDSINMEANYEAR = 9.55+60.0*(9.0+60.0*(6.0+24.0*365.0))
    ! 365 days, 6 hours, 9 minutes, 9.55 seconds

    ! Local variables
    ! First the dimensions of the slice
    integer :: NOHEIGHTS, NOLATS, NOLONS, NOLSTS, NOSZAS, NODATES
    ! Now the number of corners
    integer :: NOCORNERS
    ! Missing value to use
    real(rgr) :: MYMISSINGVALUE

    ! Now the indcies for each corner
    integer, dimension(size(heights),2) :: HEIGHTI
    integer, dimension(size(lats),2)    :: LATI
    integer, dimension(size(lons),2)    :: LONI
    integer, dimension(size(lsts),2)    :: LSTI
    integer, dimension(size(szas),2)    :: SZAI
    integer, dimension(size(dates),2)   :: DATEI

    ! Now the weights for each corner
    real(rgr), dimension(size(heights),2) :: HEIGHTW
    real(rgr), dimension(size(lats),2)    :: LATW
    real(rgr), dimension(size(lons),2)    :: LONW
    real(rgr), dimension(size(lsts),2)    :: LSTW
    real(rgr), dimension(size(szas),2)    :: SZAW
    real(rgr), dimension(size(dates),2)   :: DATEW

    ! Some extra special stuff for the dates
    real(r8), dimension(grid%noDates) :: MEANGRIDDATES
    real(r8), dimension(size(dates)) :: MODIFIEDINDATES
    integer :: YEARNUMBER

    ! Indices and loop counters
    integer :: C
    integer :: HEIGHT, LAT, LON, LST, SZA, DATE
    integer :: HEIGHTFAC, LATFAC, LONFAC, LSTFAC, SZAFAC, DATEFAC
    ! Indices of corner stuff
    integer, dimension(:), pointer :: HEIGHTC
    integer, dimension(:), pointer :: LATC
    integer, dimension(:), pointer :: LONC
    integer, dimension(:), pointer :: LSTC
    integer, dimension(:), pointer :: SZAC
    integer, dimension(:), pointer :: DATEC
    ! One result
    real(rgr) :: VAL

    ! Executable code
    myMissingValue = grid%missingValue
    if ( present ( missingValue ) ) myMissingValue = missingValue
    ! Get size of problem and check things out
    noHeights = size(heights)
    noLats = size(lats)
    noLons = size(lons)
    noLsts = size(lsts)
    noSzas = size(szas)
    noDates = size(dates)
    
    if ( any ( shape(slice) /= &
      &  (/ noHeights, noLats, noLons, noLsts, noSzas, noDates /) ) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Inappropriate size for slice' )

    ! Ensure the the input grid is 'wrapped'
    call WrapGriddedData ( grid )

    ! Now work out vertices and weights, and the number of corners
    noCorners = 1

    ! Height:
    if ( grid%noHeights == 1 ) then
      heightI = 1
      heightW(:,1) = 1.0_rgr
      heightW(:,2) = 0.0_rgr
      heightFac = 0
    else
      if ( grid%verticalCoordinate == v_is_pressure ) then
        call Hunt ( -log10(grid%heights), -log10(heights), heightI(:,1) )
        heightI(:,2) = heightI(:,1) + 1
        heightW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
          & ( -log10(grid%heights(heightI(:,2))) + log10(heights) ) / &
          & ( -log10(grid%heights(heightI(:,2))) + &
          &    log10(grid%heights(heightI(:,1))) ) ))
        heightW(:,2) = 1.0_rgr - heightW(:,1)
        
      else
        call Hunt ( grid%heights, heights, heightI(:,1) )
        heightI(:,2) = heightI(:,1) + 1
        heightW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
          & ( grid%heights(heightI(:,2)) - heights ) / &
          & ( grid%heights(heightI(:,2)) - grid%heights(heightI(:,1)) ) ))
        heightW(:,2) = 1.0_rgr - heightW(:,1)
      end if
      heightFac = noCorners
      noCorners = noCorners * 2
    end if

    ! Latitude
    if ( grid%noLats == 1 ) then
      latI = 1
      latW(:,1) = 1.0_rgr
      latW(:,2) = 0.0_rgr
      latFac = 0
    else
      call Hunt ( grid%lats, lats, latI(:,1) )
      latI(:,2) = latI(:,1) + 1
      latW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
        & ( grid%lats(latI(:,2)) - lats ) / &
        & ( grid%lats(latI(:,2)) - grid%lats(latI(:,1)) ) ))
      latW(:,2) = 1.0_rgr - latW(:,1)
      latFac = noCorners
      noCorners = noCorners * 2
    end if

    ! Longitude
    if ( grid%noLons == 1 ) then
      lonI = 1
      lonW(:,1) = 1.0_rgr
      lonW(:,2) = 0.0_rgr
      lonFac = 0
    else
      call Hunt ( grid%lons, lons, lonI(:,1) )
      lonI(:,2) = lonI(:,1) + 1
      lonW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
        & ( grid%lons(lonI(:,2)) - lons ) / &
        & ( grid%lons(lonI(:,2)) - grid%lons(lonI(:,1)) ) ))
      lonW(:,2) = 1.0_rgr - lonW(:,1)
      lonFac = noCorners
      noCorners = noCorners * 2
    end if
    
    ! Local solar time
    if ( grid%noLsts == 1 ) then
      lstI = 1
      lstW(:,1) = 1.0_rgr
      lstW(:,2) = 0.0_rgr
      lstFac = 0
    else
      call Hunt ( grid%lsts, lsts, lstI(:,1) )
      lstI(:,2) = lstI(:,1) + 1
      lstW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
        & ( grid%lsts(lstI(:,2)) - lsts ) / &
        & ( grid%lsts(lstI(:,2)) - grid%lsts(lstI(:,1)) ) ))
      lstW(:,2) = 1.0_rgr - lstW(:,1)
      lstFac = noCorners
      noCorners = noCorners * 2
    end if
    
    ! Solar zenith angle
    if ( grid%noSzas == 1 ) then
      szaI = 1
      szaW(:,1) = 1.0_rgr
      szaW(:,2) = 0.0_rgr
      szaFac = 0
    else
      call Hunt ( grid%szas, szas, szaI(:,1) )
      szaI(:,2) = szaI(:,1) + 1
      szaW(:,1) = max ( 0.0_rgr, min ( 1.0_rgr, &
        & ( grid%szas(szaI(:,2)) - szas ) / &
        & ( grid%szas(szaI(:,2)) - grid%szas(szaI(:,1)) ) ))
      szaW(:,2) = 1.0_rgr - szaW(:,1)
      szaFac = noCorners
      noCorners = noCorners * 2
    end if

    ! Dates.  Here we have to be a little clever
    meanGridDates = ( grid%dateStarts + grid%dateEnds ) / 2.0
    modifiedInDates = dates
    if ( grid%noYear ) then
      ! For 'noYear' grids, subtract the year of the dates
      ! Note, I'm not going to be terribly graceful about going over the
      ! year boundary, or for that matter leap years etc.
      yearNumber = int ( minval(modifiedInDates) / noSecondsInMeanYear )
      modifiedInDates = modifiedInDates - yearNumber * noSecondsInMeanYear
    endif
    if ( grid%noDates == 1 ) then
      dateI = 1
      dateW(:,1) = 1.0_rgr
      dateW(:,2) = 0.0_rgr
      dateFac = 0
    else
      call Hunt ( meanGridDates, modifiedInDates, dateI(:,1) )
      dateI(:,2) = dateI(:,1) + 1
      dateW(:,1) = max ( 0.0_r8, min ( 1.0_r8, &
        & ( meanGridDates(dateI(:,2)) - modifiedInDates ) / &
        & ( meanGridDates(dateI(:,2)) - meanGridDates(dateI(:,1)) ) ))
      dateW(:,2) = 1.0_rgr - dateW(:,1)
      dateFac = noCorners
      noCorners = noCorners * 2
    end if

    ! Now work out the corner information
    nullify ( heightC, latC, lonC, lstC, szaC, dateC )
    call Allocate_test ( heightC, noCorners, 'heightC', ModuleName )
    call Allocate_test ( latC, noCorners, 'latC', ModuleName )
    call Allocate_test ( lonC, noCorners, 'lonC', ModuleName )
    call Allocate_test ( lstC, noCorners, 'lstC', ModuleName )
    call Allocate_test ( szaC, noCorners, 'szaC', ModuleName )
    call Allocate_test ( dateC, noCorners, 'dateC', ModuleName )

    do c = 1, noCorners
      if ( heightFac == 0 ) then
        heightC(c) = 1
      else
        heightC(c) = merge ( 1, 2, mod((c-1)/heightFac,2) == 0 )
      end if

      if ( latFac == 0 ) then
        latC(c) = 1
      else
        latC(c) = merge ( 1, 2, mod((c-1)/latFac,2) == 0 )
      end if

      if ( lonFac == 0 ) then
        lonC(c) = 1
      else
        lonC(c) = merge ( 1, 2, mod((c-1)/lonFac,2) == 0 )
      end if

      if ( lstFac == 0 ) then
        lstC(c) = 1
      else
        lstC(c) = merge ( 1, 2, mod((c-1)/lstFac,2) == 0 )
      end if

      if ( szaFac == 0 ) then
        szaC(c) = 1
      else
        szaC(c) = merge ( 1, 2, mod((c-1)/szaFac,2) == 0 )
      end if

      if ( dateFac == 0 ) then
        dateC(c) = 1
      else
        dateC(c) = merge ( 1, 2, mod((c-1)/dateFac,2) == 0 )
      end if
    end do


    ! OK, now we do the work.
    ! We have two possibilities here.  If there are no bad datapoints
    ! then we can storm through, otherwise we have to be careful
    if ( any ( EssentiallyEqual ( grid%field, grid%missingValue ) ) ) then
      ! There are some missing values be very careful going throug the grid
      do date = 1, noDates
        do sza = 1, noSzas
          do lst = 1, noLsts
            do lon = 1, noLons
              do lat = 1, noLats
                do height = 1, noHeights
                  slice(height,lat,lon,lst,sza,date) = 0.0
                  cornerLoopMissing: do c = 1, noCorners
                    val = grid%field ( &
                      & heightI(height,heightC(c)), latI(lat,latC(c)), lonI(lon,lonC(c)), &
                      & lstI(lst,lstC(c)), szaI(sza,szaC(c)), dateI(date,dateC(c)) )
                    if ( EssentiallyEqual ( val, grid%missingValue ) ) then
                      slice(height,lat,lon,lst,sza,date) = myMissingValue
                      exit cornerLoopMissing
                    end if
                    slice(height,lat,lon,lst,sza,date) = &
                      & slice(height,lat,lon,lst,sza,date) + val * &
                      & heightW(height,heightC(c)) * latW(lat,latC(c)) * lonW(lon,lonC(c)) * &
                      & lstW(lst,lstC(c)) * szaW(sza,szaC(c)) * dateW(date,dateC(c))
                  end do cornerLoopMissing
                end do
              end do
            end do
          end do
        end do
      end do
    else
      ! Otherwise there are none, storm through
      do date = 1, noDates
        do sza = 1, noSzas
          do lst = 1, noLsts
            do lon = 1, noLons
              do lat = 1, noLats
                do height = 1, noHeights
                  slice(height,lat,lon,lst,sza,date) = 0.0
                  do c = 1, noCorners
                    val = grid%field ( &
                      & heightI(height,heightC(c)), latI(lat,latC(c)), lonI(lon,lonC(c)), &
                      & lstI(lst,lstC(c)), szaI(sza,szaC(c)), dateI(date,dateC(c)) )
                    slice(height,lat,lon,lst,sza,date) = &
                      & slice(height,lat,lon,lst,sza,date) + val * &
                      & heightW(height,heightC(c)) * latW(lat,latC(c)) * lonW(lon,lonC(c)) * &
                      & lstW(lst,lstC(c)) * szaW(sza,szaC(c)) * dateW(date,dateC(c))
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end if

    ! Tidy up
    call Deallocate_test ( heightC, 'heightC', ModuleName )
    call Deallocate_test ( latC, 'latC', ModuleName )
    call Deallocate_test ( lonC, 'lonC', ModuleName )
    call Deallocate_test ( lstC, 'lstC', ModuleName )
    call Deallocate_test ( szaC, 'szaC', ModuleName )
    call Deallocate_test ( dateC, 'dateC', ModuleName )
  end subroutine SliceGriddedData

  ! -------------------------------------------- WrapGriddedData ---
  subroutine WrapGriddedData ( GRID )
    ! Given a grid, possibly add extra points in longitude beyond +/-180
    ! and in solar time beyond 0..24 to aid in interpolations.
    ! Dummy arguments
    type ( GriddedData_T ), intent(inout) :: GRID
    ! Local variables
    real (rgr), dimension(:), pointer :: NEWLONS
    real (rgr), dimension(:), pointer :: NEWLSTS
    real (rgr), dimension(:,:,:,:,:,:), pointer :: NEWFIELD
    integer :: LOWERLON
    integer :: UPPERLON
    integer :: LOWERLST
    integer :: UPPERLST
    integer :: STATUS
    ! Executable code
    ! Don't bother with quantities that have no lon or lst variation.
    if ( grid%noLons <= 1 .and. grid%noLsts <= 1 ) return

    ! Check that this field is appropriate
    ! These integers are the number of points at or beyond our boundaries
    lowerLon = count ( grid%lons <= -180.0 )
    upperLon = count ( grid%lons >= 180.0 )
    lowerLst = count ( grid%lsts <= 0.0 )
    upperLst = count ( grid%lsts >= 24.0 )

    ! Check for unconventional grids.  We'll allow the case of 1 point at or
    ! beyond our boundaries, as that could just mean it's already been wrapped.
    ! More than that however is a problem.
    if ( lowerLon > 1 .or. upperLon > 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Gridded data disobeys EOS convention of -180..180 longitude' )
    endif
    if ( lowerLst > 1 .or. upperLst > 1 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Gridded data disobeys EOS convention of 0..24 local time' )
    endif

    ! From this point on the integers are the number of points
    ! we need to add (0 or 1)
    lowerLon = 0
    upperLon = 0
    lowerLst = 0
    upperLst = 0

    if ( grid%noLons > 1 ) then
      if ( grid%lons(1) > -180.0 ) lowerLon = 1
      if ( grid%lons(grid%noLons) < 180.0 ) upperLon = 1
    end if
    if ( grid%noLsts > 1 ) then
      if ( grid%lsts(1) > 0.0 ) lowerLst = 1
      if ( grid%lsts(grid%noLsts ) < 24.0 ) upperLst = 1
    end if

    ! If these are all zero then there's nothing to do so we might
    ! as well quit and save time
    if ( all ( (/ lowerLon, upperLon, lowerLst, upperLst /) == 0 ) ) return

    ! Create new arrays
    nullify ( newLons, newLSTs, newField )
    call Allocate_test ( newLons, grid%noLons + lowerLon + upperLon, &
      & 'newLons', ModuleName )
    call Allocate_test ( newLsts, grid%noLsts + lowerLst + upperLst, &
      & 'newLsts', ModuleName )
    allocate ( newField ( grid%noHeights, grid%noLats, &
      & grid%noLons + lowerLon + upperLon, &
      & grid%noLsts + lowerLst + upperLst, &
      & grid%noSzas, grid%noDates ), STAT=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Allocate//'newField' )

    ! Fill the 'central' part of the fields
    newLons ( 1+lowerLon : lowerLon+grid%noLons ) = grid%lons
    newLsts ( 1+lowerLst : lowerLst+grid%noLsts ) = grid%lsts
    newField ( :, :, &
      & 1+lowerLon : lowerLon+grid%noLons, &
      & 1+lowerLst : lowerLst+grid%noLsts, &
      & :, : ) = grid%field

    ! Wrap edges
    if ( lowerLon == 1 ) then
      newLons(1) = grid%lons ( grid%noLons ) - 360.0
      newField ( :, :, 1, :, :, : ) = grid%field ( :, :, grid%noLons, :, :, : )
    end if
    if ( upperLon == 1 ) then
      newLons ( grid%noLons + lowerLon + 1 ) = grid%lons ( 1 ) + 360.0
      newField ( :, :, grid%noLons + lowerLon + 1, :, :, : ) = &
        & grid%field ( :, :, 1, :, : ,: )
    end if

    if ( lowerLst == 1 ) then
      newLsts(1) = grid%lsts ( grid%noLsts ) - 24.0
      newField ( :, :, :, 1, :, : ) = grid%field ( :, :, :, grid%noLsts, :, : )
    end if
    if ( upperLst == 1 ) then
      newLsts ( grid%noLsts + lowerLst + 1 ) = grid%lsts ( 1 ) + 24.0
      newField ( :, :, :, grid%noLsts + lowerLst + 1, :, : ) = &
        & grid%field ( :, :, :, 1, : ,: )
    end if

    ! Tidy up
    call Deallocate_test ( grid%lons, 'grid%lons', ModuleName )
    call Deallocate_test ( grid%lsts, 'grid%lsts', ModuleName )
    deallocate ( grid%field, stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & MLSMSG_Deallocate//'grid%values' )

    ! Make grid use new values
    grid%noLons = grid%noLons + lowerLon + upperLon
    grid%noLsts = grid%noLsts + lowerLst + upperLst
    grid%lons => newLons
    grid%lsts => newLsts
    grid%field => newField
  end subroutine WrapGriddedData

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GriddedData

!
! $Log$
! Revision 2.33  2004/08/03 17:59:34  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.32  2004/05/19 18:54:44  vsnyder
! Remove declaration for unused variable
!
! Revision 2.31  2003/12/16 01:13:47  livesey
! Changed name in trace call
!
! Revision 2.30  2003/05/19 19:37:41  vsnyder
! Remove USE's for unreferenced names
!
! Revision 2.29  2003/05/09 01:54:58  livesey
! Got SliceGriddedData working
!
! Revision 2.28  2003/05/08 05:27:18  livesey
! Added the slicing
!
! Revision 2.27  2003/04/04 23:01:26  pwagner
! Short-curcuits dump for empty GriddedDatas
!
! Revision 2.26  2003/04/04 19:10:58  livesey
! Bug fix
!
! Revision 2.25  2003/04/04 00:09:32  livesey
! Added empty field, ConcatenateGriddedData, CopyGrid, WrapGriddedData,
! and appropriate changes to SetupNewGriddedData
!
! Revision 2.24  2003/03/01 00:22:30  pwagner
! Dump also prints missing value for that grid
!
! Revision 2.23  2003/02/28 02:26:42  livesey
! Added missingValue field
!
! Revision 2.22  2003/02/27 18:38:13  pwagner
! Removed some intent(out); Lahey takes perverse delight in resetting such to undefined
!
! Revision 2.21  2003/02/21 20:59:53  pwagner
! Tweaked dump settings
!
! Revision 2.20  2003/02/20 21:22:20  pwagner
! Changed default dump details to no for multidim arrays
!
! Revision 2.19  2003/02/19 19:13:28  pwagner
! new GriddedData_T with reduced precision
!
! Revision 2.18  2002/11/22 12:46:26  mjf
! Added nullify routine(s) to get round Sun's WS6 compiler not
! initialising derived type function results.
!
! Revision 2.17  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.16  2002/07/01 23:57:18  livesey
! Added the noYear field
!
! Revision 2.15  2002/06/27 00:09:50  vsnyder
! Don't deallocate qty%field if it's not associated.  Don't say allocate
! failed when the deallocate fails.
!
! Revision 2.14  2002/06/26 22:03:31  vsnyder
! Cosmetic changes
!
! Revision 2.13  2001/10/26 23:17:14  pwagner
! Provides a single dump module interface and details
!
! Revision 2.12  2001/09/10 23:36:56  livesey
! Expanded and tidied up.  Stuff that was in ncep_dao is now here
!
! Revision 2.11  2001/04/10 20:04:54  livesey
! Changed some NameLens to LineLens
!
! Revision 2.10  2001/03/30 00:24:40  pwagner
! Added sourceFileName
!
! Revision 2.9  2001/03/15 21:28:08  pwagner
! Now only the data type
!
! Revision 2.8  2001/03/15 00:37:13  pwagner
! Still not complete; missing gdrdfld
!
! Revision 2.7  2001/03/14 00:32:47  pwagner
! More changes--still wrong, though
!
! Revision 2.6  2001/03/10 00:33:16  pwagner
! Some corrections in ReadGriddedData
!
! Revision 2.5  2001/03/09 01:02:55  pwagner
! Fixed announce_error
!
! Revision 2.4  2001/03/08 01:08:35  pwagner
! Added announce_error
!
! Revision 2.3  2001/03/07 01:03:19  pwagner
! ReadGriddedData added
!
! Revision 2.2  2001/02/21 00:36:43  pwagner
! l3ascii_read_field now has eof as intent(out) arg
!
! Revision 2.1  2001/02/20 21:51:39  pwagner
! Functions absorbed from gridded_data_module
!
! Revision 2.0  2000/09/05 17:41:05  dcuddy
! Change revision to 2.0
!
! Revision 1.5  2000/06/20 22:19:22  lungu
! Changed DOUBBLE PRECISION to REAL (r8).
!
