! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GriddedData ! Contains the derived TYPE GriddedData_T

  use Allocate_Deallocate, only: Allocate_Test, Byte_Size, Bytes, &
    & Deallocate_Test, NoBytesAllocated, &
    & Test_Allocate, Test_Deallocate
  use HyperSlabs, only: EssentiallyEqual
  use Dates_Module, only: TAI2CCSDS
  use Dump_0, only: Dump
  use Dump_1, only: DumpDates
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: L_GeodAltitude, L_GPH, L_Eta, L_Pressure, &
    & L_Theta
  use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
  use MLSCommon, only: LineLen, NameLen, UndefinedValue
  ! r4 corresponds to sing. prec. :: same as stored in files
  ! (except for dao dimensions)
  use MLSKinds, only: RGR=>R4, R8
  use MLSMessageModule, only: MLSMSG_Error, &
    & MLSMSG_Warning, MLSMessageConfig, MLSMessage
  use MLSStringLists, only: SnipList, SwitchDetail
  use MLSStrings, only: LowerCase, ReadIntsFromChars
  use Output_m, only: OutputOptions, Blanks, Output, NewLine, &
    & ResumeOutput, SuspendOutput
  use Toggles, only: Gen, Levels, Switches, Toggle
  use Trace_m, only: Trace_Begin, Trace_End

  implicit none
  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters and data)
! GriddedData_T            The data type to hold climatology, gmao, etc.

!     (subroutines and functions)
! AddGriddedDataToDatabase What the name says
!                          
! ConcatenateGriddedData   Produce a a combination of two grids, 
!                           an earlier and a later
! ConvertFromEtaLevelGrids Convert two eta-level grids, one of them pressures,
!                           to a pressure-level Grid
! DestroyGriddedData       Deallocate all the pointer components
! DestroyGriddedDataDatabase
!                          Destroy all the elements of the database
! Diff                     Print differences between two grids
! DoGriddeddataMatch       Same shapes, are on same geolocations, etc.?
! DownsampleGriddeddata    Resample grid onto coarser grid derived by step sizes
! Dump                     Print details, values of a grid
! NullifyGriddedData       Nullify all its pointer components
! SetupNewGriddedData
!                          What the name says
! SliceGriddedData         Resample grid at supplied coords (Consider renaming it)
! WrapGriddedData          Add extra points in longitude beyond +/-180
!                           and in solar time beyond 0..24 to aid in interpolations
! === (end of toc) ===

! === (start of api) ===
! === (end of api) ===
  public :: GriddedData_T, AddgriddedDatatoDatabase, &
    & ConcatenategriddedData, Convertfrometalevelgrids, Copygrid, &
    & DestroygriddedData, DestroygriddedDataDatabase, &
    & Diff, DogriddedDatamatch, DownsamplegriddedData, Dump, &
    & NullifygriddedData, Rgr, SetupnewgriddedData, SlicegriddedData, &
    & WrapgriddedData
  
  logical, private, parameter :: MaydumpfieldValues = .true.

  interface ConcatenategriddedData
    module procedure ConcatenateGriddedData_2
    module procedure ConcatenateGriddedData_array
  end interface

  interface Diff
    module procedure DiffGriddedData
  end interface

  interface Dump
    module procedure DumpGriddedData
    module procedure DumpGriddedDatabase
  end interface

  ! These are 'enumerated types' consistent with hph's
  ! work in l3ascii_read_field
  public :: V_Is_Pressure, V_Is_Altitude, V_Is_Gph, V_Is_Theta, V_Is_Eta

  integer, parameter :: V_is_pressure = L_PRESSURE ! 1
  integer, parameter :: V_is_altitude = L_GEODALTITUDE ! v_is_pressure+1
  integer, parameter :: V_is_GPH      = L_GPH ! v_is_altitude+1
  integer, parameter :: V_is_theta    = L_THETA ! v_is_gph+1
  integer, parameter :: V_is_eta      = L_ETA ! V_is_theta+1
  
  ! If dumping gridded data, always give some details of any matching these
  character(len=*), parameter :: ALWAYSDUMPTHESE = ' ' ! 'dao,ncep,geos5'
  ! and for these automatically dumped ones, this level of detail for multi-dim
  integer, parameter :: AUTOMATICDETAILS = 0 ! 1 means dump, 0 means no
  
  ! Dump datebase memory footprint after every op
  logical, parameter, public :: DumpDBFootPrint = .true.

  ! This type reflects the format of the Level 3 ASCII files, though note that
  ! these files can store multiple quantities such as these.

  type GriddedData_T

    ! First the comment line(s) from the relevant input file
    logical :: EMPTY                    ! Set for example when file read failed
    logical :: deferReading             ! Set when to be read later
    character (len=LineLen), pointer, dimension(:) :: fileComments => NULL()

    ! Now the name, description and units information
    integer                 :: fileType = 0 ! From spec_indices, e.g. s_gridded
    character (len=LineLen) :: sourceFileName ! Input file name
    character (len=NameLen) :: quantityName ! From input file
    character (len=LineLen) :: description ! Quantity description
    character (len=LineLen) :: dimList
    character (len=LineLen) :: fieldNames
    character (len=NameLen) :: units ! Units for quantity

    ! Now define the various coordinate systems, first vertical
    integer :: verticalCoordinate ! An 'enumerated' type (e.g., V_is_...)
    integer :: noHeights         ! Number of surfaces
    real (rgr), pointer, dimension(:) :: heights  => NULL()
    ! Surfaces (e.g. pressures etc.) [noHeights]
    character (len=NameLen) :: heightsUnits ! Units for heights, e.g. 'Pa'

    ! Now the latitudinal coordinate
    logical :: equivalentLatitude       ! If set, coordinate is equivalent latitude
    logical :: noYear                   ! If set, field is for any year
    integer :: noLats                   ! Number of latitudes
    real (rgr), pointer, dimension(:) :: Lats => NULL() ! Latitudes [noLats]
    integer :: noLons                   ! Number of longitudes
    real (rgr), pointer, dimension(:) :: Lons => NULL() ! Longitudes [noLons]
    integer :: noLsts                      ! Number of local times
    real (rgr), pointer, dimension(:) :: Lsts => NULL() ! Local times [noLsts]
    integer :: noSzas                      ! Number of solar zenith angles
    real (rgr), pointer, dimension(:) :: Szas => NULL() ! Zenith angles [noSzas]
    integer :: noDates                     ! Number of dates in data
    real (r8), pointer, dimension(:) :: DateStarts => NULL()
    ! Starting dates in SDP toolkit format tai
    real (r8), pointer, dimension(:) :: DateEnds => NULL()
    ! Ending dates in SDP toolkit format tai

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
  subroutine ConcatenateGriddedData_2 ( A, B, X )
    ! This routine takes two grids A and B, B dated after A and tries
    ! to produce a third grid which is a combination of A and B
    type ( GriddedData_T ), intent(in) :: A
    type ( GriddedData_T ), intent(in) :: B
    type ( GriddedData_T ), intent(inout) :: X ! inout to let us deallocate it

    ! Local variables
    logical, parameter :: DEEBUG = .false.
    integer :: Me = -1                  ! String index for trace cacheing

    ! Executable code
    call trace_begin ( me, 'ConcatenateGriddedData_2', &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    ! First, check that the grids A and B are conformable.
    if ( (a%verticalCoordinate /= b%verticalCoordinate) .or. &
      & (a%equivalentLatitude .neqv. b%equivalentLatitude) .or. &
      & (a%noYear .neqv. b%noYear) ) then
        call output( 'Gridded Data a', advance='yes' )
        call dump( a )
        call output( 'Gridded Data b', advance='yes' )
        call dump( b )
        call outputNamedValue ( 'a%verticalCoordinate', a%verticalCoordinate )
        call outputNamedValue ( 'b%verticalCoordinate', b%verticalCoordinate )
        call outputNamedValue ( '=?', &
          & a%verticalCoordinate == b%verticalCoordinate )
        call outputNamedValue ( 'a%equivalentLatitude', a%equivalentLatitude )
        call outputNamedValue ( 'b%equivalentLatitude', b%equivalentLatitude )
        call outputNamedValue ( '=?', &
          & a%equivalentLatitude .eqv. b%equivalentLatitude )
        call outputNamedValue ( 'a%noYear', a%noYear )
        call outputNamedValue ( 'b%noYear', b%noYear )
        call outputNamedValue ( '=?', &
          & a%noYear .eqv. b%noYear )
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Grids for Concatenate are not compatible' )
    endif
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
    if ( maxval ( a%dateEnds ) > minval ( b%dateStarts ) ) then
      call Dump( a, details=0 )
      call Dump( b, details=0 )
      call output( 'Ending date of a: ', advance='no' )
      call output( maxval ( a%dateEnds ), format='(1pe16.9)', advance='yes' )
      call output( 'Starting date of b: ', advance='no' )
      call output( minval ( b%dateStarts ), format='(1pe16.9)', advance='yes' )
      call output( 'Difference: ', advance='no' )
      call output( maxval ( a%dateEnds ) - minval ( b%dateStarts ), format='(1pe16.9)', advance='yes' )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Grids for concatenation are not in correct time order' )
    endif

    ! OK, now we're ready
    call DestroyGriddedData ( X )
    call SetupNewGriddedData ( X, source=A, noDates= a%noDates + b%noDates )
    if ( DEEBUG ) then
    call output( 'source%verticalCoordinate: ', advance='no' )
    call output( a%verticalCoordinate, advance='no' )
    call blanks(3)
    call output( v_is_pressure, advance='yes' )
    call output( 'result%verticalCoordinate: ', advance='no' )
    call output( x%verticalCoordinate, advance='no' )
    call blanks(3)
    call output( v_is_pressure, advance='yes' )
    endif
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
    call trace_end ( 'ConcatenateGriddedData_2' , cond=toggle(gen) )
  end subroutine ConcatenateGriddedData_2

  subroutine ConcatenateGriddedData_array ( Database, indices, X )
    ! This routine takes an array of grids and a set of index values
    ! to produce a third grid which is a combination database elements
    ! at the index values
    type ( GriddedData_T ), dimension(:), intent(in) :: database
    integer, dimension(:), intent(in) :: indices ! index values
    type ( GriddedData_T ), intent(inout) :: X ! inout to let us deallocate it
    ! Local variables
    logical, parameter :: DEEBUG = .false.
    integer :: i
    integer, dimension(size(indices)) :: i1
    integer, dimension(size(indices)) :: i2
    integer :: index1
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: noDates
    ! Executable code
    call trace_begin ( me, 'ConcatenateGriddedData_array', &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    call DestroyGriddedData ( X )
    if ( size(indices) < 1 .or. size(database) < 1 ) go to 9
    index1 = indices(1)
    if ( index1 < 1 .or. index1 > size(Database) ) go to 9
    if ( size(indices) < 2 ) then
      call CopyGrid ( X, Database(index1) )
      go to 9
    endif
    noDates = 0
    do i=1, size(indices)
      i1(i) = noDates + 1
      noDates = noDates + Database(indices(i))%noDates
      i2(i) = noDates
    enddo
    call SetupNewGriddedData ( X, source=Database(index1), noDates= noDates )
    ! Copy over the unchanged position data
    x%heights = Database(index1)%heights
    x%lats = Database(index1)%lats
    x%lons = Database(index1)%lons
    x%lsts = Database(index1)%lsts
    x%szas = Database(index1)%szas

    x%dateStarts(1:i2(1)) = Database(index1)%dateStarts
    x%dateEnds  (1:i2(1)) = Database(index1)%dateEnds

    x%field ( :, :, :, :, :, 1:i2(1) ) = Database(index1)%field
    do i=2, size(indices)
      index1 = indices(i)
      if ( index1 < 1 .or. index1 > size(database) ) cycle
      ! First, check that the grids A and B are conformable.
      if ( Database(index1)%verticalCoordinate /= x%verticalCoordinate .or. &
        & Database(index1)%equivalentLatitude .neqv. x%equivalentLatitude &
        & .or. &
        & Database(index1)%noYear .neqv. x%noYear ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate are not compatible' )
      if ( .not. all ( EssentiallyEqual &
        & ( Database(index1)%heights, x%heights ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate do not share heights' )
      if ( .not. all ( EssentiallyEqual ( Database(index1)%lats, x%lats ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate do not share lats' )
      if ( .not. all ( EssentiallyEqual ( Database(index1)%lons, x%lons ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate do not share lons' )
      if ( .not. all ( EssentiallyEqual ( Database(index1)%lsts, x%lsts ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate do not share lsts' )
      if ( .not. all ( EssentiallyEqual ( Database(index1)%szas, x%szas ) ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate do not share szas' )
      if ( .not. EssentiallyEqual &
        & ( Database(index1)%missingValue, x%missingValue ) ) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for Concatenate have different missing values' )

      ! Check that the dates are going to work out
      if ( maxval ( x%dateEnds(:i1(i)-1) ) > &
        & minval ( Database(index1)%dateStarts ) ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Grids for concatenation are not in correct time order' )
      endif

      x%dateStarts(i1(i):i2(i)) = Database(index1)%dateStarts
      x%dateEnds(i1(i):i2(i)) = Database(index1)%dateEnds
      x%field ( :, :, :, :, :, i1(i):i2(i) ) = Database(index1)%field
    enddo
  9 call trace_end ( 'ConcatenateGriddedData_array' , cond=toggle(gen) )

  end subroutine ConcatenateGriddedData_array

  ! ---------------------------------------- ConvertFromEtaLevelGrids ------------------
  subroutine ConvertFromEtaLevelGrids ( TGrid, PGrid, NGrid, OutGrid, &
    & VGrid, ByLog )
    ! Converts two eta-level grids, one of them pressures, 
    ! to a pressure-level Grid
    use MLSNumerics, only: Interpolatevalues, Uselookuptable
    use Dump_0, only: Dump
    use VGridsDatabase, only: VGrid_T
    type ( GriddedData_T ), intent(in)  :: TGrid  ! E.g., T on eta surfaces
    type ( GriddedData_T ), intent(in)  :: PGrid  ! Pressures on eta surfaces
    type ( GriddedData_T ), intent(in)  :: NGrid  ! What surfaces to use
    type ( GriddedData_T ), intent(out) :: OutGrid ! T on pressure level
    type ( VGrid_T ), pointer           :: VGrid  ! What surfaces to use
    logical, optional, intent(in)       :: ByLog  ! Interpolate using zeta

    ! Internal variables
    integer :: iDate, iSza, iLst, iLon, iLat, iHeight
    integer :: Me = -1                  ! String index for trace cacheing
    real(rgr), dimension(TGrid%noHeights) :: pEta
    logical, parameter                    :: DEEBUG = .false.
    logical                               :: GotVgrid ! Are heights supplied by a VGrid?
    logical                               :: why
    logical                               :: ZetaSurfaces
    real (rgr), pointer, dimension(:)     :: heights  => NULL()

    ! Executable code
    call trace_begin ( me, 'ConvertFromEtaLevelGrids', &
      & cond=toggle(gen) .and. levels(gen) > 1 )
    ZetaSurfaces = .false.
    if ( present(ByLog) ) ZetaSurfaces = ByLog
    GotVgrid = associated(VGrid)

    ! GriddedData must match
    if ( .not. DoGriddeddataMatch( PGrid, TGrid ) ) then
      call DumpGriddedData( PGrid )
      call DumpGriddedData( TGrid )
      why = DoGriddeddataMatch( PGrid, TGrid, SayWhyNot=.true. )
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Gridded T,P data must match' )
    endif
    call DestroyGriddedData ( OutGrid )
    ! Copy the information over
    if ( GotVgrid ) then
      call SetupNewGriddedData ( OutGrid, source=TGrid, noHeights=VGrid%noSurfs )
      OutGrid%verticalCoordinate = VGrid%VerticalCoordinate
      OutGrid%heights(1:OutGrid%noHeights) = VGrid%Surfs(1:OutGrid%noHeights,1)
    else
      call SetupNewGriddedData ( OutGrid, source=TGrid, noHeights=NGrid%noHeights )
      OutGrid%verticalCoordinate = NGrid%VerticalCoordinate
      if ( size(OutGrid%heights) /= size(NGrid%heights) ) then
        call outputNamedValue('num heights(outGrid)', size(OutGrid%heights), advance='yes' )
        call outputNamedValue('num heights(NGrid)', size(NGrid%heights), advance='yes' )
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Grid shapes do not conform' )
        go to 9
      endif
      OutGrid%heights(1:OutGrid%noHeights) = NGrid%heights(1:OutGrid%noHeights)
    endif
    if ( PGrid%empty .or. TGrid%empty) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'Temperatures or Pressures grid was empty' )
      go to 9
    endif
    OutGrid%lats = TGrid%lats
    OutGrid%lons = TGrid%lons
    OutGrid%lsts = TGrid%lsts
    OutGrid%szas = TGrid%szas
    OutGrid%dateStarts = TGrid%dateStarts
    OutGrid%dateEnds = TGrid%dateEnds
    heights => OutGrid%heights
    ! Now we'll interpolate to the OutGrid surfaces
    do idate=1, TGrid%noDates
      do iSza=1, TGrid%noSzas
        do iLst=1, TGrid%noLsts
          do iLon=1, TGrid%noLons
            do iLat=1, TGrid%noLats
              ! But do we need to convert pressures from one system to another?
              ! Because we're in a hurry, we won't use ConvertVGrid
              ! but just do it by hand
              ! What we'll do is assume
              ! The grid we need OutGrid is in 'hPa'
              ! The eta-level pressures are either already in hPa or else in pa
              if ( lowercase(PGrid%units) == 'hpa' ) then
                pEta = PGrid%field( :, iLat, iLon, iLst, iSza, iDate )
              elseif ( lowercase(PGrid%units) == 'pa' ) then
                pEta = PGrid%field( :, iLat, iLon, iLst, iSza, iDate ) / 100
              else
                call MLSMessage(MLSMSG_Error, ModuleName, &
                  & 'unrecognized quantity for eta-level pressures '// trim(PGrid%units) )
              endif
              if ( ZetaSurfaces ) then
                call InterpolateValues( &
                & log10(pEta), TGrid%field( :, iLat, iLon, iLst, iSza, iDate ), &
                & log10(heights), OutGrid%field( :, iLat, iLon, iLst, iSza, iDate ), &
                & 'L', 'B', TGrid%missingValue )
              else
                call InterpolateValues( &
                & pEta, TGrid%field( :, iLat, iLon, iLst, iSza, iDate ), &
                & heights, OutGrid%field( :, iLat, iLon, iLst, iSza, iDate ), &
                & 'L', 'B', TGrid%missingValue )
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
    if ( DEEBUG ) then
      call outputNamedValue( 'PGrid%units', PGrid%units )
      call dump( TGrid%field( :, 1, 1, 1, 1, 1 ), 'T' )
      call dump( PGrid%field( :, 1, 1, 1, 1, 1 ), 'P' )
      call dump( heights, 'heights' )
      call dump( OutGrid%field( :, 1, 1, 1, 1, 1 ), 'T out' )
      pEta = PGrid%field( :, 1, 1, 1, 1, 1 ) / 100
      call dump( pEta, 'heights(eta)' )
      do iHeight = 1, OutGrid%noHeights
        call outputNamedValue( 'height', heights(iHeight) )
        OutGrid%field( iHeight, 1, 1, 1, 1, 1 ) = &
        & UseLookUpTable ( heights(iHeight), TGrid%field( :, 1, 1, 1, 1, 1 ), &
        & xtable=pEta, missingValue=TGrid%missingValue, options='ip' )
      enddo
      call dump( OutGrid%field( :, 1, 1, 1, 1, 1 ), 'T out (2)' )
      stop
    endif
  9 call trace_end ( 'ConvertFromEtaLevelGrids' , cond=toggle(gen) )
  end subroutine ConvertFromEtaLevelGrids

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
    ! This subroutine deallocates all the pointer components
    use Trace_M, only: Trace_begin, Trace_end

    ! Dummy argument
    type (GriddedData_T), intent(inout) :: Qty

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S ! Size in bytes of a deallocated field
    integer :: STATUS
    logical :: verbose
    integer :: Me = -1           ! String index for trace

    ! Executable code
    call trace_begin ( me, "DestroyGriddedData", &
      & cond=toggle(gen) .and. levels(gen) > 1  )

    verbose = switchDetail(switches, 'grid') > -1
    call Deallocate_test ( qty%heights, "griddedQty%heights", ModuleName )
    call Deallocate_test ( qty%lats, "griddedQty%lats", ModuleName )
    call Deallocate_test ( qty%lons, "griddedQty%lons", ModuleName )
    call Deallocate_test ( qty%lsts, "griddedQty%lsts", ModuleName )
    call Deallocate_test ( qty%szas, "griddedQty%szas", ModuleName )

    ! Now the temporal coordinates
    call Deallocate_test ( qty%dateStarts, "griddedQty%dateStarts", ModuleName )
    call Deallocate_test ( qty%dateEnds, "griddedQty%dateEnds", ModuleName )

    ! Now the data itself
    if ( associated(qty%field) ) then
      s = byte_size(qty%field)
      if ( verbose ) then
        call outputnamedValue ( 'Grid size to be destroyed', shape(qty%field) )
        call outputnamedValue ( 'Its size in MB', s/1.e6 )
      endif
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(qty%field(1,1,1,1,1,1)), addr)
      deallocate(qty%field, STAT=status)
      call test_deallocate ( status, moduleName, "griddedQty%field", s, address=addr )
    else if ( verbose ) then
      call output( 'This grid not allocated', advance='true' )
    end if
    qty%empty = .true.
    call trace_end ( "DestroyGriddedData", &
      & cond=toggle(gen) .and. levels(gen) > 1  )

  end subroutine DestroyGriddedData

  ! ----------------------------------  DestroyGriddedDataDatabase -----
  subroutine DestroyGriddedDataDatabase ( Database )
  ! This subroutine destroys a quantity template database

    use Allocate_Deallocate, only: Test_Deallocate

    use Trace_M, only: Trace_begin, Trace_end

    ! Dummy argument
    type (GriddedData_T), dimension(:), pointer :: Database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: Me = -1           ! String index for trace
    integer :: qtyIndex, s, status

    call trace_begin ( me, "DestroyGriddedDataDatabase", &
      & cond=toggle(gen) .and. levels(gen) > 1 )

    if (associated(database)) then
      do qtyIndex=1,size(database)
        call DestroyGriddedData(database(qtyIndex))
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate(database, stat=status)
      call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if
    call trace_end ( "DestroyGriddedDataDatabase", cond=toggle(gen) )

  end subroutine DestroyGriddedDataDatabase

   ! --------------------------------------------  DiffGriddedData  -----
  subroutine DiffGriddedData ( GriddedData1, GriddedData2, options )
    use Diff_1, only: Diff

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao
    type (GriddedData_T) :: GriddedData1, GriddedData2
    character(len=*), intent(in), optional :: options

    ! Local Variables
    integer :: date
    integer :: details
    character(len=4) :: dateChar

    ! Executable code
    if ( GriddedData1%empty .or. GriddedData2%empty ) then
      call output('At least one Gridded quantity was empty (perhaps the file name' &
        & // ' was wrong)', advance='yes')
      return
    endif
    details = 0
    if ( present(options) ) then
      call readIntsFromChars( options, details, ignore='-*' )
      call outputNamedValue('options', &
      & options, advance='yes')
    endif
    call output('Gridded quantity (1) name ' // &
      & GriddedData1%quantityName, advance='yes')
    call output('Gridded quantity (2) name ' // &
      & GriddedData2%quantityName, advance='yes')
    call outputNamedValue('details', &
      & details, advance='yes')
   ! May dump a 3-d slice of 6-d array
    if ( details == 1 ) then
      call diff ( &
        & GriddedData1%field(:,1,1,1,1,1), &
        & GriddedData1%quantityName, &
        & GriddedData2%field(:,1,1,1,1,1), &
        & GriddedData2%quantityName, &
        & options=options )
    elseif ( details == 2 ) then
      call diff ( &
        & GriddedData1%field(:,:,1,1,1,1), &
        & GriddedData1%quantityName, &
        & GriddedData2%field(:,:,1,1,1,1), &
        & GriddedData2%quantityName, &
        & options=options )
    elseif ( GriddedData1%noDates == 1 .and. GriddedData1%noSzas == 1 &
      & .and. GriddedData1%noLsts == 1 ) then
      call diff ( &
        & GriddedData1%field(:,:,:,1,1,1), &
        & GriddedData1%quantityName, &
        & GriddedData2%field(:,:,:,1,1,1), &
        & GriddedData2%quantityName, &
        & options=options )
    elseif ( GriddedData1%noSzas == 1 .and. GriddedData1%noLsts == 1 ) then
   ! May dump a 4-d array as noDates instances of 3-d arrays
      do date=1, GriddedData1%noDates
        write(dateChar, '(i4)') date
        call diff ( &
          & GriddedData1%field(:,:,:,1,1,date), &
          & GriddedData1%quantityName, &
          & GriddedData2%field(:,:,:,1,1,date), &
          & GriddedData2%quantityName, &
          & options=options )
      enddo
    else
      call output ( ' *(Sorry, diff_6d_... not yet coded)* ' , &
        & advance='yes')
    endif

  end subroutine DiffGriddedData

  ! ----------------------------------------  DoGriddeddataMatch  -----
  function DoGriddeddataMatch ( a, b, sayWhyNot ) result( match )
    ! Test whether two griddeddata have same shapes, are on same geolocations
    ! and so on
    ! They may have different name, description, and unit information
    ! And the field values may differ
    type(GriddedData_T), intent(in) :: a
    type(GriddedData_T), intent(in) :: b
    logical, optional, intent(in)   :: sayWhyNot
    logical                         :: match
    ! local variables
    logical                         :: verbose
    ! Executable
    verbose = .false.
    if ( present(sayWhyNot) ) verbose = sayWhyNot
    match = .true.
    if ( verbose ) then
      if ( a%empty .or. b%empty ) call explain( 'One empty, the other not' )
      if ( any( shape(a%field) /= shape(b%field) ) ) call explain( 'shapes dont match' )
      if ( a%verticalCoordinate /= b%verticalcoordinate ) call explain( 'vertical coords diff' )
      if ( a%noHeights /= b%noHeights ) call explain( 'noHeights diff' )
      if ( any( a%Heights /= b%Heights ) ) call explain( 'Heights diff' )
      if ( lowercase(a%heightsUnits) /= lowercase(b%heightsUnits) ) call explain( 'height units idff' )
      ! if ( lowercase(a%Units) /= lowercase(b%Units) ) call explain( 'One empty, the other not' )
      if ( a%noLats /= b%noLats ) call explain( 'noLats diff' )
      if ( a%noLons /= b%noLons ) call explain( 'nolens diff' )
      if ( a%noLsts /= b%noLsts ) call explain( 'noLsts diff' )
      if ( a%noSzas /= b%noSzas ) call explain( 'noSzas diff' )
      if ( a%noDates /= b%noDates ) call explain( 'nodates difft' )
      if ( a%missingValue /= b%missingValue ) call explain( 'missing values diff' )
      if ( any( a%Lats /= b%Lats ) ) call explain( 'lats diff' )
      if ( any( a%Lons /= b%Lons ) ) call explain( 'lons diff' )
      if ( any( a%Lsts /= b%Lsts ) ) call explain( 'lsts diff' )
      if ( any( a%Szas /= b%Szas ) ) call explain( 'szas diff' )
      if ( any( a%DateStarts /= b%DateStarts ) ) call explain( 'dateStarts diff' )
      if ( any( a%DateEnds /= b%DateEnds ) ) call explain( 'DateEnds diff' )
    endif
    if ( a%empty .and. b%empty ) return
    match = .false.
    if ( a%empty .or. b%empty ) return
    if ( any( shape(a%field) /= shape(b%field) ) ) return
    if ( a%verticalCoordinate /= b%verticalcoordinate ) return
    if ( a%noHeights /= b%noHeights ) return
    if ( any( a%Heights /= b%Heights ) ) return
    if ( lowercase(a%heightsUnits) /= lowercase(b%heightsUnits) ) return
    ! if ( lowercase(a%Units) /= lowercase(b%Units) ) return
    if ( a%noLats /= b%noLats ) return
    if ( a%noLons /= b%noLons ) return
    if ( a%noLsts /= b%noLsts ) return
    if ( a%noSzas /= b%noSzas ) return
    if ( a%noDates /= b%noDates ) return
    if ( a%missingValue /= b%missingValue ) return
    if ( any( a%Lats /= b%Lats ) ) return
    if ( any( a%Lons /= b%Lons ) ) return
    if ( any( a%Lsts /= b%Lsts ) ) return
    if ( any( a%Szas /= b%Szas ) ) return
    if ( any( a%DateStarts /= b%DateStarts ) ) return
    if ( any( a%DateEnds /= b%DateEnds ) ) return
    match = .true.
  contains
    subroutine explain( message )
      character(len=*), intent(in) :: message
      call output( trim(message), advance='yes' )
    end subroutine explain
  end function DoGriddeddataMatch

  ! ---------------------------------------- DownSampleGriddedData ------
  ! Given a gridded data type on a (too-)fine resolution, create a new
  ! one on a coarser mesh defined by steps through the finer mesh
  ! E.g., a latsStep of 2 means choosing every 2nd latitude
  
  ! See also SliceGriddeddata
  ! Why did we choose 2, 4, 6, .. instead of 1, 3, 5 ?
  ! Also what if number of steps does not divide mesh points
  ! into whole number?
  
  ! Suggested improvement:
  ! we already have various implementations of hyperslab in the software
  ! so just reuse one of them;
  ! i.e., pass in start, stride, and edge arrays
  subroutine DownSampleGriddedData ( Grid, Newgrid, &
    & heightsStep, latsStep, lonsStep, lstsStep, szasStep, datesStep, &
    & firstheights, firstlats, firstlons, firstlsts, firstszas, firstdates )

    type ( GriddedData_T ), intent(inout) :: Grid ! Input grid
    type ( GriddedData_T ), intent(out)   :: Newgrid ! Result
    integer, intent(in)                   :: heightsStep
    integer, intent(in)                   :: latsStep
    integer, intent(in)                   :: lonsStep
    integer, intent(in)                   :: lstsStep
    integer, intent(in)                   :: szasStep
    integer, intent(in)                   :: datesStep
    integer, intent(in)                   :: firstheights
    integer, intent(in)                   :: firstlats
    integer, intent(in)                   :: firstlons
    integer, intent(in)                   :: firstlsts
    integer, intent(in)                   :: firstszas
    integer, intent(in)                   :: firstdates
    ! Internal variables
    integer                   :: noheights
    integer                   :: nolats
    integer                   :: nolons
    integer                   :: nolsts
    integer                   :: noszas
    integer                   :: nodates
    ! Executable
    noheights = (grid%noheights - 1)/heightsStep + 1
    nolats    = (grid%nolats - 1   )/latsStep    + 1
    nolons    = (grid%nolons - 1   )/lonsStep    + 1
    nolsts    = (grid%nolsts - 1   )/lstsStep    + 1
    noszas    = (grid%noszas - 1   )/szasStep    + 1
    nodates   = (grid%nodates - 1  )/datesStep   + 1
    ! call outputNamedValue( 'grid%noheights ', grid%noheights )
    ! call outputNamedValue( 'heightsStep ', heightsStep )
    ! call outputNamedValue( 'noheights ', noheights )
    ! call outputNamedValue( 'nolats    ', nolats    )
    ! call outputNamedValue( 'nolons    ', nolons    )
    ! call outputNamedValue( 'nolsts    ', nolsts    )
    ! call outputNamedValue( 'noszas    ', noszas    )
    ! call outputNamedValue( 'nodates   ', nodates   )
    ! call dump( grid, details=1 )
    ! Set up newgrid
    call SetupNewGriddedData ( newgrid, grid, NoHeights, NoLats, &
      & NoLons, NoLsts, NoSzas, NoDates )
    ! Now fill its geo-fields
    newgrid%heights = grid%heights( heightsStep::heightsStep )
    newgrid%lats = grid%lats( firstlats::latsStep )
    newgrid%lons = grid%lons( firstlons::lonsStep )
    newgrid%lsts = grid%lsts( firstlsts::lstsStep )
    newgrid%szas = grid%szas( firstszas::szasStep )
    newgrid%dateStarts = grid%dateStarts( datesStep::datesStep )
    newgrid%dateEnds = grid%dateEnds( datesStep::datesStep )
    ! Now the multidimension gridded field values
    newgrid%field = grid%field( &
      & firstheights::heightsStep, &
      & firstlats   ::latsStep, &
      & firstlons   ::lonsStep, &
      & firstlsts   ::lstsStep, &
      & firstszas   ::szasStep, &
      & firstdates  ::datesStep &
      & )
    ! call dump( newgrid, details=1 )
  end subroutine DownSampleGriddedData

  ! ----------------------------------------  DumpGriddedDatabase  -----
  subroutine DumpGriddedDatabase ( GriddedData, Details, options )

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao

    type (GriddedData_T), dimension(:), pointer :: GriddedData 

    integer, intent(in), optional               :: DETAILS ! If -4, just show footprint
    character(len=*), intent(in), optional      :: OPTIONS
    
    ! Local Variables
    integer            :: i
    integer            :: MYDETAILS
    real               :: total

    myDetails = 0
    if ( present(details) ) myDetails = details
    if ( .not. associated(GriddedData)) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, 'Gridded database still null')

    if ( myDetails < -3 ) call SuspendOutput
    call output ( '============ Gridded Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    call output ( 'database: a priori grids: SIZE = ' )
    call output ( size(GriddedData), advance='yes' )
    if ( size(GriddedData) < 1 ) return
    total = 0.
    do i = 1, size(GriddedData)

      call outputNamedValue ( 'item number', i, advance='no' )
      if ( myDetails < -2 ) then
        call outputNamedValue('  Gridded quantity name ', &
          & trim(GriddedData(i)%quantityName), advance='no')
        call OutputNamedValue ( '  empty', GriddedData(i)%empty )
      else
        call DumpGriddedData( GriddedData(i), Details, Options )
      endif
      if ( .not. associated(GriddedData(i)%field) ) cycle
      total = total + product(shape(GriddedData(i)%field))
    end do ! i
    if ( myDetails < -3 ) call ResumeOutput
    call outputNamedValue( 'Database Total Memory Footprint', &
      & bytes(GriddedData(1)%field) * total )
  end subroutine DumpGriddedDatabase

  ! --------------------------------------------  DumpGriddedData  -----
  subroutine DumpGriddedData ( GriddedData, Details, options )
    use Dump_0, only: Dump
    use Ieee_Arithmetic, only: Ieee_Is_Finite, Ieee_Is_Nan
    use Intrinsic, only: Lit_Indices, Spec_Indices
    use String_Table, only: Display_String

    ! Imitating what dump_pointing_grid_database does, but for gridded data
    ! which may come from climatology, ncep, dao
    type (GriddedData_T) :: GriddedData 
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 0
    !                                        ! > 2 Climatology text file format
    character(len=*), intent(in), optional      :: OPTIONS
    ! Options, if present, says which of the 6 possible dimensions
    ! to include in the multi-dim dump of the gridded field
    ! E.g., if options = 'lat, height, time' then just dump the 3-d array
    !   field (:,:,1,1,1,:)
    ! Local Variables
    integer :: caseID
    character(len=8) :: ccsds
    integer :: date
    character(len=4) :: dateChar
    integer :: FIELDVALUESDETAILS
    integer :: i
    logical :: lookLikeClimatologyTxtfile
    integer :: MYDETAILS
    character(len=16) :: myOptions
    integer :: numElmnts
    integer :: numSurfs
    character(len=32) :: oldInfo

    ! Executable code
    myDetails = 0
    if ( present(details) ) myDetails = details
    fieldvaluesdetails = myDetails
    if ( index(ALWAYSDUMPTHESE, trim(GriddedData%description)) > 0 ) then
      myDetails = max( myDetails, 1 )
      fieldvaluesdetails = max(fieldvaluesdetails, AUTOMATICDETAILS)
    endif
    lookLikeClimatologyTxtfile = ( myDetails > 2 ) .and. .not. present(options)
    myOptions = ' '
    if ( present(options) ) myOptions = lowerCase(options)
    if ( myOptions == '-' ) myOptions = ' '
    if ( GriddedData%empty ) then
      call output('This Gridded quantity was empty', advance='yes')
      return
    endif
    if ( lookLikeClimatologyTxtfile ) then
      call newLine
      ! Format dump to look like Climatology text file
      call output('Field ' // trim(GriddedData%quantityName), advance='no')
      call blanks (3)
      call output(trim(GriddedData%units), advance='yes')
      ! Vertical coords
      outputOptions%arrayElmntSeparator = ','
      outputOptions%nArrayElmntsPerLine = 6
      call output('; Define vertical coordinates', advance='yes')
      call output('Pressure ', advance='no')
      call output('Explicit ', advance='no')
      call output('( ', advance='yes')
      call output(GriddedData%heights, advance='no')
      call output(') ', advance='yes')
      ! Horizontal coords
      call output('; Define horizontal coordinates etc.', advance='yes')
      call output('Latitude ', advance='no')
      call output('Explicit ', advance='no')
      call output('( ', advance='yes')
      call output(GriddedData%Lats, advance='no')
      call output(') ', advance='yes')
      call output('Longitude ', advance='no')
      call output('Explicit ', advance='no')
      call output('( ', advance='yes')
      call output(GriddedData%Lons, advance='no')
      call output(') ', advance='yes')
      if ( associated(GriddedData%Lsts) ) then
        call output('lst ', advance='no')
        call output('Explicit ', advance='no')
        call output('( ', advance='yes')
        call output(GriddedData%Lsts, advance='no')
        call output(') ', advance='yes')
      endif
      if ( associated(GriddedData%Szas) ) then
        call output('sza ', advance='no')
        call output('Explicit ', advance='no')
        call output('( ', advance='yes')
        call output(GriddedData%Szas, advance='no')
        call output(') ', advance='yes')
      endif
      outputOptions%arrayElmntSeparator = ' '
      outputOptions%nArrayElmntsPerLine = 7
      outputOptions%nBlanksBtwnElmnts   = 1
      outputOptions%sdFormatDefault     = '(1pe10.3)'
      numElmnts = GriddedData%noHeights * &
        & max(1, size(GriddedData%Lats)) * &
        & max(1, size(GriddedData%Lons)) * &
        & max(1, size(GriddedData%Lsts)) * &
        & max(1, size(GriddedData%Szas))
      do i=1, size(GriddedData%DateStarts)
        ccsds = tai2ccsds( int(GriddedData%DateStarts(i) / 86400 + 0.5) )
        !if ( i > 1 ) then
        !   date = (GriddedData%DateStarts(i) - GriddedData%DateStarts(1)) / 86400 + 1.5
        !else
        !  date = 1
        !endif
        read( ccsds(6:8), *) date
        call output('Date Single ', advance='no')
        call output(-date, advance='yes')
        call output( reshape( GriddedData%field(:,:,:,:,:,i), (/numElmnts/) ), &
          & advance='no')
        call newLine
      enddo
      outputOptions%sdFormatDefault     = '*' ! Restore default format
      outputOptions%nBlanksBtwnElmnts   = 3
      return
    endif
    call output('Gridded quantity name ' // GriddedData%quantityName, advance='yes')
    if ( associated(GriddedData%field)) then
      call OutputNamedValue ( 'shape(field)', &
        &  shape(GriddedData%field) )
      call OutputNamedValue ( 'Grid Memory Footprint', &
        &  byte_size(GriddedData%field) )
    endif
    if ( myDetails < -1 ) return
    oldInfo = MLSMessageConfig%Info
    MLSMessageConfig%Info = GriddedData%quantityName
    call output( 'fileType: ')
    if ( GriddedData%fileType < 1 ) then
      call output( '(unknown)', advance='yes' )
    else
      call display_string( spec_indices(GriddedData%fileType), advance='yes' )
    endif
    call output(' description ' // GriddedData%description, advance='yes')
    call output(' units ' // trim(GriddedData%units), advance='yes')
    call output(' missing value ', advance='no')
    call output(GriddedData%missingValue, advance='yes')

    call output ( ' ************ Geometry ********** ' ,advance='yes')

    call output ( ' Vertical coordinate = ' )
    if ( GriddedData%verticalCoordinate > 0 ) then
      call display_string ( lit_indices(GriddedData%verticalCoordinate), advance='yes' )
    else
      call output ( GriddedData%verticalCoordinate, advance='yes' )
    end if
    call output(' heights units ' // trim(GriddedData%heightsUnits), advance='yes')
    numSurfs = GriddedData%noHeights
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

    call dumpDates( GriddedData%dateStarts, 'starting dates' )
    call dumpDates( GriddedData%dateEnds, 'ending dates' )
    
    call output( 'No year in date? ' )
    call output( GriddedData%NoYear, advance='yes' )

    call output( 'Seasonal cyclic symmetry? ' )
    call output( GriddedData%NoYear .and. Griddeddata%noDates > 1, &
      & advance='yes' )

    if ( .not. all(ieee_is_finite(GriddedData%field)) ) &
      & call output ( '*** Gridded Data contains non-finite values', &
        & advance='yes' )
    if ( any(ieee_is_nan(GriddedData%field)) ) &
      & call output ( '*** Gridded Data contains nans', advance='yes' )
    if ( .not. associated(GriddedData%field) ) then
      call output ( ' Gridded Data field values empty', advance='yes' )
      return
    endif
    call outputNamedValue( 'max(values)', maxval(GriddedData%field) )
    call outputNamedValue( 'min(values)', minval(GriddedData%field) )
    call outputNamedValue( 'shape(values)', shape(GriddedData%field) )
    if ( len_trim(myOptions) > 0 ) then
      if ( GriddedData%noLons < 2 .and. index( myOptions, 'lon') > 0 ) &
        & myOptions = SnipList( myOptions, 'lon' )
      if ( GriddedData%noheights < 2 .and. index( myOptions, 'height') > 0 ) &
        & myOptions = SnipList( myOptions, 'height' )
      if ( GriddedData%nolats < 2 .and. index( myOptions, 'lat') > 0 ) &
        & myOptions = SnipList( myOptions, 'lat' )
      if ( GriddedData%noDates < 2 .and. index( myOptions, 'date') > 0 ) &
        & myOptions = SnipList( myOptions, 'date' )
      if ( GriddedData%nolsts < 2 .and. index( myOptions, 'lst') > 0 ) &
        & myOptions = SnipList( myOptions, 'lst' )
      if ( GriddedData%noszas < 2 .and. index( myOptions, 'sza') > 0 ) &
        & myOptions = SnipList( myOptions, 'sza' )
      caseID = 0
      if ( index(myOptions, 'height') > 0 ) caseID = caseID + 1
      if ( index(myOptions, 'lat') > 0 ) caseID = caseID + 2
      if ( index(myOptions, 'lon') > 0 ) caseID = caseID + 4
      if ( index(myOptions, 'lst') > 0 ) caseID = caseID + 8
      if ( index(myOptions, 'sza') > 0 ) caseID = caseID + 16
      if ( index(myOptions, 'date') > 0 ) caseID = caseID + 32
      call outputNamedValue( 'values dump options', trim(myOptions) )
      call outputNamedValue( 'caseID', caseID )
      select case (caseID)
      case ( 0 )
        call output( 'Default options for dumping values', advance='yes' )
        call outputNamedValue( 'values(1,1,..,1)', GriddedData%field(1,1,1,1,1,1) )
      case ( 1 )
        call dump ( GriddedData%field(:,1,1,1,1,1), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case ( 2 )
        call dump ( GriddedData%field(1,:,1,1,1,1), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case ( 3 )
        call dump ( GriddedData%field(:,:,1,1,1,1), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case ( 4 )
        call dump ( GriddedData%field(1,1,:,1,1,1), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case ( 7 )
        call dump ( GriddedData%field(:,:,:,1,1,1), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case ( 35 )
        ! call dump ( shape(GriddedData%field), 'shape(values)' )
        ! call dump ( shape(GriddedData%field(:,:,1,1,1,:)), 'shape(dumped values)' )
        call dump ( GriddedData%field(:,:,1,1,1,:), &
          & '    gridded field values', FillValue=GriddedData%MissingValue )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'This case not coded' )
      end select
    elseif ( MAYDUMPFIELDVALUES .and. fieldvaluesdetails > 0 ) then
      call output ( ' ************ tabulated field values ********** ' ,advance='yes')
      if ( fieldvaluesdetails < 2 ) then
      ! Do 1d slices through data
        if ( size(GriddedData%field, 1 ) > 1 ) &
          & call dump ( GriddedData%field(:,1,1,1,1,1), &
          & '    gridded field values sliced through heights =', &
          & FillValue=GriddedData%MissingValue )
        if ( size(GriddedData%field, 2 ) > 1 ) &
          & call dump ( GriddedData%field(numSurfs,:,1,1,1,1), &
          & '    gridded field values sliced through lats =', &
          & FillValue=GriddedData%MissingValue )
        if ( size(GriddedData%field, 3 ) > 1 ) &
          & call dump ( GriddedData%field(numSurfs,1,:,1,1,1), &
          & '    gridded field values sliced through lons =', &
          & FillValue=GriddedData%MissingValue )
        if ( size(GriddedData%field, 4 ) > 1 ) &
          & call dump ( GriddedData%field(1,1,1,:,1,1), &
          & '    gridded field values sliced through sol times =', &
          & FillValue=GriddedData%MissingValue )
        if ( size(GriddedData%field, 5 ) > 1 ) &
          & call dump ( GriddedData%field(1,1,1,1,:,1), &
          & '    gridded field values sliced through sol zeniths =', &
          & FillValue=GriddedData%MissingValue )
        if ( size(GriddedData%field, 6 ) > 1 ) &
          & call dump ( GriddedData%field(1,1,1,1,1,:), &
          & '    gridded field values sliced through dates =', &
          & FillValue=GriddedData%MissingValue )
     ! May dump a 3-d slice of 6-d array
      elseif ( GriddedData%noDates == 1 .and. GriddedData%noSzas == 1 &
        & .and. GriddedData%noLsts == 1 ) then
        call dump ( GriddedData%field(:,:,:,1,1,1), &
          & '    gridded field values =', &
          & FillValue=GriddedData%MissingValue )
      elseif ( GriddedData%noDates == 1 .and. GriddedData%noSzas == 1 ) then
     ! May dump a 3-d slices through this 4-d array
        call dump ( GriddedData%field(:,:,:,1,1,1), &
          & '    gridded field values (1st solar time) =' , &
          & FillValue=GriddedData%MissingValue )
        call dump ( GriddedData%field(:,:,1,:,1,1), &
          & '    gridded field values (1st longitude) =' , &
          & FillValue=GriddedData%MissingValue )
        call dump ( GriddedData%field(:,1,:,:,1,1), &
          & '    gridded field values (1st latitude) =' )
        call dump ( GriddedData%field(1,:,:,:,1,1), &
          & '    gridded field values (1st height) =' , &
          & FillValue=GriddedData%MissingValue )
      elseif ( GriddedData%noSzas == 1 .and. GriddedData%noLsts == 1 ) then
     ! May dump a 4-d array as noDates instances of 3-d arrays
        do date=1, GriddedData%noDates
          write(dateChar, '(i4)') date
          call dump ( GriddedData%field(:,:,:,1,1,date), &
            & '    gridded field values (' // &
            & trim(dateChar) // 'th solar time) =' , &
            & FillValue=GriddedData%MissingValue )
        enddo
      else
        ! No dump for 6-dimensional double arrays yet, anyway
        !     call dump ( GriddedData%field, &
        !      & '    gridded field values =' )
        call output ( ' *(Sorry, dump_6d_... not yet coded)* ' , &
          & advance='yes')
      endif
    endif
    MLSMessageConfig%Info = oldInfo
  end subroutine DumpGriddedData

  ! ----------------------------------------  SetupNewGriddedData  -----
  subroutine SetupNewGriddedData ( Qty, Source, NoHeights, NoLats, &
    & NoLons, NoLsts, NoSzas, NoDates, missingValue, verticalCoordinate, &
    & units, heightsUnits, Empty )
  ! This first routine sets up a new grid according to the user
  ! input.  This may be based partly on an already-defined source
  ! or created from scratch.
    ! Dummy arguments
    type (GriddedData_T) :: QTY ! Result
    type (GriddedData_T), optional, intent(in) :: SOURCE ! Template
    integer, optional, intent(in) :: NOHEIGHTS, NOLATS, NOLONS, NOLSTS, NOSZAS, NODATES
    integer, optional, intent(in) :: verticalCoordinate
    character(len=*), optional, intent(in) :: units
    character(len=*), optional, intent(in) :: heightsUnits
    logical, optional, intent(in) :: EMPTY
    real(rgr), optional, intent(in) :: missingValue
    ! Local parameters
    ! real(rgr), parameter :: DefaultMissingValue = undefinedValue ! -999.99
    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: status                   ! Status from allocates etc.
    logical :: myEmpty                  ! Copy of empty possibly
    logical :: verbose

    ! Executable code
    verbose = switchDetail(switches, 'grid') > -1
    qty%description    = '(unknown)'
    qty%sourceFileName = '(unknown)'
    qty%quantityName   = '(unknown)'
    qty%units          = '(unknown)'
    qty%heightsUnits   = '(unknown)'
    qty%dimList        = '(unknown)'
    qty%deferReading   = .false.

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
        qty%units = source%units
        qty%quantityName = source%quantityName
        qty%heightsUnits = source%heightsUnits
        qty%dimList = source%dimList
      else ! We have no template, setup a very bare quantity
        qty%noHeights = 1
        qty%noLats = 1
        qty%noLons = 1
        qty%noLsts = 1
        qty%noSzas = 1
        qty%noDates = 1
        qty%equivalentLatitude = .false.
        qty%noYear = .false.
        qty%missingValue = undefinedValue ! defaultMissingValue
        qty%verticalCoordinate = -1
        qty%units = '(unknown)'
        qty%heightsUnits = '(unknown)'
      endif
      
      ! Now, see if the user asked for modifications to this
      if (present(noHeights)) qty%noHeights = noHeights
      if (present(noLats)) qty%noLats = noLats
      if (present(noLons)) qty%noLons = noLons
      if (present(noLsts)) qty%noLsts = noLsts
      if (present(noSzas)) qty%noSzas = noSzas
      if (present(noDates)) qty%noDates = noDates
      if ( present ( missingValue ) ) qty%missingValue = missingValue
      if ( present ( verticalCoordinate ) ) qty%verticalCoordinate = verticalCoordinate
      if ( present ( units ) ) qty%units = units
      if ( present ( heightsunits ) ) qty%heightsunits = heightsunits
    end if

    ! First the vertical/horizontal coordinates
    if ( verbose .and. qty%noHeights > 0 ) then
      call outputNamedValue( 'qty%noHeights', qty%noHeights )
      call outputNamedValue( 'qty%noLats', qty%noLats )
      call outputNamedValue( 'qty%noLons', qty%noLons )
      call outputNamedValue( 'qty%noLsts', qty%noLsts )
      call outputNamedValue( 'qty%noSzas', qty%noSzas )
      call outputNamedValue( 'qty%noDates', qty%noDates )
      call outputNamedValue( 'Bytes before allocating vertical/horizontal coordinates', &
        & NoBytesAllocated )
    endif
    call Allocate_test ( qty%heights, qty%noHeights, "griddedQty%heights", ModuleName )
    call Allocate_test ( qty%lats, qty%noLats, "griddedQty%lats", ModuleName )
    call Allocate_test ( qty%lons, qty%noLons, "griddedQty%lons", ModuleName )
    call Allocate_test ( qty%lsts, qty%noLsts, "griddedQty%lsts", ModuleName )
    call Allocate_test ( qty%szas, qty%noSzas, "griddedQty%szas", ModuleName )

    ! Now the temporal coordinates
    call Allocate_test ( qty%dateStarts, qty%noDates, "griddedQty%dateStarts", ModuleName )
    call Allocate_test ( qty%dateEnds, qty%noDates, "griddedQty%dateEnds", ModuleName )

    ! Now the data itself
    if ( verbose ) call outputNamedValue( 'Bytes before allocating grid%field', NoBytesAllocated )
    allocate(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
      qty%noLsts, qty%noSzas, qty%noDates), STAT=status)
    if ( status == 0 ) then
      if ( size(qty%field) > 0 ) addr = transfer(c_loc(qty%field(1,1,1,1,1,1)), addr)
    end if
    call test_allocate ( status, moduleName, "griddedQty%field", (/1,1,1,1,1,1/), &
      & (/ qty%noHeights, qty%noLats, qty%noLons, qty%noLsts, qty%noSzas, &
      &    qty%noDates /), bytes(qty%field), address=addr )
    if ( verbose ) then
      call outputNamedValue( 'Bytes after allocating grid%field', NoBytesAllocated )
      call outputNamedValue( 'Size of field (elements)', product(shape(qty%field)) )
      call outputNamedValue( 'Byte Size of field', byte_size(qty%field) )
    endif
  
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
  ! 1st suggestion:
  ! Rename this procedure something more apt
  ! E.g., 'Regrid' or 'Resample'
  ! Here would be a good place to describe the purpose and method of this
  ! procedure (hint, hint)
  ! I (paw) did not write this but I can take a stab at it
  ! return an array with values of the griddeddata's field at mesh
  ! points defined by the supplied heights, lats, etc.
  ! Where and if the new mesh points do not exactly line up with the old
  ! produce a weighted average from the corners (vertices) of the smallest
  ! 6-dimensional hypercube of old mesh points surrounding each new mesh point
  ! Confusing? Maybe a picture would help ..
  ! If we had not 6 dimensions but only 2, let lowercase letters be the old mesh
  ! points and upper case letters be the new mesh points. Then we would look
  ! for points a, b, c, and d such that for point P
  !
  !   a            b
  !
  !       P
  !
  !   c            d
  ! and determine the value at P by weighting the values at a, b, c, and d
  ! according to a mysterious formula I have not taken the trouble to
  ! understand, but which, among other wrinkles, treats heights logarithmically
  ! and makes an effort to avoid allowing undefined values to exert a
  ! corrupting and baleful influence if they should be found at a, b, c, or d
  
  ! See also DownSampleGriddeddata
  subroutine SliceGriddedData ( grid, slice, &
    & heights, lats, lons, lsts, szas, dates, missingValue )

    use MLSNumerics, only: Hunt

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
    logical, parameter :: DEEBUG = .false.
    real(r8), parameter :: NOSECONDSINMEANYEAR = 9.55+60.0*(9.0+60.0*(6.0+24.0*365.0))
    ! 365 days, 6 hours, 9 minutes, 9.55 seconds

    ! Local variables
    integer :: Me = -1                  ! String index for trace cacheing
    ! Dimensions of the slice
    integer :: NOHEIGHTS, NOLATS, NOLONS, NOLSTS, NOSZAS, NODATES
    ! Number of corners
    integer :: NOCORNERS
    ! Missing value to use
    real(r8) :: MYMISSINGVALUE

    ! Now the indices for each corner
    integer, dimension(size(heights),2) :: HEIGHTI
    integer, dimension(size(lats),2)    :: LATI
    integer, dimension(size(lons),2)    :: LONI
    integer, dimension(size(lsts),2)    :: LSTI
    integer, dimension(size(szas),2)    :: SZAI
    integer, dimension(size(dates),2)   :: DATEI

    ! Now the weights for each corner
    real(r8), dimension(size(heights),2) :: HEIGHTW
    real(r8), dimension(size(lats),2)    :: LATW
    real(r8), dimension(size(lons),2)    :: LONW
    real(r8), dimension(size(lsts),2)    :: LSTW
    real(r8), dimension(size(szas),2)    :: SZAW
    real(r8), dimension(size(dates),2)   :: DATEW

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
    real(r8) :: VAL
    real(r8), dimension(size(heights)) :: subSlice

    ! Executable code
    call trace_begin ( me, 'SliceGriddedData' , cond=.false. )
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
      heightW(:,1) = 1.0_r8
      heightW(:,2) = 0.0_r8
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
      latW(:,1) = 1.0_r8
      latW(:,2) = 0.0_r8
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
      lonW(:,1) = 1.0_r8
      lonW(:,2) = 0.0_r8
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
      lstW(:,1) = 1.0_r8
      lstW(:,2) = 0.0_r8
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
      szaW(:,1) = 1.0_r8
      szaW(:,2) = 0.0_r8
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
      if ( minval(modifiedInDates) < 0. ) yearNumber = yearNumber - 1
      modifiedInDates = modifiedInDates - yearNumber * noSecondsInMeanYear
      if ( deebug ) then
        call outputNamedValue( 'yearNumber', yearNumber )
        call outputNamedValue( 'modifiedInDates', modifiedInDates(1) )
        call outputNamedValue( 'meanGridDates', meanGridDates(1) )
      endif
    endif
    if ( grid%noDates == 1 ) then
      dateI = 1
      dateW(:,1) = 1.0_r8
      dateW(:,2) = 0.0_r8
      dateFac = 0
    else
      call Hunt ( meanGridDates, modifiedInDates, dateI(:,1) )
      dateI(:,2) = dateI(:,1) + 1
      dateW(:,1) = max ( 0.0_r8, min ( 1.0_r8, &
        & ( meanGridDates(dateI(:,2)) - modifiedInDates ) / &
        & ( meanGridDates(dateI(:,2)) - meanGridDates(dateI(:,1)) ) ))
      dateW(:,2) = 1.0_r8 - dateW(:,1)
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

    if ( DEEBUG ) then
    call output('slicing .. ' , advance='yes' )
    call outputNamedValue ( 'num corners', noCorners )
    call outputNamedValue ( 'heightFac', heightFac )
    call dump( heightC, 'heightC' )
    call dump( heightI, 'heightI' )
    call dump( heightW, 'heightW' )
    call outputNamedValue ( 'latFac', latFac )
    call dump( latC, 'latC' )
    call dump( latI, 'latI' )
    call dump( latW, 'latW' )
    call outputNamedValue ( 'lonFac', lonFac )
    call dump( lonC, 'lonC' )
    call dump( lonI, 'lonI' )
    call dump( lonW, 'lonW' )
    call outputNamedValue ( 'dateFac', dateFac )
    call dump( dateC, 'dateC' )
    call dump( dateI, 'dateI' )
    call dump( dateW, 'dateW' )
    call outputNamedValue ( 'szaFac', szaFac )
    call dump( szaC, 'szaC' )
    call dump( szaI, 'szaI' )
    call dump( szaW, 'szaW' )
    call outputNamedValue ( 'lstFac', lstFac )
    call dump( lstC, 'lstC' )
    call dump( lstI, 'lstI' )
    call dump( lstW, 'lstW' )
    call outputNamedValue ( 'missing value', grid%missingValue ) 
    call outputNamedValue ( 'any missing value?', any ( EssentiallyEqual ( grid%field, grid%missingValue ) ) )
    endif
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
                subSlice = 0._r8
                do height = 1, noHeights
                  ! slice(height,lat,lon,lst,sza,date) = 0.0
                  cornerLoopMissing: do c = 1, noCorners
                    val = grid%field ( &
                      & heightI(height,heightC(c)), latI(lat,latC(c)), lonI(lon,lonC(c)), &
                      & lstI(lst,lstC(c)), szaI(sza,szaC(c)), dateI(date,dateC(c)) )
                    if ( EssentiallyEqual ( val, real(grid%missingValue, r8) ) ) then
                      ! slice(height,lat,lon,lst,sza,date) = myMissingValue
                      subSlice(height) = myMissingValue
                      exit cornerLoopMissing
                    end if
                    ! slice(height,lat,lon,lst,sza,date) = &
                      ! & slice(height,lat,lon,lst,sza,date) + val * &
                    subSlice(height) = subSlice(height) + val * &
                      & heightW(height,heightC(c)) * latW(lat,latC(c)) * lonW(lon,lonC(c)) * &
                      & lstW(lst,lstC(c)) * szaW(sza,szaC(c)) * dateW(date,dateC(c))
                  end do cornerLoopMissing
                end do
                slice(:,lat,lon,lst,sza,date) = subSlice
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
                subSlice = 0._r8
                do height = 1, noHeights
                  ! slice(height,lat,lon,lst,sza,date) = 0.0
                  do c = 1, noCorners
                    val = grid%field ( &
                      & heightI(height,heightC(c)), latI(lat,latC(c)), lonI(lon,lonC(c)), &
                      & lstI(lst,lstC(c)), szaI(sza,szaC(c)), dateI(date,dateC(c)) )
                    ! slice(height,lat,lon,lst,sza,date) = &
                    !  & slice(height,lat,lon,lst,sza,date) + val * &
                    subSlice(height) = subSlice(height) + val * &
                      & heightW(height,heightC(c)) * latW(lat,latC(c)) * lonW(lon,lonC(c)) * &
                      & lstW(lst,lstC(c)) * szaW(sza,szaC(c)) * dateW(date,dateC(c))
                    if ( height == 8 .and. lat == 1 .and. lon == 1 .and. &
                      & lst == 1 .and. sza == 1 .and. date == 1 .and. DEEBUG ) &
                      & call output ( (/ val, heightW(height,heightC(c)), &
                        & latW(lat,latC(c)), lonW(lon,lonC(c)), &
                        & lstW(lst,lstC(c)), szaW(sza,szaC(c)), dateW(date,dateC(c)) /), &
                        & advance='yes' )
                  end do
                end do
                slice(:,lat,lon,lst,sza,date) = subSlice
                if ( lat == 1 .and. lon == 1 .and. &
                  & lst == 1 .and. sza == 1 .and. date == 1 .and. DEEBUG ) &
                  & call output ( subSlice, advance='yes' )
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
    if ( DEEBUG ) call dump ( slice(:,1:10,1,:,1,1), &
        & '    sliced field values (1st longitude) =' , &
        & FillValue=Grid%MissingValue )
    call trace_end ( 'SliceGriddedData' , cond=.false. )

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
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: LOWERLON
    integer :: LOWERLST
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: S ! Size in bytes of a deallocated field
    integer :: STATUS
    integer :: UPPERLON
    integer :: UPPERLST

    ! Executable code
    call trace_begin ( me, 'WrapGriddedData' , cond=.false. )
    ! Don't bother with quantities that have no lon or lst variation.
    if ( grid%noLons <= 1 .and. grid%noLsts <= 1 ) go to 9

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
    end if

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
    if ( all ( (/ lowerLon, upperLon, lowerLst, upperLst /) == 0 ) ) go to 9

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
    addr = 0
    if ( status==0 ) then
      if ( size(newField) > 0 ) addr = transfer(c_loc(newField(1,1,1,1,1,1)), addr)
    end if
    call test_allocate ( status, moduleName, 'newField', (/1,1,1,1,1,1/), &
      & (/ grid%noHeights, grid%noLats, &
      &    grid%noLons + lowerLon + upperLon, &
      &    grid%noLsts + lowerLst + upperLst, &
      &    grid%noSzas, grid%noDates /), bytes(newField), address=addr )

    ! Fill the 'central' part of the fields
    newLons ( 1+lowerLon : lowerLon+grid%noLons ) = grid%lons
    newLsts ( 1+lowerLst : lowerLst+grid%noLsts ) = grid%lsts

    ! Wrap edges
    if ( lowerLon == 1 ) then
      newLons(1) = grid%lons ( grid%noLons ) - 360.0
    end if
    if ( upperLon == 1 ) then
      newLons ( grid%noLons + lowerLon + 1 ) = grid%lons ( 1 ) + 360.0
    end if

    if ( lowerLst == 1 ) then
      newLsts(1) = grid%lsts ( grid%noLsts ) - 24.0
    end if
    if ( upperLst == 1 ) then
      newLsts ( grid%noLsts + lowerLst + 1 ) = grid%lsts ( 1 ) + 24.0
    end if

    call wraparrays( grid%field, newfield, lowerLon, upperLon, lowerLst, upperLst, &
    & grid%noLons, grid%noLsts )
    ! Tidy up
    call Deallocate_test ( grid%lons, 'grid%lons', ModuleName )
    call Deallocate_test ( grid%lsts, 'grid%lsts', ModuleName )
    s = byte_size(grid%field)
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(grid%field(1,1,1,1,1,1)), addr)
    deallocate ( grid%field, stat=status )
    call test_deallocate ( status, moduleName, 'grid%field', s, address=addr )

    ! Make grid use new values
    grid%noLons = grid%noLons + lowerLon + upperLon
    grid%noLsts = grid%noLsts + lowerLst + upperLst
    grid%lons => newLons
    grid%lsts => newLsts
    grid%field => newField
  9 call trace_end ( 'WrapGriddedData' , cond=.false. )

  end subroutine WrapGriddedData

  ! -------- Private ---------------
  subroutine wrapArrays ( old, new, lowerLon, upperLon, lowerLst, upperLst, &
    & oldNoLons, oldNoLsts )
  ! Hide our array shenanigans to
  ! (1) prevent compilers like NAG from using a temp array and running out
  !     of memory
  ! (2) prevent cleverer compilers like Intel from checking whether they
  !     need a temp array
  ! Args
    real (rgr), dimension(:,:,:,:,:,:) :: old
    real (rgr), dimension(:,:,:,:,:,:) :: new
    integer, intent(in) :: lowerLon, upperLon, lowerLst, upperLst, &
    & oldNoLons, oldNoLsts
    new ( :, :, &
      & 1+lowerLon : lowerLon+oldNoLons, &
      & 1+lowerLst : lowerLst+oldNoLsts, &
      & :, : ) = old
    ! Wrap edges
    if ( lowerLon == 1 ) then
      new ( :, :, 1, :, :, : ) = old ( :, :, oldNoLons, :, :, : )
    end if
    if ( upperLon == 1 ) then
      new ( :, :, oldNoLons + lowerLon + 1, :, :, : ) = &
        & old ( :, :, 1, :, : ,: )
    end if

    if ( lowerLst == 1 ) then
      new ( :, :, :, 1, :, : ) = old ( :, :, :, oldNoLsts, :, : )
    end if
    if ( upperLst == 1 ) then
      new ( :, :, :, oldNoLsts + lowerLst + 1, :, : ) = &
        & old ( :, :, :, 1, : ,: )
    end if
  end subroutine wrapArrays

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module GriddedData

!
! $Log$
! Revision 2.90  2021/07/08 23:27:24  pwagner
! Correct unbalanced trace for DestroyGriddedData
!
! Revision 2.89  2020/05/29 21:54:49  pwagner
! More detailed tracing of memory usage
!
! Revision 2.88  2020/04/30 23:20:59  pwagner
! If Details lt -2 will dump only name and whether empty
!
! Revision 2.87  2020/04/27 21:34:51  pwagner
! trace_.. added to more carefully track memory usage
!
! Revision 2.86  2019/10/03 17:31:24  pwagner
! Convert from eta levels may now take a vGrid field
!
! Revision 2.85  2019/09/23 20:38:32  pwagner
! Conversion from Eta surfaces may optionally be logarithmic
!
! Revision 2.84  2019/09/05 17:50:13  pwagner
! SetUp with sources quantityName, dimList, too
!
! Revision 2.83  2017/11/03 20:01:16  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.82  2017/08/23 16:45:25  pwagner
! Fixed bug when Dates are before 1993
!
! Revision 2.81  2017/03/17 00:10:49  pwagner
! Optionally sayWhyNot if not DoGriddeddataMatch
!
! Revision 2.80  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.79  2015/03/28 01:01:03  vsnyder
! Stuff to trace allocate/deallocate addresses -- mostly commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.78  2015/01/21 19:29:15  pwagner
! More debugging if PGrid, TGrid dont match
!
! Revision 2.77  2014/09/04 23:38:45  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.76  2014/07/18 22:00:24  pwagner
! Aimed for consistency in names passed to allocate_test
!
! Revision 2.75  2014/05/29 18:28:26  pwagner
! Correctly compute memory usage
!
! Revision 2.74  2014/01/09 00:25:06  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.73  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.72  2013/08/30 03:56:01  vsnyder
! Revise use of trace_begin and trace_end
!
! Revision 2.71  2013/06/12 02:17:27  vsnyder
! UBYTES and BYTE_SIZE from Allocate_Deallocate
!
! Revision 2.70  2012/10/11 21:00:28  pwagner
! Print quantityName instead of moduleName during Dump
!
! Revision 2.69  2012/03/06 19:32:39  pwagner
! Remove more unused local variables
!
! Revision 2.68  2012/03/06 19:12:32  pwagner
! Say whether dumped grid has seasonal cyclic symmetry
!
! Revision 2.67  2012/03/01 20:01:40  pwagner
! When dumped, say whether grid is for no year
!
! Revision 2.66  2011/08/30 22:22:29  pwagner
! Destroying a Griddeddata now sets empty to TRUE so we can dump it later w/o error
!
! Revision 2.65  2011/06/29 21:37:17  pwagner
! Corrected bug in recalculating downsized dim sizes
!
! Revision 2.64  2011/06/16 23:00:09  pwagner
! Added DownSampleGriddedData
!
! Revision 2.63  2011/05/05 15:21:34  pwagner
! Added fileType field to datatype
!
! Revision 2.62  2011/04/27 17:34:25  pwagner
! dateStarts and -Ends now dumped as dates, not doubles
!
! Revision 2.61  2011/04/20 00:25:18  pwagner
! Dump now responds to options arg
!
! Revision 2.60  2010/02/09 16:24:12  pwagner
! Hide large array section assignments in subroutine to prevent NAG from creating temps
!
! Revision 2.59  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.58  2010/02/03 23:09:23  vsnyder
! Use Test_Deallocate
!
! Revision 2.57  2009/10/30 23:04:50  pwagner
! Added diff of gridded data
!
! Revision 2.56  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.55  2009/06/16 17:25:13  pwagner
! Can dump Gridded data to look like Climatology file
!
! Revision 2.54  2009/01/12 18:45:46  pwagner
! Added print statement to not_used_here
!
! Revision 2.53  2008/09/16 21:05:10  pwagner
! Improved dumps
!
! Revision 2.52  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.51  2008/01/07 21:36:33  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.50  2007/10/04 01:49:14  vsnyder
! Don't overflow the call stack
!
! Revision 2.49  2007/08/17 00:27:23  pwagner
! push more procedures onto MLSCallStack
!
! Revision 2.48  2007/08/13 17:35:24  pwagner
! Push SliceGriddedData onto new MLSSCallStack
!
! Revision 2.47  2007/07/27 00:23:57  vsnyder
! Use Test_Allocate after allocations
!
! Revision 2.46  2007/06/07 20:26:16  pwagner
! Prevents some crashes, writing non-ascii chars to stdout
!
! Revision 2.45  2007/03/23 00:10:49  pwagner
! Valiant attempts to bring two Lahey version results closer
!
! Revision 2.44  2007/01/30 21:57:44  pwagner
! Avoids a bug with empty grids in ConvertFromEtaLevelGrids
!
! Revision 2.43  2007/01/11 20:31:53  vsnyder
! Spiff up the dump
!
! Revision 2.42  2006/11/01 20:27:05  pwagner
! Hasty fix to units bug in eta-level conversion
!
! Revision 2.41  2006/06/14 23:58:46  pwagner
! Concatenate can take array of grids
!
! Revision 2.40  2006/06/13 22:10:19  pwagner
! changed interface to ConvertFromEtaLevelGrids
!
! Revision 2.39  2006/05/12 21:24:13  pwagner
! Added extra debugging statements
!
! Revision 2.38  2006/05/09 00:13:32  pwagner
! Added DoGriddeddataMatch; heightsUnits as component of GriddedData_T
!
! Revision 2.37  2006/05/02 19:00:42  pwagner
! Added ConvertFromEtaLevelGrids; preliminary, untested
!
! Revision 2.36  2006/02/10 21:21:20  pwagner
! Added V_is_eta verticalCoordinate
!
! Revision 2.35  2006/01/19 00:24:17  pwagner
! Small improvements to dump routine
!
! Revision 2.34  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
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
