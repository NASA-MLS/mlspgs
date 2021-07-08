! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module HGridsDatabase                   ! Horizontal grid information

  use Intrinsic, only: L_GeodAngle
  use MLSKinds, only: R8
  use Generate_QTM_m, only: QTM_Tree_t

  implicit none
  private

  public :: HGrid_T, HGrids_T, HGridGeolocations_T, HGridGeolocations
  public :: AddHGridtodatabase, CopyHGrid, CreateEmptyHGrid, &
    & DestroyHGridContents, DestroyHGridDatabase, Dump, FindClosestMatch, &
    & L1BGeoLocation, L1BSubsample, NullifyHGrid, &
    & TrimHGrid

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information.
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  type HGrid_T
    integer :: Name                 = 0 ! String index of name.            
    integer :: MasterCoordinate     = l_GeodAngle ! Its lit index;
    integer :: NoProfs                  ! Number of profiles in this grid  
    integer :: NoProfsLowerOverlap  = 0 ! Number of profiles in the lower overlap
    integer :: NoProfsUpperOverlap  = 0 ! Number of profiles in the upper overlap
    integer :: Type                     ! L_Explicit, L_Fixed ...
    integer :: Module                   ! index into modules database
    logical :: ForbidOverspill      = .false.   

    ! This is the maf number passing nearest to the grid point
    integer, dimension(:), pointer    :: MAF          => NULL()
    ! Now the various coordinates in the HGrid, all dimensioned (1,noProfs)
    real(r8), contiguous, pointer :: Phi(:,:)         => NULL()
    real(r8), contiguous, pointer :: GeodLat(:,:)     => NULL()
    real(r8), contiguous, pointer :: Lon (:,:)        => NULL()
    real(r8), contiguous, pointer :: Time(:,:)        => NULL()
    real(r8), contiguous, pointer :: SolarTime(:,:)   => NULL()
    real(r8), contiguous, pointer :: SolarZenith(:,:) => NULL()
    real(r8), contiguous, pointer :: LOSAngle(:,:)    => NULL()
    ! For QTM
    type(QTM_tree_t), allocatable :: QTM_Tree ! for finding things and
                                 ! representing the geolocations of vertices
                                 ! of a QTM that are within or adjacent to a
                                 ! specified polygon
  end type HGrid_T

  ! To construct an array of pointers to HGrid_T.  The reason for this
  ! is that the quantity template has a pointer to an HGrid_T.  When we
  ! add an item to the HGrid_T database, using an array of HGrid_T, all
  ! those pointers become invalid.  So now we allocate each HGrid_T
  ! individually, and grow an array of HGrids_T.
  type HGrids_T
    type(hGrid_t), pointer :: The_HGrid => NULL()
  end type HGrids_T

  ! Put here all the l1boa quantities that we don't
  ! wish to read again and again and again ....
  type HGridGeolocations_T
    double precision, dimension(:,:), pointer :: MAFStartTimeTAI => null()
    double precision, dimension(:,:), pointer :: Orbincl         => null()
    double precision, dimension(:,:), pointer :: GHzGeodAngle    => null()
    double precision, dimension(:,:), pointer :: GHzGeodAlt      => null()
    double precision, dimension(:,:), pointer :: GHzGeodLat      => null()
    double precision, dimension(:,:), pointer :: GHzLosAngle     => null()
    double precision, dimension(:,:), pointer :: GHzLon          => null()
    double precision, dimension(:,:), pointer :: GHzSolarTime    => null()
    double precision, dimension(:,:), pointer :: GHzSolarZenith  => null()
  end type HGridGeolocations_T
  
  type(HGridGeolocations_T), save :: HGridGeolocations

  interface Dump
    module procedure Dump_a_HGrid
    module procedure Dump_HGrids
  end interface

contains ! =========== Public procedures ===================================

 ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( Database, The_HGrid )

    ! This adds an item to the array of pointers to HGrid_T's.  Then it
    ! allocates an HGrid_T as the last element of the new array.  Then it
    ! does a shallow copy of the_hGrid to the allocated hGrid -- so DO NOT
    ! destroy the original one, just nullify all its pointers.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (HGrids_T), dimension(:), pointer :: Database
    type (HGrid_T), intent(in) :: The_HGrid

    ! Local variables
    type (HGrids_T) :: Item
    type (HGrids_T), dimension(:), pointer :: TempDatabase

    ! Add another pointer to HGrid_T to the HGrids_T database
    include "addItemToDatabase.f9h"

    allocate ( database(newSize)%the_hGrid, stat=status )
    call test_allocate ( status, moduleName, "database(newSize)%the_hGrid" )
    database(newSize)%the_hGrid = the_hGrid

    AddHGridToDatabase = newSize
  end function AddHGridToDatabase

  ! -------------------------------------------  copyHGrid  -----
  ! Deep copy all fields from aGrid to hGrid.
  ! Allocates all the pointer components, so a full deep copy results
  ! allowing you to destroy agrid afterwards.
  subroutine copyHGrid ( aGrid, hGrid )
    ! Just does allocates etc.

    use Allocate_Deallocate, only: Allocate_Test

    type (HGrid_T), intent(in)  :: AGRID
    type (HGrid_T), intent(out) :: HGRID
    integer                     :: stat

    ! Executable code
    hGrid%Name                 = aGrid%Name               
    hGrid%masterCoordinate     = aGrid%masterCoordinate   
    hGrid%noProfs              = aGrid%noProfs            
    hGrid%noProfsLowerOverlap  = aGrid%noProfsLowerOverlap
    hGrid%noProfsUpperOverlap  = aGrid%noProfsUpperOverlap
    hGrid%Type                 = aGrid%Type               
    hGrid%forbidOverspill      = aGrid%forbidOverspill    

    call Allocate_Test ( hGrid%maf, hGrid%noProfs, &
      & 'hGrid%maf', ModuleName )
    call Allocate_Test ( hGrid%phi, 1, hGrid%noProfs, &
      & 'hGrid%phi', ModuleName )
    call Allocate_Test ( hGrid%geodLat, 1, hGrid%noProfs, &
      & 'hGrid%geodLat', ModuleName )
    call Allocate_Test ( hGrid%lon, 1, hGrid%noProfs, &
      & 'hGrid%lon', ModuleName )
    call Allocate_Test ( hGrid%time, 1, hGrid%noProfs, &
      & 'hGrid%time', ModuleName )
    call Allocate_Test ( hGrid%solarTime, 1, hGrid%noProfs, &
      & 'hGrid%solarTime', ModuleName )
    call Allocate_Test ( hGrid%solarZenith, 1, hGrid%noProfs, &
      & 'hGrid%solarZenith', ModuleName )
    call Allocate_Test ( hGrid%losAngle, 1, hGrid%noProfs, &
      & 'hGrid%losAngle', ModuleName )

    hGrid%maf          = aGrid%maf
    hGrid%phi          = aGrid%phi        
    hGrid%geodLat      = aGrid%geodLat    
    hGrid%lon          = aGrid%lon        
    hGrid%time         = aGrid%time       
    hGrid%solarTime    = aGrid%solarTime  
    hGrid%solarZenith  = aGrid%solarZenith
    hGrid%losAngle     = aGrid%losAngle   
    if ( allocated(aGrid%QTM_tree) ) then
      allocate ( hGrid%QTM_tree, stat=stat )
      hGrid%losAngle     = aGrid%losAngle   
    endif
  end subroutine copyHGrid

  ! -------------------------------------------  CreateEmptyHGrid  -----
  subroutine CreateEmptyHGrid ( hGrid )
    ! Just does allocates etc.

    use Allocate_Deallocate, only: Allocate_Test

    type (HGrid_T), intent(inout) :: HGRID

    ! Executable code
    call Allocate_Test ( hGrid%maf, hGrid%noProfs, &
      & 'hGrid%maf', ModuleName )
    call Allocate_Test ( hGrid%phi, 1, hGrid%noProfs, &
      & 'hGrid%phi', ModuleName )
    call Allocate_Test ( hGrid%geodLat, 1, hGrid%noProfs, &
      & 'hGrid%geodLat', ModuleName )
    call Allocate_Test ( hGrid%lon, 1, hGrid%noProfs, &
      & 'hGrid%lon', ModuleName )
    call Allocate_Test ( hGrid%time, 1, hGrid%noProfs, &
      & 'hGrid%time', ModuleName )
    call Allocate_Test ( hGrid%solarTime, 1, hGrid%noProfs, &
      & 'hGrid%solarTime', ModuleName )
    call Allocate_Test ( hGrid%solarZenith, 1, hGrid%noProfs, &
      & 'hGrid%solarZenith', ModuleName )
    call Allocate_Test ( hGrid%losAngle, 1, hGrid%noProfs, &
      & 'hGrid%losAngle', ModuleName )

  end subroutine CreateEmptyHGrid

  ! ---------------------------------------------  L1BGeoLocation  -----
  ! Read the approprriate geolocations from the l1boa file.
  ! The first time through, read from the l1boa, and store in the
  ! appropriate component of HGridGeolocations.
  ! In all subsequent times, simply copy the stored values.
  
  ! We check for stored values by whether the component is associated
  subroutine L1BGeoLocation ( filedatabase, name, moduleStr, values, values2d, &
                            & root, NeverFail )

    use Allocate_Deallocate, only: Allocate_Test
    use Dump_0, only: Dump
    use HighOutput, only: OutputNamedValue
    use L1BData, only: L1BData_T, &
      & AssembleL1BQtyName, ConvertL1BData, DeallocateL1BData, ReadL1BData
    use MLSCommon, only: MLSFile_T
    use MLSFiles, only: GetMLSFileByType
    use MLSFillValues, only: IsFillValue, Monotonize
    use MLSKinds, only: Rk => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use MLSStrings, only: Lowercase, StrEq
    use Output_M, only: Output
    use Toggles, only: Gen, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    ! Args
    type (MLSFile_T), dimension(:), pointer     :: fileDatabase
    character(len=*), intent(in)                :: name
    character(len=*), intent(in)                :: moduleStr
    real(rk), dimension(:), pointer, optional   :: values
    real(rk), dimension(:,:), pointer, optional :: values2d
    integer, intent(in), optional :: Root ! For tracing
    logical, intent(in), optional :: NeverFail ! For ReadL1BData
    ! Local variables
    type (L1BData_T) :: l1bField ! L1B data
    logical, parameter                    :: DeeBug = .false.
    integer                               :: HdfVersion
    type (MLSFile_T), pointer             :: L1BFile
    character(len=64)                     :: L1bItemName
    integer                               :: Me = -1 ! String index for tracing
    logical                               :: MyNeverFail
    integer                               :: NoFreqs
    integer                               :: NoMAFs
    character(len=64)                     :: ReadItemName
    logical, parameter                    :: ShowName = .false.
    character(len=64)                     :: SnippedName
    integer                               :: Status
    logical                               :: Verbose 
    ! Executable
    if ( present(root) ) call trace_begin ( me, "L1BGeoLocation", root, &
      & cond=toggle(gen) .and. levels(gen) > 2 )
    myNeverFail = .false.
    if ( present(neverFail) ) myNeverFail = neverFail
    l1bItemName = adjustl(lowercase(name))
    ! Sometimes we're called as GHz/Name; other times as tpName
    if ( l1bItemName(1:2) == 'tp' ) l1bItemName = trim(moduleStr) // l1bItemName(3:)
    if ( l1bItemName(1:2) == 'sc' ) l1bItemName = 'sc/' // l1bItemName(3:)
    ! `1st--check if we have already read this from the l1boa file
    if ( DeeBug ) call outputnamedValue( 'name', name )
    if ( DeeBug ) call outputnamedValue( 'module name', moduleStr )
    if ( DeeBug ) call outputnamedValue( 'l1bItemName (before select case)', l1bItemName )
    snippedName = snipModuleStr(l1bItemName)
    if ( showName ) call output ( '############### ' // trim(snippedName), advance='yes' )
    select case ( snippedName )
    case ('mafstarttimetai')
      if ( associated ( HGridGeolocations%MAFStartTimeTAI ) ) then
        call recallValues( HGridGeolocations%MAFStartTimeTAI, values, values2d )
        go to 9
      endif
    case ('orbincl')
      if ( associated ( HGridGeolocations%OrbIncl ) ) then
        call recallValues( HGridGeolocations%OrbIncl, values, values2d )
        go to 9
      endif
    case ('geodangle   ')
      if ( associated ( HGridGeolocations%GHzGeodAngle ) ) then
        call recallValues( HGridGeolocations%GHzGeodAngle, values, values2d )
        go to 9
      endif
    case ('geodalt     ')
      if ( associated ( HGridGeolocations%GHzGeodAlt ) ) then
        call recallValues( HGridGeolocations%GHzGeodAlt, values, values2d )
        go to 9
      endif
    case ('geodlat     ')
      if ( associated ( HGridGeolocations%GHzGeodLat ) ) then
        call recallValues( HGridGeolocations%GHzGeodLat, values, values2d )
        go to 9
      endif
    case ('solartime   ')
      if ( associated ( HGridGeolocations%GHzSolarTime ) ) then
        call recallValues( HGridGeolocations%GHzSolarTime, values, values2d )
        go to 9
      endif
    case ('solarzenith   ')
      if ( associated ( HGridGeolocations%GHzSolarZenith ) ) then
        call recallValues( HGridGeolocations%GHzSolarZenith, values, values2d )
        go to 9
      endif
    case ('losangle   ')
      if ( associated ( HGridGeolocations%GHzlosAngle ) ) then
        call recallValues( HGridGeolocations%GHzlosAngle, values, values2d )
        go to 9
      endif
    case ('lon   ')
      if ( associated ( HGridGeolocations%GHzlon ) ) then
        call recallValues( HGridGeolocations%GHzlon, values, values2d )
        go to 9
      endif
    end select
    ! If we get here, 
    ! it is because we hadn't read and saved the proper data type
    L1BFile => GetMLSFileByType( filedatabase, content='l1boa' )
    hdfversion = L1BFile%HDFVersion
    if ( DeeBug ) call outputnamedValue( 'name', name )
    if ( DeeBug ) call outputnamedValue( '1st l1bitemname', l1bitemname )
    if ( index( lowercase(name), 'ghz') > 0 ) then
      ! readItemName = "/" // adjustl(Name)
      call output( 'We should have stopped specifying GHz by now', advance='yes' )
      call outputNamedValue( 'name', trim(Name) )
      call outputNamedValue( 'moduleStr', trim(moduleStr) )
      if ( moduleStr /= 'none' .and. moduleStr /= 'GHz' .and. &
        & len_trim(moduleStr) > 1 ) then
        readItemName = trim(moduleStr) // '/' // snipModuleStr(Name, 'GHz')
      else
        readItemName = Name
      endif
      call outputNamedValue( 'readItemName to read', trim(readItemName) )
    elseif ( index( name, 'tp') > 0 ) then
      ! l1bItemName = AssembleL1BQtyName ( 'GHz.' // Name, hdfVersion, .false. )
      if ( DeeBug ) then
        call outputNamedValue( 'AssembleL1BQtyName Name:', trim(Name) )
        call outputNamedValue( 'ModuleStr:', trim(ModuleStr) )
        call outputNamedValue( 'hdfVersion:', hdfVersion )
      endif
      readItemName = AssembleL1BQtyName ( Name(3:), hdfVersion, .true., moduleStr )
    else
      readItemName = AssembleL1BQtyName ( Name, hdfVersion, .false. )
    endif
    if ( DeeBug ) call outputnamedValue( '2nd readItemName', readItemName )
    ! We must not pad because that inserts unwanted Fill values 
    ! into the l1b data fields!
    ! Before we began using L1BGeoLocation, we used to read each
    ! chunk's worth of l1b data by setting FirstMAF and lastMAF
    ! which auttomatically set dontPad=.true.
    call ReadL1BData ( L1BFile, readItemName, l1bField, noMAFs, status, &
      & dontpad=.true., neverFail=neverFail )
    if ( myNeverFail .and. status /= 0 ) go to 9
    if ( .not. associated(l1bField%dpField) ) call ConvertL1BData( l1bField )
    if ( any(isFillValue(l1bField%dpField) ) .and. .not. & 
      & any(StrEq( (/ '/GHz/Lon ', 'scOrbIncl' /), trim(readItemname) ) ) ) then 
      call output( 'Fill values among ' // trim(readItemName), advance='yes' ) 
      call MLSMessage ( MLSMSG_Warning, trim(ModuleName) // '%L1BGeoLocation', & 
        & 'Required monotonization' ) 
      verbose = ( trim(readItemname) == '/GHz/Lon' ) & 
        & .or. ( trim(readItemname) == '/GHz/GeodAngle' ) 
      if ( verbose ) call dump( l1bField%dpField, 'lons before monotony' ) 
      call Monotonize( l1bField%dpField ) 
      if ( verbose ) call dump( l1bField%dpField, 'lons after monotony' ) 
    endif 

    noFreqs = size(l1bField%dpField,2)
    if ( present(values) ) then
      nullify(values)
      call Allocate_test ( values, noMAFs, 'values of ' // trim(name), ModuleName )
      values = l1bField%dpField(1,1,:)
    endif
    if ( present(values2d) ) then
      nullify(values2d)
      call Allocate_test ( values2d, noFreqs, noMAFs, '2d values of ' // trim(name), ModuleName )
      values2d = l1bField%dpField(1,:,:)
    endif
    ! Now before leaving, save these values for the future,
    ! eliminating the need to read the same dataset each time through
    select case ( snippedName )
    case ('mafstarttimetai')
      call Allocate_test ( HGridGeolocations%MAFStartTimeTAI, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%MAFStartTimeTAI = l1bField%dpField(1,:,:)
    case ('geodangle   ')
      call Allocate_test ( HGridGeolocations%GHzGeodAngle, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodAngle = l1bField%dpField(1,:,:)
    case ('geodalt     ')
      call Allocate_test ( HGridGeolocations%GHzGeodAlt, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodAlt = l1bField%dpField(1,:,:)
    case ('geodlat     ')
      call Allocate_test ( HGridGeolocations%GHzGeodLat, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodLat = l1bField%dpField(1,:,:)
    case ('losangle   ')
      call Allocate_test ( HGridGeolocations%GHzlosangle, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzlosangle = l1bField%dpField(1,:,:)
    case ('lon   ')
      call Allocate_test ( HGridGeolocations%GHzlon, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzlon = l1bField%dpField(1,:,:)
    case ('orbincl')
      call Allocate_test ( HGridGeolocations%orbincl, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%orbincl = l1bField%dpField(1,:,:)
    case ('solartime   ')
      call Allocate_test ( HGridGeolocations%GHzSolarTime, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzSolarTime = l1bField%dpField(1,:,:)
    case ('solarzenith   ')
      call Allocate_test ( HGridGeolocations%GHzsolarzenith, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzsolarzenith = l1bField%dpField(1,:,:)
    end select
    call deallocateL1BData ( l1bField ) ! Avoid memory leaks
9   if ( present(root) ) call trace_end ( "L1BGeoLocation", &
      & cond=toggle(gen) .and. levels(gen) > 2 )
  contains
    function snipModuleStr ( itemName, moduleGHz ) result ( bareName )
      character(len=*), intent(in)                    :: itemname
      character(len=32)                               :: bareName
      character(len=*), intent(in), optional          :: moduleGHz
      ! Internal variables
      character(len=32)                               :: modName
      integer                                         :: i
      ! Executable
      modName = moduleStr
      if ( present(moduleGHz) ) modName = moduleGHz
      i = index( itemname, trim(modName) // '/')
      if ( i < 1 ) then
        bareName = itemName
      else
        bareName = itemName( i+len_trim(modName)+1 : ) ! To account for "/"
      endif
      i = index( barename, trim(modName) )
      if ( i < 1 ) then
        bareName = bareName
      else
        bareName = bareName( i+len_trim(modName) : )
      endif
      if ( bareName(1:3) == 'sc/' ) barename = bareName(4:)
      if ( bareName(1:3) == 'sc' ) barename = bareName(3:)
    end function snipModuleStr
    subroutine recallValues( array, values, values2d )
      ! Args
      double precision, dimension(:,:)            :: array
      real(rk), dimension(:), pointer, optional   :: values
      real(rk), dimension(:,:), pointer, optional :: values2d
      ! Internal variables
      integer :: noFreqs
      integer :: noMAFs
      ! Executable
      noFreqs = size( array, 1 )
      noMAFs = size( array, 2 )
      if ( present(values) ) then
        nullify(values)
        call Allocate_test ( values, noMAFs, 'values of array', ModuleName )
        values = array(1,:)
      endif
      if ( present(values2d) ) then
        nullify(values2d)
        call Allocate_test ( values2d, noFreqs, noMAFs, '2d values of array', &
          & ModuleName )
        values2d = array(:,:)
      endif
    end subroutine recallValues
  end subroutine L1BGeoLocation

  ! ---------------------------------------------  L1BSubsample  -----
  ! Pull out just the array elements corresponding to startMAF and LastMAF
  subroutine L1BSubsample ( chunk, fullArray, fullArray2d, values, values2d )
    use Allocate_Deallocate, only: Allocate_Test
    use HighOutput, only: BeVerbose, OutputNamedValue
    use MLSCommon, only: MLSChunk_t
    use MLSKinds, only: Rk => R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    ! Args
    type(MLSChunk_t), intent(in)                 :: chunk
    real(rk), dimension(:), optional             :: fullArray
    real(rk), dimension(:,:), optional           :: fullArray2d
    real(rk), dimension(:), pointer, optional    :: values
    real(rk), dimension(:,:), pointer, optional  :: values2d
    logical                                      :: verbose
    ! Local variables
    integer :: noFreqs
    integer :: noChunkMAFs
    ! Executable
    verbose = BeVerbose( 'hgrid', 1 ) ! .or. .true.
    noChunkMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
    if ( verbose ) then
      call outputnamedValue( 'noChunkMAFs', noChunkMAFs )
      call outputnamedValue( '/chunk%firstMAFIndex+1, chunk%lastMAFIndex+1/', &
        &                    (/chunk%firstMAFIndex+1, chunk%lastMAFIndex+1/) )
    endif
    if ( present(values) ) then
      if ( .not. present(fullArray) )  &
        & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '%L1BSubsample', &
        & 'need FullArray' )
      nullify(values)
      call Allocate_test ( values, noChunkMAFs, 'values of array', ModuleName )
      values = fullArray ( chunk%firstMAFIndex+1:chunk%lastMAFIndex+1 )
    endif
    if ( present(values2d) ) then
      noFreqs = size(fullArray2d, 1)
      if ( .not. present(fullArray2d) )  &
        & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '%L1BSubsample', &
        & 'need FullArray2d' )
      nullify(values2d)
      call Allocate_test ( values2d, noFreqs, noChunkMAFs, '2d values of array', ModuleName )
      values2d = fullArray2d ( :, chunk%firstMAFIndex+1:chunk%lastMAFIndex+1 )
    endif
  end subroutine L1BSubsample

  ! -------------------------------------------------  TrimHGrid  ------
  subroutine TrimHGrid ( hGrid, side, NOTODELETE )

    use Allocate_Deallocate, only: Allocate_Test
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error

    type (HGrid_T), intent(inout) :: HGrid
    integer, intent(in) :: SIDE         ! -1 = lower, 1 = upper
    integer, intent(in) :: NOTODELETE ! How many to delete, default all
    ! Local variables
    integer :: newNoProfs
    integer, dimension(hGrid%noProfs) :: inttemp
    real(r8), dimension(hGrid%noProfs) :: temp
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
    if ( newNoProfs < 0 ) call MLSMessage ( &
        & MLSMSG_Error, ModuleName, 'Too many profiles to delete' )

    hGrid%noProfs = newNoProfs

    ! Now allocate each entry and trim it
    inttemp = hGrid%maf(:)               ! ------------------------- Maf
    call Allocate_Test ( hGrid%maf, newNoProfs, 'hGrid%maf', ModuleName)
    hGrid%maf(:) = inttemp ( first : last )
    temp = hGrid%phi(1,:)               ! ------------------------- Phi
    call Allocate_Test ( hGrid%phi, 1, newNoProfs, 'hGrid%phi', ModuleName)
    hGrid%phi(1,:) = temp ( first : last )
    temp = hGrid%geodLat(1,:)           ! ------------------------- GeodLat
    call Allocate_Test ( hGrid%geodLat, 1, newNoProfs, 'hGrid%geodLat', ModuleName)
    hGrid%geodLat(1,:) = temp ( first : last )
    temp = hGrid%lon(1,:)               ! ------------------------- Lon
    call Allocate_Test ( hGrid%lon, 1, newNoProfs, 'hGrid%lon', ModuleName)
    hGrid%lon(1,:) = temp ( first : last )
    temp = hGrid%time(1,:)              ! ------------------------- Time
    call Allocate_Test ( hGrid%time, 1, newNoProfs, 'hGrid%time', ModuleName)
    hGrid%time(1,:) = temp ( first : last )
    temp = hGrid%solarTime(1,:)         ! ------------------------- SolarTime
    call Allocate_Test ( hGrid%solarTime, 1, newNoProfs, 'hGrid%solarTime', ModuleName)
    hGrid%solarTime(1,:) = temp ( first : last )
    temp = hGrid%solarZenith(1,:)       ! ------------------------- SolarZenith
    call Allocate_Test ( hGrid%solarZenith, 1, newNoProfs, 'hGrid%solarZenith', ModuleName)
    hGrid%solarZenith(1,:) = temp ( first : last )
    temp = hGrid%losAngle(1,:)          ! ------------------------- LosAngle
    call Allocate_Test ( hGrid%losAngle, 1, newNoProfs, 'hGrid%losAngle', ModuleName)
    hGrid%losAngle(1,:) = temp ( first : last )

  end subroutine TrimHGrid
    
  ! ---------------------------------------  DestroyHGridContents  -----
  subroutine DestroyHGridContents ( hGrid )

    use Allocate_Deallocate, only: Deallocate_Test

  ! This routine destroys the information associated with an hGrid

    ! Dummy arguments
    type (HGrid_T), intent(inout) :: hGrid

    ! Executable code
    
    call deallocate_test ( hGrid%maf, 'hGrid%maf', ModuleName )
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

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

  ! This subroutine destroys a quantity template database

    ! Dummy argument
    type (HGrids_T), dimension(:), pointer :: database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: hGridIndex, s, status

    if ( associated(database) ) then
      do hGridIndex=1,SIZE(database)
        call DestroyHGridContents ( database(hGridIndex)%the_hGrid )
        s = storage_size(database(hGridIndex)%the_hGrid) / 8
        addr = 0
        if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
        deallocate ( database(hGridIndex)%the_hGrid, stat=status )
        call test_deallocate ( status, ModuleName, &
          & "database(hGridIndex)%the_hGrid", s, address=addr )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if

  end subroutine DestroyHGridDatabase

  ! ------------------------------------------------  Dump_a_HGrid  -----
  subroutine Dump_a_HGrid ( aHGRID, Details, ZOT )
    use Dates_Module, only: TAI93s2hid
    use HighOutput, only: OutputNamedValue
    use Intrinsic, only: Lit_Indices, L_QTM
    use Optional_m, only: Default
    use Output_m, only: Blanks, NewLine, Output
    use QTM_Output, only: Dump_QTM_Tree
    use String_Table, only: Display_String, IsStringInTable

    type(hGrid_T), intent(in) :: aHGRID
    logical, intent(in), optional :: ZOT      ! Dump QTM coordinates in ZOT
                                              ! Default false
    integer, intent(in), optional :: Details  ! < 1 => Don't dump QTM search tree
                                              ! > 1 => Dump sons of QTM vertices
                                              ! Default zero

    integer :: IERR
    integer :: J
    ! Executable

    IERR = 0
    call output ( 'Name = ', advance='no' )
    if ( isStringInTable(aHgrid%name) ) then
      call display_string ( aHgrid%name, ierr=ierr )
      if ( ierr /= 0 ) call output ( ' (not found in string table)' )
    else
      call output('(unknown)' )
    end if
    call newLine
    call output ( aHgrid%noProfs, before=' noProfs = ', advance='no' )
    call output ( aHgrid%noProfsLowerOverlap, before=' lowerOverlap = ', advance='no' )
    call output ( aHgrid%noProfsUpperOverlap, before=' upperOverlap = ', advance='no' )
    call output ( ' masterCoordinate = ' )
    if ( .not. isStringInTable(aHgrid%masterCoordinate, lit_indices) ) then
      call output ( '0', advance='yes' )
    else
      call display_string ( lit_indices(aHgrid%masterCoordinate), advance='yes' )
    end if
    if ( .not. isStringInTable( aHgrid%type, lit_indices ) ) then
      call outputNamedValue ( 'type', 0 )
    else
      call display_string ( lit_indices(aHgrid%type), before=' type = ' )
    end if
    if ( aHGrid%type == l_QTM ) &
      & call output ( aHGrid%QTM_tree%level, before=' tree with level = ' )
    call newLine
    if ( aHgrid%noProfs > 0 .and. default(details,0) > 0 ) then
      call output ( ' prof       phi       geodLat           lon', advance='no' )
      call output ( '      hours in day     solarTime   solarZenith', advance='no' )
      call output ( '      losAngle  nearest maf', advance='yes' )
      do j = 1, aHgrid%noProfs
        if ( &                                                      
          & ( aHgrid%noProfsLowerOverlap > 0  .and. &               
          & j == aHgrid%noProfsLowerOverlap + 1 ) &                 
          & .or. &                                                  
          & ( aHgrid%noProfsUpperOverlap > 0 .and. &                
          & aHgrid%noProfs - j == aHgrid%noProfsUpperOverlap - 1 ) &
          & ) &                                                     
          & call blanks( 80, fillChar='-', advance='yes' )          
        call output ( j, places=5 )
        call output ( aHgrid%phi(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%geodLat(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%lon(1,j), '(1x,1pg13.6)' )
        call output ( tai93s2hid( &                                       
    &         aHgrid%time(1,j), leapsec=.true. ), '(1x,1pg13.6)' )
        call output ( aHgrid%solarTime(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarZenith(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%losAngle(1,j), '(1x,1pg13.6)', advance='no' )
        call output ( aHgrid%maf(j), places=6, advance='yes' )
      end do
    end if
    if ( aHGrid%type == l_QTM .and. default(details,0) > 0 ) then
      call dump_QTM_tree ( aHGrid%QTM_tree, latLon=.not. default(zot,.false.), &
        & details=details )
    end if
  end subroutine Dump_a_HGrid

  ! ------------------------------------------------  Dump_HGrids  -----
  subroutine Dump_HGrids ( HGrids, Details, ZOT )
    use Output_M, only: Output
    type(hGrids_T), intent(in) :: HGrids(:)
    integer, intent(in), optional :: Details
    logical, intent(in), optional :: ZOT
    integer :: I
    call output ( size(hgrids), before='HGRIDS: SIZE = ', advance='yes' )
    do i = 1, size(hgrids)
      call output ( i, 4, after=': ' )
      call dump ( hGrids(i)%the_hGrid, details, ZOT )
    end do
  end subroutine Dump_HGrids

  ! ---------------------------------------- FindClosestMatch ---
  integer function FindClosestMatch ( reference, sought, instance )
    use HighOutput, only: OutputNamedValue
    use Hunt_m, only: Hunt
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    ! Given a sought quantity, 
    ! the profile in reference is found
    ! that is closest to profile 'instance' in sought quantity
    real(r8), dimension(:), intent(in) :: REFERENCE ! e.g. temperature
    real(r8), dimension(:,:), intent(in) :: SOUGHT ! e.g. ptan, radiance
    integer, intent(in) :: INSTANCE
    ! Local variables
    integer :: BESTGUESS                ! A guessed index
    integer :: FIRSTGUESS               ! A guessed index
    integer :: HIGHGUESS                ! A guessed index
    integer :: LOWGUESS                 ! A guessed index
    real(r8) :: BESTCOST                ! A cost for a guess
    real(r8) :: COST                    ! A cost for a guess
    logical, parameter                    :: DeeBug = .false.

    ! Executable code
    ! Get starting Guess from Hunt
    call Hunt ( reference, sought(1,instance), firstGuess, &
      & start=max(min(instance,size(reference)),1), &
      & allowTopValue=.true., nearest=.true. )

    ! Now check for better ones either side
    lowGuess = max ( firstGuess-1, 1 )
    highGuess = min ( firstGuess+1, size(reference) )
    bestCost = 0.0
    if ( DeeBug ) then
      call outputNamedValue( 'lowGuess', lowGuess )
      call outputNamedValue( 'highGuess', highGuess )
      call outputNamedValue( 'instance', instance )
    endif
    do firstGuess = lowGuess, highGuess
      if ( firstGuess < lbound(reference, 1) .or. firstGuess > ubound(reference, 1) ) then
        call outputNamedValue( 'bound(reference)', (/ lbound(reference), lbound(reference) /) )
        call outputNamedValue( 'firstGuess', firstGuess )
        call MLSMessage ( &
          & MLSMSG_Error, ModuleName, 'firstGuess outside range' )
      elseif ( instance < lbound(sought, 2 ) .or. instance > ubound(sought, 2 ) ) then
        call outputNamedValue( 'lbound(sought)', lbound(sought) )
        call outputNamedValue( 'ubound(sought)', ubound(sought) )
        call outputNamedValue( 'instance', instance )
        call MLSMessage ( &
          & MLSMSG_Error, ModuleName, 'instance outside range' )
      endif
      cost = sum ( abs ( reference(firstGuess) - sought(:,instance) ) )
      if ( ( firstGuess == lowGuess ) .or. ( cost < bestCost ) ) then
        bestGuess = firstGuess
        bestCost = cost
      end if
    end do

    FindClosestMatch = bestGuess
  end function FindClosestMatch

  ! ---------------------------------------- NullifyHGrid -----
  subroutine NullifyHGrid ( IntentionallyNotUsed )
    ! Given a hGrid, nullify all the pointers associated with it
    type ( HGrid_T ), intent(out) :: IntentionallyNotUsed

    ! Executable code.  Nothing is needed since intent(out) causes
    ! default initialization, which nullifies the pointers because
    ! they all have default initialization => NULL()

  end subroutine NullifyHGrid

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HGridsDatabase

! $Log$
! Revision 2.46  2021/07/08 23:30:38  pwagner
! Require -SGrid2 for verbose
!
! Revision 2.45  2020/01/27 18:30:44  pwagner
! Works wsith new AssembleL1BQtyName
!
! Revision 2.44  2020/01/09 22:22:08  pwagner
! Now wont try to Monotonize sids l1boa scOrbIncl
!
! Revision 2.43  2018/08/04 00:44:16  vsnyder
! Spiff a warning message
!
! Revision 2.42  2018/08/03 23:26:10  vsnyder
! Add Root (for tracing) and NeverFail (for debug printing) to L1BGeoLocation.
! Alphabetize local variable declarations.  Pass NeverFail to ReadL1BData.
! Make sure always to end tracing.
!
! Revision 2.41  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.40  2016/10/19 00:11:51  pwagner
! Try to avoid certain crashes QTM HGrids and ASMLS data
!
! Revision 2.39  2016/10/11 23:28:37  pwagner
! Removed stop statement; downgraded to just a warning
!
! Revision 2.38  2016/10/01 01:37:28  vsnyder
! Make QTM_Tree component of HGrid_t allocatable
!
! Revision 2.37  2016/09/23 03:08:37  vsnyder
! Get optional argument values using Optional_m
!
! Revision 2.36  2016/09/23 02:41:47  vsnyder
! Remove USE for unused use names
!
! Revision 2.35  2016/09/23 01:34:56  vsnyder
! Remove QTM_Lats; moved to QTM_Tree_t in Generate_QTM_m
!
! Revision 2.34  2016/09/13 20:09:13  vsnyder
! Get Dump_QTM from its new home in QTM_Output
!
! Revision 2.33  2016/08/23 18:34:44  pwagner
! ShowName not needed any longer
!
! Revision 2.32  2016/08/23 00:41:55  vsnyder
! Components within or adjacent to the polygon are now within the QTM_Tree_t
! structure instead of the HGrid_t structure.
!
! Revision 2.31  2016/08/20 01:11:39  vsnyder
! Cannonball polishing
!
! Revision 2.30  2016/08/12 00:35:47  pwagner
! Seems to restore tthe gold brick
!
! Revision 2.29  2016/08/09 18:16:30  pwagner
! Survives encounter with non-satellite data
!
! Revision 2.28  2016/07/28 01:34:55  vsnyder
! Remove unreferenced USE
!
! Revision 2.27  2016/05/18 01:34:37  vsnyder
! HGridsDatabase.f90
!
! Revision 2.26  2016/05/17 00:14:40  pwagner
! Tries harder not to bomb in display_string
!
! Revision 2.25  2016/04/15 22:24:43  pwagner
! Tried to fix latest cause of Fill Values among l1b fields
!
! Revision 2.24  2016/02/26 02:05:25  vsnyder
! Add QTM support
!
! Revision 2.23  2016/01/11 23:12:58  pwagner
! Skip Attempting to monotonize longitudes containing Fills
!
! Revision 2.22  2015/07/14 23:22:27  pwagner
! Prevent insertion of Fill values in L1BGeolocations when gaps occur in counterMAF
!
! Revision 2.21  2015/06/25 00:37:52  pwagner
! Stop printing low and high Guesses
!
! Revision 2.20  2015/06/19 20:35:06  pwagner
! Fixed mixup in defining noChunkMAFs
!
! Revision 2.19  2015/06/19 00:34:19  pwagner
! Many changes to speed up computing HGrid offsets
!
! Revision 2.18  2015/05/27 22:39:39  vsnyder
! Get Hunt from Hunt_m to avoid circular dependence
!
! Revision 2.17  2015/03/31 21:03:55  pwagner
! ForbidOverspill now a field of HGrid; improved Dump
!
! Revision 2.16  2015/03/28 01:03:12  vsnyder
! Added type component.  Made geolocation components contiguous so they can
! be rank-remapping pointer targets.
! Added stuff to trace allocate/deallocate addresses -- mostly commented out
! because NAG build 1017 doesn't yet allow arrays as arguments to C_LOC.
!
! Revision 2.15  2014/09/04 23:40:28  vsnyder
! More complete and accurate allocate/deallocate size tracking.
! Convert some local pointer temps to automatic.
!
! Revision 2.14  2014/08/07 22:44:28  vsnyder
! Set default initialization of MasterCoordinate to L_GeodAngle instead
! of zero.  Spiff up the dump a little bit.
!
! Revision 2.13  2014/04/24 23:50:25  pwagner
! Added masterCoordinate component
!
! Revision 2.12  2013/10/01 22:16:45  pwagner
! Added maf component to HGrid_T
!
! Revision 2.11  2013/08/12 23:47:07  vsnyder
! Default initialize more stuff in HGrid_t
!
! Revision 2.10  2013/03/01 01:04:51  pwagner
! Get R8 from MLSKinds
!
! Revision 2.9  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.8  2008/05/02 00:38:48  vsnyder
! Simplify NullifyHGrid
!
! Revision 2.7  2005/09/21 23:11:54  pwagner
! Unnecessary changes
!
! Revision 2.6  2005/08/25 20:18:05  pwagner
! Protect against crashing when dumping anonymous HGrids
!
! Revision 2.5  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/05/20 19:48:25  vsnyder
! Move Dump*HGrid here from dumper
!
! Revision 2.3  2004/05/18 01:05:06  vsnyder
! Delete unused MAFIndex and MAFCounter fields from HGrid_T
!
! Revision 2.2  2003/07/07 20:20:28  livesey
! New FindClosestMatch routine
!
! Revision 2.1  2003/06/20 19:34:45  pwagner
! Quanities now share grids stored separately in databses
!
