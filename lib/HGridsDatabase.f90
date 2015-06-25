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

  implicit none
  private

  public :: HGrid_T, HGridGeolocations_T, HGridGeolocations
  public :: addHGridtodatabase, copyHGrid, createEmptyHGrid, destroyHGridContents, &
    & destroyHGridDatabase, Dump, findClosestMatch, &
    & L1BGeoLocation, L1BSubsample, nullifyHGrid, &
    & trimHGrid

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  type HGrid_T
    integer :: Name                 = 0    ! String index of name.            
    integer :: masterCoordinate     = l_GeodAngle ! Its lit index;
    integer :: noProfs                 ! Number of profiles in this grid  
    integer :: noProfsLowerOverlap  = 0 ! Number of profiles in the lower overlap
    integer :: noProfsUpperOverlap  = 0 ! Number of profiles in the upper overlap
    integer :: Type                    ! L_Explicit, L_Fixed ...
    logical :: forbidOverspill      = .false.   

    ! This is the maf number passing nearest to the grid point
    integer, dimension(:), pointer    :: maf         => NULL()
    ! Now the various coordinates in the HGrid, all dimensioned (1,noProfs)
    real(r8), contiguous, pointer :: phi(:,:)         => NULL()
    real(r8), contiguous, pointer :: geodLat(:,:)     => NULL()
    real(r8), contiguous, pointer :: lon (:,:)        => NULL()
    real(r8), contiguous, pointer :: time(:,:)        => NULL()
    real(r8), contiguous, pointer :: solarTime(:,:)   => NULL()
    real(r8), contiguous, pointer :: solarZenith(:,:) => NULL()
    real(r8), contiguous, pointer :: losAngle(:,:)    => NULL()
  end type HGrid_T

  ! Put here all the 
  ! l1boa quantities that we don't wish to read again and again and again ..
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

  interface DUMP
    module procedure DUMP_a_HGRID
    module procedure DUMP_HGRIDS
  end interface

contains ! =========== Public procedures ===================================

 ! -----------------------------------------  AddHGridToDatabase  -----
  integer function AddHGridToDatabase ( database, item )

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    ! Dummy arguments
    type (HGrid_T), dimension(:), pointer :: database
    type (HGrid_T), intent(in) :: item

    ! Local variables
    type (HGrid_T), dimension(:), pointer :: tempDatabase

    include "addItemToDatabase.f9h"

    AddHGridToDatabase = newSize
  end function AddHGridToDatabase

  ! -------------------------------------------  copyHGrid  -----
  ! Copy all fields from aGrid to hGrid
  ! Allocates all the pointer components, so a full copy results
  ! allowing you to destroy agrid afterwards
  subroutine copyHGrid ( aGrid, hGrid )
    ! Just does allocates etc.

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST

    type (HGrid_T), intent(in)  :: AGRID
    type (HGrid_T), intent(out) :: HGRID

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
  end subroutine copyHGrid

  ! -------------------------------------------  CreateEmptyHGrid  -----
  subroutine CreateEmptyHGrid ( hGrid )
    ! Just does allocates etc.

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST

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
  subroutine L1BGeoLocation ( filedatabase, name, values, values2d )

    use allocate_deallocate, only: allocate_test
    use HighOutput, only: outputNamedValue
    use L1BData, only: deallocateL1BData, L1BData_t, readL1BData, &
      & assembleL1BQtyName
    use MLSCommon, only: MLSfile_t
    use MLSFiles, only: getMLSfilebytype
    use MLSKinds, only: rk => r8
    use MLSstrings, only: lowercase
    ! Args
    type (MLSFile_T), dimension(:), pointer     :: fileDatabase
    character(len=*), intent(in)                :: name
    real(rk), dimension(:), pointer, optional   :: values
    real(rk), dimension(:,:), pointer, optional :: values2d
    ! Local variables
    type (L1BData_T) :: l1bField ! L1B data
    integer                               :: hdfVersion
    type (MLSFile_T), pointer             :: L1BFile
    character(len=64)                     :: l1bItemName
    integer                               :: noFreqs
    integer                               :: noMAFs
    character(len=64)                     :: readItemName
    integer                               :: status
    logical, parameter                    :: DeeBug = .false.
    ! Executable
    l1bItemName = adjustl(lowercase(name))
    ! Sometimes we're called as GHz/Name; othertimes as tpName
    if ( l1bItemName(1:2) == 'tp' ) l1bItemName = 'ghz/' // l1bItemName(3:)
    if ( l1bItemName(1:2) == 'sc' ) l1bItemName = 'sc/' // l1bItemName(3:)
    ! `1st--check if we have already read this from the l1boa file
    if ( DeeBug ) call outputnamedValue( 'name', name )
    if ( DeeBug ) call outputnamedValue( 'l1bItemName (before select case)', l1bItemName )
    select case (l1bItemName)
    case ('mafstarttimetai')
      if ( associated ( HGridGeolocations%MAFStartTimeTAI ) ) then
        call recallValues( HGridGeolocations%MAFStartTimeTAI, values, values2d )
        return
      endif
    case ('sc/orbincl')
      if ( associated ( HGridGeolocations%OrbIncl ) ) then
        call recallValues( HGridGeolocations%OrbIncl, values, values2d )
        return
      endif
    case ('ghz/geodangle   ')
      if ( associated ( HGridGeolocations%GHzGeodAngle ) ) then
        call recallValues( HGridGeolocations%GHzGeodAngle, values, values2d )
        return
      endif
    case ('ghz/geodalt     ')
      if ( associated ( HGridGeolocations%GHzGeodAlt ) ) then
        call recallValues( HGridGeolocations%GHzGeodAlt, values, values2d )
        return
      endif
    case ('ghz/geodlat     ')
      if ( associated ( HGridGeolocations%GHzGeodLat ) ) then
        call recallValues( HGridGeolocations%GHzGeodLat, values, values2d )
        return
      endif
    case ('ghz/solartime   ')
      if ( associated ( HGridGeolocations%GHzSolarTime ) ) then
        call recallValues( HGridGeolocations%GHzSolarTime, values, values2d )
        return
      endif
    case ('ghz/solarzenith   ')
      if ( associated ( HGridGeolocations%GHzSolarZenith ) ) then
        call recallValues( HGridGeolocations%GHzSolarZenith, values, values2d )
        return
      endif
    case ('ghz/losangle   ')
      if ( associated ( HGridGeolocations%GHzlosAngle ) ) then
        call recallValues( HGridGeolocations%GHzlosAngle, values, values2d )
        return
      endif
    case ('ghz/lon   ')
      if ( associated ( HGridGeolocations%GHzlon ) ) then
        call recallValues( HGridGeolocations%GHzlon, values, values2d )
        return
      endif
    end select
    ! If we got here, 
    ! it was because we hadn't read and saved the proper data type
    L1BFile => GetMLSFileByType( filedatabase, content='l1boa' )
    hdfversion = L1BFile%HDFVersion
    if ( DeeBug ) call outputnamedValue( 'name', name )
    if ( DeeBug ) call outputnamedValue( '1st l1bitemname', l1bitemname )
    if ( index( lowercase(name), 'ghz') > 0 ) then
      readItemName = "/" // adjustl(Name)
    elseif ( index( name, 'tp') > 0 ) then
      ! l1bItemName = AssembleL1BQtyName ( 'GHz.' // Name, hdfVersion, .false. )
      readItemName = AssembleL1BQtyName ( 'GHz.' // Name, hdfVersion, .false. )
    else
      readItemName = AssembleL1BQtyName ( Name, hdfVersion, .false. )
    endif
    if ( DeeBug ) call outputnamedValue( '2nd readItemName', readItemName )
    call ReadL1BData ( L1BFile, readItemName, l1bField, noMAFs, status )
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
    select case (lowercase(l1bItemName))
    case ('mafstarttimetai')
      call Allocate_test ( HGridGeolocations%MAFStartTimeTAI, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%MAFStartTimeTAI = l1bField%dpField(1,:,:)
    case ('ghz/geodangle   ')
      call Allocate_test ( HGridGeolocations%GHzGeodAngle, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodAngle = l1bField%dpField(1,:,:)
    case ('ghz/geodalt     ')
      call Allocate_test ( HGridGeolocations%GHzGeodAlt, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodAlt = l1bField%dpField(1,:,:)
    case ('ghz/geodlat     ')
      call Allocate_test ( HGridGeolocations%GHzGeodLat, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzGeodLat = l1bField%dpField(1,:,:)
    case ('ghz/losangle   ')
      call Allocate_test ( HGridGeolocations%GHzlosangle, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzlosangle = l1bField%dpField(1,:,:)
    case ('ghz/lon   ')
      call Allocate_test ( HGridGeolocations%GHzlon, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzlon = l1bField%dpField(1,:,:)
    case ('sc/orbincl   ')
      call Allocate_test ( HGridGeolocations%orbincl, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%orbincl = l1bField%dpField(1,:,:)
    case ('ghz/solartime   ')
      call Allocate_test ( HGridGeolocations%GHzSolarTime, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzSolarTime = l1bField%dpField(1,:,:)
    case ('ghz/solarzenith   ')
      call Allocate_test ( HGridGeolocations%GHzsolarzenith, noFreqs, noMAFs, &
        & '2d values of ' // trim(name), ModuleName )
      HGridGeolocations%GHzsolarzenith = l1bField%dpField(1,:,:)
    end select
    call deallocateL1BData ( l1bField ) ! Avoid memory leaks
  contains
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
    use allocate_deallocate, only: allocate_test
    use MLSCommon, only: MLSChunk_t
    use MLSKinds, only: rk => r8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    ! Args
    type(MLSChunk_t), intent(in)                 :: chunk
    real(rk), dimension(:), optional             :: fullArray
    real(rk), dimension(:,:), optional           :: fullArray2d
    real(rk), dimension(:), pointer, optional    :: values
    real(rk), dimension(:,:), pointer, optional  :: values2d
    ! Local variables
    integer :: noFreqs
    integer :: noChunkMAFs
    ! Executable
    noChunkMAFs = chunk%lastMAFIndex - chunk%firstMAFIndex + 1
    if ( present(values) ) then
      if ( .not. present(fullArray) )  &
        & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // 'L1BSubsample', &
        & 'need FullArray' )
      nullify(values)
      call Allocate_test ( values, noChunkMAFs, 'values of array', ModuleName )
      values = fullArray ( chunk%firstMAFIndex+1:chunk%lastMAFIndex+1 )
    endif
    if ( present(values2d) ) then
      noFreqs = size(fullArray2d, 1)
      if ( .not. present(fullArray2d) )  &
        & call MLSMessage ( MLSMSG_Error, trim(ModuleName) // 'L1BSubsample', &
        & 'need FullArray2d' )
      nullify(values2d)
      call Allocate_test ( values2d, noFreqs, noChunkMAFs, '2d values of array', ModuleName )
      values2d = fullArray2d ( :, chunk%firstMAFIndex+1:chunk%lastMAFIndex+1 )
    endif
  end subroutine L1BSubsample

  ! -------------------------------------------------  TrimHGrid  ------
  subroutine TrimHGrid ( hGrid, side, NOTODELETE )

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_ERROR

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

    use ALLOCATE_DEALLOCATE, only: DEALLOCATE_TEST

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
    type (HGrid_T), dimension(:), pointer :: database

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: hGridIndex, s, status

    if ( associated(database) ) then
      do hGridIndex=1,SIZE(database)
        call DestroyHGridContents ( database(hGridIndex) )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, ModuleName, "database", s, address=addr )
    end if

  end subroutine DestroyHGridDatabase

  ! ------------------------------------------------  DUMP_A_HGRID  -----
  subroutine DUMP_a_HGRID ( aHGRID )
    use dates_module, only: tai93s2hid 
    use Intrinsic, only: lit_indices
    use output_m, only: blanks, newLine, output
    use string_table, only: display_string
    type(hGrid_T), intent(in) :: aHGRID
    integer :: IERR
    integer :: J
      IERR = 0
      call output ( 'Name = ', advance='no' )
      if ( aHgrid%name > 0 ) then
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
      if ( aHgrid%masterCoordinate == 0 ) then
        call output ( '0', advance='yes' )
      else
        call display_string ( lit_indices(aHgrid%masterCoordinate), advance='yes' )
      end if
      call display_string ( lit_indices(aHgrid%type), before=' type = ', advance='yes' )
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
  &           aHgrid%time(1,j), leapsec=.true. ), '(1x,1pg13.6)' )
        call output ( aHgrid%solarTime(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%solarZenith(1,j), '(1x,1pg13.6)' )
        call output ( aHgrid%losAngle(1,j), '(1x,1pg13.6)', advance='no' )
        call output ( aHgrid%maf(j), places=6, advance='yes' )
      end do
      if ( aHgrid%masterCoordinate > 0 ) then
        call output( ' Master coordinate: ', advance='no' )
        call display_string ( lit_indices(aHgrid%masterCoordinate), ierr=ierr )
        call newLine
      endif
  end subroutine DUMP_a_HGRID

  ! ------------------------------------------------  DUMP_HGRIDS  -----
  subroutine DUMP_HGRIDS ( HGRIDS )
    use OUTPUT_M, only: OUTPUT
    type(hGrid_T), intent(in) :: HGRIDS(:)
    integer :: I
    call output ( size(hgrids), before='HGRIDS: SIZE = ', advance='yes' )
    do i = 1, size(hgrids)
      call output ( i, 4, after=': ' )
      call dump ( hgrids(i) )
    end do
  end subroutine DUMP_HGRIDS

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
