! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=======================================================================================

module DirectWrite_m  ! Write l2gp and l2aux products out to files 
                      ! chunk-by-chunk

!=======================================================================================

    ! 
    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux 
    ! or hdfeos-formatted l2gp for datasets that
    ! are too big to keep all chunks stored in memory
    ! or simply take too much time doing i/o
    ! so instead write them out chunk-by-chunk

  use Allocate_Deallocate, only: Allocate_Test
  use HighOutput, only: BeVerbose, LetsDebug, OutputNamedValue, OutputTable
  use Init_Tables_Module, only: L_GeodAltitude, L_Pressure, L_Zeta, &
    & L_L2gp, L_L2aux, L_L2dgg, L_L2fwm
  use L2ParInfo, only: Parallel
  use Machine, only: USleep
  use MLSCommon, only: Interval_T, MLSFile_T, &
    & InRange
  use MLSKinds, only: Rv
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning
  use MLSFiles, only: HDFversion_4, HDFversion_5, Dump, MLS_Exists, &
    & MLS_CloseFile, MLS_OpenFile, Split_Path_Name
  use MLSFinds, only: FindFirst
  use MLSHDFEOS, only: MLS_Swath_In_File
  use MLSL2Options, only: L2Options, MLSL2Message, WriteFileAttributes
  use MLSStringLists, only: SwitchDetail
  use MLSStrings, only: LowerCase
  use Output_M, only: Blanks, NewLine, Output
  use PCFHdr, only: GlobalAttributes
  use String_Table, only: Get_String
  use Toggles, only: Switches
  use VectorsModule, only: Vector_T, VectorValue_T, Dump

  implicit none

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (data types and parameters)
! DirectData_T                    An L2AUX or L2GP data file

!     (subroutines and functions)
! AddDirectToDatabase             Adds a direct data type to a database of that type
! DestroyDirectDatabase           Deallocates all the arrays for entire database
! DirectRead                      Read the quantity from a file
! DirectWrite                     Write the quantity to a file; or
!                                   Write every quantity in a vector;
!                                   May be either l2gp or l2aux
! Dump                            Prints info on one direct or entire database
! ExpandDirectDB                  Expands a db if needed, or else
!                                   returns the matching direct
! ExpandSDNames                   Expands the array direct%sdNames if needed
! FileNameToID                    returns filename's index if found in db
! SetupNewDirect                  Allocates the arrays for a direct
! === (end of toc) ===

! === (start of api) ===
!     (user-defined types)
! DirectData_T  ( int type, int autoType, int fileIndex, int Handle, &
!    int NSDNames, char* sdNames(:), char* fileNameBase, char* fileName )

!     (subroutines and functions)
! AddDirectToDatabase AddDirectToDatabase( &
!   DirectData_T database(:), DirectData_T item )
! SetupNewDirect ( DirectData_T directData, int NsdNames )
! DestroyDirectDB ( DirectData_T database(:) )
! DirectRead ( MLSFile_T File, VectorValue_T quantity, char* qtyName, &
!    & int chunkNo, [char* options], [int rank] )
! ResizeL2AUXData ( L2AUXData_T l2aux, int newSize )
! int AddL2AUXToDatabase ( *L2AUXData_T DATABASE(:), L2AUXData_T ITEM )
! DestroyL2AUXDatabase ( *L2AUXData_T DATABASE(:) )
! Dump ( l2auxData_T L2aux(:), [char* Name], [int Details], [char* options] )
!    or Dump ( l2auxData_T L2aux, [int Details], [char* options] )
! ReadL2AUXData ( int sd_id, char* quantityname, l2auxData_T l2aux, 
!    [int firstProf], [int lastProf] )
! WriteHDF5Data ( real(r8) array(:,:,:), int l2FileHandle, int returnStatus, 
!    char* sdName )
! WriteL2AUXData ( l2auxData_T l2aux, int l2FileHandle, int returnStatus, 
!    [char* sdName], [int NoMAFS], [log WriteCounterMAF], [char* DimNames] )
! === (end of api) ===

  private
  public :: DirectData_T, &
    & AddDirectToDatabase, &
    & DestroyDirectDatabase, DirectRead, DirectWrite, Dump, &
    & ExpandDirectDB, ExpandSDNames, FileNameToID, &
    & SetupNewDirect

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface DirectRead
    module procedure DirectRead_Quantity
  end interface

  interface DirectWrite
    module procedure DirectWrite_Hdf
    module procedure DirectWrite_L2Aux
    module procedure DirectWrite_L2GP
    module procedure DirectWrite_Quantity
    module procedure DirectWriteVector_Hdf
    module procedure DirectWriteVector_L2GP
    module procedure DirectWriteVector_L2Aux
    module procedure DirectWriteVector_EveryQuantity
  end interface

  interface DUMP
    module procedure DumpDirectWrite
    module procedure DumpDirectDB
  end interface

  type DirectData_T
    integer :: type       ! l_l2aux or l_l2gp
    integer :: autoType   ! 0 if not to be filled automatically, else type
    integer :: fileIndex  ! Index into character tables
    integer :: Handle     ! Set by toolkit
    integer :: NSDNames   ! How many sources
    character(len=80), dimension(:), pointer :: sdNames => null() ! source names
    character(len=1024) :: fileNameBase ! E.g., 'H2O'
    character(len=1024) :: fileName ! E.g., '/data/../MLS..H2O...he5'
  end type DirectData_T
  
  ! logical, parameter :: countEmpty = .true.
  ! integer, parameter :: S2US  = 1000000 ! How many microseconds in a s
  ! integer, parameter :: DELAY = 1*S2US  ! How long to sleep in microseconds
  logical, parameter :: DEEBUG = .false.
  logical, parameter :: MAYWRITEPOSTOVERLAPS = .true.
  ! For Announce_Error
  integer :: ERROR
  integer, save :: lastProfTooBigWarns = 0
  integer, parameter :: MAXNUMWARNS = 40

contains ! ======================= Public Procedures =========================

  !-------------------------------------------  AddDirectToDatabase  -----
  integer function AddDirectToDatabase( DATABASE, ITEM )

    ! This function adds a directly writeabledata type to a database
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where it is put.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer :: DATABASE
    type (DirectData_T), intent(in) :: ITEM

    ! Local variables
    type (DirectData_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddDirectToDatabase = newSize
  end function AddDirectToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a directly writeable database

  subroutine DestroyDirectDatabase ( DATABASE )
    use Allocate_Deallocate, only: Deallocate_Test, Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End

    ! Dummy argument
    type (DirectData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: directIndex, s, status
    integer :: Me = -1       ! String index for trace

    call trace_begin ( me, "DestroyDirectDatabase", cond=toggle(gen) )

    if ( associated(database) ) then
      do directIndex=1, size(DATABASE)
        call deallocate_test ( database(directIndex)%sdNames, &
          "database%sdNames", moduleName )
      end do
      s = size(database) * storage_size(database) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
      deallocate ( database, stat=status )
      call test_deallocate ( status, moduleName, "database", s, address=addr )
    end if
    call trace_end ( "DestroyDirectDatabase", cond=toggle(gen) )
  end subroutine DestroyDirectDatabase

  ! ------------------------------------------ DirectRead_Quantity --------
  ! We Read the quantity from a file, chunk by chunk, overlaps and all
  ! E.g., the qtyname is 'Q', and we have written chunks 1 and 2. Thus the
  ! file layout must be
  ! /
  !   Q/
  !     1/
  !       values
  !       lons
  !       lats
  !         ...
  !     2/
  !       values
  !       lons
  !       lats
  !         ...
  
  subroutine DirectRead_Quantity ( File, quantity, qtyName, &
    & chunkNo, options, rank )

    use Dates_Module, only: Tai93s2hid
    use Dump_0, only: Dump
    use HDF5, only: H5gclose_F, H5gopen_F
    use MLSHDF5, only: GetHDF5Attribute, IsHDF5GroupPresent, &
      & GetHDF5DSRank, LoadFromHDF5DS
    use MLSStringLists, only: OptionDetail
    use MLSStrings, only: WriteIntsToChars
    ! Args:
    type (VectorValue_T), intent(inout) :: QUANTITY
    character(len=*), intent(in) :: qtyName       ! Name of qty in output file
    type(MLSFile_T)                :: File
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    character(len=*), intent(in), optional :: options
    integer, intent(in), optional :: rank

    ! Local variables
    logical :: already_there
    logical :: geolocations
    integer :: grp_id
    integer :: itsRank
    integer :: myRank
    character(len=8) :: chunkStr        ! '1', '2', ..
    integer :: returnStatus
    ! logical, parameter :: DEEBUG = .true.
    logical :: verbose

    ! executable code
    verbose = BeVerbose ( 'direct', -1 )
    myRank = 0
    if ( present(rank) ) then
      if ( rank > 0 .and. rank < 5 ) myRank = rank
      ! call outputNamedValue ( 'input rank', rank )
    endif
    geolocations = optionDetail ( options, 'g' ) /= 'no'
    
    call WriteIntsToChars ( chunkNo, chunkStr )
    chunkStr = adjustl ( chunkStr )
    ! Create or access the SD
    already_There = mls_exists( File%name ) == 0
    if ( .not. already_There ) then
      call outputNamedValue( '  Sorry, file not found', trim(File%name) )
      return
    endif
    call mls_openFile ( File, returnStatus )
    if ( returnStatus /= 0 )  then
      call Dump( File )
      call outputNamedValue( '  Sorry, unable to open file', trim(File%name) )
      return
    endif
    already_there = IsHDF5GroupPresent( File%fileID%f_id, trim(qtyName) )
    if ( .not. already_There ) then
      call Dump( File )
      call outputNamedValue( '  Sorry, quantity not found in file', trim(qtyName) )
      return
    else
      call h5GOpen_f ( File%fileID%f_id, qtyName, grp_id, returnStatus )
    endif
    File%fileID%grp_id = grp_id
    call h5GOpen_f ( File%fileID%grp_id, chunkStr, grp_id, returnStatus )
    ! Begin the Reads
    ! Values
    ! call outputNamedValue ( 'rank', myRank )
    call GetHDF5DSRank ( grp_id, 'values', itsRank )
    if ( myRank > 0 .and. myRank /= itsRank ) &
      & call MLSL2Message ( MLSMSG_Warning, ModuleName // '%DirectRead', &
      & 'ranks differ for ' // trim(qtyName), &
      & MLSFile=File )
    select case ( myRank )
    case ( 1 )
      call LoadFromHDF5DS ( grp_id, 'values', quantity%value1 )
    case ( 2 )
      call LoadFromHDF5DS ( grp_id, 'values', quantity%values )
    case ( 3 )
      call LoadFromHDF5DS ( grp_id, 'values', quantity%value3 )
    case ( 4 )
      call LoadFromHDF5DS ( grp_id, 'values', quantity%value4 )
    case default
      call LoadFromHDF5DS ( grp_id, 'values', quantity%values )
    end select
    ! Geolocations
    if ( geolocations ) then
      call readQuantityAttributes ( quantity )
      call LoadFromHDF5DS( grp_id, 'surfs      ', quantity%template%surfs       )
      if ( associated(quantity%template%Geolocation) ) &
      & call LoadFromHDF5DS( grp_id, 'Geolocation  ', quantity%template%Geolocation )
      if ( allocated(quantity%template%Phi) ) &
      & call LoadFromHDF5DS( grp_id, 'Phi          ', quantity%template%Phi         )
      call LoadFromHDF5DS( grp_id, 'GeodLat        ', quantity%template%GeodLat     )
      call LoadFromHDF5DS( grp_id, 'Lon            ', quantity%template%Lon         )
      call LoadFromHDF5DS( grp_id, 'Time           ', quantity%template%Time        )
      call LoadFromHDF5DS( grp_id, 'SolarTime      ', quantity%template%SolarTime   )
      call LoadFromHDF5DS( grp_id, 'SolarZenith    ', quantity%template%SolarZenith )
      call LoadFromHDF5DS( grp_id, 'LosAngle       ', quantity%template%LosAngle    )
      if ( associated(quantity%template%CrossAngles) ) &
      & call LoadFromHDF5DS( grp_id, 'CrossAngles    ', quantity%template%CrossAngles   )
      if ( associated(quantity%template%Frequencies) ) &
      & call LoadFromHDF5DS( grp_id, 'Frequencies    ', quantity%template%Frequencies  )
      if ( associated(quantity%template%ChanInds) ) &
      & call LoadFromHDF5DS( grp_id, 'ChanInds       ', quantity%template%ChanInds   )
      if ( verbose ) then
        call dump( quantity%template%surfs        , 'surfs        ' )
        if ( associated(quantity%template%Geolocation) ) &
        & call dump( quantity%template%Geolocation, 'Geolocation  '  )
        if ( allocated(quantity%template%Phi) ) &
        & call dump( quantity%template%Phi, 'Phi          '          )
        call dump( quantity%template%GeodLat     , 'GeodLat        ' )
        call dump( quantity%template%Lon         , 'Lon            ' )
        call dump( &
        &  tai93s2hid( quantity%template%Time, leapsec=.true. ) &
        &                                        , 'Time (hid)     ' )
        call dump( quantity%template%SolarTime   , 'SolarTime      ' )
        call dump( quantity%template%SolarZenith , 'SolarZenith    ' )
        call dump( quantity%template%LosAngle    , 'LosAngle       ' )
        if ( associated(quantity%template%CrossAngles) ) &
        & call dump( quantity%template%CrossAngles  , 'CrossAngles    ' )
        if ( associated(quantity%template%Frequencies) ) &
        & call dump( quantity%template%Frequencies, 'Frequencies    '  )
        if ( associated(quantity%template%ChanInds) ) &
        & call dump( quantity%template%ChanInds  , 'ChanInds       '  )
      endif
    endif
    ! Close everything up
    call h5GClose_f ( grp_id, returnStatus )
    call h5GClose_f ( File%fileID%grp_id, returnStatus )

    call mls_CloseFile( File )
  contains
    subroutine readQuantityAttributes ( quantity )
      ! Args
      type (VectorValue_T), intent(inout) :: QUANTITY
      ! Local variables
      ! character(len=80) :: str
      ! Executable
      ! call get_string( lit_indices(quantity%template%quantityType ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'type', str )

      call GetHDF5Attribute ( File, 'NoSurfs', quantity%template%NoSurfs )
      call GetHDF5Attribute ( File, 'NoChans', quantity%template%NoChans )
      call GetHDF5Attribute ( File, 'NoCrossTrack', quantity%template%NoCrossTrack )
      call GetHDF5Attribute ( File, 'coherent', quantity%template%coherent )
      call GetHDF5Attribute ( File, 'stacked', quantity%template%stacked )
      call GetHDF5Attribute ( File, 'regular', quantity%template%regular )
      call GetHDF5Attribute ( File, 'minorFrame', quantity%template%minorFrame )
      call GetHDF5Attribute ( File, 'majorFrame', quantity%template%majorFrame )
      call GetHDF5Attribute ( File, 'logBasis', quantity%template%logBasis )
      call GetHDF5Attribute ( File, 'minValue', quantity%template%minValue )
      call GetHDF5Attribute ( File, 'badValue', quantity%template%badValue )

      ! call get_string( lit_indices(quantity%template%unit ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'unit', str )
      ! call get_string( lit_indices(quantity%template%verticalCoordinate ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'verticalCoordinate', str )
      ! call get_string( lit_indices(quantity%template%horizontalCoordinate ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'horizontalCoordinate', str )
      ! call get_string( lit_indices(quantity%template%latitudeCoordinate ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'latitudeCoordinate', str )
      ! call get_string( lit_indices(quantity%template%frequencyCoordinate ), str, strip=.true. )
      ! call GetHDF5Attribute ( File, 'frequencyCoordinate', str )
    end subroutine readQuantityAttributes
  end subroutine DirectRead_Quantity

  ! ------------------------------------------- DirectWriteVector_EveryQuantity --------
  subroutine DirectWriteVector_EveryQuantity ( File, Vector, &
    & chunkNo, options, rank )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    ! It dooes so an entrire vector's worth of vector quantities
    use MLSStrings, only: WriteIntsToChars
    use String_Table, only: Get_String
    ! Args:
    type (Vector_T), intent(in)   :: VECTOR
    type(MLSFile_T)               :: File
    integer, intent(in)           :: CHUNKNO      ! Index into chunks
    character(len=*), intent(in)  :: options
    integer, intent(in), optional :: rank
    ! Local parameters
    integer                       :: j
    logical                       :: nameQtyByTemplate
    type (VectorValue_T), pointer :: QUANTITY
    character(len=32)             :: SDNAME       ! Name of sd in output file
    logical :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', -1 )
    nameQtyByTemplate = .not. ( index(options, 'num') > 0 )
    if ( verbose ) then
      if ( vector%name > 0 ) then
        call get_string( vector%name, sdName )
      else
        sdname = '(unknown)'
      endif
      call outputNamedValue( 'DW L2AUX vector name', trim(sdName) )
    endif
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      if ( nameQtyByTemplate ) then
        call get_string( quantity%template%name, sdname )
      else
        call writeIntsToChars ( j, sdName )
        sdName = 'Quantity ' // trim(sdName)
      endif
      call DirectWrite_Quantity ( File, quantity, sdName, &
        & chunkNo, options, rank )
    enddo
  end subroutine DirectWriteVector_EveryQuantity

  ! ------------------------------------------ DirectWriteVector_L2GP --------
  subroutine DirectWriteVector_L2GP ( L2gpFile, &
    & vector, &
    & chunkNo, createSwath, lowerOverlap, upperOverlap, maxChunkSize )

    ! Purpose:
    ! Write out all the quantities in a vector as swaths to an hdfeos file
    ! Notes and limitations:
    ! Why don't you also supply an entire vector's worth of
    ! precision, quality, etc.?
    ! Args:
    type(MLSFile_T)               :: L2GPFile
    type (Vector_T), intent(in)   :: VECTOR
    integer, intent(in)              :: chunkNo
    logical, intent(in), optional :: createSwath
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    integer, intent(in), optional :: maxChunkSize
    ! Local variables
    integer                       :: j
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: precision
    type (VectorValue_T), pointer :: quality
    type (VectorValue_T), pointer :: status
    type (VectorValue_T), pointer :: Convergence
    type (VectorValue_T), pointer :: AscDescMode
    character(len=32)             :: SDNAME       ! Name of sd in output file
    ! Executable
    nullify(precision, quality, status, convergence, AscDescMode)
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      call get_string( quantity%template%name, sdname )
      call DirectWrite_L2GP ( L2gpFile, &
        & quantity, precision, quality, status, Convergence, AscDescMode, &
        & sdName, chunkNo, createSwath, lowerOverlap, upperOverlap, &
        & maxChunkSize )
    enddo
  end subroutine DirectWriteVector_L2GP

  ! ------------------------------------------ DirectWrite_L2GP --------
  subroutine DirectWrite_L2GP ( L2gpFile, &
    & quantity, precision, quality, status, Convergence, AscDescMode, &
    & sdName, chunkNo, createSwath, lowerOverlap, upperOverlap, &
    & maxChunkSize )

    ! Purpose:
    ! Write a swath to an hdfeos file composed of 
    ! what              supplied by
    ! l2gpvalues        quantity
    ! l2gpprecision     precision
    ! quality           quality
    ! status            status
    ! Convergence       Convergence
    ! AscDescMode       AscDescMode
    use HDF, only: Dfacc_Create, Dfacc_Rdonly, Dfacc_Rdwr
    use L2GPData, only: L2GPData_T, &
      & AppendL2GPData, DestroyL2GPContents, Dump
    use ReadApriori, only: WriteAPrioriAttributes
    ! Args:
    type(MLSFile_T)               :: L2GPFile
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: precision
    type (VectorValue_T), pointer :: quality
    type (VectorValue_T), pointer :: status
    type (VectorValue_T), pointer :: Convergence
    type (VectorValue_T), pointer :: AscDescMode
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    integer, intent(in)              :: chunkNo
    logical, intent(in), optional :: createSwath
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    integer, intent(in), optional :: maxChunkSize
    ! Local variables
    logical :: DeeBug
    logical :: alreadyThere
    integer :: FIRSTINSTANCE
    integer :: GRANDTOTALINSTANCES
    type (L2GPData_T) :: l2gp
    integer :: LASTINSTANCE
    integer :: OFFSET
    character(len=8) :: overlaps ! 'lower', 'upper', or 'none'
    integer :: NOTOWRITE
    integer :: TOTALPROFS
    logical :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', -1 )
    deebug = LetsDebug ( 'direct', 0 )
    if ( verbose ) call outputNamedValue ('DW ChunkNo', ChunkNo )
    if ( present(createSwath) ) then
      alreadyThere = mls_exists( L2GPFile%name ) == 0
      if ( .not. alreadyThere .and. verbose ) &
        & call outputNamedValue( '  creating file', trim(L2GPFile%name) )
      if ( alreadyThere ) &
        & alreadyThere = mls_swath_in_file ( L2GPFile%name, trim(sdname), &
        & L2GPFile%hdfVersion )
      if ( .not. verbose ) then
        ! Remain silent
        overlaps = 'none'
      elseif ( createSwath .and. alreadyThere ) then
        call outputNamedValue( '  recreating swath', trim(sdname) )
      elseif ( createSwath ) then
        call outputNamedValue( '  creating swath', trim(sdname) )
      elseif ( .not. alreadyThere ) then
        call outputNamedValue( '  How can we add to a nonexistent swath?', trim(sdname) )
      else
        call outputNamedValue( '  adding to swath', trim(sdname) )
      endif
    endif
    ! Size the problem
    overlaps = 'none'
    if ( present(lowerOverlap) ) then
      if ( lowerOverlap ) overlaps = 'lower'
    endif
    if ( present(upperOverlap) ) then
      if ( upperOverlap ) overlaps = 'upper'
    endif
    select case ( overlaps )
    case ( 'lower' )
      firstInstance =  1
      lastInstance = quantity%template%noInstancesLowerOverlap
      ! DeeBug = .true.
    case ( 'upper' )
      firstInstance =  quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap + 1
      lastInstance = quantity%template%noInstances
      ! DeeBug = .true.
    case ( 'none' )
      firstInstance = quantity%template%noInstancesLowerOverlap + 1
      lastInstance = quantity%template%noInstances - &
        & quantity%template%noInstancesUpperOverlap
    end select
    noToWrite = lastInstance - firstInstance + 1
    offset = quantity%template%instanceOffset
    grandtotalinstances = quantity%template%grandTotalInstances
    totalProfs = grandtotalinstances
    ! Check sanity
    if ( offset + noToWrite - 1 > grandTotalInstances .and. grandTotalInstances > 0 ) then
      call output('offset: ', advance='no')
      call output(offset, advance='no')
      call output('   noToWrite: ', advance='no')
      call output(noToWrite, advance='no')
      call output('   grandTotalInstances: ', advance='no')
      call output(grandTotalInstances, advance='yes')
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
        & 'last profile > grandTotalInstances for ' // trim(sdName), &
        & MLSFile=L2GPFile )
      ! Last-ditch effort resets offset
      offset = max(0, min(offset, grandTotalInstances - noToWrite))
    endif
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSL2Message( MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile )
    ! Convert vector quantity to l2gp
    call vectorValue_to_l2gp( quantity, &
      & precision, quality, status, convergence, AscDescMode, &
      & l2gp, &
      & sdname, chunkNo, offset=0, &
      & firstInstance=firstInstance, lastInstance=lastInstance )
    ! Output the l2gp into the file
    ! if ( l2gp%name == 'Temperature-InitPtan' .and. deebug ) then
    if ( deebug ) then
      call output('firstInstance: ', advance='no')
      call output(firstInstance, advance='no')
      call output('  lastInstance: ', advance='no')
      call output(lastInstance, advance='no')
      call output('  TotalProfs: ', advance='no')
      call output(TotalProfs, advance='yes')
      call dump(l2gp, Details=-1)
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
      call dump( l2gp%chunkNumber, 'Appending chunkNumber' )
    endif
    if ( verbose ) call outputNamedValue( 'DW L2GP qty name', trim(sdName) )
    ! call usleep ( delay ) ! Should we make this parallel%delay?
    call usleep ( parallel%delay ) ! Done!
    if ( present(createSwath) ) then
      print *, 'create swath? ', createSwath
    else
      print *, 'create swath not present '
    endif
    call AppendL2GPData( l2gp, l2gpFile, &
      & sdName, offset, lastprofile=lastInstance, &
      & TotNumProfs=TotalProfs, createSwath=createSwath, &
      & maxChunkSize=maxChunkSize )
    if ( L2GPFile%access == DFACC_CREATE ) L2GPFile%access = DFACC_RDWR
    call writeAPrioriAttributes( l2gpFile, dontreplace=.true. )
    if ( switchDetail(switches, 'l2gp') > -1 ) call dump(l2gp)
    ! Clear up our temporary l2gp
    call DestroyL2GPContents(l2gp)
  end subroutine DirectWrite_L2GP

  ! ------------------------------------------- DirectWriteVector_L2Aux --------
  subroutine DirectWriteVector_L2Aux ( L2AUXFile, Vector, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options, groupName )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    ! It dooes so an entrire vector's worth of vector quantities
    use Chunks_M, only: MLSChunk_T
    use ForwardModelConfig, only: ForwardModelConfig_T
    use MLSStrings, only: WriteIntsToChars
    ! Args:
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (Vector_T), intent(in)   :: VECTOR
    type(MLSFile_T)               :: L2AUXFile
    integer, intent(in)           :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write only the 1st instance
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: groupName
    ! Local parameters
    integer                       :: j
    logical                       :: nameQtyByTemplate
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=32)             :: SDNAME       ! Name of sd in output file
    logical                       :: useGroupName
    logical                       :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', 0 )
    nameQtyByTemplate = .true.
    if ( present(options) ) nameQtyByTemplate = &
      & .not. ( index(options, 'num') > 0 )
    useGroupName = .false.
    if ( present(options) ) useGroupName = &
      & ( index(options, 'g') > 0 ) .and. present( groupName )
    nullify(precision)
    if ( verbose ) then
      if ( vector%name > 0 ) then
        call get_string( vector%name, sdName )
      else
        sdname = '(unknown)'
      endif
      call outputNamedValue( 'DW L2AUX vector name', trim(sdName) )
    endif
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      if ( nameQtyByTemplate ) then
        call get_string( quantity%template%name, sdname )
      else
        call writeIntsToChars ( j, sdName )
        sdName = 'Quantity ' // trim(sdName)
      endif
      if ( useGroupName ) sdName = trim(groupName) // '/' // sdName
      call DirectWrite_L2Aux ( L2AUXFile, quantity, precision, sdName, &
        & chunkNo, chunks, FWModelConfig, &
        & lowerOverlap, upperOverlap, single, options )
    enddo
  end subroutine DirectWriteVector_L2Aux

  ! ------------------------------------------- DirectWriteVector_Hdf --------
  subroutine DirectWriteVector_Hdf ( L2AUXFile, Vector, &
    & options, groupName )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    ! It dooes so an entrire vector's worth of vector quantities
    use MLSStrings, only: WriteIntsToChars
    ! Args:
    type (Vector_T), intent(in)   :: VECTOR
    type(MLSFile_T)               :: L2AUXFile
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: groupName
    ! Local parameters
    integer                       :: j
    logical                       :: nameQtyByTemplate
    type (VectorValue_T), pointer :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=32)             :: SDNAME       ! Name of sd in output file
    logical                       :: useGroupName
    logical                       :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', 0 )
    nameQtyByTemplate = .true.
    if ( present(options) ) nameQtyByTemplate = &
      & .not. ( index(options, 'num') > 0 )
    useGroupName = .false.
    if ( present(options) ) useGroupName = &
      & ( index(options, 'g') > 0 ) .and. present( groupName )
    nullify(precision)
    if ( verbose ) then
      if ( vector%name > 0 ) then
        call get_string( vector%name, sdName )
      else
        sdname = '(unknown)'
      endif
      call outputNamedValue( 'DW L2AUX vector name', trim(sdName) )
    endif
    do j = 1, size(vector%quantities)
      quantity => vector%quantities(j)
      if ( nameQtyByTemplate ) then
        call get_string( quantity%template%name, sdname )
      else
        call writeIntsToChars ( j, sdName )
        sdName = 'Quantity ' // trim(sdName)
      endif
      if ( useGroupName ) sdName = trim(groupName) // '/' // sdName
      call DirectWrite_Hdf ( L2AUXFile, quantity, precision, sdName, &
        & options )
    enddo
  end subroutine DirectWriteVector_Hdf

  ! ------------------------------------------- DirectWrite_Hdf --------
  subroutine DirectWrite_Hdf ( L2AUXFile, quantity, precision, sdName, &
    & options )

    ! Purpose:
    ! Write plain hdf-formatted files
    use MLSHDF5, only: SaveAsHDF5DS

    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    character(len=*), intent(in), optional :: options
    ! Local parameters
    logical :: alreadyOpen
    logical :: already_there
    ! logical, parameter :: DEEBUG = .false.
    ! logical :: deebughere
    integer :: returnStatus
    ! character(len=*), parameter :: sdDebug = "R1A:118.B1F:PT.S0.FB25-1 Core"
    logical :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', 0 )
    alreadyOpen = L2AUXFile%stillOpen
    already_There = mls_exists( L2AUXFile%name ) == 0
    if ( .not. already_There ) then
      call outputNamedValue( '  creating file', trim(L2AUXFile%name) )
    elseif ( verbose ) then
      call outputNamedValue( '  no need to recreate file', trim(L2AUXFile%name) )
      call outputNamedValue( '  already open?', alreadyOpen )
    endif
    ! deebughere = ( deebug .or. sdname == sdDebug ) .and. .false.
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2AUXFile, returnStatus)
      if ( returnStatus /= 0 ) &
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to open l2aux file', MLSFile=L2AUXFile )
    endif
    call SaveAsHDF5DS ( L2AUXFile%FileID%f_id, SDname, quantity%values )
    if ( .not. alreadyOpen )  call mls_closeFile(L2AUXFile, returnStatus)
    L2AUXFile%errorCode = returnStatus
    L2AUXFile%lastOperation = 'write'
    if ( switchDetail(switches, 'l2aux') < 0 ) return
    call dump(quantity)
    if ( associated(precision) ) call dump(precision)
  end subroutine DirectWrite_Hdf

  ! ------------------------------------------- DirectWrite_L2Aux --------
  subroutine DirectWrite_L2Aux ( L2AUXFile, quantity, precision, sdName, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options )

    ! Purpose:
    ! Write plain hdf-formatted files ala l2aux for datasets that
    ! are too big to keep all chunks stored in memory
    ! so instead write them out chunk-by-chunk
    
    ! Despite the name the routine takes vector quantities, not l2aux ones
    use Chunks_M, only: MLSChunk_T
    use ChunkDivide_M, only: ChunkDivideConfig
    use ForwardModelConfig, only: ForwardModelConfig_T
    use HDF, only: Dfacc_Rdonly

    ! Args:
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    type (VectorValue_T), intent(in) :: QUANTITY
    type (VectorValue_T), pointer :: PRECISION
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write only the 1st instance
    character(len=*), intent(in), optional :: options
    ! Local parameters
    logical :: alreadyOpen
    logical :: already_there
    logical, parameter :: DEEBUG = .false.
    logical :: deebughere
    integer :: lastMAF
    integer :: returnStatus
    character(len=*), parameter :: sdDebug = "R1A:118.B1F:PT.S0.FB25-1 Core"
    logical :: verbose
    ! Executable
    verbose = BeVerbose ( 'direct', 0 )
    alreadyOpen = L2AUXFile%stillOpen
    already_There = mls_exists( L2AUXFile%name ) == 0
    if ( .not. already_There ) then
      call outputNamedValue( '  creating file', trim(L2AUXFile%name) )
    elseif ( verbose ) then
      call outputNamedValue( '  no need to recreate file', trim(L2AUXFile%name) )
      call outputNamedValue( '  already open?', alreadyOpen )
    endif
    deebughere = ( deebug .or. sdname == sdDebug ) .and. .false.
    lastMAF = (quantity%template%instanceOffset+quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap - &
        & quantity%template%noInstancesUpperOverlap)
    if ( chunkNo == size(chunks) .and. MAYWRITEPOSTOVERLAPS .and. &
      & ChunkDivideConfig%allowPostOverlaps ) &
      & lastMAF = quantity%template%instanceOffset + &
      & quantity%template%noInstances - quantity%template%noInstancesLowerOverlap
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2AUXFile, returnStatus)
      if ( returnStatus /= 0 ) &
        call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to open l2aux file', MLSFile=L2AUXFile )
    endif
    if ( L2AUXFile%access == DFACC_RDONLY )  &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'l2aux file is rdonly', MLSFile=L2AUXFile )
    if ( (lastProfTooBigWarns <= MAXNUMWARNS) &
      & .and. &
      & (quantity%template%instanceOffset+quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap) &
      & > &
      & quantity%template%grandTotalInstances &
      & .and. &
      & quantity%template%grandTotalInstances > 0 ) then
      lastProfTooBigWarns = lastProfTooBigWarns + 1
      call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'last profile > grandTotalInstances for ' // trim(sdName), &
          & MLSFile=L2AUXFile)
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('last instance: ', advance='no')
      call output(lastMAF, advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
      if ( lastProfTooBigWarns > MAXNUMWARNS ) &
          & call MLSL2Message ( MLSMSG_Warning, ModuleName, &
          & 'Max no. of warnings reached--suppressing further ones')
    endif
    if ( deebughere ) then
      print *, 'Direct Writing to ', trim(L2AUXFile%name)
      print *, 'hdfVersion ', L2AUXFile%hdfVersion
      call output('instanceOffset: ', advance='no')
      call output(quantity%template%instanceOffset, advance='yes')
      call output('noInstances: ', advance='no')
      call output(quantity%template%noInstances, advance='yes')
      call output('noInstancesLowerOverlap: ', advance='no')
      call output(quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('noInstancesUpperOverlap: ', advance='no')
      call output(quantity%template%noInstancesUpperOverlap, advance='yes')
      call output('last profile: ', advance='no')
      call output(lastMAF, advance='yes')
      call output('grandTotalInstances: ', advance='no')
      call output(quantity%template%grandTotalInstances, advance='yes')
    endif
    select case (l2AUXFile%hdfversion)
    case (HDFVERSION_4)
      call DirectWrite_L2Aux_hdf4 ( quantity, sdName, L2AUXFile, &
        & chunkNo, chunks, &
        & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_hdf4 ( precision, &
        & trim(sdName) // 'precision', L2AUXFile, chunkNo, chunks, &
        & lowerOverlap=lowerOverlap, upperOverlap=upperOverlap )
    case (HDFVERSION_5)
      call DirectWrite_L2Aux_hdf5 ( quantity, sdName, L2AUXFile, &
        & chunkNo, chunks, &
        & FWModelConfig, lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
        & single=single, options=options )
      if ( associated(precision) ) & 
        & call DirectWrite_L2Aux_hdf5 ( precision, &
        & trim(sdName) // 'precision', L2AUXFile, chunkNo, chunks, &
        & FWModelConfig, lowerOverlap=lowerOverlap, upperOverlap=upperOverlap, &
        & single=single, options=options )
    case default
      call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unsupported hdfVersion for DirectWrite_L2Aux (currently only 4 or 5)' )
    end select
    if ( verbose ) call outputNamedValue( 'DW L2AUX qty name', trim(sdName) )
    if ( .not. alreadyOpen )  call mls_closeFile(L2AUXFile, returnStatus)
    L2AUXFile%errorCode = returnStatus
    L2AUXFile%lastOperation = 'write'
    if ( switchDetail(switches, 'l2aux') < 0 ) return
    call dump(quantity)
    if ( associated(precision) ) call dump(precision)
  end subroutine DirectWrite_L2Aux

  ! ------------------------------------------ DirectWrite_L2Aux_hdf4 --------
  subroutine DirectWrite_L2Aux_hdf4 ( quantity, sdName, L2AUXFile, &
    & chunkNo, chunks, lowerOverlap, upperOverlap )

    use Chunks_M, only: MLSChunk_T
    use HDF, only: Sfn2index, Sfselect, Sfcreate, &
      & Sfendacc, Dfnt_Float32, SfwData_F90
    use Intrinsic, only: L_None
    use MLSKinds, only: R4, R8

    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    ! integer, intent(in) :: SDNAME       ! Name of sd in output file
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap

    ! Local variables
    integer :: SDINDEX                  ! Index of sd
    integer :: SDID                     ! Handle for sd
    integer :: STATUS                   ! Status flag
    integer :: START(3)                 ! HDF array starting position
    integer :: STRIDE(3)                ! HDF array stride
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: NODIMS                   ! Also index of maf dimension
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    real (r8) :: HUGER4

    ! executable code

    hugeR4 = real ( huge(0.0_r4), r8 )

    if ( quantity%template%frequencyCoordinate == L_None ) then
      noDims = 2
    else
      noDims = 3
    end if

    ! Create or access the SD
    sdIndex = sfn2index ( L2AUXFile%fileID%f_id, trim(sdName) )
    if ( sdIndex == -1 ) then
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
      sdId  = sfCreate ( L2AUXFile%fileID%f_id, trim(sdName), DFNT_FLOAT32, &
        & noDims, sizes )
    else
      sdId = sfSelect ( L2AUXFile%fileID%f_id, sdIndex )
    end if
    if ( sdId == -1 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
     & 'Error accessing SD '//trim(sdName) // ' (hdf4)', MLSFile=L2AUXFile)

    ! What exactly will be our contribution
    stride = 1
    start = 0
    sizes = 1
    sizes(noDims) = quantity%template%noInstances - &
      & quantity%template%noInstancesLowerOverlap - &
      & quantity%template%noInstancesUpperOverlap
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    start(noDims) = quantity%template%instanceOffset

    ! Now write it out
    status = 0
    if ( sizes(noDims) > 0 ) &
      & status = SFWDATA_F90(sdId, start(1:noDims), &
      & stride(1:noDims), sizes(1:noDims), &
      & real ( max ( -hugeR4, min ( hugeR4, &
      &   quantity%values ( :, &
      &   1+quantity%template%noInstancesLowerOverlap : &
      &    quantity%template%noInstances - quantity%template%noInstancesUpperOverlap &
      &  ) ) ) ) )
    if ( status /= 0 ) then
      call announce_error (0,&
        & "Error writing SDS data " // trim(sdName) // " to l2aux file:  " )
    end if
    if ( DEEBUG ) then
      call output('noDims: ', advance='no')
      call output(noDims, advance='yes')
      call output('start: ', advance='no')
      call output(start, advance='yes')
      call output('stride: ', advance='no')
      call output(stride, advance='yes')
      call output('sizes: ', advance='no')
      call output(sizes, advance='yes')
      call output('shape(values): ', advance='no')
      call output(shape(quantity%values), advance='yes')
      call output('first instance: ', advance='no')
      call output(1+quantity%template%noInstancesLowerOverlap, advance='yes')
      call output('last instance: ', advance='no')
      call output(quantity%template%noInstances &
        & - quantity%template%noInstancesUpperOverlap, advance='yes')
    end if

    ! End access to the SD and close the file
    status = sfEndAcc ( sdId )
    if ( status == -1 ) call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Error ending access to direct write sd (hdf4)', MLSFile=L2AUXFile )

  end subroutine DirectWrite_L2Aux_hdf4

  ! ------------------------------------------ DirectWrite_L2Aux_hdf5 --------
  subroutine DirectWrite_L2Aux_hdf5 ( quantity, sdName, L2AUXFile, &
    & chunkNo, chunks, FWModelConfig, &
    & lowerOverlap, upperOverlap, single, options )

    use Chunks_M, only: MLSChunk_T
    use ChunkDivide_M, only: ChunkDivideConfig
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelSupport, only: Showfwdmodelnames
    use HDF5, only: H5gclose_F, H5gopen_F
    use Intrinsic, only: L_None
    use L2AuxData, only: L2AuxData_T, PhaseNameAttributes, &
      & DestroyL2Auxcontents, &
      & SetupnewL2Auxrecord, WriteL2Auxattributes, WriteHDF5Data
    use MLSHDF5, only: IsHDF5attributepresent, IsHDF5dspresent, &
      & MakeHDF5Attribute
    use MLSL2Timings, only: ShowTimingNames
    use PCFHdr, only: H5_WriteMLSFileattr, H5_Writeglobalattr
    use QuantityTemplates, only: WriteAttributes

    ! Args:
    type (VectorValue_T), intent(in) :: QUANTITY
    character(len=*), intent(in) :: SDNAME       ! Name of sd in output file
    type(MLSFile_T)                :: L2AUXFile
    integer, intent(in) :: CHUNKNO      ! Index into chunks
    type (MLSChunk_T), dimension(:), intent(in) :: CHUNKS
    type(ForwardModelConfig_T), dimension(:), pointer :: FWModelConfig
    logical, intent(in), optional :: lowerOverlap
    logical, intent(in), optional :: upperOverlap
    logical, intent(in), optional :: single       ! Write 1st instance only
    character(len=*), intent(in), optional :: options

    ! Local variables
    logical :: addQtyAttributes
    logical :: already_there
    ! logical :: attributes_there
    ! character(len=128) :: barename
    integer :: first_maf
    ! character(len=128), dimension(25) :: groupNames
    integer :: grp_id
    type (L2AUXData_T) :: l2aux
    integer :: last_maf
    type ( MLSChunk_T ) :: LASTCHUNK    ! The last chunk in the file
    logical :: mySingle
    ! integer :: n
    integer :: NODIMS                   ! Also index of maf dimension
    integer :: Num_qty_values
    character(len=8) :: overlaps        ! 'lower', 'upper', or 'none'
    ! character(len=1024) :: path
    integer :: returnStatus
    integer :: SIZES(3)                 ! HDF array sizes
    integer :: START(3)                 ! HDF array starting position
    integer :: total_DS_size
    logical, parameter :: MAYCOLLAPSEDIMS = .false.
    ! logical, parameter :: DEEBUG = .true.
    logical :: verbose

    ! executable code
    verbose = BeVerbose ( 'direct', 0 )
    Num_qty_values = size(quantity%values, 1)*size(quantity%values, 2)
    addQtyAttributes = .false.
    if ( present(options) ) addQtyAttributes = ( index(options, 'A') > 0 )
    if ( quantity%template%frequencyCoordinate == L_None &
      & .and. MAYCOLLAPSEDIMS) then
      noDims = 2
    else
      noDims = 3
    end if
    
    overlaps = 'none'
    if ( present(lowerOverlap) ) then
      if ( lowerOverlap ) overlaps = 'lower'
    endif
    if ( present(upperOverlap) ) then
      if ( upperOverlap ) overlaps = 'upper'
    endif
    
    mySingle = .false.
    if ( present(single) ) mySingle = single

    ! Create or access the SD
    already_There = mls_exists( L2AUXFile%name ) == 0
    if ( .not. already_There ) &
      & call outputNamedValue( '  creating file', trim(L2AUXFile%name) )
    already_there = IsHDF5DSPresent(L2AUXFile%fileID%f_id, trim(sdName))
    if ( .not. already_there ) then
      call outputNamedValue( '  creating sd', trim(sdname) )
      lastChunk = chunks(size(chunks))
      sizes(noDims) = lastChunk%lastMAFIndex - lastChunk%noMAFSUpperOverlap + 1
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps ) &
        & sizes(noDims) = lastChunk%lastMAFIndex + 1
      sizes(noDims-1) = quantity%template%noSurfs
      if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    elseif ( verbose )then
      call outputNamedValue( '  adding to sd', trim(sdname) )
    end if

    ! What exactly will be our contribution
    start = 0
    sizes = 1
    sizes(noDims-1) = quantity%template%noSurfs
    if ( noDims == 3 ) sizes(1) = quantity%template%noChans
    start(noDims) = quantity%template%instanceOffset
    select case ( overlaps )
    case ( 'lower' )
      sizes(noDims) = quantity%template%noInstancesLowerOverlap
      first_maf = 1
      last_maf = quantity%template%noInstances
    case ( 'upper' )
      sizes(noDims) = quantity%template%noInstancesUpperOverlap
      first_maf = quantity%template%noInstances &
        &       - quantity%template%noInstancesUpperOverlap + 1
      last_maf = quantity%template%noInstances
    case ( 'none' )
      sizes(noDims) = quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap - &
        & quantity%template%noInstancesUpperOverlap
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps .and. &
        & chunkNo == size(chunks) ) &
        & sizes(noDims) = quantity%template%noInstances - &
        & quantity%template%noInstancesLowerOverlap
      first_maf = 1+quantity%template%noInstancesLowerOverlap
      last_maf = quantity%template%noInstances &
        &       - quantity%template%noInstancesUpperOverlap
      if ( MAYWRITEPOSTOVERLAPS .and. ChunkDivideConfig%allowPostOverlaps .and. &
        & chunkNo == size(chunks) ) &
        & last_maf = quantity%template%noInstances
    end select
    
    if ( mySingle ) then
      ! In case we write only the 1st instance
      last_maf = first_maf
      sizes(noDims) = 1
    endif

    if ( DEEBUG ) then
      print *, 'sdname ', trim(sdName)
      print *, 'already_there ', already_there
      print *, 'noDims ', noDims
      print *, 'instanceOffset ', quantity%template%instanceOffset
      print *, 'minorFrame ', quantity%template%minorFrame
      print *, 'start ', start
      print *, 'sizes ', sizes
      print *, 'shape(quantity%values) ', shape(quantity%values)
      print *, 'first_maf ', first_maf
      print *, 'last_maf ', last_maf
    endif
    ! Make certain things will fit
    if ( noDims == 3 ) then
      total_DS_size = sizes(1)*sizes(2)*sizes(3)
      if ( DEEBUG ) then
        print *, 'total_DS_size ', total_DS_size
        print *, 'Num_qty_values ', Num_qty_values
      endif
      if ( total_DS_size > Num_qty_values ) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Number of 3d array elements to write > number stored in qty values', &
        & MLSFile=L2AUXFile )
      call WriteHDF5Data (  real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:3)) &
        & ), L2AUXFile%FileID%f_id, &
        & returnStatus, trim(sdname), already_there, &
        & start, sizes )
    else
      total_DS_size = sizes(1)*sizes(2)
      if ( DEEBUG ) then
        print *, 'total_DS_size ', total_DS_size
        print *, 'Num_qty_values ', Num_qty_values
      endif
      if ( total_DS_size > Num_qty_values ) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Number of 2d array elements to write > number stored in qty values', &
        & MLSFile=L2AUXFile )
      call WriteHDF5Data (  real( &
        &   reshape(quantity%values(:,first_maf:last_maf), sizes(1:2)) &
        & ), L2AUXFile%FileID%f_id, &
        & returnStatus, trim(sdname), already_there, &
        & start(1:2), sizes(1:2) )
    endif

    ! call mls_CloseFile(L2AUXFILE)
    ! attributes_there = IsHDF5AttributeInFile( L2AUXFile%name, 'Phase Names' )
    ! call mls_OpenFile(L2AUXFILE)
    ! attributes_there = .false.
    ! Now some attribute stuff
    if ( PHASENAMEATTRIBUTES ) then
      if ( DeeBUG ) call dump(L2AUXFile, details=1)
      call h5gopen_f(L2AUXFile%fileID%f_id, '/', grp_id, returnstatus)
      if ( returnstatus /= 0 ) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to open group "/" to write attribute', &
        & MLSFile=L2AUXFile )
      if ( associated ( FWModelConfig ) ) then
        call MakeHDF5Attribute(grp_id, &
          & 'ForwardModel Names', ShowFwdModelNames(FWModelConfig), &
          & skip_if_already_there=.false., dont_trim=.true.)
      end if
      call MakeHDF5Attribute(grp_id, &
        & 'Phase Names', showTimingNames('phases', .true.), &
        & skip_if_already_there=.false., dont_trim=.true.)
      call MakeHDF5Attribute(grp_id, &
       & 'MiscNotes', GlobalAttributes%MiscNotes, .true.)
      call h5gclose_f(grp_id, returnstatus)
    endif

    if ( already_there ) return
    ! sd-level attributes
    call SetupNewL2AUXRecord ( L2AUX, quantity%template, &
      & first_MAF, last_MAF-first_MAF+1 )
    if ( DEEBUG ) then
      call output('Writing attributes to: ', advance='no')
      call output(trim(sdName), advance='yes')
    endif
    call WriteL2AUXAttributes(L2AUXFile%fileID%f_id, l2aux, trim(sdName))
    if ( addQtyAttributes ) &
      & call WriteAttributes ( L2AUXFile%fileID%f_id, trim(sdName), quantity%template )
    ! Deallocate memory used by the l2aux
    call DestroyL2AUXContents ( l2aux )
    ! file-level attributes
    call h5_writeglobalattr(L2AUXFile%fileID%f_id, skip_if_already_there=.true.)
    if ( PHASENAMEATTRIBUTES ) then
      call h5gopen_f(L2AUXFile%fileID%f_id, '/', grp_id, returnstatus)
      if ( .not. &
        & IsHDF5AttributePresent('/', L2AUXFile%fileID%f_id, 'Section Names') ) &
        & call MakeHDF5Attribute(grp_id, &
        & 'Section Names', trim(showTimingNames('sections', .true.)), .true.)
      call h5gclose_f(grp_id, returnstatus)
    endif
    if ( WRITEFILEATTRIBUTES ) call h5_writeMLSFileAttr(L2AUXFile, skip_if_already_there=.true.)

  end subroutine DirectWrite_L2Aux_hdf5

  ! ------------------------------------------ DirectWrite_Quantity --------
  ! We write the quantity to a file, chunk by chunk, overlaps and all
  ! E.g., the qtyname is 'Q', and we have written chunks 1 and 2. Then the
  ! file layout will be
  ! /
  !   Q/
  !     1/
  !       values
  !       lons
  !       lats
  !         ...
  !     2/
  !       values
  !       lons
  !       lats
  !         ...
  
  ! Alternative use: we may instead the contents of an input file as if it
  ! contained the quantity values
  !
  subroutine DirectWrite_Quantity ( File, Quantity, QtyName, &
    & ChunkNo, Options, Rank, inputFile )

    use HDF5, only: H5gclose_F, H5gcreate_F, H5gopen_F
    use Intrinsic, only: Lit_Indices
    use MLSHDF5, only: IsHDF5grouppresent, &
      & MakeHDF5attribute, SaveasHDF5ds
    use MLSL2Options, only: MLSL2Message
    use MLSSTrings, only: Writeintstochars
    use MoreMessage, only: MoreMLSMessage => MLSMessage
    ! Args:
    type(MLSFile_T)                  :: File
    type (VectorValue_T), intent(in) :: Quantity
    character(len=*), intent(in)     :: QtyName  ! Name of qty in output file
    integer, intent(in)              :: ChunkNo  ! Index into chunks
    character(len=*), intent(in), optional :: Options
    integer, intent(in), optional    :: Rank
    character(len=*), intent(in), optional :: inputFile

    ! Local variables
    logical :: already_there
    ! logical :: attributes_there
    integer :: grp_id
    integer :: myRank
    character(len=10) :: chunkStr ! chunk number '1', '2', ..,
    integer :: returnStatus
    integer :: TrueRank ! of the quantity's value
    ! logical, parameter :: DEEBUG = .true.
    logical :: verbose

    ! executable code
    verbose = BeVerbose ( 'direct', -1 )
    myRank = 2
    if ( present(rank) ) then
      if ( rank < 5 ) myRank = rank
      ! call outputNamedValue ( 'input rank', rank )
    end if
    call writeIntsToChars ( chunkNo, chunkStr )
    chunkStr = adjustl ( chunkStr )
    ! Create or access the SD
    already_There = mls_exists( File%name ) == 0
    if ( .not. already_There ) then
      call outputNamedValue( '  creating file', trim(File%name) )
    end if
    call mls_openFile ( File, returnStatus )
    already_there = IsHDF5GroupPresent( File%fileID%f_id, trim(qtyName) )
    if ( .not. already_There ) then
      if ( verbose ) call outputNamedValue( '  creating file', trim(File%name) )
      call h5GCreate_f ( File%fileID%f_id, qtyName, grp_id, returnStatus )
      if ( returnStatus /= 0 ) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to create group ' // trim(qtyName), MLSFile=File )
      call writeQuantityAttributes ( grp_id, quantity )
    else
      call h5GOpen_f ( File%fileID%f_id, qtyName, grp_id, returnStatus )
      if ( returnStatus /= 0 ) &
        & call MLSL2Message ( MLSMSG_Error, ModuleName, &
        & 'Unable to open group ' // trim(qtyName), MLSFile=File )
    end if
    File%fileID%grp_id = grp_id
    call h5GCreate_f ( File%fileID%grp_id, chunkStr, grp_id, returnStatus )
    if ( returnStatus /= 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Unable to create group ' // trim(chunkStr), MLSFile=File )
    ! Begin the writes
    ! Values
    ! Are we getting the values from an input file?
    if ( present(inputFile) ) then
      if ( verbose ) call outputNamedValue( 'inputFile', trim(inputFile) )
      call SaveAsHDF5DS ( inputFile, grp_id, &
        & 'values', maxLineLen=4096, fromNull='@' )
    else
    
      ! call outputNamedValue ( 'rank', myRank )
      trueRank = 2 + merge(1,0,quantity%template%noChans>1) + &
               &     merge(1,0,quantity%template%noCrossTrack>1)
      if ( myRank <= 0 ) then
        myRank = trueRank
      else if ( trueRank > myRank ) then
        call MoreMLSMessage ( MLSMSG_Warning, moduleName, "Actual rank of Value %S is " // &
          & " %D but the specified rank is %D.  Are you sure this is what you want?", &
          & datum=[ quantity%template%name, trueRank, myRank] )
      end if
      select case ( myRank )
      case ( 1 )
        call SaveAsHDF5DS( grp_id, 'values', quantity%value1 )
      case ( 2 )
        call SaveAsHDF5DS( grp_id, 'values', quantity%values )
      case ( 3 )
        call SaveAsHDF5DS( grp_id, 'values', quantity%value3 )
      case ( 4 )
        call SaveAsHDF5DS( grp_id, 'values', quantity%value4 )
      case default
        call SaveAsHDF5DS( grp_id, 'values', quantity%values )
      end select
    endif
    ! Geolocations
    call SaveAsHDF5DS( grp_id, 'surfs', quantity%template%surfs )
    if ( associated(quantity%template%Geolocation) ) &
  & call SaveAsHDF5DS( grp_id, 'Geolocation', quantity%template%Geolocation )
    if ( allocated(quantity%template%Phi) ) &
  & call SaveAsHDF5DS( grp_id, 'Phi        ', quantity%template%Phi         )
    call SaveAsHDF5DS( grp_id, 'GeodLat    ', quantity%template%GeodLat     )
    call SaveAsHDF5DS( grp_id, 'Lon        ', quantity%template%Lon         )
    call SaveAsHDF5DS( grp_id, 'Time       ', quantity%template%Time        )
    call SaveAsHDF5DS( grp_id, 'SolarTime  ', quantity%template%SolarTime   )
    call SaveAsHDF5DS( grp_id, 'SolarZenith', quantity%template%SolarZenith )
    call SaveAsHDF5DS( grp_id, 'LosAngle   ', quantity%template%LosAngle    )
    if ( associated(quantity%template%CrossAngles) ) &
  & call SaveAsHDF5DS( grp_id, 'CrossAngles', quantity%template%CrossAngles )
    if ( associated(quantity%template%Frequencies) ) &
  & call SaveAsHDF5DS( grp_id, 'Frequencies', quantity%template%Frequencies )
    if ( associated(quantity%template%ChanInds) ) &
  & call SaveAsHDF5DS( grp_id, 'ChanInds   ', quantity%template%ChanInds    )
    
    ! Close everything up
    call h5GClose_f ( grp_id, returnStatus )
    if ( returnStatus /= 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Unable to close chunk group ' // trim(chunkStr), MLSFile=File )
    call h5GClose_f ( File%fileID%grp_id, returnStatus )
    if ( returnStatus /= 0 ) &
      & call MLSL2Message ( MLSMSG_Error, ModuleName, &
      & 'Unable to close qty group ' // trim(qtyName), MLSFile=File )

    call mls_CloseFile( File )
  contains
    subroutine writeQuantityAttributes ( locID, quantity )
      ! Args
      integer, intent(in) :: locID
      type (VectorValue_T), intent(in) :: QUANTITY
      ! Local variables
      character(len=80) :: str
      ! Executable
      call get_string( lit_indices(quantity%template%quantityType ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'type', str )

      call MakeHDF5Attribute ( locID, 'NoSurfs', quantity%template%NoSurfs )
      call MakeHDF5Attribute ( locID, 'NoChans', quantity%template%NoChans )
      call MakeHDF5Attribute ( locID, 'NoCrossTrack', quantity%template%NoCrossTrack )
      call MakeHDF5Attribute ( locID, 'coherent', quantity%template%coherent )
      call MakeHDF5Attribute ( locID, 'stacked', quantity%template%stacked )
      call MakeHDF5Attribute ( locID, 'regular', quantity%template%regular )
      call MakeHDF5Attribute ( locID, 'minorFrame', quantity%template%minorFrame )
      call MakeHDF5Attribute ( locID, 'majorFrame', quantity%template%majorFrame )
      call MakeHDF5Attribute ( locID, 'logBasis', quantity%template%logBasis )
      call MakeHDF5Attribute ( locID, 'minValue', quantity%template%minValue )
      call MakeHDF5Attribute ( locID, 'badValue', quantity%template%badValue )

      call get_string( lit_indices(quantity%template%unit ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'unit', str )
      call get_string( lit_indices(quantity%template%verticalCoordinate ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'verticalCoordinate', str )
      call get_string( lit_indices(quantity%template%horizontalCoordinate ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'horizontalCoordinate', str )
      call get_string( lit_indices(quantity%template%latitudeCoordinate ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'latitudeCoordinate', str )
      call get_string( lit_indices(quantity%template%frequencyCoordinate ), str, strip=.true. )
      call MakeHDF5Attribute ( locID, 'frequencyCoordinate', str )
    end subroutine writeQuantityAttributes
  end subroutine DirectWrite_Quantity

  !------------------------------------------  DumpDirectDB  -----
  subroutine DumpDirectDB ( directDB, Details )

    ! This routine dumps the DirectWrite DB

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer  :: directDB
    integer, intent(in), optional :: DETAILS

    ! Local variables
    character(len=256), dimension(:,:), pointer :: array
    integer :: i
    integer :: myDetails
    integer :: n
    ! Executable
    myDetails = 1
    if ( present(details) ) myDetails = details
    call output ( '========== DirectWrite Data Base ==========', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( .not. associated(directDB) ) then
      call output ( '**** directWrite Database empty ****', advance='yes' )
      return
    elseif ( size(directDB) < 1 ) then
      call output ( '**** directWrite Database empty ****', advance='yes' )
      return
    endif
    if ( myDetails < -1 ) then
      nullify( array )
      n = size(directDB)
      allocate( array(n+1, 2 ) )
      array(1,1) = 'name'
      array(1,2) = 'path'
      do i = 1, size(directDB)
        call split_path_name( directDB(i)%filename, array(i+1, 2), array(i+1, 1) )
      end do
      call outputTable( array, border='-', headliner='-' )
      deallocate( array )
    else
      do i = 1, size(directDB)
        call DumpDirectWrite(directDB(i), Details)
      end do
    endif
  end subroutine DumpDirectDB

  !------------------------------------------  DumpDirectWrite  -----
  subroutine DumpDirectWrite ( directWrite, Details )

    ! This routine dumps the DirectWrite DB

    ! Dummy arguments
    type (DirectData_T)  :: directWrite
    integer, intent(in), optional :: DETAILS

    ! Local variables
    integer :: i
    integer :: myDetails
    ! Executable
    myDetails = 1
    if ( present(details) ) myDetails = details

    if ( myDetails > -2 ) then
      call output ( 'File Name index: ')
      call output ( directWrite%fileIndex )
      call output ( '   File Name (base): ')
      call output ( trim(directWrite%fileNameBase), advance='yes' )
    endif
    call output ( 'Full File Name  : ')
    call output ( trim(directWrite%fileName), advance='yes' )
    if ( myDetails < -1 ) return
    call output ( 'Type  : ')
    call output ( directWrite%type, advance='no' )
    call blanks ( 3 )
    if ( directWrite%type == l_l2aux ) then
      call output ( '(l2aux)', advance='yes')
    elseif ( directWrite%type == l_l2gp ) then
      call output ( '(l2gp)', advance='yes')
    elseif ( directWrite%type == l_l2dgg ) then
      call output ( '(l2dgg)', advance='yes')
    elseif ( directWrite%type == l_l2fwm ) then
      call output ( '(l2fwm)', advance='yes')
    else
      call output ( '(unknown)', advance='yes')
    endif
    if ( directWrite%autoType < 1 ) then
      call output ( '(Is not', advance='no')
    else
      call output ( '(Is', advance='no')
    endif
    call output (' eligible to be auto-filled)', advance='yes' )
    call output ( 'Handle  : ')
    call output ( directWrite%Handle, advance='yes' )
    call output ( 'Num sci. datasets  : ')
    call output ( directWrite%NSDNames, advance='yes' )
    if ( .not. associated(directWrite%SDNames)  .or. myDetails < 1 ) then
      call output ( 's.d. names not associated', advance='yes' )
      return
    endif
    do i = 1, size(directWrite%sdNames)
      call blanks ( 3 )
      call output(trim(directWrite%SDNames(i)), advance='yes')
    end do
  end subroutine DumpDirectWrite

  !------------------------------------------  ExpandDirectDB  -----
  subroutine ExpandDirectDB ( directDB, fileNameBase, directData, isNew )

    ! This routine expands the DirectWrite DB if necessary
    ! Returning either a match or a new one

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer  :: directDB
    type (DirectData_T), pointer :: directData
    character(len=*), intent(in) :: fileNameBase
    logical, intent(out) :: isNew

    ! Local variables
    integer :: dbID
    type (DirectData_T)                      :: tempDirectData
    ! Begin executable
    isNew = .true.
    ! Do we have any database yet?
    if ( size(directDB) < 1 ) then
      ! print *, 'Setting up initial database with ', trim(fileNameBase)
      call SetupNewDirect(tempDirectData, 0)
      dbID = AddDirectToDatabase( directDB, tempDirectData )
      directData => directDB(dbID)
    else
      ! Check if the fileNameBase already there
      ! print *, 'Checking if ', trim(fileNameBase), ' is in database'
      do dbID = 1, size(directDB)
        ! print *, dbID, trim(directDB(dbID)%fileNameBase)
      enddo
      dbID = FindFirst(fileNameBase == directDB%fileNameBase)
      if ( dbID > 0 ) then
        directData => directDB(dbID)
        isNew = .false.
        ! print *, 'Recognized ', trim(fileNameBase), ' in ', trim(directData%FileName)
      else
        call SetupNewDirect(tempDirectData, 0)
        dbID = AddDirectToDatabase( directDB, tempDirectData )
        directData => directDB(dbID)
        ! print *, 'Adding to database ', trim(fileNameBase)
      endif
    endif
  end subroutine ExpandDirectDB

  !------------------------------------------  ExpandSDNames  -----
  subroutine ExpandSDNames ( directData, sdName )
    use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test

    ! This routine adds to the sdNames arrays for a DirectWrite if necessary

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    character(len=*), intent(in) :: sdName

    ! Local variables
    logical :: alreadyThere
    integer :: i
    integer :: newsize
    character(len=80), dimension(:), pointer :: sdNames
    ! Check if the sdName already there
    ! print *, 'Check if the sdName already there'
    nullify ( sdNames )
    if ( associated(directData%sdNames) ) then
      newSize=SIZE(directData%sdNames)+1
      ! alreadyThere = any(sdName == directData%sdNames)
      alreadyThere = .false.
      do i=1, size(directData%sdNames)
        alreadyThere = alreadyThere .or. (sdName == directData%sdNames(i))
      end do
      ! if ( alreadyThere ) print *, trim(sdName), ' already there'
      if ( alreadyThere ) return
    else
      ! print *, 'Allocating for ', trim(sdName)
      newSize=1
    end if
    directData%NsdNames = newsize
    ! print *, ' allocating ', newsize
    call allocate_test ( sdNames, newSize, "sdNames", moduleName )
    if ( newSize>1 ) sdNames(1:newSize-1) = directData%sdNames
    if ( associated(directData%sdNames) ) &
      call deallocate_test ( directData%sdNames, "directData%sdNames", &
        & moduleName )
    directData%sdNames => sdNames
    directData%sdNames(newSize) = sdName
  end subroutine ExpandSDNames

  !------------------------------------------  FileNameToID  -----
  function FileNameToID ( fileName, DataBase )  result(ID)

    ! Given filename, returns index; if name not found in db, returns 0

    ! Dummy arguments
    type (DirectData_T), dimension(:), pointer :: DATABASE
    character(len=*), intent(in) :: FileName
    integer                      :: ID

    ! Local variables
    ! Executable
    id = 0
    if ( .not. associated(dataBase) .or. len_trim(filename) < 1 ) return
    do id =1, size(database)
      if ( fileName == dataBase(id)%fileName ) exit
      if ( fileName == dataBase(id)%fileNameBase ) exit
    enddo
    if ( id > size(database) ) id = 0
  end function FileNameToID

  !------------------------------------------  SetupNewDirect  -----
  subroutine SetupNewDirect ( directData, NsdNames )

    ! This routine sets up the arrays for a directWrite

    ! Dummy arguments
    type (DirectData_T), intent(inout)  :: directData
    integer, intent(in) :: NsdNames            ! Dimensions

    ! Allocate the sdNames
    directData%autoType = 0
    directData%Handle = 0
    directData%Type = 0
    directData%NSDNames = NSDNames
    directData%FileNameBase = ' '
    directData%FileName = ' '
    nullify(directData%sdNames)
    if ( NSDNames <= 0 ) return
    call allocate_test ( directData%sdNames, NsdNames, "directData%sdNames", &
         & ModuleName )
  end subroutine SetupNewDirect

! =====     Private Procedures     =====================================
  ! ---------------------------------------------  vectorValue_to_l2gp  -----
  subroutine vectorValue_to_l2gp ( quantity, &
    & precision, quality, status, convergence, AscDescMode, &
    & l2gp, name, chunkNo, offset, firstInstance, lastInstance )
    use Dump_0, only: Dump
    use Intrinsic, only: L_None, Lit_Indices, L_GPH
    use L2GPData, only: AscDescModeIsField, DescendingRange, L2GPData_T, RGP, &
      & SetupNewL2GPRecord
    use QuantityTemplates, only: Dump
    ! Args:
    type (VectorValue_T), intent(in) :: quantity
    type (VectorValue_T), pointer    :: precision
    type (VectorValue_T), pointer    :: quality
    type (VectorValue_T), pointer    :: status
    type (VectorValue_T), pointer    :: convergence
    type (VectorValue_T), pointer    :: AscDescMode
    type (L2GPData_T)                :: l2gp
    character(len=*), intent(in)     :: name
    integer, intent(in)              :: chunkNo
    type (Interval_T)                :: myRange
    integer, intent(in), optional    :: offset
    integer, intent(in), optional    :: firstInstance
    integer, intent(in), optional    :: lastInstance

    ! Local variables
    logical :: DEEBUG
    integer :: firstProfile
    real(rv) :: HUGERGP
    integer :: lastProfile
    integer :: noSurfsInL2GP
    integer :: noFreqsInL2GP
    integer :: nProfiles
    logical :: profiled  ! Print grep-able info about overlaps discarded
    character(len=32)     :: QtyName
    integer :: useFirstInstance
    integer :: useLastInstance
    logical :: verbose
    ! Executable
    profiled = BeVerbose ( 'profiled', -1 )
    verbose = BeVerbose ( 'direct', -1 )
    deebug = LetsDebug ( 'direct', 0 )
    
    hugeRgp = real ( huge(0.0_rgp), rv )

    ! Work out what to do with the first and last instance information
    
    if ( present(firstInstance) ) then
      useFirstInstance = firstInstance
    else
      useFirstInstance = quantity%template%noInstancesLowerOverlap+1
    end if

    if ( present(lastInstance) ) then
      useLastInstance = lastInstance
    else
      useLastInstance = quantity%template%noInstances- &
        & quantity%template%noInstancesUpperOverlap
    end if

    ! noOutputInstances = useLastInstance-useFirstInstance+1
    ! Now create an empty L2GP record with this dimension

    if (any(quantity%template%verticalCoordinate == &
      & (/l_Pressure, l_Zeta, l_geodAltitude /) )) then
      noSurfsInL2GP = quantity%template%noSurfs
    else
      noSurfsInL2GP = 0
    end if

    if ( quantity%template%frequencyCoordinate == l_None) then
       noFreqsInL2GP = 0
    else
       noFreqsInL2GP = quantity%template%noChans
    end if
    
    if ( present(offset) ) then
      firstProfile=offset+1
    else
      firstProfile=1
    endif
    nProfiles = useLastInstance - useFirstInstance + 1
    lastProfile = firstProfile - 1 + nProfiles

    call get_string( lit_indices(quantity%template%QuantityType ), &
      & QtyName, strip=.true. )      

    if ( DEEBUG ) print *, 'noFreqsInL2GP, noSurfsInL2GP, lastProfile: ', &
      & noFreqsInL2GP, noSurfsInL2GP, lastProfile
      
    if ( profiled ) then
      call output ( '## ## ' // trim(QtyName)                , advance='no' )
      call Blanks ( 2 )
      call output ( ChunkNo                                  , advance='no' )
      call Blanks ( 2 )
      call output ( quantity%template%noInstancesLowerOverlap, advance='no' )
      call Blanks ( 2 )
      call output ( quantity%template%noInstancesUpperOverlap, advance='no' )
      call Blanks ( 2 )
      call output ( firstProfile                             , advance='no' )
      call Blanks ( 2 )
      call output ( lastProfile                              , advance='no' )
      call NewLine
    endif
    call SetupNewl2gpRecord ( l2gp, noFreqsInL2GP, noSurfsInL2GP, lastProfile )
    ! Setup the standard stuff, only pressure as it turns out.
    select case ( quantity%template%verticalCoordinate )
    case ( l_Pressure ) 
      l2gp%pressures = quantity%template%surfs(:,1)
    case ( l_Zeta ) 
      l2gp%pressures = 10.0**(-quantity%template%surfs(:,1))
    case default
      call get_string( lit_indices(quantity%template%verticalCoordinate ), &
        & l2gp%verticalCoordinate, strip=.true. )      
      if ( noSurfsInL2GP > 0 ) l2gp%pressures = quantity%template%surfs(:,1)
      ! Don't warn about quantities with no vertically-resolved levels,
      ! e.g. Column amounts
      if ( noSurfsInL2GP > 1 ) call MLSL2Message( MLSMSG_Warning, ModuleName, &
        & 'Converting qty with non-pressure vertical coordinate ' // &
        & l2gp%verticalCoordinate )
      if ( verbose ) call dump( l2gp%pressures, 'vertical coordinates' )
      if ( deebug  ) call dump( quantity%template )
    end select
    
    ! We must choose a custom FillValue for GPH because
    ! -999.99 is a possibly legitimate value
    if ( quantity%template%QuantityType == l_gph .or. &
      & index( lowercase(name), 'gph' ) > 0 ) &
      & l2gp%MissingL2GP = L2Options%GPH_MissingValue
    ! It inherits its quantity type from the quantity template
    ! l2gp%quantityType=quantity%template%quantityType
    ! Do something about frequency
    if ( associated ( quantity%template%frequencies ) ) then
      l2gp%frequency = quantity%template%frequencies
    else
      l2gp%frequency = 0.0
    end if
    ! call ExpandL2GPDataInPlace ( l2gp, &
    !  & lastProfile-firstProfile+1 )

    ! Now copy the information from the quantity to the l2gpData
    L2GP%nTimesTotal = quantity%template%grandTotalInstances

    ! name is an integer, but L2GP%name is Character data
    l2gp%nameIndex = 0
    l2gp%name = name

    ! Now fill the data, first the geolocation
    l2gp%latitude(firstProfile:lastProfile) = &
      & quantity%template%geodLat(1,useFirstInstance:useLastInstance)
    l2gp%longitude(firstProfile:lastProfile) = &
      & quantity%template%lon(1,useFirstInstance:useLastInstance)
    l2gp%solarTime(firstProfile:lastProfile) = &
      & quantity%template%solarTime(1,useFirstInstance:useLastInstance)
    l2gp%solarZenith(firstProfile:lastProfile) = &
      & quantity%template%solarZenith(1,useFirstInstance:useLastInstance)
    l2gp%losAngle(firstProfile:lastProfile) = &
      & quantity%template%losAngle(1,useFirstInstance:useLastInstance)
    ! If the quantity isn't stacked, Phi is meaningless, and has zero size.
    if ( allocated(quantity%template%phi) .and. quantity%template%stacked ) then
      l2gp%geodAngle(firstProfile:lastProfile) = &
        & quantity%template%phi(1,useFirstInstance:useLastInstance)
    else
      l2gp%geodAngle(firstProfile:lastProfile) = 0
    end if
    l2gp%time(firstProfile:lastProfile) = &
      & quantity%template%time(1,useFirstInstance:useLastInstance)
    l2gp%chunkNumber(firstProfile:lastProfile) = chunkNo
    if ( deebug ) then
      call outputNamedValue( 'firstProfile, lastProfile', (/firstProfile, lastProfile /) )
      call outputNamedValue( 'chunkNos', l2gp%chunkNumber(firstProfile:lastProfile) )
    endif

    ! Now the various data quantities.

    l2gp%l2gpValue(:,:,firstProfile:lastProfile) = &
      & reshape ( max ( -hugeRgp, min ( hugeRgp, &
      &   quantity%values(:,useFirstInstance:useLastInstance) ) ), &
      &  (/max(l2gp%nFreqs,1),max(l2gp%nLevels,1),lastProfile-firstProfile+1/))
    if (associated(precision)) then
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = &
        & reshape ( max ( -hugeRgp, min ( hugeRgp, &
        &   precision%values(:,useFirstInstance:useLastInstance) ) ), &
        &  (/max(l2gp%nFreqs,1),max(l2gp%nLevels,1),lastProfile-firstProfile+1/))
    else
      l2gp%l2gpPrecision(:,:,firstProfile:lastProfile) = 0.0
    end if
    if (associated(quality)) then
      l2gp%quality(firstProfile:lastProfile) = &
        & quality%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%quality(firstProfile:lastProfile) = 0.0
    endif
    if (associated(status)) then
      l2gp%status(firstProfile:lastProfile) = &
        & status%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%status(firstProfile:lastProfile) = 0
    endif
    if (associated(convergence)) then
      l2gp%convergence(firstProfile:lastProfile) = &
        & convergence%values(1,useFirstInstance:useLastInstance)
    else
      l2gp%convergence(firstProfile:lastProfile) = 0.0
    endif
    if ( .not. AscDescModeIsField ) then
      ! The field has been dropped from l2gp swaths
    elseif (associated(AscDescMode)) then
      ! This sets the AscDescMode field to 
      ! +1 for values of +1
      ! -1 for values of -1
      l2gp%AscDescMode(firstProfile:lastProfile) = &
        & merge( 1, -1, &
        & (AscDescMode%values(1,useFirstInstance:useLastInstance) > 0._rv) &
        & )
    elseif( DescendingRange%Bottom < DescendingRange%Top ) then
      ! Turnaround points are assigned Mode values according to
      ! 90 deg  -> descending
      ! 270 deg -> ascending
      myRange%Bottom = DescendingRange%Bottom
      myRange%Top    = DescendingRange%Top / (1. + 1.e-6)
      l2gp%AscDescMode(firstProfile:lastProfile) = &
        & merge( -1, 1, &
        & inRange( &
        &  mod( l2gp%geodAngle(firstProfile:lastProfile), 360._rgp ) &
        & , myRange ) &
        & )
      if ( DEEBUG ) then
        call Dump( (/myRange%Bottom, myRange%top/), 'myRange' )
        call Dump( mod( l2gp%geodAngle(firstProfile:lastProfile), 360._rgp ), &
          & 'Orbit Angle' )
        call Dump( l2gp%AscDescMode(firstProfile:lastProfile), 'Mode' )
      endif
    else
      l2gp%AscDescMode(firstProfile:lastProfile) = 0
    endif
    if ( DEEBUG ) then
      print *, 'Vector converted to l2gp; name: ', trim(name)
      print *, 'firstProfile, lastProfile: ', firstProfile, lastProfile
      print *, 'useFirstInstance, useLastInstance: ', useFirstInstance, useLastInstance
      call outputNamedValue ( 'lons', l2gp%longitude(firstProfile:lastProfile) )
      call outputNamedValue ( 'lats', l2gp%latitude(firstProfile:lastProfile) )
    endif
    if ( all(l2gp%chunkNumber == -999) ) &
      & call output( 'all converted chunk numbers are -999', advance='yes' )
    if ( associated(quantity%BinNumber) ) then
      allocate( l2gp%BinNumber(1:lastProfile) )
      l2gp%BinNumber(1:lastProfile) = &
        & quantity%BinNumber(useFirstInstance:useLastInstance)
      if ( DEEBUG ) call Dump( l2gp%BinNumber, 'l2gp bin numbers' )
    else
      call output( 'Bin numbers not allocated', advance='yes' )
    endif
    if ( associated(quantity%MAF) ) then
      allocate( l2gp%MAF(1:lastProfile) )
      l2gp%MAF(1:lastProfile) = quantity%MAF(useFirstInstance:useLastInstance)
      if ( DEEBUG ) call Dump( l2gp%MAF, 'l2gp MAFs' )
    else
      call output( 'MAFs not allocated', advance='yes' )
    endif
  end subroutine vectorValue_to_l2gp

  ! ---------------------------------------------  ANNOUNCE_ERROR  -----
  subroutine ANNOUNCE_ERROR ( Where, Full_message, Code, Penalty )

    use Lexer_Core, only: Print_Source
    use Tree, only: Where_At=>where

    ! Args:
    integer, intent(in) :: Where   ! Tree node where error was noticed
    character(LEN=*), intent(in) :: Full_Message
    integer, intent(in), optional :: Code    ! Code for error message
    integer, intent(in), optional :: Penalty
    integer :: myPenalty

    myPenalty = 1
    if ( present(penalty) ) myPenalty = penalty
    error = max(error,myPenalty)

    call output ( '***** At ' )
    if ( where > 0 ) then
      call print_source ( where_at(where) )
    else
      call output ( '(no lcf node available)' )
    end if
    call output ( ModuleName // ' complained: ' )

    call output ( trim(full_message), advance='yes', &
      & from_where=ModuleName )
    if ( present(code) ) then
      select case ( code )
      end select
    end  if
  end subroutine ANNOUNCE_ERROR

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module DirectWrite_m

! $Log$
! Revision 2.98  2021/07/22 23:12:54  pwagner
! Swicth setting -Sprofiled prints grepable info about chunk overlaps
!
! Revision 2.97  2021/06/10 23:47:44  pwagner
! When creating l2gps from qties, copy their BinNumber and MAF, too
!
! Revision 2.96  2020/03/04 21:27:32  pwagner
! Skip Warnings about Column amounts: already know they lack vertical coords
!
! Revision 2.95  2020/02/13 21:27:45  pwagner
! Fix errors relating to separate MissingValue for GPH
!
! Revision 2.94  2019/02/13 18:58:42  pwagner
! New GPH_MissingValue field of L2Options can now be st to a value other than default -999.99
!
! Revision 2.93  2018/11/12 23:13:26  pwagner
! Deprecated AscDescMode
!
! Revision 2.92  2018/07/27 23:17:16  pwagner
! Renamed level 2-savvy MLSMessage MLSL2Message
!
! Revision 2.91  2018/04/19 01:14:16  vsnyder
! Remove USE statements for unused names
!
! Revision 2.90  2018/04/13 00:20:42  pwagner
! Plain hdf DirectWrites and -Reads are now 'auto'
!
! Revision 2.89  2018/01/12 00:22:34  pwagner
! Reduce amount written by lowest verbose level
!
! Revision 2.88  2017/08/10 22:47:36  pwagner
! Use WriteHDF5Data from L2AuxData
!
! Revision 2.87  2017/07/27 17:02:08  pwagner
! If an sdName contains one or more '/', will attempt ot open or create nested hdf5 groups
!
! Revision 2.86  2017/02/24 19:47:52  pwagner
! Sleep time now same as parallel%delay
!
! Revision 2.85  2016/09/21 00:39:46  pwagner
! Usually dump DirectDB as a table
!
! Revision 2.84  2016/09/07 22:46:49  pwagner
! Removed unused QuantityType component from L2GPData type
!
! Revision 2.83  2016/08/09 21:39:40  pwagner
! Removed unused variables
!
! Revision 2.82  2016/05/25 00:07:18  pwagner
! To appease NAG which disliked ambiguity in generic DirectWrites
!
! Revision 2.81  2016/05/18 19:05:44  vsnyder
! Change HGrids database from an array of HGrid_T to an array of pointers
! to HGrid_T using the new type HGrids_T.
!
! Revision 2.80  2016/02/29 19:49:29  pwagner
! Usleep got from machine module instead of being an external
!
! Revision 2.79  2015/10/14 23:24:46  pwagner
! Cope better with l2gp qty with non-pressure vertical coordinate
!
! Revision 2.78  2015/10/06 17:37:05  pwagner
! Added more error checking; allow for altitude vertical coordinate
!
! Revision 2.77  2015/09/17 23:24:50  pwagner
! Passes Max chunk size for l2gp DirectWrites
!
! Revision 2.76  2015/07/29 00:29:54  vsnyder
! Convert Phi from pointer to allocatable
!
! Revision 2.75  2015/07/14 23:32:01  pwagner
! label and inputFile fields in DirectWrite
!
! Revision 2.74  2015/04/29 01:15:09  vsnyder
! Calculate TrueRank correctly, allow to use it if rank<0
!
! Revision 2.73  2015/04/25 02:10:36  vsnyder
! Spiff a message
!
! Revision 2.72  2015/04/21 17:53:03  pwagner
! May DirectRead quantity geolocations
!
! Revision 2.71  2015/04/09 22:19:59  pwagner
! We may DirectRead a quantity type
!
! Revision 2.70  2015/04/07 02:56:29  vsnyder
! Add warning about rank being wrong in DirectWrite_Quantity.  Add ability to
! detect rank if requested rank is <= 0.  Add FrequencyCoordinate.
!
! Revision 2.69  2015/03/31 21:01:57  pwagner
! rank is a new field for DirectWrite-ing quantity values as, say, rank3; now write qty attributes, too
!
! Revision 2.68  2015/03/28 02:30:58  vsnyder
! Added stuff to trace allocate/deallocate addresses.  Paul added
! DirectWriteQuantity stuff.
!
! Revision 2.67  2015/02/05 21:42:23  vsnyder
! Don't use Phi for unstacked quantities
!
! Revision 2.66  2014/12/11 21:27:04  pwagner
! Turn off more printing unless verbose or debug
!
! Revision 2.65  2014/12/10 23:02:59  pwagner
! Correct erroneous calculation of AscDescMode; now based on geodAngle
!
! Revision 2.64  2014/11/04 01:25:41  pwagner
! May switch on verbose, deebug modes with -Sdirect[n]
!
! Revision 2.63  2014/10/02 17:22:23  pwagner
! DirectWrite now responsible for writing MiscNotes for l2gp files
!
! Revision 2.62  2014/09/05 00:40:24  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.61  2014/04/10 00:47:18  pwagner
! More consistent in when to announce creating files, datasets
!
! Revision 2.60  2014/04/07 18:03:03  pwagner
! May specify AscDescMode when DirectWrite-ing swaths
!
! Revision 2.59  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.58  2013/11/20 00:56:19  pwagner
! Reduce default dump of directWriteDB to just file names
!
! Revision 2.57  2013/09/24 23:47:22  vsnyder
! Use Where instead of Source_Ref for messages
!
! Revision 2.56  2013/08/12 23:49:41  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.55  2013/05/08 20:17:35  pwagner
! Changed sleep to usleep to link with NAG
!
! Revision 2.54  2013/05/07 22:23:30  pwagner
! A workaround to prevent unexpected Missing/Fill values
!
! Revision 2.53  2012/05/10 00:48:19  pwagner
! A little less clumsy at creating files, swaths
!
! Revision 2.52  2012/03/12 17:09:44  pwagner
! Use new writeAPrioriAttributes api
!
! Revision 2.51  2012/02/24 21:19:15  pwagner
! May DirectWrite a /single instance only
!
! Revision 2.50  2011/12/15 01:52:55  pwagner
! 'num' option names quantity in entireVector DirectWrites numerically
!
! Revision 2.49  2011/11/01 21:44:31  pwagner
! May optionally add quantity Attributes via options='A' field
!
! Revision 2.48  2011/05/09 18:06:59  pwagner
! May directwrite an entire vector
!
! Revision 2.47  2010/02/04 23:12:44  vsnyder
! Remove USE or declaration for unreferenced names
!
! Revision 2.46  2010/02/04 19:10:26  pwagner
! Change condition on directwrite db from associated to zero size
!
! Revision 2.45  2010/01/08 00:11:19  pwagner
! Added ability to write MLSFile_T fields as attributes
!
! Revision 2.44  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.43  2009/04/23 23:02:52  pwagner
! May specify upperOverlap or lowerOverlap in DirectWrites
!
! Revision 2.42  2008/09/20 00:03:00  pwagner
! Added print statement to not_used_here
!
! Revision 2.41  2006/10/11 22:58:00  pwagner
! Will write convergence ratio as another quality-like field of l2gp
!
! Revision 2.40  2006/06/20 00:13:34  pwagner
! Fwm files should haves MAFs matching radiance files
!
! Revision 2.39  2006/05/11 19:39:19  pwagner
! Should not segment fault if dumping empty DB
!
! Revision 2.38  2006/04/20 23:24:12  pwagner
! More bugs squashed related to extra-range MAFs; one crashed final chunk
!
! Revision 2.37  2006/04/19 20:48:13  pwagner
! Undid most of the changes regarding extra MAFs; perhaps fixed bugs
!
! Revision 2.36  2006/04/12 22:20:42  pwagner
! Attempted workaround for hdf5-1.6.5 bugs rewriting string attributes
!
! Revision 2.35  2006/04/11 23:33:14  pwagner
! Fixed bug which added excess profiles
!
! Revision 2.34  2006/01/26 00:34:50  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.33  2005/08/19 23:28:02  pwagner
! Trying to avoid possibility of Lahey-causes memory leak
!
! Revision 2.32  2005/07/12 17:38:57  pwagner
! Writes APriori File names as an attribute to every l2gp file
!
! Revision 2.31  2005/06/22 18:57:01  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.30  2005/06/14 20:41:17  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.29  2004/11/29 21:52:41  livesey
! Bug fix for handling cases where no forward models are defined.
!
! Revision 2.28  2004/08/03 18:02:01  pwagner
! Sets fillValue for l2aux type
!
! Revision 2.27  2004/07/22 20:42:57  cvuu
! May write ForwardModel Names as file-level attributes and fix the phase names
!
! Revision 2.26  2004/06/29 18:05:26  pwagner
! May write phase, section names as file-level attributes
!
! Revision 2.25  2004/06/10 00:58:45  vsnyder
! Move FindFirst, FindNext from MLSCommon to MLSSets
!
! Revision 2.24  2004/05/19 19:16:09  vsnyder
! Move MLSChunk_t to Chunks_m
!
! Revision 2.23  2004/05/05 21:31:48  pwagner
! More debug printing
!
! Revision 2.22  2004/03/03 19:25:45  pwagner
! Fixed poorly understood prob with single chunks; more initializing in setup
!
! Revision 2.21  2004/02/19 23:54:26  pwagner
! Dumps during directWrites if l2gp or l2aux switch set
!
! Revision 2.20  2004/02/11 17:23:25  pwagner
! l2gp status an integer, not a char
!
! Revision 2.19  2004/02/10 19:30:55  pwagner
! Cures serial directWrites from writing more than one chunk at a time
!
! Revision 2.18  2004/02/05 23:38:41  pwagner
! Writes attributes to directwrite l2aux file
!
! Revision 2.17  2004/01/23 01:09:48  pwagner
! Only directwrite files entered in global settings eligible to be auto-sourced
!
! Revision 2.16  2004/01/22 06:38:26  livesey
! Typo fix
!
! Revision 2.15  2004/01/22 00:56:35  pwagner
! Fixed many bugs in auto-distribution of DirectWrites
!
! Revision 2.14  2003/12/03 17:50:54  pwagner
! L2GP tracks both nTimes (for this slave) and nTimesTotal (done by all)
!
! Revision 2.13  2003/11/14 23:38:45  pwagner
! Uses DirectWrite databse in preference to repeated calls to toolkit
!
! Revision 2.12  2003/08/28 23:51:39  livesey
! Various bug fixes and simplifications
!
! Revision 2.11  2003/08/01 20:39:34  pwagner
! Skip warning mesg about grandTotalInstances when 0
!
! Revision 2.10  2003/07/18 16:05:26  pwagner
! Stops printing warnings after 40 times
!
! Revision 2.9  2003/07/15 23:41:47  pwagner
! l2aux always rank 3 unless MAYCOLLAPSEDIMS; disabled most printing
!
! Revision 2.8  2003/07/09 21:49:53  pwagner
! Tries to figure out in advance whether to create swath
!
! Revision 2.7  2003/07/07 21:03:43  pwagner
! Removed DirectWrite_L2Aux_FileName; tries to deal sensibly with all-overlap chunks
!
! Revision 2.6  2003/07/07 17:32:30  livesey
! New approach to DirectWrite
!
! Revision 2.5  2003/07/02 00:55:27  pwagner
! Some improvements in DirectWrites of l2aux, l2gp
!
! Revision 2.4  2003/06/25 18:24:57  pwagner
! A fix to vectorValue_to_l2gp
!
! Revision 2.3  2003/06/24 23:53:27  pwagner
! Allows unassociated precisions as args to direct write
!
! Revision 2.2  2003/06/23 23:55:17  pwagner
! Added DirectData_T to keep track of data written directly
!
! Revision 2.1  2003/06/20 19:43:16  pwagner
! First commit
!
