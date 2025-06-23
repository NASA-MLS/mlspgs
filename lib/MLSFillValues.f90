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
module MLSFillValues              ! Some FillValue-related stuff
!=============================================================================

  use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
  use Ieee_Arithmetic, only: IsFinite => Ieee_Is_Finite, IsNaN => Ieee_Is_NaN, &
    & Ieee_Is_Finite, Ieee_Is_NaN
  use MLSCommon, only: Fill_Signal, Inf_Signal, NaN_Signal, UndefinedValue, &
    & Is_What_Ieee, MLSFill_T, MLSFills
  use MLSKinds ! Everything
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSFinds, only: FindFirst, FindLast
  use MLSStrings, only: Lowercase
  use Output_M, only: Blanks, Output

  implicit none

  private

  public :: AddMLSFillToDatabase, Dump
  public :: FillFunction, InfFunction, NaNFunction
  public :: FilterValues
  public :: IsFillValue
  public :: IsFinite, IsInfinite, IsNaN
  public :: Monotonize
  public :: RemoveFillValues, ReorderFillValues, ReplaceFillValues
  public :: RoundUpOrDown
  public :: WhereAreTheFills, WhereAreTheNaNs, WhereAreTheInfs

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! This module contains stuff related to:
! FillValues, NaNs, Infinities
! Rounding up or dowwn, and 
! Bridging or Monotonize-ing values that interrupt a monotonic array

! Fill values are values in an array which have been so marked as to be ignored
! Sometimes they are preassigned values, as with mls hdf[eos] datasets
! where we use -999.99
! Sometimes we give them values explicitly to mark where an agorithm failed
! or where a an excursion outside its region of validity was detected
! Sometimes we simply want to consider only non-zero values

! Could we generalize this to mark as Fill values any values above or
! below a threshold, or whose absolute value is above a threshold?

! === (start of toc) ===                                                 
!     c o n t e n t s
!     - - - - - - - -

!         Functions, operations, routines
! addMLSFillToDatabase
!                   Create or add to a database of MLSFills, to be used in 
!                   deciding whether an argument is a Fill, broadly defined
! Dump              Dump the MLSFills database
! FillFunction      Returns the Fill value
! InfFunction       Returns the Infinite value
! NaNFunction       Returns the NaN value
! FilterValues      Filters entries in two arrays
! IsFillValue       Returns true if argument is FillValue
! IsFinite          Returns true if argument is finite
! IsInfinite        Returns true if argument is infinite
! IsNaN             Returns true if argument is NaN
! Monotonize        Replace any non-monotonic elements
! RemoveFillValues  Removes FillValues from an array
!                     returning a new array smaller in size
! ReorderFillValues Reorders FillValue entries at the end of an array
! ReplaceFillValues Replaces FillValue entries in an array
!                     either by interpolating or stutter the non-Fill values
! RoundUpOrDown     Rounds an arg up or down depending on fraction
! WhereAreTheFills  Find which array elements are Fill values
! WhereAreTheInfs   Find which array elements are Inf
! WhereAreTheNaNs   Find which array elements are NaN
! === (end of toc) ===                                                   

! === (start of api) ===
! int addMLSFillToDatabase( MLSFill_T arg, [MLSFill_T yourDatabase(:)] )
! int addMLSFillToDatabase( r8 FillValue, &
!   [r8 tol], [char* condition], [MLSFill_T yourDatabase(:)] )
! Dump ( [MLSFill_T Database(:)] )
! nprec FillFunction( nprec arg )
! nprec InfFunction( nprec arg )
! nprec NaNFunction( nprec arg )
! filterValues ( a, atab, b, btab, log warn, [fillValue], [precision] )
! log IsFillValue ( a, [fillvalue], [log strict] )
! log IsFinite ( a )
! log Isinfinite ( a )
! log IsNaN ( a )
! Monotonize( values, [Period], [FillValue], [log strict] )
! RemoveFillValues ( array, FillValue, newArray, [second], [newSecond] )
! ReorderFillValues ( values, FillValue )
! ReplaceFillValues ( values, FillValue, [newValues], [newFill], [char* options] )
! nprec RoundUpOrDown( nprec value )
! WhereAreTheFills ( array, [which(:)], [int howMany], [mode], [inds(:,:)] )
! WhereAreTheInfs ( array, [which(:)], [int howMany], [mode], [inds(:,:)] )
! WhereAreTheNaNs ( array, [which(:)], [int howMany], [mode], [inds(:,:)] )
!
! These last 3 search for the targeted numerical type (Fill, Inf, or NaN)
! according to mode:
!     'any'  Finds which column or plane (corresponding to last index) has any
!     'all'  Finds which column or plane (corresponding to last index) is all
!     'ind'  Finds index pair (for rank 2) or triple (for rank 3) of each
! More Notes:
! (1) for rank r array, inds should be dimensioned (r,n) where n is the max
! number of Fills or whatever you expect
! (2) For mode = 'all' or 'any', inds is ignored, howMany is number of columns
!     or planes
! (3) For mode = 'ind', which is ignored, howMany is number of individuals
! Example: we're searching for NaNs; the array looks like the following
!
!  0   0   0 NaN
! NaN  0   0 NaN
! NaN NaN  0 NaN
!  0   0   0 NaN
!
!      Then
! mode      howMany     which        inds
! 'any'       3        1, 2, 4     (ignored)
! 'all'       1           4        (ignored)
! 'ind'       7       (ignored)  (1,4),(2,1),(2,4),(3,1),(3,2),(3,4),(4,4)
! We show the pairs of indices under inds above, but as returned
! by the routine the 2x7 array of inds will be laid out as follows
!  1  2  2  3  3  3  4
!  4  1  4  1  2  4  4
!
! === (end of api) ===                                                 
  interface addMLSFillToDatabase
    module procedure addFillValueToDatabase, addMLSFillTypeToDatabase
  end interface

  interface BridgeMissingValues
    module procedure BridgeMissingValues_1dr4, BridgeMissingValues_1dr8, BridgeMissingValues_1dint
    module procedure BridgeMissingValues_2dr4, BridgeMissingValues_2dr8, BridgeMissingValues_2dint
    module procedure BridgeMissingValues_3dr4, BridgeMissingValues_3dr8, BridgeMissingValues_3dint
  end interface

  interface Dump
    module procedure DumpMLSFillsDatabase
  end interface

  interface FillFunction
    module procedure FillFunction_REAL, FillFunction_DOUBLE
  end interface

  interface InfFunction
    module procedure InfFunction_REAL, InfFunction_DOUBLE
  end interface

  interface NaNFunction
    module procedure NaNFunction_REAL, NaNFunction_DOUBLE
  end interface

  interface FilterValues
    module procedure FilterValues_REAL, FilterValues_DOUBLE
    module procedure FilterValues_REAL_2d, FilterValues_DOUBLE_2d
    module procedure FilterValues_REAL_3d, FilterValues_DOUBLE_3d
  end interface
  
  interface IsFillValue
    module procedure IsFillValue_r4, IsFillValue_r8, IsFillValue_int
  end interface

  interface IsFinite
    module procedure IsFinite_INTEGER
  end interface
  
  interface IsInfinite
    module procedure IsInfinite_REAL, IsInfinite_DOUBLE, IsInfinite_INTEGER
  end interface
  
  interface IsNaN
    module procedure IsNaN_INTEGER
    module procedure IsNaN_CHARACTER
  end interface
  
  interface Monotonize
    module procedure Monotonize_1dr4, Monotonize_1dr8, Monotonize_1dint
    module procedure Monotonize_2dr4, Monotonize_2dr8, Monotonize_2dint
    module procedure Monotonize_3dr4, Monotonize_3dr8, Monotonize_3dint
  end interface

  interface RemoveFillValues
    module procedure RemoveFill1d_r4, RemoveFill1d_r8, RemoveFill1d_int
    ! I'm not certain it makes sense to do this except in 1-d
    module procedure RemoveFill2d_r4, RemoveFill2d_r8, RemoveFill2d_int
    module procedure RemoveFill3d_r4, RemoveFill3d_r8, RemoveFill3d_int
  end interface

  interface ReorderFillValues
    module procedure ReorderFillValues_r4, ReorderFillValues_r8, ReorderFillValues_int
  end interface

  interface ReplaceFillValues
    module procedure ReplaceFill1d_r4, ReplaceFill1d_r8, ReplaceFill1d_int
    module procedure ReplaceFill2d_r4, ReplaceFill2d_r8, ReplaceFill2d_int
    module procedure ReplaceFill3d_r4, ReplaceFill3d_r8, ReplaceFill3d_int
  end interface

  interface RoundUpOrDown
    module procedure roundUpOrDown_r4, roundUpOrDown_r8, roundUpOrDown_int
  end interface

  interface swap
    module procedure swap_int, swap_r4, swap_r8
  end interface

  interface WhereAreTheFills
    module procedure WhereAreTheFills_REAL, WhereAreTheFills_DOUBLE
    module procedure WhereAreTheFills_REAL_2d, WhereAreTheFills_DOUBLE_2d
    module procedure WhereAreTheFills_REAL_3d, WhereAreTheFills_DOUBLE_3d
    module procedure WhereAreTheFills_INTEGER
    module procedure WhereAreTheFills_INTEGER_2d, WhereAreTheFills_INTEGER_3d
  end interface
  
  interface WhereAreTheInfs
    module procedure WhereAreTheInfs_REAL, WhereAreTheInfs_DOUBLE
    module procedure WhereAreTheInfs_REAL_2d, WhereAreTheInfs_DOUBLE_2d
    module procedure WhereAreTheInfs_REAL_3d, WhereAreTheInfs_DOUBLE_3d
    module procedure WhereAreTheInfs_REAL_4d, WhereAreTheInfs_DOUBLE_4d
    module procedure WhereAreTheInfs_INTEGER
    module procedure WhereAreTheInfs_INTEGER_2d, WhereAreTheInfs_INTEGER_3d
  end interface
  
  interface WhereAreTheNaNs
    module procedure WhereAreTheNaNs_REAL, WhereAreTheNaNs_DOUBLE
    module procedure WhereAreTheNaNs_REAL_2d, WhereAreTheNaNs_DOUBLE_2d
    module procedure WhereAreTheNaNs_REAL_3d, WhereAreTheNaNs_DOUBLE_3d
    module procedure WhereAreTheNaNs_REAL_4d, WhereAreTheNaNs_DOUBLE_4d
    module procedure WhereAreTheNaNs_INTEGER
    module procedure WhereAreTheNaNs_INTEGER_2d, WhereAreTheNaNs_INTEGER_3d
  end interface
  
  interface WhereAreThey
    module procedure WhereAreThey_REAL, WhereAreThey_DOUBLE
    module procedure WhereAreThey_REAL_2d, WhereAreThey_DOUBLE_2d
    module procedure WhereAreThey_REAL_3d, WhereAreThey_DOUBLE_3d
    module procedure WhereAreThey_REAL_4d, WhereAreThey_DOUBLE_4d
    module procedure WhereAreThey_INTEGER
    module procedure WhereAreThey_INTEGER_2d, WhereAreThey_INTEGER_3d
  end interface
  
  ! This tolerance won't work w/o a little fudging 
  ! when fill values get really huge; e.g. gmao use 1.e15, gloria 1.e12
  ! The fudging will take the form of
  ! max( FILLVALUETOLERANCE, abs(FILLVALUE/100000) )
  ! where obviously 100000 is an arbitrary number
  ! Should we study this more carefully?
  real, parameter, private :: FILLVALUETOLERANCE = 0.2 ! Poss. could make it 1
  character(len=3), save   :: NaNString          = 'NaN'
  character(len=3), save   :: InfString          = 'Inf'
  logical, parameter       :: DEEBUG             = .false.
  logical, save            :: DONTINTERPOLATE    = .false. ! If true, stutter

contains

  !-------------------------------------------  AddFillValueToDatabase  -----
  integer function AddFillValueToDatabase( FillValue, &
    & tol, condition, yourDatabase )

    ! This function adds an FillValue data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where FillValue is put.

    ! Dummy arguments
    type (MLSFill_T), dimension(:), optional, pointer :: yourDatabase
    real(r8), optional, intent(in)                    :: tol
    real(r8), intent(in)                              :: FillValue
    character(len=*), optional, intent(in)            :: condition

    ! Local variables
    type (MLSFill_T)                        :: item
    ! Executable
    item%value = FillValue
    if ( present(condition) ) then
      item%condition = condition
    else
      item%condition = '='
    endif
    if ( present(tol) ) then
      item%tol = tol
    else
      item%tol = max( real(FILLVALUETOLERANCE, r8), abs(FillValue*1.d-6) )
    endif
    AddFillValueToDatabase = addMLSFillTypeToDatabase( item, yourDatabase )
  end function AddFillValueToDatabase

  !-------------------------------------------  addMLSFillTypeToDatabase  -----
  integer function addMLSFillTypeToDatabase( ITEM, yourDatabase )

    ! This function adds an MLSFill data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where MLSFill is put.

    ! Dummy arguments
    type (MLSFill_T), dimension(:), optional, pointer :: yourDatabase
    type (MLSFill_T), intent(in) :: ITEM

    ! Executable
    if ( present(yourDatabase) ) then
      addMLSFillTypeToDatabase = AddMLSFillToDefinite( item, yourDatabase )
      ! call dump( yourDatabase )
    else
      addMLSFillTypeToDatabase = AddMLSFillToDefinite( item, MLSFills )
      ! call dump( MLSFills )
    endif
  end function addMLSFillTypeToDatabase

  !-------------------------------------------  AddMLSFillToDefinite  -----
  integer function AddMLSFillToDefinite( ITEM, DATABASE )

    ! This function adds an MLSFill data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where MLSFill is put.

    ! Dummy arguments
    type (MLSFill_T), dimension(:), pointer :: DATABASE
    type (MLSFill_T), intent(in) :: ITEM

    ! Local variables
    type (MLSFill_T), dimension(:), pointer :: tempDatabase
    ! Executable
    include "addItemToDatabase.f9h" 

    AddMLSFillToDefinite = newSize
  end function AddMLSFillToDefinite

  ! ---------------------------------------------  Dump  -----
  ! This family of routines performs any Dumps appropriate
  ! to the args--in our case the MLSFill types
  subroutine DumpMLSFillsDatabase ( database )
    ! Args
    type(MLSFill_T), dimension(:), pointer, optional :: database
    ! Internal variables
    type(MLSFill_T), dimension(:), pointer :: dumpingdatabase
    integer :: i
    ! Executable
    if ( present(database) ) then
      dumpingdatabase => database
    else
      dumpingdatabase => MLSFills
    endif
    if ( .not. associated(dumpingdatabase) ) then
      call output( 'MLSFill database is not associated', advance='yes' )
      return
    endif
    call output( 'MLSFill database', advance='yes' )
    call output( '         i' )
    call blanks(12)
    call output( 'value')
    call blanks(4)
    call output( 'tolerance')
    call blanks(4)
    call output( '          condition', advance='yes' )
    do i=1, size(dumpingdatabase)
      call output( i )
      call blanks ( 4 )
      call output( dumpingdatabase(i)%value )
      call blanks ( 4 )
      call output( dumpingdatabase(i)%tol )
      call blanks ( 4 )
      call output( dumpingdatabase(i)%condition, advance='yes' )
    enddo
  end subroutine DumpMLSFillsDatabase

! -------------------------------------------------  FillFunction  -----
! Returns the fillValue
  function FillFunction_DOUBLE(arg) result(value)
    double precision, intent(in) :: arg
    real                         :: value
    value = UndefinedValue
  end function FillFunction_DOUBLE

  function FillFunction_REAL(arg) result(value)
    real, intent(in)        :: arg
    real                    :: value
    value = UndefinedValue
  end function FillFunction_REAL

! -------------------------------------------------  InfFunction  -----
! Returns the Inf Value
  function InfFunction_DOUBLE(arg) result(value)
    double precision, intent(in) :: arg
    real                         :: value
    read ( InfString, * ) value
  end function InfFunction_DOUBLE

  function InfFunction_REAL(arg) result(value)
    real, intent(in)        :: arg
    real                    :: value
    read ( InfString, * ) value
  end function InfFunction_REAL

! -------------------------------------------------  NaNFunction  -----
! Returns the NaN Value
  function NaNFunction_DOUBLE(arg) result(value)
    double precision, intent(in) :: arg
    real                         :: value
    read ( NaNString, * ) value
  end function NaNFunction_DOUBLE

  function NaNFunction_REAL(arg) result(value)
    real, intent(in)        :: arg
    real                    :: value
    read ( NaNString, * ) value
  end function NaNFunction_REAL

! -------------------------------------------------  FilterValues  -----
  subroutine filterValues_REAL ( a, ATAB, b, BTAB, warn, fillValue, precision )
    ! Return arrays filtered
    ! where corresponding precision array < 0
    ! or whose values are not finite
    ! "Filter" means offending elements set to 0. or optional fillValue
    ! Returned arrays are assigned values as appropriate
    
    ! An older version was able to reset any fillValues found to 0
    ! but that functionality has been removed
    ! It can be accomplished with ease in your code via
    ! where ( isFillValue(array) )
    !    array = 0.
    ! endwhere
    
    ! See also ReplaceFillValues
    integer, parameter :: RK = kind(0.0) ! Kind type parameter for default real
    character(*), parameter :: P = 'default reals'
    include 'FilterValues_1d.f9h'
  end subroutine filterValues_REAL

  subroutine filterValues_DOUBLE ( a, ATAB, b, BTAB, warn, fillValue, precision )
    ! Return arrays filtered etc.
    character(*), parameter :: P = 'double precision reals'
    integer, parameter :: RK = kind(0.0d0) ! Kind type parameter for double precision
    include 'FilterValues_1d.f9h'
  end subroutine filterValues_DOUBLE

  subroutine filterValues_REAL_2d(a, ATAB, b, BTAB, warn, fillValue, precision)
    ! Return arrays filtered etc.
    integer, parameter :: RK = kind(0.0) ! Kind type parameter for default real
    ! Args
    real(rk), dimension(:,:), intent(in)           :: a
    real(rk), dimension(:,:), intent(in)           :: b
    real(rk), dimension(:,:), intent(out)          :: atab
    real(rk), dimension(:,:), intent(out)          :: btab
    logical, intent(out)                           :: warn
    real(rk), optional, intent(in)                 :: fillValue
    real(rk), dimension(:,:), optional, intent(in) :: precision
    include 'FilterValues_nd.f9h'
  end subroutine filterValues_REAL_2d

  subroutine filterValues_DOUBLE_2d(a, ATAB, b, BTAB, warn, fillValue, precision)
    ! Return arrays filtered etc.
    integer, parameter :: RK = kind(0.0d0) ! Kind type parameter for double precision
    ! Args
    real(rk), dimension(:,:), intent(in)           :: a
    real(rk), dimension(:,:), intent(in)           :: b
    real(rk), dimension(:,:), intent(out)          :: atab
    real(rk), dimension(:,:), intent(out)          :: btab
    logical, intent(out)                           :: warn
    real(rk), optional, intent(in)                 :: fillValue
    real(rk), dimension(:,:), optional, intent(in) :: precision
    include 'FilterValues_nd.f9h'
  end subroutine filterValues_DOUBLE_2d

  subroutine filterValues_REAL_3d(a, ATAB, b, BTAB, warn, fillValue, precision)
    ! Return arrays filtered etc.
    integer, parameter :: RK = kind(0.0) ! Kind type parameter for default real
    ! Args
    real(rk), dimension(:,:,:), intent(in)           :: a
    real(rk), dimension(:,:,:), intent(in)           :: b
    real(rk), dimension(:,:,:), intent(out)          :: atab
    real(rk), dimension(:,:,:), intent(out)          :: btab
    logical, intent(out)                             :: warn
    real(rk), optional, intent(in)                   :: fillValue
    real(rk), dimension(:,:,:), optional, intent(in) :: precision
    include 'FilterValues_nd.f9h'
  end subroutine filterValues_REAL_3d

  subroutine filterValues_DOUBLE_3d(a, ATAB, b, BTAB, warn, fillValue, precision)
    ! Return arrays filtered etc.
    integer, parameter :: RK = kind(0.0d0) ! Kind type parameter for double precision
    ! Args
    real(rk), dimension(:,:,:), intent(in)           :: a
    real(rk), dimension(:,:,:), intent(in)           :: b
    real(rk), dimension(:,:,:), intent(out)          :: atab
    real(rk), dimension(:,:,:), intent(out)          :: btab
    logical, intent(out)                             :: warn
    real(rk), optional, intent(in)                   :: fillValue
    real(rk), dimension(:,:,:), optional, intent(in) :: precision
    include 'FilterValues_nd.f9h'
  end subroutine filterValues_DOUBLE_3d

! ------------------------------------------------- IsFillValue ---

  ! This family of routines checks to see if an arg is a fillValue
  elemental logical function IsFillValue_int ( A, FILLVALUE, STRICT )
    integer, intent(in) :: A
    integer, intent(in), optional :: FILLVALUE
    logical, optional, intent(in) :: STRICT ! true-> check FILLVALUE against arg
    ! Executable
    ! The following somewhat confusing business short-circuits
    ! checking the arg if we are strict but FillValue is not present
    IsFillValue_int = .true.
    if ( present(strict) ) IsFillValue_int = .not. ( &
      & strict .and. present(FillValue) &
      & )
    if ( .not. IsFillValue_int ) return
    if ( .not. present(FillValue) ) then
      IsFillValue_int = is_what_ieee( fill_signal, a )
    else
      IsFillValue_int = &
        & abs(a - FillValue) < 1 ! FILLVALUETOLERANCE
    endif
  end function IsFillValue_int

  elemental logical function IsFillValue_r4 ( A, FILLVALUE, STRICT )
    real(r4), intent(in) :: A
    real(r4), intent(in), optional :: FILLVALUE
    logical, optional, intent(in) :: STRICT ! true-> check FILLVALUE against arg
    ! Executable
    ! The following somewhat confusing business short-circuits
    ! checking the arg if we are strict but FillValue is not present
    IsFillValue_r4 = .true.
    if ( present(strict) ) IsFillValue_r4 = .not. ( &
      & strict .and. present(FillValue) &
      & )
    if ( .not. IsFillValue_r4 ) return
    if ( .not. present(FillValue) ) then
      IsFillValue_r4 = is_what_ieee( fill_signal, a )
    else
      IsFillValue_r4 = &
        & abs(a - FillValue) < max( FILLVALUETOLERANCE, abs(FillValue/100000) )
    endif
  end function IsFillValue_r4

  elemental logical function IsFillValue_r8 ( A, FILLVALUE, STRICT )
    real(r8), intent(in) :: A
    real(r8), intent(in), optional :: FILLVALUE
    logical, optional, intent(in) :: STRICT ! true-> check FILLVALUE against arg
    ! Executable
    ! The following somewhat confusing business short-circuits
    ! checking the arg if we are strict but FillValue is not present
    IsFillValue_r8 = .true.
    if ( present(strict) ) IsFillValue_r8 = .not. ( &
      & strict .and. present(FillValue) &
      & )
    if ( .not. IsFillValue_r8 ) return
    if ( .not. present(FillValue) ) then
      IsFillValue_r8 = is_what_ieee( fill_signal, a )
    else
      IsFillValue_r8 = &
        & abs(a - FillValue) < &
        & max( Real(FILLVALUETOLERANCE, r8), abs(FillValue/100000) )
    endif
  end function IsFillValue_r8

! ------------------------------------------------- IsFinite ---

  ! This family of routines checks to see if an arg is finite
  elemental logical function IsFinite_INTEGER ( A ) result( finite )
    integer, intent(in) :: A
    finite = .true.
  end function isfinite_INTEGER

! ------------------------------------------------- IsInfinite ---

  ! This family of routines checks to see if an arg is Infinite
  elemental logical function IsInfinite_REAL ( A ) result( Infinite )
    real, intent(in) :: A
    Infinite = .not. ( ieee_is_finite(a) .or. ieee_is_NaN(a) )
  end function isInfinite_real
  elemental logical function IsInfinite_DOUBLE ( A ) result( Infinite )
    double precision, intent(in) :: A
    Infinite = .not. ( ieee_is_finite(a) .or. ieee_is_NaN(a) )
  end function isInfinite_DOUBLE
  elemental logical function IsInfinite_INTEGER ( A ) result( Infinite )
    integer, intent(in) :: A
    Infinite = .false.
  end function isInfinite_INTEGER

! ------------------------------------------------- IsNaN ---

  ! This family of routines checks to see if an arg is NaN
  elemental logical function IsNaN_INTEGER ( A ) result( NaN )
    integer, intent(in) :: A
    NaN = .false.
  end function isNaN_INTEGER

  elemental logical function IsNaN_CHARACTER ( A ) result( NaN )
    character(len=8), intent(in) :: A
    NaN = .false.
  end function IsNaN_CHARACTER

  ! This family of subroutines makes arrays of values monotonically (increasing)
  ! May be used instead of BridgeMissing Values
  ! if Fill Values are -999.99
  ! Will also handle Periodic arrays, e.g. angles
  ! optional arg strict will insist that values[i] < values[i+1]
  subroutine Monotonize_1dint( values, Period, FillValue, strict )
    ! Args
    integer, dimension(:), intent(inout) :: values
    integer, optional, intent(in) :: Period
    integer, optional, intent(in) :: FillValue
    logical, optional, intent(in) :: strict
    ! Internal variables
    integer :: dx
    integer :: dxmin
    integer :: x1
    integer :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dint

  subroutine Monotonize_1dr4( values, Period, FillValue, strict )
    ! Args
    real(r4), dimension(:), intent(inout) :: values
    real(r4), optional, intent(in) :: Period
    real(r4), optional, intent(in) :: FillValue
    logical, optional, intent(in) :: strict
    ! Internal variables
    real(r4) :: dx
    real(r4) :: dxmin
    real(r4) :: x1
    real(r4) :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dr4

  subroutine Monotonize_1dr8( values, Period, FillValue, strict )
    ! Args
    real(r8), dimension(:), intent(inout) :: values
    real(r8), optional, intent(in) :: Period
    real(r8), optional, intent(in) :: FillValue
    logical, optional, intent(in) :: strict
    ! Internal variables
    real(r8) :: dx
    real(r8) :: dxmin
    real(r8) :: x1
    real(r8) :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dr8

  ! (But what about Period and FillValue? Don't you want them, too?)

  ! Actually, what doess it mean to be monotonic in 2d or 3d?
  ! The gradient is a vector, do you want to constrain it?
  ! Where will you need this? Will you ever need this?
  subroutine Monotonize_2dint(values)
    ! Args
    integer, dimension(:,:), intent(inout) :: values
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call Monotonize(values(i,:))
    enddo
  end subroutine Monotonize_2dint

  subroutine Monotonize_2dr4(values)
    ! Args
    real(r4), dimension(:,:), intent(inout) :: values
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call Monotonize(values(i,:))
    enddo
  end subroutine Monotonize_2dr4

  subroutine Monotonize_2dr8(values)
    ! Args
    real(r8), dimension(:,:), intent(inout) :: values
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call Monotonize(values(i,:))
    enddo
  end subroutine Monotonize_2dr8

  subroutine Monotonize_3dint(values)
    ! Args
    integer, dimension(:,:,:), intent(inout) :: values
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call Monotonize(values(i,j,:))
      enddo
    enddo
  end subroutine Monotonize_3dint

  subroutine Monotonize_3dr4(values)
    ! Args
    real(r4), dimension(:,:,:), intent(inout) :: values
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call Monotonize(values(i,j,:))
      enddo
    enddo
  end subroutine Monotonize_3dr4

  subroutine Monotonize_3dr8(values)
    ! Args
    real(r8), dimension(:,:,:), intent(inout) :: values
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call Monotonize(values(i,j,:))
      enddo
    enddo
  end subroutine Monotonize_3dr8

! -------------------------------------  RemoveFillValues  -----

  ! This family of routines removes fill values from an array
  ! returning a new, possibly smaller array
  ! If a second array of the same size is supplied
  ! it is reduced similarly by removing corresponding
  ! elements
  ! The new arrays must be large enough to accomodate non-fill elements
  ! 
  ! This should be a standard way of dealing with fill values:
  ! Say we read in a set of GMAO temperatures and their associated heights
  ! Some of the temperatures are fill values
  ! Before using the arrays to calculate tropopause pressures, we
  ! wish to remove all fill values (whose inclusion would sink the calculation)
  ! So count how many non-fill temperatures we've got, allocate new arrays
  ! to hold temperatures and heights, then call RemoveFillValues

  ! See also ReplaceFillValues

  subroutine RemoveFill1d_int ( array, FillValue, newArray, second, newSecond )
    integer, dimension(:), intent(in) :: array
    integer, intent(in) :: FillValue
    integer, dimension(:) :: newArray
    integer, dimension(:), optional, intent(in) :: second
    integer, dimension(:), optional :: newSecond
    !
    ! Local variables
    ! More local variables and executable
    include 'RemoveFillValues.f9h'
  end subroutine RemoveFill1d_int

  subroutine RemoveFill1d_r4 ( array, FillValue, newArray, second, newSecond )
    real(r4), dimension(:), intent(in) :: array
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:) :: newArray
    real(r4), dimension(:), optional, intent(in) :: second
    real(r4), dimension(:), optional :: newSecond
    !
    ! Local variables
    ! More local variables and executable
    include 'RemoveFillValues.f9h'
  end subroutine RemoveFill1d_r4

  subroutine RemoveFill1d_r8 ( array, FillValue, newArray, second, newSecond )
    real(r8), dimension(:), intent(in) :: array
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:) :: newArray
    real(r8), dimension(:), optional, intent(in) :: second
    real(r8), dimension(:), optional :: newSecond
    !
    ! Local variables
    ! More local variables and executable
    include 'RemoveFillValues.f9h'
  end subroutine RemoveFill1d_r8

  subroutine RemoveFill2d_int ( array, FillValue, newArray, second, newSecond )
    integer, dimension(:,:), intent(in) :: array
    integer, intent(in) :: FillValue
    integer, dimension(:) :: newArray
    integer, dimension(:,:), optional, intent(in) :: second
    integer, dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 2
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill2d_int

  subroutine RemoveFill2d_r4 ( array, FillValue, newArray, second, newSecond )
    real(r4), dimension(:,:), intent(in) :: array
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:) :: newArray
    real(r4), dimension(:,:), optional, intent(in) :: second
    real(r4), dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 2
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill2d_r4

  subroutine RemoveFill2d_r8 ( array, FillValue, newArray, second, newSecond )
    real(r8), dimension(:,:), intent(in) :: array
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:) :: newArray
    real(r8), dimension(:,:), optional, intent(in) :: second
    real(r8), dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 2
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill2d_r8

  subroutine RemoveFill3d_int ( array, FillValue, newArray, second, newSecond )
    integer, dimension(:,:,:), intent(in) :: array
    integer, intent(in) :: FillValue
    integer, dimension(:) :: newArray
    integer, dimension(:,:,:), optional, intent(in) :: second
    integer, dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 3
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill3d_int

  subroutine RemoveFill3d_r4 ( array, FillValue, newArray, second, newSecond )
    real(r4), dimension(:,:,:), intent(in) :: array
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:) :: newArray
    real(r4), dimension(:,:,:), optional, intent(in) :: second
    real(r4), dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 3
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill3d_r4

  subroutine RemoveFill3d_r8 ( array, FillValue, newArray, second, newSecond )
    real(r8), dimension(:,:,:), intent(in) :: array
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:) :: newArray
    real(r8), dimension(:,:,:), optional, intent(in) :: second
    real(r8), dimension(:), optional :: newSecond
    !
    ! Local variables
    integer, parameter       :: rank = 3
    integer, dimension(rank) :: shp
    ! More local variables and executable
    shp = shape(array)
    if ( .not. present(second) .or. .not. present(newSecond) ) then
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray )
    else
      call RemoveFillValues( reshape(array, (/product(shp)/) ), &
        & FillValue, newArray, reshape(second, (/product(shp)/) ), newSecond )
      endif
  end subroutine RemoveFill3d_r8

! -------------------------------------  ReorderFillValues  -----

  ! This family of routines reorders fillvalues in an array so they
  ! appear at the end
  ! E.g., given (/ 1, 4, -999, 6, -999, -1 /)
  ! returns (/ 1, 4, 6, -1, -999, -999 /)

  subroutine ReorderFillValues_r4 ( values, FillValue )
    real(r4), dimension(:), intent(inout) :: values
    real(r4), intent(in)                  :: FillValue
    !
    ! Local variables
    ! More local variables and executable
    include 'ReorderFillValues.f9h'
  end subroutine ReorderFillValues_r4

  subroutine ReorderFillValues_r8 ( values, FillValue )
    real(r8), dimension(:), intent(inout) :: values
    real(r8), intent(in)                  :: FillValue
    !
    ! Local variables
    ! More local variables and executable
    include 'ReorderFillValues.f9h'
  end subroutine ReorderFillValues_r8

  subroutine ReorderFillValues_int ( values, FillValue )
    integer, dimension(:), intent(inout) :: values
    integer, intent(in)                  :: FillValue
    !
    ! Local variables
    ! More local variables and executable
    include 'ReorderFillValues.f9h'
  end subroutine ReorderFillValues_int

! -------------------------------------  ReplaceFillValues  -----

  ! This family of routines replaces entries in an array
  ! based on whether they
  ! (1) are equal to FillValue; or
  ! (2) other criteria set by options
  ! The replacement values are supplied either by 
  ! newvalues, newFill, or according to options (e.g., you may interpolate)
  ! In short, you ensure result is positive-definite, negative definite,
  ! monotonic, or an analogous condition
  ! Note:
  ! When interpolating arrays with rank > 1, the interpolated-against
  ! index is the last one
  ! Thus we don't do true multi-dimensional interpolation
  ! If you wish to interpolate against another index, you must reshape
  
  ! This should be a standard way of dealing with fill values:
  ! Say you have an array of results, some of which are fillValues
  ! and also some of which are unphysical
  ! (Perhaps the calculations were invalid for some points)
  ! Before passing the array on for extra calculations you
  ! want to ensure all its entries will be valid (to keep next step
  ! from bombing)

  ! Can also be used to replace one set of fill values (e.g. 10^12)
  ! with another (e.g., -999.99)

  ! See also RemoveFillValues

  subroutine ReplaceFill1d_int ( values, FillValue, newValues, newFill, options )
    integer, dimension(:), intent(inout) :: values
    integer, intent(in) :: FillValue
    integer, dimension(:), optional, intent(in) :: newvalues
    integer, optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill1d_int

  subroutine ReplaceFill1d_r4 ( values, FillValue, newValues, newFill, options )
    real(r4), dimension(:), intent(inout) :: values
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:), optional, intent(in) :: newvalues
    real(r4), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill1d_r4

  subroutine ReplaceFill1d_r8 ( values, FillValue, newValues, newFill, options )
    real(r8), dimension(:), intent(inout) :: values
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:), optional, intent(in) :: newvalues
    real(r8), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill1d_r8

  subroutine ReplaceFill2d_int ( values, FillValue, newValues, newFill, options )
    integer, dimension(:, :), intent(inout) :: values
    integer, intent(in) :: FillValue
    integer, dimension(:, :), optional, intent(in) :: newvalues
    integer, optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill2d_int

  subroutine ReplaceFill2d_r4 ( values, FillValue, newValues, newFill, options )
    real(r4), dimension(:, :), intent(inout) :: values
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:, :), optional, intent(in) :: newvalues
    real(r4), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill2d_r4

  subroutine ReplaceFill2d_r8 ( values, FillValue, newValues, newFill, options )
    real(r8), dimension(:, :), intent(inout) :: values
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:, :), optional, intent(in) :: newvalues
    real(r8), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill2d_r8

  subroutine ReplaceFill3d_int ( values, FillValue, newValues, newFill, options )
    integer, dimension(:, :, :), intent(inout) :: values
    integer, intent(in) :: FillValue
    integer, dimension(:, :, :), optional, intent(in) :: newvalues
    integer, optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill3d_int

  subroutine ReplaceFill3d_r4 ( values, FillValue, newValues, newFill, options )
    real(r4), dimension(:, :, :), intent(inout) :: values
    real(r4), intent(in) :: FillValue
    real(r4), dimension(:, :, :), optional, intent(in) :: newvalues
    real(r4), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill3d_r4

  subroutine ReplaceFill3d_r8 ( values, FillValue, newValues, newFill, options )
    real(r8), dimension(:, :, :), intent(inout) :: values
    real(r8), intent(in) :: FillValue
    real(r8), dimension(:, :, :), optional, intent(in) :: newvalues
    real(r8), optional, intent(in) :: newFill
    character(len=*), optional, intent(in) :: options
    !
    ! Local variables
    ! More local variables and executable
    include 'ReplaceFillValues.f9h'
  end subroutine ReplaceFill3d_r8

! ============================================================================
  ! This family of subroutines bridges missing values by interpolation
  subroutine BridgeMissingValues_1dint(values, MissingValue)
    ! Args
    integer, dimension(:), intent(inout) :: values
    integer, intent(in), optional        :: missingValue
    ! Internal variables
    integer :: dx
    integer :: myMissingValue
    integer :: x1
    integer :: x2
    include 'BridgeMissingValues.f9h'
  end subroutine BridgeMissingValues_1dint

  subroutine BridgeMissingValues_1dr4(values, MissingValue)
    ! Args
    real(r4), dimension(:), intent(inout) :: values
    real(r4), intent(in), optional        :: missingValue
    ! Internal variables
    real(r4) :: dx
    real(r4) :: myMissingValue
    real(r4) :: x1
    real(r4) :: x2
    include 'BridgeMissingValues.f9h'
  end subroutine BridgeMissingValues_1dr4

  subroutine BridgeMissingValues_1dr8(values, MissingValue)
    ! Args
    real(r8), dimension(:), intent(inout) :: values
    real(r8), intent(in), optional        :: missingValue
    ! Internal variables
    real(r8) :: dx
    real(r8) :: myMissingValue
    real(r8) :: x1
    real(r8) :: x2
    include 'BridgeMissingValues.f9h'
  end subroutine BridgeMissingValues_1dr8

  subroutine BridgeMissingValues_2dint(values, MissingValue)
    ! Args
    integer, dimension(:,:), intent(inout) :: values
    integer, intent(in), optional          :: missingValue
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call BridgeMissingValues(values(i,:), MissingValue)
    enddo
  end subroutine BridgeMissingValues_2dint

  subroutine BridgeMissingValues_2dr4(values, MissingValue)
    ! Args
    real(r4), dimension(:,:), intent(inout) :: values
    real(r4), intent(in), optional          :: missingValue
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call BridgeMissingValues(values(i,:), MissingValue)
    enddo
  end subroutine BridgeMissingValues_2dr4

  subroutine BridgeMissingValues_2dr8(values, MissingValue)
    ! Args
    real(r8), dimension(:,:), intent(inout) :: values
    real(r8), intent(in), optional          :: missingValue
    ! Internal variables
    integer :: i
    do i=1, size(values, 1)
      call BridgeMissingValues(values(i,:), MissingValue)
    enddo
  end subroutine BridgeMissingValues_2dr8

  subroutine BridgeMissingValues_3dint(values, MissingValue)
    ! Args
    integer, dimension(:,:,:), intent(inout) :: values
    integer, intent(in), optional            :: missingValue
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call BridgeMissingValues(values(i,j,:), MissingValue)
      enddo
    enddo
  end subroutine BridgeMissingValues_3dint

  subroutine BridgeMissingValues_3dr4(values, MissingValue)
    ! Args
    real(r4), dimension(:,:,:), intent(inout) :: values
    real(r4), intent(in), optional            :: missingValue
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call BridgeMissingValues(values(i,j,:), MissingValue)
      enddo
    enddo
  end subroutine BridgeMissingValues_3dr4

  subroutine BridgeMissingValues_3dr8(values, MissingValue)
    ! Args
    real(r8), dimension(:,:,:), intent(inout) :: values
    real(r8), intent(in), optional            :: missingValue
    ! Internal variables
    integer :: i, j
    do j=1, size(values, 2)
      do i=1, size(values, 1)
        call BridgeMissingValues(values(i,j,:), MissingValue)
      enddo
    enddo
  end subroutine BridgeMissingValues_3dr8

  ! ----------------------------------  roundUpOrDown  -----
  ! This family of functions returns a result rounded up if
  ! the fractional part > .5
  ! otherwise rounded down
  elemental function roundUpOrDown_r4( value ) result ( rounded )
    integer, parameter :: RK = R4
    ! Args
    real(rk), intent(in) :: value
    real(rk) :: rounded
    ! Internal variables
    real(rk) :: frac
    integer :: sgn
    ! Executable
    rounded = int( value )
    sgn = sign(1._rk, value )
    frac = sgn * ( value - rounded )
    if ( frac > 0.5_rk ) rounded = rounded + sgn
  end function roundUpOrDown_r4

  elemental function roundUpOrDown_r8( value ) result ( rounded )
    integer, parameter :: RK = R8
    ! Args
    real(rk), intent(in) :: value
    real(rk) :: rounded
    ! Internal variables
    real(rk) :: frac
    integer :: sgn
    ! Executable
    rounded = int( value )
    sgn = sign(1._rk, value )
    frac = sgn * ( value - rounded )
    if ( frac > 0.5_rk ) rounded = rounded + sgn
  end function roundUpOrDown_r8

  elemental function roundUpOrDown_int( value ) result ( rounded )
    ! Trivial; supplied just to satisfy other generic functions
    ! Args
    integer, intent(in) :: value
    integer :: rounded
    ! Executable
    rounded = value
  end function roundUpOrDown_int

  ! ----------------------------------  WhereAreTheFills  -----
  ! This family of procedures finds which array elements are Fills
  subroutine WhereAreTheFills_DOUBLE ( array, which, howMany )
    ! Args
    double precision, dimension(:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( Fill_signal, array, which, howMany )
  end subroutine WhereAreTheFills_DOUBLE

  subroutine WhereAreTheFills_integer ( array, which, howMany )
    ! Args
    integer, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( Fill_signal, array, which, howMany )
  end subroutine WhereAreTheFills_integer

  subroutine WhereAreTheFills_REAL ( array, which, howMany )
    ! Args
    real, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( Fill_signal, array, which, howMany )
  end subroutine WhereAreTheFills_REAL

  subroutine WhereAreTheFills_DOUBLE_2d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_DOUBLE_2d

  subroutine WhereAreTheFills_integer_2d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_integer_2d

  subroutine WhereAreTheFills_REAL_2d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_REAL_2d

  subroutine WhereAreTheFills_DOUBLE_3d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_DOUBLE_3d

  subroutine WhereAreTheFills_integer_3d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_integer_3d

  subroutine WhereAreTheFills_REAL_3d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( Fill_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheFills_REAL_3d

  ! ----------------------------------  WhereAreTheInfs  -----
  ! This family of procedures finds which array elements are Infs
  subroutine WhereAreTheInfs_DOUBLE ( array, which, howMany )
    ! Args
    double precision, dimension(:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( inf_signal, array, which, howMany )
  end subroutine WhereAreTheInfs_DOUBLE

  subroutine WhereAreTheInfs_integer ( array, which, howMany )
    ! Args
    integer, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( inf_signal, array, which, howMany )
  end subroutine WhereAreTheInfs_integer

  subroutine WhereAreTheInfs_REAL ( array, which, howMany )
    ! Args
    real, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( inf_signal, array, which, howMany )
  end subroutine WhereAreTheInfs_REAL

  subroutine WhereAreTheInfs_DOUBLE_2d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_DOUBLE_2d

  subroutine WhereAreTheInfs_integer_2d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_integer_2d

  subroutine WhereAreTheInfs_REAL_2d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_REAL_2d

  subroutine WhereAreTheInfs_DOUBLE_3d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_DOUBLE_3d

  subroutine WhereAreTheInfs_DOUBLE_4d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_DOUBLE_4d

  subroutine WhereAreTheInfs_integer_3d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_integer_3d

  subroutine WhereAreTheInfs_REAL_3d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_REAL_3d

  subroutine WhereAreTheInfs_REAL_4d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( inf_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheInfs_REAL_4d

  ! ----------------------------------  WhereAreTheNaNs  -----
  ! This family of procedures finds which array elements are NaNs
  subroutine WhereAreTheNaNs_DOUBLE ( array, which, howMany )
    ! Args
    double precision, dimension(:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( NaN_signal, array, which, howMany )
  end subroutine WhereAreTheNaNs_DOUBLE

  subroutine WhereAreTheNaNs_integer ( array, which, howMany )
    ! Args
    integer, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( NaN_signal, array, which, howMany )
  end subroutine WhereAreTheNaNs_integer

  subroutine WhereAreTheNaNs_REAL ( array, which, howMany )
    ! Args
    real, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    call WhereAreThey( NaN_signal, array, which, howMany )
  end subroutine WhereAreTheNaNs_REAL

  subroutine WhereAreTheNaNs_DOUBLE_2d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_DOUBLE_2d

  subroutine WhereAreTheNaNs_integer_2d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_integer_2d

  subroutine WhereAreTheNaNs_REAL_2d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_REAL_2d

  subroutine WhereAreTheNaNs_DOUBLE_3d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:,:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_DOUBLE_3d

  subroutine WhereAreTheNaNs_DOUBLE_4d ( array, which, howMany, mode, inds )
    ! Args
    double precision, dimension(:,:,:,:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_DOUBLE_4d

  subroutine WhereAreTheNaNs_integer_3d ( array, which, howMany, mode, inds )
    ! Args
    integer, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_integer_3d

  subroutine WhereAreTheNaNs_REAL_3d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_REAL_3d

  subroutine WhereAreTheNaNs_REAL_4d ( array, which, howMany, mode, inds )
    ! Args
    real, dimension(:,:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    call WhereAreThey( NaN_signal, array, which, howMany, mode, inds )
  end subroutine WhereAreTheNaNs_REAL_4d

  ! ----------------------------------  private procedures  -----
  subroutine announce_error(message, int1, int2, dontstop)
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: int1
    integer, optional, intent(in) :: int2
    logical, optional, intent(in) :: dontstop
    logical :: keepgoing
    !
    keepgoing = .false.
    if ( present(dontstop) ) keepgoing=dontstop
    if ( .not. keepgoing ) then
      call output('*** Error in MLSFillValues module ***', advance='yes')
    endif
    call output(trim(message), advance='no')
    call blanks(3)
    if ( present(int1) ) write(*,'(i4)',advance='no') int1
    call blanks(3)
    if ( present(int2) ) write(*,'(i4)', advance='no') int2
    if ( .not. keepgoing ) stop
  end subroutine announce_error

  subroutine swap_int( first, second )
    ! Swap first and second args
    integer, intent(inout) :: first
    integer, intent(inout) :: second
    integer :: temp
    ! Executable
    temp = first
    first = second
    second = temp
  end subroutine swap_int

  subroutine swap_r4( first, second )
    ! Swap first and second args
    real(r4), intent(inout) :: first
    real(r4), intent(inout) :: second
    real(r4) :: temp
    ! Executable
    temp = first
    first = second
    second = temp
  end subroutine swap_r4

  subroutine swap_r8( first, second )
    ! Swap first and second args
    real(r8), intent(inout) :: first
    real(r8), intent(inout) :: second
    real(r8) :: temp
    ! Executable
    temp = first
    first = second
    second = temp
  end subroutine swap_r8

  ! ----------------------------------  WhereAreThey  -----
  ! This family of procedures finds which array elements are whats
  ! where what is one of our ieee signalling flags
  subroutine WhereAreThey_DOUBLE ( what, array, which, howMany )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    double precision, dimension(:), intent(in)   :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    ! Internal variables
    integer :: i
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array) < 1 ) return
    n = 0
    do i=1, size(array)
      if ( is_what_ieee(what, array(i)) ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_DOUBLE

  subroutine WhereAreThey_integer ( what, array, which, howMany )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    integer, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    ! Internal variables
    integer :: i
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array) < 1 ) return
    n = 0
    do i=1, size(array)
      if ( is_what_ieee(what, array(i)) ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_integer

  subroutine WhereAreThey_REAL ( what, array, which, howMany )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    real, dimension(:), intent(in)               :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    ! Internal variables
    integer :: i
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array) < 1 ) return
    n = 0
    do i=1, size(array)
      if ( is_what_ieee(what, array(i)) ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_REAL

  subroutine WhereAreThey_DOUBLE_2d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    double precision, dimension(:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,2) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    yes = .false.
    do i=1, size(array, 2)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 1)
          if ( .not. is_what_ieee(what, array(j,i)) ) cycle
          n = n + 1
          if ( present(inds) ) inds( 1:2, min(n,size(inds,2)) ) = (/ j, i /)
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_DOUBLE_2d

  subroutine WhereAreThey_integer_2d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    integer, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,2) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 2)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 1)
          if ( .not. is_what_ieee(what, array(j,i)) ) cycle
          n = n + 1
          if ( present(inds) ) inds( 1:2, min(n,size(inds,2)) ) = (/ j, i /)
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_integer_2d

  subroutine WhereAreThey_REAL_2d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    real, dimension(:,:), intent(in)             :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,2) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 2)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 1)
          if ( .not. is_what_ieee(what, array(j,i)) ) cycle
          n = n + 1
          if ( present(inds) ) inds( 1:2, min(n,size(inds,2)) ) = (/ j, i /)
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_REAL_2d

  subroutine WhereAreThey_DOUBLE_3d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    double precision, dimension(:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j, k
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,3) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 3)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 2)
          do k=1, size(array, 1)
            if ( .not. is_what_ieee(what, array(k, j,i)) ) cycle
            n = n + 1
          if ( present(inds) ) inds( 1:3, min(n,size(inds,2)) ) = (/ k, j, i /)
          enddo
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_DOUBLE_3d

  subroutine WhereAreThey_DOUBLE_4d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    double precision, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j, k, l
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,4) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 4)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 3)
          do k=1, size(array, 2)
            do l=1, size(array, 1)
              if ( .not. is_what_ieee(what, array(l, k, j, i)) ) cycle
              n = n + 1
            if ( present(inds) ) inds( 1:4, min(n,size(inds,2)) ) = (/ l, k, j, i /)
            enddo
          enddo
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,:,:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,:,:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_DOUBLE_4d

  subroutine WhereAreThey_integer_3d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    integer, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j, k
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,3) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 3)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 2)
          do k=1, size(array, 1)
            if ( .not. is_what_ieee(what, array(k, j,i)) ) cycle
            n = n + 1
          if ( present(inds) ) inds( 1:3, min(n,size(inds,2)) ) = (/ k, j, i /)
          enddo
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_integer_3d

  subroutine WhereAreThey_REAL_3d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    real, dimension(:,:,:), intent(in)           :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j, k
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,3) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 3)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 2)
          do k=1, size(array, 1)
            if ( .not. is_what_ieee(what, array(k, j,i)) ) cycle
            n = n + 1
          if ( present(inds) ) inds( 1:3, min(n,size(inds,2)) ) = (/ k, j, i /)
          enddo
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_REAL_3d


  subroutine WhereAreThey_real_4d ( what, array, which, howMany, mode, inds )
    ! Args
    integer, intent(in)                          :: what ! a signal flag
    real, dimension(:,:,:,:), intent(in) :: array
    integer, dimension(:), optional, intent(out) :: which
    integer, optional, intent(out)               :: howMany
    character(len=*), optional, intent(in)       :: mode ! any, ind, or all
    integer, dimension(:,:), optional, intent(out) :: inds
    ! Internal variables
    integer :: i, j, k, l
    character(len=3) :: myMode
    integer :: n
    integer :: nWhich ! size(which);  > 0 only if which is present
    logical :: yes
    ! Executable
    nWhich = 0
    if ( present(which) ) nWhich = size(which)
    if ( nWhich > 0 ) which = 0
    if ( present(howMany) ) howMany = 0
    if ( size(array,4) < 1 ) return
    myMode = 'any'
    if ( present(mode) ) myMode = mode
    n = 0
    do i=1, size(array, 4)
      select case( lowercase(myMode(1:3)) )
      case ( 'ind' )
        do j=1, size(array, 3)
          do k=1, size(array, 2)
            do l=1, size(array, 1)
              if ( .not. is_what_ieee(what, array(l, k, j, i)) ) cycle
              n = n + 1
            if ( present(inds) ) inds( 1:4, min(n,size(inds,2)) ) = (/ l, k, j, i /)
            enddo
          enddo
        enddo
      case ( 'any' )
        yes = any( is_what_ieee(what, array(:,:,:,i)) )
      case ( 'all' )
        yes = all( is_what_ieee(what, array(:,:,:,i)) )
      case default ! unrecognized mode; sorry
        yes = .false.
      end select
      if ( yes ) then
        n = n + 1
        if ( n <= nWhich ) which(n) = i
      endif
    enddo
    if ( present(howMany) ) howMany = n
  end subroutine WhereAreThey_real_4d

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSFillValues
!=============================================================================

!
! $Log$
! Revision 2.39  2018/04/19 02:00:36  vsnyder
! Compute address for allocate/deallocate tracking.  Remove USE statements for
! unused names.
!
! Revision 2.38  2017/11/03 19:56:00  pwagner
! Most array gymnastics moved to HyperSlabs module
!
! Revision 2.37  2017/11/02 00:06:50  pwagner
! rerank uses Subscripts from Array_Stuff now
!
! Revision 2.36  2017/08/25 00:18:10  pwagner
! Dropped all the internal output routines in favor of Output_m
!
! Revision 2.35  2016/01/13 00:45:44  pwagner
! May optionally stutter instead of interpolate
!
! Revision 2.34  2015/03/28 01:49:59  vsnyder
! Moved IsMonotonic to Monotone module
!
! Revision 2.33  2014/09/05 00:04:34  vsnyder
! More complete and accurate allocate/deallocate size tracking.  Remove
! unnecessary POINTER attribute from some arguments.  Remove declarations
! of unused variables.
!
! Revision 2.32  2013/09/17 22:35:13  pwagner
! Changed api of Embed, Extract arrays to match hyperslab
!
! Revision 2.31  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.30  2011/12/13 01:07:36  pwagner
! Now uses MLSFills if allocated; can add to or dump it
!
! Revision 2.29  2011/12/07 01:16:15  pwagner
! Added Bandwidth calculation
!
! Revision 2.28  2011/07/26 20:44:21  pwagner
! Added some 4d interfaces
!
! Revision 2.27  2011/07/07 00:28:27  pwagner
! Made extremum elemental
!
! Revision 2.26  2011/03/22 23:38:30  pwagner
! Rerank now public
!
! Revision 2.25  2011/01/20 01:14:31  pwagner
! Added Collapse: Turns array into lower-rank representation, collapsing last index
!
! Revision 2.24  2010/11/30 00:33:52  pwagner
! Added instance for character arg to isNaN
!
! Revision 2.23  2010/11/11 19:54:29  pwagner
! May select values by relative size
!
! Revision 2.22  2010/10/13 00:41:23  pwagner
! WhereAreThe.. can now take'ind' mode, returning indexes of each
!
! Revision 2.21  2010/09/24 23:44:48  pwagner
! removed all but integervalued isFinite, isNaN in favor of ieee_arithmetic-supplied versions
!
! Revision 2.20  2010/08/13 22:02:41  pwagner
! Renamed Repopulate; renamed decimate to Depopulate
!
! Revision 2.19  2010/02/17 22:30:53  pwagner
! Added Depopulate routines to pick out non-zeros in sparse arrays
!
! Revision 2.18  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.17  2008/09/03 20:41:17  pwagner
! Added GatherArray
!
! Revision 2.16  2008/07/10 00:12:43  pwagner
! Added extremum, halfWaves
!
! Revision 2.15  2008/06/06 22:52:21  pwagner
! EssentiallyEqual moved to MLSFillValues
!
! Revision 2.14  2008/06/04 21:42:56  pwagner
! Monotonize now takes optional arg strict
!
! Revision 2.13  2008/04/10 20:25:50  pwagner
! Montonize can take optional arg FillValue
!
! Revision 2.12  2008/03/07 01:33:47  pwagner
! Added optional arg Period to Monotonize
!
! Revision 2.11  2008/01/09 20:50:47  pwagner
! WhereAreThe.. support integer arrays
!
! Revision 2.10  2008/01/07 21:34:41  pwagner
! Added new functions to id and return ieee signals; WhereAreThe.. procedures
!
! Revision 2.9  2006/06/15 17:31:30  pwagner
! Added ReorderFillValues
!
! Revision 2.8  2006/03/15 17:32:35  pwagner
! Can removeFillValues from multi-dimensional arrays if results are rank 1
!
! Revision 2.7  2006/02/28 21:43:31  pwagner
! Improve comments regarding FilterValues
!
! Revision 2.6  2006/02/02 16:19:57  pwagner
! Should not use optional arg FillValue unless present
!
! Revision 2.5  2006/02/01 23:54:20  pwagner
! Fixed bug in integer format
!
! Revision 2.4  2006/02/01 23:42:14  pwagner
! Added RemoveFillValues
!
! Revision 2.3  2006/01/14 00:50:15  pwagner
! Added procedures to embed, extract blocs from larger arrays
!
! Revision 2.2  2005/12/23 03:10:31  vsnyder
! Make some routines more generic, using include
!
! Revision 2.1  2005/12/16 00:00:23  pwagner
! Created to hold fillValue-related stuff
