! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Diff_1

  ! Dump differences of arrays

  use HyperSlabs, only: HalfWaves
  use Dump_Options, only: AuBrick, &
    & DefaultPCTFormat, Direct, Dopts, DumpDumpOptions, DumpTableSide, &
    & NameOnEachLine, NameHasBeenPrinted, NaNs, PCTFormat, &
    & MyRatios=>Ratios, RMS, RMSFormat, Stats, Table, &
    & TheDumpBegins, TheDumpEnds, WholeArray
  use Dump_0, only: Dump, FinishLine, PrintName, PrintRMSEtc
  use Dump_1, only: DumpTable
  use HighOutput, only: OutputNamedValue
  use IEEE_Arithmetic, only: IEEE_Is_Finite
  use MLSFillValues, only: FilterValues, NaNFunction, &
    & ReorderFillValues, ReplaceFillValues, &
    & WhereAreTheInfs, WhereAreTheNaNs
  use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
  use MLSStats1, only: AllStats, HowFar, HowNear, MLSStdDev, Ratios, Reset, &
    & Stat_T
  use Output_m, only: OutputOptions, NewLine, Output
  use Time_M, only: Time_Now

  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)

!     (subroutines and functions)
! Diff_Fun           returns differences between scalars, arrays, etc.
! Diff               dump diffs between pair of arrays of numeric type
! SelfDiff           dump increments between successive array values
! === (end of toc) ===

! === (start of api) ===
! num diff_fun ( num value1, num value2, num auxvalue, char* options )
! diff ( array1, char* name1, array2, char* name2,
!      [fillvalue], [int width], [char* format],
!      [int lbound], [char* options] ) 
!       where array1, array2 can be 1, 2, or 3d arrays of 
!       ints, reals, or doubles, compatible in size and type
!       and fillValue is a scalar of the same type, if present

! The optional arguments FillValue, Width, Format, LBound and Options
! are described in the module Dump_0.
! === (end of api) ===

  public :: Diff, Diff_Fun, SelfDiff

  interface Diff            ! dump diffs between pair of n-d arrays of numeric type
    module procedure Diff_1D_Double, Diff_1D_Integer, Diff_1D_Real
    module procedure Diff_2D_Double, Diff_2D_Integer, Diff_2D_Real
    module procedure Diff_3D_Double, Diff_3D_Real
    module procedure Diff_4D_Double, Diff_4D_Real
  end interface

  interface Diff_Fun        ! return diffs between args or arrays of numeric type
    module procedure Diff_Scalar_Double, Diff_Scalar_Real
  end interface

  interface SelfDiff        ! dump increments between successive array values
    module procedure SelfDiff_Integer
    module procedure SelfDiff_Real
    module procedure SelfDiff_Double
  end interface

  ! =====     Private     ==============================================

  interface FilteredDiff    ! dump FilteredDiff_s between pair of n-d arrays of numeric type
    module procedure FilteredDiff_1D_Double, FilteredDiff_1D_Integer, FilteredDiff_1D_Real
    module procedure FilteredDiff_2D_Double, FilteredDiff_2D_Integer, FilteredDiff_2D_Real
    module procedure FilteredDiff_3D_Double, FilteredDiff_3D_Real
    module procedure FilteredDiff_4D_Double, FilteredDiff_4D_Real
  end interface

  interface UnfilteredDiff  ! dump UnfilteredDiffs between pair of n-d arrays of numeric type
    module procedure UnfilteredDiff_1D_Double, UnfilteredDiff_1D_Integer, UnfilteredDiff_1D_Real
    module procedure UnfilteredDiff_2D_Double, UnfilteredDiff_2D_Integer, UnfilteredDiff_2D_Real
    module procedure UnfilteredDiff_3D_Double, UnfilteredDiff_3D_Real
    module procedure UnfilteredDiff_4D_Double, UnfilteredDiff_4D_Real
  end interface

  ! Declared module-wide purely for convenience
  integer            :: How_Many
  character(len=16)  :: MyPCTFormat
  integer            :: NumNonFill
  logical, parameter :: ShortcutDiffs      = .false.
  logical, save      :: ThisIsADiff        = .false.
  integer, parameter :: TooManyElements    = 125*130*3500 ! Can diff l1b DACS
  integer, dimension(1024) :: Which

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

 ! -------------------------------------------------------  Diff_  -----
 ! This family of routines dumps the differences between two arrays
 ! with the same shape and numeric type
 ! Its behavior is modified by the following optional parameters
 ! FillValue   Ignore these values when computing min, max
 ! Width       Row size when dumping wholearrays
 ! Format      Output format when printing wholearray
 ! LBound      Lower bound when printing wholearray indices
 !
 ! The following optional params are now set by options
 ! Wholearray  Whether to print whole array of differences
 ! Stats       Dump number of differences found, %
 ! RMS         Dump min, max, rms
 ! Clean       Clean up after any prior dumps
  subroutine Diff_1D_Double ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: Array1(:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: Array2(:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_1D_Double

  subroutine Diff_1D_Integer ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    integer, intent(in) :: Array1(:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: Array2(:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FILLVALUE, Width, Format, LBound, Options )
    end if
  end subroutine Diff_1D_Integer

  subroutine Diff_1D_Real ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: Array1(:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: Array2(:)
    character(len=*), intent(in) :: Name2
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_1D_Real

  subroutine Diff_2D_Integer ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    integer, intent(in) :: Array1(:,:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: Array2(:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_2D_Integer

  subroutine Diff_2D_Double ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: Array1(:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: Array2(:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_2D_Double

  subroutine Diff_2D_Real ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: Array1(:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: Array2(:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_2D_Real

  subroutine Diff_3D_Double ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: Array1(:,:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: Array2(:,:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_3D_Double

  subroutine Diff_3D_Real ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: Array1(:,:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: Array2(:,:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_3D_Real

  subroutine Diff_4D_Double ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: Array1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: Array2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    call theDumpBegins ( options )
    if ( any(shape(array1) == 0) ) then
      call output( 'array sizes are 0', advance='yes' )
    else if ( size(array1,1) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '1st size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(1,:,:,:), name1, array2(1,:,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,2) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '2nd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,1,:,:), name1, array2(:,1,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,3) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '3rd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,1,:), name1, array2(:,:,1,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,4) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '4th size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,:,1), name1, array2(:,:,:,1), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( .not. present(fillValue) ) then
      if ( dopts(wholeArray)%v ) call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_4D_Double

  subroutine Diff_4D_Real ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: Array1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: Array2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    call theDumpBegins ( options )
    if ( any(shape(array1) == 0) ) then
      call output( 'array sizes are 0', advance='yes' )
    else if ( size(array1,1) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '1st size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(1,:,:,:), name1, array2(1,:,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,2) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '2nd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,1,:,:), name1, array2(:,1,:,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,3) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '3rd size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,1,:), name1, array2(:,:,1,:), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( size(array1,4) == 1 ) then
      if ( dopts(wholeArray)%v ) call output( '4th size is 1: reducing rank to 3', advance='yes' )
      call diff ( array1(:,:,:,1), name1, array2(:,:,:,1), name2, &
        & fillvalue, width, format, lbound, options )
    else if ( .not. present(fillValue) ) then
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else if ( product(shape(array1)) > tooManyElements ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'array size of ' // trim(name1) // ' too large to filter Fill values' )
      call UnfilteredDiff( Array1, Name1, Array2, Name2, &
        & Width, Format, LBound, Options )
    else
      call FilteredDiff( Array1, Name1, Array2, Name2, &
        & FillValue, Width, Format, LBound, Options )
    end if
  end subroutine Diff_4D_Real

  ! -----------------------------------------------  Diff_Scalar  -----
  ! This family of functions differences two values and returns
  ! according to options:
  ! if options contains
  ! character        meaning
  ! ---------        -------
  !     a            diff the absolute values
  !     r            divide the difference by the max abs of the 2 values
  !     f            treat auxvalue as a fillvalue: return fillvalue if either
  !                    value is fillvalue
  !     p            treat auxvalue as a period: return 
  !                    min(value1 - value2 + n*period)
  elemental function Diff_Scalar_Double ( VALUE1, VALUE2, AUXVALUE, Options ) &
    & result(d)
    ! Args
    double precision, intent(in)           :: Value1
    double precision, intent(in)           :: Value2
    double precision, optional, intent(in) :: Auxvalue
    character(len=*), optional, intent(in) :: Options
    double precision                       :: D
    ! Internal variables
    integer, parameter :: RK = kind(0.0d0)
    ! Executable
    include 'diff_scalar.f9h'
  end function Diff_Scalar_Double

  elemental function Diff_Scalar_Real ( VALUE1, VALUE2, AUXVALUE, Options ) &
    & result(d)
    ! Args
    real, intent(in)                       :: Value1
    real, intent(in)                       :: Value2
    real, optional, intent(in)             :: Auxvalue
    character(len=*), optional, intent(in) :: Options
    real                                   :: D
    ! Internal variables
    integer, parameter :: RK = kind(0.0e0)
    ! Executable
    include 'diff_scalar.f9h'
  end function Diff_Scalar_Real

  ! -----------------------------------------------  FilteredDiff  -----
  subroutine FilteredDiff_1D_Double ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: inArray1(:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: inArray2(:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    double precision, dimension(size(inArray1)) :: array1
    double precision, dimension(size(inArray2)) :: array2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_1D_Double

  subroutine FilteredDiff_1D_Integer ( IArray1, Name1, IArray2, Name2, &
    & IFillValue, Width, Format, LBound, Options )
    integer, intent(in) :: IArray1(:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: IArray2(:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: IFillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1)) :: array1
    real, dimension(size(iarray2)) :: array2
    real :: fillValue
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    fillValue = iFillValue
    call diff ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
  end subroutine FilteredDiff_1D_Integer

  subroutine FilteredDiff_1D_Real ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: inArray1(:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: inArray2(:)
    character(len=*), intent(in) :: Name2
    real, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(inArray1)) :: array1
    real, dimension(size(inArray2)) :: array2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_1D_Real

  subroutine FilteredDiff_2D_Double ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: inArray1(:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: inArray2(:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options
    !
    double precision, dimension(product(shape(inArray1))) :: array1
    double precision, dimension(product(shape(inArray2))) :: array2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_2D_Double

  subroutine FilteredDiff_2D_Integer ( IArray1, Name1, IArray2, Name2, &
    & IFillValue, Width, Format, LBound, Options )
    integer, intent(in) :: IArray1(:,:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: IArray2(:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: IFillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1,1), size(iarray1,2)) :: array1
    real, dimension(size(iarray1,1), size(iarray1,2)) :: array2
    real :: FillValue
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    fillValue = iFillValue
    call diff ( Array1, Name1, Array2, Name2, &
    & FillValue, Width, Format, LBound, Options )
  end subroutine FilteredDiff_2D_Integer

  subroutine FilteredDiff_2D_Real ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: inArray1(:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: inArray2(:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options
    !
    real, dimension(product(shape(inArray1))) :: array1
    real, dimension(product(shape(inArray2))) :: array2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_2D_Real

  subroutine FilteredDiff_3D_Double ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: inArray1(:,:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: inArray2(:,:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    double precision, dimension(product(shape(inArray1))) :: array1
    double precision, dimension(product(shape(inArray2))) :: array2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_3D_Double

  subroutine FilteredDiff_3D_Real ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: inArray1(:,:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: inArray2(:,:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    real, dimension(product(shape(inArray1))) :: array1
    real, dimension(product(shape(inArray2))) :: array2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_3D_Real

  subroutine FilteredDiff_4D_Double ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    double precision, intent(in) :: inArray1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    double precision, intent(in) :: inArray2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    double precision, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    double precision, dimension(product(shape(inArray1))) :: array1
    double precision, dimension(product(shape(inArray2))) :: array2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_4D_Double

  subroutine FilteredDiff_4D_Real ( inArray1, Name1, inArray2, Name2, &
    & FillValue, Width, Format, LBound, Options )
    real, intent(in) :: inArray1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    real, intent(in) :: inArray2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    real, intent(in):: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    real, dimension(product(shape(inArray1))) :: array1
    real, dimension(product(shape(inArray2))) :: array2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine FilteredDiff_4D_Real

 ! ---------------------------------------------  SelfDiff_Double  -----
  subroutine SelfDiff_Double ( Array, Name, &
    & FillValue, Width, Format, waves, LBound, Options )
    ! dump the running increment == ( array(i) - array(i-1) )
    double precision, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    double precision, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    double precision, dimension(size(array)-1) :: increment
    integer, dimension(100) :: lengths
    logical :: myWaves
    integer :: nWaves
    ! Executable
    myWaves = .false.
    if ( present(waves) ) myWaves = waves
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    end do
    if ( myWaves ) then
      call halfWaves( increment, lengths, nWaves )
      if ( nWaves > 0 ) then
        call dump( lengths(1:nWaves), 'half-waves of ' // Name )
      end if
      return
    end if
    call dump ( increment, Name, &
      & FillValue, Width, Format, LBound, Options )
  end subroutine SelfDiff_Double

 ! ---------------------------------------------  SelfDiff_Integer  -----
  subroutine SelfDiff_Integer ( Array, Name, &
    & FillValue, Format, Width, waves, LBound, Options )
    ! dump the running increment == ( array(i) - array(i-1) )
    integer, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    integer, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    integer, dimension(size(array)-1) :: increment
    ! Executable
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    end do
    call dump ( increment, Name, &
      & FillValue, Format, Width, LBound, Options )
  end subroutine SelfDiff_Integer

 ! ---------------------------------------------  SelfDiff_Real  -----
  subroutine SelfDiff_Real ( Array, Name, &
    & FillValue, Format, Width, waves, LBound, Options )
    ! dump the running increment == ( array(i) - array(i-1) )
    real, intent(in) :: Array(:)
    character(len=*), intent(in), optional :: Name
    real, intent(in), optional :: FillValue
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    logical, intent(in), optional :: WAVES
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options
    ! Local variables
    integer                                  :: i
    real, dimension(size(array)-1) :: increment
    ! logical, parameter :: unique = .false. 
    integer, dimension(100) :: lengths
    logical :: myWaves
    integer :: nWaves
    ! Executable
    myWaves = .false.
    if ( present(waves) ) myWaves = waves
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    end do
    if ( myWaves ) then
      call halfWaves( increment, lengths, nWaves )
      if ( nWaves > 0 ) then
        call dump( lengths(1:nWaves), 'half-waves of ' // Name )
      end if
      return
    end if
    call dump ( increment, Name, &
      & FillValue, Width, Format, LBound, Options )
  end subroutine SelfDiff_Real

  ! ---------------------------------------------  UnfilteredDiff  -----
  subroutine UnfilteredDiff_1D_Double ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: Array1(:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_1D_Double

  subroutine UnfilteredDiff_1D_Integer ( IArray1, Name1, IArray2, Name2, &
    & Width, Format, LBound, Options )
    integer, intent(in) :: IArray1(:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: IArray2(:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1)) :: array1
    real, dimension(size(iarray2)) :: array2
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    call diff ( Array1, Name1, Array2, Name2, &
      & Width=Width, Format=Format, &
      & LBound=LBound, Options=Options )
  end subroutine UnfilteredDiff_1D_Integer

  subroutine UnfilteredDiff_1D_Real ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Array1(:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_1D_Real

  subroutine UnfilteredDiff_2D_Double ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: Array1(:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_2D_Double

  subroutine UnfilteredDiff_2D_Integer ( IArray1, Name1, IArray2, Name2, &
    & Width, Format, LBound, Options )
    integer, intent(in) :: IArray1(:,:)
    character(len=*), intent(in) :: Name1
    integer, intent(in) :: IArray2(:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound ! Low bound for Array
    character(len=*), intent(in), optional :: options

    real, dimension(size(iarray1,1), size(iarray1,2)) :: array1
    real, dimension(size(iarray1,1), size(iarray1,2)) :: array2
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    call diff ( Array1, Name1, Array2, Name2, &
      & Width=Width, Format=Format, &
      & LBound=LBound, Options=Options )
  end subroutine UnfilteredDiff_2D_Integer

  subroutine UnfilteredDiff_2D_Real ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Array1(:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_2D_Real

  subroutine UnfilteredDiff_3D_Double ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: Array1(:,:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_3D_Double

  subroutine UnfilteredDiff_3D_Real ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Array1(:,:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_3D_Real

  subroutine UnfilteredDiff_4D_Double ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0d0)
    real(rk), intent(in) :: Array1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_4D_Double

  subroutine UnfilteredDiff_4D_Real ( Array1, Name1, Array2, Name2, &
    & Width, Format, LBound, Options )
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Array1(:,:,:,:)
    character(len=*), intent(in) :: Name1
    real(rk), intent(in) :: Array2(:,:,:,:)
    character(len=*), intent(in) :: Name2
    integer, intent(in), optional :: Width
    character(len=*), intent(in), optional :: Format
    integer, intent(in), optional :: LBound
    character(len=*), intent(in), optional :: options

    include "unfiltereddiff.f9h"
  end subroutine UnfilteredDiff_4D_Real

  subroutine Zdonewithdiff ( Array1, Array2 )
    ! Not an actual subroutine--here only so "make depends" knows we
    ! have a dependency on donewithdiff.f9h
    integer, parameter :: RK = kind(1.0e0)
    real(rk), intent(in) :: Array1(:,:,:,:)
    real(rk), intent(in) :: Array2(:,:,:,:)

    integer, parameter                       :: MAXPCTS = 10
    character(len=64)                        :: DiffName
    logical                                  :: alreadySbtrcted = .false.
    real(rk), dimension(3, MAXPCTS)          :: array1AtNAbs
    real(rk), dimension(3, MAXPCTS)          :: array2AtNAbs
    real(rk), dimension(3, MAXPCTS)          :: array1AtNRel
    real(rk), dimension(3, MAXPCTS)          :: array2AtNRel
    real(rk), dimension(2)                   :: exvalues
    real(rk), dimension(2)                   :: exratios
    real(rk)                                 :: fillvalue
    real(rk), dimension(MAXPCTS)             :: gaps
    type(Stat_T), dimension(MAXPCTS)         :: gapStat
    real(rk), dimension(MAXPCTS)             :: gapratios
    type(Stat_T), dimension(MAXPCTS)         :: gapRatioStat
    real(rk)                                 :: minratio
    real(rk)                                 :: maxratio
    real(rk)                                 :: medianratio
    real(rk)                                 :: meanratio
    real(rk)                                 :: minvalue
    real(rk)                                 :: maxvalue
    real(rk)                                 :: medianvalue
    real(rk)                                 :: meanvalue
    real(rk)                                 :: rmsvalue
    integer                                  :: numTot
    integer                                  :: numEqual
    real(rk)                                 :: pctDiff
    real(rk)                                 :: pctEqual
    real(rk), dimension(MAXPCTS)             :: pcts
    real(rk), dimension(MAXPCTS)             :: pctratios
    real(rk), dimension(MAXPCTS)             :: pctMaxGaps
    real(rk), dimension(MAXPCTS)             :: pctMaxGapAsRatios
    real(rk), dimension(MAXPCTS)             :: pctMeanGaps
    real(rk), dimension(MAXPCTS)             :: pctMaxRatios
    real(rk), dimension(MAXPCTS)             :: pctMaxRatioAsGaps
    real(rk), dimension(MAXPCTS)             :: pctMeanRatios
    real(rk)                                 :: refmin, refmax, refrms
    real(rk)                                 :: rmsratio
    real(rk)                                 :: stddev
    real(rk)                                 :: stddevratio
    real                                     :: t1 = 0
    real                                     :: t2 = 0
    real(rk), dimension(MAXPCTS,7)           :: TheTable
    logical, parameter                       :: DEBUG = .false.
    real(rk), dimension(MAXPCTS), parameter  :: PCTAges = &
      & (/ 99.9, 99.8, 99.7, 99.5, 99., 98., 97., 95., 90., 80. /)
    real, dimension(MAXPCTS), parameter      :: PCTFactors = &
      & (/ .01, .02, .05, .1, .2, .5, 1., 2., 5., 10. /)
    ! Executable
  contains
    include "donewithdiff.f9h"
  end subroutine Zdonewithdiff

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Diff_1

! $Log$
! Revision 2.9  2020/01/30 18:23:34  pwagner
! Increased TooManyElements; can now diff l1b DACS
!
! Revision 2.8  2017/11/03 20:02:02  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.7  2017/07/31 23:06:41  pwagner
! Removed unneeded param
!
! Revision 2.6  2017/07/19 22:48:46  pwagner
! Get PrintRMSetc from Dump_0
!
! Revision 2.5  2016/10/21 23:12:50  vsnyder
! Remove unused USE name
!
! Revision 2.4  2016/10/06 20:22:14  pwagner
! parts commom to unfiltered and filtered diffs moved to donewithdiff.f9h
!
! Revision 2.3  2016/09/09 20:34:53  pwagner
! Added Au (Gold) brick option removing some hay from the stack of statistics
!
! Revision 2.2  2016/07/30 00:09:31  pwagner
! NAG complained at link about not having a Time_Now
!
! Revision 2.1  2016/07/28 01:41:48  vsnyder
! Initial Commit
!
