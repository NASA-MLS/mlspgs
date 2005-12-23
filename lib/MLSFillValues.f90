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

  use ieee_arithmetic, only: ieee_is_finite
  use MLSCommon, only: r4, r8, DEFAULTUNDEFINEDVALUE
  use MLSKinds ! Everything
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSStrings, only: Lowercase

  implicit none

  private

  public :: FilterValues
  public :: IsFillValue, ReplaceFillValues
  public :: IsFinite
  public :: IsMonotonic, Monotonize

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! This module contains some low level numerical stuff, hunting, interpolating
  ! etc.
!     c o n t e n t s
!     - - - - - - - -

!         Functions, operations, routines
! FilterValues                 Filters entries in two arrays
! IsFillValue                  Returns true if argument is FillValue
! IsFinite                     Returns true if argument is finite
! ReplaceFillValues            Replaces FillValue entries in an array

  interface BridgeMissingValues
    module procedure BridgeMissingValues_1dr4, BridgeMissingValues_1dr8, BridgeMissingValues_1dint
    module procedure BridgeMissingValues_2dr4, BridgeMissingValues_2dr8, BridgeMissingValues_2dint
    module procedure BridgeMissingValues_3dr4, BridgeMissingValues_3dr8, BridgeMissingValues_3dint
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
    module procedure IsFinite_REAL, IsFinite_DOUBLE, IsFinite_INTEGER
  end interface
  
  interface ReplaceFillValues
    module procedure ReplaceFill1d_r4, ReplaceFill1d_r8, ReplaceFill1d_int
    module procedure ReplaceFill2d_r4, ReplaceFill2d_r8, ReplaceFill2d_int
    module procedure ReplaceFill3d_r4, ReplaceFill3d_r8, ReplaceFill3d_int
  end interface

  interface IsMonotonic
    module procedure IsMonotonic_r4, IsMonotonic_r8, IsMonotonic_int
  end interface

  interface Monotonize
    module procedure Monotonize_1dr4, Monotonize_1dr8, Monotonize_1dint
    module procedure Monotonize_2dr4, Monotonize_2dr8, Monotonize_2dint
    module procedure Monotonize_3dr4, Monotonize_3dr8, Monotonize_3dint
  end interface

  real, parameter, private :: FILLVALUETOLERANCE = 0.2 ! Poss. could make it 1

  logical, parameter ::   DEEBUG = .false.

contains

! -------------------------------------------------  FilterValues  -----
  subroutine filterValues_REAL ( a, ATAB, b, BTAB, warn, fillValue, precision )
    ! Return arrays filtered of any fillValues
    ! or where corresponding precision array < 0
    ! or whose values are not finite
    ! "Filter" means offending elements set to 0.
    ! Returned arrays are allocated and assigned values as appropriate
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
  elemental logical function IsFillValue_int ( A, FILLVALUE )
    integer, intent(in) :: A
    integer, intent(in), optional :: FILLVALUE
    integer  :: MYFILLVALUE
    myFillValue = int(DEFAULTUNDEFINEDVALUE)
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_int = &
      & abs(a - myFillValue) < 1 ! FILLVALUETOLERANCE
  end function IsFillValue_int

  elemental logical function IsFillValue_r4 ( A, FILLVALUE )
    real(r4), intent(in) :: A
    real(r4), intent(in), optional :: FILLVALUE
    real(r4)  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_r4 = &
      & abs(a - myFillValue) < FILLVALUETOLERANCE
  end function IsFillValue_r4

  elemental logical function IsFillValue_r8 ( A, FILLVALUE )
    real(r8), intent(in) :: A
    real(r8), intent(in), optional :: FILLVALUE
    real(r8)  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_r8 = &
      & abs(a - myFillValue) < Real(FILLVALUETOLERANCE, r8)
  end function IsFillValue_r8

! ------------------------------------------------- IsFinite ---

  ! This family of routines checks to see if an arg is finite
  elemental logical function IsFinite_REAL ( A ) result( finite )
    real, intent(in) :: A
    finite = ieee_is_finite(a)
  end function isfinite_real
  elemental logical function IsFinite_DOUBLE ( A ) result( finite )
    double precision, intent(in) :: A
    finite = ieee_is_finite(a)
  end function isfinite_DOUBLE
  elemental logical function IsFinite_INTEGER ( A ) result( finite )
    integer, intent(in) :: A
    finite = .true.
  end function isfinite_INTEGER

! -------------------------------------  ReplaceFillValues  -----

  ! This family of routines replaces entries in an array
  ! based on whether they
  ! (1) are equal to FillValue; or
  ! (2) other criteria set by options
  ! The replacement values are supplied either by 
  ! newvalues, newFill, or according to options (e.g., you may interpolate)
  ! Note:
  ! When interpolating arrays with rank > 1, the interpolated-against
  ! index is the last one
  ! Thus we don't do true multi-dimensional interpolation
  ! If you wish to interpolate against another index, you must reshape

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

! ------------------------------------------------- IsMonotonic ---

  ! This family of routines checks to see if an array monotonically increases
  ! By default it returns TRUE if either monotonic increasing or decreasing
  ! If optional DIRECTION is supplied and '+', TRUE means increasing
  ! if '-', TRUE means drecreasing
  ! We are strict in the sense that any "stalling" produces FALSE
  ! Thus {0, 1, 2, 2, 3} is not monotonic
  !
  ! If you would like to know instead whether an array is "never falling"
  ! or "never rising" that's a different though related condition
  logical function IsMonotonic_int ( ARRAY, DIRECTION ) result(sooDesu)
    integer, dimension(:), intent(in) :: ARRAY
    character(len=*), intent(in), optional :: DIRECTION
    ! Internal variables
    character(len=1) :: MYDIRECTION
    integer :: i
    integer, dimension(size(array)-1) :: increment
    integer :: incMax
    integer :: incMin
    integer, parameter :: ZERO = 0
    ! Executable
    myDirection = '0' ! Either direction is OK, as long as monotonic
    if ( present(direction) ) myDirection = direction
    if ( size(array) < 2 ) then
      sooDesu = .true.
      return
    endif
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    incMin = minval(increment)
    incMax = maxval(increment)
    select case (myDirection)
    case ('+') ! Must be increasing monotonic
      sooDesu = ( incMin > zero )
    case ('-') ! Must be decreasing monotonic
      sooDesu = ( incMax < zero )
    case default ! either will do
      sooDesu = ( incMin > zero ) .or. ( incMax < zero )
    end select
  end function IsMonotonic_int

  logical function IsMonotonic_r4 ( ARRAY, DIRECTION ) result(sooDesu)
    real(r4), dimension(:), intent(in) :: ARRAY
    character(len=*), intent(in), optional :: DIRECTION
    ! Internal variables
    character(len=1) :: MYDIRECTION
    integer :: i
    real(r4), dimension(size(array)-1) :: increment
    real(r4) :: incMax
    real(r4) :: incMin
    real(r4), parameter :: ZERO = 0._r4
    ! Executable
    myDirection = '0' ! Either direction is OK, as long as monotonic
    if ( present(direction) ) myDirection = direction
    if ( size(array) < 2 ) then
      sooDesu = .true.
      return
    endif
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    incMin = minval(increment)
    incMax = maxval(increment)
    select case (myDirection)
    case ('+') ! Must be increasing monotonic
      sooDesu = ( incMin > zero )
    case ('-') ! Must be decreasing monotonic
      sooDesu = ( incMax < zero )
    case default ! either will do
      sooDesu = ( incMin > zero ) .or. ( incMax < zero )
    end select
  end function IsMonotonic_r4

  logical function IsMonotonic_r8 ( ARRAY, DIRECTION ) result(sooDesu)
    real(r8), dimension(:), intent(in) :: ARRAY
    character(len=*), intent(in), optional :: DIRECTION
    ! Internal variables
    character(len=1) :: MYDIRECTION
    integer :: i
    real(r8), dimension(size(array)-1) :: increment
    real(r8) :: incMax
    real(r8) :: incMin
    real(r8), parameter :: ZERO = 0._r8
    ! Executable
    myDirection = '0' ! Either direction is OK, as long as monotonic
    if ( present(direction) ) myDirection = direction
    if ( size(array) < 2 ) then
      sooDesu = .true.
      return
    endif
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    incMin = minval(increment)
    incMax = maxval(increment)
    select case (myDirection)
    case ('+') ! Must be increasing monotonic
      sooDesu = ( incMin > zero )
    case ('-') ! Must be decreasing monotonic
      sooDesu = ( incMax < zero )
    case default ! either will do
      sooDesu = ( incMin > zero ) .or. ( incMax < zero )
    end select
  end function IsMonotonic_r8

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

  ! This family of subroutines makes arrays of values monotonically (increasing)
  subroutine Monotonize_1dint(values)
    ! Args
    integer, dimension(:), intent(inout) :: values
    ! Internal variables
    integer :: dx
    integer :: x1
    integer :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dint

  subroutine Monotonize_1dr4(values)
    ! Args
    real(r4), dimension(:), intent(inout) :: values
    ! Internal variables
    real(r4) :: dx
    real(r4) :: x1
    real(r4) :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dr4

  subroutine Monotonize_1dr8(values)
    ! Args
    real(r8), dimension(:), intent(inout) :: values
    ! Internal variables
    real(r8) :: dx
    real(r8) :: x1
    real(r8) :: x2
    include 'Monotonize.f9h'
  end subroutine Monotonize_1dr8

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
      call output('*** Error in MLSCommon module ***', advance='yes')
    endif
    call output(trim(message), advance='no')
    call blanks(3)
    if ( present(int1) ) write(*,'(i4)',advance='no') int1
    call blanks(3)
    if ( present(int2) ) write(*,'(i4)', advance='no') int2
    if ( .not. keepgoing ) stop
  end subroutine announce_error

  subroutine blanks(num, advance)
    integer, intent(in) :: num
    character(len=*), optional, intent(in) :: advance
    integer :: i
    character(len=3) :: myAdvance
    myAdvance = 'no'
    if ( present(advance) ) myAdvance = advance
    do i=1, num
      write(*, '(a1)', advance='no') ' '
    enddo
    if ( myAdvance == 'yes' ) write(*, '(a1)', advance='yes') ''
    
  end subroutine blanks

  subroutine output(str, advance)
    character(len=*), intent(in) :: str
    character(len=*), optional, intent(in) :: advance
    write(*, '(a1)', advance=advance) trim(str)
    
  end subroutine output

  ! ----------------------------------  DUMP_NAME_V_PAIRS  -----
  subroutine DUMP_NAME_V_PAIRS ( VALUES, WIDTH )
    integer, intent(in)                         :: values(:)
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    integer :: MyWidth
    character(len=24) :: myName
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          write(myName, *) 'integer # ', l, ': '
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          write(*,'(i4)', advance='no') values(l)
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSFillValues
!=============================================================================

!
! $Log$
! Revision 2.2  2005/12/23 03:10:31  vsnyder
! Make some routines more generic, using include
!
! Revision 2.1  2005/12/16 00:00:23  pwagner
! Created to hold fillValue-related stuff
!
