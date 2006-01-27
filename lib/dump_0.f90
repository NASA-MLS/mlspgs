! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module DUMP_0

! Low-level dump routines -- for some arrays of intrinsic type.
! Should handle most combinations of rank and type
! Behavior depends on optional parameters
! Actual output device determined by output_m module

  use ieee_arithmetic, only: ieee_is_finite
  use MLSCommon, only: DEFAULTUNDEFINEDVALUE
  use MLSFillValues, only : FilterValues, IsFinite
  use MLSSets, only: FindAll
  use MLSStats1, only: ALLSTATS
  use MLSStringLists, only: GetStringElement, NumStringElements
  use OUTPUT_M, only: BLANKS, NEWLINE, OUTPUT

  implicit none
  private
! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!     (parameters)
! AfterSub                 character printed between row, col id and data
! STATSONONELINE           stats, rms each printed on a single line

!     (subroutines and functions)
! DIFF                     dump diffs between pair of arrays of numeric type
! DUMP                     dump an array to output
! DUMP_NAME_V_PAIRS        dump an array of paired names and values
! SELFDIFF                 dump increments between successive array values
! === (end of toc) ===

! === (start of api) ===
! diff ( array1, char* name1, array2, char* name2,
!      [fillvalue], [log clean], [int width], [char* format],
!      [log wholearray], [log stats], [log rms], [int lbound] ) 
!       where array1, array2 can be 1, 2, or 3d arrays of 
!       ints, reals, or doubles, compatible in size and type
!       and fillValue is a scalar of the same type, if present
! dump ( array, char* name,
!      [fillvalue], [log clean], [int width], [char* format],
!      [log wholearray], [log stats], [log rms], [int lbound] ) 
!       where array can be a 1, 2, or 3d array of
!       chars, ints, reals, or doubles,
!       and fillValue is a scalar of the same type, if present
! dump ( strlist string, char* name, [char* fillvalue], [log clean] )
! dump ( log countEmpty, strlist keys, strlist values, char* name, 
!       [char* separator] )
! dump_name_v_pairs ( values, strlist names,
!      [[log clean], [char* format, [int width]] ) 
!       where values can be a 1d array of ints or reals, and
!       names is a string list of corresponding names

! in the above, a string list is a string of elements (usu. comma-separated)
! === (end of api) ===

  public :: DIFF, DUMP, DUMP_NAME_V_PAIRS, SELFDIFF

  interface DIFF        ! dump diffs between pair of n-d arrays of numeric type
    module procedure DIFF_1D_DOUBLE, DIFF_1D_INTEGER, DIFF_1D_REAL
    module procedure DIFF_2D_DOUBLE, DIFF_2D_REAL
    module procedure DIFF_3D_DOUBLE, DIFF_3D_REAL
  end interface
  interface DUMP        ! dump n-d arrays of homogeneous type
    module procedure DUMP_1D_CHAR, DUMP_1D_COMPLEX, DUMP_1D_DCOMPLEX
    module procedure DUMP_1D_DOUBLE, DUMP_1D_INTEGER
    module procedure DUMP_1D_LOGICAL, DUMP_1D_REAL
    module procedure DUMP_2D_CHAR, DUMP_2D_COMPLEX, DUMP_2D_DCOMPLEX
    module procedure DUMP_2D_DOUBLE, DUMP_2D_INTEGER
    module procedure DUMP_2D_LOGICAL, DUMP_2D_REAL
    module procedure DUMP_3D_CHAR, DUMP_3D_DOUBLE, DUMP_3D_INTEGER
    module procedure DUMP_3D_REAL
    module procedure DUMP_HASH, DUMP_STRLIST
  end interface
  interface DUMP_NAME_V_PAIRS   ! dump name-value pairs, names in string list
    module procedure DUMP_NAME_V_PAIRS_DOUBLE, DUMP_NAME_V_PAIRS_INTEGER
    module procedure DUMP_NAME_V_PAIRS_REAL
  end interface
  interface SELFDIFF       ! dump increments between successive array values
    module procedure SELFDIFF_INTEGER
    module procedure SELFDIFF_REAL
    module procedure SELFDIFF_DOUBLE
  end interface
  interface printRMSetc
    module procedure printRMSetc_DOUBLE, printRMSetc_INT, printRMSetc_REAL
  end interface
  interface SAY_FILL
    module procedure SAY_FILL_CHAR, SAY_FILL_DOUBLE, SAY_FILL_INT, SAY_FILL_REAL
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  character, public, parameter :: AfterSub = '#'
  logical, public, save ::   STATSONONELINE = .true.
  logical, parameter ::   DEEBUG = .false.
  logical :: myStats, myRMS, myWholeArray
  integer :: numNonFill, numFill
  real :: pctNonFill, pctFill

  character(*), parameter :: MyFormatDefault = '(1pg14.6)'
  character(*), parameter :: MyFormatDefaultCmplx = &
    & '(1x,"(",1pg13.6,",",1pg13.6,")")'

contains

 ! ---------------------------------------------  DIFF_1D_DOUBLE  -----
  subroutine DIFF_1D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY1(:)
    character(len=*), intent(in), optional :: NAME1
    double precision, intent(in) :: ARRAY2(:)
    character(len=*), intent(in), optional :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    double precision, dimension(size(array1)) :: filtered1
    double precision, dimension(size(array2)) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_1D_DOUBLE

 ! ---------------------------------------------  DIFF_1D_INTEGER  -----
  subroutine DIFF_1D_INTEGER ( IARRAY1, NAME1, IARRAY2, NAME2, &
    & IFILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    integer, intent(in) :: IARRAY1(:)
    character(len=*), intent(in), optional :: NAME1
    integer, intent(in) :: IARRAY2(:)
    character(len=*), intent(in), optional :: NAME2
    integer, intent(in), optional :: IFILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    real, dimension(size(iarray1)) :: array1
    real, dimension(size(iarray2)) :: array2
    real :: fillValue
    ! So we don't have to write an integer-version of allstats
    array1 = iarray1
    array2 = iarray2
    if ( present(iFillValue) ) then
      fillValue = iFillValue
    else
      fillValue = int(DEFAULTUNDEFINEDVALUE)
    endif
    call DIFF ( ARRAY1, NAME1, ARRAY2, NAME2, &
      & FILLVALUE=FILLVALUE, CLEAN=CLEAN, WIDTH=WIDTH, FORMAT=FORMAT, &
      & WHOLEARRAY=WHOLEARRAY, STATS=STATS, RMS=RMS, LBOUND=LBOUND )
  end subroutine DIFF_1D_INTEGER

 ! ---------------------------------------------  DIFF_1D_REAL  -----
  subroutine DIFF_1D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY1(:)
    character(len=*), intent(in), optional :: NAME1
    real, intent(in) :: ARRAY2(:)
    character(len=*), intent(in), optional :: NAME2
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    real, dimension(size(array1)) :: filtered1
    real, dimension(size(array2)) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_1D_REAL

  ! ---------------------------------------------  DIFF_2D_DOUBLE  -----
  subroutine DIFF_2D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
      & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in), optional :: NAME1
    double precision, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in), optional :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND
    !
    double precision, dimension(size(array1,1), size(array1,2)) :: filtered1
    double precision, dimension(size(array2,1), size(array2,2)) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_2D_DOUBLE

  ! ---------------------------------------------  DIFF_2D_REAL  -----
  subroutine DIFF_2D_REAL ( ARRAY1, NAME1, ARRAY2, Name2, &
      & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY1(:,:)
    character(len=*), intent(in), optional :: NAME1
    real, intent(in) :: ARRAY2(:,:)
    character(len=*), intent(in), optional :: NAME2
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND
    !
    real, dimension(size(array1,1), size(array1,2)) :: filtered1
    real, dimension(size(array2,1), size(array2,2)) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_2D_REAL

  ! ---------------------------------------------  DIFF_3D_DOUBLE  -----
  subroutine DIFF_3D_DOUBLE ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in), optional :: NAME1
    double precision, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in), optional :: NAME2
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    double precision, dimension(size(array1,1), size(array1,2), size(array1,3)) :: filtered1
    double precision, dimension(size(array2,1), size(array2,2), size(array2,3)) :: filtered2
    double precision :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_3D_DOUBLE

  ! ---------------------------------------------  DIFF_3D_REAL  -----
  subroutine DIFF_3D_REAL ( ARRAY1, NAME1, ARRAY2, NAME2, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY1(:,:,:)
    character(len=*), intent(in), optional :: NAME1
    real, intent(in) :: ARRAY2(:,:,:)
    character(len=*), intent(in), optional :: NAME2
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    real, dimension(size(array1,1), size(array1,2), size(array1,3)) :: filtered1
    real, dimension(size(array2,1), size(array2,2), size(array2,3)) :: filtered2
    real :: refmin, refmax, refrms
    include "diff.f9h"
  end subroutine DIFF_3D_REAL

  ! -----------------------------------------------  DUMP_1D_CHAR  -----
  subroutine DUMP_1D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN, TRIM )
    character(len=*), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    logical, intent(in), optional :: TRIM

    integer :: J, K
    integer :: LON
    logical :: MyClean
    logical :: MyTRIM
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myTrim = .false.
    if ( present(trim) ) myTrim = trim
    lon = len(array(1))
    if ( myTrim ) lon = maxval(len_trim(array))

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1)(1:lon), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), 10
        if (.not. myClean) then
          if ( any(array(j:min(j+9, size(array))) /= myFillValue) ) then
            call say_fill ( (/ j-1, size(array) /), numZeroRows, myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+9, size(array))) /= myFillValue) ) then
          do k = j, min(j+9, size(array))
              call output ( array(k)(1:lon) // ' ' )
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      call say_fill ( (/ j-10, size(array) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_1D_CHAR

  ! --------------------------------------------  DUMP_1D_COMPLEX  -----
  subroutine DUMP_1D_COMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: MyClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 3
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_COMPLEX

  ! -------------------------------------------  DUMP_1D_DCOMPLEX  -----
  subroutine DUMP_1D_DCOMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT

    logical :: MyClean
    integer :: J, K, MyWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 3
    if ( present(width) ) myWidth = width
    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          call output ( j, max(myWidth-1,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+myWidth-1, size(array))
          call output ( array(k), myFormat )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_DCOMPLEX

 ! ---------------------------------------------  DUMP_1D_DOUBLE  -----
  subroutine DUMP_1D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    integer :: Base
    integer, parameter :: DefaultWidth = 5
    integer :: J, K, MyWidth
    logical :: MyClean
    double precision :: myFillValue
    character(len=64) :: MyFormat
    integer :: NumZeroRows

    ! Executable
    myFillValue = 0.d0
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dumpstats.f9h'
    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = defaultWidth
    if ( present(width) ) myWidth = width
    myFormat = myFormatDefault
    if ( present(format) ) myFormat = format
    base = 0
    if ( present(lbound) ) base = lbound - 1

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 .and. base == 0 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
!         if (.not. myClean) then
!           call output ( j, max(defaultWidth-1,ilog10(size(array))+1) )
!           call output ( afterSub )
!         end if
!         do k = j, min(j+myWidth-1, size(array))
!           call output ( array(k), myFormat )
!         end do
        if (.not. myClean) then
          if ( any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
            call say_fill ( (/ j-1+base, size(array) /), numZeroRows, myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
          do k = j, min(j+myWidth-1, size(array))
           call output ( array(k), myFormat )
          end do
          call output ( '', advance='yes' )
        endif
      end do
      call say_fill ( (/ j-myWidth, size(array) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_1D_DOUBLE

  ! --------------------------------------------  DUMP_1D_INTEGER  -----
  subroutine DUMP_1D_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, FORMAT, WIDTH, WHOLEARRAY, STATS, RMS, LBOUND )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    integer :: Base, J, K
    logical :: MyClean
    integer :: MyWidth
    integer :: NumZeroRows

    integer :: myFillValue
    ! Executable
    myFillValue = 0
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dumpstats.f9h'
    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width
    base = 0
    if ( present(lbound) ) base = lbound - 1

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 .and. base == 0 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          if ( any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
            call say_fill ( (/ j-1+base, size(array) /), numZeroRows, myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
          do k = j, min(j+myWidth-1, size(array))
            if ( present(format) ) then
              call output ( array(k), format=format )
            else
              call output ( array(k), places=6 )
            end if
          end do
          call output ( '', advance='yes' )
        end if
      end do ! j
      call say_fill ( (/ j-myWidth, size(array) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_1D_INTEGER

  ! ----------------------------------------------  DUMP_1D_LOGICAL ----
  subroutine DUMP_1D_LOGICAL ( ARRAY, NAME, CLEAN, LBOUND )
    logical, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: LBOUND ! Low bound of Array

    integer :: Base, J, K
    logical :: myClean

    base = 0
    if ( present(lbound) ) base = lbound - 1

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 .and. base == 0 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), 34
        if (.not. myClean) then
          call output ( j+base, max(4,ilog10(size(array))+1) )
          call output ( afterSub )
        end if
        do k = j, min(j+33, size(array))
          call output ( array(k) )
        end do
        call output ( '', advance='yes' )
      end do
    end if
  end subroutine DUMP_1D_LOGICAL

  ! -----------------------------------------------  DUMP_1D_REAL  -----
  subroutine DUMP_1D_REAL ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    integer :: Base
    integer, parameter :: DefaultWidth = 5
    integer :: J, K, MyWidth
    logical :: myClean
    character(len=64) :: MyFormat
    integer :: NumZeroRows

    real :: myFillValue
    ! Executable
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dumpstats.f9h'

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = defaultWidth
    if ( present(width) ) myWidth = width
    myFormat = myFormatDefault
    if ( present(format) ) myFormat = format
    base = 0
    if ( present(lbound) ) base = lbound - 1

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 .and. base == 0 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1), myFormat, advance='yes' )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, size(array), myWidth
        if (.not. myClean) then
          if ( any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
            call say_fill ( (/ j-1+base, size(array) /), numZeroRows, myFillValue, inc=1 )
          else
            numZeroRows = numZeroRows + 1
          end if
        end if
        if ( myClean .or. any(array(j:min(j+myWidth-1, size(array))) /= myFillValue) ) then
          do k = j, min(j+myWidth-1, size(array))
           call output ( array(k), myFormat )
          end do
          call output ( '', advance='yes' )
        endif
      end do
      call say_fill ( (/ j-myWidth, size(array) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_1D_REAL

  ! -----------------------------------------------  DUMP_2D_CHAR  -----
  subroutine DUMP_2D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN, TRIM )
    character(len=*), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    logical, intent(in), optional :: TRIM

    integer :: I, J, K
    integer :: LON
    logical :: MyClean
    logical :: MyTRIM
    integer :: NumZeroRows
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myTrim = .false.
    if ( present(trim) ) myTrim = trim
    lon = len(array(1,1))
    if ( myTrim ) lon = maxval(len_trim(array))

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1)(1:lon), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, fillValue=fillValue, clean=clean, trim=trim )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), 10
          if (.not. myClean) then
            if ( any(array(i,j:min(j+9, size(array,2))) /= myFillValue) ) then
              call say_fill ( (/ i-1, size(array,1), j, size(array,2) /), &
                & numZeroRows, myFillValue, inc=1 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+9, size(array,2))) /= myFillValue) ) then
            do k = j, min(j+9, size(array,2))
                call output ( array(i,k)(1:lon) // ' ' )
            end do
            call output ( '', advance='yes' )
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size(array,1), j-10, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_CHAR

  ! --------------------------------------------  DUMP_2D_COMPLEX  -----
  subroutine DUMP_2D_COMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0e0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: myWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myWidth = 3
    if ( present(width) ) myWidth = width

    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), format=myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=myFormat )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), myWidth
            if (.not. myClean) then
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            end if
            do k = j, min(j+myWidth-1, size(array,2))
              call output ( array(i,k), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), width
            call output ( i, places=max(4,ilog10(size(array,1))+1) )      
            call output ( j, places=max(4,ilog10(size(array,2))+1) )      
            call output ( afterSub )                                      
            do k = i, min(i+width-1, size(array,1))
              call output ( array(k,j), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end if
    end if
  end subroutine DUMP_2D_COMPLEX

  ! --------------------------------------------  DUMP_2D_COMPLEX  -----
  subroutine DUMP_2D_DCOMPLEX ( ARRAY, NAME, CLEAN, WIDTH, FORMAT )
    integer, parameter :: RK = kind(0.0d0)
    complex(rk), intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH ! How many per line?
    character(len=*), optional :: FORMAT

    logical :: myClean
    integer :: I, J, K
    integer :: myWidth
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myWidth = 3
    if ( present(width) ) myWidth = width

    myFormat = MyFormatDefaultCmplx
    if ( present(format) ) myFormat = format

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), format=myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=myFormat )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), myWidth
            if (.not. myClean) then
              call output ( i, places=max(4,ilog10(size(array,1))+1) )
              call output ( j, places=max(4,ilog10(size(array,2))+1) )
              call output ( afterSub )
            end if
            do k = j, min(j+myWidth-1, size(array,2))
              call output ( array(i,k), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), width
            call output ( i, places=max(4,ilog10(size(array,1))+1) )      
            call output ( j, places=max(4,ilog10(size(array,2))+1) )      
            call output ( afterSub )                                      
            do k = i, min(i+width-1, size(array,1))
              call output ( array(k,j), format=myFormat )
            end do
            call output ( '', advance='yes' )
          end do
        end do
      end if
    end if
  end subroutine DUMP_2D_DCOMPLEX

  ! ---------------------------------------------  DUMP_2D_DOUBLE  -----
  subroutine DUMP_2D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    logical :: myClean
    integer :: I, J, K
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myFillValue = 0.0d0
    if ( present(FillValue) ) myFillValue = FillValue
    include 'dumpstats.f9h'

    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), 5
            if (.not. myClean) then
              if ( any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                  & numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
              do k = j, min(j+4, size(array,2))
                call output ( array(i,k), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
        call say_fill ( (/ i-1, size(array,1), j-5, size(array,2) /), &
          & numZeroRows, myFillValue )
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), 5
            if ( any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then  
              call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                & numZeroRows, myFillValue, inc=3 )
            else                                                            
              numZeroRows = numZeroRows + 1                                 
            end if                                                          
            if ( myClean .or. any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then
              do k = i, min(i+4, size(array,1))
                call output ( array(k,j), myFormat )
              end do
              call output ( '', advance='yes' )
            end if                                                          
          end do
        end do
      end if
      call say_fill ( (/ i-5, size(array,1), j-1, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_DOUBLE

  ! --------------------------------------------  DUMP_2D_INTEGER  -----
  subroutine DUMP_2D_INTEGER ( ARRAY, NAME, CLEAN, FORMAT, WIDTH )
    integer, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?

    integer :: I, J, K
    logical :: MyClean
    integer :: MyWidth
    integer :: NumZeroRows

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean, format=format, width=width )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), myWidth
          if (.not. myClean) then
            if ( any(array(i,j:min(j+myWidth-1, size(array,2))) /= 0) ) then
              call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                & numZeroRows, 0, inc=3 )
            else
              numZeroRows = numZeroRows + 1
            end if
          end if
          if ( myClean .or. any(array(i,j:min(j+myWidth-1, size(array,2))) /= 0) ) then
            do k = j, min(j+myWidth-1, size(array,2))
              if ( present(format) ) then
                call output ( array(i,k), format=format )
              else
                call output ( array(i,k), places=6 )
              end if
            end do
            call output ( '', advance='yes' )
          end if
        end do ! j
      end do ! i
      call say_fill ( (/ i-1, size(array,1), j-myWidth, size(array,2) /), &
        & numZeroRows, 0 )
    end if
  end subroutine DUMP_2D_INTEGER

  ! --------------------------------------------  DUMP_2D_LOGICAL  -----
  subroutine DUMP_2D_LOGICAL ( ARRAY, NAME, CLEAN )
    logical, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    logical, intent(in), optional :: CLEAN

    integer :: I, J, K
    logical :: MyClean
    integer, parameter :: MyWidth = 34

    myClean = .false.
    if ( present(clean) ) myClean = clean

    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2), myWidth
          if (.not. myClean) then
            call output ( i, places=max(4,ilog10(size(array,1))+1) )
            call output ( j, places=max(4,ilog10(size(array,2))+1) )
            call output ( afterSub )
          end if
          do k = j, min(j+myWidth-1, size(array,2))
            call output ( array(i,k) )
          end do
          call output ( '', advance='yes' )
        end do ! j
      end do ! i
    end if
  end subroutine DUMP_2D_LOGICAL

  ! -----------------------------------------------  DUMP_2D_REAL  -----
  subroutine DUMP_2D_REAL ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY(:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    logical :: myClean
    integer :: I, J, K
    integer :: NumZeroRows
    real :: myFillValue
    character(len=64) :: MyFormat

    myClean = .false.
    if ( present(clean) ) myClean = clean

    myFillValue = 0.0e0
    if ( present(FillValue) ) myFillValue = FillValue
    include 'dumpstats.f9h'

    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 ) then
      call dump ( array(:,1), name, clean=clean )
    else 
      call name_and_size ( name, myClean, size(array) )
      if ( size(array,2) >= min(5,size(array,1)) .or. myClean ) then
        if ( present(name) ) call output ( '', advance='yes' )
        do i = 1, size(array,1)
          do j = 1, size(array,2), 5
            if (.not. myClean) then
              if ( any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2) /), &
                  & numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j:min(j+4, size(array,2))) /= myFillValue) ) then
              do k = j, min(j+4, size(array,2))
                call output ( array(i,k), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
        call say_fill ( (/ i-1, size(array,1), j-5, size(array,2) /), &
          & numZeroRows, myFillValue )
      else ! Dump the transpose
        if ( present(name) ) call output ( ' ' )
        call output ( '(transposed)', advance='yes' )
        do j = 1, size(array,2)
          do i = 1, size(array,1), 5
            if ( any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then  
              call say_fill ( (/ i-1, size(array,1), j, size(array,2) /), &
                & numZeroRows, myFillValue, inc=1 )
            else                                                            
              numZeroRows = numZeroRows + 1                                 
            end if                                                          
            if ( myClean .or. any(array(i:min(i+4, size(array,1)),j) /= myFillValue) ) then
              do k = i, min(i+4, size(array,1))
                call output ( array(k,j), myFormat )
              end do
              call output ( '', advance='yes' )
            end if                                                          
          end do
        end do
      end if
      call say_fill ( (/ i-5, size(array,1), j-1, size(array,2) /), &
        & numZeroRows, myFillValue )
    end if
  end subroutine DUMP_2D_REAL

  ! -----------------------------------------------  DUMP_3D_CHAR  -----
  subroutine DUMP_3D_CHAR ( ARRAY, NAME, FILLVALUE, CLEAN, TRIM )
    character(len=*), intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    logical, intent(in), optional :: TRIM

    integer :: LON
    logical :: MyClean
    logical :: MyTRIM
    integer :: I, J, K, L
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many
    character(len=len(array)) :: myFillValue

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myTrim = .false.
    if ( present(trim) ) myTrim = trim
    lon = len(array(1,1,1))
    if ( myTrim ) lon = maxval(len_trim(array))
    call FindAll( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1)(1:lon), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, fillValue=fillValue, &
        & clean=clean, trim=trim )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, fillValue=fillValue, clean=clean, trim=trim )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 10
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+9, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+9, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+9, size(array,3))
                  call output ( array(i,j,l)(1:lon) // ' ' )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-10, size(array,3) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_3D_CHAR

  ! ---------------------------------------------  DUMP_3D_DOUBLE  -----
  subroutine DUMP_3D_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    double precision, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    double precision :: myFillValue
    character(len=64) :: myFormat

    ! Executable
    myFillValue = 0.d0
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dumpstats.f9h'

    myClean = .false.
    if ( present(clean) ) myClean = clean
    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name, clean=clean )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+4, size(array,3))
                call output ( array(i,j,l), myFormat )
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-5, size(array,3) /), numZeroRows, myFillValue )
   end if
  end subroutine DUMP_3D_DOUBLE

  ! --------------------------------------------  DUMP_3D_INTEGER  -----
  subroutine DUMP_3D_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, FORMAT, WIDTH, WHOLEARRAY, STATS, RMS, LBOUND )
    integer, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many numbers per line (10)?
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array

    integer :: I, J, K, L
    logical :: myClean
    integer :: MyWidth
    integer :: NumZeroRows
    integer, dimension(3) :: which, re_mainder
    integer :: how_many

    integer :: myFillValue
    ! Executable
    myFillValue = 0
    if ( present(FillValue) ) myFillValue=FillValue
    include 'dumpstats.f9h'
    myClean = .false.
    if ( present(clean) ) myClean = clean
    myWidth = 10
    if ( present(width) ) myWidth = width
    call FindAll( (/ size(array, 1), size(array, 2), size(array, 3)/), &
      & 1, which, how_many, re_mainder=re_mainder)

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), advance='yes' )
    else if ( how_many == 2 ) then
      call dump ( reshape(array, (/ re_mainder(1) /)), name, clean=clean, &
      & format=format )
    else if ( how_many == 1 ) then
      call dump ( reshape(array, (/ re_mainder(1), re_mainder(2) /)), &
        & name, clean=clean, format=format )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), myWidth
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+myWidth-1, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+myWidth-1, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+myWidth-1, size(array,3))
                if ( present(format) ) then
                  call output ( array(i,j,l), format=format )
                else
                  call output ( array(i,j,l), places=6 )
                end if
              end do
              call output ( '', advance='yes' )
            end if
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-myWidth, size(array,3) /), numZeroRows, myFillValue )
    end if
  end subroutine DUMP_3D_INTEGER

  ! ---------------------------------------------  DUMP_3D_REAL  -----
  subroutine DUMP_3D_REAL ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    real, intent(in) :: ARRAY(:,:,:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, optional, intent(in) :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND

    logical :: myClean
    integer :: I, J, K, L
    integer :: NumZeroRows
    real    :: myFillValue
    character(len=64) :: MyFormat

    ! Executable
    myFillValue = 0.
    if ( present(FillValue) ) myFillValue=FillValue
    myClean = .false.
    if ( present(clean) ) myClean = clean
    include 'dumpstats.f9h'

    myFormat = MyFormatDefault
    if ( present(format) ) myFormat = format

    numZeroRows = 0
    if ( size(array) == 0 ) then
      call empty ( name )
    else if ( size(array) == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( array(1,1,1), myFormat, advance='yes' )
    else if ( size(array,2) == 1 .and. size(array,3) == 1 ) then
      call dump ( array(:,1,1), name, clean=clean )
    else if ( size(array,3) == 1 ) then
      call dump ( array(:,:,1), name, fillValue=fillValue, clean=clean )
    else
      call name_and_size ( name, myClean, size(array) )
      if ( present(name) ) call output ( '', advance='yes' )
      do i = 1, size(array,1)
        do j = 1, size(array,2)
          do k = 1, size(array,3), 5
            if (.not. myClean) then
              if ( any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
                call say_fill ( (/ i, size(array,1), j-1, size(array,2), &
                  & k, size(array,3) /), numZeroRows, myFillValue, inc=3 )
              else
                numZeroRows = numZeroRows + 1
              end if
            end if
            if ( myClean .or. any(array(i,j,k:min(k+4, size(array,3))) /= myFillValue) ) then
              do l = k, min(k+4, size(array,3))
                call output ( array(i,j,l), myFormat )
              end do
              call output ( '', advance='yes' )
            endif
          end do
        end do
      end do
      call say_fill ( (/ i-1, size(array,1), j-1, size(array,2), &
        & k-5, size(array,3) /), numZeroRows, myFillValue )
   end if
  end subroutine DUMP_3D_REAL

  ! -----------------------------------------------  DUMP_HASH  -----
  subroutine DUMP_HASH ( COUNTEMPTY, KEYS, VALUES, NAME, SEPARATOR )
    ! Dumps a hash composed of two string lists: keys and values
    logical, intent(in)          :: COUNTEMPTY
    character(len=*), intent(in) :: KEYS
    character(len=*), intent(in) :: VALUES
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: SEPARATOR

    character( len=max(len(values), len(keys)) ) :: element
    integer :: J
    integer :: NumElements
    character(len=1) :: mySeparator
    character(len=*), parameter :: COMMA = ','

    mySeparator = COMMA
    if ( present(SEPARATOR) ) mySeparator = SEPARATOR

    NumElements = NumStringElements(keys, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else
      call name_and_size ( name, .false., NumElements )
      if ( present(name) ) call output ( '', advance='yes' )
      call output ( '(key)    =   (value)', advance='yes' )
      do j = 1, NumElements
        call GetStringElement( keys, element, j, countEmpty, separator)
        call output ( trim(element), advance='no' )
        call output( ' = ', advance='no' )
        call GetStringElement( values, element, j, countEmpty, separator)
        call output ( trim(element), advance='yes' )
      end do ! j
    end if
  end subroutine DUMP_HASH

  ! -----------------------------------------------  DUMP_STRLIST  -----
  subroutine DUMP_STRLIST ( STRING, NAME, FILLVALUE, CLEAN )
    ! Dumps a ','-separated list of strings
    character(len=*), intent(in) :: STRING
    character(len=*), intent(in), optional :: NAME
    character(len=*), intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN

    integer :: J
    logical :: MyClean
    integer :: NumElements
    character(len=len(string)) :: myFillValue
    character(len=*), parameter :: SEPARATOR = ','
    logical, parameter :: COUNTEMPTY = .true.

    myFillValue = ' '
    if ( present(FillValue) ) myFillValue = FillValue

    myClean = .false.
    if ( present(clean) ) myClean = clean

    NumElements = NumStringElements(string, countEmpty, &
     & separator)
    if ( NumElements == 0 ) then
      call empty ( name )
    else if ( NumElements == 1 ) then
      call name_and_size ( name, myClean, 1 )
      call output ( trim(string), advance='yes' )
    else
      call name_and_size ( name, myClean, NumElements )
      if ( present(name) ) call output ( '', advance='yes' )
      do j = 1, NumElements
        call GetStringElement(string, myFillValue, j, countEmpty, separator)
        call output ( trim(myFillValue), advance='yes' )
      end do ! j
    end if
  end subroutine DUMP_STRLIST

  ! -----------------------------------  DUMP_NAME_V_PAIRS  -----
  ! Another hash-like dump:
  ! Show names and related (numerical) values
  ! -----------------------------------  DUMP_NAME_V_PAIRS_DOUBLE  -----
  subroutine DUMP_NAME_V_PAIRS_DOUBLE ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    double precision, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'double # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_DOUBLE

  ! ----------------------------------  DUMP_NAME_V_PAIRS_INTEGER  -----
  subroutine DUMP_NAME_V_PAIRS_INTEGER ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    integer, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'integer # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format=format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_INTEGER

  ! -------------------------------------  DUMP_NAME_V_PAIRS_REAL  -----
  subroutine DUMP_NAME_V_PAIRS_REAL ( VALUES, NAMES, CLEAN, FORMAT, WIDTH )
    real, intent(in)                         :: values(:)
    character(len=*), intent(in), optional :: NAMES
    logical, intent(in), optional :: CLEAN
    character(len=*), intent(in), optional :: FORMAT
    integer, intent(in), optional :: WIDTH ! How many pairs per line (1)?
    
    integer :: J, K, L
    logical :: MyClean
    integer :: MyWidth
    character(len=24) :: myName
    myClean = .false.
    if ( present(clean) ) myClean = clean
    MyWidth = 1
    if ( present(width) ) myWidth = max(width, 1)
    if ( size(values) < 1 ) return
    if ( myClean ) call output(' ', advance='yes')
    l = 0
    do j=1, size(values), MyWidth
      do k=1, MyWidth
        call blanks(3, advance='no')
        l = l + 1
        if ( l <= size(values) ) then
          if ( present (names) ) then
            call GetStringElement(names, myName, l, .true.)
          else
            write(myName, *) 'single # ', l, ': '
          end if
          call output(myName,  advance='no')
          call blanks(3, advance='no')
          call output(values(l), format, advance='no')
        end if
      end do
      call output(' ', advance='yes')
    end do

  end subroutine DUMP_NAME_V_PAIRS_REAL

 ! ---------------------------------------------  SELFDIFF_DOUBLE  -----
  subroutine SELFDIFF_DOUBLE ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    ! dump the running increment == ( array(i) - array(i-1) )
    double precision, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    double precision, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    ! Local variables
    integer                                  :: i
    double precision, dimension(size(array)) :: increment
    ! Executable
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    call dump ( increment, NAME, &
      & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
  end subroutine SELFDIFF_DOUBLE

 ! ---------------------------------------------  SELFDIFF_INTEGER  -----
  subroutine SELFDIFF_INTEGER ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, FORMAT, WIDTH, WHOLEARRAY, STATS, RMS, LBOUND )
    ! dump the running increment == ( array(i) - array(i-1) )
    integer, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    integer, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    ! Local variables
    integer                                  :: i
    integer, dimension(size(array)) :: increment
    ! Executable
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    call dump ( increment, NAME, &
    & FILLVALUE, CLEAN, FORMAT, WIDTH, WHOLEARRAY, STATS, RMS, LBOUND )
  end subroutine SELFDIFF_INTEGER

 ! ---------------------------------------------  SELFDIFF_REAL  -----
  subroutine SELFDIFF_REAL ( ARRAY, NAME, &
    & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
    ! dump the running increment == ( array(i) - array(i-1) )
    real, intent(in) :: ARRAY(:)
    character(len=*), intent(in), optional :: NAME
    real, intent(in), optional :: FILLVALUE
    logical, intent(in), optional :: CLEAN
    integer, intent(in), optional :: WIDTH
    character(len=*), intent(in), optional :: FORMAT
    logical, intent(in), optional :: WHOLEARRAY
    logical, intent(in), optional :: STATS
    logical, intent(in), optional :: RMS
    integer, intent(in), optional :: LBOUND ! Low bound for Array
    ! Local variables
    integer                                  :: i
    real, dimension(size(array)) :: increment
    ! Executable
    if ( size(array) < 2 ) return
    do i=2, size(array)
      increment(i-1) = array(i) - array(i-1)
    enddo
    call dump ( increment, NAME, &
      & FILLVALUE, CLEAN, WIDTH, FORMAT, WHOLEARRAY, STATS, RMS, LBOUND )
  end subroutine SELFDIFF_REAL

  ! ------------------------------------------------------  Empty  -----
  subroutine Empty ( Name )
    character(len=*), intent(in), optional :: Name

    if ( present(name) ) then
      call output ( name )
      call output ( ' is ' )
    end if
    call output ( 'empty', advance='yes' )

  end subroutine Empty

  ! -----------------------------------------------------  ILOG10  -----
  integer function ILOG10(int)
    integer, intent(in) :: int
    ilog10=nint(log10(real(int)))
  end function ILOG10

  ! ----------------------------------------------  Name_And_Size  -----
  subroutine Name_And_Size ( Name, Clean, Size )
    character(len=*), intent(in), optional :: Name
    logical, intent(in) :: Clean
    integer, intent(in) :: Size

    if ( present(name) ) then
      call output ( name )
      if ( clean ) then 
        call output ( trim(" \ ") ) ! This goofiness is to outwit an incorrect
                                    ! Intel compiler.
        call output ( size )
      end if
      if ( size == 1 ) call output ( ' ' )
    end if

  end subroutine Name_And_Size

  ! ----------------------------------------------  printRMSetc  -----
  ! This family of routines prints a nicely-formatted list of min, max, etc.
  ! using output
  subroutine printRMSetc_double ( Name, min, max, rms, mean  )
    character(len=*), intent(in) :: Name
    double precision, intent(in) :: min
    double precision, intent(in) :: max
    double precision, intent(in) :: rms
    double precision, intent(in), optional :: mean
    if ( STATSONONELINE ) then
      call output ( trim(name), advance='no' )
      call output ( ' min : max, rms: ', advance='no' )
      if ( present(mean) ) call output ( ' mean: ', advance='no' )
      call output(min, advance='no')
      call output(' : ', advance='no')
      call output(max, advance='no')
      call output(', ', advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) then
        call output ( ': ', advance='no' )
        call output ( mean, advance='no' )
      end if
      call newline
    else
      call output ( trim(name), advance='no' )
      call output ( ' min      max        rms', advance='no' )
      if ( present(mean) ) call output ( '       mean ', advance='no' )
      call newline
      call output(min, advance='no')
      call output(max, advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) call output ( mean, advance='no' )
      call newline
    end if
  end subroutine printRMSetc_double

  subroutine printRMSetc_real ( Name, min, max, rms, mean  )
    character(len=*), intent(in) :: Name
    real, intent(in) :: min
    real, intent(in) :: max
    real, intent(in) :: rms
    real, intent(in), optional :: mean
    if ( STATSONONELINE ) then
      call output ( trim(name), advance='no' )
      call output ( ' min : max, rms: ', advance='no' )
      if ( present(mean) ) call output ( ' mean: ', advance='no' )
      call output(min, advance='no')
      call output(' : ', advance='no')
      call output(max, advance='no')
      call output(', ', advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) then
        call output ( ': ', advance='no' )
        call output ( mean, advance='no' )
      end if
      call newline
    else
      call output ( trim(name), advance='no' )
      call output ( ' min      max        rms', advance='no' )
      if ( present(mean) ) call output ( '       mean ', advance='no' )
      call newline
      call output(min, advance='no')
      call output(max, advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) call output ( mean, advance='no' )
      call newline
    end if
  end subroutine printRMSetc_real

  subroutine printRMSetc_int ( Name, min, max, rms, mean  )
    character(len=*), intent(in) :: Name
    integer, intent(in) :: min
    integer, intent(in) :: max
    real, intent(in) :: rms
    real, intent(in), optional :: mean
    if ( STATSONONELINE ) then
      call output ( trim(name), advance='no' )
      call output ( ' min : max, rms: ', advance='no' )
      if ( present(mean) ) call output ( ' mean: ', advance='no' )
      call output(min, advance='no')
      call output(' : ', advance='no')
      call output(max, advance='no')
      call output(', ', advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) then
        call output ( ': ', advance='no' )
        call output ( mean, advance='no' )
      end if
      call newline
    else
      call output ( trim(name), advance='no' )
      call output ( ' min      max        rms', advance='no' )
      if ( present(mean) ) call output ( '       mean ', advance='no' )
      call newline
      call output(min, advance='no')
      call output(max, advance='no')
      call output(rms, advance='no')
      if ( present(mean) ) call output ( mean, advance='no' )
      call newline
    end if
  end subroutine printRMSetc_int

  ! ----------------------------------------------  Say_Fill_Char  -----
  subroutine Say_Fill_Char ( Subs, NumZeroRows, Fill, Inc  )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    character(len=*), intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows  )
      call output ( '"', advance='no' )
      call output ( trim(fill), advance='no' )
      call output ( '" not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Char

  ! --------------------------------------------  Say_Fill_Double  -----
  subroutine Say_Fill_Double ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    double precision, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Double

  ! -----------------------------------------------  Say_Fill_Int  -----
  subroutine Say_Fill_Int ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    integer, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Int

  ! ----------------------------------------------  Say_Fill_Real  -----
  subroutine Say_Fill_Real ( Subs, NumZeroRows, Fill, Inc )
    integer, intent(in) :: Subs(:)
    integer, intent(inout) :: NumZeroRows
    real, intent(in) :: Fill
    integer, intent(in), optional :: Inc
    if ( numZeroRows /= 0 ) then
      call say_subs ( subs, numZeroRows )
      call output ( fill, advance='no' )
      call output ( ' not printed.', advance='yes' )
      numZeroRows = 0
    end if
    if ( present(inc) ) &
      & call say_subs_only ( (/ subs(:inc-1), subs(inc)+1, subs(inc+1:) /) )
  end subroutine Say_Fill_Real

  ! ----------------------------------------------------  Say_Subs -----
  subroutine Say_Subs ( Subs, NumZeroRows )
    integer, intent(in) :: Subs(:)
    integer, intent(in) :: NumZeroRows
    call say_subs_only ( subs )
    call output ( ' ' )
    call output ( numZeroRows )
    call output ( ' rows of ', advance='no' )
  end subroutine Say_Subs

  ! -----------------------------------------------  Say_Subs_Only -----
  subroutine Say_Subs_Only ( Subs )
    integer, intent(in) :: Subs(:)
    integer :: I
    do i = 1, size(subs), 2
      call output ( subs(i), places=max(4,ilog10(subs(i+1))+1) )
    end do
    call output ( afterSub )
  end subroutine Say_Subs_Only

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DUMP_0

! $Log$
! Revision 2.52  2006/01/27 01:01:37  pwagner
! Can now dump hashes
!
! Revision 2.51  2005/12/17 00:58:56  pwagner
! dumps of rms, etc. should appear uniform
!
! Revision 2.50  2005/12/16 23:25:13  pwagner
! dumpSize moved from dump0 to output_m
!
! Revision 2.49  2005/12/16 00:04:06  pwagner
! Changes to reflect new MLSFillValues module
!
! Revision 2.48  2005/11/04 18:49:02  pwagner
! Added SelfDiff procedures to dump diffs among consecutive 1d array elems
!
! Revision 2.47  2005/10/03 18:05:52  pwagner
! Allocated memory now dumped in units of MEMORY_UNITS
!
! Revision 2.46  2005/07/20 01:33:47  vsnyder
! Simplify DumpSize routines
!
! Revision 2.45  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.44  2005/05/12 20:43:54  pwagner
! Added diff routines
!
! Revision 2.43  2004/12/14 21:31:59  pwagner
! Optional trim arg added to multi-dim char dumps
!
! Revision 2.42  2004/08/16 17:09:14  pwagner
! 3d integer dumps interface redone like others
!
! Revision 2.41  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.40  2004/07/23 19:47:20  vsnyder
! Add LBOUND to dump_1d_[double,integer,real]
!
! Revision 2.39  2004/07/23 18:34:59  vsnyder
! Add LBOUND argument to Dump_1d_logical
!
! Revision 2.38  2004/07/21 19:57:21  pwagner
! New rms, wholeArray options to some dumps
!
! Revision 2.37  2004/06/11 19:04:29  pwagner
! which_ints_are_it renamed findAll, moved MLSSets.f90
!
! Revision 2.36  2004/05/27 23:25:33  pwagner
! Added stats parameter; also more 1-d get fillvalue
!
! Revision 2.35  2004/04/05 17:47:43  livesey
! Split dumpsize into real and integer versions
!
! Revision 2.34  2004/04/03 05:43:23  livesey
! Added DumpSize
!
! Revision 2.33  2004/03/30 00:44:10  vsnyder
! Remove unused variable declaration
!
! Revision 2.32  2004/02/26 21:53:31  pwagner
! Can dump ,-separated string list
!
! Revision 2.31  2004/01/21 22:02:19  vsnyder
! Don't use number of entities per line as default width for counter
!
! Revision 2.30  2003/09/19 02:00:14  vsnyder
! More about the goofy Intel compiler
!
! Revision 2.29  2003/09/15 17:43:41  livesey
! Cosmetic change for fussy (and wrong) intel compiler
!
! Revision 2.28  2003/09/06 00:48:40  vsnyder
! Specify default formats with a module parameter instead of literals.
! Change default (1x,1pg13.6) to (1pg14.6) to avoid problems with length
! calculation in output_m.
!
! Revision 2.27  2003/08/08 20:45:42  vsnyder
! Made say_fill_* generic, made them test for numZeroRows, and made them
! optionally do say_subs_only.  This simplified several dump routines.
! Added optional FORMAT arguments in several more routines.
!
! Revision 2.26  2003/07/04 02:41:33  vsnyder
! Substantial simplification by putting little things into subroutines
!
! Revision 2.25  2003/07/02 01:07:27  vsnyder
! Add complex output
!
! Revision 2.24  2003/05/21 19:20:40  vsnyder
! Start a new line after \ 1 if size==1 and clean
!
! Revision 2.23  2003/05/06 00:15:03  pwagner
! Fixed incompatibility with FilterShapes
!
! Revision 2.20.2.4  2003/05/05 23:00:05  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.20.2.3  2003/04/18 20:26:05  vsnyder
! Add WIDTH and FORMAT arguments to 1D_REAL and 1D_DOUBLE
!
! Revision 2.20.2.2  2003/03/27 23:18:33  vsnyder
! Put new-lines in better places
!
! Revision 2.20.2.1  2003/03/14 00:25:47  vsnyder
! Add Dump_2D_Logical, cosmetic changes
!
! Revision 2.20  2002/12/02 23:34:14  pwagner
! Now can dump name/value pairs
!
! Revision 2.19  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.18  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.17  2002/02/14 23:21:18  vsnyder
! Work on dumping masks
!
! Revision 2.16  2001/12/08 00:47:51  pwagner
! Added dump_1d_real for s.p. arrays
!
! Revision 2.15  2001/11/29 23:50:53  pwagner
! Added optional blase arg to dump_nd_char; fixed bug where optional
! format not passed from dump_3d_int
!
! Revision 2.14  2001/11/28 23:32:01  livesey
! Fixed bug where dump_2d_integer didn't pass format to 1d dump.
!
! Revision 2.13  2001/10/25 23:30:39  pwagner
! Improved dump_nd_double to skip rows (e.g., of zeros)
!
! Revision 2.12  2001/10/24 18:11:14  pwagner
! which_ints_are_it now works properly
!
! Revision 2.11  2001/10/23 22:40:37  pwagner
! Now dumps 1d,2d,3d char arrays and 3d ints
!
! Revision 2.10  2001/09/28 22:43:20  vsnyder
! Don't print rows of zeroes
!
! Revision 2.9  2001/09/11 22:52:32  livesey
! Added printing of sizes
!
! Revision 2.8  2001/05/11 22:44:54  vsnyder
! Print transpose of 2d-double if it would take fewer lines.  Get rid of
! double printing of "without mask"
!
! Revision 2.7  2001/05/08 20:27:24  vsnyder
! Added an optional 'format' argument in a few more places
!
! Revision 2.6  2001/05/08 17:21:02  livesey
! Added a `clean' option to the array dumps.  This omits the indices at
! the start, making it easier for other programs to read output.
!
! Revision 2.5  2001/05/03 02:12:34  vsnyder
! Insert copyright notice, clean up CVS stuff, cosmetics
!
! Revision 2.4  2001/03/10 03:39:58  vsnyder
! Improve handling of "name" if size==1 or size==0
!
! Revision 2.3  2001/03/02 01:32:08  livesey
! Handles larger arrays better
!
! Revision 2.2  2001/02/28 21:35:27  livesey
! Added dump logical 1d
!
! Revision 2.1  2000/09/13 20:38:50  vsnyder
! Initial code
!
