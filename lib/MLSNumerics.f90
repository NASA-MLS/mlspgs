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
module MLSNumerics              ! Some low level numerical stuff
!=============================================================================

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use DUMP_0, only : DUMP
  use MatrixModule_0, only: CreateBlock, M_Absent, MatrixElement_T, Sparsify
  use MLSCommon, only : DEFAULTUNDEFINEDVALUE, R4, R8, Rm, &
    & filterValues
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSStrings, only: Capitalize, lowercase

  implicit none

  private
  public :: Dump, EssentiallyEqual, Hunt, InterpolateArraySetup
  public :: InterpolateArrayTeardown, InterpolateValues
  public :: IsFillValue, ReplaceFillValues

  type, public :: Coefficients_R4
    private
    integer, pointer :: LowerInds(:) => NULL(), UpperInds(:) => NULL()
    !{ Coefficients for linear interpolation:
    !  {\tt Gap} $= g = x_{i+1}-x_i$.  $A = \frac{x_{i+1}-x}g$.
    !                                  $B = 1-A = \frac{x-x_i}g$.
    !  Coefficients for differentiation in the linear case (and of the linear
    !  terms in the spline case) are just $-1$ and $+1$, so they're not
    !  computed here.
    real(r4), pointer :: A(:) => NULL(), B(:) => NULL(), Gap(:) => NULL()
    !{ Coefficients for spline interpolation:
    !  $C = (A^3-A) \frac{g^2}6$.  $D = (B^3-B) \frac{g^2}6$.
    real(r4), pointer :: C(:) => NULL(), D(:) => NULL()
    real(r4), pointer :: Sig(:) => NULL() ! for second derivative guesser
    !{ Coefficients for spline derivatives:
    !  $E = \frac{\text{d}C}{\text{d}A} = \frac6g \frac{\text{d}C}{\text{d}x}
    !     = 3 A^2 - 1$.
    !  $F = \frac{\text{d}D}{\text{d}A} = \frac6g \frac{\text{d}D}{\text{d}x}
    !     = 3 B^2 - 1$.
    real(r4), pointer :: E(:) => NULL(), F(:) => NULL()
    !{ Coefficients for integration:
    !  $\int A \text{d}x =  \frac{x(x_{i+1}-\frac{x}2)}g$.
    !  $\int B \text{d}x = -\frac{x(x_i    -\frac{x}2)}g
    !                    = x - \int A \text{d}x$.\\
    !  $\int C \text{d}x = -\frac{g^2}6
    !                      \left( \int A \text{d}x + \frac{g A^4}4 \right)$.
    !  $\int D \text{d}x = -\frac{g^2}6
    !                      \left( \int B \text{d}x - \frac{g B^4}4 \right)$.
    real(r4), pointer :: AI(:) => NULL(), BI(:) => NULL(), &
      &                  CI(:) => NULL(), DI(:) => NULL()
    ! Stuff for extrapolation == "B"ad
    logical, pointer :: BadValue(:) => NULL()
  end type Coefficients_R4

  type, public :: Coefficients_R8
    private
    integer, pointer :: LowerInds(:) => NULL(), UpperInds(:) => NULL()
    !{ Coefficients for linear interpolation:
    !  {\tt Gap} $= g = x_{i+1}-x_i$.  $A = \frac{x_{i+1}-x}g$.
    !                                  $B = 1-A = \frac{x-x_i}g$.
    !  Coefficients for differentiation in the linear case (and of the linear
    !  terms in the spline case) are just $-1$ and $+1$, so they're not
    !  computed here.
    real(r8), pointer :: A(:) => NULL(), B(:) => NULL(), Gap(:) => NULL()
    !{ Coefficients for spline interpolation:
    !  $C = (A^3-A) \frac{g^2}6$.  $D = (B^3-B) \frac{g^2}6$.
    real(r8), pointer :: C(:) => NULL(), D(:) => NULL()
    real(r8), pointer :: Sig(:) => NULL() ! for second derivative guesser
    !{ Coefficients for spline derivatives:
    !  $E = \frac{\text{d}C}{\text{d}A} = \frac6g \frac{\text{d}C}{\text{d}x}
    !     = 3 A^2 - 1$.
    !  $F = \frac{\text{d}D}{\text{d}A} = \frac6g \frac{\text{d}D}{\text{d}x}
    !     = 3 B^2 - 1$.
    real(r8), pointer :: E(:) => NULL(), F(:) => NULL()
    !{ Coefficients for integration:
    !  $\int A \text{d}x =  \frac{x(x_{i+1}-\frac{x}2)}g$.
    !  $\int B \text{d}x = -\frac{x(x_i    -\frac{x}2)}g
    !                    = x - \int A \text{d}x$.\\
    !  $\int C \text{d}x = -\frac{g^2}6
    !                      \left( \int A \text{d}x + \frac{g A^4}4 \right)$.
    !  $\int D \text{d}x = -\frac{g^2}6
    !                      \left( \int B \text{d}x - \frac{g B^4}4 \right)$.
    real(r8), pointer :: AI(:) => NULL(), BI(:) => NULL(), &
      &                  CI(:) => NULL(), DI(:) => NULL()
    ! Stuff for extrapolation == "B"ad
    logical, pointer :: BadValue(:) => NULL()
  end type Coefficients_R8

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
! Dump                         Dump coefficients structure
! EssentiallyEqual             Returns true if two real arguments 'close enough'
!                                (See comments below for interpretation
!                                 of array versions)
! Hunt                         Finds index of item(s) in list closest to prey
! HuntArray                    Hunts for multiple items
! HuntScalar                   Hunts for just one
! InterpolateArraySetup        Compute coefficients for InterpolateUsingSetup
! InterpolateArrayTeardown     Deallocate tables created by InterpolateArraySetup
! InterpolateValues            Interpolate for new y value(s): given old (x,y), new (x), method
! IsFillValue                  Returns true if argument is FillValue
! ReplaceFillValues            Replaces FillValue entries in an array

  interface BridgeMissingValues
    module procedure BridgeMissingValues_1dr4, BridgeMissingValues_1dr8, BridgeMissingValues_1dint
    module procedure BridgeMissingValues_2dr4, BridgeMissingValues_2dr8, BridgeMissingValues_2dint
    module procedure BridgeMissingValues_3dr4, BridgeMissingValues_3dr8, BridgeMissingValues_3dint
  end interface

  interface Dump
    module procedure DumpCoefficients_r4, DumpCoefficients_r8
  end interface

  interface EssentiallyEqual
    module procedure EssentiallyEqual_r4, EssentiallyEqual_r8
    module procedure EssentiallyEqual_r4_1d, EssentiallyEqual_r8_1d
    module procedure EssentiallyEqual_r4_2d, EssentiallyEqual_r8_2d
    module procedure EssentiallyEqual_r4_3d, EssentiallyEqual_r8_3d
  end interface

  interface Hunt
    module procedure HuntArray_r8, HuntArray_r4
    module procedure HuntScalar_r8, HuntScalar_r4
  end interface

  interface InterpolateArraySetup
    module procedure InterpolateArraySetup_r4, InterpolateArraySetup_r8
  end interface

  interface InterpolateArrayTeardown
    module procedure InterpolateArrayTeardown_r4, InterpolateArrayTeardown_r8
  end interface

  interface InterpolateValues
    module procedure InterpolateArray_r4, InterpolateArray_r8
    module procedure InterpolateScalar_r4, InterpolateScalar_r8
    module procedure InterpolateUsingSetup_r4, InterpolateUsingSetup_r8
    module procedure InterpolateScalarUsingSetup_r4, InterpolateScalarUsingSetup_r8
  end interface

  interface IsFillValue
    module procedure IsFillValue_r4, IsFillValue_r8, IsFillValue_int
  end interface

  interface ReplaceFillValues
    module procedure ReplaceFill1d_r4, ReplaceFill1d_r8, ReplaceFill1d_int
    module procedure ReplaceFill2d_r4, ReplaceFill2d_r8, ReplaceFill2d_int
    module procedure ReplaceFill3d_r4, ReplaceFill3d_r8, ReplaceFill3d_int
  end interface

  real, parameter, private :: FILLVALUETOLERANCE = 0.2 ! Poss. could make it 1

contains

! ---------------------------------------------------------  Dump  -----
  subroutine DumpCoefficients_r4 ( Coeffs )
    use Dump_0, only: Dump
    type(coefficients_r4), intent(in) :: Coeffs
    include 'DumpCoefficients.f9h'
  end subroutine DumpCoefficients_r4

  subroutine DumpCoefficients_r8 ( Coeffs )
    use Dump_0, only: Dump
    type(coefficients_r8), intent(in) :: Coeffs
    include 'DumpCoefficients.f9h'
  end subroutine DumpCoefficients_r8

! ---------------------------------------------  EssentiallyEqual  -----

  ! This family of routines checks to see if two reals are essentially
  ! the same.
  elemental logical function EssentiallyEqual_r4 ( A, B )
    real(r4), intent(in) :: A
    real(r4) ,intent(in) :: B
    EssentiallyEqual_r4 = &
      & a >= nearest ( b, -1.0_r4 ) .and. a <= nearest ( b, 1.0_r4 )
  end function EssentiallyEqual_r4

  elemental logical function EssentiallyEqual_r8 ( A, B )
    real(r8), intent(in) :: A
    real(r8) ,intent(in) :: B
    EssentiallyEqual_r8 = &
      & a >= nearest ( b, -1.0_r8 ) .and. a <= nearest ( b, 1.0_r8 )
  end function EssentiallyEqual_r8

  function EssentiallyEqual_r4_1d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
  ! Warn if an element of one array is finite while the other is not
    real(r4), dimension(:), intent(in)             :: A
    real(r4), dimension(:), intent(in)             :: B
    real(r4), intent(in)                           :: fillValue
    real(r4), dimension(:), optional, intent(in)   :: precision
    logical                                        :: equal
    real(r4), dimension(size(A))                   :: atab
    real(r4), dimension(size(B))                   :: btab
    logical                                        :: warn
    equal = .false.
    call filterValues(A, ATAB, B, BTAB, warn, fillValue, precision)
    if ( .not. warn ) equal = all( &
      & a >= nearest ( b, -1.0_r4 ) .and. a <= nearest ( b, 1.0_r4 ) &
      & )
  end function EssentiallyEqual_r4_1d

  function EssentiallyEqual_r8_1d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
    real(r8), dimension(:), intent(in)             :: A
    real(r8), dimension(:), intent(in)             :: B
    real(r8), intent(in)                           :: fillValue
    real(r8), dimension(:), optional, intent(in)   :: precision
    logical                                        :: equal
    real(r8), dimension(size(A))                   :: atab
    real(r8), dimension(size(B))                   :: btab
    logical                                        :: warn
    equal = .false.
    call filterValues(A, ATAB, B, BTAB, warn, fillValue, precision)
    if ( .not. warn ) equal = all( &
      & a >= nearest ( b, -1.0_r8 ) .and. a <= nearest ( b, 1.0_r8 ) &
      & )
  end function EssentiallyEqual_r8_1d

  function EssentiallyEqual_r4_2d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
    real(r4), dimension(:,:), intent(in)             :: A
    real(r4), dimension(:,:), intent(in)             :: B
    real(r4), intent(in)                             :: fillValue
    real(r4), dimension(:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(2)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)/)) )
    else
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r4_2d

  function EssentiallyEqual_r8_2d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
    real(r8), dimension(:,:), intent(in)             :: A
    real(r8), dimension(:,:), intent(in)             :: B
    real(r8), intent(in)                             :: fillValue
    real(r8), dimension(:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(2)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)/)) )
    else
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)/)), &
        & reshape(b, (/shp(1)*shp(2)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r8_2d

  function EssentiallyEqual_r4_3d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
    real(r4), dimension(:,:,:), intent(in)             :: A
    real(r4), dimension(:,:,:), intent(in)             :: B
    real(r4), intent(in)                               :: fillValue
    real(r4), dimension(:,:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(3)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)*shp(3)/)) )
    else
      equal = EssentiallyEqual_r4_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r4_3d

  function EssentiallyEqual_r8_3d ( A, B, FillValue, Precision ) &
    & result(equal)
  ! This function is slightly different:
  ! A scalar logical testing that every element of two arrays are equal
  ! You may filter out values in either array equal to the FillValue
  ! or for which the corresponding Precision array is negative
  ! or where either element is not finite
    real(r8), dimension(:,:,:), intent(in)             :: A
    real(r8), dimension(:,:,:), intent(in)             :: B
    real(r8), intent(in)                               :: fillValue
    real(r8), dimension(:,:,:), optional, intent(in)   :: precision
    logical                                        :: equal
    ! Internal variables
    integer, dimension(3)                          :: shp
    shp =shape(a)
    if ( present(Precision) ) then
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), &
        & FillValue, reshape(Precision, (/shp(1)*shp(2)*shp(3)/)) )
    else
      equal = EssentiallyEqual_r8_1d(reshape(a, (/shp(1)*shp(2)*shp(3)/)), &
        & reshape(b, (/shp(1)*shp(2)*shp(3)/)), FillValue )
    endif
      
  end function EssentiallyEqual_r8_3d

! -------------------------------------------------  HuntArray_r4  -----

  ! This routine does the classic hunt the value kind of thing.  This does the
  ! hunt/bisect implemention a la Numerical Recipes.  List must be
  ! monotonically increasing or decreasing. There is no such requirements for
  ! values.

  subroutine HuntArray_r4 ( list, values, indices, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace )
    integer, parameter :: RK = R4

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! Result
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    include "HuntArray.f9h"
  end subroutine HuntArray_r4

! -------------------------------------------------  HuntArray_r8  -----

  subroutine HuntArray_r8 ( list, values, indices, start, allowTopValue, allowBelowValue, &
    & nearest, logSpace )
    integer, parameter :: RK = R8

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! Result
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    include "HuntArray.f9h"
  end subroutine HuntArray_r8

! ------------------------------------------------  HuntScalar_r4  -----

  ! This routine is a scalar wrapper for the above one

  subroutine HuntScalar_r4 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace )
    integer, parameter :: RK = R4

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! Resulting index
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    ! Local variables

    integer, dimension(1) :: indices ! To pass to HuntScalar

    call Hunt ( list, (/ value /), indices, start, &
      & allowTopValue, allowBelowValue, nearest )
    index = indices(1)
  end subroutine HuntScalar_r4

! ------------------------------------------------  HuntScalar_r8  -----

  subroutine HuntScalar_r8 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace )
    integer, parameter :: RK = R8

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! Resulting index
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    ! Local variables

    integer, dimension(1) :: indices ! To pass to HuntScalar

    call Hunt ( list, (/ value /), indices, start, &
      & allowTopValue, allowBelowValue, nearest )
    index = indices(1)
  end subroutine HuntScalar_r8

! ------------------------------------------  InterpolateArray_r4  -----

  ! This next subroutine is a workhorse interpolation routine, loosely based on
  ! my (Nathaniel) IDL routine of the same name.

  ! Method is one of 'L'inear, or 'S'pline
  !                                (Numerical Recipes, more later no doubt)
  ! Extrapolate is one of 'A'llow, 'C'onstant or 'B'ad

  ! Notes:
  !   oldX must be monotonically increasing or decreasing
  !   newX can be in any order
  !   one can't ask for spline interpolation with missing regions.
  !   missingRegions will probably slow the code down, as will extrapolate=B

  subroutine InterpolateArray_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld, skipNewY, IntYdX )
    integer, parameter :: RK = R4

    ! Dummy arguments
    real(rk), dimension(:), intent(IN) :: oldX
    real(rk), dimension(:,:), intent(IN) :: oldY
    real(rk), dimension(:), intent(IN) :: newX
    real(rk), dimension(:,:), intent(OUT) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    type (matrixElement_T), intent(out), optional :: dNewByDOld ! Derivatives
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    type(coefficients_r4) :: Coeffs

    include "InterpolateArray.f9h"

  end subroutine InterpolateArray_r4

! ------------------------------------------  InterpolateArray_r8  -----

  subroutine InterpolateArray_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld, skipNewY, IntYdX )
    integer, parameter :: RK = R8

    ! Dummy arguments
    real(rk), dimension(:), intent(IN) :: oldX
    real(rk), dimension(:,:), intent(IN) :: oldY
    real(rk), dimension(:), intent(IN) :: newX
    real(rk), dimension(:,:), intent(OUT) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    type (matrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    type(coefficients_r8) :: Coeffs

    include "InterpolateArray.f9h"

  end subroutine InterpolateArray_r8

! -------------------------------------  InterpolateArraySetup_r4  -----

  subroutine InterpolateArraySetup_r4 ( OldX, NewX, Method, Coeffs, &
    & Extrapolate, Width, DyByDx, dNewByDOld, IntYdX )

    integer, parameter :: RK = R4

    real(rk), intent(in) :: OldX(:), NewX(:)
    character(len=*), intent(in) :: Method
    type(coefficients_R4), intent(out) :: Coeffs
    character(len=*), intent(in), optional :: Extrapolate
    integer, intent(in), optional :: Width ! Second dimension for OldY when
                                           ! interpolations get done
    logical, optional, intent(in) :: DyByDx ! just a signal to
                                           ! compute more coeffs for splines
    type (matrixElement_T), intent(inout), optional :: dNewByDOld ! Derivatives
                                           ! inout so createBlock can clean up
                                           ! after an old one
    logical, optional, intent(in) :: IntYdX ! just a signal to
                                           ! compute more coeffs for splines

    include "InterpolateArraySetup.f9h"

  end subroutine InterpolateArraySetup_r4

! -------------------------------------  InterpolateArraySetup_r8  -----

  subroutine InterpolateArraySetup_r8 ( OldX, NewX, Method, Coeffs, &
    & Extrapolate, Width, DyByDx, dNewByDOld, IntYdX )

    integer, parameter :: RK = R8

    real(rk), intent(in) :: OldX(:), NewX(:)
    character(len=*), intent(in) :: Method
    type(coefficients_R8), intent(out) :: Coeffs
    character(len=*), intent(in), optional :: Extrapolate
    integer, intent(in), optional :: Width ! Second dimension for OldY when
                                           ! interpolations get done
    logical, optional, intent(in) :: DyByDx ! just a signal to
                                           ! compute more coeffs for splines
    type (matrixElement_T), intent(inout), optional :: dNewByDOld ! Derivatives
                                           ! inout so createBlock can clean up
                                           ! after an old one
    logical, optional, intent(in) :: IntYdX ! just a signal to
                                           ! compute more coeffs for splines

    include "InterpolateArraySetup.f9h"

  end subroutine InterpolateArraySetup_r8

! ----------------------------------  InterpolateArrayTeardown_r4  -----

  subroutine InterpolateArrayTeardown_r4 ( Coeffs )

    type(coefficients_R4), intent(inout) :: Coeffs

    include "InterpolateArrayTeardown.f9h"

  end subroutine InterpolateArrayTeardown_r4

! ----------------------------------  InterpolateArrayTeardown_r8  -----

  subroutine InterpolateArrayTeardown_r8 ( Coeffs )

    type(coefficients_R8), intent(inout) :: Coeffs

    include "InterpolateArrayTeardown.f9h"

  end subroutine InterpolateArrayTeardown_r8

! -----------------------------------------  InterpolateScalar_r4  -----

! This subroutine is a scalar wrapper for the array one.

  subroutine InterpolateScalar_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, RangeOfPeriod, skipNewY, IntYdX )
    integer, parameter :: RK = R4

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badValue
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    real(rk), dimension(2), optional, intent(in) :: rangeofperiod	  ! for periodic data
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalar.f9h"
  end subroutine InterpolateScalar_r4

! -----------------------------------------  InterpolateScalar_r8  -----

  subroutine InterpolateScalar_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, RangeOfPeriod, skipNewY, IntYdX )
    integer, parameter :: RK = R8

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badValue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    real(rk), dimension(2), optional, intent(in) :: rangeofperiod	  ! for periodic data
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalar.f9h"
  end subroutine InterpolateScalar_r8

! -------------------------------  InterpolateScalarUsingSetup_r4  -----

  subroutine InterpolateScalarUsingSetup_r4 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = R4

    ! Dummy arguments
    type(coefficients_r4), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalarUsingSetup.f9h"

  end subroutine InterpolateScalarUsingSetup_r4

! -------------------------------  InterpolateScalarUsingSetup_r8  -----

  subroutine InterpolateScalarUsingSetup_r8 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = R8

    ! Dummy arguments
    type(coefficients_r8), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateScalarUsingSetup.f9h"

  end subroutine InterpolateScalarUsingSetup_r8

! -------------------------------------  InterpolateUsingSetup_r4  -----

  subroutine InterpolateUsingSetup_r4 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = R4

    ! Dummy arguments
    type(coefficients_r4), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:,:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:,:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateUsingSetup.f9h"

  end subroutine InterpolateUsingSetup_r4

! -------------------------------------  InterpolateUsingSetup_r8  -----

  subroutine InterpolateUsingSetup_r8 ( coeffs, oldX, oldY, newX, newY, &
    & method, extrapolate, badValue, missingRegions, dyByDx, skipNewY, IntYdX )
    integer, parameter :: RK = R8

    ! Dummy arguments
    type(coefficients_r8), intent(in) :: Coeffs
    real(rk), dimension(:), intent(in) :: oldX
    real(rk), dimension(:,:), intent(in) :: oldY
    real(rk), dimension(:), intent(in) :: newX
    real(rk), dimension(:,:), intent(out) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(rk), optional, intent(in) :: badvalue
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    real(rk), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: SKIPNEWY ! Don't compute newY
    real(rk), dimension(:,:), optional, intent(out) :: IntYdX ! Antiderivative
                                              ! of Y at X

    include "InterpolateUsingSetup.f9h"

  end subroutine InterpolateUsingSetup_r8

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

! ------------------------------------------------- IsFillValue ---

  ! This family of routines checks to see if an arg is a fillValue
  ! This family of routines checks to see if an arg is a fillValue
  elemental logical function IsFillValue_int ( A, FILLVALUE )
    integer, intent(in) :: A
    integer ,intent(in), optional :: FILLVALUE
    integer  :: MYFILLVALUE
    myFillValue = int(DEFAULTUNDEFINEDVALUE)
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_int = &
      & abs(a - myFillValue) < 1 ! FILLVALUETOLERANCE
  end function IsFillValue_int

  elemental logical function IsFillValue_r4 ( A, FILLVALUE )
    real(r4), intent(in) :: A
    real(r4) ,intent(in), optional :: FILLVALUE
    real(r4)  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_r4 = &
      & abs(a - myFillValue) < FILLVALUETOLERANCE
  end function IsFillValue_r4

  elemental logical function IsFillValue_r8 ( A, FILLVALUE )
    real(r8), intent(in) :: A
    real(r8) ,intent(in), optional :: FILLVALUE
    real(r8)  :: MYFILLVALUE
    myFillValue = DEFAULTUNDEFINEDVALUE
    if ( present(fillValue) ) myFillValue = fillValue
    IsFillValue_r8 = &
      & abs(a - myFillValue) < Real(FILLVALUETOLERANCE, r8)
  end function IsFillValue_r8

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

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSNumerics
!=============================================================================

!
! $Log$
! Revision 2.37  2005/08/06 01:36:30  vsnyder
! Add Dump for interpolation coefficients
!
! Revision 2.36  2005/08/05 20:34:47  pwagner
! ReplaceFillValues can now to interpolate to bridge across MissingValues
!
! Revision 2.35  2005/08/03 16:36:46  pwagner
! antiderivatives; added replaceFillValues
!
! Revision 2.34  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.33  2005/05/12 20:47:56  pwagner
! Added new versions of EssentaillyEqual discounting Fills, NaNs
!
! Revision 2.32  2005/01/19 17:16:48  pwagner
! Now gets dump from dump_0
!
! Revision 2.31  2004/09/28 23:14:24  pwagner
! Added isFillValue function
!
! Revision 2.30  2004/09/10 23:52:29  livesey
! Added the logSpace options for Hunt (not actually needed but never
! mind).
!
! Revision 2.29  2003/09/11 23:09:18  livesey
! Added skipNewY argument
!
! Revision 2.28  2003/04/04 00:10:08  livesey
! Added EssentiallyEqual
!
! Revision 2.27  2002/11/25 18:51:19  vsnyder
! More interfaces
!
! Revision 2.26  2002/11/23 00:01:00  vsnyder
! Modify interpolation to separate setup and teardown
!
! Revision 2.25  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.24  2002/10/04 16:40:30  pwagner
! Fixed missing close quotes on include lines
!
! Revision 2.23  2002/10/04 01:48:27  vsnyder
! Get rid of a local variable
!
! Revision 2.22  2002/10/04 00:48:05  vsnyder
! Move declarations of local variables to includes.  Use generic
! InterpolateValues in InterpolateScalar routines.  Move guts of
! InterpolateScalar routines to an include file.
!
! Revision 2.21  2002/09/13 18:08:12  pwagner
! May change matrix precision rm from r8
!
! Revision 2.20  2002/09/11 17:43:38  pwagner
! Began changes needed to conform with matrix%values type move to rm from r8
!
! Revision 2.19  2002/05/24 17:00:55  livesey
! Fixed bug with interpolation with only one oldX value
!
! Revision 2.18  2001/12/01 01:03:07  livesey
! Bug fix with derivatives.
!
! Revision 2.17  2001/11/14 01:47:40  livesey
! Added nearest option to Hunt
!
! Revision 2.16  2001/10/19 23:41:43  livesey
! Fixed bug with dyByDx extrapolation
!
! Revision 2.15  2001/09/24 17:27:50  pwagner
! Removed random number things
!
! Revision 2.14  2001/09/20 20:54:55  pwagner
! Added random number stuff from MATH77
!
! Revision 2.13  2001/07/10 17:55:23  livesey
! Cosmetic changes only
!
! Revision 2.12  2001/07/06 18:44:40  dwu
! Add a feature in InterpolateScaler to handle periodic data
!
! Revision 2.11  2001/06/07 21:59:41  pwagner
! Added Copyright statement
!
! Revision 2.10  2001/05/08 23:25:58  livesey
! Added a nullify for oldSecond (how did I miss that?)
!
! Revision 2.9  2001/05/03 23:13:28  livesey
! Made Van's changes compile with NAG
!
! Revision 2.8  2001/05/03 21:54:49  vsnyder
! Added some nullify's, did some cosmetic changes
!
! Revision 2.7  2001/04/28 19:42:48  livesey
! Hunt now correctly handles cases where list is of size one.
! Changed to lower case keywords and reformatted.
!
! Revision 2.6  2001/04/28 07:05:05  livesey
! Minor bug fix in spline
!
! Revision 2.5  2001/04/11 22:43:19  vsnyder
! Let sparsify do the deallocate_test
!
! Revision 2.4  2001/03/06 00:35:23  livesey
! Missed one pointer nullification
!
! Revision 2.3  2001/03/05 01:20:36  livesey
! Nullified pointers
!
! Revision 2.2  2001/02/22 01:59:52  vsnyder
! Remove declarations for unused variables and parameters
!
! Revision 2.1  2001/02/09 00:38:55  livesey
! Various changes
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.13  2000/06/23 01:08:48  vsnyder
! Delete unused variables (except ID) to keep NAG f95 happy
!
