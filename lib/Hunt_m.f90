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
module Hunt_m
!=============================================================================

  use Diff_1, only: SelfDiff
  use Dump_0, only: Dump
  use MLSCommon, only : UndefinedValue
  use MLSFillValues, only: IsFillValue, ReplaceFillValues
  use MLSFinds, only: FindFirst, FindLast
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, &
    & MLSMessage

  implicit none

  private

  public :: Hunt, HuntBox, HuntRange, PureHunt

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------


  interface Hunt
    module procedure HuntArray_r4, HuntArray_r8
    module procedure HuntScalar_r4, HuntScalar_r8
  end interface

  interface HuntBox
    module procedure HuntBox_r4, HuntBox_r8
  end interface

  interface HuntRange
    module procedure HuntRange_int, HuntRange_r4, HuntRange_r8
  end interface

  interface purehunt
    module procedure purehunt_r4, purehunt_r8
  end interface

contains

! -------------------------------------------------  HuntArray_r4  -----

  ! This routine does the classic hunt the value kind of thing.  This
  ! does the hunt/bisect implemention a la Numerical Recipes.  List must
  ! be monotonically increasing or decreasing. There is no such
  ! requirement for values.

  subroutine HuntArray_r4 ( list, values, indices, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace, fail )
    use ieee_arithmetic, only: IEEE_Is_NaN
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! list(indices) <= values
      !                                               <= list(indices+1)
      !                                      1 <= indices < N unless
      !                                      allowTopValue or allowBelowValue
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space
    logical, optional, intent(out) :: Fail    ! True for failure

    include "HuntArray.f9h"
  end subroutine HuntArray_r4

! -------------------------------------------------  HuntArray_r8  -----

  subroutine HuntArray_r8 ( list, values, indices, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace, fail )
    use ieee_arithmetic, only: IEEE_Is_NaN
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! list(indices) <= values
      !                                               <= list(indices+1)
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space
    logical, optional, intent(out) :: Fail    ! True for failure

    include "HuntArray.f9h"
  end subroutine HuntArray_r8

! ------------------------------------------------  HuntScalar_r4  -----

  ! This routine is a scalar wrapper for the above one

  subroutine HuntScalar_r4 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace )
    integer, parameter :: RK = kind(0.0e0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! list(index) <= value <= list(index+1)
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    ! Local variables

    integer, dimension(1) :: indices ! To pass to HuntArray

    call Hunt ( list, (/ value /), indices, start, &
      & allowTopValue, allowBelowValue, nearest, logSpace )
    index = indices(1)
  end subroutine HuntScalar_r4

! ------------------------------------------------  HuntScalar_r8  -----

  subroutine HuntScalar_r8 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest, logSpace )
    integer, parameter :: RK = kind(0.0d0)

    ! Dummy arguments
    real(rk), dimension(:), intent(in) :: list ! List to search
    real(rk), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! list(index) <= value <= list(index+1)
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead
    logical, optional, intent(in) :: logSpace ! Choose nearest based on log space

    ! Local variables

    integer, dimension(1) :: indices ! To pass to HuntArray

    call Hunt ( list, (/ value /), indices, start, &
      & allowTopValue, allowBelowValue, nearest, logSpace )
    index = indices(1)
  end subroutine HuntScalar_r8

  ! ----------------------------------------------------  HuntBox  -----
  ! A binary search routine with a hunt procedure, to start from last known
  ! location (if 0 < JLO < N) or from the begining otherwise.
  subroutine HuntBox_r4 ( GRIDPOINTS, MGRIDPOINTS, COORDS, INDICES, VERTICES )
    integer, parameter :: RK = kind(0.0e0)
    include 'HuntBox.f9h'
  end subroutine HuntBox_r4

  subroutine HuntBox_r8 ( GRIDPOINTS, MGRIDPOINTS, COORDS, INDICES, VERTICES )
    integer, parameter :: RK = kind(0.0d0)
    include 'HuntBox.f9h'
  end subroutine HuntBox_r8

! ----------------------------------------------------  HuntRange  -----
! This family of subroutines searchs not for a single index
! but for a range of values within which all list elements
! lie, inclusive
! If none, return (/ 0, 0 /)
! Special interpretation: inclusive means
! if vrange(1) == vrange(2), any values of list also == vrange
! are within that range
! As with other Hunts, list must be monotonic
  subroutine HuntRange_int ( list, vrange, irange, options )
    integer, parameter :: RK = kind(0.0e0)
    ! Dummy args
    integer, dimension(:) :: list
    integer, dimension(2) :: vrange
    include 'HuntRange.f9h'
  end subroutine HuntRange_int

  subroutine HuntRange_r4 ( list, vrange, irange, options )
    integer, parameter :: RK = kind(0.0e0)
    ! Dummy args
    real(rk), dimension(:) :: list
    real(rk), dimension(2) :: vrange
    include 'HuntRange.f9h'
  end subroutine HuntRange_r4

  subroutine HuntRange_r8 ( list, vrange, irange, options )
    integer, parameter :: RK = kind(0.0d0)
    ! Dummy args
    real(rk), dimension(:) :: list
    real(rk), dimension(2) :: vrange
    include 'HuntRange.f9h'
  end subroutine HuntRange_r8

  ! ---------------------------------------------------  PureHunt  -----
  ! A binary search routine with a hunt procedure, to start from last known
  ! location (if 0 < JLO < N) or from the begining otherwise.
  pure subroutine purehunt_r4 ( ELEMENT, ARRAY, N, JLO, JHI )
    integer, parameter :: RK = kind(0.0e0)
    include 'hunt.f9h'
  end subroutine purehunt_r4

  pure subroutine purehunt_r8 ( ELEMENT, ARRAY, N, JLO, JHI )
    integer, parameter :: RK = kind(0.0d0)
    include 'hunt.f9h'
  end subroutine purehunt_r8

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Hunt_m
!=============================================================================

! $Log$
! Revision 2.2  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.1  2015/05/27 22:36:35  vsnyder
! Initial commit -- moved stuff here from MLSNumerics to avoid circularity
!
