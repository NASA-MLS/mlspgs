! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSNumerics              ! Some low level numerical stuff
  !=============================================================================

  use MLSCommon, only : r8, rm, r4
  use MLSMessageModule, only: MLSMessage,MLSMSG_Error
  use MLSStrings, only: Capitalize
  use MatrixModule_0, only: MatrixElement_T,CreateBlock_0,M_Column_Sparse,M_Absent,Sparsify
  use Allocate_Deallocate, only : Allocate_test, Deallocate_test

  implicit none

  private
  public :: Hunt, InterpolateValues

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! This module contains some low level numerical stuff, hunting, interpolating
  ! etc.
!     c o n t e n t s
!     - - - - - - - -

!         Functions, operations, routines
! Hunt                         Finds index of item(s) in list closest to prey
! HuntArray                    Hunts for multiple items
! HuntScalar                   Hunts for just one
! InterpolateValues            Interpolate for new y value(s): given old (x,y), new (x), method
! InterpolateArray             Interpolates for multiple values
! InterpolateScalar            Interpolates for just one


  interface Hunt
    module procedure HuntArray_r8, HuntArray_r4
    module procedure HuntScalar_r8, HuntScalar_r4
  end interface

  interface InterpolateValues
    module procedure InterpolateArray_r4, InterpolateArray_r8
    module procedure InterpolateScalar_r4, InterpolateScalar_r8
  end interface

contains

  ! This routine does the classic hunt the value kind of thing.  This does the
  ! hunt/bisect implemention a la Numerical Recipes.  List must be
  ! monotonically increasing or decreasing. There is no such requirements for
  ! values.

  subroutine HuntArray_r8 ( list, values, indices, start, allowTopValue, allowBelowValue, &
    & nearest )

    ! Dummy arguments
    real(r8), dimension(:), intent(in) :: list ! List to search
    real(r8), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! Result
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below

    ! Local variables
    real(r8) :: thisValue
    include "HuntArray.f9h"
  end subroutine HuntArray_r8

  subroutine HuntArray_r4 ( list, values, indices, start, allowTopValue, allowBelowValue, &
    & nearest )

    ! Dummy arguments
    real(r4), dimension(:), intent(in) :: list ! List to search
    real(r4), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! Result
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value not one below

    ! Local variables
    real(r4) :: thisValue
    include "HuntArray.f9h"
  end subroutine HuntArray_r4

  ! ---------------------------------------------------------------------------

  ! This routine is a scalar wrapper for the above one

  subroutine HuntScalar_r8 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest )

    ! Dummy arguments
    real(r8), dimension(:), intent(in) :: list ! List to search
    real(r8), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! Resulting index
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead

    ! Local variables

    real(r8), dimension(1) :: values ! To pass to HuntArray
    integer, dimension(1) :: indices ! To pass to HuntScalar

    values(1) = value
    call HuntArray_r8 ( list, values, indices, start, &
      & allowTopValue, allowBelowValue, nearest )
    index = indices(1)
  end subroutine HuntScalar_r8

  subroutine HuntScalar_r4 (list, value, index, start, allowTopValue, &
    & allowBelowValue, nearest )

    ! Dummy arguments
    real(r4), dimension(:), intent(in) :: list ! List to search
    real(r4), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! Resulting index
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0
    logical, optional, intent(in) :: nearest ! Choose nearest value instead

    ! Local variables

    real(r4), dimension(1) :: values ! To pass to HuntArray
    integer, dimension(1) :: indices ! To pass to HuntScalar

    values(1) = value
    call HuntArray_r4 ( list, values, indices, start, &
      & allowTopValue, allowBelowValue, nearest )
    index = indices(1)
  end subroutine HuntScalar_r4

  ! ---------------------------------------------------------------------------

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

  subroutine InterpolateArray_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld )

    ! Dummy arguments
    real(r8), dimension(:), intent(IN) :: oldX
    real(r8), dimension(:,:), intent(IN) :: oldY
    real(r8), dimension(:), intent(IN) :: newX
    real(r8), dimension(:,:), intent(OUT) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(r8), optional, intent(in) :: badvalue
    real(r8), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    type (matrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives

    ! Local variables
    real(r8), dimension(:),   pointer :: A
    real(r8), dimension(:,:), pointer :: AA
    real(r8), dimension(:),   pointer :: B
    real(r8), dimension(:,:), pointer :: BB
    real(r8), dimension(:),   pointer :: C
    real(r8), dimension(:,:), pointer :: CC
    real(r8), dimension(:),   pointer :: D ! Coefficients
    real(r8), dimension(:,:), pointer :: DD ! Spread coefs.
    real(r8), dimension(:),   pointer :: gap
    real(r8), dimension(:),   pointer :: gap2
    integer, dimension(:),    pointer :: lowerInds
    real(r8), dimension(:),   pointer :: maskVector
    real(r8), dimension(:,:), pointer :: oldSecond
    real(r8), dimension(:,:), pointer :: oldSecondLower
    real(r8), dimension(:,:), pointer :: oldSecondUpper
    real(r8), dimension(:,:), pointer :: oldYlower
    real(r8), dimension(:,:), pointer :: oldYupper
    real(r8), dimension(:),   pointer :: p ! For 2nd der. guess
    real(r8), dimension(:,:), pointer :: spreadGap
    real(r8), dimension(:,:), pointer :: temp ! For 2nd der. guess
    real(rm), dimension(:,:), pointer :: tempDNewByDOld ! Dense version.
    integer, dimension(:),    pointer :: upperInds
    real(r8) :: sig       ! For second derivative guesser
    real(r8) :: dyByDxFill              ! Fill value for dyByDx
    include "InterpolateArray.f9h"

end subroutine InterpolateArray_r8

  subroutine InterpolateArray_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld )

    ! Dummy arguments
    real(r4), dimension(:), intent(IN) :: oldX
    real(r4), dimension(:,:), intent(IN) :: oldY
    real(r4), dimension(:), intent(IN) :: newX
    real(r4), dimension(:,:), intent(OUT) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(r4), optional, intent(in) :: badvalue
    real(r4), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    type (matrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives

    ! Local variables
    real(r4), dimension(:),   pointer :: A
    real(r4), dimension(:,:), pointer :: AA
    real(r4), dimension(:),   pointer :: B
    real(r4), dimension(:,:), pointer :: BB
    real(r4), dimension(:),   pointer :: C
    real(r4), dimension(:,:), pointer :: CC
    real(r4), dimension(:),   pointer :: D ! Coefficients
    real(r4), dimension(:,:), pointer :: DD ! Spread coefs.
    real(r4), dimension(:),   pointer :: gap
    real(r4), dimension(:),   pointer :: gap2
    integer, dimension(:),    pointer :: lowerInds
    real(r4), dimension(:),   pointer :: maskVector
    real(r4), dimension(:,:), pointer :: oldSecond
    real(r4), dimension(:,:), pointer :: oldSecondLower
    real(r4), dimension(:,:), pointer :: oldSecondUpper
    real(r4), dimension(:,:), pointer :: oldYlower
    real(r4), dimension(:,:), pointer :: oldYupper
    real(r4), dimension(:),   pointer :: p ! For 2nd der. guess
    real(r4), dimension(:,:), pointer :: spreadGap
    real(r4), dimension(:,:), pointer :: temp ! For 2nd der. guess
    real(rm), dimension(:,:), pointer :: tempDNewByDOld ! Dense version.
    integer, dimension(:),    pointer :: upperInds
    real(r4) :: sig       ! For second derivative guesser
    real(r4) :: dyByDxFill              ! Fill value for dyByDx
    include "InterpolateArray.f9h"

end subroutine InterpolateArray_r4

! --------------------------------------------------------------------------

! This subroutine is a scalar wrapper for the first one.

subroutine InterpolateScalar_r8 ( oldX, oldY, newX, newY, method, extrapolate, &
  & badValue, missingRegions, dyByDx, RangeOfPeriod )

  ! Dummy arguments
  real(r8), dimension(:), intent(in) :: oldX
  real(r8), dimension(:), intent(in) :: oldY
  real(r8), dimension(:), intent(in) :: newX
  real(r8), dimension(:), intent(out) :: newY

  character (len=*), intent(in) :: method ! See comments above
  character (len=*), optional, intent(in) :: extrapolate ! See comments above
  real(r8), optional, intent(in) :: badValue
  real(r8), dimension(:), optional, intent(out) :: dyByDx
  real(r8), dimension(2), optional, intent(in) :: rangeofperiod	  ! for periodic data
  logical, optional, intent(in) :: missingRegions ! Allow missing regions

! local working space
  real(r8), dimension(:,:), pointer :: tempDerivative
  real(r8), dimension(size(newX), 1) :: tempResult
  real(r8), dimension(size(oldY)) :: tempY
  real(r8) period
  integer jump, j

  ! Executable code

  tempY = oldY

  if ( present(rangeofperiod) ) then
	period  = rangeofPeriod(2)-rangeofPeriod(1)
	jump = -1
	do j =1, size(oldY)-1
		if(abs(tempY(j+1)-tempY(j)) > period/2. ) jump = j 
	enddo
	if(jump /= -1) then
	   if(tempY(jump+1) > tempY(jump)) then
		tempY(jump+1:) = tempY(jump+1:) - period
	   else 
		tempY(jump+1:) = tempY(jump+1:) + period
	   end if
	end if
  end if

  nullify ( tempDerivative )

  if ( present(dyByDx) ) then
    call Allocate_Test ( tempDerivative, size(newX), 1, &
      & "tempDerivative", ModuleName )

    call InterpolateArray_r8 ( oldX, spread(tempY,2,1), newX, tempResult, method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions, dyByDx=tempDerivative )
    dyByDx = tempDerivative(:,1)

    call Deallocate_Test ( tempDerivative, "tempDerivative", ModuleName )
  else
    call InterpolateArray_r8 ( oldX, spread(tempY,2,1), newX, tempResult, method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions )
  end if
  newY = tempResult(:,1)

  if ( present(rangeofperiod) ) then
	period  = rangeofPeriod(2)-rangeofPeriod(1)
	where (newY > rangeofperiod(2)) 
	  newY = newY - period
	elsewhere (newY < rangeofperiod(1)) 
	  newY = newY + period
	end where
  end if
end subroutine InterpolateScalar_r8

subroutine InterpolateScalar_r4 ( oldX, oldY, newX, newY, method, extrapolate, &
  & badValue, missingRegions, dyByDx, RangeOfPeriod )

  ! Dummy arguments
  real(r4), dimension(:), intent(in) :: oldX
  real(r4), dimension(:), intent(in) :: oldY
  real(r4), dimension(:), intent(in) :: newX
  real(r4), dimension(:), intent(out) :: newY

  character (len=*), intent(in) :: method ! See comments above
  character (len=*), optional, intent(in) :: extrapolate ! See comments above
  real(r4), optional, intent(in) :: badValue
  real(r4), dimension(:), optional, intent(out) :: dyByDx
  real(r4), dimension(2), optional, intent(in) :: rangeofperiod	  ! for periodic data
  logical, optional, intent(in) :: missingRegions ! Allow missing regions

! local working space
  real(r4), dimension(:,:), pointer :: tempDerivative
  real(r4), dimension(size(newX), 1) :: tempResult
  real(r4), dimension(size(oldY)) :: tempY
  real(r4) period
  integer jump, j

  ! Executable code

  tempY = oldY

  if ( present(rangeofperiod) ) then
	period  = rangeofPeriod(2)-rangeofPeriod(1)
	jump = -1
	do j =1, size(oldY)-1
		if(abs(tempY(j+1)-tempY(j)) > period/2. ) jump = j 
	enddo
	if(jump /= -1) then
	   if(tempY(jump+1) > tempY(jump)) then
		tempY(jump+1:) = tempY(jump+1:) - period
	   else 
		tempY(jump+1:) = tempY(jump+1:) + period
	   end if
	end if
  end if

  nullify ( tempDerivative )

  if ( present(dyByDx) ) then
    call Allocate_Test ( tempDerivative, size(newX), 1, &
      & "tempDerivative", ModuleName )

    call InterpolateArray_r4 ( oldX, spread(tempY,2,1), newX, tempResult, method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions, dyByDx=tempDerivative )
    dyByDx = tempDerivative(:,1)

    call Deallocate_Test ( tempDerivative, "tempDerivative", ModuleName )
  else
    call InterpolateArray_r4 ( oldX, spread(tempY,2,1), newX, tempResult, method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions )
  end if
  newY = tempResult(:,1)

  if ( present(rangeofperiod) ) then
	period  = rangeofPeriod(2)-rangeofPeriod(1)
	where (newY > rangeofperiod(2)) 
	  newY = newY - period
	elsewhere (newY < rangeofperiod(1)) 
	  newY = newY + period
	end where
  end if
end subroutine InterpolateScalar_r4

!=============================================================================
end module MLSNumerics
!=============================================================================

!
! $Log$
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
