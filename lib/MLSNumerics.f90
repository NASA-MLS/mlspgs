! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSNumerics              ! Some low level numerical stuff
  !=============================================================================

  use MLSCommon, only : R8
  use MLSMessageModule, only: MLSMessage,MLSMSG_Error
  use MLSStrings, only: Capitalize
  use MatrixModule_0, only: MatrixElement_T,CreateBlock_0,M_Column_Sparse, Sparsify
  use Allocate_Deallocate, only : Allocate_test, Deallocate_test

  implicit none

  private
  public :: Hunt, InterpolateValues, srang, drang

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
! drang                        Returns a gauss. dis. dbl.pr. num. 0 mean, 1 s.d.
! srang                        Returns a gauss. dis. sngl..pr. num. 0 mean, 1 s.d.


  interface Hunt
    module procedure HuntArray
    module procedure HuntScalar
  end interface

  interface InterpolateValues
    module procedure InterpolateArray
    module procedure InterpolateScalar
  end interface

contains

  ! This routine does the classic hunt the value kind of thing.  This does the
  ! hunt/bisect implemention a la Numerical Recipes.  List must be
  ! monotonically increasing or decreasing. There is no such requirements for
  ! values.

  subroutine HuntArray ( list, values, indices, start, allowTopValue, allowBelowValue )

    ! Dummy arguments
    real(r8), dimension(:), intent(in) :: list ! List to search
    real(r8), dimension(:), intent(in) :: values ! Values to search for
    integer, dimension(:), intent(out) :: indices ! Result
    integer, optional, intent(in) :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0

    ! Local variables
    integer :: listLen, valuesLen ! Array sizes
    integer :: valueIndex       ! Loop counters
    integer :: index            ! Temporary result

    logical :: useAllowTopValue, useAllowBelowValue
    integer :: useStart
    integer :: upperLimit       ! Highest value that can be returned
    integer :: stride           ! Value to step by
    logical :: expanding        ! Whether we're expanding or reducing our search
    logical :: lowerBelow       ! Flag
    logical :: upperAbove       ! Another flag

    integer :: listDirection    ! +1 if list ascends, -1 descends
    integer :: searchDirection  ! (in index space) 
    integer :: oldSearchDirection ! Previous value of above

    real(r8) :: thisValue

    ! Executable code

    if ( present(allowTopValue) ) then
      useAllowTopValue = allowTopValue
    else
      useAllowTopValue = .false.
    end if

    if ( present(allowBelowValue) ) then
      useAllowBelowValue = allowBelowValue
    else
      useAllowBelowValue = .false.
    end if

    if ( present(start) ) then
      useStart = start
    else
      useStart = 1
    end if

    listLen = size(list)
    valuesLen = size(values)
    if ( size(indices) < valuesLen ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Result array is too small" )

    ! Try to work out the direction, also skip if there's only one value

    if ( size(list) == 1 ) then
      indices = 1
      if ( useAllowBelowValue ) then
        where ( values < list(1) )
          indices = 0
        end where
      end if
      return
    else
      if ( list(size(list)) >= list(1) ) then
        listDirection = 1
      else
        listDirection = -1
      end if
    end if

    ! Some last bits of setup before we get going.

    if ( useAllowTopValue ) then
      upperLimit = listLen
    else
      upperLimit = listLen-1
    end if

    ! Now we're ready to hit the road, loop over all the values to hunt for

    index = max(1,min(useStart,upperLimit))
    do valueIndex = 1, valuesLen
      thisValue = values(valueIndex)
      expanding = .true.
      searchDirection = 0
      stride = 1
      HuntLoop: do
        lowerBelow= (thisValue-list(index))*listDirection >= 0.0
        if ( index<listLen ) then 
          upperAbove = (list(index+1)-thisValue)*listDirection > 0.0
        else ! We're off the end of the list
          upperAbove = .true.   
        end if

        ! Now we know what the state of play is, what does it mean?

        ! First see if we've found the place
        if ( lowerBelow.and.upperAbove ) exit HuntLoop

        ! The other cases are a little more complex
        oldSearchDirection = searchDirection

        if ( lowerBelow.and. (.not. upperAbove) ) then
          ! If we're at the end, get out
          if ( index == upperLimit ) exit HuntLoop
          ! We're too low, keep looking upwards
          index = index+stride
          searchDirection = 1
        end if

        if ( (.not. lowerBelow).and.upperAbove ) then
          ! If we're at the begning, get out
          if ( index == 1 ) exit HuntLoop

          ! We're too high but not at begining, look back downwards
          index = index-stride
          searchDirection = -1
        end if

        ! Now the very first change of direction is the end of hte
        ! `expanding' phase

        if ( (searchDirection /= oldSearchDirection) .and. &
          & (oldSearchDirection /= 0) .and. (expanding) ) &
          & expanding = .false.

        if ( expanding ) then
          stride = min(stride*2,listLen/2)
        else
          stride = max(stride/2,1)
        end if

        ! Make sure we don't fall off an end

        index = min(max(index,1),upperLimit)
      end do HuntLoop

      ! Final check for off the bottom of the list

      if ( useAllowBelowValue ) then
        if ( (thisValue-list(index))*listDirection<0.0D0) index = 0
      end if

      indices(valueIndex) = index
    end do
  end subroutine HuntArray

  ! ---------------------------------------------------------------------------

  ! This routine is a scalar wrapper for the above one

  subroutine HuntScalar (list, value, index, start, allowTopValue, allowBelowValue )

    ! Dummy arguments
    real(r8), dimension(:), intent(in) :: list ! List to search
    real(r8), intent(in) :: value ! Value to search for
    integer, intent(out) :: index ! Resulting index
    integer, intent(in), optional :: start ! Optional start index
    logical, optional, intent(in) :: allowTopValue ! Can return N
    logical, optional, intent(in) :: allowBelowValue ! Can return 0

    ! Local variables

    real(R8), dimension(1) :: values ! To pass to HuntArray
    integer, dimension(1) :: indices ! To pass to HuntScalar

    values(1) = value
    call HuntArray ( list, values, indices, start, allowTopValue, allowBelowValue )
    index = indices(1)
  end subroutine HuntScalar

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

  subroutine InterpolateArray ( oldX, oldY, newX, newY, method, extrapolate, &
    & badValue, missingRegions, dyByDx, dNewByDOld )

    ! Dummy arguments
    real(R8), dimension(:), intent(IN) :: oldX
    real(R8), dimension(:,:), intent(IN) :: oldY
    real(R8), dimension(:), intent(IN) :: newX
    real(R8), dimension(:,:), intent(OUT) :: newY

    character (len=*), intent(in) :: method ! See comments above
    character (len=*), optional, intent(in) :: extrapolate ! See comments above
    real(r8), optional, intent(in) :: badvalue
    real(r8), dimension(:,:), optional, intent(out) :: dyByDx
    logical, optional, intent(in) :: missingRegions ! Allow missing regions
    type (matrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives

    ! Local variables
    integer :: noOld, noNew, width ! Dimensions
    logical :: spline              ! Flag
    logical :: useMissingRegions   ! Copy of missing regions
    integer :: ind, newInd         ! Loop counters
    logical :: computeDNewByDOld   ! Set if dNewByDOld is present

    real(R8), dimension(:),   pointer :: A
    real(R8), dimension(:,:), pointer :: AA
    real(R8), dimension(:),   pointer :: B
    real(R8), dimension(:,:), pointer :: BB
    real(R8), dimension(:),   pointer :: C
    real(R8), dimension(:,:), pointer :: CC
    real(R8), dimension(:),   pointer :: D ! Coefficients
    real(R8), dimension(:,:), pointer :: DD ! Spread coefs.
    real(R8), dimension(:),   pointer :: gap
    real(R8), dimension(:),   pointer :: gap2
    integer, dimension(:),    pointer :: lowerInds
    real(R8), dimension(:),   pointer :: maskVector
    real(R8), dimension(:,:), pointer :: oldSecond
    real(R8), dimension(:,:), pointer :: oldSecondLower
    real(R8), dimension(:,:), pointer :: oldSecondUpper
    real(R8), dimension(:,:), pointer :: oldYlower
    real(R8), dimension(:,:), pointer :: oldYupper
    real(R8), dimension(:),   pointer :: p ! For 2nd der. guess
    real(R8), dimension(:,:), pointer :: spreadGap
    real(R8), dimension(:,:), pointer :: temp ! For 2nd der. guess
    real(R8), dimension(:,:), pointer :: tempDNewByDOld ! Dense version.
    integer, dimension(:),    pointer :: upperInds
    real(R8) :: sig       ! For second derivative guesser

    character :: extrapolateMethod ! Tidy copy of extrapolate parameter

    ! Executable code

    nullify ( a, aa, b, bb, c, cc, d, dd, gap, gap2, lowerInds, maskVector, &
      &       oldSecond, oldSecondLower, oldSecondUpper, oldYlower, oldYupper, p, &
      &       spreadGap, temp, tempDNewByDOld, upperInds )

    ! Size the problem, check sanity, set up arrays etc.

    noOld=size(oldX,1)
    noNew=size(newX,1)
    width=size(oldY,2)

    spline=(Capitalize(method(1:1))=="S")

    extrapolateMethod="A"
    if ( present(extrapolate)) extrapolateMethod=Capitalize(extrapolate(1:1))

    useMissingRegions=.false.
    if ( present(missingRegions)) useMissingRegions=missingRegions

    computeDNewByDOld=present(dNewByDOld)

    if ( useMissingRegions.and.spline ) call MLSMessage &
      & ( MLSMSG_Error, ModuleName, "Cannot use missing regions with spline")

    if ( computeDNewByDOld .and. spline ) call MLSMessage &
      & ( MLSMSG_Error, ModuleName, "Cannot get dNewBydOld from spline")

    call Allocate_Test ( lowerInds, noNew, "lowerInds", ModuleName )
    call Allocate_Test ( upperInds, noNew, "upperInds", ModuleName )
    call Allocate_Test ( gap, noNew, "gap", ModuleName )
    call Allocate_Test ( A, noNew, "A", ModuleName )
    call Allocate_Test ( B, noNew, "B", ModuleName )
    call Allocate_Test ( AA, noNew, width, "AA", ModuleName )
    call Allocate_Test ( BB, noNew, width, "BB", ModuleName )
    call Allocate_Test ( oldYlower, noNew, width, "oldYlower", ModuleName )
    call Allocate_Test ( oldYupper, noNew, width, "oldYupper", ModuleName )

    ! Setup arrays needed if dyByDx is requested

    if ( present(dyByDx)) call Allocate_Test(spreadGap,noNew,width,&
      "spreadGap",ModuleName)

    ! Setup Matrix block needed if DNewByDOld is needed.
    if ( computeDNewByDOld ) then
      call CreateBlock_0 ( dNewByDOld, noNew*width, noOld*width, &
        M_Column_Sparse, NumberNonZero=2*noNew*width )
    end if

    ! Do special stuff for the case of spline, allocate arrays, find 2nd
    ! derivatives etc.

    if ( spline ) then
      call Allocate_Test ( oldSecond, noOld, width, "oldSecond", ModuleName )
      call Allocate_Test ( C, noNew, "C", ModuleName )
      call Allocate_Test ( D, noNew, "D", ModuleName )
      call Allocate_Test ( CC, noNew, width, "CC", ModuleName )
      call Allocate_Test ( DD, noNew, width, "DD", ModuleName )
      call Allocate_Test ( oldSecondlower, noNew, width, "oldSecondlower", ModuleName )
      call Allocate_Test ( oldSecondupper, noNew, width, "oldSecondupper", ModuleName )
      call Allocate_Test ( gap2, noNew, "gap2", ModuleName )
      call Allocate_Test ( temp, noOld, width, "temp", ModuleName )
      call Allocate_Test ( p, width, "p", ModuleName )

      ! Here we have to solve the a tridiagonal equation
      ! This is a straight copy of my idl code
      oldSecond(1,:) = 0.0
      temp(1,:) = 0.0
      do ind = 2, noOld-1
        sig = (oldX(ind)-oldX(ind-1))/(oldX(ind+1)-oldX(ind-1))
        p = sig*oldSecond(ind-1,:)+2.0
        oldSecond(ind,:) = (sig-1.0)/p
        temp(ind,:) = (oldY(ind+1,:)-oldY(ind,:))/(oldX(ind+1)-oldX(ind)) - &
          & (oldY(ind,:)-oldY(ind-1,:))/(oldX(ind)-oldX(ind-1))
        temp(ind,:) = (6.0*temp(ind,:)/ &
          & (oldX(ind+1)-oldX(ind-1))-sig*temp(ind-1,:))/p
      end do
      oldSecond(noOld,:) = 0.0

      ! Now do the back substitution
      do ind = noOld-1, 1, -1
        oldSecond(ind,:) = oldSecond(ind,:)*oldSecond(ind+1,:)+temp(ind,:)
      end do

      call Deallocate_test ( temp, "Temp", ModuleName ) 
      call Deallocate_test ( p, "p", ModuleName )
    end if

    ! Now we're ready to begin the real work.

    ! Clear the result array(s)
    newY = 0.0
    if ( present(dyByDx) ) dyByDx = 0.0

    ! Now hunt for the indices

    call Hunt ( oldX, newX, lowerInds )
    upperInds = lowerInds+1
    gap = oldX(upperInds)-oldX(lowerInds)
    if ( present(dyByDx) ) spreadGap = spread(gap,2,width)

    A = (oldX(upperInds)-newX)/gap

    ! If extrapolate is "C"onstant, deal with that
    if ( extrapolateMethod=="C" ) A = max(min(A,1.0_r8),0.0_r8)

    B=1.0_r8-A

    ! If extrapolate mode is "B"ad, deal with that
    if ( extrapolateMethod=="B" ) then
      call Allocate_Test ( maskVector, noNew, "maskVector", ModuleName )
      maskVector = 0.0
      where ( (A<0.0) .or. (A>1.0) )
        maskVector = badValue
        A = 0.0
        B = 0.0
      end where
      newY = spread(maskVector,2,width)
      if ( present(dyByDx) ) dyByDx = newY
      call Deallocate_Test ( maskVector, "maskVector", ModuleName )
    end if

    ! Now spread out the coefficients
    AA = spread(A,2,width)
    BB = spread(B,2,width)
    oldYlower = oldY(lowerInds,:)
    oldYupper = oldY(upperInds,:)

    ! Now worry about the missing regions flag
    if ( useMissingRegions ) then
      where( (oldYlower==badValue) .or. (oldYupper==badValue))
        newY = badValue
        AA = 0.0
        BB = 0.0
      end where
      if ( present(dyByDx) ) then
        where( (oldYlower==badValue) .or. &
          & (oldYupper==badValue) )
        dyByDx = badValue
        oldYlower = 0.0      ! Only way to guarentee bad derivative
        oldYupper  =0.0      ! But don't need to worry about spline
      end where
    end if
  end if

  ! Now do the linear interpolation calculation
  newY = newY+AA*oldYlower+BB*oldYupper
  if ( present(dyByDx) ) dyByDx = (oldYupper-oldYlower)/spreadGap

  ! Write the output derivative matrix if needed
  if ( computeDNewByDOld ) then
    ! While the matrix is ideally suited to row sparse, our storage method
    ! is column sparse, so to be lazy we'll create it full and then sparsify
    ! it.

    call Allocate_Test ( tempDNewByDOld, noNew*width, noOld*width, &
      & "tempDNewByDOld", ModuleName )
    do newInd = 1, noNew
      do ind = 1, width
        tempDNewByDOld(newInd+ind*noNew,lowerInds(newInd)+ind*noOld) = A(newInd)
        tempDNewByDOld(newInd+ind*noNew,upperInds(newInd)+ind*noOld) = B(newInd)
      end do
    end do
    call Sparsify ( tempDNewByDOld, dNewbyDOld, &
      & "tempDNewByDOld", ModuleName ) ! dNewbyDOld := tempDNewByDOld
  end if

  ! Now do the spline calculation
  if ( spline ) then
    gap2 = gap**2
    C = (A**3-A)*gap2/6.0    ! Note the extrapolate bad case is covered as..
    D = (B**3-B)*gap2/6.0    !   A=B=0.0

    ! Spread out the coefficients etc.
    CC = spread(C,2,width)
    DD = spread(D,2,width)
    oldSecondLower = oldSecond(lowerInds,:)
    oldSecondUpper = oldSecond(upperInds,:)

    newY = newY+CC*oldSecondLower+DD*oldSecondUpper
    if ( present(dyByDx)) dyByDx=dyByDx+(spreadGap/6.0)*( &
      & (3.0*BB**2-1.0)*oldSecondUpper- &
      & (3.0*AA**2-1.0)*oldSecondLower)
  end if

  call Deallocate_Test ( lowerInds, "lowerInds", ModuleName )
  call Deallocate_Test ( upperInds, "upperInds", ModuleName )
  call Deallocate_Test ( gap, "gap", ModuleName )
  call Deallocate_Test ( A, "A", ModuleName )
  call Deallocate_Test ( B, "B", ModuleName )
  call Deallocate_Test ( AA, "AA", ModuleName )
  call Deallocate_Test ( BB, "BB", ModuleName )
  call Deallocate_Test ( oldYlower, "oldYlower", ModuleName )
  call Deallocate_Test ( oldYupper, "oldYupper", ModuleName )

  if ( spline ) then
    call Deallocate_Test( oldSecond, "oldSecond", ModuleName )
    call Deallocate_Test( C, "C", ModuleName )
    call Deallocate_Test( D, "D", ModuleName )
    call Deallocate_Test( CC, "CC", ModuleName )
    call Deallocate_Test( DD, "DD", ModuleName )
    call Deallocate_Test( oldSecondlower, "oldSecondlower", ModuleName )
    call Deallocate_Test( oldSecondupper, "oldSecondupper", ModuleName )
    call Deallocate_Test( gap2, "gap2", ModuleName )
    call Deallocate_Test( temp, "temp", ModuleName )
    call Deallocate_Test( p, "p", ModuleName )
  end if
  if ( present(dyByDx) ) &
    & call Deallocate_Test ( spreadGap, "spreadGap", ModuleName )

end subroutine InterpolateArray

! --------------------------------------------------------------------------

! This subroutine is a scalar wrapper for the first one.

subroutine InterpolateScalar ( oldX, oldY, newX, newY, method, extrapolate, &
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

    call InterpolateArray ( oldX, spread(tempY,2,1), newX, tempResult, method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions, dyByDx=tempDerivative )
    dyByDx = tempDerivative(:,1)

    call Deallocate_Test ( tempDerivative, "tempDerivative", ModuleName )
  else
    call InterpolateArray ( oldX, spread(tempY,2,1), newX, tempResult, method, &
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
end subroutine InterpolateScalar

!------------------------------------------------------------------
!   Random number routines from MATH77 libraries
! ../l2:69% ls *.f
!  drang.f  ranpk1.f  ranpk2.f  srang.f
!  ../l2:71% cat *.f > stuff.sed
!  ../l2:72% sed 's/^[Cc]/\!/' stuff.sed > sed.out
! plus a small amount of subsequent editing
!------------------------------------------------------------------
      double precision function  DRANG ()
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!>> 1996-04-16 DRANG WVS SQRT(abs(TWO*log(U3))) avoids sqrt(-0.0)
!>> 1994-10-20 DRANG Krogh  Changes to use M77CON
!>> 1994-06-24 DRANG CLL Changed common to use RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-16 DRANG CLL
!>> 1991-11-26 DRANG CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 DRANG CLL Added call to RAN0, and DGFLAG in common.
!>> 1991-01-15 DRANG CLL Reordered common contents for efficiency.
!>> 1990-01-23 DRANG CLL Making names in common same in all subprogams.
!>> 1987-06-09 DRANG CLLawson  Initial code.
!
!     Returns one pseudorandom number from the Gausian (Normal)
!     distribution with zero mean and unit standard deviation.
!     Method taken from Algorithm 334, Comm. ACM, July 1968, p. 498.
!     Implemented at JPL in Univac Fortran V in 1969 by Wiley R. Bunton
!     of JPL and Stephen L. Richie of Heliodyne Corp.
!
!     Adapted to Fortran 77 for the MATH 77 library by C. L. Lawson and
!     S. Y. Chiu, JPL, April 1987, 6/9/87.
!     ------------------------------------------------------------------
!--D replaces "?": ?RANG, ?RANUA, RANC?1, RANC?2, ?PTR, ?NUMS, ?GFLAG
!     RANCD1 and RANCD2 are common blocks.
!     Uses intrinsic functions, log and sqrt.
!     Calls DRANUA to obtain an array of uniform random numbers.
!     Calls RAN0 to initialize DPTR and DGFLAG.
!     ------------------------------------------------------------------
!                        Important common variables
!
!     DPTR [integer] and DGFLAG [logical]
!
!          Will be set to DPTR = 1 and DGFLAG = .false.
!          when RAN0 is called from this subr if this
!          is the first call to RAN0.  Also set to these values when
!          RAN1 or RANPUT is called by a user.
!
!          DGFLAG will be set true on return from this subr when the
!          algorithm has values that it can save to reduce the amount
!          of computation needed the next time this function is
!          referenced.  Will be set false in the contrary case.
!
!          DPTR is the index of the last location used in the
!          common array DNUMS().  This index is counted down.
!
!     DNUMS() [floating point]  Buffer of previously computed uniform
!          random numbers.
!     ------------------------------------------------------------------
      integer M
      parameter(M = 97)
      double precision    ONE, TWO
      parameter(ONE = 1.0D0, TWO = 2.0D0)
      double precision    DNUMS(M), R, S, U3, X, XX, Y, YY
      logical FIRST
      common/RANCD2/DNUMS
      integer DPTR
      logical DGFLAG
      common/RANCD1/DPTR, DGFLAG
      save  /RANCD1/, /RANCD2/, FIRST
      save    R, X, Y
      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
         FIRST = .false.
         call RAN0
      endif
!
      if (.not. DGFLAG .or. DPTR .eq. 1) then
!
!     Use the Von Neuman rejection method for choosing a random point
!     (X,Y) in the unit circle, X**2 + Y**2 .le. 1.0.
!     Then the angle Theta = arctan(Y/X) is random, uniform in
!     (-Pi/2, Pi/2), and Phi = 2*Theta is random, uniform in (-Pi, Pi).
!     Define S = X**2 + Y**2, then
!     sin(Theta) = Y/sqrt(S),    cos(Theta) = X/sqrt(S),
!     sin(Phi) = 2*X*Y/S,    and cos(Phi) = (X**2 - Y**2)/S.
!
   10    continue
!                              Set X = random, uniform in [0., 1.]
      DPTR = DPTR - 1
      if(DPTR .eq. 0) then
         call DRANUA(DNUMS, M)
         DPTR = M
      endif
            X = DNUMS(DPTR)
!                              Set Y = random, uniform in [-1., 1.]
            DPTR = DPTR - 1
            if(DPTR .eq. 0) then
               call DRANUA(DNUMS,M)
               DPTR = M
            endif
            Y = TWO*DNUMS(DPTR) - ONE
!
            XX=X*X
            YY=Y*Y
            S=XX+YY
         if(S .gt. ONE) go to 10
!
!     Choose R randomly from Chi-squared distribution and
!     normalize with S.
!
!                              Set U3 = random, uniform in [0., 1.]
         DPTR = DPTR - 1
         if(DPTR .eq. 0) then
            call DRANUA(DNUMS,M)
            DPTR = M
         endif
         U3 = DNUMS(DPTR)
!        Changed -TWO*log(U3) to abs(TWO*log(U3)) because Lahey LF90
!        2.00 on a pentium produced -0.0 for -TWO*log(1.0), then got a
!        floating point exception on sqrt(-0.0).
         R = sqrt(abs(TWO*(log(U3))))/S
!
!                                Compute result as  R*Sin(PHI)
!
         DRANG = (XX-YY)*R
         DGFLAG = .true.
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when DGFLAG is true and DPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      DRANG = TWO*X*Y*R
      DGFLAG=.false.
      return
      end function  DRANG
      subroutine RAN1
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!                            Program unit: RANPK1
!>> 1995-11-21 RAMPK1 Krogh Removed multiple entries.
!>> 1994-06-24 CLL Reorganized common. Using RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-13 CLL Fixed error in RAN0
!>> 1991-11-26 CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 CLL Added Entry RAN0 and common variables SGFLAG,DGFLAG
!>> 1991-01-15 CLL Reordered common contents for efficiency.
!>> 1990-01-23 CLL Corrected type stmt for SNUMS in common.
!>> 1987-04-22 RANPK1 Lawson  Initial code.
!
!        This program unit, RANPK1, along with RANPK2,
!     supports random number generation.
!
!        This prog unit has entries RAN1, RAN0, and RANPUT.
!     The library user can call RAN1 to initialize random number
!     generation at a standard initial seed,
!     or call RANPUT(KSEED) to initialize random number generation
!     at a seed value provided by the integer array, KSEED().
!
!     Other higher level random number subrs call RAN0 on their first
!     time flags to be sure the package is initialized.
!
!     As a result of any of these entries this subroutine will
!     set the pointers in the COMMON arrays to 1, indicating to higher
!     level random number subprograms that these buffer arrays are
!     empty.  It also sets SGFLAG and DGFLAG to .false. to indicate to
!     Gaussian generators that they have no internal saved value.
!
!     The user can determine the appropriate dimension for the array,
!     KSEED() by first calling the entry RANSIZ in prog unit RANPK2.
!
!     The user can retrieve the current seed value by calling entry,
!     RANGET in prog unit RANPK2.  This will be the seed that will be
!     used the next time a batch of random numbers are computed.  This
!     is not necessarily the seed associated with the next number that
!     will be returned.
!     C. L. Lawson, F. T. Krogh, & S. Y. Chiu, JPL, Apr 1987.
!     ------------------------------------------------------------------
!
      integer DPTR, SPTR
      logical DGFLAG, SGFLAG
      common/RANCD1/DPTR, DGFLAG
      common/RANCS1/SPTR, SGFLAG
      save  /RANCD1/, /RANCS1/
!     ------------------------------------------------------------------
!                      For use by library users: CALL RAN1
      call RN1
      DPTR = 1
      SPTR = 1
      DGFLAG = .false.
      SGFLAG = .false.
      return
      end subroutine RAN1
!     ------------------------------------------------------------------
!                      For use by other library subprograms: CALL RAN0
      subroutine RAN0
      integer DPTR, SPTR
      logical DGFLAG, SGFLAG
      common/RANCD1/DPTR, DGFLAG
      common/RANCS1/SPTR, SGFLAG
      save  /RANCD1/, /RANCS1/
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if(FIRST) then
         FIRST = .false.
         DPTR = 1
         SPTR = 1
         DGFLAG = .false.
         SGFLAG = .false.
      end if
      return
      end subroutine RAN0
!     ------------------------------------------------------------------
!                         For use by library users: CALL RANPUT(KSEED)
      subroutine RANPUT(KSEED)
      integer KSEED(*)
      integer DPTR, SPTR
      logical DGFLAG, SGFLAG
      common/RANCD1/DPTR, DGFLAG
      common/RANCS1/SPTR, SGFLAG
      save  /RANCD1/, /RANCS1/
!
      call   RNPUT(KSEED)
      DPTR = 1
      SPTR = 1
      DGFLAG = .false.
      SGFLAG = .false.
      return
      end subroutine RANPUT
      subroutine RN1
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!>> 1997-12-17 RANPK2 Krogh Removed unreferenced labels
!>> 1995-11-22 RANPK2 Krogh Removed multiple entries.
!>> 1992-03-17 CLL Moved SAVE stmt ahead of DATA stmts.
!>> 1992-03-09 CLL Removed "save FIRST" because FIRST is not used.
!>> 1992-03-02 CLL Fix error: Set MODE = 1 in data stmt.
!>> 1991-11-21 CLL Add MODE with values 1, 2, 3, & 4.
!>> 1989-09-11 CLL Multiversion file. RANPK2 or RANPK3
!>> 1987-05-05 RANPK2 Lawson  Initial code.
!
!        This program unit, along with RANPK1, supports random number
!     generation.
!
!     The functionality of this random number package was modeled on the
!     specifications of the random number routines RANDOM and RANDOMSEED
!     in the February, 1987 working draft of the Fortran 8x language
!     standard.  This functionality remains similar to the later draft,
!     Fortran 90, S8.115, June 1990, in which the routine names have
!     been changed to RANDOM_MUMBER and RANDOM_SEED.  This should
!     facilitate replacement of use of this package by Fortran intrinsic
!     when and if Fortran 90 compilers come into widespread use.
!
!     The library user may call RANSIZ or RANGET in this prog unit,
!     or RAN1 or RANPUT in RANPK1 to obtain the functionality of
!     RANDOM_SEED of Fortran 90.  This relates to setting or fetching
!     the seed.
!        Entries RN1 and RNPUT in this prog unit should not be called by
!     library users.  They are intended only to be called from RANPK1.
!        Entry RN2 returns the value of MODE.  This is a convenience in
!     case one is interested in knowing the value of MODE the package
!     has selected.
!        Entries SRANUA and DRANUA (s.p. and d.p. respectively)
!     generate arrays of pseudorandom numbers from the uniform
!     distribution in the range [0., 1.].  These may be called by users
!     and are called by other library routines.
!        Entries SRANUS and DRANUA (s.p. and d.p. respectively)
!     generate arrays of pseudorandom numbers from the uniform
!     distribution in the range [0., 1.] and then transformed to
!     A + B*U.  These are intended for direct use by users.
!     ------------------------------------------------------------------
!              Algorithm for generating pseudorandom numbers
!                     from the uniform distribution.
!
!  The current integer value in the random integer sequence is XCUR,
!  and the next is defined mathematically by the statement:
!
!                 XCUR = mod(AFAC * XCUR,  MDIV)
!  where
!                 MDIV = m = 6_87194_76503 = 2**36 - 233
!  and
!                 AFAC = a = 612_662 = (approx.) 0.58 * 2**20
!
!  XCUR may be any integer value in the range from 1 to m-1, and all
!  integer values in this range will be generated before the sequence
!  cycles.
!
!  We call the above computational algorithm for XCUR the "short"
!  algorithm.  There is also a "long" algorithm that produces exactly
!  the same sequence of values of XCUR:
!                  Q = aint(XCUR/B)
!                  R = XCUR - Q * B
!                  XCUR = AFAC * R - C * Q
!                  do while(XCUR .lt. 0.0)
!                     XCUR = XCUR + MDIV
!                  end do
!  where B and C are constants related to MDIV and AFAC by
!            MDIV = B * AFAC + C
!  We use B = 112165 and C = 243273.  The average number of executions
!  of the statement XCUR = XCUR + MDIV is 1.09 and the maximum number of
!  executions is 3.
!
!  The largest number that must be handled in the "short" algorithm
!  is the product of AFAC  with the max value of XCUR, i.e.,
!    612_662 * 6_87194_76502 = 42_10181_19126_68324 ~= 0.58 * 2**56.
!  Thus this algorithm requires arithmetic exact to at least 56 bits.
!
!  The largest number that must be handled in the "long" algorithm
!  is the product of C with the max value of aint(XCUR/B), i.e.,
!               243273 * 612664 ~= 0.14904e12 ~= 0.54 * 2**38.
!  Thus this algorithm requires arithmetic exact to at least 38 bits.
!
!  To accommodate different compiler/computer systems this program unit
!  contains code for 3 different ways of computing the new XCUR from the
!  old XCUR, each producing exactly the same sequence of of XCUR.
!
!  Initially we have MODE = 1.  When MODE = 1 the code does tests to
!  see which of three implementation methods will be used, and sets
!  MODE = 2, 3, or 4 to indicate the choice.
!
!  Mode 2 will be used in machines such as the Cray that have at
!  least a 38 bit significand in SP arithmetic.  XCUR will be advanced
!  using the "long" algorithm in SP arithmetic.
!
!  Mode 3 will be used on machines that don't meet the Mode 2 test,
!  but can maintain at least a 56 bits exactly in computing
!  mod(AFAC*XCUR, MDIV) in DP arithmetic.  This includes VAX, UNISYS,
!  IBM 30xx, and some IEEE machines that have clever compilers that
!  keep an extended precision representation of the product AFAC*XCUR
!  within the math processor for use in the division by MDIV.  XCUR will
!  be advanced using the "short" algorithm in DP arithmetic.
!
!  Mode 4 will be used on machines that don't meet the Mode 2 or 3
!  tests, but have at least a 38 bit significand in DP arithmetic.
!  This includes IEEE machines that have not-so-clever compilers.
!  XCUR is advanced using the "long" algorithm in DP arithmetic.
!  ---------------------------------------------------------------------
!               Properties of the generated sequence.
!
!        This m is one of the prime numbers cited in
!     Table 1, Page 390, of Knuth, Seminumerical Algorithms, Second
!     edition, 1981, Addison-Wesley.
!     The prime factorization of m-1 is
!           m-1 = p1 * p2 * p3 = 2 * 43801 * 784_451
!     The complementary factors are
!           q(1) = 3_43597_38251, q(2) = 15_68902, and q(3) = 87602.
!
!     The value a is a primitive root of m as is verified by
!     computing a**q(i) mod m for i = 1,3, and finding these values are
!     not 1.  These values are m-1, 2_49653_21011, and 1_44431_31136.
!     The fact that a is a primitive root of m assures that the period
!     of the generator is m-1, i.e. starting with any integer from 1
!     through m-1, all integers in this range will be produced.
!
!     The value a has relatively large values of the measures nu and mu
!     computed for the Spectral Test as described in Knuth, pp. 89-105.
!        (Log10(nu(i)), i=2,6) = 5.4, 3.6,  2.6,  2.2,  1.8
!        (mu(i), i=2,6)        = 3.0, 3.05, 3.39, 4.55, 6.01
!     This assures that the generated sequence will have relatively low
!     autocorrelation.
!     ------------------------------------------------------------------
!                          Alternative algorithm
!
!  An alternative set of constants that has been used widely in
!  commercial and public domain software packages is
!               m = 21474_83647 = 2**31 - 1
!               a = 16807       = 7**5 = (approx.) 0.513 * 2**15
!
!  The largest product that must be handled exactly is approximately
!  0.513 * 2**46 which is approximately  0.36E14.  This is within the
!  double-precision capability of most computer systems.
!
!  The sequence can be started with any integer from 1 through m-1
!  and will generate all integers in this range.  The autocorrelation
!  properties of the whole sequence will not be as good as with the
!  larger values for m and a.
!     ------------------------------------------------------------------
!     C. L. Lawson, F. T. Krogh & S. Chiu, JPL, July 1986, Apr 1987.
!     ------------------------------------------------------------------
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!              These same parameters are also defined below in RANMOD.
      double precision X1DP
      real             X1SP
      parameter( X1DP=123456789.0D0, X1SP=123456789.0e0 )
!     ------------------------------------------------------------------
!                         Entered using CALL RN1
!              This entry should not be called by general users.
!              User should call RAN1 in RANPK1.
!
         XCURDP = X1DP
         XCURSP = X1SP
         return
         end subroutine RN1
!     ------------------------------------------------------------------
      subroutine RANSIZ(KSIZE)
         integer KSIZE
         KSIZE = 2
         return
         end subroutine RANSIZ
!     ------------------------------------------------------------------
!                         Entered using CALL RNPUT(KSEED)
!              This entry should not be called by general users.
!              User should call RANPUT(KSEED) in RANPK1.
!
      subroutine RNPUT(KSEED)
      integer KSEED(2)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      double precision MDIVDP, SCALDP
      real             MDIVSP, SCALSP
      parameter( MDIVDP=68719476503.0D0, MDIVSP=68719476503.0e0 )
      parameter( SCALDP=100000.0D0, SCALSP=100000.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      if(MODE .eq. 2) then
         XCURSP = SCALSP * real(abs(KSEED(1))) + real(abs(KSEED(2)))
         XCURSP = max(1.e0, min(XCURSP, MDIVSP - 1.0e0))
      else
!                            Here for MODE = 3 or 4
         XCURDP = SCALDP * dble(abs(KSEED(1))) + dble(abs(KSEED(2)))
         XCURDP = max(1.D0, min(XCURDP, MDIVDP - 1.0D0))
      end if
      return
      end subroutine RNPUT
!     ------------------------------------------------------------------
      subroutine RANGET(KSEED)
      integer KSEED(2)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      double precision SCALDP
      real             SCALSP
      parameter( SCALDP=100000.0D0, SCALSP=100000.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      if(MODE .eq. 2) then
         KSEED(1) = int(XCURSP/SCALSP)
         KSEED(2) = int(XCURSP - SCALSP * real(KSEED(1)))
      else
!                            Here for MODE = 3 or 4
!
         KSEED(1) = int(XCURDP/SCALDP)
         KSEED(2) = int(XCURDP - SCALDP * dble(KSEED(1)))
      end if
      return
      end subroutine RANGET
!     ------------------------------------------------------------------
      subroutine RN2(MODE1)
      integer MODE1
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      MODE1 = MODE
      return
      end subroutine RN2
!     ------------------------------------------------------------------
      subroutine SRANUA(USP, N)
      integer N
      real USP(N)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      go to (310, 320, 330, 340), MODE
  310 stop'In file RANPK2, subroutine SRANUA -- Ivalid value for MODE'
!                                         Mode 2.
  320   continue
      do 325 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  322    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 322
         end if
         USP(I) =  XCURSP/MDIVSP
  325 continue
      go to 350
!                                         Mode 3.
  330   continue
      do 335 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         USP(I) = real( XCURDP)/MDIVSP
  335 continue
      go to 350
!                                         Mode 4.
  340   continue
      do 345 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  342    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 342
         end if
         USP(I) = real( XCURDP)/MDIVSP
  345 continue
  350 continue
      return
      end subroutine SRANUA
!     ------------------------------------------------------------------
      subroutine DRANUA(UDP, N)
      integer N
      double precision UDP(N)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      go to (410, 420, 430, 440), MODE
  410 stop'In file RANPK2, subroutine DRANUA -- Ivalid value for MODE'
!                                         Mode 2.
  420   continue
      do 425 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  422    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 422
         end if
         UDP(I) = dble(XCURSP) / MDIVDP
  425 continue
      go to 450
!                                         Mode 3.
  430   continue
      do 435 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         UDP(I) =  XCURDP/MDIVDP
  435 continue
      go to 450
!                                         Mode 4.
  440   continue
      do 445 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  442    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 442
         end if
         UDP(I) =  XCURDP/MDIVDP
  445 continue
  450 continue
      return
      end subroutine DRANUA
!     ------------------------------------------------------------------
      subroutine SRANUS(USP, N, ASP, BSP)
      integer N
      real ASP, BSP, USP(N)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      go to (510, 520, 530, 540), MODE
  510 stop'In file RANPK2, subroutine SRANUS -- Ivalid value for MODE'
!                                         Mode 2.
  520   continue
      do 525 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  522    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 522
         end if
         USP(I) = ASP + BSP * XCURSP/MDIVSP
  525 continue
      go to 550
!                                         Mode 3.
  530   continue
      do 535 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         USP(I) = ASP + BSP * real( XCURDP)/MDIVSP
  535 continue
      go to 550
!                                         Mode 4.
  540   continue
      do 545 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  542    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 542
         end if
         USP(I) = ASP + BSP * real( XCURDP)/MDIVSP
  545 continue
  550 continue
      return
      end subroutine SRANUS
!     ------------------------------------------------------------------
      subroutine DRANUS(UDP, N, ADP, BDP)
      integer N
      double precision ADP, BDP, UDP(N)
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
      logical FIRST
      save FIRST
      data FIRST / .true. /
!
      if (FIRST) then
         FIRST = .false.
         call RANMOD
      end if
      go to (610, 620, 630, 640), MODE
  610 stop'In file RANPK2, subroutine DRANUS -- Ivalid value for MODE'
!                                         Mode 2.
  620   continue
      do 625 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  622    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 622
         end if
         UDP(I) = ADP + BDP * dble(XCURSP) / MDIVDP
  625 continue
      go to 650
!                                         Mode 3.
  630   continue
      do 635 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         UDP(I) = ADP + BDP * XCURDP/MDIVDP
  635 continue
      go to 650
!                                         Mode 4.
  640   continue
      do 645 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  642    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 642
         end if
         UDP(I) = ADP + BDP * XCURDP/MDIVDP
  645 continue
  650 continue
      return
      end subroutine DRANUS
!     ------------------------------------------------------------------
      subroutine RANMOD
!
!     Do tests to decide whether to set the MODE = 2, 3, or 4.
!     Outcome will depend on precision of SP and DP floating point
!     arithmetic on the host system.
!                    Common Block
      integer MODE
      real XCURSP
      double precision XCURDP
      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I, J, K
      double precision TEST(0:3,0:1), DIFF
      double precision QDP, RDP, XTSTDP
      double precision AFACDP, BFACDP, C2DP, MDIVDP, X1DP
      real             QSP, RSP, XTSTSP
      real             AFACSP, BFACSP, C2SP, MDIVSP, X1SP
      parameter( AFACDP=612662.0D0, MDIVDP=68719476503.0D0, &
     & X1DP=123456789.0D0, BFACDP=112165.0d0, C2DP=243273.0d0)
      parameter( AFACSP=612662.0e0, MDIVSP=68719476503.0e0, &
     &  X1SP=123456789.0e0, BFACSP=112165.0e0, C2SP=243273.0e0)
      logical DONE
      save DONE
!
      data DONE / .false. /
      data TEST / 24997965550.0d0, 68719476502.0d0, &
     &            68718863841.0d0, 36962132774.0d0, &
     &            43721510953.0d0,           1.0d0, &
     &                 612662.0d0, 31757343729.0d0/
!
      if (DONE) return
      DONE = .true.
      do 880 MODE = 2, 4
         do 870 J = 0,1
            XTSTDP = TEST(0,J)
            XTSTSP = real(XTSTDP)
            do 860 I = 1,3
               go to (820, 830, 840),MODE-1
!                                                Test of MODE 2
  820          continue
                  QSP = aint(XTSTSP / BFACSP)
                  RSP = XTSTSP - QSP * BFACSP
                  XTSTSP = AFACSP * RSP - C2SP * QSP
                  do 822 K = 1,3
                     if(XTSTSP .ge. 0.0e0) go to 825
                       XTSTSP = XTSTSP + MDIVSP
  822               continue
  825               continue
                  DIFF = dble(XTSTSP) - TEST(I,J)
               go to 850
!                                                Test of MODE 3
  830          continue
                  XTSTDP = mod(AFACDP * XTSTDP,  MDIVDP)
                  DIFF = XTSTDP - TEST(I,J)
               go to 850
!                                                Test of MODE 4
  840          continue
                  QDP = aint(XTSTDP / BFACDP)
                  RDP = XTSTDP - QDP * BFACDP
                  XTSTDP = AFACDP * RDP - C2DP * QDP
                  do 842 K = 1,3
                     if(XTSTDP .ge. 0.0d0) go to 845
                     XTSTDP = XTSTDP + MDIVDP
  842             continue
  845             continue
                  DIFF = XTSTDP - TEST(I,J)
  850          continue
!*              print'(1x,a,3i3,g11.3)', 'RANPK2.. MODE,J,I,DIFF=',
!*    *         MODE,J,I,DIFF                !****** For Testing ******
              print *, 'RANPK2.. MODE,J,I,DIFF'
              print *, MODE,J,I,DIFF
               if(DIFF .ne. 0.0d0) go to 880
!                            Following line ends I loop.
  860         continue
!                            Following line ends J loop.
  870      continue
!
!        Here the computations using the current value of MODE have
!        passed all tests, so we accept this value of MODE.
         XCURDP = X1DP
         XCURSP = X1SP
         return
!
!*        print*,'From RANPK2.. MODE =',MODE  !***** For Testing *****
          print*, 'From RANPK2.. MODE =',MODE  !***** For Testing *****
!                            Following line ends MODE loop.
  880 continue
!        The computations were unsuccessful for all values of MODE.
!        This means this random number package will not work on the
!        current host system.  ****** Fatal Error Stop ******
!
!      call ERMSG('RANPK2',1, 2,
!     *'This random no. code will not work on this computer system.','.')
      print *, '***RANPK2 wishes to report the following error'
      print *, 'This rnc will not work on this computer system.'
      return
      end subroutine RANMOD

      real             function  SRANG ()
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!>> 1996-04-16 SRANG WVS SQRT(abs(TWO*log(U3))) avoids sqrt(-0.0)
!>> 1994-10-20 SRANG Krogh  Changes to use M77CON
!>> 1994-06-24 SRANG CLL Changed common to use RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-16 SRANG CLL
!>> 1991-11-26 SRANG CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 SRANG CLL Added call to RAN0, and SGFLAG in common.
!>> 1991-01-15 SRANG CLL Reordered common contents for efficiency.
!>> 1990-01-23 SRANG CLL Making names in common same in all subprogams.
!>> 1987-06-09 SRANG CLLawson  Initial code.
!
!     Returns one pseudorandom number from the Gausian (Normal)
!     distribution with zero mean and unit standard deviation.
!     Method taken from Algorithm 334, Comm. ACM, July 1968, p. 498.
!     Implemented at JPL in Univac Fortran V in 1969 by Wiley R. Bunton
!     of JPL and Stephen L. Richie of Heliodyne Corp.
!
!     Adapted to Fortran 77 for the MATH 77 library by C. L. Lawson and
!     S. Y. Chiu, JPL, April 1987, 6/9/87.
!     ------------------------------------------------------------------
!--S replaces "?": ?RANG, ?RANUA, RANC?1, RANC?2, ?PTR, ?NUMS, ?GFLAG
!     RANCS1 and RANCS2 are common blocks.
!     Uses intrinsic functions, log and sqrt.
!     Calls SRANUA to obtain an array of uniform random numbers.
!     Calls RAN0 to initialize SPTR and SGFLAG.
!     ------------------------------------------------------------------
!                        Important common variables
!
!     SPTR [integer] and SGFLAG [logical]
!
!          Will be set to SPTR = 1 and SGFLAG = .false.
!          when RAN0 is called from this subr if this
!          is the first call to RAN0.  Also set to these values when
!          RAN1 or RANPUT is called by a user.
!
!          SGFLAG will be set true on return from this subr when the
!          algorithm has values that it can save to reduce the amount
!          of computation needed the next time this function is
!          referenced.  Will be set false in the contrary case.
!
!          SPTR is the index of the last location used in the
!          common array SNUMS().  This index is counted down.
!
!     SNUMS() [floating point]  Buffer of previously computed uniform
!          random numbers.
!     ------------------------------------------------------------------
      integer M
      parameter(M = 97)
      real                ONE, TWO
      parameter(ONE = 1.0E0, TWO = 2.0E0)
      real                SNUMS(M), R, S, U3, X, XX, Y, YY
      logical FIRST
      common/RANCS2/SNUMS
      integer SPTR
      logical SGFLAG
      common/RANCS1/SPTR, SGFLAG
      save  /RANCS1/, /RANCS2/, FIRST
      save    R, X, Y
      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
         FIRST = .false.
         call RAN0
      endif
!
      if (.not. SGFLAG .or. SPTR .eq. 1) then
!
!     Use the Von Neuman rejection method for choosing a random point
!     (X,Y) in the unit circle, X**2 + Y**2 .le. 1.0.
!     Then the angle Theta = arctan(Y/X) is random, uniform in
!     (-Pi/2, Pi/2), and Phi = 2*Theta is random, uniform in (-Pi, Pi).
!     Define S = X**2 + Y**2, then
!     sin(Theta) = Y/sqrt(S),    cos(Theta) = X/sqrt(S),
!     sin(Phi) = 2*X*Y/S,    and cos(Phi) = (X**2 - Y**2)/S.
!
9910    continue
!                              Set X = random, uniform in [0., 1.]
      SPTR = SPTR - 1
      if(SPTR .eq. 0) then
         call SRANUA(SNUMS, M)
         SPTR = M
      endif
            X = SNUMS(SPTR)
!                              Set Y = random, uniform in [-1., 1.]
            SPTR = SPTR - 1
            if(SPTR .eq. 0) then
               call SRANUA(SNUMS,M)
               SPTR = M
            endif
            Y = TWO*SNUMS(SPTR) - ONE
!
            XX=X*X
            YY=Y*Y
            S=XX+YY
         if(S .gt. ONE) go to 9910
!
!     Choose R randomly from Chi-squared distribution and
!     normalize with S.
!
!                              Set U3 = random, uniform in [0., 1.]
         SPTR = SPTR - 1
         if(SPTR .eq. 0) then
            call SRANUA(SNUMS,M)
            SPTR = M
         endif
         U3 = SNUMS(SPTR)
!        Changed -TWO*log(U3) to abs(TWO*log(U3)) because Lahey LF90
!        2.00 on a pentium produced -0.0 for -TWO*log(1.0), then got a
!        floating point exception on sqrt(-0.0).
         R = sqrt(abs(TWO*(log(U3))))/S
!
!                                Compute result as  R*Sin(PHI)
!
         SRANG = (XX-YY)*R
         SGFLAG = .true.
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when SGFLAG is true and SPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      SRANG = TWO*X*Y*R
      SGFLAG=.false.
      return
      end function  SRANG

!=============================================================================
end module MLSNumerics
!=============================================================================

!
! $Log$
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
