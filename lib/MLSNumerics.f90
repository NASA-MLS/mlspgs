!=============================================================================
module MLSNumerics              ! Some low level numerical stuff
  !=============================================================================

  use MLSCommon, only : R8
  use MLSMessageModule, only: MLSMessage,MLSMSG_Error
  use MLSStrings, only: Capitalize
  use MatrixModule_0, only: MatrixElement_T,CreateBlock_0,M_Column_Sparse, Sparsify
  use Allocate_Deallocate, only : Allocate_test, Deallocate_test

  implicit none

  public

  private :: Id,ModuleName
  !------------------------------- RCS Ident Info ----------------------------
  character(LEN=256) :: Id = & 
    "$Id$"
  character(LEN=*), parameter :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

!
  ! This module contains some low level numerical stuff, hunting, interpolating
  ! etc.
!

  interface Hunt
    module procedure HuntArray
    module procedure HuntScalar
  end interface

  interface InterpolateValues
    module procedure InterpolateArray
    module procedure InterpolateScalar
  end interface

  private :: InterpolateArray, InterpolateScalar, HuntArray, HuntScalar

contains

  ! This routine does the classic hunt the value kind of thing.  This does the
  ! hunt/bisect implemention a la Numerical Recipes.  List must be
  ! monotonically increasing or decreasing. There is no such requirements for
  ! values.

  subroutine HuntArray(list,values,indices,start,allowTopValue,allowBelowValue)

    ! Dummy arguments
    real(R8), dimension(:), intent(IN) :: list ! List to search
    real(R8), dimension(:), intent(IN) :: values ! Values to search for
    integer, dimension(:), intent(OUT) :: indices ! Result
    integer, optional, intent(IN) :: start ! Optional start index
    logical, optional, intent(IN) :: allowTopValue ! Can return N
    logical, optional, intent(IN) :: allowBelowValue ! Can return 0

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

    real(R8) :: thisValue

    ! Executable code

    if (present(allowTopValue)) then
      useAllowTopValue=allowTopValue
    else
      useAllowTopValue=.false.
    endif

    if (present(allowBelowValue)) then
      useAllowBelowValue=allowBelowValue
    else
      useAllowBelowValue=.false.
    endif

    if (present(start)) then
      useStart=start
    else
      useStart=1
    endif

    listLen=size(list)
    valuesLen=size(values)
    if (size(indices) < valuesLen) call MLSMessage(MLSMSG_Error,ModuleName, &
      & "Result array is too small")

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
      if (list(size(list)) >= list(1)) then
        listDirection=1
      else
        listDirection=-1
      endif
    endif

    ! Some last bits of setup before we get going.

    if (useAllowTopValue) then
      upperLimit=listLen
    else
      upperLimit=listLen-1
    endif

    ! Now we're ready to hit the road, loop over all the values to hunt for

    index=max(1,min(useStart,upperLimit))
    do valueIndex=1,valuesLen
      thisValue=values(valueIndex)
      expanding=.true.
      searchDirection=0
      stride=1
      HuntLoop: do
        lowerBelow= (thisValue-list(index))*listDirection >= 0.0D0
        if (index<listLen) then 
          upperAbove= (list(index+1)-thisValue)*listDirection > 0.0D0
        else ! We're off the end of the list
          upperAbove=.true.   
        endif

        ! Now we know what the state of play is, what does it mean?

        ! First see if we've found the place
        if (lowerBelow.and.upperAbove) exit HuntLoop

        ! The other cases are a little more complex
        oldSearchDirection=searchDirection

        if (lowerBelow.and. (.not. upperAbove)) then
          ! If we're at the end, get out
          if (index==upperLimit) exit HuntLoop
          ! We're too low, keep looking upwards
          index=index+stride
          searchDirection=1
        endif

        if ((.not. lowerBelow).and.upperAbove) then
          ! If we're at the begning, get out
          if (index==1) exit HuntLoop

          ! We're too high but not at begining, look back downwards
          index=index-stride
          searchDirection=-1
        endif

        ! Now the very first change of direction is the end of hte
        ! `expanding' phase

        if ( (searchDirection /= oldSearchDirection) .and. &
          & (oldSearchDirection /= 0) .and. (expanding)) &
          & expanding=.false.

        if (expanding) then
          stride=min(stride*2,listLen/2)
        else
          stride=max(stride/2,1)
        endif

        ! Make sure we don't fall off an end

        index=min(max(index,1),upperLimit)
      end do HuntLoop

      ! Final check for off the bottom of the list

      if (useAllowBelowValue) then
        if ((thisValue-list(index))*listDirection<0.0D0) index=0
      endif

      indices(valueIndex)=index
    enddo
  end subroutine HuntArray

  ! ---------------------------------------------------------------------------

  ! This routine is a scalar wrapper for the above one

  subroutine HuntScalar(list,value,index,start,allowTopValue,allowBelowValue)

    ! Dummy arguments
    real(R8), dimension(:), intent(IN) :: list ! List to search
    real(R8), intent(IN) :: value ! Value to search for
    integer, intent(OUT) :: index ! Resulting index
    integer, intent(IN), optional :: start ! Optional start index
    logical, optional, intent(IN) :: allowTopValue ! Can return N
    logical, optional, intent(IN) :: allowBelowValue ! Can return 0

    ! Local variables

    real(R8), dimension(1) :: values ! To pass to HuntArray
    integer, dimension(1) :: indices ! To pass to HuntScalar

    values(1)=value
    call HuntArray(list,values,indices,start,allowTopValue,allowBelowValue)
    index=indices(1)
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

  subroutine InterpolateArray(oldX,oldY,newX,newY,method,extrapolate, &
    & badValue,missingRegions,dyByDx,dNewByDOld)

    ! Dummy arguments
    real(R8), dimension(:), intent(IN) :: oldX
    real(R8), dimension(:,:), intent(IN) :: oldY
    real(R8), dimension(:), intent(IN) :: newX
    real(R8), dimension(:,:), intent(OUT) :: newY

    character (LEN=*), intent(IN) :: method ! See comments above
    character (LEN=*), optional, intent(IN) :: extrapolate ! See comments above
    real(R8), optional, intent(IN) :: badValue
    real(R8), dimension(:,:), optional, intent(OUT) :: dyByDx
    logical, optional, intent(IN) :: missingRegions ! Allow missing regions
    type (MatrixElement_T), intent(OUT), optional :: dNewByDOld ! Derivatives

    ! Local variables
    integer :: noOld,noNew,width ! Dimensions
    logical :: spline           ! Flag
    logical :: useMissingRegions ! Copy of missing regions
    integer :: ind,newInd       ! Loop counters
    logical :: computeDNewByDOld   ! Set if dNewByDOld is present

    integer, dimension(:),    pointer :: lowerInds=> NULL()
    integer, dimension(:),    pointer :: upperInds=> NULL()
    real(R8), dimension(:),   pointer :: maskVector=> NULL()
    real(R8), dimension(:),   pointer :: gap=>NULL()
    real(R8), dimension(:),   pointer :: gap2=> NULL()
    real(R8), dimension(:,:), pointer :: spreadGap=> NULL()
    real(R8), dimension(:,:), pointer :: oldSecond=> NULL()
    real(R8), dimension(:,:), pointer :: oldYupper=> NULL()
    real(R8), dimension(:,:), pointer :: oldYlower=> NULL()
    real(R8), dimension(:,:), pointer :: oldSecondLower=> NULL()
    real(R8), dimension(:,:), pointer :: oldSecondUpper=> NULL()
    real(R8), dimension(:),   pointer :: A=> NULL()
    real(R8), dimension(:),   pointer :: B=> NULL()
    real(R8), dimension(:),   pointer :: C=> NULL()
    real(R8), dimension(:),   pointer :: D=> NULL() ! Coefficients
    real(R8), dimension(:,:), pointer :: AA=> NULL()
    real(R8), dimension(:,:), pointer :: BB=> NULL()
    real(R8), dimension(:,:), pointer :: CC=> NULL()
    real(R8), dimension(:,:), pointer :: DD=> NULL() ! Spread coefs.
    real(R8), dimension(:,:), pointer :: temp=> NULL() ! For 2nd der. guess
    real(R8), dimension(:),   pointer :: p=> NULL() ! For 2nd der. guess
    real(R8), dimension(:,:), pointer :: tempDNewByDOld=> NULL() ! Dense version.
    real(R8) :: sig       ! For second derivative guesser

    character :: extrapolateMethod ! Tidy copy of extrapolate parameter

    ! Executable code

    ! Size the problem, check sanity, set up arrays etc.

    noOld=size(oldX,1)
    noNew=size(newX,1)
    width=size(oldY,2)

    spline=(Capitalize(method(1:1))=="S")

    extrapolateMethod="A"
    if (present(extrapolate)) extrapolateMethod=Capitalize(extrapolate(1:1))

    useMissingRegions=.false.
    if (present(missingRegions)) useMissingRegions=missingRegions

    computeDNewByDOld=present(dNewByDOld)

    if (useMissingRegions.and.spline) call MLSMessage(MLSMSG_Error,ModuleName, &
      & "Cannot use missing regions with spline")

    if (computeDNewByDOld .and. spline) call MLSMessage(MLSMSG_Error,ModuleName,&
      "Cannont get dNewBydOld from spline")

    call Allocate_Test(lowerInds,noNew,"lowerInds",ModuleName)
    call Allocate_Test(upperInds,noNew,"upperInds",ModuleName)
    call Allocate_Test(gap,noNew,"gap",ModuleName)
    call Allocate_Test(A,noNew,"A",ModuleName)
    call Allocate_Test(B,noNew,"B",ModuleName)
    call Allocate_Test(AA,noNew,width,"AA",ModuleName)
    call Allocate_Test(BB,noNew,width,"BB",ModuleName)
    call Allocate_Test(oldYlower,noNew,width,"oldYlower",ModuleName)
    call Allocate_Test(oldYupper,noNew,width,"oldYupper",ModuleName)

    ! Setup arrays needed if dyByDx is requested

    if (present(dyByDx)) call Allocate_Test(spreadGap,noNew,width,&
      "spreadGap",ModuleName)

    ! Setup Matrix block needed if DNewByDOld is needed.
    if (computeDNewByDOld) then
      call CreateBlock_0(dNewByDOld,noNew*width,noOld*width,M_Column_Sparse,&
        NumberNonZero=2*noNew*width)
    endif

    ! Do special stuff for the case of spline, allocate arrays, find 2nd
    ! derivatives etc.

    if (spline) then
      call Allocate_Test(oldSecond,noOld,width,"oldSecond",ModuleName)
      call Allocate_Test(C,noNew,"C",ModuleName)
      call Allocate_Test(D,noNew,"D",ModuleName)
      call Allocate_Test(CC,noNew,width,"CC",ModuleName)
      call Allocate_Test(DD,noNew,width,"DD",ModuleName)
      call Allocate_Test(oldSecondlower,noNew,width,"oldSecondlower",ModuleName)
      call Allocate_Test(oldSecondupper,noNew,width,"oldSecondupper",ModuleName)
      call Allocate_Test(gap2,noNew,"gap2",ModuleName)
      call Allocate_Test(temp,noOld,width,"temp",ModuleName)
      call Allocate_Test(p,width,"p",ModuleName)

      ! Here we have to solve the a tridiagonal equation
      ! This is a straight copy of my idl code
      oldSecond(1,:)=0.0D0
      temp(1,:)=0.0D0
      do ind=2,noOld-1
        sig=(oldX(ind)-oldX(ind-1))/(oldX(ind+1)-oldX(ind-1))
        p=sig*oldSecond(ind-1,:)+2.0D0
        oldSecond(ind,:)=(sig-1.0D0)/p
        temp(ind,:)=(oldY(ind+1,:)-oldY(ind,:))/(oldX(ind+1)-oldX(ind)) - &
          & (oldY(ind,:)-oldY(ind-1,:))/(oldX(ind)-oldX(ind-1))
        temp(ind,:)=(6.0D0*temp(ind,:)/ &
          & (oldX(ind+1)-oldX(ind-1))-sig*temp(ind-1,:))/p
      enddo
      oldSecond(noOld,:)=0.0D0

      ! Now do the back substitution
      do ind=noOld-1,1,-1
        oldSecond(ind,:)=oldSecond(ind,:)*oldSecond(ind+1,:)+temp(ind,:)
      enddo

      call Deallocate_test(temp,"Temp",ModuleName) 
      call Deallocate_test(p,"p",ModuleName)
    endif

    ! Now we're ready to begin the real work.

    ! Clear the result array(s)
    newY=0.0D0
    if (present(dyByDx)) dyByDx=0.0D0

    ! Now hunt for the indices

    call Hunt(oldX,newX,lowerInds)
    upperInds=lowerInds+1
    gap=oldX(upperInds)-oldX(lowerInds)
    if (present(dyByDx)) spreadGap=spread(gap,2,width)

    A=(oldX(upperInds)-newX)/gap

    ! If extrapolate is "C"onstant, deal with that
    if (extrapolateMethod=="C") A=max(min(A,1.0D0),0.0D0)

    B=1.0D0-A

    ! If extrapolate mode is "B"ad, deal with that
    if (extrapolateMethod=="B") then
      call Allocate_Test(maskVector,noNew,"maskVector",ModuleName)
      maskVector=0.0D0
      where ((A<0.0D0).or.(A>1.0D0))
        maskVector=badValue
        A=0.0D0
        B=0.0D0
      end where
      newY=spread(maskVector,2,width)
      if (present(dyByDx)) dyByDx=newY
      call Deallocate_Test(maskVector,"maskVector",ModuleName)
    endif

    ! Now spread out the coefficients
    AA=spread(A,2,width)
    BB=spread(B,2,width)
    oldYlower=oldY(lowerInds,:)
    oldYupper=oldY(upperInds,:)

    ! Now worry about the missing regions flag
    if (useMissingRegions) then
      where( (oldYlower==badValue) .or. (oldYupper==badValue))
        newY=badValue
        AA=0.0D0
        BB=0.0D0
      end where
      if (present(dyByDx)) then
        where( (oldYlower==badValue) .or. &
          & (oldYupper==badValue))
        dyByDx=badValue
        oldYlower=0.0      ! Only way to guarentee bad derivative
        oldYupper=0.0      ! But don't need to worry about spline
      ENDWHERE
    endif
  endif

  ! Now do the linear interpolation calculation
  newY=newY+AA*oldYlower+BB*oldYupper
  if (present(dyByDx)) dyByDx=(oldYupper-oldYlower)/spreadGap

  ! Write the output derivative matrix if needed
  if (computeDNewByDOld) then
    ! While the matrix is ideally suited to row sparse, our storage method
    ! is column sparse, so to be lazy we'll create it full and then sparsify
    ! it.
!
    call Allocate_Test(tempDNewByDOld,noNew*width,noOld*width,&
      "tempDNewByDOld",ModuleName)
    do newInd=1,noNew
      do ind=1,width
        tempDNewByDOld(newInd+ind*noNew,lowerInds(newInd)+ind*noOld)=A(newInd)
        tempDNewByDOld(newInd+ind*noNew,upperInds(newInd)+ind*noOld)=B(newInd)
      enddo
    enddo
    call Sparsify ( tempDNewByDOld, dNewbyDOld, &
      & "tempDNewByDOld", ModuleName ) ! dNewbyDOld := tempDNewByDOld
  endif

  ! Now do the spline calculation
  if (spline) then
    gap2=gap**2
    C=(A**3-A)*gap2/6.0D0    ! Note the extrapolate bad case is covered as..
    D=(B**3-B)*gap2/6.0D0    !   A=B=0.0D0

    ! Spread out the coefficients etc.
    CC=spread(C,2,width)
    DD=spread(D,2,width)
    oldSecondLower=oldSecond(lowerInds,:)
    oldSecondUpper=oldSecond(upperInds,:)

    newY=newY+CC*oldSecondLower+DD*oldSecondUpper
    if (present(dyByDx)) dyByDx=dyByDx+(spreadGap/6.0D0)*( &
      & (3.0D0*BB**2-1.0D0)*oldSecondUpper- &
      & (3.0D0*AA**2-1.0D0)*oldSecondLower)
  endif

  call Deallocate_Test(lowerInds,"lowerInds",ModuleName)
  call Deallocate_Test(upperInds,"upperInds",ModuleName)
  call Deallocate_Test(gap,"gap",ModuleName)
  call Deallocate_Test(A,"A",ModuleName)
  call Deallocate_Test(B,"B",ModuleName)
  call Deallocate_Test(AA,"AA",ModuleName)
  call Deallocate_Test(BB,"BB",ModuleName)
  call Deallocate_Test(oldYlower,"oldYlower",ModuleName)
  call Deallocate_Test(oldYupper,"oldYupper",ModuleName)

  if (spline) then
    call Deallocate_Test(oldSecond,"oldSecond",ModuleName)
    call Deallocate_Test(C,"C",ModuleName)
    call Deallocate_Test(D,"D",ModuleName)
    call Deallocate_Test(CC,"CC",ModuleName)
    call Deallocate_Test(DD,"DD",ModuleName)
    call Deallocate_Test(oldSecondlower,"oldSecondlower",ModuleName)
    call Deallocate_Test(oldSecondupper,"oldSecondupper",ModuleName)
    call Deallocate_Test(gap2,"gap2",ModuleName)
    call Deallocate_Test(temp,"temp",ModuleName)
    call Deallocate_Test(p,"p",ModuleName)
  endif
  if (present(dyByDx)) call Deallocate_Test(spreadGap,"spreadGap",ModuleName)

end subroutine InterpolateArray

! --------------------------------------------------------------------------

! This subroutine is a scalar wrapper for the first one.

subroutine InterpolateScalar(oldX,oldY,newX,newY,method,extrapolate, &
  & badValue,missingRegions,dyByDx)

  ! Dummy arguments
  real(R8), dimension(:), intent(IN) :: oldX
  real(R8), dimension(:), intent(IN) :: oldY
  real(R8), dimension(:), intent(IN) :: newX
  real(R8), dimension(:), intent(OUT) :: newY

  character (LEN=*), intent(IN) :: method ! See comments above
  character (LEN=*), optional, intent(IN) :: extrapolate ! See comments above
  real(R8), optional, intent(IN) :: badValue
  real(R8), dimension(:), optional, intent(OUT) :: dyByDx
  logical, optional, intent(IN) :: missingRegions ! Allow missing regions

  real(R8), dimension(:,:), pointer :: tempResult => NULL()
  real(R8), dimension(:,:), pointer :: tempDerivative => NULL()

  ! Executable code

  call Allocate_Test(tempResult,size(newX),1,"tempResult",ModuleName)
  if (present(dyByDx)) then
    call Allocate_Test(tempDerivative,size(newX),1,"tempDerivative",ModuleName)

    call InterpolateArray(oldX,spread(oldY,2,1),newX,tempResult,method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions, dyByDx=tempDerivative)
    dyByDx=reshape(tempDerivative,shape(newX))

    call Deallocate_Test(tempDerivative,"tempDerivative",ModuleName)
  else
    call InterpolateArray(oldX,spread(oldY,2,1),newX,tempResult,method, &
      & extrapolate=extrapolate, badValue=badValue, &
      & missingRegions=missingRegions)
  endif
  newY=reshape(tempResult,shape(newX))
  call Deallocate_test(tempResult,"tempResult",ModuleName)
end subroutine InterpolateScalar

!=============================================================================
end module MLSNumerics
!=============================================================================

!
! $Log$
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
