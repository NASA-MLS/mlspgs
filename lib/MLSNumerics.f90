!=============================================================================
MODULE MLSNumerics              ! Some low level numerical stuff
!=============================================================================

  UsE MLSCommon
  USE MLSMessageModule
  USE MLSStrings

  IMPLICIT NONE
  
  PUBLIC

  PRIVATE :: Id,ModuleName
  !------------------------------- RCS Ident Info ----------------------------
  CHARACTER(LEN=256) :: Id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

!
! This module contains some low level numerical stuff, hunting, interpolating
! etc.
!

  INTERFACE Hunt
     MODULE PROCEDURE HuntArray
     MODULE PROCEDURE HuntScalar
  END INTERFACE

  INTERFACE InterpolateValues
     MODULE PROCEDURE InterpolateArray
     MODULE PROCEDURE InterpolateScalar
  END INTERFACE

  PRIVATE :: InterpolateArray, InterpolateScalar, HuntArray, HuntScalar

CONTAINS

  ! This routine does the classic hunt the value kind of thing.  This does the
  ! hunt/bisect implemention a la Numerical Recipes.  List must be
  ! monotonically increasing or decreasing. There is no such requirements for
  ! values.

  SUBROUTINE HuntArray(list,values,indices,start,allowTopValue,allowBelowValue)

    ! Dummy arguments
    REAL(R8), DIMENSION(:), INTENT(IN) :: list ! List to search
    REAL(R8), DIMENSION(:), INTENT(IN) :: values ! Values to search for
    INTEGER, DIMENSION(:), INTENT(OUT) :: indices ! Result
    INTEGER, OPTIONAL, INTENT(IN) :: start ! Optional start index
    LOGICAL, OPTIONAL, INTENT(IN) :: allowTopValue ! Can return N
    LOGICAL, OPTIONAL, INTENT(IN) :: allowBelowValue ! Can return 0

    ! Local variables
    INTEGER :: listLen, valuesLen ! Array sizes
    INTEGER :: valueIndex       ! Loop counters
    INTEGER :: index            ! Temporary result

    LOGICAL :: useAllowTopValue, useAllowBelowValue
    INTEGER :: useStart
    INTEGER :: upperLimit       ! Highest value that can be returned
    INTEGER :: stride           ! Value to step by
    LOGICAL :: expanding        ! Whether we're expanding or reducing our search
    LOGICAL :: lowerBelow       ! Flag
    LOGICAL :: upperAbove       ! Another flag

    INTEGER :: listDirection    ! +1 if list ascends, -1 descends
    INTEGER :: searchDirection  ! (in index space) 
    INTEGER :: oldSearchDirection ! Previous value of above

    REAL(R8) :: thisValue

    ! Executable code

    IF (PRESENT(allowTopValue)) THEN
       useAllowTopValue=allowTopValue
    ELSE
       useAllowTopValue=.FALSE.
    ENDIF

    IF (PRESENT(allowBelowValue)) THEN
       useAllowBelowValue=allowBelowValue
    ELSE
       useAllowBelowValue=.FALSE.
    ENDIF

    IF (PRESENT(start)) THEN
       useStart=start
    ELSE
       useStart=1
    ENDIF

    listLen=SIZE(list)
    valuesLen=SIZE(values)
    IF (SIZE(indices) < valuesLen) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Result array is too small")

    ! Note that this next bit makes two assumptions:
    !   1 - list(1)/=list(2)
    !   2 SIZE(list) > 1
    ! We might want to be more clever with this later

    IF (list(2) >= list(1)) THEN
       listDirection=1
    ELSE
       listDirection=-1
    ENDIF

    ! Some last bits of setup before we get going.

    IF (useAllowTopValue) THEN
       upperLimit=listLen
    ELSE
       upperLimit=listLen-1
    ENDIF

    ! Now we're ready to hit the road, loop over all the values to hunt for

    index=MAX(1,MIN(useStart,upperLimit))
    DO valueIndex=1,valuesLen
       thisValue=values(valueIndex)
       expanding=.TRUE.
       searchDirection=0
       stride=1
       HuntLoop: DO
          lowerBelow= (thisValue-list(index))*listDirection >= 0.0D0
          IF (index<listLen) THEN 
             upperAbove= (list(index+1)-thisValue)*listDirection > 0.0D0
          ELSE ! We're off the end of the list
             upperAbove=.TRUE.   
          ENDIF

          ! Now we know what the state of play is, what does it mean?

          ! First see if we've found the place
          IF (lowerBelow.AND.upperAbove) EXIT HuntLoop

          ! The other cases are a little more complex
          oldSearchDirection=searchDirection

          IF (lowerBelow.AND. (.NOT. upperAbove)) THEN
             ! If we're at the end, get out
             IF (index==upperLimit) EXIT HuntLoop
             ! We're too low, keep looking upwards
             index=index+stride
             searchDirection=1
          ENDIF

          IF ((.NOT. lowerBelow).AND.upperAbove) THEN
             ! If we're at the begning, get out
             IF (index==1) EXIT HuntLoop

             ! We're too high but not at begining, look back downwards
             index=index-stride
             searchDirection=-1
          ENDIF

          ! Now the very first change of direction is the end of hte
          ! `expanding' phase

          IF ( (searchDirection /= oldSearchDirection) .AND. &
               & (oldSearchDirection /= 0) .AND. (expanding)) &
               & expanding=.FALSE.

          IF (expanding) THEN
             stride=MIN(stride*2,listLen/2)
          ELSE
             stride=MAX(stride/2,1)
          ENDIF

          ! Make sure we don't fall off an end

          index=MIN(MAX(index,1),upperLimit)
       END DO HuntLoop

       ! Final check for off the bottom of the list

       IF (useAllowBelowValue) THEN
          IF ((thisValue-list(index))*listDirection<0.0D0) index=0
       ENDIF

       indices(valueIndex)=index
    ENDDO
  END SUBROUTINE HuntArray

  ! ---------------------------------------------------------------------------

  ! This routine is a scalar wrapper for the above one

  SUBROUTINE HuntScalar(list,value,index,start,allowTopValue,allowBelowValue)
    
    ! Dummy arguments
    REAL(R8), DIMENSION(:), INTENT(IN) :: list ! List to search
    REAL(R8), INTENT(IN) :: value ! Value to search for
    INTEGER, INTENT(OUT) :: index ! Resulting index
    INTEGER, INTENT(IN), OPTIONAL :: start ! Optional start index
    LOGICAL, OPTIONAL, INTENT(IN) :: allowTopValue ! Can return N
    LOGICAL, OPTIONAL, INTENT(IN) :: allowBelowValue ! Can return 0

    ! Local variables

    REAL(R8), DIMENSION(1) :: values ! To pass to HuntArray
    INTEGER, DIMENSION(1) :: indices ! To pass to HuntScalar

    values(1)=value
    CALL HuntArray(list,values,indices,start,allowTopValue,allowBelowValue)
    index=indices(1)
  END SUBROUTINE HuntScalar

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

  SUBROUTINE InterpolateArray(oldX,oldY,newX,newY,method,extrapolate, &
       & badValue,missingRegions,dyByDx)

    ! Dummy arguments
    REAL(R8), DIMENSION(:), INTENT(IN) :: oldX
    REAL(R8), DIMENSION(:,:), INTENT(IN) :: oldY
    REAL(R8), DIMENSION(:), INTENT(IN) :: newX
    REAL(R8), DIMENSION(:,:), INTENT(OUT) :: newY

    CHARACTER (LEN=*), INTENT(IN) :: method ! See comments above
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: extrapolate ! See comments above
    REAL(R8), OPTIONAL, INTENT(IN) :: badValue
    REAL(R8), DIMENSION(:,:), OPTIONAL, INTENT(OUT) :: dyByDx
    LOGICAL, OPTIONAL, INTENT(IN) :: missingRegions ! Allow missing regions

    ! Local parameters
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_Allocate="Allocation failed for "

    ! Local variables
    INTEGER :: noOld,noNew,width ! Dimensions
    LOGICAL :: spline           ! Flag
    LOGICAL :: useMissingRegions ! Copy of missing regions
    INTEGER :: ind           ! Loop counter
    INTEGER :: status           ! From Allocate

    INTEGER, DIMENSION(:), ALLOCATABLE :: lowerInds,upperInds
    REAL(R8), DIMENSION(:),   ALLOCATABLE :: maskVector
    REAL(R8), DIMENSION(:),   ALLOCATABLE :: gap,gap2
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: spreadGap
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: oldSecond
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: oldYupper,oldYlower
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: oldSecondLower
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: oldSecondUpper
    REAL(R8), DIMENSION(:),   ALLOCATABLE :: A,B,C,D ! Coefficients
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: AA,BB,CC,DD ! Spread coefs.
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: temp ! For 2nd der. guess
    REAL(R8), DIMENSION(:),   ALLOCATABLE :: p ! For 2nd der. guess
    REAL(R8) :: sig       ! For second derivative guesser
    
    CHARACTER :: extrapolateMethod ! Tidy copy of extrapolate parameter

    ! Executable code

    ! Size the problem, check sanity, set up arrays etc.

    noOld=SIZE(oldX,1)
    noNew=SIZE(newX,1)
    width=SIZE(oldY,2)

    spline=(Capitalize(method(1:1))=="S")

    IF (PRESENT(extrapolate)) THEN
       extrapolateMethod=Capitalize(extrapolate(1:1))
    ELSE
       extrapolateMethod="A"
    ENDIF

    IF (PRESENT(missingRegions)) THEN
       useMissingRegions=missingRegions
    ELSE
       useMissingRegions=.FALSE.
    ENDIF

    IF (useMissingRegions.AND.spline) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Cannot use missing regions with spline")

    ALLOCATE(lowerInds(noNew),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lowerInds")
    ALLOCATE(upperInds(noNew),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"upperInds")
    ALLOCATE(gap(noNew),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"gap")

    ALLOCATE(A(noNew),B(noNew),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"A or B")
    ALLOCATE(AA(noNew,width),BB(noNew,width),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"AA or BB")
    ALLOCATE(oldYlower(noNew,width),oldYupper(noNew,width),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"oldYupper or oldYlower")

    ! Setup arrays needed if dyByDx is requested

    IF (PRESENT(dyByDx)) THEN
       ALLOCATE(spreadGap(noNew,width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            MLSMSG_Allocate//"spreadGap")
    ENDIF

    ! Do special stuff for the case of spline, allocate arrays, find 2nd
    ! derivatives etc.

    IF (spline) THEN
       ALLOCATE(oldSecond(noOld,width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            MLSMSG_Allocate//"oldSecond")
       ALLOCATE(C(noNew),D(noNew),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"C or D")
       ALLOCATE(CC(noNew,width),DD(noNew,width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"CC or DD")
       ALLOCATE(oldSecondLower(noNew,width),&
            & oldSecondUpper(noNew,width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"oldSecondLower or oldSecondUpper")
       ALLOCATE(gap2(noNew),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"gap2")

       ALLOCATE(temp(noOld,width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"temp")
       ALLOCATE(p(width),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"p")

       ! Here we have to solve the a tridiagonal equation
       ! This is a straight copy of my idl code
       oldSecond(1,:)=0.0D0
       temp(1,:)=0.0D0
       DO ind=2,noOld-1
          sig=(oldX(ind)-oldX(ind-1))/(oldX(ind+1)-oldX(ind-1))
          p=sig*oldSecond(ind-1,:)+2.0D0
          oldSecond(ind,:)=(sig-1.0D0)/p
          temp(ind,:)=(oldY(ind+1,:)-oldY(ind,:))/(oldX(ind+1)-oldX(ind)) - &
               & (oldY(ind,:)-oldY(ind-1,:))/(oldX(ind)-oldX(ind-1))
          temp(ind,:)=(6.0D0*temp(ind,:)/ &
               & (oldX(ind+1)-oldX(ind-1))-sig*temp(ind-1,:))/p
       ENDDO
       oldSecond(noOld,:)=0.0D0
       
       ! Now do the back substitution
       DO ind=noOld,1,-1
          oldSecond(ind,:)=oldSecond(ind,:)*oldSecond(ind+1,:)+temp(ind,:)
       ENDDO

       DEALLOCATE(temp,p, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"temp p")
    ENDIF

    ! Now we're ready to begin the real work.

    ! Clear the result array(s)
    newY=0.0D0
    IF (PRESENT(dyByDx)) dyByDx=0.0D0

    ! Now hunt for the indices

    CALL Hunt(oldX,newX,lowerInds)
    upperInds=lowerInds+1
    gap=oldX(upperInds)-oldX(lowerInds)
    IF (PRESENT(dyByDx)) spreadGap=SPREAD(gap,2,width)

    A=(oldX(upperInds)-newX)/gap

    ! If extrapolate is "C"onstant, deal with that
    IF (extrapolateMethod=="C") A=MAX(MIN(A,1.0D0),0.0D0)

    B=1.0D0-A

    ! If extrapolate mode is "B"ad, deal with that
    IF (extrapolateMethod=="B") THEN
       ALLOCATE(maskVector(noNew),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"maskVector")
       maskVector=0.0D0
       WHERE ((A<0.0D0).OR.(A>1.0D0))
          maskVector=badValue
          A=0.0D0
          B=0.0D0
       END WHERE
       newY=SPREAD(maskVector,2,width)
       IF (PRESENT(dyByDx)) dyByDx=newY
       DEALLOCATE(maskVector, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"maskVector")
    ENDIF

    ! Now spread out the coefficients
    AA=SPREAD(A,2,width)
    BB=SPREAD(B,2,width)
    oldYlower=oldY(lowerInds,:)
    oldYupper=oldY(upperInds,:)

    ! Now worry about the missing regions flag
    IF (useMissingRegions) THEN
       WHERE( (oldYlower==badValue) .OR. (oldYupper==badValue))
          newY=badValue
          AA=0.0D0
          BB=0.0D0
       END WHERE
       IF (PRESENT(dyByDx)) THEN
          WHERE( (oldYlower==badValue) .OR. &
               & (oldYupper==badValue))
             dyByDx=badValue
             oldYlower=0.0      ! Only way to guarentee bad derivative
             oldYupper=0.0      ! But don't need to worry about spline
          ENDWHERE
       ENDIF
    ENDIF

    ! Now do the linear interpolation calculation
    newY=newY+AA*oldYlower+BB*oldYupper
    IF (PRESENT(dyByDx)) dyByDx=(oldYupper-oldYlower)/spreadGap

    ! Now do the spline calculation
    IF (spline) THEN
       gap2=gap**2
       C=(A**3-A)*gap2/6.0D0    ! Note the extrapolate bad case is covered as..
       D=(B**3-B)*gap2/6.0D0    !   A=B=0.0D0

       ! Spread out the coefficients etc.
       CC=SPREAD(C,2,width)
       DD=SPREAD(D,2,width)
       oldSecondLower=oldSecond(lowerInds,:)
       oldSecondUpper=oldSecond(upperInds,:)
       
       newY=newY+CC*oldSecondLower+DD*oldSecondUpper
       IF (PRESENT(dyByDx)) dyByDx=dyByDx+(spreadGap/6.0D0)*( &
            & (3.0D0*BB**2-1.0D0)*oldSecondUpper- &
            & (3.0D0*AA**2-1.0D0)*oldSecondLower)
    ENDIF

    DEALLOCATE (lowerInds,upperInds,gap,A,B,AA,BB,oldYlower,oldYupper, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"lowerInds")
    IF (spline) DEALLOCATE(gap2,oldSecond,C,D,CC,DD, &
         & oldSecondLower,oldSecondUpper, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"gap2")
    IF (PRESENT(dyByDx)) DEALLOCATE(spreadGap, stat=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"spreadGap")

  END SUBROUTINE InterpolateArray

  ! --------------------------------------------------------------------------

  ! This subroutine is a scalar wrapper for the first one.

  SUBROUTINE InterpolateScalar(oldX,oldY,newX,newY,method,extrapolate, &
       & badValue,missingRegions,dyByDx)

    ! Dummy arguments
    REAL(R8), DIMENSION(:), INTENT(IN) :: oldX
    REAL(R8), DIMENSION(:), INTENT(IN) :: oldY
    REAL(R8), DIMENSION(:), INTENT(IN) :: newX
    REAL(R8), DIMENSION(:), INTENT(OUT) :: newY

    CHARACTER (LEN=*), INTENT(IN) :: method ! See comments above
    CHARACTER (LEN=*), OPTIONAL, INTENT(IN) :: extrapolate ! See comments above
    REAL(R8), OPTIONAL, INTENT(IN) :: badValue
    REAL(R8), DIMENSION(:), OPTIONAL, INTENT(OUT) :: dyByDx
    LOGICAL, OPTIONAL, INTENT(IN) :: missingRegions ! Allow missing regions

    ! Local parameters
    CHARACTER (LEN=*), PARAMETER :: MLSMSG_Allocate="Allocation failed for "

    ! Local variables
    INTEGER :: status

    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: tempResult
    REAL(R8), DIMENSION(:,:), ALLOCATABLE :: tempDerivative

    ! Executable code

    ALLOCATE(tempResult(SIZE(newX),1),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         MLSMSG_Allocate//"tempResult")

    IF (PRESENT(dyByDx)) THEN
       ALLOCATE(tempDerivative(SIZE(newX),1),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            MLSMSG_Allocate//"tempDerivative")

       CALL InterpolateArray(oldX,SPREAD(oldY,2,1),newX,tempResult,method, &
            & extrapolate=extrapolate, badValue=badValue, &
            & missingRegions=missingRegions, dyByDx=tempDerivative)
       dyByDx=RESHAPE(tempDerivative,SHAPE(newX))
       DEALLOCATE(tempDerivative, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"tempDerivative")
    ELSE
       CALL InterpolateArray(oldX,SPREAD(oldY,2,1),newX,tempResult,method, &
            & extrapolate=extrapolate, badValue=badValue, &
            & missingRegions=missingRegions)
    ENDIF
    newY=RESHAPE(tempResult,SHAPE(newX))
    DEALLOCATE(tempResult, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"tempResult")
  END SUBROUTINE InterpolateScalar
    
!=============================================================================
END MODULE MLSNumerics
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.13  2000/06/23 01:08:48  vsnyder
! Delete unused variables (except ID) to keep NAG f95 happy
!
