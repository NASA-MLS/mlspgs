! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE VerticalCoordinate      ! Defines enum. type for MLS vert. coords.
!=============================================================================

  USE MLSCommon
  USE MLSStrings                ! String handling routines
  USE Units

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id,ModuleName
!------------------------------- RCS Ident Info ------------------------------
CHARACTER(LEN=130) :: id = & 
   "$Id$"
CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!-----------------------------------------------------------------------------

!
! This module simply defines an enumerated type giving the various MLS vertical
! coordinates, and also supplies a routine for parsing requests.
!

INTEGER, PARAMETER :: VC_Invalid    =-1 ! Invalid coordinate system (error)
INTEGER, PARAMETER :: VC_None       =0 ! No vertical coordinate
INTEGER, PARAMETER :: VC_Pressure   =1 ! Pressure in mb
INTEGER, PARAMETER :: VC_Zeta       =2 ! -alog10(pressure/mb)
INTEGER, PARAMETER :: VC_Altitude   =3 ! Altitude in meters
INTEGER, PARAMETER :: VC_GPH        =4 ! Geopotential height in meters
INTEGER, PARAMETER :: VC_Theta      =5 ! Potential Temperature in Kelvin
INTEGER, PARAMETER :: VC_Angle      =6 ! A field of view angle (radians?)

! This routine takes a string and returns a VC_???? value corresponding to
! the name, or VC_Invalid(=-1) if it's an invalid request

CONTAINS

! This routine takes a string and returns a VC_???? value corresponding to
! the name, or VC_Invalid(=-1) if it's an invalid request

FUNCTION ParseVerticalCoordinateName(name)

  ! Dummy arguments

  CHARACTER (LEN=*), INTENT(IN) :: name ! Input name e.g. "pressure"
  INTEGER :: ParseVerticalCoordinateName

  ! Executable code
  SELECT CASE (Capitalize(name))
  CASE ('NONE')
     ParseVerticalCoordinateName=VC_None
  CASE ('PRESSURE') 
     ParseVerticalCoordinateName=VC_Pressure
  CASE ('ZETA') 
     ParseVerticalCoordinateName=VC_Zeta
  CASE ('ALTITUDE') 
     ParseVerticalCoordinateName=VC_Altitude
  CASE ('GPH') 
     ParseVerticalCoordinateName=VC_GPH
  CASE ('THETA')
     ParseVerticalCoordinateName=VC_Theta
  CASE ('ANGLE')
     ParseVerticalCoordinateName=VC_Angle
  CASE DEFAULT
     ParseVerticalCoordinateName=VC_Invalid
  END SELECT
END FUNCTION ParseVerticalCoordinateName

! ----------------------------------------------------------------------------

! This routine takes a string containing a request for vertical coordinate
! values, as for example in the vGrid specification in the L2CF and creates an
! array of said values.

SUBROUTINE ParseVertCoordSpec(line,values,family)

  ! Dummy arguments
  CHARACTER (LEN=*), INTENT(IN) :: line
  REAL (r8), DIMENSION(:), POINTER :: values
  INTEGER, INTENT(OUT) :: family

  ! Local variables
  CHARACTER (LEN=LEN(line)) :: method
  CHARACTER (LEN=LEN(line)) :: thisWord
  CHARACTER (LEN=LEN(line)) :: remainder
  CHARACTER (LEN=LEN(line)) :: newRemainder

  REAL (r8), DIMENSION(:), POINTER :: tmpValues
  REAL (r8) :: startValue,resolution,thisValue
  INTEGER :: noSurfs, noSurfsInSection, surfNo, thisFamily, status

  ! Executable code

  ! Get the first word, this is method, one of log, linear or explicit

  CALL SplitWords(line,method,remainder,delimiter=" ")

  ! Setup some stuff
  family=PHYQ_Invalid
  noSurfs=0
  method=Capitalize(method)
  IF (TRIM(method)=="EXPLICIT") THEN
     ! This is an `explicit' request.
     ExplicitLoop: DO
        ! Get this value
        IF (LEN_TRIM(remainder)==0) EXIT ExplicitLoop
        CALL ReadNumberWithUnitsFromString(remainder,thisValue,thisFamily,&
             & remainder=newRemainder,unitIsOptional=.TRUE.)
        remainder=newRemainder

        ! First make sure the family is consistent with what we have so far
        IF (thisFamily/=PHYQ_Invalid) THEN
           IF ((thisFamily/=family).AND.(family/=PHYQ_Invalid)) &
                & CALL MLSMessage(MLSMSG_Error,ModuleName,&
                & "Inconsistent unit physical quantity families")
           IF (family==PHYQ_Invalid) family=thisFamily
        ENDIF

        ! Now add this value to our list
        noSurfs=noSurfs+1
        IF (noSurfs==1) THEN
           ALLOCATE(values(1),STAT=status)
           IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                & MLSMSG_Allocate//"values")
           values(1)=thisValue
        ELSE
           ALLOCATE(tmpValues(noSurfs),STAT=status)
           IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                & MLSMSG_Allocate//"tmpValues")
           tmpValues(1:noSurfs-1)=values
           tmpValues(noSurfs)=thisValue
           DEALLOCATE(values, STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_DeAllocate//"values")
           values=>tmpValues
        ENDIF
     END DO ExplicitLoop
     IF (family==PHYQ_Invalid) CALL MLSMessage(MLSMSG_Error,ModuleName, &
          & "Must give a unit for at least one explicit coordinate")
  ELSE
     ! We're dealing with Log or Linear type requests
     ! The next word is a starting point
     CALL ReadNumberWithUnitsFromString(remainder,startValue,family,&
          & remainder=newRemainder)
     remainder=newRemainder
     ! Now only allow log for pressure

     IF ((TRIM(method)=="LOG").AND.(family/=PHYQ_Pressure)) CALL MLSMessage( &
          & MLSMSG_Error,ModuleName,"Can only use Log for pressure coordinates")
     IF ((family==PHYQ_Pressure).AND.(Capitalize(TRIM(method))/="LOG")) CALL MLSMessage( &
          & MLSMSG_Error,ModuleName, &
          & "Expecting Log not Linear for pressure coordinates")

     ! Now the remaining sections are requests for numbers of surfaces
     ! at a particular resolution.

     ! Now setup the start value as one surface

     noSurfs=1
     ALLOCATE(values(1),STAT=status)
     IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
          & MLSMSG_Allocate//"values(1)")
     values(1)=startValue

     LogLinearLoop: DO
        CALL SplitWords(remainder,thisWord,newRemainder,delimiter=" ")
        remainder=newRemainder
        CALL ReadNumberWithUnitsFromString(remainder,resolution,thisFamily,&
             & remainder=newRemainder,unitIsOptional=(TRIM(method)=="LOG"))
        remainder=newRemainder
        IF (family==PHYQ_Pressure) THEN
           IF ((thisFamily/=PHYQ_Invalid).AND.&
                &    (thisFamily/=PHYQ_Dimensionless)) &
                & CALL MLSMessage(MLSMSG_Error, ModuleName,&
                &  "Resolution must be dimensionless for pressure coordinates")
        ELSE
           IF (thisFamily /= family) CALL MLSMessage(MLSMSG_Error,ModuleName, &
                & "Inconsistent family between start value and resolution")
        ENDIF


        READ (UNIT=thisWord,FMT=*) noSurfsInSection

        ALLOCATE(tmpValues(noSurfsInSection+noSurfs),STAT=status)
        IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
             & MLSMSG_Allocate//"tmpValues")
        tmpValues(1:noSurfs)=values
        DEALLOCATE(values, STAT=status)
           IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_DeAllocate//"values")
        values=>tmpValues

        ! Now add these surfaces to those already present

        SELECT CASE (TRIM(method))
        CASE ("LINEAR")
           DO surfNo=1,noSurfsInSection
              values(noSurfs+surfNo)=values(noSurfs+surfNo-1)+resolution
           END DO
        CASE ("LOG")
           DO surfNo=1,noSurfsInSection
              values(noSurfs+surfNo)=values(noSurfs+surfNo-1)/ &
                   & 10.0**(1.0/resolution)
           ENDDO
        CASE DEFAULT
           CALL MLSMessage(MLSMSG_Error,ModuleName, &
                & "Log, Linear or Explicit expected")
        END SELECT
        noSurfs=noSurfs+noSurfsInSection

        IF (LEN_TRIM(remainder)==0) EXIT LogLinearLoop
     END DO LogLinearLoop
  ENDIF
  IF (noSurfs==0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
       & "No vertical coordinate values set")

END SUBROUTINE ParseVertCoordSpec
  
  
!=============================================================================
END MODULE VerticalCoordinate
!=============================================================================

! $Log$
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.8  2000/06/23 01:08:48  vsnyder
! Delete unused variables (except ID) to keep NAG f95 happy
!
