! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L2GPData                 ! Data types for storing L2GP data internally
!=============================================================================

  USE MLSMessageModule
  USE MLSCommon
  USE MLSStrings

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.

  INTEGER, PARAMETER :: L2GPNameLen=80
  INTEGER, PARAMETER, PRIVATE :: CCSDSLen=27

  ! This datatype is the main one, it simply defines one l2gp swath

  TYPE L2GPData_T

    CHARACTER (LEN=L2GPNameLen) :: Name ! Name for quantity to be output

    INTEGER :: noProfs           ! Total number of profiles
    INTEGER :: noSurfs           ! Total number of surfaces
    INTEGER :: noFreqs           ! Number of frequencies in breakdown

    ! Now we store the geolocation fields, first the vertical one:
    REAL, POINTER, DIMENSION(:) :: pressures ! Vertical coordinates (noSurfs)

    ! Now the horizontal geolocation information. Dimensioned (noProfs)
    REAL, POINTER, DIMENSION(:) :: latitude,longitude,solarTime, &
         & solarZenith, losAngle, geodAngle
    REAL(r8), POINTER, DIMENSION(:) :: time
    INTEGER, POINTER, DIMENSION(:) :: chunkNumber
    CHARACTER (LEN=CCSDSLen), POINTER, DIMENSION(:) :: ccsdsTime

    ! Now we store the `frequency' geolocation field

    REAL(r8), POINTER, DIMENSION(:) :: frequency
    !        dimensioned (noFreqs)

    ! Finally we store the data fields

    REAL, POINTER, DIMENSION(:,:,:) :: l2gpValue
    REAL, POINTER, DIMENSION(:,:,:) :: l2gpPrecision
    ! dimensioned (noFreqs, noSurfs, noProfs)

    CHARACTER (LEN=1), POINTER, DIMENSION(:) :: l2gpStatus
    !                (status is a reserved word in F90)
    REAL, POINTER, DIMENSION(:) :: quality
    ! Both the above dimensioned (noProfs)

    INTEGER :: accumulatedProfiles
    ! This last one is needed by Join to keep track of which profiles
    ! have been output so far.

  END TYPE L2GPData_T

CONTAINS

  !---------------------------------------------------------------------------

  ! This first routine sets up the arrays for an l2gp datatype.

  SUBROUTINE SetupNewL2GPRecord(noProfs,noSurfs,noFreqs,l2gp)

    ! Dummy arguments
    INTEGER, INTENT(IN) :: noProfs,noSurfs,noFreqs ! Dimensions
    TYPE (L2GPData_T), INTENT(OUT)  :: l2gp

    ! Local variables
    INTEGER :: status

    ALLOCATE (l2gp%pressures(noSurfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"l2gp%pressure")

    ALLOCATE (l2gp%latitude(noProfs), l2gp%longitude(noProfs), &
         & l2gp%solarTime(noProfs), l2gp%solarZenith(noProfs), &
         & l2gp%losAngle(noProfs), l2gp%geodAngle(noProfs), &
         & l2gp%chunkNumber(noProfs), l2gp%ccsdsTime(noProfs), &
         & STAT=status)
    IF (status /=0) CALL MLSMessage(MLSMSG_error,ModuleName, &
         & MLSMSG_Allocate//"l2gp horizontal coordinates")

    IF (noFreqs/=0) THEN
       ALLOCATE (l2gp%frequency(noFreqs),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//"l2gp%frequency")
    END IF

    ALLOCATE (l2gp%l2gpValue(MAX(1,noFreqs),noSurfs,noProfs), &
         & l2gp%l2gpPrecision(MAX(1,noFreqs),noSurfs,noProfs), STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"l2gp data/precision fields")

    ALLOCATE (l2gp%l2gpStatus(noProfs),l2gp%quality(noProfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"l2gp status/quality fields")

    l2gp%noProfs=noProfs
    l2gp%noSurfs=noSurfs
    l2gp%noFreqs=noFreqs
  END SUBROUTINE SetupNewL2GPRecord
    
  !---------------------------------------------------------------------------

  ! This routine deallocates all the arrays allocated above.

  SUBROUTINE DestroyL2GPContents(l2gp)

    ! Dummy arguments
    TYPE (L2GPData_T), INTENT(INOUT) :: l2gp

    ! Executable code

    DEALLOCATE(l2gp%pressures)
    DEALLOCATE(l2gp%latitude, l2gp%longitude, l2gp%solarTime, &
         & l2gp%solarZenith, l2gp%losAngle, l2gp%geodAngle, &
         & l2gp%chunkNumber, l2gp%ccsdsTime)
    DEALLOCATE(l2gp%frequency)
    DEALLOCATE(l2gp%l2gpValue,l2gp%l2gpPrecision)
    DEALLOCATE(l2gp%l2gpStatus,l2gp%quality)
    l2gp%noProfs=0
    l2gp%noSurfs=0
    l2gp%noFreqs=0
  END SUBROUTINE DestroyL2GPContents

  !---------------------------------------------------------------------------

  ! This subroutine adds an l2gp data type to a database of said types,
  ! creating a new database if it doesn't exist

  SUBROUTINE AddL2GPToDatabase(database,l2gp)

    ! Dummy arguments
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: database
    TYPE (L2GPData_T), INTENT(IN) :: l2gp

    ! Local variables
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,l2gp%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//l2gp%name)
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF

    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=l2gp
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddL2GPToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyL2GPDatabase(database)

    ! Dummy argument
    TYPE (L2GPData_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: l2gpIndex

    IF (ASSOCIATED(database)) THEN
       DO l2gpIndex=1,SIZE(database)
          CALL DestroyL2GPContents(database(l2gpIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyL2GPDatabase


!=============================================================================
END MODULE L2GPData
!=============================================================================

!
! $Log$
! Revision 1.7  1999/12/21 00:15:15  livesey
! Added accumulatedProfiles to L2GPData_T for Join
!
! Revision 1.6  1999/12/18 01:06:28  livesey
! Added USE of MLSStrings
!
! Revision 1.5  1999/12/17 21:41:31  livesey
! Added check for duplicate name
!
! Revision 1.4  1999/12/14 00:53:50  livesey
! Changed DOUBLE PRECISION to REAL(r8)
!
! Revision 1.3  1999/12/03 22:24:50  livesey
! Tidied up some of the INTENT stuff
!
! Revision 1.2  1999/12/03 21:28:56  livesey
! Renamed L2GP_T to L2GPData_T
!
! Revision 1.1  1999/12/03 19:10:34  livesey
! First version
!
!
