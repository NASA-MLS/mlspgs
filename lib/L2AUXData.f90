! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L2AUXData                 ! Data types for storing L2AUX data internally
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
  ! manipulating L2AUX data.

  INTEGER, PARAMETER :: L2AUXNameLen=80
  INTEGER, PARAMETER, PRIVATE :: CCSDSLen=27

  ! This is a set of possible values for dimension%dimensionFamily

  INTEGER, PARAMETER :: NoL2AUXDimTypes=6
  CHARACTER(LEN=12), DIMENSION(NoL2AUXDimTypes),PARAMETER::L2AUXDimNames= (/ &
       & "Channel  ", &
       & "Frequency", &
       & "MIF      ", &
       & "MAF      ", &
       & "Time     ", &
       & "GeodAngle" /)

  INTEGER, PARAMETER :: L2AUXDim_None=0
  INTEGER, PARAMETER :: L2AUXDim_Channel=1
  INTEGER, PARAMETER :: L2AUXDim_Frequency=2
  INTEGER, PARAMETER :: L2AUXDim_MIF=3
  INTEGER, PARAMETER :: L2AUXDim_MAF=4
  INTEGER, PARAMETER :: L2AUXDim_Time=5
  INTEGER, PARAMETER :: L2AUXDim_GeodAngle=6

  ! This datatype describes a dimension for an L2AUX quantity

  TYPE L2AUX_Dimension_T
     INTEGER :: noValues        ! Length of this dimension
     INTEGER :: dimensionFamily ! What is this dimension
     REAL(r8), DIMENSION(:), POINTER :: values ! (noValues)
  END TYPE L2AUX_Dimension_T

  ! This datatype describes an l2aux quantity itself.
  ! The dimensions will typically be ordered as follows:
  ! [Channel or frequency], MIF, [MAF or time or geodAngle]

  TYPE L2AUXData_T

    ! A name for the L2AUX quantity, goes into SD name
    CHARACTER (LEN=L2AUXNameLen) :: Name ! Name for quantity to be output

    ! The dimensions for the quantity
    TYPE (L2AUX_Dimension_T), DIMENSION(3) :: dimensions

    REAL(r8), POINTER, DIMENSION(:,:,:) :: values
  END TYPE L2AUXData_T

CONTAINS

  !---------------------------------------------------------------------------

  ! This first routine sets up the arrays for an l2aux datatype.
  ! The user supplies a set of three dimensionFamilies (e.g. L2AUXDim_MAF)
  ! Quantities can have upto three valid dimensions.  L2AUXDim_None can be used
  ! to indicate later dimensions are invalid.

  SUBROUTINE SetupNewL2AUXRecord(dimensionFamilies,dimSizes,l2aux)

    ! Dummy arguments
    INTEGER, DIMENSION(3), INTENT(IN) :: dimensionFamilies
    INTEGER, DIMENSION(3), INTENT(IN) :: dimSizes
    TYPE (L2AUXData_T), INTENT(OUT) :: l2aux

    ! Local variables
    INTEGER :: dimIndex
    INTEGER :: status

    ! Fill the dimensions data structure

    l2aux%dimensions%dimensionFamily=dimensionFamilies
    l2aux%dimensions%noValues=dimSizes

    ! Allocate the values for each dimension

    DO dimIndex=1,3
       IF (dimensionFamilies(dimIndex)/=L2AUXDim_None) THEN
          ALLOCATE(l2aux%dimensions(dimIndex)%values(dimSizes(dimIndex)), &
               & STAT=status)
          IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & MLSMSG_Allocate//"l2aux dimension values")
       ELSE
          l2aux%dimensions(dimIndex)%noValues=1
       END IF
    END DO

    ! Allocate the values for the data itself

    ALLOCATE(l2aux%values( &
         & l2aux%dimensions(1)%noValues, &
         & l2aux%dimensions(2)%noValues, &
         & l2aux%dimensions(3)%noValues), STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate// &
         & "l2aux values")

  END SUBROUTINE SetupNewL2AUXRecord
    
  !---------------------------------------------------------------------------

  ! This routine deallocates all the arrays allocated above.

  SUBROUTINE DestroyL2AUXContents(l2aux)

    ! Dummy arguments
    TYPE (L2AUXData_T), INTENT(INOUT) :: l2aux

    ! Executable code
    DEALLOCATE(l2aux%dimensions(1)%values, &
         & l2aux%dimensions(2)%values, &
         & l2aux%dimensions(3)%values)
    DEALLOCATE(l2aux%values)

  END SUBROUTINE DestroyL2AUXContents

  !---------------------------------------------------------------------------

  ! This subroutine expands an L2AUXData_T in place, allowing the user to
  ! add more `profiles' to it.  Note that the `profile' dimension is the last
  ! one.

  SUBROUTINE ExpandL2AUXDataInPlace(l2aux,newSize)

    ! Dummy arguments
    TYPE (L2AUXData_T), INTENT(INOUT) :: l2aux
    INTEGER, INTENT(IN) :: newSize

    ! Local variables
    INTEGER :: status           ! From ALLOCATE
    ! The following are temporary arrays for copying data around
    REAL (r8), DIMENSION(:), POINTER :: temp1D
    REAL (r8), DIMENSION(:,:,:), POINTER :: temp3D
    INTEGER :: expandingDimension
    INTEGER :: oldSize

    ! Executable code

    ! First identity which is the `last' dimension.

    expandingDimension=3
    DO WHILE( (l2Aux%dimensions(expandingDimension)%dimensionFamily==&
         & L2AUXDim_None).AND.(expandingDimension>1))
       expandingDimension=expandingDimension-1
    END DO

    ! Now see how long this is
    oldSize=l2aux%dimensions(expandingDimension)%noValues

    ! Do a sanity check
    IF (newSize<oldSize) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "This l2aux is getting smaller not bigger")

    ! Now expand this dimension
    temp1D=>l2aux%dimensions(expandingDimension)%values
    ALLOCATE(l2aux%dimensions(expandingDimension)%values(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "New dimension information")
    l2aux%dimensions(expandingDimension)%values(1:oldSize)=temp1D
    DEALLOCATE(temp1D)

    ! Now expand the data in this dimension
    temp3d=>l2aux%values
    SELECT CASE (expandingDimension)
    CASE (1)
       ALLOCATE(l2aux%values(newSize,1,1),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "New l2aux values")
       l2aux%values(1:oldSize,:,:)=temp3d
    CASE(2)
       ALLOCATE(l2aux%values(l2aux%dimensions(1)%noValues,newSize,1),&
            & STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "New l2aux values")
       l2aux%values(:,1:oldSize,:)=temp3d
    CASE(3)
       ALLOCATE(l2aux%values(l2aux%dimensions(1)%noValues, &
            &                l2aux%dimensions(2)%noValues,&
            &                newSize),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & "New l2aux values")
       l2aux%values(:,:,1:oldSize)=temp3d
    END SELECT
    DEALLOCATE(temp3d)

    ! That's it.

  END SUBROUTINE ExpandL2AUXDataInPlace

  !---------------------------------------------------------------------------

  ! This subroutine adds an l2aux data type to a database of said types,
  ! creating a new database if it doesn't exist

  SUBROUTINE AddL2AUXToDatabase(database,l2aux)

    ! Dummy arguments
    TYPE (L2AUXData_T), DIMENSION(:), POINTER :: database
    TYPE (L2AUXData_T), INTENT(IN) :: l2aux

    ! Local variables
    TYPE (L2AUXData_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,l2aux%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//l2aux%name)
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF

    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=l2aux
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddL2AUXToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyL2AUXDatabase(database)

    ! Dummy argument
    TYPE (L2AUXData_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: l2auxIndex

    IF (ASSOCIATED(database)) THEN
       DO l2auxIndex=1,SIZE(database)
          CALL DestroyL2AUXContents(database(l2auxIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyL2AUXDatabase


!=============================================================================
END MODULE L2AUXData
!=============================================================================

!
! $Log$
! Revision 1.8  2000/01/19 21:42:18  livesey
! Just tided up some comments.
!
! Revision 1.7  2000/01/07 23:53:34  livesey
! Nearly integrated, just a few tweaks.
!
! Revision 1.6  1999/12/18 01:06:28  livesey
! Added USE of MLSStrings
!
! Revision 1.5  1999/12/17 21:41:00  livesey
! Added check for duplicate name
!
! Revision 1.4  1999/12/14 00:53:17  livesey
! Changed DOUBLE PRECISION to REAL(r8)
!
! Revision 1.3  1999/12/03 22:25:57  livesey
! Tidied up some of the INTENT stuff
!
! Revision 1.2  1999/12/03 21:22:23  livesey
! Removed old log data
!
! Revision 1.1  1999/12/03 21:22:02  livesey
! First versions, modified from L2GPData module
!
