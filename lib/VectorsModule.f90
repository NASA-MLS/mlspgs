! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE VectorsModule            ! Vectors in the MLS PGS suite
!=============================================================================

  USE MLSMessageModule
  USE VectorQuantities

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module provides the simple functionality for vector quantities in the
  ! MLS Level 2 software, and related programs.

  !---------------------------------------------------------------------------

  ! This datatype describes a vector template

  INTEGER, PARAMETER :: VectorNameLen=80

  TYPE VectorTemplate_T
     
     ! First some administrative stuff
     CHARACTER (LEN=VectorNameLen) :: name ! Name for the vector template
     INTEGER :: id              ! Id code for vector (for checking purposes)

     ! Now general information about the vector

     INTEGER :: noQuantities    ! Number of quantities in the vector
     INTEGER :: noSubVectors    ! Number of subvectors in the vector
     INTEGER :: totalElements   ! Total of numbers of elements in the vector

     ! Now we describe the ordering of the subvectors in the vector

     INTEGER, DIMENSION(:), POINTER :: subVectorNoElements ! (noSubVectors)
     INTEGER, DIMENSION(:), POINTER :: subVectorFirstElements ! (noSubVectors)
     INTEGER, DIMENSION(:), POINTER :: subVectorQuantityNo ! (noSubVectors)
     INTEGER, DIMENSION(:), POINTER :: subVectorProfileNo ! (noSubVectors)

     ! Now point to the quantity templates themselves

     TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: quantities
  END TYPE VectorTemplate_T

  ! This datatype is a vector itself

  TYPE Vector_T
     CHARACTER (LEN=VectorNameLen) :: name ! Name for vector
     TYPE (VectorTemplate_T) :: template ! Copy from template database
     REAL(r8), DIMENSION(:), POINTER :: values
  END TYPE Vector_T

  ! This incrementing counter is used to set the id field for a vector template

  INTEGER, PRIVATE :: vectorTemplateCounter=0

CONTAINS

  !---------------------------------------------------------------------------

  ! This subroutine creates a vectorTemplate from a list of quantities.

  SUBROUTINE ConstructVectorTemplate(name,quantities,vectorTemplate)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: name ! Name for vector template
    TYPE (QuantityTemplate_T), DIMENSION(:), INTENT(IN) :: quantities
    TYPE (VectorTemplate_T), INTENT(OUT) :: vectorTemplate

    ! Local variables
    INTEGER :: status

    ! Executable code

    vectorTemplate%name=name
    vectorTemplate%noQuantities=SIZE(quantities)
    vectorTemplate%noSubVectors=SUM(quantities%noSubVectors)
    vectorTemplate%totalElements= &
         & SUM(quantities%noSubVectors*quantities%subVectorLen)
    
    ! Allocate some arrays

    ALLOCATE(&
         & vectorTemplate%subVectorNoElements(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorFirstElements(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorQuantityNo(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorProfileNo(vectorTemplate%noSubVectors), &
         & STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Vector information")

    ALLOCATE (vectorTemplate%quantities(vectorTemplate%noQuantities), &
         & STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Vector quantities")

    vectorTemplate%quantities=quantities

    vectorTemplate%id=vectorTemplateCounter
    vectorTemplateCounter=vectorTemplateCounter+1
    
  END SUBROUTINE ConstructVectorTemplate

  !---------------------------------------------------------------------------

  ! This routine destroys a vector template created above

  SUBROUTINE DestroyVectorTemplateInfo(vectorTemplate)

    ! Dummy arguments
    TYPE (VectorTemplate_T), INTENT(INOUT) :: vectorTemplate

    ! Executable code

    DEALLOCATE (&
         & vectorTemplate%subVectorNoElements, &
         & vectorTemplate%subVectorFirstElements, &
         & vectorTemplate%subVectorQuantityNo, &
         & vectorTemplate%subVectorProfileNo)
    DEALLOCATE(vectorTemplate%quantities)
  END SUBROUTINE DestroyVectorTemplateInfo

  !---------------------------------------------------------------------------

  ! This routine adds a vector template to a database of such templates, 
  ! creating the database if necessary.

  SUBROUTINE AddVectorTemplateToDatabase(database,vectorTemplate)

    ! Dummy arguments
    TYPE (VectorTemplate_T), DIMENSION(:), POINTER :: database
    TYPE (VectorTemplate_T), INTENT(IN) :: vectorTemplate

    ! Local variables
    TYPE (VectorTemplate_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF

    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=vectorTemplate
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddVectorTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a vector template database

  SUBROUTINE DestroyVectorTemplateDatabase(database)

    ! Dummy argument
    TYPE (VectorTemplate_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: l2gpIndex

    IF (ASSOCIATED(database)) THEN
       DO l2gpIndex=1,SIZE(database)
          CALL DestroyVectorTemplateInfo(database(l2gpIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyVectorTemplateDatabase

  ! --------------------------------------------------------------------------

  ! Now, having dealt with the templates, let's move onto the vectors
  ! This routine creates an empty vector according to a given template

  SUBROUTINE CreateVector(name,vectorTemplate,vector)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: name   ! Name for vector
    TYPE (VectorTemplate_T), INTENT(IN) :: vectorTemplate ! Template for vector
    TYPE (Vector_T), INTENT(OUT) :: vector

    ! Local variables
    INTEGER :: status

    ! Executable code

    vector%name=name
    vector%template=vectorTemplate
    ALLOCATE(vector%values(vectorTemplate%totalElements),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"Vector values")
  END SUBROUTINE CreateVector

  ! --------------------------------------------------------------------------
  
  ! This routine destroys a vector created above

  SUBROUTINE DestroyVectorInfo(vector)

    ! Dummy arguments
    TYPE (Vector_T), INTENT(IN) :: vector

    ! Executable code

    DEALLOCATE(vector%values)
  END SUBROUTINE DestroyVectorInfo

  !---------------------------------------------------------------------------

  ! This routine adds a vector template to a database of such templates, 
  ! creating the database if necessary.

  SUBROUTINE AddVectorToDatabase(database,vector)

    ! Dummy arguments
    TYPE (Vector_T), DIMENSION(:), POINTER :: database
    TYPE (Vector_T), INTENT(IN) :: vector

    ! Local variables
    TYPE (Vector_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF

    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=vector
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddVectorToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a vector database

  SUBROUTINE DestroyVectorDatabase(database)

    ! Dummy argument
    TYPE (Vector_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: l2gpIndex

    IF (ASSOCIATED(database)) THEN
       DO l2gpIndex=1,SIZE(database)
          CALL DestroyVectorInfo(database(l2gpIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyVectorDatabase

!=============================================================================
END MODULE VectorsModule
!=============================================================================

!
! $Log$
! Revision 1.1  1999/12/04 00:30:18  livesey
! First version.
!
!
