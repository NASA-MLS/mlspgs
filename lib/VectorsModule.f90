! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE VectorsModule            ! Vectors in the MLS PGS suite
!=============================================================================

  USE MLSMessageModule
  USE QuantityTemplates

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

  ! This datatype describes how the subVectors of a quantity are laid out
  ! within a vector

  TYPE QuantityLayout_T
     INTEGER, DIMENSION(:), POINTER :: entry ! (quantityTemplate%noSubVectors)
  END TYPE QuantityLayout_T

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
     INTEGER, DIMENSION(:), POINTER :: subVectorFirstElement ! (noSubVectors)
     INTEGER, DIMENSION(:), POINTER :: subVectorQuantityNo ! (noSubVectors)
     INTEGER, DIMENSION(:), POINTER :: subVectorProfileNo ! (noSubVectors)

     ! Now a `reverse' form of this information for each quantity

     TYPE (QuantityLayout_T), DIMENSION(:), POINTER :: layout ! (noQuantities)

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
  ! The default ordering is currently by quantity. Later versions may
  ! have optional parameters to request other orderings.

  SUBROUTINE ConstructVectorTemplate(name,quantities,vectorTemplate)

    ! Dummy arguments
    CHARACTER (LEN=*), INTENT(IN) :: name ! Name for vector template
    TYPE (QuantityTemplate_T), DIMENSION(:), INTENT(IN) :: quantities
    TYPE (VectorTemplate_T), INTENT(OUT) :: vectorTemplate

    ! Local variables
    INTEGER :: status,qty,subVector,startSubVector
    INTEGER :: accumulatedElements,qtySubVector,subVectorLen

    ! Executable code

    vectorTemplate%name=name
    vectorTemplate%noQuantities=SIZE(quantities)
    vectorTemplate%noSubVectors=SUM(quantities%noSubVectors)
    vectorTemplate%totalElements= &
         & SUM(quantities%noSubVectors*quantities%subVectorLen)
    
    ! Allocate some arrays

    ALLOCATE(&
         & vectorTemplate%subVectorNoElements(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorFirstElement(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorQuantityNo(vectorTemplate%noSubVectors), &
         & vectorTemplate%subVectorProfileNo(vectorTemplate%noSubVectors), &
         & vectorTemplate%layout(vectorTemplate%noQuantities), &
         & STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"Vector information")

    DO qty=1,vectorTemplate%noQuantities
       ALLOCATE(vectorTemplate%layout(qty)%entry(quantities(qty)%noSubVectors),&
            & STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//"Quantity layout information")
    END DO

    ALLOCATE (vectorTemplate%quantities(vectorTemplate%noQuantities), &
         & STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"Vector quantities")

    ! Copy quantities over

    vectorTemplate%quantities=quantities

    vectorTemplate%id=vectorTemplateCounter
    vectorTemplateCounter=vectorTemplateCounter+1

    ! Fill up the mapping information, store by quantity for this version

    startSubVector=0
    accumulatedElements=0
    DO qty=1,vectorTemplate%noQuantities
       subVectorLen=vectorTemplate%quantities(qty)%subVectorLen
       DO qtySubVector=1,vectorTemplate%quantities(qty)%noSubVectors-1
          subVector=qtySubVector+startSubVector

          vectorTemplate%layout(qty)%entry(qtySubVector)=subVector

          vectorTemplate%subVectorNoElements(subVector)=subVectorLen
          vectorTemplate%subVectorFirstElement(subVector)=accumulatedElements
          vectorTemplate%subVectorQuantityNo(subVector)=qty
          vectorTemplate%subVectorProfileNo(subVector)=qtySubVector

          accumulatedElements=accumulatedElements+subVectorLen
       END DO
    END DO
    
  END SUBROUTINE ConstructVectorTemplate

  !---------------------------------------------------------------------------

  ! This routine destroys a vector template created above

  SUBROUTINE DestroyVectorTemplateInfo(vectorTemplate)

    ! Dummy arguments
    TYPE (VectorTemplate_T), INTENT(INOUT) :: vectorTemplate

    ! Local variables
    INTEGER :: qty

    ! Executable code

    DEALLOCATE (&
         & vectorTemplate%subVectorNoElements, &
         & vectorTemplate%subVectorFirstElement, &
         & vectorTemplate%subVectorQuantityNo, &
         & vectorTemplate%subVectorProfileNo)

    DO qty=1,vectorTemplate%noQuantities
       DEALLOCATE(vectorTemplate%layout(qty)%entry)
    END DO
    DEALLOCATE(vectorTemplate%layout)
    DEALLOCATE(vectorTemplate%quantities)
    
    vectorTemplate%noQuantities=0
    vectorTemplate%noSubVectors=0
    vectorTemplate%totalElements=0
    vectorTemplate%name=''
    vectorTemplate%id=0
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
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,vectorTemplate%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//vectorTemplate%name)
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
    TYPE (Vector_T) :: vector

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
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,vector%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//vector%name)
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
    TYPE (Vector_T),  DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: l2gpIndex

    IF (ASSOCIATED(database)) THEN
       DO l2gpIndex=1,SIZE(database)
          CALL DestroyVectorInfo(database(l2gpIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyVectorDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine returns a pointer to a subVector within a vector

  FUNCTION GetSubVectorPointer(vector,quantity,profile,quantityName)

    ! Dummy arguments
    TYPE (Vector_T), INTENT(IN) :: vector
    INTEGER, INTENT(IN), OPTIONAL :: quantity
    INTEGER, INTENT(IN), OPTIONAL :: profile
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: quantityName

    ! Result
    REAL(r8), DIMENSION(:), POINTER :: GetSubVectorPointer

    ! Local variables
    INTEGER :: useQuantity,useProfile,subVector

    ! Executable code
    IF (PRESENT(quantityName)) THEN
       IF (PRESENT(quantity)) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Cannot use both quantity and quantityName in GetSubVectorPointer")
       useQuantity=LinearSearchStringArray(vector%template%quantities%name,&
            & quantityName,caseInsensitive=.TRUE.)
    ELSE
       IF (.NOT. PRESENT(quantity)) CALL MLSMessage(MLSMSG_Error,ModuleName,&
        & "Must supply either quantity or quantityName in GetSubVectorPointer")
       useQuantity=quantity
    END IF

    IF (PRESENT(profile)) THEN
       useProfile=profile
    ELSE
       useProfile=1
    ENDIF

    IF ((useQuantity<1).OR.(useQuantity>vector%template%noQuantities)) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName,"Invalid quantity request")

    IF ((useProfile<1).OR.&
         & (useProfile>vector%template%quantities(useQuantity)%noSubVectors)) &
         & CALL MLSMessage(MLSMSG_Error,ModuleName,"Invalid profile request")

    subVector=vector%template%layout(useQuantity)%entry(useProfile)
    GetSubVectorPointer=> &
         & vector%values(vector%template%subVectorFirstElement(subVector): &
         &   vector%template%subVectorFirstElement(subVector)+ &
         &   vector%template%subVectorNoElements(subVector)-1)
    
  END FUNCTION GetSubVectorPointer

  ! --------------------------------------------------------------------------

  ! This function uses the one above, and retuns a 2D array for a subVector
  ! The array defaults to (noChans,noSurfs), but if the firstIndexChannel
  ! flag is set .FALSE., the indices are (noSurfs,noChans)

  FUNCTION GetSubVectorAs2DArray(vector,quantity,profile,quantityName,&
       & firstIndexChannel)
    
    ! Dummy arguments
    TYPE (Vector_T), INTENT(IN) :: vector
    INTEGER, INTENT(IN), OPTIONAL :: quantity
    INTEGER, INTENT(IN), OPTIONAL :: profile
    CHARACTER (LEN=*), INTENT(IN), OPTIONAL :: quantityName
    LOGICAL, INTENT(IN), OPTIONAL :: firstIndexChannel

    ! Result
    REAL(r8), DIMENSION(:,:), POINTER :: GetSubVectorAs2DArray

    ! Local variables
    LOGICAL :: useFirstIndexChannel
    REAL(r8), DIMENSION(:), POINTER :: values
    REAL(r8), DIMENSION(:,:), POINTER :: tmpResult
    INTEGER :: useQuantity,useProfile,surf,chan,status
    TYPE (QuantityTemplate_T), POINTER :: qty

    ! Executable code

    IF (PRESENT(firstIndexChannel)) THEN
       useFirstIndexChannel=firstIndexChannel
    ELSE
       useFirstIndexChannel=.TRUE.
    END IF

    IF (PRESENT(quantityName)) THEN
       IF (PRESENT(quantity)) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Cannot use both quantity and quantityName in GetSubVectorPointer")
       useQuantity=LinearSearchStringArray(vector%template%quantities%name,&
            & quantityName,caseInsensitive=.TRUE.)
    ELSE
       IF (.NOT. PRESENT(quantity)) CALL MLSMessage(MLSMSG_Error,ModuleName,&
        & "Must supply either quantity or quantityName in GetSubVectorPointer")
       useQuantity=quantity
    END IF

    IF (PRESENT(profile)) THEN
       useProfile=profile
    ELSE
       useProfile=1
    ENDIF

    values=>GetSubVectorPointer(vector,useQuantity,useProfile)
    qty=>vector%template%quantities(useQuantity)

    IF (useFirstIndexChannel) THEN
       ALLOCATE(tmpResult(qty%noChans,qty%noSurfs),STAT=status)
    ELSE
       ALLOCATE(tmpResult(qty%noSurfs,qty%noChans),STAT=status)
    END IF
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "tmpResult")

    IF (useFirstIndexChannel.AND.qty%firstIndexChannel) THEN
       tmpResult=RESHAPE(values,(/qty%noChans,qty%noSurfs/))
    ELSE IF (useFirstIndexChannel.AND.(.NOT. qty%firstIndexChannel)) THEN
       tmpResult=TRANSPOSE(RESHAPE(values,(/qty%noSurfs,qty%noChans/)))
    ELSE IF ((.NOT. useFirstIndexChannel).AND.(qty%firstIndexChannel)) THEN
       tmpResult=TRANSPOSE(RESHAPE(values,(/qty%noChans,qty%noSurfs/)))
    ELSE
       tmpResult=RESHAPE(values,(/qty%noSurfs,qty%noChans/))
    END IF

    GetSubVectorAs2DArray=>tmpResult
  END FUNCTION GetSubVectorAs2DArray
       

!=============================================================================
END MODULE VectorsModule
!=============================================================================

!
! $Log$
! Revision 1.8  2000/04/14 20:27:43  vsnyder
! OOPS -- Replaced unspecified INTENT with INTENT(INOUT) -- see previous rev.
!
! Revision 1.7  2000/04/13 23:45:33  vsnyder
! Removed INTENT(IN) for the argument "vector" of DestroyVectorInfo, for
! which a component was deallocated
!
! Revision 1.6  2000/01/20 01:29:59  livesey
! Added GetSubVectorAs2DArray
!
! Revision 1.5  2000/01/19 18:36:46  livesey
! Sorted out the sub vector layout stuff.  Added the layout element of
! the template, and wrote code to deal with it.  Also wrote
! GetSubVectorPointer
!
! Revision 1.4  1999/12/17 21:43:17  livesey
! Added check for duplicate name
!
! Revision 1.3  1999/12/16 18:32:26  livesey
! Changed do to name change from VectorQuantities to QuantityTemplates
!
! Revision 1.2  1999/12/14 01:01:52  livesey
! Changed DOUBLE PRECISION to REAL(r8)
!
! Revision 1.1  1999/12/04 00:30:18  livesey
! First version.
!
!
