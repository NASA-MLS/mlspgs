! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ConstructVectorTemplates ! Construct a template for a vector
!=============================================================================

  USE MLSCommon
  USE MLSMessageModule
  USE VectorsModule
  USE QuantityTemplates
  USE MLSCF

  IMPLICIT NONE
  PUBLIC
  
  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module performs the vector template aspects of the construct task

CONTAINS

  ! This is the main routine for this module

  SUBROUTINE CreateVecTemplateFromMLSCfInfo(cfInfo,vectorTemplate,&
       & quantityTemplates)

    ! Dummy arguments
    TYPE(mlscfEntry_T), INTENT(IN) :: cfInfo
    TYPE (VectorTemplate_T), INTENT(OUT) :: vectorTemplate
    TYPE (QuantityTemplate_T), DIMENSION(:) :: quantityTemplates

    ! Local variables
    INTEGER :: keyNo,status,quantityNo
    CHARACTER (LEN=NameLen) :: name=""
    LOGICAL, DIMENSION(:), ALLOCATABLE :: selected
    TYPE (MLSCFCell_T) :: cell

    ! Executable code

    ALLOCATE(selected(SIZE(quantityTemplates)),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "selected")
    selected=.FALSE.

    ! Simply loop through the mlscf information supplied.
    !  If it's a NAME= then copy it to name, else it's a
    ! quantityTemplate name.

    DO keyNo=1,cfInfo%mlscfEntryNoKeys
       cell=cfInfo%cells(keyNo)
       IF (TRIM(cell%keyword)=="NAME") THEN
          name=cell%charValue
       ELSE
          quantityNo=LinearSearchStringArray(quantityTemplates%name,&
               & cell%charValue,caseInsensitive=.TRUE.)
          IF (quantityNo==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & "No such quantityTemplate: "//cell%charValue)
          selected(quantityNo)=.TRUE.
       ENDIF
    END DO

    IF (LEN_TRIM(name)==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "No name given for vector")
    IF (COUNT(selected)==0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "No quantities chosen for vector: "//name)

    ! Now finally construct the vector template

    CALL ConstructVectorTemplate(name,PACK(quantityTemplates,selected),&
         & vectorTemplate)

  END SUBROUTINE CreateVecTemplateFromMLSCfInfo

!=============================================================================
END MODULE ConstructVectorTemplates
!=============================================================================

!
! $Log$
! Revision 1.1  1999/12/18 03:00:45  livesey
! First version
!
!
