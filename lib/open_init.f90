

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens and reads l2cf info, NCEP, DAO and climatology(ies) file
s
!=============================================================================


  USE MLSCommon
  USE MLSMessageModule
  USE MLSSignalNomenclature
  USE VGrid

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id:"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCS$"
  !-----------------------------------------------------------------------------

 

 

  ! First we'll define some global parameters and data types.

! This datatype stores a single gridded atmospheric quantity.  For example
! temperature, if an uncertainty field is also required, this is stored in a
! separate quantity.
!
! This type reflects the format of the Level 3 ASCII files, though note that
! these files can store multiple quantities such as these.
!
TYPE GriddedData_T
   !
   ! First the comment line(s) from the relevant input file
   !
   CHARACTER (LEN=LineLen), POINTER, DIMENSION(:) :: fileComments
   !
   ! Now the name, description and units information
   !
   CHARACTER (LEN=NameLen) :: quantityName ! From input file
   CHARACTER (LEN=LineLen) :: description ! Quantity description
   CHARACTER (LEN=NameLen) :: units ! Units for quantity
   !
   ! Now define the various coordinate systems, first vertical
   !
   INTEGER :: verticalCoordinate ! An 'enumerated' type
   INTEGER :: noHeights         ! Number of surfaces
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: heights 
             ! Surfaces (e.g. pressures etc.) [noHeights]
   !
   ! Now the latitudinal coordinate
   !
   LOGICAL :: equivalentLatitude ! If set, coordinate is equivalent latitude
   INTEGER :: noLats            ! Number of latitudes
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: lats ! Latitudes [noLats]
   !
   INTEGER :: noLons            ! Number of longitudes
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: lons ! Longitudes [noLons]
   !
   INTEGER noLsts               ! Number of local times
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: lsts ! Local times [noLsts]
   !
   INTEGER noSzas               ! Number of solar zenith angles
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: szas ! Zenith angles [noSzas]
   !
   INTEGER noDates              ! Number of dates in data
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: dateStarts
      ! Starting dates in SDP toolkit format
   DOUBLE PRECISION, POINTER, DIMENSION(:) :: dateEnds
      ! Ending dates in SDP toolkit format
   !
   DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:,:,:) :: field
   !
   ! The data itself.  This is stored as
   !  [noHeights, noLats, noLons, noLsts, noSzas, noDates]
   !
END TYPE GriddedData_T




  ! --------------------------------------------------------------------------

  CONTAINS

  ! Now we have some subroutines to deal with these quantitites

  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

  SUBROUTINE SetupNewGridTemplate(qty, source)

    ! Dummy arguments
    TYPE (GriddedData_T), INTENT(OUT) :: qty ! Result

    TYPE (GriddedData_T), OPTIONAL, INTENT(IN) :: source ! Template

    ! Local variables
    INTEGER :: status           ! Status from allocates etc.

    ! Executable code

    ! First, if we have a template setup according to that
    IF (PRESENT(source)) THEN
       qty%noHeights=source%noHeights
       qty%noLats=source%noLats
       qty%noLons=source%noLons
      
    ELSE ! We have no template, setup a very bare quantity
       qty%noHeights=1
       qty%noLats=1
       qty%noLons=1

    ENDIF

    ! Now, see if the user asked for modifications to this
    IF (PRESENT(noHeights)) qty%noHeights=noHeights
    IF (PRESENT(noLats)) qty%noLats=noLats
    IF (PRESENT(noLons)) qty%noLons=noLons


    ! First the vertical coordinates

    ALLOCATE (qty%surfs(qty%noLats,noHeightsToAllocate),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"surfs")

    ! Now the horizontal coordinates

    ALLOCATE (qty%phi(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"phi")

    ALLOCATE (qty%geodLat(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"geodLat")

    ALLOCATE (qty%lon(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lon")

    ALLOCATE (qty%time(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"time")

    ALLOCATE (qty%solarTime(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"solarTime")

    ALLOCATE (qty%solarZenith(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"solarZenith")

    ALLOCATE (qty%losAngle(noLatsToAllocate,qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"losAngle")

    ! Now some other stuff to allocate
    ALLOCATE (qty%subVectorIndex(qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"subVectorIndex")

    IF (.NOT. qty%regular) THEN
       ALLOCATE (qty%surfIndex(qty%subVectorLen,qty%noHeights),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"surfIndex")

       ALLOCATE (qty%chanIndex(qty%subVectorLen,qty%noHeights),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"chanIndex")
    ENDIF

    ! Se the id field and increment the counter
    qty%id=quantityTemplateCounter
    quantityTemplateCounter=quantityTemplateCounter+1
  END SUBROUTINE SetupNewGriddedDataTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  SUBROUTINE DestroyGriddedDataTemplateContents(qty)

    ! Dummy argument
    TYPE (GriddedData_T), INTENT(INOUT) :: qty

    ! Executable code

    DEALLOCATE (qty%surfs)
    DEALLOCATE (qty%phi)
    DEALLOCATE (qty%geodLat)
    DEALLOCATE (qty%lon)
    DEALLOCATE (qty%time)
    DEALLOCATE (qty%solarTime)
    DEALLOCATE (qty%solarZenith)
    DEALLOCATE (qty%losAngle)
    DEALLOCATE (qty%subVectorIndex)
    
    IF (.NOT. qty%regular) THEN
       DEALLOCATE (qty%surfIndex)
       DEALLOCATE (qty%chanIndex)
    ENDIF

  END SUBROUTINE DestroyGriddedDataTemplateContents

  ! --------------------------------------------------------------------------

  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

  SUBROUTINE AddGriddedDataTemplateToDatabase(database,qty)

    ! Dummy arguments
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database
    TYPE (GriddedData_T), INTENT(IN) :: qty

    ! Local variables
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,qty%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//qty%name)
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF

    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=qty
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddGriddedDataTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyGriddedDataTemplateDatabase(database)

    ! Dummy argument
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: qtyIndex

    IF (ASSOCIATED(database)) THEN
       DO qtyIndex=1,SIZE(database)
          CALL DestroyGriddedDataTemplateContents(database(qtyIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyGriddedDataTemplateDatabase

!=============================================================================
END MODULE OpenInit
!=============================================================================

!
! $Log$

!

