

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE GriddedData ! Collections of subroutines to handle TYPE GriddedData_T
!=============================================================================


  USE MLSCommon
  USE MLSMessageModule
  USE MLSStrings

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
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
   REAL (R8), POINTER, DIMENSION(:) :: heights 
             ! Surfaces (e.g. pressures etc.) [noHeights]
   !
   ! Now the latitudinal coordinate
   !
   LOGICAL :: equivalentLatitude ! If set, coordinate is equivalent latitude
   INTEGER :: noLats            ! Number of latitudes
   REAL (R8), POINTER, DIMENSION(:) :: lats ! Latitudes [noLats]
   !
   INTEGER :: noLons            ! Number of longitudes
   REAL (R8), POINTER, DIMENSION(:) :: lons ! Longitudes [noLons]
   !
   INTEGER noLsts               ! Number of local times
   REAL (R8), POINTER, DIMENSION(:) :: lsts ! Local times [noLsts]
   !
   INTEGER noSzas               ! Number of solar zenith angles
   REAL (R8), POINTER, DIMENSION(:) :: szas ! Zenith angles [noSzas]
   !
   INTEGER noDates              ! Number of dates in data
   REAL (R8), POINTER, DIMENSION(:) :: dateStarts
      ! Starting dates in SDP toolkit format
   REAL (R8), POINTER, DIMENSION(:) :: dateEnds
      ! Ending dates in SDP toolkit format
   !
   REAL (R8), POINTER, DIMENSION(:,:,:,:,:,:) :: field
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

  SUBROUTINE SetupNewGridTemplate(qty, source, noHeights, noLats, noLons, noLsts, noSzas, noDates)

    ! Dummy arguments
    TYPE (GriddedData_T), INTENT(OUT) :: qty ! Result

    TYPE (GriddedData_T), OPTIONAL, INTENT(IN) :: source ! Template

    INTEGER, OPTIONAL, INTENT(IN) :: noHeights, noLats, noLons, noLsts, noSzas, noDates

    ! Local variables
    INTEGER :: status           ! Status from allocates etc.

    ! Executable code

    ! First, if we have a template setup according to that
    IF (PRESENT(source)) THEN
       qty%noHeights=source%noHeights
       qty%noLats=source%noLats
       qty%noLons=source%noLons
       qty%noLsts=source%noLsts
       qty%noSzas=source%noSzas
       qty%noDates=source%noDates

      
    ELSE ! We have no template, setup a very bare quantity
       qty%noHeights=1
       qty%noLats=1
       qty%noLons=1
       qty%noLsts=1
       qty%noSzas=1
       qty%noDates=1

    ENDIF

    ! Now, see if the user asked for modifications to this
    IF (PRESENT(noHeights)) qty%noHeights=noHeights
    IF (PRESENT(noLats)) qty%noLats=noLats
    IF (PRESENT(noLons)) qty%noLons=noLons
    IF (PRESENT(noLsts)) qty%noLsts=noLsts
    IF (PRESENT(noSzas)) qty%noSzas=noSzas
    IF (PRESENT(noDates)) qty%noDates=noDates
    ! First the vertical coordinates

    ALLOCATE (qty%heights(qty%noHeights),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"heights")

    ! Now the geolocation coordinates
    ALLOCATE (qty%lats(qty%noLats),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lats")

    ALLOCATE (qty%lons(qty%noLons),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lons")

    ALLOCATE (qty%lsts(qty%noLsts),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lsts")

    ALLOCATE (qty%szas(qty%noSzas),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"szas")

    !Now the temporal coordinates
    ALLOCATE (qty%DateStarts(qty%noDates),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"DateStarts")

    ALLOCATE (qty%DateEnds(qty%noDates),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"DateEnds")

    !Now the data itself
    ALLOCATE(qty%field(qty%noHeights, qty%noLats, qty%noLons,  &
             qty%noLsts, qty%noSzas, qty%noDates), STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"field")


  END SUBROUTINE SetupNewGridTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  SUBROUTINE DestroyGridTemplateContents(qty)

    ! Dummy argument
    TYPE (GriddedData_T), INTENT(INOUT) :: qty
    ! Local variables
    INTEGER status

    ! Executable code

    DEALLOCATE (qty%heights, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"heights")

    DEALLOCATE (qty%lats, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"lats")

    DEALLOCATE (qty%lons, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"lons")

    DEALLOCATE (qty%lsts, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"lsts")

    DEALLOCATE (qty%szas, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"szas")

    DEALLOCATE (qty%DateStarts, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"DateStarts")

    DEALLOCATE (qty%DateEnds, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"DateEnds")

    DEALLOCATE (qty%field, STAT=status)

    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"field")    


  END SUBROUTINE DestroyGridTemplateContents

  ! --------------------------------------------------------------------------

  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

  SUBROUTINE AddGridTemplateToDatabase(database,qty)

    ! Dummy arguments
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database
    TYPE (GriddedData_T), INTENT(IN) :: qty

    ! Local variables
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: tempDatabase
    INTEGER :: newSize,status

    ! Executable code

    IF (ASSOCIATED(database)) THEN
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%quantityName, qty%quantityName, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//qty%quantityName)
       newSize=SIZE(database)+1
    ELSE
       newSize=1
    ENDIF
    ALLOCATE(tempDatabase(newSize),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Allocation failed for tempDatabase")

    IF (newSize>1) tempDatabase(1:newSize-1)=database
    tempDatabase(newSize)=qty
    IF (ASSOCIATED(database)) THEN
       DEALLOCATE(database, STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"database")
    end if
    database=>tempDatabase
  END SUBROUTINE AddGridTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyGridTemplateDatabase(database)

    ! Dummy argument
    TYPE (GriddedData_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: qtyIndex, status

    IF (ASSOCIATED(database)) THEN
       DO qtyIndex=1,SIZE(database)
          CALL DestroyGridTemplateContents(database(qtyIndex))
       ENDDO
       DEALLOCATE(database, stat=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_DeAllocate//"database")
    ENDIF
  END SUBROUTINE DestroyGridTemplateDatabase

!=============================================================================
END MODULE GriddedData
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 17:41:05  dcuddy
! Change revision to 2.0
!
! Revision 1.5  2000/06/20 22:19:22  lungu
! Changed DOUBBLE PRECISION to REAL (r8).
!
