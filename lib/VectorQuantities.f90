! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
MODULE VectorQuantities         ! Quantities within vectors
!=============================================================================

  USE MLSSignalNomenclature

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: id

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  !-----------------------------------------------------------------------------

  ! This module defines the `quantities' that make up vectors and their
  ! template information.

  ! First we'll define some global parameters and data types.

  INTEGER, PARAMETER :: QtyNameLen=20

  ! This data type defines a quantity template

  TYPE QuantityTemplate_T

     ! First some administrative stuff

     CHARACTER (LEN=QtyNameLen) :: name ! Simple name for quantity
     INTEGER :: id              ! Id code for quantity (for checking stuff)

     ! Now the dimensions of this quantity

     INTEGER :: noSubVectors    ! Number of subVectors in this quantity
     INTEGER :: noSurfs         ! Number of surfaces per subvector 
     INTEGER :: noChans         ! Number of channels

     ! Now some flags describing the quantity

     LOGICAL :: coherent        ! Do subvectors have same vertical coordiantes?
     LOGICAL :: stacked         ! Are subvectors true vertical profiles?
     LOGICAL :: regular         ! Are all channels/heights represented

     ! Now the vertical coordinate

     INTEGER :: verticalCoordinate ! The vertical coordinate used

     ! Now some misc. information

     DOUBLE PRECISION :: badValue ! Value used to flag bad/missing data
     DOUBLE PRECISION :: scaleFactor ! Scale factor used when printing etc.

     ! Now, for regular quantities the number of elements each subvector is
     ! simply noSurfs*noChans.  For irregular ones it is less, but it is
     ! constant from subvector to subvector, this is that number

     INTEGER :: subVectorLen

     ! Now we define how the data in each subvector is stored whether by
     ! surface or channels (regular quantities only).

     INTEGER :: surfStride, chanStride

     ! Now we give the vertical coordinates

     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: surfs

     ! This is dimensioned (noSurfs,1) for regular quantities and (noSurfs
     !,noSubVectors) for irregular ones.

     ! Now the horizontal coordinates

     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: phi

     ! This is dimensioned (1,noSubVectors) for stacked quantities and (noSurfs
     !,noSubVectors) for unstacked ones.  These other coordinates are
     ! dimensioned in the same manner.

     DOUBLE PRECISION, DIMENSION(:,:), POINTER :: geodLat,lon,time, &
          & solarTime,solarZenith,losAngle

     ! This integer array can be used to tie subVectors to say l2gp profile
     ! numbers or to MAF indices

     INTEGER, DIMENSION(:), POINTER :: subVectorIndex

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! For quantities containing `channels' the following information may or
     ! may not be useful.

     ! Some quantities are on abritrary freqency grids, these quantities refer
     ! to those.

     INTEGER :: frequencyCoordinate ! An enumerated type, explains next fields

     DOUBLE PRECISION, DIMENSION(:), POINTER :: frequencies 
     ! List of frequencies (noChans)

     DOUBLE PRECISION :: lo     ! Local oscillator (optional)

     ! Other quantities refer to a particular MLS signal designation

     TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal ! Signal (optional)
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ! Now for irregular quantities, instead of using the stride information,
     ! we have these arrays to help us navigate around the quantity.

     INTEGER, DIMENSION(:,:), POINTER :: surfIndex
     INTEGER, DIMENSION(:,:), POINTER :: chanIndex
     ! These are actually dimensioned (subVectorLen,noSubVectors)
  END TYPE QuantityTemplate_T

  ! --------------------------------------------------------------------------

  CONTAINS

  ! Now we have some subroutines to deal with these quantitites

  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

  SUBROUTINE SetupNewQuantityTemplate(qty, source, noSubVectors, noSurfs, &
       & noChans, coherent, stacked, regular, subVectorLen, storeByChannel)

    ! Dummy arguments
    TYPE (QuantityTemplate_T), INTENT(OUT) :: qty ! Result

    TYPE (QuantityTemplate_T), OPTIONAL, INTENT(IN) :: source ! Template
    INTEGER, INTENT(IN), OPTIONAL :: noSubVectors
    INTEGER, INTENT(IN), OPTIONAL :: noSurfs
    INTEGER, INTENT(IN), OPTIONAL :: noChans
    LOGICAL, INTENT(IN), OPTIONAL :: coherent
    LOGICAL, INTENT(IN), OPTIONAL :: stacked
    LOGICAL, INTENT(IN), OPTIONAL :: regular
    INTEGER, INTENT(IN), OPTIONAL :: subVectorLen
    LOGICAL, INTENT(IN), OPTIONAL :: storeByChannel

    ! Local parameters
    CHARACTER (LEN=*), PARAMETER :: AllocateMessage= &
         & "Allocation failed for "

    ! Local variables
    INTEGER :: status           ! Status from allocates etc.
    LOGICAL :: useStoreByChannel ! Copy of store by channel
    INTEGER :: noSurfsToAllocate ! For allocations
    INTEGER :: noSubVectorsToAllocate ! For allocations

    ! Executable code

    ! First, if we have a template setup according to that
    IF (PRESENT(source)) THEN
       qty%noSubVectors=source%noSubVectors
       qty%noSurfs=source%noSurfs
       qty%noChans=source%noChans
       qty%coherent=source%coherent
       qty%stacked=source%stacked
       qty%regular=source%regular
       qty%subVectorLen=source%subVectorLen
    ELSE ! We have no template, setup a very bare quantity
       qty%noSubVectors=1
       qty%noSurfs=1
       qty%noChans=1
       qty%coherent=.TRUE.
       qty%stacked=.TRUE.
       qty%regular=.TRUE.
       qty%subVectorLen=1
    ENDIF

    ! Now, see if the user asked for modifications to this
    IF (PRESENT(noSubVectors)) qty%noSubVectors=noSubVectors
    IF (PRESENT(noSurfs)) qty%noSurfs=noSurfs
    IF (PRESENT(noChans)) qty%noChans=noChans
    IF (PRESENT(coherent)) qty%coherent=coherent
    IF (PRESENT(stacked)) qty%stacked=stacked
    IF (PRESENT(regular)) qty%regular=regular

    ! Now think about subVectorLen
    IF ((.NOT. qty%regular).AND.(PRESENT(subVectorLen))) THEN
       qty%subVectorLen=subVectorLen
    ELSE
       qty%subVectorLen=qty%noSurfs*qty%noChans
    ENDIF

    ! Deal with the storeByChannel argument
    useStoreByChannel=.FALSE.
    IF (PRESENT(storeByChannel)) useStoreByChannel=storeByChannel

    ! Now think about strides
    IF (qty%regular) THEN
       IF (useStoreByChannel) THEN
          qty%surfStride=qty%noChans
          qty%chanStride=1
       ELSE
          qty%chanStride=qty%noSurfs
          qty%surfStride=1
       ENDIF
    ELSE
       qty%surfStride=0
       qty%chanStride=0
    ENDIF

    ! Now we allocate all the arrays we're going to need

    IF (qty%coherent) THEN 
       noSubVectorsToAllocate=1
    ELSE
       noSubVectorsToAllocate=qty%noSubVectors
    ENDIF

    IF (qty%stacked) THEN
       noSurfsToAllocate=1
    ELSE
       noSurfsToAllocate=qty%noSurfs
    ENDIF

    ! First the vertical coordinates

    ALLOCATE (qty%surfs(qty%noSurfs,noSubVectorsToAllocate),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"surfs",error=.TRUE.)

    ! Now the horizontal coordinates

    ALLOCATE (qty%phi(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"phi",error=.TRUE.)

    ALLOCATE (qty%geodLat(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"geodLat",error=.TRUE.)

    ALLOCATE (qty%lon(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"lon",error=.TRUE.)

    ALLOCATE (qty%time(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"time",error=.TRUE.)

    ALLOCATE (qty%solarTime(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"solarTime",error=.TRUE.)

    ALLOCATE (qty%solarZenith(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"solarZenith",error=.TRUE.)

    ALLOCATE (qty%losAngle(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"losAngle",error=.TRUE.)

    ! Now some other stuff to allocate
    ALLOCATE (qty%subVectorIndex(qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(AllocateMessage//"subVectorIndex", &
         & error=.TRUE.)

    IF (.NOT. qty%regular) THEN
       ALLOCATE (qty%surfIndex(qty%subVectorLen,qty%noSubVectors),STAT=status)
       IF (status /= 0) CALL MLSMessage(AllocateMessage//"surfIndex", &
         & error=.TRUE.)

       ALLOCATE (qty%chanIndex(qty%subVectorLen,qty%noSubVectors),STAT=status)
       IF (status /= 0) CALL MLSMessage(AllocateMessage//"chanIndex", &
         & error=.TRUE.)
    ENDIF
  END SUBROUTINE CreateNewQuantityTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  SUBROUTINE DestroyQuantityTemplateContents(qty)

    ! Parameters
    TYPE (QuantityTemplate_T) :: qty

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
    DEALLOCATE (qty%subVectorLen)
    
    IF (.NOT. qty%regular) THEN
       DEALLOCATE (qty%surfIndex)
       DEALLOCATE (qty%chanIndex)
    ENDIF

  END SUBROUTINE DestroyQuantityTemplateContents

!=============================================================================
END MODULE VectorQuantities
!=============================================================================

!
! $Log$
!

