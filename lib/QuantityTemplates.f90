! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7407 is acknowledged.

!=============================================================================
MODULE QuantityTemplates         ! Quantities within vectors
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
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  ! This module defines the `quantities' that make up vectors and their
  ! template information.

  ! This set of integers define various families that quantities can fall into
  ! No doubt more will be added.

  INTEGER, PARAMETER :: NoQTYTypes=7
  CHARACTER (LEN=12), PARAMETER, DIMENSION(NoQtyTypes) :: &
       & QTYTypeNames= (/ &
       & "TEMPERATURE ", &
       & "VMR         ", &
       & "RADIANCE    ", &
       & "PTAN        ", &
       & "BASELINE    ", &
       & "EXTINCTION  ", &
       & "GPH         "/)
  INTEGER, PARAMETER :: QTY_Invalid=0
  INTEGER, PARAMETER :: QTY_Temperature=1
  INTEGER, PARAMETER :: QTY_Vmr=2
  INTEGER, PARAMETER :: QTY_Radiance=3
  INTEGER, PARAMETER :: QTY_Ptan=4
  INTEGER, PARAMETER :: QTY_Baseline=5
  INTEGER, PARAMETER :: QTY_Extinction=6
  INTEGER, PARAMETER :: QTY_Gph=7

  ! This set of integers defines how quantities are broken into `channels'
  ! These are known as FGrid information.  Rather than have an FGrid module,
  ! for the moment we'll just handle this information here. This may be split
  ! out in later versions of the code.

  INTEGER, PARAMETER :: NoFGTypes=5
  CHARACTER (LEN=24), PARAMETER, DIMENSION(NoFGTypes) :: &
       & FGTypeNames= (/ &
       & "No frequency dependence ", &
       & "MLS Channel             ", &
       & "Intermediate Frequecny  ", &
       & "Upper Sideband Frequency", &
       & "Lower Sideband Frequency"/)
  INTEGER, PARAMETER :: FG_Invalid=0
  INTEGER, PARAMETER :: FG_None=1
  INTEGER, PARAMETER :: FG_InstrumentChannel=2
  INTEGER, PARAMETER :: FG_IntermediateFrequency=3
  INTEGER, PARAMETER :: FG_USBFrequency=4
  INTEGER, PARAMETER :: FG_LSBFrequency=5

  ! First we'll define some global parameters and data types.

  TYPE QuantityTemplate_T

     ! First some administrative stuff

     CHARACTER (LEN=NameLen) :: name ! Simple name for quantity
     INTEGER :: id              ! Id code for quantity (for checking stuff)

     ! This integer is an enumerated type describing what kind of quantity this
     ! is. e.g. QTY_Temperature

     INTEGER :: quantityType

     ! Now the dimensions of this quantity

     INTEGER :: noSubVectors    ! Number of subVectors in this quantity
     INTEGER :: noSurfs         ! Number of surfaces per subvector 
     INTEGER :: noChans         ! Number of channels

     ! Now some flags describing the quantity

     LOGICAL :: coherent        ! Do subvectors have same vertical coordiantes?
     LOGICAL :: stacked         ! Are subvectors true vertical profiles?
     LOGICAL :: regular         ! Are all channels/heights represented

     ! This next one allows software using the vector quantities to be somewhat
     ! lazy and, for example, avoid interpolation.  Minor frame quantities are
     ! incoherent and unstacked, but may be regular or irregular.  However, not
     ! all incoherent unstacked quantities are minor frame quantities.

     LOGICAL :: minorFrame      ! Is this a minor frame quantity.

     ! This information describes how much of the data is in the overlap
     ! regions if any.

     INTEGER :: noSubVectorsLowerOverlap
     INTEGER :: noSubVectorsUpperOverlap

     ! Now the vertical coordinate

     INTEGER :: verticalCoordinate ! The vertical coordinate used

     ! Now some misc. information

     REAL(r8) :: badValue ! Value used to flag bad/missing data
     CHARACTER(LEN=NameLen) :: unit ! Unit quantity is in when scaled as below
     REAL(r8) :: scaleFactor ! Scale factor used when printing etc.

     ! Now, for regular quantities the number of elements each subvector is
     ! simply noSurfs*noChans.  For irregular ones it is less, but it is
     ! constant from subvector to subvector, this is that number

     INTEGER :: subVectorLen

     ! Now we define how the data in each subvector is stored whether by
     ! surface or channels (regular quantities only).

     LOGICAL :: firstIndexChannel

     ! Now we give the vertical coordinates

     REAL(r8), DIMENSION(:,:), POINTER :: surfs

     ! This is dimensioned (noSurfs,1) for coherent quantities and (noSurfs
     !,noSubVectors) for incoherent ones.

     ! Now the horizontal coordinates

     REAL(r8), DIMENSION(:,:), POINTER :: phi

     ! This is dimensioned (1,noSubVectors) for stacked quantities and (noSurfs
     !,noSubVectors) for unstacked ones.  These other coordinates are
     ! dimensioned in the same manner.

     REAL(r8), DIMENSION(:,:), POINTER :: geodLat,lon,time, &
          & solarTime,solarZenith,losAngle

     ! These integer arrays are used for minor frame quantities to index the
     ! major frames.

     INTEGER, DIMENSION(:), POINTER :: mafIndex,mafCounter

     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! For quantities containing `channels' the following information may or
     ! may not be useful.

     ! Some quantities are on abritrary freqency grids, these quantities refer
     ! to those.

     INTEGER :: frequencyCoordinate ! An enumerated type, e.g. FG_USBFreq
     REAL(r8), DIMENSION(:), POINTER :: frequencies 
     ! List of frequencies (noChans)

     REAL(r8) :: lo     ! Local oscillator (optional)
     
     TYPE (MLSSignal_T), DIMENSION(:), POINTER :: signal ! Signal (optional)
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     ! Some families of quantities require special additional information.
     ! This is given here if needed.

     INTEGER :: radiometerIndex ! Which radiometer does a ptan qty refer to?
     CHARACTER(LEN=NameLen) :: molecule ! What molecule does this refer to?
     
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     ! Now for irregular quantities, instead of using the firstIndexChannel
     ! information, we have these arrays to help us navigate around the
     ! quantity.

     INTEGER, DIMENSION(:,:), POINTER :: surfIndex
     INTEGER, DIMENSION(:,:), POINTER :: chanIndex
     ! These are actually dimensioned (subVectorLen,noSubVectors)
  END TYPE QuantityTemplate_T

  ! This incrementing counter is used to set the id field of a quantity template

  INTEGER, PRIVATE :: quantityTemplateCounter=0

  ! --------------------------------------------------------------------------

  CONTAINS

  ! Now we have some subroutines to deal with these quantitites

  ! This first routine sets up a new quantity template according to the user
  ! input.  This may be based on a previously supplied template (with possible
  ! modifications), or created from scratch.

  SUBROUTINE SetupNewQuantityTemplate(qty, source, noSubVectors, noSurfs, &
       & noChans, coherent, stacked, regular, subVectorLen, firstIndexChannel, &
       & minorFrame)

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
    LOGICAL, INTENT(IN), OPTIONAL :: firstIndexChannel
    LOGICAL, INTENT(IN), OPTIONAL :: minorFrame

    ! Local variables
    INTEGER :: status           ! Status from allocates etc.
    LOGICAL :: useFirstindexchannel ! Copy of store by channel
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
       qty%minorFrame=source%minorFrame
       qty%subVectorLen=source%subVectorLen
    ELSE ! We have no template, setup a very bare quantity
       qty%noSubVectors=1
       qty%noSurfs=1
       qty%noChans=1
       qty%coherent=.TRUE.
       qty%stacked=.TRUE.
       qty%regular=.TRUE.
       qty%minorFrame=.FALSE.
       qty%subVectorLen=1
    ENDIF

    ! Now, see if the user asked for modifications to this
    IF (PRESENT(noSubVectors)) qty%noSubVectors=noSubVectors
    IF (PRESENT(noSurfs)) qty%noSurfs=noSurfs
    IF (PRESENT(noChans)) qty%noChans=noChans
    IF (PRESENT(regular)) qty%regular=regular
    IF (PRESENT(minorFrame)) qty%minorFrame=minorFrame
    IF (qty%minorFrame) THEN
       IF (PRESENT(coherent)) THEN
          IF (coherent) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Minor frame quantities must be incoherent")
       ENDIF
       qty%coherent=.FALSE.
       IF (PRESENT(stacked)) THEN
          IF (stacked) CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Minor frame quantities must be unstacked")
       ENDIF
       qty%stacked=.FALSE.
    ELSE
       IF (PRESENT(coherent)) qty%coherent=coherent
       IF (PRESENT(stacked)) qty%stacked=stacked
    ENDIF

    ! Now think about subVectorLen
    IF ((.NOT. qty%regular).AND.(PRESENT(subVectorLen))) THEN
       qty%subVectorLen=subVectorLen
    ELSE
       qty%subVectorLen=qty%noSurfs*qty%noChans
    ENDIF

    ! Deal with the firstindexchannel argument
    useFirstindexchannel=.TRUE.
    IF (PRESENT(firstindexchannel)) useFirstindexchannel=firstindexchannel
    qty%firstIndexChannel=useFirstIndexChannel

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
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"surfs")

    ! Now the horizontal coordinates

    ALLOCATE (qty%phi(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"phi")

    ALLOCATE (qty%geodLat(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"geodLat")

    ALLOCATE (qty%lon(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"lon")

    ALLOCATE (qty%time(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"time")

    ALLOCATE (qty%solarTime(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"solarTime")

    ALLOCATE (qty%solarZenith(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"solarZenith")

    ALLOCATE (qty%losAngle(noSurfsToAllocate,qty%noSubVectors),STAT=status)
    IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & MLSMSG_Allocate//"losAngle")

    ! Now some other stuff to allocate
    IF (qty%minorFrame) THEN
       ALLOCATE (qty%MAFIndex(qty%noSubVectors), &
            & qty%MAFCounter(qty%noSubVectors), &
            & STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"MAFIndex/MAFCounter")
    END IF

    IF (.NOT. qty%regular) THEN
       ALLOCATE (qty%surfIndex(qty%subVectorLen,qty%noSubVectors),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"surfIndex")

       ALLOCATE (qty%chanIndex(qty%subVectorLen,qty%noSubVectors),STAT=status)
       IF (status /= 0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
            & MLSMSG_Allocate//"chanIndex")
    ENDIF

    ! Se the id field and increment the counter
    qty%id=quantityTemplateCounter
    quantityTemplateCounter=quantityTemplateCounter+1
  END SUBROUTINE SetupNewQuantityTemplate

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template

  SUBROUTINE DestroyQuantityTemplateContents(qty)

    ! Dummy argument
    TYPE (QuantityTemplate_T), INTENT(INOUT) :: qty

    ! Executable code

    DEALLOCATE (qty%surfs)
    DEALLOCATE (qty%phi)
    DEALLOCATE (qty%geodLat)
    DEALLOCATE (qty%lon)
    DEALLOCATE (qty%time)
    DEALLOCATE (qty%solarTime)
    DEALLOCATE (qty%solarZenith)
    DEALLOCATE (qty%losAngle)

    IF (qty%minorFrame) DEALLOCATE (qty%MAFIndex,qty%MAFCounter)
    
    IF (.NOT. qty%regular) THEN
       DEALLOCATE (qty%surfIndex)
       DEALLOCATE (qty%chanIndex)
    ENDIF

  END SUBROUTINE DestroyQuantityTemplateContents

  ! --------------------------------------------------------------------------

  ! This subroutine adds a quantity template to a database, or creates the
  ! database if it doesn't yet exist

  SUBROUTINE AddQuantityTemplateToDatabase(database,qty)

    ! Dummy arguments
    TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: database
    TYPE (QuantityTemplate_T), INTENT(IN) :: qty

    ! Local variables
    TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: tempDatabase
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
    IF (ASSOCIATED(database))DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddQuantityTemplateToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyQuantityTemplateDatabase(database)

    ! Dummy argument
    TYPE (QuantityTemplate_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: qtyIndex

    IF (ASSOCIATED(database)) THEN
       DO qtyIndex=1,SIZE(database)
          CALL DestroyQuantityTemplateContents(database(qtyIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyQuantityTemplateDatabase

!=============================================================================
END MODULE QuantityTemplates
!=============================================================================

!
! $Log$
! Revision 1.10  2000/05/15 22:52:35  livesey
! Typo fix.
!
! Revision 1.9  2000/01/20 21:59:28  livesey
! Replaced subVectorIndex with MAFIndex and MAFCounter
!
! Revision 1.8  2000/01/20 01:28:21  livesey
! Removed the horizontal coordinate information, and beefed up the
! frequency coordinate information.
!
! Revision 1.7  2000/01/18 21:27:00  livesey
! Added noSubVectors(Lower/Upper)Overlap, copied from HGrid or similar.
!
! Revision 1.6  2000/01/12 20:55:49  livesey
! Added minorFrame flag.
!
! Revision 1.5  2000/01/07 23:53:35  livesey
! Nearly integrated, just a few tweaks.
!
! Revision 1.4  1999/12/17 21:42:16  livesey
! Added check for duplicate name
!
! Revision 1.3  1999/12/16 23:46:35  livesey
! Added the unit field and fixed a couple of compile glitches
!
! Revision 1.2  1999/12/16 23:12:09  livesey
! Added quantityType and other support entries
!
! Revision 1.1  1999/12/16 18:31:43  livesey
! First version. Renamed from VectorQuantities
!
!
! This module was previously known as VectorQuantities.  This is it's previous
! revision history.
!
! Revision 1.9  1999/12/16 01:28:02  livesey
! Routine checkin
!
! Revision 1.8  1999/12/14 00:55:29  livesey
! Changed DOUBLE PRECISION to REAL(r8)
!
! Revision 1.7  1999/12/04 00:26:33  livesey
! Added a few comments.
!
! Revision 1.6  1999/12/03 22:27:08  livesey
! Tidied up some of the INTENT stuff
!
! Revision 1.5  1999/12/03 21:57:34  livesey
! Added the code to set the id field with an incrementing counter
!
! Revision 1.4  1999/12/01 23:01:41  livesey
! Before renaming things to upper/lower case
!
! Revision 1.3  1999/12/01 05:03:56  livesey
! Nightly checkin
!
! Revision 1.2  1999/11/30 04:03:51  livesey
! Bug fix
!
! Revision 1.1  1999/11/24 23:06:33  livesey
! First simple version.
!
!

