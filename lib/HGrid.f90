! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE HGrid                    ! Horizontal grid information
!=============================================================================

  USE MLSCommon
  USE Temporary_Types
  USE MLSMessage
  USE MLSStrings
  USE L1BData

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  ! ---------------------------------------------------------------------------

  ! This module contains datatypes and routines for handling HGrid information
  ! HGrids are the horizontal gridding information that get into vector
  ! quantities.

  ! This is the main datatype, an HGrid.

  TYPE HGrid_T
     CHARACTER (LEN=NameLen) :: name ! Name for HGrid
     INTEGER :: noProfs         ! Number of profiles in this grid

     ! Now the various coordinates in the HGrid, all dimensioned (noProfs)
     REAL(r8), DIMENSION(:), POINTER :: phi, geodLat, lon, time, &
          & solarTime,solarZenith,losAngle
   
     ! This is an array of profiles numbers to go into l2gp etc. data
     INTEGER, DIMENSION(:), POINTER :: profileIndices ! Dimensioned (noProfs)
  END TYPE HGrid_T

  CONTAINS

  ! ---------------------------------------------------------------------------

  ! This routine creates an empty HGrid

  SUBROUTINE SetupEmptyHGrid(hGrid,noProfs)
    
    ! Dummy arguments
    TYPE(HGrid_T), INTENT(OUT) :: hGrid
    INTEGER, INTENT(IN) :: noProfs

    ! Local variables
    INTEGER :: status

    ! Executable code

    hGrid%noProfs=noProfs
    ALLOCATE(hGrid%phi(noProfs), hGrid%geodLat(noProfs), hGrid%lon(noProfs), &
         & hGrid%time(noProfs), hGrid%solarTime(noProfs), &
         & hGrid%solarZenith(noProfs), hGrid%losAngle(noProfs), &
         & hGrid%profileIndices(noProfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,"hGrid information")
  END SUBROUTINE SetupEmptyHGrid


  ! --------------------------------------------------------------------------

  ! This routine creates an hGrid based on the user requests.
  
  SUBROUTINE CreateHGridFromMLSCFInfo(hGrid,cfInfo,l1bInfo,chunk)

    ! Dummy arguments
    TYPE (hGrid_T), INTENT(OUT) :: hGrid ! HGrid to create
    TYPE (L2CFEntry), INTENT(IN) :: cfInfo ! Instructions for its creation
    TYPE (L1BInfo_T), INTENT(IN) :: l1bInfo ! File handles for l1b data
    TYPE (Chunk_T), INTENT(IN) :: chunk ! This chunk

    ! Local parameters
    INTEGER, PARAMETER :: HG_Invalid=0
    INTEGER, PARAMETER :: HG_Height=1
    INTEGER, PARAMETER :: HG_Fractional=2

    INTEGER, PARAMETER :: NoL1BItemsToRead=8
    CHARACTER (LEN=15), DIMENSION(NoL1BItemsToRead), &
         &   PARAMETER :: L1bItemsToRead=&
         & (/"time           ","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpSolarTime    ","tpLosAngle     "/)
    INTEGER, PARAMETER :: TransitionToModularItems=2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz


    ! Local variables
    INTEGER :: hGridType=HG_Invalid
    INTEGER :: sampled=1.0
    INTEGER :: instrumentModule=MLSInstrumentModule_Invalid
    REAL(r8) :: fraction= -1e10, height= -1e10

    INTEGER :: keyNo            ! Entry in the l2cf line
    TYPE (L2CFCell) :: cell     ! Part of the l2cf information

    TYPE (L1BData_T) :: l1bItem ! L1B data
    REAL(r8), DIMENSION(:,:,:), POINTER :: tpGeodAngle,tpGeodAlt

    REAL(r8) :: minAngle,maxAngle,desiredAngle

    INTEGER :: noMAFs           ! Number of MAFs of L1B data read
    INTEGER :: maf,l1bItem      ! Loop counters
    INTEGER, DIMENSION(:), ALLOCATABLE :: defaultMIFs
    ! MIFs it would choose in the non over/undersampled case

    INTEGER :: status

    CHARACTER (LEN=NameLen) :: l1bItemName

    ! Executable code

    DO keyNo=1,cfInfo%l2cfEntryNoKeys
       cell=cfInfo%cells(keyNo)
       SELECT CASE(TRIM(cell%keyword))
       CASE("NAME")
          hGrid%name=cell%charValue
       CASE("TYPE")
          SELECT CASE(TRIM(Capitalize(cell%charValue)))
          CASE ("HEIGHT")
             hGridType=HG_Height
          CASE ("FRACTIONAL")
             hGridType=HG_Fractional
          CASE DEFAULT
             CALL MLSMessage(MLSMSG_Error,ModuleName,&
                  & "Unrecognised hGrid type: "//TRIM(cell%charValue))
          END SELECT
       CASE("MODULE")
          SELECT CASE(TRIM(Capitalize(cell%charValue)))
          CASE("GHZ")
             instrumentModule=MLSInstrumentModule_GHz
          CASE("THZ")
             instrumentModule=MLSInstrumentModule_THz
          CASE DEFAULT
             CALL MLSMessage(MLSMSG_Error,ModuleName,&
                  & "Unrecognised instrument module: "//TRIM(cell%charValue))
          END SELECT
       CASE("HEIGHT")
          ! Code needed here when l2cf finalized ***
       CASE("FRACTION")
          ! Code needed here when l2cf finalized ***
       CASE("SAMPLED")
          ! Code needed here when l2cf finalized ***
       CASE DEFAULT
          CALL MLSMessage(MLSMSG_Error,ModuleName,&
               & MLSMSG_Keyword//TRIM(cell%keyword))
       END SELECT
    END DO

    ! Now check the sanity of what we have

    IF (LEN_TRIM(hGrid%name)==0) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "Invalid or absent hGrid name")
    IF (hGridType==HG_Invalid) CALL MLSMessage(MLSMSG_Error,ModuleName, &
         & "No type specified for hGrid "//TRIM(hGrid%name))
    IF ((hGridType==HG_Height).AND.(height<-1e9)) CALL MLSMessage(&
         & MLSMSG_Error,ModuleName,"No height specified for hGrid "//&
         & TRIM(hGrid%name))
    IF ((hGridType==HG_Fractional).AND.(fraction<-1e9)) CALL MLSMessage(&
         & MLSMSG_Error,ModuleName,"No fraction specified for hGrid "//&
         & TRIM(hGrid%name))
    IF (instrumentModule==MLSModule_Invalid) CALL MLSMessage(&
         & MLSMSG_Error,ModuleName,"Must specifiy module=GHz/THz")

    ! This is where we will start reading the l1bdata get the name to read

    l1bItemName=MLSInstrumentModuleNames(instrumentModule)
    SELECT CASE(hGridType)
    CASE (HG_Fractional)
       l1bItemName=TRIM(l1bItemName)//"tpGeodAngle"
    CASE (HG_Height)
       l1bItemName=TRIM(l1bItemName)//"tpGeodAlt"
    END SELECT

    ! Read the data

    CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,noMAFs,l1bFlag,l1bField, &
         &firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex)
    IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_L1BRead//l1bItemName)
 
    IF (hGridType==HG_Fractional) THEN
       ! A fractional hGrid, we need to read the tangent point phi
       tpGeodAngle=>l1bField%dpField

       ALLOCATE(defaultMIFs(noMAFs),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//"defaultMIFs")

       ! Loop over the MAFs
       DO maf=1,noMAFs
          ! Think about missing data here! ***
          ! Probably need to do a pack on tpGeodAngle and then unpack on
          ! defaultMIFs

          minAngle=MINVAL(tpGeodAngle(1,:,maf))
          maxAngle=MAXVAL(tpGeodAngle(1,:,maf))

          defaultMIFs(maf)=Hunt(tpGeodAngle(1,:,maf),minAngle+&
               & fraction*(maxAngle-minAngle))
       END DO
    ELSE
       tpGeodAlt=>l1bField%dpField

       ALLOCATE(defaultMIFs(noMAFs),STAT=status)
       IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_Allocate//"defaultMIFs")

       ! Loop over the MAFs
       DO maf=1,noMAFs
          ! Think about missing data here! ***
          ! Probably need to do a pack on tpGeodAngle and then unpack on
          ! defaultMIFs

          defaultMIFs(maf)=Hunt(tpGeodAlt(1,:,maf),height)
       END DO
    ENDIF

    ! Done with this piece of l1b data for the moment
    CALL DeallocateL1BData(tpGeodAltData,flag)

    ! Now we have a default MIFs array, this is a list of the `standard'
    ! MIFs we would choose in the sampled=1 case.
    ! Work out how many profiles this is going to be

    noProfs=NINT(noMAFs*sampled)
    CALL SetupEmptyHgrid(hGrid,noProfs)

    ! Now we go through all the important geolocation quantities, read them
    ! in, interpolate them if required and store the result in the hGrid

    DO l1bItem=1,NoL1BItemsToRead
       ! Get the name of the item to read
       l1bItemName=l1bItemsToRead(l1bItem)
       IF (l1bItem>=TransitionToModularItems) l1bItemName=&
            & MLSIntrumentModuleNames(instrumentModule)//l1bItemName
       
       ! Read it from the l1boa file
       CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,noMAFsl, &
            & l1bFlag,l1bField, &
            & firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex)
       IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_L1BRead//l1bItemName)

       ! Got to here ******************************
    END DO

    
  END SUBROUTINE CreateHGridFromMLSCFInfo

  ! ---------------------------------------------------------------------------

  ! This routine destroys the information associated with an hGrid

  SUBROUTINE DestroyHGridContents(hGrid)

    ! Dummy arguments
    TYPE (HGrid_T), INTENT(OUT) :: hGrid

    ! Executable code

    DEALLOCATE(hGrid%phi, hGrid%geodLat, hGrid%lon, hGrid%time, &
         & hGrid%solarTime, hGrid%solarZenith, hGrid%losAngle, &
         & hGrid%profileIndices)

    hGrid%noProfs=0
    hGrid%name=""
  END SUBROUTINE DestroyHGridContents

  ! ---------------------------------------------------------------------------

  SUBROUTINE AddHGridToDatabase(database,hGrid)

    ! Dummy arguments
    TYPE (HGrid_T), DIMENSION(:), POINTER :: database
    TYPE (HGrid_T), INTENT(IN) :: hGrid

    ! Local variables
    TYPE (HGrid_T), DIMENSION(:), POINTER :: tempDatabase
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
    tempDatabase(newSize)=hGrid
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddHGridToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  SUBROUTINE DestroyHGridDatabase(database)

    ! Dummy argument
    TYPE (HGrid_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: hGridIndex

    IF (ASSOCIATED(database)) THEN
       DO hGridIndex=1,SIZE(database)
          CALL DestroyHGridContents(database(hGridIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyHGridDatabase

!=============================================================================
END MODULE HGrid
!=============================================================================

!
! $Log$
!
