! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE HGrid                    ! Horizontal grid information
!=============================================================================

  USE MLSCommon
  USE MLSCF
  USE MLSMessageModule
  USE MLSStrings
  USE L1BData
  USE MLSNumerics

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
     INTEGER :: noProfsLowerOverlap ! Number of profiles in the lower overlap
     INTEGER :: noProfsUpperOverlap ! Number of profiles in the upper overlap

     ! Now the various coordinates in the HGrid, all dimensioned (noProfs)
     REAL(r8), DIMENSION(:), POINTER :: phi, geodLat, lon, time, &
          & solarTime,solarZenith,losAngle
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
         & STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,"hGrid information")

    hGrid%noProfsLowerOverlap=0
    hGrid%noProfsUpperOverlap=0
  END SUBROUTINE SetupEmptyHGrid


  ! --------------------------------------------------------------------------

  ! This routine creates an hGrid based on the user requests.
  
  SUBROUTINE CreateHGridFromMLSCFInfo(hGrid,cfInfo,l1bInfo,chunk)

    ! Dummy arguments
    TYPE (hGrid_T), INTENT(OUT) :: hGrid ! HGrid to create
    TYPE (MLSCFEntry_T), INTENT(IN) :: cfInfo ! Instructions for its creation
    TYPE (L1BInfo_T), INTENT(IN) :: l1bInfo ! File handles for l1b data
    TYPE (MLSChunk_T), INTENT(IN) :: chunk ! This chunk

    ! Local parameters
    INTEGER, PARAMETER :: HG_Invalid=0
    INTEGER, PARAMETER :: HG_Height=1
    INTEGER, PARAMETER :: HG_Fractional=2

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! IF MODIFYING THIS SECTION PLEASE TAKE CARE, SEE BELOW!
    INTEGER, PARAMETER :: NoL1BItemsToRead=7
    CHARACTER (LEN=15), DIMENSION(NoL1BItemsToRead), &
         &   PARAMETER :: L1bItemsToRead=&
         & (/"MAFStartTimeTAI","tpGeodLat      ","tpLon          ",&
         &   "tpGeodAngle    ","tpSolarZenith  ","tpSolarTime    ",&
         &   "tpLosAngle     "/)
    INTEGER, PARAMETER :: TransitionToModularItems=2
    ! Entries in the above array below TransitionToModularItems are prefixed
    ! with either GHz or THz.  The layout of the above array is critically
    ! bound to the SELECT CASE(l1bItem) code below.  So TAKE CARE! when
    ! modifing it.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    ! Local variables
    INTEGER :: hGridType=HG_Invalid
    REAL(r8) :: interpolationFactor=1.0
    INTEGER :: instrumentModule=MLSInstrumentModule_Invalid
    REAL(r8) :: fraction= -1e10, height= -1e10

    INTEGER :: keyNo            ! Entry in the mlscf line
    TYPE (MLSCFCell_T) :: cell  ! Part of the mlscf information

    TYPE (L1BData_T) :: l1bField ! L1B data
    REAL(r8), DIMENSION(:,:,:), POINTER :: tpGeodAngle,tpGeodAlt

    REAL(r8) :: minAngle,maxAngle,desiredAngle

    INTEGER :: noMAFs           ! Number of MAFs of L1B data read
    INTEGER :: noProfs          ! Number of profiles in output hGrid
    INTEGER :: maf,l1bItem,prof ! Loop counters
    INTEGER, DIMENSION(:), ALLOCATABLE :: defaultMIFs
    ! MIFs it would choose in the non over/undersampled case
    REAL(r8), DIMENSION(:), ALLOCATABLE :: defaultField,interpolatedField
    REAL(r8), DIMENSION(:), ALLOCATABLE :: intermediateField

    INTEGER :: status,l1bFlag   ! Flags

    CHARACTER (LEN=NameLen) :: l1bItemName

    ! Executable code

    DO keyNo=1,cfInfo%mlscfEntryNoKeys
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
          ! Code needed here when mlscf finalized ***
       CASE("FRACTION")
          ! Code needed here when mlscf finalized ***
       CASE("INTERPOLATIONFACTOR")
          ! Code needed here when mlscf finalized ***
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
    IF (instrumentModule==MLSInstrumentModule_Invalid) CALL MLSMessage(&
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

    CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,l1bField,noMAFs,l1bFlag, &
         & firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex)
    IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_L1BRead//l1bItemName)

    ! Allocate default MIFs

    ALLOCATE(defaultMIFs(noMAFs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & MLSMSG_Allocate//"defaultMIFs")
   
    ! Work out which MIF should have the profile for each MAF.

    IF (hGridType==HG_Fractional) THEN
       ! A fractional hGrid, we need to read the tangent point phi
       tpGeodAngle=>l1bField%dpField

       ! Loop over the MAFs
       DO maf=1,noMAFs
          ! Think about missing data here! ***
          ! Probably need to do a pack on tpGeodAngle and then unpack on
          ! defaultMIFs

          minAngle=MINVAL(tpGeodAngle(1,:,maf))
          maxAngle=MAXVAL(tpGeodAngle(1,:,maf))

          CALL Hunt(tpGeodAngle(1,:,maf),minAngle+&
               & fraction*(maxAngle-minAngle),defaultMIFs(maf))
       END DO
    ELSE
       tpGeodAlt=>l1bField%dpField

       ! Loop over the MAFs
       DO maf=1,noMAFs
          ! Think about missing data here! ***
          ! Probably need to do a pack on tpGeodAngle and then unpack on
          ! defaultMIFs

          CALL Hunt(tpGeodAlt(1,:,maf),height,defaultMIFs(maf))
       END DO
    ENDIF

    ! Done with this piece of l1b data for the moment
    CALL DeallocateL1BData(l1bField,l1bFlag)

    ! Now we have a default MIFs array, this is a list of the `standard'
    ! MIFs we would choose in the interpolationFactor=1 case.
    ! Work out how many profiles this is going to be

    noProfs=NINT(noMAFs*interpolationFactor)
    CALL SetupEmptyHgrid(hGrid,noProfs)

    ! Setup some arrays
    ALLOCATE(defaultField(noMAFs),interpolatedField(noProfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         & "defaultField and/or interpolatedField")

    ! Now we go through all the important geolocation quantities, read them
    ! in, interpolate them if required and store the result in the hGrid

    DO l1bItem=1,NoL1BItemsToRead
       ! Get the name of the item to read
       l1bItemName=l1bItemsToRead(l1bItem)
       IF (l1bItem>=TransitionToModularItems) l1bItemName=&
            & MLSInstrumentModuleNames(instrumentModule)//l1bItemName
       
       ! Read it from the l1boa file
       CALL ReadL1BData(l1bInfo%l1boaid,l1bItemName,l1bField,noMAFs, &
            & l1bFlag,firstMAF=chunk%firstMafIndex, lastMAF=chunk%lastMafIndex)
       IF (l1bFlag==-1) CALL MLSMessage(MLSMSG_Error,ModuleName,&
            & MLSMSG_L1BRead//l1bItemName)

       IF (l1bItem==1) THEN     ! Do something special for time
          DO maf=1,noMAFs
             defaultField(maf)=l1bField%dpField(1,1,maf)+&
                  & (defaultMIFs(maf)-1)*0.166666666666666666666667D0
          END DO
       ELSE                     ! Otherwise this is fairly easy.
          DO maf=1,noMAFs
             defaultField(maf)=l1bField%dpField(1,defaultMIFs(maf),maf)
          END DO
       END IF

       IF (interpolationFactor==1.0) THEN
          interpolatedField=defaultField
       ELSE
          ! Some interpolation is wanted.  I'm going to hold off writing this
          ! because we certaintly don't need it for 0.1 and probably won't till
          ! 1.0.  For the sake of getting things down I'll state here what I
          ! think would be implemented.  One would simply interpolate from the
          ! defaultField to the interpolatedField, using linear or spline I
          ! imagine.  However, there are issues with roll overs for quantities
          ! such as longitude and solarTime.  This is why I have chosen to
          ! defer this piece of code.
          ! NJL - 16 December 1999
          CALL MLSMessage(MLSMSG_Error,ModuleName, &
               & "Sorry interpolation of hGrids is not yet supported")
       ENDIF

       ! Now we have to save this field in the hGrid data.  This is rather a
       ! kludgy way of doing it but this worked out the least boring way to
       ! write the code.  See the definition of L1BItemsToRead above for
       ! reference.

       SELECT CASE(l1bItem)
       CASE(1)
          hGrid%time=interpolatedField
       CASE(2)
          hGrid%geodLat=interpolatedField
       CASE(3)
          hGrid%lon=interpolatedField
       CASE(4)
          hGrid%phi=interpolatedField
       CASE(5)
          hGrid%solarZenith=interpolatedField
       CASE(6)
          hGrid%solarTime=interpolatedField
       CASE(7)
          hGrid%losAngle=interpolatedField
       END SELECT
    END DO

    DEALLOCATE (defaultMIFs)
    DEALLOCATE (defaultField,interpolatedField)

    ! This calculation may need attention! ***
    hGrid%noProfsLowerOverlap=NINT(chunk%noMAFsLowerOverlap*interpolationFactor)
    hGrid%noProfsUpperOverlap=NINT(chunk%noMAFsUpperOverlap*interpolationFactor)
    
  END SUBROUTINE CreateHGridFromMLSCFInfo

  ! ---------------------------------------------------------------------------

  ! This routine destroys the information associated with an hGrid

  SUBROUTINE DestroyHGridContents(hGrid)

    ! Dummy arguments
    TYPE (HGrid_T), INTENT(OUT) :: hGrid

    ! Executable code

    DEALLOCATE(hGrid%phi, hGrid%geodLat, hGrid%lon, hGrid%time, &
         & hGrid%solarTime, hGrid%solarZenith, hGrid%losAngle)

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
       ! Check we don't already have one of this name
       IF (LinearSearchStringArray(database%name,hGrid%name, &
            & caseInsensitive=.TRUE.)/=0) CALL MLSMessage(MLSMSG_Error,&
            & ModuleName,MLSMSG_Duplicate//hGrid%name)
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
! Revision 1.6  2000/01/18 00:14:51  livesey
! Removed profileIndices etc. No longer relevant, as Join deals with this stuff
! for l2gp quantities.
!
! Revision 1.5  2000/01/11 22:51:34  livesey
! Dealt with ramifications of change from read_parse_l2cf to MLSCF
!
! Revision 1.4  2000/01/07 23:53:34  livesey
! Nearly integrated, just a few tweaks.
!
! Revision 1.3  1999/12/17 21:39:34  livesey
! Added check for duplicate name in database.
!
! Revision 1.2  1999/12/16 18:23:20  livesey
! First version that compiles.
!
! Revision 1.1  1999/12/16 01:26:46  livesey
! Nightly checkin
!
!
