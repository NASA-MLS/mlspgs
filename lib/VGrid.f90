! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE vGrid                    ! Definitions for vGrids in vector quantities
!=============================================================================

  USE MLSCommon                 ! General constants etc.
  USE Temporary_Types           ! L2CF info etc.
  USE VerticalCoordinate        ! The various vertical coorindate systems.
  USE MLSStrings                ! String handling routines

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: ID, ModuleName
  !------------------------------- RCS Ident Info ---------------------------
  CHARACTER (LEN=130) :: Id= "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !--------------------------------------------------------------------------

  ! Define the vGrid data type.  This is used to store all the vGrid
  ! information. Note that this is only relevant for coherent quantities. 
  ! Incoherent ones deal with vGrids seperately.

  INTEGER, PARAMETER :: vGridNameLen=32
  TYPE vGrid_T
     CHARACTER (LEN=vGridNameLen) :: name ! Name for vGrid
     INTEGER :: verticalCoordinate ! Enumerated type e.g. VC_Pressure
     INTEGER :: noSurfs         ! Number of surfaces
     REAL(r8), DIMENSION(:), POINTER :: surfs  ! Array of surfaces
     ! (actually dimensioned noSurfs)
  END TYPE vGrid_T

CONTAINS

  !--------------------------------------------------------------------------

  ! This routine creates a vGrid according to user supplied information in the
  ! l2cf.

  SUBROUTINE CreateVGridFromMLSCFInfo(vGrid, cfInfo)

    ! Dummy arguments
    TYPE(vGrid_T),    INTENT(OUT) :: vGrid ! Returned vGrid
    TYPE(L2CfEntry), INTENT(IN)  :: cfInfo ! Input info. from l2cf

    ! Local parameters
    CHARACTER (LEN=*), PARAMETER :: UnitsMessage= &
         & "Inappropriate units for vertical coordinates"

    ! Local variables
    INTEGER :: keyNo            ! Entry in the l2cf line
    INTEGER :: family           ! Physical quantity family
    TYPE (L2CFCell) :: cell     ! Part of the l2cf information

    ! Executable code

    ! We will go through the information given in the l2cf and create an
    ! appropriate vGrid for it.

    vGrid%name=""
    vGrid%verticalCoordinate=VC_Invalid
    vGrid%noSurfs=0
    IF (ASSOCIATED(vGrid%surfs)) CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "vGrid%surfs already associated")

    DO keyNo=1,cfInfo%l2cfEntryNoKeys
       cell=cfInfo%cells(keyNo)
       SELECT CASE(TRIM(cell%keyword))
       CASE ("NAME")
          vGrid%name=cell%charValue
       CASE ("COORDINATE")
          vGrid%verticalCoordinate=ParseVerticalCoordinateName(cell%charValue)
       CASE ("VALUES")
          CALL ParseVerticalCoordinateSpec(cell%charValue,vGrid%surfs,family)
       CASE DEFAULT
          CALL MLSMessage(MLSMSG_Error,ModuleName,"Unexpected key "//&
               & TRIM(cell%keyword))
       END SELECT
    END DO

    ! Now check that this is a sensible vGrid, first the obvious stuff

    IF (LEN_TRIM(vGrid%name)=="") CALL MLSMessage(MLSMSG_Error,ModuleName,&
         & "Invalid or absent vGrid name")
    IF (vGrid%verticalCoordinate==VC_Invalid) CALL MLSMessage(MLSMSG_Error,&
         & ModuleName,"Invalid vertical coordinate for vGrid "//&
         & TRIM(vGrid%name))
    IF (.NOT.(ASSOCIATED(vGrid%surfs))) CALL MLS_Message(MLSMSG_Error,&
         & ModuleName,"Invalid/absent surfaces in vGrid "//TRIM(vGrid%name))

    ! Check that the given surfaces are in an appropriate unit

    SELECT CASE(vGrid%verticalCoordinate)
    CASE (VC_None)
       IF (family/=PHYQ_Dimensionless) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_Pressure)
       IF (family/=PHYQ_Pressure) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_Zeta)
       IF (family/=PHYQ_Dimensionless) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_Altitude)
       IF (family/=PHYQ_Length) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_GPH)
       IF (family/=PHYQ_Length) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_Theta)
       IF (family/=PHYQ_Temperature) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    CASE (VC_Angle)
       IF (family/=PHYQ_Angle) &
            & CALL MLSMessage(MLSMSG_Error,ModuleName,UnitsMessage)
    END SELECT

  END SUBROUTINE CreateVGridFromMLSCFInfo

  !--------------------------------------------------------------------------

  ! This routine destroys the array information created with the vGrid

  SUBROUTINE DestroyVGridInformation(vGrid)

    ! Dummy arguments

    TYPE (vGrid_T), INTENT(INOUT) :: vGrid

    ! Executable code

    vGrid%name=''
    vGrid%noSurfs=0
    vGrid%verticalCoordinate=VC_Invalid

    DEALLOCATE(vGrid%surfs)

  END SUBROUTINE DestroyVGridInformation

  !--------------------------------------------------------------------------

  ! This routine adds a vGrid to a database of vGrids, creating the database
  ! if necessary.

  SUBROUTINE AddVGridToDatabase(database,vgrid)

    ! Dummy arguments
    TYPE (VGrid_T), DIMENSION(:), POINTER :: database
    TYPE (VGrid_T), INTENT(IN) :: vgrid

    ! Local variables
    TYPE (VGrid_T), DIMENSION(:), POINTER :: tempDatabase
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
    tempDatabase(newSize)=vgrid
    DEALLOCATE(database)
    database=>tempDatabase
  END SUBROUTINE AddVGridToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a vGrid database

  SUBROUTINE DestroyVGridDatabase(database)

    ! Dummy argument
    TYPE (VGrid_T), DIMENSION(:), POINTER :: database

    ! Local variables
    INTEGER :: vgridIndex

    IF (ASSOCIATED(database)) THEN
       DO vgridIndex=1,SIZE(database)
          CALL DestroyVGridContents(database(vgridIndex))
       ENDDO
       DEALLOCATE(database)
    ENDIF
  END SUBROUTINE DestroyVGridDatabase  

END MODULE vGrid

!
! $Log$
! Revision 1.2  1999/12/14 22:56:02  livesey
! Regular commit
!
! Revision 1.1  1999/11/24 23:06:19  livesey
! First simple version
!
!
