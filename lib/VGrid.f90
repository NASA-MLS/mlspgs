! Copyright (c) 1999, California Institute of Technology. ALL RIGHTS RESERVED.
! U.S. Government sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE vGrid                    ! Definitions for vGrids in vector quantities
!=============================================================================

  USE MLSCommon                 ! General constants etc.
  USE L2CF                      ! Information on l2cf data
  USE VerticalCoordinates       ! The various vertical coorindate systems.
  USE MLSStrings                ! String handling routines

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: ID, ModuleName
  !------------------------------- RCS Ident Info ---------------------------
  CHARACTER (LEN=130) :: Id= "$Id$"
  !--------------------------------------------------------------------------

  ! Define the vGrid data type.  This is used to store all the vGrid
  ! information. Note that this is only relevant for coherent quantities. 
  ! Incoherent ones deal with vGrids seperately.

  TYPE vGrid_T
     CHARACTER (LEN=NameLen) :: name ! Name for vGrid
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

    ! Local variables

    INTEGER :: keyNo            ! Entry in the l2cf line

    ! Executable code

    ! We will go through the information given in the l2cf and create an
    ! appropriate vGrid for it.

    DO keyNo=1,cfInfo%l2cfEntryNoKeys
       SELECT 
       
    END DO

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

END MODULE vGrid

!
! $Log$
! Revision 1.1  1999/11/24 23:06:19  livesey
! First simple version
!
!
