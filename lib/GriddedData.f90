

! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE GriddedData ! Contains the derived TYPE GriddedData_T
!=============================================================================

  USE MLSCommon, only: R8, LineLen, NameLen

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id,ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  public::GriddedData_T

  public::v_is_pressure,v_is_altitude,v_is_gph,v_is_theta
	! These are 'enumerated types' consistent with hph's
	! work in l3ascii_read_field
	
	integer, parameter :: v_is_pressure = 1
	integer, parameter :: v_is_altitude = v_is_pressure+1
	integer, parameter :: v_is_gph = v_is_altitude+1
	integer, parameter :: v_is_theta = v_is_gph+1
	
!
! This type reflects the format of the Level 3 ASCII files, though note that
! these files can store multiple quantities such as these.
!
TYPE GriddedData_T
   !
   ! First the comment line(s) from the relevant input file
   !
   CHARACTER (LEN=LineLen), POINTER, DIMENSION(:) :: fileComments => NULL()
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
   REAL (R8), POINTER, DIMENSION(:) :: heights  => NULL()
             ! Surfaces (e.g. pressures etc.) [noHeights]
   !
   ! Now the latitudinal coordinate
   !
   LOGICAL :: equivalentLatitude ! If set, coordinate is equivalent latitude
   INTEGER :: noLats            ! Number of latitudes
   REAL (R8), POINTER, DIMENSION(:) :: lats => NULL() ! Latitudes [noLats]
   !
   INTEGER :: noLons            ! Number of longitudes
   REAL (R8), POINTER, DIMENSION(:) :: lons => NULL() ! Longitudes [noLons]
   !
   INTEGER noLsts               ! Number of local times
   REAL (R8), POINTER, DIMENSION(:) :: lsts => NULL() ! Local times [noLsts]
   !
   INTEGER noSzas               ! Number of solar zenith angles
   REAL (R8), POINTER, DIMENSION(:) :: szas => NULL() ! Zenith angles [noSzas]
   !
   INTEGER noDates              ! Number of dates in data
   REAL (R8), POINTER, DIMENSION(:) :: dateStarts => NULL()
      ! Starting dates in SDP toolkit format
   REAL (R8), POINTER, DIMENSION(:) :: dateEnds => NULL()
      ! Ending dates in SDP toolkit format
   !
   REAL (R8), POINTER, DIMENSION(:,:,:,:,:,:) :: field => NULL()
   !
   ! The data itself.  This is stored as
   !  [noHeights, noLats, noLons, noLsts, noSzas, noDates]
   !
END TYPE GriddedData_T


!=============================================================================
END MODULE GriddedData
!=============================================================================

!
! $Log$
! Revision 2.9  2001/03/15 21:28:08  pwagner
! Now only the data type
!
!
