

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
   CHARACTER (LEN=NameLen) :: sourceFileName ! Input file name
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
! Revision 2.10  2001/03/30 00:24:40  pwagner
! Added sourceFileName
!
! Revision 2.9  2001/03/15 21:28:08  pwagner
! Now only the data type
!
! Revision 2.8  2001/03/15 00:37:13  pwagner
! Still not complete; missing gdrdfld
!
! Revision 2.7  2001/03/14 00:32:47  pwagner
! More changes--still wrong, though
!
! Revision 2.6  2001/03/10 00:33:16  pwagner
! Some corrections in ReadGriddedData
!
! Revision 2.5  2001/03/09 01:02:55  pwagner
! Fixed announce_error
!
! Revision 2.4  2001/03/08 01:08:35  pwagner
! Added announce_error
!
! Revision 2.3  2001/03/07 01:03:19  pwagner
! ReadGriddedData added
!
! Revision 2.2  2001/02/21 00:36:43  pwagner
! l3ascii_read_field now has eof as intent(out) arg
!
! Revision 2.1  2001/02/20 21:51:39  pwagner
! Functions absorbed from gridded_data_module
!
! Revision 2.0  2000/09/05 17:41:05  dcuddy
! Change revision to 2.0
!
! Revision 1.5  2000/06/20 22:19:22  lungu
! Changed DOUBBLE PRECISION to REAL (r8).
!
