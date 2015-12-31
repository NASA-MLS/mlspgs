! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Geolocation_0
!=============================================================================

  ! Define types used in geolocation modules.

  use MLSKinds, only: RG => R8 ! Kind of REAL geolocation variables

  implicit NONE

  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

  type :: H_t           ! For horizontal (lat/lon) grids, geoc/geod unspecified
    real(rg) :: Lon     ! Degrees
    real(rg) :: Lat     ! Degrees
  end type H_t 

  type, extends(h_t) :: H_Geoc
  !  real(rg) :: Lon    ! Longitude, degrees
  !  real(rg) :: Lat    ! Geocentric latitude, degrees
  end type H_Geoc

  type, extends(h_t) :: H_Geod
  !  real(rg) :: Lon    ! Longitude, degrees
  !  real(rg) :: Lat    ! Geodetic latitude, degrees
  end type H_Geod

  type, abstract :: Lat_t ! For latitude; can't instantiate it
    real(rg) :: D       ! Degrees
  end type Lat_t

  type, extends(lat_t) :: GeocLat_t
  ! real(rg) :: D       ! Geocentric latitude, degrees
  end type GeocLat_t

  type, extends(lat_t) :: GeodLat_t
  ! real(rg) :: D       ! Geodetic latitude, degrees
  end type GeodLat_t

  type :: Lon_t         ! Longitude, degrees
    real(rg) :: D       ! Degrees
  end type Lon_t

  type :: V_t           ! For vertical grids of unspecified type
    real(rg) :: V       ! Geocentric (meters), geodetic (meters), or zeta
  end type V_t

  type, extends(v_t) :: V_geoc
  ! real(rg) :: V       ! geocentric height (meters from the Earth center)
  end type V_geoc

  type, extends(v_t) :: V_geod
  ! real(rg) :: V       ! geodetic height (meters above sea level)
  end type V_geod

  type, extends(v_t) :: V_zeta
  ! real(rg) :: V       ! Zeta -- geocentric latitude doesn't make sense if
                        ! the coordinate system is stacked
  end type V_zeta

  ! For irregular grids, or lines of sight, where H and V are expected to
  ! have the same extents
  type, extends(h_t) :: H_V_t
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
    real(rg) :: V       ! Meters or zeta
  end type H_V_t

  type, extends(h_v_t) :: H_V_Geoc
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
  ! real(rg) :: V       ! geocentric height (meters from the Earth center)
  end type H_V_Geoc

  type, extends(h_v_t) :: H_V_Geod
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
  ! real(rg) :: V       ! geodetic height (meters above sea level)
  end type H_V_Geod

  type, extends(h_v_t) :: H_V_Zeta
  ! real(rg) :: Lon     ! Degrees
  ! real(rg) :: Lat     ! Degrees
  ! real(rg) :: V       ! Zeta
  end type H_V_Zeta

  type :: ECR_t
    real(rg) :: XYZ(3)  ! Cartesian ECR, meters
  end type ECR_t

contains

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Geolocation_0

! $Log$
! Revision 2.3  2015/12/31 00:04:14  vsnyder
! Make V_t and H_V_t nonabstract
!
! Revision 2.2  2015/12/01 21:05:03  vsnyder
! Make H_t non-abstract
!
! Revision 2.1  2015/11/14 18:02:36  vsnyder
! Initial commit
!
