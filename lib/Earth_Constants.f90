! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Earth_Constants

  ! This module provides constants related to the Earth's geometry.
  ! They're not in the Geometry module so that others can use them
  ! without dragging along all the procedures therein.

  use MLSKinds, only: R8, RP

  implicit NONE
  public


  ! Earth dimensions.

  real(r8), parameter :: EarthRadA = 6378137.0_r8    ! Semi-Major axis in m
  real(r8), parameter :: EarthRadB = 6356752.3141_r8 ! Semi-Minor axis in m
  real(r8), parameter :: Earth_Axis_Ratio = EarthRadB / EarthRadA
  real(r8), parameter :: Earth_Axis_Ratio_Squared = Earth_Axis_Ratio**2
  real(r8), parameter :: Eccentricity_Sq = &
    & 1.0_r8 - Earth_Axis_Ratio_Squared ! 1 - b**2 / a**2

  ! Gravity-related terms.

  ! IUGG (International Union of Geodesy and Geodynamics) 1999 values from
  ! http://www.gfy.ku.dk/~iag/Travaux_99/sc3.htm (part III)

  ! Geocentric Earth gravitational constant, including Earth's atmosphere

! real(rp), parameter :: GM = 3.986004418e14_rp ! m^3/sec^2 (IUGG 1999 value)
  real(rp), parameter :: GM = 3.98600436e14_rp  ! m^3/sec^2

  ! Mean equatorial gravity

  real(r8), parameter :: G0 = 9.80665           ! Nominal little g ms-2
! real(r8), parameter :: G0 = 9.7803278         ! IUGG 1999 value

  ! Stokes second- and fourth-degree zonal harmonics

  real(rp), parameter :: J2 = 0.0010826256_rp   ! 1980 reference geoid value
! real(rp), parameter :: J2 = 0.0010826359_rp   ! IUGG 1999 value
  real(rp), parameter :: J4 = -.0000023709122_rp  ! 1980 reference geoid value

  ! earth rotational velocity.

  real(rp), parameter :: W = 7.292115e-05_rp    ! rad/sec

  ! Earth surface geopotential height.

  real(r8), parameter :: EarthSurfaceGPH = 6387182.265_r8 ! meters

  ! Seconds per tropical year, 1994-1998, based on orbital elements by
  ! Laskar.  See http://scienceworld.wolfram.com/astronomy/TropicalYear.html

  real(r8), parameter :: SecPerYear = 365.242190_r8 * 86400.0_r8

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Earth_Constants

! $Log$
! Revision 2.1  2016/01/05 03:14:17  vsnyder
! Moved constants from Geometry
!
