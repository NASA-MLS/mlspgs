! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Geometry

  ! This module contains some geometry routines and constants common to the
  ! forward model and the scan model.

  use MLSCommon, only: R8, RP

  implicit NONE
  private

  public :: Earth_Axis_Ratio_Squared, EarthRadA, EarthRadB, EarthSurfaceGPH
  public :: GM, G0, J2, J4, W

  public :: GeodToGeocLat

  ! Earth dimensions.

  real(r8), parameter :: EarthRadA = 6378137.0_r8    ! Major axis radius in m 
  real(r8), parameter :: EarthRadB = 6356752.3141_r8 ! Minor axis radius in m 
  real(r8), parameter :: Earth_Axis_Ratio_Squared = EarthRadA**2 / EarthRadB**2

  ! Gravity-related terms.

  real(r8), parameter :: G0 = 9.80665          ! Nominal little g ms-2
  real(rp), parameter :: GM = 3.98600436e14_rp ! m^3/sec^2

  ! These are the 1980 reference geoid values.

  real(rp), parameter :: J2 = 0.0010826256_rp
  real(rp), parameter :: J4 = -.0000023709122_rp

  ! earth rotational velocity.

  real(rp), parameter :: W = 7.292115e-05_rp   ! rad/sec

  ! Earth surface geopotential height.

  real(r8), parameter :: EarthSurfaceGPH = 6387182.265_r8 ! meters

  ! Seconds per tropical year, 1994-1998, based on orbital elements by
  ! Laskar.  See http://scienceworld.wolfram.com/astronomy/TropicalYear.html

  real(r8), parameter :: SecPerYear = 365.242190_r8 * 86400.0_r8

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ------------------------------- Subroutines and functions ----

  ! This function converts a geodetic latitude (IN DEGREES!) into a geocentric
  ! one (IN RADIANS!)
  
  real(r8) elemental function GeodToGeocLat ( geodLat )

    use Units, only: Deg2Rad, PI

    ! Arguments
    real (r8), intent(IN) :: geodLat

    ! Executable code, use special method for high latitudes.
    if ( geodLat > 89.0 ) then
      geodtoGeocLat = 0.5 * PI + Earth_Axis_Ratio_Squared * &
           (geodLat-90.0) * deg2rad
    else if ( geodLat < -89.0 ) then
      geodToGeocLat = -0.5 * PI + Earth_Axis_Ratio_Squared * &
           (geodLat+90.0) * deg2rad
    else
       geodToGeocLat=atan( (1.0/Earth_Axis_Ratio_Squared) * tan(geodLat*deg2rad))
    end if
  end function GeodToGeocLat

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Geometry

! $Log$
! Revision 2.11  2003/01/15 02:35:08  vsnyder
! Add SecPerYear, move a USE to procedure scope
!
! Revision 2.10  2003/01/10 21:55:12  vsnyder
! Move SpeedOfLight from Geometry ot Units
!
! Revision 2.9  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/10/02 21:05:39  vsnyder
! Add SpeedOfLight constant
!
! Revision 2.7  2002/09/26 20:56:53  vsnyder
! Move Earth_Axis_Ratio_Squared from heights_module to Geometry
!
! Revision 2.6  2002/09/26 20:33:19  vsnyder
! Remove PI and LN10 -- they're in Units
!
! Revision 2.5  2002/09/26 16:27:14  livesey
! Changes from Van, new constants etc.
!
! Revision 2.4  2001/03/28 19:50:19  vsnyder
! Change constants from d0 to _r8
!
! Revision 2.3  2001/03/27 18:09:11  vsnyder
! Revised CVS stuff to use CHARACTER(len=*) parameter
!
! Revision 2.2  2001/03/26 13:23:34  livesey
! Added CVS stuff
!
