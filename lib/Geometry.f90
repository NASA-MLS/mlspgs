! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Geometry

  ! This module contains some geometry routines and constants common to the
  ! forward model and the scan model.

  use MLSKinds, only: R8, RP

  implicit NONE
  private

  ! Constants
  public :: Earth_Axis_Ratio_Squared, Earth_Axis_Ratio_Squared_m1 ! a^2/b^2, a^2/b^2-1
  public :: EarthRadA, EarthRadB, EarthSurfaceGPH
  public :: GM, G0, J2, J4, SecPerYear, W, MaxRefraction

  ! Functions
  public :: GeodToGeocLat, Get_R_Eq

  ! Earth dimensions.

  real(r8), parameter :: EarthRadA = 6378137.0_r8    ! Major axis radius in m 
  real(r8), parameter :: EarthRadB = 6356752.3141_r8 ! Minor axis radius in m 
  real(r8), parameter :: Earth_Axis_Ratio_Squared = EarthRadA**2 / EarthRadB**2
  real(r8), parameter :: Earth_Axis_Ratio_Squared_m1 = &
    & (EarthRadA / EarthRadB - 1) * (EarthRadA / EarthRadB + 1)

  ! Gravity-related terms.

  real(rp), parameter :: GM = 3.98600436e14_rp ! m^3/sec^2
  real(r8), parameter :: G0 = 9.80665          ! Nominal little g ms-2

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

  ! This is the maximum amount of refraction allowed
  real(r8), parameter :: MaxRefraction = 0.0004 ! Add one to get refractive index

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------  GeodToGeocLat  -----
  real(r8) elemental function GeodToGeocLat ( geodLat )

  ! Convert a geodetic latitude (IN DEGREES!) into a geocentric one (IN RADIANS!)
  ! (IN RADIANS!)

    use Constants, only: Deg2Rad, PI

    ! Arguments
    real (r8), intent(in) :: geodLat

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

  ! ----------------------------------------------------  Get_R_eq -----
  real(rp) elemental function Get_R_Eq ( Phi, Csq ) result ( R_eq )

  !{ Given the orbit geodetic longitude {\tt Phi} = $\phi$ in radians and the
  !  square of the minor axis of the orbit plane projected Earth ellipse in
  !  meters {\tt Csq} = $R_c^2$ compute the radius in kilometers of an
  !  equivalent circular Earth tangent to the elliptical Earth and having the
  !  same radius of curvature as the elliptical Earth at $\phi$.
  !%
  ! \begin{equation*}
  ! R_{eq} = \sqrt \frac{R_a^4 \sin^2 \phi + R_c^4 \cos^2 \phi}
  !                    {R_a^2 \cos^2 \phi + R_c^2 \sin^2 \phi}
  !        = \sqrt \frac{R_a^4 - (R_a^2+R_c^2)(R_a^2-R_c^2) \cos^2 \phi}
  !                       {R_c^2 +              (R_a^2-R_c^2) \cos^2 \phi}
  ! \end{equation*}
  !%
  ! This is Equation (5.21) in the 19 August 2004 ATBD JPL D-18130.

  real(rp), intent(in) :: Phi
  real(rp), intent(in) :: Csq

  real(rp), parameter :: Earthrada_sq = earthrada ** 2
  real(rp), parameter :: Earthrada_4 = earthrada_sq ** 2

  r_eq = (earthrada_sq - csq) * COS(phi)**2
  ! Earthrad[abc] are in meters, but r_eq needs to be in km.
  r_eq = 0.001_rp * SQRT( &
    & ( earthrada_4 -(earthrada_sq + csq) * r_eq ) / &
    & ( csq + r_eq ) )

  end function Get_R_Eq

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Geometry

! $Log$
! Revision 2.19  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.18  2009/05/13 20:13:58  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.17  2008/10/08 01:11:37  vsnyder
! Add PRINT statement in not_used_here to prevent compilation cascades
!
! Revision 2.16  2008/10/08 01:07:59  vsnyder
! Add Get_R_Eq function
!
! Revision 2.15  2006/09/28 20:51:46  vsnyder
! Add Earth_Axis_Ratio_Squared_m1
!
! Revision 2.14  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.13  2003/01/16 23:13:10  livesey
! Added maxRefraction
!
! Revision 2.12  2003/01/15 02:45:50  vsnyder
! Make SecPerYear public (oops!)
!
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
