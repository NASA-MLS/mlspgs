! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Radius_of_Curvature_m

  !{ Compute the radius of curvature $R_N$ in the \emph{prime vertical}
  !  direction, at a geodetic latitude {\tt Lat} denoted by $\phi$.
  !  \begin{equation}
  !  R_N = \frac{a}{\sqrt{1-e^2 \sin^2 \phi}}
  !  \end{equation}
  !  where $a$ is the Earth's equatorial radius, $b$ is the polar radius, and
  !  $e = \sqrt{1-\frac{b^2}{a^2}}$ is the eccentricity.
  !
  !  Compute the radius of curvature $R_M$ in the plane of the meridian.
  !  \begin{equation}
  !  R_M = \frac{a (1 - e^2)}{\sqrt{(1 - e^2 \sin^2(\phi))^3}} =
  !  \frac{b^2}{a\sqrt{(1 - e^2 \sin^2(\phi))^3}} =
  !  R_N \frac{1-e^2}{1-e^2 \sin^t \phi}
  !  \end{equation}
  !
  !  The plane of the meridian and the plane in which the radius in the
  !  prime vertical direction is computed are normal to each other, and both
  !  pass through the center of the Earth.
  !
  !  Compute the mean radius of curvature $R_A$ at geodetic latitude $\phi$.
  !  \begin{equation}
  !  R_A = \frac2{\frac1{R_M} + \frac1{R_N}}
  !  \end{equation}

  implicit NONE
  private

  public :: Radius_of_Curvature_Mean
  public :: Radius_of_Curvature_Meridian
  public :: Radius_of_Curvature_Normal

  interface Radius_of_Curvature_Mean
    module procedure Radius_of_Curvature_Mean_D
    module procedure Radius_of_Curvature_Mean_S
  end interface

  interface Radius_of_Curvature_Meridian
    module procedure Radius_of_Curvature_Meridian_D
    module procedure Radius_of_Curvature_Meridian_S
  end interface

  interface Radius_of_Curvature_Normal
    module procedure Radius_of_Curvature_Normal_D
    module procedure Radius_of_Curvature_Normal_S
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  pure elemental double precision function Radius_of_Curvature_Mean_D ( Lat ) result ( R_A )
    double precision, intent(in) :: Lat ! Geodetic latitude in degrees

    R_A = 2.0d0 / ( 1.0d0 / Radius_of_Curvature_Meridian ( lat ) + &
                  & 1.0d0 / Radius_of_Curvature_Normal ( lat ) )

  end function Radius_of_Curvature_Mean_D

  pure elemental real function Radius_of_Curvature_Mean_S ( Lat ) result ( R_A )
    real, intent(in) :: Lat ! Geodetic latitude in degrees

    R_A = 2.0e0 / ( 1.0e0 / Radius_of_Curvature_Meridian ( lat ) + &
                  & 1.0e0 / Radius_of_Curvature_Normal ( lat ) )

  end function Radius_of_Curvature_Mean_S

  pure elemental double precision function Radius_of_Curvature_Meridian_D ( Lat ) result ( R_M )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq, &
      & R2 => Earth_Axis_Ratio_Squared ! b**2 / a**2 = 1 - E2
    double precision, intent(in) :: Lat ! Geodetic latitude in degrees

    R_M = ( a * r2 ) / sqrt( 1 - e2 * sin(lat*deg2rad)**2 ) ** 3

  end function Radius_of_Curvature_Meridian_D

  pure elemental real function Radius_of_Curvature_Meridian_S ( Lat ) result ( R_M )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq, &
      & R2 => Earth_Axis_Ratio_Squared ! b**2 / a**2 = 1 - E2
    real, intent(in) :: Lat ! Geodetic latitude in degrees

    R_M = ( a * r2 ) / sqrt( 1 - e2 * sin(lat*deg2rad)**2 ) ** 3

  end function Radius_of_Curvature_Meridian_S

  pure elemental double precision function Radius_of_Curvature_Normal_D ( Lat ) result ( R_N )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq
    double precision, intent(in) :: Lat ! Geodetic latitude in degrees

    R_N = a / sqrt( 1 - e2 * sin(lat*deg2rad)**2 )

  end function Radius_of_Curvature_Normal_D

  pure elemental real function Radius_of_Curvature_Normal_S ( Lat ) result ( R_N )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq
    real, intent(in) :: Lat ! Geodetic latitude in degrees

    R_N = a / sqrt( 1 - e2 * sin(lat*deg2rad)**2 )

  end function Radius_of_Curvature_Normal_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Radius_of_Curvature_m

! $Log$
! Revision 2.2  2016/03/08 02:52:32  vsnyder
! Change name of Radius_Of_Curvature to Radius_Of_Curvature_Meridian.  Add
! Radius_Of_Curvature_Normal and Radius_Of_Curvature_Mean.
!
! Revision 2.1  2016/02/24 01:15:37  vsnyder
! Initial commit
!
