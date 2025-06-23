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
  !  $e = \sqrt{1-\frac{b^2}{a^2}}$ is the eccentricity of the Earth oblate
  !  spheroid.
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
  !
  !  Compute the radius of curvature $R_\text{eq}^\oplus$ in the normal
  !  section containing the line of sight, the sub-tangent point at the
  !  Earth surface, and the gradient at the sub-tangent point.
  !  See wvs-146.
  !  \begin{equation}
  !  R^\oplus_{\text{eq}}
  !   = a' \, (1-\epsilon^2)
  !     \left( \frac{A^2+B^2}{A^2+(1-\epsilon^2)B^2} \right)^\frac32
  !   = F^2 \sqrt{a^2 W^2 -d^2}
  !      \left( \frac{A^2+B^2}{A^2 W^2 + B^2 F^2} \right) ^\frac32\,,
  !  \end{equation}
  !  where
  !  \begin{equation}
  !    a'^2 = a^2 - \frac{d^2}{W^2}\,,\,
  !    \epsilon^2 = 1 - \frac{F^2}{W^2} = \frac{A^2+B^2}{W^2} \,,
  !  \end{equation}
  !  $\vec{G}$ is the outward gradient at the sub-tangent point,
  !  $\vec{V}$ is a vector parallel to the line of sight from the instrument
  !  to the tangent point,
  !  $\vec{N} = \vec{G} \times \vec{V}$ is normal to the plane containing
  !  the line of sight, the sub-tangent point, and the gradient at the
  !  sub-tangent point, $\hat{n} = \frac{\vec{N}}{|\vec{N}|} = [ A, B, C ]^T$
  !  is the unit vector parallel to $\vec{N}$, $d = \hat{n} \cdot \vec{T}$
  !  is the distance from the center of the Earth to the plane normal to
  !  $\hat{n}$, $F^2 = (1-e^2)\, C^2$, and $W^2 = A^2 + B^2 + F^2$.

  implicit NONE
  private

  public :: Radius_of_Curvature_Mean
  public :: Radius_of_Curvature_Meridian
  public :: Radius_of_Curvature_Normal
  public :: Radius_of_Curvature_Prime_Vertical
  
  interface Radius_of_Curvature_Mean
    module procedure Radius_of_Curvature_Mean_D
    module procedure Radius_of_Curvature_Mean_S
  end interface

  interface Radius_of_Curvature_Meridian
    module procedure Radius_of_Curvature_Meridian_D
    module procedure Radius_of_Curvature_Meridian_S
  end interface

  interface Radius_of_Curvature_Prime_Vertical
    module procedure Radius_of_Curvature_Prime_Vertical_D
    module procedure Radius_of_Curvature_Prime_Vertical_S
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
                  & 1.0d0 / Radius_of_Curvature_Prime_Vertical ( lat ) )

  end function Radius_of_Curvature_Mean_D

  pure elemental real function Radius_of_Curvature_Mean_S ( Lat ) result ( R_A )
    real, intent(in) :: Lat ! Geodetic latitude in degrees

    R_A = 2.0e0 / ( 1.0e0 / Radius_of_Curvature_Meridian ( lat ) + &
                  & 1.0e0 / Radius_of_Curvature_Prime_Vertical ( lat ) )

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

  pure elemental double precision function Radius_of_Curvature_Prime_Vertical_D ( Lat ) result ( R_N )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq
    double precision, intent(in) :: Lat ! Geodetic latitude in degrees

    R_N = a / sqrt( 1 - e2 * sin(lat*deg2rad)**2 )

  end function Radius_of_Curvature_Prime_Vertical_D

  pure elemental real function Radius_of_Curvature_Prime_Vertical_S ( Lat ) result ( R_N )
    use Constants, only : Deg2Rad
    use Earth_Constants, only : A => EarthRadA, E2 => Eccentricity_Sq
    real, intent(in) :: Lat ! Geodetic latitude in degrees

    R_N = a / sqrt( 1 - e2 * sin(lat*deg2rad)**2 )

  end function Radius_of_Curvature_Prime_Vertical_S

  ! Radius of curvature in the normal section from the instrument to the
  ! tangent point
  pure elemental real(rg) function &
    & Radius_of_Curvature_Normal ( V, T ) result ( R_eq )
    use Earth_Constants, only: EarthRadA, R2 => Earth_Axis_Ratio_Squared
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: V, T
    double precision :: A, B, D, F, W ! A^2, B^2, d^2, F^2, W^2

    type(ECR_t) :: G, N ! Gradient, Unit normal

    g = ecr_t ( [ t%xyz(1), t%xyz(2), t%xyz(3) / r2 ] ) ! X, Y, Z/(1-e^2)
    n = g .crossnorm. v     ! unit normal to the plane of the normal section
    d = ( n .dot. t ) ** 2  ! d**2
    a = n%xyz(1) ** 2       ! A**2
    b = n%xyz(2) ** 2       ! B**2
    f = n%xyz(3) ** 2 / r2  ! F**2 = C**2 / ( 1-e^2 ) = C**2 / R2
    w = a + b + f
    R_eq = f * sqrt ( ( earthRadA**2 * w - d ) * &
                    & ( ( a + b ) / ( a * w + b * f ) )**3 )

  end function Radius_of_Curvature_Normal

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
! Revision 2.3  2020/01/28 00:46:21  vsnyder
! Change name of Radius_of_Curvature_Normal to Radius_of_Curvature_Prime_Vertical.
! Add new Radius_of_Curvature_Normal for normal sections that are not prime
! vertical sections, i.e., not east-west sections.
!
! Revision 2.2  2016/03/08 02:52:32  vsnyder
! Change name of Radius_Of_Curvature to Radius_Of_Curvature_Meridian.  Add
! Radius_Of_Curvature_Normal and Radius_Of_Curvature_Mean.
!
! Revision 2.1  2016/02/24 01:15:37  vsnyder
! Initial commit
!
