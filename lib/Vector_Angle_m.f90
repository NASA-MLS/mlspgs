! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Vector_Angle_m

  use Constants, only: Pi
  use Cross_m, only: Cross

  implicit NONE
  private
  public :: Vector_Angle, Vector_Angle_D, Vector_Angle_S

  !{ Compute the angle between two vectors.
  !
  ! The usual formula is
  ! \begin{equation*}
  !  \theta = \cos^{-1} z \text{ where }
  !   z = \frac{x \cdot y}{ || x || \, || y || } \,.
  ! \end{equation*}
  ! This can become severely inaccurate if $x$ and $y$ are nearly parallel
  ! or antiparallel because $\text{d} \theta / \text{d} z$ has a singularity
  ! at $z = 1$.
  !
  ! An alternative is
  ! \begin{equation*}
  !  \theta = \sin^{-1} \left( \frac { || x \times y || }
  !                                  { || x || \, || y || }
  !                     \right) \,\,\, ( x \cdot y \geq 0 ) \text{ or }
  !  \theta = \pi -\sin^{-1} \left( \frac { || x \times y || }
  !                               { || x || \, || y || }
  !                     \right) \,\,\, ( x \cdot y < 0 ) \,,
  ! \end{equation*}
  ! but this only works for three-dimensional vectors, and is inaccurate
  ! when $x$ and $y$ are nearly orthogonal because $\sin^{-1}$ has a
  ! derivative singularity at $\pm 90^\circ$.
  !
  ! A reliable but expensive formula is
  ! \begin{equation*}
  !  \theta = 2 \tan^{-1} \left(
  !   \frac{ || ( x || y || - y || x || ) || }
  !        { || ( x || y || + y || x || ) || } \right)
  ! \end{equation*}
  !
  ! See wvs-154.

  interface Vector_Angle
    module procedure Vector_Angle_D, Vector_Angle_S
  end interface Vector_Angle

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  double precision function Vector_Angle_D ( X, Y ) result ( Theta )
    double precision, intent(in) :: X(:), Y(:)
    double precision :: XN, YN, Z
    xn = norm2(x)
    yn = norm2(y)
    z = dot_product ( x, y ) / ( xn * yn )
    if ( abs(z) <= 0.9 ) then
      theta = acos(z)
    else if ( size(x) == 3 ) then
      theta = asin ( norm2(cross(x,y)) / ( xn * yn ) )
      if ( z < 0 ) theta = pi - theta
    else
      theta = atan2 ( norm2(x * yn - y * xn), norm2(x * yn + y * xn) )
    end if
  end function Vector_Angle_D

  real function Vector_Angle_S ( X, Y ) result ( Theta )
    real, intent(in) :: X(:), Y(:)
    real :: XN, YN, Z
    xn = norm2(x)
    yn = norm2(y)
    z = dot_product ( x, y ) / ( xn * yn )
    if ( abs(z) <= 0.9 ) then
      theta = acos(z)
    else if ( size(x) == 3 ) then
      theta = asin ( norm2(cross(x,y)) / ( xn * yn ) )
      if ( z < 0 ) theta = pi - theta
    else
      theta = atan2 ( norm2(x * yn - y * xn), norm2(x * yn + y * xn) )
    end if
  end function Vector_Angle_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Vector_Angle_m

! $Log$
! Revision 2.1  2019/09/17 00:53:36  vsnyder
! Initial commit
!
