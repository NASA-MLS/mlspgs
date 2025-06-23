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

  !{ Compute the angle in radians between two vectors.
  !
  ! The usual formula is
  ! \begin{equation}\label{one}
  !  \theta = \cos^{-1} z \text{ where }
  !   z = \frac{x \cdot y}{ || x || \, || y || } \,.
  ! \end{equation}
  ! This can become severely inaccurate if $x$ and $y$ are nearly parallel
  ! or antiparallel because $\text{d} \theta / \text{d} z$ has a singularity
  ! at $z = 1$.
  !
  ! An alternative is
  ! \begin{equation}\label{two}
  !  \theta = \sin^{-1} \left( \frac { || x \times y || }
  !                                  { || x || \, || y || }
  !                     \right) \,\,\, ( x \cdot y \geq 0 ) \text{ or }
  !  \theta = \pi -\sin^{-1} \left( \frac { || x \times y || }
  !                               { || x || \, || y || }
  !                     \right) \,\,\, ( x \cdot y < 0 ) \,,
  ! \end{equation}
  ! but this only works for three-dimensional vectors, and is inaccurate
  ! when $x$ and $y$ are nearly orthogonal because $\sin^{-1}$ has a
  ! derivative singularity at $\pm 90^\circ$.
  !
  ! A reliable but expensive formula is
  ! \begin{equation}\label{three}
  !  \theta = 2 \tan^{-1} \left(
  !   \frac{ || ( x || y || - y || x || ) || }
  !        { || ( x || y || + y || x || ) || } \right)
  ! \end{equation}
  !
  ! See wvs-154.
  !
  ! The errors of methods (\ref{one}) and (\ref{two}) are equal at
  ! $z = \frac12 \sqrt{2}$, where the errors are
  ! $\sqrt{2}\, \epsilon \approx 0.7071 \, \epsilon$.
  !
  ! The errors of methods (\ref{one}) and (\ref{three}) are equal at
  ! $z = \sqrt{3 - 2 \sqrt{3}}$, where the errors are
  ! $\frac1{\sqrt{4-2\sqrt{3}}} \, \epsilon \approx 1.366 \, \epsilon$.

  interface Vector_Angle
    module procedure Vector_Angle_D, Vector_Angle_S
  end interface Vector_Angle

  real(kind(pi)), parameter :: HalfPi = 0.5d0 * Pi

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!                                                        !!!!!
  !!!!! ASSUMES SIZE(X) == SIZE(Y)!                            !!!!!
  !!!!!                                                        !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision function Vector_Angle_D ( X, Y ) result ( Theta )
    double precision, intent(in) :: X(:), Y(:)
    integer :: N
    double precision, parameter :: HalfPi = 0.5d0 * Pi
    ! S1 ~ 0.7071 is the value at which arccos and arcsin methods have
    ! equal error.
    double precision, parameter :: S1 = 0.5d0 * sqrt(2.0d0)
    ! S2 ~ 0.6813 is the value at which arccos and arctan methods have
    ! equal error. The difference in the errors of the arccos method at
    ! S1 and S2 ~ 0.048 epsilon(1.0d0).
    double precision, parameter :: S2 = sqrt ( 2.0d0 * sqrt(3.0d0) - 3.0d0 )
    double precision :: XN, YN, Z
    n = size(x)
    xn = norm2(x)
    yn = norm2(y)
    if ( xn == 0.0 .or. yn == 0.0 ) then
      theta = halfpi
      return
    end if
    z = dot_product ( x, y ) / ( xn * yn )
    if ( n == 3 ) then
      if ( abs(z) <= s1 ) then
        theta = acos(z)
      else
        theta = asin ( norm2(cross(x,y)) / ( xn * yn ) )
        if ( z < 0 ) theta = pi - theta
      end if
    else if ( abs(z) <= s2 ) then
      theta = acos(z)
    else
      theta = 2.0d0 * atan2 ( norm2(x * yn - y * xn), norm2(x * yn + y * xn) )
      if ( z < 0 ) theta = pi - theta
    end if
  end function Vector_Angle_D

  real function Vector_Angle_S ( X, Y ) result ( Theta )
    real, intent(in) :: X(:), Y(:)
    integer :: N
    ! S1 ~ 0.7071 is the value at which arccos and arcsin methods have
    ! equal error.
    real, parameter :: S1 = 0.5e0 * sqrt(2.0e0)
    ! S2 ~ 0.6813 is the value at which arccos and arctan methods have
    ! equal error. The difference in the errors of the arccos method at
    ! S1 and S2 ~ 0.048 epsilon(1.0e0).
    real, parameter :: s2 = sqrt ( 2.0e0 * sqrt(3.0e0) - 3.0e0 )
    real :: XN, YN, Z
    n = size(x)
    xn = norm2(x)
    yn = norm2(y)
    if ( xn == 0.0 .or. yn == 0.0 ) then
      theta = halfpi
      return
    end if
    z = dot_product ( x, y ) / ( xn * yn )
    if ( n == 3 ) then
      if ( abs(z) <= s2 ) then
        theta = acos(z)
      else
        theta = asin ( norm2(cross(x,y)) / ( xn * yn ) )
        if ( z < 0 ) theta = pi - theta
      end if
    else if ( abs(z) <= s2 ) then
      theta = acos(z)
    else
      theta = 2.0e0 * atan2 ( norm2(x * yn - y * xn), norm2(x * yn + y * xn) )
      if ( z < 0 ) theta = pi - theta
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
! Revision 2.6  2019/09/18 22:49:01  vsnyder
! Add 'in radians' to the description
!
! Revision 2.5  2019/09/18 22:48:00  vsnyder
! Handle the case of zero norm
!
! Revision 2.4  2019/09/18 22:20:37  vsnyder
! Correct arctan method for angles > 90 degrees
!
! Revision 2.3  2019/09/17 22:20:13  vsnyder
! Correct the atan2 case
!
! Revision 2.2  2019/09/17 03:19:54  vsnyder
! More accurate switching points to minimize error
!
! Revision 2.1  2019/09/17 00:53:36  vsnyder
! Initial commit
!
