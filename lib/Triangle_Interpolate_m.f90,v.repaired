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
module Triangle_Interpolate_m
!=============================================================================

  !{ \parskip 3pt
  !  Compute the coefficients to interpolate linearly within a triangle.  The
  !  interpolation coefficients are the normalized barycentric coordinates of
  !  the point with respect to the triangle.  Consider a point $(x, y)$ within
  !  the triangle.  Then the $i^{\text{th}}$ normalized barycentric coordinate
  !  $\lambda_i$ of $(x, y)$ is the ratio of the area of the triangle
  !  $( (x, y), (x_j, y_j), (x_k, y_k) )$ to the area of the entire triangle,
  !  where $i \neq j \neq k$.
  !
  !  The coordinates $\lambda_i$ can be computed by solving the system
  !%
  !  \begin{equation*}
  !    \left( \begin{array}{ccc}
  !           x_1 & x_2 & x_3 \\
  !           y_1 & y_2 & y_3 \\
  !           1   & 1   & 1   \\ \end{array} \right)
  !    \left( \begin{array}{c}
  !           \lambda_1 \\
  !           \lambda_2 \\
  !           \lambda_3   \\ \end{array} \right)
  !    =
  !    \left( \begin{array}{c}
  !           x \\
  !           y \\
  !           1 \\ \end{array} \right)
  !  \end{equation*}
  !%
  !  Explicitly, we have
  !%
  !  \begin{equation*}
  !    \lambda_k = \frac{(y_j-y_i)(x-x_i)   - (x_j-x_i)(y-y_i)}
  !                     {(y_j-y_i)(x_k-x_i) - (x_j-x_i)(y_k-y_i)}
  !  \end{equation*}
  !%
  !  where $i$, $j$, $k$ are distinct and in $\{1,2,3\}$.  Since
  !  $\sum_{k=1}^3 \lambda_k = 1$, we only calculate $\lambda_1$ and
  !  $\lambda_2$. It doesn't matter in in which order the vertices are
  !  numbered.  If any $\lambda_k<0$ or $\lambda_k>1$, the point $(x,y)$ is
  !  outside the triangle.
  !
  !  The interpolated value from $f(x_1,y_1)$, $f(x_2,y_2)$, $f(x_3,y_3)$ to
  !  $(x, y)$ is
  !%
  !  \begin{equation*}
  !  f(x,y) = \sum_{i=1}^3 \lambda_i f(x_i,y_i)
  !  \end{equation*}
  !%
  !  The system was introduced in 1827 by August Ferdinand M{\"o}bius.

  implicit NONE
  private

  public :: Triangle_Interpolate, Triangle_Interpolate_d, Triangle_Interpolate_s

  interface Triangle_Interpolate
    module procedure Triangle_Interpolate_d, Triangle_Interpolate_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  subroutine Triangle_Interpolate_d ( X, Y, XQ, YQ, W )
    ! We assume the triangle (X, Y) is not three points on one line.
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(3), Y(3) ! Coordinates of triangle vertices
    real(rk), intent(in) :: XQ, YQ     ! Coordinates to which to interpolate
    real(rk), intent(out) :: W(3)      ! Interpolation coefficients
    real(rk) :: D                      ! 1 / det( [ ( x1 - x3 )  ( x2 - x3 ) ]
                                       !          [ ( y1 - y3 )  ( y2 - y3 ) ]
    d = 1.0 / ( ( y(2) - y(3) ) * ( x(1) - x(3) ) - &
              & ( x(2) - x(3) ) * ( y(1) - y(3) ) )
    w(1) = d * ( ( y(2) - y(3) ) * ( xq - x(3) ) + ( x(3) - x(2) ) * ( yq - y(3) ) )
    w(2) = d * ( ( y(3) - y(1) ) * ( xq - x(3) ) + ( x(1) - x(3) ) * ( yq - y(3) ) )
    w(3) = 1.0 - w(1) - w(2)
    if ( .not. any(w >= 0) ) w = -w    ! faster than all(w < 0)
  end subroutine Triangle_Interpolate_d

  subroutine Triangle_Interpolate_s ( X, Y, XQ, YQ, W )
    ! We assume the triangle (X, Y) is not three points on one line.
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X(3), Y(3) ! Coordinates of triangle vertices
    real(rk), intent(in) :: XQ, YQ     ! Coordinates to which to interpolate
    real(rk), intent(out) :: W(3)      ! Interpolation coefficients
    real(rk) :: D                      ! 1 / det( [ ( x1 - x3 )  ( x2 - x3 ) ]
                                       !          [ ( y1 - y3 )  ( y2 - y3 ) ]
    d = 1.0 / ( ( y(2) - y(3) ) * ( x(1) - x(3) ) - &
              & ( x(2) - x(3) ) * ( y(1) - y(3) ) )
    w(1) = d * ( ( y(2) - y(3) ) * ( xq - x(3) ) + ( x(3) - x(2) ) * ( yq - y(3) ) )
    w(2) = d * ( ( y(3) - y(1) ) * ( xq - x(3) ) + ( x(1) - x(3) ) * ( yq - y(3) ) )
    w(3) = 1.0 - w(1) - w(2)
    if ( .not. any(w >= 0) ) w = -w    ! faster than all(w < 0)
  end subroutine Triangle_Interpolate_s

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

end module Triangle_Interpolate_m

! $Log$
! Revision 2.4  2016/11/02 00:38:23  vsnyder
! Clearer description, faster code
!
! Revision 2.3  2015/12/30 23:53:39  vsnyder
! Correct calculation of negative weights for outside points
!
! Revision 2.2  2015/11/12 20:45:38  vsnyder
! Remove unused declarations
!
! Revision 2.1  2015/11/12 20:14:18  vsnyder
! Initial commit
!
