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
  !  the point with respect to the triangle.
  !
  !  The unnormalized barycentric coordinates $\hat\lambda_i$
  !  of a point $(x, y)$ within a triangle having vertices $(x_1, y_1)$,
  !  $(x_2, y_2)$ and $(x_3, y_3)$ are
  !%
  !  \begin{equation*}
  !  \hat\lambda_k = (y_j-y_i)(x-x_i) - (x_j-x_i)(y-y_i)
  !  \end{equation*}
  !%
  !  where $i$, $j$, $k$ are distinct and in $\{1,2,3\}$.  The signs depend
  !  upon whether the vertices are numbered in clockwise or anticlockwise
  !  order.  Normalized coordinates $\lambda_i$ are obtained by dividing the
  !  unnormalized coordinates by their sum.  Therefore, it doesn't matter in
  !  which order the vertices are numbered.  If any of the normalized
  !  barycentric coordinates are negative, the point $(x,y)$ is outside the
  !  triangle.
  !
  !  The interpolated value from $f(x_1,y_1)$, $f(x_2,y_2)$, $f(x_3,y_3)$ to
  !  $(x, y)$ is
  !%
  !  \begin{equation*}
  !  f(x,y) = \sum_{i=1}^3 \lambda_i f(x_i,y_i)
  !  \end{equation*}
  !%
  !  where $\lambda_i$ are normalized barycentric coordinates.
  !
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
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(3), Y(3) ! Coordinates of triangle vertices
    real(rk), intent(in) :: XQ, YQ     ! Coordinates to which to interpolate
    real(rk), intent(out) :: W(3)      ! Interpolation coefficients
    integer :: I, J                    ! Subscripts
    do i = 1, 3
      j = mod(i,3) + 1
      w(6-i-j) = (x(j) - x(i)) * ( yq - y(i) ) - (y(j) - y(i)) * ( xq - x(i) )
    end do
    if ( .not. any(w >= 0) ) w = -w    ! faster than all(w < 0)
    w = w / sum(abs(w))
  end subroutine Triangle_Interpolate_d

  subroutine Triangle_Interpolate_s ( X, Y, XQ, YQ, W )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X(3), Y(3) ! Coordinates of triangle vertices
    real(rk), intent(in) :: XQ, YQ     ! Coordinates to which to interpolate
    real(rk), intent(out) :: W(3)      ! Interpolation coefficients
    integer :: I, J                    ! Subscripts
    do i = 1, 3
      j = mod(i,3) + 1
      w(6-i-j) = (x(j) - x(i)) * ( yq - y(i) ) - (y(j) - y(i)) * ( xq - x(i) )
    end do
    if ( .not. any(w >= 0) ) w = -w    ! faster than all(w < 0)
    w = w / sum(abs(w))
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
! Revision 2.3  2015/12/30 23:53:39  vsnyder
! Correct calculation of negative weights for outside points
!
! Revision 2.2  2015/11/12 20:45:38  vsnyder
! Remove unused declarations
!
! Revision 2.1  2015/11/12 20:14:18  vsnyder
! Initial commit
!
