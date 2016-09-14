! Copyright 2015, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Line_and_Polygon
!=============================================================================

  implicit NONE
  private

  ! Report which edge, if any, of a polygon intersects a line segment, else
  ! report zero.

  public :: Line_Intersects_Any_Edge, Lines_Intersect, Perpendicular_Intersection

  interface Line_Intersects_Any_Edge
    module procedure Line_Intersects_Any_Edge_d, Line_Intersects_Any_Edge_s
  end interface

  interface Lines_Intersect
    module procedure Lines_Intersect_d, Lines_Intersect_s
  end interface

  interface Perpendicular_Intersection
    module procedure Perpendicular_Intersection_d, Perpendicular_Intersection_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  integer function Line_Intersects_Any_Edge_d ( X, Y, PX, PY, Ignore ) result ( I )
    ! Line intersects which edge of Polygon?  Returns zero if none.
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(2), Y(2)   ! Ends of the line
    real(rk), intent(in) :: PX(:), PY(:) ! Vertices of the polygon
    logical, intent(in), optional :: Ignore(:) ! Ignore edge (i,1+mod(i,n))

    real(rk) :: DXL, DYL     ! Deltas for line
    integer :: N             ! min(size(px),size(py))

    dxl = x(1) - x(2)
    dyl = y(1) - y(2)
    n = min(size(px),size(py))
    if ( present(ignore) ) then
      do i = 1, n - 1
        if ( .not. ignore(i) ) then
          if ( test ( px(i), px(i+1), py(i), py(i+1) ) ) return
        end if
      end do
      if ( .not. ignore(i) ) then
        if ( test ( px(n), px(1), py(n), py(1) ) ) return
      end if
    else
      do i = 1, n - 1
        if ( test ( px(i), px(i+1), py(i), py(i+1) ) ) return
      end do
      if ( test ( px(n), px(1), py(n), py(1) ) ) return
    end if

    i = 0

  contains

    logical function Test ( PX1, PX2, PY1, PY2 )
      real(rk), intent(in) :: PX1, PX2, PY1, PY2 ! Polygon edge end points
      real(rk) :: A, B       ! det(line), det(poly)
      real(rk) :: D          ! det(dLine, dPoly)
      real(rk) :: DXP, DYP   ! Deltas for polygon
      real(rk) :: TX, TY     ! D * ( Coordinates of intersection )
      dxp = px1 - px2
      dyp = py1 - py2
      d = dxl * dyp - dyl * dxp
      if ( abs(d) >= max( maxval(abs(x)), maxval(abs(y)) ) * epsilon(1.0_rk) ) then

        ! lines are not parallel
        a = x(1) * y(2) - x(2) * y(1)
        b = px1 * py2   - px2 * py1
        tx = a * dxp - b * dxl
        ty = a * dyp - b * dyl
        test = ( tx-d*x(1) >= 0.0_rk .eqv. d*x(2)-tx >= 0.0_rk ) .and. &
               ( ty-d*y(1) >= 0.0_rk .eqv. d*y(2)-ty >= 0.0_rk ) .and. &
               ( tx-d*px1  >= 0.0_rk .eqv. d*px2-tx  >= 0.0_rk ) .and. &
               ( ty-d*py1  >= 0.0_rk .eqv. d*py2-ty  >= 0.0_rk )


      else

        ! Lines are too close to parallel
        test = ( x(1) - px1 ) * ( x(1) - px2 ) <= 0.0_rk .and. &
             & ( y(1) - py1 ) * ( y(1) - py2 ) <= 0.0_rk .or. &
             & ( x(2) - px1 ) * ( x(2) - px2 ) <= 0.0_rk .and. &
             & ( y(2) - py1 ) * ( y(2) - py2 ) <= 0.0_rk .or. &
             & ( px1 - x(1) ) * ( px1 - x(2) ) <= 0.0_rk .and. &
             & ( py1 - y(1) ) * ( py1 - y(2) ) <= 0.0_rk .or. &
             & ( px2 - x(1) ) * ( px2 - x(2) ) <= 0.0_rk .and. &
             & ( py2 - y(1) ) * ( py2 - y(2) ) <= 0.0_rk

        if ( .not. test ) return ! no end point of one line is within
                                 ! end points of the other line

        ! Compute the distance from ( x(1), y(1) ) to the line
        ! defined by ( px1, py1 ) and ( px2, py2 ).
        d = abs( dxp * ( py1 - y(1) ) - dyp * ( px1 - x(1) ) ) / &
            norm2 ( [ dxp, dyp ] )
        test = d <= epsilon(1.0_rk) ! Lines are essentially collinear

      end if
    end function Test

  end function Line_Intersects_Any_Edge_d

  integer function Line_Intersects_Any_Edge_s ( X, Y, PX, PY, Ignore  ) result ( I )
    ! Does Line intersect any edge of Polygon?
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X(2), Y(2)   ! Ends of the line
    real(rk), intent(in) :: PX(:), PY(:) ! Vertices of the polygon
    logical, intent(in), optional :: Ignore(:) ! Ignore edge (i,1+mod(i,n))

    real(rk) :: DXL, DYL     ! Deltas for line
    integer :: N             ! min(size(px),size(py))

    dxl = x(1) - x(2)
    dyl = y(1) - y(2)
    n = min(size(px),size(py))
    if ( present(ignore) ) then
      do i = 1, n - 1
        if ( .not. ignore(i) ) then
          if ( test ( px(i), px(i+1), py(i), py(i+1) ) ) return
        end if
      end do
      if ( .not. ignore(i) ) then
        if ( test ( px(n), px(1), py(n), py(1) ) ) return
      end if
    else
      do i = 1, n - 1
        if ( test ( px(i), px(i+1), py(i), py(i+1) ) ) return
      end do
      if ( test ( px(n), px(1), py(n), py(1) ) ) return
    end if

    i = 0

  contains

    logical function Test ( PX1, PX2, PY1, PY2 )
      real(rk), intent(in) :: PX1, PX2, PY1, PY2 ! Polygon edge end points
      real(rk) :: A, B       ! det(line), det(poly)
      real(rk) :: D          ! det(dLine, dPoly)
      real(rk) :: DXP, DYP   ! Deltas for polygon
      real(rk) :: TX, TY     ! D * ( Coordinates of intersection )
      dxp = px1 - px2
      dyp = py1 - py2
      if ( abs(d) >= max( maxval(abs(x)), maxval(abs(y)) ) * epsilon(1.0_rk) ) then

        ! lines are not parallel
        a = x(1) * y(2) - x(2) * y(1)
        b = px1 * py2   - px2 * py1
        tx = a * dxp - b * dxl
        ty = a * dyp - b * dyl
        test = ( tx-d*x(1) >= 0.0_rk .eqv. d*x(2)-tx >= 0.0_rk ) .and. &
               ( ty-d*y(1) >= 0.0_rk .eqv. d*y(2)-ty >= 0.0_rk ) .and. &
               ( tx-d*px1  >= 0.0_rk .eqv. d*px2-tx  >= 0.0_rk ) .and. &
               ( ty-d*py1  >= 0.0_rk .eqv. d*py2-ty  >= 0.0_rk )

      else

        ! Lines are too close to parallel
        test = ( x(1) - px1 ) * ( x(1) - px2 ) <= 0.0_rk .and. &
             & ( y(1) - py1 ) * ( y(1) - py2 ) <= 0.0_rk .or. &
             & ( x(2) - px1 ) * ( x(2) - px2 ) <= 0.0_rk .and. &
             & ( y(2) - py1 ) * ( y(2) - py2 ) <= 0.0_rk .or. &
             & ( px1 - x(1) ) * ( px1 - x(2) ) <= 0.0_rk .and. &
             & ( py1 - y(1) ) * ( py1 - y(2) ) <= 0.0_rk .or. &
             & ( px2 - x(1) ) * ( px2 - x(2) ) <= 0.0_rk .and. &
             & ( py2 - y(1) ) * ( py2 - y(2) ) <= 0.0_rk

        if ( .not. test ) return ! no end point of one line is within
                                 ! end points of the other line

        ! Compute the distance from ( x(1), y(1) ) to the line
        ! defined by ( px1, py1 ) and ( px2, py2 ).
        d = abs( dxp * ( py1 - y(1) ) - dyp * ( px1 - x(1) ) ) / &
            norm2 ( [ dxp, dyp ] )
        test = d <= epsilon(1.0_rk) ! Lines are essentially collinear

      end if
    end function Test

  end function Line_Intersects_Any_Edge_s

  logical function Lines_Intersect_d ( X, Y, P, B1, B2 ) result ( Yes )
    ! Compute (P1,P2) the point where (X1,Y1) (X2,Y2) intersects
    ! (X3,Y3) (X4,Y4).  Result is false if lines are parallel and not
    ! colinear.
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(:), Y(:) ! Two points on each line
    real(rk), intent(out) :: P(:)      ! (X,Y) of intersection
    logical, intent(out), optional :: B1 ! Intersect between (x1,y1) and (x2,y2)
    logical, intent(out), optional :: B2 ! Intersect between (x3,y3) and (x4,y4)
    real(rk) :: DX12, DY12, DX34, DY34 ! Deltas for lines
    real(rk) :: D                      ! det | dx12 dx34 |
                                       !     | dy12 dy34 |
    real(rk) :: D12                    ! det | x1 x2 |
                                       !     | y1 y2 |
    real(rk) :: D34                    ! det | x3 x4 |
                                       !     | y3 y4 |
    logical :: T1, T2                  ! Local values for B1, B2
    dx12 = x(1) - x(2)
    dx34 = x(3) - x(4)
    dy12 = y(1) - y(2)
    dy34 = y(3) - y(4)
    d = dx12 * dy34 - dx34 * dy12
    d12 = x(1) * y(2) - x(2) * y(1)
    d34 = x(3) * y(4) - x(4) * y(3)
    p(1) = ( d12 * dx34 - d34 * dx12 )
    p(2) = ( d12 * dy34 - d34 * dy12 )
    yes = abs(d) >= max( maxval(abs(x)), maxval(abs(y)) ) * epsilon(1.0_rk)
    if ( yes ) then
      p = p / d
      if ( present(b1) ) b1 = (p(1)-x(1)) * (p(1)-x(2)) <= 0 .and. &
                              (p(2)-y(1)) * (p(2)-y(2)) <= 0
      if ( present(b2) ) b2 = (p(1)-x(3)) * (p(1)-x(4)) <= 0 .and. &
                              (p(2)-y(3)) * (p(2)-y(4)) <= 0
    else
      ! Lines are too close to parallel.  This will not work if d == 0.
      t1 = (p(1)-d*x(1)) * (p(1)-d*x(2)) <= 0 .and. &
           (p(2)-d*y(1)) * (p(2)-d*y(2)) <= 0
      t2 = (p(1)-d*x(3)) * (p(1)-d*x(4)) <= 0 .and. &
           (p(2)-d*y(3)) * (p(2)-d*y(4)) <= 0
      yes = t1 .or. t2
      if ( .not. yes ) return
      d = abs( dx12 * ( y(2) - y(3) ) - dy12 * ( x(2) - x(3) ) ) / &
          norm2( [ dx12, dy12 ] )
      yes = d >= epsilon(1.0_rk) ! Lines are essentially collinear
      p = [ x(1), y(1) ]
      if ( present(b1) ) b1 = t1
      if ( present(b2) ) b2 = t2
    end if
  end function Lines_Intersect_d

  logical function Lines_Intersect_s ( X, Y, P, B1, B2 ) result ( Yes )
    ! Compute (P1,P2) the point where (X1,Y1) (X2,Y2) intersects
    ! (X3,Y3) (X4,Y4).  Result is false if lines are parallel and not
    ! colinear.
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X(:), Y(:) ! Two points on each line
    real(rk), intent(out) :: P(:)      ! (X,Y) of intersection
    logical, intent(out), optional :: B1 ! Intersect between (x1,y1) and (x2,y2)
    logical, intent(out), optional :: B2 ! Intersect between (x3,y3) and (x4,y4)
    real(rk) :: DX12, DY12, DX34, DY34 ! Deltas for lines
    real(rk) :: D                      ! det | dx12 dx34 |
                                       !     | dy12 dy34 |
    real(rk) :: D12                    ! det | x1 x2 |
                                       !     | y1 y2 |
    real(rk) :: D34                    ! det | x3 x4 |
                                       !     | y3 y4 |
    logical :: T1, T2                  ! Local values for B1, B2
    dx12 = x(1) - x(2)
    dx34 = x(3) - x(4)
    dy12 = y(1) - y(2)
    dy34 = y(3) - y(4)
    d = dx12 * dy34 - dx34 * dy12
    d12 = x(1) * y(2) - x(2) * y(1)
    d34 = x(3) * y(4) - x(4) * y(3)
    p(1) = ( d12 * dx34 - d34 * dx12 )
    p(2) = ( d12 * dy34 - d34 * dy12 )
    yes = abs(d) >= max( maxval(abs(x)), maxval(abs(y)) ) * epsilon(1.0_rk)
    if ( yes ) then
      p = p / d
      if ( present(b1) ) b1 = (p(1)-x(1)) * (p(1)-x(2)) <= 0 .and. &
                              (p(2)-y(1)) * (p(2)-y(2)) <= 0
      if ( present(b2) ) b2 = (p(1)-x(3)) * (p(1)-x(4)) <= 0 .and. &
                              (p(2)-y(3)) * (p(2)-y(4)) <= 0
    else
      ! Lines are too close to parallel.  This will not work if d == 0.
      t1 = (p(1)-d*x(1)) * (p(1)-d*x(2)) <= 0 .and. &
           (p(2)-d*y(1)) * (p(2)-d*y(2)) <= 0
      t2 = (p(1)-d*x(3)) * (p(1)-d*x(4)) <= 0 .and. &
           (p(2)-d*y(3)) * (p(2)-d*y(4)) <= 0
      yes = t1 .or. t2
      if ( .not. yes ) return
      d = abs( dx12 * ( y(2) - y(3) ) - dy12 * ( x(2) - x(3) ) ) / &
          norm2( [ dx12, dy12 ] )
      yes = d >= epsilon(1.0_rk) ! Lines are essentially collinear
      p = [ x(1), y(1) ]
      if ( present(b1) ) b1 = t1
      if ( present(b2) ) b2 = t2
    end if
  end function Lines_Intersect_s

  subroutine Perpendicular_Intersection_d ( X, Y, P, Between )
    ! Calculate (P1,P2) where a line perpendicular to (X1,Y1), (X2,Y2) and
    ! passing through (X3,Y3) intersects the line (X1,Y1), (X2,Y2).
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: X(3), Y(3)
    real(rk), intent(out) :: P(2) ! The intersection
    logical, intent(out), optional :: Between ! P is between (X1,Y1) and (X2,Y2)
    real(rk) :: Det    ! x(1)*y(2) - x(2)*y(1) = det | x(1) x(2) |
                       !                     | y(1) y(2) |
    real(rk) :: DX, DY ! x(1)-x(2), y(1)-y(2)
    real(rk) :: L2     ! dx**2 + dy**2, (length of (X1,Y1), (X2,Y2) ) **2
    real(rk) :: T      ! dx*x(3) + dy*y(3)
    dx = x(1) - x(2)
    dy = y(1) - y(2)
    l2 = dx**2 + dy**2
    if ( l2 == 0.0 ) then
      p = [ x(1), y(1) ]
      if ( present(between) ) between = .true.
      return
    end if
    det = x(1)*y(2) - x(2)*y(1)
    t = dx*x(3) + dy*y(3)
    p(1) = ( dx * t - dy * det ) / l2
    p(2) = ( dy * t + dx * det ) / l2
    if ( present(between) ) &
      & between = ( p(1) - x(1) ) * ( p(1) - x(2) ) <= 0.0 .and. &
                & ( p(2) - y(1) ) * ( p(2) - y(2) ) <= 0.0
  end subroutine Perpendicular_Intersection_d

  subroutine Perpendicular_Intersection_s ( X, Y, P, Between )
    ! Calculate (P1,P2) where a line perpendicular to (X1,Y1), (X2,Y2) and
    ! passing through (X3,Y3) intersects the line (X1,Y1), (X2,Y2).
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: X(3), Y(3)
    real(rk), intent(out) :: P(2) ! The intersection
    logical, intent(out), optional :: Between ! P is between (X1,Y1) and (X2,Y2)
    real(rk) :: Det    ! x(1)*y(2) - x(2)*y(1) = det | x(1) x(2) |
                       !                     | y(1) y(2) |
    real(rk) :: DX, DY ! x(1)-x(2), y(1)-y(2)
    real(rk) :: L2     ! dx**2 + dy**2, (length of (X1,Y1), (X2,Y2) ) **2
    real(rk) :: T      ! dx*x(3) + dy*y(3)
    dx = x(1) - x(2)
    dy = y(1) - y(2)
    l2 = dx**2 + dy**2
    if ( l2 == 0.0 ) then
      p = [ x(1), y(1) ]
      if ( present(between) ) between = .true.
      return
    end if
    det = x(1)*y(2) - x(2)*y(1)
    t = dx*x(3) + dy*y(3)
    p(1) = ( dx * t - dy * det ) / l2
    p(2) = ( dy * t + dx * det ) / l2
    if ( present(between) ) &
      & between = ( p(1) - x(1) ) * ( p(1) - x(2) ) <= 0.0 .and. &
                & ( p(2) - y(1) ) * ( p(2) - y(2) ) <= 0.0
  end subroutine Perpendicular_Intersection_s

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

end module Line_and_Polygon

! $Log$
! Revision 2.2  2016/09/14 19:49:43  vsnyder
! Correct detection of collinear line segments
!
! Revision 2.1  2015/12/30 23:52:03  vsnyder
! Initial commit
!
