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
module PnPoly_m
!=============================================================================

  implicit NONE
  private

  ! Determine whether a point is inside a polygon

  ! Method:
  ! A vertical line is drawn thru the point in question. If it
  ! crosses the polygon an odd number of times, then the
  ! point is inside the polygon.

  public :: PnPoly, PnPoly_d, PnPoly_s
  public :: PnPoly2, PnPoly2_d, PnPoly2_s

  ! PnPoly is not as fast as PnPoly2, but it reliably distinguishes
  ! points on the edge or at a vertex.  PnPoly2 doesn't have a separate
  ! result value for the edge or vertex case.  It sometimes reports
  ! inside and sometimes reports them outside.

  interface PnPoly
    module procedure PnPoly_d, PnPoly_s
  end interface

  interface PnPoly2
    module procedure PnPoly2_d, PnPoly2_s
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

  pure integer function PnPoly_d ( PX, PY, XX, YY ) result ( In_or_Out )
    ! Result: -1 => (px,py) outside
    !         +1 => (px,py) inside,
    !          0 => (px,py) on an edge or at a vertex

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: PX, PY        ! Coordinates of point in question
    real(rk), intent(in) :: XX(:), YY(:)  ! Vertices of the polygon, in order

    include 'PnPoly.f9h'

  end function PnPoly_d

  pure integer function PnPoly_s ( PX, PY, XX, YY ) result ( In_or_Out )
    ! Result: -1 => (px,py) outside
    !         +1 => (px,py) inside,
    !          0 => (px,py) on an edge or at a vertex

    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: PX, PY        ! Coordinates of point in question
    real(rk), intent(in) :: XX(:), YY(:)  ! Vertices of the polygon, in order

    include 'PnPoly.f9h'

  end function PnPoly_s

  pure logical function PnPoly2_d ( PX, PY, X, Y ) result ( In )
    ! Result: "(px,py) is inside"
    ! If (PX,PY) is on the border or at a vertex, it doesn't reliably
    ! distingish whether it is inside or outside.

    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: PX, PY      ! Coordinates of point in question
    real(rk), intent(in) :: X(:), Y(:)  ! Vertices of the polygon, in order

    integer :: I, J, N

    n = min(size(x), size(y))

    in = .false.
    do i = 1, n
      j = mod(i,n) + 1
      if ( y(i) > py .neqv. y(j) > py ) then
        in = in .neqv. &
          & px < ( x(j) - x(i) ) * ( py - y(i) ) / ( y(j) - y(i) ) + x(i)
      end if
    end do

  end function PnPoly2_d

  pure logical function PnPoly2_s ( PX, PY, X, Y ) result ( In )
    ! Result: "(px,py) is inside"
    ! If (PX,PY) is on the border or at a vertex, it doesn't reliably
    ! distingish whether it is inside or outside.

    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: PX, PY      ! Coordinates of point in question
    real(rk), intent(in) :: X(:), Y(:)  ! Vertices of the polygon, in order

    integer :: I, J, N

    n = min(size(x), size(y))

    in = .false.
    do i = 1, n
      j = mod(i,n) + 1
      if ( y(i) > py .neqv. y(j) > py ) then
        in = in .neqv. &
          & px < ( x(j) - x(i) ) * ( py - y(i) ) / ( y(j) - y(i) ) + x(i)
      end if
    end do

  end function PnPoly2_s

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

end module PnPoly_m

!  !
