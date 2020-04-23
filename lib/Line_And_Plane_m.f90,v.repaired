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
module Line_And_Plane_m
!=============================================================================

  private
  public :: Line_And_Plane, Line_Reflection

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Line_And_Plane ( Plane, Line, Intersect, Intersection, S )

    ! Compute the intersections of Line with a Plane.  The plane is defined
    ! by three vertices in ECR coordinates (in any order).
    ! The Line is given by a vector to a point on it (Line(1)) and a vector
    ! along it (Line(2)).  The Intersection is given in ECR coordinates.

    ! See wvs-130 for a derivation.

    use Geolocation_0, only: Cross
    use Geolocation_0, only: ECR_t, RG

    type(ECR_t), intent(in) :: Plane(3) ! Three points in the plane
    type(ECR_t), intent(in) :: Line(2) ! The line is of the form
                                       ! Line(1) + s * Line(2), i.e., Line(1)
                                       ! is a vector to a point on the line,
                                       ! and Line(2) is a vector along the
                                       ! line.
    integer, intent(out) :: Intersect  ! -1 => no intersection,
                                       !  0 => line is in the plane,
                                       ! +1 => one intersection.
    type(ECR_t), intent(out), optional :: Intersection ! ECR coordinates of
                                       ! intersection if Intersect > 0, else
                                       ! undefined.
    real(rg), intent(out), optional :: S ! S-value of intersection if
                                       ! Intersect > 0, else undefined.

    integer :: I
    real(rg) :: L_dot_N                ! L .dot. N
    real(rg) :: MyS
    type(ECR_t) :: N                   ! Normal to the plane
    real(rg) :: Z_dot_N                ! ( P_0 - L_0 ) .dot. N

!     ifort 16.0.2.181 doesn't like .CROSS. to have operands that are both
!     expressions in parentheses.
!     n = ( plane(2) - plane(1) ) .cross. ( plane(3) - plane(1) )
    n = cross( plane(2) - plane(1), plane(3) - plane(1) )
    l_dot_n = line(2) .dot. n
    i = 1
    if ( all(plane(1)%xyz == line(1)%xyz) ) i = 2
    z_dot_n = ( plane(i) - line(1) ) .dot. n
    if ( l_dot_n == 0 ) then   ! Line and plane are parallel
      if ( z_dot_n == 0 ) then ! Line is in the plane
        intersect = 0          ! Caller chooses a point
      else                     ! Line and plane do not intersect
        intersect = -1
      end if
    else                       ! Line and plane intersect in one point
      intersect = 1
      myS = z_dot_n / l_dot_n
      if ( present(intersection) ) intersection = line(1) + myS * line(2)
      if ( present(s) ) s = myS
    end if

  end subroutine Line_And_Plane

  subroutine Line_Reflection ( Line, Normal, Reflection )
    ! Given a vector parallel to a line, and a vector parallel to an unit
    ! normal to a plane, compute a vector parallel to the reflection of
    ! the line from the plane.  For derivation, see wvs-133.
    use Geolocation_0, only: ECR_t
    type(ECR_t), intent(in) :: Line    ! A vector parallel to the line
    type(ECR_t), intent(in) :: Normal  ! An unit normal to the plane
    type(ECR_t), intent(out) :: Reflection ! A vector parallel to the
                                       ! reflection, with the same length
                                       ! as Line
    Reflection = Line - 2 * ( Line .dot. Normal ) * Normal
  end subroutine Line_Reflection

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Line_And_Plane_m

! $Log$
! Revision 2.5  2020/04/23 22:57:46  vsnyder
! Use C instead of U to compute Z_dot_N. Choose a different p_0 if C == p_0.
!
! Revision 2.4  2016/03/25 00:07:47  vsnyder
! Remove unused USE name
!
! Revision 2.3  2016/03/03 21:39:38  vsnyder
! Change name of T argument to S (because it's arc length)
!
! Revision 2.2  2016/03/03 03:05:38  vsnyder
! Add Line_Reflection
!
! Revision 2.1  2016/02/05 03:38:42  vsnyder
! Initial commit
!
