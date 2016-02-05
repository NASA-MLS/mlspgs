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
  public :: Line_And_Plane

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Line_And_Plane ( Plane, Line, Intersection, Intersect )

    ! Compute the intersections of Line with a Plane.  The plane is defined
    ! by three vertices in ECR coordinates (in any order).
    ! The Line is given by a vector to a point on it and a vector along it.
    ! The Intersection is given in ECR coordinates.

    ! See wvs-130 for a derivation.

    use Cross_m, only: Cross
    use Geolocation_0, only: ECR_t, RG

    type(ECR_t), intent(in) :: Plane(3) ! Three points in the plane
    type(ECR_t), intent(in) :: Line(2) ! The line is of the form
                                       ! Line(1) + t * Line(2), i.e., Line(1)
                                       ! is a vector to a point on the line,
                                       ! and Line(2) is a vector along the
                                       ! line.
    type(ECR_t), intent(out) :: Intersection ! ECR coordinates of intersection
                                       ! if Intersect > 0, else undefined.
    integer, intent(out) :: Intersect  ! -1 => no intersection,
                                       !  0 => line is in the plane,
                                       ! +1 => one intersection.

    real(rg) :: L_dot_N                ! L .dot. N
    real(rg) :: N(3)                   ! Normal to the plane
    real(rg) :: Z_dot_N                ! ( P_0 - L_0 ) .dot. N

    n = cross ( plane(2)%xyz - plane(1)%xyz, plane(3)%xyz - plane(1)%xyz )
    l_dot_n = dot_product ( line(2)%xyz, n )
    z_dot_n = dot_product ( plane(1)%xyz - line(2)%xyz, n )
    if ( l_dot_n == 0 ) then   ! Line and plane are parallel
      if ( z_dot_n == 0 ) then ! Line is in the plane
        intersect = 0          ! Caller chooses a point
      else                     ! Line and plane do not intersect
        intersect = -1
      end if
    else                       ! Line and plane intersect in one point
      intersection = ECR_t(line(1)%xyz + z_dot_n / l_dot_n * line(2)%xyz )
      intersect = 1
    end if

  end subroutine Line_And_Plane

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
! Revision 2.1  2016/02/05 03:38:42  vsnyder
! Initial commit
!
