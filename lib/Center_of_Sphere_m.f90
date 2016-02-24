! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Center_of_Sphere_m

  ! Determine the center of a sphere defined by three points and its
  ! radius.  There are two such spheres.  Assume the points are defined
  ! in ECR coordinates, and we wish the center of the sphere to be the
  ! one nearest the center of the Earth.

  ! For derivation see wvs-132.

  implicit NONE
  private

  public :: Center_of_Sphere, Circumcenter

  interface Center_of_Sphere
    module procedure Center_of_Sphere_3, Center_of_Sphere_3_V
    module procedure Center_of_Sphere_V, Center_of_Sphere_V_V
  end interface

  interface Circumcenter ! of three points in ECR
    module procedure Circumcenter_3, Circumcenter_3_N
    module procedure Circumcenter_V, Circumcenter_V_N
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Center_of_Sphere_3 ( P1, P2, P3, R, C )     
    ! Compute the center C of the sphere defined by {P1, P2, P3} and R that
    ! is nearest to the center of the Earth, i.e., the one for which |C| is
    ! minimum.

    ! There is an implicit assumption that the distance from {P1, P2, P3} to
    ! their circumcenter < R.  Otherwise, the circumcenter is returned.

    use Geolocation_0, only: ECR_t, RG

    type(ECR_t), intent(in) :: P1, P2, P3
    real(rg), intent(in) :: R
    type(ECR_t), intent(out) :: C

    type(ECR_t) :: N    ! A normal to the plane defined by {P1, P2, P3}
    real(rg) :: N2      ! |N|**2
    type(ECR_t) :: V    ! Vector from P3 to circumcenter of {P1, P2, P3}

    ! Compute the vector from P3 to the circumcenter, a normal to the plane
    ! defined by P1, P2, P3}, and the square of the length of the normal.
    call circumcenter ( p1, p2, p3, v, n, n2 )
    ! Now get the center of the sphere
    call center_of_sphere ( p3, v, n, n2, r, c )

  end subroutine Center_of_Sphere_3

  subroutine Center_of_Sphere_3_V ( P3, V, N, N2, R, C )
    ! Given the circumcenter V of three points {P1, P2, P3}, a normal N to
    ! the plane defined by those points, N2 = |N|^2, and a radius R, compute
    ! the center C of the sphere defined by {P1, P2, P3} and R that is
    ! nearest to the center of the Earth, i.e., the one for which |C| is
    ! minimum.

    ! There is an implicit assumption that the distance from {P1, P2, P3} to
    ! their circumcenter < R.  Otherwise, the circumcenter is returned.

    ! For derivation see wvs-132.

    use Geolocation_0, only: Dot_Product, ECR_t, RG

    type(ECR_t), intent(in) :: P3
    type(ECR_t), intent(in) :: V  ! Vector from P3 to circumcenter of {P1, P2, P3}
    type(ECR_t), intent(in) :: N  ! A normal to the plane defined by {P1, P2, P3}
    real(rg), intent(in) :: N2    ! |N|**2
    real(rg), intent(in) :: R     ! Desired radius
    type(ECR_t), intent(out) :: C ! Center of the sphere

    type(ECR_t) :: A, B ! Two solutions
    real(rg) :: T       ! Distance from circumcenter of {P1, P2, P3} to C
    real(rg) :: V2      ! |V|**2

    v2 = dot_product(v,v) ! ifort 15.0.2 doesn't like v .DOT. v

    t = ( r**2 - v2 ) / n2

    c = p3 + v   ! Start at the circumcenter, p(0) in wvs-132.
    if ( t > 0 ) then
      t = sqrt(t)
      ! We don't need A and B any more, so use them for the two solutions
      a = c + t * n
      b = c - t * n
      if ( dot_product(a,a) < dot_product(b,b) ) then
        c = a
      else
        c = b
      end if
    end if

  end subroutine Center_of_Sphere_3_V

  subroutine Center_of_Sphere_V_V ( V, X, N, N2, R, C )
    ! Given the circumcenter X of three points V, a normal N to the plane
    ! defined by V, N2 = |N|**2, and a radius R, compute the center C of
    ! the sphere defined by V and R that is nearest to the center of the
    ! Earth, i.e., the one for which |C| is minimum.
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: V(3)
    type(ECR_t), intent(in) :: X
    type(ECR_t), intent(in) :: N
    real(rg), intent(in) :: N2
    real(rg), intent(in) :: R
    type(ECR_t), intent(out) :: C
    call center_of_sphere ( v(3), x, n, n2, r, c )
  end subroutine Center_of_Sphere_V_V

  subroutine Center_of_Sphere_V ( V, R, C )
    ! Compute the center C of the sphere defined by V and R that is nearest
    ! to the center of the Earth, i.e., the one for which |C| is minimum.
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: V(3)
    real(rg), intent(in) :: R
    type(ECR_t), intent(out) :: C
    call center_of_sphere ( v(1), v(2), v(3), r, c )
  end subroutine Center_of_Sphere_V

  pure subroutine Circumcenter_3 ( P1, P2, P3, C )
    ! Compute the vector C from P3 to the circumcenter of the plane defined
    ! by {P1, P2, P3}.

    use Geolocation_0, only: Cross, Dot_Product, ECR_t, RG

    type(ECR_t), intent(in) :: P1, P2, P3
    type(ECR_t), intent(out) :: C  ! Vector from P3 to circumcenter

    type(ECR_t) :: A, B ! p1-p3, p2-p3
    real(rg) :: A2, B2  ! |A|**2 and |B|**2
    type(ECR_t) :: D    ! a2 * b - b2 * a
    type(ECR_t) :: N  ! Normal to the plane defined by {P1, P2, P3}.
    real(rg) :: N2    ! |N|**2

    a = p1 - p3
    b = p2 - p3

    a2 = dot_product(a,a) ! ifort 15.0.2 doesn't like a .DOT. a
    b2 = dot_product(b,b) ! ifort 15.0.2 doesn't like b .DOT. b

    n = cross(a,b)        ! ifort 15.0.2 doesn't like a .CROSS. b

    n2 = dot_product(n,n) ! ifort 15.0.2 doesn't like n .DOT. n

    d = a2 * b - b2 * a

    c = cross(d,n) / (2.0 / n2 ) ! ifort 15.0.2 doesn't like d .CROSS. n

  end subroutine Circumcenter_3

  pure subroutine Circumcenter_3_N ( P1, P2, P3, C, N, N2 )
    ! Compute the vector C from P3 to the circumcenter of the plane defined
    ! by {P1, P2, P3}, a vector N that is a normal to that plane, and |N|**2.

    use Geolocation_0, only: Cross, Dot_Product, ECR_t, RG

    type(ECR_t), intent(in) :: P1, P2, P3
    type(ECR_t), intent(out) :: C  ! Vector from P3 to circumcenter
    type(ECR_t), intent(out) :: N  ! Normal to the plane defined by {P1, P2, P3}.
    real(rg), intent(out) :: N2    ! |N|**2

    type(ECR_t) :: A, B ! p1-p3, p2-p3
    real(rg) :: A2, B2  ! |A|**2 and |B|**2
    type(ECR_t) :: D    ! a2 * b - b2 * a

    a = p1 - p3
    b = p2 - p3

    a2 = dot_product(a,a) ! ifort 15.0.2 doesn't like a .DOT. a
    b2 = dot_product(b,b) ! ifort 15.0.2 doesn't like b .DOT. b

    n = cross(a,b)        ! ifort 15.0.2 doesn't like a .CROSS. b

    n2 = dot_product(n,n) ! ifort 15.0.2 doesn't like n .DOT. n

    d = a2 * b - b2 * a

    c = cross(d,n) / (2.0 / n2 ) ! ifort 15.0.2 doesn't like d .CROSS. n

  end subroutine Circumcenter_3_N

  pure subroutine Circumcenter_V ( V, C )
    ! Compute the vector C from V(3) to the circumcenter of the plane defined
    ! by V, a vector N that is a normal to that plane, and |N|**2.

    use Geolocation_0, only: ECR_t

    type(ECR_t), intent(in) :: V(3)
    type(ECR_t), intent(out) :: C  ! Vector from V(3) to circumcenter
    call circumcenter ( v(1), v(2), v(3), c )
  end subroutine Circumcenter_V

  pure subroutine Circumcenter_V_N ( V, C, N, N2 )
    ! Compute the vector C from V(3) to the circumcenter of the plane defined
    ! by V, a vector N that is a normal to that plane, and |N|**2.

    use Geolocation_0, only: ECR_t, RG

    type(ECR_t), intent(in) :: V(3)
    type(ECR_t), intent(out) :: C  ! Vector from V(3) to circumcenter
    type(ECR_t), intent(out) :: N  ! Normal to the plane defined by V
    real(rg), intent(out) :: N2    ! |N|**2
    call circumcenter ( v(1), v(2), v(3), c, n, n2 )
  end subroutine Circumcenter_V_N

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Center_of_Sphere_m

! $Log$
! Revision 2.1  2016/02/24 01:19:50  vsnyder
! Initial commit
!
