! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Line_And_Ellipsoid_m

  !{ Compute the intersections of a line with an ellipsoid, or sphere. For an
  !  ellipsoid with center $\mathbf{p}_0$ and semi-minor axes $a$, $b$, and
  !  $c$, and a line given by $\mathbf{C} + t\, \mathbf{U}$, where
  !  $\mathbf{C}$ is a point on the line and $\mathbf{U}$ is a vector along
  !  the line, the intersections occur on the line at values of $t$ such that
  !%
  !  \begin{equation}
  !    (\mathbf{U}\, \mathbf{M})^T (\mathbf{U}\, \mathbf{M})\, t^2 +
  !    (\mathbf{V} \mathbf{M} )^T (\mathbf{U}\, \mathbf{M})\, t +
  !    (\mathbf{V} \mathbf{M} )^T
  !    (\mathbf{V} \mathbf{M} ) = 1
  !  \end{equation}
  !%
  !  where $\mathbf{M}$ is a diagonal matrix with the inverses of the
  !  semi-minor axes on the diagonal, and
  !  $\mathbf{V} = \mathbf{C} - \mathbf{p}_0$.
  !
  !  For the case of a sphere of radius $r$ this becomes
  !%
  !  \begin{equation}
  !    \mathbf{U} \cdot \mathbf{U} \, t^2 +
  !    \mathbf{V} \cdot \mathbf{U} \, t +
  !    \mathbf{V} \cdot \mathbf{V} = r^2 \,.
  !  \end{equation}
  !
  !  See wvs-131 for details of the derivation.

  implicit NONE
  private

  public Line_And_Ellipsoid, Line_And_Sphere

  interface Line_And_Ellipsoid
    module procedure Line_And_Earth_Geoid
    module procedure Line_And_Ellipsoid_RG
  end interface

  interface Line_And_Sphere
    module procedure Line_And_Sphere_RG
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Line_And_Earth_Geoid ( Line, Intersections )
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: Line(2) ! C + t U, two vectors
    type(ECR_t), intent(out), allocatable :: Intersections (:) ! 0..2 elements
    real(rg) :: A2, A1, A0
    type(ECR_t) :: UM, VM
    UM = ECR_t( Line(2)%xyz / [ A, A, B ] )
    VM = ECR_t( Line(1)%xyz / [ A, A, B ] )
    a2 = UM%dot(UM)       ! ifort 15.0.2.164 doesn't like UM .dot. UM
    a1 = VM%dot(UM)       ! ifort 15.0.2.164 doesn't like VM .dot. UM
    a0 = VM%dot(VM) - 1.0 ! ifort 15.0.2.164 doesn't like VM .dot. VM
    call solve ( a2, a1, a0, line, intersections )
  end subroutine Line_And_Earth_Geoid

  subroutine Line_And_Ellipsoid_RG ( Axes, Line, Intersections )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(in) :: Line(2) ! V + t U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V
    type(ECR_t), intent(out), allocatable :: Intersections (:) ! 0..2 elements
    real(rg) :: A2, A1, A0
    type(ECR_t) :: UM, VM
    UM = ECR_t( Line(2)%xyz / axes )
    VM = ECR_t( Line(1)%xyz / axes )
    a2 = UM%dot(UM)       ! ifort 15.0.2.164 doesn't like UM .dot. UM
    a1 = VM%dot(UM)       ! ifort 15.0.2.164 doesn't like VM .dot. UM
    a0 = VM%dot(VM) - 1.0 ! ifort 15.0.2.164 doesn't like VM .dot. VM
    call solve ( a2, a1, a0, line, intersections )
  end subroutine Line_And_Ellipsoid_RG

  subroutine Line_And_Sphere_RG ( Radius, Line, Intersections )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Radius     ! In same units as Center
    type(ECR_t), intent(in) :: Line(2) ! V + t U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V
    type(ECR_t), intent(out), allocatable :: Intersections (:) ! 0..2 elements
    real(rg) :: A2, A1, A0
    a2 = line(2)%dot(line(2)) ! ifort 15.0.2.164 doesn't like line(2) .dot. line(2)
    a1 = line(1)%dot(line(2))       ! ifort 15.0.2.164 doesn't like line(1) .dot. line(2)
    a0 = line(1)%dot(line(1)) - radius**2 ! ifort 15.0.2.164 doesn't like line(1) .dot. line(1)
    call solve ( a2, a1, a0, line, intersections )
  end subroutine Line_And_Sphere_RG

  ! =====     Private Procedure     ====================================

  ! Hopefully inlined above
  subroutine Solve ( A2, A1, A0, Line, Intersections )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: A2, A1, A0
    type(ECR_t), intent(in) :: Line(2)
    type(ECR_t), intent(out), allocatable :: Intersections (:) ! 0..2 elements
    real(rg) :: D
        d = a1**2 - 4.0 * a2 * a0
    if ( d < 0.0 ) then
      allocate ( intersections(0) )
    else if ( d == 0.0 ) then
      allocate ( intersections(1) )
      intersections(1) = line(1) + (-0.5 * a1 / a2) * line(2)
    else
      allocate ( intersections(2) )
      d = sqrt(d)
      intersections(1) = line(1) + (0.5 *(- a1 + d) / a2) * line(2)
      intersections(2) = line(1) + (0.5 *(- a1 - d) / a2) * line(2)
    end if
  end subroutine Solve

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Line_And_Ellipsoid_m

! $Log$
! Revision 2.1  2016/02/24 01:19:50  vsnyder
! Initial commit
!
