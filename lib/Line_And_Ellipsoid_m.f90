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
  !
  !  Compute the point on a line that is nearest to an ellipsoid.
  !  This occurs at
  !%
  !  \begin{equation}
  !   t =
  ! -\frac{\left(\mathbf{M}\mathbf{D}\right)^T
  !        \left(\mathbf{M}\mathbf{U}\right)}
  !       {\left(\mathbf{M}\mathbf{U}\right)^T
  !        \left(\mathbf{M}\mathbf{U}\right)} =
  ! -\frac{\frac{x_1 x_2}{a^2} + \frac{y_1 y_2}{b^2} + \frac{z_1 z_2}{c^2}}
  !       {\frac{x_2^2}{a^2} + \frac{y_2^2}{b^2} + \frac{z_2^2}{c^2}} \,.
  ! \end{equation}
  !
  !  where $\mathbf{D} = \mathbf{C} - \mathbf{p}(0)$ and $\mathbf{p}(0)$ is
  !  the center of the ellipsoid.
  !  See wvs-030 for details of the derivation.

  implicit NONE
  private

  public :: Ellipsoid_Gradient
  public :: Line_And_Ellipsoid, Line_And_Sphere, Line_Nearest_Ellipsoid

  interface Ellipsoid_Gradient
    module procedure Earth_Geoid_Gradient
    module procedure Ellipsoid_Gradient_RG
  end interface

  interface Line_And_Ellipsoid
    module procedure Line_And_Earth_Geoid
    module procedure Line_And_Ellipsoid_RG
  end interface

  interface Line_And_Sphere
    module procedure Line_And_Sphere_RG
  end interface

  interface Line_Nearest_Ellipsoid
    module procedure Line_Nearest_Earth_Geoid
    module procedure Line_Nearest_Ellipsoid_RG
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Earth_Geoid_Gradient ( Where, Grad )
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t
    type(ECR_t), intent(in) :: Where ! Where on the Geoid the gradient is desired
    type(ECR_t), intent(out) :: Grad ! Vector parallel to the gradient
    grad = ECR_t ( where%xyz / [ A, A, B ]**2 )
  end subroutine Earth_Geoid_Gradient

  subroutine Ellipsoid_Gradient_RG ( Axes, Center, Where, Grad )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(in) :: Center  ! Center of the ellipsoid
    type(ECR_t), intent(in) :: Where ! Where on the Geoid the gradient is desired
    type(ECR_t), intent(out) :: Grad ! Vector parallel to the gradient
    grad = ECR_t ( ( where%xyz - center%xyz ) / axes**2 )
  end subroutine Ellipsoid_Gradient_RG

  subroutine Line_And_Earth_Geoid ( Line, Intersections, T )
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: Line(2) ! C + t U, two vectors
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: T(:) ! T-values of intersections
    real(rg) :: A2, A1, A0
    type(ECR_t) :: UM, VM
    UM = ECR_t( Line(2)%xyz / [ A, A, B ] )
    VM = ECR_t( Line(1)%xyz / [ A, A, B ] )
    a2 = UM .dot. UM
    a1 = VM .dot. UM
    a0 = ( VM .dot. VM ) - 1.0
    call solve ( a2, a1, a0, line, intersections, t )
  end subroutine Line_And_Earth_Geoid

  subroutine Line_And_Ellipsoid_RG ( Axes, Center, Line, Intersections, T )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(in) :: Center  ! Center of the ellipsoid
    type(ECR_t), intent(in) :: Line(2) ! V + t U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: T(:) ! T-values of intersections
    real(rg) :: A2, A1, A0
    type(ECR_t) :: D ! Line(1) - Center
    type(ECR_t) :: UM, VM
    d = line(1) - center
    UM = ECR_t( Line(2)%xyz / axes )
    VM = ECR_t( d%xyz / axes )
    a2 = UM .dot. UM
    a1 = VM .dot. UM
    a0 = ( VM .dot. VM ) - 1.0
    call solve ( a2, a1, a0, line, intersections, t )
  end subroutine Line_And_Ellipsoid_RG

  subroutine Line_And_Sphere_RG ( Radius, Center, Line, Intersections, T )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Radius     ! In same units as Center
    type(ECR_t), intent(in) :: Center  ! Center of the sphere
    type(ECR_t), intent(in) :: Line(2) ! V + t U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: T(:) ! T-values of intersections
    real(rg) :: A2, A1, A0
    type(ECR_t) :: D ! Line(1) - Center
    d = line(1) - center
    a2 = line(2) .dot. line(2)
    a1 = d .dot. line(2)
    a0 = ( d .dot. d ) - radius**2
    call solve ( a2, a1, a0, line, intersections, t )
  end subroutine Line_And_Sphere_RG

  subroutine Line_Nearest_Earth_Geoid ( Line, T, R, Grad, Intersect )
    ! See wvs-030 for derivation.
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(in) :: Line(2) ! C + t U, two vectors
    real(rg), intent(out) :: T ! T-value of nearest point
    real(rg), intent(out) :: R ! R < 1 => C + t U is an intersection
    type(ECR_t), intent(out), optional :: Grad ! If R >= 1, gradient to the
                               ! Geoid that intersects Line(1) + t * Line(2),
                               ! else undefined
    type(ECR_t), intent(out), optional :: Intersect ! If R >= 1, where G
                               ! intersects the Geoid, i.e., the point on
                               ! the Geoid nearest to Line, else undefined
    type(ECR_t) :: N(2)        ! Normal to the ellipsoid
    type(ECR_t) :: PM, UM, VM
    real(rg), allocatable :: T_int(:) ! Intersection of N with the Geoid
    UM = ECR_t( line(2)%xyz / [ A, A, B ] )
    VM = ECR_t( line(1)%xyz / [ A, A, B ] )
    t = - ( UM .dot. VM ) / ( UM .dot. UM )
    PM = ECR_t( ( line(1)%xyz + t * line(2)%xyz ) / [ A, A, B ] )
    r = sqrt(PM .dot. PM)
    if ( r >= 1 .and. ( present(intersect) .or. present(grad) ) ) then
      n(2) = ECR_t ( PM%xyz / [ A, A, B ] ) ! Gradient to Geoid
      if ( present(grad) ) grad = n(2)
      if ( present(intersect) ) then
        n(1) = line(1) + t * line(2)
        call line_and_ellipsoid ( n, t=t_int )
        intersect = n(1) + t_int(minloc(abs(t_int),1)) * n(2)
      end if
    end if
  end subroutine Line_Nearest_Earth_Geoid

  subroutine Line_Nearest_Ellipsoid_RG ( Axes, Center, Line, T, R, Grad, Intersect )
    ! See wvs-030 for derivation.
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(in) :: Center  ! Center of the ellipsoid
    type(ECR_t), intent(in) :: Line(2) ! C + t U, two vectors
    real(rg), intent(out) :: T ! T-value of nearest point
    real(rg), intent(out) :: R ! R < 1 => C + t U is an intersection
    type(ECR_t), intent(out), optional :: Grad ! If R >= 1, gradient to the
                               ! Geoid that intersects Line(1) + t * Line(2),
                               ! else undefined
    type(ECR_t), intent(out), optional :: Intersect ! If R >= 1, where G
                               ! intersects the Geoid, i.e., the point on
                               ! the Geoid nearest to Line, else undefined
    type(ECR_t) :: D ! Line(1) - Center
    type(ECR_t) :: N(2)        ! Normal to the ellipsoid
    type(ECR_t) :: PM, UM, VM
    real(rg), allocatable :: T_int(:) ! Intersection of N with the Geoid
    d = line(1) - center
    UM = ECR_t( line(2)%xyz / axes )
    VM = ECR_t( d%xyz / axes )
    t = - ( UM .dot. VM ) / ( UM .dot. UM )
    PM = ECR_t( ( line(1)%xyz + t * line(2)%xyz ) / axes )
    r = sqrt(PM .dot. PM)
    if ( r >= 1 .and. ( present(intersect) .or. present(grad) ) ) then
      d = d + t * line(2)
      n(2) = ECR_t ( d%xyz / axes**2 ) ! Gradient to Geoid
      if ( present(grad) ) grad = n(2)
      if ( present(intersect) ) then
        n(1) = line(1) + t * line(2)
        call line_and_ellipsoid ( axes, center, n, t=t_int )
        intersect = n(1) + t_int(minloc(abs(t_int),1)) * n(2)
      end if
    end if
  end subroutine Line_Nearest_Ellipsoid_RG

  ! =====     Private Procedure     ====================================

  ! Hopefully inlined above
  subroutine Solve ( A2, A1, A0, Line, Intersections, T )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: A2, A1, A0
    type(ECR_t), intent(in) :: Line(2)
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: T(:)
    real(rg) :: D
    real(rg) :: MyT(2)
    integer :: N
    d = a1**2 - 4.0 * a2 * a0
    if ( d < 0.0 ) then
      n = 0
    else if ( d == 0.0 ) then
      n = 1
      myT(1) = -0.5 * a1 / a2
    else
      n = 2
      d = sqrt(d)
      myT = [ 0.5 *(- a1 + d) / a2, 0.5 *(- a1 - d) / a2 ]
    end if
    if ( present(intersections) ) then
      allocate ( intersections(n) )
      intersections = line(1) + myT(:n) * line(2)
    end if
    if ( present(t) ) then
      allocate ( t(n) )
      t = myT(:n)
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
! Revision 2.4  2016/03/02 02:25:06  vsnyder
! Add Ellipsoid_Gradient and Line_Nearest_Ellipsoid
!
! Revision 2.3  2016/02/26 01:59:07  vsnyder
! Add Center argument
!
! Revision 2.2  2016/02/25 21:02:24  vsnyder
! Use .DOT. operator bound to ECR_t
!
! Revision 2.1  2016/02/24 01:19:50  vsnyder
! Initial commit
!
