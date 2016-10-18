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
  !  $c$, and a line given by $\mathbf{C} + s\, \mathbf{U}$, where
  !  $\mathbf{C}$ is a point on the line and $\mathbf{U}$ is a vector along
  !  the line, the intersections occur on the line at values of $s$ such that
  !%
  !  \begin{equation}
  !    (\mathbf{U}\, \mathbf{M})^T (\mathbf{U}\, \mathbf{M})\, s^2 +
  !    (\mathbf{V} \mathbf{M} )^T (\mathbf{U}\, \mathbf{M})\, s +
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
  !    \mathbf{U} \cdot \mathbf{U} \, s^2 +
  !    \mathbf{V} \cdot \mathbf{U} \, s +
  !    \mathbf{V} \cdot \mathbf{V} = r^2 \,.
  !  \end{equation}
  !
  !  See wvs-131 for details of the derivation.
  !
  !  Compute the point on a line that is nearest to an ellipsoid.
  !  This is an approximation developed by replacing the usual equation
  !  for an ellipsoid by $(\mathbf{M} \mathbf{X})^T \mathbf{M} \mathbf{X} = r$
  !  and then solving for a value of $r$ such that the line intersects the
  !  ellipsoid in just one place.  For $|r-1|$ very small, this shouldn't be
  !  very far from the true tangent point.
  !  This occurs at
  !%
  !  \begin{equation}
  !   s =
  ! -\frac{\left(\mathbf{M}\mathbf{V}\right)^T
  !        \left(\mathbf{M}\mathbf{U}\right)}
  !       {\left(\mathbf{M}\mathbf{U}\right)^T
  !        \left(\mathbf{M}\mathbf{U}\right)} =
  ! -\frac{\frac{x_1 x_2}{a^2} + \frac{y_1 y_2}{b^2} + \frac{z_1 z_2}{c^2}}
  !       {\frac{x_2^2}{a^2} + \frac{y_2^2}{b^2} + \frac{z_2^2}{c^2}} \,.
  !  \end{equation}
  !
  !  The position on the ellipsoid that is nearest to the line is
  !%
  !  \begin{equation}
  !  \mathbf{P} = \frac{\mathbf{M}^{-T} \mathbf{M}^{-1} \mathbf{R}}{\ell} =
  !  \frac{(a^2 \mathbf{R}_x, b^2 \mathbf{R}_y, c^2 \mathbf{R}_z)^T}
  !       {\ell}
  !  \end{equation}
  !%
  !  where $\ell = \sqrt{\mathbf{R} \mathbf{M}^{-T} \mathbf{M}^{-1}
  !                            \mathbf{R}} =
  !        \sqrt{a^2 \mathbf{R}_x^2 + b^2 \mathbf{R}_y^2 + 
  !              c^2 \mathbf{R}_z^2}$, and
  !  $\mathbf{R} = \mathbf{C} - ( \mathbf{C} \cdot \mathbf{U} )
  !  \mathbf{U}$ is the vector from the origin to the point on the line
  !  that is nearest to the origin.  The tangent height $h = |\mathbf{R}| -
  !  \ell/|\mathbf{R}|$.  Let $F(\mathbf{X}) =
  !  (\mathbf{M} \mathbf{X})^T \mathbf{M} \mathbf{X} = 1$ define the
  !  ellipsoid.  Then $\nabla \mathbf{F}(\mathbf{P}) =
  !  \nabla F(\mathbf{X})|_{\mathbf{X}=\mathbf{P}} =
  !  2 \, \mathbf{M}^T \mathbf{M} \mathbf{P}$ is normal to the ellipsoid at
  !  $\mathbf{P}$, $\mathbf{T} = \mathbf{P} + h \nabla F(\mathbf{P})$ is
  !  the vector from $\mathbf{P}$ to the tangent point, and the position of
  !  the tangent point is determined by
  !  $s = ( \mathbf{T} - \mathbf{C} ) \cdot \mathbf{U}$.
  !
  !  See wvs-030 for details of the derivations.

  ! The procedures with "Earth Geoid" in their name actually refer to the
  ! Earth Reference Ellipsoid, with axes specified in the Earth_Constants
  ! module.

  implicit NONE
  private

  ! Generic identifiers
  public :: Exact_Line_Nearest_Ellipsoid, Line_And_Ellipsoid
  public :: Line_And_Sphere, Line_Nearest_Ellipsoid

  ! Specific procedures.  The suffix RG refers to the kind type parameter
  ! RG from the module Geolocation_0, meaning "Real kind for geolocations."
  public :: Exact_Line_Nearest_Earth_Geoid
  public :: Exact_Line_Nearest_Ellipsoid_RG
  public :: Line_And_Earth_Geoid
  public :: Line_And_Earth_Geoid_at_H
  public :: Line_And_Ellipsoid_RG
  public :: Line_And_Sphere_RG
  public :: Line_Nearest_Earth_Geoid
  public :: Line_Nearest_Ellipsoid_RG

  interface Exact_Line_Nearest_Ellipsoid
    module procedure Exact_Line_Nearest_Earth_Geoid
    module procedure Exact_Line_Nearest_Ellipsoid_RG
  end interface

  interface Line_And_Ellipsoid
    module procedure Line_And_Earth_Geoid
    module procedure Line_And_Earth_Geoid_at_H
    module procedure Line_And_Ellipsoid_RG
  end interface

  interface Line_And_Sphere
    module procedure Line_And_Sphere_RG
  end interface

  interface Line_Nearest_Ellipsoid
    module procedure Line_Nearest_Earth_Geoid
    module procedure Line_Nearest_Ellipsoid_RG
  end interface

  logical, private, parameter :: Debug = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Line_And_Earth_Geoid ( Line, Intersections, S )
    ! Find intersections of Line with the Earth reference ellipsoid.
    ! See wvs-131.
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: S(:) ! S-values of intersections
    real(rg) :: A2, A1, A0
    type(ECR_t) :: UM, VM
    line(2) = line(2) / line(2)%norm2()
    UM = ECR_t( line(2)%xyz / real([ A, A, B ],rg) )
    VM = ECR_t( line(1)%xyz / real([ A, A, B ],rg) )
    a2 = UM .dot. UM
    a1 = VM .dot. UM
    a0 = ( VM .dot. VM ) - 1.0
    call solve ( a2, a1, a0, line, intersections, s )
  end subroutine Line_And_Earth_Geoid

  subroutine Line_And_Earth_Geoid_at_H ( Line, H, Tol, Intersections, S )
    ! Find intersections of Line with a surface at height H above the
    ! Earth reference ellipsoid.  See wvs-134.
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, H_V_Geod, RG
    use Zero_m, only: Zero
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    real(rg), intent(in) :: H          ! Height above Earth reference ellipsoid
    real(rg), intent(in) :: Tol        ! Tolerance criterion for height
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: S(:) ! S-values of intersections
    real(rg) :: Alpha                  ! Cosine of angle between Line(2) and Grad
    real(rg) :: F1, F2                 ! For Zero, q.v.
    type(h_v_geod) :: Geod             ! Geodetic coordinates of candidate solution
    type(ECR_t) :: Grad                ! Gradient to ellipsoid
    real(rg) :: Hgt                    ! at point nearest, if no intersection
    integer :: I
    integer :: Mode(2)                 ! For Zero, q.v.
    type(ECR_t), allocatable :: MyInts(:)
    type(ECR_t) :: MyLine(2)
    real(rg), allocatable :: MyS(:)
    integer :: N, M                    ! Number of intersections
    logical :: NeedGrad
    type(ECR_t) :: Point               ! Line(1) + s * Line(2)
    real(rg) :: R                      ! Output from Line_Nearest_Ellipsoid
    real(rg) :: X1, X2                 ! For Zero, q.v.
    ! Is there an intersection with an ellipse with semi-minor axes
    ! [a+h,a+h,b+h]?
    line(2) = line(2) / line(2)%norm2()
    myLine = [line(1),line(2)]
    call line_and_ellipsoid ( real([a+h,a+h,b+h],rg), myLine, myInts, myS )
    n = size(myS)
    if ( n == 0 ) then
      if ( debug ) print '(a)', 'No intersection with [a+h,a+h,b+h]'
      ! No intersection.  How close is the closest approach?
      deallocate ( myInts, myS )
      allocate ( myInts(1), myS(1) )
      call Line_Nearest_Ellipsoid ( myLine, myS(1), r, grad, myInts(1), hgt )
      if ( hgt > h ) then
        ! Closest approach is too high, therefore no intersection.
        if ( debug ) print '(a,f8.0,a)', 'Closest approach at H =', hgt, ' is too high'
        if ( present(intersections) ) allocate ( intersections(0) )
        if ( present(s) ) allocate ( intersections(0) )
        return
      end if
      if ( debug ) print '(a,f8.0)', 'Closest approach at H =', hgt
      n = 1
      needGrad = .false.
    else
      needGrad = .true.
    end if
    if ( debug ) then
      print '(i0,a,1p,2g15.8)', n, ' intersections with [a+h,a+h,b+h] at S =', myS
      print '(a,2(" [",3g15.8,"]":))', 'Geodetic coordinates', myInts%geod()
      print '(a,2(" [",3g15.8,"]":))', 'ECR coordinates     ', myInts
    end if
    do i = 1, n
      x2 = myS(i)
      point = line(1) + x2 * line(2)
      geod = point%geod()
      f2 = h - geod%v
      if ( abs(f2) < abs(tol) ) then
        if ( debug ) print '(4(a,g17.10))', 'X2 = ', x2, ' F2 =', f2, &
          & ' H2 =', geod%v
        cycle
      end if
      if ( needGrad ) then
      ! grad = point / ( real([ A**2, A**2, B**2 ], rg) ) ! ifort 16.0.2 can't do this
        grad = ECR_t ( point%xyz / ( real([ A**2, A**2, B**2 ], rg) ) )
      end if
      needGrad = .true. ! needed next time, if there is a next time
      grad = grad / grad%norm2()
      ! Cosine of angle between line and gradient
      alpha = line(2) .DOT. grad
      if ( debug ) print '(4(a,g17.10))', 'X2 = ', x2, ' F2 =', f2, &
        & ' H2 =', geod%v, ' Alpha =', alpha
      x1 = x2 + ( h - geod%v ) / sign(max(abs(alpha),0.5_rg),alpha)
      mode(i) = 0
      do while ( mode(i) < 2 )
        point = line(1) + x1 * line(2)
        geod = point%geod()
        f1 = h - geod%v
        if ( debug ) print '(3(a,g17.10))', 'X1 = ', x1, ' F1 =', f1, ' H1 =', geod%v
        call zero ( x1, f1, x2, f2, mode(i), -abs(tol) ) ! Tolerance is on F not X
      end do
      myS(i) = x1
    end do
    m = count(mode<4) ! How many normal terminations of Zero
    if ( present(intersections) ) allocate ( intersections(m) )
    if ( present(s) ) allocate ( s(m) )
    m = 0
    do i = 1, n
      if ( mode(i) < 4 ) then
        m = m + 1
        if ( present(intersections) ) intersections(m) = line(1) + myS(i) * line(2)
        if ( present(s) ) s(m) = myS(i)
      end if
    end do
  end subroutine Line_And_Earth_Geoid_at_H

  subroutine Line_And_Ellipsoid_RG ( Axes, Line, Intersections, S, Center )
    ! Intersection of a line with arbitrary ellipsoid, not necessarily
    ! the Earth Reference Ellipsoid.  See wvs-131.
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(inout) :: Line(2) ! V + s U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V; U is made
                                       ! an unit vector here
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: S(:) ! S-values of intersections
    type(ECR_t), intent(in), optional :: Center  ! Center of the ellipsoid
    real(rg) :: A2, A1, A0
    type(ECR_t) :: D ! Line(1) - Center
    type(ECR_t) :: UM, VM
    if ( present(center) ) then
      d = line(1) - center
      VM = ECR_t( d%xyz / axes )
    else
      VM = ECR_t( line(1)%xyz / axes )
    end if
    UM = ECR_t( Line(2)%xyz / axes )
    a2 = UM .dot. UM
    a1 = VM .dot. UM
    a0 = ( VM .dot. VM ) - 1.0
    call solve ( a2, a1, a0, line, intersections, s )
  end subroutine Line_And_Ellipsoid_RG

  subroutine Line_And_Sphere_RG ( Radius, Line, Intersections, S, Center )
    ! Intersection of a line with a sphere.
    ! See wvs-131.
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Radius     ! In same units as Center
    type(ECR_t), intent(inout) :: Line(2) ! V + s U, two vectors, where V is a
                                       ! vector from the center of the ellipsoid
                                       ! to a point on the line, and U is a
                                       ! vector along the line at V; U might be
                                       ! made an unit vector by Solve
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: S(:) ! S-values of intersections
    type(ECR_t), intent(in), optional :: Center  ! Center of the sphere
    real(rg) :: A2, A1, A0
    type(ECR_t) :: D ! Line(1) - Center
    if ( present(center) ) then
      d = line(1) - center
    else
      d = line(1)
    end if
    a2 = line(2) .dot. line(2)
    a1 = d .dot. line(2)
    a0 = ( d .dot. d ) - radius**2
    call solve ( a2, a1, a0, line, intersections, s )
  end subroutine Line_And_Sphere_RG

  subroutine Exact_Line_Nearest_Earth_Geoid ( Line, S, P, H )
    ! Find the point  on a line that is nearest to the Earth
    ! geoid.  See wvs-030 for derivation.
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, RG
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    real(rg), intent(out), optional :: S  ! S-value of nearest point
    type(ECR_t), intent(out), optional :: P ! Point on Earth Geoid nearest
                                          ! to C + s U
    real(rg), intent(out), optional :: H  ! Tangent height

    real(rg), parameter :: Axes(3) = [ A, A, B ]
    real(rg) :: D      ! L / |R|
    real(rg) :: L      ! sqrt( [ A**2, A**2, B**2 ] .dot. R**2 )
    real(rg) :: MyH
    type(ECR_t) :: MyP
    type(ECR_t) :: R   ! Point on Line nearest the origin
    real(rg) :: Rnorm  ! |R|

    line(2) = line(2) / line(2)%norm2() ! U
    ! R = C - ( C . U ) U
    r = line(1) - ( line(1) .dot. line(2) ) * line(2)
    rnorm = r%norm2()
    l = sqrt ( dot_product(axes**2, r%xyz**2 ) )
    d = l / rnorm
    myH = rnorm - d
    myP = ecr_t ( axes**2 * r%xyz ) / l
    ! ( P + h * \nabla F(P) - C ) . U; grad_geoid() computes the unit gradient
    if ( present(s) ) s = ( myP + myH * myP%grad_geoid() - line(1) ) .dot. line(2)
    if ( present(h) ) h = myH
    if ( present(p) ) p = myP
  end subroutine Exact_Line_Nearest_Earth_Geoid

  subroutine Exact_Line_Nearest_Ellipsoid_RG ( Axes, Line, S, P, H )
    ! Find the point  on a line that is nearest to an arbitrary ellipsoid
    ! with its center at the origin.  See wvs-030 for derivation.
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)       ! Semi-minor axes in same units as Center
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    real(rg), intent(out), optional :: S  ! S-value of nearest point
    type(ECR_t), intent(out), optional :: P ! Point on Earth Geoid nearest
                                          ! to C + s U
    real(rg), intent(out), optional :: H  ! Tangent height

    real(rg) :: D      ! L / |R|
    real(rg) :: L      ! sqrt( [ A**2, A**2, B**2 ] .dot. R**2 )
    real(rg) :: MyH
    type(ECR_t) :: MyP
    type(ECR_t) :: R   ! Point on Line nearest the origin
    real(rg) :: Rnorm  ! |R|

    line(2)  = line(2) / line(2)%norm2()
    ! R = C - ( C . U ) U
    r = line(1) - ( line(1) .dot. line(2) ) * line(2)
    rnorm = r%norm2()
    l = sqrt ( dot_product(axes**2, r%xyz**2 ) )
    d = l / rnorm
    myH = rnorm - d
    myP = ecr_t ( axes**2 * r%xyz ) / l
    ! ( P + h * \nabla F(P) - C ) . U; grad_geoid() computes the unit gradient
    if ( present(s) ) s = ( myP + myH * myP%grad_geoid() - line(1) ) .dot. line(2)
    if ( present(h) ) h = myH
    if ( present(p) ) p = myP
  end subroutine Exact_Line_Nearest_Ellipsoid_RG

  subroutine Line_Nearest_Earth_Geoid ( Line, S, R, Grad, Intersect, H )
    ! Find a point (approximately) on a line that is nearest to the Earth
    ! geoid.  See wvs-030 for derivation.
    use Earth_Constants, only: A => EarthRadA, B => EarthRadB
    use Geolocation_0, only: ECR_t, Norm2, RG
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    real(rg), intent(out) :: S ! S-value of nearest point
    real(rg), intent(out) :: R ! R < 1 => C + s U is an intersection
    type(ECR_t), intent(out), optional :: Grad ! If R >= 1, gradient to the
                               ! Geoid that intersects Line(1) + s * Line(2),
                               ! else undefined
    type(ECR_t), intent(out), optional :: Intersect ! If R >= 1, where G
                               ! intersects the Geoid, i.e., the point on
                               ! the Geoid nearest to Line, else undefined
    real(rg), intent(out), optional :: H ! Distance at nearest point, zero if
                               ! there is an intersection ( R < 1 ).
    type(ECR_t) :: Int         ! Intersection of gradient with ellipsoid
    type(ECR_t) :: N(2)        ! Normal to the ellipsoid
    type(ECR_t) :: PM, UM, VM
    real(rg), allocatable :: S_int(:) ! Intersection of N with the Geoid
    line(2) = line(2) / line(2)%norm2()
    UM = ECR_t( line(2)%xyz / real([ A, A, B ],rg) )
    VM = ECR_t( line(1)%xyz / real([ A, A, B ],rg) )
    s = - ( UM .dot. VM ) / ( UM .dot. UM )
    PM = ECR_t( ( line(1)%xyz + s * line(2)%xyz ) / real([ A, A, B ],rg) )
    r = sqrt(PM .dot. PM)
    if ( r >= 1 .and. &
       & ( present(intersect) .or. present(grad) .or. present(h) ) ) then
      n(2) = ECR_t ( 2.0 * PM%xyz / real([ A, A, B ],rg) ) ! Gradient to Geoid
      if ( present(grad) ) grad = n(2)
      if ( present(intersect) .or. present(h) ) then
        n(1) = line(1) + s * line(2)
        call line_and_ellipsoid ( n, s=s_int )
        int = n(1) + s_int(minloc(abs(s_int),1)) * n(2)
        if ( present(intersect) ) intersect = int
        if ( present(h) ) h = norm2 ( ( line(1) + s * line(2) ) - int )
      end if
    else if ( present(h) ) then
      h = 0
    end if
  end subroutine Line_Nearest_Earth_Geoid

  subroutine Line_Nearest_Ellipsoid_RG ( Axes, Line, S, R, Grad, Intersect, Center )
    ! Find a point (approximately) on a line that is nearest to an
    ! arbitrary ellipsoid.  See wvs-030 for derivation.
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: Axes(3)    ! Semi-minor axes in same units as Center
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U is made
                                          ! an unit vector here
    real(rg), intent(out) :: S ! S-value of nearest point
    real(rg), intent(out) :: R ! R < 1 => C + s U is an intersection
    type(ECR_t), intent(out), optional :: Grad ! If R >= 1, gradient to the
                               ! Geoid that intersects Line(1) + s * Line(2),
                               ! else undefined
    type(ECR_t), intent(out), optional :: Intersect ! If R >= 1, where G
                               ! intersects the Geoid, i.e., the point on
                               ! the Geoid nearest to Line, else undefined
    type(ECR_t), intent(in), optional :: Center  ! Center of the ellipsoid
    type(ECR_t) :: D ! Line(1) - Center
    type(ECR_t) :: N(2)        ! Normal to the ellipsoid
    type(ECR_t) :: PM, UM, VM
    real(rg), allocatable :: S_int(:) ! Intersection of N with the Geoid
    if ( present(center) ) then
      d = line(1) - center
    else
      d = line(1)
    end if
    line(2) = line(2) / line(2)%norm2()
    UM = ECR_t( line(2)%xyz / axes )
    VM = ECR_t( d%xyz / axes )
    s = - ( UM .dot. VM ) / ( UM .dot. UM )
    PM = ECR_t( ( line(1)%xyz + s * line(2)%xyz ) / axes )
    r = sqrt(PM .dot. PM)
    if ( r >= 1 .and. ( present(intersect) .or. present(grad) ) ) then
      d = d + s * line(2)
      n(2) = ECR_t ( 2.0 * d%xyz / axes**2 ) ! Gradient to Geoid
      if ( present(grad) ) grad = n(2)
      if ( present(intersect) ) then
        n(1) = line(1) + s * line(2)
        call line_and_ellipsoid ( axes, n, s=s_int, center=center )
        intersect = n(1) + s_int(minloc(abs(s_int),1)) * n(2)
      end if
    end if
  end subroutine Line_Nearest_Ellipsoid_RG

  ! =====     Private Procedure     ====================================

  ! Solve a2 s^2 + 2 a1 s + a0 = 0.  Then compute Line(1) + S * Line(2).
  ! This is hopefully inlined above.
  pure subroutine Solve ( A2, A1, A0, Line, Intersections, S )
    use Geolocation_0, only: ECR_t, RG
    real(rg), intent(in) :: A2, A1, A0
    type(ECR_t), intent(inout) :: Line(2) ! C + s U, two vectors, U might be
                                          ! made an unit vector here
    type(ECR_t), intent(out), allocatable, optional :: Intersections (:) ! 0..2 elements
    real(rg), intent(out), allocatable, optional :: S(:)
    real(rg) :: D
integer :: I ! ifort 14.0.0 doesn't like elemental * operator bound to ECR_t
             ! ifort 15.0.2 and ifort 16.0.2 accept it
    real(rg) :: MyS(2)
    integer :: N
    d = a1**2 - a2 * a0
    if ( d < 0.0 ) then
      n = 0
    else if ( d == 0.0 ) then
      n = 1
      myS(1) = -a1 / a2
    else
      n = 2
      d = sqrt(d)
      myS = [ (- a1 + d) / a2, (- a1 - d) / a2 ]
    end if
    if ( present(intersections) ) then
      allocate ( intersections(n) )
        line(2) = line(2) / line(2)%norm2()
!       intersections(:n) = line(1) + myS(:n) * line(2)
do i = 1, n
intersections(i) = line(1) + myS(i) * line(2)
end do
    end if
    if ( present(s) ) then
      allocate ( s(n) )
      s = myS(:n)
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
! Revision 2.11  2016/10/18 00:33:54  vsnyder
! Make line(2) an unit vector, correct ellipsoid grad
!
! Revision 2.10  2016/09/27 00:55:17  vsnyder
! Repair computation of S in Exact_Line_Nearest...
!
! Revision 2.9  2016/05/26 19:24:28  vsnyder
! Correct height calculation in Exact_Line_Nearest_Ellipsoid
!
! Revision 2.8  2016/05/24 01:18:55  vsnyder
! Change "t" to "s" because it's path length.  Add exact non-iterative
! computation of tangent point.  Add TeXnicalities.
!
! Revision 2.7  2016/03/25 00:45:42  vsnyder
! Make specific procedures public.  Add a Debug parameter and some
! debugging output.  Add Line_And_Earth_Geoid_at_H to compute intersection
! of a line with a surface at a specified height above the Earth Reference
! Ellipsoid.  Convert the kind of the semi-axes to RG so that type-bound
! divide for ECR_t will work.  Make Center argument optional and improve
! its handling in several places.  Add H argument to
! Line_Nearest_Earth_Geoid to report height at nearest point.  Correct
! Solve subroutine.
!
! Revision 2.6  2016/03/03 21:39:38  vsnyder
! Change name of T argument to S (because it's arc length)
!
! Revision 2.5  2016/03/02 21:48:20  vsnyder
! Move ellipsoid gradient code to Geolocation_0
!
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
