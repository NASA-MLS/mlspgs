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
module Line_And_Cone_m
!=============================================================================

  private
  public :: Line_And_Cone

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Line_And_Cone ( Lat, Line, Intersections, S )

    ! Compute the intersections of Line with a cone at latitude Lat, its
    ! vertex at the center of the Earth if Lat is geocentric, and its axis
    ! collinear with the Earth's axis.
    ! The Line is given by a vector to a point on it and a vector along it.
    ! Line and Intersections are both given in ECR coordinates.

    ! See wvs-129 for a derivation.

    use Constants, only: Deg2Rad
    use Earth_Constants, only: A => EarthRadA, E2 => Eccentricity_Sq
    use Geolocation_0, only: ECR_t, GeodLat_t, Lat_t, RG

    class(lat_t), intent(in) :: Lat    ! Latitude in degrees.  Assumed to be
                                       ! geocentric unless the dynamic type
                                       ! is GeodLat_t.
    type(ECR_t), intent(in) :: Line(2) ! The line is of the form
                                       ! Line(1) + s * Line(2), i.e., Line(1)
                                       ! is a vector to a point on the line,
                                       ! and Line(2) is a vector along the
                                       ! line.
    type(ECR_t), intent(out), allocatable, optional :: Intersections(:) ! 0 <= size <= 2
    real(rg), intent(out), allocatable, optional :: S(:) ! S-value for intersections

    real(rg) :: A0, A1, A2    ! Polynomial coefficients
    real(rg) :: C(3)          ! line(1)%xyz - [ 0, 0, v ]
    real(rg) :: D             ! Discriminant: a1**2 - a0 * a2
    real(rg) :: MyS(2)        ! Root of a2 * s**2 + 2 * a1 * s + a0
    integer :: N              ! How many intersections?  0 ... 2
    real(rg) :: Sl            ! Sin(Lat) or Sin(Lat)**2
    real(rg) :: V             ! Z-coordinate of vertex, zero for geocentric lat
    real(rg) :: XYZ(3,2)

    sl = sin(deg2rad*lat%d)
    v = 0
    c = line(1)%xyz
    c(3) = c(3) - v           ! C is now C - U
    select type ( lat )
    class is ( geodLat_t )
      v = ( a * e2 * sl ) / sqrt(1 - e2 * sl**2 )
    end select
    sl = sl**2

    !{ The line is $\mathbf{C} + t\, \mathbf{U}$ where $\mathbf{C}$ is
    !  line(1) and $\mathbf{U}$ is line(2).  Intersections are where $t$
    !  is a root of $a_2\,t^2 + 2 \, a_1 \, t + a_0 = 0$, with
    !  \begin{equation*}\begin{split}
    !  a_2 = \,& u_3^2 - | \mathbf{U} |^2 \sin^2 \theta \\
    !  a_1 = \,& (c_3 - v ) \, u_3 -
    !            ( \mathbf{C} - \mathbf{V} ) \cdot \mathbf{U}\, \sin^2 \theta \\
    !  a_0 = \,& ( c_3 -v )^2 -
    !            | \mathbf{C} - \mathbf{V} |^2 \, \sin^2 \theta \\
    !  \end{split}\end{equation*}
    
    a2 = line(2)%xyz(3)**2 - sl * dot_product(line(2)%xyz,line(2)%xyz)
    a1 = c(3) * line(2)%xyz(3) - sl * dot_product(c, line(2)%xyz)
    a0 = c(3)**2 - sl * dot_product(c,c)

    if ( a2 == 0 ) then ! Zero or one intersections
      if ( a1 == 0 ) then ! No intersections
        n = 0
      else
        n = 1
        myS(1) = -a0 / a1
        xyz(:,1) = line(1)%xyz + myS(1) * line(2)%xyz
        if ( xyz(3,1) * lat%d < 0 ) n = 0 ! Wrong hemisphere
      end if
    else
      d = a1**2 - a0 * a2
      if ( d < 0 ) then ! No intersections
        n = 0
      else if ( d == 0 ) then ! One intersection
        n = 1
        myS(1) = -a1 / a2
        xyz(:,1) = line(1)%xyz + myS(1) * line(2)%xyz
        if ( xyz(3,1) * lat%d < 0 ) n = 0 ! Wrong hemisphere
      else
        n = 2
        d = sqrt(d)
        myS(2) = ( -a1 + d ) / a2
        xyz(:,2) = line(1)%xyz + myS(2) * line(2)%xyz
        if ( xyz(3,1) * lat%d < 0 ) n = 1 ! Wrong hemisphere
        myS(1) = ( -a1 - d ) / a2
        xyz(:,1) = line(1)%xyz + myS(1) * line(2)%xyz
        if ( xyz(3,1) * lat%d < 0 ) then ! Wrong hemisphere
          if ( n == 1 ) then
            n = 0
          else
            n = 1
            xyz(:,1) = xyz(:,2)
          end if
        end if
      end if
    end if
    if ( present(intersections) ) then
      allocate ( intersections(n) )
      if ( n > 0 ) intersections(1)%xyz = xyz(:,1)
      if ( n == 2 ) intersections(2)%xyz = xyz(:,2)
    end if
    if ( present(s) ) then
      allocate ( s(n) )
      s = myS(:n)
    end if

  end subroutine Line_And_Cone

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Line_And_Cone_m

! $Log$
! Revision 2.3  2017/10/31 23:44:13  vsnyder
! Add TeXnicalities
!
! Revision 2.2  2016/03/03 21:39:38  vsnyder
! Change name of T argument to S (because it's arc length)
!
! Revision 2.1  2016/02/05 03:38:42  vsnyder
! Initial commit
!
