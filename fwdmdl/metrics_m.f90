! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Metrics_m

  use Constants, only: Rad2Deg
  use MLSKinds, only: RP

  implicit NONE
  private
  public :: Height_Metrics, More_Metrics, Solve_H_Phi, Tangent_Metrics
  public :: More_Points

  ! Values for Stat argument to/from Solve_H_Phi:
  integer, parameter, public :: No_sol = 0  ! No solution
  integer, parameter, public :: NaN_sol = 1 ! No solution -- NaN
  integer, parameter, public :: Order = 2   ! No solution -- Phi out of order
  integer, parameter, public :: Complex = 3 ! Complex start for Newton iteration
                                            ! Good, Grid1, Grid2 must be last three
  integer, parameter, public :: Good = 4    ! Newton converged, tangent point, extrapolated
                                            ! Grid1 and Grid2 have to be the last two
  integer, parameter, public :: Grid1 = 5   ! Close to left grid point
  integer, parameter, public :: Grid2 = 6   ! Close to right grid point
  character(5), parameter :: nStat(no_sol:grid2) = &
    & (/ 'none ', 'NaN  ', 'order', 'cmplx', 'good ', 'grid1', 'grid2' /)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Private stuff
  real(rp), parameter :: DefaultTol = 0.001_rp ! km, for Solve_h_phi
  ! To control debugging
  logical, parameter :: Bin42 = .false. ! Write inputs on unit 42 if failure
  logical, parameter :: Debug = .false.
  logical, parameter :: NewtonDetails = .true. .and. debug
  logical, parameter :: ComplexDebug = .true. .and. debug
  ! For debugging output format:
  ! logical, parameter :: Clean = .false.
  character(len=4), parameter :: Options=' '
  real, parameter :: Ang = rad2deg ! degrees
! real, parameter :: Ang = 1.0     ! .not. degrees
  character(*), parameter :: Fmt = '( i4,i2,i4,f12.5,f10.3,a )'

contains

  ! --------------------------------------------  Tangent_Metrics  -----
                               ! Inputs
  subroutine Tangent_Metrics ( Phi_T, P_Basis, H_Ref, Tan_Ind_F, &
                               ! Outputs
    &                          H_Surf, H_Tan, &
                               ! Optional inputs
    &                          Z_Ref, Tan_press, Surf_temp, Surf_height )

  ! Compute the surface height and the tangent height at Phi_T.

    use MLSKinds, only: RP
    use Sparse_Eta_m, only: Sparse_Eta_t

    ! inputs:

    real(rp), intent(in) :: Phi_t      ! orbit projected tangent geodetic angle,
                                       ! radians
    real(rp), intent(in) :: P_basis(:) ! horizontal temperature representation
                                       ! basis, radians
    real(rp), intent(in) :: H_ref(:,:) ! heights by fine zeta and p_basis,
                                       ! referenced to equivalent circular earth
                                       ! center, km
    integer, intent(in) :: Tan_ind_f   ! Index of tangent in H_Ref, which is on
                                       ! fine Zeta grid

    ! outputs:
    real(rp), intent(out) :: H_Surf    ! Height of the pressure reference
                                       ! surface at Phi_t -- interpolated in
                                       ! Surf_Height if present(Surf_Height)
                                       ! else interpolated in row 1 of H_REF,
                                       ! referenced to equivalent circular earth
                                       ! center, km
    real(rp), intent(out) :: H_Tan     ! Tangent height above H_Surf (negative
    !                                  ! for Earth-intersecting rays),
                                       ! referenced to equivalent circular earth
                                       ! center, km

    ! optional inputs
    ! We have Z_Ref, Tan_Press and Surf_Temp if the pointing is below the
    ! surface index (where zeta Z_Ref) and we don't have Surf_Height.
    ! Otherwise we have Surf_Height.
    real(rp), optional, intent(in) :: Z_Ref        ! Zeta corresponding to H_Ref(1,:)
    real(rp), optional, intent(in) :: Tan_press    ! Tangent pressure
    real(rp), optional, intent(in) :: Surf_temp(:) ! Surface temperature at p_basis
    real(rp), optional, intent(in) :: Surf_height(:) ! Surface height in km
                                       ! above mean sea level (whatever that
                                       ! means) on P_Basis.

    type(sparse_eta_t) :: ETA_T        ! Interpolating coefficients

    ! Get interpolating coefficients (eta_t) from p_basis to phi_t
    call eta_t%eta ( p_basis, phi_t )

    if ( present(surf_height) ) then
      ! We set the surface reference at the actual surface height if we
      ! have it, and adjust r_eq and h_tan relative to this, and adjust
      ! h_ref accordingly.
      h_surf = eta_t%row_dot_vec ( 1, surf_height )
    else
      ! If we don't have the actual surface height, we set the surface
      ! reference at the input z_ref and adjust r_eq and h_tan relative
      ! to this, and adjust h_ref accordingly.
      h_surf = eta_t%row_dot_vec ( 1, h_ref(1,:) )
    end if

    ! compute the tangent height above H_surf.
    if ( present(tan_press) .and. present(surf_temp) .and. &
      &  .not. present(surf_height) ) then
      ! Earth intersecting ray. Compute GP height (km) of tangent pressure
      ! below surface. This will be negative because tan_press < z_ref.
      ! present(tan_press) requires present(surf_temp).  We don't need to
      ! subtract h_surf here because this gives km from the z_ref surface.
      h_tan = eta_t%row_dot_vec ( 1, surf_temp ) * (tan_press-z_ref)/14.8
    else
      h_tan = eta_t%row_dot_vec ( 1, h_ref(tan_ind_f,:) ) - h_surf
    end if

  end subroutine Tangent_Metrics

  ! ---------------------------------------------  Height_Metrics  -----

                            ! Inputs:
  subroutine Height_Metrics (  phi_t, tan_ind, p_basis, h_ref, h_surf, &
                            &  tan_ht_s, z_ref, r_eq,                  &
                            ! Outputs:
                            &  req_s, vert_inds, h_path, p_path,       &
                            ! Optional inputs:
                            &  h_tol, forward )

    !{ This subroutine computes {\tt h\_path} and {\tt p\_path} that define
    !  a 2-d integration path.  These points are at the intersections of the
    !  line of sight given by $H = H_\text{tan} \sec \phi$ and the heights of
    !  surfaces of constant $\zeta$ represented in piecewise linear segments
    !  by consecutive elements of each row of {\tt H\_ref}.

    ! If the ray is an Earth-intersecting ray it is assumed to reflect from
    ! the pressure reference surface, not the geometric Earth surface.

    use Constants, only: Pi, Rad2Deg
    use Dump_0, only: Dump
    use GLNP, only: NG, NGP1
    use IEEE_Arithmetic, only: IEEE_Is_NaN
    use MLSKinds, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMsg_Error, MLSMsg_Warning
    use MLSStringLists, only: SwitchDetail
    use Output_M, only: Output
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: Phi_t      ! Orbit projected tangent geodetic angle
    integer, intent(in) :: Tan_ind     ! Tangent height index, 1 = center of
                                       ! longest path
    real(rp), intent(in) :: P_basis(:) ! Horizontal temperature representation
                                       ! basis, radians
    real(rp), intent(in) :: H_ref(:,:) ! Geodetic heights (above R_Eq) by z_ref
                                       ! and p_basis, km
    real(rp), intent(in) :: H_Surf     ! Height of the pressure reference
                                       ! surface z_ref(1) above R_eq at phi_t, km
    real(rp), intent(in) :: Tan_Ht_s   ! Tangent height above H_Surf -- negative
                                       ! for Earth-intersecting ray, km
    real(rp), intent(in) :: Z_ref(:)   ! -log pressures (zetas) for which
                                       ! heights/temps are needed.  Only used
                                       ! where H/Phi iteration fails.
    real(rp), intent(in) :: R_Eq       ! Equivalent circular earth radius of
                                       ! true surface at Phi_T, km

    ! outputs:
    real(rp), intent(out) :: Req_s     ! H_Surf + R_Eq, Equivalent circular
                                       ! earth radius of pressure reference
                                       ! surface at Phi_T, km
    integer, intent(out) :: Vert_Inds(:) ! What to use in h_ref, 1:n_path
    real(rp), intent(out) :: H_Path(:) ! computed path heights, referenced to
                                       ! equivalent circular Earth center, km
    real(rp), intent(out) :: P_Path(:) ! computed path phi's, radians

    ! optional inputs
    real(rp), optional, intent(in) :: H_Tol ! Height tolerance in kilometers
                                       ! for convergence of phi/h iteration
    logical, optional, intent(in) :: Forward  ! for subsurface rays, indicates
                                       ! xi >= xi_sub (see Generate_TScat in
                                       ! FullForwardModel)

    ! Local variables.
    integer :: Do_Dumps    ! <0 = no dump, >=0 = dump, >0 = stop
    integer :: H_Phi_Dump  ! <0 = no dump, >=0 = dump, >0 = stop
    integer :: I, I1, I2, IBAD, J, K
    logical :: MyForward   ! Default .true., else Forward if it's present
    integer :: No_Bad_Fits
    integer :: No_Grid_Fits
    integer :: N_Path      ! Path length = 2*(n_vert+ngp1-tan_ind)
    integer :: N_Tan       ! Tangent index in path, (N_Path - ngp1 + 1)/2
    integer :: N_Vert      ! Size(H_ref,1)
    integer :: P_Coeffs    ! Size(P_basis)

    integer :: Stat(size(vert_inds))

    logical :: Outside     ! Solve_H_Phi converged outside the allowed area
    real(rp) :: A, B, C    ! Polynomial coefficients used to solve for phi and H
    real(rp) :: DP         ! p_basis(j+1) - p_basis(j)
    real(rp) :: DPJ0       ! p_basis(j  )-phi_t
    real(rp) :: DPJ1       ! p_basis(j+1)-phi_t
    real(rp) :: H          ! Tentative solution for H
    real(rp) :: My_H_Tol   ! Tolerance in kilometers for height convergence
    real(rp) :: P          ! Tentative solution for phi
    real(rp) :: Start      ! Starting phi for Newton iteration if > -0.5*huge(0.0_rp)
    real(rp) :: Tan_Ht     ! tan_ht_s + Req_s
    real(rp) :: Theta      ! 2.0_rp*Acos(tan_ht/req_s), angle between incident
                           ! and reflected ray for subsurface tangent height

    type(sparse_eta_t) :: ETA_T        ! Interpolating coefficients
    real(rp) :: Phi_Offset(size(vert_inds)) ! Orbit projected tangent geodetic
                           ! angle if the ray is not an earth intersecting
                           ! ray.  Otherwise the orbit projected tangent
                           ! geodetic angle of the ray between the instrument
                           ! and the intersection point, then the orbit projected
                           ! tangent geodetic angle of the reflected ray.

    real(rp) :: Phi_Sign(size(vert_inds))   ! +/- 1.0

    real(rp), parameter :: Pix2 = 2.0*pi

    ! For debugging
    logical, parameter :: TestStat = .false. ! Test whether any stat < good

!   It would be nice to do this the first time only, but the
!   retrieve command in the L2CF can now change switches
!   if ( do_dumps < 0 ) then ! First time only
      do_dumps = switchDetail(switches,'metd')
      h_phi_dump = switchDetail(switches,'hphi')
!   end if

    n_path = size(vert_inds)
    n_tan = (n_path - ng ) / 2
    n_vert = size(h_ref,1)

    my_h_tol = defaultTol ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    p_coeffs = size(p_basis)

    ! Compute equivalent earth radius req_s at phi_t of the pressure
    ! reference surface interpolated to phi_t

    req_s = r_eq + h_surf            ! Radius of z_grid(1) surface at phi_t
    tan_ht = tan_ht_s + req_s           ! From equivalent earth center

    stat = no_sol ! No solutions
    ! The tangent points and the zero-thickness layer between them are special
    stat(n_tan:n_tan+ngp1) = good

    ! p_basis and p_path are phi's in offset radians relative to
    ! phi_t, that is, the phi_t, p_basis or p_path = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t

    if ( tan_ht_s < 0.0 ) then ! Earth-intersecting ray
      myForward = .true.
      if ( present(forward) ) myForward = forward
      theta = 2.0_rp*Acos(tan_ht/req_s)
      ! Shouldn't need ever to look at phi_offset(n_tan+1:n_tan+ng)
      if ( myForward ) then
        phi_offset(1:n_tan) = phi_t + theta
        phi_offset(n_tan+1:) = phi_t
      else
        phi_offset(:n_tan+ng) = phi_t
        phi_offset(n_tan+ngp1:) = phi_t - theta
      end if
      ! This expression for p_path puts the midpoint of the path at the
      ! reflection point.
      p_path(n_tan:n_tan+ngp1) = 0.5 * ( phi_offset(n_tan) + phi_offset(n_tan+ngp1) )
      h_path(n_tan:n_tan+ngp1) = req_s
    else
      phi_offset = phi_t
      p_path(n_tan:n_tan+ngp1) = phi_offset(n_tan:n_tan+ngp1)
      h_path(n_tan:n_tan+ngp1) = tan_ht
    end if

    ! Only use the parts of the reference grids that are germane to the
    ! present path

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (tan_ind, i=1,ng), &
                &  (i, i=tan_ind,n_vert) /)

    ! sign of phi vector
    phi_sign = (/ (-1.0_rp, i=n_vert,tan_ind,-1), (1.0_rp, i=1,ng), &
               &  (+1.0_rp, i=tan_ind,n_vert) /)

    if ( h_phi_dump >= 0 ) call dumpInput ( tan_ind )

      if ( debug ) &
        & print '(a)', &
        & '   I J   K     offset   TAN_HT   HREF(K,J) HREF(K,J+1)    A        B    phi sign    C&
        &            D(K,J)        D(K,J+1)   STAT', &
        & '      P          phi    H_PATH   HREF(K,J) HREF(K,J+1)  H(TRIG)    D           N   Diffs'

    no_bad_fits = n_path - ngp1 - 1 ! Exclude coarse tangent points and fine
                                    ! points between them
    no_grid_fits = 0

    ! Get solutions outside of p_basis, if any, assuming constant H
    ! outside of p_basis.  Skip the points, if any, in the zero-thickness
    ! tangent layer.
    do i1 = 1, n_path
      if ( stat(i1) == good ) cycle ! Skip the tangent points
      k = vert_inds(i1)
      h = max(tan_ht,h_ref(k,1)+r_eq) ! or h_ref(k,1) - h_surf + req_s
      p = acos(tan_ht/h) * phi_sign(i1) + phi_offset(i1)
      if ( p >= p_basis(1) ) exit     ! phi is monotone
      p_path(i1) = p
      h_path(i1) = h
      no_bad_fits = no_bad_fits - 1
      stat(i1) = good
        if ( debug ) &
          & print fmt, i1, 0, k, ang*p_path(i1), h_path(i1), '  Below p_basis'

    end do
    do i2 = n_path, 1, -1
      if ( stat(i2) == good ) cycle ! Skip the tangent points
      k = vert_inds(i2)
      h = max(tan_ht,h_ref(k,p_coeffs)+r_eq)
      p = acos(tan_ht/h) * phi_sign(i2) + phi_offset(i2)
      if ( p <= p_basis(p_coeffs) ) exit  ! phi is monotone
      p_path(i2) = p
      h_path(i2) = h
      no_bad_fits = no_bad_fits - 1
      stat(i2) = good
        if ( debug ) &
          & print fmt, i2, p_coeffs, k, ang*p_path(i2), h_path(i2), '  Above p_basis'
    end do
      if ( debug ) then
        print fmt, ( n_tan+j, 0, vert_inds(n_tan+j), ang*p_path(n_tan+j), &
                   & h_path(n_tan+j), '  Tangent point', j = 0, ngp1 )
      end if

    !{ Let $H_{ij}$ denote an element of {\tt H\_ref}.  Remember that all
    ! heights are referenced from $R^\oplus_{\text{eq}}$.  Calculate the
    ! heights by solving for the intersections of the line of sight ($H =
    ! H_{\text{tan}} \sec \phi)$ with constant-$\zeta$ surfaces, piecewise
    ! linear segments of which are represented by two consecutive columns in
    ! each row of {\tt H\_ref}.  See wvs-048.
    !
    ! Extrapolated solutions are not accepted.

    ibad = i2 ! Index of row of first failure in a column
phi:do j = 1, p_coeffs-1
      dp = p_basis(j+1)-p_basis(j)
      a = dp * tan_ht ! Tan_Ht is height from equivalent Earth center
path: do i = i1, i2
        if ( stat(i) == good ) cycle  ! probably the tangent point
!        Don't do this; it causes trouble:
!        if ( stat(i) == grid1 ) cycle ! can't possibly be in this column
        k = vert_inds(i)
        dpj1 = p_basis(j+1)-phi_offset(i)
        dpj0 = p_basis(j)-phi_offset(i)
        b = -(h_ref(k,j+1) - h_ref(k,j)) ! h_surf cancels here
        c = -(h_ref(k,j)-h_surf)   * dpj1 &
          & +(h_ref(k,j+1)-h_surf) * dpj0 &
          & + dp * tan_ht_s ! tan_ht_s is height above H_surf
          if ( newtonDetails ) &
            & print 120, i, j, k, ang*phi_offset(i), tan_ht, &
            & h_ref(k,j)+r_eq, h_ref(k,j+1)+r_eq, &
            & a, b, int(phi_sign(i)), c, &
            & a*(1.0/cos(p_basis(j:j+1)-phi_offset(i))-1.0)+ &
            &         b*(p_basis(j:j+1)-phi_offset(i))+c, &
            & trim(nstat(stat(i)))
          120 format ( i4, i2, i4, f12.5, 2f10.3,f11.3, f10.5, 1p,g14.6, &
            & i3, 3g14.6, 1x, a )
        if ( stat(i) >= grid1 ) then
          no_grid_fits = no_grid_fits - 1
          no_bad_fits = no_bad_fits + 1
        end if
        start = -0.5*huge(0.0_rp)
        if ( i > 1 .and. i < n_path ) then
          if ( stat(i) == complex .and. stat(i-1) == good .and. stat(i+1) == good) then
            ! Complex starting point last time.  Start halfway between two good ones.
            start = 0.5 * ( p_path(i-1) + p_path(i+1) ) ! Solve_H_Phi subtracts phi_offset
            if ( newtonDetails ) &
              print '(2(a,i0),2(a,3f11.5),2(a,f11.5))', 'Phi(',i-1,':',i+1,')', &
                & p_path(i-1:i+1), ' =', p_path(i-1:i+1)*rad2deg, &
                & ' Starting phi ', start, ' = ', start*rad2deg
          end if
        end if
        call Solve_H_Phi ( p_basis(j:j+1), phi_offset(i), phi_sign(i), &
          &                h_ref(k,j:j+1), a, b, c, &
          &                tan_ht, r_eq, my_h_tol, i /= i2 .or. j /= p_coeffs-1, &
          &                start, h_path(i), p_path(i), stat(i), outside )
        if ( IEEE_Is_NaN(p_path(i)) ) then
          stat(i) = NaN_sol
          ! print *, 'i ', i, '  stat(i) ', stat(i)
        endif
        ! Test for phi out of order.
        if ( i > 1 .and. stat(i) >= good ) then
          ! print *, '- ', p_path(i), p_path(i-1), stat(i-1)
          if ( p_path(i) < p_path(i-1) .and. stat(i-1) >= good ) then
            stat(i) = order ! Phi out of order, can't be right
          end if
        end if
        if ( i < n_path .and. stat(i) >= good ) then
          ! print *, '+ ', p_path(i), p_path(i+1), stat(i+1)
          if ( any( stat(i+1) == (/ No_sol, NaN_sol /) ) )  then
            ! NAG now crashes under Oracle Linux 8 when this is the case
            stat(i+1) = NaN_sol
          elseif ( p_path(i) > p_path(i+1) .and. stat(i+1) >= good ) then
            stat(i) = order ! Phi out of order, can't be right
          end if
        end if
        if ( stat(i) >= good ) then ! still good (good < grid1 < grid2)
          no_bad_fits = no_bad_fits - 1
          if ( stat(i) >= grid1 ) no_grid_fits = no_grid_fits + 1
        end if
        if ( outside ) then
           ! phi is monotone; the remaining solutions are in the next
           ! column or even farther over, and not before the current row
           ! or the oldest "bad" row in this column.
          i1 = min(i,ibad)   
          ibad = i2
          exit
        end if
        if ( stat(i) >= good ) cycle  ! Newton iteration ended successfully
!        Don't do this; it causes trouble:
!        if ( stat(i) == grid1 ) cycle ! Can't be further to the right
        ibad=min(i,ibad)
        if ( (debug .or. complexDebug) .and. stat(i) == complex ) then ! Complex solution
            if ( debug ) &
              & print '(i4,i2,18x,f11.3,f12.3,1x,a)', &
                & i, j, h_ref(k,j)+r_eq, h_ref(k,j+1)+r_eq, 'Complex solution'
            if ( complexDebug ) then
              call output ( req_s, before='&in req_s = ' )
              call output ( h_surf, before=', H_Surf = ' )
              call output ( tan_ht_s, before=', tan_ht_s = ', after=',', advance='yes' )
              call output ( phi_t, before='  Phi_T = ' )
              call output ( rad2deg*phi_t, before=" = " )
              call output ( tan_ind, before=', tan_ind = ', advance='yes' )
              call dump ( h_ref(k,:), name='h_ref', lbound=k )
              call dump ( rad2deg*p_basis, name='p_basis (degrees)' )
            end if
        end if
      end do path ! i
    end do phi ! j

      if ( debug ) then
        print *, 'no_bad_fits =', no_bad_fits, ', no_grid+fits =', no_grid_fits
        call dump ( nStat(stat), name='Stat' )
      end if
    if ( no_bad_fits > 0 ) then
      call MLSMessage ( MLSMsg_Warning, ModuleName, &
        & "Height_Metrics failed to find H/Phi solution for some path segment" )
      if ( bin42 ) then
        rewind 42
        write ( 42 ) shape(p_path), shape(h_ref), shape(p_basis), phi_t, &
                   & tan_ind, h_surf, tan_ht_s, r_eq
        write ( 42 ) p_basis, h_ref, z_ref
        rewind 42
      end if
      if ( switchDetail(switches,'mhpx') > -1 ) then
        call dumpInput ( 1 )
        if ( .not. debug ) then
          call dump ( nStat(stat), name='Stat' )
        end if
        call dump ( h_path, 'H_Path before 1d fixup' )
        call dump ( rad2deg*p_path, 'P_Path (degrees) before 1d fixup' )
        call dump ( z_ref(tan_ind:), 'Z_Ref', lbound=tan_ind )
      end if
      ! We shouldn't get here at all, so don't worry about efficiency
      if ( stat(1) < good .or. stat(n_path) < good ) then
        call MLSMessage ( MLSMSG_Warning, moduleName, 'Resorting to 1d' )
        ! Get interpolating coefficients (eta_t) from p_basis to phi_t
        call eta_t%eta ( p_basis, phi_t )
        do i1 = 1, n_path
          if ( stat(i1) < good ) then
            ! Interpolate to tangent phi as a first guess
            h_path(i1) = eta_t%row_dot_vec ( 1, h_ref(vert_inds(i1),:) ) + r_eq
          end if
        end do
        ! Make sure H is monotone increasing away from the tangent point
        i1 = n_tan - 1
        i2 = 1
        do k = -1, 1, 2
          do j = i1, i2, k
            if ( h_path(j) <= h_path(j-k) ) then
              do i = j+k, i2, k
                if ( h_path(j-k) < h_path(i) ) exit
              end do
              if ( i == i2+k  ) then ! ran off the end
                h_path(j) = h_path(j-k) + 1.0 ! use SWAG method
              else
                h_path(j) = h_path(j-k) + (h_path(i)-h_path(j-k)) * &
                  & ( z_ref(vert_inds(j)) - z_ref(vert_inds(j-k))) / &
                  & ( z_ref(vert_inds(i)) - z_ref(vert_inds(j-k)))
              end if
            end if
          end do
          i1 = n_tan + ngp1 + 1
          i2 = n_path
        end do
        ! Just do them all to make sure they're consistent
        p_path = acos(tan_ht/h_path) * phi_sign + phi_offset
      else
        do j = n_tan-1, 2, -1
          if ( stat(j) == good ) cycle
          do i2 = j, 2, -1
            if ( stat(i2) == good ) exit
          end do
          do i1 = j, n_tan
            if ( stat(i1) == good ) exit
          end do
          h_path(j) = h_path(i1) + (h_path(i2)-h_path(i1)) * &
            & ( z_ref(vert_inds(j)) - z_ref(vert_inds(i1))) / &
            & ( z_ref(vert_inds(i2)) - z_ref(vert_inds(i1)))
          p_path(j) = acos(tan_ht/h_path(j)) * phi_sign(j) + phi_offset(j)
        end do
        do j = n_tan+ngp1+1, n_path-1
          if ( stat(j) == good ) cycle
          do i2 = j, n_path-1
            if ( stat(i2) == good ) exit
          end do
          do i1 = j, n_tan+ngp1, -1
            if ( stat(i1) == good ) exit
          end do
          h_path(j) = h_path(i1) + (h_path(i2)-h_path(i1)) * &
            & ( z_ref(vert_inds(j)) - z_ref(vert_inds(i1))) / &
            & ( z_ref(vert_inds(i2)) - z_ref(vert_inds(i1)))
          p_path(j) = acos(tan_ht/h_path(j)) * phi_sign(j) + phi_offset(j)
        end do
      end if
      if ( switchDetail(switches,'mhpx') > -1 ) then
        call dump ( h_path, name='H_Path after 1d fixup' )
        call dump ( rad2deg*p_path, name='P_Path (degrees) after 1d fixup' )
        if ( switchDetail(switches,'mhpx') > 0 ) &
          & call MLSMessage ( MLSMSG_Error, moduleName, 'Halt requested by mhpx1' )
      end if
    else if ( testStat ) then
      if ( any(stat < good) ) then
        call dump ( nStat(stat), name='Stat' )
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'No_Bad_Fits = 0 but some stat < good' )
      end if
    end if

    if ( tan_ht_s < 0.0 ) then ! Earth-intersecting ray
      if ( n_path >= 4 ) then
        ! Put the reflection point on the same side of the Earth as the rest
        ! of the path.
        if ( abs(p_path(n_tan-1) - p_path(n_tan)) > 0.9 * pi ) then
          ! Reflection point on the wrong side of the earth
          p_path(n_tan:n_tan+ngp1) = mod(p_path(n_tan:n_tan+ngp1) + &
            & sign(pi,p_path(n_tan-1) - p_path(n_tan)),pix2)
        end if
      end if
    end if

    ! Since we have solved for the intersection of sec(phi) with the heights
    ! on constant-zeta surfaces, it is impossible for the heights not to be
    ! monotone increasing away from the tangent point.

    if ( do_dumps >= 0 ) then
      call output ( tan_ht_s, before='tan_ht_s = ' )
      call output ( req_s, before=' km, r_eq+h_surf = ' )
      call output ( h_surf, before=' km, h_surf = ', after=' km', advance='yes' )
      call dump ( rad2deg*p_path, name='P_Path (degrees) before refractive correction', &
        & format='(f14.8)', options=options )
      if ( h_phi_dump < 0 .and. do_dumps >= 1 ) &
        & call dump ( rad2deg*p_basis, name='p_basis (degrees)', format='(f14.6)', options=options )
      call dump ( h_path, name='h_path (km from center of equivalent circular Earth)', &
        & format='(f14.6)', options=options )
      if ( do_dumps >= 2 ) call dump ( nStat(stat), name='Stat' )
      if ( do_dumps >= 2 ) &
        & call dump ( z_ref(tan_ind:), name='z_ref', options=options, lbound=tan_ind )
    end if

  contains

    subroutine DumpInput ( FirstRow )
      integer, intent(in) :: FirstRow
      call output ( tan_ht_s, before='Tan_Ht_s = ' )
      call output ( req_s, before=', R_Eq+H_Surf = ' )
      call dump ( rad2deg*p_basis, name='P_Basis (degrees)', format='(f14.6)', options=options )
      call output ( phi_t, before='Phi_t = ', format='(f11.8)' )
      call output ( rad2deg*phi_t, before=" = ", format='(f13.8)' )
      call output ( h_surf, before =', H_Surf = ' )
      call output ( tan_ind, before=', Tan_Ind = ', advance='yes' )
      call output ( r_eq, before='R_Eq ' )
      call dump ( z_ref(tan_ind:), name=' Z_Ref ', lbound=tan_ind )
      if ( tan_ht_s < 0.0 ) then
        call output ( n_tan, before='Phi_Offset (' )
        call output ( n_tan+ngp1, before=':', after=') (radians) ' )
        call dump ( phi_offset(n_tan:n_tan+ngp1), options='c' ) ! clean=.true.
        call output ( n_tan, before='Phi_Offset (' )
        call output ( n_tan+ngp1, before=':', after=') (degrees) ' )
        call dump ( rad2deg*phi_offset(n_tan:n_tan+ngp1), options='c' ) ! clean=.true.
      end if
      call dump ( h_ref(firstRow:n_vert,:), name='H_Ref', format='(f14.7)', &
        & options=options, lbound=firstRow )
      if ( debug ) call dump ( vert_inds, name='Vert_Inds' )
    end subroutine DumpInput

  end subroutine Height_Metrics

  ! ------------------------------------------------  Solve_H_Phi  -----

  !{Solve for an intersection of $H = H_t \sec \phi$ and the line from
  ! $(\phi_1,h_1)$ to $(\phi_2,h_2)$, where $\phi_1$ and $\phi_2$ are
  ! {\tt p\_basis(1)} and {\tt p\_basis(2)}, and $h_1$ and $h_2$ are
  ! {\tt h\_ref(1)} and {\tt h\_ref(2)}.
  !
  ! The solution is only acceptable if $\phi_1 \leq \phi \leq \phi_2$ and
  ! $(h_1 - h) (h_2-h) \leq 0$ or $|h_1 - h| \leq \tau$ or $|h_2 - h| \leq
  ! \tau$, where $\phi$ is {\tt p\_path}, $h$ is {\tt h\_path}, and $\tau$
  ! is the height tolerance {\tt tol}.

  subroutine Solve_H_Phi ( & ! inputs
    &                      p_basis, phi_offset, phi_sign, h_ref, a, b, c, &
    &                      tan_ht, r_eq, tol, inside, start, &
                             ! outputs and inouts
    &                      h_path, phi, stat, outside )

    use MLSKinds, only: RP
    use Output_m, only: Output

    real(rp), intent(in) :: P_Basis(2) ! Coordinates for H_Ref.
    real(rp), intent(in) :: Phi_Offset ! Offset from center of path.  Either the
                                       ! tangent point if the ray is not an
                                       ! Earth-intersecting ray, or the tangent
                                       ! point for the incident ray on one side,
                                       ! and the tangent point for the reflected
                                       ! ray on the other side.
    real(rp), intent(in) :: Phi_Sign   ! Which way from tangent?
    real(rp), intent(in) :: H_Ref(2)   ! Height reference, km above R_Eq.
      ! A, B and C are coefficients of a quadratic approximation to the solution
    real(rp), intent(in) :: A          ! (p_basis(2)-p_basis(1)) * tan_ht
    real(rp), intent(in) :: B          ! -(h_ref(2) - h_ref(1))
    real(rp), intent(in) :: C          ! -(h_ref(1)-h_surf) * (p_basis(2)-phi_t)
                                       ! +(h_ref(2)-h_surf) * (p_basis(1)-phi_t)
                                       ! +(p_basis(2)-p_basis(1)) * h_tan
    real(rp), intent(in) :: Tan_Ht     ! From equivalent circular Earth center,
                                       ! H_Tan + r_eq, km
    real(rp), intent(in) :: R_eq       ! From center of equivalent circular
                                       ! earth, km
    real(rp), intent(in) :: Tol        ! Height tolerance for Newton convergence,
                                       ! km
    logical, intent(in) :: Inside      ! P_Basis(2) is not the last one in the grid
    real(rp), intent(in) :: Start      ! Start here if > -0.5*huge(0.0_rp)
    real(rp), intent(inout) :: H_Path  ! H solution, inout in case there is none
    real(rp), intent(out) :: Phi       ! Phi solution
    integer, intent(inout) :: Stat     ! "good" or "grid" or "complex" or unchanged
    logical, intent(out) :: Outside    ! Newton iteration converged to a point
                                       ! outside p_basis(1:2)

    real(rp) :: D      ! B^2 - 4 A C
    real(rp) :: Deriv  ! Derivative of A (sec(p)-1 ) + B p + C or
                       ! its polynomial approximation if p < p8Tol
    real(rp) :: H      ! tan_ht * ( 1.0_rp + secM1 ) - r_eq,
                       ! ~ tan_ht * sec(p) - r_eq
    integer :: N       ! Newton iteration counter
    integer, parameter :: NMax = 10 ! Maximum number of Newton iterations
    real(rp) :: P, P2  ! Candidate solution, p^2
    real(rp) :: P_Path ! Candidate path phi solution
    real(rp) :: Secm1  ! Sec(p) - 1
    real(rp) :: Secp   ! Sec(p) = 1/Cos(p), to compute SecM1 when p > p8Tol

    ! Coefficients in expansion of Sec(phi)-1 = c2*phi^2 + c4*phi^4 ...
    ! At phi-phi_t = 0.2 radians and tan_ht = 6400 km, this gives 56 cm accuracy
    real(rp), parameter :: C2 = 0.5_rp, C4 = 5.0_rp/24, C6 = 61.0_rp/720
    real(rp), parameter :: C8 = 277.0_rp/8064
    ! Coefficients in expansion of Sec(phi)*Tan(phi) = d/dPhi(sec(phi)-1)
    real(rp), parameter :: D1 = 2*c2, D3 = 4*c4, D5 = 6*c6 ! ... 2n * c_2n
    real(rp), parameter :: D7 = 8*c8
    ! Boundary on P for using the expansion instead of sec(p) - 1
    real(rp), parameter :: P8Tol = 0.2

    ! For debugging output
    real(rp) :: DD(merge(10,0,debug))

    if ( stat == order ) stat = no_sol
    outside = .false.
    if ( start < -0.25*huge(0.0_rp) ) then
      ! Solve 1/2 a p^2 + b p + c to start
      d = b*b - 2.0 * a * c
      if ( d < 0.0 ) then ! Complex starting solution, assume 1-d pressure
        p = phi_sign*acos(tan_ht/max(tan_ht,0.5*(h_ref(1)+h_ref(2))+r_eq))
        if ( stat < good ) &
          & stat = complex ! In case Newton iteration doesn't converge
          if ( debug ) then
            call output ( p + phi_offset, format='(f11.5)', &
              & before='Complex starting point' )
            call output ( ( p + phi_offset) * rad2deg, format='(f11.5)', &
              & before=' radians = ', after=' degrees', advance='yes' )
          end if
      else
        p = (-b + phi_sign * sqrt(d) ) / a
        if ( p < p_basis(1)-phi_offset .or. p > p_basis(2)-phi_offset ) then
          p = (-b - phi_sign * sqrt(d) ) / a
          if ( p < p_basis(1)-phi_offset .or. p > p_basis(2)-phi_offset ) &
            & p = phi_sign*acos(tan_ht/max(tan_ht,0.5*(h_ref(1)+h_ref(2))+r_eq))
        end if
      end if
    else
      p = start - phi_offset
    end if
    p_path = sign(p,phi_sign) + phi_offset
    !{Use P in Newton iterations with an expansion for $\sec(\phi)$ that is
    ! higher order than two.  We don't use $\sec(\phi)-1$ for small $\phi$
    ! because it suffers cancellation.  Rather, use more terms of its Taylor
    ! series than we used for the quadratic approximation.
    do n = 1, nMax
      if ( abs(p) < P8Tol ) then ! Small phi -- sec(phi)-1 has cancellation
        p2 = p**2                             ! $ ( \phi - \phi_t )^2 $
        secM1 = p2*(c2+p2*(c4+p2*(c6+p2*c8))) ! ~ sec(p) - 1 to eighth order
        secp = secM1 + 1.0                    ! ~ sec(p) to eighth order
      else ! Large phi -- sec(phi)-1 is good, polynomial approximation is not
        secp = 1.0 / cos(p)                   ! sec(p)
        secM1 = secp - 1.0                    ! sec(p) - 1
      end if
      d = a*secM1 + b*p + c
        if ( debug ) dd(n) = d
      h_path = tan_ht * secp
      h = h_path - r_eq  ! to test against h_ref
      if ( abs(d) < tol ) then ! difference is small enough
        if ( (p_path-p_basis(1)) * (p_path-p_basis(2)) <= 0.0 .and. &
               ! Converged within bounds?
                &  ( (h-h_ref(1))*(h-h_ref(2)) <= 0.0 .or. &
               ! or the difference in reference heights is within tolerance and
               ! the current H is outside the bounds by less than tolerance
          &    abs(b) < tol .and. &
          &    ( abs(h-h_ref(1)) < tol .or. &
          &      abs(h-h_ref(2)) < tol ) ) ) then
          ! Not extrapolating in phi
          stat = good
          if ( debug ) call debug1 ( 1, n, 'GOOD', p_path, p_print=p )
          phi = p_path
          return
        end if
        if ( p_path > p_basis(2) .and. abs(a*(h-h_ref(2))) > tan_ht*abs(d) &
          & .and. inside ) then
          ! Phi is beyond the range, and H is outside by more than the
          ! move.  Go on to the next columns of H_ref and P_Basis.
          outside = .true. ! phi is monotone; the rest are in the next column
          phi = p_path ! so phi has a defined value, even if it's no good
          if ( debug ) then
            call debug1 ( 1, n, 'OUTSIDE', p_path, p_print=p )
          end if
          if ( debug ) call debug2
          return           ! or even farther over
        end if
      end if
      if ( abs(p) < P8Tol ) then
        deriv = a*p*(d1+p2*(d3+p2*(d5+p2*d7))) + b       ! a P8' + b
      else
        deriv = a * secp * sign(sqrt(secp**2-1.0),p) + b ! a sec(p) tan(p) + b
      end if
        if ( newtonDetails ) then
          call debug1 ( 1, 0, 'STEP ' // merge('P','S',abs(p) < P8Tol), p_path, p, deriv )
        end if
      p = p - d / deriv ! do the Newton step
      p_path = sign(p,phi_sign) + phi_offset
    end do ! n
    ! Newton iteration failed to converge.  Use the break point if the
    ! result is very near to it.  We'll probably find it in the next
    ! panel, but just in case....
      if ( debug ) then
        if ( .not. outside ) then
          call debug1 ( max(n-3,1), n-1, merge('CONVERGE', 'BOUNDS  ', abs(d) < tol), p_path, &
            & p_print=p )
        else
          call debug1 ( max(n-3,1), n-1, 'NONE', p_path, p_print=p )
        end if
      end if
    phi = p_path ! make sure it has a defined value, even if it's no good
    if ( stat < grid1 ) then
      if ( abs(h-h_ref(1)) <  0.1 .and. &
        & abs(p_path-p_basis(1)) < abs(p_path-p_basis(2)) ) then
        h_path = h_ref(1) + r_eq
        if ( h_path > tan_ht ) then
          phi = p_basis(1)
          stat = grid1
          if ( debug ) call debug1 ( 1, 0, 'GRID 1', phi )
        end if
      else if ( abs(h-h_ref(2)) <  0.1 .and. &
        & abs(p_path-p_basis(2)) < abs(p_path-p_basis(1)) ) then
        h_path = h_ref(2) + r_eq
        if ( h_path > tan_ht ) then
          phi = p_basis(2)
          stat = grid2
          if ( debug ) call debug1 ( 1, 0, 'GRID 2', phi )
        end if
      end if
    end if
    if ( p_path > p_basis(2) .and. inside ) then
      outside = .true. ! phi is monotone; the rest are in the next column
    end if
    if ( debug ) call debug2

  contains

    subroutine debug1 ( K, L, Why, Phi, P_print, Deriv )
      integer, intent(in) :: K, L
      character(len=*), intent(in) :: Why
      real(rp), intent(in) :: Phi
      real(rp), intent(in), optional :: P_print
      real(rp), intent(in), optional :: Deriv
      character(merge(13,0,debug)) :: OOPS
      if ( debug ) then
        if ( present(P_print) ) then
          write ( *, '(f11.5)', advance='no' ) ang*P_print
        else
          write ( *, '(11x)', advance='no' )
        end if
        oops=''
        if ( abs(cos(p)*(h+r_eq)-tan_ht) > 5.0e-4 * cos(p) ) oops='TRIG'
        if ( why /= '' ) oops(6:) = why
      100 format (f11.5,f10.3,f10.3,f11.3,f10.3,1p,g14.6,i3,11g10.2)
        write (*, 100, advance='no') ang*phi, h_path, &
          & h_ref(1)+r_eq, h_ref(2)+r_eq, &
        ! & tan_ht / cos(p), &
          & tan_ht * ( 1.0_rp + p2*((c2+p2*(c4+p2*c6))) ), & ! ~ tan_ht * sec(p)
          & d, n, dd(k:l)
        if ( present(deriv) ) then
          write ( *, '(f12.3)', advance='no' ) deriv
        else if ( k > l ) then
          write ( *, '(10x)', advance='no' )
        end if
        write (*, '(2x,a)') trim(adjustl(oops))
      end if
    end subroutine

    subroutine debug2
      call output ( phi, before='Phi ', format='(f11.5)' )
      call output ( phi*rad2deg, before=' = ', format='(f11.5)' )
      call output ( stat, before=' stat ', after=' = ' )
      call output ( nStat(stat), advance='yes' )
    end subroutine debug2

  end subroutine Solve_H_Phi

  ! -----------------------------------------------  More_Metrics  -----
  subroutine More_Metrics ( &
          ! Inputs:
          & Tan_Ind, N_Tan, T_Sv, Vert_Inds, T_Ref, dHidZij, P_Path,   &
          & Eta_P,                                                     &
          ! Outputs:
          & T_Path, dHitdZi,                                           &
          ! Optional inputs:
          & ddHidHidTl0, dHidTlm, Z_Ref,                               &
          ! Optional outputs:
          & ddHtdHtdTl0, dHitdTlm, dHtdTl0, dHtdZt,                    &
          & Eta_ZP, T_Der_Path_Flags, Tan_Phi_t )

    ! This subroutine computes metrics-related things after H_Path and
    ! P_Path are computed by Height_Metrics, and then perhaps augmented
    ! with, for example, the minimum Zeta point.

    use Dump_0, only: Dump
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use MLSStringLists, only: SwitchDetail
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Sparse_m, only: Sparse_t
    use Toggles, only: Switches

    ! inputs:

    integer, intent(in) :: Tan_Ind      ! Tangent height index, 1 = center of
                                        ! longest path
    integer, intent(in) :: N_Tan        ! Tangent index in path, usually n_path/2
    type(grids_t), intent(in) :: T_Sv   ! Temperature state vector stuff
    integer, intent(in) :: Vert_Inds(:) ! First (vertical) subscripts for
                                        ![zt]_ref  at points on the path
    real(rp), intent(in) :: T_Ref(:,:)  ! Temperatures at Z_Ref X t_sv%phi_basis
    real(rp), intent(in) :: dHidZij(:,:)! Vertical derivative at Z_Ref X t_sv%phi_basis
    real(rp), intent(in) :: P_Path(:)   ! Phi's on the path
    type(sparse_eta_t), intent(inout) :: Eta_P ! Interpolating coefficients from
                                        ! Temperature's Phi basis to P_Path.
                                        ! InOut so as not to default initialize
                                        ! everything.

    ! outputs:

    real(rp), intent(out) :: T_Path(:)  ! computed temperatures on the path                            
    real(rp), intent(out) :: dHitdZi(:) ! derivative of height wrt zeta
                                        !--may be useful in future computations

    ! optional inputs

    real(rp), optional, intent(in) :: ddHidHidTl0(:,:,:) ! second order
             ! reference temperature derivatives. This is (height, phi_basis,
             ! zeta_basis). Needed only if present(dHidTlm).
    real(rp), optional, intent(inout) :: dHidTlm(:,:,:) ! reference temperature
             ! derivatives. This gets adjusted so that at ref_h(1,@tan phi)) is
             ! 0.0 for all temperature coefficients.
             ! This is height X zeta_basis X phi_basis
    real(rp), optional, intent(in) :: Z_Ref(:)   ! -log pressures (zetas) for
             ! which derivatives are needed.  Only the parts from the tangent
             ! outward are used.  Needed only if present(dHidTlm).

    ! Optional outputs.

    real(rp), optional, intent(out), target :: ddHtdHtdTl0(:)  ! Second order
             ! derivatives of height w.r.t T_Ref at the tangent only -- used
             ! for antenna affects. Computed if present(dHidTlm).
    type(sparse_t), optional, intent(inout) :: dHitdTlm        ! Derivative of
             ! path position wrt temperature state vector
             ! (t_sv%zet_basis X t_sv%phi_basis)
    real(rp), optional, intent(out), target :: dHtdTl0(:)      ! First order derivatives
             ! of height w.r.t T_Ref at the tangent only.  Computed if
             ! present(dHidTlm).
    real(rp), optional, intent(out) :: dHtdZt          ! Height derivative wrt
             ! pressure at the tangent.  Computed if present(dHidTlm).
    type(sparse_eta_t), optional, intent(inout) :: Eta_ZP ! Interpolating
             ! coefficients for Temperature from Zeta X Phi to the path
    logical, optional, intent(out) :: T_Der_Path_Flags(:) ! An element in any
             ! column of row I of Eta_ZP is nonzero. Computed if
             ! present(dHidTlm).
    real(rp), optional, intent(out) :: Tan_Phi_t ! temperature at the tangent

    ! Local variables.

    integer :: Do_Dumps    ! <0 = no dump, >=0 = dump, >0 = stop
    integer :: I           ! Subscript, loop inductor
    integer :: N_Path      ! Path length = size(vert_inds)

!   It would be nice to do this the first time only, but the
!   retrieve command in the L2CF can now change switches
!   if ( do_dumps < 0 ) then ! First time only
      do_dumps = switchDetail(switches,'metd')
!   end if

    n_path = size(vert_inds)

    ! Interpolate Temperature (T_Ref) and the vertical height derivative
    ! (dHidZij) to the path (T_Path and dHitdZi).

    ! Compute the interpolating coefficients Eta_p for temperature.
    call eta_p%empty ! Clean out data but don't deallocate components
    call eta_p%eta_1d ( t_sv%phi_basis, p_path(:n_path), &
                      & create=.false., sorted=.true. )
    do i = 1, n_path
      ! Interpolate to t_path and dHitdZi one row at a time because
      ! t_ref(vert_inds,:) would be a big copy.
      t_path(i) = eta_p%row_dot_vec ( i, t_ref(vert_inds(i),:) )
      dHitdZi(i) = eta_p%row_dot_vec ( i, dHidZij(vert_inds(i),:) )
    end do

    ! Now for the optional tangent quantities.
    if ( n_tan <= n_path ) then
      ! Why aren't these two just t_ref(n_tan) and dHitdZi(n_tan)?
      if ( present(tan_phi_t) ) &
        & tan_phi_t = eta_p%row_dot_vec ( i, t_ref(tan_ind,:) )
      if ( present(dhtdzt) ) &
        & dhtdzt = eta_p%row_dot_vec ( i,  dHidZij(tan_ind,:) )
      ! compute temperature derivatives
      if ( present(dHidTlm) ) call Temperature_Derivatives
    end if

    if ( do_dumps >= 1 ) then
      call dump ( t_ref, name='t_ref', format='(1pg14.6)', options=options )
      call dump ( t_path(:n_path), name='T_Path', format='(1pg14.6)', options=options )
      call dump ( dHitdZi(:n_path), name='dHitdZi', format='(1pg14.6)', options=options )
      call eta_p%dump ( name='Eta_P', format='(1pg14.6)' )
      if ( do_dumps > 1 ) stop
    end if

  contains

    subroutine Temperature_Derivatives

      use Array_Stuff, only: Element_Position, Subscripts
      use Comp_Eta_DoCalc_Sparse_m, only: Comp_One_Eta_Z
      use GLNP, only: NG, NGP1
      logical :: Change                      ! Interpolator got set to zero
      real(rp), pointer :: ddHtdHtdTl0_2(:,:) ! Zeta X Path
      real(rp), pointer :: dHtdTl0_2(:,:)    ! Zeta X Path
      type(sparse_eta_t) :: Eta_Z            ! Zeta interpolating coefficients
      integer :: I, J, K, SV_P, SV_T, SV_Z   ! Loop inductors and subscripts
      integer :: P_Coeffs                    ! # of phi's, t_sv%l_p(1)
      integer :: Two_D_Bounds(2)             ! [ Z_Coeffs, P_Coeffs ]
      integer :: Two_Subs(2)                 ! [ Sv_Z, Sv_P ]
      real(rp) :: V                          ! Value to store in dHitdTlm
      integer :: WS                          ! t_sv%window_start(1)
      integer :: Z_Coeffs                    ! # of zetas, t_sv%l_z(1)
      real(rp) :: Z_Path(n_path)             ! Z_Ref(vert_inds(1:n_path))
      equivalence ( two_d_bounds(1), z_coeffs ), ( two_d_bounds(2), p_coeffs )
      equivalence ( two_subs(1), sv_z ), ( two_subs(2), sv_p )

      ! Adjust the 2d hydrostatic temperature derivative relative to the
      ! surface. Even though this is updated on every invocation, that is,
      ! with a new phi_t, it works as if the original value were updated with
      ! the current phi_t, because the interpolation represented by eta_p is
      ! linear. Thus, the effect of cumulative updates for each new phi_t are
      ! the same as starting from the original dHidTlm and updating with the
      ! latest phi_t.  The algebra is horrible, but Maple has verified this.

      p_coeffs = t_sv%l_p(1) ! Also sets two_d_bounds(2)
      z_coeffs = t_sv%l_z(1) ! Also sets two_d_bounds(1)

      do i = 1, z_coeffs
        dHidTlm(:,i,:) = dHidTlm(:,i,:) - eta_p%row_dot_vec ( n_tan, dHidTlm(1,i,:) )
      end do

      dHtdTl0_2(1:z_coeffs,1:p_coeffs) => dHtdTl0 ( 1 : z_coeffs * p_coeffs )
      ddHtdHtdTl0_2(1:z_coeffs,1:p_coeffs) => ddHtdHtdTl0 ( 1 : z_coeffs * p_coeffs )

      dHtdTl0_2 = 0
      ddHtdHtdTl0_2 = 0
      ! Now fill the places that have nonzero Phi interpolating coefficients
      ! in the tangent row.
      do i = 1, z_coeffs
        call eta_p%row_times_vec ( n_tan, dHidTlm(tan_ind,i,:), dHtdTl0_2(i,:) )
        call eta_p%row_times_vec ( n_tan, ddHidHidTl0(tan_ind,i,:), &
                                 &        ddHtdHtdTl0_2(i,:) )
      end do

      ! Compute Zeta interpolation coefficients from the temperature zeta
      ! basis to path zetas.

      z_path = z_ref(vert_inds(1:n_path))
      call comp_one_eta_z ( t_sv, n_tan, z_path, eta_z, skip=ngp1 )

      ! Compute interpolation coefficients from Zeta X Phi to the path.
      call eta_zp%empty ! clean it out but don't deallocate components
      ! Compute Eta_ZP =  Eta_Z x Eta_P^T where t_sv%deriv_flags is true
      call eta_zp%eta_nD ( eta_z, eta_p, flags=t_sv%deriv_flags )
      ! Don't need Eta_ZP at GL points between two coarse tangent points
      if ( n_path > n_tan ) eta_zp%rows(n_tan+1:n_tan+ng) = 0
      ! Compute the path temperature derivative, noting where the nonzeros are
      call dHitdTlm%empty ! Clean it out but don't deallocate components
      dHitdTlm%nRows = n_path ! usually filled by Sparse_Eta%*D
      t_der_path_flags = .false.
      do i = 1, n_path
        t_der_path_flags(i) = eta_zp%rows(i) /= 0
        ! An iterator to traverse a row of a sparse_t would be helpful
        j = eta_p%rows(i) ! Last element in the row
        if ( j /= 0 ) then
          do
            j = eta_p%e(j)%nr   ! Next element in the row
            sv_p = eta_p%e(j)%c ! Column subscript of the element
            do sv_z = 1, z_coeffs
              sv_t = element_position ( two_subs, two_d_bounds )
            ! sv_t = element_position ( [ sv_z, sv_p ], [ z_coeffs, p_coeffs ] )
              if ( t_sv%deriv_flags(sv_t) ) then
!               v = max(dHidTlm(vert_inds(i),sv_z,sv_p),0.0_rp) * eta_p%e(j)%v
                v = dHidTlm(vert_inds(i),sv_z,sv_p) * eta_p%e(j)%v
                if ( v > 0 ) call dHitdTlm%add_element ( v, i, sv_t )
              end if
            end do
            if ( j == eta_p%rows(i) ) exit ! just processed the last element
          end do
        end if
      end do

    end subroutine Temperature_Derivatives

  end subroutine More_Metrics

  ! ------------------------------------------------  More_Points  -----
                         ! Inputs:
  subroutine More_Points (  phi_t, tan_ind, p_basis, z_ref, h_ref, r_eq,  &
                         &  h_surf, h_tan, p_path,                        &
                         ! Outputs:
                         &  z_new, h_new, p_new, n_new,                   &
                         ! Optional inputs:
                         &  h_tol )

    ! Check for path crossings in h_ref from tan_ind down, which Height_Metrics
    ! doesn't do.

    use MLSKINDS, only: RP
    ! inputs:

    real(rp), intent(in) :: phi_t      ! Orbit projected tangent geodetic angle
    integer, intent(in) :: tan_ind     ! Tangent height index, 1 = center of
                                       !  longest path
    real(rp), intent(in) :: p_basis(:) ! Horizontal temperature representation
                                       !  basis
    real(rp), intent(in) :: z_ref(:)   ! Reference zetas
    real(rp), intent(in) :: h_ref(:,:) ! Heights by z_ref and p_basis
    real(rp), intent(in) :: R_eq       ! equivalent elliptical earth radius at
                                       ! H_Surf
    real(rp), intent(in) :: H_Surf     ! Height of the pressure reference
                                       ! surface
    real(rp), intent(in) :: H_Tan      ! Tangent height above H_Surf -- negative
                                       !  for Earth-intersecting ray
    real(rp), intent(in) :: p_path(:)  ! From Height_Metrics

    ! outputs:
    real(rp), intent(out) :: z_new(:)  ! computed addional zetas
    real(rp), intent(out) :: h_new(:)  ! computed heights, referenced to Earth center
    real(rp), intent(out) :: p_new(:)  ! computed phis
    integer, intent(out) :: n_new      ! How many points put into h_new, p_new

    ! optional inputs
    real(rp), optional, intent(in) :: H_Tol ! Height tolerance in kilometers
                                       !   for convergence of phi/h iteration

    ! Local variables
    real(rp) :: A, B, C
    real(rp) :: H1, H2    ! Height along line of sight at p_basis(j-1:j)
    integer :: I, J
    real(rp) :: My_H_Tol
    logical :: None
    logical :: Outside
    integer :: P_Coeffs
    real(rp) :: REQ_S     ! R_eq - H_Surf = radius at true surface
    integer :: Stat
    real(rp) :: Tan_Ht    ! H_Tan + R_eq = height above h_surf

    if ( h_tan < 0.0 ) return ! Can't be any intersections below the surface
    my_h_tol = defaultTol ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    req_s = r_eq - h_surf
    tan_ht = H_Tan + R_eq
    p_coeffs = size(p_basis)
    n_new = 0
    h2 = 0.0 ! Just so it's defined; this value is never used
    do i = min(tan_ind,ubound(z_ref,1)), 1, -1
      none = .true. ! Assume there will be no intersections
      do j = 1, size(h_ref,2)
        h1 = h2
        h2 = tan_ht / cos(p_basis(j) - phi_t) - req_s
        if ( j == 1 ) cycle ! It takes two to tango
        if ( (h1-h_ref(i,j-1)) * (h2-h_ref(i,j)) < 0.0 ) then
          ! Line of sight intersects constant-zeta surface.  Solve for where.
          n_new = n_new + 1
          a = (p_basis(j)-p_basis(j-1)) * tan_ht
          b = -(h_ref(i,j)-h_ref(i,j-1)) ! h_surf cancels here
          c = -(h_ref(i,j-1)-h_surf)*(p_basis(j  )-phi_t) &
            & +(h_ref(i,j  )-h_surf)*(p_basis(j-1)-phi_t) &
            & +(p_basis(j)-p_basis(j-1)) * h_tan
          if ( debug ) &
            & print 120, j, i, h_ref(i,j)+req_s, h_ref(i,j+1)+req_s, &
            & a, b, c
          120 format ( 4x, i2, i4, 14x, f11.3,f12.3, f9.3, 1p,g14.6, 3x, g14.6 )
          stat = no_sol
          call Solve_H_Phi ( p_basis(j-1:j), phi_t, sign(1.0_rp,phi_t-p_basis(j-1)), &
          &                h_ref(i,j-1:j), a, b, c, &
          &                tan_ht, req_s, my_h_tol, i /= tan_ind .or. j /= p_coeffs, &
          &                -0.5*huge(0.0_rp), & ! Start
          &                h_new(n_new), p_new(n_new), stat, outside )
          if ( stat >= grid1 ) then
            if ( minval(abs(h_new(n_new)- h_ref(i,j-1:j))) > my_h_tol ) then
              stat = no_sol
              n_new = n_new - 1
              cycle
            end if
          end if
          if ( stat < good .or. outside ) then
            ! What's going on?
            n_new = n_new - 1
          else
            if ( minval(abs(p_new(n_new)-p_path)) < 1.0e-4 ) then
              ! New one is too close to a grid point -- ignore it.
              stat = no_sol
              n_new = n_new - 1
            else
              none = .false. ! There is an intersection in this row of H_Ref
              z_new(n_new) = z_ref(i)
            end if
          end if
        end if
        if ( n_new == size(z_new)-1 ) return ! No room for any more
      end do ! j
      if ( none ) exit ! Heights are monotone in each column of H_Ref, so if
                       ! there aren't any intersections in this row, there
                       ! aren't any others farther down.
    end do ! i

  end subroutine More_Points

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Metrics_m

! $Log$
! Revision 2.91  2018/08/28 22:15:35  vsnyder
! Rearrange arguments because Eta_FZP is optional in Comp_Sps_Path_Sparse_No_Frq
!
! Revision 2.90  2018/08/15 01:15:47  vsnyder
! Add calculation of T_Der_Path_Flags for temperature derivatives
!
! Revision 2.89  2018/05/24 03:24:09  vsnyder
! Use sparse representation for dHitdTlm
!
! Revision 2.88  2018/05/17 02:15:45  vsnyder
! Use sparse instead of dense interpolation
!
! Revision 2.87  2018/05/14 23:37:35  vsnyder
! Change to sparse eta representation
!
! Revision 2.86  2017/09/20 01:16:22  vsnyder
! Revise some dumping, delete a redundant dump
!
! Revision 2.85  2017/09/14 19:42:24  vsnyder
! Better control over what's dumped using metd#
!
! Revision 2.84  2017/09/13 19:39:24  vsnyder
! Move nStat to module scope so more routines can use it.  Add unformatted
! output of inputs if debugging.  Add more debugging output.  Add more
! checking on phi order.  Add method to recover from complex starting point.
! Correct missing application of Phi_Sign, which caused phi to be out of
! order.  Some cannonball polishing.
!
! Revision 2.83  2017/09/08 16:44:33  pwagner
! Fixes bug that broke Platinum brick
!
! Revision 2.82  2017/08/28 20:28:08  livesey
! Changed the n,nf,np,nz elements to j,jf,...
!
! Revision 2.81  2017/07/26 20:05:22  vsnyder
! When the L4 component of Grids_t is constructed, the low bound for the
! first subscript is WindowStart, not 1.  That needs to be considered when
! looking at L4.
!
! Revision 2.80  2017/03/11 00:56:05  vsnyder
! Use list interpolation in More_Metrics, cosmetic changes
!
! Revision 2.79  2016/11/17 02:06:12  vsnyder
! Correct some LaTeX
!
! Revision 2.78  2016/11/17 01:39:45  vsnyder
! Change path variable names from xyz_grid to xyz_path
!
! Revision 2.77  2016/10/18 00:25:37  vsnyder
! Improve comments, some cannobball polishing
!
! Revision 2.76  2016/08/20 00:55:08  vsnyder
! Correct comments about units for Phi_t and P_Basis in Tangent_Metrics
!
! Revision 2.75  2015/09/22 01:59:32  vsnyder
! Correct some comments, dump Z
!
! Revision 2.74  2015/05/28 23:14:42  vsnyder
! Add units in comments about variables
!
! Revision 2.73  2013/06/12 02:32:30  vsnyder
! Cruft removal
!
! Revision 2.72  2013/05/21 23:55:14  vsnyder
! Compute No_Bad_Fits correctly, revise dump switches
!
! Revision 2.71  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.70  2013/03/30 00:06:11  vsnyder
! Put dummy argument declarations in same order as in subroutine statement
!
! Revision 2.69  2013/02/28 21:08:29  vsnyder
! Try not to access array out of bounds
!
! Revision 2.68  2011/05/09 18:02:19  pwagner
! Converted to using switchDetail
!
! Revision 2.67  2010/02/04 23:10:20  vsnyder
! Remove USE for unreferenced names
!
! Revision 2.66  2010/02/02 01:34:55  vsnyder
! Make do_calc_t intent(inout)
!
! Revision 2.65  2009/12/15 03:20:06  vsnyder
! Don't force phi_offset to -pix2..pix2 for Earth-intersecting rasy
!
! Revision 2.64  2009/12/09 21:32:03  vsnyder
! Remove ill-advised mod on p_grid
!
! Revision 2.63  2009/08/20 23:30:56  vsnyder
! Handwaving to put p_grid in -2 pi .. 2 pi and the reflection point of
! an Earth-intersecting path on the same side of the Earth as the rest of
! the path.
!
! Revision 2.62  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.61  2009/06/16 17:37:26  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.60  2009/06/13 01:14:29  vsnyder
! Extensive changes almost everywhere.  Correct many problems with
! Earth-reflecting rays.
!
! Revision 2.59  2008/06/26 00:28:16  vsnyder
! Cannonball polishing
!
! Revision 2.58  2008/05/20 00:16:40  vsnyder
! Correct some TeXnicalities
!
! Revision 2.57  2007/11/08 01:49:39  vsnyder
! Req should be Req_s in one place in More_Points
!
! Revision 2.56  2007/10/19 23:28:30  vsnyder
! Don't accept grid solutions less than Ht
!
! Revision 2.55  2007/10/11 20:17:35  vsnyder
! Accept grid solutions more readily
!
! Revision 2.54  2007/10/03 20:54:10  vsnyder
! Require phi to be monotone
!
! Revision 2.53  2007/09/07 01:37:23  vsnyder
! Spiff up some dumps
!
! Revision 2.52  2007/08/23 20:41:07  vsnyder
! Don't clobber "grid" status with "cplx" status on next panel
! Don't clobber "grid" phi with failed iteration value on next panel
! Separate "grid" into "grid1" and "grid2"; don't bother to look for
! intersection on next panel after "grid1" solution.
!
! Revision 2.51  2007/07/31 23:48:57  vsnyder
! Try to recover from failure of H/Phi convergence
!
! Revision 2.50  2007/06/29 19:34:40  vsnyder
! Make the default h/phi convergence tolerance a parameter
!
! Revision 2.49  2007/06/26 00:35:39  vsnyder
! Use column-sparse eta, make *_zxp arguments
!
! Revision 2.48  2007/06/08 22:05:48  vsnyder
! More work on min zeta
!
! Revision 2.47  2007/03/22 23:24:34  vsnyder
! Remove irrelevant USE
!
! Revision 2.46  2007/02/06 23:50:59  vsnyder
! Polish up some messages.  Use 8th order approximation for sec(phi)-1 since
! sixth order doesn't give one-meter accuracy.
!
! Revision 2.45  2007/02/06 21:09:29  vsnyder
! Don't skip unsolved segments when going to next phi
!
! Revision 2.44  2007/02/01 02:53:15  vsnyder
! Make Solve_H_Phi a subroutine, add More_Points subroutine
!
! Revision 2.43  2007/01/20 01:07:49  vsnyder
! Use Earth-intersection phi instead of tangent phi to compute tangent
! temerature derivatives, do some decrufting too.
!
! Revision 2.42  2007/01/19 02:38:53  vsnyder
! Include water in phi refractive correction
!
! Revision 2.41  2007/01/18 00:26:32  vsnyder
! Split Pure_Metrics into Tangent_Metrics and Height_Metrics
!
! Revision 2.40  2006/12/21 22:59:17  vsnyder
! Don't reference optional arguments that aren't present
!
! Revision 2.39  2006/12/20 21:22:16  vsnyder
! Split metrics into pure H-Phi calculation, and everything else, in
! preparation for inserting the minimum-Zeta point into the path.
!
! Revision 2.37  2006/12/13 02:31:35  vsnyder
! Revision 2.38  2006/12/19 02:50:35  vsnyder
! Get rid of STATUS, reference H_Grid to Earth center instead of surface,
! some cannonball polishing.
!
! Polish up a comment
!
! Revision 2.36  2006/12/09 02:25:42  vsnyder
! Size some arrays more accurately, use First:Last version of get_eta_matrix,
! remove Tangent_Temperature_Derivatives module procedure since it isn't used
!
! Revision 2.35  2006/12/08 23:57:08  vsnyder
! Revise earth-intersecting metrics
!
! Revision 2.34  2006/12/08 22:50:55  vsnyder
! Complete revision
!
