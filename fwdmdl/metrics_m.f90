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

  implicit NONE
  private
  public :: Height_Metrics, More_Metrics, Tangent_Metrics

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! --------------------------------------------  Tangent_Metrics  -----
                               ! Inputs
  subroutine Tangent_Metrics ( Phi_T, P_Basis, Z_Ref, H_Ref, Csq, &
                               ! Inout
    &                          Tan_Ind_C, NZ, &
                               ! Outputs
    &                          Req, H_Surf, H_Tan, Z_Ref_New, &
                               ! Optional inputs
    &                          Tan_press, Surf_temp, Surf_height )

  ! Compute the surface height, the tangent height, and the equivalent
  ! spherical earth radius.
  ! If the ray is an Earth-intersecting ray (H_TAN < 0) determine whether to
  ! copy Z_Ref to Z_Ref_New and add a new Zeta to it, in which case NZ =
  ! size(Z_Ref) + 1, else NZ = size(Z_Ref).

    use Geometry, only: EarthRadA
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use GLNP, only: NG, NGP1
    use Make_Z_Grid_m, only: Default_Thresh
    use MLSKinds, only: IP, RP, R8
    use Output_m, only: OUTPUT
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: phi_t      ! orbit projected tangent geodetic angle
    real(rp), intent(in) :: p_basis(:) ! horizontal temperature representation
                                       ! basis
    real(rp), intent(in) :: z_ref(:)   ! Coarse grid of -log pressures (zetas).
    real(rp), intent(in) :: h_ref(:,:) ! heights by fine zeta and p_basis
    real(rp), intent(in) :: csq        ! (minor axis of orbit plane projected
                                       ! Earth ellipse)**2

    ! inout
    integer, intent(inout) :: tan_ind_c ! Index of tangent in Z_Ref.  For an
                                       ! Earth-intersecting ray, it might be
                                       ! changed to a value such that
                                       ! Z_Ref_New(tan_ind_c) = Z_Tan.
    integer, intent(inout) :: NZ       ! Effective size of Z_Ref_New, increased
                                       ! by 1 if a new Earth-intersecting zeta
                                       ! is inserted.

    ! outputs:
    real(rp), intent(out) :: req       ! equivalent elliptical earth radius
    real(rp), intent(out) :: H_SURF    ! Height of the reference surface --
                                       ! interpolated in Surf_Height if
                                       ! present(Surf_Height) else interpolated
                                       ! in row 1 of H_REF.
    real(rp), intent(out) :: H_Tan     ! Tangent height above H_Surf (negative
    !                                  ! for Earth-intersecting rays)
    real(rp), intent(out) :: z_ref_new(:) ! Coarse grid of -log pressures (zetas).
                                       ! If the ray intersects the Earth
                                       ! surface, and the zeta of the
                                       ! intersection is sufficiently different
                                       ! from elements of Z_Ref, Z_Ref is copied
                                       ! here and a new one is added.

    ! optional inputs
    ! We have Tan_Press and Surf_Temp if the pointing is below the surface
    ! index (where zeta < -3) and we don't have Surf_Height.  Otherwise
    ! we have Surf_Height.
    real(rp), optional, intent(in) :: Tan_press    ! Tangent pressure
    real(rp), optional, intent(in) :: Surf_temp(:) ! Surface temperature at phi_basis
    real(rp), optional, intent(in) :: Surf_height(:) ! Surface height in km
                                       ! above mean sea level (whatever that
                                       ! means) on P_Basis.

    real(r8), parameter :: EarthRadA_2 = EarthRadA**2, EarthRadA_4 = EarthRadA**4

    integer :: Do_Dumps
    integer :: First, Last ! Nonzeros in Eta_T
    integer :: I
    real(rp) :: CP2, SP2   ! Cos^2 phi_t, Sin^2 phi_t
    real(rp) :: ETA_T(size(p_basis)) ! Interpolating coefficients
    real(rp) :: H_T1, H_T2 ! Used to determine where H_Tan is in H_Ref
    integer :: Tan_Ind_F   ! Tangent index in H_Ref, which is on fine Zeta grid
    real(rp) :: Z_Tan      ! Zeta at H_Tan

    do_dumps = max(index(switches,'metD'),index(switches,'metd'))

    tan_ind_f = (tan_ind_c-1) * ngp1 + 1

    ! Get interpolating coefficients (eta_t) from p_basis to phi_t
    call get_eta_sparse ( p_basis, phi_t, eta_t, first, last )

    if ( present(surf_height) ) then
      ! We set the surface reference at the actual surface height if we
      ! have it, and adjust the req and the h_tan relative to this, and
      ! adjust h_ref accordingly.
      h_surf = dot_product(surf_height(first:last), eta_t(first:last))
    else
      ! If we don't have the actual surface height, we set the surface
      ! reference at the input z_ref(1) and adjust the req and the h_tan
      ! relative to this, and adjust h_ref accordingly.
      h_surf = dot_product(h_ref(1,first:last),eta_t(first:last))
    end if

    !{ Compute equivalent earth radius (REQ) at phi\_t(1) (nearest surface)
    ! and adjust to the input z\_grid(1) using Equation (5.21) in the 19
    ! August 2004 ATBD JPL D-18130.
    !
    ! $c^2 = \frac{a^2 b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta}$.
    ! This is Equation (5.3) in the 19 August 2004 ATBD JPL D-18130.
    !
    ! $R^\oplus_{\text{eq}} \equiv H^\oplus_t =
    ! N(\phi_t) \sqrt{\sin^\phi_t + \frac{c^4}{a^4}\cos^2 \phi_t} + h_{\text{surf}} =
    ! \sqrt{\frac{a^4 \sin^2 \phi_t + c^4 \cos^2 \phi_t}
    !            {a^2 \cos^2 \phi_t + c^2 \sin^2 \phi_t}} + h_{\text{surf}}$
    !
    ! $a$ and $c$ are in meters; we want $R^\oplus_{\text{eq}}$ in kilometers.

    cp2 = cos(phi_t)**2
    sp2 = 1.0_rp - cp2
    req = 0.001_rp*sqrt((earthrada_4*sp2 + csq**2*cp2) / &
                      & (earthrada_2*cp2 + csq*sp2)) + h_surf

    ! compute the tangent height distance above H_surf.
    if ( present(tan_press) .and. .not. present(surf_height) ) then
      ! Earth intersecting ray. Compute GP height (km) of tangent pressure
      ! below surface. This will be negative because tan_press < z_ref(1).
      h_tan = dot_product(surf_temp(first:last),eta_t(first:last)) * &
        &     (tan_press-z_ref(1))/14.8
    else
      h_tan = dot_product(h_ref(tan_ind_f,first:last),eta_t(first:last)) - h_surf
    end if

    ! Determine whether the ray intersects the Earth surface
    if ( h_tan < 0.0 ) then ! Earth intersecting ray
      ! Find where h_tan is in h_ref
      h_t1 = h_tan
      do tan_ind_f = tan_ind_f, size(h_ref,1), ngp1
        h_t2 = dot_product(h_ref(tan_ind_f+ngp1,first:last),eta_t(first:last)) &
          &    - h_surf
        if ( h_t2 >= 0 ) exit
        h_t1 = h_t2
      end do
      tan_ind_c = tan_ind_f/ngp1 + 1
      ! h_t1 < 0 <= h_t2 here, so z_ref(tan_ind_c) <= z_tan <= z_ref(tan_ind_c+1)
      ! Interpolate linearly to get z_tan
      z_tan = z_ref(tan_ind_c) + h_t2/(h_t2-h_t1)
      if ( abs(z_tan-z_ref(tan_ind_c+1)) <= default_thresh ) then
        ! Tangent zeta is at tan_ind_c + 1; there is no new zeta
        tan_ind_c = tan_ind_c + 1
      else if ( abs(z_tan-z_ref(tan_ind_c)) > default_thresh ) then
        ! Add a new zeta
        z_ref_new(:tan_ind_c-1) = z_ref(:tan_ind_c-1)
        z_ref_new(tan_ind_c) = z_tan
        z_ref_new(tan_ind_c+1:nz+1) = z_ref(tan_ind_c:nz)
        nz = nz + 1
      else
        ! Tangent zeta is at tan_ind_c; nothing needs doing
      end if
      call output ( z_tan, before='Z_Tan = ', after=', ' )
    end if

    if ( do_dumps > 0 ) then
      call output ( h_tan, before='H_Tan = ' )
      call output ( h_surf, before=', H_Surf = ' )
      call output ( tan_ind_c, before=', Tan_ind_c = ' )
      call output ( tan_ind_f, before=', Tan_ind_f = ', advance='yes' )
    end if

  end subroutine Tangent_Metrics

  ! ---------------------------------------------  Height_Metrics  -----

                            ! Inputs:
  subroutine Height_Metrics (  phi_t, tan_ind, p_basis, z_ref, h_ref, req, &
                            &  h_surf, h_tan,                              &
                            ! Outputs:
                            &  h_grid, p_grid,                             &
                            ! Optional inputs:
                            &  h_tol )

    !{ This subroutine computes {\tt h\_grid} and {\tt p\_grid} that define
    !  a 2-d integration path.  These points are at the intersections of the
    !  line of sight given by $H = H_\text{tan} \sec \phi$ and the heights of
    !  surfaces of constant $\zeta$ represented in piecewise linear segments
    !  by consecutive elements of each row of {\tt H\_ref}.

    use Dump_0, only: Dump
    use Geometry, only: EarthRadA
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: IP, RP, R8
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_m, only: OUTPUT
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: phi_t      ! Orbit projected tangent geodetic angle
    integer, intent(in) :: tan_ind     ! Tangent height index, 1 = center of
    !                                     longest path
    real(rp), intent(in) :: p_basis(:) ! Horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_ref(:)   ! -log pressures (zetas) for which
    !                                     heights/temps are needed.  Only the
    !                                     parts from the tangent outward are used
    real(rp), intent(in) :: h_ref(:,:) ! Heights by z_ref and p_basis
    real(rp), intent(in) :: req        ! Equivalent elliptical earth radius
    real(rp), intent(in) :: H_Surf     ! Height at the Earth surface
    real(rp), intent(in) :: H_Tan      ! Tangent height above H_Surf -- negative
    !                                     for Earth-intersecting ray

    ! outputs:
    real(rp), intent(out) :: h_grid(:) ! computed heights, referenced to Earth center
    real(rp), intent(out) :: p_grid(:) ! computed phi's

    ! optional inputs
    real(rp), optional, intent(in) :: H_Tol ! Height tolerance in kilometers
                                       !   for convergence of phi/h iteration

    ! Local variables.
    integer :: Do_Dumps    ! 0 = no dump, >0 = dump
    integer :: Dump_Stop   ! 0 = no dump, >0 = dump/stop
    integer :: H_Phi_Dump  ! 0 = no dump, >0 = dump
    integer :: H_Phi_Stop  ! 0 = no dump, >0 = dump/stop
    integer :: First, Last ! Nonzeros in Eta_T
    integer :: I, I1, I2, J, K, N
    integer :: NO_BAD_FITS
    integer :: NO_GRID_FITS
    integer :: N_PATH      ! Path length = 2*(size(z_ref)+1-tan_ind)
    integer :: N_Tan       ! Tangent index in path, N_Path/2
    integer :: N_VERT      ! Size(z_ref)
    integer :: P_COEFFS    ! Size(P_basis)

    integer :: VERT_INDS(2*(size(z_ref)+1-tan_ind)) ! What to use in z_ref
    integer :: Stat(size(vert_inds))
    ! Values for stat:
    integer, parameter :: No_sol = 0 ! No solution
    integer, parameter :: Good = 1   ! Newton converged, tangent point, extrapolated
    integer, parameter :: Grid = 2   ! Close to a grid point

    real(rp) :: A, B, C, D ! Polynomial coefficients used to solve for phi and H
    real(rp) :: DP         ! p_basis(j+1) - p_basis(j)
    real(rp) :: DPJ0       ! p_basis(j  )-phi_t
    real(rp) :: DPJ1       ! p_basis(j+1)-phi_t
    real(rp) :: CP2        ! Cos^2 phi_t
    real(rp) :: H          ! Tentative solution for H
    real(rp) :: H_T2       ! Used when searching for surface height in h_ref
    real(rp) :: HTAN_R     ! H_Tan + req
    real(rp) :: My_H_Tol   ! Tolerance in kilometers for height convergence
    real(rp) :: P, Q       ! Tentative solutions for phi
    real(rp) :: P2         ! P**2
    real(rp) :: REQ_S      ! Req - H_Surf
    real(rp) :: SecM1      ! Taylor series for sec(phi)-1
    real(rp) :: SP2        ! Sin^2 phi_t
    real(rp) :: Z_Tan      ! Zeta at the tangent point, usually from
                           ! z_ref(tan_ind), but solved here for earth-
                           ! intersecting rays.

    real(rp) :: ETA_T(size(p_basis))
    real(rp) :: PHI_OFFSET(size(vert_inds)) ! PHI_T or a function of NEG_H_TAN
    real(rp) :: PHI_SIGN(size(vert_inds))   ! +/- 1.0

    ! Coefficients in expansion of Sec(phi)-1 = c2*phi^2 + c4*phi^4 ...
    real(rp), parameter :: C2 = 0.5_rp, C4 = 5.0_rp/24, C6 = 61.0_rp/720.0
    ! Coefficients in expansion of Sec(phi)*Tan(phi) = d/dPhi(sec(phi)-1)
    real(rp), parameter :: D1 = 2*c2, D3 = 4*c4, D5 = 6*c6 ! ... 2n * c_2n

    real(r8), parameter :: EarthRadA_2 = EarthRadA**2, EarthRadA_4 = EarthRadA**4

    ! To control debugging
    logical, parameter :: Debug = .false.
    logical, parameter :: NewtonDetails = .true. .and. debug
    ! For debugging output format:
    logical, parameter :: clean = .false.
    ! For debugging output
    real(rp) :: DD(merge(10,0,debug))
    character(merge(10,0,debug)) :: OOPS

!   It would be nice to do this the first time only, but the
!   retrieve command in the L2CF can now change switches
!   if ( do_dumps < 0 ) then ! First time only
      dump_stop = index(switches,'metD')
      do_dumps = max(dump_stop,index(switches,'metd'))
      h_phi_stop = index(switches,'hphI')
      h_phi_dump = max(h_phi_stop,index(switches,'hphi'))
!   end if

    n_path = size(vert_inds)
    n_tan = n_path / 2

    my_h_tol = 0.001_rp ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    p_coeffs = size(p_basis)
    n_vert = size(z_ref)

    ! Get interpolating coefficients (eta_t) from p_basis to phi_t
    call get_eta_sparse ( p_basis, phi_t, eta_t, first, last )

    ! p_basis and p_grid are phi's in offset radians relative to phi_t, that
    ! is, the phi_t, p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t

    phi_offset(:n_tan) = phi_t
    if ( h_tan < 0.0 ) then ! Earth-intersecting ray
      phi_offset(n_tan+1:) = phi_t-2.0_rp*Acos((req+h_tan)/req)
    else
      phi_offset(n_tan+1:) = phi_t
    end if

    req_s = req - h_surf
    htan_r = h_tan + req

    if ( h_phi_dump > 0 ) then
      call dump ( p_basis, name='p_basis', format='(1pg14.6)', clean=clean )
      call output ( phi_t, before='phi_t = ', format='(1pg14.6)' )
      call output ( h_surf, format='(f7.2)', before =', h_surf = ' )
      call output ( req, format='(f7.2)', before=', req = ', advance='yes' )
      call output ( tan_ind, before='tan_ind = ' )
      call dump ( phi_offset, name='phi_offset', clean=clean )
      call dump ( h_ref(tan_ind:n_vert,:), name='h_ref', clean=clean )
!     call dump ( z_ref(tan_ind:n_vert), name='z_ref', clean=clean )
    end if

      if ( debug ) &
        & print '(2a,f9.3,a,f7.3)', '   I J    phi  H(TAYLOR)',&
          & '  HREF(I,J) HREF(I,J+1)  H(TRIG)    H-HREF      N   Diffs'

    ! Only use the parts of the reference grids that are germane to the
    ! present path

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (i, i=tan_ind,n_vert) /)

    ! sign of phi vector
    phi_sign = (/ (-1, i=n_vert,tan_ind,-1), (+1, i=tan_ind,n_vert) /)

    stat = no_sol ! No solutions
    ! The tangent point is special
    stat(n_tan:n_tan+1) = good
    p_grid(n_tan:n_tan+1) = phi_offset(n_tan:n_tan+1)
    h_grid(n_tan:n_tan+1) = htan_r
    no_bad_fits = n_path - 2
    no_grid_fits = 0

    ! Get solutions outside of p_basis, if any, assuming constant H
    ! outside of p_basis.
    do i1 = 1, n_path
      if ( stat(i1) == good ) cycle ! Skip the tangent point
      k = vert_inds(i1)
      h = max(htan_r,h_ref(k,1)+req_s)
      p = acos(htan_r/h) * phi_sign(i1) + phi_offset(i1)
      if ( p >= p_basis(1) ) exit ! phi is monotone
      p_grid(i1) = p
      h_grid(i1) = h
      no_bad_fits = no_bad_fits - 1
      stat(i1) = good
        if ( debug ) then
          print 3, i1, 0, p_grid(i1), h_grid(i1), '  Outside p_basis'
        3 format ( i4,i2,f9.5,f9.3,a )
        end if
    end do
    j = p_coeffs
    do i2 = n_path, 1, -1
      if ( stat(i2) == good ) cycle ! Skip the tangent point
      k = vert_inds(i2)
      h = max(htan_r,h_ref(k,j)+req_s)
      p = acos(htan_r/h) * phi_sign(i2) + phi_offset(i2)
      if ( p <= p_basis(p_coeffs) ) exit ! phi is monotone
      p_grid(i2) = p
      h_grid(i2) = h
      no_bad_fits = no_bad_fits - 1
      stat(i2) = good
        if ( debug ) &
          & print 3, i2, p_coeffs, p_grid(i2), h_grid(i2), '  Outside p_basis'
    end do
      if ( debug ) then
        print 3, n_tan, 0, p_grid(n_tan), h_grid(n_tan), '  Tangent point'
        print 3, n_tan+1, 0, p_grid(n_tan+1), h_grid(n_tan+1), '  Tangent point'
      end if

    !{ Let $H_{ij}$ denote an element of {\tt H\_ref}.  Remember that all
    ! heights are referenced from $R^\oplus_{\text{eq}}$.  Calculate the
    ! heights by solving
    ! \begin{equation*}
    ! H_{ij} + (H_{i,j+1}-H_{ij})
    !                   \frac{\phi-\phi_t-(\phi_j-\phi_t)}
    !                        {\phi_{j+1}-\phi_t-(\phi_j-\phi_t)} - H_t
    ! = (H_t+R^\oplus_{\text{eq}})(\sec (\phi-\phi_t) - 1)\,.
    ! \end{equation*}
    ! 
    ! Notice that $R^\oplus_{\text{eq}}$ cancels everywhere on the left-hand
    ! side.  Using $\sec\phi \approx 1 + \frac12 \phi^2$, writing
    ! $\delta\phi=\phi-\phi_t$ and $H_{\text{tan}}=H_t+R^\oplus_{\text{eq}}$,
    ! and putting the resulting polynomial in standard form, we solve for
    ! $\delta\phi$ in
    ! 
    ! \begin{equation*}
    ! \frac12 (\phi_{j+1}-\phi_j) H_{\text{tan}} \delta\phi^2 -
    ! (H_{i,j+1}-H_{ij})\delta\phi -
    ! H_{ij} (\phi_{j+1}-\phi_t) + H_{i,j+1} (\phi_j-\phi_t) +
    ! (\phi_{j+1}-\phi_j)H_t
    ! \approx 0
    ! \end{equation*}
    !
    ! which arises from solving for the intersections of the line of sight
    ! ($H = H_{\text{tan}} \sec \phi)$ with constant-$\zeta$ surfaces,
    ! piecewise linear segments of which are represented by two consecutive
    ! columns in each row of {\tt H\_ref} (see wvs-048).  Extrapolated
    ! solutions are not accepted.  Outside of P\_Basis assume constant {\tt
    ! H\_ref}.  For I $<$ N\_Tan, try the minimum solution, otherwise try
    ! the maximum solution.

phi:do j = 1, p_coeffs-1
      dp = p_basis(j+1)-p_basis(j)
      a = dp * htan_r ! Actually 2*A
      dpj1 = p_basis(j+1)-phi_t
      dpj0 = p_basis(j)-phi_t
path: do i = i1, i2
        if ( stat(i) == good ) cycle
        k = vert_inds(i)
        b = -(h_ref(k,j+1) - h_ref(k,j)) ! h_surf cancels here
        c = -(h_ref(k,j)-h_surf)   * dpj1 &
          & +(h_ref(k,j+1)-h_surf) * dpj0 &
          & + dp * h_tan
        d = b*b - 2.0 * a * c ! Actually b^2-4ac, since a is actually 2*A
        if ( d < 0.0 ) then ! Complex solution
            if ( debug ) &
              & write (*, '(i4,i2,18x,f11.3,f12.3,1x,a)' ) &
                & i, j, h_ref(k,j)+req_s, h_ref(k,j+1)+req_s, 'Complex solution'
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            & "Complex starting value for Newton iteration shouldn't happen" )
          cycle  ! in case error level changed, and then crash somewhere else
        end if
        d = sqrt(d)
        p = (-b + phi_sign(i) * d ) / a
        !{Use P in Newton iterations with a higher-order expansion for
        ! $\sec(\phi)$.  We don't use $\sec(\phi)-1$ because it suffers
        ! cancellation for small $\phi$.  Rather, use more terms of its
        ! Taylor series than we used for the quadratic approximation.
        do n = 1, 10
          p2 = p**2
          secM1 = p2*((c2+p2*(c4+p2*c6))) ! ~ sec(p)-1
          d = a*secM1 + b*p + c
            if ( debug ) dd(n) = d
          h = htan_r * ( 1.0_rp + secM1 ) - req_s ! ~ htan_r * sec(p) - req + h_surf
          if ( abs(d) < my_h_tol ) then ! difference is small enough
            q = p + phi_offset(i)
              if ( newtonDetails ) then
                print 110, n, q, h+req_s, h_ref(k,j)+req_s, h_ref(k,j+1)+req_s,h_ref(k,j+1)-h_ref(k,j)
            110 format ( 4x,i2,f9.5,f9.3,f11.3,f12.3,11x,g14.6 )
              end if
            if ( q >= p_basis(j) .and. q <= p_basis(j+1) .and. &
                   ! Converged within bounds?
                    &  ( (h-h_ref(k,j))*(h-h_ref(k,j+1)) <= 0.0 .or. &
                   ! or the difference in reference heights is within tolerance and
                   ! the current H is outside the bounds by less than tolerance
              &    abs(b) < 0.001 .and. &
              &    ( abs(h-h_ref(k,j)) < my_h_tol .or. &
              &      abs(h-h_ref(k,j+1)) < my_h_tol ) ) ) then
              ! Not extrapolating in phi
              h_grid(i) = h + req_s
              p_grid(i) = q
                if ( debug ) then
                  oops=''
                  if ( abs(h+req_s-htan_r/cos(p)) > 5.0e-4 ) oops='      TRIG'
                100 format (i4,i2,f9.5,f9.3,f11.3,f12.3,f9.3, &
                      & 1p,g14.6,i3,11g10.2)
                  write (*, 100, advance='no') &
                    & i, j, p_grid(i), h_grid(i), &
                    & h_ref(k,j)+req_s, h_ref(k,j+1)+req_s, &
                  htan_r / cos(p), &
                  !htan_r * ( 1.0_rp + p2*((c2+p2*(c4+p2*c6))) ), & ! ~ htan_r * sec(p)
                  !a, b, c, &
                  d, n, dd(1:n)
                  write (*, '(2x,a)') trim(adjustl(oops))
                end if
              if ( stat(i) == grid ) then
                no_grid_fits = no_grid_fits - 1
              else
                no_bad_fits = no_bad_fits - 1
              end if
              stat(i) = good
              cycle path ! end the Newton iteration successfully
            end if
            if ( q > p_basis(j+1) .and. abs(h-h_ref(k,j+1)) > abs(d) &
              & .and. (i /= i2 .or. j /= p_coeffs-1) ) then
              ! Phi is beyond the range, and H is outside by more than the
              ! move.  Go on to the next columns of H_ref and P_Basis.
              i1 = i     ! phi is monotone; the rest are in the next column
              exit path  ! or even farther over
            end if
          end if
          p = p - d / (a*p*(d1+p2*(d3+p2*d5)) + b) ! do the Newton step
        end do ! n
        ! Newton iteration failed to converge.  Use the break point if the
        ! result is very near to it.  We'll probably find it in the next
        ! panel, but just in case....
        q = p + phi_offset(i)
          if ( debug ) then
            write (*,100, advance='no') i, j, q, h+req_s, h_ref(k,j)+req_s, h_ref(k,j+1)+req_s
            if ( abs(d) < 0.001_rp ) then
              write (*,'(a,i3,1p,3g10.2)') " Newton out of bounds  ", n, dd(1:min(n,3))
            else
              write (*,'(a,i3,1p,3g10.2)') " Newton didn't converge", n, dd(1:min(n,3))
            end if
          end if
        if ( stat(i) /= grid ) then
          if ( abs(h-h_ref(k,j)) <  0.1 .and. &
            & abs(q-p_basis(j)) < abs(q-p_basis(j+1)) ) then
            h_grid(i) = h_ref(k,j)+req_s
            p_grid(i) = p_basis(j)
            stat(i) = grid
            no_grid_fits = no_grid_fits + 1
            no_bad_fits = no_bad_fits - 1
          else if ( abs(h-h_ref(k,j+1)) <  0.1 .and. &
            & abs(q-p_basis(j+1)) < abs(q-p_basis(j)) ) then
            h_grid(i) = h_ref(k,j+1)+req_s
            p_grid(i) = p_basis(j+1)
            stat(i) = grid
            no_grid_fits = no_grid_fits + 1
            no_bad_fits = no_bad_fits - 1
          end if
        end if
        if ( q > p_basis(j+1) .and. (i /= i2 .or. j /= p_coeffs-1) ) then
          i1 = i    ! phi is monotone; the remaining solutions are in the next
          exit path ! column or even farther over, and not before the current row
        end if
      end do path ! i
    end do phi ! j

      if ( debug .and. no_bad_fits > 0 ) &
        & print *, 'no_bad_fits =', no_bad_fits, ', no_grid+fits =', no_grid_fits

    ! Since we have solved for the intersection of sec(phi) with the heights
    ! on constant-zeta surfaces, it is impossible for the heights not to be
    ! monotone increasing away from the tangent point.

    if ( do_dumps > 0 ) then
      call output ( h_tan, before='h_tan = ' )
      call output ( req, before=', req = ', advance='yes' )
      call dump ( p_grid(:n_path), name='p_grid before refractive correction', &
        & format='(1pg14.6)', clean=clean )
      call dump ( h_grid(:n_path), name='h_grid', format='(1pg14.6)', clean=clean )
      if ( dump_stop > 0 ) stop
    end if

  end subroutine Height_Metrics

  ! -----------------------------------------------  More_Metrics  -----
  subroutine More_Metrics ( &
          ! Inputs:
          & phi_t, tan_ind, n_tan, p_basis, z_ref, n_ref, t_ref,     &
          & h_grid, dhidzij, refract,                                &
          ! Inout:
          & p_grid,                                                  &
          ! Outputs:
          & t_grid, dhitdzi,                                         &
          ! Optional inputs:
          & ddhidhidtl0, dhidtlm, t_deriv_flag, z_basis,             &
          ! Optional outputs:
          & ddhtdhtdtl0, dhitdtlm, dhtdtl0, dhtdzt,                  &
          & do_calc_hyd, do_calc_t, eta_zxp, tan_phi_t )

    ! This subroutine computes metrics-related things after H_Grid and
    ! P_Grid are computed by Height_Metrics, and then perhaps augmented
    ! with, for example, the minimum Zeta point.

    use Dump_0, only: Dump
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: RP, IP
    use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: phi_t      ! orbit projected tangent geodetic angle
    integer(ip), intent(in) :: tan_ind ! tangent height index, 1 = center of
    !                                     longest path
    integer(ip), intent(in) :: n_tan   ! tangent index in path, usually n_path/2
    real(rp), intent(in) :: p_basis(:) ! horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_ref(:)   ! -log pressures (zetas) for which
    !                                     heights/temps are needed.  Only the
    !                                     parts from the tangent outward are used
    real(rp), intent(in) :: n_ref(:,:) ! Indices of refraction by t_phi_basis
    real(rp), intent(in) :: t_ref(:,:) ! temperatures by t_phi_basis
    real(rp), intent(in) :: h_grid(:)  ! computed heights, referenced to Earth center
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative by t_phi_basis
    logical,  intent(in) :: refract    ! compute phi refractive correction

    ! Inout:

    real(rp), intent(inout) :: p_grid(:)  ! phi's on the path

    ! outputs:

    real(rp), intent(out) :: t_grid(:) ! computed temperatures
    real(rp), intent(out) :: dhitdzi(:)! derivative of height wrt zeta
    !                                    --may be useful in future computations

    ! optional inputs

    real(rp), optional, intent(in) :: ddhidhidtl0(:,:,:) ! second order reference
    !   temperature derivatives. This is (height, phi_basis, zeta_basis).
    !   Needed only if present(dhidtlm).
    real(rp), optional, intent(inout) :: dhidtlm(:,:,:) ! reference temperature
    !   derivatives. This gets adjusted so that at ref_h(1,@tan phi)) is 0.0 for
    !   all temperature coefficients. This is height X zeta_basis X phi_basis
    logical, optional, intent(in) :: t_deriv_flag(:)  ! User's deriv. flags for
    !   Temperature. needed only if present(dhidtlm).
    real(rp), optional, intent(in) :: z_basis(:) ! vertical temperature basis
    !   Needed only if present(dhidtlm).

    ! optional outputs.

    real(rp), optional, intent(out) :: ddhtdhtdtl0(:) ! Second order 
    !          derivative at the tangent only---used for antenna affects.
    !          Computed if present(dhidtlm)
    real(rp), optional, intent(out) :: dhitdtlm(:,:)
    !                             derivative of path position wrt temperature
    !                             statevector (z_basis X phi_basis)
    real(rp), optional, intent(out) :: dhtdtl0(:)  ! First order derivative at
    !           the tangent.  Computed if present(dhidtlm)
    real(rp), optional, intent(out) :: dhtdzt      ! height derivative wrt
    !                                                pressure at the tangent
    logical, optional, intent(out) :: do_calc_hyd(:,:) ! nonzero locator for
    !           hydrostatic calculations.  Computed if present(dhidtlm)
    logical, optional, intent(out) :: do_calc_t(:,:) ! nonzero locater for
    !           temperature bases computations.  Computed if present(dhidtlm)
    real(rp), optional, intent(out) :: eta_zxp(:,:) ! eta matrix for temperature
    !           Computed if present(dhidtlm)
    real(rp), optional, intent(out) :: tan_phi_t ! temperature at the tangent

    ! Local variables.

    integer :: Do_Dumps    ! 0 = no dump, >0 = dump
    integer :: Dump_Stop   ! 0 = no dump, >0 = dump/stop
    integer :: First, Last ! Nonzeros in Eta_T
    integer :: H_Phi_Dump  ! 0 = no dump, >0 = dump
    integer :: I           ! Subscript, loop inductor
    integer :: N_PATH      ! Path length = 2*(size(z_ref)+1-tan_ind)
    integer :: N_VERT      ! size(z_ref)
    integer :: VERT_INDS(2*(size(z_ref)+1-tan_ind)) ! What to use in z_ref

    logical :: NOT_ZERO_P(size(vert_inds),size(p_basis))

    real(rp) :: ETA_P(size(vert_inds),size(p_basis))
    real(rp) :: ETA_T(size(p_basis))
    real(rp) :: N_GRID(size(vert_inds))     ! index of refraction
    real(rp) :: PHI_CORR(size(vert_inds))   ! the refractive correction
    real(rp) :: PHI_SIGN(size(vert_inds))   ! +/- 1.0

    ! For debugging output format:
    logical, parameter :: clean = .false.

!   It would be nice to do this the first time only, but the
!   retrieve command in the L2CF can now change switches
!   if ( do_dumps < 0 ) then ! First time only
      dump_stop = index(switches,'metD')
      do_dumps = max(dump_stop,index(switches,'metd'))
      h_phi_dump = index(switches,'hphi')
!   end if

    n_vert = size(z_ref)

    ! Only use the parts of the reference grids that are germane to the
    ! present path

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (i,  i=tan_ind,n_vert) /)
    n_path = size(vert_inds)

    ! sign of phi vector
    phi_sign = (/ (-1, i=n_vert,tan_ind,-1), (+1, i=tan_ind,n_vert) /)

    ! Interpolate the index of refraction (N_Ref) to the path (N_Grid)
    ! and correct phi for refraction.
    if ( refract ) then
      call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p(:,:) )
      do i = 1, n_path
        n_grid(i) = dot_product(n_ref(vert_inds(i),:), eta_p(i,:))
      end do
      call phi_refractive_correction ( n_tan, n_grid, h_grid(:n_path), phi_corr )
      p_grid(:n_path) = p_grid(:n_path) + phi_sign * phi_corr
      if ( h_phi_dump > 0 ) &
        & call dump ( phi_corr, format='(1pg14.6)', &
        &             name='Refractive correction', clean=clean )
      if ( do_dumps > 0 ) &
        & call dump ( p_grid(:n_path), name='p_grid after refractive correction', &
          & format='(1pg14.6)', clean=clean )
    end if

    ! Interpolate Temperature (T_Ref) and the vertical height derivative
    ! (dhidzij) to the path (T_Grid and dhitdzi).

    call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p, NOT_ZERO = not_zero_p )
    do i = 1, n_path
      t_grid(i) = dot_product(t_ref(vert_inds(i),:), eta_p(i,:))
      ! compute the vertical derivative grid
      dhitdzi(i) = dot_product(dhidzij(vert_inds(i),:), eta_p(i,:))
    end do

    ! now for the optional tangent quantities.
    if ( present(tan_phi_t) .or. present(dhtdzt) .or. present(dhidtlm) ) &
      & call get_eta_sparse ( p_basis, phi_t, eta_t, first, last )
    if ( present(tan_phi_t) ) &
      & tan_phi_t = dot_product(t_ref(tan_ind,first:last),eta_t(first:last))
    if ( present(dhtdzt) ) &
      & dhtdzt = dot_product(dhidzij(tan_ind,first:last),eta_t(first:last))

    ! compute tangent temperature derivatives
    if ( present(dhidtlm) ) call Tangent_Temperature_Derivatives

    if ( do_dumps > 0 ) then
      call dump ( t_grid(:n_path), name='t_grid', format='(1pg14.6)', clean=clean )
      call dump ( dhitdzi(:n_path), name='dhitdzi', format='(1pg14.6)', clean=clean )
      if ( dump_stop > 0 ) stop
    end if

  contains

    subroutine Tangent_Temperature_Derivatives
      ! This is a subroutine instead of inline so that the references
      ! to Z_Basis in the dimensions of automatic variables are only
      ! attempted if Z_Basis is present.
      real(rp) :: ETA_T2(size(vert_inds),size(z_basis))
      logical :: NOT_ZERO_T(size(vert_inds),size(z_basis))
      integer :: I, J, SV_P, SV_T, SV_Z ! Loop inductors and subscripts
      integer :: P_COEFFS    ! size(p_basis)
      integer :: Z_COEFFS    ! size(z_basis)

      ! Adjust the 2d hydrostatic relative to the surface. Even though
      ! this is updated on every invocation, that is, with a new phi_t, it
      ! works as if the original value were updated with the current
      ! phi_t, because the interpolation represented by eta_t is linear.
      ! Thus, the effect of cumulative updates for each new phi_t are the
      ! same as starting from the original dhidtlm and updating with the
      ! latest phi_t.  The algebra is horrible, but Maple has verified this.
      p_coeffs = size(p_basis)
      z_coeffs = size(z_basis)
      do i = 1, z_coeffs
        dhidtlm(:,i,:) = dhidtlm(:,i,:) - &
          & dot_product(dhidtlm(1,i,first:last),eta_t(first:last))
      end do

      j = z_coeffs * p_coeffs
      dhtdtl0 = RESHAPE(dhidtlm(tan_ind,:,:) * SPREAD(eta_t,1,z_coeffs),&
                     & (/j/))

      ddhtdhtdtl0 = RESHAPE( &
                   ddhidhidtl0(tan_ind,:,:) * SPREAD(eta_t,1,z_coeffs), &
                     & (/j/))

      ! compute the path temperature noting where the zeros are
      call get_eta_sparse ( z_basis, z_ref(vert_inds), eta_t2, NOT_ZERO = not_zero_t )

      sv_t = 0
      do sv_p = 1 , p_coeffs
        do sv_z = 1 , z_coeffs
          sv_t = sv_t + 1
          if ( t_deriv_flag(sv_t) ) then
            do_calc_t(:,sv_t) = not_zero_t(:,sv_z) .and. not_zero_p(:,sv_p)
            where ( do_calc_t(:,sv_t) )
              eta_zxp(:,sv_t) = eta_t2(:,sv_z) * eta_p(:,sv_p)
            elsewhere
              eta_zxp(:,sv_t) = 0.0
            end where
            do_calc_hyd(:,sv_t) = not_zero_p(:,sv_p) .and. &
              &                   dhidtlm(vert_inds(:),sv_z,sv_p) > 0.0_rp
            where ( do_calc_hyd(:,sv_t) )
              dhitdtlm(:,sv_t) = dhidtlm(vert_inds(:),sv_z,sv_p) * eta_p(:,sv_p)
            elsewhere
              dhitdtlm(:,sv_t) = 0.0
            end where
          else
            do_calc_t(:,sv_t) = .false.
            eta_zxp(:,sv_t) = 0.0
            do_calc_hyd(:,sv_t) = .false.
            dhitdtlm(:,sv_t) = 0.0
          end if
        end do
      end do

    end subroutine Tangent_Temperature_Derivatives

  end subroutine More_Metrics

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Metrics_m

! $Log$
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
