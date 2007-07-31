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

  use MLSKinds, only: RP

  implicit NONE
  private
  public :: Height_Metrics, More_Metrics, Solve_H_Phi, Tangent_Metrics
  public :: More_Points

  ! Values for stat argument to/from Solve_H_Phi:
  integer, parameter, public :: No_sol = 0  ! No solution
  integer, parameter, public :: Complex = 1 ! Complex start for Newton iteration
  integer, parameter, public :: Good = 2    ! Newton converged, tangent point, extrapolated
  integer, parameter, public :: Grid = 3    ! Close to a grid point

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Private stuff
  real(rp), parameter :: DefaultTol = 0.001_rp ! km, for Solve_h_phi
  ! To control debugging
  logical, parameter :: Debug = .false.
  logical, parameter :: NewtonDetails = .true. .and. debug
  logical, parameter :: ComplexDebug = .true. .and. debug
  ! For debugging output format:
  logical, parameter :: Clean = .false.

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
    use GLNP, only: NGP1
    use Make_Z_Grid_m, only: Default_Thresh ! Threshold for rejecting duplicates
    use MLSKinds, only: RP, R8
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
    real(rp) :: CP2, SP2   ! Cos^2 phi_t, Sin^2 phi_t
    real(rp) :: ETA_T(size(p_basis)) ! Interpolating coefficients
    real(rp) :: H_T1, H_T2 ! Used to determine where H_Tan is in H_Ref
    real(rp) :: R          ! H_T2/(H_T1-H_T2), used for interpolating in Z_Ref
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
      h_surf = dot_product(h_ref(1,first:last), eta_t(first:last))
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
      r = h_t2/(h_t2-h_t1)
      z_tan = z_ref(tan_ind_c) * r + (r-1.0) * z_ref(tan_ind_c+1)
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
      if ( do_dumps > 0 ) call output ( z_tan, before='Z_Tan = ', after=', ' )
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
  subroutine Height_Metrics (  phi_t, tan_ind, p_basis, h_ref, req, h_surf, &
                            &  h_tan, z_ref,                                &
                            ! Outputs:
                            &  vert_inds, h_grid, p_grid,                   &
                            ! Optional inputs:
                            &  h_tol )

    !{ This subroutine computes {\tt h\_grid} and {\tt p\_grid} that define
    !  a 2-d integration path.  These points are at the intersections of the
    !  line of sight given by $H = H_\text{tan} \sec \phi$ and the heights of
    !  surfaces of constant $\zeta$ represented in piecewise linear segments
    !  by consecutive elements of each row of {\tt H\_ref}.

    use Dump_0, only: Dump
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: RP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use Output_m, only: OUTPUT
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: Phi_t      ! Orbit projected tangent geodetic angle
    integer, intent(in) :: Tan_ind     ! Tangent height index, 1 = center of
    !                                     longest path
    real(rp), intent(in) :: P_basis(:) ! Horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: H_ref(:,:) ! Heights by z_ref and p_basis
    real(rp), intent(in) :: Req        ! Equivalent elliptical earth radius
    real(rp), intent(in) :: H_Surf     ! Height at the Earth surface
    real(rp), intent(in) :: H_Tan      ! Tangent height above H_Surf -- negative
    !                                     for Earth-intersecting ray
    real(rp), intent(in) :: Z_ref(:)   ! -log pressures (zetas) for which
    !                                     heights/temps are needed.  Only used
    !                                     where H/Phi iteration fails.

    ! outputs:
    integer, intent(out) :: Vert_Inds(:) ! What to use in h_ref, 1:n_path
    real(rp), intent(out) :: H_grid(:)   ! computed heights, referenced to Earth center
    real(rp), intent(out) :: P_grid(:)   ! computed phi's

    ! optional inputs
    real(rp), optional, intent(in) :: H_Tol ! Height tolerance in kilometers
                                       !   for convergence of phi/h iteration

    ! Local variables.
    integer :: Do_Dumps    ! 0 = no dump, >0 = dump
    integer :: Dump_Stop   ! 0 = no dump, >0 = dump/stop
    integer :: H_Phi_Dump  ! 0 = no dump, >0 = dump
    integer :: H_Phi_Stop  ! 0 = no dump, >0 = dump/stop
    integer :: I, I1, I2, IBAD, J, K
    integer :: NO_BAD_FITS
    integer :: NO_GRID_FITS
    integer :: N_PATH      ! Path length = 2*(n_vert+1-tan_ind)
    integer :: N_Tan       ! Tangent index in path, N_Path/2
    integer :: N_Vert      ! Size(H_ref,1)
    integer :: P_COEFFS    ! Size(P_basis)

    integer :: Stat(size(vert_inds))

    logical :: Outside     ! Solve_H_Phi converged outside the allowed area
    real(rp) :: A, B, C    ! Polynomial coefficients used to solve for phi and H
    real(rp) :: DP         ! p_basis(j+1) - p_basis(j)
    real(rp) :: DPJ0       ! p_basis(j  )-phi_t
    real(rp) :: DPJ1       ! p_basis(j+1)-phi_t
    real(rp) :: H          ! Tentative solution for H
    real(rp) :: HTAN_R     ! H_Tan + req
    real(rp) :: My_H_Tol   ! Tolerance in kilometers for height convergence
    real(rp) :: P          ! Tentative solution for phi
    real(rp) :: REQ_S      ! Req - H_Surf

    real(rp) :: ETA_T(size(p_basis)) ! Interpolating coefficients
    real(rp) :: PHI_OFFSET(size(vert_inds)) ! PHI_T or a function of NEG_H_TAN
    real(rp) :: PHI_SIGN(size(vert_inds))   ! +/- 1.0

    ! For debugging
    character(4), parameter :: nStat(no_sol:grid) = &
      & (/ 'none', 'cplx', 'good', 'grid' /)

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
    n_vert = size(h_ref,1)

    my_h_tol = defaultTol ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    p_coeffs = size(p_basis)

    ! p_basis and p_grid are phi's in offset radians relative to phi_t, that
    ! is, the phi_t, p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t

    phi_offset(:n_tan) = phi_t
    if ( h_tan < 0.0 ) then ! Earth-intersecting ray
      phi_offset(n_tan+1:) = phi_t - 2.0_rp*Acos((req+h_tan)/req)
    else
      phi_offset(n_tan+1:) = phi_t
    end if

    req_s = req - h_surf
    htan_r = h_tan + req

    if ( h_phi_dump > 0 ) call dumpInput ( tan_ind )

      if ( debug ) &
        & print '(a)', &
        & '   I J   K                HREF(K,J) HREF(K,J+1)    A        B    phi sign    C', &
        & '          phi    H_GRID   HREF(K,J) HREF(K,J+1)  H(TRIG)    D           N   Diffs'

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
          print 3, i1, 0, p_grid(i1), h_grid(i1), '  Below p_basis'
        3 format ( i4,i2,f9.5,f9.3,a )
        end if
    end do
    do i2 = n_path, 1, -1
      if ( stat(i2) == good ) cycle ! Skip the tangent point
      k = vert_inds(i2)
      h = max(htan_r,h_ref(k,p_coeffs)+req_s)
      p = acos(htan_r/h) * phi_sign(i2) + phi_offset(i2)
      if ( p <= p_basis(p_coeffs) ) exit ! phi is monotone
      p_grid(i2) = p
      h_grid(i2) = h
      no_bad_fits = no_bad_fits - 1
      stat(i2) = good
        if ( debug ) &
          & print 3, i2, p_coeffs, p_grid(i2), h_grid(i2), '  Above p_basis'
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

    ibad = i2 ! Index of row of first failure in a column
phi:do j = 1, p_coeffs-1
      dp = p_basis(j+1)-p_basis(j)
      a = dp * htan_r ! Actually 2*A
      dpj1 = p_basis(j+1)-phi_t
      dpj0 = p_basis(j)-phi_t
path: do i = i1, i2
        if ( stat(i) == good ) cycle ! probably the tangent point
        k = vert_inds(i)
        b = -(h_ref(k,j+1) - h_ref(k,j)) ! h_surf cancels here
        c = -(h_ref(k,j)-h_surf)   * dpj1 &
          & +(h_ref(k,j+1)-h_surf) * dpj0 &
          & + dp * h_tan
          if ( newtonDetails ) &
            & print 120, i, j, k, h_ref(k,j)+req_s, h_ref(k,j+1)+req_s, &
            & a, b, int(phi_sign(i)), c
          120 format ( i4, i2, i4, 14x, f11.3,f12.3, f9.3, 1p,g14.6, i3, g14.6 )
        if ( stat(i) == grid ) then
          no_grid_fits = no_grid_fits - 1
          no_bad_fits = no_bad_fits + 1
        end if
        call Solve_H_Phi ( p_basis(j:j+1), phi_offset(i), phi_sign(i), &
          &                h_ref(k,j:j+1), a, b, c, &
          &                htan_r, req_s, my_h_tol, i /= i2 .or. j /= p_coeffs-1, &
          &                h_grid(i), p_grid(i), stat(i), outside )
        if ( stat(i) == good ) then
          no_bad_fits = no_bad_fits - 1
        else if ( stat(i) == grid ) then
          no_grid_fits = no_grid_fits + 1
          no_bad_fits = no_bad_fits - 1
        end if
        if ( outside ) then
           ! phi is monotone; the remaining solutions are in the next
           ! column or even farther over, and not before the current row
           ! or the oldest "bad" row in this column.
          i1 = min(i,ibad)   
          ibad = i2
          exit
        end if
        if ( stat(i) == good ) cycle ! Newton iteration ended successfully
        ibad=min(i,ibad)
        if ( (debug .or. complexDebug) .and. stat(i) == complex ) then ! Complex solution
            if ( debug ) &
              & print '(i4,i2,18x,f11.3,f12.3,1x,a)', &
                & i, j, h_ref(k,j)+req_s, h_ref(k,j+1)+req_s, 'Complex solution'
            if ( complexDebug ) then
              call output ( req, before='&in Req = ' )
              call output ( h_surf, before=', H_Surf = ' )
              call output ( h_tan, before=', H_Tan = ', after=',', advance='yes' )
              call output ( phi_t, before='  Phi_T = ' )
              call output ( tan_ind, before=', tan_ind = ', advance='yes' )
              call dump ( h_ref, name='h_ref' )
              call dump ( p_basis, name='p_basis' )
            end if
        end if
      end do path ! i
    end do phi ! j

    if ( no_bad_fits > 0 ) then
        if ( debug ) then
          print *, 'no_bad_fits =', no_bad_fits, ', no_grid+fits =', no_grid_fits
          do i = 1, size(stat), 10
            print "(i4,'#',10(1x,a:))", i, nStat(stat(i:min(size(stat),i+9)))
          end do
        end if
      call dumpInput ( 1 )
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "Height_Metrics failed to find H/Phi solution for some path segment" )
      if ( index(switches,'MHPX') /= 0 ) then
        call dump ( stat, 'STAT' )
        call dump ( h_grid, 'H_Grid before 1d fixup' )
        call dump ( p_grid, 'P_Grid before 1d fixup' )
      end if
      ! We shouldn't get here at all, so don't worry about efficiency
      if ( stat(1) /= good .and. stat(1) /= grid .or. &
        &  stat(n_path) /= good .and. stat(n_path) /= grid ) then
        call MLSMessage ( MLSMSG_Warning, moduleName, 'Resorting to 1d' )
        call get_eta_sparse ( p_basis, phi_t, eta_t )
        do i1 = 1, n_path
          if ( stat(i1) /= good .and. stat(i1) /= grid ) then
            h_grid(i1) = dot_product(h_ref(vert_inds(i1),:),eta_t)
          end if
        end do
        ! Make sure H is monotone increasing away from the tangent point
        i1 = n_tan - 1
        i2 = 1
        do k = -1, 1, 2
          do j = i1, i2, k
            if ( h_grid(j) <= h_grid(j-k) ) then
              do i = j-1, 1, -1
                if ( h_grid(j-k) < h_grid(i) ) exit
              end do
              if ( i < 1 ) then
                h_grid(j) = h_grid(j-k) + 1.0 ! use SWAG method
              else
                h_grid(j) = h_grid(j-k) + (h_grid(i)-h_grid(j-k)) * &
                  & ( z_ref(vert_inds(j)) - z_ref(vert_inds(j-k))) / &
                  & ( z_ref(vert_inds(i)) - z_ref(vert_inds(j-k)))
              end if
            end if
          end do
          i1 = n_tan + 2
          i2 = n_path
        end do
      else
        do j = n_tan-1, 1, -1
          do i2 = j, 1, -1
            if ( stat(i2) == good .or. stat(i2) == grid ) exit
          end do
          do i1 = j, n_tan
            if ( stat(i1) == good .or. stat(i1) == grid ) exit
          end do
          h_grid(j) = h_grid(i1) + (h_grid(i2)-h_grid(i1)) * &
            & ( z_ref(vert_inds(j)) - z_ref(vert_inds(i1))) / &
            & ( z_ref(vert_inds(i2)) - z_ref(vert_inds(i1)))
        end do
        do j = n_tan+2, n_path
          do i2 = j, n_path
            if ( stat(i2) == good .or. stat(i2) == grid ) exit
          end do
          do i1 = j, n_tan+1, -1
            if ( stat(i1) == good .or. stat(i1) == grid ) exit
          end do
          h_grid(j) = h_grid(i1) + (h_grid(i2)-h_grid(i1)) * &
            & ( z_ref(vert_inds(j)) - z_ref(vert_inds(i1))) / &
            & ( z_ref(vert_inds(i2)) - z_ref(vert_inds(i1)))
        end do
      end if
      do i = 1, n_path
        if ( stat(i) /= good .and. stat(i) /= grid ) &
          & p_grid(i) = acos(h_grid(n_tan)/h_grid(i))
      end do
      if ( index(switches,'MHPX') /= 0 ) then
        call dump ( h_grid, name='H_Grid after 1d fixup' )
        call dump ( p_grid, name='P_Grid after 1d fixup' )
        call MLSMessage ( MLSMSG_Error, moduleName, 'Halt requested by MHPX' )
      end if
    end if

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

  contains

    subroutine DumpInput ( FirstRow )
      integer, intent(in) :: FirstRow
      call output ( h_tan, before='h_tan = ' )
      call output ( req, before=', req = ', advance='yes' )
      call dump ( p_basis, name='p_basis', format='(f14.8)', clean=clean )
      call output ( phi_t, before='phi_t = ', format='(f14.8)' )
      call output ( h_surf, before =', h_surf = ' )
      call output ( tan_ind, before=', tan_ind = ', advance='yes' )
      call dump ( phi_offset, name='phi_offset', clean=clean )
      call dump ( h_ref(firstRow:n_vert,:)+req_s, name='h_ref+req_s', format='(f14.7)',clean=clean )
    end subroutine DumpInput

  end subroutine Height_Metrics

  ! ------------------------------------------------  Solve_H_Phi  -----

  !{Solve for an intersection of $H = H_t \sec \phi$ and the line from
  ! $(\phi_1,h_1)$ to $(\phi_2,h_2)$, where $\phi_1$ and $\phi_2$ are
  ! {\tt p\_basis(1)} and {\tt p\_basis(2)}, $h_1$ and $h_2$ are {\tt
  ! h\_ref(1)} and {\tt h\_ref(2)}.
  !
  ! The solution is only acceptable if $\phi_1 \leq \phi \leq \phi_2$ and
  ! $(h_1 - h) (h_2-h) \leq 0$ or $|h_1 - h| \leq \tau$ or $|h_2 - h| \leq
  ! \tau$, where $\phi$ is {\tt p\_grid}, $h$ is {\tt h\_grid}, and $\tau$
  ! is the height tolerance.

  subroutine Solve_H_Phi ( & ! inputs
    &                      p_basis, phi_offset, phi_sign, h_ref, a, b, c, &
    &                      htan_r, req_s, tol, inside, &
                             ! outputs and inouts
    &                      h_grid, p_grid, stat, outside )

    use MLSKinds, only: RP
    use Output_m, only: OUTPUT

    real(rp), intent(in) :: P_Basis(:) ! Coordinates for H_Ref.  Actually 1:2,
                                       ! but we don't want to force copy-in if
                                       ! the actual argument isn't contiguous.
    real(rp), intent(in) :: Phi_Offset ! Offset from tangent point
    real(rp), intent(in) :: Phi_Sign   ! Which way from tangent?
    real(rp), intent(in) :: H_Ref(:)   ! Height reference.  Actually 1:2,
                                       ! but we don't want to force copy-in if
                                       ! the actual argument isn't contiguous.
    real(rp), intent(in) :: A          ! (p_basis(2)-p_basis(1)) * Htan_r
    real(rp), intent(in) :: B          ! -(h_ref(2) - h_ref(1))
    real(rp), intent(in) :: C          ! -(h_ref(1)-h_surf) * (p_basis(2)-phi_t)
                                       ! +(h_ref(2)-h_surf) * (p_basis(1)-phi_t)
                                       ! +(p_basis(2)-p_basis(1)) * h_tan
      ! A, B and C are coefficients of a quadratic approximation to the solution
    real(rp), intent(in) :: Htan_r     ! H_Tan + req
    real(rp), intent(in) :: Req_s      ! Req - H_Surf
    real(rp), intent(in) :: Tol        ! Height tolerance for Newton convergence
    logical, intent(in) :: Inside      ! P_Basis(2) is not the last one in the grid
    real(rp), intent(inout) :: H_Grid  ! H solution, inout in case there is none
    real(rp), intent(out) :: P_Grid    ! Phi solution
    integer, intent(inout) :: Stat     ! "good" or "grid" or "complex" or unchanged
    logical, intent(out) :: Outside    ! Newton iteration converged to a point
                                       ! outside p_basis(1:2)

    real(rp) :: D     ! B^2 - 4 a c
    real(rp) :: H     ! htan_r * ( 1.0_rp + secM1 ) - req_s,
                      ! ~ htan_r * sec(p) - req + h_surf
    integer :: N      ! Newton iteration counter
    integer, parameter :: NMax = 10 ! Maximum number of Newton iterations
    real(rp) :: P, P2 ! Candidate solution, p^2
    real(rp) :: Secm1 ! Sec(phi) - 1

    ! Coefficients in expansion of Sec(phi)-1 = c2*phi^2 + c4*phi^4 ...
    real(rp), parameter :: C2 = 0.5_rp, C4 = 5.0_rp/24, C6 = 61.0_rp/720.0
    real(rp), parameter :: C8 = 277.0_rp/8064
    ! Coefficients in expansion of Sec(phi)*Tan(phi) = d/dPhi(sec(phi)-1)
    real(rp), parameter :: D1 = 2*c2, D3 = 4*c4, D5 = 6*c6 ! ... 2n * c_2n
    real(rp), parameter :: D7 = 8*c8

    ! For debugging output
    real(rp) :: DD(10)

    outside = .false.
    d = b*b - 2.0 * a * c ! Actually b^2-4ac, since a is actually 2*A
    if ( d < 0.0 ) then ! Complex solution, use SWAG method
      ! Pretend D is zero
      p = -b/a
      p = min(max(p,p_basis(1)-phi_offset),p_basis(2)-phi_offset)
      stat = complex ! In case Newton iteration doesn't converge
        if ( debug ) call output ( p + phi_offset, format='(f15.5)', &
          & after='  Complex starting point', advance='yes' )
    else
      p = (-b + phi_sign * sqrt(d) ) / a
      if ( p < p_basis(1)-phi_offset ) &
        & p = 0.5 * (p_basis(1)+p_basis(2)) - phi_offset
    end if
    p_grid = p + phi_offset
    !{Use P in Newton iterations with an expansion for $\sec(\phi)$ that is
    ! higher order than two.  We don't use $\sec(\phi)-1$ because it suffers
    ! cancellation for small $\phi$.  Rather, use more terms of its Taylor
    ! series than we used for the quadratic approximation.
    do n = 1, nMax
      p2 = p**2
      secM1 = p2*(c2+p2*(c4+p2*(c6+p2*c8))) ! ~ sec(p)-1 to eighth order
      d = a*secM1 + b*p + c
        if ( debug ) dd(n) = d
      h = htan_r * ( 1.0_rp + secM1 ) - req_s ! ~ htan_r * sec(p) - req + h_surf
      if ( abs(d) < tol ) then ! difference is small enough
        if ( (p_grid-p_basis(1)) * (p_grid-p_basis(2)) <= 0.0 .and. &
               ! Converged within bounds?
                &  ( (h-h_ref(1))*(h-h_ref(2)) <= 0.0 .or. &
               ! or the difference in reference heights is within tolerance and
               ! the current H is outside the bounds by less than tolerance
          &    abs(b) < tol .and. &
          &    ( abs(h-h_ref(1)) < tol .or. &
          &      abs(h-h_ref(2)) < tol ) ) ) then
          ! Not extrapolating in phi
          h_grid = h + req_s
          stat = good
          if ( debug ) call debug1 ( 1, n )
          return
        end if
        if ( p_grid > p_basis(2) .and. abs(h-h_ref(2)) > abs(d) &
          & .and. inside ) then
          ! Phi is beyond the range, and H is outside by more than the
          ! move.  Go on to the next columns of H_ref and P_Basis.
          outside = .true. ! phi is monotone; the rest are in the next column
          if ( debug ) then
            h_grid = h + req_s
            call debug1 ( 1, n, 'OUTSIDE' )
          end if
          return           ! or even farther over
        end if
      end if
        if ( newtonDetails ) then
          h_grid = h + req_s
          call debug1 ( n, n, 'STEP' )
        end if
      p = p - d / (a*p*(d1+p2*(d3+p2*(d5+p2*d7))) + b) ! do the Newton step
      p_grid = p + phi_offset
    end do ! n
    ! Newton iteration failed to converge.  Use the break point if the
    ! result is very near to it.  We'll probably find it in the next
    ! panel, but just in case....
      if ( debug ) then
        if ( stat /= good .and. .not. outside ) then
          h_grid = h + req_s
          call debug1 ( max(n-3,1), n-1, merge('BOUNDS  ', 'CONVERGE', abs(d) < tol) )
        else
          call debug1 ( max(n-3,1), n-1 )
        end if
      end if
    if ( stat /= grid ) then
      if ( abs(h-h_ref(1)) <  0.1 .and. &
        & abs(p_grid-p_basis(1)) < abs(p_grid-p_basis(2)) ) then
        h_grid = h_ref(1)+req_s
        p_grid = p_basis(1)
        stat = grid
      else if ( abs(h-h_ref(2)) <  0.1 .and. &
        & abs(p_grid-p_basis(2)) < abs(p_grid-p_basis(1)) ) then
        h_grid = h_ref(2)+req_s
        p_grid = p_basis(2)
        stat = grid
      end if
    end if
    if ( p_grid > p_basis(2) .and. inside ) then
      outside = .true. ! phi is monotone; the rest are in the next column
      return           ! or even farther over
    end if

  contains
    subroutine debug1 ( k, l, why )
      integer, intent(in) :: k, l
      character(len=*), intent(in), optional :: why
      character(merge(13,0,debug)) :: OOPS
      if ( debug ) then
        oops=''
        if ( abs(h+req_s-htan_r/cos(p)) > 5.0e-4 ) oops='TRIG'
        if ( present(why) ) oops(6:) = why
      100 format (f15.5,f9.3,f11.3,f12.3,f9.3,1p,g14.6,i3,11g10.2)
        write (*, 100, advance='no') p_grid, h_grid, &
          & h_ref(1)+req_s, h_ref(2)+req_s, &
          & htan_r / cos(p), &
        ! & htan_r * ( 1.0_rp + p2*((c2+p2*(c4+p2*c6))) ), & ! ~ htan_r * sec(p)
          & d, n, dd(k:l)
        write (*, '(2x,a)') trim(adjustl(oops))
      end if
    end subroutine
  end subroutine Solve_H_Phi

  ! -----------------------------------------------  More_Metrics  -----
  subroutine More_Metrics ( &
          ! Inputs:
          & tan_ind, n_tan, p_basis, vert_inds, t_ref, dhidzij, p_grid, &
          ! Outputs:
          & t_grid, dhitdzi,                                            &
          ! Optional inputs:
          & ddhidhidtl0, dhidtlm, t_deriv_flag, z_basis, z_ref,         &
          ! Optional outputs:
          & ddhtdhtdtl0, dhitdtlm, dhtdtl0, dhtdzt,                     &
          & do_calc_hyd, do_calc_t, eta_zxp, nz_zxp, nnz_zxp, tan_phi_t )

    ! This subroutine computes metrics-related things after H_Grid and
    ! P_Grid are computed by Height_Metrics, and then perhaps augmented
    ! with, for example, the minimum Zeta point.

    use Dump_0, only: Dump
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse, Multiply_Eta_Column_Sparse
    use MLSKinds, only: RP
    use Toggles, only: Switches
    ! inputs:

    integer, intent(in) :: tan_ind     ! tangent height index, 1 = center of
    !                                     longest path
    integer, intent(in) :: n_tan       ! tangent index in path, usually n_path/2
    real(rp), intent(in) :: p_basis(:) ! horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: t_ref(:,:) ! temperatures at z_ref X p_basis
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative at z_ref X p_basis
    real(rp), intent(in) :: p_grid(:)  ! phi's on the path
    integer, intent(in) :: vert_inds(:) ! What to use in [zt]_ref

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
    real(rp), optional, intent(in) :: z_ref(:)   ! -log pressures (zetas) for
    !   which derivatives are needed.  Only the parts from the tangent outward
    !   are used.  Needed only if present(dhidtlm).

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
    real(rp), optional, intent(inout) :: eta_zxp(:,:) ! eta matrix for temperature
    !           Computed if present(dhidtlm)
    integer, optional, intent(inout) :: NZ_ZXP(:,:)   ! Nonzeros in Eta_zxp
    integer, optional, intent(inout) :: NNZ_ZXP(:)    ! Numbers of rows in NZ_ZXP
    real(rp), optional, intent(out) :: tan_phi_t ! temperature at the tangent

    ! Local variables.

    integer :: Do_Dumps    ! 0 = no dump, >0 = dump
    integer :: Dump_Stop   ! 0 = no dump, >0 = dump/stop
    integer :: I           ! Subscript, loop inductor
    integer :: N_PATH      ! Path length = size(vert_inds)
    integer :: N_VERT      ! size(z_ref)

    logical :: NOT_ZERO_P(size(vert_inds),size(p_basis))
    integer :: col1(size(vert_inds))        ! First nonzero in rows of Eta_P
    integer :: col2(size(vert_inds))        ! Last nonzero in rows of Eta_P
    integer :: NZ_P(size(vert_inds),size(p_basis)) ! Nonzeros in columns of Eta_P
    integer :: NNZ_P(size(p_basis))         ! Numbers of rows in NZ_P

    real(rp) :: ETA_P(size(vert_inds),size(p_basis))

!   It would be nice to do this the first time only, but the
!   retrieve command in the L2CF can now change switches
!   if ( do_dumps < 0 ) then ! First time only
      dump_stop = index(switches,'metD')
      do_dumps = max(dump_stop,index(switches,'metd'))
!   end if

    n_vert = size(t_ref,1)
    n_path = size(vert_inds)

    ! Interpolate Temperature (T_Ref) and the vertical height derivative
    ! (dhidzij) to the path (T_Grid and dhitdzi).

    nnz_p = 0
    eta_p = 0.0
    not_zero_p = .false.
    call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p, 1, n_path, &
      & nz_p, nnz_p, col1, col2, not_zero_p )
    do i = 1, n_path
      t_grid(i) = dot_product(t_ref(vert_inds(i),col1(i):col2(i)), &
        & eta_p(i,col1(i):col2(i)))
      ! compute the vertical derivative grid
      dhitdzi(i) = dot_product(dhidzij(vert_inds(i),col1(i):col2(i)), &
        & eta_p(i,col1(i):col2(i)))
    end do

    ! now for the optional tangent quantities.
    if ( present(tan_phi_t) ) &
      & tan_phi_t = dot_product(t_ref(tan_ind,col1(n_tan):col2(n_tan)), &
        & eta_p(n_tan,col1(n_tan):col2(n_tan)))
    if ( present(dhtdzt) ) &
      & dhtdzt = dot_product(dhidzij(tan_ind,col1(n_tan):col2(n_tan)), &
        & eta_p(n_tan,col1(n_tan):col2(n_tan)))

    ! compute tangent temperature derivatives
    if ( present(dhidtlm) ) call Tangent_Temperature_Derivatives ( size(z_basis) )

    if ( do_dumps > 0 ) then
      call dump ( t_grid(:n_path), name='t_grid', format='(1pg14.6)', clean=clean )
      call dump ( dhitdzi(:n_path), name='dhitdzi', format='(1pg14.6)', clean=clean )
      if ( dump_stop > 0 ) stop
    end if

  contains

    subroutine Tangent_Temperature_Derivatives ( N )
      ! This is a subroutine instead of inline so that the references
      ! to Z_Basis in the dimensions of automatic variables are only
      ! attempted if Z_Basis is present.
      integer, intent(in) :: N
      real(rp) :: ETA_T2(n_path,n)    ! n = size(z_basis) but Intel's
      integer :: NZ_T2(n_path,n)      ! Nonzeros in Eta_T2
      integer :: NNZ_T2(n)            ! Numbers of rows in NZ_T2
      integer :: I, J, SV_P, SV_T, SV_Z ! Loop inductors and subscripts
      integer :: P_COEFFS    ! size(p_basis)
      integer :: Z_COEFFS    ! size(z_basis)

      ! Adjust the 2d hydrostatic relative to the surface. Even though
      ! this is updated on every invocation, that is, with a new phi_t, it
      ! works as if the original value were updated with the current
      ! phi_t, because the interpolation represented by eta_p is linear.
      ! Thus, the effect of cumulative updates for each new phi_t are the
      ! same as starting from the original dhidtlm and updating with the
      ! latest phi_t.  The algebra is horrible, but Maple has verified this.
      p_coeffs = size(p_basis)
      z_coeffs = size(z_basis)
      do i = 1, z_coeffs
        dhidtlm(:,i,:) = dhidtlm(:,i,:) - &
          & dot_product(dhidtlm(1,i,col1(n_tan):col2(n_tan)), &
            & eta_p(n_tan,col1(n_tan):col2(n_tan)))
      end do

      j = z_coeffs * p_coeffs
      dhtdtl0 = RESHAPE(dhidtlm(tan_ind,:,:) * SPREAD(eta_p(n_tan,:),1,z_coeffs),&
                     & (/j/))

      ddhtdhtdtl0 = RESHAPE( &
                   ddhidhidtl0(tan_ind,:,:) * SPREAD(eta_p(n_tan,:),1,z_coeffs), &
                     & (/j/))

      ! compute the path temperature noting where the zeros are
      nnz_t2 = 0
      eta_t2 = 0.0
      call get_eta_sparse ( z_basis, z_ref(vert_inds), eta_t2, &
        & n_tan, 1, nz_t2, nnz_t2, .false. )
      call get_eta_sparse ( z_basis, z_ref(vert_inds), eta_t2, &
        & n_tan+1, size(vert_inds), nz_t2, nnz_t2, .true. )

      call multiply_eta_column_sparse ( eta_t2, nz_t2, nnz_t2, eta_p, nz_p, nnz_p, &
        & eta_zxp, nz_zxp, nnz_zxp, do_calc_t )

      sv_t = 0
      do sv_p = 1 , p_coeffs
        do sv_z = 1 , z_coeffs
          sv_t = sv_t + 1
          if ( t_deriv_flag(sv_t) ) then
            dhitdtlm(:,sv_t) = max(dhidtlm(vert_inds,sv_z,sv_p),0.0_rp) * eta_p(:,sv_p)
            do_calc_hyd(:,sv_t) = dhitdtlm(:,sv_t) /= 0.0_rp
          else
            do_calc_hyd(:,sv_t) = .false.
            dhitdtlm(:,sv_t) = 0.0
            eta_zxp(nz_zxp(:nnz_zxp(sv_t),sv_t),sv_t) = 0.0
            do_calc_t(nz_zxp(:nnz_zxp(sv_t),sv_t),sv_t) = .false.
            nnz_zxp(sv_t) = 0
          end if
        end do
      end do

    end subroutine Tangent_Temperature_Derivatives

  end subroutine More_Metrics

  ! ------------------------------------------------  More_Points  -----
                         ! Inputs:
  subroutine More_Points (  phi_t, tan_ind, p_basis, z_ref, h_ref, req,  &
                         &  h_surf, h_tan, p_grid,                       &
                         ! Outputs:
                         &  z_new, h_new, p_new, n_new,                  &
                         ! Optional inputs:
                         &  h_tol )

    ! Check for path crossings in h_ref from tan_ind down, which Height_Metrics
    ! doesn't do.

    use MLSKinds, only: RP
    ! inputs:

    real(rp), intent(in) :: phi_t      ! Orbit projected tangent geodetic angle
    integer, intent(in) :: tan_ind     ! Tangent height index, 1 = center of
    !                                     longest path
    real(rp), intent(in) :: p_basis(:) ! Horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_ref(:)   ! Reference zetas
    real(rp), intent(in) :: h_ref(:,:) ! Heights by z_ref and p_basis
    real(rp), intent(in) :: req        ! Equivalent elliptical earth radius
    real(rp), intent(in) :: H_Surf     ! Height at the Earth surface
    real(rp), intent(in) :: H_Tan      ! Tangent height above H_Surf -- negative
    !                                     for Earth-intersecting ray
    real(rp), intent(in) :: p_grid(:)  ! From Height_Metrics

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
    real(rp) :: Htan_r    ! H_Tan + Req
    integer :: I, J
    real(rp) :: My_H_Tol
    logical :: None
    logical :: Outside
    integer :: P_Coeffs
    real(rp) :: Phi_Offset
    real(rp) :: REQ_S     ! Req - H_Surf
    integer :: Stat

    my_h_tol = defaultTol ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    hTan_r = H_Tan + Req
    p_coeffs = size(p_basis)
    req_s = req - h_surf
    n_new = 0
    h2 = 0.0 ! Just so it's defined; this value is never used
    do i = tan_ind, 1, -1
      none = .true. ! Assume there will be no intersections
      do j = 1, size(h_ref,2)
        if ( h_tan < 0.0 .and. phi_t > p_basis(j) ) then ! Earth-intersecting ray
          phi_offset = phi_t - 2.0_rp*Acos((req+h_tan)/req)
        else
          phi_offset = phi_t
        end if
        h1 = h2
        h2 = hTan_r / cos(p_basis(j) - phi_offset) - req_s
        if ( j == 1 ) cycle ! It takes two to tango
        if ( (h1-h_ref(i,j-1)) * (h2-h_ref(i,j)) < 0.0 ) then
          ! Line of sight intersects constant-zeta surface.  Solve for where.
          n_new = n_new + 1
          a = (p_basis(j)-p_basis(j-1)) * hTan_r
          b = -(h_ref(i,j)-h_ref(i,j-1))
          c = -(h_ref(i,j-1)-h_surf)*(p_basis(j  )-phi_t) &
            & +(h_ref(i,j  )-h_surf)*(p_basis(j-1)-phi_t) &
            & +(p_basis(j)-p_basis(j-1)) * h_tan
          stat = no_sol
          call Solve_H_Phi ( p_basis(j-1:j), phi_offset, sign(1.0_rp,phi_t-p_basis(j-1)), &
          &                h_ref(i,j-1:j), a, b, c, &
          &                htan_r, req_s, my_h_tol, i /= tan_ind .or. j /= p_coeffs, &
          &                h_new(n_new), p_new(n_new), stat, outside )
          if ( stat == grid ) then
            if ( minval(abs(h_new(n_new)- h_ref(i,j-1:j))) > my_h_tol ) then
              stat = no_sol
              n_new = n_new - 1
              cycle
            end if
          end if
          if ( stat == no_sol .or. stat == complex .or. outside ) then
            ! What's going on?
            n_new = n_new - 1
          else
            if ( minval(abs(p_new(n_new)-p_grid)) < 1.0e-4 ) then
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
