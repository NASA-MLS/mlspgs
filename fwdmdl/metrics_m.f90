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
  public :: Metrics

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ----------------------------------------------------  Metrics  -----

  subroutine Metrics ( &
          ! Input:
          &  phi_t, tan_ind, p_basis, z_ref, n_ref, h_ref, t_ref,      &
          &  dhidzij, csq, refract,                                    &
          ! Output:
          &  h_grid, p_grid, t_grid, dhitdzi, req, status,             &
          ! Optional inputs:
          &  ddhidhidtl0, dhidtlm, tan_press, surf_temp, t_deriv_flag, &
          &  z_basis, h_tol,                                           &
          ! Optional outputs:
          &  ddhtdhtdtl0, dhitdtlm, dhtdtl0, dhtdzt,                   &
          &  do_calc_hyd, do_calc_t, eta_zxp, tan_phi_h, tan_phi_t )

    ! The goal of this subroutine is to return h_grid, p_grid
    ! that define a 2 d integration path

    use Dump_0, only: Dump
    use Geometry, only: EarthRadA
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: RP, IP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Output_m, only: OUTPUT
    use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction
    use Toggles, only: Switches

    ! inputs:

    real(rp), intent(in) :: phi_t      ! orbit projected tangent geodetic angle
    integer(ip), intent(in) :: tan_ind ! tangent height index, 1 = center of
    !                                     longest path
    real(rp), intent(in) :: p_basis(:) ! horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_ref(:)   ! -log pressures (zetas) for which
    !                                     heights/temps are needed.  Only the
    !                                     parts from the tangent outward are used
    real(rp), intent(in) :: n_ref(:,:) ! Indices of refraction by t_phi_basis
    real(rp), intent(in) :: h_ref(:,:) ! heights by t_phi_basis
    real(rp), intent(in) :: t_ref(:,:) ! temperatures by t_phi_basis
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative by t_phi_basis
    real(rp), intent(in) :: csq        ! (minor axis of orbit plane projected
    !                                    Earth ellipse)**2
    logical,  intent(in) :: refract    ! compute phi refractive correction
    ! outputs:
    real(rp), intent(out) :: h_grid(:) ! computed heights
    real(rp), intent(out) :: p_grid(:) ! computed phi's
    real(rp), intent(out) :: t_grid(:) ! computed temperatures
    real(rp), intent(out) :: dhitdzi(:)! derivative of height wrt zeta
    !                                    --may be useful in future computations
    real(rp), intent(out) :: req       ! equivalent elliptical earth radius
    integer, intent(out) :: Status     ! 0 => No trouble, 1 => Convergence failed,
    !                                    2 => Resorted to 1d

    ! optional inputs
    real(rp), optional, intent(in) :: ddhidhidtl0(:,:,:) ! second order reference
    !   temperature derivatives. This is (height, phi_basis, zeta_basis).
    !   Needed only if present(dhidtlm).
    real(rp), optional, intent(inout) :: dhidtlm(:,:,:) ! reference temperature
    !   derivatives. This gets adjusted so that at ref_h(1,@tan phi)) is 0.0 for
    !   all temperature coefficients. This is height X zeta_basis X phi_basis
    real(rp), optional, intent(in) :: Tan_press  ! Tangent pressure
    real(rp), optional, intent(in) :: Surf_temp(:)    ! Surface temperature at phi_basis
    logical, optional, intent(in) :: t_deriv_flag(:)  ! User's deriv. flags for
    !   Temperature. needed only if present(dhidtlm).
    real(rp), optional, intent(in) :: z_basis(:) ! vertical temperature basis
    !   Needed only if present(dhidtlm).
    real(rp), optional, intent(in) :: H_Tol      ! Height tolerance in kilometers
    !   for convergence of phi/h iteration

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
    real(rp), optional, intent(out) :: tan_phi_h ! height at the tangent
    real(rp), optional, intent(out) :: tan_phi_t ! temperature at the tangent

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
    integer :: N_VERT      ! size(z_ref)
    integer :: P_COEFFS    ! Size(P_basis)

    real(rp) :: A, B, C, D ! Polynomial coefficients used to solve for phi and H
    real(rp) :: DP         ! p_basis(j+1) - p_basis(j)
    real(rp) :: DPJ0       ! p_basis(j  )-phi_t
    real(rp) :: DPJ1       ! p_basis(j+1)-phi_t
    real(rp) :: CP2        ! Cos^2 phi_t
    real(rp) :: H          ! Tentative solution for H
    real(rp) :: H_SURF     ! Height of the reference surface -- interpolated in
                           ! row 1 of H_REF
    real(rp) :: H_T        ! Tangent height -- interpolated in row TAN_IND
                           ! of H_REF, then H_SURF is subtracted
    real(rp) :: H_TAN      ! Either H_T or NEG_H_TAN
    real(rp) :: HTAN_R     ! H_Tan + req
    real(rp) :: My_H_Tol   ! Tolerance in kilometers for height convergence
    real(rp) :: P, Q       ! Tentative solutions for phi
    real(rp) :: P2         ! P**2
    real(rp) :: REQ_S      ! Req - H_Surf
    real(rp) :: SecM1      ! Taylor series for sec(phi)-1
    real(rp) :: SP2        ! Sin^2 phi_t

    integer :: VERT_INDS(2*(size(z_ref)+1-tan_ind)) ! What to use in z_ref
    integer :: Stat(size(vert_inds))
    ! Values for stat:
    integer, parameter :: No_sol = 0 ! No solution
    integer, parameter :: Good = 1   ! Newton converged, tangent point, extrapolated
    integer, parameter :: Grid = 2   ! Close to a grid point

    logical :: NOT_ZERO_P(size(vert_inds),size(p_basis))

    real(rp) :: ETA_P(size(vert_inds),size(p_basis))
    real(rp) :: ETA_T(size(p_basis))
    real(rp) :: N_GRID(size(vert_inds))     ! index of refraction
    real(rp) :: PHI_CORR(size(vert_inds))   ! the refractive correction
    real(rp) :: PHI_OFFSET(size(vert_inds)) ! PHI_T or a function of NEG_H_TAN
    real(rp) :: PHI_SIGN(size(vert_inds))   ! +/- 1.0

    ! Coefficients in expansion of Sec(phi)-1 = c2*phi^2 + c4*phi^4 ...
    real(rp), parameter :: C2 = 0.5_rp, C4 = 5.0_rp/24, C6 = 61.0_rp/720.0
    ! Coefficients in expansion of Sec(phi)*Tan(phi) = d/dPhi(sec(phi)-1)
    real(rp), parameter :: D1 = 2*c2, D3 = 4*c4, D5 = 6*c6 ! ... 2n * c_2n

    ! For debugging output format:
    logical, parameter :: clean = .false.

    ! To control debugging
    logical, parameter :: Debug = .false.
    logical, parameter :: NewtonDetails = .true. .and. debug
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

    status = 0 ! assume it will work
    my_h_tol = 0.001_rp ! kilometers
    if ( present(h_tol) ) my_h_tol = h_tol ! H_Tol is in kilometers
    p_coeffs = size(p_basis)
    n_vert = size(z_ref)

    ! compute the tangent height vertical.
    ! For simplicity, we set the surface reference at the input z_ref(1)
    ! and adjust the req and the h_t relative to this, and adjust h_ref
    ! accordingly
    call get_eta_sparse ( p_basis, phi_t, eta_t, first, last )
    h_surf = dot_product(h_ref(1,first:last),eta_t(first:last))
    h_t = dot_product(h_ref(tan_ind,first:last),eta_t(first:last)) - h_surf
    h_tan = h_t

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
    req_s = 0.001_rp*sqrt((earthrada**4*sp2 + csq**2*cp2) / &
                        & (earthrada**2*cp2 + csq*sp2))
    req = req_s + h_surf

    ! Only use the parts of the reference grids that are germane to the
    ! present path

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (i, i=tan_ind,n_vert) /)
    n_path = size(vert_inds)
    n_tan = n_path / 2

    ! sign of phi vector
    phi_sign = (/ (-1, i=1, n_vert-tan_ind+1), (+1, i=n_vert+tan_ind, 2*n_vert) /)

    ! p_basis and p_grid are phi's in offset radians relative to phi_t, that
    ! is, the phi_t, p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t

    phi_offset(:n_path) = phi_t
    if ( present(tan_press) ) then ! Earth intersecting ray
      ! Compute GP height (in KM.) of tan. press. below surface
      h_tan = dot_product(surf_temp(first:last),eta_t(first:last))*(tan_press-z_ref(1))/14.8
      phi_offset(n_tan+1:) = phi_t-2.0_rp*Acos((req+h_tan)/req)
    end if
    htan_r = h_tan + req

    if ( h_phi_dump > 0 ) then
      call dump ( p_basis, name='p_basis', format='(1pg14.6)', clean=clean )
      call output ( phi_t, before='phi_t = ', format='(1pg14.6)' )
      call output ( h_t, format='(f7.2)', before=', h_t = ' )
      call output ( h_surf, format='(f7.2)', before =', h_surf = ' )
      call output ( req, format='(f7.2)', before=', req = ', advance='yes' )
      call output ( tan_ind, before='tan_ind = ' )
      call output ( csq, before=', csq = ', advance='yes' )
      call dump ( phi_offset, name='phi_offset', clean=clean )
      call dump ( h_ref(tan_ind:n_vert,:), name='h_ref', clean=clean )
!     call dump ( z_ref(tan_ind:n_vert), name='z_ref', clean=clean )
    end if

      if ( debug ) &
        & print '(2a,f9.3,a,f7.3)', '   I J    phi  H(TAYLOR)',&
          & '  HREF(I,J) HREF(I,J+1)  H(TRIG)    H-HREF      N   Diffs'

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
    end if

    ! The forward model wants H_grid referenced to the equivalent Earth
    ! surface instead of its center.
    h_grid(:n_path) = h_grid(:n_path) - req

    ! Interpolate Temperature (T_Ref) and the vertical height derivative
    ! (dhidzij) to the path (T_Grid and dhitdzi).

    call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p, NOT_ZERO = not_zero_p )
    do i = 1, n_path
      t_grid(i) = dot_product(t_ref(vert_inds(i),:), eta_p(i,:))
      ! compute the vertical derivative grid
      dhitdzi(i) = dot_product(dhidzij(vert_inds(i),:), eta_p(i,:))
    end do

    ! now for the optional tangent quantities.
    if ( present(tan_phi_h) ) tan_phi_h = h_tan
    if ( present(tan_phi_t) ) &
      & tan_phi_t = dot_product(t_ref(tan_ind,first:last),eta_t(first:last))
    if ( present(dhtdzt) ) &
      & dhtdzt = dot_product(dhidzij(tan_ind,first:last),eta_t(first:last))

    ! compute tangent temperature derivatives
    if ( present(dhidtlm) ) call Tangent_Temperature_Derivatives

    if ( do_dumps > 0 ) then
      call output ( h_tan, before='h_tan = ' )
      call output ( req, before=', req = ', advance='yes' )
      call dump ( p_grid(:n_path), name='p_grid', format='(1pg14.6)', clean=clean )
      call dump ( h_grid(:n_path), name='h_grid', format='(1pg14.6)', clean=clean )
      call dump ( t_grid(:n_path), name='t_grid', format='(1pg14.6)', clean=clean )
      call dump ( dhitdzi(:n_path), name='dhitdzi', format='(1pg14.6)', clean=clean )
      if ( dump_stop > 0 ) stop
    end if

  contains

    ! --------------------------  Tangent_Temperature_Derivatives  -----
    subroutine Tangent_Temperature_Derivatives

      real(rp), dimension(n_path,size(z_basis)) :: ETA_T2
      logical, dimension(n_path,size(z_basis)) :: NOT_ZERO_T
      integer :: I, J, SV_P, SV_T, SV_Z ! Loop inductors and subscripts
      integer :: P_COEFFS               ! size(p_basis)
      integer :: Z_COEFFS               ! size(z_basis)

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

  end subroutine Metrics

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
