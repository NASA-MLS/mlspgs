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
          &  phi_t, tan_ind, p_basis, z_ref, n_ref, h_ref, t_ref,     &
          &  dhidzij, csq, refract,                                   &
          ! Output:
          &  h_grid, p_grid, t_grid, dhitdzi, req, status,            &
          ! Optional inputs:
          &  ddhidhidtl0, dhidtlm, neg_h_tan, t_deriv_flag, z_basis,  &
          &  h_tol,                                                   &
          ! Optional outputs:
          &  ddhtdhtdtl0, dhitdtlm, dhtdtl0, dhtdzt,                  &
          &  do_calc_hyd, do_calc_t, eta_zxp, tan_phi_h, tan_phi_t,   &
          &  z_grid, n_path_out )

    ! The goal of this program is to return a matrix of h_grids
    ! and t_grids that define 2 d integration paths

    use Dump_0, only: Dump
    use Geometry, only: EarthRadA
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: RP, IP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Output_m, only: OUTPUT
    use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction
    use Toggles, only: Emit, Switches, Toggle

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
    ! t_grid is not computed if present(z_grid) .and. present(n_path_out)
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
    real(rp), optional, intent(in) :: neg_h_tan  ! sub earth tangent height
    logical, optional, intent(in) :: t_deriv_flag(:)  ! User's deriv. flags for
    !   Temperature. needed only if present(dhidtlm).
    real(rp), optional, intent(in) :: z_basis(:) ! vertical temperature basis
    !   Needed only if present(dhidtlm).
    real(rp), optional, intent(in) :: H_Tol      ! Height tolerance in kilometers
    !   for convergence of phi/h iteration
    !
    ! optional outputs.
    ! Most are not computed if present(z_grid) .and. present(n_path_out)
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
    real(rp), optional, intent(out) :: z_grid(:) ! If present, augment p_grid
    !           with phi's from p_basis (excluding really close values),
    !           compute corresponding h_grid, t_grid, z_grid and dhitdzi, and
    !           report the revised length in N_Path_out.
    integer, optional, intent(out) :: N_Path_out ! Goes with z_grid, q.v.

    ! Local variables.  CAN WE GET SOME COMMENTS FOR EACH OF THESE?
    integer, save :: Do_Dumps = -1  ! -1 = first time, 0 = no dump, >0 = dump
    integer, save :: Dump_Stop = -1 ! -1 = first time, 0 = no dump, >0 = dump/stop
    integer, save :: H_Phi_Dump = -1 ! -1 = first time, 0 = no dump, >0 = dump
    integer, save :: H_Phi_Stop = -1 ! -1 = first time, 0 = no dump, >0 = dump/stop
    integer :: I, I1(1)  ! I1 is for output from MAXLOC
    integer :: ITER      ! Counter for H/Phi iteration
    integer :: Min_Zeta  ! Position where minimum zeta added to z_grid
    integer :: My_Tan    ! Tangent index might change if we add phi's
    integer :: N_PATH    ! Path length = 2*(size(z_ref)+1-tan_ind)
    integer :: NO_OF_BAD_FITS
    integer :: N_VERT    ! size(z_ref)
    integer :: P_COEFFS  ! size(p_basis)

    logical :: Converge  ! H/Phi iteration has converged

    real(rp) :: CP2      ! Cos^2 phi_t
    real(rp) :: H_SURF   ! Height of the reference surface
    real(rp) :: H_T      ! Tangent height
    real(rp) :: H_TAN    ! Either H_T or NEG_H_TAN
    real(rp) :: My_H_Tol
    real(rp) :: Phi_Min  ! Phi where minimum zeta occurs
    real(rp) :: SP2      ! Sin^2 phi_t
    real(rp) :: Zeta_Min ! Minimum zeta

    integer :: VERT_INDS(2*(size(z_ref)+1-tan_ind)) ! What to use in z_ref

    logical :: MASK(size(vert_inds))
    logical :: NOT_ZERO_P(size(vert_inds),size(p_basis))

    real(rp) :: CVF_Z_GRID(size(vert_inds)) ! z_ref(vert_inds)
    real(rp) :: ETA_T(size(p_basis))
    real(rp) :: H_ZMIN                      ! Height at min zeta
    real(rp) :: H_NEW(2) ! At zeta max and another tangent pressure pt
    real(rp) :: N_GRID(size(p_grid))        ! index of refraction
    real(rp) :: OLD_HTS1(size(vert_inds))   ! Previous H in H/Phi iteration,
                                            ! from interpolation
    real(rp) :: OLD_HTS2(2)  ! Previous H in H/Phi iteration, from ACOS(phi)
    real(rp) :: PHI_CORR(size(p_grid))      ! the refractive correction
    real(rp) :: PHI_OFFSET(size(p_grid))    ! PHI_T or a function of NEG_H_TAN
    real(rp) :: PHI_SIGN(size(p_grid))      ! +/- 1.0

    real(rp), target :: ETA_P(size(p_grid),size(p_basis))

    ! Begin program

    if ( do_dumps < 0 ) then ! First time only
      dump_stop = index(switches,'metD')
      do_dumps = max(dump_stop,index(switches,'metd'))
      h_phi_stop = index(switches,'hphI')
      h_phi_dump = max(h_phi_stop,index(switches,'hphi'))
    end if

    my_h_tol = 0.01_rp
    if ( present(h_tol) ) my_h_tol = h_tol

    status = 0 ! assume it will work
    p_coeffs = size(p_basis)
    n_vert = size(z_ref)

    ! compute the tangent height vertical.
    ! For simplicity, we set the surface reference at the input z_ref(1)
    ! and adjust the req and the h_t relative to this, and adjust h_ref
    ! accordingly
    call get_eta_sparse ( p_basis, phi_t, eta_t )
    h_surf = dot_product(h_ref(1,:),eta_t)
    h_t = dot_product(h_ref(tan_ind,:),eta_t) - h_surf
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
    req = 0.001_rp*sqrt((earthrada**4*sp2 + csq**2*cp2) / &
                      & (earthrada**2*cp2 + csq*sp2)) &
        & + h_surf

    ! Only use the parts of the reference grids that are germane to the
    ! present path

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (i, i=tan_ind,n_vert) /)
    n_path = size(vert_inds)
    if ( present(z_grid) ) z_grid(:n_path) = z_ref(vert_inds)

    ! sign of phi vector
    phi_sign(:n_path) = (/ (-1, i=1, n_vert-tan_ind+1), (+1, i=n_vert+tan_ind, 2*n_vert) /)

    ! p_basis and p_grid are phi's in offset radians relative to phi_t, that
    ! is, the phi_t, p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t

    phi_offset(:n_path) = phi_t
    if ( present(neg_h_tan) ) then
      phi_offset(n_vert+1:2*n_vert) = phi_t-2.0_rp*Acos((req+neg_h_tan)/req)
      h_tan = neg_h_tan
    end if

    cvf_z_grid = z_ref(vert_inds)

    ! We will subtract req from h_grid at the end of everything.
    ! We add it in now to avoid forming req+h_grid(:n_path), which would
    ! require an array temp.
    h_grid(:n_path) = h_t + req
    h_new = h_grid(1)

    !{Estimate the phi's using the tangent height as a first guess.

    p_grid(:n_path) = phi_offset(:n_path)

    if ( .not. present(neg_h_tan) ) then
      ! tangent is above the earth surface
      ! force the del phi at the tangent to be zero
      p_grid(n_vert - tan_ind + 1) = phi_offset(n_vert - tan_ind + 1)
      p_grid(n_vert - tan_ind + 2) = phi_offset(n_vert - tan_ind + 2)
    end if

    ! Estimate some new heights.
    ! This is going to have some small surface errors away
    ! from the tangent but we ignore this issue for now.
    ! The following lines replace an extremely sparse matrix
    ! multiplication with a less memory intensive
    ! and probably faster alternative.

    ! now these parts will change with iteration therefore these steps
    ! get repeated

    !{ Interpolate $h$ onto the $\phi$ grid, then compute
    !  $\phi = \phi_t \pm \cos^{-1}
    !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right) =
    !   \phi_t \pm \cos^{-1}
    !   \left( \frac{h_t + R^\oplus_{\text{eq}}}{h+R^\oplus_{\text{eq}}} \right)$.
    !  Iterate until $h$ doesn't change too much.  This is Equation (5.24)
    !  in the 19 August 2004 ATBD JPL D-18130.

    if ( h_phi_dump > 0 ) call dump ( p_basis, name='p_basis', clean=.true. )
    iter = 0
    min_zeta = 0
    old_hts2 = h_t
    call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p(:n_path,:) )
    do

      old_hts1(:n_path) = h_grid(:n_path)
      do i = 1, n_path
        h_grid(i) = max ( h_tan, &
          &               dot_product(h_ref(vert_inds(i),:), eta_p(i,:)) - h_surf ) &
          &               + req ! We will subtract req at the end of everything
        if ( refract ) &
          & n_grid(i) = dot_product(n_ref(vert_inds(i),:), eta_p(i,:))
      end do

      call Assure_H_is_Monotone

      mask = abs(old_hts1(:n_path)-h_grid(:n_path)) > my_h_tol
      no_of_bad_fits = count(mask)
      converge = no_of_bad_fits == 0
      iter = iter + 1
      if ( iter == 20 ) exit

      ! recompute subsurface phi_s - phi_t because we need it for p_grid,
      ! which we need for h_grid, upon which we are iterating.
      if ( present(neg_h_tan) ) &
        & phi_offset(n_vert+1:2*n_vert) = phi_t - &
                  & 2.0_rp*Acos((req+neg_h_tan)/h_grid(n_vert+1:2*n_vert))

      if ( refract ) then
        call phi_refractive_correction ( n_grid(:n_path), h_grid(:n_path), phi_corr(:n_path) )

        !{ Apply the refractive correction
        !
        !  $\phi = \phi_t \pm \left[ \phi_{\text{corr}} + \cos^{-1}
        !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right)\right] =
        !   \phi_t \pm \left[ \phi_{\text{corr}} +  \cos^{-1}
        !   \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !               {h+R^\oplus_{\text{eq}}} \right)\right]$.
        p_grid(:n_path) = phi_offset(:n_path) + &
          & phi_sign(:n_path) * ( phi_corr(:n_path) + Acos((req+h_tan)/h_grid(:n_path)) )
      else
        !{ or not\dots. $\phi = \phi_t \pm \cos^{-1} \left( \frac{h_t + H^\oplus_t}{h + H^\oplus_t} \right)
        !       = \phi_t \pm \cos^{-1} \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !                                          {h + R^\oplus_{\text{eq}}} \right)$.
        !
        ! This is Equation (5.24) in the 19 August 2004 ATBD JPL D-18130.
        p_grid(:n_path) = phi_offset(:n_path) + &
          & phi_sign * Acos((req+h_tan)/h_grid(:n_path))
      end if

      if ( h_phi_dump > 0 ) then
        call dump ( (/ req+h_tan /), name='req+h_tan' )
        call dump ( h_grid(:n_path), name='h_grid', clean=.true. )
        call dump ( p_grid(:n_path), name='p_grid', clean=.true. )
      end if

      call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p(:n_path,:) )
      if ( present(z_grid) ) then
        ! Find the maximum zeta by Hermite interpolation
        do i = 1, n_path
          ! compute the height vertical derivative w.r.t. zeta grid
          dhitdzi(i) = dot_product(dhidzij(vert_inds(i),:), eta_p(i,:))
        end do
        call get_min_zeta
        if ( min_zeta > 0 ) then
          if ( abs(old_hts2(1)-h_zmin ) > my_h_tol ) converge = .false.
          old_hts2(1)=h_zmin
        end if
      end if

      if ( converge ) exit

    end do

    if ( no_of_bad_fits > 0 ) call bad_fits ( eta_p(1,:) )

    !  Compute final set of angles (Equation (5.24) again) now that we have
    !  final heights

    if ( h_phi_dump > 0 ) then
      if ( any((req+h_tan)/h_grid(:n_path) > 1.0 ) ) then
        call output ( req+h_tan, before='ACOS((req+h_tan)/h_grid(:n_path)) is going to fail, req+h_tan =' )
        call dump ( h_grid(:n_path), name='h_grid', clean=.true. )
      end if
    end if

    if ( refract ) then
        !{ Apply the refractive correction
        !
        !  $\phi = \phi_t \pm \left[ \phi_{\text{corr}} + \cos^{-1}
        !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right)\right] =
        !   \phi_t \pm \left[ \phi_{\text{corr}} +  \cos^{-1}
        !   \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !               {h+R^\oplus_{\text{eq}}} \right)\right]$.

      p_grid(:n_path) = phi_offset(:n_path) + &
        & phi_sign(:n_path) * ( phi_corr(:n_path) + Acos((req+h_tan)/h_grid(:n_path)) )
    else
        !{ or not\dots. $\phi = \phi_t \pm \cos^{-1} \left( \frac{h_t + H^\oplus_t}{h + H^\oplus_t} \right)
        !       = \phi_t \pm \cos^{-1} \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !                                          {h + R^\oplus_{\text{eq}}} \right)$.
        !
        ! This is Equation (5.24) in the 19 August 2004 ATBD JPL D-18130.

      p_grid(:n_path) = phi_offset(:n_path) + &
        & phi_sign * Acos((req+h_tan)/h_grid(:n_path))
    end if
    if ( .not. present(neg_h_tan) ) then
      p_grid(n_vert - tan_ind + 1) = phi_offset(n_vert - tan_ind + 1)
      p_grid(n_vert - tan_ind + 2) = phi_offset(n_vert - tan_ind + 2)
    end if

    h_grid(:n_path) = h_grid(:n_path) - req

    ! Must we add p_basis to p_grid?
    if ( present(z_grid) .and. present(n_path_out) ) then
      my_tan = n_path / 2
      do i = 1, size(p_basis)
        call augment_grids ( p_basis(i) )
      end do
      ! For now, we are assuming the tangent is at the midpoint of the
      ! path.  If we've added points unsymmetrically, we need to add some
      ! more points to balance things.  When the rest of the code is able
      ! to handle paths for which the tangent isn't at the midpoint, delete
      ! this.
      do while ( my_tan /= n_path / 2 )
        if ( my_tan < n_path / 2 ) then
          ! Need to add a point before the tangent to balance it
          i1 = maxloc(p_grid(2:my_tan)-p_grid(1:my_tan-1))
        else ! Need to add a point after the tangent to balance it
          i1 = maxloc(p_grid(my_tan+1:n_path)-p_grid(my_tan:n_path-1)) + my_tan
        end if
        i = i1(1) ! Because maxloc returns a rank-1 array of extent 1
        call augment_grids ( 0.5 * (p_grid(i)+p_grid(i+1)) )
      end do
      n_path_out = n_path
      return
    end if

    call get_eta_sparse ( p_basis, p_grid(:n_path), eta_p(:n_path,:), NOT_ZERO = not_zero_p )
    do i = 1, n_path
      t_grid(i) = dot_product(t_ref(vert_inds(i),:), eta_p(i,:))
      ! compute the vertical derivative grid
      dhitdzi(i) = dot_product(dhidzij(vert_inds(i),:), eta_p(i,:))
    end do

    ! now for the optional tangent quantities.
    if ( present(tan_phi_h) ) tan_phi_h = h_tan
    if ( present(tan_phi_t) ) tan_phi_t = dot_product(t_ref(tan_ind,:),eta_t)
    if ( present(dhtdzt) ) dhtdzt = dot_product(dhidzij(tan_ind,:),eta_t)

    ! compute tangent temperature derivatives
    if ( present(dhidtlm) ) call Tangent_Temperature_Derivatives

    if ( do_dumps > 0 ) then
      call dump ( h_grid(:n_path), name='h_grid', clean=.true. )
      call dump ( p_grid(:n_path), name='p_grid', clean=.true. )
      call dump ( t_grid(:n_path), name='t_grid', clean=.true. )
      call dump ( dhitdzi(:n_path), name='dhitdzi', clean=.true. )
      call output ( req, before='req \1 ', advance='yes' )
      if ( dump_stop > 0 ) stop
    end if

  contains

    ! ---------------------------------------  Assure_H_is_Monotone  -----
    subroutine Assure_H_is_Monotone
      ! Check whether H_Grid is monotone increasing away from the midpoint.
      ! If not, replace nonincreasing points by interpolation in Zeta
      ! (cvf_z_grid) using nearby ones that are increasing.
      integer :: I, I1, I2, J, M
      real :: R
      i1 = n_path/2-1
      i2 = 1
      do m = -1, 1, 2
        do i = i1, i2, m
          if ( h_grid(i) <= h_grid(i-m) ) then ! not increasing
            do j = i+m, i2, m
              if ( h_grid(j) > h_grid(i) ) exit
            end do
            if ( j /= i2+m ) then ! not at the end
              r = (cvf_z_grid(i)-cvf_z_grid(i-m)) / &
                & (cvf_z_grid(j)-cvf_z_grid(i-m))
              h_grid(i) = (1.0 - r) * h_grid(i-m) + r * h_grid(j)
            else ! at the end
              h_grid(i) = h_grid(i-m) + 1.0 ! use SWAG method
            end if
          end if
        end do
        i1 = n_path/2 + 2
        i2 = n_path
      end do
    end subroutine Assure_H_is_Monotone

    ! ----------------------------------------------  Augment_Grids  -----
    subroutine Augment_Grids ( New_Phi )

      !{ Let $\phi$ be represented by {\tt p\_basis}, $\zeta^\text{ref}$ by
      !  {\tt z\_ref}, $H^\text{ref}$ by {\tt h\_ref} and $\zeta$ by {\tt
      !  z\_grid}. Compute additional $\zeta_m = \zeta^\text{ref}_i - \frac{T_i
      !  - \sqrt{T_i^2 + 2 \Gamma Z (H_m-H^\text{ref}_{i,m})}}{\Gamma}$ at
      !  $\phi_m$ points along the path, where $H_m = \frac{H_t}{\cos \phi_m}$,
      !  $H_t$ is the tangent height, $\Gamma =
      !  \frac{T_{i+1}-T_i}{\zeta^\text{ref}_{i+1}-\zeta^\text{ref}_i}$, $Z =
      !  \frac{(T_{i+1}+T_i)(\zeta^\text{ref}_{i+1}-\zeta^\text{ref}_i)}
      !  {2(H^\text{ref}_{i+1,m}-H^\text{ref}_{i,m})}$, and $i$ is such that
      !  $H^\text{ref}_{i,m} < H_m < H^\text{ref}_{i+1,m}$.

      use MLSKinds, only: RP

      real(rp), intent(in) :: New_Phi

      integer :: J
      real(rp) :: G, R, Z

      ! Determine whether to add a phi.  Don't add it at the end.
      if ( all(abs(new_phi-p_grid(:n_path)) > 0.01) .and. &
        & new_phi > p_grid(1) .and. new_phi < p_grid(n_path) ) then
        n_path = n_path + 1
        ! Find where the new phi goes by doing an insertion sort,
        ! assuming p_grid is already sorted.
        do j = n_path, 2, -1
          if ( new_phi > p_grid(j) ) then
            h_grid(j) = h_grid(j-1)
            p_grid(j) = p_grid(j-1)
            if ( refract ) phi_corr(j) = phi_corr(j-1)
            phi_offset(j) = phi_offset(j-1)
            phi_sign(j) = phi_sign(j-1)
            t_grid(j) = t_grid(j-1)
            z_grid(j) = z_grid(j-1)
            if ( j == my_tan ) my_tan = my_tan + 1
          else
            p_grid(j) = new_phi
            r = (new_phi - p_grid(j-1)) / &
              & (p_grid(j+1) - p_grid(j-1))
            if ( refract ) then
              phi_corr(j) = phi_corr(j-1) * (1.0-r) + phi_corr(j+1) * r
              h_grid(j) = (h_tan + req) / &
                & cos(p_grid(j)-phi_offset(j)-phi_sign(j)*phi_corr(j)) - req
            else ! No refractive correction
              h_grid(j) = h_tan/cos(p_grid(j)-phi_offset(j))
            end if
            g = (t_grid(j+1)-t_grid(j-1))/(z_grid(j+1)-z_grid(j-1))
            z = ((t_grid(j+1)+t_grid(j-1))*(z_grid(j+1)-z_grid(j-1))) / &
              & (2.0_rp * (h_grid(j+1)-h_grid(j-1)))
            z_grid(j) = z_grid(j-1) - (t_grid(j-1) - &
              &            sqrt(t_grid(j-1)**2 + 2.0_rp * G * Z * (h_grid(j+1) - &
              &                 h_grid(j-1)))) / &
              &           g
            ! In case another one lands in the same area:
            t_grid(j) = t_grid(j-1) * (1.0-r) + t_grid(j+1) * r
            exit
          end if
        end do ! j = n_path, 2, -1
      end if
    end subroutine Augment_Grids

    ! -------------------------------------------------  Bad_Fits  -----
    subroutine Bad_Fits ( ETA_P_1 )

      real(rp), dimension(:), intent(out) :: ETA_P_1 ! Scratch, not used by caller
      
      integer, dimension(no_of_bad_fits) :: JUNK
      integer :: END_IND
      integer :: HI_PT
      integer :: LOW_PT
      integer :: ST_IND

      if ( toggle(emit) ) then
        call MLSMessage ( MLSMSG_Warning, moduleName, &
          & 'Full convergence not achieved, implementing an improved approximation patch' )
        status = 1
        call output ( 'pth ind, error', advance='yes' )
      end if

      ! We are going to assume that the tangent value is good.
      ! The following is an F90 specific design that is quite different
      ! from the IDL code

      junk = PACK((/(i,i=1,n_path)/),mask)
      st_ind = 1
      do end_ind = 1, no_of_bad_fits
        ! Find ranges of contiguous indicies
        if ( end_ind < no_of_bad_fits ) then
          if ( junk(end_ind+1) - junk(end_ind) < 2 ) cycle
        end if

        if ( toggle(emit) ) then
          call output ( junk(st_ind), places=7, after=', ' )
          call output ( old_hts1(junk(st_ind))-h_grid(junk(st_ind)), advance='yes' )
        end if

        if ( junk(st_ind) ==  1 .or. junk(end_ind) == n_path ) then
          call mlsMessage ( mlsmsg_warning, moduleName, 'resorting to 1d option in metrics' )
          status = 2
          ! resort to 1 d equivalent
          call get_eta_sparse ( p_basis, phi_t, eta_p_1 )
          do i = st_ind, end_ind
            h_grid(junk(i)) = &
              & dot_product(h_ref(vert_inds(junk(i)),:),eta_p_1) - h_surf + req
          end do
        else
          ! find which side of the tangent we are on
          if ( junk(st_ind) <= n_path/2 ) then ! near observer side
            low_pt = junk(end_ind) + 1
            hi_pt  = MAX(junk(st_ind) - 1, 1)
          else                                 ! far observer side
            low_pt = junk(st_ind) - 1
            hi_pt  = MIN(junk(end_ind) + 1, n_path)
          end if
          ! Correct unconverged H_Grid using linear Zeta interpolation
          h_grid(junk(st_ind):junk(end_ind)) = h_grid(low_pt) + &
             & (h_grid(hi_pt) - h_grid(low_pt)) * &
             & (cvf_z_grid(junk(st_ind):junk(end_ind))-cvf_z_grid(low_pt)) / &
             & (cvf_z_grid(hi_pt) - cvf_z_grid(low_pt))
        end if
        st_ind = end_ind + 1
        if ( h_phi_dump > 0 ) call dump ( h_grid(:n_path), 'h_grid in Bad_Fits', clean=.true. )
      end do

      call Assure_H_is_Monotone
      if ( h_phi_stop > 0 ) stop

    end subroutine Bad_Fits

    ! -----------------------------------------------  Get_Min_Zeta  -----
    subroutine Get_Min_Zeta
    ! Find Zeta_Min = the minimum zeta and Min_Zeta = where in z_grid it
    ! appears.  Min_Zeta == 0 if it's too close to a phi that's already there.
      use Hermite_Coeffs_m, only: Hermite_coeffs
      real(rp) :: Coeffs(0:3,2:n_path) ! coeffs(:,i) for [i-1,i] segment of path
      real(rp) :: D, dz_dphi(n_path)
      integer :: I
      real(rp) :: X_min, Zeta
      !{ First get {\tt dz\_dphi} =
      !  $\frac{\text{d}h}{\text{d}\phi} \frac{\text{d}\zeta}{\text{d}h} =
      !   \sec \phi \tan \phi \frac{\text{d}\zeta}{\text{d}h} =
      !   (H+R^\oplus_{\text{eq}}) \tan \phi \frac{\text{d}\zeta}{\text{d}h}$.
      if ( refract ) then
        dz_dphi = (h_grid(:n_path)+req) * &
          &  tan(p_grid(:n_path)-phi_offset(:n_path)-phi_sign(:n_path)*phi_corr(:n_path)) / &
          &  dhitdzi(:n_path)
      else
        dz_dphi = (h_grid(:n_path)+req) * &
          &  tan(p_grid(:n_path)-phi_offset(:n_path)) / &
          &  dhitdzi(:n_path)
      end if
      !{ Get the coefficients for the Hermite polynomials in each segment of
      ! the path.  Use them to compute zeta min:
      ! $\zeta^\prime = c_1 + c_2 x + c_3 x^2$, so $\zeta^\prime = 0$ when
      ! $x = \frac{-c_2 -\sqrt{c_2^2-3 c_1 c_3}}{3c_3}$ provided $c_3 \neq 0$
      ! or $x = -\frac{c_1}{2 c_2}$ if $c_3 = 0$.  If $x<0$ or $x>1$ or $x$
      ! is complex, the extremum is at a grid point.  See wvs-041r1.
      call hermite_coeffs ( z_grid(:n_path), dz_dphi, coeffs )
      zeta_min = minval(z_grid(:n_path))
      min_zeta = 0               ! Indicate it's already in the grid
      do i = 2, n_path
        if ( coeffs(3,i) == 0.0 ) then
          if ( coeffs(2,i) <= 0.0 ) cycle ! no minimum within the interval
          x_min = -0.5 * coeffs(1,i) / coeffs(2,i)
        else
          d = coeffs(2,i)**2 - 3.0*coeffs(1,i)*coeffs(3,i)
          if ( d < 0 ) cycle ! complex zero
          x_min = (-coeffs(2,i) - sqrt(d))/(3.0*coeffs(3,i))
        end if
        zeta = coeffs(0,i) + x_min * ( coeffs(1,i) + x_min * (coeffs(2,i) + &
          &      x_min * coeffs(3,i)))
        if ( zeta < zeta_min ) then
          if ( min(abs(x_min-0.005),abs(x_min-0.995)) < 0.005 ) then
            ! 0 < x_min < 0.01 or 0.99 < x_min < 1.0
            min_zeta = 0 ! It's already in the grid
          else
            zeta_min = zeta
            phi_min = x_min * p_grid(i) + (1.0-x_min) * p_grid(i-1)
            h_zmin = (h_tan + req) / &
                & cos(phi_min-phi_offset(i-1)-phi_sign(i-1)*phi_corr(i-1)) - req
            min_zeta = i
          end if
        end if
      end do ! i
    end subroutine Get_Min_Zeta

    ! --------------------------  Tangent_Temperature_Derivatives  -----
    subroutine Tangent_Temperature_Derivatives

      real(rp), dimension(n_path,size(z_basis)) :: ETA_T2
      logical, dimension(n_path,size(z_basis)) :: NOT_ZERO_T
      integer :: I, J, SV_P, SV_T, SV_Z ! Loop inductors and subscripts
      integer :: Z_COEFFS               ! size(z_basis)

      ! adjust the 2d hydrostatic relative to the surface
      z_coeffs = size(z_basis)
      do i = 1, z_coeffs
        dhidtlm(:,i,:) = dhidtlm(:,i,:) - dot_product(dhidtlm(1,i,:), eta_t)
      end do

      j = z_coeffs * p_coeffs
      dhtdtl0 = RESHAPE(dhidtlm(tan_ind,:,:) * SPREAD(eta_t,1,z_coeffs),&
                     & (/j/))

      ddhtdhtdtl0 = RESHAPE( &
                   ddhidhidtl0(tan_ind,:,:) * SPREAD(eta_t,1,z_coeffs), &
                     & (/j/))

      ! compute the path temperature noting where the zeros are
      call get_eta_sparse ( z_basis, cvf_z_grid, eta_t2, NOT_ZERO = not_zero_t )

      sv_t = 0
      do sv_p = 1 , p_coeffs
        do sv_z = 1 , z_coeffs
          sv_t = sv_t + 1
          if ( t_deriv_flag(sv_t) ) then
            do_calc_t(:,sv_t) = not_zero_t(:,sv_z) .and. not_zero_p(:,sv_p)
            where ( do_calc_t(:,sv_t) )
              eta_zxp(:,sv_t) = eta_t2(:,sv_z) * eta_p(:n_path,sv_p)
            elsewhere
              eta_zxp(:,sv_t) = 0.0
            end where
            do_calc_hyd(:,sv_t) = not_zero_p(:,sv_p) .and. &
              &                   dhidtlm(vert_inds(:),sv_z,sv_p) > 0.0_rp
            where ( do_calc_hyd(:,sv_t) )
              dhitdtlm(:,sv_t) = dhidtlm(vert_inds(:),sv_z,sv_p) * eta_p(:n_path,sv_p)
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
! Revision 2.29  2006/06/29 19:34:17  vsnyder
! Remove capability to handle more than one tangent -- substantial simplification
!
! Revision 2.28  2006/03/06 20:44:48  vsnyder
! Cannonball polishing
!
! Revision 2.27  2006/01/26 03:06:17  vsnyder
! Don't deallocate 'mask' until after it's no longer needed!
!
! Revision 2.26  2006/01/05 00:03:52  vsnyder
! Implement refractive correction for Phi
!
! Revision 2.25  2005/12/10 03:31:09  vsnyder
! Replace SUM(A*B) with Dot_Product(A,B) to avoid array temps
!
! Revision 2.24  2005/12/10 01:53:23  vsnyder
! Use get_eta_sparse_1d
!
! Revision 2.23  2005/06/22 18:08:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.22  2004/09/01 01:47:56  vsnyder
! Add status argument
!
! Revision 2.21  2004/01/23 18:44:57  bill
! problem with approximate correction patch when interpolation can't be done tentatively implimentaed
!
! Revision 2.20  2003/11/14 21:22:46  livesey
! Bug fix and output tidy up
!
! Revision 2.19  2003/06/27 23:43:34  vsnyder
! Remove unreferenced USE names
!
! Revision 2.18  2003/06/20 23:41:48  vsnyder
! Revise how compound etas are computed -- hopefully a faster method
!
! Revision 2.17  2003/05/22 20:12:32  vsnyder
! Get rid of some array temps
!
! Revision 2.16  2002/10/25 22:24:08  livesey
! Made the warning output optional
!
! Revision 2.15  2002/10/08 17:08:05  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.14  2002/10/08 14:42:57  bill
! fixed bug in non convergent estimator?
!
! Revision 2.13  2002/10/01 02:27:49  vsnyder
! Reduce number of array temps, move allocation for some out of loops.  Fix
! a bug (junk could be accessed after being deallocated).  Cosmetic changes.
!
! Revision 2.12  2002/09/26 20:14:24  vsnyder
! Get PI from Units module
!
! Revision 2.11  2002/09/25 23:35:13  vsnyder
! Simplify equivalent earth radius, insert copyright notice
!
! Revision 2.10  2002/09/06 18:16:41  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.9  2002/08/10 00:13:33  livesey
! Tiny bug fix to Bill's bug fix
!
! Revision 2.8  2002/08/10 00:08:47  bill
! fixed bad pt correction bug
!
! Revision 2.7  2002/07/05 07:52:50  zvi
! Coor. switch (phi,z) -> (z,phi)
!
! Revision 2.6  2002/06/19 11:00:35  zvi
! Some cosmetic corrections
!
! Revision 2.5  2002/06/07 04:50:47  bill
! fixes and improvements--wgr
!
! Revision 2.4  2002/02/08 00:48:09  zvi
! Restoring the t_deriv_flag code
!
! Revision 2.2  2002/01/30 01:11:21  zvi
! Fix bug in user selectable coeff. code
!
! Revision 2.1  2001/11/20 01:19:30  zvi
! Some clarification of code
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.1.2.8  2001/09/13 01:26:53  zvi
! fixing allocation bug
!
! Revision 1.1.2.7  2001/09/13 00:42:18  livesey
! Fixed memory leak
!
! Revision 1.1.2.6  2001/09/12 22:45:54  livesey
! A version that seems to work
!
! Revision 1.1.2.5  2001/09/12 22:04:07  livesey
! Interim still broken version
!
! Revision 1.1.2.4  2001/09/12 21:47:50  livesey
! More bug fixes
!
! Revision 1.1.2.3  2001/09/12 21:36:51  livesey
! Bug fix.
!
! Revision 1.1.2.2  2001/09/12 21:30:04  livesey
! Updates changes to pointers etc.
!
