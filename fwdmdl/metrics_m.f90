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
          &  dhidzij, beta, refract,                                  &
          ! Output:
          &  h_grid, p_grid, t_grid, dhitdzi, req, status,            &
          ! Optional inputs:
          &  ddhidhidtl0, dhidtlm, neg_h_tan, t_deriv_flag, z_basis,  &
          ! Optional outputs:
          &  ddhtdhtdtl0, dhitdtlm, dhtdtl0, dhtdzt,                  &
          &  do_calc_hyd, do_calc_t, eta_zxp, tan_phi_h, tan_phi_t )

    ! The goal of this program is to return a matrix of h_grids
    ! and t_grids that define 2 d integration paths

    use Dump_0, only: Dump
    use Geometry, only: EarthRadA, EarthRadB
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSKinds, only: RP, IP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Output_m, only: OUTPUT
    use Phi_Refractive_Correction_m, only: Phi_Refractive_Correction
    use Toggles, only: Emit, Switches, Toggle

    ! inputs:

    real(rp), intent(in) :: phi_t      ! orbit projected tangent geodetic angle
    integer(ip), intent(in) :: tan_ind ! tangent height index
    real(rp), intent(in) :: p_basis(:) ! horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_ref(:)   ! -log pressures (zetas) for which
    !                                     heights/temps are needed.  Only the
    !                                     parts from the tangent outward are used
    real(rp), intent(in) :: n_ref(:,:) ! Indices of refraction by t_phi_basis
    real(rp), intent(in) :: h_ref(:,:) ! heights by t_phi_basis
    real(rp), intent(in) :: t_ref(:,:) ! temperatures by t_phi_basis
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative by t_phi_basis
    real(rp), intent(in) :: beta       ! orbital incline angle (Radians)
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

    ! Keywords:
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
    !
    ! optional outputs
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

    ! NOTES
    ! p_basis and p_grid are phi's in offset radians relative to phi_t, that
    ! is the phi_t p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t
    ! compute primary hydrostatic grid
    ! tangent phi index

    ! Local variables.  CAN WE GET SOME COMMENTS FOR EACH OF THESE?
    integer, save :: Do_Dumps = -1  ! -1 = first time, 0 = no dump, >0 = dump
    integer, save :: Dump_Stop = -1 ! -1 = first time, 0 = no dump, >0 = dump/stop
    integer :: I
    integer :: ITER
    integer :: J
    integer :: N_CVF
    integer :: NO_OF_BAD_FITS
    integer :: N_VERT
    integer :: P_COEFFS

    real(rp) :: CP2      ! Cos^2 phi_t
    real(rp) :: CSQ      ! c^2
    real(rp) :: H_T
    real(rp) :: H_SURF
    real(rp) :: SP2      ! Sin^2 phi_t

    integer, dimension(2*(size(z_ref)+1-tan_ind)) :: CVF_INDS
    integer, dimension(size(cvf_inds)) :: VERT_INDS

    logical, dimension(size(cvf_inds)) :: MASK
    logical, dimension(size(cvf_inds),size(p_basis)) :: NOT_ZERO_P

! TEMPORARY until z_grid becomes an argument
real(rp), dimension(size(h_grid)) :: Z_GRID

    real(rp), dimension(size(cvf_inds)) :: CVF_ANG_OFFSET
    real(rp), dimension(size(cvf_inds)) :: CVF_H_TAN
    real(rp), dimension(size(cvf_inds)) :: CVF_SIGN
    real(rp), dimension(size(cvf_inds)) :: CVF_Z_GRID
    real(rp), dimension(size(p_basis)) :: ETA_T
    real(rp), dimension(size(cvf_inds)) :: N_GRID   ! index of refraction
    real(rp), dimension(size(cvf_inds)) :: OLD_HTS
    real(rp), dimension(size(cvf_inds)) :: PHI_CORR ! the refractive correction

    real(rp), dimension(size(cvf_inds),size(p_basis)), target :: ETA_P
    real(rp), dimension(:), pointer :: ETA_P_1 ! Target is a section of eta_p
real(rp), dimension(:,:), pointer :: ETA_PI

    ! For optional calculations

    ! Begin program

    if ( do_dumps < 0 ) then ! First time only
      dump_stop = index(switches,'metD')
      do_dumps = max(dump_stop,index(switches,'metd'))
    end if

    status = 0 ! assume it will work
    p_coeffs = size(p_basis)
    n_vert = size(z_ref)

    !{ $c^2 = \frac{a^2 b^2}{a^2 \sin^2 \beta + b^2 \cos^2 \beta}$.
    !  This is Equation (5.3) in the 19 August 2004 ATBD JPL D-18130.

    csq = (earthrada * earthradb)**2 / &
          & ((earthrada**2-earthradb**2)*sin(beta)**2 + earthradb**2)

    ! compute the tangent height vertical.
    ! For simplicity, we are going to set the
    ! surface reference at the input z_ref(1) and adjust the req and the
    ! h_t relative to this, and adjust h_ref accordingly
    call get_eta_sparse ( p_basis, phi_t, eta_t )
    h_surf = dot_product(h_ref(1,:),eta_t)
    h_t = dot_product(h_ref(tan_ind,:),eta_t) - h_surf

    !{ Compute equivalent earth radius (REQ) at phi\_t(1), nearest surface
    ! and adjust to the input z\_grid(1) using Equation (5.21) in the 19
    ! August 2004 ATBD JPL D-18130.
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

    ! construct some basic arrays
    cvf_h_tan = h_t

    ! determine condensed vector format so we don't mess with a bunch
    ! of zeros

    vert_inds = (/ (i, i=n_vert,tan_ind,-1), (i, i=tan_ind,n_vert) /)

    n_cvf = size(cvf_inds)

    ! cvf_inds = 1:n_vert-tan_ind+1 // n_vert+tan_ind:2*n_vert
    cvf_inds = (/ (i, i=1, n_vert-tan_ind+1), (i, i=n_vert+tan_ind, 2*n_vert) /)

    ! sign of phi vector
    cvf_sign = (/ (-1, i=1, n_vert-tan_ind+1), (+1, i=n_vert+tan_ind, 2*n_vert) /)

    cvf_ang_offset = phi_t
    ! compute phi_s - phi_t
    if ( present(neg_h_tan) ) then
      cvf_ang_offset(n_vert+1:2*n_vert) = phi_t-2.0_rp*Acos((req+neg_h_tan)/req)
      cvf_h_tan = neg_h_tan
    end if

    cvf_z_grid = z_ref(vert_inds)

    h_grid = h_t
 
    !{Estimate the phi's using the tangent heights as a first guess.
    !
    ! $\phi = \phi_t \pm \cos^{-1} \left( \frac{h_t + H^\oplus_t}{h + H^\oplus_t} \right)
    !       = \phi_t \pm \cos^{-1} \left( \frac{h_t + R^\oplus_{\text{eq}}}
    !                                          {h + R^\oplus_{\text{eq}}} \right)$
    !
    ! This is Equation (5.24) in the 19 August 2004 ATBD JPL D-18130.

    p_grid = cvf_ang_offset + cvf_sign * Acos((req+cvf_h_tan)/(req+h_grid))

    ! These are the tangent indicies on a cvf array
    if ( .not. present(neg_h_tan) ) then
      ! some tangents above the earth surface
      ! force the del phi at the tangent to be zero
      p_grid(n_vert - tan_ind + 1:n_vert - tan_ind + 2) = &
      & cvf_ang_offset(n_vert - tan_ind + 1:n_vert - tan_ind + 2)
    end if

    ! Estimate some new heights.
    ! This is going to have some small surface errors away
    ! from the tangent but we will ignore this issue for now.
    ! The following lines replace an extremely sparse matrix
    ! multiplication with a less memory intensive
    ! and probably faster alternative.

    ! now these parts will change with iteration therefore these steps
    ! get repeated
nullify ( eta_pi )

    !{ Interpolate $h$ onto the $\phi$ grid, then compute
    !  $\phi = \phi_t \pm \cos^{-1}
    !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right) =
    !   \phi_t \pm \cos^{-1}
    !   \left( \frac{h_t + R^\oplus_{\text{eq}}}{h+R^\oplus_{\text{eq}}} \right)$.
    !  Iterate until $h$ doesn't change too much.  This is Equation (5.24)
    !  in the 19 August 2004 ATBD JPL D-18130.

    iter = 0
    do

      old_hts = h_grid
      call get_eta_sparse ( p_basis, p_grid, eta_p )
      do i = 1, n_cvf
        h_grid(i) = max ( cvf_h_tan(i), &
          &               dot_product(h_ref(vert_inds(i),:)-h_surf, eta_p(i,:)) )
        if ( refract ) &
          & n_grid(i) = dot_product(n_ref(vert_inds(i),:), eta_p(i,:))
      end do

      iter = iter + 1
      mask = abs(old_hts-h_grid) > 0.01_rp
      no_of_bad_fits = count(mask)
      if ( no_of_bad_fits == 0 ) exit
      if ( iter == 20 ) exit

      ! recompute subsurface phi_s - phi_t because we need it for p_grid,
      ! which we need for h_grid, upon which we are iterating.
      if ( present(neg_h_tan) ) &
        & cvf_ang_offset(n_vert+1:2*n_vert) = phi_t - &
                  & 2.0_rp*Acos((req+neg_h_tan)/(req+h_grid(j)))

      if ( refract ) then
        call phi_refractive_correction ( n_grid, req+h_grid, phi_corr )

        !{ Apply the refractive correction
        !
        !  $\phi = \phi_t \pm \left[ \phi_{\text{corr}} + \cos^{-1}
        !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right)\right] =
        !   \phi_t \pm \left[ \phi_{\text{corr}} +  \cos^{-1}
        !   \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !               {h+R^\oplus_{\text{eq}}} \right)\right]$.
        p_grid = cvf_ang_offset + cvf_sign * &
          & ( phi_corr + Acos((req+cvf_h_tan)/(req+h_grid)) )
      else ! or not....
        p_grid = cvf_ang_offset + cvf_sign * Acos((req+cvf_h_tan)/(req+h_grid))
      end if

    end do

    eta_p_1 => eta_p(1,:)

    if ( no_of_bad_fits > 0 ) call bad_fits

    !  Compute final set of angles (Equation (5.24) again) now that we have
    !  final heights

    if ( refract ) then
        !{ Apply the refractive correction
        !
        !  $\phi = \phi_t \pm \left[ \phi_{\text{corr}} + \cos^{-1}
        !   \left( \frac{h_t + H^\oplus_t}{h+H^\oplus_t} \right)\right] =
        !   \phi_t \pm \left[ \phi_{\text{corr}} +  \cos^{-1}
        !   \left( \frac{h_t + R^\oplus_{\text{eq}}}
        !               {h+R^\oplus_{\text{eq}}} \right)\right]$.
      p_grid = cvf_ang_offset + cvf_sign * &
        & ( phi_corr + Acos((req+cvf_h_tan)/(req+h_grid)) )
    else ! or not...
      p_grid = cvf_ang_offset + cvf_sign * Acos((req+cvf_h_tan)/(req+h_grid))
    end if
    if ( .not. present(neg_h_tan) ) then
      p_grid(n_vert - tan_ind + 1:n_vert - tan_ind + 2) = &
        & cvf_ang_offset(n_vert - tan_ind + 1:n_vert - tan_ind + 2)
    end if

    call get_eta_sparse ( p_basis, p_grid, eta_p, NOT_ZERO = not_zero_p )
    do i = 1, n_cvf
      t_grid(i) = dot_product(t_ref(vert_inds(i),:), eta_p(i,:))
      ! compute the vertical derivative grid
      dhitdzi(i) = dot_product(dhidzij(vert_inds(i),:), eta_p(i,:))
    end do

    !{ Compute $\zeta$ along the path $= \zeta_m$ ({\tt z\_grid})
    !  $= \frac{\Gamma \zeta_0 - T_0 + \sqrt{T_0^2 + 2 \Gamma Z
    !   (H_m-H_0)}}{\Gamma}$ from the
    !  reference $\zeta = \zeta_l$ ({\tt z\_ref}) using $\phi$ along the
    !  path $= \phi_m$ ({\tt p\_basis}) and $H$ along the path $= H_m =
    !  \frac{H_t}{\cos \phi_m}$, where $\Gamma = \frac{T_1-T_0}{\zeta_m-\zeta_0}$
    !  and $Z = \frac{T_1+T_0)(\zeta_1-\zeta_0)}{2(H_1-H_0)}$ and the 0 and 1
    !  subscripts refer to positions on the $H_l$ grid that bracket $H_m$.

    ! now for the optional tangent quantities.
    if ( present(tan_phi_h) ) then
      if ( present(neg_h_tan) ) then
        tan_phi_h = neg_h_tan - h_surf
      else
        tan_phi_h = h_t
      end if
    end if
    if ( present(tan_phi_t) ) tan_phi_t = dot_product(t_ref(tan_ind,:),eta_t)
    if ( present(dhtdzt) ) dhtdzt = dot_product(dhidzij(tan_ind,:),eta_t)

    ! compute tangent temperature derivatives
    if ( present(dhidtlm) ) call Tangent_Temperature_Derivatives

    if ( do_dumps > 0 ) then
      call dump ( h_grid, name='h_grid', clean=.true. )
      call dump ( p_grid, name='p_grid', clean=.true. )
      call dump ( t_grid, name='t_grid', clean=.true. )
      call dump ( dhitdzi, name='dhitdzi', clean=.true. )
      call output ( req, before='req \1 ', advance='yes' )
      if ( dump_stop > 0 ) stop
    end if

  contains
    ! -------------------------------------------------  Bad_Fits  -----
    subroutine Bad_Fits

      integer, dimension(no_of_bad_fits) :: JUNK
      integer :: END_IND
      integer :: END_IND1
      integer :: HI_PT
      integer :: LOW_PT
      integer :: Path_Ind
      integer :: ST_IND
      integer :: TAN_PTR

      if ( toggle(emit) ) then
        call MLSMessage ( MLSMSG_Warning, moduleName, &
          & 'Full convergence not achieved, implementing an improved approximation patch' )
        status = 1
        call output ( 'pth ind, error', advance='yes' )
      end if

      ! We are going to assume that the tangent value is good.
      ! The following is an F90 specific design that is quite different
      ! from the IDL code

      ! Find ranges of contiguous indicies
      junk = PACK((/(i,i=1,n_cvf)/),mask)
      st_ind = 1
      path_ind = cvf_inds(junk(st_ind))
      do end_ind = 1, no_of_bad_fits
        if ( end_ind < no_of_bad_fits ) then
          if ( junk(end_ind+1) - junk(end_ind) < 2 ) cycle
        end if

        if ( toggle(emit) ) then
          call output ( path_ind, places=7 )
          call output ( ', ' )
          call output ( old_hts(junk(st_ind))-h_grid(junk(st_ind)), advance='yes' )
        end if

        ! find which side of the tangent we are on
        if ( path_ind < n_vert + 1 ) then ! near observer side
          low_pt = junk(end_ind) + 1
          hi_pt  = MAX(junk(st_ind) - 1,1)
        else                              ! far observer side
          low_pt = junk(st_ind) - 1
          hi_pt  = MIN(junk(end_ind) + 1,n_cvf)
        end if
        tan_ptr = 2 * n_vert
        if ( cvf_inds(junk(st_ind)) ==  1 ) then
          call mlsMessage(mlsmsg_warning,moduleName,'resorting to 1d option in metrics')
          status = 2
          ! resort to 1 d equvalent
          call get_eta_sparse ( p_basis, phi_t, eta_p_1 )
          do i = st_ind, end_ind
            h_grid(junk(i)) = dot_product(h_ref(n_vert-cvf_inds(junk(i))-1,:),eta_p_1) - h_surf
          end do
        else if ( cvf_inds(junk(end_ind)) > tan_ptr - 1 ) then
          ! calculate the path ending index.
          call mlsMessage ( mlsmsg_warning, moduleName, 'resorting to 1d option in metrics' )
          status = 2
          end_ind1 = st_ind + tan_ptr - cvf_inds(junk(st_ind))
          ! resort to 1 d equivalent
          call get_eta_sparse ( p_basis, phi_t, eta_p_1 )
          do i = st_ind, end_ind1
            h_grid(junk(i)) = dot_product(h_ref(cvf_inds(junk(i))-n_vert,:),eta_p_1) - h_surf
          end do
        else

          ! Correct
          h_grid(junk(st_ind):junk(end_ind)) = h_grid(low_pt) + &
             & (h_grid(hi_pt) - h_grid(low_pt)) * &
             & (cvf_z_grid(junk(st_ind):junk(end_ind))-cvf_z_grid(low_pt)) / &
             & (cvf_z_grid(hi_pt) - cvf_z_grid(low_pt))
        end if
        st_ind = end_ind + 1
        if ( st_ind < no_of_bad_fits ) &
          & path_ind = cvf_inds(junk(st_ind))
      end do

    end subroutine Bad_Fits

    ! --------------------------  Tangent_Temperature_Derivatives  -----
    subroutine Tangent_Temperature_Derivatives

      real(rp), dimension(n_cvf,size(z_basis)) :: ETA_T2
      logical, dimension(n_cvf,size(z_basis)) :: NOT_ZERO_T
      integer :: SV_P, SV_T, SV_Z ! Loop inductors and subscripts
      integer :: Z_COEFFS         ! size(z_basis)

      ! adjust the 2d hydrostatic relative to the surface
      z_coeffs = size(z_basis)
      do i = 1, z_coeffs
        dhidtlm(1:n_vert,i,1:p_coeffs) = dhidtlm(1:n_vert,i,1:p_coeffs) - &
          & dot_product(dhidtlm(1,i,:), eta_t)
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
              eta_zxp(:,sv_t) = eta_t2(:,sv_z) * eta_p(:,sv_p)
            elsewhere
              eta_zxp(:,sv_t) = 0.0
            end where
            do_calc_hyd(:,sv_t) = not_zero_p(:,sv_p) .and. dhidtlm(vert_inds(:),sv_z,sv_p) > 0.0_rp
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
