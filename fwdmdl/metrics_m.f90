! The goal of this program is to return a matrix of h_grids
! and t_grids which define 2 d integration paths

Module metrics_m

  USE MLSCommon, only: RP, IP
  USE Geometry, only: earthrada,earthradb,pi
  USE get_eta_matrix_m, only: get_eta_sparse
  USE Allocate_deallocate, only: Allocate_test, Deallocate_test

  USE dump_0, only: dump

  IMPLICIT NONE
  private

  public :: metrics
  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  Integer, save :: CALLEDTIMES = 0

contains

  ! ------------------------------------------------ Metrics --------------

  subroutine metrics(phi_t,tan_inds,p_basis,z_grid,h_ref,t_ref,dhidzij, &
          &  beta,t_deriv_flag,h_grid,p_grid,t_grid,dhitdzi,req,dhidtlm,&
          &  ddhidhidtl0,dhitdtlm,eta_zxp,tan_phi_h_grid,tan_phi_t_grid,&
          &  dhtdzt,dhtdtl0,ddhtdhtdtl0,neg_h_tan,z_basis,do_calc_t,    &
          &  do_calc_hyd)

    ! The main change in this routine is to interchange the z_coeffs and
    ! p_coeffs so as to make it consistent with the radiative transfer program.

    ! inputs:

    real(rp), intent(in) :: phi_t(:)   !orbit projected tangent geodetic angles
    integer(ip), intent(in) :: tan_inds(:)!tangent height indicies
    real(rp), intent(in) :: p_basis(:) !horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_grid(:)!pressures for which heights/temps are
    !                                 needed
    real(rp), intent(in) :: h_ref(:,:)! heights by t_phi_basis
    real(rp), intent(in) :: t_ref(:,:)! temperatures by t_phi_basis
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative by t_phi_basis
    real(rp), intent(in) :: beta       !orbital incline angle (Radians)
    logical,  intent(in) :: t_deriv_flag(:)  ! User's deriv. flags for Temp.
    ! outputs:
    real(rp), intent(out) :: h_grid(:)!computed heights
    real(rp), intent(out) :: p_grid(:)!computed phi's
    real(rp), intent(out) :: t_grid(:)!computed temperatures
    real(rp), intent(out) :: dhitdzi(:)!derivative of height wrt zeta
    !                              --maybe useful in future computations
    real(rp), intent(out) :: req       !equivalent elliptical earth radius

    ! Keywords:
    ! optional inputs
    real(rp), optional, intent(in) :: z_basis(:) ! vertical temperature basis
    real(rp), optional, intent(inout) :: dhidtlm(:,:,:) ! reference temperature
    !   derivatives. This gets adjusted so that at ref_h(1,@tan phi)) is 0.0 for
    !   all temperature coefficients. This is height X phi_basis X zeta_basis
    real(rp), optional, intent(in) :: ddhidhidtl0(:,:,:)! second order reference
    !   temperature derivatives. This is (height, phi_basis, zeta_basis)
    real(rp), optional, intent(in) :: neg_h_tan(:)  !sub earth tangent heights
    ! optional outputs
    real(rp), optional, intent(out) :: tan_phi_h_grid(:)!heights along the
    !                                                      tangent
    real(rp), optional, intent(out) :: tan_phi_t_grid(:)!temperature along the
    !                                                      tangent
    real(rp), optional, intent(out) :: dhtdzt(:)     !height derivative wrt
    !                                                   pressure along the tangent
    real(rp), optional, intent(out) :: dhtdtl0(:,:)!First order derivative at
    !                                                 the tangent
    real(rp), optional, intent(out) :: ddhtdhtdtl0(:,:)!second order derivative
    !                            at the tangent only---used for antenna affects
    real(rp), optional, intent(out) :: dhitdtlm(:,:)
    !                             derivative of path position wrt temperature
    !                             statevector (z_basisXphi_basis)
    real(rp), optional, intent(out) :: eta_zxp(:,:) ! eta matrix for temperature.
    logical, optional, intent(out) :: do_calc_t(:,:) ! non zero locater for
    !                                       temperature bases computations.
    logical, optional, intent(out) :: do_calc_hyd(:,:) !non zero locator for
    !                                       hydrostatic calculations.

    ! NOTES
    ! The phi basis is wholly independent of phi_t
    ! compute primary hydrostatic grid
    ! tangent phi index

    ! ===>> EVERYTHING comes out in cvf format <<===!

    ! Local variables.  CAN WE GET SOME COMMENTS FOR EACH OF THESE?
    integer :: P_COEFFS
    integer :: Z_COEFFS
    integer :: N_VERT
    integer :: N_TAN
    integer :: SS_HTAN
    integer :: SV_T
    integer :: SV_P
    integer :: SV_Z
    integer :: I
    integer :: J
    integer :: N_CVF
    integer :: ITER
    integer :: NO_OF_BAD_FITS
    integer :: ST_IND
    integer :: END_IND
    integer :: LOW_PT
    integer :: HI_PT

    real(rp) :: C
    real(rp) :: CP2
    real(rp) :: SP2
    real(rp) :: H_SURF

    logical, dimension(:), pointer :: MASK
    logical, dimension(:,:), pointer :: NOT_ZERO_P
    logical, dimension(:,:), pointer :: NOT_ZERO_T

    integer, dimension(:), pointer :: CVF_INDS
    integer, dimension(:), pointer :: PATH_INDS
    integer, dimension(:), pointer :: VERT_INDS
    integer, dimension(:), pointer :: NEAR_INDS
    integer, dimension(:), pointer :: PATH_IND
    integer, dimension(:), pointer :: JUNK
    integer, dimension(:), pointer :: FORCE_ZERO
    integer, dimension(:), pointer :: INDS
    integer, dimension(:), pointer :: TAN_IND

    real(rp), dimension(:), pointer :: ARG
    real(rp), dimension(:), pointer :: CVF_ANG_OFFSET
    real(rp), dimension(:), pointer :: CVF_H_TAN
    real(rp), dimension(:), pointer :: CVF_SIGN
    real(rp), dimension(:), pointer :: CVF_Z_GRID
    real(rp), dimension(:), pointer :: DHTDTL
    real(rp), dimension(:), pointer :: H_BETTER
    real(rp), dimension(:), pointer :: H_GRID_T
    real(rp), dimension(:), pointer :: H_GRID_TT
    real(rp), dimension(:), pointer :: OLD_HTS
    real(rp), dimension(:), pointer :: SOME_PHI

    real(rp), dimension(:,:), pointer :: ETA_P
    real(rp), dimension(:,:), pointer :: ETA_T
    real(rp), dimension(:,:), pointer :: H_TANS
    real(rp), dimension(:,:), pointer :: H_ZF
    real(rp), dimension(:,:), pointer :: M_Z_GRID
    real(rp), dimension(:,:), pointer :: PHI_TANS

    ! Begin program

    calledTimes=calledTimes + 1
    p_coeffs = size(p_basis)
    n_vert = size(z_grid)
    n_tan = size(tan_inds)

    c = earthrada*earthradb &
      / sqrt(earthrada**2*sin(beta)**2 &
      + earthradb**2*cos(beta)**2)

    nullify ( h_grid_t, eta_t, m_z_grid, h_tans, phi_tans )
    call Allocate_test ( h_grid_t, n_tan, 'h_grid_t', ModuleName )
    call Allocate_test ( eta_t, n_tan, p_coeffs, 'eta_t', ModuleName )
    call Allocate_test ( m_z_grid, 2*n_vert, n_tan, 'm_z_grid', ModuleName )
    call Allocate_test ( h_tans, 2*n_vert, n_tan, 'h_tans', ModuleName )
    call Allocate_test ( phi_tans, 2*n_vert, n_tan, 'phi_tans', ModuleName )

    ! compute the tangent height vertical
    call get_eta_sparse(p_basis,phi_t,eta_t)
    h_grid_t = sum(h_ref(tan_inds,:)*eta_t,dim=2)
    h_surf = sum(h_ref(1,:)*eta_t(1,:))

    ! compute equivalent earth radius at phi_t(1), nearest surface
    cp2 = cos(phi_t(1))**2
    sp2 = 1.0_rp - cp2
    req = 0.001_rp*sqrt((earthrada**4*sp2 + c**4*cp2)/(earthrada**2*cp2 &
      + c**2*sp2))

    ! now for the path quantities, for simplicity, we are going to set the
    ! surface reference at the inputted z_grid(1) and adjust the req and the
    ! h_grid_t relative to this
    req = req + h_surf

    ! adjust h_ref accordingly
    h_grid_t = h_grid_t - h_surf
    if(present(tan_phi_h_grid)) tan_phi_h_grid = h_grid_t
    if(present(tan_phi_t_grid)) &
      tan_phi_t_grid = sum(t_ref(tan_inds,:)*eta_t,dim=2)
    if(present(dhtdzt)) dhtdzt = sum(dhidzij(tan_inds,:)*eta_t,dim=2)

    ! compute tangent temperature derivative
    if (present(dhidtlm)) then
      z_coeffs = size(z_basis)
      nullify ( dhtdtl )
      call Allocate_test ( dhtdtl, z_coeffs, 'dhtdtl', ModuleName )
      dhtdtl = sum(RESHAPE(dhidtlm(1,:,:),(/p_coeffs,z_coeffs/)) &
        * SPREAD(RESHAPE(eta_t(1,:),(/p_coeffs/)),2,z_coeffs),dim=1)

      ! adjust the 2d hydrostatic relative to the surface
      dhidtlm = dhidtlm - SPREAD(SPREAD(dhtdtl,1,n_vert),2,p_coeffs)
      if(PRESENT(dhtdtl0)) dhtdtl0 = sum(dhidtlm(tan_inds,:,:) &
        * SPREAD(eta_t,3,z_coeffs),dim=2)
      if(PRESENT(ddhtdhtdtl0)) ddhtdhtdtl0 = &
        sum(ddhidhidtl0(tan_inds,:,:) * SPREAD(eta_t,3,z_coeffs),dim=2)
      call Deallocate_test ( dhtdtl, 'dhtdtl', ModuleName )
    endif

    call Deallocate_test ( eta_t, 'eta_t', ModuleName )

    ss_htan = 0
    if(PRESENT(neg_h_tan)) then
      ss_htan = size(neg_h_tan)
      if(PRESENT(tan_phi_h_grid)) tan_phi_h_grid(1:ss_htan)=neg_h_tan-h_surf
    endif

    ! construct some basic matricies
    h_tans = SPREAD(h_grid_t,1,2*n_vert)
    phi_tans = SPREAD(phi_t,1,2*n_vert)

    nullify ( eta_t, h_grid_tt )
    call Allocate_test ( eta_t, n_vert, max (1,n_tan-ss_htan), 'n_vert', ModuleName )
    call Allocate_test ( h_grid_tt, n_vert, 'h_grid_tt', ModuleName )

    call get_eta_sparse(z_grid(tan_inds(min(n_tan,ss_htan+1):n_tan)),z_grid, &
      eta_t)
    h_grid_tt = matmul(eta_t,h_grid_t(min(n_tan,ss_htan+1):n_tan))

    ! basic pressure grid
    m_z_grid = SPREAD(z_grid((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/)),2,n_tan)

    ! determine condensed vector format so we don't mess with a bunch
    ! of zeros
    call Deallocate_test ( h_grid_t, 'h_grid_t', ModuleName )
    call Deallocate_test ( eta_t, 'eta_t', ModuleName )

    n_cvf = 2*sum(n_vert + 1 - tan_inds)

    nullify ( old_hts, cvf_h_tan, cvf_z_grid, cvf_ang_offset, cvf_sign )
    nullify ( h_zf, path_inds, vert_inds, near_inds, cvf_inds )

    call Allocate_test ( old_hts, n_cvf, 'old_hts', ModuleName )
    call Allocate_test ( cvf_h_tan, n_cvf, 'cvf_h_tan', ModuleName )
    call Allocate_test ( cvf_z_grid, n_cvf, 'cvf_z_grid', ModuleName )
    call Allocate_test ( cvf_ang_offset, n_cvf, 'cvf_ang_offset', ModuleName )
    call Allocate_test ( cvf_sign, n_cvf, 'cvf_sign', ModuleName )
    call Allocate_test ( h_zf, n_cvf, p_coeffs, 'h_zf', ModuleName )
    call Allocate_test ( path_inds, n_cvf, 'path_inds', ModuleName )
    call Allocate_test ( vert_inds, n_cvf, 'vert_inds', ModuleName )
    call Allocate_test ( near_inds, n_cvf/2, 'near_inds', ModuleName )
    call Allocate_test ( cvf_inds, n_cvf, 'cvf_inds', ModuleName )

    cvf_inds = PACK((/(i,i=1,2*n_vert*n_tan)/),                     &
    & RESHAPE(SPREAD((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/),2,n_tan), &
    &             (/2*n_vert*n_tan/)) -                             &
    & RESHAPE(SPREAD(tan_inds,1,2*n_vert),(/2*n_vert*n_tan/)) >= 0)

    ! sign of phi matrix
    cvf_sign = (-1.0_rp)**((modulo(cvf_inds-1,2*n_vert)/n_vert)+1)

    ! Apparently you can't sequentially access indicies in a multirank array
    ! in f90 like you can in idl which is a major inconvenience
    nullify ( mask )
    call Allocate_test ( mask, 2*n_vert*n_tan, 'mask', ModuleName )
    mask = .false.
    mask(cvf_inds) = .true.
    cvf_ang_offset = PACK(phi_tans,RESHAPE(mask,(/2*n_vert,n_tan/)))

    ! compute phi_s - phi_t
    j = n_vert
    do i = 1 , ss_htan
      cvf_ang_offset(j+1:j+n_vert) = phi_t(i) &
        -2.0_rp*acos((req + neg_h_tan(i)) / req)
      h_tans(:,i) = neg_h_tan(i)
      j = j + n_vert + n_vert
    end do
    h_grid = PACK(SPREAD(h_grid_tt((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/)), &
      2,n_tan),RESHAPE(mask,(/2*n_vert,n_tan/)))
    cvf_h_tan = PACK(h_tans,RESHAPE(mask,(/2*n_vert,n_tan/)))
    cvf_z_grid = PACK(m_z_grid,RESHAPE(mask,(/2*n_vert,n_tan/)))

    call Deallocate_test ( m_z_grid, 'm_z_grid', ModuleName )
    call Deallocate_test ( h_tans, 'h_tans', ModuleName )
    call Deallocate_test ( phi_tans, 'phi_tans', ModuleName )
    call Deallocate_test ( h_grid_tt, 'h_grid_tt', ModuleName )
    call Deallocate_test ( mask, 'mask', ModuleName )

    ! These are the tangent indicies on a cvf array
    if(n_tan > ss_htan) then
      ! some tangents above the earth surface
      nullify ( force_zero )
      call Allocate_test ( force_zero, 2*(n_tan-ss_htan), 'force_zero', ModuleName )
      force_zero(1) = n_vert*(2*ss_htan + 1) - tan_inds(1+ss_htan) + 1
      force_zero(2) = force_zero(1) + 1
      j = 3
      do i = 2,n_tan - ss_htan
        force_zero(j) = force_zero(j-2)+2*n_vert &
          - tan_inds(ss_htan+i-1)-tan_inds(ss_htan+i)+2
        force_zero(j+1) = force_zero(j)+1
        j = j + 2
      end do
    endif

    ! estimate the phi's using the tangent heights as a first guess
    old_hts = h_grid
    p_grid = cvf_ang_offset + cvf_sign &
      * acos((req + cvf_h_tan)/(req + old_hts))

    ! force the del phi at the tangent to be zero
    if(n_tan > ss_htan) p_grid(force_zero) = cvf_ang_offset(force_zero)

    ! estimate some new heights
    ! This is going to have some small surface errors away
    ! from the tangent but we will ignore this issue for now.
    ! the following lines replace an extremely sparse matrix
    ! multiplication with a less memory intensive alternative
    ! and probably faster alternative
    ! convert the cvf indicies into path indicies

    path_inds = modulo((cvf_inds - 1),2*n_vert) + 1
    vert_inds = path_inds - n_vert
    near_inds = PACK((/(i,i=1,n_cvf)/),path_inds <= n_vert)
    vert_inds(near_inds) = n_vert - path_inds(near_inds) + 1

    ! h_ref with surface adjustment
    h_zf = h_ref(vert_inds,:) - h_surf

    ! now these parts will change with iteration therefore these steps
    ! get repeated
    nullify ( mask, eta_p )
    call Allocate_test ( mask, n_cvf, 'mask', ModuleName )
    call Allocate_test ( eta_p, n_cvf, p_coeffs, 'eta_p', ModuleName )

    call get_eta_sparse(p_basis,p_grid,eta_p)
    h_grid = sum(h_zf * eta_p,dim=2)
    call Deallocate_test( eta_p, 'eta_p', ModuleName )

    mask = abs(old_hts-h_grid) > 0.01_rp
    no_of_bad_fits = count(mask)

    nullify ( junk )
    call Allocate_Test ( junk , no_of_bad_fits, 'junk', ModuleName )
    junk = PACK((/(i,i=1,n_cvf)/),mask)
    where(h_grid < cvf_h_tan) 
      h_grid = cvf_h_tan
    endwhere

    ! this construct is considered more desireable than a do while
    ! according to the fortran explained book

    iter = 0
    do

      iter = iter + 1
      if (no_of_bad_fits == 0 .or. iter == 20) EXIT
      old_hts = h_grid

      ! recompute phi_s - phi_t because we have to iterate on phi_s - phi_t
      j = n_vert
      do i = 1 , ss_htan
        cvf_ang_offset(j+1:j+n_vert) = phi_t(i) &
          -2.0_rp*acos((req + neg_h_tan(i)) / (req + h_grid(j)))
        j = j + n_vert + n_vert
      end do

      nullify ( arg, some_phi, h_better, eta_p )
      call Allocate_test ( arg, no_of_bad_fits, 'arg', ModuleName )
      call Allocate_test ( some_phi, no_of_bad_fits, 'some_phi', ModuleName )
      call Allocate_test ( h_better, no_of_bad_fits, 'h_better', ModuleName )
      call Allocate_test ( eta_p, no_of_bad_fits, p_coeffs, 'eta_p', ModuleName )

      arg = (req + cvf_h_tan(junk))/ (req + old_hts(junk))
      some_phi = cvf_ang_offset(junk) + cvf_sign(junk) * acos(arg)
      call get_eta_sparse(p_basis,some_phi,eta_p)
      h_better = sum(h_zf(junk,:) * eta_p,dim=2)
      h_better = max ( cvf_h_tan(junk), h_better )
      !       where (h_better < cvf_h_tan(junk) )
      !         h_better = cvf_h_tan(junk)
      !       endwhere
      h_grid(junk) = h_better
      call Deallocate_test ( junk, 'junk', ModuleName )
      call Deallocate_test ( h_better, 'h_better', ModuleName )
      call Deallocate_test ( arg, 'arg', ModuleName )
      call Deallocate_test ( some_phi, 'some_phi', ModuleName )
      call Deallocate_test ( eta_p, 'eta_p', ModuleName )

      mask = abs(old_hts-h_grid) > 0.01_rp
      no_of_bad_fits = count(mask)

      call Deallocate_test ( junk, 'junk', ModuleName )
      call Allocate_test ( junk, no_of_bad_fits, 'junk', ModuleName )
      junk = PACK((/(i,i=1,n_cvf)/),mask)

    end do

    call deallocate_test ( mask, 'mask', ModuleName )

    if (no_of_bad_fits > 0) then

      print *,'Warning: full convergence not acheived on:'
      print *,'implementing an improved approximation patch'
      print *,'path index, tangent index, error'

      ! we are going to assume that the tangent value is good
      ! The following is an F90 specific design which is quite different
      ! from the IDL code
      
      nullify ( path_ind, tan_ind )
      call Allocate_test ( path_ind, no_of_bad_fits, 'path_ind', ModuleName )
      call Allocate_test ( tan_ind, no_of_bad_fits, 'tan_ind', ModuleName )
      path_ind = modulo(cvf_inds(junk)-1,2*n_vert) + 1
      tan_ind = 1 + (cvf_inds(junk)-1) / (2*n_vert)
      st_ind = 1
      end_ind = 1
      do i = 1 , no_of_bad_fits
        print '(i5,1x,i5,1x,f7.3)',path_ind(i),tan_ind(i), &
          old_hts(junk(i))-h_grid(junk(i))
      end do

      do

        if(end_ind > no_of_bad_fits) EXIT

        do
          if(tan_ind(end_ind) > tan_ind(st_ind)) exit
          end_ind = end_ind + 1
        end do

        end_ind = end_ind - 1

        ! find which side of the tangent we are on
        if (path_ind(st_ind) < n_vert + 1) then
          ! this is the near observer side
          low_pt = junk(end_ind) + 1
          hi_pt  = junk(st_ind) - 1
        else
          ! this is the far observer side
          low_pt = junk(st_ind) - 1
          hi_pt  = junk(end_ind) + 1
        end if

        ! Correct
        h_grid(junk(st_ind):junk(end_ind)) = h_grid(low_pt) &
          + (h_grid(hi_pt) - h_grid(low_pt)) &
          * (cvf_z_grid(junk(st_ind):junk(end_ind)) - cvf_z_grid(low_pt)) &
          / (cvf_z_grid(hi_pt) - cvf_z_grid(low_pt))
        st_ind = end_ind + 1
        end_ind = st_ind

      end do

      print *,'completed re-estimation of bad points'

      call Deallocate_test ( path_ind, 'path_ind', ModuleName )
      call Deallocate_test ( tan_ind, 'tan_ind', ModuleName )

    endif

    ! deallocate the loop iteration stuff
    call Deallocate_test ( junk, 'junk', ModuleName )
    call Deallocate_test ( old_hts, 'old_hts', ModuleName )
    call Deallocate_test ( h_zf, 'h_zf', ModuleName )
    call Deallocate_test ( path_inds, 'path_inds', ModuleName )
    call Deallocate_test ( near_inds, 'near_inds', ModuleName )
    
    nullify ( not_zero_p )
    call Allocate_test ( not_zero_p, n_cvf, p_coeffs, 'not_zero_p', ModuleName )

    ! compute final set of angles
    p_grid = cvf_ang_offset + cvf_sign  &
      * acos((req + cvf_h_tan)/(req + h_grid))
    if(n_tan > ss_htan) then
      p_grid(force_zero) = cvf_ang_offset(force_zero)
      call Deallocate_test ( force_zero, 'force_zero', ModuleName )
    endif

    ! now compute the temperature grid
    nullify ( eta_p )
    call Allocate_test ( eta_p, n_cvf, p_coeffs, 'eta_p', ModuleName )

    call get_eta_sparse(p_basis,p_grid,eta_p,NOT_ZERO = not_zero_p)
    t_grid = sum(t_ref(vert_inds,:) * eta_p,dim=2)

    ! compute the vertical derivative grid
    dhitdzi = sum(dhidzij(vert_inds,:) * eta_p,dim=2)

    ! compute the temperature derivative grid

    if (PRESENT(dhitdtlm)) then

      nullify ( inds, eta_t, not_zero_t )
      call Allocate_test ( inds, n_cvf, 'inds', ModuleName )
      call Allocate_test ( eta_t, n_cvf, z_coeffs, 'eta_t', ModuleName )
      call Allocate_test ( not_zero_t, n_cvf, z_coeffs, 'not_zero_t', ModuleName )

      ! desparate attempt to try something different
      inds = modulo(cvf_inds-1,2*n_vert) - n_vert
      where(inds >= 0) inds = inds + 1
      inds = abs(inds)
      ! compute the path temperature noting where the zeros are
      call get_eta_sparse(z_basis,cvf_z_grid,eta_t,NOT_ZERO = not_zero_t)
!
! Initialize all ..
!
      do_calc_t = .false.
      do_calc_hyd = .false.

      eta_zxp = 0.0_rp
      dhitdtlm = 0.0_rp

      sv_t = 0
      do sv_p = 1 , p_coeffs
        do sv_z = 1 , z_coeffs
          sv_t = sv_t + 1
!          if(.NOT. t_deriv_flag(sv_t)) CYCLE
          where (not_zero_p(:,sv_p) .and. not_zero_t(:,sv_z))
            do_calc_t(:,sv_t) = .true.
            eta_zxp(:,sv_t) = eta_t(:,sv_z) * eta_p(:,sv_p)
          endwhere
          where (not_zero_p(:,sv_p) .and. dhidtlm(inds(:),sv_p,sv_z) > 0.0_rp)
            do_calc_hyd(:,sv_t) = .true.
            dhitdtlm(:,sv_t) = dhidtlm(inds(:),sv_p,sv_z) * eta_p(:,sv_p)
          endwhere
        end do
      end do

      call Deallocate_test ( inds, 'inds', ModuleName )
      call Deallocate_test ( eta_t, 'eta_t', ModuleName )
      call Deallocate_test ( not_zero_t, 'not_zero_t', ModuleName )

    endif
    
    call Deallocate_test ( eta_p, 'eta_p', ModuleName )
    call Deallocate_test ( cvf_ang_offset, 'cvf_ang_offset', ModuleName )
    call Deallocate_test ( cvf_z_grid, 'cvf_z_grid', ModuleName )
    call Deallocate_test ( cvf_sign, 'cvf_sign', ModuleName )
    call Deallocate_test ( cvf_h_tan, 'cvf_h_tan', ModuleName )
    call Deallocate_test ( cvf_inds, 'cvf_inds', ModuleName )
    call Deallocate_test ( vert_inds, 'vert_inds', ModuleName )
    call Deallocate_test ( not_zero_p, 'not_zero_p', ModuleName )

  end subroutine metrics

end module metrics_m
! $Log$
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
