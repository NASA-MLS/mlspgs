! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Metrics_m

  implicit NONE
  private
  public :: Metrics

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  integer, save :: CalledTimes = 0

contains

  ! ----------------------------------------------------  Metrics  -----

  subroutine Metrics ( phi_t, tan_inds, p_basis, z_grid, h_ref, t_ref, dhidzij,&
          &  beta, t_deriv_flag, h_grid, p_grid, t_grid, dhitdzi, req, dhidtlm,&
          &  ddhidhidtl0, dhitdtlm, eta_zxp, tan_phi_h_grid, tan_phi_t_grid,   &
          &  dhtdzt, dhtdtl0, ddhtdhtdtl0, neg_h_tan, z_basis, do_calc_t,      &
          &  do_calc_hyd )

    ! The goal of this program is to return a matrix of h_grids
    ! and t_grids which define 2 d integration paths

    ! The main change in this routine is to interchange the z_coeffs and
    ! p_coeffs so as to make it consistent with the radiative transfer program.

    use Allocate_deallocate, only: Allocate_test, Deallocate_test
    use Geometry, only: EarthRadA, EarthRadB
    use Get_Eta_Matrix_m, only: Get_Eta_Sparse
    use MLSCommon, only: RP, IP
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use Output_m, only: OUTPUT
    use Toggles, only: Emit, Toggle

    ! inputs:

    real(rp), intent(in) :: phi_t(:)   !orbit projected tangent geodetic angles
    integer(ip), intent(in) :: tan_inds(:)!tangent height indicies
    real(rp), intent(in) :: p_basis(:) !horizontal temperature representation
    !                                     basis
    real(rp), intent(in) :: z_grid(:)  !pressures for which heights/temps are
    !                                     needed
    real(rp), intent(in) :: h_ref(:,:) ! heights by t_phi_basis
    real(rp), intent(in) :: t_ref(:,:) ! temperatures by t_phi_basis
    real(rp), intent(in) :: dhidzij(:,:)! vertical derivative by t_phi_basis
    real(rp), intent(in) :: beta       !orbital incline angle (Radians)
    logical,  intent(in) :: t_deriv_flag(:)  ! User's deriv. flags for Temp.
    ! outputs:
    real(rp), intent(out) :: h_grid(:) !computed heights
    real(rp), intent(out) :: p_grid(:) !computed phi's
    real(rp), intent(out) :: t_grid(:) !computed temperatures
    real(rp), intent(out) :: dhitdzi(:)!derivative of height wrt zeta
    !                                   --may be useful in future computations
    real(rp), intent(out) :: req       !equivalent elliptical earth radius

    ! Keywords:
    ! optional inputs
    real(rp), optional, intent(in) :: z_basis(:) ! vertical temperature basis
    real(rp), optional, intent(inout) :: dhidtlm(:,:,:) ! reference temperature
    !   derivatives. This gets adjusted so that at ref_h(1,@tan phi)) is 0.0 for
    !   all temperature coefficients. This is height X zeta_basis X phi_basis
    real(rp), optional, intent(in) :: ddhidhidtl0(:,:,:)! second order reference
    !   temperature derivatives. This is (height, phi_basis, zeta_basis)
    real(rp), optional, intent(in) :: neg_h_tan(:)  !sub earth tangent heights
    ! optional outputs
    real(rp), optional, intent(out) :: tan_phi_h_grid(:)!heights along the
    !                                                      tangent
    real(rp), optional, intent(out) :: tan_phi_t_grid(:)!temperature along the
    !                                                      tangent
    real(rp), optional, intent(out) :: dhtdzt(:)   ! height derivative wrt
    !                                                pressure along the tangent
    real(rp), optional, intent(out) :: dhtdtl0(:,:) ! First order derivative at
    !                                                the tangent
    real(rp), optional, intent(out) :: ddhtdhtdtl0(:,:) ! Second order 
    !          derivative at the tangent only---used for antenna affects
    !
    real(rp), optional, intent(out) :: dhitdtlm(:,:)
    !                             derivative of path position wrt temperature
    !                             statevector (z_basisXphi_basis)
    real(rp), optional, intent(out) :: eta_zxp(:,:) ! eta matrix for temperature.

    logical, optional, intent(out) :: do_calc_t(:,:) ! non zero locater for
    !                                       temperature bases computations.
    logical, optional, intent(out) :: do_calc_hyd(:,:) !non zero locator for
    !                                       hydrostatic calculations.

    ! NOTES
    ! p_basis and p_grid are phi's in offset degrees relative to phi_t, that
    ! is the phi_t p_basis or p_grid = 0.0 is phi_t.
    ! The phi basis is wholly independent of phi_t
    ! compute primary hydrostatic grid
    ! tangent phi index

    ! ===>> EVERYTHING comes out in cvf format <<===!

    ! Local variables.  CAN WE GET SOME COMMENTS FOR EACH OF THESE?
    integer :: END_IND
    integer :: END_IND1
    integer :: HI_PT
    integer :: I
    integer :: ITER
    integer :: J
    integer :: LOW_PT
    integer :: N_CVF
    integer :: NO_OF_BAD_FITS
    integer :: N_TAN
    integer :: N_VERT
    integer :: P_COEFFS
    integer :: Path_Ind
    integer :: SS_HTAN
    integer :: ST_IND
    integer :: SV_P
    integer :: SV_T
    integer :: SV_Z
    integer :: Z_COEFFS
    integer :: TAN_PTR

    real(rp) :: CP2
    real(rp) :: CSQ
    real(rp) :: H_SURF
    real(rp) :: SP2

    logical, dimension(:), pointer :: MASK
    logical, dimension(:,:), pointer :: MASK2
    logical, dimension(:,:), pointer :: NOT_ZERO_P
    logical, dimension(:,:), pointer :: NOT_ZERO_T

    integer, dimension(:), pointer :: CVF_INDS
    integer, dimension(:), pointer :: FORCE_ZERO
    integer, dimension(:), pointer :: INDS
    integer, dimension(:), pointer :: JUNK, JUN ! JUN => some_of_JUNK
    integer, dimension(:), pointer :: NEAR_INDS
    integer, dimension(:), pointer :: PATH_INDS
    integer, dimension(:), pointer :: VERT_INDS

    real(rp), dimension(:), pointer :: CVF_ANG_OFFSET
    real(rp), dimension(:), pointer :: CVF_H_TAN
    real(rp), dimension(:), pointer :: CVF_SIGN
    real(rp), dimension(:), pointer :: CVF_Z_GRID
    real(rp), dimension(:), pointer :: DHTDTL
    real(rp), dimension(size(tan_inds)) :: H_GRID_T
    real(rp), dimension(:), pointer :: H_GRID_TT
    real(rp), dimension(:), pointer :: OLD_HTS

    real(rp), dimension(:,:), pointer :: ETA_P, ET_P ! ET_P => some_of_ETA_P
    real(rp), dimension(:,:), pointer :: ETA_T
    real(rp), dimension(2*size(z_grid),size(tan_inds)) :: H_TANS
    real(rp), dimension(:,:), pointer :: H_ZF
    real(rp), dimension(2*size(z_grid),size(tan_inds)) :: M_Z_GRID
    real(rp), dimension(2*size(z_grid),size(tan_inds)) :: PHI_TANS

    ! Begin program

    calledTimes = calledTimes + 1
    p_coeffs = size(p_basis)
    n_vert = size(z_grid)
    n_tan = size(tan_inds)

    csq = (earthrada * earthradb)**2 / &
          &  ((earthrada**2-earthradb**2)*sin(beta)**2 + earthradb**2)

    nullify ( eta_t )
    call allocate_test ( eta_t, n_tan, p_coeffs, 'eta_t', ModuleName )

    ! compute the tangent height vertical
    call get_eta_sparse ( p_basis,phi_t, eta_t )
    h_grid_t = sum(h_ref(tan_inds,:)*eta_t,dim=2)
    h_surf = sum(h_ref(1,:)*eta_t(1,:))

    ! compute equivalent earth radius at phi_t(1), nearest surface
    cp2 = cos(phi_t(1))**2
    sp2 = 1.0_rp - cp2
    req = 0.001_rp*sqrt((earthrada**4*sp2 + csq**2*cp2) / &
                      & (earthrada**2*cp2 + csq*sp2))

    ! now for the path quantities. for simplicity, we are going to set the
    ! surface reference at the input z_grid(1) and adjust the req and the
    ! h_grid_t relative to this
    req = req + h_surf

    ! adjust h_ref accordingly
    h_grid_t = h_grid_t - h_surf
    if ( present(tan_phi_h_grid) ) tan_phi_h_grid = h_grid_t
    if ( present(tan_phi_t_grid) ) &
      tan_phi_t_grid = sum(t_ref(tan_inds,:)*eta_t,dim=2)
    if ( present(dhtdzt) ) dhtdzt = sum(dhidzij(tan_inds,:)*eta_t,dim=2)

    ! compute tangent temperature derivative
    if ( present(dhidtlm) ) then
      z_coeffs = size(z_basis)
      nullify ( dhtdtl )
      call allocate_test ( dhtdtl, z_coeffs, 'dhtdtl', ModuleName )
!     dhtdtl = sum(RESHAPE(dhidtlm(1,:,:),(/z_coeffs,p_coeffs/)) &
!       * SPREAD(RESHAPE(eta_t(1,:),(/p_coeffs/)),1,z_coeffs),dim=2)
      dhtdtl = sum(dhidtlm(1,:,:) * SPREAD(eta_t(1,:),1,z_coeffs),dim=2)

      ! adjust the 2d hydrostatic relative to the surface
      dhidtlm = dhidtlm - SPREAD(SPREAD(dhtdtl,1,n_vert),3,p_coeffs)
!     dhidtlm = dhidtlm - SPREAD(SPREAD(dhtdtl,1,n_vert),2,p_coeffs)

      j = z_coeffs * p_coeffs
      dhtdtl0 = RESHAPE(dhidtlm(tan_inds,:,:) * SPREAD(eta_t,2,z_coeffs),&
                     & (/n_tan, j/))

      ddhtdhtdtl0 = RESHAPE( &
                   ddhidhidtl0(tan_inds,:,:) * SPREAD(eta_t,2,z_coeffs), &
                     & (/n_tan, j/))

      call deallocate_test ( dhtdtl, 'dhtdtl', ModuleName )

    end if

    call deallocate_test ( eta_t, 'eta_t', ModuleName )

    ss_htan = 0
    if ( present(neg_h_tan) ) then
      ss_htan = size(neg_h_tan)
      if ( present(tan_phi_h_grid) ) tan_phi_h_grid(1:ss_htan)=neg_h_tan-h_surf
    end if

    ! construct some basic matrices
    h_tans = SPREAD(h_grid_t,1,2*n_vert)
    phi_tans = SPREAD(phi_t,1,2*n_vert)

    nullify ( eta_t, h_grid_tt )
    call allocate_test ( h_grid_tt, n_vert, 'h_grid_tt', ModuleName )

    call allocate_test ( eta_t, n_vert, max (1,n_tan-ss_htan), 'n_vert', ModuleName )
    call get_eta_sparse ( z_grid(tan_inds(min(n_tan,ss_htan+1):n_tan)), z_grid, &
                        & eta_t )
    h_grid_tt = matmul(eta_t,h_grid_t(min(n_tan,ss_htan+1):n_tan))
    call deallocate_test ( eta_t, 'eta_t', ModuleName )

    ! basic pressure grid
    m_z_grid = SPREAD(z_grid((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/)),2,n_tan)

    ! determine condensed vector format so we don't mess with a bunch
    ! of zeros

    n_cvf = 2*sum(n_vert + 1 - tan_inds)

    nullify ( old_hts, cvf_h_tan, cvf_z_grid, cvf_ang_offset, cvf_sign )
    nullify ( h_zf, cvf_inds )

    call allocate_test ( old_hts, n_cvf, 'old_hts', ModuleName )
    call allocate_test ( cvf_h_tan, n_cvf, 'cvf_h_tan', ModuleName )
    call allocate_test ( cvf_z_grid, n_cvf, 'cvf_z_grid', ModuleName )
    call allocate_test ( cvf_ang_offset, n_cvf, 'cvf_ang_offset', ModuleName )
    call allocate_test ( cvf_sign, n_cvf, 'cvf_sign', ModuleName )
    call allocate_test ( h_zf, n_cvf, p_coeffs, 'h_zf', ModuleName )
    call allocate_test ( cvf_inds, n_cvf, 'cvf_inds', ModuleName )

    cvf_inds = PACK((/(i,i=1,2*n_vert*n_tan)/),                     &
    & RESHAPE(SPREAD((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/),2,n_tan), &
    &             (/2*n_vert*n_tan/)) -                             &
    & RESHAPE(SPREAD(tan_inds,1,2*n_vert),(/2*n_vert*n_tan/)) >= 0)

    ! sign of phi matrix
    cvf_sign = (-1.0_rp)**((modulo(cvf_inds-1,2*n_vert)/n_vert)+1)

    ! Apparently you can't sequentially access elements of a multirank array
    ! in f90 like you can in idl, which is a major inconvenience
    nullify ( mask2 )
    call allocate_test ( mask2, 2*n_vert, n_tan, 'mask2', ModuleName )
    mask2 = .false.
    ! mask2(cvf_inds) = .true.
    do i = 1, size(cvf_inds)
      mask2(mod(cvf_inds(i)-1,2*n_vert)+1, (cvf_inds(i)-1)/(2*n_vert) + 1) &
        & = .true.
    end do
    cvf_ang_offset = PACK(phi_tans,mask2)
    h_grid = PACK(SPREAD(h_grid_tt((/(i,i=n_vert,1,-1),(i,i=1,n_vert)/)), &
                      &  2,n_tan),mask2)
    call deallocate_test ( h_grid_tt, 'h_grid_tt', ModuleName )

    ! compute phi_s - phi_t
    j = n_vert
    do i = 1 , ss_htan
      cvf_ang_offset(j+1:j+n_vert)=phi_t(i)-2.0_rp*Acos((req+neg_h_tan(i))/req)
      h_tans(:,i) = neg_h_tan(i)
      j = j + n_vert + n_vert
    end do
    cvf_h_tan = PACK(h_tans,mask2)
    cvf_z_grid = PACK(m_z_grid,mask2)

    call deallocate_test ( mask2, 'mask', ModuleName )

    ! These are the tangent indicies on a cvf array
    if ( n_tan > ss_htan ) then
      ! some tangents above the earth surface
      nullify ( force_zero )
      call allocate_test ( force_zero, 2*(n_tan-ss_htan), 'force_zero', ModuleName )
      force_zero(1) = n_vert*(2*ss_htan + 1) - tan_inds(1+ss_htan) + 1
      force_zero(2) = force_zero(1) + 1
      j = 3
      do i = 2, n_tan - ss_htan
        force_zero(j) = force_zero(j-2)+2*n_vert - &
                      & tan_inds(ss_htan+i-1)-tan_inds(ss_htan+i)+2
        force_zero(j+1) = force_zero(j)+1
        j = j + 2
      end do
    end if

    ! estimate the phi's using the tangent heights as a first guess
    old_hts = h_grid
    p_grid = cvf_ang_offset + cvf_sign * Acos((req+cvf_h_tan)/(req+old_hts))

    ! force the del phi at the tangent to be zero
    if ( n_tan > ss_htan ) p_grid(force_zero) = cvf_ang_offset(force_zero)

    ! estimate some new heights
    ! This is going to have some small surface errors away
    ! from the tangent but we will ignore this issue for now.
    ! the following lines replace an extremely sparse matrix
    ! multiplication with a less memory intensive alternative
    ! and probably faster alternative
    ! convert the cvf indicies into path indicies

    nullify ( path_inds, vert_inds, near_inds )
    call allocate_test ( vert_inds, n_cvf, 'vert_inds', ModuleName )
    call allocate_test ( path_inds, n_cvf, 'path_inds', ModuleName )
    call allocate_test ( near_inds, n_cvf/2, 'near_inds', ModuleName )

    path_inds = modulo((cvf_inds - 1),2*n_vert) + 1
    vert_inds = path_inds - n_vert
    near_inds = PACK((/(i,i=1,n_cvf)/),path_inds <= n_vert)
    vert_inds(near_inds) = n_vert - path_inds(near_inds) + 1

    call deallocate_test ( near_inds, 'near_inds', ModuleName )
    call deallocate_test ( path_inds, 'path_inds', ModuleName )

    ! h_ref with surface adjustment
    h_zf = h_ref(vert_inds,:) - h_surf

    ! now these parts will change with iteration therefore these steps
    ! get repeated
    nullify ( mask, eta_p, junk )
    call allocate_test ( junk, n_cvf, 'junk', moduleName )
    call allocate_test ( mask, n_cvf, 'mask', ModuleName )
    call allocate_test ( eta_p, n_cvf, p_coeffs, 'eta_p', ModuleName )

    call get_eta_sparse ( p_basis, p_grid, eta_p )
    h_grid = max(cvf_h_tan, sum(h_zf * eta_p,dim=2))

    ! this construct is considered more desireable than a do while
    ! according to the Fortran explained book

    iter = 0
    do

      iter = iter + 1
      mask = abs(old_hts-h_grid) > 0.01_rp
      no_of_bad_fits = count(mask)
      if ( no_of_bad_fits == 0 ) exit
      jun => junk(:no_of_bad_fits)
      jun = PACK((/(i,i=1,n_cvf)/),mask)
      if ( iter == 20 ) exit

      old_hts = h_grid
      et_p => eta_p(:no_of_bad_fits,:)

      ! recompute phi_s - phi_t because we have to iterate on phi_s - phi_t
      j = n_vert
      do i = 1 , ss_htan
        cvf_ang_offset(j+1:j+n_vert) = phi_t(i) - &
                  & 2.0_rp*Acos((req+neg_h_tan(i))/(req+h_grid(j)))
        j = j + n_vert + n_vert
      end do

      call get_eta_sparse ( p_basis, cvf_ang_offset(jun) + cvf_sign(jun) &
        & * Acos((req + cvf_h_tan(jun))/ (req + old_hts(jun))), et_p )
      h_grid(jun) = max ( cvf_h_tan(jun), sum(h_zf(jun,:) * et_p, dim=2) )

    end do

    call deallocate_test ( eta_p, 'eta_p', ModuleName )
    call deallocate_test ( mask, 'mask', ModuleName )

    if ( no_of_bad_fits > 0 ) then

      if ( toggle(emit) ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
          & 'Full convergence not achieved, implementing an improved approximation patch' )
        call output ( 'pth ind, error', advance='yes' )
      end if

      ! we are going to assume that the tangent value is good
      ! The following is an F90 specific design which is quite different
      ! from the IDL code

      ! find ranges of contigous indicies
      st_ind = 1
      path_ind = modulo(cvf_inds(jun(st_ind))-1,2*n_vert) + 1
      do end_ind = 1, no_of_bad_fits
        if ( end_ind < no_of_bad_fits ) then
          if ( jun(end_ind+1) - jun(end_ind) < 2 ) cycle
        end if

        if ( toggle(emit) ) then
          call output ( path_ind, places=7 )
          call output ( ', ' )
          call output ( old_hts(jun(st_ind))-h_grid(jun(st_ind)), advance='yes' )
        end if

        ! find which side of the tangent we are on
        if ( path_ind < n_vert + 1 ) then ! near observer side
          low_pt = jun(end_ind) + 1
          hi_pt  = MAX(jun(st_ind) - 1,1)
        else                              ! far observer side
          low_pt = jun(st_ind) - 1
          hi_pt  = MIN(jun(end_ind) + 1,n_cvf)
        end if
        tan_ptr = 2 * n_vert * ((cvf_inds(low_pt)-1) / (2*n_vert) + 1)
        IF (cvf_inds(jun(st_ind)) == tan_ptr - 2*n_vert + 1) THEN
          CALL mlsmessage(mlsmsg_warning,modulename,'resorting to 1d option in metrics')
! resort to 1 d equvalent
          CALL allocate_test(eta_p,1,p_coeffs,'eta_p',modulename)
          CALL get_eta_sparse(p_basis,(/phi_t((cvf_inds(low_pt)-1) &
               & / (2*n_vert) + 1)/), eta_p)
! this should be eta_p(1,p_coeffs)
          DO i = st_ind, end_ind
            h_grid(jun(i)) = SUM(h_ref(n_vert - MODULO(cvf_inds(jun(i))-1, &
            & 2*n_vert),:)*eta_p(1,:)) - h_surf
          enddo
          call deallocate_test ( eta_p, 'eta_p', ModuleName )
        ELSE IF (cvf_inds(jun(end_ind)) > tan_ptr - 1) THEN
! calculate the path ending index.
          CALL mlsmessage(mlsmsg_warning,modulename,'resorting to 1d option in metrics')
          end_ind1 = st_ind + tan_ptr - cvf_inds(jun(st_ind))
! resort to 1 d equvalent
          CALL allocate_test(eta_p,1,p_coeffs,'eta_p',modulename)
          CALL get_eta_sparse(p_basis,(/phi_t((cvf_inds(low_pt) - 1) &
          & / (2*n_vert) + 1)/), eta_p)
          DO i = st_ind, end_ind1
            h_grid(jun(i)) = SUM(h_ref(MODULO(cvf_inds(jun(i))-1,2*n_vert) &
                          & + 1 - n_vert,:)*eta_p(1,:)) - h_surf
          enddo
! This is in case the anomaly wraps to the next higher tangent level.
          IF (end_ind1 < end_ind) THEN
! resort to 1 d equvalent
            CALL get_eta_sparse(p_basis,(/phi_t((cvf_inds(low_pt)-1) &
                 & / (2*n_vert) + 2)/), eta_p)
! this should be eta_p(1,p_coeffs)
            DO i = end_ind1+1, end_ind
              h_grid(jun(i)) = SUM(h_ref(n_vert - MODULO(cvf_inds(jun(i))-1, &
              & 2*n_vert),:)*eta_p(1,:)) - h_surf
            enddo
          endif
          call deallocate_test ( eta_p, 'eta_p', ModuleName )
        else

          ! Correct
          h_grid(jun(st_ind):jun(end_ind)) = h_grid(low_pt) + &
             & (h_grid(hi_pt) - h_grid(low_pt)) * &
             & (cvf_z_grid(jun(st_ind):jun(end_ind))-cvf_z_grid(low_pt)) / &
             & (cvf_z_grid(hi_pt) - cvf_z_grid(low_pt))
        endif
        st_ind = end_ind + 1
        if ( st_ind < no_of_bad_fits ) &
          & path_ind = modulo(cvf_inds(jun(st_ind))-1,2*n_vert) + 1
      end do

    end if

    ! deallocate the loop iteration stuff
    call deallocate_test ( junk, 'junk', ModuleName )
    call deallocate_test ( old_hts, 'old_hts', ModuleName )
    call deallocate_test ( h_zf, 'h_zf', ModuleName )

    nullify ( not_zero_p )
    call allocate_test ( not_zero_p, n_cvf, p_coeffs, 'not_zero_p', ModuleName )

    ! compute final set of angles
    p_grid = cvf_ang_offset + cvf_sign * Acos((req+cvf_h_tan)/(req+h_grid))
    if ( n_tan > ss_htan ) then
      p_grid(force_zero) = cvf_ang_offset(force_zero)
      call deallocate_test ( force_zero, 'force_zero', ModuleName )
    end if

    ! now compute the temperature grid
    nullify ( eta_p )
    call allocate_test ( eta_p, n_cvf, p_coeffs, 'eta_p', ModuleName )

    call get_eta_sparse ( p_basis, p_grid, eta_p, NOT_ZERO = not_zero_p )
    t_grid = sum(t_ref(vert_inds,:) * eta_p,dim=2)

    ! compute the vertical derivative grid
    dhitdzi = sum(dhidzij(vert_inds,:) * eta_p,dim=2)

    ! compute the temperature derivative grid

    if ( present(dhitdtlm) ) then

      nullify ( inds, eta_t, not_zero_t )
      call allocate_test ( inds, n_cvf, 'inds', ModuleName )
      call allocate_test ( eta_t, n_cvf, z_coeffs, 'eta_t', ModuleName )
      call allocate_test ( not_zero_t, n_cvf, z_coeffs, 'not_zero_t', ModuleName )

      ! desparate attempt to try something different
      inds = modulo(cvf_inds-1,2*n_vert) - n_vert
      where ( inds >= 0 ) inds = inds + 1
      inds = abs(inds)
      ! compute the path temperature noting where the zeros are
      call get_eta_sparse ( z_basis, cvf_z_grid, eta_t, NOT_ZERO = not_zero_t )

      sv_t = 0
      do sv_p = 1 , p_coeffs
        do sv_z = 1 , z_coeffs
          sv_t = sv_t + 1
          if ( t_deriv_flag(sv_t) ) then
            do_calc_t(:,sv_t) = not_zero_t(:,sv_z) .and. not_zero_p(:,sv_p)
            where ( do_calc_t(:,sv_t) )
              eta_zxp(:,sv_t) = eta_t(:,sv_z) * eta_p(:,sv_p)
            elsewhere
              eta_zxp(:,sv_t) = 0.0
            end where
            do_calc_hyd(:,sv_t) = not_zero_p(:,sv_p) .and. dhidtlm(inds(:),sv_z,sv_p) > 0.0_rp
            where ( do_calc_hyd(:,sv_t) )
              dhitdtlm(:,sv_t) = dhidtlm(inds(:),sv_z,sv_p) * eta_p(:,sv_p)
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

      call deallocate_test ( not_zero_t, 'not_zero_t', ModuleName )
      call deallocate_test ( eta_t, 'eta_t', ModuleName )
      call deallocate_test ( inds, 'inds', ModuleName )

    end if

    call deallocate_test ( eta_p, 'eta_p', ModuleName )
    call deallocate_test ( cvf_ang_offset, 'cvf_ang_offset', ModuleName )
    call deallocate_test ( cvf_z_grid, 'cvf_z_grid', ModuleName )
    call deallocate_test ( cvf_sign, 'cvf_sign', ModuleName )
    call deallocate_test ( cvf_h_tan, 'cvf_h_tan', ModuleName )
    call deallocate_test ( cvf_inds, 'cvf_inds', ModuleName )
    call deallocate_test ( vert_inds, 'vert_inds', ModuleName )
    call deallocate_test ( not_zero_p, 'not_zero_p', ModuleName )

  end subroutine Metrics

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Metrics_m

! $Log$
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
