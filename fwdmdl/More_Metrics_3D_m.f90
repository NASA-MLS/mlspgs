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
module More_Metrics_3D_m
!=============================================================================

  implicit NONE
  private

  public :: More_Metrics_3D

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine More_Metrics_3D ( &
          ! Inputs
          & S, N_Tan, T_Ref, dHidZij, &
          ! Outputs
          & T_Path, dHitdZi, &
          ! Optional inputs
          & QTM_Tree, ddHidHidTl0, dHidTlm, T_Sv, Z_Ref, &
          ! Optional outputs
          & ddHtdHtdTl0, dHitdTlm, dHtdTl0, dHtdZt, Eta_zQT, &
          & Do_Calc_T, Do_Calc_Hyd, Tan_T )

    use Generate_QTM_m, only: QTM_Tree_t
    use Indexed_Values_m, only: Value_QTM_2D_List_t
    use Load_Sps_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use QTM_Interpolation_Weights_3D_m, only: QTM_Interpolation_Weights, &
      & S_QTM_t

    ! Inputs

    type(S_QTM_t), intent(in) :: S(:)    ! Intersections of path with constant-zeta
    integer, intent(in) :: N_Tan         ! Tangent index in S
    real(rp), intent(in) :: T_Ref(:,:)   ! Temperatures at Z_Ref X hGrid
                                         ! adjacent to the path
    real(rp), intent(in) :: dHidZij(:,:) ! Vertical derivative at Z_Ref X hGrid
                                         ! adjacent to the path

    ! Outputs

    real(rp), intent(out) :: T_Path(:)   ! Computed temperatures
    real(rp), intent(out) :: dHitdZi(:)  ! Derivative of height wrt zeta --
                                         ! may be useful in future computations

    ! Optional inputs for temperature derivatives

    type(QTM_tree_t), optional, intent(in) :: QTM_Tree
    real(rp), optional, intent(in) :: ddHidHidTl0(:,:,:) ! second order
    !          derivatives of height w.r.t T_Ref. This is an extract of
    !          (height, zeta_basis, QTM) adjacent to the path.
    !          Needed only if present(dHidTlm).
    real(rp), optional, intent(inout) :: dHidTlm(:,:,:) ! reference temperature
    !          derivatives. This gets adjusted so that at ref_h(1,@tan)) is
    !          0.0 for all temperature coefficients. This is an extract of
    !          height X zeta_basis X QTM adjacent to the path.
    type(grids_t), optional, intent(in) :: T_Sv    ! State vector temperature
    !          Needed only if present(dHidTlm).
    real(rp), optional, intent(in) :: Z_Ref(:)     ! -log pressures (zetas) for
    !          which derivatives are needed.  Only the parts from the tangent
    !          outward are used.  Needed only if present(dHidTlm).

    ! Optional outputs

    real(rp), optional, intent(out), target :: ddHtdHtdTl0(:)  ! Second
    !          order derivatives of height w.r.t T_Ref at the tangent only --
    !          used for antenna affects. Computed if present(dHidTlm).
    real(rp), optional, intent(out), target :: dHitdTlm(:,:,:)   ! Derivative of
    !          path position wrt temperature state vector
    !          (S X t_sv%zet_basis X # vertices adjacent to path)
    real(rp), optional, intent(out), target :: dHtdTl0(:)      ! First order
    !          derivatives of height w.r.t T_Ref at the tangent only.  Computed
    !          if present(dHidTlm).
    real(rp), optional, intent(out) :: dHtdZt      ! Height derivative wrt
    !          pressure at the tangent.  Computed if present(dHidTlm).
    type(Value_QTM_2D_List_t), optional, intent(out) :: Eta_zQT(:) ! Interpolation
    !          coefficients from QTM to path for temperature.
    logical, optional, intent(out) :: Do_Calc_T(:,:,:)     ! Eta_zQT /= 0 on
    !          path Z X temperature Z X temperature QTM grid
    logical, optional, intent(out) :: Do_Calc_Hyd(:,:,:)   ! dHitdTlm /= 0 on
    !          path Z X temperature Z X temperature QTM grid
    real(rp), optional, intent(out) :: Tan_T ! temperature at the tangent

    integer :: I
    integer :: IZ     ! Zeta index in Z_Ref, first subscript of T_Ref, dHidZij
    integer :: NC     ! Number of S(i)%Coeff%v to use
    integer :: N_Path ! Size(S)

    n_path = size(s)

    ! Path quantities.  These involve only horizontal interpolation because
    ! temperature and dHidZij were interpolated onto the GL grid by
    ! Two_D_Hydrostatic. The horizontal interpolation coefficients from
    ! vertices adjacent to the path onto the path, were computed and put into
    ! S%Coeff%v by Metrics_3D.
    do i = 1, n_path
      nc = s(i)%coeff%n
      iz = max(s(i)%h_ind,1) ! h_ind == zero for reflection below H_GLGrid
      t_path(i) =  dot_product ( t_ref(iz,s(i)%coeff%v(:nc)%jp),   &
                               & s(i)%coeff%v(:nc)%v )
      dHitdZi(i) = dot_product ( dHidZij(iz,s(i)%coeff%v(:nc)%jp), &
                               & s(i)%coeff%v(:nc)%v )
    end do

    ! Optional tangent quantities.
    if ( n_tan <= n_path ) then
      if ( present(tan_t) ) tan_t = t_path(n_tan)
      if ( present(dHtdZt) ) dHtdZt = dHitdZi(n_tan)
      if ( present(dHidTlm) ) call Temperature_Derivatives ( s(n_tan) )
    end if

  contains

    subroutine Temperature_Derivatives ( TP )

      type(S_QTM_t), intent(in) :: TP  ! Tangent point

      real(rp), pointer :: ddH2(:,:) ! Two-d (Z x H) pointer for ddHtdHtdTl0
      logical, pointer :: dF2(:,:)   ! Two-d (Z x H) pointer for t_sv%deriv_flags
      real(rp), pointer :: dH2(:,:)  ! Two-d (Z x H) pointer for dHtdTl0
      integer :: IC     ! Index of a coefficient in S, 1...S(ip)%nc
      integer :: IP     ! Index of vertex adjacent to path
      integer :: IS     ! Index of point on path
      integer :: IZ     ! First (vertical) subscript at tangent
      integer :: NC     ! Number of interpolating coefficients at the tangent
      integer :: P_Coeffs   ! t_sv%l_p(1)
      integer :: Sv_Z       ! Zeta subscript for t_sv%deriv_flags
      real(rp) :: W         ! Horizontal interpolating coefficient from one
                            ! profile to one path point
      integer :: Z_Coeffs   ! t_sv%l_z(1) = size(T_ref,1)

      ! Adjust the 2d hydrostatic temperature derivative relative to the
      ! surface. Even though this is updated on every invocation, that is,
      ! with a new S, it works as if the original value were updated with
      ! the current S, because the interpolation represented by S is
      ! linear. Thus, the effect of cumulative updates for each new S are
      ! the same as starting from the original dHidTlm and updating with the
      ! latest S.  The algebra is horrible, but Maple has verified this.

      p_coeffs = t_sv%l_p(1)
      z_coeffs = t_sv%l_z(1)

      nc = tp%coeff%n
      do sv_z = 1, z_coeffs
        dHidTlm(:,sv_z,:) = dHidTlm(:,sv_z,:) - &
          & dot_product ( dHidTlm(1,sv_z,tp%coeff%v(:nc)%jp), &
                        & tp%coeff%v(:nc)%v )
      end do
 
      iz = max(tp%h_ind,1) ! h_ind == zero for reflection below H_GLGrid
      dHtdTl0 = 0
      ddHtdHtdTl0 = 0
      dH2(1:z_coeffs,1:size(t_ref,2)) => dHtdTl0
      ddH2(1:z_coeffs,1:size(t_ref,2)) => ddHtdHtdTl0
      do sv_z = 1, z_coeffs
        dH2(sv_z,tp%coeff%v(:nc)%jp) = &
          & dHidTlm(iz,sv_z,tp%coeff%v(:nc)%jp) * tp%coeff%v(:nc)%v
        ddH2(sv_z,tp%coeff%v(:nc)%jp) = &
          & ddHidHidTl0(iz,sv_z,tp%coeff%v(:nc)%jp) * tp%coeff%v(:nc)%v
      end do

      dF2(1:z_coeffs,1:p_coeffs) => t_sv%deriv_flags
      dHitdTlm = 0          ! Sparse representation might be desirable some day
      do_calc_hyd = .false. ! Sparse representation might be desirable some day
      do is = 1, n_path             ! Path length
        iz = max(s(is)%h_ind,1)     ! h_ind == zero for reflection below H_GLGrid
        do ic = 1, s(ip)%coeff%n    ! # horizontal interpolation coefficients
          ip = s(ip)%coeff%v(ic)%jp ! Index among path-adjacent profiles
          w = s(is)%coeff%v(ic)%v   ! Horizontal interpolation coefficient from
                                    ! profile IP to path point IS
          do sv_z = 1, z_coeffs     ! Zeta levels of profiles
            if ( dF2(sv_z,s(is)%coeff%v(ic)%j) ) then ! Derivative flag for temperature
              dHitdTlm(is,sv_z,ip) = &
                & max(dHidTlm(iz,sv_z,ip),0.0_rp) * w
              do_calc_hyd(is,sv_z,ip) = dHitdTlm(is,sv_z,ip) /= 0
            end if
          end do
        end do
      end do

      ! Get coefficients to interpolate from temperature (zeta, adjacent-to-path)
      ! basis to the path, noting where nonzeros are, relative to the QTM
      call QTM_Interpolation_Weights ( QTM_Tree, t_sv%zet_basis, s, z_ref, &
        & eta_zQT(:size(z_ref)) )
      do_calc_t = .false.
      do is = 1, n_path          ! Path length
        iz = max(s(is)%h_ind,1)  ! h_ind == zero for reflection below H_GLGrid
        do ic = 1, eta_zQT(is)%n ! Number of nonzero coefficients
          !         path Z      temperature Z        QTM serial #
          do_calc_t(iz,eta_zQT(is)%v(ic)%jz,eta_zQT(is)%v(ic)%j) = &
            & .true. ! There are only nonzero coefficients in eta_zQT
        end do
      end do

    end subroutine Temperature_Derivatives

  end subroutine More_Metrics_3D

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module More_Metrics_3D_m

! $Log$
! Revision 2.6  2017/08/28 20:28:08  livesey
! Changed the n,nf,np,nz elements to j,jf,...
!
! Revision 2.5  2017/03/11 00:54:13  vsnyder
! Pass correct-size result array to Comp_Sps_Path_Sparse_m
!
! Revision 2.4  2016/11/23 00:12:28  vsnyder
! Use types from Indexed_Values_m.
!
! Revision 2.3  2016/11/17 01:47:40  vsnyder
! Use the Grids_t structure for temperature instead of pieces of it
!
! Revision 2.2  2016/11/03 19:11:47  vsnyder
! Inching toward 3D forward model
!
! Revision 2.1  2016/10/25 18:24:29  vsnyder
! Initial commit -- still a lot to do
!
