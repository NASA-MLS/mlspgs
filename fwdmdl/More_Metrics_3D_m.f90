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
          & QTM_Tree, ddHidHidTl0, dHidTlm, T_Deriv_Flag, Z_Basis, Z_Ref, &
          ! Optional outputs
          & ddHtdHtdTl0, dHitdTlm, dHtdTl0, dHtdZt, Eta_zQT, Tan_T )

    use Generate_QTM_m, only: QTM_Tree_t
    use Geolocation_0, only: H_V_t, RG
    use MLSKinds, only: RP
    use QTM_Interpolation_Weights_3D_m, only: MaxCoeff, &
      & QTM_Interpolation_Weights, S_QTM_t, Weight_ZQ_t

    ! Inputs

    type(S_QTM_t), intent(in) :: S(:)    ! Intersections of path with constant-zeta
    integer, intent(in) :: N_Tan         ! Tangent index in S
    real(rp), intent(in) :: T_Ref(:,:)   ! temperatures at Z_Ref X hGrid
    real(rp), intent(in) :: dHidZij(:,:) ! vertical derivative at Z_Ref X hGrid

    ! Outputs

    real(rp), intent(out) :: T_Path(:)   ! computed temperatures
    real(rp), intent(out) :: dHitdZi(:)  ! derivative of height wrt zeta
                                       ! --may be useful in future computations

    ! Optional inputs

    type(QTM_tree_t), optional, intent(in) :: QTM_Tree
    real(rp), optional, intent(in) :: ddHidHidTl0(:,:,:) ! second order
    !          derivatives of height w.r.t T_Ref. This is (height, QTM,
    !          zeta_basis). Needed only if present(dHidTlm).
    real(rp), optional, intent(inout) :: dHidTlm(:,:,:) ! reference temperature
    !          derivatives. This gets adjusted so that at ref_h(1,@tan)) is
    !          0.0 for all temperature coefficients.
    !          This is height X zeta_basis X QTM
    logical, optional, intent(in) :: T_Deriv_Flag(:)  ! User's deriv. flags for
    !          Temperature.  Needed only if present(dHidTlm).
    real(rg), optional, intent(in) :: Z_Basis(:)   ! Vertical temperature basis.
    !          Needed only if present(dHidTlm).
    real(rp), optional, intent(in) :: Z_Ref(:)     ! -log pressures (zetas) for
    !          which derivatives are needed.  Only the parts from the tangent
    !          outward are used.  Needed only if present(dHidTlm).

    ! Optional outputs

    real(rp), optional, intent(out), target :: ddHtdHtdTl0(:)  ! Second
    !          order derivatives of height w.r.t T_Ref at the tangent only --
    !          used for antenna affects. Computed if present(dHidTlm).
    real(rp), optional, intent(out) :: dHitdTlm(:,:)   ! Derivative of path
    !          position wrt temperature state vector (z_basis X QTM)
    real(rp), optional, intent(out), target :: dHtdTl0(:)      ! First order
    !          derivatives of height w.r.t T_Ref at the tangent only.  Computed
    !          if present(dHidTlm).
    real(rp), optional, intent(out) :: dHtdZt      ! Height derivative wrt
    !          pressure at the tangent.  Computed if present(dHidTlm).
    type(weight_ZQ_t), optional, intent(out) :: Eta_zQT(:) ! Interpolation
    !          coefficients from QTM to path for temperature.
    real(rp), optional, intent(out) :: Tan_T ! temperature at the tangent

    integer :: I
    integer :: IZ     ! Zeta index
    integer :: NC     ! Number of S(i)%Coeff to use
    integer :: N_Path ! Size(S)

    n_path = size(s)

    ! Path quantities.  These involve only horizontal interpolation.
    ! The interpolation coefficients were computed and put into S by
    ! Metrics_3D.
    do i = 1, n_path
      nc = s(i)%n_coeff
      iz = s(i)%h_ind
      t_path(i) =  dot_product ( t_ref(iz,s(i)%ser(:nc)),   s(i)%coeff(:nc) )
      dHitdZi(i) = dot_product ( dHidZij(iz,s(i)%ser(:nc)), s(i)%coeff(:nc) )
    end do

    ! Optional tangent quantities.
    if ( n_tan <= n_path ) then
      if ( present(tan_t) ) tan_t = t_path(n_tan)
      if ( present(dHtdZt) ) dHtdZt = dHitdZi(n_tan)

      if ( present(dHidTlm) ) &
        & call Temperature_Derivatives ( s(n_tan), size(z_basis) )
    end if

  contains

    subroutine Temperature_Derivatives ( S, Z_Coeffs )
      ! This is a subroutine instead of inline so that the references to the
      ! size of Z_Basis in the dimensions of automatic variables are only
      ! attempted if Z_Basis is present.
      type(S_QTM_t), intent(in) :: S   ! Tangent point
      integer, intent(in) :: Z_Coeffs  ! size(z_basis) = size(T_ref,1)

      real(rp), pointer :: ddH2(:,:) ! Two-d (V x H) pointer for ddHtdHtdTl0
      real(rp), pointer :: dH2(:,:)  ! Two-d (V x H) pointer for dHtdTl0
      integer :: I
      integer :: IZ     ! First (vertical) subscript at tangent
      integer :: NC     ! Number of interpolating coefficients at the tangent
      type(h_v_t) :: points(z_coeffs)

      ! Adjust the 2d hydrostatic temperature derivative relative to the
      ! surface. Even though this is updated on every invocation, that is,
      ! with a new S, it works as if the original value were updated with
      ! the current S, because the interpolation represented by S is
      ! linear. Thus, the effect of cumulative updates for each new S are
      ! the same as starting from the original dHidTlm and updating with the
      ! latest S.  The algebra is horrible, but Maple has verified this.

      nc = s%n_coeff
      do i = 1, z_coeffs
        dHidTlm(:,i,:) = dHidTlm(:,i,:) - &
          & dot_product ( dHidTlm(1,i,s%ser(:nc)), s%coeff(:nc) )
      end do

      ! Get coefficients to interpolate to temperature (zeta,QTM) basis
      call QTM_Interpolation_Weights ( QTM_Tree, z_basis, points, eta_zQT )

      iz = s%h_ind
      dHtdTl0 = 0
      ddHtdHtdTl0 = 0
      dH2(1:z_coeffs,1:size(t_ref,2)) => dHtdTl0
      ddH2(1:z_coeffs,1:size(t_ref,2)) => ddHtdHtdTl0
      do i = 1, z_coeffs
        dH2(i,s%ser(:nc)) = dHidTlm(iz,i,s%ser(:nc)) * s%coeff(:nc)
        ddH2(i,s%ser(:nc)) = ddHidHidTl0(iz,i,s%ser(:nc)) * s%coeff(:nc)
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
! Revision 2.1  2016/10/25 18:24:29  vsnyder
! Initial commit -- still a lot to do
!
