! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Phi_Refractive_Correction_m

  implicit NONE

  private
  public :: Phi_Refractive_Correction, Phi_Refractive_Correction_Up

  interface Phi_Refractive_Correction
    module procedure Phi_Refractive_Correction_GL, Phi_Refractive_Correction_T
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! -------------------------------  Phi_Refractive_Correction_GL  -----
  subroutine Phi_Refractive_Correction_GL ( Tan_Pt, N_path, H_path, dHdz_gw, &
    &                                       Z_path, Phi_Corr )

  ! Compute the refractive correction for Phi using Gauss-Legendre quadrature,
  ! in zeta coordinates.

  !{ Write Equation A.6 from Appendix A of [1] as
  !
  ! \begin{equation}
  ! \Delta\phi_i = \int_{H_i}^{H_{i+1}} \frac{\text{d}H}H \left[
  !  \frac{\mathcal{N}_t H_t}{\sqrt{\mathcal{N}^2 H^2 - \mathcal{N}_t^2 H_t^2}}
  !  - \frac{H_t}{\sqrt{H^2-H_t^2}} \right ]\,.
  ! \end{equation}
  !
  ! Substitute $a = \frac{H_t}H$,
  ! $\frac{\text{d}H}{\text{d}a} = -\frac{H_t}{a^2}$, and
  ! $b = \frac{\mathcal{N}_t}{\mathcal{N}}$, giving
  !
  ! \begin{equation}
  ! \Delta\phi_i = - \int_{a_i}^{a_{i+1}} \text{d}a \left[
  ! \frac{b}{\sqrt{1-a^2b^2}} - \frac1{\sqrt{1-a^2}} \right]
  ! \end{equation}
  !
  ! Assuming that $\mathcal{N}$ is roughly constant in the range $(a_i,a_{i+1})$,
  ! and therefore that $b$ is roughly constant in that range, the integral can be
  ! evaluated analytically, giving
  !
  ! \begin{equation}\begin{split}
  ! \Delta\phi_i \approx\,& \left. - \left[ \sin^{-1} ab - \sin^{-1} a \right]
  ! \right|_{a=a_i}^{a_{i+1}} =
  ! \left. \sin^{-1} \frac{H_t}H - \sin^{-1} \frac{\mathcal{N}_tH_t}{\mathcal{N}H}
  !  \right|_{H=H_i}^{H_{i+1}} \\
  !  =\,&
  ! \left. \chi_{\text{eq}} - \chi_{\text{eq}}^{\text{refr}} \right|_{H=H_i}^{H_{i+1}}\,,
  ! \end{split}\end{equation}
  !
  ! Unfortunately, the arcsine method doesn't work; it has serious problems
  ! near the tangent point.  So 3-point Gauss quadrature in zeta coordinates
  ! is used instead to compute the integral on the coarse grid, and then
  ! interpolation is used to estimate it on the fine grid:
  !
  ! \begin{equation}
  ! \Delta\phi_i = \frac1{H_t} \int_{\zeta_i}^{\zeta_{i+1}} \text{d}\zeta
  !  \frac{\text{d}H}{\text{d}\zeta}
  !   \left[ \frac{a^2b}{\sqrt{1-a^2b^2}} - \frac{a^2}{\sqrt{1-a^2}} \right]
  ! \end{equation}
  !
  ! \vskip 8pt
  ! [1] W.G. Read, Z. Shippony, and W.V. Snyder, {\bf Microwave Limb Sounder (MLS)
  ! Forward Model Algorithm Theoretical Basis Document}, Jet Propulsion Laboratory
  ! document JPL D-18130, 19 August 2004.

    use GLNP, only: GX, NG, NGP1 ! Gauss abscissae, # Gauss abscissae
    use MLSKinds, only: RP, R8

  ! Inputs
    integer, intent(in) :: Tan_pt      ! Tangent point index in N_Path etc.
    real(rp), intent(in) :: N_Path(:)  ! Indices of refraction - 1 on the
                                       ! entire path.  Only the fine-grid (GL)
                                       ! part is used.
    real(rp), intent(in) :: H_Path(:)  ! Equivalent circular earth heights on
                                       ! the entire path.  Only the fine-grid
                                       ! (GL) part is used.
    real(rp), intent(in) :: dHdz_gw(:) ! dH/dz times Gauss quadrature weights,
                                       ! on the entire path.  Only the
                                       ! fine-grid (GL) part is used.
    real(rp), intent(in) :: Z_path(:)  ! Zeta on the entire path.  Only the
                                       ! coarse-grid part is used.

  ! Output
    real(rp), intent(out) :: Phi_Corr(:) ! Phi refractive corrections on the
                                       ! entire path.  The ones on the coarse
                                       ! grid are the indefinite sums of
                                       ! quadratures, working outward from the
                                       ! tangent point, of the incremental Phi
                                       ! refractive corrections for each path
                                       ! segment.  The ones at the GL points
                                       ! are interpolated linearly.

    real(rp), parameter :: GX01(NG) = 0.5 * ( 1.0_r8 + gx ) ! GX on 0..1, for
                           ! interpolating Phi to the fine grid.

    real(rp) :: A          ! a**2 at I
    real(rp) :: AB         ! A*b = a**2 * b at I
    real(rp) :: DA, DAB    ! 1 - a, 1 - ab*b
    real(rp) :: B          ! b at I
    real(rp) :: HT         ! H at tangent point
    real(rp) :: II         ! Integrand, then integral, at I
    real(r8) :: NT         ! N at tangent point
    integer :: G           ! Gauss-point index, i-mng .. i-m
    integer :: I, I1, I2   ! Subscript, loop limits
    integer :: M           ! Direction away from tangent point, +/- 1
    integer :: MNG         ! M * (NG + 1)
    integer :: My_Tan      ! min(tan_pt,ubound(phi_corr,1))

    ! Assumes ubound(phi_corr,1) == ubound(n_path,1) == u_bound(h_path,1)
    my_tan = min(tan_pt,ubound(phi_corr,1))
    ht = h_path(my_tan)
    nt = 1.0_r8 + n_path(my_tan)

    ! Account for zero-thickness tangent layer
    do i1 = my_tan, my_tan+ngp1
      if ( i1 > ubound(phi_corr,1) ) exit
      phi_corr(i1) = 0.0
    end do

    i1 = my_tan - (ng + 1)
    i2 = 1

    do m = -1, 1, 2
      mng = m * (ng + 1)
      do i = i1, i2, mng
        ! Gauss quadrature on one panel
        ii = 0.0
        do g = i-mng+m, i-m, m
          a = ( ht / h_path(g) ) ** 2     ! a^2
          b = nt / (1.0_r8 + n_path(g))   ! b
          ab = a * b                      ! a^2 b
          da = 1.0 - a
          dab = 1.0 - ab*b
          if ( da > epsilon(da) .and. dab > epsilon(dab) ) &
            & ii = ii + dHdz_gw(g) * &
              & ( ab / sqrt(dab) - a / sqrt(da) )
        end do
        ! The factor of 0.5 is because GX are on -1..1
        ii = 0.5 * (z_path(i) - z_path(i-mng)) * ii / ht
        phi_corr(i) = phi_corr(i-mng) + ii
        ! Linear interpolation from coarse to fine grid
        phi_corr(i-mng+m:i-m:m) = phi_corr(i-mng) + ii * gx01
      end do ! i
      i1 = my_tan + ng + ngp1 + 1
      i2 = size(n_path)
    end do ! m

  end subroutine Phi_Refractive_Correction_GL

  ! --------------------------------  Phi_Refractive_Correction_T  -----
  subroutine Phi_Refractive_Correction_T ( Tan_Pt, N_path, H_path, Phi_Corr )

  ! Compute the refractive correction for Phi using trapezoidal quadrature,
  ! in H coordinates.

    use GLNP, only: NGP1
    use MLSKinds, only: RP, R8

  ! Inputs
    integer, intent(in) :: Tan_pt      ! Tangent point index in N_Path etc.
    real(rp), intent(in) :: N_Path(:)  ! Indices of refraction - 1 on the
                                       ! entire path.
    real(rp), intent(in) :: H_Path(:)  ! Equivalent circular earth heights on
                                       ! the entire path.

  ! Output
    real(rp), intent(out) :: Phi_Corr(:) ! Phi refractive corrections on the
                                       ! entire path.  These are indefinite
                                       ! sums of quadratures, working outward
                                       ! from the tangent point, of the
                                       ! incremental Phi refractive corrections
                                       ! for each path segment.

    real(rp) :: A, AP      ! a I, I-M
    real(rp) :: B          ! b at I
    real(rp) :: DA, DAB    ! 1 - a**2, 1 - (a*b)**2
    real(rp) :: HT         ! H at tangent point
    real(rp) :: II, IP     ! Integrand at I, I-M
    real(r8) :: NT         ! N at tangent point
    integer :: I, I1, I2   ! Subscript, loop limits
    integer :: M           ! Direction away from tangent point, +/- 1
    integer :: My_Tan      ! min(tan_pt,size(h_path))

    ! Assumes ubound(phi_corr,1) == ubound(n_path,1) == u_bound(h_path,1)
    my_tan = min(tan_pt,ubound(phi_corr,1))
    ht = h_path(my_tan)
    nt = 1.0_r8 + n_path(my_tan)

    i1 = my_tan - 1
    i2 = 1
    ! No correction at tangent points or in zero-thickness tangent layer
    phi_corr(my_tan:min(my_tan+ngp1,ubound(phi_corr,1))) = 0.0

    do m = -1, 1, 2
      ip = 0.0
      ap = 1.0
      do i = i1, i2, m
        ! trapezoidal quadrature on one panel
        a = ( ht / h_path(i) )
        b = nt / (1.0_r8 + n_path(i))
        da = 1.0 - a**2
        dab = 1.0 - (a*b)**2
        if ( da < epsilon(da) .or. dab < epsilon(dab) ) then
          ii = 0.0
        else
          ii = 1.0 / sqrt(da) - b / sqrt(dab)
        end if
        phi_corr(i) = phi_corr(i-m) + 0.5 * ( a - ap ) * ( ip + ii )
        ! Get ready for next panel
        ap = a
        ip = ii
      end do ! i
      i1 = my_tan + ngp1 + 1
      i2 = size(n_path)
    end do ! m

  end subroutine Phi_Refractive_Correction_T

  ! -------------------------------  Phi_Refractive_Correction_Up  -----
  subroutine Phi_Refractive_Correction_Up ( N_path, H_path, Phi_Corr )

  ! Compute the refractive correction for Phi using trapezoidal quadrature,
  ! in H coordinates, from the tangent point upward -- not from the zero-
  ! thickness tangent layer outward, as the other two subroutines do.

    use MLSKinds, only: RP, R8

  ! Inputs
    real(rp), intent(in) :: N_Path(:)  ! Indices of refraction - 1.
    real(rp), intent(in) :: H_Path(:)  ! Equivalent circular earth heights.

  ! Output
    real(rp), intent(out) :: Phi_Corr(:) ! Phi refractive corrections.

    real(rp) :: A, AP      ! a I, I-M
    real(rp) :: B          ! b at I
    real(rp) :: DA, DAB    ! 1 - a**2, 1 - (a*b)**2
    real(rp) :: HT         ! H at tangent point
    real(rp) :: II, IP     ! Integrand at I, I-M
    real(r8) :: NT         ! N at tangent point
    integer :: I           ! Subscript

    if ( ubound(h_path,1) < 1 ) return
    ht = h_path(1)
    nt = 1.0_r8 + n_path(1)

    phi_corr(1) = 0.0

    ip = 0.0
    ap = 1.0
    do i = 2, ubound(h_path,1)
      ! trapezoidal quadrature on one panel
      a = ( ht / h_path(i) )
      b = nt / (1.0_r8 + n_path(i))
      da = 1.0 - a**2
      dab = 1.0 - (a*b)**2
      if ( da < epsilon(da) .or. dab < epsilon(dab) ) then
        ii = 0.0
      else
        ii = 1.0 / sqrt(da) - b / sqrt(dab)
      end if
      phi_corr(i) = phi_corr(i-1) + 0.5 * ( a - ap ) * ( ip + ii )
      ! Get ready for next panel
      ap = a
      ip = ii
    end do ! i

  end subroutine Phi_Refractive_Correction_Up

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Phi_Refractive_Correction_m

! $Log$
! Revision 2.6  2013/05/18 00:34:44  vsnyder
! Insert NG fine-grid (GL) points between tangent points, thereby
! regularizing coarse-grid spacing, and reducing significantly the need
! to use c_inds to extract coarse-grid points from the composite grid.
!
! Revision 2.5  2013/02/28 21:05:48  vsnyder
! Try to cope with short paths
!
! Revision 2.4  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.3  2006/12/13 02:32:03  vsnyder
! Drag the tangent point around instead of assuming it's the middle one
!
! Revision 2.2  2006/03/06 20:45:34  vsnyder
! Avoid simgularities caused by temperature inversions
!
! Revision 2.1  2005/12/30 01:28:01  vsnyder
! Initial commit
!
