! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module O2_Abs_CS_M

  implicit NONE
  private
  public :: O2_Abs_CS, D_O2_Abs_CS_dT, Get_QN_By_Frequency

  integer :: O2_in_catalog = -1   ! Index in Spectroscopy catalog of O2 line.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

! ----------------------------------------------------  O2_Abs_CS  -----
  subroutine O2_Abs_CS ( Freq, Qn, H, Slabs, Sigma_p, Pi, Sigma_m )

! Compute the complex absorption cross section.
! Modified to use Voigt with interfered lineshape

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: IP, R8, Rk => RP
    use SpectroscopyCatalog_m, only: Lines
    use Units, only: SqrtPi

    real(r8), intent(in) :: Freq              ! Observation frequency
    integer(ip), intent(in) :: Qn(:)          ! Quantum numbers
    real(rk), intent(in) :: H                 ! Magnetic field
    type(slabs_struct), intent(in) :: Slabs   ! contains, among others:
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths

    complex(rk), intent(out) :: Sigma_p       ! Output Beta values
    complex(rk), intent(out) :: Pi
    complex(rk), intent(out) :: Sigma_m

    integer(ip) :: J, No_lines
    logical, pointer :: Polarized(:)          ! Which lines to use.  Same
                                              ! as Slabs%catalog%polarized.

    real(rk) :: F_o_v0, Z, Denomm
    real(rk) :: Slabs_nonres                  ! Nonresonant absorption
    real(r8) :: V0                            ! zero field line center frequency
    real(rk) :: Y_nonres                      ! Nonresonant ratio

    complex(rk) :: Wing

    no_lines = size(slabs%catalog%lines)
    polarized => slabs%catalog%polarized
    slabs_nonres = slabs%catalog%continuum(1)
    y_nonres = slabs%catalog%continuum(3)

    sigma_p = 0.0_rk
    pi = 0.0_rk
    sigma_m = 0.0_rk
    wing = 0.0_rk

! Do magnetic calculation even if h = 0 for all of the Zeeman-split lines

    do j = 1, no_lines

      if ( .not. polarized(j) ) cycle

      v0 = lines(slabs%catalog%lines(j))%v0
      call mag_o2_abs_cs ( qn(j), freq, v0, h, slabs%x1(j), slabs%slabs1(j), &
        &                  slabs%y(j), slabs%yi(j), slabs%v0s(j), &
        &                  sigma_p, pi, sigma_m )

!{ Fill in negative frequency part of VVW lineshape for jth line
!
!  $W = \frac12 S \frac{\nu}{\nu_0} \left[
!         \frac{y-y_i z}D +
!       i \left(\frac{z + \frac{y}{\nu_{0_s}} \left( \frac{y}{x_1} - \nu y_i \right)}{D} -
!               \frac2{x_1 \nu_{0_s}} \right) \right]$, where
!  $z = x_1 ( \nu + \nu_{0_s} )$ and $D = \sqrt{\pi} ( z^2 + y^2 )$

      f_o_v0 = freq / v0
      z = slabs%x1(j) * (slabs%v0s(j) + freq)
      denomm = sqrtPi * (z*z + slabs%y(j)*slabs%y(j))
      wing = wing + ( 0.5_rk * slabs%slabs1(j) * f_o_v0 ) * cmplx( &
        &  (slabs%y(j) - slabs%yi(j)*z) / denomm, & ! Real part
        & (z + slabs%y(j) * (slabs%y(j)/slabs%x1(j) -   & ! Imaginary part...
        &   freq*slabs%yi(j)) / slabs%v0s(j)) / denomm -  &
        &   2.0_rk/(slabs%x1(j) * slabs%v0s(j)) &
        & )

    end do

    sigma_p = sigma_p + 0.5_rk * wing
    pi = pi + wing
    sigma_m = sigma_m + 0.5_rk * wing

! Contribution from non-Zeeman-split lines is done in get_beta_path_scalar.

! Non resonant contribution is done in get_beta_path_scalar too.

  end subroutine O2_Abs_CS

! -----------------------------------------------  D_O2_Abs_CS_dT  -----
  subroutine D_O2_Abs_CS_dT ( Freq, Qn, H, Slabs, Sigma_p, Pi, Sigma_m, &
    &                         dSigma_p_dT, dPi_dT, dSigma_m_dT )

! Compute the complex absorption cross section and its temperature derivative.
! Modified to use Voigt with interfered lineshape

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use MLSCommon, only: IP, R8, Rk => RP
    use SpectroscopyCatalog_m, only: Lines
    use Units, only: SqrtPi

    real(r8), intent(in) :: Freq              ! Observation frequency
    integer(ip), intent(in) :: Qn(:)          ! Quantum numbers
    real(rk), intent(in) :: H                 ! Magnetic field
    type(slabs_struct), intent(in) :: Slabs   ! contains, among others:
!    v0s(:)         ! pressure shifted line centers
!    x1(:)          ! Doppler width
!    y(:)           ! ratio Pressure to Doppler widths
!    yi(:)          ! Interference coefficients
!    slabs1(:)      ! strengths

    complex(rk), intent(out) :: Sigma_p       ! Output Beta values
    complex(rk), intent(out) :: Pi
    complex(rk), intent(out) :: Sigma_m

    complex(rk), intent(out) :: dSigma_p_dT   ! Output dBeta_dT values
    complex(rk), intent(out) :: dPi_dT
    complex(rk), intent(out) :: dSigma_m_dT

    integer(ip) :: J, No_lines
    logical, pointer :: Polarized(:)          ! Which lines to use.  Same
                                              ! as Slabs%catalog%polarized.

    real(rk) :: Slabs_nonres                  ! Nonresonant absorption
    real(r8) :: V0                            ! zero field line center frequency
    real(rk) :: Y_nonres                      ! Nonresonant ratio

    real(rk) ::  D,  I1,  I2,  I3,  I4,  R1,  R2,  S, S1, SIGMA, Y2, Z
    real(rk) :: dD, dI1, dI2, dI3, dI4, dR1, dR2,         dZ ! 1/p dp/dT
    real(rk) :: dv0 ! 1/v0s dv0s/dT
    real(rk), parameter :: OneOvSqpi = 1.0_rk / sqrtPi
    complex(rk) :: Wing, dWing

    no_lines = size(slabs%catalog%lines)
    polarized => slabs%catalog%polarized
    slabs_nonres = slabs%catalog%continuum(1)
    y_nonres = slabs%catalog%continuum(3)

    sigma_p = 0.0_rk
    pi = 0.0_rk
    sigma_m = 0.0_rk
    wing = 0.0_rk
    dSigma_p_dT = 0.0_rk
    dPi_dT = 0.0_rk
    dSigma_m_dT = 0.0_rk
    dWing = 0.0_rk

! Do magnetic calculation even if h = 0 for all of the Zeeman-split lines

    do j = 1, no_lines

      if ( .not. polarized(j) ) cycle

      v0 = lines(slabs%catalog%lines(j))%v0
      call d_mag_o2_abs_cs_dT ( qn(j), freq, v0, h,             &
        &  slabs%x1(j),     slabs%slabs1(j),     slabs%y(j),    &
        &  slabs%yi(j),     slabs%v0s(j),                       &
        & slabs%dx1_dT(j), slabs%dslabs1_dT(j), slabs%dy_dT(j), &
        & slabs%dyi_dT(j), slabs%dv0s_dT(j),                    &
        &  sigma_p,     pi,     sigma_m,                        &
        & dSigma_p_dT, dPi_dT, dSigma_m_dT )

!{ Fill in negative frequency part of VVW lineshape for jth line
!
!  $W = \frac12 S \frac{\nu}{\nu_0} \left[
!         \frac{y-y_i z}D +
!       i \left(\frac{z + \frac{y}{\nu_{0_s}} \left( \frac{y}{x_1} - \nu y_i \right)}{D} -
!               \frac2{x_1 \nu_{0_s}} \right) \right]$, where
!  $z = x_1 \sigma$, $\sigma = \nu + \nu_{0_s}$ and
!  $D = \sqrt{\pi} ( z^2 + y^2 )$.
!
!  Expanding the terms we have $W = R_1 - R_2 + i(I_1 + I_2 - I_3 - I_4)$
!  where $R_1 = \frac12 S \frac{\nu}{\nu_0} \frac{y}D$,
!        $R_2 = \frac12 S \frac{\nu}{\nu_0} \frac{y_i z}D$,
!        $I_1 = \frac12 S \frac{\nu}{\nu_0} \frac{z}{D}$,
!        $I_2 = \frac12 S \frac{\nu}{\nu_0} \frac{y^2}{D x_1 \nu_{0_s}}$,
!        $I_3 = \frac12 S \frac{\nu}{\nu_0} \frac{\nu y y_i}{D \nu_{0_s}}$, and
!        $I_4 = \frac12 S \frac{\nu}{\nu_0} \frac2{x_1 \nu_{0_s}}$, from which
!  $\frac1{R_1}\frac{\partial R_1}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} + \frac1y \frac{\partial y}{\partial T}
!      - \frac1D \frac{\partial D}{\partial T}$,
!  $\frac1{R_2}\frac{\partial R_2}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} + \frac1{y_i} \frac{\partial y_i}{\partial T}
!      + \frac1z \frac{\partial z}{\partial T} - \frac1D \frac{\partial D}{\partial T}$,
!  $\frac1{I_1}\frac{\partial I_1}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} + \frac1z \frac{\partial z}{\partial T}
!      - \frac1D \frac{\partial D}{\partial T}$,
!  $\frac1{I_2}\frac{\partial I_2}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} + 2\frac1y \frac{\partial y}{\partial T}
!      - \frac1D \frac{\partial D}{\partial T} - \frac1{\nu_{0_s}} \frac{\partial \nu_{0_s}}{\partial T}
!      - \frac1{x_1} \frac{\partial x_1}{\partial T}$,
!  $\frac1{I_3}\frac{\partial I_3}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} + \frac1y \frac{\partial y}{\partial T}
!      + \frac1{y_i} \frac{\partial y_i}{\partial T} - \frac1D \frac{\partial D}{\partial T}
!      - \frac1{\nu_{0_s}} \frac{\partial \nu_{0_s}}{\partial T}$,
!  $\frac1{I_4}\frac{\partial I_4}{\partial T}
!      = \frac1S \frac{\partial S}{\partial T} - \frac1{x_1} \frac{\partial x_1}{\partial T}
!      - \frac1{\nu_{0_s}} \frac{\partial \nu_{0_s}}{\partial T}$,
!  $\frac1D \frac{\partial D}{\partial T}
!    = 2 \frac{ z \frac{\partial z}{\partial T} + y \frac{\partial y}{\partial T}}
!             {z^2+y^2}
!    = \frac{2 \sqrt{\pi}}D
!      \left(z^2 \frac1z \frac{\partial z}{\partial T}
!          + y^2 \frac1y \frac{\partial y}{\partial T} \right)$ and
!  $\frac1z \frac{\partial z}{\partial T}
!    = \frac1{x_1} \frac{\partial x_1}{\partial T}
!    + \frac1{\sigma} \frac{\partial \nu_{0_s}}{\partial T}$.
!  Finally, we have
!  $\frac{\partial W}{\partial T} =
!     R_1 \left( \frac1{R_1}\frac{\partial R_1}{\partial T} \right)
!   + R_2 \left( \frac1{R_1}\frac{\partial R_2}{\partial T} \right)
!   + i \left[
!     I_1 \left( \frac1{I_1}\frac{\partial I_1}{\partial T} \right)
!   + I_2 \left( \frac1{I_2}\frac{\partial I_2}{\partial T} \right)
!   + I_3 \left( \frac1{I_3}\frac{\partial I_3}{\partial T} \right)
!   + I_4 \left( \frac1{I_4}\frac{\partial I_4}{\partial T} \right)
!     \right]$.
!  The quantities $\frac{\partial \nu_{0_s}}{\partial T}$,
!                 $\frac1{x_1} \frac{\partial x_1}{\partial T}$,
!                 $\frac1{y} \frac{\partial y}{\partial T}$,
!                 $\frac1{y_i} \frac{\partial y_i}{\partial T}$ and
!                 $\frac1{S} \frac{\partial S}{\partial T}$
!  are gotten from the {\tt slabs} structure.

      s = 0.5_rk * slabs%slabs1(j) * freq / v0
      sigma = freq + slabs%v0s(j)
      z = slabs%x1(j) * sigma
      r1 = slabs%y(j)
      y2 = r1**2
      r2 = slabs%yi(j)*z
      i1 = 1.0_rk / slabs%v0s(j)
      i4 = 2.0_rk * i1 / slabs%x1(j)
      i2 = y2 * 0.5_rk * i4
      i3 = freq * r1 * slabs%yi(j) * i1
      dv0 = slabs%dv0s_dT(j) * i1 ! 1/v0s dv0s/dT
      i1 = z
      d = oneOvSqpi / (z*z + y2)

      wing = wing + s * cmplx ( d * (r1 - r2), d * (i1 + i2 - i3) - i4 )

      dz = slabs%dx1_dT(j) + slabs%dv0s_dT(j) / sigma
      dD = (2.0_rk * sqrtPi) * ( z**2 * dz + y2 * slabs%dy_dT(j) ) * d
      s1 = slabs%dslabs1_dT(j) - dD ! 1/S dS/dT - 1/D dD/dT

      dr1 = s1 + slabs%dy_dT(j)
      dr2 = s1 + slabs%dyi_dT(j) + dz

      di1 = s1 + dz
      di2 = s1 - dv0 + 2.0_rk * slabs%dy_dT(j) - slabs%dx1_dT(j)
      di3 = s1 - dv0 + slabs%dy_dT(j) + slabs%dyi_dT(j)
      di4 = slabs%dslabs1_dT(j) - slabs%dx1_dT(j) - dv0

      dWing = dWing + s * &
        &             cmplx ( d * (r1 * dr1 - r2 * dr2), &
        &                     d * (i1 * di1 + i2 * di2 - i3 * di3) - i4 * di4 )

    end do

    sigma_p = sigma_p + 0.5_rk * wing
    pi =      pi      + wing
    sigma_m = sigma_m + 0.5_rk * wing

    dsigma_p_dT = dsigma_p_dT + 0.5_rk * dWing
    dpi_dT      = dpi_dT      + dWing
    dsigma_m_dT = dsigma_m_dT + 0.5_rk * dWing

! Contribution from non-Zeeman-split lines is done in get_beta_path_scalar.

! Non resonant contribution is done in get_beta_path_scalar too.

  end subroutine D_O2_Abs_CS_dT

! ------------------------------------------------  Mag_O2_Abs_CS  -----
  subroutine Mag_O2_Abs_CS ( n, nu, v0, h, x1, s, w, y, v0s, &
    &                        sigma_p, pi, sigma_m )

! Compute the frequency dependent absorption cross section for magnetic o2.

! Other notes:
! Document refers to "MLS Spectroscopic Data Base" W. G. Read, Version 1.0
! October 19, 1990.

    use MLSCommon, only: IP, R8, Rk => RP
    use Physics, only: Bohr, G_e
    use SLabs_SW_M, only: Simple_Voigt

    integer, intent(in) :: N    ! rotational quantum number, sign indicates delta J
    real(r8), intent(in) :: Nu  ! transmission frequency in MHz
    real(r8), intent(in) :: V0  ! zero magnetic field line position
    real(rk), intent(in) :: H   ! magnetic field in Gauss
    real(r8), intent(in) :: X1  ! Doppler width factor sqrt(ln 2) / D_width
    real(r8), intent(in) :: S   ! strength factor slabs1 from slabs routine
    real(r8), intent(in) :: W   ! collision to doppler width ratio in the
                                ! document corrected for temperature and pressure
    real(r8), intent(in) :: Y   ! interference coefficient in the document
                                ! corrected for pressure and temperature
    real(r8), intent(in) :: V0S ! Pressure-shifted zero magnetic field line position

    complex(rk), intent(inout) :: Sigma_P ! absorption coefficient at unity mixing
                                ! ratio for Delta M = +1
    complex(rk), intent(inout) :: Pi ! absorption coefficient at unity mixing
                                ! ratio for Delta M = 0
    complex(rk), intent(inout) :: Sigma_M ! absorption coefficient at unity mixing
                                ! ratio for Delta M = -1
    integer(ip) :: M

    real(rk), parameter :: Kappa = Bohr * G_e ! ~ 2.8024

    real(rk) :: Del_nu, Denom1, Denom2, F_o_v0, KappaH, Nu_offst
    real(rk) :: U, V, WW, Z, Zr, Zi
    real(rk) :: X

    real(rk) :: Xi ! This is actually 0.5 \xi \frac{\nu}{\nu_0} from the
                   ! ATBD, not just \xi.

! Compute the absorption coefficient at unity mixing ratio for N transition

    if ( n == 0 ) return

    f_o_v0 = nu / v0
    del_nu = v0s - nu
    kappaH = kappa * h ! H is magnetic field in Gauss here, not Planck's constant
    x = w / ( x1 * nu )
    ww = w

    if ( n == -1 ) then

!     n = -1 => denom1 = =  n * (n - 1) == 2,
!               denom2 = -n * (4 * n * n - 1) == 3

!     denom1 =  2 ! n * (n - 1)
!     denom2 =  3 ! -n * (4 * n * n - 1)

!     m = n

! sigma_p transition

! m == n == -1 => ((m+1)*(2-n)-1) == -1
!     nu_offst = x1*(del_nu + kappa*((m+1)*(2-n)-1)*h / denom1)
      nu_offst = x1 * ( del_nu - 0.5_rk*kappaH )
      call simple_voigt ( nu_offst, ww, u, v )
! m == n == -1 and denom2 == 3 => xi == 0.5
!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
!     xi = 0.5_rk
!     z  = 0.5_rk * s * xi * f_o_v0
      z  = 0.25_rk * s * f_o_v0
      zr = z * (u - y*v)
      zi = z * (v + u * (x + y))
      sigma_p = sigma_p + cmplx(zr, zi)

!     m = -n

! sigma_m transition

! m == -n == 1 => ((m-1)*(2-n)+1) == +1
!     nu_offst = x1*(del_nu + kappa*((m-1)*(2-n)+1)*h / denom1)
      nu_offst = x1 * (del_nu + 0.5_rk*kappaH )
      call simple_voigt ( nu_offst, ww, u, v )
! m == -n == 1 and denom2 == 3 => xi == 0.5
!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
!     xi = 0.5_rk
!     z  =  0.5_rk * s * xi * f_o_v0
!     z  =  0.25_rk * s * f_o_v0
      zr = z * (u - y*v)
      zi = z * (v + u * (x + y))
      sigma_m = sigma_m + cmplx(zr, zi)

!     m = n + 1

! pi transition

! m == n + 1 == 0 => kappa*m*(2-n)*h == 0
!     nu_offst = x1*(del_nu + kappa*m*(2-n)*h / denom1)
      nu_offst = x1 * del_nu
      call simple_voigt ( nu_offst, ww, u, v )
! m == n + 1 == 0 and denom2 == 3 => xi == 1.0
!     xi = 3.0_rk * (n * n - m * m) / denom2
!     z  = 0.5_rk * s * xi * f_o_v0
!     z  = 0.5_rk * s * f_o_v0
      z = 2.0_rk * z
      zr = z * (u - y*v)
      zi = z * (v + u * (x + y))
      pi = pi + cmplx(zr, zi)

    else if ( n > 0 ) then

! Delta J = +1

      denom1 = n * (n + 1)
      denom2 = (n + 1) * (2 * n +1) * (2 * n + 3)

      do m = -n , n

! pi transition

        nu_offst = x1*(del_nu + kappaH*m*(1-n) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * ((n + 1) * (n + 1) - m * m) / denom2
        xi = 1.5_rk * ((n + 1) * (n + 1) - m * m) / denom2 * f_o_v0
        z  = s * xi
        zr = z * (u - y * v)
        zi = z * (v + u  * (w / (x1*nu) + y))
        pi = pi + cmplx(zr, zi)

! sigma_p transition

        nu_offst = x1*(del_nu + kappaH*(m*(1-n)-n) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * (n + m + 1) * (n + m + 2) / (4.0_rk * denom2)
        xi = 0.375_rk * (n + m + 1) * (n + m + 2) / denom2 * f_o_v0
        z  = s * xi
        zr = z  * (u - y * v)
        zi = z * (v + u  * (w/(x1 * nu) + y))
        sigma_p = sigma_p + cmplx(zr, zi)

! sigma_m transition

        nu_offst = x1*(del_nu + kappaH*(m*(1-n)+n) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * (n - m + 1) * (n - m + 2) / (4.0_rk * denom2)
        xi = 0.375_rk * (n - m + 1) * (n - m + 2) / denom2 * f_o_v0
        z  = s * xi
        zr = z * (u - y * v)
        zi = z * (v + u  * (x + y))
        sigma_m = sigma_m + cmplx(zr, zi)

      end do

    else ! n < -1

! Delta J = -1 n is a negative number

      denom1 =  n * (n - 1)
      denom2 = -n * (4 * n * n - 1)

      do m = n+2 , -(n+2)

! pi transition

        nu_offst = x1*(del_nu + kappaH*m*(2-n) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * (n * n - m * m) / denom2
        xi = 1.5_rk * (n * n - m * m) / denom2 * f_o_v0
        z  = s * xi
        zr = z * (u - y * v)
        zi = z * (v + u  * (x + y))
        pi = pi + cmplx(zr, zi)

! sigma_p transition

        nu_offst = x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
        xi = 0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0
        z  = s * xi
        zr = z * (u - y*v)
        zi = z * (v + u  * (x + y))
        sigma_p = sigma_p + cmplx(zr, zi)

! sigma_m transition

        nu_offst = x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1)
        call simple_voigt ( nu_offst, ww, u, v )
!       xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
        xi = 0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0
        z  = s * xi
        zr = z * (u - y * v)
        zi = z * (v + u  * (x + y))
        sigma_m = sigma_m + cmplx(zr, zi)

      end do

! Finish off special transition cases for pi and sigmas
! m = -n

      m = n

! sigma_p transition

      nu_offst = x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      xi = 0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      sigma_p = sigma_p + cmplx(zr, zi)

! m = n

      m = -n

! sigma_m transition

      nu_offst = x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      xi = 0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      sigma_m = sigma_m + cmplx(zr, zi)

! m = -(n-1)

      m = n + 1

! pi transition

      nu_offst = x1*(del_nu + kappaH*m*(2-n) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (n * n - m * m) / denom2
      xi = 1.5_rk * (n * n - m * m) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      pi = pi + cmplx(zr, zi)

! sigma_p transition

      nu_offst = x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      xi = 0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      sigma_p = sigma_p + cmplx(zr, zi)

! m = n - 1

      m = -(n + 1)

! pi transition

      nu_offst = x1*(del_nu + kappaH*m*(2-n) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (n * n - m * m) / denom2
      xi = 1.5_rk * (n * n - m * m) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      pi = pi + cmplx(zr, zi)

! sigma_m transition

      nu_offst = x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1)
      call simple_voigt ( nu_offst, ww, u, v )
!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      xi = 0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0
      z  = s * xi
      zr = z * (u - y*v)
      zi = z * (v + u  * (x + y))
      sigma_m = sigma_m + cmplx(zr, zi)

    end if

  end subroutine Mag_O2_Abs_CS

! -------------------------------------------  d_Mag_O2_Abs_CS_dT  -----
  subroutine d_Mag_O2_Abs_CS_dT ( n, nu, v0, h, x1,  s,  w,  y,  v0s, &
    &                                          dx1, ds, dw, dy, dv0s, &
    &                              sigma_p,     pi,     sigma_m,      &
    &                             dSigma_p_dT, dPi_dT, dSigma_m_dT )

! Compute the frequency dependent absorption cross section for magnetic o2
! and its temperature derivative.

! Other notes:
! Document refers to "MLS Spectroscopic Data Base" W. G. Read, Version 1.0
! October 19, 1990.

    use MLSCommon, only: IP, R8, Rk => RP
    use Physics, only: Bohr, G_e
    use SLabs_SW_M, only: D_Simple_Voigt

    integer, intent(in) :: N     ! rotational quantum number, sign indicates delta J
    real(r8), intent(in) :: Nu   ! transmission frequency in MHz
    real(r8), intent(in) :: V0   ! zero magnetic field line position
    real(rk), intent(in) :: H    ! magnetic field in Gauss
    real(r8), intent(in) :: X1   ! Doppler width factor sqrt(ln 2) / D_width
    real(r8), intent(in) :: S    ! strength factor slabs1 from slabs routine
    real(r8), intent(in) :: W    ! collision to doppler width ratio in the
                                 ! document corrected for temperature and pressure
    real(r8), intent(in) :: Y    ! interference coefficient in the document
                                 ! corrected for pressure and temperature
    real(r8), intent(in) :: V0S  ! Pressure-shifted zero magnetic field line position

    real(r8), intent(in) :: dx1  ! 1/x1 dx1/dT
    real(r8), intent(in) :: ds   ! 1/S  dS/dT
    real(r8), intent(in) :: dw   ! 1/w  dw/dT
    real(r8), intent(in) :: dy   ! 1/y  dy/dT
    real(r8), intent(in) :: dv0s ! dv0s/dT, not 1/dv0s dv0s/dT

    complex(rk), intent(inout) :: Sigma_P ! absorption coefficient at unity mixing
                                 ! ratio for Delta M = +1
    complex(rk), intent(inout) :: Pi ! absorption coefficient at unity mixing
                                 ! ratio for Delta M = 0
    complex(rk), intent(inout) :: Sigma_M ! absorption coefficient at unity mixing
                                 ! ratio for Delta M = -1

    complex(rk), intent(inout) :: dSigma_P_dT ! d sigma_p / dT
    complex(rk), intent(inout) :: dPi_dT      ! d pi      / dT
    complex(rk), intent(inout) :: dSigma_M_dT ! d sigma_m /dT
    integer(ip) :: M

    real(rk), parameter :: Kappa = Bohr * G_e ! ~ 2.8024

    real(rk) :: Del_nu, Denom1, Denom2, F_o_v0, KappaH! , Nu_offst
    real(rk) :: dNu_offst, dDel_nu
    real(rk) :: X, dX, WW, WdW

! Compute the absorption coefficient at unity mixing ratio for N transition

    if ( n == 0 ) return

    f_o_v0 = nu / v0
    del_nu = v0s - nu
    dDel_nu = dv0s
    dNu_offst = x1 * ( dx1 * del_nu + dDel_nu ) ! dx1 = 1/x1 d x1 / dT
    kappaH = kappa * h
    x = w / ( x1 * nu )
    dX = x * ( dw - dx1 ) ! dw is 1/w dw/dT, dx1 is 1/dx1 dx1/dT
    ww = w
    wdw = w * dw

    if ( n == -1 ) then

!     n = -1 => denom1 = =  n * (n - 1) == 2,
!               denom2 = -n * (4 * n * n - 1) == 3

!     denom1 =  2 ! n * (n - 1)
!     denom2 =  3 ! -n * (4 * n * n - 1)

!     m = n

! sigma_p transition

! m == n == -1 => ((m+1)*(2-n)-1) == -1
!     nu_offst = x1*(del_nu + kappa*((m+1)*(2-n)-1)*h / denom1)
! m == n == -1 and denom2 == 3 => xi == 0.5
!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      call absorption ( x1 * ( del_nu - 0.5_rk*kappaH ), &
        &               0.25_rk * f_o_v0, sigma_p, dSigma_p_dT )

!     m = -n

! sigma_m transition

! m == -n == 1 => ((m-1)*(2-n)+1) == +1
!     nu_offst = x1*(del_nu + kappa*((m-1)*(2-n)+1)*h / denom1)
! m == -n == 1 and denom2 == 3 => xi == 0.5
!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      call absorption ( x1 * (del_nu + 0.5_rk*kappaH), &
        &               0.25_rk * f_o_v0, sigma_m, dSigma_m_dT )

!     m = n + 1

! pi transition

! m == n + 1 == 0 => kappa*m*(2-n)*h == 0
!     nu_offst = x1*(del_nu + kappa*m*(2-n)*h / denom1)
! m == n + 1 == 0 and denom2 == 3 => xi == 1.0
!     xi = 3.0_rk * (n * n - m * m) / denom2
      call absorption ( x1 * del_nu, 0.5_rk * f_o_v0, pi, dPi_dT )

    else if ( n > 0 ) then

! Delta J = +1

      denom1 = n * (n + 1)
      denom2 = (n + 1) * (2 * n +1) * (2 * n + 3)

      do m = -n , n

! pi transition

!       xi = 3.0_rk * ((n + 1) * (n + 1) - m * m) / denom2
        call absorption ( x1*(del_nu + kappaH*m*(1-n) / denom1), &
          &               1.5_rk * ((n + 1) * (n + 1) - m * m) / denom2 * f_o_v0, &
          &               pi, dPi_dT )

! sigma_p transition

!       xi = 3.0_rk * (n + m + 1) * (n + m + 2) / (4.0_rk * denom2)
        call absorption ( x1*(del_nu + kappaH*(m*(1-n)-n) / denom1), &
          &               0.375_rk * (n + m + 1) * (n + m + 2) / denom2 * f_o_v0, &
          &               sigma_p, dSigma_p_dT )

! sigma_m transition

!       xi = 3.0_rk * (n - m + 1) * (n - m + 2) / (4.0_rk * denom2)
        call absorption ( x1*(del_nu + kappaH*(m*(1-n)+n) / denom1), &
          &               0.375_rk * (n - m + 1) * (n - m + 2) / denom2 * f_o_v0, &
          &               sigma_m, dSigma_m_dT )

      end do

    else ! n < -1

! Delta J = -1 n is a negative number

      denom1 =  n * (n - 1)
      denom2 = -n * (4 * n * n - 1)

      do m = n+2 , -(n+2)

! pi transition

!       xi = 3.0_rk * (n * n - m * m) / denom2
        call absorption ( x1*(del_nu + kappaH*m*(2-n) / denom1), &
          &               1.5_rk * (n * n - m * m) / denom2 * f_o_v0, pi, dPi_dT )

! sigma_p transition

!       xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
        call absorption ( x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1), &
          &               0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0, &
          &               sigma_p, dSigma_p_dT )

! sigma_m transition

!       xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
        call absorption ( x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1), &
          &               0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0, &
          &               sigma_m, dSigma_m_dT )

      end do

! Finish off special transition cases for pi and sigmas
! m = -n

      m = n

! sigma_p transition

!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      call absorption ( x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1), &
        &               0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0, &
        &               sigma_p, dSigma_p_dT )

! m = n

      m = -n

! sigma_m transition

!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      call absorption ( x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1), &
        &               0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0, &
        &               sigma_m, dSigma_m_dT )

! m = -(n-1)

      m = n + 1

! pi transition

!     xi = 3.0_rk * (n * n - m * m) / denom2
      call absorption ( x1*(del_nu + kappaH*m*(2-n) / denom1), &
        &               1.5_rk * (n * n - m * m) / denom2 * f_o_v0, pi, dPi_dT )

! sigma_p transition

!     xi = 3.0_rk * (n + m) * (n + m + 1) / (4.0_rk * denom2)
      call absorption ( x1*(del_nu + kappaH*((m+1)*(2-n)-1) / denom1), &
        &               0.375_rk * (n + m) * (n + m + 1) / denom2 * f_o_v0, &
        &               sigma_p, dSigma_p_dT )

! m = n - 1

      m = -(n + 1)

! pi transition

!     xi = 3.0_rk * (n * n - m * m) / denom2
      call absorption ( x1*(del_nu + kappaH*m*(2-n) / denom1), &
                        1.5_rk * (n * n - m * m) / denom2 * f_o_v0, pi, dPi_dT )

! sigma_m transition

!     xi = 3.0_rk * (m - n) * (m - n - 1) / (4.0_rk * denom2)
      call absorption ( x1*(del_nu + kappaH*((m-1)*(2-n)+1) / denom1), &
        &               0.375_rk * (m - n) * (m - n - 1) / denom2 * f_o_v0, &
        &               sigma_m, dSigma_m_dT )

    end if

  contains

    subroutine Absorption ( nu_offst, xi, r, dr ) ! Update r ("result"), dr

      real(rk), intent(in) :: nu_offst ! x1*(del_nu + kappa*...*h/...)
      real(rk), intent(in) :: Xi ! This is actually 0.5 \xi \frac{\nu}{\nu_0}
        !                          from the ATBD, not just \xi.
      complex(rk), intent(inout) :: R, dR

      real(rk) :: U, dU, v, dV, Z, Zr, dZr, Zi, dZi
      ! dNu_offst, S, dS, w, dW, X, dX, X1, dX1, Y, dY by host association

!{ Compute absorption
!  $A = S \xi \left[ ( u - y_i v) + i ( v + u ( x + y_i ) ) \right]$ where
!  $x = \frac{y}{x_1\nu}$, $u + i v = w(\nu_o+i y)$, and $w(z)$ is the
!  Fadeeva function, and its derivative $\frac{\partial A}{\partial T}$. 
!  Here, $y_i$ is given by the variable {\tt y}, $y$ is given by the
!  variable {\tt ww}, and $\nu_o$ is given by the variable {\tt nu\_offst}.
!
!  $\frac{\partial x}{\partial T} =
!    x \left( \frac1y \frac{\partial y}{\partial T} -
!             \frac1{x_1} \frac{\partial x_1}{\partial T} \right)$.
!
!  $\Re \frac{\partial A}{\partial T} =
!   S \xi \left( \frac1S \frac{\partial S}{\partial T} ( u - y_i v ) +
!                \frac{\partial u}{\partial T} - y_i \left(\frac{\partial v}{\partial T}
!                - v \frac1{y_i} \frac{\partial y_i}{\partial T}\right) \right)$.
!
!  $\Im \frac{\partial A}{\partial T} =
!   S \xi \left( \frac1S \frac{\partial S}{\partial T}
!     (v + u (x + y_i) ) + \frac{\partial v}{\partial T} +
!    \frac{\partial u}{\partial T} ( x + y_i ) +
!    u \left( \frac{\partial x}{\partial T} +
!    y_i \frac1{y_i}\frac{\partial y_i}{\partial T} \right) \right)$

      call d_simple_voigt ( nu_offst, ww, dNu_offst, wdw, u, v, du, dv )
      zr = u - y*v
      dZr = ds * zr + du - y*(dv + dy*v)
      z = x + y
      zi = v + u * z
      dZi = ds * zi + dv + du * z + u * ( dx + y * dy )

      z  = s * xi
      r  =  r + z * cmplx( zr,  zi)
      dr = dr + z * cmplx(dZr, dZi)
    end subroutine Absorption

  end subroutine d_Mag_O2_Abs_CS_dT

! ------------------------------------------------------  Find_O2  -----
  subroutine Find_O2
  ! Find the O2 in the spectroscopy catalog

    use Molecules, only: l_o2
    use SpectroscopyCatalog_m, only: Catalog

    do o2_in_catalog = 1, size(catalog)
      if ( catalog(o2_in_catalog)%molecule == l_o2 ) return
    end do

    o2_in_catalog = -1 ! Not found
  end subroutine Find_O2

! ------------------------------------------  Get_QN_By_Frequency  -----
  subroutine Get_QN_By_Frequency ( V0, N )

    ! Get the quantum number for the O2 line that has a center frequency
    ! nearest to V0.  The quantum number of interest is given by the
    ! difference between the fourth and second numbers in the QN field
    ! of the line specification of the spectroscopy database.

    use MLSCommon, only: IP, R8, Rk => RP
    use SpectroscopyCatalog_M, only: Catalog, Lines

    real(r8), intent(in) :: V0
    integer(ip), intent(out) :: N

    integer(ip) :: I, J, K
    real(rk) :: Q, R

    ! Find O2 line the spectroscopy catalog
    if ( o2_in_catalog < 0 ) call find_o2
    n = -1
    ! Find the O2 line that has a center frequency nearest to V0
    j = catalog(o2_in_catalog)%lines(1)
    r = abs(v0-lines(j)%v0)
    do i = 2, size(catalog(o2_in_catalog)%lines)
      k = catalog(o2_in_catalog)%lines(i)
      q = abs(v0-lines(k)%v0)
      if ( q < r ) then
        r = q
        j = k
      end if
    end do
    if ( associated(lines(j)%qn) ) then
      if ( size(lines(j)%qn) >= 4 ) n = lines(j)%qn(4) - lines(j)%qn(2)
    end if
  end subroutine Get_QN_By_Frequency


! ------------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module O2_Abs_CS_M

! $Log$
! Revision 2.11  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.10  2004/04/02 01:00:20  vsnyder
! Inching toward analytic temperature derivatives
!
! Revision 2.9  2003/08/15 00:17:45  michael
! Removed extra f_o_v0 terms from mag_o2_abs_cs.  Now zero-field polarized and scalar are consistent.
!
! Revision 2.8  2003/08/14 02:15:03  vsnyder
! Optimize n=-1 special case
!
! Revision 2.7  2003/06/03 23:58:38  vsnyder
! Cosmetic changes
!
! Revision 2.6  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.5  2003/05/17 00:30:39  pwagner
! Ousted last bogus ref to sp_o2
!
! Revision 2.4  2003/05/16 23:52:53  livesey
! Now uses molecule indices rather than spectags
!
! Revision 2.3  2003/05/16 02:45:08  vsnyder
! Removed USE's for unreferenced symbols
!
! Revision 2.2  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.1.2.5  2003/03/05 03:29:59  vsnyder
! Add in the wing calculation
!
! Revision 2.1.2.4  2003/03/01 03:15:19  vsnyder
! Finish deleting Get_QN -- form the PUBLIC statement
!
! Revision 2.1.2.3  2003/03/01 03:13:47  vsnyder
! Delete unused procedure Get_QN so we won't need o2_DBase
!
! Revision 2.1.2.2  2003/03/01 03:11:02  vsnyder
! Use 'polarized' array for size, delete nonresonant computation
!
! Revision 2.1.2.1  2003/02/27 23:35:24  vsnyder
! Process all Zeeman-split lines, and no others
!
! Revision 2.1  2003/02/03 22:55:26  vsnyder
! Initial commit
!
