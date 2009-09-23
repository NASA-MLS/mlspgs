! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Ice_Size_Distribution

  implicit NONE

  private

  public :: MH_Coeffs, MH_Distribution, MH_Distribution_Whole
  public :: MH_Distribution_Derivs

  interface MH_Coeffs
    module procedure MH_Coeffs_Only, MH_Coeffs_Derivs
  end interface

  interface MH_Distribution
    module procedure MH_Distribution_Pre, MH_Distribution_1_Pre
    module procedure MH_Distribution_2_Pre
  end interface

  interface MH_Distribution_Derivs
    module procedure MH_Distribution_Derivs_
    module procedure MH_Distribution_Derivs_1, MH_Distribution_Derivs_2
  end interface

  real, parameter, public :: D_0 = 1.0          ! um
  real, parameter, public :: Rho_ice = 0.91     ! g/cm^3

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

contains

!{ Compute
!  \newcommand{\I}{\text{IWC}}
!   \begin{equation*}\begin{split}
!    & \I_{<100} = \min[\I, 0.252(\I/\I_0)^{0.837}] \\
!    & \I_{>100} = \I - \I_{<100} \\
!    & \alpha = -4.99\times 10^{-3} + 0.0494
!                 \log_{10} (\I_{<100}/\I_0) \\
!    & \mu = (5.2 + 0.0013 T) + (0.026 - 1.2 \times 10^{-3}T)
!             \log_{10} (\I_{>100}/\I_0) \\
!    & \sigma = (0.47 + 2.1 \times 10^{-3}T) + (0.018 - 2.1 \times 10^{-4} T)
!                \log_{10} (\I_{<100}/\I_0) \\
!    & N_1 = \frac{\I_{<100} \rho_{\text{ice}} \alpha^5}{4 \pi}
!    \text{ and }
!    N_2 = \frac6{\sqrt{2 \pi^3}}
!          \frac{\I_{>100}}
!               {D_0^3 \rho_{\text{ice}} \sigma \exp( 3 \mu + 4.5 \sigma^2)}
!   \end{split}\end{equation*}
!%
!   where $D_0 = 1 \mu$m, $\rho_{\text{ice}} = 0.91$ g/cm$^3$,
!   IWC$_0$ = 1 g/m$^3$ and $T$ is the atmospheric temperature in Celsius.
!%
!  The derivatives w.r.t. IWC when $\I \leq 0.252(\I/\I_0)^{0.837}$ are
!%
! \begin{equation*}\begin{split}
! \frac{\partial\alpha}{\partial\I} =\,& \frac{-0.0494\, \I_0}{\I\, \ln 10} =
!  -0.02145414741\, \frac{\I_0}{\I} \\
! \frac1{N_1} \frac{\partial N_1}{\partial\I} =\,& \frac1{\I} +
!  \frac5{\alpha} \frac{\partial\alpha}{\partial\I} \\
! \frac1{n_1} \frac{\partial n_1}{\partial\I} =\,&
!  \frac1{N_1} \frac{\partial N_1}{\partial\I} -
!  D \frac{\partial\alpha}{\partial\I} \\
! \frac{\partial\mu}{\partial\I} =\,&
!  \frac{\partial\sigma}{\partial\I} =\,
!  \frac{\partial n_2}{\partial\I} =\, 0\\
! \end{split}\end{equation*}
!%
!  Otherwise
!%
! \begin{equation*}\begin{split}
! \frac{\partial\I_{<100}}{\partial\I} =\,&\, 0.252 \times
!  \frac{0.837}\I
!  \left(\frac{\I}{\I_0}\right)^{0.837} = \frac{0.837}{\I}\, \I_{<100}
!  \text{ and }
!  \frac{\partial\I_{>100}}{\partial\I} = 1 -
! \frac{\partial\I_{<100}}{\partial\I} \\
! \frac{\partial\alpha}{\partial\I} =\,&
!  -\frac{0.0494}{\I\, \ln 10} \,
!  \frac{\partial\I_{<100}}{\partial\I} =
!  -\frac{0.02145414741}{\I}\, \frac{\partial\I_{<100}}{\partial\I} \\
! \frac{\partial\mu}{\partial\I} =\,&
!  \frac{-1.2\times10^{-3}\, T}{\I\, \ln 10} \,
!  \frac{\partial\I_{>100}}{\partial\I} =
!  -5.2115 \times 10^{-4}\, \frac{T}{\I}\,
!   \frac{\partial\I_{>100}}{\partial\I} \\
! \frac{\partial\sigma}{\partial\I} =\,&
!  \frac{-2.1\times10^{-4} \, T}{\I\, \ln 10} \,
!  \frac{\partial\I_{>100}}{\partial\I} =
!  -9.1202 \times 10^{-5}\, \frac{T}{\I}\,
!  \frac{\partial\I_{>100}}{\partial\I} =
!  0.175\, \frac{\partial\mu}{\partial\I}\\
! \frac1{N_1} \frac{\partial N_1}{\partial\I} =\,&
!  \frac1{\I_{<100}} \, \frac{\partial\I_{<100}}{\partial\I} +
!  \frac5{\alpha} \frac{\partial\alpha}{\partial\I} =
!  \frac{0.837}\I + \frac5{\alpha} \frac{\partial\alpha}{\partial\I} \\
! \frac1{N_2} \frac{\partial N_2}{\partial\I} =\,&
!  \frac1{\I_{>100}} \, \frac{\partial\I_{>100}}{\partial\I} -
!  3 \frac{\partial\mu}{\partial\I} -
! \left( \frac1{\sigma} + 4.5\right) \frac{\partial\sigma}{\partial\I} \\
!  =\,&
!  \frac1{\I_{>100}} \, \frac{\partial\I_{>100}}{\partial\I} -
!  \left( 3.7875 + \frac{0.175}\sigma \right) \frac{\partial\mu}{\partial\I} \\
! \frac1{n_1} \frac{\partial n_1}{\partial\I} =\,&
!  \frac1{N_1} \frac{\partial N_1}{\partial\I} -
!  D \frac{\partial\alpha}{\partial\I} \\
! \frac1{n_2}\frac{\partial n_2}{\partial\I} =\,&
!  \frac1{N_2} \frac{\partial N_2}{\partial\I} + \frac\gamma\sigma \,
!  \left( \frac{\partial\mu}{\partial\I} +
!  \gamma \frac{\partial\sigma}{\partial\I} \right) \\
! \end{split}\end{equation*}
!%
! The derivatives w.r.t. temperature are
! \begin{equation*}
! \frac{\partial\alpha}{\partial T} =\,
! \frac{\partial N_1}{\partial T} =\, 0
! \end{equation*}
!%
! For $\I \leq 0.252\, (\I/\I_0)^{0.837}$, $N_2 = 0$ and therefore the remaining
! derivatives are zero.  Otherwise
!%
! \begin{equation*}\begin{split}
! \frac{\partial\mu}{\partial T} =\,&\,
!  1.3 \times 10^{-3} - 1.2 \times 10^{-3}\, \log_{10} (\I_{>100}/\I_0) \\
! \frac{\partial\sigma}{\partial T} =\,&\,
!  2.1 \times 10^{-3} -2.1 \times 10^{-4}\, \log_{10} (\I_{>100}/\I_0) \\
! \frac1{N_2} \frac{\partial N_2}{\partial T} =\,&
! -3 \frac{\partial\mu}{\partial T}
! - \left( \frac1{\sigma} + 4.5\right) \frac{\partial\sigma}{\partial T} \\
! \frac1{n_2}\frac{\partial n_2}{\partial T} =\,&
!  \frac1{N_2} \frac{\partial N_2}{\partial T} +
!  \frac{\gamma}{\sigma} \left( \frac{\partial\mu}{\partial T} +
!  \gamma \frac{\partial\sigma}{\partial T} \right) \\
! \end{split}\end{equation*}
!
!   These coefficients and their derivatives depend only upon temperature and
!   IWC, but not on particle radius.

  pure subroutine MH_Coeffs_Only ( IWC, T, Alpha, Mu, Sigma, N_1, N_2 )

    use Constants, only: Pi, SqrtPi, Sqrt2
    use MLSKinds, only: R8

    real(r8), intent(in) :: IWC  ! Water content, g/m^3, require IWC > 0
    real(r8), intent(in) :: T    ! Temperature, K

    real(r8), intent(out) :: Alpha
    real(r8), intent(out) :: Mu
    real(r8), intent(out) :: Sigma
    real(r8), intent(out) :: N_1
    real(r8), intent(out) :: N_2

    ! C = 0 in Kelvins
    real(r8), parameter :: C_0 = 273.15

    ! Normalization parameters
    real(r8), parameter :: IWC_0 = 1.0    ! g/m^3

    ! Messy constants in N_1 and N_2
    real(r8), parameter :: KN1 = 1.0 / (4.0*pi*rho_ice)
    real(r8), parameter :: KN2 = 6.0 / (sqrt2 * sqrtPi**3 * rho_ice * d_0**3)

    real(r8) :: IWC_lt, IWC_gt ! Ice water content <100, >100, resp.
    real(r8) :: log_gt         ! Log10 ( iwc_gt/iwc_0 )
    real(r8) :: TC             ! Temperature, Celsius

    tc = t - c_0

    iwc_lt = min(iwc, 0.252 * ( iwc/iwc_0 ) ** 0.837 )
    iwc_gt = iwc - iwc_lt

    alpha = max( 0.0_r8, -4.99e-3_r8 - 0.0494_r8 * log10(iwc_lt/iwc_0) )
    n_1 = kn1 * iwc_lt * alpha**5

    if ( iwc_gt <= 0.0 ) then
      mu = 0.0
      sigma = 0.0
      n_2 = 0.0
    else
      log_gt = log10(iwc_gt/iwc_0)
      mu = 5.2 + 0.0013*tc + (0.026 - 1.2e-3*tc) * log_gt
      sigma = 0.47 + 2.1e-3*tc + (0.018 - 2.1e-4*tc) * log_gt
      n_2 = iwc_gt * kn2 / ( sigma * exp(3.0*mu + 4.5 * sigma**2) )
    end if

  end subroutine MH_Coeffs_Only

  pure subroutine MH_Coeffs_Derivs ( IWC, T, Alpha, Mu, Sigma, N_1, N_2, &
    &          dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
    &          dMu_dT, dSigma_dT, dN_2_dT )

    use Constants, only: Pi, SqrtPi, Sqrt2
    use MLSKinds, only: R8

    real(r8), intent(in) :: IWC  ! Water content, g/m^3, require IWC > 0
    real(r8), intent(in) :: T    ! Temperature, K

    real(r8), intent(out) :: Alpha
    real(r8), intent(out) :: Mu
    real(r8), intent(out) :: Sigma
    real(r8), intent(out) :: N_1
    real(r8), intent(out) :: N_2

    real(r8), intent(out):: dAlpha_dIWC, dMu_dIWC, dSigma_dIWC
    real(r8), intent(out):: dN_1_dIWC, dN_2_dIWC ! Actually 1/N_1... and 1/N_2...
    real(r8), intent(out):: dMu_dT, dSigma_dT
    real(r8), intent(out):: dN_2_dT ! Actually 1/N_2...

    ! C = 0 in Kelvins
    real(r8), parameter :: C_0 = 273.15

    ! Normalization parameters
    real(r8), parameter :: IWC_0 = 1.0    ! g/m^3

    ! Messy constants in N_1 and N_2
    real(r8), parameter :: KN1 = 1.0 / (4.0*pi*rho_ice)
    real(r8), parameter :: KN2 = 6.0 / (sqrt2 * sqrtPi**3 * rho_ice * d_0**3)

    real(r8) :: dIWC_lt_dIWC, dIWC_gt_dIWC ! Derivatives
    real(r8) :: IWC_i          ! 1.0 / IWC
    real(r8) :: IWC_lt, IWC_gt ! Ice water content <100, >100, resp.
    real(r8) :: log_gt         ! Log10 ( iwc_gt/iwc_0 )
    real(r8) :: SigFac         ! 1/sigma + 9 * sigma
    real(r8) :: TC             ! Temperature, Celsius

    tc = t - c_0

    iwc_lt = min(iwc, 0.252 * ( iwc/iwc_0 ) ** 0.837 )
    iwc_gt = iwc - iwc_lt

    alpha = max( 0.0_r8, -4.99e-3_r8 - 0.0494_r8 * log10(iwc_lt/iwc_0) )
    n_1 = kn1 * iwc_lt * alpha**5
    dAlpha_dIWC = -0.02145414741 / iwc_lt

    if ( iwc_gt <= 0.0 ) then
      mu = 0.0
      sigma = 0.0
      n_2 = 0.0
      dN_1_dIWC = 1.0 / iwc + 5.0 /alpha * dAlpha_dIWC
      dMu_dIWC = 0.0
      dSigma_dIWC = 0.0
      dN_2_dIWC = 0.0
      dMu_dT = 0.0
      dSigma_dT = 0.0
      dN_2_dT = 0.0
    else
      log_gt = log10(iwc_gt/iwc_0)
      mu = 5.2 + 0.0013*tc + (0.026 - 1.2e-3*tc) * log_gt
      sigma = 0.47 + 2.1e-3*tc + (0.018 - 2.1e-4*tc) * log_gt
      n_2 = iwc_gt * kn2 / ( sigma * exp(3.0*mu + 4.5 * sigma**2) )
      iwc_i = 1.0/iwc
      dIWC_lt_dIWC = 0.837 * iwc_i * iwc_lt
      dIWC_gt_dIWC = 1.0 - dIWC_lt_dIWC
      dAlpha_dIWC = dAlpha_dIWC * dIWC_lt_dIWC
      dN_1_dIWC = 0.837 * iwc_i + 5.0 / alpha * dAlpha_dIWC
      sigFac = 1/sigma + 9*sigma
      dMu_dIWC = (0.01129 - 5.2115e-4 * tc) / iwc_gt * dIWC_gt_dIWC
      dSigma_dIWC = (7.817e-3 - 9.120e-5 * tc) / iwc_gt * dIWC_gt_dIWC
      dN_2_dIWC = dIWC_gt_dIWC / iwc_gt - 3 * dMu_dIWC - sigFac * dSigma_dIWC
      dMu_dT = 1.3e-3 - 1.2e-3 * log_gt
      dSigma_dT = 2.1e-3 - 2.1e-4 * log_gt
      dN_2_dT = -3.0 * dMu_dT - sigFac * dSigma_dT
    end if

  end subroutine MH_Coeffs_Derivs

  pure function MH_Distribution_Pre ( R, Alpha, Mu, Sigma, N_1, N_2 ) result (NR)
    use MLSKinds, only: R8

    ! Compute the MH_Distribution using precomputed coefficients that
    ! are functions only of T and IWC.

    real(r8), intent(in) :: R    ! Particle radius, um
    real(r8), intent(in) :: Alpha
    real(r8), intent(in) :: Mu
    real(r8), intent(in) :: Sigma
    real(r8), intent(in) :: N_1
    real(r8), intent(in) :: N_2

    real(r8) :: NR               ! Number density

    real(r8) :: D                ! Particle diameter, um

    d = 2.0 * r
    nr = n_1 * d * exp(-alpha*d)
!    if ( d > 100.0 .and. n_2 > 0.0 ) &
    if ( n_2 > 0.0 ) &
      nr = nr + n_2 * exp(-0.5*((log(d/d_0)-mu)/sigma)**2) / d

  end function MH_Distribution_Pre

  pure function MH_Distribution_1_Pre ( R, Alpha, N_1 ) result (NR)
    use MLSKinds, only: R8

    !{ Compute $r^2 \times$ the first term in the MH\_Distribution, \emph{viz.}
    !  $N_1 D \exp(-\alpha D)$, using coefficients that are functions only of
    !  IWC and T, precomputed by MH\_Coeffs.

    real(r8), intent(in) :: R    ! Particle radius, um
    real(r8), intent(in) :: Alpha
    real(r8), intent(in) :: N_1

    real(r8) :: NR               ! Number density


    nr = 2.0 * r**3 * n_1 * exp(-2.0*alpha*r)

  end function MH_Distribution_1_Pre

  pure function MH_Distribution_2_Pre ( R, Mu, Sigma, N_2 ) result (NR)
    use MLSKinds, only: R8

    !{ Compute $r^2 \times$ the second term in the MH\_Distribution, \emph{viz.}
    !  $N_2 \exp \left[
    !    -\frac12 \left( \frac{\log(D/D_0)-\mu}\sigma \right)^2 \right]$
    !  using coefficients that are functions only of IWC and T, precomputed by
    !  MH\_Coeffs.

    real(r8), intent(in) :: R    ! Particle radius, um
    real(r8), intent(in) :: Mu
    real(r8), intent(in) :: Sigma
    real(r8), intent(in) :: N_2

    real(r8) :: NR               ! Number density


    nr = 0.5 * r * n_2 * exp(-0.5*((log((2.0/d_0)*r)-mu)/sigma)**2)

  end function MH_Distribution_2_Pre

  pure function MH_Distribution_Whole ( R, IWC, T ) result ( NR )
    use Constants, only: Pi, SqrtPi, Sqrt2
    use MLSKinds, only: R8

    real(r8), intent(in) :: R    ! Particle radius, um
    real(r8), intent(in) :: IWC  ! Water content, g/m^3
    real(r8), intent(in) :: T    ! Temperature, K

    real(r8) :: NR               ! Number density

    ! C = 0 in Kelvins
    real(r8), parameter :: C_0 = 273.15

    ! Normalization parameters
    real(r8), parameter :: IWC_0 = 1.0    ! g/m^3

    ! Messy constants in N_1 and N_2
    real(r8), parameter :: KN1 = 1.0 / (4.0*pi*rho_ice)
    real(r8), parameter :: KN2 = 6.0 / (sqrt2 * sqrtPi**3 * rho_ice * d_0)

    real(r8) :: Alpha, D
    real(r8) :: IWC_lt, IWC_gt ! Ice water content <100, >100, resp.
    real(r8) :: Log_lt, log_gt ! Log10 ( iwc_lt/iwc_0 ), Log10 ( iwc_gt/iwc_0 )
    real(r8) :: Mu, Sigma
    real(r8) :: TC             ! Temperature, Celsius

    tc = t - c_0

    d = 2.0 * r

    iwc_lt = min(iwc, 0.252 * ( iwc/iwc_0 ) ** 0.837 )
    iwc_gt = iwc - iwc_lt

    log_lt = log10(iwc_lt/iwc_0)
    log_gt = log10(iwc_gt/iwc_0)

    alpha = max( 0.0_r8, -4.99e-3_r8 - 0.0494_r8 * log_lt )

    ! Contribution from particles < 100 um
    nr = iwc_lt * kn1 * alpha**5 * d * exp(-alpha*d)

    if ( d > 100.0 .and. iwc_gt > 0.0 ) then
      ! Contribution from particles > 100 um
      mu = 5.2 + 0.0013*tc + (0.026 - 1.2e-3*tc) * log_gt
      sigma = 0.47 + 2.1e-3*tc + (0.018 - 2.1e-4*tc) * log_gt
      nr = nr + iwc_gt * kn2 / ( sigma * exp(3.0*mu + 4.5 * sigma**2) ) * &
        &       exp(-0.5*((log(d/d_0)-mu)/sigma)**2) / d
    end if

  end function MH_Distribution_Whole

  pure subroutine MH_Distribution_Derivs_ ( &
    & R, Alpha, Mu, Sigma, N_1, N_2, &
    & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
    & dMu_dT, dSigma_dT, dN_2_dT, &
    & NR, dNR_dIWC, dNR_dT )
    use MLSKinds, only: R8

    real(r8), intent(in) :: R, Alpha, Mu, Sigma, N_1, N_2
    real(r8), intent(in) :: dAlpha_dIWC, dMu_dIWC, dSigma_dIWC
    real(r8), intent(in) :: dN_1_dIWC, dN_2_dIWC
    real(r8), intent(in) :: dMu_dT, dSigma_dT, dN_2_dT
    real(r8), intent(out) :: NR, dNR_dIWC, dNR_dT

! Derivatives dN_1/... and dN_2/... are really 1/N_1 dN_1/... and 1/N_2 dN_2/....

    real(r8) :: NR_1, NR_2       ! Number density terms n_1, n_2

    real(r8) :: D                ! Particle diameter, um
    real(r8) :: Gamma            ! (log(d/d_0)-mu)/sigma

    d = 2.0 * r
    nr_1 = n_1 * d * exp(-alpha*d)
!    if ( d > 100.0 .and. n_2 > 0.0 ) then
    dNR_dIWC = nr_1 * ( dN_1_dIWC - d * dAlpha_dIWC )
    if ( n_2 > 0.0 ) then
      gamma = (log(d/d_0)-mu)/sigma
      nr_2 = n_2 * exp(-0.5*(gamma**2) ) / d
      nr = nr_1 + nr_2
      dNR_dIWC = dNR_dIWC + &
        & nr_2 * ( dN_2_dIWC + gamma/sigma * ( dMu_dIWC + gamma * dSigma_dIWC ) )
      dNR_dT = &
        & nr_2 * ( dN_2_dT   + gamma/sigma * ( dMu_dT + gamma   * dSigma_dT   ) )
    else
      nr = nr_1
      dNR_dT = 0.0
    end if

  end subroutine MH_Distribution_Derivs_

  pure subroutine MH_Distribution_Derivs_1 ( &
    & R, Alpha, N_1, dAlpha_dIWC, dN_1_dIWC, NR, dNR_dIWC )
    use MLSKinds, only: R8

  ! Compute only the first term of the MH number distribution,
  ! and its derivatives.

    real(r8), intent(in) :: R, Alpha, N_1, dAlpha_dIWC
    real(r8), intent(in) :: dN_1_dIWC ! Actually 1/N_1 dN_1/dIWC
    real(r8), intent(out) :: NR, dNR_dIWC

    real(r8) :: D                ! Particle diameter, um

    d = 2.0 * r
    nr = n_1 * d * exp(-alpha*d)
    dNR_dIWC = nr * ( dN_1_dIWC - d * dAlpha_dIWC )

  end subroutine MH_Distribution_Derivs_1

  pure subroutine MH_Distribution_Derivs_2 ( &
    & R, Mu, Sigma, N_2, &
    & dMu_dIWC, dSigma_dIWC, dN_2_dIWC, &
    & dMu_dT, dSigma_dT, dN_2_dT, &
    & NR, dNR_dIWC, dNR_dT )
    use MLSKinds, only: R8

  ! Compute only the second term of the MH number distribution,
  ! and its derivatives.

    real(r8), intent(in) :: R, Mu, Sigma, N_2
    real(r8), intent(in) :: dMu_dIWC, dSigma_dIWC
    real(r8), intent(in) :: dN_2_dIWC ! ! Actually 1/N_2 dN_2/dIWC
    real(r8), intent(in) :: dMu_dT, dSigma_dT, dN_2_dT
    real(r8), intent(out) :: NR, dNR_dIWC, dNR_dT

! Derivatives dN_1/... and dN_2/... are really 1/N_1 dN_1/... and 1/N_2 dN_2/....

    real(r8) :: D                ! Particle diameter, um
    real(r8) :: Gamma            ! (log(d/d_0)-mu)/sigma

    d = 2.0 * r
    gamma = (log(d/d_0)-mu)/sigma

    nr = n_2 * exp(-0.5*(gamma**2) ) / d
    dNR_dIWC = &
      & nr * ( dN_2_dIWC + gamma/sigma * ( dMu_dIWC + gamma * dSigma_dIWC ) )
    dNR_dT = &
      & nr * ( dN_2_dT +   gamma/sigma * ( dMu_dT +   gamma * dSigma_dT ) )

  end subroutine MH_Distribution_Derivs_2

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Ice_Size_Distribution

! $Log$
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
