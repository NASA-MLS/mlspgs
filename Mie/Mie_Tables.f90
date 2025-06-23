! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Mie_Tables

!{ Compute $\beta_{c\_e} = 2\pi \int_0^\infty n(r) r^2 \xi_e(r)\, \text{d}D$,
!          $\beta_{c\_s} = 2\pi \int_0^\infty n(r) r^2 \xi_s(r)\, \text{d}D$,
!  where $\xi_e(r)$ and $\xi_s(r)$ are computed by Mie\_Efficiencies, q.v.,
!  their derivatives w.r.t. temperature and IWC,
!          IWC$_{\text{total}}$ =
!            $\frac83 \pi \rho_{\text{ice}} \int_0^\infty n(r)\, r^3 \text{d}D$,
!  the integrated phase function
!  $\frac{\lambda^2}{2\pi\beta_{c\_s}} \int_0^\infty n(r) p_0(\theta,r) \text{d} r$,
!  and its derivatives w.r.t. temperature and IWC.

  use Cadre_m, only: Cadre_Reverse
  use Constants, only: Deg2Rad, Pi
  use Ice_Size_Distribution, only: D_0, MH_Coeffs, MH_Distribution, &
    & MH_Distribution_Derivs, Rho_Ice
  use Mie_Efficiencies_m, only: Mie_Efficiencies, Mie_Efficiencies_Derivs
  use MLSKinds, only: R8
  use Polyfit_m, only: Polyfit
  use Phase_m, only: Phase, Phase_Deriv
  use Physics, only: SpeedOfLight
  use P1_m, only: Coeffs, Compute_P1_Derivs
  use RefractiveIndex, only: UKISub, UKISub_dT
  use WriteHDF_m, only: WriteHDF

  implicit NONE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  ! Software from Math77
  interface
    subroutine DINT1 ( A, B, Answer, Work, IOPT )
      double precision, intent(in) :: A, B
      double precision, intent(out) :: Answer
      double precision, intent(inout) :: Work(*)
      integer, intent(inout) :: IOPT(*)
    end subroutine DINT1
    subroutine DINTA ( Answer, Work, IOPT )
      double precision, intent(out) :: Answer
      double precision, intent(inout) :: Work(*)
      integer, intent(inout) :: IOPT(*)
    end subroutine DINTA
  end interface

  ! Ice water content, values are log10 of bounds in g/m^3
  real(r8) :: IWC_min = -4.0    
  real(r8) :: IWC_max = 0.0
  real(r8) :: dIWC ! = ((iwc_max - iwc_min) / (n_iwc-1))
  integer :: N_IWC = 5

  ! Temperature
  real(r8),parameter :: C_0 = 273.15_r8 ! C = 0 in Kelvins
  real(r8) :: T_Min = -100.0    ! Celsius
  real(r8) :: T_Max = -15.0     ! Celsius
  real(r8) :: dT ! = ((t_max - t_min) / (n_t-1))
  integer :: N_T = 5

  ! Phase function
  real(r8) :: Theta_min = 10.0  ! Degrees
  real(r8) :: Theta_max = 170.0 ! Degrees
  real(r8) :: dTheta ! = ((theta_max - theta_min) / (n_theta-1))
  integer :: N_Theta = 2
  logical :: Half = .false.     ! Theta covers a half circle, 0..180 degrees

  real(r8), parameter :: C = 1.0e-3 * speedOfLight ! km/s
  integer :: N_f = 1
  real(r8) :: F                 ! Frequency, GHz
  real(r8), parameter :: PI2 = 2.0*pi

  ! Integration bounds
  real(r8) :: R_min = 1.0       ! Minimum ice particle size, um
  real(r8) :: R_max = 2000.0    ! Maximum ice particle size, um
  real(r8) :: R_max_cut         ! Cutoff for R_max so nothing underflows

  integer :: N_Cut = 100        ! Maximum order for Bessel functions

  ! For Xi calculation
  complex(r8), allocatable :: A(:), B(:)  ! Intermediates for xi calculations
  real(r8), allocatable :: Beta(:,:,:,:) ! Temperature X IWC X F X 2
    ! Final dimension: 1 = Beta(c_e), 2 = Beta(c_s)
  integer, parameter :: I_c_e = 1, I_c_s = 2
  real(r8) :: Chi_fac           ! 2 * pi / lambda
  ! For Eest, the final "2" dimension is Beta(c_e), Beta(c_s)
  real(r8), allocatable :: Eest(:,:,:,:) ! Error estimate: Temperature X IWC X F X 2
  real(r8), allocatable :: EestI(:,:)    ! Error estimate for IWC: Temperature X IWC
  real(r8), allocatable :: E_P(:,:,:,:)  ! Error est: Temperature X IWC X Theta X F
  real(r8) :: Lambda            ! For printing
  ! For MaxOrd, NFunc, the final "6" dimension is Beta(c_e), Beta(c_s),
  ! d_Beta(c_e)_dIWC, d_beta(c_s)_dIWC, d_Beta(c_e)_dT, d_beta(c_s)_dT
  ! For NFunc, the first "2" dimension is for number of functions, iflag
  integer, allocatable :: MaxOrd(:,:,:,:)  ! Temperature X IWC X F X 6
  integer, allocatable :: MaxOrdI(:,:)     ! Temperature X IWC
  integer, allocatable :: NFunc(:,:,:,:,:) ! 2 X Temperature X IWC X F X 6
  integer, allocatable :: NFuncI(:,:,:)    ! 2 X Temperature X IWC
  ! For MaxOrdP, NFuncP, the final "3" dimension is Phase, dPhase_dIWC,
  ! dPhase_dT.  For NFuncP, the first "2" dimension is for number of functions,
  ! iflag
  integer, allocatable :: MaxOrdP(:,:,:,:,:)  ! Temperature X IWC X Theta X F X 3
  integer, allocatable :: NFuncP(:,:,:,:,:,:) ! 2 X Temperature X IWC X Theta X F X 3
  integer :: Ord
  complex(r8) :: RF             ! Index of refraction
  real(r8) :: Xi_e, Xi_s        ! Mie efficiencies
  real(r8) :: Xis(2)            ! Both xis
  equivalence ( xis(1), xi_e ), ( xis(2), xi_s )

  ! For ice size distribution calculation
  real(r8) :: Alpha, Mu, Sigma, N_1, N_2

  ! For total IWC = 4/3 \pi \rho_{ice} \int_0^\infty r^3 n(r) dr
  real(r8), allocatable :: IWC_Tot(:,:) ! Temperature X IWC

  ! For diameter cutoff
  real(r8) :: Cut1, Cut2   ! Where terms of size distribution underflow
  real(r8) :: LT           ! -log(tiny(1.0_r8))

  ! For the phase function
  real(r8), allocatable :: C1(:), C2(:), W(:) ! Coefficients independent of R, theta
  real(r8), allocatable :: P1(:,:)     ! 0:n_Cut X n_Theta, Legendre function P_j^1
  real(r8), allocatable :: dP1_dTheta(:,:)    ! n_Cut X n_Theta
  real(r8), allocatable :: P(:,:,:,:)  ! Phase, Temperature X IWC X Theta X F
  logical, allocatable :: Blunder(:,:) ! For testing believability, Theta X 2
  real(r8), allocatable :: Diff(:,:)   ! For testing believability, Theta X 2
  real(r8), allocatable :: Fit(:,:)    ! To replace blunders, Theta X 2
  real(r8) :: Test_Tol = 3.0           ! # stdev's for blunder testing
  integer :: Blunder_Details = 0       ! <1 = no printing
                                       ! >0 = Outliers
                                       ! >1 = avg, avg2, stdev2, maxloc, replacements
                                       ! >2 = differences
  integer :: Blunder_Order = 5         ! Order of test for blunders

  ! For DINT
  integer :: Details = 0   ! <0 = Nothing
                           ! 0 = totals only
                           ! 1 = error estimates and numbers of function values
                           ! 2 = orders of Bessel functions
                           ! 3 = coefficients
  integer :: IOPT(20)
  ! Reverse communication, output, nfunc, nfmax
  data IOPT / 0, 6, 2, 0, 10, 0, 9, 2000000, 0, 0, 10*0 /
  integer :: NFUSE                    ! Number of integrands actually used
  equivalence ( IOPT(6), NFUSE )
  integer :: NFMAX                    ! Maximum integrands per integral
  equivalence ( IOPT(8), NFMAX )
  integer, parameter :: NF_s = 100    ! Maximum number of frequencies
  real(r8) :: F_s(nf_s)               ! Frequencies
  real(r8) :: IWC                     ! Ice water content, g/m^3
  real(r8), allocatable :: IWC_s(:)   ! log10(iwc's)
  real(r8) :: R                       ! Particle radius, um
  real(r8) :: T                       ! Temperature, Kelvins
  real(r8), allocatable :: T_s(:)     ! temperatures
  real(r8) :: Theta                   ! Phase angle, radians
  real(r8), allocatable :: Theta_s(:) ! Thetas
  real(r8) :: Work(20)
  equivalence ( Work(1), R )
  logical :: Tan_u = .false.          ! Do x = tan(u) change of variable
  character(7) :: Only = ''           ! Only do this integral
  integer :: Capture = 0              ! Capture integrands on this unit

  ! For CADRE
  integer :: Level = 2 ! Print level

  ! For derivatives
  logical :: Derivs = .false.
  logical :: Diffs = .false.  ! Print differences, for checking derivs
  real(r8) :: dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC
  real(r8) :: dMu_dT, dSigma_dT, dN_2_dT
  real(r8) :: dXi_e_dT, dXi_s_dT
  real(r8) :: dXis_dT(2)
  equivalence ( dXis_dT(1), dXi_e_dT ), ( dXis_dT(2), dXi_s_dT )
  complex(r8) :: dRF_dT
  complex(r8), allocatable :: dA_dT(:), dB_dT(:)  ! For phase calculations
  ! For the next four, final dimension = 1 for Beta(c_e), 2 for Beta(c_s),
  real(r8), allocatable :: dBeta_dT(:,:,:,:)     ! Temperature X IWC X F X 2
  real(r8), allocatable :: dBeta_dIWC(:,:,:,:)   ! Temperature X IWC X F X 2
  real(r8), allocatable :: e_dBeta_dT(:,:,:,:)   ! Temperature X IWC X F X 2
  real(r8), allocatable :: e_dBeta_dIWC(:,:,:,:) ! Temperature X IWC X F X 2

  real(r8), allocatable :: diffX_IWC(:,:,:)    ! Temperature X IWC-1 X F
  real(r8), allocatable :: diffIWC(:)          ! IWC-1

  real(r8), allocatable :: dP_dIWC(:,:,:,:)    ! Phase, Temperature X IWC X Theta X F
  real(r8), allocatable :: dP_dT(:,:,:,:)      ! Phase, Temperature X IWC X Theta X F
  real(r8), allocatable :: e_dP_dIWC(:,:,:,:)  ! Phase, Temperature X IWC X Theta X F
  real(r8), allocatable :: e_dP_dT(:,:,:,:)    ! Phase, Temperature X IWC X Theta X F

  ! For output
  
  character(len=*), parameter :: esv = "es" ! For Beta_ce, Beta_cs output
  character :: ES           ! Either e or s, from esv
  logical :: HDF = .true.   ! "Output HDF else output Fortran unformatted"
  logical :: Norm = .false. ! Report quantities / total IWC too
  character(7) :: String    ! Internal write for theta
  integer :: N              ! len_trim(adjustl(string))
  logical :: Progress = .false. ! Print something after every integral
  real(r8) :: T0, T1        ! for timing of progress
  real(r8) :: OT0, OT1, OT2 ! For timing for each input request.
  logical :: WantBeta = .true., WantIWC = .true., WantP = .true.
  logical :: Warn = .false. ! Warnings from DINT
  character(len=1023) :: File = '' ! Output to this file if not blank; see HDF

  integer :: I, I_Beta, I_IWC, I_T, I_Theta, I_F, J

  namelist / in / IWC_Min, IWC_Max, N_IWC, T_Min, T_Max, N_T, &
    &             F_s, R_Min, R_Max, &
    &             Theta_Min, Theta_Max, N_Theta, Half, N_Cut, File, HDF, &
    &             IOPT, NFMAX, WORK, Warn, Level, Tan_u, Only, Capture, &
    &             Derivs, Details, Diffs, Norm, Progress, &
    &             WantBeta, WantIWC, WantP, Test_Tol, Blunder_Details, &
    &             Blunder_Order

  ! Compute diameter cutoff parameters related to machine arithmetic
!   lt = -log(tiny(1.0_r8))
!   lt = -log(tiny(1.0)) ! Clearly wrong, but this is what it was at one time
  lt = -2.0 * log(epsilon(1.0_r8)/1000.0)

  call cpu_time ( ot0 )
  do
    f_s = -1.0 ! Sentinel indicating no more frequencies
    file = ''  ! Can't write more than one data set to a single file
    read ( *, in, end=9 )

    do n_f = 0, size(f_s)-1
      if ( f_s(n_f+1) < 0.0 ) exit
    end do

    if ( allocated(beta) ) then
      if ( any(shape(p) /= (/ n_t, n_iwc, n_theta, n_f /)) .or. size(a) /= n_cut ) then
        deallocate ( IWC_s, T_s, theta_s, beta, eest, eestI, maxOrd, maxOrdI, &
          & maxOrdP, nFunc, nfuncI, nfuncP, a, b, c1, c2, w, p1, p, iwc_tot,  &
          & e_P, blunder, diff, fit )
        if ( derivs ) &
          & deallocate ( dBeta_dT, dBeta_dIWC, e_dBeta_dT, e_dBeta_dIWC, &
            & diffX_IWC, diffIWC, da_dT, db_dT, dP1_dTheta, dP_dIWC, dP_dT, &
            & e_dP_dIWC, e_dP_dT )
      end if
    end if

    if ( .not. allocated(beta) ) then
      allocate ( IWC_s(n_iwc), T_s(n_t), theta_s(n_theta), &
        &        beta(n_t, n_iwc, n_f, 2), eest(n_t, n_iwc, n_f, 2), &
        &        eestI(n_t, n_iwc), maxOrd(n_t, n_iwc, n_f, 6), &
        &        maxOrdI(n_t, n_iwc), maxOrdP(n_t, n_iwc, n_theta, n_f, 3), &
        &        nFunc(2, n_t, n_iwc, n_f, 6), nFuncI(2, n_t, n_iwc), &
        &        nFuncP(2, n_t, n_iwc, n_theta, n_f, 3), &
        &        a(n_cut), b(n_cut), c1(2:n_cut), c2(2:n_cut), w(n_cut), &
        &        p1(0:n_cut,n_Theta), dP1_dTheta(n_Cut,n_theta), &
        &        iwc_tot(n_t, n_iwc), &
        &        p(n_t, n_iwc, n_theta, n_f), e_P(n_t, n_iwc, n_theta, n_f), &
        &        blunder(n_theta,2), diff(n_theta,2), fit(n_theta,2) )
      call coeffs ( c1, c2, w ) ! For Legendre function

      if ( derivs ) &
        & allocate ( dBeta_dT(n_t, n_iwc, n_f, 2), dBeta_dIWC(n_t, n_iwc, n_f, 2), &
            &        e_dBeta_dT(n_t, n_iwc, n_f, 2), e_dBeta_dIWC(n_t, n_iwc, n_f, 2), &
            &        diffX_IWC(n_t, n_iwc-1,n_f), diffIWC(n_iwc-1), &
            &        da_dT(n_cut), db_dT(n_cut), &
            &        dP_dIWC(n_t, n_iwc, n_theta, n_f), &
            &        dP_dT(n_t, n_iwc, n_theta, n_f), &
            &        e_dP_dIWC(n_t, n_iwc, n_theta, n_f), &
            &        e_dP_dT(n_t, n_iwc, n_theta, n_f) )
    end if

    if ( details > 2 ) print '(3a)', &
      & '   T      IWC    Alpha    mu   Sigma', &
      & '  Rmax1   Rmax2    N1          N2'

    if ( n_iwc > 1 ) then
      dIWC = ((iwc_max - iwc_min) / (n_iwc-1))
    else
      dIWC = 0.0
    end if
    if ( diffs .and. n_iwc > 1 ) then
      do i_IWC = 2, n_IWC ! Loop for IWC
        diffIWC(i_iwc-1) = 10.0 ** (iwc_min + (i_iwc - 1) * dIWC) - &
          &                10.0 ** (iwc_min + (i_iwc - 2) * dIWC)
      end do
    end if
    if ( n_t > 1 ) then
      dT = (t_max - t_min) / (n_t-1)
    else
      dT = 0.0
    end if
    if ( n_Theta > 1 ) then
      dTheta = ((theta_max - theta_min) / (n_theta-1)) * deg2Rad
    else
      dTheta = 0.0
    end if

    ! Compute coordinate arrays.
    do i_iwc = 1, n_iwc
      iwc_s(i_iwc) = iwc_min + (i_iwc - 1) * dIWC
    end do
    do i_t = 1, n_t
      t_s(i_t) = t_min + c_0 + (i_t - 1) * dT
    end do
    do i_theta = 1, n_theta
      theta_s(i_theta) = theta_min * deg2Rad + (i_theta-1) * dTheta
    end do

    ! Compute Legendre functions
    do i_theta = 1, n_theta
      call compute_p1_derivs ( &
        & c1, c2, theta_s(i_theta), p1(:,i_theta), dP1_dTheta(:,i_theta) )
    end do
    maxOrd = 0
    maxOrdP = 0

    call cpu_time ( t0 )
    do i_F = 1, n_f ! Loop for frequencies
      f = f_s(i_f)
      lambda = c / f              ! Wavelength in micrometers
      chi_fac = pi2 / lambda
      do i_T = 1, n_T ! Loop for temperatures
        t = t_s(i_t)
        if ( derivs ) then
          call ukisub_dT ( f, t, rf, dRF_dT )  ! Get index of refraction for ice
        else
          call ukisub ( f, t, rf )             ! Get index of refraction for ice
        end if

        do i_IWC = 1, n_IWC ! Loop for IWC
          iwc = 10.0 ** iwc_s(i_iwc)
          ! Get radius-independent but temperature- and IWC-dependent coefficients
          if ( derivs ) then
            call mh_coeffs ( IWC, T, alpha, mu, sigma, N_1, N_2, &
      &          dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
      &          dMu_dT, dSigma_dT, dN_2_dT )
          else
            call mh_coeffs ( IWC, T, alpha, mu, sigma, N_1, N_2 )
          end if

          cut1 = 0.5 * (log(n_1)+lt)/alpha ! r**2 * N_1 underflows at cut1
          if ( n_2 <= 0.0 ) then
            r_max_cut = min(r_max,cut1)
          else
                                           ! r**2 * N_2 underflows at cut2
            cut2 = 0.5*d_0 * exp(mu+sqrt((log(n_2)+lt))*sigma)
            r_max_cut = min(r_max,cut2)
          end if

          if ( details > 2 ) then
            if ( n_2 <= 0.0 ) then
              print '(f7.2,f8.4,f7.4,14x,f8.2,9x,2g12.5)', &
                & t-c_0, iwc, alpha,            &
                & cut1,       n_1, n_2
            else
              print '(f7.2,f8.4,3f7.4,f8.2,f9.0,2g12.5)', &
                & t-c_0, iwc, alpha, mu, sigma, &
                & cut1, cut2, n_1, n_2
            end if
          end if

          if ( (wantBeta .or. wantP) .and. only == '' ) then
            do i_beta = 1, 2        ! 1 = Beta(c_e), 2 = Beta(c_s)
              es = esv(i_beta:i_beta)
              ! Integrate the whole range at once -- let DINT find peaks
              call do_dint_beta ( r_min, r_max_cut, beta(i_T,i_IWC,i_f,i_beta), &
                & eest(i_T,i_IWC,i_f,i_beta) )
              call progress_report ( 'beta_c'//es, beta(i_T,i_IWC,i_f,i_beta), &
                & eest(i_T,i_IWC,i_f,i_beta), nfunc(:,i_T,i_IWC,i_f,i_beta) )
              ! Derivatives?
              if ( derivs ) then
                ! Integrate the whole range at once -- let DINT find peaks
                call do_dint_dBeta_dIWC ( r_min, r_max_cut, dBeta_dIWC(i_T,i_IWC,i_f,i_beta), &
                  & e_dBeta_dIWC(i_T,i_IWC,i_f,i_beta) )
                call progress_report ( 'dBeta_c'//es//'_dIWC', dBeta_dIWC(i_T,i_IWC,i_f,i_beta), &
                  & e_dBeta_dIWC(i_T,i_IWC,i_f,i_beta), nfunc(:,i_T,i_IWC,i_f,2+i_beta) )
                call do_dint_dBeta_dT ( r_min, r_max_cut, dBeta_dT(i_T,i_IWC,i_f,i_beta), &
                  & e_dBeta_dT(i_T,i_IWC,i_f,i_beta) )
                call progress_report ( 'dBeta_c'//es//'_dT', dBeta_dT(i_T,i_IWC,i_f,i_beta), &
                  & e_dBeta_dT(i_T,i_IWC,i_f,i_beta), nfunc(:,i_T,i_IWC,i_f,4+i_beta) )
              end if
            end do ! Beta c_e or Beta c_s
          end if
          if ( wantIWC ) then
            if ( i_f == 1 ) then
              call do_dint_ice ( r_min, r_max_cut, iwc_tot(i_T,i_IWC), &
                & eestI(i_T,i_IWC) )
              call progress_report ( 'IWC', iwc_tot(i_T,i_IWC), eestI(i_T,i_IWC), &
                & nfuncI(:,i_T,i_IWC) )
            end if
          end if

          if ( wantP .or. only(1:1) /= '' ) then
            if ( only(1:1) == '' .or. only == 'P' ) then
              do i_Theta = 1, n_Theta
                theta = theta_s(i_theta)
                call do_dint_phase ( r_min, r_max_cut, p(i_T,i_IWC,i_Theta,i_f), &
                  & e_p(i_T,i_IWC,i_Theta,i_f), only )
                call progress_report ( 'P', p(i_T,i_IWC,i_Theta,i_f), &
                  & e_p(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,1), sub3=i_theta )
                if ( iopt(1) > 0 ) then
                  call do_cadre_phase ( r_min, r_max_cut, p(i_T,i_IWC,i_Theta,i_f), &
                    & e_p(i_T,i_IWC,i_Theta,i_f) )
                  call progress_report ( 'P', p(i_T,i_IWC,i_Theta,i_f), &
                    & e_p(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,1), sub3=i_theta )
                end if
              end do ! i_theta
            end if
            if ( only(1:1) == '' ) then
              ! Test whether the results are a smooth function of theta.
              call blunder_test ( p, half, i_T, i_IWC, i_f, "P", 1, &
                                & blunder(:,1), diff(:,1), fit(:,1) )
              do i_Theta = 1, n_Theta
                if ( blunder(i_theta,1) ) then
                  theta = theta_s(i_theta)
                  call do_cadre_phase ( r_min, r_max_cut, p(i_T,i_IWC,i_Theta,i_f), &
                    & e_p(i_T,i_IWC,i_Theta,i_f) )
                  call progress_report ( 'P', p(i_T,i_IWC,i_Theta,i_f), &
                    & e_p(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,2), sub3=i_theta )
                end if
              end do
              do j = 2, 9
                if ( .not. any(blunder(:,1)) ) exit
                call blunder_test ( p, half, i_T, i_IWC, i_f, "P", j, &
                                  & blunder(:,2), diff(:,2), fit(:,2) )
                call blunder_check ( p, blunder, diff, fit, &
                                   & i_T, i_IWC, i_f, "P", j )
                blunder(:,1)=blunder(:,2); diff(:,1)=diff(:,2); fit(:,1)=fit(:,2)
              end do
            end if
            if ( derivs ) then
              if ( only(1:1) == '' .or. only == 'dP_dIWC' ) then
                do i_Theta = 1, n_Theta
                  theta = theta_s(i_theta)
                  call do_dint_dPhase_dIWC ( r_min, r_max_cut, dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                    & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f), only )
                  call progress_report ( 'dP_dIWC', dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                    & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,2), sub3=i_theta )
                  if ( iopt(1) > 0 ) then
                    call do_cadre_dPhase_dIWC ( r_min, r_max_cut, dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f) )
                    call progress_report ( 'dP_dIWC', dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,2), sub3=i_theta )
                  end if
                end do
              end if
              if ( only(1:1) == '' ) then
                call blunder_test ( dP_dIWC, half, i_T, i_IWC, i_f, "dP_dIWC", 1, &
                                  & blunder(:,1), diff(:,1), fit(:,1) )
                do i_Theta = 1, n_Theta
                  if ( blunder(i_theta,1) ) then
                    call do_cadre_dPhase_dIWC ( r_min, r_max_cut, dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f) )
                    call progress_report ( 'dP_dIWC', dP_dIWC(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dIWC(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,2), sub3=i_theta )
                  end if
                end do
                do j = 2, 9
                  if ( .not. any(blunder(:,1)) ) exit
                  call blunder_test ( dP_dIWC, half, i_T, i_IWC, i_f, "dP_dIWC", j, &
                                    & blunder(:,2), diff(:,2), fit(:,2) )
                  call blunder_check ( dP_dIWC, blunder, diff, fit, &
                                     & i_T, i_IWC, i_f, "dP_dIWC", j )
                  blunder(:,1)=blunder(:,2); diff(:,1)=diff(:,2); fit(:,1)=fit(:,2)
                end do
              end if
              if ( only(1:1) == '' .or. only == 'dP_dT' ) then
                do i_Theta = 1, n_Theta
                  theta = theta_s(i_theta)
                  call do_dint_dPhase_dT ( r_min, r_max_cut, dP_dT(i_T,i_IWC,i_Theta,i_f), &
                    & e_dP_dT(i_T,i_IWC,i_Theta,i_f), only )
                  call progress_report ( 'dP_dT', dP_dT(i_T,i_IWC,i_Theta,i_f), &
                    & e_dP_dT(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,3), sub3=i_theta )
                  if ( iopt(1) > 0 ) then
                    call do_cadre_dPhase_dT ( r_min, r_max_cut, dP_dT(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dT(i_T,i_IWC,i_Theta,i_f) )
                    call progress_report ( 'dP_dT', dP_dT(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dT(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,3), sub3=i_theta )
                  end if                  
                end do
              end if
              if ( only(1:1) == '' ) then
                call blunder_test ( dP_dT, half, i_T, i_IWC, i_f, "dP_dT", 1, &
                                  & blunder(:,1), diff(:,1), fit(:,1) )
                do i_Theta = 1, n_Theta
                  if ( blunder(i_theta,1) ) then
                    call do_cadre_dPhase_dT ( r_min, r_max_cut, dP_dT(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dT(i_T,i_IWC,i_Theta,i_f) )
                    call progress_report ( 'dP_dT', dP_dT(i_T,i_IWC,i_Theta,i_f), &
                      & e_dP_dT(i_T,i_IWC,i_Theta,i_f), nfuncP(:,i_T,i_IWC,i_Theta,i_f,3), sub3=i_theta )
                  end if
                end do
                do j = 2, 9
                  if ( .not. any(blunder(:,1)) ) exit
                  call blunder_test ( dP_dT, half, i_T, i_IWC, i_f, "dP_dT", j, &
                                    & blunder(:,2), diff(:,2), fit(:,2) )
                  call blunder_check ( dP_dT, blunder, diff, fit, &
                                     & i_T, i_IWC, i_f, "dP_dT", j )
                  blunder(:,1)=blunder(:,2); diff(:,1)=diff(:,2); fit(:,1)=fit(:,2)
                end do
              end if
            end if
          end if ! wantP
        end do ! IWC
      end do ! T
      call cpu_time ( ot2 )
      write ( *, 7 ) ot2 - ot0, "so far"
    end do ! F

! Report the results

    if ( details >= 0 ) call printResults

    if ( file /= '' ) then
      if ( hdf ) then
        if ( derivs ) then
          call writeHDF ( &
            & File, R_max, R_min, N_cut, IWC_s, T_s, Theta_s, F_s(:n_f), &
            & Beta, Eest, nFunc, MaxOrd, P, E_P, nFuncP, MaxOrdP, &
            & wantBeta, wantIWC, wantP, &
            & dBeta_dIWC, E_dBeta_dIWC, dBeta_dT, E_dBeta_dT, &
            & dP_dIWC, E_dP_dIWC, dP_dT, E_dP_dT )
        else
          call writeHDF ( &
            & File, R_max, R_min, N_cut, IWC_s, T_s, Theta_s, F_s(:n_f), &
            & Beta, Eest, nFunc, MaxOrd, P, E_P, nFuncP, MaxOrdP, &
            & wantBeta, wantIWC, wantP )
        end if
      else
        open ( 10, file=trim(file), form='unformatted' )
        write ( 10 ) n_f, n_iwc, n_t, n_theta, r_min, r_max, n_cut, derivs
        if ( wantIWC ) &
          & write ( 10 ) iwc_s, t_s, theta_s, f_s
        if ( wantBeta ) then
          write ( 10 ) beta(:,:,:,i_c_e), eest(:,:,:,i_c_e), nFunc(:,:,:,:,i_c_e), maxOrd(:,:,:,i_c_e)
          write ( 10 ) beta(:,:,:,i_c_s), eest(:,:,:,i_c_s), nFunc(:,:,:,:,i_c_s), maxOrd(:,:,:,i_c_s)
        end if
        if ( wantP ) &
          & write ( 10 ) p, e_p, nFuncP(:,:,:,:,:,1), maxOrdP(:,:,:,:,1)
        if ( derivs ) then
          if ( wantBeta ) then
            write ( 10 ) dBeta_dIWC(:,:,:,i_c_e),  e_dBeta_dIWC(:,:,:,i_c_e), nFunc(:,:,:,:,3), maxOrd(:,:,:,3)
            write ( 10 ) dBeta_dIWC(:,:,:,i_c_s),  e_dBeta_dIWC(:,:,:,i_c_s), nFunc(:,:,:,:,4), maxOrd(:,:,:,4)
            write ( 10 ) dBeta_dT(:,:,:,i_c_e),  e_dBeta_dT(:,:,:,i_c_e), nFunc(:,:,:,:,5), maxOrd(:,:,:,5)
            write ( 10 ) dBeta_dT(:,:,:,i_c_s),  e_dBeta_dT(:,:,:,i_c_s), nFunc(:,:,:,:,6), maxOrd(:,:,:,6)
          end if
          if ( wantP ) then
            write ( 10 ) dP_dIWC, e_dP_dIWC, nFuncP(:,:,:,:,:,2), maxOrdP(:,:,:,:,2)
            write ( 10 ) dP_dT, e_dP_dT, nFuncP(:,:,:,:,:,3), maxOrdP(:,:,:,:,3)
          end if
        end if
        close ( 10 )
      end if
    end if

    call cpu_time ( ot1 )
    write ( *, 7 ) ot1 - ot0
7   format ( 'Used ', f9.2, ' CPU seconds', :, 1x, a )
    ot0 = ot1

  end do ! input
9 continue

contains

  subroutine Blunder_Check ( To_Test, Blunder, Diff, Fit, i_T, i_IWC, i_f, &
    &                        What, Iter )
    ! Where blunder(i_theta) check whether
    ! abs(diff(i_theta,1)) < abs(diff(i_theta,2)).
    ! If so, DINT did a better job than CADRE, so move Save(i_theta) back to
    ! To_Test(i_T,i_IWC,i_theta,i_f) and announce it.
    real(r8), intent(inout) :: To_Test(:,:,:,:)
    logical, intent(in) :: Blunder(:,:) ! n_theta X 2
    real(r8), intent(in) :: Diff(:,:)   ! n_theta X 2
    real(r8), intent(in) :: Fit(:,:)    ! n_theta X 2
    integer, intent(in) :: i_T, i_IWC, i_f
    character(len=*), intent(in) :: What
    integer, intent(in) :: Iter
    integer :: I_Theta

    if ( any(blunder(:,2)) ) then
      write ( *, '(a)' ) repeat("*",72)
      do i_theta = 1, size(blunder,1)
        if ( blunder(i_theta,2) ) then ! Cadre didn't work either
          if ( .not. blunder(i_theta,1) ) then
            write ( *, 666 ) "Spurious blunder report after CADRE for ", &
            & what, i_T, i_IWC, i_theta, i_f
            cycle
          end if
          write ( *, 666, advance="no" ) "CADRE and DINT both failed for ", &
            & what, i_T, i_IWC, i_theta, i_f
666       format ( a, a, "(", i0, 3(",",i0), ")" )
          if ( iter > 2 ) then
            to_test(i_T, i_IWC, i_theta, i_f) = fit(i_theta,1)
            write ( *, '(a)' ) ", using LS with previous LS result"
          else if ( abs(diff(i_theta,1)) < abs(diff(i_theta,2)) ) then
            to_test(i_T, i_IWC, i_theta, i_f) = fit(i_theta,1)
            write ( *, '(a)' ) ", using LS with previous DINT result"
          else
            to_test(i_T, i_IWC, i_theta, i_f) = fit(i_theta,2)
            write ( *, '(a)' ) ", using LS with CADRE result"
          end if
        end if
      end do
    end if
    write ( *, '(a)' ) repeat("*",72)
  end subroutine Blunder_Check

  subroutine Blunder_Test ( To_Test, Half, i_T, i_IWC, i_f, What, Which, &
    &                       Blunder, Diff, Fit )
    ! Test whether To_Test is a smooth function of theta.  Form differences
    ! relative to To_Test, then look for consecutive relative differences
    ! greater than Test_Tol in magnitude but with opposite signs.
    real(r8), intent(in) :: To_Test(:,:,:,:)
    logical, intent(in) :: Half  ! Theta_s cover 0..180
    integer, intent(in) :: i_T, i_IWC, i_f
    character(len=*), intent(in) :: What
    integer, intent(in) :: Which ! 1 = DINT, 2 = CADRE
    logical, intent(out) :: Blunder(:)
    real(r8), intent(out) :: Diff(:)
    real(r8), intent(out) :: Fit(:)
    integer :: I_Theta, I, J, M, N, S
    integer :: Order
    real(r8), dimension(1-blunder_order:size(to_test,3)+blunder_order) :: &
              & My_Test, X, My_Fit, My_Diff
    real(r8) :: D_Theta ! Theta stepsize -- theta_s are evenly spaced
    real(r8) :: Stdev   ! of the fit

    m = 1
    n = size(to_test,3)
    s = n

    d_theta = theta_s(2) - theta_s(1) ! theta_s are evenly spaced
    x(m:n) = theta_s
    do i = 0, 1-blunder_order, -1
      x(i) = x(i+1) - d_theta
    end do
    do i = n + 1, n + blunder_order
      x(i) = x(i-1) + d_theta
    end do

    my_test(m:n) = to_test(i_T,i_IWC,:,i_f)
    if ( half ) then ! extend cyclically
      my_test(1-blunder_order:0) = to_test(i_T,i_IWC,blunder_order+1:2:-1,i_f)
      my_test(n+1:n+blunder_order) = to_test(i_T,i_IWC,n-1:n-blunder_order:-1,i_f)
    else             ! extend with constant
      my_test(1-blunder_order:0) = to_test(i_T,i_IWC,1,i_f)
      my_test(n+1:n+blunder_order) = to_test(i_T,i_IWC,n,i_f)
    end if
    m = 1-blunder_order
    n = n+blunder_order

    ! Do least-squares fit to polynomial
    call polyfit ( x, my_test, blunder_order, order, stdev, my_fit, my_diff )
    fit(1:s) = my_fit(1:s)
    diff(1:s) = my_diff(1:s)

    ! Compute blunders
    blunder = abs(diff) > test_tol * stdev

    if ( blunder_details > 0 ) then
      write ( *, '(a,1pg15.6,a,i0)' ) 'Standard deviation of residual =', stdev, &
        & ', Order of fit = ', order
      if ( blunder_details > 1 ) then
        write ( *, '(a)' ) 'Fitted result'
        do i = m, n, 5
          write ( *, '(i4,": ",1p,5g13.6)' ) i, (my_fit(j), j=i,min(i+4,n))
        end do
      end if
    end if

    if ( any(blunder(1:s)) ) then
      write ( *, '(a)' ) repeat("*",72)
      do i_theta = 1, s
        if ( blunder(i_theta) ) then
          write ( *, 666 ) "Blunder ", which , what, i_T, i_IWC, i_theta-1, i_theta+1, i_f, &
            &              my_test(i_theta-1:i_theta+1) !, cut(i_theta)
666       format ( a, i1, ": ", a, "(", i0, 2(",", i0), ":", i0, ",", i0, ") = ", &
                 & 1pg13.6, 2(", ", 1pg13.6): ", Cut = ", g13.6 : &
                 & ", Avg, Avg2, Stdev2 = ", 1p,3g13.6 )
          if ( blunder_details > 0 ) &
            & write ( *, 666 ) "Diff     ", which, what, i_T, i_IWC, i_theta-1, i_theta+1, i_f, &
              &              my_diff(i_theta-1:i_theta+1)
        end if
      end do ! i_theta
      write ( *, '(a)' ) repeat("*",72)
    end if
  end subroutine Blunder_Test

  subroutine Do_dint_beta ( R_min, R_max, Answer, Error )
  !{ Do an integration to get beta using both terms of the number
  !  distribution: $\int_0^\infty r^2 \xi(r) \left\{
  !    N_1 2 r \exp(-\alpha 2 r) +
  !    \frac{N_2}{2 r} \exp \left[ -\frac12
  !      \left( \frac{\log(2 r/D_0) - \mu}{\sigma} \right)^2 \right]
  !    \right\} \text{d} r$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    ! Start the quadrature
    call dint1 ( r_min, r_max, answer, work, iopt )
    ! Evaluate the integrand
    do
      call dinta ( answer, work, iopt )
      if ( iopt(1) /= 0 ) exit
      call mie_efficiencies ( rf, chi_fac * r, xi_e, xi_s, a, b, ord )
      answer = xis(i_beta) * r * r * mh_distribution(r,alpha,mu,sigma,n_1,n_2)
      maxOrd(i_T,i_IWC,i_f,i_beta) = max(maxOrd(i_T,i_IWC,i_f,i_beta),ord)
    end do
    ! Finished the quadrature.  Integral is in Answer.
    answer = pi2 * answer
    error = pi2 * work(1)
    nFunc(1,i_T,i_IWC,i_f,i_beta) = nfuse
    nFunc(2,i_T,i_IWC,i_f,i_beta) = -iopt(1)
    if ( iopt(1) > 0 ) call errorReport ( "Beta", iopt(1), answer, error )
  end subroutine Do_dint_beta

  subroutine Do_dint_dBeta_dIWC ( R_min, R_max, Answer, Error )
  !{ Do an integration to get the derivative of beta w.r.t.\ IWC:\\
  !  $\int_0^\infty r^2 \left( \xi(r) \frac{\partial n(r)}{\partial \text{IWC}} +
  !    \frac{\partial \xi(r)}{\partial \text{IWC}} n(r) \right) \text{d} r$,
  !  where $n(r) = N_1(r) + N_2(r)$.
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8) :: NR, dNR_dIWC, dNR_dT
    ! Start the quadrature
    call dint1 ( r_min, r_max, answer, work, iopt )
    ! Evaluate the integrand
    do
      call dinta ( answer, work, iopt )
      if ( iopt(1) /= 0 ) exit
      call mie_efficiencies ( rf, chi_fac * r, xi_e, xi_s, a, b, ord )
      call mh_distribution_derivs ( r, alpha, mu, sigma, n_1, n_2, &
        & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
        & dMu_dT, dSigma_dT, dN_2_dT, &
        & nr, dNR_dIWC, dNR_dT )
      answer = r * r * dNR_dIWC * xis(i_beta)
      maxOrd(i_T,i_IWC,i_f,i_beta+2) = max(maxOrd(i_T,i_IWC,i_f,i_beta+2),ord)
    end do
    ! Finished the quadrature.  Integral is in Answer.
    answer = pi2 * answer
    error = pi2 * work(1)
    nFunc(1,i_T,i_IWC,i_f,i_beta+2) = nfuse
    nFunc(2,i_T,i_IWC,i_f,i_beta+2) = -iopt(1)
    if ( iopt(1) > 0 ) call errorReport ( "dBeta_dIWC", iopt(1), answer, error )
  end subroutine Do_dint_dBeta_dIWC

  subroutine Do_dint_dBeta_dT ( R_min, R_max, Answer, Error )
  !{ Do an integration to get the derivative of beta w.r.t.\ $T$:
  !  $\int_0^\infty r^2 \left( \xi(r) \frac{\partial n(r)}{\partial T} +
  !    \frac{\partial \xi(r)}{\partial T} n(r) \right) \text{d} r$,
  !  where $n(r) = N_1(r) + N_2(r)$.
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8) :: NR, dNR_dIWC, dNR_dT
    ! Start the quadrature
    call dint1 ( r_min, r_max, answer, work, iopt )
    ! Evaluate the integrand
    do
      call dinta ( answer, work, iopt )
      if ( iopt(1) /= 0 ) exit
      call mie_efficiencies_derivs ( rf, chi_fac * r, dRF_dT, xi_e, xi_s, &
        & dXi_e_dT, dXi_s_dT, a, b, dA_dT, dB_dT, ord )
      call mh_distribution_derivs ( r, alpha, mu, sigma, n_1, n_2, &
        & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
        & dMu_dT, dSigma_dT, dN_2_dT, &
        & nr, dNR_dIWC, dNR_dT )
      answer = r * r * ( nr * dXis_dT(i_beta) + dNR_dT * xis(i_beta) )
      maxOrd(i_T,i_IWC,i_f,i_beta+4) = max(maxOrd(i_T,i_IWC,i_f,i_beta+4),ord)
    end do
    ! Finished the quadrature.  Integral is in Answer.
    answer = pi2 * answer
    error = pi2 * work(1)
    nFunc(1,i_T,i_IWC,i_f,i_beta+4) = nfuse
    nFunc(2,i_T,i_IWC,i_f,i_beta+4) = -iopt(1)
    if ( iopt(1) > 0 ) call errorReport ( "dBeta_dT", iopt(1), answer, error )
  end subroutine Do_dint_dBeta_dT

  subroutine Do_dint_ice ( R_min, R_max, Answer, Error )
  !{ Do an integration to get IWC\_total using both terms of the number
  !  distribution:\\
  ! $\frac43 \pi \rho_{\text{ice}} \int_0^\infty r^3 \left\{
  !    N_1 2 r \exp(-\alpha 2 r) +
  !    \frac{N_2}{2 r} \exp \left[ -\frac12
  !      \left( \frac{\log(2 r/D_0) - \mu}{\sigma} \right)^2 \right]
  !    \right\} \text{d} r$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8), parameter :: Pi83_r = 8.0/3.0 * Pi * rho_ice
    ! Start the quadrature
    call dint1 ( r_min, r_max, answer, work, iopt )
    ! Evaluate the integrand
    do
      call dinta ( answer, work, iopt )
      if ( iopt(1) /= 0 ) exit
      answer = r * r * r * mh_distribution(r,alpha,mu,sigma,n_1,n_2)
    end do
    ! Finished the quadrature.  Integral is in Answer.
    answer = Pi83_r * answer
    error = Pi83_r * work(1)
    nFuncI(1,i_T,i_IWC) = nfuse
    nFuncI(2,i_T,i_IWC) = -iopt(1)
    if ( iopt(1) > 0 ) call errorReport ( "IWC", iopt(1), answer, error )
  end subroutine Do_dint_ice

  subroutine Do_Cadre_phase ( R_min, R_max, Answer, Error )
  !{ Do the integration to get the integrated phase function
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}} \int_0^\infty n(r) p_0(\theta,r) \text{d} r$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8) :: tan_r ! Tan(r)
    iopt(1) = 0
    nfuse = 0
    if ( .not. tan_u ) then
      do
        call cadre_reverse ( answer, r_min, r_max, 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
      ! Evaluate the integrand
        call mie_efficiencies ( rf, chi_fac * answer, xi_e, xi_s, a, b, ord )
        answer = mh_distribution(answer,alpha,mu,sigma,n_1,n_2) * &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta))
        maxOrdP(i_T,i_IWC,i_theta,i_f,1) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,1),ord)
      end do
    else
      do
        call cadre_reverse ( answer, atan(r_min), atan(r_max), 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        tan_r = tan(answer)
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
      ! Evaluate the integrand
        call mie_efficiencies ( rf, chi_fac * tan_r, xi_e, xi_s, a, b, ord )
        answer = mh_distribution(tan_r,alpha,mu,sigma,n_1,n_2) * &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta)) * &
          & (1.0 + tan_r**2)
        maxOrdP(i_T,i_IWC,i_theta,i_f,1) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,1),ord)
      end do
    end if
    iopt(1) = iopt(1) + 1000
    ! Finished the quadrature.  Integral is in Answer, error is in work(1).
    nFuncP(1,i_T,i_IWC,i_theta,i_f,1) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,1) = iopt(1)
    answer = lambda**2 / (pi2*beta(i_T,i_IWC,i_f,i_c_s)) * answer
    !{ Error in $P$ depends on work(1) and error in $\beta_{c\_s}$.
    ! Suppose $x$ and $y$ are functions with errors $e$ and $f$ respectively.
    ! Let $a=e/x$ and $b=f/y$.  Then the error in $x/y$ is
    ! $\frac{x(1-a)}{y(1-b)} -\frac{x}y = \frac{x}y \frac{a+b}{1-b}$.  Substituting
    ! $e$ and $f$ we have $\frac1{|y|} \frac{|x|f + |y|e}{|y|-f}$.  Neglecting
    ! $f$ w.r.t.\ $y$, we have $\frac{x}y \left( f + \frac{x}y e \right)$.
    error = answer * ( work(1) + answer * eest(i_T,i_IWC,i_f,i_c_s) )
    if ( iopt(1) < 1004 ) then
      nFuncP(2,i_T,i_IWC,i_theta,i_f,1) = -iopt(1)
      return
    end if
    call errorReport ( "P", iopt(1), answer, error, i_theta )
  end subroutine Do_Cadre_phase

  subroutine Do_dint_phase ( R_min, R_max, Answer, Error, Only )
  !{ Do the integration to get the integrated phase function
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}} \int_0^\infty n(r) p_0(\theta,r) \text{d} r$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    character(*), intent(in) :: Only ! Only do the integral if present
    real(r8) :: tan_r ! Tan(r)
    if ( .not. tan_u ) then
    ! Start the quadrature
      call dint1 ( r_min, r_max, answer, work, iopt )
      ! Evaluate the integrand
      do
        call dinta ( answer, work, iopt )
        if ( iopt(1) /= 0 ) exit
        call mie_efficiencies ( rf, chi_fac * r, xi_e, xi_s, a, b, ord )
        answer = mh_distribution(r,alpha,mu,sigma,n_1,n_2) * &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta))
        if ( capture > 0 ) &
          & write ( capture, 666 ) work(1), answer, mh_distribution(r,alpha,mu,sigma,n_1,n_2), &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta))
666     format (f8.2,1p,5g15.6)
        maxOrdP(i_T,i_IWC,i_theta,i_f,1) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,1),ord)
      end do
    else
      call dint1 ( atan(r_min), atan(r_max), answer, work, iopt )
      do
        call dinta ( answer, work, iopt )
        tan_r = tan(r)
        if ( iopt(1) /= 0 ) exit
        call mie_efficiencies ( rf, chi_fac * tan_r, xi_e, xi_s, a, b, ord )
        answer = mh_distribution(tan_r,alpha,mu,sigma,n_1,n_2) * &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta)) * &
          & (1.0 + tan_r**2)
        if ( capture > 0 ) &
          & write ( capture, 666 ) tan_r, answer, mh_distribution(tan_r,alpha,mu,sigma,n_1,n_2), &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta)), &
          & 1.0 + tan_r**2, r
        maxOrdP(i_T,i_IWC,i_theta,i_f,1) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,1),ord)
      end do
    end if
    ! Finished the quadrature.  Integral is in Answer, error is in work(1).
    if ( capture > 0 ) write ( capture, 666 ) -1.0, answer, work(1)
    nFuncP(1,i_T,i_IWC,i_theta,i_f,1) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,1) = -iopt(1)
    if ( iopt(1) <= 0 ) then
      if ( only(1:1) /= '' ) then
        error = work(1)
        return
      end if
      answer = lambda**2 / (pi2*beta(i_T,i_IWC,i_f,i_c_s)) * answer
      !{ Error in $P$ depends on work(1) and error in $\beta_{c\_s}$.
      ! Suppose $x$ and $y$ are functions with errors $e$ and $f$ respectively.
      ! Let $a=e/x$ and $b=f/y$.  Then the error in $x/y$ is
      ! $\frac{x(1-a)}{y(1-b)} -\frac{x}y = \frac{x}y \frac{a+b}{1-b}$.  Substituting
      ! $e$ and $f$ we have $\frac1{|y|} \frac{|x|f + |y|e}{|y|-f}$.  Neglecting
      ! $f$ w.r.t.\ $y$, we have $\frac{x}y \left( f + \frac{x}y e \right)$.
      error = answer * ( work(1) + answer * eest(i_T,i_IWC,i_f,i_c_s) )
      return
    end if
    call errorReport ( "P", iopt(1), answer, work(1), i_theta )
  end subroutine Do_dint_phase

  subroutine Do_Cadre_dPhase_dIWC ( R_min, R_max, Answer, Error )
  !{ Do the integration to get the derivative of the integrated phase function
  ! w.r.t.~IWC.\\
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}}
  !  \int_0^\infty p_0(\theta,r) \frac{\partial n(r)}{\partial IWC}
  !          \text{d} r
  !  -\frac{P(\theta)}{\beta_{c\_s}} \frac{\partial \beta_{c\_s}}{\partial IWC}$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8) :: Numer        ! of outer integral, then error in numer of answer
    real(r8) :: NR, dNR_dIWC, dNR_dT
    real(r8) :: tan_r ! Tan(r)
    iopt(1) = 0
    nfuse = 0
    if ( .not. tan_u ) then
      do
        call cadre_reverse ( answer, r_min, r_max, 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
        ! Evaluate the integrand
        call mie_efficiencies ( rf, chi_fac * answer, xi_e, xi_s, a, b, ord )
        call mh_distribution_derivs ( answer, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        answer = dNR_dIWC * phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta) )
        maxOrdP(i_T,i_IWC,i_theta,i_f,2) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,2),ord)
      end do
    else
      do
        call cadre_reverse ( answer, atan(r_min), atan(r_max), 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        tan_r = tan(answer)
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
        ! Evaluate the integrand
        call mie_efficiencies ( rf, chi_fac * tan_r, xi_e, xi_s, a, b, ord )
        call mh_distribution_derivs ( tan_r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        answer = dNR_dIWC * phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta)) * &
               & ( 1.0 + tan_r**2 )
        maxOrdP(i_T,i_IWC,i_theta,i_f,2) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,2),ord)
      end do
    end if
    iopt(1) = iopt(1) + 1000
    ! Finished the quadrature.  Integral is in Answer, error is in work(1).
    nFuncP(1,i_T,i_IWC,i_theta,i_f,2) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,2) = iopt(1)
    answer = ( lambda**2 / pi2 * answer - &
           &   p(i_T,i_IWC,i_Theta,i_f) * dBeta_dIWC(i_T,i_IWC,i_f,i_c_s) ) / &
           & beta(i_T,i_IWC,i_f,i_c_s)
    !{ Final error depends on work(1), error in $\beta_{c\_s}$, error in
    ! $P(\theta)$ and error in $\frac{\partial \beta_{c\_s}}{\partial
    ! \text{IWC}}$. Let $x$ and $y$ be quantities with errors $e$ and $f$. 
    ! Then neglecting $ef$ the error in $xy$ is $xf + ye$. See {\tt
    ! do\_dint\_phase} for the error in a quotient.
    numer = work(1) + & ! Error in numerator of answer
          &   e_p(i_T,i_IWC,i_Theta,i_f) * abs(dBeta_dIWC(i_T,i_IWC,i_f,i_c_s)) + &
          &   abs(p(i_T,i_IWC,i_Theta,i_f)) * e_dBeta_dIWC(i_T,i_IWC,i_f,i_c_s)
    error = abs(answer) * ( numer + abs(answer) * eest(i_T,i_IWC,i_f,i_c_s) )
    if ( iopt(1) < 1004 ) then
      nFuncP(2,i_T,i_IWC,i_theta,i_f,1) = -iopt(1)
      return
    end if
    call errorReport ( "dP_dIWC", iopt(1), answer, error, i_theta )
  end subroutine Do_Cadre_dPhase_dIWC

  subroutine Do_dint_dPhase_dIWC ( R_min, R_max, Answer, Error, Only )
  !{ Do the integration to get the derivative of the integrated phase function
  ! w.r.t.~IWC.\\
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}}
  !  \int_0^\infty p_0(\theta,r) \frac{\partial n(r)}{\partial IWC}
  !          \text{d} r
  !  -\frac{P(\theta)}{\beta_{c\_s}} \frac{\partial \beta_{c\_s}}{\partial IWC}$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    character(*), intent(in) :: Only ! Only do the integral if present
    real(r8) :: Numer        ! of outer integral, then error in numer of answer
    real(r8) :: NR, dNR_dIWC, dNR_dT
    real(r8) :: Tan_r ! tan(r)
    if ( .not. tan_u ) then
      ! Start the quadrature
      call dint1 ( r_min, r_max, answer, work, iopt )
      ! Evaluate the integrand
      do
        call dinta ( answer, work, iopt )
        if ( iopt(1) /= 0 ) exit
        call mie_efficiencies ( rf, chi_fac * r, xi_e, xi_s, a, b, ord )
        call mh_distribution_derivs ( r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        answer = dNR_dIWC * phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta) )
        if ( capture > 0 ) &
          & write ( capture, '(f8.2,1p,3g15.6)' ) work(1), answer, dNR_dIWC, &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta))
        maxOrdP(i_T,i_IWC,i_theta,i_f,2) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,2),ord)
      end do
    else
      ! Start the quadrature
      call dint1 ( atan(r_min), atan(r_max), answer, work, iopt )
      ! Evaluate the integrand
      do
        call dinta ( answer, work, iopt )
        if ( iopt(1) /= 0 ) exit
        tan_r = tan(r)
        call mie_efficiencies ( rf, chi_fac * tan_r, xi_e, xi_s, a, b, ord )
        call mh_distribution_derivs ( tan_r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        answer = dNR_dIWC * phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta)) * &
               & ( 1.0 + tan_r**2 )
        if ( capture > 0 ) &
          & write ( capture, '(f8.2,1p,3g15.6)' ) work(1), answer, dNR_dIWC, &
          & phase(theta,a(:ord),b,c1,c2,w,p1(:,i_theta),dP1_dTheta(:,i_theta))
        maxOrdP(i_T,i_IWC,i_theta,i_f,2) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,2),ord)
      end do
    end if
    ! Finished the quadrature.  Integral is in Answer.
    nFuncP(1,i_T,i_IWC,i_theta,i_f,2) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,2) = -iopt(1)
    if ( iopt(1) <= 0 ) then
      if ( only(1:1) /= '' ) then
        error = work(1)
        return
      end if
      answer = ( lambda**2 / pi2 * answer - &
             &   p(i_T,i_IWC,i_Theta,i_f) * dBeta_dIWC(i_T,i_IWC,i_f,i_c_s) ) / &
             & beta(i_T,i_IWC,i_f,i_c_s)
      !{ Final error depends on work(1), error in $\beta_{c\_s}$, error in
      ! $P(\theta)$ and error in $\frac{\partial \beta_{c\_s}}{\partial
      ! \text{IWC}}$. Let $x$ and $y$ be quantities with errors $e$ and $f$. 
      ! Then neglecting $ef$ the error in $xy$ is $xf + ye$. See {\tt
      ! do\_dint\_phase} for the error in a quotient.
      numer = work(1) + & ! Error in numerator of answer
            &   e_p(i_T,i_IWC,i_Theta,i_f) * abs(dBeta_dIWC(i_T,i_IWC,i_f,i_c_s)) + &
            &   abs(p(i_T,i_IWC,i_Theta,i_f)) * e_dBeta_dIWC(i_T,i_IWC,i_f,i_c_s)
      error = abs(answer) * ( numer + abs(answer) * eest(i_T,i_IWC,i_f,i_c_s) )
      return
    end if
    call errorReport ( "dP_dIWC", iopt(1), answer, error, i_theta )
  end subroutine Do_dint_dPhase_dIWC

  subroutine Do_Cadre_dPhase_dT ( R_min, R_max, Answer, Error )
  !{ Do the integration to get the derivative of the integrated phase function
  ! w.r.t.~T.\\
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}}
  !  \int_0^\infty n(r) \frac{\partial p_0(\theta,r)}{\partial T} +
  !                \frac{\partial n(r)}{\partial T} p_0(\theta,r)
  !                \, \text{d} r
  !  -\frac{P(\theta)}{\beta_{c\_s}} \frac{\partial \beta_{c\_s}}{\partial T}$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    real(r8) :: Numer          ! Error in numerator of answer
    real(r8) :: P0, dP0_dT
    real(r8) :: NR, dNR_dIWC, dNR_dT
    real(r8) :: Tan_r ! tan(r)
    iopt(1) = 0
    nfuse = 0
    if ( .not. tan_u ) then
      do
        call cadre_reverse ( answer, r_min, r_max, 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
        ! Evaluate the integrand
        call mie_efficiencies_derivs ( rf, chi_fac * answer, dRF_dT, &
          & xi_e, xi_s, dXi_e_dT, dXi_s_dT, a, b, dA_dT, dB_dT, ord )
        maxOrdP(i_T,i_IWC,i_theta,i_f,3) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,3),ord)
        call mh_distribution_derivs ( answer, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        call phase_deriv ( theta, a(:ord), b, dA_dT, dB_dT, c1, c2, w, &
          & p0, dp0_dT, p1(:,i_theta), dP1_dTheta(:,i_theta) )
        answer = nr * dp0_dT + dNR_dT * p0
      end do
    else
      do
        call cadre_reverse ( answer, atan(r_min), atan(r_max), 1.0e-13_r8, 1.0e-13_r8, &
          & level, work(1), iopt(1) )
        if ( iopt(1) > 0 ) exit
        tan_r = tan(tan_r)
        nfuse = nfuse + 1
        if ( nfuse > nfmax ) then
          iopt(1) = 6
          exit
        end if
        ! Evaluate the integrand
        call mie_efficiencies_derivs ( rf, chi_fac * tan_r, dRF_dT, &
          & xi_e, xi_s, dXi_e_dT, dXi_s_dT, a, b, dA_dT, dB_dT, ord )
        maxOrdP(i_T,i_IWC,i_theta,i_f,3) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,3),ord)
        call mh_distribution_derivs ( tan_r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        call phase_deriv ( theta, a(:ord), b, dA_dT, dB_dT, c1, c2, w, &
          & p0, dp0_dT, p1(:,i_theta), dP1_dTheta(:,i_theta) )
        answer = nr * dp0_dT + dNR_dT * p0 * ( 1.0 + tan_r**2 )
      end do
    end if
    iopt(1) = iopt(1) + 1000
    ! Finished the quadrature.  Integral is in Answer, error is in work(1).
    nFuncP(1,i_T,i_IWC,i_theta,i_f,3) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,3) = iopt(1)
    answer = ( lambda**2 /pi2 * answer - &
           &   p(i_T,i_IWC,i_Theta,i_f) * dBeta_dT(i_T,i_IWC,i_f,i_c_s) ) / &
           & beta(i_T,i_IWC,i_f,i_c_s)
    !{ Final error depends on work(1), error in $\beta_{c\_s}$, error in
    ! $P(\theta)$ and error in $\frac{\partial \beta_{c\_s}}{\partial T}$.
    ! Let $x$ and $y$ be quantities with errors $e$ and $f$.  Then neglecting
    ! $ef$ the error in $xy$ is $xf + ye$. See {\tt do\_dint\_phase} for the
    ! error in a quotient.
    numer = work(1) + & ! Error in numerator of answer
          &   e_p(i_T,i_IWC,i_Theta,i_f) * abs(dBeta_dT(i_T,i_IWC,i_f,i_c_s)) + &
          &   abs(p(i_T,i_IWC,i_Theta,i_f)) * e_dBeta_dT(i_T,i_IWC,i_f,i_c_s)
    error = abs(answer) * ( numer + abs(answer) * eest(i_T,i_IWC,i_f,i_c_s) )
    if ( iopt(1) < 1004 ) then
      nFuncP(2,i_T,i_IWC,i_theta,i_f,1) = -iopt(1)
      return
    end if
    call errorReport ( "dP_dT", iopt(1), answer, error, i_theta )
  end subroutine Do_Cadre_dPhase_dT

  subroutine Do_dint_dPhase_dT ( R_min, R_max, Answer, Error, Only )
  !{ Do the integration to get the derivative of the integrated phase function
  ! w.r.t.~T.\\
  ! $\frac{\lambda^2}{2\pi\beta_{c\_s}}
  !  \int_0^\infty n(r) \frac{\partial p_0(\theta,r)}{\partial T} +
  !                \frac{\partial n(r)}{\partial T} p_0(\theta,r)
  !                \, \text{d} r
  !  -\frac{P(\theta)}{\beta_{c\_s}} \frac{\partial \beta_{c\_s}}{\partial T}$
    real(r8), intent(in) :: R_min, R_max
    real(r8), intent(inout) :: Answer
    real(r8), intent(out) :: Error
    character(*), intent(in) :: Only ! Only do the integral if present
    real(r8) :: Numer          ! Error in numerator of answer
    real(r8) :: P0, dP0_dT
    real(r8) :: NR, dNR_dIWC, dNR_dT
    real(r8) :: Tan_r ! tan(r)
    if ( .not. tan_u ) then
      ! Start the quadrature
      call dint1 ( r_min, r_max, answer, work, iopt )
      ! Evaluate the integrand
      do
        call dinta ( answer, work, iopt )
        if ( iopt(1) /= 0 ) exit
        call mie_efficiencies_derivs ( rf, chi_fac * r, dRF_dT, &
          & xi_e, xi_s, dXi_e_dT, dXi_s_dT, a, b, dA_dT, dB_dT, ord )
        maxOrdP(i_T,i_IWC,i_theta,i_f,3) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,3),ord)
        call mh_distribution_derivs ( r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        call phase_deriv ( theta, a(:ord), b, dA_dT, dB_dT, c1, c2, w, &
          & p0, dp0_dT, p1(:,i_theta), dP1_dTheta(:,i_theta) )
        answer = nr * dp0_dT + dNR_dT * p0
        if ( capture > 0 ) &
          & write ( capture, '(f8.2,1p,5g15.6)' ) work(1), answer, nr, dp0_dT, &
          & dNR_dT, p0
      end do
    else
      ! Start the quadrature
      call dint1 ( atan(r_min), atan(r_max), answer, work, iopt )
      ! Evaluate the integrand
      do
        call dinta ( answer, work, iopt )
        if ( iopt(1) /= 0 ) exit
        tan_r = tan(r)
        call mie_efficiencies_derivs ( rf, chi_fac * tan_r, dRF_dT, &
          & xi_e, xi_s, dXi_e_dT, dXi_s_dT, a, b, dA_dT, dB_dT, ord )
        maxOrdP(i_T,i_IWC,i_theta,i_f,3) = max(maxOrdP(i_T,i_IWC,i_theta,i_f,3),ord)
        call mh_distribution_derivs ( tan_r, alpha, mu, sigma, n_1, n_2, &
          & dAlpha_dIWC, dMu_dIWC, dSigma_dIWC, dN_1_dIWC, dN_2_dIWC, &
          & dMu_dT, dSigma_dT, dN_2_dT, &
          & nr, dNR_dIWC, dNR_dT )
        call phase_deriv ( theta, a(:ord), b, dA_dT, dB_dT, c1, c2, w, &
          & p0, dp0_dT, p1(:,i_theta), dP1_dTheta(:,i_theta) )
        answer = nr * dp0_dT + dNR_dT * p0 * ( 1.0 + tan_r**2 )
        if ( capture > 0 ) &
          & write ( capture, '(f8.2,1p,5g15.6)' ) work(1), answer, nr, dp0_dT, &
          & dNR_dT, p0
      end do
    end if
    ! Finished the quadrature.  Integral is in Answer.
    nFuncP(1,i_T,i_IWC,i_theta,i_f,3) = nfuse
    nFuncP(2,i_T,i_IWC,i_theta,i_f,3) = -iopt(1)
    if ( iopt(1) <= 0 ) then
      if ( only(1:1) /= '' ) then
        error = work(1)
        return
      end if
      answer = ( lambda**2 /pi2 * answer - &
             &   p(i_T,i_IWC,i_Theta,i_f) * dBeta_dT(i_T,i_IWC,i_f,i_c_s) ) / &
             & beta(i_T,i_IWC,i_f,i_c_s)
      !{ Final error depends on work(1), error in $\beta_{c\_s}$, error in
      ! $P(\theta)$ and error in $\frac{\partial \beta_{c\_s}}{\partial T}$.
      ! Let $x$ and $y$ be quantities with errors $e$ and $f$.  Then neglecting
      ! $ef$ the error in $xy$ is $xf + ye$. See {\tt do\_dint\_phase} for the
      ! error in a quotient.
      numer = work(1) + & ! Error in numerator of answer
            &   e_p(i_T,i_IWC,i_Theta,i_f) * abs(dBeta_dT(i_T,i_IWC,i_f,i_c_s)) + &
            &   abs(p(i_T,i_IWC,i_Theta,i_f)) * e_dBeta_dT(i_T,i_IWC,i_f,i_c_s)
      error = abs(answer) * ( numer + abs(answer) * eest(i_T,i_IWC,i_f,i_c_s) )
      return
    end if
    call errorReport ( "dP_dT", iopt(1), answer, error, i_theta )
  end subroutine Do_dint_dPhase_dT

  subroutine ErrorReport ( What, Iflag, Answer, Error, I_Theta )
    character(*), intent(in) :: What
    integer, intent(in) :: Iflag
    real(r8), intent(out) :: Answer, Error
    integer, intent(in), optional :: I_Theta
    write ( *, '(a)' ) repeat("*",72)
    write ( *, '(a,"(")', advance='no' ) trim(what)
    write ( *, '(i0,",",i0)', advance='no' ) i_t, i_iwc
    if ( present(i_theta) ) write ( *, '(",",i0)', advance='no' ) i_theta
    write ( *, '(",",i0,")")', advance='no' ) i_f
    select case ( iflag )
    case ( -5 )
      write ( *, '(a)' ) ": Incorrect usage: WORK too small"
    case ( -6, 4 )
      write ( *, '(a)' ) ": Incorrect usage: Bad value for an element of IOPT"
    case ( -7, 5 )
      write ( *, '(a)' ) ": Too many integrand values required"
    case ( -8, -9 )
      write ( *, '(a,i0)' ) ": Non integrable singularity in dimension ", 10+Iflag
    case ( 6 )
      write ( *, '(a)' ) ": Non integrable singularity"
    case ( 1004 )
      write ( *, '(a)' ) ": CADRE ran out of stack space"
    case ( 1005 )
      write ( *, '(a)' ) ": CADRE required too small a subinterval"
    case default
      write ( *, '(a,i0,a)' ) ": Why was IFLAG = ", iflag, " produced?"
    end select
    write ( *, '(a)' ) repeat("*",72)
    answer = huge(answer)
    error = huge(error)
  end subroutine ErrorReport

  subroutine PrintResults
    write ( *, 1 ) f_s(:n_f)
1   format ( 'Frequencies = ', 8f7.1, (10f7.1) )

    write ( *, 2 ) t_min, t_max, n_t, iwc_min, iwc_max, n_iwc, &
      &            r_min, r_max, theta_min, theta_max, n_theta, n_cut
2   format ( 'T = ', f5.0, ' : ', f5.0, ' (', i0, &
           & ') Log10 IWC = ', f6.2, ' : ', f6.2, ' (', i0, ')' / &
           & 'R = ', f6.1, ' : ', f6.1, ' Theta = ', f7.3, ' : ', f7.3, &
           & ' (', i0, ')', ' N_Cut = ', i0 )

    if ( wantIWC ) &
      & call report_2 ( 'IWC_total', iwc_tot, .false., eestI, nFuncI )
    if ( wantBeta ) then
      call report ( 'Beta(c_e)', beta(:,:,:,i_c_e), norm, eest(:,:,:,i_c_e), &
        & nFunc(:,:,:,:,i_c_e), maxOrd(:,:,:,i_c_e) )
      call report ( 'Beta(c_s)', beta(:,:,:,i_c_s), norm, eest(:,:,:,i_c_s), &
        & nFunc(:,:,:,:,i_c_s), maxOrd(:,:,:,i_c_s) )
      if ( derivs ) then  
        call report ( 'dBeta(c_e)/dIWC', dBeta_dIWC(:,:,:,i_c_e), norm, &
          & e_dBeta_dIWC(:,:,:,i_c_e), nFunc(:,:,:,:,3), maxOrd(:,:,:,3) )
        if ( diffs .and. n_iwc > 1 ) then
          diffX_IWC = beta(:,2:n_iwc,:,i_c_e)-beta(:,:n_iwc-1,:,i_c_e)
          do i_IWC = 2, n_iwc
            diffX_IWC(:,i_iwc-1,:) = diffX_IWC(:,i_iwc-1,:) / diffIWC(i_iwc-1)
          end do
          call report ( 'Diff Beta(c_e) / diff IWC', diffX_IWC, .false. )
        end if
        call report ( 'dBeta(c_s)/dIWC', dBeta_dIWC(:,:,:,i_c_s), norm, &
          & e_dBeta_dIWC(:,:,:,i_c_s), nFunc(:,:,:,:,4), maxOrd(:,:,:,4) )
        if ( diffs .and. n_iwc > 1 ) then
          diffX_IWC = beta(:,2:n_iwc,:,i_c_s)-beta(:,:n_iwc-1,:,i_c_s)
          do i_IWC = 2, n_iwc
            diffX_IWC(:,i_iwc-1,:) = diffX_IWC(:,i_iwc-1,:) / diffIWC(i_iwc-1)
          end do
          call report ( 'Diff Beta(c_s) / diff IWC', diffX_IWC, .false. )
        end if
        call report ( 'dBeta(c_e)/dT', dBeta_dT(:,:,:,i_c_e), norm, &
          & e_dBeta_dT(:,:,:,i_c_e), nFunc(:,:,:,:,5), maxOrd(:,:,:,5) )
        if ( diffs .and. n_t > 1 ) call report ( 'Diff Beta(c_e) / Diff T', &
          & (beta(2:n_t,:,:,i_c_e)-beta(:n_t-1,:,:,i_c_e)) / dT, .false. )
        call report ( 'dBeta(c_s)/dT', dBeta_dT(:,:,:,i_c_s), norm, &
          & e_dBeta_dT(:,:,:,i_c_s), nFunc(:,:,:,:,6), maxOrd(:,:,:,6) )
        if ( diffs .and. n_t > 1 ) call report ( 'Diff Beta(c_s) / Diff T', &
              & (beta(2:n_t,:,:,i_c_s)-beta(:n_t-1,:,:,i_c_s)) / dT, .false. )
      end if ! Derivs
    end if ! wantBeta

    if ( wantP ) then
      do i_Theta = 1, n_Theta
        theta = theta_s(i_theta)
        write ( string, '(f7.3)' ) theta/deg2Rad
        string = adjustl(string)
        n = len_trim(string)                     ! trim trailing blanks
        n = verify(string(:n), '0', back=.true.) ! Trim trailing zeros
        n = verify(string(:n), '.', back=.true.) ! Trim trailing decimal point
        call report ( "P(" // string(:n) // ")", &
          & p(:,:,i_theta,:), norm, e_p(:,:,i_theta,:), nFuncP(:,:,:,i_theta,:,1), &
          & maxOrdP(:,:,i_theta,:,1) )
        if ( derivs ) then
          call report ( "dP(" // string(:n) // ")/dIWC", &
          & dP_dIWC(:,:,i_theta,:), norm, e_dP_dIWC(:,:,i_theta,:), &
          & nFuncP(:,:,:,i_theta,:,2), maxOrdP(:,:,i_theta,:,2) )
          if ( diffs .and. n_iwc > 1 ) then
            diffX_IWC = p(:,2:n_iwc,i_theta,:) - p(:,:n_iwc-1,i_theta,:)
            do i_IWC = 2, n_iwc
              diffX_IWC(:,i_iwc-1,:) = diffX_IWC(:,i_iwc-1,:) / diffIWC(i_iwc-1)
            end do
            call report ( "Diff P(" // string(:n) // ")/diff IWC", &
              & diffX_IWC, .false. )
          end if
          call report ( "dP(" // string(:n) // ")/dT", &
            & dP_dT(:,:,i_theta,:), norm, e_dP_dT(:,:,i_theta,:), &
            & nFuncP(:,:,:,i_theta,:,3), maxOrdP(:,:,i_theta,:,3) )
          if ( diffs .and. n_t > 1 ) call report ( "Diff P(" // string(:n) // ")/diff T", &
                & (p(2:n_t,:,i_theta,:) - p(:n_t-1,:,i_theta,:)) / dT, .false. )
        end if
      end do ! Theta
    end if ! wantP
  end subroutine PrintResults

  subroutine Progress_Report ( What, Value, Error, NFunc, Sub3 )
    character(*), intent(in) :: What
    real(r8), intent(in) :: Value
    real(r8), intent(in) :: Error
    integer, intent(in) :: NFunc(:)
    integer, intent(in), optional :: Sub3
    integer, save :: Lines = 0
    if ( .not. progress ) return
    if ( lines == 0 ) write ( *, 1 )
1   format ( 'What', t22, '   Integral     Error        NFunc  *  Time    T(C) log IWC  Theta    F' )
    call cpu_time ( t1 )
    if ( present(sub3) ) then
      write ( *, 2 ) what, i_t, i_iwc, sub3, i_f, value, error, nfunc(1), -nfunc(2), &
        &      t1-t0, t_s(i_t)-c_0, iwc_s(i_iwc), theta_s(sub3)/deg2rad, f_s(i_f)
2     format ( a, '(', 3(i0,','), i0, ')', t22, 1pg13.6, g13.6, i8, i3, 0p, f6.2, &
        &      f9.3, f7.2, 0p, f7.2, f8.2 )
    else
      write ( *, 3 ) what, i_t, i_iwc,       i_f, value, error, nfunc(1), -nfunc(2), &
        &      t1-t0, t_s(i_t)-c_0, iwc_s(i_iwc),                        f_s(i_f)
3     format ( a, '(', 2(i0,','), i0, ')', t22, 1pg13.6, g13.6, i8, i3, 0p, f6.2, &
        &      f9.3, f7.2, 0p, 7x,   f8.2 )
    end if
    lines = mod(lines + 1, 50)
    t0 = t1
  end subroutine Progress_Report

  subroutine Report_2 ( Title, Value, Norm, EEst, NFunc, MaxOrd )
    ! For two-D values
    character(len=*), intent(in) :: Title
    real(r8), intent(in) :: Value(:,:)
    logical, intent(in) :: Norm   ! Report also Value/IWC_Tot
    real(r8), intent(in), optional :: EEst(:,:)
    integer, intent(in), optional :: NFunc(:,:,:)
    integer, intent(in), optional :: MaxOrd(:,:) ! Ignored if NFunc is absent
    integer :: I, J
1   format ( a, " \", i0 )
2   format ( 1p, 10g13.5 )
3   format ( 10(i7,i2,i4) )
4   format ( 10(i8,i5) )

    write ( *, 1 ) title, size(value)
    do i = 1, size(value,2)
      write ( *, 2 ) value(:,i)
    end do
    if ( norm .and. wantIWC ) then
      write ( *, 1 ) title // ' / IWC total', size(value)
      do i = 1, size(iwc_tot,2)
        write ( *, 2 ) value(:,i) / iwc_tot(:,i)
      end do
    end if
    if ( details > 0 ) then
      if ( present(eest) ) then
        write ( *, 1 ) 'Error estimate', size(value)
        do i = 1, size(eest,2); write ( *, 2 ) eest(:,i); end do
      end if
      if ( present(nFunc) ) then
        if ( present(maxOrd) ) then
          write ( *, 1 ) &
            & 'Number of integrands, -flag, maximum order of Bessel functions', &
            &  3*size(value)
          do i = 1, size(nFunc,3)
            write ( *, 3 ) ( nFunc(:,j,i), maxOrd(j,i), j = 1, size(nFunc,2) )
          end do
        else
          write ( *, 1 ) 'Number of integrands, -flag', 2*size(value)
          do i = 1, size(nFunc,3)
            write ( *, 4 ) nFunc(:,:,i)
          end do
        end if
      end if
    end if
  end subroutine Report_2

  subroutine Report ( Title, Value, Norm, EEst, NFunc, MaxOrd )
    ! For three-D values
    character(len=*), intent(in) :: Title
    real(r8), intent(in) :: Value(:,:,:)
    logical, intent(in) :: Norm   ! Report also Value/IWC_Tot
    real(r8), intent(in), optional :: EEst(:,:,:)
    ! Either both of NFunc and MaxOrd, or neither, is assumed
    integer, intent(in), optional :: NFunc(:,:,:,:)
    integer, intent(in), optional :: MaxOrd(:,:,:)
    integer :: I, I_F, J
1   format ( a, ", F = ", f6.1, " \", i0 )
2   format ( a, a, " \", i0 )
3   format ( 1p, 10g13.5 )
4   format ( a, " \", i0 )
5   format ( 10(i7,i2,i4) )

    if ( size(value,1) == 1 ) then
      if ( size(value,2) == 1 ) then
        write ( *, 2 ) title, ' F ', size(value)
        write ( *, 3 ) value
        if ( norm .and. wantIWC ) then
          write ( *, 2 ) title, ' / IWC total', size(value)
          write ( *, 3 ) value(1,1,:) / iwc_tot(1,1)
        end if
        if ( details > 0 ) then
          if ( present(eest) ) then
            write ( *, 4 ) 'Error estimate', size(value)
            write ( *, 3 ) eest
          end if
          if ( present(nFunc) ) then
            write ( *, 4 ) &
              & 'Number of integrands, -flag, maximum order of Bessel functions', &
              & 3*size(value)
            write ( *, 5 ) (nFunc(:,1,1,j), maxOrd(1,1,j), j=1, size(nFunc,4))
          end if
        end if
      else
        write ( *, 2 ) title, ' IWC X F ', size(value)
        do i = 1, size(value,3)
          write ( *, 3 ) value(1,:,i)
        end do
        if ( norm .and. wantIWC ) then
          write ( *, 2 ) title, ' / IWC total', size(value)
          do i = 1, size(iwc_tot,2)
            write ( *, 3 ) value(1,:,i) / iwc_tot(:,i)
          end do
        end if
        if ( details > 0 ) then
          if ( present(eest) ) then
            write ( *, 4 ) 'Error estimate', size(value)
            do i = 1, size(eest,3); write ( *, 3 ) eest(1,:,i); end do
          end if
          if ( present(nFunc) ) then
            write ( *, 4 ) &
              & 'Number of integrands, -flag, maximum order of Bessel functions', &
              & 3*size(value)
            do i = 1, size(nFunc,4)
              write ( *, 5 ) (nFunc(:,1,j,i), maxord(1,j,i), j=1,size(nFunc,3))
            end do
          end if
        end if
      end if
    else if ( size(value,2) == 1 ) then
      write ( *, 2 ) title, ' T X F ', size(value)
      do i = 1, size(value,3)
        write ( *, 3 ) value(:,1,i)
      end do
      if ( norm .and. wantIWC ) then
        write ( *, 2 ) title, ' / IWC total', size(value)
        do i = 1, size(iwc_tot,2)
          write ( *, 3 ) value(:,1,i) / iwc_tot(:,i)
        end do
      end if
      if ( details > 0 ) then
        if ( present(eest) ) then
          write ( *, 4 ) 'Error estimate', size(value)
          do i = 1, size(eest,3); write ( *, 3 ) eest(:,1,i); end do
        end if
        if ( present(nFunc) ) then
          write ( *, 4 ) &
            & 'Number of integrands, -flag, maximum order of Bessel functions', &
            & 3*size(value)
          do i = 1, size(nFunc,4)
            write ( *, 5 ) (nFunc(:,j,1,i), maxord(j,1,i), j=1, size(nFunc,2))
          end do
        end if
      end if
    else
      do i_f = 1, size(value,3)
        write ( *, 1 ) title, f_s(i_f), size(value(:,:,1))
        do i = 1, size(value,2)
          write ( *, 3 ) value(:,i,i_f)
        end do
        if ( norm .and. wantIWC ) then
          write ( *, 2 ) title, ' / IWC total', size(value(:,:,1))
          do i = 1, size(iwc_tot,2)
            write ( *, 3 ) value(:,i,i_f) / iwc_tot(:,i)
          end do
        end if
        if ( details > 0 ) then
          if ( present(eest) ) then
            write ( *, 4 ) 'Error estimate', size(value(:,:,1))
            do i = 1, size(eest,2); write ( *, 3 ) eest(:,i,i_f); end do
          end if
          if ( present(nFunc) ) then
            write ( *, 4 ) &
              & 'Number of integrands, -flag, maximum order of Bessel functions', &
              & 3*size(value(:,:,1))
            do i = 1, size(nFunc,3)
              write ( *, 5 ) (nFunc(:,j,i,i_f), maxOrd(j,i,i_f), j=1, size(nFunc,2))
            end do
          end if
        end if
      end do ! i_f
    end if
  end subroutine Report

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end program Mie_Tables

! $Log$
! Revision 1.4  2009/07/01 20:18:35  vsnyder
! Correct some LaTeX
!
! Revision 1.3  2008/06/05 02:20:09  vsnyder
! Added HDF output, added explicit frequencies
!
! Revision 1.2  2008/05/22 01:56:31  vsnyder
! List of frequencies, HDF output
!
! Revision 1.1  2008/04/19 01:15:27  vsnyder
! Initial commit
!
