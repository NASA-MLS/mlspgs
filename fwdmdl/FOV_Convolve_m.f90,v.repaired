! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module FOV_Convolve_m

  ! Antenna field of view convolutions

  use AntennaPatterns_m, only: AntennaPattern_T
  use MLSKinds, only: R8
  use MLSNumerics, only: Coefficients

  implicit NONE
  private
  public :: FOV_Convolve_Setup, FOV_Convolve_1D, FOV_Convolve_2d
  public :: FOV_Convolve_Temp_Derivs, FOV_Convolve_Teardown
  public :: FOV_Convolve_Temp_Derivs_Normalization
  public :: Dump, Dump_Convolve_Support
  public :: AntennaPattern_T, Coefficients ! for full f95 compatibility

  interface Dump
    module procedure Dump_Convolve_Support
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! FFT-related parameters
  integer, parameter, public :: pwr=12, no_fft=2**pwr, ffth = no_fft / 2

  ! Results of setting up for convolution
  type, public :: Convolve_Support_T
    real(r8), dimension(no_fft) :: angles ! Basis for FFT
    type(antennaPattern_t), pointer :: AntennaPattern ! For temperature derivs
    type(coefficients(r8)) :: Coeffs_1 ! for chi_in-init_angle -> angles(ffth:no_fft)
    type(coefficients(r8)) :: Coeffs_2 ! for angles(ffth:no_fft)-ang_step -> chi_out-init_angle
    real(r8), pointer :: del_chi_in(:) => NULL(), del_chi_out(:) => NULL()
    real(r8) :: Init_angle
    real(r8), dimension(no_fft) :: p, dp  ! From antenna pattern
  end type Convolve_Support_T

contains

  ! --------------------------------------  Dump_Convolve_Support  -----
  subroutine Dump_Convolve_Support ( cs, Name, Details )
    use Dump_0, only: DUMP
    use MLSNumerics, only: Dump
    use Output_m, only: OUTPUT

    type (convolve_support_t), intent(in) :: cs
    character(len=*), optional, intent(in) :: Name
    integer, optional, intent(in) :: Details ! 0 => no coeffs (default)
    integer :: My_Details

    if ( present(name) ) then
      call output ( name, advance='yes' )
    else
      call output ( 'Convolve support:', advance='yes' )
    end if
    call dump ( cs%angles, name="Angles" )
    my_Details = 0
    if ( present(details) ) my_Details = details
    if ( my_Details > 0 ) then
      call dump ( cs%coeffs_1, name='Coefficients 1' )
      call dump ( cs%coeffs_2, name='Coefficients 2' )
    end if
    call dump ( cs%del_chi_in, name="Del_Chi_In" )
    call dump ( cs%del_chi_out, name="Del_Chi_Out" )
    call output ( cs%init_angle, before='Init angle: ', advance='yes' )
    call dump ( cs%p, name="P" )
    call dump ( cs%dp, name="dP" )

  end subroutine Dump_Convolve_Support

  ! -----------------------------------------  FOV_Convolve_Setup  -----
  subroutine FOV_Convolve_Setup ( AntennaPattern, Chi_in, Chi_out, &
    & Convolve_Support, R_eq, R_sc, Earth_frac, Do_dRad_dx, Do_Scan_Avg )

    use Allocate_Deallocate, only: Allocate_test
    use AntennaPatterns_m, only: AntennaPattern_T
    use MLSKinds, only: Rp
    use MLSNumerics, only: InterpolateArraySetup

    ! Inputs
    type(antennapattern_t), intent(in), target :: AntennaPattern

    real(rp), intent(in) :: Chi_in(:)  ! input pointing angles radians
    real(rp), intent(in) :: Chi_out(:) ! output pointing angles radians

    real(rp), intent(in) :: r_eq ! equivalent earth radius
    real(rp), intent(in) :: r_sc ! spacecraft radius from equivalent center
    real(rp), optional, intent(in) :: earth_frac ! fraction of earth in total
    !                                   filled-out pattern
    ! req, rsc and earth_frac are non critical parameters and don't
    ! really need to be supplied externally. They are used to partition the
    ! full FFT field between earth and space components.

    logical, optional, intent(in) :: Do_dRad_dx ! "dRad_dx will be convolved"
    logical, optional, intent(in) :: Do_Scan_Avg ! "Scan averaging will be done."

    ! Output
    type(Convolve_Support_T), intent(out) :: Convolve_Support

    integer :: AAAPN, I
    real(r8) :: AAAP_step, Ang_step
    real(r8) :: E_frac, R_ratio
    logical :: My_Do_dRad_dx, My_Do_Scan_Avg

    e_frac = 0.185
    if ( present(earth_frac)) e_frac = 0.5 * earth_frac
    r_ratio = r_eq / r_sc ! earth radius / spacecraft radius

    my_Do_dRad_dx = .false.
    if ( present(do_dRad_dx) ) my_Do_dRad_dx = do_DRad_dx
    my_Do_Scan_Avg = .false.
    if ( present(do_Scan_Avg) ) my_Do_Scan_Avg = do_Scan_Avg
    ! load up the antenna pattern

    convolve_support%antennaPattern => antennaPattern
    aaap_step = antennaPattern%lambda
    ang_step = 1.0_r8 / (no_fft * aaap_step)
    aaapn = min(no_fft,size(antennaPattern%aaap))
    convolve_support%p(1:aaapn) = antennaPattern%aaap(1:aaapn)
    convolve_support%p(aaapn+1:) = 0.0_r8
    ! p are really complex numbers masquerading as real ones

    ! construct the angles

    convolve_support%angles = (/(i*ang_step,i=-ffth,no_fft-ffth-1)/)
    convolve_support%init_angle = &
      & asin((r_ratio - e_frac*sqrt(1.0_r8-r_ratio**2)/aaap_step))

    call allocate_test ( convolve_support%del_chi_in, size(chi_in), 'Del_Chi_In', moduleName )
    call allocate_test ( convolve_support%del_chi_out, size(chi_out), 'Del_Chi_In', moduleName )

    convolve_support%del_chi_in = chi_in - convolve_support%init_angle
    convolve_support%del_chi_out = chi_out - convolve_support%init_angle

    ! set up for interpolations

    call interpolateArraySetup ( convolve_support%del_chi_in, &
      & convolve_support%angles(ffth:no_fft), &
      & METHOD='S', coeffs=convolve_support%coeffs_1, EXTRAPOLATE='C' )
    call interpolateArraySetup ( convolve_support%angles(ffth-1:no_fft-1), &
      & convolve_support%del_chi_out, &
      & METHOD='S', coeffs=convolve_support%coeffs_2, EXTRAPOLATE='C', &
      & dyByDx=my_Do_dRad_dx .and. .not. my_Do_Scan_Avg, intYdX = do_Scan_Avg )

  end subroutine FOV_Convolve_Setup

  ! --------------------------------------------  FOV_Convolve_1D  -----
  subroutine FOV_Convolve_1D ( Convolve_Support, &
    & Rad_In, MIF_Times, DeadTime, Rad_Out, DRad_dx_out, Rad_FFT_out )

    use MLSKinds, only: Rp, Rv
    use MLSNumerics, only: InterpolateValues
    use ScanAverage_m, only: ScanAverage

    type(convolve_support_t), intent(in) :: Convolve_Support
    real(rp), intent(in) :: Rad_In(:)   ! input radiances, I
    real(rv), pointer :: MIF_Times(:)   ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(out) :: rad_out(:) ! output radiances, IA = I*G
    real(rp), optional, intent(out) :: drad_dx_out(:) ! output derivative
!                                         of radiance wrt to Chi_out
    real(r8), optional, intent(out), target :: Rad_FFT_out(:) ! FFT(I).  Temp derivs need it

    integer :: I, J
    real(r8), dimension(ffth+1), target :: My_rad_FFT
    real(r8), dimension(no_fft) :: rad_fft1
    real(r8), dimension(:), pointer :: rad_fft

    rad_fft => my_rad_fft(1:ffth+1)
    if ( present(rad_fft_out) ) rad_fft => rad_fft_out(1:ffth+1)

    ! interpolate from input grid to FFT angles

    call interpolateValues ( convolve_support%coeffs_1, &
      & convolve_support%del_chi_in, rad_in, &
      & convolve_support%angles(ffth:no_fft), rad_fft, &
      & METHOD='S', EXTRAPOLATE='C' )

    ! might need this if the signs coming out wrong is a problem

    rad_fft = rad_fft(ffth+1:1:-1)

    ! take cosine transform of interpolated input array
    ! For DTCST the coefficients come out in order, and they're
    ! twice the coefficients from DRFT1

    ! rad_fft = FFT(I)
    call dtcst_t ( rad_fft, 'a' )
    rad_fft = 0.5 * rad_fft

    ! apply convolution theorem

    ! p = FFT(G), where G is the convolution kernel.
    ! rad_fft1 = FFT(I) FFT(G)

    ! Handle first and last coefficients first
    rad_fft1(1) = rad_fft(1) * convolve_support%p(1)
    rad_fft1(2) = rad_fft(ffth+1) * convolve_support%p(2)
    j = 1
    do i = 3, no_fft - 1, 2
      j = j + 1
      rad_fft1(i)   = rad_fft(j) * convolve_support%p(i)
      rad_fft1(i+1) = rad_fft(j) * convolve_support%p(i+1)
    end do

    call drft1_t ( rad_fft1, 's' )

    ! interpolate from FFT angles to output grid

    if ( associated(MIF_Times) ) then
      call scanAverage ( MIF_Times, deadTime(1,1), &
        & real(convolve_support%angles(ffth-1:no_fft-1),rp), &
        & convolve_support%del_chi_out, real(rad_fft1(ffth:no_fft),rp), rad_out )
    else
      call interpolateValues ( convolve_support%coeffs_2, &
        & convolve_support%angles(ffth-1:no_fft-1), &
        & rad_fft1(ffth:no_fft), convolve_support%del_chi_out, rad_out, &
        & METHOD='S', EXTRAPOLATE='C', dyByDx=drad_dx_out )
    end if

  end subroutine FOV_Convolve_1D

  ! --------------------------------------------  FOV_Convolve_2d  -----
  subroutine FOV_Convolve_2d ( Convolve_Support, dI_df, MIF_Times, DeadTime, &
    & dI_df_flag, dRad_df_out )

    use MLSKinds, only: Rp, Rv

    ! Inputs

    type(convolve_support_t), intent(in) :: Convolve_Support
    real(rp), intent(in) :: di_df(:,:) ! mixing ratio derivatives or any
    !                                 parameter where a simple convolution
    !                                 will suffice
    real(rv), pointer :: MIF_Times(:)   ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    logical, optional, intent(in) :: di_df_flag(:) ! Flag to indicate which of
    !                                 di_df are to be calculated.  Deafult true.

    ! outputs

    real(rp), intent(out) :: drad_df_out(:,:) ! output radiance
    !                                 derivatives for input di_df.

    integer :: I, N_Coeffs

    ! nominally the mixing ratio derivatives but can be used for any
    ! quantity requiring a simple convolution.

    n_coeffs = size(di_df,dim=2)

    do i = 1, n_coeffs
      if ( present(di_df_flag) ) then
        if ( .not. di_df_flag(i) ) then
          drad_df_out(:,i) = 0.0
          cycle
        end if
      end if

      call fov_convolve_1d ( convolve_support, di_df(:,i), MIF_Times, DeadTime, &
        & drad_df_out(:,i) )

    end do

  end subroutine FOV_Convolve_2d

  ! -----------------------------------  FOV_Convolve_Temp_Derivs  -----
  subroutine FOV_Convolve_Temp_Derivs ( Convolve_Support, Rad_In, &
    & Rad_FFT, Surf_Angle, MIF_Times, DeadTime, dI_dT, dx_dT, ddx_dxdT, &
    & dx_dT_out, di_dT_flag, dRad_dT_out )

    use MLSKinds, only: Rp, Rv, R8
    use MLSNumerics, only: Coefficients, Hunt, &
      & InterpolateArraySetup, InterpolateArrayTeardown, InterpolateValues
    use ScanAverage_m, only: ScanAverage

    ! inputs

    type(convolve_support_t), intent(in) :: Convolve_Support

    real(rp), intent(in) :: Rad_In(:)  ! input radiances
    real(r8), intent(in) :: Rad_FFT(:) ! Convolved radiances on FFT grid
    real(rp), intent(in) :: Surf_Angle ! An angle (radians) that defines the
    !                       Earth surface.
    real(rv), pointer :: MIF_Times(:)  ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:) ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                       temperature on chi_in
    real(rp), intent(in) :: dx_dT(:,:) ! derivative of angle wrt
    !                       temperature on chi_in
    real(rp), intent(in) :: ddx_dxdT(:,:) ! 2nd derivative wrt angle and
    !                       temperature on chi_in
    real(rp), intent(in) :: dx_dT_out(:,:) ! derivative of angle wrt
    !                       temperature on chi_out
    logical, optional, intent(in) :: dI_dT_flag(:) ! Indicates whether to
    !                       accumulate di_dT.  Assumed true if absent.

    ! outputs

    real(rp), optional, intent(out) :: dRad_dT_out(:,:) ! output radiance
    !                       derivatives wrt temperature.

    type(coefficients(r8)) :: Coeffs_t ! for chi_in-init_angle -> angles(ffth+zero_out_s+1:no_fft)
    integer :: AAAPN, I, J, K, N_Coeffs, Zero_out_s, Zero_out_t
    real(r8), dimension(no_fft) :: dp, rad_fft1, rad_fft2, rad_fft3
    real(r8) :: drad_dT_temp(size(convolve_support%del_chi_out))

    ! Set up for interpolations.  First find the surface dimension
    call hunt ( convolve_support%angles(ffth:no_fft), &
      & surf_angle-convolve_support%init_angle, zero_out_s )
    call hunt ( convolve_support%angles(ffth:no_fft), &
      & convolve_support%del_chi_in(SIZE(convolve_support%del_chi_in)), &
      & zero_out_t )
    call interpolateArraySetup ( convolve_support%del_chi_in, &
      & convolve_support%angles(ffth+zero_out_s+1:no_fft), &
      & METHOD='S', coeffs=coeffs_t, EXTRAPOLATE='C' )

    ! temperature derivatives calculation
    ! compute the antenna derivative function

    n_coeffs = size(di_dT,dim=2)

    ! third term first (its fft is coefficient independent)
    ! apply convolution theorem

    aaapn = min(no_fft,size(convolve_support%antennaPattern%aaap))
    ! Derivative of antenna pattern
    dp(1:aaapn) = convolve_support%antennaPattern%d1aap(1:aaapn)
    dp(aaapn+1:) = 0.0_r8
    ! dp are really complex numbers masquerading as real ones

    rad_fft1(1:2) = 0.0_rp
    j = 1
    do i = 3, no_fft-1, 2
      j = j + 1
      rad_fft1(i)   = rad_fft(j) * dp(i)
      rad_fft1(i+1) = rad_fft(j) * dp(i+1)
    end do

    call drft1_t ( rad_fft1, 's' )

    ! interpolate to output grid

    call interpolateValues ( convolve_support%coeffs_2, &
      & convolve_support%angles(ffth-1:no_fft-1), &
      & rad_fft1(ffth:no_fft), convolve_support%del_chi_out, drad_dT_temp, &
      & METHOD='S', EXTRAPOLATE='C' )

    do i = 1, n_coeffs

      if ( present(di_dT_flag) ) then
        if ( .not. di_dT_flag(i) ) cycle
      end if

    ! estimate the error compensation

      k = maxloc(ddx_dxdT(:,i),1)

    ! do rad_in * ddx_dxdT piece

      call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
        &  (rad_in-rad_in(k)) * ddx_dxdT(:,i), &
        &  convolve_support%angles(ffth+zero_out_s+1:no_fft), &
        &  rad_fft2(ffth+zero_out_s+1:no_fft), METHOD='S',  &
        &  EXTRAPOLATE='C' )

    ! zero out the subsurface stuff

      rad_fft2(ffth:ffth+zero_out_s) = 0.0_rp

    ! add in di_dT part

      call interpolateValues ( convolve_support%coeffs_1, &
        & convolve_support%del_chi_in, di_dT(:, i), &
        & convolve_support%angles(ffth:no_fft), rad_fft1(ffth:no_fft), &
        & METHOD='S', EXTRAPOLATE='C' )

      rad_fft2(ffth:no_fft) = rad_fft2(ffth:no_fft) + rad_fft1(ffth:no_fft)

    ! zero out this array above toa

      rad_fft2(ffth+zero_out_t + 1:no_fft) = 0.0_rp

      if ( any(rad_fft2(ffth+zero_out_s+1:ffth+zero_out_t) /= 0.0) ) then

    ! resymetrize

        rad_fft2(1:ffth-1) = rad_fft2(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the radiances
    ! identically to the prototype code

        rad_fft2 = cshift(rad_fft2,-1)

    ! take cosine transform of rad_in * ddx_dxdT + di_dT array
    ! Coefficients from DTCST come out in order, and are twice the
    ! coefficients from DRFT1

        call dtcst_t ( rad_fft2(1:ffth+1), 'a' )

      else

        rad_fft2(1:ffth+1) = 0.0

      end if

    ! do the rad_in * dx_dT term

      call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
        &   (rad_in-rad_in(k))*dx_dT(:,i), &
        &   convolve_support%angles(ffth+zero_out_s+1:no_fft), &
        &   rad_fft1(ffth+zero_out_s+1:no_fft), METHOD='S', EXTRAPOLATE='C' )

    ! zero out array below surf_angle

      rad_fft1(ffth:ffth+zero_out_s) = 0.0_rp

    ! resymetrize

      rad_fft1(1:ffth-1) = rad_fft1(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the radiances
    ! identically to the prototype code

      rad_fft1 = cshift(rad_fft1,-1)

    ! take fft of rad_in * ddx_dxdT + di_dT + rad_in * dx_dT array

      call dtcst_t ( rad_fft1(1:ffth+1), 'a' )

    ! Rearrange rad_fft2 from dtcst order to drft1 order


    ! apply convolution theorem

      rad_fft3(1) = 0.5 * rad_fft2(1) * convolve_support%p(1)
      rad_fft3(2) = 0.5 * rad_fft2(ffth+1) * convolve_support%p(2)
      k = 1
      do j = 3, no_fft-1, 2
        k = k + 1
        rad_fft3(j+1) = 0.5 * ( rad_fft2(k) * convolve_support%p(j+1) - rad_fft1(k) * dp(j+1) )
        rad_fft3(j)   = 0.5 * ( rad_fft2(k) * convolve_support%p(j)   - rad_fft1(k) * dp(j)   )
      end do

    ! interplolate to chi_out

      call drft1_t ( rad_fft3, 's' )

      if ( associated(MIF_Times) ) then
        call scanAverage ( MIF_Times, deadTime(1,1), &
          & real(convolve_support%angles(ffth-1:no_fft-1),rp), &
          & convolve_support%del_chi_out, real(rad_fft3(ffth:no_fft),rp), &
          & drad_dT_out(:, i) )
      else
        call interpolateValues ( convolve_support%coeffs_2, &
          & convolve_support%angles(ffth-1:no_fft-1), &
          & rad_fft3(ffth:no_fft), convolve_support%del_chi_out, drad_dT_out(:, i), &
          & METHOD='S', EXTRAPOLATE='C' )
      end if

    ! compute final result

      drad_dT_out(:,i) = drad_dT_out(:,i) + dx_dT_out(:,i)*drad_dT_temp

    end do               ! On i = 1, n_coeffs

    call interpolateArrayTeardown ( coeffs_t )

  end subroutine FOV_Convolve_Temp_Derivs

  ! -----------------------------------  FOV_Convolve_Temp_Derivs_Normalization  -----
  ! FOV Convolution of Temperature derivatives (with normalization)
  subroutine FOV_Convolve_Temp_Derivs_Normalization ( Convolve_Support, Rad_Diff, &
    & Rad_Diff_FFT, Surf_Angle, MIF_Times, DeadTime, dI_dT, dx_dT, ddx_dxdT, &
    & dx_dT_out, di_dT_flag, dRad_dT_out )

    use MLSKinds, only: Rp, Rv, R8
!   use output_m, only: outputNamedValue                       ! IGOR
    use MLSNumerics, only: Coefficients, Hunt, &
      & InterpolateArraySetup, InterpolateArrayTeardown, InterpolateValues
    use ScanAverage_m, only: ScanAverage

    ! inputs

    type(convolve_support_t), intent(in) :: Convolve_Support

    real(rp), intent(in) :: Rad_Diff(:)  ! input radiances differences I-IA
    real(r8), intent(in) :: Rad_Diff_FFT(:) ! FFT(I-IA)
    real(rp), intent(in) :: Surf_Angle ! An angle (radians) that defines the
    !                       Earth surface.
    real(rv), pointer :: MIF_Times(:)  ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:) ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                       temperature on chi_in
    real(rp), intent(in) :: dx_dT(:,:) ! derivative of angle wrt
    !                       temperature on chi_in
    real(rp), intent(in) :: ddx_dxdT(:,:) ! 2nd derivative wrt angle and
    !                       temperature on chi_in
    real(rp), intent(in) :: dx_dT_out(:,:) ! derivative of angle wrt
    !                       temperature on chi_out
    logical, optional, intent(in) :: dI_dT_flag(:) ! Indicates whether to
    !                       accumulate di_dT.  Assumed true if absent.

    ! outputs

    real(rp), optional, intent(out) :: dRad_dT_out(:,:) ! output radiance
    !                       derivatives wrt temperature.

    type(coefficients(r8)) :: Coeffs_t ! for chi_in-init_angle -> angles(ffth+zero_out_s+1:no_fft)
    integer :: AAAPN, I, J, K, N_Coeffs, Zero_out_s, Zero_out_t
    real(r8), dimension(no_fft) :: dp, rad_fft1, rad_fft2, rad_fft3
    real(r8) :: drad_dT_temp(size(convolve_support%del_chi_out))


    ! Set up for interpolations.  First find the surface dimension
    call hunt ( convolve_support%angles(ffth:no_fft), &
      & surf_angle-convolve_support%init_angle, zero_out_s )
    call hunt ( convolve_support%angles(ffth:no_fft), &
      & convolve_support%del_chi_in(SIZE(convolve_support%del_chi_in)), &
      & zero_out_t )
    call interpolateArraySetup ( convolve_support%del_chi_in, &
      & convolve_support%angles(ffth+zero_out_s+1:no_fft), &
      & METHOD='S', coeffs=coeffs_t, EXTRAPOLATE='C' )

    ! temperature derivatives calculation
    ! compute the antenna derivative function

    n_coeffs = size(di_dT,dim=2)

    ! third term first (its fft is coefficient independent)
    ! apply convolution theorem

    aaapn = min(no_fft,size(convolve_support%antennaPattern%aaap))
    ! Derivative of antenna pattern
    dp(1:aaapn) = convolve_support%antennaPattern%d1aap(1:aaapn)
    dp(aaapn+1:) = 0.0_r8
    ! dp are really complex numbers masquerading as real ones

    
    !{  dx$\_$dT $= \displaystyle \frac{d\epsilon}{dT} = \frac{\tan \epsilon}{H} \frac{dH}{dT}$, \\
    ! dx$\_$dT$\_$out $= \displaystyle \frac{d\epsilon_t}{dT} = \frac{\tan \epsilon_t}{H_t} \frac{dH_t}{dT}$, \\
    ! ddx$\_$dxdT $= \displaystyle \frac{d^2 \epsilon}{d\epsilon dT} = \frac{2 +
    ! \tan^2 \epsilon}{H} \frac{dH}{dT} + \frac{\eta}{T}$.

    !{ Only consider rad$\_$diff$\_$fft $(= I-I^A)$, instead of rad$\_$in $(=I)$, throughout. \\
    !  p $= FT(G(\epsilon))$, \ \ dp $= FT \left( \displaystyle \frac{dG(\epsilon)}{d\epsilon} \right)$, \\
    !  rad$\_$fft1 $= FT(I-I^A) \cdot FT \left( \displaystyle \frac{dG(\epsilon)}{d\epsilon} \right)$.

    rad_fft1(1:2) = 0.0_rp
    j = 1
    do i = 3, no_fft-1, 2
      j = j + 1
      rad_fft1(i)   = rad_diff_fft(j) * dp(i)
      rad_fft1(i+1) = rad_diff_fft(j) * dp(i+1)
    end do

    !{ rad$\_$fft1 = $IFT \left( FT(I-I^A) \cdot FT \left( \displaystyle \frac{dG(\epsilon)}{d\epsilon} \right) \right)$ 
    ! $= (I-I^A) \displaystyle \frac{dG(\epsilon)}{d\epsilon}$  \\
    !  Mode = 's' - Systhesis (i.e. ifft) \\
    !  drad_dT_temp $= (I-I^A) \displaystyle \frac{dG(\epsilon)}{d\epsilon}$   (interpolated)
   
    call drft1_t ( rad_fft1, 's' )

    ! interpolate to output grid: drad_dT_temp = interpolate(rad_fft1)

    call interpolateValues ( convolve_support%coeffs_2, &
      & convolve_support%angles(ffth-1:no_fft-1), &
      & rad_fft1(ffth:no_fft), convolve_support%del_chi_out, drad_dT_temp, &
      & METHOD='S', EXTRAPOLATE='C' )

    do i = 1, n_coeffs

      if ( present(di_dT_flag) ) then
        if ( .not. di_dT_flag(i) ) cycle
      end if

    ! estimate the error compensation

      k = maxloc(ddx_dxdT(:,i),1)

    !{ rad$\_$fft2 = $(I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$   (interpolated)

     ! call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
     !   &  (rad_in-rad_in(k)) * ddx_dxdT(:,i), &
     !   &  convolve_support%angles(ffth+zero_out_s+1:no_fft), &
     !   &  rad_fft2(ffth+zero_out_s+1:no_fft), METHOD='S',  &
     !   &  EXTRAPOLATE='C' )			                  ! IGOR

      call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
        &  rad_diff * ddx_dxdT(:,i), &
        &  convolve_support%angles(ffth+zero_out_s+1:no_fft), &
        &  rad_fft2(ffth+zero_out_s+1:no_fft), METHOD='S',  &
        &  EXTRAPOLATE='C' )                                      ! IGOR

    ! zero out the subsurface stuff

      rad_fft2(ffth:ffth+zero_out_s) = 0.0_rp

    !{ rad$\_$fft1 = $\displaystyle \frac{\partial I}{\partial T}$   (interpolated)\\
    !  rad$\_$fft2 = $(I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$ 
    ! $+ \displaystyle \frac{\partial I}{\partial T}$

      call interpolateValues ( convolve_support%coeffs_1, &
        & convolve_support%del_chi_in, di_dT(:, i), &
        & convolve_support%angles(ffth:no_fft), rad_fft1(ffth:no_fft), &
        & METHOD='S', EXTRAPOLATE='C' )

      rad_fft2(ffth:no_fft) = rad_fft2(ffth:no_fft) + rad_fft1(ffth:no_fft)

    ! zero out this array above toa

      rad_fft2(ffth+zero_out_t + 1:no_fft) = 0.0_rp

      if ( any(rad_fft2(ffth+zero_out_s+1:ffth+zero_out_t) /= 0.0) ) then

    ! resymetrize

        rad_fft2(1:ffth-1) = rad_fft2(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the radiances
    ! identically to the prototype code

        rad_fft2 = cshift(rad_fft2,-1)

    !{ rad$\_$fft2 $= FT \bigg( (I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$
    ! $+ \displaystyle \frac{\partial I}{\partial T} \bigg)$

    ! take cosine transform of rad_diff * ddx_dxdT + di_dT array        ! IGOR
    ! Coefficients from DTCST come out in order, and are twice the
    ! coefficients from DRFT1

        call dtcst_t ( rad_fft2(1:ffth+1), 'a' )

      else

        rad_fft2(1:ffth+1) = 0.0

      end if

    !{ rad$\_$fft1 = $(I-I^A) \displaystyle \frac{d\epsilon}{dT}$   (interpolated)

!      call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
!        &   (rad_in-rad_in(k))*dx_dT(:,i), &
!        &   convolve_support%angles(ffth+zero_out_s+1:no_fft), &
!        &   rad_fft1(ffth+zero_out_s+1:no_fft), METHOD='S', EXTRAPOLATE='C' )    ! IGOR

                                                                      ! IGOR
      call interpolateValues ( coeffs_t, convolve_support%del_chi_in, &
        &   rad_diff*dx_dT(:,i), &
        &   convolve_support%angles(ffth+zero_out_s+1:no_fft), &
        &   rad_fft1(ffth+zero_out_s+1:no_fft), METHOD='S', EXTRAPOLATE='C' )

    ! zero out array below surf_angle

      rad_fft1(ffth:ffth+zero_out_s) = 0.0_rp

    ! resymetrize

      rad_fft1(1:ffth-1) = rad_fft1(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the radiances
    ! identically to the prototype code

      rad_fft1 = cshift(rad_fft1,-1)

    ! take fft of rad_diff * ddx_dxdT + di_dT + rad_diff * dx_dT array     ! IGOR

    !{ rad$\_$fft1 $= FT\left( (I-I^A) \displaystyle \frac{d\epsilon}{dT} \right)$
 
      call dtcst_t ( rad_fft1(1:ffth+1), 'a' )

    ! Rearrange rad_fft2 from dtcst order to drft1 order


    ! apply convolution theorem

    !{ rad$\_$fft3 $= FT \bigg( (I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$
    ! $+ \displaystyle \frac{\partial I}{\partial T} \bigg) \cdot FT\big(G(\epsilon)\big)$
    ! $- FT \left( (I-I^A) \displaystyle \frac{d\epsilon}{dT} \right) $
    ! $\cdot FT \left( \displaystyle \frac{dG(\epsilon)}{d\epsilon} \right)$

      rad_fft3(1) = 0.5 * rad_fft2(1) * convolve_support%p(1)
      rad_fft3(2) = 0.5 * rad_fft2(ffth+1) * convolve_support%p(2)
      k = 1
      do j = 3, no_fft-1, 2
        k = k + 1
        rad_fft3(j+1) = 0.5 * ( rad_fft2(k) * convolve_support%p(j+1) - rad_fft1(k) * dp(j+1) )
        rad_fft3(j)   = 0.5 * ( rad_fft2(k) * convolve_support%p(j)   - rad_fft1(k) * dp(j)   )
      end do

    ! interplolate to chi_out

    !{ rad$\_$fft3 $= \bigg( (I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$
    ! $+ \displaystyle \frac{\partial I}{\partial T} \bigg) \cdot G(\epsilon)$
    ! $- (I-I^A) \displaystyle \frac{d\epsilon}{dT} $
    ! $\cdot \displaystyle \frac{dG(\epsilon)}{d\epsilon}$

      call drft1_t ( rad_fft3, 's' )

    !{ drad$\_$dT$\_$out $= \bigg( (I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$
    ! $+ \displaystyle \frac{\partial I}{\partial T} \bigg) \cdot G(\epsilon)$
    ! $- (I-I^A) \displaystyle \frac{d\epsilon}{dT} $
    ! $\cdot \displaystyle \frac{dG(\epsilon)}{d\epsilon}$   (interpolated)

      if ( associated(MIF_Times) ) then
        call scanAverage ( MIF_Times, deadTime(1,1), &
          & real(convolve_support%angles(ffth-1:no_fft-1),rp), &
          & convolve_support%del_chi_out, real(rad_fft3(ffth:no_fft),rp), &
          & drad_dT_out(:, i) )
      else
        call interpolateValues ( convolve_support%coeffs_2, &
          & convolve_support%angles(ffth-1:no_fft-1), &
          & rad_fft3(ffth:no_fft), convolve_support%del_chi_out, drad_dT_out(:, i), &
          & METHOD='S', EXTRAPOLATE='C' )
      end if

    ! compute final result

    !{ drad$\_$dT$\_$out $= \bigg( (I-I^A) \displaystyle \frac{d^2\epsilon}{d\epsilon dT}$
    ! $+ \displaystyle \frac{\partial I}{\partial T} \bigg) \cdot G(\epsilon)$
    ! $- (I-I^A) \displaystyle \frac{d\epsilon}{dT} $
    ! $\cdot \displaystyle \frac{dG(\epsilon)}{d\epsilon}$
    ! $+ \displaystyle \frac{d\epsilon_t}{dT} \cdot (I-I^A) \displaystyle \frac{dG(\epsilon)}{d\epsilon}$

      drad_dT_out(:,i) = drad_dT_out(:,i) + dx_dT_out(:,i)*drad_dT_temp

      ! IGOR
      !sum_p = sum(convolve_support%p)
      !call outputNamedValue( 'sum_p = sum(convolve_support%p)', sum_p )
      !call outputNamedValue( 'size(convolve_support%p)', size(convolve_support%p) )

    end do               ! On i = 1, n_coeffs

    call interpolateArrayTeardown ( coeffs_t )

  end subroutine FOV_Convolve_Temp_Derivs_Normalization

  ! --------------------------------------  FOV_Convolve_Teardown  -----
  subroutine FOV_Convolve_Teardown ( Convolve_Support )
  ! Destroy the stuff in Convolve_Support
    use Allocate_Deallocate, only: Deallocate_Test
    use MLSNumerics, only: InterpolateArrayTeardown

    type(convolve_support_t), intent(inout) :: Convolve_Support

    call deallocate_test ( convolve_support%del_chi_in, 'Del_Chi_In', moduleName )
    call deallocate_test ( convolve_support%del_chi_out, 'Del_Chi_In', moduleName )

    call interpolateArrayTeardown ( convolve_support%coeffs_1 )
    call interpolateArrayTeardown ( convolve_support%coeffs_2 )

  end subroutine FOV_Convolve_Teardown

! =====     Private procedures  ========================================

  subroutine DRFT1_T ( A, Mode )
  ! Call DRFT1 and test its status flag
    use DFFT_M, only: DRFT1
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use SineTables_m, only: CreateSineTable, DestroySineTable, &
      & LogSize_SineTable_R8, SineTable_R8
    real(r8) :: A(:)
    character :: Mode
    call createSineTable ( pwr - 2 )
    call drft1 ( a, mode, pwr, logSize_SineTable_R8, sineTable_R8 )
    if ( logSize_SineTable_R8 == -2 ) then
      call DestroySineTable
      call MLSMessage ( MLSMSG_Error, ModuleName, "Error in drft1" )
    end if
  end subroutine DRFT1_T

  subroutine DTCST_T ( A, Mode )
  ! Call DTCST and test its status flag
    use DFFT_M, only: DTCST
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use SineTables_m, only: CreateSineTable, DestroySineTable, &
      & LogSize_SineTable_R8, SineTable_R8
    real(r8), intent(inout) :: A(:)
    character, intent(in) :: Mode
    integer, parameter :: DTCST_SIZE(1) = (/ pwr - 1 /) ! Avoid a run-time temp
    call createSineTable ( pwr - 2 )
    call dtcst ( a, 'C', mode, dtcst_size, 1, logSize_SineTable_R8, sineTable_R8 )
    if ( logSize_SineTable_R8 == -2 ) then
      call DestroySineTable
      call MLSMessage ( MLSMSG_Error, ModuleName, "Error in dtcst" )
    end if
  end subroutine DTCST_T

  ! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module FOV_Convolve_m

! $Log$
! Revision 2.19  2018/05/14 23:43:51  vsnyder
! Move Hessians stuff to Hessians_m
!
! Revision 2.18  2017/10/31 23:49:35  vsnyder
! Make Coefficients a parameterized type
!
! Revision 2.17  2013/06/12 02:19:37  vsnyder
! Cruft removal
!
! Revision 2.16  2012/07/07 00:14:33  vsnyder
! Shorten some comments to avoid gripes about long lines
!
! Revision 2.15  2012/07/06 21:30:00  yanovsky
! Added FOV_Convolve_Temp_Derivs_Normalization subroutine that computes
! normalized Temperature derivatives
!
! Revision 2.14  2011/03/23 23:49:20  vsnyder
! Make some array bounds explicit, to avoid a bounds violation if the
! actual argument is bigger than needed for the FFT.
!
! Revision 2.13  2011/03/23 23:45:32  vsnyder
! This log entry is bogus.  Check in again to get the right one.
! FOV_Convolve_m.f90
!
! Revision 2.12  2011/03/02 02:05:59  vsnyder
! Make AntennaPattern_T, Coefficients public, for F95 compatibility
!
! Revision 2.11  2010/07/19 18:21:09  yanovsky
! Add FOV_Convolve_3d subroutine
!
! Revision 2.10  2010/02/05 03:19:05  vsnyder
! Remove USE for unreferenced name
!
! Revision 2.9  2009/12/15 03:17:07  vsnyder
! Get kinds from MLSKinds instead of MLSCommon
!
! Revision 2.8  2009/11/17 23:38:11  vsnyder
! Add Dump, make R_eq, R_sc nonoptional in FOV_Convolve_Setup
!
! Revision 2.7  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2007/06/25 20:34:15  vsnyder
! Use DTCST instead of DFFT for even transforms
!
! Revision 2.5  2007/05/23 22:39:03  vsnyder
! Make sure drad_df_out gets defined
!
! Revision 2.4  2005/08/06 01:40:45  vsnyder
! ScanAverage doesn't need coeffs
!
! Revision 2.3  2005/08/03 18:03:20  vsnyder
! Scan averaging
!
! Revision 2.2  2005/07/08 00:12:11  vsnyder
! Get Rad_FFT from Convolve_Radiance to Convolve_Temperature_Deriv
!
! Revision 2.1  2005/07/06 02:16:34  vsnyder
! Initial commit, replacing fov_convolve_m.f90
!
