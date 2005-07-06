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
  use MLSCommon, only: R8
  use MLSNumerics, only: Coefficients => Coefficients_r8

  implicit NONE
  private
  public :: FOV_Convolve_Setup, FOV_Convolve_1D, FOV_Convolve_2d
  public :: FOV_Convolve_Temp_Derivs, FOV_Convolve_Teardown

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! FFT-related parameters
  integer, parameter, private :: pwr=12, no_fft=2**pwr, ffth = no_fft / 2

  ! Results of setting up for convolution
  type, public :: Convolve_Support_T
    real(r8), dimension(no_fft) :: angles ! Basis for FFT
    type(antennaPattern_t), pointer :: AntennaPattern ! For temperature derivs
    type(coefficients) :: Coeffs_1 ! for chi_in-init_angle -> angles(ffth:no_fft)
    type(coefficients) :: Coeffs_2 ! for angles(ffth:no_fft)-ang_step -> chi_out-init_angle
    real(r8), pointer :: del_chi_in(:) => NULL(), del_chi_out(:) => NULL()
    real(r8) :: Init_angle
    real(r8), dimension(no_fft) :: p, dp  ! From antenna pattern
  end type Convolve_Support_T

contains

  ! -----------------------------------------  FOV_Convolve_Setup  -----
  subroutine FOV_Convolve_Setup ( AntennaPattern, Chi_in, Chi_out, &
    & Convolve_Support, Req, Rsc, Earth_frac, Do_dRad_dx )

    use Allocate_Deallocate, only: Allocate_test
    use AntennaPatterns_m, only: AntennaPattern_T
    use MLSCommon, only: Rp
    use MLSNumerics, only: InterpolateArraySetup

    ! Inputs
    type(antennapattern_t), intent(in), target :: AntennaPattern

    real(rp), intent(in) :: Chi_in(:)  ! input pointing angles radians
    real(rp), intent(in) :: Chi_out(:) ! output pointing angles radians

    real(rp), optional, intent(in) :: req ! equivalent earth radius
    real(rp), optional, intent(in) :: rsc ! spacecraft radius
    real(rp), optional, intent(in) :: earth_frac ! fraction of earth in total
    !                                   filled-out pattern
    ! req, rsc and earth_frac are non critical parameters and don't
    ! really need to be supplied externally. They are used to partition the
    ! full FFT field between earth and space components.

    logical, optional, intent(in) :: Do_dRad_dx ! "dRad_dx will be convolved"

    ! Output
    type(Convolve_Support_T), intent(out) :: Convolve_Support

    integer :: AAAPN, I
    real(r8) :: AAAP_step, Ang_step
    real(r8) :: E_frac, R_eq, R_sc, R_ratio
    logical :: My_Do_dRad_dx

    r_eq = 6371.0_rp
    r_sc = r_eq + 705.0_rp
    e_frac = 0.185
    if ( present(req)) r_eq = req
    if ( present(rsc)) r_sc = rsc
    if ( present(earth_frac)) e_frac = 0.5 * earth_frac
    r_ratio = r_eq / r_sc ! earth radius / spacecraft radius

    my_Do_dRad_dx = .false.
    if ( present(do_dRad_dx) ) my_Do_dRad_dx = do_DRad_dx

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
      & dyByDx=my_Do_dRad_dx )

  end subroutine FOV_Convolve_Setup

  ! --------------------------------------------  FOV_Convolve_1D  -----
  subroutine FOV_Convolve_1D ( Convolve_Support, &
    & Rad_In, Rad_Out, DRad_dx_out )

    use MLSCommon, only: Rp
    use MLSNumerics, only: InterpolateValues

    type(convolve_support_t), intent(in) :: Convolve_Support
    real(rp), intent(in) :: Rad_In(:)   ! input radiances
    real(rp), intent(out) :: rad_out(:) ! output radiances
    real(rp), optional, intent(out) :: drad_dx_out(:) ! output derivative
!                                      of radiance wrt to Chi_out

    integer :: I
    real(r8), dimension(no_fft) :: rad_fft, rad_fft1

    ! interpolate from input grid to FFT angles

    call interpolateValues ( convolve_support%coeffs_1, &
      & convolve_support%del_chi_in, rad_in, &
      & convolve_support%angles(ffth:no_fft), rad_fft(ffth:no_fft), &
      & METHOD='S', EXTRAPOLATE='C' )

    ! mirror reflect this

    rad_fft(1:ffth-1) = rad_fft(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the
    ! radiances identically to the prototype code

    rad_fft = cshift(rad_fft,-1)

    ! take fft of interpolated input array

    call drft1_t ( rad_fft, 'a' )

    ! apply convolution theorem

    rad_fft1(1:2) = rad_fft(1:2) * convolve_support%p(1:2)
    do i = 3, no_fft - 1, 2
      rad_fft1(i)   = rad_fft(i) * convolve_support%p(i)
      rad_fft1(i+1) = rad_fft(i) * convolve_support%p(i+1)
    end do

    call drft1_t ( rad_fft1, 's' )

    ! interpolate from FFT angles to output grid

    call interpolateValues ( convolve_support%coeffs_2, &
      & convolve_support%angles(ffth-1:no_fft-1), &
      & rad_fft1(ffth:no_fft), convolve_support%del_chi_out, rad_out, &
      & METHOD='S', EXTRAPOLATE='C', dyByDx=drad_dx_out )

  end subroutine FOV_Convolve_1D

  ! --------------------------------------------  FOV_Convolve_2d  -----
  subroutine FOV_Convolve_2d ( Convolve_Support, dI_df, dI_df_flag, &
    & dRad_df_out )

    use MLSCommon, only: Rp

    ! Inputs

    type(convolve_support_t), intent(in) :: Convolve_Support
    real(rp), intent(in) :: di_df(:,:) ! mixing ratio derivatives or any
    !                                 parameter where a simple convolution
    !                                 will suffice
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
        if ( .not. di_df_flag(i) ) cycle
      end if

      call fov_convolve_1d ( convolve_support, di_df(:,i), drad_df_out(:,i) )

    end do

  end subroutine FOV_Convolve_2d

  ! -----------------------------------  FOV_Convolve_Temp_Derivs  -----
  subroutine FOV_Convolve_Temp_Derivs ( Convolve_Support, &
    & Rad_In, Surf_Angle, dI_dT, dx_dT, ddx_dxdT, dx_dT_out, di_dT_flag, &
    & dRad_dT_out )

    use MLSCommon, only: Rp
    use MLSNumerics, only: Coefficients=>Coefficients_r8, Hunt, &
      & InterpolateArraySetup, InterpolateArrayTeardown, InterpolateValues

    ! inputs

    type(convolve_support_t), intent(in) :: Convolve_Support

    real(rp), intent(in) :: Rad_In(:)  ! input radiances
    real(rp), intent(in) :: Surf_Angle ! An angle (radians) that defines the
    !                       Earth surface.
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

    type(coefficients) :: Coeffs_t ! for chi_in-init_angle -> angles(ffth+zero_out_s+1:no_fft)
    integer :: AAAPN, I, J, K, N_Coeffs, Zero_out_s, Zero_out_t
    real(r8), dimension(no_fft) :: dp, rad_fft, rad_fft1
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
    do i = 3, no_fft-1, 2
      rad_fft1(i)   = rad_fft(i) * dp(i)
      rad_fft1(i+1) = rad_fft(i) * dp(i+1)
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
        &  rad_fft(ffth+zero_out_s+1:no_fft), METHOD='S',  &
        &  EXTRAPOLATE='C' )

    ! zero out the subsurface stuff

      rad_fft(ffth:ffth+zero_out_s) = 0.0_rp


    ! add in di_dT part

      call interpolateValues ( convolve_support%coeffs_1, &
        & convolve_support%del_chi_in, di_dT(:, i), &
        & convolve_support%angles(ffth:no_fft), rad_fft1(ffth:no_fft), &
        & METHOD='S', EXTRAPOLATE='C' )

      rad_fft(ffth:no_fft) = rad_fft(ffth:no_fft) + rad_fft1(ffth:no_fft)

    ! zero out this array above toa

      rad_fft(ffth+zero_out_t + 1:no_fft) = 0.0_rp

    ! resymetrize

      rad_fft(1:ffth-1) = rad_fft(no_fft-1:no_fft-ffth+1:-1)

    ! I don't know if this step is truly necessary but it rephases the radiances
    ! identically to the prototype code

      rad_fft = cshift(rad_fft,-1)

    ! take fft of rad_in * ddx_dxdT + di_dT array

      call drft1_t ( rad_fft, 'a' )

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

      call drft1_t ( rad_fft1, 'a' )

    ! apply convolution theorem

      rad_fft(1:2) = rad_fft(1:2) * convolve_support%p(1:2)
      do j = 3, no_fft-1, 2
        rad_fft(j+1) = rad_fft(j) * convolve_support%p(j+1) - rad_fft1(j) * dp(j+1)
        rad_fft(j)   = rad_fft(j) * convolve_support%p(j)   - rad_fft1(j) * dp(j)
      end do

    ! interplolate to chi_out

      call drft1_t ( rad_fft, 's' )

      call interpolateValues ( convolve_support%coeffs_2, &
        & convolve_support%angles(ffth-1:no_fft-1), &
        & rad_fft(ffth:no_fft), convolve_support%del_chi_out, drad_dT_out(:, i), &
        & METHOD='S', EXTRAPOLATE='C' )

    ! compute final result

      drad_dT_out(:,i) = drad_dT_out(:,i) + dx_dT_out(:,i)*drad_dT_temp

    end do               ! On i = 1, n_coeffs

    call interpolateArrayTeardown ( coeffs_t )

  end subroutine FOV_Convolve_Temp_Derivs

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

  ! ----------------------------------------------  not_used_here  -----
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module FOV_Convolve_m

! $Log$
! Revision 2.1  2005/07/06 02:16:34  vsnyder
! Initial commit, replacing fov_convolve_m.f90
!
