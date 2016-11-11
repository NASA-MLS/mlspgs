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
module Convolution_m
!=============================================================================

  implicit NONE
  private

  public :: Convolution, Convolution_Setup

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ------------------------------------------------  Convolution  -----
  subroutine Convolution ( DH_DZ_OUT, DX_DH_OUT, DX_DT, DXDT_Surface, &
                         & DXDT_TAN, D2X_DXDT, &
                         & EarthRadC_sq, Est_ScGeocAlt, FirstSignal, FmStat, &
                         & FwdModelConf, FwdModelExtra, FwdModelIn, FwdModelOut, &
                         & Grids_f, Grids_n, Grids_tmp, Grids_v, Grids_w, &
                         & H_Atmos, K_Atmos, K_Spect_DN, K_Spect_DV, K_Spect_DW, &
                         & K_Temp, L1BMIF_TAI, MIFDeadTime, PTan, PTan_Der, &
                         & Ptg_Angles, Radiances, Sideband, S_T, Surf_Angle, &
                         & Tan_Chi_Out, Tan_Phi, Temp, ThisSideband, &
                         & Jacobian, ExtraJacobian, Hessian )

    ! Convolution if needed, or interpolation to ptan ----------------

    use ANTENNAPATTERNS_M, only: ANTENNAPATTERNS
    use INTRINSIC, only: L_ELEVOFFSET, L_LIMBSIDEBANDFRACTION
    use CONSTANTS, only: Deg2Rad
    use CONVOLVE_ALL_M, only: CONVOLVE_RADIANCE, CONVOLVE_TEMPERATURE_DERIV, &
      & CONVOLVE_OTHER_DERIV, CONVOLVE_OTHER_SECOND_DERIV, INTERPOLATE_RADIANCE, &
      & INTERPOLATE_TEMPERATURE_DERIV, INTERPOLATE_OTHER_DERIV
!     use CONVOLVE_ALL_M, only: CONVOLVE_TEMPERATURE_DERIV_NORMALIZATION, &
!       & CONVOLVE_RADIANCE_NORMALIZATION
    use DUMP_0, only: DUMP
    use ForwardModelConfig, only: Beta_Group_T, Channels_T, ForwardModelConfig_T
    use ForwardModelIntermediate, only: B_Ptg_Angles, ForwardModelStatus_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use FOV_CONVOLVE_M, only: CONVOLVE_SUPPORT_T, FOV_CONVOLVE_SETUP, &
      & FOV_CONVOLVE_TEARDOWN, NO_FFT
    use Geometry, only: Get_R_EQ
    use HessianModule_1, only: Hessian_T
    use Intrinsic, only: L_Radiance
    use Load_SPS_Data_M, only: Grids_T
    use MatrixModule_1, only: Matrix_T
    use MLSKinds, only: RP, RV, R4, R8
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_WARNING
    use MLSNUMERICS, only: Average, Coefficients => Coefficients_R8, &
      & InterpolateArraySetup, InterpolateArrayTeardown
    use MLSSignals_M, only: AreSignalsSuperset, Signal_T
    use MLSStringLists, only: SwitchDetail
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use VectorsModule, only: Vector_T, VectorValue_T

    real(rp), intent(in), contiguous :: DH_DZ_OUT(:) ! (ptan%template%nosurfs)
    real(rp), intent(in), contiguous :: DX_DH_OUT(:) ! (ptan%template%nosurfs)
    real(rp), intent(in), contiguous :: DX_DT(:,:) ! (no_tan_hts,s_t*sv_t_len)     ! (No_tan_hts, nz*np)
    real(rp), intent(in), contiguous :: DXDT_Surface(:,:) ! (1,s_t*sv_t_len)
    real(rp), intent(in), contiguous :: DXDT_TAN(:,:) ! (ptan%template%nosurfs,s_t*sv_t_len)
    real(rp), intent(in), contiguous :: D2X_DXDT(:,:) ! (no_tan_hts,s_t*sv_t_len)    ! (No_tan_hts, nz*np)
    type(vectorValue_t), intent(in) :: EarthRadC_sq  ! (minor axis of orbit plane projected Earth ellipse)**2
    real(rp), intent(in), contiguous :: Est_ScGeocAlt(:) ! Est S/C geocentric altitude /m
    type (Signal_T), intent(in)  :: FirstSignal ! The first signal we're dealing with
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type (Grids_T), intent(in) :: Grids_f   ! All the coordinates for VMR
    type (Grids_T), intent(in) :: Grids_n   ! All the spectroscopy(N) coordinates
    type (Grids_T), intent(in) :: Grids_tmp ! All the coordinates for TEMP
    type (Grids_T), intent(in) :: Grids_v   ! All the spectroscopy(V) coordinates
    type (Grids_T), intent(in) :: Grids_w   ! All the spectroscopy(W) coordinates
    real(r4), intent(in) :: H_Atmos(:,:,:,:) ! (noUsedChannels,no_tan_hts,s_h*size(grids_f%values),s_h*size(grids_f%values))
    real(r4), intent(in) :: K_Atmos(:,:,:) ! (noUsedChannels,no_tan_hts,s_a*size(grids_f%values))
    real(r4), intent(in) :: K_Spect_DN(:,:,:) ! (noUsedChannels,no_tan_hts,s_td*size(grids_n%values))
    real(r4), intent(in) :: K_Spect_DV(:,:,:) ! (noUsedChannels,no_tan_hts,s_lc*size(grids_v%values))
    real(r4), intent(in) :: K_Spect_DW(:,:,:) ! (noUsedChannels,no_tan_hts,s_lw*size(grids_w%values))
    real(r4), intent(in) :: K_Temp(:,:,:) ! (noUsedChannels,no_tan_hts,s_t*sv_t_len)
    real(rv), intent(in), pointer :: L1BMIF_TAI(:,:)   ! MIF Times
    real(rv), intent(in), pointer :: MIFDeadTime(:,:)  ! Not collecting data
    type (VectorValue_T), intent(in) :: PTan   ! Tangent pressure component of state vector
    logical, intent(in) :: PTan_Der       ! Compute derivative w.r.t. PTan
    real(rp), intent(inout), contiguous :: Ptg_Angles(:)
    real(rp), intent(in), contiguous :: Radiances(:,:) ! (noUsedChannels,no_tan_hts)
    integer, intent(in) :: Sideband       ! Either zero or from firstSignal
    integer, intent(in) :: S_T            ! Multiplier for temp derivative sizes, 0 or 1
    real(rp), intent(in) :: Surf_Angle(1)
    real(rp), intent(in), contiguous :: Tan_Chi_Out(:)
    real(rp), intent(in), contiguous :: Tan_Phi(:)
    type (VectorValue_T), intent(in) :: Temp   ! Temperature component of state vector
    integer, intent(in) :: ThisSideband   ! -1 = LSB, +1 = USB

    type(matrix_T), intent(inout), optional :: Jacobian
    type(matrix_T), intent(inout), optional :: ExtraJacobian
    type(hessian_T), intent(inout), optional :: Hessian

    logical :: Atmos_Der         ! Compute derivative w.r.t. VMR
    logical :: Atmos_Second_Der  ! Compute 2nd derivative w.r.t. VMR
    type (Beta_Group_T), pointer :: Beta_Group(:) ! from FwdModelConf%Beta_Group
    type (Channels_T), pointer :: Channels(:) ! from FwdModelConf%Channels
    integer ChanInd, Channel, I, J, SigInd, Ptg_I
    real(rp) :: DELTAPTG         ! Used for patching the pointings
    integer :: MAF
    integer :: Me = -1           ! String index for trace
    integer :: MINSUPERSET       ! Min. value of superset > 0
    integer :: No_Tan_Hts        ! Number of tangent heights
    integer :: NoUsedChannels    ! Number of channels used
    logical :: PATCHEDAPTG       ! Used in patching the pointings
    logical :: Print_Ptg         ! For debugging, from -Sptg
    integer :: PTG_J             ! Loop counters for patching the pointings
    real(r8) :: RAD_FFT(s_t*no_fft) ! FFT(I)             ! IGOR
!     real(r8) :: RAD_DIFF_FFT(s_t*no_fft) ! FFT(I-IA)   ! IGOR
    logical :: Spect_Der, Spect_Der_Center, Spect_Der_Width, Spect_Der_Width_TDep
    integer :: SUPERSET          ! Output from AreSignalsSuperset
    real(rp) :: THISELEV         ! An elevation offset
    real(rp) :: THISFRACTION     ! A sideband fraction
    logical :: Update            ! Just update radiances etc.
    integer :: WHICHPATTERN      ! Index of antenna pattern
    type (Coefficients) :: Coeffs ! For interpolation
    type (Convolve_Support_T) :: Convolve_Support
    type (VectorValue_T), pointer :: ELEVOFFSET       ! Elevation offset
    type (VectorValue_T), pointer :: SIDEBANDFRACTION ! The sideband fraction to use
    logical :: Temp_Der          ! Compute derivative w.r.t. Temp
    type (VectorValue_T), pointer :: THISRADIANCE     ! A radiance vector quantity

    call trace_begin ( me, 'ForwardModel.Convolution', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

    atmos_der = present ( jacobian ) .and. FwdModelConf%atmos_der
    atmos_second_der = present (hessian ) .and. FwdModelConf%atmos_second_der

    spect_der = present ( jacobian ) .and. FwdModelConf%spect_der
    spect_der_center = spect_der .and. size(fwdModelConf%lineCenter) > 0
    spect_der_width = spect_der .and. size(fwdModelConf%lineWidth) > 0
    spect_der_width_TDep = spect_der .and. size(fwdModelConf%lineWidth_TDep) > 0

    temp_der = present ( jacobian ) .and. FwdModelConf%temp_der

    beta_group => fwdModelConf%beta_group
    channels => fwdModelConf%channels
    MAF = fmstat%MAF
    no_tan_hts = size(ptg_angles)
    noUsedChannels = size(channels)

    ! Check that the angles are in the correct order.  If they
    ! are not it means (give or take some approximations in the
    ! horizontal according to Bill), that the rays crossed over
    ! between the tangent point and the spacecraft.  One could dream
    ! up all sorts of elegant schemes to get around that problem, but
    ! it's simplest just to bail out (and is certainly preferable to
    ! the infinite loop in the convolution (Hunt on angles) that
    ! results otherwise).

    print_Ptg = switchDetail(switches, 'ptg') > -1
    if ( print_Ptg ) &
      & call Dump ( ptg_angles, 'ptg_angles (before any patch)', format='(1PG22.17)' )

    ! This code is needed to ensure that the ptg_angles are monotonic
    ! (and not flat even)
    deltaPtg = 1e-3              ! Some starting value
    patchedAPtg = .false.
    do ptg_i = 2, no_tan_hts
      if ( ptg_angles(ptg_i) <= ptg_angles(ptg_i-1) ) then
        patchedAPtg = .true.
        ! This one is at or below its predecessor, find the next one above
        ! If there is no next one above just use the previous spacing
        do ptg_j = ptg_i + 1, no_tan_hts
          if ( ptg_angles(ptg_j) > ptg_angles(ptg_i-1) ) then
            ! Found one above.  Work out spacing to fill in with
            deltaPtg = ( ptg_angles(ptg_j) - ptg_angles(ptg_i-1) ) / ( ptg_j - ptg_i + 1 )
            exit
          end if
        end do
        ptg_angles(ptg_i) = ptg_angles(ptg_i-1) + deltaPtg
      else
        ! This value is above the previous one so compute a delta from it
        ! to use if needed later
        deltaPtg = ptg_angles(ptg_i) - ptg_angles(ptg_i-1)
      end if
    end do

    if ( patchedAPtg ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & 'Had to patch some out-of-order ptg_angles' )
      fmStat%flags = ior(fmStat%flags,b_ptg_angles)
      if ( print_Ptg ) &
        & call dump ( ptg_angles, 'ptg_angles (after patching)', format='(1PG22.17)' )
    end if

    ! Work out which antenna patterns we're going to need ------------------
    do i = 1, noUsedChannels
      channel = channels(i)%used
      chanInd = channel + 1 - channels(i)%origin
      sigInd = channels(i)%signal
      ! Get the radiance
      thisRadiance =>  &
        GetQuantityForForwardModel (fwdModelOut, quantityType=l_radiance, &
        & signal=fwdModelConf%signals(sigInd)%index, sideband=sideband,   &
        & config=fwdModelConf )
      ! Get the sideband fraction if we need to
      if ( firstSignal%sideband == 0 .or. fwdModelConf%forceSidebandFraction ) then
        sidebandFraction => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_limbSidebandFraction, &
          & signal=fwdModelConf%signals(sigInd)%index, &
          & sideband=thisSideband, config=fwdModelConf )
        thisFraction = sidebandFraction%values(chanInd,1)
      else                  ! Otherwise, want just unfolded signal
        thisFraction = 1.0
      end if
      ! Get the elevation offset
      elevOffset => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_elevOffset, signal=fwdModelConf%signals(sigInd)%index, &
        & sideband=thisSideband, config=fwdModelConf )
      thisElev = elevOffset%values(chanInd,1) * deg2Rad

      ! Here comes the Convolution (or not) codes
      update = ( thisSideband /= fwdModelConf%sidebandStart )

      if ( FwdModelConf%do_conv ) then

        whichPattern = -1
        minSuperset = huge(0)
        do j = 1, size(antennaPatterns)
          superset = AreSignalsSuperset ( antennaPatterns(j)%signals, &
            & fwdModelConf%signals(sigInd:sigInd), sideband=thisSideband, &
            & channel=channel )
          if ( superset >= 0 .and. superset <= minSuperset ) then
            minSuperset = superset
            whichPattern = j
          end if
        end do
        if ( whichPattern < 0 ) call Announce_Error ( &
          & "No matching antenna patterns." )

        call fov_convolve_setup ( antennaPatterns(whichPattern), ptg_angles, &
          & tan_chi_out-thisElev, convolve_support, &
          & get_r_eq(average(tan_phi), &
          &          average(earthradc_sq%values(:,MAF))), &  ! Average r_eq
          & average(0.001_rp*est_scgeocalt), &                ! Average alt
          & do_dRad_dx=ptan_der, do_Scan_Avg=fwdModelConf%scanAverage )

        if ( temp_der ) then

          ! To turn on/off Temperature Derivatives Normalization:
          ! comment/uncomment two function calls:

          call convolve_radiance ( convolve_support, maf, chanInd, &
            & radiances(i,:), thisFraction, update, ptan, thisRadiance, &
            & L1BMIF_TAI, MIFDeadTime, &
            & Jacobian, fmStat%rows, dh_dz_out, dx_dh_out, ptan_der, rad_FFT )    ! IGOR

          !call convolve_radiance_normalization ( convolve_support, maf, chanInd, &
          !  & radiances(i,:), thisFraction, update, ptan, thisRadiance, &
          !  & L1BMIF_TAI, MIFDeadTime, &
          !  & Jacobian, fmStat%rows, dh_dz_out, dx_dh_out, ptan_der, rad_FFT, &
          !  & radiances_diff(i,:), rad_diff_FFT )                                  ! IGOR

          call convolve_temperature_deriv ( convolve_support, maf, chanInd, &
            & radiances(i,:), rad_fft, thisFraction, update, thisRadiance, &
            & temp, grids_tmp, surf_angle(1), L1BMIF_TAI, MIFDeadTime, &
            & real(k_temp(i,:,:),kind=rp), &
            & dx_dt, d2x_dxdt, dxdt_tan, dxdt_surface, Jacobian, fmStat%rows )    ! IGOR

          !call convolve_temperature_deriv_normalization ( convolve_support, maf, chanInd, &
          !  & radiances_diff(i,:), rad_diff_FFT, thisFraction, update, thisRadiance, &
          !  & temp, grids_tmp, surf_angle(1), L1BMIF_TAI, MIFDeadTime, &
          !  & real(k_temp(i,:,:),kind=rp), &
          !  & dx_dt, d2x_dxdt, dxdt_tan, dxdt_surface, Jacobian, fmStat%rows )     ! IGOR

        else ! No temperature derivative
          call convolve_radiance ( convolve_support, maf, chanInd, &
            & radiances(i,:), thisFraction, update, ptan, thisRadiance, &
            & L1BMIF_TAI, MIFDeadTime, &
            & Jacobian, fmStat%rows, dh_dz_out, dx_dh_out, ptan_der )
        end if

        if ( atmos_der ) then
          call convolve_other_deriv ( convolve_support, maf, chanInd, &
          & thisFraction, update, thisRadiance, beta_group%qty, Grids_f, &
          & L1BMIF_TAI, MIFDeadTime, real(k_atmos(i,:,:),kind=rp), &
          & Jacobian, fmStat%rows, extraJacobian )

          if (atmos_second_der ) then
            call convolve_other_second_deriv ( convolve_support, maf, chanInd, &
            & thisFraction, update, thisRadiance, beta_group%qty, Grids_f, &
            & L1BMIF_TAI, MIFDeadTime, real(h_atmos(i,:,:,:),kind=rp), &
            & Hessian, fmStat%rows )
          end if
        end if

        if ( spect_der_center ) &
          & call convolve_other_deriv ( convolve_support, maf, chanInd, &
            & thisFraction, update, thisRadiance, &
            & fwdModelConf%lineCenter%qty, grids_v, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dv(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & extraJacobian )
        if ( spect_der_Width ) &
          & call convolve_other_deriv ( convolve_support, maf, chanInd, &
            & thisFraction, update, thisRadiance, &
            & fwdModelConf%lineWidth%qty, grids_w, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dw(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & extraJacobian )
        if ( spect_der_Width_TDep ) &
          & call convolve_other_deriv ( convolve_support, maf, chanInd, &
            & thisFraction, update, thisRadiance, &
            & fwdModelConf%lineWidth_TDep%qty, grids_n, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dn(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & extraJacobian )

        call fov_convolve_teardown ( convolve_support )

      else          ! No convolution needed ..

        call interpolateArraySetup ( ptg_angles, tan_chi_out-thisElev, &
          & method='S', extrapolate='C', coeffs=coeffs, &
          & dyByDx=ptan_der.or.fwdModelConf%scanAverage )

        call interpolate_radiance ( coeffs, maf, chanInd, ptg_angles, &
          & radiances(i,:), thisFraction, update, ptan, tan_chi_out-thisElev, &
          & thisRadiance, L1BMIF_TAI, MIFDeadTime, Jacobian, fmStat%rows, &
          & dh_dz_out, dx_dh_out, ptan_der )

        if ( temp_der ) &
          & call interpolate_temperature_deriv ( coeffs, maf, chanInd, &
            & ptg_angles, thisFraction, update, tan_chi_out-thisElev,  &
            & thisRadiance, temp, grids_tmp, L1BMIF_TAI, MIFDeadTime, &
            & real(k_temp(i,:,:),kind=rp), Jacobian, fmStat%rows )

        if ( atmos_der ) &
          call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
            & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
            & beta_group%qty, grids_f, L1BMIF_TAI, MIFDeadTime, &
            & real(k_atmos(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & linear=.true., extraJacobian=extraJacobian )

        if ( spect_der_center ) &
          & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
            & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
            & fwdModelConf%lineCenter%qty, grids_v, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dv(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & linear=.true., extraJacobian=extraJacobian )

        if ( spect_der_Width ) &
          & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
            & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
            & fwdModelConf%lineWidth%qty, grids_w, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dw(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & linear=.true., extraJacobian=extraJacobian )

        if ( spect_der_Width_TDep ) &
          & call interpolate_other_deriv ( coeffs, maf, chanInd, ptg_angles, &
            & thisFraction, update, tan_chi_out-thisElev, thisRadiance, &
            & fwdModelConf%lineWidth_TDep%qty, grids_n, L1BMIF_TAI, MIFDeadTime, &
            & real(k_spect_dn(i,:,:),kind=rp), Jacobian, fmStat%rows, &
            & linear=.true., extraJacobian=extraJacobian )

        call interpolateArrayTeardown ( coeffs )

      end if ! Convolve or interpolate

    end do                            ! Channel loop

    call trace_end ( 'ForwardModel.Convolution', &
      & cond=toggle(emit) .and. levels(emit) > 2 )

  contains

  ! .............................................  Announce_Error  .....
    subroutine Announce_Error ( Message )
    ! Announce Message using MLSMessage.  Include the configuration name
    ! in the message
      use MLSMessageModule, only: MLSMSG_Error
      use MoreMessage, only: MLSMessage
      character(len=*), intent(in) :: Message
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & "With config(%S): " // message, datum=fwdModelConf%name )
    end subroutine Announce_Error

  end subroutine Convolution

  ! ----------------------------------------  Convolution_Setup  -----
  subroutine Convolution_Setup ( DH_DZ_OUT, DX_DH_OUT, DXDT_Surface, &
                               & DXDT_TAN, EarthRadC_sq, Est_ScGeocAlt, &
                               & FwdModelConf, FwdModelExtra, FwdModelIn, &
                               & Grids_f, Grids_tmp, &
                               & L1BMIF_TAI, MAF, MIFDeadTime, &
                               & PhiTan, PTan, RefGPH, SCGeocAlt, &
                               & Surf_Angle, Tan_Chi_Out, Tan_Phi, TAN_Press, &
                               & WindowFinish, WindowStart )
  ! set up output pointing angles ------------------------------------

    use CONSTANTS, only: Deg2Rad
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Geometry, only: Get_R_EQ
    use Get_Chi_Out_m, only: Get_Chi_Out
    use Intrinsic, only: L_MIFDEADTIME, L_L1BMIF_TAI
    use Load_SPS_Data_M, only: Grids_T
    use MLSKinds, only: RP, RV
    use Molecules, only: L_H2O
    use Toggles, only: Emit, Levels, Toggle
    use Trace_M, only: Trace_Begin, Trace_End
    use VectorsModule, only: Vector_T, VectorValue_T

    real(rp), intent(out), contiguous :: DH_DZ_OUT(:) ! (ptan%template%nosurfs)
    real(rp), intent(out), contiguous :: DX_DH_OUT(:) ! (ptan%template%nosurfs)
    real(rp), intent(out), contiguous :: DXDT_Surface(:,:) ! (1,s_t*sv_t_len)
    real(rp), intent(out), contiguous :: DXDT_TAN(:,:) ! (ptan%template%nosurfs,s_t*sv_t_len)
    type(VectorValue_T), intent(in) :: EarthRadC_sq ! (minor axis of orbit plane projected Earth ellipse)**2
    real(rp), intent(in), contiguous :: Est_ScGeocAlt(:) ! (no_tan_hts) ! Est S/C geocentric altitude /m
    type(forwardModelConfig_T), intent(in) :: FwdModelConf
    type(vector_T), intent(in) :: FwdModelExtra, FwdModelIn
    type (Grids_T), intent(in) :: Grids_f      ! All the coordinates for VMR
    type (Grids_T), intent(in) :: Grids_tmp    ! All the coordinates for TEMP
    real(rv), intent(out), pointer :: L1BMIF_TAI(:,:)   ! MIF Times
    integer, intent(in) :: MAF
    real(rv), intent(out), pointer :: MIFDeadTime(:,:)  ! Not collecting data
    type (VectorValue_T), intent(in) :: PhiTan ! Tangent geodAngle component of state vector
    type (VectorValue_T), intent(in) :: PTan   ! Tangent pressure component of state vector
    type (VectorValue_T), intent(in) :: RefGPH ! Reference geopotential height
    type (VectorValue_T), intent(in) :: SCGeocAlt     ! S/C geocentric altitude /m
    real(rp), intent(out) :: Surf_Angle(1)
    real(rp), intent(out), contiguous :: Tan_Chi_Out(:)
    real(rp), intent(in), contiguous :: Tan_Phi(:) ! Radians
    real(rp), intent(in) :: TAN_Press(:)       ! Pressures corresponding to Z_PSIG
    integer, intent(in) :: WindowFinish        ! End of temperature `window'
    integer, intent(in) :: WindowStart         ! Start of temperature `window'

    real(rp) :: D2XDXDT_SURFACE(1,size(dxdt_surface,2))
    real(rp) :: D2XDXDT_TAN(size(dxdt_tan,1),size(dxdt_tan,2))
    integer :: H2O_Ind
    real(rp) :: One_dhdz(1), One_dxdh(1)
    real(rp) :: REQ_OUT(phitan%template%nosurfs)

    integer :: Me = -1                    ! String index for trace
    type (VectorValue_T), pointer :: WORK ! Temporary stuff

    call trace_begin ( me, 'ForwardModel.Convolution_Setup', &
      & cond=toggle(emit)  .and. levels(emit) > 0 )

    h2o_ind = grids_f%s_ind(l_h2o)

    ! Compute equivalent earth radius

    ! Although this is the same mathematical formula as used in metrics,
    ! the phi used here is different: These are on MIFs, not hypothetical
    ! pointings to the desired tangent zetas.  Therefore, we can't use
    ! these values in metrics, or where its output r_eq value is
    ! used.

    req_out = get_R_eq ( phitan%values(:,maf)*Deg2Rad, earthradc_sq%values(:,maf) )

    ! Temperature's windowStart:windowFinish are correct here.
    ! RefGPH and Temperature have the same horizontal basis.
    ! Grids_F is only needed for H2O, for calculating refractive index.
    ! SCgeocAlt and RefGPH are in meters, but Get_Chi_Out wants them in km.
    ! This is only used for convolution, which is done for both sidebands.
    call get_chi_out ( ptan%values(:,maf), phitan%values(:,maf)*deg2rad,  &
       & scGeocAlt%values(:,maf), Grids_tmp,                              &
       & spread(refGPH%template%surfs(1,1),1,windowFinish-windowStart+1), &
       & refGPH%values(1,windowStart:windowFinish), 0.0_rp,               &
       & req_out, grids_f, h2o_ind, tan_chi_out, dh_dz_out, dx_dh_out,    &
       & dxdt_tan=dxdt_tan, d2xdxdt_tan=d2xdxdt_tan )

    ! This is a lazy way to get the surface angle
    ! Temperature's windowStart:windowFinish are correct here.
    ! refGPH and temperature have the same horizontal basis.
    ! Grids_tmp is only needed for H2O, for calculating refractive index.
    ! Est_scgeocalt and RefGPH are in meters, but Get_Chi_Out wants them in km.
    ! This is only used for convolution, which is done for both sidebands.
    call get_chi_out ( tan_press(1:1), tan_phi(1:1),                      &
       & est_scgeocalt(1:1), Grids_tmp,                                   &
       & spread(refGPH%template%surfs(1,1),1,windowFinish-windowStart+1), &
       & refGPH%values(1,windowStart:windowFinish), 0.0_rp,               &
       & req_out(1:1), grids_f, h2o_ind, surf_angle, one_dhdz, one_dxdh,  &
       & dxdt_tan=dxdt_surface, d2xdxdt_tan=d2xdxdt_surface )

    ! This is only used for convolution, which is done for both sidebands.
    if ( fwdModelConf%scanAverage ) then
      work => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_l1bMIF_TAI, config=fwdModelConf )
      l1bMIF_TAI => work%values
      work => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_MIFDeadTime, config=fwdModelConf )
      MIFDeadTime => work%values ! Only the (1,1) element is used.
    else
      nullify ( l1bMIF_TAI, MIFDeadTime )
    end if

    call trace_end ( 'ForwardModel.Convolution_Setup', &
      & cond=toggle(emit) .and. levels(emit) > 0 )

  end subroutine Convolution_Setup

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Convolution_m

! $Log$
! Revision 2.8  2016/11/11 01:54:36  vsnyder
! Call Get_Chi_Out with ScGeocAlt and RefGPH in m instead of km, because
! Get_Chi_Out does the conversion now (to avoid an array temp).  Use
! SPREAD intrinsic function instead of an array constructor to create an
! array of all-the-same reference surface values.
!
! Revision 2.7  2016/10/24 22:13:42  vsnyder
! Eliminate orbit inclination because get_chi_out no longer needs it
!
! Revision 2.6  2016/08/20 00:53:48  vsnyder
! Correct a comment about units for Ten_Phi
!
! Revision 2.5  2016/07/28 00:44:55  vsnyder
! Remove unused variable declaration
!
! Revision 2.4  2016/06/03 23:49:50  vsnyder
! Change EarthRadC_sq from a real scalar to a vector quantity.  Use average
! of tangent phi / average EarthRacC_sq to do FOV_Convolve_Setup.
!
! Revision 2.3  2016/05/10 00:02:35  vsnyder
! Remove unused variable
!
! Revision 2.2  2016/05/02 23:32:15  vsnyder
! Remove unused dummy arguments from Convolution_Setup
!
! Revision 2.1  2016/04/21 01:59:12  vsnyder
! Initial commit
!
