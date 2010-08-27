! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module TScat_Support_m

  ! Provide types and procedures for use during calculation and use of TScat
  ! tables by the Full Forward Model.

  use MLSKinds, only: RP, R8
  implicit NONE
  private
  public :: Get_TScat, Interpolate_P_to_theta_e, TScat_Linear_t
  public :: Get_TScat_Setup, TScat_Gen_Setup

  ! Type for work area for Get_TScat.  The linear model is used to fill this
  ! at frequency Frq.  Then TScat (and its derivatives) are interpolated in
  ! frequency from one or two of these to the desired frequency.  These are
  ! filled here by Get_TScat and used here by Get_TScat, and not filled or
  ! used anywhere else.  The caller is expected, however, to allocate the
  ! array components.  Only the part germane to a particular path is used,
  ! so the caller can allocate them to the length of the longest coarse path.

  type :: TScat_Linear_t
    real(r8) :: Frq           ! Frequency at which TScat was evaluated
    logical :: Got = .false.  ! "there are data here"
    real(rp), pointer :: TScat_Path(:) => NULL()    ! TScat on the path
    real(rp), pointer :: K_Temp_TScat(:) => NULL()  ! dTScat / dT, T-SV length
    real(rp), pointer :: K_Atmos_TScat(:) => NULL() ! dTScat / dVMR, vmr-SV length
  end type

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  ! --------------------------------------------------  Get_TScat  -----
  subroutine Get_TScat ( TScat_Conf, Frq, Eta_T, Eta_F, Atmos_Der, Temp_Der, &
    &                    L2PC, xP, DeltaX, XstarPTAN,                        &
    &                    W, TScat_Path, K_Temp_TScat, K_Atmos_TScat )

  ! Get TScat and its derivatives along the integration path using
  ! the TScat tables and the same strategy as the linear forward model.

    use ForwardModelConfig, only: ForwardModelConfig_t
    use L2PC_m, only: L2PC_t
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: Hunt
    use Molecules, only: L_CLOUD_A
    use VectorsModule, only: VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(in) :: TScat_Conf
    real(r8), intent(in) :: Frq        ! Frequency at which to evaluate TScat
    real(rp), intent(in) :: Eta_T(:,:) ! Interpolation factors from temperature grid
                                       ! to path, path length X SV length
    real(rp), intent(in) :: Eta_F(:,:) ! Interpolation factors from VMR grid
                                       ! to path, path length X SV length
    logical, intent(in) :: Atmos_Der   ! Want dTScat / dIWC
    logical, intent(in) :: Temp_Der    ! Want dTScat / dT
    type(L2PC_t), intent(in) :: L2PC   ! The selected L2PC
    type(vector_T), intent(inout) :: XP     ! Same form as xStar, contents as x
    type(vector_T), intent(inout) :: DELTAX ! xp-xStar
    type(vectorValue_T), intent(inout) :: XSTARPTAN ! Tangent pressure in l2pc
    type(TScat_Linear_t), intent(inout) :: W(2) ! Work area, see type def above
    real(rp), intent(out) :: TScat_Path(:) ! TScat on the path defined by the
                                       ! interpolation factors in Eta_T, at
                                       ! frequency Frq.
    real(rp), intent(out) :: K_Temp_TScat(:)  ! dTScat / dT, T-SV length
    real(rp), intent(out) :: K_Atmos_TScat(:) ! dTScat / dVMR, vmr-SV length

    real(r8) :: F(2 )  ! Interpolation factors
    integer :: FrqIdx  ! Index of frq in TScat_Conf%signals(1)%frequencies
    integer :: I
    integer :: N       ! path length
    logical :: Need(2) ! Need to evaluate at new frequencies
    type(TScat_Linear_t) :: Wtemp ! for swapping

    n = size(TScat_Path)

    ! Decide for what frequencies we want data in W.
    call hunt ( TScat_Conf%signals(1)%frequencies, frq, frqIdx, &
              & allowTopValue=.true., allowBelowValue=.true. )

    ! Of the frequencies for which we want data in W, determine which if
    ! any we already have and which we need to compute.
    if ( frqIdx == 0 ) then
      if ( w(1)%got ) &
        & w(1)%got = w(1)%frq == TScat_Conf%signals(1)%frequencies(1)
      need(1) = .not. w(1)%got
      need(2) = .false.
    else if ( frqIdx < size(TScat_Conf%signals(1)%frequencies) ) then
      if ( w(1)%got ) then
        w(1)%got = w(1)%frq == TScat_Conf%signals(1)%frequencies(frqIdx)
        if ( .not. w(1)%got .and. w(2)%got ) then
          ! Use w(2) if it now has the frequency we should use for w(1)
          if ( w(2)%frq == TScat_Conf%signals(1)%frequencies(frqIdx) ) then
            ! Swap w(1) and w(2)
            wTemp = w(1)
            w(1) = w(2)
            w(2) = wTemp
          end if
        end if
      end if
      if ( w(2)%got ) &
        & w(2)%got = w(1)%frq == TScat_Conf%signals(1)%frequencies(frqIdx+1)
      need = .not. w%got
    else
      if ( w(2)%got ) &
        & w(2)%got = w(1)%frq == TScat_Conf%signals(1)%frequencies(frqIdx)
      need(1) = .false.
      need(2) = .not. w(2)%got
    end if

    do i = 1, 2
      if ( need(i) ) then
        ! Evaluate TScat and derivatives the frequency in w(i)
      end if
    end do

    ! If we have TScat at two frequencies, interpolate TScat radiance and
    ! derivatives to FRQ, else use the result from the one we have.

    if ( w(1)%got .and. w(2)%got ) then
      ! Interpolate the results from the linear model results in W(1:2)
      f(1) = ( w(2)%frq - frq ) / ( w(2)%frq - w(1)%frq )
      f(2) = ( frq - w(1)%frq ) / ( w(2)%frq - w(1)%frq )
      TScat_path = f(1)*w(1)%TScat_path(:n) + f(2)*w(2)%TScat_path(:n)
      if ( temp_der ) k_temp_TScat = f(1)*w(1)%k_temp_TScat(:n) + &
                                     f(2)*w(2)%k_temp_TScat(:n)
      if ( atmos_der ) k_atmos_TScat = f(1)*w(1)%k_atmos_TScat(:n) + &
                                          f(2)*w(2)%k_atmos_TScat(:n)
    else if ( w(1)%got ) then ! use linear model results in W(1)
      TScat_path = w(1)%TScat_path(:n)
      if ( temp_der ) k_temp_TScat = w(1)%k_temp_TScat(:n)
      if ( atmos_der ) k_atmos_TScat = w(1)%k_atmos_TScat(:n)
    else if ( w(2)%got ) then ! use linear model results in W(2)
      TScat_path = w(2)%TScat_path(:n)
      if ( temp_der ) k_temp_TScat = w(2)%k_temp_TScat(:n)
      if ( atmos_der ) k_atmos_TScat = w(2)%k_atmos_TScat(:n)
    else
      ! Can we get here?
      call MLSMessage ( MLSMSG_Error, moduleName, "No linear-model TScat results" )
    end if

  end subroutine Get_TScat

  ! --------------------------------------------  Get_TScat_Setup  -----
  subroutine Get_TScat_Setup ( TScat_Conf, FmConf, FwdModelIn, FwdModelExtra, &
    & Radiance, Sideband, MAF, &
    & L2PC, xP, DeltaX, XstarPTAN )

    ! Create an ersatz ForwardModelConfig to drive a linearized
    ! computation of TScat radiances during clear-sky path
    ! integration.

    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: L_PTAN, L_Radiance
    use L2PC_m, only: L2PC_t, L2PCDatabase, PopulateL2PCbin
    use L2PCBins_m, only: SelectL2PCBins
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: GetVectorQuantityByType, VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(out) :: TScat_Conf
    type(forwardModelConfig_T), intent(in) :: FMconf
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type (VectorValue_T), intent(in) :: RADIANCE ! The radiance we're after
    integer, intent(in) :: SIDEBAND       ! Which sideband (see below)
    integer, intent(in) :: MAF            ! MAF index
    type(L2PC_t), pointer :: L2PC         ! The selected L2PC
    type(vector_T), intent(out) :: XP     ! Same form as xStar, contents as x
    type(vector_T), intent(out) :: DELTAX ! xp-xStar
    type(vectorValue_T), pointer :: XSTARPTAN ! Tangent pressure in l2pc

    integer, dimension(-1:1) :: L2PCBINS ! Result
    type(vectorValue_T), pointer :: RADINL2PC ! The relevant radiance part of yStar
    integer :: Stat

    TScat_conf = fmConf
    TScat_conf%lockBins = .false.
    TScat_conf%molecules => fmConf%TScatMolecules
    TScat_conf%moleculeDerivatives => fmConf%TScatMoleculeDerivatives
    TScat_conf%phiWindow = 0.0
    TScat_conf%xStar = 0
    TScat_conf%yStar = 0

    allocate ( TScat_conf%signals(1), stat=stat )
    call test_allocate ( stat, moduleName, "TScat_Conf%signals", (/1/), (/1/), -1 )

    ! Use the same signals as FmConf%signals(1), but with all channels turned on
    TScat_conf%signals(1) = FmConf%signals(1)
    nullify ( TScat_conf%signals(1)%channels )
    call allocate_test ( TScat_conf%signals(1)%channels, &
      & size(FmConf%signals(1)%channels), "TScat_conf%signals(1)%channels", &
      & moduleName )
    TScat_conf%signals(1)%channels = .true.

    ! Determine the best bin to use
    call selectL2PCBins ( TScat_conf, fwdModelIn, fwdModelExtra, &
      & radiance, sideband, maf, &
      & l2pcBins )

    ! Load the bin
    call populateL2PCBin ( l2pcBins(sideband), ignoreHessian=.true. )

    ! Get the L2PC from the bin
    l2pc => l2pcDatabase(l2pcBins(sideband))

    radInl2pc => GetVectorQuantityByType ( &
      & l2pc%j%row%vec, quantityType=l_radiance, &
      & signal=TScat_conf%signals(1)%index, sideband=sideband )
    if ( radInL2PC%template%noChans /= radiance%template%noChans ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Channel dimension in l2pc not same as in measurements" )

    xStarPtan => GetVectorQuantityByType ( l2pc%j%col%vec, &
      & quantityType = l_ptan )

  end subroutine Get_TScat_Setup

  ! -----------------------------------------  Get_TScat_Teardown  -----
  subroutine Get_TScat_Teardown ( TScat_Conf )

    ! Clean up the ersatz ForwardModelConfig created by Get_TScat_Setup

    use Allocate_Deallocate, only: Deallocate_Test, Test_DeAllocate
    use ForwardModelConfig, only: ForwardModelConfig_t
    type(forwardModelConfig_T), intent(out) :: TScat_Conf
    integer :: Stat

    call deallocate_test ( TScat_Conf%signals(1)%channels, &
      & "TScat_Conf%signals%channels", moduleName )
    deallocate ( TScat_Conf%signals, stat=stat )
    call test_deallocate ( stat, moduleName, "TScat_Conf%signals", -1 )

  end subroutine Get_TScat_Teardown

  ! ---------------------------------  Interpolate_P_to_theta_e  -----
  subroutine Interpolate_P_to_theta_e ( P, Eta_t_iwc, Theta_e, Beg_Pos_Theta, &
    &                                   Xis, Coeffs, P_on_Xi )

    ! Interpolate P or dP_dT or dP_dIWC from Theta to Xi, for TScat
    ! generation, using periodic spline interpolation.

    use MLSKinds, only: Rp
    use MLSNumerics, only: Coefficients_r8, InterpolateValues

    real(rp), intent(in) :: P(0:,0:,:) ! P or dP_dT or dP_dIWC on T x IWC x Theta
    real(rp), intent(in) :: Eta_T_IWC(0:,0:) ! 2x2, to interpolate to T, IWC
    real(rp), intent(in) :: Theta_e(:) ! Theta extended to negative values
    integer, intent(in) :: Beg_Pos_Theta ! 1 or 2, depending on whether
                                       ! theta_e(1) is nonzero or zero
    real(rp), intent(in) :: Xis(:) ! Angles on which radiative transfer calculated
    type(coefficients_r8), intent(in) :: Coeffs ! To interpolate from Theta_e to Xis
    real(rp), intent(out) :: P_on_Xi(:) ! Interpolated result

    real(rp) :: P_T_IWC(size(p,3)) ! P interpolated to T and IWC
    real(rp) :: P_on_Theta_e(2*size(p_t_iwc)+1-beg_pos_theta) ! P * sin(abs(theta))
                             ! interpolated to scattering point IWC and T for
                             ! each theta (including negative theta), intermediary
                             ! to getting P_On_Xi

    integer :: I

    ! Interpolate P to scattering point IWC and temperature and
    ! multiply that by sin(theta)
    do i = 1, size(p_t_iwc)
      p_t_iwc(i) = &
        & sum(eta_t_iwc(0:1,0:1) * p(0:1,0:1,i) ) * abs(sin(theta_e(i)))
    end do

    ! The phase function is symmetric in theta, but the radiances
    ! are not.  Unfold p_t_iwc to negative theta.
    p_on_theta_e(1:size(p_t_iwc)) = p_t_iwc(size(p_t_iwc):1:-1)
    p_on_theta_e(size(p_t_iwc)+1:2*size(p_t_iwc)+1-beg_pos_theta) = &
      & p_t_iwc(beg_pos_theta:)

    ! Now interpolate P_On_Theta_e to P_On_Xi using a periodic spline
    call interpolateValues ( coeffs, theta_e, p_on_theta_e, &
      & xis, p_on_xi, method='S', extrapolate='P' )

  end subroutine Interpolate_P_to_theta_e

  ! --------------------------------------------  TScat_Gen_Setup  -----
  subroutine TScat_Gen_Setup ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    &                      FwdModelOut, NoUsedChannels, Temp, Sideband, &
    &                      IWC, ScatteringAngles )

    ! Get ready for computation of TScat table

    use ForwardModelConfig, only: Channels_T, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_ScatteringAngle, L_TSCAT, L_VMR
    use Molecules, only: L_CLOUD_A
    use Read_Mie_m, only: dP_dT, F_s, P
    use VectorsModule, only: Vector_T, VectorValue_T

    type (forwardModelConfig_t), intent(in) :: FwdModelConf
    type (vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type (vector_T), intent(in) :: FwdModelOut  ! Radiances, etc.
    type (vectorValue_T), intent(in) :: TEMP ! Temperature component of state vector
    integer, intent(in) :: NoUsedChannels, Sideband
    type (vectorValue_T), pointer :: IWC  ! IWC at scattering points
    type (vectorValue_T), pointer :: ScatteringAngles ! for TScat computation

    logical :: Error
    integer :: Q ! Loop index
    type (Channels_T), pointer, dimension(:) :: Channels
    type (VectorValue_T), pointer :: TScat, TScat2
    if ( .not. associated(p) .or. .not. associated(F_s) .or. &
      &  .not. associated(dP_dT) ) &
      & call announce_error ( 'TScat table computation requires Mie tables.' )
    if ( FwdModelConf%do_conv ) call announce_error ( &
      & 'Convolution and TScat computation are incompatible.' )
    if ( FwdModelConf%incl_cld ) call announce_error ( &
      & 'Cloud modeling and TScat computation are incompatible.' )
    if ( FwdModelConf%polarized ) call announce_error ( &
      & 'Cannot compute TScat for the polarized model.' )
    if ( FwdModelConf%refract ) call announce_error ( &
      & 'Refractive correction and TScat computation are incompatible.' )
    if ( FwdModelConf%spect_der ) call announce_error ( &
      & 'Spectroscopy derivatives and TScat computation are incompatible.' )
    channels => fwdModelConf%channels
    ! Get IWC field
    iwc => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra,  &
      & quantityType=l_vmr, molecule=l_cloud_a, config=fwdModelConf )
    ! Make sure all the TScat quantities for the selected signals have the
    ! same grids
    TScat => GetQuantityForForwardModel ( &
        & fwdModelOut, quantityType=l_TScat, &
        & signal=fwdModelConf%signals(channels(1)%signal)%index, &
        & sideband=sideband )
    if ( .not. TScat%template%coherent .or. &
      &  .not. TScat%template%stacked ) call announce_error ( &
        & 'TScat coordinates must be stacked and coherent' )
    TScat%values = 0.0
    error = .false.
    do q = 2, noUsedChannels
      TScat2 => GetQuantityForForwardModel ( &
        & fwdModelOut, quantityType=l_TScat, &
        & signal=fwdModelConf%signals(channels(q)%signal)%index, &
        & sideband=sideband )
      error = .not. TScat2%template%coherent .or. &
        &     .not. TScat2%template%stacked .or. &
        &     TScat%template%noChans /= TScat2%template%noChans .or. &
        &     TScat%template%hGridIndex /= TScat2%template%hGridIndex .or. &
        &     TScat%template%vGridIndex /= TScat2%template%vGridIndex
      if ( error ) exit
      TScat2%values = 0.0
    end do
    if ( error ) call announce_error ( "TScat quantities must all have the same grids." )
    ! Make sure IWC and temperature are coherent and stacked
    if ( .not. temp%template%stacked .or. .not. temp%template%coherent .or. &
       & .not.  iwc%template%stacked .or. .not.  iwc%template%coherent ) &
       & call announce_error ( "IWC and temperature must be coherent and stacked." )
    ! Make sure IWC and temperature have the same grids as TScat
    if ( TScat%template%hGridIndex /= temp%template%hGridIndex .or. &
       & TScat%template%vGridIndex /= temp%template%vGridIndex .or. &
       & TScat%template%hGridIndex /=  iwc%template%hGridIndex .or. &
       & TScat%template%vGridIndex /=  iwc%template%vGridIndex ) &
       & call announce_error ( "TScat, IWC and temperature must have the same grids." )
    scatteringAngles => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_scatteringAngle, config=fwdModelConf )

  contains

  ! .............................................  Announce_Error  .....
    subroutine Announce_Error ( Message )
    ! Announce Message using MLSMessage.  Include the configuration name
    ! in the message.  This is the same as in FullForwardModel.
      use MLSMessageModule, only: MLSMessage, MLSMSG_Error
      use String_Table, only: Get_String, String_Length
      character(len=*), intent(in) :: Message
      integer, parameter :: C = len('With config(')
      integer :: L
      character(511) :: Work ! Should be plenty of room
      l = string_length(fwdModelConf%name)
      work(:c) = 'With config('
      call get_string ( fwdModelConf%name, work(c+1:) )
      l = c + l + 1
      work(l:l+3) = '): '
      l = l + 3
      work(l+1:l+len_trim(message)) = message
      call MLSMessage ( MLSMSG_Error, moduleName, work(:l+len_trim(message)) )
    end subroutine Announce_Error

  end subroutine TScat_Gen_Setup

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module TScat_Support_m

! $Log$
! Revision 2.1  2010/08/27 23:39:13  vsnyder
! Initial commit -- far from working
!
