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
  public :: Get_TScat, Interpolate_P_to_theta_e
  public :: Get_TScat_Setup, TScat_Gen_Setup

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  ! --------------------------------------------------  Get_TScat  -----
  subroutine Get_TScat ( TScat_Conf, State, MAF, Phitan, Frq, Z_path,        &
    &                    Phi_path, Tan_PT, Grids_tmp, Eta_T, Grids_f, Eta_F, &
    &                    L2PC, X, dX, TScat, RadInL2PC, Grids_TScat,         &
    &                    Grids_L2PC_X, XstarPTAN, TScat_Path, dTScat )

  ! Get TScat and its derivatives along the integration path using
  ! the TScat tables and the same strategy as the linear forward model.

    use Comp_Eta_Docalc_No_Frq_m, only: Comp_Eta_Docalc_No_Frq, Comp_Eta_fzp
    use Comp_Sps_Path_Frq_m, only: Comp_One_Path_Frq
    use ForwardModelConfig, only: ForwardModelConfig_t
    use HGridsDatabase, only: FindClosestMatch
    use L2PC_m, only: L2PC_t
    use Load_Sps_Data_m, only: Grids_T, Fill_Grids_2
    use ManipulateVectorQuantities, only: FindOneClosestInstance
    use MatrixModule_1, only: MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP, &
      & FINDBLOCK, CREATEBLOCK
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use Molecules, only: L_CLOUD_A
    use String_Table, only: Get_String
    use VectorsModule, only: GetVectorQuantityByType, &
      & VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(in) :: TScat_Conf
    type(vector_t), intent(in) :: State ! Forward model's state vector.
    integer, intent(in) :: MAF          ! Major frame index
    type(vectorValue_T), intent(in) :: Phitan ! Limb tangent geodAngle
    real(r8), intent(in) :: Frq         ! Frequency at which to evaluate TScat
    real(rp), intent(in) :: Z_Path(:)   ! Zeta along the path
    real(rp), intent(in) :: Phi_Path(:) ! Phi along the path
    integer, intent(in) :: Tan_PT       ! Index of tangent point in Z_Path
    type (Grids_T), intent(in) :: Grids_tmp ! All the coordinates for TEMP
    real(rp), intent(in) :: Eta_T(:,:)  ! Interpolation factors from temperature grid
                                        ! to path, path length X SV length
    type (Grids_T), intent(in) :: Grids_f   ! All the coordinates for VMR
    real(rp), intent(in) :: Eta_F(:,:)  ! Interpolation factors from VMR grid
                                        ! to path, path length X SV length
    type(L2PC_t), intent(in) :: L2PC    ! The selected L2PC
    type(vector_T), intent(inout) :: X  ! Quantities from State that correspond
                                        ! to quantities from L2PC%col%vec. Work
                                        ! space to avoid cloning; (inout) to
                                        ! avoid nullifying pointer components.
    type(vector_T), intent(inout) :: dX ! X - L2PC%col%vec. Work space to avoid
                                        ! cloning; (inout) to avoid nullifying
                                        ! pointer components.
    type(vector_T), intent(inout) :: TScat ! Same form as L2PC%row%vec. Work
                                        ! space to avoid cloning; (inout) to
                                        ! avoid nullifying pointer components.
    integer, intent(in) :: RadInL2PC    ! Index of qty in TScat, X, dX, L2PC%J
    type(grids_t), intent(inout) :: Grids_TScat ! All the grids for TScat
    type(grids_t), intent(inout) :: Grids_L2PC_X ! All the grids for L2PC%j%col%vec
    type(vectorValue_T), intent(inout) :: XSTARPTAN ! Tangent pressure in l2pc
    real(rp), intent(out) :: TScat_Path(:) ! TScat on the path.
    real(rp), intent(out) :: dTScat(:,:) ! dTScat / etc.
                                        ! path X size(L2PC%col%vec%quantities)

    integer :: Center      ! of TQ instances
    integer :: ClosestInstance ! Instance of SQ nearest to instance of TScat
    integer :: DeltaInstance   ! Distance from ClosestInstance
    real(r8) :: DeltaPhi
    logical :: Do_Calc_FZP(size(z_path), max(size(grids_TScat%values), size(grids_L2PC_x%values)))
    logical :: Do_Calc_ZP(size(z_path), max(grids_TScat%p_len,grids_L2PC_x%p_len))
    real(rp) :: Eta_FZP(size(z_path), max(size(grids_TScat%values), size(grids_L2PC_x%values)))
    real(rp) :: Eta_ZP(size(z_path), max(grids_TScat%p_len,grids_L2PC_x%p_len))
    integer :: F_Inda, F_Indb
    integer :: I           ! Loop index and subscript
    integer :: Instance    ! For which instance in State do we want TScat?
    integer :: N           ! path length
    logical :: Not_Zero_f(grids_L2PC_x%l_f(ubound(grids_L2PC_x%l_f,1)))
    integer :: Qty         ! Loop index and subscript
    type(vectorValue_t), pointer :: RQ ! TScat radiance quantity
    integer :: RQProfile
    type(vectorValue_t), pointer :: SQ ! State vector quantity
    integer :: SQInstance  ! Instance of SQ
    integer :: SV_f, SV_zp ! Loop inductor and subscript
    type(vectorValue_t), pointer :: TQ ! L2PC TScat state vector quantity
    integer :: TQInstance  ! Instance of TQ
    integer :: TSProfile   ! L2PC TScat profile
    integer :: V_Ind
    integer :: W_Inda, W_Indb
    character(len=80) :: Word ! From the string table

    n = size(TScat_Path)

    rq => TScat%quantities(radInL2PC)
    do rqProfile = 1, size(rq%template%phi)   ! Profile in TScat rad qty
      do qty = 1, size(l2pc%j%col%vec%quantities,1) ! TScat state quantities
        ! Get quantity in State corresponding to quantity in L2PC TScat's
        ! column vector
        tq => l2pc%j%col%vec%quantities(qty)
        sq => getVectorQuantityByType( State, quantityType=tq%template%quantityType )
        if ( .not. associated(sq) ) then ! Complain
          call get_string ( tq%template%name, word, strip=.true. )
          call MLSMessage ( MLSMSG_Error, ModuleName, &
            &  "No quantity in state vector to match L2PC TScat quantity "//trim(word) )
        end if
        ! Compute dX using instances of sq most nearly corresponding to
        ! instances of tq.
        closestInstance = findClosestMatch( sq%template%phi(1,:), &
          & rq%template%phi, rqProfile )
        center = tq%template%noInstances/2 + 1
        do tqInstance = 1, tq%template%noInstances
          deltaInstance = tqInstance - center
          phiWindowLoop: do
            deltaPhi = tq%template%phi(1,center+deltaInstance) - &
              & tq%template%phi(1,center)
            if ( deltaPhi > TScat_conf%phiWindow ) then
              deltaInstance = deltaInstance - 1
            else if ( deltaPhi < -TScat_conf%phiWindow ) then
              deltaInstance = deltaInstance + 1
            else
              exit phiWindowLoop
            end if
          end do phiWindowLoop
          sqInstance = max(1, min( closestInstance + deltaInstance, &
                                   sq%template%noInstances ) )
          ! Fill the tqInstance part of dX
          dX%quantities(qty)%values(:,tqInstance) = sq%values(:,sqInstance) - &
            &                                       tq%values(:,tqInstance)
          ! Interpolate Jacobian to phitan here?
        end do ! tqInstance
      end do ! qty
    end do ! rqProfile

    ! Evaluate TScat%quantities(radInL2PC) using first-order approximation.
    ! TScat(out) = TScat(L2PC) + J * dX
    TScat%quantities(radInL2PC)%values = l2pc%j%row%vec%quantities(radInL2PC)%values
    call MultiplyMatrixVectorNoT ( l2pc%J, dX, TScat, update = .true., row=radInL2PC )

    ! Put values into Grids_TScat; grids were put there by Get_TScat_Setup.
    call Fill_Grids_2 ( grids_TScat, i, TScat%quantities(radInL2PC), &
      & setDerivFlags = .true. )

    ! Get eta_zp
    call comp_eta_docalc_no_frq ( grids_TScat, z_path, phi_path, eta_zp, &
      & do_calc_zp, tan_pt = tan_pt )
    ! Get TScat_path.  sideband=0 means absolute frequency, not I.F.
    call comp_one_path_frq ( grids_TScat, frq, eta_zp, do_calc_zp, &
      & TScat_path, do_calc_fzp, eta_fzp, sideband=0 )

    ! We have TScat on the path.  Now get the derivatives.

    call comp_eta_fzp ( grids_L2PC_x, frq, eta_zp, do_calc_zp, 0, &
                      & eta_fzp, not_zero_f, do_calc_fzp )

    dTScat = 0.0
    f_inda = 0
    w_inda = 0

    do qty = 1, size(l2pc%j%col%vec%quantities)
      f_indb = grids_L2PC_x%l_f(qty)
      v_ind = grids_L2PC_x%l_v(qty)
      w_indb = w_inda + (Grids_L2PC_x%l_z(qty) - Grids_L2PC_x%l_z(qty-1)) * &
                        (Grids_L2PC_x%l_p(qty) - Grids_L2PC_x%l_p(qty-1))

      do sv_zp = w_inda + 1, w_indb
        do sv_f = f_inda + 1, f_indb
          v_ind = v_ind + 1
          if ( not_zero_f(sv_f) ) then
            dTscat(:,qty) = dTscat(:,qty) + &
                          & Grids_L2PC_x%values(v_ind) * eta_fzp(:,v_ind)
          else
            dTscat(:,qty) = 0.0
          end if
        end do
      end do
    end do

  end subroutine Get_TScat

  ! --------------------------------------------  Get_TScat_Setup  -----
  subroutine Get_TScat_Setup ( TScat_Conf, FmConf, FwdModelIn, FwdModelExtra, &
    &                          Radiance, Sideband, MAF, Phitan, &
    &                          L2PC, X, dX, TScat, RadInL2PC, Grids_TScat, &
    &                          Grids_L2PC_X, XstarPTAN )

    ! Create an ersatz ForwardModelConfig to drive a linearized
    ! computation of TScat radiances during clear-sky path
    ! integration.

    use Allocate_Deallocate, only: Allocate_Test, Test_Allocate
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Intrinsic, only: L_PTAN, L_Radiance
    use Load_Sps_Data_m, only: Grids_T, Load_One_Item_Grid, Load_Grid_From_Vector
    use L2PC_m, only: L2PC_t, L2PCDatabase, PopulateL2PCbin
    use L2PCBins_m, only: SelectL2PCBins
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use VectorsModule, only: CloneVector, GetVectorQuantityByType, &
      & GetVectorQuantityIndexByType, VectorValue_T, Vector_T

    type(forwardModelConfig_T), intent(out) :: TScat_Conf
    type(forwardModelConfig_T), intent(in) :: FMconf
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type (VectorValue_T), intent(in) :: RADIANCE ! The radiance we're after
    integer, intent(in) :: SIDEBAND       ! Which sideband (see below)
    integer, intent(in) :: MAF            ! MAF index
    type(vectorValue_T), intent(in) :: Phitan ! Limb tangent geodAngle
    type(L2PC_t), pointer :: L2PC         ! The selected L2PC
    type(vector_T), intent(out) :: X      ! Same form as L2PC%col%vec
    type(vector_T), intent(out) :: dX     ! Same form as L2PC%col%vec
    type(vector_T), intent(out) :: TScat  ! Same form as L2PC%row%vec
    integer, intent(out) :: RadInL2PC     ! Qty in L2PC for TScat_conf%signals(1)
    type(grids_t), intent(inout) :: Grids_TScat ! All the grids for TScat
    type(grids_t), intent(inout) :: Grids_L2PC_X ! All the grids for L2PC%j%col%vec
    type(vectorValue_T), pointer :: XSTARPTAN ! Tangent pressure in L2PC

    integer, dimension(-1:1) :: L2PCBINS ! Result
    integer :: Stat

    TScat_conf = fmConf
    TScat_conf%lockBins = .false.
    TScat_conf%molecules => fmConf%TScatMolecules
    TScat_conf%moleculeDerivatives => fmConf%TScatMoleculeDerivatives
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

    radInl2pc = GetVectorQuantityIndexByType ( &
      & l2pc%j%row%vec, quantityType=l_radiance, &
      & signal=TScat_conf%signals(1)%index, sideband=sideband )
    if ( l2pc%j%row%vec%quantities(radInL2PC)%template%noChans /= radiance%template%noChans ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Channel dimension in l2pc not same as in measurements" )

    xStarPtan => GetVectorQuantityByType ( l2pc%j%col%vec, &
      & quantityType = l_ptan )

    call cloneVector ( x, l2pc%j%col%vec, vectorNameText='__x' )
    call cloneVector ( dX, x, vectorNameText='__dX' )
    call cloneVector ( TScat, l2pc%j%row%vec, vectorNameText='__TScat' )

    ! Create and fill grids for TScat radiance
    call load_one_item_grid ( grids_TScat, TScat%quantities(radInl2pc), &
      & phitan, MAF, TScat_Conf, setDerivFlags = .true. )

    ! Create and fill grids for L2PC state vector.
    call load_grid_from_vector ( grids_L2PC_x, L2PC%j%col%vec, phitan, MAF, &
      & TScat_Conf, setDerivFlags=spread(.true.,1,size(L2PC%j%col%vec%quantities)) )
    ! Save some space
    call allocate_test ( grids_L2PC_x%values, 0, 'Grids_L2PC_x%values', moduleName )

  end subroutine Get_TScat_Setup

  ! -----------------------------------------  Get_TScat_Teardown  -----
  subroutine Get_TScat_Teardown ( TScat_Conf, X, dX, TScat, Grids_TScat, &
                                & Grids_L2PC_X )

    ! Clean up the ersatz ForwardModelConfig created by Get_TScat_Setup

    use Allocate_Deallocate, only: Deallocate_Test, Test_DeAllocate
    use ForwardModelConfig, only: ForwardModelConfig_t
    use Load_Sps_Data_m, only: DestroyGrids_t, Grids_T
    use VectorsModule, only: DestroyVectorInfo, Vector_t

    type(forwardModelConfig_T), intent(out) :: TScat_Conf
    type(vector_t), intent(inout) :: X      ! Same form as L2PC%col%vec
    type(vector_t), intent(inout) :: dX     ! Same form as L2PC%col%vec
    type(vector_t), intent(inout) :: TScat  ! Same form as L2PC%row%vec.
    type(grids_t), intent(inout) :: Grids_TScat, Grids_L2PC_X

    integer :: Stat

    call deallocate_test ( TScat_Conf%signals(1)%channels, &
      & "TScat_Conf%signals%channels", moduleName )
    deallocate ( TScat_Conf%signals, stat=stat )
    call test_deallocate ( stat, moduleName, "TScat_Conf%signals", -1 )

    call destroyVectorInfo ( x )
    call destroyVectorInfo ( dX )
    call destroyVectorInfo ( TScat ) 

    call destroyGrids_t ( grids_TScat )
    call destroyGrids_t ( grids_L2PC_X )

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
  subroutine TScat_Gen_Setup ( FwdModelConf, FwdModelIn, FwdModelExtra,     &
    &                          FwdModelOut, NoUsedChannels, Temp, Sideband, &
    &                          IWC, ScatteringAngles )

    ! Get ready for computation of TScat table

    use ForwardModelConfig, only: Channels_T, ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_ScatteringAngle, L_TSCAT, L_VMR
    use Molecules, only: L_CLOUD_A
    use Read_Mie_m, only: dP_dT, F_s, P
    use VectorsModule, only: Vector_T, VectorValue_T

    ! Inputs:
    type (forwardModelConfig_t), intent(in) :: FwdModelConf
    type (vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type (vector_T), intent(in) :: FwdModelOut  ! Radiances, etc.
    type (vectorValue_T), intent(in) :: TEMP ! Temperature component of state vector
    integer, intent(in) :: NoUsedChannels, Sideband
    ! Outputs:
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
! Revision 2.2  2010/09/25 01:09:01  vsnyder
! Inching along
!
! Revision 2.1  2010/08/27 23:39:13  vsnyder
! Initial commit -- far from working
!
