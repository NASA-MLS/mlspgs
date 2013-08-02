! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Convolve_All_m

  ! Subprograms CONVOLVE_... add the effects of antenna smearing.
  ! Call FOV_Convolve_Setup before them, and FOV_Convolve_Teardown after them.

  ! Subprograms INTERPOLATE_... simply interpolate from the calculation grid
  ! to the output grid.  Call InterpolateArraySetup with METHOD='S',
  ! EXTRAPOLATE='C', dyByDx=ptan_der before them, and InterpolateArrayTeardown
  ! after them.

  implicit NONE
  private
  public :: CONVOLVE_OTHER_DERIV, CONVOLVE_RADIANCE, CONVOLVE_TEMPERATURE_DERIV
  public :: CONVOLVE_RADIANCE_NORMALIZATION, CONVOLVE_TEMPERATURE_DERIV_NORMALIZATION
  public :: CONVOLVE_OTHER_SECOND_DERIV
  public :: INTERPOLATE_RADIANCE, INTERPOLATE_OTHER_DERIV
  public :: INTERPOLATE_TEMPERATURE_DERIV
  public :: STORE_OTHER_DERIV, STORE_TEMPERATURE_DERIV

  interface Store_Other_Deriv
    module procedure Store_Other_Deriv_1D, Store_Other_Deriv_2D
  end interface

  interface Store_Temperature_Deriv
    module procedure Store_Temperature_Deriv_1D, Store_Temperature_Deriv_2D
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  ! ---------------------------------  Convolve_Radiance  -----
  subroutine Convolve_Radiance ( Convolve_Support, MAF, Channel, Rad_In, &
           & SbRatio, Update, Ptan, Radiance, MIF_Times, DeadTime, &
           & Jacobian, RowFlags, dh_dz_out, dx_dh_out, ptan_Der, Rad_FFT )

  ! Convolve the radiance, and maybe dRadiance/dPtan, with the antenna pattern

    use FOV_Convolve_m, only: Convolve_Support_t, FOV_Convolve_1d
    use MatrixModule_1, only: FINDBLOCK, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use VectorsModule, only: VectorValue_T

    ! Required inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Rad_In(:)   ! input radiances
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update       ! "add to radiance, don't overwrite"
    type(vectorvalue_t), intent(in) :: PTAN ! Only for some indices
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.

    ! Output
    type (VectorValue_T), intent(inout) :: RADIANCE   ! Output radiances

    ! Optional if PTan derivative processed.  All or none.
    type (Matrix_t), optional, intent(inout) :: Jacobian
    logical, optional, intent(inout) :: rowFlags(:) ! Flag to calling code
    real(rp), optional, intent(in) :: dh_dz_out(:) ! dh/dz on the output pointing grid
    real(rp), optional, intent(in) :: dx_dh_out(:) ! dx/dh on the output pointing grid
    logical, optional, intent(in) :: Ptan_der ! "Process PTAN derivatives"

    ! Temperature derivatives need this
    real(r8), intent(out), optional :: Rad_FFT(:) ! Convolved radiances on FFT grid

    logical :: my_ptan_der
    integer :: Col, NoChans, NoPtan, Row
    real(rv), pointer :: MIF_Times_for_MAF(:)
    real(r8) :: SRad(ptan%template%noSurfs), di_dx(ptan%template%noSurfs)

    my_ptan_der = .false.
    if ( present ( ptan_der ) ) my_ptan_der = ptan_der

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans
    noPtan = ptan%template%nosurfs

    ! Convolve the radiances (and maybe the derivative w.r.t. PTan):

    if ( my_ptan_der ) then ! Convolve radiance and get di_dx
      call fov_convolve_1d ( convolve_support, rad_in, MIF_Times_for_MAF, DeadTime, &
        & SRad, dRad_dx_out=di_dx, rad_fft_out=rad_fft )
    else
      call fov_convolve_1d ( convolve_support, rad_in, MIF_Times_for_MAF, DeadTime, &
        & SRad, rad_fft_out=rad_fft )
    end if

    ! Load the Radiance values into the Radiance structure:

    call loadVectorValue ( sRad, Radiance%values(channel::noChans,maf), &
      & sbRatio, update )

    if ( my_ptan_der ) then ! Jacobian better be present!

      ! First, find index location in Jacobian and set the derivative flags

      row = FindBlock( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .TRUE.

      ! Compute dI/dPtan using the chain rule:
      SRad = di_dx * dx_dh_out * dh_dz_out

      col = FindBlock ( Jacobian%col, ptan%index, maf )
      call getBandedBlock ( jacobian, row, col, noChans, noPtan )
      call loadMatrixValue ( sRad, Jacobian%block(row,col)%values(channel::noChans,1), &
        & sbRatio, update )

    end if ! present(jacobian) .and. my_ptan_der

  end subroutine Convolve_Radiance

  ! ---------------------------------  Convolve_Temperature_Deriv  -----
  subroutine Convolve_Temperature_Deriv ( Convolve_Support, MAF, Channel, &
           & Rad_In, Rad_FFT, SbRatio, Update, Radiance, Temp, Grids_Tmp, &
           & surf_angle, MIF_Times, DeadTime, di_dT, dx_dT, d2x_dxdT, &
           & dxdt_tan, dxdt_surface, Jacobian, RowFlags )

    use Fov_Convolve_m, only: Convolve_Support_T, &
      & FOV_Convolve_Temp_Derivs
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Rad_In(:)  ! input radiances
    real(r8), intent(in) :: Rad_FFT(:) ! convolved radiance on FFT grid
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type (vectorvalue_t), intent(in) :: TEMP     ! Only for some indices
    type (Grids_T), intent(in) :: Grids_Tmp      ! Temperature's grids, etc.
    real(rp), intent(in) :: Surf_angle ! An angle that defines the
    !                      Earth surface.
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: dx_dT(:,:) ! derivative of angle wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: d2x_dxdT(:,:) ! 2nd derivative wrt angle and
    !                      temperature on chi_in
    real(rp), intent(in) :: dxdT_tan(:,:) ! derivative of angle wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: dxdT_surface(:,:) ! derivative of angle
    !                      wrt temperature at the surface

    ! Outputs
    type (Matrix_t), intent(inout) :: Jacobian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Local variables
    integer :: Col, JF, K, NoChans, NoPtan, N_T_Zeta, Row, SV_I
    real(r8) :: dRad_dT_out(size(convolve_support%del_chi_out), &
                & temp%template%noSurfs*( &
                  & grids_tmp%windowFinish(1) - grids_tmp%windowStart(1) + 1))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    ! Start here

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans
    noPtan = size(convolve_support%del_chi_out)
    n_t_zeta = temp%template%noSurfs

    call fov_convolve_temp_derivs ( convolve_support, rad_in, &
      & rad_fft, surf_angle, MIF_Times_for_MAF, DeadTime, dI_dT, dx_dT, &
      & d2x_dxdT, dxdt_tan - SPREAD(dxdt_surface(1,:),1,noPtan), &
      & grids_tmp%deriv_flags, dRad_dT_out )

    ! Load the Temp. derivative values into the Jacobian
    ! First, find index location in Jacobian and set the derivative flags

    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    sv_i = 0
    do jf = grids_tmp%windowStart(1), grids_tmp%windowFinish(1)

      col = FindBlock ( Jacobian%col, temp%index, jf )
      call getFullBlock ( jacobian, row, col, 'temperature' )

      do k = 1, n_t_zeta

        ! Load derivatives for this (zeta & phi) if needed :

        sv_i = sv_i + 1
        if ( grids_tmp%deriv_flags(sv_i) ) &
          & call loadMatrixValue ( dRad_dT_out(:,sv_i), &
            & jacobian%block(row,col)%values(channel::noChans,k), sbRatio, &
            & update )
      end do

    end do

  end subroutine Convolve_Temperature_Deriv

  ! ---------------------------------  Convolve_Radiance_Normalization  -----
  ! Convolves Radiance, similar to Convolve_Radiance subroutine.
  ! This implemenation also computes Rad_Diff and Rad_Diff_FFT that 
  ! are required for Convolve_Temperature_Deriv_Normalization subroutine.

  subroutine Convolve_Radiance_Normalization ( Convolve_Support, MAF, Channel, Rad_In, &
           & SbRatio, Update, Ptan, Radiance, MIF_Times, DeadTime, &
           & Jacobian, RowFlags, dh_dz_out, dx_dh_out, ptan_Der, Rad_FFT, Rad_Diff, Rad_Diff_FFT )

  ! Convolve the radiance, and maybe dRadiance/dPtan, with the antenna pattern

    use FOV_Convolve_m, only: Convolve_Support_t, FOV_Convolve_1d
    use MatrixModule_1, only: FINDBLOCK, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use MLSNumerics, only: InterpolateValues    ! IGOR
    use VectorsModule, only: VectorValue_T

    ! Required inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Rad_In(:)   ! input radiances, I
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update       ! "add to radiance, don't overwrite"
    type(vectorvalue_t), intent(in) :: PTAN ! Only for some indices
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.

    ! Output
    type (VectorValue_T), intent(inout) :: RADIANCE   ! Output radiances, IA = I*G

    ! Optional if PTan derivative processed.  All or none.
    type (Matrix_t), optional, intent(inout) :: Jacobian
    logical, optional, intent(inout) :: rowFlags(:) ! Flag to calling code
    real(rp), optional, intent(in) :: dh_dz_out(:) ! dh/dz on the output pointing grid
    real(rp), optional, intent(in) :: dx_dh_out(:) ! dx/dh on the output pointing grid
    logical, optional, intent(in) :: Ptan_der ! "Process PTAN derivatives"

    ! Temperature derivatives need this    IGOR - in new version, rad_fft will not be needed,
    !        but rad_diff (= I-IA) and rad_diff_FFT (= FFT(I-IA)) will be 
    real(r8),             intent(out)  , optional :: Rad_FFT(:)      ! FFT(I), on FFT grid
    real(rp),             intent(out)  , optional :: Rad_Diff(:)     ! Output radiance differences, I-IA    ! IGOR
    real(r8),             intent(out)  , optional :: Rad_Diff_FFT(:) ! FFT(Rad_Diff) = FFT(I-IA)            ! IGOR

    logical :: my_ptan_der
    integer :: Col, NoChans, NoPtan, Row
    real(rv), pointer :: MIF_Times_for_MAF(:)
    real(r8) :: SRad(ptan%template%noSurfs), di_dx(ptan%template%noSurfs)
    real(r8) :: SRad_out(size(rad_in))                             ! IGOR
    real(r8) :: SRad_temp(ptan%template%noSurfs)                   ! IGOR

    my_ptan_der = .false.
    if ( present ( ptan_der ) ) my_ptan_der = ptan_der

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans
    noPtan = ptan%template%nosurfs

    ! Convolve the radiances (and maybe the derivative w.r.t. PTan):

    if ( my_ptan_der ) then ! Convolve radiance and get di_dx
      call fov_convolve_1d ( convolve_support, rad_in, MIF_Times_for_MAF, DeadTime, &
        & SRad, dRad_dx_out=di_dx, rad_fft_out=rad_fft )
    else
      call fov_convolve_1d ( convolve_support, rad_in, MIF_Times_for_MAF, DeadTime, &
        & SRad, rad_fft_out=rad_fft )
    end if

    ! Load the Radiance values into the Radiance structure:

    call loadVectorValue ( sRad, Radiance%values(channel::noChans,maf), &
      & sbRatio, update )

    ! IGOR
    call InterpolateValues ( convolve_support%del_chi_out, sRad, convolve_support%del_chi_in, sRad_out, &
                                 & METHOD='S', extrapolate='C' )

    ! I_diff = I - IA            ! IGOR:
    
    rad_diff = rad_in - sRad_out      ! IGOR

    ! IGOR:
    if ( my_ptan_der ) then ! Get FFT(I-IA)
      call fov_convolve_1d ( convolve_support, rad_diff, MIF_Times_for_MAF, DeadTime, &
        & sRad_temp, rad_fft_out=Rad_Diff_FFT )
    end if

    ! Load the I-IA values into the Radiance_Diff structure:

    !call loadVectorValue ( rad_diff, Radiance_diff%values(channel::noChans,maf), &
    !  & sbRatio, update )          ! IGOR


    if ( my_ptan_der ) then ! Jacobian better be present!

      ! First, find index location in Jacobian and set the derivative flags

      row = FindBlock( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .TRUE.

      ! Compute dI/dPtan using the chain rule:
      SRad = di_dx * dx_dh_out * dh_dz_out

      col = FindBlock ( Jacobian%col, ptan%index, maf )
      call getBandedBlock ( jacobian, row, col, noChans, noPtan )
      call loadMatrixValue ( sRad, Jacobian%block(row,col)%values(channel::noChans,1), &
        & sbRatio, update )

    end if ! present(jacobian) .and. my_ptan_der

  end subroutine Convolve_Radiance_Normalization

  ! ---------------------------------  Convolve_Temperature_Deriv_Normalization  -----
  ! Performs convolution of Temperature derivatives (with normalization)
  subroutine Convolve_Temperature_Deriv_Normalization ( Convolve_Support, MAF, Channel, &
           & Rad_Diff, Rad_Diff_FFT, SbRatio, Update, Radiance, Temp, Grids_Tmp, &
           & surf_angle, MIF_Times, DeadTime, di_dT, dx_dT, d2x_dxdT, &
           & dxdt_tan, dxdt_surface, Jacobian, RowFlags )

    use Fov_Convolve_m, only: Convolve_Support_T, &
      & FOV_Convolve_Temp_Derivs, FOV_Convolve_Temp_Derivs_Normalization
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Rad_Diff(:)  ! Input radiance differences, I-IA        ! IGOR
    real(r8), intent(in) :: Rad_Diff_FFT(:) ! FFT(Rad_Diff) = FFT(I-IA)            ! IGOR
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type (vectorvalue_t), intent(in) :: TEMP     ! Only for some indices
    type (Grids_T), intent(in) :: Grids_Tmp      ! Temperature's grids, etc.
    real(rp), intent(in) :: Surf_angle ! An angle that defines the
    !                      Earth surface.
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: dx_dT(:,:) ! derivative of angle wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: d2x_dxdT(:,:) ! 2nd derivative wrt angle and
    !                      temperature on chi_in
    real(rp), intent(in) :: dxdT_tan(:,:) ! derivative of angle wrt
    !                      temperature on chi_in
    real(rp), intent(in) :: dxdT_surface(:,:) ! derivative of angle
    !                      wrt temperature at the surface

    ! Outputs
    type (Matrix_t), intent(inout) :: Jacobian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Local variables
    integer :: Col, JF, K, NoChans, NoPtan, N_T_Zeta, Row, SV_I
    real(r8) :: dRad_dT_out(size(convolve_support%del_chi_out), &
                & temp%template%noSurfs*( &
                  & grids_tmp%windowFinish(1) - grids_tmp%windowStart(1) + 1))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    ! Start here

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans
    noPtan = size(convolve_support%del_chi_out)
    n_t_zeta = temp%template%noSurfs

    !call fov_convolve_temp_derivs ( convolve_support, rad_in, &
    !  & rad_fft, surf_angle, MIF_Times_for_MAF, DeadTime, dI_dT, dx_dT, &
    !  & d2x_dxdT, dxdt_tan - SPREAD(dxdt_surface(1,:),1,noPtan), &
    !  & grids_tmp%deriv_flags, dRad_dT_out )                          ! IGOR

    call fov_convolve_temp_derivs_normalization ( convolve_support, rad_diff, &
      & rad_diff_fft, surf_angle, MIF_Times_for_MAF, DeadTime, dI_dT, dx_dT, &
      & d2x_dxdT, dxdt_tan - SPREAD(dxdt_surface(1,:),1,noPtan), &
      & grids_tmp%deriv_flags, dRad_dT_out )   ! IGOR


    ! Load the Temp. derivative values into the Jacobian
    ! First, find index location in Jacobian and set the derivative flags

    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    sv_i = 0
    do jf = grids_tmp%windowStart(1), grids_tmp%windowFinish(1)

      col = FindBlock ( Jacobian%col, temp%index, jf )
      call getFullBlock ( jacobian, row, col, 'temperature' )

      do k = 1, n_t_zeta

        ! Load derivatives for this (zeta & phi) if needed :

        sv_i = sv_i + 1
        if ( grids_tmp%deriv_flags(sv_i) ) &
          & call loadMatrixValue ( dRad_dT_out(:,sv_i), &
            & jacobian%block(row,col)%values(channel::noChans,k), sbRatio, &
            & update )
      end do

    end do

  end subroutine Convolve_Temperature_Deriv_Normalization

  ! ----------------------------------------  Convolve_OtherDeriv  -----
  subroutine Convolve_Other_Deriv ( Convolve_Support, MAF, Channel, &
             & SbRatio, Update, Radiance, Qtys, Grids_f, MIF_Times, &
             & DeadTime, dI_df, Jacobian, RowFlags, ExtraJacobian, &
             & Derivs )

    use ForwardModelConfig, only: QtyStuff_T
    use Fov_Convolve_m, only: Convolve_Support_T, &
      & FOV_Convolve_2d
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_df(:,:) ! mixing ratio derivatives or any
    !                                    parameter for which a simple
    !                                    convolution will suffice

    ! Outputs
    type (Matrix_t), intent(inout), target :: Jacobian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Optional output
    type (Matrix_t), intent(inout), optional, target :: ExtraJacobian

    ! Optional input
    logical, intent(in), optional :: Derivs(:) ! Same size as Qtys; if not
                                               ! present, do all columns.  Use
                                               ! in case grids_f%deriv_flags
                                               ! hasn't been set.

    ! Local variables
    integer :: Row, Col
    integer :: JF
    integer :: K
    integer :: NFZ    ! # frequency elements times # z elements in a specie
    integer :: NoChans
    integer :: SPS_I  ! species index
    integer :: SV_F   ! index of final element in a state vector
    real(r8) :: drad_df_out(size(convolve_support%del_chi_out), &
      &                     size(di_df,dim=2))
    real(rv), pointer :: MIF_Times_for_MAF(:)
    type (Matrix_t), pointer :: MyJacobian

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans

    ! do the convolution

    call fov_convolve_2d ( convolve_support, di_df, MIF_Times_for_MAF, DeadTime, &
      & grids_f%deriv_flags, drad_df_out )

    ! load derivatives into jacobian
    ! First, find index location in Jacobian and set the derivative flags

    row = FindBlock( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    do sps_i = 1, size(qtys)

      if ( present(derivs) ) then
        if ( .not. derivs(sps_i) ) cycle
      end if

      if ( qtys(sps_i)%foundInFirst ) then
        myJacobian => jacobian
      else
        if ( .not. present(extraJacobian) ) cycle
        myJacobian => extraJacobian
      end if

      sv_f = grids_f%l_v(sps_i-1)
      nfz = (Grids_f%l_f(sps_i) - Grids_f%l_f(sps_i-1)) * &
          & (Grids_f%l_z(sps_i) - Grids_f%l_z(sps_i-1))

      do jf = Grids_f%windowStart(sps_i), Grids_f%windowfinish(sps_i)

        col = FindBlock ( myJacobian%col, qtys(sps_i)%qty%index, jf)
        call getFullBlock ( myJacobian, row, col, 'atmospheric' )

        do k = 1, nfz

          ! load derivatives for this (zeta & phi) if needed:

          sv_f = sv_f + 1
          if ( Grids_f%deriv_flags(sv_f) ) &
            & call loadMatrixValue ( drad_df_out(:,sv_f), &
              & myJacobian%block(row,col)%values(channel::noChans,k), sbRatio, &
              & update )
        end do

      end do

    end do

  end subroutine Convolve_Other_Deriv

  ! -------------------------------------  Convolve_Other_Second_Deriv -----
  subroutine Convolve_Other_Second_Deriv ( Convolve_Support, MAF, Channel, &
             & SbRatio, Update, Radiance, Qtys, Grids_f, &
             & MIF_Times, DeadTime, d2I_df2, Hessian, RowFlags )
    use Dump_0, only: DUMP
    use ForwardModelConfig, only: QtyStuff_T
    use Fov_Convolve_m, only: Convolve_Support_T, &
      & FOV_Convolve_3d
    use HessianModule_0, only: DUMP
    use HessianModule_1, only: HESSIAN_T
    use Load_sps_data_m, only: Grids_T, Dump
    use MLSFillValues, only: isNaN
    use MLSKinds, only: R8, RP, RV
    use MatrixModule_1, only: FINDBLOCK
    use MLSStringLists, only: switchDetail
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use output_m, only: outputNamedValue, resumeOutput, suspendOutput
    use TOGGLES, only: switches
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(r8), intent(in) :: SbRatio    ! Sideband ratio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: d2I_df2(:,:,:) ! mixing ratio derivatives or any
    !                                    parameter for which a simple
    !                                    convolution will suffice

    ! Outputs
    type (Hessian_t), intent(inout) :: Hessian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Local variables
    integer :: Row, Col1, Col2
    integer :: err
    integer :: JF_I, JF_J
    integer :: KI, KJ
    integer :: NFZ_I, NFZ_J ! # frequency elements times # z elements in species
    integer :: NoChans
    integer :: SPS_I, SPS_J ! species indices
    integer :: q, r         ! indices of final elements in state vectors
    real(r8) :: d2rad_df2_out( size(convolve_support%del_chi_out), &
                             & size(d2i_df2,dim=2), size(d2i_df2,dim=3))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans

    ! do the convolution
    if ( any(isNaN(d2i_df2)) ) then
      call dump( d2i_df2, 'd2i_df2', options='-H' )
      call MLSMessage( MLSMSG_Warning, ModuleName // 'Convolve_Other_Second_Deriv', &
        & 'NaNs found in d2i_df2' )
    endif

    call fov_convolve_3d ( convolve_support, d2i_df2, MIF_Times_for_MAF, &
                         & DeadTime, grids_f%deriv_flags, d2rad_df2_out )

    ! load second derivatives into Hessian
    ! First, find index location in Hessian and set the derivative flags

    row = FindBlock( Hessian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    if ( switchDetail( switches, 'hnan', options='-fc' ) > -1 ) then
      call outputNamedValue( 'MAF', MAF )
      call outputNamedValue( 'Channel', Channel )
      call outputNamedValue( 'row', row )
      call outputNamedValue( 'noChans', noChans )
      call outputNamedValue( 'size(qtys)', size(qtys) )
      if ( .not. update) then
        call dump( Grids_f, details=2 )
      endif
    endif

    do sps_i = 1, size(qtys)

      if ( .not. qtys(sps_i)%foundInFirst ) cycle

      q = grids_f%l_v(sps_i-1)
      nfz_i = (Grids_f%l_f(sps_i) - Grids_f%l_f(sps_i-1)) * &
            & (Grids_f%l_z(sps_i) - Grids_f%l_z(sps_i-1))

      do jf_i = Grids_f%windowStart(sps_i), Grids_f%windowfinish(sps_i)

        col1 = FindBlock ( Hessian%col, qtys(sps_i)%qty%index, jf_i)

        do ki = 1, nfz_i

          q = q + 1

          do sps_j = 1, size(qtys)

            if ( .not. qtys(sps_j)%foundInFirst ) cycle

            r = grids_f%l_v(sps_j-1)
            nfz_j = (Grids_f%l_f(sps_j) - Grids_f%l_f(sps_j-1)) * &
                  & (Grids_f%l_z(sps_j) - Grids_f%l_z(sps_j-1))

            do jf_j = Grids_f%windowStart(sps_j), Grids_f%windowfinish(sps_j)

              col2 = FindBlock ( Hessian%col, qtys(sps_j)%qty%index, jf_j)
        
              call getFullBlock_Hessian ( hessian, row, col1, col2, 'atmospheric' )

              do kj = 1, nfz_j
            
                r = r + 1

                ! load derivatives for this (zeta & phi) if needed:

                if ( Grids_f%deriv_flags(q) .and. Grids_f%deriv_flags(r) ) then

                  call loadMatrixValue ( d2rad_df2_out(:,q,r), &
                  & hessian%block(row,col1,col2)%values(channel::noChans,ki,kj), &
                  & sbRatio, update, ERR )

                  select case ( err ) 
                  case (1)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out bigger than hessian slice' )
                  case (2)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out smaller than hessian slice' )
                  case (4)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'NaNs found in d2rad_df2_out' )
                  case (5)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out bigger than hessian slice; NaNs, too' )
                  case (6)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out smaller than hessian slice; NaNs, too' )
                  case default
                  end select
                endif
              end do
              if ( switchDetail( switches, 'hnan', options='-fc' ) > -1) then
                call suspendOutput
                call dump( hessian%block(row,col1,col2), details=-3 )
                call resumeOutput
              endif
   
            end do
        
          end do

        end do

      end do

    end do

  end subroutine Convolve_Other_Second_Deriv

  ! ---------------------------------------  Interpolate_Radiance  -----
  subroutine Interpolate_Radiance ( Coeffs, MAF, Channel, Chi_In, Rad_In, &
           & SbRatio, Update, Ptan, Chi_Out, Radiance, MIF_Times, DeadTime, &
           & Jacobian, RowFlags, dh_dz_out, dx_dh_out, ptan_Der )

    ! Interpolate the radiance from Chi_In to Chi_Out, and maybe dI/dPTan too.

    use MatrixModule_1, only: FINDBLOCK, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use MLSNumerics, only: Coefficients => Coefficients_r8, InterpolateValues
    use ScanAverage_m, only: ScanAverage
    use VectorsModule, only: VectorValue_T

    ! Required inputs
    type(coefficients), intent(in) :: Coeffs
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Chi_In(:)  ! input pointing angles radians
    real(rp), intent(in) :: Rad_In(:)   ! input radiances
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update       ! "add to radiance, don't overwrite"
    type(vectorvalue_t), intent(in) :: PTAN ! Only for some indices
    real(rp), intent(in) :: Chi_Out(:) ! output pointing angles radians

    ! output
    type (VectorValue_T), intent(inout) :: RADIANCE   ! Output radiances
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.

    ! Optional if PTan derivative processed.  All or none.
    type (Matrix_t), optional, intent(inout) :: Jacobian
    logical, optional, intent(inout) :: rowFlags(:) ! Flag to calling code
    real(rp), optional, intent(in) :: dh_dz_out(:) ! dh/dz on the output pointing grid
    real(rp), optional, intent(in) :: dx_dh_out(:) ! dx/dh on the output pointing grid
    logical, optional, intent(in) :: Ptan_der ! "Process PTAN derivatives"

    logical :: my_ptan_der
    integer :: Col, NoChans, NoPtan, Row
    real(r8), dimension(ptan%template%noSurfs) :: SRad, dI_dx, Rad_Out
    real(rv), pointer :: MIF_Times_for_MAF(:)

    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)

    my_ptan_der = .false.
    if ( present ( ptan_der ) ) my_ptan_der = ptan_der

    noPtan = ptan%template%noSurfs
    noChans = radiance%template%noChans

    if ( my_ptan_der ) then
      row = FindBlock ( Jacobian%row, radiance%index, maf )
      rowFlags(row) = .TRUE.
      col = FindBlock ( Jacobian%col, ptan%index, maf )

      if ( associated(MIF_Times) ) then
        call scanAverage ( MIF_Times_for_MAF, deadTime(1,1), &
          & chi_in, chi_out, rad_in, rad_out, dY_dX_out=di_dx )
      else
        call InterpolateValues ( coeffs, chi_in, rad_in, chi_out, rad_out, &
                               & METHOD='S', extrapolate='C', dyByDx=di_dx )
      end if

      ! Use the chain rule to compute dI/dPtan on the output grid:

      SRad = di_dx * dx_dh_out * dh_dz_out
      call getBandedBlock ( jacobian, row, col, noChans, noPtan )
      call loadMatrixValue ( sRad, Jacobian%block(row,col)%values(channel::noChans,1), &
                           & sbRatio, update )

    else

        if ( associated(MIF_Times) ) then
          call scanAverage ( MIF_Times_for_MAF, deadTime(1,1), &
            & chi_in, chi_out, rad_in, rad_out )
        else
          call InterpolateValues ( coeffs, chi_in, rad_in, chi_out, rad_out, &
                                 & METHOD='S', extrapolate='C' )
        end if

    end if

    ! Load the Radiance values into the Radiance structure:

    call loadVectorValue ( rad_out, Radiance%values(channel::noChans,maf), &
                         & sbRatio, update )

  end subroutine Interpolate_Radiance

  ! ------------------------------  Interpolate_Temperature_Deriv  -----
  subroutine Interpolate_Temperature_Deriv ( Coeffs, MAF, Channel, Chi_In, &
           & SbRatio, Update, Chi_Out, Radiance, Temp, Grids_Tmp, &
           & MIF_Times, DeadTime, di_dT, Jacobian, RowFlags )

    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use MLSNumerics, only: Coefficients => Coefficients_r8, InterpolateValues
    use ScanAverage_m, only: ScanAverage
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(coefficients), intent(in) :: Coeffs
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Chi_In(:)  ! input pointing angles radians
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    real(rp), intent(in) :: Chi_Out(:) ! output pointing angles radians
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type (vectorvalue_t), intent(in) :: TEMP     ! Only for some indices
    type (Grids_T), intent(in) :: Grids_Tmp      ! Temperature's grids, etc.
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                      temperature on chi_in

    ! Outputs
    type (Matrix_t), intent(inout) :: Jacobian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Local variables
    integer :: Col, JF, JZ, K, NoChans, Row, SV_I
    real(r8) :: dRad_dT_in(size(chi_in)), dRad_dT_out(size(chi_out))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    ! Start here

    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)

    sv_i = 0
    k = size(chi_in)
    noChans = radiance%template%noChans

    row = FindBlock ( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    do jf = grids_tmp%windowStart(1), grids_tmp%windowFinish(1)

      col = FindBlock ( Jacobian%col, temp%index, jf )
      call getFullBlock ( jacobian, row, col, 'temperature' )

      do jz = 1, temp%template%noSurfs

        ! Check if derivatives are needed for this (zeta & phi) :

        sv_i = sv_i + 1
        if ( grids_tmp%deriv_flags(sv_i) ) then
          dRad_dT_in = di_dt(1:k,sv_i)
          if ( associated(MIF_Times) ) then
            call scanAverage ( MIF_Times_for_MAF, deadTime(1,1), &
              & chi_in, chi_out, dRad_dT_in, dRad_dT_out )
          else
            call InterpolateValues ( coeffs, chi_in, dRad_dT_in, &
                                   & chi_out, dRad_dT_out, &
                                   & method = 'S', extrapolate = 'C' )
          end if
          call loadMatrixValue ( dRad_dT_out, &
            & jacobian%block(row,col)%values(channel::noChans,jz), sbRatio, &
            & update )
        end if

      end do

    end do

  end subroutine Interpolate_Temperature_Deriv

  ! ------------------------------------  Interpolate_Other_Deriv  -----
  subroutine Interpolate_Other_Deriv ( Coeffs, MAF, Channel, Chi_In, &
             & SbRatio, Update, Chi_Out, Radiance, Qtys, Grids_f, &
             & MIF_Times, DeadTime, dI_df, Jacobian, RowFlags, Linear, &
             & Derivs, ExtraJacobian )

    use ForwardModelConfig, only: QtyStuff_T
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: R8, RP, RV
    use MLSNumerics, only: Coefficients => Coefficients_r8, InterpolateValues
    use ScanAverage_m, only: ScanAverage
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(coefficients), intent(in) :: Coeffs
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    real(rp), intent(in) :: Chi_In(:)  ! input pointing angles radians
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    real(rp), intent(in) :: Chi_Out(:) ! output pointing angles radians
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: dI_df(:,:) ! mixing ratio derivatives or any
    !                                    parameter for which a simple
    !                                    convolution will suffice

    ! Outputs
    type (Matrix_t), intent(inout), target :: Jacobian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Optional inputs
    logical, intent(in), optional :: Linear ! "Use linear interpolation"
                                            ! Deafult true.
    logical, intent(in), optional :: Derivs(:) ! Same size as Qtys; if not
                                               ! present, do all columns.  Use
                                               ! in case grids_f%deriv_flags
                                               ! hasn't been set.
    ! Optional outputs
    type (Matrix_t), intent(inout), optional, target :: ExtraJacobian

    ! Local variables
    integer :: Col, IS, JF, K, NFZ, NoChans, Row, SV_F
    type (Matrix_t), pointer :: MyJacobian
    logical :: MyLinear

    real(r8) :: dRad_df_in(size(chi_in)), dRad_df_out(size(chi_out))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)

    myLinear = .true.
    if ( present(linear) ) myLinear = linear

    noChans = Radiance%template%noChans

    row = FindBlock ( Jacobian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    do is = 1, size(qtys)

      if ( present(derivs) ) then
        if ( .not. derivs(is) ) cycle
      end if

      if ( qtys(is)%foundInFirst ) then
        myJacobian => jacobian
      else
        if ( .not. present(extraJacobian) ) cycle
        myJacobian => extraJacobian
      end if

      sv_f = grids_f%l_v(is-1)
      nfz = (Grids_f%l_f(is) - Grids_f%l_f(is-1)) * &
          & (Grids_f%l_z(is) - Grids_f%l_z(is-1))

      do jf = Grids_f%windowStart(is), Grids_f%windowfinish(is)

        col = FindBlock ( myJacobian%col, qtys(is)%qty%index, jf)
        call getFullBlock ( myJacobian, row, col, 'atmospheric' )

        do k = 1, nfz

          ! Check if derivatives are needed for this (zeta & phi) :

          sv_f = sv_f + 1
          if ( Grids_f%deriv_flags(sv_f) ) then
            dRad_df_in = di_df(:,sv_f)
            if ( associated(MIF_Times) ) then
              call scanAverage ( MIF_Times_for_MAF, deadTime(1,1), &
                & chi_in, chi_out, dRad_df_in, dRad_df_out )
            else if ( myLinear ) then
              call InterpolateValues ( chi_in, dRad_df_in, &
                                     & chi_out, dRad_df_out, &
                                     & method = 'L', extrapolate = 'C' )
            else
              call InterpolateValues ( coeffs, chi_in, dRad_df_in, &
                                     & chi_out, dRad_df_out, &
                                     & method = 'S', extrapolate = 'C' )
            end if
            call loadMatrixValue ( dRad_df_out, &
              & myJacobian%block(row,col)%values(channel::noChans,k), &
              & sbRatio, update )
          end if

        end do

      end do

    end do

  end subroutine Interpolate_Other_Deriv

  ! ---------------------------------  Store_Temperature_Deriv_1D  -----
  subroutine Store_Temperature_Deriv_1d ( MAF, Row_0, Radiance, Temp, &
                                        & Grids_Tmp, di_dT, Jacobian )

    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: RP
    use VectorsModule, only: VectorValue_T

    ! Inputs
    integer, intent(in) :: MAF
    integer, intent(in) :: Row_0 ! within the block, (channel-1)*noChans+zeta
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type (Vectorvalue_t), intent(in) :: TEMP     ! Only for some indices
    type (Grids_T), intent(in) :: Grids_Tmp      ! Temperature's grids, etc.
    real(rp), intent(in) :: dI_dT(:) ! derivative of radiance in Row_0 wrt
    !                                  temperature state vector coordinates

    ! Outputs
    type (Matrix_t), intent(inout) :: Jacobian

    ! Local variables
    integer :: Col, JF, JZ, Row_1, SV_I

    ! Start here

    sv_i = 0

    row_1 = FindBlock ( Jacobian%row, radiance%index, maf )

    do jf = grids_tmp%windowStart(1), grids_tmp%windowFinish(1)

      col = FindBlock ( Jacobian%col, temp%index, jf )
      call getFullBlock ( jacobian, row_1, col, 'temperature' )

      do jz = 1, temp%template%noSurfs

        ! Check if derivatives are needed for this (zeta & phi) :

        sv_i = sv_i + 1
        if ( grids_tmp%deriv_flags(sv_i) ) &
          & jacobian%block(row_1,col)%values(row_0,jz) = dI_dT(sv_i)

      end do

    end do

  end subroutine Store_Temperature_Deriv_1D

  ! ------------------------------------  Store_Temperature_Deriv_2D  -----
  subroutine Store_Temperature_Deriv_2D ( MAF, Channel, Radiance, Temp, &
                                        & Grids_Tmp, di_dT, Jacobian )

    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: RP
    use VectorsModule, only: VectorValue_T

    ! Inputs
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type (vectorvalue_t), intent(in) :: TEMP     ! Only for some indices
    type (Grids_T), intent(in) :: Grids_Tmp      ! Temperature's grids, etc.
    real(rp), intent(in) :: dI_dT(:,:) ! derivative of radiance wrt
    !                      temperature on chi_in

    ! Outputs
    type (Matrix_t), intent(inout) :: Jacobian

    ! Local variables
    integer :: Col, JF, JZ, NoChans, Row, SV_I

    ! Start here

    sv_i = 0
    noChans = radiance%template%noChans

    row = FindBlock ( Jacobian%row, radiance%index, maf )

    do jf = grids_tmp%windowStart(1), grids_tmp%windowFinish(1)

      col = FindBlock ( Jacobian%col, temp%index, jf )
      call getFullBlock ( jacobian, row, col, 'temperature' )

      do jz = 1, temp%template%noSurfs

        ! Check if derivatives are needed for this (zeta & phi) :

        sv_i = sv_i + 1
        if ( grids_tmp%deriv_flags(sv_i) ) &
          & jacobian%block(row,col)%values(channel::noChans,jz) = dI_dT(:,sv_i)

      end do

    end do

  end subroutine Store_Temperature_Deriv_2D

  ! ---------------------------------------  Store_Other_Deriv_1D  -----
  subroutine Store_Other_Deriv_1D ( MAF, Row_0, Radiance, Qtys, Grids_f, &
                                  & dI_df, Jacobian, ExtraJacobian, Derivs )

    use ForwardModelConfig, only: QtyStuff_T
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: RP
    use VectorsModule, only: VectorValue_T

    ! Inputs
    integer, intent(in) :: MAF
    integer, intent(in) :: Row_0 ! within the block, (channel-1)*noChans+zeta
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rp), intent(in) :: dI_df(:) ! mixing ratio derivatives or any
    !                                  parameter for which a simple
    !                                  convolution will suffice

    ! Output
    type (Matrix_t), intent(inout), target :: Jacobian

    ! Optional output
    type (Matrix_t), intent(inout), optional, target :: ExtraJacobian

    ! Optional input
    logical, intent(in), optional :: Derivs(:) ! Same size as Qtys; if not
                                               ! present, do all columns.  Use
                                               ! in case grids_f%deriv_flags
                                               ! hasn't been set.

    ! Local variables
    integer :: Col, IS, JF, K, NFZ, Row_1, SV_F
    type (Matrix_t), pointer :: MyJacobian

    row_1 = FindBlock ( Jacobian%row, radiance%index, maf )

    do is = 1, size(qtys)

      if ( present(derivs) ) then
        if ( .not. derivs(is) ) cycle
      end if

      if ( qtys(is)%foundInFirst ) then
        myJacobian => jacobian
      else
        if ( .not. present(extraJacobian) ) cycle
        myJacobian => extraJacobian
      end if

      sv_f = grids_f%l_v(is-1)
      nfz = (Grids_f%l_f(is) - Grids_f%l_f(is-1)) * &
          & (Grids_f%l_z(is) - Grids_f%l_z(is-1))

      do jf = Grids_f%windowStart(is), Grids_f%windowfinish(is)

        col = FindBlock ( myJacobian%col, qtys(is)%qty%index, jf)
        call getFullBlock ( myJacobian, row_1, col, 'atmospheric' )

        do k = 1, nfz

          ! Check if derivatives are needed for this (zeta & phi) :

          sv_f = sv_f + 1
          if ( Grids_f%deriv_flags(sv_f) ) &
            & myJacobian%block(row_1,col)%values(row_0,k) = di_df(sv_f)

        end do

      end do

    end do

  end subroutine Store_Other_Deriv_1D

  ! ---------------------------------------  Store_Other_Deriv_2D  -----
  subroutine Store_Other_Deriv_2D ( MAF, Channel, Radiance, Qtys, Grids_f, &
                                  & dI_df, Jacobian, ExtraJacobian, Derivs )

    use ForwardModelConfig, only: QtyStuff_T
    use Load_sps_data_m, only: Grids_T
    use MatrixModule_1, only: FINDBLOCK, GetFullBlock, MATRIX_T
    use MLSKinds, only: RP
    use VectorsModule, only: VectorValue_T

    ! Inputs
    integer, intent(in) :: MAF
    integer, intent(in) :: CHANNEL
    type (VectorValue_T), intent(in) :: RADIANCE ! Only for some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rp), intent(in) :: dI_df(:,:) ! mixing ratio derivatives or any
    !                                    parameter for which a simple
    !                                    convolution will suffice

    ! Outputs
    type (Matrix_t), intent(inout), target :: Jacobian

    ! Optional outputs
    type (Matrix_t), intent(inout), optional, target :: ExtraJacobian

    ! Optional input
    logical, intent(in), optional :: Derivs(:) ! Same size as Qtys; if not
                                               ! present, do all columns.  Use
                                               ! in case grids_f%deriv_flags
                                               ! hasn't been set.

    ! Local variables
    integer :: Col, IS, JF, K, NFZ, NoChans, Row, SV_F
    type (Matrix_t), pointer :: MyJacobian

    noChans = Radiance%template%noChans

    row = FindBlock ( Jacobian%row, radiance%index, maf )

    do is = 1, size(qtys)

      if ( present(derivs) ) then
        if ( .not. derivs(is) ) cycle
      end if

      if ( qtys(is)%foundInFirst ) then
        myJacobian => jacobian
      else
        if ( .not. present(extraJacobian) ) cycle
        myJacobian => extraJacobian
      end if

      sv_f = grids_f%l_v(is-1)
      nfz = (Grids_f%l_f(is) - Grids_f%l_f(is-1)) * &
          & (Grids_f%l_z(is) - Grids_f%l_z(is-1))

      do jf = Grids_f%windowStart(is), Grids_f%windowfinish(is)

        col = FindBlock ( myJacobian%col, qtys(is)%qty%index, jf)
        call getFullBlock ( myJacobian, row, col, 'atmospheric' )

        do k = 1, nfz

          ! Check if derivatives are needed for this (zeta & phi) :

          sv_f = sv_f + 1
          if ( Grids_f%deriv_flags(sv_f) ) &
            & myJacobian%block(row,col)%values(channel::noChans,k) = di_df(:,sv_f)

        end do

      end do

    end do

  end subroutine Store_Other_Deriv_2D

! =====     Private Procedures     =====================================

  subroutine GetBandedBlock ( Jacobian, Row, Col, NoChans, NoPtan )
    use MatrixModule_0, only: M_ABSENT, M_BANDED, CHECKFORSIMPLEBANDEDLAYOUT
    use MatrixModule_1, only: CREATEBLOCK, MATRIX_T
     use MLSKinds, only: RM
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use String_Table, only: Get_String
    type (Matrix_t), intent(inout) :: Jacobian
    integer, intent(in) :: Row, Col, NoChans, NoPtan
    character(len=63) :: ForWhom
    if ( jacobian%name /= 0 ) then
      call get_string ( jacobian%name, forWhom )
      forWhom = trim(forWhom) // " in GetBandedBlock"
    else
      forWhom = "GetBandedBlock"
    end if
    select case (jacobian%block(Row,col)%kind)
      case (m_absent)
        call CreateBlock ( Jacobian, row, col, m_banded, noPtan*noChans, &
                         & bandHeight=noChans, init=0.0_rm, forWhom=forWhom )
      case (m_banded)
        call CheckForSimpleBandedLayout ( jacobian%block(row,col), noChans, &
          & 'd[radiance]/d[ptan] in convolution' )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Wrong matrix block type for ptan derivative' )
    end select
  end subroutine GetBandedBlock

  subroutine GetFullBlock_Hessian ( Hessian, Row, Col1, Col2, What )
    use HessianModule_0, only: H_ABSENT, H_FULL, RH
    use HessianModule_1, only: CREATEBLOCK, HESSIAN_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    type (Hessian_t), intent(inout) :: Hessian
    integer, intent(in) :: Row, Col1, Col2
    character(len=*), intent(in) :: What

    select case ( Hessian%block(row,col1,col2)%kind )
      case ( h_absent )
        ! If profiling shows we're spending too much time filling the block,
        ! we could send in a signal and only fill for the channels for which
        ! we are not computing results.
        call CreateBlock ( Hessian, row, col1, col2, h_full, inittuples=0, fill=0.0_rh )
      case ( h_full )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Wrong matrix block type for ' // what // ' second derivative matrix' )
    end select

  end subroutine GetFullBlock_Hessian


  subroutine LoadMatrixValue ( In, Out, SbRatio, Update, err )
    ! If Update add SbRatio*In to Out else assign SbRatio*In to Out.
    ! Identical to LoadVectorValue except for kind of Out.
    use MLSFillValues, only: isNaN
    use MLSKinds, only: R8, RM
    real(r8), intent(in) :: In(:)
    real(rm), intent(inout) :: Out(:)
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update
    integer, optional, intent(out) :: err
    integer :: myErr ! poss. values: 0, 1, 2, 4, 5, 6
    integer :: N
    n = size(in)
    myErr = 0
    if ( n > size(out) ) then
      myErr = 1
    endif
    if ( n < size(out) ) then
      myErr = myErr + 2
    endif
    if ( any(isNaN(in)) ) then
      myErr = myErr + 4
    endif
    if ( present(err) ) err = myErr
    if ( any(myErr == (/1, 4, 5, 6/)) ) return
    if ( update ) then
      out(:n) = out(:n) + sbRatio * in
    else
      out(:n) = sbRatio * in
    end if
  end subroutine LoadMatrixValue

  subroutine LoadVectorValue ( In, Out, SbRatio, Update )
    ! If Update add SbRatio*In to Out else assign SbRatio*In to Out.
    ! Identical to LoadMatrixValue except for kind of Out.
    use MLSKinds, only: R8, RV
    real(r8), intent(in) :: In(:)
    real(rv), intent(inout) :: Out(:)
    real(r8), intent(in) :: SbRatio
    logical, intent(in) :: Update
    integer :: N
    n = size(in)
    if ( update ) then
      out(:n) = out(:n) + sbRatio * in
    else
      out(:n) = sbRatio * in
    end if
  end subroutine LoadVectorValue

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Convolve_All_m

! $Log$
! Revision 2.23  2013/08/02 01:24:06  vsnyder
! Add ExtraJacobian to compute derivatives not in state vector
!
! Revision 2.22  2013/06/12 02:19:16  vsnyder
! Cruft removal
!
! Revision 2.21  2012/07/10 04:10:59  vsnyder
! Remove noChans from Store_*_Deriv_1D because it's not used
!
! Revision 2.20  2012/07/07 00:14:33  vsnyder
! Shorten some comments to avoid gripes about long lines
!
! Revision 2.19  2012/07/06 21:30:31  yanovsky
! Added Convolve_Radiance_Normalization and
! Convolve_Temperature_Deriv_Normalization subroutines that compute
! normalized Temperature derivatives
!
! Revision 2.18  2011/12/17 00:35:29  vsnyder
! Move GetFullBlock to MatrixModule_1
!
! Revision 2.17  2011/08/20 00:44:32  vsnyder
! Remove unused USE statements
!
! Revision 2.16  2010/11/09 01:01:25  vsnyder
! Get RH from HessianModule_0 instead of RM from MLSKinds.  Delete unused
! declarations.  Remove explicit copy-out for loadMatrixValue (should have
! been copy-in/copy-out anyway).  Don't get strings that aren't needed.
!
! Revision 2.15  2010/10/13 22:18:22  pwagner
! Intermediate steps in eliminating NaNs from analytic Hessians
!
! Revision 2.14  2010/08/27 05:53:23  yanovsky
! Fixed the order of do loops in Convolve_Other_Second_Deriv subroutine.  Added comments.
!
! Revision 2.13  2010/07/27 01:12:56  vsnyder
! Fold some lines the compiler said were too long
!
! Revision 2.12  2010/07/18 23:39:19  yanovsky
! Add Convolve_Other_Second_Deriv and GetFullBlock_Hessian
!
! Revision 2.11  2010/02/05 03:18:11  vsnyder
! Remove USE for unreferenced names
!
! Revision 2.10  2009/12/22 03:23:05  vsnyder
! Add Store_Other_Deriv, Store_Temp_Deriv versions
!
! Revision 2.9  2009/12/15 03:17:07  vsnyder
! Get kinds from MLSKinds instead of MLSCommon
!
! Revision 2.8  2009/11/17 23:40:08  vsnyder
! Add Store_*_Deriv routines
!
! Revision 2.7  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.6  2008/08/27 19:56:51  vsnyder
! Add PRINT to not_used_here
!
! Revision 2.5  2006/08/02 19:55:33  vsnyder
! Tell CreateBlock that Convolve creates it, for leak tracking
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
! Revision 2.1  2005/07/06 02:16:54  vsnyder
! Initial commit, replacing convolve_all_m.f90 and no_conv_at_all.f90
!
