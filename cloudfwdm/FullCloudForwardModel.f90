! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullCloudForwardModel

  ! This module will contain the full cloud forward model.  Currently it only
  ! contains a wrapper, but eventually will contain the top level routine for
  ! the model itself.

  use Allocate_deallocate, only: Allocate_test, Deallocate_test
  use MLSCommon, only: r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use VectorsModule, only: GETVECTORQUANTITYBYTYPE, VECTOR_T, VECTORVALUE_T
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, &
    & FORWARDMODELSTATUS_T
  use MatrixModule_1, only: MATRIX_T
  use Intrinsic, only: L_TEMPERATURE, L_PTAN, L_VMR, L_GPH
  use ManipulateVectorQuantities, only: FindClosestInstances

  implicit none
  private

  public :: FullCloudForwardModelWrapper

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Local parameters ---------------------------------------------------------

  character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

  
contains ! ============= Public Procedures ==========================

  subroutine FullCloudForwardModelWrapper ( ForwardModelConfig, FwdModelIn, FwdModelExtra, &
        FwdModelOut, Ifm, fmStat, Jacobian )
    ! This is simply a wrapper for Jonathan's full up forward model.  Hopefully
    ! we won't need to wrap it in the long run, but will call it directly.
    
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Local variables

    type (VectorValue_T), pointer :: GPH ! Geopotential height sv qty
    type (VectorValue_T), pointer :: PTAN ! Tangent pressure
    type (VectorValue_T), pointer :: RADIANCE ! Radiance quantity
    type (VectorValue_T), pointer :: TEMP ! Temperature state vector quantity
    type (VectorValue_T), pointer :: VMR ! One vmr quantity

    type (Signal_T) :: signal           ! The signal we're working on.

    real(r8), dimension(:,:), pointer :: VMRARRAY ! The vmr's

    integer :: i                        ! Loop counter
    integer :: MAF                      ! The major frame we're working on.
    integer :: VMRINST                  ! Instance index
    integer :: INSTANCE                 ! Relevant instance for temperature
    integer :: GPHINST                  ! Relevant instance for GPH

    integer, dimension(:), pointer :: closestInstances ! From find routine

    real(r8), dimension(:,:), pointer :: A_CLEARSKYRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDINDUCEDRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_CLOUDSENSITIVITY
    real(r8), dimension(:,:), pointer :: A_EFFECTIVEOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_MASSMEANDIAMETER
    real(r8), dimension(:,:), pointer :: A_TOTALOPTICALDEPTH

    real(r8), dimension(:), pointer :: FREQUENCIES

    ! Executable code

    ! Check the configuration we've been given, and deduce some stuff
    if ( size ( forwardModelConfig%signals ) /= 1 ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Cannot call the full cloud forward model with multiple signals' )
    signal = forwardModelConfig%signal(1)
    maf = fmStat%maf

    ! For the moment make it only single sideband
    if ( signal%sideband == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
      & 'Only single sidebands allowed in FullForwardCloudModel for now' )
    call Allocate_test ( frequencies, count ( signal%channels ), &
      & 'frequencies', ModuleName )
    frequencies = signal%lo + signal%sideband * ( signal%centerFrequency + &
      & pack ( signal%frequencies, signal%channels ) )

    ! Get the quantities we need from thevectors
    radiance => GetVectorQuantityByType ( fwdModelOut, &
      & quantityType=l_radiance, signal=signal%index, sideband=signal%sideband )
    cloudInducedRadiance => GetVectorQuantityByType ( fwdModelOut, &
      & quantityType=l_cloudInducedRadiance, &
      & signal=signal%index, sideband=signal%sideband )

    ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_ptan, radiometer=signal%radiometer )
    temp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_temperature )
    gph => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_gph )
    ! Do the vmr's one by one later on.

    ! Make sure the quantities we have are OK
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(gph, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )

    ! Set up some temporary quantities
    call Allocate_test ( closestInstances, radiance%template%noInstances, &
      & 'closestInstances', ModuleName )      

    ! Assemble the vmr array
    call allocate_test ( vmrArray, &
      & size(forwardModelConfig%molecules), temp%template%noSurfs, &
      & 'vmrArray', ModuleName )

    do i = 1, size(forwardModelConfig%molecules)
      vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_vmr, molecule=forwardModelConfig%molecules(i) )
      if (.not. ValidateVectorQuantity ( vmr,  stacked=.true., coherent=.true., &
        & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'vmr' )

      call FindClosestInstances ( vmr, radiance, closestInstances )
      vmrInst = closestInstances(maf)
      call InterpolateValues ( &
        & vmr%template%surfs, &         ! Old X
        & vmr%values(:,vmrInst), &         ! Old Y
        & temp%template%surfs, &        ! New X
        & vmrArray(i,:), &              ! New Y
        & method='Linear', extrapolate='Clamp' )
    end do

    ! Work out the closest instances for the other quantities
    call FindClosestInstances ( temp, radiance, closestInstances )
    instance = closestInstances(maf)
    noLayers = temp%template%noSurfs - 1

    ! Make temporary arrays for Jonathan's stuff
    call Allocate_test ( a_clearSkyRadiance, &
      & radiance%template%noSurfs, noFreqs, &
      & 'a_clearSkyRadiance', ModuleName )
    call Allocate_test ( a_cloudInducedRadiance, &
      & radiance%template%noSurfs, noFreqs, &
      & 'a_cloudInducedRadiance', ModuleName )
    call Allocate_test ( a_effectiveOpticalDepth, &
      & radiance%template%noSurfs, noFreqs, &
      & 'a_effectiveOpticalDepth', ModuleName )
    call Allocate_test ( a_cloudSensitivity, &
      & radiance%template%noSurfs, noFreqs, &
      & 'a_cloudSensitivity', ModuleName )
    call Allocate_test ( a_totalOpticalDepth, &
      & temp%template%noSurfs, noFreqs, &
      & 'a_totalOpticalDepth', ModuleName )
    call Allocate_test ( a_cloudOpticalDepth, &
      & temp%template%noSurfs, noFreqs, &
      & 'a_cloudOpticalDepth', ModuleName )
    call Allocate_test ( a_massMeanDiameter, &
      & 2, temp%template%noSurfs, &
      & 'a_massMeanDiameterModuleName', ModuleName )
    
    ! Now call Jonathan's f77 code

    call CloudForwardModel ( &
      & frequency/1e3_r8, &
      & 10.0^(-temp%template%surfs), &
      & gph%values(:, instance), &
      & temp%values(:,instance), &
      & vmrArray, &
      & 10.0^(-ptan%values(:,maf)), &
      & size(frequency), &
      & temp%template%noSurfs, &
      & radiance%template%noSurfs, &
      & WCscale%values(1,wcScaleInst), &
      & cloudType%values(1,instance), &
      & PSDtype%values(1,instance), &
      & surfaceType%values(1, instance), &
      & forwardModelConfig%cloud_der, &
      & saturatedBelowClouds%values(1, instance), &
      & cloudHeight%values(1, instance), &
      & a_clearSkyRadiance, &
      & a_cloudInducedRadiance, &
      & a_totalOpticalDepth, &
      & a_cloudOpticalDepth, &
      & layerThickness%value( 1:noLayers, instance ), &
      & a_massMeanDiameter, &
      & a_effectiveOpticalDepth, &
      & cloudSensitivity )

    effectiveOpticalDeptmassMeanDiameterIce%values(:,instance) = 0.0_r8

    ! Now store results in relevant vectors
    ! Vectors are stored ( noChannels * noSurfaces, noInstances ), so transpose
    ! and reshape temporary variables to be in the right form.
    ! First the minor frame stuff
    radiance%values ( :, maf) = &
      & reshape ( transpose(a_clearSkyRadiance), &
      & (/clearSkyRadiance%template%instanceLen/) )
    cloudInducedRadiance%values ( :, maf ) = &
      & reshape ( transpose(a_cloudInducedRadiance), &
      & (/cloudInducedRadiance%template%instanceLen/) )
    effectiveOpticalDepth%values ( :, maf ) = &
      & reshape ( transpose(a_effectiveOpticalDepth), &
      & (/effectiveOpticalDepth%template%instanceLen/) )
    cloudSensitivity%values ( :, maf ) = &
      & reshape ( transpose(a_cloudSensitivity), &
      & (/cloudSensitivity%template%instanceLen/) )

    ! For layer (noTempSurfs-1) stuff make sure all are zero to start, then do rest
    cloudOpticalDepth%values(:,instance) =     0.0_r8
    layerThickness%value( :, instance) =       0.0_r8
    massMeanDiameterIce%values(:,instance) =   0.0_r8
    massMeanDiameterWater%values(:,instance) = 0.0_r8
    totalOpticalDepth%values(:,instance) =     0.0_r8
    
    cloudOpticalDepth%values ( 1:noLayers, instance ) = &
      & reshape ( transpose(a_cloudOpticalDepth), (/noLayers*noFreqs/) )
    massMeanDiameterIce%values ( 1:noLayers, instance ) = &
      & reshape ( transpose(a_massMeanDiameterIce), (/noLayers*noFreqs/) )
    massMeanDiameterWater%values ( 1:noLayers, instance ) = &
      & reshape ( transpose(a_massMeanDiameterWater), (/noLayers*noFreqs/) )
    totalOpticalDepth%values ( 1:noLayers, instance ) = &
      & reshape ( transpose(a_totalOpticalDepth), (/noLayers*noFreqs/) )

    ! Remove temporary quantities
    call Deallocate_test ( a_massMeanDiameter, 'a_massMeanDiameterModuleName', ModuleName )
    call Deallocate_test ( a_cloudOpticalDepth, 'a_cloudOpticalDepth', ModuleName )
    call Deallocate_test ( a_totalOpticalDepth, 'a_totalOpticalDepth', ModuleName )
    call Deallocate_test ( a_cloudSensitivity, 'a_cloudSensitivity', ModuleName )
    call Deallocate_test ( a_effectiveOpticalDepth, 'a_effectiveOpticalDepth', ModuleName )
    call Deallocate_test ( a_cloudInducedRadiance, 'a_cloudInducedRadiance', ModuleName )
    call Deallocate_test ( a_clearSkyRadiance, 'a_clearSkyRadiance', ModuleName )
    call Deallocate_test ( vmrArray, 'vmrArray', ModuleName )
    call Deallocate_test ( closestInstances, 'closestInstances', ModuleName )

  end subroutine FullCloudForwardModelWrapper

end module FullCloudForwardModel

