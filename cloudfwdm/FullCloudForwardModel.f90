
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FullCloudForwardModel

! -------------------------------------------------------------------------
! THIS MODULE CONTAINS THE FULL CLOUD FORWARD MODEL  
! Jonathan Jiang, Paul Wagner, Jul 16, 2001 
! -------------------------------------------------------------------------

  use Allocate_deallocate, only: Allocate_test, Deallocate_test
  use Hdf,      only: DFACC_READ, DFACC_CREATE
  use HDFEOS,   only: SWOPEN,     SWCLOSE
  use L2GPData, only: L2GPData_T, ReadL2GPData, WriteL2GPData
  use MLSCommon,only: NameLen,    FileNameLen, r8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSSignals_m, only: SIGNAL_T
  use MatrixModule_1, only: MATRIX_T
  use ManipulateVectorQuantities, only: FindClosestInstances
  use MLSNumerics, only: InterpolateValues
  use Molecules, only: L_H2O, L_O3
  use VectorsModule, only: GETVECTORQUANTITYBYTYPE,                          &
                         & VECTOR_T, VECTORVALUE_T,                          &
                         & VALIDATEVECTORQUANTITY
  
! -----------------------------------------------------------------------
! THE FOLLOWING IS MODIFICATIONS FOR THE CLOUD FORWARD MODEL PARAMETERS
! -----------------------------------------------------------------------

  use ForwardModelConfig, only: FORWARDMODELCONFIG_T   
  use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T,            &
                                    & FORWARDMODELSTATUS_T

! ----------------------------------------------------------
! DEFINE INTRINSIC CONSTANTS NEEDED BY Init_Tables_Module
! I'm not sure anything else is needed !
! ----------------------------------------------------------

  use Intrinsic, only: L_TEMPERATURE,L_PTAN,L_VMR,L_GPH,L_RADIANCE,L_NONE,   &
                     & L_CLOUDINDUCEDRADIANCE,                               &
                     & L_EFFECTIVEOPTICALDEPTH,                              &
                     & L_CLOUDRADSENSITIVITY,                                &
                     & L_TOTALEXTINCTION,                                    &
                     & L_CLOUDEXTINCTION,                                    &
                     & L_MASSMEANDIAMETERICE,                                & 
                     & L_MASSMEANDIAMETERWATER,                              &
                     & L_SURFACETYPE,                                        &
                     & L_SIZEDISTRIBUTION,                                   &
                     & L_TNGTGEOCALT,                                        &
                     & L_EARTHRADIUS,                                      &
                     & L_CLOUDICE,                                           &
                     & L_CLOUDWATER 

  implicit none
  private

  public :: FullCloudForwardModelWrapper

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Local parameters ---------------------------------------------------------

  character, parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

         ! ---------------------------------------------------------------------
contains ! THIS SUBPROGRAM CONTAINS THE WRAPPER ROUTINE FOR CALLING THE FULL
         ! CLOUD FORWARD MODEL
         ! ---------------------------------------------------------------------
  subroutine FullCloudForwardModelWrapper ( ForwardModelConfig, FwdModelIn,  &
                                            FwdModelExtra, FwdModelOut, Ifm, &
                                            fmStat, Jacobian                 )  
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT                ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM      ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT         ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Local variables
    type (VectorValue_T), pointer :: CLOUDICE                   ! Profiles
    type (VectorValue_T), pointer :: CLOUDWATER                 ! Profiles
    type (VectorValue_T), pointer :: CLOUDEXTINCTION            ! Profiles
    type (VectorValue_T), pointer :: CLOUDINDUCEDRADIANCE       ! Like radiance
    type (VectorValue_T), pointer :: CLOUDRADSENSITIVITY        ! Like radiance
    type (VectorValue_T), pointer :: EFFECTIVEOPTICALDEPTH      ! Quantity
    type (VectorValue_T), pointer :: GPH                        ! Geop height
    type (VectorValue_T), pointer :: MASSMEANDIAMETERICE        ! Quantity
    type (VectorValue_T), pointer :: MASSMEANDIAMETERWATER      ! Quantity
    type (VectorValue_T), pointer :: PTAN                       ! Tgt pressure
    type (VectorValue_T), pointer :: RADIANCE                   ! Quantity
    type (VectorValue_T), pointer :: SIZEDISTRIBUTION           ! Integer really
    type (VectorValue_T), pointer :: EARTHRADIUS              ! Scalar 
    type (VectorValue_T), pointer :: SURFACETYPE                ! Integer really
    type (VectorValue_T), pointer :: TEMP                       ! Temperature 
    type (VectorValue_T), pointer :: TOTALEXTINCTION            ! Profile
    type (VectorValue_T), pointer :: VMR                        ! Quantity

    type (Signal_T) :: signal               ! I don't know how to define this!

    real(r8), dimension(:,:), pointer :: VMRARRAY               ! The vmr's

    integer :: i                        ! Loop counter
    integer :: ivmr
    integer :: NQ1
    integer :: NQ2
    integer :: MAF                      ! The major frame 
    integer :: VMRINST                  ! Instance index
    integer :: INSTANCE                 ! Relevant instance for temperature
    integer :: GPHINST                  ! Relevant instance for GPH
    integer :: NOLAYERS                 ! temp.noSurfs - 1
    integer :: NOFREQS                  ! Number of frequencies
    integer :: NOsurf                   ! Number of pressure levels
    integer :: status                   ! allocation status 
    
    integer :: quantity_type, L_quantity_type       ! added on Jul 13

    integer, dimension(:), pointer :: closestInstances 

    real(r8), dimension(:,:), pointer :: A_CLEARSKYRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDINDUCEDRADIANCE
    real(r8), dimension(:,:), pointer :: A_CLOUDEXTINCTION
    real(r8), dimension(:,:), pointer :: A_CLOUDRADSENSITIVITY
    real(r8), dimension(:,:), pointer :: A_EFFECTIVEOPTICALDEPTH
    real(r8), dimension(:,:), pointer :: A_MASSMEANDIAMETER
    real(r8), dimension(:,:), pointer :: A_TOTALEXTINCTION

    real(r8), dimension(:), pointer :: FREQUENCIES
    real(r8), dimension(:,:), allocatable :: WC

    logical :: Got(2) = .false.  

    nullify( CLOUDICE, CLOUDWATER, CLOUDEXTINCTION, CLOUDINDUCEDRADIANCE,    &
             CLOUDRADSENSITIVITY, EFFECTIVEOPTICALDEPTH, GPH,                &
             MASSMEANDIAMETERICE, MASSMEANDIAMETERWATER, PTAN,               &
             RADIANCE, SIZEDISTRIBUTION, EARTHRADIUS, SURFACETYPE,         &
             TEMP, TOTALEXTINCTION, VMR, VMRARRAY,closestInstances,          &
             A_CLEARSKYRADIANCE, A_CLOUDINDUCEDRADIANCE,                     &
             A_CLOUDEXTINCTION, A_CLOUDRADSENSITIVITY,                       &
             A_EFFECTIVEOPTICALDEPTH, A_MASSMEANDIAMETER,                    &
             A_TOTALEXTINCTION, FREQUENCIES )

    ! Check the model configuration 
    if ( size ( forwardModelConfig%signals ) /= 1 )                          &
      & call MLSMessage ( MLSMSG_Error, ModuleName,                          &
      & 'Cannot call the full cloud forward model with multiple signals' )
    signal = forwardModelConfig%signals(1)
    maf = fmStat%maf

    ! For the moment make it only single sideband
    if ( signal%sideband == 0 ) call MLSMessage ( MLSMSG_Error, ModuleName,  &
      & 'Only single sidebands allowed in FullForwardCloudModel for now' )
    call Allocate_test ( frequencies, count ( signal%channels ),             &
      & 'frequencies', ModuleName )
    frequencies = signal%lo + signal%sideband * ( signal%centerFrequency +   &
      & pack ( signal%frequencies, signal%channels ) )
    noFreqs = size (frequencies)

    ! Get the quantities we need from the vectors

    ! Outputs
    do quantity_type = 1, fwdModelOut%template%noQuantities
      l_quantity_type = fwdModelOut%quantities(quantity_type)%template%quantityType
      print*,'quantity_type: ', 'outputs', quantity_type
      print*,'l_quantity_type: ', 'outputs', l_quantity_type
      select case (l_quantity_type)
        case (l_radiance) 
          radiance => GetVectorQuantityByType ( fwdModelOut,                 &
          & quantityType=l_radiance,                                         &
          & signal=signal%index, sideband=signal%sideband )
        case (l_cloudInducedRadiance)
          cloudInducedRadiance => GetVectorQuantityByType ( fwdModelOut,     &
          & quantityType=l_cloudInducedRadiance,                             &
          & signal=signal%index, sideband=signal%sideband )
        case (l_cloudExtinction)
          cloudExtinction => GetVectorQuantityByType ( fwdModelOut,          &
          & quantityType=l_cloudExtinction,                                  &
          & signal=signal%index, sideband=signal%sideband )
        case (l_cloudRADSensitivity)
          cloudRADSensitivity => GetVectorQuantityByType ( fwdModelOut,      &
          & quantityType=l_cloudRADSensitivity,                              &
          & signal=signal%index, sideband=signal%sideband )
        case (l_totalExtinction)
          totalExtinction => GetVectorQuantityByType ( fwdModelOut,          &
          & quantityType=l_totalExtinction,                                  &
          & signal=signal%index, sideband=signal%sideband )
        case (l_effectiveOpticalDepth)
          effectiveOpticalDepth => GetVectorQuantityByType ( fwdModelOut,    &
          & quantityType=l_effectiveOpticalDepth,                            &
          & signal=signal%index, sideband=signal%sideband )
        case (l_massMeanDiameterIce)
          massMeanDiameterIce => GetVectorQuantityByType ( fwdModelOut,      &
          & quantityType=l_massMeanDiameterIce )
        case (l_massMeanDiameterWater)
          massMeanDiameterWater => GetVectorQuantityByType ( fwdModelOut,    &
          & quantityType=l_massMeanDiameterWater )
        case default
          print*, 'l_radiance: ', l_radiance
          print*, 'l_cloudInducedRadiance: ', l_cloudInducedRadiance
          print*, 'l_cloudextinction: ', l_cloudextinction
          print*, 'l_massmeandiameterice: ', l_massmeandiameterice
          print*, 'l_cloudRADSensitivity: ', l_cloudRADSensitivity
          print*, 'l_totalExtinction: ', l_totalExtinction
          print*, 'l_effectiveOpticalDepth: ', l_effectiveOpticalDepth
          print*, 'l_massMeanDiameterWater: ', l_massMeanDiameterWater
          print*, 'l_quantity_type: ', l_quantity_type

          call MLSMessage ( MLSMSG_Error, ModuleName,                        &
                            'Did not understand output l_quantity_types')
      end select
    enddo

!-------------------------------------------------------------------

    ! Inputs
    NQ1 = fwdModelIn%template%noQuantities
    NQ2 = fwdModelExtra%template%noQuantities
    do quantity_type = 1, NQ1+NQ2
      if (quantity_type .le. NQ1) &
      l_quantity_type = fwdModelIn%quantities(quantity_type)%template%quantityType
      if (quantity_type .gt. NQ1) &
      l_quantity_type = fwdModelExtra%quantities(quantity_type-NQ1)%template%quantityType
      print*,'quantity_type: ', 'inputs', quantity_type
      print*,'l_quantity_type: ', 'inputs', l_quantity_type

      select case (l_quantity_type)
        case (l_ptan)
          ptan => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,      &
          & quantityType=l_ptan, instrumentModule = radiance%template%instrumentModule)
        case (l_temperature)
          temp => GetVectorQuantityByType ( fwdModelIn,  fwdModelExtra,      &
          & quantityType=l_temperature )
        case (l_gph)
          gph => GetVectorQuantityByType ( fwdModelIn,  fwdModelExtra,       &
          & quantityType=l_gph )
        case (l_cloudIce)
          cloudIce => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,   &
          & quantityType=l_cloudIce )
        case (l_cloudWater)
          cloudWater => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_cloudWater )
        case (l_surfaceType)
          surfaceType => GetVectorQuantityByType ( fwdModelIn,fwdModelExtra, &
          & quantityType=l_surfaceType )
        case (l_sizeDistribution)
          sizeDistribution=>GetVectorQuantityByType(fwdModelIn,fwdModelExtra, &
          & quantityType=l_sizeDistribution )
        case (l_earthradius)
          earthradius=>GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_earthradius ) 
        case (l_vmr)
!          need to do nothing, will be treated below.
        case default
          call MLSMessage ( MLSMSG_Error, ModuleName,                        &
                            'Did not understand Input l_quantity_types')
      end select
    enddo

    ! Make sure the quantities we have are OK
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,    &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(gph, stacked=.true., coherent=.true.,  &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,    &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true.,               &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,    &
      & ModuleName, InvalidQuantity//'ptan' )

    ! Set up some temporary quantities
    call Allocate_test ( closestInstances, radiance%template%noInstances,    &
      & 'closestInstances', ModuleName )      

    ! Assemble the vmr array

    if ( size(forwardModelConfig%molecules) .lt. 2 ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, 'Not enough molecules' )
    endif
    call allocate_test ( vmrArray,                                           &
      & size(forwardModelConfig%molecules), temp%template%noSurfs,           &
      & 'vmrArray', ModuleName )

    ivmr=0
    do i = 1, size(forwardModelConfig%molecules)
      select case (forwardModelConfig%molecules(i))
        case(L_H2O)
          ivmr=1
        case(L_O3)
          ivmr=2
        case default
          ivmr=0
      end select
      if(ivmr==0) then
        cycle
      endif
      print*, 'i: ', i, 'about to get vmr for molecule of i'
      vmr => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra,            &
        & quantityType=l_vmr, molecule=forwardModelConfig%molecules(i) )
      if (.not.ValidateVectorQuantity( vmr, stacked=.true., coherent=.true., &
        & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error,  &
        & ModuleName, InvalidQuantity//'vmr' )
      print*, 'i: ', i, 'got vmr for molecule of i'

      call FindClosestInstances ( vmr, radiance, closestInstances )
      vmrInst = closestInstances(maf)

      call InterpolateValues ( &
        & vmr%template%surfs(:,1), &    ! Old X
        & vmr%values(:,vmrInst), &      ! Old Y
        & temp%template%surfs(:,1), &   ! New X
        & vmrArray(i,:), &              ! New Y
        & 'Linear', extrapolate='Clamp' )
      Got(ivmr)=.true.
    end do

    if ( .not. got(1) .or. .not. got(2) ) then
      call MLSMessage( MLSMSG_Error, ModuleName,                             &
                      'Missing the required molecules' )
    endif

    ! Work out the closest instances for the other quantities
    call FindClosestInstances ( temp, radiance, closestInstances )
    instance = closestInstances(maf)
    noLayers = temp%template%noSurfs - 1

    ! Make temporary arrays for the cloud forward model
    call Allocate_test ( a_clearSkyRadiance,                                 &
      & radiance%template%noSurfs, noFreqs,                                  &
      & 'a_clearSkyRadiance', ModuleName )
    call Allocate_test ( a_cloudInducedRadiance,                             &
      & radiance%template%noSurfs, noFreqs,                                  &
      & 'a_cloudInducedRadiance', ModuleName )
    call Allocate_test ( a_effectiveOpticalDepth,                            &
      & radiance%template%noSurfs, noFreqs,                                  &
      & 'a_effectiveOpticalDepth', ModuleName )
    call Allocate_test ( a_cloudRADSensitivity,                              &
      & radiance%template%noSurfs, noFreqs,                                  &
      & 'a_cloudRADSensitivity', ModuleName )
    call Allocate_test ( a_totalExtinction,                                  &
      & temp%template%noSurfs, noFreqs,                                      &
      & 'a_totalExtinction', ModuleName )
    call Allocate_test ( a_cloudExtinction,                                  &
      & temp%template%noSurfs, noFreqs,                                      &
      & 'a_cloudExtinction', ModuleName )
    call Allocate_test ( a_massMeanDiameter,                                 &
      & 2, temp%template%noSurfs,                                            &
      & 'a_massMeanDiameter', ModuleName )
    
    ! Now call the full CloudForwardModel code

    NOsurf=temp%template%noSurfs
    if (Nosurf /= GPH%template%nosurfs) then
      call MLSMessage ( MLSMSG_Error, ModuleName,                            &
      & 'number of levels in gph does not match no of levels in temp' )
     else if (radiance%template%nosurfs /= ptan%template%nosurfs) then
        call MLSMessage ( MLSMSG_Error, ModuleName,                          &
        & 'number of levels in radiance does not match no of levels in ptan' )
    endif

    allocate ( WC(2,NOsurf), STAT=status )

    WC (1,:) = CloudIce%values(:,instance)
    WC (2,:) = CloudWater%values(:,instance)

    call CloudForwardModel (                                                 &
      & noFreqs,                                                             &
      & noSurf,                                                              & 
      & radiance%template%noSurfs,                                           &
      & size( ForwardModelConfig%molecules),                                 &
      & ForwardModelConfig%no_cloud_species,                                 &
      & ForwardModelConfig%no_model_surfs,                                   &
      & frequencies/1e3_r8,                                                  &
      & 10.0**(-temp%template%surfs),                                        &
      & gph%values(:, instance),                                             &
      & temp%values(:,instance),                                             &
      & vmrArray,                                                            &
      & WC,                                                                  &
      & sizeDistribution%values(1,instance),                                 &
      & 10.0**(-ptan%values(:,maf)),                                         &
      & earthradius%values(1,1),                                           &
      & surfaceType%values(1, instance),                                     &
      & forwardModelConfig%cloud_der,                                        &
      & forwardModelConfig%cloud_width,                                      &
      & a_clearSkyRadiance,                                                  &
      & a_cloudInducedRadiance,                                              &
      & a_totalExtinction,                                                   &
      & a_cloudExtinction,                                                   &
      & a_massMeanDiameter,                                                  &
      & a_effectiveOpticalDepth,                                             &
      & cloudRADSensitivity,                                                 &
      & forwardModelConfig%NUM_SCATTERING_ANGLES,                            &  
      & forwardModelConfig%NUM_AZIMUTH_ANGLES,                               &
      & forwardModelConfig%NUM_AB_TERMS,                                     &
      & forwardModelConfig%NUM_SIZE_BINS )

    deallocate (WC, stat=status)

    ! Now store results in relevant vectors
    ! Vectors are stored (noChannels*noSurfaces, noInstances), so transpose
    ! and reshape temporary variables to be in the right form.
    ! First the minor frame stuff
    radiance%values ( :, maf) =                                              &
      & reshape ( transpose(a_clearSkyRadiance),                             &
      & (/radiance%template%instanceLen/) )
    cloudInducedRadiance%values ( :, maf ) =                                 &
      & reshape ( transpose(a_cloudInducedRadiance),                         &
      & (/cloudInducedRadiance%template%instanceLen/) )
    effectiveOpticalDepth%values ( :, maf ) =                                &
      & reshape ( transpose(a_effectiveOpticalDepth),                        &
      & (/effectiveOpticalDepth%template%instanceLen/) )
    cloudRADSensitivity%values ( :, maf ) =                                  &
      & reshape ( transpose(a_cloudRADSensitivity),                          &
      & (/cloudRADSensitivity%template%instanceLen/) )

 ! For layer(noTempSurfs-1) stuff make sure all are zero to start, then do rest
    cloudExtinction%values(:,instance) =       0.0_r8
    massMeanDiameterIce%values(:,instance) =   0.0_r8
    massMeanDiameterWater%values(:,instance) = 0.0_r8
    totalExtinction%values(:,instance) =       0.0_r8
    
    cloudExtinction%values ( 1:noLayers, instance ) =                        &
      & reshape ( transpose(a_cloudExtinction), (/noLayers*noFreqs/) )
    massMeanDiameterIce%values (1:noLayers,instance)=a_massMeanDiameter(1,:)
    massMeanDiameterWater%values(1:noLayers,instance)=a_massMeanDiameter(2,:)
    totalExtinction%values ( 1:noLayers, instance ) =                        &
      & reshape ( transpose(a_totalExtinction), (/noLayers*noFreqs/) )

    ! Remove temporary quantities
    call Deallocate_test ( a_massMeanDiameter,                               &
                          'a_massMeanDiameterModuleName', ModuleName )
    call Deallocate_test ( a_cloudExtinction,                                &
                          'a_cloudExtinction',            ModuleName )
    call Deallocate_test ( a_totalExtinction,                                &
                          'a_totalExtinction',            ModuleName )
    call Deallocate_test ( a_cloudRADSensitivity,                            &
                          'a_cloudRADSensitivity',        ModuleName )
    call Deallocate_test ( a_effectiveOpticalDepth,                          &
                          'a_effectiveOpticalDepth',      ModuleName )
    call Deallocate_test ( a_cloudInducedRadiance,                           &
                          'a_cloudInducedRadiance',       ModuleName )
    call Deallocate_test ( a_clearSkyRadiance,                               &
                          'a_clearSkyRadiance',           ModuleName )
    call Deallocate_test ( vmrArray,                                         &
                          'vmrArray',                     ModuleName )
    call Deallocate_test ( closestInstances,                                 &
                          'closestInstances',             ModuleName )

  end subroutine FullCloudForwardModelWrapper

end module FullCloudForwardModel

