! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module SwitchingMirrorModel_m

  ! This model contains a special forward model for dealing with
  ! the 'radiance at the switching mirror' stuff

  implicit none
  private
  public :: SwitchingMirrorModel

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ====================================== SwitchingMirrorModel =====

  subroutine SwitchingMirrorModel ( FwdModelConf, FwdModelIn, FwdModelExtra, &
    & FwdModelOut, oldIFM, fmStat, jacobian )
    ! Imports
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use D_STAT_TEMP_M, only: STAT_TEMP
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use ForwardModelIntermediate, only: FORWARDMODELINTERMEDIATE_T, FORWARDMODELSTATUS_T
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, GETVECTORQUANTITYBYTYPE
    use MatrixModule_1, only: MATRIX_T
    use Intrinsic, only: L_RADIANCE, L_REFLREFL, L_REFLSPILL, L_REFLTEMP, &
      & L_REFLTRANS, L_STRAYRADIANCE, L_PRIMARY, L_SECONDARY, L_TERTIARY, L_COMPLETE, &
      & L_CALSIDEBANDFRACTION, L_THZ, LIT_INDICES
    use MLSSignals_m, only: SIGNALS, SIGNAL_T, GETSIDEBANDLOOP, MODULES
    use MLSCommon, only: R8

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: fwdModelConf
    type(vector_T), intent(in) ::  FwdModelIn, FwdModelExtra
    type(vector_T), intent(inout) :: FwdModelOut  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: oldIfm ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FmStat ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: Jacobian

    ! Local parameters
    integer, parameter :: NOREFLECTORS = 3
    integer, dimension(noReflectors), parameter :: REFLECTORS = &
      & (/ l_primary, l_secondary, l_tertiary /)

    ! Local variables
    integer :: MAF                      ! Which MAF are we doing
    integer :: SIGINDEX                 ! Which signal in our list are we doing
    integer :: CHAN                     ! Which channel are we doing
    integer :: NOCHANS                  ! Number of channels for this signal
    integer :: SIDEBAND                 ! Which sideband are we considering
    integer :: REFLECTOR                ! Which reflector are we doing
    integer :: SIDEBANDSTART            ! Loop start
    integer :: SIDEBANDSTOP             ! Loop end
    integer :: SIDEBANDSTEP             ! Loop step
    integer :: THZMODULE                ! Index of THz module

    real(r8), dimension(noReflectors) :: ACLPRODUCT ! An intermediate term in the calculation
    real(r8), dimension(:), pointer :: OFFSET ! The offset to add for each channel
    real(r8), dimension(:), pointer :: EMISSION ! Thermal emission term
    real(r8), dimension(:,:), pointer :: DELTATRANS ! Difference in transmissions

    type(Signal_T), pointer :: SIGNAL   ! One signal to deal with
    type(VectorValue_T), pointer :: RADIANCE    ! One radiance to fill
    type(VectorValue_T), pointer :: REFLREFL    ! Reflectivity for one reflector
    type(VectorValue_T), pointer :: REFLSPILL   ! Spill over for one reflector
    type(VectorValue_T), pointer :: REFLTEMP    ! Temperature for one reflector
    type(VectorValue_T), pointer :: REFLTRANS   ! Transmission for one reflector
    type(VectorValue_T), pointer :: STRAYRADIANCE ! Stray radiance term
    type(VectorValue_T), pointer :: TOTALTRANS  ! Transmission of the complete system
    type(VectorValue_T), pointer :: PRIMARYTRANS  ! Transmission of the primary
    type(VectorValue_T), pointer :: SIDEBANDFRACTION  ! Transmission of the primary

    ! Executable code
    ! Don't do this model if we've not been told to
    if ( .not. fwdModelConf%switchingMirror ) return
    ! Setup various low level stuff
    nullify ( offset, deltaTrans, emission )
    maf = fmStat%maf

    ! Setup an intermediate result
    do reflector = noReflectors, 1, -1
      reflRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
        & quantityType=l_reflRefl, reflector=reflectors(reflector) )
      aclProduct ( reflector ) = reflRefl%values(1,1)
      if ( reflector < noReflectors ) aclProduct ( reflector ) = &
        & aclProduct ( reflector ) * aclProduct ( reflector+1 )
    end do

    ! Loop over the signals this forward model is considering
    signalLoop: do sigIndex = 1, size(fwdModelConf%signals)
      ! What signal is this?
      signal => fwdModelConf%signals(sigIndex)
      ! Don't do this model at all for the THz
      if ( modules(signal%instrumentModule)%name == lit_indices(l_THz) ) cycle signalLoop

      ! Look for the radiance quantity for this signal
      radiance => GetVectorQuantityByType ( fwdModelOut, quantityType=l_radiance, &
        & signal=signal%index, sideband=signal%sideband )
      ! Setup some work arrays
      noChans = size ( signal%frequencies )
      call Allocate_Test ( emission, noChans, 'emission', ModuleName )
      call Allocate_Test ( offset, noChans, 'offset', ModuleName )
      call Allocate_Test ( deltaTrans, noChans, noReflectors, 'deltaTrans', ModuleName )
      
      ! Loop over the sidebands
      call GetSidebandLoop ( signal%index, signal%sideband, .true., &
        & sidebandStart, sidebandStop, sidebandStep )
      do sideband = sidebandStart, sidebandStop, sidebandStep
        ! Setup another intermediate quantity
        do reflector = 1, noReflectors
          reflTrans => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_reflTrans, reflector=reflectors(reflector), &
            & signal=signal%index, sideband=sideband )
          deltaTrans ( :, reflector ) = reflTrans%values(:,1)
        end do
        do reflector = 1, noReflectors - 1
          deltaTrans ( :, reflector ) = deltaTrans ( :, reflector+1 ) - &
            & deltaTrans ( :, reflector )
        end do
        deltaTrans ( :, noReflectors ) = 1.0 - deltaTrans ( :, noReflectors )
        
        ! Now get some state vector quantities for this signal
        totalTrans => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_reflTrans, reflector=l_complete, &
          & signal=signal%index, sideband=sideband )
        strayRadiance => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_strayRadiance, &
          & signal=signal%index, sideband=sideband )
        primaryTrans => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_reflTrans, reflector=l_primary, &
          & signal=signal%index, sideband=sideband )
        sidebandFraction => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
          & quantityType=l_calSidebandFraction, &
          & signal=signal%index, sideband=sideband )

        ! Do the reflector independent stuff first
        offset = aclProduct ( 1 ) * primaryTrans%values(:,1) * &
          & ( 1.0 - totalTrans%values(:,1) ) * strayRadiance%values(:,maf)
        
        ! Now loop over the reflectors
        do reflector = 1, noReflectors
          ! Get more relevant quantities
          reflRefl => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_reflRefl, reflector=reflectors(reflector) )
          reflTemp => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_reflTemp, reflector=reflectors(reflector) )
          reflSpill => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_reflSpill, reflector=reflectors(reflector), &
            & signal=signal%index, sideband=sideband )
          reflTrans => GetVectorQuantityByType ( fwdModelIn, fwdModelExtra, &
            & quantityType=l_reflTrans, reflector=reflectors(reflector), &
            & signal=signal%index, sideband=sideband )
          ! Do the actual calculation
          emission = stat_temp ( reflTemp%values(1,1), &
            & signal%lo + sideband * &
            & ( signal%direction*signal%frequencies+signal%centerFrequency ) )
          offset = offset + aclProduct ( reflector ) * &
            & ( ( 1 - reflrefl%values(1,1) ) * reflTrans%values(:,1) * emission + &
            &   deltaTrans(:,reflector) * reflSpill%values(:,maf) )
        end do                          ! End loop over reflectors
        ! Finally multiply by the sideband fraction
        offset = offset * sidebandFraction%values(:,1)
        ! Now adjust the appropriate channels
        channelLoop: do chan = 1, noChans
          ! Perhaps skip this channel
          if ( associated ( signal%channels ) ) then
            if ( .not. signal%channels ( chan ) ) cycle channelLoop
          end if
          radiance%values ( chan:radiance%template%instanceLen:noChans, &
            & maf ) = &
            & radiance%values ( chan:radiance%template%instanceLen:noChans, &
            & maf ) + offset ( chan )
        end do channelLoop            ! End loop over channels
      end do                            ! End loop over sidebands
      call Deallocate_test ( emission, 'emission', ModuleName )
      call Deallocate_test ( deltaTrans, 'deltaTrans', ModuleName )
      call Deallocate_test ( offset, 'offset', ModuleName )
    end do signalLoop                   ! End loop over signals
  end subroutine SwitchingMirrorModel

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here
end module SwitchingMirrorModel_m

! $Log$
! Revision 2.1  2003/05/29 16:46:06  livesey
! First version
!
