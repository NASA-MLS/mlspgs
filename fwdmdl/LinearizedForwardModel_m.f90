! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module LinearizedForwardModel_m

  ! Approximate a forward model with a first-order Taylor's series.

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use Intrinsic, only: LIT_INDICES
  use L2PC_m, only: L2PCDATABASE
  use ForwardModelConfig, only: FORWARDMODELCONFIG_T
  use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T, &
    & FORWARDMODELINTERMEDIATE_T
  use ManipulateVectorQuantities, only: FINDCLOSESTINSTANCES
  use MatrixModule_1, only:  MATRIX_T, MULTIPLYMATRIXVECTORNOT, DUMP
  use MLSCommon, only: r8
  use MLSSignals_m, only: Signal_T
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ERROR
  use MLSNumerics, only: INTERPOLATEVALUES
  use VectorsModule, only: assignment(=), OPERATOR(-), OPERATOR(+), &
    & CLONEVECTOR, CONSTRUCTVECTORTEMPLATE, COPYVECTOR, CREATEVECTOR,&
    & DESTROYVECTORINFO, GETVECTORQUANTITYBYTYPE, VECTOR_T, &
    & VECTORVALUE_T, VECTORTEMPLATE_T, DUMP
  use QuantityTemplates, only: QuantityTemplate_T
  use String_Table, only: Display_String
  use Intrinsic, only: L_RADIANCE, L_TEMPERATURE, L_PTAN,L_VMR

  implicit none
  private
  public :: LinearizedForwardModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
       "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  !-----------------------------------------------------------------------------

contains ! =====     Public Procedures     =============================

  ! -------------------------------------  LinearizedForwardModel  -----
  subroutine LinearizedForwardModel ( ForwardModelConfig, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, Ifm, fmStat, Jacobian )
    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FORWARDMODELCONFIG
    type(vector_T), intent(in) ::  FWDMODELIN, FwdModelExtra
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN

    ! Local variables
    integer :: MAF                      ! Major frame to do
    integer :: CENTERPROF               ! Center profile corresponding to MAF
    integer :: CLOSESTINSTANCE          ! A closest profile to this MAF
    integer :: INSTANCE                 ! Loop index
    integer :: L2PCBIN                  ! Which l2pc to use
    integer :: QTYIND                   ! Loop index for main loop
    integer :: XSTARINSTANCE            ! Loop counter
    integer :: XINSTANCE                ! Instance in x corresponding to xStarInstance

    logical :: DODERIVATIVES            ! Flag
    logical :: FOUNDINFIRST             ! Flag for state quantities

    integer, dimension(:), pointer :: closestInstances ! Of qty to radiance

    ! The `prime' quantities are important.
    ! - yPrime (yP) is contains one maf of the relevant radiances, but on the
    !   l2pc's vertical coordinates
    ! - xPrime (xP) is like xStar but contains state values
    ! - kPrime is the jacobian for these two.

    real (r8), dimension(:,:), pointer :: yPmapped ! Remapped values of yP
    real (r8), dimension(:,:), pointer :: resultMapped ! Remapped values of result
    real (r8), dimension(:,:), pointer :: dyByDX ! Raw dRad/dPtan

    type(vector_T) :: XP                ! Same form as xStar, contents as x
    type(vector_T) :: YP                ! Same form as yStar,=kstar*(xp-xStar)
    type(vector_T) :: DELTAX            ! xp-xStar
    type(Signal_T) :: signal            ! A signal

    type(Matrix_T), pointer :: L2PC       ! The l2pc to use
    type(VectorValue_T), pointer :: RADIANCE ! The radiance quantity to fill
    type(VectorValue_T), pointer :: RADINL2PC ! The relevant radiance part of yStar
    type(VectorValue_T), pointer :: STATEQ ! A state vector quantity
    type(VectorValue_T), pointer :: PTAN ! Tangent pressure quantity
    type(VectorValue_T), pointer :: XSTARPTAN ! Tangent pressure in l2pc
    type(VectorValue_T), pointer :: L2PCQ ! A quantity in the l2pc
    ! Executable code

    nullify ( yPmapped, resultMapped, dyByDX, closestInstances )

    ! Identify the band/maf we're looking for

    maf = fmStat%maf
    print*,'Linear model doing maf:',maf
    if ( size ( forwardModelConfig%signals ) /= 1 ) call MLSMessage ( &
      & MLSMSG_Error, ModuleName, &
      & 'Can only have one signal for linearized models')

    signal = forwardModelConfig%signals(1)
    if ( signal%sideband /= 0 ) call MLSMessage(MLSMSG_Error, ModuleName,&
      & "Linearized forward model can only do folded signals for now")
    radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal= signal%index, sideband=signal%sideband )

    ! Here in later versions, we'd identify which l2pc bin to use
    ! It has to go here as we need to know yStar%noSurfs for the next
    ! part.
    ! For the moment choose the first
    l2pcBin = 1
    l2pc => l2pcDatabase(l2pcBin)

    radInl2pc => GetVectorQuantityByType (l2pc%row%vec, quantityType=l_radiance, &
      & signal=signal%index, sideband=signal%sideband )

    ! Now we loop over the quantities in the l2pc file and construct an xPrime
    ! for them
    call cloneVector ( xP, l2pc%col%vec )
    call cloneVector ( deltaX, xP )     ! Note sets values to 0.0

    ! Set up some other stuff before main loop
    call Allocate_test ( closestInstances, radiance%template%noInstances, &
      & 'closestInstances', ModuleName )

    ! -------- Main loop over xStar quantities -------------------------------
    do qtyInd = 1, size ( l2pc%col%vec%quantities )

      ! Identify this quantity in xStar
      l2pcQ => l2pc%col%vec%quantities(qtyInd)
      
      ! Now see if we are wanting to deal with this
      if ( l2pcQ%template%quantityType == l_ptan ) cycle ! Get this from interpolation
      if ( l2pcQ%template%quantityType == l_vmr ) then
        if ( .not. any &
          & (l2pcQ%template%molecule == forwardModelConfig%molecules)) cycle
      end if

      ! Identify this quantity in x
      select case ( l2pcQ%template%quantityType )
      case ( l_temperature )
        stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
          & quantityType = l_temperature, &
          & foundInFirst = foundInFirst, noError=.true. )
      case ( l_vmr )
        stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
          & quantityType = l_vmr, &
          & molecule = l2pcQ%template%molecule, &
          & foundInFirst = foundInFirst, noError=.true. )
      case default
        stateQ => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra,&
          & quantityType = l2pcQ%template%quantityType, &
          & foundInFirst = foundInFirst, noError=.true. )
      end select

      ! If it's not in the state vector, perhaps make a fuss
      if ( .not. associated(stateQ) ) then
        if ( any(l2pcQ%template%quantityType == &
          & (/ l_temperature, l_vmr /))) &
          &   call MLSMessage ( MLSMSG_Error, ModuleName, &
          &  "Temperature or vmr quantity absent from state")
      end if

      ! Now check that surfaces are the same
      if (stateQ%template%noSurfs /= l2pcQ%template%noSurfs) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "State vector surfaces not same as l2pc" )
      ! The l2pc writer insisted that things were on zeta coords or none
      ! at all, so check surfaces OK.
      if ( any (abs (stateQ%template%surfs-l2pcQ%template%surfs) > 0.01)) &
        & call MLSMessage(MLSMSG_Error,ModuleName,&
        & 'State vector surfaces not same as l2pc')

      ! Check that no. chans. is the same (will it ever /=1!!?)
      if (stateQ%template%noChans /= l2pcQ%template%noChans) &
        & call MLSMessage ( MLSMSG_Error, ModuleName, &
        & "State vector channels not same as l2pc" )

      call FindClosestInstances ( stateQ, radiance, closestInstances)
      closestInstance=closestInstances(maf)
      ! Loop over profiles
      do xStarInstance = 1, l2pcQ%template%noInstances
        ! Identify this instance in state
        xInstance = closestInstance - &
          & l2pcQ%template%noInstances/2 + &
          & xStarInstance - 1
        xInstance = max ( 1, min ( xInstance, stateQ%template%noInstances ) )
        ! Fill this part of xP
        xP%quantities(qtyInd)%values(:,xStarInstance) =&
          & stateQ%values(:,xInstance)

        ! Do we need derivatives for this?
        doDerivatives = .true.
        if ( .not. foundInFirst ) doDerivatives = .false.
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_temperature) ) &
          & doDerivatives = forwardModelConfig%temp_der
        if (doDerivatives .and. (l2pcQ%template%quantityType==l_vmr) ) then
          doDerivatives = forwardModelConfig%atmos_der
          if ( doDerivatives .and. .not. any (l2pcQ%template%molecule == &
            & pack(ForwardModelConfig%molecules, &
            &      ForwardModelConfig%moleculeDerivatives))) &
            & doDerivatives = .false.
        end if
        
        ! If so, interpolate this block of kStar and place in jacobian
        if (doDerivatives) then
          !??????????? Code here.
        endif
      end do
      ! End loop over profiles
      ! Compute this instance of deltaX
      deltaX%quantities(qtyInd)%values = &
        & xP%quantities(qtyInd)%values - &
        & l2pc%col%vec%quantities(qtyInd)%values

      !???? SPECIAL TESTS, NJL
      if ( l2pcQ%template%quantityType == l_Temperature ) then
        deltaX%quantities(qtyInd)%values(max(1,l2pcQ%template%noSurfs-maf):,:) = 0.0
      endif

    end do                              ! End loop over quantities

    ! Now compute yP

    call cloneVector( yp, l2pc%row%vec )
    call MultiplyMatrixVectorNoT ( l2pc, deltaX, yP )

    yP = yP + l2pc%row%vec

    ! Now we interpolate yP to ptan

    ptan => GetVectorQuantityByType ( FwdModelIn, FwdModelExtra, &
      & quantityType = l_ptan, &
      & instrumentModule = radiance%template%instrumentModule )

    xStarPtan => GetVectorQuantityByType ( l2pc%col%vec, &
      & quantityType = l_ptan )

    call allocate_test( ypMapped, radInL2PC%template%noSurfs, &
      & radInL2PC%template%noChans, 'ypMapped', ModuleName )
    call allocate_test( resultMapped, radiance%template%noSurfs, &
      & radInL2PC%template%noChans, 'resultMapped', ModuleName )
    call allocate_test( dyByDx, radiance%template%noSurfs, &
      & radInL2PC%template%noChans, 'dyByDX', ModuleName )

    yPmapped = transpose ( &
      & reshape ( yp%quantities(1)%values(:,1), &
      &           (/radInL2PC%template%noChans, radInL2PC%template%noSurfs/) ) )

    call InterpolateValues ( &
      & xStarPtan%values(:,1), &        ! OldX
      & yPmapped, &                     ! OldY
      & ptan%values(:,maf), &           ! NewX
      & resultMapped, &                 ! NewY
      & 'Spline', &                     ! use spline
      & extrapolate='Constant', &       ! dont extrapolate, clamp
      & dyByDx=dyByDx )

    ! For the moment, overwrite, may make it an addition later for unfolded
    ! cases

    radiance%values(:,maf) = reshape(transpose(resultMapped),&
      & (/radiance%template%instanceLen/))

    ! Put ptan derivatives here. !?????????

    call DestroyVectorInfo ( xP )
    call DestroyVectorInfo ( deltaX )
    call Deallocate_test ( closestInstances, 'closestInstances', ModuleName )
    call Deallocate_test ( dyByDx, 'dyByDx', ModuleName)
    call Deallocate_test ( resultMapped, 'resultMapped', ModuleName )
    call Deallocate_test ( ypMapped, 'ypMapped', ModuleName )

    if ( maf == radiance%template%noInstances ) fmStat%finished = .true.

  end subroutine LinearizedForwardModel
end module LinearizedForwardModel_m

! $Log$
! Revision 1.5  2001/04/28 04:40:56  livesey
! Another interim version.
!
! Revision 1.4  2001/04/28 01:54:54  livesey
! Interim, non working but compiling version
!
! Revision 1.3  2001/04/27 17:35:42  livesey
! An interim version for some testing.
!
! Revision 1.2  2001/04/26 23:54:32  livesey
! Interim version
!
! Revision 1.1  2001/04/25 20:49:27  vsnyder
! Initial commit
!
