! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PolarLinearModel_m

  ! This is a special instance of the linearized forward model which is
  ! designed to handle polarized case.  It actually calls the
  ! regular linearized forward model three times as part of its computation

  implicit none
  private
  public :: PolarLinearModel

  !------------------------------- RCS Ident Info ------------------------------
  character(len=*), parameter :: IdParm = & 
    "$Id$"
  character(len=len(idParm)) :: Id = idParm
  character(len=*), parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here 
  !-----------------------------------------------------------------------------


contains ! =====     Public Procedures     =============================

  subroutine PolarLinearModel ( fmConf, FwdModelIn, FwdModelExtra,&
    & FwdModelOut, Ifm, fmStat, Jacobian, vectors )

    ! Import stuff
    use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use VectorsModule, only: VECTOR_T, VECTORVALUE_T, ADDTOVECTOR, &
      & CLEARVECTOR, DESTROYVECTORINFO, GETVECTORQUANTITYBYTYPE, CLONEVECTOR
    use ForwardModelConfig, only: FORWARDMODELCONFIG_T
    use ForwardModelIntermediate, only: FORWARDMODELSTATUS_T, &
      & FORWARDMODELINTERMEDIATE_T
    use Intrinsic, only: L_LINEAR, L_FULL, L_RADIANCE, L_FIELDAZIMUTH
    use FullForwardModel_m, only: FULLFORWARDMODEL
    use L2PC_m, only: BINSELECTORS, DEFAULTSELECTOR_FIELDAZIMUTH
    use LinearizedForwardModel_m, only: LINEARIZEDFORWARDMODEL
    use MatrixModule_1, only: MATRIX_T, COPYMATRIX, ADDTOMATRIX, CLEARMATRIX, &
      & CREATEEMPTYMATRIX, DESTROYMATRIX
    use MLSCommon, only: R8
    use MLSSignals_M, only: SIGNAL_T
    use Units, only: DEG2RAD
    use ForwardModelVectorTools, only: GETQUANTITYFORFORWARDMODEL
    use ManipulateVectorQuantities, only: FINDONECLOSESTINSTANCE

    ! Dummy arguments
    type(forwardModelConfig_T), intent(inout) :: FMCONF
    type(vector_T), intent(in) ::  FWDMODELIN
    type(vector_T), intent(in) ::  FWDMODELEXTRA
    type(vector_T), intent(inout) :: FWDMODELOUT  ! Radiances, etc.
    type(forwardModelIntermediate_T), intent(inout) :: IFM ! Workspace
    type(forwardModelStatus_t), intent(inout) :: FMSTAT ! Reverse comm. stuff
    type(matrix_T), intent(inout), optional :: JACOBIAN
    type(vector_T), dimension(:), pointer, optional :: VECTORS

    ! Local parameters
    ! Local variables
    type(forwardModelConfig_T) :: THISCONFIG
    type(Vector_T), target :: RADIANCECONTRIBUTION
    type(Matrix_T), target :: JACOBIANCONTRIBUTION
    type(VectorValue_T), pointer :: PHI
    type(VectorValue_T), pointer :: RADIANCE
    type(VectorValue_T), pointer :: FIELDAZIMUTH

    integer :: PASS                     ! Which 'pass' of the linear model are we on
    integer :: MAF                      ! Which MAF are we computing for
    integer :: INSTANCE                 ! Instance of field azimuth
    real(r8) :: SIN2PHI, COS2PHI        ! What they say!
    real(r8) :: SCALING                 ! How much to multiply this contribution by
    real(r8), dimension(:,:), pointer :: SAVEPHIVALUES
    ! A store for the original values of phi

    ! Setup our temporary stuff
    if ( present(jacobian) ) then
      call CreateEmptyMatrix ( jacobianContribution, 0, &
        & jacobian%row%vec, jacobian%col%vec, &
        & .not. jacobian%row%instFirst, .not. jacobian%col%instFirst, &
        & 'jacobianContribution' )
    end if
    call CloneVector ( radianceContribution, fwdModelOut )

    maf = fmStat%maf

    ! This forward model configuration is somewhat crafty.   We create our
    ! own configuration we pass to the linear model 3 times to form
    ! our composite radiances.  To do this we patch in the default bin
    ! selector for the field azimuth and then hack around at it's value.

    thisConfig = fmConf
    nullify ( thisConfig%binSelectors )
    call Allocate_test ( thisConfig%binSelectors, size(fmConf%binSelectors)+1, &
      & 'thisConfig%binSelectors', ModuleName )
    
    thisConfig%binSelectors(1) = DefaultSelector_FieldAzimuth
    thisConfig%binSelectors(2:) = fmConf%binSelectors

    ! Get stuff from various vectors
    radiance => GetVectorQuantityByType (fwdModelOut, quantityType=l_radiance, &
      & signal=fmConf%signals(1)%index, sideband=fmConf%signals(1)%sideband )
    phi => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
      & quantityType=l_fieldAzimuth, config=fmConf )
    instance = FindOneClosestInstance ( phi, radiance, maf )

    ! Create terms for scalings
    sin2phi = sin ( 2*phi%values(1,instance) * deg2rad )
    cos2phi = cos ( 2*phi%values(1,instance) * deg2rad )

    ! Store the values of phi away for the moment
    nullify ( savePhiValues )
    call Allocate_Test ( savePhiValues, phi%template%instanceLen, &
      & phi%template%noInstances, 'savePhiValues', ModuleName )
    savePhiValues = phi%values

    ! Now loop over the three calls we're going to make
    do pass = 1, 3
      ! Get the appropriate values of phi and scalings
      select case ( pass )
      case ( 1 ) ! phi=0 case
        scaling = ( 1.0_r8 - sin2phi + cos2phi ) / 2.0_r8
        phi%values = 0.0_r8
      case ( 2 ) ! phi=pi/4 case
        scaling = sin2phi
        phi%values = 45.0_r8
      case ( 3 ) ! phi=pi/2 case
        scaling = ( 1.0_r8 - sin2phi - cos2phi ) / 2.0_r8
        phi%values = 90.0_r8
      end select

      ! Do this one's contribution
      call ClearVector ( radianceContribution )
      if ( present ( jacobian ) ) then
        call ClearMatrix ( jacobianContribution )
        call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
          & radianceContribution, ifm, fmStat, jacobianContribution, vectors )
      else
        call LinearizedForwardModel ( thisConfig, fwdModelIn, fwdModelExtra, &
          & radianceContribution, ifm, fmStat, vectors=vectors )
      end if
      call AddToVector ( fwdModelOut, radianceContribution, scaling )
      if ( present ( jacobian ) ) &
        & call AddToMatrix ( jacobian, jacobianContribution, scaling )
    end do

    ! Put the original phi back
    phi%values = savePhiValues

    ! Tidy up
    call Deallocate_test ( savePhiValues, 'savePhiValues', ModuleName )
    if ( present(jacobian) )call DestroyMatrix ( jacobianContribution )
    call DestroyVectorInfo ( radianceContribution )
    call Deallocate_test ( thisConfig%binSelectors, 'thisConfig%binSelectors', ModuleName )

  end subroutine PolarLinearModel

  ! ----------------------------------------------------------------------------

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PolarLinearModel_m
