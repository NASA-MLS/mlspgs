! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ScanModelModule          ! Scan model and associated calculations

  ! This module defines the `scan model' for the MLS Level 2 software.  It also
  ! contains functionality such as a pressure guesser etc, which may be useful
  ! for other programs.

  ! Note that currently the scan model is 1D, later versions will be two
  ! dimensional.  I will try to place comments in the code to indicate where 2D
  ! aware stuff needs to be placed.

  ! Also note that there is currently no height offset term in the state
  ! vector.  This was never used in UMLS V5, and so it is not clear that it
  ! will ever be required.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Constants, only: Deg2rad, Ln10, Pi
  use ForwardModelConfig, only: Forwardmodelconfig_T
  use ForwardModelIntermediate, only: Forwardmodelstatus_T
  use ForwardModelVectortools, only: Getquantityforforwardmodel
  use Geometry, only: EarthRadA, EarthRadB, EarthSurfaceGPH, GeodToGeocLat, &
    & G0, Gm, J2, J4, Omega => W
  use Init_Tables_Module, only: L_Refgph, L_Zeta
  use Intrinsic, only: L_Heightoffset, L_None, L_Ptan, L_Scanresidual, &
    & L_Temperature, L_Tngtgeocalt, L_Vmr, L_Phitan, L_Orbitinclination, &
    & Phyq_Length
  use ManipulateVectorQuantities, only: FindClosestInstances, &
    & Findinstancewindow
  use MatrixModule_0, only: DestroyBlock, MatrixElement_T, M_Absent, &
    & M_Full, UpdateDiagonal
  use MatrixModule_1, only: CreateBlock, FindBlock, Matrix_T, &
    & CreateEmptyMatrix, DestroyMatrix, ClearMatrix, RM
  use MLSKinds, only: R8, RP, RV
  use MLSMessagemodule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
  use MLSNumerics, Only : Hunt, InterpolateValues
  use MLSStringLists, only: SwitchDetail
  use Molecules, only: L_H2o
  use Output_M, only: Output
  use QuantityTemplates, only: Dump
  use Refraction_M, only: MaxRefraction, Refractive_Index, &
    & Refractive_Index_Deriv, Refractive_Index_F
  use Toggles, only: Emit, Toggle, Switches
  use Trace_M, only: Trace_Begin, Trace_End
  use VectorsModule, only: ValidateVectorquantity, &
    & Vector_T, VectorTemplate_T, VectorValue_T, CreateVector, &
    & ConstructVectorTemplate, DestroyVectorinfo

  implicit none

  private

  public :: DestroyForwardModelIntermediate, DumpInstanceWindows, &
    & GetBasisGPH, GetGPHPrecision, &
    & GetHydrostaticTangentPressure, Get2DHydrostaticTangentPressure,      &
    & ScanForwardModel, TwoDScanForwardModel

  type, private :: ForwardModelIntermediate_T

    ! These ones are for the scan model
    real (r8), dimension(:,:),             pointer :: BasisGph=>NULL()
    real (r8), dimension(:),               pointer :: R=>NULL()
    real (r8), dimension(:,:),             pointer :: RT=>NULL()
    integer :: BelowRef                 ! T. basis at or below refGPH

  end type ForwardModelIntermediate_T

  ! Workspace type stuff for the scan model
  type (ForwardModelIntermediate_T), private, save :: IFM

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  ! Define various constants etc.
  ! Some terms to do with refraction
  
  real (r8), parameter :: MaxPressure = 1400.0 ! /mb Don't allow very large pressures
  real (r8), parameter :: MinZeta = log10(maxPressure)

contains ! =============== Subroutines and functions ==========================

  ! -----------------------------  DestroyForwardModelIntemediate  -----
  subroutine DestroyForwardModelIntermediate

    ! Exectuable code

    call Deallocate_test ( ifm%basisGPH, 'basisGPH', ModuleName )
    call Deallocate_test ( ifm%RT, 'RT', ModuleName )
    call Deallocate_test ( ifm%R, 'R', ModuleName )

  end subroutine DestroyForwardModelIntermediate

  ! ----------------------------------------  DumpInstanceWindows  -----
    ! Dump the instance windows (e.g., in temperature) 
    ! that a 2d forward model must consider
    ! for the range of mafs [maf1, maf2]
  subroutine DumpInstanceWindows ( STATE, EXTRA, MAF1, MAF2, FMCONF, DETAILS )
    use Dump_0, only: Dump
    use ForwardModelConfig, only: Dump
    use MLSStrings, only: WriteIntsToChars
    use HighOutput, only: HeadLine, OutputNamedValue
    ! Args
    type (Vector_T), intent(in) :: STATE ! The state vector
    type (Vector_T), intent(in) :: EXTRA ! Other stuff in the state vector
    integer, intent(in)         :: MAF1
    integer, intent(in)         :: MAF2
    type (ForwardModelConfig_T), intent(in) :: FMCONF
    integer, intent(in)         :: DETAILS ! dump FMConf if > 0
    ! Internal variables
    character(len=4) :: CHARS
    character(len=16) :: HEADING
    integer :: INSTLOW
    integer :: INSTHIGH
    type (VectorValue_T), pointer :: PHITAN ! Phitan component of state
    type (VectorValue_T), pointer :: TEMP   ! Temperature component of state
    integer :: WINDOWSTART_T         ! first instance for temperature
    integer :: WINDOWFINISH_T        ! last instance for temperature
    ! Executable
    if ( fmConf%instrumentModule > 0 ) then
      phitan => getQuantityForForwardModel ( state, extra, &
        & quantityType=l_phitan, instrumentModule=fmConf%instrumentModule, &
        & config=fmConf )
    else
      phitan => getQuantityForForwardModel ( state, extra, &
        & quantityType=l_phitan, config=fmConf )
    end if
    temp => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_temperature, config=fmConf )
    call FindInstanceWindow( temp, phitan, maf1, FMConf%phiWindow, &
      & FMConf%windowUnits, windowstart_t, windowfinish_t )
    instLow = min( windowstart_t, windowfinish_t )
    instHigh = max( windowstart_t, windowfinish_t )
    call FindInstanceWindow( temp, phitan, maf2, FMConf%phiWindow, &
      & FMConf%windowUnits, windowstart_t, windowfinish_t )
    instLow = min( windowstart_t, instLow )
    instHigh = max( instHigh, windowfinish_t )
    call writeintstochars( maf1, chars )
    heading = adjustl(chars)
    call writeintstochars( maf2, chars )
    heading = trim(heading) // ',' // adjustl(chars)
    call headline( 'Instance Window for mafs in range: '// trim(heading), &
      & fillChar='-', before='*', after='*' )
    call dump( FMConf, details=details-1 )
    call outputNamedValue ( 'Instance Window', [ instLow, instHigh ] )
  end subroutine DumpInstanceWindows

  ! ------------------------------------------------  GetBasisGPH  -----
  subroutine GetBasisGPH ( temp, refGPH, gph, R, RT, belowRef )
    ! This function calculates gph based on
    ! temperature temp
    ! reference gph refgph
    !
    ! It also marks with FillValues any levels affected
    ! by FillValues among the Temperatures
    use Dump_0, only: Dump
    use HighOutput, only: OutputNamedValue
    use MLSCommon, only: UndefinedValue
    use MLSFillValues, only: IsFillValue
    use MLSFinds, only: FindFirst
    use Physics, only: Boltz
    use Trace_M, only: Trace_Begin, Trace_End

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: TEMP ! The temperature field
    type (VectorValue_T), intent(IN) :: REFGPH ! The reference gph field
    real (r8), dimension(:,:), intent(OUT) :: GPH ! Result (temp%noSurfs,temp%noInstances)
    real (r8), dimension(:), pointer, optional :: R ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer, optional :: RT ! R*T (noSurfs,noInstances)
    integer, intent(out), optional :: BELOWREF ! Level in temperature basis
  
    ! Some terms to do with the gas 'constant'

!    real (r8), parameter :: GASM0 = 25.34314957d-3 ! Constant
!    real (r8), parameter :: GASM1 = 2.89644d-3 ! Linear term
!    real (r8), parameter :: GASM2 = -0.579d-3 ! Quadratic term
!    real (r8), parameter :: GASR0 = 8.31441 ! `Standard' gas constant
    REAL(r8), parameter :: c_mass = 0.02_rp ! mass reduction coefficient
    REAL(r8), parameter :: basiscutoff = 2.5_rp ! zeta where mass reduction begins

    ! Local variables, many automatic arrays

    real (r8), dimension(:), pointer :: MYR      ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer :: MYRT   ! R*T (noSurfs,noInstances)

    real (r8), dimension(temp%template%noSurfs) :: LOGP ! -log10 pressure
    real (r8), dimension(temp%template%noSurfs) :: MODIFIEDBASIS ! noSurfs
    real (r8), dimension(temp%template%noInstances) :: CURRENTREFGPH ! From 1st calc
    real (r8), dimension(temp%template%noInstances) :: CORRECTION ! To apply to gph
    real (r8), dimension(temp%template%noInstances) :: DELTAGEOPOT ! noInstances

    integer :: Me = -1                  ! String index for trace cacheing

    integer :: MYBELOWREF               ! Result of a hunt
    real (r8) :: ABOVEREFWEIGHT         ! Interpolation weight
    
    integer :: INSTANCE                 ! Loop counter
    integer :: SURF                     ! Loop counter
    
!    real (r8) :: BASISCUTOFF            ! Threshold level for gas constant
    real (r8) :: REFLOGP                ! Log p of pressure reference surface
    real (r8) :: BASISGAP               ! Space between adjacent surfaces
    logical, parameter :: DEEBUG = .false.

    call trace_begin ( me, 'GetBasisGPH', cond=.false. )
    nullify ( myR, myRT )
    ! Check that we get the right kinds of quantities
    if ( ( .not. ValidateVectorQuantity( temp,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) .or. &
      &  ( .not. ValidateVectorQuantity( refGPH,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) ) &
      & call MLSMessage(MLSMSG_Error,ModuleName,&
        & 'Inappropriate temp/refGPH quantity' )

    ! Now allocate intermediate quantities
    call allocate_test ( myR, temp%template%noSurfs, 'myR', ModuleName )
    call allocate_test ( myRT, temp%template%noSurfs, &
      & temp%template%noInstances, 'myR', ModuleName )

    ! Now the main calculation. This is two parts.  First we compute a
    ! geopotential height grid with an arbitrary offset of H=0 at the lowest
    ! basis point.

    ! To do this we get a gas constant for all the temperature basis points
    logP= temp%template%surfs(:,1)

!    basisCutoff= -gasM1 / (2*gasM2) ! Set a threshold value
    modifiedBasis = max ( logP, basisCutoff ) ! Either logP or this threshold
!    myR = gasR0 / ( gasM0 + gasM1*modifiedBasis + gasM2*modifiedBasis**2 )
    myR = boltz * (1.0_r8 + c_mass*(modifiedBasis - basiscutoff)**2)   

    ! Compute R*T for each point, avoid spread intrinsic to save memory
    ! and cpu time.
    do instance = 1, temp%template%noInstances
      myRT(:,instance) = myR * temp%values(:,instance)
    end do

    gph(1,:) = 0.0
    do surf = 2, temp%template%noSurfs
!     deltaGeopot = (ln10/ (2*g0) ) * &
!         & ( myRT(surf,:) + myRT(surf-1,:) ) * &
!         & ( logP(surf) - logP(surf-1) )
      deltaGeopot = ( myRT(surf,:) + myRT(surf-1,:) ) * &
         & ( logP(surf) - logP(surf-1) ) / (2 * g0)
      gph(surf,:) = gph(surf-1,:) + deltaGeopot
    end do

    ! Now we need to correct for the reference geopotential, find the layer the
    ! reference surface is within.

    refLogP = refGPH%template%surfs(1,1)
    call Hunt ( logP, refLogP, myBelowRef )

    ! Get weights
    basisGap = logP(myBelowRef+1) - logP(myBelowRef)
    aboveRefWeight = ( refLogP - logP(myBelowRef) )/basisGap

    ! Forbid extrapolation
    aboveRefWeight = max ( min ( aboveRefWeight, 1.0D0 ), 0.0D0 )

    ! Get geopotential at the reference surface from our intermediate result
!    currentRefGPH = gph(myBelowRef,:) + ((basisGap*ln10)/(2*g0))* &
!      & ( myRT(myBelowRef,:) * aboveRefWeight * (2-aboveRefWeight) + &
!      &   myRT(myBelowRef+1,:) * (aboveRefWeight**2))
    currentRefGPH = gph(myBelowRef,:) + ((basisGap)/(2 * g0))* &
      & ( myRT(myBelowRef,:) * aboveRefWeight * (2-aboveRefWeight) + &
      &   myRT(myBelowRef+1,:) * (aboveRefWeight**2))

    ! Now make the correction, again avoid spread to save time/memory,
    ! Also convert to geopotential height.
    correction = refGPH%values(1,:) - currentRefGPH
    do surf = 1, temp%template%noSurfs         
       gph(surf,:) = gph(surf,:) + correction
    end do
    ! Put back intermediate data if wanted.
    if ( present(R) ) then
      R=>myR
    else
      call deallocate_test ( myR, 'myR', moduleName )
    end if

    if ( present(rt) ) then
      rt=>myRT
    else
      call deallocate_test ( myRT, 'myRT', moduleName )
    end if

    if ( present(belowRef) ) belowRef = myBelowRef
    if ( any(IsFillValue(temp%values) ) ) then
      if ( DEEBUG ) call outputNamedValue ( 'myBelowRef', myBelowRef )
      ! We will reset gph to FillValue for any ..
      do instance=1, temp%template%noInstances
        if ( DEEBUG ) call outputNamedValue ( 'instance', instance )
        if ( DEEBUG ) call dump( temp%values(:, instance), 'Temperatures' )
        ! heights above the 1st neg T starting at the refGPH and going up
        surf = FindFirst ( isFillValue(temp%values(myBelowRef:, instance) ) )
        if ( DEEBUG ) call outputNamedValue ( 'upper -999.99 surf', surf )
        if ( surf > 0 ) gph(myBelowRef+surf-1:,instance) = UndefinedValue
        ! heights below the 1st neg T starting at the refGPH and going down
        surf = FindFirst ( isFillValue(temp%values(myBelowRef:1:-1, instance) ) )
        if ( DEEBUG ) call outputNamedValue ( 'lower -999.99 surf', surf )
        if ( surf > 0 ) gph(myBelowRef-surf+1:1:-1,instance) = UndefinedValue
      enddo
    endif
         
    call trace_end ( 'GetBasisGPH', cond=.false. )
    ! That's it  
  end subroutine GetBasisGPH

  ! --------------------------------------------  GetGPHPrecision  -----
  subroutine GetGPHPrecision ( tempPrec, refGPHPrec, &
    & gphPrec )
    ! This function takes a state vector, containing one and only one
    ! temperature and reference geopotential height precisions, and
    ! returns the GPH precision.
    use HighOutput, only: OutputNamedValue
    use Physics, only: Boltz
    use Trace_M, only: Trace_Begin, Trace_End

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: TEMPPREC ! The temperature precision
                                                 ! field
    type (VectorValue_T), intent(IN) :: REFGPHPREC ! The reference gph
                                                   ! precision field
    real (r8), dimension(:,:), intent(OUT) :: GPHPREC ! Result (temp%noSurfs,temp%noInstances)
  
    ! Some terms to do with the gas 'constant'

!    real (r8), parameter :: GASM0 = 25.34314957d-3 ! Constant
!    real (r8), parameter :: GASM1 = 2.89644d-3 ! Linear term
!    real (r8), parameter :: GASM2 = -0.579d-3 ! Quadratic term
!    real (r8), parameter :: GASR0 = 8.31441 ! `Standard' gas constant
    real(r8), parameter :: c_mass = 0.02_rp ! mass reduction coefficient
    real(r8), parameter :: basiscutoff = 2.5_rp ! zeta where mass reduction begins

    ! Local variables, many automatic arrays

    real (r8), dimension(tempPrec%template%noSurfs) :: MYR ! Gas constant

    real (r8), dimension(tempPrec%template%noSurfs) :: LOGP ! -log10 pressure
    real (r8), dimension(tempPrec%template%noSurfs) :: MODIFIEDBASIS ! noSurfs
!   real (r8), dimension(tempPrec%template%noInstances) :: CURRENTREFGPH ! From 1st calc
!   real (r8), dimension(tempPrec%template%noInstances) :: CORRECTION ! To apply to gph

    ! The derivatives are effectively 1D at each instance
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: DMyRT_dT   ! d(R*T)/dT
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: DCurrentRefGPH_dT ! From 1st calc
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: Dcorrection_dT ! To apply to dgph_dT
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: DDeltaGeoPot_dT ! 
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: DGPH_dT ! 
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: DGPH_dRefGPH ! 

    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: GPHPrec2 ! squared PGH precision
    real (r8), dimension(tempPrec%template%noSurfs, &
      & tempPrec%template%noInstances) :: GPHPrec2a ! squared PGH precision (alt)
    ! We usee this next to avoid segment faulting due to bug in NAG compiler
    real (r8), dimension(tempPrec%template%noSurfs) :: TPrecVals ! T Precision

    integer :: Me = -1                  ! String index for trace cacheing

    integer :: MyBelowRef               ! Result of a hunt
    real (r8) :: AboveRefWeight         ! Interpolation weight
    
    integer :: Instance                 ! Loop counter
    integer :: Surf                     ! Loop counter
    
!    real (r8) :: BASISCUTOFF            ! Threshold level for gas constant
    real (r8) :: RefLogP                ! Log p of pressure reference surface
    real (r8) :: BasisGap               ! Space between adjacent surfaces
    integer, dimension( max(1,size(GPHPrec, 1)), max(1,size(GPHPrec, 2)) ) &
      &       :: ItsSign                ! < 0 if tempprec or refgphprec are
    integer, save :: TimesHere = 0
    logical :: DeeBug

    TimesHere = TimesHere + 1
    DeeBug = .false. ! ( TimesHere > 1 )
    if ( DeeBug ) call output( 'Now in GetGPHPrecision', advance='yes' )
    call trace_begin ( me, 'GetGPHPrecision', cond=.false. )
    ! Check that we get the right kinds of quantities
    if ( ( .not. ValidateVectorQuantity( tempPrec,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) .or. &
      &  ( .not. ValidateVectorQuantity( refGPHPrec,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) ) &
      & call MLSMessage(MLSMSG_Error,ModuleName,&
        & 'Inappropriate tempPrec/refGPHPrec quantity' )

    ! Now the main calculation. This is two parts.  First we compute a
    ! geopotential height grid with an arbitrary offset of H=0 at the lowest
    ! basis point.

    ! To do this we get a gas constant for all the temperature basis points
    if ( DeeBug )  then
      call Dump( tempPrec%template )
      call Dump( refGPHPrec%template )
      call output( 'Now the main calculation', advance='yes' )
    end if
    logP= tempPrec%template%surfs(:,1)

!    basisCutoff= -gasM1 / (2*gasM2) ! Set a threshold value
    modifiedBasis = max ( logP, basisCutoff ) ! Either logP or this threshold
!    myR = gasR0 / ( gasM0 + gasM1*modifiedBasis + gasM2*modifiedBasis**2 )
    myR = boltz * (1.0_r8 + c_mass*(modifiedBasis - basiscutoff)**2)   

    ! Compute T derivative of R*T for each point.
    dmyRT_dT = 0.0
    do surf = 1, tempPrec%template%noSurfs
      dmyRT_dT(surf,surf,:) = myR(surf)
    end do

    dgph_dT = 0.0
    do surf = 2, tempPrec%template%noSurfs
!       ddeltaGeopot_dT = (ln10/ (2*g0) ) * &
!         & ( dmyRT_dT(surf,:,:) + dmyRT_dT(surf-1,:,:) ) * &
!         & ( logP(surf) - logP(surf-1) )
       ddeltaGeopot_dT = ( dmyRT_dT(surf,:,:) + dmyRT_dT(surf-1,:,:) ) * &
         & ( logP(surf) - logP(surf-1) ) / (2 * g0)
       dgph_dT(surf,:,:) = dgph_dT(surf-1,:,:) + ddeltaGeopot_dT
    end do

    ! Now we need to correct for the reference geopotential, find the layer the
    ! reference surface is within.

    refLogP = refGPHPrec%template%surfs(1,1)
    call Hunt ( logP, refLogP, myBelowRef )

    ! Get weights
    basisGap = logP(myBelowRef+1) - logP(myBelowRef)
    aboveRefWeight = ( refLogP - logP(myBelowRef) )/basisGap

    ! Forbid extrapolation
    aboveRefWeight = max ( min ( aboveRefWeight, 1.0D0 ), 0.0D0 )

    ! Get derivative of the geopotential at the reference surface
    ! from our intermediate result
!    dcurrentRefGPH_dT = dgph_dT(myBelowRef,:,:) + ((basisGap*ln10)/(2*g0))* &
!      & ( dmyRT_dT(myBelowRef,:,:) * aboveRefWeight * (2-aboveRefWeight) + &
!      &   dmyRT_dT(myBelowRef+1,:,:) * (aboveRefWeight**2))
    dcurrentRefGPH_dT = dgph_dT(myBelowRef,:,:) + ((basisGap)/(2*g0))* &
      & ( dmyRT_dT(myBelowRef,:,:) * aboveRefWeight * (2-aboveRefWeight) + &
      &   dmyRT_dT(myBelowRef+1,:,:) * (aboveRefWeight**2))

    ! Now make the correction, again avoid spread to save time/memory,
    ! Also convert to geopotential height.
    dcorrection_dT = - dcurrentRefGPH_dT
    do surf = 1, tempPrec%template%noSurfs         
       dgph_dT(surf,:,:) = dgph_dT(surf,:,:) + dcorrection_dT
    end do

    ! The refGPH derivative is easy!
    dgph_drefGPH = 1.0
    if ( DeeBug ) then
      call output( 'Transform the Temp, refGPH precision to GPH precision', advance='yes' )
      call outputNamedValue ( 'shape(tempPrec)', shape(tempPrec%values) )
      call outputNamedValue ( 'shape(refGPHPrec)', shape(refGPHPrec%values) )
      call outputNamedValue ( 'shape(GPHPrec)', shape(GPHPrec) )
      call outputNamedValue ( 'shape(ITSSIGN)', shape(ITSSIGN) )
      call outputNamedValue ( 'shape(GPHPrec2)', shape(GPHPrec2) )
      call outputNamedValue ( 'shape(dgph_dT)', shape(dgph_dT) )
      call outputNamedValue ( 'shape(dgph_drefGPH)', shape(dgph_drefGPH) )
      call outputNamedValue ( 'noInstances', tempPrec%template%noInstances )
      call outputNamedValue ( 'noSurfs', tempPrec%template%noSurfs )
      call outputNamedValue ( 'myBelowRef', myBelowRef )
      call Dump( tempPrec%values, 'temperature precision' )
    end if
    !  Transform the temperature, refGPH precision to GPH precision
    GPHPrec2 = 0.0
    ItsSign  = 1
    do instance = 1, tempPrec%template%noInstances
      TPrecVals = tempPrec%values(:,instance)
      do surf = 1, tempPrec%template%noSurfs
        GPHPrec2(:,instance) = GPHPrec2(:,instance) + &
          & ( dgph_dT(:,surf,instance) * TPrecVals(surf) )**2
      end do
      GPHPrec2(:,instance) = GPHPrec2(:,instance) + &
        & ( dgph_drefGPH(:,instance) * refGPHPrec%values(1,instance) )**2
      if ( DeeBug ) then
        call outputnamedValue ( 'instance', instance )
        call outputnamedValue ( 'refGPHPrec%values(1,instance)', refGPHPrec%values(1,instance) )
        call outputnamedValue ( ' all( TPrecVals >= 0._rv',  all( TPrecVals >= 0._rv) )
      end if
      if ( refGPHPrec%values(1,instance) < 0._rv ) then
        ItsSign(:, instance) = -1
      else if ( all( TPrecVals >= 0._rv ) ) then
        ! ITSS already 1
      else if ( refGPHPrec%values(1,instance) >= 0._rv ) then
        if ( DeeBug ) then
          call outputnamedValue ( 'myBelowRef', myBelowRef )
          call outputnamedValue ( 'tempPrec%template%noSurfs', tempPrec%template%noSurfs )
          call outputnamedValue ( 'any( TPrecVals < 0._rv )', any( TPrecVals < 0._rv ) )
          call outputnamedValue ( 'shape( TPrecVals(myBelowRef:myBelowRef))', &
            & shape( TPrecVals(myBelowRef:myBelowRef)) )
        end if
        if ( any( TPrecVals < 0._rv ) ) then
          if ( DeeBug ) then
            call outputnamedValue( 'shape', shape(TPrecVals(1:myBelowRef)))
            call outputnamedValue( 'value', TPrecVals(1:myBelowRef))
          end if
          do surf = 1, myBelowRef
            ! if ( any( tempPrec%values(surf:myBelowRef, instance) < 0._rv ) ) &
            !    & ITSSIGN(surf, instance) = -1
            ITSSIGN(surf, instance) = sign ( 1.d0, minval(TPrecVals(surf:myBelowRef)) )
          end do
          do surf = myBelowRef+1, tempPrec%template%noSurfs
            if ( any( TPrecVals(myBelowRef:surf) < 0._rv ) ) &
              ITSSIGN(surf, instance) = -1
          end do
        end if
      else
        if ( DeeBug ) then
          call output( 'in else clause', advance='yes' )
          call outputnamedValue ( 'myBelowRef', myBelowRef )
          call outputnamedValue ( 'tempPrec%template%noSurfs', tempPrec%template%noSurfs )
        end if
        do surf = 1, tempPrec%template%noSurfs
          ! On which side of the reference surface are we?
          if ( surf < myBelowRef ) then
            if ( any( TPrecVals(surf:myBelowRef) < 0._rv ) ) &
              & ITSSIGN(surf, instance) = -1
          else if ( surf == myBelowRef ) then
            if ( TPrecVals(surf) < 0._rv ) &
              & ITSSIGN(surf, instance) = -1
          else if ( &
            & any( TPrecVals(myBelowRef:surf) < 0._rv ) ) then
            ITSSIGN(surf, instance) = -1
          end if
        end do
      end if
    end do

    ! if ( TimesHere > 1 ) then
    !  call output( '2nd time, so quitting', advance='yes' )
    !  stop
    ! end if
    if ( DeeBug ) call output( 'Finish the calculation', advance='yes' )
    do surf = 1, tempPrec%template%noSurfs
      GPHPrec2a(surf,:) = sum ( &
        & dgph_dT(surf,:,:) * tempPrec%values(:,:) * &
        & tempPrec%values(:,:) * dgph_dT(surf,:,:), &
        & dim=1 )
      GPHPrec2a(surf,:) = GPHPrec2a(surf,:) + &
        & dgph_drefGPH(surf,:) * refGPHPrec%values(1,:) * &
        & refGPHPrec%values(1,:) * dgph_drefGPH(surf,:)
    end do
    GPHPrec = ITSSIGN * sqrt ( GPHPrec2 )
     
    call trace_end ( 'GetGPHPrecision', cond=.false. )
    ! That's it  
  end subroutine GetGPHPrecision

  ! ----------------------------  Get2DHydroStaticTangentPressure  -----
  subroutine Get2DHydrostaticTangentPressure ( ptan, temp, refGPH, h2o, &
    & orbIncl, phiTan, geocAlt, maxIterations, phiWindow, phiWindowUnits, &
    & chunkNo )
    ! This is a new pressure guesser routine which works by simply running the
    ! new 2D scan model 'backwards'
    ! Dummy arguments

    use Dump_0, only: Dump
    use Output_M, only: Output
    use QuantityTemplates, only: QuantityTemplate_T

    type (VectorValue_T), intent(inout) :: PTAN ! Tangent pressure sv. component
    type (VectorValue_T), intent(in) :: Temp ! Temperature
    type (VectorValue_T), intent(in) :: RefGPH ! Reference GPH
    type (VectorValue_T), intent(in) :: H2O ! H2O
    type (VectorValue_T), intent(in) :: Orbincl ! Inclination
    type (VectorValue_T), intent(in) :: PhiTan ! Tangent phi
    type (VectorValue_T), intent(in) :: GeocAlt ! L1B Tp Geoc alt.
    integer, intent(IN) :: MaxIterations ! Number of iterations to use
    real(r8), intent(in) :: PhiWindow(2)   ! For 2D or not
    integer, intent(in) :: PhiWindowUnits
    integer, intent(in), optional :: ChunkNo

    ! Local parameters
    integer, parameter :: NoQtys = 8

    ! Local variables
    type (QuantityTemplate_T) :: ScanResidual ! Scan residual
    type (QuantityTemplate_T) :: MyQTs(NoQtys) ! All our quantities
    type (Vector_T) :: State            ! Gets ptan
    type (Vector_T) :: Extra            ! Gets temp,refGPH,h2o and geocAlt
    type (Vector_T) :: Residual         ! Gets the scan residual
    type (VectorTemplate_T) :: StateTemplate ! Gets ptan
    type (VectorTemplate_T) :: ExtraTemplate ! Gets temp,refGPH,h2o and geocAlt
    type (VectorTemplate_T) :: ResidualTemplate ! Gets the scan residual
    type (Matrix_T) :: Jacobian         ! dScanresidual/dPtan

    type (ForwardModelConfig_T) :: FMConf
    type (ForwardModelStatus_T) :: FMStat

    integer :: I                        ! Iteration counter
    integer :: MAF                      ! MAF counter
    integer :: Me = -1                  ! String index for trace
    real (r8), dimension(ptan%template%noSurfs, ptan%template%noInstances) :: &
      EARTHRADIUS ! Earth radius (minor frame)

    ! Executable code
    ! Get a really simple first guess, assuming a uniform log scale height of
    ! 16km

    call trace_begin (me, 'Get2DHydrostaticTangentPressure', cond=toggle(emit) )

    earthRadius = earthRadA*earthRadB/sqrt( &
      & (earthRadA**2-earthRadB**2)*sin(GeodToGeocLat(ptan%template%geodLat))**2 + earthRadB**2)
    ptan%values = -3.0+(geocAlt%values - earthRadius)/16.0e3
    if ( any(ptan%values < -4.0 .or. ptan%values > 10.0) ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'Why is PTan weird?  Module turned off?  PTan guess set to -3.0' )
      ptan%values = -3.0
    end if
    if ( .not. associated(orbIncl%values) ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'orbIncl not associated; non-satellite data?' )
      call trace_end ( 'Get2DHydrostaticTangentPressure', cond=toggle(emit) )
      return
    else if ( size(orbIncl%values) < 1 ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & 'orbIncl array size 0; non-satellite data?' )
      call trace_end ( 'Get2DHydrostaticTangentPressure', cond=toggle(emit) )
      return
    end if
    
    if ( switchDetail ( switches, 'pguess' ) > -1 ) then
      call dump ( ptan%values, 'Initial ptan guess' )
      if ( switchDetail ( switches, 'pguess' ) > 1 ) then
        call dump ( geocAlt%values, name='GeocAlt' )
        call dump ( earthRadius, name='EarthRadius' )
      end if
    end if

    fmConf%phiWindow = phiWindow
    fmConf%windowUnits = phiWindowUnits
    fmConf%instrumentModule = ptan%template%instrumentModule
    fmConf%differentialScan = .false.

    call Allocate_test ( fmStat%rows, ptan%template%noInstances, &
      & 'fmStat%rows', ModuleName )
    fmStat%rows = .false.

    ! Construct a scan residual quantity by hand
    scanResidual = ptan%template
    scanResidual%quantityType = l_scanResidual
    scanResidual%unit = PHYQ_Length
    ! Now create an array of all our quantity templates
    myQTs = [ ptan%template, temp%template, refGPH%template, &
      & h2o%template, orbIncl%template, phiTan%template, &
      & geocAlt%template, scanResidual ]
    ! Now create the vector templates
    ! state: ptan
    call ConstructVectorTemplate ( 0, myQTs, [1], stateTemplate, forWhom=moduleName )
    ! extra: temp, refGPH, h2o, geocAlt
    call ConstructVectorTemplate ( 0, myQTs, [2,3,4,5,6,7],  extraTemplate, forWhom=moduleName )
    ! residual
    call ConstructVectorTemplate ( 0, myQTs, [8], residualTemplate, forWhom=moduleName )

    ! Now create the vectors
    state = CreateVector ( 0, stateTemplate, myQTs )
    extra = CreateVector ( 0, extraTemplate, myQTs )
    residual = CreateVector ( 0, residualTemplate, myQTs )

    ! Fill the values
    state%quantities(1)%values = ptan%values

    extra%quantities(1)%values = temp%values
    extra%quantities(2)%values = refGPH%values
    extra%quantities(3)%values = h2o%values
    extra%quantities(4)%values = orbIncl%values
    extra%quantities(5)%values = phiTan%values
    extra%quantities(6)%values = geocAlt%values

    ! Create the matrix
    call CreateEmptyMatrix ( jacobian, 0, residual, state )

    ! Now do the iterations
    do i = 1, maxIterations
      ! Get residual for all mafs
      do maf = 1, ptan%template%noInstances
        fmStat%maf = maf
        call TwoDScanForwardModel ( fmConf, state, extra, residual, &
          & fmStat, jacobian, chunkNo )
        state%quantities(1)%values(:,maf) = state%quantities(1)%values(:,maf) - &
          & residual%quantities(1)%values(:,maf) / &
          &   jacobian%block(maf,maf)%values(:,1)
        call ClearMatrix ( jacobian )
      end do
      if ( switchDetail ( switches, 'pguess' ) > -1 ) then
        call output ( 'Pressure guesser iteration ' )
        call output ( i )
        call output ( ' mean, min abs, max abs residual:  ', advance='yes' )
        call output ( sum(residual%quantities(1)%values) / &
          & (ptan%template%noInstances*ptan%template%noSurfs) )
        call output ( minval(abs(residual%quantities(1)%values)), before=', ' )
        call output ( maxval(abs(residual%quantities(1)%values)), before=', ', &
          & advance='yes' )
      end if
    end do

    ! Put the result back in state
    ptan%values = state%quantities(1)%values

    if ( switchDetail (switches, 'pguess' ) > -1 ) &
      & call dump ( ptan%values, 'Final ptan guess' )

    ! Destroy our vectors
    call DestroyMatrix ( jacobian )
    call DestroyVectorInfo ( state )
    call DestroyVectorInfo ( extra )
    call DestroyVectorInfo ( residual )
    
    call Deallocate_test ( fmStat%rows, 'fmStat%rows', ModuleName )

    call trace_end ( 'Get2DHydrostaticTangentPressure', cond=toggle(emit) )

  end subroutine Get2DHydrostaticTangentPressure

  ! ------------------------------  GetHydroStaticTangentPressure  -----
  subroutine GetHydrostaticTangentPressure ( ptan, temp, refGPH, h2o, geocAlt, &
    maxIterations )
    ! This routine is a pressure `guesser'.  It works by comparing the
    ! geopotential heights estimated from the level 1 heights with those based on
    ! a hydrostatic calculation.  The tangent point pressure is adjusted over
    ! several iterations until a good match is achieved.

    ! Dummy arguments
    ! Temp and H2O have the target attribute so we can take a pointer to
    ! their %template%surfs components.
    type (VectorValue_T), intent(inout) :: PTAN ! Tangent pressure sv. component
    type (VectorValue_T), intent(in), target :: Temp ! Temperature
    type (VectorValue_T), intent(in) :: RefGPH       ! Reference GPH
    type (VectorValue_T), intent(in), target :: H2O  ! H2O
    type (VectorValue_T), intent(in) :: GeocAlt      ! L1B Tp Geoc alt.
    integer, intent(IN) :: MaxIterations             ! Number of iterations to use

    ! Local variables

    real (r8), dimension(:,:), pointer :: RT           ! rt=R*T
    real (r8), target :: BasisGPH(temp%template%noSurfs, temp%template%noInstances)

    real (r8) :: ACoeff(ptan%template%noSurfs)           ! Quadratic term
    real (r8) :: BasisLower(ptan%template%noSurfs)       ! For temperature
    real (r8) :: BasisSpacing(ptan%template%noSurfs)     ! For temperature
    real (r8) :: BasisUpper(ptan%template%noSurfs)       ! For temperature
    real (r8) :: BCoeff(ptan%template%noSurfs)           ! Quadratic term
    real (r8) :: CCoeff(ptan%template%noSurfs)           ! Quadratic term
    real (r8) :: DeltaRT(ptan%template%noSurfs) 
    real (r8) :: EarthRadius(ptan%template%noSurfs)
    real (r8) :: GeocLat(ptan%template%noSurfs)
    real (r8) :: GeometricGPH(ptan%template%noSurfs)
    real (r8) :: N(ptan%template%noSurfs)                ! Refractive index
    real (r8) :: P2(ptan%template%noSurfs)
    real (r8) :: P4(ptan%template%noSurfs)
    real (r8) :: PointingH2O(ptan%template%noSurfs)      ! t.p. h2o
    real (r8) :: PointingPres(ptan%template%noSurfs)     ! t.p. press
    real (r8) :: PointingTemp(ptan%template%noSurfs)     ! t.p. temp.
    real (r8) :: Ratio2(ptan%template%noSurfs)           ! minor frame
    real (r8) :: Ratio4(ptan%template%noSurfs)           ! minor frame
    real (r8) :: RefractedGeocAlt(ptan%template%noSurfs) ! minor frame
    real (r8) :: RTLower(ptan%template%noSurfs)
    real (r8) :: RTUpper(ptan%template%noSurfs)
    real (r8) :: S2(ptan%template%noSurfs)

    real (r8), dimension(:), pointer :: H2OBasis         ! H2O basis
    real (r8), dimension(:), pointer :: H2OVals          ! Part of h2o
    real (r8), dimension(:), pointer :: PTANVals         ! Part of ptan
    real (r8), dimension(:), pointer :: TempBasis        ! Temp basis
    real (r8), dimension(:), pointer :: TempVals         ! Part of temp
    real (r8), dimension(:), pointer :: ThisBasisGPH     ! Part of basisGPH

    integer :: ClosestTempProfiles(ptan%template%noInstances)
    integer :: ClosestH2OProfiles(ptan%template%noInstances)
    integer :: Lower(ptan%template%noSurfs) ! index into temperature profile
    integer :: Upper(ptan%template%noSurfs) ! index into temperature profile

    integer :: Iteration                ! Loop stuff
    integer :: MAF                      ! Loop counter
    integer :: Me = -1                  ! String index for trace

    ! Executable code

    call trace_begin ( me, 'GetHydrostaticTangentPressure', cond=toggle(emit) )

    ! Check that we get the right kinds of quantities
    if ( ( .not. ValidateVectorQuantity( temp,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) .or. &
      &  ( .not. ValidateVectorQuantity( refGPH,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) .or. &
      &  ( .not. ValidateVectorQuantity( h2o,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=[ l_Zeta ]) ) ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
      &   'Inappropriate temp/refGPH/h2o quantity' )

    ! Now precompute the basis geopotential height
    nullify ( rt ) ! Allocated in GetBasisGPH
    call GetBasisGPH ( temp, refGPH, basisGPH, rt=rt )

    ! Find the closest temperature and h2o profiles to each MAF.  Note that
    ! this makes this calculation 1D
    call FindClosestInstances ( temp, ptan, closestTempProfiles )
    call FindClosestInstances ( h2o, ptan, closestH2OProfiles )

    ! Rather than try to be too clever and do this all with 2D arrays, we'll
    ! loop over major frame.

    tempBasis => temp%template%surfs(:,1)
    h2oBasis => h2o%template%surfs(:,1)

    do maf = 1, ptan%template%noInstances
      ! Setup pointers for this maf
      ptanVals => ptan%values(:,maf)
      tempVals => temp%values(:,closestTempProfiles(maf))
      h2oVals => h2o%values(:,closestH2OProfiles(maf))
      thisBasisGPH => basisGPH(:,closestTempProfiles(maf))

      ! Get a really simple first guess, assuming a uniform log scale height of
      ! 16km
      geocLat=GeodToGeocLat(ptan%template%geodLat(:,maf))
      s2=sin(geocLat)**2
      earthRadius = earthRadA*earthRadB/sqrt(&
        & (earthRadA**2-earthRadB**2)*s2 + earthRadB**2)

      ptanVals = -3.0+(geocAlt%values(:,maf) - &
        & earthRadius)/16e3 ! Do a better job later !???
      ! Precompute some spherical geometry terms, these are
      ! particularly worthy of note. We're doing the geometricGPH
      ! in a 2D manner (apart from refraction effects).  The hydrostatic GPH is
      ! being done 2D.
      p2=0.5*(3*s2-1)
      p4=0.125*(35*(s2**2)-30*s2+3)

      ! Now we have an iteration loop where we try to fit the tangent pressure
      do iteration = 1, maxIterations
        ! The first stage is to refract the tangent point altitudes, for this we
        ! need the refractive index, which is dependent on temperature, pressure
        ! and water vapor concentration for the given altitude.

        ! Note here an assumption that temp and h2o are coherent.
        call InterpolateValues ( tempBasis, tempVals, ptanVals, &
          pointingTemp, "Linear", extrapolate='Clamp' )
        call InterpolateValues ( h2oBasis, h2oVals, ptanVals, &
          pointingH2O, "Linear", extrapolate='Clamp' )

        pointingH2O = max(pointingH2O, 0.0_r8)
        where ( geocAlt%values(:,maf) /= geocAlt%template%badValue )

          pointingPres = min(10**(-ptanVals), maxPressure)

          n = refractive_index_f ( pointingPres, pointingTemp, pointingH2O )

          n = min(n,maxRefraction) ! Forbid stupidly large values of n

          ! Now we're going to compute refracted geocentric altitudes
          refractedGeocAlt=geocAlt%values(:,maf)/(1.0+n)

          ! Now convert these to geopotential heights
          ratio2=(earthRadA/refractedGeocAlt)**2
          ratio4=ratio2**2

          geometricGPH = -((GM/g0)/refractedGeocAlt)*&
            &              (1-j2*p2*ratio2- j4*p4*ratio4) - &
            ((omega**2/(2*g0))*(refractedGeocAlt*cos(geocLat))**2)+earthSurfaceGPH
        elsewhere
          n = 0.0
        end where

        ! Now, we're effectively going to compare this with a hydrostatic
        ! calculation.

        call Hunt ( thisBasisGPH, geometricGPH, lower )
        upper=lower+1
        basisLower= tempBasis(lower)
        basisUpper= tempBasis(upper)
        basisSpacing= basisUpper - basisLower
        rtLower = rt(lower, closestTempProfiles(maf))
        rtUpper = rt(upper, closestTempProfiles(maf))
        deltaRT = rtUpper - rtLower

        where ( geocAlt%values(:,maf) == geocAlt%template%badValue )
          ptanVals = ptan%template%badValue

          ! Below the bottom
        elsewhere ( geometricGPH < thisBasisGPH(1) )
          ptanVals = tempBasis(1) - (thisBasisGPH(1) - geometricGPH) * &
            & ((g0 / ln10) / rtLower)

          ! Above the top
        elsewhere ( geometricGPH >= thisBasisGPH(temp%template%noSurfs) )
          ptanVals = tempBasis(temp%template%noSurfs) + &
            & (geometricGPH-thisBasisGPH(upper)) * &
            & ((g0 / ln10) / rtUpper)

          ! Everywhere else
        elsewhere
          aCoeff = deltaRT
          bCoeff = 2*rtLower
          cCoeff = -2*(g0/ln10)*(geometricGPH-thisBasisGPH(lower)) / &
            & (basisSpacing)
          ! Let ptan contain upperWeight (ptan-basisLower)/basisSpacing first
          ptanVals = 2*cCoeff/( -bCoeff - &
            & sqrt(max(bCoeff**2-4*aCoeff*cCoeff,0.0_r8)))
          ! The max is here just in case of slips.
          ! Now let it contain ptan
          ptanVals = basisLower + ( ptanVals * basisSpacing )
        end where
      end do                            ! End iteration loop
    end do                              ! Major frame loop

    call Deallocate_Test ( rt, "rt", ModuleName ) ! Was allocated in GetBasisGPH

    call trace_end ( 'GetHydrostaticTangentPressure', cond=toggle(emit) )

  end subroutine GetHydrostaticTangentPressure

  ! -------------------------------------------  ScanForwardModel  -----
  subroutine ScanForwardModel ( fmConf, state, extra, &
    & fwmOut, fmStat, jacobian )
    ! This is the main `scan model' for the module. It compares altitude reported
    ! geometrically, as converted into geopotential height, with a hydrostatic
    ! equivalent.

    ! Note that while this model takes multi-MAF quantities as input, it only
    ! considers one MAF at a time (fmStat%maf).

    ! Dummy arguments
    type (ForwardModelConfig_T), intent(in) :: FMConf ! Configuration options
    type (Vector_T), intent(in) :: State ! The state vector
    type (Vector_T), intent(in) :: Extra ! Other stuff in the state vector
    type (Vector_T), intent(inout) :: FWMOut ! Output vector, residual filled
    type (ForwardModelStatus_T), intent(inout) :: FMStat ! Which maf etc.
    type (Matrix_T), intent(inout), optional :: Jacobian ! The derivative matrix

    ! Local parameters

    character(len=*), parameter :: InvalidQuantity = "Invalid vector quantity for "

    ! Local variables

    integer :: Me = -1                  ! String index for trace
    integer :: NoMAFS                   ! Dimension
    integer :: NoMIFS                   ! Dimension
    integer :: NoTemps                  ! Dimension

    logical :: H2OInState               ! Set if H2O in state, not extra
    logical :: HeightOffsetInState      ! Set if heightOffset in state, not extra
    logical :: PTANInState              ! Set if ptan in state, not extra
    logical :: REFGPHInState            ! Set if refGPH in state not extra
    logical :: TempInState              ! Set if temp in state not extra

    type (VectorValue_T), pointer :: H2O      ! H2O component of state
    type (VectorValue_T), pointer :: HeightOffset ! Height offset component of state
    type (VectorValue_T), pointer :: L1Alt    ! Tangent point altitude
    type (VectorValue_T), pointer :: PTAN     ! Ptan component of state
    type (VectorValue_T), pointer :: RefGPH   ! Ref gph component of state
    type (VectorValue_T), pointer :: Residual ! Resulting component of fwmOut
    type (VectorValue_T), pointer :: Temp     ! Temperature component of state

    ! Executable code -----------------------
    call MLSMessage( MLSMSG_Error, ModuleName,&
      & 'The 1d Scan model is broken; how could you get here?' )
    call trace_begin ( me, 'ScanForwardModel', cond=toggle(emit) )

    ! Identify the vector quantities from state/extra
    temp => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_temperature, config=fmConf, foundInFirst=tempInState )
    ptan => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_ptan, instrumentModule=fmConf%instrumentModule,&
      & config=fmConf, foundInFirst=ptanInState )
    h2o => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_vmr, molecule=l_h2o, &
      & config=fmConf, foundInFirst=h2oInState, noError=.true.)
    refgph => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_refGPH, config=fmConf, foundInFirst=refGPHInState )
    l1Alt => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_tngtGeocAlt, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf )
    heightOffset => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_heightOffset, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf, noError=.true., foundInFirst=heightOffsetInState )

    ! Identify the vector quantities from fwmOut
    residual => getQuantityForForwardModel ( fwmOut, &
      & quantityType=l_scanResidual, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf )

    ! Now check that they make sense
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=[ l_none ]) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=[ l_none ]) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )
    if ( .not. ValidateVectorQuantity(h2o, stacked=.true., coherent=.true., &
      & frequencyCoordinate=[ l_none ]) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'h2o' )
    if ( .not. ValidateVectorQuantity(refGPH, stacked=.true., coherent=.true., &
      & frequencyCoordinate=[ l_none ]) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'refGPH' )
    if ( .not. ValidateVectorQuantity(l1Alt, minorFrame=.true., &
      & frequencyCoordinate=[ l_none ]) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'l1Alt' )
    if ( associated(heightOffset) ) then
      if ( .not. ValidateVectorQuantity(heightOffset, minorFrame=.false., &
        & frequencyCoordinate=[ l_none ],&
        & noInstances=[1] ) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'heightOffset' )
    end if

    ! Now setup some standard dimensions etc.
    noMIFs = ptan%template%noSurfs
    noMAFs = ptan%template%noInstances
    noTemps = temp%template%noSurfs
    call scanForwardModelAuto
    call trace_end ( 'ScanForwardModel', cond=toggle(emit) )

  contains

    subroutine ScanForwardModelAuto

      ! Local variables

      integer :: Col                      ! Block col in jacobian
      integer :: I                        ! Loop index
      integer :: J                        ! Loop index
      integer :: Lower                    ! Index into T profile
      integer :: MAF                      ! Major frame to consider
      integer :: Row                      ! Block row in jacobian
      real (r8), target :: Zero=0.0_r8    ! Need this if no heightOffset
      real (r8), pointer :: HTOff         ! HeightOffset%values or zero
      real (r8) :: UsePTAN                ! A tangent pressure
      real (r8) :: Basis                  ! Index
      real (r8) :: BasisMinus             ! Index
      real (r8) :: BasisPlus              ! Index

      integer, dimension(noMifs) :: PointTempLayer ! Pointing layer
      integer, dimension(noMifs) :: PointH2OLayer  ! Pointing layer
      integer, dimension(noMafs) :: ClosestTempProfiles
      integer, dimension(noMafs) :: ClosestH2OProfiles

      real (r8), dimension(noMifs) :: DH2OByDPtan ! Derivative
      real (r8), dimension(noMifs) :: DHydrosGPHByDPtan ! Derivative
      real (r8), dimension(noMifs) :: DL1GPHByDHeightOffset ! Derivative
      real (r8), dimension(noMifs) :: DL1GPHByDL1RefrGeocAlt ! Derivative
      real (r8), dimension(noMifs) :: DL1GPHByDPtan ! Derivative
      real (r8), dimension(noMifs) :: DL1GPHByDtempLower ! Derivative
      real (r8), dimension(noMifs) :: DL1GPHByDtempUpper ! Derivative
      real (r8), dimension(noMifs) :: DL1RefrGeocAltByDHeightOffset ! Derivative
      real (r8), dimension(noMifs) :: DL1refrGeocAltByDn ! Derivative
      real (r8), dimension(noMifs) :: DL1RefrGeocAltByDptan ! Derivative
      real (r8), dimension(noMifs) :: DL1RefrGeocAltByDtempLower ! Derivative
      real (r8), dimension(noMifs) :: DL1RefrGeocAltByDtempUpper ! Derivative
      real (r8), dimension(noMifs) :: dnByDH2O  ! Derivative
      real (r8), dimension(noMifs) :: DnByDPtan ! Derivative
      real (r8), dimension(noMifs) :: DnByDTemp ! Derivative
      real (r8), dimension(noMifs) :: DnByDTempLower ! Derivative
      real (r8), dimension(noMifs) :: DnByDTempUpper ! Derivative
      real (r8), dimension(noMifs) :: DTempByDPtan ! Derivative
      real (r8), dimension(noMifs) :: GeocLat ! TP GeocLat
      real (r8), dimension(noMifs) :: H2OBasisGap ! Basis spacing
      real (r8), dimension(noMifs) :: H2OLowerWeight ! Weight
      real (r8), dimension(noMifs) :: H2OUpperWeight ! Weight    
      real (r8), dimension(noMifs) :: HydrosGPH ! GPH from hydrostatic calc.
      real (r8), dimension(noMifs) :: L1GPH ! Geometric geopotential
      real (r8), dimension(noMifs) :: L1RefrGeocAlt ! Refracted geocentric alt
      real (r8), dimension(noMifs) :: N ! Refractive index - 1
      real (r8), dimension(noMifs) :: P2 ! Polynomial term
      real (r8), dimension(noMifs) :: P4 ! Polynomial term
      real (r8), dimension(noMifs) :: PointingH2O ! tangent h2o.
      real (r8), dimension(noMifs) :: PointingPres ! tangent pressure.
      real (r8), dimension(noMifs) :: PointingTemp ! tangent temp.
      real (r8), dimension(noMifs) :: Ratio2 ! For geopotential calculation
      real (r8), dimension(noMifs) :: Ratio4 ! For geopotential calculation
      real (r8), dimension(noMifs) :: S2 ! sin^2 geocLat
      real (r8), dimension(noMifs) :: TempBasisGap ! Basis spacing
      real (r8), dimension(noMifs) :: TempLowerWeight ! Weight
      real (r8), dimension(noMifs) :: TempUpperWeight ! Weight

      ! These are all dimensioned noMIFs
      real (r8), dimension(:), pointer :: BasisGPH ! From ifm
      real (r8), dimension(:), pointer :: H2OBasis ! => h2o%template%surfs
      real (r8), dimension(:), pointer :: H2OVals ! => h2o%values
      real (r8), dimension(:), pointer :: PtanVals ! => ptan%values
      real (r8), dimension(:), pointer :: RT ! Gas constant*T
      real (r8), dimension(:), pointer :: TempBasis ! => temp%template%surfs
      real (r8), dimension(:), pointer :: TempVals ! => temp%values

      real (r8), dimension(noMifs+1, noTemps) :: A ! Derivative array,
                                          ! extra row for refGPH
      real (r8), dimension(noMifs, noTemps) :: DHydrosGPHBydTemp ! Derivative array

      type (MatrixElement_T), pointer :: Block ! A matrix block

      maf = fmStat%maf

      ! This could maybe be store in ifm, but for the moment, we'll do it each
      ! time
      call FindClosestInstances ( temp, ptan, closestTempProfiles )
      call FindClosestInstances ( h2o, ptan, closestH2OProfiles )

      ! Make some pointers to save time
      tempVals => temp%values (:,closestTempProfiles(maf) )
      h2oVals => h2o%values(:,closestH2OProfiles(maf) )

      tempBasis => temp%template%surfs(:,1)
      h2oBasis => h2o%template%surfs(:,1)
      ptanVals => ptan%values(:,maf)
      if ( associated( heightOffset ) ) then
        htOff => heightOffset%values(1,1)
      else
        htOff => zero
      end if

      ! Do some preamble calcaulations
      geocLat = GeodToGeocLat( ptan%template%geodLat(:,maf) )
      ! Get terms for geopotential expression
      s2 = (sin(geocLat))**2
      p2 = 0.5 * (3*s2-1) ! Polynomial terms
      p4 = 0.125 * (35*(s2**2) - 30*s2 + 3)

      ! The first part of the calculation is the geometric calculation. This is in
      ! three stages.

      ! First we convert the level 1 geodetic heights to geocentric ones.
      ! Then we refract these heights given our knowledge of temperature and pressure
      ! Lastly we convert these heights to geopotentials.
      ! For each of these stages we need to calculate derivatives


      ! We have to refract the L1altitudes. In order to do that we need the
      ! refractive index, for which we need pressure and the temperature corresponding
      ! to that pressure.

      ! First find the layer each pointing is in, both for T and H2O, and
      ! thence calculate the refractive indicies


      ! Find temp and h2o layers
      call Hunt ( tempBasis, ptanVals, pointTempLayer )
      call Hunt ( h2oBasis,  ptanVals, pointH2OLayer )

      ! Calculate lower and upper weights
      tempBasisGap = tempBasis(pointTempLayer+1) - tempBasis(pointTempLayer)
      h2oBasisGap = h2oBasis(pointH2OLayer+1) - h2oBasis(pointH2OLayer)

      ! Compute weights
      tempUpperWeight = (ptanVals-tempBasis(pointTempLayer))/tempBasisGap
      tempLowerWeight = 1.0 - tempUpperWeight
      h2oUpperWeight = (ptanVals-h2oBasis(pointH2oLayer))/h2oBasisGap
      h2oLowerWeight = 1.0 - h2oUpperWeight

      ! Constrain to constant so we get no extrapolation
      tempUpperWeight = max ( 0.0_r8, min( 1.0_r8, tempUpperWeight ) )
      tempLowerWeight = max ( 0.0_r8, min( 1.0_r8, tempLowerWeight ) )
      h2oUpperWeight = max ( 0.0_r8, min( 1.0_r8, h2oUpperWeight ) )
      h2oLowerWeight = max ( 0.0_r8, min( 1.0_r8, h2oLowerWeight ) )

      ! Now compute derivatives wrt ptan for later chain rule use.
      where (ptanVals>=tempBasis(1) .and. ptanVals<tempBasis(noTemps))
        dTempByDPTan = ( tempVals(pointTempLayer+1) - tempVals(pointTempLayer) ) / &
          & tempBasisGap
      elsewhere
        dTempByDPTan = 0.0_r8
      end where

      where (ptanVals>=h2oBasis(1) .and. ptanVals<h2oBasis(h2o%template%noSurfs))
        dH2OByDPTan = ( h2oVals(pointH2oLayer+1) - h2oVals(pointH2oLayer) ) / &
          & h2oBasisGap
      elsewhere
        dH2OByDPTan = 0.0_r8
      end where

      ! Now get interpolated temperature and H2O
      pointingTemp = tempVals(pointTempLayer) * tempLowerWeight + &
        &            tempVals(pointTempLayer+1) * tempUpperWeight 
      pointingH2O = h2oVals(pointH2OLayer) * h2oLowerWeight + &
        &           h2oVals(pointH2OLayer+1) * h2oUpperWeight

      ! Get pointing pressure in mb
      pointingPres = min(10.0_r8**(-ptanVals), maxPressure)

      ! Now get the refractive index and its derivatives w.r.t. temperature
      ! and H2O.
      call refractive_index_deriv ( pointingPres, pointingTemp, pointingH2O, &
        & n, dnByDTemp, dnByDH2O )

      where ( n > maxRefraction )
        ! where n is too big limit it and set derivatives to zero
        n = maxRefraction
        dNByDPtan = 0.0
        dNByDTempLower = 0.0
        dNByDTempUpper = 0.0
      elsewhere ! ( n <= maxRefraction )
        !{ Compute derivative of the refractive index w.r.t.\ PTan ($\zeta$).\\
        !  $\frac{\text{d} n}{\text{d} \zeta} =
        !   \frac{\text{d} n}{\text{d} P}
        !   \frac{\text{d} P}{\text{d} \zeta} =
        !   \left( \frac{n}P +
        !   \frac{\text{d} n}{\text{d} T} \frac{\text{d} T}{\text{d} P} +
        !   \frac{\text{d} n}{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}
        !   \frac{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}{\text{d} P} \right)
        !   \frac{\text{d} P}{\text{d} \zeta} =
        !   \frac{n}P \frac{\text{d} P}{\text{d} \zeta} +
        !   \frac{\text{d} n}{\text{d} T} \frac{\text{d} T}{\text{d} \zeta} +
        !   \frac{\text{d} n}{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}
        !   \frac{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}{\text{d} \zeta}$. \\
        !  Since $\frac{\text{d} P}{\text{d} \zeta} = - P \ln 10$, we have
        !  $\frac{\text{d} n}{\text{d} \zeta} = -n \ln 10 +
        !   \frac{\text{d} n}{\text{d} T} \frac{\text{d} T}{\text{d} \zeta} +
        !   \frac{\text{d} n}{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}
        !   \frac{\text{d} \mathbf{f}_{\text{H}_2\text{O}}}{\text{d} \zeta}$ 
        
        dnBydPTan = -n * ln10 + dnByDTemp * dTempByDPTan + dnByDH2O * dH2OByDPTan
        ! Compute derivatives w.r.t. temperature at lower and upper layers.
        dNByDTempLower = tempLowerWeight * dnByDTemp
        dNByDTempLower = tempUpperWeight * dnByDTemp
      end where

      ! Set derivatives w.r.t. pressure to zero for high pressure cases
    ! where ( 10**(-ptanVals) > maxPressure ) dNByDPtan=0.0
      where ( ptanVals < minZeta ) dNByDPtan=0.0

      ! Now from this calculate the refracted geocentric altitudes
      l1RefrGeocAlt= ( l1Alt%values(:,maf) + htOff ) / (1.0+n)

      ! And calculate their derivatives, again ignore H2O
      dl1RefrGeocAltByDHeightOffset = 1.0 / (1.0+n)
      dl1RefrGeocAltByDn = -l1RefrGeocAlt / (1.0+n)

      dl1RefrGeocAltByDPtan = dL1RefrGeocAltByDn * dNByDPtan

      dl1RefrGeocAltByDTempLower = dl1RefrGeocAltByDn * dNByDTempLower
      dl1RefrGeocAltByDTempUpper = dl1RefrGeocAltByDn * dNByDTempUpper

      ! Now we want to convert these to geopotentials
      ratio2=(earthRadA/l1RefrGeocAlt)**2
      ratio4=ratio2**2

      l1GPH = -((GM/g0)/l1RefrGeocAlt) * &
        & (1-J2*P2*ratio2-J4*P4*ratio4) - &
        & ((omega**2)/(2*g0)*(l1RefrGeocAlt*cos(geocLat))**2) + earthSurfaceGPH

      dL1GPHByDL1RefrGeocAlt= ( (GM/(l1RefrGeocAlt**2)) * &
        & (1-3*J2*P2*ratio2-5*J4*P4*ratio4) - &
        & omega**2*l1RefrGeocAlt*(cos(geocLat))**2 ) / g0

      ! Now get the combined derivatives for the l1Geopt's
      dL1GPHByDHeightOffset = dL1GPHByDL1RefrGeocAlt * &
        & dL1RefrGeocAltByDHeightOffset
      dL1GPHByDPtan = dL1GPHByDL1RefrGeocAlt * &
        & dL1RefrGeocAltByDPtan
      dL1GPHByDTempLower = dL1GPHByDL1RefrGeocAlt * &
        & dL1RefrGeocAltByDTempLower
      dL1GPHByDTempUpper = dL1GPHByDL1RefrGeocAlt * &
        & dL1RefrGeocAltByDTempUpper

      ! --------------------------------------------------------
      ! That's the end of the geometric geopotential calculation
      ! Now we're onto the hydrostatic calculation
      ! --------------------------------------------------------

      ! If this is the very first call, setup the hydrostatic temperature
      ! grid.
      if ( fmStat%newScanHydros ) then
        call allocate_test ( ifm%basisGPH, noTemps, &
          & temp%template%noInstances, 'ifm%gph', ModuleName )
        call allocate_test ( ifm%RT,  noTemps, &
          & temp%template%noInstances, 'ifm%gph', ModuleName )
        call allocate_test ( ifm%R, noTemps, &
          & 'ifm%gph', ModuleName )
        call GetBasisGPH ( temp, refGPH, ifm%basisGPH, &
          & ifm%R, ifm%RT, belowRef=ifm%belowRef )
        fmStat%newScanHydros = .false.
      end if

      rt => ifm%rt(:, closestTempProfiles(maf) )
      basisGPH => ifm%basisGPH (:, closestTempProfiles(maf) )

      ! Compute the basisGPH for each minor frame
      where (ptanVals>=tempBasis(1) .and. ptanVals<tempBasis(noTemps))
        hydrosGPH = basisGPH ( pointTempLayer ) + &
          & (rt(pointTempLayer)*tempUpperWeight*(2-tempUpperWeight) + &
          &  rt(pointTempLayer+1)*(tempUpperWeight**2))*(ln10*tempBasisGap/(2*g0))
      end where

      where ( ptanVals < tempBasis(1) )
        hydrosGPH = basisGPH ( pointTempLayer ) + &
          & rt(pointTempLayer) * (ln10/g0) * &
          & (ptanVals-tempBasis(pointTempLayer))
      end where

      where ( ptanVals > tempBasis(noTemps) )
        hydrosGPH = basisGPH (pointTempLayer+1) + &
          & rt(pointTempLayer+1) * (ln10/g0) * &
          & (ptanVals-tempBasis(pointTempLayer+1))
      end where

      ! Now we need dHydrosGPHByDT. To do this we calculate a matrix A containing
      ! the integrated basis functions.  There is an extra dimension in this
      ! to account for the reference surface which is then globally subtracted at
      ! the end.

      ! Unlike the code above, I think I'll do this in a loop to make life easier

      do i = 1, noMIFs + 1
        if ( i == noMIFs+1 ) then
          usePtan = refGPH%template%surfs(1,1)
          lower = ifm%belowRef
        else
          usePtan = ptanVals(i)
          lower = pointTempLayer(i)
        end if

        ! Now do the main bulk of the temperature bases
        do j = 1, noTemps

          ! Loop over temperature basis surfaces
          A(i,j) = 0.0_r8          ! Set to zero to start

          ! Now do nothing for bad ptans
          if ( usePtan /= ptan%template%badValue ) then
            basis = tempBasis(j)
            basisPlus =  tempBasis( min(j+1,noTemps))
            basisMinus = tempBasis( max(j-1,1) )

            ! Consider the extremes (j==1, j==noTemps) cases first

            ! Below bottom of basis
            if ( (j==1) .and. (usePtan < tempBasis(j)) ) then
              A(i,j) = 2*(usePtan-tempBasis(j))

              ! Above top
            else if ( (j==noTemps) .and. (usePtan > tempBasis(j)) ) then
              A(i,j) = 2*(usePtan-tempBasis(j))

              ! Basis triangles completely below ptan
            else if ( j < lower ) then
              A(i,j) = basisPlus - basisMinus

              ! Basis triangles with ptan in the top half (not when ptan above top)
            else if ( j == lower .and. usePtan < tempBasis(noTemps) ) then
              A(i,j) = ( basis - basisMinus ) + &
                & ( usePtan - basis ) * &
                & ( 2 * basisPlus - basis - usePtan ) / &
                & (basisPlus - basis)

              ! Basis triangles with ptan in the bottom half (not when ptan below bottom)
            else if ( j == lower+1 .and. usePtan >= tempBasis(1) ) then
              A(i,j) = ( ( usePtan - basisMinus )**2 ) / &
                & (basis - basisMinus)

            end if                         ! Basis triangles completely above ptan, no impact
          end if                           ! Good ptan
        end do                             ! Loop over temp basis
      end do                               ! Loop over ptans

      ! Now use this to get temperature derivatives
      do i = 1, noMIFs
        if ( ptanVals(i) /= ptan%template%badValue ) then
          do j = 1, noTemps
            dHydrosGPHByDTemp(i,j) = ( A(i,j) - A(noMIFs+1,j) ) * &
              & (ifm%R(j)*ln10) / (2*g0)
          end do
        end if
      end do

      ! Now get dHydrosGPHByDPtan
      do i= 1, noMIFs
        if ( ptanVals(i) /= ptan%template%badValue ) then
          dHydrosGPHByDPtan(i)= (ln10/g0) * ( &
            & rt(pointTempLayer(i)) * tempLowerWeight(i)+ &
            & rt(pointTempLayer(i)+1) * tempUpperWeight(i) )
        else
          dHydrosGPHByDPtan(i)=0.0
        end if
      end do

      ! -----------------------------------------------------------------------
      ! Now we calculate the residual, which is defined as the difference of the
      ! l1GPH and hydrosGPH. We also need it's derivatives wrt temp and ptan.
      ! First calculate residual. Also calculate the derivatives.
      ! -----------------------------------------------------------------------

      ! Compute residual (or differential residual)
      if ( fmConf%differentialScan ) then
        residual%values ( 1 : noMIFs-1, maf ) = &
          & ( l1GPH(2:noMIFs) - l1GPH(1:noMIFs-1) ) -&
          & ( hydrosGPH(2:noMIFs) - hydrosGPH(1:noMIFs-1) )
        residual%values( noMIFs, maf ) = 0.0
      else
        residual%values(:,maf) = l1GPH - hydrosGPH
      end if

      ! Find row in jacobian
      if ( present ( jacobian ) ) then
        row = FindBlock ( jacobian%row, residual%index, maf )
        fmStat%rows(row) = .true.

        ! Store the ptan derivatives
        if ( ptanInState ) then
          col = FindBlock ( jacobian%col, ptan%index, maf )
          block => jacobian%block(row,col)
          if ( block%kind /= M_Absent ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Found a prexisting d(residual)/d(ptan), removing' )
            call DestroyBlock ( block )
          end if
          if ( fmConf%differentialScan ) then
            call updateDiagonal ( block, [ &
              & ( dL1GPHByDPtan(2:noMIFs) - dL1GPHByDPtan(1:noMIFs-1) ) - &
              & ( dHydrosGPHByDPtan(2:noMIFs) - dHydrosGPHByDPtan(1:noMIFs-1) ) - &
              & 0.0 ] )
          else
            call updateDiagonal ( block, dL1GPHByDPtan-dHydrosGPHByDPtan )
          end if
        end if

        ! Store refGPH derivatives
        if ( refGPHInState ) then
          col = FindBlock ( jacobian%col, refGPH%index, maf )
          block => jacobian%block(row,col)
          if ( fmConf%differentialScan ) then
            call DestroyBlock ( block )
          else
            if ( block%kind /= M_Absent ) then
              call MLSMessage ( MLSMSG_Warning, ModuleName, &
                & 'Found a prexisting d(residual)/d(refGPH), removing' )
              call DestroyBlock ( block )
            end if
            call CreateBlock ( jacobian, row, col, M_Full )
            block%values = -1.0
          end if
        end if

        ! Store heightOffset derivatives if any
        if ( associated (heightOffset) .and. heightOffsetInState ) then
          col = FindBlock ( jacobian%col, heightOffset%index, maf )
          block => jacobian%block(row,col)
          if ( block%kind /= M_Absent ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Found a prexisting d(residual)/d(heightOffset), removing' )
            call DestroyBlock ( block )
          end if
          call CreateBlock ( jacobian, row, col, M_Full )
          if  ( fmConf%differentialScan ) then
            block%values(:,1) = [ &
              & dL1GPHByDHeightOffset(2:noMIFs) - dL1GPHByDHeightOffset(1:noMIFs-1), &
              & 0.0_r8 ]
          else
            block%values(:,1) = dL1GPHByDHeightOffset
          end if
        end if

        ! Now the temperature derivatives
        if ( tempInState ) then
          col = FindBlock ( jacobian%col, temp%index, maf )
          block => jacobian%block(row,col)
          if ( block%kind /= M_Absent ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Found a prexisting d(residual)/d(temp), removing' )
            call DestroyBlock ( block )
          end if
          call CreateBlock ( jacobian, row, col, M_Full )
          if ( fmConf%differentialScan ) then
            ! ------------- Diffential model
            block%values(noMIFs,:) = 0.0
            block%values(1:noMIFs-1,:) = - ( dHydrosGPHByDTemp(2:noMIFs,:) - &
              & dHydrosGPHByDTemp(1:noMIFs-1,:) )
            do i = 2, noMIFs
              block%values(i-1,pointTempLayer(i)) = &
                & block%values(i-1,pointTempLayer(i)) + dL1GPHByDTempLower(i)
              block%values(i-1,pointTempLayer(i)+1) = &
                & block%values(i-1,pointTempLayer(i)+1) + dL1GPHByDTempUpper(i)
            end do
            do i = 1, noMIFs-1
              block%values(i,pointTempLayer(i)) = &
                & block%values(i,pointTempLayer(i)) - dL1GPHByDTempLower(i)
              block%values(i,pointTempLayer(i)+1) = &
                & block%values(i,pointTempLayer(i)+1) - dL1GPHByDTempUpper(i)
            end do
          else
            ! ------------- Non differential model
            block%values = - dHydrosGPHByDTemp
            do i = 1, noMIFs
              block%values(i,pointTempLayer(i)) = &
                & block%values(i,pointTempLayer(i)) + dL1GPHByDTempLower(i)
              block%values(i,pointTempLayer(i)+1) = &
                & block%values(i,pointTempLayer(i)+1) + dL1GPHByDTempUpper(i)
            end do
          end if
        end if

      end if

    end subroutine ScanForwardModelAuto

  end subroutine ScanForwardModel
  ! ---------------------------------------  TwodScanForwardModel  -----
  subroutine TwoDScanForwardModel ( FmConf, State, Extra, FwmOut, &
                                  & FmStat, Jacobian, ChunkNo )

use dump_0, only: dump
    use Array_Stuff, only: Subscripts
    use Check_QTM_m, only: Check_QTM
    use Get_Species_Data_M, only:  Get_Species_Data
    use HGridsDatabase, only: HGrid_T
    use MatrixModule_0, only: MatrixElement_t
    use Output_M, only: Blanks, Output
    use Physics, only: Boltz
    use Piq_Int_M, only: Piq_Int
    use Sparse_Eta_m, only: Sparse_Eta_t
    use Tangent_Qty_m, only: Tangent_Qty
    use Time_M, only: Time_Now

  ! This is a two D version of ScanForwardModel
  ! inputs
    type (ForwardModelConfig_T), intent(inout) :: FmConf ! Configuration options
    type (Vector_T), intent(in) :: State ! The state vector
    type (Vector_T), intent(in) :: Extra ! Other stuff in the state vector
  ! outputs
    type (Vector_T), intent(inout) :: FwmOut ! Output vector, residual filled
    type (ForwardModelStatus_T), intent(inout) :: FmStat ! Which maf etc.
    type (Matrix_T), intent(inout), optional :: Jacobian ! The derivative matrix
    integer, intent(in), optional :: ChunkNo

    ! To make an array of pointers to MatrixElement_t blocks
    type :: Blocks_t
      type(matrixElement_t), pointer :: B
    end type Blocks_t

  ! local variables
    type (HGrid_t), pointer :: QTM_HGrid        ! QTM with finest resolution
    type (VectorValue_T), pointer :: OrbIncline ! Orbital inclination
    type (VectorValue_T), pointer :: Temp       ! Temperature component of state
    type (VectorValue_T), pointer :: PTAN       ! Ptan component of state
    type (VectorValue_T), pointer :: PhiTAN     ! Phitan component of state
    type (VectorValue_T), pointer :: H2O        ! H2O component of state
    type (VectorValue_T), pointer :: RefGPH     ! Ref gph component of state
    type (VectorValue_T), pointer :: L1Alt      ! Tangent point altitude
    type (VectorValue_T), pointer :: Residual   ! Resulting component of fwmOut
    type (MatrixElement_T), pointer :: Block    ! A matrix block

    logical :: TempInState           ! Set if Temp in state not extra
    logical :: RefGPHInState         ! Set if RefGPH in state not extra
    logical :: H2OInState            ! Set if H2O in state, not extra
    logical :: PTANInState           ! Set if PTan in state, not extra
    logical :: UsingQTM              ! Set if Temp and all species are QTM

    integer :: E                     ! Index of Sparse_t E component
    integer :: Me = -1               ! String index for trace
    integer :: R                     ! Row index from Sparse_t E(e) component
    integer :: Row                   ! Block row in jacobian
    integer :: Col                   ! Block col in jacobian
    integer :: Sv_Z                  ! height basis index
    integer :: Sv_P                  ! horizontal basis index
    integer :: Windowstart_t         ! first instance for temperature
    integer :: Windowfinish_t        ! last instance for temperature
    integer :: Windowstart_h2o       ! first instance for water vapor
    integer :: Windowfinish_h2o      ! last instance for water vapor

    real(rp) :: z_surf
    real(rp) :: surf_temp
    real(rp) :: surf_refr_indx(1)

    type(sparse_eta_t) :: Eta_P_H2O     ! H2O horizontal
    type(sparse_eta_t) :: Eta_P_T       ! Temperature horizontal
    type(sparse_eta_t) :: Eta_ZXP_T     ! Temperature zeta X horizontal
    type(sparse_eta_t) :: Eta_ZXP_H2O   ! Water zeta X horizontal

    real :: T0, T1, T2                     ! For timing
    logical, parameter :: always_timing = .false.  ! if worried about NAG taking so long
    logical :: Timing  ! if worried about NAG taking so long
    logical, parameter :: Total_Times = .true.

integer :: Times = 0
times = times + 1
    call trace_begin ( me, 'TwoDScanForwardModel, MAF=', index=fmstat%maf, &
      & cond=toggle(emit) ) ! set by -f command-line switch
  ! Time this
    call time_now ( t0 )
    call time_now ( t1 )
    Timing = always_timing
  !  if ( present(chunkNo) ) &
  !    & Timing = always_timing .or. &
  !    & chunkNo == 4 .or. chunkNo == 12 .or. chunkNo == 25 .or.  &
  !    & chunkNo == 50 .or. chunkNo == 167
    if ( timing ) &
      & call output('beginning timing for 2d scan forward model', advance='yes')
  ! Identify the vector quantities from state/extra
    orbIncline => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_orbitInclination, config=fmConf )
    temp => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_temperature, config=fmConf, foundInFirst=tempInState )
    ptan => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_ptan, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf, foundInFirst=ptanInState )
    phitan => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_phitan, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf)
    h2o => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_vmr, molecule=l_h2o, &
      & config=fmConf, foundInFirst=h2oInState, noError=.true.)
    refgph => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_refGPH, config=fmConf, foundInFirst=refGPHInState )
    l1Alt => getQuantityForForwardModel ( state, extra, &
      & quantityType=l_tngtGeocAlt, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf)
  ! Identify the vector quantities from fwmOut
    residual => getQuantityForForwardModel ( fwmOut, &
      & quantityType=l_scanResidual, instrumentModule=fmConf%instrumentModule, &
      & config=fmConf)
    if ( timing ) call sayTime ( 'Getting vector quantities' )
    call time_now ( t1 )

    ! Get pointers to the temperature and other quantities from the state
    ! vector into the configuration. This has to be done AFTER
    ! DeriveFromForwardModelConfig, which is invoked from ForwardModelWrappers.

    call get_species_data ( fmConf, state, extra )

    ! Check whether temperature and all the mixing ratios either all have QTM
    ! HGrid, or none do.  If so, find the QTM HGrid with the finest resolution.
    ! (We currently require that if they're QTM they all have the same grid.)

    call check_QTM ( fmConf, QTM_HGrid, usingQTM )

    ! Get windows
    if ( usingQTM ) then
      windowstart_t = 1
      windowfinish_t = size(temp%template%the_hGrid%QTM_tree%geo_in)
      windowstart_h2o = 1
      windowfinish_h2o = size(h2o%template%the_hGrid%QTM_tree%geo_in)
    else
      call FindInstanceWindow ( temp, phitan, fmStat%maf, FMConf%phiWindow, &
        & FMConf%windowUnits, windowstart_t, windowfinish_t )
      call FindInstanceWindow( h2o, phitan, fmStat%maf, FMConf%phiWindow, &
        & FMConf%windowUnits, windowstart_h2o, windowfinish_h2o )
    end if

    if ( timing ) call sayTime ( 'Finding instance windows' )
    call time_now ( t1 )

    block ! So that local arrays are automatic instead of allocatable
      type(blocks_t) :: Blocks(windowstart_t:windowfinish_t)
      real(rp) :: Coslat2(ptan%template%nosurfs)
      real(rp) :: dGPHdR(ptan%template%nosurfs)
      real(rp) :: Earth_radius(ptan%template%nosurfs)
      real(rp) :: Eff_earth_radius(ptan%template%nosurfs)
      real(rp) :: Eta_piqxp( ptan%template%nosurfs, temp%template%nosurfs, &
                           & windowfinish_t - windowstart_t + 1 )
      real(rp) :: G_ref(ptan%template%nosurfs)
      real(rp) :: L1AltRefr(ptan%template%nosurfs)
      real(rp) :: L1RefAlt(ptan%template%nosurfs)
      real(rp) :: Mass_corr(ptan%template%nosurfs)
      real(rp) :: P2(ptan%template%nosurfs)
      real(rp) :: P4(ptan%template%nosurfs)
      real(rp) :: piq(ptan%template%nosurfs,temp%template%nosurfs)
      real(rp) :: Ratio2_gph(ptan%template%nosurfs)
      real(rp) :: Ratio2(ptan%template%nosurfs)
      real(rp) :: Ratio4_gph(ptan%template%nosurfs)
      real(rp) :: Ratio4(ptan%template%nosurfs)
      real(rp) :: RefGeomAlt_denom(ptan%template%nosurfs)
      real(rp) :: Refgph_Surfs(ptan%template%nosurfs)
      real(rp) :: Tan_H2O(ptan%template%nosurfs)
      real(rp) :: Tan_Refr_Indx(ptan%template%nosurfs)
      real(rp) :: Tan_temp(ptan%template%nosurfs)
      real(rp) :: Work(ptan%template%nosurfs)

      ! Geometry calculations
      block ! So that some even-more local arrays are automatic
        real(rp) :: earthradc(ptan%template%nosurfs)  ! square of minor axis of
    !                       earth ellipsoid in orbit plane projected system
        real(rp) :: red_phi_t(ptan%template%nosurfs)
        real(rp) :: sinbeta(ptan%template%nosurfs)
        real(rp) :: sinphi2(ptan%template%nosurfs)
        real(rp) :: cosphi2(ptan%template%nosurfs)
        real(rp) :: geoclats(ptan%template%nosurfs)
        real(rp) :: sinlat2(ptan%template%nosurfs)

        ! Sometimes, orbincline may not have the same number of surfaces as
        ! ptan, e.g. for non-satellite data
        if ( size(orbincline%values, 1) < ptan%template%noSurfs ) then
          sinbeta = sin(deg2rad*orbincline%values(1,fmStat%maf))
        else
          sinbeta = sin(deg2rad*orbincline%values(1:ptan%template%noSurfs,fmStat%maf))
        end if
        earthradc = (earthrada*earthradb)**2 / &
          & ( (earthRada**2 - earthRadb**2) * sinbeta**2 + earthRadb**2) ! in meters
      ! rephase the phi
        red_phi_t = modulo(deg2rad*phitan%values(:,fmStat%maf),2.0_rp*Pi)
        where ( 0.5_rp*Pi < red_phi_t .AND. red_phi_t <= 1.5_rp*Pi )
          red_phi_t = Pi - red_phi_t
        elsewhere ( red_phi_t > 1.5_rp*Pi )
          red_phi_t = red_phi_t - 2.0_rp*Pi
        end where
        if ( timing ) call sayTime ( 'rephasing phi' )
        call time_now ( t1 )
      ! compute sin^2(phi) and cos^2(phi)
        sinphi2 = sin(red_phi_t)**2
        cosphi2 = 1.0_rp - sinphi2
      ! convert phitan into geocentric latitude
        geoclats = asin(earthradc * sin(red_phi_t) * sinbeta / &
          & sqrt(earthrada**4*cosphi2 + earthradc**2*sinphi2))
        sinlat2 = sin(geoclats)**2
        coslat2 = 1.0_rp - sinlat2
        p2=0.5_rp * (3.0_rp*sinlat2 - 1.0_rp)
        p4=0.125_rp * (35.0_rp*sinlat2**2 - 30.0_rp*sinlat2 + 3.0_rp)
      ! compute the local gravitational acceleration at the surface
        earth_radius = sqrt((earthrada**4*cosphi2 + earthradc**2*sinphi2) &
                     &    / (earthrada**2*cosphi2 + earthradc*sinphi2))
        ratio2=(earthRadA/earth_radius)**2
        ratio4=ratio2**2
        g_ref =  GM * (1.0_rp - 3.0_rp*j2*p2*ratio2 - 5.0_rp*j4*p4*ratio4) &
          & / earth_radius**2 - omega**2*earth_radius*coslat2
      ! get the effective earth radius
        eff_earth_radius = 2.0_rp * g_ref / (2.0_rp * gm * (1.0_rp-6.0_rp*j2*p2 &
          & * ratio2 - 15.0_rp*j4*p4*ratio4) / earth_radius**3 &
          & + omega**2*coslat2)
      end block

      ! Get Eta functions, and H2O and temperature tangent profiles.
      call tangent_qty ( fmConf, usingQTM, fmStat%MAF, ptan, phiTan, &
                       & h2o, windowStart_h2o, windowFinish_h2o, tan_h2o, &
                       & eta_zxp_h2o, eta_p_h2o, logLin=.true. )
      call tangent_qty ( fmConf, usingQTM, fmStat%MAF, ptan, phiTan, &
                       & temp, windowStart_t, windowFinish_t, tan_temp, &
                       & eta_zxp_t, eta_p_t )
    ! compute refractive index at the MIF tangent points
      call refractive_index ( 10.0_rp**(-max(-10.0_rp,ptan%values(:,fmStat%maf))), &
                            & tan_temp, tan_refr_indx, h2o_path = tan_h2o )
      if ( timing ) call sayTime ( 'getting etas etc.' )
      call time_now ( t1 )

    ! construct piq integral
      call piq_int ( ptan%values(:,fmStat%maf), temp%template%surfs(:,1), &
        & refGPH%template%surfs(1,1), piq, Z_MASS = 2.5_rp, C_MASS = 0.02_rp)
      if ( timing ) call sayTime ( 'constructing piq integral' )
      call time_now ( t1 )

    ! Convert level 1 reference geopotential height into geometric altitude.
    ! This assumes that the reference ellipsoid is equivalent to a reference
    ! geopotential height of 0 meters.
    ! Assume RefGPH and Temperature have the same horizontal basis.
      call eta_p_t%sparse_dot_vec ( &
        & refgph%values(1,windowstart_t:windowfinish_t), refgph_surfs )
      refgeomalt_denom = g_ref*eff_earth_radius - g0*refgph_surfs
      l1refalt = g_ref*eff_earth_radius**2 / refgeomalt_denom - eff_earth_radius &
        & + earth_radius
      eta_piqxp = 0
      do e = 1, eta_p_t%ne
        r = eta_p_t%e(e)%r
        eta_piqxp(r,:,eta_p_t%e(e)%c) = piq(r,:) * eta_p_t%e(e)%v
      end do
      if ( timing ) call sayTime ( 'constructing geometric altitude' )
      call time_now ( t1 )

      ! Compute surface pressure. Compute temperature at surface phi.
      ! Compute H2O at surface phi and surface pressure.
      block
        type(sparse_eta_t) :: Eta_At_One_Zeta
        type(sparse_eta_t) :: Eta_At_Z_Surf
        real(rp) :: Temp_At_Surf_Phi(temp%template%nosurfs)
        real(rp) :: Surf_H2O(h2o%template%nosurfs)

        call eta_p_t%row_dot_matrix ( 1, &
          & temp%values(:,windowstart_t:windowfinish_t), temp_at_surf_phi )
        z_surf = z_surface( temp%template%surfs(:,1),temp_at_surf_phi, g_ref(1), &
                          & refgph_surfs(1), refGPH%template%surfs(1,1), boltz )
      ! compute temperature at Z_Surf
        call eta_at_z_surf%eta_1d ( temp%template%surfs(:,1), [z_surf] )
        surf_temp = eta_at_z_surf%row_dot_vec ( 1, temp_at_surf_phi )
        ! Compute H2O at Z_Surf
        ! Logarithmic interpolation:
        ! First surf_h2o = eta_p_h2o .dot. log(max(h2o%values,1.0e-9_rp)),
        ! then h2o_path = exp(eta_at_z_surf .dot. surf_h2o)
        ! surf_h2o = 0 ! not necessary because most elements won't be used
        if ( eta_p_h2o%rows(1) /= 0 ) then
          call eta_at_z_surf%eta_1d ( h2o%template%surfs(:,1), [z_surf], &
            & empty=.true. )
          ! Surf_H2O is log(H2O at Z_Surf).  It will be an operand in a dot
          ! product with Eta_At_Z_Surf, so we only need it where Eta_At_Z_Surf
          ! has a nonzero column -- and there are only two of those.
          ! If TangentQuantity is changed to send the basis boundaries to
          ! Eta_p%Eta_1D, DO NOT put windowStart_h2o:windowFinish_h2o as a
          ! section subscript on the second dimension of H2O%Values!
          call eta_p_h2o%row_dot_matrix ( 1, &
            & h2o%values(:,windowStart_h2o:windowFinish_h2o), surf_h2o, &
            & logOnly = .true., rows=eta_at_z_surf )
          ! Compute refractive index at Z_Surf
          call refractive_index(10.0_rp**([ -z_surf ]), [ surf_temp ], &
            & surf_refr_indx, &
            & h2o_path = [ exp(eta_at_z_surf%row_dot_vec(1,surf_h2o)) ] )
        else
          ! Compute refractive index at Z_Surf
          call refractive_index(10.0_rp**([ -z_surf ]), [ surf_temp ], &
            & surf_refr_indx )
        end if
      end block
    ! set all refr indicies below the surface to the surface value
      where(ptan%values(:,fmStat%maf) < z_surf) 
        tan_refr_indx = surf_refr_indx(1)
        tan_temp = surf_temp
      end where
    ! compute l1 geopotential
      l1altrefr = l1alt%values(:,fmStat%maf) / (1.0_rp + tan_refr_indx)
      ratio2 = (earthRadA/l1altrefr)**2
      ratio4 = ratio2**2
      ratio2_gph = (earthRadA/l1refalt)**2
!????? Bill says ratio4_gph ought to be ratio2_gph**2 ?????
      ratio4_gph = ratio2**2
    ! do a simple mass correction
      mass_corr = 1.0_rp
    !  where(ptan%values(:,fmStat%maf) > 2.5) mass_corr = 1.0_rp &
    !    & / (0.875_rp + 0.1_rp*ptan%values(:,fmStat%maf) &
    !    & - 0.02_rp*ptan%values(:,fmStat%maf)**2)
    ! This is a reasonable approximation to the above.
      where(ptan%values(:,fmStat%maf) > 2.5) mass_corr = 1.0_rp &
         & + 0.02_rp*(ptan%values(:,fmStat%maf) - 2.5)**2
      if ( timing ) call sayTime ( 'refractive index, l1 geopotential' )
      call time_now ( t1 )
    ! forward model calculation
      do row = 1, ptan%template%nosurfs
        work(row) = sum ( temp%values(:,windowstart_t:windowfinish_t) * &
                        & eta_piqxp(row,:,:) ) * (boltz/g0)
      end do
      residual%values(:,fmStat%maf) = &
      & ( GM/g0) * ((1.0_rp - j2*p2*ratio2 - j4*p4*ratio4)  / l1altrefr &
      &             - (1.0_rp - j2*p2*ratio2_gph - j4*p4*ratio4_gph) / &
      &               l1refalt) &
      & + (omega**2/(2*g0)) * coslat2 * (l1altrefr - l1refalt) * &
      &   (l1altrefr + l1refalt) &
      & + work
      if ( timing ) call sayTime ( 'forward model calculation' )
      call time_now ( t1 )
      if ( fmConf%differentialScan ) residual%values(:,fmStat%maf) = &
        & EOSHIFT(residual%values(:,fmStat%maf), 1, &
        & residual%values(ptan%template%nosurfs,fmStat%maf)) &
        & - residual%values(:,fmStat%maf)
    ! derivatives--for simplicity we are ignoring h2o contributions
      if ( present ( jacobian ) ) then
        ! Find row in jacobian
        row = FindBlock ( jacobian%row, residual%index, fmStat%maf )
        fmStat%rows(row) = .true.
    ! compute dudr
        dgphdr =  -( (GM/g0) / l1altrefr**2) * (1.0_rp - 3.0_rp*j2*p2*ratio2 &
        & - 5.0_rp*j4*p4*ratio4) + (omega**2/g0)*l1altrefr*coslat2
    ! Store the ptan derivatives
        if ( ptanInState ) then
          col = FindBlock ( jacobian%col, ptan%index, fmStat%maf )
          block => jacobian%block(row,col)
          if ( block%kind /= M_Absent ) then
            call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Found a prexisting d(residual)/d(ptan), removing' )
            call DestroyBlock ( block )
          end if
          block
            real(rp) :: dScandz(ptan%template%nosurfs)
            !{ $\frac{\text{d Scan}}{\text{d} z} = 
            !  \frac{k}{g_0} \, \eta_{\zeta \times p}^T \cdot T \,\times\, $
            !  mass\_corr
            call eta_zxp_t%sparse_dot_vec ( &
              & temp%values(:,windowstart_t:windowfinish_t), dScandz )
            dscandz = ( boltz / g0 ) * mass_corr * dscandz
            where ( ptan%values(:,fmStat%maf) > z_surf ) dscandz = dscandz + &
              & dgphdr*l1altrefr*tan_refr_indx*ln10 / (1.0_rp + tan_refr_indx)
              ! UpdateDiagonal creates a banded block
            if ( fmConf%differentialScan ) then
              call updateDiagonal ( block, EOSHIFT(dscandz, 1, &
              & dscandz(ptan%template%nosurfs)) - dscandz)
            else
              call updateDiagonal ( block, dscandz )
            end if
! call dump ( dscandz, 'PTan derivatives' )
          end block
        end if
    ! Store refGPH derivatives
        if ( refGPHInState ) then
          ! Find/create necessary Jacobian blocks
          do sv_p = windowstart_t, windowfinish_t
            col = FindBlock ( jacobian%col, refGPH%index, sv_p )
            blocks(sv_p)%b => jacobian%block(row,col)
            if ( blocks(sv_p)%b%kind /= M_Absent ) then
              call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'Found a prexisting d(residual)/d(refGPH), removing' )
              call DestroyBlock ( blocks(sv_p)%b )
            end if
            if ( any(eta_p_t%e%r == sv_p ) ) &
              call CreateBlock ( jacobian, row, col, M_Full, init=0.0_rm )
          end do
          work = (GM * (1.0_rp - 3.0_rp*j2*p2*ratio2_gph &
               & - 5.0_rp*j4*p4*ratio4_gph) / l1refalt**2 &
               & + omega**2*l1refalt*coslat2) * &
               & (l1refalt+eff_earth_radius - earth_radius) / &
               & refgeomalt_denom
          !{ The dot product $\eta_p^T \cdot $ work is done ``the hard way''
          !  because the column index in eta_p_t\%e(e)\%c needs to be offset
          !  by WindowStart_t, and it gives the block index, not a subscript.
          do e = 1, eta_p_t%ne
            r = eta_p_t%e(e)%r
            sv_p = eta_p_t%e(e)%c + windowStart_t - 1
            blocks(sv_p)%b%values(r,1) =  blocks(sv_p)%b%values(r,1) + &
              & work(r) * eta_p_t%e(e)%v
          end do
          if ( fmConf%differentialScan ) then
    ! ------------- Differential model
            do sv_p = windowStart_t, windowFinish_t
              blocks(sv_p)%b%values = EOSHIFT(blocks(sv_p)%b%values, &
                & SPREAD(1, 1, ptan%template%nosurfs), &
                & RESHAPE(blocks(sv_p)%b%values(ptan%template%nosurfs,:), &
                & [ temp%template%nosurfs ]), dim=2) - blocks(sv_p)%b%values
            end do
          end if
! do sv_p = windowstart_t, windowfinish_t
! print '(a,i0)', 'SV_P ', sv_p
! call dump ( blocks(sv_p)%b%values, 'refGPH derivatives' )
! end do
          if ( timing ) call sayTime ( 'differential model' )
          call time_now ( t1 )
        end if
    ! Now the temperature derivatives
        if ( tempInState ) then
          block
            integer :: C     ! Column of Eta_ZXP_T element
            integer :: C2(2) ! Temperature's Zeta, Phi subscripts,
                             ! Sub-column subscripts of PIQXP, second and third
                             ! subscripts of PIQXP3
            ! Find/create necessary Jacobian blocks
            do sv_p = windowstart_t, windowfinish_t
              col = FindBlock ( jacobian%col, temp%index, sv_p )
              blocks(sv_p)%b => jacobian%block(row,col)
              if ( blocks(sv_p)%b%kind /= M_Absent ) then
                call MLSMessage ( MLSMSG_Warning, ModuleName, &
                  & 'Found a prexisting d(residual)/d(temp), removing' )
                call DestroyBlock ( blocks(sv_p)%b )
              end if
              if ( any(eta_zxp_t%e%r == sv_p) ) &
                call CreateBlock ( jacobian, row, col, M_Full, init=0.0_rm )
              blocks(sv_p)%b%values = (boltz/g0) * eta_piqxp(:,:,sv_p-windowstart_t+1)
            end do
            work = dgphdr * l1altrefr * tan_refr_indx / &
                 & ((1.0_rp + tan_refr_indx) * tan_temp)
            !{ The dot product $\eta_{\zeta \times p}^T \cdot$ work is done
            !  ``the hard way'' because the column index in eta_p_t\%e(e)\%c is
            !  actually the array-order index of a two-dimensional array, one of
            !  which needs to be offset by WindowStart_t, and gives the block
            !  index, not a subscript.
            do e = 1, eta_zxp_t%ne
              r = eta_zxp_t%e(e)%r
              c = eta_zxp_t%e(e)%c
              c2 = subscripts( c, [ temp%template%nosurfs, windowfinish_t ], &
                             &    [ 1, windowstart_t ] )
              blocks(c2(2))%b%values(r,c2(1)) = blocks(c2(2))%b%values(r,c2(1)) + &
                & work(r) * eta_zxp_t%e(e)%v
            end do
            if ( fmConf%differentialScan ) then
    ! ------------- Differential model
              do sv_p = windowStart_t, windowFinish_t
                blocks(sv_p)%b%values = EOSHIFT(blocks(sv_p)%b%values, &
                & SPREAD(1, 1, ptan%template%nosurfs), &
                & RESHAPE(blocks(sv_p)%b%values(ptan%template%nosurfs,:), &
                & [ temp%template%nosurfs ]), dim=2) - blocks(sv_p)%b%values
              end do
            end if
          end block
          if ( timing ) call sayTime ( 'temperature derivatives' )
          call time_now ( t1 )
        end if
      end if
    end block

    if ( timing ) call sayTime ( 'deallocating the rest' )
    if ( timing ) then
      t1 = t0
      call sayTime ( 'all of 2d scan forward model' )
    end if

    call trace_end ( 'TwoDScanForwardModel MAF=', fmStat%maf, &
      & cond=toggle(emit) )

  contains
    ! ..................................................  SayTime  .....
    subroutine SayTime ( What )
      character(len=*), intent(in) :: What
      call time_now ( t2 )
      if ( total_times ) then
        call output ( "Total time = " )
        call output ( dble(t2), advance = 'no' )
        call blanks ( 4, advance = 'no' )
      end if
      call output ( "Timing for " // what // " = " )
      call output ( dble(t2 - t1), advance = 'yes' )
    end subroutine SayTime

    ! ................................................  Z_surface  .....
    real(rp) function Z_surface ( z_basis, t_values, g_ref, h_ref, z_ref, boltz, &
      & threshold, maxiterations )
! finds the surface pressure
      use Sparse_Eta_m, only: Sparse_Eta_t
! inputs:
      real(rp), intent(in) :: Z_basis(:) ! zeta basis for temperature
      real(rp), intent(in) :: T_values(:) ! temperature coefficients
      real(rp), intent(in) :: G_ref      ! reference acceloration at the surface
      real(rp), intent(in) :: H_ref      ! reference geopotential
      real(rp), intent(in) :: Z_ref      ! reference zeta for h_ref
      real(rp), intent(in) :: Boltz      ! k*ln10/m in your favorite units
! note that the units of h_ref, g_ref, and boltz must be consistent
      real(rp), optional, intent(in) :: Threshold ! zeta convergence criteria
      integer, optional, intent(in) :: MaxIterations ! maximum number of
!                                  iterations
! This is a one dimensional calculation where g_ref is assumed to apply
! throughout the entire vertical range of t_basis.
! internals
      integer :: Iter
      integer :: MaxIter

      type(sparse_eta_t) :: Eta
      real(rp) :: GHB                    ! g_ref * h_ref / boltz
      real(rp) :: PIQ(1,size(z_basis))
      real(rp) :: Thresh
      real(rp) :: Z_old(1)

! begin code
      maxiter = 10
      thresh = 0.0001_rp
      if ( present(maxiterations) ) maxiter = maxiterations
      if ( present(threshold) ) thresh = threshold

! intitial guess
      ghb = g_ref * h_ref / boltz
      z_old = z_ref - ghb / t_values(1)

      iter = 0
      do
        call piq_int ( z_old, z_basis, z_ref, piq )
        call eta%eta ( z_basis, z_old, empty=.true. )
! correct
        z_surface = z_old(1) - (ghb + dot_product(piq(1,:),t_values)) / &
                             & eta%row_dot_vec(1,t_values)

        iter = iter + 1
        if ( abs(z_surface - z_old(1)) < thresh .or. iter >= maxiter ) exit
        z_old = z_surface
      end do
      if (Timing) call output( 'Num iterations in z_surface: ', advance='no')
      if (Timing) call output( iter, advance='yes')
    end function Z_surface

  end subroutine TwoDScanForwardModel

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ScanModelModule

! $Log$
! Revision 2.95  2020/05/15 21:33:41  pwagner
! Fixed long-standing in sending basis boundaries
!
! Revision 2.94  2018/12/13 01:45:28  vsnyder
! More QTM work -- hopefully the last of it
!
! Revision 2.93  2018/12/04 23:22:33  vsnyder
! Add QTM horizontal interpolation.  Spell outputNamedValue correctly.
!
! Revision 2.92  2018/11/30 00:13:28  pwagner
! Resets gph to Fill Values where Temperatures do, too
!
! Revision 2.91  2018/11/01 00:43:21  vsnyder
! Make sure there are no uninitialized elements of Jacobian blocks
!
! Revision 2.90  2018/10/30 20:59:14  vsnyder
! Completely revised.  Use sparse interpolators.  Convert allocatable
! variables to automatic variables.
!
! Revision 2.89  2017/08/10 22:44:12  pwagner
! Exit witth eror message if buggy 1s Scan Model called

! Revision 2.88  2016/08/12 00:31:25  pwagner
! Seems to restore tthe gold brick

! Revision 2.87  2016/08/09 18:46:16  pwagner
! Survives encounter with non-satellite data

! Revision 2.86  2016/01/23 02:59:18  vsnyder
! Get MexRefraction from refraction_m, not geometry.  Get refractive index
! derivatives from refraction_m instead of computing them here.  Add LaTeX.

! Revision 2.85  2015/09/24 22:06:36  pwagner
! Code around segment fault due to NAG bug

! Revision 2.84  2015/08/25 17:35:56  vsnyder
! PhiWindow is a tuple, with the first element specifying the angles or
! number of profiles/MAFs before the tangent point, and the second
! specifying the angles or number after.

! Revision 2.83  2015/06/19 21:15:50  pwagner
! Maneuvers around calculation of ITSS; buggy NAG sometimes segment faults here

! Revision 2.82  2015/06/04 03:14:14  vsnyder
! Make Surfs component of quantity template allocatable

! Revision 2.81  2014/07/18 23:17:59  pwagner
! Aimed for consistency in names passed to allocate_test

! Revision 2.80  2014/01/09 00:30:24  pwagner
! Some procedures formerly in output_m now got from highOutput

! Revision 2.79  2013/08/31 02:29:12  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end

! Revision 2.78  2013/08/30 02:45:47  vsnyder
! Revise calls to trace_begin and trace_end

! Revision 2.77  2013/06/12 02:38:33  vsnyder
! Cruft removal

! Revision 2.76  2013/03/01 01:11:10  pwagner
! Added DumpInstanceWindows

! Revision 2.75  2012/06/06 20:13:37  vsnyder
! Avoid overflow if ptan is bogus because module is turned off

! Revision 2.74  2012/04/20 01:54:20  vsnyder
! Remove GeocLat from Get2DHydrostaticTangentPressure, make EarthRadius
! automatic.  Add some dumps.  Print a warning if the initial pressure
! guess is <-4 or >10.  Probably still gets overflows if the module is turned
! off.  More tracing.

! Revision 2.73  2011/05/09 18:25:12  pwagner
! Converted to using switchDetail

! Revision 2.72  2009/06/23 18:46:18  pwagner
! Prevent Intel from optimizing ident string away

! Revision 2.71  2009/05/13 20:41:55  vsnyder
! Get constants from Constants, kinds from MLSKinds

! Revision 2.70  2007/08/20 22:06:08  pwagner
! Two procedures now push their names onto MLSCallStack

! Revision 2.69  2007/07/25 20:08:46  vsnyder
! Delete declaration for unused variable

! Revision 2.68  2007/06/29 19:32:07  vsnyder
! Make ForwardModelIntermediate_t private to ScanModelModule

! Revision 2.67  2006/12/13 01:33:08  vsnyder
! Use slightly faster get_eta_sparse

! Revision 2.66  2006/08/05 02:36:06  vsnyder
! Delete unused symbols

! Revision 2.65  2006/08/05 02:12:27  vsnyder
! Add ForWhom argument to ConstructVectorTemplate

! Revision 2.64  2006/08/02 20:10:56  vsnyder
! Use deallocate_test for myR and myRT, for leak tracking

! Revision 2.63  2006/04/03 20:23:46  pwagner
! Preserve correct sign for GPH Precision

! Revision 2.62  2005/06/22 18:57:02  pwagner
! Reworded Copyright statement, moved rcs id

! Revision 2.61  2004/03/20 04:05:23  vsnyder
! Moved Boltz from units to physics

! Revision 2.60  2003/08/15 23:58:20  vsnyder
! Get PHYQ_... directly from Intrinsic instead of indirectly via Units

! Revision 2.59  2003/03/19 19:24:10  pwagner
! Prints timing info when necessary; speedup of TwoDScanForwardModel despite terrible NAG memory management

! Revision 2.58  2003/02/18 23:59:28  livesey
! Added phiWindow stuff for hydrostatic fill

! Revision 2.57  2003/02/13 01:15:39  bill
! made reference geopotential height calculation constants the same as those used in the hydrostatic calcs

! Revision 2.56  2003/02/12 22:46:49  bill
! fixed mass correction integrations

! Revision 2.55  2003/01/31 22:30:56  livesey
! Got rid of a print statement

! Revision 2.54  2003/01/26 04:43:15  livesey
! Changed to handle units for phiWindow

! Revision 2.53  2003/01/16 23:13:35  livesey
! Now gets MaxRefraction (smaller value) from Geometry

! Revision 2.52  2002/11/22 12:22:41  mjf
! Use getQuantityForForwardModel instead of GetVectorQuantityByType.

! Revision 2.51  2002/10/25 22:24:54  livesey
! Two bug fixes with the h2o stuff

! Revision 2.50  2002/10/16 20:13:55  mjf
! Added GetGPHPrecision, based on GetBasisGPH.

! Revision 2.49  2002/10/08 17:36:22  pwagner
! Added idents to survive zealous Lahey optimizer

! Revision 2.48  2002/09/27 01:47:45  vsnyder
! Move some USEs from module scope to procedure scope.  Move parameters from
! module scope to procedure scope.  Simplify iteration in z_surface.
! Convert some variables from pointers to automatic arrays.  Cosmetic changes.

! Revision 2.47  2002/09/27 00:01:16  vsnyder
! Remove unused variables and dead calculations

! Revision 2.46  2002/09/26 23:48:21  vsnyder
! Avoid computing sine and cosine of the same angle

! Revision 2.45  2002/09/26 20:38:19  vsnyder
! Get some constants from Geometry and Units instead of declaring them, cosmetics

! Revision 2.44  2002/06/27 00:20:43  livesey
! Typo!

! Revision 2.43  2002/06/27 00:19:26  livesey
! Fixed another log h2o

! Revision 2.42  2002/06/26 23:37:27  livesey
! Tidied up dumps.

! Revision 2.41  2002/06/26 23:29:52  bill
! fixed residual calculation bug--wgr

! Revision 2.40  2002/06/26 20:59:25  livesey
! Fixed bug in 2d pressure guesser, was ignoring very first guess

! Revision 2.39  2002/06/26 19:18:53  livesey
! Merge of Bills fixes and my diagnostics

! Revision 2.38  2002/06/26 18:44:12  bill
! added a feature to limit the index of refraction to its surface value--wgr

! Revision 2.37  2002/06/26 01:26:10  livesey
! Added 2D pressure guesser

! Revision 2.36  2002/06/25 22:17:09  bill
! doesn't store blocks of zeros--wgr

! Revision 2.35  2002/06/25 21:55:49  bill
! first debugged? version of 2d scan module--wgr

! Revision 2.34  2002/06/25 17:03:30  bill
! fixed p4 calc in pressure guesser--wgr

! Revision 2.33  2002/06/25 14:55:15  bill
! work in progress--wgr

! Revision 2.32  2002/06/25 00:01:32  bill
! work in progress--wgr

! Revision 2.31  2002/06/24 18:27:02  livesey
! Debugging

! Revision 2.30  2002/06/24 17:55:41  bill
! added a two d scan model subroutine--wgr

! Revision 2.29  2002/06/05 00:00:29  livesey
! Typo fixed

! Revision 2.28  2002/02/13 00:08:58  livesey
! Added differential scan model

! Revision 2.27  2001/12/06 23:44:57  livesey
! Moved Omega into units so L1BOASim can use it.

! Revision 2.26  2001/11/03 01:33:47  livesey
! Changed the togle from gen to emit

! Revision 2.25  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models

! Revision 2.24  2001/06/29 23:03:28  livesey
! Whoops, bad parameter set.

! Revision 2.23  2001/06/19 22:43:23  pwagner
! Eliminated l_temperature from things got from init_tables_module

! Revision 2.22  2001/06/04 22:42:36  livesey
! Gets belowRef from intermediate

! Revision 2.21  2001/05/10 00:46:49  livesey
! Changed else where to elsewhere (NAG problem?)

! Revision 2.20  2001/05/09 17:43:44  vsnyder
! Use UpdateDiagonal to store PTAN derivatives; cosmetic changes

! Revision 2.19  2001/05/08 19:43:57  livesey
! Pretty stable version.  Pressure now converges better.  Caught a few
! more gremlins in d(residual)/dT.

! Revision 2.18  2001/05/08 03:05:05  livesey
! Working version, pressure guesser a little unstable in extreme cases

! Revision 2.17  2001/05/05 00:03:30  livesey
! Close to working.  One last part of the temperature derivatives that disagrees with
! numerical calculation.

! Revision 2.16  2001/05/04 05:41:26  livesey
! Added some more conditions for derivative calculation

! Revision 2.15  2001/05/04 05:05:46  livesey
! The scan calculation works in ScanModel, haven't tested the
! derivatives yet though.

! Revision 2.14  2001/05/03 23:42:57  livesey
! Closer to working

! Revision 2.13  2001/05/03 23:09:52  livesey
! First version of scan model.  Compiles, but probably doesn't run.

! Revision 2.12  2001/05/03 20:34:08  vsnyder
! Cosmetic changes

! Revision 2.11  2001/04/21 00:52:24  livesey
! Fixed memory leak.

! Revision 2.10  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.

! Revision 2.9  2001/04/03 19:03:07  vsnyder
! Get L_H2O from Molecules -- can't get it from Intrinsic after the order of
! Molecules and Intrinsic was reversed.

! Revision 2.8  2001/03/20 21:43:41  livesey
! Moved geodtogeoc lat to geometry in lib

! Revision 2.7  2001/03/15 23:31:15  livesey
! Made it default private.

! Revision 2.6  2001/03/07 22:42:38  livesey
! Got pressure guesser working

! Revision 2.5  2001/03/06 00:34:54  livesey
! Pretty good version

! Revision 2.4  2001/03/05 01:20:26  livesey
! Interim version

! Revision 2.3  2001/02/28 17:35:36  livesey
! Added RCS log stuff

