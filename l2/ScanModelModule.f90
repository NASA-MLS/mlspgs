! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

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

  use Allocate_Deallocate, only: ALLOCATE_TEST, DEALLOCATE_TEST
  use Dump_0, only: DUMP
  use ForwardModelConfig, only: ForwardModelConfig_T
  use ForwardModelIntermediate, only: ForwardModelIntermediate_T, &
    & ForwardModelStatus_T
  use Geometry, only: EARTHRADA, EARTHRADB, GEODTOGEOCLAT, LN10, PI
  use Init_Tables_Module, only: L_HEIGHT, L_REFGPH, L_ZETA
  use intrinsic, only: L_HEIGHTOFFSET, L_NONE, L_PTAN, L_SCANRESIDUAL, &
    & L_TEMPERATURE, L_TNGTGEOCALT, L_VMR
  use ManipulateVectorQuantities, only: FINDCLOSESTINSTANCES
  use MatrixModule_0, only: DESTROYBLOCK, MATRIXELEMENT_T, M_ABSENT, M_BANDED, &
    & M_FULL, UpdateDiagonal
  use MatrixModule_1, only: CREATEBLOCK, FINDBLOCK, MATRIX_T
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMESSAGE, MLSMSG_ALLOCATE, MLSMSG_DEALLOCATE, &
       MLSMSG_ERROR, MLSMSG_WARNING
  use MLSNumerics, only : HUNT, INTERPOLATEVALUES
  use Molecules, only: L_H2O
  use Output_M, only: OUTPUT
  use Toggles, only: GEN, TOGGLE
  use Trace_M, only: TRACE_BEGIN, TRACE_END
  use VectorsModule, only : GETVECTORQUANTITYBYTYPE, VALIDATEVECTORQUANTITY, &
    & VECTOR_T, VECTORVALUE_T

  implicit none

  private

  public :: GetBasisGPH, GetHydrostaticTangentPressure, Omega, ScanForwardModel

  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=130), private :: Id = &
    & "$Id$"
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Define various constants etc.

  ! First some constants to do with the earths dimensions and rotation

  real (r8), parameter :: OMEGA = 7.292115D-5 ! Earth's angular velocity (s-1)
  real (r8), parameter :: EARTHSURFACEGPH = 6387182.265D0 ! GPH at earth's surface m

  ! Now some other constants to do with geopotentials and GPHs
  
  real (r8), parameter :: J2=0.0010826256D0 ! 2nd coefficient
  real (r8), parameter :: J4=-0.0000016165D0 ! 4th coeficient.
  real (r8), parameter :: GM=3.98600436D14 ! Big G times earth mass (m2s-2)
  real (r8), parameter :: G0=9.80665    ! Nominal little g ms-2
  
  ! Now some terms to do with refraction
  
  real (r8), parameter :: REFRATERM = 0.0000776D0 ! First term
  real (r8), parameter :: REFRBTERM = 4810.0D0 ! Second term
  real (r8), parameter :: MAXREFRACTION = 10.0!1e-2 ! Don't allow stupidly large n.
  real (r8), parameter :: MAXPRESSURE = 1400.0 ! /mb Don't allow very large pressures
  
  ! Now some terms to do with the gas 'constant'
  
  real (r8), parameter :: GASM0 = 25.34314957d-3 ! Constant
  real (r8), parameter :: GASM1 = 2.89644d-3 ! Linear term
  real (r8), parameter :: GASM2 = -0.579d-3 ! Quadratic term
  real (r8), parameter :: GASR0 = 8.31441 ! `Standard' gas constant

contains ! =============== Subroutines and functions ==========================

  ! ----------------------------------------------- GetBasisGPH ---------------
  subroutine GetBasisGPH(temp,refGPH,gph,R,RT,belowRef)
    ! This function takes a state vector, containing one and only one temperature
    ! and reference geopotential height quantity, and returns

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: TEMP ! The temperature field
    type (VectorValue_T), intent(IN) :: REFGPH ! The reference gph field
    real (r8), dimension(:,:), intent(OUT) :: GPH ! Result (temp%noSurfs,temp%noInstances)
    real (r8), dimension(:), pointer, optional :: R ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer, optional :: RT ! R*T (noSurfs,noInstances)
    integer, intent(out), optional :: BELOWREF ! Level in temperature basis

    ! Local variables, many automatic arrays

    real (r8), dimension(:), pointer :: MYR      ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer :: MYRT   ! R*T (noSurfs,noInstances)

    real (r8), dimension(temp%template%noSurfs) :: LOGP ! -log10 pressure
    real (r8), dimension(temp%template%noSurfs) :: MODIFIEDBASIS ! noSurfs
    real (r8), dimension(temp%template%noInstances) :: CURRENTREFGPH ! From 1st calc
    real (r8), dimension(temp%template%noInstances) :: CORRECTION ! To apply to gph
    real (r8), dimension(temp%template%noInstances) :: DELTAGEOPOT ! noInstances

    integer :: MYBELOWREF                 ! Result of a hunt
    real (r8) :: ABOVEREFWEIGHT         ! Interpolation weight
    
    integer :: INSTANCE                 ! Loop counter
    integer :: SURF                     ! Loop counter
    integer :: STATUS                   ! Flag
    
    real (r8) :: BASISCUTOFF            ! Threshold level for gas constant
    real (r8) :: REFLOGP                ! Log p of pressure reference surface
    real (r8) :: BASISGAP               ! Space between adjacent surfaces

    nullify ( myR, myRT )
    ! Check that we get the right kinds of quantities
    if ( ( .not. ValidateVectorQuantity( temp,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=(/l_Zeta/)) ) .or. &
      &  ( .not. ValidateVectorQuantity( refGPH,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=(/l_Zeta/)) ) ) &
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

    basisCutoff= -gasM1 / (2*gasM2) ! Set a threshold value
    modifiedBasis = max ( logP, basisCutoff ) ! Either logP or this threshold
    myR = gasR0 / ( gasM0 + gasM1*modifiedBasis + gasM2*modifiedBasis**2 )

    ! Compute R*T for each point, avoid spread intrinsic to save memory
    ! and cpu time.
    do instance = 1, temp%template%noInstances
       myRT(:,instance) = myR * temp%values(:,instance)
    end do

    gph(1,:) = 0.0
    do surf = 2, temp%template%noSurfs
       deltaGeopot = (ln10/ (2*g0) ) * &
         & ( myRT(surf,:) + myRT(surf-1,:) ) * &
         & ( logP(surf) - logP(surf-1) )
       gph(surf,:) = gph(surf-1,:) + deltaGeopot
    end do

    ! Now we need to correct for the reference geopotential, find the layer the
    ! reference surface is within.

    refLogP = refGPH%template%surfs(1,1)
    call Hunt( logP, refLogP, myBelowRef )

    ! Get weights
    basisGap = logP(myBelowRef+1) - logP(myBelowRef)
    aboveRefWeight = ( refLogP - logP(myBelowRef) )/basisGap

    ! Forbid extrapolation
    aboveRefWeight = max ( min ( aboveRefWeight, 1.0D0 ), 0.0D0 )

    ! Get geopotential at the reference surface from our intermediate result
    currentRefGPH = gph(myBelowRef,:) + ((basisGap*ln10)/(2*g0))* &
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
       deallocate ( myR )
    end if

    if ( present(rt) ) then
       rt=>myRT
    else
       deallocate ( myRT )
    end if

    if ( present(belowRef) ) belowRef = myBelowRef
         
    ! That's it  
  end subroutine GetBasisGPH

  ! ---------------------------------- GetHydroStaticTangentPressure ----------
  subroutine GetHydrostaticTangentPressure ( ptan, temp, refGPH, h2o, geocAlt,&
    maxIterations )
    ! This routine is a pressure `guesser'.  It works by comparing the
    ! geopotential heights estimated from the level 1 heights with those based on
    ! a hydrostatic calculation.  The tangent point pressure is adjusted over
    ! several iterations until a good match is achieved.

    ! Dummy arguments
    type (VectorValue_T), intent(inout) :: PTAN ! Tangent pressure sv. component
    type (VectorValue_T), intent(in) :: TEMP    ! Temperature
    type (VectorValue_T), intent(in) :: REFGPH  ! Reference GPH
    type (VectorValue_T), intent(in) :: H2O     ! H2O
    type (VectorValue_T), intent(in) :: GEOCALT ! L1B Tp Geoc alt.
    integer, intent(IN) :: MAXITERATIONS        ! Number of iterations to use

    ! Local variables

    real (r8), dimension(:,:), pointer :: RT           ! rt=R*T
    real (r8), dimension(:,:), pointer :: BASISGPH     ! temp(noSurfs,noInstances)

    real (r8), dimension(:), pointer :: ACOEFF         ! Quadratic term
    real (r8), dimension(:), pointer :: BASISLOWER     ! For temperature
    real (r8), dimension(:), pointer :: BASISSPACING   ! For temperature
    real (r8), dimension(:), pointer :: BASISUPPER     ! For temperature
    real (r8), dimension(:), pointer :: BCOEFF         ! Quadratic term
    real (r8), dimension(:), pointer :: CCOEFF         ! Quadratic term
    real (r8), dimension(:), pointer :: DELTART 
    real (r8), dimension(:), pointer :: EARTHRADIUS
    real (r8), dimension(:), pointer :: GEOCLAT
    real (r8), dimension(:), pointer :: GEOMETRICGPH
    real (r8), dimension(:), pointer :: H2OBASIS ! Points to h2o basis
    real (r8), dimension(:), pointer :: H2OVALS ! Points to part of h2o
    real (r8), dimension(:), pointer :: N              ! Refractive index
    real (r8), dimension(:), pointer :: P2
    real (r8), dimension(:), pointer :: P4
    real (r8), dimension(:), pointer :: POINTINGH2O    ! t.p. h2o
    real (r8), dimension(:), pointer :: POINTINGPRES ! t.p. press
    real (r8), dimension(:), pointer :: POINTINGTEMP   ! t.p. temp.
    real (r8), dimension(:), pointer :: PTANVALS       ! Points to part of ptan
    real (r8), dimension(:), pointer :: RATIO2         ! minor frame
    real (r8), dimension(:), pointer :: RATIO4         ! minor frame
    real (r8), dimension(:), pointer :: REFRACTEDGEOCALT ! minor frame
    real (r8), dimension(:), pointer :: RTLOWER
    real (r8), dimension(:), pointer :: RTUPPER
    real (r8), dimension(:), pointer :: S2
    real (r8), dimension(:), pointer :: TEMPBASIS ! Points to temp basis
    real (r8), dimension(:), pointer :: TEMPVALS ! Points to part of temp
    real (r8), dimension(:), pointer :: THISBASISGPH ! Points to part of basisGPH

    integer, dimension(:), pointer :: CLOSESTTEMPPROFILES
    integer, dimension(:), pointer :: CLOSESTH2OPROFILES
    integer, dimension(:), pointer :: LOWER ! index into temperature profile
    integer, dimension(:), pointer :: UPPER ! index into temperature profile

    integer :: iteration                ! Loop stuff
    integer :: maf                      ! Loop counter

    ! Executable code

    if ( toggle(gen) ) call trace_begin ('GetHydrostaticTangentPressure' )

    nullify ( aCoeff, basisGPH, basisLower, basisSpacing, basisUpper, &
      & bCoeff, cCoeff, closestH2OProfiles, closestTempProfiles, deltaRT, &
      & earthRadius, geocLat, geometricGPH, lower, n, &
      & pointingH2O, pointingPres, pointingTemp, ratio2, ratio4, &
      & refractedGeocAlt, rt, rtLower, rtUpper, upper, s2, p2, p4 )

    ! Check that we get the right kinds of quantities
    if ( ( .not. ValidateVectorQuantity( temp,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=(/l_Zeta/)) ) .or. &
      &  ( .not. ValidateVectorQuantity( refGPH,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=(/l_Zeta/)) ) .or. &
      &  ( .not. ValidateVectorQuantity( h2o,&
      &            coherent=.true., &
      &            stacked=.true., &
      &            regular=.true., &
      &            verticalCoordinate=(/l_Zeta/)) ) ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
      &   'Inappropriate temp/refGPH/h2o quantity' )

    ! Allocate temporary arrays
    call Allocate_Test ( basisGPH, temp%template%noSurfs, temp%template%noInstances, &
      & "basisGPH", ModuleName )
    call Allocate_Test ( earthRadius, ptan%template%noSurfs, &
      "earthRadius", ModuleName )
    call Allocate_Test ( geocLat, ptan%template%noSurfs, &
      "geocLat", ModuleName )
    call Allocate_Test ( s2, ptan%template%noSurfs, "s2", ModuleName )
    call Allocate_Test ( p2, ptan%template%noSurfs, "p2", ModuleName )
    call Allocate_Test ( p4, ptan%template%noSurfs, "p4", ModuleName )
    call Allocate_Test ( aCoeff, ptan%template%noSurfs, "aCoeff", ModuleName )
    call Allocate_Test ( rtLower, ptan%template%noSurfs, "rtLower", ModuleName )
    call Allocate_Test ( rtUpper, ptan%template%noSurfs, "rtUpper", ModuleName )
    call Allocate_Test ( basisLower, ptan%template%noSurfs, "basisLower", ModuleName )
    call Allocate_Test ( basisUpper, ptan%template%noSurfs, "basisUpper", ModuleName )
    call Allocate_Test ( basisSpacing, ptan%template%noSurfs, "basisSpacing", ModuleName )
    call Allocate_Test ( bCoeff, ptan%template%noSurfs, "bCoeff", ModuleName )
    call Allocate_Test ( cCoeff, ptan%template%noSurfs, "cCoeff", ModuleName )
    call Allocate_Test ( deltaRT, ptan%template%noSurfs, "deltaRT", ModuleName )
    call Allocate_Test ( geometricGPH, ptan%template%noSurfs, &
      "geometricGPH", ModuleName )
    call Allocate_Test ( n, ptan%template%noSurfs, "n", ModuleName )
    call Allocate_Test ( pointingH2O, ptan%template%noSurfs, "pointingH2O", ModuleName )
    call Allocate_Test ( pointingpres, ptan%template%noSurfs, "pointingpres", ModuleName )
    call Allocate_Test ( pointingTemp, ptan%template%noSurfs, "pointingTemp", ModuleName )
    call Allocate_Test ( ratio2, ptan%template%noSurfs, "ratio2", ModuleName )
    call Allocate_Test ( ratio4, ptan%template%noSurfs, "ratio4", ModuleName )
    call Allocate_Test ( refractedGeocAlt, ptan%template%noSurfs, &
      & "refractedGeocAlt", ModuleName )
    call Allocate_Test ( closestTempProfiles, ptan%template%noInstances, &
      "closestTempProfiles", ModuleName )
    call Allocate_Test ( closestH2OProfiles, ptan%template%noInstances, &
      "closestH2OProfiles", ModuleName )
    call Allocate_Test ( lower, ptan%template%noSurfs, "lower", ModuleName )
    call Allocate_Test ( upper, ptan%template%noSurfs, "upper", ModuleName )

    ! Now precompute the basis geopotential height
    call GetBasisGPH ( temp, refGPH, basisGPH, rt=rt )

    ! Find the closest temperature and h2o profiles to each MAF.  Note that
    ! this makes this calculation 1D
    call FindClosestInstances ( temp, ptan, closestTempProfiles )
    call FindClosestInstances ( h2o, ptan, closestH2OProfiles )

    ! Rather than try to be too clever and do this all with 2D arrays, we'll
    ! loop over major frame.

    tempBasis => temp%template%surfs(:,1)
    h2oBasis => h2o%template%surfs(:,1)

    do maf=1,ptan%template%noInstances
      ! Setup pointers for this maf
      ptanVals => ptan%values(:,maf)
      tempVals => temp%values(:,closestTempProfiles(maf))
      h2oVals => h2o%values(:,closestH2OProfiles(maf))
      thisBasisGPH => basisGPH(:,closestTempProfiles(maf))

      ! Get a really simple first guess, assuming a uniform log scale height of
      ! 16km
      geocLat=GeodToGeocLat(ptan%template%geodLat(:,maf))
      earthRadius= earthRadA*earthRadB/sqrt(&
        & (earthRadA*sin(geocLat))**2+ &
        & (earthRadB*cos(geocLat))**2)

      ptanVals = -3.0+(geocAlt%values(:,maf) - &
        & earthRadius)/16e3 ! Do a better job later !???
      ! Precompute some spherical geometry terms, these are
      ! particularly worthy of note. We're doing the geometricGPH
      ! in a 2D manner (apart from refraction effects).  The hydrostatic GPH is
      ! being done 2D.
      s2=sin(geocLat)**2
      p2=0.5*(3*s2-1)
      p4=0.125*(35*(s2**2)-30*s2+30)

      ! Now we have an iteration loop where we try to fit the tangent pressure
      do iteration=1,maxIterations
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

          pointingpres = min(10**(-ptanVals), maxPressure)

          n=( (refrATerm*pointingpres) / pointingTemp ) * &
            & (1.0 + refrBTerm*pointingH2O/pointingTemp)

          n=min(n,maxRefraction) ! Forbid stupidly large values of n           

          ! Now we're going to compute refracted geocentric altitudes
          refractedGeocAlt=geocAlt%values(:,maf)/(1.0+n)

          ! Now convert these to geopotential heights
          ratio2=(earthRadA/refractedGeocAlt)**2
          ratio4=ratio2**2

          geometricGPH= -(GM/(refractedGeocAlt*g0))*&
            &              (1-j2*p2*ratio2- j4*p4*ratio4) - &
            ((omega*refractedGeocAlt*cos(geocLat))**2)/(2*g0)+earthSurfaceGPH
        elsewhere
          n=0.0
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
            & (g0 / (ln10*rtLower))

          ! Above the top
        elsewhere ( geometricGPH >= thisBasisGPH(temp%template%noSurfs) )
          ptanVals = tempBasis(temp%template%noSurfs) + &
            & (geometricGPH-thisBasisGPH(upper)) * &
            & (g0 / (ln10*rtUpper))

          ! Everywhere else
        elsewhere
          aCoeff = deltaRT
          bCoeff = 2*rtLower
          cCoeff = -2*g0*(geometricGPH-thisBasisGPH(lower)) / &
            & (basisSpacing*ln10)
          ! Let ptan contain upperWeight (ptan-basisLower)/basisSpacing first
          ptanVals = 2*cCoeff/( -bCoeff - &
            & sqrt(max(bCoeff**2-4*aCoeff*cCoeff,0.0_r8)))
          ! The max is here just in case of slips.
          ! Now let it contain ptan
          ptanVals = basisLower + ( ptanVals * basisSpacing )
        end where
      end do                            ! End iteration loop
    end do                              ! Major frame loop

    ! Deallocate all the various things
    call Deallocate_Test ( basisGPH, "basisGPH", ModuleName )
    call Deallocate_Test ( aCoeff, "aCoeff", ModuleName )
    call Deallocate_Test ( rtLower, "rtLower", ModuleName )
    call Deallocate_Test ( rtUpper, "rtUpper", ModuleName )
    call Deallocate_Test ( basisLower, "basisLower", ModuleName )
    call Deallocate_Test ( basisUpper, "basisUpper", ModuleName )
    call Deallocate_Test ( basisSpacing, "basisSpacing", ModuleName )
    call Deallocate_Test ( bCoeff, "bCoeff", ModuleName )
    call Deallocate_Test ( cCoeff, "cCoeff", ModuleName )
    call Deallocate_Test ( deltaRT, "deltaRT", ModuleName )
    call Deallocate_Test ( earthRadius, "earthRadisu", ModuleName )
    call Deallocate_Test ( geocLat, "geocLat", ModuleName )
    call Deallocate_Test ( p2, "p2", ModuleName )
    call Deallocate_Test ( p4, "p4", ModuleName )
    call Deallocate_Test ( s2, "s2", ModuleName )
    call Deallocate_Test ( geometricGPH, "geometricGPH", ModuleName )
    call Deallocate_Test ( n, "n", ModuleName )
    call Deallocate_Test ( pointingH2O, "pointingH2O", ModuleName )
    call Deallocate_Test ( pointingpres, "pointingpres", ModuleName )
    call Deallocate_Test ( pointingTemp, "pointingTemp", ModuleName )
    call Deallocate_Test ( ratio2, "ratio2", ModuleName )
    call Deallocate_Test ( ratio4, "ratio4", ModuleName )
    call Deallocate_Test ( refractedGeocAlt, "refractedGeocAlt", ModuleName )

    call Deallocate_Test ( closestTempProfiles, "closestTempProfiles", ModuleName )
    call Deallocate_Test ( closestH2OProfiles, "closestH2OProfiles", ModuleName )
    call Deallocate_Test ( lower, "lower", ModuleName )
    call Deallocate_Test ( upper, "upper", ModuleName )

    call Deallocate_Test ( rt, "rt", ModuleName ) ! Was allocated in GetBasisGPH

    if ( toggle(gen) ) call trace_end ('GetHydrostaticTangentPressure' )

  end subroutine GetHydrostaticTangentPressure

  ! ---------------------------------------- ScanForwardModel --------------
  subroutine ScanForwardModel ( fmConf, state, extra, &
    & fwmOut, ifm, fmStat, jacobian )
    ! This is the main `scan model' for the module. It compares altitude reported
    ! geometrically, as converted into geopotential height, with a hydrostatic
    ! equivalent.

    ! Note that while this model takes multi-MAF quantities as input, it only
    ! considers one MAF at a time (fmStat%maf).

    ! Dummy arguments
    type (ForwardModelConfig_T), intent(in) :: FMCONF ! Configuration options
    type (Vector_T), intent(in) :: STATE ! The state vector
    type (Vector_T), intent(in) :: EXTRA ! Other stuff in the state vector
    type (Vector_T), intent(inout) :: FWMOUT ! Output vector, residual filled
    type (ForwardModelIntermediate_T), intent(inout) :: IFM ! Workspace type stuff
    type (ForwardModelStatus_T), intent(inout) :: FMSTAT ! Which maf etc.
    type (Matrix_T), intent(inout), optional :: JACOBIAN ! The derivative matrix

    ! Local parameters

    character(len=*), parameter :: INVALIDQUANTITY = "Invalid vector quantity for "

    ! Local variables

    integer :: COL                      ! Block col in jacobian
    integer :: H2OINST                  ! H2O Instance for this MAF
    integer :: I                        ! Loop index
    integer :: J                        ! Loop index
    integer :: LOWER                    ! Index into T profile
    integer :: MAF                      ! Major frame to consider
    integer :: NOMAFS                   ! Dimension
    integer :: NOMIFS                   ! Dimension
    integer :: NOTEMPS                  ! Dimension
    integer :: ROW                      ! Block row in jacobian
    integer :: TEMPINST                 ! Temperature instance for this MAF
    integer :: UPPER                    ! Index into T profile

    integer, dimension(:), pointer :: POINTTEMPLAYER ! Pointing layer
    integer, dimension(:), pointer :: POINTH2OLAYER ! Pointing layer
    integer, dimension(:), pointer :: CLOSESTTEMPPROFILES
    integer, dimension(:), pointer :: CLOSESTH2OPROFILES

    real (r8), target :: ZERO=0.0_r8    ! Need this if no heightOffset
    real (r8), pointer :: HTOFF         ! HeightOffset%values or zero
    real (r8) :: USEPTAN                ! A tangent pressure
    real (r8) :: UPPERWEIGHT            ! A weight
    real (r8) :: BASIS                    ! Index
    real (r8) :: BASISMINUS               ! Index
    real (r8) :: BASISPLUS                ! Index

    ! All these will be dimensioned noMIFs.
    real (r8), dimension(:), pointer :: BASISGPH ! From ifm
    real (r8), dimension(:), pointer :: DH2OBYDPTAN ! Derivative
    real (r8), dimension(:), pointer :: DHYDROSGPHBYDPTAN ! Derivative
    real (r8), dimension(:), pointer :: DL1GPHBYDHEIGHTOFFSET ! Derivaite
    real (r8), dimension(:), pointer :: DL1GPHBYDL1REFRGEOCALT ! Derivaite
    real (r8), dimension(:), pointer :: DL1GPHBYDPTAN ! Derivaite
    real (r8), dimension(:), pointer :: DL1GPHBYDTEMPLOWER ! Derivaite
    real (r8), dimension(:), pointer :: DL1GPHBYDTEMPUPPER ! Derivaite
    real (r8), dimension(:), pointer :: DL1REFRGEOCALTBYDHEIGHTOFFSET ! Derivaite
    real (r8), dimension(:), pointer :: DL1REFRGEOCALTBYDN ! Derivaite
    real (r8), dimension(:), pointer :: DL1REFRGEOCALTBYDPTAN ! Derivaite
    real (r8), dimension(:), pointer :: DL1REFRGEOCALTBYDTEMPLOWER ! Derivaite
    real (r8), dimension(:), pointer :: DL1REFRGEOCALTBYDTEMPUPPER ! Derivaite
    real (r8), dimension(:), pointer :: DNBYDPTAN ! Derivaite
    real (r8), dimension(:), pointer :: DNBYDTEMPLOWER ! Derivaite
    real (r8), dimension(:), pointer :: DNBYDTEMPUPPER ! Derivaite
    real (r8), dimension(:), pointer :: DTEMPBYDPTAN ! Derivative
    real (r8), dimension(:), pointer :: GEOCLAT ! TP GeocLat
    real (r8), dimension(:), pointer :: H2OBASIS ! => h2o%template%surfs
    real (r8), dimension(:), pointer :: H2OBASISGAP ! Basis spacing
    real (r8), dimension(:), pointer :: H2OLOWERWEIGHT ! Weight
    real (r8), dimension(:), pointer :: H2OUPPERWEIGHT ! Weight    
    real (r8), dimension(:), pointer :: H2OVALS ! => h2o%values
    real (r8), dimension(:), pointer :: HYDROSGPH ! GPH from hydrostatic calc.
    real (r8), dimension(:), pointer :: L1GPH ! Geometric geopotential
    real (r8), dimension(:), pointer :: L1REFRGEOCALT ! Refracted geocentric alt
    real (r8), dimension(:), pointer :: N ! Refractive index - 1
    real (r8), dimension(:), pointer :: P2 ! Polynomial term
    real (r8), dimension(:), pointer :: P4 ! Polynomial term
    real (r8), dimension(:), pointer :: POINTINGH2O ! tangent h2o.
    real (r8), dimension(:), pointer :: POINTINGPRES ! tangent pressure.
    real (r8), dimension(:), pointer :: POINTINGTEMP ! tangent temp.
    real (r8), dimension(:), pointer :: POVERTSQUARED ! As name implies
    real (r8), dimension(:), pointer :: PTANVALS ! => ptan%values
    real (r8), dimension(:), pointer :: RATIO2 ! For geopotential calculation
    real (r8), dimension(:), pointer :: RATIO4 ! For geopotential calculation
    real (r8), dimension(:), pointer :: RT ! Gas constant*T
    real (r8), dimension(:), pointer :: S2 ! sin^2 geocLat
    real (r8), dimension(:), pointer :: TEMPBASIS ! => temp%template%surfs
    real (r8), dimension(:), pointer :: TEMPBASISGAP ! Basis spacing
    real (r8), dimension(:), pointer :: TEMPLOWERWEIGHT ! Weight
    real (r8), dimension(:), pointer :: TEMPUPPERWEIGHT ! Weight
    real (r8), dimension(:), pointer :: TEMPVALS ! => temp%values

    real (r8), dimension(:,:), pointer :: A ! Derivative array
    real (r8), dimension(:,:), pointer :: DHYDROSGPHBYDTEMP ! Derivative array

    logical :: TEMPINSTATE              ! Set if temp in state not extra
    logical :: REFGPHINSTATE            ! Set if refGPH in state not extra
    logical :: H2OINSTATE               ! Set if H2O in state, not extra
    logical :: PTANINSTATE              ! Set if ptan in state, not extra
    logical :: HEIGHTOFFSETINSTATE      ! Set if heightOffset in state, not extra

    type (VectorValue_T), pointer :: TEMP ! Temperature component of state
    type (VectorValue_T), pointer :: PTAN ! Ptan component of state
    type (VectorValue_T), pointer :: H2O ! H2O component of state
    type (VectorValue_T), pointer :: REFGPH ! Ref gph component of state
    type (VectorValue_T), pointer :: L1ALT ! Tangent point altitude
    type (VectorValue_T), pointer :: RESIDUAL ! Resulting component of fwmOut
    type (VectorValue_T), pointer :: HEIGHTOFFSET ! Height offset component of state

    type (MatrixElement_T), pointer :: BLOCK ! A matrix block

    ! Executable code -----------------------

    if ( toggle ( gen ) ) call trace_begin ( 'ScanForwardModel' )

    ! Identify the vector quantities from state/extra
    temp => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_temperature, foundInFirst=tempInState )
    ptan => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_ptan, instrumentModule=fmConf%instrumentModule,&
      & foundInFirst=ptanInState )
    h2o => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_vmr, molecule=l_h2o, &
      & foundInFirst=h2oInState, noError=.true.)
    refgph => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_refGPH, foundInFirst=refGPHInState )
    l1Alt => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_tngtGeocAlt, instrumentModule=fmConf%instrumentModule )
    heightOffset => GetVectorQuantityByType ( state, extra, &
      & quantityType=l_heightOffset, instrumentModule=fmConf%instrumentModule, &
      & noError=.true., foundInFirst=heightOffsetInState )

    ! Identify the vector quantities from fwmOut
    residual => GetVectorQuantityByType ( fwmOut, &
      & quantityType=l_scanResidual, instrumentModule=fmConf%instrumentModule )

    ! Now check that they make sense
    if ( .not. ValidateVectorQuantity(temp, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'temperature' )
    if ( .not. ValidateVectorQuantity(ptan, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'ptan' )
    if ( .not. ValidateVectorQuantity(h2o, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'h2o' )
    if ( .not. ValidateVectorQuantity(refGPH, stacked=.true., coherent=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'refGPH' )
    if ( .not. ValidateVectorQuantity(l1Alt, minorFrame=.true., &
      & frequencyCoordinate=(/l_none/)) ) call MLSMessage ( MLSMSG_Error, &
      & ModuleName, InvalidQuantity//'l1Alt' )
    if ( associated(heightOffset) ) then
      if ( .not. ValidateVectorQuantity(heightOffset, minorFrame=.false., &
        & frequencyCoordinate=(/l_none/),&
        & noInstances=(/1/) ) ) call MLSMessage ( MLSMSG_Error, &
        & ModuleName, InvalidQuantity//'heightOffset' )
    end if

    ! Now setup some standard dimensions etc.
    maf = fmStat%maf
    noMIFs = ptan%template%noSurfs
    noMAFs = ptan%template%noInstances
    noTemps = temp%template%noSurfs

    nullify ( pointTempLayer, pointH2OLayer )
    nullify ( closestTempProfiles, closestH2OProfiles )

    nullify ( dh2obydptan, dhydrosgphbydptan, dl1gphbydheightoffset )
    nullify ( dl1gphbydl1refrgeocalt, dl1gphbydptan, dl1gphbydtemplower )
    nullify ( dl1gphbydtempupper,  dl1refrgeocaltbydheightoffset )
    nullify ( dl1refrgeocaltbydn, dl1refrgeocaltbydptan )
    nullify ( dl1refrgeocaltbydtemplower, dl1refrgeocaltbydtempupper )
    nullify ( dnbydptan, dnbydtemplower, dnbydtempupper, dtempbydptan )
    nullify ( geoclat, h2obasisgap, h2olowerweight )
    nullify ( h2oupperweight, hydrosgph, l1gph, l1refrgeocalt )
    nullify ( n, p2, p4, pointingh2o, pointingpres, pointingtemp )
    nullify ( povertsquared, ratio2, ratio4 )
    nullify ( rt, s2, tempbasisgap )
    nullify ( templowerweight, tempupperweight )

    nullify ( A, dHydrosGPHByDTemp )

    ! Now set up pointer arrays (there are a lot of these)
    call allocate_test ( pointTempLayer, noMIFs,&
      & 'pointTempLayer', ModuleName )
    call allocate_test ( pointH2OLayer, noMIFs, &
      & 'pointH2OLayer', ModuleName )
    call allocate_test ( closestTempProfiles, noMAFs,&
      & 'closestTempProfiles', ModuleName )
    call allocate_test ( closestH2OProfiles, noMAFs,&
      & 'closestH2OProfiles', ModuleName )

    call allocate_test ( dH2OByDPtan, noMIFs, &
      & 'dH2OByDPtan', ModuleName )
    call allocate_test ( dHydrosGPHByDPtan, noMIFs, &
      & 'dHydrosGPHByDPtan', ModuleName )
    call allocate_test ( dL1GPHByDHeightOffset, noMIFs, &
      & 'dL1GPHByDHeightOffset', ModuleName )
    call allocate_test ( dL1GPHByDL1RefrGeocAlt, noMIFs, &
      & 'dL1GPHByDL1RefrGeocAlt', ModuleName )
    call allocate_test ( dL1GPHByDPtan, noMIFs, &
      & 'dL1GPHByDPtan', ModuleName )
    call allocate_test ( dL1GPHByDTempLower, noMIFs, &
      & 'dL1GPHByDTempLower', ModuleName )
    call allocate_test ( dL1GPHByDTempUpper, noMIFs, &
      & 'dL1GPHByDTempLower', ModuleName )
    call allocate_test ( dL1RefrGeocAltByDHeightOffset, noMIFs, &
      & 'dL1RefrGeocAltByDHeightOffset', ModuleName )
    call allocate_test ( dL1RefrGeocAltByDN, noMIFs, &
      & 'dL1RefrGeocAltByDN', ModuleName )
    call allocate_test ( dL1RefrGeocAltByDPtan, noMIFs, &
      & 'dL1RefrGeocAltByDPtan', ModuleName )
    call allocate_test ( dL1RefrGeocAltByDTempLower, noMIFs, &
      & 'dL1RefrGeocAltByDTempLower', ModuleName )
    call allocate_test ( dL1RefrGeocAltByDTempUpper, noMIFs, &
      & 'dL1RefrGeocAltByDTempUpper', ModuleName )
    call allocate_test ( dNByDPtan, noMIFs, &
      & 'dNByDPtan', ModuleName )
    call allocate_test ( dNByDTempLower, noMIFs, &
      & 'dNByDTempLower', ModuleName )
    call allocate_test ( dNByDTempUpper, noMIFs, &
      & 'dNByDTempUpper', ModuleName )
    call allocate_test ( dTempByDPtan, noMIFs, &
      & 'dTempByDPtan', ModuleName )
    call allocate_test ( geocLat, noMIFs, &
      & 'geocLat', ModuleName )
    call allocate_test ( h2oBasisGap, noMIFs, &
      & 'h2oBasisGap', ModuleName )
    call allocate_test ( h2oLowerWeight, noMIFs, &
      & 'h2oLowerWeight', ModuleName )
    call allocate_test ( h2oUpperWeight, noMIFs, &
      & 'h2oUpperWeight', ModuleName )
    call allocate_test ( hydrosGPH, noMIFs, &
      & 'hydrosGPH', ModuleName )
    call allocate_test ( l1GPH, noMIFs, &
      & 'l1GPH', ModuleName )
    call allocate_test ( l1RefrGeocAlt, noMIFs, &
      & 'l1RefrGeocAlt', ModuleName )
    call allocate_test ( n, noMIFs, &
      & 'n', ModuleName )
    call allocate_test ( p2, noMIFs, &
      & 'p2', ModuleName )
    call allocate_test ( p4, noMIFs, &
      & 'p4', ModuleName )
    call allocate_test ( pointingH2O, noMIFs, &
      & 'pointingH2O', ModuleName )
    call allocate_test ( pointingPres, noMIFs, &
      & 'pointingPres', ModuleName )
    call allocate_test ( pointingTemp, noMIFs, &
      & 'pointingTemp', ModuleName )
    call allocate_test ( pOverTSquared, noMIFs, &
      & 'pOverTSquared', ModuleName )
    call allocate_test ( ratio2, noMIFs, &
      & 'ratio2', ModuleName )
    call allocate_test ( ratio4, noMIFs, &
      & 'ratio4', ModuleName )
    call allocate_test ( s2, noMIFs, &
      & 's2', ModuleName )
    call allocate_test ( tempBasisGap, noMIFs, &
      & 'tempBasisGap', ModuleName )
    call allocate_test ( tempLowerWeight, noMIFs, &
      & 'tempLowerWeight', ModuleName )
    call allocate_test ( tempUpperWeight, noMIFs, &
      & 'tempUpperWeight', ModuleName )

    call allocate_test ( A, noMIFs+1, noTemps, & ! Extra for refGPH
      & 'A', ModuleName )
    call allocate_test ( dHydrosGPHByDTemp, noMIFs, noTemps, &
      & 'dHydrosGPHByDTemp', ModuleName )

    ! This could maybe be store in ifm, but for the moment, we'll do it each
    ! time
    call FindClosestInstances ( temp, ptan, closestTempProfiles )
    call FindClosestInstances ( h2o, ptan, closestH2OProfiles )

    ! Make some pointers to save time
    tempInst = closestTempProfiles ( maf )
    h2oInst = closestH2OProfiles ( maf )
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
    p4 = 0.125 * (35*(s2**2) - 30*s2 + 30)

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
    endwhere

    where (ptanVals>=h2oBasis(1) .and. ptanVals<h2oBasis(h2o%template%noSurfs))
      dH2oByDPTan = ( h2oVals(pointH2oLayer+1) - h2oVals(pointH2oLayer) ) / &
        & h2oBasisGap
    elsewhere
      dH2oByDPTan = 0.0_r8
    endwhere
  
    ! Now get interpolated temperature and H2O
    pointingTemp = tempVals(pointTempLayer) * tempLowerWeight + &
      &            tempVals(pointTempLayer+1) * tempUpperWeight 
    pointingH2O = h2oVals(pointH2OLayer) * h2oLowerWeight + &
      &           h2oVals(pointH2OLayer+1) * h2oUpperWeight

    ! Get pointing pressure in mb
    pointingPres = min(10.0_r8**(-ptanVals), maxPressure)

    ! Now get the refractive index n
    n=refrATerm/(pointingTemp/pointingPres) * &
      & (1.0+refrBterm*pointingH2O/pointingTemp)

    ! Get its derivatives wrt temp, and ptan (ignore h2o)
    pOverTSquared=pointingPres/(pointingTemp**2)

    dNByDPtan=n*(-LN10-dTempByDPtan/pointingTemp)+ &
      & (refrBTerm*refrATerm*pOverTSquared) * &
      & (dH2OByDPtan-(pointingH2O/pointingTemp)*dTempByDPtan)

    dNByDTempLower=-(tempLowerWeight*refrATerm) * &
      & pOverTSquared*(1+2*refrBTerm*pointingH2O/pointingTemp)
    dNByDTempUpper=-(tempUpperWeight*refrATerm) * &
      & pOverTSquared*(1+2*refrBTerm*pointingH2O/pointingTemp)

    ! Now if n is too big limit it and set derivatives to zero
    where ( n > maxRefraction )
      n=maxRefraction
      dNByDTempLower = 0.0
      dNByDTempUpper = 0.0
      dNByDPtan = 0.0
    end where

    ! Do similar things for high pressure cases
    where ( 10.0_r8**(-ptanVals) > maxPressure )
      dNByDPtan=0.0
    end where

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

    l1GPH = -(GM/(l1RefrGeocAlt*g0)) * &
      & (1-J2*P2*ratio2-J4*P4*ratio4) - &
      & ((omega*l1RefrGeocAlt*cos(geocLat))**2)/(2*g0) + earthSurfaceGPH

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
        & rt(pointTempLayer) * ln10 * &
        & (ptanVals-tempBasis(pointTempLayer))/g0
    end where

    where ( ptanVals > tempBasis(noTemps) )
      hydrosGPH = basisGPH (pointTempLayer+1) + &
        & rt(pointTempLayer+1) * ln10 * &
        & (ptanVals-tempBasis(pointTempLayer+1))/g0
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
        upper = ifm%belowRef + 1
      else
        usePtan = ptanVals(i)
        lower = pointTempLayer(i)
        upper = lower+1
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
          if ( j /= noTemps ) then
            upperWeight = (usePtan - basis) / (basisPlus - basis)
          else
            upperWeight = 1.0
          end if
          
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

    ! Compute residual
    residual%values(:,maf) = l1GPH - hydrosGPH

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
        call updateDiagonal ( block, dL1GPHByDPtan-dHydrosGPHByDPtan )
      end if

      ! Store refGPH derivatives
      if ( refGPHInState ) then
        col = FindBlock ( jacobian%col, refGPH%index, maf )
        block => jacobian%block(row,col)
        if ( block%kind /= M_Absent ) then
          call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Found a prexisting d(residual)/d(refGPH), removing' )
          call DestroyBlock ( block )
        end if
        call CreateBlock ( jacobian, row, col, M_Full )
        block%values = -1.0
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
        block%values(:,1) = dL1GPHByDHeightOffset
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
        block%values = - dHydrosGPHByDTemp
        do i = 1, noMIFs
          block%values(i,pointTempLayer(i)) = &
            & block%values(i,pointTempLayer(i)) + dL1GPHByDTempLower(i)
          block%values(i,pointTempLayer(i)+1) = &
            & block%values(i,pointTempLayer(i)+1) + dL1GPHByDTempUpper(i)
        end do
      end if
      
    end if

    ! Now deallocate pointers
    call deallocate_test ( pointTempLayer, 'pointTempLayer', ModuleName )
    call deallocate_test ( pointH2OLayer, 'pointH2OLayer', ModuleName )
    call deallocate_test ( closestTempProfiles, 'closestTempProfiles', ModuleName )
    call deallocate_test ( closestH2OProfiles, 'closestH2OProfiles', ModuleName )

    call deallocate_test ( dH2OByDPtan, 'dH2OByDPtan', ModuleName )
    call deallocate_test ( dHydrosGPHByDPtan, 'dHydrosGPHByDPtan', ModuleName )
    call deallocate_test ( dL1GPHByDHeightOffset, 'dL1GPHByDHeightOffset', ModuleName )
    call deallocate_test ( dL1GPHByDL1RefrGeocAlt, 'dL1GPHByDL1RefrGeocAlt', ModuleName )
    call deallocate_test ( dL1GPHByDPtan, 'dL1GPHByDPtan', ModuleName )
    call deallocate_test ( dL1GPHByDTempLower, 'dL1GPHByDTempLower', ModuleName )
    call deallocate_test ( dL1GPHByDTempUpper, 'dL1GPHByDTempLower', ModuleName )
    call deallocate_test ( dL1RefrGeocAltByDHeightOffset,&
      & 'dL1RefrGeocAltByDHeightOffset', ModuleName )
    call deallocate_test ( dL1RefrGeocAltByDN,&
      & 'dL1RefrGeocAltByDN', ModuleName )
    call deallocate_test ( dL1RefrGeocAltByDPtan,&
      & 'dL1RefrGeocAltByDPtan', ModuleName )
    call deallocate_test ( dL1RefrGeocAltByDTempLower,&
      & 'dL1RefrGeocAltByDTempLower', ModuleName )
    call deallocate_test ( dL1RefrGeocAltByDTempUpper, &
      & 'dL1RefrGeocAltByDTempUpper', ModuleName )
    call deallocate_test ( dNByDPtan, 'dNByDPtan', ModuleName )
    call deallocate_test ( dNByDTempLower, 'dNByDTempLower', ModuleName )
    call deallocate_test ( dNByDTempUpper, 'dNByDTempUpper', ModuleName )
    call deallocate_test ( dTempByDPtan, 'dTempByDPtan', ModuleName )
    call deallocate_test ( geocLat, 'geocLat', ModuleName )
    call deallocate_test ( h2oBasisGap, 'h2oBasisGap', ModuleName )
    call deallocate_test ( h2oLowerWeight, 'h2oLowerWeight', ModuleName )
    call deallocate_test ( h2oUpperWeight, 'h2oUpperWeight', ModuleName )
    call deallocate_test ( hydrosGPH, 'hydrosGPH', ModuleName )
    call deallocate_test ( l1GPH, 'l1GPH', ModuleName )
    call deallocate_test ( l1RefrGeocAlt, 'l1RefrGeocAlt', ModuleName )
    call deallocate_test ( n, 'n', ModuleName )
    call deallocate_test ( p2, 'p2', ModuleName )
    call deallocate_test ( p4, 'p4', ModuleName )
    call deallocate_test ( pointingH2O, 'pointingH2O', ModuleName )
    call deallocate_test ( pointingPres, 'pointingPres', ModuleName )
    call deallocate_test ( pointingTemp, 'pointingTemp', ModuleName )
    call deallocate_test ( pOverTSquared, 'pOverTSquared', ModuleName )
    call deallocate_test ( ratio2, 'ratio2', ModuleName )
    call deallocate_test ( ratio4, 'ratio4', ModuleName )
    call deallocate_test ( s2, 's2', ModuleName )
    call deallocate_test ( tempBasisGap, 'tempBasisGap', ModuleName )
    call deallocate_test ( tempLowerWeight, 'tempLowerWeight', ModuleName )
    call deallocate_test ( tempUpperWeight, 'tempUpperWeight', ModuleName )

    call deallocate_test ( A, 'A', ModuleName )
    call deallocate_test ( dHydrosGPHByDTemp, 'dHydrosGPHByDTemp', ModuleName )

    if ( toggle ( gen ) ) call trace_end ( 'ScanForwardModel' )

  end subroutine ScanForwardModel

end module ScanModelModule

! $Log$
! Revision 2.25  2001/10/02 16:49:56  livesey
! Removed fmStat%finished and change loop ordering in forward models
!
! Revision 2.24  2001/06/29 23:03:28  livesey
! Whoops, bad parameter set.
!
! Revision 2.23  2001/06/19 22:43:23  pwagner
! Eliminated l_temperature from things got from init_tables_module
!
! Revision 2.22  2001/06/04 22:42:36  livesey
! Gets belowRef from intermediate
!
! Revision 2.21  2001/05/10 00:46:49  livesey
! Changed else where to elsewhere (NAG problem?)
!
! Revision 2.20  2001/05/09 17:43:44  vsnyder
! Use UpdateDiagonal to store PTAN derivatives; cosmetic changes
!
! Revision 2.19  2001/05/08 19:43:57  livesey
! Pretty stable version.  Pressure now converges better.  Caught a few
! more gremlins in d(residual)/dT.
!
! Revision 2.18  2001/05/08 03:05:05  livesey
! Working version, pressure guesser a little unstable in extreme cases
!
! Revision 2.17  2001/05/05 00:03:30  livesey
! Close to working.  One last part of the temperature derivatives that disagrees with
! numerical calculation.
!
! Revision 2.16  2001/05/04 05:41:26  livesey
! Added some more conditions for derivative calculation
!
! Revision 2.15  2001/05/04 05:05:46  livesey
! The scan calculation works in ScanModel, haven't tested the
! derivatives yet though.
!
! Revision 2.14  2001/05/03 23:42:57  livesey
! Closer to working
!
! Revision 2.13  2001/05/03 23:09:52  livesey
! First version of scan model.  Compiles, but probably doesn't run.
!
! Revision 2.12  2001/05/03 20:34:08  vsnyder
! Cosmetic changes
!
! Revision 2.11  2001/04/21 00:52:24  livesey
! Fixed memory leak.
!
! Revision 2.10  2001/04/10 22:27:47  vsnyder
! Nullify explicitly instead of with <initialization> so as not to give
! pointers the SAVE attribute.  <initialization> is NOT executed on each
! entry to a procedure.
!
! Revision 2.9  2001/04/03 19:03:07  vsnyder
! Get L_H2O from Molecules -- can't get it from Intrinsic after the order of
! Molecules and Intrinsic was reversed.
!
! Revision 2.8  2001/03/20 21:43:41  livesey
! Moved geodtogeoc lat to geometry in lib
!
! Revision 2.7  2001/03/15 23:31:15  livesey
! Made it default private.
!
! Revision 2.6  2001/03/07 22:42:38  livesey
! Got pressure guesser working
!
! Revision 2.5  2001/03/06 00:34:54  livesey
! Pretty good version
!
! Revision 2.4  2001/03/05 01:20:26  livesey
! Interim version
!
! Revision 2.3  2001/02/28 17:35:36  livesey
! Added RCS log stuff
!
