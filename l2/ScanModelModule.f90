! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ScanModelModule          ! Scan model and associated calculations
!=============================================================================

  ! This module defines the `scan model' for the MLS Level 2 software.  It also
  ! contains functionality such as a pressure guesser etc, which may be useful
  ! for other programs.

  ! Note that currently the scan model is 1D, later versions will be two
  ! dimensional.  I will try to place comments in the code to indicate where 2D
  ! aware stuff needs to be placed.

  ! Also note that there is currently no height offset term in the state
  ! vector.  This was never used in UMLS V5, and so it is not clear that it
  ! will ever be required.

  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, &
       MLSMSG_Error
  USE VectorsModule, ONLY : Vector_T, VectorValue_T, GetVectorQuantityByType, &
    &  ValidateVectorQuantity
  USE Init_Tables_Module, ONLY: L_Zeta, L_Temperature, L_RefGPH
  USE ManipulateVectorQuantities, ONLY: FindClosestInstances
  USE MLSNumerics, ONLY : Hunt, InterpolateValues
  USE Intrinsic, ONLY: L_H2O, L_TEMPERATURE, L_VMR, L_PTAN, L_TNGTGEODALT, L_TNGTGEOCALT
  USE Allocate_Deallocate, ONLY: Allocate_Test, Deallocate_Test

  USE MLSCommon, ONLY: r8

  IMPLICIT NONE

  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=130), PRIVATE :: Id = &
    & "$Id$"
  CHARACTER (LEN=*), PARAMETER, PRIVATE :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Define various constants etc.

  ! First some constants to do with the earths dimensions and rotation

  REAL (r8), PARAMETER :: earthRadA=6378.137D0 ! Major axis radius in km
  REAL (r8), PARAMETER :: earthRadB=6356.7523141D0 ! Minor axis radius in km
  REAL (r8), PARAMETER :: omega=7.292115D-5 ! Earth's angular velocity (s-1)

  ! Now some other constants to do with geopotentials and GPHs
  
  REAL (r8), PARAMETER :: j2=0.0010826256D0, j4=-0.0000016165D0 ! 2nd and 4th coefs.
  REAL (r8), PARAMETER :: GM=3.98600436D5 ! Big G times earth mass (km2s-2)
  REAL (r8), PARAMETER :: g0=9.80665D-3 ! Nominal little g
  
  ! Now some terms to do with refraction
  
  REAL (r8), PARAMETER :: refrATerm=0.0000776D0, refrBterm=4810.0D0
  REAL (r8), PARAMETER :: maxRefraction=1.0 ! Don't allow stupidly large n.
  
  ! Now some terms to do with the gas 'constant'
  
  REAL (r8), PARAMETER :: gasM0=25.34314957d-3, gasM1=2.896d-3, gasM2=-0.579d-3
  REAL (r8), PARAMETER :: gasR0=8.31432d-6
    
  ! Just in case this constant isn't defined everywhere else.
  
  REAL (r8), PARAMETER :: PI=3.1415926535897931159979635D0
  REAL (r8), PARAMETER :: LN10=2.302585124969482421875D0
  
CONTAINS ! --------------- Subroutines and functions --------------------------

  ! This function converts a geodetic latitude (IN DEGREES!) into a geocentric
  ! one (IN RADIANS!)

  REAL(r8) ELEMENTAL FUNCTION GeodToGeocLat(geodLat)

    ! Arguments
    REAL (r8), INTENT(IN) :: geodLat

    ! Executable code, use special method for high latitudes.
    IF (ABS(geodLat).GT.89.0) THEN
       IF (geodLat.GT.0.0) THEN
          geodtoGeocLat=PI/2.0+((earthRadA/earthRadB)**2)*&
               (geodLat-90.0)*PI/180.0
       ELSE
          geodToGeocLat=-PI/2.0+((earthRadA/earthRadB)**2)*&
               (geodLat+90.0)*PI/180.0
       ENDIF
    ELSE
       geodToGeocLat=ATAN(((earthRadB/earthRadA)**2)*TAN(geodLat*PI/180.0))
    ENDIF
  END FUNCTION GeodToGeocLat
  
  ! --------------------------------------------------------------------------

  ! This function takes a state vector, containing one and only one temperature
  ! and reference geopotential height quantity, and returns

  SUBROUTINE GetBasisGPH(state,gph,R,art)

    ! Dummy arguments
    TYPE (Vector_T), INTENT(IN) :: state
    REAL (r8), DIMENSION(:,:), POINTER :: gph ! Result (temp%noSurfs,temp%noInstances)
    REAL (r8), DIMENSION(:), POINTER, OPTIONAL :: R ! Gas constant, noSurfs
    REAL (r8), DIMENSION(:,:), POINTER, OPTIONAL :: art ! ln10*R*T (noSurfs,noInstances)

    ! Local variables
    TYPE (VectorValue_T), POINTER :: temp
    TYPE (VectorValue_T), POINTER :: refGPH

    REAL (r8), DIMENSION(:), ALLOCATABLE :: logP ! -log10 pressure, noSurfs
    REAL (r8), DIMENSION(:), ALLOCATABLE :: modifiedBasis ! noSurfs
    REAL (r8), DIMENSION(:), POINTER :: myR ! Gas constant, noSurfs
    REAL (r8), DIMENSION(:,:), POINTER :: myArt ! ln10*R*T (noSurfs,noInstances)
    REAL (r8), DIMENSION(:), ALLOCATABLE :: currentRefGeopot ! (noInstances)
    REAL (r8), DIMENSION(:), ALLOCATABLE :: correction ! (noInstances)

    REAL (r8), DIMENSION(:), ALLOCATABLE :: deltaGeopot ! noInstances

    INTEGER :: belowRef
    REAL (r8) :: aboveRefWeight ! a weight for an interpolation
    
    INTEGER :: instance,surf ! Loop counters
    INTEGER :: status ! Flag
    
    REAL (r8) :: basisCutoff    ! Threshold level for gas constant
    REAL (r8) :: refLogP        ! Log p of pressure reference surface
    REAL (r8) :: basisGap       ! Space between adjacent surfaces

    ! First we identify the temperature and refGPH quantities in the vector
    ! Note that no return if no such quantity found.
    temp=>GetVectorQuantityByType(state, l_temperature)
    refGPH=>GetVectorQuantityByType(state, l_refGPH)

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

    ! Now allocate a result array
    IF (ASSOCIATED(gph)) CALL Deallocate_test(gph,"GPH in GetBasisGPH",ModuleName)
    CALL Allocate_Test(gph,temp%template%noSurfs,temp%template%noInstances,&
         "GPH in GetBasisGPH",ModuleName)

    ! Now allocate intermediate quantities
    ALLOCATE(logP(temp%template%noSurfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "logP in GetBasisGPH")

    ALLOCATE(modifiedBasis(temp%template%noSurfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "modifiedBasis in GetBasisGPH")

    ALLOCATE(myR(temp%template%noSurfs),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "myR in GetBasisGPH")

    ALLOCATE(myArt(temp%template%noSurfs,temp%template%noInstances),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "myArt in GetBasisGPH")

    ALLOCATE(deltaGeopot(temp%template%noInstances),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "deltaGeopot in GetBasisGPH")

    ALLOCATE(currentRefGeopot(temp%template%noInstances),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "currentRefGeopot in GetBasisGPH")

    ALLOCATE(correction(temp%template%noInstances),STAT=status)
    IF (status/=0) CALL MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "correction in GetBasisGPH")

    ! Now the main calculation. This is two parts.  First we compute a
    ! geopotential height grid with an arbitrary offset of H=0 at the lowest
    ! basis point.

    ! To do this we get a gas constant for all the temperature basis points
    logP= temp%template%surfs(:,1)

    basisCutoff=-gasM1/(2*gasM2) ! Set a threshold value
    modifiedBasis=MAX(logP,basisCutoff) ! Either logP or this threshold
    myR=gasR0/(gasM0+gasM1*modifiedBasis+gasM2*modifiedBasis**2)
    
    ! Compute ln10*R*T for each point, avoid spread intrinsic to save memory
    ! and cpu time.
    DO instance=1,temp%template%noInstances
       myArt(:,instance)=LN10*myR*temp%values(:,instance)
    END DO

    gph(1,:)=0.0
    DO surf=2,temp%template%noSurfs
       deltaGeopot=(myArt(surf,:)+myArt(surf-1,:))*(logP(surf)-logP(surf-1))/2.0
       gph(surf,:)=gph(surf-1,:)+deltaGeopot
    END DO

    ! Now we need to correct for the reference geopotential, find the layer the
    ! reference surface is within.

    refLogP=refGPH%template%surfs(1,1)
    CALL Hunt(logP,refLogP,belowRef)

    ! Get weights
    basisGap=logP(belowRef+1)-logP(belowRef)
    aboveRefWeight=(refLogP-logP(belowRef))/basisGap

    ! Forbid extrapolation
    aboveRefWeight=MAX(MIN(aboveRefWeight,1.0D0),0.0D0)

    ! Get geopotential at the reference surface from our intermediate result
    currentRefGeopot=gph(belowRef,:)+( &
         myArt(belowRef,:)*aboveRefWeight*(2-aboveRefWeight)+&
         myArt(belowRef+1,:)*(aboveRefWeight**2))*basisGap/2.0

    ! Now make the correction, again avoid spread to save time/memory,
    ! Also convert to geopotential height.
    correction=refGPH%values(1,:)-currentRefGeopot/g0
    DO surf=1,temp%template%noSurfs         
       gph(surf,:)=gph(surf,:)/g0+correction
    END DO

    ! Put back intermediate data if wanted.
    IF (PRESENT(R)) THEN
       R=>myR
    ELSE
       DEALLOCATE(myR)
    ENDIF

    IF (PRESENT(art)) THEN
       art=>myArt
    ELSE
       DEALLOCATE(myArt)
    ENDIF
         
    ! That's it  
  END SUBROUTINE GetBasisGPH

  ! ----------------------------------------------------------------------------

  ! This routine is a pressure `guesser'.  It works by comparing the
  ! geopotential heights estimated from the level 1 heights with those based on
  ! a hydrostatic calculation.  The tangent point pressure is adjusted over
  ! several iterations until a good match is achieved.

  SUBROUTINE GetHydrostaticTangentPressure(state,measurements,radiometer,&
       noIterations)

    ! Dummy arguments
    TYPE (Vector_T), INTENT(INOUT) :: state ! Contains temp, refGPH and ptan
    TYPE (Vector_T), INTENT(IN) :: measurements ! Contains l1 altitudes
    INTEGER, INTENT(IN) :: radiometer ! L_THz or L_GHZ
    INTEGER, INTENT(IN), OPTIONAL :: noIterations ! Number of iterations to use

    ! Local variables
    TYPE (VectorValue_T), POINTER :: geodAlt
    TYPE (VectorValue_T), POINTER :: geocAlt
    TYPE (VectorValue_T), POINTER :: temp
    TYPE (VectorValue_T), POINTER :: refGPH
    TYPE (VectorValue_T), POINTER :: h2o
    TYPE (VectorValue_T), POINTER :: ptan

    REAL (r8), DIMENSION(:,:), POINTER :: basisGPH ! temp(noSurfs,noInstances)
    REAL (r8), DIMENSION(:,:), POINTER :: art ! art=ln10*R*T

    REAL (r8), DIMENSION(:), POINTER :: artLower
    REAL (r8), DIMENSION(:), POINTER :: basisSpacing ! For temperature
    REAL (r8), DIMENSION(:), POINTER :: cCoeff ! Quadratic term
    REAL (r8), DIMENSION(:), POINTER :: deltaArt 
    REAL (r8), DIMENSION(:), POINTER :: geometricGeopotential
    REAL (r8), DIMENSION(:), POINTER :: n ! Refractive index
    REAL (r8), DIMENSION(:), POINTER :: pointingH2O ! t.p. h2o
    REAL (r8), DIMENSION(:), POINTER :: pointingTemp ! t.p. temp.
    REAL (r8), DIMENSION(:), POINTER :: ratio2 ! minor frame
    REAL (r8), DIMENSION(:), POINTER :: ratio4 ! minor frame
    REAL (r8), DIMENSION(:), POINTER :: refractedGeocAlt ! minor frame

    INTEGER, DIMENSION(:), POINTER :: closestTempProfiles
    INTEGER, DIMENSION(:), POINTER :: closestH2OProfiles
    INTEGER, DIMENSION(:), POINTER :: lower ! index into temperature profile

    REAL (r8) :: geocLat        ! Geocentric latitude
    REAL (r8) :: s2, s4, p2, p4 ! Polynomial terms

    INTEGER :: myNoIterations,iteration ! Loop stuff
    INTEGER :: maf              ! Loop counter

    ! Executable code

    ! Set default values

    myNoIterations=4
    IF (PRESENT(noIterations)) myNoIterations=noIterations

    ! The first stage is to identify the vector quantities involved.
    geodAlt=> GetVectorQuantityByType(state, l_tngtGeodAlt, &
         radiometer=radiometer)
    geocAlt=> GetVectorQuantityByType(state, l_tngtGeocAlt, &
         radiometer=radiometer)
    temp=> GetVectorQuantityByType(state, l_temperature)
    refGPH=> GetVectorQuantityByType(state, l_refGPH)
    h2o=> GetVectorQuantityByType(state, l_vmr, molecule=l_h2o)
    ptan=> GetVectorQuantityByType(state, l_ptan, radiometer= &
         radiometer)

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

    if ( ( .not. ValidateVectorQuantity( geodAlt, minorFrame=.true.) ) .or. &
      &  ( .not. ValidateVectorQuantity( geocAlt, minorFrame=.true.) ) .or. &
      &  ( .not. ValidateVectorQuantity( ptan,    minorFrame=.true.) ) ) &
      & call MLSMessage(MLSMSG_Error,ModuleName, &
      &   'Inapporpriate geodAlt/geocAlt/ptan quantity' )

    ! Allocate temporary arrays
    CALL Allocate_Test(artLower,ptan%template%noSurfs,"artLower",ModuleName)
    CALL Allocate_Test(basisSpacing,ptan%template%noSurfs,"basisSpacing",ModuleName)
    CALL Allocate_Test(cCoeff,ptan%template%noSurfs,"cCoeff",ModuleName)
    CALL Allocate_Test(deltaArt,ptan%template%noSurfs,"deltaArt",ModuleName)
    CALL Allocate_Test(geometricGeopotential,ptan%template%noSurfs, &
         "geometricGeopotential",ModuleName)
    CALL Allocate_Test(n,ptan%template%noSurfs,"n",ModuleName)
    CALL Allocate_Test(pointingH2O,ptan%template%noSurfs,"pointingH2O",ModuleName)
    CALL Allocate_Test(pointingTemp,ptan%template%noSurfs,"pointingTemp",ModuleName)
    CALL Allocate_Test(ratio2,ptan%template%noSurfs,"ratio2",ModuleName)
    CALL Allocate_Test(ratio4,ptan%template%noSurfs,"ratio4",ModuleName)
    CALL Allocate_Test(refractedGeocAlt,ptan%template%noSurfs,"refractedGeocAlt",ModuleName)
    
    CALL Allocate_Test(closestTempProfiles,ptan%template%noInstances,&
         "closestTempProfiles", ModuleName)
    CALL Allocate_Test(closestH2OProfiles,ptan%template%noInstances,&
         "closestH2OProfiles", ModuleName)
    CALL Allocate_Test(lower,ptan%template%noSurfs,"lower",ModuleName)

    ! Now precompute the basis geopotential height
    CALL GetBasisGPH(state,basisGPH,art=art)

    ! Get a really simple first guess, assuming a uniform log scale height of
    ! 16km
    ptan%values=3.0-geodAlt%values/16e3

    ! Find the closest temperature and h2o profiles to each MAF.  Note that
    ! this makes this calculation 1D
    CALL FindClosestInstances(temp,ptan,closestTempProfiles)
    CALL FindClosestInstances(h2o,ptan,closestH2OProfiles)

    ! Rather than try to be too clever and do this all with 2D arrays, we'll
    ! loop over major frame.

    DO maf=1,ptan%template%noInstances

       ! Precompute some spherical geometry terms, these are
       ! particularly worthy of note.  There is a 1D approximation
       ! the same one we made in early versions of UMLSV5 and had to fix.
       ! Make sure we update this one when we go 2D!
       geocLat=GeodToGeocLat(temp%template%geodLat(1, &
            closestTempProfiles(maf)))
       s2=SIN(geocLat)**2
       p2=0.5*(3*s2-1)
       p4=0.125*(35*(s2**2)-30*s2+30)

       ! Now we have an iteration loop where we try to fit the tangent pressure
       DO iteration=1,myNoIterations

          ! The first stage is to refract the tangent point altitudes, for this we
          ! need the refractive index, which is dependent on temperature, pressure
          ! and water vapor concentration for the given altitude.
          
          ! Note here an assumption that temp and h2o are coherent.
          CALL InterpolateValues( &
               temp%template%surfs(:,1), &
               temp%values(:,closestTempProfiles(maf)), &
               ptan%values(:,maf), &
               pointingTemp,&
               "Linear")
          CALL InterpolateValues( &
               h2o%template%surfs(:,1), &
               h2o%values(:,closestH2OProfiles(maf)), &
               ptan%values(:,maf), &
               pointingH2O,&
               "Linear")

          WHERE(geocAlt%values(:,maf) /= geocAlt%template%badValue)

             n=(refrATerm/pointingTemp*(10**ptan%values(:,maf)))*&
                  (1+refrBTerm*pointingH2O/pointingTemp)

             n=MIN(n,maxRefraction) ! Forbid stupidly large values of n             
             ! Now we're going to compute refracted geocentric altitudes
             refractedGeocAlt=geocAlt%values(:,maf)/(1.0+n)

             ! Now convert these to geopotential heights
             ratio2=(earthRadA/refractedGeocAlt)**2
             ratio4=ratio2**2

             geometricGeopotential= -(GM/refractedGeocAlt)*(1-j2*p2*ratio2- &
                  j4*p4*ratio4)- &
                  0.5*((omega*refractedGeocAlt*COS(geocLat))**2)
          ENDWHERE

          ! Now, we're effectively going to compare this with a hydrostatic
          ! calculation.

          CALL Hunt(temp%template%surfs(:,1), &
               ptan%values(:,maf),lower)
          basisSpacing= temp%template%surfs(lower+1,1)- &
               temp%template%surfs(lower,1)
          deltaArt=art(closestTempProfiles(maf),lower+1)-&
               art(closestTempProfiles(maf),lower)
          artLower=art(closestTempProfiles(maf),lower)

          WHERE (geocAlt%values(:,maf) /= geocAlt%template%badValue)
             cCoeff=2*(geometricGeopotential-&
                  basisGPH(closestTempProfiles(maf),lower)) &
                  /basisSpacing
             ptan%values(:,maf)=temp%template%surfs(lower,1)+ &
                  basisSpacing*cCoeff/ &
                  (artLower*SQRT(ABS(artLower**2+deltaArt*cCoeff)))
          ELSEWHERE
             ptan%values(:,maf)=ptan%template%badValue
          END WHERE

       END DO                   ! End iteration loop
    END DO                      ! Major frame loop

    ! Deallocate all the various things
    CALL Deallocate_Test(closestTempProfiles,"closestTempProfiles", ModuleName)
    CALL Deallocate_Test(closestH2OProfiles,"closestH2OProfiles", ModuleName)
    CALL Deallocate_Test(pointingTemp,"pointingTemp",ModuleName)
    CALL Deallocate_Test(pointingH2O,"pointingH2O",ModuleName)
    CALL Deallocate_Test(n,"n",ModuleName)
    CALL Deallocate_Test(lower,"lower",ModuleName)

  END SUBROUTINE GetHydrostaticTangentPressure

END MODULE ScanModelModule

! $Log$
! Revision 2.3  2001/02/28 17:35:36  livesey
! Added RCS log stuff
!
