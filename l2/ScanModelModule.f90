! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module ScanModelModule          ! Scan model and associated calculations
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

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: DUMP
  use Geometry, only: earthRadA,earthRadB, PI, LN10, GeodToGeocLat
  use Init_Tables_Module, only: L_Zeta, L_Temperature, L_RefGPH
  use intrinsic, only: L_TEMPERATURE, L_VMR, L_PTAN, L_TNGTGEOCALT
  use ManipulateVectorQuantities, only: FindClosestInstances
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Deallocate, &
       MLSMSG_Error
  use MLSNumerics, only : Hunt, InterpolateValues
  use Molecules, only: L_H2O
  use Output_M, only: OUTPUT
  use VectorsModule, only : Vector_T, VectorValue_T, GetVectorQuantityByType, &
    &  ValidateVectorQuantity

  use MLSCommon, only: r8

  implicit none

  private

  public :: GetBasisGPH, GetHydrostaticTangentPressure, Omega

  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=130), private :: Id = &
    & "$Id$"
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

  ! Define various constants etc.

  ! First some constants to do with the earths dimensions and rotation

  real (r8), parameter :: omega=7.292115D-5 ! Earth's angular velocity (s-1)
  real (r8), parameter :: earthSurfaceGPH=6387182.265D0 ! GPH at earth's surface m

  ! Now some other constants to do with geopotentials and GPHs
  
  real (r8), parameter :: j2=0.0010826256D0, j4=-0.0000016165D0 ! 2nd and 4th coefs.
  real (r8), parameter :: GM=3.98600436D14 ! Big G times earth mass (m2s-2)
  real (r8), parameter :: g0=9.80665    ! Nominal little g ms-2
  
  ! Now some terms to do with refraction
  
  real (r8), parameter :: refrATerm=0.0000776D0, refrBterm=4810.0D0
  real (r8), parameter :: maxRefraction=1.0 ! Don't allow stupidly large n.
  
  ! Now some terms to do with the gas 'constant'
  
  real (r8), parameter :: gasM0=25.34314957d-3, gasM1=2.89644d-3, gasM2=-0.579d-3
  real (r8), parameter :: gasR0=8.31441

contains ! --------------- Subroutines and functions --------------------------


  ! This function takes a state vector, containing one and only one temperature
  ! and reference geopotential height quantity, and returns

  subroutine GetBasisGPH(temp,refGPH,gph,R,RT)

    ! Dummy arguments
    type (VectorValue_T), intent(IN) :: temp
    type (VectorValue_T), intent(IN) :: refGPH
    real (r8), dimension(:,:), intent(OUT) :: gph ! Result (temp%noSurfs,temp%noInstances)
    real (r8), dimension(:), pointer, optional :: R ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer, optional :: RT ! R*T (noSurfs,noInstances)

    ! Local variables

    real (r8), dimension(:), allocatable :: logP ! -log10 pressure, noSurfs
    real (r8), dimension(:), allocatable :: modifiedBasis ! noSurfs
    real (r8), dimension(:), pointer :: myR      ! Gas constant, noSurfs
    real (r8), dimension(:,:), pointer :: myRT   ! R*T (noSurfs,noInstances)
    real (r8), dimension(:), allocatable :: currentRefGPH ! (noInstances)
    real (r8), dimension(:), allocatable :: correction ! (noInstances)

    real (r8), dimension(:), allocatable :: deltaGeopot ! noInstances

    integer :: belowRef
    real (r8) :: aboveRefWeight ! a weight for an interpolation
    
    integer :: instance,surf ! Loop counters
    integer :: status ! Flag
    
    real (r8) :: basisCutoff    ! Threshold level for gas constant
    real (r8) :: refLogP        ! Log p of pressure reference surface
    real (r8) :: basisGap       ! Space between adjacent surfaces

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
    allocate(logP(temp%template%noSurfs),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "logP in GetBasisGPH")

    allocate(modifiedBasis(temp%template%noSurfs),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "modifiedBasis in GetBasisGPH")

    allocate(myR(temp%template%noSurfs),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "myR in GetBasisGPH")

    allocate(myRT(temp%template%noSurfs,temp%template%noInstances),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "myRT in GetBasisGPH")

    allocate(deltaGeopot(temp%template%noInstances),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "deltaGeopot in GetBasisGPH")

    allocate(currentRefGPH(temp%template%noInstances),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "currentRefGPH in GetBasisGPH")

    allocate(correction(temp%template%noInstances),STAT=status)
    if (status/=0) call MLSMessage(MLSMSG_Error,ModuleName,MLSMSG_Allocate//&
         "correction in GetBasisGPH")

    ! Now the main calculation. This is two parts.  First we compute a
    ! geopotential height grid with an arbitrary offset of H=0 at the lowest
    ! basis point.

    ! To do this we get a gas constant for all the temperature basis points
    logP= temp%template%surfs(:,1)

    basisCutoff=-gasM1/(2*gasM2) ! Set a threshold value
    modifiedBasis=max(logP,basisCutoff) ! Either logP or this threshold
    myR=gasR0/(gasM0+gasM1*modifiedBasis+gasM2*modifiedBasis**2)

    ! Compute R*T for each point, avoid spread intrinsic to save memory
    ! and cpu time.
    do instance=1,temp%template%noInstances
       myRT(:,instance)=myR*temp%values(:,instance)
    end do

    gph(1,:)=0.0
    do surf=2,temp%template%noSurfs
       deltaGeopot=(ln10/(2*g0)) * &
         & (myRT(surf,:)+myRT(surf-1,:)) * &
         & (logP(surf)-logP(surf-1))
       gph(surf,:)=gph(surf-1,:)+deltaGeopot
    end do

    ! Now we need to correct for the reference geopotential, find the layer the
    ! reference surface is within.

    refLogP=refGPH%template%surfs(1,1)
    call Hunt(logP,refLogP,belowRef)

    ! Get weights
    basisGap=logP(belowRef+1)-logP(belowRef)
    aboveRefWeight=(refLogP-logP(belowRef))/basisGap

    ! Forbid extrapolation
    aboveRefWeight=max(min(aboveRefWeight,1.0D0),0.0D0)

    ! Get geopotential at the reference surface from our intermediate result
    currentRefGPH=gph(belowRef,:)+((basisGap*ln10)/(2*g0))* &
      & (myRT(belowRef,:)*aboveRefWeight*(2-aboveRefWeight)+ &
      &  myRT(belowRef+1,:)*(aboveRefWeight**2))

    ! Now make the correction, again avoid spread to save time/memory,
    ! Also convert to geopotential height.
    correction=refGPH%values(1,:)-currentRefGPH
    do surf=1,temp%template%noSurfs         
       gph(surf,:)=gph(surf,:)+correction
    end do
    ! Put back intermediate data if wanted.
    if (present(R)) then
       R=>myR
    else
       deallocate(myR)
    endif

    if (present(rt)) then
       rt=>myRT
    else
       deallocate(myRT)
    endif
         
    ! That's it  
  end subroutine GetBasisGPH

  ! ----------------------------------------------------------------------------

  ! This routine is a pressure `guesser'.  It works by comparing the
  ! geopotential heights estimated from the level 1 heights with those based on
  ! a hydrostatic calculation.  The tangent point pressure is adjusted over
  ! several iterations until a good match is achieved.

  subroutine GetHydrostaticTangentPressure ( ptan, temp, refGPH, h2o, geocAlt,&
       maxIterations )

    ! Dummy arguments
    type (VectorValue_T), intent(INOUT) :: PTAN
    type (VectorValue_T), intent(IN) :: TEMP
    type (VectorValue_T), intent(IN) :: REFGPH
    type (VectorValue_T), intent(IN) :: H2O
    type (VectorValue_T), intent(IN) :: GEOCALT
    integer, intent(IN) :: maxIterations ! Number of iterations to use

    ! Local variables

    real (r8), dimension(:,:), pointer :: rt           ! rt=R*T
    real (r8), dimension(:,:), pointer :: basisGPH     ! temp(noSurfs,noInstances)
    real (r8), dimension(:,:), pointer :: earthRadius
    real (r8), dimension(:,:), pointer :: geocLat

    real (r8), dimension(:), pointer :: aCoeff         ! Quadratic term
    real (r8), dimension(:), pointer :: rtLower
    real (r8), dimension(:), pointer :: rtUpper
    real (r8), dimension(:), pointer :: basisLower     ! For temperature
    real (r8), dimension(:), pointer :: basisUpper     ! For temperature
    real (r8), dimension(:), pointer :: basisSpacing   ! For temperature
    real (r8), dimension(:), pointer :: bCoeff         ! Quadratic term
    real (r8), dimension(:), pointer :: cCoeff         ! Quadratic term
    real (r8), dimension(:), pointer :: deltaRT 
    real (r8), dimension(:), pointer :: geometricGPH
    real (r8), dimension(:), pointer :: n              ! Refractive index
    real (r8), dimension(:), pointer :: pointingH2O    ! t.p. h2o
    real (r8), dimension(:), pointer :: pointingTemp   ! t.p. temp.
    real (r8), dimension(:), pointer :: ratio2         ! minor frame
    real (r8), dimension(:), pointer :: ratio4         ! minor frame
    real (r8), dimension(:), pointer :: refractedGeocAlt ! minor frame

    integer, dimension(:), pointer :: closestTempProfiles
    integer, dimension(:), pointer :: closestH2OProfiles
    integer, dimension(:), pointer :: lower ! index into temperature profile
    integer, dimension(:), pointer :: upper ! index into temperature profile

    real (r8) :: s2, s4, p2, p4 ! Polynomial terms

    real (r8) :: geocLat1               ! Geocentric latitude
    integer :: iteration                ! Loop stuff
    integer :: maf                      ! Loop counter

    ! Executable code

    nullify ( rt, basisGPH, earthRadius, geocLat, aCoeff, rtLower, rtUpper, &
      & basisLower, basisUpper, basisSpacing, bCoeff, cCoeff, deltaRT, &
      & geometricGPH, n, pointingH2O, pointingTemp, ratio2, ratio4, &
      & refractedGeocAlt, closestTempProfiles, closestH2OProfiles, lower, &
      & upper )
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
    call Allocate_Test(basisGPH,temp%template%noSurfs,temp%template%noInstances,&
      & "basisGPH", ModuleName)
    call Allocate_Test(earthRadius,ptan%template%noSurfs,ptan%template%noInstances, &
         "earthRadius",ModuleName)
    call Allocate_Test(geocLat,ptan%template%noSurfs,ptan%template%noInstances, &
         "geocLat",ModuleName)
    call Allocate_Test(aCoeff,ptan%template%noSurfs,"aCoeff",ModuleName)
    call Allocate_Test(rtLower,ptan%template%noSurfs,"rtLower",ModuleName)
    call Allocate_Test(rtUpper,ptan%template%noSurfs,"rtUpper",ModuleName)
    call Allocate_Test(basisLower,ptan%template%noSurfs,"basisLower",ModuleName)
    call Allocate_Test(basisUpper,ptan%template%noSurfs,"basisUpper",ModuleName)
    call Allocate_Test(basisSpacing,ptan%template%noSurfs,"basisSpacing",ModuleName)
    call Allocate_Test(bCoeff,ptan%template%noSurfs,"bCoeff",ModuleName)
    call Allocate_Test(cCoeff,ptan%template%noSurfs,"cCoeff",ModuleName)
    call Allocate_Test(deltaRT,ptan%template%noSurfs,"deltaRT",ModuleName)
    call Allocate_Test(geometricGPH,ptan%template%noSurfs, &
         "geometricGPH",ModuleName)
    call Allocate_Test(n,ptan%template%noSurfs,"n",ModuleName)
    call Allocate_Test(pointingH2O,ptan%template%noSurfs,"pointingH2O",ModuleName)
    call Allocate_Test(pointingTemp,ptan%template%noSurfs,"pointingTemp",ModuleName)
    call Allocate_Test(ratio2,ptan%template%noSurfs,"ratio2",ModuleName)
    call Allocate_Test(ratio4,ptan%template%noSurfs,"ratio4",ModuleName)
    call Allocate_Test(refractedGeocAlt,ptan%template%noSurfs,"refractedGeocAlt",ModuleName)
    
    call Allocate_Test(closestTempProfiles,ptan%template%noInstances,&
         "closestTempProfiles", ModuleName)
    call Allocate_Test(closestH2OProfiles,ptan%template%noInstances,&
         "closestH2OProfiles", ModuleName)
    call Allocate_Test(lower,ptan%template%noSurfs,"lower",ModuleName)
    call Allocate_Test(upper,ptan%template%noSurfs,"upper",ModuleName)

    ! Now precompute the basis geopotential height
    call GetBasisGPH(temp,refGPH,basisGPH,rt=rt)

    ! Get a really simple first guess, assuming a uniform log scale height of
    ! 16km
    geocLat=GeodToGeocLat(ptan%template%geodLat)
    earthRadius= earthRadA*earthRadB/sqrt(&
      & (earthRadA*sin(geocLat))**2+ &
      & (earthRadB*cos(geocLat))**2)

    ptan%values= -3.0+(geocAlt%values-earthRadius)/16e3 ! Do a better job later !???

    ! Find the closest temperature and h2o profiles to each MAF.  Note that
    ! this makes this calculation 1D
    call FindClosestInstances(temp,ptan,closestTempProfiles)
    call FindClosestInstances(h2o,ptan,closestH2OProfiles)

    ! Rather than try to be too clever and do this all with 2D arrays, we'll
    ! loop over major frame.

    do maf=1,ptan%template%noInstances
       ! Precompute some spherical geometry terms, these are
       ! particularly worthy of note.  There is a 1D approximation
       ! the same one we made in early versions of UMLSV5 and had to fix.
       ! Make sure we update this one when we go 2D!
       geocLat1=GeodToGeocLat(temp%template%geodLat(1, &
            closestTempProfiles(maf)))
       s2=sin(geocLat1)**2
       p2=0.5*(3*s2-1)
       p4=0.125*(35*(s2**2)-30*s2+30)

       ! Now we have an iteration loop where we try to fit the tangent pressure
       do iteration=1,maxIterations
          ! The first stage is to refract the tangent point altitudes, for this we
          ! need the refractive index, which is dependent on temperature, pressure
          ! and water vapor concentration for the given altitude.
          
          ! Note here an assumption that temp and h2o are coherent.
          call InterpolateValues( temp%template%surfs(:,1), &
               temp%values(:,closestTempProfiles(maf)), ptan%values(:,maf), &
               pointingTemp, "Linear", extrapolate='C')
          call InterpolateValues( h2o%template%surfs(:,1), &
               h2o%values(:,closestH2OProfiles(maf)), ptan%values(:,maf), &
               pointingH2O, "Linear", extrapolate='C')

          pointingH2O = max(pointingH2O, 0.0_r8)

          where(geocAlt%values(:,maf) /= geocAlt%template%badValue)
            
            n=( (refrATerm*(10**(-ptan%values(:,maf)))) / pointingTemp ) * &
              & (1.0 + refrBTerm*pointingH2O/pointingTemp)
            
            n=min(n,maxRefraction) ! Forbid stupidly large values of n             
            ! Now we're going to compute refracted geocentric altitudes
            refractedGeocAlt=geocAlt%values(:,maf)/(1.0+n)
            
            ! Now convert these to geopotential heights
            ratio2=(earthRadA/refractedGeocAlt)**2
            ratio4=ratio2**2
            
            geometricGPH= -(GM/(refractedGeocAlt*g0))*&
              &              (1-j2*p2*ratio2- j4*p4*ratio4)- &
              ((omega*refractedGeocAlt*cos(geocLat1))**2)/(2*g0)+earthSurfaceGPH
          elsewhere
            n=0.0
          end where

          ! Now, we're effectively going to compare this with a hydrostatic
          ! calculation.

          call Hunt(temp%template%surfs(:,1), ptan%values(:,maf),lower)
          upper=lower+1
          basisLower= temp%template%surfs(lower,1)
          basisUpper= temp%template%surfs(upper,1)
          basisSpacing= basisUpper - basisLower
          rtLower = rt(lower, closestTempProfiles(maf))
          rtUpper = rt(upper, closestTempProfiles(maf))
          deltaRT = rtUpper - rtLower

          where (geocAlt%values(:,maf) /= geocAlt%template%badValue)
            aCoeff = deltaRT
            bCoeff = 2*rtLower
            cCoeff = -2*g0*(geometricGPH-basisGPH(lower,closestTempProfiles(maf))) /  &
              & (basisSpacing*ln10)
            ! Let ptan contain upperWeight (ptan-basisLower)/basisSpacing first
            ptan%values(:,maf)=2*cCoeff/( -bCoeff - &
              & sqrt(max(bCoeff**2-4*aCoeff*cCoeff,0.0_r8)))
            ! The max is here just in case of slips.
            ! Now let it contain ptan
            ptan%values(:,maf) = basisLower+(ptan%values(:,maf) * basisSpacing)
          elsewhere
            ptan%values(:,maf)=ptan%template%badValue
          end where
       end do                   ! End iteration loop
    end do                      ! Major frame loop

    ! Deallocate all the various things
    call Deallocate_Test(basisGPH,"basisGPH", ModuleName)
    call Deallocate_Test(aCoeff,"aCoeff",ModuleName)
    call Deallocate_Test(rtLower,"rtLower",ModuleName)
    call Deallocate_Test(rtUpper,"rtUpper",ModuleName)
    call Deallocate_Test(basisLower,"basisLower",ModuleName)
    call Deallocate_Test(basisUpper,"basisUpper",ModuleName)
    call Deallocate_Test(basisSpacing,"basisSpacing",ModuleName)
    call Deallocate_Test(bCoeff,"bCoeff",ModuleName)
    call Deallocate_Test(cCoeff,"cCoeff",ModuleName)
    call Deallocate_Test(deltaRT,"deltaRT",ModuleName)
    call Deallocate_Test(earthRadius,"earthRadisu",ModuleName)
    call Deallocate_Test(geocLat,"geocLat",ModuleName)
    call Deallocate_Test(geometricGPH,"geometricGPH",ModuleName)
    call Deallocate_Test(n,"n",ModuleName)
    call Deallocate_Test(pointingH2O,"pointingH2O",ModuleName)
    call Deallocate_Test(pointingTemp,"pointingTemp",ModuleName)
    call Deallocate_Test(ratio2,"ratio2",ModuleName)
    call Deallocate_Test(ratio4,"ratio4",ModuleName)
    call Deallocate_Test(refractedGeocAlt,"refractedGeocAlt",ModuleName)
    
    call Deallocate_Test(closestTempProfiles,"closestTempProfiles", ModuleName)
    call Deallocate_Test(closestH2OProfiles,"closestH2OProfiles", ModuleName)
    call Deallocate_Test(lower,"lower",ModuleName)
    call Deallocate_Test(upper,"upper",ModuleName)

    call Deallocate_Test(rt,"rt",ModuleName) ! Was allocated in GetBasisGPH

  end subroutine GetHydrostaticTangentPressure

end module ScanModelModule

! $Log$
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
