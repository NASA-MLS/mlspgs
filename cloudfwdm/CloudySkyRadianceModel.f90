! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module CloudySkyRadianceModel

! -------------------------------------------------------------------------  
! CLOUDY-SKY MICROWAVE RADIATIVE TRANSFER MODEL      (VERSION 1.0)  
! -------------------------------------------------------------------------
      use AntennaPatterns_m,       only: AntennaPattern_T
      use ClearSkyModule,          only: CLEAR_SKY
      use CloudySkyModule,         only: CLOUDY_SKY
      use FOV_CONVOLVE_M,          only: FOV_CONVOLVE ! , FOV_CONVOLVE_OLD
      use MLSCommon,               only: r8
      use MLSMessageModule,        only: MLSMessage, MLSMSG_Error
      use MLSNumerics,             only: INTERPOLATEVALUES
      use ModelInput,              only: MODEL_ATMOS
      use ModelOutput,             only: SENSITIVITY
      use RadiativeTransferModule, only: RADXFER
      use ScatteringAngle,         only: ANGLE
      use SpectroscopyCatalog_m,   only: CATALOG_T
      use Tmp,                     only: GET_TAN_PRESS
      USE Intrinsic,               only: l_clear,                  &
                                      &  l_clear_0RH,              &
                                      &  l_clear_110RH_below_top,  &
                                      &  l_clear_lowest_0_110RH,   &
                                      &  l_cloudy_110RH_below_top, &
                                      &  l_cloudy_110RH_in_cloud,  &
                                      &  l_cloudy_nearside_only,   &
                                      &  l_none

      use Toggles, only: Emit, Levels, Switches, Toggle
      use Trace_M, only: Trace_begin, Trace_end

      IMPLICIT NONE
      private
      public :: CloudForwardModel

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------

contains 

  !----------------------------------------  CloudForwardModel  -----
      SUBROUTINE CloudForwardModel (doChannel, NF, NZ, NT, NS, N,      &
             &   NZmodel,                                              &
             &   FREQUENCY, PRESSURE, HEIGHT, TEMPERATURE, VMRin,      &
             &   WCin, IPSDin, ZZT, RE, ISURF, Cloud_der, I_Saturation, IFOV,   &
             &   Bill_data,                                            &
             &   h_obs, elev_offset, AntennaPattern,                   &
             &   TB0, DTcir, Trans, BETA, BETAc, Dm, TAUeff, SS,       &
             &   NU, NUA, NAB, NR, Slevl, noS, Catalog, LosVel )

!============================================================================C
!   >>>>>>>>> FULL CLOUD FORWARD MODEL FOR MICROWAVE LIMB SOUNDER >>>>>>>>   C
!----------------------------------------------------------------------------C
!                                                                            C
!     THIS PROGRAM IS USED IN EOS/AURA MLS LEVEL 2 DATA PROCESSING SOFTWARE  C
!     TO SIMULATE CLOUD INDUCED RADIANCES AND TO COMPUTE THE CLOUD RADIANCE  C
!     SENSITIVITY.  IT CAN ALSO BE USED IN CLEAR-SKY CONDITIONS AS A CLEAR-  C
!     SKY FORWARD MODEL.                                                     C
!                                                                            C
!     JONATHAN H. JIANG                                                      C
!     -- MAY 18, 2001: FIRST F77 WORKING VERSION.                            C
!     -- JUNE 9, 2001: FIRST F90 VERSION.                                    C
!     -- AUG  6, 2001: ADDED TRANS FUNCTION FOR CLOUD RETRIEVAL              C
!     -- AUG 18, 2001: ADDED FIELD OF VIEW AVERAGING                         C
!     -- SEP 19, 2001: MODIFIED F95 VERSION USE MODULES                      C
!     -- SEP 25, 2001: ADDED EFFECT OF ATMOSPHERIC REFRACTION                C 
!     -- Dec  1, 2001: ADDED OPTION TO USE BILL'S CLEAR SKY MODEL            C  
!----------------------------------------------------------------------------C
!                                                                            C
!     <<< INPUT PARAMETERS >>>                                               C
!     ------------------------                                               C
!                                                                            C
!     0. DEMENSION:                                                          C
!     -------------                                                          C
!     NF:                 -> Number of Frequencies.                          C
!     NZ:                 -> Number of Pressure Levels.                      C
!     NT:                 -> Number of Tangent Pressures.                    C
!     NS:                 -> Number of Chemical Species.                     C
!     N:                  -> Number of Cloud Species.                        C
!                                                                            C
!     1. ATMOSPHERIC PROFILES:                                               C
!     ------------------------                                               C
!     FREQUENCY  (NF)     -> Frequency (GHz).                                C 
!     PRESSURE   (NZ)     -> Pressure (hPa).                                 C
!     HEIGHT     (NZ)     -> Geopotential Height (hPa).                      C
!     TEMPERATURE(NZ)     -> Temperature (K).                                C 
!     VMRin   (NS,NZ)     -> NS=1: H2O Volumn Mixing Ratios (ppv).           C
!                            NS=2: O3 Volumn Mixing Ratio (ppv).             C
!                            NS=3: N2O Volumn Mixing Ratio (ppv).            C
!                            NS=4: HNO3 Volumn Mixing Ratio (ppv).           C
!                                                                            C
!     2. CLOUD PARAMETERS:                                                   C
!     --------------------                                                   C
!     WCin     (N,NZ)     -> N=1: Cloud Ice Water Content (g/m3).            C
!                            N=2: Cloud Liquid Water Content (g/m3).         C
!     IPSDin     (NZ)     -> Particle Size Distribution Flag.                C
!                                                                            C
!     3. OTHER PARAMETERS:                                                   C
!     --------------------                                                   C
!     ZT         (NT)     -> Tangent Pressure (hPa).                         C
!     RE                  -> Radius of Earth (m).                            C
!     ISURF               -> Surface Types (Default: 0).                     C
!     Cloud_der           -> Switch for Sensitivity Calculation.             C
!           none                                                             C
!           iwc_low_height                                                   C
!           iwc_high_height                                                  C
!     I_saturation        -> Cloud Control Switch.                           C
!           ICON = 0                if WCin==0 and i_saturation is cloudy    C
!           ICON = I_Saturation     otherwise                                C
!     -----------------------------------------------                        C
!                                                                            C
!     >>> OUTPUT PARAMETERS <<<                                              C
!     -------------------------                                              C
!                                                                            C
!     4. STANDARD OUTPUTS:                                                   C
!     --------------------                                                   C
!     TB0     (NT,NF)     -> Background Clear-sky Radiance (K).              C 
!     DRcir   (NT,NF)     -> Cloud Induced Radiance (K).                     C
!     TAUeff  (NT,NF)     -> Effective Cloud Optical Depth.                  C
!     SS      (NT,NF)     -> Cloud Radiance Sensitivity (K).                 C
!     BETA    (NZ,NF)     -> Total Extinction Profile (m-1).                 C
!     BETAc   (NZ,NF)     -> Cloud Extinction Profile (m-1).                 C
!     Trans   (noS,NT,NF) -> Clear Sky Transmittance Function                C
!     Dm      (N,NZ)      -> Mass-Mean-Diameters (micron).                   C 
!                            Note:1=Ice,2=Liquid.                            C 
!                                                                            C
!     -----------------------------------------------                        C
!                                                                            C
!     ((( INTERNAL MODEL PARAMETERS )))                                      C
!                                                                            C
!     5. INTERNAL MODEL PARAMETERS                                           C
!     ----------------------------                                           C
!     NZmodel = 640       -> Number of pressure levels.                      C
!     NU      = 16        -> Number of scattering angles.                    C
!     NUA     = 8         -> Number of azimuth angles.                       C
!     NAB     = 50        -> Maximum number of truncation terms.             C
!     NR      = 40        -> Number of particle size bins.                   C 
!     --------------------------------------------------------------         C
!                                                                            C
!     FREQUENCY RANGE: 1-3000GHz                                             C
!                                                                            C
!     VERSION: 1.0, MAY 18, 2001                                             C
!     MICROWAVE ATMOSPHERIC SCIENCE TEAM                                     C
!     JET PROPULSION LABORATORY                                              C
!     CALIFORNIA INSTITUTE OF TECHNOLOGY                                     C
!     4800 OAK GROVE DRIVE                                                   C
!     PASADENA, CA 91109-8099                                                C
!                                                                            C
!     EMAIL: JONATHAN@MLS.JPL.NASA.GOV                                       C
!     PHONE: (818) 354-7135                                                  C
!     FAX:   (818) 393-5065                                                  C
!============================================================================C

!---------------------------------------
!     INPUT PARAMETERS (INPUTS FROM L2)        ! -- INTERFACE AEA -- ! 
!---------------------------------------
      INTEGER, intent(in) :: NF                ! NUMBER OF FREQUENCIES
      INTEGER, intent(in) :: NZ                ! NUMBER OF PRESSURE LEVELS
      INTEGER, intent(in) :: NT                ! NUMBER OF TANGENT HEIGHTS
      INTEGER, intent(in) :: NS                ! NUMBER OF CHEMICAL SPECIES
      INTEGER, intent(in) :: N                 ! NUMBER OF CLOUD SPECIES
      INTEGER, intent(in) :: NZmodel           ! NUMBER OF INTERNAL MODEL LEVELS
      INTEGER, intent(in) :: noS               ! NUMBER OF S GRID

      INTEGER, intent(in) :: IPSDin(NZ)        ! SIZE-DISTRIBUTION FLAG
                                               ! IILL     (I=ICE, L=LIQUID)
                                               ! 1000:     I->MH, L->GAMMA 
                                               ! 1100:     I->LIU-CURRY
                                               ! 2000-3900:I->MODIFIED GAMMA
                                               !              WITH VARIOUS De,
                                               !              alpha
                                               ! 4000-5900:I->KNOLLENBERG WITH 
                                               !              VARIOUS b1
                                               ! 6000:     I->PSD FOR PSC

      INTEGER, intent(in) :: ISURF             ! SURFACE TYPE
                                               ! 0 = SIMPLE MODEL
                                               ! 1 = LAND
                                               ! 2 = SEA 

      INTEGER, intent(in) :: Cloud_der         ! SENSITIVITY SWITCH 

      INTEGER, intent(in) :: I_Saturation      ! Clear-Cloudy-sky CONTROL SWITCH 

      LOGICAL, intent(in) :: IFOV              ! FIELD OF VIEW AVERAGING SWITCH
                                               ! False = OFF
                                               ! True = ON

      REAL(r8), intent(in) :: FREQUENCY(NF)    ! FREQUENCIES (GHz)
      REAL(r8), intent(in) :: PRESSURE(NZ)     ! PRESSURE LEVEL
      REAL(r8), intent(in) :: HEIGHT(NZ)       ! PRESSURE HEIGHT
      REAL(r8), intent(in) :: TEMPERATURE(NZ)  ! ATMOSPHERIC TEMPERATURE

      REAL(r8), intent(in) :: VMRin(NS,NZ)     ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3 VOLUME MIXING RATIO

      REAL(r8), intent(in) :: WCin(N,NZ)       ! CLOUD WATER CONTENT
                                               ! N=1: ICE; N=2: LIQUID
      REAL(r8), intent(in) :: RE               ! EARTH RADIUS
      REAL(r8), intent(in) :: Slevl(noS)       ! Sgrid levels
      REAL(r8), intent(in) :: LosVel           ! Line of sight velocity

      LOGICAL, intent(in) :: doChannel(NF)     ! do only true channels

!--------------------------------------
!     OUTPUT PARAMETERS (OUTPUT TO L2)        
!--------------------------------------

      REAL(r8) :: TB0(NT,NF)                   ! CLEAR-SKY TB AT ZT
      REAL(r8) :: DTcir(NT,NF)                 ! CLOUD-INDUCED RADIANCE
      REAL(r8) :: TAUeff(NT,NF)                ! CLOUD EFFECTIVE OPTICAL DEPTH
      REAL(r8) :: SS(NT,NF)                    ! CLOUD RADIANCE SENSITIVITY
                                               ! (NT+1) FOR ZENITH LOOKING) 

      REAL(r8) :: Trans(noS,NT,NF)             ! Clear Trans Func
      REAL(r8) :: BETA(NZ,NF)                  ! TOTAL OPTICAL DEPTH
      REAL(r8) :: BETAc(NZ,NF)                 ! CLOUDY OPTICAL DEPTH
 
      REAL(r8) :: Dm(N,NZ)                     ! MASS-MEAN-DIAMETER

                                               ! -- END OF INTERFACE AREA -- !

!-----------------------------------------------------------------------------
                                               
!-------------------------------
!     INTERNAL MODEL PARAMETERS                ! -- INTERNAL AREA -- !
!-------------------------------

      INTEGER :: MULTI
      integer, parameter :: Nsub=5             !Below surface tangent grids
      
      INTEGER :: NU                            ! NO. OF SCATTERING ANGLES
      INTEGER :: NUA                           ! NO. OF SCAT. AZIMUTH ANGLES
      INTEGER :: NIWC                          ! NO. OF ICE WATER CONTENTS

      INTEGER :: NAB                           ! MAX. NO. OF 
      INTEGER :: NR                            ! NUMBER OF SIZE BINS
      INTEGER :: NABR(NR)                      ! TRUNCATION FOR A, B

      INTEGER :: LORS                          ! SURFACE TYPE: LAND OR SEA

      REAL(r8) :: U(NU)                        ! COSINE OF SCATTERING ANGLES
      REAL(r8) :: DU(NU)                       ! DELTA U
      REAL(r8) :: UA(NUA)                      ! COSINE OF SCAT AZIMUTH ANGLES
      REAL(r8) :: THETA(NU)                    ! SCATTERING ANGLES
      REAL(r8) :: PHI(NUA)                     ! SCATTERING AZIMUTH ANGLES
      REAL(r8) :: UI(NU,NU,NUA)                ! COSINE OF INCIDENT TB ANGLES
      REAL(r8) :: THETAI(NU,NU,NUA)            ! ANGLES FOR INCIDENT TB
      REAL(r8) :: PHH(N,NU,NZmodel-1)          ! PHASE FUNCTION 
      REAL(r8) :: tau_wetAll(NZmodel-1)            ! TOTAL OPTICAL DEPTH for 100%RH
      REAL(r8) :: tau_dry(NZmodel-1)            ! TOTAL OPTICAL DEPTH for 0%RH
      REAL(r8) :: W0(N,NZmodel-1)              ! SINGLE SCATTERING ALBEDO
      
      REAL(r8) :: S                            ! SALINITY
      REAL(r8) :: SWIND                        ! SEA SURFACE WIND
      REAL(r8) :: TS                           ! SURFACE TEMPERATURE (K)
      REAL(r8) :: RS(NU/2)                     ! SURFACE REFLECTIVITY
      REAL(r8) :: RC0(3)                       ! GAS ABS.SCAT.EXT COEFFS.
      REAL(r8) :: RC_TOT(3)                    ! TOTAL ABS.SCAT.EXT COEFFS.
      REAL(r8) :: Z(NZmodel-1)                 ! MODEL LAYER THICKNESS (m)
      REAL(r8) :: TAU0(NZmodel-1)              ! CLEAR-SKY DEPTH
      REAL(r8) :: TEMP(NZmodel-1)              ! MEAN LAYER TEMPERATURE (K)

      REAL(r8) :: TT(NZmodel/8+Nsub,NZmodel)   ! CLOUDY-SKY TB AT TANGENT 
                                               ! HEIGHT ZT (LAST INDEX FOR 
                                               ! ZENITH LOOKING)

      REAL(r8) :: TT0(NZmodel/8+Nsub,NZmodel)  ! CLEAR-SKY TB AT TANGENT
                                               ! HEIGHT ZT

!---------------------------------------------
!     INTERNAL ATMOSPHERIC PROFILE PARAMETERS
!---------------------------------------------

      INTEGER :: IPSD (NZmodel)

      REAL(r8) :: WC (N,NZmodel)
      REAL(r8) :: YP (NZmodel)
      REAL(r8) :: YZ (NZmodel)
      REAL(r8) :: YT (NZmodel)
      REAL(r8) :: YQ (NZmodel)                 ! H2O VOLUME MIXING RATIO
      REAL(r8) :: VMR(NS-1,NZmodel)            ! 1=O3 VOLUME MIXING RATIO
                                               ! 2=N2O VOLUME MIXING RATIO
                                               ! 3=HNO3 Volumn Mixing Ratio
      REAL(r8) :: DDm(N,NZmodel-1)                         

!----------------------------
!     CLOUD MODEL PARAMETERS
!----------------------------

      REAL(r8) :: CWC                          ! CLOUD WATER CONTENT (g/m3)
      REAL(r8) :: CDEPTH(N)                    ! CLOUD OPTICAL DEPTH
      REAL(r8) :: DEPTH                        ! TOTAL OPTICAL DEPTH
      REAL(r8) :: tau(NZmodel-1)            ! TOTAL EXTINCTION
                                               ! (CLEAR+CLOUD)
      REAL(r8) :: tau_clear(NZmodel-1)      ! EXTINCTION (clear part)
      REAL(r8) :: tau_wetCld(NZmodel-1)     ! EXTINCTION (clear part)
                                               ! assuming 110%RHi
                                               ! below 100mb or cloud top
      REAL(r8) :: tau_fewCld(NZmodel-1)     ! EXTINCTION (clear part)
                                               ! assuming 110%RHi
                                               ! only inside cloud
      REAL(r8) :: tauc(NZmodel-1)           ! EXTINCTION (cloud part)
      
!---------------------------
!     WORK SPACE PARAMETERS
!---------------------------

      INTEGER :: I
      INTEGER :: J
      INTEGER :: K
      INTEGER :: IFR 
      INTEGER :: ILYR
      INTEGER :: IL
      INTEGER :: ISPI
      INTEGER :: IIWC
      INTEGER :: ICLD_TOP                      ! cloud top index
      INTEGER :: I100_TOP                      ! 100 mb index
      INTEGER :: n1, n2

      REAL(r8) :: DMA
      REAL(r8) :: RATIO
      REAL(r8) :: PH0(N,NU,NZmodel-1)
      REAL(r8) :: W00(N,NZmodel-1)  
      REAL(r8) :: P11(NU)
      REAL(r8) :: RC11(3)
      REAL(r8) :: RC_TMP(N,3)
      REAL(r8) :: CHK_CLD(NZmodel)                        

      REAL(r8) :: ZZT(NT)                      ! TANGENT HEIGHTS (meters)

      REAL(r8) :: ZZT1(NZmodel/8-1+Nsub)       ! TANGENT HEIGHTS (meters) FOR 
                                               ! (a subset OF YZ)     | TO ZZT
                                               ! THE RESULT WILL BE INTERPOLATED

      REAL(r8) :: ZPT1(NZmodel/8-1+Nsub)       ! TANGENT PRESSURE (mb) FOR 
                                               ! (a subset OF YP)
                                               ! THIS IS ASSOCIATED WITH ZZT1

      REAL(r8) :: ZTT1(NZmodel/8-1+Nsub)       ! TEMPERATURE CORRESPONDING TO 
                                               ! (a subset OF YT)
                                               ! THIS IS ASSOCIATED WITH ZZT1

      REAL(r8) :: ZVT1(NZmodel/8-1+Nsub) 
      REAL(r8) :: ZNT1(NZmodel/8-1+Nsub) 

      REAL(r8) :: ptg_angle(NZmodel/8-1+Nsub)  ! POINTING ANGLES CORRESPONDING 

      type(antennaPattern_T), intent(in) :: AntennaPattern
      Type(Catalog_T), INTENT(IN) :: Catalog(:)

      integer :: FFT_pts
      REAL(r8) :: schi

      Real(r8), dimension( NZmodel/8-1+Nsub ) :: RAD0, RAD, SRad0, SRad

      REAL(r8) :: PH1(NU)                      ! SINGLE PARTICLE PHASE FUNCTION
      REAL(r8) :: P(NAB,NU)                    ! LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: DP(NAB,NU)                   ! Delt LEGENDRE POLYNOMIALS l=1
      
      REAL(r8) :: R(NR)                        ! PARTICLE RADIUS
      REAL(r8) :: RN(NR)                       ! NUMBER OF PARTICLES IN EACH BIN
      REAL(r8) :: BC(3,NR)                     ! SINGLE PARTICLE ABS/SCAT/EXT 
                                               ! COEFFS
      REAL(r8) :: DZ(NZ-1)

      REAL(r8) :: RT

      REAL(r8) :: h_obs, Rs_eq, elev_offset

      ! Constants used in formula for antenna smearing
      REAL(r8), PARAMETER :: CONST1 = 0.0000776_r8
      REAL(r8), PARAMETER :: CONST2 = 4810.0_r8

      Logical :: Bill_data
      Logical :: isClear                        ! clear-sky indicator
      Logical :: doSensitivity                  ! do Sensitivity ?

      COMPLEX(r8) A(NR,NAB),B(NR,NAB)          ! MIE COEFFICIENCIES

!-------------------------------------------------------------------------

!---------------<<<<<<<<<<<<< START EXCUTION >>>>>>>>>>>>-----------------

      if ( toggle(emit) ) call trace_begin ( 'Cloudy Sky Forward Model' )

!      CALL HEADER(1)

! initialization of TB0, DTcir, Trans, BETA, BETAc, Dm, TAUeff, SS

      tt0 = 0._r8
      tt  = 0._r8
      Tb0 = 0._r8
      TAU0 = 0.0_r8
      TAU=0.0_r8
      tau_wetAll=0.0_r8
      DTcir = 0._r8
      trans = 0._r8
      Betac = 0._r8
      Beta = 0._r8
      Dm = 0._r8
      Taueff = 0._r8
      SS = 0._r8
      fft_pts =0_r8
      RS=0.0_r8
      RC0=0.0_r8
      RC_TMP=0.0_r8
      RC_tot=0.0_r8
      chk_cld = 0._r8
         
      !Save time if there is no cloud in cloudy-sky requests
      isClear = .false.
      if(maxval(WCin) == 0) isClear = .true.
      doSensitivity = .true.
      if(cloud_der == l_none) doSensitivity = .false.

!print*,maxval(WCin)
!=========================================================================
!                    >>>>>> CHECK MODEL-INPUT <<<<<<< 
!-------------------------------------------------------------------------
!     CHECK IF THE INPUT PROFILE MATCHS THE MODEL INTERNAL GRID. IF NOT, 
!     INTERPOLATE IT TO INTERNAL GRID. ALSO SET TANGENT PRESSURE (hPa) 
!     TO TANGENT HEIGHT (km)
!=========================================================================

      if ( toggle(emit) .and. levels(emit) > 4 ) then
           call trace_begin ( 'Build Internal Grids' )
      endif

      CALL MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMRin,NZ,NS,N,   &
           &           WCin,IPSDin,YP,YZ,YT,YQ,VMR,WC,NZmodel,CHK_CLD,IPSD)
 
!     TAKE FIRST NT ELEMENT OF YZ AS TANGENT Z

      IF (NZmodel-1 .LE. NT) THEN
         call MLSMessage ( MLSMSG_Error, ModuleName, &
           &   'Too many tangent heights' )
      ENDIF

! ----------------------------------
! DEFINE INTERNAL TANGENT PRESSURES 
! ----------------------------------

      MULTI=NZmodel/8-1+Nsub

      n1=0
      Do i=1,NZmodel
      IF(YZ(i) .le. 28000.0_r8) n1=n1+1
      enddo
      n1=n1*2/(multi-nsub)
      n2=0
      Do i=1,NZmodel
      IF(YZ(i) .gt. 28000.0_r8) n2=n2+1
      enddo
      n2=n2*2/(multi-nsub)
     
! --The following sets first 5 tangent heights below the surface

      DO I=1, Multi
         if (i .le. Nsub) zzt1(i) = -20000.0_r8 + i*2000.0_r8 
         if (i .gt. Nsub) then
            if(i .le. (multi-Nsub)/2+Nsub) zzt1(i) = YZ( (I-Nsub)*n1-n1+1 )
            if(i .gt. (multi-Nsub)/2+Nsub) zzt1(i) = &
          &                   YZ( (I-Nsub)*n2-n2+1+((multi-Nsub)/2)*(n1-n2) )
         endif

         ZPT1(I) = 0._r8    
         ZTT1(I) = 0._r8
         ZVT1(I) = 0._r8                

      ENDDO

      IF (IFOV) THEN         
         CALL GET_TAN_PRESS ( YP, YZ, YT, YQ, NZmodel,        &
                ZPT1, ZZT1, ZTT1, ZVT1, Multi )

      ENDIF

!-----------------------------------------------
!     INITIALIZE SCATTERING AND INCIDENT ANGLES 
!-----------------------------------------------

      CALL ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

!----------------------------------
!     SURFACE INFORMATION
!----------------------------------

      LORS  = ISURF          

      TS    = 288._r8
      S     = 35._r8
      SWIND = 0._r8

      RATIO=1._r8
      NIWC=1                ! FAST SENSITIVITY CALCULATION

!------------------------------------------
!     PERFORM FULL SENSITIVITY CALCULATION
!------------------------------------------

      iwc_loop: do IIWC=1, NIWC
!=========================================================================
!                   >>>>>>> CLEAR-SKY MODULE <<<<<<<<
!-------------------------------------------------------------------------
!     COMPUTE CLEAR-SKY ABSORPTION COEFFICIENTS, INCLUDING DRY, WET 
!     CONTINUMA AND LINE EMISSIONS.   
!=========================================================================
      if ( toggle(emit) .and. levels(emit) > 4 ) then
         call Trace_End ( 'Internal Grids' )
      end if


      frequency_loop: do IFR=1, NF

      IF ( doChannel(IFR) ) then
               
         if ( toggle(emit) .and. levels(emit) > 4 ) then
            call Trace_Begin ( 'ClearSky' )
         end if

         CALL CLEAR_SKY(NZmodel-1,NU,TS,S,LORS,SWIND,                   &
              &         YZ,YP,YT,YQ,VMR,NS,                             &
              &         FREQUENCY(IFR),RS,U,TEMP,Z,TAU0,tau_wetAll,     &
              &         tau_dry,Catalog, Bill_data, LosVel, i_saturation ) 
!         CALL HEADER(3)

       ! find various layer indices         
         ICLD_TOP = 0
         I100_TOP = 0
         tau_wetCld=TAU0                       ! Initialize to TAU0
         tau_fewCld = tau0

         DO IL=1, NZmodel-1
            IF(CHK_CLD(IL) .NE. 0.)THEN
               ICLD_TOP=IL                     ! FIND INDEX FOR CLOUD-TOP 
               tau_fewCld(IL)=tau_wetAll(IL)   ! 100% SATURATION INSIDE CLOUD 
            ENDIF
            IF(YP(IL) .GE. 100._r8) THEN       ! IF BELOW 100MB
               I100_TOP=IL                     ! FIND INDEX FOR 100MB
            ENDIF
         ENDDO

         DO IL=1,MAX(ICLD_TOP,I100_TOP)        
              tau_wetCld(IL)=tau_wetAll(IL)    ! 110%rhi BELOW cloud top or 100mb
         ENDDO                              
         ! tau_wetCld will be used for transmission function calculation

         ! the following stuffs are only used in cloudy cases
            PHH      = 0._r8     ! phase function
            DDm      = 0._r8     ! mass-mean diameter
            W0       = 0._r8     ! single scattering albedo
            PH0      = 0._r8 
            W00      = 0._r8
            tau_clear = 0._r8

            if ( toggle(emit) .and. levels(emit) > 4 ) then
               call Trace_End ( 'ClearSky' )
            end if

       IF (.not. isClear) then 
      
            select case (i_saturation)
            case (l_cloudy_110RH_below_top) 
               tau_clear = tau_wetCld
            case (l_cloudy_110RH_in_cloud) 
               tau_clear = tau_fewCld
            case default
               tau_clear = tau_fewCld
            end select
               
         model_layer_loop: do ILYR=1, NZmodel-1
 
            RC0(1)=tau_clear(ILYR)/Z(ILYR)  ! GAS ABSORPTION COEFFICIENT
            RC0(2)=0._r8
            RC0(3)=RC0(1)                   ! CLEAR-SKY EXTINCTION COEFFICIENT

            RC_TOT =RC0                     ! initialized to clear-sky coeffs.
            RC_TMP =0._r8                   ! tmp for cloud coeffs
            cdepth = 0._r8

            DO ISPI=1,N
               CWC = RATIO*WC(ISPI,ILYR)
               IF(CWC .ge. 1.e-9_r8) then           
              
               if ( toggle(emit) .and. levels(emit) > 4 ) then
                  call Trace_Begin ( 'Clouds' )
               end if

                  CALL CLOUDY_SKY ( ISPI,CWC,TEMP(ILYR),FREQUENCY(IFR),  &
                       &          NU,U,DU,P11,RC11,IPSD(ILYR),DMA,       &
                       &          PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR)

                  PHH(ISPI,:,ILYR)=P11            ! INTERGRATED PHASE FUNCTION
                  CDEPTH(ISPI)=RC11(3)*Z(ILYR)
                  RC_TMP(ISPI,:)=RC11             ! VOLUME EXT/SCAT/ABS COEFFS
                  DDm(ISPI,ILYR)=DMA              ! MASS-MEAN-DIAMETER

               if ( toggle(emit) .and. levels(emit) > 4 ) then
                  call Trace_End ( 'Clouds' )
               end if

               ENDIF
                              
               DO J=1,3                           ! ADD cloud COEFFICIENTS
                  RC_TOT(J)=RC_TOT(J)+RC_TMP(ISPI,J)
               ENDDO
            
            ENDDO

            DO ISPI=1,N
               W0(ISPI,ILYR) = min( 1._r8, RC_TMP(ISPI,2)/RC_TOT(3) )
            ENDDO

            DEPTH=RC_TOT(3)*Z(ILYR)

            tau(ILYR) = max(0._r8, DEPTH )
            tauc(ILYR)= max(0._r8, CDEPTH(1) )     ! only consider scattering coeff
         enddo model_layer_loop
       Endif

! call radiative transfer calculations
      
       if ( toggle(emit) .and. levels(emit) > 4 ) then
            call Trace_Begin ( 'Radiative Transfer' )
       end if

       select case ( i_saturation )
       case ( l_clear )
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
       case ( l_clear_110RH_below_top )
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_wetCld,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
       case ( l_clear_0RH )
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_dry,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
       case ( l_clear_lowest_0_110RH )
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_dry,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_wetCld,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT,0,RE)                          !CLEAR-SKY
         TT0 = min(TT, TT0)
         TT  = TT0	   ! so that dTcir=0

       case (l_cloudy_110RH_in_cloud)
         ! this is for the sids/retrieval calculations (saturation only inside clouds)
         ! wc is supplied by model if it's the retrieval mode
         ! wc is supplied from IWC file if it's the sids mode
         IF (isClear) THEN
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
         ELSE
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_clear,&
              &   RS,TS,FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &   UI,UA,TT0,0,RE)                          !CLEAR-SKY
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PHH,MULTI,ZZT1,W0,TAU,RS,TS,&
              &  FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,         &
              &  UI,UA,TT,1,RE)                            !CLOUDY-SKY
         ENDIF
       case (l_cloudy_110RH_below_top)
         ! this is for the sids/retrieval calculations (saturation below cloud top)
         ! wc is supplied by model if it's the retrieval mode
         ! wc is supplied from IWC file if it's the sids mode
         IF (isClear) THEN
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
         ELSE
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_clear,&
              &   RS,TS,FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &   UI,UA,TT0,0,RE)                          !CLEAR-SKY
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PHH,MULTI,ZZT1,W0,TAU,RS,TS,&
              &  FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,         &
              &  UI,UA,TT,1,RE)                            !CLOUDY-SKY
         ENDIF
       case (l_cloudy_nearside_only)
         ! this is for the sids/retrieval calculations (cloud on the nearside of the tgt pt)
         ! wc is supplied by model if it's the retrieval mode
         ! wc is supplied from IWC file if it's the sids mode
         IF (isClear) THEN
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT  = TT0	   ! so that dTcir=0
         ELSE
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,tau_clear,&
              &   RS,TS,FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &   UI,UA,TT0,0,RE)                          !CLEAR-SKY
         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PHH,MULTI,ZZT1,W0,TAU,RS,TS,&
              &  FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,         &
              &  UI,UA,TT,3,RE)                            !CLOUDY-SKY
         ENDIF
       case default
         call MLSMessage(MLSMSG_Error, ModuleName,'invalid i_saturation')
       end select
 
       if ( toggle(emit) .and. levels(emit) > 4 ) then
          call Trace_End ( 'Radiative Transfer' )
       end if


!========================================================================
!     FOV case
!========================================================================
       IF (IFOV) THEN       ! **** BEGIN FOV AVERAGING ****

       if ( toggle(emit) .and. levels(emit) > 4 ) then
          call Trace_Begin ( 'FOV Averaging' )
       end if

!----------------------------------------------------------------------------
!    >>>>>> ADDS THE EFFECTS OF ANTENNA SMEARING TO THE RADIANCE <<<<<<
!----------------------------------------------------------------------------
!        Compute the effect of air refraction. The following method is based 
!        on the discussion with Bill Read, personal communication (09/25/01)    
!        (1+znt1) is air refractive index.   
!
!        NOTE IF TANGENT HEIGHT ZZT1 IS BELOW SURFACE, USE SURFACE VALUE OF
!        ZPT1, ZTT1 AND ZVT1 TO COMPUTE REFRACTIVE INDEX   
!----------------------------------------------------------------------

         K=0
         DO I=1, Multi
            IF (ZZT1(I) .LE. 0._r8) K=I
         END DO

         DO I=1, Multi                      
            IF (ZZT1(I) .LE. 0._r8) THEN
               znt1(I) = const1 * zpt1(K+1)/ztt1(K+1)
               znt1(I) = znt1(I)*(1.0_r8 + const2*zvt1(K+1)/ztt1(K+1)) 
            ELSE IF (ZZT1(I) .GT. 0._r8) THEN
               znt1(I) = const1 * zpt1(I)/ztt1(I)
               znt1(I) = znt1(I)*(1.0_r8 + const2*zvt1(I)/ztt1(I))
            END IF
         END DO
    
!---------------------------------------------------------------------------
!	 FIRST COMPUTE THE POINTING ANGLES (ptg_angle) 
!---------------------------------------------------------------------------
!	 Rs_eq = h_obs + 38.9014 * Sin(2.0*(phi_tan - 51.6814 * deg2rad)) 

	 Rs_eq = h_obs

         schi = 0.0_r8
         RT = 0.0_r8
         ptg_angle = 0.0_r8

  	    DO I = 1, Multi
            
            If (ZZT1(I) .LT. 0._r8) then
               RT= ( ZZT1(I) + RE )
               schi = (1+znt1(I)) * ( RT/RE ) * (ZZT1(I) + RE) / Rs_eq    
            else if (ZZT1(I) .GE. 0._r8) then
               schi = (1+znt1(I)) * (ZZT1(I) + RE) / Rs_eq 
            Endif

   	    IF(ABS(schi) > 1.0) THEN
              call MLSMessage ( MLSMSG_Error, ModuleName, &
                   & 'Pointing angle computation error: arg > 1.0 in ArcSin(arg)' )
            END IF
    	    ptg_angle(i) = Asin(schi) + elev_offset

  	    END DO

! ----------------------------------------------------------------
! 	 THEN DO THE FIELD OF VIEW AVERAGING
! ----------------------------------------------------------------

         RAD0=0.0_r8
         RAD0(1:Multi)=TT0(1:Multi,NZmodel)

         RAD=0.0_r8
         RAD(1:Multi)=TT(1:Multi,NZmodel)

         SRad0=0.0_r8
         SRad =0.0_r8

         call fov_convolve ( AntennaPattern, ptg_angle, RAD0, ptg_angle, SRad0 )
         call fov_convolve ( AntennaPattern, ptg_angle, RAD,  ptg_angle, SRad )

! -----------------------------------------------------------------------------

         ! CLEAR-SKY BACKGROUND
         CALL INTERPOLATEVALUES(ZZT1,SRad0(:),ZZT,TB0(:,IFR),method='Linear')

         ! CLOUD-INDUCED RADIANCE
         CALL INTERPOLATEVALUES(ZZT1,SRad(:)-SRad0(:),ZZT,DTcir(:,IFR), &
              &                 method='Linear')

         if ( toggle(emit) .and. levels(emit) > 4 ) then
              call Trace_End ( 'FOV Averaging' )
         end if

!========================================================================
!    non-FOV case
!========================================================================

       ELSE      

         ! CLEAR-SKY BACKGROUND
         CALL INTERPOLATEVALUES(ZZT1, TT0(:,NZmodel), ZZT, TB0(:,IFR), &
          & method='Linear')

         ! CLOUD-INDUCED RADIANCE
         CALL INTERPOLATEVALUES(ZZT1, TT(:,NZmodel)-TT0(:,NZmodel), ZZT, &
          & DTcir(:,IFR), &
              &                 method='Linear')

       ENDIF                            ! END OF FOV 

       IF(doSensitivity) THEN

         CALL SENSITIVITY (DTcir(:,IFR),ZZT,NT,YP,YZ,NZmodel,PRESSURE,NZ, &
              &      tau,tauc,tau_clear,TAUeff(:,IFR),SS(:,IFR), &
              &      Trans(:,:,IFR), BETA(:,IFR), BETAc(:,IFR), DDm, Dm, &
              &      Z, DZ, N,RE, noS, Slevl)     ! COMPUTE SENSITIVITY
              
              !for small dTcir at low tangent heights
               DO I=1,NT
               IF(DTcir(I,IFR) .LT. 0. .AND. DTcir(I,IFR) .GT. -3.) &
               & SS(I,IFR) = -50.
               ENDDO
       ENDIF

      END IF                                 ! END IF doCHANNEL
      enddo frequency_loop
      enddo iwc_loop                         ! 

!      CALL HEADER(5)

!-----------------------<<<<<<<<<<<<< END >>>>>>>>>>>>------------------------C

!      RETURN
!      END

      if ( toggle(emit) ) call trace_end ( 'Cloudy Sky Forward Model' )

      END SUBROUTINE CloudForwardModel

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module CloudySkyRadianceModel

! $Log$
! Revision 1.64  2004/01/08 00:40:50  jonathan
! fix bug
!
! Revision 1.63  2004/01/08 00:36:57  livesey
! Fixed some of the tracing
!
! Revision 1.62  2004/01/08 00:24:58  jonathan
! add tracing signals
!
! Revision 1.61  2003/07/22 23:49:52  jonathan
! change dimension for RAD0 etc to NZmodel/8-1+Nsub
!
! Revision 1.60  2003/05/14 23:12:29  dwu
! fix a dimension problem and modify the structure
!
! Revision 1.59  2003/05/07 23:05:00  jonathan
! some clean-up and cosmetic changes
!
! Revision 1.58  2003/04/10 20:26:21  dwu
! make i_saturation and cloud_der as a verbal argument
!
! Revision 1.57  2003/04/05 17:30:27  dwu
! make icon=0 for cloudy-sky request if there is no cloud
!
! Revision 1.56  2003/04/03 22:37:15  dwu
! fix a couple of bugs plus some cleanups
!
! Revision 1.55  2003/04/02 21:47:33  dwu
! initialize RC_TMP inside layer loop
!
! Revision 1.54  2003/04/02 20:00:04  dwu
! some clearup and replace ifov with do_conv
!
! Revision 1.53  2003/02/03 19:26:11  dwu
! fix a bug
!
! Revision 1.52  2003/02/03 19:25:00  dwu
! fix a bug
!
! Revision 1.51  2003/02/02 01:45:17  dwu
! some minor changes
!
! Revision 1.50  2003/01/31 23:17:15  dwu
! rename tau_wet to tau_wetAll
!
! Revision 1.49  2003/01/31 23:07:16  dwu
! mask 100%rhi only below 100mb
!
! Revision 1.48  2003/01/30 23:48:17  dwu
! add icon=-4 for min Tb and some clearups
!
! Revision 1.47  2003/01/23 00:19:09  pwagner
! Some cosmetic only (or so I hope) changes
!
! Revision 1.46  2003/01/16 18:39:41  pwagner
! Removed some unused variables
!
! Revision 1.45  2003/01/12 04:49:09  dwu
! change icon=4 to icon=-1 for better definition
!
! Revision 1.44  2003/01/12 03:52:22  dwu
! add ICON=4 for calculating clear-sky radiance limits assuming 100% RHi below 100hPa
!
! Revision 1.43  2003/01/09 21:08:23  dwu
! drop orbI
!
! Revision 1.42  2002/12/18 16:09:57  jonathan
! add orbI
!
! Revision 1.41  2002/12/09 17:34:04  jonathan
! *** empty log message ***
!
! Revision 1.40  2002/11/06 18:19:37  jonathan
! some minor changes
!
! Revision 1.39  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.38  2002/10/03 22:02:33  vsnyder
! Get Deg2Rad from Units instead of l2pc_file_parameters
!
! Revision 1.37  2002/09/27 22:41:09  jonathan
! change old FOV to new FOV routine
!
! Revision 1.36  2002/08/22 00:13:58  jonathan
! upgrade to include more molecules
!
! Revision 1.35  2002/08/19 22:22:03  jonathan
! debug stuff
!
! Revision 1.34  2002/08/09 22:21:42  jonathan
! new internal tangent heights
!
! Revision 1.33  2002/08/08 22:46:15  jonathan
! newly improved version
!
! Revision 1.32  2002/06/17 17:49:41  bill
! fix tangent height bug
!
! Revision 1.30  2001/11/16 00:41:08  jonathan
! add losVel
!
! Revision 1.29  2001/11/15 23:52:20  jonathan
! add default_spectroscopy
!
! Revision 1.28  2001/11/09 18:07:09  jonathan
! add spectra catalog
!
! Revision 1.27  2001/10/30 21:14:01  dwu
! change order of delTau100 and Tau0 initialization
!
! Revision 1.26  2001/10/30 05:55:35  dwu
! make clear sky dTcir =0
!
! Revision 1.25  2001/10/26 00:01:45  jonathan
! fixed a bug, change from w0*0._r8 etc to w00, ph0
!
! Revision 1.24  2001/10/25 16:10:37  dwu
! fix a bug
!
! Revision 1.23  2001/10/24 23:15:47  dwu
! some cleanup and correction
!
! Revision 1.22  2001/10/24 17:30:49  jonathan
! some minor changes
!
! Revision 1.21  2001/10/19 19:33:47  dwu
! change DDm dimension from NH to NH-1
!
! Revision 1.20  2001/10/19 19:30:36  dwu
! initialize output quantities
!
! Revision 1.19  2001/10/19 19:16:39  dwu
! change Nz-1 to NZ
!
! Revision 1.18  2001/10/12 22:40:07  dwu
! fix a bug in delTAU100
!
! Revision 1.17  2001/10/11 22:20:02  dwu
! make tangent height coming from outside
!
! Revision 1.16  2001/10/11 20:02:46  jonathan
! *** empty log message ***
!
! Revision 1.15  2001/10/09 22:12:10  jonathan
! *** empty log message ***
!
! Revision 1.14  2001/10/04 23:35:05  dwu
! *** empty log message ***
!
! Revision 1.13  2001/09/27 15:52:06  jonathan
! minor changes
!
! Revision 1.12  2001/09/26 17:04:37  jonathan
! added Air Refraction Correction, Jonathan
!
! Revision 1.11  2001/09/24 23:51:35  jonathan
! minor changes
!
! Revision 1.10  2001/09/21 15:51:37  jonathan
! modified F95 version
!

