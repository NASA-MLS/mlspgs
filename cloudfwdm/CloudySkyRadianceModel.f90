
! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
 
module CloudySkyRadianceModel

! -------------------------------------------------------------------------  
! MICROWAVE LIMB RADIATIVE TRANSFER MODEL FOR CLOUDY ATMOSPHERES  
! -------------------------------------------------------------------------
      use AntennaPatterns_m, only: AntennaPattern_T
      use ClearSkyModule, only: CLEAR_SKY
      use CloudySkyModule, only: CLOUDY_SKY
      use SpectroscopyCatalog_m, only: CATALOG_T
      use DCSPLINE_DER_M, only: CSPLINE_DER
      use FOV_CONVOLVE_M, only: FOV_CONVOLVE
      use HYDROSTATIC_INTRP, only: GET_PRESSURES
      use L2PC_FILE_PARAMETERS, only: DEG2RAD
      use ModelInput, only: MODEL_ATMOS
      use ModelOutput, only: SENSITIVITY
      use MLSCommon, only: r8
      use MLSNumerics, only: INTERPOLATEVALUES
      use PrtMsg, only: HEADER
      use RadiativeTransferModule, only: RADXFER
      use ScatteringAngle, only: ANGLE
      use Tmp, only: GET_TAN_PRESS

      IMPLICIT NONE
      private
      public :: CloudForwardModel

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------

contains 

      SUBROUTINE CloudForwardModel (doChannel, NF, NZ, NT, NS, N,      &
             &   NZmodel,                                              &
             &   FREQUENCY, PRESSURE, HEIGHT, TEMPERATURE, VMRin,      &
             &   WCin, IPSDin, ZT, ZZT, RE, ISURF, ISWI, ICON, IFOV,   &
             &   Bill_data,                                            &
             &   phi_tan, h_obs, elev_offset, AntennaPattern,          &
             &   TB0, DTcir, Trans, BETA, BETAc, Dm, TAUeff, SS,       &
             &   NU, NUA, NAB, NR, Slevl, noS, Catalog )

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
!----------------------------------------------------------------------------C
!                                                                            C
!     <<< INPUT PARAMETERS >>>                                               C
!     ------------------------                                               C
!                                                                            C
!     0. DEMENSION:                                                          C
!     -------------                                                          C
!     NF:             -> Number of Frequencies.                              C
!     NZ:             -> Number of Pressure Levels.                          C
!                        -------------------------------------------------   C
!               NOTE:    The internal model grid resuires NZmodel=640,       C
!                        where pressure levels are defined as following:     C
!                            NZmodel = 640                                   C
!                            Ptop    = 80./16.-3.                            C
!                            Pbottom =-ALOG10(PRESSURE(lowest_level))        C
!                            dP      =(Ptop-Pbottom)/NZmodel                 C
!                            ZH(I)   =Pbottom+(I-1)*dP,   I=1,NZmodel        C
!                            P(I)    = 10**(-ZH(I)) ,        I=1,NZmodel     C
!                        -------------------------------------------------   C
!     NT:             -> Number of Tangent Pressures.                        C
!     NS:             -> Number of Chemical Species.                         C
!     N:              -> Number of Cloud Species.                            C
!                                                                            C
!     1. ATMOSPHERIC PROFILES:                                               C
!     ------------------------                                               C
!     FREQUENCY  (NF) -> Frequency (GHz).                                    C 
!     PRESSURE   (NZ) -> Pressure (hPa).                                     C
!     HEIGHT     (NZ) -> Geopotential Height (hPa).                          C
!     TEMPERATURE(NZ) -> Temperature (K).                                    C 
!     VMRin   (NS,NZ) -> NS=1: H2O Volumn Mixing Ratios (ppm).               C
!                        NS=2: O3 Volumn Mixing Ratio (ppm).                 C
!                                                                            C
!     2. CLOUD PARAMETERS:                                                   C
!     --------------------                                                   C
!     WCin     (N,NZ) -> N=1: Cloud Ice Water Content (g/m3).                C
!                        N=2: Cloud Liquid Water Content (g/m3).             C
!     IPSDin     (NZ) -> Particle Size Distribution Flag.                    C
!                                                                            C
!     3. OTHER PARAMETERS:                                                   C
!     --------------------                                                   C
!     ZT         (NT) -> Tangent Pressure (hPa).                             C
!     RE              -> Radius of Earth (m).                                C
!     ISURE           -> Surface Types (Default: 0).                         C
!     ISWI            -> Switch for Sensitivity Calculation (Default: 0).    C
!     ICON            -> Cloud Control Switch (Default: 2).                  C
!     -----------------------------------------------                        C
!                                                                            C
!     >>> OUTPUT PARAMETERS <<<                                              C
!     -------------------------                                              C
!                                                                            C
!     4. STANDARD OUTPUTS:                                                   C
!     --------------------                                                   C
!     TB0     (NT,NF) -> Background Clear-sky Radiance (K).                  C 
!     DRcir   (NT,NF) -> Cloud Induced Radiance (K).                         C
!     TAUeff  (NT,NF) -> Effective Cloud Optical Depth.                      C
!     SS      (NT,NF) -> Cloud Radiance Sensitivity (K).                     C
!     BETA  (NZ,NF) -> Total Extinction Profile (m-1).                       C
!     BETAc (NZ,NF) -> Cloud Extinction Profile (m-1).                       C
!     Trans (noS,NT,NF) -> Clear Sky Transmittance Function                  C
!     Dm     (N,NZ) -> Mass-Mean-Diameters (micron).                         C 
!                        Note:1=Ice,2=Liquid.                                C 
!                                                                            C
!     -----------------------------------------------                        C
!                                                                            C
!     ((( INTERNAL MODEL PARAMETERS )))                                      C
!                                                                            C
!     5. INTERNAL MODEL PARAMETERS                                           C
!     ----------------------------                                           C
!     NZmodel = 640   -> Number of pressure levels.                          C
!     NU      = 16    -> Number of scattering angles.                        C
!     NUA     = 8     -> Number of azimuth angles.                           C
!     NAB     = 50    -> Maximum number of truncation terms.                 C
!     NR      = 40    -> Number of particle size bins.                       C 
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
      INTEGER :: NF                            ! NUMBER OF FREQUENCIES
      LOGICAL :: doChannel(NF)                 ! do only true channels
      INTEGER :: NZ                            ! NUMBER OF PRESSURE LEVELS
      INTEGER :: NT                            ! NUMBER OF TANGENT HEIGHTS
      INTEGER :: NS                            ! NUMBER OF CHEMICAL SPECIES
      INTEGER :: N                             ! NUMBER OF CLOUD SPECIES
      INTEGER :: NZmodel                       ! NUMBER OF INTERNAL MODEL LEVELS
      INTEGER :: MULTI
      INTEGER :: noS                           ! NUMBER OF S GRID

      INTEGER :: IPSDin(NZ)                    ! SIZE-DISTRIBUTION FLAG
                                               ! IILL     (I=ICE, L=LIQUID)
                                               ! 1000:     I->MH, L->GAMMA 
                                               ! 1100:     I->LIU-CURRY
                                               ! 2000-3900:I->MODIFIED GAMMA
                                               !              WITH VARIOUS De,
                                               !              alpha
                                               ! 4000-5900:I->KNOLLENBERG WITH 
                                               !              VARIOUS b1
                                               ! 6000:     I->PSD FOR PSC

      INTEGER :: ISURF                         ! SURFACE TYPE
                                               ! 0 = SIMPLE MODEL
                                               ! 1 = LAND
                                               ! 2 = SEA 
  
      INTEGER :: ISWI                          ! SENSITIVITY SWITCH
                                               ! 0 = OFF
                                               ! 1 = ON

      INTEGER :: ICON                          ! CONTROL SWITCH
                                               ! 0 = CLEAR-SKY
                                               ! 1 = 100% R.H. BELOW CLOUD
                                               ! 2 = DEFAULT
                                               ! 3 = NEAR SIDE CLOUD ONLY

      INTEGER :: IFOV                          ! FIELD OF VIEW AVERAGING SWITCH
                                               ! 0 = OFF
                                               ! 1 = ON

      REAL(r8) :: FREQUENCY(NF)                ! FREQUENCIES (GHz)
      REAL(r8) :: PRESSURE(NZ)                 ! PRESSURE LEVEL
      REAL(r8) :: HEIGHT(NZ)                   ! PRESSURE HEIGHT
      REAL(r8) :: TEMPERATURE(NZ)              ! ATMOSPHERIC TEMPERATURE
      REAL(r8) :: VMRin(NS,NZ)                 ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3 VOLUME MIXING RATIO

      REAL(r8) :: WCin(N,NZ)                   ! CLOUD WATER CONTENT
                                               ! N=1: ICE; N=2: LIQUID
      REAL(r8) :: ZT(NT)                       ! TANGENT PRESSURE
      REAL(r8) :: RE                           ! EARTH RADIUS
      REAL(r8) :: Slevl(noS)                   ! Sgrid levels

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

      REAL :: PI
      PARAMETER (PI=3.1415926)
      
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
      REAL(r8) :: TAU(NZmodel-1)               ! TOTAL OPTICAL DEPTH
      REAL(r8) :: TAU100(NZmodel-1)            ! TOTAL OPTICAL DEPTH AT 100%RH
      REAL(r8) :: W0(N,NZmodel-1)              ! SINGLE SCATTERING ALBEDO
      
      REAL(r8) :: S                            ! SALINITY
      REAL(r8) :: SWIND                        ! SEA SURFACE WIND
      REAL(r8) :: TS                           ! SURFACE TEMPERATURE (K)
      REAL(r8) :: RS(NU/2)                     ! SURFACE REFLECTIVITY
      REAL(r8) :: RC0(3)                       ! GAS ABS.SCAT.EXT COEFFS.
      REAL(r8) :: RC_TOT(3)                    ! TOTAL ABS.SCAT.EXT COEFFS.
      REAL(r8) :: Z(NZmodel-1)                 ! MODEL LAYER THICKNESS (m)
      REAL(r8) :: TAU0(NZmodel-1)              ! CLEAR-SKY OPTICAL DEPTH
      REAL(r8) :: TEMP(NZmodel-1)              ! MEAN LAYER TEMPERATURE (K)

      REAL(r8) :: TT(NZmodel/8,NZmodel)        ! CLOUDY-SKY TB AT TANGENT 
                                               ! HEIGHT ZT (LAST INDEX FOR 
                                               ! ZENITH LOOKING)
      REAL(r8) :: TT0(NZmodel/8,NZmodel)       ! CLEAR-SKY TB AT TANGENT
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
      REAL(r8) :: VMR(NS,NZmodel)              ! 1=O3 VOLUME MIXING RATIO
      REAL(r8) :: DDm(N,NZmodel-1)                         

!----------------------------
!     CLOUD MODEL PARAMETERS
!----------------------------

      REAL(r8) :: CWC                          ! CLOUD WATER CONTENT (g/m3)
      REAL(r8) :: CDEPTH(N)                    ! CLOUD OPTICAL DEPTH
      REAL(r8) :: DEPTH                        ! TOTAL OPTICAL DEPTH
                                               ! (CLEAR+CLOUD)

      REAL(r8) :: delTAU100(NZmodel-1)         ! TOTAL AIR EXTINCTION for 100%RHi
      REAL(r8) :: delTAU(NZmodel-1)            ! TOTAL EXTINCTION
      REAL(r8) :: delTAUc(NZmodel-1)           ! CLOUDY-SKY EXTINCTION
      
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
      INTEGER :: MY_NIWC
      INTEGER :: L

      REAL(r8) :: DMA
      REAL(r8) :: RATIO
      REAL(r8) :: PH0(N,NU,NZmodel-1)
      REAL(r8) :: W00(N,NZmodel-1)  
      REAL(r8) :: P11(NU)
      REAL(r8) :: RC11(3)
      REAL(r8) :: RC_TMP(N,3)
      REAL(r8) :: CHK_CLD(NZmodel)                        
      REAL(r8) :: ZZT(NT)                      ! TANGENT HEIGHTS (meters)

      REAL(r8) :: ZZT1(NZmodel/8-1)            ! TANGENT HEIGHTS (meters) FOR CALCULATION 
                                               ! (a subset OF YZ)
                                               ! THE RESULT WILL BE INTERPOLATED TO ZZT

      REAL(r8) :: ZPT1(NZmodel/8-1)            ! TANGENT PRESSURE (mb) FOR CALCULATION 
                                               ! (a subset OF YP)
                                               ! THIS IS ASSOCIATED WITH ZZT1

      REAL(r8) :: ZTT1(NZmodel/8-1)            ! TEMPERATURE CORRESPONDING TO TANGENT PRESSURE 
                                               ! (a subset OF YT)
                                               ! THIS IS ASSOCIATED WITH ZZT1

      REAL(r8) :: ZVT1(NZmodel/8-1) 
      REAL(r8) :: ZNT1(NZmodel/8-1) 

      REAL(r8) :: ptg_angle(NZmodel/8-1)       ! POINTING ANGLES CORRESPONDING TO ZZT1

      type(antennaPattern_T), intent(in) :: AntennaPattern
      Type(Catalog_T), INTENT(IN) :: Catalog(:)
      integer :: FFT_INDEX(size(antennaPattern%aaap))
      integer :: FFT_pts, ntr, ier, is, ktr, Ptg_i, brkpt
      REAL(r8) :: schi, center_angle, Q
      Real(r8) :: dTB0_dZT(NT,NF), dDTcir_dZT(NT,NF)

      Real(r8), dimension(size(fft_index)) :: FFT_ANGLES, FFT_PRESS, RAD0, RAD

      REAL(r8) :: PH1(NU)                      ! SINGLE PARTICLE PHASE FUNCTION
      REAL(r8) :: P(NAB,NU)                    ! LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: DP(NAB,NU)                   ! Delt LEGENDRE POLYNOMIALS l=1
      
      REAL(r8) :: R(NR)                        ! PARTICLE RADIUS
      REAL(r8) :: RN(NR)                       ! NUMBER OF PARTICLES IN EACH BIN
      REAL(r8) :: BC(3,NR)                     ! SINGLE PARTICLE ABS/SCAT/EXT 
                                               ! COEFFS
      REAL(r8) :: DZ(NZ-1)

      REAL(r8) :: RT

      REAL(r8) :: phi_tan, h_obs, Rs_eq, elev_offset, pp, ti, rr,freq

      REAL(r8), PARAMETER :: CONST1 = 0.0000776_r8
      REAL(r8), PARAMETER :: CONST2 = 4810.0_r8
      Logical :: Bill_data 

      COMPLEX(r8) A(NR,NAB),B(NR,NAB)          ! MIE COEFFICIENCIES

!-------------------------------------------------------------------------


!---------------<<<<<<<<<<<<< START EXCUTION >>>>>>>>>>>>-----------------

!      CALL HEADER(1)

! initialization of TB0, DTcir, Trans, BETA, BETAc, Dm, TAUeff, SS
      Tb0 = 0._r8
      DTcir = 0._r8
      trans = 0._r8
      Betac = 0._r8
      Beta = 0._r8
      Dm = 0._r8
      Taueff = 0._r8
      SS = 0._r8
      
!=========================================================================
!                    >>>>>> CHECK MODEL-INPUT <<<<<<< 
!-------------------------------------------------------------------------
!     CHECK IF THE INPUT PROFILE MATCHS THE MODEL INTERNAL GRID; 
!     SET TANGENT PRESSURE (hPa) TO TANGENT HEIGHT (km)
!=========================================================================

      CALL MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMRin,NZ,NS,N,   &
           &           WCin,IPSDin,                                 &
           &           YP,YZ,YT,YQ,VMR,WC,NZmodel,CHK_CLD,IPSD,     &
           &           ZT,ZZT,NT)
 
!     TAKE FIRST NT ELEMENT OF YZ AS TANGENT Z

      IF (NZmodel-1 .LE. NT) THEN
        PRINT*,'TOO MANY TANGENT HEIGHTS'
        STOP
      ENDIF

! ----------------------------------
! DEFINE INTERNAL TANGENT PRESSURES 
! ----------------------------------

      MULTI=NZmodel/8-1

!         DO I=1, Multi
!            ZZT1(I)=YZ(I*8-7)
!            ZPT1(I)=YP(I*8-7)
!            ZTT1(I)=YT(I*8-7)
!            ZVT1(I)=YQ(I*8-7)      
!         ENDDO
!
! ---------------The following sets first 9 tangent heights below zero

      DO I=1, Multi
         if (i .le. 10) zzt1(i) = -20000.0_r8 + i*2000.0_r8 
         if (i .gt. 10) zzt1(i) = 1000._r8*i-10000._r8
         ZPT1(I) = 0._r8    
         ZTT1(I) = 0._r8
         ZVT1(I) = 0._r8                
      ENDDO

      IF (IFOV .EQ. 1) THEN         
               CALL GET_TAN_PRESS ( YP, YZ, YT, YQ, NZmodel,        &
                    &               ZPT1, ZZT1, ZTT1, ZVT1, Multi )

      ENDIF

!      DO I=1, Multi
!       print*,  ZZT1(I), ZPT1(I), ZTT1(I), ZVT1(I)                
!      ENDDO
!      stop

!-----------------------------------------------
!     INITIALIZE SCATTERING AND INCIDENT ANGLES 
!-----------------------------------------------

      CALL ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

!----------------------------------
!     SURFACE INFORMATION
!----------------------------------

      LORS  = ISURF          
!      lors=0

      TS    = 288._r8
      S     = 35._r8
      SWIND = 0._r8
      NIWC  = 10

!      IF (ISWI .EQ. 0) THEN                  
         RATIO=1._r8
         MY_NIWC=1                ! SKIP FULL SENSITIVITY CALCULATION
!      ELSE
!          MY_NIWC=NIWC
!      ENDIF

!------------------------------------------
!     PERFORM FULL SENSITIVITY CALCULATION
!------------------------------------------

      DO 3000 IIWC=1, MY_NIWC    ! START OF IWC LOOP
!         IF (ISWI .EQ. 0) THEN
            RATIO=1._r8
!         ELSE
!            RATIO = 10.*(IIWC-1)**2*0.004+1.0E-9_r8
!         ENDIF

!=========================================================================
!                   >>>>>>> CLEAR-SKY MODULE <<<<<<<<
!-------------------------------------------------------------------------
!     COMPUTE CLEAR-SKY ABSORPTION COEFFICIENTS, INCLUDING DRY, WET 
!     CONTINUMA AND LINE EMISSIONS.   
!=========================================================================

      DO 2000 IFR=1, NF
      IF ( doChannel(IFR) ) then
        
         CALL CLEAR_SKY(NZmodel-1,NU,TS,S,LORS,SWIND,           &
              &         YZ,YP,YT,YQ,VMR,NS,                     &
              &         FREQUENCY(IFR),RS,U,TEMP,TAU0,Z,TAU100, &
              &         Catalog, Bill_data ) 

!         CALL HEADER(3)

!-----------------------------------------------------------------------------
!        ASSUME 100% SATURATION IN CLOUD LAYER
! 	 N.B.	ICON=0 is Clear-Sky only 
!		ICON=1 is for 100%RH inside and below Cloud
!		ICON=2 is default for 100%RH inside cloud ONLY
!		ICON=3 is for near-side cloud only 
!-----------------------------------------------------------------------------

         ICLD_TOP = 0
         I100_TOP = 0

         DO IL=1, NZmodel-1  
            IF(CHK_CLD(IL) .NE. 0.)THEN
               ICLD_TOP=IL                    ! FIND INDEX FOR CLOUD-TOP 
               TAU0(IL)=TAU100(IL)            ! 100% SATURATION INSIDE CLOUD 
            ENDIF
            IF(YP(IL) .GE. 100._r8) THEN      ! IF BELOW 100MB                    
               I100_TOP=IL                    ! FIND INDEX FOR 100MB
            ENDIF
         ENDDO

!        delTAU100 will be used for transmission function calculation
         delTAU100=TAU0                       ! Initialize to TAU0,                                         
         if (ICON .ne. 0) then                ! only for cloudy retrieval cases
           DO IL=1,MAX(ICLD_TOP,I100_TOP)        
              delTAU100(IL)=TAU100(IL)        ! MASK 100% SATURATION BELOW
           ENDDO                                ! 100MB or cloud top
         endif

         IF (ICON .EQ. 1) THEN
            DO IL=1, ICLD_TOP                 
               TAU0(IL)=TAU100(IL)            ! 100% SATURATION BELOW CLOUD
            ENDDO
         ENDIF
!----------------------------------------------------------------------------------

            PHH      = 0._r8     ! phase function
            DDm      = 0._r8     ! mass-mean diameter
            W0       = 0._r8     ! single scattering albedo

            PH0      = 0._r8 
            W00      = 0._r8 

         DO 1000 ILYR=1, NZmodel-1            ! START OF MODEL LAYER LOOP:   
 
            RC0(1)=TAU0(ILYR)/Z(ILYR)         ! GAS ABSORPTION COEFFICIENT
            RC0(2)=0._r8
            RC0(3)=RC0(1)                     ! CLEAR-SKY EXTINCTION COEFFICIENT

            RC_TOT   = 0._r8     ! total extinction
            RC11     = 0._r8     ! cloud abs. scat. ext. coefficients
            RC_TMP   = 0._r8     ! cloud abs. scat. ext. coefficients (for ice/water)
            CDEPTH   = 0._r8     ! cloud optical depth

            DEPTH  = 0._r8

            DO ISPI=1,N
            CWC = RATIO*WC(ISPI,ILYR)
            IF(CWC .ne. 0._r8 .and. ICON .ne. 0) then           
            CWC = MAX(1.E-9_r8,CWC)
              
!=================================================
!    >>>>>>>>> CLOUDY-SKY MODULE <<<<<<<<<<<
!=================================================

               CALL CLOUDY_SKY ( ISPI,CWC,TEMP(ILYR),FREQUENCY(IFR),  &
                       &          NU,U,DU,P11,RC11,IPSD(ILYR),DMA,    &
                       &          PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR)

               PHH(ISPI,:,ILYR)=P11          ! INTERGRATED PHASE FUNCTION
               CDEPTH(ISPI)=RC11(3)*Z(ILYR)
               RC_TMP(ISPI,:)=RC11        ! VOLUME EXT/SCAT/ABS COEFFS
               DDm(ISPI,ILYR)=DMA                     ! MASS-MEAN-DIAMETER
            ENDIF
            ENDDO
                              
            DO J=1,3                               ! ADD all COEFFICIENTS
               RC_TOT(J)=RC0(J)+RC_TMP(1,J)+RC_TMP(2,J)
            ENDDO

            DO ISPI=1,N
               W0(ISPI,ILYR)=RC_TMP(ISPI,2)/RC_TOT(3) ! SINGLE SCATTERING ALBEDO   
               IF(W0(ISPI,ILYR) .GT. 1.) THEN         
                  W0(ISPI,ILYR)=1.
               ENDIF
            ENDDO

            TAU(ILYR)=RC_TOT(3)*Z(ILYR)
            DEPTH=RC_TOT(3)*Z(ILYR)

            delTAU(ILYR) = max(0._r8, DEPTH )
            delTAUc(ILYR)= max(0._r8, CDEPTH(1) )

 1000    CONTINUE                         ! END OF MODEL LAYER LOOP

!==================================================
!    >>>>>>> RADIATIVE TRANSFER MODULE <<<<<<<<<<
!==================================================

!         CALL HEADER(4)

         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,MULTI,ZZT1,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,0,RE)                          !CLEAR-SKY
         TT = TT0	! so that dTcir=0

         IF(ICON .GE. 1) THEN                               

           CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PHH,MULTI,ZZT1,W0,TAU,RS,TS,&
                &  FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,         &
                &  UI,UA,TT,ICON,RE)                            !CLOUDY-SKY

         ENDIF


         IF (IFOV .EQ. 1) THEN       ! **** BEGIN FOV AVERAGING ****

! ==========================================================================
!    >>>>>> ADDS THE EFFECTS OF ANTENNA SMEARING TO THE RADIANCE <<<<<<
! ==========================================================================

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
            IF (ZZT1(I) .LT. 0._r8) K=I
         END DO

         DO I=1, Multi                      
            IF (ZZT1(I) .LT. 0._r8) THEN
               znt1(I) = const1 * zpt1(K+1)/ztt1(K+1)
               znt1(I) = znt1(I)*(1.0_r8 + const2*zvt1(K+1)/ztt1(K+1)) 
            ELSE IF (ZZT1(I) .GE. 0._r8) THEN
               znt1(I) = const1 * zpt1(I)/ztt1(I)
               znt1(I) = znt1(I)*(1.0_r8 + const2*zvt1(I)/ztt1(I))                
            END IF
         END DO
    
!---------------------------------------------------------------------------
!	 FIRST COMPUTE THE POINTING ANGLES (ptg_angle) 
!---------------------------------------------------------------------------
!	 Rs_eq = h_obs + 38.9014 * Sin(2.0*(phi_tan - 51.6814 * deg2rad)) 

	 Rs_eq = h_obs

  	 DO I = 1, Multi
            
            If (ZZT1(I) .LT. 0._r8) then
               RT= MIN( (ZZT1(I)+RE), RE)
               schi = (1+znt1(I)) * ( RT/RE ) * (ZZT1(I) + RE) / Rs_eq    
            else if (ZZT1(I) .GE. 0._r8) then
               schi = (1+znt1(I)) * (ZZT1(I) + RE) / Rs_eq 
            Endif

   	    IF(ABS(schi) > 1.0) THEN
      	       PRINT *,'*** ERROR IN COMPUTING POINTING ANGLES'
               PRINT *,'    arg > 1.0 in ArcSin(arg) ..'
               STOP
            END IF
    	    ptg_angle(i) = Asin(schi) + elev_offset
!            print*, real(ptg_angle(i)), real(TT0(i,NZmodel))
  	 END DO

	 center_angle = ptg_angle(1)

! ----------------------------------------------------------------
! 	 THEN DO THE FIELD OF VIEW AVERAGING
! ----------------------------------------------------------------

         Ier = 0
         ntr = size(antennaPattern%aaap)

	 fft_pts = nint(log(real(size(AntennaPattern%aaap)))/log(2.0))

         RAD0=0.0_r8
         RAD0(1:Multi)=TT0(1:Multi,NZmodel)

         RAD=0.0_r8
         RAD(1:Multi)=TT(1:Multi,NZmodel)

         fft_angles=0.0
         fft_angles(1:Multi) = ptg_angle(1:Multi)    ! Multi = No. of tangent heights
	
         Call fov_convolve ( fft_angles, RAD0, center_angle, 1, Multi,   &
              &              fft_pts, AntennaPattern, Ier )
              if ( Ier /= 0) then
	         print*,'error in FOV CONV'
	         stop
	      endif

         Call fov_convolve ( fft_angles, RAD, center_angle, 1, Multi,   &
              &              fft_pts, AntennaPattern, Ier )
              if ( Ier /= 0) then
	         print*,'error in FOV CONV'
	         stop
	      endif

              !------------------------------------------------------------------
              !  TT   ___                                     RAD  ___|___
              !          \                                        /   |   \
              !           \          after fov_convolve ==>      /    |    \
              !            \                                    /  -  | +   \
              !             \___                            ___/      |      \___
              !-------------------------------------------------------------------


         !  Get 'ntr' pressures associated with the fft_angles:
         Call get_pressures ( 'a', ptg_angle, ZTT1, -log10(ZPT1), Multi,     &
              &                    fft_angles, fft_press, Ntr, Ier )
              if ( Ier /= 0) then
	         print*,'error in get_pressures'
	         stop
	      endif

         ! Make sure the fft_press array is MONOTONICALY increasing:
         is = 1
         do while (is < Ntr-1  .and.  fft_press(is) >= fft_press(is+1)) 
            is = is + 1                                                  
         end do                                                         
          ! There was an error in zvi's code, in which he wrote is<Ntr. This will cause 
          ! subscript fft_press(is+1) out of range when run the program after compile 
          ! with -C, as "Subscript 1 of FFT_PRESS (value 1025) is out of range (1:1024)"
          ! UPDATE: zvi and I have fixed above error on 08/29/01

         Ktr = 1
         Rad0(Ktr) = Rad0(is)
         Rad(Ktr) = Rad(is)
         fft_index(Ktr) = is
         fft_press(Ktr) = fft_press(is)

         do ptg_i = is+1, Ntr
           q = fft_press(ptg_i)
           if ( q > fft_press(Ktr)) then
             Ktr = Ktr + 1
             fft_press(Ktr) = q
             Rad0(Ktr) = Rad0(ptg_i)
             Rad(Ktr) = Rad(ptg_i)
             fft_index(Ktr) = ptg_i
           end if
         end do
 
! ----------------------------------------------------------------------------
!        INTERPOLATE THE OUTPUT VALUES RAD0 TO TB0 WITH RESPECT TO L2 TANGENT
!        PRESSURES ZT (i.e. ptan).
!        (ALSO STORE THE RADIANCE DERIVATIVES WITH RESPECT TO ZT)
! ----------------------------------------------------------------------------

         Call Cspline_der ( fft_press, -log10(ZT), RAD0, TB0(:,IFR), dTB0_dZT(:,IFR), Ktr, NT )
         Call Cspline_der ( fft_press, -log10(ZT), RAD-RAD0, DTcir(:,IFR), dDTcir_dZT(:,IFR), Ktr, NT )

! -----------------------------------------------------------------------------
         END IF     ! **** END OF FOV AVERAGING ****

!====================================
!    >>>>>>> MODEL-OUTPUT <<<<<<<<<
!====================================

         IF (IFOV .EQ. 0) THEN       

         ! CLEAR-SKY BACKGROUND
         CALL INTERPOLATEVALUES(ZZT1,TT0(:,NZmodel),ZZT,TB0(:,IFR),method='Linear')

         ! CLOUD-INDUCED RADIANCE
         CALL INTERPOLATEVALUES(ZZT1,TT(:,NZmodel)-TT0(:,NZmodel),ZZT,DTcir(:,IFR), &
              &                 method='Linear')

         ENDIF

!         if (iswi == 0) &
         CALL SENSITIVITY (DTcir(:,IFR),ZZT,NT,YP,YZ,NZmodel,PRESSURE,NZ, &
              &            delTAU,delTAUc,delTAU100,TAUeff(:,IFR),SS(:,IFR), &
              &            Trans(:,:,IFR), BETA(:,IFR), BETAc(:,IFR), DDm, Dm, Z, DZ, &
              &            N,RE, noS, Slevl) ! COMPUTE SENSITIVITY

      END IF                                 ! END OF DO CHANNEL

 2000 CONTINUE                               ! END OF FREQUENCY LOOP   

 3000 CONTINUE                               ! END OF IWC LOOP

!      CALL HEADER(5)

!-----------------------<<<<<<<<<<<<< END >>>>>>>>>>>>------------------------C

!      RETURN
!      END

      END SUBROUTINE CloudForwardModel

end module CloudySkyRadianceModel

! $Log$
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

