      SUBROUTINE CloudForwardModel (NF, NZ, NT, NS, N, NZmodel,   &
             &   FREQUENCY, PRESSURE, HEIGHT, TEMPERATURE, VMRin, &
             &   WCin, IPSDin,                                    &
             &   ZT, RE, ISURF, ISWI, ICON,                       &
             &   TB0, DTcir, BETA, BETAc, Dm, TAUeff, SS,         &
             &   NU, NUA, NAB, NR)

!============================================================================C
!   >>>>>>>>> FULL CLOUD FORWARD MODEL FOR MICROWAVE LIMB SOUNDER >>>>>>>>   C
!----------------------------------------------------------------------------C
!                                                                            C
!     THIS PROGRAM IS USED IN LEVEL 2 DATA PROCESSING TO SIMULATE CLOUD      C
!     INDUCED RADIANCES AND CLOUD RADIANCE SENSITIVITY.  IT CAN ALSO BE      C
!     IN CLEAR-SKY CONDITION AS A CLEAR-SKY FORWARD MODEL, THAT WILL BE      C
!     IN CLOUD FLAGING PROCESS. IN THE CLOUD FLAGING CASE, THIS PROGRAM      C
!     MUST BE CALLED TWICE TO COMPUTE CLEAR-SKY RADIANCES IN BOTH DRY &      C
!     WET (100% RELATIVE HUMIDITY) CONDITIONS.                               C
!                                                                            C
!     JONATHAN H. JIANG                                                      C
!     -- MAY 18, 2001: FIRST WORKING VERSION.                                C
!     -- JUNE 9, 2001: ELIMINATE INTERNAL GRID SO THAT THE INPUT/OUTPUT      C
!                      GRID ARE THE SAME AS THE INTERNAL GRID. HOWEVER,      C
!                      THE NUMBER OF MODEL PARAMETER NEEDED TO PASSE TO      C
!                      LEVEL 2 WAS THEREFORE INCREASED.                      C
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
!               NOTE:    The internal model grid resuires NZ=640, and        C
!                        PRESSURE levels are defined as the following:       C
!                            NZ=640                                          C
!                            Ptop = 40./16.-3.                               C
!                            Pbottom=-ALOG10(PRESSURE(lowest-level))         C
!                            dP=(Ptop-Pbottom)/NZ                            C
!                            ZH(I)=Pbottom+(I-1)*dP,      I=1,NZ             C
!                            PRESSURE(I) = 10**(-ZH(I)) , I=1,NZ             C
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
!     BETA  (NZ-1,NF) -> Total Extinction Profile (m-1).                     C
!     BETAc (NZ-1,NF) -> Cloud Extinction Profile (m-1).                     C
!     Dm     (N,NZ-1) -> Mass-Mean-Diameters (micron). Note:1=Ice,2=Liquid.  C
!                                                                            C
!     -----------------------------------------------                        C
!                                                                            C
!     ((( INTERNAL MODEL PARAMETERS )))                                      C
!                                                                            C
!     5. INTERNAL MODEL PARAMETERS                                           C
!     ----------------------------                                           C
!     NU  = 16        -> Number of scattering angles.                        C
!     NUA = 8         -> Number of azimuth angles.                           C
!     NAB = 50        -> Maximum number of truncation terms.                 C
!     NR  = 40        -> Number of particle size bins.                       C 
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

      use MLSCommon, only: r8
      IMPLICIT NONE

!---------------------------------------
!     INPUT PARAMETERS (INPUTS FROM L2)        ! -- INTERFACE AEA -- ! 
!---------------------------------------

      INTEGER :: NF                            ! NUMBER OF FREQUENCIES
      INTEGER :: NZ                            ! NUMBER OF PRESSURE LEVELS
      INTEGER :: NT                            ! NUMBER OF TANGENT HEIGHTS
      INTEGER :: NS                            ! NUMBER OF CHEMICAL SPECIES
      INTEGER :: N                             ! NUMBER OF CLOUD SPECIES
      INTEGER :: NZmodel                       ! NUMBER OF INTERNAL MODEL LEVELS


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
                                               ! 1 = CLEAR-SKY, 100% R.H. 
                                               !     BELOW 100hPa
                                               ! 2 = DEFAULT
                                               ! 3 = NEAR SIDE CLOUD ONLY

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

!--------------------------------------
!     OUTPUT PARAMETERS (OUTPUT TO L2)        
!--------------------------------------

      REAL(r8) :: TB0(NT,NF)                   ! CLEAR-SKY TB AT ZT
      REAL(r8) :: DTcir(NT,NF)                 ! CLOUD-INDUCED RADIANCE

      REAL(r8) :: TAUeff(NT,NF)                ! CLOUD EFFECTIVE OPTICAL DEPTH
      REAL(r8) :: SS(NT,NF)                    ! CLOUD RADIANCE SENSITIVITY
                                               ! (NT+1) FOR ZENITH LOOKING) 

      REAL(r8) :: BETA(NZ-1,NF)                ! TOTAL OPTICAL DEPTH
      REAL(r8) :: BETAc(NZ-1,NF)               ! CLOUDY OPTICAL DEPTH

      REAL(r8) :: Dm(N,NZ-1)                   ! MASS-MEAN-DIAMETER

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
      REAL(r8) :: RC(N,3)                      ! CLOUD ABS.SCAT.EXT COEFFS.
      REAL(r8) :: RC_TOT(3)                    ! TOTAL ABS.SCAT.EXT COEFFS.
      REAL(r8) :: Z(NZmodel-1)                 ! MODEL LAYER THICKNESS (m)
      REAL(r8) :: TAU0(NZmodel-1)              ! CLEAR-SKY OPTICAL DEPTH
      REAL(r8) :: TEMP(NZmodel-1)              ! MEAN LAYER TEMPERATURE (K)

      REAL(r8) :: TT(NT+1,NZmodel)             ! CLOUDY-SKY TB AT TANGENT 
                                               ! HEIGHT ZT (LAST INDEX FOR 
                                               ! ZENITH LOOKING)
      REAL(r8) :: TT0(NT+1,NZmodel)            ! CLEAR-SKY TB AT TANGENT
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
      REAL(r8) :: DDm(N,NZmodel)                         

!----------------------------
!     CLOUD MODEL PARAMETERS
!----------------------------

      REAL(r8) :: CWC                          ! CLOUD WATER CONTENT (g/m3)
      REAL(r8) :: CDEPTH(N)                    ! CLOUD OPTICAL DEPTH
      REAL(r8) :: DEPTH                        ! TOTAL OPTICAL DEPTH
                                               ! (CLEAR+CLOUD)

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
      INTEGER :: ICLD_TOP
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
      REAL(r8) :: ZZT(NT)
      REAL(r8) :: PH1(NU)                      ! SINGLE PARTICLE PHASE FUNCTION
      REAL(r8) :: P(NAB,NU)                    ! LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: DP(NAB,NU)                   ! Delt LEGENDRE POLYNOMIALS l=1
      
      REAL(r8) :: R(NR)                        ! PARTICLE RADIUS
      REAL(r8) :: RN(NR)                       ! NUMBER OF PARTICLES IN EACH BIN
      REAL(r8) :: BC(3,NR)                     ! SINGLE PARTICLE ABS/SCAT/EXT 
                                               ! COEFFS

      REAL(r8) :: DZ(NZ-1)
      COMPLEX(r8) A(NR,NAB),B(NR,NAB)          ! MIE COEFFICIENCIES

!---------------<<<<<<<<<<<<< START EXCUTION >>>>>>>>>>>>-------------------C

      CALL HEADER(1)
      RE = 6370000._r8

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

      IF (ISWI .EQ. 0) THEN                  
         RATIO=1.
         MY_NIWC=1                ! SKIP FULL SENSITIVITY CALCULATION
      ELSE
          MY_NIWC=NIWC
      ENDIF

!------------------------------------------
!     PERFORM FULL SENSITIVITY CALCULATION
!------------------------------------------

      DO 3000 IIWC=1, MY_NIWC    ! START OF IWC LOOP
         IF (ISWI .EQ. 0) THEN
            RATIO=1.
         ELSE
            RATIO = 10.*(IIWC-1)**2*0.004+1.0E-9
         ENDIF

!=========================================================================
!                   >>>>>>> CLEAR-SKY MODULE <<<<<<<<
!-------------------------------------------------------------------------
!     COMPUTE CLEAR-SKY ABSORPTION COEFFICIENTS, INCLUDING DRY, WET 
!     CONTINUMA AND LINE EMISSIONS.   
!=========================================================================

       DO 2000 IFR=1, NF

         CALL CLEAR_SKY(NZmodel-1,NU,TS,S,LORS,SWIND,           &
              &         YZ,YP,YT,YQ,VMR,NS,                     &
              &         FREQUENCY(IFR),RS,U,TEMP,TAU0,Z,TAU100) 

         CALL HEADER(3)

!-----------------------------------------------------
!        ASSUME 100% SATURATION IN CLOUD LAYER
!-----------------------------------------------------

         DO IL=1, NZmodel-1                   ! 100% SATURATION INSIDE CLOUD 
            IF(CHK_CLD(IL) .NE. 0.)THEN
               ICLD_TOP=IL
               IF(YZ(IL) .LT. 20.)THEN
                  TAU0(IL)=TAU100(IL)
               ENDIF
            ENDIF
         ENDDO

         IF (ICON .EQ. 1) THEN
            DO IL=1,ICLD_TOP             ! 100% SATURATION BELOW CLOUD 
               TAU0(IL)=TAU100(IL)       
            ENDDO
         ENDIF   

!--------------------------------------------------------

         DO 1000 ILYR=1, NZmodel-1             ! START OF MODEL LAYER LOOP:   
 
            RC0(1)=TAU0(ILYR)/Z(ILYR)     ! GAS ABSORPTION COEFFICIENT
            RC0(2)=0.
            RC0(3)=RC0(1)                 ! CLEAR-SKY EXTINCTION COEFFICIENT

            DO J=1,3
               RC_TOT(J) = 0.
               RC11(J)=0.
            ENDDO

            DO ISPI=1,2
               DO J=1,3
                  RC(ISPI,J)=0.
                  RC_TMP(ISPI,J)=0.
               ENDDO
               CDEPTH(ISPI) = 0.
               DO K=1,NU
                  PHH(ISPI,K,ILYR) = 0.
                  PH0(ISPI,K,ILYR)=0.
               ENDDO
               DDm(ISPI,ILYR)=0.
               W0(ISPI,ILYR)=0.
               W00(ISPI,ILYR)=0.
            ENDDO

            DEPTH  = 0.
            CWC = 1.E-9_r8

            IF(CHK_CLD(ILYR) .NE. 0.) THEN 

               DO ISPI=1,N
                  CWC = RATIO*WC(ISPI,ILYR) 
                  CWC = MAX(1.E-9_r8,CWC)
              
!=================================================
!    >>>>>>>>> CLOUDY-SKY MODULE <<<<<<<<<<<
!=================================================
!                  write(11,*) 'cloud'
!                  write(11,*) ISPI,CWC,TEMP(ILYR),FREQUENCY(IFR)

                  CALL CLOUDY_SKY ( ISPI,CWC,TEMP(ILYR),FREQUENCY(IFR),  &
                       &          NU,U,DU,P11,RC11,IPSD(ILYR),DMA,       &
                       &          PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR)

!                  write(11,*) PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR

                  DO K=1,NU
                     PHH(ISPI,K,ILYR)=P11(K)      ! INTERGRATED PHASE FUNCTION
                  ENDDO
                  CDEPTH(ISPI)=RC11(3)*Z(ILYR)
                  
                  DO J=1,3
                     RC_TMP(ISPI,J)=RC11(J)       ! VOLUME EXT/SCAT/ABS COEFFS
                  ENDDO
               ENDDO
            ENDIF
                              
            DO J=1,3                    ! ADD CLEAR-SKY COEFFICIENTS
               RC_TOT(J)=RC0(J)+RC_TMP(1,J)+RC_TMP(2,J)
            ENDDO

            DO ISPI=1,N        
                                        
               W0(ISPI,ILYR)=RC_TMP(ISPI,2)/RC_TOT(3) ! SINGLE SCATTERING ALBEDO   
               IF(W0(ISPI,ILYR) .GT. 1.) THEN         
                  W0(ISPI,ILYR)=1.
               ENDIF
               DDm(ISPI,ILYR)=DMA                     ! MASS-MEAN-DIAMETER
            ENDDO

            TAU(ILYR)=RC_TOT(3)*Z(ILYR)
            DEPTH=RC_TOT(3)*Z(ILYR)

            delTAU(ILYR) = max(0._r8, DEPTH )
            delTAUc(ILYR)= max(0._r8, CDEPTH(1) )

 1000    CONTINUE                         ! END OF MODEL LAYER LOOP

!==================================================
!    >>>>>>> RADIATIVE TRANSFER MODULE <<<<<<<<<<
!==================================================

         CALL HEADER(4)

         CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PH0,NT,ZZT,W00,TAU0,RS,TS,&
              &     FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,        &
              &     UI,UA,TT0,NT,ICON,RE)                          !CLEAR-SKY

         IF(ICON .GT. 1) THEN                               

           CALL RADXFER(NZmodel-1,NU,NUA,U,DU,PHH,NT,ZZT,W0,TAU,RS,TS,&
                &  FREQUENCY(IFR),YZ,TEMP,N,THETA,THETAI,PHI,         &
                &  UI,UA,TT,NT,ICON,RE)                            !CLOUDY-SKY

         ENDIF

!====================================
!    >>>>>>> MODEL-OUTPUT <<<<<<<<<
!====================================

         DO I=1,NT
            TB0(I,IFR)=TT0(I,NZmodel)                  ! CLEAR-SKY BACKGROUND
            DTcir(I,IFR)=TT(I,NZmodel)-TT0(I,NZmodel)  ! CLOUD-INDUCED RADIANCE
         ENDDO

         CALL SENSITIVITY (DTcir,ZZT,NT,YP,YZ,NZmodel,PRESSURE,NZ, &
              &            delTAU,delTAUc,TAUeff,SS,               &
              &            BETA, BETAc, DDm, Dm, Z, DZ,            &
              &            N,NF,IFR,ISWI,RE) ! COMPUTE SENSITIVITY

 2000 CONTINUE                               ! END OF FREQUENCY LOOP   

 3000 CONTINUE                               ! END OF IWC LOOP

      CALL HEADER(5)

!-----------------------<<<<<<<<<<<<< END >>>>>>>>>>>>------------------------C

      RETURN
      END

! $Log: CloudForwardModel.f90,v      








