
      SUBROUTINE CloudForwardModel (NF, NZ, NT, NS,
     1           FREQUENCY,PRESSURE,HEIGHT,TEMPERATURE,VMRin,IWC,LWC,
     1           ZT, RE, IPSD, ISURF, ISWI, ICON,
     2           TB0,DTcir,BETA,BETAc,Dm,TAUeff,SS )  

C===========================================================================C
C   >>>>>>>> CLOUDY-SKY FORWARD MODEL FOR MICROWAVE LIMB SOUNDER <<<<<<<<   C
C---------------------------------------------------------------------------C
C                                                                           C
C     DESCRIPTION OF MODEL PARAMETERS           NOTES                       C
C     -------------------------------           -----                       C
C     1. INPUT:  FREQUENCY  (NF)                                            C 
C                PRESSURE   (NZ)    . The internal model resuires NZ=640    C
C                HEIGHT     (NZ)      and the top model layer is defined    C
C                TEMPERATURE(NZ)      at PRESSURE(640)=10^(-(40./16.-3.)).  C 
C                VMRin      (NS,NZ) . Volume mixing ratios (N:1=H2O,2=O3)   C
C                IWC        (NZ)                                            C
C                LWC        (NZ)                                            C
C                ZT         (NT)    . Tangent pressure                      C
C                RE                 . Earth radius                          C
C                IPSD,ISURE,ISWI,   . User determined size-distribution,    C
C                ICON                 surface-type, sensitivity calcula-    C
C                                     tion and cloud control switches.      C 
C                                     Default: (0,0,0,2)                    C
C                                                                           C
C     2. OUTPUT: TB0     (NT,NF)    . Background clear-sky radiance         C 
C                DRcir   (NT,NF)    . Cloud induced radiance                C
C                BETA    (NZ-1,NF)  . Total extinction profile              C
C                BETAc   (NZ-1,NF)  . Cloud extinction profile              C
C                TAUeff  (NT,NF)    . Effective cloud optical depth         C
C                SS      (NT,NF)    . Cloud radiance sensitivity            C
C                Dm      (2,NZ-1)   . Mass-mean-diameters (1=ice,2=water)   C
C     -------------------------------                                       C
C     FREQUENCY RANGE: 1-3000GHz                                            C
C                                                                           C
C     VERSION: 1.0, MAY 18, 2001                                            C
C     MICROWAVE ATMOSPHERIC SCIENCE TEAM                                    C
C     JET PROPULSION LABORATORY                                             C
C     CALIFORNIA INSTITUTE OF TECHNOLOGY                                    C
C     4800 OAK GROVE DRIVE                                                  C
C     PASADENA, CA 91109-8099                                               C
C                                                                           C
C     EMAIL: JONATHAN@MLS.JPL.NASA.GOV                                      C
C     PHONE: (818) 354-7135                                                 C
C     FAX:   (818) 393-5065                                                 C
C===========================================================================C

      IMPLICIT NONE

C---------------------------------------
C     INPUT PARAMETERS (INPUTS FROM L2)        ! -- INTERFACE AEA -- ! 
C---------------------------------------

      INTEGER NF                               ! NUMBER OF FREQUENCIES
      INTEGER NZ                               ! NUMBER OF PRESSURE LEVELS
      INTEGER NT                               ! NUMBER OF TANGENT HEIGHTS
      INTEGER NS                               ! NUMBER OF SPECIES

      INTEGER IPSD                             ! SIZE-DISTRIBUTION TYPE
                                               ! 0 = MH          
                                               ! 1 = LIU-CURRY
                                               ! 2 = PSD FOR LIQUID WATER
                                               ! 3 = PSD FOR PSC
                                               ! 10-19 MODIFIED GAMMA WITH 
                                               !       VARIOUS Deff, Alpha
                                               ! 20-29 KNOLLENBERG WITH 
                                               !       VARIOUS b1

      INTEGER ISURF                            ! SURFACE TYPE
                                               ! 0 = SIMPLE MODEL
                                               ! 1 = LAND
                                               ! 2 = SEA 
  
      INTEGER ISWI                             ! SENSITIVITY SWITCH
                                               ! 0 = NO SENSITIVITY CALCULATION
                                               ! 1 = COMPUTE SENSITIVITY

      INTEGER ICON                             ! CLOUD CONTROL SWITCH
                                               ! 0 = NO CLOUD
                                               ! 1 = NEAR SIDE CLOUD ONLY
                                               ! 2 = DEFAULT

      REAL FREQUENCY(NF)                       ! FREQUENCIES (GHz)
      REAL PRESSURE(NZ)                        ! PRESSURE LEVEL
      REAL HEIGHT(NZ)                          ! PRESSURE HEIGHT
      REAL TEMPERATURE(NZ)                     ! ATMOSPHERIC TEMPERATURE
      REAL VMRin(NS,NZ)                        ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3 VOLUME MIXING RATIO

      REAL IWC(NZ)                             ! ICE WATER CONTENT
      REAL LWC(NZ)                             ! LIQUID WATER CONTENT

      REAL ZT(NT)                              ! TANGENT PRESSURE
      REAL*8 RE                                ! EARTH RADIUS

C--------------------------------------
C     OUTPUT PARAMETERS (OUTPUT TO L2)        
C--------------------------------------

      REAL TB0(NT,NF)                          ! CLEAR-SKY TB AT ZT
      REAL DTcir(NT,NF)                        ! CLOUD-INDUCED RADIANCE

      REAL TAUeff(NT,NF)                       ! CLOUD EFFECTIVE OPTICAL DEPTH
      REAL SS(NT,NF)                           ! CLOUD RADIANCE SENSITIVITY
                                               ! (NT+1) FOR ZENITH LOOKING) 

      REAL BETA(NZ-1,NF)                       ! TOTAL OPTICAL DEPTH
      REAL BETAc(NZ-1,NF)                      ! CLOUDY OPTICAL DEPTH

      REAL Dm(2,NZ-1)                          ! MASS-MEAN-DIAMETER

                                               ! -- END OF INTERFACE AREA -- !

C-----------------------------------------------------------------------------
                                               
C-------------------------------
C     INTERNAL MODEL PARAMETERS                ! -- INTERNAL AREA -- !
C-------------------------------

      REAL PI
      PARAMETER (PI=3.1415926)
      
      INTEGER L                                ! NO. OF ATMOS. LAYERS
      INTEGER NFREQ                            ! NO. OF FREQ READ IN
      INTEGER NU                               ! NO. OF SCATTERING ANGLES
      INTEGER NUA                              ! NO. OF SCAT. AZIMUTH ANGLES
      INTEGER NIWC                             ! NO. OF ICE WATER CONTENTS
      INTEGER N                                ! NO. OF AEROSOL SPECIES

      INTEGER MAXFS                            ! MAX. NO. OF FREQUENCIES
      INTEGER MNT                              ! MAX. NO. OF TANGENT HEIGHTS 

      PARAMETER (L=639,NU=16,NUA=8,N=2,NIWC=10,
     >           MAXFS=500,MNT=500)            ! THESE ARE THE FIXED INTERNAL 
                                               ! PARAMETERS. WE SET THESE 
                                               ! CONSTANTS HERE TO REDUCE THE 
                                               ! NUMBER OF PARAMETERS PASSED 
                                               ! THROUGH L2

      INTEGER LORS                             ! SURFACE TYPE: LAND OR SEA

      REAL U(NU)                               ! COSINE OF SCATTERING ANGLES
      REAL DU(NU)                              ! DELTA U
      REAL UA(NUA)                             ! COSINE OF SCAT AZIMUTH ANGLES
      REAL THETA(NU)                           ! SCATTERING ANGLES
      REAL PHI(NUA)                            ! SCATTERING AZIMUTH ANGLES
      REAL UI(NU,NU,NUA)                       ! COSINE OF INCIDENT TB ANGLES
      REAL THETAI(NU,NU,NUA)                   ! ANGLES FOR INCIDENT TB
      REAL PHH(N,NU,L)                         ! PHASE FUNCTION 
      REAL TAU(L)                              ! TOTAL OPTICAL DEPTH
      REAL TAU100(L)                           ! TOTAL OPTICAL DEPTH AT 100%RH
      REAL W0(N,L)                             ! SINGLE SCATTERING ALBEDO
      
      REAL S                                   ! SALINITY
      REAL SWIND                               ! SEA SURFACE WIND
      REAL TS                                  ! SURFACE TEMPERATURE (K)
      REAL RS(NU/2)                            ! SURFACE REFLECTIVITY
      REAL RC0(3)                              ! GAS ABS.SCAT.EXT COEFFS.
      REAL RC(N,3)                             ! CLOUD ABS.SCAT.EXT COEFFS.
      REAL RC_TOT(3)                           ! TOTAL ABS.SCAT.EXT COEFFS.
      REAL Z(L)                                ! MODEL LAYER THICKNESS (m)
      REAL TAU0(L)                             ! CLEAR-SKY OPTICAL DEPTH
      REAL TEMP(L)                             ! MEAN LAYER TEMPERATURE (K)

      REAL TT(MNT,L+1)                         ! CLOUDY-SKY TB AT TANGENT 
                                               ! HEIGHT ZT (LAST INDEX FOR 
                                               ! ZENITH LOOKING)
      REAL TT0(MNT,L+1)                        ! CLEAR-SKY TB AT TANGENT
                                               ! HEIGHT ZT

      REAL FREQ(MAXFS)                         ! FREQUENCIES (GHz)
  
C------------------------------------
C     ATMOSPHERIC PROFILE PARAMETERS
C------------------------------------

      REAL YZ(L+1)                             ! PRESSURE HEIGHT (m)
      REAL YP(L+1)                             ! PRESSURE (hPa)
      REAL YT(L+1)                             ! TEMPERATURE PROFILE
      REAL YQ(L+1)                             ! H2O VOLUME MIXING RATIO
      REAL VMR(5,L+1)                          ! 1=O3 VOLUME MIXING RATIO
                               
C----------------------------
C     CLOUD MODEL PARAMETERS
C----------------------------

      REAL WC(N,L+1)                           ! CLOUD WATER CONTENT (g/m3)
                                               ! (N=1 ICE, N=2 WATER)
      REAL CWC                                 ! CLOUD WATER CONTENT (g/m3)

      REAL CDEPTH(N)                           ! CLOUD OPTICAL DEPTH
      REAL DEPTH                               ! TOTAL OPTICAL DEPTH
                                               ! (CLEAR+CLOUD)
      REAL DDm(N,L)                            ! MASS-MEAN-DIAMETER

      REAL delTAU(L)                           ! TOTAL EXTINCTION
      REAL delTAUc(L)                          ! CLOUDY-SKY EXTINCTION
      
C---------------------------
C     WORK SPACE PARAMETERS
C---------------------------

      INTEGER I, J, K, IFR, ILYR, IL,ISPI,IIWC, ICLD_TOP,MY_NIWC

      REAL HT,DMA,RATIO
      REAL PH0(N,NU,L),W00(N,L)  
      REAL P11(NU), RC11(3),RC_TMP(N,3)
      REAL CHK_CLD(L+1)                        
      REAL ZZT(MNT)


C---------------<<<<<<<<<<<<< START EXCUTION >>>>>>>>>>>>-------------------C

      CALL HEADER(1)

C=========================================================================
C                    >>>>>> CHECK MODEL-INPUT <<<<<<< 
C-------------------------------------------------------------------------
C     CHECK IF THE INPUT PROFILE MATCHS THE MODEL INTERNAL GRID; 
C     SET TANGENT PRESSURE (hPa) TO TANGENT HEIGHT (km)
C=========================================================================

      CALL MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMRin,NZ,NS,  
     >                  YP,YZ,YT,YQ,VMR,L+1,ZT,ZZT,NT) 
                                                       
      CALL CLOUD_CHECK(PRESSURE,IWC,LWC,NZ,WC,CHK_CLD,L+1)     

C-----------------------------------------------
C     INITIALIZE SCATTERING AND INCIDENT ANGLES 
C-----------------------------------------------

      CALL ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

C----------------------------------
C     SURFACE INFORMATION
C----------------------------------

      LORS  = ISURF          
      TS    = 288.
      S     = 35.
      SWIND = 0.

      IF (ISWI .EQ. 0) THEN                  
         RATIO=1.
         MY_NIWC=1                ! SKIP FULL SENSITIVITY CALCULATION
      ELSE
          MY_NIWC=NIWC
      endif

C------------------------------------------
C     PERFORM FULL SENSITIVITY CALCULATION
C------------------------------------------

      DO 3000 IIWC=1, MY_NIWC
         if (ISWI .EQ. 0) THEN
            ratio=1.
         else
      RATIO = 10.*(IIWC-1)**2*0.004+1.0E-9
      endif
C=========================================================================
C                   >>>>>>> CLEAR-SKY MODULE <<<<<<<<
C-------------------------------------------------------------------------
C     COMPUTE CLEAR-SKY ABSORPTION COEFFICIENTS, INCLUDING DRY, WET 
C     CONTINUMA AND LINE EMISSIONS.   
C=========================================================================

      NFREQ=NF

      DO 2000 IFR=1, NFREQ      ! START OF FREQUENCY LOOP

         FREQ(IFR)=FREQUENCY(IFR)

         CALL CLEAR_SKY(L,NU,TS,S,LORS,SWIND,YZ,YP,YT,YQ,VMR,
     >                     FREQ(IFR),RS,U,TEMP,TAU0,Z,TAU100) 

         CALL HEADER(3)

C-----------------------------------------------------
C        ASSUME 100% SATURATION IN CLOUD LAYER
C-----------------------------------------------------

         DO IL=1,L                        ! 100% SATURATION INSIDE CLOUD 
            IF(CHK_CLD(IL) .NE. 0.)THEN
               ICLD_TOP=IL
               IF(YZ(IL) .LT. 20.)THEN
                  TAU0(IL)=TAU100(IL)
               ENDIF
            ENDIF
         ENDDO

C--------------------------------------------------------

         DO 1000 ILYR=1,L                 ! START OF MODEL LAYER LOOP:   
 
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
               W0(ISPI,ILYR)=0.
               W00(ISPI,ILYR)=0.
            ENDDO

            DEPTH  = 0.
            CWC = 1.E-9

            IF(CHK_CLD(ILYR) .NE. 0.) THEN 

               DO ISPI=1,N
                  CWC = WC(ISPI,ILYR) 
                  CWC = MAX(1.E-9,CWC)
              
C=================================================
C    >>>>>>>>> CLOUDY-SKY MODULE <<<<<<<<<<<
C=================================================
       
                  CALL CLOUDY_SKY(ISPI,CWC,TEMP(ILYR),FREQ(IFR),
     >                            NU,U,DU,P11,RC11,IPSD,DMA)

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

            delTAU(ILYR) = DEPTH
            delTAUc(ILYR)= CDEPTH(1)

 1000    CONTINUE                         ! END OF MODEL LAYER LOOP

C==================================================
C    >>>>>>> RADIATIVE TRANSFER MODULE <<<<<<<<<<
C==================================================

         CALL HEADER(4)

         CALL RADXFER(L,NU,NUA,U,DU,PH0,NT,ZZT,W00,TAU0,RS,TS,FREQ(IFR),
     >                YZ,TEMP,N,THETA,THETAI,PHI,UI,UA,TT0,MNT,ICON,RE)  !CLEAR-SKY

         IF(ICON .NE. 0) THEN                                          

           CALL RADXFER(L,NU,NUA,U,DU,PHH,NT,ZZT,W0,TAU,RS,TS,FREQ(IFR),
     >                YZ,TEMP,N,THETA,THETAI,PHI,UI,UA,TT,MNT,ICON,RE)   !CLOUDY-SKY

         ENDIF

C====================================
C    >>>>>>> MODEL-OUTPUT <<<<<<<<<
C====================================

         DO I=1,NT
            TB0(I,IFR)=TT0(I,L+1)                  ! CLEAR-SKY BACKGROUND      
            DTcir(I,IFR)=TT(I,L+1)-TT0(I,L+1)      ! CLOUD-INDUCED RADIANCE
         ENDDO

         CALL SENSITIVITY (DTcir,ZZT,NT,YP,YZ,L+1,PRESSURE,NZ,
     >                     delTAU,delTAUc,BETA,BETAc,TAUeff,SS,
     >                     DDm,Dm,N,NF,IFR,ISWI,RE) ! COMPUTE SENSITIVITY

 2000 CONTINUE                               ! END OF FREQUENCY LOOP   

c      IF (ISWI .EQ. 0) THEN
c         GOTO 200
c      ENDIF

 3000 CONTINUE

      CALL HEADER(5)


C-----------------------<<<<<<<<<<<<< END >>>>>>>>>>>>------------------------C

      RETURN
      END

! $Log: CloudForwardModel.f,v      








