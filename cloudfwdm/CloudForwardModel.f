
      SUBROUTINE CloudForwardModel (FREQUENCY,PRESSURE,HEIGHT,
     1            TEMPERATURE,VMRin,ZT, NF,NZ,NT,
     2            WCscale,ITYPE,IPSD,ISURF,ISWI,IRHB,CHT,
     3            TB0,DTcir,DTAU,DTAUc,DZ,Dm,TAUeff,SS)  

C===========================================================================C
C   >>>>>>>>>> CLOUDY-SKY MICROWAVE RADIATIVE TRANSFER PROGRAM <<<<<<<<<<   C
C---------------------------------------------------------------------------C
C                                                                           C
C     DESCRIPTION OF MODEL PARAMETERS                                       C
C     -------------------------------                                       C
C     1. INPUT: FREQUENCY, L2 ATMOSPHERIC PROFILES (INDICED AT PRESSURE     C 
C                LEVELS), AND TANGENT HEIGHTS.                              C
C     2. USER DEFINED CLOUD/MODEL PARAMETERS. DEFAULT: 0.,0,0,0,0,0,0.      C
C     3. OUTPUT: TB0, DTcir:  CLEAR-SKY, AND CLOUD-INDUCED RADIANCES FOR    C 
C                             EACH TANGENT HEIGHT ZT;                       C
C                DTAU, DTAUc: TOTAL, AND CLOUD OPTICAL DEPTHS FOR EACH      C
C                             ATMOSPHERIC LAYER DZ;                         C
C                TAUeff, SS:  EFFECTIVE CLOUD OPTICAL DEPTH, AND CLOUD      C
C                             RADIANCE SENSITIVITY AT TANGENT HEIGHT ZT;    C
C                Dm:          MASS-MEAN-DIAMETER OF CLOUD HYDROMETERS.      C
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

      INTEGER ITYPE                            ! CLOUD MODEL TYPE
                                               ! 0 = CLEAR-SKY    
                                               ! 1 = DEEP-CONVECTIVE
                                               ! 2 = FRONTAL
                                               ! 3 = ANVIL
                                               ! 4 = THIN-LAYER
                                               ! 101-104 NEAR SIDE CLOUDS ONLY 

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
     
      INTEGER IRHB                             ! RELITIVE HUMIDITY BELOW CLOUD
                                               ! 0 = < 100% R.H. BELOW CLOUD 
                                               ! 1 =   100% R.H. BELOW CLOUD

      REAL CHT                                 ! CLOUD PROFILE/HEIGHT SHIFT
                                               ! 0.0 = DEFAULT CLOUD PROFILE
                                               ! 1.0 = MOVE CLOUD 1KM UP
                                               ! -1. = MOVE CLOUD 1KM DOWN

      REAL WCscale                             ! CLOUD WATER SCALING FACTOR
                                               ! (DEFAULT MAX IWC/LWC = 1.)

      REAL FREQUENCY(NF)                       ! FREQUENCIES (GHz)
      REAL PRESSURE(NZ)                        ! PRESSURE LEVEL
      REAL HEIGHT(NZ)                          ! PRESSURE HEIGHT
      REAL TEMPERATURE(NZ)                     ! ATMOSPHERIC TEMPERATURE

      REAL VMRin(5,NZ)                         ! 1=H2O VOLUME MIXING RATIO
                                               ! 2=O3 VOLUME MIXING RATIO

      REAL ZT(NT)                              ! TANGENT HEIGHT

C--------------------------------------
C     OUTPUT PARAMETERS (OUTPUT TO L2)        
C--------------------------------------

      REAL TB0(NT+1,NF)                        ! CLEAR-SKY TB AT ZT
      REAL DTcir(NT+1,NF)                      ! CLOUD-INDUCED RADIANCE

      REAL TAUeff(NT+1,NF)                     ! CLOUD EFFECTIVE OPTICAL DEPTH
      REAL SS(NT+1,NF)                         ! CLOUD RADIANCE SENSITIVITY
                                               ! ZT(NT+1) FOR ZENITH LOOKING 

      REAL DTAU(NZ-1,NF)                       ! TOTAL OPTICAL DEPTH
      REAL DTAUc(NZ-1,NF)                      ! CLOUDY OPTICAL DEPTH
      REAL DZ(NZ-1)                            ! LAYER THICKNESS

      REAL Dm(2,NZ-1)                          ! MASS-MEAN-DIAMETER

                                               ! -- END OF INTERFACE AREA -- !
C-----------------------------------------------------------------------------
                                               
C-------------------------------
C     INTERNAL MODEL PARAMETERS                ! -- INTERNAL AREA -- !
C-------------------------------

      REAL PI
      REAL*8 RE                                ! EARTH RADIUS
      PARAMETER (PI=3.1415926, RE=6370.D3)
      
      INTEGER L                                ! NO. OF ATMOS. LAYERS
      INTEGER NFREQ                            ! NO. OF FREQ READ IN
      INTEGER NU                               ! NO. OF SCATTERING ANGLES
      INTEGER NUA                              ! NO. OF SCAT. AZIMUTH ANGLES
      INTEGER NIWC                             ! NO. OF ICE WATER CONTENTS
      INTEGER N                                ! NO. OF AEROSOL SPECIES

      INTEGER MAXFS                            ! MAX. NO. OF FREQUENCIES
      INTEGER MNT                              ! MAX. NO. OF TANGENT HEIGHTS 

      PARAMETER (L=639,NU=16,NUA=8,MAXFS=500,MNT=500,N=2)

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

      REAL TT(MNT+1,L+1)                       ! CLOUDY-SKY TB AT TANGENT 
                                               ! HEIGHT ZT (LAST INDEX FOR 
                                               ! ZENITH LOOKING)
      REAL TT0(MNT+1,L+1)                      ! CLEAR-SKY TB AT TANGENT
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

      INTEGER I, J, K, IFR, ILYR, IL,ISPI,IIWC, ICLD_TOP

      REAL HT,DMA
      REAL PH0(N,NU,L),W00(N,L)  
      REAL P11(NU), RC11(3),RC_TMP(N,3)
      REAL CHK_CLD(L+1)                        


C---------------<<<<<<<<<<<<< START EXCUTION >>>>>>>>>>>>-------------------C

      CALL HEADER(1)

C==============================
C   >>>>>> MODEL-INPUT <<<<<<< 
C==============================

      CALL MODEL_ATMOS(PRESSURE,HEIGHT,TEMPERATURE,VMRin,NZ,  
     >                  YP,YZ,YT,YQ,VMR,L+1)         ! INTERPOLATE NZ TO L+1
      CALL CLOUD_MODEL(ITYPE,CHT,WC,CHK_CLD,YZ,L+1)  ! DEFINE CLOUD PROFILES

C-----------------------------------------------
C     INITIALIZE SCATTERING AND INCIDENT ANGLES 
C-----------------------------------------------

      CALL ANGLE(THETA,U,DU,NU,PHI,UA,NUA,UI,THETAI)

C----------------------------------
C        SURFACE INFORMATION
C----------------------------------

         LORS  = ISURF          
         TS    = 288.
         S     = 35.
         SWIND = 0.

C=======================================
C    >>>>>>> CLEAR-SKY MODULE <<<<<<<<
C=======================================

      NFREQ=NF

      DO 9000 IFR=1, NFREQ                  ! START OF FREQUENCY LOOP

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

         IF (IRHB .EQ. 1) THEN
            DO IL=1,ICLD_TOP              ! ASSUME 100% SATURATION BELOW CLOUD 
               TAU0(IL)=TAU100(IL)        ! (OPTIONAL FOR CONVECTIVE CLOUD)
            ENDDO
         ENDIF   

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
                  CWC = WCscale*WC(ISPI,ILYR) 
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

         CALL RADXFER(L,NU,NUA,U,DU,PH0,NT,ZT,W00,TAU0,RS,TS,FREQ(IFR),
     >                YZ,TEMP,N,THETA,THETAI,PHI,UI,UA,TT0,MNT,0)   ! CLEAR-SKY

         IF(ITYPE .NE. 0) THEN                                      ! CLOUDY-SKY
            CALL RADXFER(L,NU,NUA,U,DU,PHH,NT,ZT,W0,TAU,RS,TS,FREQ(IFR),
     >                   YZ,TEMP,N,THETA,THETAI,PHI,UI,UA,TT,MNT,ITYPE) 

         ENDIF

C====================================
C    >>>>>>> MODEL-OUTPUT <<<<<<<<<
C====================================

         DO I=1,NT+1
            TB0(I,IFR)=TT0(I,L+1)                  ! CLEAR-SKY BACKGROUND           
            DTcir(I,IFR)=TT(I,L+1)-TT0(I,L+1)      ! CLOUD-INDUCED RADIANCE
         ENDDO

         CALL SENSITIVITY (DTcir,ZT,NT,YP,YZ,L+1,PRESSURE,NZ,
     >                     delTAU,delTAUc,DTAU,DTAUc,TAUeff,SS,
     >                     DDm,Dm,Z,DZ,N,NF,IFR,ISWI)  ! SENSITIVITY PARAMETERS

 9000 CONTINUE                               ! END OF FREQUENCY LOOP   

      CALL HEADER(5)


C-----------------------<<<<<<<<<<<<< END >>>>>>>>>>>>------------------------C

      RETURN
      END

     








