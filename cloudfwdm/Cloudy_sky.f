
      SUBROUTINE CLOUDY_SKY(ISPI,CWC,TEMP,F,NU,U,DU,PH,RC,IPSD,Dm,
     >                      PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR)

C=======================================================
C     >>>>>>>>CLOUDY-SKY RADIATION SCHEME<<<<<<<<<<
C
C     CALCULATE CLOUD SCATTERING AND EXTINCTION
C     LATEST UPDATE: J.JIANG, MAY 20, 2001
C=======================================================
      
      IMPLICIT NONE
      REAL PI
      PARAMETER (PI=3.1415926)
      REAL F                                ! FREQUENCY IN GHz
      REAL WL                               ! WAVELENGTH IN METERS

C---------------------------------------
C     CLOUD PARAMETERS
C---------------------------------------

      INTEGER ISPI                          ! CLOUD TYPE (1:ICE,2:WATER)

      INTEGER IPSD                          ! SIZE-DISTRIBUTION FLAG
                                            ! IILL     (I=ICE, L=LIQUID)
                                            ! 1000:     I->MH, L->GAMMA 
                                            ! 1100:     I->LIU-CURRY
                                            ! 2000-3900:I->MODIFIED GAMMA
                                            !              WITH VARIOUS De,
                                            !              alpha
                                            ! 4000-5900:I->KNOLLENBERG WITH 
                                            !              VARIOUS b1
                                            ! 6000:     I->PSD FOR PSC
      
      INTEGER NR                            ! NO OF PARTICLE SIZE BINS
      INTEGER NAB                           ! MAXIMUM NO OF A, B TERMS
      INTEGER NU                            ! DIMENSION FOR U

      INTEGER NABR(NR)                      ! TRUNCATION FOR A, B AT R

      REAL U(NU)                            ! COSINE OF SCATTERING ANGLES
      REAL DU(NU)                           ! DELTA U
      REAL CWC                              ! ICE WATER CONTENT
      REAL R(NR)                            ! PARTICLE RADIUS
      REAL RN(NR)                           ! NUMBER OF PARTICLES IN EACH BIN
      REAL TEMP                             ! TEMPERATURE (K)
      REAL BC(3,NR)                         ! SINGLE PARTICLE ABS,SCAT,EXT COEFFS
      REAL RC(3)                            ! VOLUME ABS,SCAT,EXT COEFFS
      REAL Dm                               ! MASS-MEAN-DIAMETER

      COMPLEX A(NR,NAB),B(NR,NAB)           ! MIE COEFFICIENCIES

      REAL PH(NU)                           ! INTERGRATED PHASE FUNCTION
      REAL PH1(NU)                         ! SINGLE PARTICLE PHASE FUNCTION

      COMPLEX S1,S2                         ! SINGLE PARTICLE AMPLITUDE FUNCTION

      REAL P(NAB,NU)                       ! LEGENDRE POLYNOMIALS l=1
      REAL DP(NAB,NU)                      ! Delt LEGENDRE POLYNOMIALS l=1
      

C---------------------------------------
C     WORK SPACES
C---------------------------------------
      INTEGER I,J,K
      REAL W1,SUM,DD,US

C-----------------------------------------------------------------------------

      WL=0.3/F
      
      IF (ISPI .EQ. 1) THEN
         DD=2000./NR**2                     ! DEFINE SIZE BINS FOR ICE CLOUD
      ELSE IF (ISPI .EQ. 2) THEN
         DD=200./NR**2                      ! DEFINE SIZE BINS FOR WATER CLOUD
      ENDIF

      DO I=1,NR
         R(I)=1+I**2*DD
         RN(I)=0.
      ENDDO

      CALL PFSETUP(NAB,P,DP,U,NU)       
      CALL DRP_SIZE(ISPI,R,RN,NR,CWC,TEMP,IPSD,Dm)
      CALL MIECOEFF(ISPI,F,TEMP,NR,R,A,B,NAB,NABR,BC)

      DO I=1,3
         RC(I)=0.
      ENDDO

      DO I=1,NR
         RC(2)=RC(2)+BC(2,I)*RN(I)*PI*R(I)*R(I)*1.E-12   ! SCATTERING
         RC(3)=RC(3)+BC(3,I)*RN(I)*PI*R(I)*R(I)*1.E-12   ! EXTINCTION
      ENDDO
      RC(1)=RC(3)-RC(2)                                  ! ABSORPTION

      DO K=1,NU
         PH(K)=0.
      ENDDO

      DO 200 I=1,NR

         W1=PI*R(I)*R(I)*1.E-12*RN(I)*BC(2,I)
         
         IF(W1 .GT. 0.) THEN
            SUM = 0.

            DO 100 K=1,NU

               US=SQRT(1-U(K)*U(K))
               S1=CMPLX(0.)
               S2=CMPLX(0.)

               DO J=1,NABR(I)
                  S1=S1+(2*J+1.)/J/(J+1.)*(A(I,J)*DP(J,K)
     >               +B(I,J)*P(J,K)/US)
                  S2=S2+(2*J+1.)/J/(J+1.)*(A(I,J)*P(J,K)/US
     >               +B(I,J)*DP(J,K))
               ENDDO
               PH1(K)=CABS(S1)**2+CABS(S2)**2
               SUM=SUM+PH1(K)*DU(K)

 100        ENDDO
   
            SUM=SUM*2           ! FACTOR 2 ACCOUNTS FOR HALF PLANE

            DO K=1,NU
               PH(K)=PH(K)+PH1(K)*W1/RC(2)/SUM
            ENDDO
         ENDIF
   
 200  CONTINUE   

      RETURN
      END

! $Log: Cloudy_sky.f,v      



