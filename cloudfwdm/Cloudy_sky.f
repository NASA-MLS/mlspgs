
      SUBROUTINE CLOUDY_SKY(ISPI,CWC,TEMP,F,NU,U,DU,PH,RC,IPSD,Dm)

C=======================================================
C     >>>>>>>>CLOUDY-SKY RADIATION SCHEME<<<<<<<<<<
C
C     CALCULATE CLOUD SCATTERING AND EXTINCTION
C     LATEST UPDATE: J.JIANG, MAY 20, 2001
C=======================================================
      
      IMPLICIT NONE
      REAL PI
      PARAMETER (PI=3.1415926)
      REAL*8 RE
      PARAMETER (RE=6370.d3)
      REAL F                                ! FREQUENCY IN GHz
      REAL WL                               ! WAVELENGTH IN METERS

C---------------------------------------
C     CLOUD PARAMETERS
C---------------------------------------

      INTEGER ISPI                          ! CLOUD TYPE (1:ICE,2:WATER)

      INTEGER IPSD                          ! SIZE-DISTRIBUTION TYPE
                                            ! 0 = USER DEFINED
                                            ! 1 = MH
                                            ! 2 = GAMMA (dme=150)
                                            ! 3 = LIU-CURRY 
                                            ! 4 = KNOLLENGBERG (b1=0.05)
                                            ! N.B. IF ISPD=0,2,4 USER MAY
                                            !      DEFINE NEW DISTRIBUTION
                                            !      OR CHANGE PARAMETERS 
                                            !      (dme,b1) IN drp_size.f  
      
      INTEGER NR                            ! NO OF PARTICLE SIZE BINS
      INTEGER NAB                           ! MAXIMUM NO OF A, B TERMS
      INTEGER NU                            ! DIMENSION FOR U
      INTEGER NU0                           ! MAXIMUM DIMENSION FOR U    

      PARAMETER (NR=40,NAB=50,NU0=256)

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
      REAL PH1(NU0)                         ! SINGLE PARTICLE PHASE FUNCTION

      COMPLEX S1,S2                         ! SINGLE PARTICLE AMPLITUDE FUNCTION

      REAL P(NAB,NU0)                       ! LEGENDRE POLYNOMIALS l=1
      REAL DP(NAB,NU0)                      ! Delt LEGENDRE POLYNOMIALS l=1
      

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

      CALL PFSETUP(NAB,P,DP,U,NU,NU0)       
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



