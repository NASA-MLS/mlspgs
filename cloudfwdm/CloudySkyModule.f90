! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module CloudySkyModule

! -------------------------------------------------------------------------  
! MLS CLOUD RADIANCE MODEL
! -------------------------------------------------------------------------

      use MieTheory, only: MieCoeff
      use MLSCommon, only: r8      
      use PhaseFunction, only: pfsetup
      use SizeDistribution, only: DRP_SIZE
      IMPLICIT NONE
      Private
      Public :: CLOUDY_SKY, CLOUD_MODEL

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE CLOUDY_SKY(ISPI,CWC,TEMP,F,NU,U,DU,PH,RC,IPSD,Dm, &
                 &          PH1,NAB,P,DP,NR,R,RN,BC,A,B,NABR)

!=======================================================
!     >>>>>>>>CLOUDY-SKY RADIATION SCHEME<<<<<<<<<<
!
!     CALCULATE CLOUD SCATTERING AND EXTINCTION
!     LATEST UPDATE: J.JIANG, MAY 20, 2001
!=======================================================

      REAL :: PI
      PARAMETER (PI=3.1415926)
      REAL(r8) :: F                            ! FREQUENCY IN GHz
      REAL(r8) :: WL                           ! WAVELENGTH IN METERS

!---------------------------------------
!     CLOUD PARAMETERS
!---------------------------------------

      INTEGER :: ISPI                          ! CLOUD TYPE (1:ICE,2:WATER)

      INTEGER :: IPSD                          ! SIZE-DISTRIBUTION FLAG
                                               ! IILL     (I=ICE, L=LIQUID)
                                               ! 1000:     I->MH, L->GAMMA 
                                               ! 1100:     I->LIU-CURRY
                                               ! 2000-3900:I->MODIFIED GAMMA
                                               !              WITH VARIOUS De,
                                               !              alpha
                                               ! 4000-5900:I->KNOLLENBERG WITH 
                                               !              VARIOUS b1
                                               ! 6000:     I->PSD FOR PSC
      
      INTEGER :: NR                            ! NO OF PARTICLE SIZE BINS
      INTEGER :: NAB                           ! MAXIMUM NO OF A, B TERMS
      INTEGER :: NU                            ! DIMENSION FOR U

      INTEGER :: NABR(NR)                      ! TRUNCATION FOR A, B AT R

      REAL(r8) :: U(NU)                        ! COSINE OF SCATTERING ANGLES
      REAL(r8) :: DU(NU)                       ! DELTA U
      REAL(r8) :: CWC                          ! ICE WATER CONTENT
      REAL(r8) :: R(NR)                        ! PARTICLE RADIUS
      REAL(r8) :: RN(NR)                       ! NUMBER OF PARTICLES IN EACH BIN
      REAL(r8) :: TEMP                         ! TEMPERATURE (K)
      REAL(r8) :: BC(3,NR)                     ! SINGLE PARTICLE ABS,SCAT,EXT COEFFS
      REAL(r8) :: RC(3)                        ! VOLUME ABS,SCAT,EXT COEFFS
      REAL(r8) :: Dm                           ! MASS-MEAN-DIAMETER

      COMPLEX(r8) :: A(NR,NAB),B(NR,NAB)       ! MIE COEFFICIENCIES

      REAL(r8) :: PH(NU)                       ! INTERGRATED PHASE FUNCTION
      REAL(r8) :: PH1(NU)                      ! SINGLE PARTICLE PHASE FUNCTION

      COMPLEX(r8) :: S1,S2                     ! SINGLE PARTICLE AMPLITUDE FUNCTION

      REAL(r8) :: P(NAB,NU)                    ! LEGENDRE POLYNOMIALS l=1
      REAL(r8) :: DP(NAB,NU)                   ! Delt LEGENDRE POLYNOMIALS l=1
      
!---------------------------------------
!     WORK SPACES
!---------------------------------------
      INTEGER :: I,J,K
      REAL(r8) :: W1,SUM,DD,US

!-----------------------------------------------------------------------------

      WL=0.3/F
      
      IF (ISPI .EQ. 1) THEN

        ! DEFINE SIZE BINS FOR ICE CLOUD
          DD=(2000._r8/NR**2)/(F/200._r8)  
!          DD=2000._r8/NR**2

      ELSE IF (ISPI .EQ. 2) THEN

        ! DEFINE SIZE BINS FOR WATER CLOUD
          DD=200._r8/NR**2                   

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
                  S1=S1+(2*J+1.)/J/(J+1.)*(A(I,J)*DP(J,K)      &
     &               +B(I,J)*P(J,K)/US)
                  S2=S2+(2*J+1.)/J/(J+1.)*(A(I,J)*P(J,K)/US    &
     &               +B(I,J)*DP(J,K))
               ENDDO
               PH1(K)=ABS(S1)**2+ABS(S2)**2
               SUM=SUM+PH1(K)*DU(K)

 100        ENDDO
   
            SUM=SUM*2           ! FACTOR 2 ACCOUNTS FOR HALF PLANE

            DO K=1,NU
               PH(K)=PH(K)+PH1(K)*W1/RC(2)/SUM
            ENDDO
         ENDIF
   
 200  CONTINUE   

      END SUBROUTINE CLOUDY_SKY


      SUBROUTINE CLOUD_MODEL(ITYPE,CHT,YZ,NH,WC)

!=========================================================================C
!  DEFINE VERTICAL PROFILES OF CLOUD ICE-WATER-CONTENT                    C
!  J.JIANG -05/18/2001                                                    C
!          -10/05/2001, MODIFIED TO FIT CLOUD RETREVIAL REQUIREMENTS      C
!=========================================================================C

      CHARACTER :: ITYPE

      INTEGER :: NH,I

      REAL(r8) :: YZ(NH)
      REAL(r8) :: WC(2,NH)

      REAL(r8) :: CLD_TOP
      REAL(r8) :: CLD_BASE
      REAL(r8) :: UPPER_LAG
      REAL(r8) :: LOWER_LAG
      REAL(r8) :: HT, CHT, WCscale

!--------------------------------------------------------------------

      IF (CHT .EQ. 0._r8) THEN
         HT = 16000. + CHT*1000.
      ELSE
         HT=CHT
      ENDIF

      WCscale=0.1

      IF(ITYPE .EQ. 'C') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: DEEP-CONVECTIVE SYSTEM '

         CLD_TOP   = HT + 1000.         ! CONVECTIVE CLOUD-TOP
         CLD_BASE  = HT - 11600.        ! CONVECTIVE CLOUD-BASE
         UPPER_LAG = HT - 8000.
         LOWER_LAG = HT - 10000.

         DO I=1,NH

            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/5000.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF

            ! LIQUID WATER CLOUD
            
            IF (YZ(I).GT.5000. .AND. YZ(I).LT. 10000.) THEN
               WC(2,I) = 0.1*0.
            ENDIF
         ENDDO


      ELSE IF (ITYPE .EQ. 'F') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: FRONTAL SYSTEM '

         CLD_TOP   = HT - 5000.         ! FRONTAL CLOUD-TOP
         CLD_BASE  = HT - 11000.        ! FRONTAL CLOUD-BASE
         UPPER_LAG = HT - 7000.
         LOWER_LAG = HT - 9000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/500.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 'A') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: ANVILS'

         CLD_TOP   = HT - 4000.         !ANVIL CLOUD-TOP
         CLD_BASE  = HT - 11000.        !ANVIL CLOUD-BASE
         UPPER_LAG = HT - 5000.
         LOWER_LAG = HT - 6000.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-UPPER_LAG)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(LOWER_LAG-YZ(I))/600.)
            ENDIF
            WC(2,I) = 0.0
         ENDDO

      ELSE IF (ITYPE .EQ. 'T') THEN

!         PRINT*,' '
!         PRINT*,' -> CLOUD-TYPE: THIN-LAYER CIRRUS'

         CLD_TOP   = HT + 500.         !1KM THICK CLOUD LAYER
         CLD_BASE  = HT - 500.
         UPPER_LAG = HT + 200.
         LOWER_LAG = HT - 200.

         DO I=1,NH
            IF (YZ(I) .GT. CLD_TOP .OR. YZ(I) .LT. CLD_BASE) THEN
               WC(1,I) = 0.0
            ELSE IF (YZ(I).LE.UPPER_LAG .AND. YZ(I).GE.LOWER_LAG) THEN
               WC(1,I) = 1.0
            ELSE IF (YZ(I).GT.UPPER_LAG .AND. YZ(I).LE.CLD_TOP) THEN
               WC(1,I) = 1.0*EXP(-(YZ(I)-HT)/500.)
            ELSE IF (YZ(I).LT.LOWER_LAG .AND. YZ(I).GE.CLD_BASE) THEN
               WC(1,I) = 1.0*EXP(-(HT-YZ(I))/500.)
            ENDIF
            WC(2,I)=0.0
         ENDDO
 
      ELSE

        PRINT*, 'NO CLOUDS'
       
        DO I=1,NH
            WC(1,I)=0.0
            WC(2,I)=0.0
        ENDDO

      ENDIF


        DO I=1,NH
            WC(1,I)=WCscale*WC(1,I)
            WC(2,I)=WCscale*WC(2,I)
        ENDDO


!------------------------------------------------------

      END SUBROUTINE CLOUD_MODEL


end module CloudySkyModule

! $Log$
! Revision 1.8  2002/04/15 22:22:40  jonathan
! check bug
!
! Revision 1.7  2002/04/15 20:08:04  jonathan
! add frequency-dependent size bin
!
! Revision 1.6  2001/10/08 21:40:21  jonathan
! add cloud_model sub
!
! Revision 1.5  2001/09/21 15:51:37  jonathan
! modified F95 version
!

