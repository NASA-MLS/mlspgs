
      SUBROUTINE GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2,
     >                    MOL,NMOL,NCNT,T,PB,F,RH,VMR,ABSC)

C==============================================================
C      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
C      LATEST UPDATE: J.JIANG, MAY 18, 2001
C==============================================================

      IMPLICIT NONE
      INCLUDE 'spectra.h'
      REAL QTP(3)                      ! TEMPERATURE ON WHICH QLG ARE GIVEN
      DATA QTP /300.0,225.0,150.0/

      REAL QRAT                        ! INTERPOLATED QLG RATIO AT A GIVEN TEMP
      REAL F                           ! FREQUENCY IN GHz
      REAL T                           ! TEMPERATURE (K)
      REAL P                           ! DRY AIR PARTIAL PRESSURE (hPa)
      REAL PB                          ! TOTAL AIR PRESSURE (hPa)
      REAL VP                          ! VAPOR PARTIAL PRESSURE (hPa)
      REAL VMR(5)                      ! MINOR SPECIES 1-O3
      REAL VMR_H2O                     ! H2O VOLUME MIXING RATIO
      REAL VMR_O2                      ! O2 VOLUME MIXING RATIO
      REAL B                           ! BETA (1/m/ppv)
      REAL RH                          ! Relative Humidity
      REAL ABSC                        ! ABSORPTION COEFFICIENT (1/m)
                                       ! 1-O2
                                       ! 2-H2O
                                       ! 3-O_18_O
                                       ! 4-H2O_18

      REAL CONT_1,CONT_2,CONT_3        ! CONTINUUM ABSORPTION COEFFICIENTS
      REAL SC_CONST
      REAL SD,G0                       ! DEBY CONTRIBUTION (LIEBE 1989)

      REAL PI
      PARAMETER (PI=3.1415926)
      REAL ZP,YY,TT,TWTH0,DWTH0        ! WORKING SPACE
      INTEGER I,J,IMOL                

      REAL*8 MYSHAPE,FF
      EXTERNAL MYSHAPE

C------------------------------------------------------------------------

      IF(RH .NE. 100.) THEN
         VMR_H2O = RH                  ! PH IS WATER VAPOR MIXING RATIO
         VP=VMR_H2O*PB                 ! VP IS VAPOR PRESSURE, PB IS TOTAL
         P=PB-VP                       ! PRESSURE, P IS DRY-AIR PRESSURE
      ELSE IF (RH .EQ. 100.) THEN
         CALL RHtoEV(PB,T,100.,VP)     ! RH HERE IS 100% RELATIVE HUMIDITY 
         P = PB-VP
         VMR_H2O = VP/P
      END IF

      VMR_O2 = 0.209476
      FF     = F*1000.                 ! CONVERT F TO MHz
      ZP     = -ALOG10(P)              ! CONVERT P(hPa) TO ZP
      TT     = 300./T

      ABSC=0.
      I=0

C-------------------------------------
C     LINE EMISSION
C-------------------------------------
      DO IMOL=1,NMOL
         B=0.
         DO J=1,NCNT(IMOL)
            I=I+1

            IF(IMOL.GE.3.AND.ABS(V0(i)-FF).GT.10000.) GOTO 100 ! SAVE CPU

            IF(T .LE. QTP(2)) THEN
c               QRAT = (QLG(3,I)-QLG(1,I))*(T-QTP(2))/(QTP(3)-QTP(2))
               QRAT = (QLG(2,I)-QLG(1,I))+(QLG(3,I)-QLG(2,I))*(T-QTP(2))/(QTP(3)-QTP(2))
            ELSE
               QRAT = (QLG(2,I)-QLG(1,I))*(T-QTP(1))/(QTP(2)-QTP(1))
            ENDIF

            YY=10.**(-ZP)*(DELTA(I)*TT**N1(I)+GAMMA(I)*TT**N2(I))

C-------------------------------------
C           PRESSURE WIDTH
C-------------------------------------

            TWTH0=WTH(I)*10.**(-ZP)*TT**NTH(I)

C-------------------------------------
C           DOPPLER WIDTH
C-------------------------------------

            IF(IMOL.EQ.1) DWTH0 = 3.58e-7*SQRT(T/32.)*FF
            IF(IMOL.EQ.2) DWTH0 = 3.58e-7*SQRT(T/18.)*FF
            IF(IMOL.EQ.3) DWTH0 = 3.58e-7*SQRT(T/34.)*FF
            IF(IMOL.EQ.4) DWTH0 = 3.58e-7*SQRT(T/20.)*FF
            IF(IMOL.EQ.5) DWTH0 = 3.58e-7*SQRT(T/48.)*FF
            
C--------------------------------------------------
C           WHITING'S APPROXIMATION FOR THE VOIGT
C--------------------------------------------------

            TWTH0 = TWTH0/2+SQRT(TWTH0**2/4+DWTH0**2)

C---------------------------------------------------------------
C	TO SAVE COMPUTING TIME CHECK IF FREQ IS CLOSE TO O3 LINES
C----------------------------------------------------------------

            B = B+2.30549d09*10.0**(IST(I) - QRAT + GSE(I)*
     >		(1.d0/300-1.d0/T)/1.600386 - ZP) * PI * 
     >		 MYSHAPE(V0(i),FF,YY,TWTH0) * FF/V0(I) *
     >		(1.d0 - dEXP(-V0(i)/(20836.7d0*T)))/
     >		(1.d0 - dEXP(-V0(i)/(20836.7d0*300.0)))/T 

 100     CONTINUE

         ENDDO   

         B=AMAX1(0.,B)                 ! AVOID NEGATIVE LINE SHAPE     

         IF(IMOL .EQ. 1) ABSC = ABSC + B*VMR_O2                   ! O2      
         IF(IMOL .EQ. 2) ABSC = ABSC + B*VMR_H2O                  ! H2O
         IF(IMOL .EQ. 3) ABSC = ABSC + B*VMR_O2*0.00409524        ! O_18_O
         IF(IMOL .EQ. 4) ABSC = ABSC + B*VMR_H2O*0.00204          ! H2O_18
         IF(IMOL .EQ. 5) ABSC = ABSC + B*VMR(1)*1e-6              ! O3

      ENDDO  

C----------------------------------------
C     DRY CONTINUUM (ATBD APPENDIX F)
C----------------------------------------

      B=(7.7e-10*EXP(-1.5e-3*(F/30)**2)+1.e-13*(60**2+(F/30)**2)
     >     *EXP(-1.e-4*(F/30)**2))*TT**1.7

      ABSC = ABSC + B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5

C     CONT_1=1.4e-10*(1-1.2e-5*F**1.5)     ! LIEBE 1989
C     CONT_1 = 1.4e-12/(1+1.93e-5*F**1.5)  ! LIEBE 1993
C
C      CONT_1 = 1.4e-12/(1+1.93e-5*F**1.5)*1.28*1.17  ! MATCH THE UARS FIT
C
C      B=CONT_1 * F * (P/10.)**2*TT**3.5
C      ABSC=ABSC + B*.182*F/4.343 
C------------------------------------------
C     DEBYE TERM FROM LIEBE
C------------------------------------------
      SD = 6.14e-4*P/10.*TT**2
      G0 = 5.6e-3*(1+1.1*VMR_H2O)*P/10.*TT
      B = SD*F/(G0 * (1.+ (F/G0)**2))		! B is Npp in LIEBE
      ABSC = ABSC + B*.182*F/4.343              ! CONVERTED TO 1/km 

C----------------------------------------
C     WET CONTINUUM 
C----------------------------------------
C     CONT_1 = 1.28e-15 	! BILL'S VALIDATION PAPER
      CONT_1 = 7.53e-16 
      CONT_2 = 4.20
      CONT_3 = 0.00
      SC_CONST = CONT_1 * FF**2 * EXP(-CONT_3 * FF**2)
      B = SC_CONST * 10.**(-2.0*ZP)* TT**CONT_2
C--------------------------------------------------------------
C     THE SECOND TERM FOR H2O-H2O COLLISION (GODON ET AL 1991)
C--------------------------------------------------------------
      B = B+SC_CONST*1.595e-8/8.592e-10*TT**(6.18)*10.**(-2*ZP)*VMR_H2O
      ABSC = ABSC + B * VMR_H2O
      ABSC = ABSC * 1.e-3       ! converted from 1/km to 1/m

      RETURN
      END

C--------------------------------------------------------------
C     DEFINE LINE-SHAPE FUNCTIONS
C--------------------------------------------------------------
	real*8 function myshape(v0,ff,yy,twth0)
	implicit none
        real pi
        parameter (pi=3.1415926)
	real*8 voffm		! frequency difference
	real*8 voffp		! frequency difference
	real*8 ff		! frequency in MHz
	real*8 v0		! line frequency in MHz
	real twth0		! line width in MHz/mb
	real yy 		! interference coeff.
	real*8 twthm, twthp
c
	  voffm = v0 - ff
	  voffp = v0 + ff
	  twthm = twth0 - voffm*yy*1.d0
	  twthp = twth0 - voffp*yy*1.d0
c... Gross function
c	myshape = 4.*ff*v0*twth0/pi/((voffm*voffp)**2+4*(ff*twthm)**2)
c	return

c... Van Vleck-Weisskopf
 	myshape = twthm/(voffm**2 + twth0**2) + 
     >      twthp/(voffp**2 + twth0**2)
 	myshape = myshape * (ff/v0)/pi
c...
	return
	end

! $Log: get_beta.f,v      



