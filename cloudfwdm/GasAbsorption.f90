! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GasAbsorption

! -------------------------------------------------------------------------  
! COMPUTE ATMOSPHERIC GASES ABSORPTION COEFFICIENTS
! -------------------------------------------------------------------------

      use MLSCommon, only: r8
      IMPLICIT NONE
      Private
      Public :: GET_BETA

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
 !---------------------------------------------------------------------------
      
contains

      SUBROUTINE GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2, &
                 &        MOL,NMOL,NCNT,T,PB,F,RH,VMR,ABSC,NS)

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      LATEST UPDATE: J.JIANG, MAY 18, 2001
!==============================================================

      INCLUDE 'spectra.f9h' 

      INTEGER :: NS
      REAL(r8) :: QTP(3)                      ! TEMPERATURE ON WHICH QLG ARE GIVEN
      DATA QTP /300.0,225.0,150.0/

      REAL(r8) :: QRAT                        ! INTERPOLATED QLG RATIO AT A GIVEN TEMP
      REAL(r8) :: F                           ! FREQUENCY IN GHz
      REAL(r8) :: T                           ! TEMPERATURE (K)
      REAL(r8) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa)
      REAL(r8) :: PB                          ! TOTAL AIR PRESSURE (hPa)
      REAL(r8) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)
      REAL(r8) :: VMR(NS)                      ! MINOR SPECIES 1-O3
      REAL(r8) :: VMR_H2O                     ! H2O VOLUME MIXING RATIO
      REAL(r8) :: VMR_O2                      ! O2 VOLUME MIXING RATIO
      REAL(r8) :: B                           ! BETA (1/m/ppv)
      REAL(r8) :: RH                          ! Relative Humidity
      REAL(r8) :: ABSC                        ! ABSORPTION COEFFICIENT (1/m)
                                              ! 1-O2
                                              ! 2-H2O
                                              ! 3-O_18_O
                                              ! 4-H2O_18

      REAL(r8) :: CONT_1,CONT_2,CONT_3        ! CONTINUUM ABSORPTION COEFFICIENTS
      REAL(r8) :: SC_CONST
      REAL(r8) :: SD,G0                       ! DEBY CONTRIBUTION (LIEBE 1989)

      REAL :: PI
      PARAMETER (PI=3.1415926)
      REAL(r8) :: ZP,YY,TT,TWTH0,DWTH0        ! WORKING SPACE
      INTEGER :: I, J, IMOL

!      REAL(r8) :: MYSHAPE, FF, MYRESULT
!      EXTERNAL MYSHAPE

       REAL(r8) :: FF

!------------------------------------------------------------------------

      IF(RH .NE. 100._r8) THEN
         VMR_H2O = RH                  ! PH IS WATER VAPOR MIXING RATIO
         VP=VMR_H2O*PB                 ! VP IS VAPOR PRESSURE, PB IS TOTAL
         P=PB-VP                       ! PRESSURE, P IS DRY-AIR PRESSURE
      ELSE IF (RH .EQ. 100._r8) THEN
!         CALL RHtoEV(PB,T,100._r8,VP)     ! RH HERE IS 100% RELATIVE HUMIDITY 
         CALL RHtoEV(T,100._r8,VP)     ! RH HERE IS 100% RELATIVE HUMIDITY 
         P = PB-VP
         VMR_H2O = VP/(max(1.e-9_r8, P))
      END IF

      VMR_O2 = 0.209476_r8
      FF     = F*1000.                 ! CONVERT F TO MHz
      ZP     = -LOG10( max(1.e-9_r8,P) ) ! CONVERT P(hPa) TO ZP
      TT     = 300._r8/T

      ABSC=0.
      I=0

!      print*,'check'

!-------------------------------------
!     LINE EMISSION
!-------------------------------------
      DO IMOL=1,NMOL
         B=0.
         DO J=1,NCNT(IMOL)
            I=I+1

            IF(IMOL.GE.3.AND.ABS(V0(i)-FF).GT.10000.) GOTO 100 ! SAVE CPU

            IF(T .LE. QTP(2)) THEN
!               QRAT = (QLG(3,I)-QLG(1,I))*(T-QTP(2))/(QTP(3)-QTP(2))
               QRAT = (QLG(2,I)-QLG(1,I))+(QLG(3,I)-QLG(2,I))*(T-QTP(2))/(QTP(3)-QTP(2))
            ELSE
               QRAT = (QLG(2,I)-QLG(1,I))*(T-QTP(1))/(QTP(2)-QTP(1))
            ENDIF

            YY=10.**(-ZP)*(DELTA(I)*TT**N1(I)+GAMMA(I)*TT**N2(I))

!-------------------------------------
!           PRESSURE WIDTH
!-------------------------------------

            TWTH0=WTH(I)*10.**(-ZP)*TT**NTH(I)

!-------------------------------------
!           DOPPLER WIDTH
!-------------------------------------

            IF(IMOL.EQ.1) DWTH0 = 3.58e-7*SQRT(T/32.)*FF
            IF(IMOL.EQ.2) DWTH0 = 3.58e-7*SQRT(T/18.)*FF
            IF(IMOL.EQ.3) DWTH0 = 3.58e-7*SQRT(T/34.)*FF
            IF(IMOL.EQ.4) DWTH0 = 3.58e-7*SQRT(T/20.)*FF
            IF(IMOL.EQ.5) DWTH0 = 3.58e-7*SQRT(T/48.)*FF
            
!--------------------------------------------------
!           WHITING'S APPROXIMATION FOR THE VOIGT
!--------------------------------------------------

            TWTH0 = TWTH0/2+SQRT(TWTH0**2/4+DWTH0**2)

!---------------------------------------------------------------
!	TO SAVE COMPUTING TIME CHECK IF FREQ IS CLOSE TO O3 LINES
!----------------------------------------------------------------
!            print*,IST(I),GSE(I),V0(I), FF, YY, TWTH0, T
!            print*,ZP, MYSHAPE(V0(i),FF,YY,TWTH0)   

            B = B+2.30549e09_r8*10.0**(IST(I)-QRAT+GSE(I)*   &
     &		(1._r8/300._r8 - 1._r8/T)/1.600386_r8 - ZP) * PI * &
     &		 MYSHAPE(V0(i),FF,YY,TWTH0) * FF/V0(I) *  &
     &		(1._r8 - EXP(-V0(i)/(20836.7_r8*T)))/           &
     &		(1._r8 - EXP(-V0(i)/(20836.7_r8*300.0_r8)))/T 

 100     CONTINUE

         ENDDO   

         B=MAX(0._r8,B)                 ! AVOID NEGATIVE LINE SHAPE     

         IF(IMOL .EQ. 1) ABSC = ABSC + B*VMR_O2                   ! O2      
         IF(IMOL .EQ. 2) ABSC = ABSC + B*VMR_H2O                  ! H2O
         IF(IMOL .EQ. 3) ABSC = ABSC + B*VMR_O2*0.00409524        ! O_18_O
         IF(IMOL .EQ. 4) ABSC = ABSC + B*VMR_H2O*0.00204          ! H2O_18
         IF(IMOL .EQ. 5) ABSC = ABSC + B*VMR(1)                   ! O3

      ENDDO  

!----------------------------------------
!     DRY CONTINUUM (ATBD APPENDIX F)
!----------------------------------------

      B=(7.7e-10*EXP(-1.5e-3*(F/30)**2)+1.e-13*(60**2+(F/30)**2)  &
     &     *EXP(-1.e-4*(F/30)**2))*TT**1.7

      ABSC = ABSC + B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5

!     CONT_1=1.4e-10*(1-1.2e-5*F**1.5)     ! LIEBE 1989
!     CONT_1 = 1.4e-12/(1+1.93e-5*F**1.5)  ! LIEBE 1993
!
!      CONT_1 = 1.4e-12/(1+1.93e-5*F**1.5)*1.28*1.17  ! MATCH THE UARS FIT
!
!      B=CONT_1 * F * (P/10.)**2*TT**3.5
!      ABSC=ABSC + B*.182*F/4.343 
!------------------------------------------
!     DEBYE TERM FROM LIEBE
!------------------------------------------
      SD = 6.14e-4*P/10.*TT**2
      G0 = 5.6e-3*(1+1.1*VMR_H2O)*P/10.*TT
      B = SD*F/(G0 * (1.+ (F/G0)**2))		! B is Npp in LIEBE
      ABSC = ABSC + B*.182*F/4.343              ! CONVERTED TO 1/km 

!----------------------------------------
!     WET CONTINUUM 
!----------------------------------------
!     CONT_1 = 1.28e-15 	! BILL'S VALIDATION PAPER
      CONT_1 = 7.53e-16 
      CONT_2 = 4.20
      CONT_3 = 0.00
      SC_CONST = CONT_1 * FF**2 * EXP(-CONT_3 * FF**2)
      B = SC_CONST * 10.**(-2.0*ZP)* TT**CONT_2
!--------------------------------------------------------------
!     THE SECOND TERM FOR H2O-H2O COLLISION (GODON ET AL 1991)
!--------------------------------------------------------------
      B = B+SC_CONST*1.595e-8/8.592e-10*TT**(6.18)*10.**(-2*ZP)*VMR_H2O
      ABSC = ABSC + B * VMR_H2O
      ABSC = ABSC * 1.e-3       ! converted from 1/km to 1/m

      END SUBROUTINE GET_BETA

!--------------------------------------------------------------
!     DEFINE LINE-SHAPE FUNCTIONS
!--------------------------------------------------------------
!J	real(r8) function myshape(v0,ff,yy,twth0)

        real(r8) function myshape(v0,ff,yy,twth0) result (myresult)
!        use MLSCommon, only: r8        
	implicit none
!        real(r8):: myresult     !J

        real :: pi
        parameter (pi=3.1415926)
	real(r8) ::  voffm		! frequency difference
	real(r8) ::  voffp		! frequency difference
	real(r8) ::  ff		        ! frequency in MHz
	real(r8) ::  v0		        ! line frequency in MHz
	real(r8) ::  twth0		! line width in MHz/mb
	real(r8) ::  yy 		! interference coeff.
	real(r8) ::  twthm, twthp
!
	  voffm = v0 - ff
	  voffp = v0 + ff
	  twthm = twth0 - voffm*yy*1.
	  twthp = twth0 - voffp*yy*1.
!... Gross function
!	myshape = 4.*ff*v0*twth0/pi/((voffm*voffp)**2+4*(ff*twthm)**2)
!	return

!... Van Vleck-Weisskopf
!J 	myshape = twthm/(voffm**2 + twth0**2) + &
!J     &      twthp/(voffp**2 + twth0**2)
!J 	myshape = myshape * (ff/v0)/pi

  	myresult = twthm/(voffm**2 + twth0**2) + &
       &      twthp/(voffp**2 + twth0**2)
 	 myresult = myresult * (ff/v0)/pi

       return
       end function myshape

end module GasAbsorption

! $Log: GasAbsorption.f90,v      







