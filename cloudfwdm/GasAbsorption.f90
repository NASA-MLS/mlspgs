! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GasAbsorption

! -------------------------------------------------------------------------  
! COMPUTE ATMOSPHERIC GASES ABSORPTION COEFFICIENTS
! -------------------------------------------------------------------------

      use MLSCommon, only: r8 
      use WaterVapor, only: RHtoEV
      use RHIFromH2O, only: RHIFromH2O_Factor
      IMPLICIT NONE
      Private
      Public :: GET_BETA

 !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm =                          &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName=                       &
    "$RCSfile$"
  private :: not_used_here 
 !---------------------------------------------------------------------------
      
      REAL, private, parameter :: PI = 3.1415926

contains

! -------------------------------------------------  GET_BETA  -----
      SUBROUTINE GET_BETA(QLG,V0,GSE,IST,WTH,NTH,DELTA,N1,GAMMA,N2, &
                 &        NMOL,NCNT,T,PB,F,RH,VMR,ABSC,NS )
!                &        MOL,NMOL,NCNT,T,PB,F,RH,VMR,ABSC,NS )

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      LATEST UPDATE: J.JIANG, MAY 18, 2001
!==============================================================

      INCLUDE 'spectra.f9h' 

      INTEGER :: NS
      REAL :: QTP(3)                          ! TEMPERATURE ON WHICH QLG 
      DATA QTP /300.0,225.0,150.0/            ! ARE GIVEN
                                              !                 | A GIVEN TEMP
      REAL(r8) :: QRAT                        ! INTERPOLATED QLG RATIO AT 
      REAL(r8) :: F                           ! FREQUENCY IN GHz
      REAL(r8) :: T                           ! TEMPERATURE (K)
      REAL(r8) :: P                           ! DRY AIR PARTIAL PRESSURE (hPa)
      REAL(r8) :: PB                          ! TOTAL AIR PRESSURE (hPa)
      REAL(r8) :: VP                          ! VAPOR PARTIAL PRESSURE (hPa)
      REAL(r8) :: VMR(NS-1)                   ! MINOR SPECIES 1=O3, 2=N2O, 3=HNO3
      REAL(r8) :: VMR_H2O                     ! H2O VOLUME MIXING RATIO
      REAL(r8) :: VMR_O2                      ! O2 VOLUME MIXING RATIO
      REAL(r8) :: B                           ! BETA (1/m/ppv)
      REAL(r8) :: RH                          ! Relative Humidity if values > 1
					      ! or H2O vmr if values < 1
      REAL(r8) :: ABSC                        ! ABSORPTION COEFFICIENT (1/m)
                                              ! 1-O2
                                              ! 2-H2O
                                              ! 3-O_18_O
                                              ! 4-H2O_18
                                              ! 5-O3
                                              ! 6-N2O
                                              ! 7-HNO3

      REAL(r8) :: CONT_1,CONT_2,CONT_3        ! CONTINUUM ABSORPTION COEFFICIENTS
      REAL(r8) :: SC_CONST
      REAL(r8) :: SD,G0                       ! DEBY CONTRIBUTION (LIEBE 1989)
      REAL(r8) :: PS, NPS                     ! PRESSURE SHIFT PARAMETERS
     
!      REAL :: PI
!      PARAMETER (PI=3.1415926)

      REAL(r8) :: ZP,YY,TT,TWTH0,DWTH0        ! WORKING SPACE
      INTEGER :: I, J, IMOL

!      REAL(r8) :: MYSHAPE, FF, MYRESULT
!      EXTERNAL MYSHAPE

       REAL(r8) :: FF

!------------------------------------------------------------------------

      IF(RH .LT. 1._r8) THEN	       ! RH IS WATER VAPOR MIXING RATIO

         VMR_H2O = RH                  
         VP=VMR_H2O*PB                 ! VP IS VAPOR PRESSURE, PB IS TOTAL
         P=PB-VP                       ! PRESSURE, P IS DRY-AIR PRESSURE

      ELSE IF (RH .GE. 1._r8) THEN     ! RH HERE IS RELATIVE HUMIDITY 
!         CALL RHtoEV(T,RH,VP)    
!         P = PB-VP
!         VMR_H2O = VP/(max(1.e-9_r8, P))
         VMR_H2O = RHIFromH2O_Factor (T, -log10(PB), 0, .true.)*RH
         !                                           *
         ! * optional 0 will return mixing ratio in parts per 1
         ! * optional 6 will return mixing ratio in units of ppmv
         VP=PB*VMR_H2O
         P =PB-VP              !approximate dry air pressure
      END IF

         VMR_H2O = max(1.e-9_r8, vmr_h2o)

      VMR_O2 = 0.209476_r8
      FF     = F*1000.                   ! CONVERT F TO MHz
      ZP     = -LOG10( max(1.e-9_r8,P) ) ! CONVERT P(hPa) TO ZP
      TT     = 300._r8/T

      ABSC=0._r8
      I=0

!-------------------------------------
!     LINE EMISSION
!-------------------------------------
      DO IMOL=1,NMOL

         B=0._r8
         
         DO J=1,NCNT(IMOL)
            I=I+1

            IF(IMOL.GE.3.AND.ABS(V0(i)-FF).GT.10000.)  cycle ! GOTO 100 SAVE CPU
            
            PS  = 0.00_r8
            NPS = 0.00_r8

         ! Determine pressure shift parameters
            select case (IMOL)
            case (1)                    ! O2
              IF ( V0(i) .eq. 118750.3410_r8 ) THEN
                PS   = -0.140_r8
                NPS  = 1.36_r8
              END IF
            case (2)                    ! H2O
              IF ( V0(i) .eq. 183310.0910_r8) THEN 
                PS   = -0.160_r8
                NPS  = 1.375_r8
              ENDIF
            case (4)                    ! H2O_18
              IF ( V0(i) .eq. 203407.5200_r8 ) THEN
                PS   = -0.160_r8
                NPS  = 1.375_r8 
              END IF
            end select
            v01(i) = 0.0_r8                     | Pressure Shift effects
            V01(i) = V0(i) + PS * P * (TT**NPS) ! Include Hugh Pumphrey's

            IF(T .LE. QTP(2)) THEN
               QRAT = (QLG(2,I)-QLG(1,I))+&
                & (QLG(3,I)-QLG(2,I))*(T-QTP(2))/(QTP(3)-QTP(2))
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
            select case (IMOL)
            case (1)
              DWTH0 = 3.58e-7*SQRT(T/32.)*FF         ! O2
            case (2)
              DWTH0 = 3.58e-7*SQRT(T/18.)*FF         ! H2O
            case (3)
              DWTH0 = 3.58e-7*SQRT(T/34.)*FF         ! O_18_0
            case (4)
              DWTH0 = 3.58e-7*SQRT(T/20.)*FF         ! H2O_18
            case (5)
              DWTH0 = 3.58e-7*SQRT(T/48.)*FF         ! O3
            case (6)
              DWTH0 = 3.58e-7*SQRT(T/44.)*FF         ! N2O
            case (7)
              DWTH0 = 3.58e-7*SQRT(T/63.)*FF         ! HNO3
            end select

!--------------------------------------------------
!           WHITING'S APPROXIMATION FOR THE VOIGT
!--------------------------------------------------

            TWTH0 = TWTH0/2+SQRT(TWTH0**2/4+DWTH0**2)

!---------------------------------------------------------------
!	TO SAVE COMPUTING TIME CHECK IF FREQ IS CLOSE TO O3 LINES
!----------------------------------------------------------------

            B = B+2.30549e09_r8*10.0**(IST(I)-QRAT+GSE(I)*   &
     &		(1._r8/300._r8 - 1._r8/T)/1.600386_r8 - ZP) * PI * &
     &		 MYSHAPE(V01(i),FF,YY,TWTH0) * FF/V01(I) *  &
     &		(1._r8 - EXP(-V01(i)/(20836.7_r8*T)))/           &
     &		(1._r8 - EXP(-V01(i)/(20836.7_r8*300.0_r8)))/T 

         ENDDO   

         B=MAX(0._r8,B)                 ! AVOID NEGATIVE LINE SHAPE     

            select case (IMOL)
            case (1)
              ABSC = ABSC + B*VMR_O2                   ! O2
            case (2)
              ABSC = ABSC + B*VMR_H2O                  ! H2O
            case (3)
              ABSC = ABSC + B*VMR_O2*0.00409524        ! O_18_O
            case (4)
              ABSC = ABSC + B*VMR_H2O*0.00204          ! H2O_18
            case (5)
              ABSC = ABSC + B*VMR(1)                   ! O3
            case (6)
              ABSC = ABSC + B*VMR(2)                   ! N2O
            end select

      ENDDO  

!----------------------------------------
!     DRY CONTINUUM (ATBD APPENDIX F)
!----------------------------------------

      B=(7.7e-10*EXP(-1.5e-3*(F/30)**2)+1.e-13*(60**2+(F/30)**2)  &
     &     *EXP(-1.e-4*(F/30)**2))*TT**1.7
!---------------------------
! from ATBD N2-N2 continuum                                     
!--------------------------
       ABSC = ABSC + B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5/1.8 !best fit to B2 & B6
!      ABSC = ABSC + B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5

!==============================================================================
! from ATBD N2-N2 continuum, 0.84 is the best fit to f15 (with N2 O2) sids at
! 640GHz. The Debye term is coded differently between the two models, which
! affects R2 and R3 mostly (check it later) 
!       ABSC = ABSC + 0.84*B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5

! The factor 1.8 is fix to match B2U of Bill's FWM without [o2, o2]
!      ABSC = ABSC + B*0.65*(P/1013.)**2*TT**2*(F/30)**2*1.e5/1.8  
!===============================================================================
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
!      CONT_1 = 1.28e-15 	! BILL'S VALIDATION PAPER
!      CONT_1 = 1.37e-15         ! UARS 203GHz  v5 

!      the following cont_1 is adjusted based on bill's 
      ! H2O_R2  continuum=[ 2.52300e-16, 3.62800 ]

! difference with Bill's FWM f15 Band2U~6K, B10L~1K

! case R1
      if(abs(f-127.) .lt. 14.) CONT_1 = 7.53e-16     ! wu's version
! case R2
      if(abs(f-192.) .lt. 16.) CONT_1 = 7.53e-16/1.3 ! best fit to R2
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

! -------------------------------------------------  myshape  -----
!--------------------------------------------------------------
!     DEFINE LINE-SHAPE FUNCTIONS
!--------------------------------------------------------------

        real(r8) function myshape(v0,ff,yy,twth0) result (myresult)

	real(r8) ::  voffm		! frequency difference
	real(r8) ::  voffp		! frequency difference
	real(r8) ::  ff		        ! frequency in MHz
	real(r8) ::  v0		        ! line frequency in MHz
	real(r8) ::  twth0		! line width in MHz/mb
	real(r8) ::  yy 		! interference coeff.
	real(r8) ::  twthm, twthp

	  voffm = 0.0_r8
	  voffp = 0.0_r8
	  twthm = 0.0_r8
	  twthp = 0.0_r8
          myresult = 0.0_r8

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

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GasAbsorption

! $Log$
! Revision 1.21  2003/04/03 17:13:40  dwu
! changes: N2O mass is 44; HNO3 mass is 63
!
! Revision 1.20  2003/04/03 16:49:18  dwu
! fix a bug: N2O mass should be 30; and add HNO3 in absorption calculation
!
! Revision 1.19  2003/02/04 19:25:52  jonathan
! Fix bug in using RHIFromH2O
!
! Revision 1.18  2003/02/03 20:28:45  dwu
! make wet continuum dependent on radiometer frequency
!
! Revision 1.17  2003/02/01 06:43:16  dwu
! some fixes
!
! Revision 1.16  2003/01/30 18:03:20  jonathan
! switch to use Bill's RHIFromH2O
!
! Revision 1.15  2003/01/23 00:19:09  pwagner
! Some cosmetic only (or so I hope) changes
!
! Revision 1.14  2002/12/18 16:11:07  jonathan
! minor changes
!
! Revision 1.13  2002/11/06 19:08:08  jonathan
! best fit to b2 and b6
!
! Revision 1.12  2002/11/06 18:20:03  jonathan
! add N2O
!
! Revision 1.11  2002/10/08 17:08:07  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.10  2002/08/22 00:14:42  jonathan
! upgrade to include more molecules
!
! Revision 1.9  2002/08/19 22:22:03  jonathan
! debug stuff
!
! Revision 1.8  2002/08/08 22:46:44  jonathan
! newly improved version
!
! Revision 1.7  2002/05/06 22:35:22  jonathan
! add Bill's version for H2O
!
! Revision 1.6  2002/04/30 18:15:22  jonathan
! change CONT_1
!
! Revision 1.5  2001/11/14 00:40:26  jonathan
! minor changes
!
! Revision 1.4  2001/09/21 15:51:37  jonathan
! modified F95 version
!
