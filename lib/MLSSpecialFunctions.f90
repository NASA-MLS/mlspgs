! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSSpecialFunctions              ! Some special functions
!=============================================================================

  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none

  private
  public :: gamma
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface gamma
    module procedure sgamma, dgamma
  end interface
contains

      double precision function DGAMMA(X)
!>> 1996-03-30 DGAMMA Krogh  Added external statement.
!>> 1994-10-20 DGAMMA Krogh  Changes to use M77CON
!>> 1991-10-21 DGAMMA CLL Eliminated DGAM1 as a separate subroutine.
!>> 1991-01-16 DGAMMA Lawson  Replaced D2MACH/R2MACH with DGAM1.
!>> 1985-08-02 DGAMMA Lawson  Initial code.
!--D replaces "?": ?GAMMA, ?ERM1, ?ERV1
!
! ----------------------------------------------------------------------
!
!  THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A double precision
!      ARGUMENT X. PERMITS NEGATIVE AS WELL AS POSITIVE X. NOTE
!      THAT THE GAMMA FUNCTION HAS POLES AT ZERO AND AT NEGATIVE
!      ARGUMENTS. COMPUTATION IS BASED ON AN ALGORITHN OUTLINED IN
!      W.J.CODY, 'AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!      FUNCTIONS', LECTURE NOTES IN MATHEMATICS, 506, NUMERICAL ANALYSIS
!      DUNDEE, 1975, G. A. WATSON (ED.),SPRINGER VERLAG, BERLIN, 1976.
!      THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!      FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS. COEFFICIENTS
!      FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!      THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM HART, ET. AL.,
!      COMPUTER APPROXIMATIONS, WILEY AND SONS, NEW YORK, 1968.
!      LOWER ORDER APPROXIMATIONS CAN BE SUBSTITUTED FOR THESE ON
!      MACHINES WITH LESS PRECISE ARITHMETIC.
!
!  Designed & programmed by W.J.CODY, Argonne National Lab.,1982.
!  Minor changes for the JPL library by C.L.LAWSON & S.CHAN,JPL,1983.
!
!***********************************************************************
!
!-- Begin mask code changes
!  EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
!  EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!           1.0 + EPS .GT. 1.0  (EPS = [D/R]1MACH(4).)
!  XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER.
!           XINF = [D/R]1MACH(2).
!  XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!           both XMININ and 1/XMININ are representable.  This is
!           the larger of [D/R]1MACH(1) and 1.0 / [D/R]1MACH(2).
!  XGBIG  - A value such that    Gamma(XGBIG) = 0.875 * XINF.
!           (Computed and used in [D/S]GAMMA.)
!  XLBIG  - A value such that LogGamma(XLBIG) = 0.875 * XINF.
!           (Computed and used in [D/S]LGAMA.)
!
!      Values of XINF, XGBIG, and XLBIG for some machines:
!
!        XINF              XGBIG     XLBIG       Machines
!
!  2**127  = 0.170e39      34.81  0.180e37     Vax SP & DP; Unisys SP
!  2**128  = 0.340e39      35.00  0.358e37     IEEE SP
!  2**252  = 0.723e76      57.54  0.376e74     IBM30xx DP
!  2**1023 = 0.899e308    171.46  0.112e306    Unisys DP
!  2**1024 = 0.180e309    171.60  0.2216e306   IEEE DP
!  2**1070 = 0.126e323    177.78  0.1501e320   CDC/7600 SP
!  2**8191 = 0.550e2466   966.94  0.8464e2462  Cray SP & DP
!-- End mask code changes
!
!***********************************************************************
!
!  ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR. THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!  AUTHOR:W. J. CODY
!         APPLIED MATHMATICS DIVISION
!         ARGONNE NATIONAL LABORATORY
!         ARGONNE, IL 60439
!
!  LATEST MODIFICATION by Cody: MAY 18, 1982
!
!     ------------------------------------------------------------------
      double precision C(7), CONST, DEL, EPS, F, FACT, FP, HALF
      double precision ONE,P(8), PI,Q(8), RES,C1
      double precision SUM,TEMP, TWELVE,TWO
      double precision X,X1, X2, XGBIG,XDEN,XINF,XMININ,XNUM
      double precision Y,Y1,YSQ,Z,ZERO
      integer I,J,N
      logical PARITY
!
      save EPS, XGBIG, XMININ, XINF
!
      parameter( ONE = 1.0D0, HALF = 0.5D0, TWO = 2.0d0)
      parameter( ZERO = 0.0D0, TWELVE = 12.0D0)
!
!                      C1 = LOG base e of SQRT(2*PI)
!
      parameter( C1 = 0.9189385332046727417803297D0)
!
      parameter( PI = 3.1415926535897932384626434D0)
!
      data XINF/0.0D0/
!
! ----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
! ----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,&
            -3.79804256470945635097577D+2,6.29331155312818442661052D+2,&
            8.66966202790413211295064D+2,-3.14512729688483675254357D+4,&
            -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,&
           -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,&
             2.25381184209801510330112D+4,4.75584627752788110767815D+3,&
           -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
! ----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
! ----------------------------------------------------------------------
      DATA C/-1.910444077728D-03,8.4171387781295D-04,&
          -5.952379913043012D-04,7.93650793500350248D-04,&
          -2.777777777777681622553D-03,8.333333333333333331554247D-02,&
           5.7083835261D-03/
! ----------------------------------------------------------------------
!
      IF (XINF .EQ. ZERO) THEN
        EPS = epsilon(x)   ! D1MACH(4)
        XINF = huge(x)     ! D1MACH(2)
         ! if(D1MACH(1) * D1MACH(2) .ge. ONE) then
         !   XMININ = D1MACH(1)
         ! else
         !   XMININ = ONE / D1MACH(2)
         ! endif
         XMININ = tiny(x)
!
!                         Compute XGBIG
!
!        XGBIG will satisfy Gamma(XGBIG) = 0.875 * XINF.
!        Use a Newton iteration and the following approximation for
!        the gamma function:
!        log(gamma(x)) ~ (x - .5)*log(x) - x + 0.5 * log(2 * PI)
!
         TEMP = log(0.875d0 * XINF)
         CONST = HALF * log(TWO * PI) - TEMP
         X1 = TEMP * 0.34d0
         do 40 J=1,7
            F = (X1-HALF)*log(X1) - X1 + CONST
            FP = ((X1-HALF)/X1)  + log(X1) - ONE
            DEL = -F/FP
            X2 = X1+DEL
            if(abs(DEL) .lt. 0.5d-5 * X2) go to 45
            X1 = X2
   40    continue
   45    continue
         XGBIG = X2
      END IF
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .GT. ZERO) GO TO 200
! ----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE OR ZERO
! ----------------------------------------------------------------------
      Y = -X
      J = INT(Y)
      RES = Y - dble(J)
      IF (RES .EQ. ZERO) GO TO 700
      IF (J .NE. (J/2)*2) PARITY = .TRUE.
      FACT = -PI / sin(PI*RES)
      Y = Y + ONE
! ----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
! ----------------------------------------------------------------------
  200 IF (Y .LT. EPS) GO TO 650
      IF (Y .GE. TWELVE) GO TO 300
      Y1 = Y
      IF (Y .GE. ONE) GO TO 210
! ----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
! ----------------------------------------------------------------------
      Z = Y
      Y = Y + ONE
      GO TO 250
! ----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
! ----------------------------------------------------------------------
  210 N= int(Y) - 1
      Y = Y - dble(N)
      Z = Y - ONE
! ----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
! ----------------------------------------------------------------------
  250 XNUM = ZERO
      XDEN = ONE
      DO 260 I = 1, 8
         XNUM = (XNUM + P(I)) * Z
         XDEN = XDEN * Z + Q(I)
  260 CONTINUE
      RES = (XNUM / XDEN + HALF) + HALF
      IF (Y .EQ. Y1) GO TO 900
      IF (Y1 .GT. Y) GO TO 280
! ----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE 0.0 .LT. ARGUMENT .LT. 1.0
! ----------------------------------------------------------------------
      RES = RES / Y1
      GO TO 900
! ----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE 2.0 .LT. 12.0
! ----------------------------------------------------------------------
  280 DO 290 I = 1, N
         RES = RES * Y
         Y = Y + ONE
  290 CONTINUE
      GO TO 900
! ----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0
! ----------------------------------------------------------------------
  300 IF (Y .GT. XGBIG) GO TO 720
      YSQ = Y * Y
      SUM = C(7)
      DO 350 I = 1, 6
         SUM = SUM / YSQ + C(I)
  350 CONTINUE
      SUM = ((SUM/Y + C1) - Y) + (Y - HALF) * log(Y)
      RES = exp(SUM)
      GO TO 900
! ----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
! ----------------------------------------------------------------------
  650 IF (Y .LT. XMININ) GO TO 740
      RES = ONE / Y
      GO TO 900
! ----------------------------------------------------------------------
!  RETURN FOR SINGULARITIES,EXTREME ARGUMENTS, ETC.
! ----------------------------------------------------------------------
  700 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "POLE AT 0 AND NEG INTEGERS" )
      GO TO 780
!
  720 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "arg SO LARGE VALUE WOULD OVERFLOW" )
      GO TO 780
!
  740 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "arg too near a sinularity--VALUE WOULD OVERFLOW" )
!
  780 DGAMMA = XINF
      GO TO 950
! ----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
! ----------------------------------------------------------------------
  900 IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
      DGAMMA = RES
  950 RETURN
    end function DGAMMA
      real             function SGAMMA(X)
!>> 1996-03-30 SGAMMA Krogh  Added external statement.
!>> 1994-10-20 SGAMMA Krogh  Changes to use M77CON
!>> 1991-10-21 SGAMMA CLL Eliminated DGAM1 as a separate subroutine.
!>> 1991-01-16 SGAMMA Lawson  Replaced D2MACH/R2MACH with DGAM1.
!>> 1985-08-02 SGAMMA Lawson  Initial code.
!--S replaces "?": ?GAMMA, ?ERM1, ?ERV1
!
! ----------------------------------------------------------------------
!
!  THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A double precision
!      ARGUMENT X. PERMITS NEGATIVE AS WELL AS POSITIVE X. NOTE
!      THAT THE GAMMA FUNCTION HAS POLES AT ZERO AND AT NEGATIVE
!      ARGUMENTS. COMPUTATION IS BASED ON AN ALGORITHN OUTLINED IN
!      W.J.CODY, 'AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!      FUNCTIONS', LECTURE NOTES IN MATHEMATICS, 506, NUMERICAL ANALYSIS
!      DUNDEE, 1975, G. A. WATSON (ED.),SPRINGER VERLAG, BERLIN, 1976.
!      THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!      FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS. COEFFICIENTS
!      FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!      THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM HART, ET. AL.,
!      COMPUTER APPROXIMATIONS, WILEY AND SONS, NEW YORK, 1968.
!      LOWER ORDER APPROXIMATIONS CAN BE SUBSTITUTED FOR THESE ON
!      MACHINES WITH LESS PRECISE ARITHMETIC.
!
!  Designed & programmed by W.J.CODY, Argonne National Lab.,1982.
!  Minor changes for the JPL library by C.L.LAWSON & S.CHAN,JPL,1983.
!
!***********************************************************************
!
!-- Begin mask code changes
!  EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
!  EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!           1.0 + EPS .GT. 1.0  (EPS = [D/R]1MACH(4).)
!  XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER.
!           XINF = [D/R]1MACH(2).
!  XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!           both XMININ and 1/XMININ are representable.  This is
!           the larger of [D/R]1MACH(1) and 1.0 / [D/R]1MACH(2).
!  XGBIG  - A value such that    Gamma(XGBIG) = 0.875 * XINF.
!           (Computed and used in [D/S]GAMMA.)
!  XLBIG  - A value such that LogGamma(XLBIG) = 0.875 * XINF.
!           (Computed and used in [D/S]LGAMA.)
!
!      Values of XINF, XGBIG, and XLBIG for some machines:
!
!        XINF              XGBIG     XLBIG       Machines
!
!  2**127  = 0.170e39      34.81  0.180e37     Vax SP & DP; Unisys SP
!  2**128  = 0.340e39      35.00  0.358e37     IEEE SP
!  2**252  = 0.723e76      57.54  0.376e74     IBM30xx DP
!  2**1023 = 0.899e308    171.46  0.112e306    Unisys DP
!  2**1024 = 0.180e309    171.60  0.2216e306   IEEE DP
!  2**1070 = 0.126e323    177.78  0.1501e320   CDC/7600 SP
!  2**8191 = 0.550e2466   966.94  0.8464e2462  Cray SP & DP
!-- End mask code changes
!
!***********************************************************************
!
!  ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR. THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!  AUTHOR:W. J. CODY
!         APPLIED MATHMATICS DIVISION
!         ARGONNE NATIONAL LABORATORY
!         ARGONNE, IL 60439
!
!  LATEST MODIFICATION by Cody: MAY 18, 1982
!
!     ------------------------------------------------------------------
      real             C(7), CONST, DEL, EPS, F, FACT, FP, HALF
      real             ONE,P(8), PI,Q(8), RES,C1
      real             SUM,TEMP, TWELVE,TWO
      real             X,X1, X2, XGBIG,XDEN,XINF,XMININ,XNUM
      real             Y,Y1,YSQ,Z,ZERO
      integer I,J,N
      logical PARITY
!
      save EPS, XGBIG, XMININ, XINF
!
      parameter( ONE = 1.0E0, HALF = 0.5E0, TWO = 2.0e0)
      parameter( ZERO = 0.0E0, TWELVE = 12.0E0)
!
!                      C1 = LOG base e of SQRT(2*PI)
!
      parameter( C1 = 0.9189385332046727417803297E0)
!
      parameter( PI = 3.1415926535897932384626434E0)
!
      data XINF/0.0E0/
!
! ----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
! ----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,&
            -3.79804256470945635097577E+2,6.29331155312818442661052E+2,&
            8.66966202790413211295064E+2,-3.14512729688483675254357E+4,&
            -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,&
           -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,&
             2.25381184209801510330112E+4,4.75584627752788110767815E+3,&
           -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
! ----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
! ----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,&
          -5.952379913043012E-04,7.93650793500350248E-04,&
          -2.777777777777681622553E-03,8.333333333333333331554247E-02,&
           5.7083835261E-03/
! ----------------------------------------------------------------------
!
      IF (XINF .EQ. ZERO) THEN
        EPS = epsilon(x)   ! R1MACH(4)
        XINF = huge(x)     ! R1MACH(2)
         ! if(R1MACH(1) * R1MACH(2) .ge. ONE) then
         !   XMININ = R1MACH(1)
         !  else
         !   XMININ = ONE / R1MACH(2)
         ! endif
         XMININ = tiny(x)
!
!                         Compute XGBIG
!
!        XGBIG will satisfy Gamma(XGBIG) = 0.875 * XINF.
!        Use a Newton iteration and the following approximation for
!        the gamma function:
!        log(gamma(x)) ~ (x - .5)*log(x) - x + 0.5 * log(2 * PI)
!
         TEMP = log(0.875e0 * XINF)
         CONST = HALF * log(TWO * PI) - TEMP
         X1 = TEMP * 0.34e0
         do 40 J=1,7
            F = (X1-HALF)*log(X1) - X1 + CONST
            FP = ((X1-HALF)/X1)  + log(X1) - ONE
            DEL = -F/FP
            X2 = X1+DEL
            if(abs(DEL) .lt. 0.5e-5 * X2) go to 45
            X1 = X2
   40    continue
   45    continue
         XGBIG = X2
      END IF
      PARITY = .FALSE.
      FACT = ONE
      N = 0
      Y = X
      IF (Y .GT. ZERO) GO TO 200
! ----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE OR ZERO
! ----------------------------------------------------------------------
      Y = -X
      J = INT(Y)
      RES = Y - real(J)
      IF (RES .EQ. ZERO) GO TO 700
      IF (J .NE. (J/2)*2) PARITY = .TRUE.
      FACT = -PI / sin(PI*RES)
      Y = Y + ONE
! ----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
! ----------------------------------------------------------------------
  200 IF (Y .LT. EPS) GO TO 650
      IF (Y .GE. TWELVE) GO TO 300
      Y1 = Y
      IF (Y .GE. ONE) GO TO 210
! ----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
! ----------------------------------------------------------------------
      Z = Y
      Y = Y + ONE
      GO TO 250
! ----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
! ----------------------------------------------------------------------
  210 N= int(Y) - 1
      Y = Y - real(N)
      Z = Y - ONE
! ----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
! ----------------------------------------------------------------------
  250 XNUM = ZERO
      XDEN = ONE
      DO 260 I = 1, 8
         XNUM = (XNUM + P(I)) * Z
         XDEN = XDEN * Z + Q(I)
  260 CONTINUE
      RES = (XNUM / XDEN + HALF) + HALF
      IF (Y .EQ. Y1) GO TO 900
      IF (Y1 .GT. Y) GO TO 280
! ----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE 0.0 .LT. ARGUMENT .LT. 1.0
! ----------------------------------------------------------------------
      RES = RES / Y1
      GO TO 900
! ----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE 2.0 .LT. 12.0
! ----------------------------------------------------------------------
  280 DO 290 I = 1, N
         RES = RES * Y
         Y = Y + ONE
  290 CONTINUE
      GO TO 900
! ----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0
! ----------------------------------------------------------------------
  300 IF (Y .GT. XGBIG) GO TO 720
      YSQ = Y * Y
      SUM = C(7)
      DO 350 I = 1, 6
         SUM = SUM / YSQ + C(I)
  350 CONTINUE
      SUM = ((SUM/Y + C1) - Y) + (Y - HALF) * log(Y)
      RES = exp(SUM)
      GO TO 900
! ----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
! ----------------------------------------------------------------------
  650 IF (Y .LT. XMININ) GO TO 740
      RES = ONE / Y
      GO TO 900
! ----------------------------------------------------------------------
!  RETURN FOR SINGULARITIES,EXTREME ARGUMENTS, ETC.
! ----------------------------------------------------------------------
  700 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "POLE AT 0 AND NEG INTEGERS" )
      GO TO 780
!
  720 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "arg SO LARGE VALUE WOULD OVERFLOW" )
      GO TO 780
!
  740 CALL MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "arg too near a sinularity--VALUE WOULD OVERFLOW" )
!
  780 SGAMMA = XINF
      GO TO 950
! ----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
! ----------------------------------------------------------------------
  900 IF (PARITY) RES = -RES
      IF (FACT .NE. ONE) RES = FACT / RES
      SGAMMA = RES
  950 RETURN
      end function sgamma
!=============================================================================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MLSSpecialFunctions
!=============================================================================

!
! $Log$
! Revision 2.1  2004/06/02 21:03:16  pwagner
! First commit
!
