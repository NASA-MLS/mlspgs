module Cadre_m

! The original version of this procedure was by Carl de Boor.  It was
! written in Fortran 66.  That version was converted to Fortran 90 free
! form, control constructs eliminated some GOTO statements, DO
! constructs were better structured, all variables were declared,
! intents were specified for dummy arguments, Hollerith edit descriptors
! were eliminated, and a reverse communication version was created, by
! Van Snyder in July 2009.

  implicit NONE
  private
  public :: Cadre, Cadre_Reverse

  integer, parameter :: RK = kind(0.0d0)

!---------------------------- RCS Module Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

      FUNCTION CADRE (F, A, B, AERR, RERR, LEVEL, ERROR, IFLAG)
!  THIS FUNCTION RETURNS AN ESTIMATE *CADRE* FOR THE NUMBER
!     INT = INTEGRAL OF *F*(X) FROM *A* TO *B*
!  WHICH HOPEFULLY SATISFIES
!     ABS(INT - *CADRE*) <= max(*AERR*, *RERR* TIMES ABS(INT)).
!     THE PROGRAM USES CAUTIOUS ADAPTIVE ROMBERG EXTRAPOLATION.
!  IN THIS SCHEME, THE INTEGRAL IS CALCULATED AS THE SUM OF INTEGRALS
!  OVER SUITABLY SMALL SUBINTERVALS. ON EACH SUBINTERVAL, AN ESTIMATE
!  *VINT*, WITH ESTIMATED ABSOLUTE ERROR *ERRER*, IS FOUND BY CAUTIOUS
!  ROMBERG EXTRAPOLATION. IF *ERRER* IS SMALL ENOUGH, *VINT* IS ACCEPT
!  ED AND ADDED TO *CADRE*, AND *ERRER* IS ADDED TO *ERROR*. OTHERWISE
!  THE SUBINTERVAL IS HALVED, AND EACH HALF IS CONSIDERED SEPARATELY,
!  INFORMATION ABOUT THE OTHER HALF BEING TEMPORARILY STACKED.

        interface
          function F ( X )
            integer, parameter :: RK = kind(0.0d0)
            real(rk) :: F
            real(rk), intent(in) :: X
          end function F
        end interface
        real(rk), intent(in) :: A, B
        real(rk), intent(in) :: Aerr, Rerr
        integer, intent(in) :: Level
        real(rk), intent(out) :: Error
        integer, intent(inout) :: Iflag
        real(rk) :: Cadre

!                   *****   INPUT   *****
!
!  F       THE NAME OF A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM
!          THIS NAME MUST APPEAR IN THE CALLING PROGRAM IN AN
!          EXTERNAL STATEMENT.
!  A,B     THE TWO ENDPOINTS OF THE INTERVAL OF INTEGRATION
!  AERR
!  RERR    DESIRED ABSOLUTE AND RELATIVE ERROR IN THE ANSWER
!  LEVEL   AN INTEGER INDICATING DESIRED LEVEL OF PRINTOUT
! <= 1,  NO PRINTOUT,
!           =  2,  FAILURE MESSAGE (IF ANY), AND LIST OF SINGULAR-
!                  ITIES ENCOUNTERED (IF ANY),
!           =  3,  IN ADDITION, ALL SUBINTERVALS CONSIDERED ARE LISTED
!                  TOGETHER WITH THE KIND OF REGULAR BEHAVIOUR FOUND
!                  (IF ANY),
!           =  4,  IN ADDITION, ALL RATIOS CONSIDERED ARE LISTED AS IS
!                  INFO ON WHICH DECISION PROCEDURE IS BASED,
! >= 5,  IN ADDITION, ALL T-TABLES ARE LISTED.
!
!                   *****   OUTPUT  *****
!
!  CADRE   ESTIMATE OF THE INTEGRAL, RETURNED VIA THE FUNCTION CALL,
!  ERROR   ESTIMATED BOUND ON THE ABSOLUTE ERROR OF THE NUMBER *CADRE*
!  IFLAG   AN INTEGER BETWEEN 1 AND 5 INDICATING WHAT DIFFICULTIES
!          WERE MET WITH, SPECIFICALLY
!          = 1, ALL IS WELL,
!          = 2, ONE OR MORE SINGULARITIES WERE SUCCESSFULLY HANDLED,
!          = 3, IN SOME SUBINTERVAL(S), THE ESTIMATE *VINT* WAS ACCEPT
!               ED MERELY BECAUSE *ERRER* WAS SMALL, EVEN THOUGH NO
!               REGULAR BEHAVIOUR COULD BE RECOGNIZED,
!          = 4, FAILURE, OVERFLOW OF STACK *TS* (THIS HAS NEVER HAPPEN
!               ED,  - SO FAR),
!          = 5, FAILURE, TOO SMALL A SUBINTERVAL IS REQUIRED. THIS MAY
!               BE DUE TO TOO MUCH NOISE IN THE FUNCTION (RELATIVE TO
!               THE GIVEN ERROR REQUIREMENTS) OR DUE TO A PLAIN ORNERY
!               INTEGRAND.
!  A VERY CAUTIOUS MAN WOULD ACCEPT *CADRE* ONLY IF IFLAG IS 1 OR 2.
!  THE MERELY REASONABLE MAN WOULD KEEP THE FAITH EVEN IF IFLAG IS 3.
!  THE ADVENTUROUS MAN IS QUITE OFTEN RIGHT IN ACCEPTING *CADRE*
!  EVEN IF IFLAG IS 4 OR 5.
!
!     *****   LIST OF MAJOR VARIABLES   *****
!
!  CUREST  BEST ESTIMATE SO FAR FOR
!             INT - (INTEGRAL OVER CURRENTLY CONSIDERED SUBINTERVAL).
!  FNSIZE  MAXIMUM AVERAGE FUNCTION SIZE SO FAR ENCOUNTERED,
!  ERRR    RELATIVE ERROR REQUIREMENT USED. DERIVED FROM INPUT *RERR*
!          AND CHOSEN TO LIE BETWEEN .1 AND 10 TIMES *TOLMCH*,
!  ERRA    = ABS(*AERR*)
!  STAGE   (MORE OR LESS) EQUAL TO 2 TO THE -(*ISTAGE*)
!  THESE FIVE QUANTITIES ARE USED IN THE DETERMINATION OF THE LOCAL
!  ERROR REQUIREMENT.
!
!  STEPMN  MINIMUM SUBINTERVAL LENGTH PERMITTED,
!  TS      STACK OF VALUES OF F(X) SO FAR COMPUTED BUT NOT YET
!          SUCCESSFULLY USED,
!  ISTAGE  AN INTEGER INDICATING THE HEIGHT OF THE STACK OF INTERVALS
!          YET TO BE PROCESSED.
!
!     *****   LIST OF PARAMETERS   *****
!
!  TOLMCH  DEPENDS ON THE LENGTH OF FLOATING POINT MANTISSA. SHOULD BE
!          ABOUT 1.E-7 FOR 27 BINARY BIT MANTISSA AND
!          ABOUT 1.E-13 FOR 48 BINARY BIT MANTISSA.
!  AITLOW  SHOULD BE SOMEWHAT GREATER THAN 1.
!  H2TOL,
!  AITTOL,
!  JUMPTL  TOLERANCES USED IN THE DECISION PROCESS TO RECOGNIZE
!          H**2 CONVERGENCE, X**ALPHA TYPE CONVERGENCE, OR
!          JUMP-TYPE CONVERGENCE OF THE TRAPEZOID SUMS.
!  MAXTS,
!  MAXTBL,
!  MXSTGE  ARE THE THREE DIFFERENT UPPER LIMITS FOR THE DIMENSION OF
!          THE VARIOUS ARRAYS.
!
!        *****   PROGRAM LAYOUT   *****
!
!          INITIALIZATION
!  5,6     BEGIN WORK ON NEXT SUBINTERVAL
!   9-14   GET NEXT TRAPEZOID SUM
!  15-19   GET RATIOS. PRELIMINARY DECISION PROCEDURE.
!  20-     ESTIMATE *VINT* ASSUMING SMOOTH INTEGRAND
!  30-     ESTIMATE *VINT* ASSUMING INTEGRAND HAS X**ALPHA TYPE
!          SINGULARITY
!  40-     NO LUCK WITH THIS TRAPEZOID SUM. GET NEXT ONE OR GET OUT.
!  50-     ESTIMATE *VINT* ASSUMING INTEGRAND HAS JUMP
!  60-     ESTIMATE *VINT* ASSUMING INTEGRAND IS STRAIGHT LINE
!  70-     ESTIMATE *VINT* ASSUMING VARIATION IN INTEGRAND
!          IS MOSTLY NOISE.
!  80-     INTEGRATION OVER CURRENT SUBINTERVAL SUCCESSFUL.
!          SET UP NEXT SUBINTERVAL, IF ANY, OR RETURN.
!  90-     INTEGRATION OVER CURRENT SUBINTERVAL NOT SUCCESSFUL.
!          MARK CURRENT SUBINTERVAL FOR SUBDIVISION AND SET UP
!          NEXT SUBINTERVAL.
!  900-    FAILURE.
!
      integer, parameter :: MAXTS = 2049
      integer, parameter :: MAXTBL = 10
      integer, parameter :: MXSTGE = 30
      real(rk) :: T (MAXTBL, MAXTBL), R (MAXTBL), AIT (MAXTBL),         &
     & DIF (MAXTBL), TS (MAXTS), IBEGS (MXSTGE), BEGIN (MXSTGE),        &
     & FINIS (MXSTGE), EST (MXSTGE)
      real(rk) :: LENGTH
      LOGICAL H2CONV, AITKEN, RIGHT, REGLAR, REGLSV (MXSTGE)
      real(rk), parameter :: JUMPTL = 0.01_rk
      real(rk), parameter :: TOLMCH = epsilon(1.0_rk)
      real(rk), parameter :: AITLOW = 1.1_rk
      real(rk), parameter :: H2TOL = 0.15_rk, AITTOL = 0.1_rk
!     "Randomly chosen" points to check for a straight line
      real(rk), parameter :: RN (4) =                                   &
     & (/ 0.71420053_rk, 0.34662815_rk, 0.843751_rk, 0.12633046_rk /)
      real(rk), parameter :: ALG4O2 = 0.3010299956639795_rk ! log10(2)
      DATA AIT (1) / 0.0_rk /
      real(rk) :: ABSI, ALPHA, ASTEP, BEG, CUREST, DIFF, END, ERGOAL
      real(rk) :: ERRA, ERRER, ERRR, FBEG, FBEG2, FEND, FEXTM1, FEXTRP
      real(rk) :: FN, FNSIZE, HOVN, H2NEXT, H2TFEX, PREVER, SING
      real(rk) :: SINGNX, SLOPE, STAGE, STEP, STEPMN, SUM, SUMABS
      real(rk) :: TABS, TABTLM, VINT
      integer :: I, IBEG, IEND, II, III, ISTAGE, ISTEP, ISTEP2, IT, J, L
      integer :: LM1, N, NNLEFT, N2

      CADRE = 0.0_rk
      ERROR = 0.0_rk
      IFLAG = 1
      LENGTH = ABS (B - A)
      IF (LENGTH == 0.0_rk) RETURN
      ERRR = min (0.1_rk, max (ABS (RERR), 10.0_rk * TOLMCH) )
      ERRA = ABS (AERR)
      STEPMN = max (LENGTH / 2.0_rk**MXSTGE,                            &
     &              max (LENGTH, ABS (A),  ABS (B) ) * TOLMCH)
      STAGE = 0.5_rk
      ISTAGE = 1
      CUREST = 0.0_rk
      FNSIZE = 0.0_rk
      PREVER = 0.0_rk
      REGLAR = .FALSE.
!
!  THE GIVEN INTERVAL OF INTEGRATION IS THE FIRST INTERVAL CONSIDERED.
      BEG = A
      FBEG = 0.5_rk * F (BEG)
      TS (1) = FBEG
      IBEG = 1
      END = B
      FEND = 0.5_rk * F (END)
      TS (2) = FEND
      IEND = 2
!
    5 RIGHT = .FALSE.
!
!  INVESTIGATION OF A PARTICULAR SUBINTERVAL BEGINS AT THIS POINT.
!           *****   MAJOR VARIABLES   *****
!  BEG,
!  END     ENDPOINTS OF THE CURRENT INTERVAL
!  FBEG,
!  FEND    ONE HALF THE VALUE OF F(X) AT THE ENDPOINTS
!  STEP    SIGNED LENGTH OF CURRENT SUBINTERVAL
!  ISTAGE  HEIGHT OF CURRENT SUBINTERVAL IN STACK OF SUBINTERVALS
!          YET TO BE DONE
!  RIGHT   A LOGICAL VARIABLE INDICATING WHETHER CURRENT SUBINTERVAL
!          IS RIGHT HALF OF PREVIOUS SUBINTERVAL. NEEDED IN 80FF AND
!          90FF TO DECIDE WHAT INTERVAL TO LOOK AT NEXT.
!  TS(I), I=IBEG,...,IEND, CONTAINS THE FUNCTION VALUES FOR THIS
!          SUBINTERVAL SO FAR COMPUTED. SPECIFICALLY,
!            TS(I) = F(BEG + (I-IBEG)/(IEND-IBEG)*STEP), ALL I
!          EXCEPT THAT TS(IBEG) = FBEG, TS(IEND) = FEND
!  REGLAR  A LOGICAL VARIABLE INDICATING WHETHER OR NOT THE CURRENT
!          SUBINTERVAL IS REGULAR (SEE NOTES)
!  H2CONV  A LOGICAL VARIABLE INDICATING WHETHER H**2 CONVERGENCE OF
!          THE TRAPEZOID SUMS FOR THIS INTERVAL IS RECOGNIZED,
!  AITKEN  A LOGICAL VARIABLE INDICATING WHETHER CONVERGENCE OF RATIOS
!          FOR THIS SUBINTERVAL IS RECOGNIZED
!  T       CONTAINS THE FIRST *L* ROWS OF THE ROMBERG T-TABLE FOR THIS
!          SUBINTERVAL IN ITS LOWER TRIANGULAR PART. SPECIFICALLY,
!            T(I,1) = TRAPEZOID SUM (WITHOUT THE FACTOR *STEP*)
!                     ON 2**(I-1) + 1 EQUISPACED POINTS, I=1,...,L,
!            T(I,J+1) = T(I,J) + (T(I,J)-T(I-1,J))/(4**J - 1),
!                       J=2,...,I-1, I=2,...,L.
!          FURTHER, THE STRICTLY UPPER TRIANGULAR PART OF T CONTAINS
!          THE RATIOS FOR THE VARIOUS COLUMNS OF THE T-TABLE.
!          SPECIFICALLY,
!            T(J,I) = (T(I,J)-T(I-1,J))/(T(I+1,J)-T(I,J)),
!                     I=J+1,...,L-1,  J=1,...,L-2.
!          FINALLY, THE LAST OR L-TH COLUMN CONTAINS
!            T(J,L) = T(L,J) - T(L-1,J), J=1,...,L-1.
    6 STEP = END - BEG
      ASTEP = ABS (STEP)
      IF (ASTEP < STEPMN) GOTO 950
      IF (LEVEL >= 3) write (*, 609) BEG, STEP, ISTAGE
  609 FORMAT(" BEG,STEP ",2E16.8,I5)
      T (1, 1) = FBEG + FEND
      TABS = ABS (FBEG) + ABS (FEND)
      L = 1
      N = 1
      H2CONV = .FALSE.
      AITKEN = .FALSE.
      GOTO 10
!
    9 IF (LEVEL >= 4) write (*, 692) L, T (1, LM1)
   10 LM1 = L
      L = L + 1
!
!CALCULATE THE NEXT TRAPEZOID SUM, T(L,1), WHICH IS BASED ON
!  *N2*  + 1 EQUISPACED POINTS. HERE,  N2 = N*2 = 2**(L-1) .
      N2 = N * 2
      FN = N2
      ISTEP = (IEND-IBEG) / N
      IF (ISTEP <= 1) then
        II = IEND
        IEND = IEND+N
        IF (IEND > MAXTS) GOTO 900
        HOVN = STEP / FN
        III = IEND
        DO I = 1, N2, 2
           TS (III) = TS (II)
           TS (III - 1) = F (END-FLOAT (I) * HOVN)
           III = III - 2
           II = II - 1
        end do
        ISTEP = 2
      end if
      ISTEP2 = IBEG + ISTEP / 2
      SUM = 0.0_rk
      SUMABS = 0.0_rk
      DO I = ISTEP2, IEND, ISTEP
         SUM = SUM + TS (I)
         SUMABS = SUMABS + ABS (TS (I) )
      end do
      T (L, 1) = 0.5_rk * T (L - 1, 1) + SUM / FN
      TABS = 0.5_rk * TABS + SUMABS / FN
      ABSI = ASTEP * TABS
      N = N2
!
!  GET PRELIMINARY VALUE FOR *VINT* FROM LAST TRAPEZOID SUM AND
!  UPDATE THE ERROR REQUIREMENT *ERGOAL* FOR THIS SUBINTERVAL.
!  THE ERROR REQUIREMENT IS NOT PRORATED ACCORDING TO THE LENGTH OF
!  THE CURRENT SUBINTERVAL RELATIVE TO THE INTERVAL OF INTEGRATION,
!  BUT ACCORDING TO THE HEIGHT *ISTAGE* OF THE CURRENT SUBINTERVAL
!  IN THE STACK OF SUBINTERVALS YET TO BE DONE.
!  THIS PROCEDURE IS NOT BACKED BY ANY RIGOROUS ARGUMENT, BUT
!  SEEMS TO WORK.
      IT = 1
      VINT = STEP * T (L, 1)
      TABTLM = TABS * TOLMCH
      FNSIZE = max (FNSIZE, ABS (T (L, 1) ) )
      ERGOAL = max (ASTEP * TOLMCH * FNSIZE,                            &
     &              STAGE * max (ERRA,  ERRR * ABS (CUREST + VINT) ) )
!
!COMPLETE ROW L AND COLUMN L OF *T* ARRAY.
      FEXTRP = 1.0_rk
      DO I = 1, LM1
         FEXTRP = FEXTRP * 4.0_rk
         T (I, L) = T (L, I) - T (L - 1, I)
         T (L, I + 1) = T (L, I) + T (I, L) / (FEXTRP - 1.)
      end do
      ERRER = ASTEP * ABS (T (1, L) )
!---------------------------------------------------------------------
!---  PRELIMINARY DECISION PROCEDURE  --------------------------------
!
!  IF L = 2 AND T(2,1) = T(1,1), GO TO 60 TO FOLLOW UP THE IMPRESSION
!  THAT INTEGRAND IS STRAIGHT LINE.
      IF ( L <= 2 ) then
         IF (ABS (T (1, 2) ) <= TABTLM) GOTO 60
         GOTO 10
      end if
!
!CALCULATE NEXT RATIOS FOR COLUMNS 1,...,L-2 OF T-TABLE
!  RATIO IS SET TO ZERO IF DIFFERENCE IN LAST TWO ENTRIES OF
!  COLUMN IS ABOUT ZERO.
   15 DO I = 2, LM1
         DIFF = 0.0_rk
         IF (ABS (T (I - 1, L) ) > TABTLM) DIFF = T (I - 1, LM1)        &
     &                                          / T (I - 1, L)
         T (I - 1, LM1) = DIFF
      end do
!
!  T(1,LM1) IS THE RATIO DERIVED FROM LAST THREE TRAPEZOID SUMS, I.E.,
!     T(1,LM1) = (T(L-1,1)-T(L-2,1))/(T(L,1)-T(L-1,1)) .
!  IF THIS RATIO IS ABOUT 4, DO ROMBERG EXTRAPOLATION.
!  IF THIS RATIO IS ZERO, I.E., IF LAST TVO TRAPEZOID SUMS ARE ABOUT
!     EQUAL, MAYBE DO NOISE CHECK.
!  IF THIS RATIO IS ABOUT 2 IN ABSOLUTE VALUE, GO TO 50 WITH THE
!     BELIEF THAT INTEGRAND HAS JUMP DISCONTINUITY.
!  IF THIS RATIO IS, AT LEAST, ABOUT EQUAL TO THE PREVIOUS RATIO, THEN
!     THE INTEGRAND MAY WELL HAVE A NICE INTEGRABLE SINGULARITY.
!     GO TO 30 TO FOLLOW UP THIS HUNCH.
      IF (ABS (4.0_rk - T (1, LM1) ) > H2TOL) then
         IF (T (1, LM1) /= 0.0_rk) then
            IF (ABS (2.0_rk - ABS (T (1, LM1) ) ) < JUMPTL) GOTO 50
            IF (L == 3) GOTO 9
            H2CONV = .FALSE.
            IF (ABS((T(1, LM1) - T(1, L - 2) ) / T(1, LM1) ) <= AITTOL) &
     &      GOTO 30
!  AT THIS POINT, NO REGULAR BEHAVIOUR WAS DETECTED.
!  IF CURRENT SUBINTERVAL IS NOT REGULAR AND ONLY FOUR TRAPEZOID SUMS
!  WERE COMPUTED SO FAR, TRY ONE MORE TRAPEZOID SUM.
!  IF, AT LEAST, LAST TWO TRAPEZOID SUMS ARE ABOUT EQUAL, THEN
!  FAILURE TO RECOGNIZE REGULAR BEHAVIOUR MAY WELL BE DUE TO NOISE.
!  GO TO 70 TO CHECK THIS OUT.
!  OTHERWISE, GO TO 90 FOR FURTHER SUBDIVISION.
!
            IF ( .not. REGLAR ) then
               IF (L == 4) GOTO 9
            end if
         end if
         IF (ERRER <= ERGOAL) GOTO 70
         IF (LEVEL >= 4) write (*, 692) L, T (1, LM1)
         GOTO 91
      end if
!----------------------------------------------------------------------
!CAUTIOUS ROMBERG EXTRAPOLATION  --------------------------------------
!
!  THE CURRENT, OR L-TH, ROW OF THE ROMBERG T-TABLE HAS L ENTRIES.
!  FOR J=1,...,L-2, THE ESTIMATE
!             STEP*T(L,J+1)
!  IS BELIEVED TO HAVE ITS ERROR BOUNDED BY
!    ABS(STEP*(T(L,J)-T(L-1,J))/(4**J - 1))
!  IF THE LAST RATIO
!      T(J,LM1) = (T(L-1,J)-T(L-2,J))/(T(L,J)-T(L-1,J))
!  FOR COLUMN J OF THE T-TABLE IS ABOUT 4**J.
!  THE FOLLOWING IS A SLIGHTLY RELAXED EXPRESSION OF THIS BELIEF.
      IF (LEVEL >= 4) write (*, 619) L, T (1, LM1)
  619 FORMAT(I5,E16.8,5X,"H2CONV")
      IF ( .not. H2CONV ) then
         AITKEN = .FALSE.
         H2CONV = .TRUE.
         IF (LEVEL >= 3) write (*, 620) L
  620    FORMAT(" H2 CONVERGENCE AT ROW",I3)
      end if
   21 FEXTRP = 4.0_rk
      do
         IT = IT + 1
         VINT = STEP * T (L, IT)
         ERRER = ABS (STEP / (FEXTRP - 1.0_rk) * T (IT - 1, L) )
         IF (ERRER <= ERGOAL) GOTO 80
         IF (IT == LM1) GOTO 40
         IF (T (IT, LM1) /= 0.0_rk) then
            IF (T (IT, LM1) <= FEXTRP) GOTO 40
            IF (ABS (T (IT, LM1) / 4. - FEXTRP) / FEXTRP < AITTOL)      &
     &        FEXTRP = FEXTRP * 4.0_rk
         end if
      end do
!---------------------------------------------------------------------
!--- INTEGRAND MAY HAVE X**ALPHA TYPE SINGULARITY---------------------
!  RESULTING IN A RATIO OF *SING*  =  2**(ALPHA + 1)
   30 IF (LEVEL >= 4) write (*, 629) L, T (1, LM1)
  629 FORMAT(I5,E16.8,5X,"AITKEN")
      IF (T (1, LM1) < AITLOW) GOTO 91
      IF ( .not. AITKEN ) then
         H2CONV = .FALSE.
         AITKEN = .TRUE.
         IF (LEVEL >= 3) write (*, 630) L
  630    FORMAT(" AITKEN AT ROW",I3)
      end if
      FEXTRP = T (L - 2, LM1)
      IF (FEXTRP > 4.5_rk) GOTO 21
      IF (FEXTRP < AITLOW) GOTO 91
      IF (ABS (FEXTRP - T (L - 3, LM1) ) / T (1, LM1) > H2TOL) GOTO 91
      IF (LEVEL >= 3) write (*, 631) FEXTRP
  631 FORMAT(" RATIO",F12.8)
      SING = FEXTRP
      FEXTM1 = FEXTRP - 1.0_rk
      DO I = 2, L
         AIT (I) = T (I, 1) + (T (I, 1) - T (I - 1, 1) ) / FEXTM1
         R (I) = T (1, I - 1)
         DIF (I) = AIT (I) - AIT (I - 1)
      end do
      IT = 2
      do
         VINT = STEP * AIT (L)
         IF (LEVEL >= 5) then
            write (*, 632) (R (I + 1), I = IT, LM1)
            write (*, 632) (AIT (I), I = IT, L)
            write (*, 632) (DIF (I + 1), I = IT, LM1)
  632       FORMAT(1X,8E15.8)
         end if
         ERRER = ERRER / FEXTM1
         IF ( ERRER <= ERGOAL ) then
            ALPHA = LOG10 (SING) / ALG4O2 - 1.0_rk
            IF (LEVEL >= 2) write (*, 633) ALPHA, BEG, END
  633       FORMAT(12X,"INTEGRAND SHOWS SINGULAR BEHAVIOR OF TYPE",     &
     &             " X**(",F4.2,") BETWEEN",E15.8," AND",E15.8)
            IFLAG = MAX (IFLAG, 2)
            GOTO 80
         end if
         IT = IT + 1
         IF (IT == LM1) exit
         IF ( IT <= 3 ) then
            H2NEXT = 4.0_rk
            SINGNX = 2.0_rk * SING
         end if
         IF ( H2NEXT >= SINGNX ) then
           FEXTRP = SINGNX
           SINGNX = 2.0_rk * SINGNX
         else
           FEXTRP = H2NEXT
           H2NEXT = 4.0_rk * H2NEXT
         end if
         DO I = IT, LM1
            R (I + 1) = 0.0_rk
            IF (ABS (DIF (I + 1) ) > TABTLM) R (I + 1) = DIF (I) / DIF  &
     &         (I + 1)
         end do
         IF (LEVEL >= 4) write (*, 638) FEXTRP, R (L - 1), R (L)
  638    FORMAT(" FEXTRP + RATIOS",3E15.8)
         H2TFEX = - H2TOL * FEXTRP
         IF (R (L) - FEXTRP < H2TFEX) exit
         IF (R (L - 1) - FEXTRP < H2TFEX) exit
         ERRER = ASTEP * ABS (DIF (L) )
         FEXTM1 = FEXTRP - 1.0_rk
         DO I = IT, L
            AIT (I) = AIT (I) + DIF (I) / FEXTM1
            DIF (I) = AIT (I) - AIT (I - 1)
         end do
      end do
!CURRENT TRAPEZOID SUM AND RESULTING EXTRAPOLATED VALUES DID NOT GIVE
!  A SMALL ENOUGH *ERRER*.
!  IF LESS THAN FIVE TRAPEZOID SUMS WERE COMPUTED SO FAR, TRY NEXT
!     TRAPEZOID SUM.
!  OTHERWISE, DECIDE WHETHER TO GO ON OR TO SUBDIVIDE AS FOLLOWS.
!  WITH T(L,IT) GIVING THE CURRENTLY BEST ESTIMATE, GIVE UP ON DEVELOP
!  ING THE T-TABLE FURTHER IF  L > IT+2, I.E., IF EXTRAPOLATION
!  DID NOT GO VERY FAR INTO THE T-TABLE.
!  FURTHER, GIVE UP IF REDUCTION IN *ERRER* AT THE CURRENT RATE
!  DOES NOT PREDICT AN *ERRER* LESS THAN *ERGOAL* BY THE TIME
!  *MAXTBL* TRAPEZOID SUMS HAVE BEEN COMPUTED.
!  ---NOTE---
!  HAVING  PREVER < ERRER  IS AN ALMOST CERTAIN SIGN OF BEGINNING
!  TROUBLE WITH NOISE IN THE FUNCTION VALUES. HENCE,
!  A WATCH FOR, AND CONTROL OF, NOISE SHOULD BEGIN HERE.
   40 FEXTRP = max (PREVER / ERRER, AITLOW)
      PREVER = ERRER
      IF (L < 5) GOTO 10
      IF (LEVEL >= 3) write (*, 641) ERRER, ERGOAL, FEXTRP, IT
  641 FORMAT(" ERRER,ERGOAL,FEXTRP,IT",2E15.8,E14.5,I3)
      IF (L - IT > 2.AND.ISTAGE < MXSTGE) GOTO 90
      IF (ERRER / FEXTRP** (MAXTBL - L) < ERGOAL) GOTO 10
      GOTO 90
!---------------------------------------------------------------------
!----    INTEGRAND HAS JUMP (SEE NOTES)  -----------------------------
   50 IF (LEVEL >= 4) write (*, 649) L, T (1, LM1)
  649 FORMAT(I5,E16.8,5X,"JUMP")
      IF (ERRER > ERGOAL) GOTO 90
!  NOTE THAT  2*FN = 2**L
      DIFF = ABS (T (1, L) ) * 2. * FN
      IF (LEVEL >= 2) write (*, 650) DIFF, BEG, END
  650 FORMAT(12X,"INTEGRAND SEEMS TO HAVE JUMP OF SIZE",E13.6,          &
     &       " BETWEEN",E15.8," AND",E15.8)
      GOTO 80
!---------------------------------------------------------------------
!----  INTEGRAND IS STRAIGHT LINE  -----------------------------------
!  TEST THIS ASSUMPTION BY COMPARING THE VALUE OF THE INTEGRAND AT
!  FOUR 'RANDOMLY CHOSEN' POINTS WITH THE VALUE OF THE STRAIGHT LINE
!  INTERPOLATING THE INTEGRAND AT THE TWO ENDPOINTS OF THE SUB-
!  INTERVAL. IF TEST IS PASSED, ACCEPT *VINT* .
   60 IF (LEVEL >= 4) write (*, 660) L
  660 FORMAT(I5,21X,"STRAIGHT LINE")
      SLOPE = (FEND-FBEG) * 2.0_rk
      FBEG2 = FBEG * 2.0_rk
      DO I = 1, 4
         DIFF = ABS (F (BEG + RN (I) * STEP) - FBEG2 - RN (I) * SLOPE)
         IF (DIFF > TABTLM) GOTO 72
      end do
      IF (LEVEL >= 3) write (*, 667) BEG, END
  667 FORMAT(27X,"INTEGRAND SEEMS TO BE STRAIGHT LINE BETWEEN",         &
     &       E15.8," AND",E15.8)
      GOTO 80
!---------------------------------------------------------------------
!------  NOISE MAY BE DOMINANT FEATURE  ------------------------------
!  ESTIMATE NOISE LEVEL BY COMPARING THE VALUE OF THE INTEGRAND AT
!  FOUR 'RANDOMLY CHOSEN' POINTS WITH THE VALUE OF THE STRAIGHT LINE
!  INTERPOLATING THE INTEGRAND AT THE TWO ENDPOINTS. IF SMALL
!  ENOUGH, ACCEPT *VINT*.
   70 IF (LEVEL >= 4) write (*, 670) L, T (1, LM1)
  670 FORMAT(I5,E16.8,5X,"NOISE")
      SLOPE = (FEND-FBEG) * 2.0_rk
      FBEG2 = FBEG * 2.0_rk
      I = 1
   71 DIFF = ABS (F (BEG + RN (I) * STEP) - FBEG2 - RN (I) * SLOPE)
   72 ERRER = max (ERRER, ASTEP * DIFF)
      IF (ERRER > ERGOAL) GOTO 91
      I = I + 1
      IF (I <= 4) GOTO 71
      IF (LEVEL >= 3) write (*, 671) BEG, END
  671 FORMAT(" NOISE BETWEEN ",E15.8," AND",E15.8)
      IFLAG = 3
!---------------------------------------------------------------------
!--- INTEGRATION OVER CURRENT SUBINTERVAL SUCCESSFUL -----------------
!  ADD *VINT* TO *CADRE* AND *ERRER* TO *ERROR*, THEN SET UP NEXT
!  SUBINTERVAL, IF ANY.
   80 CADRE = CADRE+VINT
      ERROR = ERROR + ERRER
      IF (LEVEL >= 3) then
         IF (LEVEL >= 5) then
            DO I = 1, L
              write (*, 692) I, (T (I, J), J = 1, L)
            end do
         end if
         write (*, 682) VINT, ERRER, L, IT
  682    FORMAT(" INTEGRAL IS",E16.8,", ERROR",E15.8,"  FROM T(",       &
     &          I0,",",I0,")")
      end if
!
      IF ( .not. RIGHT ) then
         ISTAGE = ISTAGE-1
         IF (ISTAGE == 0) RETURN
!
         REGLAR = REGLSV (ISTAGE)
         BEG = BEGIN (ISTAGE)
         END = FINIS (ISTAGE)
         CUREST = CUREST - EST (ISTAGE+1) + VINT
         IEND = IBEG - 1
         FEND = TS (IEND)
         IBEG = IBEGS (ISTAGE)
         GOTO 94
      end if
      CUREST = CUREST + VINT
      STAGE = STAGE * 2.0_rk
      IEND = IBEG
      IBEG = IBEGS (ISTAGE)
      END = BEG
      BEG = BEGIN (ISTAGE)
      FEND = FBEG
      FBEG = TS (IBEG)
      GOTO 5
!---------------------------------------------------------------------
!--- INTEGRATION OVER CURRENT SUBINTERVAL IS UNSUCCESSFUL-------------
!  MARK SUBINTERVAL FOR FURTHER SUBDIVISION. SET UP NEXT SUBINTERVAL.
   90 REGLAR = .TRUE.
   91 IF (ISTAGE == MXSTGE) GOTO 950
      IF (LEVEL >= 5) then
         DO I = 1, L
            write (*, 692) I, (T (I, J), J = 1, L)
         end do
  692    FORMAT(I5,7E16.8/3E16.8)
      end if
      IF (RIGHT) GOTO 95
      REGLSV (ISTAGE+1) = REGLAR
      BEGIN (ISTAGE) = BEG
      IBEGS (ISTAGE) = IBEG
      STAGE = 0.5_rk * STAGE
   94 RIGHT = .TRUE.
      BEG = 0.5_rk * (BEG + END)
      IBEG = (IBEG + IEND) / 2
      TS (IBEG) = 0.5_rk * TS (IBEG)
      FBEG = TS (IBEG)
      GOTO 6
   95 NNLEFT = IBEG - IBEGS (ISTAGE)
      IF ( IEND+NNLEFT < MAXTS ) then
         III = IBEGS (ISTAGE)
         II = IEND
         DO I = III, IBEG
            II = II + 1
            TS (II) = TS (I)
         end do
         DO I = IBEG, II
            TS (III) = TS (I)
            III = III + 1
         end do
         IEND = IEND+1
         IBEG = IEND-NNLEFT
         FEND = FBEG
         FBEG = TS (IBEG)
         FINIS (ISTAGE) = END
         END = BEG
         BEG = BEGIN (ISTAGE)
         BEGIN (ISTAGE) = END
         REGLSV (ISTAGE) = REGLAR
         ISTAGE = ISTAGE+1
         REGLAR = REGLSV (ISTAGE)
         EST (ISTAGE) = VINT
         CUREST = CUREST + EST (ISTAGE)
         GOTO 5
      end if
!---------------------------------------------------------------------
!--- FAILURE TO HANDLE GIVEN INTEGRATION PROBLEM----------------------
  900 IF (LEVEL >= 2) write (*, 6900) BEG, END
 6900 FORMAT(" TOO MANY FUNCTION EVALUATIONS AROUND"/                   &
     &       10X,E15.8," AND",E15.8)
      IFLAG = 4
      GOTO 999
  950 IFLAG = 5
      IF ( LEVEL >= 2 ) then
         IF ( LEVEL >= 5 ) then
            DO I = 1, L
               write (*, 692) I, (T (I, J), J = 1, L)
            end do
         end if
         write (*, 6959) BEG, END
 6959    FORMAT(12X,"INTEGRAND SHOWS SINGULAR BEHAVIOUR OF "            &
     &          ,"UNKNOWN TYPE BETWEEN",E15.8," AND",E15.8)
      end if
  999 CADRE = CUREST + VINT
      RETURN
      END FUNCTION CADRE

      subroutine CADRE_Reverse (Answer, A, B, AERR, RERR, LEVEL, ERROR, &
     &                          IFLAG)
!  THIS Subroutine RETURNS AN ESTIMATE *Answer* FOR THE NUMBER
!     INT = INTEGRAL OF *F*(X) FROM *A* TO *B*
!  WHICH HOPEFULLY SATISFIES
!     ABS(INT - *Answer**) <= max(*AERR*, *RERR* TIMES ABS(INT)).
!     THE PROGRAM USES CAUTIOUS ADAPTIVE ROMBERG EXTRAPOLATION.
!  IN THIS SCHEME, THE INTEGRAL IS CALCULATED AS THE SUM OF INTEGRALS
!  OVER SUITABLY SMALL SUBINTERVALS. ON EACH SUBINTERVAL, AN ESTIMATE
!  *VINT*, WITH ESTIMATED ABSOLUTE ERROR *ERRER*, IS FOUND BY CAUTIOUS
!  ROMBERG EXTRAPOLATION. IF *ERRER* IS SMALL ENOUGH, *VINT* IS ACCEPT
!  ED AND ADDED TO *Answer*, AND *ERRER* IS ADDED TO *ERROR*. OTHERWISE
!  THE SUBINTERVAL IS HALVED, AND EACH HALF IS CONSIDERED SEPARATELY,
!  INFORMATION ABOUT THE OTHER HALF BEING TEMPORARILY STACKED.
!
!  USAGE:
!    iflag = 0
!    do
!      call cadre_reverse ( Answer, A, B, Aerr, Rerr, Level, Error, Iflag )
!      if ( iflag > 0 ) exit
!      ! Evaluate integrand at Answer and put the result in Answer
!    end do

        real(rk), intent(inout) :: Answer
        real(rk), intent(in) :: A, B
        real(rk), intent(in) :: Aerr, Rerr
        integer, intent(in) :: Level
        real(rk), intent(inout) :: Error ! inout so it's not undefined
                                         ! at every re-entry
        integer, intent(inout) :: Iflag

!                   *****   INPUT   *****
!
!  Answer  Integrand value
!  A,B     THE TWO ENDPOINTS OF THE INTERVAL OF INTEGRATION
!  AERR
!  RERR    DESIRED ABSOLUTE AND RELATIVE ERROR IN THE ANSWER
!  LEVEL   AN INTEGER INDICATING DESIRED LEVEL OF PRINTOUT
! <= 1,  NO PRINTOUT,
!           =  2,  FAILURE MESSAGE (IF ANY), AND LIST OF SINGULAR-
!                  ITIES ENCOUNTERED (IF ANY),
!           =  3,  IN ADDITION, ALL SUBINTERVALS CONSIDERED ARE LISTED
!                  TOGETHER WITH THE KIND OF REGULAR BEHAVIOUR FOUND
!                  (IF ANY),
!           =  4,  IN ADDITION, ALL RATIOS CONSIDERED ARE LISTED AS IS
!                  INFO ON WHICH DECISION PROCEDURE IS BASED,
! >= 5,  IN ADDITION, ALL T-TABLES ARE LISTED.
!  IFLAG  must initially be zero.  Do not change it otherwise.
!
!                   *****   OUTPUT  *****
!
!  Answer  If IFLAG <= 0, abscissa at which to evaluate integrand
!          If IFLAG > 0,  ESTIMATE OF THE INTEGRAL
!  ERROR   ESTIMATED BOUND ON THE ABSOLUTE ERROR OF THE NUMBER *Answer*
!  IFLAG   AN INTEGER BETWEEN 1 AND 5 INDICATING WHAT DIFFICULTIES
!          WERE MET WITH, SPECIFICALLY
!          < 0, Evaluate integrand at Answer, put the value in Answer
!          = 1, ALL IS WELL,
!          = 2, ONE OR MORE SINGULARITIES WERE SUCCESSFULLY HANDLED,
!          = 3, IN SOME SUBINTERVAL(S), THE ESTIMATE *VINT* WAS ACCEPT
!               ED MERELY BECAUSE *ERRER* WAS SMALL, EVEN THOUGH NO
!               REGULAR BEHAVIOUR COULD BE RECOGNIZED,
!          = 4, FAILURE, OVERFLOW OF STACK *TS* (THIS HAS NEVER HAPPEN
!               ED,  - SO FAR),
!          = 5, FAILURE, TOO SMALL A SUBINTERVAL IS REQUIRED. THIS MAY
!               BE DUE TO TOO MUCH NOISE IN THE FUNCTION (RELATIVE TO
!               THE GIVEN ERROR REQUIREMENTS) OR DUE TO A PLAIN ORNERY
!               INTEGRAND.
!  A VERY CAUTIOUS MAN WOULD ACCEPT *Answer* ONLY IF IFLAG IS 1 OR 2.
!  THE MERELY REASONABLE MAN WOULD KEEP THE FAITH EVEN IF IFLAG IS 3.
!  THE ADVENTUROUS MAN IS QUITE OFTEN RIGHT IN ACCEPTING *Answer*
!  EVEN IF IFLAG IS 4 OR 5.
!
!     *****   LIST OF MAJOR VARIABLES   *****
!
!  CUREST  BEST ESTIMATE SO FAR FOR
!             INT - (INTEGRAL OVER CURRENTLY CONSIDERED SUBINTERVAL).
!  FNSIZE  MAXIMUM AVERAGE FUNCTION SIZE SO FAR ENCOUNTERED,
!  ERRR    RELATIVE ERROR REQUIREMENT USED. DERIVED FROM INPUT *RERR*
!          AND CHOSEN TO LIE BETWEEN .1 AND 10 TIMES *TOLMCH*,
!  ERRA    = ABS(*AERR*)
!  STAGE   (MORE OR LESS) EQUAL TO 2 TO THE -(*ISTAGE*)
!  THESE FIVE QUANTITIES ARE USED IN THE DETERMINATION OF THE LOCAL
!  ERROR REQUIREMENT.
!
!  STEPMN  MINIMUM SUBINTERVAL LENGTH PERMITTED,
!  TS      STACK OF VALUES OF F(X) SO FAR COMPUTED BUT NOT YET
!          SUCCESSFULLY USED,
!  ISTAGE  AN INTEGER INDICATING THE HEIGHT OF THE STACK OF INTERVALS
!          YET TO BE PROCESSED.
!
!     *****   LIST OF PARAMETERS   *****
!
!  TOLMCH  DEPENDS ON THE LENGTH OF FLOATING POINT MANTISSA. SHOULD BE
!          ABOUT 1.E-7 FOR 27 BINARY BIT MANTISSA AND
!          ABOUT 1.E-13 FOR 48 BINARY BIT MANTISSA.
!  AITLOW  SHOULD BE SOMEWHAT GREATER THAN 1.
!  H2TOL,
!  AITTOL,
!  JUMPTL  TOLERANCES USED IN THE DECISION PROCESS TO RECOGNIZE
!          H**2 CONVERGENCE, X**ALPHA TYPE CONVERGENCE, OR
!          JUMP-TYPE CONVERGENCE OF THE TRAPEZOID SUMS.
!  MAXTS,
!  MAXTBL,
!  MXSTGE  ARE THE THREE DIFFERENT UPPER LIMITS FOR THE DIMENSION OF
!          THE VARIOUS ARRAYS.
!
!        *****   PROGRAM LAYOUT   *****
!
!          INITIALIZATION
!  5,6     BEGIN WORK ON NEXT SUBINTERVAL
!   9-14   GET NEXT TRAPEZOID SUM
!  15-19   GET RATIOS. PRELIMINARY DECISION PROCEDURE.
!  20-     ESTIMATE *VINT* ASSUMING SMOOTH INTEGRAND
!  30-     ESTIMATE *VINT* ASSUMING INTEGRAND HAS X**ALPHA TYPE
!          SINGULARITY
!  40-     NO LUCK WITH THIS TRAPEZOID SUM. GET NEXT ONE OR GET OUT.
!  50-     ESTIMATE *VINT* ASSUMING INTEGRAND HAS JUMP
!  60-     ESTIMATE *VINT* ASSUMING INTEGRAND IS STRAIGHT LINE
!  70-     ESTIMATE *VINT* ASSUMING VARIATION IN INTEGRAND
!          IS MOSTLY NOISE.
!  80-     INTEGRATION OVER CURRENT SUBINTERVAL SUCCESSFUL.
!          SET UP NEXT SUBINTERVAL, IF ANY, OR RETURN.
!  90-     INTEGRATION OVER CURRENT SUBINTERVAL NOT SUCCESSFUL.
!          MARK CURRENT SUBINTERVAL FOR SUBDIVISION AND SET UP
!          NEXT SUBINTERVAL.
!  900-    FAILURE.
!
      integer, parameter :: MAXTS = 2049
      integer, parameter :: MAXTBL = 10
      integer, parameter :: MXSTGE = 30
      real(rk) :: T (MAXTBL, MAXTBL), R (MAXTBL), AIT (MAXTBL),         &
     & DIF (MAXTBL), TS (MAXTS), IBEGS (MXSTGE), BEGIN (MXSTGE),        &
     & FINIS (MXSTGE), EST (MXSTGE)
      real(rk) :: LENGTH
      LOGICAL H2CONV, AITKEN, RIGHT, REGLAR, REGLSV (MXSTGE)
      real(rk), parameter :: JUMPTL = 0.01_rk
      real(rk), parameter :: TOLMCH = epsilon(1.0_rk)
      real(rk), parameter :: AITLOW = 1.1_rk
      real(rk), parameter :: H2TOL = 0.15_rk, AITTOL = 0.1_rk
!     "Randomly chosen" points to check for a straight line
      real(rk), parameter :: RN (4) =                                   &
     & (/ 0.71420053_rk, 0.34662815_rk, 0.843751_rk, 0.12633046_rk /)
      real(rk), parameter :: ALG4O2 = 0.3010299956639795_rk ! log10(2)
      DATA AIT (1) / 0._rk /
      real(rk) :: ABSI, ALPHA, ASTEP, BEG, CUREST, DIFF, END, ERGOAL
      real(rk) :: ERRA, ERRER, ERRR, FBEG, FBEG2, FEND, FEXTM1, FEXTRP
      real(rk) :: FN, FNSIZE, HOVN, H2NEXT, H2TFEX, PREVER, SING
      real(rk) :: SINGNX, SLOPE, STAGE, STEP, STEPMN, SUM, SUMABS
      real(rk) :: TABS, TABTLM, VINT
      real(rk) :: RESULT
      integer :: Final_Iflag ! Status detected during integration
      integer :: I, IBEG, IEND, II, III, ISTAGE, ISTEP, ISTEP2, IT, J, L
      integer :: LM1, N, NNLEFT, N2
      save

      if ( iflag > 0 ) then
        write ( *, 1000 ) iflag
1000    format ( " Improper usage: Cadre_Reverse entered with IFLAG = ",&
     &           i0, " > 0" )
        answer = -huge(0.0_rk)
        error = -huge(0.0_rk)
        return
      end if
      if ( iflag < 0 ) go to ( 1001, 1002, 1003, 1004, 1005 ), -iflag

      answer = 0.0_rk
      result = 0.0_rk
      ERROR = 0.0_rk
      FINAL_IFLAG = 1
      IFLAG = 1
      LENGTH = ABS (B - A)
      IF (LENGTH == 0.0_rk) RETURN
      ERRR = min (0.1_rk, max (ABS (RERR), 10.0_rk * TOLMCH) )
      ERRA = ABS (AERR)
      STEPMN = max (LENGTH / 2.0_rk**MXSTGE,                            &
     &              max (LENGTH, ABS (A), ABS (B) ) * TOLMCH)
      STAGE = 0.5_rk
      ISTAGE = 1
      CUREST = 0.0_rk
      FNSIZE = 0.0_rk
      PREVER = 0.0_rk
      REGLAR = .FALSE.
!
!  THE GIVEN INTERVAL OF INTEGRATION IS THE FIRST INTERVAL CONSIDERED.
      BEG = A
      answer = beg
      iflag = -1
      return
 1001 fbeg = 0.5_rk * answer
      TS (1) = FBEG
      IBEG = 1
      END = B
      answer = end
      iflag = -2
      return
 1002 fend = 0.5_rk * answer
      TS (2) = FEND
      IEND = 2
!
    5 RIGHT = .FALSE.
!
!  INVESTIGATION OF A PARTICULAR SUBINTERVAL BEGINS AT THIS POINT.
!           *****   MAJOR VARIABLES   *****
!  BEG,
!  END     ENDPOINTS OF THE CURRENT INTERVAL
!  FBEG,
!  FEND    ONE HALF THE VALUE OF F(X) AT THE ENDPOINTS
!  STEP    SIGNED LENGTH OF CURRENT SUBINTERVAL
!  ISTAGE  HEIGHT OF CURRENT SUBINTERVAL IN STACK OF SUBINTERVALS
!          YET TO BE DONE
!  RIGHT   A LOGICAL VARIABLE INDICATING WHETHER CURRENT SUBINTERVAL
!          IS RIGHT HALF OF PREVIOUS SUBINTERVAL. NEEDED IN 80FF AND
!          90FF TO DECIDE WHAT INTERVAL TO LOOK AT NEXT.
!  TS(I), I=IBEG,...,IEND, CONTAINS THE FUNCTION VALUES FOR THIS
!          SUBINTERVAL SO FAR COMPUTED. SPECIFICALLY,
!            TS(I) = F(BEG + (I-IBEG)/(IEND-IBEG)*STEP), ALL I
!          EXCEPT THAT TS(IBEG) = FBEG, TS(IEND) = FEND
!  REGLAR  A LOGICAL VARIABLE INDICATING WHETHER OR NOT THE CURRENT
!          SUBINTERVAL IS REGULAR (SEE NOTES)
!  H2CONV  A LOGICAL VARIABLE INDICATING WHETHER H**2 CONVERGENCE OF
!          THE TRAPEZOID SUMS FOR THIS INTERVAL IS RECOGNIZED,
!  AITKEN  A LOGICAL VARIABLE INDICATING WHETHER CONVERGENCE OF RATIOS
!          FOR THIS SUBINTERVAL IS RECOGNIZED
!  T       CONTAINS THE FIRST *L* ROWS OF THE ROMBERG T-TABLE FOR THIS
!          SUBINTERVAL IN ITS LOWER TRIANGULAR PART. SPECIFICALLY,
!            T(I,1) = TRAPEZOID SUM (WITHOUT THE FACTOR *STEP*)
!                     ON 2**(I-1) + 1 EQUISPACED POINTS, I=1,...,L,
!            T(I,J+1) = T(I,J) + (T(I,J)-T(I-1,J))/(4**J - 1),
!                       J=2,...,I-1, I=2,...,L.
!          FURTHER, THE STRICTLY UPPER TRIANGULAR PART OF T CONTAINS
!          THE RATIOS FOR THE VARIOUS COLUMNS OF THE T-TABLE.
!          SPECIFICALLY,
!            T(J,I) = (T(I,J)-T(I-1,J))/(T(I+1,J)-T(I,J)),
!                     I=J+1,...,L-1,  J=1,...,L-2.
!          FINALLY, THE LAST OR L-TH COLUMN CONTAINS
!            T(J,L) = T(L,J) - T(L-1,J), J=1,...,L-1.
    6 STEP = END - BEG
      ASTEP = ABS (STEP)
      IF (ASTEP < STEPMN) GOTO 950
      IF (LEVEL >= 3) write (*, 609) BEG, STEP, ISTAGE
  609 FORMAT(" BEG,STEP ",2E16.8,I5)
      T (1, 1) = FBEG + FEND
      TABS = ABS (FBEG) + ABS (FEND)
      L = 1
      N = 1
      H2CONV = .FALSE.
      AITKEN = .FALSE.
      GOTO 10
!
    9 IF (LEVEL >= 4) write (*, 692) L, T (1, LM1)
   10 LM1 = L
      L = L + 1
!
!CALCULATE THE NEXT TRAPEZOID SUM, T(L,1), WHICH IS BASED ON
!  *N2*  + 1 EQUISPACED POINTS. HERE,  N2 = N*2 = 2**(L-1) .
      N2 = N * 2
      FN = N2
      ISTEP = (IEND-IBEG) / N
      IF (ISTEP > 1) GOTO 12
        II = IEND
        IEND = IEND+N
        IF (IEND > MAXTS) GOTO 900
        HOVN = STEP / FN
        III = IEND
        i = 1
   11   if ( i > n2 ) go to 1111
           TS (III) = TS (II)
           answer = END-FLOAT (I) * HOVN
           iflag = -3
           return
 1003      ts (iii - 1) = answer
           III = III - 2
           II = II - 1
           i = i + 2
        go to 11
 1111   continue
        ISTEP = 2
   12 ISTEP2 = IBEG + ISTEP / 2
      SUM = 0.0_rk
      SUMABS = 0.0_rk
      DO I = ISTEP2, IEND, ISTEP
         SUM = SUM + TS (I)
         SUMABS = SUMABS + ABS (TS (I) )
      end do
      T (L, 1) = 0.5_rk * T (L - 1, 1) + SUM / FN
      TABS = 0.5_rk * TABS + SUMABS / FN
      ABSI = ASTEP * TABS
      N = N2
!
!  GET PRELIMINARY VALUE FOR *VINT* FROM LAST TRAPEZOID SUM AND
!  UPDATE THE ERROR REQUIREMENT *ERGOAL* FOR THIS SUBINTERVAL.
!  THE ERROR REQUIREMENT IS NOT PRORATED ACCORDING TO THE LENGTH OF
!  THE CURRENT SUBINTERVAL RELATIVE TO THE INTERVAL OF INTEGRATION,
!  BUT ACCORDING TO THE HEIGHT *ISTAGE* OF THE CURRENT SUBINTERVAL
!  IN THE STACK OF SUBINTERVALS YET TO BE DONE.
!  THIS PROCEDURE IS NOT BACKED BY ANY RIGOROUS ARGUMENT, BUT
!  SEEMS TO WORK.
      IT = 1
      VINT = STEP * T (L, 1)
      TABTLM = TABS * TOLMCH
      FNSIZE = max (FNSIZE, ABS (T (L, 1) ) )
      ERGOAL = max (ASTEP * TOLMCH * FNSIZE,                            &
     &              STAGE * max (ERRA, ERRR * ABS (CUREST + VINT) ) )
!
!COMPLETE ROW L AND COLUMN L OF *T* ARRAY.
      FEXTRP = 1.0_rk
      DO I = 1, LM1
         FEXTRP = FEXTRP * 4.0_rk
         T (I, L) = T (L, I) - T (L - 1, I)
         T (L, I + 1) = T (L, I) + T (I, L) / (FEXTRP - 1.)
      end do
      ERRER = ASTEP * ABS (T (1, L) )
!---------------------------------------------------------------------
!---  PRELIMINARY DECISION PROCEDURE  --------------------------------
!
!  IF L = 2 AND T(2,1) = T(1,1), GO TO 60 TO FOLLOW UP THE IMPRESSION
!  THAT INTEGRAND IS STRAIGHT LINE.
      IF ( L <= 2 ) then
         IF (ABS (T (1, 2) ) <= TABTLM) GOTO 60
         GOTO 10
      end if
!
!CALCULATE NEXT RATIOS FOR COLUMNS 1,...,L-2 OF T-TABLE
!  RATIO IS SET TO ZERO IF DIFFERENCE IN LAST TWO ENTRIES OF
!  COLUMN IS ABOUT ZERO.
   15 DO I = 2, LM1
         DIFF = 0.0_rk
         IF (ABS (T (I - 1, L) ) > TABTLM) DIFF = T (I - 1, LM1)        &
     &                                          / T (I - 1, L)
         T (I - 1, LM1) = DIFF
      end do
!
!  T(1,LM1) IS THE RATIO DERIVED FROM LAST THREE TRAPEZOID SUMS, I.E.,
!     T(1,LM1) = (T(L-1,1)-T(L-2,1))/(T(L,1)-T(L-1,1)) .
!  IF THIS RATIO IS ABOUT 4, DO ROMBERG EXTRAPOLATION.
!  IF THIS RATIO IS ZERO, I.E., IF LAST TVO TRAPEZOID SUMS ARE ABOUT
!     EQUAL, MAYBE DO NOISE CHECK.
!  IF THIS RATIO IS ABOUT 2 IN ABSOLUTE VALUE, GO TO 50 WITH THE
!     BELIEF THAT INTEGRAND HAS JUMP DISCONTINUITY.
!  IF THIS RATIO IS, AT LEAST, ABOUT EQUAL TO THE PREVIOUS RATIO, THEN
!     THE INTEGRAND MAY WELL HAVE A NICE INTEGRABLE SINGULARITY.
!     GO TO 30 TO FOLLOW UP THIS HUNCH.
      IF (ABS (4.0_rk - T (1, LM1) ) > H2TOL) then
         IF (T (1, LM1) /= 0.0_rk) then
            IF (ABS (2.0_rk - ABS (T (1, LM1) ) ) < JUMPTL) GOTO 50
            IF (L == 3) GOTO 9
            H2CONV = .FALSE.
            IF (ABS((T(1, LM1) - T(1, L - 2) ) / T(1, LM1) ) <= AITTOL) &
     &      GOTO 30
!  AT THIS POINT, NO REGULAR BEHAVIOUR WAS DETECTED.
!  IF CURRENT SUBINTERVAL IS NOT REGULAR AND ONLY FOUR TRAPEZOID SUMS
!  WERE COMPUTED SO FAR, TRY ONE MORE TRAPEZOID SUM.
!  IF, AT LEAST, LAST TWO TRAPEZOID SUMS ARE ABOUT EQUAL, THEN
!  FAILURE TO RECOGNIZE REGULAR BEHAVIOUR MAY WELL BE DUE TO NOISE.
!  GO TO 70 TO CHECK THIS OUT.
!  OTHERWISE, GO TO 90 FOR FURTHER SUBDIVISION.
!
            IF ( .not. REGLAR ) then
               IF (L == 4) GOTO 9
            end if
         end if
         IF (ERRER <= ERGOAL) GOTO 70
         IF (LEVEL >= 4) write (*, 692) L, T (1, LM1)
         GOTO 91
      end if
!----------------------------------------------------------------------
!CAUTIOUS ROMBERG EXTRAPOLATION  --------------------------------------
!
!  THE CURRENT, OR L-TH, ROW OF THE ROMBERG T-TABLE HAS L ENTRIES.
!  FOR J=1,...,L-2, THE ESTIMATE
!             STEP*T(L,J+1)
!  IS BELIEVED TO HAVE ITS ERROR BOUNDED BY
!    ABS(STEP*(T(L,J)-T(L-1,J))/(4**J - 1))
!  IF THE LAST RATIO
!      T(J,LM1) = (T(L-1,J)-T(L-2,J))/(T(L,J)-T(L-1,J))
!  FOR COLUMN J OF THE T-TABLE IS ABOUT 4**J.
!  THE FOLLOWING IS A SLIGHTLY RELAXED EXPRESSION OF THIS BELIEF.
      IF (LEVEL >= 4) write (*, 619) L, T (1, LM1)
  619 FORMAT(I5,E16.8,5X,"H2CONV")
      IF ( .not. H2CONV ) then
         AITKEN = .FALSE.
         H2CONV = .TRUE.
         IF (LEVEL >= 3) write (*, 620) L
  620    FORMAT(" H2 CONVERGENCE AT ROW",I3)
      end if
   21 FEXTRP = 4.0_rk
      do
         IT = IT + 1
         VINT = STEP * T (L, IT)
         ERRER = ABS (STEP / (FEXTRP - 1.) * T (IT - 1, L) )
         IF (ERRER <= ERGOAL) GOTO 80
         IF (IT == LM1) GOTO 40
         IF (T (IT, LM1) /= 0.0_rk) then
            IF (T (IT, LM1) <= FEXTRP) GOTO 40
            IF (ABS (T (IT, LM1) / 4. - FEXTRP) / FEXTRP < AITTOL)      &
     &         FEXTRP = FEXTRP * 4.0_rk
         end if
      end do
!---------------------------------------------------------------------
!--- INTEGRAND MAY HAVE X**ALPHA TYPE SINGULARITY---------------------
!  RESULTING IN A RATIO OF *SING*  =  2**(ALPHA + 1)
   30 IF (LEVEL >= 4) write (*, 629) L, T (1, LM1)
  629 FORMAT(I5,E16.8,5X,"AITKEN")
      IF (T (1, LM1) < AITLOW) GOTO 91
      IF ( .not. AITKEN ) then
         H2CONV = .FALSE.
         AITKEN = .TRUE.
         IF (LEVEL >= 3) write (*, 630) L
  630    FORMAT(" AITKEN AT ROW",I3)
      end if
      FEXTRP = T (L - 2, LM1)
      IF (FEXTRP > 4.5_rk) GOTO 21
      IF (FEXTRP < AITLOW) GOTO 91
      IF (ABS (FEXTRP - T (L - 3, LM1) ) / T (1, LM1) > H2TOL) GOTO 91
      IF (LEVEL >= 3) write (*, 631) FEXTRP
  631 FORMAT(" RATIO",F12.8)
      SING = FEXTRP
      FEXTM1 = FEXTRP - 1.0_rk
      DO I = 2, L
         AIT (I) = T (I, 1) + (T (I, 1) - T (I - 1, 1) ) / FEXTM1
         R (I) = T (1, I - 1)
         DIF (I) = AIT (I) - AIT (I - 1)
      end do
      IT = 2
      do
         VINT = STEP * AIT (L)
         IF (LEVEL >= 5) then
            write (*, 632) (R (I + 1), I = IT, LM1)
            write (*, 632) (AIT (I), I = IT, L)
            write (*, 632) (DIF (I + 1), I = IT, LM1)
  632       FORMAT(1X,8E15.8)
         end if
         ERRER = ERRER / FEXTM1
         IF ( ERRER <= ERGOAL ) then
            ALPHA = LOG10 (SING) / ALG4O2 - 1.0_rk
            IF (LEVEL >= 2) write (*, 633) ALPHA, BEG, END
  633       FORMAT(12X,"INTEGRAND SHOWS SINGULAR BEHAVIOR OF TYPE",     &
     &             " X**(",F4.2,") BETWEEN",E15.8," AND",E15.8)
            FINAL_IFLAG = MAX (FINAL_IFLAG, 2)
            GOTO 80
         end if
         IT = IT + 1
         IF (IT == LM1) exit
         IF ( IT <= 3 ) then
            H2NEXT = 4.0_rk
            SINGNX = 2.0_rk * SING
         end if
         IF ( H2NEXT <= SINGNX ) then
            FEXTRP = SINGNX
            SINGNX = 2.0_rk * SINGNX
         else
            FEXTRP = H2NEXT
            H2NEXT = 4.0_rk * H2NEXT
         end if
         DO I = IT, LM1
            R (I + 1) = 0.0_rk
            IF (ABS (DIF (I + 1) ) > TABTLM) R (I + 1) = DIF (I) / DIF  &
     &        (I + 1)
         end do
         IF (LEVEL >= 4) write (*, 638) FEXTRP, R (L - 1), R (L)
  638    FORMAT(" FEXTRP + RATIOS",3E15.8)
         H2TFEX = - H2TOL * FEXTRP
         IF (R (L) - FEXTRP < H2TFEX) exit
         IF (R (L - 1) - FEXTRP < H2TFEX) exit
         ERRER = ASTEP * ABS (DIF (L) )
         FEXTM1 = FEXTRP - 1.0_rk
         DO I = IT, L
            AIT (I) = AIT (I) + DIF (I) / FEXTM1
            DIF (I) = AIT (I) - AIT (I - 1)
         end do
      end do
!CURRENT TRAPEZOID SUM AND RESULTING EXTRAPOLATED VALUES DID NOT GIVE
!  A SMALL ENOUGH *ERRER*.
!  IF LESS THAN FIVE TRAPEZOID SUMS WERE COMPUTED SO FAR, TRY NEXT
!     TRAPEZOID SUM.
!  OTHERWISE, DECIDE WHETHER TO GO ON OR TO SUBDIVIDE AS FOLLOWS.
!  WITH T(L,IT) GIVING THE CURRENTLY BEST ESTIMATE, GIVE UP ON DEVELOP
!  ING THE T-TABLE FURTHER IF  L > IT+2, I.E., IF EXTRAPOLATION
!  DID NOT GO VERY FAR INTO THE T-TABLE.
!  FURTHER, GIVE UP IF REDUCTION IN *ERRER* AT THE CURRENT RATE
!  DOES NOT PREDICT AN *ERRER* LESS THAN *ERGOAL* BY THE TIME
!  *MAXTBL* TRAPEZOID SUMS HAVE BEEN COMPUTED.
!  ---NOTE---
!  HAVING  PREVER < ERRER  IS AN ALMOST CERTAIN SIGN OF BEGINNING
!  TROUBLE WITH NOISE IN THE FUNCTION VALUES. HENCE,
!  A WATCH FOR, AND CONTROL OF, NOISE SHOULD BEGIN HERE.
   40 FEXTRP = max (PREVER / ERRER, AITLOW)
      PREVER = ERRER
      IF (L < 5) GOTO 10
      IF (LEVEL >= 3) write (*, 641) ERRER, ERGOAL, FEXTRP, IT
  641 FORMAT(" ERRER,ERGOAL,FEXTRP,IT",2E15.8,E14.5,I3)
      IF (L - IT > 2.AND.ISTAGE < MXSTGE) GOTO 90
      IF (ERRER / FEXTRP** (MAXTBL - L) < ERGOAL) GOTO 10
      GOTO 90
!---------------------------------------------------------------------
!----    INTEGRAND HAS JUMP (SEE NOTES)  -----------------------------
   50 IF (LEVEL >= 4) write (*, 649) L, T (1, LM1)
  649 FORMAT(I5,E16.8,5X,"JUMP")
      IF (ERRER > ERGOAL) GOTO 90
!  NOTE THAT  2*FN = 2**L
      DIFF = ABS (T (1, L) ) * 2.0_rk * FN
      IF (LEVEL >= 2) write (*, 650) DIFF, BEG, END
  650 FORMAT(12X,"INTEGRAND SEEMS TO HAVE JUMP OF SIZE",E13.6,          &
     &       " BETWEEN",E15.8," AND",E15.8)
      GOTO 80
!---------------------------------------------------------------------
!----  INTEGRAND IS STRAIGHT LINE  -----------------------------------
!  TEST THIS ASSUMPTION BY COMPARING THE VALUE OF THE INTEGRAND AT
!  FOUR 'RANDOMLY CHOSEN' POINTS WITH THE VALUE OF THE STRAIGHT LINE
!  INTERPOLATING THE INTEGRAND AT THE TWO ENDPOINTS OF THE SUB-
!  INTERVAL. IF TEST IS PASSED, ACCEPT *VINT* .
   60 IF (LEVEL >= 4) write (*, 660) L
  660 FORMAT(I5,21X,"STRAIGHT LINE")
      SLOPE = (FEND-FBEG) * 2.0_rk
      FBEG2 = FBEG * 2.0_rk
      i = 1
   61 if ( i > 4 ) go to 1161
         answer = BEG + RN (I) * STEP
         iflag = -4
         return
 1004    diff = abs ( answer - FBEG2 - RN (I) * SLOPE)
         IF (DIFF > TABTLM) GOTO 72
         i = i + 1
         go to 61
 1161 continue
      IF (LEVEL >= 3) write (*, 667) BEG, END
  667 FORMAT(27X,"INTEGRAND SEEMS TO BE STRAIGHT LINE BETWEEN",         &
     &       E15.8," AND",E15.8)
      GOTO 80
!---------------------------------------------------------------------
!------  NOISE MAY BE DOMINANT FEATURE  ------------------------------
!  ESTIMATE NOISE LEVEL BY COMPARING THE VALUE OF THE INTEGRAND AT
!  FOUR 'RANDOMLY CHOSEN' POINTS WITH THE VALUE OF THE STRAIGHT LINE
!  INTERPOLATING THE INTEGRAND AT THE TWO ENDPOINTS. IF SMALL
!  ENOUGH, ACCEPT *VINT*.
   70 IF (LEVEL >= 4) write (*, 670) L, T (1, LM1)
  670 FORMAT(I5,E16.8,5X,"NOISE")
      SLOPE = (FEND-FBEG) * 2.0_rk
      FBEG2 = FBEG * 2.0_rk
      I = 1
   71 answer = BEG + RN (I) * STEP
      iflag = -5
      return
 1005 DIFF = abs( answer - FBEG2 - RN (I) * SLOPE)
   72 ERRER = max (ERRER, ASTEP * DIFF)
      IF (ERRER > ERGOAL) GOTO 91
      I = I + 1
      IF (I <= 4) GOTO 71
      IF (LEVEL >= 3) write (*, 671) BEG, END
  671 FORMAT(" NOISE BETWEEN ",E15.8," AND",E15.8)
      FINAL_IFLAG = 3
!---------------------------------------------------------------------
!--- INTEGRATION OVER CURRENT SUBINTERVAL SUCCESSFUL -----------------
!  ADD *VINT* TO *result* AND *ERRER* TO *ERROR*, THEN SET UP NEXT
!  SUBINTERVAL, IF ANY.
   80 result = result+VINT
      ERROR = ERROR + ERRER
      IF (LEVEL >= 3) then
        IF (LEVEL >= 5) then
           DO I = 1, L                                
             write (*, 692) I, (T (I, J), J = 1, L)   
           end do                                     
        end if
        write (*, 682) VINT, ERRER, L, IT          
  682   FORMAT(" INTEGRAL IS",E16.8,", ERROR",E15.8,"  FROM T(",        &
     &         I0,",",I0,")")
      end if
!
      IF ( .not. RIGHT ) then
         ISTAGE = ISTAGE-1
         IF (ISTAGE == 0) then
           answer = result
           iflag = final_iflag
           RETURN
         end if
!
         REGLAR = REGLSV (ISTAGE)
         BEG = BEGIN (ISTAGE)
         END = FINIS (ISTAGE)
         CUREST = CUREST - EST (ISTAGE+1) + VINT
         IEND = IBEG - 1
         FEND = TS (IEND)
         IBEG = IBEGS (ISTAGE)
         GOTO 94
      end if
      CUREST = CUREST + VINT
      STAGE = STAGE * 2.0_rk
      IEND = IBEG
      IBEG = IBEGS (ISTAGE)
      END = BEG
      BEG = BEGIN (ISTAGE)
      FEND = FBEG
      FBEG = TS (IBEG)
      GOTO 5
!---------------------------------------------------------------------
!--- INTEGRATION OVER CURRENT SUBINTERVAL IS UNSUCCESSFUL-------------
!  MARK SUBINTERVAL FOR FURTHER SUBDIVISION. SET UP NEXT SUBINTERVAL.
   90 REGLAR = .TRUE.
   91 IF (ISTAGE == MXSTGE) GOTO 950
      IF ( LEVEL >= 5 ) then
         DO I = 1, L
            write (*, 692) I, (T (I, J), J = 1, L)
         end do
  692    FORMAT(I5,7E16.8/3E16.8)
      end if
      IF (RIGHT) GOTO 95
      REGLSV (ISTAGE+1) = REGLAR
      BEGIN (ISTAGE) = BEG
      IBEGS (ISTAGE) = IBEG
      STAGE = 0.5_rk * STAGE
   94 RIGHT = .TRUE.
      BEG = 0.5_rk * (BEG + END)
      IBEG = (IBEG + IEND) / 2
      TS (IBEG) = 0.5_rk * TS (IBEG)
      FBEG = TS (IBEG)
      GOTO 6
   95 NNLEFT = IBEG - IBEGS (ISTAGE)
      IF ( IEND+NNLEFT < MAXTS ) then
         III = IBEGS (ISTAGE)
         II = IEND
         DO I = III, IBEG
            II = II + 1
            TS (II) = TS (I)
         end do
         DO I = IBEG, II
            TS (III) = TS (I)
            III = III + 1
         end do
         IEND = IEND+1
         IBEG = IEND-NNLEFT
         FEND = FBEG
         FBEG = TS (IBEG)
         FINIS (ISTAGE) = END
         END = BEG
         BEG = BEGIN (ISTAGE)
         BEGIN (ISTAGE) = END
         REGLSV (ISTAGE) = REGLAR
         ISTAGE = ISTAGE+1
         REGLAR = REGLSV (ISTAGE)
         EST (ISTAGE) = VINT
         CUREST = CUREST + EST (ISTAGE)
         GOTO 5
      end if
!---------------------------------------------------------------------
!--- FAILURE TO HANDLE GIVEN INTEGRATION PROBLEM----------------------
  900 IF (LEVEL >= 2) write (*, 6900) BEG, END
 6900 FORMAT(" TOO MANY FUNCTION EVALUATIONS AROUND"/                   &
     &       10X,E15.8," AND",E15.8)
      FINAL_IFLAG = 4
      GOTO 999
  950 FINAL_IFLAG = 5
      IF ( LEVEL >= 2 ) then
         IF ( LEVEL >= 5 ) then
            DO I = 1, L
               write (*, 692) I, (T (I, J), J = 1, L)
            end do
         end if
         write (*, 6959) BEG, END
 6959    FORMAT(12X,"INTEGRAND SHOWS SINGULAR BEHAVIOUR OF "            &
     &          ,"UNKNOWN TYPE BETWEEN",E15.8," AND",E15.8)
      end if
  999 result = CUREST + VINT
      answer = result
      iflag = final_iflag
      RETURN
      END subroutine CADRE_Reverse

      logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
      character (len=*), parameter :: IdParm = &
           "$Id$"
      character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
        not_used_here = (id(1:1) == ModuleName(1:1))
        print *, Id ! .mod files sometimes change if PRINT is added
      end function not_used_here

end module Cadre_m

! $Log$
