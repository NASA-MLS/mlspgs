! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MathUtils ! Math utilities
!=============================================================================

  USE ERMSG_M, ONLY: DERM1, DERV1, ERFIN, ERMSG, SERM1, SERV1

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  CONTAINS

  ! --------------------------------------------------------  AMACH  -----
  SUBROUTINE AMACH (MODE, I, I1, R1, D1)
!>> 1996-03-30 AMACH  Krogh   Added external statement.
!>> 1994-10-26 AMACH  Krogh   Changes to use M77CON
!>> 1994-09-23 AMACH  Snyder  Add VAX G parameters
!>> 1994-06-21 AMACH  Snyder  Compute only round-off and u-flow at first
!>> 1994-05-25 AMACH  Snyder  Added an option to compute at run time.
!>> 1992-04-07 AMACH  Oken    Removed ^Z at EOF (error found by VAX comp
!>> 1992-02-20 AMACH  Snyder  Added Cray-YMP stuff, q.v.
!>> 1990-06-11 AMACH  Snyder  Added Apollo DN-10000 stuff, q.v.
!>> 1990-12-14 AMACH  Lawson  Changed to eliminate ENTRY statements.
!>> 1990-08-21 AMACH  Krogh   No test was getting done for bad machine.
!>> 1990-02-28 AMACH  Krogh   Correct missing DOUBLE PRECISION AMSUB1
!>> 1989-08-14 AMACH  Krogh   Parameterized everything -- Massive change
!>> 1989-03-30 AMACH  Snyder  Correct missing "/" line 921
!>> 1989-01-30 AMACH  Snyder  Incorporate more constants from NETLIB.
!>> 1988-05-19 AMACH  Lawson  Initial code.
! File AMACH.FOR contains user-callable functions I1MACH, D1MACH, and
! R1MACH, plus second-level subroutines AMACH, AMTEST, and AMSUB1.
! Appropriate lines must be switched between comment and non-comment
! status when this code is moved to a different computer system.
!     These changes can be done with any text editor, however the "c++"
! lines permit automation of the change using the M77CON processor.
! Note that when the M77CON processor activates a line it shifts
! Columns 2-72 to 1-71 and puts a blank in Column 72.  When it inactiv-
! ates a line it shifts Columns 1-71 to 2-72 and puts a C in Column 1.
!     The possible choices using M77CON (don't include parenthetical
!     comments) are:
!      c++ CURRENT HAS SYS = IEEE
!      c++ CURRENT HAS SYS = ALPHA_D3
!      c++ CURRENT HAS SYS = AMDAHL
!      c++ CURRENT HAS SYS = APOLLO_10000
!      c++ CURRENT HAS SYS = BUR1700
!      c++ CURRENT HAS SYS = BUR5700
!      c++ CURRENT HAS SYS = BUR67_7700
!      c++ CURRENT HAS SYS = CDC60_7000
!      c++ CURRENT HAS SYS = CONVEXC_1
!      c++ CURRENT HAS SYS = CRAY1
!      c++ CURRENT HAS SYS = CRAY1_SD (Sngl prec.arith. used for dble.)
!      c++ CURRENT HAS SYS = CRAY1_64 (64 bit integers)
!      c++ CURRENT HAS SYS = CRAY1_SD_64 (64 bit int, SP used for DP)
!      c++ CURRENT HAS SYS = CRAY_T3D
!      c++ CURRENT HAS SYS = CRAY_YMP
!      c++ CURRENT HAS SYS = CRAY_YMP_SD (Sngl prec. used for dble.)
!      c++ CURRENT HAS SYS = DG_S2000
!      c++ CURRENT HAS SYS = HARRIS220
!      c++ CURRENT HAS SYS = HON600_6000
!      c++ CURRENT HAS SYS = HON_DPS_8_70
!      c++ CURRENT HAS SYS = HP700Q
!      c++ CURRENT HAS SYS = IBM360_370
!      c++ CURRENT HAS SYS = INTERDATA_8_32
!      c++ CURRENT HAS SYS = PDP10_KA
!      c++ CURRENT HAS SYS = PDP10_KB
!      c++ CURRENT HAS SYS = PDP11
!      c++ CURRENT HAS SYS = PRIME50
!      c++ CURRENT HAS SYS = SEQ_BAL_8000
!      c++ CURRENT HAS SYS = UNIVAC
!      c++ CURRENT HAS SYS = VAX
!      c++ CURRENT HAS SYS = VAX_G
!     The current choice is:
!++ CURRENT HAS SYS = IEEE
!
!     One can also select whether floating point constants are created
!     by the compiler or created at run time.  The choices using M77CON
!     are:
!      c++ CURRENT HAS HOW = COMPILER
!      c++ CURRENT HAS HOW = RUN
!     The current choice is:
!++ CURRENT HAS HOW = COMPILER
!
!     If the constants are created at run time, and they fail the run-
!     time check for reasonableness, they are re-created assuming IEEE.
!     If they still fail, the program stops.
!
!  I/O UNIT NUMBERS:
!
!    IM1 = I1MACH( 1) = THE STANDARD INPUT UNIT.
!    IM2 = I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    IM3 = I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    IM4 = I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!
!  WORDS:
!
!    IM5 = I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    IM6 = I1MACH( 6) = THE NUMBER OF CHARACTERS/INTEGER STORAGE UNIT.
!
!  INTEGERS:
!
!    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
!
!               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
!
!    IM7 = I1MACH( 7) = A, THE BASE.
!    IM8 = I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    IM9 = I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS:
!
!    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
!    BASE-B FORM
!
!               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
!               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
!
!    IM10 = I1MACH(10) = B, THE BASE.
!
!  SINGLE-PRECISION:
!
!    IM11 = I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    IM12 = I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    IM13 = I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION:
!
!    IM14 = I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    IM15 = I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    IM16 = I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
!  CONVERSION FROM FUNCTIONAL TO STRUCTURAL FLOATING POINT CONSTANTS
!
!    IM17 = CONSTANT SUCH THAT IM14 + IM17 = ACTUAL NUMBER OF BASE-B
!           DIGITS IN DOUBLE PRECISION, USED FOR CHECKING THAT CORRECT
!           VERSION OF THIS PROGRAM IS INSTALLED.  (SEE DEFINITION OF
!           DM6, AND THE USE OF DM6 IN CALLING AMTEST.)
!
!  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
!  THE DESIRED SET OF PARAMETER STATEMENTS SHOULD BE ACTIVATED BY
!  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
!  IM1 - IM4 SHOULD BE CHECKED FOR CONSISTENCY
!  WITH THE LOCAL OPERATING SYSTEM.
!     -----------------------------------------------------------------
!     Original design and code due to P. A. Fox, A. D. Hall, and
!     N. L. Schryer, Bell Laboratories.  See ACM TOMS, 4,(1978),177-188.
!     Adapted to Univac 1100 by Kris Stewart, JPL, 7/30/81.
!     Adapted for the JPL MATH77 library by C. L. Lawson and F. T. Krogh
!     Sept, 1987.
!     1989-08-14 AMACH  Krogh   Parameterized everything. Major changes.
!     1990 Dec. CLL reorganized code to avoid using ENTRY statements
!     for functions of different types.  Also added save statements.
!     -----------------------------------------------------------------
!     On the first call to this function, tests are done to verify that
!     IM10 and IM14 are not grossly wrong for the host environment.
!     This gives some protection against using the wrong version of this
!     subprogram.
!     -----------------------------------------------------------------
      INTEGER MODE, I, I1
      REAL R1
      DOUBLE PRECISION D1, TEST
!
      INTEGER IMACH(17)
      INTEGER IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10, IM11, &
                  IM12, IM13, IM14, IM15, IM16, IM17
!++ Code for HOW=RUN is INACTIVE
!      integer IEEE
!      integer ID1, ID2, ID3, ID4, ID5, ID6, ID7, ID8, ID10, ID11,
!     1   ID12, ID13, ID14, ID15, ID16, ID17
!++ Code for (HOW=RUN) | SYS=IEEE is ACTIVE
      INTEGER IE1, IE2, IE3, IE4, IE5, IE6, IE7, IE8, IE10, IE11, &
         IE12, IE13, IE14, IE15, IE16, IE17
!++ end
      REAL             RMACH(5), RM1, RM2, RM3, RM4, RM5, &
                       RMA, RMB, RBASE
      DOUBLE PRECISION DMACH(5), DM1, DM2, DM3, DM4, DM5, DM6, &
                       DMA, DMB, DBASE
      SAVE TEST, IMACH, RMACH, DMACH
!     -----------------------------------------------------------------
!     Machine constants for IEEE standard binary floating-point
!     processors.  This includes PC's and work-stations using the
!     Intel 8087, 80287, 80387, ... processors or the
!     Motorola 68881, 68882, ... processors.
!     Note:  We are setting the "most negative exponent" (IMACH(12) and
!     IMACH(15)) to be the exponent of the smallest normalized number.
!     An IEEE processor actually handles smaller numbers before
!     underflowing, however these "unnormalized" numbers have
!     diminished precision.
!
!++ Code for (HOW=RUN) | SYS=IEEE is ACTIVE
!     Parameters for IEEE when generating at run time:
      PARAMETER (IE1 =5, IE2 =6, IE3 =7, IE4 =6)
      PARAMETER (IE5 =32, IE6 =4, IE7 =2, IE8 =31)
      PARAMETER (IE10 =2, IE11 =24, IE12 =-125, IE13 =128)
      PARAMETER (IE14 =53, IE15 =-1021, IE16 =1024, IE17=0)
!++ Code for SYS = IEEE is ACTIVE
      PARAMETER (IM1 = IE1, IM2 = IE2, IM3 = IE3, IM4 = IE4)
      PARAMETER (IM5 = IE5, IM6 = IE6, IM7 = IE7, IM8 = IE8)
      PARAMETER (IM10 = IE10, IM11 = IE11, IM12 = IE12, IM13 = IE13)
      PARAMETER (IM14 = IE14, IM15 = IE15, IM16 = IE16, IM17 = IE17)
!     -----------------------------------------------------------------
!++ Code for SYS = ALPHA_D3 is INACTIVE
!c     MACHINE CONSTANTS for the VAX/VMS F and D-3 format for Alpha

!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-127, IM16 =127, IM17=0)
!++ end
!     -----------------------------------------------------------------
!++ Code for HOW = RUN is INACTIVE
!c     MACHINE CONSTANTS for the VAX/VMS F and D-3 format for Alpha

!      PARAMETER (ID1 =5, ID2 =6, ID3 =7, ID4 =6)
!      PARAMETER (ID5 =32, ID6 =4, ID7 =2, ID8 =31)
!      PARAMETER (ID10 =2, ID11 =24, ID12 =-127, ID13 =127)
!      PARAMETER (ID14 =53, ID15 =-127, ID16 =127, ID17=0)
!++ end
!     -----------------------------------------------------------------
!++ Code for SYS = AMDAHL is INACTIVE
!C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!      -----------------------------------------------------------------
!++ Code for SYS = APOLLO_10000 is INACTIVE
!c     MACHINE CONSTANTS FOR APOLLO DN_10000 MACHINES.
!c     The only difference from IEEE is IM13.  This difference has
!c     nothing to do with the arithmetic or representation used by the
!c     machine.  It is caused by a bug in the compiler:  The right-hand
!c     side of RM2 (below) is apparently evaluated in double precision.
!c     When the compiler is ready to store the resulting value into its
!c     internal data structures, it compares it to an incorrect value
!c     of the overflow limit.  It appears the incorrect value has the
!c     correct exponent, but the fraction is 1.5 instead of 2-2**(-p),
!c     where p is the precision in bits.  You can get the correct result
!c     by changing IM13 to 128, changing RM2 from a parameter to a
!c     variable, and changing the parameter statement that assigns a
!c     value to RM2 into an ordinary assignment statement.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17 =0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR1700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
!C
!      PARAMETER (IM1 =7, IM2 =2, IM3 =2, IM4 =2)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =33)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-256, IM13 =255)
!      PARAMETER (IM14 =60, IM15 =-256, IM16 =255, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR5700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
!      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
!      PARAMETER (IM14 =26, IM15 =-50, IM16 =76, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = BUR67_7700 is INACTIVE
!C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =48, IM6 =6, IM7 =2, IM8 =39)
!      PARAMETER (IM10 =8, IM11 =13, IM12 =-50, IM13 =76)
!      PARAMETER (IM14 =26, IM15 =-32754, IM16 =32780, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CDC60_7000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =60, IM6 =10, IM7 =2, IM8 =48)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-929, IM13 =1070)
!      PARAMETER (IM14 =94, IM15 =-929, IM16 =1069, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CONVEXC_1 is INACTIVE
!C     MACHINE CONSTANTS FOR CONVEX C-1.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-1024, IM16 =1023, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY_T3D is INACTIVE
!c     Machine constants for Cray T3D.  IEEE double for both precisions.
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =53, IM12 =-1021, IM13 =1024)
!      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY_YMP is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY YMP
!C     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
!C     the underflow exponent (IM12 and IM15) is -8189, but these values
!C     don't seem to work in cf77:  the underflow limit underflows, and
!C     the overflow limit overflows when using Cray's values.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
!      PARAMETER (IM14 =94, IM15 =-8188, IM16 =8189, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY_YMP_SD is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY YMP
!C     Cray claims the overflow exponent (IM13 and IM16) is 8189, and
!C     the underflow exponent (IM12 and IM15) is -8189, but these
!C     values don't seem to work in cf77:  the underflow limit under-
!C     flows, and the overflow limit overflows when using Cray's values.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8188, IM13 =8189)
!      PARAMETER (IM14 =47, IM15 =-8188, IM16 =8189, IM17=1)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1_SD is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
!C     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =46)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1_64 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =94, IM15 =-8099, IM16 =8190, IM17=2)
!C     -----------------------------------------------------------------
!++ Code for SYS = CRAY1_SD_64 is INACTIVE
!C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3, WHEN DOUBLE
!C     PRECISION IS TO USE SINGLE PRECISION ARITHMETIC.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =102, IM4 =6)
!      PARAMETER (IM5 =64, IM6 =8, IM7 =2, IM8 =63)
!      PARAMETER (IM10 =2, IM11 =47, IM12 =-8189, IM13 =8190)
!      PARAMETER (IM14 =47, IM15 =-8189, IM16 =8190, IM17=1)
!C     -----------------------------------------------------------------
!++ Code for SYS = DG_S2000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
!C
!      PARAMETER (IM1 =11, IM2 =12, IM3 =8, IM4 =10)
!      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HARRIS220 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HARRIS 220, SLASH 6, SLASH 7.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =0, IM4 =6)
!      PARAMETER (IM5 =24, IM6 =3, IM7 =2, IM8 =23)
!      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =38, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HON600_6000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =6, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HON_DPS_8_70 is INACTIVE
!C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =43, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =63, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = HP700Q is INACTIVE
!c     Machine constants for HP-700 using the +autodblpad option,
!c     which automatically increases DOUBLE PRECISION to REAL*16, and
!c     REAL to DOUBLE PRECISION.
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =53, IM12 =-1021, IM13 =1024)
!      PARAMETER (IM14 = 113, IM15 = -16381, IM16 = 16384, IM17 = 0)
!C     -----------------------------------------------------------------
!++ Code for SYS = IBM360_370 is INACTIVE
!C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
!C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.

!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =63)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =63, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = INTERDATA_8_32 is INACTIVE
!C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
!C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =6, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =16, IM11 =6, IM12 =-64, IM13 =62)
!      PARAMETER (IM14 =14, IM15 =-64, IM16 =62, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP10_KA is INACTIVE
!C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =54, IM15 =-101, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP10_KB is INACTIVE
!C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
!C
!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =5, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =62, IM15 =-128, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PDP11 is INACTIVE
!C     MACHINE CONSTANTS FOR PDP-11 FORTRAN'S SUPPORTING
!C     16-BIT INTEGER ARITHMETIC.

!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =16, IM6 =2, IM7 =2, IM8 =15)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = PRIME50 is INACTIVE
!C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
!C     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
!C     SUPPLIED BY IGOR BRAY.

!      PARAMETER (IM1 =1, IM2 =1, IM3 =2, IM4 =1)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =23, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =47, IM15 =-32895, IM16 =32637, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = SEQ_BAL_8000 is INACTIVE
!C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
!C
!      PARAMETER (IM1 =0, IM2 =0, IM3 =7, IM4 =0)
!      PARAMETER (IM5 =32, IM6 =1, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-125, IM13 =128)
!      PARAMETER (IM14 =53, IM15 =-1021, IM16 =1024, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = UNIVAC is INACTIVE
!C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!C
!C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 1
!C     WHICH IS APPROPRIATE FOR THE UNIVAC-FTN SYSTEM.
!C     IF YOU HAVE THE UNIVAC-FOR SYSTEM, SET IT TO 7.
!C     IM6 = 4 for FTN (4 chars per word), 6 for FOR (6 chars per word).
!c
!      PARAMETER (IM1 =5, IM2 =6, IM3 =1, IM4 =6)
!      PARAMETER (IM5 =36, IM6 =4, IM7 =2, IM8 =35)
!      PARAMETER (IM10 =2, IM11 =27, IM12 =-128, IM13 =127)
!      PARAMETER (IM14 =60, IM15 =-1024, IM16 =1023, IM17=0)
!C     -----------------------------------------------------------------
!++ Code for SYS = VAX is INACTIVE
!c     MACHINE CONSTANTS for the VAX/VMS F and D formats
!c     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.

!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =56, IM15 =-127, IM16 =127, IM17=0)
!++ end
!     -----------------------------------------------------------------
!++ Code for SYS = VAX_G is INACTIVE
!c     MACHINE CONSTANTS for the VAX/VMS F and G formats
!c     and for PDP-11 FORTRAN SUPPORTING 32-BIT INTEGER ARITHMETIC.

!      PARAMETER (IM1 =5, IM2 =6, IM3 =7, IM4 =6)
!      PARAMETER (IM5 =32, IM6 =4, IM7 =2, IM8 =31)
!      PARAMETER (IM10 =2, IM11 =24, IM12 =-127, IM13 =127)
!      PARAMETER (IM14 =53, IM15 =-1023, IM16 =1023, IM17=0)
!++ end
!     -----------------------------------------------------------------
!
!
! Real parameters
!
!  RM1 = R1MACH(1) = B**(EMIN-1), The smallest positive number, i.e.,
!                    the underflow limit.
!  RM2 = R1MACH(2) = B**EMAX*(1 - B**(-T)), The largest number, i.e.,
!                    the overflow limit.
!  RM3 = R1MACH(3) = B**(-T), The smallest relative spacing, i.e., the
!                    difference between 1.0 and the next smaller number.
!  RM4 = R1MACH(4) = B**(1-T), The largest relative spacing, i.e., the
!                     difference between 1.0 and the next larger number.
!  RM5 = R1MACH(5) = LOG10(B).  When B = 2 this value is
!              Log10(2) = 0.30102_99956_63981_19521_37388_94724
!
! Parameter RMA and RMB are selected so that for values of the base =
! 2, 8, 16, 10, RMA has the values 1, 3, 4, 0, and RMB has the values 0,
! 0, 0, 1.  These values are used in computing RM5.
! $$$$ Note that if other bases are to be supported, the calculation of
! $$$$ RMA and RMB will have to be generalized.
!
!++   Code for HOW = COMPILER is ACTIVE
      PARAMETER (IM9 = 2 * (2**(IM8-1) - 1) + 1)
      PARAMETER (RMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 + &
          12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (RMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (RBASE = IM10)
!
!     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
!     some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and RM2 using
!     these subterfuges so it will be clear we're computing the REAL
!     and DOUBLE PRECISION characteristics in the same way.
      PARAMETER (RM1 = (RBASE**(IM12/2)) * (RBASE**(IM12-IM12/2-1)))
      PARAMETER (RM2 = RBASE**(IM13-IM11) * ((RBASE**IM11 - RBASE) &
                     + (RBASE - 1.0E0)))
      PARAMETER (RM3 = RBASE**(-IM11))
      PARAMETER (RM4 = RBASE**(1-IM11))
!     PARAMETER (RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB)
      PARAMETER (RM5 = RMA*0.301029995663981195213738894724E0+RMB)
!
! Double precision parameters -- (Defined like the real ones.)
!
      PARAMETER (DMA = ((IM10 - 10) * (-3 + ((IM10 - 2) * (-77 + &
          12 * (IM10 - 8))) / 14)) / 24)
      PARAMETER (DMB = ((IM10 - 2) * (IM10 - 8) * (16 - IM10)) / 96)
      PARAMETER (DBASE = IM10)
!
!     Weird subterfuges below are NECESSARY to compute DM1 and DM2 on
!     some systems.  DON'T SIMPLIFY THEM.
      PARAMETER (DM1 = (DBASE**(IM15/2)) * (DBASE**(IM15-IM15/2-1)))
      PARAMETER (DM2 = DBASE**(IM16-IM14) * ((DBASE**IM14 - DBASE) &
                     + (DBASE - 1.0D0)))
      PARAMETER (DM3 = DBASE**(-IM14))
      PARAMETER (DM4 = DBASE**(1-IM14))
!     PARAMETER (DM5 = DMA *
!    1 0.30102 99956 63981 19521 37388 94724 49302 67681 89881 46211 D0
!    2 + DMB)
!
      PARAMETER (DM5 = DMA* &
       0.30102999566398119521373889472449302676818988146211D0 + DMB)
! DM6 and TEST are used in checking that the correct constants have
! been selected.
      PARAMETER (DM6 = DBASE**(-IM14-IM17))
!++   END
      DATA TEST / 0.D0 /
!
!     DATA IMACH / IM1, IM2, IM3, IM4, IM5, IM6, IM7, IM8, IM9, IM10,
!    1   IM11, IM12, IM13, IM14, IM15, IM16 /
!     DATA RMACH / RM1, RM2, RM3, RM4, RM5 /
!     DATA DMACH / DM1, DM2, DM3, DM4, DM5 /
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (TEST .EQ. 0.0D0) THEN
!         IM9 = 2 * (2**(IM8-1) - 1) + 1
         IMACH(1) = IM1
         IMACH(2) = IM2
         IMACH(3) = IM3
         IMACH(4) = IM4
         IMACH(5) = IM5
         IMACH(6) = IM6
         IMACH(7) = IM7
         IMACH(8) = IM8
         IMACH(10) = IM10
         IMACH(11) = IM11
         IMACH(12) = IM12
         IMACH(13) = IM13
         IMACH(14) = IM14
         IMACH(15) = IM15
         IMACH(16) = IM16
         IMACH(17) = IM17
!++   Code for HOW = RUN is INACTIVE
!         IEEE = 0
!100      continue
!      DBASE = IMACH(10)
!C
!C     Weird subterfuge below is NECESSARY to compute DM1 on
!C     some systems.  DON'T SIMPLIFY IT.
!      DM1=(DBASE**(IMACH(15)/2)) * (DBASE**(IMACH(15)-IMACH(15)/2-1))
!C DM6 and TEST are used in checking that the correct constants have
!C been selected.
!      DM6 = DBASE**(-IMACH(14)-IMACH(17))
!++   end
         CALL AMTEST (TEST, DM6)
         IF (dm1 .EQ. 0.0d0 .OR. test .EQ. 0.0d0) THEN
!++   Code for HOW = RUN is INACTIVE
!           if (IEEE .eq. 0) then
!              IEEE = 1
!              IMACH(1) = IE1
!              IMACH(2) = IE2
!              IMACH(3) = IE3
!              IMACH(4) = IE4
!              IMACH(5) = IE5
!              IMACH(6) = IE6
!              IMACH(7) = IE7
!              IMACH(8) = IE8
!              IMACH(10) = IE10
!              IMACH(11) = IE11
!              IMACH(12) = IE12
!              IMACH(13) = IE13
!              IMACH(14) = IE14
!              IMACH(15) = IE15
!              IMACH(16) = IE16
!              IMACH(17) = IE17
!              go to 100
!           end if
!           if (IEEE .eq. 1) then
!              IEEE = 2
!              IMACH(1) = ID1
!              IMACH(2) = ID2
!              IMACH(3) = ID3
!              IMACH(4) = ID4
!              IMACH(5) = ID5
!              IMACH(6) = ID6
!              IMACH(7) = ID7
!              IMACH(8) = ID8
!              IMACH(10) = ID10
!              IMACH(11) = ID11
!              IMACH(12) = ID12
!              IMACH(13) = ID13
!              IMACH(14) = ID14
!              IMACH(15) = ID15
!              IMACH(16) = ID16
!              IMACH(17) = ID17
!              go to 100
!           end if
!++   END
            PRINT*,'AMACH has bad parameters for current environment.'
            STOP
         END IF
!++   Code for HOW = RUN is INACTIVE
!         IM9 = 2 * (2**(IMACH(8)-1) - 1) + 1
!         RMA = ((IMACH(10) - 10) * (-3 + ((IMACH(10) - 2) * (-77 +
!     1       12 * (IMACH(10) - 8))) / 14)) / 24
!         RMB = ((IMACH(10)-2) * (IMACH(10)-8) * (16-IMACH(10)))/96
!         RBASE = IMACH(10)
!C
!C        Weird subterfuges below are NECESSARY to compute DM1 and DM2
!C        on some systems.  DON'T SIMPLIFY THEM.  We compute RM1 and
!C        RM2 using these subterfuges so it will be clear we're
!C        computing the REAL and DOUBLE PRECISION characteristics in
!c        the same way.
!         RM1=(RBASE**(IMACH(12)/2))*(RBASE**(IMACH(12)-IMACH(12)/2-1))
!         RM2 = RBASE**(IMACH(13)-IMACH(11))*((RBASE**IMACH(11) - RBASE)
!     1                  + (RBASE - 1.0E0))
!         RM3 = RBASE**(-IMACH(11))
!         RM4 = RBASE**(1-IMACH(11))
!c        RM5 = RMA*0.30102 99956 63981 19521 37388 94724E0+RMB
!         RM5 = RMA*0.301029995663981195213738894724E0+RMB
!C
!C Double precision parameters -- (Defined like the real ones.)
!C
!         DMA = ((IMACH(10) - 10) * (-3 + ((IMACH(10) - 2) * (-77 +
!     1       12 * (IMACH(10) - 8))) / 14)) / 24
!         DMB = ((IMACH(10)-2) * (IMACH(10)-8) * (16-IMACH(10)))/96
!C
!C        Weird subterfuge below is NECESSARY to compute DM2 on
!C        some systems.  DON'T SIMPLIFY IT.
!         DM2 = DBASE**(IMACH(16)-IMACH(14))*((DBASE**IMACH(14) - DBASE)
!     1                  + (DBASE - 1.0D0))
!         DM3 = DBASE**(-IMACH(14))
!         DM4 = DBASE**(1-IMACH(14))
!c        DM5 = DMA*0.30102 99956 63981 19521 37388 94724D0+DMB
!         DM5 = DMA*0.301029995663981195213738894724D0+DMB
!++   END
         IMACH(9) = IM9
         RMACH(1) = RM1
         RMACH(2) = RM2
         RMACH(3) = RM3
         RMACH(4) = RM4
         RMACH(5) = RM5
         DMACH(1) = DM1
         DMACH(2) = DM2
         DMACH(3) = DM3
         DMACH(4) = DM4
         DMACH(5) = DM5
      END IF

      IF (MODE .EQ. 0) THEN
         I1=IMACH(I)
      ELSE IF (MODE .EQ. 1) THEN
         R1=RMACH(I)
!                                  Here we assume MODE = 2.
      ELSE
         D1=DMACH(I)
      END IF
      RETURN
  END SUBROUTINE AMACH
  ! -------------------------------------------------------  I1MACH  -----
  INTEGER FUNCTION I1MACH(I)
      INTEGER I, I1
      REAL R1
      DOUBLE PRECISION D1
      IF (I .LT. 1  .OR.  I .GT. 16) THEN
         PRINT*,'I1MACH.. Bad argument: I =',I
         STOP 'I1MACH error'
      END IF
      CALL AMACH (0, I, I1, R1, D1)
      I1MACH = I1
      RETURN
  END FUNCTION I1MACH

  ! -------------------------------------------------------  R1MACH  -----
  REAL FUNCTION R1MACH(I)
      INTEGER I, I1
      REAL R1
      DOUBLE PRECISION D1
      IF (I .LT. 1  .OR.  I .GT. 5) THEN
         PRINT*,'R1MACH.. Bad argument: I = ',I
         STOP 'R1MACH error'
      END IF
      CALL AMACH (1, I, I1, R1, D1)
      R1MACH = R1
      RETURN
  END FUNCTION R1MACH

  ! -------------------------------------------------------  D1MACH  -----
  DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I, I1
      REAL R1
      DOUBLE PRECISION D1
      IF (I .LT. 1  .OR.  I .GT. 5) THEN
         PRINT*,'D1MACH.. Bad argument: I = ',I
         STOP 'D1MACH error'
      END IF
      CALL AMACH (2, I, I1, R1, D1)
      D1MACH = D1
      RETURN
  END FUNCTION D1MACH

  ! -------------------------------------------------------  AMTEST  -----
  SUBROUTINE AMTEST (TEST, D6)
! Verifies that D6 is an appropriate value for DM6.
! Returns TEST = D6 + D6 - 1, .ne. 0 if D6 is an appropriate value for
! DM6, else returns TEST = 0.  The caller uses TEST = 0 as a signal to
! try again with IEEE settings (unless that's already been done).

      DOUBLE PRECISION D6, TEST
      TEST = AMSUB1(1.D0 + D6)
!
! The comparison with 1.875E0*D6 in the line below is to guard
! against the possibility that TEST is > 0 as a result of rounding
! up in the addition of D6 to 1.
!
      IF ((TEST .EQ. 0.D0) .OR. (TEST .GT. 1.875D0*D6)) THEN
         TEST = (D6 + D6) + 1.D0
         IF (AMSUB1(TEST) .NE. 0.D0) RETURN
      END IF
      test = 0.0d0
  END SUBROUTINE AMTEST

  ! -------------------------------------------------------  AMSUB1  -----
  DOUBLE PRECISION FUNCTION AMSUB1 (TEST1)
      DOUBLE PRECISION TEST1
!     Returns the value of TEST1 - 1.
      AMSUB1 = TEST1 - 1.0D0
      RETURN
  END FUNCTION AMSUB1

  ! --------------------------------------------------------  DERFI  -----
  DOUBLE PRECISION FUNCTION DERFI (X)
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!--D replaces "?": ?ERFI, ?ERFCI, ?ERFIX, ?ERM1
!>> 1996-06-18 DERFI Krogh  Changes to use .C. and C%%. J not changed.
!>> 1996-03-30 DERFI Krogh  Added external statements.
!>> 1995-11-28 DERFI Krogh  Removed multiple entries.
!>> 1995-11-03 DERFI Krogh  Removed blanks in numbers for C conversion.
!>> 1994-10-20 DERFI Krogh  Changes to use M77CON
!>> 1994-04-20 DERFI CLL Edited type stmts to make DP & SP files similar
!>> 1987-10-29 DERFI Snyder  Initial code.
!
!     For -1.0 .LT. X .LT. 1.0 calculate the inverse of the error
!     function.  That is, X = ERF(ERFI).
!
!     For 0.0 .LT. X .LT. 2.0 calculate the inverse of the
!     complementary error function.  that is, X = ERFC(ERFCI).  This
!     calculation is carried out by invoking the alternate entry *ERFCI.
!
!     If X is out of range, program execution is terminated by calling
!     the error message processor.
!
!     This subprogram uses approximations due to A. Strecok from
!     Mathematics of Computation 22, (1968) pp 144-158.
!

      DOUBLE PRECISION X
      DOUBLE PRECISION ARG, D(6), FSIGN, S
      INTEGER J
!
!     *****     Parameters     *****************************************
!
! MAX...  is the position in C of the last coefficient of a Chebyshev
!         polynomial expansion.
! MIN...  is the position in C of the first coefficient of a Chebyshev
!         polynomial expansion.
! NC      is the upper dimension of the array of coefficients.
! NDELTA  is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         0.9975 .LT. X .LE. 1-5.0D-16
! NLAMDA  is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         0.8 .LT. X .LE. 0.9975.
! NMU     is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         5.0D-16 .GT. 1-X .GE. 1.D-300.
! NXI     is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate DERFCI(X)/X in the
!         range 0.0 .LE. X .LE. 0.8.
!
!
!     *****     External References     ********************************
!
! D1MACH   Provides the round-off level.  Used to calculate the number
!          of coefficients to retain in each Chebyshev expansion.
! DERM1    Prints an error message and stops if X .LE. -1.0 or
!          X .GE. 1.0 (ERFI) or X .LE. 0.0 or X .GE. 2.0 (ERFCI).
! LOG      Calculates the natural logarithm.
! SQRT     Calculates the square root.
!
!
!     *****     Local Variables      ***********************************
!
! ARG     If ERFI or ERFCI is being approximated by a Chebyshev
!         expansion then ARG is the argument of ERFI or the argument
!         that would be used if ERFCI(X) were computed as ERFC(1-X),
!         that is, ARG = X if ERFI is being computed, or ARG = 1-X if
!         ERFCI is being computed.  If ERFI or ERFCI is being computed
!         using the initial approximation ERFI=SQRT(-LOG((1-X)*(1+X))),
!         then ARG is that initial approximation.
! C       contains the coefficients of polynomial expansions.  They are
!         stored in C in the order DELTA(0..37), LAMDA(0..26),
!         MU(0..25), XI(0..38).
! D       are used to scale the argument of the Chebyshev polynomial
!         expansion in the range 1.D-300 .LT. 1-X .LT. 0.2.
! DELTA   are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 0.9975 .LT. X .LE. 1-5.0D-16.
! FIRST   is a logical SAVE variable indicating whether it is necessary
!         to calculate the number of coefficients to use for each
!         Chebyshev expansion.
! FSIGN   is X or 1.0 - X.  It is used to remember the sign to be
!         assigned to the function value.
! I, J    are used as indices.
! IMIN    is the minimum index of a coefficient in the Chebyshev
!         polynomial expansion to be used.
! JIX     is an array containing MINXI, MAXXI, MINLAM, MAXLAM, MINDEL,
!         MAXDEL, MINMU, MAXMU in locations -1..6
! LAMDA   are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 0.8 .LT. X .LE. 0.9975.
! MU      are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 5.0D-16 .GT. 1-X .GE. 1.D-300.
! S2      is 2.0 * S.
! S       is the argument of the Chebyshev polynomial expansion.
! W1..W3  are adjacent elements of the recurrence used to evaluate the
!         Chebyshev polynomial expansion.
! XI      are coefficients of the Chebyshev polynomial expansion of
!         ERFC(X)/X for 0.0 .LE. X .LE. 0.8.
!
      DATA D /-1.548813042373261659512742D0, &
               2.565490123147816151928163D0, &
              -.5594576313298323225436913D0, &
               2.287915716263357638965891D0, &
              -9.199992358830151031278420D0, &
               2.794990820124599493768426D0/
!
      fsign = x
      arg = ABS(x)
      IF (arg.LT.0.0d0 .OR. arg.GE.1.0d0)THEN
         CALL derm1 ('DERFI',1,2,'Argument out of range','X',x,'.')
!     In case the error level is shifted to zero by the caller:
         derfi = 0.0d0
         RETURN
      END IF
      IF (arg.EQ.0.0d0) THEN
         derfi = 0.0d0
         RETURN
      END IF
      IF (arg.LE.0.8d0) THEN
         s = 3.125d0*arg*arg - 1.0d0
         j = -1
      ELSE
         IF (arg.LE.0.9975d0) THEN
            j = 1
         ELSE
            j = 3
         END IF
         arg = SQRT(-LOG((1.0d0-arg)*(1.0d0+arg)))
         s = d(j)*arg + d(j+1)
      END IF
      DERFI = SIGN(arg*DERFIX(s, j), fsign)
      RETURN
  END FUNCTION DERFI

  ! -------------------------------------------------------  DERFCI  -----
  DOUBLE PRECISION FUNCTION DERFCI(X)
!     Calculate the inverse of the complementary error function.

      DOUBLE PRECISION X
      DOUBLE PRECISION ARG, D(6), FSIGN, S
      INTEGER J
      DATA D /-1.548813042373261659512742D0, &
               2.565490123147816151928163D0, &
              -.5594576313298323225436913D0, &
               2.287915716263357638965891D0, &
              -9.199992358830151031278420D0, &
               2.794990820124599493768426D0/
!
!     Decide which approximation to use, and calculate the argument of
!     the Chebyshev polynomial expansion.
!
      IF (x.LE.0.0d0 .OR. x.GE.2.0d0) THEN
         CALL derm1('DERFCI',1,2,'Argument out of range','X',x,'.')
!     In case the error level is shifted to zero by the caller:
         derfci = 0.0d0
      END IF
      IF (x.EQ.1.0d0) THEN
         derfci = 0.0d0
         RETURN
      END IF
      fsign = 1.0d0 - x
      arg = ABS(fsign)
      IF (arg.LE.0.8d0) THEN
         s = 3.125d0*arg*arg - 1.0d0
         j = -1
      ELSE
         arg = 2.0d0 - x
         IF (x.LT.1.0d0) THEN
            s = x
         ELSE
            s = arg
         END IF
         arg = SQRT(-LOG(x*arg))
         IF (s.LT.5.0d-16) THEN
            j = 5
            s = d(5)/SQRT(arg) + d(6)
         ELSE
            IF (s.GE.0.0025d0) THEN
               j = 1
            ELSE IF (s.GE.5.0d-16) THEN
               j = 3
            END IF
            s = d(j)*arg + d(j+1)
         END IF
      END IF
      DERFCI = SIGN(arg*DERFIX(s, j), fsign)
      RETURN
  END FUNCTION DERFCI

  ! -------------------------------------------------------  DERFIX  -----
  DOUBLE PRECISION FUNCTION DERFIX(S, J)
!             Subroutine where most of calculations are done.

      INTEGER MAXDEL, MAXLAM, MAXMU, MAXXI, MINDEL, MINLAM
      INTEGER MINMU, MINXI, NC, NDELTA, NLAMDA, NMU, NXI
      PARAMETER (MINDEL = 0)
      PARAMETER (NDELTA = 37)
      PARAMETER (MAXDEL = MINDEL + NDELTA)
      PARAMETER (MINLAM = MAXDEL + 1)
      PARAMETER (NLAMDA = 26)
      PARAMETER (MAXLAM = MINLAM + NLAMDA)
      PARAMETER (MINMU = MAXLAM + 1)
      PARAMETER (NMU = 25)
      PARAMETER (MAXMU = MINMU + NMU)
      PARAMETER (MINXI = MAXMU + 1)
      PARAMETER (NXI = 38)
      PARAMETER (MAXXI = MINXI + NXI)
      PARAMETER (NC = MAXXI)
      DOUBLE PRECISION C(0:NC), DELTA(0:NDELTA)
      LOGICAL FIRST
      SAVE FIRST
      INTEGER I, J, JIX(-1:6)
      SAVE JIX
      INTEGER IMIN
      DOUBLE PRECISION LAMDA(0:NLAMDA), MU(0:NMU), S, S2
      DOUBLE PRECISION W1, W2, W3
      DOUBLE PRECISION XI(0:NXI)
!
!     *****     Equivalence Statements     *****************************
!
      EQUIVALENCE (C(MINDEL),DELTA(0))
      EQUIVALENCE (C(MINLAM),LAMDA(0))
      EQUIVALENCE (C(MINMU),MU(0))
      EQUIVALENCE (C(MINXI),XI(0))
!
!     *****     Data Statements     ************************************
!
!     DELTA(J), J = 0..NDELTA
!
!++ With first index 0, save data by elements if ~.C.
      DATA DELTA(0) /  .9566797090204925274526373D0 /
      DATA DELTA(1) / -.0231070043090649036999908D0 /
      DATA DELTA(2) / -.0043742360975084077333218D0 /
      DATA DELTA(3) / -.0005765034226511854809364D0 /
      DATA DELTA(4) / -.0000109610223070923931242D0 /
      DATA DELTA(5) /  .0000251085470246442787982D0 /
      DATA DELTA(6) /  .0000105623360679477511955D0 /
      DATA DELTA(7) /  .0000027544123300306391503D0 /
      DATA DELTA(8) /  .0000004324844983283380689D0 /
      DATA DELTA(9) /   -.0000000205303366552086916D0 /
      DATA DELTA(10) / -.0000000438915366654316784D0 /
      DATA DELTA(11) / -.0000000176840095080881795D0 /
      DATA DELTA(12) / -.0000000039912890280463420D0 /
      DATA DELTA(13) / -.0000000001869324124559212D0 /
      DATA DELTA(14) /  .0000000002729227396746077D0 /
      DATA DELTA(15) /  .0000000001328172131565497D0 /
      DATA DELTA(16) /  .0000000000318342484482286D0 /
      DATA DELTA(17) /  .0000000000016700607751926D0 /
      DATA DELTA(18) / -.0000000000020364649611537D0 /
      DATA DELTA(19) / -.0000000000009648468127965D0 /
      DATA DELTA(20) / -.0000000000002195672778128D0 /
      DATA DELTA(21) / -.0000000000000095689813014D0 /
      DATA DELTA(22) /  .0000000000000137032572230D0 /
      DATA DELTA(23) /  .0000000000000062538505417D0 /
      DATA DELTA(24) /  .0000000000000014584615266D0 /
      DATA DELTA(25) /  .0000000000000001078123993D0 /
      DATA DELTA(26) / -.0000000000000000709229988D0 /
      DATA DELTA(27) / -.0000000000000000391411775D0 /
      DATA DELTA(28) / -.0000000000000000111659209D0 /
      DATA DELTA(29) / -.0000000000000000015770366D0 /
      DATA DELTA(30) /  .0000000000000000002853149D0 /
      DATA DELTA(31) /  .0000000000000000002716662D0 /
      DATA DELTA(32) /  .0000000000000000000957770D0 /
      DATA DELTA(33) /  .0000000000000000000176835D0 /
      DATA DELTA(34) / -.0000000000000000000009828D0 /
      DATA DELTA(35) / -.0000000000000000000020464D0 /
      DATA DELTA(36) / -.0000000000000000000008020D0 /
      DATA DELTA(37) / -.0000000000000000000001650D0 /
!
!     LAMDA(J), J = 0..NLAMDA
!
!++ With first index 0, save data by elements if ~.C.
      DATA LAMDA(0) /  .9121588034175537733059200D0 /
      DATA LAMDA(1) / -.0162662818676636958546661D0 /
      DATA LAMDA(2) /  .0004335564729494453650589D0 /
      DATA LAMDA(3) /  .0002144385700744592065205D0 /
      DATA LAMDA(4) /  .0000026257510757648130176D0 /
      DATA LAMDA(5) / -.0000030210910501037969912D0 /
      DATA LAMDA(6) / -.0000000124060618367572157D0 /
      DATA LAMDA(7) /  .0000000624066092999917380D0 /
      DATA LAMDA(8) / -.0000000005401247900957858D0 /
      DATA LAMDA(9) /   -.0000000014232078975315910D0 /
      DATA LAMDA(10) /  .0000000000343840281955305D0 /
      DATA LAMDA(11) /  .0000000000335848703900138D0 /
      DATA LAMDA(12) / -.0000000000014584288516512D0 /
      DATA LAMDA(13) / -.0000000000008102174258833D0 /
      DATA LAMDA(14) /  .0000000000000525324085874D0 /
      DATA LAMDA(15) /  .0000000000000197115408612D0 /
      DATA LAMDA(16) / -.0000000000000017494333828D0 /
      DATA LAMDA(17) / -.0000000000000004800596619D0 /
      DATA LAMDA(18) /  .0000000000000000557302987D0 /
      DATA LAMDA(19) /  .0000000000000000116326054D0 /
      DATA LAMDA(20) / -.0000000000000000017262489D0 /
      DATA LAMDA(21) / -.0000000000000000002784973D0 /
      DATA LAMDA(22) /  .0000000000000000000524481D0 /
      DATA LAMDA(23) /  .0000000000000000000065270D0 /
      DATA LAMDA(24) / -.0000000000000000000015707D0 /
      DATA LAMDA(25) / -.0000000000000000000001475D0 /
      DATA LAMDA(26) /  .0000000000000000000000450D0 /
!
!     MU(J), J = 0..NMU
!
!++ With first index 0, save data by elements if ~.C.
      DATA MU(0) /  .9885750640661893136460358D0 /
      DATA MU(1) /  .0108577051845994776160281D0 /
      DATA MU(2) / -.0017511651027627952494825D0 /
      DATA MU(3) /  .0000211969932065633437984D0 /
      DATA MU(4) /  .0000156648714042435087911D0 /
      DATA MU(5) / -.0000005190416869103124261D0 /
      DATA MU(6) / -.0000000371357897426717780D0 /
      DATA MU(7) /  .0000000012174308662357429D0 /
      DATA MU(8) / -.0000000001768115526613442D0 /
      DATA MU(9) /   -.0000000000119372182556161D0 /
      DATA MU(10) /  .0000000000003802505358299D0 /
      DATA MU(11) / -.0000000000000660188322362D0 /
      DATA MU(12) / -.0000000000000087917055170D0 /
      DATA MU(13) / -.0000000000000003506869329D0 /
      DATA MU(14) / -.0000000000000000697221497D0 /
      DATA MU(15) / -.0000000000000000109567941D0 /
      DATA MU(16) / -.0000000000000000011536390D0 /
      DATA MU(17) / -.0000000000000000001326235D0 /
      DATA MU(18) / -.0000000000000000000263938D0 /
      DATA MU(19) /  .0000000000000000000005341D0 /
      DATA MU(20) / -.0000000000000000000022610D0 /
      DATA MU(21) /  .0000000000000000000009552D0 /
      DATA MU(22) / -.0000000000000000000005250D0 /
      DATA MU(23) /  .0000000000000000000002487D0 /
      DATA MU(24) / -.0000000000000000000001134D0 /
      DATA MU(25) /  .0000000000000000000000420D0 /
!
!     XI(J), J = 0..NXI
!
!++ With first index 0, save data by elements if ~.C.
      DATA XI(0) /  .9928853766189408231495800D0 /
      DATA XI(1) /  .1204675161431044864647846D0 /
      DATA XI(2) /  .0160781993420999447257039D0 /
      DATA XI(3) /  .0026867044371623158279591D0 /
      DATA XI(4) /  .0004996347302357262947170D0 /
      DATA XI(5) /  .0000988982185991204409911D0 /
      DATA XI(6) /  .0000203918127639944337340D0 /
      DATA XI(7) /  .0000043272716177354218758D0 /
      DATA XI(8) /  .0000009380814128593406758D0 /
      DATA XI(9) /  .0000002067347208683427411D0 /
      DATA XI(10) /  .0000000461596991054300078D0 /
      DATA XI(11) /  .0000000104166797027146217D0 /
      DATA XI(12) /  .0000000023715009995921222D0 /
      DATA XI(13) /  .0000000005439284068471390D0 /
      DATA XI(14) /  .0000000001255489864097987D0 /
      DATA XI(15) /  .0000000000291381803663201D0 /
      DATA XI(16) /  .0000000000067949421808797D0 /
      DATA XI(17) /  .0000000000015912343331469D0 /
      DATA XI(18) /  .0000000000003740250585245D0 /
      DATA XI(19) /  .0000000000000882087762421D0 /
      DATA XI(20) /  .0000000000000208650897725D0 /
      DATA XI(21) /  .0000000000000049488041039D0 /
      DATA XI(22) /  .0000000000000011766394740D0 /
      DATA XI(23) /  .0000000000000002803855725D0 /
      DATA XI(24) /  .0000000000000000669506638D0 /
      DATA XI(25) /  .0000000000000000160165495D0 /
      DATA XI(26) /  .0000000000000000038382583D0 /
      DATA XI(27) /  .0000000000000000009212851D0 /
      DATA XI(28) /  .0000000000000000002214615D0 /
      DATA XI(29) /  .0000000000000000000533091D0 /
      DATA XI(30) /  .0000000000000000000128488D0 /
      DATA XI(31) /  .0000000000000000000031006D0 /
      DATA XI(32) /  .0000000000000000000007491D0 /
      DATA XI(33) /  .0000000000000000000001812D0 /
      DATA XI(34) /  .0000000000000000000000439D0 /
      DATA XI(35) /  .0000000000000000000000106D0 /
      DATA XI(36) /  .0000000000000000000000026D0 /
      DATA XI(37) /  .0000000000000000000000006D0 /
      DATA XI(38) /  .0000000000000000000000002D0 /
!
      DATA FIRST /.TRUE./
!
      DATA JIX /MINXI, MAXXI, MINLAM, MAXLAM, MINDEL, MAXDEL, &
                MINMU, MAXMU/
!
!     *****     Procedures     *****************************************
!
!     Decide which approximation to use, and calculate the argument of
!     the Chebyshev polynomial expansion.
!
!
!     If this is the first call, calculate the degree of each expansion.
!
      IF (first) THEN
         first = .FALSE.
         s2 = 0.5d0*d1mach(3)
         DO 120 imin = -1, 5, 2
            DO 110 i = jix(imin), jix(imin+1)
               IF (ABS(c(i)).LT.s2) THEN
                  jix(imin+1) = i
                  go to 120
               END IF
110         CONTINUE
120      CONTINUE
      END IF
!
!     Evaluate the Chebyshev polynomial expansion.
!
      s2 = s + s
      w1 = 0.0d0
      w2 = 0.0d0
      imin = jix(j)
      i = jix(j+1)
200      w3 = w2
         w2 = w1
         w1 = (s2*w2 - w3) + c(i)
         i = i - 1
         IF (i.GT.imin) go to 200
      derfix = (s*w1 - w2) + c(imin)
      RETURN
  END FUNCTION DERFIX

  ! --------------------------------------------------------  SERFI  -----
  REAL             FUNCTION SERFI (X)
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.
!--S replaces "?": ?ERFI, ?ERFCI, ?ERFIX, ?ERM1
!>> 1996-06-18 SERFI Krogh  Changes to use .C. and C%%. J not changed.
!>> 1996-03-30 SERFI Krogh  Added external statements.
!>> 1995-11-28 SERFI Krogh  Removed multiple entries.
!>> 1995-11-03 SERFI Krogh  Removed blanks in numbers for C conversion.
!>> 1994-10-20 SERFI Krogh  Changes to use M77CON
!>> 1994-04-20 SERFI CLL Edited type stmts to make DP & SP files similar
!>> 1987-10-29 SERFI Snyder  Initial code.
!
!     For -1.0 .LT. X .LT. 1.0 calculate the inverse of the error
!     function.  That is, X = ERF(ERFI).
!
!     For 0.0 .LT. X .LT. 2.0 calculate the inverse of the
!     complementary error function.  that is, X = ERFC(ERFCI).  This
!     calculation is carried out by invoking the alternate entry *ERFCI.
!
!     If X is out of range, program execution is terminated by calling
!     the error message processor.
!
!     This subprogram uses approximations due to A. Strecok from
!     Mathematics of Computation 22, (1968) pp 144-158.
!

      REAL             X
      REAL             ARG, D(6), FSIGN, S
      INTEGER J
!
!     *****     Parameters     *****************************************
!
! MAX...  is the position in C of the last coefficient of a Chebyshev
!         polynomial expansion.
! MIN...  is the position in C of the first coefficient of a Chebyshev
!         polynomial expansion.
! NC      is the upper dimension of the array of coefficients.
! NDELTA  is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         0.9975 .LT. X .LE. 1-5.0E-16
! NLAMDA  is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         0.8 .LT. X .LE. 0.9975.
! NMU     is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate R(X) in the range
!         5.0E-16 .GT. 1-X .GE. 1.E-300.
! NXI     is the number of coefficients of the Chebyshev polynomial
!         expansion used to approximate SERFCI(X)/X in the
!         range 0.0 .LE. X .LE. 0.8.
!
!
!     *****     External References     ********************************
!
! R1MACH   Provides the round-off level.  Used to calculate the number
!          of coefficients to retain in each Chebyshev expansion.
! SERM1    Prints an error message and stops if X .LE. -1.0 or
!          X .GE. 1.0 (ERFI) or X .LE. 0.0 or X .GE. 2.0 (ERFCI).
! LOG      Calculates the natural logarithm.
! SQRT     Calculates the square root.
!
!
!     *****     Local Variables      ***********************************
!
! ARG     If ERFI or ERFCI is being approximated by a Chebyshev
!         expansion then ARG is the argument of ERFI or the argument
!         that would be used if ERFCI(X) were computed as ERFC(1-X),
!         that is, ARG = X if ERFI is being computed, or ARG = 1-X if
!         ERFCI is being computed.  If ERFI or ERFCI is being computed
!         using the initial approximation ERFI=SQRT(-LOG((1-X)*(1+X))),
!         then ARG is that initial approximation.
! C       contains the coefficients of polynomial expansions.  They are
!         stored in C in the order DELTA(0..37), LAMDA(0..26),
!         MU(0..25), XI(0..38).
! D       are used to scale the argument of the Chebyshev polynomial
!         expansion in the range 1.E-300 .LT. 1-X .LT. 0.2.
! DELTA   are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 0.9975 .LT. X .LE. 1-5.0E-16.
! FIRST   is a logical SAVE variable indicating whether it is necessary
!         to calculate the number of coefficients to use for each
!         Chebyshev expansion.
! FSIGN   is X or 1.0 - X.  It is used to remember the sign to be
!         assigned to the function value.
! I, J    are used as indices.
! IMIN    is the minimum index of a coefficient in the Chebyshev
!         polynomial expansion to be used.
! JIX     is an array containing MINXI, MAXXI, MINLAM, MAXLAM, MINDEL,
!         MAXDEL, MINMU, MAXMU in locations -1..6
! LAMDA   are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 0.8 .LT. X .LE. 0.9975.
! MU      are coefficients of the Chebyshev polynomial expansion of R(X)
!         for 5.0E-16 .GT. 1-X .GE. 1.E-300.
! S2      is 2.0 * S.
! S       is the argument of the Chebyshev polynomial expansion.
! W1..W3  are adjacent elements of the recurrence used to evaluate the
!         Chebyshev polynomial expansion.
! XI      are coefficients of the Chebyshev polynomial expansion of
!         ERFC(X)/X for 0.0 .LE. X .LE. 0.8.
!
      DATA D /-1.548813042373261659512742E0, &
               2.565490123147816151928163E0, &
              -.5594576313298323225436913E0, &
               2.287915716263357638965891E0, &
              -9.199992358830151031278420E0, &
               2.794990820124599493768426E0/
!
      fsign = x
      arg = ABS(x)
      IF (arg.LT.0.0e0 .OR. arg.GE.1.0e0)THEN
         CALL serm1 ('SERFI',1,2,'Argument out of range','X',x,'.')
!     In case the error level is shifted to zero by the caller:
         serfi = 0.0e0
         RETURN
      END IF
      IF (arg.EQ.0.0e0) THEN
         serfi = 0.0e0
         RETURN
      END IF
      IF (arg.LE.0.8e0) THEN
         s = 3.125e0*arg*arg - 1.0e0
         j = -1
      ELSE
         IF (arg.LE.0.9975e0) THEN
            j = 1
         ELSE
            j = 3
         END IF
         arg = SQRT(-LOG((1.0e0-arg)*(1.0e0+arg)))
         s = d(j)*arg + d(j+1)
      END IF
      SERFI = SIGN(arg*SERFIX(s, j), fsign)
      RETURN
    END FUNCTION SERFI

  ! -------------------------------------------------------  SERFCI  -----
  REAL             FUNCTION SERFCI(X)
!     Calculate the inverse of the complementary error function.
!

      REAL             X
      REAL             ARG, D(6), FSIGN, S
      INTEGER J
      DATA D /-1.548813042373261659512742E0, &
               2.565490123147816151928163E0, &
              -.5594576313298323225436913E0, &
               2.287915716263357638965891E0, &
              -9.199992358830151031278420E0, &
               2.794990820124599493768426E0/
!
!     Decide which approximation to use, and calculate the argument of
!     the Chebyshev polynomial expansion.
!
      IF (x.LE.0.0e0 .OR. x.GE.2.0e0) THEN
         CALL serm1('SERFCI',1,2,'Argument out of range','X',x,'.')
!     In case the error level is shifted to zero by the caller:
         serfci = 0.0e0
      END IF
      IF (x.EQ.1.0e0) THEN
         serfci = 0.0e0
         RETURN
      END IF
      fsign = 1.0e0 - x
      arg = ABS(fsign)
      IF (arg.LE.0.8e0) THEN
         s = 3.125e0*arg*arg - 1.0e0
         j = -1
      ELSE
         arg = 2.0e0 - x
         IF (x.LT.1.0e0) THEN
            s = x
         ELSE
            s = arg
         END IF
         arg = SQRT(-LOG(x*arg))
         IF (s.LT.5.0e-16) THEN
            j = 5
            s = d(5)/SQRT(arg) + d(6)
         ELSE
            IF (s.GE.0.0025e0) THEN
               j = 1
            ELSE IF (s.GE.5.0e-16) THEN
               j = 3
            END IF
            s = d(j)*arg + d(j+1)
         END IF
      END IF
      SERFCI = SIGN(arg*SERFIX(s, j), fsign)
      RETURN
  END FUNCTION SERFCI

  ! -------------------------------------------------------  SERFIX  -----
  REAL             FUNCTION SERFIX(S, J)
!             Subroutine where most of calculations are done.

      INTEGER MAXDEL, MAXLAM, MAXMU, MAXXI, MINDEL, MINLAM
      INTEGER MINMU, MINXI, NC, NDELTA, NLAMDA, NMU, NXI
      PARAMETER (MINDEL = 0)
      PARAMETER (NDELTA = 37)
      PARAMETER (MAXDEL = MINDEL + NDELTA)
      PARAMETER (MINLAM = MAXDEL + 1)
      PARAMETER (NLAMDA = 26)
      PARAMETER (MAXLAM = MINLAM + NLAMDA)
      PARAMETER (MINMU = MAXLAM + 1)
      PARAMETER (NMU = 25)
      PARAMETER (MAXMU = MINMU + NMU)
      PARAMETER (MINXI = MAXMU + 1)
      PARAMETER (NXI = 38)
      PARAMETER (MAXXI = MINXI + NXI)
      PARAMETER (NC = MAXXI)

      REAL             C(0:NC), DELTA(0:NDELTA)
      LOGICAL FIRST
      SAVE FIRST
      INTEGER I, J, JIX(-1:6)
      SAVE JIX
      INTEGER IMIN
      REAL             LAMDA(0:NLAMDA), MU(0:NMU), S, S2
      REAL             W1, W2, W3
      REAL             XI(0:NXI)
!
!     *****     Equivalence Statements     *****************************
!
      EQUIVALENCE (C(MINDEL),DELTA(0))
      EQUIVALENCE (C(MINLAM),LAMDA(0))
      EQUIVALENCE (C(MINMU),MU(0))
      EQUIVALENCE (C(MINXI),XI(0))
!
!     *****     Data Statements     ************************************
!
!     DELTA(J), J = 0..NDELTA
!
!++ With first index 0, save data by elements if ~.C.
      DATA DELTA(0) /  .9566797090204925274526373E0 /
      DATA DELTA(1) / -.0231070043090649036999908E0 /
      DATA DELTA(2) / -.0043742360975084077333218E0 /
      DATA DELTA(3) / -.0005765034226511854809364E0 /
      DATA DELTA(4) / -.0000109610223070923931242E0 /
      DATA DELTA(5) /  .0000251085470246442787982E0 /
      DATA DELTA(6) /  .0000105623360679477511955E0 /
      DATA DELTA(7) /  .0000027544123300306391503E0 /
      DATA DELTA(8) /  .0000004324844983283380689E0 /
      DATA DELTA(9) /   -.0000000205303366552086916E0 /
      DATA DELTA(10) / -.0000000438915366654316784E0 /
      DATA DELTA(11) / -.0000000176840095080881795E0 /
      DATA DELTA(12) / -.0000000039912890280463420E0 /
      DATA DELTA(13) / -.0000000001869324124559212E0 /
      DATA DELTA(14) /  .0000000002729227396746077E0 /
      DATA DELTA(15) /  .0000000001328172131565497E0 /
      DATA DELTA(16) /  .0000000000318342484482286E0 /
      DATA DELTA(17) /  .0000000000016700607751926E0 /
      DATA DELTA(18) / -.0000000000020364649611537E0 /
      DATA DELTA(19) / -.0000000000009648468127965E0 /
      DATA DELTA(20) / -.0000000000002195672778128E0 /
      DATA DELTA(21) / -.0000000000000095689813014E0 /
      DATA DELTA(22) /  .0000000000000137032572230E0 /
      DATA DELTA(23) /  .0000000000000062538505417E0 /
      DATA DELTA(24) /  .0000000000000014584615266E0 /
      DATA DELTA(25) /  .0000000000000001078123993E0 /
      DATA DELTA(26) / -.0000000000000000709229988E0 /
      DATA DELTA(27) / -.0000000000000000391411775E0 /
      DATA DELTA(28) / -.0000000000000000111659209E0 /
      DATA DELTA(29) / -.0000000000000000015770366E0 /
      DATA DELTA(30) /  .0000000000000000002853149E0 /
      DATA DELTA(31) /  .0000000000000000002716662E0 /
      DATA DELTA(32) /  .0000000000000000000957770E0 /
      DATA DELTA(33) /  .0000000000000000000176835E0 /
      DATA DELTA(34) / -.0000000000000000000009828E0 /
      DATA DELTA(35) / -.0000000000000000000020464E0 /
      DATA DELTA(36) / -.0000000000000000000008020E0 /
      DATA DELTA(37) / -.0000000000000000000001650E0 /
!
!     LAMDA(J), J = 0..NLAMDA
!
!++ With first index 0, save data by elements if ~.C.
      DATA LAMDA(0) /  .9121588034175537733059200E0 /
      DATA LAMDA(1) / -.0162662818676636958546661E0 /
      DATA LAMDA(2) /  .0004335564729494453650589E0 /
      DATA LAMDA(3) /  .0002144385700744592065205E0 /
      DATA LAMDA(4) /  .0000026257510757648130176E0 /
      DATA LAMDA(5) / -.0000030210910501037969912E0 /
      DATA LAMDA(6) / -.0000000124060618367572157E0 /
      DATA LAMDA(7) /  .0000000624066092999917380E0 /
      DATA LAMDA(8) / -.0000000005401247900957858E0 /
      DATA LAMDA(9) /   -.0000000014232078975315910E0 /
      DATA LAMDA(10) /  .0000000000343840281955305E0 /
      DATA LAMDA(11) /  .0000000000335848703900138E0 /
      DATA LAMDA(12) / -.0000000000014584288516512E0 /
      DATA LAMDA(13) / -.0000000000008102174258833E0 /
      DATA LAMDA(14) /  .0000000000000525324085874E0 /
      DATA LAMDA(15) /  .0000000000000197115408612E0 /
      DATA LAMDA(16) / -.0000000000000017494333828E0 /
      DATA LAMDA(17) / -.0000000000000004800596619E0 /
      DATA LAMDA(18) /  .0000000000000000557302987E0 /
      DATA LAMDA(19) /  .0000000000000000116326054E0 /
      DATA LAMDA(20) / -.0000000000000000017262489E0 /
      DATA LAMDA(21) / -.0000000000000000002784973E0 /
      DATA LAMDA(22) /  .0000000000000000000524481E0 /
      DATA LAMDA(23) /  .0000000000000000000065270E0 /
      DATA LAMDA(24) / -.0000000000000000000015707E0 /
      DATA LAMDA(25) / -.0000000000000000000001475E0 /
      DATA LAMDA(26) /  .0000000000000000000000450E0 /
!
!     MU(J), J = 0..NMU
!
!++ With first index 0, save data by elements if ~.C.
      DATA MU(0) /  .9885750640661893136460358E0 /
      DATA MU(1) /  .0108577051845994776160281E0 /
      DATA MU(2) / -.0017511651027627952494825E0 /
      DATA MU(3) /  .0000211969932065633437984E0 /
      DATA MU(4) /  .0000156648714042435087911E0 /
      DATA MU(5) / -.0000005190416869103124261E0 /
      DATA MU(6) / -.0000000371357897426717780E0 /
      DATA MU(7) /  .0000000012174308662357429E0 /
      DATA MU(8) / -.0000000001768115526613442E0 /
      DATA MU(9) /   -.0000000000119372182556161E0 /
      DATA MU(10) /  .0000000000003802505358299E0 /
      DATA MU(11) / -.0000000000000660188322362E0 /
      DATA MU(12) / -.0000000000000087917055170E0 /
      DATA MU(13) / -.0000000000000003506869329E0 /
      DATA MU(14) / -.0000000000000000697221497E0 /
      DATA MU(15) / -.0000000000000000109567941E0 /
      DATA MU(16) / -.0000000000000000011536390E0 /
      DATA MU(17) / -.0000000000000000001326235E0 /
      DATA MU(18) / -.0000000000000000000263938E0 /
      DATA MU(19) /  .0000000000000000000005341E0 /
      DATA MU(20) / -.0000000000000000000022610E0 /
      DATA MU(21) /  .0000000000000000000009552E0 /
      DATA MU(22) / -.0000000000000000000005250E0 /
      DATA MU(23) /  .0000000000000000000002487E0 /
      DATA MU(24) / -.0000000000000000000001134E0 /
      DATA MU(25) /  .0000000000000000000000420E0 /
!
!     XI(J), J = 0..NXI
!
!++ With first index 0, save data by elements if ~.C.
      DATA XI(0) /  .9928853766189408231495800E0 /
      DATA XI(1) /  .1204675161431044864647846E0 /
      DATA XI(2) /  .0160781993420999447257039E0 /
      DATA XI(3) /  .0026867044371623158279591E0 /
      DATA XI(4) /  .0004996347302357262947170E0 /
      DATA XI(5) /  .0000988982185991204409911E0 /
      DATA XI(6) /  .0000203918127639944337340E0 /
      DATA XI(7) /  .0000043272716177354218758E0 /
      DATA XI(8) /  .0000009380814128593406758E0 /
      DATA XI(9) /  .0000002067347208683427411E0 /
      DATA XI(10) /  .0000000461596991054300078E0 /
      DATA XI(11) /  .0000000104166797027146217E0 /
      DATA XI(12) /  .0000000023715009995921222E0 /
      DATA XI(13) /  .0000000005439284068471390E0 /
      DATA XI(14) /  .0000000001255489864097987E0 /
      DATA XI(15) /  .0000000000291381803663201E0 /
      DATA XI(16) /  .0000000000067949421808797E0 /
      DATA XI(17) /  .0000000000015912343331469E0 /
      DATA XI(18) /  .0000000000003740250585245E0 /
      DATA XI(19) /  .0000000000000882087762421E0 /
      DATA XI(20) /  .0000000000000208650897725E0 /
      DATA XI(21) /  .0000000000000049488041039E0 /
      DATA XI(22) /  .0000000000000011766394740E0 /
      DATA XI(23) /  .0000000000000002803855725E0 /
      DATA XI(24) /  .0000000000000000669506638E0 /
      DATA XI(25) /  .0000000000000000160165495E0 /
      DATA XI(26) /  .0000000000000000038382583E0 /
      DATA XI(27) /  .0000000000000000009212851E0 /
      DATA XI(28) /  .0000000000000000002214615E0 /
      DATA XI(29) /  .0000000000000000000533091E0 /
      DATA XI(30) /  .0000000000000000000128488E0 /
      DATA XI(31) /  .0000000000000000000031006E0 /
      DATA XI(32) /  .0000000000000000000007491E0 /
      DATA XI(33) /  .0000000000000000000001812E0 /
      DATA XI(34) /  .0000000000000000000000439E0 /
      DATA XI(35) /  .0000000000000000000000106E0 /
      DATA XI(36) /  .0000000000000000000000026E0 /
      DATA XI(37) /  .0000000000000000000000006E0 /
      DATA XI(38) /  .0000000000000000000000002E0 /
!
      DATA FIRST /.TRUE./
!
      DATA JIX /MINXI, MAXXI, MINLAM, MAXLAM, MINDEL, MAXDEL, &
                MINMU, MAXMU/
!
!     *****     Procedures     *****************************************
!
!     Decide which approximation to use, and calculate the argument of
!     the Chebyshev polynomial expansion.
!
!
!     If this is the first call, calculate the degree of each expansion.
!
      IF (first) THEN
         first = .FALSE.
         s2 = 0.5e0*r1mach(3)
         DO 120 imin = -1, 5, 2
            DO 110 i = jix(imin), jix(imin+1)
               IF (ABS(c(i)).LT.s2) THEN
                  jix(imin+1) = i
                  go to 120
               END IF
110         CONTINUE
120      CONTINUE
      END IF
!
!     Evaluate the Chebyshev polynomial expansion.
!
      s2 = s + s
      w1 = 0.0e0
      w2 = 0.0e0
      imin = jix(j)
      i = jix(j+1)
200      w3 = w2
         w2 = w1
         w1 = (s2*w2 - w3) + c(i)
         i = i - 1
         IF (i.GT.imin) go to 200
      serfix = (s*w1 - w2) + c(imin)
      RETURN
  END FUNCTION SERFIX

!=============================================================================
END MODULE MathUtils
!=============================================================================

! $Log$
! Revision 2.1  2002/03/29 20:22:37  perun
! Version 1.0 commit
!
