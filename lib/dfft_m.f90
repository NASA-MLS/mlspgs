! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DFFT_M

  implicit NONE
  private
  public :: DRFT1, DFFT

  !---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: IdParm = &
    & "$Id$"
  character(len=len(idparm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
       & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

!=====================================================================
!=================== Coversion of the f77 fft.f to f90 ===============
!  The JPL Double Precision FFT Package.
!---------------------------------------------------------------------

      SUBROUTINE DRFT1 (A,MODE,M,MS,S)

!>> 1994-11-11 DRFT1  Krogh   Declared all vars.
!>> 1994-10-20 DRFT1 Krogh  Changes to use M77CON
!>> 1989-05-07 DRFT1 FTK & CLL
!>> 1989-04-21 FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.

!     This subroutine computes Fourier transforms of real data using
!     the Cooley-Tukey fast Fourier transform.

!     Variables in the calling sequence have the following types
      DOUBLE PRECISION A(*), S(*)
      INTEGER          M, MS
      CHARACTER        MODE

!     Programmed by Fred T. Krogh at the Jet Propulsion Laboratory,
!     Pasadena, Calif.   August 1, 1969.
!     Revised for portability by Krogh -- January 22, 1988

!     In describing the usage the following notation is used
!     N = 2 ** M
!     W = EXP(2*PI*I/N), where I = SQRT(-1) and PI = 3.14159...

!     The usage is as follows

! A() is an array of function values if one is doing Fourier analysis,
!  and is an array of Fourier coefficients if one is doing Fourier
!  synthesis.  In our description here we assume that A is a real
!  array with dimension A(N) when A contains the real data, X, and
!  that A is a complex array with dimension A(N/2) when A contains
!  complex Fourier coefficients, C.  (C(k) for k > N/2 need not be
!  saved since for 0 < k < N, C(N-k) = conjugate of C(k).  It is
!  assumed that the imaginary part of a complex number is stored
!  in the cell immediately following its real part, except that
!  A(1) = C(0), and A(2) = C(N/2).  This is possible since these
!  values of C are real and by doing this both X and C require the
!  same storage in A. Of course the storage required for A can be
!  reserved by the user in any way that works.

! MODE  Selects Synthesis or Analysis.
!  If MODE = 'A' or 'a', do Fourier analysis, which amounts to setting
!  C(k) = sum for j=0 to N-1 of X(j)*T(M,j,k), for k = 0, N/2
!  with  T(M,j,k) = (1/N) * W ** (-j*k).
!  If MODE = 'S' or 's', do Fourier synthesis, which amounts to setting
!  X(j) = sum for k=0 to N-1 of C(k)*T(M,j,k), for j = 0, N - 1
!  with  T(M,j,k) = W ** (j*k)
!  (Recall that C(N-k) is the conjugate of C(k).)

! M is used to indicate N = 2**M, the number of real points.  The
!  number of points must satisfy 1 .le. N .le. 2**21.
!  M = 0 gives an immediate return.

! MS gives the state of the sine table.  If MS > 0, there are NT =
!    2 ** (MS-2) good entries in the sine table.  On the initial call,
!    MS must be set to 0, and when the return is made, it will be set
!    to M, which is the value of MS required for computing a
!    transform of size N.  If MS = -1, the sine table will be computed
!    as for the case MS = 0, and THEN a return to the user will be made
!    with MS set as before, but no transform will be computed.  This
!    option is useful if the user would like access to the sine table
!    before computing the FFT.
!    On detected errors the error message subrs are called and
!    execution stops.  If the user overrides the stop to cause
!    continuation, THEN this subr will return with MS = -2.

! S() is a vector, S(j) = sin(pi*j/2*NT)), j = 1, 2, ..., NT-1, where
!  NT is defined in the description of MS above.  S is computed by the
!  subroutine if M .gt. MS.  (If S is altered, set MS=0 so that S
!  is recomputed.)

!     ------------------------------------------------------------------
!                Notes on COMMON, PARAMETER's, and local variables

!     MMAX is the largest value allowed for M
!     The dimension of KE must be at least as large as MMAX-1.
!     The named common CDFFTC is used for communication between this
!     subroutine and the subroutine DFFT which computes a one
!     dimensional complex Fourier transform and computes the sine table.
!     The use of the variables in CDFFTC is contained in the listing
!     of DFFT.

!     ANAL = .TRUE. when doing Fourier analysis, and .false. otherwise.

!     N1 = 2 ** M
!     ------------------------------------------------------------------
!--D replaces "?": ?RFT1, ?FFT, C?FFTC
!     Both versions use ERMSG, IERM1
!     and need ERFIN, IERV1
!     ------------------------------------------------------------------
      INTEGER MMAX
      INTEGER I, II, II1, II2, IR, IR1, IR2
      INTEGER J, JDIF, JJ
      INTEGER K1, K1N, KN2
      INTEGER L
      INTEGER MA, MSI
      INTEGER N1, N1P, KEDIM
 
      DOUBLE PRECISION FN, HALF
      DOUBLE PRECISION T, TI, TT, TTI, TWO, WI, WR
 
      LOGICAL ANAL
 
      PARAMETER (TWO = 2.0D0)
      PARAMETER (HALF = 0.5D0)
      EQUIVALENCE (ILAST, N1)
! Common variables
      PARAMETER (KEDIM=20)
      LOGICAL NEEDST
      INTEGER MT, NT, MM, KS, ILAST, KE(KEDIM), KEE(KEDIM+1)
! Note that KEE(1) is equivalent to ILAST.
      EQUIVALENCE (KE(1), KEE(2))
      COMMON /CDFFTC/ NEEDST, MT, NT, MM, KS, ILAST, KE
      SAVE /CDFFTC/
      PARAMETER (MMAX = KEDIM+1)
!     ------------------------------------------------------------------

      if ( MODE .eq. 'A' .or. MODE .eq. 'a' ) then
         ANAL = .true.
      else if ( MODE .eq. 'S' .or. MODE .eq. 's' ) then
         ANAL = .false.
      else
         Print *,'** Fatal error in DRFT1 **'
         Print *,'   Bad MODE.  MODE = ',MODE
         MS = -2
         return
      end if
      MA = M
      MSI = MS
      NEEDST = MA .GT. MSI
      if ( NEEDST ) then
         if ( MA .GT. MMAX .or. MA .lt. 0 ) then
!                               Fatal error, default is to stop in IERM1
            Print *,'** Fatal error in DRFT1 **'
            Print *,'   Require (0 .le. M .le. 21), M=',M
            MS = -2
            return
         end if
         MS = MA
         MT = MA - 2
         call DFFT (A, A, S)
         if ( MSI .EQ. -1) return
      end if
      if ( MA .NE. 0 ) then
         MM = MA - 1
         N1 = 2 ** MA
         N1P = N1 + 2
         KN2 = N1 / 2
         JDIF = (4 * NT) / N1
         KS = 2
         if ( ANAL ) then
!                               Set flags for Fourier analysis
            IR = 2
            II = 1
            FN = HALF / DBLE(N1)
!           Doing Fourier analysis, so multiply by 2 ** M
            do 10 I = 1, N1
               A(I) = A(I) * FN
   10       continue
         else
!                              Set flags for Fourier synthesis
            IR = 1
            II = 2
            go to 40
         end if
 
!                              Compute complex Fourier transform
   20    do 30 L = 1, MM
            KEE(L+1) = KEE(L) / 2
   30    continue

         call DFFT (A(IR), A(II), S)
!                              End of computing complex transform

         if ( .NOT. ANAL) return

!        Beginning of calculations relating coefficients of real data
!        with coefficients of associated complex data

!        Special case --  K1 = 0
   40    T = A(1) + A(2)
         TI = A(1) - A(2)
         if ( ANAL ) then
            T = TWO * T
            TI = TWO * TI
         end if
         A(1) = T
         A(2) = TI
         if ( MM .GT. 0 ) then
!                           Special kase -- K1 = N1 / 4
            A(KN2+1) = TWO * A(KN2+1)
            A(KN2+2) = -TWO * A(KN2+2)
            if ( MM .GT. 1 ) then
               J = 0
               do 50 K1 = 3, KN2, 2
                  K1N = N1P - K1
                  if ( ANAL ) then
                     IR1 = K1N
                     IR2 = K1
                  else
                     IR1 = K1
                     IR2 = K1N
                  end if
                  II2 = IR2 + 1
                  II1 = IR1 + 1
                  J = J + JDIF
                  WI = S(J)
                  JJ = NT - J
                  WR = S(JJ)
                  T = A(IR1) - A(IR2)
                  TI = A(II1) + A(II2)
                  TT = T * WI + TI * WR
                  TTI = T * WR - TI * WI
                  T = A(IR1) + A(IR2)
                  TI = A(II1) - A(II2)
                  A(IR1) = T - TT
                  A(IR2) = T + TT
                  A(II1) = TTI + TI
                  A(II2) = TTI - TI
   50          continue
            end if
         end if
         if ( .NOT. ANAL) go to 20
      end if

      return

      END SUBROUTINE DRFT1

!==========================================================================

      SUBROUTINE DFFT(AR, AI, S)

!>> 1994-11-11 DFFT  Krogh   Declared all vars.
!>> 1994-10-20 DFFT  Krogh  Changes to use M77CON
!>> 1989-08-15 DFFT  FTK -- Parameter constants given to more precision.
!>> 1989-04-21 DFFT  FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  All rights reserved.  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged..

!     This subroutine computes one dimensional complex Fourier
!     transforms using the Cooley-Tukey algorithm. It is used by a
!     number of subroutines which compute Fourier transforms.

!     Programmed by Fred T. Krogh at the Jet Propulsion Laboratory,
!     Pasadena, Calif.   August 1, 1969.
!     Revised by Krogh at JPL -- January 19, 1988 -- For portability

      DOUBLE PRECISION AR(*), AI(*), S(*)
!     Minimum dimensions are AR(ILAST), AI(ILAST), S(NT-1), where ILAST
!     and NT are defined in the common block below.

! AR and AI give the arrays used to hold the real and imaginary parts of
!     the Fourier coefficients and data.
! S   contains the sine table required in the calculations.

!     -----------------------------------------------------------------
!                Notes on COMMON, PARAMETER's, and local variables

!     SPI4 = SIN(PI/4) = 0.5 * SQRT(2)
!     PI4 = PI / 4
!     THE DIMENSION OF KE MUST BE AS LARGE AS THE MAXIMUM VALUE
!     PERMITTED FOR MM.

!     WHERE
!     NEEDST= .TRUE. if the sine table must be computed.
!     MT    = base 2 log(NT)
!     NT    = number of entries in the sine table
!     MM    = base 2 log(number of complex fourier coefficients)
!     KS    = distance in memory between successive points.  The i-th
!             coefficient, a(i) is given by AR((I+1)*KS)+AI((I+1)*KS)*
!             sqrt(-1), i=0, 1, ..., (2**MM)-1.
!     ILAST = KS * 2**MM
!     KE(L) = KS * 2**(MM-L), L=1, 2, ..., MM
!     -----------------------------------------------------------------
!--D replaces "?": ?FFT, C?FFTC
!     -----------------------------------------------------------------
      INTEGER I, I1, I2, I3, IJ
      INTEGER J, JDIF, JGO, JI, JI2, JJ, JR
      INTEGER K, KSI
      INTEGER L, L1, L4, LJ, LL
      INTEGER MTC
 
      DOUBLE PRECISION HALF, PI4, SPI4
      DOUBLE PRECISION T, T1, T1I, T2, T2I, T3, T3I, THETA
      DOUBLE PRECISION TI, TP, TP1, TP1I, TPI
      DOUBLE PRECISION WI1, WI2, WI3, WR1, WR2, WR3
 
      PARAMETER (SPI4 = 0.7071067811865475244008443621048490D0)
      PARAMETER (PI4 = 0.7853981633974483096156608458198757D0)
      PARAMETER (HALF = 0.5D0)
 
! Common variables
      INTEGER KEDIM
      PARAMETER (KEDIM=20)
      LOGICAL NEEDST
      INTEGER MT, NT, MM, KS, ILAST, KE(KEDIM), KEE(KEDIM+1)
! Note that KEE(1) is equivalent to ILAST.
      EQUIVALENCE (KE(1), KEE(2))
      COMMON /CDFFTC/ NEEDST, MT, NT, MM, KS, ILAST, KE
      SAVE /CDFFTC/
!     -----------------------------------------------------------------

      if ( NEEDST ) then
!                      Compute the sine table
         NEEDST = .FALSE.
         MTC = MT
         NT = 2**MTC
         if ( MTC .GT. 0 ) then
!                            SET FOR L=1
            J = NT
            JJ = J / 2
            S(JJ) = SPI4
            if ( MTC .GT. 1 ) then
               THETA = PI4
               do 20 L = 2, MTC
                  THETA = HALF * THETA
                  K = J
                  J = JJ
                  JJ = JJ / 2
!                       At this point J = 2**(MT-L+1) and JJ = 2**(MT-L)
                  S(JJ) = SIN(THETA)
                  L1 = NT - JJ
                  S(L1) = COS(THETA)
                  LL = NT - K
                  if ( LL .GE. J ) then
                     do 10 I = J, LL, J
                        I1 = NT - I
                        I2 = I + JJ
                        S(I2) = S(I) * S(L1) + S(JJ) * S(I1)
   10                continue
                  end if
   20          continue
            end if
         end if
      else
!                      Compute the transform
!                           Scramble A

         IJ = 1
         JI = 1
         L1 = KS
         if ( MM .GT. 1 ) then
   30       IJ = IJ + L1
            do 40 I = 1, MM
               JI = JI + KE(I)
               KE(I) = -KE(I)
               if ( KE(I) .LT. 0 ) then
                  if ( IJ .LT. JI ) then
!                       INTERCHANGE THE IJ-TH COEFFICIENT WITH THE JI-TH
                     T = AR(IJ)
                     AR(IJ) = AR(JI)
                     AR(JI) = T
                     T = AI(IJ)
                     AI(IJ) = AI(JI)
                     AI(JI) = T
                  end if
                  go to 30
               end if
   40       continue
         end if
!                          END OF SCRAMBLING A
         JDIF = NT
         if ( MOD(MM, 2) .NE. 0 ) then
!                    SPECIAL CASE -- L = 1,  MM ODD  (RADIX 2 ALGORITHM)
            L1 = L1 + L1
            do 50 I = 1, ILAST, L1
               KSI = I + KS
               T = AR(I)
               AR(I) = T + AR(KSI)
               AR(KSI) = T - AR(KSI)
               T = AI(I)
               AI(I) = T + AI(KSI)
               AI(KSI) = T - AI(KSI)
   50       continue
            JDIF = JDIF / 2
         end if

         do 140 L = 2, MM, 2
            L4 = 4 * L1
            LJ = L1 / 2
            J = 0
            JI = 0

!           ASSIGN 70 TO JGO
            JGO = 70

!           START OF I LOOP  (RADIX 4 ALGORITHM)
   60       IJ = J + 1
            do 120 I = IJ, ILAST, L4
               I1 = I + L1
               I2 = I1 + L1
               I3 = I2 + L1

!              go to JGO, (70, 80, 90)

               if ( JGO.eq.70 ) then
                 go to 70
               else if ( JGO.eq.80 ) then
                 go to 80
               else if ( JGO.eq.90 ) then
                 go to 90
               end if

!              SPECIAL CASE -- J = 0
   70          T = AR(I) + AR(I1)
               T1 = AR(I) - AR(I1)
               TI = AI(I) + AI(I1)
               T1I = AI(I) - AI(I1)
               T2 = AR(I2) + AR(I3)
               T3 = AR(I2) - AR(I3)
               T2I = AI(I2) + AI(I3)
               T3I = AI(I2) - AI(I3)
               go to 110

!              SPECIAL CASE -- J = L1 / 2
   80          T2 = SPI4 * AR(I2)
               T2I = SPI4 * AI(I2)
               T3 = SPI4 * AR(I3)
               T3I = SPI4 * AI(I3)
               TP = T2 - T2I
               TPI = T2 + T2I
               TP1 = -T3 - T3I
               TP1I = T3 - T3I
               T = AR(I) - AI(I1)
               T1 = AR(I) + AI(I1)
               TI = AI(I) + AR(I1)
               T1I = AI(I) - AR(I1)
               go to 100

!              USUAL CASE -- J .NE. 0  AND  J .NE. L1 / 2

!              WRK AND WIK (K = 1, 2, 3) ARE THE REAL AND IMAGINARY PART
!              RESP. OF EXP(SQRT(-1) * PI * K*(2**(-L-MOD(MM, 2)))*J/KS)
!                             =EXP(SQRT(-1) * PI * K * (J / (4 * L1)))

   90          T2 = WR2 * AR(I1) - WI2 * AI(I1)
               T2I = WI2 * AR(I1) + WR2 * AI(I1)
               T = AR(I) + T2
               T1 = AR(I) - T2
               TI = AI(I) + T2I
               T1I = AI(I) - T2I
               TP = WR1 * AR(I2) - WI1 * AI(I2)
               TPI = WI1 * AR(I2) + WR1 * AI(I2)
               TP1 = WR3 * AR(I3) - WI3 * AI(I3)
               TP1I = WI3 * AR(I3) + WR3 * AI(I3)
  100          T2 = TP + TP1
               T3 = TP - TP1
               T2I = TPI + TP1I
               T3I = TPI - TP1I
  110          AR(I) = T + T2
               AI(I) = TI + T2I
               AR(I1) = T1 - T3I
               AI(I1) = T1I + T3
               AR(I2) = T - T2
               AI(I2) = TI - T2I
               AR(I3) = T1 + T3I
               AI(I3) = T1I - T3
  120       continue
!           END OF I LOOP

            if ( J .LT. LJ ) then
               if ( J .EQ. 0 ) then
!                 ASSIGN 90 TO JGO
                  JGO = 90
                  J = KS
               else
                  J = L1 - J
!                 COMPUTE WR-S AND WI-S FOR J REPLACED BY L1 - J
                  T = WI1
                  WI1 = WR1
                  WR1 = T
                  WR2 = -WR2
                  T = -WI3
                  WI3 = -WR3
                  WR3 = T
                  go to 60
               end if
            else if ( J .EQ. LJ ) then
               go to 130
            else
               J = L1 - J + KS
            end if
            if ( J .LT. LJ ) then
!                              COMPUTE WR-S AND WI-S
               JI = JI + JDIF
               JR = NT - JI
               WR1 = S(JR)
               WI1 = S(JI)
               JI2 = JI + JI
               WI2 = S(JI2)
               JR = NT - JI2
               WR2 = S(JR)
               JI2 = JI + JI2
               if ( JI2 .LE. NT ) then
                  JR = NT - JI2
                  WR3 = S(JR)
                  WI3 = S(JI2)
                  go to 60
               end if
               JR = JI2 - NT
               JI2 = NT - JR
               WI3 = S(JI2)
               WR3 = -S(JR)
               go to 60
            else if ( J .EQ. LJ ) then
!                                    SET FOR J = L1 / 2
!              ASSIGN 80 TO JGO
               JGO = 80
               go to 60
            end if
!           END OF COMPUTING WR-S AND WI-S

  130       L1 = L4
            JDIF = JDIF / 4
  140    continue
!        END OF L LOOP
      end if

      return
 
      END SUBROUTINE DFFT

!=====================================================================
end module DFFT_M

! $Log$
! Revision 2.3  2002/10/01 22:03:55  pwagner
! Fixed RCS Ident Block
!
! Revision 2.2  2002/10/01 20:06:00  bwknosp
! Added Id and RCS Info
!
! Revision 2.1  2002/09/06 22:31:51  vsnyder
! Initial commit
!
