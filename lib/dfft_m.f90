! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module DFFT_M

  implicit NONE
  private
  public :: DFFT, DRFT1, DTCST, FFT
  public :: InitSineTable, RFT1, TCST

  interface FFT
    module procedure DFFT
  end interface FFT

  interface InitSineTable
    module procedure DInitSineTable
  end interface InitSineTable

  interface RFT1
    module procedure DRFT1
  end interface RFT1

  interface TCST
    module procedure DTCST
  end interface TCST

  !---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: IdParm = &
    & "$Id$"
  character(len=len(idparm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
    & "$RCSfile$"
  private :: not_used_here
  !---------------------------------------------------------------------------

  integer, parameter :: RK = kind(0.0d0) ! Real Kind type parameter

  real(rk), parameter :: HALF = 0.5_rk
  ! PI4 = Pi/4
  real(rk), parameter :: PI4 = 0.7853981633974483096156608458198757_rk
  ! SPI4 = Sin(Pi/4) = Sqrt(2)/2
  real(rk), parameter :: SPI4 = .7071067811865475244008443621048490_rk


  ! THE DIMENSION OF KE MUST BE AS LARGE AS THE MAXIMUM VALUE
  ! PERMITTED FOR MM.

  ! NEEDST= .TRUE. if the sine table must be computed.
  ! MT    = base 2 log(NT)
  ! NT    = number of entries in the sine table
  ! MM    = base 2 log(number of complex fourier coefficients)
  ! KS    = distance in memory between successive points.  The i-th
  !         coefficient, a(i) is given by AR((I+1)*KS)+AI((I+1)*KS)*
  !         sqrt(-1), i=0, 1, ..., (2**MM)-1.
  ! ILAST = KS * 2**MM
  ! KE(L) = KS * 2**(MM-L), L=1, 2, ..., MM
  integer, parameter :: KEDIM = 30
  logical, save :: NEEDST = .true.
  integer, save :: MT, N1, NT, MM, KS, ILAST, KE(KEDIM), KEE(KEDIM+1)
  equivalence (KE(1), KEE(2)), (ILAST,KEE(1),N1)

contains

!=====================================================================
!=================== Conversion of the f77 fft.f to f90 ==============
!  The JPL real(rk) FFT Package.
!---------------------------------------------------------------------

      subroutine DRFT1 ( A, MODE, M, MS, S )

!>> 1994-11-11 DRFT1  Krogh   Declared all vars.
!>> 1994-10-20 DRFT1 Krogh  Changes to use M77CON
!>> 1989-05-07 DRFT1 FTK & CLL
!>> 1989-04-21 FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.

!     This subroutine computes Fourier transforms of real data using
!     the Cooley-Tukey fast Fourier transform.

      use ERMSG_M, only: ERM1, ERMSG

!     Variables in the calling sequence have the following types
      real(rk), intent(inout) :: A(:), S(:)
      integer, intent(in)     :: M
      integer, intent(inout)  :: MS
      character, intent(in)   :: MODE

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
!  if MODE = 'A' or 'a', do Fourier analysis, which amounts to setting
!  C(k) = sum for j=0 to N-1 of X(j)*T(M,j,k), for k = 0, N/2
!  with  T(M,j,k) = (1/N) * W ** (-j*k).
!  if MODE = 'S' or 's', do Fourier synthesis, which amounts to setting
!  X(j) = sum for k=0 to N-1 of C(k)*T(M,j,k), for j = 0, N - 1
!  with  T(M,j,k) = W ** (j*k)
!  (Recall that C(N-k) is the conjugate of C(k).)

! M is used to indicate N = 2**M, the number of real points.  The
!  number of points must satisfy 1 <= N <= 2**21.
!  M = 0 gives an immediate return.

! MS gives the state of the sine table.  if MS > 0, there are NT =
!    2 ** (MS-2 ) good entries in the sine table.  On the initial call,
!    MS must be set to 0, and when the return is made, it will be set
!    to M, which is the value of MS required for computing a
!    transform of size N.  if MS = -1, the sine table will be computed
!    as for the case MS = 0, and then a return to the user will be made
!    with MS set as before, but no transform will be computed.  This
!    option is useful if the user would like access to the sine table
!    before computing the FFT.
!    On detected errors the error message subrs are called and
!    execution stops.  if the user overrides the stop to cause
!    continuation, then this subr will return with MS = -2.

! S() is a vector, S(j) = sin(pi*j/2*NT)), j = 1, 2, ..., NT-1, where
!  NT is defined in the description of MS above.  S is computed by the
!  subroutine if M > MS.  (if S is altered, set MS=0 so that S
!  is recomputed.)

!     ------------------------------------------------------------------
!                Notes on local variables

!     MMAX is the largest value allowed for M
!     The dimension of KE must be at least as large as MMAX-1.

!     ANAL = .TRUE. when doing Fourier analysis, and .false. otherwise.

!     N1 = 2 ** M
!     ------------------------------------------------------------------
!--D replaces "?": ?RFT1, ?FFT, C?FFTC
!     Both versions use ERMSG, ERM1
!     and need ERFIN, IERV1
!     ------------------------------------------------------------------
      integer, parameter :: MMAX = KEDIM + 1
      INTEGER II, II1, II2, IR, IR1, IR2
      INTEGER J, JDIF, JJ
      INTEGER K1, K1N, KN2
      INTEGER L
      INTEGER MA, MSI
      INTEGER N1P

      real(rk) FN
      real(rk) T, TI, TT, TTI, TWO, WI, WR
      CHARACTER(LEN=19) :: MSG1 = 'Bad MODE, MODE =  .'
      LOGICAL ANAL

      PARAMETER (TWO = 2.0_rk)
!     ------------------------------------------------------------------

      if ( MODE == 'A' .or. MODE == 'a' ) then
         ANAL = .true.
      else if ( MODE == 'S' .or. MODE == 's' ) then
         ANAL = .false.
      else
         msg1(18:18) = mode
         call ermsg ( 'DRFT1', 1, 2, msg1, '.' )
         MS = -2
         return
      end if
      MA = M
      MSI = MS
      NEEDST = MA > MSI
      if ( NEEDST ) then
         if ( MA > MMAX .or. MA < 0 ) then
!                               Fatal error, default is to stop in ERM1
            call erm1 ( 'DRFT1', 1, 2, 'Require (0 <= M <= 21)', 'M', m, '.' )
            MS = -2
            return
         end if
         MS = MA
         MT = MA - 2
         call fft (A, A, S)
         if ( MSI == -1 ) return
      end if
      if ( MA /= 0 ) then
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
            a(1:n1) = a(1:n1) * fn
         else
!                              Set flags for Fourier synthesis
            IR = 1
            II = 2
            go to 40
         end if

!                              Compute complex Fourier transform
   20    do L = 1, MM
            KEE(L+1) = KEE(L) / 2
         end do ! L

         call fft (A(IR:), A(II:), S)
!                              End of computing complex transform

         if ( .NOT. ANAL ) return

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
         if ( MM > 0 ) then
!                           Special kase -- K1 = N1 / 4
            A(KN2+1) = TWO * A(KN2+1)
            A(KN2+2) = -TWO * A(KN2+2)
            if ( MM > 1 ) then
               J = 0
               do K1 = 3, KN2, 2
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
               end do ! K1
            end if
         end if
         if ( .NOT. ANAL ) go to 20
      end if

      return

      end subroutine DRFT1

!==========================================================================

      subroutine DFFT(AR, AI, S)

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

      real(rk), intent(inout) :: AR(:), AI(:), S(:)
!     Minimum dimensions are AR(ILAST), AI(ILAST), S(NT-1), where ILAST
!     and NT are module variables.

! AR and AI give the arrays used to hold the real and imaginary parts of
!     the Fourier coefficients and data.
! S   contains the sine table required in the calculations.

!     -----------------------------------------------------------------
!--D replaces "?": ?FFT, C?FFTC
!     -----------------------------------------------------------------
      integer I, I1, I2, I3, IJ
      integer J, JDIF, JI, JI2, JR
      integer KSI
      integer L, L1, L4, LJ

      real(rk) T, T1, T1I, T2, T2I, T3, T3I
      real(rk) TI, TP, TP1, TP1I, TPI
      real(rk) WI1, WI2, WI3, WR1, WR2, WR3
      logical SPCASE

!     -----------------------------------------------------------------

      if (NEEDST) then
!                      Compute the sine table
         call initSineTable ( s, mt )
      else
!                      Compute the transform
!                           Scramble A
         IJ = 1
         JI = 1
         L1 = KS
         if (MM > 1) then
   30       IJ = IJ + L1
            do I = 1, MM
               JI = JI + KE(I)
               KE(I) = -KE(I)
               if (KE(I) < 0) then
                  if (IJ < JI) then
!                       Interchange the IJ-th coefficient with the JI-th
                     T = AR(IJ)
                     AR(IJ) = AR(JI)
                     AR(JI) = T
                     T = AI(IJ)
                     AI(IJ) = AI(JI)
                     AI(JI) = T
                  end if
                  go to 30
               end if
            end do
         end if
!                          End of scrambling A
         JDIF = NT
         if (MOD(MM, 2) /= 0) then
!                    Special case -- L = 1,  MM odd  (radix 2 algorithm)
            L1 = L1 + L1
            do I = 1, ILAST, L1
               KSI = I + KS
               T = AR(I)
               AR(I) = T + AR(KSI)
               AR(KSI) = T - AR(KSI)
               T = AI(I)
               AI(I) = T + AI(KSI)
               AI(KSI) = T - AI(KSI)
            end do
            JDIF = JDIF / 2
         end if
!
         do L = 2, MM, 2
            L4 = 4 * L1
            LJ = L1 / 2
            SPCASE = .TRUE.
            J = 0
            JI = 0
!
!           Start of I loop  (radix 4 algorithm)
ij_l:       do
              IJ = J + 1
              do I = IJ, ILAST, L4
                 I1 = I + L1
                 I2 = I1 + L1
                 I3 = I2 + L1
                 if (SPCASE) then
                    if (J == 0) then
!                                     Special case -- J = 0
                       T = AR(I) + AR(I1)
                       T1 = AR(I) - AR(I1)
                       TI = AI(I) + AI(I1)
                       T1I = AI(I) - AI(I1)
                       T2 = AR(I2) + AR(I3)
                       T3 = AR(I2) - AR(I3)
                       T2I = AI(I2) + AI(I3)
                       T3I = AI(I2) - AI(I3)
                       go to 110
                    end if
!                                     Special case -- J = L1 / 2
                    T2 = SPI4 * AR(I2)
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
                 else
!
!                Usual case -- J /= 0  and  J /= L1 / 2
!
!                WRK and WIK (K = 1, 2, 3) are the real and imaginary part
!                resp. of exp(sqrt(-1) * PI * K*(2**(-L-MOD(MM, 2)))*J/KS)
!                               = exp(sqrt(-1) * PI * K * (J / (4 * L1)))
!
                    T2 = WR2 * AR(I1) - WI2 * AI(I1)
                    T2I = WI2 * AR(I1) + WR2 * AI(I1)
                    T = AR(I) + T2
                    T1 = AR(I) - T2
                    TI = AI(I) + T2I
                    T1I = AI(I) - T2I
                    TP = WR1 * AR(I2) - WI1 * AI(I2)
                    TPI = WI1 * AR(I2) + WR1 * AI(I2)
                    TP1 = WR3 * AR(I3) - WI3 * AI(I3)
                    TP1I = WI3 * AR(I3) + WR3 * AI(I3)
                 end if
                 T2 = TP + TP1
                 T3 = TP - TP1
                 T2I = TPI + TP1I
                 T3I = TPI - TP1I
  110            AR(I) = T + T2
                 AI(I) = TI + T2I
                 AR(I1) = T1 - T3I
                 AI(I1) = T1I + T3
                 AR(I2) = T - T2
                 AI(I2) = TI - T2I
                 AR(I3) = T1 + T3I
                 AI(I3) = T1I - T3
              end do ! I
!             End of I loop
!
              if (J <= LJ) then
                 if (SPCASE) then
                    if (J /= 0) &
            exit ij_l
                    J = KS
                    SPCASE = .FALSE.
                 else
                    J = L1 - J
!                   Compute WR-S and WI-S for J replaced by L1 - J
                    T = WI1
                    WI1 = WR1
                    WR1 = T
                    WR2 = -WR2
                    T = -WI3
                    WI3 = -WR3
                    WR3 = T
            cycle ij_l
                 end if
              else
                 J = L1 - J + KS
              end if
              if (J < LJ) then
!                                Compute WR-S and WI-S
                 JI = JI + JDIF
                 JR = NT - JI
                 WR1 = S(JR)
                 WI1 = S(JI)
                 JI2 = JI + JI
                 WI2 = S(JI2)
                 JR = NT - JI2
                 WR2 = S(JR)
                 JI2 = JI + JI2
                 if (JI2 <= NT) then
                    JR = NT - JI2
                    WR3 = S(JR)
                    WI3 = S(JI2)
            cycle ij_l
                 end if
                 JR = JI2 - NT
                 JI2 = NT - JR
                 WI3 = S(JI2)
                 WR3 = -S(JR)
            cycle ij_l
              else if (J == LJ) then
                 SPCASE = .TRUE.
            cycle ij_l
              end if
            exit ij_l
            end do ij_l
!           End of computing WR-S and WI-S
!
            L1 = L4
            JDIF = JDIF / 4
         end do ! L
!        End of L loop
      end if

      end subroutine DFFT

!=====================================================================

      subroutine DTCST (A, TCS, MODE, M, ND, MS, S)
!>> 1997-03-31 DTCST Krogh  Increased KEDIM, more sine table checks.
!>> 1996-01-23 DTCST Krogh  Changes to simplify conversion to C.
!>> 1994-11-11 DTCST Krogh  Declared all vars.
!>> 1994-10-20 DTCST Krogh  Changes to use M77CON
!>> 1989-06-16 DTCST FTK Fix error message on MODE, and TCS.
!>> 1989-06-05 WVS Change length of MODE and TCS from (ND) to (*)
!>> 1989-05-08 FTK & CLL
!>> 1989-04-21 FTK & CLL
!     .  Copyright (C) 1989, California Institute of Technology.
!     .  U. S. Government sponsorship under
!     .  NASA contract NAS7-918 is acknowledged.

!     This subroutine computes trigonometirc (sine-cosine), sine, or
!     cosine transforms of real data in up to 6 dimensions using the
!     Cooley-Tukey fast Fourier transform.

      use ERMSG_M, only: ERM1, ERMSG

!     Variables in the calling sequence have the following types
      real(rk), intent(inout)      :: A(:), S(:)
      integer, intent(in)          :: ND, M(ND)
      integer, intent(inout)       :: MS
      character(len=*), intent(in) :: TCS, MODE

!     Programmed by Fred F. Krogh at the Jet Propulsion Laboratory,
!     Pasadena, Calif.   August 1, 1969.
!     Revised for portability by Krogh -- January 29, 1988

!     Values for A, TCS, MODE, M, ND, and MS must be specified before
!     calling the subroutine.

!     In describing the usage the following notation is used
!     N(K) = 2 ** M(K)
!     MA = M(1) + M(2) + ... + M(ND)
!     NA = 2 ** MA

!     MTCS(K) = M(K)     TCS(K:K) = 'T'
!             = M(K)+1   otherwise

!     MX = MAX(MTCS(1), MTCS(2), ..., MTCS(ND))
!     NX = 2 ** MX

!     T(L,j,k) is defined differently for different values of TCS(L)

!       if TCS(L:L) = 'T' and MODE(L:L) = 'S', T(L,j,k)
!         =1/2                     if k = 0
!         =(1/2)*(-1)**j           if k = 1
!         =COS(j*k*PI/N(L))        if k IS EVEN  (k = 2, 4, ..., N(L)-2)
!         =SIN(j*(k-1)*PI/N(L))    if k IS ODD   (k = 3, 5, ..., N(L)-1)
!                    and if MODE(L:L) = 'A', T(L,j,k)
!         = (4/N) * (value of T(L,k,j) defined above)   if j<2
!         = (2/N) * (value of T(L,k,j) defined above)   Otherwise

!       if TCS(L:L) = 'C' and MODE(L:L) = 'S', T(L,j,k)
!         =1/2                     if k = 0
!         =COS(j*k*PI/N(L))           k = 1, 2, ..., N(L)-1
!         =(1/2)*(-1)**j           if k = N(L)
!                    and if MODE(L:L) = 'A', T(L,j,k)
!         = (2/N) * (value of T(L,j,k) defined above)

!       if TCS(L:L) = 'S' and MODE(L:L) = 'S', T(L,j,k)
!         =SIN(j*k*PI/N(L))           k = 0, 1, ..., N(L)-1
!                    and if MODE(L:L) = 'A', T(L,j,k)
!         = (2/N) * (value of T(L,j,k) defined above)

!     D(L) = N(L)     if TCS(L:L) /= 'C'
!          = N(L)+1   if TCS(L:L) = 'C'

!     The usage is as follows

! A() on input is an array of function values if one is doing Fourier
!   analysis, and is an array of Fourier coefficients if one is doing
!   Fourier synthesis.  On output, these are reversed.  In either case
!   A() is a real array with dimension A(D(1), D(2), ..., D(ND)).

! TCS  is a character variable of length ND.  The k-th character must be
!   'T' or 't' to select the general Trigonometric transform, or
!   'C' or 'c' to select the Cosine transform, or
!   'S' or 's' to select the Sine transform.
!     See the description of T(L,j,k) and M above.

! MODE  A character variable of length ND.  The k-th character must be
!   'A' or 'a' to select Analysis in the k-th dimension, or
!   'S' or 's' to select Synthesis in the k-th dimension.
!   One may be doing analysis, MODE(k:k) = 'A', with respect to one
!   dimension and synthesis, MODE(k:k) = 'S', with respect to
!   another.  A(j1+1, j2+1, ..., jND+1) is replaced by the sum over
!   0 <= k1 <= D(1)-1, 0 <= k2 <= D(2)-1, ..., 0 <= kND <=
!   D(ND)-1, of A(k1+1, k2+1, ..., kND+1) * T(1, k1, j1) * T(2, k2, j2)
!   ... * T(ND, kND, jND), 0 <= j1 <= D(1)-1, ..., 0 <= jND <=
!   D(ND)-1.

! M() is a vector used to indicate N(k) = 2**M(k)).  The number of
!   points in the k-th variable is then given by D(k) (see above).  M
!   must be specified in such a way that 0 < M(k) < 22
!   for k = 1, 2, ..., ND.

! ND is the dimension of (i.e. the number of subscripts in) the
!    array A.  ND must satisfy 1 <= ND <= 6.

! MS gives the state of the sine table.  if MS > 0, there are NT =
!    2 ** (MS-2 ) good entries in the sine table.  On the initial call,
!    MS must be set to 0, and when the return is made, it will be set
!    to MX, which is the value of MS required for computing a
!    transform of size N.  if MS = -1, the sine table will be computed
!    as for the case MS = 0, and then a return to the user will be made
!    with MS set as before, but no transform will be computed.  This
!    option is useful if the user would like access to the sine table
!    before computing the FFT.
!    On detected errors the error message subrs are called and
!    execution stops.  if the user overrides the stop to cause
!    continuation, then this subr will return with MS = -2.

! S() is a vector, S(j) = sin(pi*j/2*NT)), j = 1, 2, ..., NT-1, where
!  NT is defined in the description of MS above.  S is computed by the
!  subroutine if MX > MS.  (if S is altered, set MS=0 so that S
!  is recomputed.)
!     -----------------------------------------------------------------
!                Notes on local variables

!     NDMAX = the maximum value for ND, and MAXMX = the maximum
!     permitted for MTCS(1), ..., MTCS(ND)

!     NF(1) = 1, NF(K+1) = NF(K) * D(K) K = 1, 2, ..., ND

!     MU is used in the process of eliminating transforms with respect
!     to the first subscript of transforms with TCS(:) = 'S'.
!     (This is only necessary if ND > 1.)

!     The input character variable TCS is mapped to the internal
!     integer array ITCS() by mapping 'T' to 1, 'C' to 2, 'S' to 3.

!     -----------------------------------------------------------------
!--D replaces "?": ?TCST, ?FFT, C?FFTC
!     Both versions use ERM1
!     and need ERFIN, IERV1
!     -----------------------------------------------------------------
      integer, parameter :: NDMAX = 6
      integer, parameter :: MAXMX = KEDIM+1

      integer I, I1, II, IR, ITCS(NDMAX), ITCSK
      integer J, JDIF, JJ, JK
      integer K, KDR
      integer KII, KIN, KK, KKI, KKL, KKN
      integer L
      integer MA, MI, MMAX, MSI, MU(NDMAX)
      integer N, NDD, NDIV, NF(NDMAX+1)
      integer NI, NI1, NI2, NI2I, NTOT2

      character :: MSG1*13 = 'MODE(K:K) = ?', MSG2*12 = 'TCS(K:K) = ?'

      real(rk) FN, SUM, T, T1, TP, TPI, TS, TS1, WI, WR

      real(rk), parameter :: ZERO = 0.0_rk, ONE = 1.0_rk, TWO = 2.0_rk
      real(rk), parameter :: FOUR = 4.0_rk


!     -----------------------------------------------------------------

      NDD = ND
      if ( (NDD <= 0) .OR. (NDD > NDMAX) ) then
!                               FATAL ERROR, DEFAULT IS TO STOP IN ERM1
         call ERM1 ('DTCST', 1, 2, 'BAD ND', 'ND', ND, '.')
         MS = -2
         return
      end if
      MA = 0
      MMAX = 0
      NDIV = 1
! Every element in the array A is divided by NDIV before computing
! the transform.  The value computed for NDIV depends on whether
! one is doing analysis or synthesis and on the type of
! transform being computed.
      do K = 1, NDD
         MM = M(K)
         if ( (MM < 0) .OR. (MM > MAXMX) ) then
           !                         Fatal error, default is to stop in ERM1
           call ERM1 ('DTCST', 4, 2, 'Require 0 <= max(M(K)) <= 31', 'M', MM, '.')
           MS = -2
           return
         end if
         MA = MA + MM
         N = 2 ** MM
         if ( MODE(K:K) == 'A' .or. MODE(K:K) == 'a' ) then
            NDIV = NDIV * N
         else if ( MODE(K:K) == 'S' .or. MODE(K:K) == 's' ) then
            NDIV = 2 * NDIV
         else
            MSG1(13:13) = MODE(K:K)
            call ERM1 ('DTCST',2,2, MSG1, 'for K =',  K, '.')
            MS = -2
            return
         end if
         ITCSK = (index('TtCcSs', TCS(K:K)) + 1) / 2
         ITCS(K) = ITCSK
         if ( ITCSK == 0 ) then
            MSG2(12:12) = TCS(K:K)
            call ERM1 ('DTCST', 3, 2, MSG2, 'for K =', K, '.')
            return
         end if
         NF(1) = 1
         if ( ITCSK >= 2 ) then
            if ( ITCSK == 2) N = N + 1
            NDIV = 2 * NDIV
            MM = MM + 1
         end if
         NF(K+1) = NF(K) * N
         MMAX = max(MMAX, MM)
      end do ! K

      MSI = MS
      NEEDST = MMAX > MSI

      if ( .NOT. NEEDST ) then
!  Check internal parameters to catch certain user errors.
         if ( MT < KEDIM ) then
            if ( MMAX <= MT + 2 ) then
!              Skip sine table computation if all appears O.K.
               if ( MT <= 0 ) go to 15
               if ( abs(S(NT/2) - SPI4) <= 1.D-7 ) go to 15
            end if
         end if
         NEEDST = .true.
         call ERMSG('DTCST', 3, 1, 'Invalid sine table (re)computed', '.')
      end if
      MS = MMAX
      MT = MMAX - 2
      call fft (A, A, S)
      if ( MSI == -1 ) return
!                   All setup for computation now
   15 NTOT2 = NF(NDD+1)

      FN = ONE / DBLE(NDIV)
!     Divide every element of A by NDIV
      a(1:ntot2) = a(1:ntot2) * fn


!     Beginning of loop for computing multiple sum
k_l:  do K = 1, NDD
         ITCSK = ITCS(K)
         MI = M(K)
         MM = MI - 1
         if ( MODE(K:K) == 'A' .or. MODE(K:K) == 'a') MI = -MI
         KDR = NF(K)
         KS = KDR + KDR
         ILAST = NF(K+1)
         if ( ITCSK == 2) ILAST = ILAST - KDR
         do L = 1, MM
            KEE(L+1) = KEE(L) / 2
         end do ! L

         I = 1
         J = NDD
j_l:     do
           do L = 1, J
              MU(L) = 0
              if ( (L /= K) .AND. (ITCS(L) > 2) ) then
!             Skip the part of the array left empty by the sine transform
                 MU(L) = NF(L)
                 I = I + NF(L)
              end if
           end do ! L

ii_l:      do
           ! Compute indices associated with the current value of I (and K)
             I1 = I + KDR
             NI1 = I + NF(K+1)
             if ( ITCSK == 2) NI1 = NI1 - KDR
             NI = NI1 - KDR
             NI2 = (NI1 + I) / 2
             NI2I = NI2 + KDR
             if ( ITCSK /= 1 ) then
               !   Doing a cosine or a sine transform -- set MI = 0 and do
               !   calculations not required for sine-cosine transforms
               MI = 0
               J = NI
               SUM = A(I1)
               T = A(J)
               do
                  JK = J - KS
                  if ( JK < I1 ) exit
                  SUM = SUM + A(J)
                  A(J) = A(JK) - A(J)
                  J = JK
               end do
               if ( ITCSK /= 2 ) then
                 !                 Calculations for the sine transform
                 A(I) = TWO * A(I1)
                 A(I1) = -TWO * T
                 if ( MM == 0 ) go to 90
                 T = TWO * A(NI2)
                 A(NI2) = -TWO * A(NI2I)
                 A(NI2I) = T
                 go to 90
               end if
               !                  Set for cosine transform
               A(I1) = A(NI1)
             end if
             if ( MM == 0 ) go to 90
             if ( MI < 0 ) then
             ! Set for Fourier analysis
               II = I
               IR = II + KDR
               go to 110
             end if
             ! Begin calculations for the sine-cosine transform
   80        A(NI2) = TWO * A(NI2)
             A(NI2I) = TWO * A(NI2I)
   90        T = A(I) + A(I1)
             A(I1) = A(I) - A(I1)
             if ( MI < 0 ) then
                if ( ITCSK == 1 ) then
                   A(I1) = TWO * A(I1)
                   T = TWO * T
                end if
             end if
             A(I) = T
             J = 0
             JDIF = 2 ** (MT - MM + 1)
             KKL = KE(1) - KDR
             if ( MM > 1 ) then
                do KK = KS, KKL, KS
                   KKI = I + KK
                   KII = KKI + KDR
                   KKN = NI1 - KK
                   KIN = KKN + KDR
                   J = J + JDIF
                   WI = S(J)
                   JJ = NT - J
                   WR = S(JJ)
                   T = A(KKI) + A(KKN)
                   TS = A(KKN) - A(KKI)
                   T1 = A(KIN) - A(KII)
                   TS1 = A(KIN) + A(KII)
                   if ( ITCSK > 2 ) then
!                             The sine-cosine transform must be computed
!                             differently in the case of the sine transform
!                             because the input data is stored differently.
                      TP = WR * T - WI * T1
                      TPI = WI * T + WR * T1
                      A(KKI) = TP - TS1
                      A(KKN) = -TP - TS1
                      A(KII) = TPI + TS
                      A(KIN) = TPI - TS
                   else
                      TP = WR * TS1 + WI * TS
                      TPI = WI * TS1 - WR * TS
                      A(KKI) = T + TP
                      A(KKN) = T - TP
                      A(KII) = T1 + TPI
                      A(KIN) = TPI - T1
                   end if
                end do ! KK
             else if ( MM == 0 ) then
                go to 120
             end if
!            end of computing sine-cosine transform

             if ( MI < 0 ) go to 140
             IR = I
             II = IR + KDR

!            Compute a one-dimensional complex Fourier transform
  110        call fft (A(IR:), A(II:), S)
             if ( MI < 0 ) go to 80
             if ( MI /= 0 ) go to 140
  120        if ( ITCSK == 1 ) go to 140
             if ( ITCSK == 3 ) then
                A(I) = ZERO
             else
!                    Compute first and last elements of the cosine transform
                SUM = FOUR * SUM
                T = TWO * A(I)
                A(NI1) = T - SUM
                A(I) = T + SUM
                if ( MM >= 0) A(NI2) = TWO * A(NI2)
             end if
             if ( MM > 0 ) then
!                   Extra calculations required by sine and cosine transform
                J = 0
                JDIF = JDIF / 2
                do KK = KDR, KKL, KDR
                   KKI = I + KK
                   KKN = NI1 - KK
                   J = J + JDIF
                   WI = TWO * S(J)
                   T = A(KKI) + A(KKN)
                   TS = A(KKI) - A(KKN)
                   if ( ITCSK /= 2 ) then
                      T = T / WI
                   else
                      TS = TS / WI
                   end if
                   A(KKI) = T + TS
                   A(KKN) = T - TS
                end do ! KK
             end if

!            Logic for deciding which one-dimensional transform to do next
  140        J = 0
             do
                J = J + 1
                if ( J == K ) then
                   J = J + 1
                   I = I + NF(J) - NF(J-1)
                end if
                if ( J > NDD ) cycle k_l
                MU(J) = MU(J) + NF(J)
                if ( MU(J) < NF(J+1) ) exit
             end do
             I = I + NF(1)
             J = J - 1
             if ( j /= 0 ) exit ii_l
           end do ii_l
         end do j_l

      end do K_L

      return


      end subroutine DTCST

!=====================================================================

  subroutine DInitSineTable ( S, LogSize )
  ! Initialize a sine table S of size 2**LogSize - 1
    real(rk), intent(out) :: S(:)
    integer, intent(in) :: LogSize

    integer :: I, I1, I2, J, JJ, K, L, L1, LL
    real(rk) :: THETA
    mt = logSize
    needst = .true.
    nt = 2**mt
    if (mt > 0) then
!                       Set for L = 1
      j = nt
      jj = j / 2
      s(jj) = spi4
      if (mt > 1) then
        theta = pi4
        do l = 2, mt
          theta = half * theta
          k = j
          j = jj
          jj = jj / 2
!               At this point J = 2**(MT-L+1) and JJ = 2**(MT-L)
          s(jj) = sin(theta)
          l1 = nt - jj
          s(l1) = cos(theta)
          ll = nt - k
          if (ll >= j) then
            do i = j, ll, j
              i1 = nt - i
              i2 = i + jj
              s(i2) = s(i) * s(l1) + s(jj) * s(i1)
            end do ! I
          end if
        end do ! L
      end if
    end if

  end subroutine DInitSineTable 
   
!=====================================================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module DFFT_M

! $Log$
! Revision 2.9  2004/02/12 02:19:12  vsnyder
! Move DInitSineTable out of DFFT into a separate procedure, more f90-isms
!
! Revision 2.8  2004/01/29 19:54:35  vsnyder
! Remove unused variable
!
! Revision 2.7  2004/01/29 01:57:10  vsnyder
! Replace DFFT subroutine with version from Math77.6.1
!
! Revision 2.6  2004/01/16 01:50:01  vsnyder
! Remove superfluous declaration for I from DRFT1
!
! Revision 2.5  2003/07/15 22:41:19  vsnyder
! F90-ify stuff, add DTCST
!
! Revision 2.4  2002/10/08 00:09:08  pwagner
! Added idents to survive zealous Lahey optimizer

! Revision 2.3  2002/10/01 22:03:55  pwagner
! Fixed RCS Ident Block

! Revision 2.2  2002/10/01 20:06:00  bwknosp
! Added Id and RCS Info

! Revision 2.1  2002/09/06 22:31:51  vsnyder
! Initial commit

