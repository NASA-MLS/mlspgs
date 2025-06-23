! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module MLSNumbers              ! Some number theoretic datatypes, procedures
!=============================================================================

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use Dump_0, only: Dump
  use HighOutput, only: OutputNamedValue, StyledOutput
  use MLSFinds, only: FindFirst, FindLast
  use MLSSets, only: FindIntersection, Intersection, Union
  use MLSStringLists, only: ReadIntsFromList
  use MLSStrings, only: ReadNumsFromChars, Reverse_Trim, WriteIntsToChars
  use MLSStrings_0, only: StringContainsForbiddens
  use Optional_m, only: Default
  use Output_m, only: Output

  implicit none

  private
!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!         data types
! CompositeNum_T           A (composite) number expressed as prime factors
!                            and their powers; must be > 0
! RationalNum_T            The ratio between two composite numbers
! VeryBigInt_T             A signed integer stored as the string of its digits
! VeryBigReal_T            A signed real stored as the string of its digits and
!                            a separate signed exp
!
!         Functions, operations, routines
!     (In the following, Num_T, or A and B, can be any of the above data types)
! Add                      A + B
!                                                            n
! BinomialCoef             Compute the binomial coefficient ( )
!                                                            m
! BigPrime                 nth prime (for n up to 10^8)
! Copy                     Copy one Num_T into a second
! Create                   Represent an integer as a Num_T or Convert between
!                            different Num_Ts
! Destroy                  Deallocate the arrays in a Num_T 
! Divide                   Divide two Num_T; c/d = a/b
!                            (removing common factors so c, d have none)
! Dump                     Dump the components in a Num_T
! Estimate                 Express a Num_T as its log (base 10)
! Evaluate                 Express a Num_T as an integer
! Factorial                Compute n!, returned as a composite
! Fibonacci                Compute F[n], the nth number in the Fibonacci sequence
!                            1, 1, 2, 3, 5, .., (F[n-1]+F[n-2])
! FindAllPrimes            Compute all the prime numbers the supplied array
!                            can hold
! GreatestCommonDivisor    Find the gcd of two composite nums
! inRange                  Is the int represented by the composite num in range?
!                           (processor-dependent)
! IntFactorial             Compute n!, returned as an integer
! Invert                   Interchange a rational's numerator and denominator
! LeastCommonMultiple      Find the lcm of two composite nums
! LogBinomialCoef          Compute the logarithm of binomial coefficient
! LogFactorial             Compute logarithm of n!, returned as a d.p. float
! Multiply                 A * B
! Power                    A ^ n
! Reduce                   Reduce numerator and denomincator of a rational 
!                            number by dividing out their gcd
! Subtract                 A - B
! AreRelativelyPrime       gcd(A,B) = 1 ?
! IsGreaterThan            A > B ?
! IsEqual                  A = B?
! IsLessThan               A < B ?
! IsOne                    A = 1 ?
! IsZero                   A = 0 ?
! IsPrime                  Is arg prime?
! NextPrime                Next prime > arg
! Prime                    nth prime (if n < 1230)
! PrimeFactors             break arg into its prime factors
! Trim                     trim off non-essential 0s representing VeryBigInt
! === (end of toc) ===

! === (start of api) ===
! Add ( Num_T A, Num_T B, Num_T C )
! log AreRelativelyPrime( compNum A, compNum B )
! int BigPrime( int n, [int primes(:)] )
! int BinomialCoef ( int n, int m )
! Copy ( Num_T A, Num_T C )
! Create ( int n, Num_T C )
! Destroy ( Num_T C )
! Divide ( compNum A, compNum B, compNum C, compNum D )
! Dump ( Num_T C, [char* options] )
! Estimate ( real logn, Num_T C )
! Evaluate ( int n, Num_T C )
! type(compNum) factorial ( int n, [int k] )
! int Fibonacci( int n, int f1, int f2 )
! FindAllPrimes( int primes(:) )
! GreatestCommonDivisor ( compNum A, compNum B, compNum C )
! int IntFactorial ( int n, [int k] )
! int Fibonacci ( int n, [int F1], [int F2] )
! log InRange( compNum C )
! Invert ( rationalNum R )
! log IsEqual( Num_T A, Num_T B )
! log IsEqual( int n, Num_T C )
! log IsGreaterThan( Num_T A, Num_T B )
! log IsLessThan( Num_T A, Num_T B )
! log IsPrime( compNum C )
! LeastCommonMultiple ( compNum A, compNum B, compNum C )
! dble LogBinomialCoef ( int n, int m )
! dble LogFactorial ( int n, [int k] )
! Multiply ( Num_T A, Num_T B, Num_T C )
! Power ( n, Num_T C )
! Reduce ( compNum A, compNum B, compNum C )
! Reduce ( rationalNum R )
! log IsPrime( int n, [int primes(:)] )
! int NextPrime( int n, [int primes(:)] )
! int Prime( int n )
! int PrimeFactors( int n, int[:] factors, [int[:] powers] )
! Subtract ( Num_T A, Num_T B, Num_T C )
!
! Note: Initially, this module uses only prime numbers below 10000 efficiently
! and correctly
! They account for the first 1229 such primes.
! For larger prime numbers, the function prime poops out, returning -999.
! Instead use BigPrime. It may be slow but it can handle primes below 10^8.
! For still larger primes, or to precompute enough primes to efficiently
! use them later, BigPrime, NextPrime, and isPrime can accept
! an optional array primes. You can use FindAllPrimes to precompute them.
! === (end of api) ===

  public :: BigPrime, FindAllPrimes, IsPrime, NextPrime, Prime, PrimeFactors
  public :: Fibonacci
  public :: CompositeNum_T, RationalNum_T, VeryBigInt_T, VeryBigReal_T
  public :: Copy, Create, Destroy, Dump, Estimate, Evaluate, InRange
  public :: AreRelativelyPrime, IsEqual, IsGreaterThan, IsLessThan, &
    & IsOne, IsZero
  public :: EstimateCompositeNum, EvaluateCompositeNum
  public :: GreatestCommonDivisor, LeastCommonMultiple
  public :: Add, Subtract, Multiply, Divide, Power, Reduce, Invert
  public :: BinomialCoef, CompBinomialCoef, LogBinomialCoef
  public :: Factorial, IntFactorial, LogFactorial
  public :: IntTimesVeryBigInt, Multiply_very

  interface Add
    module procedure Add_composite, Add_rational, Add_very
  end interface

  interface Copy
    module procedure Copy_composite, Copy_rational
  end interface

  interface Create
    module procedure Create_composite, Create_rational, Create_very, &
      & CreateVeryFromComposite
  end interface

  interface Destroy
    module procedure Destroy_composite, Destroy_rational
  end interface

  ! -- We must first decide:
  ! -- how to handle integer division:
  !       (1) Produce a Rational? A truncated integer?
  !       (2) Integer and a remainder?
  ! -- Dividing Rationals seems straightforward
  ! interface Divide
  !  module procedure Divide_composite, Divide_rational, Divide_very
  ! end interface

  interface Dump
    module procedure Dump_composite, Dump_rational, Dump_Very
  end interface

  interface Estimate
    module procedure EstimateCompositeNum_single, EstimateCompositeNum_double
    module procedure EstimateRationalNum_single, EstimateRationalNum_double
    module procedure EstimateVeryBigInt
  end interface

  interface EstimateCompositeNum
    module procedure EstimateCompositeNum_single, EstimateCompositeNum_double
  end interface

  interface EstimateRationalNum
    module procedure EstimateRationalNum_single, EstimateRationalNum_double
  end interface

  interface Evaluate
    module procedure EvaluateCompositeNum, EvaluateVeryBigInt
  end interface

  interface inRange
    module procedure inRange_comp, inRange_very
  end interface

  interface isEqual
    module procedure isEqual_int, isEqual_composite, isEqual_rational, &
      & isEqual_very
  end interface

  interface isGreaterThan
    module procedure isGreaterThan_composite, isGreaterThan_rational, &
      & isGreaterThan_very
  end interface

  interface isLessThan
    module procedure isLessThan_composite, isLessThan_rational, isLessThan_very
  end interface

  interface isOne
    module procedure IsOne_composite, isOne_very
  end interface

  interface isPrime
    module procedure isPrime_int, isPrime_composite
  end interface

  interface isZero
    module procedure isZero_very
  end interface

  interface Multiply
    module procedure IntTimesVeryBigInt, &
      & Multiply_composite, Multiply_rational, Multiply_very
  end interface

  interface Power
    module procedure Power_composite, Power_rational, Power_Very
  end interface

  interface Reduce
    module procedure Reduce_composite, Reduce_rational
  end interface

  interface Subtract
    module procedure Subtract_composite, Subtract_rational, Subtract_Very
  end interface

  ! Should we move these datatype and their procedures 
  ! to a separate module? Say MLSExtendedNumbers?
  
  ! ------------------------------------------------------------
  ! This datatype represents a positive integer as a composite 
  ! of its prime factors
  ! I.e., n = Product ( factors[k] ^ powers[k] )
  ! where each factors[k] is a prime number, and powers[k] > 0
  ! Why might this be useful?
  ! (1) Represents even very large non-prime integers without loss of precision
  ! (2) Some operations become very easy, e.g. lcm, gcd, isPrime, power

  ! Beware of calling EvaluateCompositeNum without first checking the size of its
  ! log by a call to inRange; if it is outside the range we
  ! can represent, EvaluateCompositeNum returns "-1"

  ! The Estimate procedures provide an alternate way to get at the value,
  ! returning its logarithm as either a single- or double-precision float

  ! "1" is represnted by c%factors remaining unassociated
  type CompositeNum_T
    integer, dimension(:), allocatable  :: factors
    integer, dimension(:), allocatable  :: powers
  end type
  ! ------------------------------------------------------------

  ! ------------------------------------------------------------
  ! This datatype represents a positive non-integer as a ratio of two integers
  type RationalNum_T
    type(CompositeNum_T)                :: numerator
    type(CompositeNum_T)                :: denominator
  end type
  ! ------------------------------------------------------------

  ! ------------------------------------------------------------
  ! Let's represent very big signed integers as the sequence of characters
  ! they would occupy if we printed each Value via 
  !    print *, Reverse_trim(Value)
  ! Note that to avoid having to do many Reverse_Trims
  ! we store the digits in reverse order.
  ! So the integer 72310 is stored as '01327'
  ! Beware! 
  !   '0' passed to IsZero returns false--should be ' '
  !   '000000' > '1' (because it's a longer string)
  !   '10' /= 1
  ! To prevent these bugs and gotchas, apply the module's own Trim function
  ! How many characters can VeryBigInt_T hold?
  integer, parameter :: MaxNumVeryCharS = 128

  type VeryBigInt_T
    character(len=MaxNumVeryCharS)      :: cval     = ' '
    logical                             :: negative = .false.
  end type
  ! Here is where a parameterized type might have been handy
  ! Note also how extravagant this type is with memory--a packed
  ! version might easily hold the range 0-99 in each position
  ! instead of just 0-9. That's twice as efficient.
  ! ------------------------------------------------------------

  ! ------------------------------------------------------------
  ! These are basically 
  !      (+ or -) verybigint * 10^exponent
  ! E.g., '4315' and exponent=5 represents
  !     (.5134 * 10^5) or 51340.0
  type VeryBigReal_T
    character(len=MaxNumVeryCharS)      :: cval     = ' '
    logical                             :: negative = .false.
    integer                             :: exponent = 0 ! Can be signed
  end type
  ! Here is where a type extension of VeryBigInt_T might have been handy
  ! ------------------------------------------------------------

  ! How many unique factors can a CompositeNum_T hold?
  integer, parameter :: MAXNUMFACTORS = 128

  ! ------------------------------------------------------------
  ! How many Fibonacci numbers do we track?
  integer, parameter :: MAXNUMFIBS = 12
  integer, dimension(MAXNUMFIBS), save :: fibnumbers = (/ &
    & 1 ,           1 ,           2 ,           3 ,           5 , &
    & 8 ,          13 ,          21 ,          34 ,          55 , &
    & 89 ,         144 /)

  ! How many primes do we track?
  integer, parameter :: MAXNUMPRIMES = 1229
  integer, dimension(MAXNUMPRIMES), save :: primenumbers = -999

  ! How many factorials do we track? How large can we calculate
  integer, parameter :: MAXNUMFACS = 12
  integer, parameter :: MAXNINFACTORIAL = 12
  integer, dimension(MAXNUMFACS), save :: Factorials = (/ &
    & 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600  &
    & /)

contains

  ! --------------- Add_composite --------------------
  ! Add two composite nums; c = a + b
  subroutine Add_composite( A, B, C, sbtrct )
    ! Args
    type(CompositeNum_T), intent(in)  :: A
    type(CompositeNum_T), intent(in)  :: B
    type(CompositeNum_T), intent(out) :: C
    logical, optional, intent(in)     :: sbtrct ! For internal use: A - B?
    ! Local variables
    type(CompositeNum_T)              :: APB
    type(CompositeNum_T)              :: AOG ! A/gcd
    type(CompositeNum_T)              :: BOG ! B/gcd
    type(CompositeNum_T)              :: gcd ! B/gcd
    integer                           :: na
    integer                           :: nb
    integer                           :: napb
    integer                           :: signB
    ! Executable
    signB = 1
    if ( present(sbtrct) ) signB = merge(-1, 1, sbtrct)
    ! 1st
    ! Find the gcd
    call GreatestCommonDivisor( A, B, gcd, AOG, BOG )
    ! Form the sum AOG+BOG
    call EvaluateCompositeNum( nA, AOG )
    call EvaluateCompositeNum( nB, BOG )
    napb = abs(nA + signB*nB)
    ! A + B = (nA + nB)*gcd
    call Create ( napb, APB )
    call Multiply ( APB, gcd, C )
  end subroutine Add_composite

  ! --------------- Add_rational --------------------
  ! Add two rational nums; cr = ar * br; Reduce result
  ! Plan: strive to express the ratios with a common denominator
  !   a       c       ad + bc      a LCM/b + c LCM/d
  !   -   +   -    =     -    =            -
  !   b       d         bd                LCM
  subroutine Add_rational( AR, BR, CR, sbtrct )
    ! Args
    type(RationalNum_T), intent(in)  :: AR
    type(RationalNum_T), intent(in)  :: BR
    type(RationalNum_T), intent(out) :: CR
    logical, optional, intent(in)    :: sbtrct
    ! Local variables
    type(CompositeNum_T)             :: A
    type(CompositeNum_T)             :: B
    type(CompositeNum_T)             :: C
    type(CompositeNum_T)             :: D
    type(CompositeNum_T)             :: LCMOB  ! LCM/B
    type(CompositeNum_T)             :: LCMOD  ! LCM/D
    type(CompositeNum_T)             :: LCM
    type(CompositeNum_T)             :: T1
    type(CompositeNum_T)             :: T2
    type(CompositeNum_T)             :: T3
    ! integer                          :: signB
    ! Executable
    ! signB = 1
    ! if ( present(sbtrct) ) signB = merge(-1, 1, sbtrct)
    call Copy ( AR%numerator  , A ) ! AR = a/b
    call Copy ( AR%denominator, B )
    call Copy ( BR%numerator  , C ) ! BR = c/d
    call Copy ( BR%denominator, D )
    ! 1st: find lcm of the two denominators
    call LeastCommonMultiple ( B, D, LCM, LCMOB, LCMOD )
    ! Form a LCM/b and c LCM/d
    call Multiply ( A, LCMOB, T1 )
    call Multiply ( C, LCMOD, T2 )
    call Add( T1, T2, T3, sbtrct )
    ! Express T3/LCM as A/B
    call Destroy ( A )
    call Destroy ( B )
    call Divide( T3, LCM, A, B )
    ! Now A/B is the rational CR
    call Copy ( A, CR%numerator )
    call Copy ( B, CR%denominator )
  end subroutine Add_rational

  ! --------------- Add_Very --------------------
  ! Add two Very nums; c = a + b
  subroutine Add_Very( A, B, C, sbtrct )
    ! Args
    type(VeryBigInt_T), intent(in), target  :: A
    type(VeryBigInt_T), intent(in), target  :: B
    type(VeryBigInt_T), intent(out) :: C
    logical, optional, intent(in)   :: sbtrct ! For internal use: A - B?
    ! Local variables
    type(VeryBigInt_T), pointer     :: alpha
    type(VeryBigInt_T), pointer     :: beta
    type(VeryBigInt_T)              :: gamma
    character(len=3)                :: row
    character(len=1)                :: na     
    character(len=1)                :: nb     
    logical                         :: negate
    character(len=1)                :: op
    logical                         :: subtractB
    logical, parameter              :: verbose = .false.
    ! Executable
    subtractB  = .false.
    if ( present(sbtrct) ) subtractB = sbtrct
    ! 1st
    ! Consider the table where we show what operations will be needed
    !             sbtract?  F             T          row
    ! A>0.B>0           |A| + |B|     |A| - |B|      +.+
    ! A>0.B<0           |A| - |B|     |A| + |B|      +.-
    ! A<0.B>0           |B| - |A|    -(|A| + |B|)    -.+
    ! A<0.B<0          -(|A| + |B|)   |B| - |A|      -.-
    ! 
    ! So we must choose
    ! (a) alpha, beta: which comes 1st, A or B?
    ! (b) op: which operation do we perform, addition or subtraction?
    ! (c) negate: do we negate the whole result?
    
    ! Some useful defaults
    alpha => A
    beta => B
    negate = .false.
    ! Format the row as 'n.m' where n is 1 if pos, 0 if neg
    na = merge ('+', '-', .not. A%negative )
    nb = merge ('+', '-', .not. B%negative )
    row = na // '.' // nb
    select case(row)
    case ('+.+')
      if ( .not. subtractB ) then
        op = '+'
      else
                                       op = '-'
      endif
    case ('+.-')
      if ( .not. subtractB ) then
        op = '-'
      else
                                       op = '+'
      endif
    case ('-.+')
      if ( .not. subtractB ) then
        op = '-'
        alpha => B
        beta => A
      else
                                       op = '+'
                                       negate=.true.
      endif
    case ('-.-')
      if ( .not. subtractB ) then
        op = '+'
        negate = .true.
      else
                                       op = '-'
                                       alpha => B
                                       beta => A
      endif
    case default
      ! How did we end up here?
      print *, 'Oops--unlikely value for row in Add_very: ', row
      stop
    end select
    if ( verbose ) then
      print *, 'op: ', op
      print *, 'negate: ', negate
      print *, 'row: ', row
    endif
    select case (op)
    case ('+')
      call Add_us
    case ('-')
      call Subtract_us
    end select
    C = Trim_very(gamma) ! In case any trailing 0's were created during op
    if ( negate ) C%negative = .not. gamma%negative
  contains
    subroutine Add_us
      ! Do the tedious, time-consuming task of adding 2 multi-char strings
      ! Not optimized for speed--an obvious improvement would deal with
      ! say blocks of 4 characters at a time
      type(VeryBigInt_T)              :: aa
      type(VeryBigInt_T)              :: bb
      ! type(VeryBigInt_T)              :: cc
      integer                         :: carry ! 1 if w carry into next position
      integer                         :: ia, ib, ic ! int values at position
      integer                         :: k     ! position in string
      !
      if ( isZero(alpha) ) then
        gamma = beta
        return
      elseif ( isZero(beta) ) then
        gamma = alpha
        return
      endif
      ! Swap aa and bb if |alpha| < |beta| so aa is longest string
      if ( isLessThan(alpha, beta) ) then
        bb%cval =  alpha%cval 
        aa%cval =  beta%cval 
      else
        aa%cval =  alpha%cval 
        bb%cval =  beta%cval 
      endif
      carry = 0
      do k=1, len_trim(aa%cval) ! The longer string
        ! If blank, ReadNumsFromChars will return '-999'
        call ReadNumsFromChars ( aa%cval(k:k), ia )
        call ReadNumsFromChars ( bb%cval(k:k), ib ) 
        ib = max(ib, 0) ! Turns '-999' -> '0'
        ic = carry + ia + ib
            if ( verbose ) print *, 'ia, ib, ic: ', ia, ib, ic
        if ( ic > 9 ) then
          ! We must carry 10 into the next position
          ic = ic - 10
          carry = 1 ! For next position
        else
          carry = 0 ! For next position
        endif
        call WriteIntsToChars ( ic, gamma%cval(k:k) )
      enddo
      ! Do we have a carry?
      if ( carry > 0 ) call WriteIntsToChars ( 1, gamma%cval(k:k) )
    end subroutine Add_us
    subroutine Subtract_us
      ! Do the tedious, time-consuming task of subtracting 2 multi-char strings
      ! Not optimized for speed--an obvious improvement would deal with
      ! say blocks of 4 characters at a time
      type(VeryBigInt_T)              :: aa
      type(VeryBigInt_T)              :: bb
      ! type(VeryBigInt_T)              :: cc
      integer                         :: borrow ! 1 if we borrow from next position
      integer                         :: ia, ib, ic ! int values at position
      integer                         :: k     ! position in string
      logical                         :: swapped ! If Alpha < Beta
      !
      if ( isEqual(alpha, beta) ) then
        return ! assume gamma enters world initilized as ' '
      elseif ( isZero(alpha) ) then
        gamma%cval = beta%cval
        gamma%negative = .not. beta%negative
        return
      elseif ( isZero(beta) ) then
        gamma = alpha
        return
      endif
      ! Swap aa and bb if |alpha| < |beta| so aa is longest string
      if ( isLessThan(alpha, beta, absval=.true.) ) then
        bb%cval =  alpha%cval 
        aa%cval =  beta%cval 
        swapped = .true.
      else
        aa%cval =  alpha%cval 
        bb%cval =  beta%cval 
        swapped = .false.
      endif
      borrow = 0
      do k=1, len_trim(aa%cval) ! The longer string
        ! If blank, ReadNumsFromChars will return '-999'
        call ReadNumsFromChars ( aa%cval(k:k), ia )
        call ReadNumsFromChars ( bb%cval(k:k), ib ) 
        ib = max(ib, 0) ! Turns '-999' -> '0'
        ic = ia - ib - borrow
        if ( ic < 0 ) then
          ! We must borrow 10 from the next position
          ic = ic + 10
          borrow = 1 ! For next position
        else
          borrow = 0 ! For next position
        endif
        call WriteIntsToChars ( ic, gamma%cval(k:k) )
      enddo
      gamma%negative = swapped
    end subroutine Subtract_us
  end subroutine Add_Very
  
  function IsGreaterThan_Very ( A, B, absval ) result ( so )
    ! args
    type(VeryBigInt_T), intent(in)          :: A
    type(VeryBigInt_T), intent(in)          :: B
    logical, optional, intent(in)           :: absval ! |A| > |B|?
    logical                                 :: so
    ! Executable
    so = isLessThan( A, B, grtr=.true., absval=absval )
  end function IsGreaterThan_Very

  function IsEqual_Very ( A, B ) result ( so )
    ! args
    type(VeryBigInt_T), intent(in)          :: A
    type(VeryBigInt_T), intent(in)          :: B
    logical                                 :: so
    ! Executable
    so = (trim(A%cval) == trim(B%cval)) .and. (A%negative .eqv. B%negative)
    ! If A and B are both zero, we don't care if either claims to be negative
    if ( isZero(A) .and. isZero(B) ) so = .true.
  end function IsEqual_Very

  function IsLessThan_Very ( A, B, grtr, absval ) result ( so )
  ! A < B ?

  !         *   *   *
  ! A trick! If grtr is present, and true, we will
  ! return whether -A < -B
  ! (which is the same as A > B)
  ! If absval is present and true, we will compare absolute values only
    ! args
    type(VeryBigInt_T), intent(in), target  :: A
    type(VeryBigInt_T), intent(in), target  :: B
    logical, optional, intent(in)           :: grtr ! A > B?
    logical, optional, intent(in)           :: absval ! |A| < |B|?
    logical                                 :: so
    ! Internal variables
    type(VeryBigInt_T), pointer             :: alpha
    type(VeryBigInt_T), pointer             :: beta
    type(VeryBigInt_T), pointer             :: gamma
    type(VeryBigInt_T), target              :: AA
    type(VeryBigInt_T), target              :: BB
    integer                                 :: ia, ib, k
    logical                                 :: absolute
    logical                                 :: reversed
    ! Executable
    so = .false.
    if ( IsEqual(A,B) ) return
    reversed = .false.
    if ( present(grtr) ) reversed = grtr
    absolute = .false.
    if ( present(absval) ) absolute = absval
    AA = A
    BB = B
    ! Here's how we use absolute values
    if ( absolute ) then
      AA%negative = .false.
      BB%negative = .false.
    endif
    if ( .not. reversed ) then
      alpha => AA
      beta => BB
    else
      alpha => BB
      beta => AA
    endif
    if ( alpha%negative .neqv. beta%negative ) then
      ! One of alpha, beta is < 0; so is true if it's alpha
      so = alpha%negative
      return
    elseif ( alpha%negative .and. beta%negative ) then
      ! Must swap alpha and beta because we'll compare absolute values below
      gamma => alpha
      alpha => beta
      beta => gamma
    endif
    ! Two more special cases
    if ( IsZero(alpha) ) then
      so = .true.
      return
    elseif ( IsZero(beta) ) then
      so = .false.
      return
    endif
    ! Before we begin the painful character-by-character comparison, 
    ! check on lengths: the longer the greater
    if ( len_trim(alpha%cval) /= len_trim(beta%cval) ) then
      so = (len_trim(alpha%cval) < len_trim(beta%cval) )
      return
    endif
    ! Now for the painful business
    ! We'll start with the most significant digit, which
    ! is the final because that's how we decided to implement this type
    do k = len_trim(alpha%cval), 1, -1
      call ReadNumsFromChars ( alpha%cval(k:k), ia )
      call ReadNumsFromChars ( beta%cval(k:k), ib )
      if ( ia < ib ) then
        so = .true.
        return
      elseif ( ib < ia ) then
        so = .false.
        return
      endif
    enddo
    ! If we reached this point, the two must be equal--but w ruled that out first
    call output( '***How are A and B equal? ' // &
      & 'We checked already in IsLessThan_Very', advance='yes' )
    call output( 'A: ' // trim(A%cval), advance='yes' )
    call output( 'B: ' // trim(B%cval), advance='yes' )
    stop
  end function IsLessThan_Very

  ! --------------------- IsZero ---------------
  function IsZero_Very ( arg ) result ( so )
    ! args
    ! Assumes you did any necessary Trim_very before calling us
    type(VeryBigInt_T), intent(in)          :: arg
    logical                                 :: so
    ! Executable
    so = (len_trim(arg%cval) < 1)
  end function IsZero_Very

  ! --------------- AreRelativelyPrime --------------------
  ! Are the two numbers relatively prime?
  ! (1) gcd(A,B) == 1?
  logical function AreRelativelyPrime( A, B )
    ! Args
    type(CompositeNum_T), intent(in)       :: A, B
    ! Local variables
    integer, dimension(:), allocatable     :: C1, C2
    integer                                :: n
    ! Executable
    AreRelativelyPrime = .true.
    if ( isOne ( a ) .or. &
      & isOne ( b ) ) return ! gcd(any,"1") == "1"
    call allocate_test( C1, size(A%factors), 'C1', &
      & ModuleName // 'AreRelativelyPrime C1' )
    call allocate_test( C2, size(A%factors), 'C2', &
      & ModuleName // 'AreRelativelyPrime C2' )
    call FindIntersection ( A%factors, B%factors, C1, C2, n )
    AreRelativelyPrime = ( n == 0 )
  end function AreRelativelyPrime

  ! --------------- BigPrime --------------------
  function BigPrime( n, primes ) result( next )
    ! Returns the next prime number greater than the arg n
    ! Method:
    ! Check if n < MAXNUMPRIMES
    ! If n < MAXNUMPRIMES use prime(n)
    ! otherwise, loop over nextprime however many times you need
    ! (This may be agonizingly slow so please be patient!)
    ! Args
    integer :: next
    integer, intent(in) :: n ! E.g., if n=100 returns 101 which is next prime > 100
    integer, optional, dimension(:), intent(in) :: primes ! instead of primenumbrs
    integer :: i
    integer :: nprimes
    ! Executable
    if ( present(primes) ) then
      nprimes = size(primes)
      if ( n <= nprimes ) then
        next = primes(n)
        return
      end if
      next = primes(nprimes)
      do i=nprimes+1, n
        next = nextprime( next+2, primes )
      enddo
      return
    endif
    if ( n <= MAXNUMPRIMES ) then
      next = prime(n)
      return
    end if
    next = prime(MAXNUMPRIMES)
    do i=MAXNUMPRIMES+1, n
      next = nextprime(next+2)
    enddo
  end function BigPrime

  ! --------------- BinomialCoef --------------------
  !                               n
  ! Form the integer bin. coef   ( )
  !                               m
  ! We must calculate n! / ( m! (n-m)! )
  !
  ! Warning--if the result is too large to fit within an integer, we may get -1
  ! See also CompBinomialCoef
  function BinomialCoef( n, m ) result ( c )
    ! Args
    integer, intent(in)               :: n
    integer, intent(in)               :: m
    integer                           :: c
    ! Local variables
    type(CompositeNum_T)              :: A
    ! Executable
    c = 1
    if ( n < 2 .or. m < 1 .or. (n-m) < 1 ) return ! m = 0 and n = m both => 1
    A = CompBinomialCoef( n, m )
    call EvaluateCompositeNum( c, A )
    call Destroy( A )
  end function BinomialCoef

  ! --------------- CompBinomialCoef --------------------
  !                               n
  ! Form the Composite result of ( )
  !                               m
  ! We must calculate n! / ( m! (n-m)! )
  ! See also BinomialCoef
  function CompBinomialCoef( n, m ) result ( C )
    ! Args
    integer, intent(in)               :: n
    integer, intent(in)               :: m
    type(CompositeNum_T)              :: C
    ! Local variables
    type(CompositeNum_T)              :: A
    type(CompositeNum_T)              :: B
    type(CompositeNum_T)              :: D
    ! Executable
    call Create( 1, C ) ! Default value is 1, not 0
    if ( n < 2 .or. m < 1 .or. (n-m) < 1 ) return ! m = 0 and n = m both => 1
    call Destroy_composite( C )
    ! 1st--multiply m! * (n-m)!
    A = factorial (m)
    C = factorial (n-m)
    call Multiply ( A, C, B )
    call Destroy_composite( A )
    call Destroy_composite( C )
    ! 2nd--divide n! by this product
    A = factorial (n)
    call Divide ( A, B, C, D )
    call Destroy_composite( A )
    call Destroy_composite( B )
    call Destroy_composite( D ) ! But we are guarenteed that D = 1
  end function CompBinomialCoef

  ! --------------- Copy_composite --------------------
  ! Copy one Composite "a" into a second "c"
  subroutine Copy_composite( A, C )
    ! Args
    type(CompositeNum_T), intent(in)     :: A
    type(CompositeNum_T), intent(out)    :: C
    ! Local variables
    integer                              :: nFactors
    ! Executable
    if ( isOne ( a ) ) return ! If "1" return "1"
    nFactors = size ( a%factors )
    if ( nFactors < 1 ) return
    call allocate_test( C%factors, nFactors, 'factors', &
      & ModuleName // 'Copy_composite' )
    call allocate_test( C%powers, nFactors, 'powers', &
      & ModuleName // 'Copy_composite' )
    C%factors = a%factors
    C%powers  = a%powers
  end subroutine Copy_composite

  ! --------------- Copy_rational --------------------
  ! Copy one rational "a" into a second "c"
  subroutine Copy_rational( A, C )
    ! Args
    type(RationalNum_T), intent(in)     :: A
    type(RationalNum_T), intent(out)    :: C
    ! Local variables
    type(CompositeNum_T)                :: numerator
    type(CompositeNum_T)                :: denominator
    ! Executable

    call Copy( A%numerator, numerator )
    call Copy( A%denominator, denominator )
    C = RationalNum_T( numerator, denominator )
  end subroutine Copy_rational

  ! --------------- Create_composite --------------------
  ! Form a Composite from an input integer
  subroutine Create_composite( n, C )
    ! Args
    integer, intent(in)                  :: n
    type(CompositeNum_T), intent(out)    :: C
    ! Local variables
    integer, dimension(MAXNUMFACTORS) :: factors
    integer, dimension(MAXNUMFACTORS) :: powers
    integer                           :: nFactors
    ! Executable
    if ( n < 2 ) return
    nFactors = primeFactors ( n, factors, powers )
    call allocate_test( C%factors, nFactors, 'factors', &
      & ModuleName // 'Create_composite' )
    call allocate_test( C%powers, nFactors, 'powers', &
      & ModuleName // 'Create_composite' )
    C%factors = factors(1:nFactors)
    C%powers  = powers (1:nFactors)
  end subroutine Create_composite

  ! --------------- Create_rational --------------------
  ! Form a rational from two input integers
  subroutine Create_rational( numerator, denominator, R )
    ! Args
    integer, intent(in)                  :: numerator
    integer, intent(in)                  :: denominator
    type(RationalNum_T), intent(out)     :: R
    ! Local variables
    type(CompositeNum_T)                 :: A
    type(CompositeNum_T)                 :: B
    ! Executable
    call Create( numerator, A )
    call Create( denominator, B )
    R = RationalNum_T( A, B )
    call Reduce ( R )
  end subroutine Create_rational

  ! --------------- Create_very --------------------
  ! Form a VeryBigInt from an input integer
  ! See also other Creates
  subroutine Create_very( n, C )
    ! Args
    integer, intent(in)                  :: n
    type(VeryBigInt_T), intent(out)      :: C
    ! Local variables
    ! Executable
    if ( n < 0 ) then
      call WriteIntsToChars ( abs(n), C%cval )
      C%negative  = .true.
    elseif ( n > 0 ) then
      call WriteIntsToChars ( n, C%cval )
      C%negative  = .false.
    else
      C%cval = ' '
      C%negative  = .false.
      return
    endif
    C%cval = reverse_trim(C%cval)
  end subroutine Create_very

  ! --------------- CreateVeryFromComposite --------------------
  ! Create a VeryBigInt From a Composite as an output integer

  ! See also EstimateCompositeNum
  subroutine CreateVeryFromComposite( C, V )
    ! Args
    type(CompositeNum_T), intent(in)  :: C
    type(VeryBigInt_T), intent(out)   :: V
    ! Local variables
    integer                           :: i
    integer                           :: nFactors
    type(VeryBigInt_T)                :: F, T ! temps
    ! Executable
    call Create(1, V ) ! C defaults to 1, not 0
    if ( isOne ( c ) ) return
    nFactors = size ( C%factors )
    if ( nFactors < 1 ) return
    do i=1, nFactors
      call Create ( C%factors(i), F ) ! This is the ith factor
      call Power( C%powers(i), F ) ! Now raised to the nth power
      T = V
      call Multiply( F, T, V ) ! the ith contribution to C
    end do
  end subroutine CreateVeryFromComposite


  ! --------------- Destroy_composite --------------------
  ! Destroy a Composite, deallocating its components
  ! In effect reducing it to "1"
  !
  ! 2 common uses of ths procedure:
  ! (1) Reclaiming the memory from a composite type before reallocating it
  ! (2) Setting it to "1"
  !
  ! The use of allocatables relieves you of having to Destroy
  ! temporary data types
  subroutine Destroy_composite( C )
    ! Args
    type(CompositeNum_T), intent(inout)    :: C
    ! Executable
    call deallocate_test( C%factors, 'factors', &
      & ModuleName // 'Create_composite' )
    call deallocate_test( C%powers, 'powers', &
      & ModuleName // 'Create_composite' )
  end subroutine Destroy_composite

  ! --------------- Destroy_rational --------------------
  ! Destroy a rational, deallocating its components
  ! In effect reucing it to "1"
  !
  ! 2 common uses of ths procedure:
  ! (1) Reclaiming the memory from a rational type before reallocating it
  ! (2) Setting it to "1"
  !
  ! The use of allocatables relieves you of having to Destroy
  ! temporary data types
  subroutine Destroy_rational( R )
    ! Args
    type(RationalNum_T), intent(inout)    :: R
    ! Executable
    call Destroy( R%numerator )
    call Destroy( R%denominator )
  end subroutine Destroy_rational

  ! --------------- Divide --------------------
  ! Divide two composite nums; c/d = a/b
  ! where c and d have no common factors
  !
  ! It does not create a rational_t, however:
  ! to do that you must 
  !          call Create ( int numerator, int denominator, rationalNum R )
  subroutine Divide( A, B, C, D )
    ! Args
    type(CompositeNum_T), intent(in)  :: A
    type(CompositeNum_T), intent(in)  :: B
    type(CompositeNum_T), intent(out) :: C
    type(CompositeNum_T), intent(out) :: D
    ! Local variables
    type(CompositeNum_T)              :: G ! GCM of A, B
    ! Executable
    if ( isOne ( A ) ) then
      call Copy_composite ( b, d ) ! a is 1, so c is 1
    else if ( isOne ( B ) ) then
      call Copy_composite ( a, c ) ! b is 1 so d is 1
    else
      call GreatestCommonDivisor( A, B, G )
      call Reduce( A, G, C )
      call Reduce( B, G, D )
    end if
  end subroutine Divide

  ! --------------- Dump_composite --------------------
  ! Dump a Composite, dumping its components
  subroutine Dump_composite( C, options )
  ! Optionally, its estimate, too
    ! Args
    type(CompositeNum_T), intent(in)       :: C
    character(len=*), optional, intent(in) :: options
    double precision                       :: logN
    ! Executable
    call StyledOutput( 'Dump of composite number', &
      & Default( options, '--Headline') )
    if ( isOne ( c ) ) then
      call output( '1 (its factors are unassociated)', advance='yes' )
    else
      call Dump( C%factors, 'prime factors' )
      call Dump( C%powers, 'their powers' )
    end if
    if ( .not. present(options) ) return
    if ( index(options, 'v') > 0 ) then
      call EstimateCompositeNum ( logN, C )
      call outputNamedValue( 'log(value)', logN )
    endif
  end subroutine Dump_composite

  ! --------------- Dump_rational --------------------
  ! Dump a Rational, dumping its components
  ! Optionally, its estimate, too
  subroutine Dump_rational( R, options )
    ! Args
    type(RationalNum_T), intent(in)        :: R
    character(len=*), optional, intent(in) :: options
    double precision                       :: logN
    ! Executable
    call StyledOutput( 'Dump of rational number; first the numerator', &
      & Default( options, '--Headline') )
    call Dump( R%numerator )
    call StyledOutput( 'next the denominator', Default( options, '--Headline') )
    call Dump( R%denominator )
    if ( .not. present(options) ) return
    if ( index(options, 'v') > 0 ) then
      call EstimateRationalNum ( logN, R )
      call outputNamedValue( 'log(value)', logN )
    endif
  end subroutine Dump_rational

  subroutine Dump_very( v, options )
    ! Args
    type(VeryBigInt_T), intent(in)         :: v
    character(len=*), optional, intent(in) :: options
    ! Executable
    call StyledOutput( 'Dump of very big int', &
      & Default( options, '--Headline') )
    call output ( reverse_trim(v%cval), advance='yes' )
    call output ( 'negative? ', advance='no' )
    call output ( v%negative, advance='yes' )
  end subroutine Dump_very

  ! --------------- EstimateCompositeNum --------------------
  ! Return a single or -double precision estimate
  ! of a composite number as the value of its logarithm
  !
  ! Why as its log? Is that choice so intuitively obvious?
  subroutine EstimateCompositeNum_single( logN, C )
    ! Args
    real, intent(out)                 :: logN
    type(CompositeNum_T), intent(in)  :: C
    ! Local variables
    integer                           :: i
    integer                           :: nFactors
    ! Executable
    logN = 0 ! Default value is 0
    if ( isOne ( c ) ) return
    nFactors = size ( C%factors )
    if ( nFactors < 1 ) return
    ! 1st--Must know if n would be outside range
    logN = 0.
    do i=1, nFactors
      logN = logN + c%powers(i)*log10(1.*C%factors(i))
    end do
  end subroutine EstimateCompositeNum_single

  subroutine EstimateCompositeNum_double( logN, C )
    ! Args
    double precision, intent(out)     :: logN
    type(CompositeNum_T), intent(in)  :: C
    ! Local variables
    integer                           :: i
    integer                           :: nFactors
    ! Executable
    logN = 0 ! Default value is 0
    if ( isOne ( c ) ) return
    nFactors = size ( C%factors )
    if ( nFactors < 1 ) return
    ! 1st--Must know if n would be outside range
    logN = 0.d0
    do i=1, nFactors
      logN = logN + c%powers(i)*log10(1.d0*C%factors(i))
    end do
  end subroutine EstimateCompositeNum_double

  ! --------------- EstimateRationalNum --------------------
  ! Return a single or -double precision estimate
  ! of a rational number in the value of its logarithm
  ! Note that the log(x/y) = log(x) - log(y)
  subroutine EstimateRationalNum_single( logN, R )
    ! Args
    real, intent(out)                 :: logN
    type(RationalNum_T), intent(in)   :: R
    ! Local variables
    real                              :: logNumerator
    real                              :: logDenominator
    ! Executable
    call EstimateCompositeNum ( logNumerator,   R%numerator )
    call EstimateCompositeNum ( logDenominator, R%denominator )
    logN = logNumerator - logDenominator
  end subroutine EstimateRationalNum_single

  subroutine EstimateRationalNum_double( logN, R )
    ! Args
    double precision, intent(out)     :: logN
    type(RationalNum_T), intent(in)   :: R
    ! Local variables
    double precision                  :: logNumerator
    double precision                  :: logDenominator
    ! Executable
    call EstimateCompositeNum ( logNumerator,   R%numerator )
    call EstimateCompositeNum ( logDenominator, R%denominator )
    logN = logNumerator - logDenominator
  end subroutine EstimateRationalNum_double

  subroutine EstimateVeryBigInt( logN, C )
    ! Return 0 if 0, 1, or < 0
    ! Args
    ! Method:
    ! Treat the cval as a sequence of digits after he decimal point
    ! instead of before
    ! Note that they're in reverse order (Who thought up that trick?)
    ! See also EvaluateVeryBigInt
    double precision, intent(out)     :: logN
    type(VeryBigInt_T), intent(in)    :: C
    ! Local variables
    double precision                  :: colMultiplier
    integer                           :: ic
    integer                           :: k
    integer                           :: nTerms
    ! Executable
    logN = 0 ! Default value is 0
    if ( isZero ( c ) .or. isOne ( c ) .or. C%negative ) return
    nTerms = len_trim ( C%cval )
    if ( nTerms < 1 ) return ! Should have caught this earlier, but still ..
    colMultiplier = 1.
    do k=nTerms, 1, -1
      colMultiplier = 0.1*colMultiplier
      call ReadNumsFromChars ( c%cval(k:k), ic )
      logN = logN + ic*colMultiplier
    end do
    logN = log10(logN)
    logN = nTerms + logN
  end subroutine EstimateVeryBigInt

  ! --------------- EvaluateCompositeNum --------------------
  ! Evaluate a Composite as an output integer
  ! Warning-- the result may be outiside the range of integers
  ! representable on the machine; if it is outside,
  ! we return 1.

  ! See also EstimateCompositeNum
  subroutine EvaluateCompositeNum( n, C )
    ! Args
    integer, intent(out)              :: n
    type(CompositeNum_T), intent(in)  :: C
    ! Local variables
    integer                           :: i
    integer                           :: nFactors
    ! Executable
    n = 1 ! Default value is 1, not 0
    if ( isOne ( c ) ) return
    nFactors = size ( C%factors )
    if ( nFactors < 1 ) return
    ! 1st--Must know if n would be outside range
    if ( .not. inRange( C ) ) then
      ! print *, 'logN > range(n) ', logN, range(n)
      n = -1 ! The value returned if n would be > Huge
      return
    end if
    do i=1, nFactors
      n = n * C%factors(i)**C%powers(i)
    end do
  end subroutine EvaluateCompositeNum

  ! --------------- EvaluateVeryBigInt --------------------
  ! Evaluate a verybigint as an output integer
  ! Warning-- the result may be outiside the range of integers
  ! representable on the machine; if it is outside,
  ! we return 0. (Or maybe just bomb? What if InRange is wrong?)

  ! See also EstimateVeryBigInt
  subroutine EvaluateVeryBigInt( n, C )
    ! Args
    integer, intent(out)              :: n
    type(VeryBigInt_T), intent(in)    :: C
    ! Executable
    n = 0 ! Default value is 0
    if ( isZero ( c ) .or. .not. InRange(C) ) return
    call ReadNumsFromChars ( reverse_trim(C%cval), n )
    if ( n == Stringcontainsforbiddens ) n = 0
    if ( C%negative ) n = -n
  end subroutine EvaluateVeryBigInt

  ! --------------- Factorial --------------------
  ! Form the Composite result of n!.
  ! Optionally, may find double-factorial, by setting k=2
  ! E.g., 6!! = 6*4*2. 7!! = 7*5*3

  ! See also IntFactorial, EstimateCompositeNum
  function Factorial( n, k ) result ( C )
    ! Args
    integer, intent(in)               :: n
    integer, optional, intent(in)     :: k ! In case we need double-factorial
    type(CompositeNum_T)              :: C
    ! Local variables
    type(CompositeNum_T)              :: A
    type(CompositeNum_T)              :: B
    integer                           :: i
    integer                           :: kstep
    ! Executable
    kstep = 1
    if ( present(k) ) kstep = max(k, 1)
    call Create( 1, C ) ! Default value is 1, not 0
    if ( n < 2 ) return ! 0! and 1! both = 1
    i = n
    do
      Call Copy( C, B )
      call Destroy_composite( C )
      call Create( i, A )
      call Multiply_composite( A, B, C )
      i = i - kstep
      ! Housekeeping
      call Destroy_composite( A )
      call Destroy_composite( B )
      if ( i < 2 ) exit ! No point in multiplying by 1
    end do
  end function Factorial

  ! --------------- IntFactorial --------------------
  ! Find the integer result of n!.
  ! If n is too large, returns -1
  ! Optionally, may find double-factorial, by setting k=2
  ! E.g., 6!! = 6*4*2. 7!! = 7*5*3

  ! See also Factorial
  recursive function IntFactorial( n, k ) result ( f )
    ! Args
    integer, intent(in)               :: n
    integer, optional, intent(in)     :: k ! In case we need double-factorial
    integer                           :: f
    ! Internal variables
    integer                           :: kstep
    ! Executable
    kstep = 1
    if ( present(k) ) kstep = max(k, 1)
    f = 1
    if ( n < 2 ) return
    if ( n <= MAXNUMFACS ) then
      if ( kstep == 1 ) then
        f = Factorials(n)
      else
        f = n * IntFactorial( n-kstep, k )
      end if
    else if ( n <= MAXNINFACTORIAL-kstep ) then
      f = n * IntFactorial( n-kstep, k )
    else
      f = -1
    end if

  end function IntFactorial

  ! --------------- Fibonacci --------------------
  ! Find the inth number in the Fibonacci sequence
  ! 1, 1, 2, 3, 5, ..
  ! using the fact that
  ! F[n] = F[n-1] + F[n-2]
  ! 
  ! Optionally, use the generalized Fibonacci sequence
  ! where you supply F[1] and F[2]
  recursive function Fibonacci( n, f1, f2 ) result ( f )
    ! Args
    integer, intent(in)               :: n
    integer, optional, intent(in)     :: f1
    integer, optional, intent(in)     :: f2
    integer                           :: f
    ! Internal variables
    integer                            :: fnm2 ! F[n-2]
    integer                            :: fnm1 ! F[n-1]
    integer                            :: i
    ! Executable
    if ( present(f1) .and. present(f2) ) then
      ! Generalized Fibonacci sequence
      if ( n < 2 ) then
        f = f1
      elseif ( n == 2 ) then
        f = f2
      else
        fnm2 = f1
        fnm1 = f2
        do i=3, n
          ! Fib's rule
          f = fnm1 + fnm2
          ! Get ready for the next iteration
          fnm2 = fnm1
          fnm1 = f
        enddo
      endif
      return
    endif
    f = 1
    if ( n < 3 ) return
    if ( n <= MAXNUMFIBS ) then
      f = fibnumbers(n)
     elseif( n > 46 ) then
       f = -1
     else
       ! This form is a million times slower (no exaggeration)
       ! f = Fibonacci(n-1) + Fibonacci(n-2)
       f = Fibonacci( n, 1, 1 )
     endif
  end function Fibonacci

  ! --------------- FindAllPrimes --------------------
  ! Precompute an array of the first n primes where
  ! n is the size of the array that is to hold them.
  subroutine FindAllPrimes( primes )
    ! Args
    integer, dimension(:), intent(out)        :: primes
    ! Local variables
    integer                                   :: i
    integer                                   :: n
    ! Executable
    n = size(primes)
    if ( n <= MAXNUMPRIMES ) then
      primes = primenumbers(1:n)
      return
    endif
   ! Now Fill primes array
   do i=1, MAXNUMPRIMES
     primes(i) = primenumbers( i )
   enddo
   do i=MAXNUMPRIMES+1, n
     primes(i) = Bigprime( i, primes(:i-1) )
   enddo
  end subroutine FindAllPrimes

  ! --------------- GreatestCommonDivisor --------------------
  ! Find the gcd of two composite nums
  ! Method:
  ! (1) Find the factors in the set
  !       {c%factors[i]} = {a%factors[i]} /\ {b%factors[i]}
  ! (2) For each c%factors(i), choose the smaller of the corresponding a or b
  !      power
  ! Note that if either A or B is "one", C will be "one", too
  ! by a side-effect of declaring C with intent(out)
  subroutine GreatestCommonDivisor( A, B, C, AOC, BOC )
    ! Args
    type(CompositeNum_T), intent(in)                 :: A
    type(CompositeNum_T), intent(in)                 :: B
    type(CompositeNum_T), intent(out)                :: C
    type(CompositeNum_T), optional, intent(out)      :: AOC ! A/C
    type(CompositeNum_T), optional, intent(out)      :: BOC ! B/C
    ! Local variables
    integer                                          :: i
    integer                                          :: j
    integer                                          :: na
    integer                                          :: nb
    ! Executable
    ! Check whether  or B is "one"
    if (  allocated(a%factors) .and. allocated(b%factors) ) then
      ! c%factors = Intersection( a%factors, b%factors )
      ! c%powers =  Intersection( a%powers, b%powers )
      c = CompositeNum_T( &
        & Intersection( a%factors, b%factors ), &
        & Intersection( a%factors, b%factors ) &
        & )
      do i=1, size(c%factors)
        na = 0
        nb = 0
        j = FindFirst( a%factors, c%factors(i) )
        na = a%powers(j)
        j = FindFirst( b%factors, c%factors(i) )
        nb = b%powers(j)
        c%powers(i) = min(na, nb)
      end do
    end if
    if ( present(AOC) ) call Reduce_Composite( A, C, AOC )
    if ( present(BOC) ) call Reduce_Composite( B, C, BOC )
  end subroutine GreatestCommonDivisor

  ! --------------- InRange --------------------
  ! Is the integer represented by the composite argument in range?
  ! i.e., not too big
  logical function InRange_comp( C )
    ! Args
    type(CompositeNum_T), intent(in) :: C
    ! Local variables
    real                             :: logN
    integer                          :: n
    ! Executable
    InRange_comp = .true.
    if ( isOne ( c ) ) return ! "1" is not too big
    call EstimateCompositeNum ( logN, C )
    InRange_comp = .not. ( logN > real(range(n)) )
  end function InRange_comp

  logical function InRange_very( C )
    ! Args
    type(VeryBigInt_T), intent(in) :: C
    ! Executable
    InRange_very = len_trim(C%cval) <= digits(1)
  end function InRange_very

  ! --------------- Invert --------------------
  ! Invert a rational number by interchanging numerator and denominator
  subroutine Invert( R )
    ! Args
    type(RationalNum_T), intent(inout) :: R
    ! Local variables
    type(CompositeNum_T)              :: C ! R%numerator
    type(CompositeNum_T)              :: D ! R%denominator
    ! Executable
    call Copy( R%numerator, C )
    call Copy( R%denominator, D )
    call Destroy( R )
    R = RationalNum_T( D, C )
  end subroutine Invert

  ! --------------- IsEqual_composite --------------------
  ! Are the two numbers equal?
  ! (1) size(c*factors) == 1
  ! (2) c%powers(1) == 1
  logical function IsEqual_composite( A, B )
    ! Args
    type(CompositeNum_T), intent(in) :: A, B
    isEqual_composite = .true.
    if ( isOne ( a ) .and. &
      & isOne ( b ) ) return ! "1" == "1"
    isEqual_composite = .false.
    if ( isOne ( a ) .or. &
      & isOne ( b ) ) return ! "1" is not Equal to "n"
    if ( size(a%factors) /= size(b%factors) ) return
    isEqual_composite = all(a%factors == b%factors) .and. &
      & all(a%powers == b%powers)
  end function IsEqual_composite

  ! --------------- IsEqual_int --------------------
  ! Are the two numbers equal?
  ! (1) size(c*factors) == 1
  ! (2) c%powers(1) == 1
  logical function isEqual_int( n, A )
    ! Args
    integer, intent(in)              :: n
    type(CompositeNum_T), intent(in) :: A
    ! Internal variables
    integer                          :: nA
    ! Executable
    if ( n < 1 ) then
      isEqual_int = .false.
      return
    else if ( .not. allocated(a%factors) ) then
      isEqual_int = ( n == 1 )
      return
    end if
    call EvaluateCompositeNum( nA, A )
    isEqual_int = ( n == nA )
  end function isEqual_int

  ! --------------- IsEqual_rational --------------------
  ! Are the two numbers equal?
  ! Having been reduced, are numerators and denominators equal?
  logical function IsEqual_rational( A, B )
    ! Args
    type(RationalNum_T), intent(in) :: A, B
    ! Internal variables
    type(RationalNum_T)             :: C, D
    ! Executable
    call Copy( A, C )
    call Copy( B, D )
    Call Reduce( C )
    Call Reduce( D )
    isEqual_rational = isEqual( C%numerator, D%numerator ) .and. &
      & isequal( C%denominator, D%denominator )
  end function IsEqual_rational

  ! --------------- IsGreaterThan --------------------
  logical function IsGreaterThan_composite( A, B )
    type(CompositeNum_T), intent(in) :: A, B
    IsGreaterThan_composite = IsLessThan_composite (A, B, grtr=.true. )
  end function IsGreaterThan_composite

  logical function IsGreaterThan_rational( A, B )
    type(rationalNum_T), intent(in) :: A, B
    IsGreaterThan_rational = IsLessThan_rational (A, B, grtr=.true. )
  end function IsGreaterThan_rational

  ! --------------- IsOne --------------------
  logical function IsOne_composite( A )
    type(CompositeNum_T), intent(in) :: A
    IsOne_composite = ( .not. allocated(A%factors) )
  end function IsOne_composite

  function IsOne_Very ( arg ) result ( so )
    ! args
    type(VeryBigInt_T), intent(in)          :: arg
    logical                                 :: so
    ! Executable
    so = (trim(arg%cval) == '1' .and. .not. arg%negative)
  end function IsOne_Very

  ! --------------- IsLessThan --------------------
  ! A < B ?

  !         *   *   *
  ! A trick! If grtr is present, and true, we will
  ! return whether -A < -B
  ! (which is the same as A > B)
  logical function IsLessThan_composite( A, B, grtr )
    ! Args
    type(CompositeNum_T), intent(in) :: A, B
    logical, optional, intent(in)    :: grtr ! A > B?
    ! Internal variables
    double precision                 :: D_A, D_B
    integer                          :: signB
    ! Executable
    signB = 1
    if ( present(grtr) ) signB = merge(-1, 1, grtr)
    IsLessThan_composite = .false.
    if ( isEqual ( a, b ) ) return
    call EstimateCompositeNum_double( D_A, A )
    call EstimateCompositeNum_double( D_B, B )
    IsLessThan_composite = ( signB*D_A < signB*D_B )
  end function IsLessThan_composite

  logical function IsLessThan_rational( A, B, grtr )
    ! Args
    type(RationalNum_T), intent(in)  :: A, B
    logical, optional, intent(in)    :: grtr ! A > B?
    ! Internal variables
    double precision                 :: D_A, D_B
    integer                          :: signB
    ! Executable
    signB = 1
    if ( present(grtr) ) signB = merge(-1, 1, grtr)
    IsLessThan_rational = .false.
    if ( isEqual ( a, b ) ) return
    call EstimateRationalNum_double( D_A, A )
    call EstimateRationalNum_double( D_B, B )
    IsLessThan_rational = ( signB*D_A < signB*D_B )
  end function IsLessThan_rational

  ! --------------- IsPrime_composite --------------------
  ! A version of the function for composite arguments
  ! Returns FALSE unless both
  ! (1) size(c*factors) == 1
  ! (2) c%powers(1) == 1
  ! See how easy that was?
  logical function IsPrime_composite( C )
    ! Args
    type(CompositeNum_T), intent(in) :: C
    isPrime_composite = .false.
    if ( isOne ( c ) ) return ! "1" is not prime
    if ( size(c%factors) /= 1 ) return
    if ( c%powers(1) == 1 ) isPrime_composite = .true.
  end function IsPrime_composite

  ! --------------- LeastCommonMultiple --------------------
  ! Find the lcm of two composite nums
  ! Method:
  ! (1) Find the factors in the set
  !       {c%factors[i]} = {a%factors[i]} U {b%factors[i]}
  ! (2) For each c%factors(i), choose the larger of the corresponding a or b
  !      power
  subroutine LeastCommonMultiple( A, B, C, COA, COB )
    ! Args
    type(CompositeNum_T), intent(in)                   :: A
    type(CompositeNum_T), intent(in)                   :: B
    type(CompositeNum_T), intent(out)                  :: C
    type(CompositeNum_T), optional, intent(out)        :: COA  ! C/A
    type(CompositeNum_T), optional, intent(out)        :: COB  ! C/A
    ! Local variables
    integer                           :: i
    integer                           :: j
    integer                           :: na
    integer                           :: nb
    ! Executable
    if ( isOne ( a ) ) then
      call Copy_composite ( b, c )
    else if ( isOne ( b ) ) then
      call Copy_composite ( a, c )
    else
      ! c%factors = Union( a%factors, b%factors )
      ! c%powers =  Union( a%powers, b%powers )
      c = CompositeNum_T( &
        & Union( a%factors, b%factors ), &
        & Union( a%factors, b%factors ) &
        & )
      do i=1, size(c%factors)
        na = 0
        nb = 0
        j = FindFirst( a%factors, c%factors(i) )
        if ( j > 0 ) na = a%powers(j)
        j = FindFirst( b%factors, c%factors(i) )
        if ( j > 0 ) nb = b%powers(j)
        c%powers(i) = max(na, nb)
      end do
    end if
    if ( present(COA) ) call Reduce ( C, A, COA )
    if ( present(COB) ) call Reduce ( C, B, COB )
  end subroutine LeastCommonMultiple

  ! --------------- LogBinomialCoef --------------------
  !                                   n
  ! Form the d.p. log of bin. coef   ( )
  !                                   m
  ! We must calculate log10( n! / ( m! (n-m)! ) )
  !
  ! See also CompBinomialCoef
  function LogBinomialCoef( n, m ) result ( c )
    ! Args
    integer, intent(in)               :: n
    integer, intent(in)               :: m
    double precision                  :: c
    ! Local variables
    type(CompositeNum_T)              :: A
    ! Executable
    c = 0.d0
    if ( n < 2 .or. m < 1 .or. (n-m) < 1 ) return ! m = 0 and n = m both => log10(1)
    A = CompBinomialCoef( n, m )
    call EstimateCompositeNum( c, A )
  end function LogBinomialCoef

  ! --------------- LogFactorial --------------------
  ! Form the d.p. log of n!
  !
  ! See also CompFactorial
  function LogFactorial( n, k ) result ( c )
    ! Args
    integer, intent(in)               :: n
    integer, optional, intent(in)     :: k
    double precision                  :: c
    ! Local variables
    type(CompositeNum_T)              :: A
    ! Executable
    c = 0.d0
    if ( n < 2  ) return ! 0! and 1! are both 1, so logs are 0
    A = Factorial( n, k )
    call EstimateCompositeNum( c, A )
  end function LogFactorial

  ! --------------- Multiply_composite --------------------
  ! Multiply two composite nums; c = a * b
  subroutine Multiply_composite( A, B, C )
    ! Args
    type(CompositeNum_T), intent(in)  :: A
    type(CompositeNum_T), intent(in)  :: B
    type(CompositeNum_T), intent(out) :: C
    ! Local variables
    integer                           :: i
    integer                           :: j
    integer                           :: na
    integer                           :: nb
    ! Executable
    if ( isOne ( a ) ) then
      call Copy_composite ( b, c )
    else if ( isOne ( b ) ) then
      call Copy_composite ( a, c )
    else
      c = CompositeNum_T( &
        & Union( a%factors, b%factors ), &
        & Union( a%factors, b%factors ) &
        & )
      do i=1, size(c%factors)
        na = 0
        nb = 0
        j = FindFirst( a%factors, c%factors(i) )
        if ( j > 0 ) na = a%powers(j)
        j = FindFirst( b%factors, c%factors(i) )
        if ( j > 0 ) nb = b%powers(j)
        c%powers(i) = na + nb
      end do
    end if
  end subroutine Multiply_composite

  ! --------------- Multiply_rational --------------------
  ! Multiply two rational nums; c = a * b; Reduce result
  subroutine Multiply_rational( A, B, C )
    ! Args
    type(RationalNum_T), intent(in)  :: A
    type(RationalNum_T), intent(in)  :: B
    type(RationalNum_T), intent(out) :: C
    ! Local variables
    type(CompositeNum_T)             :: N
    type(CompositeNum_T)             :: D
    ! Executable
    call Multiply( A%numerator, B%numerator, N )
    call Multiply( A%denominator, B%denominator, D )
    C = RationalNum_T ( N, D )
    call Reduce ( C )
  end subroutine Multiply_rational

  subroutine intTimesVeryBigInt ( int, very, product )
    ! Compute the product of an int Times a veryBigInt
    ! Args
    integer, intent(in)             :: int
    type(VeryBigInt_T), intent(in)  :: very
    type(VeryBigInt_T), intent(out) :: product
    integer                         :: carry ! 1 if w carry into next position
    integer                         :: ia, ib, ic ! int values at position
    integer                         :: k     ! position in string
    !
    call Create (0, product)
    if ( isZero(very) .or. int==0 ) then
      return
    elseif ( int == 1 ) then
      product = very
      return
    elseif ( isOne(very) ) then
      call Create (int, product)
      return
    endif
    carry = 0
    do k=1, len_trim(very%cval)
      call ReadNumsFromChars ( very%cval(k:k), ia )
      ic = carry + int*ia
      if ( ic > 9 ) then
        ! We must carry ic/10 into the next position
        carry = ic/10
        ic = ic - 10*carry
      else
        carry = 0 ! For next position
      endif
      call WriteIntsToChars ( ic, product%cval(k:k) )
    enddo
    ! Do we have a carry?
    if ( carry > 0 ) call WriteIntsToChars ( carry, product%cval(k:k) )
    product%negative = very%negative
  end subroutine intTimesVeryBigInt

  ! --------------- Multiply_very --------------------
  ! Multiply two very big ints; c = a * b
  subroutine Multiply_very( A, B, C )
    ! Args
    type(VeryBigInt_T), intent(in)  :: A
    type(VeryBigInt_T), intent(in)  :: B
    type(VeryBigInt_T), intent(out) :: C
    ! Local variables
    integer                           :: k
    integer                           :: ib
    type(VeryBigInt_T)                :: C1, C2
    type(VeryBigInt_T)                :: LeadingZeros ! Just store "trailing" 0s
    logical, parameter                :: verbose = .false.
    ! Executable
    call Create (0, c)
    if ( isZero ( a ) .or. isZero( b ) ) then
      return
    else if ( isOne ( a ) ) then
      c = b
    else if ( isOne ( b ) ) then
      c = a
    endif
    call Create (0, c1)
    call Create (0, c2)
    call Create (0, LeadingZeros)
    do k=1, len_trim(b%cval)
      call ReadNumsFromChars ( b%cval(k:k), ib )
      call intTimesVeryBigInt ( ib, a, c2 )
      if ( len_trim(LeadingZeros%cval) > 0 ) then
        c2%cval = &
          & trim(LeadingZeros%cval) // c2%cval
        LeadingZeros%cval = trim(LeadingZeros%cval) // '0'
      else
        LeadingZeros%cval = '0'
      endif
      if ( verbose ) then
        print *, 'ib: ', ib
        print *, 'c1%cval: ', reverse_trim(c1%cval)
        print *, 'c2%cval: ', reverse_trim(c2%cval)
        print *, 'LeadingZeros%cval: ', reverse_trim(LeadingZeros%cval)
      endif
      call Add ( c2, c1, c )
      c1 = c
    enddo
    ! The product is negative only if the factors are of opposite signs
    c%negative = ( a%negative .neqv. b%negative )
  end subroutine Multiply_very

  ! --------------- Power_Composite --------------------
  ! Set the Composite number "c" to the "n"th power
  ! If n < 1, returns "1"
  ! if n == 1, returns C
  ! otherwise just multiply each power by n
  subroutine Power_Composite( n, C )
    ! Args
    integer, intent(in)                  :: n
    type(CompositeNum_T), intent(inout)  :: C
    ! Executable
    if (  isOne ( c ) ) return ! In case it is 1
    if ( n < 1 ) then
      ! Return 1
      call Destroy_composite( C )
    else if ( n > 1 ) then
      c%powers = n*c%powers
    end if
  end subroutine Power_Composite

  ! --------------- Power_Rational --------------------
  ! Set the Rational number "R" to the "n"th power
  ! where n can be any integer, even < 1
  subroutine Power_Rational( n, R )
    ! Args
    integer, intent(in)                  :: n
    type(RationalNum_T), intent(inout)   :: R
    ! Executable
    if ( n < 0 ) then
      ! R^(-k) = (1/R)^k
      call Invert ( R )
      call Power( n, R%numerator )
      call Power( n, R%denominator )
    else if ( n == 0 ) then
      call Destroy( R ) ! anything to the "0" is "1"
    else if ( n > 1 ) then
      call Power( n, R%numerator )
      call Power( n, R%denominator )
    end if
  end subroutine Power_Rational

  ! --------------- Power_Very --------------------
  ! Set the VeryBigInt "c" to the "n"th power
  ! If n < 1, returns "1"
  ! if n == 1, returns C
  ! otherwise just multiply c by itself as many times a n requires
  subroutine Power_Very( n, C )
    ! Args
    integer, intent(in)                  :: n
    type(VeryBigInt_T), intent(inout)    :: C
    ! Local variables
    type(VeryBigInt_T)                   :: T
    ! Executable
    if (  isOne ( c ) .or. isZero( c ) ) return ! In case it is 1
    if ( mod(n, 2) == 0 ) c%negative = .false.
    if ( n < 1 ) then
      c%cval = '1'
    elseif ( n > 1 ) then
      t = c
      call c_to_the_nth ( n )
    end if
  contains
    recursive subroutine c_to_the_nth ( n )
      integer, intent(in)                  :: n
      type(VeryBigInt_T)                   :: Q ! A temp for c
      ! print *, 'In c_to_the_nth with n=', n
      ! call Dump( t )
      ! call Dump( c )
      q = c
      call Multiply( t, q, c )
      ! print *, '*** ', reverse_trim(c%cval)
      if ( n < 3 ) return
      call c_to_the_nth ( n-1 )
    end subroutine c_to_the_nth
  end subroutine Power_Very

  ! --------------- Reduce_Composite --------------------
  ! Divide two composite nums; c = a / b
  ! knowing each factor in b is also in a
  subroutine Reduce_Composite( A, B, C )
    ! Args
    type(CompositeNum_T), intent(in)  :: A
    type(CompositeNum_T), intent(in)  :: B
    type(CompositeNum_T), intent(out) :: C
    ! Local variables
    integer                           :: i
    integer                           :: j
    integer                           :: k
    integer                           :: nFactors
    integer, dimension(:), pointer    :: factors
    integer, dimension(:), pointer    :: powers
    ! Executable
    if ( isEqual ( A, B ) ) return
    call Copy_composite( A, C )
    c%factors = 0
    nFactors = size(a%factors)
    k = 0 ! This is index into c's arrays
    do i=1, nFactors
      j = FindFirst( b%factors, a%factors(i) )
      if ( j < 1 ) then
        ! Not a factor in common, so a keeps its full power
        k = k + 1
        c%factors(k) = a%factors(i)
        c%powers(k)  = a%powers(i)
      else if ( a%powers(i) > b%powers(j) ) then
        ! A factor in common, so a suffers a reduction
        k = k + 1
        c%factors(k) = a%factors(i)
        c%powers(k)  = a%powers(i) - b%powers(j)
      else ! if ( a%powers(i) == b%powers(j) )
        ! Both factor and power in common, so it is eliminated totally
        ! no operation
      end if
    end do
    nFactors = FindLast( c%factors /= 0 )
    nullify ( factors, powers )
    call allocate_test( factors, nFactors, 'factors', &
      & ModuleName // 'ReduceCompositeNum' )
    call allocate_test( powers, nFactors, 'powers', &
      & ModuleName // 'ReduceCompositeNum' )
    factors = c%factors(1:nFactors)
    powers  = c%powers(1:nFactors)
    call Destroy_composite( c )
    c = CompositeNum_T( factors, powers )
    ! call deallocate_test( factors, 'factors', &
    !  & ModuleName // 'ReduceCompositeNum' )
    ! call deallocate_test( powers, 'powers', &
    !  & ModuleName // 'ReduceCompositeNum' )
  end subroutine Reduce_Composite

  ! --------------- Reduce_Rational --------------------
  ! Reduce a rational number by dividing numerator and denominator
  ! by their gcd
  subroutine Reduce_Rational( R )
    ! Args
    type(RationalNum_T), intent(inout) :: R
    ! Local variables
    type(CompositeNum_T)              :: gcd
    type(CompositeNum_T)              :: C ! R%numerator / gcd
    type(CompositeNum_T)              :: D ! R%denominator / gcd
    ! Executable
    call GreatestCommonDivisor( R%numerator, R%denominator, gcd )
    if ( isOne ( gcd ) ) return
    call Reduce( R%numerator, gcd, C )
    call Reduce( R%denominator, gcd, D )
    call Destroy( R )
    R = RationalNum_T( C, D )
  end subroutine Reduce_Rational

  ! --------- IsPrime ------------------
  logical function IsPrime_int( n, primes )
    ! Return TRUE if arg is prime, FALSE if not
    ! Method:
    ! if n is < M (largest of stored array) primenumbers,
    ! just check if n is an element of the array
    ! Otherwise, check if arg is divisible by any of them
    ! Obvious bug:
    ! If n is too large (i.e., n > M^2) it might be divisible
    ! by a prime too large to be in the array
    ! Args
    integer, intent(in) :: n
    integer, optional, dimension(:), intent(in) :: primes ! instead of primenumbrs
    integer, parameter                          :: nBins = 2000
    ! Executable
    isPrime_int = .false.
    if ( n < 2 ) return
    if ( present(primes) ) then
      call myPoorCode( n, primes )
    else
      call myPoorCode( n, primenumbers )
    endif
  contains
    subroutine myPoorCode( n, primenumbers )
    ! Args
    integer, intent(in) :: n
    integer, dimension(:), intent(in) :: primenumbers
    ! Internal variables
    integer :: kM    ! p_k^2 < n for k = 1 .. kM
    integer :: M     ! size(primenumbers)
    integer :: sqrtn
    ! Due to poor coding practices, we must initialize the primenumbers array
    ! by attempting to access it via the prime function
    kM = prime(1)
    M = size(primenumbers)
    ! print *, 'maxval(primenumbers): ', maxval(primenumbers)
    isPrime_int = .false.
    if ( n < primenumbers(M)+1 ) then
      ! kM = FindFirst( primenumbers > n )
      kM = FindFirstHere( primenumbers, n, nBins )
      if ( kM < 2 ) return
      isPrime_int = any(n == primenumbers(1:kM-1))
      return
    end if
    sqrtn = sqrt(n * 1.0)
    ! kM = FindFirst( primenumbers > sqrtn )
    kM = FindFirstHere( primenumbers, sqrtn, nBins )
    if ( kM < 1 ) kM = size(primenumbers)
    ! print *, 'kM: ', kM
    isPrime_int = all( mod(n, primenumbers(:kM)) > 0 )
    end subroutine myPoorCode

    integer function FindFirstHere ( array, probe, nBins )
    ! Find the first array element > probe
    ! Exploit the fact that the arrays is monotonically increaing
    ! by subdividing it into nBins
    ! and choosing the bin over which to search
    !
    ! Note: Someday we should look into adding nBins as an optional arg
    ! to the FindFirst function in the MLSFinds module
    ! Args
    integer, dimension(:), intent(in) :: array
    integer, intent(in)               :: probe
    integer, intent(in)               :: nBins
    ! Internal variables
    integer                           :: BinSize
    integer                           :: iBin
    integer                           :: k1
    integer                           :: k2
    integer                           :: M
    ! Executable
    M = size(array)
    ! Don't let bins become too small
    BinSize =  M/nBins
    if ( BinSize < 5 ) then
      FindFirstHere = FindFirst( array > probe )
    else
      FindFirstHere = 0
      k2 = 0
      do iBin=1, nBins
        k1 = k2 + 1
        k2 = min(M, k2+BinSize)
        ! Have we reached the proper bin yet?
        if ( array(k2) > probe ) then
          ! Apply an offset stemming from the array section k1:k2
          FindFirstHere = k1 - 1 + FindFirst( array(k1:k2) > probe )
          return
        endif
      enddo
    endif
    end function FindFirstHere
  end function IsPrime_int

  function NextPrime( n, primes ) result( next )
    ! Returns the next prime number greater than the arg n
    ! Method:
    ! Starting with n, we'll check each integer until we find one that is prime
    ! Args
    integer, intent(in) :: n ! E.g., if n=100 returns 101 which is next prime > 100
    integer, optional, dimension(:), intent(in) :: primes ! instead of primenumbrs
    integer :: next
    ! Executable
    next = 2
    if ( n < 3 ) return
    ! The following trick starts at the next odd integer in case n is even
    ! (Because no even integer > 2 is prime)
    next = n + ( 1 - mod(n, 2) )
    do
      if ( isPrime( next, primes ) ) return
      next = next + 2
    end do
  end function NextPrime

  integer function Prime(n)
    ! Returns the nth prime number
    ! Args
    integer, intent(in) :: n ! E.g., if n=1 returns 2 which is first prime
    prime = -999
    if ( n < 1 .or. n > MAXNUMPRIMES ) return
    if ( primenumbers(1) > 0 ) then
      prime = primenumbers(n)
      return
    end if
    ! Initializing
    ! We need to build array of primenumbers
    call appendValues(primenumbers, '    2     3     5     7    11    13    17    19    23    29')
    call appendValues(primenumbers, '   31    37    41    43    47    53    59    61    67    71')
    call appendValues(primenumbers, '   73    79    83    89    97   101   103   107   109   113')
    call appendValues(primenumbers, '  127   131   137   139   149   151   157   163   167   173')
    call appendValues(primenumbers, '  179   181   191   193   197   199   211   223   227   229')
    call appendValues(primenumbers, '  233   239   241   251   257   263   269   271   277   281')
    call appendValues(primenumbers, '  283   293   307   311   313   317   331   337   347   349')
    call appendValues(primenumbers, '  353   359   367   373   379   383   389   397   401   409')
    call appendValues(primenumbers, '  419   421   431   433   439   443   449   457   461   463')
    call appendValues(primenumbers, '  467   479   487   491   499   503   509   521   523   541')
    call appendValues(primenumbers, '  547   557   563   569   571   577   587   593   599   601')
    call appendValues(primenumbers, '  607   613   617   619   631   641   643   647   653   659')
    call appendValues(primenumbers, '  661   673   677   683   691   701   709   719   727   733')
    call appendValues(primenumbers, '  739   743   751   757   761   769   773   787   797   809')
    call appendValues(primenumbers, '  811   821   823   827   829   839   853   857   859   863')
    call appendValues(primenumbers, '  877   881   883   887   907   911   919   929   937   941')
    call appendValues(primenumbers, '  947   953   967   971   977   983   991   997  1009  1013')
    call appendValues(primenumbers, ' 1019  1021  1031  1033  1039  1049  1051  1061  1063  1069')
    call appendValues(primenumbers, ' 1087  1091  1093  1097  1103  1109  1117  1123  1129  1151')
    call appendValues(primenumbers, ' 1153  1163  1171  1181  1187  1193  1201  1213  1217  1223')
    call appendValues(primenumbers, ' 1229  1231  1237  1249  1259  1277  1279  1283  1289  1291')
    call appendValues(primenumbers, ' 1297  1301  1303  1307  1319  1321  1327  1361  1367  1373')
    call appendValues(primenumbers, ' 1381  1399  1409  1423  1427  1429  1433  1439  1447  1451')
    call appendValues(primenumbers, ' 1453  1459  1471  1481  1483  1487  1489  1493  1499  1511')
    call appendValues(primenumbers, ' 1523  1531  1543  1549  1553  1559  1567  1571  1579  1583')
    call appendValues(primenumbers, ' 1597  1601  1607  1609  1613  1619  1621  1627  1637  1657')
    call appendValues(primenumbers, ' 1663  1667  1669  1693  1697  1699  1709  1721  1723  1733')
    call appendValues(primenumbers, ' 1741  1747  1753  1759  1777  1783  1787  1789  1801  1811')
    call appendValues(primenumbers, ' 1823  1831  1847  1861  1867  1871  1873  1877  1879  1889')
    call appendValues(primenumbers, ' 1901  1907  1913  1931  1933  1949  1951  1973  1979  1987')
    call appendValues(primenumbers, ' 1993  1997  1999  2003  2011  2017  2027  2029  2039  2053')
    call appendValues(primenumbers, ' 2063  2069  2081  2083  2087  2089  2099  2111  2113  2129')
    call appendValues(primenumbers, ' 2131  2137  2141  2143  2153  2161  2179  2203  2207  2213')
    call appendValues(primenumbers, ' 2221  2237  2239  2243  2251  2267  2269  2273  2281  2287')
    call appendValues(primenumbers, ' 2293  2297  2309  2311  2333  2339  2341  2347  2351  2357')
    call appendValues(primenumbers, ' 2371  2377  2381  2383  2389  2393  2399  2411  2417  2423')
    call appendValues(primenumbers, ' 2437  2441  2447  2459  2467  2473  2477  2503  2521  2531')
    call appendValues(primenumbers, ' 2539  2543  2549  2551  2557  2579  2591  2593  2609  2617')
    call appendValues(primenumbers, ' 2621  2633  2647  2657  2659  2663  2671  2677  2683  2687')
    call appendValues(primenumbers, ' 2689  2693  2699  2707  2711  2713  2719  2729  2731  2741')
    call appendValues(primenumbers, ' 2749  2753  2767  2777  2789  2791  2797  2801  2803  2819')
    call appendValues(primenumbers, ' 2833  2837  2843  2851  2857  2861  2879  2887  2897  2903')
    call appendValues(primenumbers, ' 2909  2917  2927  2939  2953  2957  2963  2969  2971  2999')
    call appendValues(primenumbers, ' 3001  3011  3019  3023  3037  3041  3049  3061  3067  3079')
    call appendValues(primenumbers, ' 3083  3089  3109  3119  3121  3137  3163  3167  3169  3181')
    call appendValues(primenumbers, ' 3187  3191  3203  3209  3217  3221  3229  3251  3253  3257')
    call appendValues(primenumbers, ' 3259  3271  3299  3301  3307  3313  3319  3323  3329  3331')
    call appendValues(primenumbers, ' 3343  3347  3359  3361  3371  3373  3389  3391  3407  3413')
    call appendValues(primenumbers, ' 3433  3449  3457  3461  3463  3467  3469  3491  3499  3511')
    call appendValues(primenumbers, ' 3517  3527  3529  3533  3539  3541  3547  3557  3559  3571')
    call appendValues(primenumbers, ' 3581  3583  3593  3607  3613  3617  3623  3631  3637  3643')
    call appendValues(primenumbers, ' 3659  3671  3673  3677  3691  3697  3701  3709  3719  3727')
    call appendValues(primenumbers, ' 3733  3739  3761  3767  3769  3779  3793  3797  3803  3821')
    call appendValues(primenumbers, ' 3823  3833  3847  3851  3853  3863  3877  3881  3889  3907')
    call appendValues(primenumbers, ' 3911  3917  3919  3923  3929  3931  3943  3947  3967  3989')
    call appendValues(primenumbers, ' 4001  4003  4007  4013  4019  4021  4027  4049  4051  4057')
    call appendValues(primenumbers, ' 4073  4079  4091  4093  4099  4111  4127  4129  4133  4139')
    call appendValues(primenumbers, ' 4153  4157  4159  4177  4201  4211  4217  4219  4229  4231')
    call appendValues(primenumbers, ' 4241  4243  4253  4259  4261  4271  4273  4283  4289  4297')
    call appendValues(primenumbers, ' 4327  4337  4339  4349  4357  4363  4373  4391  4397  4409')
    call appendValues(primenumbers, ' 4421  4423  4441  4447  4451  4457  4463  4481  4483  4493')
    call appendValues(primenumbers, ' 4507  4513  4517  4519  4523  4547  4549  4561  4567  4583')
    call appendValues(primenumbers, ' 4591  4597  4603  4621  4637  4639  4643  4649  4651  4657')
    call appendValues(primenumbers, ' 4663  4673  4679  4691  4703  4721  4723  4729  4733  4751')
    call appendValues(primenumbers, ' 4759  4783  4787  4789  4793  4799  4801  4813  4817  4831')
    call appendValues(primenumbers, ' 4861  4871  4877  4889  4903  4909  4919  4931  4933  4937')
    call appendValues(primenumbers, ' 4943  4951  4957  4967  4969  4973  4987  4993  4999  5003')
    call appendValues(primenumbers, ' 5009  5011  5021  5023  5039  5051  5059  5077  5081  5087')
    call appendValues(primenumbers, ' 5099  5101  5107  5113  5119  5147  5153  5167  5171  5179')
    call appendValues(primenumbers, ' 5189  5197  5209  5227  5231  5233  5237  5261  5273  5279')
    call appendValues(primenumbers, ' 5281  5297  5303  5309  5323  5333  5347  5351  5381  5387')
    call appendValues(primenumbers, ' 5393  5399  5407  5413  5417  5419  5431  5437  5441  5443')
    call appendValues(primenumbers, ' 5449  5471  5477  5479  5483  5501  5503  5507  5519  5521')
    call appendValues(primenumbers, ' 5527  5531  5557  5563  5569  5573  5581  5591  5623  5639')
    call appendValues(primenumbers, ' 5641  5647  5651  5653  5657  5659  5669  5683  5689  5693')
    call appendValues(primenumbers, ' 5701  5711  5717  5737  5741  5743  5749  5779  5783  5791')
    call appendValues(primenumbers, ' 5801  5807  5813  5821  5827  5839  5843  5849  5851  5857')
    call appendValues(primenumbers, ' 5861  5867  5869  5879  5881  5897  5903  5923  5927  5939')
    call appendValues(primenumbers, ' 5953  5981  5987  6007  6011  6029  6037  6043  6047  6053')
    call appendValues(primenumbers, ' 6067  6073  6079  6089  6091  6101  6113  6121  6131  6133')
    call appendValues(primenumbers, ' 6143  6151  6163  6173  6197  6199  6203  6211  6217  6221')
    call appendValues(primenumbers, ' 6229  6247  6257  6263  6269  6271  6277  6287  6299  6301')
    call appendValues(primenumbers, ' 6311  6317  6323  6329  6337  6343  6353  6359  6361  6367')
    call appendValues(primenumbers, ' 6373  6379  6389  6397  6421  6427  6449  6451  6469  6473')
    call appendValues(primenumbers, ' 6481  6491  6521  6529  6547  6551  6553  6563  6569  6571')
    call appendValues(primenumbers, ' 6577  6581  6599  6607  6619  6637  6653  6659  6661  6673')
    call appendValues(primenumbers, ' 6679  6689  6691  6701  6703  6709  6719  6733  6737  6761')
    call appendValues(primenumbers, ' 6763  6779  6781  6791  6793  6803  6823  6827  6829  6833')
    call appendValues(primenumbers, ' 6841  6857  6863  6869  6871  6883  6899  6907  6911  6917')
    call appendValues(primenumbers, ' 6947  6949  6959  6961  6967  6971  6977  6983  6991  6997')
    call appendValues(primenumbers, ' 7001  7013  7019  7027  7039  7043  7057  7069  7079  7103')
    call appendValues(primenumbers, ' 7109  7121  7127  7129  7151  7159  7177  7187  7193  7207')
    call appendValues(primenumbers, ' 7211  7213  7219  7229  7237  7243  7247  7253  7283  7297')
    call appendValues(primenumbers, ' 7307  7309  7321  7331  7333  7349  7351  7369  7393  7411')
    call appendValues(primenumbers, ' 7417  7433  7451  7457  7459  7477  7481  7487  7489  7499')
    call appendValues(primenumbers, ' 7507  7517  7523  7529  7537  7541  7547  7549  7559  7561')
    call appendValues(primenumbers, ' 7573  7577  7583  7589  7591  7603  7607  7621  7639  7643')
    call appendValues(primenumbers, ' 7649  7669  7673  7681  7687  7691  7699  7703  7717  7723')
    call appendValues(primenumbers, ' 7727  7741  7753  7757  7759  7789  7793  7817  7823  7829')
    call appendValues(primenumbers, ' 7841  7853  7867  7873  7877  7879  7883  7901  7907  7919')
    call appendValues(primenumbers, ' 7927  7933  7937  7949  7951  7963  7993  8009  8011  8017')
    call appendValues(primenumbers, ' 8039  8053  8059  8069  8081  8087  8089  8093  8101  8111')
    call appendValues(primenumbers, ' 8117  8123  8147  8161  8167  8171  8179  8191  8209  8219')
    call appendValues(primenumbers, ' 8221  8231  8233  8237  8243  8263  8269  8273  8287  8291')
    call appendValues(primenumbers, ' 8293  8297  8311  8317  8329  8353  8363  8369  8377  8387')
    call appendValues(primenumbers, ' 8389  8419  8423  8429  8431  8443  8447  8461  8467  8501')
    call appendValues(primenumbers, ' 8513  8521  8527  8537  8539  8543  8563  8573  8581  8597')
    call appendValues(primenumbers, ' 8599  8609  8623  8627  8629  8641  8647  8663  8669  8677')
    call appendValues(primenumbers, ' 8681  8689  8693  8699  8707  8713  8719  8731  8737  8741')
    call appendValues(primenumbers, ' 8747  8753  8761  8779  8783  8803  8807  8819  8821  8831')
    call appendValues(primenumbers, ' 8837  8839  8849  8861  8863  8867  8887  8893  8923  8929')
    call appendValues(primenumbers, ' 8933  8941  8951  8963  8969  8971  8999  9001  9007  9011')
    call appendValues(primenumbers, ' 9013  9029  9041  9043  9049  9059  9067  9091  9103  9109')
    call appendValues(primenumbers, ' 9127  9133  9137  9151  9157  9161  9173  9181  9187  9199')
    call appendValues(primenumbers, ' 9203  9209  9221  9227  9239  9241  9257  9277  9281  9283')
    call appendValues(primenumbers, ' 9293  9311  9319  9323  9337  9341  9343  9349  9371  9377')
    call appendValues(primenumbers, ' 9391  9397  9403  9413  9419  9421  9431  9433  9437  9439')
    call appendValues(primenumbers, ' 9461  9463  9467  9473  9479  9491  9497  9511  9521  9533')
    call appendValues(primenumbers, ' 9539  9547  9551  9587  9601  9613  9619  9623  9629  9631')
    call appendValues(primenumbers, ' 9643  9649  9661  9677  9679  9689  9697  9719  9721  9733')
    call appendValues(primenumbers, ' 9739  9743  9749  9767  9769  9781  9787  9791  9803  9811')
    call appendValues(primenumbers, ' 9817  9829  9833  9839  9851  9857  9859  9871  9883  9887')
    call appendValues(primenumbers, ' 9901  9907  9923  9929  9931  9941  9949  9967  9973')

    prime = primenumbers(n)
  end function Prime

  function PrimeFactors( n, factors, powers ) result(nFactors)
    ! Break n into its prime factors and, optionally, their powers
    ! E.g., called with 200, returns 2 (the number of prime factors)
    ! along with:
    ! factors = (/ 2, 5 /)
    ! powers  = (/ 3, 2 /)
    ! for 200 = (2^3) (5^2)
    ! Args
    integer, intent(in)                              :: n
    integer, dimension(:), intent(out)               :: factors
    integer, dimension(:), optional, intent(out)     :: powers
    integer                                          :: nFactors
    ! Internal variables
    logical :: addOneTok
    integer :: i ! check prime(i)
    integer :: k ! counts which factor
    integer :: m
    integer :: primo
    ! Executable
    nFactors = 0
    factors = 1
    if ( present(powers) ) powers = 0
    if ( isPrime(n) .or. n < 2 ) then
      factors(1) = n
      if ( present(powers) ) powers(1) = 1
      nFactors = 1
      return
    end if
    m = n
    k = 1
    do i = 1, n
      if ( m < 2 .or. k > size(factors) ) then
        return
      end if
      ! Now check to see if prime(i) divides m
      if ( prime(i) < 0 ) then
        primo = nextPrime(primo+1)
      else
        primo = prime(i)
      end if
      ! print *, 'm, primo', m, primo
      if ( primo < 0 ) return
      addOneTok = .false.
      do
        if ( mod(m, primo) /= 0 .or. m < primo ) exit
        addOneTok = .true.
        factors(k) = primo
        if ( present(powers) ) powers(k) = powers(k) + 1
        nFactors = k
        m = m / primo
      end do
      if ( addOneTok ) k = k + 1
    end do
  end function PrimeFactors

  function PrimeIndex(n) result(i)
    ! Returns the index i of prime number prime(i) greater than the arg n
    ! Method:
    ! Starting with n, we'll check each integer until we find one that is prime
    ! Args
    integer, intent(in) :: n ! E.g., if n=100 returns 101 which is next prime > 100
    integer :: i
    ! Executable
    ! First, some intializing
    if ( primenumbers(1) < 0 ) i = prime(n)
    i = findFirst( primenumbers > n )
  end function PrimeIndex

  ! --------------- Subtract_composite --------------------
  ! Subtract two composite nums; c = a - b
  ! If a < b, returns b - a
  subroutine Subtract_composite( A, B, C )
    ! Args
    type(CompositeNum_T), intent(in)  :: A
    type(CompositeNum_T), intent(in)  :: B
    type(CompositeNum_T), intent(out) :: C
    call Add_composite( A, B, C, sbtrct=.true. )
  end subroutine Subtract_composite

  ! --------------- Subtract_rational --------------------
  ! Subtract two rational nums; c = a - b
  ! If a < b, returns b - a
  subroutine Subtract_rational( AR, BR, CR )
    ! Args
    type(RationalNum_T), intent(in)  :: AR
    type(RationalNum_T), intent(in)  :: BR
    type(RationalNum_T), intent(out) :: CR
    call Add_rational( AR, BR, CR, sbtrct=.true. )
  end subroutine Subtract_rational

  ! --------------- Subtract_very --------------------
  ! Subtract two rational nums; c = a - b
  ! If a < b, returns b - a
  subroutine Subtract_very( A, B, C )
    ! Args
    type(VeryBigInt_T), intent(in)  :: A
    type(VeryBigInt_T), intent(in)  :: B
    type(VeryBigInt_T), intent(out) :: C
    call Add_Very( A, B, C, sbtrct=.true. )
  end subroutine Subtract_very

  ! --------------------- Trim ---------------
  ! Will we need this for the VeryBigReal_T, too?
  function Trim_Very ( arg ) result ( trimmed )
    ! args
    type(VeryBigInt_T), intent(in)          :: arg
    type(VeryBigInt_T)                      :: trimmed
    ! Internal variables
    integer :: k
    ! Executable
    if ( isZero(arg) ) return ! Because each trimmed springs into existence as ' '
    trimmed = arg
    ! Start at the end and work backward, trimming off '0' chars until
    ! we hit a non-zero, at which point we're done.
    do k=len_trim(arg%cval), 1, -1
      if ( trimmed%cval(k:k) == '0' ) then
        trimmed%cval(k:k) = ' '
      else
        return
      endif
    enddo
  end function Trim_Very

  ! ------------------- Private procedures -------------------------
  subroutine AppendValues( array, chars )
    ! Append new values to end of array where new values
    ! are encoded by chars
    ! Args
    integer, dimension(:), intent(inout) :: array
    character(len=*), intent(in)         :: chars
    ! Internal variables
    integer :: error, nBegin, nEnd, nRead
    integer, dimension(10) :: ints
    ! Executable
    nBegin = FindFirst( array, -999 )
    if ( nBegin < 1 ) return
    call readIntsFromList( chars, ints, error )
    nRead = FindLast( ints /= -999 )
    if ( nRead < 1 ) return
    nEnd = min( size(array), nBegin+nRead-1 )
    array(nBegin:nEnd) = ints(:nEnd-nBegin+1)
  end subroutine AppendValues
  
!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MLSNumbers
!=============================================================================

!
! $Log$
! Revision 2.9  2019/10/09 22:33:22  pwagner
! Added two new numeric types and their procedures: VeryBig for ints and reals
!
! Revision 2.8  2019/09/19 16:01:02  pwagner
! Added AreRelativelyPrime function
!
! Revision 2.7  2019/08/08 17:41:09  pwagner
! Added procedures for -Num types to Add, Subtract, islessThan, and isGreaterThan
!
! Revision 2.6  2017/12/06 01:02:58  pwagner
! Added BigPrime, FindAllPrimes; improved performance
!
! Revision 2.5  2017/10/17 23:40:40  pwagner
! Added Fibonacci
!
! Revision 2.4  2017/03/10 00:39:19  vsnyder
! Make components of CompositeNum_T allocatable
!
! Revision 2.3  2016/07/28 01:35:51  vsnyder
! Remove unused variable declaration
!
! Revision 2.2  2016/03/05 00:17:18  pwagner
! Added Factorial, BinaryCoef functions
!
! Revision 2.1  2015/06/11 22:59:05  pwagner
! First commit
!
