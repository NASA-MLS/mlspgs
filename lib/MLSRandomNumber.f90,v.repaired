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
module MLSRandomNumber              ! Some random number-generating things
  !=============================================================================

  ! use Machine, only: USleep
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error

  implicit none

  private
  public :: Srang, Drang
  public :: MLS_lottery,  MLS_random_number, MLS_random_seed, MLS_scramble
  
  ! Calculate random numbers via MATH77 routines from ran_pack (if true)
  ! or else via the f95 intrinsic (if false)
  logical, public, save  ::                   MATH77_RAN_PACK = .true.

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  integer, private, parameter          ::      haystack =        97
  real, private, dimension(haystack), save  :: harvest
  integer, private, save               ::      lastclock =       0
  integer, private, save               ::      lastcount =       0
  integer, parameter                   ::      MaxSeed1  =     400

  ! The following recapture the effects of ranpk common blocs
  real, private, save  ::                      Xcursp =          123456789.0e0
  double precision, private, save  ::          Xcurdp =          123456789.0D0
  integer, private, save  ::                   Mode
  integer, private, save  ::                   Dptr =            1
  logical, private, save  ::                   Dgflag =          .false.
  integer, private, save  ::                   Sptr =            1
  logical, private, save  ::                   Sgflag =          .false.
  logical, private, save  ::                   First =           .true.
  logical, private, save  ::                   was_last_time_m77 = .true.

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -

!      settable parameter
! MATH77_RAN_PACK   Whether to use MATH77 routines (if TRUE) or intrinsic f95

!      functions and subroutines
! Drang             gauss. distribution: 0 mean, 1 s.d. (double)
! Srang             gauss. distribution: 0 mean, 1 s.d. (single)
! MLS_Lottery       random integer(s) in supplied range
! MLS_Random_number uniform distribution in interval [0,1]
! MLS_Random_seed   put or set or size seed
! MLS_Scramble      randomly permute the elements of an integer array

! Modified to use f95 intrinsic random_number instead of ranpk1 and ranpk2
! Unfortunately, using NAG f95, this had side-effect of making test case
! results worse than when MATH77 libraries were used.
! Revisiting this issue later would be justified if this module
! became required for regular processing.
! However, its current use is for testing whether diagnostic products
! were correctly coded

! Notes: 
! (1) Discovered why ranpk1 and ranpk2 ignored putting a new seed
! Patched using the Lastcount mechanism (is there something better?)
! (2) The latest NAG results are not obviously worse than ranpk1
! === (end of toc) ===

! === (start of api) ===
! MLS_Lottery ( int array_arg(:), range(2), [log verbose], [log unique] )      
! MLS_Random_number ( int array_arg(:) )
! MLS_Random_seed ( [int ssize], [int pput(:)], [int gget(:), &
!   [int new_seed(:)], [log verbose] )      
! MLS_Scramble ( int array_arg(:), [log verbose] )
! real   Srang( )
! double Drang( )

! === (end of api) ===

contains

!     -------------- mls_lottery -----------------------
    ! Method:                                                 
    ! We'll reap an array of random floats and for each      
    ! (i) multiply by 10,                                     
    ! (ii) take the fractional part                           
    ! (iii) convert to base N                                 
    ! (iv) multiply by N                                      
    ! (v) use the leading coefficient of the fractional part  
    !     i.e., first "decimal"                               
    ! Optionally, insist each element is unique, i.e. not repeated
    subroutine MLS_Lottery( array_arg, range, verbose, unique )
      use Dump_0, only: Dump
      use MLSStringLists, only: Readintsfromlist
      use MLSStrings, only: Writeintasbasen
      ! Formal arguments
      integer, intent(out)   :: array_arg(:)
      integer, intent(in)    :: range(2)
      logical, optional, intent(in) :: verbose
      logical, optional, intent(in) :: unique
      ! Internal variables
      integer, parameter     :: ArgSizeMult = 3 ! How big should this really be?
      real                   :: real_args(ArgSizeMult*size(array_arg))
      integer                :: arg
      integer                :: i
      integer                :: k
      integer, dimension(10) :: ints
      integer                :: sizeArgs
      logical                :: myUnique
      logical                :: myVerbose
      integer                :: N
      integer                :: p
      integer, dimension(630):: new_seed = 0
      character(len=16)      :: strs
      ! Executable
      myVerbose = .false.
      if ( present(verbose) ) myVerbose = verbose
      myUnique = .false.
      if ( present(unique) ) myUnique = unique
      ! We use a longer array if we require the result to feature
      ! unique elements so there is less chance of omitting one
      ! during the reaping process
      if ( myUnique ) then
        sizeArgs = ArgSizeMult*size(array_arg)
      else
        sizeArgs = size(array_arg)
      endif
      N = range(2) - range(1) + 1
      if ( N <= 1 ) then
        if ( myVerbose ) print *, 'N too small, just ', N
        array_arg = range(1)
        return
      endif
      if ( myVerbose ) print *, 'Getting new seed'
      ! call USleep ( FIVE )
      call mls_random_seed( new_seed=new_seed, verbose=verbose )
      call mls_random_seed( gget=new_seed, verbose=verbose )
      if ( myVerbose ) then
        call dump( new_seed, 'new_seed' )
        print *, 'xcurdp ', xcurdp
      endif
      call mls_random_number( real_args(1:sizeArgs) )
      if ( myVerbose ) then
        call dump( real_args(1:sizeArgs), 'real_args' )
        print *, 'Mod ', Mode
      endif
      k = 0
      do i=1, sizeArgs
        real_args(i) = 10*real_args(i)
        real_args(i) = real_args(i) - floor(real_args(i))
        ! Now we don't actually convert a fraction to base N
        ! Instead we keep the leftmost 4 digits by
        ! multiplying by N^4 and taking the integer part
        arg = real_args(i)*N**4
        p = convert( arg )
        if ( k < 1 ) then
          k = k + 1
          array_arg(k) = p
        elseif ( myUnique ) then
          if ( any( p == (/ array_arg(1:k) /) ) ) cycle
          k = k + 1
          array_arg(k) = p
        else
          array_arg(i) = p
        endif
      enddo
      if ( .not. myUnique .or. k >= size(array_arg) ) return
      ! We must "fill in" the rest of array_arg
      do
        if ( k >= size(array_arg) ) return
        if ( (k + 1) == size(array_arg) ) then
          ! We are looking for the integer p in the range 1 .. N 
          ! and yet not already found among array_arg[j] j = 1.. k
          ! Method:
          ! Now we know  the arithmetic series 
          ! Sum(i, i=1 .. k+1) = (k+2) * (k+1)/2
          ! so 
          ! array_arg[1] + array_arg[2] + .. array_arg[k] + array_arg[k+1] =
          ! (k+2) * (k+1)/2
          array_arg(k + 1) = ((k+2) * (k+1))/2 - &
            & Sum(array_arg(1:k))
          return
        endif
        if ( myVerbose ) print *, 'Needing another random reaping'
        call mls_random_number( real_args )
        do i=1, sizeArgs
          real_args(i) = 10*real_args(i)
          real_args(i) = real_args(i) - floor(real_args(i))
          ! Now we don't actually convert a fraction to base N
          ! Instead we keep the leftmost 4 digits by
          ! multiplying by N^4 and taking the integer part
          arg = real_args(i)*N**4
          p = convert( arg )
         if ( any( p == (/ array_arg(1:k) /) ) ) cycle
         k = k + 1
         array_arg(k) = p
        enddo
      enddo
    contains
      function convert ( arg ) result ( p )
        integer, intent(in)          :: arg
        integer                      :: p
        ! Convert to baseN
        call writeIntAsBaseN ( arg, N, strs )
        ! Multiplying by N and then taking the leading coefficient
        ! means the same thing as taking the next-leading one, i.e.
        ! the silver medalist
        call ReadIntsFromList( strs, ints )
        p = max( ints(2), 0 ) + range(1)
      end function convert
    end subroutine mls_lottery

!     -------------- mls_random_number -----------------------
      subroutine MLS_Random_number( array_arg )
      ! Formal arguments
      real, intent(inout) :: array_arg(:)

      if ( MATH77_RAN_PACK ) then
        call SRANUA(array_arg, size(array_arg))
      else
        call random_number(array_arg)
      endif
      was_last_time_m77 = MATH77_RAN_PACK
      return
      end subroutine MLS_Random_number

!     -------------- mls_random_seed -----------------------
      subroutine MLS_Random_Seed( ssize, pput, gget, new_seed, verbose )
      ! Formal arguments
      ! (Esssentially stammerings of random_zeed's)
      integer, optional, intent(out)   :: ssize
      integer, optional, intent(in)    :: pput(:)
      integer, optional, intent(out)   :: gget(:)
      integer, optional, intent(inout) :: new_seed(:)
      logical, optional, intent(in)    :: verbose
      ! Local variables
      integer               ::          count
      real                  ::          seconds
      integer, dimension(8) ::          values
      logical               ::          verboseHere
      verboseHere = .false.
      if ( present(verbose) ) verboseHere = verbose
      ! stop
      
      if ( present(ssize) ) then
        if ( MATH77_RAN_PACK ) then
          call RANSIZ(ssize)
        else
          call random_seed(size=ssize)
        endif
      elseif ( present(pput) ) then
        if ( MATH77_RAN_PACK ) then
          call RANPUT(pput(1:2))
        else
          call random_seed(put=pput(1:))
        endif
      elseif ( present(gget) ) then
        if ( MATH77_RAN_PACK ) then
          call RANGET(gget(1:2))
        else
          call random_seed(get=gget(1:))
        endif
      elseif ( present(new_seed) ) then
        ! We try hard to avoid choosing
        ! the same random seed over and over
        
        ! We tried system_clock; success was mixed
        ! call system_clock( count=count )
        call cpu_time( seconds )
        count = 1000*seconds ! This should be tthe number of milliseconds
        call date_and_time(values=values)
        if ( count == lastclock ) then
          count = lastcount + 1
        else
          count = max( count, lastcount + 1 )
          lastclock = count
        endif
        lastcount = count
        new_seed(1) = mod( count, MaxSeed1 )
        new_seed(2) = MaxSeed1*values(8) + values(1)
        ! new_seed(2) = mod( new_seed(2), MaxSeed1 )
        ! if ( mod(new_seed(2), 2) < 1 ) new_seed(1) = new_seed(1) + 1
        if ( count == 0 ) new_seed(1) = new_seed(2) - 60*values(6)
        ! NAG has sometimes complained when any element of new_seed was 0
        ! (Why does it need an array of length 630, anyway? Intel makes do
        ! with just 2))
        ! if ( size(new_seed) > 2 ) then
        !  new_seed(3:) = max( new_seed(3:), new_seed(1) )
        ! endif
        if ( verboseHere ) then
          print *, 'count, values ', count, values(1), values(8)
          print *, 'new_seed(1,2) ', new_seed(1:2)
        endif
        if ( MATH77_RAN_PACK ) then
          call RANPUT(new_seed(1:2))
        else
          call random_seed(put=new_seed(1:))
        endif
      else
        if ( MATH77_RAN_PACK ) then
          call RAN1
        else
          call random_seed()
        endif
      endif
      return
      end subroutine MLS_Random_seed

!     -------------- mls_scramble -----------------------
    ! Scramble array_arg "randomly"
    ! Method:                                                 
    ! Relying on mls_lottery, create array of random ints lots[i] in range 1 .. N
    ! outputs[i] = inputs[lots[i]], i=1 ..
    subroutine MLS_Scramble( array_arg, verbose )
      ! Formal arguments
      integer, intent(inout) :: array_arg(:)
      logical, optional, intent(in) :: verbose
      ! Internal variables
      integer                :: lots(size(array_arg))
      integer                :: temp_args(size(array_arg))
      integer                :: ind
      integer                :: i
      integer                :: N
      ! Executable
      N = size(array_arg)
      if ( N <= 1 )  return
      call mls_lottery( lots, (/1, N/), verbose, unique=.true. )
      temp_args = array_arg
      do i=1, N
       array_arg(i) = temp_args(lots(i))
      enddo
    end subroutine MLS_Scramble

!------------------------------------------------------------------
!   Random number routines from MATH77 libraries
! ../l2:69% ls *.f
!  drang.f  ranpk1.f  ranpk2.f  srang.f
!  ../l2:71% cat *.f > stuff.sed
!  ../l2:72% sed 's/^[Cc]/\!/' stuff.sed > sed.out
! plus a small amount of subsequent editing
!------------------------------------------------------------------

!     ---------------------------- drang ---------------------------
      double precision function  DRANG ()
!>> 1996-04-16 DRANG WVS SQRT(abs(TWO*log(U3))) avoids sqrt(-0.0)
!>> 1994-10-20 DRANG Krogh  Changes to use M77CON
!>> 1994-06-24 DRANG CLL Changed common to use RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-16 DRANG CLL
!>> 1991-11-26 DRANG CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 DRANG CLL Added call to RAN0, and DGFLAG in common.
!>> 1991-01-15 DRANG CLL Reordered common contents for efficiency.
!>> 1990-01-23 DRANG CLL Making names in common same in all subprogams.
!>> 1987-06-09 DRANG CLLawson  Initial code.
!
!     Returns one pseudorandom number from the Gausian (Normal)
!     distribution with zero mean and unit standard deviation.
!     Method taken from Algorithm 334, Comm. ACM, July 1968, p. 498.
!     Implemented at JPL in Univac Fortran V in 1969 by Wiley R. Bunton
!     of JPL and Stephen L. Richie of Heliodyne Corp.
!
!     Adapted to Fortran 77 for the MATH 77 library by C. L. Lawson and
!     S. Y. Chiu, JPL, April 1987, 6/9/87.
!     ------------------------------------------------------------------
!--D replaces "?": ?RANG, ?RANUA, RANC?1, RANC?2, ?PTR, ?NUMS, ?GFLAG
!     RANCD1 and RANCD2 are common blocks.
!     Uses intrinsic functions, log and sqrt.
!     Calls DRANUA to obtain an array of uniform random numbers.
!     Calls RAN0 to initialize DPTR and DGFLAG.
!     ------------------------------------------------------------------
!                        Important common variables
!
!     DPTR [integer] and DGFLAG [logical]
!
!          Will be set to DPTR = 1 and DGFLAG = .false.
!          when RAN0 is called from this subr if this
!          is the first call to RAN0.  Also set to these values when
!          RAN1 or RANPUT is called by a user.
!
!          DGFLAG will be set true on return from this subr when the
!          algorithm has values that it can save to reduce the amount
!          of computation needed the next time this function is
!          referenced.  Will be set false in the contrary case.
!
!          DPTR is the index of the last location used in the
!          common array harvest().  This index is counted down.
!
!     harvest() [floating point]  Buffer of previously computed uniform
!          random numbers.
!     ------------------------------------------------------------------
      integer M
      parameter(M = haystack)
      double precision    ONE, TWO
      parameter(ONE = 1.0D0, TWO = 2.0D0)
      double precision    R, S, U3, X, XX, Y, YY
!      logical FIRST
!      integer DPTR
!      logical DGFLAG
!      common/RANCD1/DPTR, DGFLAG
!      save  DPTR, DGFLAG, FIRST
      save    R, X, Y
!      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
!         FIRST = .false.
         FIRST = MATH77_RAN_PACK
         DPTR = 1
         DGFLAG = .true.
      endif
!
      if (.not. DGFLAG .or. DPTR .eq. 1) then
!
!     Use the Von Neuman rejection method for choosing a random point
!     (X,Y) in the unit circle, X**2 + Y**2 .le. 1.0.
!     Then the angle Theta = arctan(Y/X) is random, uniform in
!     (-Pi/2, Pi/2), and Phi = 2*Theta is random, uniform in (-Pi, Pi).
!     Define S = X**2 + Y**2, then
!     sin(Theta) = Y/sqrt(S),    cos(Theta) = X/sqrt(S),
!     sin(Phi) = 2*X*Y/S,    and cos(Phi) = (X**2 - Y**2)/S.
!
     do
!   10    continue
!                              Set X = random, uniform in [0., 1.]
      DPTR = DPTR - 1
      if(DPTR .eq. 0) then
           call mls_random_number(harvest)
         DPTR = M
      endif
            X = harvest(DPTR)
!                              Set Y = random, uniform in [-1., 1.]
            DPTR = DPTR - 1
            if(DPTR .eq. 0) then
                 call mls_random_number(harvest)
               DPTR = M
            endif
            Y = TWO*harvest(DPTR) - ONE
!
            XX=X*X
            YY=Y*Y
            S=XX+YY
         if(S <= ONE) exit
!         if(S .gt. ONE) go to 10
      enddo
!
!     Choose R randomly from Chi-squared distribution and
!     normalize with S.
!
!                              Set U3 = random, uniform in [0., 1.]
         DPTR = DPTR - 1
         if(DPTR .eq. 0) then
              call mls_random_number(harvest)
            DPTR = M
         endif
         U3 = harvest(DPTR)
!        Changed -TWO*log(U3) to abs(TWO*log(U3)) because Lahey LF90
!        2.00 on a pentium produced -0.0 for -TWO*log(1.0), then got a
!        floating point exception on sqrt(-0.0).
         R = sqrt(abs(TWO*(log(U3))))/S
!
!                                Compute result as  R*Sin(PHI)
!
!         DRANG = (XX-YY)*R
         DRANG = (X-Y)*(X+Y)*R
         DGFLAG = .true.
!         was_last_time_m77 = MATH77_RAN_PACK  ! Already done in mls_random_number
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when DGFLAG is true and DPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      DRANG = TWO*X*Y*R
      DGFLAG=.false.
!         was_last_time_m77 = MATH77_RAN_PACK  ! Already done in mls_random_number
      return
      end function  DRANG

      real             function  SRANG ()
!>> 1996-04-16 SRANG WVS SQRT(abs(TWO*log(U3))) avoids sqrt(-0.0)
!>> 1994-10-20 SRANG Krogh  Changes to use M77CON
!>> 1994-06-24 SRANG CLL Changed common to use RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-16 SRANG CLL
!>> 1991-11-26 SRANG CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 SRANG CLL Added call to RAN0, and SGFLAG in common.
!>> 1991-01-15 SRANG CLL Reordered common contents for efficiency.
!>> 1990-01-23 SRANG CLL Making names in common same in all subprogams.
!>> 1987-06-09 SRANG CLLawson  Initial code.
!
!     Returns one pseudorandom number from the Gausian (Normal)
!     distribution with zero mean and unit standard deviation.
!     Method taken from Algorithm 334, Comm. ACM, July 1968, p. 498.
!     Implemented at JPL in Univac Fortran V in 1969 by Wiley R. Bunton
!     of JPL and Stephen L. Richie of Heliodyne Corp.
!
!     Adapted to Fortran 77 for the MATH 77 library by C. L. Lawson and
!     S. Y. Chiu, JPL, April 1987, 6/9/87.
!     ------------------------------------------------------------------
!--S replaces "?": ?RANG, ?RANUA, RANC?1, RANC?2, ?PTR, ?NUMS, ?GFLAG
!     RANCS1 and RANCS2 are common blocks.
!     Uses intrinsic functions, log and sqrt.
!     Calls SRANUA to obtain an array of uniform random numbers.
!     Calls RAN0 to initialize SPTR and SGFLAG.
!     ------------------------------------------------------------------
!                        Important common variables
!
!     SPTR [integer] and SGFLAG [logical]
!
!          Will be set to SPTR = 1 and SGFLAG = .false.
!          when RAN0 is called from this subr if this
!          is the first call to RAN0.  Also set to these values when
!          RAN1 or RANPUT is called by a user.
!
!          SGFLAG will be set true on return from this subr when the
!          algorithm has values that it can save to reduce the amount
!          of computation needed the next time this function is
!          referenced.  Will be set false in the contrary case.
!
!          SPTR is the index of the last location used in the
!          common array harvest().  This index is counted down.
!
!     harvest() [floating point]  Buffer of previously computed uniform
!          random numbers.
!     ------------------------------------------------------------------
      integer M
      parameter(M = haystack)
!      integer i
      real                ONE, TWO
      parameter(ONE = 1.0E0, TWO = 2.0E0)
      real                R, S, U3, X, XX, Y, YY
!      logical FIRST
!      integer SPTR
!      logical SGFLAG
!      common/RANCS1/SPTR, SGFLAG
!      save  SPTR, SGFLAG, FIRST
      save    R, X, Y
!      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
!         FIRST = .false.
         FIRST = MATH77_RAN_PACK
         SPTR = 1
         SGFLAG = .false.
      endif
!
      if (.not. SGFLAG .or. SPTR .eq. 1) then
!
!     Use the Von Neuman rejection method for choosing a random point
!     (X,Y) in the unit circle, X**2 + Y**2 .le. 1.0.
!     Then the angle Theta = arctan(Y/X) is random, uniform in
!     (-Pi/2, Pi/2), and Phi = 2*Theta is random, uniform in (-Pi, Pi).
!     Define S = X**2 + Y**2, then
!     sin(Theta) = Y/sqrt(S),    cos(Theta) = X/sqrt(S),
!     sin(Phi) = 2*X*Y/S,    and cos(Phi) = (X**2 - Y**2)/S.
!
   do
!9910    continue
!                              Set X = random, uniform in [0., 1.]
      SPTR = SPTR - 1
      if(SPTR .eq. 0) then
           call mls_random_number(harvest)
         SPTR = M
      endif
            X = harvest(SPTR)
!                              Set Y = random, uniform in [-1., 1.]
            SPTR = SPTR - 1
            if(SPTR .eq. 0) then
                 call mls_random_number(harvest)
               SPTR = M
            endif
            Y = TWO*harvest(SPTR) - ONE
!
            XX=X*X
            YY=Y*Y
            S=XX+YY
         if(S <= ONE) exit
!         if(S .gt. ONE) go to 9910
    enddo
!
!     Choose R randomly from Chi-squared distribution and
!     normalize with S.
!
!                              Set U3 = random, uniform in [0., 1.]
         SPTR = SPTR - 1
         if(SPTR .eq. 0) then
              call mls_random_number(harvest)
            SPTR = M
         endif
         U3 = harvest(SPTR)
!        Changed -TWO*log(U3) to abs(TWO*log(U3)) because Lahey LF90
!        2.00 on a pentium produced -0.0 for -TWO*log(1.0), then got a
!        floating point exception on sqrt(-0.0).
         if ( S <= 0.0 ) then
          call MLSMessage ( MLSMSG_Error, ModuleName // '%srang', &
          & "Illegal S: X, Y " )
         elseif ( U3 == 0.0 ) then
          call MLSMessage ( MLSMSG_Error, ModuleName // '%srang', &
          & "Illegal U3, harvest " )
         endif
         R = sqrt(abs(TWO*(log(U3))))/S
!
!                                Compute result as  R*Sin(PHI)
!
!         SRANG = (XX-YY)*R
         SRANG = (X-Y)*(X+Y)*R
!         was_last_time_m77 = MATH77_RAN_PACK  ! Already done in mls_random_number
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when SGFLAG is true and SPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      SRANG = TWO*X*Y*R
      SGFLAG=.false.
!         was_last_time_m77 = MATH77_RAN_PACK  ! Already done in mls_random_number
      return
      end function  SRANG
      
! ------------------ Private Procedures --------------------------

!>> 1997-12-17 RANPK2 Krogh Removed unreferenced labels
!>> 1995-11-22 RANPK2 Krogh Removed multiple entries.
!>> 1992-03-17 CLL Moved SAVE stmt ahead of DATA stmts.
!>> 1992-03-09 CLL Removed "save FIRST" because FIRST is not used.
!>> 1992-03-02 CLL Fix error: Set MODE = 1 in data stmt.
!>> 1991-11-21 CLL Add MODE with values 1, 2, 3, & 4.
!>> 1989-09-11 CLL Multiversion file. RANPK2 or RANPK3
!>> 1987-05-05 RANPK2 Lawson  Initial code.
!
!        This program unit, along with RANPK1, supports random number
!     generation.
!
!     The functionality of this random number package was modeled on the
!     specifications of the random number routines RANDOM and RANDOMSEED
!     in the February, 1987 working draft of the Fortran 8x language
!     standard.  This functionality remains similar to the later draft,
!     Fortran 90, S8.115, June 1990, in which the routine names have
!     been changed to RANDOM_MUMBER and RANDOM_SEED.  This should
!     facilitate replacement of use of this package by Fortran intrinsic
!     when and if Fortran 90 compilers come into widespread use.
!
!     The library user may call RANSIZ or RANGET in this prog unit,
!     or RAN1 or RANPUT in RANPK1 to obtain the functionality of
!     RANDOM_SEED of Fortran 90.  This relates to setting or fetching
!     the seed.
!        Entries RN1 and RNPUT in this prog unit should not be called by
!     library users.  They are intended only to be called from RANPK1.
!        Entry RN2 returns the value of MODE.  This is a convenience in
!     case one is interested in knowing the value of MODE the package
!     has selected.
!        Entries SRANUA and DRANUA (s.p. and d.p. respectively)
!     generate arrays of pseudorandom numbers from the uniform
!     distribution in the range [0., 1.].  These may be called by users
!     and are called by other library routines.
!        Entries SRANUS and DRANUA (s.p. and d.p. respectively)
!     generate arrays of pseudorandom numbers from the uniform
!     distribution in the range [0., 1.] and then transformed to
!     A + B*U.  These are intended for direct use by users.
!     ------------------------------------------------------------------
!              Algorithm for generating pseudorandom numbers
!                     from the uniform distribution.
!
!  The current integer value in the random integer sequence is XCUR,
!  and the next is defined mathematically by the statement:
!
!                 XCUR = mod(AFAC * XCUR,  MDIV)
!  where
!                 MDIV = m = 6_87194_76503 = 2**36 - 233
!  and
!                 AFAC = a = 612_662 = (approx.) 0.58 * 2**20
!
!  XCUR may be any integer value in the range from 1 to m-1, and all
!  integer values in this range will be generated before the sequence
!  cycles.
!
!  We call the above computational algorithm for XCUR the "short"
!  algorithm.  There is also a "long" algorithm that produces exactly
!  the same sequence of values of XCUR:
!                  Q = aint(XCUR/B)
!                  R = XCUR - Q * B
!                  XCUR = AFAC * R - C * Q
!                  do while(XCUR .lt. 0.0)
!                     XCUR = XCUR + MDIV
!                  end do
!  where B and C are constants related to MDIV and AFAC by
!            MDIV = B * AFAC + C
!  We use B = 112165 and C = 243273.  The average number of executions
!  of the statement XCUR = XCUR + MDIV is 1.09 and the maximum number of
!  executions is 3.
!
!  The largest number that must be handled in the "short" algorithm
!  is the product of AFAC  with the max value of XCUR, i.e.,
!    612_662 * 6_87194_76502 = 42_10181_19126_68324 ~= 0.58 * 2**56.
!  Thus this algorithm requires arithmetic exact to at least 56 bits.
!
!  The largest number that must be handled in the "long" algorithm
!  is the product of C with the max value of aint(XCUR/B), i.e.,
!               243273 * 612664 ~= 0.14904e12 ~= 0.54 * 2**38.
!  Thus this algorithm requires arithmetic exact to at least 38 bits.
!
!  To accommodate different compiler/computer systems this program unit
!  contains code for 3 different ways of computing the new XCUR from the
!  old XCUR, each producing exactly the same sequence of of XCUR.
!
!  Initially we have MODE = 1.  When MODE = 1 the code does tests to
!  see which of three implementation methods will be used, and sets
!  MODE = 2, 3, or 4 to indicate the choice.
!
!  Mode 2 will be used in machines such as the Cray that have at
!  least a 38 bit significand in SP arithmetic.  XCUR will be advanced
!  using the "long" algorithm in SP arithmetic.
!
!  Mode 3 will be used on machines that don't meet the Mode 2 test,
!  but can maintain at least a 56 bits exactly in computing
!  mod(AFAC*XCUR, MDIV) in DP arithmetic.  This includes VAX, UNISYS,
!  IBM 30xx, and some IEEE machines that have clever compilers that
!  keep an extended precision representation of the product AFAC*XCUR
!  within the math processor for use in the division by MDIV.  XCUR will
!  be advanced using the "short" algorithm in DP arithmetic.
!
!  Mode 4 will be used on machines that don't meet the Mode 2 or 3
!  tests, but have at least a 38 bit significand in DP arithmetic.
!  This includes IEEE machines that have not-so-clever compilers.
!  XCUR is advanced using the "long" algorithm in DP arithmetic.
!  ---------------------------------------------------------------------
!               Properties of the generated sequence.
!
!        This m is one of the prime numbers cited in
!     Table 1, Page 390, of Knuth, Seminumerical Algorithms, Second
!     edition, 1981, Addison-Wesley.
!     The prime factorization of m-1 is
!           m-1 = p1 * p2 * p3 = 2 * 43801 * 784_451
!     The complementary factors are
!           q(1) = 3_43597_38251, q(2) = 15_68902, and q(3) = 87602.
!
!     The value a is a primitive root of m as is verified by
!     computing a**q(i) mod m for i = 1,3, and finding these values are
!     not 1.  These values are m-1, 2_49653_21011, and 1_44431_31136.
!     The fact that a is a primitive root of m assures that the period
!     of the generator is m-1, i.e. starting with any integer from 1
!     through m-1, all integers in this range will be produced.
!
!     The value a has relatively large values of the measures nu and mu
!     computed for the Spectral Test as described in Knuth, pp. 89-105.
!        (Log10(nu(i)), i=2,6) = 5.4, 3.6,  2.6,  2.2,  1.8
!        (mu(i), i=2,6)        = 3.0, 3.05, 3.39, 4.55, 6.01
!     This assures that the generated sequence will have relatively low
!     autocorrelation.
!     ------------------------------------------------------------------
!                          Alternative algorithm
!
!  An alternative set of constants that has been used widely in
!  commercial and public domain software packages is
!               m = 21474_83647 = 2**31 - 1
!               a = 16807       = 7**5 = (approx.) 0.513 * 2**15
!
!  The largest product that must be handled exactly is approximately
!  0.513 * 2**46 which is approximately  0.36E14.  This is within the
!  double-precision capability of most computer systems.
!
!  The sequence can be started with any integer from 1 through m-1
!  and will generate all integers in this range.  The autocorrelation
!  properties of the whole sequence will not be as good as with the
!  larger values for m and a.
!     ------------------------------------------------------------------
!     C. L. Lawson, F. T. Krogh & S. Chiu, JPL, July 1986, Apr 1987.
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
      subroutine SRANUA(USP, N)
      integer N
      real USP(N)
!                    Common Block
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
!      logical FIRST
!      save FIRST
!      data FIRST / .true. /
!
      if (FIRST .or. .not. was_last_time_m77) then
         FIRST = .false.
         call RANMOD
      end if
!      go to (310, 320, 330, 340), MODE
    select case (MODE)
!  310 stop'In file RANPK2, subroutine SRANUA -- Ivalid value for MODE'
!                                         Mode 2.
!  310  continue
!    case 1 
!          call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "invalid value for MODE in subroutine SRANUA" )
!  320   continue
    case (2)
      do 325 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  322    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 322
         end if
         USP(I) =  XCURSP/MDIVSP
  325 continue
!      go to 350
!                                         Mode 3.
!  330   continue
    case (3)
      do 335 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         USP(I) = real( XCURDP)/MDIVSP
  335 continue
!      go to 350
!                                         Mode 4.
!  340   continue
    case (4)
      do 345 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  342    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 342
         end if
         USP(I) = real( XCURDP)/MDIVSP
  345 continue
!  350 continue
    case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "invalid value for MODE in subroutine SRANUA" )
    end select
      return
      end subroutine SRANUA
!     ------------------------------------------------------------------
      subroutine DRANUA(UDP, N)
      integer N
      double precision UDP(N)
!                    Common Block
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I
      double precision AFACDP, BFACDP, C2DP, MDIVDP, QDP, RDP
      real             AFACSP, BFACSP, C2SP, MDIVSP, QSP, RSP
      parameter( AFACDP=612662.0D0, BFACDP=112165.0d0, &
     &   C2DP=243273.0d0, MDIVDP=68719476503.0D0 )
      parameter( AFACSP=612662.0e0, BFACSP=112165.0e0, &
     &   C2SP=243273.0e0, MDIVSP=68719476503.0e0 )
!      logical FIRST
!      save FIRST
!      data FIRST / .true. /
!
      if (FIRST .or. .not. was_last_time_m77) then
         FIRST = .false.
         call RANMOD
      end if
!      go to (410, 420, 430, 440), MODE
    select case (MODE)
!  410 stop'In file RANPK2, subroutine DRANUA -- Ivalid value for MODE'
!                                         Mode 2.
!  410   continue
!          call MLSMessage ( MLSMSG_Error, ModuleName, &
!          & "invalid value for MODE in subroutine DRANUA" )
!  420   continue
    case (2)
      do 425 I = 1,N
         QSP = aint(XCURSP / BFACSP)
         RSP = XCURSP - QSP * BFACSP
         XCURSP = AFACSP * RSP - C2SP * QSP
  422    if(XCURSP .lt. 0.0e0) then
            XCURSP = XCURSP + MDIVSP
            go to 422
         end if
         UDP(I) = dble(XCURSP) / MDIVDP
  425 continue
!      go to 450
!                                         Mode 3.
!  430   continue
    case (3)
      do 435 I = 1,N
         XCURDP = mod(AFACDP * XCURDP,  MDIVDP)
         UDP(I) =  XCURDP/MDIVDP
  435 continue
!      go to 450
!                                         Mode 4.
!  440   continue
    case (4)
      do 445 I = 1,N
         QDP = aint(XCURDP / BFACDP)
         RDP = XCURDP - QDP * BFACDP
         XCURDP = AFACDP * RDP - C2DP * QDP
  442    if(XCURDP .lt. 0.0d0) then
            XCURDP = XCURDP + MDIVDP
            go to 442
         end if
         UDP(I) =  XCURDP/MDIVDP
  445 continue
!  450 continue
    case default
          call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "invalid value for MODE in subroutine DRANUA" )
    end select
      return
      end subroutine DRANUA

!     ------------------------------------------------------------------
      subroutine RANMOD
!
!     Do tests to decide whether to set the MODE = 2, 3, or 4.
!     Outcome will depend on precision of SP and DP floating point
!     arithmetic on the host system.
!                    Common Block
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!
      integer I, J, K
      double precision TEST(0:3,0:1), DIFF
      double precision QDP, RDP, XTSTDP
      double precision AFACDP, BFACDP, C2DP, MDIVDP, X1DP
      real             QSP, RSP, XTSTSP
      real             AFACSP, BFACSP, C2SP, MDIVSP, X1SP
      parameter( AFACDP=612662.0D0, MDIVDP=68719476503.0D0, &
     & X1DP=123456789.0D0, BFACDP=112165.0d0, C2DP=243273.0d0)
      parameter( AFACSP=612662.0e0, MDIVSP=68719476503.0e0, &
     &  X1SP=123456789.0e0, BFACSP=112165.0e0, C2SP=243273.0e0)
      logical DONE
      save DONE
!
      data DONE / .false. /
      data TEST / 24997965550.0d0, 68719476502.0d0, &
     &            68718863841.0d0, 36962132774.0d0, &
     &            43721510953.0d0,           1.0d0, &
     &                 612662.0d0, 31757343729.0d0/
!
      if (DONE) return
      DONE = .true.
      do 880 MODE = 2, 4
         do 870 J = 0,1
            XTSTDP = TEST(0,J)
            XTSTSP = real(XTSTDP)
            do 860 I = 1,3
!               go to (820, 830, 840),MODE-1
                select case (MODE)
!                                                Test of MODE 2
!  820          continue
                case (2)
                  QSP = aint(XTSTSP / BFACSP)
                  RSP = XTSTSP - QSP * BFACSP
                  XTSTSP = AFACSP * RSP - C2SP * QSP
                  do 822 K = 1,3
                     if(XTSTSP .ge. 0.0e0) go to 825
                       XTSTSP = XTSTSP + MDIVSP
  822               continue
  825               continue
                  DIFF = dble(XTSTSP) - TEST(I,J)
!               go to 850
!                                                Test of MODE 3
!  830          continue
                case (3)
                  XTSTDP = mod(AFACDP * XTSTDP,  MDIVDP)
                  DIFF = XTSTDP - TEST(I,J)
!               go to 850
!                                                Test of MODE 4
!  840          continue
                case (4)
                  QDP = aint(XTSTDP / BFACDP)
                  RDP = XTSTDP - QDP * BFACDP
                  XTSTDP = AFACDP * RDP - C2DP * QDP
                  do 842 K = 1,3
                     if(XTSTDP .ge. 0.0d0) go to 845
                     XTSTDP = XTSTDP + MDIVDP
  842             continue
  845             continue
                  DIFF = XTSTDP - TEST(I,J)
                case default
                      call MLSMessage ( MLSMSG_Error, ModuleName, &
                      & "invalid value for MODE in subroutine RANMOD" )
                end select
!  850          continue
               if(DIFF .ne. 0.0d0) go to 880
!                            Following line ends I loop.
  860         continue
!                            Following line ends J loop.
  870      continue
!
!        Here the computations using the current value of MODE have
!        passed all tests, so we accept this value of MODE.
         XCURDP = X1DP
         XCURSP = X1SP
         return
!
!                            Following line ends MODE loop.
  880 continue
!        The computations were unsuccessful for all values of MODE.
!        This means this random number package will not work on the
!        current host system.  ****** Fatal Error Stop ******
!
!      call ERMSG('RANPK2',1, 2,
!     *'This random no. code will not work on this computer system.','.')
          call MLSMessage ( MLSMSG_Error, ModuleName // '%ranmod', &
          & "This rnc will not work on this computer system " )
      end subroutine RANMOD

!     ------------------------------------------------------------------
!                         For use by library users: CALL RANPUT(KSEED)
      subroutine RANPUT(KSEED)
      integer KSEED(:)     ! Switch to assumed-shape array
!      integer KSEED(*)
!      integer DPTR, SPTR
!      logical DGFLAG, SGFLAG
!      common/RANCD1/DPTR, DGFLAG
!      common/RANCS1/SPTR, SGFLAG
!      save  /RANCD1/, /RANCS1/
!
      call   RNPUT(KSEED)
      DPTR = 1
      SPTR = 1
      DGFLAG = .false.
      SGFLAG = .false.
      return
      end subroutine RANPUT

!     ------------------------------------------------------------------
!                         Entered using CALL RNPUT(KSEED)
!              This entry should not be called by general users.
!              User should call RANPUT(KSEED) in RANPK1.
!
      subroutine RNPUT(KSEED)
      integer KSEED(:)     ! Switch to assumed-shape array
!      integer KSEED(2)
!                    Common Block
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!
      double precision MDIVDP, SCALDP
      real             MDIVSP, SCALSP
      parameter( MDIVDP=68719476503.0D0, MDIVSP=68719476503.0e0 )
      parameter( SCALDP=100000.0D0, SCALSP=100000.0e0 )
!      logical FIRST
!      save FIRST
!      data FIRST / .true. /
!
      if (FIRST .or. .not. was_last_time_m77) then
         FIRST = .false.
         call RANMOD
      end if
      if(MODE .eq. 2) then
         XCURSP = SCALSP * real(abs(KSEED(1))) + real(abs(KSEED(2)))
         XCURSP = max(1.e0, min(XCURSP, MDIVSP - 1.0e0))
      else
!                            Here for MODE = 3 or 4
         XCURDP = SCALDP * dble(abs(KSEED(1))) + dble(abs(KSEED(2)))
         XCURDP = max(1.D0, min(XCURDP, MDIVDP - 1.0D0))
      end if
      return
      end subroutine RNPUT
!     ------------------------------------------------------------------
      subroutine RANGET(KSEED)
      integer KSEED(:)     ! Switch to assumed-shape array
!      integer KSEED(2)
!                    Common Block
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!
      double precision SCALDP
      real             SCALSP
      parameter( SCALDP=100000.0D0, SCALSP=100000.0e0 )
!      logical FIRST
!      save FIRST
!      data FIRST / .true. /
!
      if (FIRST .or. .not. was_last_time_m77) then
         FIRST = .false.
         call RANMOD
      end if
      if(MODE .eq. 2) then
         KSEED(1) = int(XCURSP/SCALSP)
         KSEED(2) = int(XCURSP - SCALSP * real(KSEED(1)))
      else
!                            Here for MODE = 3 or 4
!
         KSEED(1) = int(XCURDP/SCALDP)
         KSEED(2) = int(XCURDP - SCALDP * dble(KSEED(1)))
      end if
      return
      end subroutine RANGET

!     ------------------------------------------------------------------
      subroutine RANSIZ(KSIZE)
         integer KSIZE
         KSIZE = 2
         return
         end subroutine RANSIZ

      subroutine RAN1
!                            Program unit: RANPK1
!>> 1995-11-21 RAMPK1 Krogh Removed multiple entries.
!>> 1994-06-24 CLL Reorganized common. Using RANC[D/S]1 & RANC[D/S]2.
!>> 1992-03-13 CLL Fixed error in RAN0
!>> 1991-11-26 CLL Reorganized common. Using RANCM[A/D/S].
!>> 1991-11-22 CLL Added Entry RAN0 and common variables SGFLAG,DGFLAG
!>> 1991-01-15 CLL Reordered common contents for efficiency.
!>> 1990-01-23 CLL Corrected type stmt for SNUMS in common.
!>> 1987-04-22 RANPK1 Lawson  Initial code.
!
!        This program unit, RANPK1, along with RANPK2,
!     supports random number generation.
!
!        This prog unit has entries RAN1, RAN0, and RANPUT.
!     The library user can call RAN1 to initialize random number
!     generation at a standard initial seed,
!     or call RANPUT(KSEED) to initialize random number generation
!     at a seed value provided by the integer array, KSEED().
!
!     Other higher level random number subrs call RAN0 on their first
!     time flags to be sure the package is initialized.
!
!     As a result of any of these entries this subroutine will
!     set the pointers in the COMMON arrays to 1, indicating to higher
!     level random number subprograms that these buffer arrays are
!     empty.  It also sets SGFLAG and DGFLAG to .false. to indicate to
!     Gaussian generators that they have no internal saved value.
!
!     The user can determine the appropriate dimension for the array,
!     KSEED() by first calling the entry RANSIZ in prog unit RANPK2.
!
!     The user can retrieve the current seed value by calling entry,
!     RANGET in prog unit RANPK2.  This will be the seed that will be
!     used the next time a batch of random numbers are computed.  This
!     is not necessarily the seed associated with the next number that
!     will be returned.
!     C. L. Lawson, F. T. Krogh, & S. Y. Chiu, JPL, Apr 1987.
!     ------------------------------------------------------------------
!
!      integer DPTR, SPTR
!      logical DGFLAG, SGFLAG
!      common/RANCD1/DPTR, DGFLAG
!      common/RANCS1/SPTR, SGFLAG
!      save  /RANCD1/, /RANCS1/
!     ------------------------------------------------------------------
!                      For use by library users: CALL RAN1
      call RN1
      DPTR = 1
      SPTR = 1
      DGFLAG = .false.
      SGFLAG = .false.
      return
      end subroutine RAN1

      subroutine RN1
!      integer MODE
!      real XCURSP
!      double precision XCURDP
!      common / RANCOM / XCURDP, XCURSP, MODE
!              These same parameters are also defined below in RANMOD.
      double precision X1DP
      real             X1SP
      parameter( X1DP=123456789.0D0, X1SP=123456789.0e0 )
!     ------------------------------------------------------------------
!                         Entered using CALL RN1
!              This entry should not be called by general users.
!              User should call RAN1 in RANPK1.
!
         XCURDP = X1DP
         XCURSP = X1SP
         return
         end subroutine RN1

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

end module MLSRandomNumber
!=============================================================================

!
! $Log$
! Revision 2.15  2019/03/21 23:47:34  pwagner
! Repaired errors in nearly every procedure; choosing a new seed now works properly
!
! Revision 2.14  2019/03/18 22:07:31  pwagner
! Repaired error where a new seed had no effect
!
! Revision 2.13  2014/08/06 23:18:34  vsnyder
! Remove declaration of unused local variable
!
! Revision 2.12  2013/08/01 20:18:38  pwagner
! Added mls_scramble and mls_lottery for randomized permutation and integer sequence
!
! Revision 2.11  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.10  2005/06/22 17:25:50  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.9  2002/10/08 00:09:12  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.8  2002/01/09 23:45:20  pwagner
! Removed vistiges of print statements
!
! Revision 2.7  2001/10/18 23:30:59  pwagner
! Works even if change to/from intrinsic
!
! Revision 2.6  2001/10/17 23:38:32  pwagner
! MATH77_ran_pack now publicly settable
!
! Revision 2.5  2001/10/15 23:49:58  pwagner
! Added back some MATH77 stuff, mls_random number and seed
!
! Revision 2.4  2001/09/27 16:39:24  pwagner
! saves harvest array
!
! Revision 2.3  2001/09/24 23:04:11  pwagner
! Now uses intrinsic random_number rather than MATH77 ranpk
!
! Revision 2.2  2001/09/24 17:27:07  pwagner
! Fixed blunder
!
! Revision 2.1  2001/09/24 17:22:08  pwagner
! First commit
!
