! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module MLSRandomNumber              ! Some random number-generating things
  !=============================================================================

  use MLSCommon, only : R8
  use MLSMessageModule, only: MLSMessage,MLSMSG_Error

!------------------------------------------------------------------
!   Random number routines from MATH77 libraries
! ../l2:69% ls *.f
!  drang.f  ranpk1.f  ranpk2.f  srang.f
!  ../l2:71% cat *.f > stuff.sed
!  ../l2:72% sed 's/^[Cc]/\!/' stuff.sed > sed.out
! plus a small amount of subsequent editing
!------------------------------------------------------------------
  implicit none

  private
  public :: srang, drang

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  integer, private, parameter          :: haystack=97
  real, private, dimension(haystack)   :: harvest

!     c o n t e n t s
!     - - - - - - - -

! drang             gauss. distribution: 0 mean, 1 s.d. (double)
! srang             gauss. distribution: 0 mean, 1 s.d. (single)


contains

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
      logical FIRST
      integer DPTR
      logical DGFLAG
!      common/RANCD1/DPTR, DGFLAG
      save  DPTR, DGFLAG, FIRST
      save    R, X, Y
      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
         FIRST = .false.
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
           call random_number(harvest)
         DPTR = M
      endif
            X = harvest(DPTR)
!                              Set Y = random, uniform in [-1., 1.]
            DPTR = DPTR - 1
            if(DPTR .eq. 0) then
                 call random_number(harvest)
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
              call random_number(harvest)
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
         DRANG = (XX-YY)*R
         DGFLAG = .true.
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when DGFLAG is true and DPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      DRANG = TWO*X*Y*R
      DGFLAG=.false.
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
      integer i
      real                ONE, TWO
      parameter(ONE = 1.0E0, TWO = 2.0E0)
      real                R, S, U3, X, XX, Y, YY
      logical FIRST
      integer SPTR
      logical SGFLAG
!      common/RANCS1/SPTR, SGFLAG
      save  SPTR, SGFLAG, FIRST
      save    R, X, Y
      data  FIRST/.true./
!     ------------------------------------------------------------------
      if(FIRST) then
         FIRST = .false.
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
           call random_number(harvest)
         SPTR = M
      endif
            X = harvest(SPTR)
!                              Set Y = random, uniform in [-1., 1.]
            SPTR = SPTR - 1
            if(SPTR .eq. 0) then
                 call random_number(harvest)
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
              call random_number(harvest)
            SPTR = M
         endif
         U3 = harvest(SPTR)
!        Changed -TWO*log(U3) to abs(TWO*log(U3)) because Lahey LF90
!        2.00 on a pentium produced -0.0 for -TWO*log(1.0), then got a
!        floating point exception on sqrt(-0.0).
         if ( S <= 0.0 ) then
           print *, 'Illegal S: X, Y ', X, Y
           stop 
         elseif ( U3 == 0.0 ) then
           print *, 'Illegal U3: harvest ', (i, harvest(i), i=1, SPTR)
           stop
         endif
         R = sqrt(abs(TWO*(log(U3))))/S
!
!                                Compute result as  R*Sin(PHI)
!
         SRANG = (XX-YY)*R
         SGFLAG = .true.
         return
      endif
!     -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!        Come here when SGFLAG is true and SPTR is not 1.
!
!                                Compute result as  R*Cos(PHI)
      SRANG = TWO*X*Y*R
      SGFLAG=.false.
      return
      end function  SRANG

!=============================================================================
end module MLSRandomNumber
!=============================================================================

!
! $Log$
! Revision 2.3  2001/09/24 23:04:11  pwagner
! Now uses intrinsic random_number rather than MATH77 ranpk
!
! Revision 2.2  2001/09/24 17:27:07  pwagner
! Fixed blunder
!
! Revision 2.1  2001/09/24 17:22:08  pwagner
! First commit
!
