! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
module TIME_M
!=============================================================================

! Compute either CPU time, in arbitrary units, or wall-clock time, in
! seconds since midnight.

  use MLSStrings, only: yyyymmdd_to_dai
  implicit NONE
  private

  public :: TIME_NOW, TIME_NOW_D, TIME_NOW_S, USE_WALL_CLOCK, &
   & WAIT, WAIT_LOOP_LIMITS, RETRY, INIT_RETRY, &
   & TRY_AGAIN, RETRY_SUCCESS, TOO_MANY_RETRIES, TIME_DIVISOR

  interface TIME_NOW
    module procedure TIME_NOW_D
    module procedure TIME_NOW_S
  end interface

  interface WAIT
    module procedure WAIT_D
    module procedure WAIT_S
  end interface

  integer, dimension(3), save :: STARTTIME = -1  ! Reset on first call to time_now
  integer, save :: TIME_DIVISOR = 1  ! Divide by this before returning result
  logical, save :: USE_WALL_CLOCK = .false.
  integer, save :: WAIT_LOOP_LIMITS = 100
  integer, save :: TRY_NUMBER
  integer, save :: RETRY_SUCCESSFUL_RESULT, RETRY_FAILED_RESULT
  integer, parameter :: RETRY_SUCCESSFUL_DEFAULT = 0
  integer, parameter :: RETRY_FAILED_DEFAULT = 999
  integer, parameter :: TRY_AGAIN = 1
  integer, parameter :: RETRY_SUCCESS = TRY_AGAIN - 1
  integer, parameter :: TOO_MANY_RETRIES = RETRY_SUCCESS - 1
  real, save    :: INIT_T0

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine TIME_NOW_D ( T )
    integer, parameter :: RK = kind(0.0d0)
    include "time_now.f9h"
  end subroutine TIME_NOW_D

  subroutine TIME_NOW_S ( T )
    integer, parameter :: RK = kind(0.0e0)
    include "time_now.f9h"
  end subroutine TIME_NOW_S

  subroutine WAIT_D ( T, ERR )
    integer, parameter :: RK = kind(0.0d0)
    real(rk), intent(in) :: T
    include "wait.f9h"
  end subroutine WAIT_D

  subroutine WAIT_S ( T, ERR )
    integer, parameter :: RK = kind(0.0e0)
    real(rk), intent(in) :: T
    include "wait.f9h"
  end subroutine WAIT_S

  subroutine INIT_RETRY ( SUCCESSFUL_RESULT, FAILED_RESULT )
    integer, intent(in), optional :: SUCCESSFUL_RESULT
    integer, intent(in), optional :: FAILED_RESULT
    ! As long as RETRY_FAILED_RESULT is default, try for success
    ! If not default, will try for any value != RETRY_FAILED_RESULT
    RETRY_FAILED_RESULT = RETRY_FAILED_DEFAULT  
    if ( present(SUCCESSFUL_RESULT) ) then
      RETRY_SUCCESSFUL_RESULT = SUCCESSFUL_RESULT
    elseif ( present(FAILED_RESULT) ) then
      RETRY_FAILED_RESULT = FAILED_RESULT
    else
      RETRY_SUCCESSFUL_RESULT = RETRY_SUCCESSFUL_DEFAULT
    endif
    
    call time_now (INIT_T0)
    TRY_NUMBER = 0
  end subroutine INIT_RETRY

  function RETRY ( TRIAL_VALUE, DELAY, MAX_RETRIES, MAX_RETRYING_TIME )
    integer, intent(in)           :: TRIAL_VALUE
    integer, intent(in), optional :: MAX_RETRIES
    real, intent(in), optional    :: DELAY
    real, intent(in), optional    :: MAX_RETRYING_TIME
    integer                       :: RETRY
    real                          :: T1
    ! Return one of three values
    ! TRY_AGAIN if TRIAL_VALUE not yet successful (or still failure)
    ! RETRY_SUCCESS if TRIAL_VALUE successful (or no longer failure)
    ! TOO_MANY_RETRIES if MAX_RETRIES or MAX_RETRYING_TIME exceeded
    
    ! Usage:
    ! Assume you want to keep calling home until you get a successful result "0"
    ! making a call each 2 seconds, and giving up after 100 tries
    ! call init_retry(SUCCESSFUL_RESULT=0)
    ! do
    !    call home(result)
    !    shall_i = retry(result, delay=2.0, max_retries=100)
    !    if ( shall_i /= try_again) exit
    ! enddo
    ! if ( shall_i /= RETRY_SUCCESS ) 
    !    call exception_handler ( shall_i )
    ! endif
    
    ! This may be useful if you can't trust the file system to cooperate
    ! well with the hdfeos library: if a gdclose is immediately
    ! followed by a gdopen, the bare gdopen may fail because the file
    ! hasn't been fully released

    ! Looking for a successful result?
    if ( RETRY_FAILED_RESULT == RETRY_FAILED_DEFAULT ) then
    ! Have we succeeded yet?
      if ( TRIAL_VALUE == RETRY_SUCCESSFUL_RESULT ) then
        retry = RETRY_SUCCESS
        return
      endif
    else
    ! or looking to avoid failure
    ! have we avoided it?
      if ( TRIAL_VALUE /= RETRY_FAILED_RESULT ) then
        retry = RETRY_SUCCESS
        return
      endif
    endif

    ! Have we tried too many times?
    if ( present(MAX_RETRIES) ) then
      if ( TRY_NUMBER > MAX_RETRIES ) then
          retry = TOO_MANY_RETRIES
          return
      endif
    endif
    ! For too long?
    if ( present(MAX_RETRYING_TIME) ) then
      call time_now (t1)
      if ( t1-init_t0 > MAX_RETRYING_TIME ) then
          retry = TOO_MANY_RETRIES
          return
      endif
    endif
    call wait (delay)
    TRY_NUMBER = TRY_NUMBER + 1
    retry = TRY_AGAIN
  end function RETRY

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module TIME_M

!$Log$
!Revision 2.4  2003/12/05 00:53:26  pwagner
!Added possible TIME_DIVISOR to scale results from time_now
!
!Revision 2.3  2002/10/08 00:09:15  pwagner
!Added idents to survive zealous Lahey optimizer
!
!Revision 2.2  2002/08/27 23:05:25  pwagner
!Added wait and retry
!
!Revision 2.1  2001/11/09 22:45:30  vsnyder
!Initial commit
!
