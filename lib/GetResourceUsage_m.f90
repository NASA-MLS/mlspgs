! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GetResourceUsage_m

  use ISO_C_Binding, only: C_Int, C_Long

  implicit NONE
  private
  public :: GetPID, GetResourceUsage, rusage_t, TicksPerSecond

  ! These are system dependent, and not specified by the C standard:
  integer, parameter :: suSeconds_T = c_long, Time_T = c_long

  ! I used this to get their values in Scientific Linux 6:
  ! #include <stdio.h>
  ! #include <sys/types.h>
  ! 
  ! void main()
  ! {
  !   __time_t tv_sec;
  !   __suseconds_t tv_usec;
  ! 
  !   printf ( "sizeof(tv_sec) = %d, sizeof(tv_usec) = %d\n", 
  !     sizeof(tv_sec), sizeof(tv_usec) );
  ! 
  ! }
  
  type, bind(c) :: TimeVal
    integer(time_t) :: tv_sec       ! Seconds
    integer(suSeconds_t) :: tv_usec ! Microseconds
  end type TimeVal

  type, bind(c) :: rusage_t
    type(timeval) :: ru_utime       ! user time used
    type(timeval) :: ru_stime       ! system time used
    integer(c_long) :: ru_maxrss    ! maximum resident set size
    integer(c_long) :: ru_ixrss     ! integral shared memory size
    integer(c_long) :: ru_idrss     ! integral unshared data size
    integer(c_long) :: ru_isrss     ! integral unshared stack size
    integer(c_long) :: ru_minflt    ! page reclaims
    integer(c_long) :: ru_majflt    ! page faults
    integer(c_long) :: ru_nswap     ! swaps
    integer(c_long) :: ru_inblock   ! block input operations
    integer(c_long) :: ru_oublock   ! block output operations
    integer(c_long) :: ru_msgsnd    ! messages sent
    integer(c_long) :: ru_msgrcv    ! messages received
    integer(c_long) :: ru_nsignals  ! signals received
    integer(c_long) :: ru_nvcsw     ! voluntary context switches
    integer(c_long) :: ru_nivcsw    ! involuntary context switches
  end type rusage_t

  interface

    integer(c_long) function GetPID() bind(C,name='GetPID')
      import
    end function GetPID

    integer(c_int) function GetResourceUsage ( usage ) &
      & bind(C,name='GetResourceUsage')
      import
      type(rusage_t), intent(out) :: usage
    end function GetResourceUsage

    integer(c_long) function TicksPerSecond ( ) &
      & bind(C,name='TicksPerSecond')
      import
    end function TicksPerSecond

  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!----------------------------------------------------------------------
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, not_used_here ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module GetResourceUsage_m

! $Log$
! Revision 2.2  2014/08/19 23:13:29  vsnyder
! Add GetPID and TicksPerSecond
!
! Revision 2.1  2014/08/19 16:54:42  vsnyder
! Initial commit
!
