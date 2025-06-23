! Copyright 2018, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
module Time_Config_M
!=============================================================================

! Compute or print either CPU time, in arbitrary units, or wall-clock time, in
! seconds since midnight.  And some other time-related actions.

  use Dates_Module, only: Yyyymmdd_To_Dai
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!   datatypes
! SayTime_config_t        say time configuration type
! Time_config_t           time configuration type

!   subroutines and functions
! Time_Now                 Returns time according to time_config: as
!                          (1) Arbitrary units if cpu time
!                          (2) s since midnight if wall clock time
!                          (3) s since 1st call to Time_Now (so 1st returns 0)
! === (end of toc) ===
! === (start of api) ===
! Time_Now ( float t, [int invalues(8)] )
! === (end of api) ===

  public :: Time_Now, Time_Now_d, Time_Now_s
  interface Time_Now
    module procedure Time_Now_d
    module procedure Time_Now_s
  end interface

  ! There are two time configurations:
  ! This is the sayTime configuration type used wheen you call sayTime()
  public :: sayTime_config_t
  type sayTime_config_t
    character(len=1)     :: TimingUnits = 's'
    real                 :: StartT1     = 0.
    logical              :: SayUnits    = .false.
    character(len=16)    :: Preamble    = ''  ! Instead of 'Timing for' at start
    character(len=16)    :: Coda        = ''      ! Print at end of line
  end type sayTime_config_t

  ! This is the time configuration type used when you call Time_Now()
  public :: time_config_t
  type time_config_t
    ! Divide by time_divisor before returning Time_Now or waiting
    integer, dimension(8) :: Starttime                   = -1  ! Reset on first call to Time_Now
    integer               :: Startdaysoff                = 0
    integer               :: Time_divisor                = 1  
    logical               :: Use_wall_clock              = .false.
    integer               :: Wait_time                   = 100 ! How many s to sleep waiting for event
    logical               :: WallClockIsElapsedFromStart = .true. ! return 0 for 1st call
  end type time_config_t
  type(time_config_t), public, save :: time_config
  type(sayTime_config_t), public, save :: sayTime_config

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Time_Now_D ( T, Invalues )
    integer, parameter :: RK = kind(0.0d0)
    integer, dimension(8), intent(in), optional :: INVALUES
    include "time_now.f9h"
  end subroutine Time_Now_D

  subroutine Time_Now_S ( T, Invalues )
    integer, parameter :: RK = kind(0.0e0)
    integer, dimension(8), intent(in), optional :: INVALUES
    include "time_now.f9h"
  end subroutine Time_Now_S

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Time_Config_M

!$Log$
!Revision 2.1  2018/12/11 01:17:33  pwagner
!First commit
!
