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
module Time_M
!=============================================================================

! Compute or print either CPU time, in arbitrary units, or wall-clock time, in
! seconds since midnight.  And some other time-related actions.

  use Dates_module, only: YYYYMMDD_to_DAI
  use IO_Stuff, only: Pause_Stdin => Pause
  use Machine, only: USleep
  use Optional_m, only: Default
  use PrintIt_m, only: MLSMSG_Warning, PrintItOut
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!   datatypes
! Retry_config_t          retry configuration type
! SayTime_config_t        say time configuration type
! Time_config_t           time configuration type

!   subroutines and functions
! Begin                   announce the time of day at Beginning
! Dump                    dump either of the time of day configurations
! Finish                  announce the time of day and consumed CPU time at Finish
! Ms_to_hms               convert millisec to (hour,minute,sec)

! Init_retry              Initialize the retry mechanism
! Pause                   wait for user input, either via stdin or InputFileName
! Retry                   Engage retry mechanism that monitors iterations,
!                         exception handling, or synchronizing events
!     Example:
!     Assume you want to keep calling home until you get a successful result "0"
!     making a call each 2 seconds, and giving up after 100 tries
!     call init_retry(SUCCESSFUL_RESULT=0)
!     do
!        call home(result)
!        shall_i = retry(result, delay=2.0, max_retries=100)
!        if ( shall_i /= try_again) exit
!     enddo
!     if ( shall_i /= RETRY_SUCCESS ) 
!        call exception_handler ( shall_i )
!     endif
! SayTime                  print a stepwise and cumulative timing line
! ConfigureSayTime         set a start time, units, etc. to be used to say timings times
! Set_time_config          Explicitly set start time in time_config_t
!                            formatted as array of integer values returned 
!                            by intrinsic date_and_time
!                          (/ year,month,day,t-t_utc,hour,minutes,s,ms /)
!                          E.g., input (/ -1, i=1,8 /) to make next call
!                          to time_now return 0
! Time_now                 Returns time according to time_config: as
!                          (1) Arbitrary units if cpu time
!                          (2) s since midnight if wall clock time
!                          (3) s since 1st call to time_now (so 1st returns 0)
! Wait                     Wait for supplied interval to elapse
! Wait_for_event           Wait for an event to occur
! === (end of toc) ===

! === (start of api) ===
! Begin ( char* string )
! Finish ( char* string )
! Init_retry ( [int successful_result], [int failed_result] )
! Int retry ( int trial_value, [real delay], [int max_retries], &
!   [real max_retrying_time] )
! Ms_to_hms ( int ms, int hrs, int mins, int secs )
! Pause ( char* mesg , [char* InputFileName )] )
! SayTime ( [char* Operation], [real t1], [real t2], [log cumulative] )
! ConfigureSayTime ( [real t1], [char timingUnits], [log reset], &
!     [log showTimingUnits], [char* Preamble], [char* Coda] )
! Set_time_config ( real values(8) )
! Time_now ( float t, [int invalues(8)] )
! Wait ( float t, [int err] )
! Wait_for_event ( integer function event, int id, [int err] )
! Wait_for_event ( logical function event, int id )
! Wait_for_event ( logical function event, int id(:), theId )
!
! Note:
! float means either real or double precision types
! === (end of api) ===
  public :: Begin, dump, Finish, &
   & init_retry, ms_to_hms, retry, Pause, retry_success, &
   & sayTime, configureSayTime, set_time_config, &
   & time_now, time_now_d, time_now_s, &
   & too_many_retries, try_again, &
   & wait, wait_for_event

  interface Dump
    module procedure Dump_SayTime
    module procedure Dump_TimeConfig
  end interface

  interface time_now
    module procedure time_now_d
    module procedure time_now_s
  end interface

  interface wait
    module procedure wait_d
    module procedure wait_s
  end interface

  interface wait_for_event
    module procedure wait_for_event_int
    module procedure wait_for_event_log
    module procedure wait_for_events
  end interface

  ! This is the retry configuration type
  public :: retry_config_t
  type retry_config_t
    integer :: Try_number
    integer :: Successful_result
    integer :: Failed_result
    real    :: Init_t0
  end type retry_config_t

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

  ! This is the time configuration type used when you call time_now()
  public :: time_config_t
  type time_config_t
    ! Divide by time_divisor before returning time_now or waiting
    integer, dimension(8) :: Starttime                   = -1  ! Reset on first call to time_now
    integer               :: Startdaysoff                = 0
    integer               :: Time_divisor                = 1  
    logical               :: Use_wall_clock              = .false.
    integer               :: Wait_time                   = 100 ! How many s to sleep waiting for event
    logical               :: WallClockIsElapsedFromStart = .true. ! return 0 for 1st call
  end type time_config_t

  integer, parameter :: Microspers         = 1000000 ! How many micros in a s
  integer, parameter :: Failed_default     = 999
  integer, parameter :: Try_again          = 1
  integer, parameter :: Retry_success      = TRY_AGAIN - 1
  integer, parameter :: Successful_default = 0
  integer, parameter :: Too_many_retries   = RETRY_SUCCESS - 1
  !                                              so first value is 0.0
  type(retry_config_t), public, save :: retry_config
  type(time_config_t), public, save :: time_config
  type(sayTime_config_t), public, save :: sayTime_config

  double precision, private, save :: Start_CPU_time
  integer         , private, save :: Start_WallClockSeconds = 0

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Begin ( SHOW )
    ! Announce the time of day of starting; set start time
    ! for calculating elapsed cpu when calling Finish
    use HighOutput, only: Output_Date_and_Time
    character(len=*), intent(in) :: SHOW
    double precision :: dt2
    call cpu_time ( start_CPU_time )
    call time_Now( dt2 )
    Start_WallClockSeconds = dt2
    call Output_Date_and_Time ( msg=show )
  end subroutine Begin

  subroutine Dump_SayTime ( config )
    ! Dump the SayTime config
    use HighOutput, only: headLine, outputNamedValue
    type(sayTime_config_t), intent(in) :: config
    call headLine ( 'SayTime Configuration' )
    call outputNamedValue ( 'timing units', config%TimingUnits )
    call outputNamedValue ( 'start t1    ', config%startT1     )
    call outputNamedValue ( 'say units   ', config%sayUnits    )
    call outputNamedValue ( 'preamble    ', config%Preamble    )
    call outputNamedValue ( 'coda        ', config%Coda        )
  end subroutine Dump_SayTime

  subroutine Dump_TimeConfig ( config )
    ! Dump the SayTime config
    use HighOutput, only: headLine, outputNamedValue
    type(time_config_t), intent(in) :: config
    call headLine ( 'Time Configuration' )
    call outputNamedValue ( 'start (date)            ', config%starttime(1:3)              )
    call outputNamedValue ( 'start (time)            ', config%starttime(5:7)              )
    call outputNamedValue ( 'days after 2001-01-01   ', config%startdaysoff                )
    call outputNamedValue ( 'time divisor            ', config%time_divisor                )
    call outputNamedValue ( 'wall clock?             ', config%use_wall_clock              )
    call outputNamedValue ( 'event wait time         ', config%wait_time                   )
    call outputNamedValue ( 'start wall clock from 0?', config%wallClockIsElapsedFromStart )
  end subroutine Dump_TimeConfig

  subroutine Finish ( SHOW )
    ! Announce the time of day and CPU time of Finishing
    use HighOutput, only: Output_Date_and_Time
    character(len=*), intent(in) :: SHOW
    double precision :: Finish_CPU_time
    double precision :: Finish_WallClockSeconds
    call cpu_time ( Finish_CPU_time )
    call time_now ( Finish_WallClockSeconds )
    call Output_Date_and_Time ( msg=show, &
      & CPU_seconds = Finish_CPU_time - start_CPU_time, &
      & WallClock_Seconds = int(Finish_WallClockSeconds-Start_WallClockSeconds) )
  end subroutine Finish

  subroutine Init_retry ( successful_result, failed_result )
    integer, intent(in), optional :: SUCCESSFUL_RESULT
    integer, intent(in), optional :: FAILED_RESULT
    ! As long as retry_config%FAILED_RESULT is default, try for success
    ! If not default, will try for any value != retry_config%FAILED_RESULT
    ! Therefore, there's no reason to input both args
    if ( present(successful_result) .and. present(failed_result) ) then
      call PrintItOut( &
        & "both args supplied in call to init_retry--only FAILED_RESULT effective", &
        & MLSMSG_Warning )
    endif
    retry_config%failed_result = FAILED_DEFAULT  
    retry_config%successful_result = SUCCESSFUL_DEFAULT
    if ( present(successful_result) ) &
      & retry_config%successful_result = SUCCESSFUL_RESULT
    if ( present(failed_result) ) &
      & retry_config%failed_result = failed_result
    
    call time_now (retry_config%init_t0)
    retry_config%try_number = 0
  end subroutine Init_retry

  subroutine MS_to_HMS (ms, hrs, mins, secs)

    ! convert millisecs to hours, minutes and seconds

    integer, intent (in) :: ms
    integer, intent (out) :: hrs, mins, secs

    hrs = ms / 3600000
    mins = mod (ms, 3600000) / 60000
    secs = mod (mod (ms, 3600000), 60000) / 1000

  end subroutine MS_to_HMS

  subroutine Pause ( mesg, InputFileName, Prompts )
    use MLSCommon, only: Filenamelen

    ! Wait to read mesg from InputFileName
    ! or stdin (if InputFileName not present)

    character(len=8), intent (out)                        :: mesg            
    character(len=*), intent (in), optional               :: InputFileName   
    character(len=*), intent (in), dimension(:), optional :: Prompts
    ! Internal variables
    character(len=Filenamelen), save                      :: FileName
    integer                                               :: i
    character(len=80)                                     :: myMesg
    integer                                               :: OldWaitTime
    ! Executable
    OldWaitTime = Time_Config%wait_time
    Time_Config%wait_time = 5
    if ( present(InputFileName) ) then
      print *, '(P a u s e d .. w a i t i n g   f o r: ', trim(InputFileName)
      FileName = InputFileName
      call Wait_then_read_mesg ( myMesg )
    else
      call Pause_stdin( mesg, Prompts )
    endif
    mesg = myMesg
    Time_Config%wait_time = OldWaitTime
  contains
    subroutine Wait_then_read_mesg ( mesg )
      character(len=*), intent(out) :: mesg
      integer                       :: unitnum, status
      call wait_for_event( File_There_Now, 1 )
      open ( newunit=unitnum, form='formatted', &
        & file=trim(FileName), status='old', iostat=status )
      read ( unitnum, '(a80)' ) MyMesg
      mesg = MyMesg
      close ( unitnum )
    end subroutine Wait_then_read_mesg

    function File_There_Now( id ) result ( now )
      ! Returns
      ! 0 if filename exists
      ! 1 otherwise
      integer, intent(in)          :: id ! Ignored
      logical                      :: exist
      integer                      :: now
      ! Use-ing MLSFiles might create a circular dependence
      ! now = mls_exists ( trim(FileName) )
      if ( len_trim(FileName) < 1 ) then
        print *, 'Was FileName ever defined?'
        stop
      endif
      inquire( file=trim(FileName), exist=exist )
      now = merge( 0, 1, exist )
      ! print *, trim(FileName) // ':', exist
    end function File_There_Now
  end subroutine Pause

  function Retry ( Trial_value, delay, max_retries, max_retrying_time )
    integer, intent(in)           :: Trial_value
    integer, intent(in), optional :: Max_retries
    real, intent(in), optional    :: Delay
    real, intent(in), optional    :: Max_retrying_time
    integer                       :: Retry
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
    
    ! For other uses, this can hide the grubby details of an event loop

    ! Looking for a successful result?
    if ( retry_config%failed_result == FAILED_DEFAULT ) then
    ! Have we succeeded yet?
      if ( TRIAL_VALUE == retry_config%successful_result ) then
        retry = RETRY_SUCCESS
        return
      endif
    else
    ! or looking to avoid failure
    ! have we avoided it?
      if ( TRIAL_VALUE /= retry_config%failed_result ) then
        retry = RETRY_SUCCESS
        return
      endif
    endif

    ! Have we tried too many times?
    if ( present(MAX_RETRIES) ) then
      if ( retry_config%try_number > MAX_RETRIES ) then
          retry = TOO_MANY_RETRIES
          return
      endif
    endif
    ! For too long?
    if ( present(MAX_RETRYING_TIME) ) then
      call time_now (t1)
      if ( t1 - retry_config%init_t0 > MAX_RETRYING_TIME ) then
          retry = TOO_MANY_RETRIES
          return
      endif
    endif
    call wait (delay)
    retry_config%try_number = retry_config%try_number + 1
    retry = TRY_AGAIN
  end function Retry

  ! ----------------------------------------------  sayTime  -----
  subroutine sayTime ( Operation, t1, t2, cumulative )
  use Output_m, only: Blanks, NewLine, Output
  ! Print stepwise and cumulative time usage for sequence of operations
  ! Or show time since start of run (if t1 is 0)
  ! Skip printing total time if cumulative is both present and FALSE
    character(len=*), optional, intent(in) :: Operation
    real, optional, intent(in)             :: t1       ! Time at last op
    real, optional, intent(out)            :: t2       ! Time at current op
    logical, optional, intent(in)          :: cumulative ! defaults to TRUE
    ! Internal variables
    double precision                       :: dt1 
    double precision                       :: dt2 
    logical                                :: myCumulative ! local
    integer                                :: timeDivisor
    character(len=*), parameter            :: TIMEFORMSMALL = '(F10.2)'
    character(len=*), parameter            :: TIMEFORMBIG = '(1PE10.2)'
    character(len=len(timeformbig))        :: timeForm
    ! Executable
    myCumulative = Default( cumulative, .true. )
    select case ( sayTime_config%TimingUnits )
    case ( 's' )
      timeDivisor = 1
    case ( 'm' )
      timeDivisor = 60
    case ( 'h' )
      timeDivisor = 3600
    end select
    dt1 = Default( t1, sayTime_config%startT1 )
    call time_now ( dt2 )
    if ( present(t2) ) t2 = dt2
    ! call output( (/ StartT1*1.d0, dt1, dt2 /), advance='yes' )
    if ( dt2/timeDivisor < 0.5d0 .or. dt2/timeDivisor > 99999.99d0 ) then
      TIMEFORM = TIMEFORMBIG
    else
      TIMEFORM = TIMEFORMSMALL
    endif
    if ( len_trim(sayTime_config%Preamble) < 1 ) then
      call output ( "Timing ", advance='no' )
    else
      call output ( trim(sayTime_config%Preamble) // " ", advance='no' )
    endif
    if ( present(operation) ) then
      call output ( "for " // trim(Operation), advance='no'  )
    endif
    call output ( (dt2-dt1)/timeDivisor, FORMAT=TIMEFORM, advance='no' )
    if ( sayTime_config%sayUnits ) &
      & call output( trim(sayTime_config%timingUnits), advance='no' )
    if ( .not. present(t1) ) then
      if ( len_trim(sayTime_config%Coda) > 0 ) &
        & call output( trim(sayTime_config%Coda), advance='no' )
      call newLine
      sayTime_config%startT1 = dt2
      return
    endif
    ! Cumulative
    if ( myCumulative ) then
      call blanks ( 3 )
      call output ( "Total ", advance='no' )
      call output ( (dt2-sayTime_config%startT1)/timeDivisor, FORMAT=TIMEFORM, advance='no' )
      if ( sayTime_config%sayUnits ) &
        & call output( trim(sayTime_config%timingUnits), advance='no' )
    endif
    if ( len_trim(sayTime_config%Coda) > 0 ) &
      & call output( trim(sayTime_config%Coda), advance='no' )
    call newLine
  end subroutine sayTime

  ! ----------------------------------------------  configureSayTime  -----
  subroutine configureSayTime ( t1, units, reset, showTimingUnits, Preamble, Coda )
  ! set start time or timing units to be used by sayTime
    real, optional, intent(in)             :: t1
    character(len=*), optional, intent(in) :: units
    logical, optional, intent(in)          :: reset
    logical, optional, intent(in)          :: showTimingUnits
    character(len=*), optional, intent(in) :: Preamble
    character(len=*), optional, intent(in) :: Coda
    ! Local variables
    logical                                :: myReset
    ! Executable
    myReset = Default( reset, .false. )
    if ( present(t1) ) then
      sayTime_config%startT1 = t1
    elseif ( myReset ) then
      call time_now ( sayTime_config%startT1 )
    endif
    if ( present(units) ) sayTime_config%timingUnits = units
    if ( present(showTimingUnits) ) sayTime_config%sayUnits = showTimingUnits
    if ( present(Preamble) ) sayTime_config%Preamble = Preamble
    if ( present(Coda) ) sayTime_config%Coda = Coda
  end subroutine configureSayTime

  subroutine set_time_config ( VALUES )
    integer, dimension(8), intent(in) :: VALUES
    character (len=8) :: date   ! in yyyymmdd format
    character (len=4) :: yyyy
    character ( len=2) :: mm, dd
    time_config%starttime = values
    write(yyyy, '(i4)') values(1)
    write(mm, '(i2)') values(2)
    write(dd, '(i2)') values(3)
    date = yyyy // mm // dd
    call yyyymmdd_to_dai(date, time_config%startdaysoff)
    ! print *, 'values: ', values
    ! print *, 'yyyy: ', yyyy
    ! print *, 'mm: ', mm
    ! print *, 'dd: ', dd
    ! print *, 'date: ', date
    ! print *, 'startdaysoff: ', startdaysoff
  end subroutine set_time_config

  subroutine TIME_NOW_D ( T, INVALUES )
    integer, parameter :: RK = kind(0.0d0)
    integer, dimension(8), intent(in), optional :: INVALUES
    include "time_now.f9h"
  end subroutine TIME_NOW_D

  subroutine TIME_NOW_S ( T, INVALUES )
    integer, parameter :: RK = kind(0.0e0)
    integer, dimension(8), intent(in), optional :: INVALUES
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

  ! ------------ wait_for_event ---------------
  ! This family waits for an "event" to occur. The occurrence is signaled
  ! by either
  ! (a) an integer-valued function returning 0 for sucess; or
  ! (b) a logical-valued function returning true
  ! The specific event can be picked out by the argument "id"
  subroutine wait_for_event_int ( event, id, err )
    integer, optional, intent(out)     :: err
    integer, intent(in)                :: id
    interface
      integer function event ( id )
        integer, intent(in) :: id
      end function event
    end interface
    do
      call Usleep( int(time_config%wait_time*MICROSPERS/time_config%time_divisor) )
      select case ( event(id) )
      case ( SUCCESSFUL_DEFAULT )
        exit
      case ( FAILED_DEFAULT )
        if ( present(err) ) err = 1
        exit
      ! case default
        ! nothing--just keep looping
      end select
    enddo
  end subroutine wait_for_event_int

  subroutine wait_for_event_log ( event, id )
    integer, intent(in)                :: id
    interface
      logical function event ( id )
        integer, intent(in) :: id
      end function event
    end interface
    do
      call Usleep( int(time_config%wait_time*MICROSPERS/time_config%time_divisor) )
      if ( event(id) ) exit
    enddo
  end subroutine wait_for_event_log

  subroutine wait_for_events ( event, id, theID )
    ! This lets you wait for any one of a possible number of events
    integer, dimension(:), intent(in)   :: id
    integer, intent(out)                :: theID
    interface
      logical function event ( id )
        integer, intent(in) :: id
      end function event
    end interface
    ! Internal variables
    ! integer                             :: i
    ! Executable
    do
      call Usleep( int(time_config%wait_time*MICROSPERS/time_config%time_divisor) )
      idLoop: do theID=1, size(id)
        if ( event(id(theID)) ) exit idLoop
      enddo idLoop
      if ( theID <= size(id) ) exit
    enddo
  end subroutine wait_for_events

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Time_M

!$Log$
!Revision 2.21  2018/10/25 22:43:51  pwagner
!Added Pause command with optional InputFilename arg
!
!Revision 2.20  2017/01/11 23:25:02  pwagner
!May now Dump the two time configs
!
!Revision 2.19  2016/11/15 19:26:13  pwagner
!May write distinctive Coda at end of line when sayTime; can track both CPU and wallClockSeconds
!
!Revision 2.18  2016/11/03 20:54:51  pwagner
!Added sayTime and configureSayTime
!
!Revision 2.17  2016/02/29 19:48:12  pwagner
!Usleep got from machine module instead of being an external
!
!Revision 2.16  2015/07/14 23:12:20  pwagner
!Added family of routines to wait for events
!
!Revision 2.15  2014/12/09 00:26:01  vsnyder
!Add MS_to_HMS
!
!Revision 2.14  2014/01/09 00:24:29  pwagner
!Some procedures formerly in output_m now got from highOutput
!
!Revision 2.13  2013/08/28 00:37:14  pwagner
!Moved more stuff from MLSMessage down to PrintIt module
!
!Revision 2.12  2013/06/12 02:15:56  vsnyder
!Cruft removal
!
!Revision 2.11  2012/04/20 01:27:53  vsnyder
!Add Begin, Finish, some cannonball polishing
!
!Revision 2.10  2009/06/23 18:25:44  pwagner
!Prevent Intel from optimizing ident string away
!
!Revision 2.9  2009/01/12 18:45:06  pwagner
!Improved cofiguration; wait means sleep
!
!Revision 2.8  2005/09/22 23:36:34  pwagner
!time_config and retry_config now hold config settings
!
!Revision 2.7  2005/06/22 17:25:51  pwagner
!Reworded Copyright statement, moved rcs id
!
!Revision 2.6  2004/08/04 23:19:02  pwagner
!Much moved from MLSStrings to MLSStringLists
!
!Revision 2.5  2003/12/07 23:12:14  pwagner
!wall clock is now time elapsed from start; may set start time
!
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
