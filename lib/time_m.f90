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

  use Dates_Module, only: Yyyymmdd_To_Dai
  use Optional_m, only: Default
  use Time_Config_m, only: SayTime_Config, SayTime_Config_T, &
    & Time_Config, Time_Config_T, Time_Now, Time_Now_D, Time_Now_S
  implicit none
  private

! === (start of toc) ===
!     c o n t e n t s
!     - - - - - - - -
!   datatypes
! SayTime_config_t        say time configuration type
! Time_config_t           time configuration type

!   subroutines and functions
! Begin                   announce the time of day at Beginning
! Dump                    dump either of the time of day configurations
! Finish                  announce the time of day and consumed CPU time at Finish
! Ms_to_hms               convert millisec to (hour,minute,sec)

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
! === (end of toc) ===

! === (start of api) ===
! Begin ( char* string )
! Finish ( char* string )
! Ms_to_hms ( int ms, int hrs, int mins, int secs )
! SayTime ( [char* Operation], [real t1], [real t2], [log cumulative] )
! ConfigureSayTime ( [real t1], [char timingUnits], [log reset], &
!     [log showTimingUnits], [char* Preamble], [char* Coda] )
! Set_time_config ( real values(8) )
! Time_now ( float t, [int invalues(8)] )
! Note:
! float means either real or double precision types
! === (end of api) ===
  public :: Begin, dump, Finish, ms_to_hms, sayTime, &
   & configureSayTime, SayTime_Config, set_time_config, time_config, &
   & time_now, time_now_d, time_now_s

  interface Dump
    module procedure Dump_SayTime
    module procedure Dump_TimeConfig
  end interface

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
  subroutine MS_to_HMS (ms, hrs, mins, secs)

    ! convert millisecs to hours, minutes and seconds

    integer, intent (in) :: ms
    integer, intent (out) :: hrs, mins, secs

    hrs = ms / 3600000
    mins = mod (ms, 3600000) / 60000
    secs = mod (mod (ms, 3600000), 60000) / 1000

  end subroutine MS_to_HMS

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
!Revision 2.22  2018/12/11 01:16:53  pwagner
!Subdivided now among Time_config, wait_m, and smaller time_m
!
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
