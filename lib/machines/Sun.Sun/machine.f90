! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE

  ! This file added by Hugh Pumphrey in order to get something or 
  ! other to compile with Sun's own Fortran 95 compiler. 
  ! It is not complete yet: I have not looked into how Sun f95 
  ! handles command line args as I didn't need them. 


  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

   ! According to Sun's own docs, e.g. 806-7987.pdf, Sun supplies
   ! getarg as an intrinsic library routine
  public :: GETARG
  interface
    subroutine GETARG ( ARGNUM, ARGVAL )
      integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
      character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    end subroutine GETARG
  end interface

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private IO_ERROR_
  public :: SHELL_COMMAND
  public :: MLS_DISABLE_AUTOGC, MLS_GC_NOW, MLS_HOWMANY_GC, MLS_CONFIG_GC

  logical, public, save :: NEVERCRASH = .true. ! Change to false for testing

  private :: NOT_USED_HERE

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  ! -------------------------------------------------  CRASH_BURN  -----
  subroutine CRASH_BURN ( CODE, MSG, OPT )
  ! Print CODE and MSG if they're present.
  ! Try to make the processor crash, so it produces a walkback if it
  ! supports such a thing.
  ! DON'T put an optional argument corresponding to OPT.  It may be
  ! referenced on some systems in an attempt to make the processor crash.
  ! STOP as a last resort.
  ! Compile everything with -gline to get a walkback -- even with -O.
  
  ! Crash is provoked by one of available methods:
  ! (a) rewind illegal unit number (default)
  ! (b) acos(2.)
  ! You may select which method to use by setting msg = '!x (rest of message)'
  ! where x is the method
    integer, intent(in), optional :: CODE
    character(len=*), intent(in), optional :: MSG
    integer, intent(in), optional :: OPT(:)
    integer :: I = -huge(0)
    real :: arg
    character(len=*), parameter :: DEFAULTMETHOD = 'b' ! 'a' or 'b'
    character(len=1) :: method
    if ( present(code) ) print *, 'In Sun MACHINE%CRASH_BURN, CODE =', code
    if ( present(msg) ) print *, 'In Sun MACHINE%CRASH_BURN, MSG = ', trim(msg)
    method = DEFAULTMETHOD
    if ( NEVERCRASH ) then
      method = 's' ! Just stop, don't crash nor burn
    elseif ( present(msg) ) then
      if ( msg(1:1) == '!' ) method = adjustl(msg(2:))
    endif
    select case (method)
    case (' ', 'a', 'A')
      rewind i
    case ('b', 'B')
      arg = 1.
      if ( .not. present(opt) ) arg = arg + 1.
      arg = acos(arg)
      if ( arg > 0. ) print *, 'arg ', arg
    case default
      ! When you just want to land safely; e.g., a sips production run
    end select
    stop
  end subroutine CRASH_BURN

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.
    character(LEN=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(LEN=*), intent(in), optional :: FILE

    integer :: L
!   character(LEN=127) :: MSG

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      if ( l /= 0 ) write (*,*) file(:l)
    end if
!   call iostat_msg (iostat, msg)       ! Lahey intrinsic
!   write (*,*) msg(:len_trim(msg))    ! Print the error message
    write (*,*) 'Error status code =', iostat
    return
  end subroutine IO_ERROR_

  subroutine exit_with_status (istatus)
    integer,intent(in)::istatus

    call exit(istatus)
  end subroutine exit_with_status

  subroutine SHELL_COMMAND ( Command, Status, Error )
  ! Submit a character variable to the system as a shell command.
  ! Based on Sun's own documentation :: 806-7985.pdf
  ! Fortran Library Reference July 2001 Revision A
  ! (Based on past adversities, that's no guarantee)

    character(len=*), intent(in) :: Command  ! The command
    integer, intent(out), optional :: Status ! Its status, if the system
                                        !  has such a concept, else zero
    integer, intent(out), optional :: Error  ! Status of the routine to submit
                                        ! the command, if the system has
                                        ! such a concept, else zero

    integer :: MyStatus

    interface
      integer function system ( CH )
        character(len=*), intent(in) :: CH
      end function system
    end interface

    myStatus = system(command)
    if ( present(error) ) error = 0
    if ( present(status) ) status = myStatus
  end subroutine SHELL_COMMAND

  ! ----------------------------------------------
  ! The following are merely introduced to satisfy 
  ! NAG call interfaces to f90_gc
  ! Not yet functional
  subroutine MLS_CONFIG_GC ( EXPAND, FREQUENCY, RETRIES, SILENT )
  ! Configures various parameters affecting garbage collection
   logical dont_EXPAND, SILENT_GC
   integer full_frequency, max_retries
    logical, optional, intent(in) :: expand    ! Autoexpand heap?
    integer, optional, intent(in) :: frequency ! How many incremental colls. betw. fulls
    integer, optional, intent(in) :: retries   ! How many attempts before giving up
    logical, optional, intent(in) :: silent    ! Quash report on each collection?
    if ( present(expand) ) dont_expand = .not. expand
    if ( present(frequency) ) full_frequency = frequency
    if ( present(retries) ) max_retries = retries
    if ( present(silent) ) silent_gc = silent
  end subroutine MLS_CONFIG_GC

  subroutine MLS_DISABLE_AUTOGC ( Which )
  ! Turns automatic garbage collection on/off
    character(len=*), intent(in) :: Which  ! 'On' or 'Off'
    logical dont_gc
    if ( Which == 'On' .or. Which == 'ON' &
      & .or. Which == 'on' ) then
      DONT_GC = .false.
    else
      DONT_GC = .true.
    endif
  end subroutine MLS_DISABLE_AUTOGC

  subroutine MLS_GC_NOW
  ! Manually collects garbage when called
  !    CALL GCOLLECT
  end subroutine MLS_GC_NOW

  integer function MLS_HOWMANY_GC()
  ! Returns how many garbage collections have been performed
    MLS_HOWMANY_GC = 0 ! NCOLLECTIONS()
  end function MLS_HOWMANY_GC

  ! ----------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MACHINE

! $Log$
! Revision 1.4  2002/02/05 00:38:39  pwagner
! Added NAG-like interfaces for garbage collection; not functional (yet)
!
! Revision 1.3  2002/01/31 19:16:28  pwagner
! Brought up-to-date with shell_command using library as appropriate; untested
!
! Revision 1.2  2001/09/13 19:15:24  pwagner
! Added getarg interface (but does it work?)
!
! Revision 1.1  2001/06/28 15:12:45  pumphrey
! Machine.f90 for Sun's f95 compiler added. It is not complete yet. (HCP)
