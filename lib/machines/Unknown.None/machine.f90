! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE
! This file must be preprocessed through a makefile
! containing sed commands to snip out offending lines
! (delimited by BAD_MATCH) depending upon whether to build 
! (1) a version allowing garbage collection (-gc on the link line); or
! (2) a version not allowing garbage collection
!
! The recipe for getting each of the versions above is
!     if  (1) set BAD_MATCH='no -gc'
! else if (2) set BAD_MATCH='-gc'
!
!     sed "/Start $BAD_MATCH/,/End $BAD_MATCH/ d" machine.f90 \
!       > my_machine.f90
! and then just compile my_machine.f90

!---------- Start -gc section
!           requires -gc among LDOPTS
!   (the following lines automatically deleted for version (2))
  use F90_GC, only: GCOLLECT, NCOLLECTIONS, &
   & DONT_EXPAND, DONT_GC, FULL_FREQUENCY, MAX_RETRIES, SILENT_GC
!---------- End -gc section
  use F90_IOSTAT				! everything; see iostat_msg_NAG
  use F90_UNIX_ENV, only: IARGC, NAG_GETARG => GETARG
  ! Exit and return an integer status to the invoking process
  use F90_UNIX_PROC, only: EXIT_WITH_STATUS => EXIT, SYSTEM
  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

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
    if ( present(code) ) print *, 'In NAG MACHINE%CRASH_BURN, CODE =', code
    if ( present(msg) ) print *, 'In NAG MACHINE%CRASH_BURN, MSG = ', trim(msg)
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

  subroutine GETARG ( ARGNUM, ARGVAL )
    integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
    character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    integer :: STATUS
    call nag_getarg ( argnum, argval, errno = status )
    if ( status /= 0 ) argval = ' '
  end subroutine GETARG

  subroutine SHELL_COMMAND ( Command, Status, Error )
  ! Submit a character variable to the system as a shell command.

    character(len=*), intent(in) :: Command  ! The command
    integer, intent(out), optional :: Status ! Its status, if the system
                                        !  has such a concept, else zero
    integer, intent(out), optional :: Error  ! Status of the routine to submit
                                        ! the command, if the system has
                                        ! such a concept, else zero

    integer :: MyError, MyStatus

    call system ( command, myStatus, myError)
    if ( present(error) ) error = myError
    if ( present(status) ) status = myStatus
  end subroutine SHELL_COMMAND

!---------- Start -gc section
!           requires -gc among LDOPTS
!   (the following lines automatically deleted for version (2))
  subroutine MLS_CONFIG_GC ( EXPAND, FREQUENCY, RETRIES, SILENT )
  ! Configures various parameters affecting garbage collection

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
    if ( Which == 'On' .or. Which == 'ON' &
      & .or. Which == 'on' ) then
      DONT_GC = .false.
    else
      DONT_GC = .true.
    endif
  end subroutine MLS_DISABLE_AUTOGC

  subroutine MLS_GC_NOW
  ! Manually collects garbage when called

    CALL GCOLLECT
  end subroutine MLS_GC_NOW

  integer function MLS_HOWMANY_GC()
  ! Returns how many garbage collections have been performed

    MLS_HOWMANY_GC = NCOLLECTIONS()
  end function MLS_HOWMANY_GC
!---------- End -gc section
!---------- Start no -gc section
!           forbids -gc among LDOPTS
!   (the following lines automatically deleted for version (1))
  ! ----------------------------------------------
  ! The following are merely introduced to satisfy 
  ! NAG call interfaces to f90_gc
  ! Because we're assuming the link statement will lack -gc
  ! we have non-functional substitutes
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
!---------- End no -gc section

  ! ----------------------------------------------  not_used_here  -----
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module MACHINE

! $Log$
! Revision 1.4  2002/02/08 18:31:50  pwagner
! Now stores -gc and no-gc versions in same machine.f90
!
! Revision 1.3  2002/02/05 00:38:55  pwagner
! Added garbage collection features
!
! Revision 1.2  2002/01/31 19:16:28  pwagner
! Brought up-to-date with shell_command using library as appropriate; untested
!
! Revision 1.1  2001/01/13 00:29:44  pwagner
! moved to lib/machines/MLSCONFG/machine.f90
