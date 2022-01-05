! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module MACHINE
  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

   ! According to Sun's own docs, e.g. 806-7987.pdf, Sun supplies
   ! getarg as an intrinsic library routine
  public :: GETARG
!   interface
!     subroutine GETARG ( ARGNUM, ARGVAL )
!       integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
!       character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
!     end subroutine GETARG
!   end interface

  interface getarg
    module procedure getarg_internal
  end interface
  private :: getarg_internal

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private IO_ERROR_
  public :: SHELL_COMMAND
  public :: MLS_DISABLE_AUTOGC, MLS_GC_NOW, MLS_HOWMANY_GC, MLS_CONFIG_GC

  ! private :: USleep
  ! These supply sleep and usleep functions from the c library
  ! which we call as procedures, thowing away their return values
  interface
    subroutine usleep ( i ) bind (c)
      ! 4.3BSD, POSIX.1-2001.  POSIX.1-2001 declares this function obsolete;
      ! use nanosleep(2) instead.  POSIX.1-2008 removes the specification of
      !  usleep().    
      use, intrinsic          :: iso_c_binding
      integer, value      :: i
    end subroutine usleep
    subroutine sleep ( i ) bind (c)
    use, intrinsic          :: iso_c_binding
      integer, value      :: i
    end subroutine sleep
  end interface

  logical, public, save :: NEVERCRASH = .true. ! Change to false for testing


!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
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
    if ( present(code) ) print *, 'In g95 MACHINE%CRASH_BURN, CODE =', code
    if ( present(msg) ) print *, 'In g95 MACHINE%CRASH_BURN, MSG = ', trim(msg)
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
  
  subroutine getarg_internal ( ARGNUM, ARGVAL )
    integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
    character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    call get_command_argument (ARGNUM, ARGVAL)
  end subroutine getarg_internal

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.

    character(len=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(len=*), intent(in), optional :: FILE

    integer :: L

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      write (*,*) file(:l)
    end if
    write (*,*) 'Error status code =', iostat
    return
  end subroutine IO_ERROR_

  subroutine exit_with_status (istatus)
    integer,intent(in)::istatus

    call exit(istatus)
  end subroutine exit_with_status

  subroutine SHELL_COMMAND ( Command, Status, Error )
  ! Submit a character variable to the system as a shell command.
  ! Based on Absoft's own documentation :: SupportLibrary.pdf
  ! Absoft Support Libraries (1991-1998)
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

    print*,"system not done yet for IFC"
    myStatus = 0
    !myStatus = system(command)
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
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module MACHINE

! $Log$
! Revision 1.6  2009/06/23 19:58:52  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 1.5  2009/06/23 18:25:44  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 1.4  2005/06/22 20:26:22  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.3  2005/03/04 21:41:16  pwagner
! Expedient until system command available
!
! Revision 1.2  2004/08/19 00:13:12  pwagner
! Added crash_burn to provoke crash with walkback
!
! Revision 1.1  2003/04/25 23:00:49  pwagner
! FIest commit
!
