! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE
  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

  public :: GETARG
  interface
    subroutine GETARG ( ARGNUM, ARGVAL )
      integer, intent(in) :: ARGNUM  ! 0 = command name, 1 = first arg, etc.
      character(len=*), intent(out) :: ARGVAL   ! Blank if argnum out-of-range
    end subroutine GETARG
  end interface

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private :: IO_ERROR_

  logical, public, save :: NEVERCRASH = .true. ! Change to false for testing

  private :: NOT_USED_HERE

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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
    if ( present(code) ) print *, 'In ifc MACHINE%CRASH_BURN, CODE =', code
    if ( present(msg) ) print *, 'In ifc MACHINE%CRASH_BURN, MSG = ', trim(msg)
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

  subroutine EXIT_WITH_STATUS ( STATUS )
  ! Exit and return STATUS to the invoking process
    integer, intent(in) :: STATUS
    write(*,*)"Fnord. IFC exiting. Should use status ",status,&
      " but I have not implemented it yet."
    stop
  end subroutine EXIT_WITH_STATUS

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.

    character(len=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(len=*), intent(in), optional :: FILE

    integer :: L
    character(len=127) :: MSG           ! From the Lahey IOSTAT_MSG intrinsic

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      write (*,*) file(:l)
    end if
!    call iostat_msg (iostat, msg)       ! Lahey intrinsic
    write(*,*) "Fnord: IFC message handling not done yet"
    write (*,*) msg(:len_trim(msg))     ! Print the error message
    write (*,*) 'Error status code =', iostat
    return
  end subroutine IO_ERROR_

  subroutine SHELL_COMMAND ( Command, Status, Error )
  ! Submit a character variable to the system as a shell command.

    character(len=*), intent(in) :: Command  ! The command
    integer, intent(out), optional :: Status ! Its status, if the system
                                        !  has such a concept, else zero
    integer, intent(out), optional :: Error  ! Status of the routine to submit
                                        ! the command, if the system has
                                        ! such a concept, else zero

    integer :: MyError, MyStatus
    print*,"system not done yet for IFC"
    !call system ( command, myStatus, myError)
    !if ( present(error) ) error = myError
    !if ( present(status) ) status = myStatus
  end subroutine SHELL_COMMAND

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
! Revision 1.3  2003/04/01 22:53:02  pwagner
! Chucked stub getarg--the v7.0 ifc has getarg
!
! Revision 1.2  2002/05/02 09:49:55  hcp
! Just syncing up
!
! Revision 1.1  2001/09/17 10:45:18  pumphrey
! Added machine.f90 file for Intel ifc fortran
!
! Revision 1.3  2001/07/25 19:36:18  vsnyder
! Added an interface for GETARG
!
! Revision 1.2  2001/05/04 23:25:10  vsnyder
! Added Exit_With_Status routine
!
! Revision 1.1  2001/01/13 00:29:44  pwagner
! moved to lib/machines/MLSCONFG/machine.f90
!
! Revision 1.1  2000/10/19 17:41:17  pwagner
! first commit
!
! Revision 2.1  2000/10/09 22:16:14  vsnyder
! Moved machine.f90 from l2 to lib
!
! Revision 2.0  2000/09/05 18:57:42  ahanzel
! Changing file revision to 2.0.
