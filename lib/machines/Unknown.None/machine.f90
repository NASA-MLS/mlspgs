! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE
  use F90_GC, only: GCOLLECT, NCOLLECTIONS, &
   & DONT_EXPAND, DONT_GC, FULL_FREQUENCY, MAX_RETRIES, SILENT_GC
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

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

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
  ! Retuyrns how many garbage collections have been performed

    MLS_HOWMANY_GC = NCOLLECTIONS()
  end function MLS_HOWMANY_GC

end module MACHINE

! $Log$
! Revision 1.2  2002/01/31 19:16:28  pwagner
! Brought up-to-date with shell_command using library as appropriate; untested
!
! Revision 1.1  2001/01/13 00:29:44  pwagner
! moved to lib/machines/MLSCONFG/machine.f90
