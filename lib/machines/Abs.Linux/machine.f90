! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module MACHINE
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

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine IO_ERROR_ ( MESSAGE, IOSTAT, FILE )
  ! Print MESSAGE and FILE, and then do something reasonable with IOSTAT.

    character(len=*), intent(in) :: MESSAGE
    integer, intent(in) :: IOSTAT
    character(len=*), intent(in), optional :: FILE

    integer :: L
!    character(len=127) :: MSG           ! From the Lahey IOSTAT_MSG intrinsic

    write (*,*) message(:len_trim(message))
    if ( present(file) ) then
      l = len_trim(file)
      write (*,*) file(:l)
    end if
!    call iostat_msg (iostat, msg)       ! Lahey intrinsic
!    write (*,*) msg(:len_trim(msg))     ! Print the error message
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

    myStatus = system(command)
    if ( present(error) ) error = 0
    if ( present(status) ) status = myStatus
  end subroutine SHELL_COMMAND

end module MACHINE

! $Log$
! Revision 1.2  2001/11/20 00:08:47  pwagner
! Small changes to work better--or go farther before failing
!
! Revision 1.1  2001/03/08 00:29:41  pwagner
! First Commit
!
