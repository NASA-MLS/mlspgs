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

  public :: SHELL_COMMAND

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine EXIT_WITH_STATUS ( STATUS )
  ! Exit and return STATUS to the invoking process
    integer, intent(in) :: STATUS
    external :: SETRCD
    call setrcd ( status )
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
    call iostat_msg (iostat, msg)       ! Lahey intrinsic
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
! Revision 1.4  2002/01/30 00:23:31  vsnyder
! Add Shell_Command subroutine
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
