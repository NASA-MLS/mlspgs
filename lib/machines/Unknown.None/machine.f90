module MACHINE
  use F90_UNIX_ENV, only: IARGC, NAG_GETARG => GETARG
  implicit none

  character(LEN=2) :: END_LINE = ' ' // char(10)
  character(LEN=1) :: FILSEP = '/'      ! '/' for Unix, '\' for DOS or NT
  integer, parameter :: HP = 0          ! Offset for first argument for GETARG

  interface IO_ERROR; module procedure IO_ERROR_; end interface
  private IO_ERROR_

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

end module MACHINE

! $Log$
! Revision 1.1  2000/10/21 00:16:25  pwagner
! First commit
!
! Revision 1.1  2000/10/19 17:40:52  pwagner
! first commit
!
! Revision 1.2  2000/10/12 22:54:12  vsnyder
! Correct a commented-out line that may get commented-in for another
! computer/os/compiler combination
!
! Revision 1.1  2000/10/12 22:21:11  vsnyder
! Change directory name from NAG to pclinuxNAG
!
! Revision 1.3  2000/10/09 22:15:55  vsnyder
! Moved machine.f90 from l2 to lib
!
! Revision 2.0  2000/09/05 18:58:04  ahanzel
! Changing file revision to 2.0.
