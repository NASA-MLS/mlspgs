module MACHINE
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

end module MACHINE

! $Log$
! Revision 1.1  2000/10/19 17:41:17  pwagner
! first commit
!
! Revision 2.1  2000/10/09 22:16:14  vsnyder
! Moved machine.f90 from l2 to lib
!
! Revision 2.0  2000/09/05 18:57:42  ahanzel
! Changing file revision to 2.0.
