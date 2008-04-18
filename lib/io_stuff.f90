! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module IO_STUFF

! Useful stuff for I/O

  implicit NONE

  private
  public :: GET_LUN
  public :: READ_TEXTFILE

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface read_textfile
    module procedure read_textfile_arr, read_textfile_sca
  end interface

contains

! ================================================     GET_LUN     =====

  subroutine GET_LUN ( LUN, MSG )
  ! Find a Fortran logical unit number that's not in use.
    integer, intent(out) :: LUN          ! The logical unit number
    logical, intent(in), optional :: MSG ! Print failure message if absent or .true.
    logical :: EXIST, OPENED             ! Used to inquire about the unit
    do lun = 20, 100
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) return
    end do
    lun = -1
    if ( present(msg) ) then
      if ( .not. msg ) return
    end if
    write(*,*) 'IO_STUFF%GET_LUN-E- Unable to get a logical unit number'
    return
  end subroutine GET_LUN

  !------------------ read_textfile
  ! Notes and limitations:
  ! formatted io
  ! No line should be longer than 512 characters
  ! (To get around that limitation supply optional arg maxLineLen)
  subroutine READ_TEXTFILE_arr ( File, string, nLines, maxLineLen )
  ! read a textfile into string array, one line per element
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), dimension(:), intent(inout) :: string    ! its contents
    integer, optional, intent(out) :: nLines ! num lines read
    integer, optional, intent(in) :: maxLineLen
    ! Internal variables
    integer :: lun
    integer :: recrd
    integer :: status
    character(len=len(string)) :: tempStr
    character(len=8) :: xfmt
    character(len=6) :: xlen
    ! What format do we use for reading each line?
    xfmt = '(a512)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i6)' ) maxLineLen
     if ( len(string) < maxLineLen ) write( xlen, '(i6)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(a' // trim(adjustl(xlen)) // ')'
    endif
    ! Try to read the textfile
    ! Don't change unread elements
    ! string = " " 
    if ( present(nLines) ) nLines = 0
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE_ARR-E- Unable to open textfile'
      return
    endif
    recrd = 0
    do
      read( UNIT=lun, fmt=xfmt, IOSTAT=status ) tempStr
      if ( status /= 0 ) exit
      recrd = recrd + 1
      string(recrd) = tempStr
    enddo
    if ( present(nLines) ) nLines = recrd
    close( UNIT=lun, iostat=status )
  end subroutine READ_TEXTFILE_arr

  subroutine READ_TEXTFILE_sca ( File, string, maxLineLen )
  ! read a textfile into sa single tring
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), intent(out) :: string    ! its contents
    integer, optional, intent(in) :: maxLineLen
    ! Internal variables
    integer :: lun
    integer :: recrd
    integer :: status
    character(len=len(string)) :: tempStr
    character(len=8) :: xfmt
    character(len=6) :: xlen
    ! What format do we use for reading each line?
    xfmt = '(a512)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i6)' ) maxLineLen
     if ( index(xlen, '*') < 1 ) xfmt = '(a' // trim(adjustl(xlen)) // ')'
    endif
    ! Try to read the textfile
    string = " "
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE_SCA-E- Unable to open textfile'
      return
    endif
    read( UNIT=lun, fmt=xfmt, IOSTAT=status ) string
    do
      read( UNIT=lun, fmt=xfmt, IOSTAT=status ) tempStr
      if ( status /= 0 ) exit
      string = trim(string) // achar(13) // tempStr
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine READ_TEXTFILE_sca

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module IO_STUFF

! $Log$
! Revision 2.7  2008/04/18 16:28:26  pwagner
! Now works properly with NAG, Lahey, and Intel
!
! Revision 2.6  2008/03/11 00:09:11  pwagner
! Added read_textfile; should work for more compilers
!
! Revision 2.5  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/05/19 23:00:18  vsnyder
! Add optional MSG argument
!
! Revision 2.3  2002/10/08 00:09:10  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.2  2001/04/26 02:39:11  vsnyder
! Fix up CVS stuff
!
! Revision 2.1  2000/10/11 18:33:24  vsnyder
! Move from lib/cf_parser to lib; insert copyright notice
!
! Revision 2.0  2000/09/05 17:41:50  dcuddy
! Change revision to 2.0
!
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
