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

! Useful, low-level stuff for mostly formatted I/O
  use MLSFinds, only: findFirstCharacter => FindFirst, &
    &                 findFirstSubstring => FindFirst

  implicit none

  private
  public :: get_lun
  public :: get_nLines
  public :: read_stdin
  public :: read_textFile
  public :: truncate_textFile
  public :: write_textFile

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (subroutines and functions)
! get_lun           Find a Fortran logical unit number that's not in use.
! get_nLines        Find how many lines are in a text file
! read_stdin        Read standard input into characters scalar or array
! read_textfile     Read contents of a textfile into characters scalar or array
! truncate_textFile          Delete contents of a text file
! write_textfile    Write characters scalar or array out to a textfile
! === (end of toc) ===

! === (start of api) ===
! get_lun( int lun, [log msg] )
! get_nLines( char* File, int nLines, [int maxLineLen] )
! read_stdin( str string, [int maxLineLen], [int nLines] )
! read_textfile( char* File, str string, [int maxLineLen], [int nLines] )
! write_textfile( char* File, str string, [int maxLineLen], [int nLines] )
! truncate_textFile( str string )
! str can be any of
! character(len=*)                 a scalar character string of any length
! character(len=*), dimension(:)   a 1d character array of any length
! character(len=*), dimension(:,:) a 2d character array of any length
!
! Notes:
! (1) Before reading a textfile into an array, make sure the array is allocated
! and of sufficient size
! (2) You could do this by first using get_nLines to discover the size you need, 
! and then allocating the array
! (3) By design, reading a textfile into a string will not change any unread
! elements (so you can prefill with nulls
! (4) get_nLines returns 0 for an empty file, and -1 if the file doesn't exist
! or can't be opened
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface read_stdin
    module procedure read_stdin_arr, read_stdin_arr2d, read_stdin_sca
  end interface

  interface read_textfile
    module procedure read_textfile_arr, read_textfile_arr2d, read_textfile_sca
  end interface

  interface write_textfile
    module procedure write_textfile_arr, write_textfile_arr2d, write_textfile_sca
  end interface

  ! The only legal unit numbers that files may be assigned
  ! for use by Fortran opens, closes, reads and writes
  integer, parameter :: bottom_unit_num = 1
  integer, parameter :: top_unit_num    = 99

  integer, parameter :: maxStrLen       = 1024

contains

! ================================================     GET_LUN     =====

  subroutine GET_LUN ( LUN, MSG, BOTTOM, TOP )
    ! Find a Fortran logical unit number that's not in use.
    ! Args
    integer, intent(out)          :: LUN    ! The logical unit number
    logical, intent(in), optional :: MSG ! Print failure message? (default: T)
    integer, intent(in), optional :: BOTTOM ! Where to begin
    integer, intent(in), optional :: TOP    ! Where to end
    ! Internal variables
    logical :: EXIST, OPENED             ! Used to inquire about the unit
    integer :: myBottom
    integer :: myTop
    ! Executable
    myBottom = bottom_unit_num
    myTop    = top_unit_num
    if ( present(Bottom) ) myBottom = Bottom
    if ( present(Top) ) myTop = Top
    do lun = myBottom, myTop
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

  !------------------ get_nLines
  ! Notes and limitations:
  ! formatted io
  ! No line should be longer than maxStrLen
  ! Returns 0 if file is of zero length
  ! Returns -1 if file doesn't exist or can't be opened
  ! (To get around that limitation supply optional arg maxLineLen)
  subroutine get_nLines ( File, nLines, maxLineLen )
  ! read a textfile into string array, one line per element
    character(len=*), intent(in)  :: File ! its path and name
    integer, intent(out)          :: nLines ! num lines read
    integer, optional, intent(in) :: maxLineLen
    ! Internal variables
    integer :: lun
    integer :: status
    character(len=1), dimension(maxStrLen) :: nullArray
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! Executable
    nLines = -1 ! The value returned if the file doesn't exist
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i8)' ) maxLineLen
     if ( maxStrLen < maxLineLen ) write( xlen, '(i8)' ) maxStrLen
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    else
     write( xlen, '(i8)' ) maxStrLen
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    endif
    ! Try to read the textfile
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      ! write(*,*) 'IO_STUFF%get_nLines-E- Unable to open textfile ' // &
      !  & trim(File)
      return
    endif
    ! print *, 'xfmt: ', xfmt
    nLines = 0
    do
      status = 0
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) nullArray
500   status = -1
50    if ( status /= 0 ) exit
      nLines = nLines + 1
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine get_nLines

  !------------------ read_stdin
  ! Notes and limitations:
  ! Won't change unread elements (so you can prefill with nulls)
  ! formatted io
  ! No line should be longer than len(string)
  ! (To get around that limitation supply optional arg maxLineLen)
  ! Even then certain compilers impose limitations
  ! E.g., NAG can't read a line longer than 1024 from stdin
  subroutine READ_stdin_arr ( string, maxLineLen, nLines )
  ! read stdin into string array, one line per element
    character(len=*), dimension(:), intent(inout) :: string    ! its contents
    integer, optional, intent(in) :: maxLineLen
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    integer :: lun = 5
    integer :: pos
    integer :: recrd
    integer :: status
    character(len=1), dimension(len(string)) :: cArray
    character(len=1), dimension(len(string)) :: nullArray
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i8)' ) maxLineLen
     if ( len(string) < maxLineLen ) write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    else
     write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    endif
    ! Try to read the stdin
    if ( present(nLines) ) nLines = 0
    recrd = 0
    ! print *, 'xfmt: ', xfmt
    do
      status = 0
      call null_fill_1d( nullArray )
      cArray = string( min(recrd+1, size(string)) )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) nullArray
500   status = -1
50    if ( status /= 0 ) exit
      do pos=1, len(string) - 1
        if ( any(nullArray(pos:pos+1) == achar(0)) ) exit
      enddo
      pos = max(pos, 2)
      cArray(1:pos-1) = nullArray(1:pos-1)
      ! print *, cArray
      recrd = min(recrd+1, size(string))
      string(recrd) = transfer( cArray, string(recrd) )
    enddo
    if ( present(nLines) ) nLines = recrd
  end subroutine READ_stdin_arr

  subroutine READ_stdin_arr2d ( chars, LineLen, nLines )
  ! read stdin into a 2d char array, one line per row
  ! leaving unread elements unchanged
  ! (So you can prefill with nulls)
    character(len=1), dimension(:,:), intent(inout) :: chars    ! its contents
    integer, optional, intent(out) :: LineLen ! max line length read
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    character(len=1), dimension(size(chars,2)) :: cArray
    integer :: lun = 5
    integer :: N ! max line length so far
    integer :: col
    integer :: recrd
    integer :: status
    character(len=12) :: xfmt
    character(len=8) :: xlen
    N = 0
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default
    write( xlen, '(i8)' ) size(chars,2)
    if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    ! Try to read stdin
    if ( present(nLines) ) nLines = 0
    recrd = 1
    ! print *, 'xfmt: ', xfmt
    do
      status = 0
      call null_fill_1d( cArray )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      do col=1, size(chars, 2)
        if( cArray(col) == achar(0) ) then
          N = max(N, col-1)
          exit
        else
          chars(recrd, col) = cArray(col)
        endif
      enddo
      if ( col > size(chars, 2) ) N = size(chars, 2)
      recrd = min(recrd + 1, size(chars, 1))
    enddo
    if ( present(nLines) ) nLines = recrd
    if ( present(LineLen) ) LineLen = N
  end subroutine READ_stdin_arr2d

  subroutine READ_stdin_sca ( string, maxLineLen, nLines )
  ! read stdin into a single string
    character(len=*), intent(inout) :: string    ! its contents
    integer, optional, intent(in) :: maxLineLen
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    character(len=1), dimension(len(string)) :: cArray
    integer :: i
    integer :: lun = 5
    integer :: pos
    integer :: recrd
    integer :: status
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! Executable
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i8)' ) maxLineLen
     if ( len(string) < maxLineLen ) write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    else
     write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    endif
    ! Try to read the stdin
    if ( present(nLines) ) nLines = 0
    recrd = 0
    ! print *, 'xfmt: ', xfmt
    i = 0
    do
      status = 0
      call null_fill_1d( carray )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      ! print *, cArray
      recrd = recrd + 1
      oneLine: do pos=1, len(string) - 1
        if ( any(carray(pos:pos+1) == achar(0)) ) exit oneLine
        i = min(i + 1, len(string))
        string(i:i) = carray(pos)
       enddo oneLine
      i = min(i + 1, len(string))
      string(i:i) = achar(13)
    enddo
    if ( present(nLines) ) nLines = recrd
  end subroutine READ_stdin_sca

  !------------------ read_textfile
  ! Notes and limitations:
  ! Won't change unread elements (so you can prefill with nulls)
  ! formatted io
  ! No line should be longer than len(string)
  ! (To get around that limitation supply optional arg maxLineLen)
  subroutine READ_TEXTFILE_arr ( File, string, maxLineLen, nLines )
  ! read a textfile into string array, one line per element
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), dimension(:), intent(inout) :: string    ! its contents
    integer, optional, intent(in) :: maxLineLen
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    integer :: lun
    integer :: pos
    integer :: recrd
    integer :: status
    character(len=1), dimension(len(string)) :: cArray
    character(len=1), dimension(len(string)) :: nullArray
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i8)' ) maxLineLen
     if ( len(string) < maxLineLen ) write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    else
     write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    endif
    ! Try to read the textfile
    if ( present(nLines) ) nLines = 0
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE_ARR-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    recrd = 0
    ! print *, 'xfmt: ', xfmt
    do
      status = 0
      call null_fill_1d( nullArray )
      cArray = string( min(recrd+1, size(string)) )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) nullArray
500   status = -1
50    if ( status /= 0 ) exit
      do pos=1, len(string) - 1
        if ( any(nullArray(pos:pos+1) == achar(0)) ) exit
      enddo
      pos = max(pos, 2)
      cArray(1:pos-1) = nullArray(1:pos-1)
      ! print *, cArray
      recrd = min(recrd+1, size(string))
      string(recrd) = transfer( cArray, string(recrd) )
    enddo
    if ( present(nLines) ) nLines = recrd
    close( UNIT=lun, iostat=status )
  end subroutine READ_TEXTFILE_arr

  subroutine READ_TEXTFILE_arr2d ( File, chars, LineLen, nLines )
  ! read a textfile into a 2d char array, one line per row
  ! leaving unread elements unchanged
  ! (So you can prefill with nulls)
    character(len=*), intent(in)  :: File ! its path and name
    character(len=1), dimension(:,:), intent(inout) :: chars    ! its contents
    integer, optional, intent(out) :: LineLen ! max line length read
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    character(len=1), dimension(size(chars,2)) :: cArray
    integer :: lun
    integer :: N ! max line length so far
    integer :: col
    integer :: recrd
    integer :: status
    character(len=12) :: xfmt
    character(len=8) :: xlen
    N = 0
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default
    write( xlen, '(i8)' ) size(chars,2)
    if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    ! Try to read the textfile
    if ( present(nLines) ) nLines = 0
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE_ARR2D-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    recrd = 1
    ! print *, 'xfmt: ', xfmt
    do
      status = 0
      call null_fill_1d( cArray )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      do col=1, size(chars, 2)
        if( cArray(col) == achar(0) ) then
          N = max(N, col-1)
          exit
        else
          chars(recrd, col) = cArray(col)
        endif
      enddo
      if ( col > size(chars, 2) ) N = size(chars, 2)
      recrd = min(recrd + 1, size(chars, 1))
    enddo
    if ( present(nLines) ) nLines = recrd
    if ( present(LineLen) ) LineLen = N
    close( UNIT=lun, iostat=status )
  end subroutine READ_TEXTFILE_arr2d

  subroutine READ_TEXTFILE_sca ( File, string, maxLineLen, nLines )
  ! read a textfile into a single string
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), intent(inout) :: string    ! its contents
    integer, optional, intent(in) :: maxLineLen
    integer, optional, intent(out) :: nLines ! num lines read
    ! Internal variables
    character(len=1), dimension(len(string)) :: cArray
    integer :: i
    integer :: lun
    integer :: pos
    integer :: recrd
    integer :: status
    character(len=12) :: xfmt
    character(len=8) :: xlen
    ! Executable
    ! What format do we use for reading each line?
    xfmt = '(128a1)' ! This is the default; if lines are larger supply maxLineLen
    if ( present(maxLineLen) ) then
     write( xlen, '(i8)' ) maxLineLen
     if ( len(string) < maxLineLen ) write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    else
     write( xlen, '(i8)' ) len(string)
     if ( index(xlen, '*') < 1 ) xfmt = '(' // trim(adjustl(xlen)) // 'a1)'
    endif
    ! Try to read the textfile
    if ( present(nLines) ) nLines = 0
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='old', iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%READ_TEXTFILE_ARR-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    recrd = 0
    ! print *, 'xfmt: ', xfmt
    i = 0
    do
      status = 0
      call null_fill_1d( carray )
      read( UNIT=lun, fmt=xfmt, eor=50, end=500, err=50, advance='no' ) cArray
500   status = -1
50    if ( status /= 0 ) exit
      ! print *, cArray
      recrd = recrd + 1
      oneLine: do pos=1, len(string) - 1
        if ( any(carray(pos:pos+1) == achar(0)) ) exit oneLine
        i = min(i + 1, len(string))
        string(i:i) = carray(pos)
       enddo oneLine
      i = min(i + 1, len(string))
      string(i:i) = achar(13)
    enddo
    if ( present(nLines) ) nLines = recrd
  end subroutine READ_TEXTFILE_sca
  
  !------------------ truncate_textFile
  subroutine truncate_textFile( filename )
    character(len=*), intent(in) :: filename
    integer :: unitnum
    call get_lun( unitnum )
    open( unit=unitnum, file=filename, form='formatted', status='replace' )
    close( unitnum )
  end subroutine truncate_textFile

  !------------------ write_textfile
  ! We assume line feeds are already in string
  subroutine write_TEXTFILE_arr ( File, string )
  ! write a string array out to a textfile, one line per element
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), dimension(:), intent(in) :: string    ! its contents
    ! Internal variables
    integer :: i, n
    integer :: lun
    integer :: status
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for writing each line?
    ! Try to write the textfile
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='unknown', access='sequential', &
      & recl=size(string)*len(string(1)) + 1, iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%write_TEXTFILE_ARR-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    do i=1, size(string)
      n = FindFirstSubString( string(i), achar(0) )
      if ( n < 2 ) then
        write ( lun, '(a)', advance='yes' ) ''
      else
        write ( lun, '(a)', advance='yes' ) string(i)(:n-1)
      endif
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine write_TEXTFILE_arr

  subroutine write_TEXTFILE_arr2d ( File, chars )
  ! write a 2-d array out to a textfile, one line per row
    character(len=*), intent(in)  :: File ! its path and name
    character(len=1), dimension(:,:), intent(in) :: chars    ! its contents
    ! Internal variables
    integer :: i, n
    integer :: lun
    integer :: status
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for writeing each line?
    ! Try to write the textfile
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='unknown', access='sequential', &
      & recl=size(chars) + 1, iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%write_TEXTFILE_arr2d-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    do i=1, size(chars,1)
      n = FindFirstCharacter( chars(i,:), achar(0) )
      if ( n < 2 ) then
        write ( lun, '(a)', advance='yes' ) ''
      else
        write ( lun, '(a)', advance='yes' ) chars(i,:n-1)
      endif
    enddo
    close( UNIT=lun, iostat=status )
  end subroutine write_TEXTFILE_arr2d

  subroutine write_TEXTFILE_sca ( File, string )
  ! write a textfile into string array, one line per element
    character(len=*), intent(in)  :: File ! its path and name
    character(len=*), intent(in) :: string    ! its contents
    ! Internal variables
    integer :: lun
    integer :: status
    ! print *, 'Name of textfile: ', trim(File)
    ! What format do we use for writeing each line?
    ! Try to write the textfile
    call GET_LUN ( LUN )
    open(UNIT=lun, form='formatted', &
      & file=trim(File), status='unknown', access='sequential', &
      & recl=len(string) + 1, iostat=status )
    if ( status /= 0 ) then
      write(*,*) 'IO_STUFF%write_TEXTFILE_sca-E- Unable to open textfile ' // &
        & trim(File)
      return
    endif
    write ( lun, '(a)', advance='no' ) string
    close( UNIT=lun, iostat=status )
  end subroutine write_TEXTFILE_sca

!------------ Private procedures
  subroutine null_fill_1d( array, nullChar )
    ! Fill array with null chars
    ! Args
    character(len=*), dimension(:), intent(out) :: array
    character(len=1), optional, intent(in)      :: nullChar
    ! Internal variables
    character(len=1) :: myNull
    integer :: col
    integer :: pos
    ! Executable
    myNull = achar(0)
    if ( present(nullChar) ) myNull = nullChar
    do col=1, size(array)
      do pos=1, len(array(1))
        array(col)(pos:pos) = myNull
      enddo
    enddo
  end subroutine null_fill_1d

  subroutine null_fill_2d( array, nullChar )
    ! Fill array with null chars
    ! Args
    character(len=*), dimension(:,:), intent(out) :: array
    character(len=1), optional, intent(in)      :: nullChar
    ! Internal variables
    integer :: col
    ! Executable
    do col=1, size(array,2)
      call null_fill_1d( array(:, col), nullChar )
    enddo
  end subroutine null_fill_2d

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module IO_STUFF

! $Log$
! Revision 2.21  2015/07/14 23:10:56  pwagner
! Added a routine to truncate_textFile
!
! Revision 2.20  2014/07/31 20:19:08  pwagner
! Improved comments; get_nLines returns 0 for an empty file, and -1 if cant open
!
! Revision 2.19  2014/06/20 20:25:46  pwagner
! Added get_nLines
!
! Revision 2.18  2013/08/20 00:29:36  pwagner
! Print name of text file that we choke on
!
! Revision 2.17  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.16  2013/04/13 00:17:40  pwagner
! Fixed typos in commets; removed unused variables
!
! Revision 2.15  2013/04/12 00:01:07  pwagner
! Added write_textfile like existing read_ routines
!
! Revision 2.14  2012/08/14 00:22:09  pwagner
! get_lun can take optional Bottom, Top args
!
! Revision 2.13  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.12  2010/01/25 23:51:11  pwagner
! Added routines to read stdin into string variables
!
! Revision 2.11  2009/06/30 15:21:21  pwagner
! Changed intent to prevent READ_TEXTFILE_sca from leaving undefineds in string
!
! Revision 2.10  2009/06/23 18:25:43  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.9  2008/12/02 23:08:28  pwagner
! Added a print to not_used_here
!
! Revision 2.8  2008/05/02 00:02:47  pwagner
! Less efficient but more failthful
!
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
