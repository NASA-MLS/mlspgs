module IO

  implicit NONE
  private

  integer, public :: INUNIT = 10   ! Input unit
  integer, public :: OUNIT         ! Output unit = PRUNIT or TBUNIT
  integer, public :: PRUNIT = 11   ! Unit for printing
  integer, public :: TBUNIT = 12   ! Unit for writing tables

  public NEW_LINE, OUTPUT     ! Public procedures

contains

  subroutine NEW_LINE ( PRINT_IT, LINE_BUF, LINE_LEN, LINE_NO, LINE_PTR )

  ! Read a new line from INUNIT into LINE_BUF.  Print it on PRUNIT if
  ! PRINT_IT is non-zero.  Set LINE_LEN to LEN_TRIM( LINE_BUF ).
  ! Increment LINE_NO.  Set LINE_PTR to 1.

    integer, intent(in) :: PRINT_IT
    character(len=*), intent(out) :: LINE_BUF
    integer, intent(out) :: LINE_LEN
    integer, intent(inout) :: LINE_NO
    integer, intent(inout) :: LINE_PTR

    integer :: IOSTAT

    read ( inunit, '(a)', iostat = iostat ) line_buf
    if ( iostat < 0 ) then
      write ( prunit, * ) ' Unexpected end-of-file, &G put in input.'
      line_buf = '&G'
    end if
    line_no = line_no + 1
    line_ptr = 1
    line_len = len_trim(line_buf)
    if ( print_it /= 0 ) write ( prunit, 10 ) line_no, line_buf(:line_len)
10  format ( i6, '.   ', a )
    return

  end subroutine NEW_LINE

  subroutine OUTPUT ( LINE )

  ! Write LINE on OUNIT.  Fill LINE with blanks.

    character(len=*), intent(inout) :: LINE

    write ( ounit, '(a)' ) trim(line)
    line = ' '
    return
  end subroutine OUTPUT

end module IO
