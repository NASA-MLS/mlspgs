module OUTPUT_M

  implicit NONE
  private

  integer, save, public :: LINE_WIDTH = 120 ! Not used here, but a convenient
                                        ! place to store it
  integer, save, public :: PRUNIT = -1  ! "printer" unit, * if < 0.

  public :: BLANKS, OUTPUT
  interface OUTPUT
    module procedure output_char, output_char_array, output_double
    module procedure output_integer
  end interface

!---------------------------- RCS Ident Info -------------------------------
  character (len=256), private :: Id = &
       "$Id$"
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

  subroutine BLANKS ( N_BLANKS, ADVANCE )
  ! Output N_BLANKS blanks to PRUNIT.
    integer, intent(in) :: N_BLANKS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=3) :: ADV
    character(len=*), parameter :: B = &
    '                                                                    '
    integer :: I    ! Blanks to write in next WRITE statement
    integer :: N    ! Blanks remaining to write
    character(len=3) :: MY_ADV
    my_adv = 'no'
    if ( present(advance) ) then; my_adv = advance; end if
    n = n_blanks
    adv = 'no'
    do
      i = min(n,len(b))
      n = n - i
      if ( n == 0 ) adv = my_adv
      if ( prunit < 0 ) then
        write ( *, '(a)', advance=adv ) b(:i)
      else
        write ( prunit, '(a)', advance=adv ) b(:i)
      end if
      if ( n == 0 ) exit
    end do
    return
  end subroutine BLANKS

  subroutine OUTPUT_CHAR ( CHARS, ADVANCE )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS
    character(len=*), intent(in), optional :: ADVANCE
    character(len=3) :: MY_ADV
    my_adv = 'no'
    if ( present(advance) ) then; my_adv = advance; end if
    if ( prunit < 0 ) then
      write ( *, '(a)', advance=my_adv ) chars
    else
      write ( prunit, '(a)', advance=my_adv ) chars
    end if
  end subroutine OUTPUT_CHAR

  subroutine OUTPUT_CHAR_ARRAY ( CHARS, ADVANCE )
  ! Output CHARS to PRUNIT.
    character(len=*), intent(in) :: CHARS(:)
    character(len=*), intent(in), optional :: ADVANCE
    integer :: I ! loop inductor
    do i = 1, size(chars)
      if ( prunit < 0 ) then
        write ( *, '(a)', advance='no' ) chars(i:i)
      else
        write ( prunit, '(a)', advance='no' ) chars(i:i)
      end if
    end do
    if ( present(advance) ) then
      if ( prunit < 0 ) then
        write ( *, '(a)', advance=advance )
      else
        write ( prunit, '(a)', advance=advance )
      end if
    end if
  end subroutine OUTPUT_CHAR_ARRAY

  subroutine OUTPUT_DOUBLE ( VALUE, FORMAT, ADVANCE )
  ! Output "double" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    double precision, intent(in) :: VALUE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: ADVANCE
    integer :: I, J, K
    character(len=30) :: LINE
    character(len=3) :: MY_ADV
    my_adv = 'no'
    if ( present(advance) ) then; my_adv = advance; end if
    if ( present(format) ) then
      if ( prunit < 0 ) then
        write ( *, format, advance=my_adv ) value
      else
        write ( prunit, format, advance=my_adv ) value
      end if
    else
      write ( line, * ) value
      i = index(line,'.')
      j = scan(line(i:),'DdEe ') + i - 1
      if ( i /= 0 ) then
        if ( j == i ) j = len(line)
        i = i + 2
        k = j
        do while ( j > i )
          j = j - 1
          if ( line(j:j) /= '0' .and. line(j:j) /= ' ') exit
        end do
        line(j:) = line(k:)
      end if
      line = adjustl(line)
      k = len_trim(line)
      if ( prunit < 0 ) then
        write ( *, '(a)', advance=my_adv ) line(:k)
      else
        write ( prunit, '(a)', advance=my_adv ) line(:k)
      end if
    end if
  end subroutine OUTPUT_DOUBLE

  subroutine OUTPUT_INTEGER ( INT, PLACES, ADVANCE )
  ! Output INT to PRUNIT using at most PLACES (default zero) places
    integer, intent(in) :: INT
    integer, intent(in), optional :: PLACES
    character(len=*), intent(in), optional :: ADVANCE
    integer :: I
    character(len=6) :: LINE
    character(len=3) :: MY_ADV
    integer :: MY_PLACES
    my_places = 0
    if ( present(places) ) then; my_places = places; end if
    my_adv = 'no'
    if ( present(advance) ) then; my_adv = advance; end if
    write ( line, '(i6)' ) int
    i = max( 1, min(len(line)+1-my_places, index(line,' ',back=.true.)+1) )
    if ( prunit < 0 ) then
      write ( *, '(a)', advance=my_adv ) line(i:)
    else
      write ( prunit, '(a)', advance=my_adv ) line(i:)
    end if
    return
  end subroutine OUTPUT_INTEGER

  subroutine OUTPUT_SINGLE ( VALUE, FORMAT, ADVANCE )
  ! Output "SINGLE" to "prunit" using * format, trimmed of insignificant
  ! trailing zeroes, and trimmed of blanks at both ends.
    real, intent(in) :: VALUE
    character(len=*), intent(in), optional :: FORMAT
    character(len=*), intent(in), optional :: ADVANCE
    integer :: I, J, K
    character(len=30) :: LINE
    character(len=3) :: MY_ADV
    my_adv = 'no'
    if ( present(advance) ) then; my_adv = advance; end if
    if ( present(format) ) then
      if ( prunit < 0 ) then
        write ( *, format, advance=my_adv ) value
      else
        write ( prunit, format, advance=my_adv ) value
      end if
    else
      write ( line, * ) value
      i = index(line,'.')
      j = scan(line(i:),'DdEe ') + i - 1
      if ( i /= 0 ) then
        if ( j == i ) j = len(line)
        i = i + 2
        k = j
        do while ( j > i )
          j = j - 1
          if ( line(j:j) /= '0' .and. line(j:j) /= ' ') exit
        end do
        line(j:) = line(k:)
      end if
      line = adjustl(line)
      k = len_trim(line)
      if ( prunit < 0 ) then
        write ( *, '(a)', advance=my_adv ) line(:k)
      else
        write ( prunit, '(a)', advance=my_adv ) line(:k)
      end if
    end if
  end subroutine OUTPUT_SINGLE
end module OUTPUT_M

! $Log$
! Revision 1.1  2000/07/06 01:43:12  vsnyder
! Initial check-in
!
