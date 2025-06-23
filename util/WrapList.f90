! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program WrapList

  ! Given a list, one per line, wrap the list onto Fortran lines, with
  ! a comma after each element but the last, and & at the end of each
  ! line.

  ! The command line arguments are the left and right margins for
  ! wrapping.

  ! Input is from stdin, and output is to stdout.

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  logical :: AmpBefore = .false., NoFinalComma = .false.
  integer :: First = 0, I, L = 0, R = 0, Status
  character(255) :: In, Out

  i = 0
  do
    i = i + 1
    call get_command_argument ( i, in ) ! Might be left margin
    if ( in(1:1) /= '-' ) exit  ! Yup, it is
    if ( in(1:2) == '-f' ) then
      if ( in(3:) == '' ) then
        i = i + 1
        call get_command_argument ( i, in(3:) )
      end if
      read ( in(3:), *, iostat=status ) first
      if ( status /= 0 ) call usage
    else if ( in(1:2) == '-b' ) then
      ampBefore = .true.
    else if ( in(1:2) == '-n' ) then
      noFinalComma = .true.
    else
      call usage
    end if
  end do
  call get_command_argument ( i+1, out ) ! Right margin
  if ( in /= '' ) read ( in, *, iostat=status ) l
  if ( status /= 0 ) call usage
  if ( out /= '' ) read ( out, *, iostat=status ) r
  if ( status /= 0 ) call usage
  if ( l < 1 .or. r-l < 40 ) call usage
  if ( first <= 0 ) first = l

  out = ''
  i = first
  if ( ampBefore .and. first == l ) then
    out(i:i) = '&'
    i = i + 2
  end if
  do
    read ( *, '(a)', end=9 ) in
    in = adjustl(in)
    if ( in(1:1) == '' ) cycle ! Ignore blank lines
    if ( i + len_trim(in) + 3 >= r ) then
      write ( *, '(a)' ) trim(out) // ' &'
      out = ''
      i = l
      if ( ampBefore ) then
        out(i:i) = '&'
        i = i + 2
      end if
    end if
    out(i:) = trim(in) // ','
    i = len_trim(out) + 2
  end do
9 continue
  if ( noFinalComma ) then
    l = len_trim(out)
    if ( out(l:l) == ',' ) out(l:l) = ''
  else
    out = trim(out) // ' &'
  end if
  write ( *, '(a)' ) trim(out)

contains

  subroutine Usage
    call get_command_argument ( 0, in )
    print '(a)', 'Usage: ' // trim(in) // ' [options] leftMargin rightMargin'
    print '(a)', '  where leftMargin > 0 and rightMargin - leftMargin > 39'
    print '(a)', '  Options: -f[ ]# => First left margin, else the same as the others'
    print '(a)', '           -b     => Ampersand before each line, except the first'
    print '(a)', '                     if first margin /= left margin'
    print '(a)', '           -n     => No comma and ampersand on last line'
    print '(a)', '           -"anything else", or a missing margin'
    print '(a)', '                  => This explanation.'
    print '(a)'
    print '(a)', '  Suppose you have a list of USE names that you normally want'
    print '(a)', '  to start in column 5, but the colon after ONLY in the USE'
    print '(a)', '  statement is in column 18 so you want the list on that line'
    print '(a)', '  to start in column 20. Suppose you want to wrap the end at'
    print '(a)', '  column 74, and you want ampersands on each line after the'
    print '(a)', '  first one.  You might use'
    print '(a)'
    print '(a)', '  ' // trim(in) // ' -f20 -b 5 78'
    stop
  end subroutine Usage

end program WrapList

! $Log$
! Revision 1.3  2013/07/19 00:29:54  vsnyder
! Add -f, -b and -n options
!
! Revision 1.2  2013/07/18 22:56:16  vsnyder
! Initial commit
!
