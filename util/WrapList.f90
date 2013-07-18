program WrapList

  ! Given a list, one per line, wrap the list onto Fortran lines, with
  ! a comma after each element but the last, and & at the end of each
  ! line.

  ! The command line arguments are the left and right margins for
  ! wrapping.

  ! Input is from stdin, and output is to stdout.

  integer :: I, L, R
  character(255) :: In, Out

  l = 0
  r = 0
  call get_command_argument ( 1, in ) ! Left margin
  call get_command_argument ( 2, out ) ! Right margin
  if ( in /= '' ) read ( in, * ) l
  if ( out /= '' ) read ( out, * ) r
  if ( l < 1 .or. r-l < 40 ) then
    call get_command_argument ( 0, in )
    print '(a)', 'Usage: ' // trim(in) // ' leftMargin rightMargin'
    print '(a)', '  where leftMargin > 0 and rightMargin - leftMargin > 39'
    stop
  end if

  out = ''
  i = l
  do
    read ( *, '(a)', end=9 ) in
    in = adjustl(in)
    if ( in(1:1) == '' ) cycle ! Ignore blank lines
    if ( i + len_trim(in) + 3 >= r ) then
      write ( *, '(a)' ) trim(out) // ' &'
      out = ''
      i = l
    end if
    out(i:) = trim(in) // ','
    i = len_trim(out) + 2
  end do
9 continue
  write ( *, '(a)' ) trim(out) // ' &'

end program WrapList
