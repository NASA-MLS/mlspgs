program WordSplit

! Scan some lines of text.

! Put each word (ending with comma, blank, ')' or newline) onto a separate
! line.

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

  integer :: I, L, S
  character(len=132) LINE

o:do
    read ( *, '(a)', iostat=s ) line
    if ( s /= 0 ) exit
    do
      line = adjustl(line)
      if ( line(1:1) == '!' .or. line(1:1) == ' ' .or. line(:i) == ')' ) &
  cycle o
      i = scan(line,', )')
      if ( line(:i) /= '&' ) write ( *, '(a)') line(:i-1)
      line(:i) = ''
    end do
  end do o

end program WordSplit

! $Log$
