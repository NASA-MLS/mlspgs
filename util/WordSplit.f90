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
    if ( s /= 0 ) &
  exit o
    do
      line = adjustl(line)
      if ( line(1:1) == '!' .or. line(1:1) == ' ' .or. line(1:1) == ')' ) &
  cycle o
      i = scan(line,', )&!')
      if ( line(:i) /= '&' .and. i /= 1 ) write ( *, '(a)') line(:i-1)
      line(:i) = ''
    end do
  end do o

end program WordSplit

! $Log$
! Revision 1.3  2004/01/26 20:49:10  vsnyder
! Allow a word to end with !
!
! Revision 1.2  2004/01/26 20:40:09  vsnyder
! Allow a word to end with ampersand
!
! Revision 1.1  2004/01/26 20:35:05  vsnyder
! Initial commit
!
