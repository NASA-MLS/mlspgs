program Convert_Spectroscopy

! Convert the spectroscopy database from the "old" format to a format
! suitable for inclusion in the L2CF.

  use Output_M, only: Blanks, Output, PrUnit

  integer, parameter :: R8 = kind(0.0d0)

  character(len=127) :: File       ! File names (from stdin)
  integer :: HowMuch               ! Current line width
  integer :: I                     ! Loop inductor, subscript
  character(len=80) :: Line        ! of input
  integer, parameter :: Lun = 51   ! To read the "old" format database
  character(len=13) :: Name        ! Of a catalog item
  integer :: NL                    ! Number of lines
  real(r8) :: Qlog(3)
  integer :: Spectag

  ! Line parameters:
  Real(r8) :: DELTA
  Real(r8) :: EL
  Real(r8) :: GAMMA
  Real(r8) :: N
  Real(r8) :: N1
  Real(r8) :: N2
  Real(r8) :: PS
  Real(r8) :: STR
  Real(r8) :: V0
  Real(r8) :: W

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    "$Id$"
  character(len=len(idparm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

  print *, 'Enter "old" format database filename: '
  read ( *, '(a)', end=99 ) File

  open ( lun, file=File, status='old', action='read' )

  print *, 'Enter "new" format database filename: '
  read ( *, '(a)', end=99 ) File
  prunit = lun+1
  open ( prunit, file=File, action='write' )

  do
    read ( lun, '(a)' ) line
    line = adjustl(line)
    if ( line(:7) == 'END_CAT' ) exit
    if ( line(1:1) == '#' ) cycle
    read ( line, * ) name, spectag, nl, qlog
    do i = 1, len(name)
      if ( name(i:i) == '-' ) name(i:i) = '_'
      if ( name(i:i) >= 'A' .and. name(i:i) <= 'Z' ) &
        & name(i:i) = achar(iachar(name(i:i)) + iachar('a') - iachar('A'))
    end do
    do i = 1, nl
      read ( lun, * ) v0, el, str, w, ps, n, delta, n1, gamma, n2
      call blanks ( 2 )
      call output ( trim(name) )
      call output ( '_' )
      call output ( i )
      call output ( ': line, v0= ' )
      call output ( v0 )
      call output ( ' MHz, el= ' )
      call output ( el )
      call output ( ', str= ' )
      call output ( str )
      call output ( ', w= ' )
      call output ( w )
      call output ( ', $', advance='yes' )
      call blanks ( len_trim(name) + 8 )
      call output ( 'ps= ' )
      call output ( ps )
      call output ( ', n= ' )
      call output ( n )
      call output ( ', delta= ' )
      call output ( delta )
      call output ( ', n1= ' )
      call output ( n1 )
      call output ( ', gamma= ' )
      call output ( gamma )
      call output ( ', n2= ')
      call output ( n2, advance='yes' )
    end do
    call output ( 'spectra, molecule=' )
    call output ( trim(name) )
    call output ( ', Qlog=[' )
    do i = 1, 3
      call output ( qlog(i) )
      if ( i < 3 ) call output ( ', ' )
    end do
    if ( nl <= 1 ) then
      call output ( '], lines=' )
    else
      call output ( '], $', advance='yes' )
      call blanks ( len('spectra, ') )
      call output ( 'lines=' )
    end if
    if ( nl > 1 ) call output ( '[' )
    howMuch = len('spectra, lines=[')
    do i = 1, nl
      howMuch = howMuch + len_trim(name) + 5
      if ( howMuch > 80 ) then
        call output ( '$', advance='yes' )
        howMuch = len('spectra, lines=[')
        call blanks ( howMuch )
        howMuch = howMuch + len_trim(name) + 5
      end if
      call output ( trim(name) )
      call output ( '_' )
      call output ( i )
      if ( i < nl ) call output ( ', ' )
    end do
    if ( nl > 1 ) call output ( ']' )
    call output ( '', advance='yes' )
  end do
99 close ( lun )
end program Convert_Spectroscopy

! $Log$
! Revision 1.2  2001/04/04 02:06:15  vsnyder
! Repair CVS variables
!
! Revision 1.1  2001/04/04 02:05:10  vsnyder
! Initial commit
!
