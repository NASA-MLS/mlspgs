! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program CopyElevMap

  ! Copy the GTOPO30 elevation map that was reduced to 10' resolution,
  ! ftp://www-surf.larc.nasa.gov/pub/surf/digelev_gtopo30.map, and
  ! further reduce its resolution to 1 degree using MAXVAL (not average).

  ! The input file is a 2160x1080 grid of two-byte big-endian signed
  ! integers, with the prime meridian and the NORTH pole at the origin.

  ! The output file is a grid of 360x180 default reals with the prime
  ! meridian and the SOUTH pole at the origin.  The coordinates of the
  ! center of the 1x1 degree square represented by the (i,j) element are
  ! i-0.5 degrees lon, j-90.5 degrees lat.

  ! The units are meters above mean sea level, whatever that means.

  ! The command-line arguments are the input file and the output file.
  ! If they're absent, the program prompts for them to be input from
  ! standard input.

  implicit NONE

  character(10) :: Date ! From Date_And_Time
  character :: Format = 'L' ! L = L3Ascii output, D = direct access output
  integer :: I, J
  integer :: IMap(2160,1080)
  character(127) :: InFile, OutFile, TheProgram
  character(2) :: Map(2160,1080)
  real :: MapR(360,180)
  integer :: Recl

!---------------------------- RCS Info -------------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  call getarg ( 0, theProgram )
! Get options and file names
  i = 0
  do
    i = i + 1
    call getarg ( i, inFile )
    if ( inFile(1:1) /= '-' ) exit ! No more options
    if ( inFile(2:2) == 'd' ) then
      format = 'D'
    else if ( infile(2:2) == 'l' ) then
      format = 'L'
    else
      print '(a,a,a)', 'Usage: ', trim(theProgram), ' [options] [[inFile] outFile]'
      print '(a)', '  Options: -d[irect]: direct access output file'
      print '(a)', '           -l[3ascii]: L3Ascii output file'
      print '(a)', '           -anything else: this message'
      stop
    end if
  end do
  call getarg ( i+1, outFile )

  if ( inFile == '' ) then
    print *, 'Enter file name: '
    read '(a)', inFile
  end if

  if ( outFile == '' ) then
    print *, 'Enter output file name: '
    read '(a)', outFile
  end if

  inquire ( iolength=recl ) map
  open ( 10, file=inFile, access='direct', form='unformatted', status='old', &
    & action='read', recl=recl )
  read ( 10, rec=1 ) map
  close ( 10 )

  imap = ichar(map(:,:)(1:1))*256 + ichar(map(:,:)(2:2))
  where ( imap > 32767 ) imap = imap - 65536

  do i = 1, 180
    do j = 1, 360
      mapr(j,i) = maxval(imap((j-1)*6+1:(j-1)*6+6,(180-i)*6+1:(180-i)*6+6))
    end do ! j
  end do ! i

  if ( format == 'D' ) then ! Direct access output
    inquire ( iolength=recl ) mapr
    open ( 11, file=outFile, access='direct', form='unformatted', recl=recl )
    write ( 11, rec=1 ) mapr
    close ( 11 )
  else if ( format == 'L' ) then ! L3Ascii output
    open ( 11, file=outFile, access='sequential', form='formatted' )
    call date_and_time ( date=date )
    date = date(1:4) // '-' // date(5:6) // '-' // date(7:8) ! yyyy-mm-dd
    write ( 11, * ) '; Written by ' // trim(theProgram) // ' on ' // date
    write ( 11, * ) '; Input file: ' // trim(inFile)
    write ( 11, * ) 'Field Surface Elevation Kilometers'
    write ( 11, * ) 'Latitude linear -89.5 180 1.0'
    write ( 11, * ) 'Longitude linear 0.5 360 1.0'
    write ( 11, * ) '; Data in latitude-major order'
    write ( 11, * ) 'date CCSDSA ' // date
    write ( 11, '(10f7.3)') transpose(mapr) / 1000.0
  end if

contains

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end program CopyElevMap

! $Log$
! Revision 1.1  2007/01/10 00:58:14  vsnyder
! Initial commit
!
