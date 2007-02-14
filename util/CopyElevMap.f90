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

  character(127) :: InFile, OutFile
  integer :: I, J
  integer :: IMap(2160,1080)
  character(2) :: Map(2160,1080)
  real :: MapR(360,180)
  integer :: Recl

!---------------------------- RCS Info -------------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  call getarg ( 1, inFile )
  if ( inFile == '' ) then
    print *, 'Enter file name: '
    read '(a)', inFile
  end if

  call getarg ( 2, outFile )
  if ( outFile == '' ) then
    print *, 'Enter output file name: '
    read '(a)', outFile
  end if

  inquire ( iolength=recl ) map
  open ( 10, file=inFile, access='direct', form='unformatted', status='old', &
    & recl=recl )
  read ( 10, rec=1 ) map
  close ( 10 )

  imap = ichar(map(:,:)(1:1))*256 + ichar(map(:,:)(2:2))
  where ( imap > 32767 ) imap = imap - 65536

  do i = 1, 180
    do j = 1, 360
      mapr(j,i) = maxval(imap((j-1)*6+1:(j-1)*6+6,(180-i)*6+1:(180-i)*6+6))
    end do ! j
  end do ! i

! print '(10i6)', int(mapr)

  inquire ( iolength=recl ) mapr
  open ( 11, file=outFile, access='direct', form='unformatted', recl=recl )
  write ( 11, rec=1 ) mapr
  close ( 11 )

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
