! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Spreadsheet

! Make a spreadsheet of the molecules and bands in a PFA database.

! Input is produced by h5ls and read from standard input.

! Output is produced on standard output.

! The first argument on the command line, if any, is printed as a heading.

! Example:
!   h5ls PFAData_DACS_v2-0-6.h5 > PFAData_DACS_v2-0-6.toc
!   ./Spreadsheet "PFA tables in PFAData_DACS_v2-0-6" < PFAData_DACS_v2-0-6.toc

!------------------------------ RCS Ident Info -------------------------------
  ! "$RCSfile$"
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = IdParm
!------------------------------------------------------------------------------

  character(80) :: Head, Line

  integer, parameter :: MaxMol=500, MaxBand=100
  character(20) :: Molecules(MaxMol), Bands(MaxBand)
  character(4*maxBand+4) :: Fmts
  character :: Sheet(MaxMol,MaxBand)

  integer :: I
  integer :: NBand, Band, LenBand ! Number, one of them, length of longest string
  integer :: NMol, Mol, LenMol    ! Number, one of them, length of longest string

  call getarg ( 1, head )

  if ( head == "--version" ) then
    print *, id
    stop
  end if

  if ( head == "--help" ) then
    call getarg ( 0, head )
    print *, trim(head) // " [heading] < input from h5ls > output"
    stop
  end if

  sheet = ''

  ! Collect the unique molecules and bands in the file, and mark SHEET
  nBand = 0; nMol = 0
  lenBand = 0; lenMol = 0
  do
    read ( *, '(a)', end=1 ) line
    i = index(line,'%')
    if ( i == 0 ) cycle
    molecules(nMol+1) = line(:i-1)
    mol = lookup(molecules(:nMol+1))
    if ( mol > nMol ) nMol = nMol + 1
    lenMol = max(lenMol,i-1)
    i = index(line,'.B')
    if ( i == 0 ) cycle
    line = line(i+1:)
    i = index(line,':')
    if ( i == 0 ) cycle
    bands(nBand+1) = line(:i-1)
    band = lookup(bands(:nBand+1))
    if ( band > nBand ) nBand = nBand + 1
    lenBand = max(lenBand,i-1)
    sheet(mol,band) = "*"
  end do
1 continue

  ! Construct a format to print the filename and print it
  write ( fmts, 2 ) lenMol+1
2 format ( "(",i0,"x,a)" )
  print trim(fmts)//")", trim(head)
  ! Construct a format to print the bands and print them
  write ( fmts, 3 ) lenMol+1, ( len_trim(bands(i))+1, i=1, nBand )
3 format ( "(",i0,"x",500(:",a",i0,) )
  print trim(fmts)//")", bands(:nBand)
  ! Construct a format to print the sheet and print it
  write ( fmts, 4 ) lenMol+3, ( len_trim(bands(i)), i=1, nBand )
4 format ( "(a",i0,500(:",a1,",i0,"x") )
  print trim(fmts)//")", (molecules(i), sheet(i,:nBand), i = 1, nMol)

contains

  integer function Lookup ( Array )
    character(*), intent(in) :: Array(:)
    ! Return where array(N) is found
    integer :: L, N
    n = size(array,1)
    do lookup = 1, n-1
      if ( array(lookup) == array(n) ) return
    end do
  end function Lookup

end program Spreadsheet

! $Log$
