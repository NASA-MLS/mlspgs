program ConvertPFA

  ! Usage: ./convertPFA infile outfile

  ! Convert PFAData array files from formatted to unformatted.

  ! Emit a PFAData command:

  ! pfaData, file=<stuff before signal from outfile>$<stuff after signal from outfile>, $
  !   temperatures=pfaTgrid, vGrid=pfaVgrid, signal=<signal from the file>

  ! The input file format is

  ! Frequency averaged cross section file for <molecule list>
  ! Signal, <signal>
  ! vgrid, type=Logarithmic,coordinate=Zeta,Start=<start>,formula = [<nPress>:<perDecade>]
  ! tgrid, type=Logarithmic,coordinate=logT,Start=<start>,formula = [<nTemps>:<perDecade>]
  ! Velocity Linearization <Velocity Linearization>
  ! ln(Absorption (km^-1)) data(logT, logp)
  ! <array of nTemps x nPress absorptions>
  ! Dln(Absorption(km^-1)/Dwc(hPa/MHz) data(logT, logp)
  ! <array of nTemps x nPress Dln(Absorption(km^-1)/Dwc(hPa/MHz)>
  ! Dln(Absorption(km^-1)/Dnc data(logT, logp)
  ! <array of nTemps x nPress Dln(Absorption(km^-1)/Dnc data(logT, logp) >
  ! Dln(Absorption(km^-1)/Dnu(MHz) data(logT, logp)
  ! <array of nTemps x nPress Dln(Absorption(km^-1)/Dnu(MHz) >

  ! The output file format is

  ! nTemps, nPress, nMol, velLin, len_trim(signal), trim(signal)

  ! absorption(:nTemps,:nPress), dAbsDwc(:nTemps,:nPress),
  ! dAbsDnc(:nTemps,:nPress), dAbsDnu(:nTemps,:nPress)

  ! { len_trim(molecules(i)), trim(molecules(i)) } one record for each i = 1 : nMol

  use Machine ! May need to get GETARG from here; otherwise, nothing is used

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  character :: Ch
  character(len=255) :: Infile, Line, Outfile, Signal*127
  character(len=31) :: Molecules(255) = ' '
  integer :: I, nTemps, nPress, nMol
  real, allocatable :: Absorption(:,:), dAbsDwc(:,:), dAbsDnc(:,:), dAbsDnu(:,:)
  real :: VelLin

  call getarg ( 1, line )
  if ( line == '' ) then
    call getarg ( 0, line )
    write ( *, '(a)' ) 'Usage: ' // trim(line) // ' nTemps nPress infile outfile'
    stop
  end if
  call getarg ( 1, inFile )
  open ( 10, file=infile, form='formatted', status='old' )
  call getarg ( 2, outFile )
  open ( 11, file=outfile, form='unformatted' )

  ! The Molecules
  read ( 10, '(a)', advance='no' ) line(:41) ! The 'Frequency averaged ... part'
  nMol = 0
  do
    do ! skip blanks
      read ( 10, '(a)', advance='no', eor=9 ) ch
      if ( ch /= '' ) exit
    end do
    nMol = nMol + 1
    i = 1
    molecules(nMol)(1:1) = ch
    do ! read the rest of the symbol
      read ( 10, '(a)', advance='no', eor=9 ) ch
      if ( ch == '' ) exit
      if ( ch == '-' ) ch = '_'
      i = i + 1
      molecules(nMol)(i:i) = ch
    end do
  end do
9 continue

  ! The Signal
  read ( 10, '(a)' ) line
  signal = adjustl(line(9:))
  i = index(outFile,trim(signal))
  if ( i == 0 ) then
    write ( *, * ) 'The output file name ', trim(outfile)
    write ( *, * ) 'does not contain the signal ', trim(signal)
    stop
  end if

  write ( *, '(a)' ) 'pfaData, file="' // outfile(:i-1) // '$' // &
    & trim(outfile(i+len_trim(signal):)) // '", $' 

  write ( *, '(a)' ) '  temperatures=pfaTgrid, vGrid=pfaVgrid, signal="' // &
    & trim(signal) // '"'

  ! The vGrid line
  read ( 10, '(a)' ) line
  i = index(line,'[')
  line = line(i+1:)
  i = index(line,':')
  read ( line(:i-1), * ) nPress

  ! The tGrid line
  read ( 10, '(a)' ) line
  i = index(line,'[')
  line = line(i+1:)
  i = index(line,':')
  read ( line(:i-1), * ) nTemps

  allocate ( absorption(nTemps,nPress), dAbsDwc(nTemps,nPress), &
    &        dAbsDnc(nTemps,nPress), dAbsDnu(nTemps,nPress) )

  ! The Velocity Linearization line
  read ( 10, '(a)' ) line
  line = adjustl(line(23:))
  read ( line, * ) velLin

  ! The "ln(Absorption (km^-1)) data(logT, logp)" line (ignore it)
  read ( 10, * ) line
  ! The "ln(Absorption (km^-1)) data(logT, logp)" data
  read ( 10, * ) absorption

  ! The "Dln(Absorption(km^-1)/Dwc(hPa/MHz) data(logT, logp)" line (ignore it)
  read ( 10, * ) line
  ! Th "Dln(Absorption(km^-1)/Dwc(hPa/MHz) data(logT, logp)" data
  read ( 10, * ) dAbsDwc

  ! The "Dln(Absorption(km^-1)/Dnc data(logT, logp)" line (ignore it)
  read ( 10, * ) line
  ! The "Dln(Absorption(km^-1)/Dnc data(logT, logp)" data
  read ( 10, * ) dAbsDnc

  ! The "Dln(Absorption(km^-1)/Dnu(MHz) data(logT, logp)" line (ignore it)
  read ( 10, * ) line
  ! The "Dln(Absorption(km^-1)/Dnu(MHz) data(logT, logp)" data
  read ( 10, * ) dAbsDnu

  ! Write the output
  write ( 11 ) nTemps, nPress, nMol, velLin, len_trim(signal), trim(signal)

  write ( 11 ) absorption, dAbsDwc, dAbsDnc, dAbsDnu

  do i = 1, nMol
    write ( 11 ) len_trim(molecules(i)), trim(molecules(i))
  end do

  ! All done
  close ( 10 )
  close ( 11 )

end program ConvertPFA

! $Log$
! Revision 1.4  2004/07/08 20:59:03  vsnyder
! Cannonball polishing
!
! Revision 1.3  2004/06/17 01:04:37  vsnyder
! Added 'use Machine' in case it's needed for GETARG
!
! Revision 1.2  2004/06/17 00:42:00  vsnyder
! Repaired a comment
!
! Revision 1.1  2004/06/17 00:39:58  vsnyder
! Initial commit
!
