program DumpPFAData

  ! Dump an unformatted PFAData file in the same format as a formatted file
  ! that could have been used to create it.  See ConvertPFAData for formats.

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  real, allocatable :: Absorption(:,:), dAbsDwc(:,:), dAbsDnc(:,:), dAbsDnu(:,:)
  character(255) :: FileName
  character(31) :: Format
  integer :: I, L, Lun = -1
  character(31), allocatable :: Molecules(:)
  integer :: nTemps, nPress, nMol
  character(127) :: Signal = ''
  real :: VelLin

  call getarg ( 1, fileName )
  open ( 10, file=trim(fileName), form='unformatted' )
  call getarg ( 2, fileName )
  if ( fileName /= '' ) then
    lun = 11
    open ( 11, file=trim(fileName) )
  end if

  read ( 10 ) nTemps, nPress, nMol, velLin, l, signal(:l)

  allocate ( absorption(nTemps,nPress), dAbsDwc(nTemps,nPress), &
    &        dAbsDnc(nTemps,nPress), dAbsDnu(nTemps,nPress), molecules(nMol) )

  read ( 10 ) absorption, dAbsDwc, dAbsDnc, dAbsDnu

  molecules = ''
  do i = 1, nMol
    read ( 10 ) l, molecules(i)(:l)
  end do

  write ( format, '(a,i0,a)' ) '(1p,', nTemps, 'e12.4)'
  if ( lun > 0 ) then
    write ( lun, '(a)', advance='no' ) 'Frequency averaged cross section file for'
    do i = 1, nMol
      write ( lun, '(1x,a)', advance='no' ) trim(molecules(i))
    end do
    write ( lun, '()' )
    write ( lun, '(a,a)' ) 'Signal, ', trim(signal)
    write ( lun, '(a,i0,a)' ) 'vgrid, type=Logarithmic,coordinate=Zeta,Start=???hpa,formula = [', &
      & nPress, ':???]'
    write ( lun, '(a,i0,a)' ) 'tgrid, type=Logarithmic,coordinate=logT,Start=???K,formula = [', &
      & nTemps, ':???]'
    write ( lun, '(a,f6.2,a)' ) 'Velocity Linearization', velLin, ' km/sec'
    write ( lun, '(a)' ) 'ln(Absorption (km^-1)) data(logT, logp)'
    write ( lun, format ) absorption
    write ( lun, '(a)' ) 'Dln(Absorption(km^-1)/Dwc(hPa/MHz) data(logT, logp)'
    write ( lun, format ) dAbsDwc
    write ( lun, '(a)' ) 'Dln(Absorption(km^-1)/Dnc data(logT, logp)'
    write ( lun, format ) dAbsDnc
    write ( lun, '(a)' ) 'Dln(Absorption(km^-1)/Dnu(MHz) data(logT, logp)'
    write ( lun, format ) dAbsDnu
    close ( lun )
  else
    write ( *, '(a)', advance='no' ) 'Frequency averaged cross section file for'
    do i = 1, nMol
      write ( *, '(1x,a)', advance='no' ) trim(molecules(i))
    end do
    write ( *, '()' )
    write ( *, '(a,a)' ) 'Signal, ', trim(signal)
    write ( *, '(a,i0,a)' ) 'vgrid, type=Logarithmic,coordinate=Zeta,Start=???hpa,formula = [', &
      & nPress, ':???]'
    write ( *, '(a,i0,a)' ) 'tgrid, type=Logarithmic,coordinate=logT,Start=???K,formula = [', &
      & nTemps, ':???]'
    write ( *, '(a,f6.2,a)' ) 'Velocity Linearization', velLin, ' km/sec'
    write ( *, '(a)' ) 'ln(Absorption (km^-1)) data(logT, logp)'
    write ( *, format ) absorption
    write ( *, '(a)' ) 'Dln(Absorption(km^-1)/Dwc(hPa/MHz) data(logT, logp)'
    write ( *, format ) dAbsDwc
    write ( *, '(a)' ) 'Dln(Absorption(km^-1)/Dnc data(logT, logp)'
    write ( *, format ) dAbsDnc
    write ( *, '(a)' ) 'Dln(Absorption(km^-1)/Dnu(MHz) data(logT, logp)'
    write ( *, format ) dAbsDnu
  end if

end program DumpPFAData

! $Log$
