! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program DumpPFAData

  ! Dump an unformatted PFAData file in the same format as a formatted file
  ! that could have been used to create it.  See ConvertPFAData for formats.

  use ISO_FORTRAN_ENV, only: Output_Unit
  use Machine ! May need to get GETARG from here; otherwise, nothing is used

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  real, allocatable :: Absorption(:,:), dAbsDwc(:,:), dAbsDnc(:,:), dAbsDnu(:,:)
  character(255) :: FileName
  character(31) :: Format
  integer :: I, L, Lun = Output_Unit
  character(31), allocatable :: Molecules(:)
  integer :: nTemps, nPress, nMol
  character(127) :: Signal = ''
  real :: TStart, TStep, VelLin, VStart, VStep

  call getarg ( 1, fileName )
  open ( 10, file=trim(fileName), form='unformatted' )
  call getarg ( 2, fileName )
  if ( fileName /= '' ) then
    lun = 11
    open ( 11, file=trim(fileName) )
  end if

  read ( 10 ) nTemps, nPress, nMol, velLin, l, signal(:l), &
    & vStart, vStep, tStart, tStep

  allocate ( absorption(nTemps,nPress), dAbsDwc(nTemps,nPress), &
    &        dAbsDnc(nTemps,nPress), dAbsDnu(nTemps,nPress), molecules(nMol) )

  read ( 10 ) absorption, dAbsDwc, dAbsDnc, dAbsDnu

  molecules = ''
  do i = 1, nMol
    read ( 10 ) l, molecules(i)(:l)
  end do

  write ( format, '(a,i0,a)' ) '(1p,', nTemps, 'e12.4)'
  write ( lun, '(a)', advance='no' ) 'Frequency averaged cross section file for'
  write ( lun, '(99(1x,a))' ) (trim(molecules(i)), i = 1, nMol)
  write ( lun, '()' )
  write ( lun, '(a,a)' ) 'Signal, ', trim(signal)
  write ( lun, '(3a,i0,3a)' ) 'vgrid, type=Logarithmic,coordinate=Zeta,Start=', &
      & trim(writeIt(vStart)),'hpa,formula = [', nPress, ':', trim(writeIt(vStep)), ']'
  write ( lun, '(3a,i0,3a)' ) 'tgrid, type=Logarithmic,coordinate=logT,Start=', &
      & trim(writeIt(tStart)),'K,formula = [', nTemps, ':', trim(writeIt(tStep)), ']'
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

contains

  character(15) function WriteIt ( Step )
    real, intent(in) :: Step
    if ( step == anint(step) ) then
      write ( writeIt, '(i0)' ) nint(step)
    else
      write ( writeIt, '(f0.3)' ) step
    end if
  end function WriteIt

end program DumpPFAData

! $Log$
! Revision 1.3  2004/10/05 23:03:51  vsnyder
! Dump grid start and step instead of ???
!
! Revision 1.2  2004/06/17 01:05:02  vsnyder
! Added 'use Machine' in case it's needed for GETARG
!
! Revision 1.1  2004/06/17 00:39:58  vsnyder
! Initial commit
!
