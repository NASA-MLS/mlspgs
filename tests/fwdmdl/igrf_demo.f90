! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

program IGRF_DEMO

! This is an alternative "demonstration driver" program for the
! IGRF model by Dr. Dieter Bilitza, (301)513-1664
!               GSFC, NSSDC, Code 933, Greenbelt, MD 20771, USA.
!               BILITZA%NSSDCA.SPAN@DFTNIC.BITNET

! Van Snyder, 818/354-6271
! Jet Propulsion Laboratory
! 4800 Oak Grove Drive, Mail Stop 183-701
! Pasadena, CA 91109-8099 USA
! van.snyder@jpl.nasa.gov

  use IGRF_INT, only: FELDCOF, FELDG, FINDB0, SHELLG
  use UNITS, only: Rad2Deg
  implicit NONE

  character(len=4) :: Edition = '2000'  ! Edition of model

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
!---------------------------------------------------------------------------

  real :: BABS                     ! Magnitude of field
  real :: BBOT(4) = (/ - 90.0, +   0.0, 0.00000, 1940.0 /) ! begin bottom
  real :: BDEL                     ! Tolerance for FINDB0
  real :: BDOWN                    ! Downward component of field
  real :: BEAST                    ! Eastward component of field
  real :: BEQ                      ! Output from FINDB0
  real :: BEQU                     ! Either BEQ or DIMO / (XL * XL * XL)
  integer :: BNORM = 0             ! B field normalization?
                                   ! 0 = absolute (Gauss), 1 = B/B0
  character(len=7), parameter :: BNORMT(0:1) = (/ 'B/Gauss', '   B/B0' /)
  real :: BNORTH                   ! Northward component of field
  real :: DEC                      ! Geomagnetic declination in degress
  real :: DIMO                     ! Output from FELDCOF
  real :: DIP                      ! Geomagnetic inclination in degrees
  real :: ETOP(4) = (/ + 90.0, + 360.0, 30000.0, 2010.0 /) ! end top
  logical :: GOOD                  ! All the input is good?
  integer :: I                     ! Loop inductor, subscript, temp
  integer :: ICODE                 ! Status code from SHELLG
  integer :: IOSTAT                ! I/O status
  character(len=4), parameter :: ITEXT(4) = (/ 'Lat ', 'Lon ', 'H/Km', 'Year' /)
  integer :: IVAR = 3              ! Which variable to step?
                                   ! 1 = lat, 2 = lon, 3 = height, 4 = time
  logical :: NEED_COEFF            ! Need to call FELDCOF
  logical :: VAL                   ! Output from FINDB0
  real :: VARB(4) = (/ + 45.1, + 293.1,   100.0, 1985.5 /) ! begin bounds
  real :: VARE(4) = (/ + 90.0, + 360.0,  1000.0, 2010.0 /) ! end bounds
  real :: VARS(4) = (/ +  1.0, +   1.0, + 100.0, +  0.1 /) ! steps
  real :: XL                       ! Output from SHELLG
  real :: XVAR(4)                  ! Independent variables current values

  namelist /in/ ivar, varb, vare, vars, bnorm

o:do
    ! Get new input values
    good = .false.
    do
      write ( *, 10 ) ivar, bbot, varb, vare, etop, vars, bnorm
10    format ( 'Current values:'/ &
        &      'Variable to step (IVAR):', i2/ &
        &      '   IVAR values:       1         2            3         4'/  &
        &      '   IVAR meaning:   Lat(deg)  Lon(deg)  Height(Km)  Time(Yr)'/ &
        &      'Minimum values:    ', f8.1,f10.1,f12.1,f10.1/ &
        &      'Begin values(VARB):', f8.1,f10.1,f12.1,f10.1/ &
        &      'End values(VARE):  ', f8.1,f10.1,f12.1,f10.1/ &
        &      'Maximum values:    ', f8.1,f10.1,f12.1,f10.1/ &
        &      'Step values(VARS): ', f8.1,f10.1,f12.1,f10.1/ &
        &      'Non-stepped variables stay at Begin value'/   &
        &      'Output normalization(BNORM):', i2,', 0 => B (Gauss), 1 => B/B0'// &
        &      'Enter new values for IVAR, VARB(:), VARE(:), VARS(:) or BNORM in the'/ &
        &      'following format'/ &
        &      ' &in name=values, name=values... /'/  &
        &      'Input is not case sensitive.  Use as many lines as needed.  Values not'/ &
        &      'needing change need not be entered.  Whole arrays do not need subscripts.'/ &
        &      'If only one element of an array value is to be changed, ', &
        &      'use "name(index)".'/ &
        &      'If a part of an array value is to be changed, use ', &
        &      '"name(start:end)".' )

      read ( *, in, iostat = iostat )
      if ( iostat < 0 ) exit o ! end of file
      if ( iostat > 0 ) then   ! error
        write ( *, * ) '*** Some name is wrong or does not have the correct ', &
          &            'type of associated value'
        write ( *, * ) 'Press enter to continue'
        read ( *, *, iostat=iostat )
        if ( iostat < 0 ) exit o ! end of file
        cycle
      end if
      good = .true.
      if ( ivar < 1 .or. ivar > 4 ) then
        write ( *, * ) '*** IVAR is not in range 1 to 4'
        good = .false.
      end if
      do i = 1, 4
        if ( varb(i) < bbot(i) .or. varb(i) > etop(i) ) then
          write ( *, * ) '*** VARB(', i, ') is not in the range', &
            & bbot(i), ' to', etop(i)
          good = .false.
        end if
        if ( vare(i) < bbot(i) .or. vare(i) > etop(i) ) then
          write ( *, * ) '*** VARE(', i, ') is not in the range', &
            & bbot(i), ' to', etop(i)
          good = .false.
        end if
      end do
      if ( bnorm /= 0 .and. bnorm /= 1 ) then
        write ( *, * ) '*** BNORM is not 0 or 1'
        good = .false.
      end if
      if ( good ) exit
        write ( *, * ) 'Press enter to continue'
        read ( *, *, iostat=iostat )
        if ( iostat < 0 ) exit o ! end of file
        cycle
    end do

    ! Do the calculations
    write ( *, 20 ) itext(ivar), bnormt(bnorm)
20  format (// 5x, a4, '   DIMO  ', a7, ' B-NORTH  B-EAST  B-DOWN ', &
      &      '   DIP    DEC  L-VALUE C')
    xvar = varb
    need_coeff = .true.
    do
      if ( need_coeff ) call feldcof ( xvar(4), dimo )
      need_coeff = ivar == 4
      call feldg ( xvar(1:3), bnorth, beast, bdown, babs )
      call shellg ( xvar(1:3), dimo, xl, icode )
      dip = asin (bdown / babs) * Rad2Deg
      dec = asin (beast / sqrt (beast * beast + bnorth * bnorth) ) * Rad2Deg
      if ( abs (icode) > 9 ) then
        write (*, *) ' ICODE=', icode,' is set to 2'
        icode = 2
      end if
      if ( bnorm == 0 ) then
        write (*, 30) xvar(ivar), dimo, babs, &
          & bnorth, beast, bdown, dip, dec, xl, icode
30      format (f9.2,f8.4,f8.5,3(1x,f7.5),2f7.1,f8.3,i3)
      else
        bequ = dimo / (xl * xl * xl)
        if ( icode == 1 ) then
          bdel = 1.0e-3
          call findb0 ( 0.05, bdel, val, beq )
          if ( val ) bequ = beq
        end if
        write (*, 40) xvar(ivar), dimo, min(babs / bequ, 9999.999), &
          & bnorth, beast, bdown, dip, dec, xl, icode
40      format (f9.2,f8.4,f8.3,3(1x,f7.5),2f7.1,f8.3,i3)
      end if
      if ( (xvar(ivar) - vare(ivar)) * sign(1.0, vars(ivar)) >= 0.0 ) exit
      xvar(ivar) = xvar(ivar) + vars(ivar)
    end do
    write (*, 50) edition, xvar
50  format (1X,'------- International Geomagnetic Reference Field', &
     &          ' --- Edition ', a4, ' -------'/' LATI=',F7.1,'  LONGI=',F6.1, &
     &    '  I   DIMO is Dipol   I   C=1  L and B0 correct'/               &
     &    '  ALT=',F7.1,'   YEAR=',F6.1,'  I  Moment in Gauss',            &
     &    '  I    =2  wrong,  =3  approx.'/1X,74('-'))
    if (xvar(3) > 5000.0) &
      & write (*, *) &
      & ' !! REMINDER: this field model does not include external sources !!'
    if ( (xvar(4) < 1945.0) .or. (xvar(4) > 2005.0) ) &
      & write (*, *) &
      & ' !! REMINDER: Recommended time period is 1945 to 2005 !!'
    write ( *, '(/)' )
  end do o

contains

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end program IGRF_DEMO

! $Log$
