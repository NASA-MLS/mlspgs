! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module IGRF_INT

! Model for magnetic field.

! This version has the coefficients built in, instead of being read from
! a file.

! Based on the package by
!               Dr. Dieter Bilitza, (301)513-1664
!               GSFC, NSSDC, Code 933, Greenbelt, MD 20771, USA.
!               BILITZA%NSSDCA.SPAN@DFTNIC.BITNET

  use UNITS, only: Rad2Deg, Deg2Rad
  use GEOMETRY, only: EarthRadA, EarthRadB

  implicit NONE
  private
  public :: FELDCOF, FELDC, FELDG, FINDB0, IGRF_SUB, SHELLC, SHELLG, TO_CART

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    &  "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    &  "$RCSfile$"
  private :: not_used_here
!---------------------------------------------------------------------------

! *****     Private Variables     **************************************

! AQUAD   Square of major half axis for earth ellipsoid
! BQUAD   Square of minor half axis for earth ellipsoid

  real, parameter :: AQUAD = (EarthRadA/1000.0)**2 ! Km**2
  real, parameter :: BQUAD = (EarthRadB/1000.0)**2 ! Km**2

! ERAD    Earth radius for normalization of cartesian
!         coordinates (6371.2 Km) as recommended by the International
!         Astronomical Union

  real, save ::              ERAD = 6371.2

! G(M)    Normalized field coefficients (see FELDCOF)
!         Size = NMAX*(NMAX+2)+1

  real, save, allocatable :: G(:)

! H       ??? (See FELDI)

  real, save, allocatable :: H(:)

! NMAX    Maximum order of spherical harmonics

  integer, save ::           NMAX

! SP      Dipole oriented coordinates from SHELLG; P(1,*),
!         P(2,*), P(3,*) closest to magnetic equator

  real, save ::              SP(3)

! U       Used to convert between geographic and geomagnetic coordinates

  real, parameter :: U(3,3) = reshape ( (/ +0.3511737,-0.9148385,-0.1993679, &
    &                                      +0.9335804,+0.3583680,+0.0000000, &
    &                                      +0.0714471,-0.1861260,+0.9799247 /), &
    &                                   (/ 3,3 /) )

! *****     Public Procedures     **************************************

contains

  ! ---------------------------------------------------  IGRF_SUB  -----
  subroutine IGRF_SUB ( WHERE, YEAR, XL, ICODE, DIP, DEC )

  !----------------------------------------------------------------
  !   Input:
  !	WHERE   (1) Geodetic latitude in degrees,
  !	        (2) Geodetic longitude in degrees,
  !	        (3) Height in km
  !	YEAR    Decimal year (year+month/12.0-0.5 or year+day-of-year/365
  !		Or 366 if leap year)
  !   Output:
  !     XL      L value                                  
  !     ICODE   =1  L is correct; =2  L is not correct;  
  !             =3  an approximation is used             
  !     DIP     Geomagnetic inclination in degrees       
  !     DEC     Geomagnetic declination in degress       
  !----------------------------------------------------------------

    real, intent(in) :: WHERE(3), YEAR
    real, intent(out) :: XL, DIP, DEC
    integer, intent(out) :: ICODE

    real ::           BABS, BDOWN, BEAST, BNORTH, DIMO

!----------------CALCULATE PROFILES-----------------------------------

    call feldcof ( year, dimo )
    call feldg ( where, bnorth, beast, bdown, babs )
    call shellg ( where, dimo, xl, icode )
    dip = asin(bdown/babs) * rad2Deg
    dec = asin(beast/sqrt(beast*beast+bnorth*bnorth)) * rad2Deg

  end subroutine IGRF_SUB

  ! ------------------------------------------------------  FINDB0  -----
  subroutine FINDB0 ( STPS, BDEL, VALUE, BEQU, RR0 )
!--------------------------------------------------------------------
! Find smallest magnetic field strength on field line

! Input:   STPS   Step size for field line tracing
!          BDEL   Required accuracy  = [ B(LAST) - BEQU ] / BEQU
!                 B(LAST)  Is field strength before BEQU

!          Module variable
!          SP     Dipole oriented coordinates from SHELLG; P(1,*),
!                 P(2,*), P(3,*) closest to magnetic equator

! Output:  VALUE  =.FALSE. if BEQU is not minimal value on field line
!          BEQU   Magnetic field strength at magnetic equator
!          RR0    Equatorial radius normalized to earth radius
!          BDEL   Final achieved accuracy
!--------------------------------------------------------------------

    real, intent(in) :: STPS
    real, intent(inout) :: BDEL
    logical, intent(out) :: VALUE
    real, intent(out) :: BEQU
    real, intent(out), optional :: RR0

    real ::           B, BMIN, BDELTA, BOLD, BQ1, BQ2, BQ3, P(8,4)
    real ::           R2, R3, ROLD, STEP, STEP12
    integer ::        IRUN

    step = stps
    irun = 0
    do
      irun = irun+1
      if ( irun > 5 ) then
        value = .false.
        exit
      end if
!*********************FIRST THREE POINTS
      p(1:3,2) = sp(1:3)
      call stoer_start ( p, step, bq1, bq2, bq3, r2, r3 )
!******************INITIALIZATION
      step12 = step/12.
      value = .true.
      bmin = huge(1.0)
      bold = huge(1.0)
!******************CORRECTOR (FIELD LINE TRACING)
      do
        p(1,3) = p(1,2) + step12*(5.*p(4,3)+8.*p(4,2)-p(4,1))
        p(2,3) = p(2,2) + step12*(5.*p(5,3)+8.*p(5,2)-p(5,1))
!******************PREDICTOR (FIELD LINE TRACING)
        p(1,4) = p(1,3) + step12*(23.*p(4,3) - 16.*p(4,2)+5.*p(4,1))
        p(2,4) = p(2,3) + step12*(23.*p(5,3) - 16.*p(5,2)+5.*p(5,1))
        p(3,4) = p(3,3) + step
        call stoer ( p(:,4), bq3, r3 )
        p(1:8,1:3) = p(1:8,2:4)
        b = sqrt(bq3)
        if ( b < bmin) bmin = b
        if ( b > bold ) exit
        bold = b
        rold = 1./r3
        sp = p(1:3,4)
      end do
      if ( bold /= bmin ) value = .false.
      bdelta = (b-bold)/bold
      if ( bdelta <= bdel ) exit
      step = step/10.
    end do
    rr0 = rold
    bequ = bold
    bdel = bdelta
    if ( present(rr0) ) rr0 = rold
  end subroutine FINDB0

  ! -----------------------------------------------------  SHELLG  -----
  subroutine SHELLG ( WHERE, DIMO, FL, ICODE, B0 )
  !--------------------------------------------------------------------
  ! Calculates L-value for specified geodetic coordinates, altitude
  ! and geomagnetic field model.
  !--------------------------------------------------------------------
  !  Input:  WHERE  (1) Geodetic latitude in degrees (north),
  !                 (2) Geodetic longitude in degrees (east)
  !                 (3) Altitude in km above sea level
  !          DIMO   Dipole moment in Gauss (normalized to earth radius)

  !          Module variable
  !          H(144)     Field model coefficients adjusted for SHELLG
  !-----------------------------------------------------------------------
  !  Output: FL           L-value
  !          ICODE        =1 Normal completion
  !                       =2 Unphysical conjugate point (FL meaningless)
  !                       =3 Shell parameter greater than limit up to
  !                          which accurate calculation is required;
  !                          approximation is used.
  !          B0           Magnetic field strength in gauss
  !-----------------------------------------------------------------------

    real, intent(in) :: WHERE(3), DIMO
    real, intent(out) :: FL
    integer, intent(out), optional :: ICODE
    real, intent(out), optional :: B0

    real ::           COSPHI, CT
    real ::           SINPHI, ST, XI(3)

    call to_cart ( where, xi, ct, st, cosphi, sinphi )
    call shellc ( xi, dimo, fl, icode, b0 )
  end subroutine SHELLG

  ! -----------------------------------------------------  SHELLC  -----
  subroutine SHELLC ( XI, DIMO, FL, ICODE, B0 )
  !--------------------------------------------------------------------
  ! Calculate L-value for specified cartesian coordinates
  ! and geomagnetic field model.
  !-----------------------------------------------------------------------
  ! REF: G. KLUGE, European Space Operations Center, Internal note
  !      NO. 67, 1970.
  !      G. KLUGE, Computer Physics Communications 3, 31-35, 1972
  !--------------------------------------------------------------------
  !  Input:
  !          XI(3)  Cartesian coordinates in earth radii (ERAD km)     
  !                  X-axis pointing to equator at 0 longitude           
  !                  Y-axis pointing to equator at 90 long.              
  !                  Z-axis pointing to north pole                       
  !          DIMO   Dipole moment in Gauss (normalized to earth radius)

  !          Module variable
  !          H(144)     Field model coefficients adjusted for SHELLG
  !-----------------------------------------------------------------------
  !  Output: FL           L-value
  !          ICODE        =1 Normal completion
  !                       =2 Unphysical conjugate point (FL meaningless)
  !                       =3 Shell parameter greater than limit up to
  !                          which accurate calculation is required;
  !                          approximation is used.
  !          B0           Magnetic field strength in Gauss
  !-----------------------------------------------------------------------
    real, intent(in) :: XI(3), DIMO
    real, intent(out) :: FL
    integer, intent(out), optional :: ICODE
    real, intent(out), optional :: B0

    real ::           ARG1, ARG2, BEQU, BQ1, BQ2, BQ3, C0, C1, C2, C3
    real ::           D0, D1, D2, DIMOB0, E0, E1, E2, FF, FI, GG, HLI
    integer ::        IEQU, MyCode, N
    real ::           MyB0, ORADIK, OTERM, P(8,100)
    real ::           R2, R3, R3H, R, RADIK, RQ

!-- RMIN, RMAX are boundaries for identification of ICODE=2 and 3
    real, parameter :: RMIN = 0.05, RMAX = 1.01

!-- STEP is step size for field line tracing
!-- STEQ is step size for integration
    real ::           STEP = 0.20, STEP2, STEP12, STEQ = 0.03, STP
    real ::           T, TERM, XX, Z, ZQ

    bequ = huge(1.0e0)
    rq = 1./(xi(1)*xi(1)+xi(2)*xi(2)+xi(3)*xi(3))
    r3h = sqrt(rq*sqrt(rq))
    p(1,2) = (xi(1)*u(1,1)+xi(2)*u(2,1)+xi(3)*u(3,1))*r3h
    p(2,2) = (xi(1)*u(1,2)+xi(2)*u(2,2)             )*r3h
    p(3,2) = (xi(1)*u(1,3)+xi(2)*u(2,3)+xi(3)*u(3,3))*rq
!*****FIRST THREE POINTS OF FIELD LINE
    call stoer_start ( p, step, bq1, bq2, bq3, r2, r3 )
    myb0 = sqrt(bq2)
!*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
    if ( bq1 < bequ ) then
      bequ = bq1
      iequ = 1
    end if
    if ( bq2 < bequ ) then
      bequ = bq2
      iequ = 2
    end if
    if ( bq3 < bequ ) then
      bequ = bq3
      iequ = 3
    end if
!*****INITIALIZATION OF INTEGRATION LOOPS
    step12 = step/12.
    step2 = step+step
    steq = sign(steq,step)
    fi = 0.
    myCode = 1
    oradik = 0.
    oterm = 0.
    stp = r2*steq
    z = p(3,2)+stp
    stp = stp/0.75
    p(8,1) = step2*(p(1,1)*p(4,1)+p(2,1)*p(5,1))
    p(8,2) = step2*(p(1,2)*p(4,2)+p(2,2)*p(5,2))
!*****MAIN LOOP (FIELD LINE TRACING)
o:  do n = 3, size(p,2)-1
!*****CORRECTOR (FIELD LINE TRACING)
      p(1,n) = p(1,n-1) + step12*(5.*p(4,n)+8.*p(4,n-1)-p(4,n-2))
      p(2,n) = p(2,n-1) + step12*(5.*p(5,n)+8.*p(5,n-1)-p(5,n-2))
!*****PREPARE EXPANSION COEFFICIENTS FOR INTERPOLATION
!*****OF SLOWLY VARYING QUANTITIES
      p(8,n) = step2*(p(1,n)*p(4,n)+p(2,n)*p(5,n))
      c0 = p(1,n-1)**2 + p(2,n-1)**2
      c1 = p(8,n-1)
      c2 = (p(8,n) - p(8,n-2))*0.25
      c3 = (p(8,n) + p(8,n-2)-c1-c1)/6.0
      d0 = p(6,n-1)
      d1 = (p(6,n) - p(6,n-2))*0.5
      d2 = (p(6,n) + p(6,n-2)-d0-d0)*0.5
      e0 = p(7,n-1)
      e1 = (p(7,n) - p(7,n-2))*0.5
      e2 = (p(7,n) + p(7,n-2)-e0-e0)*0.5
!*****INNER LOOP (FOR QUADRATURE)
      do
        t = (z-p(3,n-1))/step
        if ( t > 1.) exit
        hli = 0.5*(((c3*t+c2)*t+c1)*t+c0)
        zq = z*z
        r = hli + sqrt(hli*hli+zq)
        if ( r <= rmin) then
          !*****APPROXIMATION FOR HIGH VALUES OF L.
          myCode = 3
          t = -p(3,n-1)/step
          fl = 1./(abs(((c3*t+c2)*t+c1)*t+c0)+1e-15)
          return
        end if
        rq = r*r
        ff = sqrt(1.+3.*zq/rq)
        radik = myb0-((d2*t+d1)*t+d0)*r*rq*ff
        if ( r > rmax) then
          myCode = 2
          radik = radik-12.*(r-rmax)**2
        end if
        if ( radik+radik <= oradik) exit o
        term = sqrt(radik)*ff*((e2*t+e1)*t+e0)/(rq+zq)
        fi = fi + stp*(oterm+term)
        oradik = radik
        oterm = term
        stp = r*steq
        z = z+stp
      end do
!*****PREDICTOR (FIELD LINE TRACING)
      p(1,n+1) = p(1,n) + step12*(23.*p(4,n)-16.*p(4,n-1)+5.*p(4,n-2))
      p(2,n+1) = p(2,n) + step12*(23.*p(5,n)-16.*p(5,n-1)+5.*p(5,n-2))
      p(3,n+1) = p(3,n) + step
      call stoer ( p(:,n+1), bq3, r3 )
!*****SEARCH FOR LOWEST MAGNETIC FIELD STRENGTH
      if ( bq3 < bequ ) then
        iequ = n + 1
        bequ = bq3
        end if
    end do o
    if ( iequ < 2) iequ = 2
    sp = p(1:3,iequ-1)
    if ( oradik >= 1e-15) fi = fi+stp/0.75*oterm*oradik/(oradik-radik)

!-- The minimal allowable value of FI was changed from 1E-15 to 1E-12,
!-- because 1E-38 is the minimal allowable arg. for LOG in our envir.
!-- D. Bilitza, Nov 87.

    fi = 0.5*abs(fi)/sqrt(myb0) + 1.0e-12

!*****Compute L from B and I.  Same as CARMEL in INVAR.

!-- Correct dipole moment is used here. D. Bilitza, Nov 87.

    dimob0 = dimo/myb0
    arg1 = log(fi)
    arg2 = log(dimob0)
    xx = 3*arg1-arg2
    if ( xx > 23.0) then
      gg = xx - 3.0460681e0
    else if ( xx > 11.7) then
      gg = (((((2.8212095e-8*xx-3.8049276e-6)*xx+2.170224e-4)*xx- &
        & 6.7310339e-3 )*xx+1.2038224e-1)*xx-1.8461796e-1)*xx+2.0007187e0
    else if ( xx > +3.0) then
      gg = ((((((((6.3271665e-10*xx-3.958306e-8)*xx+9.9766148e-07)*xx- &
        & 1.2531932e-5)*xx+7.9451313e-5)*xx-3.2077032e-4)*xx+2.1680398e-3)* &
        & xx+1.2817956e-2)*xx+4.3510529e-1)*xx+6.222355e-1
    else if ( xx > -3.0) then
      gg = ((((((((2.6047023e-10*xx+2.3028767e-9)*xx-2.1997983e-8)*xx- &
        & 5.3977642e-7)*xx-3.3408822e-6)*xx+3.8379917e-5)*xx+1.1784234e-3)* &
        & xx+1.4492441e-2)*xx+4.3352788e-1)*xx+6.228644e-1
    else if ( xx > -22.) then
      gg = ((((((((-8.1537735e-14*xx+8.3232531e-13)*xx+1.0066362e-9)*xx+ &
        & 8.1048663e-8)*xx+3.2916354e-6)*xx+8.2711096e-5)*xx+1.3714667e-3)* &
        & xx+1.5017245e-2)*xx+4.3432642e-1)*xx+6.2337691e-1
    else
      gg = 3.33338e-1*xx+3.0062102e-1
    end if
    fl = exp(log((1.0+exp(gg))*dimob0)/3.0)
    if ( present(icode) ) icode = mycode
    if ( present(b0) ) b0 = myB0

  end subroutine SHELLC

  ! ------------------------------------------------------  FELDG  -----
  subroutine FELDG ( WHERE, BNORTH, BEAST, BDOWN, BABS )
  !-------------------------------------------------------------------
  ! Calculate earth magnetic field from spherical harmonics model
  ! REF: G. KLUGE, European Space Operations Centre, Internal note 61,
  !      1970.
  !--------------------------------------------------------------------
  !  Input:       WHERE  (1) Geodetic latitude in degrees (north),
  !                      (2) Geodetic longitude in degrees (east),
  !                      (3) Altitude in km above sea level

  !------------------------------------------------------------------------
  !  Output: BABS   Magnetic field strength in Gauss
  !          BNORTH, BEAST, BDOWN   Components of the field with respect
  !                 to the local geodetic coordinate system, with axis
  !                 pointing in the tangential plane to the north, east
  !                 and downward.
  !-----------------------------------------------------------------------
    real, intent(in) :: WHERE(3)
    real, intent(out) :: BNORTH, BEAST, BDOWN, BABS

    real ::           B(3), BRHO, CP, CT, SP, ST, V(3)

    call to_cart ( where, v, ct, st, cp, sp )
    call feldc ( v, b )
    babs = sqrt(b(1)*b(1)+b(2)*b(2)+b(3)*b(3))   
    beast = b(2)*cp-b(1)*sp                      
    brho = b(2)*sp+b(1)*cp                       
    bnorth = b(3)*st-brho*ct                     
    bdown = -b(3)*ct-brho*st                     
  end subroutine FELDG

  ! ------------------------------------------------------  FELDC  -----
  subroutine FELDC ( V, B )

  !  Input:       V(3)  Cartesian coordinates in earth radii, km (ERAD)
  !                       X-axis pointing to equator at 0 longitude
  !                       Y-axis pointing to equator at 90 long.
  !                       Z-axis pointing to north pole
  !  Output:      B(3)  Components of the magnetic field with respect to
  !                     cartesian coordinates

    real, intent(in) :: V(3)
    real, intent(out) :: B(3)
    real :: RQ, S, T, XI(3)

    rq = 1./(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
    xi = v*rq
    call feldi ( xi ) ! Fills H
    s = .5*h(1)+2.*(h(2)*xi(3)+h(3)*xi(1)+h(4)*xi(2))
    t = (rq+rq)*sqrt(rq)
    b(1) = t*(h(3)-s*v(1))
    b(2) = t*(h(4)-s*v(2))
    b(3) = t*(h(2)-s*v(3))
  end subroutine FELDC

  ! ----------------------------------------------------  FELDCOF  -----
  subroutine FELDCOF ( YEAR, DIMO )
  !------------------------------------------------------------------------
  !  Determine coefficients and dipole moment from IGRF models

  !  Input:       YEAR  Decimal year for which geomagnetic field is to  
  !                     be calculated                                   
  !  Output:      DIMO  Geomagnetic dipole moment in GAUSS (normalized  
  !                     to earth's radius) at the time (year)           
  !-----------------------------------------------------------------------
    real, intent(in) :: YEAR
    real, intent(out), optional :: DIMO

    ! GHS     From coefficients file.  Created by Remake_GH program using
    !         the original GETSHC subroutine.

    include 'gh.f9h'

    real, allocatable ::         GHA(:)
    integer ::      I
!   IS is .true. for Schmidt normalization, .false. for Gauss normalization
    logical, parameter :: IS = .true.
    integer ::      IYEA, J, L, M, N, NUMYE
    real ::         SQRT2
    double precision :: F, F0, X

    numye = size(ghs) - 1
    iyea  =  int(year/5.)*5
    l  =  max(1.0,min(real(numye),(iyea - ghs(1)%year)/5.0 + 1.0))
    m = max(ghs(l)%ngh, ghs(l+1)%ngh)
    if ( allocated(g) ) then
      if ( size(g) < m+1 ) deallocate ( g )
    end if
    if ( allocated(h) ) then
      if ( size(h) < m+1 ) deallocate ( h )
    end if
    allocate ( gha(m) )
    if ( .not. allocated(g) ) allocate ( g(m+1) )
    if ( .not. allocated(h) ) allocate ( h(m+1) )
!-- Determine IGRF coefficients for year
    if ( l < numye ) then
      call intershc ( year, ghs(l), ghs(l+1), nmax, gha )
    else
      call extrashc ( year, ghs(l), ghs(l+1), nmax, gha )
    end if
!-- Determine magnetic dipole moment and coeffiecients G
    f0 = 0.0d0
    do j = 1, 3
      f = gha(j) * 1.0d-5
      f0 = f0 + f * f
    end do
    if ( present(dimo) ) dimo = sqrt(f0)

    g(1) = 0.0
    i = 2
    f0 = 1.0d-5
    if ( is ) f0 = -f0
    sqrt2 = sqrt(2.)

    do n = 1, nmax
      x = n
      f0 = f0 * x * x / (4.0d0 * x - 2.d0)
      if ( is ) f0 = f0 * (2.0d0 * x - 1.0d0) / x
      f = f0 * 0.5d0
      if ( is ) f = f * sqrt2
      g(i) = gha(i-1) * f0
      i = i + 1
      do m = 1 , n
        f = f * (x + m) / (x - m + 1.d0)
        if ( is ) f = f * sqrt((x - m + 1.0d0) / (x + m))
        g(i) = gha(i-1) * f
        g(i+1) = gha(i) * f
        i = i + 2
      end do
    end do
    deallocate ( gha )

  contains

    ! .................................................  INTERSHC  .....
    subroutine INTERSHC ( DATE, GH1, GH2, NMAX, GH )

    ! ===============================================================
    !
    !       Version 1.01
    !
    !       Interpolates linearly, in time, between two spherical
    !       harmonic models.
    !
    !       Input:
    !           DATE  - Date of resulting model (in decimal year)
    !           GH1   - G and H structure, type GH_T, including
    !              YEAR - year epoch for parameters
    !              NMAX - Maximum degree and order
    !              GH   - Schmidt quasi-normal internal spherical
    !                   harmonic coefficients of earlier model
    !           GH2   - Same as GH1, but for a later year
    !
    !       Output:
    !           GH    - Coefficients of resulting model
    !           NMAX  - Maximum degree and order of resulting model
    !
    !       A. Zunde
    !       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
    !
    ! ===============================================================

      real, intent(in) :: DATE
      type(gh_t), intent(in) :: GH1, GH2
      real, intent(out) :: GH(:)
      integer, intent(out) :: NMAX

    ! ---------------------------------------------------------------
    !       The coefficients (GH) of the resulting model, at date
    !       DATE, are computed by linearly interpolating between the
    !       coefficients of the earlier model (GH1), at date DTE1,
    !       and those of the later model (GH2), at date DTE2. If one
    !       model is smaller than the other, the interpolation is
    !       performed with the missing coefficients assumed to be 0.
    ! ---------------------------------------------------------------

      real :: FACTOR
      integer :: K, L

      factor = (date - gh1%year) / (gh2%year - gh1%year)

      if ( gh1%nmax == gh2%nmax ) then
        k = gh1%nmax * (gh1%nmax + 2)
        nmax = gh1%nmax
      else if ( gh1%nmax > gh2%nmax ) then
        k = gh2%nmax * (gh2%nmax + 2)
        l = gh1%nmax * (gh1%nmax + 2)
        gh(k+1:l) = gh1%gh(k+1:l) + factor * (-gh1%gh(k+1:l))
        nmax = gh1%nmax
      else
        k = gh1%nmax * (gh1%nmax + 2)
        l = gh2%nmax * (gh2%nmax + 2)
        gh(k+1:l) = factor * gh2%gh(k+1:l)
        nmax = gh2%nmax
      end if

      gh(1:k) = gh1%gh(1:k) + factor * (gh2%gh(1:k) - gh1%gh(1:k))

    end subroutine INTERSHC

    ! .................................................  EXTRASHC  .....
    subroutine EXTRASHC ( DATE, GH1, GH2, NMAX, GH )

    ! ===============================================================
    !
    !       Version 1.01
    !
    !       Extrapolates linearly a spherical harmonic model with a
    !       rate-of-change model.
    !
    !       Input:
    !           DATE  - Date of resulting model (in decimal year)
    !           GH1   - G and H structure, type GH_T, including
    !              YEAR - year epoch for parameters
    !              NMAX - Maximum degree and order
    !              GH   - Schmidt quasi-normal internal spherical
    !                   harmonic coefficients of earlier model
    !           GH2   - Same as GH1, but for a later year
    !
    !       Output:
    !           GH    - Coefficients of resulting model
    !           NMAX  - Maximum degree and order of resulting model
    !
    !       A. Zunde
    !       USGS, MS 964, Box 25046 Federal Center, Denver, CO  80225
    !
    ! ===============================================================

      real, intent(in) :: DATE
      type(gh_t), intent(in) :: GH1, GH2
      real, intent(out) :: GH(:)
      integer, intent(out) :: NMAX

    ! ---------------------------------------------------------------
    !       The coefficients (GH) of the resulting model, at date
    !       DATE, are computed by linearly extrapolating the coef-
    !       ficients of the base model (GH1), at date DTE1, using
    !       those of the rate-of-change model (GH2), at date DTE2. If
    !       one model is smaller than the other, the extrapolation is
    !       performed with the missing coefficients assumed to be 0.
    ! ---------------------------------------------------------------

      real :: FACTOR
      integer :: K, L

      factor = (date - gh1%year)

      if ( gh1%nmax == gh2%nmax ) then
        k = gh1%nmax * (gh1%nmax + 2)
        nmax = gh1%nmax
      else if ( gh1%nmax > gh2%nmax ) then
        k = gh2%nmax * (gh2%nmax + 2)
        l = gh1%nmax * (gh1%nmax + 2)
        gh(k+1:l) = gh1%gh(k+1:l)
        nmax = gh1%nmax
      else
        k = gh1%nmax * (gh1%nmax + 2)
        l = gh2%nmax * (gh2%nmax + 2)
        gh(k+1:l) = factor * gh2%gh(k+1:l)
        nmax = gh2%nmax
      end if

      gh(1:k) = gh1%gh(1:k) + factor * gh2%gh(1:k)

    end subroutine EXTRASHC

  end subroutine FELDCOF

  ! ----------------------------------------------------  To_Cart  -----
  subroutine To_Cart ( WHERE, CART, CT, ST, CP, SP )
  ! Convert Geodetic latitude and longitude, and altitude, to Cartesian
    real, intent(in) :: WHERE(3)        ! Latitude, Longitude, Altitude
    real, intent(out) :: CART(3)        ! X, Y, Z
    real, intent(out), optional :: CT, ST, CP, SP ! Cosine, Sine of Lat, Lon
    real :: D, MyCt, MySt, MyCp, MySp, RHO, RLAT, RLON

    rlat = where(1)*deg2Rad
    myct = sin(rlat)
    myst = cos(rlat)
    d = sqrt(aquad-(aquad-bquad)*myct*myct)
    rlon = where(2)*deg2Rad
    mycp = cos(rlon)
    mysp = sin(rlon)
    cart(3) = (where(3)+bquad/d)*myct/erad
    rho = (where(3)+aquad/d)*myst/erad
    cart(1) = rho*mycp
    cart(2) = rho*mysp
    if ( present(ct) ) ct = myct
    if ( present(st) ) st = myst
    if ( present(cp) ) cp = mycp
    if ( present(sp) ) sp = mysp

  end subroutine To_Cart

! *****     Private Procedures     *************************************

  ! ------------------------------------------------------  FELDI  -----
  subroutine FELDI ( XI )
  ! Get H from G -- used for L computation

  !        Module variables
  !             H       ???
  !             NMAX    Maximum order of spherical harmonicS
  !             G(M)    Normalized field coefficients (see FELDCOF)
  !                     M=NMAX*(NMAX+2)

    real, intent(in) :: XI(3) ! Cartesian coordinates

    real ::           F
    integer ::        I, IH, IHMAX, IL, IMAX, K, LAST, M
    real ::           X, Y, Z

    ihmax = nmax*nmax+1
    last = ihmax+nmax+nmax
    imax = nmax+nmax-1
    h(ihmax:last) = g(ihmax:last)
    do k = 1, 3, 2
      i = imax
      ih = ihmax
      do
        il = ih-i
        f = 2./float(i-k+2)
        x = xi(1)*f
        y = xi(2)*f
        z = xi(3)*(f+f)
        i = i - 2
        if ( i >= 1 ) then
          do m = 3, i, 2
            h(il+m+1) = g(il+m+1)+z*h(ih+m+1)+x*(h(ih+m+3)-h(ih+m-1)) &
              &                              -y*(h(ih+m+2)+h(ih+m-2))
            h(il+m) = g(il+m)+z*h(ih+m)+x*(h(ih+m+2)-h(ih+m-2)) &
              &                        +y*(h(ih+m+3)+h(ih+m-1))
          end do
          h(il+2) = g(il+2)+z*h(ih+2)+x*h(ih+4)-y*(h(ih+3)+h(ih))
          h(il+1) = g(il+1)+z*h(ih+1)+y*h(ih+4)+x*(h(ih+3)-h(ih))
        end if
        h(il) = g(il)+z*h(ih)+2.*(x*h(ih+1)+y*h(ih+2))
        ih = il
        if ( i < k ) exit
      end do
    end do
  end subroutine FELDI

  ! ------------------------------------------------------  STOER  -----
  subroutine STOER ( P, BQ, R )
!*******************************************************************
!* Field line tracing for FINDB0 and SHELLG                        *
!*******************************************************************
    real, intent(inout) :: P(:) ! p(1:3) are in, p(4:7) are out
    real, intent(out) :: BQ, R

    ! *XM,*YM,*ZM  ARE GEOMAGNETIC CARTESIAN INVERSE CO-ORDINATES
    real ::           DR, DSQ, DX, DXM, DY, DYM, DZ, DZM, FLI, Q, RQ, WR
    real ::           XI(3), XM, YM, ZM

    zm = p(3)
    fli = p(1)*p(1)+p(2)*p(2)+1e-15
    r = 0.5*(fli+sqrt(fli*fli+(zm+zm)**2))
    rq = r*r
    wr = sqrt(r)
    xm = p(1)*wr
    ym = p(2)*wr
!*****TRANSFORM TO GEOGRAPHIC CO-ORDINATE SYSTEM
    xi(1) = xm*u(1,1) + ym*u(1,2) + zm*u(1,3)
    xi(2) = xm*u(2,1) + ym*u(2,2) + zm*u(2,3)
    xi(3) = xm*u(3,1)             + zm*u(3,3)
!*****COMPUTE DERIVATIVES
    call feldi ( xi ) ! Fills H
    q = h(1)/rq
    dx = h(3) + h(3) + q*xi(1)
    dy = h(4) + h(4) + q*xi(2)
    dz = h(2) + h(2) + q*xi(3)
!*****TRANSFORM BACK TO GEOMAGNETIC CO-ORDINATE SYSTEM
    dxm = u(1,1)*dx + u(2,1)*dy + u(3,1)*dz
    dym = u(1,2)*dx + u(2,2)*dy
    dzm = u(1,3)*dx + u(2,3)*dy + u(3,3)*dz
    dr = (xm*dxm+ym*dym+zm*dzm)/r
!*****FORM SLOWLY VARYING EXPRESSIONS
    p(4) = (wr*dxm-0.5*p(1)*dr)/(r*dzm)
    p(5) = (wr*dym-0.5*p(2)*dr)/(r*dzm)
    dsq = rq*(dxm*dxm+dym*dym+dzm*dzm)
    bq = dsq*rq*rq
    p(6) = sqrt(dsq/(rq+3.*zm*zm))
    p(7) = p(6)*(rq+zm*zm)/(rq*dzm)

  end subroutine STOER

  ! ------------------------------------------------  STOER_START  -----
  subroutine STOER_START ( P, STEP, BQ1, BQ2, BQ3, R2, R3 )
    real, intent(inout) :: P(:,:), STEP
    real, intent(out) :: BQ1, BQ2, BQ3, R2, R3
    integer :: I
    real :: R1, ZZ

    step = -sign(step,p(3,2))
    call stoer ( p(:,2), bq2, r2 )
    p(1,3) = p(1,2) + 0.5*step*p(4,2)
    p(2,3) = p(2,2) + 0.5*step*p(5,2)
    p(3,3) = p(3,2) + 0.5*step
    call stoer ( p(:,3), bq3, r3 )
    p(1,1) = p(1,2) - step*(2.*p(4,2)-p(4,3))
    p(2,1) = p(2,2) - step*(2.*p(5,2)-p(5,3))
    p(3,1) = p(3,2) - step
    call stoer ( p(:,1), bq1, r1 )
    p(1,3) = p(1,2) + step*(20.*p(4,3)-3.*p(4,2)+p(4,1))/18.
    p(2,3) = p(2,2) + step*(20.*p(5,3)-3.*p(5,2)+p(5,1))/18.
    p(3,3) = p(3,2) + step
    call stoer ( p(:,3), bq3, r3 )
    !*****INVERT SENSE IF REQUIRED
    if ( bq3 >  bq1) then
      step = -step
      r3 = r1
      bq3 = bq1
      do i = 1,7
        zz = p(i,1)
        p(i,1) = p(i,3)
        p(i,3) = zz
      end do
    end if
  end subroutine STOER_START

  ! ----------------------------------------------  not_used_here  -----
  logical function not_used_here()
  ! This just makes sure ID and ModuleName get used, so eager optimizers
  ! don't throw them away.  It's never referenced.
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module IGRF_INT

! $Log$
! Revision 2.1  2003/01/14 21:37:31  vsnyder
! Initial conversion from Bilitza's F77 code
!
