module Voigt_M

  ! Single-Line Absorption Software

  use MLSCommon, only: RP
  use Units, only: Pi, SqrtPi

  implicit NONE

  private

  ! Routines to compute Fadeeva/Voigt/Lorentz:
  public :: Real_Simple_Voigt, D_Real_Simple_Voigt,                     &
         &  Simple_Voigt, D_Simple_Voigt,                               &
         &  RLorentz, CLorentz, RVoigth2, CVoigth2, RVoigth6, CVoigth6, &
         &  RHui6, CHui6, RDrayson, CDrayson, Taylor

  real(rp), parameter :: OneOvSPi = 1.0_rp / sqrtPi  ! 1.0/Sqrt(Pi)

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
    & "$Id$"
  character ( len=len(idParm)) :: Id = idParm
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! ------------------------------------------  Real_Simple_Voigt  -----
  elemental subroutine Real_Simple_Voigt ( x, y, u )

! simple REAL(Voigt) function

! inputs

    real(rp), intent(in) :: x ! doppler width, del frequency ratio
    real(rp), intent(in) :: y ! doppler width, collision width ratio

! outputs

    real(rp), intent(out) :: u ! real part of Voigt

! internals

    integer :: I, J, MAXJ, N
    real(rp) :: DELY, DX, D1, D2, F, FD, IS, IT, Q(5), R(5), RS, RT, TR, TI
    real(rp) :: W, XA, X2, Y2, Y4

    !   Note that the dimensions of R and Q are 2 less than the a and b
    !   coefficients because we start the r and q coeffients at 2.

    ! For Hui
    real(rp), parameter :: a(7) = (/ &
       &     122.607931777104326_rp, 214.382388694706425_rp, &
       &     181.928533092181549_rp,  93.155580458138441_rp, &
       &      30.180142196210589_rp,   5.912626209773153_rp, &
       &       0.564189583562615_rp /)

    real(rp), parameter :: b(7) = (/ &
       &     122.607931773875350_rp, 352.730625110963558_rp, &
       &     457.334478783897737_rp, 348.703917719495792_rp, &
       &     170.354001821091472_rp,  53.992906912940207_rp, &
       &      10.479857114260399_rp /)

    ! For Drayson
    real(rp), parameter :: TwoOvSPi = 2.0_rp * oneOvSPi  ! 2.0/Sqrt(Pi)

    real(rp), parameter :: H = 0.2_rp
    real(rp), parameter :: HN(26) = (/(h*(i-1),i=1,26)/)
    real(rp), parameter :: RI(15) = (/(-i/2.0_rp,i=1,15)/)

    real(rp), parameter :: Dawson(26) = (/  &
       &   0.000000000000000000000_rp, 0.194751033368028049654_rp, &
       &   0.359943481934888104273_rp, 0.474763203662977930602_rp, &
       &   0.532101707056365429017_rp, 0.538079506912768419134_rp, &
       &   0.507273496407739614173_rp, 0.456507237526897257242_rp, &
       &   0.399939894323081412623_rp, 0.346772769114872245155_rp, &
       &   0.301340388923791966033_rp, 0.264510759950831957658_rp, &
       &   0.235313055663842576224_rp, 0.212165124242499004111_rp, &
       &   0.193550723859366792343_rp, 0.178271030610558287342_rp, &
       &   0.165461999878675203167_rp, 0.154524057736963452532_rp, &
       &   0.145041773054088859351_rp, 0.136721221674636496320_rp, &
       &   0.129348001236005115591_rp, 0.122760816006522922545_rp, &
       &   0.116835039953297254075_rp, 0.111472268532125307267_rp, &
       &   0.106593431283281074400_rp, 0.102134074424276835438_rp /)

    real(rp), parameter :: FDer1(26) = 1.0_rp-2.0_rp*hn*dawson
    real(rp), parameter :: FDer2(26) = (hn*FDer1+Dawson)/ri(2)
    real(rp), parameter :: FDer3(26) = (hn*FDer2+FDer1)/ri(3)
    real(rp), parameter :: FDer4(26) = (hn*FDer3+FDer2)/ri(4)
    real(rp), parameter :: FDer5(26) = (hn*FDer4+FDer3)/ri(5)
    real(rp), parameter :: FDer6(26) = (hn*FDer5+FDer4)/ri(6)

    ! For 2-pt Gauss-Hermite (we divide GW by Pi here)
    real(rp), parameter :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
    real(rp), parameter :: gw = 0.5_rp / sqrtPi

    ! For 3-pt Gauss-Hermite (we divide GW3 by Pi here)
    real(rp), parameter :: gx3 = 1.22474487139158904909864203735_rp ! sqrt(6)/2
    real(rp), parameter :: gx3_sq = gx3*gx3, gx3_2 = 2.0_rp * gx3
    real(rp), parameter :: gw3(2) = (/ 2.0_rp / 3.0_rp / sqrtPi , & ! Midpoint
                                    &  1.0_rp / 3.0_rp / sqrtPi /)  ! * 2, at GX3

    ! For 6-pt Gauss-Hermite
    real(rp), parameter :: gx6(3) = (/ 4.36077411927616508679e-1_rp, &
                                       1.33584907401369694971_rp,    &
                                       2.35060497367449222283_rp /)
    real(rp), parameter :: gw6(3) = (/ 7.24629595224392524092e-1_rp, &
                                       1.57067320322856643916e-1_rp, &
                                       4.53000990550884564086e-3_rp /) / Pi

    real(rp), parameter :: TwoThirds = 2.0_rp / 3.0_rp
    real(rp), parameter :: XL=5.2_rp, YL=0.05_rp, YH=0.6_rp, DYDX=(yh-yl)/xl

    ! For Taylor series;
    real(rp), parameter :: TwoOvSqpi = 2.0_rp / sqrtPi

    real(rp), parameter :: Q0 = TwoOvSqpi, Q1 = 2.0 * TwoOvSqpi / 3.0
    real(rp), parameter :: Q2 = 4.0 * TwoOvSqpi / 15.0, Q3 = 2.0 * q2

! This is sorted in likely occurance of each case

    xa = ABS(x)
    y2 = y*y

! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed

    if ( y + TwoThirds*xa > 100.0_rp ) then

      ! Here x is sqrt(ln2)*delnu / wd and y is sqrt(ln2)*wc / wd
      ! u = rlorentz ( x, y )
      ! This is the first term of the asymptotic expansion
      !{ $w(z) \sim \frac{i}{\sqrt{\pi}z}
      !   \left(1 + \sum_{m=1}^\infty \frac{1 \cdot 3 \dots (2m-1)}{(2z^2)^m}\right)$.

      u = OneOvSPi * y / (y2 + x*x)

    else if ( y + 0.6875_rp * xa > 11.0_rp ) then

      ! Drayson's quick 2pt hermite quadrature (essentially a lorentz)

      ! u = rvoigth2(xa,y)
!     u = gw * y * (1.0_rp/(y2 + (xa-gx)**2) + &
!       &           1.0_rp/(y2 + (xa+gx)**2))

      ! 3pt hermite quadrature
!     u = y * ( gw3(1) / (y2+x*x) + &
!       &       gw3(2) * ( 1.0_rp / ( y2 + (xa-gx3)**2 ) +
!       &                  1.0_rp / ( y2 + (xa+gx3)**2 ) ) )
!     except we eliminate one divide by combining the last two terms

      f = y2 + x * x
      d1 = f + gx3_sq
      d2 = xa * gx3_2
      u = y * ( gw3(1) / f + gw3(2) * d1 / (( d1 + d2 ) * (d1 - d2 )) )

    else if ( y > 0.6_rp .OR. y > yl + dydx*xa ) then

! Intermediate region

      ! u = rhui6(xa,y)

      !   Fill the r, q coefficients with the recursion
      !   r(0)=1.0, r(1) = y q(0) = 0.0, q(1) = -x

      r(1) = (y - xa) * (y + xa)
      q(1) = -2.0_rp*xa*y
      do i = 2, 5
        r(i) = xa * q(i-1) +  y * r(i-1)
        q(i) = y  * q(i-1) - xa * r(i-1)
      end do
      rs = a(1) + a(2)*y  + dot_product(a(3:),r)
      rt = b(1) + b(2)*y  + dot_product(b(3:),r) + xa * q(5) +  y * r(5)
      is =      - a(2)*xa + dot_product(a(3:),q)
      it =      - b(2)*xa + dot_product(b(3:),q) +  y * q(5) - xa * r(5)

      u = (rs*rt + is*it) / (rt*rt + it*it)

    else if ( xa > 5.2_rp ) then

! small y large x limit

      ! u = rvoigth6(xa,y)

      ! Real Voigt region IV 6pt GH integration

      u = y * ( gw6(1) * (1.0_rp/(y2 + (xa-gx6(1))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(1))**2)) + &
                gw6(2) * (1.0_rp/(y2 + (xa-gx6(2))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(2))**2)) + &
                gw6(3) * (1.0_rp/(y2 + (xa-gx6(3))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(3))**2)) )
    else

! Near line center where Doppler dominates Pressure broadening

      ! u = rdrayson(xa,y)

      x2 = xa * xa
      if ( x2 + y2 > 0.0036_rp ) then

        !******** Region I. Compute Dawson's function at x from Taylor series

        if ( y <= 1.0e-12_rp ) then
          u = exp(-x2)
        else

          !  Taylor series expansion about y = 0.0

          j = xa / h
          n = 1 + min(j, 25)
          dx = xa - hn(n)
          f = 0.0_rp
          if ( xa > 0.05*h) &
             f = (((((fDer6(n)*dx + fDer5(n))*dx + fDer4(n))*dx  + &
                   &  fDer3(n))*dx + fDer2(n))*dx + fDer1(n))*dx + Dawson(n)
          dely = -y
          fd = 1.0_rp - 2.0_rp * xa * f
          w = fd * dely
          j = 5.0_rp + (12.5_rp - xa) * 0.8_rp * y
          maxj = min(j, 14)
          do j = 2, maxj, 2
            f  = (xa*fd + f) / ri(j)
            fd  = (xa*f + fd) / ri(j+1)
            dely = -y2 * dely
            w = w + fd * dely
          end do

          u = exp(y2-xa*xa)*cos(2.0_rp*xa*y) + TwoOvSpi*w
        end if
      else
        ! Really to the origin, where the above seems to have a bug.
        ! Use three terms of a Taylor series

        ! (x2,y2) = -z**2
        x2 = y2 - x2
        y2 = -2.0_rp * x * y
        y4 = y2 * y2

        tr = q0 - q2 * y4 + ( q1 + q2 * x2 ) * x2
        ti = ( q1 + q3 * x2 ) * y2
        u = ( 1.0_rp + 0.5_rp * x2 ) * x2 + 1.0 - 0.5 * y4 - tr * y - ti * x
      end if
    end if

  end subroutine Real_Simple_Voigt

  ! -----------------------------------------------  Simple_Voigt  -----
  elemental subroutine Simple_Voigt ( x, y, u, v, du, dv )

! Simple Voigt function, also known as Fadeeva function.  exp(-z^2)*erfc(-I*z).

! inputs

    real(rp), intent(in) :: x ! doppler width, delta frequency ratio
    real(rp), intent(in) :: y ! doppler width, collision width ratio

! outputs

    real(rp), intent(out) :: u            ! real part of Voigt
    real(rp), optional, intent(out) :: v  ! imaginary part of Voigt
    real(rp), optional, intent(out) :: du ! real part of Voigt derivative
    real(rp), optional, intent(out) :: dv ! imaginary part of Voigt derivative

! internals

    ! For asymptotic expansion
    real(rp), parameter :: A = OneOvSpi, B = 1.5 * a

    ! For 2-pt Gauss-Hermite
!   real(rp), parameter :: GX = 0.70710678118655_rp ! 1.0/sqrt(2.0)
!   real(rp), parameter :: GW = 0.5_rp / sqrtPi

    real(rp), parameter :: XL=5.2_rp, YL=0.05_rp, YH=0.6_rp, DYDX=(yh-yl)/xl
    real(rp) :: P, Q, R, S, VV, XA, X2, Y2
!   real(rp) :: DENOM1, DENOM2, XM, XP
    complex(rp) :: UV

! This is sorted in likely occurance of each case

    xa = ABS(x)

! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed

    if ( y + 0.666666*xa > 100.0_rp ) then

!{    clorentz is one term of the asymptotic expansion
!     $w(z) \sim \frac{i}{z\sqrt{\pi}}
!       \left(1 + \sum_{m=1}^\infty \frac{1 \cdot 3 \dots (2m-1)}{(2z^2)^m}\right)$.
!     This exactly cancels in the derivative
!     $w^\prime(z) = \frac{2i}{\sqrt\pi} - 2 z w(z)$, and is therefore not
!     accurate enough for derivative computations.

      if ( present(du) .or. present(dv) ) then
        !{ Compute the derivative first using two terms of the asymptotic
        ! expansion, then compute w(z) from it, to avoid cancellation
        !%
        ! Taking the first three terms,
        ! $w(z) \approx \frac{i}{\sqrt{\pi}z}
        !   \left( 1 + \frac1{2z^2}\right)$.
        !%
        ! Applying $w^\prime(z) = \frac{2i}{\sqrt{\pi}} - 2 z w(z)$ we have
        ! $w^\prime(z) \approx -\frac{i}{\sqrt{\pi}z^2}$.
        !%
        ! Substituting $z^{-2} = x_2 + i y_2$ and using $a = \frac1{\sqrt{\pi}}$
        ! gives
        !%
        ! $w^\prime(z) \approx a y_2 - i a x_2$.
        !%
        ! This is four multiplies and three adds cheaper than in the next
        ! region inward.

        r = 1.0_rp / ( x*x + y*y )  ! 1 / |z|
        s = -y * r                  ! Im(1/z)
        r = x * r                   ! Re(1/z)

        x2 = ( r - s ) * ( r + s )  ! Re(1/z^2)
        y2 = 2.0 * r * s            ! Im(1/z^2)

        p = a * y2                  ! Re(w')
        q = -a * x2                 ! Im(w')

        !{ Writing $z^{-1} = r + i s$ and using
        ! $w(z) = \frac1z \left( \frac{i}{\sqrt\pi} - \frac12 w^\prime(z)\right)$
        ! we have
        !
        ! $w(z) \approx \frac12\left( -r p + s q \right) - a s +
        !       i \left[\frac12\left( -r q - s p \right) + a r \right]$.

        u = 0.5 * ( s * q - r * p ) - s * a ! Re(w)
        if ( present(v) ) v = sign(r * a - 0.5 * ( r * q + s * p ), x) ! Im(w)
        if ( present(du) ) du = p
        if ( present(dv) ) dv = q

      else ! no derivatives; use one term.  This isn't quite as accurate.

        ! uv = clorentz(xa,y) ! i/(\sqrt\pi z)
        r = oneOvSpi / ( x*x + y*y )
        u = y * r
        v = x * r

        ! If we need something more accurate, use two terms:
        ! $w(z) ~ \frac{i}{\sqrt\pi} \left( 1 + \frac1{2 z^2} )$.

!       r = 1.0_rp / ( x*x + y*y )  ! 1 / |z|
!       s = -y * r                  ! Im(1/z)
!       r = x * r                   ! Re(1/z)

!       x2 = 1.0_rp + 0.5_rp * ( r - s ) * ( r + s )  ! 1 + Re(1/(2 z^2))
!       y2 = r * s                  ! Im(1/(2 z^2))

!       u = oneOvSpi * ( r * x2 - s * y2 )
!       v = oneOvSpi * ( r * y2 + s * x2 )

      end if

    else if ( y + 0.6875_rp * xa > 11.0_rp ) then

      !{Three terms of the asymptotic expansion get seven digits correct.
      ! Compute the derivative first, then get $w(z)$ from it, because the
      ! first term of the asymptotic expansion cancels in the derivative:
      !%
      ! $w(z) \sim \frac{i}{\sqrt{\pi}z}
      !   \left( 1 + \sum_{m=1}^\infty \frac{1 \cdot 3 \dots (2m-1)}
      !                                      {(2 z^2)^m} \right)$.
      !%
      ! Taking the first three terms,
      ! $w(z) \approx \frac{i}{\sqrt{\pi}z}
      !   \left( 1 + \frac1{2z^2} + \frac3{4z^4} \right)$.
      !%
      ! Applying $w^\prime(z) = \frac{2i}{\sqrt{\pi}} - 2 z w(z)$ we have
      ! $w^\prime(z) \approx
      !  - \frac{i}{\sqrt{\pi}z^2} \left( 1 + \frac3{2z^2} \right)$.
      !%
      ! Substituting $z^{-2} = x_2 + i y_2$ and using $a = \frac1{\sqrt{\pi}}$
      ! and $b = \frac3{2\sqrt{\pi}}$ gives
      !
      ! $w^\prime(z) \approx b x_2 y_2 + ( a + b x_2 ) y_2 +
      !             i \left( -( a + b x_2 ) x_2 + b y_2^2 \right)$.
      !%
      ! Writing the common subexpressions $u = b y_2$ and $v = a + b x_2$
      ! we finally have $w^\prime(z) \approx u x_2 + v y_2 + i ( -v x_2 + u y_2 )$.

      r = 1.0_rp / ( x*x + y*y )  ! 1 / |z|
      s = -y * r                  ! Im(1/z)
      r = x * r                   ! Re(1/z)

      x2 = ( r - s ) * ( r + s )  ! Re(1/z^2)
      y2 = 2.0 * r * s            ! Im(1/z^2)

      u = b * y2        ! $\frac{3 y2}{2 \sqrt{\pi}}
      vv = a + b * x2   ! $\frac1{\sqrt{\pi}} + $\frac{3 x2}{2 \sqrt{\pi}}
      p = x2*u + y2*vv  ! Re(w')
      q = -x2*vv + y2*u ! Im(w')

      !{ Writing $z^{-1} = r + i s$ and using
      ! $w(z) = \frac1z \left( \frac{i}{\sqrt\pi} - \frac12 w^\prime(z)\right)$
      ! we have
      !
      ! $w(z) \approx \frac12\left( -r u + s v \right) - a s +
      !       i \left[\frac12\left( -r v - s u \right) + a r \right]$.

      u = 0.5 * ( s * q - r * p ) - s * a ! Re(w)
      if ( present(v) ) v = sign(r * a - 0.5 * ( r * q + s * p ), x) ! Im(w)
      if ( present(du) ) du = p
      if ( present(dv) ) dv = q

      ! Drayson's quick 2pt hermite integral (essentially a lorentz)

      ! uv = cvoigth2(xa,y)
!     xm = xa - gx
!     xp = xa + gx
!     y2 = y**2
!     denom1 = gw/(y2 + xm**2)
!     denom2 = gw/(y2 + xp**2)
!     u = y*(denom1+denom2)
!     if ( present(v) ) v = sign(xm*denom1 + xp*denom2, x)

    else

      if ( y > 0.6_rp .OR. y > yl + dydx*xa ) then

! Intermediate region

        uv = chui6(xa,y)

      else if ( xa > 5.2_rp ) then

! small y large x limit

        uv = cvoigth6(xa,y)

      else if ( x*x + y*y > 0.0036 ) then

! Near line center where Doppler dominates Pressure broadening

        uv = cdrayson(xa,y)

      else

! Very close to the line center, where cdrayson seems to have a bug

        uv = taylor(xa,y)

      end if

      u = real(uv,kind=rp)
      vv = sign(aimag(uv),x)
      if ( present(v) ) v = vv

      ! w' = 2*I/sqrt(pi) - 2*z*w
      if ( present(du) ) du = 2.0_rp * ( y*vv - x*u )
      if ( present(dv) ) dv = 2.0_rp * ( OneOvSpi - x*vv - y*u )
    end if

  end subroutine Simple_Voigt

  ! ---------------------------------------------------  RLorentz  -----
    real(rp) pure function RLorentz ( x, y )

  ! Real Lorentz

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

  ! Internals

!{ This is one term of the asymptotic expansion
!     $w(z) \sim \frac{i}{z\sqrt{\pi}}
!       \left(1 + \sum_{m=1}^\infty \frac{1 \cdot 3 \dots (2m-1)}{(2z^2)^m}\right)$.

    rlorentz = OneOvSPi * y / (y*y + x*x)

  end function RLorentz

  ! ---------------------------------------------------  CLorentz  -----
  complex(rp) elemental function CLorentz(x,y)

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    real(rp) :: denom

!{ This is one term of the asymptotic expansion
!     $w(z) \sim \frac{i}{z\sqrt{\pi}}
!       \left(1 + \sum_{m=1}^\infty \frac{1 \cdot 3 \dots (2m-1)}{(2z^2)^m}\right)$.

    denom = OneOvSPi / (x*x + y*y)
    clorentz = CMPLX(y*denom,x*denom,KIND=rp)

  end function CLorentz

  ! ---------------------------------------------------  RVoigth2  -----
  real(rp) elemental function RVoigth2 ( x, y )

! Real Voigt region IV 2pt GL integration

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    real(rp), parameter :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
    real(rp), parameter :: gw = 0.88622692545277_rp / Pi
    real(rp) :: y2

    y2 = y**2
    rvoigth2 = gw * y * (1.0_rp/(y2 + (x-gx)**2) + &
             &           1.0_rp/(y2 + (x+gx)**2))

  end function RVoigth2

  ! ---------------------------------------------------  CVoigth2  -----
  complex(rp) elemental function CVoigth2 ( x, y )

! Voigt region  2pt GL integration

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    real(rp), parameter :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
    real(rp), parameter :: gw = 0.88622692545277_rp / Pi
    real(rp) :: denom1,denom2,xm,xp,y2

    xm = x-gx
    xp = x+gx
    y2 = y**2
    denom1 = gw/(y2 + xm**2)
    denom2 = gw/(y2 + xp**2)
    cvoigth2 = CMPLX(y*(denom1+denom2),xm*denom1 + xp*denom2,KIND=rp)

  end function CVoigth2

  ! ---------------------------------------------------  RVoigth6  -----
  real(rp) elemental function RVoigth6 ( x, y )

! Real Voigt region IV 6pt GL integration

    use Units, only: Pi
    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    integer, parameter :: n = 3
    real(rp), parameter :: gx(n) = (/ 4.36077411927616508271e-1_rp, &
                                      1.33584907401369694957_rp,    &
                                      2.35060497367449222280_rp /)
    real(rp), parameter :: gw(n) = (/ 7.24629595224392524608e-1_rp, &
                                      1.57067320322856644036e-1_rp, &
                                      4.53000990550884564224e-3_rp /) / Pi
    real(rp) :: y2
    integer :: i

    rvoigth6 = 0.0
    y2 = y**2
    do i = 1 , n
      rvoigth6 = rvoigth6 + gw(i) * y * (1.0_rp/(y2 + (x-gx(i))**2) + &
                                     &   1.0_rp/(y2 + (x+gx(i))**2))
    end do

  end function RVoigth6

  ! ---------------------------------------------------  CVoigth6  -----
  complex(rp) elemental function CVoigth6 ( x, y )

! Voigt region IV 6pt GL integration

    use Units, only: Pi
    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    integer, parameter :: n = 3
    real(rp), parameter :: gx(n) = (/ 4.36077411927616508271e-1_rp, &
                                      1.33584907401369694957_rp,    &
                                      2.35060497367449222280_rp /)
    real(rp), parameter :: gw(n) = (/ 7.24629595224392524608e-1_rp, &
                                      1.57067320322856644036e-1_rp, &
                                      4.53000990550884564224e-3_rp /) / Pi
    integer :: i
    real(rp) :: denom1,denom2,xp,xm,y2

    cvoigth6 = cmplx(0.0_rp,0.0_rp)
    y2 = y**2
    do i = 1 , n
      xm = x-gx(i)
      xp = x+gx(i)
      denom1 = gw(i)/(y2 + xm**2)
      denom2 = gw(i)/(y2 + xp**2)
      cvoigth6 = cvoigth6+cmplx(y*(denom1+denom2), xm*denom1+xp*denom2, KIND=rp)
    end do

    end function CVoigth6

  ! ------------------------------------------------------  RHui6  -----
  real(rp) elemental function RHui6 ( x, y )

! Voigt region II Hui polynomial

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internal stuff

    real(rp), parameter :: a(7) = (/ &
       &     122.607931777104326_rp, 214.382388694706425_rp, &
       &     181.928533092181549_rp,  93.155580458138441_rp, &
       &      30.180142196210589_rp,   5.912626209773153_rp, &
       &       0.564189583562615_rp /)

    real(rp), parameter :: b(7) = (/ &
       &     122.607931773875350_rp, 352.730625110963558_rp, &
       &     457.334478783897737_rp, 348.703917719495792_rp, &
       &     170.354001821091472_rp,  53.992906912940207_rp, &
       &      10.479857114260399_rp /)

    integer :: i
    real :: rs,rt,is,it,r(5),q(5)

!   Note that this dimension is 2 less than the a and b coefficients
!   because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
!   q(0) = 0.0, q(1) = -x
!   Fill the r, q coefficients with this recursion

    r(1) = y**2 - x**2
    q(1) = -2.0_rp*x*y
    do i = 2 , 5
      r(i) = x * q(i-1) + y * r(i-1)
      q(i) = y * q(i-1) - x * r(i-1)
    end do
    rs = a(1) + a(2)*y + dot_product(a(3:),r)
    rt = b(1) + b(2)*y + dot_product(b(3:),r) + x * q(5) + y * r(5)
    is =      - a(2)*x + dot_product(a(3:),q)
    it =      - b(2)*x + dot_product(b(3:),q) + y * q(5) - x * r(5)

    rhui6 = (rs*rt + is*it) / (rt*rt + it*it)

!   FYI ihui6 = (is*rt - rs*it) / (rt*rt + it*it)

  end function rhui6

  ! ------------------------------------------------------  CHui6  -----
  complex(rp) elemental function CHui6 ( x, y )

! Voigt region II Hui polynomial

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

! Internal stuff

    real(rp), parameter :: a(7) = (/ &
       &     122.607931777104326_rp, 214.382388694706425_rp, &
       &     181.928533092181549_rp,  93.155580458138441_rp, &
       &      30.180142196210589_rp,   5.912626209773153_rp, &
       &       0.564189583562615_rp /)

    real(rp), parameter :: b(7) = (/ &
       &     122.607931773875350_rp, 352.730625110963558_rp, &
       &     457.334478783897737_rp, 348.703917719495792_rp, &
       &     170.354001821091472_rp,  53.992906912940207_rp, &
       &      10.479857114260399_rp /)

    integer :: i
    real :: rs,rt,is,it,r(5),q(5),denom

! Note that this dimension is 2 less than the a and b coefficients
! because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
! q(0) = 0.0, q(1) = -x
! Fill the r, q coefficients with this recursion

    r(1) = y**2 - x**2
    q(1) = -2.0_rp*x*y
    do i = 2, 5
      r(i) = x * q(i-1) + y * r(i-1)
      q(i) = y * q(i-1) - x * r(i-1)
    end do

    rs = a(1) + a(2)*y + dot_product(a(3:),r)
    rt = b(1) + b(2)*y + dot_product(b(3:),r) + x * q(5) + y * r(5)
    is =      - a(2)*x + dot_product(a(3:),q)
    it =      - b(2)*x + dot_product(b(3:),q) + y * q(5) - x * r(5)
    denom = 1.0_rp / (rt*rt + it*it)
    chui6 = cmplx((rs*rt+is*it)*denom, (is*rt-rs*it)*denom,kind=rp)

  end function chui6

  ! ---------------------------------------------------  RDrayson  -----
  real(rp) elemental function RDrayson ( x, y )

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

! Internal stuff

    integer :: I, J, MAXJ, N
    real(rp), parameter :: TwoOvSPi = 2.0_rp * oneOvSPi  ! 2.0/Sqrt(Pi)

    real(rp), parameter :: H = 0.2_rp
    real(rp), parameter :: HN(26) = (/(h*(i-1),i=1,26)/)
    real(rp), parameter :: RI(15) = (/(-i/2.0_rp,i=1,15)/)

    real(rp), parameter :: Dawson(26) = (/  &
       &   0.000000000000000000000_rp, 0.194751033368028049654_rp, &
       &   0.359943481934888104273_rp, 0.474763203662977930602_rp, &
       &   0.532101707056365429017_rp, 0.538079506912768419134_rp, &
       &   0.507273496407739614173_rp, 0.456507237526897257242_rp, &
       &   0.399939894323081412623_rp, 0.346772769114872245155_rp, &
       &   0.301340388923791966033_rp, 0.264510759950831957658_rp, &
       &   0.235313055663842576224_rp, 0.212165124242499004111_rp, &
       &   0.193550723859366792343_rp, 0.178271030610558287342_rp, &
       &   0.165461999878675203167_rp, 0.154524057736963452532_rp, &
       &   0.145041773054088859351_rp, 0.136721221674636496320_rp, &
       &   0.129348001236005115591_rp, 0.122760816006522922545_rp, &
       &   0.116835039953297254075_rp, 0.111472268532125307267_rp, &
       &   0.106593431283281074400_rp, 0.102134074424276835438_rp /)

    real(rp), parameter :: FDer1(26) = 1.0_rp-2.0_rp*hn*dawson
    real(rp), parameter :: FDer2(26) = (hn*FDer1+Dawson)/ri(2)
    real(rp), parameter :: FDer3(26) = (hn*FDer2+FDer1)/ri(3)
    real(rp), parameter :: FDer4(26) = (hn*FDer3+FDer2)/ri(4)
    real(rp), parameter :: FDer5(26) = (hn*FDer4+FDer3)/ri(5)
    real(rp), parameter :: FDer6(26) = (hn*FDer5+FDer4)/ri(6)

! internal stuff

    real(rp) :: Y2, DX, F, FD, W, dely

!******** Region I. Compute Dawson's function at x from Taylor series

    if ( y <= 1.0e-12_rp ) then
      rdrayson = exp(-X*X)
      return
    end if

!  Taylor series expansion about y = 0.0

    y2 = y*y
    j = x / h
    n = 1 + min(j, 25)
    dx = x - hn(n)
    f = 0.0_rp
    if ( x > 0.05*h) &
       f = (((((fDer6(n)*dx + fDer5(n))*dx + fDer4(n))*dx  + &
             &  fDer3(n))*dx + fDer2(n))*dx + fDer1(n))*dx + Dawson(n)
    dely = -y
    fd = 1.0_rp - 2.0_rp * x * f
    w = fd * dely
    j = 5.0_rp + (12.5_rp - x) * 0.8_rp * y
    maxj = min(j, 14)
    do j = 2, maxj,2
      f  = (x*fd + f) / ri(j)
      fd  = (x*f + fd) / ri(j+1)
      dely = -y2 * dely
      w = w + fd * dely
    end do

    rdrayson = exp(y2-x*x)*cos(2.0_rp*x*y) + TwoOvSpi*w

  end function RDrayson

  ! ---------------------------------------------------  CDrayson  -----
  complex(rp) elemental function CDrayson ( x, y )

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

! Internal stuff

    integer :: I, J, MAXJ, N
    real(rp), parameter :: TwoOvSPi = 2.0_rp * oneOvSPi  ! 2.0/Sqrt(Pi)

! This is a way to set the mesh points without having to use a SAVE statement

    real(rp), parameter :: H = 0.2_rp
    real(rp), parameter :: HN(26) = (/(H*(I-1),I=1,26)/)
    real(rp), parameter :: RI(15) = (/(-I/2.0_rp,I=1,15)/)

    real(rp), parameter :: Dawson(26) = (/  &
       &   0.000000000000000000000_rp, 0.194751033368028049654_rp, &
       &   0.359943481934888104273_rp, 0.474763203662977930602_rp, &
       &   0.532101707056365429017_rp, 0.538079506912768419134_rp, &
       &   0.507273496407739614173_rp, 0.456507237526897257242_rp, &
       &   0.399939894323081412623_rp, 0.346772769114872245155_rp, &
       &   0.301340388923791966033_rp, 0.264510759950831957658_rp, &
       &   0.235313055663842576224_rp, 0.212165124242499004111_rp, &
       &   0.193550723859366792343_rp, 0.178271030610558287342_rp, &
       &   0.165461999878675203167_rp, 0.154524057736963452532_rp, &
       &   0.145041773054088859351_rp, 0.136721221674636496320_rp, &
       &   0.129348001236005115591_rp, 0.122760816006522922545_rp, &
       &   0.116835039953297254075_rp, 0.111472268532125307267_rp, &
       &   0.106593431283281074400_rp, 0.102134074424276835438_rp /)

    real(rp), parameter :: FDer1(26) = 1.0_rp-2.0_rp*hn*Dawson
    real(rp), parameter :: FDer2(26) = (hn*fDer1+Dawson)/ri(2)
    real(rp), parameter :: FDer3(26) = (hn*fDer2+fDer1)/ri(3)
    real(rp), parameter :: FDer4(26) = (hn*fDer3+fDer2)/ri(4)
    real(rp), parameter :: FDer5(26) = (hn*fDer4+fDer3)/ri(5)
    real(rp), parameter :: FDer6(26) = (hn*fDer5+fDer4)/ri(6)

! internal stuff

    real(rp) :: DX, F, FD,wr,wi,dely,twoxy

!******** Region I. Compute Dawson's function at x from Taylor series

    j = x / h
    n = 1 + min(j, 25)
    dx = x - hn(n)
    f = 0.0_rp
    if ( x > 0.05*h) &
       f = (((((fder6(n)*dx + fder5(n))*dx + fder4(n))*dx  + &
             &  fder3(n))*dx + fder2(n))*dx + fder1(n))*dx + dawson(n)
    if ( y <= 1.0e-12_rp ) then
      cdrayson = cmplx(exp(-x*x),twoovspi * f,kind=rp)
      return
    end if

!  Taylor series expansion about y = 0.0

    wr = exp(y*y-x*x)
    twoxy = 2.0_rp*x*y
    cdrayson = cmplx(wr*cos(twoxy),-wr*sin(twoxy),kind=rp)
    dely = -twoovspi*y
    fd = 1.0_rp - 2.0_rp * x * f
    wi = twoovspi*f
    wr = fd*dely
    j = 5.0 + 10.0*y - 0.4*twoxy
    maxj = min(j, 14)
    do j = 2, maxj,2
      f  = (x*fd + f) / ri(j)
      fd = (x*f + fd) / ri(j+1)
      dely = y * dely
      wi = wi + f*dely
      dely = -y*dely
      wr = wr + fd*dely
    end do

    cdrayson = cdrayson + cmplx(wr,wi,kind=rp)

  end function CDrayson

  ! -----------------------------------------------------  Taylor  -----
  complex(rp) elemental function Taylor ( x, y )

  !{ $w(z) = \sum_{n=0}^\infty \frac{(iz)^n}{\Gamma(\frac{n}2+1)}$.
  !  Separating even and odd terms, we have
  !  $w(z) = \sum_{n=0}^\infty \frac{(-z^2)^n}{n!} +
  !          \frac{2 i z}{\sqrt{\pi}}
  !            \sum_{n=0}^\infty \frac{(-z^2)^n}{1\cdot3\cdot\cdot\cdot(2n+1)}$
  !  The even terms are $e^{-z^2}$; the odd ones are erf$(iz)$.  We only use
  !  terms up to second order in each of the even and odd series, so this
  !  approximation gets seven digits only for $x, y \leq 0.06$.

    use MLSCommon, only: RP
    use Units, only: SqrtPi
    real(rp), intent(in) :: X, Y
    real(rp) :: U, V
    real(rp) :: X2, Y2, Y4, TR, TI
    real(rp), parameter :: TwoOvSqpi = 2.0_rp / sqrtPi

    real(rp), parameter :: Q0 = TwoOvSqpi, Q1 = 2.0 * TwoOvSqpi / 3.0
    real(rp), parameter :: Q2 = 4.0 * TwoOvSqpi / 15.0, Q3 = 2.0 * q2

    ! (x2,y2) = -z**2
    x2 = y*y - x*x
    y2 = -2.0_rp * x * y
    y4 = y2 * y2

    tr = q0 - q2 * y4 + ( q1 + q2 * x2 ) * x2
    ti = ( q1 + q3 * x2 ) * y2
    u = ( 1.0_rp + 0.5_rp * x2 ) * x2 + 1.0 - 0.5 * y4 - tr * y - ti * x
    v = ( x2 + 1.0_rp ) * y2 + x * tr - y * ti
    taylor = cmplx(u,v)

  end function Taylor

!{ \newpage

  ! ----------------------------------------  D_Real_Simple_Voigt  -----
  subroutine D_Real_Simple_Voigt ( x, y, dx, dy, u, du )
  ! Compute the real part of Fadeeva = Voigt (without Lorentz) and
  ! the real part of the derivative (which isn't just the derivative
  ! of Voigt).

!{ The Fadeeva function $w(z)$, where $z = x + i y$, can be written as
!  $V(x,y) + i L(x,y)$, where $V(x,y)$ is the Voigt function and
!  $L(x,y)$ is the Lorentz function. 
!
!  From 7.1.20 in {\bf Handbook of Mathematical Functions} by Abramowitz
!  and Stegun (National Bureau of Standards Applied Math Series 55) we have
!  $w^{\prime}(z) = \frac{2i}{\sqrt{\pi}} - 2 z w(z)$.\\
!  $\Re \, \frac{\partial w(z(t))}{\partial T} =
!   2 \left[ \left( -x V(x,y) + y L(x,y) \right) \frac{\partial x}{\partial T} +
!            \left( x L(x,y) + y V(x,y) +\frac1{\sqrt{\pi}} \right)
!             \frac{\partial y}{\partial T} \right]$.

    real(rp), intent(in) :: x, y, dx, dy
    real(rp), intent(out) :: u, du

    real(rp) :: v ! Lorentz function

    call simple_voigt ( x, y, u, v )
    du = 2.0_rp * ( (-x * u + y * v) * dx + (x * v + y * u - oneOvSpi) * dy )

  end subroutine D_Real_Simple_Voigt

  ! ---------------------------------------------  D_Simple_Voigt  -----
  subroutine D_Simple_Voigt ( x, y, dx, dy, u, v, du, dv )
  ! Compute Fadeeva = Voigt & Lorentz and its derivative

!{ Let $w^\prime(z) = a + i b$ and $z = x + i y$, and assume $x$ and $y$
!  are functions of some parameter, say $T$.  Then $\frac{d w(z)}{d T} =
!  \frac{d w(z)}{d z} \frac{d z}{d T} = a ^\prime - b y^\prime +
!  i ( a y^\prime + b x^\prime)$.

    real(rp), intent(in) :: x, y, dx, dy
    real(rp), intent(out) :: u, v, du, dv

    real(rp) :: a, b

    call simple_voigt ( x, y, u, v, a, b )
    du = a * dx - b * dy
    dv = a * dy + b * dx

  end subroutine D_Simple_Voigt

!=====================================================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, ModuleName(1:1)
  end function not_used_here

end module Voigt_M

! $Log$
! Revision 2.2  2004/04/26 22:07:12  vsnyder
! Reinstate 1-term asymptotic expansion in non-derivative case in Simple_Voigt
!
! Revision 2.1  2004/04/24 02:13:48  vsnyder
! Separate from slabs_sw_m, polish up some regions
!
! Revision 2.33  2004/04/20 00:48:06  vsnyder
! Only use Taylor really close to the origin
!
! Revision 2.32  2004/04/19 21:02:19  vsnyder
! Use Taylor instead of CDrayson near the origin
!
! Revision 2.31  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.30  2004/04/06 23:40:21  vsnyder
! Do slabs_prep where derivatives not requested in get_gl_slabs_arrays
! instead of doing nothing.
!
! Revision 2.29  2004/04/02 00:59:24  vsnyder
! Get catalog from slabs structure
!
! Revision 2.28  2004/03/30 02:25:08  vsnyder
! Comment out call to slabs_prep_dt until n1==0 etc are worked out
!
! Revision 2.27  2004/03/27 03:35:27  vsnyder
! Add pointer to catalog in slabs_struct.  Use it so as not to need to drag
! line centers and line widths around.  Write slabs_lines and slabswint_lines
! to get sum of beta over all lines; put slabs_struct instead of its components
! in the calling sequence.
!
! Revision 2.26  2004/03/20 03:17:44  vsnyder
! Steps along the way toward analytic temperature derivatives
!
! Revision 2.24  2003/07/09 22:46:24  vsnyder
! Futzing
!
! Revision 2.23  2003/07/08 00:09:18  vsnyder
! Inlined several functions
!
! Revision 2.22  2003/07/04 02:49:03  vsnyder
! Simplify interface to Get_GL_Slabs_Arrays
!
! Revision 2.21  2003/06/18 14:45:00  bill
! added subsetting feature for T-ders
!
! Revision 2.20  2003/06/13 21:28:20  bill
! fixed/improved some bugs with line shape derivative computations
!
! Revision 2.19  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.18  2003/05/16 23:53:05  livesey
! Now uses molecule indices rather than spectags
!
! Revision 2.17  2003/05/09 19:25:31  vsnyder
! Expect T+DT instead of T and DT separately in Get_GL_Slabs_Arrays
!
! Revision 2.16  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.15.2.2  2003/02/27 00:57:20  vsnyder
! Cosmetic changes, get rid of declared but unused variables
!
! Revision 2.15.2.1  2003/02/13 17:29:26  bill
! abs coeff obeys detailed balance
!
! Revision 2.15  2003/01/16 19:41:42  jonathan
! tested version: in 1D case, compute only each element along the LOS path before tangent point, and fill otherside accordingly
!
! Revision 2.14  2003/01/16 19:08:43  jonathan
! testing
!
! Revision 2.13  2003/01/16 18:50:20  jonathan
! For 1D FWM compute first element along the LOS path then fill other grid points with value of the first grid point
!
! Revision 2.12  2003/01/16 18:04:12  jonathan
! add Do_1D option to get_gl_slabs_arrays
!
! Revision 2.11  2003/01/10 21:55:26  vsnyder
! Move SpeedOfLight from Geometry to Units
!
! Revision 2.10  2002/12/20 20:22:59  vsnyder
! Cosmetic changes
!
! Revision 2.9  2002/12/03 00:34:23  vsnyder
! Test optional argument presence before using them
!
! Revision 2.8  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.7  2002/10/02 21:06:03  vsnyder
! Get SpeedOfLight from Geometry module
!
! Revision 2.6  2002/09/12 23:00:04  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.5  2002/08/05 17:51:15  jonathan
! debug
!
! Revision 2.4  2001/12/14 23:43:44  zvi
! Modification for Grouping concept
!
! Revision 2.3  2001/11/30 01:18:11  zvi
! Correcting a minor bug
!
! Revision 2.1  2001/10/17 22:01:00  zvi
! Eliminate computation of: ns
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.4.2.4  2001/09/12 21:38:54  zvi
! Added CVS stuff
!
! Revision 1.4.2.3  2001/09/12 00:05:44  livesey
! Corrected sign of velocity correction
!
! Revision 1.4.2.2  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2001/01/31 18:12:06  Z.Shippony
! Initial conversion to Fortran 90
