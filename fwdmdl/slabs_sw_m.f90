module SLABS_SW_M
  use MLSCommon, only: I4, R4, R8
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
!
 Subroutine VoigtH6(x,y,u,v)
!
!****************************************************************************
! Computes the complex Voigt function: (i/pi)*integral from - to + infinity *
! of : exp(-t*t)/(z -t) dt, where z=(x+iy), x, y >= 0. Using Gauss-Hermite  *
! Quadrature. The real part is u, the integral of:                          *
!                                                                           *
!                 (y/pi)*Exp(-t*t)/(y*y + (x-t)*(x-t)) dt                   *
!                                                                           *
! The imaginary part is v, the integral of:                                 *
!                                                                           *
!                 (1/pi)*(x-t)*Exp(-t*t)/(y*y + (x-t)*(x-t)) dt             *
!                                                                           *
!****************************************************************************
!
  Real(r8), INTENT(IN)  :: x,y
  Real(r8), INTENT(OUT) :: u,v
!
  Real(r8), Parameter :: Pi_i = 1.0_r8/3.1415926535897932385_r8   ! 1/Pi
!
! These are the 6-point-Gauss-Hermite abscissa (X-axis) values in [-1,1]:
! Note: The roots are both positive and negative !
!
  Integer(i4), Parameter :: N = 3
  Real(r8), Parameter :: Gx(N) = (/ &
     & 4.36077411927616508271e-1_r8, 1.33584907401369694957_r8, &
     & 2.35060497367449222280_r8 /)
!
  Real(r8), Parameter :: Gw(N) = (/ &
     & 7.24629595224392524608e-1_r8, 1.57067320322856644036e-1_r8, &
     & 4.53000990550884564224e-3_r8 /)

  Integer(i4) :: i
  Real(r8) :: y2, Sumv, Sumu, t, xpt, xmt, fm, fp
!
  y2 = y * y
  Sumu = 0.0_r8
  Sumv = 0.0_r8
!
  do i = 1, N
    t = Gx(i)
    xpt = x + t
    xmt = x - t
    fm = 1.0_r8 / (y2 + xmt * xmt)
    fp = 1.0_r8 / (y2 + xpt * xpt)
    Sumu = Sumu + Gw(i) * (fm + fp)
    Sumv = Sumv + Gw(i) * (xmt * fm + xpt * fp)
  end do
!
  u = y * Sumu * Pi_i
  v = Sumv * Pi_i
!
  Return
 End Subroutine VoigtH6
!
!-----------------------------------------------------------------------
!
 Subroutine Drayson(Xin,Yin,Uo,Vo)
!
  Real(r8),INTENT(IN)  :: Xin,Yin
  Real(r8),INTENT(OUT) :: Uo,Vo
!
  Real(r8), Parameter :: TwoOvSPi = 1.1283791670955125739_r8  ! 2.0/Sqrt(Pi)
!
!***********************************************************************
! Computes the voigt function: y/pi*integral from - to + infinity of   *
! exp(-t*t)/(y*y + (x-t)*(x-t)) dt, for x, y >= 0.  This algorithm was *
! obtained from "Rapid Computation of the Voigt Profile" by S. R.      *
! Drayson in the J. Quant. Spectrosc. Radiat. Transfer, Vol. 16,       *
! pp. 611-614, Pergamon Press 1976.                                    *
!                                                                      *
!               x = sqrt(ln 2) * (v - v0) / aD                         *
!               y = sqrt(ln 2) * aL / aD                               *
!                                                                      *
!  where v is the wave number, v0 is the line center wave number, aL   *
! is the Lorentzian line half-width and aD is the Doppler line         *
! half-width.                                                          *
!                                                                      *
!              voigt(0,0) = 1.                                         *
!              voigt(0,y) = 1 / (sqrt(pi) * y)   for large y           *
!                                                                      *
!***********************************************************************
!
  LOGICAL, SAVE :: FIRST = .True.

  REAL(r4), SAVE :: B(22)
  REAL(r4), SAVE :: RI(15)
  REAL(r4), SAVE :: HN(25), D0(25), D1(25), D2(25), D3(25), D4(25)
!
  REAL(r4), Parameter :: H = 0.201

  REAL(r4), Parameter :: C(21) = &
    & (/ 0.7093602E-7, -0.2518434E-6, 0.8566874E-6, -0.2787638E-5 &
    &  , 0.8660740E-5, -0.2565551E-4, 0.7228775E-4, -0.1933631E-3 &
    &  , 0.4899520E-3, -0.1173267E-2, 0.2648762E-2, -0.5623190E-2 &
    &  , 0.1119601E-1, -0.2084976E-1, 0.3621573E-1, -0.5851412E-1 &
    &  , 0.8770816E-1, -0.121664, 0.15584, -0.184, 0.2 /)

  REAL(r4) :: XL, YL, UU, VV, Y2, DX, CO, F, FD, FR, FI, WI, WR, VOIGT
!
  INTEGER :: I, J, MAXJ, N
!
  Data B(1), B(2)/0.0, 0.7093602E-7/
!
  IF(FIRST) THEN
!
    FIRST = .False.
!
!******* Initialize Region I. Compute Dawson's function at mesh points
!
    DO I = 1, 15
      RI(I) = -Real(I) / 2.0
    END DO
!
    DO I = 1, 25
      HN(I) = H * (I - 0.5)
      CO = 4.0 * HN(I) * HN(I) / 25.0 - 2.0
      DO J = 2, 21
        B(J+1) = CO * B(J) - B(J-1) + C(J)
      END DO
      D0(I) = HN(I) * (B(22) - B(21)) / 5.0
      D1(I) = 1.0 - 2.0 * HN(I) * D0(I)
      D2(I) = (HN(I) * D1(I) + D0(I)) / RI(2)
      D3(I) = (HN(I) * D2(I) + D1(I)) / RI(3)
      D4(I) = (HN(I) * D3(I) + D2(I)) / RI(4)
    END DO
!
  ENDIF
!
  XL = Xin
  YL = Yin
!
!******** Region I. Compute Dawson's function at x from Taylor series
!
  J = XL / H
  N = MIN0(J, 24)
  DX = XL - HN(N+1)
  F = (((D4(N+1)*DX + D3(N+1))*DX + D2(N+1))*DX + D1(N+1))*DX + D0(N+1)
!
  FD = 1.0 - 2.0 * XL * F
!
!  Taylor series expansion about y = 0.0
!
  Y2=  YL*YL
  WI = EXP(Y2 - XL*XL)
  DX = -2.0*XL*YL
  WR = WI*COS(DX)
  WI = WI*SIN(DX)
  FI = F
  UU = -YL
  FR = UU*FD
  J = 5.0 + (12.5-XL)*0.8*YL
  MAXJ = MIN0(J, 14)
  DO J = 2, MAXJ, 2
     F  = (XL*FD + F) / RI(J)
     FD = (XL*F + FD) / RI(J+1)
     UU = UU * YL
     FI = FI + F*UU
     UU = -UU * YL
     FR = FR + FD*UU
  END DO

  VOIGT = WR + TwoOvSpi * FR
  VV    = WI + TwoOvSPi * FI
!
  Uo = VOIGT
  Vo = VV
!
  RETURN

END Subroutine Drayson
!
!---------------------------------------------------------------------
!
 Subroutine Hui6(x,y,u,v)
!
!***********************************************************************
! Computes the Voigt function: y/pi*integral from - to + infinity of   *
! exp(-t*t)/(y*y + (x-t)*(x-t)) dt, for x, y >= 0.  This algorithm was *
! obtained from "Rapid computation of the Voigt and complex error      *
! function"  by: A.K. Hui, B.H. Aramstrong and A.A. Wray               *
! J. Quant. Spectrosc. Radiat. Transfer, Vol 19, pp. 506-516,          *
! Pergamon Press 1978.                                                 *
!                                                                      *
!               x = sqrt(ln 2) * (v - v0) / aD                         *
!               y = sqrt(ln 2) * aL / aD                               *
!                                                                      *
!  where v is the wave number, v0 is the line center wave number, aL   *
! is the Lorentzian line half-width and aD is the Doppler line         *
! half-width.                                                          *
!                                                                      *
!***********************************************************************
!
  Real(r8), INTENT(IN)  :: x,y
  Real(r8), INTENT(OUT) :: u,v

  Complex(r8) :: z, w, wa, wb
!
  Real(r8), Parameter :: a(7) = (/ &
     &     122.607931777104326_r8, 214.382388694706425_r8, &
     &     181.928533092181549_r8,  93.155580458138441_r8, &
     &      30.180142196210589_r8,   5.912626209773153_r8, &
     &       0.564189583562615_r8 /)
!
  Real(r8), Parameter :: b(7) = (/ &
     &     122.607931773875350_r8, 352.730625110963558_r8, &
     &     457.334478783897737_r8, 348.703917719495792_r8, &
     &     170.354001821091472_r8,  53.992906912940207_r8, &
     &      10.479857114260399_r8 /)
!
  z = Cmplx(y,-x)
  wa = a(1)+z*(a(2)+z*(a(3)+z*(a(4)+z*(a(5)+z*(a(6)+z*a(7))))))
  wb = b(1)+z*(b(2)+z*(b(3)+z*(b(4)+z*(b(5)+z*(b(6)+z*(b(7)+z))))))
  w = wa / wb
  u = Real(w)
  v = AImag(w)
!
  Return
 End Subroutine Hui6
!
!---------------------------------------------------------------------
!
 Subroutine Z_Slabs(x,y,u,v)
!
!***********************************************************************
!                                                                      *
! Computes the Voigt function: Integral from - to + infinity of:       *
!                                                                      *
!     u = (y/Pi)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt        (Real(W(z)))    *
!                                                                      *
!   and:                                                               *
!                                                                      *
!     v = (1/Pi)*(x-t)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt  (Image(W(z)))   *
!                                                                      *
!   Here:                                                              *
!               x = sqrt(ln 2) * (v - v0) / aD    (x >= 0.0)           *
!               y = sqrt(ln 2) * aL / aD          (y >= 0.0)           *
!                                                                      *
!   Where v is the wave number, v0 is the line center wave number, aL  *
! is the Lorentzian line half-width and aD is the Doppler line         *
! half-width.                                                          *
!                                                                      *
!***********************************************************************
!
  Real(r8), INTENT(IN)  :: x,y
  Real(r8), INTENT(OUT) :: u,v

  Real(r8) :: xa, yh
!
  xa = abs(x)
  if(xa+y > 6.0_r8) then           ! Region 4, Gauss-Hermite 6 points
    Call VoigtH6(xa,y,u,v)
  else
    yh = 1.0e-4*xa*xa
    if(y <= yh) then
      if(xa > 5.0_r8) then
        Call VoigtH6(xa,y,u,v)    ! Region 4, Gauss-Hermite 6 points
      else
        Call Drayson(xa,y,u,v)    ! Region 1, Dowson+Taylor
      endif
    else                          ! Region 2, Hui (p=6)
      Call Hui6(xa,y,u,v)
    endif
  endif
!
  if(x < 0.0_r8) v = -v
!
  Return
 End Subroutine Z_Slabs
!
!---------------------------------------------------------------------
! This function will compute a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
!
! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
!
 Subroutine Slabs_prep(t,m,v0,el,w,ps,p,n,i,q,delta,gamma,n1,n2, &
                    &  v0s,x1,y,yi,slabs1,dslabs1)
!
  implicit none
!
! inputs:
!
  Real(r8), INTENT(IN) :: t        ! Temperature K
  Real(r8), INTENT(IN) :: m        ! Molecular mass amu
  Real(r8), INTENT(IN) :: v0       ! Line center frequency MHz
  Real(r8), INTENT(IN) :: el       ! Lower state energy cm-1
  Real(r8), INTENT(IN) :: w        ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
  Real(r8), INTENT(IN) :: ps       ! Pressure shift parameter in MHz/mbar
  Real(r8), INTENT(IN) :: p        ! Pressure mbar
  Real(r8), INTENT(IN) :: n        ! Temperature power dependence of w
  Real(r8), INTENT(IN) :: i        ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
  Real(r8), INTENT(IN) :: q(3)     ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
  Real(r8), INTENT(IN) :: delta    ! Delta interference coefficient at 300K 1/mb
  Real(r8), INTENT(IN) :: gamma    ! Gamma               "
  Real(r8), INTENT(IN) :: n1       ! Temperature dependency of delta
  Real(r8), INTENT(IN) :: n2       ! Temperature dependency of gamma
!
! outputs:
!
  Real(r8), INTENT(OUT) :: v0s     ! Pressure shifted line position
  Real(r8), INTENT(OUT) :: x1      ! Sqrt(Ln(2))/Doppler half width MHz
  Real(r8), INTENT(OUT) :: y       ! Sqrt(Ln(2))*collision width /
                                   !             doppler width
  Real(r8), INTENT(OUT) :: yi      ! Interference contribution
  Real(r8), INTENT(OUT) :: slabs1  ! Frequency independent piece of slabs
  Real(r8), INTENT(OUT) :: dslabs1 ! Derivative of slabs1 w.r.t. v0
!
!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.
!
! Internal constants:
!
!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K
!  sqrtln2  - sqrt(ln(2))
!  loge     - log10(e)
!
  Real(r8), Parameter :: i2abs = 3.402136078e9_r8
  Real(r8), Parameter :: dc = 3.58117369e-7_r8
  Real(r8), Parameter :: boltzcm = 0.6950387_r8
  Real(r8), Parameter :: boltzmhz = 20836.74_r8
  Real(r8), Parameter :: sqrtln2 = 8.32554611157698e-1_r8
  Real(r8), Parameter :: loge = 4.34294481903251828e-1_r8
  Real(r8), Parameter :: oned300 = 1.0_r8/300.0_r8
!
  Real(r8), parameter :: tl1 = 1.76091259055681e-1_r8     ! Log10(225/150)
  Real(r8), parameter :: tl2 = 1.24938736608300e-1_r8     ! Log10(300/225)
!
! Internal data:
!
  Real(r8) :: Wd, Q_Log, betae, betav, t3t, onedt, ns, r, e1, e2, de1, &
              de2, g, s, ds
!
! The action begins here
!
  onedt = 1.0_r8 / t
  t3t = 300.0_r8 * onedt
  yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
  if (t < 225.0_r8) then
    r = (q(2)-q(3))/tl1
    Q_Log = q(2)-q(1)+r*Log10(t/225.0_r8)
  else
    r = (q(1)-q(2))/tl2
    Q_Log = r*Log10(t/300.0_r8)
  endif
!
  ns = 0.25_r8 + 1.5_r8 * n
  v0s = v0 + ps * p * (t3t**ns)
!
  betae = el / boltzcm
  betav = v0s / boltzmhz
  Wd = v0s * Sqrt(t/m) * dc
  x1 = sqrtln2 / Wd
  y = x1 * w * p * (t3t**n)
  g = i - Q_Log + loge *  betae * (oned300 - onedt)
  r = (i2abs * p * (10.0**g)) / (t * Wd)
  e1 = exp(-betav*onedt)
  e2 = exp(-betav*oned300)
  de1 = -e1*onedt/boltzmhz
  de2 = -e2*oned300/boltzmhz
  g = 1.0_r8 - e2
  s = (1.0_r8 - e1) / g
  ds = (-de1*g+(1.0_r8-e1)*de2)/(g*g)
  slabs1 = r * s
  dslabs1 = r * ds
!
  Return
 End Subroutine Slabs_prep
!
!---------------------------------------------------------------------
! Slabs_prep_wder: ** ORIGINALLY: Subroutine Slabs_prep()
!
! This function will compute a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
!
! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! ** CHANGED: Jan/5/00  To Include derivatives of x1,y & slabs w.r.t. Nu0
!
 Subroutine Slabs_prep_wder(t,m,v0,el,w,ps,p,n,i,q,delta,gamma,n1,n2, &
     &      v0s,x1,y,yi,slabs1,dx1_dv0,dy_dv0,dslabs1_dv0)
!
! inputs:
!
  Real(r8), INTENT(IN) :: t        ! Temperature K
  Real(r8), INTENT(IN) :: m        ! Molecular mass amu
  Real(r8), INTENT(IN) :: v0       ! Line center frequency MHz
  Real(r8), INTENT(IN) :: el       ! Lower state energy cm-1
  Real(r8), INTENT(IN) :: w        ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
  Real(r8), INTENT(IN) :: ps       ! Pressure shift parameter in MHz/mbar
  Real(r8), INTENT(IN) :: p        ! Pressure mbar
  Real(r8), INTENT(IN) :: n        ! Temperature power dependence of w
  Real(r8), INTENT(IN) :: i        ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
  Real(r8), INTENT(IN) :: q(3)     ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
  Real(r8), INTENT(IN) :: delta    ! Delta interference coefficient at 300K 1/mb
  Real(r8), INTENT(IN) :: gamma    ! Gamma               "
  Real(r8), INTENT(IN) :: n1       ! Temperature dependency of delta
  Real(r8), INTENT(IN) :: n2       ! Temperature dependency of gamma
!
! outputs:
!
  Real(r8), INTENT(OUT) :: v0s       ! Pressure shifted line position
  Real(r8), INTENT(OUT) :: x1        ! Sqrt(Ln(2))/Doppler half width MHz
  Real(r8), INTENT(OUT) :: y         ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
  Real(r8), INTENT(OUT) :: yi        ! Interference contribution
  Real(r8), INTENT(OUT) :: slabs1    ! Frequency independent piece of slabs
!
  Real(r8), INTENT(OUT) :: dx1_dv0       ! Derivative of x1 w.r.t. v0
  Real(r8), INTENT(OUT) :: dy_dv0        ! Derivative of y w.r.t. v0
  Real(r8), INTENT(OUT) :: dslabs1_dv0   ! Derivative of slabs1 w.r.t. v0
!
!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.
!
! Internal constants:
!
!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K
!  sqrtln2  - sqrt(ln(2))
!  loge     - log10(e)
!
  Real(r8), Parameter :: i2abs = 3.402136078e9_r8
  Real(r8), Parameter :: dc = 3.58117369e-7_r8
  Real(r8), Parameter :: boltzcm = 0.6950387_r8
  Real(r8), Parameter :: boltzmhz = 20836.74_r8
  Real(r8), Parameter :: sqrtln2 = 8.32554611157698e-1_r8
  Real(r8), Parameter :: loge = 4.34294481903251828e-1_r8
  Real(r8), Parameter :: oned300 = 1.0_r8/300.0_r8
!
  Real(r8), parameter :: tl1 = 1.76091259055681e-1_r8     ! Log10(225/150)
  Real(r8), parameter :: tl2 = 1.24938736608300e-1_r8     ! Log10(300/225)
!
! Internal data:
!
  Real(r8) Wd, Q_Log, betae, betav, t3t, onedt, ns, r, e1, e2, &
     &     de1, de2, g, s, dWd_dv0, dr_dv0, ds_dv0
!
! The action begins here
!
  onedt = 1.0_r8 / t
  t3t = 300.0_r8 * onedt
  yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
  if (t < 225.0_r8) then
    r = (q(2)-q(3))/tl1
    Q_Log = q(2)-q(1)+r*Log10(t/225.0_r8)
  else
    r = (q(1)-q(2))/tl2
    Q_Log = r*Log10(t/300.0_r8)
  endif
!
  ns = 0.25_r8 + 1.5_r8 * n
  v0s = v0 + ps * p * (t3t**ns)
!
  betae = el / boltzcm
  betav = v0s / boltzmhz
  Wd = v0s * Dsqrt(t/m) * dc
  x1 = sqrtln2 / Wd
  y = x1 * w * p * (t3t**n)
  g = i - Q_Log + loge *  betae * (oned300 - onedt)
  r = (i2abs * p * (10.0**g)) / (t * Wd)
  e1 = exp(-betav*onedt)
  e2 = exp(-betav*oned300)
  de1 = -e1*onedt/boltzmhz
  de2 = -e2*oned300/boltzmhz
  g = 1.0_r8 - e2
  s = (1.0_r8 - e1) / g
  ds_dv0 = (-de1*g+(1.0_r8-e1)*de2)/(g*g)
  slabs1 = r * s
!
  dWd_dv0 = Wd/v0s
  dr_dv0 = -r*dWd_dv0/Wd
!
  dx1_dv0 = -x1*dWd_dv0/Wd
  dy_dv0 = y*dx1_dv0/x1
  dslabs1_dv0 = r * ds_dv0 + s * dr_dv0
!
  Return
!
 End Subroutine Slabs_prep_wder
!
!---------------------------------------------------------------------
!
 Real(r8) Function Slabswint(dNu,v0s,x1,slabs1,y,yi)
!
  Real(r8), INTENT(IN) :: dNu, v0s, x1, slabs1, y, yi
!
  Real(r8), Parameter :: sqrt_pi_i = 1.0_r8/1.7724538509055160273_r8
!
!  Note: dNu = v - v0s
!
! If the molecular transition and temperature have not changed but
! frequency has enter here.
!
! inputs: dNu , x1 , slabs1 , y, v0s, yi
! output: slabswint (slab with interference)

  Real(r8) :: x, u, v, q, p, z, r, w
!
  x = x1 * dNu
  Call Z_Slabs(x,y,u,v)
!
!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)
!
  q = 1.0_r8 + dNu / v0s
  p = x1 * (2.0_r8 * v0s + dNu)
  z = (y - p * yi)/(p * p + y * y)
  r = z * sqrt_pi_i + yi * v
  w = (u + r) * q * q
  Slabswint = slabs1 *  w
!
  Return
 End Function Slabswint
!
!---------------------------------------------------------------------
! Computes the voigt function and its first derivative with respect
! to spectaral parameters: w, n & Nu0
!
! NOTE: Before calling this routine, the user needs to call slabs_prep_wder()
!       routine to compute dx1_dv0,dy_dv0 and dslabs1_dNu0
!
 Subroutine dvoigt_spectral(dNu,Nu0,x1,yi,y,w,t,slabs1,dx1_dv0, &
     &      dy_dv0,dslabs1_dNu0,SwI,dSwI_dw,dSwI_dn,dSwI_dNu0)
!
  Real(r8), Parameter :: twovspi = 1.1283791670955125739_r8   ! 2.0/Sqrt(Pi)
  Real(r8), Parameter :: oneovspi = 0.5_r8 * twovspi
!
  Real(r8), INTENT(IN) :: dNu, Nu0, x1, yi, y, w, t, slabs1, &
                          dslabs1_dNu0, dx1_dv0, dy_dv0

  Real(r8), INTENT(OUT) :: SwI,dSwI_dw,dSwI_dn,dSwI_dNu0
!
  Real(r8) :: x, u, v, du_dx, du_dy, dv_dx, dv_dy, q, q2, b, g, z, r
!
  Real(r8) :: dq_dv0, dx_dv0, du_dv0, dv_dv0, db_dv0, dg_dv0, dz_dv0, &
              dr_dv0, dvvw_dv0, vvw
!
  x = x1 * dNu
  Call Z_Slabs(x,y,u,v)
!
  du_dx = -2.0_r8 * (x * u - y * v)
  du_dy =  2.0_r8 * (x * v + y * u) - twovspi
!
  dv_dx = -du_dy         ! Cauchy-Riemann equation
  dv_dy =  du_dx         ! Cauchy-Riemann equation
!
!  Van Vleck - Wieskopf (VVW) line shape with Voigt
!
  q = 1.0_r8 + dNu / Nu0
  q2 = q * q
!
  b = x1 * (2.0_r8 * Nu0 + dNu)
  g = b * b + y * y
  z = (y - b * yi) / g
  r = z * oneovspi + yi * v
  vvw = (u + r) * q2
  SwI = slabs1 * vvw
!
! Compute the derivative of SwI w.r.t. w
!
  dSwI_dw = slabs1 * (y / w) * (du_dy +          &
     &      q2 * oneovspi * (b*b-y*y) / (g*g) +  &
     &      q2 * yi * du_dx)
!
! Compute the derivative of SwI w.r.t. n
!
  dSwI_dn = slabs1 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)
!
! Finaly, compute the derivative of SwI w.r.t. Nu0
!
! ***** Analytically *****
!
  dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
  dx_dv0 = dNu * dx1_dv0 - x1
  du_dv0 = du_dx * dx_dv0 + du_dy * dy_dv0
  dv_dv0 = dv_dx * dx_dv0 + dv_dy * dy_dv0
  db_dv0 = (2.0_r8 * Nu0 + dNu)*dx1_dv0 + x1
  dg_dv0 = 2.0_r8 * (b*db_dv0+y*dy_dv0)
  dz_dv0 = (dy_dv0-yi*db_dv0-z*dg_dv0)/g
  dr_dv0 = dz_dv0*oneovspi+yi*dv_dv0
  dvvw_dv0 = (du_dv0+dr_dv0)*q2 + 2.0_r8*q*dq_dv0*(u+r)
  dSwI_dNu0 = dslabs1_dNu0*vvw + slabs1*dvvw_dv0
!
  Return
 End Subroutine dvoigt_spectral
!
End module SLABS_SW_M
! $Log$
! Revision 1.1  2001/01/31 18:12:06  Z.Shippony
! Initial conversion to Fortran 90
