module SLABS
  use GH06, only: GW, GX
  use MLSCommon, only: R4, R8
  implicit NONE
  private
  public :: SLABS_PREP

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

contains

!---------------------------------------------------------------------
! This subroutine will compute a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
!
! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! =============================================     SLABS_PREP     =====
  Subroutine Slabs_prep(t, m, v0, el, w, ps, p, n, i, q, delta, gamma, &
 &                      n1, n2, v0s, x1, y, yi, slabs1, dslabs1)
!
! inputs:
!
    Real(r8), intent(in) :: t      ! Temperature K
    Real(r8), intent(in) :: m      ! Molecular mass amu
    Real(r8), intent(in) :: v0     ! Line center frequency MHz
    Real(r8), intent(in) :: el     ! Lower state energy cm-1
    Real(r8), intent(in) :: w      ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
    Real(r8), intent(in) :: ps     ! Pressure shift parameter in MHz/mbar
    Real(r8), intent(in) :: p      ! Pressure mbar
    Real(r8), intent(in) :: n      ! Temperature power dependence of w
    Real(r8), intent(in) :: i      ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
    Real(r8), intent(in) :: q(3)   ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
    Real(r8), intent(in) :: delta  ! Delta interference coefficient at
                                   ! 300K 1/mb
    Real(r8), intent(in) :: gamma  ! Gamma               "
    Real(r8), intent(in) :: n1     ! Temperature dependency of delta
    Real(r8), intent(in) :: n2     ! Temperature dependency of gamma
!
! outputs:
!
    Real(r8), intent(out) :: v0s      ! Pressure shifted line position
    Real(r8), intent(out) :: x1       ! Sqrt(Ln(2))/Doppler half width MHz
    Real(r8), intent(out) :: y        ! Sqrt(Ln(2))*collision width/doppler
                                      ! width
    Real(r8), intent(out) :: yi       ! Interference contribution
    Real(r8), intent(out) :: slabs1   ! Frequency independent piece of slabs
    Real(r8), intent(out) :: dslabs1  ! Derivative of slabs1 w.r.t. v0
!
!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.
!
! Internal constants:
!
!                          Boltzmann constant cm-1/K:
    Real(r8), parameter :: boltzcm = 0.6950387_r8
!                          Boltzmann constant MHz/K:
    Real(r8), parameter :: boltzmhz = 20836.74_r8
!                          sqrt(amu/K) used to calculate doppler width:
    Real(r8), parameter :: dc = 3.58117369e-7_r8
!                          converts intensity into absorption:
    Real(r8), parameter :: i2abs = 3.402136078e9_r8
    Real(r8), parameter :: loge = 4.34294481903251828e-1_r8  ! log10(e)
    Real(r8), parameter :: oned300 = 1.0_r8/300.0_r8         ! 1.0 / 300.0
    Real(r8), parameter :: sqrtln2 = 8.32554611157698e-1_r8  ! sqrt(ln(2))
!
! Internal data:
!
    Real(r8) :: betae, betav, de1, de2, ds, e1, e2, g, ns
    Real(r8) :: onedt, Wd, Q_Log, r, s, t3t
!
! The action begins here
!
    onedt = 1.0_r8 / t
    t3t = 300.0_r8 * onedt
    yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
    ns = 0.25_r8 + 1.5_r8 * n
    v0s = v0 + ps * p * (t3t**ns)
!
    betae = el / boltzcm
    betav = v0s / boltzmhz
    Wd = v0s * sqrt(t/m) * dc
    x1 = sqrtln2 / Wd
    y = x1 * w * p * (t3t**n)
    g = i - Q_Log(q,t) + loge *  betae * (oned300 - onedt)
    r = (i2abs * p * (10.0_r8**g)) / (t * Wd)
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
! ==================================================     Q_LOG     =====
  real(r8) Function Q_LOG ( Q0, T )
!
!   LOG OF PARTITION FUNCTION NORMALIZED TO Q(T=300K)
!
!      Returns the base 10 logarithm of the ratio of partition function
!              at temperature T to that at temperature 300K
!
!      Linear interpolation in log(T) and log(Q) is used.  JPL catalog
!              accuracy is maintained over temperature range 150-300K.
!
!      Inputs: Q0 = 3-element array of base 10 logarithms of partition
!                   function at T=300, 225, 150K from JPL catalog.
!                   [Poynter and Pickett, Appl. Opt. 24, 2235, July 1985]
!
!              T  = temperature [K]
!
!      Written by J.Waters.   JPL  2 June 1986
!****************************************************************************
!
    real(r8), intent(in) :: Q0(3), T
!
    real(r8) :: SLOPE, TLOG    
    real(r8), parameter :: tlog0(3) = & ! Log(T) for T = 300.0, 225.0, 150.0:
      (/ 2.47712125471966d0, 2.35218251811136d0, 2.17609125905568d0 /)
!
    tlog = dlog10(T)
!
    if (tlog < tlog0(2)) then
!
      slope = (q0(2)-q0(3))/(tlog0(2)-tlog0(3))
      Q_Log = q0(2)-q0(1)+slope*(tlog-tlog0(2))
!
    else
!
      slope = (q0(1)-q0(2))/(tlog0(1)-tlog0(2))
      Q_Log = slope*(tlog-tlog0(1))
!
    end if
!
    Return
  End Function Q_LOG
! ==============================================     SLABSWINT     =====
  real(r8) Function SLABSWINT ( DNU, V0S, X1, SLABS1, Y, YI )
! Compute slab with interference
    real(r8), intent(in) :: DNU, V0S, X1, SLABS1, Y, YI
    real(r8), Parameter :: sqrt_pi_i = 1.0_r8 / 1.7724538509055160273_r8
!
    real(r8) :: P, Q, R, U, V, W, X, Z
!
!  Note: dNu = v - v0s
!
! If the molecular transition and temperature have not changed but
! frequency has enter here.
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
  End Function SLABSWINT
! ================================================     Z_SLABS     =====
  Subroutine Z_SLABS ( X, Y, U, V )
!
!***********************************************************************
!                                                                      *
! Computes the Voigt function: Integral from - to + infinity of:       *
!                                                                      *
!     u = (y/Pi)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt        (Real(W(z)))    *
!                                                                      *
!   and:                                                               *
!                                                                      *
!     v = (1/Pi)*(x-t)*Exp(-t*t)/(y*y+(x-t)*(x-t)) dt  (Imag(W(z)))    *
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
    real(r8), intent(in) :: X, Y
    real(r8), intent(out) :: U, V
!
    real(r8) :: XA, YH
!
    xa = abs(x)
    if (xa+y > 6.0_r8) then             ! Region 4, Gauss-Hermite 6 points
      Call VoigtH6(xa,y,u,v)
    else
      yh = 1.0e-4_r8*xa*xa
      if (y <= yh) then
        if (xa > 5.0_r8) then
          Call VoigtH6(xa,y,u,v)       ! Region 4, Gauss-Hermite 6 points
        else
          Call Drayson(xa,y,u,v)       ! Region 1, Dowson+Taylor
        end if
      else                             ! Region 2, Hui (p=6)
        Call Hui6(xa,y,u,v)
      end if
    end if
!
    if (x < 0.0_r8) v = -v
!
    Return
  End Subroutine Z_SLABS
! ===============================================     VOIGTH6     =====
  Subroutine VOIGTH6 ( X, Y, U, V )
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
    real(r8), intent(in) :: X, Y
    real(r8), intent(out) :: U, V
!
    real(r8) :: FM, FP
    integer :: I
    real(r8), parameter :: PI_I = 1.0_r8 / 3.1415926535897932385_r8
    real(r8) :: SUMU, SUMV, T, XMT, XPT, Y2
!
    y2 = y * y
    Sumu = 0.0
    Sumv = 0.0
!
!  Note: The roots are both positive and negative !
!
    do i = 1, size(Gx)
      t = Gx(i)
      xpt = x + t
      xmt = x - t
      fm = 1.0 / (y2 + xmt * xmt)
      fp = 1.0 / (y2 + xpt * xpt)
      Sumu = Sumu + Gw(i) * (fm + fp)
      Sumv = Sumv + Gw(i) * (xmt * fm + xpt * fp)
    end do
!
    u = y * Sumu * Pi_i
    v = Sumv * Pi_i
!
    Return
  End Subroutine VOIGTH6
! ================================================     DRAYSON     =====
  Subroutine DRAYSON ( XIN, YIN, UO, VO )
!
    real(r8), intent(in) :: XIN, YIN
    real(r8), intent(out) :: UO, VO
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
    logical, save :: FIRST = .true.
    real(r4), save :: B(22)
    real(r4), parameter :: C(21) = &
    & (/ 0.7093602E-7, -0.2518434E-6, 0.8566874E-6, -0.2787638E-5 &
    &  , 0.8660740E-5, -0.2565551E-4, 0.7228775E-4, -0.1933631E-3 &
    &  , 0.4899520E-3, -0.1173267E-2, 0.2648762E-2, -0.5623190E-2 &
    &  , 0.1119601E-1, -0.2084976E-1, 0.3621573E-1, -0.5851412E-1 &
    &  , 0.8770816E-1, -0.121664, 0.15584, -0.184, 0.2 /)
    real(r4) :: CO
    real(r4), save :: D0(25), D1(25), D2(25), D3(25), D4(25)
    real(r4) :: DX
    real(r4) :: F, FD, FI, FR
    real(r4), parameter :: H = 0.201
    real(r4), save :: HN(25)
    integer :: I, J, N
    real(r4), save :: RI(15)
    real(r8), parameter :: TwoOvSPi = 1.1283791670955125739_r8 ! 2.0/Sqrt(Pi)
    real(r4) :: UU, VOIGT, VV, WI, WR
    real(r4) :: XL, YL,  Y2
!
    data B(1), B(2)/0.0, 0.7093602E-7/
!
    if (first) then
!
      first = .false.
!
!******* Initialize Region I. Compute Dawson's function at mesh points
!
      do i = 1, 15
        ri(i) = -i / 2.0
      end do
!
      do i = 1, 25
        hn(i) = h * (i - 0.5)
        co = 4.0 * hn(i) * hn(i) / 25.0 - 2.0
        do j = 2, 21
          b(j+1) = co * b(j) - b(j-1) + c(j)
        end do
        d0(i) = hn(i) * (b(22) - b(21)) / 5.0
        d1(i) = 1.0 - 2.0 * hn(i) * d0(i)
        d2(i) = (hn(i) *d1(i) + d0(i)) / ri(2)
        d3(i) = (hn(i) *d2(i) + d1(i)) / ri(3)
        d4(i) = (hn(i) *d3(i) + d2(i)) / ri(4)
      end do
!
    end if
!
    xl = xin
    yl = yin
!
!******** Region I. Compute Dawson's function at x from Taylor series
!
    j = xl / h
    n = min(j, 24)
    dx = xl - hn(n+1)
    f = (((d4(n+1)*dx + d3(n+1))*dx + d2(n+1))*dx + d1(n+1)) &
   &    * dx + d0(n+1)
!
    fd = 1.0 - 2.0 * xl * f
!
!  Taylor series expansion about y = 0.0
!
    y2=  yl*yl
    wi = exp(y2 - xl*xl)
    dx = -2.0*xl*yl
    wr = wi*cos(dx)
    wi = wi*sin(dx)
    fi = f
    uu = -yl
    fr = uu*fd
    j = 5.0 + (12.5-xl)*0.8*yl
    do j = 2, min(j, 14), 2
       f  = (xl*fd + f) / ri(j)
       fd = (xl*f + fd) / ri(j+1)
       uu = uu * yl
       fi = fi + f*uu
       uu = -uu * yl
       fr = fr + fd*uu
    end do
    voigt = wr + twoovspi * fr
    vv    = wi + twoovspi * fi
!
    uo = voigt
    vo = vv
!
    return
  end Subroutine DRAYSON
! ===================================================     HUI6     =====
  Subroutine HUI6 ( X, Y, U, V )
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
    real(r8), intent(in) :: X, Y
    real(r8), intent(out) :: U, V

    real(r8), parameter :: A(7) = &
    & (/ 122.607931777104326_r8, 214.382388694706425_r8, &
    &    181.928533092181549_r8,  93.155580458138441_r8, &
    &     30.180142196210589_r8,   5.912626209773153_r8, &
    &      0.564189583562615_r8 /)

    real(r8), parameter :: B(7) = &
    & (/ 122.607931773875350_r8, 352.730625110963558_r8, &
    &    457.334478783897737_r8, 348.703917719495792_r8, &
    &    170.354001821091472_r8,  53.992906912940207_r8, &
    &     10.479857114260399_r8 /)

    Complex(r8) :: W, WA, WB, Z
!
    z = cmplx(y,-x)
    wa = a(1)+z*(a(2)+z*(a(3)+z*(a(4)+z*(a(5)+z*(a(6)+z*a(7))))))
    wb = b(1)+z*(b(2)+z*(b(3)+z*(b(4)+z*(b(5)+z*(b(6)+z*(b(7)+z))))))
    w = wa / wb
    u = Real(w)
    v = aImag(w)
!
    Return
  End Subroutine HUI6
! ========================================     SLABS_PREP_WDER     =====
! ** ORIGINALLY: Subroutine Slabs_prep()
!
! This function will compute a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
!
! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! ** CHANGED: Jan/5/00  To Include derivatives of x1,y & slabs w.r.t. Nu0
!
  Subroutine SLABS_PREP_WDER(T, M, V0, EL, W, PS, P, N, I, Q, DELTA, GAMMA, &
  &                          N1, N2, V0S, X1, Y, YI, SLABS1, DX1_DV0, DY_DV0, DSLABS1_DV0)
!
! inputs:
!
    Real(r8), intent(in) :: T     ! Temperature K
    Real(r8), intent(in) :: M     ! Molecular mass amu
    Real(r8), intent(in) :: V0    ! Line center frequency MHz
    Real(r8), intent(in) :: EL    ! Lower state energy cm-1
    Real(r8), intent(in) :: W     ! Collision broadening parameter
                                  ! MHz/mbar at 300 K
    Real(r8), intent(in) :: PS    ! Pressure shift parameter in MHz/mbar
    Real(r8), intent(in) :: P     ! Pressure mbar
    Real(r8), intent(in) :: N     ! Temperature power dependence of w
    Real(r8), intent(in) :: I     ! Integrated spectral intensity
                                  ! Log(nm**2 MHz) at 300 K
    Real(r8), intent(in) :: Q(3)  ! Logarithm of the partition function
                                  ! At 300 , 225 , and 150 K
    Real(r8), intent(in) :: DELTA ! Delta interference coefficient at 300K 1/mb
    Real(r8), intent(in) :: GAMMA ! Gamma               "
    Real(r8), intent(in) :: N1    ! Temperature dependency of delta
    Real(r8), intent(in) :: N2    ! Temperature dependency of gamma
!
! outputs:
!
    Real(r8), intent(out) :: v0s         ! Pressure shifted line position
    Real(r8), intent(out) :: x1          ! Sqrt(Ln(2))/Doppler half width MHz
    Real(r8), intent(out) :: y           ! Sqrt(Ln(2))*collision width/doppler
                                         ! width
    Real(r8), intent(out) :: yi          ! Interference contribution
    Real(r8), intent(out) :: slabs1      ! Frequency independent piece of slabs
!
    Real(r8), intent(out) :: dx1_dv0     ! Derivative of x1 w.r.t. v0
    Real(r8), intent(out) :: dy_dv0      ! Derivative of y w.r.t. v0
    Real(r8), intent(out) :: dslabs1_dv0 ! Derivative of slabs1 w.r.t. v0
!
!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.
!
! Internal constants:
!
! boltzmann constant cm-1/K:
    Real(r8), parameter :: boltzcm = 0.6950387_r8
! boltzmann constant MHz/K:
    Real(r8), parameter :: boltzmhz = 20836.74_r8
! sqrt(amu/K) used to calculate doppler width:
    Real(r8), parameter :: dc = 3.58117369e-7_r8
! converts intensity into absorption:
    Real(r8), parameter :: i2abs = 3.402136078e9_r8
    Real(r8), parameter :: loge = 4.34294481903251828e-1_r8 ! log10(e)
    Real(r8), parameter :: oned300= 1.0_r8 / 300.0_r8       ! 1.0 / 300.0
    Real(r8), parameter :: sqrtln2 = 8.32554611157698e-1_r8 ! sqrt(ln(2))
!
! Internal data:
!
    Real(r8) :: BETAE, BETAV, DE1, DE2, DR_DV0, DS_DV0, DWD_DV0
    Real(r8) :: E1, E2, G, NS, ONEDT, Q_LOG, R, S, T3T, WD
!
! The action begins here
!
    onedt = 1.0_r8 / t
    t3t = 300.0_r8 * onedt
    yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
    ns = 0.25_r8 + 1.5_r8 * n
    v0s = v0 + ps * p * (t3t**ns)
!
    betae = el / boltzcm
    betav = v0s / boltzmhz
    Wd = v0s * sqrt(t/m) * dc
    x1 = sqrtln2 / Wd
    y = x1 * w * p * (t3t**n)
    g = i - Q_Log(q,t) + loge *  betae * (oned300 - onedt)
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
  End Subroutine SLABS_PREP_WDER
! ========================================     DVOIGT_SPECTRAL     =====
! Computes the voigt function and its first derivative with respect
! to spectaral parameters: w, n & Nu0
!
! NOTE: Before calling this routine, the user needs to call slabs_prep_wder()
!       routine to compute dx1_dv0,dy_dv0 and dslabs1_dNu0
!
  Subroutine DVOIGT_SPECTRAL ( DNU, NU0, X1, YI, Y, W, T, SLABS1, DX1_DV0, &
 &           DY_DV0, DSLABS1_DNU0, SWI, DSWI_DW, DSWI_DN, DSWI_DNU0 )

    real(r8), intent(in) :: DNU
    real(r8), intent(in) :: NU0
    real(r8), intent(in) :: X1
    real(r8), intent(in) :: YI
    real(r8), intent(in) :: Y
    real(r8), intent(in) :: W
    real(r8), intent(in) :: T
    real(r8), intent(in) :: SLABS1
    real(r8), intent(in) :: DX1_DV0
    real(r8), intent(in) :: DY_DV0
    real(r8), intent(in) :: DSLABS1_DNU0
    real(r8), intent(out) :: SWI
    real(r8), intent(out) :: DSWI_DW
    real(r8), intent(out) :: DSWI_DN
    real(r8), intent(out) :: DSWI_DNU0
!
    real(r8), parameter :: twovspi = 1.1283791670955125739_r8 ! 2.0/Sqrt(Pi)
    real(r8), parameter :: oneovspi = 0.5_r8 * twovspi

    real(r8) :: B, DB_DV0, DG_DV0, DQ_DV0, DR_DV0, DU_DV0, DU_DX, DU_DY
    real(r8) :: DV_DV0, DV_DX, DV_DY, DVVW_DV0, DX_DV0, DZ_DV0
    real(r8) :: G, Q, Q2, R, U, V, VVW, X, Z

    x = x1 * dNu
    Call Z_Slabs(x,y,u,v)
!
    du_dx = -2.0d0 * (x * u - y * v)
    du_dy =  2.0d0 * (x * v + y * u) - twovspi
!
    dv_dx = -du_dy         ! Cauchy-Riemann equation
    dv_dy =  du_dx         ! Cauchy-Riemann equation
!
!  Van Vleck - Wieskopf (VVW) line shape with Voigt
!
    q = 1.0d0 + dNu / Nu0
    q2 = q * q
!
    b = x1 * (2.0d0 * Nu0 + dNu)
    g = b * b + y * y
    z = (y - b * yi) / g
    r = z * oneovspi + yi * v
    vvw = (u + r) * q2
    SwI = slabs1 * vvw
!
! Compute the derivative of SwI w.r.t. w
!
    dSwI_dw = slabs1 * (y / w) * (du_dy +                 &
   &                q2 * oneovspi * (b*b-y*y) / (g*g) +   &
   &                q2 * yi * du_dx)
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
    db_dv0 = (2.0d0 * Nu0 + dNu)*dx1_dv0 + x1
    dg_dv0 = 2.0d0 * (b*db_dv0+y*dy_dv0)
    dz_dv0 = (dy_dv0-yi*db_dv0-z*dg_dv0)/g
    dr_dv0 = dz_dv0*oneovspi+yi*dv_dv0
    dvvw_dv0 = (du_dv0+dr_dv0)*q2 + 2.0d0*q*dq_dv0*(u+r)
    dSwI_dNu0 = dslabs1_dNu0*vvw + slabs1*dvvw_dv0
!
    Return
  End Subroutine DVOIGT_SPECTRAL
end module SLABS

! $Log$
