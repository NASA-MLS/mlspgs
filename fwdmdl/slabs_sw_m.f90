module SLABS_SW_M

  use MLSCommon, only: R8, RP, IP

  use SpectroscopyCatalog_m, only: CATALOG_T, Lines

  Implicit NONE

  Private
  Public :: dvoigt_spectral, slabs, slabswint, voigt_lorentz, &
        &  real_simple_voigt, simple_voigt, rlorentz, rvoigth2, &
        &  cvoigth2, rvoigth6, cvoigth6, rhui6, chui6, rdrayson, &
        &  cdrayson, slabs_prep, slabs_prep_arrays, get_gl_slabs_arrays

  REAL(rp), PARAMETER :: OneOvSPi = 0.56418958354775628695_rp  ! 1.0/Sqrt(Pi)
!
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
CONTAINS
!---------------------------------------------------------------------------
! Computes the voigt function and its first derivative with respect
! to spectaral parameters: w, n & Nu0
!
! NOTE: Before calling this routine, the user needs to call slabs_prep_wder()
!       routine to compute dslabs1_dNu0
!
 Subroutine dvoigt_spectral(dNu,Nu0,x1,yi,y,w,t,slabs1,SwI, &
                         &  dslabs1_dNu0,dSwI_dw,dSwI_dn,dSwI_dNu0)
!
  REAL(rp), INTENT(IN) :: dnu, nu0, x1, yi, y, w, t, slabs1
  Real(rp), INTENT(IN), OPTIONAL :: dslabs1_dNu0

  Real(rp), INTENT(OUT) :: SwI, dSwI_dw,dSwI_dn,dSwI_dNu0
!
  REAL(rp) :: x, u, v, du_dx, du_dy, dv_dx, dv_dy, q, q2, b, g, z, r
!
  Real(rp) :: dq_dv0, dx_dv0, du_dv0, dv_dv0, db_dv0, dg_dv0, dz_dv0, &
              dr_dv0, dvvw_dv0, vvw
!
  x = x1 * dNu
  Call simple_voigt(x,y,u,v)
!
!  Van Vleck - Wieskopf (VVW) line shape with Voigt
!
  q = 1.0_rp + dNu / Nu0
  q2 = q * q
!
  b = x1 * (2.0_r8 * Nu0 + dNu)
  g = b * b + y * y
  z = (y - b * yi) / g
  r = z * OneOvSPi + yi * v
  vvw = (u + r) * q2
  SwI = slabs1 * vvw
!
  du_dx = 2.0_rp * (y * v - x * u)
  du_dy = 2.0_rp * (y * u + x * v - OneOvSPi)
!
  dv_dx = -du_dy         ! Cauchy-Riemann equation
  dv_dy =  du_dx         ! Cauchy-Riemann equation
!
! Compute the derivative of SwI w.r.t. w
!
  dSwI_dw = q2 * slabs1* (y/w) * (du_dy + yi*du_dx + &
                              &   OneOvSPi*((1.0_rp-2.0_rp*z*y)/g))
!
! Compute the derivative of SwI w.r.t. n
!
  dSwI_dn = q2 * slabs1 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)
!
! Finaly, compute the derivative of SwI w.r.t. Nu0
!
! ***** Analytically *****
!
  dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
  dx_dv0 = -x1
  du_dv0 = du_dx * dx_dv0
  dv_dv0 = dv_dx * dx_dv0
  db_dv0 = x1
  dg_dv0 = 2.0_rp * b*db_dv0
  dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
  dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
  dvvw_dv0 = (du_dv0+dr_dv0)*q2 + 2.0_rp*q*dq_dv0*(u+r)
  dSwI_dNu0 = dslabs1_dNu0*vvw + slabs1*dvvw_dv0
!
  Return
 End Subroutine dvoigt_spectral
!
!---------------------------------------------------------------------
!
 Real(rp) Function Slabs(dNu,v0s,x1,slabs1,y)
!
  REAL(rp), INTENT(IN) :: dNu, v0s, x1, slabs1, y
!
!  Note: dNu = v - v0s
!
! If the molecular transition and temperature have not changed but
! frequency has enter here.
!
! inputs: dNu , x1 , slabs1 , y, v0s, yi
! output: slabs
!
  REAL(rp) :: u
!
  Call real_simple_voigt(x1*dNu,y,u)
!
!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)
!
  Slabs = slabs1 * (1.0_rp + dNu / v0s)**2 &
            * (u + OneOvSPi*y/((x1*(2.0_rp*v0s+dNu))**2 + y*y))
!
  Return
 End Function Slabs
!
!---------------------------------------------------------------------
!
 Real(rp) Function Slabswint(dNu,v0s,x1,slabs1,y,yi)
!
  REAL(rp), INTENT(IN) :: dNu, v0s, x1, slabs1, y, yi
!
!  Note: dNu = v - v0s
!
! If the molecular transition and temperature have not changed but
! frequency has enter here.
!
! inputs: dNu , x1 , slabs1 , y, v0s, yi
! output: slabswint (slab with interference)

  Real(rp) :: x, u, q, p, z, y2, w
!
  x = x1 * dNu
  Call real_simple_voigt(x,y,u)
!
!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)
!
  q = (1.0_rp + dNu / v0s)**2
  p = x1 * (2.0_rp * v0s + dNu)
  y2 = y*y
  z = OneOvSPi*((y - p*yi)/(p*p + y2) + yi*x/(x*x+y2))
  w = (u + z) * q
  Slabswint = slabs1 *  w
!
  Return
 End Function Slabswint
!
!---------------------------------------------------------------------------
! Computes the Voigt/Lorentz function and its first derivative with respect
! to spectaral parameters: w, n & Nu0
!
! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0
!
 Subroutine Voigt_Lorentz(dNu, Nu0, x1, yi, y, w, t, slabs1, VL,  &
                      &   dslabs1_dNu0, dVL_dw, dVL_dn, dVL_dNu0)
!
  Real(rp), INTENT(IN) :: dNu, Nu0, x1, yi, y, w, t, slabs1, dslabs1_dNu0

  Real(rp), INTENT(OUT) :: VL, dVL_dw, dVL_dn, dVL_dNu0
!
  Real(rp) :: xj, zj, q, y2, q2, u, v, Sum, up1, up2, dn1, dn2, dup1, &
   &          dup2, ddn1, ddn2, dy_dw, dy_dn, dq_dNu0, dSum_dw, dSum_dn, &
   &          dSum_dNu0
!
  q = 1.0_rp + dNu / Nu0
  q2 = q * q
!
  y2 = y * y
  xj = x1 * dNu
  zj = x1 * (2.0_r8 * Nu0 + dNu)
  dn1 = zj * zj + y2
  up1 = y - zj * yi
!
! Van Vleck - Wieskopf (VVW) line shape with Voigt
!
  Call simple_voigt(xj,y,u,v)
  dup1 = up1 * OneOvSPi / dn1 + yi * v
  ddn2 = u + dup1
  Sum = ddn2 / OneOvSPi
  VL = slabs1 * ddn2 * q2            ! This is the Voigt + VVW correction
!
  dn2 = xj * xj + y2
  up2 = y - yi * xj
!
  dy_dw = y / w
  dup1 = dy_dw
  ddn1 = 2.0 * y * dy_dw
  dup2 = dy_dw
  ddn2 = 2.0 * y * dy_dw
!
  dSum_dw = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
 &          (dn2*dup2-up2*ddn2)/(dn2*dn2)
!
  dVL_dw = OneOvSPi * slabs1 * q2 * dSum_dw
!
  dy_dn = y * Log(300.0/t)
  dup1 = dy_dn
  ddn1 = 2.0 * y * dy_dn
  dup2 = dy_dn
  ddn2 = 2.0 * y * dy_dn
  dSum_dn = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
 &          (dn2*dup2-up2*ddn2)/(dn2*dn2)

  dVL_dn = OneOvSPi * slabs1 * q2 * dSum_dn
!
  dup2 =  yi * x1               !  x1 = -dxj_dNu0
  ddn2 = -2.0 * xj * x1         !  x1 = -dxj_dNu0
  dSum_dNu0 = (dn2*dup2-up2*ddn2)/(dn2*dn2)
  dq_dNu0 = -(dNu+Nu0)/(Nu0*Nu0)
!
  dVL_dNu0 = OneOvSPi * (dslabs1_dNu0 * q2 * Sum          + &
            &         2.0 * slabs1 * q * dq_dNu0 * Sum + &
            &         slabs1 * q2 * dSum_dNu0)
!
  Return
 End Subroutine Voigt_Lorentz
!
!---------------------------------------------------------------------------
! simple REAL(Voigt) function
!
  ELEMENTAL SUBROUTINE real_simple_voigt(x,y,u)
!
! inputs
!
  REAL(rp), INTENT(in) :: x ! doppler width, del frequency ratio
  REAL(rp), INTENT(in) :: y ! doppler width, collision width ratio
!
! outputs
!
  REAL(rp), INTENT(out) :: u ! real part of Voigt
!
! internals
!
  REAL(rp), PARAMETER :: xl=5.2_rp, yl=0.05_rp, yh=0.6_rp, dydx=(yh-yl)/xl
  REAL(rp) :: xa
!
! This is sorted in likely occurance of each case
!
  xa = ABS(x)
!
! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed
!
  IF(y + 0.666666*xa > 100.0_rp) THEN
!
    u = rlorentz(xa,y)
!
  ELSE IF(y + 0.6875_rp * xa > 11.0_rp) THEN
!
! Drayson's quick 2pt hermite integral (essentially a lorentz)
!
    u = rvoigth2(xa,y)
!
  ELSE IF(y > 0.6_rp .OR. y > yl + dydx*xa) THEN
!
! Intermediate region
!
    u = rhui6(xa,y)
!
  ELSE IF(xa > 5.2_rp) THEN
!
! small y large x limit
!
    u = rvoigth6(xa,y)
!
  ELSE
!
! Near line center where Doppler dominates Pressure broadening
!
    u = rdrayson(xa,y)
!
  ENDIF
!
  END SUBROUTINE real_simple_voigt
!
!---------------------------------------------------------------------------
! simple Voigt function
!
  ELEMENTAL SUBROUTINE simple_voigt(x,y,u,v)
!
! inputs
!
  REAL(rp), INTENT(in) :: x ! doppler width, del frequency ratio
  REAL(rp), INTENT(in) :: y ! doppler width, collision width ratio
!
! outputs
!
  REAL(rp), INTENT(out) :: u ! real part of Voigt
  REAL(rp), OPTIONAL, INTENT(out) :: v ! imaginary part of Voigt
!
! internals
!
  REAL(rp), PARAMETER :: xl=5.2_rp, yl=0.05_rp, yh=0.6_rp, dydx=(yh-yl)/xl
  REAL(rp) :: xa
  COMPLEX(rp) :: uv
!
! This is sorted in likely occurance of each case
!
  xa = ABS(x)
!
! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed
!
  IF(y + 0.666666*xa > 100.0_rp) THEN
!
    uv = clorentz(xa,y)
!
  ELSE IF(y + 0.6875_rp * xa > 11.0_rp) THEN
!
! Drayson's quick 2pt hermite integral (essentially a lorentz)
!
    uv = cvoigth2(xa,y)
!
  ELSE IF(y > 0.6_rp .OR. y > yl + dydx*xa) THEN
!
! Intermediate region
!
    uv = chui6(xa,y)
!
  ELSE IF(xa > 5.2_rp) THEN
!
! small y large x limit
!
    uv = cvoigth6(xa,y)
!
  ELSE
!
! Near line center where Doppler dominates Pressure broadening
!
    uv = cdrayson(xa,y)
!
  ENDIF
!
  u = REAL(uv,KIND=rp)
  v = SIGN(AIMAG(uv),x)
!
  END SUBROUTINE simple_voigt
!
!---------------------------------------------------------------------------
! Real Lorentz
!
  REAL(rp) ELEMENTAL FUNCTION rlorentz(x,y)
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  rlorentz = OneOvSPi * y / (y*y + x*x)
!
  END FUNCTION rlorentz
!
!---------------------------------------------------------------------------
  COMPLEX(rp) ELEMENTAL FUNCTION clorentz(x,y)
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  REAL(rp) :: denom
!
  denom = OneOvSPi / (x*x + y*y)
  clorentz = CMPLX(y*denom,x*denom,KIND=rp)
!
  END FUNCTION clorentz
!
!---------------------------------------------------------------------------
! Real Voigt region IV 2pt GL integration
!
  REAL(rp) ELEMENTAL FUNCTION rvoigth2(x,y)
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  REAL(rp), PARAMETER :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
  REAL(rp), PARAMETER :: gw = 0.28209479177388_rp
  REAL(rp) :: y2
!
  y2 = y*y
  rvoigth2 = gw * y * (1.0_rp/(y2 + (x-gx)**2) + &
           &           1.0_rp/(y2 + (x+gx)**2))
!
  END FUNCTION rvoigth2
!
!---------------------------------------------------------------------------
! Voigt region  2pt GL integration
!
  COMPLEX(rp) ELEMENTAL FUNCTION cvoigth2(x,y)
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu) / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  REAL(rp), PARAMETER :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
  REAL(rp), PARAMETER :: gw = 0.28209479177388_rp
  REAL(rp) :: denom1,denom2,xm,xp,y2
!
  xm = x-gx
  xp = x+gx
  y2 = y**2
  denom1 = gw/(y2 + xm**2)
  denom2 = gw/(y2 + xp**2)
  cvoigth2 = CMPLX(y*(denom1+denom2),xm*denom1 + xp*denom2,KIND=rp)
!
  END FUNCTION cvoigth2
!
!---------------------------------------------------------------------------
! Real Voigt region IV 6pt GL integration
!
  REAL(rp) ELEMENTAL FUNCTION rvoigth6(x,y)
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  INTEGER, PARAMETER :: n = 3
  REAL(rp), PARAMETER :: Pi = 3.1415926535897932385_rp
  REAL(rp), PARAMETER :: gx(n) = (/ 4.36077411927616508271e-1_rp, &
                                    1.33584907401369694957_rp,    &
                                    2.35060497367449222280_rp /)
  REAL(rp), PARAMETER :: gw(n) = (/ 7.24629595224392524608e-1_rp, &
                                    1.57067320322856644036e-1_rp, &
                                    4.53000990550884564224e-3_rp /) / Pi
  REAL(rp) :: y2
  INTEGER :: i
!
  rvoigth6 = 0.0
  y2 = y**2
  DO i = 1 , n
    rvoigth6 = rvoigth6 + gw(i) * y * (1.0_rp/(y2 + (x-gx(i))**2) + &
                                   &   1.0_rp/(y2 + (x+gx(i))**2))
  ENDDO
!
  END FUNCTION rvoigth6
!
!---------------------------------------------------------------------------
! Voigt region IV 6pt GL integration
!
  COMPLEX(rp) ELEMENTAL FUNCTION cvoigth6(x,y)
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu) / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internals
!
  INTEGER, PARAMETER :: n = 3
  REAL(rp), PARAMETER :: Pi = 3.1415926535897932385_rp
  REAL(rp), PARAMETER :: gx(n) = (/ 4.36077411927616508271e-1_rp, &
                                    1.33584907401369694957_rp,    &
                                    2.35060497367449222280_rp /)
  REAL(rp), PARAMETER :: gw(n) = (/ 7.24629595224392524608e-1_rp, &
                                    1.57067320322856644036e-1_rp, &
                                    4.53000990550884564224e-3_rp /) / Pi
  INTEGER :: i
  REAL(rp) :: denom1,denom2,xp,xm
!
  cvoigth6 = CMPLX(0.0_rp,0.0_rp)
  DO i = 1 , n
    xm = x-gx(i)
    xp = x+gx(i)
    denom1 = gw(i)/(y**2 + xm**2)
    denom2 = gw(i)/(y**2 + xp**2)
    cvoigth6 = cvoigth6+CMPLX(y*(denom1+denom2),xm*denom1+xp*denom2,KIND=rp)
  ENDDO
!
  END FUNCTION cvoigth6
!
!---------------------------------------------------------------------------
  REAL(rp) ELEMENTAL FUNCTION rhui6(x,y)
!
! Voigt region II Hui polynomial
! This too looks complicated to split into reals and imaginaries so I
! will do the complex arithmetic here
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internal stuff
!
  REAL(rp), PARAMETER :: a(7) = (/ &
     &     122.607931777104326_rp, 214.382388694706425_rp, &
     &     181.928533092181549_rp,  93.155580458138441_rp, &
     &      30.180142196210589_rp,   5.912626209773153_rp, &
     &       0.564189583562615_rp /)
!
  REAL(rp), PARAMETER :: b(7) = (/ &
     &     122.607931773875350_rp, 352.730625110963558_rp, &
     &     457.334478783897737_rp, 348.703917719495792_rp, &
     &     170.354001821091472_rp,  53.992906912940207_rp, &
     &      10.479857114260399_rp /)
!
  INTEGER :: i
  REAL :: rs,rt,is,it,r(5),q(5)
!
! note that this dimension is 2 less than the a and b coefficients
! because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
! q(0) = 0.0, q(1) = -x
! fill the r, q coefficients with this recursion
!
  r(1) = y**2 - x**2
  q(1) = -2.0_rp*x*y
  DO i = 2 , 5
    r(i) = x * q(i-1) + y * r(i-1)
    q(i) = y * q(i-1) - x * r(i-1)
  ENDDO
  rs = a(1) + a(2)*y + SUM(a(3:)*r)
  rt = b(1) + b(2)*y + SUM(b(3:)*r) + x * q(5) + y * r(5)
  is =      - a(2)*x + SUM(a(3:)*q)
  it =      - b(2)*x + SUM(b(3:)*q) + y * q(5) - x * r(5)
!
  rhui6 = (rs*rt + is*it) / (rt*rt + it*it)
!
! FYI ihui6 = (is*rt - rs*it) / (rt*rt + it*it)
!
  END FUNCTION rhui6
!
!---------------------------------------------------------------------------
  COMPLEX(rp) ELEMENTAL FUNCTION chui6(x,y)
!
! Voigt region II Hui polynomial
! This too looks complicated to split into reals and imaginaries so I
! will do the complex arithmetic here
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu) / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internal stuff
!
  REAL(rp), PARAMETER :: a(7) = (/ &
     &     122.607931777104326_rp, 214.382388694706425_rp, &
     &     181.928533092181549_rp,  93.155580458138441_rp, &
     &      30.180142196210589_rp,   5.912626209773153_rp, &
     &       0.564189583562615_rp /)
!
  REAL(rp), PARAMETER :: b(7) = (/ &
     &     122.607931773875350_rp, 352.730625110963558_rp, &
     &     457.334478783897737_rp, 348.703917719495792_rp, &
     &     170.354001821091472_rp,  53.992906912940207_rp, &
     &      10.479857114260399_rp /)
!
  INTEGER :: i
  REAL :: rs,rt,is,it,r(5),q(5),denom
!
! note that this dimension is 2 less than the a and b coefficients
! because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
! q(0) = 0.0, q(1) = -x
! fill the r, q coefficients with this recursion
!
  r(1) = y**2 - x**2
  q(1) = -2.0_rp*x*y
  DO i = 2, 5
    r(i) = x * q(i-1) + y * r(i-1)
    q(i) = y * q(i-1) - x * r(i-1)
  ENDDO
!
  rs = a(1) + a(2)*y + SUM(a(3:)*r)
  rt = b(1) + b(2)*y + SUM(b(3:)*r) + x * q(5) + y * r(5)
  is =      - a(2)*x + SUM(a(3:)*q)
  it =      - b(2)*x + SUM(b(3:)*q) + y * q(5) - x * r(5)
  denom = 1.0_rp / (rt*rt + it*it)
  chui6 = CMPLX((rs*rt+is*it)*denom, (is*rt-rs*it)*denom,KIND=rp)
!
  END FUNCTION chui6
!
!---------------------------------------------------------------------------
  REAL(rp) ELEMENTAL FUNCTION rdrayson(x,y)
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu) / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internal stuff
!
  INTEGER :: I, J, MAXJ, N
  REAL(rp), PARAMETER :: TwoOvSPi = 1.1283791670955125739_rp  ! 2.0/Sqrt(Pi)
!
! This is a way to set the mesh points without having to use a SAVE statement
!
  REAL(rp), Parameter :: H = 0.2_rp
  REAL(rp), PARAMETER :: HN(26) = (/(H*(I-1),I=1,26)/)
  REAL(rp), PARAMETER :: RI(15) = (/(-I/2.0_rp,I=1,15)/)
!
  REAL(rp), Parameter :: Dawson(26) = (/  &
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
!
  REAL(rp), PARAMETER :: FDer1(26) = &
           & (/(1.0_rp-2.0_rp*HN(I)*Dawson(I),I=1,26)/)
  REAL(rp), PARAMETER :: FDer2(26) = &
           & (/((HN(I)*FDer1(I)+Dawson(I))/RI(2),I=1,26)/)
  REAL(rp), PARAMETER :: FDer3(26) = &
           & (/((HN(I)*FDer2(I)+FDer1(I))/RI(3),I=1,26)/)
  REAL(rp), PARAMETER :: FDer4(26) = &
           & (/((HN(I)*FDer3(I)+FDer2(I))/RI(4),I=1,26)/)
  REAL(rp), PARAMETER :: FDer5(26) = &
           & (/((HN(I)*FDer4(I)+FDer3(I))/RI(5),I=1,26)/)
  REAL(rp), PARAMETER :: FDer6(26) = &
           & (/((HN(I)*FDer5(I)+FDer4(I))/RI(6),I=1,26)/)
!
! internal stuff
!
  REAL(rp) :: Y2, DX, F, FD, W, dely
  Y2 = y*y
!
!******** Region I. Compute Dawson's function at x from Taylor series
!
  J = X / H
  N = 1 + MIN(J, 25)
  DX = X - HN(N)
  F = 0.0_rp
  IF(X > 0.05*H) &
     F = (((((FDer6(N)*DX + FDer5(N))*DX + FDer4(N))*DX  + &
           &  FDer3(N))*DX + FDer2(N))*DX + FDer1(N))*DX + Dawson(N)
  IF(Y <= 1.0e-12_rp) THEN
    rdrayson = EXP(-X*X)
    RETURN
  ENDIF
!
!  Taylor series expansion about y = 0.0
!
  dely = -Y
  FD = 1.0_rp - 2.0_rp * X * F
  w = FD * dely
  J = 5.0_rp + (12.5_rp - X) * 0.8_rp * Y
  MAXJ = MIN(J, 14)
  DO J = 2, MAXJ,2
    F  = (X*FD + F) / RI(J)
    FD  = (X*F + FD) / RI(J+1)
    dely = -y2 * dely
    w = w + FD * dely
  END DO
!
  rdrayson = EXP(y2-x*x)*COS(2.0_rp*x*y) + TwoOvSpi*w
!
  END FUNCTION rdrayson
!
!---------------------------------------------------------------------------
  COMPLEX(rp) ELEMENTAL FUNCTION cdrayson(x,y)
!
  REAL(rp), INTENT(in) :: x ! sqrt(ln2)*delnu) / wd
  REAL(rp), INTENT(in) :: y ! sqrt(ln2)*wc / wd
!
! Internal stuff
!
  INTEGER :: I, J, MAXJ, N
  REAL(rp), PARAMETER :: TwoOvSPi = 1.1283791670955125739_rp  ! 2.0/Sqrt(Pi)
!
! This is a way to set the mesh points without having to use a SAVE statement
!
  REAL(rp), Parameter :: H = 0.2_rp
  REAL(rp), PARAMETER :: HN(26) = (/(H*(I-1),I=1,26)/)
  REAL(rp), PARAMETER :: RI(15) = (/(-I/2.0_rp,I=1,15)/)
!
  REAL(rp), Parameter :: Dawson(26) = (/  &
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
!
  REAL(rp), PARAMETER :: FDer1(26) = &
           & (/(1.0_rp-2.0_rp*HN(I)*Dawson(I),I=1,26)/)
  REAL(rp), PARAMETER :: FDer2(26) = &
           & (/((HN(I)*FDer1(I)+Dawson(I))/RI(2),I=1,26)/)
  REAL(rp), PARAMETER :: FDer3(26) = &
           & (/((HN(I)*FDer2(I)+FDer1(I))/RI(3),I=1,26)/)
  REAL(rp), PARAMETER :: FDer4(26) = &
           & (/((HN(I)*FDer3(I)+FDer2(I))/RI(4),I=1,26)/)
  REAL(rp), PARAMETER :: FDer5(26) = &
           & (/((HN(I)*FDer4(I)+FDer3(I))/RI(5),I=1,26)/)
  REAL(rp), PARAMETER :: FDer6(26) = &
           & (/((HN(I)*FDer5(I)+FDer4(I))/RI(6),I=1,26)/)
!
! internal stuff
!
  REAL(rp) :: DX, F, FD,wr,wi,dely,twoxy
!
!******** Region I. Compute Dawson's function at x from Taylor series
!
  J = X / H
  N = 1 + MIN(J, 25)
  DX = X - HN(N)
  F = 0.0_rp
  IF(X.gt.0.05*H) &
     F = (((((FDer6(N)*DX + FDer5(N))*DX + FDer4(N))*DX  + &
           &  FDer3(N))*DX + FDer2(N))*DX + FDer1(N))*DX + Dawson(N)
  IF(Y <= 1.0e-12_rp) THEN
    cdrayson = CMPLX(EXP(-X*X),TwoOvSPi * F,KIND=rp)
    RETURN
  ENDIF
!
!  Taylor series expansion about y = 0.0
!
  wr = EXP(y*y-x*x)
  twoxy = 2.0_rp*x*y
  cdrayson = CMPLX(wr*COS(twoxy),-wr*SIN(twoxy),KIND=rp)
  dely = -TwoOvSpi*Y
  FD = 1.0_rp - 2.0_rp * X * F
  wi = TwoOvSpi*F
  wr = FD*dely
  j = 5.0 + 10.0*y - 0.4*twoxy
  MAXJ = MIN(J, 14)
  DO J = 2, MAXJ,2
    F  = (X*FD + F) / RI(J)
    FD = (X*F + FD) / RI(J+1)
    dely = y * dely
    wi = wi + f*dely
    dely = -y*dely
    wr = wr + FD*dely
  END DO
!
  cdrayson = cdrayson + CMPLX(wr,wi,KIND=rp)
!
  END FUNCTION cdrayson
!
!---------------------------------------------------------------------
! This function will compute a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
!
! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
!
 Subroutine Slabs_prep(t,m,v0,el,w,ps,p,n,ns,i,q,delta,gamma,n1,n2, &
                    &  v0s,x1,y,yi,slabs1,dslabs1)
!
  implicit none
!
! inputs:
!
  Real(rp), INTENT(IN) :: t        ! Temperature K
  Real(rp), INTENT(IN) :: m        ! Molecular mass amu
  Real(r8), INTENT(IN) :: v0       ! Line center frequency MHz
  Real(rp), INTENT(IN) :: el       ! Lower state energy cm-1
  Real(rp), INTENT(IN) :: w        ! Collision broadening parameter
                                   ! MHz/mbar at 300 K
  Real(rp), INTENT(IN) :: ps       ! Pressure shift parameter in MHz/mbar
  Real(rp), INTENT(IN) :: p        ! Pressure mbar
  Real(rp), INTENT(IN) :: n        ! Temperature power dependence of w
  Real(rp), INTENT(IN) :: ns       ! Temperature power dependence of ps
  Real(rp), INTENT(IN) :: i        ! Integrated spectral intensity
                                   ! Log(nm**2 MHz) at 300 K
  Real(rp), INTENT(IN) :: q(3)     ! Logarithm of the partition function
                                   ! At 300 , 225 , and 150 K
  Real(rp), INTENT(IN) :: delta    ! Delta interference coefficient at 300K 1/mb
  Real(rp), INTENT(IN) :: gamma    ! Gamma               "
  Real(rp), INTENT(IN) :: n1       ! Temperature dependency of delta
  Real(rp), INTENT(IN) :: n2       ! Temperature dependency of gamma
!
! outputs:
!
  Real(r8), INTENT(OUT) :: v0s     ! Pressure shifted line position
  Real(rp), INTENT(OUT) :: x1      ! Sqrt(Ln(2))/Doppler half width MHz
  Real(rp), INTENT(OUT) :: y       ! Sqrt(Ln(2))*collision width /
                                   !             doppler width
  Real(rp), INTENT(OUT) :: yi      ! Interference contribution
  Real(rp), INTENT(OUT) :: slabs1  ! Frequency independent piece of slabs
  Real(rp), INTENT(OUT) :: dslabs1 ! Derivative of slabs1 w.r.t. v0
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
  Real(rp), Parameter :: i2abs = 3.402136078e9_rp
  Real(rp), Parameter :: dc = 3.58117369e-7_rp
  Real(rp), Parameter :: boltzcm = 0.6950387_rp
  Real(rp), Parameter :: boltzmhz = 20836.74_rp
  Real(rp), Parameter :: sqrtln2 = 8.32554611157698e-1_rp
  Real(rp), Parameter :: loge = 4.34294481903251828e-1_rp
  Real(rp), Parameter :: oned300 = 1.0_rp/300.0_rp
!
  Real(rp), parameter :: tl1 = 1.76091259055681e-1_rp     ! Log10(225/150)
  Real(rp), parameter :: tl2 = 1.24938736608300e-1_rp     ! Log10(300/225)
!
! Internal data:
!
  Real(rp) :: Wd, Q_Log, betae, betav, t3t, onedt, r, e1, e2, de1, &
              de2, g, s, ds
!
! The action begins here
!
  onedt = 1.0_rp / t
  t3t = 300.0_rp * onedt
  yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))
!
  if (t < 225.0_rp) then
    r = (q(2)-q(3))/tl1
    Q_Log = q(2)-q(1)+r*Log10(t/225.0_rp)
  else
    r = (q(1)-q(2))/tl2
    Q_Log = r*Log10(t/300.0_rp)
  endif
!
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
  g = 1.0_rp - e2
  s = (1.0_rp - e1) / g
  ds = (-de1*g+(1.0_rp-e1)*de2)/(g*g)
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
 Subroutine Slabs_prep_wder(t,m,v0,el,w,ps,p,n,ns,i,q,delta,gamma,n1,n2, &
     &      v0s,x1,y,yi,slabs1,dx1_dv0,dy_dv0,dslabs1_dv0)
!
! inputs:
!
  Real(rp), INTENT(IN) :: t       ! Temperature K
  Real(rp), INTENT(IN) :: m       ! Molecular mass amu
  Real(r8), INTENT(IN) :: v0      ! Line center frequency MHz
  Real(rp), INTENT(IN) :: el      ! Lower state energy cm-1
  Real(rp), INTENT(IN) :: w       ! Collision broadening parameter
                                  ! MHz/mbar at 300 K
  Real(rp), INTENT(IN) :: ps      ! Pressure shift parameter in MHz/mbar
  Real(rp), INTENT(IN) :: p       ! Pressure mbar
  Real(rp), INTENT(IN) :: n       ! Temperature power dependence of w
  Real(rp), INTENT(IN) :: ns       ! Temperature power dependence of ps
  Real(rp), INTENT(IN) :: i       ! Integrated spectral intensity
                                  ! Log(nm**2 MHz) at 300 K
  Real(rp), INTENT(IN) :: q(3)    ! Logarithm of the partition function
                                  ! At 300 , 225 , and 150 K
  Real(rp), INTENT(IN) :: delta   ! Delta interference coefficient at 300K 1/mb
  Real(rp), INTENT(IN) :: gamma   ! Gamma               "
  Real(rp), INTENT(IN) :: n1      ! Temperature dependency of delta
  Real(rp), INTENT(IN) :: n2      ! Temperature dependency of gamma
!
! outputs:
!
  Real(r8), INTENT(OUT) :: v0s       ! Pressure shifted line position
  Real(rp), INTENT(OUT) :: x1        ! Sqrt(Ln(2))/Doppler half width MHz
  Real(rp), INTENT(OUT) :: y         ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
  Real(rp), INTENT(OUT) :: yi        ! Interference contribution
  Real(rp), INTENT(OUT) :: slabs1    ! Frequency independent piece of slabs
!
  Real(rp), INTENT(OUT) :: dx1_dv0       ! Derivative of x1 w.r.t. v0
  Real(rp), INTENT(OUT) :: dy_dv0        ! Derivative of y w.r.t. v0
  Real(rp), INTENT(OUT) :: dslabs1_dv0   ! Derivative of slabs1 w.r.t. v0
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
  Real(r8) Wd, Q_Log, betae, betav, t3t, onedt, r, e1, e2, &
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
! --------------------------------  slabs_prep_arrays   -----
!
 Subroutine Slabs_Prep_Arrays(Spectag,nl,t,p,mass,Qlog,Catalog,v0s,x1,y, &
                           &  yi,slabs1,dslabs1_dv0)
!
 type(catalog_T) :: Catalog

 Integer(ip), intent(in) :: Spectag, nl

 Real(rp), intent(in) :: t, p, mass,Qlog(:)

 Real(r8), intent(out) :: v0s(:)
 Real(rp), intent(out) :: x1(:),y(:),yi(:),slabs1(:),dslabs1_dv0(:)
!
  Integer :: j, k

  if(Spectag==18999 .or. Spectag==28964 .or. Spectag==28965) Return
!
! Check for anything but liquid water and dry air:
!
  do j = 1, nl
!
! Prepare the temperature weighted coefficients:
!
    k = Catalog%Lines(j)
    Call Slabs_prep(t,mass,Lines(k)%V0,Lines(k)%EL,Lines(k)%W,        &
      &  Lines(k)%PS, p, Lines(k)%N,Lines(k)%NS,Lines(k)%STR,         &
      &  Qlog,Lines(k)%DELTA,Lines(k)%GAMMA,Lines(k)%N1,Lines(k)%N2,  &
      &  v0s(j),x1(j),y(j),yi(j),slabs1(j),dslabs1_dv0(j))
!
  end do
!
 End Subroutine Slabs_Prep_Arrays
!
!------------------------------------------------------------
  SUBROUTINE get_gl_slabs_arrays(Catalog,p_path,t_path,vel_z,gl_slabs, &
                             &   no_ele,dt)
!
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
!
  Type(Catalog_T), dimension(:), intent(in) :: Catalog

 Integer(ip), INTENT(IN) :: no_ele
!
 REAL(rp), INTENT(IN) :: p_path(:) ! Pressure in hPa or mbar
 REAL(rp), INTENT(IN) :: t_path(:)
!
 Real(rp), INTENT(IN) :: vel_z, dt
!
 Type (slabs_struct), POINTER :: gl_slabs(:,:)
!
!  ----------------
!  Local variables:
!  ----------------
!
 Real(rp), PARAMETER :: c = 299792.4583_rp     ! Speed of Light Km./Sec.
!
 Integer :: nl,i,j,no_sps,spectag
!
 Real(rp) :: mass, Vel_z_correction, Qlog(3)
!
! Begin code:
!
   no_sps = Size(Catalog)
!
   Vel_z_correction = 1.0_rp - vel_z / c
!
   DO i = 1, no_sps
!
     Spectag = Catalog(i)%spec_tag
     mass = Real(Spectag) / 1000.0_rp
!
     nl = Size(Catalog(i)%Lines)
     gl_slabs(1:no_ele,i)%no_lines = nl

     Qlog(1:3) = Catalog(i)%QLOG(1:3)
!
     do j = 1, no_ele
!
       Call Slabs_Prep_Arrays(Spectag,nl,t_path(j)+dt,p_path(j),mass,Qlog, &
         &  Catalog(i),gl_slabs(j,i)%v0s,gl_slabs(j,i)%x1,gl_slabs(j,i)%y, &
         &  gl_slabs(j,i)%yi,gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dslabs1_dv0)
!
!  Apply velocity corrections:
!
       gl_slabs(j,i)%v0s = gl_slabs(j,i)%v0s * Vel_z_correction
!
     end do
!
   END DO              ! On i
!
   Return
!
 END SUBROUTINE get_gl_slabs_arrays
!
!=====================================================================

End module SLABS_SW_M
! $Log$
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
