module SLABS_SW_M

  use MLSCommon, only: R8, RP, IP

  use SpectroscopyCatalog_m, only: CATALOG_T, Lines

  implicit NONE

  private
  public :: DVoigt_Spectral, Slabs, Slabswint, Voigt_Lorentz, &
        &  Real_Simple_Voigt, Simple_Voigt, RLorentz, CLorentz, RVoigth2, &
        &  CVoigth2, RVoigth6, CVoigth6, RHui6, CHui6, RDrayson, &
        &  CDrayson, Slabs_Prep, Slabs_Prep_Arrays, Get_GL_Slabs_Arrays

  real(rp), parameter :: OneOvSPi = 0.56418958354775628695_rp  ! 1.0/Sqrt(Pi)

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

  ! --------------------------------------------  dVoigt_spectral  -----
  SUBROUTINE dVoigt_spectral ( dNu, Nu0, x1, yi, y, w, t, tanh1, slabs1, SwI, &
                         &  dslabs1_dNu0, dSwI_dw, dSwI_dn, dSwI_dNu0 )

! Compute the Voigt function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep_wder()
!       routine to compute dslabs1_dNu0

! NOTE: In here and in all other routines in this module, 
!       tanh1 = tanh(nu*expa / 2.0)

    REAL(r8), INTENT(in) :: nu0
    REAL(rp), INTENT(in) :: dnu, x1, yi, y, w, t, tanh1, slabs1             
    real(rp), intent(in), optional :: dslabs1_dNu0                        

    real(rp), intent(out) :: SwI, dSwI_dw,dSwI_dn,dSwI_dNu0               

    real(rp) :: x, u, v, du_dx, du_dy, dv_dx, dv_dy, q, b, g, z, r    

!    real(rp) :: dq_dv0
    real(rp) :: dx_dv0, du_dv0, dv_dv0, db_dv0, dg_dv0, dz_dv0, & 
                dr_dv0, dvvw_dv0, vvw, slabs2

    x = x1 * dNu                                                          
    call simple_voigt(x,y,u,v)  

!  Van Vleck - Wieskopf (VVW) line shape with Voigt

    q = 1.0_rp + dNu / Nu0

    b = x1 * (2.0_r8 * Nu0 + dNu)
    g = b * b + y * y
    z = (y - b * yi) / g
    r = z * OneOvSPi + yi * v
    vvw = (u + r) * q
    slabs2 = slabs1 * tanh1
    SwI = slabs2 * vvw

    du_dx = 2.0_rp * (y * v - x * u)
    du_dy = 2.0_rp * (y * u + x * v - OneOvSPi)

    dv_dx = -du_dy         ! Cauchy-Riemann equation
    dv_dy =  du_dx         ! Cauchy-Riemann equation

! Compute the derivative of SwI w.r.t. w

    dSwI_dw = q * slabs2* (y/w) * (du_dy + yi*du_dx + &
                                &   OneOvSPi*((1.0_rp-2.0_rp*z*y)/g))

! Compute the derivative of SwI w.r.t. n

    dSwI_dn = q * slabs2 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)

! Finaly, compute the derivative of SwI w.r.t. Nu0

! ***** Analytically *****

!    dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
    dx_dv0 = -x1
    du_dv0 = du_dx * dx_dv0
    dv_dv0 = dv_dx * dx_dv0
    db_dv0 = x1
    dg_dv0 = 2.0_rp * b*db_dv0
    dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
    dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
!    dvvw_dv0 = (du_dv0+dr_dv0)*q + dq_dv0*(u+r)
    dvvw_dv0 = (du_dv0+dr_dv0)*q
    if ( present(dslabs1_dNu0) ) dSwI_dNu0 = dslabs1_dNu0*vvw + slabs2*dvvw_dv0

  end subroutine dVoigt_spectral

  ! ------------------------------------------------------  Slabs  -----
  REAL(rp) FUNCTION Slabs ( dNu, v0s, x1, tanh1, slabs1, y )

    REAL(r8), INTENT(in) :: v0s
    REAL(rp), INTENT(in) :: dNu, x1, tanh1, slabs1, y

!  Note: dNu = v - v0s

! If the molecular transition and temperature have not changed but
! frequency has enter here.

! inputs: dNu , x1 , tanh1, slabs1 , y, v0s, yi
! output: slabs

    real(rp) :: u

    call real_simple_voigt(x1*dNu,y,u)

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

    Slabs = slabs1 * (1.0_rp + dNu / v0s) * tanh1 &
              * (u + OneOvSPi*y/((x1*(2.0_rp*v0s+dNu))**2 + y*y))

  end function Slabs

  ! --------------------------------------------------  Slabswint  -----
  REAL(rp) FUNCTION Slabswint ( dNu, v0s, x1, tanh1, slabs1, y, yi )

    REAL(r8), INTENT(in) :: v0s
    REAL(rp), INTENT(in) :: dNu, x1, tanh1, slabs1, y, yi

!  Note: dNu = v - v0s

! If the molecular transition and temperature have not changed but
! frequency has enter here.

! inputs: dNu , x1 , slabs1 , y, v0s, yi
! output: slabswint (slab with interference)

    real(rp) :: x, u, q, p, z, y2, w

    x = x1 * dNu
    call real_simple_voigt(x,y,u)

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

    q = (1.0_rp + dNu / v0s)
    p = x1 * (2.0_rp * v0s + dNu)
    y2 = y*y
    z = OneOvSPi*((y - p*yi)/(p*p + y2) + yi*x/(x*x+y2))
    w = (u + z) * q
    Slabswint = slabs1 * tanh1 *  w

  end function Slabswint

  ! ----------------------------------------------  Voigt_Lorentz  -----

  subroutine Voigt_Lorentz ( dNu,  Nu0,  x1,  yi,  y,  w,  t,  tanh1, slabs1,  &
                         & VL, dslabs1_dNu0,  dVL_dw,  dVL_dn,  dVL_dNu0 )

! Compute the Voigt/Lorentz function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

    real(r8), intent(in) :: nu0
    real(rp), intent(in) :: dNu, x1, yi, y, w, t, tanh1, slabs1, dslabs1_dNu0

    real(rp), intent(out) :: VL, dVL_dw, dVL_dn, dVL_dNu0

    real(rp) :: xj, zj, q, y2, u, v, up1, up2, dn1, dn2, dup1, &
     &          dup2, ddn1, ddn2, dy_dw, dy_dn, dSum_dw, dSum_dn, &
     &          dSum_dNu0, slabs2
!    real(rp) :: dq_dNu0, Sum

    q = 1.0_rp + dNu / Nu0

    y2 = y * y
    xj = x1 * dNu
    zj = x1 * (2.0_r8 * Nu0 + dNu)
    dn1 = zj * zj + y2
    up1 = y - zj * yi

! Van Vleck - Wieskopf (VVW) line shape with Voigt

    call simple_voigt ( xj, y, u, v )
    dup1 = up1 * OneOvSPi / dn1 + yi * v
    ddn2 = u + dup1
!    Sum = ddn2 / OneOvSPi
    slabs2 = slabs1 * tanh1
    VL = slabs2 * ddn2 * q            ! This is the Voigt + VVW correction

    dn2 = xj * xj + y2
    up2 = y - yi * xj

    dy_dw = y / w
    dup1 = dy_dw
    ddn1 = 2.0 * y * dy_dw
    dup2 = dy_dw
    ddn2 = 2.0 * y * dy_dw

    dSum_dw = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
   &          (dn2*dup2-up2*ddn2)/(dn2*dn2)

    dVL_dw = OneOvSPi * slabs2 * q * dSum_dw

    dy_dn = y * Log(300.0/t)
    dup1 = dy_dn
    ddn1 = 2.0 * y * dy_dn
    dup2 = dy_dn
    ddn2 = 2.0 * y * dy_dn
    dSum_dn = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
   &          (dn2*dup2-up2*ddn2)/(dn2*dn2)

    dVL_dn = OneOvSPi * slabs2 * q * dSum_dn

    dup2 =  yi * x1               !  x1 = -dxj_dNu0
    ddn2 = -2.0 * xj * x1         !  x1 = -dxj_dNu0
    dSum_dNu0 = (dn2*dup2-up2*ddn2)/(dn2*dn2)
!    dq_dNu0 = -(dNu+Nu0)/(Nu0*Nu0)

!    dVL_dNu0 = OneOvSPi * (dslabs1_dNu0 * q * Sum          + &
!              &         slabs2 * dq_dNu0 * Sum + &
!              &         slabs2 * q * dSum_dNu0)

    dVL_dNu0 = OneOvSPi * slabs2 * q * dSum_dNu0

  end subroutine Voigt_Lorentz

  ! ------------------------------------------  Real_Simple_Voigt  -----
  elemental subroutine Real_Simple_Voigt ( x, y, u )

! simple REAL(Voigt) function

! inputs

    real(rp), intent(in) :: x ! doppler width, del frequency ratio
    real(rp), intent(in) :: y ! doppler width, collision width ratio

! outputs

    real(rp), intent(out) :: u ! real part of Voigt

! internals

    real(rp), parameter :: xl=5.2_rp, yl=0.05_rp, yh=0.6_rp, dydx=(yh-yl)/xl
    real(rp) :: xa

! This is sorted in likely occurance of each case

    xa = ABS(x)

! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed

    if ( y + 0.666666*xa > 100.0_rp ) then

      u = rlorentz(xa,y)

    else if ( y + 0.6875_rp * xa > 11.0_rp ) then

! Drayson's quick 2pt hermite integral (essentially a lorentz)

      u = rvoigth2(xa,y)

    else if ( y > 0.6_rp .OR. y > yl + dydx*xa ) then

! Intermediate region

      u = rhui6(xa,y)

    else if ( xa > 5.2_rp ) then

! small y large x limit

      u = rvoigth6(xa,y)

    else

! Near line center where Doppler dominates Pressure broadening

      u = rdrayson(xa,y)

    end if

  end subroutine Real_Simple_Voigt

  ! -----------------------------------------------  Simple_Voigt  -----
  elemental subroutine Simple_Voigt ( x, y, u, v )

! simple Voigt function

! inputs

    real(rp), intent(in) :: x ! doppler width, del frequency ratio
    real(rp), intent(in) :: y ! doppler width, collision width ratio

! outputs

    real(rp), intent(out) :: u ! real part of Voigt
    real(rp), optional, intent(out) :: v ! imaginary part of Voigt

! internals

    real(rp), parameter :: xl=5.2_rp, yl=0.05_rp, yh=0.6_rp, dydx=(yh-yl)/xl
    real(rp) :: xa
    complex(rp) :: uv

! This is sorted in likely occurance of each case

    xa = ABS(x)

! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed

!    if ( y + 0.666666*xa > 100.0_rp ) then
! NOTE: clorentz is not accurate enough for spectral derivative
!       computations. This may be something for Van S. to investigate
!       later.

!      uv = clorentz(xa,y)

!    else if ( y + 0.6875_rp * xa > 11.0_rp ) then

    if ( y + 0.6875_rp * xa > 11.0_rp ) then

! Drayson's quick 2pt hermite integral (essentially a lorentz)

      uv = cvoigth2(xa,y)

    else if ( y > 0.6_rp .OR. y > yl + dydx*xa ) then

! Intermediate region

      uv = chui6(xa,y)

    else if ( xa > 5.2_rp ) then

! small y large x limit

      uv = cvoigth6(xa,y)

    else

! Near line center where Doppler dominates Pressure broadening

      uv = cdrayson(xa,y)

    end if

    u = REAL(uv,KIND=rp)
    if ( present(v) ) v = SIGN(AIMAG(uv),x)

  end subroutine Simple_Voigt

  ! ---------------------------------------------------  RLorentz  -----
  real(rp) elemental function RLorentz ( x, y )

! Real Lorentz

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

! Internals

    rlorentz = OneOvSPi * y / (y*y + x*x)

  end function RLorentz

  ! ---------------------------------------------------  CLorentz  -----
  complex(rp) elemental function CLorentz(x,y)

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    real(rp) :: denom

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
    real(rp), parameter :: gw = 0.28209479177388_rp
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
    real(rp), parameter :: gw = 0.28209479177388_rp
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

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    integer, parameter :: n = 3
    real(rp), parameter :: Pi = 3.1415926535897932385_rp
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

    real(rp), intent(in) :: x ! sqrt(ln2)*delnu) / wd
    real(rp), intent(in) :: y ! sqrt(ln2)*wc / wd

!   Internals

    integer, parameter :: n = 3
    real(rp), parameter :: Pi = 3.1415926535897932385_rp
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
! This looks too complicated to split into reals and imaginaries so I
! will do the complex arithmetic here

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

!   note that this dimension is 2 less than the a and b coefficients
!   because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
!   q(0) = 0.0, q(1) = -x
!   fill the r, q coefficients with this recursion

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
! This too looks complicated to split into reals and imaginaries so I
! will do the complex arithmetic here

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

! note that this dimension is 2 less than the a and b coefficients
! because we start the r and q coeffients at 2. r(0)=1.0,r(1) = y
! q(0) = 0.0, q(1) = -x
! fill the r, q coefficients with this recursion

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
    real(rp), parameter :: TwoOvSPi = 1.1283791670955125739_rp  ! 2.0/Sqrt(Pi)

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
    y2 = y*y

!******** Region I. Compute Dawson's function at x from Taylor series

    j = x / h
    n = 1 + min(j, 25)
    dx = x - hn(n)
    f = 0.0_rp
    if ( x > 0.05*h) &
       f = (((((fDer6(n)*dx + fDer5(n))*dx + fDer4(n))*dx  + &
             &  fDer3(n))*dx + fDer2(n))*dx + fDer1(n))*dx + Dawson(n)
    if ( y <= 1.0e-12_rp ) then
      rdrayson = exp(-X*X)
      return
    end if

!  Taylor series expansion about y = 0.0

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
    real(rp), parameter :: TwoOvSPi = 1.1283791670955125739_rp  ! 2.0/Sqrt(Pi)

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
    if ( x.gt.0.05*h) &
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

  ! -------------------------------------------------  Slabs_prep  -----
  Subroutine Slabs_prep ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                      &   n1, n2, &
                      &   v0s, x1, y, yi, slabs1, dslabs1 )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects

! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(rp), intent(in) :: M        ! Molecular mass amu
    real(r8), intent(in) :: V0       ! Line center frequency MHz
    real(rp), intent(in) :: El       ! Lower state energy cm-1
    real(rp), intent(in) :: W        ! Collision broadening parameter
                                     ! MHz/mbar at 300 K
    real(rp), intent(in) :: Ps       ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P        ! Pressure mbar
    real(rp), intent(in) :: N        ! Temperature power dependence of w
    real(rp), intent(in) :: Ns       ! Temperature power dependence of ps
    real(rp), intent(in) :: I        ! Integrated spectral intensity
                                     ! Log(nm**2 MHz) at 300 K
    real(rp), intent(in) :: Q(3)     ! Logarithm of the partition function
                                     ! At 300 , 225 , and 150 K
    real(rp), intent(in) :: Delta    ! Delta interference coefficient at 300K 1/mb
    real(rp), intent(in) :: Gamma    ! Gamma               "
    real(rp), intent(in) :: N1       ! Temperature dependency of delta
    real(rp), intent(in) :: N2       ! Temperature dependency of gamma

! outputs:

    real(r8), intent(out) :: V0s     ! Pressure shifted line position
    real(rp), intent(out) :: X1      ! Sqrt(Ln(2))/Doppler half width MHz
    real(rp), intent(out) :: Y       ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
    real(rp), intent(out) :: Yi      ! Interference contribution
    real(rp), intent(out) :: Slabs1  ! Frequency independent piece of slabs
    real(rp), intent(out) :: Dslabs1 ! Derivative of slabs1 w.r.t. v0

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K
!  sqrtln2  - sqrt(ln(2))
!  loge     - log10(e)

    real(rp), parameter :: I2abs = 3.402136078e9_rp
    real(rp), parameter :: Dc = 3.58117369e-7_rp
    real(rp), parameter :: Boltzcm = 0.6950387_rp
    real(rp), parameter :: Boltzmhz = 20836.74_rp
    real(rp), parameter :: Sqrtln2 = 8.32554611157698e-1_rp
    real(rp), parameter :: Loge = 4.34294481903251828e-1_rp
    real(rp), parameter :: Oned300 = 1.0_rp/300.0_rp

    real(rp), parameter :: LT2 = 2.35218251811136_rp      ! Log10(225)
    real(rp), parameter :: LT3 = 2.47712125471966_rp      ! Log10(300)

    real(rp), parameter :: Tl1 = 0.176091259055681_rp     ! Log10(225/150)
    real(rp), parameter :: Tl2 = 0.124938736608300_rp     ! Log10(300/225)

! Internal data:

    real(rp) :: Wd, Q_Log, betae, betav, t3t, onedt

! The action begins here

   onedt = 1.0_rp / t
   t3t = 300.0_rp * onedt
   yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))

   if ( t < 225.0_rp ) then
     Q_Log = q(2)-q(1)+(q(2)-q(3))/tl1*(Log10(t)-lt2)
   else
     Q_Log =           (q(1)-q(2))/tl2*(Log10(t)-lt3)
   end if

   v0s = v0 + ps * p * (t3t**ns)

   betae = el / boltzcm
   betav = v0s / boltzmhz
   Wd = v0s * Sqrt(t/m) * dc
   x1 = sqrtln2 / Wd
   y = x1 * w * p * (t3t**n)
   slabs1 = i2abs * p * 10.0**(i - Q_Log + loge *  betae * (oned300 - onedt)) &
        & * (1.0_rp + EXP(-betav*onedt)) &
        & / (t * Wd * (1.0_rp - EXP(-betav*oned300)))
   dslabs1 = 0.0_rp

 end subroutine Slabs_prep

  !  -------------------------------------------  Slabs_prep_wder  -----
  Subroutine Slabs_prep_wder ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, &
                            &  gamma, n1, n2,  &
                            &  v0s, x1, y, yi, slabs1, &
                            &  dx1_dv0, dy_dv0, dslabs1_dv0 )

! Slabs_prep_wder: ** ORIGINALLY: Subroutine Slabs_prep()

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! ** CHANGED: Jan/5/00  To Include derivatives of x1,y & slabs w.r.t. Nu0
! ** This routine is obselete and probably can be deleted.

! inputs:

    real(rp), intent(in) :: T       ! Temperature K
    real(rp), intent(in) :: M       ! Molecular mass amu
    real(r8), intent(in) :: V0      ! Line center frequency MHz
    real(rp), intent(in) :: El      ! Lower state energy cm-1
    real(rp), intent(in) :: W       ! Collision broadening parameter
                                    ! MHz/mbar at 300 K
    real(rp), intent(in) :: Ps      ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P       ! Pressure mbar
    real(rp), intent(in) :: N       ! Temperature power dependence of w
    real(rp), intent(in) :: Ns       ! Temperature power dependence of ps
    real(rp), intent(in) :: I       ! Integrated spectral intensity
                                    ! Log(nm**2 MHz) at 300 K
    real(rp), intent(in) :: Q(3)    ! Logarithm of the partition function
                                    ! At 300 , 225 , and 150 K
    real(rp), intent(in) :: Delta   ! Delta interference coefficient at 300K 1/mb
    real(rp), intent(in) :: Gamma   ! Gamma               "
    real(rp), intent(in) :: N1      ! Temperature dependency of delta
    real(rp), intent(in) :: N2      ! Temperature dependency of gamma

! outputs:

    real(r8), intent(out) :: V0s       ! Pressure shifted line position
    real(rp), intent(out) :: X1        ! Sqrt(Ln(2))/Doppler half width MHz
    real(rp), intent(out) :: Y         ! Sqrt(Ln(2))*collision width /
                                       !             doppler width
    real(rp), intent(out) :: Yi        ! Interference contribution
    real(rp), intent(out) :: Slabs1    ! Frequency independent piece of slabs

    real(rp), intent(out) :: Dx1_dv0       ! Derivative of x1 w.r.t. v0
    real(rp), intent(out) :: Dy_dv0        ! Derivative of y w.r.t. v0
    real(rp), intent(out) :: Dslabs1_dv0   ! Derivative of slabs1 w.r.t. v0

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K
!  sqrtln2  - sqrt(ln(2))
!  loge     - log10(e)

    real(r8), parameter :: I2abs = 3.402136078e9_r8
    real(r8), parameter :: Dc = 3.58117369e-7_r8
    real(r8), parameter :: Boltzcm = 0.6950387_r8
    real(r8), parameter :: Boltzmhz = 20836.74_r8
    real(r8), parameter :: Sqrtln2 = 8.32554611157698e-1_r8
    real(r8), parameter :: Loge = 4.34294481903251828e-1_r8
    real(r8), parameter :: Oned300 = 1.0_r8/300.0_r8

    real(r8), parameter :: Tl1 = 1.76091259055681e-1_r8     ! Log10(225/150)
    real(r8), parameter :: Tl2 = 1.24938736608300e-1_r8     ! Log10(300/225)

! Internal data:

    real(r8) :: Wd, Q_Log, betae, betav, t3t, onedt, r

! The action begins here

    onedt = 1.0_r8 / t
    t3t = 300.0_r8 * onedt
    yi = p * (delta*(t3t**n1) + gamma*(t3t**n2))

    if ( t < 225.0_r8 ) then
      r = (q(2)-q(3))/tl1
      Q_Log = q(2)-q(1)+r*Log10(t/225.0_r8)
    else
      r = (q(1)-q(2))/tl2
      Q_Log = r*Log10(t/300.0_r8)
    end if

    v0s = v0 + ps * p * (t3t**ns)

    betae = el / boltzcm
    betav = v0s / boltzmhz
    Wd = v0s * Dsqrt(t/m) * dc
    x1 = sqrtln2 / Wd
    y = x1 * w * p * (t3t**n)
    slabs1 = i2abs * p * 10.0**(i - Q_Log + loge *  betae * (oned300-onedt)) &
        &  * (1.0_rp + EXP(-betav*onedt)) &
        &  / (t * Wd * (1.0_rp - EXP(-betav*oned300)))
    dx1_dv0 = 0.0_rp
    dy_dv0 = 0.0_rp
    dslabs1_dv0 = 0.0_rp

 end subroutine Slabs_prep_wder

  ! -----------------------------------------  Slabs_Prep_Arrays   -----
  Subroutine Slabs_Prep_Arrays ( molecule, nl, t, p, mass, Qlog, Catalog, &
                               & v0s, x1, y, yi, slabs1, dslabs1_dv0 )

    use Molecules, only: L_Extinction

    type(catalog_T) :: Catalog

    integer(ip), intent(in) :: molecule, nl

    real(rp), intent(in) :: t, p, mass,Qlog(:)

    real(r8), intent(out) :: v0s(:)
    REAL(rp), INTENT(out) :: x1(:),y(:),yi(:),slabs1(:),dslabs1_dv0(:)

    integer :: j, k

    if ( any ( molecule == (/ l_extinction /) ) ) return

! Check for anything but dry air or extinction:

    do j = 1, nl

! Prepare the temperature weighted coefficients:

      k = Catalog%Lines(j)
      call Slabs_prep ( t, mass, Lines(k)%V0, Lines(k)%EL, Lines(k)%W,      &
        &  Lines(k)%PS, p, Lines(k)%N, Lines(k)%NS, Lines(k)%STR,           &
        &  Qlog, Lines(k)%DELTA, Lines(k)%GAMMA, Lines(k)%N1, Lines(k)%N2,  &
        &  v0s(j), x1(j), y(j), yi(j), slabs1(j), dslabs1_dv0(j) )

    end do

  end subroutine Slabs_Prep_Arrays

  ! ----------------------------------------  Get_GL_Slabs_Arrays  -----
  subroutine Get_GL_Slabs_Arrays ( Catalog, P_path, T_path, Vel_z, GL_slabs, &
                             &     No_ele, Do_1D )

    use Units, only: SpeedOfLight
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT

    type(Catalog_T), dimension(:), intent(in) :: Catalog

    integer(ip), intent(in) :: no_ele

    real(rp), intent(in) :: p_path(:) ! Pressure in hPa or mbar
    real(rp), intent(in) :: t_path(:)

    real(rp), intent(in) :: vel_z

    type (slabs_struct), pointer :: gl_slabs(:,:)

!  ----------------
!  Local variables:
!  ----------------

    real(rp), parameter :: c = speedOfLight/1000.0_rp ! Speed of Light Km./Sec.

    integer :: nl,i,j,n_sps, k

    real(rp) :: vel_z_correction, Qlog(3)

    Logical :: Do_1D

! Begin code:

    n_sps = Size(Catalog)

    Vel_z_correction = 1.0_rp - vel_z / c

    do i = 1, n_sps

      nl = Size(Catalog(i)%Lines)
      gl_slabs(1:no_ele,i)%no_lines = nl

      Qlog(1:3) = Catalog(i)%QLOG(1:3)

      if ( .not. Do_1D ) then

        do j = 1, no_ele

          call Slabs_Prep_Arrays ( catalog(i)%molecule, nl, t_path(j),&
            &  p_path(j), catalog(i)%mass, Qlog, &
            &  Catalog(i), gl_slabs(j,i)%v0s, gl_slabs(j,i)%x1, gl_slabs(j,i)%y, &
            &  gl_slabs(j,i)%yi,gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dslabs1_dv0 )

!  Apply velocity corrections:

          gl_slabs(j,i)%v0s = gl_slabs(j,i)%v0s * Vel_z_correction

        end do

      else

        ! compute each element along the LOS path before tangent point

        do j = 1, no_ele/2

          call Slabs_Prep_Arrays ( catalog(i)%molecule, nl, t_path(j), p_path(j), &
            & catalog(i)%mass, Qlog, &
            & Catalog(i), gl_slabs(j,i)%v0s, gl_slabs(j,i)%x1, gl_slabs(j,i)%y, &
            & gl_slabs(j,i)%yi,gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dslabs1_dv0 )

          gl_slabs(j,i)%v0s = gl_slabs(j,i)%v0s * Vel_z_correction

        end do
        
        ! fill in grid points on other side with above value
        
        do j = no_ele, no_ele/2+1, -1
          
          k = no_ele - j + 1
          gl_slabs(j,i)%v0s         = gl_slabs(k,i)%v0s
          gl_slabs(j,i)%x1          = gl_slabs(k,i)%x1
          gl_slabs(j,i)%y           = gl_slabs(k,i)%y
          gl_slabs(j,i)%yi          = gl_slabs(k,i)%yi 
          gl_slabs(j,i)%slabs1      = gl_slabs(k,i)%slabs1 
          gl_slabs(j,i)%dslabs1_dv0 = gl_slabs(k,i)%dslabs1_dv0
        
        end do
        
      end if

    end do              ! On i

  end subroutine Get_GL_Slabs_Arrays

!=====================================================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SLABS_SW_M

! $Log$
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
