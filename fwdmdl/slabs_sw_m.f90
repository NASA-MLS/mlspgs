module SLABS_SW_M

  ! Single-Line Absorption Software

  use MLSCommon, only: R8, RP
  use Units, only: SqrtPi

  implicit NONE

  private

  ! Routines to compute Betas and their derivatives, and to get ready to do so:
  public :: Get_GL_Slabs_Arrays,                                          &
         &  Slabs, Slabs_dT, Slabs_Lines, Slabs_Lines_dT, Slabs_Prep,     &
         &  Slabs_Prep_Struct, Slabs_Prep_dT, Slabswint, Slabswint_dT,    &
         &  Slabswint_Lines, Slabswint_Lines_dT, Voigt_Lorentz,           &
         &  DVoigt_Spectral, DVoigt_Spectral_Lines

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

  ! --------------------------------------------  dVoigt_spectral  -----
  subroutine dVoigt_spectral ( dNu, Nu0, x1, yi, y, w, t, tanh1, slabs1, SwI, &
                         &  dslabs1_dNu0, dSwI_dw, dSwI_dn, dSwI_dNu0 )

! Compute the Voigt function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

! NOTE: In here and in all other routines in this module, 
!       tanh1 = tanh(h * nu / ( 2.0 * k * T ) )

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dnu, nu0
    real(rp), intent(in) :: x1, yi, y, w, t, tanh1, slabs1             
    real(rp), intent(in), optional :: dslabs1_dNu0                        

    real(rp), intent(out) :: SwI, dSwI_dw,dSwI_dn,dSwI_dNu0               

    real(rp) :: x, u, v, du_dx, du_dy, dv_dx, dv_dy, q, b, g, z, r    

    real(rp) :: dx_dv0, du_dv0, dv_dv0, vvw, slabs2

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
!    db_dv0 = x1
!    dg_dv0 = 2.0_rp * b*db_dv0
!    dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
!    dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
!    dvvw_dv0 = (du_dv0+dr_dv0)*q + dq_dv0*(u+r)
!    dvvw_dv0 = (du_dv0+dr_dv0)*q
!    if (present(dslabs1_dNu0)) dSwI_dNu0 = dslabs1_dNu0*vvw + slabs2*dvvw_dv0
    if (present(dslabs1_dNu0)) dSwI_dNu0 = swi * (dslabs1_dNu0/slabs1 &
                   - 1.0_r8/Nu0) + slabs2*q*(du_dv0 + yi * dv_dv0)

  end subroutine dVoigt_spectral

  ! --------------------------------------  dVoigt_spectral_Lines  -----
  subroutine dVoigt_spectral_Lines ( dNu, Slabs, t, tanh1, &
                                  &  SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 )

! Compute the sums of the Voigt function and its first derivatives with respect
! to spectral parameters: w, n & Nu0 for all lines in Slabs.

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

! NOTE: In here and in all other routines in this module, 
!       tanh1 = tanh( h * nu / ( 2.0 * k * T ) )

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dnu
    type(slabs_struct), intent(in) :: Slabs
    real(rp), intent(in) :: T
    real(rp), intent(in) :: Tanh1

    real(rp), intent(out) :: SwI, dSwI_dw, dSwI_dn, dSwI_dNu0               

    integer, pointer :: CatLines(:)  ! slabs%catalog%lines
    integer :: L  ! Line index

    real(r8) :: Nu0
    real(rp) :: x1, yi, y, w, slabs1             
    real(rp) :: dslabs1_dNu0                        

    real(rp) :: x, u, v, SwI1, du_dx, du_dy, dv_dx, dv_dy, q, b, g, z, r    

    real(rp) :: dx_dv0, du_dv0, dv_dv0, vvw, slabs2

    catLines => slabs%catalog%lines
    SwI = 0.0_rp
    dSwI_dw = 0.0_rp
    dSwI_dn = 0.0_rp
    dSwI_dNu0 = 0.0_rp     
    do l = 1, size(catLines)

      nu0 = slabs%v0s(l)
      x1 = slabs%x1(l)
      yi = slabs%yi(l)
      y = slabs%y(l)
      w = lines(catLines(l))%w
      slabs1 = slabs%slabs1(l)
      dslabs1_dNu0 = slabs%dslabs1_dv0(l)

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
      SwI1 = slabs2 * vvw

      du_dx = 2.0_rp * (y * v - x * u)
      du_dy = 2.0_rp * (y * u + x * v - OneOvSPi)

      dv_dx = -du_dy         ! Cauchy-Riemann equation
      dv_dy =  du_dx         ! Cauchy-Riemann equation

  ! Compute the derivative of SwI w.r.t. w

      dSwI_dw = dSwI_dw + q * slabs2* (y/w) * (du_dy + yi*du_dx + &
                                  &   OneOvSPi*((1.0_rp-2.0_rp*z*y)/g))

  ! Compute the derivative of SwI w.r.t. n

      dSwI_dn = dSwI_dn + q * slabs2 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)

  ! Finaly, compute the derivative of SwI w.r.t. Nu0

  ! ***** Analytically *****

  !    dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
      dx_dv0 = -x1
      du_dv0 = du_dx * dx_dv0
      dv_dv0 = dv_dx * dx_dv0
  !    db_dv0 = x1
  !    dg_dv0 = 2.0_rp * b*db_dv0
  !    dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
  !    dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
  !    dvvw_dv0 = (du_dv0+dr_dv0)*q + dq_dv0*(u+r)
  !    dvvw_dv0 = (du_dv0+dr_dv0)*q
  !    dSwI_dNu0 = dslabs1_dNu0*vvw + slabs2*dvvw_dv0
      dSwI_dNu0 = dSwI_dNu0 + swi1 * (dslabs1_dNu0/slabs1 &
                     - 1.0_r8/Nu0) + slabs2*q*(du_dv0 + yi * dv_dv0)

      SwI = SwI + SwI1

    end do

  end subroutine dVoigt_spectral_Lines

  ! ------------------------------------------------------  Slabs  -----
  real(rp) function Slabs ( Nu, v0, v0s, x1, tanh1, slabs1, y )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    real(r8), intent(in) :: v0    ! Line center frequency
    real(r8), intent(in) :: v0s   ! Pressure-shifted line center frequency
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: slabs1
    real(rp), intent(in) :: y

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: u

    call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta, b = x_1 \sigma$, and
!  $D = \frac1{b^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + V(b,y) \right)$.  $b$ is always huge, so we approximate
!  $V(b,y)$ with one term of an asymptotic expansion, \emph{viz.}
!  $V(b,y) \sim \frac{y D}{\sqrt{\pi}}$.

    Slabs = slabs1 * real(nu / v0, rp) * tanh1 * &
      & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

  end function Slabs

  ! ----------------------------------------------------  Slabs_dT  -----
  subroutine Slabs_dT ( Nu, v0, v0s, x1, tanh1, slabs1, y, &
    &                           dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, dy_dT, &
    &                   Slabs, dSlabs_dT )

  ! Compute single-line absorption and its derivative w.r.t. temperature.

    use Voigt_m, only: D_Real_Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(out) :: Slabs, dSlabs_dT

    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu-v0s
    real(rp) :: Du      ! du/dT
    real(rp) :: Da      ! d(x1*delta)
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: Sigma   ! Nu+v0s
    real(rp) :: SigmaX1 ! sigma * x1
    real(rp) :: U       ! Voigt
    real(rp) :: Y2      ! y**2

    delta = Nu - v0s
    da = x1 * ( delta * dx1_dT - dv0s_dT ) ! remember, dx1_dT = 1/x1 dx1 / dT
    call D_Real_Simple_Voigt ( x1*delta, y, da, y*dy_dT, u, du )

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta, b = x_1 \sigma$, and
!  $D = \frac1{b^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + V(b,y) \right)$.  $b$ is always huge, so we approximate
!  $V(b,y)$ with one term of an asymptotic expansion, \emph{viz.}
!  $V(b,y) \sim \frac{y D}{\sqrt{\pi}}$.

    sigma = nu + v0s
    sigmaX1 = sigma * x1
    y2 = y*y
    d = 1.0_rp / ( sigmaX1**2 + y2 )
    sb = slabs1 * real(nu / v0,rp) * tanh1
    sa = sb * u
    sb = sb * OneOvSPi * y * d
    Slabs = sa + sb

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as $V(a,y) +
!  i L(a,y)$, where $V(a,y)$ is the Voigt function ({\tt u} above) and $L(a,y)$
!  is the Lorentz function.  All we want for $S$ is its real part, so we don't
!  need $L(a,y)$ for $S$.  For $\frac{\partial S}{\partial T}$ we want the real
!  part of the derivative; this requires the real part of the derivative of
!  $w(z)$, not just Voigt.
!
!  Write $S = S_a + S_b$ where
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) V(a,y)$ and
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{y D}{\sqrt{\pi}}$.\\
!  Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$ and\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!    \frac1y \frac{\partial y}{\partial T} - 
!    2 D \left( x_1 \sigma \left( x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!     \sigma \frac{\partial x_1}{\partial T} \right) +
!     y \frac{\partial y}{\partial T} \right)$.\\
!    Notice that the first two terms of
!    $\frac1{S_a}\frac{\partial S_a}{\partial T}$ and
!    $\frac1{S_b}\frac{\partial S_b}{\partial T}$ are the same.\\
!  For $\Re \frac{\partial w(z)}{\partial T}$ we need
!  $\frac{\partial a}{\partial T} = -x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!   \delta \frac{\partial x_1}{\partial T}$ and $\frac{\partial y}{\partial T}$.

    c = dSlabs1_dT + dtanh_dT
    dSlabs_dT = sa * ( c + du / u ) + &
      &         sb * ( c + dy_dT - 2.0_rp * D * &
      &           ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

  end subroutine Slabs_dT

  ! ------------------------------------------------  Slabs_Lines  -----
  function Slabs_Lines ( Nu, Slabs, tanh1, NoPolarized ) result ( Beta )

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    real(rp) :: Beta ! Output

! If the molecular transition and temperature have not changed but
! frequency has, enter here to calculate sum of Beta for all lines in SLABS.

    integer, pointer :: CatLines(:)  ! slabs%catalog%lines
    integer :: L      ! Line Index
    real(rp) :: U     ! Voigt = real part of Fadeeva
    real(r8) :: V0S   ! Pressure-shifted line center
    real(r8) :: X1, Y ! Doppler width, ratio Pressure to Doppler widths

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, and
!  $D = \frac1{\sigma^2 x_1^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + \frac{y D}{\sqrt{\pi}} \right)$.

    beta = 0.0_rp
    catLines => slabs%catalog%lines
    if ( .not. noPolarized ) then
      do l = 1, size(catLines)
        v0s = slabs%v0s(l)
        x1 = slabs%x1(l)
        y = slabs%y(l)
        call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

        beta = beta + slabs%slabs1(l) * &
          &           real(nu / lines(catLines(l))%v0, rp) * tanh1 * &
          & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

      end do
    else
      do l = 1, size(slabs%v0s)
        if ( slabs%catalog%polarized(l) ) cycle
        v0s = slabs%v0s(l)
        x1 = slabs%x1(l)
        y = slabs%y(l)
        call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

        beta = beta + slabs%slabs1(l) * &
          &           real(nu / lines(catLines(l))%v0, rp) * tanh1 * &
          & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

      end do
    end if

  end function Slabs_Lines

  ! ---------------------------------------------  Slabs_Lines_dT  -----
  subroutine Slabs_Lines_dT ( Nu, Slabs, Tanh1, dTanh_dT, &
    &                         Beta, dBeta_dT, NoPolarized )

  ! Compute single-line absorption and its derivative w.r.t. temperature
  ! for all lines in the Slabs structure.

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: D_Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: Tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: dTanh_dT ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT

    real(rp), intent(out) :: Beta, dBeta_dT

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    integer, pointer :: CatLines(:)  ! slabs%catalog%lines
    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu-v0s
    real(rp) :: Du      ! du/dT
    real(rp) :: Da      ! d(x1*delta)
    real(rp) :: Dv0s_dT ! dv0s / dT
    real(rp) :: Dx1_dT  ! 1/x1 dx1/dT
    real(rp) :: Dy_dT   ! 1/y dy/dT
    integer :: L        ! Line index
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: Sigma   ! Nu+v0s
    real(rp) :: SigmaX1 ! sigma * x1
    real(rp) :: U       ! Voigt
    real(r8) :: V0S     ! Pressure-shifted line center
    real(rp) :: X1, Y   ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Y2      ! y**2

! See Slabs_dT for TeXnicalities

    catLines => slabs%catalog%lines
    Beta = 0.0_rp
    dBeta_dT = 0.0_rp

    do l = 1, size(catLines)

      if ( noPolarized ) then
        if ( slabs%catalog%polarized(l) ) cycle
      end if

      v0s = slabs%v0s(l)
      x1 = slabs%x1(l)
      y = slabs%y(l)
      dv0s_dT = slabs%dv0s_dT(l)
      dx1_dT = slabs%dx1_dT(l)
      dy_dT =slabs%dy_dT(l)
      delta = Nu - v0s
      da = x1 * ( delta * dx1_dT - dv0s_dT )
      call D_Real_Simple_Voigt ( x1*delta, y, da, y*dy_dT, u, du )
      sigma = nu + v0s
      sigmaX1 = sigma * x1
      y2 = y * y
      d = 1.0_rp / ( sigmaX1**2 + y2 )
      sb = slabs%slabs1(l) * real(nu / lines(catLines(l))%v0,rp) * tanh1
      sa = sb * u
      sb = sb * OneOvSPi * y * d
      beta = beta + sa + sb

      c = slabs%dSlabs1_dT(l) + dtanh_dT

      dBeta_dT = dBeta_dT + sa * ( c + du / u ) &
        &                 + sb * ( c + dy_dT - 2.0_rp * d * &            
        &                          ( sigmaX1 * ( x1 * dv0s_dT + &        
        &                            sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

    end do

  end subroutine Slabs_Lines_dT

  ! --------------------------------------------------  Slabswint  -----
  real(rp) function Slabswint ( Nu, v0, v0s, x1, tanh1, slabs1, y, yi )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    real(r8), intent(in) :: v0    ! Line center frequency
    real(r8), intent(in) :: v0s   ! Pressure-shifted line center frequency
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: slabs1
    real(rp), intent(in) :: y
    real(rp), intent(in) :: yi

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: a, sigmaX1, u, y2

    a = x1 * real(nu-v0s,rp)
    call real_simple_voigt ( a, y, u )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$,
!  $D_1 = \frac1{\sigma^2 x_1^2 + y^2}$ and $D_2 = \frac1{a^2 + y^2}$.
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!  \tanh\left(\frac{h \nu}{2 k T}\right)
!   \left( V(a,y) + \frac{(y - \sigma x_1 y_i) D_1}{\sqrt{\pi}}
!     + \frac{y_i a D_2}{\sqrt{\pi}} \right)$.

    sigmaX1 = x1 * (nu + v0s)
    y2 = y*y
    Slabswint = slabs1 * real(Nu / v0, rp) * tanh1 * &
      & (u + OneOvSPi*((y - sigmaX1*yi)/(sigmaX1*sigmaX1 + y2) + yi*a/(a*a+y2)))

  end function Slabswint

  ! -----------------------------------------------  Slabswint_dT  -----
  subroutine Slabswint_dT ( Nu, v0, v0s, x1, tanh1, slabs1, y, yi, &
    &                              dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, &
    &                              dy_dT, dyi_dT, &
    &                              Slabswint, dSlabs_dT )

    use Voigt_m, only: D_Real_Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y, yi
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(in) :: dyi_dT     ! 1/yi dyi / dT

    real(rp), intent(out) :: Slabswint, dSlabs_dT

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: A               ! x1 * delta
    real(rp) :: C               ! Terms common to the parts of dSlabs_dT
    real(rp) :: D1              ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: D2              ! 1 / (a**2 + y**2)
    real(rp) :: DD1, DD2        ! 1/D1 d(D1)/dT, 1/D2 d(D2)/dT, 
    real(rp) :: Da              ! da / dT
    real(rp) :: Delta           ! Nu-v0s
    real(rp) :: DU              ! du/dT
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of SlabsWint
    real(rp) :: Sigma           ! Nu+v0s
    real(rp) :: SigmaX1         ! sigma * x1
    real(rp) :: U               ! Voigt
    real(rp) :: Y2              ! Y**2

    delta = nu - v0s
    a = x1 * delta
    da = x1 * ( delta * dx1_dT - dv0s_dT ) ! remember, dx1_dT = 1/x1 dx1 / dT
    call D_Real_Simple_Voigt ( a, y, da, y*dy_dT, u, du )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$,
!  $D_1 = \frac1{\sigma^2 x_1^2 + y^2}$ and $D_2 = \frac1{a^2 + y^2}$.
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!  \tanh\left(\frac{h \nu}{2 k T}\right)
!   \left( V(a,y) + \frac{(y - \sigma x_1 y_i) D_1}{\sqrt{\pi}}
!     + \frac{y_i a D_2}{\sqrt{\pi}} \right)$.

    sigma = nu + v0s
    sigmaX1 = sigma * x1
    y2 = y * y
    d1 = 1.0_rp / ( sigmaX1**2 + y2 )
    d2 = 1.0_rp / ( a * a + y2 )
    c = slabs1 * real(nu / v0,rp) * tanh1
    sa = c * u
    c = c * OneOvSPi
    sb = c * y * d1
    sc = c * sigmaX1 * yi * d1
    sd = c * yi * d2
    Slabswint = sa + sb - sc + sd * a

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as $V(a,y) +
!  i L(a,y)$, where $V(a,y)$ is the Voigt function ({\tt u} above) and $L(a,y)$
!  is the Lorentz function.  All we want for $S$ is its real part, so we don't
!  need $L(a,y)$ for $S$. For $\frac{\partial S}{\partial T}$ we want the real
!  part of the derivative; this requires the real part of the derivative of
!  $w(z)$, not just Voigt.
!
!  Write $S = S_a + S_b - S_c + S_d$ where
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) V(a,y)$,
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{y D_1}{\sqrt{\pi}}$,
!  $S_c = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{\sigma x_1 y_i D_1}{\sqrt{\pi}}$, and
!  $S_d = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{a y_i D_2}{\sqrt{\pi}}$.\\
!  Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T} - \frac{\partial S_c}{\partial T} +
!   \frac{\partial S_d}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$,\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{D_1}\frac{\partial D_1}{\partial T} +
!   \frac1y \frac{\partial y}{\partial T}$,\\
!  $\frac1{S_c}\frac{\partial S_c}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{D_1}\frac{\partial D_1}{\partial T} +
!   \frac1{y_i}\frac{\partial y_i}{\partial T} +
!   \frac1{x_1}\frac{\partial x_1}{\partial T} +
!   \frac1{\sigma}\frac{\partial \nu_{0_s}}{\partial T}$, and\\
!  $\frac1{S_d}\frac{\partial S_d}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{D_2}\frac{\partial D_2}{\partial T} +
!   \frac1{y_i}\frac{\partial y_i}{\partial T} +
!   \frac1{a}\frac{\partial a}{\partial T}$
!   , where\\
!  $\frac1{D_1}\frac{\partial D_1}{\partial T} = -2 D_1 \left( x_1 \sigma \left( x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!     \sigma \frac{\partial x_1}{\partial T} \right) +
!     y \frac{\partial y}{\partial T} \right)$ and
!  $\frac1{D_2}\frac{\partial D_2}{\partial T} = -2 D_2 \left( a \frac{\partial a}{\partial T} +
!     y \frac{\partial y}{\partial T} \right)$.\\
!  Notice that the first two terms of
!    $\frac1{S_a}\frac{\partial S_a}{\partial T}$,
!    $\frac1{S_b}\frac{\partial S_b}{\partial T}$,
!    $\frac1{S_c}\frac{\partial S_c}{\partial T}$ and
!    $\frac1{S_d}\frac{\partial S_d}{\partial T}$ are the same,
!  the third terms of
!    $\frac1{S_b}\frac{\partial S_b}{\partial T}$ and
!    $\frac1{S_c}\frac{\partial S_c}{\partial T}$ are the same, and
!  the fourth terms of
!    $\frac1{S_c}\frac{\partial S_c}{\partial T}$ and
!    $\frac1{S_d}\frac{\partial S_d}{\partial T}$ are the same.\\
!  For $\Re \frac{\partial w(z)}{\partial T}$ we need
!  $\frac{\partial a}{\partial T} = -x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!   \delta \frac{\partial x_1}{\partial T}$ and $\frac{\partial y}{\partial T}$.

    c = dSlabs1_dT + dtanh_dT
    dd1 = -2.0_rp * d1 * ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
    dd2 = -2.0_rp * d2 * ( a * da + y2 * dy_dT )
    dSlabs_dT = sa * ( c + du / u ) &
      &       + sb * ( c + dd1 + dy_dT ) &
      &       - sc * ( c + dd1 + dyi_dT + dx1_dT + dv0s_dT / sigma ) &
      &       + sd * ( a * ( c + dd2 + dyi_dT ) + da )

  end subroutine Slabswint_dT

  ! --------------------------------------------  Slabswint_Lines  -----
  function Slabswint_Lines ( Nu, Slabs, tanh1, NoPolarized ) result ( Beta )

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    real(rp) :: Beta ! Output

! If the molecular transition and temperature have not changed but
! frequency has, enter here to calculate sum of Beta for all lines in SLABS.

    real(rp) :: A        ! First argument for real_simple_voigt
    integer, pointer :: CatLines(:)  ! slabs%catalog%lines
    integer :: L         ! Line Index
    real(rp) :: SigmaX1  ! (nu + nu0s) * x1
    real(rp) :: U        ! Voigt = real part of Fadeeva
    real(r8) :: V0S      ! Pressure-shifted line center
    real(rp) :: X1, Y    ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Yi       ! Interference coefficient
    real(rp) :: Y2       ! y**2

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$,
!  $D_1 = \frac1{\sigma^2 x_1^2 + y^2}$ and $D_2 = \frac1{a^2 + y^2}$.
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!  \tanh\left(\frac{h \nu}{2 k T}\right)
!   \left( V(a,y) + \frac{(y - \sigma x_1 y_i) D_1}{\sqrt{\pi}}
!     + \frac{y_i a D_2}{\sqrt{\pi}} \right)$.

    beta = 0.0_rp
    catLines => slabs%catalog%lines
    do l = 1, size(slabs%v0s)
      if ( noPolarized .and. slabs%catalog%polarized(l) ) cycle
      v0s = slabs%v0s(l)
      x1 = slabs%x1(l)
      y = slabs%y(l)
      yi = slabs%yi(l)
      a = x1 * real(nu-v0s,rp)
      call real_simple_voigt ( a, y, u )

      sigmaX1 = x1 * (nu + v0s)
      y2 = y*y
      if ( abs(yi) > 1.0e-6_rp ) then ! Include interference effect
        beta = beta + slabs%slabs1(l) * &
          &           real(Nu / lines(catLines(l))%v0, rp) * tanh1 * &
          & (u + OneOvSPi*((y - sigmaX1*yi)/(sigmaX1*sigmaX1 + y2) + yi*a/(a*a+y2)))
      else
        beta = beta + slabs%slabs1(l) * &
          &           real(Nu / lines(catLines(l))%v0, rp) * tanh1 * &
          & (u + OneOvSPi*(y/(sigmaX1*sigmaX1 + y2)))
      end if
    end do

  end function Slabswint_Lines

  ! -----------------------------------------  Slabswint_Lines_dT  -----
  subroutine Slabswint_Lines_dT ( Nu, Slabs, tanh1, dTanh_dT, &
    &                             Beta, dBeta_dT, NoPolarized )

  ! Compute single-line absorption and its derivative w.r.t. temperature,
  ! with interference, for all lines in the Slabs structure.

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: D_Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: Tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: dTanh_dT ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT

    real(rp), intent(out) :: Beta, dBeta_dT

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: A               ! x1 * delta
    real(rp) :: C1, C2          ! Common terms
    integer, pointer :: CatLines(:)  ! slabs%catalog%lines
    real(rp) :: D1              ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: D2              ! 1 / (a**2 + y**2)
    real(rp) :: DD1, DD2        ! 1/D1 d(D1)/dT, 1/D2 d(D2)/dT, 
    real(rp) :: Da              ! da / dT
    real(rp) :: Delta           ! Nu-v0s
    real(rp) :: DU              ! du / dT
    real(rp) :: Dv0s_dT         ! dv0s / dT
    real(rp) :: Dx1_dT          ! 1/x1 dx1/dT
    real(rp) :: Dy_dT           ! 1/y dy/dT
    real(rp) :: Dyi_dT          ! 1/yi dyi/dT
    integer :: L                ! Line index
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of SlabsWint
    real(rp) :: Sigma           ! Nu+v0s
    real(rp) :: SigmaX1         ! sigma * x1
    real(rp) :: U               ! Voigt
    real(r8) :: V0S             ! Pressure-shifted line center
    real(rp) :: X1, Y           ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Yi              ! Interference coefficient
    real(rp) :: Y2              ! Y**2

    ! See Slabswint_dT for TeXnicalities.

    catLines => slabs%catalog%lines
    beta = 0.0_rp
    dBeta_dT = 0.0_rp
    do l = 1, size(catLines)

      if ( noPolarized .and. slabs%catalog%polarized(l) ) cycle

      v0s = slabs%v0s(l)
      x1 = slabs%x1(l)
      y = slabs%y(l)
      yi = slabs%yi(l)
      dv0s_dT = slabs%dv0s_dT(l)
      dx1_dT = slabs%dx1_dT(l)
      dy_dT = slabs%dy_dT(l)
      dyi_dT = slabs%dyi_dT(l)
      delta = nu - v0s
      a = x1 * delta
      da = x1 * ( delta * dx1_dT - dv0s_dT )
      call D_Real_Simple_Voigt ( a, y, da, y*dy_dT, u, du )

      sigma = nu + v0s
      sigmaX1 = sigma * x1
      y2 = y * y
      d1 = 1.0_rp / ( sigmaX1**2 + y2 )
      d2 = 1.0_rp / ( a * a + y2 )
      c1 = slabs%slabs1(l) * real(nu / lines(catLines(l))%v0,rp) * &
        & tanh1
      sa = c1 * u
      c1 = c1 * OneOvSPi
      sb = c1 * y * d1
      c2 = slabs%dSlabs1_dT(l) + dtanh_dT
      dd1 = -2.0_rp * d1 * ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
      if ( abs(yi) > 1.0e-6_rp ) then
        sc = c1 * sigmaX1 * yi * d1
        sd = c1 * yi * d2
        beta = beta + sa + sb - sc + sd * a

        dd2 = -2.0_rp * d2 * ( a * da + y2 * dy_dT )

        dBeta_dT = dBeta_dT &
          &      + sa * ( c2 + du / u ) &
          &      + sb * ( c2 + dd1 + dy_dT ) &
          &      - sc * ( c2 + dd1 + dyi_dT + dx1_dT + dv0s_dT / sigma ) &
          &      + sd * ( a * (c2 + dd2 + dyi_dT) + da )
      else
        beta = beta + sa + sb
        dBeta_dT = dBeta_dT + sa * ( c2 + du / u )
      end if

    end do

  end subroutine Slabswint_Lines_dT

  ! ----------------------------------------------  Voigt_Lorentz  -----

  subroutine Voigt_Lorentz ( dNu,  Nu0,  x1,  yi,  y,  w,  t,  tanh1, slabs1,  &
                         &   VL, dslabs1_dNu0,  dVL_dw,  dVL_dn,  dVL_dNu0 )

! Compute the Voigt/Lorentz function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dNu, nu0
    real(rp), intent(in) :: x1, yi, y, w, t, tanh1, slabs1, dslabs1_dNu0

    real(rp), intent(out) :: VL, dVL_dw, dVL_dn, dVL_dNu0

    real(rp) :: xj, zj, q, y2, u, v, up1, up2, dn1, dn2, dup1, &
     &          dup2, ddn1, ddn2, dy_dw, dy_dn, dSum_dw, dSum_dn, &
     &          slabs2

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
    slabs2 = slabs1 * tanh1
    VL = slabs2 * ddn2 * q            ! This is the Voigt + VVW correction

    dn2 = xj * xj + y2
    up2 = y + yi * xj

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

!    dup2 =  yi * x1               !  x1 = -dxj_dNu0
!    ddn2 = -2.0 * xj * x1         !  x1 = -dxj_dNu0
!    dSum_dNu0 = (dn2*dup2-up2*ddn2)/(dn2*dn2)
!    dq_dNu0 = -(dNu+Nu0)/(Nu0*Nu0)

!    dVL_dNu0 = OneOvSPi * (dslabs1_dNu0 * q * Sum          + &
!              &         slabs2 * dq_dNu0 * Sum + &
!              &         slabs2 * q * dSum_dNu0)

!    dVL_dNu0 = OneOvSPi * slabs2 * q * dSum_dNu0
    dVL_dNu0 = VL * (dslabs1_dNu0 / slabs1 - 1.0_rp / Nu0) &
           & + OneOvSPi * slabs2 * q * ( &
           &  (-yi*x1*(xj*xj+y2)+2.0_rp*(y+yi*xj)*xj*x1)/(xj*xj+y2)**2 &
           & + (yi*x1*  dn1     -2.0_rp*   up1   *x1*zj)/dn1**2)

  end subroutine Voigt_Lorentz

  ! -------------------------------------------------  Slabs_prep  -----
  subroutine Slabs_prep ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                      &   n1, n2, velCor, useYi, &
                      &   v0s, x1, y, yi, slabs1, dslabs1 )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! >>2004-03-18 WV Snyder Use 1 - exp(-v0/300/Boltzmhz) in denominator instead
!                        of 1 - exp(-v0s/300/Boltzmhz).

    use Physics, only: H_OVER_K, k, SpeedOfLight
    use Units, only: Ln10, Sqrtln2, SqrtPi

! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(r8), intent(in) :: M        ! Molecular mass amu
    real(r8), intent(in) :: V0       ! Line center frequency MHz
    real(r8), intent(in) :: El       ! Lower state energy cm-1
    real(r8), intent(in) :: W        ! Collision broadening parameter
                                     ! MHz/mbar at 300 K
    real(r8), intent(in) :: Ps       ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P        ! Pressure mbar
    real(r8), intent(in) :: N        ! Temperature power dependence of w
    real(r8), intent(in) :: Ns       ! Temperature power dependence of ps
    real(r8), intent(in) :: I        ! Integrated spectral intensity
                                     ! Log(nm**2 MHz) at 300 K
    real(r8), intent(in) :: Q(3)     ! Logarithm of the partition function
                                     ! At 300 , 225 , and 150 K
    real(r8), intent(in) :: Delta    ! Delta interference coefficient at 300K 1/mb
    real(r8), intent(in) :: Gamma    ! Gamma               "
    real(r8), intent(in) :: N1       ! Temperature dependency of delta
    real(r8), intent(in) :: N2       ! Temperature dependency of gamma
    real(r8), intent(in) :: VelCor   ! Doppler velocity correction term
    logical, intent(in) :: UseYi     ! delta + gamma > 0.0

! outputs:

    real(r8), intent(out) :: V0s     ! Pressure shifted line position
    real(r8), intent(out) :: X1      ! Sqrt(Ln(2))/Doppler half width MHz
    real(r8), intent(out) :: Y       ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
    real(r8), intent(out) :: Yi      ! Interference contribution
    real(r8), intent(out) :: Slabs1  ! Frequency independent piece of slabs
    real(r8), intent(out) :: Dslabs1 ! Derivative of slabs1 w.r.t. v0

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K
!  sqrtln2  - sqrt(ln(2))

    real(rp), parameter :: I2abs = sqrtln2 / ( sqrtPi * 1.0e13 * k )
!   real(rp), parameter :: I2abs = 3.402155052e9_rp ! using above constants
!   real(rp), parameter :: I2abs = 3.402136078e9_rp ! Zvi's original value
    real(rp), parameter :: Dc = 3.58116514e-7_rp ! sqrt(1000 k ln 4 avogadro) / c
!   real(rp), parameter :: Dc = 3.58117369e-7_rp ! Zvi's original value
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Oned300 = 1.0_rp/300.0_rp

    real(rp), parameter :: LT2 = 2.35218251811136_rp      ! Log10(225)
    real(rp), parameter :: LT3 = 2.47712125471966_rp      ! Log10(300)

    real(rp), parameter :: Tl1 = 0.176091259055681_rp     ! Log10(225/150)
    real(rp), parameter :: Tl2 = 0.124938736608300_rp     ! Log10(300/225)

! Internal data:

    real(rp) :: betae, betav, expd, expn, log_T, onedt
    real(rp) :: Q_Log_a, Q_Log_b, t3t, Wd, z1, z2

! The action begins here

    onedt = 1.0_rp / t
    log_T = log(t)
    t3t = lt3 * ln10 - log_T ! log(300/T)

!{ $y_i = p \left( \delta \left( \frac{300}T \right)^{n_1} +
!                  \gamma \left( \frac{300}T \right)^{n_2} \right)$.

    yi = 0.0
    if ( useYi ) yi = p * (delta*exp(n1*t3t) + gamma*exp(n2*t3t))

!{ $\nu_{0_s} = v_c \left[ \nu_0 + p_s p \left( \frac{300}T \right) ^{n_s} \right]$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = v_c$.

    v0s = velCor * ( v0 + ps * p * exp(ns*t3t) )

!{ $w_d = \nu_0 d_c \sqrt{\frac{T}M}$.

    Wd = v0 * dc * Sqrt(t/m)

!{ $x_1 = \frac{\sqrt{\ln 2}}{w_d}$.
!  $\frac{\partial x_1}{\partial \nu_0} = -\frac{x_1}{\nu_0}$.

    x1 = real(sqrtln2,rp) / Wd

!{ $y = x_1 w p \left( \frac{300}T \right) ^n$.
!  $\frac{\partial y}{\partial \nu_0} = -\frac{y}{\nu_0}$.

    y = x1 * w * p * exp(n*t3t)

!{ $\beta_e = 10^{-4} \frac{h}k c\, e_l$. $\beta_v = \frac{h}k \nu_{0_s}$.

    betae = el / boltzcm
    betav = v0s / boltzmhz

    if ( t < 225.0_rp ) then
      q_log_b = (q(2)-q(3)) / tl1
      q_log_a = q(2) - q(1) - lt2 * q_log_b
    else
      q_log_b = (q(1)-q(2)) / tl2
      q_log_a = -lt3 * q_log_b
    end if
    q_log_b = -1.0 - q_log_b ! This is more useful below

    expd = EXP(-v0*(oned300/boltzmhz)) ! H
    expn = EXP(-betav*onedt) ! G
    z1 = 1.0 + expn          ! 1 + G
    z2 = 1.0 - expd          ! 1 - H

!{ $S = \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T w_d \left(1 - e^{-\frac{h \nu_0}{300 k}}\right)} =
!  I_2 \, p \, e^{(i-a) \ln 10 -(b+1) \log T +
!    \beta_e \left( \frac1{300} - \frac1T \right)}
!  \frac{(1+G)} {w_d (1-H)}$, where
!  $G = e^{-\frac{\beta_v}T}$, $H = e^{-\frac{h \nu_0}{300 k}}$ and
!  $S =$ {\tt slabs1}.

    slabs1 = i2abs * p * &
      & exp((i - q_log_a)*ln10 + q_log_b * log_T + betae * (oned300 - onedt)) &
      & * z1 / (Wd * z2)

!{ $\frac{\partial S}{\partial \nu_0} =
!   -S \left[ \left( \frac{G_1}T \frac{\partial \nu_{0_s}}{\partial \nu_0}
!   + \frac{H_1}{300} \right) \frac{h}k + \frac1{\nu_0} \right]$ where
!   $G_1 = \frac{G}{1 + G}$ and $H_1 = \frac{H}{1 - H}$.

    dslabs1 = -slabs1 * ( ( expn / (t * z1) * velCor + expd / (300.0_rp * z2)) / &
      & boltzmhz + 1.0 / v0 )

  end subroutine Slabs_prep

  ! ------------------------------------------  Slabs_prep_struct  -----
  subroutine Slabs_prep_struct ( T, P, Catalog, VelCor, Derivs, Slabs )
  ! Fill all the fields of the Slabs structure

    use L2PC_PFA_STRUCTURES, only: Slabs_Struct
    use SpectroscopyCatalog_m, only: Catalog_T, Lines

    ! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(rp), intent(in) :: P        ! Pressure
    type(catalog_t), intent(in) :: Catalog ! The spectroscopy
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction term, 
                                     ! 1 - losVel / C
    logical, intent(in) :: Derivs    ! "Setup for derivative calculations"

    ! output:

    type(slabs_struct), intent(inout) :: Slabs ! inout so as not to clobber
                                     ! pointer associations

    integer :: I ! A loop index
    integer :: L ! Index in the Lines array

    slabs%useYi = .false.
    do i = 1, size(catalog%lines)
      slabs%dx1_dv0(i) = 0.0
      slabs%dy_dv0(i) = 0.0
      l = catalog%lines(i)
      slabs%useYi = slabs%useYi .or. lines(l)%useYi
      if ( derivs ) then
        call slabs_prep_dT ( t, catalog%mass, &
          & lines(l)%v0, lines(l)%el, lines(l)%w, lines(l)%ps, p, &
          & lines(l)%n, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%v0s(i), slabs%x1(i), slabs%y(i), &
          & slabs%yi(i), slabs%slabs1(i), &
          & slabs%dslabs1_dv0(i), &
          & slabs%dv0s_dT(i), slabs%dx1_dT(i), &
          & slabs%dy_dT(i), slabs%dyi_dT(i), &
          & slabs%dslabs1_dT(i) )
      else
        call slabs_prep ( t, catalog%mass, &
          & lines(l)%v0, lines(l)%el, lines(l)%w, lines(l)%ps, p, &
          & lines(l)%n, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%v0s(i), slabs%x1(i), slabs%y(i), &
          & slabs%yi(i), slabs%slabs1(i), &
          & slabs%dslabs1_dv0(i) )
      end if
    end do ! i = 1, size(catalog%lines)

  end subroutine Slabs_prep_struct

  ! ----------------------------------------------  Slabs_prep_DT  -----
  subroutine Slabs_prep_dT ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                         &   n1, n2, velCor, useYi, &
                         &   v0s, x1, y, yi, slabs1, dslabs1_dv0, &
                         &   dv0s_dT, dx1_dT, dy_dT, dyi_dT, dslabs1_dT )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
! Compute the derivatives with respect to temperature, too.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects

    use Physics, only: H_OVER_K, k, SpeedOfLight
    use Units, only: Ln10, Sqrtln2, SqrtPi

! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(r8), intent(in) :: M        ! Molecular mass amu
    real(r8), intent(in) :: V0       ! Line center frequency MHz
    real(r8), intent(in) :: El       ! Lower state energy cm-1
    real(r8), intent(in) :: W        ! Collision broadening parameter
                                     ! MHz/mbar at 300 K
    real(r8), intent(in) :: Ps       ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P        ! Pressure mbar
    real(r8), intent(in) :: N        ! Temperature power dependence of w
    real(r8), intent(in) :: Ns       ! Temperature power dependence of ps
    real(r8), intent(in) :: I        ! Integrated spectral intensity
                                     ! Log(nm**2 MHz) at 300 K
    real(r8), intent(in) :: Q(3)     ! Logarithm of the partition function
                                     ! At 300 , 225 , and 150 K
    real(r8), intent(in) :: Delta    ! Delta interference coefficient at 300K 1/mb
    real(r8), intent(in) :: Gamma    ! Gamma               "
    real(r8), intent(in) :: N1       ! Temperature dependency of delta
    real(r8), intent(in) :: N2       ! Temperature dependency of gamma
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction term
    logical, intent(in) :: UseYi     ! delta + gamma > 0

! outputs:

    real(r8), intent(out) :: V0s     ! Pressure shifted line position
    real(r8), intent(out) :: X1      ! Sqrt(Ln(2))/Doppler half width MHz
    real(r8), intent(out) :: Y       ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
    real(r8), intent(out) :: Yi      ! Interference contribution
    real(r8), intent(out) :: Slabs1  ! Frequency independent piece of slabs
    real(r8), intent(out) :: Dslabs1_dv0 ! Derivative of slabs1 w.r.t. v0

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Derivatives with respect to temperature:

    real(r8), intent(out) :: dv0s_dT    ! dv0s/dT
    real(r8), intent(out) :: dx1_dT     ! 1/x1 dx1/dT
    real(r8), intent(out) :: dy_dT      ! 1/y dy/dT
    real(r8), intent(out) :: dyi_dT     ! 1/yi dyi/dT
    real(r8), intent(out) :: dslabs1_dT ! 1/slabs1 dslabs1/dT

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K = k/h
!  sqrtln2  - sqrt(ln(2))

    real(rp), parameter :: I2abs = sqrtln2 / ( sqrtPi * 1.0e13 * k )
!   real(rp), parameter :: I2abs = 3.402155052e9_rp ! using above constants
!   real(rp), parameter :: I2abs = 3.402136078e9_rp ! Zvi's original value
    real(rp), parameter :: Dc = 3.58116514e-7_rp ! sqrt(1000 k ln 4 avogadro) / c
!   real(rp), parameter :: Dc = 3.58117369e-7_rp ! Zvi's original value
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Oned300 = 1.0_rp/300.0_rp

    real(rp), parameter :: LT2 = 2.35218251811136_rp      ! Log10(225)
    real(rp), parameter :: LT3 = 2.47712125471966_rp      ! Log10(300)

    real(rp), parameter :: Tl1 = 0.176091259055681_rp     ! Log10(225/150)
    real(rp), parameter :: Tl2 = 0.124938736608300_rp     ! Log10(300/225)

! Internal data:

    real(rp) :: Betae, Betav, dBetav_dT, onedt, expn, expd
    real(rp) :: Q_Log_a, Q_Log_b ! Q_Log = a + b * log10(T)
    real(rp) :: Log_T         ! log(t)
    real(rp) :: t3t           ! log(300/t) = log(300) - log(t)
    real(rp) :: Wd, DWd_dT    ! dWd_dT is actually -dWd_dT/Wd
    real(rp) :: Z1, Z2 ! Temps

! The action begins here

    onedt = 1.0_rp / t
    log_t = log(t)
    t3t = lt3 * ln10 - log_t  ! log(300/T)

!{ $y_i = p \left( \delta \left( \frac{300}T \right)^{n_1} +
!                  \gamma \left( \frac{300}T \right)^{n_2} \right)$.
!  $\frac{\partial y_i}{\partial T} = -\frac{p}T \left(
!                    n_1 \delta \left( \frac{300}T \right)^{n_1} +
!                    n_2 \gamma \left( \frac{300}T \right)^{n_2} \right)$.

    if ( useYi ) then
      z1 = delta*exp(n1*t3t)
      z2 = gamma*exp(n2*t3t)
      yi = ( z1 + z2 )
      dyi_dT = -onedt * ( n1 * z1 + n2 * z2 ) / yi ! 1/yi dyi/dT
      yi = p * yi
    else ! yi == 0.0
      ! If yi == 0.0, dyi_dT will necessarily be zero.  The 1/yi cancels
      ! a yi in a numerator where dyi_dT is used, so 0.0 is the correct
      ! result.  We don't need a fancy l'Hospital argument to justify it.
      yi = 0.0
      dyi_dT = 0.0
    end if

!{ $\nu_{0_s} = v_c \left[ \nu_0 + p_s p \left( \frac{300}T \right)^{n_s} \right]$.
!  $\frac{\partial \nu_{0_s}}{\partial T} = \frac{-n_s}T ( \nu_{0_s} - v_c \nu_0 )$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = v_c$.

    if ( ps /= 0.0_r8 ) then
      v0s = velCor * ps * p * exp(ns*t3t)
      dv0s_dT = -ns * v0s * onedt
      v0s = velCor * v0 + v0s
    else
      v0s = velCor * v0
      dv0s_dT = 0.0
    end if

!{ $w_d = \nu_0 d_c \sqrt{\frac{T}M}$.  The $\nu_0$ term should
!  really be $\nu$.  We approximate $\nu$ by $\nu_0$ so that we can use
!  this routine outside the frequency loop.  Thus
!  $-\frac1{w_d}\frac{\partial w_d}{\partial T} = 
!   - \frac1{2 T}$.
!  $-\frac1{w_d}\frac{\partial w_d}{\partial T}$ is what's actually
!  useful later.

    Wd = v0 * dc * sqrt(t/m)
    dWd_dT = - 0.5 * onedt ! Actually -dWd_dT/Wd

!{ $x_1 = \frac{\sqrt{\ln 2}}{w_d}$.
!  $\frac1{x_1}\frac{\partial x_1}{\partial T} =
!   -\frac1{w_d} \frac{\partial w_d}{\partial T}
!   = -\frac1{2T}$.
!  We don't calculate $x$ = $x_1 ( \nu - \nu_{0_s} )$ here because it
!  depends on frequency.  Here's $\frac1x \frac{\partial x}{\partial T} =
!  \frac1{x_1}\frac{\partial x_1}{\partial T} - \frac1{\nu - \nu_{0_s}}
!  \frac{\partial \nu_{0_s}}{\partial T}$ anyway, for reference.
!  $\frac{\partial x_1}{\partial \nu_0} = -\frac{x_1}{\nu_0}$.

    x1 = real(sqrtln2,rp) / Wd
    dx1_dT = dWd_dT ! 1/x1 dx1/dT

!{ $y = x_1 w p \left( \frac{300}T \right)^n$.
!  $\frac1y \frac{\partial y}{\partial T} =
!    \left( \frac1{x_1} \frac{\partial x_1}{\partial T} - \frac{n}T \right)
!    = -\frac1{2T}(1+2n)$.
!  $\frac{\partial y}{\partial \nu_0} = -\frac{y}{\nu_0}$.

    y = x1 * w * p * exp(n*t3t)
    dy_dT = ( dx1_dT - n * onedt ) ! 1/y dy/dT

    if ( t < 225.0_rp ) then
      q_log_b = (q(2)-q(3)) / tl1
      q_log_a = q(2) - q(1) - lt2 * q_log_b
    else
      q_log_b = (q(1)-q(2)) / tl2
      q_log_a = -lt3 * q_log_b
    end if
    q_log_b = -q_log_b - 1.0 ! This is what's interesting later

!{ $\beta_e = 10^{-4} \frac{h}k c\, e_l$. $\beta_v = \frac{h}k \nu_{0_s}$.
!  $\frac{\partial \beta_v}{\partial T} =
!    \frac{h}k \frac{\partial \nu_{0_s}}{\partial T}$.

    betae = el / boltzcm
    betav = v0s / boltzmhz ! should not be velocity corrected
    dBetav_dT = dv0s_dT / boltzmhz

!{ Write {\tt slabs1} $= \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T w_d \left(1 - e^{\frac{\beta_v}{300}}\right)}$
!  as $S = f \frac{T^{-b-1} e^{-\frac{\beta_e}T} (1+G)}
!                                    {w_d}$, where
!  $f = \frac{I_2 \, p \, e^{(i-a) \ln 10 + \frac{\beta_e}{300}}}{1-H}$
!  is independent of T,
!  $H = e^{-\frac{h \nu_0}{300 k}}$ and $G = e^{-\frac{\beta_v}T}$.  Then
!  $\frac1S \frac{\partial S}{\partial T} =
!    \frac1T \left( -b -1 -G_1 \frac{\partial \beta_v}{\partial T} +
!     \frac1T \left[ \beta_e + G_1 \beta_v \right ] \right )
!    - \frac1{w_d}\frac{\partial w_d}{\partial T}$, where
!  $G_1 = \frac{G}{1+G}$.

    expd = EXP(-v0*(oned300/boltzmhz)) ! H
    expn = EXP(-betav*onedt)   ! G
    z1 = 1.0 + expn            ! 1 + G
    z2 = 1.0 - expd            ! 1 - H

    ! This is rearranged to reduce the number of references to "exp".
    slabs1 = i2abs * p * z1 / ( wd * z2 ) * &
      & exp((i-q_log_a)*ln10 + betae*(oned300 -onedt) + q_log_b * log_t )

    z1 = expn / z1             ! G1
    ! 1/slabs1 dslabs1/dT:
    dslabs1_dT = onedt * ( q_log_b - z1 * dBetav_dT + &
      & onedt * ( betae + z1 * betav ) ) + &
      & dWd_dT ! Remember dWd_dT is really -dWd_dT/Wd

!{ $\frac{\partial S}{\partial \nu_0} =
!   -S \left[ \left( \frac{G_1}T \frac{\partial \nu_{0_s}}{\partial \nu_0}
!    + \frac{H_1}{300} \right) \frac{h}k + \frac1{\nu_0} \right] $
!   where $H_1 = \frac{H}{1-H}$.

    dslabs1_dv0 = -slabs1 * ( (z1 * onedt * velCor + expd / z2 * oned300) / boltzmhz + &
      & 1.0 / v0 )

  end subroutine Slabs_prep_dT

!{\newpage

  ! ----------------------------------------  Get_GL_Slabs_Arrays  -----
  subroutine Get_GL_Slabs_Arrays ( P_path, T_path, Vel_z, GL_Slabs, &
                             &     Do_1D, t_der_flags )

    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
    use Molecules, only: L_Extinction
    use Physics, only: SpeedOfLight
    use SpectroscopyCatalog_m, only: Catalog_T, Lines

    real(rp), intent(in) :: p_path(:) ! Pressure in hPa or mbar
    real(rp), intent(in) :: t_path(:)

    real(rp), intent(in) :: vel_z     ! Meters per second

    ! GL_Slabs needs to have been created by L2PC_pfa_structures%AllocateSlabs
    type (slabs_struct), intent(inout) :: GL_Slabs(:,:)

    logical, intent(in) :: Do_1D
    logical, intent(in), optional :: t_der_flags(:) ! do derivatives if present

!  ----------------
!  Local variables:
!  ----------------

    type(Catalog_T), pointer :: Catalog
    integer :: i, j, k, n, n_sps, nl, no_ele
    logical :: Temp_Der

    real(rp) :: vel_z_correction

! Begin code:

    no_ele = size(p_path)
    n_sps = size(gl_slabs,2)
    n = no_ele
    if ( Do_1D ) n = n / 2

    ! opposite sign convention here from ATBD
    Vel_z_correction = 1.0_rp - vel_z / speedOfLight

    do i = 1, n_sps

      catalog => gl_slabs(1,i)%catalog ! gl_slabs(:,i)%catalog are all the same
      gl_slabs(:,i)%useYi = any(lines(catalog%lines)%useYi)

      nl = Size(catalog%Lines)
      if ( nl == 0 ) cycle

      if ( catalog%molecule == l_extinction ) cycle

      do j = 1, n
        temp_der = present(t_der_flags)
        if ( temp_der ) temp_der = t_der_flags(j)
        call slabs_prep_struct ( t_path(j), p_path(j), catalog, &
          &                      Vel_z_correction, temp_der, gl_slabs(j,i) )
      end do ! j = 1, n

      if ( Do_1D ) then
        ! fill in grid points on other side with above value
        do j = no_ele, no_ele/2+1, -1
          k = no_ele - j + 1
          gl_slabs(j,i)%v0s         = gl_slabs(k,i)%v0s
          gl_slabs(j,i)%x1          = gl_slabs(k,i)%x1
          gl_slabs(j,i)%y           = gl_slabs(k,i)%y
          gl_slabs(j,i)%yi          = gl_slabs(k,i)%yi 
          gl_slabs(j,i)%slabs1      = gl_slabs(k,i)%slabs1 
          gl_slabs(j,i)%dslabs1_dv0 = gl_slabs(k,i)%dslabs1_dv0
          gl_slabs(j,i)%dx1_dv0     = gl_slabs(k,i)%dx1_dv0
          gl_slabs(j,i)%dy_dv0      = gl_slabs(k,i)%dy_dv0
        end do ! j = no_ele, no_ele/2+1, -1

        if ( present(t_der_flags) ) then
          do j = no_ele, no_ele/2+1, -1
            if ( t_der_flags(j) ) then ! do derivative stuff
              k = no_ele - j + 1
              gl_slabs(j,i)%dv0s_dT    = gl_slabs(k,i)%dv0s_dT
              gl_slabs(j,i)%dx1_dT     = gl_slabs(k,i)%dx1_dT
              gl_slabs(j,i)%dy_dT      = gl_slabs(k,i)%dy_dT
              gl_slabs(j,i)%dyi_dT     = gl_slabs(k,i)%dyi_dT
              gl_slabs(j,i)%dslabs1_dT = gl_slabs(k,i)%dslabs1_dT
            end if
          end do ! j = no_ele, no_ele/2+1, -1
        end if
      end if

    end do              ! On i = 1, n_sps

  end subroutine Get_GL_Slabs_Arrays

!=====================================================================

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module SLABS_SW_M

! $Log$
! Revision 2.41  2004/12/13 20:55:36  vsnyder
! Make Slabs_Prep and Slabs_Prep_dT public.  Add Slabs_Prep_Struct.  Polish
! some TeXnicalities.  Revise Slabswint_dT and Slabswint_Lines_dT not to
! compute interference if |yi| < 1.0e-6.  Added UseYi argument to Slabs_Prep
! and Slabs_Prep_dT.  Use Slabs_Prep_Struct from Get_GL_Slabs_Arrays.
!
! Revision 2.40  2004/09/23 20:08:47  vsnyder
! Finish correcting divide by zero
!
! Revision 2.39  2004/09/16 22:16:21  vsnyder
! Avoid dividing by zero in Slabswint_dT also
!
! Revision 2.38  2004/09/16 20:24:23  vsnyder
! Avoid dividing by zero in Slabswing_Lines_dT
!
! Revision 2.37  2004/09/01 01:14:48  vsnyder
! Correct 'not_used_here' routine
!
! Revision 2.36  2004/08/05 20:59:32  vsnyder
! Don't do any calculations for gl_slabs with no lines
!
! Revision 2.35  2004/05/11 02:52:43  vsnyder
! Remove USE for Pi, which isn't referenced
!
! Revision 2.34  2004/04/24 02:26:54  vsnyder
! Move Voigt stuff to its own module
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
