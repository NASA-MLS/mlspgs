module SLABS_SW_M

  use MLSCommon, only: R8, RP, IP
  use SpectroscopyCatalog_m, only: CATALOG_T, Lines
  use Units, only: Pi, SqrtPi

  implicit NONE

  private
  public :: DVoigt_Spectral, Slabs, Slabswint, Voigt_Lorentz, &
        &  Real_Simple_Voigt, Simple_Voigt, RLorentz, CLorentz, RVoigth2, &
        &  CVoigth2, RVoigth6, CVoigth6, RHui6, CHui6, RDrayson, &
        &  CDrayson, Slabs_Prep, Slabs_Prep_Arrays, Get_GL_Slabs_Arrays

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
!       tanh1 = tanh(nu*expa / 2.0)

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

  ! ------------------------------------------------------  Slabs  -----
  real(rp) function Slabs ( Nu, v0, v0s, x1, tanh1, slabs1, y )

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

    call real_simple_voigt(x1*real(nu-v0s,rp),y,u)

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, and
!  $D = \frac1{\sigma^2 x_1^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh(\frac{h \nu}{2 k T}) \left( V(a,y) + \frac{y D}{\sqrt{\pi}} \right)$.

    Slabs = slabs1 * real(nu / v0, rp) * tanh1 * &
      & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

  end function Slabs

  ! ----------------------------------------------------  Slabs_dT  -----
  subroutine Slabs_dT ( Nu, v0, v0s, x1, tanh1, slabs1, y, &
    &                           dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, dy_dT, &
    &                   Slabs, dSlabs_dT )

  ! Compute single-line absorption and its derivative w.r.t. temperature.

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t ) (tanh(...) - 1/tanh(...) )
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
    call D_Real_Simple_Voigt ( x1*delta, y, da, dy_dT, u, du )

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, and
!  $D = \frac1{\sigma^2 x_1^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh(\frac{h \nu}{2 k T}) \left( V(a,y) + \frac{y D}{\sqrt{\pi}} \right)$.

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
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) V(a,y)$ and
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) \frac{y D}{\sqrt{\pi}}$.\\
!  Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$ and\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
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

  ! --------------------------------------------------  Slabswint  -----
  real(rp) function Slabswint ( Nu, v0, v0s, x1, tanh1, slabs1, y, yi )

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
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T})
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

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y, yi
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t ) (tanh(...) - 1/tanh(...) )
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
    call D_Real_Simple_Voigt ( a, y, da, dy_dT, u, du )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$,
!  $D_1 = \frac1{\sigma^2 x_1^2 + y^2}$ and $D_2 = \frac1{a^2 + y^2}$.
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T})
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
    sd = c * a * yi * d2
    Slabswint = sa + sb - sc + sd

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as $V(a,y) +
!  i L(a,y)$, where $V(a,y)$ is the Voigt function ({\tt u} above) and $L(a,y)$
!  is the Lorentz function.  All we want for $S$ is its real part, so we don't
!  need $L(a,y)$ for $S$. For $\frac{\partial S}{\partial T}$ we want the real
!  part of the derivative; this requires the real part of the derivative of
!  $w(z)$, not just Voigt.
!
!  Write $S = S_a + S_b - S_c + S_d$ where
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) V(a,y)$,
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) \frac{y D_1}{\sqrt{\pi}}$,
!  $S_c = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) \frac{\sigma x_1 y_i D_1}{\sqrt{\pi}}$, and
!  $S_d = S_1 \frac{\nu}{\nu_0} \tanh(\frac{h \nu}{2 k T}) \frac{a y_i D_2}{\sqrt{\pi}}$.\\
!  Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T} - \frac{\partial S_c}{\partial T} +
!   \frac{\partial S_d}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$,\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
!   \frac1{D_1}\frac{\partial D_1}{\partial T} +
!   \frac1y \frac{\partial y}{\partial T}$,\\
!  $\frac1{S_c}\frac{\partial S_c}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
!   \frac1{D_1}\frac{\partial D_1}{\partial T} +
!   \frac1{y_i}\frac{\partial y_i}{\partial T} +
!   \frac1{x_1}\frac{\partial x_1}{\partial T} +
!   \frac1{\sigma}\frac{\partial \nu_{0_s}}{\partial T}$, and\\
!  $\frac1{S_d}\frac{\partial S_d}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh(\frac{h \nu}{2 k T})}\frac{\partial}{\partial T} \tanh(\frac{h \nu}{2 k T}) +
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
    dSlabs_dT = sa * ( c + du / u ) + &
      &         sb * ( c + dd1 + dy_dT ) - &
      &         sc * ( c + dd1 + dyi_dT + dx1_dT + dv0s_dT / sigma ) + &
      &         sd * ( c + dd2 + dyi_dT + da / a )

  end subroutine Slabswint_dT

  ! ----------------------------------------------  Voigt_Lorentz  -----

  subroutine Voigt_Lorentz ( dNu,  Nu0,  x1,  yi,  y,  w,  t,  tanh1, slabs1,  &
                         &   VL, dslabs1_dNu0,  dVL_dw,  dVL_dn,  dVL_dNu0 )

! Compute the Voigt/Lorentz function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

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
    real :: DELY, DX, F, FD, IS, IT, Q(5), R(5), RS, RT, W, XA, Y2

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

    ! For 2-pt Gauss-Hermite
    real(rp), parameter :: gx = 0.70710678118655_rp ! 1.0/sqrt(2.0)
    real(rp), parameter :: gw = 0.88622692545277_rp / Pi

    ! For 6-pt Gauss-Hermite
    real(rp), parameter :: gx6(3) = (/ 4.36077411927616508271e-1_rp, &
                                       1.33584907401369694957_rp,    &
                                       2.35060497367449222280_rp /)
    real(rp), parameter :: gw6(3) = (/ 7.24629595224392524608e-1_rp, &
                                       1.57067320322856644036e-1_rp, &
                                       4.53000990550884564224e-3_rp /) / Pi

    real(rp), parameter :: TwoThirds = 2.0_rp / 3.0_rp
    real(rp), parameter :: XL=5.2_rp, YL=0.05_rp, YH=0.6_rp, DYDX=(yh-yl)/xl

! This is sorted in likely occurance of each case

    xa = ABS(x)

! I am assuming that the OR are evaluated sequentially until the first
! true is found. Also routines are ordered according to speed

    if ( y + TwoThirds*xa > 100.0_rp ) then

      ! Here x is sqrt(ln2)*delnu / wd and y is sqrt(ln2)*wc / wd
      ! u = rlorentz ( x, y )
      u = OneOvSPi * y / (y*y + x*x)

    else if ( y + 0.6875_rp * xa > 11.0_rp ) then

! Drayson's quick 2pt hermite integral (essentially a lorentz)

      ! u = rvoigth2(xa,y)
      y2 = y**2
      u = gw * y * (1.0_rp/(y2 + (xa-gx)**2) + &
        &           1.0_rp/(y2 + (xa+gx)**2))

    else if ( y > 0.6_rp .OR. y > yl + dydx*xa ) then

! Intermediate region

      ! u = rhui6(xa,y)

      !   Fill the r, q coefficients with the recursion
      !   r(0)=1.0, r(1) = y q(0) = 0.0, q(1) = -x

      r(1) = y**2 - xa**2
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

      y2 = y**2

      u = y * ( gw6(1) * (1.0_rp/(y2 + (xa-gx6(1))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(1))**2)) + &
                gw6(2) * (1.0_rp/(y2 + (xa-gx6(2))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(2))**2)) + &
                gw6(3) * (1.0_rp/(y2 + (xa-gx6(3))**2) + &
                     &    1.0_rp/(y2 + (xa+gx6(3))**2)) )
    else

! Near line center where Doppler dominates Pressure broadening

      ! u = rdrayson(xa,y)

      !******** Region I. Compute Dawson's function at x from Taylor series

      if ( y <= 1.0e-12_rp ) then
        u = exp(-xa*xa)
      else

        !  Taylor series expansion about y = 0.0

        y2 = y*y
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

    ! For 2-pt Gauss-Hermite
    real(rp), parameter :: GX = 0.70710678118655_rp ! 1.0/sqrt(2.0)
    real(rp), parameter :: GW = 0.88622692545277_rp / Pi

    real(rp), parameter :: XL=5.2_rp, YL=0.05_rp, YH=0.6_rp, DYDX=(yh-yl)/xl
    real(rp) :: DENOM1, DENOM2, XA, XM, XP, Y2
    complex(rp) :: UV

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

      ! uv = cvoigth2(xa,y)
      xm = xa - gx
      xp = xa + gx
      y2 = y**2
      denom1 = gw/(y2 + xm**2)
      denom2 = gw/(y2 + xp**2)
      uv = CMPLX(y*(denom1+denom2), xm*denom1 + xp*denom2, KIND=rp)

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

    u = real(uv,kind=rp)
    if ( present(v) ) v = sign(aimag(uv),x)

  end subroutine Simple_Voigt

  ! ---------------------------------------------------  RLorentz  -----
    real(rp) pure function RLorentz ( x, y )

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

!{ \newpage

  ! ----------------------------------------  D_Real_Simple_Voigt  -----
  subroutine D_Real_Simple_Voigt ( x, y, dx, dy, u, du )
  ! Compute the real part of Fadeeva = Voigt (without Lorentz) and
  ! the real part of the derivative (which isn't just the derivative
  ! of Voigt).

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as
!  $V(a,y) + i L(a,y)$, where $V(a,y)$ is the Voigt function and
!  $L(a,y)$ is the Lorentz function. 
!
!  From 7.1.20 in {\bf Handbook of Mathematical Functions} by Abramowitz
!  and Stegun (National Bureau of Standards Applied Math Series 55) we have
!  $w^{\prime}(z) = \frac{2i}{\sqrt{\pi}} - 2 z w(z)$.\\
!  $\Re \, \frac{\partial w(z(t))}{\partial t} =
!   2 \left[ \left( -x V(x,y) + y L(x,y) \right) \frac{\partial x}{\partial t} +
!            \left( x L(x,y) + y V(x,y) -\frac1{\sqrt{\pi}} \right)
!             \frac{\partial y}{\partial t} \right]$.

    real(rp), intent(in) :: x, y, dx, dy
    real(rp), intent(out) :: u, du

    real(rp) :: v ! Lorentz function

    call simple_voigt ( x, y, u, v )
    du = 2.0_rp * ( (-x * u + y * v) * dx + (x * v + y * u) * dy )

  end subroutine D_Real_Simple_Voigt

  ! ---------------------------------------------  D_Simple_Voigt  -----
  subroutine D_Simple_Voigt ( x, y, dx, dy, u, v, du, dv )
  ! Compute Fadeeva = Voigt and Lorentz its derivative

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as
!  $V(a,y) + i L(a,y)$, where $V(a,y)$ is the Voigt function and
!  $L(a,y)$ is the Lorentz function. 
!
!  From 7.1.20 in {\bf Handbook of Mathematical Functions} by Abramowitz
!  and Stegun (National Bureau of Standards Applied Math Series 55) we have
!  $w^{\prime}(z) = \frac{2i}{\sqrt{\pi}} - 2 z w(z)$.\\
!  \begin{equation*}\begin{split}
!  \frac{\partial w(z(t))}{\partial t}
!   =& 2 \left[ \left( -x V(x,y) + y L(x,y) \right) \frac{\partial x}{\partial t} +
!            \left( x L(x,y) + y V(x,y) -\frac1{\sqrt{\pi}} \right)
!             \frac{\partial y}{\partial t} \right]\\
!   +& 2 i \left [ - \left( x L(x,y) + y V(x,y) -\frac1{\sqrt{\pi}} \right)
!             \frac{\partial x}{\partial t} +
!            \left( - x V(x,y) + y L(x,y) \right) \frac{\partial y}{\partial t} \right]\\
!   =& a \frac{\partial x}{\partial t} + b \frac{\partial y}{\partial t}
!    + i \left( -b \frac{\partial x}{\partial t} + a \frac{\partial y}{\partial t} \right).
!  \end{split}\end{equation*}

    real(rp), intent(in) :: x, y, dx, dy
    real(rp), intent(out) :: u, v, du, dv

    real :: a, b

    call simple_voigt ( x, y, u, v )
    a = -x * u + y * v
    b = x * v + y * u
    du = 2.0_rp * ( a * dx + b * dy )
    dv = 2.0_rp * ( -b * dx + a * dy )

  end subroutine D_Simple_Voigt

  ! -------------------------------------------------  Slabs_prep  -----
  subroutine Slabs_prep ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                      &   n1, n2, &
                      &   v0s, x1, y, yi, slabs1, dslabs1 )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! >>2004-03-18 WV Snyder Use 1 - exp(-v0/300/Boltzmhz) in denominator instead
!                        of 1 - exp(-v0s/300/Boltzmhz).

    use Physics, only: H_OVER_K, SpeedOfLight
    use Units, only: Ln10

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

    real(rp), parameter :: I2abs = 3.402136078e9_rp ! sqrt(ln(2)/pi)*1.0e-13/k
    real(rp), parameter :: Dc = 3.58117369e-7_rp ! sqrt(1000 k ln 4 avogadro) / c
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Sqrtln2 = 8.32554611157698e-1_rp
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

    yi = p * (delta*exp(n1*t3t) + gamma*exp(n2*t3t))

!{ $\nu_{0_s} = \nu_0 + p_s p \left( \frac{300}T \right) ^{n_s}$.

    v0s = v0 + ps * p * exp(ns*t3t)

!{ $\omega_d = \nu_0 d_c \sqrt{\frac{T}M}$.

    Wd = v0 * Sqrt(t/m) * dc

!{ $x_1 = \frac{\sqrt{\ln 2}}{\omega_d}$.

    x1 = sqrtln2 / Wd

!{ $y = x_1 w p \left( \frac{300}T \right) ^n$.

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

    expn = EXP(-betav*onedt)
    z1 = 1.0 + expn
    expd = EXP(-v0*(oned300/boltzmhz))
    z2 = 1.0 - expd

!{ $S = \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T \omega_d \left(1 - e^{-\frac{h \nu_0}{300 k}}\right)} =
!  I_2 \, p \, e^{(i-a) \ln 10 -(b+1) \log T +
!    \beta_e \left( \frac1{300} - \frac1T \right)}
!  \frac{(1+G)} {\omega_d (1-H)}$, where
!  $G = e^{-\frac{\beta_v}T}$, $H = e^{-\frac{h \nu_0}{300 k}}$ and
!  $S =$ {\tt slabs1}.

    slabs1 = i2abs * p * &
      & exp((i - q_log_a)*ln10 + q_log_b * log_T + betae * (oned300 - onedt)) &
      & * z1 / (Wd * z2)

!{ $\frac{\partial S}{\partial \nu_0} =
!   -S \left( \frac{G_1}T + \frac{H_1}{300} \right) \frac{h}k$ where
!   $G_1 = 1 + G$ and $H_1 = 1 - H$.

    dslabs1 = -slabs1 * (expn / (t * z1) + expd / (300.0_rp * z2)) / &
      & boltzmhz

  end subroutine Slabs_prep

  ! ----------------------------------------------  Slabs_prep_DT  -----
  subroutine Slabs_prep_dT ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                         &   n1, n2, &
                         &   v0s, x1, y, yi, slabs1, dslabs1_dv0, &
                         &   dv0s_dT, dx1_dT, dy_dT, dyi_dT, dslabs1_dT )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
! Compute the derivatives with respect to temperature, too.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects

    use Physics, only: H_OVER_K, SpeedOfLight
    use Units, only: Ln10

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
    real(rp), intent(out) :: Dslabs1_dv0 ! Derivative of slabs1 w.r.t. v0

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Derivatives with respect to temperature:

    real(r8), intent(out) :: dv0s_dT    ! dv0s/dT
    real(rp), intent(out) :: dx1_dT     ! 1/x1 dx1/dT
    real(rp), intent(out) :: dy_dT      ! 1/y dy/dT
    real(rp), intent(out) :: dyi_dT     ! 1/yi dyi/dT
    real(rp), intent(out) :: dslabs1_dT ! 1/slabs1 dslabs1/dT

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K = k/h
!  sqrtln2  - sqrt(ln(2))

    real(rp), parameter :: I2abs = 3.402136078e9_rp ! sqrt(ln(2)/pi)*1.0e-13/k
    real(rp), parameter :: Dc = 3.58117369e-7_rp
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Sqrtln2 = 8.32554611157698e-1_rp
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

    z1 = delta*exp(n1*t3t)
    z2 = gamma*exp(n2*t3t)
    yi = ( z1 + z2 )
    dyi_dT = -onedt * ( n1 * z1 + n2 * z2 ) / yi ! 1/yi dyi/dT
    yi = p * yi

!{ $\nu_{0_s} = \nu_0 + p_s p \left( \frac{300}T \right)^{n_s}$.
!  $\frac{\partial \nu_{0_s}}{\partial T} = \frac{-n_s}T ( \nu_{0_s} - \nu_0 )$.

    v0s = ps * p * exp(ns*t3t)
    dv0s_dT = -ns * v0s * onedt
    v0s = v0 + v0s

!{ $\omega_d = \nu_0 d_c \sqrt{\frac{T}M}$.  The $\nu_0$ term should
!  really be $\nu$.  We approximate $\nu$ by $\nu_0$ so that we can use
!  this routine outside the frequency loop.  Thus
!  $-\frac1{\omega_d}\frac{\partial \omega_d}{\partial T} = 
!   - \frac1{2 T}$.
!  $-\frac1{\omega_d}\frac{\partial \omega_d}{\partial T}$ is what's actually
!  useful later.

    Wd = v0 * dc * sqrt(t/m)
    dWd_dT = - 0.5 * onedt ! Actually -dWd_dT/Wd

!{ $x_1 = \frac{\sqrt{\ln 2}}{\omega_d}$.
!  $\frac1{x_1}\frac{\partial x_1}{\partial T} =
!   -\frac1{\omega_d} \frac{\partial \omega_d}{\partial T}
!   = -\frac1{2T}$.

    x1 = sqrtln2 / Wd
    dx1_dT = dWd_dT ! 1/x1 dx1/dT

!{ $y = x_1 w p \left( \frac{300}T \right)^n$.
!  $\frac1y \frac{\partial y}{\partial T} =
!    \left( \frac1{x_1} \frac{\partial x_1}{\partial T} - \frac{n}T \right)
!    = -\frac1{2T}(1+2n)$.

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
    betav = v0s / boltzmhz
    dBetav_dT = dv0s_dT / boltzmhz

!{ Write {\tt slabs1} $= \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T \omega_d \left(1 - e^{\frac{\beta_v}{300}}\right)}$
!  as $S = f \frac{T^{-b-1} e^{-\frac{\beta_e}T} (1+G)}
!                                    {\omega_d (1-H)}$, where
!  $f = I_2 \, p \, e^{(i-a) \ln 10 + \frac{\beta_e}{300}}$,
!  $G = e^{-\frac{\beta_v}T}$ and $H = e^{-\frac{\beta_v}{300}}$.  Then
!  $\frac1S \frac{\partial S}{\partial T} =
!    \frac1T \left( -b -1 -G_1 \frac{\partial \beta_v}{\partial T} +
!     \frac1T \left[ \beta_e + G_1 \beta_v \right ] \right )
!    - \frac1{\omega_d}\frac{\partial \omega_d}{\partial T}
!    - \frac{H_1}{300} \frac{\partial \beta_v}{\partial T}$, where
!  $G_1 = \frac{G}{1+G}$ and $H_1 = \frac{H}{1-H}$.

    expd = EXP(-betav*oned300) ! H
    expn = EXP(-betav*onedt)   ! G
    z1 = 1.0 + expn            ! 1 + G
    z2 = 1.0 - expd            ! 1 - H

    ! This is rearranged to reduce the number of references to "exp".
    slabs1 = i2abs * p * z1 / ( wd * z2 ) * &
      & exp((i-q_log_a)*ln10 + betae*(oned300 -onedt) + q_log_b * log_t )

    z1 = expn / z1             ! G1
    z2 = oned300 * expd / z2   ! H1 / 300
    ! 1/slabs1 dslabs1/dT:
    dslabs1_dT = onedt * ( q_log_b - z1 * dBetav_dT + &
      & onedt * ( betae + z1 * betav ) ) + &
      & dWd_dT - z2 * dbetav_dT ! Remember dWd_dT is really -dWd_dT/Wd

!{ $\frac{\partial S}{\partial \nu_0} =
!   -S \left( \frac{G_1}T + \frac{H_1}{300} \right) \frac{h}k$.

    dslabs1_dv0 = -slabs1 * (z1 * onedt + z2) / boltzmhz

  end subroutine Slabs_prep_dT

!{ \newpage

  ! -----------------------------------------  Slabs_Prep_Arrays   -----
  subroutine Slabs_Prep_Arrays ( molecule, nl, t, p, mass, Qlog, Catalog, &
                               & v0s, x1, y, yi, slabs1, dslabs1_dv0 )

    use Molecules, only: L_Extinction

    type(catalog_T) :: Catalog

    integer(ip), intent(in) :: molecule, nl

    real(rp), intent(in) :: t, p, mass,Qlog(:)

    real(r8), intent(out) :: v0s(:)
    real(rp), intent(out) :: x1(:),y(:),yi(:),slabs1(:),dslabs1_dv0(:)

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
                             &     Do_1D, t_der_flags )

    use Physics, only: SpeedOfLight
    use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT

    type(Catalog_T), dimension(:), intent(in) :: Catalog

    real(rp), intent(in) :: p_path(:) ! Pressure in hPa or mbar
    real(rp), intent(in) :: t_path(:)

    real(rp), intent(in) :: vel_z
    logical, intent(in) :: Do_1D
    logical, intent(in) :: t_der_flags(:)

    type (slabs_struct) :: gl_slabs(:,:)

!  ----------------
!  Local variables:
!  ----------------

    real(rp), parameter :: c = speedOfLight/1000.0_rp ! Speed of Light Km./Sec.

    integer :: i, j, k, n_sps, nl, no_ele

    real(rp) :: Qlog(3), vel_z_correction

! Begin code:

    no_ele = size(p_path)
    n_sps = Size(catalog)

    Vel_z_correction = 1.0_rp - vel_z / c

    do i = 1, n_sps

      nl = Size(Catalog(i)%Lines)

      Qlog(1:3) = Catalog(i)%QLOG(1:3)

      if ( .not. Do_1D ) then

        do j = 1, no_ele
          if (.not. t_der_flags(j)) cycle
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
          if (.not. t_der_flags(j)) cycle
          call Slabs_Prep_Arrays ( catalog(i)%molecule, nl, t_path(j), p_path(j), &
            & catalog(i)%mass, Qlog, &
            & Catalog(i), gl_slabs(j,i)%v0s, gl_slabs(j,i)%x1, gl_slabs(j,i)%y, &
            & gl_slabs(j,i)%yi,gl_slabs(j,i)%slabs1,gl_slabs(j,i)%dslabs1_dv0 )

          gl_slabs(j,i)%v0s = gl_slabs(j,i)%v0s * Vel_z_correction

        end do
        
        ! fill in grid points on other side with above value
        
        do j = no_ele, no_ele/2+1, -1
          if (.not. t_der_flags(j)) cycle          
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
