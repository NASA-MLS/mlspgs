!
module REFRACTION_M
  use MLSCommon, only: I4, R4, R8
  use L2PCDim, only: NLVL, Nptg
  use D_LINTRP_M, only: LINTRP
  use D_SOLVE_QUAD_M, only: SOLVE_QUAD
  Implicit None
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
contains

!-------------------------------------------------------------------

SUBROUTINE refractive_index(h2o_vmr,h2o_press,no_h2o_press,z_grid, &
                            t_grid,n_grid,h2o_grid,n_lvls)

! This routine computes the refractive index as a function of altitude
! The returned value has one subtracted from it

!  ===============================================================
!  Declaration of variables for sub-program: refractive_index
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls
Integer(i4), INTENT(IN) :: no_h2o_press
Real(r8), INTENT(IN) :: h2o_vmr(:), h2o_press(:), z_grid(:), t_grid(:)

Real(r8), INTENT(OUT) :: n_grid(:)
Real(r8), INTENT(IN OUT) :: h2o_grid(:)

!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: ht_i
Real(r8) :: r, p, v, ti, pv

! Derive the water function by linear interpolation of the basis coefficients

  CALL lintrp(h2o_press,z_grid,h2o_vmr,h2o_grid,no_h2o_press,n_lvls)

! Compute the relative refractive index function.
! ***************** Caution ******************************
! ADD: 1.000 TO GET THE ABSOLUTE REFRACTIVE INDEX !!!!

  pv = h2o_grid(1)
  DO ht_i = 1, n_lvls
    v = h2o_grid(ht_i)
    IF(v > 0.0) THEN
      pv = v
    ELSE
      v = pv
      h2o_grid(ht_i) = pv
    END IF
    ti = 1.0_r8 / t_grid(ht_i)
    p = 10.0_r8**(-z_grid(ht_i))
    r = 7.76D-5 * p * ti * (1.0 + 4810.0_r8 * v * ti)
    n_grid(ht_i) = r
  END DO

  RETURN
END SUBROUTINE refractive_index

!---------------------------------------------------------------

SUBROUTINE refraction_correction(n_lvls, no_conv_hts, h_grid, n_grid, &
                                 conv_hts, E_rad, ref_corr)

!  ===============================================================
!  Declaration of variables for sub-program: refraction_correction
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls, no_conv_hts

Real(r8), INTENT(IN) :: conv_hts(:), n_grid(:), h_grid(:), E_rad
Real(r8), INTENT(OUT) :: ref_corr(:,:)
!  ----------------------
!  PARAMETER Declaration:
!  ----------------------
Integer(i4), PARAMETER :: n = 4
Integer(i4), PARAMETER :: n2g = 2*n
!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, l, m, is

Real(r8), save :: xv(n2g)
Real(r8), save :: wv(n2g) = -1.0

Real(r8) :: ngrid(nlvl), hgrid(nlvl), hn, pr, hneg, hpos, q, r, hsk, hsk1, &
            an1, an2, dndh, n1, n2, h_tan, n_tan, nt, ht2, ntht2, xm, ym, &
            an, dh, h, s1, s2, h1, h2, x1, x2, eps, sum

Real(r8) :: gx(n) = &
     (/ 1.83434642495649804939D-1, 5.25532409916328985818D-1, &
        7.96666477413626739592D-1, 9.60289856497536231684D-1/)

Real(r8) :: gw(n) = &
     (/ 3.62683783378361982965D-1, 3.13706645877887287338D-1, &
        2.22381034453374470540D-1, 1.01228536290376259153D-1/)

! For first time in, Load the Gauss-Legendre abscissae & weights

  IF(wv(1) <= 0.0) THEN
    j = 0
    DO k = 1, n
      h = gx(k)
      r = gw(k)
      j = j + 1
      m = n + j
      i = n+1-j
      xv(m) = h
      wv(m) = r
      xv(i) = -h
      wv(i) = r
    END DO
  END IF

! ref_corr(1:nlvl,1:nptg) = 1.0
  DO i = 1, nptg
    ref_corr(1:nlvl,i) = 1.0
  END DO

  ngrid(1:n_lvls) = n_grid(1:n_lvls)
  DO k = 1, n_lvls
    hgrid(k) = h_grid(k) + E_rad
  END DO

  k = n_lvls
  x1 = ngrid(n_lvls-1)
  x2 = ngrid(n_lvls)
  DO WHILE(k < nlvl)
    k = k + 1
    xm = x2 * (x2 / x1)
    IF(xm >= 1.0D-13) THEN
      x1 = x2
      x2 = xm
    ELSE
      x1 = x2
      xm = x2
    END IF
    ngrid(k) = xm
  END DO

!  First approximation:

  j = 1
  is = -1
  DO i = 1, no_conv_hts
    h_tan = conv_hts(i) + E_rad
    IF(h_tan >= h_grid(1) .AND. is < 1) is = i
    ht2 = h_tan * h_tan
    DO WHILE(hgrid(j) <= h_tan .AND. j < n_lvls)
      j = j + 1
    END DO
    IF(j == 1) THEN
      l = 2
      ref_corr(1,i) = 1.0
      n_tan = ngrid(1) + ngrid(2)
      hsk1 = hgrid(2) * hgrid(2)
    ELSE
      l = j
      m = MAX(1,j-1)
      ref_corr(m,i) = 1.0
      hsk1 = hgrid(j) * hgrid(j)
      CALL exp_lintrp(hgrid(m),hgrid(j),h_tan,ngrid(m),ngrid(j),an)
      n_tan = ngrid(j) + an
    END IF
    DO k = l, n_lvls-1
      r = 0.0
      hsk = hsk1
      hsk1 = hgrid(k+1) * hgrid(k+1)
      an = n_tan - ngrid(k+2) - ngrid(k+1)
      q = 2.0_r8 * SQRT(hsk-ht2) * SQRT(hsk1-ht2)
      IF(q > 0.0) r = ht2 * an / q
      ref_corr(k,i) = 1.0_r8 + r
    END DO
  END DO

!  Modified approximation: Integral{dx/(1+n+h*dn/dh)},
!     Where: x = Sqrt((n*h)^2 - (nt*htan)^2)

  j = 1
  DO i = is, no_conv_hts

    h_tan = conv_hts(i) + E_rad

    DO WHILE(hgrid(j) <= h_tan .AND. j < n_lvls)
      j = j + 1
    END DO

    m = MAX(1,j-1)
    CALL exp_lintrp(hgrid(m),hgrid(j),h_tan,ngrid(m),ngrid(j),an)

    nt = 1.0_r8 + an
    ht2 = h_tan * h_tan
    ntht2 = nt * nt * ht2

    s2 = 0.0
    h2 = hgrid(m)
    r = h2*h2 - ht2
    IF(r > 1.0D-8) s2 = SQRT(r)

    x2 = 0.0
    an2 = ngrid(m)
    n2 = 1.0_r8 + an2
    q = n2 * h2
    r = q * q - ntht2
    IF(r > 1.0D-13) x2 = SQRT(r)

    pr = 1.0D2
    DO  k = m, n_lvls - 1

      h1 = h2
      h2 = hgrid(k+1)

      s1 = s2
      IF(h2 <= h_tan) CYCLE
      s2 = SQRT(h2*h2-ht2)
      IF(s2 <= s1) CYCLE

      an1 = an2
      an2 = ngrid(k+1)

      n1 = n2
      n2 = 1.0_r8 + an2

      x1 = x2
      q = n2 * h2
      r = q * q - ntht2
      IF(r < 0.0) CYCLE
      x2 = SQRT(r)

      dh = h2 - h1
      eps = LOG(an2/an1)/dh

      an = n2 * h2 - n1 * h1
      IF(an >= 0.0) THEN
        hpos = h2
      ELSE
        hpos = h1
      END IF

      hneg = h1 + h2 - hpos

      an = 2.0
      sum = 0.0
      xm = 0.5_r8 * (x2 + x1)
      ym = 0.5_r8 * (x2 - x1)
      DO l = 1, n2g
        q = xm + ym * xv(l)
        hn = SQRT(q*q + ntht2)
        CALL solve_hn(hn,h1,an1,eps,hneg,hpos,h,an)
        dndh = eps * an
        r = ym / (1.0_r8 + an + h * dndh)
        sum = sum + r * wv(l)
      END DO

      an = sum / (s2 - s1)
      IF(an >= 1.0_r8 .AND. an <= pr) ref_corr(k,i) = an
      pr = ref_corr(k,i)

    END DO

  END DO

  RETURN
END SUBROUTINE refraction_correction

!
!---------------------------------------------------------------------
! Exponential interpolation subroutine
!
      Subroutine Exp_Lintrp(h1,h2,h,a1,a2,a)
!
      Real*8 h1,h2,h,a1,a2,a,r,q
!
      if(abs(h2-h1).lt.1.0d-8) then
        a = 0.50 * (a1 + a2)
        return
      else if(abs(h-h1).lt.5.0d-9) then
        a = a1
        return
      else if(abs(h-h2).lt.5.0d-9) then
        a = a2
        return
      endif
!
      r = -1.0
      q = (h - h1) / (h2 - h1)
      if(abs(a1).gt.1.0d-11) r = a2 / a1
!
      if(r .gt. 0.0d0) then
!
        a = a1 * Exp(q*Log(r))
!
      else
!
! This is a sign change transition revert to linear interpolation
!
        a = a1 + (a2 - a1) * q
!
      endif
!
      Return
      End Subroutine Exp_Lintrp
!-------------------------------------------------------------------
! Solve the equation h*(1+N(h)) = Hn, where N(h) is an exponential:
!    N(h) = an1*exp(eps*(h-h1))

SUBROUTINE solve_hn(hn,h1,an1,eps,rneg,rpos,vh,vn)

REAL(r8), INTENT(IN) :: hn, h1, an1, eps, rneg, rpos
REAL(r8), INTENT(IN OUT)  :: vh
REAL(r8), INTENT(OUT)     :: vn

INTEGER(i4) :: i
REAL(r8) :: hneg, hpos, hmid, a, b, c, r1, r2, vh0, f0, f, af, df, vn0

hneg = rneg
hpos = rpos
hmid = 0.5_r8 * (hpos + hneg)

! First guess - Use first two Taylor series terms of:  exp(eps*(h-h1))
! Then solve quadratic eqation in h: (1.0+an1*exp(eps*(h-h1)))*h-Hn = 0

  a = an1 * eps
  b = 1.0_r8 + an1 * (1.0_r8 - eps * h1)
  c = -hn
  CALL solve_quad(a,b,c,r1,r2)     ! r2 <= r1, use the smaller root
  IF(r2 < MIN(hpos,hneg) .OR. r2 > MAX(hpos,hneg)) THEN
    IF(vn < 1.0) THEN
      vh0 = vh
    ELSE
      vh0 = hmid
    END IF
  ELSE
    vh0 = r2
  END IF

  vn0 = an1 * EXP(eps*(vh0-h1))
  f0 = (1.0_r8 + vn0) * vh0 - hn
  IF(f0 > 0.0) THEN
    hpos = vh0
  ELSE
    hneg = vh0
  END IF

  i = 0
  af = 1.0_r8
  DO WHILE(i < 10 .AND. af > 1.0D-7)
    i = i + 1
    df = 1.0_r8 + vn0 * (1.0_r8 + eps * vh0)
    vh = vh0 - f0 / df
    hmid = 0.5_r8 * (hpos + hneg)
    IF(vh < MIN(hpos,hneg) .OR. vh > MAX(hpos,hneg)) vh = hmid
    vn = an1 * EXP(eps*(vh-h1))
    f  = (1.0_r8 + vn) * vh - hn
    IF(f > 0.0) THEN
      hpos = vh
    ELSE
      hneg = vh
    END IF
    f0 = f
    vh0 = vh
    vn0 = vn
    af = ABS(f)
  END DO

  IF(af > 1.0D-7) THEN
    PRINT *,'** WARNING: From subroutine: SOLVE_HN **'
    PRINT *,'   Solving: h*(1+N(h)) - Hn = 0'
    PRINT *,'   Failed to converge after ',i,' iterations !!'
    PRINT *,'   Last value: h*(1+N(h))-Hn =',f
  END IF

  RETURN
END SUBROUTINE solve_hn

end module REFRACTION_M
! $Log$
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
