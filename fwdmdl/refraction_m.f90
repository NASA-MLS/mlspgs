module REFRACTION_M
  use MLSCommon, only: I4, R4, R8
  use L2PCDIM, only: NLVL, N2lvl, Nptg
  use D_LINTRP_M, only: LINTRP
  use D_SOLVE_QUAD_M, only: SOLVE_QUAD
  use TWO_D_POLATE_M, only: TWO_D_POLATE
  use PATH_ENTITIES_M, only: PATH_INDEX, PATH_VECTOR
  use GL6P, only: NG, GX, GW
  Implicit None
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
 "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
 "$RCSfile$"
!---------------------------------------------------------------------------
contains

!-------------------------------------------------------------------

SUBROUTINE refractive_index(h2o_vmr,h2o_press,h2o_phi,no_h2o_press, &
           no_h2o_phi,ndx_path,z_path,t_path,phi_path,n_path,       &
           wet,no_tan_hts)

! This routine computes the refractive index as a function of altitude
! and phi. The returned value has one subtracted from it

!  ===============================================================
!  Declaration of variables for sub-program: refractive_index
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Logical, INTENT(IN) :: wet
Integer(i4), INTENT(IN) :: no_h2o_press,no_h2o_phi,no_tan_hts

Real(r8), INTENT(IN) :: h2o_vmr(:,:),h2o_press(:),h2o_phi(:)

Type(path_index),  INTENT(IN)  :: ndx_path(:)
Type(path_vector),  INTENT(IN) :: z_path(:), t_path(:), phi_path(:)

Type(path_vector), INTENT(OUT) :: n_path(:)

!  ----------------
!  Local variables:
!  ----------------
Integer(i4) :: ht_i, ptg_i, j, ier
Real(r8) :: r, p, v, z, ti, phi

! Derive the water function by linear interpolation of the basis coefficients
! and compute the relative refractive index function.
! ***************** Caution ******************************
! ADD: 1.000 TO GET THE ABSOLUTE REFRACTIVE INDEX !!!!

  ptg_i = no_tan_hts
  DEALLOCATE(n_path(ptg_i)%values,STAT=j)

  DO ptg_i = 1, no_tan_hts
    j = ndx_path(ptg_i)%total_number_of_elements
    ALLOCATE(n_path(ptg_i)%values(j), STAT=ier)
    if(Ier /= 0) then
      Print *,'** Allocation error in "refractive_index" routine !'
      Stop
    endif
    n_path(ptg_i)%values(1:j) = 0.0
    DO ht_i = 1, j
      z = z_path(ptg_i)%values(ht_i)
      if(z < -5.0) CYCLE
      ti = 1.0_r8 / t_path(ptg_i)%values(ht_i)
      p = 10.0_r8**(-z)
      r = 7.76D-5 * p * ti
      if(wet) then
        phi = phi_path(ptg_i)%values(ht_i)
        Call TWO_D_POLATE (h2o_press,h2o_vmr,no_h2o_press,h2o_phi, &
                           no_h2o_phi,z,phi,v)
        r = r * (1.0_r8 + 4810.0_r8 * v * ti)
      endif
      n_path(ptg_i)%values(ht_i) = r
    END DO
  END DO

  RETURN
END SUBROUTINE refractive_index

!---------------------------------------------------------------

SUBROUTINE refraction_correction(no_tan_hts, tan_hts, h_path, n_path, &
                      ndx_path, E_rad, ref_corr)

!  ===============================================================
!  Declaration of variables for sub-program: refraction_correction
!  ===============================================================

!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: no_tan_hts

Real(r8), INTENT(IN)  :: tan_hts(:), E_rad

Type(path_index),  INTENT(IN) :: ndx_path(:)
Type(path_vector), INTENT(IN) :: h_path(:), n_path(:)

Real(r8), INTENT(OUT) :: ref_corr(:,:)       ! (N2lvl,Nptg)

!  ----------------------
!  Local variables:
!  ----------------
Integer(i4) :: i, j, k, m, ptg_i, no_ele, Ngp1, count, revc, brkpt

Real(r8) :: ngrid(Nlvl), hgrid(Nlvl), refcor(Nlvl), h_tan, ht2, h

  Ngp1 = Ng + 1

  DO ptg_i = 1, no_tan_hts-1

    brkpt = ndx_path(ptg_i)%break_point_index
    no_ele = ndx_path(ptg_i)%total_number_of_elements

    h_tan = tan_hts(ptg_i) + E_rad
    ht2 = h_tan * h_tan
!
    hgrid(1:Nlvl) = -1.0
    ngrid(1:Nlvl) =  1.0
    ref_corr(1:N2lvl,ptg_i) = 1.0

    k = 0
    j = 1 - Ngp1
    revc = (brkpt + Ng) / Ngp1
    do i = 1, revc
      j = j + Ngp1
      h = h_path(ptg_i)%values(j) + E_rad
      if(h >= h_tan) k = k + 1
    end do

    count = 0
    j = 1 - Ngp1
    do i = 1, revc
      j = j + Ngp1
      h = h_path(ptg_i)%values(j) + E_rad
      if(h >= h_tan) then
        count = count + 1
        hgrid(k+1-count) = h
        ngrid(k+1-count) = n_path(ptg_i)%values(j)
      endif
    end do

    m = count
    CALL create_ref_corr(count, h_tan, ht2, hgrid, ngrid, refcor)
    do i = 1, count
      ref_corr(i,ptg_i) = refcor(count+1-i)
    end do

    hgrid(1:Nlvl) = -1.0
    ngrid(1:Nlvl) =  1.0

    count = 0
    j = brkpt+1-Ngp1
    revc = (no_ele-brkpt+Ng) / Ngp1
    do i = 1, revc
      j = j + Ngp1
      h = h_path(ptg_i)%values(j) + E_rad
      if(h >= h_tan) then
        count = count + 1
        hgrid(count) = h
        ngrid(count) = n_path(ptg_i)%values(j)
      endif
    end do

    CALL create_ref_corr(count, h_tan, ht2, hgrid, ngrid, refcor)
    do i = 1, count
      m = m + 1
      ref_corr(m,ptg_i) = refcor(i)
    end do
!
  END DO

  k = no_tan_hts
  ref_corr(1:,k) = ref_corr(1:,k-1)

  RETURN
END SUBROUTINE refraction_correction

!---------------------------------------------------------------

SUBROUTINE create_ref_corr(n_lvls, ht, ht2, hgrid, ngrid, refcor)

!  ===============================================================
!  Declaration of variables for sub-program: refraction_correction
!  ===============================================================
!  ---------------------------
!  Calling sequence variables:
!  ---------------------------
Integer(i4), INTENT(IN) :: n_lvls

Real(r8), INTENT(IN) :: ht, ht2
Real(r8), INTENT(INOUT) :: ngrid(:), hgrid(:)
Real(r8), INTENT(OUT) :: refcor(:)
!  ----------------------
!  Local variables:
!  ----------------
Integer(i4) :: j, k, l, m

Real(r8) :: hn, pr, hneg, hpos, q, r, hsk, hsk1, an1, an2, dndh, n1, n2, &
            ngrd_tan, nt, ntht2, xm, ym, an, dh, h, s1, s2, h1, h2, x1,  &
            x2, eps, sum

! For first time in, Load the Gauss-Legendre abscissae & weights

  refcor(1:nlvl) = 1.0

  k = n_lvls
  x1 = ngrid(n_lvls-1)
  x2 = ngrid(n_lvls)
  DO WHILE(k < nlvl .AND. x1 > 1.0d-13)
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
  DO WHILE(hgrid(j) <= ht .AND. j < n_lvls)
    j = j + 1
  END DO

  IF(j == 1) THEN
    l = 2
    refcor(1) = 1.0
    ngrd_tan = ngrid(1) + ngrid(2)
    hsk1 = hgrid(2) * hgrid(2)
  ELSE
    l = j
    m = MAX(1,j-1)
    refcor(m) = 1.0
    hsk1 = hgrid(j) * hgrid(j)
    CALL exp_lintrp(hgrid(m),hgrid(j),ht,ngrid(m),ngrid(j),an)
    ngrd_tan = ngrid(j) + an
  END IF
  DO k = l, n_lvls-1
    r = 0.0
    hsk = hsk1
    hsk1 = hgrid(k+1) * hgrid(k+1)
    an = ngrd_tan - ngrid(k+2) - ngrid(k+1)
    q = 2.0_r8 * SQRT(hsk-ht2) * SQRT(hsk1-ht2)
    IF(q > 0.0) r = ht2 * an / q
    refcor(k) = 1.0_r8 + r
  END DO

  q = refcor(n_lvls-1)
  refcor(n_lvls:nlvl) = q

  if(ht < 0.0) Return

! For derivation of the code below, please see: "FWD Model" paper,
! Page 16, Eqn. 26 & 27

!  Modified approximation: Integral{dx/(1+n+h*dn/dh)},
!     Where: x = Sqrt((n*h)^2 - (nt*htan)^2)
!
  j = 1
  DO WHILE(hgrid(j) <= ht .AND. j < n_lvls)
    j = j + 1
  END DO

  m = MAX(1,j-1)
  CALL exp_lintrp(hgrid(m),hgrid(j),ht,ngrid(m),ngrid(j),an)

  nt = 1.0_r8 + an
  ntht2 = nt * nt * ht2

  s2 = 0.0
  h2 = hgrid(m)
  r = h2*h2 - ht2
  IF(r > 1.0D-8) s2 = SQRT(r)

  x2 = 0.0_r8
  an2 = ngrid(m)
  n2 = 1.0_r8 + an2
  q = n2 * h2
  r = q * q - ntht2
  IF(r > 1.0D-13) x2 = SQRT(r)

  pr = 100.0
  DO  k = m, n_lvls - 1

    h1 = h2
    h2 = hgrid(k+1)

    s1 = s2
    IF(h2 <= ht) CYCLE
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
    sum = 0.0_r8
    xm = 0.5_r8 * (x2 + x1)
    ym = 0.5_r8 * (x2 - x1)
    DO l = 1, ng
      q = xm + ym * Gx(l)
      hn = SQRT(q*q + ntht2)
      CALL solve_hn(hn,h1,an1,eps,hneg,hpos,h,an)
      dndh = eps * an
      r = ym / (1.0_r8 + an + h * dndh)
      sum = sum + r * Gw(l)
    END DO

    an = sum / (s2 - s1)
    IF(an >= 1.0_r8 .AND. an <= pr) refcor(k) = an
    pr = refcor(k)

  END DO

  q = refcor(n_lvls-1)
  refcor(n_lvls:nlvl) = q

  RETURN
END SUBROUTINE create_ref_corr

!---------------------------------------------------------------------
! Exponential interpolation subroutine
!
      Subroutine Exp_Lintrp(h1,h2,h,a1,a2,a)
!
      Real(r8) :: h1,h2,h,a1,a2,a,r,q
!
      if(abs(h2-h1).lt.1.0e-8) then
        a = 0.50 * (a1 + a2)
        return
      else if(abs(h-h1).lt.5.0e-9) then
        a = a1
        return
      else if(abs(h-h2).lt.5.0e-9) then
        a = a2
        return
      endif
!
      r = -1.0
      q = (h - h1) / (h2 - h1)
      if(abs(a1).gt.1.0e-11) r = a2 / a1
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
!    N(h) = an1*Exp(eps*(h-h1))

SUBROUTINE solve_hn(hn,h1,an1,eps,rneg,rpos,vh,vn)

REAL(r8), INTENT(IN) :: hn, h1, an1, eps, rneg, rpos
REAL(r8), INTENT(IN OUT)  :: vh
REAL(r8), INTENT(OUT)     :: vn

INTEGER(i4) :: i
REAL(r8) :: hneg, hpos, hmid, a, b, c, r1, r2, vh0, f0, f, af, df, vn0

hneg = rneg
hpos = rpos
hmid = 0.5_r8 * (hpos + hneg)

! First guess - Use first two Taylor series terms of:  Exp(eps*(h-h1))
! Then solve quadratic eqation in h: (1.0+an1*Exp(eps*(h-h1)))*h-Hn = 0

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
