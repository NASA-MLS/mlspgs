module GENERIC_DELTA_INTEGRAL_M
  use MLSCOmmon, only: I4, R8
  use GL6P, only: GW, NG
  use ELLIPSE, only: HT, HT2, PS, ROC
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use PATH_ENTITIES_M, only: PATH_VECTOR
  implicit NONE
  private
  public :: GENERIC_DELTA_INTEGRAL
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
  "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes a general form for the delta Integral.
!  Using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrates in ZETA Space !
!
  Subroutine generic_delta_integral(mid, brkpt, no_ele, z_path, h_path,   &
 &           phi_path, dhdz_path, N_lvls, ref_corr, integrand, s_z_basis, &
 &           s_phi_basis, s_nz, s_np, iz, ip, fq, delta, Ier )
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, NO_ELE, S_NZ, S_NP, &
   &                           IZ, IP
!
    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH
!
    Real(r8), intent(in) :: REF_CORR(:), INTEGRAND(:), FQ

    Real(r8), intent(in) :: S_Z_BASIS(:), S_PHI_BASIS(:)

    Real(r8), intent(out) :: DELTA(:)

    Integer(i4), intent(out) :: IER

    Integer(i4) :: K, MP, NGP1, H_I, HEND

    Real(r8) :: H_GL(Ng), GW_DHDZ(Ng), INTEGRAND_GL(Ng), VETA(Ng)

    Real(r8) :: INTEGRAND_ZS, ETANP_SING, DS_DH, ETAP, ETAZ, FS, ZH, &
   &            ZL, ZS, HD, HH, HL, HTAN2, PH, PL, RC, SA, SB, SING
!
    Ier = 0
    htan2 = ht2
!
!  Initialize all arrays:
!
    integrand_GL = 0.0
!
    hend = 2 * N_lvls - 1
    delta(1:hend+1) = 0.0
!
! First, do the right hand side of the ray path:
!
    mp = 1
    ps = -1.0
    Ngp1 = Ng + 1
    zh = z_path%values(mp)
    hh = h_path%values(mp)
    ph = phi_path%values(mp)

    hd = hh + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do h_i = 1, mid
!
      mp = mp + Ngp1
      if (mp > brkpt) EXIT
!
      hl = hh
      hh = h_path%values(mp)
      if (hh < ht-0.001) EXIT
!
      zl = zh
      zh = z_path%values(mp)
!
      pl = ph
      ph = phi_path%values(mp)
!
      sa = sb
      hd = hh + RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      if (abs(sa-sb) < 0.05) EXIT
!
      rc = ref_corr(h_i+1)
!
      call define_gl_grid_entities
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
      zs = zh
      fs = ph
      etanp_sing = 0.0
      Call get_one_eta(zs,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        Call get_one_eta(fs,s_phi_basis,s_np,ip,etap)
        etanp_sing = etaz * etap
      end if
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      k = mp - Ngp1
      integrand_zs = integrand(mp)
      sing = integrand_zs * etanp_sing
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      integrand_GL(1:Ng) = integrand(k+1:k+Ng)
!
      Call gauss_legendre
!
    end do
!
! Second, do the left hand side of the ray path:
!
    ps = 1.0
    h_i = mid
    mp = brkpt + 1
    hh = h_path%values(mp)
    do while (hh < ht-0.0001)
      h_i = h_i + 1
      mp = mp + Ngp1
      hh = h_path%values(mp)
    end do
!
    zh = z_path%values(mp)
    ph = phi_path%values(mp)

    hd = hh + RoC
    sb = Sqrt(abs(hd*hd-htan2))
!
    do while (h_i < hend)
!
      mp = mp + Ngp1
      if(mp > no_ele) Return
!
      hl = hh
      hh = h_path%values(mp)
!
      zl = zh
      zh = z_path%values(mp)
!
      pl = ph
      ph = phi_path%values(mp)
!
      sa = sb
      hd = hh + RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      rc = ref_corr(h_i)
!
      Call define_gl_grid_entities
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
      zs = zl
      fs = pl
      etanp_sing = 0.0
      Call get_one_eta(zs,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        Call get_one_eta(fs,s_phi_basis,s_np,ip,etap)
        etanp_sing = etaz * etap
      end if
!
      k = mp - Ngp1
      Integrand_zs = integrand(k)
      sing = integrand_zs * etanp_sing
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      integrand_GL(1:Ng) = integrand(k+1:k+Ng)
!
      Call gauss_legendre
!
      h_i = h_i + 1
!
    end do
!
    Return
  ! *****     Internal procedures     **********************************
  contains
!
! --------------------------------     DEFINE_GL_GRID_ENTITIES     -----
    Subroutine DEFINE_GL_GRID_ENTITIES
!
    Integer :: i, j
    Real(r8) :: aym, q, v, z, phi
!
! Define the various GL grid entities for this sub-interval:
!
      aym = 0.5 * abs(zh - zl)
!
      j = mp - Ngp1
      do i = 1, Ng
        j = j + 1
        z = z_path%values(j)
        phi = phi_path%values(j)
        h_GL(i) = h_path%values(j)
        Gw_dHdZ(i) = Gw(i) * dhdz_path%values(j) * aym
        Call get_one_eta(z,s_z_basis,s_nz,iz,v)
        Call get_one_eta(phi,s_phi_basis,s_np,ip,q)
        veta(i) = v * q
      end do
    End subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    Subroutine GAUSS_LEGENDRE
!
    Integer :: i
    Real(r8) :: Sum, fv, q
!
! Now compute the Gauss-Legendre quadrature, subtructing the 'singularities
! factors':
!
      Sum = 0.0
!
      do i = 1, Ng
!
! Compute the "Hydrostatic" contribution to the derivative:
!
        hd = h_GL(i) + RoC
        ds_dh = hd / Sqrt(hd*hd-htan2)
!
        q = integrand_GL(i) * veta(i)
!
! The final integrand:
!
        Sum = Sum + (q - sing) * ds_dh * Gw_dHdz(i)
!
      end do
!
! Now add the 'singularities factors' back in, multiplied by their
! respective analytical integral:
!
      fv = Sum + sing * abs(sb-sa)
!
! And Finally - define the delta:
!
      delta(h_i) = fv * rc * fq             ! for (iz,ip)
!
    End subroutine GAUSS_LEGENDRE
  End Subroutine generic_delta_integral
End module GENERIC_DELTA_INTEGRAL_M
! $Log$
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:16  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
