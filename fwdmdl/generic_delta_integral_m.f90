! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module GENERIC_DELTA_INTEGRAL_M
  use MLSCOmmon, only: I4, R8
  use GLNP, only: GW, NG
  use ELLIPSE_M, only: ELLIPSE
  use D_GET_ONE_ETA_M, only: GET_ONE_ETA
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_INT_VECTOR_2D
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
 &           s_phi_basis,s_nz,s_np,iz,ip,fq,elvar,midval_delta,midval_ndx,&
 &           no_midval_ndx,gl_ndx,no_gl_ndx,Sps_zeta_loop,Sps_phi_loop,   &
 &           check_sps_loop,delta,Ier)
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, NO_ELE, S_NZ, S_NP, &
   &             IZ, IP, no_midval_ndx,no_gl_ndx, check_sps_loop
!
   Integer(i4), intent(in) :: midval_ndx(:,:),gl_ndx(:,:)
!
    Type(path_vector), intent(in) :: Z_PATH, H_PATH, PHI_PATH, &
   &                                 DHDZ_PATH

    Type(path_int_vector_2d), intent(in) :: Sps_phi_loop, Sps_zeta_loop
!
    Real(r8), intent(in) :: midval_delta(:)
    Real(r8), intent(in) :: REF_CORR(:), INTEGRAND(:), FQ
    Real(r8), intent(in) :: S_Z_BASIS(:), S_PHI_BASIS(:)
!
    Type(ELLIPSE), intent(in out) :: elvar

    Real(r8), intent(out) :: DELTA(:)

    Integer(i4), intent(out) :: IER

    Integer(i4) :: J, K, L, MP, NGP1, H_I, HEND, NZP

    Real(r8) :: H_GL(Ng), GW_DHDZ(Ng), INTEGRAND_GL(Ng), VETA(Ng)

    Real(r8) :: INTEGRAND_ZS, ETANP_SING, DS_DH, ETAP, ETAZ, FS, ZH, ZL, &
   &            ZS, HD, HH, HL, HTAN2, PH, PL, RC, SA, SB, SING, DS
!
    Ier = 0
!
!  Initialize the delta array:
!
    hend = 2 * N_lvls - 1
    delta(1:hend+1) = 0.0
!
! Load the mid_value delta into array: delta
!
    if(no_midval_ndx > 0) then
      do j = 1, no_midval_ndx
        h_i = midval_ndx(j,1)
        delta(h_i) = midval_delta(j)
      end do
    endif
!
    if(no_gl_ndx < 1) Return
!
! Now, do the GL deltas:
! First, do the right hand side of the ray path:
!
    k = 0
    integrand_GL = 0.0

    Ngp1 = Ng + 1
    elvar%ps = -1.0
    htan2 = elvar%ht2
    nzp = s_nz * s_np * check_sps_loop
!
    do j = 1, no_gl_ndx
!
      mp = gl_ndx(j,2)
      if (mp >= brkpt) EXIT
!
      h_i = gl_ndx(j,1)
      if(h_i >= mid) EXIT
!
      if(nzp > 1) then
        if(Sps_phi_loop%values(h_i,1) > ip  .or. &
       &   Sps_phi_loop%values(h_i,2) < ip  .or. &
       &   Sps_zeta_loop%values(h_i,1) > iz .or. &
       &   Sps_zeta_loop%values(h_i,2) < iz) CYCLE
      endif
!
      zl = z_path%values(mp)
      hl = h_path%values(mp)
      pl = phi_path%values(mp)
!
      hd = hl + elvar%RoC
      sa = Sqrt(abs(hd*hd-htan2))
!
      l = mp + Ngp1
      zh = z_path%values(l)
      hh = h_path%values(l)
      ph = phi_path%values(l)
!
      hd = hh + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
      ds = abs(sa-sb)
!
      if (ds < 0.05) EXIT
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
      integrand_zs = integrand(l)
      sing = integrand_zs * etanp_sing
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      integrand_GL(1:Ng) = integrand(mp+1:mp+Ng)
!
      k = j
      Call gauss_legendre
!
    end do

    if(k == no_gl_ndx) Return
!
    j = k
    l = gl_ndx(j,2)
    do
      j = j + 1
      if(j >= no_gl_ndx) EXIT
      mp = gl_ndx(j,2)
      if(mp - l == 1) EXIT
      l = mp
      k = j
    end do
!
! Second, do the left hand side of the ray path:
!
    do j = k+1, no_gl_ndx
!
      mp = gl_ndx(j,2)
      l = mp + Ngp1
      if (l > no_ele) EXIT
!
      h_i = gl_ndx(j,1)
      if(h_i >= hend) Return
!
      if(nzp > 1) then
        if(Sps_phi_loop%values(h_i,1) > ip  .or. &
       &   Sps_phi_loop%values(h_i,2) < ip  .or. &
       &   Sps_zeta_loop%values(h_i,1) > iz .or. &
       &   Sps_zeta_loop%values(h_i,2) < iz) CYCLE
      endif
!
      zl = z_path%values(mp)
      hl = h_path%values(mp)
      pl = phi_path%values(mp)
!
      hd = hl + elvar%RoC
      sa = Sqrt(abs(hd*hd-htan2))
!
      zh = z_path%values(l)
      hh = h_path%values(l)
      ph = phi_path%values(l)
!
      hd = hh + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
      ds = abs(sa-sb)
!
      if (ds < 0.05) EXIT
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
      integrand_zs = integrand(mp)
      sing = integrand_zs * etanp_sing
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      integrand_GL(1:Ng) = integrand(mp+1:mp+Ng)
!
      Call gauss_legendre
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
      do i = 1, Ng
        j = mp + i
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
        hd = h_GL(i) + elvar%RoC
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
      fv = Sum + sing * ds
!
! And Finally - define the delta:
!
      delta(h_i) = fv * rc * fq             ! for (iz,ip)
!
    End subroutine GAUSS_LEGENDRE
  End Subroutine generic_delta_integral
End module GENERIC_DELTA_INTEGRAL_M
! $Log$
! Revision 1.5  2001/06/07 23:30:34  pwagner
! Added Copyright statement
!
! Revision 1.4  2001/03/30 20:28:21  zvi
! General fix-up to get rid of COMMON BLOCK (ELLIPSE)
!
! Revision 1.3  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/06/21 21:56:16  zvi
! First version D.P.
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
