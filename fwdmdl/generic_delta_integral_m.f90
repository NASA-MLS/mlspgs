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
 &           delta,Ier)
!
    Integer(i4), intent(in) :: N_LVLS, MID, BRKPT, NO_ELE, S_NZ, S_NP, &
   &             IZ, IP, no_midval_ndx,no_gl_ndx
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

    Integer(i4) :: J, K, L, MP, NGP1, H_I, HEND

    Real(r8) :: H_GL(Ng),GW_DHDZ(Ng),INTEGRAND_GL(Ng),DS_DH(Ng),VETA(Ng)

    Real(r8) :: ETAP, ETAZ, ZH, ZL, HD, HH, HL, HTAN2, PH, PL, RC, SA, &
             &  SB, SING, DS
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
!
    do j = 1, no_gl_ndx
!
      mp = gl_ndx(j,2)
      if (mp >= brkpt) EXIT
!
      h_i = gl_ndx(j,1)
      if(h_i >= mid) EXIT
!
      if(Sps_phi_loop%values(h_i,1) > ip  .or. &
     &   Sps_phi_loop%values(h_i,2) < ip  .or. &
     &   Sps_zeta_loop%values(h_i,1) > iz .or. &
     &   Sps_zeta_loop%values(h_i,2) < iz) CYCLE
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
      k = j
      rc = ref_corr(h_i+1)
!
      call define_gl_grid_entities
!
!  Compute the value of the integrand at the starting point of the interval
!  (This is done in order to eliminate the singularities. We call these
!  the 'singularities factors')
!
      sing = 0.0
      Call get_one_eta(zh,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        Call get_one_eta(ph,s_phi_basis,s_np,ip,etap)
        sing = etaz * etap * integrand(l)
      end if
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval:
!
      integrand_GL(1:Ng) = integrand(mp+1:mp+Ng)
!
      Call gauss_legendre
!
    end do

    if(k == no_gl_ndx) Return
!
    j = max(1,k)
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
      if(Sps_phi_loop%values(h_i,1) > ip  .or. &
     &   Sps_phi_loop%values(h_i,2) < ip  .or. &
     &   Sps_zeta_loop%values(h_i,1) > iz .or. &
     &   Sps_zeta_loop%values(h_i,2) < iz) CYCLE
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
      sing = 0.0
      Call get_one_eta(zl,s_z_basis,s_nz,iz,etaz)
      if (etaz > 0.0) then
        Call get_one_eta(pl,s_phi_basis,s_np,ip,etap)
        sing = etaz * etap * integrand(mp)
      end if
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
    Real(r8) :: aym, q, v
!
! Define the various GL grid entities for this sub-interval:
!
      aym = 0.5 * abs(zh - zl)
      h_GL(1:Ng) = h_path%values(mp+1:mp+Ng) + elvar%RoC
      DS_DH(1:Ng) = h_GL(1:Ng) / Sqrt(h_GL(1:Ng)**2-htan2)
      Gw_dHdZ(1:Ng) = Gw(1:Ng) * dhdz_path%values(mp+1:mp+Ng) * aym
!
      do i = 1, Ng
        j = mp + i
        Call get_one_eta(z_path%values(j),s_z_basis,s_nz,iz,v)
        Call get_one_eta(phi_path%values(j),s_phi_basis,s_np,ip,q)
        veta(i) = v * q
      end do

    End subroutine DEFINE_GL_GRID_ENTITIES
! -----------------------------------------     GAUSS_LEGENDRE     -----
    Subroutine GAUSS_LEGENDRE
!
    Real(r8) :: fv, q
!
! Now compute the Gauss-Legendre quadrature, subtructing the 'singularities
! factors':
!
      q = SUM((integrand_GL(1:)*veta(1:)-sing)*DS_DH(1:)*Gw_dHdz(1:))
!
! Now add the 'singularities factors' back in, multiplied by their
! respective analytical integral:
!
      fv = q + sing * ds
!
! And Finally - define the delta:
!
      delta(h_i) = fv * rc * fq             ! for (iz,ip)
!
    End subroutine GAUSS_LEGENDRE

  End Subroutine generic_delta_integral

End module GENERIC_DELTA_INTEGRAL_M
! $Log$
! Revision 1.7  2001/06/22 04:05:28  zvi
! Fixing indexing error
!
! Revision 1.6  2001/06/21 13:07:08  zvi
! Speed enhancement MAJOR update
!
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
