! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module RAD_DELTA_INTEGRAL_M
  use MLSCOmmon, only: I4, R8
  use GLNP, only: GW, NG
  use ELLIPSE_M, only: ELLIPSE
  use PATH_ENTITIES_M, only: PATH_VECTOR
  implicit NONE
  private
  public :: RAD_DELTA_INTEGRAL
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
  Subroutine rad_delta_integral(mid,brkpt,no_ele,z_path,h_path,dhdz_path, &
          &  N_lvls,ref_corr,integrand,elvar,gl_ndx,no_gl_ndx,delta,Ier)
!
  Integer(i4), intent(in) :: N_lvls, mid, brkpt, no_ele, no_gl_ndx
!
  Integer(i4), intent(in) :: gl_ndx(:,:)
!
  Type(path_vector), intent(in) :: Z_PATH, H_PATH, DHDZ_PATH

  Real(r8), intent(in) :: REF_CORR(:), INTEGRAND(:)
!
  Type(ELLIPSE), intent(in out) :: elvar

  Real(r8), intent(out) :: DELTA(:)

  Integer(i4), intent(out) :: IER

  Integer(i4) :: J, K, L, MP, NGP1, H_I, HEND

  Real(r8) :: DS_DH(Ng), H_GL(Ng), GW_DHDZ(Ng), INTEGRAND_GL(Ng)

  Real(r8) :: HD, HTAN2, RC, SA, SB, SING, DS, AYM, q, fv
!
    Ier = 0
!
!  Initialize the delta array:
!
    hend = 2 * N_lvls - 1
    delta(1:hend+1) = 0.0
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
      hd = h_path%values(mp) + elvar%RoC
      sa = Sqrt(abs(hd*hd-htan2))
!
      l = mp + Ngp1
      hd = h_path%values(l) + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
      ds = abs(sa-sb)
!
      if (ds < 0.05) EXIT
!
      k = j
      rc = ref_corr(h_i+1)
!
      h_GL(1:Ng) = h_path%values(mp+1:mp+Ng) + elvar%RoC
      DS_DH(1:Ng) = h_GL(1:Ng) / Sqrt(h_GL(1:Ng)**2 - htan2)
      aym = 0.5 * abs(z_path%values(l)-z_path%values(mp))
      Gw_dHdZ(1:Ng) = Gw(1:Ng) * dhdz_path%values(mp+1:mp+Ng) * aym
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval,
! subtructing the value of the integrand at the starting point of the 
! interval (This is done in order to eliminate the singularities. We call 
! these the 'singularities factors')
!
      integrand_GL(1:Ng) = integrand(mp+1:mp+Ng) - integrand(l)
!
! Now compute the Gauss-Legendre quadrature:
!
      q = SUM(integrand_GL(1:)*DS_DH(1:)*Gw_dHdz(1:))

! Now add the 'singularities factors' back in, multiplied by their
! respective analytical integral (delta_s):
!
      fv = q + integrand(l) * ds
!
! And Finally - define the delta:
!
      delta(h_i) = fv * rc             ! for (iz,ip)
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
      hd = h_path%values(mp) + elvar%RoC
      sa = Sqrt(abs(hd*hd-htan2))
!
      hd = h_path%values(l) + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
      ds = abs(sa-sb)
!
      if (ds < 0.05) EXIT
!
      rc = ref_corr(h_i)
!
      h_GL(1:Ng) = h_path%values(mp+1:mp+Ng) + elvar%RoC
      DS_DH(1:Ng) = h_GL(1:Ng) / Sqrt(h_GL(1:Ng)**2 - htan2)
      aym = 0.5 * abs(z_path%values(l)-z_path%values(mp))
      Gw_dHdZ(1:Ng) = Gw(1:Ng) * dhdz_path%values(mp+1:mp+Ng) * aym
!
! Define integrand on the Gauss-Legendre grid for the current sub-interval,
! subtructing the value of the integrand at the starting point of the 
! interval (This is done in order to eliminate the singularities. We call 
! these the 'singularities factors')
!
      integrand_GL(1:Ng) = integrand(mp+1:mp+Ng) - integrand(mp)
!
! Now compute the Gauss-Legendre quadrature:
!
      q = SUM(integrand_GL(1:)*DS_DH(1:)*Gw_dHdz(1:))

! Now add the 'singularities factors' back in, multiplied by their
! respective analytical integral (delta_s):
!
      fv = q + integrand(mp) * ds
!
! And Finally - define the delta:
!
      delta(h_i) = fv * rc             ! for (iz,ip)
!
    end do
!
    Return

  End Subroutine rad_delta_integral

End module RAD_DELTA_INTEGRAL_M
! $Log$
! Revision 1.0  2001/07/02 17:45:00  zvi
! First version
