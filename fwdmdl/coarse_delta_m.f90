! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module COARSE_DELTA_M
  use GLNP, only: NG
  use ELLIPSE_M, only: ELLIPSE
  use MLSCommon, only: I4, R8
  use PATH_ENTITIES_M, only: PATH_VECTOR, PATH_BETA, PATH_INT_VECTOR_2D
  implicit NONE
  private
  public :: COARSE_DELTA
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------
!  This routine computes the d_delta function. Integration is done
!  using the Gauss-Legendre method.
!
!  ** NOTE: This routine integrate in ZETA Space !
!
  Subroutine COARSE_DELTA(mid,brkpt,no_ele,h_path,beta_path,spsfunc_path,&
          &  n_sps,N_lvls,ref_corr,elvar,midval_delta,Ier)
!
    Integer(i4), intent(in) :: N_SPS, N_LVLS, MID, BRKPT, NO_ELE

    Real(r8), intent(in) :: REF_CORR(:)

    Real(r8), intent(out) :: midval_delta(:,:)     ! (N2lvl,Nsps)

    Integer(i4), intent(out) :: IER

    Type(ELLIPSE), intent(in out) :: elvar

    Type(path_beta), intent(in) :: BETA_PATH(:)      ! (Nsps)

    Type(path_vector), intent(in) :: H_PATH
    Type(path_vector), intent(in) :: SPSFUNC_PATH(:)
!
    Integer(i4) :: J, MP, H_I, HEND, Ngp1
    Real(r8) :: HD, HH, HL, HTAN2, RC, SA, SB, FH, FL
!
! -----     Executable statements     ----------------------------------
!
    Ier = 0
    Ngp1 = Ng + 1
    htan2 = elvar%ht2
    hend = 2 * N_lvls - 1
!
    do j = 1, n_sps
!
! Initialize some arrays:
!
      midval_delta(1:hend+1,j) = 0.0
!
! First, do the right hand side of the ray path:
!
      mp = 1
      elvar%ps = -1.0
      hh = h_path%values(mp)
      fh = spsfunc_path(j)%values(mp) * beta_path(j)%values(mp)
!
      hd = hh + elvar%RoC
      sb = Sqrt(abs(hd*hd-htan2))
!
      do h_i = 1, mid
!
        mp = mp + Ngp1
        if (mp > brkpt) EXIT
!
        hl = hh
        hh = h_path%values(mp)
        if (hh < elvar%ht-0.001) EXIT
!
        fl = fh
        fh = spsfunc_path(j)%values(mp) * beta_path(j)%values(mp)
!
        rc = ref_corr(h_i+1) * (fl + fh) / 2.0
!
        sa = sb
        hd = hh + elvar%RoC
        sb = Sqrt(abs(hd*hd-htan2))
!
        if (abs(sa-sb) < 0.05) EXIT
!
        midval_delta(h_i,j) = rc * abs(sa-sb)
!
      end do
!
! Second, do the left hand side of the ray path:
!
      h_i = mid
      mp = brkpt + 1
      elvar%ps = 1.0
      hh = h_path%values(mp)
      do while (hh < elvar%ht-0.0001)
        h_i = h_i + 1
        mp = mp + Ngp1
        hh = h_path%values(mp)
      end do
!
      fh = spsfunc_path(j)%values(mp) * beta_path(j)%values(mp)
!
      hd = hh + elvar%RoC
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
        fl = fh
        fh = spsfunc_path(j)%values(mp) * beta_path(j)%values(mp)

        rc = ref_corr(h_i) * (fl + fh) / 2.0
!
        sa = sb
        hd = hh + elvar%RoC
        sb = Sqrt(abs(hd*hd-htan2))
!
        midval_delta(h_i,j) = rc * abs(sa-sb)
!
        h_i = h_i + 1
!
      end do
!
    end do
!
    Return

 End Subroutine COARSE_DELTA

End module COARSE_DELTA_M
! $Log$
! Revision 1.0  2001/06/20 22:19:36  Z.Shippony
! Initial release
