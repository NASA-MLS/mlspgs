! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module PATH_CONTRIB_M
  use MLSCommon, only: I4, R8
  use GLNP, only: NG
  implicit NONE
  private
  public :: PATH_CONTRIB
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!----------------------------------------------------------------------
! This routine computes the contributions (along the path) of each interval
! of the (coarse) pre-selected integration grid

 Subroutine Path_contrib(tau, brkpt, no_ele, tol, mid, midval_ndx, &
         &  no_midval_ndx, gl_ndx, no_gl_ndx, Ier)
!
    Integer(i4), intent(in) :: BRKPT, NO_ELE

    Integer(i4), intent(out) :: Ier
    Integer(i4), intent(out) :: no_midval_ndx,no_gl_ndx,mid
    Integer(i4), intent(out) :: midval_ndx(:,:), gl_ndx(:,:)

    Real(r8), intent(in) :: tol

    Real(r8), intent(in) :: tau(:)

    Integer(i4) :: Ngp1, mp, h_i, k

    Real(r8) :: r
!
    ier = 0
    Ngp1 = Ng + 1
!
    no_gl_ndx = 0
    no_midval_ndx = 0
!
    h_i = 0
    mp = 1 - Ngp1
    do
      mp = mp + Ngp1
      if(mp > brkpt) EXIT
      h_i = h_i + 1
      r = 250.0*abs(tau(h_i+1)-tau(h_i))
      if(r > tol .or. brkpt - mp <= 2*Ngp1) then
        no_gl_ndx = no_gl_ndx + 1
        gl_ndx(no_gl_ndx,1) = h_i
        gl_ndx(no_gl_ndx,2) = mp
      else
        no_midval_ndx = no_midval_ndx + 1
        midval_ndx(no_midval_ndx,1) = h_i
        midval_ndx(no_midval_ndx,2) = mp
      endif
    end do
!
    mid = h_i
    h_i = mid - 1
    mp = brkpt + 1 - Ngp1
    do
      mp = mp + Ngp1
      if(mp >= no_ele) EXIT
      h_i = h_i + 1
      r = 250.0*abs(tau(h_i+1)-tau(h_i))
      if(r > tol .or. mp - brkpt <= 2*Ngp1) then
        no_gl_ndx = no_gl_ndx + 1
        gl_ndx(no_gl_ndx,1) = h_i
        gl_ndx(no_gl_ndx,2) = mp
      else
        no_midval_ndx = no_midval_ndx + 1
        midval_ndx(no_midval_ndx,1) = h_i
        midval_ndx(no_midval_ndx,2) = mp
      endif
    end do
!
    k = no_gl_ndx + 1
    gl_ndx(k,1) = gl_ndx(k-1,1) + 1
    gl_ndx(k,2) = gl_ndx(k-1,2) + Ngp1

    if(no_midval_ndx > 0) then
      k = no_midval_ndx + 1
      midval_ndx(k,1) = midval_ndx(k-1,1) + 1
      midval_ndx(k,2) = midval_ndx(k-1,2) + Ngp1
    endif
!
    Return
!
 End Subroutine PATH_CONTRIB

End module PATH_CONTRIB_M
! $Log$
! Revision 1.0  2001/06/07 23:30:34  Z.Shippony
! Initial release
