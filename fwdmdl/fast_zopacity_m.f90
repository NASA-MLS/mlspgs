! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module FAST_ZOPACITY_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: FAST_ZOPACITY
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
     "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
     "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This subroutine computes the incremental opacity matrices.
!
  Subroutine FAST_ZOPACITY (n_sps,Ngp1,brkpt,no_ele,delta,del_opacity)
!
    integer(i4), intent(in) :: Ngp1,brkpt,no_ele
    integer(i4), intent(in) :: N_SPS
!
    real(r8), intent(in) :: DELTA(:,:)
!
    real(r8), intent(out) :: del_opacity(:)

    integer(i4) :: h_i,i,m,IndxR
!
! Cover all heights
!
! First, right hand side of the ray path:
!
    h_i = 0
    m = 1 - Ngp1
    IndxR = brkpt / Ngp1
    do i = 1, IndxR
      m = m + Ngp1
      if(m > brkpt) EXIT
      h_i = h_i + 1
      call sum_incremental_opacity
    end do
!
! Second, left hand side of the ray path:
!
    m = brkpt + 1
    do i = 1, IndxR
      m = m + Ngp1
      if(m > no_ele) EXIT
      h_i = h_i + 1
      call sum_incremental_opacity
    end do
!
    Return
!
  contains
! --------------------------------     SUM_INCREMENTAL_OPACITY     -----
!
    subroutine SUM_INCREMENTAL_OPACITY
!
      Integer(i4) :: J
!
      Real(r8) :: SUM
!
! Compute the total incremental opacity for the segment
!
      Sum = 0.0
      do j = 1, n_sps                  ! Run over the species list
        Sum = Sum + delta(h_i,j)
      end do
!
      del_opacity(h_i) = Sum
!
    end subroutine SUM_INCREMENTAL_OPACITY
!
  End subroutine FAST_ZOPACITY
!
end module FAST_ZOPACITY_M
! $Log$
! Revision 1.4  2001/03/31 23:40:55  zvi
! Eliminate l2pcdim (dimension parameters) move to allocatable ..
!
! Revision 1.3  2001/03/29 08:51:01  zvi
! Changing the (*) toi (:) everywhere
!
! Revision 1.2  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
