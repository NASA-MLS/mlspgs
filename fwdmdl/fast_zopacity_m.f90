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
  Subroutine FAST_ZOPACITY (sps_tbl, n_sps, Ngp1, N2lvl, brkpt, &
 &                          no_ele, delta, del_opacity)
!
    integer(i4), intent(in) :: Ngp1,brkpt,no_ele,N2lvl
    integer(i4), intent(in) :: SPS_TBL(*), N_SPS
!
    real(r8), intent(in) :: DELTA(N2lvl,*)
!
    real(r8), intent(out) :: del_opacity(*)

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
      Integer(i4) :: J, SPS_I
!
      Real(r8) :: SUM
!
! Compute the total incremental opacity for the segment
!
      Sum = 0.0
      do sps_i = 1, n_sps             ! Run over the species list
        j = sps_tbl(sps_i)
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
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
