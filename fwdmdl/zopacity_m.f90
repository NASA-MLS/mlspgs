module ZOPACITY_M
  use MLSCommon, only: I4, R8
  implicit NONE
  private
  public :: ZOPACITY
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This subroutine computes the incremental opacity matrices.
!
  Subroutine ZOPACITY ( mr_f, n_coeffs, n_sps, Ngp1, &
 &           N2lvl, Nc, mnp, no_phi_f, brkpt, no_ele, delta, del_opacity)
!
    integer(i4), intent(in) :: NC, MNP,Ngp1,brkpt,no_ele
    real(r8), intent(in) :: MR_F(:,:,:)
    integer(i4), intent(in) :: N_COEFFS(:), N_SPS
    integer(i4), intent(in) :: N2lvl, NO_PHI_F(:)
    real(r8), intent(in) :: DELTA(:,:,:,:)
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
  contains
! --------------------------------     SUM_INCREMENTAL_OPACITY     -----
    subroutine SUM_INCREMENTAL_OPACITY
      integer(i4) :: FI, J, K, NF, NO_C
    Real(r8) :: SUM
!
! initialize sums
!
      Sum = 0.0
!
      do j = 1, n_sps                  ! Run over the species list
        nf = no_phi_f(j)
        no_c = n_coeffs(j)
        do fi = 1, nf                  ! Loop over Phi's
          do k = 1, no_c               ! Run over the coefficients
!
! Compute the total incremental opacity for the segment
!
            Sum = Sum + mr_f(k,fi,j) * delta(h_i,k,fi,j)
!
          end do
        end do
      end do
!
      del_opacity(h_i) = Sum
    end subroutine SUM_INCREMENTAL_OPACITY
  End subroutine ZOPACITY
end module ZOPACITY_M
! $Log$
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
