module ZOPACITY_M
  use MLSCommon, only: I4, R4, R8
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
  Subroutine ZOPACITY ( mr_f, sps_tbl, n_coeffs, n_sps, N_lvls, &
 &           N2lvl, Nc, mnp, no_phi_f, delta, del_opacity, IndxR, IndxL)
!
    integer(i4), intent(in) :: NC, MNP
    real(r8), intent(in) :: MR_F(Nc,mnp,*)
    integer(i4), intent(in) :: SPS_TBL(*), N_COEFFS(*), N_SPS, N_lvls
    integer(i4), intent(in) :: N2lvl, NO_PHI_F(*)
    real(r8), intent(in) :: DELTA(N2lvl,Nc,mnp,*)
    real(r8), intent(out) :: del_opacity(*)
    integer(i4), intent(in) :: IndxR, IndxL
    integer(i4) :: H_I
!
! Cover all heights
!
! First, right hand side of the ray path:
!
    do h_i = 1, IndxR
      call sum_incremental_opacity
    end do
!
! Second, left hand side of the ray path:
!
    do h_i = IndxL - 1, 2 * N_lvls - 2
      call sum_incremental_opacity
    end do
!
    Return
  contains
! --------------------------------     SUM_INCREMENTAL_OPACITY     -----
    subroutine SUM_INCREMENTAL_OPACITY
      integer(i4) :: FI, J, K, NF, NO_C, SPS_I
    Real(r8) :: SUM
!
! initialize sums
!
      Sum = 0.0
!
      do sps_i = 1, n_sps             ! Run over the species list
        j = sps_tbl(sps_i)
        nf = no_phi_f(j)
        no_c = n_coeffs(j)
        do fi = 1, nf                 ! Loop over Phi's
          do k = 1, no_c              ! Run over the coefficients
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
! Revision 1.1  2000/05/04 18:12:06  vsnyder
! Initial conversion to Fortran 90
!
