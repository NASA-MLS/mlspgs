module GET_ETA_M
  use MLSCommon, only: Ip, Rp
  implicit NONE
  private
  public :: GET_ETA
!---------------------------- RCS Ident Info -------------------------------
  private ID, ModuleName
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------
contains
! This subroutine gets the eta function for temperature
!
  Subroutine GET_ETA ( GRID, PEAKS, NO_GRID, NO_PEAKS, ETA )

    Real(rp), intent(in) :: GRID(:), PEAKS(:)
    Integer(ip), intent(in) :: NO_GRID, NO_PEAKS

    Real(rp), intent(out) :: ETA(:,:)
!
    Integer(ip) :: I, J
    Real(rp) :: R
!
! The first coefficient is one for all values of grid below peaks(1)
! until grid = peaks(1),then it ramps down in the usual triangular sense
! i is the independent variable grid index and j is the coefficient index
!
    eta = 0.0_rp
    i = 1
!
! first basis calculation
!
    WHERE(grid <= peaks(1)) eta(:,1) = 1.0_rp
!
! Normal triangular function for j=2 to j=no_peaks-1
!
    DO j = 2, no_peaks
      r = peaks(j) - peaks(j-1)
      WHERE(peaks(j-1) < grid .AND. grid <= peaks(j))
        eta(:,j-1) = (peaks(j  ) - grid) / r
        eta(:,j  ) = (grid - peaks(j-1)) / r
      ENDWHERE
    END DO
!
! last basis calculation
!
    WHERE(peaks(no_peaks) < grid) eta(:,no_peaks) = 1.0_rp

    Return
  End Subroutine GET_ETA
end module GET_ETA_M
! $Log$
! Revision 1.8.2.2  2001/09/13 11:18:22  zvi
! get the correct eta
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
