module GET_ETA_M
  use MLSCommon, only: I4, R4
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
  Subroutine GET_ETA ( X, PEAKS, NO_X, NO_PEAKS, Nlvl, ETA )
    Real(r4), intent(in) :: X(*), PEAKS(*)
    Integer(i4), intent(in) :: NO_X, NO_PEAKS, Nlvl
    Real(r4), intent(out) :: ETA(Nlvl,*)
!
    Integer(i4) :: I, J
    Real(r4) :: R
!
! The first coefficient is one for all values of x below peaks(1)
! until x = peaks(1),then it ramps down in the usual triangular sense
! i is the independent variable x index and j is the coefficient index
!
    i = 1
    do while (i <= no_x .and. peaks(1) > x(i))
      eta(i,1) = 1.0
      eta(i,2:no_peaks) = 0.0
      i = i + 1
    end do
!
    if (i > no_x) Return
!
! Normal triangular function for j=2 to j=no_peaks-1
!
    j = 2
    do while (i <= no_x .and. x(i) < peaks(no_peaks))
      do while (peaks(j) < x(i) .and. j < no_peaks)
        j = j + 1
      end do
      eta(i,1:no_peaks) = 0.0
      r = peaks(j) - peaks(j-1)
      eta(i,j-1) = (peaks(j  ) - x(i)) / r
      eta(i,j  ) = (x(i) - peaks(j-1)) / r
      i = i + 1
    end do
!
    if (i > no_x) Return
!
! The no_peaks coefficient ramps up as a triangle until x = peaks(no_peaks)
! then afterwards it is equal to one
!
    do while (i <= no_x)
      eta(i,1:no_peaks-1) = 0.0
      eta(i,no_peaks) = 1.0
      i = i + 1
    end do
!
    Return
  End Subroutine GET_ETA
end module GET_ETA_M

! $Log$
