module GET_ETA_M
  use MLSCommon, only: I4, R8
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
    Real(r8), intent(in) :: X(*), PEAKS(*)
    Integer(i4), intent(in) :: NO_X, NO_PEAKS, Nlvl
    Real(r8), intent(out) :: ETA(Nlvl,*)
!
    Integer(i4) :: I, J
    Real(r8) :: R
!
! The first coefficient is one for all values of x below peaks(1)
! until x = peaks(1),then it ramps down in the usual triangular sense
! i is the independent variable x index and j is the coefficient index
!
    i = 1
    do while (i <= no_x .and. peaks(1) > x(i))
      eta(i,1) = 1.0
      eta(i,2:no_peaks) = 0.0
      if ( i == no_x ) exit
      i = i + 1
    end do
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
      if ( i == no_x ) exit
      i = i + 1
    end do
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
! Revision 1.4  2001/01/31 01:08:48  zvi
! New version of forward model
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
