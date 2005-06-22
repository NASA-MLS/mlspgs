! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module GET_ETA_M

  use MLSCommon, only: Ip, Rp
  implicit NONE
  private
  public :: GET_ETA

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

! This subroutine gets the eta function for temperature

  subroutine GET_ETA ( GRID, PEAKS, ETA )

    real(rp), intent(in) :: GRID(:), PEAKS(:)

    real(rp), intent(out) :: ETA(:,:)

    integer(ip) :: I, J,  NO_PEAKS
    real(rp) :: R

! The first coefficient is one for all values of grid below peaks(1)
! until grid = peaks(1),then it ramps down in the usual triangular sense
! i is the independent variable grid index and j is the coefficient index

    no_peaks = size(peaks)
    eta = 0.0_rp
    i = 1

! first basis calculation

    where ( grid <= peaks(1) ) eta(:,1) = 1.0_rp

! Normal triangular function for j=2 to j=no_peaks-1

    do j = 2, no_peaks
      r = peaks(j) - peaks(j-1)
      where (peaks(j-1) < grid .AND. grid <= peaks(j) )
        eta(:,j-1) = (peaks(j  ) - grid) / r
        eta(:,j  ) = (grid - peaks(j-1)) / r
      endwhere
    end do

! last basis calculation

    where ( peaks(no_peaks) < grid ) eta(:,no_peaks) = 1.0_rp

    return
  end subroutine GET_ETA

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module GET_ETA_M

! $Log$
! Revision 2.2  2003/09/16 00:21:07  vsnyder
! Remove unused dummy arguments
!
! Revision 2.1  2002/10/08 17:08:04  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.8.2.2  2001/09/13 11:18:22  zvi
! get the correct eta
!
! Revision 1.1  2000/05/04 18:12:05  vsnyder
! Initial conversion to Fortran 90
!
