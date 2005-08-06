! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ScanAverage_m

  implicit NONE
  private

  public ScanAverage

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface ScanAverage
    module procedure ScanAverage_1d, ScanAverage_2d
  end interface

contains

  ! ---------------------------------------------  ScanAverage_1d  -----
  subroutine ScanAverage_1d ( MIF_Times, DeadTime, Chi_In, Chi_Out, Y_in, &
    &                         Y_Out, DY_DX_Out )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume MIF_Times are midway between T1 and T2-deadTime, where T1 and
    ! T2 are the beginning and ending times for MIF integration.

    use MLSCommon, only: RP, RV
    use MLSNumerics, only: InterpolateValues

    real(rv), intent(in) :: MIF_Times(:)
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_In(:)   ! For Y_In
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  Same size as MIF_Times
    real(rp), intent(in) :: Y_In(:)     ! Quantity to be scan averaged
                                        ! size(Chi_in)
    real(rp), intent(out) :: Y_Out(:)   ! Scan averaged quantity
                                        ! size(Chi_out)
    real(rp), intent(out), optional :: DY_DX_Out(:)

    real(rp), dimension(size(chi_out)) :: DYBDX, DYEDX ! \int YB dX, \int YE dX
!   real(rp), dimension(size(chi_out)) :: IntYBDX, IntYEDX ! \int YB dX, \int YE dX
    real(rp), dimension(size(chi_out)) :: TB, TE ! Bounding times
    real(rp), dimension(size(chi_out)) :: XB, XE ! Bounding angles
    real(rp), dimension(size(chi_out)) :: YB, YE ! Y at XB and XE

    call scanAverage_prep ( MIF_Times, deadTime, chi_out, xb, xe )

    ! Scan smooth to get the output
    call InterpolateValues ( chi_in, y_in, xb, yb, &
                           & METHOD='S', extrapolate='C', dYbYdX=dYBdX ) ! intYdx=intYbDx )
    call InterpolateValues ( chi_in, y_in, xe, ye, &
                           & METHOD='S', extrapolate='C', dYbYdX=dYEdX ) ! intYdx=intYeDx )
    y_out = 0.5 * ( yb + ye ) + (xe - xb ) * ( dyEdX - dyBdx ) / 12.0
!   y_out = intYeDx - intYbDx
    if ( present(dY_dX_out) ) dY_dX_out = ye - yb ! \int y' dx = y

  end subroutine ScanAverage_1d

  ! ---------------------------------------------  ScanAverage_1d  -----
  subroutine ScanAverage_2d ( MIF_Times, DeadTime, Chi_In, Chi_Out, Y_in, &
    &                         Y_Out, DY_DX_Out )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume MIF_Times are midway between T1 and T2-deadTime, where T1 and
    ! T2 are the beginning and ending times for MIF integration.

    use MLSCommon, only: RP, RV
    use MLSNumerics, only: Coefficients => Coefficients_r8, InterpolateValues

    real(rv), intent(in) :: MIF_Times(:)
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_In(:)   ! For Y_In
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  Same size as MIF_Times
    real(rp), intent(in) :: Y_In(:,:)   ! Quantity to be scan averaged
                                        ! size(Chi_in) x ???
    real(rp), intent(out) :: Y_Out(:,:) ! Scan averaged quantity
                                        ! size(Chi_out) x size(Y_in,2)
    real(rp), intent(out), optional :: DY_DX_Out(:,:)

    real(rp), dimension(size(chi_out)) :: DYBDX, DYEDX ! \int YB dX, \int YE dX
!   real(rp), dimension(size(chi_out)) :: IntYBDX, IntYEDX ! \int YB dX, \int YE dX
    real(rp), dimension(size(chi_out)) :: R      ! Time or angle difference ratios
    real(rp), dimension(size(chi_out)) :: TB, TE ! Bounding times
    real(rp), dimension(size(chi_out)) :: XB, XE ! Bounding angles
    real(rp), dimension(size(chi_out)) :: YB, YE ! Y at XB and XE

    integer :: I     ! Subscript and loop inductor

    call scanAverage_prep ( MIF_Times, deadTime, chi_out, xb, xe )

    ! Scan smooth to get the output
    r = ( xe - xb ) / 12.0
    do i = 1, size(y_out,2)
      call InterpolateValues ( chi_in, y_in(:,i), xb, yb, &
                               & METHOD='S', extrapolate='C', dYbYdX=dYBdX ) ! intYdx=intYbDx )
      call InterpolateValues ( chi_in, y_in(:,i), xe, ye, &
                               & METHOD='S', extrapolate='C', dYbYdX=dYEdX ) ! intYdx=intYeDx )
      y_out(:,i) = 0.5 * ( yb + ye ) + r * ( dyEdX - dyBdx )
!     y_out(:,i) = intYeDx - intYbDx
      if ( present(dY_dX_out) ) dY_dX_out(:,i) = ye - yb ! \int y' dx = y
    end do

  end subroutine ScanAverage_2d

  ! ---------------------------------------------  ScanAverage_Prep  -----
  subroutine ScanAverage_Prep ( MIF_Times, DeadTime, Chi_Out, XB, XE )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume MIF_Times are midway between T1 and T2-deadTime, where T1 and
    ! T2 are the beginning and ending times for MIF integration.

    use MLSCommon, only: RP, RV
    real(rv), intent(in) :: MIF_Times(:)
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  No larger than MIF_Times
    real(rp), intent(out) :: XB(:), XE(:) ! Bounding angles

    real(rp), dimension(1:size(chi_out)-4) :: A  ! Acceleration
    real(rp), dimension(2:size(chi_out)-5) :: B  ! Acceleration's moment
    real(rp), dimension(size(chi_out)-1) :: R    ! Time or angle difference ratios
    real(rp), dimension(size(chi_out)) :: TB, TE ! Bounding times

    real(rp) :: AMAX ! Maxval(A)
    integer :: I     ! Subscript and loop inductor
    integer :: M     ! minloc(abs(a)+abs(b))
    integer :: N     ! Size(chi_Out)

    n = size(chi_out)

    ! Get the times at the MIF starts and ends
    tb(1) = ( 3.0 * MIF_times(1) - MIF_times(2) + deadTime ) * 0.5
    te(1) = 2.0 * MIF_Times(1) - tb(1)
    do i = 2, n
      tb(i) = 2.0 * MIF_times(i-1) - tb(i-1) + deadTime
      te(i) = 2.0 * MIF_times(i) - tb(i)
    end do

    ! Get the bounding angles.
    ! The accuracy depends critically upon getting a good initial condition.
    ! We choose a place where the scan acceleration is zero or velocity is
    ! most constant.
    ! Compute acceleration:
    a = ((chi_out(5:n)-chi_out(4:n-1)) / ((MIF_times(5:n)-MIF_times(4:n-1)) * &
      &                                   (MIF_times(4:n-1)-MIF_times(3:n-2))) + &
      &  (chi_out(4:n-1)-2.0*chi_out(3:n-2)+chi_out(2:n-3)) / &
      &     ((MIF_times(4:n-1)-MIF_times(3:n-2)) * &
      &      (MIF_times(3:n-2)-MIF_times(2:n-3))) - &
      &  (chi_out(2:n-3)-chi_out(1:n-4)) / ((MIF_times(3:n-2)-MIF_times(2:n-3)) * &
      &  (MIF_times(2:n-3)-MIF_times(1:n-4)))) / 4.0
    amax = maxval(a)
    ! Compute first moment of acceleration
    b = ((a(3:n-4)-a(2:n-5)) / (MIF_times(5:n-2)-MIF_times(4:n-3)) + &
      & (a(2:n-5)-a(1:n-6)) / (MIF_times(4:n-3)-MIF_times(3:n-4)))*0.5
    ! Find where acceleration + moment is minimum.  We add 2 to the index
    ! because acceleration is a second difference.
    m = minloc(abs(a(2:n-5))+abs(b),1) + 2

    ! Estimate an initial value for XB(m)
    xb(m) = chi_out(m-1) + (chi_out(m)-chi_out(m-1)) * (tb(m)-MIF_times(m-1)) / &
      &                  (MIF_times(m)-MIF_times(m-1))

    ! Get time difference ratios
    r = (tb(2:n)-tb(1:n-1)) / (MIF_times(1:n-1)-tb(1:n-1))

    ! Get XB for MIFs above M
    do i = m+1, n
      xb(i) = xb(i-1) + ( chi_out(i-1) - xb(i-1) ) * r(i-1)
    end do
    ! Get XB for MIFs below B
    do i = m-1, 1, -1
      xb(i) = ( xb(i+1) - r(i) * chi_out(i) ) / ( 1.0 - r(i) )
    end do

    ! Get XE
    xe(1:n-1) = xb(1:n-1) + ( xb(2:n) - xb(1:n-1) ) * ( te(1:n-1) - tb(1:n-1) ) / &
      &                                               ( tb(2:n)   - tb(1:n-1) )
    ! With an estimated extrapolated XE(n)
    xe(n) = xb(n) + (( chi_out(n) - xb(n) ) * ( te(n) - tb(n) ) / &
      &                                       ( MIF_times(n) - tb(n)))

    ! We could have replaced this with a simple interpolation from
    ! MIF_times, Chi_out to XB, XE at TB, TE, but that does not
    ! guarantee a constant velocity between Xb and XE, which is a
    ! requirement for the scan smoothing integration algorithm.  In
    ! practice it probably doesn't matter.

  end subroutine ScanAverage_Prep

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ScanAverage_m

! $Log$
! Revision 2.2  2005/08/06 01:41:10  vsnyder
! Get rid of coeffs, use trapezoid rule with endpoint corrections
!
! Revision 2.1  2005/08/03 02:31:08  vsnyder
! Initial commit
!
