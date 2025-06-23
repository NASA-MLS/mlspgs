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
  use MLSKinds, only: RP, RV

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

  integer, parameter :: NG = 8 ! Order of Gauss quadrature
  real(rp), parameter :: GX(ng) = (/ &          ! Gauss abscissae
    & -.960289856497536231684_rp, -.796666477413626739592_rp, &
    & -.525532409916328985818_rp, -.183434642495649804939_rp, &
    & 0.183434642495649804939_rp, 0.525532409916328985818_rp, &
    & 0.796666477413626739592_rp, 0.960289856497536231684_rp /)
  real(rp), parameter :: GW(ng) = 0.5_rp * (/ & ! Gauss weights
    & 0.101228536290376259153_rp, 0.222381034453374470544_rp, &
    & 0.313706645877887287338_rp, 0.362683783378361982965_rp, &
    & 0.362683783378361982965_rp, 0.313706645877887287338_rp, &
    & 0.222381034453374470544_rp, 0.101228536290376259153_rp /)

contains

  ! ---------------------------------------------  ScanAverage_1d  -----
  subroutine ScanAverage_1d ( MIF_Times, DeadTime, Chi_In, Chi_Out, Y_in, &
    &                         Y_Out, DY_DX_Out )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume the MIFs begin at MIF_Times + deadTime.

    use MLSFILLVALUES, only: ISNAN
    use MLSMESSAGEMODULE, only: MLSMESSAGE, MLSMSG_WARNING
    use MLSNUMERICS, only: INTERPOLATEVALUES

    real(rv), intent(in) :: MIF_Times(:)   ! MIF starts here + deadTime
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_In(:)   ! For Y_In
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  Same size as MIF_Times
    real(rp), intent(in) :: Y_In(:)     ! Quantity to be scan averaged
                                        ! size(Chi_in)
    real(rp), intent(out) :: Y_Out(:)   ! Scan averaged quantity
                                        ! size(Chi_out)
    real(rp), intent(out), optional :: DY_DX_Out(:)

    integer :: I ! Subscript and loop inductor

    real(rp), dimension(ng*size(chi_out)) :: X     ! Quadrature abscissae
    real(rp), dimension(ng*size(chi_out)) :: Y, DY ! Quadrature ordinates

    call scanAverage_prep ( MIF_Times, deadTime, chi_out, x )

    if ( present(dY_dX_out) ) then
      ! Interpolate from (Chi_In,Y_in) to (x,y) and to (x,dy)
      call interpolateValues ( chi_in, y_in, x, y, &
                             & METHOD='S', extrapolate='C', dYbYdX=dY )
      if ( any(isNaN(y)) ) &
        & call MLSMessage( MLSMSG_Warning, ModuleName // 'ScanAverage_1D', &
        & 'NaNs returned by 1st interpolation' )
      ! Integrate panels of y and dy
      do i = 1, size(y_out)
        y_out(i) = dot_product(y(1+(i-1)*ng:i*ng),gw)
        dy_dx_out(i) = dot_product(dy(1+(i-1)*ng:i*ng),gw)
      end do
    else
      ! Interpolate from (Chi_In,Y_in) to (x,y)
      call interpolateValues ( chi_in, y_in, x, y, &
                             & METHOD='S', extrapolate='C' )
      if ( any(isNaN(y)) ) &
        & call MLSMessage( MLSMSG_Warning, ModuleName // 'ScanAverage_1D', &
        & 'NaNs returned by 2nd interpolation' )
      ! Integrate panels of y
      do i = 1, size(y_out)
        y_out(i) = dot_product(y(1+(i-1)*ng:i*ng),gw)
      end do
    end if

  end subroutine ScanAverage_1d

  ! ---------------------------------------------  ScanAverage_2d  -----
  subroutine ScanAverage_2d ( MIF_Times, DeadTime, Chi_In, Chi_Out, Y_in, &
    &                         Y_Out, DY_DX_Out )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume MIF_Times are midway between T1 and T2-deadTime, where T1 and
    ! T2 are the beginning and ending times for MIF integration.

    use MLSNUMERICS, only: Coefficients, INTERPOLATEARRAYSETUP, &
      & INTERPOLATEARRAYTEARDOWN, INTERPOLATEVALUES

    real(rv), intent(in) :: MIF_Times(:) ! MIF starts here + deadTime
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_In(:)   ! For Y_In
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  Same size as MIF_Times
    real(rp), intent(in) :: Y_In(:,:)   ! Quantity to be scan averaged
                                        ! size(Chi_in) x ???
    real(rp), intent(out) :: Y_Out(:,:) ! Scan averaged quantity
                                        ! size(Chi_out) x size(Y_in,2)
    real(rp), intent(out), optional :: DY_DX_Out(:,:)

    type(coefficients(rp)) :: Coeffs        ! to interpolate from Chi_In to X

    integer :: I, J ! Subscripts and loop inductors

    real(rp), dimension(ng*size(chi_out)) :: X     ! Quadrature abscissae
    real(rp), dimension(ng*size(chi_out)) :: Y, DY ! Quadrature ordinates

    call scanAverage_prep ( MIF_Times, deadTime, chi_out, x )

    call interpolateArraySetup ( chi_in, x, 'S', coeffs, extrapolate='C', &
      & DyByDx=present(dY_dX_out) )

    do j = 1, size(y_out,2)
      if ( present(dY_dX_out) ) then
        ! Interpolate from (Chi_In,Y_in) to (x,y) and to (x,dy)
        call interpolateValues ( coeffs, chi_in, y_in(:,j), x, y, &
                               & METHOD='S', extrapolate='C', dYbYdX=dY )
        ! Integrate panels of y and dy
        do i = 1, size(y_out)
          y_out(i,j) = dot_product(y(1+(i-1)*ng:i*ng),gw)
          dy_dx_out(i,j) = dot_product(dy(1+(i-1)*ng:i*ng),gw)
        end do
      else
        ! Interpolate from (Chi_In,Y_in) to (x,y)
        call interpolateValues ( coeffs, chi_in, y_in(:,j), x, y, &
                               & METHOD='S', extrapolate='C' )
        ! Integrate panels of y
        do i = 1, size(y_out)
          y_out(i,j) = dot_product(y(1+(i-1)*ng:i*ng),gw)
        end do
      end if
    end do

    call interpolateArrayTeardown ( coeffs )

  end subroutine ScanAverage_2d

  ! ---------------------------------------------  ScanAverage_Prep  -----
  subroutine ScanAverage_Prep ( MIF_Times, DeadTime, Chi_Out, X, Pos1P )
    ! Average over each MIF.  Assume linear motion during the MIF.
    ! Assume MIF_Times are midway between T1 and T2-deadTime, where T1 and
    ! T2 are the beginning and ending times for MIF integration.

    use ALLOCATE_DEALLOCATE, only: ALLOCATE_TEST, DEALLOCATE_TEST
    use MLSNUMERICS, only: INTERPOLATEVALUES
    use SORT_M, only: SORTP

    real(rv), intent(in) :: MIF_Times(:) ! MIF starts here + deadTime
    real(rv), intent(in) :: DeadTime ! How much of each MIF is not collecting data
    real(rp), intent(in) :: Chi_Out(:)  ! For Y_Out.  No larger than MIF_Times
    real(rp), intent(out) :: X(:)       ! Quadrature abscissae
    real(rp), intent(in), optional :: Pos1P(:) ! Positions (Chi) at MIF_Times.
                                        ! If present, Chi_Out is assumed to be
                                        ! Pos2P.  Size(Chi_Out).

    integer :: I     ! Subscript and loop inductor
    integer :: N     ! Size(chi_Out)

    real(rp), dimension(ng*size(chi_out)) :: T     ! Quadrature abscissae
    real(rp), dimension(size(chi_out)) :: TB, TE   ! Bounding times
    real(rp), dimension(size(chi_out)) :: TBAR     ! Time at Chi_Out

    ! If Pos1P is present
    integer, dimension(:), pointer :: I_POS
    real(rp), dimension(:), pointer :: X_POS, T_POS

    n = size(chi_out)
    ! Integration start times for each MIF
    tb = MIF_Times(:n) + deadTime
    ! Integration end times for each MIF
    te(1:n-1) = MIF_Times(2:n)
    te(n) = 2.0*MIF_Times(n) - MIF_Times(n-1)
    ! Times corresponding to Chi_Out
    tbar(1:n-1) = 0.5 * ( MIF_Times(1:n-1) + MIF_Times(2:n) )
    tbar(n) = 0.5 * ( MIF_Times(n) + te(n) )
    ! Times at quadrature abscissae
    do i = 1, n
      t(1+(i-1)*ng:i*ng) = 0.5 * ( tb(i) + te(i) + ( te(i) - tb(i) ) * gx )
    end do

    if ( present(pos1P) ) then
      call allocate_test ( i_pos, 2*size(chi_out), 'I_Pos', moduleName )
      call allocate_test ( t_pos, 2*size(chi_out), 'T_Pos', moduleName )
      call allocate_test ( x_pos, 2*size(chi_out), 'X_Pos', moduleName )
      x_pos(1:n) = chi_out
      x_pos(n+1:) = pos1P
      t_pos(1:n) = tbar
      t_pos(n+1:) = MIF_Times
      call sortp ( t_pos, 1, 2*n, i_pos )
      call interpolateValues ( t_pos(i_pos), x_pos(i_pos), t, x, &
                             & METHOD='S', extrapolate='C' )
      call deallocate_test ( i_pos, 'I_Pos', moduleName )
      call deallocate_test ( t_pos, 'T_Pos', moduleName )
      call deallocate_test ( x_pos, 'X_Pos', moduleName )
    else
      ! Interpolate from (tbar,Chi_Out) to (t,x)
      call interpolateValues ( tbar, Chi_Out, t, x, &
                             & METHOD='S', extrapolate='C' )
    end if

  end subroutine ScanAverage_Prep

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ScanAverage_m

! $Log$
! Revision 2.8  2017/10/31 23:49:36  vsnyder
! Make Coefficients a parameterized type
!
! Revision 2.7  2011/05/09 17:54:29  pwagner
! Warns of NaNs retruned by interpolations
!
! Revision 2.6  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.5  2008/05/02 00:46:56  vsnyder
! Remove unused symbols
!
! Revision 2.4  2007/11/29 23:28:22  vsnyder
! Repair some comments
!
! Revision 2.3  2005/09/03 01:21:59  vsnyder
! Completely revised
!
! Revision 2.2  2005/08/06 01:41:10  vsnyder
! Get rid of coeffs, use trapezoid rule with endpoint corrections
!
! Revision 2.1  2005/08/03 02:31:08  vsnyder
! Initial commit
!
