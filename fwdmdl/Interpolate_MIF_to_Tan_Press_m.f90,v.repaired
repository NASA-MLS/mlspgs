! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Interpolate_MIF_to_Tan_Press_m

  ! Interpolate minor frame quantities from PTan to tangent pressures to be
  ! used for ray tracing

  implicit NONE
  private
  public :: Get_Lines_of_Sight, Interpolate_MIF_to_Tan_Press

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  subroutine Get_Lines_of_Sight ( MAF, PTan, Tan_Press, Q_LOS, Path )

    ! Interpolate lines of sight from the PTan pressures for minor frame Q_LOS
    ! to Tan_Press pressures for Path.  Probably only used for QTM.
    ! This only computes Path%Lines(:,1), not Path%Lines(:,2), i.e., only from
    ! the reference to the tangent or Earth-surface intersection.

    use Insertion_Sort_m, only: Insertion_Sort
    use Path_Representation_m, only: Path_t
    use MLSKinds, only: RK => RP
    use MLSNumerics, only: Coefficients, InterpolateArraySetup, &
      & InterpolateArrayTeardown, InterpolateValues
    use VectorsModule, only: VectorValue_T

    ! Inputs
    integer, intent(in) :: MAF                    ! MAF under consideration
    type(vectorValue_T), intent(in) :: PTan       ! Tangent pressure (minor
                                                  ! frame), zeta
    real(rk), intent(in) :: Tan_Press(:)          ! Tangent pressure levels
                                                  ! where the output quantities
                                                  ! are desired, zeta.
    type(Path_t), intent(in):: Q_LOS(:)           ! Minor frame LOS as C + s U
                                                  ! C is ScECR_MIF, U is ECRtoFOV

    ! Output
    type(path_t), intent(out) :: Path(:)          ! Interpolated to Tan_Press

    ! Local variables
    type(coefficients(rk)) :: Coeffs
    integer :: I, J
    integer :: P(ptan%template%noSurfs) ! Permutation result of sorting Z_MIF
    real(rk) :: Q_LOS_P(size(Q_LOS))              ! From Q_LOS(p)%lines(i,1)%xyz
    real(rk) :: Z_MIF(0:ptan%template%noSurfs)    ! Zetas for sorting MIFs
                                                  ! 0'th is a sentinel

    ! Since the interpolateValues routine needs the OldX array to be sorted
    ! we have to sort ptan%values and re-arrange phitan%values, scgeocalt%values
    ! and losvel%values accordingly

    z_mif(0) = -0.5 * huge(1.0_rk)    ! Sentinel to simplify sort. Use 0.5*Huge
                                      ! instead of Huge in case some compiler
                                      ! uses subtraction to compare.
    z_mif(1:) = ptan%values(1:,maf)   ! Zeta

    ! Sort Z_MIF.
    ! Use insertion sort since it might be nearly in order already.
    ! Compute the permutation P of Z_MIF that put it into order, to permute
    ! Q_LOS%value3, from which we interpolate.
    call insertion_sort ( z_mif, p )

    call interpolateArraySetup ( z_mif(1:), tan_press, METHOD = 'L', &
      COEFFS=coeffs, EXTRAPOLATE='A' )

    do j = 1, 3 ! XYZ
      do i = 1, 2 ! 1 = C, 2 = U
        Q_LOS_P = Q_LOS(p)%lines(i,1)%xyz(j)
        ! Extrapolate for subsurface values, if necessary.
        call interpolateValues ( coeffs, z_mif(1:), Q_LOS_P, &
                               & tan_press, path%lines(i,1)%xyz(j), &
                               & METHOD = 'L', EXTRAPOLATE='A' )
      end do
    end do

    call interpolateArrayTeardown ( coeffs )

  end subroutine Get_Lines_of_Sight

  subroutine Interpolate_MIF_to_Tan_Press ( Nlvl, MAF, PTan, &
                              & PhiTan,  ScGeocAlt,     LOSVel, &
                              & Tan_Press, &
                              & Tan_Phi, Est_ScGeocAlt, Est_LOS_Vel, &
                              & Q_EarthRadC_sq, Q_TanHt, &
                              & EarthRadC,      TanHt, &
                              & Tan_Pt_Geod )

    ! Interpolate minor frame quantities PhiTan, ScGeocAlt, and LOSVel, at
    ! pressure levels given by PTan, to Tan_Phi, SC_Geoc_Alt and LOS Velocity,
    ! at pressure levels given by Tan_Press.  If Q_... are present, interpolate
    ! them too.  If Tan_Pt_Geod is present, interpolate it too.

    use Geolocation_0, only: H_V_Geod, Lon_T
    use Insertion_Sort_m, only: Insertion_Sort
    use MLSKinds, only: RK => RP
    use MLSNumerics, only: InterpolateValues
    use Constants, only: Deg2Rad
    use VectorsModule, only: VectorValue_T

  ! Inputs
    integer, intent(in) :: NLvl                   ! Size of integration grid
    integer, intent(in) :: MAF                    ! MAF under consideration
    type(vectorValue_T), intent(in) :: PTan       ! Tangent pressure (minor
                                                  ! frame), zeta
    type(vectorValue_T), intent(in) :: PhiTan     ! Tangent geodAngle (minor
                                                  ! frame), degrees
    type(vectorValue_T), intent(in) :: ScGeocAlt  ! S/C geocentric altitude
                                                  ! (minor frame), meters
    type(vectorValue_T), intent(in) :: LOSVel     ! Line of sight velocity
                                                  ! (minor frame), m/s
    real(rk), intent(in) :: Tan_Press(:)          ! Tangent pressure levels
                                                  ! where the output quantities
                                                  ! are desired, zeta.

  ! Outputs, on pressure levels given by Tan_Press
    real(rk), intent(out) :: Tan_Phi(:)       ! Radians
    real(rk), intent(out) :: Est_ScGeocAlt(:) ! Est S/C geocentric altitude /m
    real(rk), intent(out) :: Est_LOS_Vel(:)

  ! Optional minor frame inputs, for QTM-based model
    type(vectorValue_T), intent(in), optional :: Q_EarthRadC_sq ! square of orbit-
                                                  ! plane projected ellipse, m^2
    type(vectorValue_t), intent(in), optional :: Q_TanHt   ! Geodetic height,
                                                  ! above plane-projected
                                                  ! ellipse, meters
  ! Optional outputs, interpolated to Tan_Press, for QTM-based model
    real(rk), intent(out), optional :: EarthRadC(:) ! Plane-projected
                                                  ! minor axis**2, m**2
    real(rk), intent(out), optional :: TanHt(:)   ! Tangent height, meters
    type(H_V_Geod), intent(out), optional :: Tan_Pt_Geod(:) ! Tangent point
                                                  ! coordinates at Tan_Press
                                                  ! (lon,lat,geod ht)
                                                  ! degrees, meters

  ! Local variables
    integer :: I, Sub
!   real(rk), dimension(merge(size(tan_press),0,present(tan_pt_Geod))) :: Lat, Lon
    real(rk), dimension(size(tan_press)) :: Lat, Lon
    integer :: P(ptan%template%noSurfs) ! Permutation result of sorting Z_MIF
    real(rk) :: Z_MIF(0:ptan%template%noSurfs)    ! Zetas for sorting MIFs
                                                  ! 0'th is a sentinel

    sub = size(tan_press) - nlvl ! # subsurface levels = SurfaceTangentIndex-1

    ! Use first MIF value for subsurface values
    tan_phi(1:sub) = phitan%values(1,maf)
    est_scgeocalt(1:sub) = scGeocalt%values(1,maf)
    est_los_vel(1:sub) = losvel%values(1,maf)

    ! Since the interpolateValues routine needs the OldX array to be sorted
    ! we have to sort ptan%values and re-arrange phitan%values, scgeocalt%values
    ! and losvel%values accordingly

    z_mif(0) = -0.5 * huge(1.0_rk)    ! Sentinel to simplify sort. Use 0.5*Huge
                                      ! instead of Huge in case some compiler
                                      ! uses subtraction to compare.
    z_mif(1:) = ptan%values(1:,maf)   ! Zeta

    ! Sort Z_MIF.
    ! Use insertion sort since it might be nearly in order already.
    ! Compute the permutation P of Z_MIF that put it into order, to permute
    ! other minor-frame quantities from which we interpolate.
    call insertion_sort ( z_mif, p )

    ! Interpolate from minor frame zetas to pointing grid zetas.
    ! Tangent phi.  This actually isn't interesting for QTM.
    call interpolateValues ( z_mif(1:),         phitan%values(p,maf), &
                           & tan_press(sub+1:), tan_phi(sub+1:), &
                           & METHOD = 'L', EXTRAPOLATE='C' )
    tan_phi = tan_phi * deg2rad
    ! Geocentric altitude
    call interpolateValues ( z_mif(1:),         scgeocalt%values(p,maf), &
                           & tan_press(sub+1:), est_scgeocalt(sub+1:), &
                           & METHOD='L', EXTRAPOLATE='C' )
    ! LOS velocity
    call interpolateValues ( z_mif(1:),         losvel%values(p,maf), &
                           & tan_press(sub+1:), est_los_vel(sub+1:), &
                           & METHOD='L', EXTRAPOLATE='C' )

    if ( present(Q_EarthRadC_sq) ) then
      earthRadC(:sub) = Q_EarthRadC_sq%values(1,maf)
      call interpolateValues ( z_mif(1:),         Q_EarthRadC_sq%values(p,maf), &
                             & tan_press(sub+1:), earthRadC(sub+1:), &
                             & METHOD = 'L', EXTRAPOLATE='C' )
    end if

    if ( present(Q_TanHt) ) then
      ! Extrapolate subsurface values if necessary
      call interpolateValues ( z_mif(1:), Q_TanHt%values(p,maf), &
                             & tan_press, tanHt, &
                             & METHOD = 'L', EXTRAPOLATE='A' )
      if ( present(tan_pt_Geod) ) then
        call interpolateValues ( z_mif(1:), phitan%template%geodLat(p,maf), &
                               & tan_press, lat, &
                               & METHOD = 'L', EXTRAPOLATE='A' )
        call interpolateValues ( z_mif(1:), phitan%template%lon(p,maf), &
                               & tan_press ,lon, &
                               & METHOD = 'L', EXTRAPOLATE='A' )
        ! Subsurface values all the same.  These might also be changed to MIF
        ! some day.
        do i = 1, size(tan_pt_geod)
          tan_pt_geod(i) = h_v_geod ( lon_t(lon(i)), lat(i), tanHt(i) )
        end do
      end if
    end if

  end subroutine Interpolate_MIF_to_Tan_Press

! ------------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Interpolate_MIF_to_Tan_Press_m

! $Log$
! Revision 2.4  2017/11/01 00:10:02  vsnyder
! Use setup/teardown with coefficients
!
! Revision 2.3  2017/10/31 17:36:15  vsnyder
! Change QTM path from vector quantity to Path_t
!
! Revision 2.2  2016/11/11 02:01:48  vsnyder
! Add Get_Lines_of_Sight to compute LOS to Tan_Press levels from MIF LOS.
! Eliminate LOS from computations done by Interpolate_MIF_to_Tan_Press.
! Interpolate or extrapolate subsurface rays to subsurface pressures, instead
! of using the first MIF direction.
!
! Revision 2.1  2016/11/09 00:24:27  vsnyder
! Moved from FullForwardModel
!
