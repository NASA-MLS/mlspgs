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
  public :: Interpolate_MIF_to_Tan_Press

!------------------------------ RCS Ident Info -------------------------------
  character (len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
!------------------------------------------------------------------------------

contains

  subroutine Interpolate_MIF_to_Tan_Press ( Nlvl, MAF, PTan, &
                              & PhiTan,  ScGeocAlt,     LOSVel, &
                              & Tan_Press, &
                              & Tan_Phi, Est_ScGeocAlt, Est_LOS_Vel, &
                              & Q_EarthRadC_sq, Q_Incline, Q_LOS, Q_TanHt, &
                              & EarthRadC,      Incline,   LOS,   TanHt, &
                              & Tan_Pt_Geod )

  ! Interpolate minor frame quantities PhiTan, ScGeocAlt, and LOSVel, at
  ! pressure levels given by PTan, to Tan_Phi, SC_Geoc_Alt and LOS Velocity, at
  ! pressure levels given by Tan_Press.  If Q_... are present, interpolate
  ! them too.

    use Geolocation_0, only: ECR_T, H_V_Geod, Lon_T
    use Insertion_Sort_m, only: Insertion_Sort
    use MLSKinds, only: RK => RP
    use MLSNumerics, only: InterpolateValues
    use Constants, only: Deg2Rad
    use VectorsModule, only: VectorValue_T

    implicit NONE

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
    type(vectorValue_t), intent(in), optional :: Q_Incline ! Degrees, geocentric
    type(vectorValue_t), intent(in), optional :: Q_LOS     ! LOS as C + s U
    type(vectorValue_t), intent(in), optional :: Q_TanHt   ! Geodetic height,
                                                  ! above plane-projected
                                                  ! ellipse, meters
  ! Optional outputs, interpolated to Tan_Press, for QTM-based model
    real(rk), intent(out), optional :: EarthRadC(:) ! Plane-projected
                                                  ! minor axis**2, m**2
    real(rk), intent(out), optional :: Incline(:) ! Inclination of the plane
                                                  ! defined by LOS and the
                                                  ! Earth center, degrees.
    type(ECR_t), intent(out), optional :: LOS(:,:)! (1,:) are points on LOS,
                                                  ! (2,:) are unit directions.
    real(rk), intent(out), optional :: TanHt(:)   ! Tangent height, meters
    type(H_V_Geod), intent(out), optional :: Tan_Pt_Geod(:) ! Tangent point
                                                  ! coordinates at Tan_Press
                                                  ! (lon,lat,geod ht)
                                                  ! degrees, meters

  ! Local variables
    integer :: I, J, K, SUB
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

    k = ptan%template%noSurfs

    z_mif(0) = -0.5 * huge(1.0_rk)    ! Sentinel to simplify sort. Use 0.5*Huge
                                      ! instead of Huge in case some compiler
                                      ! uses subtraction to compare.
    z_mif(1:k) = ptan%values(1:k,maf) ! Zeta

    ! Sort Z_MIF.
    ! Use insertion sort since it might be nearly in order already.
    ! Compute the permutation P of Z_MIF that put it into order, to permute
    ! other minor-frame quantities from which we interpolate.
    call insertion_sort ( z_mif, p )

    ! Interpolate from minor frame zetas to pointing grid zetas.
    ! Tangent phi.  This actually isn't interesting for QTM.
    call interpolateValues ( z_mif(1:k),        phitan%values(p,maf), &
                           & tan_press(sub+1:), tan_phi(sub+1:), &
                           & METHOD = 'L', EXTRAPOLATE='C' )
    tan_phi = tan_phi * deg2rad
    ! Geocentric altitude
    call interpolateValues ( z_mif(1:k),        scgeocalt%values(p,maf), &
                           & tan_press(sub+1:), est_scgeocalt(sub+1:), &
                           & METHOD='L', EXTRAPOLATE='C' )
    ! LOS velocity
    call interpolateValues ( z_mif(1:k),        losvel%values(p,maf), &
                           & tan_press(sub+1:), est_los_vel(sub+1:), &
                           & METHOD='L', EXTRAPOLATE='C' )

    if ( present(Q_EarthRadC_sq) ) then
      earthRadC(:sub) = Q_EarthRadC_sq%values(1,maf)
      call interpolateValues ( z_mif(1:k),        Q_EarthRadC_sq%values(p,maf), &
                             & tan_press(sub+1:), earthRadC(sub+1:), &
                             & METHOD = 'L', EXTRAPOLATE='C' )
    end if

    if ( present(Q_Incline) ) then
      incline(:sub) = Q_Incline%values(1,maf)
      call interpolateValues ( z_mif(1:k),        Q_Incline%values(p,maf), &
                             & tan_press(sub+1:), incline(sub+1:), &
                             & METHOD = 'L', EXTRAPOLATE='C' )
    end if

    ! ==========================================================================
    ! !!!!! This doesn't create LOS for Earth-intersecting rays correctly. !!!!!
    ! !!!!! I'm not sure what to do here because we don't know which MIF   !!!!!
    ! !!!!! zetas are for intersecting rays.  Even if we did, we don't     !!!!!
    ! !!!!! have zetas for subsurface tangent points to which to           !!!!!
    ! !!!!! interpolate them.                                              !!!!!
    ! ==========================================================================
    if ( present(Q_LOS) ) then
      do i = 1, 2 ! 1 = C, 2 = U
        ! Use the first MAF value for subsurface values, if any.
        ! ==================================================================
        ! !!!!! This is wrong.  LOS(1,:) ought to be the vector from   !!!!!
        ! !!!!! the instrument to the intersection point, and LOS(2,:) !!!!!
        ! !!!!! ought to be the vector that's a reflection of the      !!!!!
        ! !!!!! direction to the intersection, but we don't know what  !!!!!
        ! !!!!! direction that is.  Tangent_Quantities does compute    !!!!!
        ! !!!!! MIF LOS correctly, but how do we interpolate that to   !!!!!
        ! !!!!! pointing-grid LOS without subsurface MIF zetas?        !!!!!
        ! !!!!! Eventually, we could do the radiative transfer on MIF  !!!!!
        ! !!!!! rays instead of aiming at an hypothetical subsurface   !!!!!
        ! !!!!! tangent pressure.  This would affect many places,      !!!!!
        ! !!!!! setting up convolution.                                !!!!!
        ! ==================================================================
        LOS(i,:sub) = ECR_t(Q_LOS%value3(3*i-2:3*i,1,maf))
        ! Now interpolate rays that don't intersect the surface from MIF zetas
        ! to pointing-grid zetas.
        do j = 1, 3 ! XYZ
          call interpolateValues ( z_mif(1:k),        Q_LOS%value3(3*i-3+j,p,maf), &
                                 & tan_press(sub+1:), LOS(i,sub+1:)%xyz(j), &
                                 & METHOD = 'L', EXTRAPOLATE='C' )
        end do
      end do
    end if

    if ( present(Q_TanHt) ) then
      tanHt(:sub) = Q_TanHt%values(1,maf)
      call interpolateValues ( z_mif(1:k),        Q_TanHt%values(p,maf), &
                             & tan_press(sub+1:), tanHt(sub+1:), &
                             & METHOD = 'L', EXTRAPOLATE='C' )
      if ( present(tan_pt_Geod) ) then
        call interpolateValues ( z_mif(1:k),        phitan%template%geodLat(p,maf), &
                               & tan_press(sub+1:), lat(sub+1:), &
                               & METHOD = 'L', EXTRAPOLATE='C' )
        call interpolateValues ( z_mif(1:k),        phitan%template%lon(p,maf), &
                               & tan_press(sub+1:), lon(sub+1:), &
                               & METHOD = 'L', EXTRAPOLATE='C' )
        ! Subsurface values all the same.  These might also be changed to MIF
        ! some day.
        tan_pt_geod(:sub) = h_v_geod ( lon_t(lon(1)), lat(1), tanHt(1) )
        do i = 1, nlvl
          tan_pt_geod(sub+i) = h_v_geod ( lon_t(lon(sub+i)), lat(sub+i), tanHt(sub+i) )
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
! Revision 2.1  2016/11/09 00:24:27  vsnyder
! Moved from FullForwardModel
!
