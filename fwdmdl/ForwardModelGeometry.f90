! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module ForwardModelGeometry

  ! This module contains procedures to compute various geometry-related
  ! quantities for the forward model.

  implicit none
  private

  public :: Compute_Viewing_Plane
  public :: Get_Quantity_XYZ
  public :: Viewing_Azimuth, Viewing_Azimuth_Qty, Viewing_Azimuth_Vec

  ! If abs(viewing_azimuth(...)) < Azimuth_Tol, assuming the viewing
  ! plane is the orbit plane.
  real, parameter, public :: Azimuth_Tol = 0.01

  interface Viewing_Azimuth
    module procedure Viewing_Azimuth_Qty, Viewing_Azimuth_Vec
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  subroutine Compute_Viewing_Plane ( Qty, TpGeocAlt, ScVelECR, FirstMAF, LastMAF )
    ! Compute Geocentric (!) latitude, longitude, and height for Qty for the
    ! specified MAFs (or all MAFs).  If Qty's NoCrossTrack field is <= 1, or the
    ! Viewing_Azimuth is zero (or close to zero), nothing is done here. Qty is
    ! assumed to be stacked and coherent, i.e., the GeodLat and Lon components
    ! will be the same for all surfaces in each instance, and they will be
    ! assumed to have been computed in the orbit plane by their HGrid.

    ! Otherwise, Qty cannot be stacked, should be coherent, and the GeodLat
    ! and Lon values are computed here for each surface and instance.
    
    ! In either case, the Surfs component is specified by a VGrid with a
    ! vertical coordinate of GeocAltitude (in meters).

    use Constants, only: Deg2Rad, Rad2Deg
    use Cross_m, only: Cross
    use Geometry, only: GeodToGeocLat, To_XYZ, XYZ_To_Geod
    use Intrinsic, only: L_Geocentric, L_Geodetic, L_GeodAltitude
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSNumerics, only: InterpolateValues
    use Monotone, only: Get_Monotone
    use QuantityTemplates, only: QuantitiesAreCompatible, RT
    use Rotation_m, only: Rotate_3d
    use VectorsModule, only: VectorValue_t

    type(vectorValue_t), intent(inout) :: Qty    ! To get geolocation computed
    type(vectorValue_t), intent(in) :: TpGeocAlt ! Tangent geolocation, not alt
    type(vectorValue_t), intent(in) :: ScVelECR  ! For SC geolocation, not vel
    integer, intent(in), optional :: FirstMAF    ! Default 1
    integer, intent(in), optional :: LastMAF     ! Default Qty%Template%NoInstances;
                                                 ! not used if FirstMAF is absent.

    real(rt) :: Alt(min(tpGeocAlt%template%noSurfs, &
                   & scVelECR%template%noSurfs))! Geocentric altitude on LOS
    real(rt) :: Geod(3)                         ! Lat, Lon, Height from XYZ_To_Geod
    integer :: I                                ! Cross-track angle index
    real(rt) :: Lat(size(alt))                  ! Geocentric (!) latitude on LOS
    real(rt) :: Lon(size(alt))                  ! Longitude on LOS
    integer :: MAF, MAF1, MAFn                  ! MAF indices
    integer :: MIF                              ! MIF index
    integer :: M1, M2                           ! Equivalenced to StartStop
    real(rt) :: N(3)                            ! Normal to SC-TP plane
    real(rt) :: R(3)                            ! TP rotated in TP-SC plane
    integer :: S                                ! Select which way to pack Alt
    real(rt) :: SC_XYZ(3)                       ! Spacecraft position
    real(rt) :: Sec_x                           ! Secant of cross-track angle
    integer :: StartStop(2)                     ! Ends of monotone altitude seq
    real(rt) :: TP_XYZ(3)                       ! Tangent point position

    equivalence ( StartStop(1), M1 ), ( StartStop(2), M2 )

    if ( qty%template%noCrossTrack <= 1 ) return ! Nothing needs doing

    MAF1 = 1; MAFn = qty%template%NoInstances
    if ( present(firstMAF) ) then
      MAF1 = firstMAF
      if ( present(lastMAF) ) MAFn = max(MAF1,lastMAF)
    end if

    if ( abs(viewing_azimuth ( tpGeocAlt, scVelECR, 1, MAF1 )) < azimuth_tol ) &
      & return ! Nothing needs doing

    if ( qty%template%stacked ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        "Cross-track viewing quantity cannot be stacked" )
      return ! Hopefully we don't get here
    end if

    if ( .not. quantitiesAreCompatible ( tpGeocAlt%template, ScVelECR%template, &
                                       & differentTypeOK=.true. ) ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        "Quantities are not compatible in Compute_Viewing_Plane" )
      return ! Hopefully we don't get here
    end if

    do MAF = MAF1, MAFn
      do i = 1, qty%template%noCrossTrack
        sec_x = 1.0 / cos(qty%template%crossAngles(i)*deg2rad)
        do MIF = 1, size(alt,1)
          ! Get SC and TP positions for MIF MIF
          ! Notice that SC_XYZ, TP_XYZ, and N do not depend upon I.  We could
          ! exchange the loop order, and hoist them out of the inner loop, at
          ! the expense of making Alt, Lat, and Lon two-dimensional arrays.
          ! The result of geodToGeocLat is in RADIANS!  To_XYZ wants DEGREES!
          sc_xyz = to_xyz ( geodToGeocLat(ScVelECR%template%geodLat(MIF,MAF))*rad2Deg, &
                          & ScVelECR%template%lon(MIF,MAF) )
          tp_xyz = to_xyz ( geodToGeocLat(tpGeocAlt%template%geodLat(MIF,MAF))*rad2Deg, &
                          & tpGeocAlt%template%lon(MIF,MAF) )
          ! Compute normal to SC-TP plane (doesn't need to be unit normal).
          ! The tangent-point position vector will be rotated in the SC-TP plane
          ! (about this vector) toward or away from the spacecraft position.
          n = cross(sc_xyz, tp_xyz )
          ! Rotate TP toward SC in TP-SC plane (about N) by CrossAngles(i) degrees
          call rotate_3d ( tp_xyz, qty%template%crossAngles(i) * deg2rad, n, r )
          lat(MIF) = asin(r(3)) * rad2deg
          lon(MIF) = atan2(r(2),r(1)) * rad2deg
          ! Compute the geocentric altitude at the rotated position.
          alt(MIF) = tpGeocAlt%value3(1,MIF,MAF) * sec_x
          ! Convert to geodetic coordinates if necessary
          if ( qty%template%verticalCoordinate == l_geodAltitude ) then
            geod = xyz_to_geod ( r * alt(MIF) )
            lat(MIF) = geod(1) * rad2deg
            lon(MIF) = geod(2) * rad2deg
            alt(MIF) = geod(3)
          end if
        end do
        ! Interpolate Latitude and Longitude at the rotated position using
        ! linear interpolation on Altitude at the rotated position to altitude
        ! in Qty.  Only use the monotone part of the altitude.
        startStop = get_monotone ( alt )
        call interpolateValues ( alt(m1:m2), lat(m1:m2), &
          & qty%template%surfs(:,1), qty%template%geodLat(:,i), &
          & method='L', extrapolate='C' )
        call interpolateValues ( alt(m1:m2), lon(m1:m2), &
          & qty%template%surfs(:,1), qty%template%lon(:,i), &
          & method='L', extrapolate='C' )
      end do
    end do
    qty%template%latitudeCoordinate = l_geocentric
    if ( qty%template%verticalCoordinate == l_geodAltitude ) &
      & qty%template%latitudeCoordinate = l_geodetic
  end subroutine Compute_Viewing_Plane

  function Get_Quantity_XYZ ( FwdModelIn, FwdModelExtra, QtyType, Config, &
    & MIF, MAF, InstrumentModule ) result (XYZ)
    ! Get the geocentric coordinates, as an unit vector in ECR, for the
    ! specified quantity and specified MIF (if provided) or MIF #1, and the
    ! specified MAF (if provided) or MAF #1.
    use Constants, only: Rad2Deg
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Geometry, only: GeodToGeocLat, To_XYZ
    use QuantityTemplates, only: RT
    use VectorsModule, only: Vector_t, VectorValue_t

    type(vector_t), intent(in) :: FwdModelIn
    type(vector_t), intent(in), optional :: FwdModelExtra
    integer, intent(in) :: QtyType
    type(forwardModelConfig_t), intent(in), optional :: Config
    integer, intent(in), optional :: MIF
    integer, intent(in), optional :: MAF
    integer, intent(in), optional :: InstrumentModule
    real(rt) :: XYZ(3)

    integer :: MyMAF, MyMIF
    type(vectorValue_t), pointer :: Qty

    myMAF = 1
    if ( present(MAF) ) myMAF = MAF
    myMIF = 1
    if ( present(MIF) ) myMIF = MIF

    qty => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
        & quantityType=qtyType, config=config, &
        & instrumentModule=InstrumentModule )
    xyz = to_xyz ( geodToGeocLat(qty%template%geodLat(myMIF,myMAF))*rad2Deg, &
                 & qty%template%lon(myMIF,myMAF) )

  end function Get_Quantity_XYZ

  function Viewing_Azimuth_Qty ( TpGeocAlt, ScVelECR, MIF, MAF ) &
    & result ( Viewing_Azimuth_Rad )
    ! Compute the angle (in Radians!) between the normal to the plane defined
    ! by the spacecraft position and velocity vectors, and the normal to the 
    ! plane defined by PTan and the spacecraft position, for the specified
    ! MIF (if provided) or MIF #1, and the the specified MAF (if provided)
    ! or MAF #1.
    use Constants, only: Rad2Deg
    use Cross_m, only: Cross
    use Geometry, only: GeodToGeocLat, To_XYZ
    use QuantityTemplates, only: RT
    use VectorsModule, only: VectorValue_t

    type(vectorValue_t), intent(in) :: TpGeocAlt   ! Tangent point, for geolocation
    type(vectorValue_t), intent(in) :: ScVelECR    ! For geolocation and velocity
    integer, intent(in), optional :: MIF
    integer, intent(in), optional :: MAF
    real(rt) :: Viewing_Azimuth_Rad

    integer :: MyMAF, MyMIF
    real(rt) :: SCVel_XYZ(3)                 ! Unit vector in velocity direction
    real(rt) :: SC_N(3)                      ! Normal to Velocity x Position plane
    real(rt) :: SC_XYZ(3)                    ! Spacecraft position, ECR
    real(rt) :: TP_N(3)                      ! Normal to Spacecraft X TP plane
    real(rt) :: TP_XYZ(3)                    ! Tangent point, ECR

    myMAF = 1
    if ( present(MAF) ) myMAF = MAF
    myMIF = 1
    if ( present(MIF) ) myMIF = MIF

    ! Get position and velocity vectors
    sc_xyz = to_xyz ( geodToGeocLat(scVelECR%template%geodLat(myMIF,myMAF))*rad2Deg, &
                    & scVelECR%template%lon(myMIF,myMAF) )
    scVel_XYZ = scVelECR%value4(:,myMIF,myMAF,1)
    scVel_XYZ = scVel_XYZ / norm2(scVel_XYZ)
    tp_xyz = to_xyz ( geodToGeocLat(tpGeocAlt%template%geodLat(myMIF,myMAF))*rad2Deg, &
                    & tpGeocAlt%template%lon(myMIF,myMAF) )

    ! Get unit normals to the Specraft Velocity X Spacecraft Position and
    ! Spacecraft Position X Tangent point planes.
    sc_n = cross(scVel_xyz,sc_xyz,norm=.true.)
    tp_n = cross(sc_xyz,tp_xyz,norm=.true.)

    ! Compute the angle (Radians!) between the planes' normals
    viewing_azimuth_rad = acos(dot_product(sc_n,tp_n))

  end function Viewing_Azimuth_Qty

  function Viewing_Azimuth_Vec ( FwdModelIn, FwdModelExtra, Config, MIF, MAF ) &
    & result ( Viewing_Azimuth_Rad )
    ! Compute the angle (in Radians!) between the normal to the plane defined
    ! by the spacecraft position and velocity vectors, and the normal to the 
    ! plane defined by PTan and the spacecraft position, for the specified
    ! MIF (if provided) or MIF #1, and the the specified MAF (if provided)
    ! or MAF #1.
    use ForwardModelConfig, only: ForwardModelConfig_T
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Intrinsic, only: L_ScVelECR, L_TngtGeocAlt
    use QuantityTemplates, only: RT
    use VectorsModule, only: Vector_t, VectorValue_t

    type(vector_t), intent(in) :: FwdModelIn
    type(vector_t), intent(in), optional :: FwdModelExtra
    type(forwardModelConfig_t), intent(in), optional :: Config
    integer, intent(in), optional :: MIF
    integer, intent(in), optional :: MAF
    real(rt) :: Viewing_Azimuth_Rad

    integer :: MyMAF, MyMIF
    type(vectorValue_t), pointer :: SCVelECR ! Spacecraft Velocity, ECR
    type(vectorValue_t), pointer :: TpGeocAlt ! Tangent point

    myMAF = 1
    if ( present(MAF) ) myMAF = MAF
    myMIF = 1
    if ( present(MIF) ) myMIF = MIF

    ! Get quantities
    scVelECR => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
             & quantityType=l_scVelECR, config=config )
    tpGeocAlt => GetQuantityForForwardModel ( fwdModelIn, fwdModelExtra, &
             & quantityType=l_tngtGeocAlt, config=config )

    viewing_azimuth_rad = viewing_azimuth ( tpGeocAlt, scVelECR, MIF, MAF )

  end function Viewing_Azimuth_Vec

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module ForwardModelGeometry

! $Log$
! Revision 2.1  2015/03/29 18:48:46  vsnyder
! Initial commit
!
