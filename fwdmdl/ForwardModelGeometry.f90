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

  public :: Compute_Viewing_Plane, Height_From_Zeta
  public :: Viewing_Azimuth, Viewing_Azimuth_Qty, Viewing_Azimuth_Vec

  ! If abs(viewing_azimuth(...)) < Azimuth_Tol, assuming the viewing
  ! plane is the orbit plane.
  real, parameter, public :: Azimuth_Tol = 0.05 ! 3 degrees

  interface Viewing_Azimuth
    module procedure Viewing_Azimuth_Qty, Viewing_Azimuth_Vec
  end interface

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures ==========================

  subroutine Compute_Viewing_Plane ( Qty, ScVelECR, TpGeocAlt, GPHQuantity, &
                                   & XYZs )

  ! Compute Geocentric (!) latitude, longitude, and height for Qty.

  ! We assume that the quantity types and the units of their values have been
  ! checked, that the relationship between the magnetic field's vertical
  ! coordinate and the existence of the other quantities has been checked, that
  ! ScVelECR and TpGeocAlt are either both present or both absent, and that the
  ! relationship between the existence, number, and values of cross-angles and
  ! the existence of the other quantities has been checked.  See
  ! UsingMagneticModel in FillUtils_1.

    use Allocate_Deallocate, only: Allocate_Test
    use Constants, only: Deg2Rad, Rad2Deg
    use Cross_m, only: Cross
    use Dump_0, only: Dump
    use Geometry, only: GeodToGeocAlt, GeodToGeocLat, To_XYZ, XYZ_To_Geod
    use Intrinsic, only: L_Geocentric, L_Geodetic, L_GeodAltitude, L_Zeta
    use MLSNumerics, only: InterpolateExtrapolate_d
    use MLSStringlists, only: SwitchDetail
    use Monotone, only: Longest_Monotone_Subsequence
    use QuantityTemplates, only: CreateGeolocationFields, &
      & DestroyGeolocationFields, Dump, PointQuantityToHGrid, RT
    use Quantity_Geometry, only: XYZ
    use Rotation_m, only: Rotate_3d
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_t

    type(vectorValue_t), intent(inout) :: Qty   ! To get geolocation computed
    type(vectorValue_t), intent(in), optional :: ScVelECR  ! For SC geolocation,
                                                ! not velocity
    type(vectorValue_t), intent(in), optional :: TpGeocAlt ! Tangent geolocation,
                                                ! and altitude
    type(vectorValue_t), intent(in), optional :: GPHQuantity ! Geopotential height,
                                                ! a proxy for altitude and used
                                                ! to calculate height from zeta
    real(rt), intent(out) :: XYZs(:,:,:,:)      ! Geocentric Cartesian coordinates
                                                ! of Qty%Value4.
                                                ! 3 X surfs X instances X cross

    logical :: Across                           ! Viewing across the orbit track
    real(rt), allocatable :: Alt(:)             ! Geocentric altitude on LOS
    integer :: C                                ! Cross-track angle index
    integer, allocatable :: Dec(:)              ! Decreasing sequence indices
    integer :: Detail                           ! Dump level
    real(rt) :: Geod(3)                         ! Geodetic Lat, Lon (degrees),
                                                ! Altitude (meters)
    real(rt), allocatable :: Heights(:,:)       ! Meters, geocentric altitude
    integer, allocatable :: Inc(:)              ! Increasing sequence indices
    integer :: Inst                             ! Second subscript for lat, lon
    integer :: InstOr1                          ! Second subscript for surfs
    integer :: Me = -1                          ! String index, for tracing
    integer :: MIF                              ! MIF index
    real(rt) :: N(3)                            ! Normal to SC-TP plane
    real(rt), allocatable :: R(:,:)             ! TP rotated in TP-SC plane
    integer :: S                                ! Surface index in magnetic field
    real(rt) :: SC_XYZ(3)                       ! Spacecraft position
    real(rt) :: Sec_x                           ! Secant of cross-track angle
    integer, allocatable :: Seq(:)              ! Increasing or decreasing
                                                ! sequence indices, => Inc or Dec
    integer :: SurfOr1                          ! First subscript for Lat, Lon
    real(rt) :: TP_XYZ(3)                       ! Tangent point position

    call trace_begin ( me, 'Compute_Viewing_Plane', &
      & cond=toggle(emit) .and. levels(emit) > 1 ) ! set by -f command-line switch

    across = associated(qty%template%crossAngles)
    if ( across ) &
      & across = abs(viewing_azimuth ( tpGeocAlt, scVelECR, 1 )) > azimuth_tol

    detail = switchDetail ( switches, 'plane' )

    if ( qty%template%verticalCoordinate == l_zeta ) then
      ! Get the geodetic altitude corresponding to each zeta in the
      ! magnetic field quantity.
      ! Hopefully, we got here from FillUtils_1%UsingMagneticModel, where
      ! it is verified that if the magnetic field's vertical coordinate
      ! is zeta, then GphQuantity is provided, and that GphQuantity has
      ! the same number of instances, at the same geolocations, as the
      ! magnetic field quantity.
      call height_from_zeta ( qty, gphQuantity, heights )
    else ! the geocentric or geodetic heights are in the magnetic field quantity
      ! Eventually, this explicit allocation will not be necessary.
      ! As of 2015-04-22, ifort 15.0.2.164 requires a command-line option to
      ! make the standard-conforming automatic allocation work.
      call allocate_test ( heights, size(qty%template%surfs,1), size(qty%template%surfs,2), &
        & moduleName, 'Heights' )
      heights = qty%template%surfs ! heights allocated automatically someday
      if ( qty%template%verticalCoordinate == l_geodAltitude ) then
        ! We know that if the vertical coordinate is zeta, then tpGeocAlt is
        ! not present, and we will therefore get XYZs from the quantity template
        ! using Quantity_Geometry%XYZ.  This wants geodetic altitudes, so don't
        ! convert them to geocentric heights.
        ! Convert geodetic altitude to geocentric height
        do inst = 1, size(heights,2)
          do s = 1, size(heights,1)
            surfOr1 = merge(1,s,qty%template%stacked)
            heights(s,inst) = geodToGeocAlt ( [ qty%template%geodLat(surfOr1,inst), &
                                              & qty%template%lon(surfOr1,inst), &
                                              & heights(s,inst) ] )
          end do
        end do
      end if
    end if

    if ( detail > 2 ) call dump ( heights, name='Heights' )

    if ( present(tpGeocAlt) ) then
      ! Make the magnetic field quantity incoherent and unstacked.
      ! Re-create the latitude and longitude arrays and re-compute their values.
      qty%template%coherent = .false.
      qty%template%stacked = .false.
      call createGeolocationFields ( qty%template, qty%template%noSurfs, 'MagneticField' )

      ! Hopefully, we got here from FillUtils_1%UsingMagneticModel, where
      ! it is verified that if TpGeocAlt is provided, then so is ScVelECR,
      ! and that they have the same numbers of instances and surfaces.
      call allocate_test ( Alt, tpGeocAlt%template%noSurfs, moduleName, 'Alt' )
      call allocate_test ( r, 3, tpGeocAlt%template%noSurfs, moduleName, 'R' )

      do inst = 1, qty%template%NoInstances
        instOr1 = merge(1,inst,qty%template%coherent)
        ! Use only the monotone part of the tangent-point altitudes
        call longest_monotone_subsequence ( tpGeocAlt%value3(1,:,inst), inc )
        call longest_monotone_subsequence ( tpGeocAlt%value3(1,:,inst), dec, -1 )
        if ( size(inc) >= size(dec) ) then
          call move_alloc ( inc, seq )
          deallocate ( dec )
        else
          call move_alloc ( dec, seq )
          deallocate ( inc )
        end if
        do c = 1, qty%template%noCrossTrack
          if ( associated(qty%template%crossAngles) ) then
            sec_x = 1.0 / cos(qty%template%crossAngles(c)*deg2rad)
          else
            sec_x = 1.0
          end if
          do MIF = 1, size(seq)
            ! Get SC and TP positions for seq(MIF) in geocentric Cartesian
            ! coordinates. Notice that SC_XYZ, TP_XYZ, and N do not depend upon C. 
            ! We could exchange the loop order, and hoist them out of the inner
            ! loop, at the expense of making Alt and R two-dimensional arrays.
            ! The result of geodToGeocLat is in RADIANS!  To_XYZ wants DEGREES!
            sc_xyz = to_xyz ( geodToGeocLat(ScVelECR%template%geodLat(seq(MIF),inst))*rad2Deg, &
                            & ScVelECR%template%lon(seq(MIF),inst) )
            tp_xyz = to_xyz ( geodToGeocLat(tpGeocAlt%template%geodLat(seq(MIF),inst))*rad2Deg, &
                            & tpGeocAlt%template%lon(seq(MIF),inst) )
            ! Compute normal to SC-TP plane (doesn't need to be unit normal).
            n = cross(sc_xyz, tp_xyz )
            ! Rotate TP toward or away from SC in the TP-SC plane (i.e., about N)
            ! by CrossAngles(c) degrees, giving R in geocentric Cartesian
            ! coordinates.  TP_xyz is a unit vector, so R is also.
            if ( associated(qty%template%crossAngles) ) then
              call rotate_3d ( tp_xyz, qty%template%crossAngles(c) * deg2rad, n, &
                             & r(:3,seq(MIF)) )
            else ! The cross-angle is assumed to be zero if there isn't one.
              r(:3,seq(MIF)) = tp_xyz
            end if
            ! Compute the geocentric altitude at the rotated position, and
            ! extend R to that length.
            alt(seq(MIF)) = tpGeocAlt%value3(1,seq(MIF),inst) * sec_x
            r(:3,seq(MIF)) = r(:3,seq(MIF)) * alt(seq(MIF))
          end do ! MIF
          ! Interpolate Cartesian ECR coordinates, in meters, at the rotated
          ! position using linear interpolation and extrapolation on Altitude
          ! at the rotated position to the altitude in Qty.  Only use the
          ! monotone part of the rotated-position altitudes.  Extrapolate
          ! outside the range of Alt using the average slope of R(i,:).
          ! This might not be exactly correct in geodetic coordinates, but it
          ! ought to be close enough.
          call interpolateExtrapolate_d ( alt(seq), r(:,seq), &
            & heights(:,instOr1), xyzs(:,:,instOr1,c), second=.true. )
          if ( qty%template%verticalCoordinate == l_geodAltitude .or. &
             & qty%template%verticalCoordinate == l_zeta ) then
            qty%template%latitudeCoordinate = l_geodetic
            ! Get geodetic latitude
            do s = 1, merge(1,qty%template%noSurfs,qty%template%stacked)
              geod = xyz_to_geod ( xyzs(:,s,instOr1,c) )
              call qty%template%putLat3 ( s, inst, c, geod(1) * rad2deg )
            end do
          else ! qty%template%verticalCoordinate == l_geocAltitude
            ! Get geocentric latitude
            qty%template%latitudeCoordinate = l_geocentric
            do s = 1, merge(1,qty%template%noSurfs,qty%template%stacked)
              call qty%template%putLat3 ( s, inst, c, &
                & asin(xyzs(3,s,instOr1,c) / norm2(xyzs(:,s,instOr1,c),1)) * rad2deg )
            end do
          end if
          ! Get longitude
          do s = 1, size(qty%template%lon,1)
            call qty%template%putLon3 ( s, inst, c, &
              & atan2(xyzs(2,s,instOr1,c),xyzs(1,s,instOr1,c)) * rad2deg )
          end do
        end do ! c = 1, qty%template%noCrossTrack
        deallocate ( seq )
      end do ! inst = 1, qty%template%NoInstances
    else ! Use the geolocations already in the magnetic field quantity.
      ! Since we have no TpGeocAlt and presumably no ScVelECR, the viewing
      ! plane must be the orbit plane.  It's easy to compute XYZ from the
      ! geolocations in the magnetic field quantity template.
      call destroyGeolocationFields ( qty%template )
      qty%template%stacked = .true.
      qty%template%coherent = .true.
      call pointQuantityToHGrid ( qty%template )
      xyzs(:,:,:,1) = xyz(qty%template, heights)
    end if

    if ( detail > 2 ) call dump ( xyzs, name='XYZs' )

    if ( detail > -1 ) call dump ( qty%template, details=detail, &
      & what='In Compute_Viewing_Plane' )

    call trace_end ( 'Compute_Viewing_Plane', &
      & cond=toggle(emit) .and. levels(emit) > 1 )

  end subroutine Compute_Viewing_Plane

  subroutine Height_From_Zeta ( Qty, GPHQuantity, Heights )
    ! Compute an array of geopotential heights in meters that correspond to
    ! Qty%Surfs in zeta.  These are roughly similar to geodetic altitudes, but
    ! depend upon the temperature profile by way of zeta's dependence upon
    ! temperature.  Then convert them to geodetic altitudes.

    use Allocate_Deallocate, only: Allocate_Test
    use Intrinsic, only: L_GeodAltitude, L_Zeta
    use Heights_Module, only: GPH_to_Geom
    use Intrinsic, only: L_Zeta
    use MLSMessageModule, only : MLSMessage, MLSMSG_Error
    use MLSNumerics, only: InterpolateValues
    use QuantityTemplates, only: RT
    use Toggles, only: Emit, Levels, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_t

    type (VectorValue_T), intent(in) :: Qty
    type (VectorValue_T), intent(in) :: GPHQuantity
    real(rt), intent(out), allocatable :: Heights(:,:) ! meters, geodetic
      ! altitude, same shape as Qty%Surfs (NoSurfs x NoInstances or 1)

    integer :: IQ, IG, JQ, NQ, NG, SurfOr1
    integer :: Me = -1                          ! String index, for tracing

    call trace_begin ( me, 'Height_From_Zeta', &
      & cond=toggle(emit) .and. levels(emit) > 1 ) ! set by -f command-line switch

    if ( qty%template%verticalCoordinate /= l_zeta ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Height requested for zeta, but vertical coordinate is not zeta' )
      return ! Hopefully, we don't get here
    end if

    if ( gphQuantity%template%verticalCoordinate /= l_zeta ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'GPHQuantity vertical coordinate is not zeta' )
      return ! Hopefully, we don't get here
    end if

    nq = size(qty%template%surfs,1)
    ng = size(gphQuantity%template%surfs,1)

    if ( nq > 1 .and. ng > 1 .and. ng /= nq .or. nq == 1 .and. ng > 1 ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Numbers of instances of Qty and GHPQuantity not compatible ' // &
        & 'in Height_From_Zeta' )
      return ! Hopefully, we don't get here
    end if

    call allocate_test ( heights, nq, size(qty%template%surfs,2), moduleName, &
                       & 'Heights' )

    do iq = 1, nq
      ig = merge(iq,1,ng>1)
      call interpolateValues ( gphQuantity%template%surfs(ig,:), &
                             & gphQuantity%values(ig,:), &
                             & qty%template%surfs(iq,:), heights(iq,:), &
                             & method='L' )
    end do

    if ( qty%template%verticalCoordinate == l_geodAltitude ) then
      do jq = 1, size(qty%template%surfs,2)
        surfOr1 = merge(1,jq,qty%template%stacked)
        do iq = 1, nq
          heights(iq,jq) = gph_to_geom ( heights(iq,jq), &
                                       & lat_geod=qty%template%geodLat(surfOr1,jq) )
        end do
      end do
    else ! qty%template%verticalCoordinate == l_geocAltitude
      do jq = 1, size(qty%template%surfs,2)
        surfOr1 = merge(1,jq,qty%template%stacked)
        do iq = 1, nq
          heights(iq,jq) = gph_to_geom ( heights(iq,jq), &
                                       & lat_geoc=qty%template%geodLat(surfOr1,jq) )
        end do
      end do
    end if

    call trace_end ( 'Height_From_Zeta', &
      & cond=toggle(emit) .and. levels(emit) > 1 )

  end subroutine Height_From_Zeta

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
    use MLSMessageModule, only: MLSMessage, MLSMSG_Warning
    use MoreMessage, only: MLSMessage
    use VectorsModule, only: VectorValue_t

    type(vectorValue_t), intent(in) :: TpGeocAlt   ! Tangent point, for geolocation
    type(vectorValue_t), intent(in) :: ScVelECR    ! For geolocation and velocity
    integer, intent(in), optional :: MIF
    integer, intent(in), optional :: MAF
    real(rt) :: Viewing_Azimuth_Rad

    integer :: MyMAF, MyMIF
    real(rt) :: N1, N2                       ! Norm of a 3-vector
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
    n1 = norm2(sc_xyz)
    scVel_XYZ = scVelECR%value4(:,myMIF,myMAF,1)
    n2 = norm2(scVel_XYZ)
    if ( n1 < 1.0 .or. n2 < 1.0 ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "Either position or velocity in %S quantity apparently not filled.", &
        & scVelECR%template%name )
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "Zero used for viewing azimuth.  Hope that's OK." )
      viewing_azimuth_rad = 0.0
      return
    end if
    scVel_XYZ = scVel_XYZ / n2
    tp_xyz = to_xyz ( geodToGeocLat(tpGeocAlt%template%geodLat(myMIF,myMAF))*rad2Deg, &
                    & tpGeocAlt%template%lon(myMIF,myMAF) )

    ! Get unit normals to the Specraft Velocity X Spacecraft Position and
    ! Spacecraft Position X Tangent point planes.
    sc_n = cross(scVel_xyz,sc_xyz,norm=.true.)
    tp_n = cross(sc_xyz,tp_xyz,norm=.true.)

    ! Compute the angle (Radians!) between the planes' normals
    viewing_azimuth_rad = acos(abs(dot_product(sc_n,tp_n)))

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

    type(vectorValue_t), pointer :: SCVelECR ! Spacecraft Velocity, ECR
    type(vectorValue_t), pointer :: TpGeocAlt ! Tangent point

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
! Revision 2.5  2015/06/04 03:13:40  vsnyder
! Make Surfs component of quantity template allocatable
!
! Revision 2.4  2015/06/03 00:00:32  vsnyder
! Use type-bound procedures to do rank-3 reference and update for latitude
! and longitude, instead of using rank-remapped pointers.
!
! Revision 2.3  2015/05/28 23:08:58  vsnyder
! Eliminate RefMIF, check that ScVelECR geolocation is filled
!
! Revision 2.2  2015/04/29 02:08:39  vsnyder
! Don't convert heights to geocentric if vertical is zeta
!
! Revision 2.1  2015/03/29 18:48:46  vsnyder
! Initial commit
!
