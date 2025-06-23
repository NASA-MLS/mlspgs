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
                                   & XYZs, Regular, ReferenceMIF )

  ! Compute Geocentric (!) latitude, longitude, and height for Qty.

  ! We assume that the quantity types and the units of their values have been
  ! checked, that the relationship between the Qty's vertical coordinate and
  ! the existence of the other quantities has been checked, that ScVelECR and
  ! TpGeocAlt are either both present or both absent, and that the relationship
  ! between the existence, number, and values of cross-angles and the existence
  ! of the other quantities has been checked.  See UsingMagneticModel in
  ! FillUtils_1.

    use Allocate_Deallocate, only: Allocate_Test
    use Constants, only: Deg2Rad, Rad2Deg
    use Cross_m, only: Cross
    use Dump_0, only: Dump
    use Geometry, only: GeodToGeocAlt, GeodToGeocLat, To_XYZ, XYZ_To_Geod
    use Intrinsic, only: L_Geocentric, L_Geodetic, L_GeodAltitude, L_Zeta
    use MLSNumerics, only: InterpolateExtrapolate
    use MLSStringlists, only: SwitchDetail
    use Monotone, only: Longest_Monotone_Subsequence
    use Output_m, only: Output
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
    logical, intent(in), optional :: Regular    ! coherent and stacked
    integer, intent(in), optional :: ReferenceMIF

    logical :: Across                     ! Viewing across the orbit track
    real(rt), allocatable :: Alt(:)       ! Geocentric altitude on LOS
    integer :: C                          ! Cross-track angle index
    integer, allocatable :: Dec(:)        ! Decreasing sequence indices
    integer :: Detail                     ! Dump level
    real(rt) :: Geod(3)                   ! Geodetic Lat, Lon (degrees),
                                          ! Altitude (meters)
    real(rt), allocatable :: Heights(:,:) ! Meters, geocentric altitude
    integer, allocatable :: Inc(:)        ! Increasing sequence indices
    integer :: Inst                       ! Second subscript for lat, lon
    integer :: InstOr1                    ! Second subscript for surfs
    integer :: Me = -1                    ! String index, for tracing
    integer :: MIF                        ! MIF index
    logical :: MyRegular                  ! Coherent and stacked
    real(rt) :: N(3)                      ! Normal to SC-TP plane
    integer :: P                          ! Index in crossAngles at which to
                                          ! compute phi -- the one with the
                                          ! smallest absolute value.
    real(rt), allocatable :: In_Phi(:,:)  ! from qty%template%phi
    real(rt), allocatable :: R(:,:)       ! TP rotated in TP-SC plane
    integer :: S                          ! Surface index in Qty
    real(rt) :: SC_XYZ(3)                 ! Spacecraft position
    real(rt) :: Sec_x                     ! Secant of cross-track angle
    integer, allocatable :: Seq(:)        ! Increasing or decreasing
                                          ! sequence indices, => Inc or Dec
    integer :: SurfOr1                    ! First subscript for Lat, Lon
    real(rt) :: TP_XYZ(3)                 ! Tangent point position
    real(rt) :: V(3)                      ! ECR Geolocation

    call trace_begin ( me, 'Compute_Viewing_Plane', &
      & cond=toggle(emit) .and. levels(emit) > 1 ) ! set by -f command-line switch

    myRegular = .false.
    if ( present(regular) ) myRegular = regular

    across = any(qty%template%crossAngles /= 0.0)
    if ( across .and. present(tpGeocAlt) ) &
        & across = abs(viewing_azimuth ( tpGeocAlt, scVelECR, 1 )) > azimuth_tol

    detail = switchDetail ( switches, 'plane' )

    if ( qty%template%verticalCoordinate == l_zeta ) then
      ! Get the geodetic altitude corresponding to each zeta in Qty.
      ! Hopefully, we got here from FillUtils_1%UsingMagneticModel, where it
      ! is verified that if Qty's vertical coordinate is zeta, that
      ! GphQuantity is provided, and that GphQuantity has the same number of
      ! instances, at the same geolocations, as Qty.
      call height_from_zeta ( qty, gphQuantity, heights )
    else ! the geocentric or geodetic heights are in Qty
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

    if ( detail > 3 ) then
      call dump ( heights, name='Heights' )
      if ( detail > 4 ) call dump ( qty%template, details=detail, &
        & what='Early in Compute_Viewing_Plane' )
    end if

    if ( present(tpGeocAlt) .and. .not. myRegular ) then
      ! Make Qty incoherent and unstacked.
      ! Re-create the latitude and longitude arrays and re-compute their values.
      qty%template%coherent = .false.
      qty%template%stacked = .false.
      call createGeolocationFields ( qty%template, qty%template%noSurfs, 'MagneticField' )

      ! Re-allocate phi to (noSurfs, noInstances).  If Qty was originally
      ! coherent, it was (noSurfs,1).  If it was stacked and coherent, it was
      ! (1,1).
      call allocate_test ( qty%template%phi, qty%template%noSurfs, &
        & qty%template%noInstances, moduleName, 'Phi' )

      ! Hopefully, we got here from FillUtils_1%UsingMagneticModel, where
      ! it is verified that if TpGeocAlt is provided, then so is ScVelECR,
      ! and that they have the same numbers of instances and surfaces.
      call allocate_test ( Alt, tpGeocAlt%template%noSurfs, moduleName, 'Alt' )
      call allocate_test ( r, 3, tpGeocAlt%template%noSurfs, moduleName, 'R' )

      p = minloc(abs(qty%template%crossAngles),1)
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
          sec_x = 1.0 / cos(qty%template%crossAngles(c)*deg2rad)
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
            if ( across ) then
              call rotate_3d ( tp_xyz, qty%template%crossAngles(c) * deg2rad, n, &
                             & r(:3,seq(MIF)) )
            else ! The cross-angle is zero or viewing is in the orbit plane.
              r(:3,seq(MIF)) = tp_xyz
            end if
            ! Compute the geocentric altitude at the rotated position, and
            ! extend R to that length.
            alt(seq(MIF)) = tpGeocAlt%value3(1,seq(MIF),inst) * sec_x
            r(:3,seq(MIF)) = r(:3,seq(MIF)) * alt(seq(MIF))
          end do ! MIF
          ! Compute Phi for the cross-angle nearest zero by interpolating
          ! phi from tpGeocAlt in height.  Compute it only for one
          ! instance if the quantity is coherent.
          if ( c == p .and. instOr1 == inst ) then
            call interpolateExtrapolate ( alt(seq), &
              & tpGeocAlt%template%phi(seq,instOr1), &
              & heights(:,instOr1), qty%template%phi(:,instOr1), &
              & second=.false. )
          end if
          ! Interpolate Cartesian ECR coordinates, in meters, at the rotated
          ! position using linear interpolation and extrapolation on Altitude
          ! at the rotated position to the altitude in Qty.  Only use the
          ! monotone part of the rotated-position altitudes.  Extrapolate
          ! outside the range of Alt using the average slope of R(i,:).
          ! This might not be exactly correct in geodetic coordinates, but it
          ! ought to be close enough.
          call interpolateExtrapolate ( alt(seq), r(:,seq), &
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
          if ( detail > 3 ) then
            call output ( c, before='R(' )
            call dump ( r(:,seq), name=')' )
          end if
        end do ! c = 1, qty%template%noCrossTrack
        deallocate ( seq )
      end do ! inst = 1, qty%template%NoInstances
    else if ( present(tpGeocAlt) .and. myRegular ) then
      ! Create a 3D stacked and coherent (i.e., regular rectangular) grid.
      ! It's then easy to compute XYZ from the geolocations in Qty's template.
      call destroyGeolocationFields ( qty%template )
      qty%template%stacked = .true.
      qty%template%coherent = .true.
      if ( .not. across ) then
        ! Viewing plane is the orbit plane.  Assume no cross angles.
        call pointQuantityToHGrid ( qty%template )
        xyzs(:,:,:,1) = xyz(qty%template, heights)
      else
        call allocate_test ( in_phi, 1, qty%template%noInstances, &
          & moduleName, 'In_Phi' )
        in_phi = qty%template%phi
        call createGeolocationFields ( qty%template, 1, 'Qty', Phi=.true. )
        MIF = tpGeocAlt%template%noSurfs / 2
        if ( present(referenceMIF) ) MIF = referenceMIF
        do inst = 1, qty%template%NoInstances
          ! Compute the normal to the plane defined by tpGeocAlt and scVelECR
          ! at MIF
          if ( detail > 2 ) then
            call output ( MIF, before='Height at reference MIF ' )
            call output (tpGeocAlt%template%surfs(MIF,inst), before=' = ', advance='yes' )
          end if
          sc_xyz = to_xyz ( ScVelECR%template%geocLat(MIF,inst), &
                          & ScVelECR%template%lon(MIF,inst) )
          tp_xyz = to_xyz ( tpGeocAlt%template%geocLat(MIF,inst), &
                          & tpGeocAlt%template%lon(MIF,inst) )
          tp_xyz = tp_xyz / norm2(tp_xyz)
          ! Compute unit normal to SC-TP plane.  The ECR coordinate of one node is
          ! is (n(2), -n(1), 0).  Use that node as the zero for phi.
          n = cross ( sc_xyz, tp_xyz, norm=.true. )
          do c = 1, qty%template%noCrossTrack
            ! Compute geolocations and a coordinate in the viewing plane that
            ! is an angle at the center of the earth between the above-computed
            ! node and the tangent point position, rotated in the plane
            ! containing the tangent point and spacecraft at the reference MIF.
            call rotate_3d ( tp_xyz, qty%template%crossAngles(c) * deg2rad, n, v )
            call qty%template%putGeocLat3 ( 1, inst, c, asin(v(3)) * rad2deg )
            call qty%template%putLon3 ( 1, inst, c, atan2(v(2),v(1)) * rad2deg )
!             Don't do this!  It measures phi along the great circle containing
!             the tangent point and spacecraft, from one of its equator crossings,
!             but the forward model wants phi to come from the hGrid.
!             call qty%template%putPhi3 ( 1, inst, c, &
!               & acos(dot_product(v(1:2),[n(2),-n(1)])) * rad2deg )
            call qty%template%putPhi3 ( 1, inst, c, &
              & in_phi(1,inst) + qty%template%crossAngles(c) ) 
            do s = 1, size(heights,1)
              xyzs(:,s,inst,c) = xyz(qty%template, s, inst, heights(s,inst),c )
            end do
          end do ! c = 1, qty%template%noCrossTrack
        end do ! inst = 1, qty%template%NoInstances
      end if
      ! Use the geolocations in Qty. 
    else ! Use the geolocations already in Qty.
      ! Since we have no TpGeocAlt and presumably no ScVelECR, the viewing
      ! plane must be the orbit plane.  It's easy to compute XYZ from the
      ! geolocations in Qty.
      call destroyGeolocationFields ( qty%template )
      qty%template%stacked = .true.
      qty%template%coherent = .true.
      call pointQuantityToHGrid ( qty%template )
      xyzs(:,:,:,1) = xyz(qty%template, heights)
    end if

    if ( detail > 1 ) call dump ( xyzs, name='XYZs' )

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
    use MLSStringLists, only: SwitchDetail
    use Output_m, only: Output
    use QuantityTemplates, only: RT
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    use MoreMessage, only: MLSMessage
    use Toggles, only: Switches
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

    if ( myMIF < 1 .or. myMIF > size(scVelECR%template%geodLat,1) ) &
      call MLSMessage ( MLSMSG_Error, moduleName, 'MIF out of range' )
    if ( myMAF < 1 .or. myMAF > size(scVelECR%template%geodLat,2) ) &
      call MLSMessage ( MLSMSG_Error, moduleName, 'MAF out of range' )

    ! Get position and velocity vectors
    sc_xyz = to_xyz ( geodToGeocLat(scVelECR%template%geodLat(myMIF,myMAF))*rad2Deg, &
                    & scVelECR%template%lon(myMIF,myMAF) )
    n1 = norm2(sc_xyz)
    scVel_XYZ = scVelECR%value4(:,myMIF,myMAF,1)
    n2 = norm2(scVel_XYZ)
    if ( n1 < 0.5 .or. n2 < 0.5 ) then
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "Either position or velocity in %S quantity apparently not filled.", &
        & scVelECR%template%name )
      call MLSMessage ( MLSMSG_Warning, moduleName, &
        & "Zero used for viewing azimuth.  Hope that's OK." )
      viewing_azimuth_rad = 0.0
      go to 9
    end if
    ! sc_xyz is already a unit vector.
    scVel_XYZ = scVel_XYZ / n2
    tp_xyz = to_xyz ( geodToGeocLat(tpGeocAlt%template%geodLat(myMIF,myMAF))*rad2Deg, &
                    & tpGeocAlt%template%lon(myMIF,myMAF) )

    ! Get unit normals to the Specraft Velocity X Spacecraft Position and
    ! Spacecraft Position X Tangent point planes.
    sc_n = cross(scVel_xyz,sc_xyz,norm=.true.)
    tp_n = cross(sc_xyz,tp_xyz,norm=.true.)

    ! Compute the angle (Radians!) between the planes' normals
    viewing_azimuth_rad = acos(abs(dot_product(sc_n,tp_n)))

9   continue
    if ( switchDetail(switches,'Azimuth') > -1 ) &
      & call output ( rad2deg*viewing_azimuth_rad, before='Azimuth ', &
        & after=' degrees', advance='yes' )

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
! Revision 2.11  2015/10/29 00:56:06  vsnyder
! Improve dump of referenceMIF
!
! Revision 2.10  2015/09/23 22:42:02  vsnyder
! Use type-bound procedures to access geocentric latitude from the quantity
! template.  This corrects a mistake in the imputed latitude of the
! generated viewing plane, which was stored as geocentric latitude
! notwithstanding the value of the template's LatitudeCoordinate.
!
! Revision 2.9  2015/09/22 23:22:20  vsnyder
! Add a reference MIF index; add option to require stacked and coherent
!
! Revision 2.8  2015/08/25 18:40:45  vsnyder
! Check MIF and MAF range.  Use 0.5 instead of 1.0 as threshold for unfilled
! scVelECR position or value.  The position is calculated from lat and lon
! (there is no height), and its norm sometimes comes out a tiny bit less than
! 1.0.
!
! Revision 2.7  2015/07/29 00:28:15  vsnyder
! Compute Phi
!
! Revision 2.6  2015/07/27 22:30:12  vsnyder
! Use nonzero crossAngles instead of associated crossAngles to set 'across'
!
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
