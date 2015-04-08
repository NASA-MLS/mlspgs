! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Magnetic_Field_Quantity

  implicit NONE
  private
  public :: Get_Magnetic_Field_Quantity
  public :: Get_Magnetic_Field_Qty_Geoc, Get_Magnetic_Field_Qty_XYZ

  interface Get_Magnetic_Field_Quantity
    module procedure Get_Magnetic_Field_Qty_Geoc
    module procedure Get_Magnetic_Field_Qty_XYZ
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Get_Magnetic_Field_Qty_Geoc ( Qty, GeocAltitudeQuantity, &
    & ScVelQuantity, MAF )

  ! Compute the magnetic field at the geolocations given by Qty (lat,lon,surfs).
  ! Determine the viewing plane defined by GeocAltitudeQuantity and
  ! ScVelQuantity. If that is not the same as the orbit plane, compute the
  ! magnetic field at across-track positions given by
  ! Qty%Template%CrossAngles. Put the values in Qty%Value4.

    use ForwardModelGeometry, only: Compute_Viewing_Plane
    use Geometry, only: ERad, SecPerYear, To_Cart, To_XYZ
    use IGRF_Int, only: Feldc, FeldCof
    use Intrinsic, only: L_GeocAltitude, L_GeodAltitude, L_Gauss
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use MLSStringLists, only: SwitchDetail
    use QuantityTemplates, only: Epoch, RT
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: Dump, VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    type (vectorValue_t), intent(in) :: GeocAltitudeQuantity
    type (vectorValue_t), intent(in) :: ScVelQuantity
    integer, intent(in), optional :: MAF ! to use for GPH quantity

    real :: B(3)                      ! Magnetic field
    integer :: Cross                  ! Index for cross-track angle
    integer, save :: Details = -10    ! From switchDetails('gmag')
    integer :: Instance               ! Loop counter
    integer :: InstanceOr1            ! Index
    integer :: MAF1, MAFn             ! MAFs to do
    integer :: Me = -1                ! String index for tracing
    character(len=8), save :: Options = ''
    integer :: Surf                   ! Loop counter
    integer :: SurfOr1                ! Index
    real    :: XYZ(3)                 ! ECR equivalent of (lat,lon,height),
                                      ! first an unit vector, then in units of
                                      ! ERad.

    call trace_begin ( me, 'Get_Magnetic_Field_Qty_Geoc', &
      & cond=levels(emit) > 1 )

    if ( details < -1 ) then ! only once
      details = switchDetail(switches,'gmag')
      if ( switchDetail(switches,'clean') > -1 ) options = '-c'
    end if

    if ( qty%template%verticalCoordinate /= l_geocAltitude .and. &
       & qty%template%verticalCoordinate /= l_geodAltitude ) then
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Magnetic field vertical coordinate is not geodetic or geocentric altitude' )
      return ! Hopefully, we don't get here
    end if

    MAF1 = 1; MAFn = qty%template%noInstances
    if ( present(MAF) ) then
      MAF1 = MAF; MAFn = MAF
    end if
    call compute_viewing_plane ( Qty, GeocAltitudeQuantity, ScVelQuantity, &
      MAF1, MAFn )

    ! Assume the time is close enough to constant that one call to
    ! FELDCOF is accurate enough.

    call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

    do instance = MAF1, MAFn
      instanceOr1 = merge ( 1, instance, qty%template%coherent )
      do surf = 1, qty%template%noSurfs
        surfOr1 = merge ( 1, surf, qty%template%stacked )
        do cross = 1, qty%template%noCrossTrack
          if ( qty%template%verticalCoordinate == l_geocAltitude ) then
            ! Convert lat(deg), lon(deg) to an unit vector in ECR.
            xyz = to_xyz ( qty%template%geodLat3(surfOr1,instance,cross), &
                         & qty%template%lon3(surfOr1,instance,cross) )
            ! Make the length of XYZ the geocentric altitude in units of ERad.
            xyz = xyz * qty%template%surfs(surf,instanceOr1) / ( 1.0e3 * ERad )
          else ! Assume the vertical coordinate is geodetic altitude in meters.
            ! To_Cart wants altitude in km above sea level.
            ! Qty%Template%Surfs is in meters.
            ! It produces a vector with a length in units of the Earth radius.
            call to_cart ( real( [ qty%template%geodLat3(surfOr1,instance,cross), &
                                 & qty%template%lon3(surfOr1,instance,cross), &
                                 & qty%template%surfs(surf,instanceOr1) / 1.0e3 ] ), &
                                 & xyz )
          end if
          ! Compute the field at and w.r.t. cartesian coordinates in XYZ.
          call feldc ( xyz, b )
          ! The first dimension of value3 is field components, not channels.
          qty%value4 ( 1:3, surf, instance, cross ) = b
        end do
      end do
    end do
    if ( details > -1 ) call dump ( qty, details=details, options=options )
    qty%template%unit = l_gauss

    call trace_end ( 'Get_Magnetic_Field_Qty_Geoc', cond=levels(emit) > 1 )

  end subroutine Get_Magnetic_Field_Qty_Geoc

  subroutine Get_Magnetic_Field_Qty_XYZ ( Qty, XYZ )

  ! Compute the magnetic field at the geolocations given by XYZ, and
  ! put the values in Qty%Value3.

    use Geometry, only: SecPerYear
    use IGRF_Int, only: Feldc, FeldCof
    use MLSStringLists, only: SwitchDetail
    use QuantityTemplates, only: Epoch
    use Toggles, only: Emit, Levels, Switches, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: Dump, VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    real, intent(in) :: XYZ(3,size(qty%value3,2),size(qty%value3,3))

    real :: B(3)                      ! Magnetic field
    integer, save :: Details = -10    ! From switchDetails('gmag')
    integer :: INSTANCE               ! Loop counter
    integer :: Me = -1                ! String index for tracing
    character(len=8), save :: Options = ''
    integer :: SURF                   ! Loop counter
    integer :: SURFOR1                ! Index

    call trace_begin ( me, 'Get_Magnetic_Field_Qty_XYZ', &
      & cond=levels(emit) > 1 )

    if ( details < -1 ) then ! only once
      details = switchDetail(switches,'gmag')
      if ( switchDetail(switches,'clean') > -1 ) options = '-c'
    end if

    ! Assume the time is close enough to constant that one call to
    ! FELDCOF is accurate enough.

    call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

    do instance = 1, qty%template%noInstances
      do surf = 1, qty%template%noSurfs
        if ( qty%template%stacked ) then
          surfOr1 = 1
        else
          surfOr1 = surf
        end if
        ! Compute the field at and w.r.t. cartesian coordinates in XYZ.
        call feldc ( xyz(1:3,surfOr1,instance), b )
        ! The first dimension of value3 is field components, not channels.
        qty%value3 ( 1:3, surf, instance) = b
      end do
    end do

    if ( details > -1 ) call dump ( qty, details=details, options=options )

    call trace_end ( 'Get_Magnetic_Field_Qty_XYZ', cond=levels(emit) > 1 )

  end subroutine Get_Magnetic_Field_Qty_XYZ

  ! ----------------------------------------------  not_used_here  -----
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Magnetic_Field_Quantity

! $Log$
! Revision 2.5  2015/04/08 01:55:13  vsnyder
! Check for geocentric or geodetic altitude vertical coordinate
!
! Revision 2.4  2015/03/28 02:13:06  vsnyder
! Add cross-track magnetic field grids
!
! Revision 2.3  2014/01/11 01:28:53  vsnyder
! Decruftification
!
! Revision 2.2  2013/08/16 02:31:25  vsnyder
! Add MAF argument, get To_Cart from Geometry
!
! Revision 2.1  2013/08/13 02:22:41  vsnyder
! Initial commit
!
