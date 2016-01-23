! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Quantity_Geometry

  ! This module contains functions to compute Cartesian coordinates from
  ! geolocations in quantity templates.

  ! The computations depend upon whether the coordinates are geodetic
  ! or geocentric.

  implicit NONE
  private

  public :: XYZ, XYZ_ECR_1, XYZ_ECR_All

  interface XYZ
    module procedure XYZ_ECR_1, XYZ_ECR_All
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  function XYZ_ECR_1 ( Qty, Surf, Inst, Height, Cross ) result ( XYZ )
    ! Get Cartesian ECR coordinates, in meters, for one surface and instance.

    use Constants, only: Rad2Deg
    use Geometry, only: GeocToGeodLat, GeodToECRm, GeodToGeocLat, To_XYZ
    use Intrinsic, only: L_Geocentric, L_GeocAltitude, L_Geodetic, &
      & L_GeodAltitude, L_Zeta
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_t, RT

    type(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst
    real(rt), intent(in), optional :: Height ! Geodetic meters, used iff zeta
    integer, intent(in), optional :: Cross
    real(rt) :: XYZ(3)

    integer :: C
    integer :: InstOr1
    real(rt) :: Lat
    integer :: SurfOr1

    c = 1
    if ( present(cross) ) c = cross

    instOr1 = merge ( 1, inst, qty%coherent )
    surfOr1 = merge ( 1, surf, qty%stacked )

    select case ( qty%verticalCoordinate )
    case ( l_geocAltitude )
      if ( qty%latitudeCoordinate == l_geodetic ) then
        lat = geodToGeocLat(qty%geodLat3(surfOr1,inst,c))*rad2deg
      else
        lat = qty%geodLat3(surfOr1,inst,c)
      end if
      ! Get a unit vector in the right direction
      xyz = to_xyz ( lat, qty%lon3(surfOr1,inst,c) )
      ! Extend it to the correct altitude
      xyz = xyz * qty%surfs(surf,instOr1)
    case ( l_geodAltitude )
      if ( qty%latitudeCoordinate == l_geocentric ) then
        lat = geocToGeodLat(qty%geodLat3(surfOr1,inst,c))
      else
        lat = qty%geodLat3(surfOr1,inst,c)
      end if
      xyz = geodToECRm ( [ lat, qty%lon3(surfOr1,inst,c), &
                       &   qty%surfs(surf,instOr1) ] )
    case ( l_zeta )
      if ( .not. present(height) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Vertical coordinate is zeta but height is not provided' )
      if ( qty%latitudeCoordinate == l_geocentric ) then
        lat = geocToGeodLat(qty%geodLat3(surfOr1,inst,c))
      else
        lat = qty%geodLat3(surfOr1,inst,c)
      end if
      xyz = geodToECRm ( [ lat, qty%lon3(surfOr1,inst,c), height ] )
    case default
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Vertical coordinate other than geocentric, geodetic, or zeta not supported' )
    end select

  end function XYZ_ECR_1

  function XYZ_ECR_All ( Qty, Heights, Cross  ) result ( XYZ)
    ! Get Cartesian ECR coordinates, in meters, for all surfaces and instances.

    use Allocate_Deallocate, only: Allocate_Test
    use Constants, only: Rad2Deg
    use Geometry, only: GeocToGeodLat, GeodToGeocLat, GeodToECRm, To_XYZ
    use Intrinsic, only: L_Geocentric, L_GeocAltitude, L_Geodetic, &
      & L_GeodAltitude, L_Zeta
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_t, RT

    type(quantityTemplate_t), intent(in) :: Qty
    real(rt), allocatable :: XYZ(:,:,:) ! 3 x noSurfs x noInstances
    real(rt), intent(in), optional :: Heights(:,:) ! Geodetic meters, used iff zeta
    integer, intent(in), optional :: Cross

    integer :: C
    integer :: Inst, InstOr1
    real(rt) :: Lat
    integer :: Surf, SurfOr1

    c = 1
    if ( present(cross) ) c = cross

    call allocate_test ( xyz, 3, qty%noSurfs, qty%noInstances, moduleName, 'XYZ' )

    select case ( qty%verticalCoordinate )
    case ( l_geocAltitude )
      do inst = 1, qty%noInstances
        instOr1 = merge ( 1, inst, qty%coherent )
        do surf = 1, qty%noSurfs
          surfOr1 = merge ( 1, surf, qty%stacked )
          if ( qty%latitudeCoordinate == l_geodetic ) then
            lat = geodToGeocLat(qty%geodLat3(surfOr1,inst,c))*rad2deg
          else
            lat = qty%geodLat3(surfOr1,inst,c)
          end if
          ! Get a unit vector in the right direction
          xyz(:,surf,inst) = to_xyz ( lat, qty%lon3(surfOr1,inst,c) )
          ! Extend it to the correct altitude
          xyz(:,surf,inst) = xyz(:,surf,inst) * qty%surfs(surf,instOr1)
        end do
      end do
    case ( l_geodAltitude )
      do inst = 1, qty%noInstances
        instOr1 = merge ( 1, inst, qty%coherent )
        do surf = 1, qty%noSurfs
          surfOr1 = merge ( 1, surf, qty%stacked )
          if ( qty%latitudeCoordinate == l_geocentric ) then
            lat = geocToGeodLat(qty%geodLat3(surfOr1,inst,c))
          else
            lat = qty%geodLat3(surfOr1,inst,c)
          end if
          xyz(:,surf,inst) = geodToECRm ( [ lat, qty%lon3(surfOr1,inst,c), &
                                        &   qty%surfs(surf,instOr1) ] )
        end do
      end do
    case ( l_zeta )
      if ( .not. present(heights) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Vertical coordinate is zeta but heights are not provided' )
      do inst = 1, qty%noInstances
        instOr1 = merge ( 1, inst, qty%coherent )
        do surf = 1, qty%noSurfs
          surfOr1 = merge ( 1, surf, qty%stacked )
          if ( qty%latitudeCoordinate == l_geocentric ) then
            lat = geocToGeodLat(qty%geodLat3(surfOr1,inst,c))
          else
            lat = qty%geodLat3(surfOr1,inst,c)
          end if
          xyz(:,surf,inst) = geodToECRm ( [ lat, qty%lon3(surfOr1,inst,c), &
                                        &   heights(surf,instOr1) ] )
        end do
      end do
    case default
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Vertical coordinate other than geocentric, geodetic, or zeta not supported' )
    end select

  end function XYZ_ECR_All

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Quantity_Geometry

! $Log$
! Revision 2.5  2016/01/23 02:48:47  vsnyder
! Replace To_Cart with GeodToECRm
!
! Revision 2.4  2015/09/22 23:13:19  vsnyder
! Add a cross-track coordinate index
!
! Revision 2.3  2015/04/30 02:55:35  vsnyder
! Allow geocentric altitude/latitude to coexist with geodetic latitude/altitude
!
! Revision 2.2  2015/04/29 01:21:21  vsnyder
! Add ability to compute XYZ if vertical is zeta
!
! Revision 2.1  2015/04/25 02:04:57  vsnyder
! Initial commit
!
