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

  function XYZ_ECR_1 ( Qty, Surf, Inst, Height ) result ( XYZ )
    ! Get Cartesian ECR coordinates, in meters, for one surface and instance.

    use Geometry, only: To_Cart, To_XYZ
    use Intrinsic, only: L_Geocentric, L_GeocAltitude, L_Geodetic, &
      & L_GeodAltitude, L_Zeta
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_t, RT

    type(quantityTemplate_t), intent(in) :: Qty
    integer, intent(in) :: Surf, Inst
    real(rt), intent(in), optional :: Height ! Geodetic meters, used iff zeta
    real(rt) :: XYZ(3)

    integer :: InstOr1
    integer :: SurfOr1

    instOr1 = merge ( 1, inst, qty%coherent )
    surfOr1 = merge ( 1, surf, qty%stacked )

    select case ( qty%verticalCoordinate )
    case ( l_geocAltitude )
      if ( qty%latitudeCoordinate == l_geodetic ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Geocentric altitude with geodetic latitude not supported')
      ! Get a unit vector in the right direction
      xyz = to_xyz ( qty%geodLat(surfOr1,inst), qty%lon(surfOr1,inst) )
      ! Extend it to the correct altitude
      xyz = xyz * qty%surfs(surf,instOr1)
    case ( l_geodAltitude )
      if ( qty%latitudeCoordinate == l_geocentric ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Geodetic altitude with geocentric latitude not supported')
      call to_cart ( [ qty%geodLat(surfOr1,inst), qty%lon(surfOr1,inst), &
                   &   qty%surfs(surf,instOr1)/1000.0 ], xyz, km=.true. )
      xyz = xyz * 1000.0_rt ! Convert km to m
    case ( l_zeta )
      if ( .not. present(height) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Vertical coordinate is zeta but height is not provided' )
      call to_cart ( [ qty%geodLat(surfOr1,inst), qty%lon(surfOr1,inst), &
                   &   height/1000.0 ], xyz, km=.true. )
    case default
      call MLSMessage ( MLSMSG_Error, moduleName, &
        & 'Vertical coordinate other than geocentric, geodetic, or zeta not supported' )
    end select

  end function XYZ_ECR_1

  function XYZ_ECR_All ( Qty, Heights  ) result ( XYZ)
    ! Get Cartesian ECR coordinates, in meters, for all surfaces and instances.

    use Allocate_Deallocate, only: Test_Allocate
    use Geometry, only: To_Cart, To_XYZ
    use Intrinsic, only: L_Geocentric, L_GeocAltitude, L_Geodetic, &
      & L_GeodAltitude, L_Zeta
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    use QuantityTemplates, only: QuantityTemplate_t, RT

    type(quantityTemplate_t), intent(in) :: Qty
    real(rt), allocatable :: XYZ(:,:,:) ! 3 x noSurfs x noInstances
    real(rt), intent(in), optional :: Heights(:,:) ! Geodetic meters, used iff zeta

    integer :: Inst, InstOr1
    integer :: Surf, SurfOr1
    integer :: Stat

    instOr1 = size(qty%surfs,2)    ! merge ( 1, qty%noInstances, qty%coherent )
    surfOr1 = size(qty%geodLat,1)  ! merge ( 1, qty%noSurfs, qty%stacked )

    allocate ( xyz(3,qty%noSurfs, qty%noInstances), stat=stat )
    call test_allocate ( stat, moduleName, 'XYZ' )

    select case ( qty%verticalCoordinate )
    case ( l_geocAltitude )
      if ( qty%latitudeCoordinate == l_geodetic ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Geocentric altitude with geodetic latitude not supported')
      do inst = 1, qty%noInstances
        instOr1 = merge ( 1, inst, qty%coherent )
        do surf = 1, qty%noSurfs
          surfOr1 = merge ( 1, surf, qty%stacked )
          ! Get a unit vector in the right direction
          xyz(:,surf,inst) = to_xyz ( qty%geodLat(surfOr1,inst), qty%lon(surfOr1,inst) )
          ! Extend it to the correct altitude
          xyz(:,surf,inst) = xyz(:,surf,inst) * qty%surfs(surf,instOr1)
        end do
      end do
    case ( l_geodAltitude )
      if ( qty%latitudeCoordinate == l_geocentric ) &
        call MLSMessage ( MLSMSG_Error, moduleName, &
          & 'Geodetic altitude with geocentric latitude not supported')
      do inst = 1, qty%noInstances
        instOr1 = merge ( 1, inst, qty%coherent )
        do surf = 1, qty%noSurfs
          surfOr1 = merge ( 1, surf, qty%stacked )
          call to_cart ( [ qty%geodLat(surfOr1,inst), &
                       &   qty%lon(surfOr1,inst), &
                       &   qty%surfs(surf,instOr1)/1000.0 ], &
                       & xyz(:,surf,inst), km=.true. )
          xyz(:,surf,inst) = xyz(:,surf,inst) * 1000.0_rt ! Convert km to m
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
          call to_cart ( [ qty%geodLat(surfOr1,inst), &
                       &   qty%lon(surfOr1,inst), &
                       &   heights(surf,instOr1)/1000.0 ], &
                       & xyz(:,surf,inst), km=.true. )
          xyz(:,surf,inst) = xyz(:,surf,inst) * 1000.0_rt ! Convert km to m
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
! Revision 2.2  2015/04/29 01:21:21  vsnyder
! Add ability to compute XYZ if vertical is zeta
!
! Revision 2.1  2015/04/25 02:04:57  vsnyder
! Initial commit
!
