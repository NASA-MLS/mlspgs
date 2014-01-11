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
  public :: Get_Magnetic_Field_Qty_GPH, Get_Magnetic_Field_Qty_XYZ

  interface Get_Magnetic_Field_Quantity
    module procedure Get_Magnetic_Field_Qty_GPH
    module procedure Get_Magnetic_Field_Qty_XYZ
  end interface

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Get_Magnetic_Field_Qty_GPH ( Qty, GPHQty, MAF )

  ! Compute the magnetic field at the geolocations given by Qty (lat,lon)
  ! and GPHQty (altitude), and put the values in Qty%Value3.

    use Geometry, only: SecPerYear, To_Cart
    use IGRF_INT, only: FELDC, FELDCOF
    use MLSStringLists, only: SwitchDetail
    use QuantityTemplates, only: Epoch
    use Toggles, only: Switches
    use VectorsModule, only: Dump, VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    type (vectorValue_t), intent(inout) :: GPHQty
    integer, intent(in), optional :: MAF ! to use for GPH quantity

    real :: B(3)                      ! Magnetic field
    integer, save :: Details = -10    ! From switchDetails('mag')
    integer :: GPHInstance
    integer :: INSTANCE               ! Loop counter
    character(len=8), save :: options = ''
    integer :: SURF                   ! Loop counter
    integer :: SURFOR1                ! Index
    real    :: XYZ(3)                 ! lat, lon, height for to_cart

    if ( details < -1 ) then ! only once
      details = switchDetail(switches,'mag')
      if ( switchDetail(switches,'clean') > -1 ) options = '-c'
    end if

    ! Assume the time is close enough to constant that one call to
    ! FELDCOF is accurate enough.

    call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

    do instance = 1, qty%template%noInstances
      gphInstance = instance
      if ( present(MAF) ) gphInstance = MAF
      do surf = 1, qty%template%noSurfs
        if ( qty%template%stacked ) then
          surfOr1 = 1
        else
          surfOr1 = surf
        end if
        ! Convert (/ lat(deg), lon(deg), height(km) /) to cartesian (e-radius)
        call to_cart ( real( (/ qty%template%geodLat(surfOr1,instance), &
          &                     qty%template%lon(surfOr1,instance), &
          &                     gphQty%values(surf,gphInstance)*1.0e-3 /) ), xyz )
        ! Compute the field at and w.r.t. cartesian coordinates in XYZ.
        call feldc ( xyz, b )
        ! The first dimension of value3 is field components, not channels.
        qty%value3 ( 1:3, surf, instance) = b
      end do
    end do

    if ( details > -1 ) call dump ( qty, details=details, options=options )

  end subroutine Get_Magnetic_Field_Qty_GPH

  subroutine Get_Magnetic_Field_Qty_XYZ ( Qty, XYZ )

  ! Compute the magnetic field at the geolocations given by XYZ, and
  ! put the values in Qty%Value3.

    use Geometry, only: SecPerYear
    use IGRF_INT, only: FELDC, FELDCOF
    use MLSStringLists, only: SwitchDetail
    use QuantityTemplates, only: Epoch
    use Toggles, only: Switches
    use VectorsModule, only: Dump, VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    real, intent(in) :: XYZ(3,size(qty%value3,2),size(qty%value3,3))

    real :: B(3)                      ! Magnetic field
    integer, save :: Details = -10    ! From switchDetails('mag')
    integer :: INSTANCE               ! Loop counter
    character(len=8), save :: options = ''
    integer :: SURF                   ! Loop counter
    integer :: SURFOR1                ! Index

    if ( details < -1 ) then ! only once
      details = switchDetail(switches,'mag')
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
! Revision 2.3  2014/01/11 01:28:53  vsnyder
! Decruftification
!
! Revision 2.2  2013/08/16 02:31:25  vsnyder
! Add MAF argument, get To_Cart from Geometry
!
! Revision 2.1  2013/08/13 02:22:41  vsnyder
! Initial commit
!
