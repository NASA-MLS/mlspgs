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

  subroutine Get_Magnetic_Field_Qty_Geoc ( Qty, ScVelQuantity, &
    & GeocAltitudeQuantity, GPHQuantity, Regular, ReferenceMIF )

  ! Compute the magnetic field at the geolocations given by Qty (lat,lon,surfs).
  ! Compute the orbit plane determined by the velocity and geolocation in
  ! ScVelQuantity. Determine the viewing plane defined by GeocAltitudeQuantity
  ! and ScVelQuantity. If that is not the same as the orbit plane, compute the
  ! magnetic field at across-track positions given by Qty%Template%CrossAngles.
  ! Put the values in Qty%Value4.

  ! We assume that the quantity types have been checked, that the relationship
  ! between the magnetic field's vertical coordinate and the existence of the
  ! other quantities has been checked, and that the relationship between the
  ! existence, number, and values of cross-angles and the existence of the other
  ! quantities has been checked.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate
    use ForwardModelGeometry, only: Compute_Viewing_Plane
    use QuantityTemplates, only: RT
    use Toggles, only: Emit, Levels
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    type (vectorValue_t), intent(in), optional :: ScVelQuantity
    type (vectorValue_t), intent(in), optional :: GeocAltitudeQuantity ! tangent
    type (vectorValue_t), intent(in), optional :: GPHQuantity
    logical, intent(in), optional :: Regular ! coherent and stacked, default false
    integer, intent(in), optional :: ReferenceMIF ! Used only if Regular

    integer :: Me = -1                 ! String index for tracing
    integer :: Stat                    ! From allocate or deallocate
    real(rt), allocatable :: XYZs(:,:,:,:) ! Geocentric Cartesian coordinates
                                       ! in viewing plane.
                                       ! 3 X instances X surfs X cross angles

    call trace_begin ( me, 'Get_Magnetic_Field_Qty_Geoc', &
      & cond=levels(emit) > 1 )

    allocate ( xyzs ( 3, &                     ! XYZ
                    & size(qty%value4,2),   &  ! NoSurfs
                    & size(qty%value4,3),   &  ! NoInstances or 1
                    & size(qty%value4,4) ), &  ! NoCrossAngles
                    & stat=stat )
    call test_allocate ( stat, moduleName, 'XYZs' )
    
    call compute_viewing_plane ( Qty, ScVelQuantity, GeocAltitudeQuantity, &
      gphQuantity, XYZs, regular, referenceMIF )

    call get_magnetic_field_qty_XYZ ( Qty, XYZs )

    deallocate ( xyzs, stat=stat )
    call test_deallocate ( stat, moduleName, 'XYZs' )

    call trace_end ( 'Get_Magnetic_Field_Qty_Geoc', cond=levels(emit) > 1 )

  end subroutine Get_Magnetic_Field_Qty_Geoc

  subroutine Get_Magnetic_Field_Qty_XYZ ( Qty, XYZs )

  ! Compute the magnetic field at the geolocations given by XYZs, and
  ! put the values in Qty%Value3.

    use Geometry, only: SecPerYear
    use IGRF_Int, only: Feldcm, FeldCof
    use Intrinsic, only: PHYQ_Gauss
    use MLSStringLists, only: SwitchDetail
    use QuantityTemplates, only: Epoch, RT
    use Toggles, only: Emit, Levels, Switches
    use Trace_m, only: Trace_Begin, Trace_End
    use VectorsModule, only: Dump, VectorValue_T

    type (vectorValue_t), intent(inout) :: Qty
    ! Cartesian ECR coordinates, meters:
    real(rt), intent(in) :: XYZs(3,size(qty%value4,2),size(qty%value4,3),size(qty%value4,4))

    integer :: Across                  ! Loop counter for cross-track angles
    real :: B(3)                       ! Magnetic field
    integer, save :: Details = -10     ! From switchDetails('gMag')
    integer :: Inst                    ! Loop counter for instances
    integer :: Me = -1                 ! String index for tracing
    character(len=8), save :: Options = ''
    integer :: Surf                    ! Loop counter for surfaces

    call trace_begin ( me, 'Get_Magnetic_Field_Qty_XYZ', &
      & cond=levels(emit) > 1 )

    ! Assume the time is close enough to constant that one call to
    ! FELDCOF is accurate enough.

    call feldcof ( real(qty%template%time(1,1)/secPerYear + epoch) )

    do across = 1, size(xyzs,4)
      do inst = 1, size(xyzs,3)
        do surf = 1, size(xyzs,2)
          ! Compute the field at and w.r.t. cartesian coordinates in xyzs.
          call feldcm ( real(xyzs(:3,surf,inst,across)), b )
          ! The first dimension of value4 is field components, not channels.
          qty%value4 ( 1:3, surf, inst, across ) = b
        end do
      end do
    end do
    qty%template%unit = phyq_gauss

    if ( details < -1 ) then ! only once
      details = switchDetail(switches,'gMag')
      if ( switchDetail(switches,'clean') > -1 ) options = '-c'
    end if
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
! Revision 2.12  2016/01/23 02:51:54  vsnyder
! Use procedures that want coordinates in meters instead of ERad
!
! Revision 2.11  2015/10/27 23:03:30  vsnyder
! Change dump switch from gmag to gMag because mag causes the forward model
! to dump the magnetic field along the integration path.
!
! Revision 2.10  2015/09/22 23:20:38  vsnyder
! Add a reference MIF number
!
! Revision 2.9  2015/07/29 00:28:36  vsnyder
! Repair a comment
!
! Revision 2.8  2015/07/23 23:44:47  vsnyder
! qty%template%unit should be phyq_gauss, not l_gauss
!
! Revision 2.7  2015/05/28 23:15:18  vsnyder
! Fiddle with tracing
!
! Revision 2.6  2015/04/29 01:23:08  vsnyder
! Coordinates now given as XYZ, not taken from the magnetic field quantity
!
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
