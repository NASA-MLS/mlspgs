! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ManipulateVectorQuantities ! Various routines for manipulating vectors

  ! This modules contains routines needed for manipulating vectors.

  use MLSMessageModule, only: MLSMessage, MLSMSG_Error
  use MLSCommon, only: r8, rv
  use MLSNumerics, only: Hunt
  use VectorsModule, only: VectorValue_T
  use Intrinsic, only: L_PHITAN, L_CHANNEL, L_NONE, PHYQ_ANGLE, PHYQ_PROFILES

  implicit none

  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256), private :: Id = &
    & "$Id$"
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

  private

  public :: AnyGoodDataInQty, FindClosestInstances, FindOneClosestInstance, &
    & FindInstanceWindow, DoHGridsMatch, DoVGridsMatch, DoFGridsMatch, &
    & DoQtysDescribeSameThing

contains

  ! --------------------------------------- AnyGoodDataInQty --------------
  logical function AnyGoodDataInQty ( a, a_precision )
    ! Returns true if any of the mask != char(0)
    ! Returns true if any of the precision values >= 0 (if present)
    type ( VectorValue_T ), intent(in) :: a           ! Precision of qty
    type ( VectorValue_T ), intent(in), optional :: a_precision ! Precision of qty

    ! Executable code
    AnyGoodDataInQty = .false.
    if ( present(a_precision) ) then
      if ( .not. associated ( a_precision%values ) ) return
      AnyGoodDataInQty = any(a_precision%values >= 0._rv)
      return
    endif
    AnyGoodDataInQty = .true.
    if ( .not. associated ( a%mask ) ) return
    AnyGoodDataInQty = any(a%mask == char(0) )
  end function AnyGoodDataInQty

  ! ------------------------------ FindClosestInstances -----------------
  subroutine FindClosestInstances ( referenceQuantity, soughtQuantity,&
    referenceIndices )
    ! This subroutine is similar to FindOneClosestInstance (and calls it in
    ! fact), and finds an array of closest instances

    ! Dummy arguments
    type (VectorValue_T), intent(in) :: referenceQuantity ! e.g. temperature
    type (VectorValue_T), intent(in) :: soughtQuantity ! e.g. ptan, radiance
    integer, dimension(soughtQuantity%template%noInstances), &
      & intent(out) :: referenceIndices

    ! Local variables
    integer :: instance

    ! Executable code
    do instance = 1, soughtQuantity%template%noInstances
      referenceIndices(instance) = FindOneClosestInstance ( &
        & referenceQuantity, soughtQuantity, instance )
    end do
  end subroutine FindClosestInstances

  ! ---------------------------------------- FindOneClosestInstance -----
  integer function FindOneClosestInstance ( referenceQuantity, &
    soughtQuantity, instance, useValue )
    use HGridsDatabase, only: FINDCLOSESTMATCH
    ! This returns the instance index into a stacked quantity for the
    ! instance 'closest' to the given instance in an unstacked one
    type (VectorValue_T), intent(in) :: referenceQuantity ! e.g. temperature
    type (VectorValue_T), intent(in) :: soughtQuantity ! e.g. ptan, radiance
    integer, intent(in) :: instance
    logical, intent(in), optional :: USEVALUE ! For phiTan as sought quantity

    ! Local variables
    logical :: MYUSEVALUE
    real (r8), dimension(:,:), pointer :: SEEK ! The thing to look for

    ! Executable code
    ! We'll skip the error checking we could do at this point, for speed.

    ! First we'll do a hunt to get ourselves in the right area.  Might as
    ! well start looking where we think that will be.
    myUseValue = .false.
    if ( present(useValue) ) myUseValue = useValue

    if ( myUseValue ) then
      if ( soughtQuantity%template%quantityType == l_phiTan ) then
        seek => soughtQuantity%values
      else
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'Cannot use useValue option for non phiTan quantities' )
      end if
    else
      seek => soughtQuantity%template%phi
    end if
    ! Call FindClosestMatch to do the work
    FindOneClosestInstance = FindClosestMatch ( referenceQuantity%template%phi(1,:), &
      & seek, instance )
  end function FindOneClosestInstance

  ! --------------------------------------- FindInstanceWindow ---------
  subroutine FindInstanceWindow ( quantity, phiTan, maf, phiWindow, &
    & windowUnits, windowStart, windowFinish )
    ! This returns the start end end of a window into a quantity such as
    ! temperature for a given instance of a minor frame quantity
    type (VectorValue_T), intent(in) :: QUANTITY ! Quantity e.g. temperature
    type (VectorValue_T), intent(in) :: PHITAN ! Phitan information
    integer, intent(in) :: MAF          ! Major frame sought
    real (r8), intent(in) :: PHIWINDOW  ! Window size input
    integer, intent(in) :: WINDOWUNITS
    integer, intent(out) :: WINDOWSTART ! Output window start
    integer, intent(out) :: WINDOWFINISH ! Output window finish

    ! Local variables
    integer :: CLOSESTINSTANCE
    real(r8) :: PHIMIN, PHIMAX           ! Limiting values of phi for this MAF

    ! Executable code
    if ( phiWindow == 0.0 ) then
      ! Just return closest instances
      closestInstance = FindOneClosestInstance ( quantity, phiTan, maf, &
        & useValue=.true. )
      windowStart = closestInstance
      windowFinish = closestInstance
    else if ( windowUnits == PHYQ_Profiles ) then
      ! Return n profiles either side of the closest instance
      closestInstance = FindOneClosestInstance ( quantity, phiTan, maf, &
        & useValue=.true. )
      windowStart = max ( 1, closestInstance - nint ( (phiWindow-1)/2 ) )
      windowFinish = min ( quantity%template%noInstances, &
        & closestInstance + nint ( (phiWindow-1)/2 ) )
    else if ( windowUnits == PHYQ_Angle ) then
      phiMin = minval ( phiTan%values(:,maf) ) - phiWindow/2.0
      phiMax = maxval ( phiTan%values(:,maf) ) + phiWindow/2.0
      call Hunt ( quantity%template%phi(1,:), phiMin, windowStart, &
        & allowTopValue=.true. )
      call Hunt ( quantity%template%phi(1,:), phiMax, windowFinish, &
        & allowTopValue=.true. )
      windowStart = max ( 1, windowStart - 1 )
      windowFinish = min ( quantity%template%noInstances, windowFinish + 1 )
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Invalid units for window specification' )
    end if
  end subroutine FindInstanceWindow

  ! --------------------------------------- DoHGridsMatch --------------
  logical function DoHGridsMatch ( a, b, spacingOnly )
    ! Returns true if quantities have same hGrid information
    type (VectorValue_T), intent(in) :: A ! First quantity
    type (VectorValue_T), intent(in) :: B ! Second quantity
    logical, optional, intent(in) :: SPACINGONLY

    ! Local parameters
    logical :: MYSPACINGONLY
    real (r8), parameter :: PHITOL = 0.01 ! Tolerance in angle
    real (r8) :: MINA, MINB, MAXA, MAXB ! Information on a and b

    ! Executable code
    DoHGridsMatch = .false.

    mySpacingOnly = .false.
    if ( present ( spacingOnly ) ) mySpacingOnly = spacingOnly

    if ( a%template%stacked .neqv. b%template%stacked ) return
    if ( mySpacingOnly .and. .not. a%template%stacked ) return

    if ( .not. mySpacingOnly ) then
      if ( a%template%noInstances /= b%template%noInstances ) return
      if ( any(abs(a%template%phi - &
        &          b%template%phi) > PhiTol) ) return
      DoHGridsMatch = .true.
    else
      ! Here we default to true
      doHGridsMatch = .true.
      if ( a%template%noInstances == 1 ) return
      if ( b%template%noInstances == 1 ) return
      mina = minval ( &
        & a%template%phi(1,2:a%template%noInstances) - &
        & a%template%phi(1,1:a%template%noInstances-1) )
      minb = minval ( &
        & b%template%phi(1,2:b%template%noInstances) - &
        & b%template%phi(1,1:b%template%noInstances-1) )
      maxa = maxval ( &
        & a%template%phi(1,2:a%template%noInstances) - &
        & a%template%phi(1,1:a%template%noInstances-1) )
      maxb = minval ( &
        & b%template%phi(1,2:b%template%noInstances) - &
        & b%template%phi(1,1:b%template%noInstances-1) )
      doHGridsMatch = all ( (/ &
        & maxa-mina, maxb-minb, abs(maxa-maxb), abs(mina-minb) /) < phiTol )
    end if

  end function DoHGridsMatch

  ! --------------------------------------- DoVGridsMatch --------------
  logical function DoVGridsMatch ( a, b )
    ! Returns true if quantities have same vGrid information
    type (VectorValue_T), intent(in) :: A ! First quantity
    type (VectorValue_T), intent(in) :: B ! Second quantity

    ! Local parameters
    real (r8), parameter :: ZTOL = 0.01 ! Tolerance in whatever coordinate

    ! Executable code
    DoVGridsMatch = .false.
    if ( a%template%noSurfs /= b%template%noSurfs ) return
    if ( a%template%verticalCoordinate /= &
      &  b%template%verticalCoordinate ) return
    if ( a%template%coherent .neqv. b%template%coherent ) return
    if ( a%template%regular .neqv. b%template%regular ) return
    if ( ( .not. a%template%coherent) .and. &
      &  ( a%template%noInstances /= b%template%noInstances ) ) return
    if ( any(abs(a%template%surfs - &
      &          b%template%surfs) > zTol) ) return
    if (.not. a%template%regular ) then
      if ( any(a%template%surfIndex /= b%template%surfIndex) .or. &
        &  any(a%template%chanIndex /= b%template%chanIndex) ) return
    end if

    DoVGridsMatch = .true.
  end function DoVGridsMatch

  ! --------------------------------------- DoVGridsMatch --------------
  logical function DoFGridsMatch ( a, b )
    ! Returns true if the quantities have the same fGrid information
    type ( VectorValue_T ), intent(in) :: A ! First quantity
    type ( VectorValue_T ), intent(in) :: B ! Second quantity

    ! Local paramterts
    real (r8), parameter :: FTOL = 1.0e-3 ! 1 kHz

    ! Executable code
    DoFGridsMatch = .false.
    if ( a%template%frequencyCoordinate /= b%template%frequencyCoordinate ) return
    if ( a%template%noChans /= b%template%noChans ) return
    select case ( a%template%frequencyCoordinate )
    case ( l_none )
    case ( l_channel )
      if ( a%template%signal /= b%template%signal ) return
      if ( a%template%sideband /= b%template%sideband ) return
    case default
      if ( .not. associated ( a%template%frequencies ) .or. &
        & .not. associated ( b%template%frequencies ) ) return
      if ( any ( shape(a%template%frequencies) /= &
        & shape(b%template%frequencies) ) ) return
      if ( any ( abs ( a%template%frequencies - &
        & b%template%frequencies ) > fTol ) ) return
    end select
    DoFGridsMatch = .true.
  end function DoFGridsMatch

  ! ---------------------------------- DoQtysDescribeSameThing ----
  logical function DoQtysDescribeSameThing ( a, b )
    ! Returns true if the quantities describe the same geophysical
    ! parameter, albeit at a different resolution perhaps
    type ( VectorValue_T ), intent(in) :: A ! First quantity
    type ( VectorValue_T ), intent(in) :: B ! Second quantity

    DoQtysDescribeSameThing = .false.
    if ( a%template%quantityType /= b%template%quantityType ) return
    if ( a%template%logBasis .neqv. b%template%logBasis ) return
    if ( a%template%verticalCoordinate /= b%template%verticalCoordinate ) return
    if ( a%template%unit /= b%template%unit ) return
    if ( a%template%signal /= b%template%signal ) return
    if ( a%template%sideband /= b%template%sideband ) return
    if ( a%template%instrumentModule /= b%template%instrumentModule ) return
    if ( a%template%radiometer /= b%template%radiometer ) return
    if ( a%template%molecule /= b%template%molecule ) return
    DoQtysDescribeSameThing = .true.

  end function DoQtysDescribeSameThing

  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module ManipulateVectorQuantities
  
! $Log$
! Revision 2.23  2003/07/07 20:21:34  livesey
! Now uses the FindClosestMatch function
!
! Revision 2.22  2003/01/26 04:42:20  livesey
! Added handling of profiles/angle units for phiWindow
!
! Revision 2.21  2002/11/22 01:07:13  vsnyder
! Delete USE'd but unreferenced symbols
!
! Revision 2.20  2002/10/08 00:09:11  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.19  2002/09/19 00:30:36  pwagner
! Added AnyGoodDataInQty; set DoQtysDescribeSameThing at its start
!
! Revision 2.18  2002/09/10 20:47:44  livesey
! Added DoQtysDescribeSameThing
!
! Revision 2.17  2002/08/23 01:24:18  livesey
! Added DoFGridsMatch
!
! Revision 2.16  2002/07/25 08:43:19  mjf
! Initialised DoHGridsMatch to .false. at start of function.
!
! Revision 2.15  2002/07/17 06:01:27  livesey
! Fixed bugs in DoH/VGrids match
!
! Revision 2.14  2002/06/12 16:53:32  livesey
! Tidied up some public/private stuff
!
! Revision 2.13  2002/06/12 16:50:39  livesey
! Added findInstanceWindow
!
! Revision 2.12  2002/02/06 01:32:58  livesey
! Rewrote FindOneClosestInstance and FindClosestInstances to reflect the
! way in which they are mostly called, and to fix a bug.
!
! Revision 2.11  2001/11/08 01:05:06  livesey
! Fixed a minor sort of bug in FindOneClosestQuantity
!
! Revision 2.10  2001/09/14 18:02:52  livesey
! Bug fix in FindOneClosestInstance and FindClosestInstances.
! Will probably come back to these and rewrite them some time.
!
! Revision 2.9  2001/09/11 01:27:27  livesey
! Bug fixes
!
! Revision 2.8  2001/09/09 21:17:30  livesey
! Imported FindOneClosestInstance from branch
!
! Revision 2.7.2.1  2001/09/08 23:46:40  livesey
! Added FindOneClosestInstance
!
! Revision 2.7  2001/05/11 00:03:41  livesey
! Fixed but with DoVGridsMatch
!
! Revision 2.6  2001/05/10 23:29:27  livesey
! Added DoHGridsMatch and DoVGridsMatch
!
! Revision 2.5  2001/03/08 02:21:08  livesey
! Fixed bug, wasn't setting minloc!
!
! Revision 2.4  2001/03/02 01:31:36  livesey
! Regular commit
!
! Revision 2.3  2001/02/27 17:18:20  livesey
! Moved ValidateVectorQuantity into vectors module
!
