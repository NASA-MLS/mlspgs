! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ManipulateVectorQuantities ! Various routines for manipulating vectors

  ! This modules contains routines needed for manipulating vectors.

  use MLSMessageModule, only: MLSMessage,MLSMSG_Error,MLSMSG_Allocate,MLSMSG_Deallocate
  use MLSCommon, only: r8
  use MLSNumerics, only: Hunt
  use VectorsModule, only: VectorValue_T
  use Dump_0, only: Dump
  use Output_m, only: Output

  implicit none

  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256), private :: Id = &
    & "$Id$"
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

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
    soughtQuantity, instance )
    ! This returns the instance index into a stacked quantity for the
    ! instance 'closest' to the given instance in an unstacked one
    type (VectorValue_T), intent(in) :: referenceQuantity ! e.g. temperature
    type (VectorValue_T), intent(in) :: soughtQuantity ! e.g. ptan, radiance
    integer, intent(in) :: instance

    ! Local variables
    integer :: FIRSTGUESS               ! The result of hunt
    integer :: LOWGUESS                 ! A profile below firstGuess
    integer :: HIGHGUESS                ! A profile above firstGuess
    integer :: BESTGUESS                ! The result
    real (r8) :: COST                   ! A cost function
    real (r8) :: BESTCOST               ! The best cost function

    ! Executable code
    ! We'll skip the error checking we could do at this point, for speed.

    ! First we'll do a hunt to get ourselves in the right area.  Might as
    ! well start looking where we think that will be.
    call Hunt ( referenceQuantity%template%phi(1,:), &
      & soughtQuantity%template%phi(1,instance), firstGuess, &
      & start=max(min(instance,referenceQuantity%template%noInstances),1), &
      & allowTopValue=.true., nearest=.true. )

    ! Now check the ones either side
    lowGuess = max ( firstGuess-1, 1 )
    highGuess = min ( firstGuess+1, referenceQuantity%template%noInstances )
    bestCost = 0.0
    do firstGuess = lowGuess, highGuess
      cost = sum ( abs ( &
        & referenceQuantity%template%phi(1,firstGuess) - &
        & soughtQuantity%template%phi(:,instance) ) )
      if ( ( firstGuess == lowGuess ) .or. ( cost < bestCost ) ) then
        bestGuess = firstGuess
        bestCost = cost
      end if
    end do
    FindOneClosestInstance = bestGuess
  end function FindOneClosestInstance

  ! --------------------------------------- DoHGridsMatch --------------
  logical function DoHGridsMatch ( a, b )
    ! Returns true if quantities have same hGrid information
    type (VectorValue_T), intent(in) :: A ! First quantity
    type (VectorValue_T), intent(in) :: B ! Second quantity

    ! Local parameters
    real (r8), parameter :: PHITOL = 0.01 ! Tolerance in angle

    ! Executable code
    DoHGridsMatch = .false.
    if ( a%template%noInstances /= b%template%noInstances ) return

    if ( any(abs(a%template%phi - &
      &          b%template%phi) > PhiTol) ) return

    DoHGridsMatch = .true.
  end function DoHGridsMatch

  ! --------------------------------------- DoVGridsMatch --------------
  logical function DoVGridsMatch ( a, b )
    ! Returns true if quantities have same hGrid information
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

end module ManipulateVectorQuantities
  
! $Log$
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
