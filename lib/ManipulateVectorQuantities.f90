! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module ManipulateVectorQuantities ! Various routines for manipulating vectors

  ! This modules contains routines needed for manipulating vectors.

  use MLSMessageModule, only: MLSMessage,MLSMSG_Error,MLSMSG_Allocate,MLSMSG_Deallocate
  use MLSCommon, only: r8
  use MLSNumerics, only: Hunt
  use VectorsModule, only: VectorValue_T

  implicit none

  !---------------------------- RCS Ident Info -------------------------------
  character (LEN=256), private :: Id = &
    & "$Id$"
  character (LEN=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------
  
  contains


    ! ------------------------------ FindClosestInstances -----------------
    subroutine FindClosestInstances(referenceQuantity,soughtQuantity,&
         referenceIndices)
      ! This subroutine finds the instance (i.e. profile) in a stacked
      ! quantity that is closest to an instance of another quantity which is
      ! typically unstacked or incoherent (e.g. a minor frame quantity such as
      ! ptan).
      
      ! Close is defined in terms of the geodetic orbit angle.  The algorithm
      ! implemented makes several (obvious) assumptions.  Please consult the code
      ! below.
      
      ! Eventually (v1.0) this routine should become unnecessary, as a full 2D
      ! interpolation will be performed.

      ! Dummy arguments
      type (VectorValue_T), intent(IN) :: referenceQuantity ! e.g. temperature
      type (VectorValue_T), intent(IN) :: soughtQuantity ! e.g. ptan, radiance
      integer, dimension(soughtQuantity%template%noInstances), &
           intent(OUT) :: referenceIndices ! Result
      
      ! Local variables
      real (r8), dimension(-1:1) :: costs
      integer :: soughtInstance ! Loop counter
      integer :: instanceOffset ! Loop counter
      integer :: referenceIndex ! Index into reference quantity

      integer, dimension(1) :: minlocResult

      ! Executable code

      ! First check the obvious
      if (.not. referenceQuantity%template%stacked) &
           call MLSMessage(MLSMSG_Error,ModuleName,&
           'Reference quantity must be stacked')
      
      ! First we're going to look for the instance within the reference
      ! quantity that starts below the current one.

      call Hunt(referenceQuantity%template%phi(1,:), &
        & soughtQuantity%template%phi(1,:), referenceIndices, &
        & allowTopValue=.true.)

      ! Now we refine these by looking at the instances found above and the
      ! ones above and below, and choosing the one that is closest over the
      ! entire vertical range.

      do soughtInstance=1,soughtQuantity%template%noInstances
         do instanceOffset= -1,1 ! Look below, at and above
            ! Look into reference quantity, make sure we don't fall off end
            referenceIndex=referenceIndices(soughtInstance)+instanceOffset
            referenceIndex=min(max(referenceIndex,1),&
                 referenceQuantity%template%noInstances)

            ! Assess cost for these
            costs(instanceOffset)=sum(abs(&
                 referenceQuantity%template%phi(1,referenceIndex)-&
                 soughtQuantity%template%phi(:,soughtInstance)))
         end do
         minLocResult=minloc(costs)
         ! Choose best
         referenceIndices(soughtInstance)=referenceIndices(soughtInstance)+&
              minlocResult(1)-2 ! Correct for the -1:1 indexing
      end do

      ! Now, again don't fall off the ends.
      referenceIndices=min(max(referenceIndices,1),&
           referenceQuantity%template%noInstances)

    end subroutine FindClosestInstances

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
