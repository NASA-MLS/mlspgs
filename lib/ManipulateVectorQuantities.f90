! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE ManipulateVectorQuantities ! Interpolate in coordinate spaces
!=============================================================================

  ! This module is used to interpolate one quantity to the coordinate system of
  ! another, optionally producing derivative matrices.  This is done in
  ! phases, typically by frequency (if applicable), then vertically and finally
  ! horizontally. Each stage can output a derivative matrix_1 for the result
  ! with respect to the input values.

  ! I may add an interpolation to ptan here too at a later stage as that will
  ! be useful.

  USE MLSMessageModule, ONLY: MLSMessage,MLSMSG_Error,MLSMSG_Allocate,MLSMSG_Deallocate
  USE MLSCommon, ONLY: r8
  USE MLSNumerics, ONLY: Hunt
  USE VectorsModule, ONLY: VectorValue_T

  IMPLICIT NONE

  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256), PRIVATE :: Id = &
    & "$Id$"
  CHARACTER (LEN=*), PARAMETER, PRIVATE :: ModuleName= &
    & "$RCSfile$"
  !---------------------------------------------------------------------------
  
  CONTAINS

    ! This subroutine finds the instance (i.e. profile) in a stacked
    ! quantity that is closest to an instance of another quantity which is
    ! typically unstacked or incoherent (e.g. a minor frame quantity such as
    ! ptan).

    ! Close is defined in terms of the geodetic orbit angle.  The algorithm
    ! implemented makes several (obvious) assumptions.  Please consult the code
    ! below.

    ! Eventually (v1.0) this routine should become unnecessary, as a full 2D
    ! interpolation will be performed.

    SUBROUTINE FindClosestInstances(referenceQuantity,soughtQuantity,&
         referenceIndices)

      ! Dummy arguments
      TYPE (VectorValue_T), INTENT(IN) :: referenceQuantity ! e.g. temperature
      TYPE (VectorValue_T), INTENT(IN) :: soughtQuantity ! e.g. ptan, radiance
      INTEGER, DIMENSION(soughtQuantity%template%noInstances), &
           INTENT(OUT) :: referenceIndices ! Result
      
      ! Local variables
      REAL (r8), DIMENSION(-1:1) :: costs
      INTEGER :: soughtInstance ! Loop counter
      INTEGER :: instanceOffset ! Loop counter
      INTEGER :: referenceIndex ! Index into reference quantity

      INTEGER, DIMENSION(1) :: minlocResult

      ! Executable code

      ! First check the obvious
      IF (.NOT. referenceQuantity%template%stacked) &
           CALL MLSMessage(MLSMSG_Error,ModuleName,&
           'Reference quantity must be stacked')
      
      ! First we're going to look for the instance within the reference
      ! quantity that starts below the current one.

      CALL Hunt(referenceQuantity%template%phi(1,:), &
        & soughtQuantity%template%phi(1,:), referenceIndices, &
        & allowTopValue=.TRUE.)

      ! Now we refine these by looking at the instances found above and the
      ! ones above and below, and choosing the one that is closest over the
      ! entire vertical range.

      DO soughtInstance=1,soughtQuantity%template%noInstances
         DO instanceOffset= -1,1 ! Look below, at and above
            ! Look into reference quantity, make sure we don't fall off end
            referenceIndex=referenceIndices(soughtInstance)+instanceOffset
            referenceIndex=MIN(MAX(referenceIndex,1),&
                 referenceQuantity%template%noInstances)

            ! Assess cost for these
            costs(instanceOffset)=SUM(ABS(&
                 referenceQuantity%template%phi(1,referenceIndex)-&
                 soughtQuantity%template%phi(:,soughtInstance)))
         END DO
         minLocResult=MINLOC(costs)
         ! Choose best
         referenceIndices(soughtInstance)=referenceIndices(soughtInstance)+&
              minlocResult(1)-2 ! Correct for the -1:1 indexing
      END DO

      ! Now, again don't fall off the ends.
      referenceIndices=MIN(MAX(referenceIndices,1),&
           referenceQuantity%template%noInstances)

    END SUBROUTINE FindClosestInstances
    
  END MODULE ManipulateVectorQuantities
  
! $Log$
! Revision 2.5  2001/03/08 02:21:08  livesey
! Fixed bug, wasn't setting minloc!
!
! Revision 2.4  2001/03/02 01:31:36  livesey
! Regular commit
!
! Revision 2.3  2001/02/27 17:18:20  livesey
! Moved ValidateVectorQuantity into vectors module
!
