! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Calib  ! Calibrate L0 to L1 Radiance
!=============================================================================

  USE MLSCommon
  USE MLSL1Common
  USE MLSL1Rad
  USE MLSMessageModule

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE Calibrate (l1bInfo, l0Sci, L1Brad, success)
!=============================================================================

    ! Arguments
    
    TYPE (L1BInfo_T), INTENT (IN):: l1bInfo  ! File handles etc. for L1B dataset
    TYPE (L0Sci_T), INTENT (IN) :: l0Sci(*)  ! Lvl0 science data for 1 MAF
    TYPE (Radiance_T), DIMENSION (:), POINTER :: L1Brad
    LOGICAL, INTENT (OUT) :: success   ! successful write

    ! Local

    INTEGER :: i

    DO i = 1, SIZE (L1Brad)

       L1Brad(i)%value = 0.0
       L1Brad(i)%precision = -1.0

    END DO

    success  = .TRUE.  ! In the future, may use integer status instead

    RETURN

  END Subroutine Calibrate

!=============================================================================
END MODULE Calib
!=============================================================================

!
! $Log$
! Revision 1.1  2000/02/08 21:25:56  perun
! Initial version
!
!
