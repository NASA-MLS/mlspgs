! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE MLSL2Common              ! Common data types for the MLSL2 program
!=============================================================================

  USE MLSCommon
  USE MLSMessageModule

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------

  ! This module simply contains data types that are common to the MLSL2
  ! program. Many such data types are `owned' by specific modules. For example
  ! the L2CF datatype is clearly under the remit of an L2CF module.  However,
  ! There are some datatypes that are passed between modules fairly freely, and
  ! not really `owned' by any particular component.  These are detailed in this
  ! module.

  !---------------------------------------------------------------------------
  
  ! Type DataProcessingRange_T has been superseded by TAI93_Range_T in
  ! MLSCommon, so this module is now just a placeholder for future type
  ! definitions.

  ! --------------------------------------------------------------------------

!=============================================================================
END MODULE MLSL2Common
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:57:03  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.7  2000/05/16 19:57:25  livesey
! Removed GetMIFTimesFromMAFTimes, belongs in ConstructQuantityTemplates
!

