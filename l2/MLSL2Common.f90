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

CONTAINS

  ! -------------------------------------------- FindFirst --------------
  integer function FindFirst ( condition )
    ! Find the first logical in the array that is true
    logical, dimension(:), intent(in) :: CONDITION

    ! Local variables
    integer :: I                        ! Loop counter

    ! Executable code
    FindFirst = -1
    do i = 1, size(condition)
      if ( condition(i) ) then
        FindFirst = i
        return
      end if
    end do
  end function FindFirst

!=============================================================================
!=============================================================================
END MODULE MLSL2Common
!=============================================================================

!
! $Log$
! Revision 2.1  2001/06/22 05:09:04  livesey
! Moved FindFirst in from L2Parallel
!
! Revision 2.0  2000/09/05 18:57:03  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.7  2000/05/16 19:57:25  livesey
! Removed GetMIFTimesFromMAFTimes, belongs in ConstructQuantityTemplates
!

