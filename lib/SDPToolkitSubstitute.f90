! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
!MODULE SDPToolkitSubstitute     ! Substitute for the essential toolkit routines
!=============================================================================

!  IMPLICIT NONE

!  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
!  CHARACTER (LEN=256) :: Id = &
!       "$Id$"
!  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------


  ! This module gives substitutes for the essential toolkit routines used by
  ! library code used in both the toolkit and non toolkit environment.

!CONTAINS

  FUNCTION PGS_SMF_GenerateStatusReport(message)
    CHARACTER (LEN=*), INTENT(IN) :: message
    INTEGER :: PGS_SMF_GenerateStatusReport

    PRINT*,message
    PGS_SMF_GenerateStatusReport=0
  END FUNCTION PGS_SMF_GenerateStatusReport

!=============================================================================
!END MODULE SDPToolkitSubstitute
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 17:41:07  dcuddy
! Change revision to 2.0
!
! Revision 1.3  1999/12/03 00:14:18  livesey
! Nightly checkin
!
