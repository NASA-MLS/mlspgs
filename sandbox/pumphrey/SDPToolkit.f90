! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SDPToolkit     ! Substitute for the essential toolkit routines
!=============================================================================

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
  !---------------------------------------------------------------------------
  ! This module gives substitutes for the essential toolkit routines used by
  ! library code used in both the toolkit and non toolkit environment.
  ! This is HCPs personal version, which has the same module name and  
  ! filename as the true SDP Toolkit module. I think this deserves a 
  ! separate directory from mlspgs/lib, so one can tell the compiler to look
  ! there first for source files

CONTAINS

  FUNCTION PGS_SMF_GenerateStatusReport(message)
    CHARACTER (LEN=*), INTENT(IN) :: message
    INTEGER :: PGS_SMF_GenerateStatusReport

    PRINT*,message
    PGS_SMF_GenerateStatusReport=0
  END FUNCTION PGS_SMF_GenerateStatusReport

!=============================================================================
END MODULE SDPToolkit
!=============================================================================

