! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1T      ! MLS Level 1 software for the THz module
!============================================================================

  USE MLSL1Common, ONLY : L1ProgType, THzType
  USE OpenInit, ONLY : OpenAndInitialize
  USE SortQualifyTHz, ONLY : SortAndQualifyTHz
  USE THzCalibration, ONLY : CalibrateTHz
  USE THzRadiances, ONLY : ProcessLimbData
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit
  USE Close_Files, ONLY: CloseFiles

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2
  LOGICAL :: more_data = .TRUE.

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 THz processing.")

  L1ProgType = THzType

  CALL OpenAndInitialize

  CALL SortAndQualifyTHz

  DO

     CALL CalibrateTHz (more_data)

     CALL ProcessLimbData

     IF (.NOT. more_data) EXIT   !all done

  ENDDO

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 THz data processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)


!=============================================================================
END PROGRAM MLSL1T
!=============================================================================

! $Log$
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
