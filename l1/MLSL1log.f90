! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1log       ! MLS Level 1 software to produce data log
!============================================================================

  USE MLSL1Common, ONLY: L1ProgType, LogType
  USE OpenInitLog, ONLY: OpenAndInitializeLog
  USE L1LogUtils, ONLY: ExamineData, LogStatus
  USE CalibWeightsFlags, ONLY: DetermineWeightsFlags
  USE Close_Files, ONLY: CloseFiles
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 log processing.")

  L1ProgType = LogType

  CALL OpenAndInitializeLog

  CALL ExamineData

  CALL DetermineWeightsFlags

  CALL CloseFiles

  CALL LogStatus

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 log data processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1log
!=============================================================================

! $Log$
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
