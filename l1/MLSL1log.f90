! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

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
! Revision 2.3  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
