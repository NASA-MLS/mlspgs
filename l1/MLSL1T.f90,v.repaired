! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!============================================================================
PROGRAM MLSL1T      ! MLS Level 1 software for the THz module
!============================================================================

  USE MLSL1RunConfig, ONLY : MLSL1Executable ! name of current executable
  USE MLSL1Common, ONLY : L1ProgType, THzType
  USE OpenInit, ONLY : OpenAndInitialize
  USE SortQualifyTHz, ONLY : SortAndQualifyTHz
  USE THzCalibration, ONLY : CalibrateTHz
  USE THzRadiances, ONLY : ProcessLimbData
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit
  USE Close_Files, ONLY: CloseFiles

  IMPLICIT NONE

!---------------------------- RCS Ident Info ------------------------------
  CHARACTER (len=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2
  character(len=256):: msg

  MLSL1Executable='mlsl1t'
  msg="Start EOS MLS Level 1 THz processing."
  print *,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  PRINT *,TRIM(msg)
  print *,"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  CALL MLSMessage (MLSMSG_Info, ModuleName, TRIM(msg))

  L1ProgType = THzType

  CALL OpenAndInitialize

  CALL SortAndQualifyTHz

  CALL CalibrateTHz

  CALL ProcessLimbData

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 THz data processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)


!=============================================================================
END PROGRAM MLSL1T
!=============================================================================

! $Log$
! Revision 2.4  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.3.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.3  2006/03/24 15:13:35  perun
! Remove DO processing for CalibrateTHz and ProcessLimbData
!
! Revision 2.2  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
