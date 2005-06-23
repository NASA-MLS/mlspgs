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

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

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
! Revision 2.2  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
