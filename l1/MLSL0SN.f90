! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!==============================================================================
PROGRAM MLSL0SN  ! MLS Level 0 software to update Band Switch Network database
!==============================================================================

  USE MLSL1RunConfig, ONLY: MLSL1Executable ! name of current executable
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit
  USE BandSwitchesUpdate, ONLY: OpenInitBsw, ExamineSciData, CloseFiles

  IMPLICIT NONE

!---------------------------- RCS Ident Info ------------------------------
  CHARACTER (len=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2

  MLSL1Executable='mlsl0sn'
  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 0 Band Switches database processing.")

  CALL OpenInitBsw

  CALL ExamineSciData

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
    & "EOS MLS Level 0 Band Switches database updated successfully")
  PRINT *,"MLSL0SN: EOS MLS Level 0 Band Switches database updated successfully"

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)


!=============================================================================
END PROGRAM MLSL0SN
!=============================================================================

! $Log$
! Revision 2.2  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.1.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.1  2006/03/24 15:09:21  perun
! Initial release of main program to update Band Switch Network database
!
