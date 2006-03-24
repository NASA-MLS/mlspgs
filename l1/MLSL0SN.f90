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

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 0 Band Switches database processing.")

  CALL OpenInitBsw

  CALL ExamineSciData

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
    & "EOS MLS Level 0 Band Switches database updated successfully")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)


!=============================================================================
END PROGRAM MLSL0SN
!=============================================================================

! $Log$
! Revision 2.1  2006/03/24 15:09:21  perun
! Initial release of main program to update Band Switch Network database
!
