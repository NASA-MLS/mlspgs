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
PROGRAM MLSL1BOA       ! MLS Level 1 software for L1BOA data
!============================================================================

  USE OpenInitBOA, ONLY : OpenAndInitBOA
  USE L1BOAutils, ONLY : OutputL1BOA
  USE CloseBOA, ONLY : CloseBOAfile
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit

  IMPLICIT NONE

!---------------------------- RCS Ident Info ------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
!---------------------------------------------------------------------------

  LOGICAL :: more_data
  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 BOA processing.")

  CALL OpenAndInitBOA

  DO

     CALL OutputL1BOA (more_data)

     IF (.NOT. more_data) EXIT

  ENDDO

  CALL CloseBOAfile

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 BOA processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1BOA
!=============================================================================

! $Log$
! Revision 2.2  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.1  2003/10/24 19:38:36  perun
! Version 1.3 commit
!
!
