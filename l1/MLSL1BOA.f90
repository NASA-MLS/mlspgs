! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1BOA       ! MLS Level 1 software for L1BOA data
!============================================================================

  USE OpenInitBOA, ONLY : OpenAndInitBOA
  USE L1BOAutils, ONLY : OutputL1BOA
  USE CloseBOA, ONLY : CloseBOAfile
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

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
! Revision 2.1  2003/10/24 19:38:36  perun
! Version 1.3 commit
!
!
