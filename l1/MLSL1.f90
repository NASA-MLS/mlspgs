! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1       ! MLS Level 1 software
!============================================================================

  USE OpenInit, ONLY : OpenAndInitialize
  USE SortQualify, ONLY : SortAndQualify
  USE Calibration, ONLY : Calibrate
  USE Radiances, ONLY : CalcLimbRads
  USE L1BOutUtils, ONLY : OutputL1Bdata
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit
  USE Close_Files, ONLY: CloseFiles

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  LOGICAL :: more_data, do_calib
  integer, parameter :: NORMAL_EXIT_STATUS = 2

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 processing.")

  CALL OpenAndInitialize

  DO

     CALL SortAndQualify (more_data, do_calib)

     IF (do_calib) THEN

        CALL Calibrate

        CALL CalcLimbRads

        CALL OutputL1Bdata

     ENDIF

     IF (.NOT. more_data) EXIT

  ENDDO

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 data processing successfully completed!")

  call MLSMessageExit(NORMAL_EXIT_STATUS)
!=============================================================================
END PROGRAM MLSL1
!=============================================================================

! $Log$
! Revision 2.2  2002/04/03 21:42:00  pwagner
! Sets status on normal exit to 2
!
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
