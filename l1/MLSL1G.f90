! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1G       ! MLS Level 1 software for the GHz module
!============================================================================

  USE OpenInit, ONLY : OpenAndInitialize
  USE SortQualify, ONLY : SortAndQualify
  USE Calibration, ONLY : Calibrate
  USE Radiances, ONLY : CalcLimbRads
  USE L1BOutUtils, ONLY : OutputL1Bdata
  USE Close_Files, ONLY : CloseFiles
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  LOGICAL :: more_data, do_calib
  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 GHz processing.")

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
       & "EOS MLS Level 1 GHz data processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1G
!=============================================================================

! $Log$
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
