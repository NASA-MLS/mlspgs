! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!============================================================================
PROGRAM MLSL1       ! MLS Level 1 software
!============================================================================

  USE OpenInit, ONLY : OpenAndInitialize
  USE SortQualify, ONLY : SortAndQualify
  USE Calibration, ONLY : Calibrate
  USE Radiances, ONLY : CalcLimbRads
  USE L1BOutUtils, ONLY : OutputL1Bdata
  USE Close_Files, ONLY: CloseFiles
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMessageExit

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: NORMAL_EXIT_STATUS = 2
  LOGICAL :: more_data, do_calib

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Start EOS MLS Level 1 data processing.")

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

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1
!=============================================================================

! $Log$
! Revision 2.5  2002/11/13 15:26:42  perun
! Restore consistent coding style
!
! Revision 2.4  2002/11/12 21:50:03  jdone
! Remove obsolete variables from previous.
!
! Revision 2.3  2002/11/07 21:35:33  jdone
! Added HDF4/HDF5 switch.
!
! Revision 2.2  2002/04/03 21:42:00  pwagner
! Sets status on normal exit to 2
!
! Revision 2.1  2001/02/23 18:26:11  perun
! Version 0.5 commit
!
