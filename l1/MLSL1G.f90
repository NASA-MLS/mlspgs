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
PROGRAM MLSL1G       ! MLS Level 1 software for the GHz module
!============================================================================

  USE OpenInit, ONLY: OpenAndInitialize
  USE SortQualify, ONLY: SortAndQualify
  USE Calibration, ONLY: Calibrate
  USE Radiances, ONLY: CalcLimbRads
  USE L1BOutUtils, ONLY: OutputL1Bdata
  USE SpectralBaseline, ONLY: UpdateBaselines
  USE DACsUtils, ONLY: FinalizeDACSdata
  USE GHzBaseline, ONLY: LatBinRads, OutputBaselinedRads
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

        CALL LatBinRads

     ENDIF

     IF (.NOT. more_data) EXIT

  ENDDO

  CALL FinalizeDACSdata

  CALL UpdateBaselines

  CALL OutputBaselinedRads

  CALL CloseFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 GHz data processing successfully completed!")

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1G
!=============================================================================

! $Log$
! Revision 2.5  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.4  2004/12/01 17:10:30  perun
! Add call to FinalizeDACSdata
!
! Revision 2.3  2004/11/10 15:35:40  perun
! Add call to UpdateBaselines
!
! Revision 2.2  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
