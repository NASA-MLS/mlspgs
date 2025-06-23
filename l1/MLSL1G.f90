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

  USE MLSL1RunConfig, ONLY: MLSL1Executable ! name of current executable
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
  USE MLSL1Debug, ONLY: closeMLSL1DebugFiles
  USE MLSL1Config, ONLY: L1Config
  use MLSStrings, only: WriteIntsToChars

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
  CHARACTER (len=255) :: ComVecsFile
  character (len=16) :: nWallsChars



  MLSL1Executable='mlsl1g'
  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "============= Start EOS MLS Level 1 GHz processing ================ ")
  PRINT *,"============= Start EOS MLS Level 1 GHz processing ================ "

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
  ! Print informative message if we encountered unusually many walls
  if ( L1Config%Output%NumAttenuationWalls >= &
    & L1Config%Output%MaxAttenuationWalls ) then
    call WriteIntsToChars ( L1Config%Output%NumAttenuationWalls, nWallsChars )
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "Unusual number of walls on this date: " // trim(nWallsChars) )
    PRINT *,"Unusual number of walls on this date: " // trim(nWallsChars)
  endif

  CALL FinalizeDACSdata

  CALL UpdateBaselines

  CALL OutputBaselinedRads

  CALL CloseFiles

  CALL closeMLSL1DebugFiles

  CALL MLSMessage (MLSMSG_Info, ModuleName, &
       & "EOS MLS Level 1 GHz data processing successfully completed!")

  !<whd:debugging>
  close(12)
  !</whd:debugging>

  CALL MLSMessageExit (NORMAL_EXIT_STATUS)

!=============================================================================
END PROGRAM MLSL1G
!=============================================================================

! $Log$
! Revision 2.7  2023/06/06 22:34:58  pwagner
! Note if unusually many attenuation walls
!
! Revision 2.6  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.5.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
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
