! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L1LogUtils
!=============================================================================

  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSL1Common, ONLY: L1BFileInfo, R8
  USE MLSL1Config, ONLY: L1Config
  USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ExamineData, LogStatus

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE BeginEnd_T
     INTEGER :: EngMAFno(2), SciMAFno(2)
     INTEGER :: OrbitNo(2)
     INTEGER :: TotalMAFcount(2)
     REAL(R8) :: EngTAI(2), SciTAI(2)
  END TYPE BeginEnd_T

  TYPE (BeginEnd_T) :: BeginEnd
  INTEGER :: EngMAFs, SciMAFs
  INTEGER :: EngGaps, SciGaps

  CHARACTER(LEN=27) :: asciiUTC(2)
  INTEGER :: unit  ! for writing log file
  INTEGER :: eng_warns, eng_errs, sci_warns, sci_errs
  INTEGER :: MAF_dif, PGS_stat
  REAL(R8) :: last_TAI
  REAL :: MAF_dur
  TYPE (TAI93_Range_T) :: TAI_range
  
CONTAINS

!=============================================================================
  SUBROUTINE ExamineEngData
!=============================================================================

    USE EngUtils, ONLY: NextEngMAF
    USE EngTbls, ONLY: EngMAF

    LOGICAL :: first = .TRUE.
    LOGICAL :: more_data = .TRUE.

    PRINT *, 'Examining eng data...'
    WRITE (unit, *) ''
    WRITE (unit, *) '################ Engineering data scan ###############'
    WRITE (unit, *) ''

    EngMAFs = 0
    EngGaps = 0
    DO

       CALL NextEngMAF (more_data)

       IF (.NOT. more_data) EXIT

       IF (first) THEN

! Catch up to the Sci MAFno:

          DO
             IF (EngMAF%MAFno == BeginEnd%SciMAFno(1)) EXIT
             CALL NextEngMAF (more_data)
          ENDDO

          BeginEnd%EngMAFno(1) = EngMAF%MAFno
          BeginEnd%EngTAI(1) = EngMAF%secTAI
          BeginEnd%TotalMAFcount(1) = EngMAF%TotalMAF
          last_TAI = EngMAF%secTAI
          first = .FALSE.
       ENDIF

! Check for data gaps:

       MAF_dif = INT ((EngMAF%secTAI - last_TAI) / MAF_dur + 0.1)
       IF (MAF_dif > 1) THEN
          EngGaps = EngGaps + 1
          Eng_Warns = Eng_Warns + 1
          WRITE (unit, *) '##### WARNING! Data Gap:'
          WRITE (unit, *) ''
          WRITE (unit, *) 'MAFs missing: ', (MAF_dif-1)
          WRITE (unit, *) 'MAFno gap: ', BeginEnd%EngMAFno(2), EngMAF%MAFno
          PGS_stat = PGS_TD_TAItoUTC (last_TAI, asciiUTC(1))
          PGS_stat = PGS_TD_TAItoUTC (EngMAF%secTAI, asciiUTC(2))
          WRITE (unit, *) 'UTC gap: ', asciiUTC(1)//' to '//asciiUTC(2)
       ENDIF
       last_TAI = EngMAF%secTAI

       EngMAFs = EngMAFs + 1
       BeginEnd%EngMAFno(2) = EngMAF%MAFno
       BeginEnd%EngTAI(2) = EngMAF%secTAI
       BeginEnd%TotalMAFcount(2) = EngMAF%TotalMAF

       more_data = ABS (EngMAF%secTAI - TAI_range%endTime) > (1.5*MAF_dur)
       IF (.NOT. more_data) more_data = (EngMAF%MAFno /= BeginEnd%SciMAFno(2))
       IF (.NOT. more_data) EXIT
   ENDDO

   IF (Eng_Warns == 0 .AND. Eng_Errs == 0) &
        WRITE (unit, *) '##### No Warnings and no Errors ###'

  END SUBROUTINE ExamineEngData

!=============================================================================
  SUBROUTINE ExamineSciData
!=============================================================================

    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: NextSciMAF

    LOGICAL :: first = .TRUE.
    LOGICAL :: more_data = .TRUE.

    PRINT *, 'Examining sci data...'
    WRITE (unit, *) ''
    WRITE (unit, *) '################## Science data scan #################'
    WRITE (unit, *) ''

    SciMAFs = 0
    SciGaps = 0
    DO

       DO
          CALL NextSciMAF (more_data)
          IF (.NOT. more_data .OR. SciMAF(0)%secTAI >= TAI_range%startTime) EXIT
       ENDDO

       IF (more_data) more_data = SciMAF(0)%secTAI <= TAI_range%endTime
       IF (.NOT. more_data) EXIT

       IF (first) THEN
          BeginEnd%SciMAFno(1) = SciMAF(0)%MAFno
          BeginEnd%SciTAI(1) = SciMAF(0)%secTAI
          last_TAI = SciMAF(0)%secTAI
          first = .FALSE.
       ENDIF

! Check for data gaps:

       MAF_dif = INT ((SciMAF(0)%secTAI - last_TAI) / MAF_dur + 0.1)
       IF (MAF_dif > 1) THEN
          SciGaps = SciGaps + 1
          Sci_Warns = Sci_Warns + 1
          WRITE (unit, *) '##### WARNING! Data Gap:'
          WRITE (unit, *) ''
          WRITE (unit, *) 'MAFs missing: ', (MAF_dif-1)
          WRITE (unit, *) 'MAFno gap: ', BeginEnd%SciMAFno(2), SciMAF(0)%MAFno
          PGS_stat = PGS_TD_TAItoUTC (last_TAI, asciiUTC(1))
          PGS_stat = PGS_TD_TAItoUTC (SciMAF(0)%secTAI, asciiUTC(2))
          WRITE (unit, *) 'UTC gap: ', asciiUTC(1)//' to '//asciiUTC(2)
       ENDIF
       last_TAI = SciMAF(0)%secTAI

       SciMAFs = SciMAFs + 1
       BeginEnd%SciMAFno(2) = SciMAF(0)%MAFno
       BeginEnd%SciTAI(2) = SciMAF(0)%secTAI

    ENDDO

    IF (Sci_Warns == 0 .AND. Sci_Errs == 0) &
         WRITE (unit, *) '##### No Warnings and no Errors ###'

  END SUBROUTINE ExamineSciData

!=============================================================================
  SUBROUTINE OutputLogSummary
!=============================================================================
    
    WRITE (unit, *) ''
    WRITE (unit, *) '################ Science data summary ################'
    WRITE (unit, *) ''
    WRITE (unit, *) 'SciMAFs: ', SciMAFs
    WRITE (unit, *) 'SciMAFnos: ', BeginEnd%SciMAFno
    WRITE (unit, *) 'SciTAI: ', BeginEnd%SciTAI(1), ' to ', BeginEnd%SciTAI(2)
    PGS_stat = PGS_TD_TAItoUTC (BeginEnd%SciTAI(1), asciiUTC(1))
    PGS_stat = PGS_TD_TAItoUTC (BeginEnd%SciTAI(2), asciiUTC(2))
    WRITE (unit, *) 'SciUTC: ', asciiUTC(1)//' to '//asciiUTC(2)
    WRITE (unit, *) 'SciGaps: ', SciGaps
    WRITE (unit, *) 'Sci_Warns, Sci_Errs: ', Sci_Warns, Sci_Errs

    WRITE (unit, *) ''
    WRITE (unit, *) '############## Engineering data summary ##############'
    WRITE (unit, *) ''
    WRITE (unit, *) 'EngMAFs: ', EngMAFs
    WRITE (unit, *) 'EngMAFnos: ', BeginEnd%EngMAFno
    WRITE (unit, *) 'EngTAI: ', BeginEnd%EngTAI(1), ' to ', BeginEnd%EngTAI(2)
    PGS_stat = PGS_TD_TAItoUTC (BeginEnd%EngTAI(1), asciiUTC(1))
    PGS_stat = PGS_TD_TAItoUTC (BeginEnd%EngTAI(2), asciiUTC(2))
    WRITE (unit, *) 'EngUTC: ', asciiUTC(1)//' to '//asciiUTC(2)
    WRITE (unit, *) 'EngTotalMAF: ', BeginEnd%TotalMAFcount
    WRITE (unit, *) 'EngGaps: ', EngGaps
    WRITE (unit, *) 'Eng_Warns, Eng_Errs: ', Eng_Warns, Eng_Errs

  END SUBROUTINE OutputLogSummary

!=============================================================================
  SUBROUTINE ExamineData
!=============================================================================

    USE InitPCFs, ONLY: L1PCF

    CHARACTER(LEN=80) :: message
    INTEGER :: window_MAFs

    unit = L1BFileInfo%LogId

    WRITE (unit, *) ''
    WRITE (unit, *) '################## Begin MLSL1log ####################'
    WRITE (unit, *) ''
    WRITE (unit, *) 'PCF filename: '//TRIM (L1PCF%PCF_filename)
    WRITE (unit, *) 'L1CF filename: '//TRIM (L1PCF%L1CF_filename)
    WRITE (unit, *) 'Input Start/End UTC: '// &
         TRIM (L1PCF%StartUTC)//' to '//TRIM (L1PCF%EndUTC)

! Determine TAI/MAF range to read:

    window_MAFs = L1Config%Calib%CalWindow
    MAF_dur = L1Config%Calib%MIF_duration * L1Config%Calib%MIFsPerMAF
    TAI_range = L1Config%Input_TAI
    TAI_range%startTime = TAI_range%startTime - (window_MAFs/2 * MAF_dur)
    TAI_range%endTime = TAI_range%endTime + ((window_MAFs/2 - 1)* MAF_dur)

! Init number of warnings and errors

    eng_warns = 0
    eng_errs = 0
    sci_warns = 0
    sci_errs = 0

! Must examine science data first!

    CALL ExamineSciData

    CALL ExamineEngData

    CALL OutputLogSummary

    IF (eng_warns > 0 .OR. sci_warns > 0) THEN
       WRITE (unit, *) ''
       WRITE (message, *) '############ Warnings:', (sci_warns+eng_warns)
       WRITE (unit, '(A)') message
    ENDIF

    IF (eng_errs > 0 .OR. sci_errs > 0) THEN
       WRITE (unit, *) ''
       WRITE (message, *) '############ Errors:', (sci_errs+eng_errs)
       WRITE (unit, '(A)') message
    ENDIF

    WRITE (unit, *) ''
    WRITE (unit, *) '################### End MLSL1log #####################'

  END SUBROUTINE ExamineData

!=============================================================================
  SUBROUTINE LogStatus

    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Warning
    USE Machine, ONLY: Shell_Command

    IF (eng_warns > 0 .OR. sci_warns > 0) THEN
       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
            "EOS MLS Level 1 log file " //TRIM(L1BFileInfo%LogFilename)//&
            " reports warnings!")
    ENDIF

    IF (eng_errs > 0 .OR. sci_errs > 0) THEN

! send email on failure (next version!):

       !CALL Shell_Command ("/usr/lib/sendmail -t <mail.vsp")

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "EOS MLS Level 1 log file " //TRIM(L1BFileInfo%LogFilename)//&
            " reports errors!")
    ENDIF

  END SUBROUTINE LogStatus

!=============================================================================
END MODULE L1LogUtils
!=============================================================================

! $Log$
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
