! Copyright 2006, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE L1LogUtils
!=============================================================================

  USE MLSCommon, ONLY: TAI93_Range_T
  USE MLSL1Common, ONLY: L1BFileInfo, R8, Chan_R_T
  !USE MLSL1Config, ONLY: L1Config
  USE MLSL1Config, ONLY: L1Config, GetL1Config
  USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ExamineData, LogStatus, EngMAFs, SciMAFs, MAF_dur, MIF_dur

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

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
  INTEGER :: unit ! for writing log file
  INTEGER :: eng_warns, eng_errs, sci_warns, sci_errs
  INTEGER :: MAF_dif, PGS_stat, last_MIF
  REAL(R8) :: last_TAI, SciEngDiff
  REAL :: MAF_dur, MIF_dur
  TYPE (TAI93_Range_T) :: TAI_range

  TYPE (Chan_R_T) :: Atten_cnts = &  ! FB, MB, WF, DACS
       Chan_R_T (HUGE(1.0), HUGE(1.0), HUGE(1.0), HUGE(1.0))
  
CONTAINS

!=============================================================================
  SUBROUTINE ExamineEngData
!=============================================================================

    ! <whd>Read the L0 PDS science data in the processing range and write out a
    ! temporary file (PCF id=921, nominally named engMAF.dat, a fortran binary
    ! file in sequential mode) which will be read later by
    ! CalibWeightsFlags::ProcessMAFdata.</whd>

    use EngUtils, ONLY: NextEngMAF
    use EngTbls, ONLY: EngMAF, EngPkt
    !use MLSL1Common, ONLY: MaxdataGaps, MaxErroneousCounterMAFs
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Warning, MLSMSG_Info
    USE HighOutput, ONLY: outputNamedValue
    
    logical :: first = .TRUE.
    logical :: more_data = .TRUE.
    integer :: counterMAF = 0
    integer, save :: dataGaps = 0
    integer, save :: erroneousCounterMAFs = 0
    character(len=512) :: msg

    WRITE (unit, *) 'MaxDataGaps = ', L1Config%Calib%MaxDataGaps
    WRITE (unit, *) 'MaxErroneousCounterMAFs = ', L1Config%Calib%MaxErroneousCounterMAFs
    WRITE (unit, *) 'DiffBeginEndEng = ', L1Config%Calib%DiffBeginEndEng
    
    
    print *, 'Examining eng data...'
    write (unit, *) ''
    write (unit, *) '################ Engineering data scan ###############'
    WRITE (unit, *) ''

    EngMAFs = 0
    EngGaps = 0

    DO
       CALL NextEngMAF (more_data)

       IF (.NOT. more_data) EXIT

       IF (first) THEN

! Catch up to the Sci MAFno:

          ! loop until we find the MAF no we want
          DO
             IF (EngMAF%MAFno == BeginEnd%SciMAFno(1)) EXIT
             CALL NextEngMAF (more_data)
             IF (.NOT. more_data) EXIT
          ENDDO

          BeginEnd%EngMAFno(1) = EngMAF%MAFno
          BeginEnd%EngTAI(1) = EngMAF%secTAI
          BeginEnd%TotalMAFcount(1) = EngMAF%TotalMAF
          counterMAF = EngMAF%TotalMAF - 1   ! previous count for comparisons
          last_TAI = EngMAF%secTAI
          first = .FALSE.
       ENDIF

       ! Save data for further processing:
       ! PCF id=921, nominally named engMAF_tmp.dat
       WRITE (L1BFileInfo%EngMAF_unit) EngMAF, EngPkt
 
! Check for data gaps:

       MAF_dif = INT ((EngMAF%secTAI - last_TAI) / MAF_dur + 0.1)
       IF (MAF_dif > 1) THEN
          EngGaps = EngGaps + 1
          Eng_Warns = Eng_Warns + 1
          dataGaps  = dataGaps + 1
          WRITE (unit, *) '##### WARNING! Data Gap:'
          WRITE (unit, *) 'MAFs missing: ', (MAF_dif-1)
          WRITE (unit, *) 'MAFno gap: ', BeginEnd%EngMAFno(2), EngMAF%MAFno
          PGS_stat = PGS_TD_TAItoUTC (last_TAI, asciiUTC(1))
          PGS_stat = PGS_TD_TAItoUTC (EngMAF%secTAI, asciiUTC(2))
          WRITE (unit, *) 'UTC gap: ', asciiUTC(1)//' to '//asciiUTC(2)
          WRITE (unit, *) ''
          !if ( dataGaps > MaxdataGaps ) &
	  if ( dataGaps > L1Config%Calib%MaxDataGaps ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
                     & 'Too many data gaps--must quit now' )
       ENDIF
       last_TAI = EngMAF%secTAI

! Check counter MAF

       IF (counterMAF /= (EngMAF%TotalMAF - 1)) THEN
          Eng_Warns = Eng_Warns + 1
          erroneousCounterMAFs  = erroneousCounterMAFs + 1
          WRITE (unit, *) '##### WARNING! Counter MAF incorrect:'
          WRITE (unit, *) 'Counter MAFs: ', counterMAF, EngMAF%TotalMAF
          PGS_stat = PGS_TD_TAItoUTC (EngMAF%secTAI, asciiUTC(1))
          WRITE (unit, *) 'UTC: ', asciiUTC(1)
          WRITE (unit, *) ''
          if ( erroneousCounterMAFs > L1Config%Calib%MaxErroneousCounterMAFs ) &
            & call MLSMessage ( MLSMSG_Error, ModuleName, &
                     & 'Too many bad counterMAFs--must quit now' )
       ENDIF
       counterMAF = EngMAF%TotalMAF

       EngMAFs = EngMAFs + 1
       BeginEnd%EngMAFno(2) = EngMAF%MAFno
       BeginEnd%EngTAI(2) = EngMAF%secTAI
       BeginEnd%TotalMAFcount(2) = EngMAF%TotalMAF


       ! <vp comment> Check for utterly bogus times
       ! signaled by mismatch between Engineering and Science times
       ! If they are more than 1 week apart, exit with error status
       ! </vp comment>

       ! <whd comment> By the time we get to this point, *all* of the
       ! science MAFs have been read (ExamineSciData is called before
       ! this routine and it loops through all the science data), so
       ! BeginEnd%SciTAI(1) is the start time of the *first* MAF and
       ! BeginEnd%SciTAI(2) is the start time of the *last* MAF of the
       ! run. However, in this section we're looping through the
       ! Engineering MAFs, so BeginEnd%EngTAI(1) is the begin time of
       ! the first Engineering MAF and, at each interation,
       ! BeginEnd%EngTAI(2) the start time for the MAF that's just
       ! been read. Therefore for the first MAF,
       ! BeginEnd%EngTAI(1)=BeginEnd%EngTAI(2) and the `end'
       ! comparisons are between the last Science MAF and the *first*
       ! Eng MAF. This doesn't seem to me to be a very useful
       ! comparison, but it's how the code is written, and I'm not
       ! going to change it now. (If someone does change this, please
       ! modify this comment accordingly)
       !
       ! Therefore,L1Config%Calib%DiffBeginEndEng (which controls how
       ! large that difference is allowed to be and which is in the
       ! L1CF file read during the startup of the Level1 run) should
       ! be large enough to allow for such a comparison during any
       ! normal run. Normally, the runs are for 24 hours and L1 adds
       ! some sloop, so I suggest setting this parameter in the L1CF
       ! code to 25 hours (90000 seconds) Anyone who runs for more
       ! than 24 hours will have to fashion their own L1CF file.</whd
       ! comment>
       SciEngDiff=ABS(BeginEnd%SciTAI(1) - BeginEnd%EngTAI(1))
       IF ( SciEngDiff > L1Config%Calib%DiffBeginEndEng ) THEN 
         CALL MLSMessage ( MLSMSG_Info, ModuleName, &
         & "start time difference between science and engineering data too big")
         ! The messages sent by OutputNamedVariable go to STDOUT
         CALL OutputNamedValue(ModuleName,&
              & "start time difference between science and engineering data too big")

         msg=ModuleName//"L1CF PARAMETER DiffBeginEndEng"
         CALL OutputNamedValue(msg,L1Config%Calib%DiffBeginEndEng)

         msg=ModuleName//"BeginEnd%SciTAI(1)"
         CALL OutputNamedValue(msg,BeginEnd%SciTAI(1))

         msg=ModuleName//"BeginEnd%EngTAI(1)"
         CALL OutputNamedValue(msg,BeginEnd%EngTAI(1))

          msg=moduleName//"Start Time Difference"
          CALL OutputNamedValue(msg,SciEngDiff)
          
          CALL MLSMessage ( MLSMSG_Error, ModuleName, "Aborting!")
       ENDIF
       SciEngDiff=ABS(BeginEnd%SciTAI(2) - BeginEnd%EngTAI(2))
       IF ( SciEngDiff > L1Config%Calib%DiffBeginEndEng ) THEN 
          CALL MLSMessage ( MLSMSG_Info, ModuleName, &
            "end time difference between science and engineering data too big")
         ! The messages sent by OutputNamedVariable go to STDOUT
          CALL OutputNamedValue(ModuleName,&
               & "end time difference between science and engineering data too big")

          msg=ModuleName//"L1CF PARAMETER DiffBeginEndEng"
          CALL OutputNamedValue(msg,L1Config%Calib%DiffBeginEndEng)


          msg=ModuleName//"BeginEnd%SciTAI(2)"
          CALL OutputNamedValue(msg,BeginEnd%SciTAI(2))

          msg=ModuleName//"BeginEnd%EngTAI(2)"
          CALL OutputNamedValue(msg,BeginEnd%EngTAI(2))

          msg=moduleName//"End Time Difference"
          CALL OutputNamedValue(msg,SciEngDiff)


          CALL MLSMessage ( MLSMSG_Error, ModuleName, "Aborting!")
        ENDIF

       more_data = (EngMAF%secTAI <= BeginEnd%SciTAI(2))
       IF (.NOT. more_data) more_data = (EngMAF%MAFno /= BeginEnd%SciMAFno(2))
       IF (.NOT. more_data) EXIT
   ENDDO

   ENDFILE L1BFileInfo%EngMAF_unit

   IF (Eng_Warns == 0 .AND. Eng_Errs == 0) &
        WRITE (unit, *) '##### No Warnings and no Errors ###'

  END SUBROUTINE ExamineEngData

!=============================================================================
  SUBROUTINE ExamineSciData
!=============================================================================


    ! Read the L0 PDS science data in the processing range and write out a
    ! temporary file (PCF id=920, nominally named sciMAF.dat, a fortran binary
    ! in sequential mode) which will be read later by
    ! CalibWeightsFlags::ProcessMAFdata.

    USE L0_sci_tbls, ONLY: SciMAF
    USE SciUtils, ONLY: NextSciMAF
    USE MLSL1Common, ONLY: FBnum, MBnum, WFnum, DACSnum, MaxMIFs

    INTEGER :: i, mindx,lastGoodMAF
    LOGICAL :: first = .TRUE.
    LOGICAL :: more_data = .TRUE.
    LOGICAL :: doneAttens = .FALSE.
    LOGICAL :: MIFmask(MaxMIFs)
    INTEGER, PARAMETER :: MIFindx(MaxMIFs) = (/ (i, i=1, MaxMIFs) /)

    PRINT *, trim('Examining sci data...')
    WRITE (unit, *) ''
    WRITE (unit, *) '################## Science data scan #################'
    WRITE (unit, *) ''

    SciMAFs = 0
    SciGaps = 0
    DO

       DO
          CALL NextSciMAF (more_data)
          IF (.NOT. more_data) EXIT
          ! <whd> what is this first test doing? </whd>
          IF (ALL(SciMAF%scAngleG < 0.0) .OR. (SciMAF(0)%scAngleG > 0.0)) THEN
             IF (SciMAF(0)%secTAI >= TAI_range%startTime) EXIT
          ENDIF
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
          WRITE (unit, *) 'MAFs missing: ', (MAF_dif-1)
          WRITE (unit, *) 'MAFno gap: ', BeginEnd%SciMAFno(2), SciMAF(0)%MAFno
          PGS_stat = PGS_TD_TAItoUTC (last_TAI, asciiUTC(1))
          PGS_stat = PGS_TD_TAItoUTC (SciMAF(0)%secTAI, asciiUTC(2))
          WRITE (unit, *) 'UTC gap: ', asciiUTC(1)//' to '//asciiUTC(2)
          WRITE (unit, *) ''
       ENDIF

       DO i = 0, last_MIF
         IF (SciMAF(i)%MAFNo /= -1) THEN
            lastGoodMAF=SciMAF(i)%MAFNo
         ENDIF
          IF (SciMAF(i)%MIFno < 0) THEN
             ! It turns out that sciMAF(i)%MAFno will also == -1 here.
             SciGaps = SciGaps + 1
             Sci_Warns = Sci_Warns + 1
             WRITE (unit, *) '##### WARNING! Data Gap:'
             WRITE (unit, *) 'MAF no:',lastGoodMAF,', MIF missing: ', i
             PGS_stat = PGS_TD_TAItoUTC (SciMAF(0)%secTAI+MIF_dur*i,asciiUTC(2))
             WRITE (unit, *) 'UTC: ', asciiUTC(2)
             WRITE (unit, *) ''
             SciMAF(i)%CRC_good = .FALSE.   ! Mark as bad
             SciMAF(i)%FB = 0.0             ! clear the counts
             SciMAF(i)%MB = 0.0
             SciMAF(i)%WF = 0.0
             SciMAF(i)%DACS = 0.0
          ENDIF
       ENDDO

! Save data for further processing:
       !PCF id=920, nominally named sciMAF.dat
       WRITE (L1BFileInfo%SciMAF_unit) SciMAF

       last_TAI = SciMAF(0)%secTAI

       SciMAFs = SciMAFs + 1
       BeginEnd%SciMAFno(2) = SciMAF(0)%MAFno
       BeginEnd%SciTAI(2) = SciMAF(0)%secTAI

       IF (ANY (SciMAF%AttenMaxed)) THEN
          CALL StoreMinAttenCnts (SciMAF)
          doneAttens = .TRUE.
       ENDIF

! Check for Attenuation changes:

       DO i = 1, FBNUM
          IF (ANY (SciMAF%DeltaAtten%FB(i))) THEN
             MIFmask = .FALSE.
             WHERE (SciMAF%DeltaAtten%FB(i))
                MIFmask = .TRUE.
             ENDWHERE
             mindx = MINVAL (MIFindx, MIFmask)
             Sci_Warns = Sci_Warns + 1
             WRITE (unit, *) '##### WARNING! Attenuation Change:'
             WRITE (unit, *) 'FBno: ', i
             WRITE (unit, *) 'MAFno, MIFno: ', SciMAF(0)%MAFno, (mindx-1)
             PGS_stat = PGS_TD_TAItoUTC (SciMAF(mindx-1)%secTAI, asciiUTC(1))
             WRITE (unit, *) 'UTC: ', asciiUTC(1)
             WRITE (unit, *) ''
          ENDIF
       ENDDO

       DO i = 1, MBNUM
          IF (ANY (SciMAF%DeltaAtten%MB(i))) THEN
             MIFmask = .FALSE.
             WHERE (SciMAF%DeltaAtten%MB(i))
                MIFmask = .TRUE.
             ENDWHERE
             mindx = MINVAL (MIFindx, MIFmask)
             Sci_Warns = Sci_Warns + 1
             WRITE (unit, *) '##### WARNING! Attenuation Change:'
             WRITE (unit, *) 'MBno: ', i
             WRITE (unit, *) 'MAFno, MIFno: ', SciMAF(0)%MAFno, (mindx-1)
             PGS_stat = PGS_TD_TAItoUTC (SciMAF(mindx-1)%secTAI, asciiUTC(1))
             WRITE (unit, *) 'UTC: ', asciiUTC(1)
             WRITE (unit, *) ''
          ENDIF
       ENDDO

       DO i = 1, WFNUM
          IF (ANY (SciMAF%DeltaAtten%WF(i))) THEN
             MIFmask = .FALSE.
             WHERE (SciMAF%DeltaAtten%WF(i))
                MIFmask = .TRUE.
             ENDWHERE
             mindx = MINVAL (MIFindx, MIFmask)
             Sci_Warns = Sci_Warns + 1
             WRITE (unit, *) '##### WARNING! Attenuation Change:'
             WRITE (unit, *) 'WFno: ', i
             WRITE (unit, *) 'MAFno, MIFno: ', SciMAF(0)%MAFno, (mindx-1)
             PGS_stat = PGS_TD_TAItoUTC (SciMAF(mindx-1)%secTAI, asciiUTC(1))
             WRITE (unit, *) 'UTC: ', asciiUTC(1)
             WRITE (unit, *) ''
          ENDIF
       ENDDO

       DO i = 1, DACSNUM
          IF (ANY (SciMAF%DeltaAtten%DACS(i))) THEN
             MIFmask = .FALSE.
             WHERE (SciMAF%DeltaAtten%DACS(i))
                MIFmask = .TRUE.
             ENDWHERE
             mindx = MINVAL (MIFindx, MIFmask)
             Sci_Warns = Sci_Warns + 1
             WRITE (unit, *) '##### WARNING! Attenuation Change:'
             WRITE (unit, *) 'DACSno: ', i
             WRITE (unit, *) 'MAFno, MIFno: ', SciMAF(0)%MAFno, (mindx-1)
             PGS_stat = PGS_TD_TAItoUTC (SciMAF(mindx-1)%secTAI, asciiUTC(1))
             WRITE (unit, *) 'UTC: ', asciiUTC(1)
             WRITE (unit, *) ''
          ENDIF
       ENDDO

    ENDDO

    ENDFILE L1BFileInfo%SciMAF_unit

    IF (Sci_Warns == 0 .AND. Sci_Errs == 0) &
         WRITE (unit, *) '##### No Warnings and no Errors ###'

    IF (doneAttens) CALL SaveDefaultZeroCnts   ! Save to file, if any done

  END SUBROUTINE ExamineSciData

!=============================================================================
  SUBROUTINE StoreMinAttenCnts (SciMAF)
!=============================================================================

    USE MLSL1Common, ONLY: MaxMIFs, FBnum, FBchans, MBnum, MBchans, WFnum, &
         WFchans
    USE L0_sci_tbls, ONLY: Sci_pkt_T

    TYPE (Sci_pkt_T), DIMENSION(0:), INTENT (IN) :: SciMAF

    INTEGER :: i, chan, MIFno
    INTEGER, PARAMETER :: lastMIF = (MaxMIFs - 1)

    DO MIFno = 0, lastMIF
       DO i = 1, FBnum   ! FB
          IF (SciMAF(MIFno)%MaxAtten%FB(i)) THEN
             DO chan = 1, FBchans
                Atten_cnts%FB(chan,i) = &
                     MIN (Atten_cnts%FB(chan,i), REAL(SciMAF(MIFno)%FB(chan,i)))
             ENDDO
          ENDIF
       ENDDO
       DO i = 1, MBnum   ! MB
          IF (SciMAF(MIFno)%MaxAtten%MB(i)) THEN
             DO chan = 1, MBchans
                Atten_cnts%MB(chan,i) = &
                     MIN (Atten_cnts%MB(chan,i), REAL(SciMAF(MIFno)%MB(chan,i)))
             ENDDO
          ENDIF
       ENDDO
       DO i = 1, WFnum   ! WF
          IF (SciMAF(MIFno)%MaxAtten%WF(i)) THEN
             DO chan = 1, WFchans
                Atten_cnts%WF(chan,i) = &
                     MIN (Atten_cnts%WF(chan,i), REAL(SciMAF(MIFno)%WF(chan,i)))
             ENDDO
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE StoreMinAttenCnts

!=============================================================================
  SUBROUTINE SaveDefaultZeroCnts
!=============================================================================

    USE MLSL1Common, ONLY: FBnum, FBchans, MBnum, MBchans, WFnum, WFchans, &
         deflt_zero
    USE MLSPCF1, ONLY : mlspcf_defltzeros_end
    USE SDPToolkit, ONLY: PGS_PC_GetReference, PGSd_PC_FILE_PATH_MAX
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info

    CHARACTER (LEN=PGSd_PC_FILE_PATH_MAX) :: PhysicalFilename
    INTEGER :: bank, chan, i, ios, returnStatus, unit, version
    REAL(R8) :: mid_TAI
    CHARACTER (LEN=*), DIMENSION(10), PARAMETER :: Header = (/ &
"                                                                            ",&
"This section is the header comments section for the default zeros file.     ",&
"                                                                            ",&
"Use a #DATA (starting in column 1) to indicate the beginning of the data.   ",&
"                                                                            ",&
"All lines before the #DATA are comments.                                    ",&
"                                                                            ",&
"Data begins after 1 comment line (beginning with a '#').  Data is free form.",&
"The data order is FB 1-19, MB 1-5 followed by WF 1-3.                       ",&
"                                                                            " &
    /)
    REAL, PARAMETER :: HUGE_F = HUGE(1.0)

    INTEGER, EXTERNAL :: PGS_IO_Gen_Track_LUN

! Open output default zeros file:

    version = 1
    returnStatus = PGS_PC_getReference (mlspcf_defltzeros_end, version, &
          & PhysicalFilename)
    returnStatus = PGS_IO_Gen_Track_LUN (unit, 0)
    OPEN (unit=unit, file=PhysicalFilename, status="REPLACE", &
         FORM="FORMATTED", ACCESS="SEQUENTIAL", iostat=ios)

    IF (ios /= 0) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open Output Default Zeros file: " // PhysicalFilename)
    ENDIF
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened Output Default Zeros file: " // PhysicalFilename)

    mid_TAI = TAI_range%StartTime + (TAI_range%EndTime - TAI_range%StartTime)/2
    returnStatus = PGS_TD_TAItoUTC (mid_TAI, asciiUTC(1))

    WRITE (unit, '(A)') (TRIM(Header(i)),i=1,SIZE(Header))
    WRITE (unit, '(A)') 'Produced for data of '//asciiUTC(1)(3:4)//'/'// &
         asciiUTC(1)(6:7)//'/'//asciiUTC(1)(9:10)
    WRITE (unit, '(/,A)') '#DATA'

    DO bank = 1, FBnum
       DO chan = 1, FBchans
          IF (Atten_cnts%FB(chan,bank) == HUGE_F) Atten_cnts%FB(chan,bank) = &
               deflt_zero%FB(chan,bank)
       ENDDO
       WRITE (unit, "(/,'# FB ', i2,/)") bank
       WRITE (unit, "(13I5)") NINT(Atten_cnts%FB(:,bank))
    ENDDO

    DO bank = 1, MBnum
       DO chan = 1, MBchans
          IF (Atten_cnts%MB(chan,bank) == HUGE_F) Atten_cnts%MB(chan,bank) = &
               deflt_zero%MB(chan,bank)
       ENDDO
       WRITE (unit, "(/,'# MB ', i1,/)") bank
       WRITE (unit, "(11I5)") NINT(Atten_cnts%MB(:,bank))
    ENDDO

    DO bank = 1, WFnum
       DO chan = 1, WFchans
          IF (Atten_cnts%WF(chan,bank) == HUGE_F) Atten_cnts%WF(chan,bank) = &
               deflt_zero%WF(chan,bank)
       ENDDO
       WRITE (unit, "(/,'# WF ', i1,/)") bank
       WRITE (unit, "(4I5)") NINT(Atten_cnts%WF(:,bank))
    ENDDO

    CLOSE (unit)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Closed Output Default Zeros file: " // PhysicalFilename)

  END SUBROUTINE SaveDefaultZeroCnts

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
    IF (SciMAFs > EngMAFs) THEN
       WRITE (unit, *) "L0 input warning-> SciMAFs greater than EngMAFs!"
    ENDIF

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

    CHARACTER(LEN=80) :: message
    INTEGER :: window_MAFs

    ! <whd> L1BFileInfo%LogId points to a file opened in
    ! OpenInit. that file is defined in the PCF file with PCF ID =
    ! 30006. Nominally it's named
    ! MLS-Aura_L1BLOG_<version>_<auraday>.txt </whd>

    unit = L1BFileInfo%LogId

! Determine TAI/MAF range to read:

    window_MAFs = L1Config%Calib%CalWindow
    MAF_dur = L1Config%Calib%MIF_duration * L1Config%Calib%MIFsPerMAF
    MIF_dur = L1Config%Calib%MIF_duration
    TAI_range = L1Config%Input_TAI
    TAI_range%startTime = TAI_range%startTime - ((window_MAFs/2 + 1) * MAF_dur)
    TAI_range%endTime = TAI_range%endTime + (window_MAFs/2 * MAF_dur)
    last_MIF = L1Config%Calib%MIFsPerMAF - 1

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

    IF (eng_warns > 0 .OR. sci_warns > 0) THEN
       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
            "EOS MLS Level 1 log file " //TRIM(L1BFileInfo%LogFilename)//&
            " reports warnings!")
    ENDIF

    IF (SciMAFs > EngMAFs) THEN

       CALL MLSMessage (MLSMSG_Warning, ModuleName, &
            "L0 input warning-> SciMAFs greater than EngMAFs!")
    ENDIF

    IF (eng_errs > 0 .OR. sci_errs > 0) THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "EOS MLS Level 1 log file " //TRIM(L1BFileInfo%LogFilename)//&
            " reports errors!")
    ENDIF

  END SUBROUTINE LogStatus

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE L1LogUtils
!=============================================================================

! $Log$
! Revision 2.24  2018/04/09 22:13:36  whdaffer
! Documentation and reportage
!
! Revision 2.23  2017/02/14 17:55:14  whdaffer
! Resolved conflict. Don't quite know how it came about, but it did.
!
! Revision 2.22  2017/01/05 21:24:12  whdaffer
! A few more tweaks to the last modification.
!
! Revision 2.21  2017/01/05 21:10:03  whdaffer
! Modified ExamineEngData to give some more useful information when
! failures arising from start/end time comparisons larger than
! DiffBeginEndEng occur. Also greatly expanded comments to explain these
! comparisons
!
! Revision 2.20  2016/05/10 20:41:23  mmadatya
! To get the error-checking parameters from the l1 configuration file instead of them being hard-coded into the source code
!
! Revision 2.19  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.18  2016/02/12 20:04:48  pwagner
! Prevent bad level 0 files from causing englog to fill all diskspace
!
! Revision 2.17.4.3  2016/03/04 21:43:26  whdaffer
! Merged PW's Trunk changes into my branch
!
! Revision 2.17.4.2  2016/03/03 18:53:08  whdaffer
! Some comments
!
! Revision 2.17.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.17  2011/01/27 15:35:51  perun
! Only check TAI time for good scAngleG or for first empty MAF
!
! Revision 2.16  2008/02/25 17:18:40  perun
! Clear counts for missing MIF and mark flag as bad.
!
! Revision 2.15  2006/09/26 16:02:00  perun
! Make MAF_dur and MIF_dur public
!
! Revision 2.14  2006/04/05 18:10:15  perun
! Remove unused variables
!
! Revision 2.13  2006/03/24 15:08:45  perun
! Removed startup message to log file
!
! Revision 2.12  2005/10/14 18:41:41  perun
! Expanded start and end times to scan data
!
! Revision 2.11  2005/08/24 15:51:29  perun
! Check for and report missing MIFs in science data
!
! Revision 2.10  2005/08/11 19:05:03  perun
! Change SciMAFs greater than EngMAFs from error to warning
!
! Revision 2.9  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.8  2005/01/25 15:45:21  perun
! Compare SciMAFs against EngMAFs and report error if greater
!
! Revision 2.7  2004/08/12 13:51:50  perun
! Version 1.44 commit
!
! Revision 2.6  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.5  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.4  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.3  2003/09/02 17:10:52  perun
! Change examine eng data exit condition
!
! Revision 2.2  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.1  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
!
