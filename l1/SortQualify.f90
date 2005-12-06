! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE SortQualify ! Sort and qualify the L0 data
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, MaxMIFs, OA_counterMAF, &
       OA_counterIndex
  USE L0_sci_tbls, ONLY: SciMAF
  USE EngTbls, ONLY: EngMAF
  USE Calibration, ONLY: CalWin, MAFdata_T, UpdateCalVectors, WeightsFlags_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SortAndQualify

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (MAFdata_T), POINTER :: CurMAFdata

CONTAINS

!=============================================================================
  SUBROUTINE UpdateCalWindow (more_data, CalWinFull)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, MIFsGHz
    USE TkL1B, ONLY: GHz_GeodAlt, GHz_GeodLat, GHz_BO_stat
    USE L1BOutUtils, ONLY: OutputL1BOA

    LOGICAL, INTENT (OUT) :: more_data
    LOGICAL, INTENT (OUT) :: CalWinFull

    INTEGER :: sci_MAFno, dif_MAFno, i, ios, windx
    INTEGER, SAVE :: prev_MAFno
    TYPE (MAFdata_T) :: EmptyMAFdata, MAFdata
    REAL(r8), SAVE :: prev_secTAI
    REAL :: MAF_dur, MIF_dur
    INTEGER :: nom_MIFs
    TYPE (WeightsFlags_T) :: WeightsFlags
    INTEGER, PARAMETER :: last_GHz_indx = MIFsGHz - 1

    nom_MIFs = L1Config%Calib%MIFsPerMAF

!! Initialize empty MAF to "D"iscard:

    EmptyMAFdata%SciPkt%GHz_sw_pos = "D"
    EmptyMAFdata%SciPkt%THz_sw_pos = "D"
    DO i = 0, (SIZE (EmptyMAFdata%ChanType) - 1)
       EmptyMAFdata%ChanType(i)%FB = "D"
       EmptyMAFdata%ChanType(i)%MB = "D"
       EmptyMAFdata%ChanType(i)%WF = "D"
       EmptyMAFdata%ChanType(i)%DACS = "D"
    ENDDO
    EmptyMAFdata%CalType = .FALSE.   ! Not a calibration type MAF
    EmptyMAFdata%EMAF%MIFsPerMAF = nom_MIFs
    EmptyMAFdata%WeightsFlags%recomp_MAF = .TRUE.

!! Get the next MAF's worth of data:

    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) WeightsFlags
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) EngMAF
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) SciMAF
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

! Do L1BOA for all good data times:

    IF (SciMAF(0)%secTAI >= L1Config%Input_TAI%startTime .AND. &
         SciMAF(0)%secTAI<= L1Config%Input_TAI%endTime) THEN

       MAFdata%SciPkt = SciMAF
       MAFdata%EMAF = EngMAF

!! Update MAFinfo from current MAF

       MAFinfo%startTAI = SciMAF(0)%secTAI
       MAFinfo%MIFsPerMAF = Nom_MIFs  ! need a fixed size for L1BOA uses
       MAFinfo%MIF_dur = MIF_dur
       OA_counterMAF(OA_counterIndex) = MAFdata%EMAF%TotalMAF
       OA_counterIndex = OA_counterIndex + 1  ! for next entry

       CALL OutputL1BOA (MAFdata)

    ENDIF

    sci_MAFno = SciMAF(0)%MAFno

PRINT *, "SCI/ENG MAF: ", sci_MAFno, EngMAF%MAFno

    MIF_dur = (SUM(SciMAF(1:100)%secTAI - SciMAF(0:99)%secTAI)/100.0)
    MAF_dur = MIF_dur * nom_MIFs   !Nominal duration of MAF

    IF (CalWin%current > 0) THEN
       dif_MAFno = sci_MAFno - prev_MAFno
       IF (dif_MAFno < 0) THEN      ! Rolled over
          dif_MAFno = NINT (((SciMAF(0)%secTAI - prev_secTAI - &
               4.0 * MIF_dur) / MAF_dur) + 0.5)
       ENDIF
    ELSE
       dif_MAFno = 1
    ENDIF

    prev_MAFno = sci_MAFno
    prev_secTAI = SciMAF(0)%secTAI
    IF (CalWin%current /= CalWin%size) THEN
       CalWin%current = CalWin%current + dif_MAFno
       IF (CalWin%current > CalWin%size) THEN  ! Beyond the end of the window
          CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, &
               (CalWin%current-CalWin%size), EmptyMAFdata)
          CalWin%current = CalWin%size
       ENDIF
    ELSE
       CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, dif_MAFno, EmptyMAFdata)
    ENDIF

    CurMAFdata => CalWin%MAFdata(CalWin%current)

    CurMAFdata%SciPkt = SciMAF
    CurMAFdata%EMAF = EngMAF
    CurMAFdata%last_MIF = EngMAF%MIFsPerMAF - 1
    CurMAFdata%BandSwitch =  &
     CurMAFdata%SciPkt(CurMAFdata%last_MIF)%BandSwitch
    CalWin%MAFdata(1)%start_index = 0
    CalWin%MAFdata(1)%end_index = CalWin%MAFdata(1)%last_MIF
    DO windx = 2, CalWin%current
       CalWin%MAFdata(windx)%start_index = CalWin%MAFdata(windx-1)%end_index + 1
       CalWin%MAFdata(windx)%end_index = CalWin%MAFdata(windx)%start_index + &
            CalWin%MAFdata(windx)%last_MIF
    ENDDO
    CurMAFdata%WeightsFlags = WeightsFlags

! Save GHz Alt, Lat for use in baseline:

    CurMAFdata%SciPkt(0:last_GHz_indx)%altG = GHz_GeodAlt
    CurMAFdata%SciPkt(0:last_GHz_indx)%latG = GHz_GeodLat

! Save GHz BO_stat for further tests:

    CurMAFdata%BO_stat = GHz_BO_stat

!! Update MAFinfo from central MAF

    MAFinfo%startTAI = CalWin%MAFdata(CalWin%central)%SciPkt(0)%secTAI
    MAFinfo%MIFsPerMAF = Nom_MIFs  ! need a fixed size for L1BOA uses
    MAFinfo%MIF_dur = MIF_dur
    MAFinfo%integTime = MIF_dur - L1Config%Calib%MIF_DeadTime

!! Determine if CalWin is full

    IF (CalWin%current == CalWin%size .AND. &
         MAFinfo%startTAI >= L1Config%Input_TAI%startTime) THEN
       CalWinFull = .TRUE.
    ELSE
       CalWinFull = .FALSE.
    ENDIF

    more_data =  MAFinfo%startTAI <= L1Config%Input_TAI%endTime

  END SUBROUTINE UpdateCalWindow

!=============================================================================
  SUBROUTINE QualifyCurrentMAF
!=============================================================================

    USE MLSL1Config, ONLY: GHz_seq, GHz_seq_use, THz_seq, THz_seq_use
    USE MLSL1Common, ONLY: SwitchBank
    USE MLSL1Rad, ONLY: BandToBanks
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC
    USE OutputL1B_DataTypes, ONLY: lenG
    USE MLSStrings, ONLY: Capitalize
    USE BrightObjects_m, ONLY: Test_BO_stat, BO_Match
 
!! Qualify the Current MAF data in the cal window

    CHARACTER(len=1) :: GHz_sw_pos, THz_sw_pos
    CHARACTER(len=80) :: msg
    CHARACTER(len=27) :: asciiUTC
    INTEGER :: MIF, last_MIF, i, band(2), bank(2), n, bno, stat, sw
    LOGICAL :: bandMask(150)   ! large enough for the maximum MIFs
    CHARACTER(len=1), PARAMETER :: TargetType = "T" !! Primary target type
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)

    CurMAFdata => CalWin%MAFdata(CalWin%current)
PRINT *, 'Data:', CurMAFdata%SciPkt%GHz_sw_pos

!! Initialize MIF Rad Precision signs:

    CurMAFdata%MIFprecSign = 1.0

    DO MIF = 0, MaxMIFno !! Check each packet

!! Initialize to "U"ndefined:

       CurMAFdata%ChanType(MIF)%FB = undefined
       CurMAFdata%ChanType(MIF)%MB = undefined
       CurMAFdata%ChanType(MIF)%WF = undefined
       CurMAFdata%ChanType(MIF)%DACS = undefined

!! Rule #1: Data quality information:

       IF (.NOT. CurMAFdata%SciPkt(MIF)%CRC_good) THEN

          !! Discard all

          CurMAFdata%ChanType(MIF)%FB = discard
          CurMAFdata%ChanType(MIF)%MB = discard
          CurMAFdata%ChanType(MIF)%WF = discard
          CurMAFdata%ChanType(MIF)%DACS = discard
          CYCLE                              !! All done for this packet
       ENDIF

!! Rule #2: User Input qualifications:

       !! Set the appropriate user input channels to "D"iscard

!! Rule #3: Disqualify based on engineering ("OFF" state, out-of-lock, etc.)

       !! Set the appropriate channels to "D"iscard

!! Rule #4: Check for "Z"ero data

       !! Set the appropriate bands to a "wall"

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%FB .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%FB)
          CurMAFdata%BankWall%FB = .TRUE.
       ENDWHERE

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%MB .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%MB)
          CurMAFdata%BankWall%MB = .TRUE.
       ENDWHERE

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%WF .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%WF)
          CurMAFdata%BankWall%WF = .TRUE.
       ENDWHERE

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%DACS .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%DACS)
          CurMAFdata%BankWall%DACS = .TRUE.
       ENDWHERE

!! Rule #5: Set Switching Mirror position

       !! Save current sw pos:
       
       GHz_sw_pos = CurMAFdata%SciPkt(MIF)%GHz_sw_pos
       IF (GHz_seq_use == match) GHz_sw_pos = Capitalize (GHz_sw_pos)
       THz_sw_pos = CurMAFdata%SciPkt(MIF)%THz_sw_pos

       !! Possibly use sequence from configuration:

       IF (GHz_seq_use == match) THEN           ! Type Match
          IF (GHz_sw_pos /= GHz_seq(MIF)) THEN  ! "D"iscard if not a match
             GHz_sw_pos = discard
          ENDIF
       ELSE IF (GHz_seq_use == override) THEN   ! Type Override
          IF (GHz_sw_pos /= discard) THEN       ! Override if not a discard
             GHz_sw_pos = GHz_seq(MIF)
          ENDIF
       ENDIF

       IF (THz_seq_use == match) THEN           ! Type Match
          IF (THz_sw_pos /= THz_seq(MIF)) THEN  ! "D"iscard if not a match
             THz_sw_pos = discard
          ENDIF
       ELSE IF (THz_seq_use == override) THEN   ! Type Override
          IF (THz_sw_pos /= discard) THEN       ! Override if not a discard
             THz_sw_pos = THz_seq(MIF)
          ENDIF
       ENDIF

       !! GHz Module:

       !! Filterbanks 1 through 14

       WHERE (CurMAFdata%ChanType(MIF)%FB(:,1:14) == undefined)
          CurMAFdata%ChanType(MIF)%FB(:,1:14) = GHz_sw_pos
          !! NOTE: The THz module could be using FB 12!
       END WHERE
       WHERE (CurMAFdata%ChanType(MIF)%MB == undefined)
          CurMAFdata%ChanType(MIF)%MB = GHz_sw_pos
       END WHERE
       WHERE (CurMAFdata%ChanType(MIF)%WF == undefined)
          CurMAFdata%ChanType(MIF)%WF = GHz_sw_pos
       END WHERE
       WHERE (CurMAFdata%ChanType(MIF)%DACS == undefined)
          CurMAFdata%ChanType(MIF)%DACS = GHz_sw_pos
       END WHERE

       !! Upcase "T"arget types for matching

       IF (GHz_seq_use == match) THEN           ! Type Match

          WHERE (CurMAFdata%ChanType(MIF)%FB(:,1:14) == "t")
             CurMAFdata%ChanType(MIF)%FB(:,1:14) = "T"
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%MB == "t")
             CurMAFdata%ChanType(MIF)%MB = "T"
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%WF == "t")
             CurMAFdata%ChanType(MIF)%WF = "T"
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%DACS == "t")
             CurMAFdata%ChanType(MIF)%DACS = "T"
          END WHERE
       ENDIF

       CurMAFdata%SciPkt(MIF)%GHz_sw_pos = GHz_sw_pos

    ENDDO

! Check if MAF data contains calibration data ("S"pace and "T"arget views):

    IF (ANY (CurMAFdata%SciPkt%GHz_sw_pos == "S" .AND. &
         ANY (CurMAFdata%SciPkt%GHz_sw_pos == "T"))) THEN
       CurMAFdata%CalType = .TRUE.
    ELSE
       CurMAFdata%CalType = .FALSE.
    ENDIF

! Check for any walls flagged in MAF: 

    IF (ANY (CurMAFdata%BankWall%FB) .OR. ANY (CurMAFdata%BankWall%MB) .OR. &
         ANY (CurMAFdata%BankWall%WF) .OR. ANY (CurMAFdata%BankWall%DACS)) THEN
       n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       msg = 'Attenuation change'
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
    ENDIF

PRINT *, 'Sort:', CurMAFdata%ChanType(0:149)%FB(1,1)

!! Rule #6: Discard based on other qualifications such as commanded "W"alls

!! Check for bright objects in Space FOV
!! Will decide later how to handle Space Temperature

    IF (ANY (BTEST (CurMAFdata%BO_stat, 0))) THEN    ! Bit 0 for Space FOV
       msg = 'MOON in Space View'
       n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       CALL MLSMessage (MLSMSG_Warning, ModuleName, TRIM(msg)//' at '//asciiUTC)
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
       CurMAFdata%BankWall%FB = .TRUE.
       CurMAFdata%BankWall%MB = .TRUE.
       CurMAFdata%BankWall%WF = .TRUE.
       CurMAFdata%BankWall%DACS = .TRUE.
    ENDIF

!! Check for bright objects in Limb FOV and mark as "D"iscards via precisions

    CALL Test_BO_stat (CurMAFdata%BO_stat)

    DO n = 1, BO_Match%Num
       msg = TRIM(BO_Match%Name(n)) // ' in Limb View'
       stat = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       CALL MLSMessage (MLSMSG_Warning, ModuleName, TRIM(msg)//' at '//asciiUTC)
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       WHERE (BO_Match%InFOV(:,n))
          CurMAFdata%MIFprecSign = BO_Match%PrecScale(n)   ! sign of precision
       ENDWHERE
    ENDDO

!! Initialize Wall MIFs to beginning of MAF

    CurMAFdata%WallMIF%FB = 0
    CurMAFdata%WallMIF%MB = 0
    CurMAFdata%WallMIF%WF = 0
    CurMAFdata%WallMIF%DACS = 0

!! Check DACS 1 for switch change:

    last_MIF = CurMAFdata%last_MIF
    IF (ANY (CurMAFdata%SciPkt(1:last_MIF)%BandSwitch(1) /= &
         CurMAFdata%SciPkt(0)%BandSwitch(1))) THEN
       CurMAFdata%BankWall%DACS(1) = .TRUE.
       n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       msg = 'DACS 1 switch change'
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
    ENDIF

!! Check FBs for switch changes:

    DO i = 2, 5
       bandMask = .FALSE.
       IF (ANY (CurMAFdata%SciPkt(1:last_MIF)%BandSwitch(i) /= &
            CurMAFdata%SciPkt(0)%BandSwitch(i))) THEN
          CurMAFdata%BankWall%FB(SwitchBank(i)) = .TRUE.
          n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
          WRITE (msg, '("FB switch ", i1, " change")') i
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
          WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
PRINT *, 'switch MAF: ', CurMAFdata%SciPkt(0)%MAFno
          WHERE (CurMAFdata%SciPkt%BandSwitch(i) > 0)
             bandMask = .TRUE.  ! contains real band numbers
          ENDWHERE
          band(1) = MINVAL (CurMAFdata%SciPkt%BandSwitch(i), bandMask)
          band(2) = MAXVAL (CurMAFdata%SciPkt%BandSwitch(i), bandMask)
          DO n = 1, 2
             CALL BandToBanks (band(n), bank)
             DO bno = 1, 2
                IF (bank(bno) /= SwitchBank(i)) THEN   !current already done
                   IF (ANY (SwitchBank(2:) == bank(bno))) THEN
                      DO sw = 2, 5
                         IF (CurMAFdata%SciPkt(0)%BandSwitch(sw) == band(n)) &
                              THEN
                            CurMAFdata%BankWall%FB(bank(bno)) = .TRUE.
                            EXIT
                         ENDIF
                      ENDDO
                   ELSE
                      CurMAFdata%BankWall%FB(bank(bno)) = .TRUE.
                  ENDIF
                ENDIF
             ENDDO
          ENDDO
       ENDIF
    ENDDO

  END SUBROUTINE QualifyCurrentMAF

!=============================================================================
  SUBROUTINE QualifyWindow
!=============================================================================

  USE MLSL1Common, ONLY: MBnum, WFnum, DACSnum, GHzNum

!! Qualify the calibration window for spikes, walls, etc.

    INTEGER, PARAMETER :: MaxWin = 10
    INTEGER :: i, indx(MaxWin) = (/ (i, i=1, MaxWin) /), wallindx(MaxWin)
    INTEGER :: bank, MIF, minwall, maxwall
    INTEGER :: cal_range(2), central, current

!! Update start/end indexes of each MAF in the calibration window

    current = CalWin%current

    cal_range(1) = CalWin%MAFdata(1)%start_index
    cal_range(2) = CalWin%MAFdata(current)%end_index
    central = CalWin%central

    DO bank = 1, GHzNum

       wallindx = 0

       WHERE (CalWin%MAFdata%BankWall%FB(bank))
          wallindx(1:current) = indx(1:current)
       END WHERE

       IF (ANY (wallindx(1:current) /= 0)) THEN

          minwall = MINVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%FB(bank))
          maxwall = MAXVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%FB(bank))
          DO i = minwall, maxwall
             CalWin%MAFdata(i)%BankWall%FB(bank) = .TRUE.
             DO MIF = 0, CalWin%MAFdata(i)%last_MIF
                IF (MIF >= CalWin%MAFdata(i)%WallMIF%FB(bank)) THEN
                   CalWin%MAFdata(i)%ChanType(MIF)%FB(:,bank) = "D"
                ENDIF
             ENDDO
          ENDDO
          IF (minwall >= central) THEN
             CalWin%MAFdata(central)%BankCalInd%FB(bank) = (/ cal_range(1), &
                  (CalWin%MAFdata(minwall-1)%end_index + &
                  CalWin%MAFdata(minwall)%WallMIF%FB(bank)) /)
          ELSE IF (maxwall < current) THEN
             CalWin%MAFdata(central)%BankCalInd%FB(bank) = (/ &
                  CalWin%MAFdata(maxwall+1)%start_index, &
                  CalWin%MAFdata(current)%end_index /)
          ELSE
             CalWin%MAFdata(central)%BankCalInd%FB(bank) = (/ 0, 0 /)
          ENDIF
PRINT *, 'bank, wall: ', bank, CalWin%MAFdata%BankWall%FB(bank)
       ELSE
          CalWin%MAFdata(central)%BankCalInd%FB(bank) = cal_range
       ENDIF
    ENDDO

    DO bank = 1, MBnum
       wallindx = 0
       WHERE (CalWin%MAFdata%BankWall%MB(bank))
          wallindx(1:current) = indx(1:current)
       END WHERE
       IF (ANY (wallindx(1:current) /= 0)) THEN
          minwall = MINVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%MB(bank))
          maxwall = MAXVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%MB(bank))
          DO i = minwall, maxwall
             CalWin%MAFdata(i)%BankWall%MB(bank) = .TRUE.
             DO MIF = 0, CalWin%MAFdata(i)%last_MIF
                IF (MIF >= CalWin%MAFdata(i)%WallMIF%MB(bank)) THEN
                   CalWin%MAFdata(i)%ChanType(MIF)%MB(:,bank) = "D"
                ENDIF
             ENDDO
          ENDDO
          IF (minwall >= central) THEN
             CalWin%MAFdata(central)%BankCalInd%MB(bank) = (/ cal_range(1), &
                  (CalWin%MAFdata(minwall-1)%end_index + &
                  CalWin%MAFdata(minwall)%WallMIF%MB(bank)) /)
          ELSE IF (maxwall < current) THEN
             CalWin%MAFdata(central)%BankCalInd%MB(bank) = (/ &
                  CalWin%MAFdata(maxwall+1)%start_index, &
                  CalWin%MAFdata(current)%end_index /)
          ELSE
             CalWin%MAFdata(central)%BankCalInd%MB(bank) = (/ 0, 0 /)
          ENDIF
       ELSE
          CalWin%MAFdata(central)%BankCalInd%MB(bank) = cal_range
       ENDIF
    ENDDO

    DO bank = 1, WFnum
       wallindx = 0
       WHERE (CalWin%MAFdata%BankWall%WF(bank))
          wallindx(1:current) = indx(1:current)
       END WHERE
       IF (ANY (wallindx(1:current) /= 0)) THEN
          minwall = MINVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%WF(bank))
          maxwall = MAXVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%WF(bank))
          DO i = minwall, maxwall
             CalWin%MAFdata(i)%BankWall%WF(bank) = .TRUE.
             DO MIF = 0, CalWin%MAFdata(i)%last_MIF
                IF (MIF >= CalWin%MAFdata(i)%WallMIF%WF(bank)) THEN
                   CalWin%MAFdata(i)%ChanType(MIF)%WF(:,bank) = "D"
                ENDIF
             ENDDO
          ENDDO
          IF (minwall >= central) THEN
             CalWin%MAFdata(central)%BankCalInd%WF(bank) = (/ cal_range(1), &
                  (CalWin%MAFdata(minwall-1)%end_index + &
                  CalWin%MAFdata(minwall)%WallMIF%WF(bank)) /)
         ELSE IF (maxwall < current) THEN
             CalWin%MAFdata(central)%BankCalInd%WF(bank) = (/ &
                  CalWin%MAFdata(maxwall+1)%start_index, &
                  CalWin%MAFdata(current)%end_index /)
          ELSE
             CalWin%MAFdata(central)%BankCalInd%WF(bank) = (/ 0, 0 /)
          ENDIF
       ELSE
          CalWin%MAFdata(central)%BankCalInd%WF(bank) = cal_range
       ENDIF
    ENDDO

    DO bank = 1, DACSnum
       wallindx = 0
       WHERE (CalWin%MAFdata%BankWall%DACS(bank))
          wallindx(1:current) = indx(1:current)
       END WHERE
       IF (ANY (wallindx(1:current) /= 0)) THEN
          minwall = MINVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%DACS(bank))
          maxwall = MAXVAL (wallindx(1:current), &
               CalWin%MAFdata%BankWall%DACS(bank))
          DO i = minwall, maxwall
             CalWin%MAFdata(i)%BankWall%DACS(bank) = .TRUE.
             DO MIF = 0, CalWin%MAFdata(i)%last_MIF
                IF (MIF >= CalWin%MAFdata(i)%WallMIF%DACS(bank)) THEN
                   CalWin%MAFdata(i)%ChanType(MIF)%DACS(:,bank) = "D"
                ENDIF
             ENDDO
          ENDDO
          IF (minwall >= central) THEN
             CalWin%MAFdata(central)%BankCalInd%DACS(bank) = (/ cal_range(1), &
                  (CalWin%MAFdata(minwall-1)%end_index + &
                  CalWin%MAFdata(minwall)%WallMIF%DACS(bank)) /)
          ELSE IF (maxwall < current) THEN
             CalWin%MAFdata(central)%BankCalInd%DACS(bank) = (/ &
                  CalWin%MAFdata(maxwall+1)%start_index, &
                  CalWin%MAFdata(current)%end_index /)
          ELSE
             CalWin%MAFdata(central)%BankCalInd%DACS(bank) = (/ 0, 0 /)
          ENDIF
PRINT *, 'DACS, wall: ', bank, CalWin%MAFdata%BankWall%DACS(bank)
       ELSE
          CalWin%MAFdata(central)%BankCalInd%DACS(bank) = cal_range
       ENDIF
    ENDDO

  END SUBROUTINE QualifyWindow

!=============================================================================
  SUBROUTINE SortAndQualify (more_data, do_calib)
!=============================================================================

    LOGICAL, INTENT (OUT) :: more_data
    LOGICAL, INTENT (OUT) :: do_calib

    LOGICAL :: CalWinFull = .FALSE.

    do_calib = .FALSE.

    CALL UpdateCalWindow (more_data, CalWinFull)

    IF (.NOT. more_data) RETURN

    CALL QualifyCurrentMAF

    IF (CalWinFull) CALL QualifyWindow

    CALL UpdateCalVectors

    do_calib = CalWinFull  !! Calibrate if calib window is full

  END SUBROUTINE SortAndQualify

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE SortQualify
!=============================================================================

! $Log$
! Revision 2.18  2005/12/06 19:29:22  perun
! Removed call to Flag_Bright_Objects and added testing BO_stat
!
! Revision 2.17  2005/10/14 18:42:57  perun
! Calculate MIF_dur from science data rather than using CF input
!
! Revision 2.16  2005/10/10 14:29:35  perun
! Correct sorting of target types and determine if MAF is calibration type
!
! Revision 2.15  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.14  2005/01/28 17:06:18  perun
! Get GHz FOV bright object flags
!
! Revision 2.13  2004/11/10 15:43:51  perun
! Add MIFprecSign to deal with Bright Objects; save GHz Alt, Lat for use in
! baseline; save OA_counterMAF
!
! Revision 2.12  2004/08/12 18:44:48  perun
! Moved "Wall" messages outside MIF loop
!
! Revision 2.11  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.10  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.9  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.8  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.7  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.6  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.5  2002/04/04 20:39:15  perun
! Corrected wallindx test.
!
! Revision 2.4  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.3  2001/09/10 16:18:25  perun
! Added CalMAFdata from Calibration module
!
! Revision 2.2  2001/03/05 19:54:41  perun
! Check TAI against input TAI
!
! Revision 2.1  2001/02/23 20:57:10  perun
! Version 0.5 commit
!
