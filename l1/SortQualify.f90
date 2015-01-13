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
MODULE SortQualify ! Sort and qualify the L0 data
!=============================================================================

  USE MLSKinds, ONLY: R8
  USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, MaxMIFs, OA_counterMAF, &
       OA_counterIndex, ChanLogical_T, BankInt_T, FBnum, MBnum, WFnum, &
       DACSnum, GHzNum
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

  TYPE (MAFdata_T), POINTER :: CurMAFdata => NULL()

  TYPE (BankInt_T) :: BankWallCnt
  INTEGER, PARAMETER :: BankWallSize = 2    ! Initialize size of Wall
  REAL, PARAMETER :: BaselineAlt = 80.0e3   !!! TEMPORARY (will add to CF?)

CONTAINS

!=============================================================================
  SUBROUTINE UpdateCalWindow (more_data, CalWinFull)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config, MIFsGHz
    USE TkL1B, ONLY: GHz_GeodAlt, GHz_GeodLat, GHz_BO_stat, scGeodAngle
    USE L1BOutUtils, ONLY: OutputL1BOA

    LOGICAL, INTENT (OUT) :: more_data
    LOGICAL, INTENT (OUT) :: CalWinFull

    INTEGER :: sci_MAFno, dif_MAFno, bankno, i, ios, windx
    INTEGER, SAVE :: prev_MAFno
    TYPE (MAFdata_T) :: EmptyMAFdata, MAFdata
    REAL(r8), SAVE :: prev_secTAI
    REAL :: MAF_dur, MIF_dur
    INTEGER :: nom_MIFs
    TYPE (WeightsFlags_T) :: WeightsFlags
    TYPE (ChanLogical_T) :: EmptyLimbAltFlag(0:(MaxMIFs-1))
    LOGICAL, SAVE :: InitBankWall = .TRUE.
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

    DO i = 0, (SIZE (EmptyLimbAltFlag) - 1)
       EmptyLimbAltFlag(i)%FB = .FALSE.
       EmptyLimbAltFlag(i)%MB = .FALSE.
       EmptyLimbAltFlag(i)%WF = .FALSE.
       EmptyLimbAltFlag(i)%DACS = .FALSE.
    ENDDO

! Take care of BankWalls:

    IF (InitBankWall) THEN   ! Init Wall counters
       BankWallCnt%FB = 0
       BankWallCnt%MB = 0
       BankWallCnt%WF = 0
       BankWallCnt%DACS = 0
       InitBankWall = .FALSE.
    ENDIF

    DO i = 1, FBnum
       IF (BankWallCnt%FB(i) == 0) THEN
          EmptyMAFdata%BankWall%FB(i) = .FALSE.
       ELSE
          EmptyMAFdata%BankWall%FB(i) = .TRUE.
          BankWallCnt%FB(i) = BankWallCnt%FB(i) - 1    ! decrement until 0
       ENDIF
    ENDDO
    DO i = 1, MBnum
       IF (BankWallCnt%MB(i) == 0) THEN
          EmptyMAFdata%BankWall%MB(i) = .FALSE.
       ELSE
          EmptyMAFdata%BankWall%MB(i) = .TRUE.
          BankWallCnt%MB(i) = BankWallCnt%MB(i) - 1    ! decrement until 0
       ENDIF
    ENDDO
    DO i = 1, WFnum
       IF (BankWallCnt%WF(i) == 0) THEN
          EmptyMAFdata%BankWall%WF(i) = .FALSE.
       ELSE
          EmptyMAFdata%BankWall%WF(i) = .TRUE.
          BankWallCnt%WF(i) = BankWallCnt%WF(i) - 1    ! decrement until 0
       ENDIF
    ENDDO
    DO i = 1, DACSnum
       IF (BankWallCnt%DACS(i) == 0) THEN
          EmptyMAFdata%BankWall%DACS(i) = .FALSE.
       ELSE
          EmptyMAFdata%BankWall%DACS(i) = .TRUE.
          BankWallCnt%DACS(i) = BankWallCnt%DACS(i) - 1    ! decrement until 0
       ENDIF
    ENDDO

! Clear Limb Alts numbers and minimum CalFlags :

    EmptyMAFdata%LimbAltNo%FB = 0
    EmptyMAFdata%LimbAltNo%MB = 0
    EmptyMAFdata%LimbAltNo%WF = 0
    EmptyMAFdata%LimbAltNo%DACS = 0

    EmptyMAFdata%MinCalFlag%FB = .FALSE.
    EmptyMAFdata%MinCalFlag%MB = .FALSE.
    EmptyMAFdata%MinCalFlag%WF = .FALSE.
    EmptyMAFdata%MinCalFlag%DACS = .FALSE.

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

    MIF_dur = L1Config%Calib%MIF_duration
    MAF_dur = MIF_dur * nom_MIFs   !Nominal duration of MAF

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
          CalWin%LimbAltFlag = EOSHIFT (CalWin%LimbAltFlag, &
              (CalWin%current-CalWin%size), EmptyLimbAltFlag, dim=2)
          CalWin%current = CalWin%size
       ENDIF
    ELSE
       CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, dif_MAFno, EmptyMAFdata)
       CalWin%LimbAltFlag = EOSHIFT (CalWin%LimbAltFlag, dif_MAFno, &
            EmptyLimbAltFlag, dim=2)
    ENDIF

    DO i = 1, CalWin%current
       CalWin%MAFdata(i)%LimbAltFlag(0:) => CalWin%LimbAltFlag(0:,i)
    ENDDO

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
!    CurMAFdata%LimbAltFlag => CalWin%LimbAltFlag(0:,CalWin%current)

! Save GHz BO_stat for further tests:

    CurMAFdata%BO_stat = GHz_BO_stat

! Save scGeodAngle for calculating stray radiances:

    CurMAFdata%scGeodAngle = scGeodAngle

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

    USE MLSL1Config, ONLY: GHz_seq, GHz_seq_use, THz_seq, THz_seq_use, L1Config
    USE MLSL1Common, ONLY: SwitchBank, MaxMIFs, L1BFileInfo, MBnum, WFnum, &
         DACSnum, GHzNum
    USE MLSL1Rad, ONLY: BandToBanks
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC
    USE MLSStrings, ONLY: Capitalize
    USE BrightObjects_m, ONLY: Test_BO_stat, BO_Match
    USE BandTbls, ONLY: BandAlt
    USE DACsUtils, ONLY: TPz
 
!! Qualify the Current MAF data in the cal window

    CHARACTER(len=1) :: GHz_sw_pos, THz_sw_pos
    CHARACTER(len=80) :: msg
    CHARACTER(len=27) :: asciiUTC
    INTEGER :: MIF, last_MIF, i, band(2), bandno, bank(2), bankno, n, bno, &
         ngood, stat, sw, sMIF, tMIF, CalDif
    LOGICAL :: bandMask(MaxMIFs)
    CHARACTER(len=1), PARAMETER :: TargetType = "T" !! Primary target type
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)
    REAL :: swFac(0:MaxMIFno)
    REAL(r8) :: TP_ana(0:MaxMIFno), TP_dig(0:MaxMIFno)

    LOGICAL :: TPisDig = .FALSE.

    INTEGER, PARAMETER :: MinCalDif = 100   ! Minimum calibration dif (T - S)

    CurMAFdata => CalWin%MAFdata(CalWin%current)

PRINT *, 'Data:', CurMAFdata%SciPkt%GHz_sw_pos

    TPisDig = L1Config%Calib%TPdigital

    sMIF = -1; tMIF = -1              ! init to nothing so far

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
          !IF (GHz_sw_pos /= discard) THEN       ! Override if not a discard
             GHz_sw_pos = GHz_seq(MIF)
          !ENDIF
       ENDIF

       IF (THz_seq_use == match) THEN           ! Type Match
          IF (THz_sw_pos /= THz_seq(MIF)) THEN  ! "D"iscard if not a match
             THz_sw_pos = discard
          ENDIF
       ELSE IF (THz_seq_use == override) THEN   ! Type Override
          !IF (THz_sw_pos /= discard) THEN       ! Override if not a discard
             THz_sw_pos = THz_seq(MIF)
          !ENDIF
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

! Save last known S and T MIF nos;

       IF (GHz_sw_pos == "S" .and. sMIF < 0) THEN
          sMIF = MIF
       ELSE IF (GHz_sw_pos == "T" .and. tMIF < 0) THEN
          tMIF = MIF
       ENDIF

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

    DO i = 1, FBnum
       IF (ANY(CurMAFdata%SciPkt%MaxAtten%FB(i) .OR. &
            CurMAFdata%SciPkt%DeltaAtten%FB(i))) THEN
          BankWallCnt%FB(i) = BankWallSize
       ENDIF
    ENDDO
    CurMAFdata%SciPkt%GHz_sw_pos = CurMAFdata%ChanType(0:MaxMIFs-1)%FB(1,1)
PRINT *, 'Sort:', CurMAFdata%SciPkt%GHz_sw_pos

! Scale DACS data based on appropriate TP values:

   swFac = 1.0        ! init to non-discard MIFs
   WHERE (CurMAFdata%SciPkt%GHz_sw_pos == discard)   ! Don't use discards
      swFac = 0.0
   ENDWHERE
   ngood = COUNT (swFac == 1.0)

   DO bankno = 1, DACSnum
      TP_ana = CurMAFdata%SciPkt%TP(bankno) * swFac
      TP_dig = (CurMAFdata%SciPkt%TPdigP(bankno) + &
           CurMAFdata%SciPkt%TPdigN(bankno)) * swFac
      IF (TPisDig) THEN
         DO MIF = 0, MaxMIFno
            CurMAFdata%SciPkt(MIF)%DACS(:,bankno) = &
                 CurMAFdata%SciPkt(MIF)%DACS(:,bankno) * TP_dig(MIF)
         ENDDO
      ELSE
!!$         TPz = ((SUM (TP_dig*TP_dig) / ngood * SUM (TP_ana) / ngood) - &
!!$              (SUM (TP_dig*TP_ana) / ngood * SUM (TP_dig) / ngood)) / &
!!$              (SUM (TP_dig*TP_dig) / ngood - (SUM (TP_dig) / ngood)**2)
         DO MIF = 0, MaxMIFno
            CurMAFdata%SciPkt(MIF)%DACS(:,bankno) = &
                 CurMAFdata%SciPkt(MIF)%DACS(:,bankno) * (TP_ana(MIF) - &
                 TPz(bankno))
         ENDDO
      ENDIF

   ENDDO

!! Rule #6: Discard based on other qualifications such as commanded "W"alls

! Check if bank is "OFF", i.e. Space and Target cnts close together:

   IF (sMIF >= 0 .AND. tMIF >= 0) THEN

      DO bankno = 1, GHzNum
         CalDif = CurMAFdata%SciPkt(tMIF)%FB(13,bankno) - &
              CurMAFdata%SciPkt(sMIF)%FB(13,bankno)
         IF (ABS (CalDif) < MinCalDif) THEN
            DO MIF = 0, MaxMIFno
               CurMAFdata%ChanType(MIF)%FB(:,bankno) = discard ! Clear channels
               CurMAFdata%SciPkt(MIF)%FB(:,bankno) = 0
            ENDDO
            CurMAFdata%BankWall%FB(bankno) = .TRUE.
         ENDIF
      ENDDO

      DO bankno = 1, MBnum
         CalDif = CurMAFdata%SciPkt(tMIF)%MB(5,bankno) - &
              CurMAFdata%SciPkt(sMIF)%MB(5,bankno)
         IF (ABS (CalDif) < MinCalDif) THEN
            DO MIF = 0, MaxMIFno
               CurMAFdata%ChanType(MIF)%MB(:,bankno) = discard ! Clear channels
               CurMAFdata%SciPkt(MIF)%MB(:,bankno) = 0
            ENDDO
            CurMAFdata%BankWall%MB(bankno) = .TRUE.
         ENDIF
      ENDDO

      DO bankno = 1, WFnum
         CalDif = CurMAFdata%SciPkt(tMIF)%WF(2,bankno) - &
              CurMAFdata%SciPkt(sMIF)%WF(2,bankno)
         IF (ABS (CalDif) < MinCalDif) THEN
            DO MIF = 0, MaxMIFno
               CurMAFdata%ChanType(MIF)%WF(:,bankno) = discard ! Clear channels
               CurMAFdata%SciPkt(MIF)%WF(:,bankno) = 0
            ENDDO
            CurMAFdata%BankWall%WF(bankno) = .TRUE.
         ENDIF
      ENDDO

      DO bankno = 1, DACSnum
         CalDif = CurMAFdata%SciPkt(tMIF)%DACS(1,bankno) - &
              CurMAFdata%SciPkt(sMIF)%DACS(1,bankno)
         IF (ABS (CalDif) < MinCalDif) THEN
            DO MIF = 0, MaxMIFno
               CurMAFdata%ChanType(MIF)%DACS(:,bankno) = discard ! Clear channels
               CurMAFdata%SciPkt(MIF)%DACS(:,bankno) = 0
            ENDDO
            CurMAFdata%BankWall%DACS(bankno) = .TRUE.
         ENDIF
      ENDDO
   ENDIF

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
       WHERE (BO_Match%InFOV(:,n) .AND. CurMAFdata%MIFprecSign > 0.0)
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

! Set Limb alt flags

    DO MIF = 0, last_MIF !! Check each packet

       DO bankno = 1, GHzNum

          SELECT CASE (bankno)
          CASE (3)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(2)
          CASE (8)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(3)
          CASE (12)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(4)
          CASE default
             bandno = bankno
          END SELECT

          IF (bandno < 0) CYCLE     ! Switch is changing

          WHERE (CurMAFdata%SciPkt(MIF)%altG > BandAlt(bandno)%Meters .AND. &
               (CurMAFdata%MIFprecSign(MIF) > 0.0))
             CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno) = .TRUE.
             WHERE (CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno))
                CurMAFdata%LimbAltNo%FB(:,bankno) = &
                     CurMAFdata%LimbAltNo%FB(:,bankno) + 1
             ENDWHERE
          ELSEWHERE
             CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno) = .FALSE.
          ENDWHERE
 
          CurMAFdata%LimbAltIndx%FB(:,bankno) = BandAlt(bandno)%indx

       ENDDO

       DO bankno = 1, MBnum
          bandno = bankno + 26

          WHERE (CurMAFdata%SciPkt(MIF)%altG > BandAlt(bandno)%Meters .AND. &
               (CurMAFdata%MIFprecSign(MIF) > 0.0))
             CurMAFdata%LimbAltFlag(MIF)%MB(:,bankno) = .TRUE.
             WHERE (CurMAFdata%LimbAltFlag(MIF)%MB(:,bankno))
                CurMAFdata%LimbAltNo%MB(:,bankno) = &
                     CurMAFdata%LimbAltNo%MB(:,bankno) + 1
             ENDWHERE
          ELSEWHERE
             CurMAFdata%LimbAltFlag(MIF)%MB(:,bankno) = .FALSE.
          ENDWHERE
          CurMAFdata%LimbAltIndx%MB(:,bankno) = BandAlt(bandno)%indx

       ENDDO

       DO bankno = 1, WFnum
          bandno = bankno + 31

          WHERE (CurMAFdata%SciPkt(MIF)%altG > BandAlt(bandno)%Meters .AND. &
               (CurMAFdata%MIFprecSign(MIF) > 0.0))
             CurMAFdata%LimbAltFlag(MIF)%WF(:,bankno) = .TRUE.
             WHERE (CurMAFdata%LimbAltFlag(MIF)%WF(:,bankno))
                CurMAFdata%LimbAltNo%WF(:,bankno) = &
                     CurMAFdata%LimbAltNo%WF(:,bankno) + 1
             ENDWHERE
          ELSEWHERE
             CurMAFdata%LimbAltFlag(MIF)%WF(:,bankno) = .FALSE.
          ENDWHERE
          CurMAFdata%LimbAltIndx%WF(:,bankno) = BandAlt(bandno)%indx

       ENDDO

       DO bankno = 1, DACSnum
          IF (bankno == 4) THEN
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(1)
          ELSE
             bandno = bankno + 21
          ENDIF

          IF (bandno < 0) CYCLE     ! Switch is changing

          WHERE (CurMAFdata%SciPkt(MIF)%altG > BandAlt(bandno)%Meters .AND. &
               (CurMAFdata%MIFprecSign(MIF) > 0.0))
             CurMAFdata%LimbAltFlag(MIF)%DACS(:,bankno) = .TRUE.
             WHERE (CurMAFdata%LimbAltFlag(MIF)%DACS(:,bankno))
                CurMAFdata%LimbAltNo%DACS(:,bankno) = &
                     CurMAFdata%LimbAltNo%DACS(:,bankno) + 1
             ENDWHERE
          ELSEWHERE
             CurMAFdata%LimbAltFlag(MIF)%DACS(:,bankno) = .FALSE.
          ENDWHERE
          CurMAFdata%LimbAltIndx%DACS(:,bankno) = BandAlt(bandno)%indx

       ENDDO
    ENDDO

  END SUBROUTINE QualifyCurrentMAF

!=============================================================================
  SUBROUTINE QualifyWindow
!=============================================================================

  USE MLSL1Common, ONLY: MBnum, WFnum, DACSnum, GHzNum, FBchans, MBchans, &
       WFchans, DACSchans
  USE MLSL1Config, ONLY: L1Config

!! Qualify the calibration window for spikes, walls, etc.

    INTEGER, PARAMETER :: MaxWin = 10
    INTEGER :: i, indx(MaxWin) = (/ (i, i=1, MaxWin) /), wallindx(MaxWin)
    INTEGER :: bank, MIF, minwall, maxwall, mincals
    INTEGER :: cal_range(2), central, current

!! Update start/end indexes of each MAF in the calibration window

    current = CalWin%current

    cal_range(1) = CalWin%MAFdata(1)%start_index
    cal_range(2) = CalWin%MAFdata(current)%end_index
    central = CalWin%central
    mincals = L1Config%Calib%MinSpaceLimbs

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

       CalWin%MAFdata(current)%MinCalFlag%FB(:,bank) = .FALSE.
       DO i = 1, FBchans
          IF (ALL (CalWin%MAFdata%LimbAltNo%FB(i,bank) > mincals)) THEN
             CalWin%MAFdata%MinCalFlag%FB(i,bank) = .TRUE.
          ENDIF
       ENDDO

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

       CalWin%MAFdata(current)%MinCalFlag%MB(:,bank) = .FALSE.
       DO i = 1, MBchans
          IF (ALL (CalWin%MAFdata%LimbAltNo%MB(i,bank) > mincals)) THEN
             CalWin%MAFdata%MinCalFlag%MB(i,bank) = .TRUE.
          ENDIF
       ENDDO

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

       CalWin%MAFdata(current)%MinCalFlag%WF(:,bank) = .FALSE.
       DO i = 1, WFchans
          IF (ALL (CalWin%MAFdata%LimbAltNo%WF(i,bank) > mincals)) THEN
             CalWin%MAFdata%MinCalFlag%WF(i,bank) = .TRUE.
          ENDIF
       ENDDO

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
 
       CalWin%MAFdata(current)%MinCalFlag%DACS(:,bank) = .FALSE.
       DO i = 1, DACSchans
          IF (ALL (CalWin%MAFdata%LimbAltNo%DACS(i,bank) > mincals)) THEN
             CalWin%MAFdata%MinCalFlag%DACS(i,bank) = .TRUE.
          ENDIF
       ENDDO

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
! Revision 2.32  2015/01/13 18:42:46  pwagner
! Changed lower bounds on pointer to match LimbAltFlag
!
! Revision 2.31  2013/07/12 15:14:45  perun
! Lowered MinCalDiff to 100 to handle tighter R4 to determine when OFF
!
! Revision 2.30  2009/10/06 16:03:48  perun
! Handle unknown switch position
!
! Revision 2.29  2009/07/24 16:38:43  perun
! Assign MIF_dur before using
!
! Revision 2.28  2009/05/13 20:33:05  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.27  2008/03/28 18:21:22  perun
! Remove debug write statements.
!
! Revision 2.26  2008/03/27 14:40:50  perun
! Calculate calibration dif only when there are both space and target MIFs
!
! Revision 2.25  2008/01/16 19:11:45  perun
! Corrected saving the sorted GHz_sw_pos
!
! Revision 2.24  2007/06/21 21:05:02  perun
! Change mininum calibration dif to 300
!
! Revision 2.23  2006/09/26 16:03:04  perun
! Increase Wall size to 3 MAFs and determine when filter bank is in OFF state
!
! Revision 2.22  2006/08/02 18:58:17  perun
! Remove calculation of TPz for current MAF and use daily TPz data from RADD
!
! Revision 2.21  2006/06/14 13:49:04  perun
! Scale DACS data based on appropriate TP (digital or analog)
!
! Revision 2.20  2006/04/05 18:09:44  perun
! Remove unused variables
!
! Revision 2.19  2006/03/24 15:18:20  perun
! Add sorting based on limb altitudes and remove criteria regarding overriding "discard" MIFs
!
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
