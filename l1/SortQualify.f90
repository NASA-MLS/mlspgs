! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SortQualify ! Sort and qualify the L0 data
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: L1BFileInfo, MAFinfo, MaxMIFs
  USE L0_sci_tbls, ONLY: SciMAF
  USE EngTbls, ONLY: EngMAF
  USE Calibration, ONLY: CalWin, MAFdata_T, UpdateCalVectors, WeightsFlags_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: SortAndQualify

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (MAFdata_T), POINTER :: CurMAFdata

CONTAINS

!=============================================================================
  SUBROUTINE UpdateCalWindow (more_data, CalWinFull)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE TkL1B, ONLY: Flag_Bright_Objects, LOG_ARR1_PTR_T
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
    TYPE (LOG_ARR1_PTR_T) :: Space_BO_Flag(2), Limb_BO_Flag(2)

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

       CALL OutputL1BOA (MAFdata)

    ENDIF

    sci_MAFno = SciMAF(0)%MAFno

print *, "SCI/ENG MAF: ", sci_MAFno, EngMAF%MAFno

    MIF_dur = L1Config%Calib%MIF_duration
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

!! Check for Bright Objects in FOVs

    Space_BO_flag(1)%ptr => CurMAFdata%SpaceView%MoonInFOV
    Space_BO_flag(2)%ptr => CurMAFdata%SpaceView%VenusInFOV
    Limb_BO_flag(1)%ptr => CurMAFdata%LimbView%MoonInFOV
    Limb_BO_flag(2)%ptr => CurMAFdata%LimbView%VenusInFOV

    CALL Flag_Bright_Objects (CurMAFdata%SciPkt%secTAI, &
         CurMAFdata%SciPkt%scAngleG, Space_BO_flag, Limb_BO_flag)

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
 
!! Qualify the Current MAF data in the cal window

    CHARACTER(len=1) :: GHz_sw_pos, THz_sw_pos
    CHARACTER(len=80) :: msg
    CHARACTER(len=27) :: asciiUTC
    INTEGER :: MIF, last_MIF, i, band(2), bank(2), n, bno, sw
    LOGICAL :: bandMask(150)   ! large enough for the maximum MIFs
    LOGICAL :: MoonInSpaceView, VenusInSpaceView, &
               MoonInLimbView, VenusInLimbView
    CHARACTER(len=1), PARAMETER :: TargetType = "T" !! Primary target type
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)

    CurMAFdata => CalWin%MAFdata(CalWin%current)

!! Initialize nominal interpolation flags

    CurMAFdata%Nominal%FB = .TRUE.
    CurMAFdata%Nominal%MB = .TRUE.
    CurMAFdata%Nominal%WF = .TRUE.
    CurMAFdata%Nominal%DACS = .TRUE.
print *, 'Data:', CurMAFdata%SciPkt%GHz_sw_pos

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

       IF (ANY (CurMAFdata%BankWall%FB) .OR. ANY (CurMAFdata%BankWall%MB) .OR. &
         ANY (CurMAFdata%BankWall%WF) .OR. ANY (CurMAFdata%BankWall%DACS)) THEN
          n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
          msg = 'Attenuation change'
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
          WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
       ENDIF

!! Rule #5: Set Switching Mirror position

       !! Save current sw pos:
       
       GHz_sw_pos = CurMAFdata%SciPkt(MIF)%GHz_sw_pos
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

       !! "D"iscard incorrect Target Type:   !TEST!!!

       IF (GHz_sw_pos /= TargetType) THEN    ! Need to set TargetType in CF
          WHERE (CurMAFdata%ChanType(MIF)%FB(:,1:14) == "T")
             CurMAFdata%ChanType(MIF)%FB(:,1:14) = discard
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%MB == "T")
             CurMAFdata%ChanType(MIF)%MB = discard
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%WF == "T")
             CurMAFdata%ChanType(MIF)%WF = discard
          END WHERE
          WHERE (CurMAFdata%ChanType(MIF)%DACS == "T")
             CurMAFdata%ChanType(MIF)%DACS = discard
          END WHERE
       ENDIF

       CurMAFdata%SciPkt(MIF)%GHz_sw_pos = GHz_sw_pos

    ENDDO

print *, 'Sort:', CurMAFdata%ChanType(0:149)%FB(1,1)

!! Rule #6: Discard based on other qualifications such as commanded "W"alls

!! Check for bright objects in Space FOV and mark as "W"alls

    MoonInSpaceView = ANY (CurMAFdata%SpaceView%MoonInFOV)
    IF (MoonInSpaceView) msg = 'Moon in Space View'
    VenusInSpaceView = ANY (CurMAFdata%SpaceView%VenusInFOV)
    IF (VenusInSpaceView) msg = 'Venus in Space View'
    IF (MoonInSpaceView .OR. VenusInSpaceView) THEN   ! Mark all banks as Walls
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

!! Check for bright objects in Limb FOV and mark as "D"iscards

    MoonInLimbView = ANY (CurMAFdata%LimbView%MoonInFOV)
    IF (MoonInLimbView) msg = 'Moon in Limb View'
    VenusInLimbView = ANY (CurMAFdata%LimbView%VenusInFOV)
    IF (VenusInLimbView) msg = 'Venus in Limb View' 
    IF (MoonInLimbView .OR. VenusInLimbView) THEN   ! Discard
       n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       CALL MLSMessage (MLSMSG_Warning, ModuleName, TRIM(msg)//' at '//asciiUTC)
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       DO i = 1, lenG
          IF (CurMAFdata%LimbView%MoonInFOV(i) .OR. &
               CurMAFdata%LimbView%VenusInFOV(i)) THEN
             CurMAFdata%ChanType(i)%FB = discard
             CurMAFdata%ChanType(i)%MB = discard
             CurMAFdata%ChanType(i)%WF = discard
             CurMAFdata%ChanType(i)%DACS = discard
          ENDIF
       ENDDO
    ENDIF

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
print *, 'switch MAF: ', CurMAFdata%SciPkt(0)%MAFno
          WHERE (CurMAFdata%SciPkt%BandSwitch(i) > 0)
             bandMask = .TRUE.  ! contains real band numbers
          ENDWHERE
          band(1) = minval (CurMAFdata%SciPkt%BandSwitch(i), bandMask)
          band(2) = maxval (CurMAFdata%SciPkt%BandSwitch(i), bandMask)
          DO n = 1, 2
             CALL BandToBanks (band(n), bank)
             DO bno = 1, 2
                IF (bank(bno) /= SwitchBank(i)) THEN   !current already done
                   IF (ANY (SwitchBank(2:) == bank(bno))) THEN
                      DO sw = 2, 5
                         IF (CurMAFdata%SciPkt(0)%BandSwitch(sw) == band(n)) &
                              THEN
                            CurMAFdata%BankWall%FB(bank(bno)) = .TRUE.
                            CurMAFdata%Nominal%FB(bank(bno)) = .FALSE.
                            EXIT
                         ENDIF
                      ENDDO
                   ELSE
                      CurMAFdata%BankWall%FB(bank(bno)) = .TRUE.
                      CurMAFdata%Nominal%FB(bank(bno)) = .FALSE.
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

          minwall = minval (wallindx(1:current), &
               CalWin%MAFdata%BankWall%FB(bank))
          maxwall = maxval (wallindx(1:current), &
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
print *, 'bank, wall: ', bank, CalWin%MAFdata%BankWall%FB(bank)
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
          minwall = minval (wallindx(1:current), &
               CalWin%MAFdata%BankWall%MB(bank))
          maxwall = maxval (wallindx(1:current), &
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
          minwall = minval (wallindx(1:current), &
               CalWin%MAFdata%BankWall%WF(bank))
          maxwall = maxval (wallindx(1:current), &
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
          minwall = minval (wallindx(1:current), &
               CalWin%MAFdata%BankWall%DACS(bank))
          maxwall = maxval (wallindx(1:current), &
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
print *, 'DACS, wall: ', bank, CalWin%MAFdata%BankWall%DACS(bank)
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
END MODULE SortQualify
!=============================================================================

! $Log$
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
