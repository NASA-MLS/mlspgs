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
  USE dates_module, ONLY: tai93s2hid
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info, &
       MLSMSG_Warning

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
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC


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
    INTEGER :: n
    CHARACTER(len=27) :: asciiUTC
    character(len=256) :: msg
    logical :: Deebug
    Deebug = L1Config%Output%DebugUpdateCalWindow

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
    EmptyMAFdata%CalType = .FALSE.   ! Assume this MAF has no calibrations
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

    !!<whd> When a BankWall caused by a change in attenuation is
    !! encountered in QualifyCurrentMAF, a counter is set to 2 for each
    !! channel for which the attenuation changed. This block of code
    !! counts down that counter until we reach 0. Basically, we skip 2
    !! MAFs after an attenuation change </whd>
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

    !! <whd> L1BFileInfo%MAF_data_unit points to a temporary data file written
    !! by MLSL1log which contains the engineering data, the science data and the
    !! flags used in calculating calibration weights. The default name is
    !! MAF_data_tmp.dat (PCF id=922)
    !! </whd>
    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) WeightsFlags
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) EngMAF
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

    READ (unit=L1BFileInfo%MAF_data_unit, iostat=ios) SciMAF
    more_data = (ios == 0)
    IF (.NOT. more_data) RETURN    !! Nothing more to do

    ! <whd> Anything in L1Config comes from the l1cf file </whd>
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

    n = PGS_TD_TAItoUTC (EngMAF%secTAI, asciiUTC)
    if ( Deebug ) then
      WRITE(msg,'("Sci/Eng MAF: ",i4,"/",i4,", UTC: ",a27)') &
           &      sci_MAFno, EngMAF%MAFno, asciiUTC
      print *,TRIM(msg)
      CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
    endif


    IF (CalWin%current > 0) THEN
      dif_MAFno = sci_MAFno - prev_MAFno
      IF (dif_MAFno < 0) THEN      ! Rolled over
        dif_MAFno = NINT (((SciMAF(0)%secTAI - prev_secTAI - &
             4.0 * MIF_dur) / MAF_dur) + 0.5)
      ENDIF
    ELSE
      ! IF ((sci_MAFno - prev_MAFno) > 1) THEN
      !    print *,'Non consecutive MAFs! and no rollover!'
      !    print *,'Pre/This MAF Number = ',prev_MAFno, sci_MAFno
      ! ENDIF
      dif_MAFno = 1
    ENDIF

    prev_MAFno = sci_MAFno
    prev_secTAI = SciMAF(0)%secTAI

    ! <whd> 
    !
    ! Shift the already accumulated data to the left and add more data on
    ! the end. The ammount of the shift depends on possible
    ! discontinuities. It'll be 1 MAF if everything is working correctly, but
    ! more if there's been a discontinuity.
    ! 
    ! </whd>

    IF (CalWin%current /= CalWin%size) THEN
      ! Not full yet!
      CalWin%current = CalWin%current + dif_MAFno
      IF (CalWin%current > CalWin%size) THEN  
        ! Beyond the end of the window. Shift previous MAFs leftward and add
        ! the new data on the end
        CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, &
             (CalWin%current-CalWin%size), EmptyMAFdata)
        CalWin%LimbAltFlag = EOSHIFT (CalWin%LimbAltFlag, &
             (CalWin%current-CalWin%size), EmptyLimbAltFlag, dim=2)
        CalWin%current = CalWin%size
      ENDIF
    ELSE
      ! Drop the first MAF, add new MAF to end.
      CalWin%MAFdata = EOSHIFT (CalWin%MAFdata, dif_MAFno, EmptyMAFdata)
      CalWin%LimbAltFlag = EOSHIFT (CalWin%LimbAltFlag, dif_MAFno, &
           EmptyLimbAltFlag, dim=2)
    ENDIF

    ! <whd> point LimbAlt flag to new locations, given the (possibly) shifted
    ! data. LimbAltFlag is used when doing 'space in limb port' calculations.
    ! </whd>
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

    ! Reset start/end indices
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

    !! Update MAFinfo from central MAF (for L1BOA)

    MAFinfo%startTAI = CalWin%MAFdata(CalWin%central)%SciPkt(0)%secTAI
    MAFinfo%MIFsPerMAF = Nom_MIFs  ! need a fixed size for L1BOA uses
    MAFinfo%MIF_dur = MIF_dur
    MAFinfo%integTime = MIF_dur - L1Config%Calib%MIF_DeadTime

    !! Determine if CalWin is full (needs to be between start/stop times and have
    !! WinMAFs number of MAFs

    IF (CalWin%current == CalWin%size .AND. &
         MAFinfo%startTAI >= L1Config%Input_TAI%startTime) THEN
      n = PGS_TD_TAItoUTC (MAFinfo%startTAI, asciiUTC)
      if ( DeeBug ) then
        write(msg,*) &
             & 'Calibration window full: time of central MAF = '//asciiUTC
        print *,trim(msg)
        CALL MLSMessage(MLSMSG_Info,ModuleName,trim(msg))
        write(msg,*) 'With integration time = ',MAFInfo%integTime
        PRINT *,TRIM(msg)
        CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
      endif

      CalWinFull = .TRUE.
    ELSE
      CalWinFull = .FALSE.
    ENDIF

    more_data =  MAFinfo%startTAI <= L1Config%Input_TAI%endTime

  END SUBROUTINE UpdateCalWindow

!=============================================================================
  SUBROUTINE QualifyCurrentMAF
!=============================================================================


    !! Qualify the Current MAF data in the cal window

    USE MLSL1Config, ONLY: GHz_seq, GHz_seq_use, THz_seq, THz_seq_use, L1Config
    USE MLSL1Common, ONLY: SwitchBank, MaxMIFs, L1BFileInfo, MBnum, WFnum, &
         DACSnum, GHzNum, MAFinfo
    USE MLSL1Rad, ONLY: BandToBanks
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Warning,MLSMSG_Info
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC
    USE MLSStrings, ONLY: Capitalize
    USE BrightObjects_m, ONLY: Test_BO_stat, BO_Match
    USE BandTbls, ONLY: BandAlt
    USE DACsUtils, ONLY: TPz
 
    !! Qualify the Current MAF data in the cal window

    CHARACTER(len=1) :: GHz_sw_pos, THz_sw_pos
    CHARACTER(len=256) :: msg
    CHARACTER(len=27) :: asciiUTC
    INTEGER :: MIF, last_MIF, i, band(2), bandno, bank(2), bankno, n, bno, &
         ngood, stat, sw, sMIF, tMIF, CalDif
    LOGICAL :: bandMask(MaxMIFs) ! true when a change has occured during a MIF
    CHARACTER(len=1), PARAMETER :: TargetType = "T" !! Primary target type
    CHARACTER(len=1), PARAMETER :: discard = "D"
    CHARACTER(len=1), PARAMETER :: match = "M"
    CHARACTER(len=1), PARAMETER :: override = "O"
    CHARACTER(len=1), PARAMETER :: undefined = "U"
    INTEGER, PARAMETER :: MaxMIFno = (MaxMIFs - 1)
    REAL :: swFac(0:MaxMIFno)
    REAL(r8) :: TP_ana(0:MaxMIFno), TP_dig(0:MaxMIFno)

    LOGICAL :: TPisDig = .FALSE.,FoundWall=.FALSE.

    INTEGER, PARAMETER :: MinCalDif = 100   ! Minimum calibration dif (T - S)
    real, parameter    :: OutOfLockScale = -1.0
    character, dimension(2)         :: GM07
    ! Added later to cope with out-of-lock in R2
    integer, parameter              :: BadFirstByte = 13
    integer, parameter              :: BadSecondByte = 202
    integer, parameter              :: GoodFirstByte = 12
    integer, parameter              :: GoodSecondByte = 194

    CurMAFdata => CalWin%MAFdata(CalWin%current)

    
    n = PGS_TD_TAItoUTC (CurMAFData%SciPkt(0)%secTAI, asciiUTC)
    write(msg,'(a,a27,1x,150a1)')'GHz SW Pos Data : ', &
         &                     asciiUTC, &
         &                     CurMAFdata%SciPkt%GHz_sw_pos
    PRINT *, TRIM(msg)
    call MLSMessage(MLSMSG_Info,ModuleName, TRIM(msg))

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

       CurMAFdata%RnPrecSign(MIF,:) = 1. ! Defaults to 1.0
       !! Rule #2: User Input qualifications:

          !! Set the appropriate user input channels to "D"iscard
          ! (paw notes that ..)
          ! Apparently the level 1 omitted Rule #2
          ! in versions prior to ???????????

       !! Rule #3: Disqualify based on engineering ("OFF" state, out-of-lock, etc.)

          !! Set the appropriate channels to "D"iscard
          ! (paw notes that ..)
          ! Apparently the level 1 omitted Rule #3
          ! in versions prior to v6.0x where x > 1
          ! Don't know what "D"iscard is supposed to mean, but
          ! we're going to set the radiance precisions negative
          ! for the affected MIFs if GM07 is out-of-lock
        GM07(1) = CurMAFdata%SciPkt(MIF)%PLL_DN(3:3)
        GM07(2) = CurMAFdata%SciPkt(MIF)%PLL_DN(4:4)
        if ( ichar(GM07(1)) /= GoodFirstByte ) then
          CurMAFdata%RnPrecSign(MIF,2) = OutOfLockScale
        ! elseif ( ichar(GM07(2)) == BadSecondByte ) then
        elseif ( ichar(GM07(2)) /= GoodSecondByte ) then
          CurMAFdata%RnPrecSign(MIF,2) = OutOfLockScale
        endif

       !! Rule #4: Check for "Z"ero data

         !! Set the appropriate bands to a "wall". 

          ! (paw notes that ..)
          ! Having completely ignored Rules #2 and #3, here
          ! the level 1 programmer goes to town with Rule #4
         !! If something has changed that prevents a calibration, we
         !! set a 'wall' that marks this data up to the next 'wall' as
         !! unusable. At the moment, the only thing that can cause
         !! this is the attenuation set at the maximum value of 63
         !! (see MaxAtten, SciUtils.f90::789) or if the attenuation
         !! has changed since the last MAF, of a bright object in the
         !! limb view or some particular banks of particular bands are
         !! 'off'. The following is the test on attenuation.

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%FB .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%FB)
          CurMAFdata%BankWall%FB = .TRUE.
       ENDWHERE
       IF (ANY(CurMAFdata%BankWall%FB)) THEN 
          foundWall=.TRUE.
       ENDIF

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%MB .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%MB)
          CurMAFdata%BankWall%MB = .TRUE.
       ENDWHERE
       IF (ANY(CurMAFdata%BankWall%MB)) THEN 
          foundWall=.TRUE.
       ENDIF

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%WF .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%WF)
          CurMAFdata%BankWall%WF = .TRUE.
       ENDWHERE
       IF (ANY(CurMAFdata%BankWall%WF)) THEN 
          foundWall=.TRUE.
       ENDIF

       WHERE (CurMAFdata%SciPkt(MIF)%MaxAtten%DACS .OR. &
            CurMAFdata%SciPkt(MIF)%DeltaAtten%DACS)
          CurMAFdata%BankWall%DACS = .TRUE.
       ENDWHERE
       IF (ANY(CurMAFdata%BankWall%DACS)) THEN 
          foundWall=.TRUE.
       ENDIF

       IF (foundWall) THEN
          L1Config%Output%NumAttenuationWalls = &
            & L1Config%Output%NumAttenuationWalls + 1
          if ( L1Config%Output%NumAttenuationWalls < &
            & L1Config%Output%MaxAttenuationWalls ) then
            msg='Found Attenuation WALL at MAF time: ' // asciiUTC
            PRINT *,TRIM(msg)
            CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
          elseif ( L1Config%Output%NumAttenuationWalls == &
            & L1Config%Output%MaxAttenuationWalls ) then
            msg='Too many Attenuation WALLS to note all on this date'
            PRINT *,TRIM(msg)
            CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
          endif
       ENDIF

       !! Rule #5: Set Switching Mirror position. 
       
       !! <whd> Here we compare the predicted with the values as read
       !! from telemetry.</whd>

       !! Save current sw pos:

       !! GHz_seq_use is read from the l1cf file, so is a *nominal* value. 
       
       GHz_sw_pos = CurMAFdata%SciPkt(MIF)%GHz_sw_pos
       IF (GHz_seq_use == match) GHz_sw_pos = Capitalize (GHz_sw_pos)
       THz_sw_pos = CurMAFdata%SciPkt(MIF)%THz_sw_pos

       !! Possibly use sequence from configuration:
       
       !! <whd> Here we compare `predicted' (GHz_seq) with observed
       !! (GHz_sw_pos) and change the latter on the basis of the value
       !! of GHz_seq_use. If GHz_seq_use=='override' then we'll use
       !! whatever is in the telemetry.
       !!
       !! IF GHz_seq_use == 'match', the telemetry has to match with
       !! what's in the calibration section of the L1CF file.
       !!
       !! One cause of mismatch could be that the GSM_Theta value that
       !! tells where the GHz Switching mirror is pointing could be in
       !! the 'discard' range, instead of in 'L'imb, 'S'pace, 'T'arget
       !! or 't'arget for some MIF.  One sees this in moontrack runs
       !! because there the ghz mirror is stopped and moved to the
       !! limb port asynchronously with the commanding, so while it
       !! *should* have been an 'S' or 'T', it's moving between any of
       !! the three recognized positions and so gets set to a 'D.' But
       !! that's what it will be set to as a result of the
       !! disagreement with what's predicted by GHz_seq, so this case
       !! is really a no-op.
       !!
       !! The only non-trivial case is when the predicted value is 'D'
       !! and the observed value is != 'D'. Then that MIF of
       !! GHz_sw_pos (and thence CurMAFdata%SciPkt(MIF)%GHz_sw_pos)
       !! will be changed to 'D' because of the mismatch
       !!
       !! This type of disagreement will be communicated in a message
       !! below. The first type will not, because
       !! CurMAFdata%SciPkt(MIF)%GHz_sw_pos will end up equaling
       !! CurMAFdata%ChanType%FB(1,1) because it won't have changed,
       !! despite the disagreement with the predicted behavior.
       !!
       !! </whd>

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
       !
       ! <whd:comment> 
       ! 15:19 are the THz slots in the FB filter bank. 
       ! GHz_sw_pos is determined from the telemetry.
       ! </whd:comment>
       WHERE (CurMAFdata%ChanType(MIF)%FB(:,1:14) == undefined)
          CurMAFdata%ChanType(MIF)%FB(:,1:14) = GHz_sw_pos
          !! vp NOTE: The THz module could be using FB 12!

          ! <whd:comment> 
          ! FB12 can take its input from GHz switch 4. One input to switch 4 is
          ! Band 20, P/T from R5V out of the THz *hardward* module.
          ! </whd:comment> 

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

       !! <whd> GHz_seq_use is set from the l1cf file. It tells what to do when
       !! comparing the expected position of the GHz switching mirror (stored in
       !! L1Config%Calib%GHz_seq) with what's found in telemetry (stored in GHz_sw_pos)
       !! </whd>
       IF (GHz_seq_use == match) THEN           ! Type Match

          ! <whd:comment> 

          ! type `match' means that the telemetry must match what is specified
          ! in the Calibration section of the l1cf file for this run (See
          ! L1Config%Calib%GHz_seq_use). That section specifies the position of
          ! the GHz switching mirror on a MIF by MIF basis and tells which MIFs
          ! are in which category ('L'imb, 'S'pace, 'T'arget and secondary
          ! 't'arget and 'D'iscard, for MIFs that aren't to be used.
          !
          ! 't' signifies that the switching mirror is pointing at the
          ! 'secondary' target and 'T' that it's pointing at the primary
          ! target. The values used in setting GHz_sw_pos (and hence,
          ! ChanType(MIF)%<whatever>. GHz_sw_pos is set in
          ! CalibWeightsFlags::ProcessMAF on a MIF by MIF basis based on what
          ! the GSM_theta angle reads, using ranges as defined by the variables
          ! GHz_SwMir_Range_{A,B,B_2}, THz_SwMir_Range defined in MLSL1Common.
          !
          ! However, there is some problem (which is probably only in my
          ! understanding, but I thought I'd put a note in here anyway) with
          ! which is the 'primary' and 'secondary' target which I'm still
          ! working out. Dominick says (using his terminology, the ranges for
          ! the GSM theta values are ...
          !
          ! 'L'imb     : ~149.5 (L1 calls this the 'L' range)
          ! 'S'pace    : ~329.5 (L1 calls this the 'S' range)
          ! Calibration: 59.5   (L1 calls this the 'T', primary target)
          ! Ambient    : 239.5  (L1 calls this the 't', secondary target, but Dom
          !                      says is never used. However, the telemetry never
          !                      shows an angle around 59.5, only 149, 329 and 239
          !
          ! As you'll see in the next bit of code, when the range reads that
          ! it's 't', it peremptorily sets it to 'T'. I don't quite know why it
          ! does this. I'm told that we've never used the secondary target.
          ! There are no .l1cf files that indicate that any MIF should be using
          ! the secondary target, and even if there were, this code would undo
          ! that. It really appears that the ranges for 'T' and 't' have been
          ! incorrectly entered in MLSL1Common, but rather than fix it there,
          ! they've chosen to 'undo' the effect here by renaming any positions
          ! identified as looking at the 'secondary' target as actually looking
          ! at the primary target.
          !
          ! Presumably that means that, at no time, will the GSM every be
          ! identified as pointing at the 'T' primary target.
          !
          ! </whd:comment>
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

! Save last known 'S'pace and 'T'arget MIF nos;

       IF (GHz_sw_pos == "S" .and. sMIF < 0) THEN
          sMIF = MIF
       ELSE IF (GHz_sw_pos == "T" .and. tMIF < 0) THEN
          tMIF = MIF
       ENDIF

    ENDDO ! loop over MIFs in this MAF

! Check if MAF data contains calibration data (requires at least 1 "S"pace and "T"arget view):

    IF (ANY (CurMAFdata%SciPkt%GHz_sw_pos == "S" .AND. &
         ANY (CurMAFdata%SciPkt%GHz_sw_pos == "T"))) THEN
       CurMAFdata%CalType = .TRUE.
    ELSE
       CurMAFdata%CalType = .FALSE.
    ENDIF

! Check for any walls flagged in MAF: 

    ! L1BFileInfo%LogId points to a file nominally named
    ! MLS-Aura_L1BLOG_<version>_<auraday>.txt in the directory where all output
    ! goes except the STDOUT from T.Sh (if still using that mechanism to run L1)
    IF (ANY (CurMAFdata%BankWall%FB) .OR. ANY (CurMAFdata%BankWall%MB) .OR. &
         ANY (CurMAFdata%BankWall%WF) .OR. ANY (CurMAFdata%BankWall%DACS)) THEN
       n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
       msg = 'Attenuation change'
       WRITE (L1BFileInfo%LogId, *) ''
       WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
       WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC
       PRINT *,TRIM(msg)//' at MAF UTC '//asciiUTC
    ENDIF

    !<whd> BankWallSize=2. In updateCalWindow BankWallCnt%<whatever>(x) is
    !decremented and the wall cleared when it reaches 0. Why this is done only
    !for FB and not MB, WF or DACS is unclear to me </whd>
    
    DO i = 1, FBnum
       IF (ANY(CurMAFdata%SciPkt%MaxAtten%FB(i) .OR. &
            CurMAFdata%SciPkt%DeltaAtten%FB(i))) THEN
         ! If, for this MAF, the bank is at maximum attenuation, or the
         ! attenuation changed, we can't use this data so we 'wall' it off.  I
         ! guess it's set to 2 MAFs because that's the amount of time it takes
         ! for the instrument to settle down after the attenuation is changed.
          BankWallCnt%FB(i) = BankWallSize 
       ENDIF
    ENDDO

    !<whd:comment>
    !
    ! Here we set GHz_sw_pos (which is the location of the GHz switching mirror
    ! in telemetry set by MLSL1log in CalibWeightsFlags::ProcessMAF on the basis
    ! of telemetry processing of SwMirPos(GSM_theta...))  equal
    ! ChanType()%FB. But THAT QUANTITY was just set on the basis of GHz_sw_pos
    ! above! (see lines containing the string
    ! 'CurMAFdata%ChanType(MIF)%FB(:,1:14)'
    !
    ! Curiouser and Curioser!
    !
    ! How this has anything to do with 'sort', I have *no* idea! All that
    ! happens to CurMAFdata%ChanType(...) is that 't's get converted to 'T'. And
    ! why that happens I don't know either.
    !
    !</whd:comment>

    DO MIF=0,MaxMIFs-1 
      IF (CurMAFdata%SciPkt(MIF)%GHz_sw_pos /= CurMAFdata%ChanType(MIF)%FB(1,1)) THEN 
        WRITE(msg,'("GHZ_sw_pos disagreement at MIF ",i3,",Before/After: ",a1,"/",a1 )') &
             & MIF,CurMAFdata%SciPkt(MIF)%GHz_sw_pos, CurMAFdata%ChanType(MIF)%FB(1,1)
        print *,TRIM(msg)
        CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))
      ENDIF
    END DO 
    CurMAFdata%SciPkt%GHz_sw_pos = CurMAFdata%ChanType(0:MaxMIFs-1)%FB(1,1)
    WRITE(msg,'("After GHz_sw_pos processing: GHz sw pos= ",150a1)') &
         &      CurMAFdata%SciPkt%GHz_sw_pos
    PRINT *, TRIM(msg)
    CALL MLSMessage(MLSMSG_Info,ModuleName,TRIM(msg))


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

      ! <whd>: Commenting on the following code... So, apparently if the space
      ! and target views are `close' (less than 100 counts) together for any
      ! 'bank' in FB 13, MB 5, WF 2 or DACCS 1, the filter bank is 'off' and a
      ! "Wall" is declared for all banks. </whd>

      DO bankno = 1, GHzNum
         CalDif = CurMAFdata%SciPkt(tMIF)%FB(13,bankno) - &
              CurMAFdata%SciPkt(sMIF)%FB(13,bankno)
         IF (ABS (CalDif) < MinCalDif) THEN
            DO MIF = 0, MaxMIFno
               CurMAFdata%ChanType(MIF)%FB(:,bankno) = discard ! Clear channels
               CurMAFdata%SciPkt(MIF)%FB(:,bankno) = 0
            ENDDO
            CurMAFdata%BankWall%FB(bankno) = .TRUE.
            n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
            ! msg = 'FB13: bank off at MAF UTC '//asciiUTC
            write(msg,'(a,i3,a,i2,a,a27)') &
                 & 'FB13: bank off, MAF:',&
                 & CurMAFdata%sciPkt(0)%MAFno,&
                 & ', bank_no: ',bankno,',UTC: ', &
                 & TRIM(asciiUTC)
            WRITE (L1BFileInfo%LogId, *) TRIM(msg)
            PRINT *,TRIM(msg)
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
            n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
            write(msg,'(a,i3,a,i2,a,a27)') &
                 & 'MB5: bank off, MAF:',&
                 & CurMAFdata%sciPkt(0)%MAFno,&
                 & ', bank_no: ',bankno,',UTC: ', &
                 & TRIM(asciiUTC) 

            !msg = 'MB5: bank off at MAF UTC '//asciiUTC

            WRITE (L1BFileInfo%LogId, *) TRIM(msg)
            PRINT *,TRIM(msg)

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
            write(msg,'(a,i3,a,i2,a,a27)') &
                 & 'WF2: bank off, MAF:',&
                 & CurMAFdata%sciPkt(0)%MAFno,&
                 & ', bank_no: ',bankno,',UTC: ', &
                 & TRIM(asciiUTC)

            !msg = 'WF2: bank off at MAF UTC '//asciiUTC
            WRITE (L1BFileInfo%LogId, *) TRIM(msg)
            PRINT *,TRIM(msg)

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
            ! msg = 'DACS1: bank off at MAF UTC '//asciiUTC

            write(msg,'(a,i3,a,i2,a,a27)') &
                 & 'DACS1: bank off, MAF:',&
                 & CurMAFdata%sciPkt(0)%MAFno,&
                 & ', bank_no: ',bankno,',UTC: ', &
                 & TRIM(asciiUTC)

            WRITE (L1BFileInfo%LogId, *) TRIM(msg)
            PRINT *,TRIM(msg)
         ENDIF
      ENDDO

   ENDIF ! come from if {s,t}MIF > 0

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

!! <whd> Following code doesn't actually do this. It marks the MIFprecSign
!! instead. Wonder if that causes any problems </whd>


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
    
    ! <whd:comment> 
    !
    ! Wall MIFs: If anything changes in the instrument that would make
    ! calibration impossible, the software declares a 'wall'. At the moment (Wed
    ! Jun 17 2015), Walls are declared in the following circumstances: if the
    ! Attenuation changes or it's maxed out, if the switch changes, if the Moon
    ! or some other bright object is in the limb or space port. The telemetry
    ! captures only changes in attenuation, not the attenuation value itself.
    !
    ! </whd:comment>

    CurMAFdata%WallMIF%FB = 0
    CurMAFdata%WallMIF%MB = 0
    CurMAFdata%WallMIF%WF = 0
    CurMAFdata%WallMIF%DACS = 0

!! Check DACS 1 for switch change:

    ! <whd:comment> 
    !
    ! Looking for changes in the GHz switch which directs some of the bands to
    ! certain sections of the filter banks. Normally, the 'switch network' is
    ! off, Dom says it's only turned on to be commanded and then immediately
    ! turned off. Ongoing processing captures changes in the switch network and
    ! stores it in the file BandSwitches.tbl file (PCF id=913, nominally
    ! BandSwitches.tbl) which is captured by the SIPS in normal Level 1
    ! processing and passed back to the SCF for our offline processing.
    !
    ! Switch 1 is always the DACS (hard wired into the hardware)
    !
    ! </whd:comment>

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
    
    DO i = 2, 5 ! check switches 2 - 5
      bandMask = .FALSE.
      IF (ANY (CurMAFdata%SciPkt(1:last_MIF)%BandSwitch(i) /= &
           CurMAFdata%SciPkt(0)%BandSwitch(i))) THEN
         CurMAFdata%BankWall%FB(SwitchBank(i)) = .TRUE.

         n = PGS_TD_TAItoUTC (CurMAFdata%SciPkt(0)%secTAI, asciiUTC)
         WRITE (msg, '("FB switch ", i1, " change")') i
         WRITE (L1BFileInfo%LogId, *) ''
         WRITE (L1BFileInfo%LogId, *) TRIM(msg)//' at MAF UTC '//asciiUTC
         WRITE (L1BFileInfo%LogId, *) 'WALL event at MAF UTC '//asciiUTC

         PRINT *, 'FB switch change at MAF/Time: ', & 
              & CurMAFdata%SciPkt(0)%MAFno,asciiUTC
         PRINT *, 'WALL event MAF/Time: ', &
              &  CurMAFdata%SciPkt(0)%MAFno,asciiUTC

         WHERE (CurMAFdata%SciPkt%BandSwitch(i) > 0)
            bandMask = .TRUE.  ! These MIFs contains real numbers, i.e. they have
            ! switch (i) 'on'
         ENDWHERE

         ! <whd> find GHz switch positions for the min/max MIFs that have real
         ! numbers in them (i.e. where the switch is 'on') Keep in mind that
         ! BandSwitch is only read from the telemetry when its on. Dom tells me
         ! that normally it's off and reports telem==0 and the values in this
         ! part of the SciPkt user type is actually read from the BandSwitches
         ! file in MLSL1log</whd>

         ! find the min/max bands that are on.
         band(1) = MINVAL (CurMAFdata%SciPkt%BandSwitch(i), bandMask)
         band(2) = MAXVAL (CurMAFdata%SciPkt%BandSwitch(i), bandMask)

         !<whd> Man, I really don't understand what this loop is trying to do!
         !This code seems to assume that BandSwitch has at most 2 possible
         !values for the whole MAF. But if it changed, it's conceivable it
         !could have changed more than once. Maybe it's a fact of the telemetry
         !that only one change can be captured in a MAF?!? Or perhaps it can
         !only be commanded once in a MAF? Have to ask Dom about this.
         !
         ! So this loops over the two 'bands'</whd>

         ! <whd>The upshot is that it declares a bank wall for those bands which have
         ! changed, according to the bandswitches data.</whd>

         DO n = 1, 2

           CALL BandToBanks (band(n), bank) ! find which banks this band is going to.
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
           ENDDO ! Loop over band going to this bank?
         ENDDO ! Loop over min/max MIFs with change
      ENDIF ! check if ANY MIFs had a change in switch
    ENDDO

! Set Limb alt flags

    DO MIF = 0, last_MIF !! Check each packet

      ! <whd>
      ! Another bit of impenetrable code!
      !
      ! This bit of code has the mapping of band to switch to bank.
      !
      ! Look at the diagram named "Aura Microwave Limb Sounders (MLS) -- switch
      ! network configuration". The nominal value of the GHz switch network
      ! (Thu May 28 2015), stored in the bandSwitch variable(s), is
      ! [25,3,8,21,15]. These are the *bands* that go through the switches,
      ! i.e. switch 1 is passing band 25, etc. So, find band 25 in switch
      ! 1. (it's labeled B25D: D=DACS). Now follow the line from the 'out' in
      ! the switch 1 and you'll see that it ends at DACS-1. Similarly, switch 2,
      ! B3F goes to FB25-3; switch 3: B8F goes to FB25-8; switch 4, B21F goes
      ! to FB235-12 and, finally, band 15 goes to FB25-15. 
      !
      ! This loop goes down the filter box (FB) for 19 'banks'. When 'bank'
      ! equals 3, it sets `bandno' to BandSwitch(2) because filter bank 3 takes
      ! it's input frome whatever is being passed through switch 2. Likewise
      ! with switches 3 (band 8 goes to FB-8 through sw 3) and 4 (band 21 goes to
      ! FB-12 through switch 4).
      ! 
      ! Switch 1 is excluded, I guess, because it's hardcoded to DACS-1. I don't
      ! quite know why switch 5 is left out, but it happens that switch 5 is
      ! band 15 which goes to filter bank 15. I certainly hope that's not why
      ! it's left out.
      
      ! Other filter banks don't go through the GHz switch, so 'bandno' just
      ! gets set to 'bankno' because that's where the non-switched bands go, to
      ! FB-xx where 'xx' = `bandno'
      !
      ! </whd>
       
       DO bankno = 1, GHzNum

          SELECT CASE (bankno) ! FB25-[bankno]
          CASE (3)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(2) ! input to that FB num
          CASE (8)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(3)
          CASE (12)
             bandno = CurMAFdata%SciPkt(MIF)%BandSwitch(4)
          CASE default
             bandno = bankno
          END SELECT

          IF (bandno < 0) CYCLE ! Switch is changing, can't do anything when
                                ! that's happening

          ! <whd> LimbAltFlag marks as TRUE the bands/banks for each MIF whose
          ! altitude is an instance of 'Space view in Limb port' (as defined in
          ! the file BandAlt.tbl (PCF id=912) and where precision > 0.
          ! LimbAltNo counts how many MIFs in a MAF satisfy this requirement and
          ! LimbAltIndx%{...}  stores an index into BandAlts::MinAlts which has
          ! minimum altitude for the band/bankno for which 'space in limb port'
          ! condition is true. (I think). Similarly for MB, WF and DACS
          ! </whd>
          
          WHERE (CurMAFdata%SciPkt(MIF)%altG > BandAlt(bandno)%Meters .AND. &
               (CurMAFdata%MIFprecSign(MIF) > 0.0))
             CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno) = .TRUE.
             WHERE (CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno))
                !count # MIFs in MAF above min altitude 
                CurMAFdata%LimbAltNo%FB(:,bankno) = &
                     CurMAFdata%LimbAltNo%FB(:,bankno) + 1
             ENDWHERE
          ELSEWHERE
             CurMAFdata%LimbAltFlag(MIF)%FB(:,bankno) = .FALSE.
          ENDWHERE
          ! %indx marks which locations have acceptible altitudes???
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
       WFchans, DACSchans, FileNameLen
  USE MLSL1Config, ONLY: L1Config
  USE SDPToolkit, ONLY: PGS_TD_TAItoUTC




  !! Qualify the calibration window for spikes, walls, etc.

    INTEGER, PARAMETER :: MaxWin = 20
    INTEGER :: i, indx(MaxWin) = (/ (i, i=1, MaxWin) /), wallindx(MaxWin)
    INTEGER :: bank, MIF, minwall, maxwall, mincals
    INTEGER :: cal_range(2), central, current, n

    CHARACTER(len=27) :: asciiUTC
    CHARACTER(len=FileNameLen) :: msg 


!! Update start/end indexes of each MAF in the calibration window

    current = CalWin%current

    !! <whd> 
    !! cal_range is the MIFs between the start of the WinMAFs cal
    !! window until the end of the current window, so it can be no
    !! longer than maxMIFs*WinMAFs (6*150 for nominal runs, 10*150 for
    !! moontrack runs)
    !! </whd>

    cal_range(1) = CalWin%MAFdata(1)%start_index
    cal_range(2) = CalWin%MAFdata(current)%end_index
    central = CalWin%central
    mincals = L1Config%Calib%MinSpaceLimbs

    DO bank = 1, GHzNum

       wallindx = 0 ! A MAF base quantity

       WHERE (CalWin%MAFdata%BankWall%FB(bank))
          wallindx(1:current) = indx(1:current)
       END WHERE

       IF (ANY (wallindx(1:current) /= 0)) THEN
         !<whd>Find the first/last MAFs that have a Wall declared.</whd>
         minwall = MINVAL (wallindx(1:current), &
              CalWin%MAFdata%BankWall%FB(bank))
         maxwall = MAXVAL (wallindx(1:current), &
              CalWin%MAFdata%BankWall%FB(bank))

         DO i = minwall, maxwall
           !! <whd> And mark all MAFs between first and last as Walls</whd>
           CalWin%MAFdata(i)%BankWall%FB(bank) = .TRUE.
           DO MIF = 0, CalWin%MAFdata(i)%last_MIF
             
             ! <whd> And mark all MIFs > the beginning MIF with a wall
             ! as 'D'. However, I don't think this code ever comes
             ! into play, except in the trivial case where WallMIF==0,
             ! because I can't find any place where
             ! CalWin%MAFdata(...)%WallMIF%<whatever>(bank) is *ever*
             ! set except the initialization that happens above.
             !
             ! So, the net effect is just to make all MIFs as "D". I suspect
             ! this is a bit of nascent code that never got fully fleshed out.
             ! </whd>
             
             IF (MIF >= CalWin%MAFdata(i)%WallMIF%FB(bank)) THEN
               CalWin%MAFdata(i)%ChanType(MIF)%FB(:,bank) = "D"
             ENDIF
           ENDDO ! loop over MIFs inside these MAFs
         ENDDO ! loop from MAF=minWall to MAF=maxWall
         
         ! <whd> 
         ! Don't quite understand this bit. 
         !
         ! BankCalInd%{whatever}(bank)=(beginning,end) of usable data
         ! in the [0,maxMIF] buffer for each bank. This bank is
         ! 150*WinMAFs long and the CalWin%MAF(*).{start,end}_index
         ! point to the beginning/end of each MAF in the calibration
         ! window. 
         !
         ! There never seems to be a case where WallMIF equals
         ! anything other than 0, so the effect of this code is the
         ! following.  
         !
         ! IF the first wall is >= central MAF, set this range to
         ! (beginning_of_cal_window,end_index_of_maf_before_wall). If
         ! the last wall is less than `current' MAF (which should
         ! always be the end of the Cal window, because this routine
         ! isn't called except when the Calibration window is full),
         ! the set the usable range to the beginning of the next MAF
         ! to the end of the window. If minwall < central and maxwall
         ! == current, mark the whole cal window bad.
         !
         ! </whd>
         
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
         n = PGS_TD_TAItoUTC (CalWin%MAFdata(central)%EMAF%secTAI, asciiUTC)
         ! PRINT *, 'FB: UTC, bank, wall, : ', asciiUTC,&
         !      &       bank, CalWin%MAFdata%BankWall%FB(bank)
         WRITE (msg, '("FB UTC(central), bank, wall ", a23, 2x, i2, 2x, 10(L1,:,1x))') asciiUTC,&
              &       bank, CalWin%MAFdata%BankWall%FB(bank)
         
         PRINT *,TRIM(msg)
         WRITE(msg,&
              & '("FB: Walls start at MAF: ", i4, ", end at ",i4,", bank:",i3)') &
              &  CalWin%MAFdata(minwall)%SciPkt(0)%MAFno, &
              &CalWin%MAFdata(maxwall)%SciPkt(0)%MAFno, bank

         Call MLSMessage(MLSMSG_Warning,ModuleName,TRIM(msg))
         PRINT *,TRIM(msg)

         CALL MLSMessage(MLSMSG_Info, ModuleName, TRIM(msg))
         
       ELSE
          CalWin%MAFdata(central)%BankCalInd%FB(bank) = cal_range
       ENDIF

       CalWin%MAFdata(current)%MinCalFlag%FB(:,bank) = .FALSE.
       DO i = 1, FBchans
          IF (ALL (CalWin%MAFdata%LimbAltNo%FB(i,bank) > mincals)) THEN
             CalWin%MAFdata%MinCalFlag%FB(i,bank) = .TRUE.
          ENDIF
       ENDDO
    ENDDO ! over GHz channels (1-14)

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

    ENDDO ! 1,MBNum

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
! Revision 2.35  2024/10/10 20:13:51  pwagner
! Detect out-of-lock for R2 and then set RnPrecSign
!
! Revision 2.34  2023/06/06 22:33:51  pwagner
! Reduce routine printing
!
! Revision 2.33  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.32.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
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
