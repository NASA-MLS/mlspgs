! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SciUtils ! L0 science utilities
!=============================================================================

  USE L0_sci_tbls, ONLY: sci_type, l0_sci1, l0_sci2, sci1_T1_fmt, sci2_T1_fmt, &
       sci1_T2_fmt, sci2_T2_fmt, sci1_T1, sci2_T1, sci1_T2, sci2_T2, sci_cptr, &
       Sci_pkt, SciMAF, type_sci1, THzSciMAF, DACS_pkt, DACS_MAF, Atten_T
  USE L0Utils, ONLY: ReadL0Sci
  USE MLSL1Utils, ONLY: BigEndianStr, ExtractBigEndians, SwapBytes, QNan, &
       Finite
  USE DACsUtils, ONLY: ExtractDACSdata, UncompressDACSdata, ProcessDACSdata
  USE MLSL1Common, ONLY: FBnum, MBnum, WFnum, WFchans, deg24, BankLogical_T, &
       DACSchans, BandSwitch, L1ProgType, THzType, MaxMIFs, NumBands

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: NextSciMAF, SwMirPos, GetScAngles

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  TYPE (BankLogical_T) :: MaxAtten = & ! initialize to NOT Max Atten
       BankLogical_T (.FALSE., .FALSE., .FALSE., .FALSE.)

  CONTAINS

!=============================================================================
  FUNCTION GetSciPkt () RESULT (OK)
!=============================================================================

    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info
    USE ERMSG_M, ONLY: ermset
    USE THzUtils, ONLY: ConvertLLO

    LOGICAL OK

    INTEGER :: i, ios, j, returnStatus
    INTEGER :: tindex   ! Type index: 1 = Type I, 2 = Type II/III

    !! Pointers to counts:

    INTEGER, DIMENSION(:), POINTER :: FBcnts => NULL()
    INTEGER, DIMENSION(:), POINTER :: MBcnts => NULL()
    INTEGER, DIMENSION(:), POINTER :: WFcnts => NULL()
    CHARACTER(len=1), PARAMETER :: sw_good = CHAR(7)
    CHARACTER(len=2) :: WFdat(WFchans)   ! wide filter raw data buffer
    CHARACTER (LEN=1024) :: scipkt(2)
    CHARACTER (LEN=24) :: DN, sw_msg

    TYPE (Atten_T) :: Attenuation
    LOGICAL :: AttenMaxed

    INTEGER :: DACS_C_K(DACSchans)
    INTEGER :: DACS_i1, DACS_i2
    INTEGER :: D(4), TP, DIO, LO, Zlag

    INTEGER, PARAMETER :: type_I = 1
    INTEGER, PARAMETER :: type_II = 3
    INTEGER, PARAMETER :: type_III = 5

    INTEGER, EXTERNAL :: PGS_TD_EOSPMGIRDtoTAI

    CALL ReadL0Sci (scipkt, OK)

    l0_sci1 = scipkt(1)
    l0_sci2 = scipkt(2)

    IF (.NOT. OK) RETURN

    sci_type = ICHAR(type_sci1)  ! science packet #1 type

    IF (sci_type == type_I) THEN

       READ (UNIT=l0_sci1, FMT=sci1_T1_fmt, iostat=ios) sci1_T1
       READ (UNIT=l0_sci2, FMT=sci2_T1_fmt, iostat=ios) sci2_T1
       tindex = 1

    ELSE IF (sci_type == type_II .OR. sci_type == type_III) THEN

       READ (UNIT=l0_sci1, FMT=sci1_T2_fmt, iostat=ios) sci1_T2
       READ (UNIT=l0_sci2, FMT=sci2_T2_fmt, iostat=ios) sci2_T2
       tindex = 2

       !! DACS indexes based on science type format

       IF (sci_type == type_II) THEN
          DACS_i1 = 1
          DACS_i2 = 2
       ELSE
          DACS_i1 = 3
          DACS_i2 = 4
       ENDIF

    ELSE

       PRINT *, 'sci_type: ', sci_type   !! Not correct type!
       OK = .FALSE.
       RETURN

    ENDIF

!! Check for good checksums (LATER!!!):

    Sci_pkt%CRC_good = .TRUE.

!! Convert to angles:

    DN = sci_cptr(tindex)%THz_sw
    IF (DN(1:1) /= sw_good) DN = DN(4:)   ! shift and try
    IF (DN(1:1) == sw_good) THEN
       Sci_pkt%TSSA_pos(1)  = deg24 * BigEndianStr (DN(3:5))
       Sci_pkt%TSSA_pos(2)  = deg24 * BigEndianStr (DN(6:8))
    ELSE
       Sci_pkt%TSSA_pos(:) = QNan()
    ENDIF
    Sci_pkt%THz_sw_pos = SwMirPos ("T", Sci_pkt%TSSA_pos)

    DN = sci_cptr(tindex)%GHz_sw
    IF (DN(1:1) == sw_good) THEN
       Sci_pkt%GSA_pos(1)  = deg24 * BigEndianStr (DN(3:5))
       Sci_pkt%GSA_pos(2)  = deg24 * BigEndianStr (DN(6:8))
    ELSE
       Sci_pkt%GSA_pos(:) = QNan()
    ENDIF
    Sci_pkt%GHz_sw_pos = SwMirPos ("G", Sci_pkt%GSA_pos)

!! Get GHz scanning angles

    Sci_pkt%APE_pos(:) = QNan()
    Sci_pkt%ASA_pos(:) = QNan()
    DN = sci_cptr(tindex)%GHz_ant_scan
    CALL Get_APE_ASA_pos (DN, Sci_pkt%APE_pos, Sci_pkt%ASA_pos)

!! Filter bank switches

    Sci_pkt%GSN = ICHAR (sci_cptr(tindex)%GSN%ptr)
    DO i = 1, 4
       CALL Band_switch (i, Sci_pkt%GSN(i), Sci_pkt%BandSwitch(i))
    ENDDO
    Sci_pkt%THzSw = ICHAR (sci_cptr(tindex)%THzSw%ptr)
    CALL Band_switch (5, Sci_pkt%THzSw, Sci_pkt%BandSwitch(5))

    DO i = 1, 5
       IF (Sci_pkt%BandSwitch(i) == 0) THEN
          Sci_pkt%BandSwitch(i) = BandSwitch(i)  ! use current value
       ELSE IF ((Sci_pkt%BandSwitch(i) > 0) .AND. &
            (BandSwitch(i) /= Sci_pkt%BandSwitch(i))) THEN
          BandSwitch(i) = Sci_pkt%BandSwitch(i)  ! new current value
          WRITE (sw_msg, '("S", i1, " switching to Band ", i2)')i, bandswitch(i)
          PRINT *, sw_msg
          CALL MLSMessage (MLSMSG_Info, ModuleName, sw_msg)
       ENDIF
    ENDDO

!! Convert raw data:

    Sci_pkt%MAFno = BigEndianStr (sci_cptr(tindex)%MAF(1)%ptr)
    Sci_pkt%MIFno = BigEndianStr (sci_cptr(tindex)%MIF(1)%ptr)
    Sci_pkt%Orbit = BigEndianStr (sci_cptr(tindex)%orbit(1)%ptr)

! Get TAI93 time

    returnstatus = PGS_TD_EOSPMGIRDtoTAI (scipkt(1)(8:15), Sci_pkt%secTAI)

!! LLO data

    Sci_pkt%LLO_DN = sci_cptr(tindex)%LLO_DN

    CALL ConvertLLO (Sci_pkt%LLO_DN, Sci_pkt%LLO_EU)

!! PLL data

    CALL Get_PLL_DN (scipkt, sci_type, Sci_pkt%PLL_DN)

!! FB data

    DO i = 1, FBnum
       CALL ExtractBigEndians (sci_cptr(tindex)%FB(i)%ptr, FBcnts)
       Sci_pkt%FB(:,i) = FBcnts
    ENDDO

!! Attenuation readings

    Attenuation%RIU = ICHAR (sci_cptr(tindex)%Attenuation%ptr(1))
    Attenuation%Addr = BigEndianStr ( &
         sci_cptr(tindex)%Attenuation%ptr(2) // &
         sci_cptr(tindex)%Attenuation%ptr(3))
    Attenuation%Mask = ICHAR (sci_cptr(tindex)%Attenuation%ptr(4))
    Attenuation%Value = ICHAR (sci_cptr(tindex)%Attenuation%ptr(5))

!! Determine latest attenuations

    CALL DetermineAttens (Attenuation%RIU, Attenuation%Addr, Attenuation%Value,&
         Sci_pkt%BandSwitch, MaxAtten, AttenMaxed)

    Sci_pkt%MaxAtten = MaxAtten   ! Use latest attenuation flags
    Sci_pkt%AttenMaxed = AttenMaxed

! Nothing more if THz

    IF (L1ProgType == THzType) RETURN

    DO i = 1, MBnum
       CALL ExtractBigEndians (sci_cptr(tindex)%MB(i)%ptr, MBcnts)
       Sci_pkt%MB(:,i) = MBcnts
    ENDDO

    DO i = 1, WFnum
       DO j = 1, WFchans
          WFdat(j) = sci_cptr(tindex)%WF(i)%ptr(j*2:j*2+1)
          CALL SwapBytes (WFdat(j), WFdat(j))
       ENDDO
       CALL ExtractBigEndians (WFdat, WFcnts)
       Sci_pkt%WF(:,i) = WFcnts
    ENDDO

!! DACS data

    CALL ermset (-4)      ! turn OFF error messages

    Sci_pkt%DACS = 0.0
    DACS_pkt%D = 0

    IF (sci_type == type_I) THEN    ! compressed DACS

       DACS_pkt%Compressed = .TRUE.
       DO i = 1, 4

          CALL UncompressDACSdata (sci_cptr(tindex)%DACS(i)%ptr, &
               DACS_C_K(1:128), D, TP, DIO, LO, Zlag)

          DACS_pkt%C_K(:,i) = DACS_C_K
          DACS_pkt%D(:,i) = D
          DACS_pkt%TP(i) = TP
          DACS_pkt%DIO(i) = DIO
          DACS_pkt%LO(i) = LO
          DACS_pkt%Zlag(i) = Zlag

       ENDDO

    ELSE

       DACS_pkt%Compressed = .FALSE.
       DO i = DACS_i1, DACS_i2

          CALL ExtractDACSdata (sci_cptr(tindex)%DACS(i)%ptr, DACS_C_K, D, TP, &
               DIO, LO)

          DACS_pkt%C_K(:,i) = DACS_C_K
          DACS_pkt%D(:,i) = D
          DACS_pkt%TP(i) = TP
          DACS_pkt%DIO(i) = DIO
          DACS_pkt%LO(i) = LO

       ENDDO

    ENDIF

  END FUNCTION GetSciPkt

!=============================================================================
  SUBROUTINE NextSciMAF (more_data)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE MLSL1Common, ONLY : L1ProgType, LogType

    !! Get the next MAF's science data

    LOGICAL, INTENT (OUT) :: more_data

    INTEGER, PARAMETER :: no_data = -1   ! no data is available
    INTEGER, SAVE :: prev_MAF = no_data
    INTEGER :: last_MIF, i
    REAL :: pos1(2)
    CHARACTER(len=1) :: mir_pos

    more_data = .TRUE.

    !! Initialize MAF/MIF counters to indicate no data

    SciMAF%MAFno = no_data
    SciMAF%MIFno = no_data
    last_MIF = L1Config%Calib%MIFsPerMAF -1
    DO i = 0, (MaxMIFs - 1); DACS_MAF(i)%D = 0; ENDDO

    !! Initialize CRC flags to good

    SciMAF%CRC_good = .TRUE.

    !! Save previously read packet (if available):

    IF (prev_MAF /= no_data) THEN
       SciMAF(Sci_pkt%MIFno) = Sci_pkt
       DACS_MAF(Sci_pkt%MIFno) = DACS_pkt
       IF (L1ProgType == THzType) THEN
          CALL Save_THz_pkt (SciMAF(Sci_pkt%MIFno), THzSciMAF(Sci_pkt%MIFno))
       ENDIF
    ENDIF

    DO

       IF (GetSciPkt()) THEN

          IF ((prev_MAF /= no_data) .AND. (Sci_pkt%MAFno /= prev_MAF)) THEN
             prev_MAF = Sci_pkt%MAFno
             IF (SciMAF(0)%MIFno /= no_data .AND. &
                  SciMAF(last_MIF)%MIFno /= no_data) THEN
                EXIT      ! Already got a full MAF's worth
             ELSE
                SciMAF%MAFno = no_data
                SciMAF%MIFno = no_data
             ENDIF
          ENDIF

          SciMAF(Sci_pkt%MIFno) = Sci_pkt  ! save current packet
          DACS_MAF(Sci_pkt%MIFno) = DACS_pkt  ! save DACS data

          IF (L1ProgType == THzType) &
           CALL Save_THz_pkt (SciMAF(Sci_pkt%MIFno), THzSciMAF(Sci_pkt%MIFno))

! Shift mechanism encoder readings to corresponding MIF

          IF (Sci_pkt%MIFno > 0) THEN
             SciMAF(Sci_pkt%MIFno-1)%GSA_pos = &
                  SciMAF(Sci_pkt%MIFno)%GSA_pos
             SciMAF(Sci_pkt%MIFno)%GSA_pos = QNan()
             SciMAF(Sci_pkt%MIFno-1)%APE_pos(2) = &
                  SciMAF(Sci_pkt%MIFno)%APE_pos(2)
             SciMAF(Sci_pkt%MIFno)%APE_pos(2) = QNan()
             SciMAF(Sci_pkt%MIFno-1)%ASA_pos(2) = &
                  SciMAF(Sci_pkt%MIFno)%ASA_pos(2)
             SciMAF(Sci_pkt%MIFno)%ASA_pos(2) = QNan()
             SciMAF(Sci_pkt%MIFno-1)%GHz_sw_pos = &
                  SciMAF(Sci_pkt%MIFno)%GHz_sw_pos
             SciMAF(Sci_pkt%MIFno)%GHz_sw_pos = "D"
             SciMAF(Sci_pkt%MIFno-1)%TSSA_pos = &
                  SciMAF(Sci_pkt%MIFno)%TSSA_pos
             SciMAF(Sci_pkt%MIFno)%TSSA_pos = QNan()

             IF (L1ProgType == THzType) THEN            ! THz if needed
                THzSciMAF(Sci_pkt%MIFno-1)%TSSA_pos = &
                     THzSciMAF(Sci_pkt%MIFno)%TSSA_pos
                THzSciMAF(Sci_pkt%MIFno)%TSSA_pos = QNan()
                THzSciMAF(Sci_pkt%MIFno-1)%SwMirPos = &
                     THzSciMAF(Sci_pkt%MIFno)%SwMirPos
                THzSciMAF(Sci_pkt%MIFno)%SwMirPos = "D"
             ENDIF

          ENDIF

          IF (Sci_pkt%MIFno > 1) THEN
             SciMAF(Sci_pkt%MIFno-2)%APE_pos(1) = &
                  SciMAF(Sci_pkt%MIFno)%APE_pos(1)
             SciMAF(Sci_pkt%MIFno)%APE_pos(1) = QNan()
             SciMAF(Sci_pkt%MIFno-2)%ASA_pos(1) = &
                  SciMAF(Sci_pkt%MIFno)%ASA_pos(1)
             SciMAF(Sci_pkt%MIFno)%ASA_pos(1) = QNan()
           ENDIF

          prev_MAF = Sci_pkt%MAFno

       ELSE

          more_data = .FALSE.
          prev_MAF = no_data
          EXIT  ! Nothing more is available

       ENDIF

    ENDDO

! Compare previous pos1 against current MIF pos1/2:

    IF (L1ProgType == THzType) THEN
       DO i = 1, (MaxMIFs - 1)
          IF (THzSciMAF(i)%SwMirPos /= "D") THEN
             pos1(1) = THzSciMAF(i-1)%TSSA_pos(1)  ! previous pos1
             pos1(2) = pos1(1)   ! need 2 angles to figure sw mir pos
             IF (Finite (pos1(1))) THEN
                mir_pos = SwMirPos ("T", pos1)
                IF (mir_pos /= THzSciMAF(i)%SwMirPos) &
                     THzSciMAF(i)%SwMirPos = "D"
             ENDIF
          ENDIF
       ENDDO
    ENDIF

! Process DACS data for the current MAF

    IF (L1ProgType == LogType) CALL ProcessDACSdata

! Get scAngles to be used for L1BOA

    CALL GetScAngles

  END SUBROUTINE NextSciMAF

!=============================================================================
  FUNCTION SwMirPos (sw_module, angle) RESULT (sw_pos)
!=============================================================================

    USE MLSL1Common, ONLY: GHz_SwMir_Range_A, GHz_SwMir_Range_B, &
       THz_SwMir_Range,  Discard_SwMir_Range, SwMir_Range_T
    USE EngTbls, ONLY: EngMAF

    CHARACTER(len=1), INTENT(IN) :: sw_module
    REAL, INTENT(IN) :: angle(2)
    CHARACTER(len=1) :: sw_pos

    INTEGER :: n
    TYPE (SwMir_Range_T), DIMENSION(:), POINTER :: SwMir_Range

    sw_pos = "D"                 ! Initialize to "Discard"
    IF (.NOT. Finite (angle(1)) .OR. .NOT. Finite (angle(2))) RETURN

    IF (sw_module == "G") THEN   ! GigaHertz Module
       IF (EngMAF%GSM_Side == "A") THEN
          SwMir_Range => GHz_SwMir_Range_A
       ELSE IF (EngMAF%GSM_Side == "B") THEN
          SwMir_Range => GHz_SwMir_Range_B
       ELSE
          SwMir_Range => Discard_SwMir_Range
       ENDIF
    ELSE                         ! TeraHertz Module
       SwMir_Range => THz_SwMir_Range
    ENDIF

    DO n = 1, SIZE (SwMir_Range)
       IF (angle(1) > SwMir_Range(n)%low_angle .AND. &
            angle(1) < SwMir_Range(n)%high_angle .AND. &
            angle(2) >  SwMir_Range(n)%low_angle .AND. &
            angle(2) <  SwMir_Range(n)%high_angle) THEN
          sw_pos = SwMir_Range(n)%pos
          EXIT
       ENDIF
    ENDDO

  END FUNCTION SwMirPos

!=============================================================================
  SUBROUTINE GetScAngles
!=============================================================================

    USE L0_sci_tbls, ONLY: APE2_dflt, TSE2_dflt
    USE MLSL1Config, ONLY: L1Config
    USE EngTbls, ONLY: EngMAF

    INTEGER :: i
    REAL :: APE, TSE
    REAL, PARAMETER :: APE_eps = 27.7340
    REAL, PARAMETER :: APE_B_A = 0.734
    REAL, PARAMETER :: TSE_eps = 26.301

    DO i = 0, (MaxMIFs - 1)
       IF (L1Config%Globals%SimOA) THEN
          APE = APE2_dflt(i)
          TSE = TSE2_dflt(i)
       ELSE
          APE = SciMAF(i)%APE_pos(2)
          IF (EngMAF%ASE_Side == "B") APE = APE - APE_B_A  ! Adjust for "B" side
          TSE = SciMAF(i)%TSSA_pos(2)
       ENDIF
       IF (APE >= 0.0) THEN
          SciMAF(i)%scAngleG = MOD ((APE + APE_eps), 360.0)
       ELSE
          SciMAF(i)%scAngleG = -999.9
       ENDIF
       IF (TSE >= 0.0) THEN
          SciMAF(i)%scAngleT = MOD ((TSE + TSE_eps), 360.0)
       ELSE
          SciMAF(i)%scAngleT = -999.9
       ENDIF
    ENDDO

  END SUBROUTINE GetScAngles

!=============================================================================
  SUBROUTINE Get_APE_ASA_pos (DN, APE_pos, ASA_pos)
!=============================================================================

    CHARACTER (LEN=24) :: DN
    REAL :: APE_pos(2), ASA_pos(2)

    INTEGER :: iape, iasa
    REAL, PARAMETER :: aaa_frac = 3706880.0

! Check for leading commands:

    SELECT CASE (ICHAR(DN(1:1)))

! possible leading commands (z'0e', z'0f', z'18', z'19', z'1a', z'1b', z'22'

       CASE (14, 15, 24, 25, 26, 27, 34)
          DN = DN(6:)   ! drop the first 5 bytes

    END SELECT

    IF (DN(1:1) /= CHAR(7) .AND. DN(1:1) /= CHAR(9)) DN = DN(4:)  ! shift

    iape = 0   ! indicate not yet available

    IF (DN(1:1) == CHAR(9)) THEN
       iape = 3
    ELSE IF (DN(10:10) == CHAR(9)) THEN
       iape = 12
    ENDIF

    IF (iape > 0) THEN
       APE_pos(1)  = deg24 * BigEndianStr (DN(iape:iape+2))
       APE_pos(2)  = deg24 * BigEndianStr (DN(iape+3:iape+5))
    ENDIF

    iasa = 0   ! indicate not yet available

    IF (DN(1:1) == CHAR(7)) THEN
       iasa = 3
    ELSE IF (DN(10:10) == CHAR(7)) THEN
       iasa = 12
    ENDIF

    IF (iasa > 0) THEN
       ASA_pos(1)  = BigEndianStr (DN(iasa:iasa+2)) / aaa_frac
       ASA_pos(2)  = BigEndianStr (DN(iasa+3:iasa+5)) / aaa_frac
    ENDIF

  END SUBROUTINE Get_APE_ASA_pos

!=============================================================================
  SUBROUTINE Get_PLL_DN (scipkt, sci_type, DN)
!=============================================================================

    CHARACTER (LEN=*), DIMENSION(:) :: scipkt
    CHARACTER (LEN=*) :: DN
    INTEGER :: sci_type

    IF (sci_type == 1) THEN
       DN(1:54) = scipkt(1)(601:654)
       DN(55:56) = scipkt(2)(941:942)
       DN(57:58) = scipkt(1)(123:124)       ! SM01
       DN(59:60) = scipkt(1)(225:226)       ! SM02
       DN(61:62) = scipkt(1)(499:500)       ! SM05
       DN(63:64) = scipkt(2)(95:96)         ! SM07
       DN(65:66) = scipkt(2)(169:170)       ! SM08
       DN(67:68) = scipkt(2)(315:316)       ! SM10
       DN(69:69)    = scipkt(1)(973:973)    ! DAC1
       DN(70:70)    = scipkt(2)(835:835)    ! DAC2
    ELSE
       DN(1:54) = scipkt(1)(527:580)
       DN(55:56) = scipkt(2)(956:957)
       DN(57:58) = scipkt(1)(123:124)       ! SM01
       DN(59:60) = scipkt(1)(225:226)       ! SM02
       DN(61:62) = scipkt(2)(95:96)         ! SM05
       DN(63:64) = scipkt(2)(169:170)       ! SM07
       DN(65:66) = scipkt(2)(243:244)       ! SM08
       DN(67:68) = scipkt(2)(389:390)       ! SM10
       DN(69:70) = CHAR(0)//CHAR(0)         ! clear DACs
       IF (sci_type == 3) THEN
          DN(69:69)    = scipkt(1)(983:983) ! DAC1
       ELSE
          DN(70:70)    = scipkt(1)(983:983) ! DAC2
       ENDIF
    ENDIF

  END SUBROUTINE Get_PLL_DN

!=============================================================================
  SUBROUTINE Band_switch (switch, sw_val, band)
!=============================================================================

    INTEGER :: band, sw_val, switch

    INTEGER :: i, pos

    INTEGER :: pos_no(6), FF
    DATA pos_no / z'FE', z'FD', z'FB', z'F7', z'EF', z'DF' /
    DATA FF / z'FF' /

    INTEGER, PARAMETER :: sw_bands(6,5) = RESHAPE ((/ &
         25, 26,  0,  0,  0,  0, &                 ! switch #1 bands
          2,  8,  3,  5,  6,  9, &                 ! switch #2 bands
          7, 13,  1,  8,  4, 21, &                 ! switch #3 bands
         10, 11, 14, 12, 20, 21, &                 ! switch #4 bands
         15, 16, 17, 18, 19, 20 /), &              ! switch #5 bands
          (/ 6, 5 /))                              ! final shape

    IF (sw_val == 0) THEN     ! no change from previous
       band = 0
       RETURN
    ENDIF

    IF (sw_val == FF) THEN  ! in process of changing
       band = -1
       PRINT *, 'switching bands for switch ', switch
       RETURN
    ENDIF

    DO i = 1, 6
       IF (sw_val == pos_no(i)) EXIT
    ENDDO
    pos = i
    IF (pos > 6) THEN
       band = 0
       RETURN
    ENDIF
    band = sw_bands (pos, switch)

  END SUBROUTINE Band_switch

!=============================================================================
  SUBROUTINE DetermineAttens (RIU, Addr, Val, BandSwitch, MaxAtten, AttenMaxed)
!=============================================================================

    USE L0_sci_tbls, ONLY: BandAtten
    USE MLSL1Common, ONLY: SwitchBank

    INTEGER, INTENT (IN) :: RIU, Addr, Val, BandSwitch(5)
    TYPE (BankLogical_T), INTENT (OUT) :: MaxAtten
    LOGICAL, INTENT (OUT) :: AttenMaxed

    INTEGER :: i, n, nMatch, nBanks, swBanks, MatchBand, BankIndx(2)
    LOGICAL :: BandMask(NumBands), SwitchMask(5), IsSwitchFB

    INTEGER, PARAMETER :: BandIndx(NumBands) = (/ (i,i=1,NumBands) /)
    INTEGER, PARAMETER :: DACS_indx(22:26) = (/ 4, 2, 3, 1, 1 /)
    INTEGER, PARAMETER :: MaxAttenVal = 63   ! Value for maximum attenuation

    AttenMaxed = (Val == MaxAttenVal)
    BandMask = .FALSE.  ! Nothing matches yet

    WHERE (BandAtten%RIU == RIU .AND. BandAtten%Addr == Addr)
       BandMask = .TRUE.          ! only for matched case(s)!
    ENDWHERE
    nMatch = COUNT (BandMask)

    DO i = 1, nMatch    ! 0, 1 or 2 possible

       BankIndx = 0     ! No indexes (yet)
       IF (i == 1) THEN
          MatchBand = MinVal (BandIndx, BandMask)
       ELSE
          MatchBand = MaxVal (BandIndx, BandMask)
       ENDIF

       SELECT CASE (MatchBand)

          CASE (1:21)     ! FBs

             IsSwitchFB = ANY (SwitchBank(2:5) == MatchBand) ! Is a switch FB
             SwitchMask = .FALSE.
             WHERE (BandSwitch == MatchBand)
                SwitchMask = .TRUE.
             END WHERE
             swBanks = COUNT (SwitchMask)
             IF (swBanks > 0) THEN
                BankIndx(1) = SwitchBank(MinVal (BandIndx(1:5), SwitchMask))
                BankIndx(2) = SwitchBank(MaxVal (BandIndx(1:5), SwitchMask))
                IF (IsSwitchFB) THEN
                   nBanks = swBanks
                ELSE
                   IF (MatchBand < 20) THEN
                      nBanks = 2
                      BankIndx(1) = MatchBand
                   ELSE
                      nBanks = swBanks
                   ENDIF
                ENDIF
             ELSE    ! Not in switch readbacks
                IF (.NOT. IsSwitchFB .AND. MatchBand < 20) THEN
                   nBanks = 1   ! Is just one bare FB
                   BankIndx(1) = MatchBand
                ELSE
                   nBanks = 0   ! FB not hooked up!
                ENDIF
             ENDIF

             DO n = 1, nBanks
                MaxAtten%FB(BankIndx(n)) = AttenMaxed
             ENDDO

          CASE (22:26)    ! DACS

             BankIndx = DACS_indx(MatchBand)
             IF (MatchBand < 25) THEN
                nBanks = 1
             ELSE
                IF (MatchBand == BandSwitch(1)) THEN
                   nBanks = 1
                ELSE
                   nBanks = 0
                ENDIF
             ENDIF

             DO n = 1, nBanks
                MaxAtten%DACS(BankIndx(n)) = AttenMaxed
             ENDDO

          CASE (27:31)    ! MBs

             nBanks = 1
             BankIndx = MatchBand - 26

             DO n = 1, nBanks
                MaxAtten%MB(BankIndx(n)) = AttenMaxed
             ENDDO

          CASE (32:34)    ! WFs

             nBanks = 1
             BankIndx = MatchBand - 31

             DO n = 1, nBanks
                MaxAtten%WF(BankIndx(n)) = AttenMaxed
             ENDDO

       END SELECT

    ENDDO

  END SUBROUTINE DetermineAttens

!=============================================================================
  SUBROUTINE Save_THz_pkt (Sci_pkt, THz_Sci_pkt)
!=============================================================================

    USE L0_sci_tbls, ONLY: Sci_pkt_T, THz_Sci_pkt_T
    USE THzUtils, ONLY: LLO_Bias
    USE MLSL1Common, ONLY: Deflt_zero

    TYPE (Sci_pkt_T), INTENT (IN) :: Sci_pkt
    TYPE (THz_Sci_pkt_T), INTENT (OUT) :: THz_Sci_pkt

    INTEGER :: bank

    THz_Sci_pkt%secTAI = Sci_pkt%secTAI
    THz_Sci_pkt%MAFno = Sci_pkt%MAFno
    THz_Sci_pkt%MIFno = Sci_pkt%MIFno
    THz_Sci_pkt%Orbit = Sci_pkt%Orbit
    THz_Sci_pkt%FB(:,1:5) = Sci_pkt%FB(:,15:)
    DO bank = 1, 5
       IF (ANY (THz_Sci_pkt%FB(:,bank) > 0)) THz_Sci_pkt%FB(:,bank) = &
            THz_Sci_pkt%FB(:,bank) - Deflt_zero%FB(:,bank+14)
    ENDDO
    THz_Sci_pkt%FB(:,6) = Sci_pkt%FB(:,12)
    IF (ANY (THz_Sci_pkt%FB(:,6) > 0)) THz_Sci_pkt%FB(:,6) = &
         THz_Sci_pkt%FB(:,6) - Deflt_zero%FB(:,12)
    THz_Sci_pkt%MaxAtten(1:5) = Sci_pkt%MaxAtten%FB(15:19)
    THz_Sci_pkt%MaxAtten(6) = Sci_pkt%MaxAtten%FB(12)
    THz_Sci_pkt%AttenMaxed = (ANY (THz_Sci_pkt%MaxAtten))
    THz_Sci_pkt%TSSA_pos = Sci_pkt%TSSA_pos
    THz_Sci_pkt%SwMirPos = Sci_pkt%THz_sw_pos
    THz_Sci_pkt%LLO_bias = LLO_Bias (Sci_pkt%LLO_DN, Sci_pkt%MIFno)
    THz_Sci_pkt%BandSwitch = Sci_pkt%BandSwitch(4:5)  !Only need sw #4 & #5
    THz_Sci_pkt%CRC_good = Sci_pkt%CRC_good

  END SUBROUTINE Save_THz_pkt

END MODULE SciUtils

! $Log$
! Revision 2.7  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.6  2003/09/15 17:15:54  perun
! Version 1.3 commit
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.3  2002/10/03 17:43:53  jdone
! moved parameter statement to data statement for LF/NAG compatitibility
!
! $Log$
! Revision 2.7  2004/01/09 17:46:23  perun
! Version 1.4 commit
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:56:11  perun
! Version 0.5 commit
!
