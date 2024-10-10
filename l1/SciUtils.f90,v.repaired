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
MODULE SciUtils ! L0 science utilities
!=============================================================================

  USE L0_sci_tbls, ONLY: sci_type, l0_sci1, l0_sci2, sci1_T1_fmt, sci2_T1_fmt, &
       sci1_T2_fmt, sci2_T2_fmt, sci1_T1, sci2_T1, sci1_T2, sci2_T2, sci_cptr, &
       Sci_pkt, SciMAF, type_sci1, THzSciMAF, DACS_pkt, DACS_MAF, Atten_T
  USE L0Utils, ONLY: ReadL0Sci, CheckSum
  USE MLSL1Utils, ONLY: BigEndianStr, ExtractBigEndians, SwapBytes, QNan, &
       Finite
  USE DACsUtils, ONLY: ExtractDACSdata, UncompressDACSdata, ProcessDACSdata
  USE MLSL1Common, ONLY: FBnum, MBnum, WFnum, WFchans, deg24, BankLogical_T, &
       DACSchans, BandSwitch, L1ProgType, THzType, MaxMIFs, BankInt_T

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: NextSciMAF, SwMirPos, GetScAngles

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

  TYPE (BankLogical_T) :: MaxAtten = & ! initialize to NOT Max Atten
       BankLogical_T (.FALSE., .FALSE., .FALSE., .FALSE.)
  TYPE (BankLogical_T) :: DeltaAtten = & ! initialize to NOT Delta Atten
       BankLogical_T (.FALSE., .FALSE., .FALSE., .FALSE.)
  TYPE (BankInt_T) :: BankAtten = & ! initialize to NO bank Atten reading
       BankInt_T (0, 0, 0, 0)

  CONTAINS

!=============================================================================
  FUNCTION GetSciPkt () RESULT (OK)
!=============================================================================

    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info
    USE MLSL1Common, ONLY: L1BFileInfo
    USE ERMSG_M, ONLY: ermset
    USE THzUtils, ONLY: ConvertLLO
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

    LOGICAL OK

    INTEGER :: i, ios, j, n, returnStatus
    INTEGER :: tindex   ! Type index: 1 = Type I, 2 = Type II/III

    !! Pointers to counts:

    INTEGER, DIMENSION(:), POINTER :: FBcnts => NULL()
    INTEGER, DIMENSION(:), POINTER :: MBcnts => NULL()
    INTEGER, DIMENSION(:), POINTER :: WFcnts => NULL()
    CHARACTER(len=1), PARAMETER :: sw_good = CHAR(7)
    CHARACTER(len=2) :: WFdat(WFchans)   ! wide filter raw data buffer
    CHARACTER (LEN=1024) :: scipkt(2)
    CHARACTER (LEN=24) :: DN
    CHARACTER (LEN=80) :: msg
    CHARACTER(len=27) :: asciiUTC

    TYPE (Atten_T) :: Attenuation
    LOGICAL :: AttenMaxed

    INTEGER :: DACS_C_K(DACSchans)
    INTEGER :: DACS_i1, DACS_i2
    INTEGER :: D(4), TP, DIO, LO, Zlag, IDN(512)

    INTEGER, PARAMETER :: type_I = 1
    INTEGER, PARAMETER :: type_II = 3
    INTEGER, PARAMETER :: type_III = 5
    logical :: show
    integer, save :: mifTotal = -1

    INTEGER, EXTERNAL :: PGS_TD_EOSPMGIRDtoTAI

    CALL ReadL0Sci (scipkt, OK)
    mifTotal = mifTotal + 1

    l0_sci1 = scipkt(1)
    l0_sci2 = scipkt(2)

    sci_type = ICHAR(type_sci1)  ! science packet #1 type

    ! attempt to convert the packet, even though it may not be 'OK'. We may be
    ! able to get more info about the packet
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
       !OK = .FALSE.
       !RETURN
       READ (UNIT=l0_sci1, FMT=sci1_T1_fmt, iostat=ios) sci1_T1
       READ (UNIT=l0_sci2, FMT=sci2_T1_fmt, iostat=ios) sci2_T1
       tindex = 1
       CALL MLSMessage (MLSMSG_Info, ModuleName, "Incorrect sci type!")

    ENDIF

    !! Convert raw data:

    Sci_pkt%MAFno = BigEndianStr (sci_cptr(tindex)%MAF(1)%ptr)
    Sci_pkt%MIFno = BigEndianStr (sci_cptr(tindex)%MIF(1)%ptr)
    Sci_pkt%Orbit = BigEndianStr (sci_cptr(tindex)%orbit(1)%ptr)

    IF (.NOT. OK) THEN 
       print *, 'Science packet was not OK'
       print *, "but check messages, could be we're out of the processing window!"
       RETURN
    ENDIF



! Get TAI93 time

    returnstatus = PGS_TD_EOSPMGIRDtoTAI (scipkt(1)(8:15), Sci_pkt%secTAI)
    Sci_pkt%secTAI = Sci_pkt%secTAI - 0.25  ! Adjust to actual time taken (whd:
    !Why?)  answer: 0.25 seconds is about 1.5 MIFs. Dominick thinks the time
    !reported is actually this amount of time away from the actual start of the
    !MIF

!! Check for bad checksums:

    Sci_pkt%CRC_good = .TRUE.
    DO i = 1, 2
       DO j = 1, 512
          IDN(j) = BigEndianStr (scipkt(i)((j*2-1):j*2))
       ENDDO
       IF (CheckSum (IDN, 511) /= IDN(512)) THEN
          Sci_pkt%CRC_good = .FALSE.
          n = PGS_TD_TAItoUTC (Sci_pkt%secTAI, asciiUTC)
          WRITE (msg, &
               '("Bad Checksum: Sci Pkt ", i1, ", MIF: ", i3, ", UTC: ", A27)')&
               i, Sci_pkt%MIFno, asciiUTC
          PRINT *, TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) '### Info: '//TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
          CALL MLSMessage (MLSMSG_Info, ModuleName, TRIM(msg))
          RETURN   ! Can't do any more
       ENDIF
    ENDDO

    !! Convert various positions counts to angles:


    !! THz_sw is the THz Single Scanning Mirror position information (THz
    !! doesn't have a switching mirrog like the GHz module

    DN = sci_cptr(tindex)%THz_sw
    IF (DN(1:1) /= sw_good) DN = DN(4:)   ! shift and try
    IF (DN(1:1) == sw_good) THEN
       Sci_pkt%TSSM_pos(1)  = deg24 * BigEndianStr (DN(3:5))
       Sci_pkt%TSSM_pos(2)  = deg24 * BigEndianStr (DN(6:8))
    ELSE
       Sci_pkt%TSSM_pos(:) = QNan()
    ENDIF

    !    Sci_pkt%THz_sw_pos = SwMirPos ("T", Sci_pkt%TSSM_pos)
    ! GHz_sw contains the GHz Switching mirror positions
    DN = sci_cptr(tindex)%GHz_sw
    IF (DN(1:1) == sw_good) THEN
       Sci_pkt%GSM_pos(1)  = deg24 * BigEndianStr (DN(3:5))
       Sci_pkt%GSM_pos(2)  = deg24 * BigEndianStr (DN(6:8))
    ELSE
       Sci_pkt%GSM_pos(:) = QNan()
    ENDIF

    !    Sci_pkt%GHz_sw_pos = SwMirPos ("G", Sci_pkt%GSM_pos)

    !! Get GHz scanning angles

    ! This is the GHz mirror encoder, giving the position of the Limb scanning
    ! mirror

    Sci_pkt%APE_pos(:) = QNan()
    Sci_pkt%ASA_pos(:) = QNan()
    DN = sci_cptr(tindex)%GHz_ant_scan
    if ( DN(1:1) == achar(0) ) then
       n = PGS_TD_TAItoUTC (Sci_pkt%secTAI, asciiUTC)
      print *, 'About to fail during Get_APE_ASA_pos at MAF ', Sci_pkt%MAFno
      print *, 'Packet UTC: ',asciiUTC
    endif
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
          WRITE (msg, '("S", i1, " switching to Band ", i2)') i, bandswitch(i)
          PRINT *, TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) '### Info: '//TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
          CALL MLSMessage (MLSMSG_Info, ModuleName, TRIM(msg))
       ENDIF
    ENDDO

!! (laser Local oscillator) LLO data (sometimes you see)GLO

    Sci_pkt%LLO_DN = sci_cptr(tindex)%LLO_DN

    CALL ConvertLLO (Sci_pkt%LLO_DN, Sci_pkt%LLO_EU)

!! (Phase Lock Loop Oscillator) PLL data
    show = (Sci_pkt%MIFNo == 0)
    ! Because the day's first Sci_pkt%MIFNo may be offset from 0
    if ( mifTotal == 0 ) mifTotal = Sci_pkt%MIFNo

    CALL Get_PLL_DN (scipkt, sci_type, Sci_pkt%PLL_DN, show=show, mifTotal=mifTotal)

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
         Sci_pkt%BandSwitch, MaxAtten, AttenMaxed, DeltaAtten)

    Sci_pkt%MaxAtten = MaxAtten   ! Use latest attenuation flags
    Sci_pkt%AttenMaxed = AttenMaxed
    Sci_pkt%DeltaAtten = DeltaAtten

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
    USE MLSL1Common, ONLY : L1ProgType, LogType, MAFinfo
    USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

    !! Get the next MAF's science data

    LOGICAL, INTENT (OUT) :: more_data

    INTEGER, PARAMETER :: no_data = -1   ! no data is available
    INTEGER, SAVE :: prev_MAF = no_data, prev_last_MIF = no_data
    INTEGER :: i, last_MIF, MIFno, m
    REAL :: pos1, pos2, dif
    REAL :: APE_pos(0:(MaxMIFs-1),2), ASA_pos(0:(MaxMIFs-1),2)
    REAL :: GSM_pos(0:(MaxMIFs-1),2), TSSM_pos(0:(MaxMIFs-1),2)
    REAL :: APE_theta(0:(MaxMIFs-1)), GSM_theta(0:(MaxMIFs-1)), &
         TSSM_theta(0:(MaxMIFs-1))
    REAL :: APE_pos_P(0:(MaxMIFs-1),2), TSSM_pos_P(0:(MaxMIFs-1),2)

    REAL, PARAMETER :: pos_offset = 4.44e-03   ! Offset between MIF 0, Pos 2
                                               ! and MIF 1, Pos 1
    CHARACTER (LEN=80) :: msg
    CHARACTER(len=27) :: asciiUTC

    more_data = .TRUE.

    !! Initialize MAF/MIF counters to indicate no data

    SciMAF%MAFno = no_data
    SciMAF%MIFno = no_data
    THzSciMAF%MAFno = no_data
    THzSciMAF%MIFno = no_data
    last_MIF = L1Config%Calib%MIFsPerMAF - 1
    DO i = 0, (MaxMIFs - 1)
     DACS_MAF(i)%D = 0
    ENDDO
    APE_pos = QNan(); ASA_pos = QNan(); GSM_pos = QNan(); TSSM_pos = QNan()
    SciMAF%secTAI = QNan()

    !! Initialize CRC flags to good

    SciMAF%CRC_good = .TRUE.

    !! Save previously read packet (if available):

    IF (prev_MAF /= no_data) THEN
       MIFno = Sci_pkt%MIFno
       SciMAF(MIFno) = Sci_pkt
       DACS_MAF(MIFno) = DACS_pkt
       IF (L1ProgType == THzType) THEN
          CALL Save_THz_pkt (SciMAF(MIFno), THzSciMAF(MIFno))
       ENDIF
       APE_pos(MIFno,:) = SciMAF(MIFno)%APE_pos
       ASA_pos(MIFno,:) = SciMAF(MIFno)%ASA_pos
       GSM_pos(MIFno,:) = SciMAF(MIFno)%GSM_pos
       TSSM_pos(MIFno,:) = SciMAF(MIFno)%TSSM_pos

! Adjust time in MIF 0, if not available
!
! whd: That comment doesn't make sense! This code runs only when MIF != 0!

       IF (MIFno /= 0) THEN
          SciMAF(0)%secTAI = SciMAF(MIFno)%secTAI - MIFno * &
               L1Config%Calib%MIF_duration
          SciMAF(0)%MAFno = SciMAF(MIFno)%MAFno
       ENDIF

    ENDIF

    DO

       IF (GetSciPkt()) THEN

          MIFno = Sci_pkt%MIFno
          ! print *, 'mif :', MIFno
          IF ((prev_MAF /= no_data) .AND. (Sci_pkt%MAFno /= prev_MAF)) THEN
             prev_MAF = Sci_pkt%MAFno
             EXIT                      ! Previous MAF is done
          ENDIF

          SciMAF(MIFno) = Sci_pkt  ! save current packet
          DACS_MAF(MIFno) = DACS_pkt  ! save DACS data

          IF (L1ProgType == THzType) &
           CALL Save_THz_pkt (SciMAF(MIFno), THzSciMAF(MIFno))

! Save mechanism positions for further processing:

          APE_pos(MIFno,:) = SciMAF(MIFno)%APE_pos
          ASA_pos(MIFno,:) = SciMAF(MIFno)%ASA_pos
          GSM_pos(MIFno,:) = SciMAF(MIFno)%GSM_pos
          TSSM_pos(MIFno,:) = SciMAF(MIFno)%TSSM_pos

          prev_MAF = Sci_pkt%MAFno
          prev_last_MIF = MIFno

       ELSE

          more_data = .FALSE.
          prev_MAF = no_data
          EXIT  ! Nothing more is available

       ENDIF

    ENDDO

    ! Seems backwards to me, what if SciMAF(m+1) == NaN? 
    DO m = 0, 2    ! Make sure time is available in first 3 MIFs!
       IF (.NOT. Finite (SciMAF(m)%secTAI)) SciMAF(m)%secTAI = &
            SciMAF(m+1)%secTAI - L1Config%Calib%MIF_duration
    ENDDO

    ! Shouldn't it be
    ! DO m = 2,0,-1    ! Make sure time is available in first 3 MIFs!
    !    IF (.NOT. Finite (SciMAF(m)%secTAI)) SciMAF(m)%secTAI = &
    !         SciMAF(m+1)%secTAI - L1Config%Calib%MIF_duration
    ! ENDDO


! Determine pointing mechanisms angles:

    APE_theta = QNan()
    GSM_theta = QNan()
    TSSM_theta = QNan()
    APE_pos_P = QNan()
    TSSM_pos_P = QNan()

    !do m = 3, (MaxMIFS - 1), 3
       !APE_pos(m,:) = QNan()    ! for testing!!!!
    !enddo
    DO m = 0, (MaxMIFs - 2)

       IF (m < (MaxMIFs - 3)) THEN   ! APE pos
          pos1 = APE_pos(m+2,1)
          pos2 = APE_pos(m+1,2)
          IF (.NOT. Finite (pos1)) THEN
             APE_theta(m) = pos2
             IF (m == 0 .AND. (.NOT. Finite (pos2))) THEN  ! No good pos data!
                APE_theta(m) = APE_pos(0,2)        ! Use good MIF 0, pos 2
                IF (prev_last_MIF /= 148) THEN     ! Nominal previous MAF number
                   APE_theta(m) = APE_theta(m) + pos_offset
                ENDIF
                IF (Finite (APE_theta(m))) THEN
                   i = PGS_TD_TAItoUTC (SciMAF(0)%secTAI, asciiUTC)
                   WRITE (msg, &
                    '("Adjusting MIF 0 pointing using offset at UTC: ", A27)') &
                    asciiUTC
                   CALL MLSMessage (MLSMSG_Info, ModuleName, msg)
                ENDIF
             ENDIF
          ELSE IF (.NOT. Finite (pos2)) THEN
             APE_theta(m) = pos1
          ELSE
             dif = MOD ((pos1 + APE_pos(m+1,1) - 2.0 * pos2 + 900.0), 360.0)
             APE_theta(m) = pos2 + 0.25 * (dif - 180.0)
          ENDIF
       ENDIF
       
       ! THz Single Scanning Mirror (?)

       pos1 = TSSM_pos(m,1)
       pos2 = TSSM_pos(m+1,2)
       IF (.NOT. Finite (pos1)) THEN
          TSSM_theta(m) = pos2
       ELSE IF (.NOT. Finite (pos2)) THEN
          TSSM_theta(m) = pos1
       ELSE
          dif = MOD ((pos1 + TSSM_pos(m+1,1) - 2.0 * pos2 + 900.0), 360.0)
          TSSM_theta(m) = pos2 + 0.25 * (dif - 180.0)
       ENDIF

       ! GHz switching Mirror
       pos1 = GSM_pos(m,1)
       pos2 = GSM_pos(m+1,2)

       !<whd> why GSM_pos(m,1) and GSM_pos(m+1,2)</!whd>

       IF (.NOT. Finite (pos1)) THEN
          GSM_theta(m) = pos2
       ELSE IF (.NOT. Finite (pos2)) THEN
          GSM_theta(m) = pos1
       ELSE
          dif = MOD ((pos1 + GSM_pos(m+1,1) - 2.0 * pos2 + 900.0), 360.0)
          GSM_theta(m) = pos2 + 0.25 * (dif - 180.0)
       ENDIF

       ! Save Pos1 and Pos2 Prime to save in the L1BOA file:

       TSSM_pos_P(m,1) = TSSM_pos(m,1)
       IF (m == 1 .OR. m > 2) THEN
          APE_pos_P((m-1),2) = APE_pos(m,2)
          TSSM_pos_P((m-1),2) = TSSM_pos(m,2)
          APE_pos_P((m-1),1) = APE_pos(m,1)
       ELSE IF (m == 2) THEN
          APE_pos_P(1,1) = 2 * APE_pos(1,2) - APE_pos(1,1) 
          TSSM_pos_P(2,1) = 2 * TSSM_pos(3,2) - TSSM_pos(3,1)
          APE_pos_P(1,2) = 2 * APE_pos(3,1) - APE_pos(3,2)
          TSSM_pos_P(1,2) = 2 * TSSM_pos(1,1) - TSSM_pos(1,2)
       ENDIF

    ENDDO
    APE_theta(MaxMIFs-2) = APE_pos(MaxMIFs-1,2)
    APE_theta(MaxMIFs-1) = APE_pos(MaxMIFs-1,1)
    TSSM_theta(MaxMIFs-1) = TSSM_pos(MaxMIFs-1,1)
    GSM_theta(MaxMIFs-1) = GSM_pos(MaxMIFs-1,1)

    DO m = 0, (MaxMIFs - 1)
       SciMAF(m)%APE_theta = APE_theta(m)
       SciMAF(m)%GSM_theta = GSM_theta(m)
       SciMAF(m)%GHz_sw_pos = SwMirPos ("G", GSM_theta(m))
       SciMAF(m)%THz_sw_pos = SwMirPos ("T", TSSM_theta(m))
       SciMAF(m)%TSSM_theta = TSSM_theta(m)
       SciMAF(m)%APE_pos_P = APE_pos_P(m,:)
       SciMAF(m)%TSSM_pos_P = TSSM_pos_P(m,:)
       THzSciMAF(m)%TSSM_theta = TSSM_theta(m)
       THzSciMAF(m)%SwMirPos = SwMirPos ("T", TSSM_theta(m))
    ENDDO

! Process DACS data for the current MAF

    IF (L1ProgType == LogType) CALL ProcessDACSdata

! Get scAngles to be used for L1BOA

    CALL GetScAngles

    IF (L1ProgType == THzType) THEN
       DO m = 0, (MaxMIFs - 1)
          THzSciMAF(m)%scAngle = SciMAF(m)%scAngleT
       ENDDO
    ENDIF

  END SUBROUTINE NextSciMAF

!=============================================================================
  FUNCTION SwMirPos (sw_module, angle) RESULT (sw_pos)
!=============================================================================

    ! returns a 1 letter code telling whether this mirror position is to be
    ! 'D'iscarded or not. If not, the possible returns are 'S' (space target)',
    ! 'T' (calibration target)', 'L' (limb view)

    USE MLSL1Common, ONLY: GHz_SwMir_Range_A, GHz_SwMir_Range_B, r8, &
     GHz_SwMir_Range_B_2, THz_SwMir_Range,  Discard_SwMir_Range, SwMir_Range_T
    USE EngTbls, ONLY: EngMAF

    CHARACTER(len=1), INTENT(IN) :: sw_module
    REAL, INTENT(IN) :: angle
    CHARACTER(len=1) :: sw_pos

    INTEGER :: n
    TYPE (SwMir_Range_T), DIMENSION(:), POINTER :: SwMir_Range

    ! 2011 DOY 153 19:10:00 adjust time. 

    ! <whd> 
    ! 
    ! I think this `TAI' time, is really TAI93. It works out to
    ! 2011-153T19:07:53. Dom tells me that we went to nominal on the new B side
    ! table at 2011-153T19:18, so maybe this time is just to set things up to
    ! make sure we're ready. I'm not going to change it, I'll just use it as is
    ! here and in my IDL code.
    !
    ! </whd>

    REAL(r8), PARAMETER :: B_TAI = 5.8112628e+08+69000d00

    sw_pos = "D"                 ! Initialize to "Discard"

    IF (.NOT. Finite (angle)) RETURN

    
    IF (sw_module == "G") THEN   

       ! GigaHertz Module Determine which ranges to use. The ranges have changed
       ! over the course of the mission. See MLSL1Common (~line 215) for
       ! definitions

       IF (EngMAF%GSM_Side == "A") THEN
          SwMir_Range => GHz_SwMir_Range_A
       ELSE IF (EngMAF%GSM_Side == "B") THEN
          IF (SciMAF(0)%secTAI > B_TAI) THEN   ! use MIF 0 TAI time
             SwMir_Range => GHz_SwMir_Range_B_2
          ELSE
             SwMir_Range => GHz_SwMir_Range_B
          ENDIF
       ELSE
          SwMir_Range => Discard_SwMir_Range
       ENDIF
    ELSE                         ! TeraHertz Module
       SwMir_Range => THz_SwMir_Range
    ENDIF

    DO n = 1, SIZE (SwMir_Range)
       IF (angle > SwMir_Range(n)%low_angle .AND. &
            angle < SwMir_Range(n)%high_angle) THEN
          sw_pos = SwMir_Range(n)%pos
          EXIT
       ENDIF
    ENDDO

  END FUNCTION SwMirPos

!=============================================================================
  SUBROUTINE GetScAngles
!=============================================================================

    USE L0_sci_tbls, ONLY: APE_theta_dflt, TSSM_theta_dflt
    USE MLSL1Config, ONLY: L1Config
    USE EngTbls, ONLY: EngMAF
    use MLSFillValues, only: Isnan

    INTEGER :: i
    REAL :: APE, TSSM
    REAL, PARAMETER :: APE_eps = 27.7340
    REAL, PARAMETER :: APE_B_A = 0.734
    REAL, PARAMETER :: TSSM_eps = 26.232   !26.301 before 9/20/04

!     print *, 'In GetScAngles'
!     print *, 'shape(APE_theta_dflt) ', shape(APE_theta_dflt)
!     print *, 'shape(TSSM_theta_dflt) ', shape(TSSM_theta_dflt)
!     print *, 'shape(SciMAF) ', shape(SciMAF)
!     print *, 'lbound,ubound(APE_theta_dflt) ', lbound(APE_theta_dflt), ubound(APE_theta_dflt)
!     print *, 'lbound,ubound(SciMAF) ', lbound(SciMAF), ubound(SciMAF)
!     do i = 0, (MaxMIFs - 1)
!           print *, 'i, SciMAF(i)%APE_theta ', i, SciMAF(i)%APE_theta
!           print *, 'i, SciMAF(i)%TSSM_theta ', i, SciMAF(i)%TSSM_theta
!     enddo
!     stop
    DO i = 0, (MaxMIFs - 1)
       IF (L1Config%Globals%SimOA) THEN
          APE = APE_theta_dflt(i)
          TSSM = TSSM_theta_dflt(i)
!           print *, 'SimOA ', i, APE
       ELSE
          APE = SciMAF(i)%APE_theta
          IF (EngMAF%ASE_Side == "B") APE = APE - APE_B_A  ! Adjust for "B" side
          TSSM = SciMAF(i)%TSSM_theta
!           print *, 'not SimOA ', i, APE
       ENDIF
       if ( isNaN(APE) ) then
         APE = -999.9
         SciMAF(i)%scAngleG = -999.9
       endif
       if ( isNaN(TSSM) ) then
         TSSM = -999.9
         SciMAF(i)%scAngleT = -999.9
       endif
       IF (APE >= 0.0) THEN
          SciMAF(i)%scAngleG = MOD ((APE + APE_eps), 360.0)
       ELSE
          SciMAF(i)%scAngleG = -999.9
       ENDIF
       IF (TSSM >= 0.0) THEN
          SciMAF(i)%scAngleT = MOD ((TSSM + TSSM_eps), 360.0)
       ELSE
          SciMAF(i)%scAngleT = -999.9
       ENDIF
    ENDDO

  END SUBROUTINE GetScAngles

!=============================================================================
  SUBROUTINE Get_APE_ASA_pos (DN, APE_pos, ASA_pos)
!=============================================================================

    ! APE = Antenna Position Electronics 
    !
    ! ASA = Antenna Switch A<something>. Might also go by the name ASE (Antenna
    ! Scan Electronics), see type Eng_MAF_T defined in EngTbls.f90


    CHARACTER (LEN=24) :: DN
    REAL :: APE_pos(2), ASA_pos(2)

    INTEGER :: iape, iasa
    REAL, PARAMETER :: aaa_frac = 3706880.0
    integer, parameter :: NULLIS_ZERO = 0

! Check if DN is all or partly NULLs
    if ( ichar(DN(1:1)) == NULLIS_ZERO ) then
      print *, 'Get_APE_ASA_pos failed 1st char in DN is NULL'
      return
    endif

! Check for leading commands:

    SELECT CASE (ICHAR(DN(1:1)))

! possible leading commands (z'0e', z'0f', z'18', z'19', z'1a', z'1b', z'22') ==
! 14, 15, 24, 25, 26, 27, 34 (so what is a 'leading command?')

       CASE (14, 15, 24, 25, 26, 27, 34)
          DN = DN(6:)   ! drop the first 5 bytes

    END SELECT

    ! The codes that indicate that this is the GHz antenna position counts is 7
    ! or 9. If either, the data begins at postion 4
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
    
    if ( iape < 1 .and. iasa < 1 ) THEN 
       print *, 'Get_APE_ASA_pos failed with DN:', DN
    ENDIF

  END SUBROUTINE Get_APE_ASA_pos

!=============================================================================
  SUBROUTINE Get_PLL_DN (scipkt, sci_type, DN, show, mifTotal)
!=============================================================================

    CHARACTER (LEN=*), DIMENSION(:) :: scipkt
    CHARACTER (LEN=*) :: DN
    INTEGER :: sci_type
    logical, optional, intent(in)   :: show
    integer, optional, intent(in)   :: mifTotal
    character, dimension(2)         :: GM07
    logical                         :: myShow
    integer, parameter              :: BadFirstByte = 13
    integer, parameter              :: BadSecondByte = 202
    integer, parameter              :: GoodFirstByte = 12
    integer, parameter              :: GoodSecondByte = 194

    ! An example of how not to write software
    ! Comes with no comments, explanation, or apologies
    ! Sigh.
    ! Here's what we leaarned from an idl snippet
     !   help, /structure, pll
     !   ** Structure <1e31ce8>, 20 tags, length=70, data length=70, refs=1:
     !      GM05            BYTE        14
     !      GM06            BYTE        14
     !      GM07            BYTE      Array[2]
     !      GM08            BYTE      Array[2]
     !      GM09            BYTE      Array[2]
     !      GM10            BYTE      Array[13]
     !      GM11            BYTE      Array[13]
     !      GM12            BYTE       255
     !      GM13            BYTE      Array[13]
     !      GM14            BYTE       243
     !      GM15            BYTE      Array[5]
     !      TM03            BYTE      Array[2]
     !      SM01            BYTE      Array[2]
     !      SM02            BYTE      Array[2]
     !      SM05            BYTE      Array[2]
     !      SM07            BYTE      Array[2]
     !      SM08            BYTE      Array[2]
     !      SM10            BYTE      Array[2]
     !      DAC1            BYTE        13
     !      DAC2            BYTE        13
    ! Do they at least line up?
    ! The offset for SM01 from the above is
    ! 2 + 3*2 + 2*13 + 1 + 13 + 1 + 5 + 2 = 56
    ! So that's why SM01 below begins at 57. Mann, that's a pain.
    ! OK so if we want to look at GM07, we'll want
    !  DN(3:4) = scipkt(1)(603:604)         ! GM07
    myshow = .false.
    if ( present(show) ) myShow = show
    
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
    ! if ( .not. myShow ) return
    GM07(1) = DN(3:3)
    GM07(2) = DN(4:4)
    ! print *, 'sci_type: ', sci_type, ' ', DN
    if ( .false. ) then
      print *, mifTotal, ' GM07: ', ichar(GM07)
    ! elseif ( ichar(GM07(1)) == BadFirstByte ) then
    elseif ( ichar(GM07(1)) /= GoodFirstByte ) then
      print *, mifTotal, ' Bad first byte GM07: ', ichar(GM07(1))
    ! elseif ( ichar(GM07(2)) == BadSecondByte ) then
    elseif ( ichar(GM07(2)) /= GoodSecondByte ) then
      print *, mifTotal, ' Bad second byte GM07: ', ichar(GM07(2))
    endif

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
  SUBROUTINE DetermineAttens (RIU, Addr, Val, BandSwitch, MaxAtten, &
       AttenMaxed, DeltaAtten)
!=============================================================================

    USE L0_sci_tbls, ONLY: BandAtten
    USE MLSL1Common, ONLY: SwitchBank

    INTEGER, INTENT (IN) :: RIU, Addr, Val, BandSwitch(5)
    TYPE (BankLogical_T), INTENT (OUT) :: MaxAtten, DeltaAtten
    LOGICAL, INTENT (OUT) :: AttenMaxed

    INTEGER :: i, n, nMatch, nBanks, swBanks, MatchBand(2), BankIndx(2), tblIndx
    LOGICAL :: BandMask(SIZE(BandAtten)), SwitchMask(5), IsSwitchFB

    INTEGER, PARAMETER :: BandIndx(SIZE(BandAtten)) = &
         (/ (i, i=1, SIZE(BandAtten)) /)
    INTEGER, PARAMETER :: DACS_indx(22:26) = (/ 4, 2, 3, 1, 1 /)
    INTEGER, PARAMETER :: MaxAttenVal = 63   ! Value for maximum attenuation

    DeltaAtten = BankLogical_T (.FALSE., .FALSE., .FALSE., .FALSE.)
    AttenMaxed = (Val == MaxAttenVal)
    BandMask = .FALSE.  ! Nothing matches yet


    !BandAtten loaded by table in L0_sci_tbls:~382. If BandAtten(tblIndx)%Link
    !!= 0 then there are two banks (channels?, ??) to look at

    !find the indices of records in BandAtten that match the input RIU and
    !Addr.
    WHERE (BandAtten%RIU == RIU .AND. BandAtten%Addr == Addr)
       BandMask = .TRUE.          ! only for matched case!
    ENDWHERE
    IF (COUNT (BandMask) /= 1) RETURN ! nothing to do, if none match

    ! Now find that records
    tblIndx = MAXVAL (BandIndx, BandMask)

    ! extract the Band and Link fields
    MatchBand = (/ BandAtten(tblIndx)%Band, BandAtten(tblIndx)%Link /)
    ! if %Link != 0 we need to check 2 bands (some are 'linked')
    IF (MatchBand(2) /= 0) THEN
       nMatch = 2
    ELSE
       nMatch = 1
    ENDIF

    DO i = 1, nMatch    ! 1 or 2 possible

       BankIndx = 0     ! No indexes (yet)

       SELECT CASE (MatchBand(i))

          CASE (1:21)     ! FBs

             ! If the input Band goes through the GHz switch
             IsSwitchFB = ANY (SwitchBank(2:5) == MatchBand(i)) ! Is a switch FB
             SwitchMask = .FALSE.
             WHERE (BandSwitch == MatchBand(i))
                SwitchMask = .TRUE. ! == T where the input Band is going thru GHz switch?
             END WHERE
             swBanks = COUNT (SwitchMask)
             IF (swBanks > 0) THEN
                BankIndx(1) = SwitchBank(MINVAL (BandIndx(1:5), SwitchMask))
                BankIndx(2) = SwitchBank(MAXVAL (BandIndx(1:5), SwitchMask))
                IF (IsSwitchFB) THEN
                   nBanks = swBanks
                ELSE
                   IF (MatchBand(i) < 20) THEN
                      nBanks = 2
                      BankIndx(1) = MatchBand(i)
                   ELSE
                      nBanks = swBanks
                   ENDIF
                ENDIF
             ELSE    ! Not in switch readbacks
                IF (.NOT. IsSwitchFB .AND. MatchBand(i) < 20) THEN
                   nBanks = 1   ! Is just one bare FB
                   BankIndx(1) = MatchBand(i)
                ELSE
                   nBanks = 0   ! FB not hooked up!
                ENDIF
             ENDIF

             DO n = 1, nBanks
                MaxAtten%FB(BankIndx(n)) = AttenMaxed
                DeltaAtten%FB(BankIndx(n)) = &
                     (BankAtten%FB(BankIndx(n)) /= Val)
                BankAtten%FB(BankIndx(n)) = Val
             ENDDO

          CASE (22:26)    ! DACS

             BankIndx = DACS_indx(MatchBand(i))
             IF (MatchBand(i) < 25) THEN
                nBanks = 1
             ELSE
                IF (MatchBand(i) == BandSwitch(1)) THEN
                   nBanks = 1
                ELSE
                   nBanks = 0
                ENDIF
             ENDIF

             DO n = 1, nBanks
                MaxAtten%DACS(BankIndx(n)) = AttenMaxed
                DeltaAtten%DACS(BankIndx(n)) = &
                     (BankAtten%DACS(BankIndx(n)) /= Val)
                BankAtten%DACS(BankIndx(n)) = Val
             ENDDO

          CASE (27:31)    ! MBs

             nBanks = 1
             BankIndx = MatchBand(i) - 26

             DO n = 1, nBanks
                MaxAtten%MB(BankIndx(n)) = AttenMaxed
                DeltaAtten%MB(BankIndx(n)) = &
                     (BankAtten%MB(BankIndx(n)) /= Val)
                BankAtten%MB(BankIndx(n)) = Val
             ENDDO

          CASE (32:34)    ! WFs

             nBanks = 1
             BankIndx = MatchBand(i) - 31

             DO n = 1, nBanks
                MaxAtten%WF(BankIndx(n)) = AttenMaxed
                DeltaAtten%WF(BankIndx(n)) = &
                     (BankAtten%WF(BankIndx(n)) /= Val)
                BankAtten%WF(BankIndx(n)) = Val
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
    THz_Sci_pkt%TSSM_pos = Sci_pkt%TSSM_pos
    THz_Sci_pkt%TSSM_theta = Sci_pkt%TSSM_theta
    THz_Sci_pkt%SwMirPos = Sci_pkt%THz_sw_pos
    THz_Sci_pkt%LLO_bias = LLO_Bias (Sci_pkt%LLO_DN)
    THz_Sci_pkt%BandSwitch = Sci_pkt%BandSwitch(4:5)  !Only need sw #4 & #5
    THz_Sci_pkt%CRC_good = Sci_pkt%CRC_good

  END SUBROUTINE Save_THz_pkt

  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE SciUtils

! $Log$
! Revision 2.24  2024/10/10 20:18:35  pwagner
! Improve comments; print GM07 when out-of-lock
!
! Revision 2.23  2023/06/06 22:40:55  pwagner
! Try to avoid printing binary NULL chars to stdout
!
! Revision 2.22  2022/11/08 23:47:59  pwagner
! Workaround for NAG signaling NaN in scAngleG and T
!
! Revision 2.21  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.20.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.20  2015/01/27 18:58:51  pwagner
! Print warning messages if NaNs will be left in OA products
!
! Revision 2.19  2011/06/07 18:56:03  perun
! Select Side B switching mirror table based on data time starting at 2011 DOY 153 19:10:00.
!
! Revision 2.18  2011/02/16 17:07:40  perun
! Estimate MIF 0 pointing when MIF 1, Position 2 is missing.
!
! Revision 2.17  2008/01/15 19:56:31  perun
! Allow unrecognized format type to pass through with a warning.
!
! Revision 2.16  2006/08/02 18:57:51  perun
! Allow MAF buffer to fill without all expected MIFs and fillin MIF 0 time data, if necessary
!
! Revision 2.15  2006/04/05 18:09:02  perun
! Remove unused variables
!
! Revision 2.14  2006/03/24 15:17:40  perun
! Initialize THzSciMAF MAF and MIF numbers to no data
!
! Revision 2.13  2005/08/24 15:53:18  perun
! Save pos1/pos2 prime data for the L1BOA file
!
! Revision 2.12  2005/08/11 19:06:02  perun
! Write bad checksum and band switching messages to log file
!
! Revision 2.11  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.10  2004/11/10 15:34:06  perun
! New value for TSSM_eps; test checksum sci value for correctness
!
! Revision 2.9  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.8  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
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
! Revision 2.24  2024/10/10 20:18:35  pwagner
! Improve comments; print GM07 when out-of-lock
!
! Revision 2.23  2023/06/06 22:40:55  pwagner
! Try to avoid printing binary NULL chars to stdout
!
! Revision 2.22  2022/11/08 23:47:59  pwagner
! Workaround for NAG signaling NaN in scAngleG and T
!
! Revision 2.21  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.20.4.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.20  2015/01/27 18:58:51  pwagner
! Print warning messages if NaNs will be left in OA products
!
! Revision 2.19  2011/06/07 18:56:03  perun
! Select Side B switching mirror table based on data time starting at 2011 DOY 153 19:10:00.
!
! Revision 2.18  2011/02/16 17:07:40  perun
! Estimate MIF 0 pointing when MIF 1, Position 2 is missing.
!
! Revision 2.17  2008/01/15 19:56:31  perun
! Allow unrecognized format type to pass through with a warning.
!
! Revision 2.16  2006/08/02 18:57:51  perun
! Allow MAF buffer to fill without all expected MIFs and fillin MIF 0 time data, if necessary
!
! Revision 2.15  2006/04/05 18:09:02  perun
! Remove unused variables
!
! Revision 2.14  2006/03/24 15:17:40  perun
! Initialize THzSciMAF MAF and MIF numbers to no data
!
! Revision 2.13  2005/08/24 15:53:18  perun
! Save pos1/pos2 prime data for the L1BOA file
!
! Revision 2.12  2005/08/11 19:06:02  perun
! Write bad checksum and band switching messages to log file
!
! Revision 2.11  2005/06/23 18:41:36  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.10  2004/11/10 15:34:06  perun
! New value for TSSM_eps; test checksum sci value for correctness
!
! Revision 2.9  2004/08/12 13:51:51  perun
! Version 1.44 commit
!
! Revision 2.8  2004/05/14 15:59:11  perun
! Version 1.43 commit
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
