! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SciUtils ! L0 science utilities
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE L0_sci_tbls, ONLY: sci_type, l0_sci1, l0_sci2, sci1_T1_fmt, sci2_T1_fmt, &
       sci1_T2_fmt, sci2_T2_fmt, sci1_T1, sci2_T1, sci1_T2, sci2_T2, &
       sci_cptr, Sci_pkt, SciMAF, type_sci1
  USE L0Utils, ONLY: ReadL0Sci
  USE MLSL1Utils, ONLY: BigEndianStr, ExtractBigEndians, SwapBytes, QNan, &
       Finite
  USE DACsUtils, ONLY: ExtractDACSdata, UncompressDACSdata, ProcessDACSdata
  USE MLSL1Common, ONLY: FBnum, MBnum, WFnum, WFchans, deg24, DACSnum, &
       DACSchans, BandSwitch

  IMPLICIT NONE

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

  CONTAINS

!=============================================================================
  FUNCTION GetSciPkt () RESULT (OK)
!=============================================================================

    USE ERMSG_M, ONLY: ermset

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

    INTEGER :: DACS_K(DACSchans)
    INTEGER :: DACS_i1, DACS_i2
    INTEGER :: D(4), Dtot, TP, DIO, LO, nchans

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

!! Convert to angles:

    IF (sci_cptr(tindex)%THz_sw(1:1) == sw_good) THEN
       Sci_pkt%THz_sw_angle(1)  = deg24 * &
            BigEndianStr (sci_cptr(tindex)%THz_sw(3:5))
       Sci_pkt%THz_sw_angle(2)  = deg24 * &
            BigEndianStr (sci_cptr(tindex)%THz_sw(6:8))
    ELSE
       Sci_pkt%THz_sw_angle(:) = QNan()
    ENDIF
    Sci_pkt%THz_sw_pos = SwMirPos ("T", Sci_pkt%THz_sw_angle)

    IF (sci_cptr(tindex)%GHz_sw(1:1) == sw_good) THEN
       Sci_pkt%GHz_sw_angle(1)  = deg24 * &
            BigEndianStr (sci_cptr(tindex)%GHz_sw(3:5))
       Sci_pkt%GHz_sw_angle(2)  = deg24 * &
            BigEndianStr (sci_cptr(tindex)%GHz_sw(6:8))
    ELSE
       Sci_pkt%GHz_sw_angle(:) = QNan()
    ENDIF
    Sci_pkt%GHz_sw_pos = SwMirPos ("G", Sci_pkt%GHz_sw_angle)

!! Filter bank switches

    Sci_pkt%GSN = ICHAR (sci_cptr(tindex)%GSN%ptr)
    DO i = 1, 4
       CALL Band_switch (i, Sci_pkt%GSN(i), Sci_pkt%BandSwitch(i))
    ENDDO
    Sci_pkt%THzSw = ICHAR (sci_cptr(tindex)%THzSw%ptr)
    CALL Band_switch (i, Sci_pkt%THzSw, Sci_pkt%BandSwitch(5))

    DO i = 1, 5
       IF (Sci_pkt%BandSwitch(i) == 0) THEN
          Sci_pkt%BandSwitch(i) = BandSwitch(i)  ! use current value
       ELSE IF ((Sci_pkt%BandSwitch(i) > 0) .AND. &
            (BandSwitch(i) /= Sci_pkt%BandSwitch(i))) THEN
          BandSwitch(i) = Sci_pkt%BandSwitch(i)  ! new current value
       ENDIF
    ENDDO

!! Convert raw data:

    Sci_pkt%MAFno = BigEndianStr (sci_cptr(tindex)%MAF(1)%ptr)
    Sci_pkt%MIFno = BigEndianStr (sci_cptr(tindex)%MIF(1)%ptr)
    Sci_pkt%Orbit = BigEndianStr (sci_cptr(tindex)%orbit(1)%ptr)

! Get TAI93 time

    returnstatus = PGS_TD_EOSPMGIRDtoTAI (scipkt(1)(8:15), Sci_pkt%secTAI)

    DO i = 1, FBnum
       CALL ExtractBigEndians (sci_cptr(tindex)%FB(i)%ptr, FBcnts)
       Sci_pkt%FB(:,i) = FBcnts

    ENDDO

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

    IF (sci_type == type_I) THEN    ! compressed DACS

       DO i = 1, 4
          CALL UncompressDACSdata (sci_cptr(tindex)%DACS(i)%ptr)
       ENDDO
!print *, "compressed DACS"
    ELSE

       Sci_pkt%DACS = 0.0

       DO i = DACS_i1, DACS_i2

          CALL ExtractDACSdata (sci_cptr(tindex)%DACS(i)%ptr, DACS_K, D, TP, &
               DIO, LO, nchans)

          CALL ProcessDACSdata (D, DACS_K, nchans, TP, Sci_pkt%DACS(:,i))

       ENDDO

    ENDIF

!! LLO data

    Sci_pkt%LLO_data = sci_cptr(tindex)%LLO_data%ptr
!print '(16(1x,Z2.2))', ICHAR(Sci_pkt%LLO_data(1:16))

!! Check for good checksums (LATER!!!):

    Sci_pkt%CRC_good = .TRUE.

    OK = .TRUE.

  END FUNCTION GetSciPkt

!=============================================================================
  SUBROUTINE NextSciMAF (more_data)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    !! Get the next MAF's science data

    LOGICAL, INTENT (OUT) :: more_data

    INTEGER, PARAMETER :: no_data = -1   ! no data is available
    INTEGER, SAVE :: prev_MAF = no_data

    INTEGER :: last_MIF !! = 147   !! Nominal last MIFno minus 1 (TEST!!!)
    more_data = .TRUE.

    !! Initialize MAF/MIF counters to indicate no data

    SciMAF%MAFno = no_data
    SciMAF%MIFno = no_data
    last_MIF = L1Config%Calib%MIFsPerMAF -1

    !! Initialize CRC flags to not good

    SciMAF%CRC_good = .FALSE.

    !! Save previously read packet (if available):

    IF (prev_MAF /= no_data) THEN
       SciMAF(Sci_pkt%MIFno) = Sci_pkt
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

          prev_MAF = Sci_pkt%MAFno

       ELSE

          more_data = .FALSE.
          prev_MAF = no_data
          EXIT  ! Nothing more is available

       ENDIF
    ENDDO

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
       IF (angle(1) >= SwMir_Range(n)%low_angle .AND. &
            angle(1) <= SwMir_Range(n)%high_angle .AND. &
            angle(2) >=  SwMir_Range(n)%low_angle .AND. &
            angle(2) <=  SwMir_Range(n)%high_angle) THEN
          sw_pos = SwMir_Range(n)%pos
          EXIT
       ENDIF
    ENDDO

  END FUNCTION SwMirPos

  SUBROUTINE Band_switch (switch, sw_val, band)

    INTEGER :: band, sw_val, switch

    INTEGER :: i, pos

    INTEGER, PARAMETER :: pos_no(6) = (/ &
         z'FE', z'FD', z'FB', z'F7', z'EF', z'DF' /)

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

    IF (sw_val == z'FF') THEN  ! in process of changing
       band = -1
       print *, 'switching bands for switch ', switch
       RETURN
    ENDIF

    DO i = 1, 6
       IF (sw_val == pos_no(i)) EXIT
    ENDDO
    pos = i
    IF (pos > 6) THEN
       ! print *, 'illegal switch position!'
       ! print *, 'switch, sw_val: ', switch, sw_val
       ! stop
       band = 0
       RETURN
    ENDIF
    band = sw_bands (pos, switch)

  END SUBROUTINE Band_switch

END MODULE SciUtils

! $Log$
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:56:11  perun
! Version 0.5 commit
!
