! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE SciUtils ! L0 science utilities
!=============================================================================

  USE L0_sci_tbls, ONLY: l0_sci1, l0_sci2, sci_type, sci1_uc_fmt, sci2_uc_fmt, &
       sci1_uc, sci2_uc, sci_cptr, MAF_uc_cptr, MIF_uc_cptr, orbit_uc_cptr, &
       fb_uc_cptr, mb_uc_cptr, GHz_ant_scan_uc_cptr, GHz_sw_uc_cptr, sci_pkt, &
       THz_sw_uc_cptr, SciMAF, FBnum, MBnum, type_sci1
  USE L0Utils, ONLY: ReadL0Sci
  USE MLSL1Utils, ONLY: BigEndianStr, ExtractBigEndians, SwapBytes, QNan, &
       Finite
  USE MLSL1Common, ONLY: deg24, GHz_SwMir_range, THz_SwMir_range, &
       SwMir_Range_T

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

    LOGICAL OK

    INTEGER :: i, ios, returnStatus
    INTEGER, SAVE :: old_type = -1

    !! Pointers to counts:

    INTEGER, DIMENSION(:), POINTER :: fbcnts => NULL()
    INTEGER, DIMENSION(:), POINTER :: mbcnts => NULL()
    CHARACTER(len=1), PARAMETER :: sw_good = CHAR(7)

    CHARACTER (LEN=1024) :: scipkt(2)

    INTEGER, EXTERNAL :: PGS_TD_EOSPMGIRDtoTAI

    CALL ReadL0Sci (scipkt, OK)
    l0_sci1 = scipkt(1)
    l0_sci2 = scipkt(2)

    IF (.NOT. OK) RETURN

    sci_type = ICHAR(type_sci1)  ! science packet #1 type

    IF (sci_type == 3) THEN      ! uncompressed format(s)

       READ (UNIT=l0_sci1, FMT=sci1_uc_fmt, iostat=ios) sci1_uc
       READ (UNIT=l0_sci2, FMT=sci2_uc_fmt, iostat=ios) sci2_uc

    ENDIF

    IF (sci_type /= old_type) THEN

       old_type = sci_type      ! Keep old type for future tests

       IF (sci_type == 3) THEN  ! sci pkt 1, type ii (uncompressed)

          !! Save uncompressed pointers:

          sci_cptr%MAF = MAF_uc_cptr
          sci_cptr%MIF = MIF_uc_cptr
          sci_cptr%orbit = orbit_uc_cptr

          sci_cptr%FB = FB_uc_cptr
          sci_cptr%MB = MB_uc_cptr

          sci_cptr%GHz_ant_scan => GHz_ant_scan_uc_cptr
          sci_cptr%GHz_sw => GHz_sw_uc_cptr
          sci_cptr%THz_sw => THz_sw_uc_cptr

       ENDIF

    ENDIF

!! Put some of the data into the correct order and convert to angles:

    CALL SwapBytes (sci_cptr%THz_sw, sci_cptr%THz_sw)
    IF (sci_cptr%THz_sw(1:1) == sw_good) THEN
       Sci_pkt%THz_sw_angle(1)  = deg24 * BigEndianStr (sci_cptr%THz_sw(3:5))
       Sci_pkt%THz_sw_angle(2)  = deg24 * BigEndianStr (sci_cptr%THz_sw(6:8))
    ELSE
       Sci_pkt%THz_sw_angle(:) = QNan()
    ENDIF
    Sci_pkt%THz_sw_pos = SwMirPos ("T", Sci_pkt%THz_sw_angle)

    CALL SwapBytes (sci_cptr%GHz_sw, sci_cptr%GHz_sw)
    IF (sci_cptr%GHz_sw(1:1) == sw_good) THEN
       Sci_pkt%GHz_sw_angle(1)  = deg24 * BigEndianStr (sci_cptr%GHz_sw(3:5))
       Sci_pkt%GHz_sw_angle(2)  = deg24 * BigEndianStr (sci_cptr%GHz_sw(6:8))
    ELSE
       Sci_pkt%GHz_sw_angle(:) = QNan()
    ENDIF
    Sci_pkt%GHz_sw_pos = SwMirPos ("G", Sci_pkt%GHz_sw_angle)

!! Convert raw data:

    Sci_pkt%MAFno = BigEndianStr (sci_cptr%MAF(1)%ptr)
    Sci_pkt%MIFno = BigEndianStr (sci_cptr%MIF(1)%ptr)
    Sci_pkt%Orbit = BigEndianStr (sci_cptr%orbit(1)%ptr)

! Get TAI93 time

    returnstatus = PGS_TD_EOSPMGIRDtoTAI (scipkt(1)(8:15), Sci_pkt%secTAI)

    DO i = 1, FBnum
       CALL ExtractBigEndians (sci_cptr%FB(i)%ptr, fbcnts)
       Sci_pkt%FB(:,i) = fbcnts
    ENDDO

    DO i = 1, MBnum
       CALL ExtractBigEndians (sci_cptr%MB(i)%ptr, mbcnts)
       Sci_pkt%MB(:,i) = mbcnts
    ENDDO

!! Check for good checksums (LATER!!!):

    Sci_pkt%CRC_good = .TRUE.

    OK = .TRUE.

  END FUNCTION GetSciPkt

!=============================================================================
  SUBROUTINE NextSciMAF (more_data)
!=============================================================================

    !! Get the next MAF's science data

    LOGICAL, INTENT (OUT) :: more_data

    INTEGER, PARAMETER :: no_data = -1   ! no data is available
    INTEGER, SAVE :: prev_MAF = no_data

    INTEGER :: last_MIFno = 147   !! Nominal last MIFno minus 1 (TEST!!!)
    more_data = .TRUE.

    !! Initialize MAF/MIF counters to indicate no data

    SciMAF%MAFno = no_data
    SciMAF%MIFno = no_data

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
                  SciMAF(last_MIFno)%MIFno /= no_data) THEN
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

    CHARACTER(len=1), INTENT(IN) :: sw_module
    REAL, INTENT(IN) :: angle(2)
    CHARACTER(len=2) :: sw_pos

    INTEGER :: n
    TYPE (SwMir_Range_T), DIMENSION(:), POINTER :: SwMir_Range

    sw_pos = "D"                 ! Initialize to "Discard"
    IF (.NOT. Finite (angle(1)) .OR. .NOT. Finite (angle(2))) RETURN

    IF (sw_module == "G") THEN   ! GigaHertz Module
       SwMir_Range => GHz_SwMir_range
    ELSE                         ! TeraHertz Module
       SwMir_Range => THz_SwMir_range
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

END MODULE SciUtils

! $Log$
! Revision 2.1  2001/02/23 20:56:11  perun
! Version 0.5 commit
!
