! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L0Utils ! Utilities to read L0 data
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: L0FileInfo
  USE SDPToolkit, ONLY: PGS_S_SUCCESS, PGSIO_W_L0_END_OF_VIRTUAL_DS, &
       PGSIO_M_L0_HEADER_CHANGED, PGS_PC_GetReference
  USE OpenInit, ONLY: OpenL0File
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMSG_Warning

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ReadL0Sci, ReadL0Eng, CheckSum

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  FUNCTION ReadL0Packet (fileHandle, buf_len, pkt_buf, TAI93, EOD) &
       & RESULT (ret_len)
!=============================================================================

    INTEGER, INTENT (IN) :: fileHandle
    INTEGER, INTENT (IN) :: buf_len   ! sized large enough for largest packet!
    CHARACTER (LEN=*), INTENT (OUT) :: pkt_buf
    REAL(r8), INTENT (OUT) :: TAI93
    LOGICAL :: EOD    ! End-Of-Data flag

    INTEGER :: ret_len   ! actual packet length read

    INTEGER :: returnstatus

    INTEGER, EXTERNAL :: PGS_IO_L0_Getpacket
    INTEGER, EXTERNAL :: PGS_TD_EOSPMGIRDtoTAI

    ret_len = 0     ! Nothing read yet
    EOD = .FALSE.   ! Not at the end-of-data

    returnstatus = PGS_IO_L0_Getpacket (fileHandle, buf_len, pkt_buf)

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       IF (returnstatus == PGSIO_W_L0_END_OF_VIRTUAL_DS) THEN
          EOD = .TRUE.
       ELSE IF (returnStatus /= PGSIO_M_L0_HEADER_CHANGED) THEN
          RETURN   ! Nothing more to do
       ENDIF
    ENDIF

    ! Get TAI93 time

    returnstatus = PGS_TD_EOSPMGIRDtoTAI (pkt_buf(8:15), TAI93)

    ! return length of CCSDS data length plus header

    ret_len = ICHAR (pkt_buf(5:5)) * 256 + ICHAR (pkt_buf(6:6)) + 7

  END FUNCTION ReadL0Packet

!=============================================================================
  SUBROUTINE ReadL0Sci (SciPkt, OK)
!=============================================================================

    USE MLSL1Config, ONLY: L1Config

    CHARACTER(LEN=*), DIMENSION(:) :: SciPkt
    LOGICAL :: OK

    INTEGER :: returnStatus, version
    REAL(r8) :: TAI93(2)
    INTEGER :: ret_len, sindx, MIF(2)
    LOGICAL :: EOD
    CHARACTER(len=132) :: filename

    INTEGER, PARAMETER :: SciLen = 1024
    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    MIF = -1
    TAI93 = 0.0
    sindx = 1
    DO

       ret_len = ReadL0Packet (L0FileInfo%sci_unit(sindx), SciLen, &
            SciPkt(sindx), TAI93(sindx), EOD)

       IF (TAI93(sindx) > L1Config%Expanded_TAI%endTime) THEN
          OK = .FALSE.
          RETURN
       ENDIF

       IF (EOD) THEN    ! May need to do this elsewhere or pass flag here
          EOD = .FALSE.       
          returnStatus = PGS_IO_L0_Close (L0FileInfo%sci_unit(sindx))

          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & 'Closed L0 Science file: '//L0FileInfo%SciFileName(sindx))

          L0FileInfo%sci_pcf(sindx) = L0FileInfo%sci_pcf(sindx) + 1 ! Next entry
          version = 1
          returnStatus = PGS_PC_getReference (L0FileInfo%sci_pcf(sindx), &
               version, filename)
          IF (returnStatus /= PGS_S_SUCCESS) THEN
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & 'No next L0 Sci File')
             OK = .FALSE.
             RETURN
          ENDIF

          CALL OpenL0File (L0FileInfo%sci_pcf(sindx), &
               L0FileInfo%sci_unit(sindx), L0FileInfo%SciFilename(sindx), &
               "Science")
       ENDIF

       MIF(sindx) = ICHAR (SciPkt(sindx)(17:17))

       IF (MIF(1) == MIF(2)) EXIT

       IF (TAI93(2) < TAI93(1)) THEN
          sindx = 2
       ELSE
          sindx = 1
       ENDIF

    ENDDO

    OK = .NOT. EOD   ! TEST!!!

  END SUBROUTINE ReadL0Sci

!=============================================================================
  SUBROUTINE ReadL0Eng (engpkt, MAFno, TotalMAF, MIFsPerMAF, MAFtime, data_Ok, &
       more_data)
!=============================================================================

    USE MLSL1Utils, ONLY: BigEndianStr
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

    CHARACTER(LEN=*), DIMENSION(:) :: engpkt
    INTEGER :: MAFno, TotalMAF, MIFsPerMAF
    REAL(r8) :: MAFtime
    LOGICAL :: more_data, data_OK

    INTEGER :: i, n, returnStatus, MAF(6), version
    REAL(r8) :: TAI93, engtime
    INTEGER :: IDN(128), ret_len
    LOGICAL :: EOD
    CHARACTER(len=27) :: asciiUTC
    CHARACTER(len=132) :: filename, msg

! For Version 3.0 and above:

    INTEGER, PARAMETER :: MAF_offset(6) = (/ 61, 253, 253, 17, 19, 17 /)

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    more_data = .TRUE.
    data_OK = .TRUE.
    MAF = -1
    i = 1   ! start with eng packet #1
    DO

       DO

          ret_len = ReadL0Packet (L0FileInfo%eng_unit(i), LEN(engpkt(i)), &
               engpkt(i), TAI93, EOD)

          IF (i == 1) THEN
             engtime = TAI93 ! Save time from packet #1 for comparisons
             MAFtime = engtime
          ENDIF

          IF (EOD) THEN    ! May need to do this elsewhere or pass flag here

             returnStatus = PGS_IO_L0_Close (L0FileInfo%eng_unit(i))

             CALL MLSMessage (MLSMSG_Info, ModuleName, &
                  & 'Closed L0 Engineering file: '//L0FileInfo%EngFileName(i))

             L0FileInfo%eng_pcf(i) = L0FileInfo%eng_pcf(i) + 1  ! Next entry
             version = 1
             returnStatus = PGS_PC_getReference (L0FileInfo%eng_pcf(i), &
                  version, filename)
             IF (returnStatus /= PGS_S_SUCCESS) THEN
                CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                     & 'No next L0 Eng File')
                more_data = .FALSE.
                RETURN
             ENDIF

             CALL OpenL0File (L0FileInfo%eng_pcf(i), L0FileInfo%eng_unit(i), &
                  L0FileInfo%EngFilename(i), "Engineering")
          ENDIF

         IF (TAI93 >= engtime) EXIT  ! Put in time order for packets 1-6

       ENDDO

       MAF(i) = BigEndianStr (EngPkt(i)(MAF_offset(i):MAF_offset(i)+1))

       IF (i == 1) THEN
          DO n = (i+1), 6
             IF (MAF(n) /= MAF(n-1)) EXIT
          ENDDO
       ELSE
          IF (MAF(i) == MAF(1)) THEN
             DO n = (i+1), 6
                IF (MAF(n) /= MAF(n-1)) EXIT
             ENDDO
          ELSE
             n = 1    ! Start with the next MAF
          ENDIF
       ENDIF
       i = n
       IF (i > 6) EXIT

    ENDDO

    MIFsPerMAF = BigEndianStr (EngPkt(6)(244:245))
    TotalMAF = BigEndianStr (EngPkt(6)(246:249))
    MAFno = MAF(1)

! Test checksum value

    DO i = 1, 6
       DO n = 1, 256
          IDN(n) = BigEndianStr (engpkt(i)((n*2-1):n*2))
       ENDDO
       IF (CheckSum (IDN, 127) /= IDN(128)) THEN
          data_OK = .FALSE.
          n = PGS_TD_TAItoUTC (engtime, asciiUTC)
          WRITE (msg, &
               '("Bad Checksum: Eng Pkt ", i1, ", UTC: ", A27)') i, asciiUTC
          PRINT *, TRIM(msg)
          CALL MLSMessage (MLSMSG_Info, ModuleName, TRIM(msg))
          RETURN   ! Can't do any more
       ENDIF
   ENDDO

  END SUBROUTINE ReadL0Eng

!=============================================================================
  FUNCTION CheckSum (DN, n) RESULT (csum)
!=============================================================================

    INTEGER :: n
    INTEGER :: DN(n)   ! contains 16-bit integer data

    INTEGER :: csum, i, seed
    DATA seed /z'A300'/

    csum = seed
    DO i = 1, n
       csum = IEOR (csum, DN(i))
    ENDDO

  END FUNCTION CheckSum

!=============================================================================
END MODULE L0Utils
!=============================================================================

! $Log$
! Revision 2.6  2004/11/10 15:33:11  perun
! Test checksum eng value for correctness
!
! Revision 2.5  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
! Revision 2.4  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
! Revision 2.3  2002/09/12 16:18:39  perun
! Align science packets to same MIF and MAF
!
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:48:47  perun
! Version 0.5 commit
!
