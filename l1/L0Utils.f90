! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE L0Utils ! Utilities to read L0 data
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: L0FileInfo
  USE SDPToolkit, ONLY: PGS_S_SUCCESS, PGSIO_W_L0_END_OF_VIRTUAL_DS, &
       PGSIO_M_L0_HEADER_CHANGED, PGS_PC_GetReference
  USE OpenInit, ONLY: OpenL0File
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMSG_Warning, MLSMSG_Error

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
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

  SUBROUTINE ReadL0Sci (scipkt, OK)

    USE MLSL1Config, ONLY: L1Config

    CHARACTER(LEN=*), DIMENSION(:) :: scipkt
    LOGICAL :: OK

    INTEGER :: i, returnStatus, version
    REAL(r8) :: TAI93
    INTEGER :: ret_len
    LOGICAL :: EOD
    CHARACTER(len=132) :: filename

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    DO i = 1, 2

       ret_len = ReadL0Packet (L0FileInfo%sci_unit(i), LEN(scipkt(i)), &
            scipkt(i), TAI93, EOD)

       IF (TAI93 > L1Config%Expanded_TAI%endTime) THEN
          OK = .FALSE.
          RETURN
       ENDIF

       IF (EOD) THEN    ! May need to do this elsewhere or pass flag here
          EOD = .FALSE.       
          returnStatus = PGS_IO_L0_Close (L0FileInfo%sci_unit(i))

          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & 'Closed L0 Science file: '//L0FileInfo%SciFileName(i))

          L0FileInfo%sci_pcf(i) = L0FileInfo%sci_pcf(i) + 1  ! Next entry
          version = 1
          returnStatus = PGS_PC_getReference (L0FileInfo%sci_pcf(i), version, &
               filename)
          IF (returnStatus /= PGS_S_SUCCESS) THEN
             CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                  & 'No next L0 Sci File')
             OK = .FALSE.
             RETURN
          ENDIF

          CALL OpenL0File (L0FileInfo%sci_pcf(i), L0FileInfo%sci_unit(i), &
               L0FileInfo%SciFilename(i), "Science")
       ENDIF

    ENDDO

    OK = .NOT. EOD   ! TEST!!!

  END SUBROUTINE ReadL0Sci

  SUBROUTINE ReadL0Eng (engpkt, MAFno, TotalMAF, MIFsPerMAF, OK)

    USE MLSL1Utils, ONLY: BigEndianStr

    CHARACTER(LEN=*), DIMENSION(:) :: engpkt
    INTEGER :: MAFno, TotalMAF, MIFsPerMAF
    LOGICAL :: OK

    INTEGER :: i, returnStatus, MAF(6), version
    REAL(r8) :: TAI93, engtime
    INTEGER :: ret_len
    LOGICAL :: EOD
    CHARACTER(len=132) :: filename

! For Version 3.0:

    INTEGER, PARAMETER :: MAF_offset(6) = (/ 61, 253, 253, 17, 19, 17 /)

! For Version 2.4:

    ! INTEGER, PARAMETER :: MAF_offset(6) = (/ 61, 241, 252, 19, 19, 17 /)

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    DO i = 1, 6

       DO

          ret_len = ReadL0Packet (L0FileInfo%eng_unit(i), LEN(engpkt(i)), &
               engpkt(i), TAI93, EOD)

          IF (i == 1) engtime = TAI93   ! Save time from packet #1 for comparisons

          IF (EOD) THEN    ! May need to do this elsewhere or pass flag here

             returnStatus = PGS_IO_L0_Close (L0FileInfo%eng_unit(i))

             CALL MLSMessage (MLSMSG_Info, ModuleName, &
                  & 'Closed L0 Engineering file: '//L0FileInfo%EngFileName(i))

             L0FileInfo%eng_pcf(i) = L0FileInfo%eng_pcf(i) + 1  ! Next entry
             version = 1
             returnStatus = PGS_PC_getReference (L0FileInfo%eng_pcf(i), version, &
                  filename)
             IF (returnStatus /= PGS_S_SUCCESS) THEN
                CALL MLSMessage (MLSMSG_Warning, ModuleName, &
                     & 'No next L0 Eng File')
                OK = .FALSE.
                RETURN
             ENDIF

             CALL OpenL0File (L0FileInfo%eng_pcf(i), L0FileInfo%eng_unit(i), &
                  L0FileInfo%EngFilename(i), "Engineering")
          ENDIF

          IF (TAI93 >= engtime) EXIT  ! Put in time order for packets 1-6

       ENDDO

       MAF(i) = BigEndianStr (EngPkt(i)(MAF_offset(i):MAF_offset(i)+1))

    ENDDO

    MIFsPerMAF = BigEndianStr (EngPkt(6)(244:245))
    TotalMAF = BigEndianStr (EngPkt(6)(246:249))

    OK = .TRUE.

    MAFno = MAF(1)

    DO i = 2, 6
       IF (MAF(i) /= MAFno) THEN  ! All must be the same MAF number
          OK = .FALSE.
          EXIT
       ENDIF
    ENDDO

  END SUBROUTINE ReadL0Eng

!=============================================================================
END MODULE L0Utils
!=============================================================================

! $Log$
! Revision 2.2  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.1  2001/02/23 20:48:47  perun
! Version 0.5 commit
!
