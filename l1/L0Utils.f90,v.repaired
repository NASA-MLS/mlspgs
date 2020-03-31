! Copyright 2008, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!=============================================================================
MODULE L0Utils ! Utilities to read L0 data
!=============================================================================

  USE MLSCommon, ONLY: r8
  USE MLSL1Common, ONLY: L0FileInfo
  USE SDPToolkit, ONLY: PGS_S_SUCCESS, PGSIO_W_L0_END_OF_VIRTUAL_DS, &
       PGSIO_M_L0_HEADER_CHANGED, PGS_PC_GetReference,PGSd_PC_FILE_PATH_MAX
  USE OpenInit, ONLY: OpenL0File
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Info, MLSMSG_Warning

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: ReadL0Sci, ReadL0Eng, CheckSum

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

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
    INTEGER :: ret_len, sindx
    LOGICAL :: EOD
    CHARACTER(len=132) :: filename

    INTEGER, PARAMETER :: SciLen = 1024
    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    TAI93 = 0.0
    sindx = 1
    DO

       ret_len = ReadL0Packet (L0FileInfo%sci_unit(sindx), SciLen, &
            SciPkt(sindx), TAI93(sindx), EOD)

       IF (TAI93(sindx) > L1Config%Expanded_TAI%endTime) THEN
          OK = .FALSE.
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & 'After time range: TAI time > L1Config%Expanded_TAI%endTime')
          ! so it gets sent to STDOUT too, otherwise we have to look in two
          ! files for the complete story
          print *,'After time range: TAI time > L1Config%Expanded_TAI%endTime'
          print *,'returning OK=FALSE'
          RETURN
       ENDIF

       IF (TAI93(sindx) < L1Config%Expanded_TAI%startTime) THEN
          OK = .FALSE.
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & 'Before time range: TAI time < L1Config%Expanded_TAI%startTime')
          ! so it gets sent to STDOUT too, otherwise we have to look in two
          ! files for the complete story

          PRINT *,'Before time range: TAI time < L1Config%Expanded_TAI%startTime'
          print *,'returning OK=TRUE nevertheless'
          ! RETURN
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
             print *,'No next L0 Sci File! returning OK=FALSE'
             OK = .FALSE.
             RETURN
          ENDIF

          CALL OpenL0File (L0FileInfo%sci_pcf(sindx), &
               L0FileInfo%sci_unit(sindx), L0FileInfo%SciFilename(sindx), &
               "Science")
       ENDIF

       IF (ret_len /= SciLen) THEN   ! incorrect length
          CALL MLSMessage (MLSMSG_Warning, ModuleName, &
               'Incorrect Sci packet length')
          CYCLE                      ! try again
       ENDIF

       IF (TAI93(1) == TAI93(2)) EXIT  ! Same packet time, so done for now

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
    USE MLSL1Common, ONLY: L1BFileInfo,FileNameLen
    USE SDPToolkit, ONLY: PGS_TD_TAItoUTC

    CHARACTER(LEN=*), DIMENSION(:) :: engpkt
    INTEGER :: MAFno, TotalMAF, MIFsPerMAF
    REAL(r8) :: MAFtime
    LOGICAL :: more_data, data_OK

    INTEGER :: i, n, returnStatus, MAF(6), version
    REAL(r8) :: TAI93(6), engtime
    INTEGER :: IDN(128), ret_len
    LOGICAL :: EOD
    CHARACTER(len=27) :: asciiUTC
    CHARACTER(len=FileNameLen) :: filename, msg


! For Version 3.0 and above:
    !offsets into packet giving the MAF number for this packet
    INTEGER, PARAMETER :: MAF_offset(6) = (/ 61, 253, 253, 17, 19, 17 /)

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    more_data = .TRUE.
    data_OK = .TRUE.
    MAF = -1
    TAI93 = -1.0     ! Init all times to less than engtime
    engtime = 0.0    ! No engtime, yet
    i = 1   ! start with eng packet #1
    DO

       IF (TAI93(i) < engtime) THEN

          ret_len = ReadL0Packet (L0FileInfo%eng_unit(i), LEN(engpkt(i)), &
               engpkt(i), TAI93(i), EOD)

          IF (engtime == 0.0) engtime = TAI93(i)

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

       ENDIF

       MAF(i) = BigEndianStr (EngPkt(i)(MAF_offset(i):MAF_offset(i)+1))

       If (TAI93(i) > engtime) THEN
          engtime = TAI93(i)              ! consider "new" time to match
          i = 1                           ! start with packet #1 again
       ELSEIF (TAI93(i) == engtime) THEN
          i = i + 1                       ! try next packet
       ENDIF

       IF (i > 6) EXIT

    ENDDO

    MIFsPerMAF = BigEndianStr (EngPkt(6)(244:245))
    TotalMAF = BigEndianStr (EngPkt(6)(246:249))
    MAFno = MAF(1)
    MAFtime = engtime

! Test checksum value

    DO i = 1, 6
       DO n = 1, 128
          IDN(n) = BigEndianStr (engpkt(i)((n*2-1):n*2))
       ENDDO
       IF (CheckSum (IDN, 127) /= IDN(128)) THEN
          data_OK = .FALSE.
          n = PGS_TD_TAItoUTC (engtime, asciiUTC)
          WRITE (msg, &
               '("Bad Checksum: Eng Pkt ", i1, ",MAF: ",i4,' // &
               '", UTC: ", A27)') i, MAFno, asciiUTC


          PRINT *, TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
          WRITE (L1BFileInfo%LogId, *) '### Info: '//TRIM(msg)
          WRITE (L1BFileInfo%LogId, *) ''
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
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE L0Utils
!=============================================================================

! $Log$
! Revision 2.17  2020/03/31 20:47:23  pwagner
! Permit timepackets prior to start time; restores MAFs lost by v5.00
!
! Revision 2.16  2018/04/17 16:19:26  whdaffer
! line 265, syntax error for NAG in write statement
!
! Revision 2.15  2018/04/10 18:14:21  whdaffer
! Gave output field a width, to mollify NAG compiler
!
! Revision 2.14  2018/04/09 22:13:21  whdaffer
! Reportage
!
! Revision 2.13  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.12.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.12  2008/03/18 17:20:14  perun
! Align Sci packets based on time, not MIF number.
!
! Revision 2.11  2008/02/14 15:00:05  perun
! Changed ReadL0Eng to accommodate missing APID data.
!
! Revision 2.10  2006/08/02 18:54:58  perun
! Warn if science packet of incorrect length is read and continue processing
!
! Revision 2.9  2005/08/11 19:03:11  perun
! Write bad checksum message to log file
!
! Revision 2.8  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.7  2005/02/28 17:15:10  perun
! Corrected loop for testing checksum
!
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
