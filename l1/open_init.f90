! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE OpenInit ! Opens input L0 files and output L1 files
!=============================================================================

  USE SDPToolkit
  USE MLSCommon
  USE MLSL1Common
  USE MLSL1Rad
  USE MLSMessageModule
  USE MLSPCF
  USE Hdf
  USE OutputL1B
  USE MLSSignalNomenclature

  IMPLICIT NONE
  PUBLIC

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = & 
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE OpenAndInitialize (processingRange, l0Info, l1bInfo, L1Brad)
!=============================================================================

    ! Arguments
    
    TYPE (TAI93_Range_T), INTENT (OUT) :: processingRange !Data processing range
    TYPE (L0Info_T), INTENT (OUT) :: l0Info !Lvl0 file info
    TYPE (L1BInfo_T), INTENT (OUT):: l1bInfo  !File handles etc. for L1B dataset
    TYPE (Radiance_T), DIMENSION (:), POINTER :: L1Brad

    INTEGER :: pgs_io_l0_open

    INTEGER, PARAMETER :: sci_pcf = 21000  ! temporary!!!
    INTEGER, PARAMETER :: nomen_pcf = 323  ! temporary!!!

    INTEGER :: returnStatus, error, unit
    INTEGER :: FileHandle, sd_id, ifl1, version
    CHARACTER (LEN=132) :: PhysicalFilename

    ! Open L0 science file (for Version 0.5 and above, 3 files will be needed)

    version = 1
    returnStatus = PGS_PC_getReference (sci_pcf, version, &
         & PhysicalFilename)

    IF (returnStatus /= PGS_S_SUCCESS) THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Error finding L0 Science file "//PhysicalFilename)

    END IF

    returnStatus = pgs_io_l0_open (sci_pcf, PGSd_EOS_PM, l0Info%sci, &
         & processingRange%startTime, processingRange%endTime)

    IF (returnStatus /= PGS_S_SUCCESS) THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Error opening L0 Science file "//PhysicalFilename)

    END IF

    l0Info%SciFileName = PhysicalFileName
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & "Opened L0 Science file "//PhysicalFilename)

    ifl1 = 0

    ! Open L1 RAD files
    ! Get the l1 file name from the PCF

    ! Will do only 1 for now!

    ALLOCATE (l1bInfo%L1BRADIDs(1), STAT=error)
    ALLOCATE (l1bInfo%L1BRADFileNames(1), STAT=error)

    DO FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_start

       version = 1
       returnStatus = Pgs_pc_getReference (FileHandle, version, &
            & PhysicalFilename)

       IF (returnStatus == PGS_S_SUCCESS) THEN

          ! Open the HDF file and initialize the SD interface

          sd_id = sfstart (PhysicalFilename, DFACC_CREATE)
          IF (sd_id == -1) THEN
             CALL MLSMessage (MLSMSG_Error, &
                  & ModuleName, "Error creating L1RAD file"//PhysicalFilename)
          ELSE
             ifl1 = ifl1 + 1
             l1bInfo%L1BRADIDs(ifl1) = sd_id
             l1bInfo%L1BRADFileNames(ifl1) = PhysicalFilename
             CALL MLSMessage (MLSMSG_Info, &
                  & ModuleName, "Opened L1RAD file"//PhysicalFilename)
          ENDIF
       END IF
    END DO ! FileHandle = mlspcf_l1b_rad_start, mlspcf_l1b_rad_end

    IF (ifl1 == 0)THEN

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find any L1BRAD file entries to create!")

    END IF

    ! Open L1BOA File

    version = 1
    returnStatus = Pgs_pc_getReference (mlspcf_l1b_oa_start, version, &
         & PhysicalFilename)

    IF (returnStatus == PGS_S_SUCCESS) THEN

       ! Open the HDF file and initialize the SD interface

       sd_id = sfstart (PhysicalFilename, DFACC_CREATE)
       IF (sd_id == -1) THEN
          CALL MLSMessage (MLSMSG_Error, &
               & ModuleName, "Error creating L1BOA file"//PhysicalFilename)
       ELSE
          CALL MLSMessage (MLSMSG_Info, &
               & ModuleName, "Opened L1BOA file"//PhysicalFilename)
          l1bInfo%L1BOAID = sd_id
          l1bInfo%L1BOAFileName = PhysicalFilename
       END IF

    ELSE

       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not find L1BOA file entry")

    END IF

    ! Open and Read nomenclature file

    version = 1
    returnStatus = PGS_IO_Gen_openF (nomen_pcf, PGSd_IO_Gen_RSeqFrm, 0, unit, &
         & version)
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not open nomenclature file")
    END IF

    CALL ReadSignalsDatabase (unit)

    returnStatus = PGS_IO_Gen_closeF (unit)  ! close unit
    IF (returnstatus /= PGS_S_SUCCESS) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & "Could not close nomenclature file")
    END IF

    ! Initialize the radiance vectors

    CALL InitRad (L1Brad)

    ! Define the SD structures in the output files

    CALL OutputL1B_create(l1bInfo)

    RETURN

  END Subroutine OpenAndInitialize

!=============================================================================
  FUNCTION ReadL0Packet (fileHandle, buf_len, pkt_buf, TAI93) RESULT (ret_len)
!=============================================================================

    INTEGER, INTENT (IN) :: fileHandle
    INTEGER, INTENT (IN) :: buf_len   ! sized large enough for largest packet!
    CHARACTER (LEN=*), INTENT (OUT) :: pkt_buf
    DOUBLE PRECISION, INTENT (OUT) :: TAI93

    INTEGER :: ret_len   ! actual packet length read

    INTEGER :: returnstatus
    INTEGER :: ccsds_len
    LOGICAL, SAVE :: more = .TRUE.  ! more to read next call

    INTEGER :: pgs_io_l0_getpacket
    INTEGER :: pgs_td_eospmtotai

    ret_len = 0   ! Nothing read yet

    IF (more) THEN
       returnstatus = pgs_io_l0_getpacket (fileHandle, buf_len, pkt_buf)
    ELSE
       RETURN
    END IF

    ! return length of CCSDS data length plus header

    ccsds_len = ichar(pkt_buf(5:5)) * 256 + ichar (pkt_buf(6:6)) + 7

    IF (returnstatus /= PGS_S_SUCCESS) THEN
       IF (returnstatus == PGSIO_W_L0_END_OF_VIRTUAL_DS) THEN
          ret_len = ccsds_len
          more = .FALSE.   ! all done reading for now
       ELSE   ! Eventually will handle multiple files
       END IF
    ELSE
       ret_len = ccsds_len
    END IF

    returnstatus = PGS_TD_EOSPMtoTAI (pkt_buf(8:15), TAI93)

  END FUNCTION ReadL0Packet

!=============================================================================
END MODULE OpenInit
!=============================================================================

!
! $Log$
! Revision 2.0  2000/09/05 18:55:15  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.3  2000/02/16 21:41:52  perun
! Opened and read Signals Database
!

