! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Close_files ! Close the production files
!=============================================================================

  USE MLSL1Common, ONLY: L0FileInfo, L1BFileInfo
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info
  USE Hdf, ONLY: sfend
  USE WriteMetaL1

  IMPLICIT NONE

  PRIVATE :: Id, ModuleName
  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CloseFiles  ! Close production Level 1 files
!=============================================================================

    USE OpenInit, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr

    INTEGER :: i, returnStatus

    INTEGER, EXTERNAL :: PGS_IO_L0_Close

    ! Close L0 Science Files:

    DO i = 1, 2
       
       returnStatus = PGS_IO_L0_Close (L0FileInfo%sci_unit(i))

       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L0 Science file: '//L0FileInfo%SciFileName(i))

    ENDDO

    ! Close L0 Engineering Files:

    DO i = 1, 6
       
       returnStatus = PGS_IO_L0_Close (L0FileInfo%eng_unit(i))

       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L0 Engineering file: '//L0FileInfo%EngFileName(i))

    ENDDO

    ! Write Annotations and Close L1RAD D file

    CALL WritePCF2Hdr (l1BFileInfo%RADDFileName, anTextPCF)

    CALL WritePCF2Hdr (l1BFileInfo%RADDFileName, anTextCF)

    returnStatus = sfend (L1BFileInfo%RADDid)
    IF (returnStatus == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & 'Failed to close L1BRAD D file: '//l1BFileInfo%RADDFileName)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)

    ! Write Annotations and Close L1RAD F file

    CALL WritePCF2Hdr (l1BFileInfo%RADFFileName, anTextPCF)

    CALL WritePCF2Hdr (l1BFileInfo%RADFFileName, anTextCF)

    returnStatus = sfend (L1BFileInfo%RADFid)
    IF (returnStatus == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & 'Failed to close L1BRAD F file: '//l1BFileInfo%RADFFileName)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD F file: '//L1BFileInfo%RADFFileName)

    ! Write Annotations and Close L1BOA file

    CALL WritePCF2Hdr (l1BFileInfo%OAFileName, anTextPCF)

    CALL WritePCF2Hdr (l1BFileInfo%OAFileName, anTextCF)

    returnStatus = sfend (L1BFileInfo%OAid)
    IF (returnStatus == -1) THEN
       CALL MLSMessage (MLSMSG_Error, ModuleName, &
            & 'Failed to close L1BOA file: '//l1BFileInfo%OAFileName)
    ENDIF

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

    CALL WriteMetaData

  END Subroutine CloseFiles

!=============================================================================
END MODULE Close_files
!=============================================================================
! $Log$
! Revision 2.4  2001/03/12 19:40:11  perun
! Write PCF and CF contents as annotations
!
! Revision 2.2  2001/03/05 22:35:13  perun
! Corrected L0_Close call
!
! Revision 2.1  2001/02/23 18:51:39  perun
! Version 0.5 commit
!
