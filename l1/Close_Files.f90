!
! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
!
!=============================================================================
MODULE Close_files ! Close the production files
!=============================================================================

  USE MLSL1Common, ONLY: L0FileInfo, L1BFileInfo
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info
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

    USE MLSL1Config, ONLY: L1Config
    USE OpenInit, ONLY: antextPCF, antextCF
    USE PCFHdr, ONLY: WritePCF2Hdr
    USE SDPToolkit, ONLY: PGS_IO_Gen_closeF
    USE MLSFiles, ONLY: mls_closeFile, HDFVERSION_4, HDFVERSION_5  

    INTEGER :: i, returnStatus, error, hdfVersion

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

! Close HDF files

    hdfVersion = L1Config%Output%HDFversion

    ! Write Annotations and Close L1RAD D file

    IF (hdfVersion == HDFVERSION_4) THEN
       CALL WritePCF2Hdr (l1BFileInfo%RADDFileName, anTextPCF)
       CALL WritePCF2Hdr (l1BFileInfo%RADDFileName, anTextCF)
    ENDIF

    CALL mls_closeFile (L1BFileInfo%RADDid, hdfVersion)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)

    ! Write Annotations and Close L1RAD F file

    IF (hdfVersion == HDFVERSION_4) THEN
       CALL WritePCF2Hdr (l1BFileInfo%RADFFileName, anTextPCF)
       CALL WritePCF2Hdr (l1BFileInfo%RADFFileName, anTextCF)
    ENDIF

    CALL mls_closeFile (L1BFileInfo%RADFid, hdfVersion)
    
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD F file: '//L1BFileInfo%RADFFileName)

    ! Write Annotations and Close L1BOA file

    IF (hdfVersion == HDFVERSION_4) THEN
       CALL WritePCF2Hdr (l1BFileInfo%OAFileName, anTextPCF)
       CALL WritePCF2Hdr (l1BFileInfo%OAFileName, anTextCF)
    ENDIF

    CALL mls_closeFile (L1BFileInfo%OAid, hdfVersion)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

    CALL WriteMetaData

    returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%EngId)
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
          CALL MLSMessage (MLSMSG_ERROR, ModuleName, &
               "Calling PGS_IO_Gen_CloseF failed." )
    ENDIF          

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BENG file: '//l1BFileInfo%EngFileName)

    returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%DiagId)
    IF (returnStatus /= PGS_S_SUCCESS) THEN 
          CALL MLSMessage (MLSMSG_ERROR, ModuleName, &
               "Calling PGS_IO_Gen_CloseF failed." )
    ENDIF          

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BDIAG file: '//l1BFileInfo%DiagFileName)

!! Close the HDF 5 Fortran Interface

    IF (HDFversion == HDFVERSION_5) THEN
       CALL h5close_f (error)
       IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
            "Fortran HDF 5 API error on closing.")
    ENDIF
       
  END SUBROUTINE CloseFiles

!=============================================================================
END MODULE Close_files
!=============================================================================
! $Log$
! Revision 2.8  2002/11/20 17:50:27  perun
! Reinstated writing PCF annotations for HDF 4 files
!
! Revision 2.7  2002/11/19 21:33:23  perun
! Use HDFversion instead of HDFVersionString
!
! Revision 2.6  2002/11/07 21:54:54  jdone
! Added HDF4/HDF5 switch.
!
! Revision 2.5  2002/03/29 20:18:34  perun
! Version 1.0 commit
!
! Revision 2.4  2001/03/12 19:40:11  perun
! Write PCF and CF contents as annotations
!
! Revision 2.2  2001/03/05 22:35:13  perun
! Corrected L0_Close call
!
! Revision 2.1  2001/02/23 18:51:39  perun
! Version 0.5 commit
!
