! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Close_files ! Close the production files
!=============================================================================

  USE MLSL1Common, ONLY: L0FileInfo, L1BFileInfo, L1ProgType, THzType, LogType
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info
  USE WriteMetaL1

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CloseFiles

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CloseFiles   ! Close production Level 1 files
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE SDPToolkit, ONLY: PGS_IO_Gen_closeF
    USE MLSFiles, ONLY: MLS_closeFile, HDFVERSION_5
    USE L1BOutUtils, ONLY: WriteHdrAnnots

    INTEGER :: i, returnStatus, error, HDFversion

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

    IF (L1ProgType == LogType) THEN

       returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%LogId)

       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BLOG file: '//l1BFileInfo%LogFileName)

       RETURN

    ENDIF

! Close HDF files

    HDFversion = L1Config%Output%HDFversion

    IF (L1ProgType == THzType) THEN

       ! Write Hdr Annotations and Close L1RAD T file

       CALL WriteHdrAnnots (L1BFileInfo%RADTFileName, L1BFileInfo%RADTid, &
            HDFversion)
       CALL MLS_closeFile (L1BFileInfo%RADTid, HDFversion=HDFversion)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BRAD T file: '//L1BFileInfo%RADTFileName)

       CALL WriteMetaData (IsTHz=.TRUE.)

       RETURN

    ENDIF

    ! Write Hdr Annotations and Close L1RAD D file

    CALL WriteHdrAnnots (L1BFileInfo%RADDFileName, L1BFileInfo%RADDid, &
         HDFversion)
    CALL MLS_closeFile (L1BFileInfo%RADDid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)

    ! Write Hdr Annotations and Close L1RAD F file

    CALL WriteHdrAnnots (L1BFileInfo%RADFFileName, L1BFileInfo%RADFid, &
         HDFversion)
    CALL MLS_closeFile (L1BFileInfo%RADFid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD F file: '//L1BFileInfo%RADFFileName)

    ! Write Hdr Annotations and Close L1BOA file

    CALL WriteHdrAnnots (L1BFileInfo%OAFileName, L1BFileInfo%OAid, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%OAid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

    CALL WriteMetaData

    returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%EngId)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BENG file: '//l1BFileInfo%EngFileName)

    returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%DiagId)

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
! Revision 2.10  2003/02/28 19:11:35  pwagner
! In accord with new MLS_closeFile interface
!
! Revision 2.9  2003/01/31 18:13:34  perun
! Version 1.1 commit
!
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
