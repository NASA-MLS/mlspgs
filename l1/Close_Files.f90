! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE Close_files ! Close the production files
!=============================================================================

  USE MLSL1Common, ONLY: L0FileInfo, L1BFileInfo, L1ProgType, THzType, &
       LogType, HDFversion
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

    USE SDPToolkit, ONLY: PGS_IO_Gen_closeF
    USE MLSFiles, ONLY: MLS_closeFile
    USE L1BOutUtils, ONLY: WriteHdrAnnots
    USE HDF5, ONLY: H5gClose_f, H5gOpen_f
    USE MLSHDF5, ONLY: MakeHDF5Attribute, mls_h5close
    USE Orbit, ONLY: OrbitNumber, OrbPeriod
!    USE H5LIB

    CHARACTER(LEN=132) :: filename
    INTEGER :: i, returnStatus, error, grp_id
    LOGICAL :: opened

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

! Close Eng/Sci MAF data unit:

    INQUIRE (unit=L1BFileInfo%MAF_data_unit, name=filename, opened=opened)
    IF (opened) THEN
       CLOSE (L1BFileInfo%MAF_data_unit)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed MAF_data file: '//filename)
    ENDIF

    IF (L1ProgType == LogType) THEN

       returnStatus = PGS_IO_Gen_CloseF (l1BFileInfo%LogId)

       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BLOG file: '//l1BFileInfo%LogFileName)

       INQUIRE (unit=L1BFileInfo%engMAF_unit, name=filename)
       CLOSE (L1BFileInfo%engMAF_unit, status="DELETE")
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed engMAF file: '//filename)

       INQUIRE (unit=L1BFileInfo%sciMAF_unit, name=filename)
       CLOSE (L1BFileInfo%sciMAF_unit, status="DELETE")
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed sciMAF file: '//filename)

       RETURN

    ENDIF

! Close HDF files

    IF (L1ProgType == THzType) THEN

       ! Write Hdr Annotations and Close L1RAD T file

       CALL WriteHdrAnnots (L1BFileInfo%RADTFileName, HDFversion)
       CALL MLS_closeFile (L1BFileInfo%RADTid, HDFversion=HDFversion)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BRAD T file: '//L1BFileInfo%RADTFileName)

       CALL WriteMetaData (IsTHz=.TRUE.)

       RETURN

    ENDIF

    ! Write Hdr Annotations and Close L1RAD D file

    CALL WriteHdrAnnots (L1BFileInfo%RADDFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%RADDid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)

    ! Write Hdr Annotations and Close L1RAD F file

    CALL WriteHdrAnnots (L1BFileInfo%RADGFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%RADGid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD G file: '//L1BFileInfo%RADGFileName)

    ! Write Hdr Annotations and Close L1BOA file

    CALL H5gOpen_f (L1BFileInfo%OAid, '/', grp_id, returnStatus)
    CALL MakeHDF5Attribute (grp_id, 'OrbitNumber', OrbitNumber, .true.)
    CALL MakeHDF5Attribute (grp_id, 'OrbitPeriod', OrbPeriod, .true.)
    CALL H5gClose_f (grp_id, returnStatus)

    CALL WriteHdrAnnots (L1BFileInfo%OAFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%OAid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

    ! Write Hdr Annotations and Close L1B Diag file

    CALL WriteHdrAnnots (L1BFileInfo%DiagFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%DiagId, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1B Diag file: '//L1BFileInfo%DiagFileName)

    CALL WriteMetaData

    INQUIRE (unit=L1BFileInfo%EngId, name=filename, opened=opened)
    IF (opened) THEN
       CLOSE (L1BFileInfo%EngId)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BENG file: '//filename)
    ENDIF

!! Close the HDF 5 Fortran Interface

    CALL mls_h5close (error)
    IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         "Fortran HDF 5 API error on closing.")
    
  END SUBROUTINE CloseFiles

!=============================================================================
END MODULE Close_files
!=============================================================================
! $Log$
! Revision 2.14  2004/05/06 21:59:23  pwagner
! Uses mls_h5open/close
!
! Revision 2.13  2004/01/09 17:46:22  perun
! Version 1.4 commit
!
! Revision 2.12  2003/09/15 17:15:53  perun
! Version 1.3 commit
!
! Revision 2.11  2003/08/15 14:25:04  perun
! Version 1.2 commit
!
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
