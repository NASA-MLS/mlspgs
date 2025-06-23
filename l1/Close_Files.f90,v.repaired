! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

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

!---------------------------- RCS Module Info ------------------------------
  CHARACTER (len=*), PRIVATE, PARAMETER :: ModuleName= &
       "$RCSfile$"
  PRIVATE :: not_used_here 
!---------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CloseFiles   ! Close production Level 1 files
!=============================================================================

    USE SDPToolkit, ONLY: PGS_IO_Gen_closeF
    USE MLSFiles, ONLY: MLS_closeFile
    USE L1BOutUtils, ONLY: WriteHdrAnnots
    USE OutputL1B, ONLY: OutputL1B_diags, OutputL1B_Chi2
    USE HDF5, ONLY: H5gClose_f, H5gOpen_f
    USE MLSHDF5, ONLY: MakeHDF5Attribute, MLS_h5close
    USE Orbit, ONLY: OrbitNumber, OrbPeriod
    USE PCFHdr, ONLY: h5_writeglobalattr, CreatePCFAnnotation, WritePCF2Hdr
    USE BrightObjects_m, ONLY: BO_name, BO_Angle_GHz, BO_Angle_THz
    USE DACsUtils, ONLY: TPz
    USE INTRINSIC, ONLY: l_hdf
    USE MLSL1Common, ONLY: FileNameLen

    CHARACTER (LEN=FileNameLen) :: filename, msg


    INTEGER :: i, returnStatus, error, grp_id
    LOGICAL :: exist
    LOGICAL :: opened

    CHARACTER (LEN=1), POINTER :: leapsecText(:), utcpoleText(:)
    INTEGER, PARAMETER :: leapsec_pcf = 10301, utcpole_pcf = 10401
    logical, parameter :: DeeBug = .false.

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
       rewind ( L1BFileInfo%engMAF_unit )
       INQUIRE (unit=L1BFileInfo%engMAF_unit, name=filename, exist=exist, opened=opened)
       if ( DeeBug ) then
         print *, 'filename: ', trim(filename)
         print *, 'exist: ', exist
         print *, 'opened: ', opened
         print *, 'engMAF_unit: ', L1BFileInfo%engMAF_unit
       endif
       CLOSE (L1BFileInfo%engMAF_unit, status="KEEP", &
         & iostat=returnStatus, iomsg=msg )
       if ( DeeBug ) then
         print *, 'msg: ', trim(msg)
         print *, 'status: ', returnStatus
       endif
       OPEN (unit=L1BFileInfo%EngMAF_unit, file=Filename, &
         status="REPLACE", FORM="UNFORMATTED", ACCESS="SEQUENTIAL", iostat=returnStatus)
       INQUIRE (unit=L1BFileInfo%engMAF_unit, name=filename, exist=exist, opened=opened)
       if ( DeeBug ) then
         print *, '(after reopening) filename: ', trim(filename)
         print *, 'status: ', returnStatus		
         print *, 'exist: ', exist
         print *, 'opened: ', opened
         print *, 'engMAF_unit: ', L1BFileInfo%engMAF_unit
       endif
       CLOSE (L1BFileInfo%engMAF_unit, status="DELETE", &
         & iostat=returnStatus, iomsg=msg )
       if ( DeeBug ) then
         print *, 'msg: ', trim(msg)
         print *, 'status: ', returnStatus
       endif
       ! stop
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed engMAF file: '//filename)

       INQUIRE (unit=L1BFileInfo%sciMAF_unit, name=filename)
       CLOSE (L1BFileInfo%sciMAF_unit, status="DELETE")
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed sciMAF file: '//filename)

    ! Write TPz Annotations and Close L1BRADD file.  

       ! <whd> I can't find where the L1BRADD file is opened when running
       ! MLSL1log, so I think that, for MLSL1log, L1BFileInfo%RADDid will
       ! *always* == 0, and I can't find TPz in the output in the L1BRADD files
       ! that I've run.  </whd>

       IF (L1BFileInfo%RADDid /= 0) THEN
          CALL H5gOpen_f (L1BFileInfo%RADDid, '/', grp_id, returnStatus)
          CALL MakeHDF5Attribute (grp_id, 'TPz', TPz, .TRUE.)
          CALL MLS_closeFile (L1BFileInfo%RADDid, HDFversion=HDFversion)
          CALL MLSMessage (MLSMSG_Info, ModuleName, &
               & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)
       ENDIF

       RETURN

    ENDIF

    WRITE (L1BFileInfo%LogId, *) ''
    IF (L1ProgType == THzType) THEN
       WRITE (L1BFileInfo%LogId, *) &
            '##################### End MLSL1T ####################'
    ELSE
       WRITE (L1BFileInfo%LogId, *) &
            '##################### End MLSL1G ####################'
    ENDIF
    WRITE (L1BFileInfo%LogId, *) ''

    CLOSE (l1BFileInfo%LogId)

    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BLOG file: '//l1BFileInfo%LogFileName)

! Close HDF files

    IF (L1ProgType == THzType) THEN

       CALL OutputL1B_Chi2 (L1BFileInfo%RADTid)   ! Write default Chi2

       ! Write Hdr Annotations and Close L1RAD T file

       CALL WriteHdrAnnots (L1BFileInfo%RADTFileName, HDFversion)
       CALL MLS_closeFile (L1BFileInfo%RADTid, HDFversion=HDFversion)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BRAD T file: '//L1BFileInfo%RADTFileName)

       ! Write Hdr Annotations and Close L1B THz Diag file

       CALL WriteHdrAnnots (L1BFileInfo%DiagTFileName, HDFversion)
       CALL h5_writeglobalattr (L1BFileInfo%DiagTid, &
            skip_if_already_there=.FALSE.)
       CALL MLS_closeFile (L1BFileInfo%DiagTId, HDFversion=HDFversion)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1B Diag THz file: '//L1BFileInfo%DiagTFileName)

       CALL WriteMetaData (IsTHz=.TRUE.)

       RETURN

    ENDIF

    IF (L1BFileInfo%RADDid /= 0) THEN
       CALL OutputL1B_Chi2 (L1BFileInfo%RADDid)   ! Write default Chi2

    ! Write Hdr Annotations and Close L1RAD D file

       CALL WriteHdrAnnots (L1BFileInfo%RADDFileName, HDFversion)
       CALL MLS_closeFile (L1BFileInfo%RADDid, HDFversion=HDFversion)
       CALL MLSMessage (MLSMSG_Info, ModuleName, &
            & 'Closed L1BRAD D file: '//L1BFileInfo%RADDFileName)
    ENDIF

    CALL OutputL1B_Chi2 (L1BFileInfo%RADGid)   ! Write default Chi2

    ! Write Hdr Annotations and Close L1RAD G file

    CALL WriteHdrAnnots (L1BFileInfo%RADGFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%RADGid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BRAD G file: '//L1BFileInfo%RADGFileName)

    ! Write Hdr Annotations and Close L1BOA file

    !! Get annotation from leapsec and utcpole files

    CALL CreatePCFAnnotation (leapsec_pcf, leapsecText)

    CALL CreatePCFAnnotation (utcpole_pcf, utcpoleText)

    CALL H5gOpen_f (L1BFileInfo%OAid, '/', grp_id, returnStatus)

    CALL WritePCF2Hdr (L1BFileInfo%OAFileName, leapsecText, &
         hdfVersion=HDFversion, fileType=l_hdf, name='/leapsec')

    CALL WritePCF2Hdr (L1BFileInfo%OAFileName, utcpoleText, &
         hdfVersion=HDFversion, fileType=l_hdf, name='/utcpole')
 
    ! Orbit attributes:

    CALL MakeHDF5Attribute (grp_id, 'OrbitNumber', OrbitNumber, .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'OrbitPeriod', OrbPeriod, .TRUE.)

    ! Bright Object attributes:

    CALL MakeHDF5Attribute (grp_id, 'BO_name', BO_name, .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'BO_Angle_GHz', BO_Angle_GHz, .TRUE.)
    CALL MakeHDF5Attribute (grp_id, 'BO_Angle_THz', BO_Angle_THz, .TRUE.)
    CALL H5gClose_f (grp_id, returnStatus)

    CALL WriteHdrAnnots (L1BFileInfo%OAFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%OAid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

    ! Write default zeros to L1B Diag file

    CALL OutputL1B_diags (L1BFileInfo%Diagid, Zeros=.TRUE.)

    ! Write Hdr Annotations and Close L1B Diag file

    CALL WriteHdrAnnots (L1BFileInfo%DiagFileName, HDFversion)
    CALL h5_writeglobalattr (L1BFileInfo%Diagid, &
         skip_if_already_there=.FALSE.)
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

    CALL MLS_h5close (error)
    IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         "Fortran HDF 5 API error on closing.")
    
  END SUBROUTINE CloseFiles

!=============================================================================
  LOGICAL FUNCTION not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (len=*), PARAMETER :: IdParm = &
       "$Id$"
  CHARACTER (len=LEN(idParm)), SAVE :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
  END FUNCTION not_used_here
END MODULE Close_files
!=============================================================================
! $Log$
! Revision 2.24  2017/01/06 23:42:00  pwagner
! Workaround for bug in NAG 6.1 regarding L1BFileInfo%engMAF_unit
!
! Revision 2.23  2016/03/15 22:17:59  whdaffer
! Merged whd-rel-1-0 back onto main branch. Most changes
! are to comments, but there's some modification to Calibration.f90
! and MLSL1Common to support some new modules: MLSL1Debug and SnoopMLSL1.
!
! Revision 2.22.6.1  2015/10/09 10:21:38  whdaffer
! checkin of continuing work on branch whd-rel-1-0
!
! Revision 2.22  2008/02/25 20:11:38  perun
! Added leapsec and utcpole file contents to L1BOA output file.
!
! Revision 2.21  2007/06/21 20:58:46  perun
! Only write to RADD file if DACS calibration is enabled
!
! Revision 2.20  2006/08/02 18:52:46  perun
! Write TPz Annotations to RADD file
!
! Revision 2.19  2005/12/06 19:22:56  perun
! Output Bright Object attributes to BOA file
!
! Revision 2.18  2005/06/23 18:41:35  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.17  2004/11/10 15:37:51  perun
! Output and close default Chi2 file
!
! Revision 2.16  2004/08/12 13:51:49  perun
! Version 1.44 commit
!
! Revision 2.15  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
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
