! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!=============================================================================
MODULE CloseBOA ! Close the BOA file
!=============================================================================

  USE MLSL1Common, ONLY: L1BFileInfo, HDFversion
  USE MLSMessageModule, ONLY: MLSMessage, MLSMSG_Error, MLSMSG_Info

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: CloseBOAfile

  !------------------------------- RCS Ident Info ------------------------------
  CHARACTER(LEN=130) :: id = &
       "$Id$"
  CHARACTER(LEN=*), PARAMETER :: ModuleName="$RCSfile$"
  !-----------------------------------------------------------------------------

CONTAINS

!=============================================================================
  SUBROUTINE CloseBOAfile   ! Close BOA file
!=============================================================================

    USE MLSL1Config, ONLY: L1Config
    USE MLSFiles, ONLY: MLS_closeFile
    USE L1BOAutils, ONLY: WriteHdrAnnots
    USE Orbit, ONLY: OrbitNumber, OrbPeriod
    USE HDF5, ONLY: H5gClose_f, H5gOpen_f
    USE MLSHDF5, ONLY: MakeHDF5Attribute, MLS_h5close

    INTEGER :: returnStatus, error, grp_id

! Close HDF file

    ! Write Hdr Annotations and Close L1BOA file

    CALL H5gOpen_f (L1BFileInfo%OAid, '/', grp_id, returnStatus)
    CALL MakeHDF5Attribute (grp_id, 'OrbitNumber', OrbitNumber, .true.)
    CALL MakeHDF5Attribute (grp_id, 'OrbitPeriod', OrbPeriod, .true.)
    CALL H5gClose_f (grp_id, returnStatus)

    CALL WriteHdrAnnots (L1BFileInfo%OAFileName, HDFversion)
    CALL MLS_closeFile (L1BFileInfo%OAid, HDFversion=HDFversion)
    CALL MLSMessage (MLSMSG_Info, ModuleName, &
         & 'Closed L1BOA file: '//L1BFileInfo%OAFileName)

!! Close the HDF 5 Fortran Interface

    CALL MLS_h5close (error)
    IF (error /= 0) CALL MLSMessage (MLSMSG_Error, ModuleName, &
         "Fortran HDF 5 API error on closing.")
    
  END SUBROUTINE CloseBOAfile

!=============================================================================
END MODULE CloseBOA
!=============================================================================
! $Log$
! Revision 2.4  2004/05/14 15:59:11  perun
! Version 1.43 commit
!
! Revision 2.3  2004/05/06 21:59:23  pwagner
! Uses mls_h5open/close
!
! Revision 2.2  2004/01/09 20:02:57  perun
! Update BOA to HDF 5
!
! Revision 2.1  2003/10/24 19:38:36  perun
! Version 1.3 commit
!
!
