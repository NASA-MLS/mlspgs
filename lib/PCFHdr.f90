
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE PCFHdr
!===============================================================================

   USE Hdf
   USE MLSCommon
   USE MLSMessageModule
   USE SDPToolkit
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

! Contents:

! Subroutines -- CreatePCFAnnotation
!                WritePCF2Hdr

! Remarks:  This module contains subroutines for writing the PCF as an annotation
! to HDF files.

CONTAINS

!--------------------------------------------------------
   SUBROUTINE CreatePCFAnnotation (mlspcfN_pcf_start, anText)
!--------------------------------------------------------

! Brief description of subroutine
! This subroutine stores the PCF as an annotation for writing to file headers.

! Arguments

      CHARACTER (LEN=1), POINTER :: anText(:)

      INTEGER :: mlspcfN_pcf_start

! Parameters

! Functions

      INTEGER, EXTERNAL :: Pgs_pc_getFileSize

! Variables

      CHARACTER (LEN=10) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: err, ios, pcfHandle, returnStatus, size, version

! Get the size of the PCF

      version = 1
      returnStatus = Pgs_pc_getFileSize(mlspcfN_pcf_start, version, size)

      ALLOCATE(anText(size), STAT=err)
      IF ( err /= 0 ) THEN
         msr = MLSMSG_Allocate // ' anText PCF array.'
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Open the PCF for reading

      version = 1
      returnStatus = Pgs_io_gen_openF (mlspcfN_pcf_start, PGSd_IO_Gen_RDirUnf, &
                                       size, pcfHandle, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Read the PCF text into the CHAR anText variable

      READ(UNIT=pcfHandle, REC=1, IOSTAT=ios) anText

! Close the PCF

      returnStatus = Pgs_io_gen_closeF (pcfHandle)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

!------------------------------------
   END SUBROUTINE CreatePCFAnnotation
!------------------------------------

!----------------------------------------
   SUBROUTINE WritePCF2Hdr (file, anText)
!----------------------------------------

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS file as an annotation.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file 

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: afEnd, afEndAccess, afFCreate, afStart, afWriteAnn
      INTEGER, EXTERNAL :: hClose, hOpen

! Variables

      CHARACTER (LEN=480) :: msr

      INTEGER :: anID, annID, fileID, status

! Open the HDF-EOS file for writing

      fileID = hOpen(file, DFACC_WRITE, 0)
      IF (fileID == -1) THEN
         msr = MLSMSG_Fileopen // file
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      ENDIF

! Initialize the AN interface

      anID = afStart(fileID)
      IF (anID == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                              &initialize the AN interface.')

! Create a file annotation

      annID = afFCreate(anID, AN_FILE_DESC)
      IF (annID == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                               &create the file annotation.')

! Write the PCF as an annotation to the file

      status = afWriteAnn( annID, anText, SIZE(anText) )

! Terminate access to the annotation

      status = afEndAccess(annID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                    &terminate access to the file annotation.')

! Terminate access to the AN interface

      status = afEnd(anID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                       &terminate access to the AN interface.')

! Close the HDF file

      status = hClose(fileID)
      IF (status == -1) CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
                                                       &close the HDF file.')

!-----------------------------
   END SUBROUTINE WritePCF2Hdr
!-----------------------------

!================
END MODULE PCFHdr
!================

!# $Log$
!# Revision 2.1  2001/03/09 21:10:32  nakamura
!# Routines for writing the PCF to an HDF file as an annotation.
!#
!#
