
! Copyright (c) 2001, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE PCFHdr
!===============================================================================

   USE Hdf
  use LEXER_CORE, only: PRINT_SOURCE
   USE MLSCommon
  USE output_m, only: output
   USE SDPToolkit
  use TREE, only: DUMP_TREE_NODE, SOURCE_REF
   IMPLICIT NONE
   PUBLIC

   PRIVATE :: ID, ModuleName

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName="$RCSfile$"
!----------------------------------------------------------

  integer, private :: ERROR

! Contents:

! Subroutines -- CreatePCFAnnotation
!                WritePCF2Hdr

! Remarks:  This module contains subroutines for writing the PCF as an annotation
! to HDF files.

CONTAINS

!------------------------------------------------------------
   SUBROUTINE CreatePCFAnnotation (mlspcfN_pcf_start, anText)
!------------------------------------------------------------

! Brief description of subroutine
! This subroutine stores the PCF as an annotation for writing to file headers.

! Arguments

      INTEGER, INTENT(IN) :: mlspcfN_pcf_start

      CHARACTER (LEN=1), POINTER :: anText(:)

! Parameters

! Functions

      INTEGER, EXTERNAL :: Pgs_pc_getFileSize

! Variables

      CHARACTER (LEN=10) :: mnemonic
      CHARACTER (LEN=480) :: msg, msr

      INTEGER :: err, ios, pcfHandle, returnStatus, size, version
		integer, parameter :: DEFAULTANTEXTSIZE=1024
		character(len=*), parameter :: DEFAULTANTEXT= &
		& 'PCF file not found--check it has the right PCF number(900)'

! Get the size of the PCF

    error = 0
      version = 1
      returnStatus = Pgs_pc_getFileSize(mlspcfN_pcf_start, version, size)

      IF (returnStatus /= PGS_S_SUCCESS) THEN
!         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
!         msr = mnemonic // ':  ' // msg
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
			call announce_error(0, &
			& 'Failed to find PCF file--check its PCF number')
			size = DEFAULTANTEXTSIZE
      	ALLOCATE(anText(size), STAT=err)
			anText(:len(DEFAULTANTEXT)) = DEFAULTANTEXT(:len(DEFAULTANTEXT))
			return
      ENDIF

      ALLOCATE(anText(size), STAT=err)
      IF ( err /= 0 ) THEN
 !        msr = MLSMSG_Allocate // ' anText PCF array.'
 !        CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
			call announce_error(0, &
			& 'Failed to allocate anText for storing PCF contents')
      ENDIF

! Open the PCF for reading

      version = 1
      returnStatus = Pgs_io_gen_openF (mlspcfN_pcf_start, PGSd_IO_Gen_RDirUnf, &
                                       size, pcfHandle, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
!         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
!         msr = mnemonic // ':  ' // msg
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
			call announce_error(0, &
			& 'Failed to open PCF file--check its PCF number or name')
      ENDIF

! Read the PCF text into the CHAR anText variable

      READ(UNIT=pcfHandle, REC=1, IOSTAT=ios) anText

! Close the PCF

      returnStatus = Pgs_io_gen_closeF (pcfHandle)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
!         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
!         msr = mnemonic // ':  ' // msg
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
			call announce_error(0, &
			& 'Failed to close PCF file')
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

    error = 0

! Open the HDF-EOS file for writing

      fileID = hOpen(file, DFACC_WRITE, 0)
      IF (fileID == -1) THEN
!         msr = MLSMSG_Fileopen // file
!         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
			call announce_error(0, &
			& 'Failed to open hdf-eos file: ' // file)
      ENDIF

! Initialize the AN interface

      anID = afStart(fileID)
      IF (anID == -1) then
!			CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
!                                              &initialize the AN interface.')
			call announce_error(0, &
			& 'Failed to initialize AN interface')
		endif
! Create a file annotation

      annID = afFCreate(anID, AN_FILE_DESC)
      IF (annID == -1) then
!			 CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
!                                               &create the file annotation.')
			call announce_error(0, &
			& 'Failed to create file annotation')
		endif
! Write the PCF as an annotation to the file

      status = afWriteAnn( annID, anText, SIZE(anText) )

! Terminate access to the annotation

      status = afEndAccess(annID)
      IF (status == -1) then
!			CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
!                                    &terminate access to the file annotation.')
			call announce_error(0, &
			& 'Failed to terminate access to annotation')
		endif
! Terminate access to the AN interface

      status = afEnd(anID)
      IF (status == -1) then
!			CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
!                                       &terminate access to the AN interface.')
			call announce_error(0, &
			& 'Failed to find terminate access to AN interface')
		endif
! Close the HDF file

      status = hClose(fileID)
      IF (status == -1) then
!			CALL MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
!                                                       &close the HDF file.')
			call announce_error(0, &
			& 'Failed to close HDF file')
		endif
!-----------------------------
   END SUBROUTINE WritePCF2Hdr
!-----------------------------

  ! ------------------------------------------------  announce_error  -----
  subroutine announce_error ( lcf_where, full_message, use_toolkit, &
  & error_number )
  
   ! Arguments
	
	integer, intent(in)    :: lcf_where
	character(LEN=*), intent(in)    :: full_message
	logical, intent(in), optional :: use_toolkit
	integer, intent(in), optional    :: error_number

	! Local
  logical :: just_print_it
  logical, parameter :: default_output_by_toolkit = .true.
	
	if(present(use_toolkit)) then
		just_print_it = use_toolkit
	elseif(default_output_by_toolkit) then
		just_print_it = .false.
	else
		just_print_it = .true.
	endif
	
	if(.not. just_print_it) then
    error = max(error,1)
    call output ( '***** At ' )

	if(lcf_where > 0) then
	    call print_source ( source_ref(lcf_where) )
		else
    call output ( '(no lcf node available)' )
		endif

    call output ( ': ' )
    call output ( "The " );
	if(lcf_where > 0) then
    call dump_tree_node ( lcf_where, 0 )
		else
    call output ( '(no lcf tree available)' )
		endif

		CALL output("Caused the following error:", advance='yes', &
		& from_where=ModuleName)
		CALL output(trim(full_message), advance='yes', &
		& from_where=ModuleName)
		if(present(error_number)) then
			CALL output('error number ', advance='no')
			CALL output(error_number, places=9, advance='yes')
		endif
	else
		print*, '***Error in module ', ModuleName
		print*, trim(full_message)
		if(present(error_number)) then
			print*, 'error number ', error_number
		endif
	endif

!===========================
  end subroutine announce_error
!===========================

!================
END MODULE PCFHdr
!================

!# $Log$
!# Revision 2.4  2001/04/04 22:23:03  pwagner
!# Added announce_error; attempt recovery if PCF file not found
!#
!# Revision 2.3  2001/04/04 19:37:25  vsnyder
!# Try to produce an error message if there's no entry for the PCF in the PCF
!#
!# Revision 2.2  2001/03/09 21:32:45  nakamura
!# Added INTENT(IN) for pcf number arg.
!#
!# Revision 2.1  2001/03/09 21:10:32  nakamura
!# Routines for writing the PCF to an HDF file as an annotation.
!#
!#
