! Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
MODULE PCFHdr
!===============================================================================
! Unfortunately from the point of view of its name, which fails to clearly 
! indicate it, this module contains ways to annotate, write
! attributes, and otherwise do miscellaneous things to mls product files.
! It might have been better named PCFHdrAndGlobalAttributes
! or split off global attribute stuff into a separate module
   USE Hdf, only: DFACC_WRITE, AN_FILE_DESC
!  USE Hdf, only: hOpen, afStart, afFCreate, afWriteAnn, afEndAccess, &
!    & afEnd, hClose
   USE MLSCommon, only: i4, r4, r8, FileNameLen
   USE MLSFiles, only: GetPCFromRef, HDFVERSION_4, HDFVERSION_5
   USE MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
     & MLSMSG_FILEOPEN
   USE SDPToolkit, only: PGSD_PC_UREF_LENGTH_MAX, PGS_S_SUCCESS, &
     & PGSD_MET_GROUP_NAME_L, PGS_IO_GEN_CLOSEF, PGS_IO_GEN_OPENF, &
     & PGSD_IO_GEN_RDIRUNF
   IMPLICIT NONE
   PUBLIC :: GlobalAttributes_T, &
     & CreatePCFAnnotation, gd_writeglobalattr, h5_writeglobalattr, &
     & InputInputPointer, &
     & sw_writeglobalattr, WritePCF2Hdr, WriteInputPointer
   private

  !---------------------------- RCS Ident Info -------------------------------
  character(len=*), private, parameter :: IdParm = &
    & "$Id$"
  character(len=len(idparm)), private :: Id = idParm
  character(len=*), private, parameter :: ModuleName = &
       & "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

! Contents:

! Subroutines -- CreatePCFAnnotation
!                WritePCF2Hdr
!                WriteInputPointer
!                InputInputPointer
!                gd_writeglobalattr
!                sw_writeglobalattr

! Remarks:  This module contains subroutines for writing the PCF as an 
! annotation to HDF files. (obsolete)
! It also contains the two routines that prepare and write the input pointer
! to metadata
! and the routines for writing attributes to the various files requiring them
! in particular swath, grid, and plain hdf5
  integer, parameter, public :: INPUTPTR_STRING_LENGTH = PGSd_PC_UREF_LENGTH_MAX
  integer, parameter, public :: GA_VALUE_LENGTH = 40

   ! May get some of these from MLSLibOptions? 
  type GlobalAttributes_T
    character(len=GA_VALUE_LENGTH) :: InstrumentName = 'MLS Aura'
    character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''
    character(len=GA_VALUE_LENGTH) :: InputVersion = ''
    character(len=GA_VALUE_LENGTH) :: PGEVersion = ''
    character(len=GA_VALUE_LENGTH) :: StartUTC = ''
    character(len=GA_VALUE_LENGTH) :: EndUTC = ''
    integer :: GranuleMonth                  = 0
    integer :: GranuleDay                    = 0
    integer :: GranuleYear                   = 0
  end type GlobalAttributes_T

  ! This variable describes the global attributes
  type (GlobalAttributes_T), public, save :: GlobalAttributes
  ! Use this in case hdfVersion omitted from call to WritePCF2Hdr
  ! E.g., in level 3 prior to conversion
  integer, public, save            :: PCFHDR_DEFAULT_HDFVERSION = HDFVERSION_4

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

!------------------------------------------------------------
   SUBROUTINE gd_writeglobalattr (gridID)
!------------------------------------------------------------

!     use HDF5, only: H5T_NATIVE_CHARACTER
      use HDFEOS5, only: HE5T_NATIVE_SCHAR
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf-eos5 grid

! Arguments

      INTEGER, INTENT(IN) :: gridID
      integer, external ::   he5_GDwrattr
! Internal variables
      integer :: status
! Executable
      status = he5_GDwrattr(gridID, &
       & 'InstrumentName', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InstrumentName)
      status = he5_GDwrattr(gridID, &
       & 'ProcessLevel', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%ProcessLevel)
      status = he5_GDwrattr(gridID, &
       & 'InputVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InputVersion)
      status = he5_GDwrattr(gridID, &
       & 'PGEVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%PGEVersion)
      status = he5_GDwrattr(gridID, &
       & 'StartUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%StartUTC)
      status = he5_GDwrattr(gridID, &
       & 'EndUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%EndUTC)
!------------------------------------------------------------
   END SUBROUTINE gd_writeglobalattr
!------------------------------------------------------------

!------------------------------------------------------------
   SUBROUTINE h5_writeglobalattr (fileID)
!------------------------------------------------------------

      use HDF5, only: h5gclose_f, h5gopen_f
      USE MLSHDF5, only: MakeHDF5Attribute
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf5-formatted file
! It does so at the root '/' group level

! Arguments

      INTEGER, INTENT(IN) :: fileID
! Local variables
      integer :: grp_id
      integer :: status
! Executable
      call h5gopen_f(fileID, '/', grp_id, status)
      call MakeHDF5Attribute(grp_id, &
       & 'InstrumentName', GlobalAttributes%InstrumentName, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'ProcessLevel', GlobalAttributes%ProcessLevel, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'InputVersion', GlobalAttributes%InputVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'PGEVersion', GlobalAttributes%PGEVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'StartUTC', GlobalAttributes%StartUTC, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'EndUTC', GlobalAttributes%EndUTC, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleMonth', GlobalAttributes%GranuleMonth, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleDay', GlobalAttributes%GranuleDay, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleYear', GlobalAttributes%GranuleYear, .true.)
      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   END SUBROUTINE h5_writeglobalattr
!------------------------------------------------------------

!----------------------------------------
   SUBROUTINE InputInputPointer (urefs, fileIDArray, fileNameArray, &
     & PCBottom, PCTop )
!----------------------------------------
!  Prepare Input for WriteInputPointer consisting of universal refs
!  This can be done for an array of fileids or of file names or both
!  that were input during the run; e.g., l1brads for level 2
!  Arguments
    CHARACTER (LEN=INPUTPTR_STRING_LENGTH), intent(out)    :: urefs(:)
    CHARACTER (LEN=*), dimension(:), intent(in), optional  :: fileNameArray
    integer, dimension(:), intent(in), optional            :: fileIDArray
    integer,  intent(IN), optional                 :: PCBottom, PCTop
!  Local variables
   integer  :: i
   integer  :: ref
   integer  :: returnStatus
   integer  :: thePC
   integer  :: version
   character(len=INPUTPTR_STRING_LENGTH) :: sval
   INTEGER, EXTERNAL :: pgs_pc_getUniversalRef

 ! Executable
   if (size(urefs) < 1) return
   urefs = ' '
   ref = 0
   if ( size(fileIDArray) > 0 ) then
      DO i = 0, size(fileIDArray)
        version = 1
        if ( fileIDArray(i) > 0 ) then
         returnStatus = pgs_pc_getUniversalRef(fileIDArray(i), &
           & version, sval)
        else
         returnStatus = PGS_S_SUCCESS + 1
        endif
        IF (returnStatus == PGS_S_SUCCESS) THEN 
           ref = min(ref+1, size(urefs))
           urefs(ref) = sval                     
        ENDIF                                   
      ENDDO
   endif

   if ( size(fileNameArray) > 0 ) then
      DO i = 0, size(fileNameArray)
        version = 1
        thePC = GetPCFromRef(trim(fileNameArray(i)), PCBottom, PCTop, &
          & .true., returnStatus, version)
        if ( thePC > 0 ) then
         version = 1
         returnStatus = pgs_pc_getUniversalRef(thePC, &
           & version, sval)
        else
         returnStatus = PGS_S_SUCCESS + 1
        endif
        IF (returnStatus == PGS_S_SUCCESS) THEN 
           ref = min(ref+1, size(urefs))
           urefs(ref) = sval                     
        ENDIF                                   
      ENDDO
   endif

!------------------------------------
   END SUBROUTINE InputInputPointer
!------------------------------------

!------------------------------------------------------------
   SUBROUTINE sw_writeglobalattr (swathID)
!------------------------------------------------------------

!     use HDF5, only: H5T_NATIVE_CHARACTER
      use HDFEOS5, only: HE5T_NATIVE_SCHAR
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf-eos5 swath

! Arguments

      INTEGER, INTENT(IN) :: swathID
      integer, external ::   he5_SWwrattr
! Internal variables
      integer :: status
! Executable
      status = he5_SWwrattr(swathID, &
       & 'InstrumentName', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InstrumentName)
      status = he5_SWwrattr(swathID, &
       & 'ProcessLevel', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%ProcessLevel)
      status = he5_SWwrattr(swathID, &
       & 'InputVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InputVersion)
      status = he5_SWwrattr(swathID, &
       & 'PGEVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%PGEVersion)
      status = he5_SWwrattr(swathID, &
       & 'StartUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%StartUTC)
      status = he5_SWwrattr(swathID, &
       & 'EndUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%EndUTC)
!------------------------------------------------------------
   END SUBROUTINE sw_writeglobalattr
!------------------------------------------------------------

!----------------------------------------
   FUNCTION WriteInputPointer (groups, attrName, inpt)
!----------------------------------------

!  Write InputPointer metadata
!  Moved here to hide inconsistency of arguments from NAGging inquiry
!  in the style of Enron's offshore limited partnerships

!  Arguments
    character (len = PGSd_MET_GROUP_NAME_L) :: Groups
    character (len=132), intent(in) :: Attrname
    CHARACTER (LEN=INPUTPTR_STRING_LENGTH), intent(in)  :: inpt(:)

    integer             :: WriteInputPointer
    integer, external   :: pgs_met_setAttr_s

!   Executable statements
       WriteInputPointer = pgs_met_setAttr_s(groups, attrName, inpt)

!------------------------------------
   END FUNCTION WriteInputPointer
!------------------------------------

!----------------------------------------
   SUBROUTINE WritePCF2Hdr (file, anText, hdfVersion)
!----------------------------------------
      use HDF5, only: H5F_ACC_RDWR_F, &
        & h5fopen_f, h5fclose_f

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS file as an annotation.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file 

      CHARACTER (LEN=1), POINTER              :: anText(:)
      integer, intent(in), optional           :: hdfVersion

! Parameters
      integer :: my_hdfVersion
      integer :: fileID, status
! Executable
      my_hdfVersion = PCFHDR_DEFAULT_HDFVERSION
      if ( present(hdfVersion) ) my_hdfVersion = hdfVersion
      select case(my_hdfVersion)
      case (HDFVERSION_4)
        call WritePCF2Hdr_hdf4 (file, anText)
      case (HDFVERSION_5)
        call h5fopen_f(trim(file), H5F_ACC_RDWR_F, fileID, status)
        if ( status /= PGS_S_SUCCESS) &
          & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Error opening hdf5 file for annotating with PCF' )
!       call WritePCF2Hdr_hdf5 (file, anText)
        call WritePCF2Hdr_hdf5 (fileID, anText)
        call h5fclose_f(fileID, status)
        if ( status /= PGS_S_SUCCESS) &
          & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Error closing hdf5 file for annotating with PCF' )
      case default
         CALL MLSMessage(MLSMSG_Error, ModuleName, &
           & 'Unrecognized hdfVersion in WritePCF2Hdr')
      end select
!-----------------------------
   END SUBROUTINE WritePCF2Hdr
!-----------------------------

!----------------------------------------
   SUBROUTINE WritePCF2Hdr_hdf4 (file, anText)
!----------------------------------------

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS4 file as an annotation.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file 

      CHARACTER (LEN=1), POINTER              :: anText(:)

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
   END SUBROUTINE WritePCF2Hdr_hdf4
!-----------------------------

!----------------------------------------
!  SUBROUTINE WritePCF2Hdr_hdf5 (file, anText)
   SUBROUTINE WritePCF2Hdr_hdf5 (fileID, anText)
!----------------------------------------

!     use HDF5, only: H5F_ACC_RDWR_F, &
!       & h5fopen_f, h5fclose_f, h5gclose_f, h5gopen_f
      use HDF5, only: h5gclose_f, h5gopen_f
      USE MLSHDF5, only: MakeHDF5Attribute
! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS5 file as an attribute.
! It does so at the root '/' group level, treating the file as if
! it were a plain hdf5 file
! This may result in the dread hybrid file syndrome which may
! confuse some hdf-eos5 readers (not tested)

! Arguments

!      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file 

      CHARACTER (LEN=1), POINTER              :: anText(:)

! Local variables
      integer :: fileID
      integer :: grp_id
      integer :: status
! Executable
!      call h5fopen_f(trim(file), H5F_ACC_RDWR_F, fileID, status)
      ! if ( status /= PGS_S_SUCCESS) &
      !  & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      !  & 'Error opening hdf5 file for annotating with PCF' )
      call h5gopen_f(fileID, '/', grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error opening hdf5 file root group for annotating with PCF' )
      call MakeHDF5Attribute(grp_id, &
       & 'PCF file text', anText, .true.)
      call h5gclose_f(grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error closing hdf5 file root group for annotating with PCF' )
      ! call h5fclose_f(fileID, status)
      ! if ( status /= PGS_S_SUCCESS) &
      !  & CALL MLSMessage(MLSMSG_Error, ModuleName, &
      !  & 'Error closing hdf5 file for annotating with PCF' )
!-----------------------------
   END SUBROUTINE WritePCF2Hdr_hdf5
!-----------------------------

!================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PCFHdr
!================

!# $Log$
!# Revision 2.9  2003/01/30 00:57:24  pwagner
!# Added much new stuff for global attributes with hdf5
!#
!# Revision 2.8  2002/10/08 00:09:13  pwagner
!# Added idents to survive zealous Lahey optimizer
!#
!# Revision 2.7  2002/10/01 22:03:54  pwagner
!# Fixed RCS Ident Block
!#
!# Revision 2.6  2002/08/29 16:54:44  pwagner
!# Added WriteInputPointer
!#
!# Revision 2.5  2001/04/06 16:54:57  pwagner
!# Reverting to version 2.2
!#
!# Revision 2.2  2001/03/09 21:32:45  nakamura
!# Added INTENT(IN) for pcf number arg.
!#
!# Revision 2.1  2001/03/09 21:10:32  nakamura
!# Routines for writing the PCF to an HDF file as an annotation.
!#
!#
