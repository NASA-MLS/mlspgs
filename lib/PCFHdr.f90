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
   USE Hdf, only: DFACC_RDWR, DFACC_WRITE, AN_FILE_DESC
!  USE Hdf, only: hOpen, afStart, afFCreate, afWriteAnn, afEndAccess, &
!    & afEnd, hClose
   USE MLSCommon, only: i4, r4, r8, FileNameLen, NameLen
   USE MLSFiles, only: GetPCFromRef, HDFVERSION_4, HDFVERSION_5, &
     & MLS_IO_GEN_OPENF, MLS_IO_GEN_CLOSEF
   USE MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_Error, &
     & MLSMSG_Warning, MLSMSG_DeAllocate, MLSMSG_FILEOPEN
   use MLSStrings, only: utc_to_yyyymmdd, lowerCase
   USE SDPToolkit, only: PGSD_PC_UREF_LENGTH_MAX, PGS_S_SUCCESS, &
     & PGSD_MET_GROUP_NAME_L, PGS_IO_GEN_CLOSEF, PGS_IO_GEN_OPENF, &
     & PGSD_IO_GEN_RDIRUNF, &
     & PGS_TD_ASCIITIME_ATOB, PGS_TD_ASCIITIME_BTOA, &
     & UseSDPToolkit
   IMPLICIT NONE
   PUBLIC :: GlobalAttributes_T, &
     & FillTAI93Attribute, &
     & CreatePCFAnnotation,  &
     & h5_writeglobalattr, he5_writeglobalattr, &
     & InputInputPointer, &
     & WritePCF2Hdr, WriteInputPointer
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
!                h5_writeglobalattr
!                he5_writeglobalattr
!                sw_writeglobalatt
! Remarks:  This module contains subroutines for writing the PCF as an 
! annotation to HDF files. (obsolete)
! It also contains the two routines that prepare and write the input pointer
! to metadata
! and the routines for writing attributes to the various files requiring them
! in particular hdfeos5 and plain hdf5
  integer, parameter, public :: INPUTPTR_STRING_LENGTH = PGSd_PC_UREF_LENGTH_MAX
  integer, parameter, public :: GA_VALUE_LENGTH = 40
  integer, parameter, public :: UTC_A_VALUE_LENGTH = 27
  integer, parameter, public :: UTC_B_VALUE_LENGTH = 25
  character(len=*), parameter, private :: PCFATTRIBUTENAME = 'PCF'

   ! May get some of these from MLSLibOptions? 
  type GlobalAttributes_T
    character(len=GA_VALUE_LENGTH) :: InstrumentName = 'MLS Aura'
    character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''
    character(len=GA_VALUE_LENGTH) :: InputVersion = ''  ! may drop eventually
    character(len=GA_VALUE_LENGTH) :: PGEVersion = ''
    character(len=GA_VALUE_LENGTH) :: StartUTC = ''
    character(len=GA_VALUE_LENGTH) :: EndUTC = ''
    ! character(len=GA_VALUE_LENGTH) :: GranuleMonth = ''
    ! character(len=GA_VALUE_LENGTH) :: GranuleDay   = ''
    ! character(len=GA_VALUE_LENGTH) :: GranuleYear  = ''
    integer :: GranuleMonth                  = 0
    integer :: GranuleDay                    = 0
    integer :: GranuleYear                   = 0
    real(r8) :: TAI93At0zOfGranule           = 0.d0
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

!----------------------------------------
   SUBROUTINE FillTAI93Attribute (LeapSecFileName)
!----------------------------------------

    use SDPTOOLKIT, only: mls_utctotai, pgs_td_utctotai
    use MLSSTRINGS, only: utc_to_yyyymmdd
!  Fill the TAI93 component of the global attribute based on the
!  StartUTC component

!  Arguments
    character(LEN=*), optional :: LeapSecFileName
    character(len=NameLen)     :: start_time_string
!  Local variables
    integer             :: returnStatus
    character(len=16)   :: year, month, day
!   Executable statements
    if ( GlobalAttributes%StartUTC /= ' ' ) then
      call utc_to_yyyymmdd(GlobalAttributes%StartUTC, returnStatus, &
        & year, month, day, utcAt0z=start_time_string)
      if ( present(LeapSecFileName) ) then
        returnStatus = mls_utctotai(trim(LeapSecFileName), start_time_string, &
        & GlobalAttributes%TAI93At0zOfGranule )
        if ( returnStatus /= 0 ) then
          CALL MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'Unable to get convert utc to tai using leapsecfile ' &
            & // trim(LeapSecFileName) // ' Error number ' &
            & // trim(int_to_char(returnStatus)) &
            & )
        end if
      else
        returnStatus = pgs_td_utctotai (start_time_string, &
        & GlobalAttributes%TAI93At0zOfGranule )
        if ( returnStatus /= 0 ) then
          CALL MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'Unable to get convert utc to tai using toolbox; Error number ' &
            & // trim(int_to_char(returnStatus)) &
            & )
        end if
      endif
    endif

!------------------------------------
   END SUBROUTINE FillTAI93Attribute
!------------------------------------

!------------------------------------------------------------
   SUBROUTINE gd_writeglobalattr (gridID)
!------------------------------------------------------------

      use HDFEOS5, only: HE5T_NATIVE_SCHAR, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE
      ! use he5_gdapi, only: he5_GDwrattr    ! Not coded yet
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
!     status = he5_GDwrattr(gridID, &
!      & 'InputVersion', HE5T_NATIVE_SCHAR, 1, &
!      &  GlobalAttributes%InputVersion)
      status = he5_GDwrattr(gridID, &
       & 'PGEVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%PGEVersion)
      status = he5_GDwrattr(gridID, &
       & 'StartUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%StartUTC)
      status = he5_GDwrattr(gridID, &
       & 'EndUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%EndUTC)
! >       status = he5_GDwrattr(gridID, &
! >        & 'GranuleMonth', HE5T_NATIVE_INT, 1, &
! >        &  GlobalAttributes%GranuleMonth)
! >       status = he5_GDwrattr(gridID, &
! >        & 'GranuleDay', HE5T_NATIVE_INT, 1, &
! >        &  GlobalAttributes%GranuleDay)
! >       status = he5_GDwrattr(gridID, &
! >        & 'GranuleYear', HE5T_NATIVE_INT, 1, &
! >        &  GlobalAttributes%GranuleYear)
! >       status = he5_GDwrattr(gridID, &
! >        & 'TAI93At0zOfGranule', HE5T_NATIVE_DOUBLE, 1, &
! >        &  GlobalAttributes%TAI93At0zOfGranule )
!------------------------------------------------------------
   END SUBROUTINE gd_writeglobalattr
!------------------------------------------------------------

!------------------------------------------------------------
   SUBROUTINE h5_writeglobalattr (fileID, skip_if_already_there)
!------------------------------------------------------------

      use HDF5, only: h5gclose_f, h5gopen_f
      USE MLSHDF5, only: IsHDF5AttributePresent, MakeHDF5Attribute
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf5-formatted file
! It does so at the root '/' group level

! Arguments

      INTEGER, INTENT(IN) :: fileID
      logical, intent(in), optional :: skip_if_already_there
! Local variables
      integer :: grp_id
      integer :: status
      logical :: my_skip

      ! Executable code
      my_skip = .false.
      if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
      if ( my_skip ) then
        if ( IsHDF5AttributePresent('/', fileID, 'InstrumentName') ) &
          & return
      endif
      call h5gopen_f(fileID, '/', grp_id, status)
      call MakeHDF5Attribute(grp_id, &
       & 'InstrumentName', GlobalAttributes%InstrumentName, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'ProcessLevel', GlobalAttributes%ProcessLevel, .true.)
!     call MakeHDF5Attribute(grp_id, &
!      & 'InputVersion', GlobalAttributes%InputVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'PGEVersion', GlobalAttributes%PGEVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'StartUTC', GlobalAttributes%StartUTC, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'EndUTC', GlobalAttributes%EndUTC, .true.)
      if ( GlobalAttributes%GranuleDay < 1 ) return
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleMonth', GranuleMonth_fun() , .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleDay', GranuleDay_fun(), .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleDayOfYear', GranuleDayOfYear_fun(), .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'GranuleYear', GlobalAttributes%GranuleYear, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'TAI93At0zOfGranule', GlobalAttributes%TAI93At0zOfGranule, .true.)
      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   END SUBROUTINE h5_writeglobalattr
!------------------------------------------------------------

!------------------------------------------------------------
   SUBROUTINE he5_writeglobalattr (fileID)
!------------------------------------------------------------

!     use HDF5, only: H5T_NATIVE_CHARACTER
      use HDFEOS5, only: HE5T_NATIVE_SCHAR, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE
      use he5_swapi, only: he5_EHwrglatt
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf-eos5 file

! Arguments

      INTEGER, INTENT(IN) :: fileID
      ! integer, external ::   he5_EHwrglatt
! Internal variables
      integer :: status
! Executable
      status = he5_EHwrglatt(fileID, &
       & 'InstrumentName', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InstrumentName)
      status = he5_EHwrglatt(fileID, &
       & 'ProcessLevel', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%ProcessLevel)
!     status = he5_EHwrglatt(fileID, &
!      & 'InputVersion', HE5T_NATIVE_SCHAR, 1, &
!      &  GlobalAttributes%InputVersion)
      status = he5_EHwrglatt(fileID, &
       & 'PGEVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%PGEVersion)
      status = he5_EHwrglatt(fileID, &
       & 'StartUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%StartUTC)
      status = he5_EHwrglatt(fileID, &
       & 'EndUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%EndUTC)
      ! if ( GlobalAttributes%GranuleDay == ' ') return
      ! if ( GlobalAttributes%GranuleMonth == ' ') &
      if ( GlobalAttributes%GranuleDay < 1 ) return
      status = he5_EHwrglatt(fileID, &
       & 'GranuleMonth', HE5T_NATIVE_INT, 1, &
       &  (/ GranuleMonth_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleDay', HE5T_NATIVE_INT, 1, &
       &  (/ GranuleDay_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleDayOfYear', HE5T_NATIVE_INT, 1, &
       &  (/ GranuleDayOfYear_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleYear', HE5T_NATIVE_INT, 1, &
       &  (/ GlobalAttributes%GranuleYear/) )
      status = he5_EHwrglatt(fileID, &
       & 'TAI93At0zOfGranule', HE5T_NATIVE_DOUBLE, 1, &
       &  (/ GlobalAttributes%TAI93At0zOfGranule/) )
!------------------------------------------------------------
   END SUBROUTINE he5_writeglobalattr
!------------------------------------------------------------

!----------------------------------------
   SUBROUTINE InputInputPointer (urefs, fileIDArray, fileNameArray, &
     & PCBottom, PCTop )
!----------------------------------------
!  Prepare Input for WriteInputPointer consisting of universal refs
!  This can be done for an array of fileids or of file names or both
!  that were input during the run; e.g., l1brads for level 2
!  If no universal refs are found, as always happens with my test cases,
!  then put the file names as returned by pgs_pc_getReference in their place
!  for the fileIDArray, or else the fileNameArray itself if appropriate
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
   INTEGER, EXTERNAL :: pgs_pc_getUniversalRef, pgs_pc_getReference

 ! Executable
   if (size(urefs) < 1) return
   urefs = ' '
   ref = 0
   if ( size(fileIDArray) > 0 ) then
      DO i = 1, size(fileIDArray)
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
        elseif(fileIDArray(i) > 0) then
          CALL MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'Unable to get Universal Ref for file ID ' &
            & // trim(int_to_char(fileIDArray(i))) // ' Error number ' &
            & // trim(int_to_char(returnStatus)) &
            & )
          returnStatus = pgs_pc_getReference(fileIDArray(i), &
           & version, sval)
          IF (returnStatus == PGS_S_SUCCESS) THEN 
            ref = min(ref+1, size(urefs))
            urefs(ref) = sval
          else
            CALL MLSMessage(MLSMSG_Error, ModuleName, &
              & 'Unable to get even a file Ref for file ID ' &
              & // trim(int_to_char(fileIDArray(i))) // ' Error number ' &
              & // trim(int_to_char(returnStatus)) &
              & )
          endif
        ENDIF                                   
      ENDDO
   endif

   if ( size(fileNameArray) > 0 ) then
      DO i = 1, size(fileNameArray)
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
        else                                   
           ref = min(ref+1, size(urefs))
           urefs(ref) = fileNameArray(i)                     
        ENDIF                                   
      ENDDO
   endif

!------------------------------------
   END SUBROUTINE InputInputPointer
!------------------------------------

!------------------------------------------------------------
   SUBROUTINE sw_writeglobalattr (swathID)
!------------------------------------------------------------

      use HDFEOS5, only: HE5T_NATIVE_SCHAR, HE5T_NATIVE_INT, HE5T_NATIVE_DOUBLE
      use HE5_SWAPI, only: he5_SWwrattr
! Brief description of subroutine
! This subroutine writes the global attributes for an hdf-eos5 swath

! Arguments

      INTEGER, INTENT(IN) :: swathID
!     integer, external ::   he5_SWwrattr
! Internal variables
      integer :: status
! Executable
      status = he5_SWwrattr(swathID, &
       & 'InstrumentName', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%InstrumentName)
      status = he5_SWwrattr(swathID, &
       & 'ProcessLevel', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%ProcessLevel)
!     status = he5_SWwrattr(swathID, &
!      & 'InputVersion', HE5T_NATIVE_SCHAR, 1, &
!      &  GlobalAttributes%InputVersion)
      status = he5_SWwrattr(swathID, &
       & 'PGEVersion', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%PGEVersion)
      status = he5_SWwrattr(swathID, &
       & 'StartUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%StartUTC)
      status = he5_SWwrattr(swathID, &
       & 'EndUTC', HE5T_NATIVE_SCHAR, 1, &
       &  GlobalAttributes%EndUTC)
      if ( GlobalAttributes%GranuleDay < 1 ) return
      if ( GlobalAttributes%GranuleMonth > 0 ) then
        status = he5_SWwrattr(swathID, &
         & 'GranuleMonth', HE5T_NATIVE_INT, 1, &
         &  (/ GlobalAttributes%GranuleMonth/) )
        status = he5_SWwrattr(swathID, &
         & 'GranuleDay', HE5T_NATIVE_INT, 1, &
         &  (/ GlobalAttributes%GranuleDay/) )
      else
        status = he5_SWwrattr(swathID, &
         & 'GranuleDayOfYear', HE5T_NATIVE_INT, 1, &
         &  (/ GlobalAttributes%GranuleDay/) )
      endif
      status = he5_SWwrattr(swathID, &
       & 'GranuleYear', HE5T_NATIVE_INT, 1, &
       &  (/ GlobalAttributes%GranuleYear/) )
      status = he5_SWwrattr(swathID, &
       & 'TAI93At0zOfGranule', HE5T_NATIVE_DOUBLE, 1, &
       &  (/ GlobalAttributes%TAI93At0zOfGranule/) )
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
   SUBROUTINE WritePCF2Hdr (file, anText, hdfVersion, fileType)
!----------------------------------------
      use HDF5, only: H5F_ACC_RDWR_F, &
        & h5fopen_f, h5fclose_f

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS file as an annotation.

! Arguments

      CHARACTER (LEN=FileNameLen), INTENT(IN) :: file 

      CHARACTER (LEN=1), POINTER              :: anText(:)
      integer, intent(in), optional           :: hdfVersion
      character(len=2), intent(in), optional  :: fileType ! 'sw', 'gd', 'hdf'
      ! logical, intent(in), optional         :: isHDFEOS

! Parameters
      integer :: fileID
      integer :: my_hdfVersion
      ! logical :: myisHDFEOS
      integer :: record_length
      integer :: status
      character (len=2) :: the_type
! Executable
      my_hdfVersion = PCFHDR_DEFAULT_HDFVERSION
      if ( present(hdfVersion) ) my_hdfVersion = hdfVersion
      ! myisHDFEOS = .false.
      ! if ( present(isHDFEOS) ) myisHDFEOS = isHDFEOS
      the_type = 'hd'
      if ( present(fileType) ) the_type = lowercase(fileType)
      select case(my_hdfVersion)
      case (HDFVERSION_4)
        call WritePCF2Hdr_hdf4 (file, anText)
      case (HDFVERSION_5)
        if ( the_type == 'sw' ) then
          fileID = mls_io_gen_openF('swopen', .TRUE., status, &
           & record_length, DFACC_RDWR, FileName=trim(file), &
           & hdfVersion=hdfVersion, debugOption=.false. )
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error opening hdfeos5 swath file for annotating with PCF' )
          call WritePCF2Hdr_hdfeos5 (fileID, anText)
          status = mls_io_gen_closeF('swclose', fileID, &
            & hdfVersion=hdfVersion)
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error closing hdfeos5 swath file for annotating with PCF' )
        elseif ( the_type == 'gd' ) then
          fileID = mls_io_gen_openF('gdopen', .TRUE., status, &
           & record_length, DFACC_RDWR, FileName=trim(file), &
           & hdfVersion=hdfVersion, debugOption=.false. )
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error opening hdfeos5 grid file for annotating with PCF' )
          call WritePCF2Hdr_hdfeos5 (fileID, anText)
          status = mls_io_gen_closeF('gdclose', fileID, &
            & hdfVersion=hdfVersion)
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error closing hdfeos5 grid file for annotating with PCF' )
        else
          call h5fopen_f(trim(file), H5F_ACC_RDWR_F, fileID, status)
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error opening hdf5 file for annotating with PCF' )
          ! call WritePCF2Hdr_hdf5 (file, anText)
          call WritePCF2Hdr_hdf5 (fileID, anText)
          call h5fclose_f(fileID, status)
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error closing hdf5 file for annotating with PCF' )
        endif
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
   SUBROUTINE WritePCF2Hdr_hdf5 (fileID, anText)
!----------------------------------------

      use HDF5, only: h5gclose_f, h5gopen_f
      USE MLSHDF5, only: MakeHDF5Attribute, SaveAsHDF5DS
! Brief description of subroutine
! This subroutine writes the PCF into an HDF5 file as 
! (1) a datset if MAKEDATASET is TRUE
! (2) an attribute if MAKEATTRIBUTE is TRUE
! It does (2) at the root '/' group level, treating the file as if
! it were a plain hdf5 file
! If the file were in fact an hdfeos5 file
! this may result in the dread hybrid file syndrome which may
! confuse some hdf-eos5 readers (not tested)

! Arguments

      CHARACTER (LEN=1), POINTER              :: anText(:)

! Local variables
      integer :: fileID
      integer :: grp_id
      integer :: status
      CHARACTER (LEN=1), dimension(:), POINTER              :: an40
      integer :: how_big
      logical, parameter :: USELENGTHONECHARS = .true.
      logical, parameter :: MAKEDATASET = .true.
      logical, parameter :: MAKEATTRIBUTE = .true.
      character(len=size(anText)) :: anScalar
      ! Executable
      if ( MAKEDATASET ) then
        anScalar = transfer(anText, anScalar)
        call SaveAsHDF5DS ( fileID, '/PCF', anScalar )
      endif
      if ( .not. MAKEATTRIBUTE ) return
      call h5gopen_f(fileID, '/', grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error opening hdf5 file root group for annotating with PCF' )
      if ( USELENGTHONECHARS ) then
        call MakeHDF5Attribute(grp_id, &
         & PCFATTRIBUTENAME, anText, .true.)
         ! & 'PCF file text', anText, .true.)
      else
        ! Find how big an40 must be to hold anText
        how_big = 1 + (size(anText)-1)/40
        allocate(an40(how_big), stat=status)
        if ( status /= 0 ) &
          & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & MLSMSG_Allocate // 'an40 for annotating hdfeos5 PCF' )
        an40 = ' '
        ! Do some nonsense here
        call MakeHDF5Attribute(grp_id, &
         & 'PCF file text', an40, .true.)
        deallocate(an40, stat=status)
        if ( status /= PGS_S_SUCCESS) &
          & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & MLSMSG_DeAllocate // 'an40 annotating hdf5 with PCF' )
      endif
      call h5gclose_f(grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error closing hdf5 file root group for annotating with PCF' )
!-----------------------------
   END SUBROUTINE WritePCF2Hdr_hdf5
!-----------------------------

!----------------------------------------
   SUBROUTINE WritePCF2Hdr_hdfeos5 (fileID, anText)
!----------------------------------------

      use HDFEOS5, only: HE5T_NATIVE_SCHAR
      use he5_swapi, only: he5_EHwrglatt
! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS5 file as an attribute.
! It does so as file level attributes
! It does not do so as a dataset because that would create a hybrid file
! For unclear reasons the attributes are stored as an array of 40-length chars
! Unless the parameter USELENGTHONECHARS is true
! Arguments

      CHARACTER (LEN=1), POINTER              :: anText(:)
      integer, intent(in)                     :: fileID

! Local variables
      integer :: divisor
      integer :: grp_id
      integer :: i
      integer :: remainder
      integer :: status
      CHARACTER (LEN=40), POINTER             :: an40(:)
      integer :: how_big
      logical, parameter :: USELENGTHONECHARS = .true.
! Executable
      if ( USELENGTHONECHARS ) then
        status = he5_EHwrglatt(fileID, &
         & PCFATTRIBUTENAME, HE5T_NATIVE_SCHAR, size(anText), anText)
         ! & 'PCF', HE5T_NATIVE_SCHAR, size(anText), anText)
        if ( status /= PGS_S_SUCCESS) &
          & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Error annotating with PCF using length one chars' )
        return
      endif
      ! Find how big an40 must be to hold anText
      how_big = 1 + (size(anText)-1)/40
      allocate(an40(how_big), stat=status)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & MLSMSG_Allocate // 'an40 for annotating hdfeos5 PCF' )
      an40 = ' '
      ! Do some nonsense here
      ! (This is clever, but is it necessary?)
      do i=1, size(anText)
        divisor = (i-1) / 40
        remainder = mod(i-1, 40)
        an40(divisor+1)(remainder+1:remainder+1) = anText(i)
      enddo
      ! print *, 'size(anText) ', size(anText)
      ! print *, 'how_big ', how_big
      ! print *, 'an40 '
      ! do i=1, how_big
      !  print *, trim(an40(i))
      ! enddo
      status = he5_EHwrglatt(fileID, &
       & PCFATTRIBUTENAME, HE5T_NATIVE_SCHAR, how_big, an40)
       ! & 'PCF', HE5T_NATIVE_SCHAR, how_big, an40)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error annotating with PCF using length 40 chars' )
      deallocate(an40, stat=status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & MLSMSG_DeAllocate // 'an40 annotating hdfeos5 with PCF' )
!-----------------------------
   END SUBROUTINE WritePCF2Hdr_hdfeos5
!-----------------------------

  function GranuleDayOfYear_fun () result (dayOfYear)
    ! Arguments
    integer :: month
    ! Local variables
    integer :: year
    integer :: dayOfYear
    integer :: status
    character (len=UTC_A_VALUE_LENGTH) :: asciiutc_a
    character (len=UTC_B_VALUE_LENGTH) :: asciiutc_b
    ! Executable
    if ( GlobalAttributes%GranuleMonth <= 0 .or. .not. UseSDPToolkit ) then
      dayOfYear = GlobalAttributes%GranuleDay
    else
      asciiutc_a = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_atob(asciiutc_a, asciiutc_b)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc A to B formats')
      call utc_to_yyyymmdd(asciiutc_b, status, &
        & year, month, dayOfYear) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to extract year, day of year from utc B format')
    endif
  end function GranuleDayOfYear_fun

  function GranuleDay_fun () result (day)
    ! Arguments
    integer :: month
    ! Local variables
    integer :: year
    integer :: day
    integer :: status
    character (len=UTC_A_VALUE_LENGTH) :: asciiutc_a
    character (len=UTC_B_VALUE_LENGTH) :: asciiutc_b
    ! Executable
    if ( GlobalAttributes%GranuleMonth > 0 .or. .not. UseSDPToolkit ) then
      day = GlobalAttributes%GranuleDay
    else
      asciiutc_b = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_btoa(asciiutc_b, asciiutc_a)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc B to A formats')
      call utc_to_yyyymmdd(asciiutc_a, status, &
        & year, month, day) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to extract year, month, day from utc A format')
    endif
  end function GranuleDay_fun

  function GranuleMonth_fun () result (month)
    ! Arguments
    integer :: month
    ! Local variables
    integer :: year
    integer :: day
    integer :: status
    character (len=UTC_A_VALUE_LENGTH) :: asciiutc_a
    character (len=UTC_B_VALUE_LENGTH) :: asciiutc_b
    ! Executable
    if ( GlobalAttributes%GranuleMonth > 0 .or. .not. UseSDPToolkit ) then
      month = GlobalAttributes%GranuleMonth
    else
      asciiutc_b = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_btoa(asciiutc_b, asciiutc_a)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc B to A formats')
      call utc_to_yyyymmdd(asciiutc_a, status, &
        & year, month, day) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to extract year, month, day from utc A format')
    endif
  end function GranuleMonth_fun

  function int_to_char (int) result (chars)
    ! Arguments
    integer , intent(in) :: int
    character(len=16) :: chars
    ! Executable
    write(chars, * ) int
    chars = adjustl(chars)
  end function int_to_char

!================
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module PCFHdr
!================

!# $Log$
!# Revision 2.17  2003/03/19 23:57:56  jdone
!# changed length of string fileType from wildcard to 2
!#
!# Revision 2.16  2003/03/11 00:20:14  pwagner
!# can WritePCF2Hdr for swath, grid, or hdf files
!#
!# Revision 2.15  2003/03/07 00:37:24  pwagner
!# Write GranuleDay -Month and -Year even if StartUTC in format yyy-ddd
!#
!# Revision 2.14  2003/02/27 21:52:48  pwagner
!# Added FillTAI93Attribute; tweaks to PCF as attribute; unsatisfactory for hdfeos5
!#
!# Revision 2.13  2003/02/12 21:47:05  pwagner
!# Optional skip_if_already_there arg to h5_writeglobalattr
!#
!# Revision 2.12  2003/02/10 22:07:23  pwagner
!# Granule Month, Day, Year now ints; he5 global attributes leave file pure hdfeos5
!#
!# Revision 2.11  2003/02/08 00:32:11  pwagner
!# Gets he5_EHwrglatt from he5_swapi
!#
!# Revision 2.10  2003/02/06 00:30:19  pwagner
!# Added h5_ and he5_writeglobalattr
!#
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
