! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. user has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!===============================================================================
module PCFHdr
!===============================================================================
! Dsepite an unfortunate choice of name, which fails to clearly 
! indicate it, this module contains ways to annotate, write
! attributes, and otherwise do miscellaneous things to mls product files.
! It might have been better named PCFHdrAndGlobalAttributes
! or split off global attribute stuff into a separate module

! See Background below for background and suggestion

  use Dates_Module, only: Utc_To_Date, Utc_To_Yyyymmdd, Utcform, &
    & Yyyymmdd_To_Doy, YyyyDoy_To_Mmdd
  use HDF, only: Dfacc_Rdwr, Dfacc_Write, An_File_Desc
  use HighOutput, only: OutputNamedValue
  use Intrinsic, only: L_HDFeos, L_HDF, L_Swath
  use MLSCommon, only: Filenamelen, MLSFile_T, Filenamelen, Namelen
  use MLSFiles, only: GetPCFromRef, HDFversion_4, HDFversion_5, &
    & IniTializeMLSFile, MLS_CloseFile, MLS_OpenFile, MLS_OpenFile, MLS_CloseFile
  use MLSKinds, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Error, MLSMSG_Fileopen, &
    & MLSMSG_Warning
  use MLSStrings, only: Lowercase
  use Output_M, only: Output
  use Sdptoolkit, only: Pgsd_Pc_Uref_Length_Max, Pgs_S_Success, &
    & Pgsd_Met_Group_Name_L, Pgs_Io_Gen_Closef, Pgs_Io_Gen_Openf, &
    & Pgsd_Io_Gen_Rdirunf, Pgs_Pc_Getreference, &
    & Pgs_Td_Asciitime_Atob, Pgs_Td_Asciitime_Btoa, &
    & Max_Orbits
  implicit none
  public :: GlobalAttributes_T, SomeGlobalAttributes_T, &
    & CreatePCFAnnotation, DumpGlobalAttributes, &
    & FillTAI93Attribute, &
    & GranuleDay, GranuleDayOfYear, GranuleMonth, GranuleYear, &
    & H5_ReadGlobalAttr, H5_ReadMLSFileAttr, He5_ReadMLSFileAttr, &
    & H5_WriteGlobalAttr, He5_WriteglobalAttr, He5_ReadGlobalAttr, &
    & H5_WriteMLSFileAttr, He5_WriteMLSFileAttr, &
    & InputInputPointer, ProcessLevelFun, SomeToGlobalAttributes, &
    & WriteInputPointer, &
    & WriteLeapSecHDFEOSAttr, WriteLeapSecHDF5DS, WritePCF2Hdr, &
    & WriteUTCPoleHDFEOSAttr, WriteUTCPoleHDF5DS
   private

! === (start of toc) ===                                                 
!     c o n t e n t s                                                    
!     - - - - - - - -                                                    

!     (data types and parameters)
! Inputptr_string_length   String length used by Inputpointer procedures
! Ga_value_length          String length used by GlobalAttributes_T
! MiscNotesLength          String length used by MiscNotes field
! Utc_a_value_length       String length used to encode utc version 'a'
! Utc_b_value_length       String length used to encode utc version 'b'
! GlobalAttributes         Which attributes to write to product files
! SomeGlobalAttributes     Which of them can be read from an opts file

!     (subroutines and functions)
! CreatePCFAnnotation      Read the PCF file into a character array
! DumpGlobalAttributes     Dumps the global attributes
! FillTAI93Attribute       Fill the TAI93 component of the global attribute 
!                           based on theStartUTC component
! GranuleDay               Return the granule's day in its month
! GranuleDayOfYear         Return the granule's day in its year
! GranuleMonth             Return the granule's month
! GranuleYear              Return the granule's year
! H5_readglobalattr        Reads the global attributes from an hdf5-formatted file
! H5_writeglobalattr       Writes the global attributes to an hdf5-formatted file
! He5_writeglobalattr      Writes the global attributes to an hdfeos5-formatted file
! He5_readglobalattr       Reads the global attributes from an hdfeos5-formatted file
! InputInputpointer        Prepare Input for WriteInputpointer
! ProcessLevelFun          Return the Current ProcessLevel (Default: 'L2')
! SomeToGlobalAttributes   Copy some of the attributes read into the global ones
! Sw_writeglobalattr       Writes the global attributes for an hdfeos5 swath
! WriteInputpointer        Write Inputpointer metadata
! WriteLeapSecHDFEOSAttr   Write contents of leapsec file as hdfeos5 attribute
! WriteLeapSecHDF5DS       Write contents of leapsec file as hdf5 dataset
! WritePCF2Hdr             Write the PCF into an HDF or HDF-EOS file
! WriteutcPoleHDFEOSAttr   Write contents of utcPole file as hdfeos5 attribute
! WriteutcPoleHDF5DS       Write contents of utcPole file as hdf5 dataset
! === (end of toc) ===      
                                             
! === (start of api) ===
! CreatePCFAnnotation ( int mlspcfN_pcf_start, char* anText(:) )
! DumpGlobalAttributes     Dumps the global attributes
! FillTAI93Attribute       Fill the TAI93 component of the global attribute 
!                           based on theStartUTC component
! GranuleDay               Return the granule's day in its month
! GranuleDayOfYear         Return the granule's day in its year
! GranuleMonth             Return the granule's month
! GranuleYear              Return the granule's year
! H5_readglobalattr        Reads the global attributes from an hdf5-formatted file
! H5_writeglobalattr       Writes the global attributes to an hdf5-formatted file
! He5_writeglobalattr      Writes the global attributes to an hdfeos5-formatted file
! He5_readglobalattr       Reads the global attributes from an hdfeos5-formatted file
! InputInputpointer        Prepare Input for WriteInputpointer
! Sw_writeglobalattr       Writes the global attributes for an hdfeos5 swath
! WriteInputpointer        Write Inputpointer metadata
! WriteLeapSecHDFEOSAttr   Write contents of leapsec file as hdfeos5 attribute
! WriteLeapSecHDF5DS       Write contents of leapsec file as hdf5 dataset
! WritePCF2Hdr             Write the PCF into an HDF or HDF-EOS file
! WriteutcPoleHDFEOSAttr   Write contents of utcPole file as hdfeos5 attribute
! WriteutcPoleHDF5DS       Write contents of utcPole file as hdf5 dataset
! === (end of api) ===

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Contents:

! Background:  This module still contains subroutines for writing the PCF as an 
! annotation to HDF files, a use now mostly obsolete.
! Thus its inappropriate name.

! It also contains 
! (a) two routines that prepare and write the input pointer
! to metadata, still used by Level 1
! (b) routines for dumping, reading and writing file-level attributes
! in particular hdfeos5 and plain hdf5
! (c) functions returning the day, month, etc. of the data date, aka granule

! Suggestion:
! Create a brand new module, named something like GlobalAttributes.f90
! and move the non-obsolete procedures and datatypes there.
  integer, parameter, public :: INPUTPTR_STRING_LENGTH = PGSd_PC_UREF_LENGTH_MAX
  integer, parameter, public :: GA_VALUE_LENGTH = 40
  integer, parameter, public :: MiscNotesLENGTH = 4096
  integer, parameter, public :: UTC_A_VALUE_LENGTH = 27
  integer, parameter, public :: UTC_B_VALUE_LENGTH = 25
  character(len=*), parameter, private :: PCFATTRIBUTENAME = 'PCF'
  character(len=*), parameter, private :: PCFPATHNAME = '/PCF'
  character(len=*), parameter, private :: HDFEOSINPTPTRVALUE = 'Found at ' // &
    & '/HDFEOS/ADDITIONAL/FILE_ATTRIBUTES/PCF'
  character(len=*), parameter, private :: HDFINPTPTRVALUE = 'Found at ' // &
    & '/PCF'
  character(len=*), parameter, private :: DEFAULTPROCESSLEVEL = 'L2'

  ! No MAF number can ever be this big (as in L1BData module)
  integer, parameter :: BIGGESTMAFCTR = huge(0)/2

   ! May we get some of these from a new module, e.g. MLSLibOptions? 
  type SomeGlobalAttributes_T
    character(len=GA_VALUE_LENGTH) :: InstrumentName = 'MLS Aura'
    character(len=GA_VALUE_LENGTH) :: ProcessLevel   = ''
    character(len=GA_VALUE_LENGTH) :: HostName       = ''  ! E.g. 'lightspeed'
    character(len=GA_VALUE_LENGTH) :: PGEVersion     = ''
    character(len=MiscNotesLENGTH) :: MiscNotes      = ''
    character(len=GA_VALUE_LENGTH) :: StartUTC       = ''
    character(len=GA_VALUE_LENGTH) :: EndUTC         = ''
    character(len=GA_VALUE_LENGTH) :: DOI            = '' ! E.g., '10.5083/AURA/MLS/DATA201'
    character(len=GA_VALUE_LENGTH) :: productionLoc  = ' '
  end type SomeGlobalAttributes_T

  ! Note that we repeat some of the above
  ! Sad but necessary because we want to be able to read some
  ! global attributes from a text file
  ! yet Fortran in its mighty wisdom forbids io if  any list item is a pointer
  type GlobalAttributes_T
    character(len=FileNameLen)     :: FileAttributesCopiedFrom = ' '  
    integer                        :: OrbNum(max_orbits)              
    real(r8)                       :: OrbPeriod(max_orbits)           
    integer, pointer, dimension(:,:) :: OrbNumDays => Null()
    real(r8), pointer, dimension(:,:) :: OrbPeriodDays => Null()
    character(len=GA_VALUE_LENGTH) :: InstrumentName = 'MLS Aura'
    character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''
    character(len=GA_VALUE_LENGTH) :: HostName = ''  ! E.g. 'lightspeed'
    character(len=GA_VALUE_LENGTH) :: PGEVersion = ''
    character(len=MiscNotesLENGTH) :: MiscNotes = ''
    character(len=GA_VALUE_LENGTH) :: StartUTC = ''
    character(len=GA_VALUE_LENGTH) :: EndUTC = ''
    character(len=GA_VALUE_LENGTH) :: DOI = '' ! E.g., '10.5083/AURA/MLS/DATA201'
    character(len=GA_VALUE_LENGTH) :: productionLoc = ' '
    character(len=GA_VALUE_LENGTH) :: geos5_type = ' '
    ! Note that the next 3 may refer to either of 2 formats
    ! (a) yyyy-mm-dd if GranuleMonth > 0
    ! (b) yyyy-Doy if GranuleMonth < 0
    integer                        :: GranuleMonth                  = 0 ! < 0 if day is Doy
    integer                        :: GranuleDay                    = 0 ! 0 < Day < 32 of month > 0
    integer                        :: GranuleYear                   = 0
    integer                        :: DayOfYear                     = 0

    real(r8)                       :: TAI93At0zOfGranule            = 0.d0
    integer                        :: FirstMAFCtr                   = BIGGESTMAFCTR
    integer                        :: LastMAFCtr                    = 0
    character(len=GA_VALUE_LENGTH) :: PhaseNames                    = ' '
    character(len=GA_VALUE_LENGTH) :: ForwardModelNames             = ''
    ! character(len=GA_VALUE_LENGTH) :: ProductionLocation            = ''
    integer                        :: NumCompletedChunks            = 0
    integer                        :: NumFailedChunks               = 0
    character(len=GA_VALUE_LENGTH) :: FailedChunks                  = ''
    character(len=GA_VALUE_LENGTH) :: FailedMachines                = ''
    character(len=GA_VALUE_LENGTH) :: FailedMsgs                    = ''
  end type GlobalAttributes_T

  ! This variable describes the global attributes
  type (SomeGlobalAttributes_T), public, save :: SomeGlobalAttributes
  type (GlobalAttributes_T), public, save :: GlobalAttributes
  ! use this in case hdfVersion omitted from call to WritePCF2Hdr
  ! E.g., in level 3 prior to conversion
  integer, public, save            :: PCFHDR_DEFAULT_HDFVERSION = HDFVERSION_5
  logical, parameter               :: DEBUG = .false.

  interface DumpGlobalAttributes
    module procedure DumpGlobalAttributes_global, dumpGlobalAttributes_local
  end interface
  
  interface GranuleDay
    module procedure GranuleDay_Fun
  end interface
  
  interface GranuleDayOfYear
    module procedure GranuleDayOfYear_fun
  end interface
  
  interface GranuleMonth
    module procedure GranuleMonth_Fun
  end interface
  
  interface GranuleYear
    module procedure GranuleYear_Fun
  end interface
  
  interface h5_writeglobalattr
    module procedure h5_writeglobalattr_fileID, h5_writeglobalattr_MLSFile
  end interface
  
  interface he5_writeglobalattr
    module procedure he5_writeglobalattr_fileID, he5_writeglobalattr_MLSFile
  end interface
  
contains

!--------------------------------------------------------------
   subroutine CreatePCFAnnotation ( mlspcfN_pcf_start, anText )
!--------------------------------------------------------------

      use Allocate_Deallocate, only: Allocate_Test

! Brief description of subroutine
! This subroutine stores the PCF as an annotation for writing to file headers.

! Arguments

      integer, intent(in) :: mlspcfN_pcf_start

      character (len=1), pointer :: anText(:)

! Parameters

! Functions

      integer, external :: Pgs_pc_getFileSize

! Variables

      character (len=10) :: mnemonic
      character (len=480) :: msg, msr

      integer :: ios, pcfHandle, returnStatus, size, version

! Get the size of the PCF

      version = 1
      returnStatus = Pgs_pc_getFileSize(mlspcfN_pcf_start, version, size)
      nullify( anText )
      call allocate_test ( anText, size, "anText", moduleName )

! Open the PCF for reading

      version = 1
      returnStatus = Pgs_io_gen_openF (mlspcfN_pcf_start, PGSd_IO_Gen_RDirUnf, &
                                       size, pcfHandle, version)
      IF (returnStatus /= PGS_S_SUCCESS) THEN
         call Pgs_smf_getMsg(returnStatus, mnemonic, msg)
         msr = mnemonic // ':  ' // msg
         CALL MLSMessage(MLSMSG_Error, ModuleName, msr)
      END IF

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
   end subroutine CreatePCFAnnotation
!------------------------------------

!----------------------------------------
   subroutine DumpGlobalAttributes_global
     call DumpGlobalAttributes_local ( GlobalAttributes )
   end subroutine DumpGlobalAttributes_global
!----------------------------------------
   
   subroutine DumpGlobalAttributes_local ( GlobalAttributes )
     use Dates_Module, only: MonthName
     use HighOutput, only: AddRow, AddRow_Divider, AddRow_Header, &
       & OutputTable, StartTable
     type(GlobalAttributes_T), intent(in)       :: GlobalAttributes
     ! Local variables
     integer                                    :: DayofYear, month
     character(len=16)                          :: month_name
     ! Executable
     DayOfYear = GranuleDayOfYear_fun()
     month_name = 'unknown'
     month = abs(GlobalAttributes%GranuleMonth )
     if ( month > 0 .and. month < 13) month_name = MonthName(month)
     call startTable
     call addRow_header ( 'Global attributes', 'c' )
     call addRow_divider ( '-' )
     call addRow ('Start UTC                    ', trim( GlobalAttributes%StartUTC        ) )
     call addRow ('End UTC                      ', trim( GlobalAttributes%EndUTC          )) 
     call addRow ('Granule day of year          ', DayOfYear )
     call addRow ('Granule month                ', trim(Month_Name) )
     call addRow ('Granule day                  ', GlobalAttributes%GranuleDay  )
     call addRow ('Granule year                 ', GlobalAttributes%GranuleYear )
     call addRow ('FirstMAFCtr                  ', GlobalAttributes%FirstMAFCtr )
     call addRow ('LastMAFCtr                   ', GlobalAttributes%LastMAFCtr )
     call addRow ('NumCompletedChunks           ', GlobalAttributes%NumCompletedChunks )
     call addRow ('NumFailedChunks              ', GlobalAttributes%NumFailedChunks )
     call addRow ('FailedMsgs                   ', trim( GlobalAttributes%FailedMsgs          ), &
       & BlocLen=38, options='-w' )
     call addRow ('FailedChunks                 ', trim( GlobalAttributes%FailedChunks        ), &
       & BlocLen=38, options='-w' )
     call addRow ('FailedMachines               ', trim( GlobalAttributes%FailedMachines      ), &
       & BlocLen=38, options='-w' )
     call addRow ('FileAttributesCopiedFrom     ', trim( GlobalAttributes%FileAttributesCopiedFrom ) )
     call addRow ('Equator crossing time (tai93)', GlobalAttributes%TAI93At0zOfGranule )
     call addRow ('ProductionLocation           ', trim( GlobalAttributes%productionLoc ))
     call addRow ('HostName                     ', trim( GlobalAttributes%HostName        ))
     call addRow ('InstrumentName               ', trim( GlobalAttributes%InstrumentName  ) )
     call addRow ('Process level                ', trim( GlobalAttributes%ProcessLevel    ) )
     call addRow ('PGE version                  ', trim( GlobalAttributes%PGEVersion      ) )
     call addRow ('DOI                          ', trim( GlobalAttributes%DOI             ) )
     call addRow ('Misc Notes                   ', trim( GlobalAttributes%MiscNotes       ), &
       & BlocLen=38, options='-w', wrappingChar=',' )
     call addRow ('DOI                          ', trim( GlobalAttributes%DOI             ))
     call outputTable ( sep='|', border='-' )
     call outputNamedValue('Orbit numbers       ', GlobalAttributes%orbNum)
     call outputNamedValue('Orbit Periods       ', GlobalAttributes%orbPeriod)
   end subroutine DumpGlobalAttributes_local

!----------------------------------------
   subroutine FillTAI93Attribute ( LeapSecFileName )
!----------------------------------------

    use SDPToolkit, only: MLS_UTCToTAI, PGS_TD_UTCToTAI
!  Fill the TAI93 component of the global attribute based on the
!  StartUTC component

!  Arguments
    character(len=*), optional :: LeapSecFileName
    character(len=NameLen)     :: start_time_string
!  Local variables
    integer             :: returnStatus
    character(len=16)   :: date
!   Executable statements
    if ( GlobalAttributes%StartUTC /= ' ' ) then
      call utc_to_date(GlobalAttributes%StartUTC, returnStatus, &
        & date, utcAt0z=start_time_string)
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
    if ( DEBUG ) call outputnamedValue ( 'FillTAI93Attribute', &
      & GlobalAttributes%TAI93At0zOfGranule )

!------------------------------------
   end subroutine FillTAI93Attribute
!------------------------------------

  ! ------------ GranuleDayOfYear ---------------
  ! These functions return day, month, etc. of the current data date
  function GranuleDayOfYear_fun () result (dayOfYear)
    ! Arguments
    integer :: month
    ! Local variables
    integer :: year
    integer :: dayOfYear
    integer :: status
    character (len=UTC_A_VALUE_LENGTH) :: asciiutc_a
    character (len=UTC_B_VALUE_LENGTH) :: asciiutc_b
    character(len=1) :: whatUTCForm
    ! Executable
    dayofYear = -999
    if ( GlobalAttributes%GranuleYear <= 0 ) then
      whatUTCForm = utcForm(GlobalAttributes%StartUTC)
      if ( DEBUG ) then
        call outputNamedValue( 'utc form', whatUTCForm )
        call outputNamedValue( 'StartUTC', trim(GlobalAttributes%StartUTC) )
      endif
      select case (whatUTCForm)
      case ('a')
        asciiutc_a = GlobalAttributes%StartUTC
        status = pgs_td_asciitime_atob(asciiutc_a, asciiutc_b)
        if ( status /= 0 ) then
          call outputNamedValue( 'StartUTC', GlobalAttributes%StartUTC )
          CALL MLSMessage(MLSMSG_Error, ModuleName, &
          & 'Unable to convert utc A to B formats')
        endif
      case ('b')
        asciiutc_b = GlobalAttributes%StartUTC
      case default
        asciiutc_b = GlobalAttributes%StartUTC
      end select
      call utc_to_yyyymmdd(asciiutc_b, status, &
        & year, month, dayOfYear) 
      if ( DEBUG ) call outputNamedValue( 'asciiutc_b', trim(asciiutc_b) )
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to extract year, day of year from utc B format')
    elseif ( GlobalAttributes%GranuleMonth <= 0 ) then
      dayOfYear = GlobalAttributes%GranuleDay
    else
      call yyyymmdd_to_Doy ( GlobalAttributes%GranuleYear, &
        & GlobalAttributes%GranuleMonth, &
        & GlobalAttributes%GranuleDay, dayOfYear )
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
    day = 0
    if ( GlobalAttributes%GranuleYear <= 0 ) then
      asciiutc_b = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_btoa(asciiutc_b, asciiutc_a)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc B to A formats')
      call utc_to_yyyymmdd(asciiutc_a, status, &
        & year, month, day) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to extract year, month, day from utc A format')
    elseif ( GlobalAttributes%GranuleMonth > 0 ) then
      day = GlobalAttributes%GranuleDay
    else
      call yyyyDoy_to_mmdd ( GlobalAttributes%GranuleYear, &
        & Month, Day, &
        & GlobalAttributes%GranuleDay )
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
    month = 0
    if ( GlobalAttributes%GranuleYear <= 0 ) then
      asciiutc_b = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_btoa(asciiutc_b, asciiutc_a)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc B to A formats')
      call utc_to_yyyymmdd(asciiutc_a, status, &
        & year, month, day) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to extract year, month, day from utc A format')
    elseif ( GlobalAttributes%GranuleMonth > 0 ) then
      month = GlobalAttributes%GranuleMonth
    else
      call yyyyDoy_to_mmdd ( GlobalAttributes%GranuleYear, &
        & Month, Day, &
        & GlobalAttributes%GranuleDay )
    endif
  end function GranuleMonth_fun

  function GranuleYear_fun () result (Year)
    ! Arguments
    integer :: Year
    ! Local variables
    integer :: month
    integer :: day
    integer :: status
    character (len=UTC_A_VALUE_LENGTH) :: asciiutc_a
    character (len=UTC_B_VALUE_LENGTH) :: asciiutc_b
    ! Executable
    Year = 0
    if ( GlobalAttributes%GranuleYear <= 0 ) then
      asciiutc_b = GlobalAttributes%StartUTC
      status = pgs_td_asciitime_btoa(asciiutc_b, asciiutc_a)
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to convert utc B to A formats')
      call utc_to_yyyymmdd(asciiutc_a, status, &
        & year, Year, day) 
      if ( status /= 0 ) &
        & CALL MLSMessage(MLSMSG_Warning, ModuleName, &
        & 'Unable to extract year, Year, day from utc A format')
    else
      Year = GlobalAttributes%GranuleYear
    endif
  end function GranuleYear_fun

!------------------------------------------------------------
   subroutine h5_readMLSFileAttr ( MLSFile )
!------------------------------------------------------------

      use HDF5, only: H5GClose_f, H5GOpen_f
      use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributePresent
! Brief description of subroutine
! This subroutine reads the components of an MLSFile_t 
! as attributes from an hdf5-formatted file
! It does so at the root '/' group level

! Arguments

      type(MLSFile_T)       :: MLSFile
! Local variables
      integer :: fileID
      integer :: grp_id
      integer :: status
      integer, dimension(2) :: ints

      ! Executable code
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      fileID = MLSFile%FileID%f_id
      if ( IsHDF5AttributePresent('/', fileID, 'ShortName') ) &
        & return
      call h5gopen_f(fileID, '/', grp_id, status)
      MLSFile%FileID%grp_id = grp_id
      MLSFile%FileID%sd_id = 0
      call GetHDF5Attribute(MLSFile, &
       & 'content', MLSFile%content)
      call GetHDF5Attribute(MLSFile, &
       & 'lastOperation', MLSFile%lastOperation)
      call GetHDF5Attribute(MLSFile, &
       & 'name', MLSFile%name)
      call GetHDF5Attribute(MLSFile, &
       & 'ShortName', MLSFile%ShortName)
      call GetHDF5Attribute(MLSFile, &
       & 'typeStr', MLSFile%typeStr)
      call GetHDF5Attribute(MLSFile, &
       & 'type', MLSFile%type)
      call GetHDF5Attribute(MLSFile, &
       & 'access', MLSFile%access)
      call GetHDF5Attribute(MLSFile, &
       & 'HDFVersion', MLSFile%HDFVersion)
      call GetHDF5Attribute(MLSFile, &
       & 'PCFID', MLSFile%PCFId)
      call GetHDF5Attribute(MLSFile, &
       & 'recordLength', MLSFile%recordLength)
      call GetHDF5Attribute(MLSFile, &
       & 'errorCode', MLSFile%errorCode)
      call GetHDF5Attribute( MLSFile, &
       & 'PCFIDRange', ints )
      MLSFile%PCFIdRange%Bottom = ints(1)
      MLSFile%PCFIdRange%Top    = ints(2)
      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   end subroutine h5_readMLSFileAttr
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine he5_readMLSFileAttr (MLSFile)
!------------------------------------------------------------

    use MLSHDFEOS, only: HE5_EHRDGlAtt
! Brief description of subroutine
! This subroutine reads the components of an MLSFile_t 
! as attributes from an hdfeos5-formatted file

! Arguments

      type(MLSFile_T)       :: MLSFile
! Local variables
      integer :: fileID, status
      integer, dimension(2) :: ints

      ! Executable code
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      fileID = MLSFile%FileID%f_id
      status = HE5_EHRDGLATT( fileID, &
       & 'content', &
       &  MLSFile%content )
      status = HE5_EHRDGLATT( fileID, &
       & 'lastOperation', &
       &  MLSFile%lastOperation )
      status = HE5_EHRDGLATT( fileID, &
       & 'name', &
       &  MLSFile%name )
      status = HE5_EHRDGLATT( fileID, &
       & 'ShortName', &
       &  MLSFile%ShortName )
      status = HE5_EHRDGLATT( fileID, &
       & 'typeStr', &
       &  MLSFile%typeStr )
      status = he5_EHrdglatt( fileID, &
       & 'type', MLSFile%type )
      status = he5_EHrdglatt( fileID, &
       & 'access', MLSFile%access )
      status = he5_EHrdglatt( fileID, &
       & 'HDFVersion', MLSFile%HDFVersion )
      status = he5_EHrdglatt( fileID, &
       & 'PCFID', MLSFile%PCFId )
      status = he5_EHrdglatt( fileID, &
       & 'recordlength', MLSFile%recordlength )
      status = he5_EHrdglatt( fileID, &
       & 'errorCode', MLSFile%errorCode )
      status = he5_EHrdglatt( fileID, &
       & 'PCFIDRange', ints )
      MLSFile%PCFIDRange%Bottom = ints(1)
      MLSFile%PCFIDRange%Top    = ints(2)

!------------------------------------------------------------
   end subroutine he5_readMLSFileAttr
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine h5_ReadGlobalAttr ( MLSFile )
!------------------------------------------------------------

      use HDF5, only: H5GClose_f, H5GOpen_f
      use MLSHDF5, only: GetHDF5Attribute, IsHDF5AttributePresent
! Brief description of subroutine
! Reads the global attributes from an hdf5-formatted file
! It finds them at the root '/' group level

! Arguments

      type(MLSFile_T)       :: MLSFile
! Local variables
      logical :: copiedFromAnEarlerFile
      integer :: fileID
      integer :: grp_id
      integer :: status
      integer, dimension(2) :: ints

      ! Executable code
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      fileID = MLSFile%FileID%f_id
      if ( IsHDF5AttributePresent( '/', fileID, 'InstrumentName') ) &
        & return
      copiedFromAnEarlerFile = &
        & IsHDF5AttributePresent( '/', fileID, 'FileAttributesCopiedFrom' )
      call h5gopen_f(fileID, '/', grp_id, status)
      MLSFile%FileID%grp_id = grp_id

      call GetHDF5Attribute( MLSFile, 'InstrumentName        ', GlobalAttributes%InstrumentName        )
      call GetHDF5Attribute( MLSFile, 'ProcessLevel          ', GlobalAttributes%ProcessLevel          )
      call GetHDF5Attribute( MLSFile, 'HostName              ', GlobalAttributes%HostName              )
      call GetHDF5Attribute( MLSFile, 'PGEVersion            ', GlobalAttributes%PGEVersion            )
      call GetHDF5Attribute( MLSFile, 'StartUTC              ', GlobalAttributes%StartUTC              )
      call GetHDF5Attribute( MLSFile, 'EndUTC                ', GlobalAttributes%EndUTC                )
      call GetHDF5Attribute( MLSFile, 'ProductionLoc         ', GlobalAttributes%productionLoc         )
      call GetHDF5Attribute( MLSFile, 'GranuleMonth          ', GlobalAttributes%GranuleMonth          )
      call GetHDF5Attribute( MLSFile, 'GranuleDay            ', GlobalAttributes%GranuleDay            )
      call GetHDF5Attribute( MLSFile, 'GranuleYear           ', GlobalAttributes%GranuleYear           )
      call GetHDF5Attribute( MLSFile, 'TAI93At0zOfGranule    ', GlobalAttributes%TAI93At0zOfGranule    )
      call GetHDF5Attribute( MLSFile, 'FirstMAFCtr           ', GlobalAttributes%FirstMAFCtr           )
      call GetHDF5Attribute( MLSFile, 'LastMAFCtr            ', GlobalAttributes%LastMAFCtr            )
      call GetHDF5Attribute( MLSFile, 'PhaseNames            ', GlobalAttributes%PhaseNames            )
      call GetHDF5Attribute( MLSFile, 'ForwardModelNames     ', GlobalAttributes%ForwardModelNames     )
      call GetHDF5Attribute( MLSFile, 'MiscNotes             ', GlobalAttributes%MiscNotes             )
      call GetHDF5Attribute( MLSFile, 'identifier_product_doi', GlobalAttributes%DOI                   )
!      call GetHDF5Attribute( MLSFile, 'ProductionLocation    ', GlobalAttributes%ProductionLocation    )
      call GetHDF5Attribute( MLSFile, 'NumCompletedChunks    ', GlobalAttributes%NumCompletedChunks    )
      call GetHDF5Attribute( MLSFile, 'NumFailedChunks       ', GlobalAttributes%NumFailedChunks       )
      call GetHDF5Attribute( MLSFile, 'FailedChunks          ', GlobalAttributes%FailedChunks          )
      call GetHDF5Attribute( MLSFile, 'FailedMachines        ', GlobalAttributes%FailedMachines        )
      call GetHDF5Attribute( MLSFile, 'FailedMsgs            ', GlobalAttributes%FailedMsgs            )

      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   end subroutine h5_ReadGlobalAttr
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine h5_writeglobalattr_fileID ( fileID, skip_if_already_there, DOI )
!------------------------------------------------------------

      use HDF5, only: H5GClose_f, H5GOpen_f
      use MLSHDF5, only: MakeHDF5Attribute, IsHDF5AttributePresent
      ! Brief description of subroutine
      ! This subroutine writes the global attributes for an hdf5-formatted file
      ! It does so at the root '/' group level

      ! Arguments

      integer, intent(in) :: fileID
      logical, intent(in), optional :: skip_if_already_there
      logical, intent(in), optional :: doi
      ! Local variables
      logical, parameter :: DeeBug = .false.
      integer :: grp_id
      integer :: status
      logical :: my_skip
      logical :: myDOI
      logical, parameter :: WRITE_ORBIT = .false.
      character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''

      ! Executable code
      if ( deebug ) then
        call output( 'Writing global attributes', advance='yes' )
        call dumpGlobalAttributes
      endif
      myDOI = len_trim(GlobalAttributes%FileAttributesCopiedFrom) > 0 ! .false.
      if ( present(DOI) ) myDOI=DOI
      my_skip = .false.
      if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
      if ( my_skip ) then
        if ( IsHDF5AttributePresent('/', fileID, 'InstrumentName') ) &
          & return
      endif
      call h5gopen_f(fileID, '/', grp_id, status)
      if (WRITE_ORBIT) then
          call MakeHDF5Attribute(grp_id, &
             & 'OrbitNumber', GlobalAttributes%OrbNum, .true.)
          call MakeHDF5Attribute(grp_id, &
             & 'OrbitPeriod', GlobalAttributes%OrbPeriod, .true.)
      end if
      call MakeHDF5Attribute(grp_id, &
        & 'InstrumentName', GlobalAttributes%InstrumentName, .true.)
      if ( len_trim(GlobalAttributes%HostName) > 0 ) then
        call MakeHDF5Attribute(grp_id, &
          & 'HostName', GlobalAttributes%HostName, .true.)
      else
        call MakeHDF5Attribute(grp_id, &
          & 'HostName', GlobalAttributes%ProductionLoc, .true.)
      endif
      ProcessLevel = ProcessLevelFun()
      call MakeHDF5Attribute(grp_id, &
       & 'ProcessLevel', ProcessLevel, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'PGEVersion', GlobalAttributes%PGEVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'StartUTC', GlobalAttributes%StartUTC, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'EndUTC', GlobalAttributes%EndUTC, .true.)
      if ( GlobalAttributes%GranuleDay < 1 &
        & .and. &
        & (len_trim(GlobalAttributes%FileAttributesCopiedFrom) < 1 ) ) return
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
      if ( lowercase(ProcessLevel(1:2)) == 'l2' ) then
        call MakeHDF5Attribute(grp_id, &
         & 'FirstMAF', GlobalAttributes%FirstMAFCtr, .true.)
        call MakeHDF5Attribute(grp_id, &
         & 'LastMAF', GlobalAttributes%LastMAFCtr, .true.)
      endif
      if ( len_trim(GlobalAttributes%FileAttributesCopiedFrom) > 0 ) then
        call MakeHDF5Attribute(grp_id, &
         & 'FileAttributesCopiedFrom', GlobalAttributes%FileAttributesCopiedFrom, .true.)
      endif
      ! We don't write these here (the master would insert an erroneous value)
      ! Instead we write them during DirectWrite operations
      ! call MakeHDF5Attribute(grp_id, &
      !  & 'MiscNotes', GlobalAttributes%MiscNotes, .true.)
      if ( len_trim(GlobalAttributes%DOI) > 0 .and. myDOI ) &
        & call MakeHDF5Attribute(grp_id, &
        & 'identifier_product_doi', GlobalAttributes%DOI, .false.)
      if ( len_trim(GlobalAttributes%productionLoc) > 0 .and. myDOI ) &
        & call MakeHDF5Attribute(grp_id, &
        & 'ProductionLocation', GlobalAttributes%productionLoc, .false.)
      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   end subroutine h5_writeglobalattr_fileID
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine h5_writeglobalattr_MLSFile ( MLSFile, skip_if_already_there, DOI )
!------------------------------------------------------------

      ! Brief description of subroutine
      ! This subroutine writes the global attributes for an hdf5-formatted file
      ! It does so at the root '/' group level

      ! Arguments

      type(MLSFile_T)       :: MLSFile
      logical, intent(in), optional :: skip_if_already_there
      logical, intent(in), optional :: doi
      ! Local variables
      logical :: alreadyOpen
      integer :: returnStatus
      ! Executable
      alreadyOpen = MLSFile%stillOpen
      if ( .not. alreadyOpen ) then
        call mls_openFile( MLSFile, returnStatus )
        if ( returnStatus /= 0 ) &
          call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'Unable to open hdf file', MLSFile=MLSFile )
      endif
      call h5_writeglobalattr_fileID ( MLSFile%fileID%f_id, &
        & skip_if_already_there, DOI )
      if ( .not. alreadyOpen ) call mls_closeFile( MLSFile, returnStatus )
   end subroutine h5_writeglobalattr_MLSFile

!------------------------------------------------------------
   subroutine h5_writeMLSFileAttr ( MLSFile, skip_if_already_there )
!------------------------------------------------------------

      use HDF5, only:  H5GClose_f, H5GOpen_f
      use MLSHDF5, only: ISHDF5AttributePresent, makeHDF5Attribute
      use HighOutput, only: BeVerbose
      ! Brief description of subroutine
      ! This subroutine writes the components of an MLSFile_t 
      ! as attributes for an hdf5-formatted file
      ! It does so at the root '/' group level

      ! Warning: Don't confuse this with h5_writeGlobalAttr
      ! which writes the run's global attributes (StartUTC, PGEVersion, etc.)

      ! Arguments

      type(MLSFile_T)       :: MLSFile
      logical, intent(in), optional :: skip_if_already_there
      ! Local variables
      integer :: fileID
      ! logical, parameter :: DeeBug = .false.
      logical :: DeeBug
      integer :: grp_id
      integer :: status
      logical :: my_skip

      ! Executable code
      Deebug = BeVerbose ( 'file', 0 )
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      if ( deebug ) then
        call output( 'Writing file attributes as global attributes', advance='yes' )
        call dumpGlobalAttributes
      endif
      my_skip = .false.
      if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
      fileID = MLSFile%FileID%f_id
      if ( my_skip ) then
        if ( IsHDF5AttributePresent('/', fileID, 'ShortName') ) &
          & return
      endif
      call h5gopen_f(fileID, '/', grp_id, status)
      call MakeHDF5Attribute(grp_id, &
       & 'content', MLSFile%content, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'lastOperation', MLSFile%lastOperation, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'name', MLSFile%name, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'ShortName', MLSFile%ShortName, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'typeStr', MLSFile%typeStr, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'type', MLSFile%type, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'access', MLSFile%access, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'HDFVersion', MLSFile%HDFVersion, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'PCFID', MLSFile%PCFId, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'recordLength', MLSFile%recordLength, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'errorCode', MLSFile%errorCode, .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'PCFIDRange', (/MLSFile%PCFIdRange%Bottom, MLSFile%PCFIdRange%Top/), .true.)
      call MakeHDF5Attribute(grp_id, &
       & 'FileID', &
       & (/MLSFile%FileID%f_id, MLSFile%FileID%grp_id, MLSFile%FileID%sd_id/), &
       & .true.)
      call h5gclose_f(grp_id, status)

!------------------------------------------------------------
   end subroutine h5_writeMLSFileAttr
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine he5_writeglobalattr_MLSFile ( MLSFile, dayNum, DOI, skip_if_already_there )
!------------------------------------------------------------

! Brief description of subroutine
! This subroutine writes the global attributes for an hdfeos5 file

! Arguments

      type(MLSFile_T)       :: MLSFile
      integer, intent(in), optional :: dayNum
      logical, intent(in), optional :: doi
      logical, intent(in), optional :: skip_if_already_there
      ! Local variables
      logical :: alreadyOpen
      integer :: returnStatus
      ! Executable
      alreadyOpen = MLSFile%stillOpen
      if ( .not. alreadyOpen ) then
        call mls_openFile( MLSFile, returnStatus )
        if ( returnStatus /= 0 ) &
          call MLSMessage( MLSMSG_Error, ModuleName, &
          & 'Unable to open hdfeos file', MLSFile=MLSFile )
      endif
      call he5_writeglobalattr_FileID ( MLSFile%fileID%f_id, &
        & dayNum, DOI, skip_if_already_there )
      if ( .not. alreadyOpen ) call mls_closeFile( MLSFile, returnStatus )
   end subroutine he5_writeglobalattr_MLSFile

!------------------------------------------------------------
   subroutine he5_writeglobalattr_FileID ( fileID, dayNum, DOI, skip_if_already_there )
!------------------------------------------------------------

    use HDFEOS5, only: HE5T_Native_Int, &
      & HE5T_Native_Double, MLS_CharType
    use MLSHDFEOS, only: HE5_EHWRGlAtt, HSize, MLS_EhwrGlAtt, MLS_IsGlAtt
! Brief description of subroutine
! This subroutine writes the global attributes for an hdfeos5 file

! Arguments

      integer, intent(in) :: fileID
      integer, intent(in), optional :: dayNum
      logical, intent(in), optional :: doi
      logical, intent(in), optional :: skip_if_already_there
! Internal variables
      integer :: status
      character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''
      logical :: myDOI
      logical, parameter :: DeeBug = .false.
      logical :: my_skip
! Executable
      myDOI = .false.
      if ( present(DOI) ) myDOI=DOI
      my_skip = .false.
      if ( present(skip_if_already_there) ) my_skip=skip_if_already_there
      if ( deebug ) then
        call output( 'Writing global attributes', advance='yes' )
        call dumpGlobalAttributes
      endif
      if ( len_trim(GlobalAttributes%DOI) > 0 .and. myDOI ) &
       & status = mls_EHwrglatt(fileID, &
       & 'identifier_product_doi', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%DOI)
      if ( len_trim(GlobalAttributes%productionLoc) > 0 .and. myDOI ) &
       & status = mls_EHwrglatt(fileID, &
       & 'ProductionLocation', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%productionLoc)
      if ( my_skip ) then
        if ( MLS_IsGlAtt ( fileID, 'OrbitNumber' ) ) &
          & return
      endif
      if (present(dayNum)) then
         status = he5_EHwrglatt(fileID, &
            & 'OrbitNumber', HE5T_NATIVE_INT, hsize(max_orbits), &
            &  GlobalAttributes%OrbNumDays(:,dayNum))
         status = he5_EHwrglatt(fileID, &
            & 'OrbitPeriod', HE5T_NATIVE_DOUBLE, hsize(max_orbits), &
            &  GlobalAttributes%OrbPeriodDays(:,dayNum))
      else
         status = he5_EHwrglatt(fileID, &
            & 'OrbitNumber', HE5T_NATIVE_INT, hsize(max_orbits), &
            &  GlobalAttributes%OrbNum)
         status = he5_EHwrglatt(fileID, &
            & 'OrbitPeriod', HE5T_NATIVE_DOUBLE, hsize(max_orbits), &
            &  GlobalAttributes%OrbPeriod)
      end if
      status = mls_EHwrglatt(fileID, &
       & 'InstrumentName', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%InstrumentName)
      if ( len_trim(GlobalAttributes%HostName) > 0 ) then
        status = mls_EHwrglatt(fileID, &
         & 'HostName', MLS_CHARTYPE, 1, &
         &  GlobalAttributes%HostName)
      else
        status = mls_EHwrglatt(fileID, &
         & 'HostName', MLS_CHARTYPE, 1, &
         &  GlobalAttributes%ProductionLoc)
      endif
      ProcessLevel = ProcessLevelFun()
      status = mls_EHwrglatt(fileID, &
       & 'ProcessLevel', MLS_CHARTYPE, 1, &
       &  ProcessLevel)
      status = mls_EHwrglatt(fileID, &
       & 'PGEVersion', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%PGEVersion)
      status = mls_EHwrglatt(fileID, &
       & 'StartUTC', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%StartUTC)
      status = mls_EHwrglatt(fileID, &
       & 'EndUTC', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%EndUTC)
      ! if ( GlobalAttributes%GranuleDay == ' ') return
      ! if ( GlobalAttributes%GranuleMonth == ' ') &
      if ( GlobalAttributes%GranuleDay < 1 ) return
      status = he5_EHwrglatt(fileID, &
       & 'GranuleMonth', HE5T_NATIVE_INT, hsize(1), &
       &  (/ GranuleMonth_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleDay', HE5T_NATIVE_INT, hsize(1), &
       &  (/ GranuleDay_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleDayOfYear', HE5T_NATIVE_INT, hsize(1), &
       &  (/ GranuleDayOfYear_fun() /) )
      status = he5_EHwrglatt(fileID, &
       & 'GranuleYear', HE5T_NATIVE_INT, hsize(1), &
       &  (/ GlobalAttributes%GranuleYear/) )
      status = he5_EHwrglatt(fileID, &
       & 'TAI93At0zOfGranule', HE5T_NATIVE_DOUBLE, hsize(1), &
       &  (/ GlobalAttributes%TAI93At0zOfGranule/) )
      if ( lowercase(ProcessLevel(1:2)) == 'l2' ) then
        status = he5_EHwrglatt(fileID, &
         & 'FirstMAF', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%FirstMAFCtr /) )
        status = he5_EHwrglatt(fileID, &
         & 'LastMAF', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%LastMAFCtr /) )
      endif
      ! if ( DEBUG ) call outputNamedValue( 'GlobalAttributes%MiscNotes: ', GlobalAttributes%MiscNotes )
      ! We don't write these here (the master would insert an erroneous value)
      ! Instead we write them during DirectWrite operations
      ! status = mls_EHwrglatt(fileID, &
      ! & 'MiscNotes', MLS_CHARTYPE, 1, &
      !  &  GlobalAttributes%MiscNotes)
!------------------------------------------------------------
   end subroutine he5_writeglobalattr_FileID
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine he5_writeMLSFileAttr ( MLSFile )
!------------------------------------------------------------

   use HDFeos5, only: He5t_Native_Int, MLS_Chartype
   use MLSHDFeos, only: He5_Ehwrglatt, Hsize, MLS_Ehwrglatt
! Brief description of subroutine
! This subroutine writes the components of an MLSFile_t 
! as attributes for an hdfeos5-formatted file

! Arguments

      ! integer, intent(in) :: fileID
      type(MLSFile_T)       :: MLSFile
! Local variables
      integer :: fileID
      integer :: status

      ! Executable code
      if ( .not. MLSFile%stillOpen ) then
        call MLS_OpenFile( MLSFile )
      endif
      fileID = MLSFile%FileID%f_id
      status = mls_EHwrglatt(fileID, &
       & 'content', MLS_CHARTYPE, 1, &
       &  MLSFile%content)
      status = mls_EHwrglatt(fileID, &
       & 'lastOperation', MLS_CHARTYPE, 1, &
       &  MLSFile%lastOperation)
      status = mls_EHwrglatt(fileID, &
       & 'name', MLS_CHARTYPE, 1, &
       &  MLSFile%name)
      status = mls_EHwrglatt(fileID, &
       & 'ShortName', MLS_CHARTYPE, 1, &
       &  MLSFile%ShortName)
      status = mls_EHwrglatt(fileID, &
       & 'typeStr', MLS_CHARTYPE, 1, &
       &  MLSFile%typeStr)
      status = he5_EHwrglatt(fileID, &
       & 'type', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%type /) )
      status = he5_EHwrglatt(fileID, &
       & 'access', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%access /) )
      status = he5_EHwrglatt(fileID, &
       & 'HDFVersion', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%HDFVersion /) )
      status = he5_EHwrglatt(fileID, &
       & 'PCFID', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%PCFId /) )
      status = he5_EHwrglatt(fileID, &
       & 'recordlength', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%recordlength /) )
      status = he5_EHwrglatt(fileID, &
       & 'errorCode', HE5T_NATIVE_INT, hsize(1), &
       &  (/ MLSFile%errorCode /) )
      status = he5_EHwrglatt(fileID, &
       & 'PCFIDRange', HE5T_NATIVE_INT, hsize(2), &
       &  (/ MLSFile%PCFIDRange%Bottom, MLSFile%PCFIDRange%Top /) )
      status = he5_EHwrglatt(fileID, &
       & 'FileID', HE5T_NATIVE_INT, hsize(3), &
       &  (/MLSFile%FileID%f_id, MLSFile%FileID%grp_id, MLSFile%FileID%sd_id/) )

!------------------------------------------------------------
   end subroutine he5_writeMLSFileAttr
!------------------------------------------------------------

!------------------------------------------------------------
   subroutine he5_readglobalattr ( fileID, gAttributes, &
     & ProcessLevel, DayofYear, TAI93At0zOfGranule, &
     & HostName, MiscNotes, DOI, returnStatus )
!------------------------------------------------------------

     use MLSHDFeos, only: Maxdlistlength, He5_Ehrdglatt, MLS_IsGlatt
! Brief description of subroutine
! This subroutine reads the global attributes from an hdf-eos5 file

! Arguments

      integer, intent(in)                      :: fileID
      type (GlobalAttributes_T), intent(inout) :: gAttributes
      character(len=*), intent(out), optional  :: ProcessLevel
      integer, intent(out), optional           :: DayofYear
      double precision, optional, intent(out)  :: TAI93At0zOfGranule
      character(len=*), intent(out), optional  :: HostName
      character(len=*), intent(out), optional  :: MiscNotes
      character(len=*), intent(out), optional  :: DOI
      integer, optional, intent(out)           :: returnStatus
! Internal variables
      integer :: status
      integer, dimension(1) :: ibuf
      real(r8), dimension(1) :: dbuf
      ! character(len=MAXDLISTLENGTH) :: attrList
      ! integer :: listSize
! Executable
      if ( DEBUG ) call output( 'Reading global attributes', advance='yes' )
      !status = he5_EHinqglatts(fileID, attrList, listSize)
      !status = 0
      ! if ( status /= 0 ) then
      if ( present(returnStatus) ) returnStatus = 1
      if ( .not. MLS_IsGlatt ( FileID, 'OrbitNumber') ) then
        call MLSMessage(MLSMSG_Warning, ModuleName, &
            & 'No global attributes in file: ' )
        return
      endif
      status = he5_EHrdglatt(fileID, &
         & 'OrbitNumber', &
         &  gAttributes%OrbNum )
      status = he5_EHrdglatt(fileID, &
         & 'OrbitPeriod', &
         &  gAttributes%OrbPeriod )
      status = he5_EHrdglatt(fileID, &
       & 'InstrumentName', &
       &  gAttributes%InstrumentName)
      if ( present(HostName) ) then
       status = he5_EHrdglatt(fileID, &
        & 'HostName', &
        &  HostName)
      endif
      status = he5_EHrdglatt(fileID, &
       & 'ProcessLevel', &
       &  gattributes%ProcessLevel)
      if ( present(ProcessLevel) ) ProcessLevel = gattributes%ProcessLevel
      status = he5_EHrdglatt(fileID, &
       & 'PGEVersion', &
       &  gAttributes%PGEVersion)
      if ( present(MiscNotes) ) then
       status = he5_EHrdglatt(fileID, &
       & 'MiscNotes', &
       &  MiscNotes)
       if ( DEBUG ) call outputNamedValue('Misc Notes (read) ', trim(MiscNotes) )
      endif
      status = he5_EHrdglatt(fileID, &
       & 'StartUTC', &
       &  gAttributes%StartUTC)
      status = he5_EHrdglatt(fileID, &
       & 'EndUTC', &
       &  gAttributes%EndUTC)
      status = he5_EHrdglatt(fileID, &
       & 'GranuleMonth', &
       &  ibuf  )
      gAttributes%GranuleMonth = ibuf(1)
      status = he5_EHrdglatt(fileID, &
       & 'GranuleDay', &
       &  ibuf  )
      gAttributes%GranuleDay = ibuf(1)
      status = he5_EHrdglatt(fileID, &
       & 'GranuleDayOfYear', &
       &  ibuf  )
      if ( present(DayofYear) ) DayOfYear=ibuf(1)
      status = he5_EHrdglatt(fileID, &
       & 'GranuleYear', &
       &  ibuf  )
      gAttributes%GranuleYear = ibuf(1)
      if ( present(TAI93At0zOfGranule) ) then
        status = he5_EHrdglatt(fileID, &
         & 'TAI93At0zOfGranule', dbuf )
        TAI93At0zOfGranule = dbuf(1)
      endif
      if ( lowercase(gAttributes%ProcessLevel(1:2)) == 'l2' ) then
        status = he5_EHrdglatt(fileID, &
         & 'FirstMAF', &
         &  ibuf  )
        gAttributes%FirstMAFCtr = ibuf(1)
        status = he5_EHrdglatt(fileID, &
         & 'LastMAF', &
         &  ibuf  )
        gAttributes%LastMAFCtr = ibuf(1)
      endif
      if ( present(DOI) ) then
       status = he5_EHrdglatt(fileID, &
       & 'identifier_product_doi', &
       &  DOI)
       if ( DEBUG ) call outputNamedValue('identifier_product_doi ', trim(DOI) )
      endif
      ! if ( present(returnStatus) ) returnStatus = status
      if ( DEBUG ) call dumpGlobalAttributes
      if ( present(returnStatus) ) returnStatus = 0
!------------------------------------------------------------
   end subroutine he5_readglobalattr
!------------------------------------------------------------

!----------------------------------------
   subroutine InputInputpointer (urefs, fileIDArray, fileNameArray, &
     & PCBottom, PCTop )
!----------------------------------------
!  Prepare Input for WriteInputpointer consisting of universal refs
!  This can be done for an array of fileids or of file names or both
!  that were input during the run; e.g., l1brads for level 2
!  If no universal refs are found, as always happens with my test cases,
!  then put the file names as returned by pgs_pc_getReference in their place
!  for the fileIDArray, or else the fileNameArray itself if appropriate
!  Arguments
    character (len=INPUTPTR_STRING_LENGTH), intent(out)    :: urefs(:)
    character (len=*), dimension(:), intent(in), optional  :: fileNameArray
    integer, dimension(:), intent(in), optional            :: fileIDArray
    integer,  intent(in), optional                 :: PCBottom, PCTop
!  Local variables
   integer  :: i
   integer  :: ref
   integer  :: returnStatus
   integer  :: thePC
   integer  :: version
   character(len=INPUTPTR_STRING_LENGTH) :: sval
   integer, external :: pgs_pc_getUniversalRef, pgs_pc_getReference

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
   end subroutine InputInputpointer
!------------------------------------

!------------------------------------------------------------
  subroutine SomeToGlobalAttributes
    ! Copy SomeGlobalAttribues into the whole enchilada
    ! (But only where Some is non-blank)
    call CopyIfNonBlank ( SomeGlobalAttributes%InstrumentName  , GlobalAttributes%InstrumentName  )
    call CopyIfNonBlank ( SomeGlobalAttributes%ProcessLevel    , GlobalAttributes%ProcessLevel    )
    call CopyIfNonBlank ( SomeGlobalAttributes%HostName        , GlobalAttributes%HostName        )
    call CopyIfNonBlank ( SomeGlobalAttributes%PGEVersion      , GlobalAttributes%PGEVersion      )
    call CopyIfNonBlank ( SomeGlobalAttributes%MiscNotes       , GlobalAttributes%MiscNotes       )
    call CopyIfNonBlank ( SomeGlobalAttributes%StartUTC        , GlobalAttributes%StartUTC        )
    call CopyIfNonBlank ( SomeGlobalAttributes%EndUTC          , GlobalAttributes%EndUTC          )
    call CopyIfNonBlank ( SomeGlobalAttributes%DOI             , GlobalAttributes%DOI             )
    call CopyIfNonBlank ( SomeGlobalAttributes%productionLoc   , GlobalAttributes%productionLoc   )
  contains
    subroutine CopyIfNonBlank( FromMe, ToMe )
      character(len=*), intent(in)    :: FromMe
      character(len=*), intent(inout) :: ToMe
      !
      if ( len_trim(FromMe) > 0 ) ToMe = FromMe
    end subroutine CopyIfNonBlank
  end subroutine SomeToGlobalAttributes

!------------------------------------------------------------
   subroutine sw_writeglobalattr (swathID)
!------------------------------------------------------------

     use HDFeos5, only: He5t_Native_Int, He5t_Native_Double, MLS_Chartype
     use He5_Swapi, only: He5_Swwrattr
     use MLSHDFeos, only: Hsize, MLS_Swwrattr
! Brief description of subroutine
! This subroutine writes the global attributes for an hdfeos5 swath

! Arguments

      integer, intent(in) :: swathID
!     integer, external ::   he5_SWwrattr
! Internal variables
      logical, parameter :: DeeBug = .false.
      integer :: status
      character(len=GA_VALUE_LENGTH) :: ProcessLevel = ''
! Executable
      if ( deeBug ) then
        call output( 'Writing global attributes', advance='yes' )
        call dumpGlobalAttributes
      endif
      status = he5_SWwrattr(swathID, &
       & 'OrbitNumber', HE5T_NATIVE_INT, hsize(max_orbits), &
       &  GlobalAttributes%OrbNum)
      status = he5_SWwrattr(swathID, &
       & 'OrbitPeriod', HE5T_NATIVE_DOUBLE, hsize(max_orbits), &
       &  GlobalAttributes%OrbPeriod)
      status = he5_SWwrattr(swathID, &
       & 'InstrumentName', MLS_CHARTYPE, hsize(1), &
       &  GlobalAttributes%InstrumentName)
      if ( len_trim(GlobalAttributes%HostName) > 0 ) then
        status = he5_SWwrattr(swathID, &
         & 'HostName', MLS_CHARTYPE, hsize(1), &
         &  GlobalAttributes%HostName)
      else
        status = he5_SWwrattr(swathID, &
         & 'HostName', MLS_CHARTYPE, hsize(1), &
         &  GlobalAttributes%ProductionLoc)
      endif
      ProcessLevel = ProcessLevelFun()
      status = mls_SWwrattr(swathID, &
       & 'ProcessLevel', MLS_CHARTYPE, 1, &
       &  ProcessLevel)
      status = mls_SWwrattr(swathID, &
       & 'PGEVersion', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%PGEVersion)
      status = mls_SWwrattr(swathID, &
       & 'StartUTC', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%StartUTC)
      status = mls_SWwrattr(swathID, &
       & 'EndUTC', MLS_CHARTYPE, 1, &
       &  GlobalAttributes%EndUTC)
      if ( GlobalAttributes%GranuleDay < 1 ) return
      if ( GlobalAttributes%GranuleMonth > 0 ) then
        status = he5_SWwrattr(swathID, &
         & 'GranuleMonth', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%GranuleMonth/) )
        status = he5_SWwrattr(swathID, &
         & 'GranuleDay', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%GranuleDay/) )
      else
        status = he5_SWwrattr(swathID, &
         & 'GranuleDayOfYear', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%GranuleDay/) )
      endif
      status = he5_SWwrattr(swathID, &
       & 'GranuleYear', HE5T_NATIVE_INT, hsize(1), &
       &  (/ GlobalAttributes%GranuleYear/) )
      status = he5_SWwrattr(swathID, &
       & 'TAI93At0zOfGranule', HE5T_NATIVE_DOUBLE, hsize(1), &
       &  (/ GlobalAttributes%TAI93At0zOfGranule/) )
      if ( lowercase(ProcessLevel(1:2)) == 'l2' ) then
        status = he5_SWwrattr(swathID, &
         & 'FirstMAF', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%FirstMAFCtr/) )
        status = he5_SWwrattr(swathID, &
         & 'LastMAF', HE5T_NATIVE_INT, hsize(1), &
         &  (/ GlobalAttributes%LastMAFCtr/) )
      endif
      ! We don't write these here (the master would insert an erroneous value)
      ! Instead we write them during DirectWrite operations
      ! status = he5_SWwrattr(swathID, &
      ! & 'MiscNotes', MLS_CHARTYPE, hsize(1), &
      !  &  GlobalAttributes%MiscNotes)
!------------------------------------------------------------
   end subroutine sw_writeglobalattr
!------------------------------------------------------------

!----------------------------------------
   function WriteInputpointer (groups, attrName, inpt, fileType)
!----------------------------------------

!  Write Inputpointer metadata

!  Arguments
    character (len = PGSd_MET_GROUP_NAME_L) :: Groups
    character (len=*), intent(in) :: Attrname
    character (len=INPUTPTR_STRING_LENGTH), intent(in), optional  :: inpt(:)
    ! character(len=*), intent(in), optional :: fileType   ! 'hdfeos', 'hdf', 'sw' or ..
    integer, intent(in), optional :: fileType   ! l_swath, l_hdf, ...

    integer             :: WriteInputpointer
    integer, external   :: pgs_met_setAttr_s
    ! character (len=6) :: the_type
    integer :: the_type

!   Executable statements
    if ( present(inpt) ) then
       WriteInputpointer = pgs_met_setAttr_s(groups, attrName, inpt)
       return
    endif
    the_type = l_hdf
    if ( present(fileType) ) the_type = fileType
    select case(the_type)
    case (l_hdf)
      WriteInputpointer = pgs_met_setAttr_s (groups, attrName, &
        &  (/HDFINPTPTRVALUE/) )
    case default
      WriteInputpointer = pgs_met_setAttr_s (groups, attrName, &
        &  (/HDFEOSINPTPTRVALUE/) )
    end select      
    
!------------------------------------
   end function WriteInputpointer
!------------------------------------

   subroutine WriteLeapSecHDFEOSAttr (fileID)
     ! Write contents of leapsec file as hdfeos5 attribute to file
    use MLSHDFEOS, only: MLS_EHWrGlAtt
     ! Args
     integer, intent(in) :: fileID
     ! Internal variables
     character(len=FileNameLen) :: LeapSecFile
     integer, parameter :: PCFid = 10301
     integer :: status
     integer :: version
     ! Executable
     version = 1
     Status = Pgs_pc_getReference( PCFid, version, &
       & LeapSecFile )
     if ( Status /= PGS_S_SUCCESS ) then
       call outputNamedValue( 'status', status )
       CALL MLSMessage( MLSMSG_Warning, ModuleName, &
         & 'Unable to get path, file name for Leap sec file using PCFid' )
       return
     end if
     status = MLS_EHWRGLATT ( trim(leapSecFile), FILEID, &
       & 'leap seconds' )
   end subroutine WriteLeapSecHDFEOSAttr

   subroutine WriteLeapSecHDF5DS (fileID)
     ! Write contents of leapsec file as hdf5 dataset to file
    use MLSHDF5, only: SaveAsHDF5DS
     ! Args
     integer, intent(in) :: fileID
     ! Internal variables
     character(len=FileNameLen) :: LeapSecFile
     integer, parameter :: PCFid = 10301
     integer :: status
     integer :: version
     ! Executable
     version = 1
     Status = Pgs_pc_getReference( PCFid, version, &
       & LeapSecFile )
     if ( Status /= PGS_S_SUCCESS ) then
       call outputNamedValue( 'status', status )
       CALL MLSMessage( MLSMSG_Warning, ModuleName, &
         & 'Unable to get path, file name for Leap sec file using PCFid' )
       return
     end if
     call SaveAsHDF5DS ( trim(leapSecFile), FILEID, &
       & 'leap seconds' )
     
   end subroutine WriteLeapSecHDF5DS

!----------------------------------------
   subroutine WritePCF2Hdr (file, anText, hdfVersion, fileType, name)
!----------------------------------------
      use HDF5, only: H5F_AcC_Rdwr_F, &
        & H5FOpen_F, H5FClosE_F

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS file as an annotation.

! Arguments

      character (len=FileNameLen), intent(in) :: file 

      character (len=1), pointer              :: anText(:)
      integer, intent(in), optional           :: hdfVersion
      ! character(len=*), intent(in), optional  :: fileType ! 'sw', 'gd', 'hdf'
      integer, intent(in), optional  :: fileType ! l_swath, l_hdf, ..
      character(len=*), intent(in), optional :: name
      ! logical, intent(in), optional         :: isHDFEOS

! Parameters
      integer :: fileID
      type(MLSFile_T)                :: MLSFile
      integer :: my_hdfVersion
      ! logical :: myisHDFEOS
    ! integer :: record_length
      integer :: status
      ! character (len=2) :: the_type
      integer :: the_type
! Executable
      my_hdfVersion = PCFHDR_DEFAULT_HDFVERSION
      if ( present(hdfVersion) ) my_hdfVersion = hdfVersion
      the_type = l_hdf
      if ( present(fileType) ) the_type = fileType
      status = InitializeMLSFile ( MLSFile, type=the_Type, access=DFACC_RDWR, &
       & name=trim(file), HDFVersion=my_hdfVersion )
      select case(my_hdfVersion)
      case (HDFVERSION_4)
        call WritePCF2Hdr_hdf4 (file, anText)
      case (HDFVERSION_5)
        if ( the_type == l_swath ) then
          ! fileID = mls_io_gen_openF(l_swath, .TRUE., status, &
          ! & record_length, DFACC_RDWR, FileName=trim(file), &
          ! & hdfVersion=hdfVersion, debugOption=.false. )
          call MLS_OpenFile( MLSFile )
          fileID = MLSFile%FileID%f_id
          ! if ( status /= PGS_S_SUCCESS) &
          !  & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          !  & 'Error opening hdfeos5 swath file for annotating with PCF' )
          call WritePCF2Hdr_hdfeos5 (fileID, anText)
          ! status = mls_io_gen_closeF(l_swath, fileID, &
          !  & hdfVersion=hdfVersion)
          call MLS_CloseFile( MLSFile )
          ! if ( status /= PGS_S_SUCCESS) &
          !  & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          !  & 'Error closing hdfeos5 swath file for annotating with PCF' )
        elseif ( the_type == l_hdfeos ) then
          ! fileID = mls_io_gen_openF(l_grid, .TRUE., status, &
          !  & record_length, DFACC_RDWR, FileName=trim(file), &
          !  & hdfVersion=hdfVersion, debugOption=.false. )
          call MLS_OpenFile( MLSFile )
          fileID = MLSFile%FileID%f_id
          ! if ( status /= PGS_S_SUCCESS) &
          !   & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          !   & 'Error opening hdfeos5 grid file for annotating with PCF' )
          call WritePCF2Hdr_hdfeos5 (fileID, anText)
          ! status = mls_io_gen_closeF(l_grid, fileID, &
          !  & hdfVersion=hdfVersion)
          call MLS_CloseFile ( MLSFile )
          ! if ( status /= PGS_S_SUCCESS) &
          !  & CALL MLSMessage(MLSMSG_Error, ModuleName, &
          !  & 'Error closing hdfeos5 grid file for annotating with PCF' )
        else
          call h5fopen_f(trim(file), H5F_ACC_RDWR_F, fileID, status)
          if ( status /= PGS_S_SUCCESS) &
            & CALL MLSMessage(MLSMSG_Error, ModuleName, &
            & 'Error opening hdf5 file for annotating with PCF' )
          call WritePCF2Hdr_hdf5 (fileID, anText, name)
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
   end subroutine WritePCF2Hdr
!-----------------------------

!----------------------------------------
   subroutine WritePCF2Hdr_hdf4 (file, anText)
!----------------------------------------

! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS4 file as an annotation.

! Arguments

      character (len=FileNameLen), intent(in) :: file 

      character (len=1), pointer              :: anText(:)

! Parameters

! Functions

      integer, external :: afEnd, afEndAccess, afFCreate, afStart, afWriteAnn
      integer, external :: hClose, hOpen

! Variables

      character (len=480) :: msr

      integer :: anID, annID, fileID, status

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
   end subroutine WritePCF2Hdr_hdf4
!-----------------------------

!----------------------------------------
   subroutine WritePCF2Hdr_hdf5 (fileID, anText, name)
!----------------------------------------

     use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
     use HDF5, only: H5gclose_F, H5gopen_F
     use MLSHDF5, only: IsHDF5dspresent, MakeHDF5attribute, SaveasHDF5ds
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

      integer :: fileID
      character (len=1), pointer              :: anText(:)
      character(len=*), intent(in), optional :: name

! Local variables
      integer :: grp_id
      integer :: status
      character (len=1), dimension(:), pointer              :: an40
      integer :: how_big
      logical, parameter :: useLENGTHONECHARS = .true.
      logical, parameter :: MAKEDATASET = .true.
      logical, parameter :: MAKEATTRIBUTE = .not. MAKEDATASET
      character(len=size(anText)) :: anScalar
      character(len=80) :: myPCFPATHNAME
      ! Executable
      myPCFPATHNAME = PCFPATHNAME
      if ( present(name) ) myPCFPATHNAME = name
      if ( MAKEDATASET ) then
        if ( IsHDF5DSPresent(fileID, trim(myPCFPATHNAME) ) ) return
        anScalar = transfer(anText, anScalar)
        call SaveAsHDF5DS ( fileID, trim(myPCFPATHNAME), anScalar )
      endif
      if ( .not. MAKEATTRIBUTE ) return
      myPCFPATHNAME = PCFATTRIBUTENAME
      if ( present(name) ) myPCFPATHNAME = name
      call h5gopen_f(fileID, '/', grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error opening hdf5 file root group for annotating with PCF' )
      if ( useLENGTHONECHARS ) then
        call MakeHDF5Attribute(grp_id, &
         & trim(myPCFPATHNAME), anText, .true.)
         ! & 'PCF file text', anText, .true.)
      else
        ! Find how big an40 must be to hold anText
        how_big = 1 + (size(anText)-1)/40
        call allocate_test ( an40, how_big, "an40", moduleName )
        an40 = ' '
        ! Do some nonsense here
        call MakeHDF5Attribute(grp_id, &
         & trim(myPCFPATHNAME), an40, .true.)
        call deallocate_test ( an40, "an40", moduleName )
      end if
      call h5gclose_f(grp_id, status)
      if ( status /= PGS_S_SUCCESS) &
        & CALL MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Error closing hdf5 file root group for annotating with PCF' )
!-----------------------------
   end subroutine WritePCF2Hdr_hdf5
!-----------------------------

!----------------------------------------
   subroutine WritePCF2Hdr_hdfeos5 (fileID, anText)
!----------------------------------------

     use HDFeos5, only: MLS_Chartype
     use MLSHDFeos, only: Hsize, He5_Ehwrglatt
! Brief description of subroutine
! This subroutine writes the PCF into an HDF-EOS5 file as an attribute.
! It does so as file level attributes
! It does not do so as a dataset because that would create a hybrid file
! For unclear reasons the attributes are stored as an array of 40-length chars
! Unless the parameter useLENGTHONECHARS is true
! Arguments

! If the PCF file if bigger than 40,000 characters, need to break it down into 
! multiple attribute names such as PCF1, PCF2, etc... up to PCF9
! Do this to solve the probem of L3 PCF files 

      character (len=1), pointer              :: anText(:)
      integer, intent(in)                     :: fileID

! Local variables
      integer :: firstChar, lastChar, blockNumber
      integer :: status
      integer, parameter :: maxheadersize = 40000
      character (len=3) :: blockChar 

! Executable

      blockNumber = 1
      firstChar = 1

      do
        lastChar = min(firstChar-1+maxheadersize, size(anText))
        ! i1 is for PCF filesize < 400,000 
       write( blockChar, '(i1)') blockNumber
       if ( blockNumber > 9 ) write( blockChar, '(i2)') blockNumber

        status = he5_EHwrglatt(fileID, 'PCF'//TRIM(ADJUSTL(blockChar)), &
                & MLS_CHARTYPE, hsize(lastChar-firstChar+1), &
                & anText(firstChar:lastChar))
        if ( status /= PGS_S_SUCCESS) then
           CALL MLSMessage(MLSMSG_Error, ModuleName, &
                & 'Error annotating with PCF')
           return
        end if
        blockNumber = blockNumber + 1
        firstChar = firstChar + maxheadersize
        if ( firstChar > size(anText)) return
      end do

!-----------------------------
   end subroutine WritePCF2Hdr_hdfeos5
!-----------------------------

   subroutine WriteutcPoleHDFEOSAttr (fileID)
     ! Write contents of utcPole file as hdfeos5 attribute to file
    use MLSHDFEOS, only: MLS_EHWrGlAtt
     ! Args
     integer, intent(in) :: fileID
     ! Internal variables
     character(len=FileNameLen) :: utcPoleFile
     integer, parameter :: PCFid = 10401
     integer :: status
     integer :: version
     ! Executable
     version = 1
     Status = Pgs_pc_getReference( PCFid, version, &
       & utcPoleFile )
     if ( Status /= PGS_S_SUCCESS ) then
       call outputNamedValue( 'status', status )
       CALL MLSMessage( MLSMSG_Warning, ModuleName, &
         & 'Unable to get path, file name for utc pole file using PCFid' )
       return
     end if
     status = MLS_EHWRGLATT ( trim(utcPoleFile), FILEID, &
       & 'utc pole' )
   end subroutine WriteutcPoleHDFEOSAttr

   subroutine WriteutcPoleHDF5DS (fileID)
     ! Write contents of leapsec file as hdf5 dataset to file
    use MLSHDF5, only: SaveAsHDF5DS
     ! Args
     integer, intent(in) :: fileID
     ! Internal variables
     character(len=FileNameLen) :: utcPoleFile
     integer, parameter :: PCFid = 10401
     integer :: status
     integer :: version
     ! Executable
     version = 1
     Status = Pgs_pc_getReference( PCFid, version, &
       & utcPoleFile )
     if ( Status /= PGS_S_SUCCESS ) then
       call outputNamedValue( 'status', status )
       CALL MLSMessage( MLSMSG_Warning, ModuleName, &
         & 'Unable to get path, file name for utc pole file using PCFid' )
       return
     end if
     call SaveAsHDF5DS ( trim(utcPoleFile), FILEID, &
       & 'utc pole' )
     
   end subroutine WriteutcPoleHDF5DS

!----------- Internal Procedures --------------
  function int_to_char (int) result (chars)
    ! Arguments
    integer , intent(in) :: int
    character(len=16) :: chars
    ! Executable
    write(chars, * ) int
    chars = adjustl(chars)
  end function int_to_char

  function ProcessLevelFun () result (ProcessLevel)
    ! Take '1', '2', '3' and return 'L1', 'L2', etc.
    ! Leaves 'L*' unchanged
    ! Arguments
    character(len=GA_VALUE_LENGTH) :: ProcessLevel
    ! Executable
    if ( GlobalAttributes%ProcessLevel == ' ' ) then
      ProcessLevel = DEFAULTPROCESSLEVEL   ! 'unknown'
    elseif ( GlobalAttributes%ProcessLevel(1:1) /= 'L' ) then
      ProcessLevel = 'L' // trim(GlobalAttributes%ProcessLevel)
    else
      ProcessLevel = GlobalAttributes%ProcessLevel
    endif
  end function ProcessLevelFun

!================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module PCFHdr
!================

!# $Log$
!# Revision 2.79  2022/12/20 23:32:07  pwagner
!# Added two more file attributes
!#
!# Revision 2.78  2020/03/20 23:02:10  pwagner
!# Attempt to relieve confusion between ProductionLocation and ProductionLoc; may optionally read DOI, etc.
!#
!# Revision 2.77  2020/03/04 21:22:20  pwagner
!# Make ProcessLevelFun public
!#
!# Revision 2.76  2019/10/30 20:08:40  pwagner
!# Dumping Misc Notes spans lines properly
!#
!# Revision 2.75  2019/07/09 23:05:59  pwagner
!# Use Table cells to DumpGlobalAttributes
!#
!# Revision 2.74  2019/01/29 21:47:37  pwagner
!# Uses MLS_IsGlatt; light housekeeping
!#
!# Revision 2.73  2018/05/31 22:43:56  pwagner
!# Changed name of attribute to identifier_product_doi to please DAAC
!#
!# Revision 2.72  2018/05/31 18:06:49  pwagner
!# We can now copy over globalattributes
!#
!# Revision 2.71  2018/05/22 23:04:41  pwagner
!# Correct initial value of FailedChunks; may dump user-supplied GlobalAttribute_T
!#
!# Revision 2.70  2018/04/19 00:47:54  vsnyder
!# Remove USE statements for unused names
!#
!# Revision 2.69  2017/07/25 22:29:35  pwagner
!# Made 4 functions public that return granules day, month, dayofyear, year
!#
!# Revision 2.68  2016/07/26 17:44:31  pwagner
!# Deebug may be turned on by switch; print clearer message
!#
!# Revision 2.67  2015/08/12 20:36:38  pwagner
!# Added some missing global attrs; added h5_readGlobalAttr
!#
!# Revision 2.66  2015/02/27 23:57:21  pwagner
!# Take pains to ensure HostName is not blank
!#
!# Revision 2.65  2015/01/21 19:28:17  pwagner
!# Prevent unassociated anText error
!#
!# Revision 2.64  2014/10/02 17:20:46  pwagner
!# Avoid letting master write erroneous MiscNotes
!#
!# Revision 2.63  2014/09/05 00:15:47  vsnyder
!# More complete and accurate allocate/deallocate size tracking.  Remove
!# USE for unreferenced names.
!#
!# Revision 2.62  2014/04/14 16:59:02  pwagner
!# he5_writeglobalattr may skip writing if already there
!#
!# Revision 2.61  2014/04/02 23:04:06  pwagner
!# Removed redundant open_ and close_MLSFile
!#
!# Revision 2.60  2014/03/27 23:59:16  pwagner
!# he5_writeglobalattr now generic; DOI optional arg controls wwhether to write DOI, prodLoc
!#
!# Revision 2.59  2014/03/26 17:43:38  pwagner
!# Added ProductionLocation, identifier_product_DOI to attributes
!#
!# Revision 2.58  2014/01/09 00:24:29  pwagner
!# Some procedures formerly in output_m now got from highOutput
!#
!# Revision 2.57  2012/09/18 18:49:27  pwagner
!# Reduced severity when unable to convert utc; capitalize names in use statements
!#
!# Revision 2.56  2011/07/12 22:35:44  honghanh
!# Change l_grid to l_hdfeos
!#
!# Revision 2.55  2010/11/10 02:01:57  pwagner
!# Extra output if DEBUG
!#
!# Revision 2.54  2010/02/04 23:08:00  vsnyder
!# Remove use or declaration for unused names
!#
!# Revision 2.53  2010/01/15 01:12:37  pwagner
!# Added routines to read MLSFile_T components
!#
!# Revision 2.52  2010/01/11 18:36:07  pwagner
!# Changed attribute names to 'PCFID..'
!#
!# Revision 2.51  2010/01/08 00:10:46  pwagner
!# Added ability to write MLSFile_T fields as attributes
!#
!# Revision 2.50  2009/10/05 23:37:59  pwagner
!# Moved use mlshdfeos statements from module scope to speedup Lahey; this is the last time we do that
!#
!# Revision 2.49  2009/09/29 23:36:28  pwagner
!# Changes needed by 64-bit build
!#
!# Revision 2.48  2009/08/26 16:33:58  pwagner
!# Workaround for buggy hdfeos function he5_EHinqglatts; added dumpGlobalAttributes
!#
!# Revision 2.47  2009/06/23 18:25:42  pwagner
!# Prevent Intel from optimizing ident string away
!#
!# Revision 2.46  2008/12/02 23:10:30  pwagner
!# mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!#
!# Revision 2.45  2008/04/25 22:52:47  pwagner
!# Remove unused 'use ..'
!#
!# Revision 2.44  2008/03/07 01:37:36  pwagner
!# Relies on MLSHDFEOS to subdivide long textfiles into attributes
!#
!# Revision 2.43  2008/02/22 21:32:43  pwagner
!# Can now write leapsecfile, utcpole files as attrs, DS
!#
!# Revision 2.42  2007/10/04 20:43:36  vsnyder
!# Remove unused symbols
!#
!# Revision 2.41  2007/06/21 00:49:52  vsnyder
!# Remove tabs, which are not part of the Fortran standard
!#
!# Revision 2.40  2007/01/12 00:27:10  pwagner
!# New HostName global attribute written to product files
!#
!# Revision 2.39  2005/09/22 23:34:31  pwagner
!# date conversion procedures and functions all moved into dates module
!#
!# Revision 2.38  2005/08/15 20:38:39  pwagner
!# FirstMAf, LastMAF global attributes written for level 2 files only
!#
!# Revision 2.37  2005/07/12 17:14:22  pwagner
!# Global attribute MiscNotes added; InputVersion dropped
!#
!# Revision 2.36  2005/06/22 17:25:50  pwagner
!# Reworded Copyright statement, moved rcs id
!#
!# Revision 2.35  2005/06/14 20:35:24  pwagner
!# Changed interface to mls_io_gen functions
!#
!# Revision 2.34  2005/02/03 19:06:01  pwagner
!# utc_to_date used to find 0 crossing
!#
!# Revision 2.33  2004/08/16 17:05:38  pwagner
!# First,LastMAFCtr now global attributes
!#
!# Revision 2.32  2004/08/04 23:19:01  pwagner
!# Much moved from MLSStrings to MLSStringLists
!#
!# Revision 2.31  2004/03/24 23:53:02  pwagner
!# Switched from HE5T_NATIVE_SCHAR to MLS_CHARTYPE
!#
!# Revision 2.30  2004/02/26 22:01:04  pwagner
!# Acts more gracefully if l2gp file lacks global attributes
!#
!# Revision 2.29  2004/02/13 00:17:12  pwagner
!# New stuff for reading swath attributes
!#
!# Revision 2.28  2003/10/30 00:03:02  pwagner
!# Prepends 'L' to GlobalAttributes%ProcessLevel if necessary
!#
!# Revision 2.27  2003/10/28 00:39:00  pwagner
!# Fixed bug where character-vlaued attributes were only 1 char long
!#
!# Revision 2.26  2003/09/15 17:13:57  cvuu
!# Optional writing Orbit info to global attribute for h5_writeglobal
!#
!# Revision 2.25  2003/09/12 16:39:05  cvuu
!# Add attributes OrbitNumber and OrbitPeriod in the global attributes
!#
!# Revision 2.24  2003/08/15 20:41:50  pwagner
!# Wont try to write another /PCF dataset if already there
!#
!# Revision 2.23  2003/08/11 17:40:42  cvuu
!# Write anotating with multiple attributes PCF to handle big size of PCF file for hdfeos5
!#
!# Revision 2.22  2003/07/07 23:46:54  pwagner
!# Changed in interfaces to make filetype a lit_name
!#
!# Revision 2.21  2003/06/11 19:33:33  pwagner
!# No longer tries to write hdf5 pcf as both ds and attr
!#
!# Revision 2.20  2003/05/30 23:47:00  pwagner
!# Now standardized to write PCF, Inputpointer for all levels
!#
!# Revision 2.19  2003/04/11 23:33:13  pwagner
!# Gets he5_EHwrglatt from new MLSHDFEOS module
!#
!# Revision 2.18  2003/03/20 01:27:00  jdone
!# variable length string FileType used in WritePCF2Hdr
!#
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
!# Added WriteInputpointer
!#
!# Revision 2.5  2001/04/06 16:54:57  pwagner
!# Reverting to version 2.2
!#
!# Revision 2.2  2001/03/09 21:32:45  nakamura
!# Added intent(in) for pcf number arg.
!#
!# Revision 2.1  2001/03/09 21:10:32  nakamura
!# Routines for writing the PCF to an HDF file as an annotation.
!#
!#
