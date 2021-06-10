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
module L2GPData                 ! Creation, manipulation and I/O for L2GP Data
!=============================================================================
  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use HyperSlabs, only: ExtractArray, GatherBloc
  use BitStuff, only: DumpBitNames
  use Constants, only: Deg2Rad
  use Diff_1, only: Diff, Diff_Fun
  use Dump_Options, only: CrashAtBeginning, StatsOnOneLine, NameOnEachLine, &
    & NameHasbeenPrinted
  use Dump_0, only: Dump
  use Dump_1, only: Dump
  use Geometry, only: EarthRadA
  use HDF, only: Dfacc_Rdonly, Dfacc_Read, Dfacc_Create, Dfacc_Rdwr, &
    & Dfnt_Float32, Dfnt_Int32, Dfnt_Float64
  use HighOutput, only: BeVerbose, OutputNamedValue, StyledOutput
  use Intrinsic ! "units" Type Literals, Beginning With L
  use MLSCommon, only: DefaultUndefinedValue, Interval_T, &
    & MLSFile_T, L2MetaData_T, MLS_HyperStart, UndefinedIntegerValue
  use MLSFiles, only: FileNotFound, &
    & HDFVersion_4, HDFversion_5, WildcardHDFversion, &
    & Dump, InitializeMLSFile, MLS_CloseFile, MLS_Exists, MLS_OpenFile, &
    & MLS_HDF_Version, MLS_Inqswath
  use MLSKinds, only: R4, R8, Rt
  use MLSFillValues, only: IsfillValue, ReplacefillValues
  use MLSHDFEOS, only: HSize
  use MLSMessageModule, only: MLSMSG_Error, MLSMSG_Warning, MLSMessage
  use MLSNumerics, only: FindInRange
  use MLSFinds, only: FindFirst, FindLast, FindUnique
  use MLSSets, only: FindIntersection, Intersection
  use MLSStats1, only: MLSMin
  use MLSStrings, only: Capitalize, Lowercase
  use MLSStringLists, only: ExtractSubstring, &
    & GetHashElement, GetStringElement, GetUniqueList, &
    & List2Array, NumStringElements, RemoveListFromList, ReplaceSubstring, &
    & StringElementNum, SwitchDetail
  use Output_M, only: Blanks, Output, ResumeOutput, SuspendOutput
  use String_Table, only: Display_String
  use Time_M, only: SayTime, ConfigureSayTime, Time_Now
  use Trace_M, only: Trace_Begin, Trace_End

  implicit none

  private
  public :: L2GPData_T
  public :: AddL2GPToDatabase, AppendL2GPData, &
    & CompactL2GPRecord, ContractL2GPRecord, ConvertL2GPToQuantity, &
    & CpHE5GlobalAttrs, CpL2GPData, CpL2GPDataToAttribute, &
    & DestroyL2GPContents, DestroyL2GPDatabase, &
    & Diff, DiffRange, Dump, DumpRange, &
    & ExpandL2GPDataInPlace, ExtractL2GPRecord, FilterL2GP, &
    & IsL2GPSetup, OutputL2GP_Attributes, &
    & ReadL2GPData, RepairL2GP, SetupNewL2GPRecord, WriteL2GPData

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  interface ContractL2GPRecord
    module procedure ContractL2GPRecord_opt, ContractL2GPRecord_names
  end interface

  interface Diff
    module procedure DiffL2GPData
    module procedure DiffL2GPData_Chunks
    module procedure DiffL2GPData_Levels
    module procedure DiffL2GPFiles_MLSFile
    module procedure DiffL2GPFiles_Name
  end interface

  interface Diffrange
    module procedure DiffL2GPData_RANGES
  end interface

  interface Diffstats
    module procedure DiffStatsInt
  end interface

  interface Dump
    module procedure Dump_L2GP
    module procedure Dump_L2GP_Chunks
    module procedure Dump_L2GP_DataBase
    module procedure DumpL2GP_attributes_hdf5
  end interface

  interface DumpRange
    module procedure DumpL2GPData_Ranges
  end interface

  interface ExtractL2GPRecord
    module procedure ExtractL2GPRecord_range, ExtractL2GPRecord_which
  end interface

  interface ReadL2GPData
    module procedure ReadL2GPData_fileID
    module procedure ReadL2GPData_fileName
    module procedure ReadL2GPData_MLSFile
  end interface

  interface AppendL2GPData
    module procedure AppendL2GPData_fileName
    module procedure AppendL2GPData_MLSFile
  end interface

  interface CpL2GPData
    module procedure CpL2GPData_fileID
    module procedure CpL2GPData_fileName
    module procedure CpL2GPData_MLSFile
  end interface
  
  interface OutputL2GP_Attributes
    module procedure OutputL2GP_Attributes_MF
  end interface

  interface RepairL2GP
    module procedure RepairL2GP_L2GP
    module procedure RepairL2GP_HGrid
  end interface

  interface WriteL2GPData
    module procedure WriteL2GPData_fileID
    module procedure WriteL2GPData_MLSFile
  end interface

  ! This module defines datatypes and gives basic routines for storing and
  ! manipulating L2GP data.
  ! It is prepared to handle io for both file versions: hdfeos2 and hdfeos5
  ! For NetCDF4 interfaces, see NCL2GP.f90

!     c o n t e n t s
!     - - - - - - - -

! L2GPData_T              The l2gp data type; holds all data for one swath

! AddL2GPToDatabase       Adds an l2gp data type to a database
! AppendL2GPData          Appends L2GP onto end of existing swath file
! CompactL2GPRecord       Reduce a L2GP from an existing L2GP
!                           to just the non-Fill instances
! ContractL2GPRecord      Subsample an existing L2GP using geolocation ranges
! ConvertL2GPToQuantity   Convert an L2GP data type into a Vector Quantity
! CpHE5GlobalAttrs        Copies global attributes from one l2gp file to another
! CpL2GPData              Copies swaths from one l2gp file to another
! CpL2GPDataToAttribute   Copies swathvalues from one l2gp file to another's
!                           file-level attributes
! DestroyL2GPContents     Deallocates all the arrays allocated for an L2GP
! DestroyL2GPDatabase     Destroys an L2GP database
! Diff                    Shows differences between two swaths
! Dump                    Reveals info about an L2GP or a database of L2GP
!                            to any desired level of detail
! ExpandL2GPDataInPlace   Adds more profiles to an existing L2GP
! ExtractL2GPRecord       Subsample an existing L2GP using index ranges
! FilterL2GP              Set Status to Crash wherever Geolocations contain
!                            FillValues
! IsL2GPSetUp             Returns TRUE if all L2GP arrays allocated
! OutputL2GP_Attributes   Write all attributes relevant to an L2GP Data type
! ReadL2GPData            Reads L2GP from existing swath file
! RepairL2GP              Replaces fillValues in one L2GP with either
!                         (1) corresponding values from another L2GP; or
!                         (2) corresponding values from an HGrid
! SetupNewL2GPRecord      Allocates arrays for a new L2GP
! WriteL2GPData           Writes an L2GP into a swath file

! The meaning of options has replaced the older logical arguments
! if the options is present and contains the following characters:
!   character         meaning
!      ---            -------
!       c              repair geolcations of bad chunks
!       d              treat differing geolcations as bad chunks
!       f              filter profiles with FillValues by marking as crashed
!       g              skip diff or dump of geolcations
!       i              ignore bad chunks
!       m              silent (mute)
!       r              rms (or repair where appropriate)
!       s              stats
!       v              verbose
  ! Assume L2GP files w/o explicit hdfVersion field are this
  ! 4 corresponds to hdf4, 5 to hdf5 in L2GP, L2AUX, etc. 
  integer, parameter :: L2GPDefault_HDFVersion = HDFVersion_5
  integer, parameter :: DangerWillRobinson = 513 ! status for a crashed profile

  ! r4 corresponds to sing. prec. :: same as stored in files
  integer, public, parameter :: rgp = r4

  integer, public, parameter :: MAXCHUNKTIMES = 120
  integer, public, parameter :: MAXNUMSWATHPERFILE = 300
  ! How long may the list of swath names grow (~80 x max num. of swaths/file)
  integer, public, parameter :: MAXSWATHNAMESBUFSIZE = 80*MAXNUMSWATHPERFILE

  ! TRUE means we can avoid using unlimited dimension and its time penalty
  logical, public            :: AVOIDUNLIMITEDDIMS         = .true.
  logical, public            :: MUSTGUARDAGIANSTHDFEOSBUG  = .false.
  logical, public            :: WRITEMASTERSFILEATTRIBUTES = .true.

  integer, parameter, public :: CHARATTRLEN = 255   ! was GA_VALUE_LENGTH
  integer, parameter, public :: L2GPNameLen = 80
  integer, parameter, public :: NumDataFields = 5
  integer, parameter, public :: NumGeolocFields = 10
  integer, parameter, public :: MAXFNFIELDS = NumGeolocFields + NumDataFields + 4
  integer, parameter, public :: MAXNLEVELS = 100
  integer, parameter, public :: MAXNUMTIMES = 10000

   ! The following are the current data fields
   character (len=*), parameter, public :: DATA_FIELD1 = 'L2gpValue'
   character (len=*), parameter, public :: DATA_FIELD2 = 'L2gpPrecision'
   ! character (len=*), parameter, public :: DATA_FIELD3 = 'Status'
   ! character (len=*), parameter, public :: DATA_FIELD4 = 'Quality'
   ! character (len=*), parameter, public :: DATA_FIELD5 = 'Convergence'
   character (len=*), parameter, public :: DATA_FIELDS = &
     & 'L2gpValue,L2gpPrecision,Quality,Status,Convergence'

   ! The following are the current geolocation fields
   ! character (len=*), parameter, public :: GEO_FIELD1 = 'Latitude'
   ! character (len=*), parameter, public :: GEO_FIELD2 = 'Longitude'
   ! character (len=*), parameter, public :: GEO_FIELD3 = 'Time'
   ! character (len=*), parameter, public :: GEO_FIELD4 = 'LocalSolarTime'
   ! character (len=*), parameter, public :: GEO_FIELD5 = 'SolarZenithAngle'
   ! character (len=*), parameter, public :: GEO_FIELD6 = 'LineOfSightAngle'
   ! character (len=*), parameter, public :: GEO_FIELD7 = 'OrbitGeodeticAngle'
   ! character (len=*), parameter, public :: GEO_FIELD8 = 'ChunkNumber'
   ! character (len=*), parameter, public :: GEO_FIELD9 = 'Pressure'
   ! character (len=*), parameter, public :: GEO_FIELD10= 'Frequency'
   character (len=*), parameter, public :: GEO_FIELDS = &
     & 'Latitude,Longitude,LocalSolarTime,SolarZenithAngle,LineOfSightAngle' // &
     & ',OrbitGeodeticAngle,Pressure,Time,Frequency,ChunkNumber'
   ! The following are HIRDLS geolocation fields
   ! character (len=*), parameter, public :: HGEO_FIELDS = &
   !  & 'Latitude,Longitude,LocalSolarTime,SolarZenithAngle,LineOfSightAngle' // &
   !  & ',OrbitGeodeticAngle,Time'

   ! The following are the dimension names according to the mls spec
   character (len=*), parameter, public :: DIM_NAME1 = 'nTimes'
   ! character (len=*), parameter, public :: DIM_NAME2 = 'nLevels'
   ! character (len=*), parameter, public :: DIM_NAME3 = 'nFreqs'
   ! An alternate name for DIM_NAME3 used by HIRDLS
   ! character (len=*), parameter, public :: HRD_DIM_NAME3 = 'nChans'
   character (len=*), parameter, public :: DIM_NAME12 = 'nLevels,nTimes'
   character (len=*), parameter, public :: DIM_NAME123 = 'nFreqs,nLevels,nTimes'
   ! These are for the new max_dimlist parameter added to  SWdefgfld.
   ! this one is for non-extendible dimensions
   character (len=*), parameter, public :: MAX_DIML = ' '
   character (len=*), parameter, public :: UNLIM = 'Unlim'
   ! This is for cases where the time dimension is extendible
   character (len=*), parameter, public :: MAX_DIML1 = UNLIM
   character (len=*), parameter, public :: MAX_DIML12 = 'nLevels,Unlim'
   character (len=*), parameter, public :: MAX_DIML123 = 'nFreqs,nLevels,Unlim'

!   INTEGER,PARAMETER::CHUNKFREQS=13,CHUNKLEVELS=17,CHUNKTIMES=9,CHUNK4=1
   
   ! character (len=*), parameter, public :: DEFAULTMAXDIM = UNLIM
   ! integer, parameter, public :: DEFAULT_CHUNKRANK = 1
   ! integer, dimension(7), parameter, public :: DEFAULT_CHUNKDIMS = (/ 1,1,1,1,1,1,1 /)
   integer, parameter, public :: HDFE_NOMERGE = 0       ! don't merge
  
  ! So far, the nameIndex component of the main data type is never set
  logical, parameter, public :: NAMEINDEXEVERSET = .false.  
  
  ! Do you want to write file and swath attributes when you append values
  logical, parameter, public :: APPENDSWRITEATTRIBUTES = .true.  
  
  ! In what range of orbit angles do we set the mode to "descending"
  type(interval_T), public, parameter :: DescendingRange = interval_T( 90._rt, 270._rt )

  ! Do you want to pre-fill arrays with MissingValue by default?
  logical, parameter :: ALWAYSFILLWITHMISSINGVALUE = .true.  
  
  ! Deprecated, demoted, doomed
  logical, parameter, public :: AscDescModeIsField = .false.  

  ! The following pair, keys and hash, tell what value the units attribute
  ! will be given when writing column abundances in hdfeos5-formatted files
  ! You should override these if you need to do so
  character(len=255), public, save :: COL_SPECIES_KEYS = &
    & 'h2o,hno3,o3,hcl,clo,co,n2o,oh,rhi,so2,' &
    & // &
    & 'ho2,bro,hocl,hcn,iwc,ch3cn'
  character(len=255), public, save :: COL_SPECIES_HASH = &
    & 'molcm2,molcm2,DU,molcm2,molcm2,molcm2,molcm2,molcm2,molcm2,molcm2,' &
    & // &
    & 'molcm2,molcm2,molcm2,molcm2,molcm2,molcm2'
  character(len=16), dimension(10), parameter :: StatusBitNames = (/ &
    & 'do not use      ', 'beware          ', 'inform          ', &
    & 'post-processed  ', 'high clouds     ', 'low clouds      ', &
    & 'no meteorology  ', 'abandoned chunk ', 'too few radiance', &
    & 'mlsl2 crashed   ' /)
!      123456789012345678901234567890123456789012345678901234567890

  ! This datatype is the main one, it simply defines one l2gp swath
  ! It is used for 
  ! (1) normal swaths, which have geolocations along nTimes
  !     vertical coordinates called "surfaces" dimensioned by pressures
  !     and at nTimes horizontal "instances" dimensioned by 
  !     latitudes and longitudes
  ! (2) column abundances which are  integrated amounts and therefore 
  !     the vertical coordinate is suppressed; instead we set nLevels=1
  !     and store the tropopause in pressures(1)
  type L2GPData_T

     ! First some variables we use for book keeping

     character (len=L2GPNameLen) :: name ! Typically the swath name.
     integer :: nameIndex       ! Used by the parser to keep track of the data
     ! integer :: QUANTITYTYPE = 0   ! E.g., l_temperature; (Where is this used?)

     ! Now the dimensions of the data

     integer :: nTimes          ! Number of profiles in this slave
     integer :: nTimesTotal     ! Total number of profiles
     integer :: nLevels         ! Total number of surfaces (==1 for col. abund)
     integer :: nFreqs          ! Number of frequencies in breakdown

     ! Now we store the geolocation fields, first the vertical one:
     ! (The following has the tropopause if the swath is a column abundance)
     real (rgp), pointer, dimension(:) :: pressures=>NULL() ! Vertical coords (nLevels)

     ! Now the horizontal geolocation information. Dimensioned (nTimes)
     real (rgp), pointer, dimension(:) :: latitude => NULL()
     real (rgp), pointer, dimension(:) :: longitude => NULL()
     real (rgp), pointer, dimension(:) :: solarTime => NULL()
     real (rgp), pointer, dimension(:) :: solarZenith => NULL()
     real (rgp), pointer, dimension(:) :: losAngle => NULL()
     real (rgp), pointer, dimension(:) :: geodAngle => NULL()
     real (r8 ), pointer, dimension(:) :: time => NULL()   ! dble prec.

     integer, pointer, dimension(:) :: chunkNumber=>NULL()

     ! Now we store the `frequency' geolocation field

     real (rgp), pointer, dimension(:) :: frequency=>NULL()
     !        dimensioned (nFreqs)

     ! Finally we store the data fields

     real (rgp), pointer, dimension(:,:,:) :: l2gpValue=>NULL()
     real (rgp), pointer, dimension(:,:,:) :: l2gpPrecision=>NULL()
     ! dimensioned (nFreqs, nLevels, nTimes)

     ! Now we've changed our minds: status will be a 4-byte integer
     integer, pointer, dimension(:)    :: status=>NULL()
     !                (status is a reserved word in F90)
     real (rgp), pointer, dimension(:) :: quality=>NULL()
     real (rgp), pointer, dimension(:) :: convergence=>NULL()
     
     ! These fields applies only to products of a neural network model
     integer, pointer, dimension(:)    :: BinNumber=>NULL()
     integer, pointer, dimension(:)    :: MAF=>NULL()
     ! We have deprecated this field
     integer, pointer, dimension(:)    :: AscDescMode=>NULL()
     ! All the above dimensioned (nTimes)

     ! These are the fill/missing values for l2gpValue
     real (rgp)                        :: MissingL2GP  = DefaultUndefinedValue
     ! These are the fill/missing values for all other arrays except status
     real (rgp)                        :: MissingValue = DefaultUndefinedValue
     integer                           :: MissingStatus = DangerWillRobinson
    ! Vertical coordinate
    character(len=8) :: verticalCoordinate ! E.g. 'Pressure', or 'Theta'
    ! integer :: verticalCoordinate ! The vertical coordinate used.  These
                                  ! are l_lits of the type t_VGridCoord
                                  ! defined in Init_Tables_Module.
  ! The following final method caused the sids tests
  ! in the nightly gold brick to bomb. A symptom of a deeper problem?
  ! contains
    ! final :: DestroyL2GPContents
  end type L2GPData_T

  ! Print debugging stuff?
  logical, parameter :: DEEBUG = .false.  
  logical, parameter ::SWATHLEVELMISSINGVALUE = .false. ! Make it swath attr?
  real :: t2

contains ! =====     Public Procedures     =============================

  !-------------------------------------------------------------

  ! ---------------------- AppendL2GPData_fileName  ---------------------------

  subroutine AppendL2GPData_fileName( l2gp, fileName, &
    & swathname, offset, lastProfile, TotNumProfs, hdfVersion, &
    & createSwath, maxchunksize )
    !------------------------------------------------------------------------

    ! Given a file name,
    ! This routine does an append operation (see AppendL2GPData_fileID)
    ! If the file doesn't exist yet, hopefully, it'll create it

    ! Arguments

    character (len=*), intent(in) :: fileName ! Name of file to append to
    type( l2GPData_T ), intent(inout) :: l2gp ! What to append
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    integer,intent(in),optional   :: offset
    integer, intent(in), optional ::lastProfile
    integer,intent(in),optional   :: TotNumProfs
    integer, optional, intent(in) :: hdfVersion
    logical, intent(in), optional :: createSwath
    integer, optional, intent(in) :: maxchunksize
    ! logical, optional, intent(in) :: clean   ! Not implemented yet

    ! Local
    integer                       :: file_access    
    logical                       :: file_exists    
    type(MLSFile_T)               :: MLSFile        
    logical                       :: myClean        
    integer                       :: status         
    integer                       :: the_hdfVersion 
    
    ! Executable code
    myClean = .false.
    ! if ( present(clean) ) myClean = clean
    file_exists = ( mls_exists(trim(FileName)) == 0 )
    the_hdfVersion = L2GPDEFAULT_HDFVERSION
    if ( present(hdfVersion) ) the_hdfVersion = hdfVersion
    if ( the_hdfVersion == WILDCARDHDFVERSION ) then


      ! Does the file exist, yet?
      if ( .not. file_exists ) then
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File not found; make sure the name and path are correct' &
          & // trim(fileName) )
      endif
      the_hdfVersion = mls_hdf_version(FileName, hdfVersion)
      if ( the_hdfVersion == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File not found; make sure the name and path are correct' &
          & // trim(fileName) )
    endif
    if ( file_exists ) then
      file_access = DFACC_RDWR
    else
      file_access = DFACC_CREATE
    endif
    status = InitializeMLSFile ( MLSFile, type=l_swath, access=file_access, &
     & name=trim(fileName), content = 'l2gp', HDFVersion=the_hdfVersion )
    ! call MLS_OpenFile( MLSFile )
    ! L2FileHandle = MLSFile%FileID%f_id
    call AppendL2GPData_MLSfile( l2gp, MLSFile, swathname, &
      & offset, lastProfile=lastProfile, totNumProfs=totNumProfs, &
      & createSwath=createSwath, &
      & maxchunksize=maxchunksize )
    ! call MLS_CloseFile( MLSFile )
  end subroutine AppendL2GPData_fileName
  ! ---------------------- AppendL2GPData_MLSFile  ---------------------------

  subroutine AppendL2GPData_MLSfile( l2gp, l2gpFile, &
    & swathName, offset, lastProfile, TotNumProfs, &
    & createSwath, maxchunksize )
    ! sticks l2gp into the swath swathName in the file pointed at by
    ! l2FileHandle,starting at the profile number "offset" (First profile
    ! in the file has offset==0). If this runs off the end of the swath, 
    ! it is lengthened automagically. 
    ! This call has been altered recently, so that it can be used to create
    ! a swath as well as adding to one. 

    use MLSHDFEOS, only: MLS_Swath_In_File

    ! Arguments

    type(MLSFile_T)                :: L2GPFile

    ! This is a L2GPData_T structure containing all the data to be written
    type (L2GPData_T), intent(INOUT) :: l2gp
    ! This is the name the swath is given in the file. By default it is
    ! the name contained in l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    ! This (offset) is the point in the swath at which the data is written. 
    ! First profile in the file has offset==0. If the swath in the file is 
    ! shorter than offset + ( num of profiles in l2gp) then it grows by magic
    integer, intent(in), optional::offset
    ! TotNumProfs is a new argument. It seems only to be used if we are 
    ! creating a swath, rather than adding to one. In that case I guess
    ! it is the total number of profiles in the swath created. I also 
    ! guess that this is done so that we can avoid growing and re-growing 
    ! the swath.
    integer, intent(in), optional ::TotNumProfs
    integer, intent(in), optional ::lastProfile
    logical, intent(in), optional :: createSwath
    integer, optional, intent(in) :: maxchunksize
    ! Local
    integer :: actual_ntimes
    logical :: alreadyOpen
    type (L2GPData_T) :: largerl2gp
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: myLastProfile
    character (len=L2GPNameLen) :: myswathName
    logical :: notUnlimited
    integer :: numProfs
    integer :: status
    logical :: swath_exists
    logical :: timing
    real :: tFile ! How long have we been fooling with this file?
    type (L2GPData_T) :: totall2gp
    ! logical, parameter :: DEEBUG = .false.

    ! Executable code
    call trace_begin ( me, 'AppendL2GPData_MLSFile', cond=.false. )
    ! call Dump ( l2gp%chunkNumber, 'Appending l2gp%chunkNumber' )
    call time_now ( tFile )
    call configureSayTime ( tFile )
    status = 0
    timing = DEEBUG .or. BeVerbose( 'l2gp', 2 )

    if (present(lastProfile)) then
      myLastProfile = lastProfile 
    elseif (present(TotNumProfs)) then
      myLastProfile = TotNumProfs 
    else
      myLastProfile = L2GP%nTimesTotal
    endif
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      if ( DEEBUG ) print *, 'Needed to open file ', trim(L2GPFile%name)
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    if ( timing ) then
      call sayTime( 'Opening file ' // trim(L2GPFile%name), tFile, t2 )
      tFile = t2
    endif
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile)
    myswathName = l2gp%name
    if ( present(swathName) ) myswathName = swathName
    
    swath_exists = mls_swath_in_file( L2GPFile%name, myswathName, &
      & L2GPFile%HdfVersion )

    ! if ( timing ) call outputNamedValue( 'tFile', tFile )
    if ( timing ) then
      call sayTime( 'Checking that swath exists', tFile, t2 )
      tFile = t2
    endif
    if ( swath_exists .and. &
      & ( MUSTGUARDAGIANSTHDFEOSBUG .or. .not. present(maxChunkSize) ) &
      & ) then
      if(DEEBUG) print *, 'Swath already exists--guarding against hdfeos bug'
      ! We must guard against a bug in HDF-EOS that prevents appending more
      ! profiles to a swath than already exist
      call ReadL2GPData( L2GPFile, myswathName, totall2gp, numProfs, &
        & ReadData=.false. )
      if ( l2gp%nTimes > numProfs ) then
        if(DEEBUG) call outputNamedValue( 'current', totall2gp%nTimes )
        if(DEEBUG) call outputNamedValue( 'needed', l2gp%nTimes )
        call DestroyL2GPContents ( totall2gp )
        call ReadL2GPData( L2GPFile, myswathName, totall2gp, numProfs )
        call ExpandL2GPDataInPlace( totall2gp, l2gp%nTimes )
        if(DEEBUG) call outputNamedValue( 'expanded', totall2gp%nTimes )
        call ExpandL2GPDataInFile( L2GPFile, myswathName, totall2gp )
        if(DEEBUG) call outputNamedValue( 'AppendL2GPData expanded', totall2gp%nTimes )
      elseif(DEEBUG) then
         print *, 'Swath already has enough profiles'
      endif
      call DestroyL2GPContents ( totall2gp )
      if ( timing ) then
        call sayTime( 'Guarding against HDFEOS error', tFile, t2 )
        tFile = t2
      endif
    elseif ( swath_exists ) then
      if(DEEBUG) print *, 'OK, swath already exists, but HDFEOS bug no longer scares us'
    else
      ! Must create swath in file w/o disturbing other swaths
      if(DEEBUG) print *, 'Must create swath'
      if(DEEBUG) print *, 'Will have ', myLastProfile, ' profiles'
      if(DEEBUG) print *, 'instead of ', l2gp%nTimes, ' profiles'
      actual_ntimes = l2gp%nTimes
      ! if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
      select case (L2GPFile%hdfVersion)
      case (HDFVERSION_4)
        ! Currently force unlimited, remove the .false. .and. to allow limited
        if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
        call OutputL2GP_createFile_MF (l2gp, L2GPFile, &
          & myswathName, notUnlimited=.false. .and. present(totNumProfs))
        l2gp%nTimes = actual_ntimes
      case (HDFVERSION_5)
        ! By default allow limited; 
        ! may force unlimited by setting avoidUnlimitedDims to FALSE
        notUnlimited = ( avoidUnlimitedDims .and. present(totNumProfs) )
        ! if ( present(TotNumProfs) ) l2gp%nTimes = TotNumProfs
        if ( present(maxchunksize) ) then
          if ( maxchunksize > l2gp%nTimes ) then
            call SetupNewL2GPRecord ( largerl2gp, proto=l2gp, &
              & nTimes=maxchunksize )
            call OutputL2GP_createFile_MF ( largerl2gp, L2GPFile, &
              & myswathName, notUnlimited=notUnlimited )
            call DestroyL2GPContents ( largerl2gp )
          else
            call OutputL2GP_createFile_MF ( l2gp, L2GPFile, &
              & myswathName, notUnlimited=notUnlimited )
          endif
        else
          call OutputL2GP_createFile_MF ( l2gp, L2GPFile, &
            & myswathName, notUnlimited=notUnlimited )
        endif
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
         & 'Illegal hdf version in AppendL2GPData_fileID', MLSFile=L2GPFile)
      end select
      ! l2gp%nTimes = actual_ntimes
    endif

    if(DEEBUG) then
      if ( present(offset) ) then
        print*,"offset=",offset,"myLastProfile=",myLastProfile,&
        "size(l2gp%l2gpValue,3)=",size(l2gp%l2gpValue,3)
      else
        print*,"no offset; myLastProfile=",myLastProfile,&
        "size(l2gp%l2gpValue,3)=",size(l2gp%l2gpValue,3)
      endif
    endif
    if ( size(l2gp%l2gpValue,3) == 0 ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, &
        & "No profiles in this chunk", MLSFile=L2GPFile )

    else
      ! actual_ntimes = l2gp%nTimes
      ! l2gp%nTimes = max(myLastProfile - offset + 1, 1)
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'all chunk numbers are -999', advance='yes' )
      call OutputL2GP_writeGeo_MF (l2gp, l2GPFile, &
        & myswathName, offset)
      if ( timing ) then
        call sayTime( 'Writing geolocations', tFile, t2 )
        tFile = t2
      endif
      call OutputL2GP_writeData_MF (l2gp, l2GPFile, &
        & myswathName, offset)
      if ( timing ) then
        call sayTime( 'Writing data', tFile, t2 )
        tFile = t2
      endif
      select case ( L2GPFile%HDFVersion )
      case ( HDFVERSION_4 )
      case ( HDFVERSION_5 )
        if ( .not. swath_exists .and. APPENDSWRITEATTRIBUTES) then
          call OutputL2GP_Attributes_MF (l2gp, l2GPFile, swathName)
          call SetL2GP_aliases_MF (l2gp, l2GPFile, swathName)
          if ( timing ) then
            call sayTime( 'Writing attributes', tFile, t2 )
            tFile = t2
          endif
        end if
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & "Unrecognized hdfVersion passed to AppendL2GPData", MLSFile=L2GPFile )
      end select
      ! l2gp%nTimes = actual_ntimes
    end if

    if ( .not. alreadyOpen )  call mls_closeFile(L2GPFile, Status)
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'append'
    call trace_end ( 'AppendL2GPData_MLSFile', cond=.false. )
  end subroutine AppendL2GPData_MLSFile

  !------------------------------------------  CompactL2GPRecord  -----
  ! Trim all the flabby FillValues from an L2GP record producing a slim
  ! athletic record suitable for Olympic tryouts

  subroutine CompactL2GPRecord ( ol2gp, l2gp )
    ! Dummy arguments
    type (L2GPData_T), intent(in)   ::     ol2gp
    type (L2GPData_T), intent(out)  ::     l2gp
    ! Local variables
    integer, parameter              :: MaxFreqs  = 150
    integer, parameter              :: MaxLevels = 150
    integer, parameter              :: MaxTimes  = 7200
    integer, dimension(MaxFreqs )   :: whichFreqs  ! which freqs        
    integer, dimension(MaxLevels)   :: whichLevels ! which levels       
    integer, dimension(MaxTimes )   :: whichTimes  ! which times        
    integer                         :: useFreqs  ! how many freqs
    integer                         :: useLevels ! how many levels
    integer                         :: useTimes  ! how many times
    integer                         :: k
    ! Executable
    useFreqs  = 0
    useLevels = 0
    useTimes  = 0
    ! Discover which values are non-Fill
    if (associated(ol2gp%frequency) ) then
      do k=1, size(ol2gp%frequency)
        if ( .not. isFillValue( ol2gp%frequency(k) ) ) then
          useFreqs = useFreqs + 1
          whichFreqs(useFreqs) = k
        endif
      enddo
    endif
    if (associated(ol2gp%pressures) ) then
      do k=1, size(ol2gp%pressures)
        if ( .not. isFillValue( ol2gp%pressures(k) ) ) then
          useLevels = useLevels + 1
          whichLevels(useLevels) = k
        endif
      enddo
    endif
    if (associated(ol2gp%time) ) then
      do k=1, size(ol2gp%time)
        if ( .not. isFillValue( ol2gp%time(k) ) ) then
          useTimes = useTimes + 1
          whichTimes(useTimes) = k
        endif
      enddo
    endif
    if ( DeeBug ) then
      call outputnamedValue ( 'useFreqs ', useFreqs  )
      call outputnamedValue ( 'useLevels', useLevels )
      call outputnamedValue ( 'useTimes ', useTimes  )
    endif
    call ExtractL2GPRecord ( ol2gp, l2gp, &
      & whichFreqs, whichLevels, whichTimes, &
      & useFreqs, useLevels, useTimes )
  end subroutine CompactL2GPRecord

  !------------------------------------------  ContractL2GPRecord  -----
  ! Gather a reduced copy of an original l2gp record
  ! picking out a set of subscripts for freqs, level2, times
  ! based on input range(s) for latitude, pressure, etc.

  ! A trick: if the lowBound > hiBound then the sense is reversed
  ! i.e. find subscripts outside corresponding range

  ! Note that this does not do the more useful task:
  ! Reposition an l2gp record onto a new set of geolocations
  ! via multi-dimensional or repeated interpolation

  subroutine ContractL2GPRecord_names ( ol2gp, l2gp, &
    & geoBoxNames, geoBoxLowBound, geoBoxHiBound )
    ! Dummy arguments
    type (L2GPData_T), intent(in)   ::     ol2gp
    type (L2GPData_T), intent(out)  ::     l2gp
    character(len=*), intent(in)  ::       geoBoxNames
    real(rgp), dimension(:), intent(in) :: geoBoxLowBound  ! range
    real(rgp), dimension(:), intent(in) :: geoBoxHiBound  ! range
    ! Internal variables
    logical, parameter :: countEmpty = .true.
    integer :: elem
    real(rgp), dimension(2) :: pressures  ! range
    real(rgp), dimension(2) :: latitudes  ! range
    real(rgp), dimension(2) :: longitudes  ! range
    real(r8), dimension(2) :: times  ! range
    ! Executable
    pressures(1)  = ol2gp%MissingValue
    latitudes(1)  = ol2gp%MissingValue
    longitudes(1) = ol2gp%MissingValue
    times(1)      = ol2gp%MissingValue
    pressures(2)  = ol2gp%MissingValue
    latitudes(2)  = ol2gp%MissingValue
    longitudes(2) = ol2gp%MissingValue
    times(2)      = ol2gp%MissingValue
    elem = stringElementNum( lowerCase(geoBoxNames), 'pressure', countEmpty )
    if ( elem > 0 ) then
      pressures(1)   = geoBoxLowBound(elem)
      pressures(2)   = geoBoxHiBound(elem)
    endif
    elem = stringElementNum( lowerCase(geoBoxNames), 'latitude', countEmpty )
    if ( elem > 0 ) then
      latitudes(1)   = geoBoxLowBound(elem)
      latitudes(2)   = geoBoxHiBound(elem)
    endif
    elem = stringElementNum( lowerCase(geoBoxNames), 'longitude', countEmpty )
    if ( elem > 0 ) then
      longitudes(1)   = geoBoxLowBound(elem)
      longitudes(2)   = geoBoxHiBound(elem)
    endif
    elem = stringElementNum( lowerCase(geoBoxNames), 'time', countEmpty )
    if ( elem > 0 ) then
      times(1)   = geoBoxLowBound(elem)
      times(2)   = geoBoxHiBound(elem)
    endif
    call ContractL2GPRecord( ol2gp, l2gp, pressures, latitudes, &
    & longitudes, times )
  end subroutine ContractL2GPRecord_names

  subroutine ContractL2GPRecord_opt ( ol2gp, l2gp, pressures, latitudes, &
    & longitudes, intimes, hoursInDay, chunks )
    ! Dummy arguments
    type (L2GPData_T), intent(in)   :: ol2gp
    type (L2GPData_T), intent(out)  :: l2gp
    real(rgp), dimension(2), intent(in), optional :: pressures  ! range
    real(rgp), dimension(2), intent(in), optional :: latitudes  ! range
    real(rgp), dimension(2), intent(in), optional :: longitudes ! range
    real(r8), dimension(2), intent(in), optional  :: intimes    ! range
    real(rgp), dimension(2), intent(in), optional :: hoursInDay ! range
    integer, dimension(2), intent(in), optional   :: chunks     ! range
    ! Local variables
    integer :: i, n
    integer, dimension(:), allocatable :: intrsctn
    integer, dimension(4000) :: tempTimes  ! subscripts
    real(r8), dimension(2)  :: times
    integer :: useFreqs  
    integer :: useLevels 
    integer :: useTimes  
    integer, dimension(300) :: whichFreqs ! subscripts
    integer, dimension(60) :: whichLevels ! subscripts
    integer, dimension(4000) :: whichTimes  ! subscripts

    ! logical, parameter :: DEEBug = .true.
    
    ! Executable
    do i=1, max( ol2gp%nFreqs, 1 )
      whichFreqs(i)=i
    enddo
    do i=1, max( ol2gp%nTimes, 1 )
      whichTimes(i)=i
    enddo
    do i=1, max( ol2gp%nLevels, 1 )
      whichLevels(i)=i
    enddo
    useFreqs = max( ol2gp%nFreqs, 1 )
    useTimes = ol2gp%nTimes
    useLevels = max(ol2gp%nLevels, 1)
    if ( present(pressures) ) then
      if ( any(IsFillValue(pressures, ol2gp%MissingValue)) ) then
        ! we do nothing
      elseif ( pressures(2) >= pressures(1) ) then
        call FindInRange( ol2gp%pressures, pressures, whichLevels, useLevels )
      else
        ! Reverse sense (i.e., find outside of range)
        call FindInRange( ol2gp%pressures, pressures, whichLevels, useLevels, options='-r' )
      endif
    endif
    if ( present(latitudes) ) then
      if ( any(IsFillValue(latitudes, ol2gp%MissingValue)) ) then
        ! we do nothing
      else
        if ( latitudes(2) >= latitudes(1) ) then
          call FindInRange( ol2gp%latitude, latitudes, tempTimes, n )
        else
          ! Reverse sense (i.e., find outside of range)
          call FindInRange( ol2gp%latitude, latitudes, tempTimes, n, options='-r' )
        endif
        intrsctn = Intersection( whichTimes(1:useTimes), tempTimes(1:n) )
        useTimes = size(intrsctn)
        whichTimes(1:useTimes) = intrsctn
        if ( DeeBug ) then
          call dump( latitudes, 'latitudes' )
          call outputNamedValue( 'n', n )
          call dump( tempTimes(1:n), 'tempTimes' )
          call dump( intrsctn(1:useTimes), 'intrsctn' )
        endif
        call deallocate_test( intrsctn, 'intersection with lats', ModuleName )
      endif
    endif
    if ( present(longitudes) ) then
      if ( any(IsFillValue(longitudes, ol2gp%MissingValue)) ) then
        ! we do nothing
      else
        if ( longitudes(2) >= longitudes(1) ) then
          call FindInRange( ol2gp%longitude, longitudes, tempTimes, n )
        else
          ! Reverse sense (i.e., find outside of range)
          call FindInRange( ol2gp%longitude, longitudes, tempTimes, n, options='-r' )
        endif
        intrsctn = Intersection( whichTimes(1:useTimes), tempTimes(1:n) )
        useTimes = size(intrsctn)
        whichTimes(1:useTimes) = intrsctn
        if ( DeeBug ) then
          call dump( longitudes, 'longitudes' )
          call outputNamedValue( 'n', n )
          call dump( tempTimes(1:n), 'tempTimes' )
          call dump( intrsctn(1:useTimes), 'intrsctn' )
        endif
        call deallocate_test( intrsctn, 'intersection with lons', ModuleName )
      endif
    endif
    ! Time may be represented as either TAI (s after midnight Jan 1 1993)
    if ( present(intimes) ) then
      times = intimes
    elseif ( present(hoursInDay) ) then
      times = hoursInDayToTime( hoursInDay, ol2gp )
    endif
    if ( present(intimes) .or. present(hoursInDay) ) then
      if ( any(IsFillValue(times, real(ol2gp%MissingValue, r8))) ) then
        ! we do nothing
      else
        if ( times(2) >= times(1) ) then
          call FindInRange( ol2gp%time, times, tempTimes, n )
        else
          ! Reverse sense (i.e., find outside of range)
          call FindInRange( ol2gp%time, times, tempTimes, n, options='-r' )
        endif
      endif
    elseif ( present(chunks) ) then
      call FindInRange( ol2gp%chunkNumber, chunks, tempTimes, n )
    endif
    if ( present(intimes) .or. present(hoursInDay) .or. present(chunks) ) then
      intrsctn = Intersection( whichTimes(1:useTimes), tempTimes(1:n) )
      useTimes = size(intrsctn)
      whichTimes(1:useTimes) = intrsctn
      if ( DeeBug ) then
        call dump( times, 'times' )
        call outputNamedValue( 'n', n )
        call dump( tempTimes(1:n), 'tempTimes' )
        call dump( intrsctn(1:useTimes), 'intrsctn' )
      endif
      call deallocate_test( intrsctn, 'intersection with times', ModuleName )
    endif
    call ExtractL2GPRecord ( ol2gp, l2gp, whichFreqs, whichLevels, whichTimes, &
    & useFreqs, useLevels, useTimes )
  end subroutine ContractL2GPRecord_opt

  !------------------------------------------  ConvertL2GPToQuantity  -----
  ! Convert an existing l2gp data type to an equivalent
  ! vector quantity

  subroutine ConvertL2GPToQuantity ( l2gp, Quantity )
    use QuantityTemplates, only: SetUpNewQuantityTemplate
    use VectorsModule, only: VectorValue_T, CreateVectorValue
    ! Dummy arguments
    type (L2GPData_T), intent(in)   ::     l2gp
    type (VectorValue_T), intent(out)  ::  Quantity
    ! Internal variables
    ! Executable
    call SetupNewQuantityTemplate ( Quantity%template, l2gp%nTimes, l2gp%nLevels, &
      & 1, .true., .true., .true. )
    call CreateVectorValue ( Quantity, 'L2GP'  )
    if ( associated(l2gp%l2gpValue) ) Quantity%values = l2gp%l2gpValue(1,:,:)
    ! Don't bother copying any of the gelocations if l2gp%latitudes not associated
    if ( .not. associated(l2gp%latitude) ) return
    Quantity%template%geodLat(1,:)        = l2gp%latitude
    Quantity%template%lon(1,:)            = l2gp%longitude
    Quantity%template%time(1,:)           = l2gp%time
    Quantity%template%losAngle(1,:)       = l2gp%losAngle
    Quantity%template%solarTime(1,:)      = l2gp%solarTime
    Quantity%template%solarZenith(1,:)    = l2gp%solarZenith
    ! If SetupNewQuantityTemplate edxplained clearly and accurately
    ! about the 2 indexes of the surfs array, this if bloc wouldn't be needed.
    ! However, all the chaff about regular, coherent, and stacked
    ! possibilities, with no effort at explaining them,
    ! leaves us no recourse but what might be called
    !   d e f e n s i v e   p r o g r a m m i n g
    if ( size(Quantity%template%surfs,1) == l2gp%nLevels ) then
      Quantity%template%surfs(:,1)          = l2gp%pressures
    elseif ( size(Quantity%template%surfs,2) == l2gp%nLevels ) then
      Quantity%template%surfs(1,:)          = l2gp%pressures
    ! else
      ! call MLSMessage ( MLSMSG_Error, ModuleName, &
      !   & 'ConvertL2GPToQuantity failed to cope with l2gp pressures' )
    endif
  end subroutine ConvertL2GPToQuantity

  ! ---------------------- cpHE5GlobalAttrs  ---------------------------

  subroutine cpHE5GlobalAttrs( File1Handle, File2Handle, status )
  ! Copy global attributes from file 1 to file 2
    use PCfhdr, only: GlobalAttributes_T, GlobalAttributes, &
      & DumpGlobalAttributes, HE5_Readglobalattr, HE5_Writeglobalattr
    ! Args
    integer, intent(in)                          :: File1Handle
    integer, intent(in)                          :: File2Handle
    integer, intent(out)                         :: status
    ! Internal variables
    type(GlobalAttributes_T)                     :: gAttributes        
    type(GlobalAttributes_T)                     :: gAttributesOriginal
    character(len=40)                            :: ProcessLevel       
    double precision                             :: TAI93At0zOfGranule 
    integer                                      :: DayofYear          
    ! Executable
    if ( DEEBUG ) then
      call output( 'Before reading global attributes', advance='yes' )
      call outputNamedValue( 'StartUTC', trim(GlobalAttributes%StartUTC) )
    endif
    call he5_readglobalattr ( File1Handle, gAttributes, &
      & ProcessLevel, DayofYear, TAI93At0zOfGranule, returnStatus=status )
    if ( status == 0 ) then
      if ( DEEBUG ) then
        call output ( '(Global Attributes read) ', advance='yes')
        call dumpGlobalAttributes
      endif
      if ( DEEBUG ) then
        call output( 'After reading global attributes', advance='yes' )
        call outputNamedValue( 'StartUTC', trim(GlobalAttributes%StartUTC) )
      endif
      ! Unfortunately, he5_writeglobalattr writes class-level data
      ! so must save original version so can copy new data into it 
      gAttributes%HostName = GlobalAttributes%HostName
      gAttributes%ProductionLoc = GlobalAttributes%ProductionLoc
      gAttributesOriginal = GlobalAttributes
      GlobalAttributes = gAttributes
      ! Why must we do this? Why do things not simply work out properly?
      GlobalAttributes%TAI93At0zOfGranule = TAI93At0zOfGranule
      if ( DEEBUG ) then
        call outputNamedValue('Misc Notes (written) ', trim(GlobalAttributes%MiscNotes) )
        call dumpGlobalAttributes
      endif
      call he5_writeglobalattr (File2Handle)
      ! Before leaving must restore original data back
      GlobalAttributes = gAttributesOriginal
    endif
  end subroutine cpHE5GlobalAttrs

  ! ---------------------- cpL2GPData_fileID  ---------------------------

  subroutine cpL2GPData_fileID( l2metaData, file1, file2, swathList, &
    & hdfVersion1, hdfVersion2, notUnlimited, rename, ReadData, &
    & HGrid, rFreqs, rLevels, rTimes, options )
    !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies swathList from 1 to 2
    ! and, depending on options, makes some repairs or applies a sanity filter
    ! which guarantees that any profiles with Fill values among the
    ! geolocations are marked with MissingStatus status

    use HGridsDatabase, only: HGrid_T
    ! Arguments

    type (L2Metadata_T) :: l2metaData
    integer, intent(in)           :: file1 ! handle of file 1
    integer, intent(in)           :: file2 ! handle of file 2
    character (len=*), intent(in) :: swathList ! copy only these; no wildcard
    integer, intent(in)           :: hdfVersion1
    integer, intent(in)           :: hdfVersion2
    logical, optional, intent(in) :: notUnlimited
    logical, optional, intent(in) :: ReadData
    character (len=*), optional, intent(in) :: rename
    type (HGrid_T), optional, intent(in)    :: HGrid
    integer, dimension(2), intent(in), optional :: rFreqs  ! subscript range
    integer, dimension(2), intent(in), optional :: rLevels ! subscript range
    integer, dimension(2), intent(in), optional :: rTimes  ! subscript range
    character (len=*), optional, intent(in) :: options ! E.g., '-v'

    ! Local variables
!     integer :: chunk
    logical, parameter            :: countEmpty = .true.
    logical :: filter
    integer :: i
    type (L2GPData_T) :: l2gp
    type (L2GPData_T) :: reducedl2gp
    integer, dimension(2) :: myFreqs  ! subscript range
    integer, dimension(2) :: myLevels ! subscript range
    integer, dimension(2) :: myTimes  ! subscript range
    character (len=8) :: myOptions
    integer :: noSwaths
    logical :: renameSwaths
    logical :: repair
    character (len=L2GPNameLen) :: swath
    character (len=L2GPNameLen) :: swath2
    logical :: verbose
    
    ! Executable code
    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = ( index(myOptions, 'v') > 0 )
    repair  = ( index(myOptions, 'r') > 0 )
    filter  = ( index(myOptions, 'f') > 0 )
    renameSwaths  = present(rename)
    if ( renameSwaths ) renameSwaths = ( rename /= ' ' )
    if ( repair .and. .not. present(HGrid) ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'cpL2GPFile must be given HGrid to repair L2GPData' )
    if ( .not. repair .and. present(HGrid) ) &
      & call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'cpL2GPFile was given HGrid but no order to repair L2GPData' )
    noSwaths = NumStringElements( trim(swathList), countEmpty )
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths to cp to file--unable to count swaths in ' // trim(swathList) )
    endif
    if ( verbose ) call dump(swathlist, 'swath names')
    if ( renameSwaths ) then
      noSwaths = min(noSwaths, NumStringElements(trim(rename), countEmpty))
      if ( verbose ) call dump(rename, 'swath names (copied)')
    endif
    if ( noSwaths < 1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'No swaths cp to file--unable to count swaths in ' // trim(rename) )
    endif
    ! Loop over swaths in file 1
    do i = 1, noSwaths
      call GetStringElement (trim(swathList), swath, i, countEmpty )
      if ( renameSwaths ) then
        call GetStringElement (trim(rename), swath2, i, countEmpty )
      else
        swath2 = swath
      endif
      ! Allocate and fill l2gp
      if ( DEEBUG ) print *, 'Reading swath from file: ', trim(swath)
      call ReadL2GPData ( file1, trim(swath), l2gp, &
           & hdfVersion=hdfVersion1, ReadData=ReadData )
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'all chunk numbers are -999', advance='yes' )
      if ( DEEBUG ) then
        print *, 'Writing swath to file: ', trim(swath)
        print *, 'l2gp%nFreqs:  ', l2gp%nFreqs
        print *, 'l2gp%nLevels: ', l2gp%nLevels
        print *, 'l2gp%nTimes:  ', l2gp%nTimes
        print *, 'shape(l2gp%l2gpvalue):  ', shape(l2gp%l2gpvalue)
        print *, 'l2gp%MissingL2GP:  ', l2gp%MissingL2GP
        print *, 'l2gp%Missingvalue:  ', l2gp%Missingvalue
      endif
      ! call dump ( l2gp%chunkNumber, 'chunkNumber as Read' )
      ! chunk = FindFirst( l2gp%chunkNumber /= -999 )
      ! call outputNamedValue ( 'index of non-Fill chunkNumber', chunk )
      ! Possibly repair or filter l2gp
      if ( repair ) call RepairL2GP( l2gp, HGrid, options=options )
      if ( filter ) call FilterL2GP( l2gp )
      if ( all(l2gp%chunkNumber == -999) ) &
        & call output( 'After repair and filtering, all chunk numbers are -999', advance='yes' )
      ! Possibly extract a reduced l2pg
      if ( present(rFreqs) .or. present(rLevels) .or. present(rTimes) ) then
        if ( present(rFreqs) ) then
          myFreqs = rFreqs
        else
          myFreqs(1) = 1
          myFreqs(2) = l2gp%nFreqs
        endif
        if ( present(rLevels) ) then
          myLevels = rLevels
        else
          myLevels(1) = 1
          myLevels(2) = l2gp%nLevels
        endif
        if ( present(rTimes) ) then
          myTimes = rTimes
        else
          myTimes(1) = 1
          myTimes(2) = l2gp%nTimes
        endif
        if ( DEEBUG ) call Dump( l2gp%BinNumber, 'Full l2gp BinNumber' )
        if ( DEEBUG ) call Dump( l2gp%MAF      , 'Full l2gp MAF' )
        call ExtractL2GPRecord ( l2gp, reducedl2gp, myFreqs, myLevels, myTimes )
        if ( DEEBUG ) call Dump( reducedl2gp%BinNumber, 'Reduced BinNumber' )
        if ( DEEBUG ) call Dump( reducedl2gp%MAF      , 'Reduced l2gp MAF' )
        call WriteL2GPData(reducedl2gp, file2, trim(swath2), &
          & hdfVersion=hdfVersion2, &
          & notUnlimited=notUnlimited)
        call DestroyL2GPContents ( reducedl2gp )
      else
      ! Write the filled l2gp to file2
        ! call dump ( l2gp%chunkNumber, 'chunkNumber as Written' )
        ! chunk = FindFirst( l2gp%chunkNumber /= -999 )
        ! call outputNamedValue ( 'index of non-Fill chunkNumber', chunk )
        call WriteL2GPData(l2gp, file2, trim(swath2), hdfVersion=hdfVersion2, &
          & notUnlimited=notUnlimited)
      endif

      l2metaData%minLat = mlsmin( l2gp%latitude, l2gp%MissingValue )
      l2metaData%maxLat = maxval(l2gp%latitude)
      l2metaData%minLon = mlsmin( l2gp%longitude, l2gp%MissingValue )
      l2metaData%maxLon = maxval(l2gp%longitude)

      call DestroyL2GPContents ( l2gp )
    enddo
       
  end subroutine cpL2GPData_fileID

  ! ---------------------- cpL2GPData_fileName  ---------------------------

  subroutine cpL2GPData_fileName(l2metaData, file1, file2, &
    & create2, hdfVersion1, hdfVersion2, swathList, rename, exclude, &
    & notUnlimited, andGlAttributes, ReadData, HGrid, &
    & rFreqs, rLevels, rTimes, options)
    !------------------------------------------------------------------------

    ! Given file names file1 and file2,
    ! This routine copies all the l2gpdata from 1 to 2
    ! (see cpL2GPData_fileID)
    ! If file2 doesn't exist yet, or if create2 is TRUE, it'll create it
    ! Optionally repairs l2gpdata

    use HGridsDatabase, only: HGrid_T
    ! Arguments

    type (L2Metadata_T) :: l2metaData
    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    logical, optional, intent(in) :: create2
    logical, optional, intent(in) :: andGlAttributes
    integer, optional, intent(in) :: hdfVersion1
    integer, optional, intent(in) :: hdfVersion2         !              is wild
    character (len=*), optional, intent(in) :: swathList ! Copy only these; '*'
    character (len=*), optional, intent(in) :: rename ! But rename them these
    character (len=*), optional, intent(in) :: exclude ! Don't copy these
    logical, optional, intent(in)           :: notUnlimited
    logical, optional, intent(in)           :: ReadData
    type (HGrid_T), optional, intent(in)    :: HGrid
    integer, dimension(2), intent(in), optional :: rFreqs  ! subscript range
    integer, dimension(2), intent(in), optional :: rLevels ! subscript range
    integer, dimension(2), intent(in), optional :: rTimes  ! subscript range
    character (len=*), optional, intent(in) :: options ! E.g., '-v'

    ! Local
    logical :: allSwaths
    logical, parameter            :: countEmpty = .true.
    ! logical, parameter :: DEEBUG = .TRUE.
    integer :: File1Handle
    integer :: File2Handle
    logical :: file_exists
    integer :: file_access
    integer :: listsize
    type(MLSFile_T)                :: MLSFile1, MLSFile2
    logical :: myandGlAttributes
    character (len=MAXSWATHNAMESBUFSIZE) :: mySwathList
    character (len=MAXSWATHNAMESBUFSIZE) :: myrename
    integer :: noSwaths
    integer :: status
    integer :: the_hdfVersion1
    integer :: the_hdfVersion2
    
    ! Executable code
    the_hdfVersion1 = L2GPDEFAULT_HDFVERSION
    if ( present(hdfVersion1) ) the_hdfVersion1 = hdfVersion1
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, trim(ModuleName) // '/cpL2GPData', &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    if ( the_hdfVersion1 == WILDCARDHDFVERSION ) then
      the_hdfVersion1 = mls_hdf_version(File1, the_hdfVersion1)
      if ( the_hdfVersion1 == FILENOTFOUND ) &
        call MLSMessage ( MLSMSG_Error, ModuleName, &
          & 'File 1 not found; make sure the name and path are correct' &
          & // trim(file1) )
    endif
    the_hdfVersion2 = the_hdfVersion1  ! Defaults to same hdf version as file 1
    if ( present(hdfVersion2) ) the_hdfVersion2 = hdfVersion2

    allSwaths = .not. present(swathList)
    if ( present(swathList) ) allSwaths = (swathList == '*')
    ! if ( present(swathList) ) then
    if ( .not. allSwaths ) then
      if ( DEEBUG ) then
        noSwaths = mls_InqSwath ( file1, mySwathList, listSize, &
           & hdfVersion=the_hdfVersion1)
        print *, 'swathList you requested to cp: ', trim(swathList)
        print *, 'mls_InqSwath finds: ', trim(mySwathList)
      endif
      mySwathList = swathList
    else
      noSwaths = mls_InqSwath ( file1, mySwathList, listSize, &
           & hdfVersion=the_hdfVersion1)
    endif
    status = InitializeMLSFile ( MLSFile1, type=l_swath, access=DFACC_READ, &
     & name=trim(file1), HDFVersion=the_hdfVersion1 )
    call MLS_OpenFile( MLSFile1 )
    File1Handle = MLSFile1%FileID%f_id

    file_exists = ( mls_exists(trim(File2)) == 0 )
    if ( file_exists ) then
      file_access = DFACC_RDWR
    else
      file_access = DFACC_CREATE
    endif
    if ( present(create2) ) then
      if ( create2 ) then
        file_access = DFACC_CREATE
      elseif ( .not. file_exists ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
         & "L2gp file: " // trim(File2) // ' for cp-ing not yet existing')
        file_access = DFACC_RDWR
      else
        file_access = DFACC_RDWR
      endif
    endif
    if ( DEEBUG ) then
      print *, 'About to open file2: ', trim(file2)
      print *, 'file_access: ', file_access
      print *, 'hdfVersion: ', the_hdfVersion2
    endif
    myandGlAttributes = (file_access == DFACC_CREATE)
    if ( present(andGlAttributes) ) &
      & myandGlAttributes = myandGlAttributes .and. andGlAttributes
    myandGlAttributes = myandGlAttributes .and. &
      & (the_hdfVersion1 == HDFVERSION_5) .and. (the_hdfVersion2 == HDFVERSION_5)
    status = InitializeMLSFile ( MLSFile2, type=l_swath, access=file_access, &
     & name=trim(file2), HDFVersion=the_hdfVersion2 )
    call MLS_OpenFile( MLSFile2 )
    File2Handle = MLSFile2%FileID%f_id
    if ( DEEBUG ) then
      print *, 'About to cp from file1 to file2: ', trim(file1), trim(file2)
      print *, trim(mySwathList)
    endif
    ! Maybe copy global attributes, too
    if ( myandGlAttributes ) then
      call cpHE5GlobalAttrs ( File1Handle, File2Handle, status )
      if ( status /= 0 ) &
        & call output ( '(Global Attributes missing) ' // trim(file1), advance='yes')
    endif

    if ( present(exclude) ) then
      call RemoveListFromList( mySwathList, myrename, exclude )
      noSwaths = NumStringElements( trim(myrename), countEmpty )
      if ( noSwaths < 1 ) then
        call MLSMessage ( MLSMSG_Warning, ModuleName, &
              & 'No swaths after excluding ' )
        call dump( trim(exclude) , 'exclude these' )
        call dump( trim(mySwathList) , 'original SwathList' )
        call dump( trim(myrename) , 'remaining' )
        if ( DEEBUG ) then
          call output( 'Lets do it elem-by-elm', advance='yes' )
          call RemoveListFromList( mySwathList, myrename, exclude, options='-v' )
        endif
        call MLS_CloseFile( MLSFile1 )
        call MLS_CloseFile( MLSFile2 )
        return
      endif
      mySwathList = myrename
    endif
    call cpL2GPData_fileID(l2metaData, File1Handle, File2Handle, &
      & mySwathList, the_hdfVersion1, the_hdfVersion2, &
      & notUnlimited=notUnlimited, rename=rename, &
      & ReadData=ReadData, HGrid=HGrid, &
      & rFreqs=rFreqs, rLevels=rlevels, rTimes=rTimes, options=options )
    if ( DEEBUG ) print *, 'About to close File1Handle: ', File1Handle
    call MLS_CloseFile( MLSFile1 )
    if ( DEEBUG ) print *, 'About to close File2Handle: ', File2Handle
    call MLS_CloseFile( MLSFile2 )
  end subroutine cpL2GPData_fileName

  ! ---------------------- cpL2GPData_MLSFile  ---------------------------


  subroutine cpL2GPData_MLSFile(l2metaData, L2GPfile1, L2GPfile2, &
    & create2, swathList, rename, exclude, &
    & notUnlimited, andGlAttributes, ReadData, HGrid, &
    & rFreqs, rLevels, rTimes, options)
    !------------------------------------------------------------------------

    ! Given MLSFiles L2GPfile1 and L2GPfile2,
    ! This routine copies all the l2gpdata from 1 to 2
    ! (see cpL2GPData_fileID)
    ! If L2GPfile2 doesn't exist yet, or if create2 is TRUE, it'll create it
    ! Optionally repairs l2gpdata
    use HGridsDatabase, only: HGrid_T
    ! Arguments

    type (L2Metadata_T) :: l2metaData
    type(MLSFile_T)               :: L2GPfile1 ! file 1
    type(MLSFile_T)               :: L2GPfile2 ! file 2
    logical, optional, intent(in) :: create2
    logical, optional, intent(in) :: andGlAttributes
    character (len=*), optional, intent(in) :: swathList ! Copy only these; '*'
    character (len=*), optional, intent(in) :: rename ! But rename them these
    character (len=*), optional, intent(in) :: exclude ! Don't copy these
    logical, optional, intent(in)           :: notUnlimited
    logical, optional, intent(in)           :: ReadData
    type (HGrid_T), optional, intent(in)    :: HGrid
    integer, dimension(2), intent(in), optional :: rFreqs  ! subscript range
    integer, dimension(2), intent(in), optional :: rLevels ! subscript range
    integer, dimension(2), intent(in), optional :: rTimes  ! subscript range
    character (len=*), optional, intent(in) :: options ! E.g., '-v'
    ! Local
    logical :: L1alreadyOpen
    logical :: L2alreadyOpen
    logical :: myCreate2
    character(len=MAXSWATHNAMESBUFSIZE) :: mySwathList
    integer :: originalAccess
    integer :: status
    !
    myCreate2 = .false.
    if ( present(create2) ) myCreate2 = create2
    mySwathList = '*'
    if ( present(swathList) ) mySwathList = swathList
    originalAccess = L2GPFile2%access
    status = 0
    L1alreadyOpen = L2GPFile1%stillOpen
    if ( L1alreadyOpen ) then
      call mls_closeFile(L2GPFile1, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to close l2gp file to cp from', MLSFile=L2GPFile1)
    endif

    status = 0
    L2alreadyOpen = L2GPFile2%stillOpen
    if ( L2alreadyOpen ) then
      if ( myCreate2 ) L2GPFile2%access = DFACC_CREATE
      call mls_closeFile(L2GPFile2, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to close l2gp file to cp to', MLSFile=L2GPFile2)
    endif
    
    call cpL2GPData_fileName(l2metaData, L2GPFile1%Name, L2GPFile2%Name, &
      & create2, L2GPFile1%hdfVersion, L2GPFile2%hdfVersion, &
      & SwathList, rename, exclude, &
      & notUnlimited, andGLAttributes, ReadData, HGrid, &
      & rFreqs, rLevels, rTimes, options)

    if ( L1alreadyOpen )  call mls_openFile(L2GPFile1, Status)
    L2GPFile1%errorCode = status
    L2GPFile1%lastOperation = 'read'

    if ( L2alreadyOpen )  call mls_openFile(L2GPFile2, Status)
    L2GPFile2%errorCode = status
    L2GPFile2%lastOperation = 'write'
    L2GPFile2%access = originalAccess

  end subroutine cpL2GPData_MLSFile

  ! ---------------------- cpL2GPDataToAttribute  ---------------------------

  subroutine cpL2GPDataToAttribute( L2GPfile1, L2GPfile2, &
    & swathname, attrname )
    !------------------------------------------------------------------------
    use HDFEOS5, only: HE5T_Native_Real
    use MLSHDFEOS, only: HE5_EhwrGlAtt, HSize

    ! Given MLSFiles L2GPfile1 and L2GPfile2,
    ! This routine copies the l2gpdata named swathname from 1 to the
    ! file level attribute named attrname in 2
    ! Arguments

    type (L2Metadata_T) :: l2metaData
    type(MLSFile_T)               :: L2GPfile1 ! file 1
    type(MLSFile_T)               :: L2GPfile2 ! file 2
    character (len=*), optional, intent(in) :: swathname ! Name to copy
    character (len=*), optional, intent(in) :: attrname ! But rename it this
    ! Local
    type (L2GPData_T) :: l2gp
    integer :: status
    !
    call ReadL2GPData ( L2GPfile1, trim(swathname), l2gp )
    if ( .not. L2GPfile2%stillOpen ) call MLS_OpenFile( L2GPfile2 )
    status = he5_ehwrglatt( L2GPfile2%fileID%f_id, &
            & trim(attrname), HE5T_NATIVE_REAL, hsize(l2gp%nTimes), &
            &  l2gp%l2gpValue(1,1,:) )
    call MLS_CloseFile( L2GPfile2 )
  end subroutine cpL2GPDataToAttribute

  ! ------------------------------------------ DiffL2GPData_CHUNKS ------------

  subroutine DiffL2GPData_CHUNKS ( fullL2gp1, fullL2gp2, Chunks, &
    & Details, options, fields, numDiffs )
    ! Diff selected chunks of an l2gp
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          FULLL2GP1
    type (l2gpData_T), intent(in) ::          FULLL2GP2
    integer, intent(in), dimension(:) :: Chunks ! Which chunks to diff
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: fields  ! diff only these fields
    integer, intent(out), optional :: numDiffs  ! how many diffs

    ! Local variables
    integer :: chunk
    ! logical, parameter :: DEEBUG = .true.
    integer :: i
    integer, dimension(2) :: irange
    type (l2gpData_T) ::          L2GP1, L2gp2
    logical :: mySilent
    ! Executable
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0
    if ( DEEBUG ) then
      call output( 'mySilent: ', advance='no' )
      call output( mySilent, advance='yes' )
      call output( 'size(chunks): ', advance='no' )
      call output( size(chunks), advance='yes' )
      call output( chunks, advance='yes' )
      call output( 'fullL2gp1%chunkNumber: ', advance='no' )
      call output( fullL2gp1%chunkNumber, advance='yes' )
    endif
    do i=1, size(Chunks)
      chunk = Chunks(i)
      irange(1) = FindFirst( fullL2gp1%chunkNumber, chunk )
      irange(2) = Findlast( fullL2gp1%chunkNumber, chunk )
      ! HuntRange does not work correctly when any chunk Numbers are -999
      ! (Missing Value)
      if ( any( irange == 0 ) ) cycle
      if ( .not. mySilent ) then
        call output ( ' - - - Chunk number:', advance='no')
        call output ( chunk, advance='no')
        call output ( ' - - -', advance='yes')
      endif
      call ExtractL2GPRecord ( fullL2gp1, l2gp1, rTimes=irange )
      call ExtractL2GPRecord ( fullL2gp2, l2gp2, rTimes=irange )
      call Diff ( L2gp1, L2gp2, &
         & Details, options, fields, &
         & numDiffs=numDiffs )
      call DestroyL2GPContents( l2gp1 )
      call DestroyL2GPContents( l2gp2 )
    enddo
  end subroutine DiffL2GPData_CHUNKS

  ! ------------------------------------------ DiffL2GPData_LEVELS ------------

  subroutine DiffL2GPData_LEVELS ( fullL2gp1, fullL2gp2, Pressures, Chunks, &
    & Details, options, fields, numDiffs )
    ! DDiff selected levels of an l2gp
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          FULLL2GP1
    type (l2gpData_T), intent(in) ::          FULLL2GP2
    real(rgp), intent(in), dimension(:) :: pressures ! Which levels to diff
    integer, intent(in), optional, dimension(:) :: Chunks ! Which chunks to diff
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: fields  ! diff only these fields
    integer, intent(out), optional :: numDiffs  ! how many diffs

    ! Local variables
    real(rgp) :: pressure
    ! logical, parameter :: DEEBUG = .true.
    integer :: i
    integer, dimension(2) :: irange
    type (l2gpData_T) ::          L2GP1, L2gp2
    logical :: mySilent
    ! Executable
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0
    if ( DEEBUG ) then
      call output( 'mySilent: ', advance='no' )
      call output( mySilent, advance='yes' )
      call output( 'size(pressures): ', advance='no' )
      call output( size(pressures), advance='yes' )
      call output( pressures, advance='yes' )
      call output( 'fullL2gp1%pressures: ', advance='no' )
      call output( fullL2gp1%pressures, advance='yes' )
    endif
    do i=1, size(pressures)
      pressure = pressures(i)
      ! The dubious "1.005" factor is to prevent missing a desired level
      ! due to unlucky roundoff
      irange(1) = FindFirst( fullL2gp1%pressures <= pressure*1.005 )
      irange(2) = irange(1)
      if ( any( irange == 0 ) ) cycle
      if ( .not. mySilent ) then
        call output ( ' - - - pressure:', advance='no')
        call output ( fullL2gp1%pressures(irange(1)), advance='no')
        call output ( ' - - -', advance='yes')
      endif
      call ExtractL2GPRecord ( fullL2gp1, l2gp1, rLevels=irange )
      call ExtractL2GPRecord ( fullL2gp2, l2gp2, rLevels=irange )
      if ( present(chunks) ) then
        call DiffL2GPData_Chunks ( L2gp1, L2gp2, chunks, &
          & Details, options, fields, numDiffs=numDiffs )
      else
        call Diff ( L2gp1, L2gp2, &
          & Details, options, fields, numDiffs=numDiffs )
      endif
      call DestroyL2GPContents( l2gp1 )
      call DestroyL2GPContents( l2gp2 )
    enddo
  end subroutine DiffL2GPData_LEVELS

  ! ------------------------------------------ DiffL2GPData_RANGES ------------

  subroutine DiffL2GPData_RANGES ( fullL2gp1, fullL2gp2, &
    & geoBoxNames, geoBoxLowBound, geoBoxHiBound, &
    & Pressures, Latitudes, Longitudes, Times, Chunks, &
    & Details, options, fields, numDiffs )
    ! Diff 2 l2gps according to prescribed ranges of pressure, longitude, etc.
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          FULLL2GP1
    type (l2gpData_T), intent(in) ::          FULLL2GP2
    character(len=*), intent(in) , optional ::       geoBoxNames
    real(rgp), dimension(:), intent(in), optional :: geoBoxLowBound  ! range
    real(rgp), dimension(:), intent(in), optional :: geoBoxHiBound  ! range
    real(rgp), intent(in), dimension(2), optional :: pressures ! range
    real(rgp), intent(in), dimension(2), optional :: latitudes ! range
    real(rgp), intent(in), dimension(2), optional :: longitudes ! range
    real(r8), intent(in), dimension(2), optional :: times ! range
    integer, intent(in), optional, dimension(:) :: Chunks ! Which chunks to diff
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: fields  ! diff only these fields
    integer, intent(out), optional :: numDiffs  ! how many diffs

    ! Local variables
    ! logical, parameter :: DEEBUG = .true.
    type (l2gpData_T) ::          L2GP1, L2gp2
    logical :: mySilent
    ! Executable
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0
    if ( DEEBUG ) then
      call output( 'mySilent: ', advance='no' )
      call output( mySilent, advance='yes' )
      if ( present(pressures) ) &
        & call dump( pressures, 'pressures' )
      if ( present(latitudes) ) &
        & call dump( latitudes, 'latitudes' )
      if ( present(longitudes) ) &
        & call dump( longitudes, 'longitudes' )
      call output( 'fullL2gp1%pressures: ', advance='no' )
      call output( fullL2gp1%pressures, advance='yes' )
    endif
    if ( present(geoBoxNames) ) then
      call ContractL2GPRecord ( fullL2gp1, l2gp1, &
        & geoBoxNames, geoBoxLowBound, geoBoxHiBound )
      call ContractL2GPRecord ( fullL2gp2, l2gp2, &
        & geoBoxNames, geoBoxLowBound, geoBoxHiBound )
    else
      call ContractL2GPRecord ( fullL2gp1, l2gp1, pressures=pressures, &
        & latitudes=latitudes, longitudes=longitudes, intimes=times )
      call ContractL2GPRecord ( fullL2gp2, l2gp2, pressures=pressures, &
        & latitudes=latitudes, longitudes=longitudes, intimes=times )
    endif
    if ( present(chunks) ) then
      call DiffL2GPData_Chunks ( L2gp1, L2gp2, chunks, &
        & Details, options, fields, numDiffs=numDiffs )
    else
      call Diff ( L2gp1, L2gp2, &
        & Details, options, fields, numDiffs=numDiffs )
    endif
    call DestroyL2GPContents( l2gp1 )
    call DestroyL2GPContents( l2gp2 )
  end subroutine DiffL2GPData_RANGES

  ! ---------------------- DiffL2GPData  ---------------------------
  subroutine DiffL2GPData ( L2gp1, L2gp2, &
    & Details, options, fields, matchTimes, numDiffs )
    ! Show diff between l2gp1 and l2gp2 down to level of Details
    
    ! Note:
    ! by default, print arrays of diffs, but not stats nor rms
    ! If either stats or rms, just print those, but not arrays
    ! To print stats and arrays, both, turn on wholeArrays and stats
    ! To print stats, rms, and arrays, both, turn on wholeArrays, rms, stats
    
    ! (Is this too complicated? Should we make wholeArrays always on
    !  by default unless explicitly turned off?)
    ! Dummy arguments
    type (l2gpData_T), intent(inout) ::          L2GP1
    type (l2gpData_T), intent(inout) ::          L2GP2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: fields  ! diff only these fields
    logical, intent(in), optional :: matchTimes
    integer, intent(out), optional :: numDiffs  ! how many diffs
    ! Local variables
    integer :: how_many
    logical :: myMatchTimes
    logical :: mySilent
    type (L2GPData_T) :: cl2gp1, cl2gp2 ! Compacted l2gps
    type (L2GPData_T) :: tl2gp1, tl2gp2 ! Synchronized l2gps
    integer, dimension(MAXNUMTIMES) :: which1, which2
    integer :: i
    ! Executable
    myMatchTimes = .false.
    if ( present(matchTimes) ) mymatchTimes = matchTimes
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0
    if ( present(numDiffs) ) numDiffs = 0
    if ( DEEBUG ) call outputNamedValue( 'match times?', myMatchTimes )
    if ( myMatchTimes .or. L2gp1%nTimes /= L2gp2%nTimes ) then
      call CompactL2GPRecord ( L2gp1, cL2gp1 )
      call CompactL2GPRecord ( L2gp2, cL2gp2 )
      ! 1st, find which profile times match (to within 0.1 s)
      call FindIntersection( cl2gp1%time, cl2gp2%time, which1, which2, how_many, &
        & tol=1.d-1 )
      if ( .not. mySilent ) then
        ! call dump( cl2gp1 )
        ! call dump( cl2gp2 )
        call OutputNamedValue ( 'how many times matched', how_many )
        call dump( which1(1:how_many), 'which1' )
        call dump( which2(1:how_many), 'which2' )
        do i = 1, how_many
          write(*,*) timeToHoursInDay ( cl2gp1%time(which1(i)) ), &
            & timeToHoursInDay ( cl2gp2%time(which2(i)) )
        enddo
      endif
      if ( how_many == 0 ) then
        ! 2nd chance--try to match Geod. Ang. to within 1/2 deg.
        if ( .not. mySilent ) call output( &
          & 'No matching times found in l2gp1, l2gp2', advance='yes' )
        call FindIntersection( cl2gp1%GeodAngle, cl2gp2%GeodAngle, which1, which2, how_many, &
          & tol=0.5 )
        if ( how_many == 0 ) then
          if ( .not. mySilent ) call output( &
            & 'No matching geod. angs found in l2gp1, l2gp2', advance='yes' )
          return
        endif
      endif
      call SetupNewL2GPRecord ( tl2gp1, proto=cl2gp1, which=which1(1:how_many) )
      call SetupNewL2GPRecord ( tl2gp2, proto=cl2gp2, which=which2(1:how_many) )
      if ( DEEBUG ) print *, 'About to enter ..atLast having matched times'
      call DiffThis ( tL2gp1, tL2gp2, &
      & Details, options, fields, numDiffs )
      call DestroyL2GPContents( cL2gp1 )
      call DestroyL2GPContents( cL2gp2 )
      call DestroyL2GPContents( tL2gp1 )
      call DestroyL2GPContents( tL2gp2 )
    else
      if ( DEEBUG ) print *, 'Diffing as we found them'
      call DiffThis ( L2gp1, L2gp2, &
      & Details, options, fields, numDiffs )
    endif
  end subroutine DiffL2GPData

  ! ---------------------- DiffThis  ---------------------------
  subroutine DiffThis ( L2gp1, L2gp2, &
    & Details, options, fields, numDiffs )
    ! Show diff between l2gp1 and l2gp2 down to level of Details
    
    ! Note:
    ! by default, print arrays of diffs, but not stats nor rms
    ! If either stats or rms, just print those, but not arrays
    ! To print stats and arrays, both, turn on wholeArrays and stats
    ! To print stats, rms, and arrays, both, turn on wholeArrays, rms, stats
    
    ! (Is this too complicated? Should we make wholeArrays always on
    !  by default unless explicitly turned off?)
    ! Dummy arguments
    type (l2gpData_T), intent(inout) ::          L2GP1
    type (l2gpData_T), intent(inout) ::          L2GP2
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    character(len=*), intent(in), optional :: options
    character(len=*), intent(in), optional :: fields  ! diff only these fields
    integer, intent(out), optional :: numDiffs  ! how many diffs
    ! Local variables
    logical :: badChunks              ! Are any chunks bad?
    logical :: badChunksDifferent     ! Are any chunks bad in one l2gp not bad in the other?
    integer :: badInstances
    character(len=*), parameter :: DEFAULTFIELDS = &
      & 'pressures, latitude, longitude, solarTime, solarZenith,' // &
      & 'losAngle, geodAngle, time, chunkNumber, frequency,'  // &
      & 'l2gpvalue, l2gpPrecision, status, quality, convergence, AscDescMode '
    ! logical, parameter :: DEEBUG = .true.
    logical :: diffGeosMeanBadChunks
    logical :: fillsInL2GP1
    logical :: fillsInL2GP2
    logical :: goldbrick
    integer :: instance
    type (l2gpData_T) ::          L2GP2Temp
    integer :: MYDETAILS
    character(len=len(DEFAULTFIELDS)) :: myFields
    integer :: myNumDiffs
    logical :: mySilent
    logical :: myVerbose
    logical :: ShapesDontMatch
    logical :: skipGeos
    logical :: ForceThrough
    ! Executable code
    nameOnEachLine = ' '
    if ( .not. isL2GPSetUp(l2gp1) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, & 
        'l2gp1 not yet allocated in diff')
      return
    elseif ( .not. isL2GPSetUp(l2gp2) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, & 
        'l2gp2 not yet allocated in diff')
      return
    endif
    myDetails = 1
    if ( present(details) ) myDetails = details
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0
    if ( mySilent ) call suspendOutput
    myVerbose = .false.
    if ( present(options) ) myVerbose = index(options, 'v') > 0
    goldbrick = .false.
    if ( present(options) ) goldbrick = index(options, '@') > 0
    ForceThrough = .false.
    if ( present(options) ) ForceThrough = index(options, 'F') > 0

    myNumDiffs = 0

    ShapesDontMatch = .false.
    badChunks = .false.
    badInstances = 0
    diffGeosMeanBadChunks = .false.
    myFields = DEFAULTFIELDS
    skipGeos = .false.
    if ( present(fields) ) myFields = fields
                                                                              
    if ( present ( options ) ) then                                           
      if ( DEEBUG ) call output( 'options: ' // trim(options), advance='yes' )
      badChunks = (index(options, 'i' ) > 0) .and. badChunks                  
      diffGeosMeanBadChunks = ( index(options, 'd') > 0 )                     
      skipGeos = ( index(options, 'g') > 0 )                                  
    else                                                                      
      badChunks = .false.                                                     
    endif                                                                     

    if ( DEEBUG ) then                                                        
      call output( 'myDetails: ', advance='no' )
      call output( myDetails, advance='yes' )
      call output( 'mySilent: ', advance='no' )
      call output( mySilent, advance='yes' )
      call output( 'myFields: ', advance='no' )
      call output( myFields, advance='yes' )   
      call output( 'skipGeos: ', advance='no' )
      call output( skipGeos, advance='yes' )   
    endif
    if ( myFields == '*' .or. lowercase(myFields) == 'all' ) &
      & myFields = DEFAULTFIELDS
    if ( DEEBUG ) then
      call output( 'myFields: ', advance='no' )
      call output( myFields, advance='yes' )
    endif
    if ( trim(l2gp1%name) /= trim(l2gp2%name) ) then
      call output('(1) name: ' // trim(l2gp1%name), advance='yes')
      call output('(2) name: ' // trim(l2gp2%name), advance='yes')
      myNumDiffs = myNumDiffs + 1
    endif
    if ( myDetails < -1 ) then
      call doneHere
      return
    endif
    if ( L2gp1%MissingL2GP /= L2gp2%MissingL2GP ) then
      call output(' (1) MissingL2GP = ', advance='no')
      call output(L2gp1%Missingl2GP, advance='yes')
      call output(' (2) Missingl2GP = ', advance='no')
      call output(L2gp2%MissingL2GP, advance='yes')
      myNumDiffs = myNumDiffs + 1
    endif
    if ( L2gp1%MissingValue /= L2gp2%MissingValue ) then
      call output(' (1) MissingValue = ', advance='no')
      call output(L2gp1%MissingValue, advance='yes')
      call output(' (2) MissingValue = ', advance='no')
      call output(L2gp2%MissingValue, advance='yes')
      myNumDiffs = myNumDiffs + 1
    endif
    if ( L2gp1%nTimes /= L2gp2%nTimes ) then
      call output(' (1) nTimes = ', advance='no')
      call output(L2gp1%nTimes, advance='yes')
      call output(' (2) nTimes = ', advance='no')
      call output(L2gp2%nTimes, advance='yes')
      ShapesDontMatch = .true.
      myNumDiffs = myNumDiffs + 1
    endif
    if ( L2gp1%nLevels /= L2gp2%nLevels ) then
      call output(' (1) nLevels = ', advance='no')
      call output(L2gp1%nLevels, advance='yes')
      call output(' (2) nLevels = ', advance='no')
      call output(L2gp2%nLevels, advance='yes')
      ShapesDontMatch = .true.
      myNumDiffs = myNumDiffs + 1
    endif
    if ( L2gp1%nFreqs /= L2gp2%nFreqs ) then
      call output(' (1) nFreqs = ', advance='no')
      call output(L2gp1%nFreqs, advance='yes')
      call output(' (2) nFreqs = ', advance='no')
      call output(L2gp2%nFreqs, advance='yes')
      ShapesDontMatch = .true.
      myNumDiffs = myNumDiffs + 1
    endif
    if ( myDetails < 0 )  then
      if ( myVerbose ) &
        & call output ( '(Data shapes and swath names match)', advance='yes' )
      call doneHere
      return
    endif
    if ( ShapesDontMatch  ) then
      if ( ForceThrough ) then
        call output('Some fields wont be diffed because shapes dont match',&
          & advance='yes')
      else
        call output('Skipping further details because shapes dont match',&
          & advance='yes')
        call doneHere
        return
      endif
    endif

    badChunks = badChunks .or. &
      & ( &
      & any( mod ( l2gp1%status, 2 ) == 1 ) &
      & .or. &
      & any( mod ( l2gp2%status, 2 ) == 1 ) &
      & .or. &
      & any( IsFillValue ( l2gp1%l2gpValue, l2gp1%MissingL2GP ) ) &
      & .or. &
      & any( IsFillValue ( l2gp2%l2gpValue, l2gp2%MissingL2GP ) ) &
      & )
    ! OK, we'll try what you suggest
    call SetupNewL2GPRecord ( l2gp2Temp, proto=l2gp2 )
    l2gp2Temp%status = l2gp2%status
    if ( diffGeosMeanBadChunks ) then
      do instance=1, L2gp1%nTimes
        if ( l2gp1%geodAngle(instance) /= l2gp2%geodAngle(instance) &
          & .or. &
          & l2gp1%solarZenith(instance) /= l2gp2%solarZenith(instance) &
          & .or. &
          & l2gp1%time(instance) /= l2gp2%time(instance) ) &
          & l2gp2Temp%status(instance) = -1
      enddo
    endif
    badChunksDifferent = .false.
    if ( badChunks ) then
      do instance=1, L2gp1%nTimes
        fillsInL2GP1 = any( &
          & IsFillValue ( l2gp1%l2gpValue(:,:,instance), l2gp1%MissingL2GP ) &
          & )
        fillsInL2GP2 = any( &
          & IsFillValue ( l2gp2%l2gpValue(:,:,instance), l2gp2%MissingL2GP ) &
          & )
        if  ( l2gp2Temp%status(instance) < 0 &
          & .or. &
          & mod(l2gp1%status(instance), 2) == 1 &
          & .or. &
          & mod(l2gp2%status(instance), 2) == 1 &
          & .or. &
          & fillsInL2GP1 &
          & .or. &
          & fillsInL2GP2 &
          & ) then
          badChunksDifferent = badChunksDifferent .or. &
            mod(l2gp1%status(instance), 2) /= mod(l2gp2%status(instance), 2) &
            & .or. &
            & fillsInL2GP1 .neqv. fillsInL2GP2
          l2gp2Temp%l2gpValue(:,:,instance) = l2gp1%l2gpValue(:,:,instance)
          l2gp2Temp%l2gpPrecision(:,:,instance) = l2gp1%l2gpPrecision(:,:,instance)
          l2gp2Temp%status(instance) = l2gp1%status(instance)
          l2gp2Temp%quality(instance) = l2gp1%quality(instance)
          l2gp2Temp%convergence(instance) = l2gp1%convergence(instance)
          if ( AscDescModeIsField ) &
            & l2gp2Temp%AscDescMode(instance) = l2gp1%AscDescMode(instance)
          ! Also geolocations will be reset to avoid displaying bogus diffs
          l2gp2Temp%latitude(instance) = l2gp1%latitude(instance)
          l2gp2Temp%longitude(instance) = l2gp1%longitude(instance)
          l2gp2Temp%solarTime(instance) = l2gp1%solarTime(instance)
          l2gp2Temp%chunkNumber(instance) = l2gp1%chunkNumber(instance)
          l2gp2Temp%geodAngle(instance) = l2gp1%geodAngle(instance)
          l2gp2Temp%solarZenith(instance) = l2gp1%solarZenith(instance)
          l2gp2Temp%time(instance) = l2gp1%time(instance)
          l2gp2Temp%losAngle(instance) = l2gp1%losAngle(instance)
          badInstances = badInstances + 1
        endif
      enddo
      if ( badChunksDifferent ) then
        call output('Number of bad instances of l2gp2 reset to l2gp1 ', advance='no')
        call output(badInstances, advance='yes')
      endif
    endif

    if ( .not. skipGeos ) then
      ! call output( 'calling diffGeoLocations', advance='yes' )
      ! call outputNamedValue( 'nameOnEachLine', trim(nameOnEachLine) )
      call diffGeoLocations( l2gp1, l2gp2Temp )
    endif
    if ( myDetails < 1 .or. ShapesDontMatch )  then
      call doneHere
      return
    endif

    NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%values', ' ', statsOnOneLine )
    if ( any(l2gp1%l2gpValue /= l2gp2Temp%l2gpValue) .and. &
      & SwitchDetail(lowercase(myFields), 'l2gpvalue', '-fc') > -1 ) then
      call diff ( l2gp1%l2gpValue, 'l2gp%l2gpValue', &
        &         l2gp2Temp%l2gpValue, ' ', &
        & options=options, fillValue=l2gp1%MissingL2GP )
      myNumDiffs = myNumDiffs + count( l2gp1%l2gpValue /= l2gp2Temp%l2gpValue )
    elseif ( all(l2gp1%l2gpValue == l2gp2Temp%l2gpValue) .and. &
      & SwitchDetail(lowercase(myFields), 'l2gpvalue', '-fc') > -1 .and. myVerbose ) then
      call output('(values fields equal)', advance='yes')
    endif

    NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%Precisions', ' ', statsOnOneLine )
    if ( any(l2gp1%l2gpPrecision /= l2gp2Temp%l2gpPrecision) .and. &
      & SwitchDetail(lowercase(myFields), 'l2gpprecision', '-fc') > -1 ) then
      call diff ( l2gp1%l2gpPrecision, 'l2gpPrecision', &
        &         l2gp2Temp%l2gpPrecision, ' ', &
        & options=options, fillValue=l2gp1%MissingValue )
      myNumDiffs = myNumDiffs + count( l2gp1%l2gpPrecision /= l2gp2Temp%l2gpPrecision )
    elseif ( all(l2gp1%l2gpPrecision == l2gp2Temp%l2gpPrecision) .and. &
      & SwitchDetail(lowercase(myFields), 'l2gpprecision', '-fc') > -1 .and. myVerbose ) then
      call output('(precision fields equal)', advance='yes')
    endif

    NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%status', ' ', statsOnOneLine )
    if ( any(l2gp1%status /= l2gp2Temp%status) .and. &
      & SwitchDetail(lowercase(myFields), 'status', '-fc') > -1 ) then
      call diff ( l2gp1%status, 'status', &
        &         l2gp2Temp%status, ' ', &
        & options=options )
      myNumDiffs = myNumDiffs + count( l2gp1%status /= l2gp2Temp%status )
    elseif ( all(l2gp1%status == l2gp2Temp%status) .and. &
      & SwitchDetail(lowercase(myFields), 'status', '-fc') > -1 .and. myVerbose ) then
      call output('(status fields equal)', advance='yes')
    endif

    NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%quality', ' ', statsOnOneLine )
    if ( any(l2gp1%quality /= l2gp2Temp%quality) .and. &
      & SwitchDetail(lowercase(myFields), 'quality', '-fc') > -1 ) then
      call diff ( l2gp1%quality, 'quality', &
        &         l2gp2Temp%quality, ' ', &
        & options=options )
      myNumDiffs = myNumDiffs + count( l2gp1%quality /= l2gp2Temp%quality )
    elseif ( all(l2gp1%quality == l2gp2Temp%quality) .and. &
      & SwitchDetail(lowercase(myFields), 'quality', '-fc') > -1 .and. myVerbose ) then
      call output('(quality fields equal)', advance='yes')
    endif

    NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%convergence', ' ', statsOnOneLine )
    if ( any(l2gp1%convergence /= l2gp2Temp%convergence) .and. &
      & SwitchDetail(lowercase(myFields), 'convergence', '-fc') > -1 ) then
      call diff ( l2gp1%convergence, 'convergence', &
        &         l2gp2Temp%convergence, ' ', &
        & options=options )
      myNumDiffs = myNumDiffs + count( l2gp1%convergence /= l2gp2Temp%convergence )
    elseif ( all(l2gp1%convergence == l2gp2Temp%convergence) .and. &
      & SwitchDetail(lowercase(myFields), 'convergence', '-fc') > -1 .and. myVerbose ) then
      call output('(convergence fields equal)', advance='yes')
    endif
    call doneHere
    return

  contains
    subroutine diffGeoLocations( l2gp1, l2gp2 )
      ! Args
      type(L2GPData_T) :: l2gp1
      type(L2GPData_T) :: l2gp2
      ! Executable
      if ( all(shape(l2gp1%pressures) == shape(l2gp2%pressures)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%pressure', ' ', statsOnOneLine )
        if ( any(l2gp1%pressures /= l2gp2%pressures) .and. &
          & SwitchDetail(lowercase(myFields), 'pressure', '-fc') > -1 ) then
            call diff ( l2gp1%pressures, 'l2gp%pressures', &
              &         l2gp2%pressures, ' ', &
              & options=options )
          myNumDiffs = myNumDiffs + count( l2gp1%pressures /= l2gp2%pressures )
        endif
      endif
      if ( all(shape(l2gp1%latitude) == shape(l2gp2%latitude)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%latitude', ' ', statsOnOneLine )
        if ( any(l2gp1%latitude /= l2gp2%latitude) .and. &
          & SwitchDetail(lowercase(myFields), 'lat', '-fc') > -1 ) then
            call diff ( l2gp1%latitude, 'latitude', &
              &         l2gp2%latitude, ' ', &
              & options=options )
          myNumDiffs = myNumDiffs + count( l2gp1%latitude /= l2gp2%latitude )
        endif
      endif
      if ( all(shape(l2gp1%longitude) == shape(l2gp2%longitude)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%longitude', ' ', statsOnOneLine )
        if ( any(l2gp1%longitude /= l2gp2%longitude) .and. &
          & SwitchDetail(lowercase(myFields), 'lon', '-fc') > -1 ) then
          if ( goldbrick ) then
            call diff ( l2gp1%longitude, 'longitude', &
              &         l2gp2%longitude, ' ', &
              & options=options )
          else
            call dump ( &
              & diff_fun( l2gp1%longitude, l2gp2%longitude, &
              &         auxvalue=360._rgp, &
              &         options=AddOntoOptions('p',options) ), &
              & 'longitude', &
              & options=options )
          endif
          myNumDiffs = myNumDiffs + count( l2gp1%longitude /= l2gp2%longitude )
        endif
      endif
      if ( all(shape(l2gp1%solarTime) == shape(l2gp2%solarTime)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%solartime', ' ', statsOnOneLine )
        if ( any(l2gp1%solarTime /= l2gp2%solarTime) .and. &
          & SwitchDetail(lowercase(myFields), 'solartime', '-fc') > -1 ) then
          if ( goldbrick ) then
            call diff ( l2gp1%solarTime, 'solarTime', &
              &         l2gp2%solarTime, ' ', &
              & options=options )
          else
            call dump ( &
              & diff_fun( l2gp1%solarTime, l2gp2%solarTime, &
              &         auxvalue=24._rgp, &
              &         options=AddOntoOptions('p',options) ), &
              & 'solarTime', &
              & options=options )
          endif
          myNumDiffs = myNumDiffs + count( l2gp1%solarTime /= l2gp2%solarTime )
        endif
      endif
      if ( all(shape(l2gp1%solarZenith) == shape(l2gp2%solarZenith)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%solarzenith', ' ', statsOnOneLine )
        CrashAtBeginning = .true.
        if ( any(l2gp1%solarZenith /= l2gp2%solarZenith) .and. &
          & SwitchDetail(lowercase(myFields), 'solarzenith', '-fc') > -1 ) then
            call diff ( l2gp1%solarZenith, 'solarZenith', &
              &         l2gp2%solarZenith, ' ', &
              & options=trim(options) )
          badChunks = .true.
          myNumDiffs = myNumDiffs + count( l2gp1%solarZenith /= l2gp2%solarZenith )
        endif
      endif
      if ( all(shape(l2gp1%losAngle) == shape(l2gp2%losAngle)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%losangle', ' ', statsOnOneLine )
        if ( any(l2gp1%losAngle /= l2gp2%losAngle) .and. &
          & SwitchDetail(lowercase(myFields), 'losangle', '-fc') > -1 ) then
            call diff ( l2gp1%losAngle, 'losAngle', &
              &         l2gp2%losAngle, ' ', &
              & options=options )
          myNumDiffs = myNumDiffs + count( l2gp1%losAngle /= l2gp2%losAngle )
        endif
      endif
      if ( all(shape(l2gp1%geodAngle) == shape(l2gp2%geodAngle)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%geodangle', ' ', statsOnOneLine )
        if ( any(l2gp1%geodAngle /= l2gp2%geodAngle) .and. &
          & SwitchDetail(lowercase(myFields), 'geodangle', '-fc') > -1 ) then
            call diff ( l2gp1%geodAngle, 'geodAngle', &
              &         l2gp2%geodAngle, ' ', &
              & options=options )
          ! badChunks = .true.
          myNumDiffs = myNumDiffs + count( l2gp1%geodAngle /= l2gp2%geodAngle )
        endif
      endif
      if ( all(shape(l2gp1%time) == shape(l2gp2%time)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%time', ' ', statsOnOneLine )
        if ( any(l2gp1%time /= l2gp2%time) .and. &
          & SwitchDetail(lowercase(myFields), 'time', '-fc') > -1 ) then
            call diff ( l2gp1%time, ' time', &
              &         l2gp2%time, ' ', &
              & options=options )
          ! badChunks = .true.
          myNumDiffs = myNumDiffs + count( l2gp1%time /= l2gp2%time )
        endif
      endif
      if ( all(shape(l2gp1%chunkNumber) == shape(l2gp2%chunkNumber)) ) then
        NameOnEachLine = merge( trim(trim(l2gp1%name)) // '%chunknumber', ' ', statsOnOneLine )
        if ( any(l2gp1%chunkNumber /= l2gp2%chunkNumber) .and. &
          & SwitchDetail(lowercase(myFields), 'chunknumber', '-fc') > -1 ) then
          call diff ( l2gp1%chunkNumber, 'chunkNumber', &
            &         l2gp2%chunkNumber, ' ', &
              & options=options )
          myNumDiffs = myNumDiffs + count( l2gp1%chunkNumber /= l2gp2%chunkNumber )
        endif
      endif
      
      if ( associated(l2gp1%frequency) .and.  associated(l2gp2%frequency)) then
        if ( all(shape(l2gp1%frequency) == shape(l2gp2%frequency)) ) then
          if ( any(l2gp1%frequency /= l2gp2%frequency) .and. &
            & SwitchDetail(lowercase(myFields), 'freq', '-fc') > -1 ) then
            call diff ( l2gp1%frequency, 'frequency', &
              &         l2gp2%frequency, ' ', &
              & options=options )
          myNumDiffs = myNumDiffs + count( l2gp1%frequency /= l2gp2%frequency )
          endif
        endif
      endif
    end subroutine diffGeoLocations

    subroutine doneHere
      ! Housekeeping
      ! Done with the temporary, so deallocate it before returning
      call DestroyL2GPContents ( L2GP2Temp )
      
      call resumeOutput
      if ( present(numDiffs) ) numDiffs = myNumDiffs
      nameOnEachLine = ' '
    end subroutine doneHere
    
    function AddOntoOptions( addOn, options) result(newOptions)
      character(len=*), intent(in)           :: AddOn
      character(len=*), optional, intent(in) :: options
      character(len=16)                      :: newOptions
      newOptions = AddOn
      if ( .not. present(options) ) return
      if ( len_trim(options) < 1 ) return
      newOptions = trim(options) // addOn
    end function AddOntoOptions
  end subroutine DiffThis
    
  ! ------------------------------------------ DiffL2GPFiles_MLSFile ------------
  subroutine DiffL2GPFiles_MLSFile ( L2GPFile1, L2GPFile2, &
    & geoBoxNames, geoBoxLowBound, geoBoxHiBound, pressures, chunks, &
    & Details, options, &
    & swList, showMissing, fields, force, swaths1, swaths2, &
    & matchTimes, numDiffs )
    ! Show diff between swaths in file1 and file2 down to level of Details
    ! Dummy arguments
    type(MLSFile_T)               :: L2GPfile1 ! file 1
    type(MLSFile_T)               :: L2GPfile2 ! file 2
    real, intent(in), dimension(:), optional :: pressures ! Which heights to dump
    integer, intent(in), dimension(:), optional :: Chunks ! Which chunks to dump
    character(len=*), intent(in) , optional ::       geoBoxNames
    real(rgp), dimension(:), intent(in), optional :: geoBoxLowBound  ! range
    real(rgp), dimension(:), intent(in), optional :: geoBoxHiBound  ! range
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    !
    ! The following parameters, if present, will override Details
    character(len=*), intent(in), optional :: options
    ! swList currently used only if showMissing is TRUE
    character (len=*), optional, intent(in) :: swList
    logical, intent(in), optional :: showMissing   ! if TRUE, just show which
    character(len=*), intent(in), optional :: fields  ! only these fields
    character(len=*), intent(in), optional :: swaths1  ! only these swaths
    character(len=*), intent(in), optional :: swaths2  ! only these swaths
    logical, intent(in), optional :: FORCE ! Force diff even if swathnames differ
    logical, intent(in), optional :: matchTimes  ! only matching profile times
    integer, intent(out), optional :: numDiffs  ! how many diffs
    ! Local                                         swaths are missing from other
    logical :: L1alreadyOpen
    logical :: L2alreadyOpen
    integer :: status
    !
    L1alreadyOpen = L2GPFile1%stillOpen
    if ( L1alreadyOpen ) then
      call mls_closeFile(L2GPFile1, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to close l2gp file to diff between', MLSFile=L2GPFile1)
    endif
    !
    L2alreadyOpen = L2GPFile2%stillOpen
    if ( L2alreadyOpen ) then
      call mls_closeFile(L2GPFile2, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to close l2gp file to diff between', MLSFile=L2GPFile2)
    endif
    
    if ( present(geoBoxNames) ) then
      call DiffL2GPFiles_Name ( L2GPFile1%Name, L2GPFile2%Name, &
      & geoBoxNames=geoBoxNames, &
      & geoBoxLowBound=geoBoxLowBound, geoBoxHiBound=geoBoxHiBound, &
      & Details=Details, options=options, &
      & swList=swList, showMissing=showMissing, &
      & fields=fields, force=force, swaths1=swaths1, swaths2=swaths2, &
      & matchTimes=matchTimes, &
      & numDiffs=numDiffs )
    elseif ( present(chunks) .and. present(pressures) ) then
      call DiffL2GPFiles_Name ( L2GPFile1%Name, L2GPFile2%Name, &
      & pressures=pressures, &
      & chunks=chunks, Details=Details, options=options, &
      & swList=swList, showMissing=showMissing, &
      & fields=fields, force=force, swaths1=swaths1, swaths2=swaths2, &
      & matchTimes=matchTimes, &
      & numDiffs=numDiffs )
    elseif ( present(chunks) ) then
      call DiffL2GPFiles_Name ( L2GPFile1%Name, L2GPFile2%Name, &
      & chunks=chunks, Details=Details, options=options, &
      & swList=swList, showMissing=showMissing, &
      & fields=fields, force=force, swaths1=swaths1, swaths2=swaths2, &
      & matchTimes=matchTimes, &
      & numDiffs=numDiffs )
    elseif ( present(pressures) ) then
      call DiffL2GPFiles_Name ( L2GPFile1%Name, L2GPFile2%Name, &
      & pressures=pressures, Details=Details, options=options, &
      & swList=swList, showMissing=showMissing, &
      & fields=fields, force=force, swaths1=swaths1, swaths2=swaths2, &
      & matchTimes=matchTimes, &
      & numDiffs=numDiffs )
    else
      call DiffL2GPFiles_Name ( L2GPFile1%Name, L2GPFile2%Name, &
      & Details=Details, options=options, &
      & swList=swList, showMissing=showMissing, &
      & fields=fields, force=force, swaths1=swaths1, swaths2=swaths2, &
      & matchTimes=matchTimes, &
      & numDiffs=numDiffs )
    endif

    if ( L1alreadyOpen )  call mls_openFile(L2GPFile1, Status)
    L2GPFile1%errorCode = status
    L2GPFile1%lastOperation = 'read'

    if ( L2alreadyOpen )  call mls_openFile(L2GPFile2, Status)
    L2GPFile2%errorCode = status
    L2GPFile2%lastOperation = 'read'

  end subroutine DiffL2GPFiles_MLSFile
    
  ! ------------------------------------------ DiffL2GPFiles_Name ------------
  subroutine DiffL2GPFiles_Name ( file1, file2, &
    & geoBoxNames, geoBoxLowBound, geoBoxHiBound, pressures, chunks, &
    & Details, options, &
    & swList, showMissing, fields, force, swaths1, swaths2, matchTimes, &
    &  numDiffs )
    ! Show diff between swaths in file1 and file2 down to level of Details
    ! Dummy arguments
    character (len=*), intent(in) :: file1 ! Name of file 1
    character (len=*), intent(in) :: file2 ! Name of file 2
    real, intent(in), dimension(:), optional :: pressures ! Which heights to diff
    integer, intent(in), dimension(:), optional :: Chunks ! Which chunks to diff
    character(len=*), intent(in) , optional ::       geoBoxNames
    real(rgp), dimension(:), intent(in), optional :: geoBoxLowBound  ! range
    real(rgp), dimension(:), intent(in), optional :: geoBoxHiBound  ! range
    integer, intent(in), optional :: DETAILS ! <=0 => Don't diff data fields
    !                                        ! -1 Skip even geolocation fields
    !                                        ! -2 Skip all but name
    !                                        ! >0 Diff even data fields
    !                                        ! Default 1
    !
    ! The following parameters, if present, will override Details
    character(len=*), intent(in), optional :: options
    ! swList currently used only if showMissing is TRUE
    character (len=*), optional, intent(in) :: swList
    logical, intent(in), optional :: showMissing   ! if TRUE, just show which
    character(len=*), intent(in), optional :: fields  ! only these fields
    character(len=*), intent(in), optional :: swaths1  ! only these swaths
    character(len=*), intent(in), optional :: swaths2  ! only these swaths
    logical, intent(in), optional :: FORCE ! Force diff even if swathnames differ
    logical, intent(in), optional :: matchTimes  ! only matching profile times
    integer, intent(out), optional :: numDiffs  ! how many diffs
    ! Local                                         swaths are missing from other
    logical, parameter            :: countEmpty = .true.
    logical :: file_exists
    integer :: File1Handle
    integer :: File2Handle
    integer :: i
    type (L2GPData_T) :: l2gp1
    type (L2GPData_T) :: l2gp2
    integer :: listsize
    type(MLSFile_T)                :: MLSFile1, MLSFile2
    logical :: myForce
    integer :: myNumDiffs
    logical :: myShowMissing
    logical :: mySilent
    integer :: noSwaths
    integer :: noSwaths2
    integer :: noUnique
    integer :: status
    character (len=L2GPNameLen) :: swath
    character (len=L2GPNameLen) :: swath2
    character (len=MAXSWATHNAMESBUFSIZE) :: swathList1
    character (len=MAXSWATHNAMESBUFSIZE) :: swathList2
    character (len=MAXSWATHNAMESBUFSIZE) :: swathUnique
    integer :: the_hdfVersion1
    integer :: the_hdfVersion2
    ! Executable code
    myNumDiffs = 0
    myShowMissing = .false.
    if ( present(showMissing) ) myShowMissing=showMissing
    file_exists = ( mls_exists(trim(File1)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 1 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    myForce = .false.
    if ( present(force) ) myForce = force
    mySilent = .false.
    if ( present(options) ) mySilent = index(options, 'm') > 0

    the_hdfVersion1 = mls_hdf_version(File1)
    the_hdfVersion2 = mls_hdf_version(File2)
    file_exists = ( mls_exists(trim(File2)) == 0 )
    if ( .not. file_exists ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File 2 not found; make sure the name and path are correct' &
        & // trim(file1) )
    endif
    the_hdfVersion2 = mls_hdf_version(File2)
    noSwaths = mls_InqSwath ( file1, swathList1, listSize, &
         & hdfVersion=the_hdfVersion1)
    noSwaths2 = mls_InqSwath ( file2, swathList2, listSize, &
         & hdfVersion=the_hdfVersion2)
    if ( present(swaths1) ) then
      if ( swaths1 /= 'all' .and. swaths1 /= '*' ) then
        swathList1 = swaths1
        noSwaths = NumStringElements( swaths1, countEmpty=.true. )
        swathList2 = swathList1
      endif
    endif
    if ( present(swaths2) ) then
      if ( swaths2 /= 'all' .and. swaths2 /= '*' ) then
        swathList2 = swaths2
        noSwaths2 = NumStringElements( swaths2, countEmpty=.true. )
      endif
    endif
    ! Are we merely to point out which swaths are missing from the other file?
    if ( myShowMissing ) then
      if ( present(swList) ) then
        call GetUniqueList( swList, swathUnique, noUnique, &
          & str2=trim(swathList1), options='-e' )
        if ( noUnique > 0 ) then
          call output('swaths missing from ' // trim(File1), advance='yes')
          call output(trim(swathUnique), advance='yes')
        else
          call output('All swaths in ' // trim(File1), advance='yes')
        endif
        call GetUniqueList( swList, swathUnique, noUnique, &
          & str2=trim(swathList2), options='-e' )
        if ( noUnique > 0 ) then
          call output('swaths missing from ' // trim(File2), advance='yes')
          call output(trim(swathUnique), advance='yes')
        else
          call output('All swaths in ' // trim(File2), advance='yes')
        endif
      else
        call output('Comparing swaths in ' // trim(File1), advance='no')
        call output(' with ' // trim(File2), advance='yes')
        call GetUniqueList( trim(swathList1), swathUnique, noUnique, &
          & str2=trim(swathList2), options='-e' )
        if ( noUnique > 0 ) then
          call output('swaths only in ' // trim(File1), advance='yes')
          call output(trim(swathUnique), advance='yes')
        endif
        call GetUniqueList(trim(swathList2), swathUnique, noUnique, &
          & str2=trim(swathList1), options='-e' )
        if ( noUnique > 0 ) then
          call output('swaths only in ' // trim(File2), advance='yes')
          call output(trim(swathUnique), advance='yes')
        endif
      endif
      return
    endif
    status = InitializeMLSFile ( MLSFile1, type=l_swath, access=DFACC_READ, &
     & name=trim(file1), HDFVersion=the_hdfVersion1 )
    call MLS_OpenFile( MLSFile1 )
    File1Handle = MLSFile1%FileID%f_id
    status = InitializeMLSFile ( MLSFile2, type=l_swath, access=DFACC_READ, &
     & name=trim(file2), HDFVersion=the_hdfVersion2 )
    call MLS_OpenFile( MLSFile2 )
    File2Handle = MLSFile2%FileID%f_id

    ! Loop over swaths in file 1

    do i = 1, noSwaths
      call GetStringElement (trim(swathList1), swath, i, countEmpty )
      if ( len_trim(swath) < 1 ) then
        if ( .not. mySilent ) call output('(Ignoring blank swath name in ' // &
          &  trim(File1), advance='yes')
        cycle
      endif
      if ( swathList2 /= swathList1 ) then
        call GetStringElement (trim(swathList2), swath2, i, countEmpty )
      else
        swath2 = swath
      endif
      status = stringElementNum(swathList2, trim(swath2), countEmpty)
      if ( status < 1 ) then
        if ( .not. myForce ) then
          if ( .not. mySilent ) call output('Swath ' // trim(swath) // ' not found in ' // &
            & trim(File2), advance='yes')
          cycle
        elseif ( .not. mySilent) then
          call GetStringElement (trim(swathList2), swath2, i, countEmpty )
          call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
          call output( '---- swath(1) name: ' // trim(swath) // ' ----', advance='yes')
          call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
          call output( '     swath(2) name: ' // trim(swath2)// '     ', advance='yes')
          call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
        endif
      elseif( .not. mySilent ) then
        call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
        call output( '---- swath name: ' // trim(swath) // ' ----', advance='yes')
        call blanks( 22+len_trim(swath), fillChar='-', advance='yes')
      endif
      call ReadL2GPData ( File1Handle, trim(swath), l2gp1, &
           & hdfVersion=the_hdfVersion1 )
      call ReadL2GPData ( File2Handle, trim(swath2), l2gp2, &
           & hdfVersion=the_hdfVersion2 )
      if ( present(geoBoxNames) ) then
        call DiffL2GPData_RANGES( l2gp1, l2gp2, &
        & geoBoxNames, geoBoxLowBound, geoBoxHiBound, chunks=chunks, &
        & details=details, options=options, fields=fields, &
        & numDiffs=numDiffs )
      elseif ( present(chunks) .and. present(pressures) ) then
        call Diff( l2gp1, l2gp2, pressures=pressures, chunks=chunks, &
        & details=details, options=options, fields=fields, &
        & numDiffs=numDiffs )
      elseif ( present(chunks) ) then
        call Diff( l2gp1, l2gp2, chunks=chunks, &
        & details=details, options=options, fields=fields, &
        & numDiffs=numDiffs )
      elseif ( present(pressures) ) then
        call Diff( l2gp1, l2gp2, pressures=pressures, &
        & details=details, options=options, fields=fields, &
        & numDiffs=numDiffs )
      else
        call Diff( l2gp1, l2gp2, &
        & details=details, options=options, fields=fields, &
        & numDiffs=numDiffs )
      endif

      if ( present(numDiffs) ) myNumDiffs = myNumDiffs + numDiffs
      call DestroyL2GPContents ( l2gp1 )
      call DestroyL2GPContents ( l2gp2 )
    enddo
    call MLS_CloseFile( MLSFile1 )
    call MLS_CloseFile( MLSFile2 )
    if ( present(numDiffs) ) numDiffs = myNumDiffs
  end subroutine DiffL2GPFiles_Name
    
  ! ------------------------------------------ DiffStatsInt ------------
  subroutine DiffStatsInt ( array1, array2, arrayName )
    ! Show diff between arrays 1 and 2: num diff num same, pctages
    integer, dimension(:), intent(in) :: array1
    integer, dimension(:), intent(in) :: array2
    character(len=*), intent(in) :: arrayName
    ! Internal variables
    integer :: numSame
    integer :: numDiff
    real :: pctSame
    real :: pctDiff
    ! Executable
    call output('Differences between (1) and (2) integer arrays ' &
      & // trim(arrayName), advance='yes')
    if ( size(array1) /= size(array2) ) then
      call output('(arrays (1) and (2) different lengths)', advance='yes')
      return
    endif
    numDiff = count(array1 /= array2)
    numSame = size(array1) - numDiff
    pctDiff = 100*numDiff/(0.+numDiff+numSame)
    pctSame = 100. - pctDiff
    call output(numSame, advance = 'no')
    call blanks(2, advance='no')
    call output('same', advance = 'no')
    call blanks(4, advance='no')
    call output(numDiff, advance = 'no')
    call blanks(2, advance='no')
    call output('different', advance = 'no')
    call blanks(4, advance='no')
    call output(pctSame, advance = 'no')
    call output('%  same', advance = 'no')
    call blanks(4, advance='no')
    call output(pctDiff, advance = 'no')
    call output('%  different', advance = 'yes')
  end subroutine DiffStatsInt

  ! ------------------------------------------ DUMP_L2GP ------------
  ! This family of routines Dumps l2gp datatype(s) according to
  ! their number (singly or an array of them)
  ! whole or subsetted
  ! the data or its metadata (i.e., attributes)
  ! ------------------------------------------ DUMP_L2GP_DATABASE ------------
  subroutine DUMP_L2GP_DATABASE ( L2gp, &
    & Name, ColumnsOnly, Details, Fields, options )

    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP(:)
    character(len=*), intent(in), optional :: Name
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS
    character(len=*), intent(in), optional :: fields ! ,-separated list of names
    character(len=*), intent(in), optional :: options ! Passed dumping arrays

    ! Local variables
    integer :: i
    call output ( '============ L2GP Data Base ============', advance='yes' )
    call output ( ' ', advance='yes' )
    if ( present(name) ) then
      call StyledOutput ( 'L2GP Database name: ' // trim(Name), options )
    endif
    if ( size(l2gp) < 1 ) then
      call output ( '**** L2GP Database empty ****', advance='yes' )
      return
    endif
    do i = 1, size(l2gp)
      call dump( l2gp(i), ColumnsOnly, Details, Fields, options=options )
    end do
      
  end subroutine DUMP_L2GP_DATABASE

  ! ------------------------------------------ DUMP_L2GP_CHUNKS ------------
  subroutine DUMP_L2GP_CHUNKS ( fullL2gp, Chunks, &
    & ColumnsOnly, Details, Fields, Width, options )
    ! Dump selected chunks of an l2gp
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          FULLL2GP
    integer, intent(in), dimension(:) :: Chunks ! Which chunks to dump
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: fields ! ,-separated list of names
    integer, intent(in), optional :: width   ! width of each dumped line
    character(len=*), intent(in), optional :: options ! Passed dumping arrays
    ! Local variables
    integer :: chunk
    integer :: i
    integer, dimension(2) :: irange
    type (l2gpData_T) ::          L2GP
    ! Executable
    do i=1, size(Chunks)
      chunk = Chunks(i)
      irange(1) = FindFirst( fullL2gp%chunkNumber, chunk )
      irange(2) = Findlast( fullL2gp%chunkNumber, chunk )
      ! HuntRange does not work correctly when any chunk Numbers are -999
      ! (Missing Value)
      ! call HuntRange( fullL2gp%chunkNumber, (/ chunk, chunk /), irange )
      if ( any( irange == 0 ) ) cycle
      call output ( ' - - - Chunk number:', advance='no')
      call output ( chunk, advance='no')
      call output ( ' - - -', advance='yes')
      call ExtractL2GPRecord ( fullL2gp, l2gp, rTimes=irange )
      call Dump ( L2gp, ColumnsOnly, Details, Fields, width, options )
      call DestroyL2GPContents( l2gp )
    enddo
  end subroutine DUMP_L2GP_CHUNKS

  ! ------------------------------------------ DUMPL2GPData_RANGES ------------
  subroutine DUMPL2GPData_RANGES ( fullL2gp, &
    & geoBoxNames, geoBoxLowBound, geoBoxHiBound, &
    & Pressures, Latitudes, Longitudes, Times, Chunks, ColumnsOnly, &
    & Details, fields, Width, options )
    ! DUMP an l2gp according to prescribed ranges of pressure, longitude, etc.
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          FULLL2GP
    character(len=*), intent(in) , optional ::       geoBoxNames
    real(rgp), dimension(:), intent(in), optional :: geoBoxLowBound  ! range
    real(rgp), dimension(:), intent(in), optional :: geoBoxHiBound  ! range
    real(rgp), intent(in), dimension(2), optional :: pressures ! range
    real(rgp), intent(in), dimension(2), optional :: latitudes ! range
    real(rgp), intent(in), dimension(2), optional :: longitudes ! range
    real(r8), intent(in), dimension(2), optional :: times ! range
    integer, intent(in), optional, dimension(:) :: Chunks ! Which chunks to dump
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: fields ! ,-separated list of names
    integer, intent(in), optional :: width   ! width of each dumped line
    character(len=*), intent(in), optional :: options ! Passed dumping arrays
    ! Local variables
    ! logical, parameter :: DEEBUG = .true.
    type (l2gpData_T) ::          L2GP
    ! Executable
    if ( DEEBUG ) then
      if ( present(pressures) ) &
        & call dump( pressures, 'pressures' )
      if ( present(latitudes) ) &
        & call dump( latitudes, 'latitudes' )
      if ( present(longitudes) ) &
        & call dump( longitudes, 'longitudes' )
      call output( 'fullL2gp%pressures: ', advance='no' )
      call output( fullL2gp%pressures, advance='yes' )
    endif
    if ( present(geoBoxNames) ) then
      call ContractL2GPRecord ( fullL2gp, l2gp, &
        & geoBoxNames, geoBoxLowBound, geoBoxHiBound )
    else
      call ContractL2GPRecord ( fullL2gp, l2gp, pressures=pressures, &
        & latitudes=latitudes, longitudes=longitudes, intimes=times )
    endif
    if ( present(chunks) ) then
      call DUMP_L2GP_Chunks ( L2gp, Chunks, &
        & ColumnsOnly, Details, Fields, width, options )
    else
      call DUMP ( L2gp, ColumnsOnly, Details, Fields, width, options )
    endif
    call DestroyL2GPContents( l2gp )
  end subroutine DUMPL2GPData_RANGES

  ! ------------------------------------------ DUMP_L2GP ------------
  subroutine Dump_L2GP ( L2gp, ColumnsOnly, Details, Fields, Width, options )
    ! Dump an l2gp
    ! Either according to level of detail set by Details
    ! or else just those fields named in Fields
    ! Note: 'solarTime' among Fields will also provoke dumping 'time'
    ! Dummy arguments
    type (l2gpData_T), intent(in) ::          L2GP
    logical, intent(in), optional ::          ColumnsOnly ! if true, dump only with columns
    integer, intent(in), optional :: DETAILS ! <=0 => Don't dump multidim arrays
    !                                        ! -1 Skip even 1-d arrays
    !                                        ! -2 Skip all but name
    !                                        ! >0 Dump even multi-dim arrays
    !                                        ! Default 1
    character(len=*), intent(in), optional :: fields ! ,-separated list of names
    integer, intent(in), optional :: width   ! width of each dumped line
    character(len=*), intent(in), optional :: options ! Passed dumping arrays

    ! Local variables
    real(r8) :: FillValue
    real(rgp) :: FillValueGP
    real(rgp), dimension(:), pointer :: hoursInDay
    integer :: i
    integer :: ierr
    logical :: myColumnsOnly
    integer :: MYDETAILS
    character(len=Len(DATA_FIELDS)+Len(GEO_FIELDS)) :: myFields
    integer :: nUnique
    logical :: skipGeos
    real(rgp), dimension(:), pointer :: spacing
    integer, dimension(1000) :: uniqueVals

    ! Executable code
    myDetails = 1
    if ( present(details) ) myDetails = details
    myFields = ' '
    if ( present(fields) ) then
      myFields = fields
      if ( .not. present(details) ) myDetails = 2
    endif
    
    if( present(ColumnsOnly)) then
      myColumnsOnly = ColumnsOnly
    else
      myColumnsOnly = .false.
    endif

    if ( myColumnsOnly .and. l2gp%nLevels > 1 ) return
    skipGeos = .false.
    if ( present(options) ) skipGeos = ( index(options, 'g') > 0 )
    FillValue = real(l2gp%MissingValue, r8)
    FillValueGP = l2gp%MissingValue
    if ( showMe(.true., myFields, 'swathname') ) then
      call StyledOutput ( 'L2GP Data: (swath name) ' // trim(l2gp%name), &
        & options='--Banner' )
      if ( NAMEINDEXEVERSET ) then
        call output ( ', (parser name) ')
        if(l2gp%nameIndex > 0) then
          call display_string ( l2gp%nameIndex, advance='yes', IERR=ierr)
          if ( ierr /= 0 ) call output ( '(not found in string table)', &
           & advance='yes')
        else
          call output ( '(the nameIndex was 0) ', advance='yes')
        endif
      else
        call output ( ' ', advance='yes')
      endif
    endif

    if ( showMe(myDetails > -2, myFields, 'ntimes') ) then
      call output ( 'nTimes: ')
      call output ( l2gp%nTimes, 5)
      call output ( '  nTimesTotal: ')
      call output ( l2gp%nTimesTotal, 5)
      call output ( '  nLevels: ')
      call output ( l2gp%nLevels, 3)
      call output ( '  nFreqs: ')
      call output ( l2gp%nFreqs, 3, advance='yes')
      call output ( 'Fill/Missing L2GP Values: ')
      call output ( l2gp%MissingL2GP, advance='yes')
      call output ( 'Fill/Missing Values (Others): ')
      call output ( l2gp%MissingValue, advance='yes')
      call output ( 'Fill/Missing Status Field: ')
      call output ( l2gp%MissingStatus, advance='yes')
     endif
    
    if ( .not. skipGeos ) then
      if ( showMe(myDetails > -1, myFields, 'pressure') ) &
        & call dump ( l2gp%pressures, trim(l2gp%verticalCoordinate) // 's:', &
        & FillValue=FillValueGP, width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'latitude') ) &
        & call dump ( l2gp%latitude, 'Latitude:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'longitude') ) &
        & call dump ( l2gp%longitude, 'Longitude:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'solartime') ) &
        & call dump ( l2gp%solarTime, 'SolarTime:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'solarzenith') ) &
        & call dump ( l2gp%solarZenith, 'SolarZenith:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'LOSAngle') ) &
        & call dump ( l2gp%losAngle, 'LOSAngle:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'geodAngle') ) &
        & call dump ( l2gp%geodAngle, 'geodAngle:', FillValue=FillValueGP, &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'time') ) then
        nullify( hoursInDay )
        call dump ( l2gp%time, 'Time:', FillValue=FillValue, &
          & width=width, options=options )
        call allocate_test ( hoursInDay, l2gp%nTimes, 'hoursInDay', ModuleName )
        where ( l2gp%time /= FillValue )
          hoursInDay = timeToHoursInDay ( l2gp%time )
        elsewhere
          hoursInDay = FillValueGP
        end where
        call dump ( hoursInDay, 'Hours In Day:', FillValue=FillValueGP, &
          & width=width, options=options )
        call deallocate_test ( hoursInDay, 'hoursInDay', ModuleName )
      endif

      if ( showMe(myDetails > -1, myFields, 'spacing') ) then
        nullify( spacing )
        call allocate_test ( spacing, l2gp%nTimes-1, 'spacing', ModuleName )
        spacing = 0.
        do i=1, l2gp%nTimes - 1
          spacing(i) = EarthRadA * sqrt ( &
            & (sin(Deg2Rad*l2gp%latitude(i+1)) - sin(Deg2Rad*l2gp%latitude(i)))**2 &
            & + &
            & (cos(Deg2Rad*l2gp%latitude(i+1))*cos(Deg2Rad*l2gp%longitude(i+1)) &
            & - cos(Deg2Rad*l2gp%latitude(i))*cos(Deg2Rad*l2gp%longitude(i)))**2 &
            & + &
            & (cos(Deg2Rad*l2gp%latitude(i+1))*sin(Deg2Rad*l2gp%longitude(i+1)) &
            & - cos(Deg2Rad*l2gp%latitude(i))*sin(Deg2Rad*l2gp%longitude(i)))**2 &
            & )
        enddo
        call dump ( spacing, 'spacing (in m):', FillValue=FillValueGP, &
          & width=width, options=options )
        call deallocate_test ( spacing, &
          & 'spacing', ModuleName )
      endif

      if ( showMe(myDetails > -1, myFields, 'chunkNumber') ) &
        & call dump ( l2gp%chunkNumber, 'ChunkNumber:', &
        & width=width, options=options )

      if ( showMe(myDetails > -1, myFields, 'pressure') .and. &
        & associated(l2gp%frequency) ) &
        & call dump ( l2gp%frequency, 'Frequencies:', FillValue=FillValueGP, &
        & width=width, options=options )

    endif      

    if ( showMe(myDetails > 0, myFields, 'l2gpvalue') ) then
      if ( l2gp%nFreqs < 2 ) then
        call dump ( l2gp%l2gpValue(1,:,:), 'L2GPValue:', &
          & FillValue=l2gp%MissingL2GP, options=options )
      else
        call dump ( l2gp%l2gpValue, 'L2GPValue:', &
          & FillValue=l2gp%MissingL2GP, options=options )
      endif
    endif
      
    if ( showMe(myDetails > 0, myFields, 'l2gpprecision') ) then
      if ( l2gp%nFreqs < 2 ) then
        call dump ( real(l2gp%l2gpprecision(1,:,:), r8), 'L2GPPrecision:', &
          & FillValue=FillValue, options=options )
      else
        call dump ( real(l2gp%l2gpprecision, r8), 'L2GPprecision:', &
          & FillValue=FillValue, options=options )
      endif
    endif
      
    if ( showMe(myDetails > 0, myFields, 'status') ) then
      call dump ( l2gp%status, 'Status:', options=options )
      call FindUnique( l2gp%status, uniqueVals, nUnique )
      call outputNamedValue(' Number unique values', nUnique )
      do i=1, nUnique
        call DumpBitNames( uniqueVals(i), StatusBitNames )
      enddo
      
    endif      
    if ( showMe(myDetails > 0, myFields, 'quality') ) &
      & call dump ( l2gp%quality, 'Quality:', &
        & FillValue=FillValueGP, options=options )
      
    if ( showMe(myDetails > 0, myFields, 'convergence') ) &
      & call dump ( l2gp%convergence, 'Convergence:', &
        & FillValue=FillValueGP, options=options )
      
    if ( showMe(myDetails > 0, myFields, 'AscDescMode ') .and. AscDescModeIsField ) &
      & call dump ( l2gp%AscDescMode, 'AscDescMode :', &
        & width=width, options=options )
      
    if ( showMe(myDetails > 0, myFields, 'BinNumber') .and. associated(l2gp%BinNumber) ) &
      & call dump ( l2gp%BinNumber, 'BinNumber:', &
        & options=options )
      
    if ( showMe(myDetails > 0, myFields, 'MAF') .and. associated(l2gp%MAF) ) &
      & call dump ( l2gp%MAF, 'MAF:', &
        & options=options )
      
  contains
    logical function showMe(detailsOK, fields, field)
      ! Determine whether this field should be dumped or not
      ! depending on whether you're going by details or by picking fields
      logical, intent(in) :: detailsOK
      character(len=*), intent(in) :: fields
      character(len=*), intent(in) :: field
      !
      if ( len_trim(fields) < 1 ) then
        showMe = detailsOK
      else
        showMe = ( &
          & switchDetail(LowerCase(fields), LowerCase(trim(field)), '-fc' ) &
          & > -1 )
      endif
    end function showMe
  end subroutine Dump_L2GP
    
  !----------------------------------------  DumpL2GP_attributes_hdf5  -----
  subroutine DumpL2GP_attributes_hdf5( l2FileHandle, l2gp, swathName )
  use HE5_SWAPI, only: HE5_SWRdattr, HE5_SWRdlattr
  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach, HE5_EHRdglatt
  use PCFHdr, only:  GlobalAttributes_T, HE5_ReadGlobalAttr, &
    & DumpGlobalAttributes
    ! Brief description of subroutine
    ! This subroutine dumps the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    integer, intent(in) :: l2FileHandle ! From swopen
    type( L2GPData_T ), intent(inout) :: l2gp
    character (len=*), optional, intent(in) :: swathName ! Defaults->l2gp%name

    ! Variables
    type (GlobalAttributes_T)         :: gAttributes
    character(len=255)                :: ProcessLevel
    integer                           :: DayofYear
    double precision                  :: TAI93At0zOfGranule
    real(rgp), dimension(MAXNLEVELS)  :: pressures

    character (len=132) :: name     ! Either swathName or l2gp%name
    ! The following pair of string list encode the Units attribute
    ! corresponding to each Title attribute; e.g., the Units for Latitude is deg
    character (len=*), parameter :: GeolocationTitles = &
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,' // &
      & 'LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,' // &
      & 'deg,deg,NoUnits,hPa,GHz'
    character (len=*), parameter :: GeoUniqueFieldDefinition = &
      & 'HMT,HMT,AS,HMT,HMT,' // &
      & 'M,M,M,AS,M'   ! These are abbreviated values
    character (len=*), parameter :: UniqueFieldDefKeys = &
      & 'HM,HMT,MT,AS,M'
    character (len=*), parameter :: UniqueFieldDefValues = &
      & 'HIRDLS-MLS-Shared,HIRDLS-MLS-TES-Shared,MLS-TES-Shared,' // &
      & 'Aura-Shared,MLS-Specific'  ! Expanded values
    ! The following associate UniqueFieldDefs with species names
    character (len=*), parameter :: Species = &
      & 'Temperature,BrO,CH3CN,CO,ClO,GPH,HCl,HCN,H2O,H2O2,' // &
      & 'HNO3,HOCl,HO2,N2,N2O,OH,O2,O3,RHI,SO2'
    character (len=*), parameter :: SpUniqueFieldDefinition = &
      & 'HMT,M,M,MT,M,M,M,M,HMT,M,' // &
      & 'HMT,M,M,M,HM,M,M,HMT,M,M'   ! These are abbreviated values

    integer :: field
    logical :: isColumnAmt
    integer :: status
    integer :: swid
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    character(len=CHARATTRLEN) :: field_name
    character(len=CHARATTRLEN) :: units_name
    character(len=CHARATTRLEN) :: uniqueness
    character(len=CHARATTRLEN) :: species_name
    character(len=CHARATTRLEN) :: abbr_uniq_fdef
    character(len=CHARATTRLEN) :: expnd_uniq_fdef
    real(rgp), dimension(1) :: MissingValue
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    call outputNamedValue ( 'L2GP Attributes: (swath name) ', trim(name), &
      & options = '-H' )
    if ( .not. any( rgp /= (/ r4, r8 /) ) ) then
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8' )
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    call output ( '(Global Attributes) ', advance='yes' )
    call he5_readglobalattr( l2FileHandle, gAttributes, &
     & ProcessLevel, DayofYear, TAI93At0zOfGranule, returnStatus=status )
    if ( status /= 0 ) then
      call output ('No global attributes found in file', advance='yes')
    else
      ! We stopped reading MiscNotes in the PCFHdr module to avoid
      ! letting the master task step on the value each slave had carefully
      ! written there
      ! Therefore, we must read it here , now
      ! (Sigh)
      status = he5_EHrdglatt( l2FileHandle, &
      & 'MiscNotes', &
      &  gAttributes%MiscNotes )
      status = he5_EHrdglatt( l2FileHandle, &
      & 'ProductionLocation', &
      &  gAttributes%ProductionLoc )
      status = he5_EHrdglatt( l2FileHandle, &
      & 'HostName', &
      &  gAttributes%HostName )
      status = he5_EHrdglatt( l2FileHandle, &
      & 'identifier_product_doi', &
      &  gAttributes%DOI )
      call dumpGlobalAttributes ( gAttributes ) 
      call dump_int ( DayOfYear, 'Granule day of year:' )
      call dump_r8 ( TAI93At0zOfGranule, 'Equator crossing time (tai93):' )
    endif
    swid = mls_SWattach (l2FileHandle, name, hdfVersion=HDFVERSION_5)
    if ( swid == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to attach swath ' // trim(name))
    end if
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    call output ( '(Swath Attributes) ', advance='yes')
    status = he5_swrdattr(swid, 'Pressure', pressures)
    call dump ( pressures, 'Vertical coordinates:' )
    status = he5_swrdattr(swid, 'VerticalCoordinate', field_name)
    call dump_chars ( field_name, 'Vertical coordinates type:' )
    status = he5_swrdlattr( swid, 'L2gpValue', 'MissingValue', MissingValue )
    call dump ( MissingValue, 'MissingValues (L2GPValues only):' )
    status = he5_swrdlattr( swid, 'L2gpPrecision', 'MissingValue', MissingValue )
    call dump ( MissingValue, 'MissingValues (other fields):' )
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
    call output ( '(Geolocation Attributes) ', advance='yes')
    do field=1, NumGeolocFields
      ! Take care not to write attributes to "missing fields"
      if ( trim(theTitles(field)) == 'Frequency' &
        & .and. l2gp%nFreqs < 1 ) then
        field_name = ''
      elseif ( trim(theTitles(field)) == 'Pressure' &
        & .and. l2gp%nLevels < 1 ) then
        field_name = ''
      else
        call GetHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        status = he5_swrdlattr(swid, trim(theTitles(field)), 'Title', field_name)
        status = he5_swrdlattr(swid, trim(theTitles(field)), 'Units', units_name)
        call dump_chars ( field_name, 'Field title:' )
        call dump_chars ( units_name, 'Units:' )

        status = he5_swrdlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', uniqueness)
        call dump_chars ( uniqueness, 'Unique field definition:' )
      endif
    enddo
    !   - -   D a t a   A t t r i b u t e s   - -
    call output ( '(Data Attributes) ', advance='yes')
    field_name = Name
    species_name = name
    isColumnAmt = ( index(species_name, 'Column') > 0 )
    if ( isColumnAmt ) then
      call ExtractSubString(Name, species_name, 'Column', 'wmo')
    endif
    call GetHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(lowercase(species_name)), &
      & abbr_uniq_fdef, .false.)
    call GetHashElement (UniqueFieldDefKeys, &
      & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
      & expnd_uniq_fdef, .false.)
    select case (trim(lowercase(species_name)))
    case ('temperature')
      units_name = 'K'
    case ('gph')
      units_name = 'm'
    case ('rhi')
      units_name = '%rhi'
    case default
      units_name = 'vmr'
    end select
    if ( isColumnAmt ) units_name = 'DU'
    status = he5_swrdlattr(swid, 'L2gpValue', 'Title', field_name)
    status = he5_swrdlattr(swid, 'L2gpValue', 'Units',  units_name)
    status = he5_swrdlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    status = he5_swrdlattr(swid, 'L2gpPrecision', 'Title', field_name)
    status = he5_swrdlattr(swid, 'L2gpPrecision', 'Units', units_name)
    status = he5_swrdlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )

    ! ('Status' data field newly written)
    status = he5_swrdlattr(swid, 'Status', 'Title', field_name)
    status = he5_swrdlattr(swid, 'Status', 'Units', units_name)
    status = he5_swrdlattr(swid, 'Status', &
      & 'UniqueFieldDefinition', uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    
    status = he5_swrdlattr(swid, 'Quality', 'Title', field_name)
    status = he5_swrdlattr(swid, 'Quality', 'Units', units_name)
    status = he5_swrdlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', uniqueness)
    call dump_chars ( field_name, 'Field title:' )
    call dump_chars ( units_name, 'Units:' )
    call dump_chars ( uniqueness, 'Unique field definition:' )
    
    status = mls_SWdetach(swid, hdfVersion=HDFVERSION_5)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface' )
    end if
    contains
    subroutine dump_chars(value, name)
      ! arguments
      character(len=*), intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(trim(value), advance='yes')
    end subroutine dump_chars
    subroutine dump_int(value, name)
      ! arguments
      integer, intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(value, advance='yes')
    end subroutine dump_int
    subroutine dump_r8(value, name)
      ! arguments
      real(r8), intent(in) :: value
      character(len=*), intent(in) :: name
      ! Executable
      call output(trim(name) // ' ', advance='no')
      call output(value, advance='yes')
    end subroutine dump_r8

  !-------------------------------------
  end subroutine DumpL2GP_attributes_hdf5
  !-------------------------------------

  !------------------------------------------  ExtractL2GPRecord  -----
  ! Extract a reduced copy of an original l2gp record
  ! picking out a range of subscripts for freqs, level2, times
  ! Omitted ranges will be same range as original l2gp

  ! Note that this does not do the more useful task:
  ! Reposition an l2gp record onto a new set of geolocations
  ! via multi-dimensional or repeated interpolation

  ! Nor does it (yet) do another potentially useful task:
  ! Extract the range of profiles within a box of times, lats, lons
  ! and a range of levels within a range of pressures

  ! These last tasks aren't too difficult, since we could define
  ! myTimes and myLevels based on HuntRange.

  ! A tricky point: invalid ranges (i.e. upper limit < 1) will
  ! be treated as if they were omitted

  subroutine ExtractL2GPRecord_range ( ol2gp, l2gp, rFreqs, rLevels, rTimes )
    ! Dummy arguments
    type (L2GPData_T), intent(in)   :: ol2gp
    type (L2GPData_T), intent(out)  :: l2gp
    integer, dimension(2), intent(in), optional :: rFreqs  ! subscript range
    integer, dimension(2), intent(in), optional :: rLevels ! subscript range
    integer, dimension(2), intent(in), optional :: rTimes  ! subscript range
    ! Local variables
    integer :: i
    integer, dimension(2) :: myFreqs  ! subscript range
    integer, dimension(2) :: myLevels ! subscript range
    integer, dimension(2) :: myTimes  ! subscript range

    integer, dimension(3) :: start
    integer, dimension(3) :: count
    integer, dimension(3) :: block
    integer, dimension(3) :: stride
    logical :: useFreqs  
    logical :: useLevels 
    logical :: useTimes  
    
    ! Executable
    block = 1
    stride = 1
    useFreqs   = .false.
    useLevels  = .false.
    useTimes   = .false.
    if ( present(rFreqs) ) then
      useFreqs = rFreqs(2) > 0
    endif
    if ( present(rLevels) ) then
      useLevels = rLevels(2) > 0
    endif
    if ( present(rTimes) ) then
      useTimes = rTimes(2) > 0
    endif
    if ( useFreqs ) then
      myFreqs = rFreqs
    else
      myFreqs(1) = 1
      myFreqs(2) = ol2gp%nFreqs
    endif
    if ( useLevels ) then
      myLevels = rLevels
    else
      myLevels(1) = 1
      myLevels(2) = ol2gp%nLevels
    endif
    if ( useTimes ) then
      myTimes = rTimes
    else
      myTimes(1) = 1
      myTimes(2) = ol2gp%nTimes
    endif
    
    call SetupNewL2GPRecord ( l2gp, &
      & ( myFreqs(2) - myFreqs(1) + 1 ), &
      & ( myLevels(2) - myLevels(1) + 1 ), &
      & ( myTimes(2) - myTimes(1) + 1 ) )

    l2gp%name               = ol2gp%name    
    l2gp%nameIndex          = ol2gp%nameIndex    
    ! l2gp%quantitytype       = ol2gp%quantitytype 
    l2gp%MissingL2GP        = ol2gp%MissingL2GP
    l2gp%MissingValue       = ol2gp%MissingValue 
    l2gp%MissingStatus      = ol2gp%MissingStatus
    l2gp%verticalCoordinate = ol2gp%verticalCoordinate    
    ! Now fill the actual arrays
    ! We'll need to have sane ranges for freqs and levels even if the
    ! we supplied 0 originally
    myFreqs(2) = max( myFreqs(2), myFreqs(1) )
    myLevels(2) = max( myLevels(2), myLevels(1) )
    ! convert to hyperslab params start, count, stride, block
    start(1) = myFreqs(1) + MLS_HyperStart -1 
    start(2) = mylevels(1) + MLS_HyperStart -1 
    start(3) = myTimes(1) + MLS_HyperStart -1 
    count(1) = myFreqs(2) - myFreqs(1) + 1
    count(2) = mylevels(2) - mylevels(1) + 1
    count(3) = myTimes(2) - myTimes(1) + 1
    if ( ol2gp%nFreqs < 1 ) then
      start(1) = MLS_HyperStart
      count(1) = 1
    endif
    call ExtractArray ( l2gp%pressures    , ol2gp%pressures    , start(2:2), count(2:2), stride(2:2), block(2:2) )
    call ExtractArray ( l2gp%frequency    , ol2gp%frequency    , start(1:1), count(1:1), stride(1:1), block(1:1) )
    call ExtractArray ( l2gp%latitude     , ol2gp%latitude     , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%longitude    , ol2gp%longitude    , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%solarTime    , ol2gp%solarTime    , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%solarZenith  , ol2gp%solarZenith  , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%losAngle     , ol2gp%losAngle     , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%geodAngle    , ol2gp%geodAngle    , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%time         , ol2gp%time         , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%chunkNumber  , ol2gp%chunkNumber  , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%l2gpValue    , ol2gp%l2gpValue    , start, count, stride, block )
    call ExtractArray ( l2gp%l2gpPrecision, ol2gp%l2gpPrecision, start, count, stride, block )
    call ExtractArray ( l2gp%status       , ol2gp%status       , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%quality      , ol2gp%quality      , start(3:3), count(3:3), stride(3:3), block(3:3) )
    call ExtractArray ( l2gp%convergence  , ol2gp%convergence  , start(3:3), count(3:3), stride(3:3), block(3:3) )
    ! We may implement this field as either an integer or a char
    ! call ExtractArray ( l2gp%AscDescMode    , ol2gp%AscDescMode    , start(2:2), count(2:2), stride(2:2), block(2:2) )
    if ( AscDescModeIsField ) then
      do i=1, count(3)
        l2gp%AscDescMode(i) = ol2gp%AscDescMode(myTimes(1) + i - 1)
      enddo
    endif
    if ( associated(ol2gp%BinNumber) ) then
      allocate(l2gp%BinNumber(1:count(3)))
      do i=1, count(3)
        l2gp%BinNumber(i) = ol2gp%BinNumber(myTimes(1) + i - 1)
      enddo
    endif
    if ( associated(ol2gp%MAF) ) then
      allocate(l2gp%MAF(1:count(3)))
      do i=1, count(3)
        l2gp%MAF(i) = ol2gp%MAF(myTimes(1) + i - 1)
      enddo
    endif
  end subroutine ExtractL2GPRecord_range

  subroutine ExtractL2GPRecord_which ( ol2gp, l2gp, &
    & whichFreqsIn, whichLevels, whichTimes, &
    & useFreqs, useLevels, useTimes )
    ! Like the above but with sets of which indexes to choose
    ! instead of their ranges
    ! Dummy arguments
    type (L2GPData_T), intent(in)   :: ol2gp
    type (L2GPData_T), intent(out)  :: l2gp
    integer, dimension(:), intent(in), target :: whichFreqsIn! which freqs
    integer, dimension(:), intent(in) :: whichLevels ! which levels
    integer, dimension(:), intent(in) :: whichTimes  ! which times
    integer, intent(in)               :: useFreqs    ! how many freqs
    integer, intent(in)               :: useLevels   ! how many levels
    integer, intent(in)               :: useTimes    ! how many times
    ! Local variables
    integer :: i

    integer, dimension(3) :: start
    integer, dimension(3) :: count
    integer, dimension(3) :: block
    integer, dimension(3) :: stride
    integer               :: numFreqs  ! how many freqs
    integer               :: numLevels ! how many levels
    integer               :: numTimes  ! how many times
    integer, dimension(:), pointer :: whichFreqs
    
    ! Executable
    call SetupNewL2GPRecord ( l2gp, &
      & ol2gp%nFreqs, &
      & useLevels, &
      & useTimes )
    
    numFreqs                = max ( useFreqs , 1 )
    numLevels               = max ( useLevels, 1 )
    numTimes                = max ( useTimes , 1 )
    l2gp%name               = ol2gp%name    
    l2gp%nameIndex          = ol2gp%nameIndex    
    ! l2gp%quantitytype       = ol2gp%quantitytype 
    l2gp%MissingL2GP        = ol2gp%MissingL2GP
    l2gp%MissingValue       = ol2gp%MissingValue 
    l2gp%MissingStatus      = ol2gp%MissingStatus
    l2gp%verticalCoordinate = ol2gp%verticalCoordinate    
    if ( usefreqs > 0 ) then
      whichFreqs => whichFreqsIn
    else
      allocate( whichFreqs(1) )
      whichFreqs = 1
    endif
    ! Now fill the actual arrays
    if ( DeeBug ) then
      call outputNamedValue( 'UseFreqs', useFreqs )
      call outputNamedValue( 'UseLevels', useLevels )
      call outputNamedValue( 'UseTimes', useTimes )
      call dump( whichFreqs(1:useFreqs), 'whichFreqs' )
      call dump( whichLevels(1:useLevels), 'whichLevels' )
      call dump( whichTimes(1:useTimes), 'whichTimes' )
    endif
    call GatherBloc ( l2gp%pressures    , ol2gp%pressures    , &
      & whichLevels(1:numLevels) )
    call GatherBloc ( l2gp%frequency    , ol2gp%frequency    , &
      & whichFreqs(1:numFreqs) )
    call GatherBloc ( l2gp%latitude     , ol2gp%latitude     , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%longitude    , ol2gp%longitude    , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%solarTime    , ol2gp%solarTime    , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%solarZenith  , ol2gp%solarZenith  , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%losAngle     , ol2gp%losAngle     , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%geodAngle    , ol2gp%geodAngle    , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%time         , ol2gp%time         , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%chunkNumber  , ol2gp%chunkNumber  , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%l2gpValue    , ol2gp%l2gpValue    , &
      & whichFreqs(1:numFreqs), whichLevels(1:numLevels), whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%l2gpPrecision, ol2gp%l2gpPrecision, &
      & whichFreqs(1:numFreqs), whichLevels(1:numLevels), whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%status       , ol2gp%status       , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%quality      , ol2gp%quality      , &
      & whichTimes(1:numTimes) )
    call GatherBloc ( l2gp%convergence  , ol2gp%convergence  , &
      & whichTimes(1:numTimes) )
    if ( AscDescModeIsField ) then
      do i=1, numTimes
        l2gp%AscDescMode(i) = ol2gp%AscDescMode(whichTimes(i))
      enddo
    endif
    if ( associated(ol2gp%BinNumber) ) then
      allocate(l2gp%BinNumber(1:count(3)))
      do i=1, numTimes
        l2gp%BinNumber(i) = ol2gp%BinNumber(whichTimes(i))
      enddo
    endif
    if ( associated(ol2gp%MAF) ) then
      allocate(l2gp%MAF(1:numTimes))
      do i=1, numTimes
        l2gp%MAF(i) = ol2gp%MAF(whichTimes(i))
      enddo
    endif
    if ( usefreqs == 0 ) then
      deallocate( whichFreqs )
    endif
  end subroutine ExtractL2GPRecord_which

  !------------------------------------------  IsL2GPSetUp  -----
  function IsL2GPSetUp ( l2gp ) result( sooDesu )

    ! return TRUE only if all arrays allocated

    ! Dummy arguments
    type (L2GPData_T), intent(in)  :: l2gp
    logical :: sooDesu
    ! Executable
    sooDesu = .true.
    sooDesu = sooDesu .and. ( size(l2gp%pressures       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%frequency       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%latitude        )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%longitude       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%solarTime       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%solarZenith     )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%losAngle        )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%geodAngle       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%time            )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%chunkNumber     )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%l2gpValue       )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%l2gpPrecision   )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%status          )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%quality         )     > 0 )
    sooDesu = sooDesu .and. ( size(l2gp%convergence     )     > 0 )
    if ( AscDescModeIsField ) &
      & sooDesu = sooDesu .and. ( size(l2gp%AscDescMode       )     > 0 )
  end function IsL2GPSetUp

  !------------------------------------------  SetupNewL2GPRecord  -----
  subroutine SetupNewL2GPRecord ( l2gp, nFreqs, nLevels, nTimes, nTimesTotal, &
    & FillIn, proto, which )

    ! This routine sets up the arrays for an l2gp datatype.
    ! It preassigns their values to MissingValue (if FillIn is TRUE)
    ! or to proptype's values (if supplied)

    ! Dummy arguments
    type (L2GPData_T), intent(inout)  :: l2gp
    integer, intent(in), optional :: nFreqs              ! Dimensions
    integer, intent(in), optional :: nLevels             ! Dimensions
    integer, intent(in), optional :: nTimes              ! Dimensions
    integer, intent(in), optional :: nTimesTotal         ! Dimensions
    logical, intent(in), optional :: FillIn    ! Fill with MissingValue
    type (L2GPData_T), intent(in), optional  :: proto    ! prototype
    integer, dimension(:), optional, intent(in) :: which ! which profiles to fill from

    ! Local variables
    logical :: myFillIn
    integer :: nn
    integer :: useNFreqs, useNLevels, useNTimes, useNTimesTotal

    if ( present(nFreqs) ) then
       useNFreqs=nFreqs
    elseif ( present(proto) ) then
       useNFreqs=proto%nFreqs
    else
       useNFreqs=0
    end if

    if ( present(nLevels) ) then
       useNLevels=nLevels
    elseif ( present(proto) ) then
       useNLevels=proto%nLevels
    else
       useNLevels=0
    end if

    if ( present(nTimes) ) then
       useNTimes=nTimes
    elseif ( present(proto) ) then
      if ( present(which) ) then
       useNTimes=size(which)
      else
       useNTimes=proto%nTimes
      endif
    else
       useNTimes=0              ! Default to empty l2gp
    end if

    if ( present(nTimesTotal) ) then
       useNTimesTotal=nTimesTotal
    elseif ( present(proto) ) then
       useNTimesTotal=proto%nTimesTotal
    else
       useNTimesTotal=0              ! Default to empty l2gp
    end if

    if ( present(FillIn) ) then
      myFillIn = FillIn
    elseif ( present(proto) ) then
     myFillIn=.FALSE. ! fill in with values of proto
    else
      myFillIn = ALWAYSFILLWITHMISSINGVALUE
    endif

    ! Sanity
    useNTimesTotal = max(useNTimesTotal, useNTimes)
    if ( useNLevels > MAXNLEVELS ) &
      & call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Too many levels in SetUpNewL2GPRecord')
    ! Store the dimensionality

    l2gp%nTimes = useNTimes
    l2gp%nTimesTotal = useNTimesTotal
    l2gp%nLevels = useNLevels
    l2gp%nFreqs = useNFreqs

    ! But allocate to at least one for times, freqs
 
    useNLevels=MAX(useNLevels,1)
    useNFreqs=MAX(useNFreqs,1)    

    ! Allocate the vertical coordinate
    l2gp%verticalCoordinate = 'Pressure'
    call allocate_test ( l2gp%pressures, useNLevels, "l2gp%pressures", &
         & ModuleName )

    ! Allocate the frequency coordinate
    
    call allocate_test ( l2gp%frequency, useNFreqs, "l2gp%frequency", ModuleName)

    ! Allocate the horizontal coordinates

    call allocate_test(l2gp%latitude,   useNTimes, "l2gp%latitude",   ModuleName)
    call allocate_test(l2gp%longitude,  useNTimes, "l2gp%longitude",  ModuleName)
    call allocate_test(l2gp%solarTime,  useNTimes, "l2gp%solarTime",  ModuleName)
    call allocate_test(l2gp%solarZenith,useNTimes, "l2gp%solarZenith",ModuleName)
    call allocate_test(l2gp%losAngle,   useNTimes, "l2gp%losAngle",   ModuleName)
    call allocate_test(l2gp%geodAngle,  useNTimes, "l2gp%geodAngle",  ModuleName)
    call allocate_test(l2gp%time,       useNTimes, "l2gp%time",       ModuleName)
    call allocate_test(l2gp%chunkNumber,useNTimes, "l2gp%chunkNumber",ModuleName)

    ! Allocate the data fields

    call allocate_test(l2gp%l2gpValue,useNFreqs,useNLevels,&
      & useNTimes,"l2gp%l2gpValue", ModuleName)
    call allocate_test(l2gp%l2gpPrecision,useNFreqs,useNLevels,&
      & useNTimes,"l2gp%l2gpPrecision", ModuleName)

    call allocate_test(l2gp%status, useNTimes,"l2gp%status", ModuleName)
    call allocate_test(l2gp%quality,useNTimes,"l2gp%quality",ModuleName)
    call allocate_test(l2gp%convergence,useNTimes,"l2gp%convergence",ModuleName)
    call allocate_test(l2gp%AscDescMode,useNTimes,"l2gp%AscDescMode",ModuleName)
    if ( myFillIn ) then
      l2gp%pressures    = l2gp%MissingValue
      l2gp%frequency    = l2gp%MissingValue
      l2gp%latitude     = l2gp%MissingValue
      l2gp%longitude    = l2gp%MissingValue
      l2gp%solarTime    = l2gp%MissingValue
      l2gp%solarZenith  = l2gp%MissingValue
      l2gp%losAngle     = l2gp%MissingValue
      l2gp%geodAngle    = l2gp%MissingValue
      l2gp%time         = l2gp%MissingValue
      l2gp%chunkNumber  = UndefinedIntegerValue
      l2gp%l2gpValue    = l2gp%MissingL2GP
      l2gp%l2gpPrecision= l2gp%MissingValue
      l2gp%status       = l2gp%MissingStatus ! l2gp%MissingValue
      l2gp%quality      = l2gp%MissingValue
      l2gp%convergence  = l2gp%MissingValue
      if ( AscDescModeIsField ) &
        & l2gp%AscDescMode    = 0
    elseif ( present(proto) ) then
      nn                = proto%nTimes
      l2gp%name         = proto%name
      l2gp%nameIndex    = proto%nameIndex
      ! l2gp%QUANTITYTYPE = proto%QUANTITYTYPE
      l2gp%MissingStatus= proto%MissingStatus
      l2gp%MissingL2GP  = proto%MissingL2GP
      l2gp%MissingValue = proto%MissingValue
      l2gp%pressures    = proto%pressures    
      l2gp%frequency    = proto%frequency
      if ( useNTimes > proto%nTimes ) then
        l2gp%latitude       (1:nn) = proto%latitude     
        l2gp%longitude      (1:nn) = proto%longitude    
        l2gp%solarTime      (1:nn) = proto%solarTime    
        l2gp%solarZenith    (1:nn) = proto%solarZenith  
        l2gp%losAngle       (1:nn) = proto%losAngle     
        l2gp%geodAngle      (1:nn) = proto%geodAngle    
        l2gp%time           (1:nn) = proto%time         
        l2gp%chunkNumber    (1:nn) = proto%chunkNumber  
        l2gp%l2gpValue      (:,:,1:nn) = proto%l2gpValue    
        l2gp%l2gpPrecision  (:,:,1:nn) = proto%l2gpPrecision
        l2gp%status         (1:nn) = proto%status       
        l2gp%quality        (1:nn) = proto%quality      
        l2gp%convergence    (1:nn) = proto%convergence  
        if ( AscDescModeIsField ) &
          & l2gp%AscDescMode      (1:nn) = proto%AscDescMode  
        if ( associated(proto%BinNumber) ) then
          allocate( l2gp%BinNumber(useNTimes) )
          l2gp%BinNumber    (1:nn) = proto%BinNumber 
        endif
        if ( associated(proto%MAF) ) then
          allocate( l2gp%MAF(useNTimes) )
          l2gp%MAF    (1:nn) = proto%MAF 
        endif
      elseif ( present( which ) ) then
        l2gp%latitude     = proto%latitude     (which)
        l2gp%longitude    = proto%longitude    (which)
        l2gp%solarTime    = proto%solarTime    (which)
        l2gp%solarZenith  = proto%solarZenith  (which)
        l2gp%losAngle     = proto%losAngle     (which) 
        l2gp%geodAngle    = proto%geodAngle    (which) 
        l2gp%time         = proto%time         (which) 
        l2gp%chunkNumber  = proto%chunkNumber  (which) 
        l2gp%l2gpValue    = proto%l2gpValue    (:,:,which)
        l2gp%l2gpPrecision= proto%l2gpPrecision(:,:,which)
        l2gp%status       = proto%status       (which)
        l2gp%quality      = proto%quality      (which)
        l2gp%convergence  = proto%convergence  (which)
        if ( AscDescModeIsField ) &
          & l2gp%AscDescMode    = proto%AscDescMode  (which)
      else
        l2gp%latitude     = proto%latitude     
        l2gp%longitude    = proto%longitude    
        l2gp%solarTime    = proto%solarTime    
        l2gp%solarZenith  = proto%solarZenith
        l2gp%losAngle     = proto%losAngle     
        l2gp%geodAngle    = proto%geodAngle    
        l2gp%time         = proto%time         
        l2gp%chunkNumber  = proto%chunkNumber  
        l2gp%l2gpValue    = proto%l2gpValue    
        l2gp%l2gpPrecision= proto%l2gpPrecision
        l2gp%status       = proto%status       
        l2gp%quality      = proto%quality      
        l2gp%convergence  = proto%convergence  
        if ( AscDescModeIsField ) &
          & l2gp%AscDescMode    = proto%AscDescMode  
        if ( associated(proto%BinNumber) ) then
          allocate( l2gp%BinNumber(useNTimes) )
          l2gp%BinNumber  = proto%BinNumber 
        endif
        if ( associated(proto%MAF) ) then
          allocate( l2gp%MAF(useNTimes) )
          l2gp%MAF        = proto%MAF 
        endif
      endif
    endif
  end subroutine SetupNewL2GPRecord

  !-----------------------------------------  DestroyL2GPContents  -----
  subroutine DestroyL2GPContents ( L2GP )

    ! This routine deallocates all the arrays allocated above.

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP

    ! Executable code

    call deallocate_test ( l2gp%pressures,     "l2gp%pressures",     ModuleName )
    call deallocate_test ( l2gp%latitude,      "l2gp%latitude",      ModuleName )
    call deallocate_test ( l2gp%longitude,     "l2gp%longitude",     ModuleName )
    call deallocate_test ( l2gp%solarTime,     "l2gp%solarTime",     ModuleName )
    call deallocate_test ( l2gp%solarZenith,   "l2gp%solarZenith",   ModuleName )
    call deallocate_test ( l2gp%losAngle,      "l2gp%losAngle",      ModuleName )
    call deallocate_test ( l2gp%geodAngle,     "l2gp%geodAngle",     ModuleName )
    call deallocate_test ( l2gp%chunkNumber,   "l2gp%chunkNumber",   ModuleName )
    call deallocate_test ( l2gp%time,          "l2gp%time",          ModuleName )
    call deallocate_test ( l2gp%frequency,     "l2gp%frequency",     ModuleName )
    call deallocate_test ( l2gp%l2gpValue,     "l2gp%l2gpValue",     ModuleName )
    call deallocate_test ( l2gp%l2gpPrecision, "l2gp%l2gpPrecision", ModuleName )
    call deallocate_test ( l2gp%status,        "l2gp%status",        ModuleName )
    call deallocate_test ( l2gp%quality,       "l2gp%quality",       ModuleName )
    call deallocate_test ( l2gp%convergence,   "l2gp%convergence",   ModuleName )
    call deallocate_test ( l2gp%AscDescMode,   "l2gp%AscDescMode",   ModuleName )
    call deallocate_test ( l2gp%BinNumber,     "l2gp%BinNumber",     ModuleName )
    call deallocate_test ( l2gp%MAF,           "l2gp%MAF        ",   ModuleName )
    l2gp%nTimes = 0
    l2gp%nTimesTotal = 0
    l2gp%nLevels = 0
    l2gp%nFreqs = 0

  end subroutine DestroyL2GPContents

  !-------------------------------------- ExpandL2GPDataInFile ---
  
  subroutine ExpandL2GPDataInFile( L2GPFile, swathName, l2gp )
    ! Repeated rewrite swathname until its numProfs >= l2gp%nTimes
    ! Args
    type(MLSFile_T) :: L2GPFile
    character(len=*), intent(in) :: swathName
    type (L2GPData_T), intent(inout) :: l2gp
    ! Internal variables
    logical :: done
    integer :: M
    integer :: ReadM
    integer, parameter :: MAXPASSES = 40
    integer, parameter :: METHOD = 1
    integer :: numProfs
    integer :: passes
    type (L2GPData_T) :: smallerl2gp, readl2gp
    ! Executable
    if ( L2GPFile%hdfVersion /= HDFVERSION_5 ) return
    M = l2gp%nTimes
    if ( M < 2 ) return
    if ( method == 1 ) then
      call ReadL2GPData( L2GPFile, swathName, readl2gp, numProfs, &
        & readData=.false. )
      ReadM = readl2gp%nTimes
      if(DEEBUG) call outputNamedValue( 'numProfs', numProfs )
      if(DEEBUG) call outputNamedValue( 'readl2gp%nTimes', readl2gp%nTimes )
      call DestroyL2GPContents ( readl2gp )
      if ( ReadM >= M ) return
      call ExtractL2GPRecord ( l2gp, smallerl2gp, rTimes=(/ M, M /) )
      call OutputL2GP_writeGeo_MF ( smallerl2gp, l2GPFile, &
        & swathName, offset=M-1 )
      call OutputL2GP_writeData_MF ( smallerl2gp, l2GPFile, &
        & swathName, offset=M-1 )
      call DestroyL2GPContents ( smallerl2gp )
      return
    endif
    done = .false.
    do passes=1, MAXPASSES
      call ReadL2GPData( L2GPFile, swathName, readl2gp, numProfs )
      ReadM = readl2gp%nTimes
      if(DEEBUG) call outputNamedValue( 'numProfs', numProfs )
      if(DEEBUG) call outputNamedValue( 'readl2gp%nTimes', readl2gp%nTimes )
      done = ( ReadM >= M )
      call DestroyL2GPContents ( readl2gp )
      if ( done ) return
      call ExtractL2GPRecord ( l2gp, smallerl2gp, rTimes=(/ ReadM+1, M /) )
      call OutputL2GP_writeGeo_MF ( smallerl2gp, l2GPFile, &
        & swathName, offset=ReadM )
      call OutputL2GP_writeData_MF ( smallerl2gp, l2GPFile, &
        & swathName, offset=ReadM )
      call DestroyL2GPContents ( smallerl2gp )
    enddo
    call MLSMessage ( MLSMSG_Error, ModuleName // '/ExpandL2GPDataInFile', &
      & "Failed to fit l2gp in file", MLSFile=L2GPFile )
  end subroutine ExpandL2GPDataInFile
  
  !---------------------------------------  ExpandL2GPDataInPlace  -----
  
  subroutine ExpandL2GPDataInPlace ( l2gp, newNTimes )

    ! This subroutine expands an L2GPData_T in place allowing the user to 
    ! (1) add more profiles to it; or
    ! 

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: l2gp
    integer, optional, intent(in)    :: newNTimes

    ! Local variables
    type (L2GPData_T) :: tempL2gp       ! For copying data around
    integer :: myNTimes
    integer :: tmpNFreqs, tmpNLevels

    ! Executable code

   if(present(newNTimes)) then
      myNTimes = newNTimes
   else
      myNTimes = l2gp%nTimes
   endif

    ! First do a sanity check


    if ( myNTimes<l2gp%nTimes ) then
          call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "The number of profiles requested is fewer than those already present" )
    endif

    tempL2gp = l2gp ! Copy the pointers to the old information

    ! Now recreate l2gp with the new size.
    ! First, nullify all of the pointers in l2gp, so that a deallocate_test
    ! won't delete them.  After all, we just went to the trouble to preserve
    ! them in TempL2GP!
    nullify ( l2gp%pressures, l2gp%latitude, l2gp%longitude, l2gp%solarTime, &
      & l2gp%solarZenith, l2gp%losAngle, l2gp%losAngle, l2gp%geodAngle, &
      & l2gp%chunkNumber, l2gp%time, l2gp%frequency, l2gp%l2gpValue, &
      & l2gp%l2gpPrecision, l2gp%status, l2gp%quality, l2gp%convergence, &
      & l2gp%AscDescMode )
    
    tmpNFreqs = l2gp%nFreqs
    tmpNLevels = l2gp%nLevels
    call SetupNewL2GPRecord( l2gp, nFreqs=tmpNFreqs, nLevels=tmpNLevels, &
      & nTimes=myNTimes, FillIn = .true. )

    ! Don't forget the `global' stuff
    l2gp%pressures=templ2gp%pressures
    l2gp%frequency=templ2gp%frequency
    l2gp%verticalCoordinate = templ2gp%verticalCoordinate

    ! Now go through the parameters one by one, and copy the previous contents
    l2gp%latitude(1:templ2gp%nTimes) = templ2gp%latitude(1:templ2gp%nTimes)
    l2gp%longitude(1:templ2gp%nTimes) = templ2gp%longitude(1:templ2gp%nTimes)
    l2gp%solarTime(1:templ2gp%nTimes) = templ2gp%solarTime(1:templ2gp%nTimes)
    l2gp%solarZenith(1:templ2gp%nTimes) = templ2gp%solarZenith(1:templ2gp%nTimes)
    l2gp%losAngle(1:templ2gp%nTimes) = templ2gp%losAngle(1:templ2gp%nTimes)
    l2gp%geodAngle(1:templ2gp%nTimes) = templ2gp%geodAngle(1:templ2gp%nTimes)
    l2gp%time(1:templ2gp%nTimes) = templ2gp%time(1:templ2gp%nTimes)
    l2gp%chunkNumber(1:templ2gp%nTimes) = templ2gp%chunkNumber(1:templ2gp%nTimes)

    l2gp%l2gpValue(:,:,1:templ2gp%nTimes) = templ2gp%l2gpValue(:,:,1:templ2gp%nTimes)
    l2gp%l2gpPrecision(:,:,1:templ2gp%nTimes) = &
         templ2gp%l2gpPrecision(:,:,1:templ2gp%nTimes)
    
    l2gp%status(1:templ2gp%nTimes) = templ2gp%status(1:templ2gp%nTimes)
    l2gp%quality(1:templ2gp%nTimes) = templ2gp%quality(1:templ2gp%nTimes)
    l2gp%convergence(1:templ2gp%nTimes) = templ2gp%convergence(1:templ2gp%nTimes)
    if ( AscDescModeIsField ) &
      & l2gp%AscDescMode(1:templ2gp%nTimes) = templ2gp%AscDescMode(1:templ2gp%nTimes)
    
    ! Deallocate the old arrays
    call DestroyL2GPContents(templ2gp)

  end subroutine ExpandL2GPDataInPlace

  !-------------------------------------------  AddL2GPToDatabase  -----
  integer function AddL2GPToDatabase( DATABASE, ITEM )

    ! This function adds an l2gp data type to a database of said types,
    ! creating a new database if it doesn't exist.  The result value is
    ! the size -- where L2gp is put.

    use Allocate_Deallocate, only: Test_Allocate, Test_Deallocate

    ! Dummy arguments
    type (l2gpdata_t), dimension(:), pointer :: DATABASE
    type (l2gpdata_t), intent(in) :: ITEM

    ! Local variables
    type (l2GPData_T), dimension(:), pointer :: tempDatabase
    !This include causes real trouble if you are compiling in a different 
    !directory.
    include "addItemToDatabase.f9h" 

    AddL2GPToDatabase = newSize
  end function AddL2GPToDatabase

  ! --------------------------------------------------------------------------

  ! This subroutine destroys a quantity template database

  subroutine DestroyL2GPDatabase ( DATABASE )

    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    use Toggles, only: Gen, Toggle
    use Trace_m, only: Trace_Begin, Trace_End
    ! Dummy argument
    type (l2GPData_T), dimension(:), pointer :: DATABASE

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: l2gpIndex, s, status
    integer :: Me = -1       ! String index for trace

    call trace_begin ( me, "DestroyL2GPDatabase", cond=toggle(gen) )

    if ( associated(database) ) then
       do l2gpIndex = 1, SIZE(database)
          call DestroyL2GPContents ( database(l2gpIndex) )
       end do
       s = size(database) * storage_size(database) / 8
       addr = 0
       if ( s > 0 ) addr = transfer(c_loc(database(1)), addr)
       deallocate ( database, stat=status )
       call test_deallocate ( status, ModuleName, "database", s, Address=addr )
    end if
    call trace_end ( "DestroyL2GPDatabase", cond=toggle(gen) )
  end subroutine DestroyL2GPDatabase

  ! ---------------------- ReadL2GPData_fileID  -----------------------------

  subroutine ReadL2GPData_fileID(L2FileHandle, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion, HMOT, ReadData)
    !------------------------------------------------------------------------

    ! Given a file handle,
    ! This routine reads an L2GP file, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! if present, hdfVersion must be one of HDFVERSION_4, HDFVERSION_5

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    integer, intent(in) :: L2FileHandle ! Returned by swopen
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    integer, optional, intent(in) :: hdfVersion
    character, optional, intent(in) :: hmot   ! 'H' 'M'(def) 'O' 'T'
    logical, optional, intent(in) :: ReadData

    ! Local
    integer :: myhdfVersion
    integer :: status
    type( MLSFile_T ) :: l2gpFile
    
    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif
    status = InitializeMLSFile(l2gpFile, type=l_swath, access=DFACC_RDONLY, &
      & content='l2gp', name='unknown', hdfVersion=myhdfVersion)
    l2gpFile%FileID%f_id = l2FileHandle
    l2gpFile%stillOpen = .true.
    call ReadL2GPData(l2gpFile, swathname, l2gp, numProfs, &
       firstProf, lastProf, HMOT, ReadData)
  end subroutine ReadL2GPData_fileID

  ! ---------------------- ReadL2GPData_fileName  -----------------------------

  subroutine ReadL2GPData_fileName(fileName, swathname, l2gp, numProfs, &
       firstProf, lastProf, hdfVersion, HMOT, ReadData)
    !------------------------------------------------------------------------

    ! Given a file name,
    ! This routine reads an L2GP file, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! hdfVersion may be WILDCARDHDFVERSION

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    character (len=*), intent(in) :: fileName ! Name of swath
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    integer, optional, intent(in) :: hdfVersion
    character, optional, intent(in) :: hmot   ! 'H' 'M'(def) 'O' 'T'
    logical, optional, intent(in) :: ReadData

    ! Local
    integer :: L2FileHandle
    type(MLSFile_T)                :: MLSFile
    integer :: status
    integer :: the_hdfVersion
    
    ! Executable code
    the_hdfVersion = mls_hdf_version(FileName, hdfVersion)
    if ( the_hdfVersion == FILENOTFOUND ) &
      call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'File not found; make sure the name and path are correct' &
        & // trim(fileName) )

    status = InitializeMLSFile ( MLSFile, type=l_swath, access=DFACC_READ, &
     & name=trim(fileName), HDFVersion=the_hdfVersion )
    call MLS_OpenFile( MLSFile )
    L2FileHandle = MLSFile%FileID%f_id
    call ReadL2GPData_fileID(L2FileHandle, swathname, l2gp, numProfs=numProfs, &
       & firstProf=firstProf, lastProf=lastProf, hdfVersion=the_hdfVersion, &
       & hmot=hmot, ReadData=ReadData)
    call MLS_CloseFile( MLSFile )
  end subroutine ReadL2GPData_fileName

  ! ---------------------- ReadL2GPData_MLSFile  -----------------------------

  subroutine ReadL2GPData_MLSFile(L2GPFile, swathname, l2gp, numProfs, &
       firstProf, lastProf, HMOT, ReadData)
    !------------------------------------------------------------------------

    ! Given a file,
    ! This routine reads an L2GP structure, in either hdfVersion,
    ! returning a filled data structure and the !
    ! number of profiles read.
    ! if present, hdfVersion must be one of HDFVERSION_4, HDFVERSION_5

    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    type(MLSFile_T)                :: L2GPFile
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( l2GPData_T ), intent(out) :: l2gp ! Result
    integer, intent(out), optional :: numProfs ! Number actually read
    character, optional, intent(in) :: hmot   ! 'H' 'M'(def) 'O' 'T'
    logical, optional, intent(in) :: ReadData

    ! Local
    logical :: alreadyOpen
    character :: my_hmot
    integer :: status
    
    ! Executable code
    my_hmot = 'M'
    if ( present(hmot)) my_hmot = Capitalize(hmot)
    ! Check for valid hmot
    select case (my_hmot)
    case ('H')
    case ('M')
    case ('O')
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unable to Read OMI L2GPData", MLSFile=L2GPFile )
    case ('T')
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unable to Read TES L2GPData", MLSFile=L2GPFile )
    case default
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized instrument key passed to ReadL2GPData: "// my_hmot, &
      & MLSFile=L2GPFile )
    end select
    status = 0
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    if (L2GPFile%hdfVersion == HDFVERSION_4) then
      call ReadL2GPData_MF_hdf(L2GPFile, swathname, l2gp, my_hmot,&
        & numProfs, firstProf, lastProf, ReadData)
    elseif (L2GPFile%hdfVersion /= HDFVERSION_5) then
      call MLSMessage ( MLSMSG_Error, ModuleName, &
      & "Unrecognized hdfVersion passed to ReadL2GPData", MLSFile=L2GPFile )
    else
      call ReadL2GPData_MF_hdf(L2GPFile, swathname, l2gp, my_hmot,&
        & numProfs, firstProf, lastProf, ReadData)
    endif
    if ( .not. alreadyOpen )  call mls_closeFile(L2GPFile, Status)
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'read'
  end subroutine ReadL2GPData_MLSFile

  ! ------------------- ReadL2GPData_MF_hdf ----------------

  subroutine ReadL2GPData_MF_hdf(L2GPFile, swathname, l2gp, HMOT, &
    & numProfs, firstProf, lastProf, ReadData)
  use HDFEOS, only: SWInqdims
  use HDFEOS5, only: HE5_SWInqdims, HE5_SWInqdflds, HE5_SWFldinfo
  use HE5_SWAPI, only: HE5_SWRdlattr
  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach, MLS_SWDiminfo, MLS_SWRdfld
  use MLSStringLists, only: IsInList
  use HDF5, only: Size_T
    !------------------------------------------------------------------------

    ! This routine reads an L2GP file, returning a filled data structure and the !
    ! number of profiles read.

    ! All the ReadConvergence harrumphing is because convergence is newly added
    ! (Some older files won't have this field)
    ! Arguments

    character (len=*), intent(in) :: swathname ! Name of swath
    type(MLSFile_T)                :: L2GPFile
    integer, intent(in), optional :: firstProf, lastProf ! Defaults to first and last
    type( L2GPData_T ), intent(out) :: l2gp ! Result
    character, intent(in) :: HMOT   ! 'H' 'M'(def) 'O' 'T'
    integer, intent(out),optional :: numProfs ! Number actually read
    logical, optional, intent(in) :: ReadData

    ! Local Parameters
    character (len=*), parameter :: MLSMSG_INPUT = 'Error in input argument '
    integer, parameter           :: MAXDIMSIZE = 4
    ! Local Variables
    character (len=80) :: DF_Name
    character (len=80) :: DF_Precision
    character (len=80) :: dimlist
    character (len=80) :: fieldlist
    character (len=80) :: list
    character (len=80) :: maxdimlist
    character (len=8)  :: maxdimName
    integer :: hdfVersion
    integer :: rank
    integer, dimension(MAXFNFIELDS) :: ranks
    integer, dimension(MAXFNFIELDS) :: types
    integer, dimension(7) :: numberType
    integer, dimension(7) :: flddims
    integer(kind=size_t), dimension(7) :: hflddims
    character (len=480) :: msr

    integer, dimension(MAXDIMSIZE) :: dims
    integer :: first, freq, lev, nDims, nFlds, size, swid, status
    integer(kind=size_t), dimension(MAXDIMSIZE) :: hdims
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: nFreqs, nLevels, nTimes, nFreqsOr1, nLevelsOr1, myNumProfs
    integer, dimension(3) :: start, stride, edge
    logical :: firstCheck, lastCheck

    real(r4), pointer, dimension(:) :: REALFREQ
    real(r4), pointer, dimension(:) :: REALSURF
    real(r4), pointer, dimension(:) :: REALPROF
    real(r4), pointer, dimension(:,:,:) :: REAL3
    logical :: deeBugHere
    logical :: dontfail
    logical :: ReadingConvergence
    logical :: ReadingAscDescMode
    logical :: ReadingData
    logical :: ReadingBinNumber
    logical :: ReadingMAF
    real(rgp), dimension(1) :: MissingValue
    ! Executable code
    call trace_begin ( me,  'ReadL2GPData_MF_hdf', cond=.false. )
    deeBugHere = DEEBUG ! .or. .true.
    nullify ( realFreq, realSurf, realProf, real3 )
    hdfVersion = L2GPFile%hdfVersion
    flddims = 0
    ! Don't fail when trying to read an mls-specific field 
    ! if the file is from another Aura instrument
    dontfail = (HMOT /= 'M')
    ReadingData = .true.
    if ( present(ReadData) ) ReadingData = ReadData
    ReadingConvergence = .false.
    ReadingAscDescMode = .false.
    ReadingBinNumber   = .false.
    ReadingMAF         = .false.
    ! Attach to the swath for reading
    l2gp%Name = swathname
    ! We have suffered surprises when hefeos character fields 
    ! were not initialized
    dimlist = ' '
    fieldlist = ' '
    list = ' '
    
    select case (HMOT)
    case ('H')
      swid = mls_SWattach( L2GPFile, 'HIRDLS' )
      DF_Name = TRIM(l2gp%Name)
      DF_Precision = TRIM(l2gp%Name) // 'Precision'
      l2gp%MissingL2GP = -999.  ! This is a HIRDLS-specific setting
    case ('M')
      swid = mls_SWattach( L2GPFile, l2gp%Name )
      DF_Name = DATA_FIELD1
      DF_Precision = DATA_FIELD2
      if ( deeBugHere ) print *, 'DF_NAME: ',DF_NAME
      if ( deeBugHere ) print *, 'DF_Precision: ',DF_Precision
      ! Here we read the MissingValue attribute
      ! If we can't for any reason, it retains its default value
      if ( hdfVersion == HDFVERSION_5 ) then
        status = &
          & he5_swrdlattr( swid, 'L2gpValue', 'MissingValue', MissingValue )
        if ( status /= 0 ) then
          call MLSMessage( MLSMSG_Warning, ModuleName, &
         &'Failed to read Missing Value attribute for ' &
         & // trim(swathname), MLSFile=L2GPFile )
        else
          l2gp%MissingL2GP = MissingValue(1)
        endif
      endif
    case default
    end select
    if (swid == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
         &'Failed to attach to hdfeos2/5 swath interface for reading' &
         & // trim(swathname), MLSFile=L2GPFile)

    ! Get dimension information

    lev = 0
    freq = 0

    if( hdfVersion == HDFVERSION_4 ) then
      nDims = swinqdims(swid, list, dims)
    else
      nDims = HE5_SWinqdims(swid, list, hdims)
      nDims = min( nDims, MAXDIMSIZE )
      dims(1:nDims) = hdims(1:nDims)                                  
      if ( nDims < MAXDIMSIZE ) dims(nDims+1:) = 0 ! Just to make sure they're defined, not junk
      nFlds = HE5_SWinqdflds( swid, fieldlist, ranks, types )
      ReadingConvergence = isInList( lowerCase(fieldList), 'convergence', '-fc' )
      ReadingAscDescMode = isInList( lowerCase(fieldList), 'AscDescMode', '-fc' )
      ReadingBinNumber   = isInList( lowerCase(fieldList), 'BinNumber'  , '-fc' )
      ReadingMAF         = isInList( lowerCase(fieldList), 'MAF'        , '-fc' )
    endif
    if ( deeBugHere ) print *, 'HMOT: ', HMOT
    if ( deeBugHere ) print *, 'swathName: ', l2gp%name
    if ( deeBugHere ) print *, 'dimlist: ', trim(list)
    if ( deeBugHere ) print *, 'ndims: ', ndims
    if ( deeBugHere ) print *, 'DF_NAME: ',DF_NAME
    if ( deeBugHere ) print *, 'DF_Precision: ',DF_Precision
    if ( deeBugHere ) print *, 'DATA_FIELD1: ', DATA_FIELD1
    if ( deeBugHere ) print *, 'dims: ', dims
    if ( deeBugHere ) print *, 'swathName: ', l2gp%name
    if ( deeBugHere ) print *, 'fieldlist: ', trim(fieldlist)
    if ( deeBugHere ) print *, 'ReadBin:   ', ReadingBinNumber
    if ( deeBugHere ) print *, 'ReadMAF:   ', ReadingMAF      
    if (nDims == -1) call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'Failed to get dimension information on hdfeos5 swath ' // &
      & trim(swathname), MLSFile=L2GPFile)
    if ( index(list,'nLevels') /= 0 ) lev = 1
    if ( index(list,'Freq') /= 0 ) freq = 1

    size = mls_swdiminfo(swid, 'nTimes', hdfVersion=hdfVersion)
    if ( deeBugHere ) print *, 'size--1st try: ', size
    if ( hdfVersion == HDFVERSION_5 ) then
      ! This will be wrong if timeIsUnlim .eq. .TRUE. . 
      ! HE5_SWdiminfo returns 1 instead of the right answer.
      status = HE5_swfldinfo(swid, trim(DF_Name), rank, hflddims, &
        & numberType, dimlist, maxdimlist)
      call GetHashElement (dimlist, &
       & maxdimlist, 'nTimes', &
       & maxDimName, .false.)
      if ( deeBugHere ) print *, 'flddims: ', flddims
      if ( deeBugHere ) print *, 'dimlist: ', trim(dimlist)
      if ( deeBugHere ) print *, 'maxdimlist: ', trim(maxdimlist)
      if ( maxDimName == 'Unlim' ) then
        status = StringElementNum(dimlist, 'nTimes', .false.)
        if ( status > 0 .and. status <= rank ) then
          flddims(1:rank) = hflddims(1:rank)
          size = max(size, flddims(status))
        endif
      endif
    endif
    l2gp%nTimes = size
    l2gp%nTimesTotal = size
    nTimes=size

    if (lev == 0) then
       nLevels = 0
    else
      size = mls_swdiminfo(swid, 'nLevels', hdfVersion=hdfVersion)
       nLevels = size

    endif

    if ( freq == 1 .and. HMOT == 'M') then
      size = mls_swdiminfo(swid, 'nFreqs', hdfVersion=hdfVersion)
       nFreqs = size
    elseif ( freq == 1 .and. HMOT == 'H') then
      size = mls_swdiminfo(swid, 'nChans', hdfVersion=hdfVersion)
      nFreqs = size
    else
       nFreqs = 0
    endif

    ! Check optional input arguments

    firstCheck = present(firstProf)
    lastCheck = present(lastProf)

    if (firstCheck) then
      ! Note that if time is an umlimited dimension, HDF-EOS won't 
      ! nTimes is wrong.
       if ( (firstProf >= l2gp%nTimes) &
         .or. (firstProf < 0) ) then
          msr = MLSMSG_INPUT // 'firstProf'
          call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2GPFile)
       else
          first = firstProf
       endif

    else

       first = 0

    endif

    if (lastCheck) then

       if (lastProf < first) then
          msr = MLSMSG_INPUT // 'lastProf'
          call MLSMessage(MLSMSG_Error, ModuleName, msr, MLSFile=L2GPFile)
       endif

       ! If user has supplied "last" _and_ time is unlimited, we have
       ! to believe the user about how many profiles there are.
       ! This is _crap_ and is a temporary workaround. 
!       if(timeIsUnlim) then
!         myNumProfs = lastProf - first + 1
!         nTimes=lastprof-first+1
!         l2gp%nTimes=nTimes
!       endif
       if (lastProf >= nTimes) then
          myNumProfs = nTimes - first
       else
          myNumProfs = lastProf - first + 1
       endif

    else

       myNumProfs = l2gp%nTimes - first

    endif

    ! Allocate result
    if ( deeBugHere ) then
      print *, 'nFreqs: ', nFreqs
      print *, 'nLevels: ', nLevels
      print *, 'mynumProfs: ', mynumProfs
    endif
    call SetupNewL2GPRecord (l2gp, nFreqs=nFreqs, nLevels=nLevels, &
      &  nTimes=mynumProfs, FillIn = .true.)

    ! Skip reading any of the data?
    if ( .not. readingData ) then
      status = mls_SWdetach(swid, hdfVersion=hdfVersion)
      if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
           &detach from swath interface after reading.', MLSFile=L2GPFile)
      if (present(numProfs)) numProfs=myNumProfs
      call trace_end (  'ReadL2GPData_MF_hdf', cond=.false. )
      if(DEEBUG) call outputNamedValue( 'no data read; nTimes', myNumProfs )
      return
    endif
    ! Allocate temporary arrays

    nFreqsOr1=max(nFreqs,1)
    nLevelsOr1=max(nLevels, 1)

    call Allocate_test ( realProf, myNumProfs, 'realProf', ModuleName )
    call Allocate_test ( realSurf, l2gp%nLevels, 'realSurf', ModuleName )
    call Allocate_test ( realFreq, l2gp%nFreqs, 'realFreq', ModuleName )
    call Allocate_test ( real3, nFreqsOr1, nLevelsOr1, myNumProfs, 'real3', ModuleName )

    ! Read the horizontal geolocation fields

    start(1) = 0
    start(2) = 0
    start(3) = first
    stride = 1
    edge(1) = nFreqsOr1
    edge(2) = nLevelsOr1
    edge(3) = myNumProfs
    status=0

    status = mls_SWrdfld(swid, 'Latitude', start(3:3), stride(3:3), &
      edge(3:3), realProf, hdfVersion=hdfVersion)
    l2gp%latitude = realProf

    status = mls_SWrdfld(swid, 'Longitude', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=hdfVersion)
    l2gp%longitude = realProf

    status = mls_SWrdfld(swid, 'Time', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%time, hdfVersion=hdfVersion)

    status = mls_SWrdfld(swid, 'LocalSolarTime', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=hdfVersion)
    l2gp%solarTime = realProf

    status = mls_SWrdfld(swid, 'SolarZenithAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=hdfVersion)
    l2gp%solarZenith = realProf

    ! These next 3 are MLS-specific
    if ( HMOT == 'M' ) then
    status = mls_SWrdfld(swid, 'LineOfSightAngle', start(3:3), stride(3:3), edge(3:3),&
      &    realProf, hdfVersion=hdfVersion, dontfail=dontfail)
    l2gp%losAngle = realProf

    status = mls_SWrdfld(swid, 'OrbitGeodeticAngle', start(3:3), stride(3:3), edge(3:3),&
      &   realProf, hdfVersion=hdfVersion, dontfail=dontfail)
    l2gp%geodAngle = realProf

    status = mls_SWrdfld(swid, 'ChunkNumber', start(3:3), stride(3:3), edge(3:3),&
      &    l2gp%chunkNumber, hdfVersion=hdfVersion, dontfail=dontfail)

    endif
    ! Read the pressures vertical geolocation field, if it exists

    if (lev /= 0) then
       status = mls_SWrdfld(swid,'Pressure',start(2:2),stride(2:2), edge(2:2),&
         & realSurf, hdfVersion=hdfVersion, dontfail=dontfail)
       l2gp%pressures = realSurf
    endif

    ! Read the frequency geolocation field, if it exists

    if (freq == 1) then

       edge(1) = l2gp%nFreqs

       status = mls_SWrdfld(swid,'Frequency',start(1:1),stride(1:1),edge(1:1),&
         & realFreq, hdfVersion=hdfVersion, dontfail=dontfail)
       l2gp%frequency = realFreq

    endif

    ! Read the data fields that may have 1-3 dimensions

    if ( freq == 1) then

       status = mls_SWrdfld(swid, trim(DF_Name), start, stride, edge, real3, &
         & hdfVersion=hdfVersion)
       l2gp%l2gpValue = real3

       status = mls_SWrdfld(swid, trim(DF_Precision), start, stride, edge, real3, &
         & hdfVersion=hdfVersion)
       l2gp%l2gpPrecision = real3

    else if ( lev == 1) then

       status = mls_SWrdfld( swid, trim(DF_Name), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:), hdfVersion=hdfVersion )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(2:3), stride(2:3), &
            edge(2:3), real3(1,:,:), hdfVersion=hdfVersion )
       l2gp%l2gpPrecision = real3

    else

       status = mls_SWrdfld(swid,trim(DF_Name),start(3:3),stride(3:3),edge(3:3),&
         &   real3(1,1,:), hdfVersion=hdfVersion )
       l2gp%l2gpValue = real3

       status = mls_SWrdfld( swid, trim(DF_Precision), start(3:3), stride(3:3), edge(3:3), &
            real3(1,1,:), hdfVersion=hdfVersion )
       l2gp%l2gpPrecision = real3

    endif
    
    ! Read the data fields that are 1-dimensional

    l2gp%status = l2gp%MissingStatus ! l2gp%MissingValue ! So it has a value.
    status = mls_swrdfld( swid, 'Status',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%status, hdfVersion=hdfVersion, dontfail=.true. )

    if ( HMOT == 'M' ) then
      status = mls_SWrdfld(swid, 'Quality', start(3:3), stride(3:3),&
        edge(3:3),realProf, hdfVersion=hdfVersion, dontfail=dontfail)
      l2gp%quality = realProf
    endif

    l2gp%Convergence = l2gp%MissingValue ! l2gp%MissingValue ! So it has a value.
    if ( ReadingConvergence ) &
      & status = mls_swrdfld( swid, 'Convergence',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%convergence, hdfVersion=hdfVersion, dontfail=.true. )

    if ( AscDescModeIsField ) then
      l2gp%AscDescMode = 0 ! l2gp%MissingValue ! So it has a value.
      if ( ReadingAscDescMode ) &
        & status = mls_swrdfld( swid, 'AscDescMode ',start(3:3),stride(3:3),edge(3:3),&
        & l2gp%AscDescMode, hdfVersion=hdfVersion, dontfail=.true. )
    endif

    if ( ReadingBinNumber ) then
      allocate( l2gp%BinNumber(myNumProfs) )
      status = mls_swrdfld( swid, 'BinNumber',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%BinNumber, hdfVersion=hdfVersion, dontfail=.true. )
    endif

    if ( ReadingMAF ) then
      allocate( l2gp%MAF(myNumProfs) )
      status = mls_swrdfld( swid, 'MAF',start(3:3),stride(3:3),edge(3:3),&
      & l2gp%MAF, hdfVersion=hdfVersion, dontfail=.true. )
    endif

    ! Deallocate local variables
    call Deallocate_test ( realProf, 'realProf', ModuleName )
    call Deallocate_test ( realSurf, 'realSurf', ModuleName )
    call Deallocate_test ( realFreq, 'realFreq', ModuleName )
    call Deallocate_test ( real3, 'real3', ModuleName )

    !  After reading, detach from HE5_SWath interface
    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if (status == -1) call MLSMessage(MLSMSG_Error, ModuleName, 'Failed to &
         &detach from swath interface after reading.', MLSFile=L2GPFile)
    !print*," leaving ReadL2GPData_hdf: first/last/read=",&
    !  firstprof,lastprof,myNumProfs
    ! Set numProfs if wanted
    if (present(numProfs)) numProfs=myNumProfs
    call trace_end (  'ReadL2GPData_MF_hdf', cond=.false. )

  end subroutine ReadL2GPData_MF_hdf

  !----------------------------------------  writeL2GPData_fileID  -----
  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  subroutine writeL2GPData_fileID(l2gp, l2FileHandle, swathName, hdfVersion, &
    & notUnlimited)

    ! Arguments

    integer, intent(in) :: l2FileHandle ! From swopen
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    integer, optional, intent(in) :: hdfVersion
    logical, optional, intent(in) :: notUnlimited
    ! Exectuable code

    ! Local
    integer :: myhdfVersion
    integer :: status
    type( MLSFile_T ) :: l2gpFile

    ! Executable code
    if (present(hdfVersion)) then
      myhdfVersion = hdfVersion
    else
      myhdfVersion = L2GPDEFAULT_HDFVERSION
    endif

    status = InitializeMLSFile(l2gpFile, type=l_swath, access=DFACC_RDWR, &
      & content='l2gp', name='unknown', hdfVersion=myhdfVersion)
    l2gpFile%FileID%f_id = l2FileHandle
    l2gpFile%stillOpen = .true.
    call WriteL2GPData( l2gp, l2gpFile, swathName, notUnlimited )

  end subroutine writeL2GPData_fileID

  ! --------------------------------------  OutputL2GP_createFile_MF  -----
  subroutine OutputL2GP_createFile_MF (l2gp, L2GPFile, &
    & swathName, nLevels, notUnlimited, compressTimes)

  use HDFEOS5, only: HE5S_Unlimited_F
  use MLSHDFEOS, only: MLS_SWDetach, &
    & MLS_Swcreate, MLS_DFldsetup, MLS_GFldsetup, MLS_SWDefdim
    ! Brief description of subroutine
    ! This subroutine sets up the structural definitions in an empty L2GP file.

    ! Arguments
    type(MLSFile_T)                :: L2GPFile
    type( L2GPData_T ), intent(inout) :: l2gp
    character (len=*), optional, intent(in) :: swathName ! Defaults to l2gp%swathName
    integer, optional, intent(in) :: nLevels
    logical, optional, intent(in) :: notUnlimited   !               as nTimes
    logical, optional, intent(in) :: compressTimes  ! don't store nTimesTotal

    ! Variables

    character (len=132) :: NAME   ! From l2gp%name
    character (len=32) :: MYDIM1, MYDIM12, MYDIM123

    ! THESE ARE HDF5 CHUNKS, _NOT_ MLS ALONG-TRACK PROCESSING CHUNKS 
    integer,dimension(7)::CHUNK_DIMS
    integer::CHUNK_RANK
    integer::CHUNKTIMES,CHUNKFREQS,CHUNKLEVELS
    integer :: hdfVersion

    integer :: SWID, STATUS
    logical :: myNotUnlimited
    logical :: mycompressTimes
    logical :: deebughere
    ! Executable
    deebughere = DEEBUG ! .or. .TRUE.
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    myNotUnlimited = .false.
    if ( present ( notUnlimited ) ) myNotUnlimited = notUnlimited
    mycompressTimes = .false.
    if ( present (compressTimes ) ) mycompressTimes = compressTimes
    hdfVersion = L2GPFile%hdfVersion
    
    ! Work out the chunking
    if ( myNotUnlimited ) then
      chunkTimes = max ( min ( MAXCHUNKTIMES, l2gp%nTimes ), 1 )
    else
      chunktimes = MAXCHUNKTIMES      ! was 1
    endif
    chunkfreqs = max ( l2gp%nFreqs, 1)
    if(present(nLevels))then
       chunklevels = nLevels
    else
      chunklevels = min(l2gp%nLevels, 500)     ! was .., 5)
      chunklevels = max(chunklevels, 1)
    endif
    
    ! Create the swath within the file

    if ( deebughere )  print *, 'About to sw_create ', TRIM(name)
    if ( deebughere )  print *, 'myNotUnlimited ', myNotUnlimited
    if ( deebughere )  print *, 'l2gp%MissingL2GP ', l2gp%MissingL2GP
    if ( deebughere )  print *, 'l2gp%MissingValue ', l2gp%MissingValue
    ! swid = mls_SWcreate(L2GPFile%FileID%f_id, trim(name), &
    swid = mls_SWcreate(L2GPFile, trim(name) )
    if ( swid == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to create swath ' // TRIM(name) &
        & // ' (maybe has the same name as another swath in this file?)', &
            & MLSFile=L2GPFile )
    end if

    ! Define dimensions

    if ( L2GPFile%hdfVersion == HDFVERSION_5 .and. .not. myNotUnlimited ) then
      ! Defining special "unlimited dimension called UNLIM
      ! print*,"Defined Unlim with size", HE5S_UNLIMITED_f
      status = mls_swdefdim(swid, UNLIM, HE5S_UNLIMITED_F, &
        & hdfVersion=hdfVersion)
    endif

    if ( myNotUnlimited ) then
      myDim1 = DIM_NAME1
      myDim12 = DIM_NAME12
      myDim123 = DIM_NAME123
    else
      myDim1 = MAX_DIML1
      myDim12 = MAX_DIML12
      myDim123 = MAX_DIML123
    endif
    if ( deebughere ) then
      print *, 'myDim1 ', myDim1
      print *, 'myDim12 ', myDim12
      print *, 'myDim123 ', myDim123
      print *, 'nTimes ', l2gp%nTimes
      print *, 'nTimesTotal ', l2gp%nTimesTotal
      print *, 'nLevels ', l2gp%nLevels
      print *, 'nFreqs ', l2gp%nFreqs
    endif
    if ( mycompressTimes ) then
      if ( deebughere ) call outputNamedValue ('nTimes',  max(l2gp%nTimes,1) )
      status = mls_swdefdim(swid, 'nTimes', max(l2gp%nTimes,1), &
        & hdfVersion=hdfVersion)
    else
      if ( deebughere ) call outputNamedValue ('nTimes',  max(l2gp%nTimesTotal,1) )
      status = mls_swdefdim(swid, 'nTimes', max(l2gp%nTimesTotal,1), &
        & hdfVersion=hdfVersion)
    endif
    status = mls_swdefdim(swid, 'nTimesTotal', max(l2gp%nTimesTotal,1), &
      & hdfVersion=hdfVersion)

    if ( l2gp%nLevels > 0 ) then
      status = mls_swdefdim(swid, 'nLevels', l2gp%nLevels, &
        & hdfVersion=hdfVersion)
    end if

    if ( l2gp%nFreqs > 0 ) then
      status = mls_swdefdim(swid, 'nFreqs', l2gp%nFreqs, &
        & hdfVersion=hdfVersion)
    end if

    ! Define horizontal geolocation fields using above dimensions

    chunk_rank=1
    chunk_dims=1
    chunk_dims(1)=CHUNKTIMES
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_gfldsetup(swid, 'Latitude', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'Longitude', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'Time', 'nTimes', MYDIM1, &
      & DFNT_FLOAT64, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, dFill=real(l2gp%MissingValue, r8))

    status = mls_gfldsetup(swid, 'LocalSolarTime', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'SolarZenithAngle', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'LineOfSightAngle', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'OrbitGeodeticAngle', 'nTimes', MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    status = mls_gfldsetup(swid, 'ChunkNumber', 'nTimes', MYDIM1, &
      & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, iFill=UndefinedIntegerValue)

    if ( l2gp%nLevels > 0 ) then

      chunk_rank=1
      chunk_dims(1)=CHUNKLEVELS
      status = mls_gfldsetup(swid, 'Pressure', 'nLevels', MAX_DIML, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)
    end if

    if ( l2gp%nFreqs > 0 ) then

      chunk_rank=1
      chunk_dims(1)=CHUNKFREQS
      status = mls_gfldsetup(swid, 'Frequency', 'nFreqs', MAX_DIML, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)
    end if

    ! Define data fields using above dimensions

    if ( (l2gp%nFreqs > 0) .and. (l2gp%nLevels > 0) ) then
       chunk_rank=3
       chunk_dims(1:3)=(/ CHUNKFREQS,CHUNKLEVELS,CHUNKTIMES /)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nFreqs,nLevels,nTimes', &
      & MYDIM123, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nFreqs,nLevels,nTimes', &
      & MYDIM123, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    else if ( l2gp%nLevels > 0 ) then
       chunk_rank=2
       chunk_dims(1:7)=(/ CHUNKLEVELS,CHUNKTIMES,37,38,39,47,49/)

      status = mls_dfldsetup(swid, 'L2gpValue', 'nLevels,nTimes', &
      & MYDIM12, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nLevels,nTimes', &
      & MYDIM12, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    else
       chunk_rank=1
       chunk_dims(1)=CHUNKTIMES

      status = mls_dfldsetup(swid, 'L2gpValue', 'nTimes', &
      & MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingL2GP)

      status = mls_dfldsetup(swid, 'L2gpPrecision', 'nTimes', &
      & MYDIM1, &
      & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    end if

    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_dfldsetup(swid, 'Status', 'nTimes', &
    & MYDIM1, &
    & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=hdfVersion, iFill=l2gp%MissingStatus)

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    if ( deebughere )  print *, 'chunk_dims ', chunk_dims(1:chunk_rank)
    status = mls_dfldsetup(swid, 'Quality', 'nTimes', &
    & MYDIM1, &
    & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)

    chunk_rank=1
    chunk_dims(1)=CHUNKTIMES
    status = mls_dfldsetup(swid, 'Convergence', 'nTimes', &
    & MYDIM1, &
    & DFNT_FLOAT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
    & hdfVersion=hdfVersion, rFill=l2gp%MissingValue)
    if ( AscDescModeIsField ) then
      chunk_rank=1
      chunk_dims(1)=CHUNKTIMES
      status = mls_dfldsetup(swid, 'AscDescMode ', 'nTimes', &
      & MYDIM1, &
      & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, iFill=int(l2gp%MissingValue) )
    endif

    if ( associated(l2gp%BinNumber) ) then
      chunk_rank=1
      chunk_dims(1)=CHUNKTIMES
      status = mls_dfldsetup(swid, 'BinNumber ', 'nTimes', &
      & MYDIM1, &
      & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, iFill=int(l2gp%MissingValue) )
    endif

    if ( associated(l2gp%MAF) ) then
      chunk_rank=1
      chunk_dims(1)=CHUNKTIMES
      status = mls_dfldsetup(swid, 'MAF ', 'nTimes', &
      & MYDIM1, &
      & DFNT_INT32, HDFE_NOMERGE, chunk_rank, chunk_dims, &
      & hdfVersion=hdfVersion, iFill=int(l2gp%MissingValue) )
    endif

    ! Detach from the HE5_SWath interface.This stores the swath info within the
    ! file and must be done before writing or reading data to or from the
    ! swath. (May be un-necessary for HDF5 -- test program works OK without.)
    ! 
    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Error, ModuleName, &
            & 'Failed to detach from swath interface after definition.', &
            & MLSFile=L2GPFile )
    end if

  !--------------------------------------
  end subroutine OutputL2GP_createFile_MF
  !--------------------------------------

  !-----------------------------------------  OutputL2GP_writeGeo_MF  -----
  subroutine OutputL2GP_writeGeo_MF (l2gp, L2GPFile, &
    & swathName,offset)

  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach, MLS_SWWrfld
    ! Brief description of subroutine
    ! This subroutine writes the geolocation fields to an L2GP output file.

    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name
    integer,intent(in),optional::offset

    ! Variables

    character (len=132) :: name ! Either swathName or l2gp%name
!     integer :: chunk
    integer :: status, swid,myOffset
    integer :: start(2), stride(2), edge(2)
    integer :: hdfVersion

    ! Begin
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    hdfVersion = L2GPFile%hdfVersion

    swid = mls_SWattach (L2GPFile, name)

    ! Write data to the fields
    edge = 0
    stride = 1
    start = myOffset ! Please do not set to zero
    edge(1) = l2gp%nTimes

    if (DEEBUG) then
      call output( 'Writing geolocations', advance='yes' )
      call outputNamedValue( 'stride, start, edge', (/ stride(1), start(1), edge(1) /) )
      call outputNamedValue( 'shape(Geod. Angle', shape(l2gp%geodAngle) )
    endif

    status = mls_SWwrfld(swid, 'Latitude', start, stride, edge, &
         real(l2gp%latitude), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'Longitude', start, stride, edge, &
         real(l2gp%longitude), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'Time', start, stride, edge, &
         l2gp%time, hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'LocalSolarTime', start, stride, edge, &
        real(l2gp%solarTime), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'SolarZenithAngle', start, stride, edge, &
         real(l2gp%solarZenith), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'LineOfSightAngle', start, stride, edge, &
         real(l2gp%losAngle), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'OrbitGeodeticAngle', start, stride, edge, &
         real(l2gp%geodAngle), hdfVersion=hdfVersion)

    ! call outputNamedValue ( 'Writing chunk number ', l2gp%chunkNumber )
    ! chunk = FindFirst( l2gp%chunkNumber /= -999 )
    ! call outputNamedValue ( 'index of non-Fill chunkNumber', chunk )
    status = mls_SWwrfld(swid, 'ChunkNumber', start, stride, edge, &
         l2gp%chunkNumber, hdfVersion=hdfVersion)

    if ( l2gp%nLevels > 0 ) then
       edge(1) = l2gp%nLevels
       start(1)=0 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Pressure', start(1:1), stride(1:1), edge(1:1), &
            real(l2gp%pressures), hdfVersion=hdfVersion)
    end if

    if ( l2gp%nFreqs > 0 ) then
       edge(1) = l2gp%nFreqs
       start(1)=0 ! needed because offset may have made this /=0
       status = mls_SWwrfld(swid, 'Frequency', start, stride, edge, &
            real(l2gp%frequency), hdfVersion=hdfVersion)
    end if

    ! Detach from the swath interface.  

    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface', MLSFile=L2GPFile )
    end if

  !------------------------------------
  end subroutine OutputL2GP_writeGeo_MF
  !------------------------------------

  !----------------------------------------  OutputL2GP_writeData_MF  -----
  subroutine OutputL2GP_writeData_MF(l2gp, L2GPFile, &
    & swathName,offset)

  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach, MLS_SWWrfld
    ! Brief description of subroutine
    ! This subroutine writes the data fields to an L2GP output file.
    ! For now, you have to write all of l2gp, but you can choose to write
    ! it at some offset into the file
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name
    integer,intent(in),optional::offset
    ! Parameters
    ! logical, parameter :: DEEBUG = .false.

    ! Variables

    character (len=132) :: name     ! Either swathName or l2gp%name

    integer :: status,myOffset
    integer :: start(3), stride(3), edge(3)
    integer :: swid
    integer :: hdfVersion
    real :: tFile ! How long have we been fooling with this file?
    

    ! Begin
    if (present(offset)) then
       myOffset=offset
    else
       myOffset=0
    endif

    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    hdfVersion = L2GPFile%hdfVersion

    start = 0
    stride = 1
    start(3)= myOffset ! Please do not set to zero
    edge(1) = l2gp%nFreqs
    edge(2) = l2gp%nLevels
    edge(3) = l2gp%nTimes
    call time_now ( tFile )
    if (DEEBUG) then
      call output( 'Writing data now ' // trim(name), advance='yes' )
      call outputNamedValue( 'stride, start, edge', (/ stride(3), start(3), edge(3) /) )
      call outputNamedValue( 'shape(Geod. Angle', shape(l2gp%geodAngle) )
    endif
    swid = mls_SWattach (L2GPFile, name)
    if ( DEEBUG ) then
      call sayTime( 'Attaching swath ' // trim(name), tFile, t2 )
      tFile = t2
    endif
    if ( l2gp%nFreqs > 0 ) then
       ! Value and Precision are 3-D fields
       if (DEEBUG) print *, 'start, stride, edge ', start, stride, edge
       ! What the devil???
       ! Why aren't we writing a 3-d array here??
       status = mls_SWwrfld(swid, 'L2gpValue', start, stride, edge, &
            & reshape(l2gp%l2gpValue, edge), &
            & hdfVersion=hdfVersion )
       status = mls_SWwrfld(swid, 'L2gpPrecision', start, stride, edge, &
            & reshape(l2gp%l2gpPrecision, edge), &
            & hdfVersion=hdfVersion )

       if ( DEEBUG ) then
         call sayTime( 'Witing 3d values and precision', tFile, t2 )
         tFile = t2
       endif
    else if ( l2gp%nLevels > 0 ) then
       ! Value and Precision are 2-D fields
      
       if (DEEBUG) print *, 'start, stride, edge: ', start(2:3), stride(2:3), edge(2:3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(2:3), stride(2:3), &
            edge(2:3), l2gp%l2gpValue(1,:,:), hdfVersion=hdfVersion)
       if ( DEEBUG ) then
         call sayTime( 'Witing 2d values', tFile, t2 )
         tFile = t2
       endif
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(2:3), stride(2:3), &
            edge(2:3), l2gp%l2gpPrecision(1,:,:), hdfVersion=hdfVersion)
       if ( DEEBUG ) then
         call sayTime( 'Witing 2d precision', tFile, t2 )
         tFile = t2
       endif
    else

       ! Value and Precision are 1-D fields
       if (DEEBUG) print *, 'start, stride, edge: ', start(3), stride(3), edge(3)
       status = mls_SWwrfld( swid, 'L2gpValue', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpValue(1,1,:) ), hdfVersion=hdfVersion)
       if ( DEEBUG ) then
         call sayTime( 'Witing 1d values', tFile, t2 )
         tFile = t2
       endif
       status = mls_SWwrfld( swid, 'L2gpPrecision', start(3:3), stride(3:3), edge(3:3), &
            real(l2gp%l2gpPrecision(1,1,:) ), hdfVersion=hdfVersion)
       if ( DEEBUG ) then
         call sayTime( 'Witing 1d precision', tFile, t2 )
         tFile = t2
       endif
    end if

    ! 1-D status & quality fields

    status = mls_swwrfld(swid, 'Status', start(3:3), stride(3:3), edge(3:3), &
       &   l2gp%status, hdfVersion=hdfVersion, dontfail=.true.)

    status = mls_SWwrfld(swid, 'Quality', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%quality), hdfVersion=hdfVersion)

    status = mls_SWwrfld(swid, 'Convergence', start(3:3), stride(3:3), edge(3:3), &
         real(l2gp%convergence), hdfVersion=hdfVersion)

    if ( AscDescModeIsField ) then
      status = mls_SWwrfld(swid, 'AscDescMode ', start(3:3), stride(3:3), edge(3:3), &
         l2gp%AscDescMode, hdfVersion=hdfVersion)
    endif

    if ( associated(l2gp%BinNumber) ) then
      status = mls_SWwrfld(swid, 'BinNumber ', start(3:3), stride(3:3), edge(3:3), &
         l2gp%BinNumber, hdfVersion=hdfVersion)
    endif

    if ( associated(l2gp%MAF) ) then
      status = mls_SWwrfld(swid, 'MAF ', start(3:3), stride(3:3), edge(3:3), &
         l2gp%MAF, hdfVersion=hdfVersion)
    endif

    if ( DEEBUG ) then
      call sayTime( 'Witing 1d status, Quality, Asc/Desc', tFile, t2 )
      tFile = t2
    endif
    !     Detach from the swath interface.
    status = mls_SWdetach(swid, hdfVersion=hdfVersion)
    if(DEEBUG) print *, 'Detached from swid ', swid
    if(DEEBUG) print *, 'file handle ', L2GPFile%FileID%f_id
    if(DEEBUG) print *, 'status ', status
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach from swath interface', MLSFile=L2GPFile )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_writeData_MF
  !-------------------------------------

  !----------------------------------------  OutputL2GP_attributes_MF  -----
  subroutine OutputL2GP_attributes_MF(l2gp, L2GPFile, swathName)

  use HDFEOS5, only: HE5T_Native_Int, HE5T_Native_Real, HE5T_Native_Double, &
    & MLS_Chartype
  use HE5_SWAPI, only: HE5_SWWrattr, HE5_SWWrlattr
  use MLSHDFEOS, only: MLS_ISGlatt, &
    & MLS_SWAttach, MLS_SWDetach, MLS_SWWrattr, MLS_SWWrlattr
  use PCFHDR, only:  HE5_Writeglobalattr
    ! Brief description of subroutine
    ! This subroutine writes the attributes for an l2gp
    ! These include
    ! Global attributes which are file level
    ! Swath level
    ! Geolocation field attributes
    ! Data field attributes
    ! Arguments

    type( L2GPData_T ), intent(inout) :: l2gp
    type(MLSFile_T)                :: L2GPFile
    character (len=*), intent(in), optional :: swathName ! Defaults->l2gp%name

    ! Parameters
    character (len=*), parameter :: NOUNITS = 'NoUnits'

    ! The following pair of string list encode the Units attribute
    ! corresponding to each Title attribute; e.g., the Units for Latitude is deg
    character (len=*), parameter :: GeolocationTitles = &
      & 'Latitude,Longitude,Time,LocalSolarTime,SolarZenithAngle,' // &
      & 'LineOfSightAngle,OrbitGeodeticAngle,ChunkNumber,Pressure,Frequency'
    character (len=*), parameter :: GeolocationUnits = &
      & 'deg,deg,s,h,deg,' // &
      & 'deg,deg,NoUnits,hPa,GHz'
    character (len=*), parameter :: GeoUniqueFieldDefinition = &
      & 'HMT,HMT,AS,HMT,HMT,' // &
      & 'M,M,M,AS,M'   ! These are abbreviated values
    character (len=*), parameter :: UniqueFieldDefKeys = &
      & 'HM,HMT,MT,AS,M'
    character (len=*), parameter :: UniqueFieldDefValues = &
      & 'HIRDLS-MLS-Shared,HIRDLS-MLS-TES-Shared,MLS-TES-Shared,' // &
      & 'Aura-Shared,MLS-Specific'  ! Expanded values
    ! The following associate UniqueFieldDefs with species names
    character (len=*), parameter :: Species = &
      & 'Temperature,BrO,CH3CN,CO,ClO,GPH,HCl,HCN,H2O,H2O2,' // &
      & 'HNO3,HOCl,HO2,N2,N2O,OH,O2,O3,RHI,SO2'
    character (len=*), parameter :: SpUniqueFieldDefinition = &
      & 'HMT,M,M,MT,M,M,M,M,HMT,M,' // &
      & 'HMT,M,M,M,HM,M,M,HMT,M,M'   ! These are abbreviated values
    ! logical, parameter :: DEEBUG = .true.

    ! Variables
    character (len=132) :: name     ! Either swathName or l2gp%name
    character(len=CHARATTRLEN) :: abbr_uniq_fdef
    character(len=CHARATTRLEN) :: expnd_uniq_fdef
    integer :: field
    character(len=CHARATTRLEN) :: field_name
    logical :: isColumnAmt
    logical :: isTPPressure
    integer :: rgp_type
    character(len=CHARATTRLEN) :: species_name ! Always lower case
    integer :: status
    integer :: swid
    character(len=CHARATTRLEN) :: temp_name
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theTitles
    character(len=CHARATTRLEN), dimension(NumGeolocFields) :: theUnits
    character(len=CHARATTRLEN) :: units_name
    ! Begin
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    if ( DEEBUG ) print *, 'About to wr attrs to: ', trim(name)
    if ( rgp == r4 ) then
      rgp_type = HE5T_NATIVE_REAL
    elseif ( rgp == r8 ) then
      rgp_type = HE5T_NATIVE_DOUBLE
    else
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'Attributes have unrecognized numeric data type; should be r4 or r8', &
        & MLSFile=L2GPFile)
    endif
    call List2Array(GeolocationTitles, theTitles, .true.)
    call List2Array(GeolocationUnits, theUnits, .true.)

    ! - -   G l o b a l   A t t r i b u t e s   - -
    if ( WRITEMASTERSFILEATTRIBUTES .and. &
      & .not. MLS_ISGLATT(L2GPFile%fileID%f_id, 'StartUTC') ) &
      & call he5_writeglobalattr(L2GPFile%fileID%f_id)

    swid = mls_SWattach (L2GPFile, name)
    
    !   - -   S w a t h   A t t r i b u t e s   - -
    status = he5_swwrattr(swid, trim(l2gp%verticalCoordinate), &
      & rgp_type, hsize(size(l2gp%pressures)), &
      & l2gp%pressures)
    field_name = l2gp%verticalCoordinate ! 'Pressure'
    status = mls_swwrattr(swid, 'VerticalCoordinate', MLS_CHARTYPE, 1, &
      & field_name)
    if ( DeeBug ) print *, 'Missing L2GP: ', real(l2gp%MissingL2GP)
    if ( DeeBug ) print *, 'Missing value: ', real(l2gp%MissingValue)
    if ( SWATHLEVELMISSINGVALUE ) &
      & status = he5_swwrattr(swid, 'MissingValue', rgp_type, hsize(1), &
      & (/ real(l2gp%MissingL2GP, rgp) /) )
    
    !   - -   G e o l o c a t i o n   A t t r i b u t e s   - -
    if ( DEEBUG ) print *, 'About to wr loc attrs to: ', trim(name)
    do field=1, NumGeolocFields
      ! Take care not to write attributes to "missing fields"
      if ( trim(theTitles(field)) == 'Frequency' &
        & .and. l2gp%nFreqs < 1 ) then
        field_name = ''
      elseif ( trim(theTitles(field)) == 'Pressure' &
        & .and. l2gp%nLevels < 1 ) then
        field_name = ''
      else
        field_name = theTitles(field)
        if ( trim(theTitles(field)) == 'Pressure' ) &
          & field_name = l2gp%verticalCoordinate
        call GetHashElement (GeolocationTitles, &
          & GeoUniqueFieldDefinition, trim(theTitles(field)), &
          & abbr_uniq_fdef, .false.)
        call GetHashElement (UniqueFieldDefKeys, &
          & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
          & expnd_uniq_fdef, .false.)
        if ( DEEBUG ) print *, 'Field Title ', trim(theTitles(field))
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Title', &
          & MLS_CHARTYPE, 1, theTitles(field))
        if ( DEEBUG ) print *, 'Units ', trim(theUnits(field))
        status = mls_swwrlattr(swid, trim(theTitles(field)), 'Units', &
          & MLS_CHARTYPE, 1, theUnits(field))

        if ( trim(theTitles(field)) == 'Time' ) then
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & HE5T_NATIVE_DOUBLE, hsize(1), (/ real(l2gp%MissingValue, r8) /) )
        elseif ( trim(theTitles(field)) == 'ChunkNumber' ) then
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & HE5T_NATIVE_INT, hsize(1), (/ UndefinedIntegerValue /) )
        else
          status = he5_swwrlattr(swid, trim(theTitles(field)), &
            & 'MissingValue', &
            & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
        endif
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to write local attribute to " // trim(theTitles(field)) )
        if ( DEEBUG ) print *, 'Uniquefielddef ', trim(expnd_uniq_fdef)
        status = mls_swwrlattr(swid, trim(theTitles(field)), &
          & 'UniqueFieldDefinition', &
          & MLS_CHARTYPE, 1, expnd_uniq_fdef)
        ! print *, 'status : ', status
        if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, ModuleName, &
         & "Unable to write local attribute to " // trim(theTitles(field)) )
      endif
    enddo
    if ( DEEBUG ) print *, 'Data'
    !   - -   D a t a   A t t r i b u t e s   - -
    ! call GetQuantityAttributes ( l2gp%quantityType, &
    !  & units_name, expnd_uniq_fdef)
    field_name = Name
    species_name = lowercase(name)
    ! The following special cases are handled by crude, despicable hacks
    ! It would be better to check on the quantity type, but paw hasn't
    ! succeeded in getting that to work properly and reliably
    isColumnAmt = ( index(species_name, 'column') > 0 )
    isTPPressure = ( index(species_name, 'tpp') > 0 )
    if ( isColumnAmt ) then
      ! This next failed sometimes
      temp_name = species_name
      call ExtractSubString(Name, species_name, 'column', 'wmo')
      if ( species_name == ' ' ) species_name=temp_name
      ! So we'll try another tack:
      temp_name = species_name
      call ReplaceSubString ( temp_name, species_name, 'column', '' )
      if ( species_name == ' ' ) species_name=temp_name
      temp_name = species_name
      call GetStringElement( temp_name, species_name, &
        & 1, countEmpty=.true., inseparator='-' )
      if ( species_name == ' ' ) species_name=temp_name
    else
      temp_name = species_name
      call GetStringElement( temp_name, species_name, &
        & 1, countEmpty=.true., inseparator='-' )
      if ( species_name == ' ' ) species_name=temp_name
    endif
    ! Hopefully by now we've turned species_name into one of the species
    if ( DEEBUG ) &
      & print *, 'full name: ', trim(name), ' Species: ', trim(species_name)
    call GetHashElement (lowercase(Species), &
      & SpUniqueFieldDefinition, trim(species_name), &
      & abbr_uniq_fdef, .false.)
    call GetHashElement (UniqueFieldDefKeys, &
      & UniqueFieldDefValues, trim(abbr_uniq_fdef), &
      & expnd_uniq_fdef, .false.)
    if ( expnd_uniq_fdef == '' .or. expnd_uniq_fdef == ',' ) &
      & expnd_uniq_fdef = 'MLS-Specific'
    select case (trim(lowercase(species_name)))
    case ('temperature')
      units_name = 'K'
    case ('gph')
      units_name = 'm'
    case ('geoheight')
      units_name = 'km'
    case ('rhi')
      units_name = '%rhi'
    case ('iwc', 'iwp')
      units_name = 'g/m^3'
    case default
      units_name = 'vmr'
    end select
    if ( isColumnAmt ) then
      ! units_name = 'DU'
      if ( DEEBUG ) &
        & call dump( .true., col_species_keys, col_species_hash, &
        & 'column species units' )
      call GetHashElement (col_species_keys, &
      & col_species_hash, trim(lowercase(species_name)), &
      & units_name, .true.)
      if ( DEEBUG ) &
      & print *, 'units_name: ', trim(units_name)
      if ( units_name == ',' ) units_name = 'molcm2'
    endif
    if ( isTPPressure ) units_name = 'hPa'
    if ( DEEBUG ) print *, 'Title ', trim(field_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'Title', &
      & MLS_CHARTYPE, 1, field_name)
    if ( DEEBUG ) print *, 'Units ', trim(units_name)
    status = mls_swwrlattr(swid, 'L2gpValue', 'Units', &
      & MLS_CHARTYPE, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpValue', 'MissingValue', &
      & rgp_type, hsize(1), (/ real(l2gp%MissingL2GP, rgp) /) )
    if ( DEEBUG ) print *, 'Title ', trim(expnd_uniq_fdef)
    status = mls_swwrlattr(swid, 'L2gpValue', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, expnd_uniq_fdef)
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Precision')
    status = mls_swwrlattr(swid, 'L2gpPrecision', 'Units', &
      & MLS_CHARTYPE, 1, units_name)
    status = he5_swwrlattr(swid, 'L2gpPrecision', 'MissingValue', &
      & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'L2gpPrecision', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, expnd_uniq_fdef)

    ! ('Status' data field newly written)
    status = mls_swwrlattr(swid, 'Status', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Status')
    status = mls_swwrlattr(swid, 'Status', 'Units', &
      & MLS_CHARTYPE, 1, NOUNITS)
    status = he5_swwrlattr(swid, 'Status', 'MissingValue', &
      & HE5T_NATIVE_INT, hsize(1), (/ l2gp%MissingStatus /) )
    status = mls_swwrlattr(swid, 'Status', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, 'MLS-Specific')
    
    status = mls_swwrlattr(swid, 'Quality', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Quality')
    status = mls_swwrlattr(swid, 'Quality', 'Units', &
      & MLS_CHARTYPE, 1, NOUNITS)
    status = he5_swwrlattr(swid, 'Quality', 'MissingValue', &
      & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'Quality', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, 'MLS-Specific')
    
    status = mls_swwrlattr(swid, 'Convergence', 'Title', &
      & MLS_CHARTYPE, 1, trim(field_name)//'Convergence')
    status = mls_swwrlattr(swid, 'Convergence', 'Units', &
      & MLS_CHARTYPE, 1, NOUNITS)
    status = he5_swwrlattr(swid, 'Convergence', 'MissingValue', &
      & rgp_type, hsize(1), (/ real(l2gp%MissingValue, rgp) /) )
    status = mls_swwrlattr(swid, 'Convergence', &
      & 'UniqueFieldDefinition', &
      & MLS_CHARTYPE, 1, 'MLS-Specific')
    if ( AscDescModeIsField ) then
      status = mls_swwrlattr(swid, 'AscDescMode', 'Title', &
        & MLS_CHARTYPE, 1, trim(field_name)//'AscDescMode ')
      status = mls_swwrlattr(swid, 'AscDescMode', 'Units', &
        & MLS_CHARTYPE, 1, NOUNITS)
      status = he5_swwrlattr(swid, 'AscDescMode', 'MissingValue', &
        & HE5T_NATIVE_INT, hsize(1), (/ 0 /) )
      status = mls_swwrlattr(swid, 'AscDescMode', &
        & 'UniqueFieldDefinition', &
        & MLS_CHARTYPE, 1, 'MLS-Specific')
      endif
    
    status = mls_SWdetach(swid, hdfVersion=HDFVERSION_5)
    if ( status == -1 ) then
       call MLSMessage ( MLSMSG_Warning, ModuleName, &
            & 'Failed to detach  from swath interface', MLSFile=L2GPFile )
    end if


  !-------------------------------------
  end subroutine OutputL2GP_attributes_MF
  !-------------------------------------

  !-----------------------------------------  FilterL2GP  -----
  subroutine FilterL2GP ( L2GP, fields, options )

    ! This routine filters an  l2gp by setting its status to crashed
    ! wherever the gelolocations contain fillvalues
    
    ! Thew optional arg fields will dtermine which fields this applies to:
    ! value                  fields
    ! -----                  ------
    ! '*' (or missing)    all
    ! 'latitude,time'     latitude, time
    ! 'l2gpvalue'         l2gpvalue
    ! 'value'             l2gpvalue, l2gpPrecision, status, quality
    ! 'geolocation'       geolocations (all except 'value')

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP
    character(len=*), optional, intent(in) :: fields
    character(len=*), optional, intent(in) :: options ! E.g., '-v'
    
    ! Internal variables
    character(len=128) :: myFields
    character (len=8) :: myOptions

    ! Executable code
    myFields = GEO_FIELDS ! '*' ! would mean values or geolocations
    if ( present(fields) ) myFields = lowercase(fields)
    myOptions = ' '
    if ( present(options) ) myOptions = options
    ! verbose = ( index(myOptions, 'v') > 0 )
    select case (trim(myFields))
    case ('value')
      myFields = lowercase(DATA_FIELDS) ! 'l2gpvalue, l2gpPrecision, status, quality'
    case ('geolocation')
      myFields = lowercase(GEO_FIELDS) ! // ',pressure,time,frequency,chunknumber'
    case default
      ! Not a special case
    end select

    if ( SwitchDetail(myFields, 'latitude', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%latitude) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'longitude', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%longitude) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'solartime', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%solartime) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'solarzenith', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%solarzenith) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'losangle', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%losangle) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'geodangle', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%geodangle) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'time', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%time) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'l2gpvalue', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%l2gpvalue(1,1,:)) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif
    if ( SwitchDetail(myFields, 'l2gpprecision', '-wfc') > -1 ) then
      where ( isFillValue(l2gp%l2gpprecision(1,1,:)) )
         l2gp%status = l2gp%MissingStatus  !DANGERWILLROBINSON
      endwhere
    endif

  end subroutine FilterL2GP
    
  !-----------------------------------------  RepairL2GP_L2GP  -----
  subroutine RepairL2GP_L2GP ( L2GP1, L2GP2, fields, options )

    ! This routine repairs l2gp1 using values from l2gp2
    ! wherever the first has fillvalues
    
    ! Thew optional arg fields will dtermine which fields this applies to:
    ! value                  fields
    ! -----                  ------
    ! '*' (or missing)    all
    ! 'latitude,time'     latitude, time
    ! 'l2gpvalue'         l2gpvalue
    ! 'value'             l2gpvalue, l2gpPrecision, status, quality
    ! 'geolocation'       geolocations (all except 'value')

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP1
    type (L2GPData_T), intent(in)    :: L2GP2
    character(len=*), optional, intent(in) :: fields
    character(len=*), optional, intent(in) :: options ! E.g., '-v'
    
    ! Internal variables
    character(len=128) :: myFields
    character (len=8) :: myOptions

    ! Executable code
    myFields = '*'
    if ( present(fields) ) myFields = lowercase(fields)
    myOptions = ' '
    if ( present(options) ) myOptions = options
    ! verbose = ( index(myOptions, 'v') > 0 )
    select case (trim(myFields))
    case ('value')
      myFields = lowercase(DATA_FIELDS) ! 'l2gpvalue, l2gpPrecision, status, quality'
    case ('geolocation')
      myFields = lowercase(GEO_FIELDS) ! // ',pressure,time,frequency,chunknumber'
    case default
      ! Not a special case
    end select

    if ( SwitchDetail(myFields, 'pressure', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%pressures, l2gp1%MissingValue, l2gp2%pressures )
    endif
    if ( SwitchDetail(myFields, 'latitude', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%latitude, l2gp1%MissingValue, l2gp2%latitude )
    endif
    if ( SwitchDetail(myFields, 'longitude', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%longitude, l2gp1%MissingValue, l2gp2%longitude )
    endif
    if ( SwitchDetail(myFields, 'solartime', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%solartime, l2gp1%MissingValue, l2gp2%solartime )
    endif
    if ( SwitchDetail(myFields, 'solarzenith', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%solarzenith, l2gp1%MissingValue, l2gp2%solarzenith )
    endif
    if ( SwitchDetail(myFields, 'LineOfSightAngle', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%losangle, l2gp1%MissingValue, l2gp2%losangle )
    endif
    if ( SwitchDetail(myFields, 'OrbitGeodeticAngle', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%geodangle, l2gp1%MissingValue, l2gp2%geodangle )
    endif
    if ( SwitchDetail(myFields, 'time', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%time, real(l2gp1%MissingValue, r8), &
        & l2gp2%time )
    endif
    if ( SwitchDetail(myFields, 'chunknumber', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%chunknumber, UndefinedIntegerValue, &
        & l2gp2%chunknumber )
    endif
    if ( SwitchDetail(myFields, 'frequency', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%frequency, l2gp1%MissingValue, l2gp2%frequency )
    endif
    if ( SwitchDetail(myFields, 'l2gpvalue', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%l2gpvalue, l2gp1%MissingL2GP, l2gp2%l2gpvalue )
    endif
    if ( SwitchDetail(myFields, 'l2gpprecision', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%l2gpprecision, l2gp1%MissingValue, l2gp2%l2gpprecision )
    endif
    if ( SwitchDetail(myFields, 'status', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%status, l2gp1%MissingStatus, l2gp2%status )
    endif
    if ( SwitchDetail(myFields, 'quality', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%quality, l2gp1%MissingValue, l2gp2%quality )
    endif
    if ( SwitchDetail(myFields, 'convergence', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp1%convergence, l2gp1%MissingValue, l2gp2%convergence )
    endif

  end subroutine RepairL2GP_L2GP

  !-----------------------------------------  RepairL2GP_HGrid  -----
  subroutine RepairL2GP_HGrid ( L2GP, HGrid, fields, offset, options )
    use HGridsDatabase, only: HGrid_T
    ! This routine repairs l2gp1 using values from HGrid
    ! wherever the first has fillvalues
    
    ! The optional arg fields will determine which fields this applies to:
    ! value                  fields
    ! -----                  ------
    ! '*' (or missing)    all
    ! 'latitude,time'     latitude, time
    ! 'geolocation'       geolocations (all)
    
    ! If the optional arg options contains
    ! value       effect
    ! -----       ------
    !   v        verbose
    !   d        DEEBUG
    !   c        forcibly repair any geolocations whose chunk number is -999

    ! Dummy arguments
    type (L2GPData_T), intent(inout) :: L2GP
    type (HGrid_T), intent(in)       :: HGrid
    character(len=*), optional, intent(in) :: fields
    integer, optional, intent(in) :: offset ! If HGrid profs offset from l2gp    
    character (len=*), optional, intent(in) :: options ! E.g., '-v'
    ! Internal variables
    ! logical, parameter :: DEEBUG = .false.
    logical            :: force
    integer            :: HGp1 ! Effective starting HGrid profile
    integer            :: HGpn ! Effective ending HGrid profile
    character(len=128) :: myFields
    character (len=8)  :: myOptions
    logical            :: verbose
    ! Executable code
    myFields = '*'
    if ( present(fields) ) myFields = lowercase(fields)
    myOptions = ' '
    if ( present(options) ) myOptions = options
    verbose = ( index(myOptions, 'v') > 0 ) .or. DEEBUG
    if ( present(offset) ) then
      HGp1 = offset
    else
      HGp1 = 1 + HGrid%noProfsLowerOverlap ! 2  ! No longer assume we're not offset; was 1
    endif
    force = ( index(myOptions, 'c') > 0 )
    ! Check that we've calculated array sizes correectly
    ! In particular, useable profiles matched HGrid and l2gp
    if ( size(hgrid%geodLat(1,HGp1:)) < size(l2gp%latitude) ) then
      call output('shape l2gp%latitude: ')
      call output(shape(l2gp%latitude), advance='yes')
      call output('shape HGrid%geodLat: ')
      call output(shape(HGrid%geodLat), advance='yes')
      call output('HGp1: ')
      call output(HGp1, advance='yes')
      call dump(hgrid%geodLat(1,HGp1:), 'hgrid%geodLat(1,HGp1:)')
      call dump(l2gp%latitude, 'l2gp%latitude')
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & 'HGrid and l2gp size mismatch: HGrid too small')
    elseif ( size(hgrid%geodLat(1,HGp1:)) > size(l2gp%latitude) ) then
      call MLSMessage ( MLSMSG_Warning, ModuleName, & 
        & 'HGrid and l2gp size mismatch: HGrid needed to be trimmed; done')
      call output('shape l2gp%latitude: ')
      call output(shape(l2gp%latitude), advance='yes')
      call output('shape HGrid%geodLat: ')
      call output(shape(HGrid%geodLat), advance='yes')
      verbose = .true.
    endif
    HGpn = size(l2gp%latitude) + HGp1 - 1
    if ( verbose ) then
      call output('1st HGrid angle: ')
      call output(hgrid%phi(1,HGp1), advance='yes')
      call output('1st l2gp angle: ')
      call output(l2gp%geodAngle(1), advance='yes')
      call output('last HGrid angle: ')
      call output(hgrid%phi(1,HGpn), advance='yes')
      call output('last l2gp angle: ')
      call output( l2gp%geodAngle( size(l2gp%geodAngle) ), advance='yes' )
      call output('HGp1: ')
      call output(HGp1, advance='yes')
      call output('HGpn: ')
      call output(HGpn, advance='yes')
      call output('shape l2gp%latitude: ')
      call output(shape(l2gp%latitude), advance='yes')
      call output('shape HGrid%phi(1,HGp1:HGpn): ')
      call output(shape(HGrid%phi(1,HGp1:HGpn)), advance='yes')
    endif
    
    if ( force ) then
      where ( l2gp%chunkNumber == -999 )
        l2gp%latitude      = l2gp%MissingValue
        l2gp%longitude     = l2gp%MissingValue
        l2gp%solartime     = l2gp%MissingValue
        l2gp%solarzenith   = l2gp%MissingValue
        l2gp%losangle      = l2gp%MissingValue
        l2gp%geodAngle     = l2gp%MissingValue
        l2gp%time          = l2gp%MissingValue
      endwhere
    endif

    select case (trim(myFields))
    case ('geolocation')
      myFields = lowercase(GEO_FIELDS) ! // ',pressure,time,frequency,chunknumber'
    case default
      ! Not a special case
    end select

    if ( DEEBUG ) print *, 'shape l2gp%latitude: ', shape(l2gp%latitude)
    if ( DEEBUG ) print *, 'shape HGrid%geodLat: ', shape(HGrid%geodLat)
    
    if ( SwitchDetail(myFields, 'latitude', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%latitude, l2gp%MissingValue, &
        & real(hgrid%geodLat(1,HGp1:HGpn), rgp) )
    endif
    if ( DEEBUG ) print *, 'shape l2gp%longitude: ', shape(l2gp%longitude)
    if ( DEEBUG ) print *, 'shape HGrid%lon: ', shape(HGrid%lon)
    if ( SwitchDetail(myFields, 'longitude', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%longitude, l2gp%MissingValue, &
        & real(hgrid%lon(1,HGp1:HGpn), rgp) )
    endif
    if ( SwitchDetail(myFields, 'solartime', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%solartime, l2gp%MissingValue, &
        & real(hgrid%solartime(1,HGp1:HGpn), rgp) )
    endif
    if ( SwitchDetail(myFields, 'solarzenith', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%solarzenith, l2gp%MissingValue, &
        & real(hgrid%solarzenith(1,HGp1:HGpn), rgp) )
    endif
    if ( SwitchDetail(myFields, 'LineOfSightAngle', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%losangle, l2gp%MissingValue, &
        & real(hgrid%losangle(1,HGp1:HGpn), rgp) )
    endif
    if ( SwitchDetail(myFields, 'OrbitGeodeticAngle', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%geodAngle, l2gp%MissingValue, &
        & real(hgrid%phi(1,HGp1:HGpn), rgp) )
    endif
    if ( SwitchDetail(myFields, 'time', '-wfc') > -1 ) then
      call ReplaceFillValues ( l2gp%time, real(l2gp%MissingValue, r8), &
        & hgrid%time(1,HGp1:HGpn) )
    endif

  end subroutine RepairL2GP_HGrid

  !----------------------------------------  SetL2GP_aliases_MF  -----
  subroutine SetL2GP_aliases_MF(l2gp, L2GPFile, swathName)

  use HDFEOS5, only: HE5_SWSetalias
  use MLSHDFEOS, only: MLS_SWAttach, MLS_SWDetach
  use SDPToolkit, only: PGS_S_Success
    ! Arguments
    type(MLSFile_T)                :: L2GPFile
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    
    ! Brief description of subroutine
    ! This subroutine creates an alias for each of the two data fields
    ! TYPE2FIELDNAME and TYPE2PRECISIONNAME
    ! Local variables
    character (len=132) :: name     ! Either swathName or l2gp%name
    character(len=*), parameter :: TYPE2FIELDNAME = 'L2gpValue'
    character(len=*), parameter :: TYPE2PRECISIONNAME = 'L2gpPrecision'
    integer :: returnStatus
    integer :: sw_id
    ! Executable
    if (present(swathName)) then
       name=swathName
    else
       name=l2gp%name
    endif
    sw_id = mls_swattach(L2GPFile, trim(name))
    if ( sw_id < 1 ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in attaching swath for setting alias.", MLSFile=L2GPFile )
    end if
    returnStatus = he5_SWsetalias(sw_id, TYPE2FIELDNAME, trim(name))
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in setting alias from " // TYPE2FIELDNAME // &
          & ' to ' // trim(name), MLSFile=L2GPFile )
    end if
    returnStatus = he5_SWsetalias(sw_id, TYPE2PRECISIONNAME, &
     & trim(name) // 'Precision')
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in setting alias from " // TYPE2PRECISIONNAME // &
          & ' to ' // trim(name) // 'Precision', MLSFile=L2GPFile )
    end if
    returnStatus = mls_SWdetach(sw_id, hdfVersion=HDFVERSION_5)
    if ( returnStatus /= PGS_S_SUCCESS ) then 
      call MLSMessage ( MLSMSG_Error, ModuleName, & 
        & "Error in detaching swath for setting alias.", MLSFile=L2GPFile )
    end if
  !-------------------------------------
  end subroutine SetL2GP_aliases_MF
  !-------------------------------------

  !----------------------------------------  writeL2GPData_MLSFile  -----
  ! This subroutine is an amalgamation of the last three
  ! Should be renamed CreateAndWriteL2GPData
  subroutine writeL2GPData_MLSFile(l2gp, L2GPFile, swathName, &
    & notUnlimited)

    ! Arguments

    type(MLSFile_T)                :: L2GPFile
    type (L2GPData_T), intent(INOUT) :: l2gp
    character (len=*), optional, intent(in) ::swathName!default->l2gp%swathName
    logical, optional, intent(in) :: notUnlimited
    ! Local
    logical :: alreadyOpen
    integer :: Me = -1                  ! String index for trace cacheing
    integer :: status
    ! Executable code

    call trace_begin ( me, 'writeL2GPData_MLSFile', cond=.false. )
    status = 0
    alreadyOpen = L2GPFile%stillOpen
    if ( .not. alreadyOpen ) then
      call mls_openFile(L2GPFile, Status)
      if ( Status /= 0 ) &
        call MLSMessage(MLSMSG_Error, ModuleName, &
        & 'Unable to open l2gp file', MLSFile=L2GPFile)
    endif
    if ( L2GPFile%access == DFACC_RDONLY )  &
      & call MLSMessage(MLSMSG_Error, ModuleName, &
      & 'l2gp file is rdonly', MLSFile=L2GPFile)
    call OutputL2GP_createFile_MF (l2gp, L2GPFile, &
      & swathName, notUnlimited=notUnlimited)
    call OutputL2GP_writeGeo_MF (l2gp, L2GPFile, &
      & swathName)
    call OutputL2GP_writeData_MF (l2gp, L2GPFile, &
      & swathName)
    if (L2GPFile%hdfVersion == HDFVERSION_5) then
      if ( DEEBUG ) print *, 'Outputting attributes'
      call OutputL2GP_attributes_MF (l2gp, L2GPFile, swathName)
      if ( DEEBUG ) print *, 'Setting aliases'
      call SetL2GP_aliases_MF (l2gp, L2GPFile, swathName)
    endif

    if ( .not. alreadyOpen )  call mls_closeFile(L2GPFile, Status)
    L2GPFile%errorCode = status
    L2GPFile%lastOperation = 'write'
    call trace_end ( 'writeL2GPData_MLSFile', cond=.false. )
  end subroutine writeL2GPData_MLSFile
  
 ! The next two functions ignore leap seconds
  elemental function hoursInDayToTime( hid, l2gp ) result( time )
    ! Given rgp hid, return r8 time as tai(s)
    ! Args:
    real(rgp), intent(in)            :: hid ! hours in day
    type (L2GPData_T), intent(in)   :: l2gp
    real(r8)                        :: time
    ! Executable
    time = 3600._r8*hid + l2gp%time(1)
  end function hoursInDayToTime

  elemental function TimeToHoursInDay( time, l2gp ) result( hid )
    ! Given r8 time as tai(s), return rgp hid
    use Dates_Module, only: TAI93S2HID
    ! Args:
    real(r8), intent(in)                      :: time
    type (L2GPData_T), optional, intent(in)   :: l2gp
    real(rgp)                                 :: hid ! hours in day
    ! Local variables
    ! Executable
    if ( present(l2gp) ) then
      hid = (time - l2gp%time(1)) / 3600
    else
      hid = tai93s2hid( time )
    endif
  end function TimeToHoursInDay

!------------------------- SayTimeHere ---------------------
  subroutine SayTimeHere ( What, startTime )
    character(len=*), intent(in) :: What
    real, intent(in)             :: startTime
    call time_now ( t2 )
    call output ( "Timing for " // what // " = " )
    call output ( dble(t2 - startTime), advance = 'yes' )
  end subroutine SayTimeHere

!=============================================================================
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

!=============================================================================
end module L2GPData
!=============================================================================

!
! $Log$
! Revision 2.249  2021/05/27 23:39:09  pwagner
! Added 'F' option to print diffs for compatible datasets
!
! Revision 2.248  2021/04/15 22:43:22  pwagner
! Now uses MLS_HyperStart
!
! Revision 2.247  2020/07/17 16:16:14  pwagner
! Avoid reading hdf5 MissingValue attribute from hdf4 files
!
! Revision 2.246  2020/06/30 23:20:18  pwagner
! ConvertL2GPToQuantity copies pressures, hopefulyy w/o crashing
!
! Revision 2.245  2020/03/20 23:04:13  pwagner
! Made consistent with new he5_readglobalattr api
!
! Revision 2.244  2020/03/04 21:24:48  pwagner
! Make stuff public needed by NCL2GPData; preFill status with l2gp%MissingStatus
!
! Revision 2.243  2020/02/13 21:27:08  pwagner
! Fix errors relating to separate MissingValue for GPH
!
! Revision 2.242  2019/10/30 20:09:37  pwagner
! Prevents an array bounds error caught by NAG
!
! Revision 2.241  2019/10/21 23:20:14  pwagner
! Converts l2gp to quantity even if geolocations not associated
!
! Revision 2.240  2019/05/13 23:32:55  pwagner
! Use UndefinedIntegerValue to prevent failure taking an int of -1.e15 in case of GPH
!
! Revision 2.239  2019/01/29 21:45:47  pwagner
! Initializes some hdfeos character fields to prevent bleed-thru
!
! Revision 2.238  2018/11/12 23:11:12  pwagner
! Deprecated AscDescMode
!
! Revision 2.237  2018/08/01 22:14:27  pwagner
! Added CompactL2GPRecord; corrected errors in ExtractL2GPRecord and Diff-ing with matchTimes
!
! Revision 2.236  2018/05/31 22:47:45  pwagner
! Read ProductionLocation, HostName, identifier_product_doi when dumping global attrs
!
! Revision 2.235  2018/05/22 23:41:14  pwagner
! Use dumpGlobalAttributes when Dumping L2GP Attributes
!
! Revision 2.234  2018/04/24 18:25:11  pwagner
! Commented out the final method; it caused gold brick crashes
!
! Revision 2.233  2018/04/19 00:50:43  vsnyder
! Remove USE statements for unused names, add final subroutine
!
! Revision 2.232  2018/02/28 19:54:26  pwagner
! Tries to prevent unwanted names on each line
!
! Revision 2.231  2018/02/03 00:24:49  pwagner
! Correct Status Bit names; add post-processed
!
! Revision 2.230  2017/12/18 16:41:09  mmadatya
! Added unit for new quantity geoHeight
!
! Revision 2.229  2017/11/30 21:15:53  pwagner
! Avoid writing MissingValue to metadata
!
! Revision 2.228  2017/11/03 19:59:08  pwagner
! Most array gymnastics moved from MLSFillValues to HyperSlabs module
!
! Revision 2.227  2017/07/19 22:51:50  pwagner
! Improved appearance and consistency when Diffing L2GPData types; esp. if gold brick
!
! Revision 2.226  2017/03/10 00:40:19  vsnyder
! Make intrsctn allocatable
!
! Revision 2.225  2016/11/03 22:39:14  pwagner
! Begin transition to the Time_m implmentation of sayTime
!
! Revision 2.224  2016/09/07 22:46:21  pwagner
! Removed unused QuantityType component from L2GPData type
!
! Revision 2.223  2016/07/28 01:42:27  vsnyder
! Refactoring dump and diff
!
! Revision 2.222  2016/06/13 17:59:16  pwagner
! Added cpHE5GlobalAttrs
!
! Revision 2.221  2016/03/23 00:20:09  pwagner
! DiffL2GPData now able to print name on each line
!
! Revision 2.220  2015/10/14 23:18:04  pwagner
! RepairL2GP may optionally repair any geolocations whose chunk number is -999
!
! Revision 2.219  2015/10/13 23:48:28  pwagner
! Warn if all chunkNumbers are FillValues; exit with error if asked to create a swath with preexisting name
!
! Revision 2.218  2015/10/06 00:20:53  pwagner
! Removed unused code
!
! Revision 2.217  2015/09/17 22:57:09  pwagner
! Guards against hdfeos bug in AppendL2GPData only if max chunk size not present
!
! Revision 2.216  2015/09/03 20:26:59  pwagner
! verbose shows timings in AppendL2GPData
!
! Revision 2.215  2015/08/29 00:42:00  vsnyder
! Don't copy undefined data, even if it won't be used
!
! Revision 2.214  2015/07/14 23:13:21  pwagner
! Added computation of interprofile geolocation spacings
!
! Revision 2.213  2015/06/04 17:05:35  pwagner
! Fixed bug causing integer overflow in integer type conversion
!
! Revision 2.212  2015/04/28 16:21:11  pwagner
! Fixed some things in Diffing
!
! Revision 2.211  2015/03/28 01:09:22  vsnyder
! Made DescendingRange a parameter.
! Added stuff to trace allocate/deallocate addresses.
!
! Revision 2.210  2015/02/27 23:57:55  pwagner
! Take pains to ensure HostName is not blank
!
! Revision 2.209  2015/02/18 00:28:01  pwagner
! Corrected DumpL2GP_attributes_hdf5 to include linefeeds, MiscNotes
!
! Revision 2.208  2014/12/11 21:25:57  pwagner
! Make Asc/Desc Mode -999 for crashed chunks
!
! Revision 2.207  2014/12/10 21:27:49  pwagner
! Commented-out unused stuff
!
! Revision 2.206  2014/12/09 01:26:52  pwagner
! Define and initialize orbit angle range for Descending mode
!
! Revision 2.205  2014/10/01 00:02:06  pwagner
! Fixed another bug in ExtractL2GPRecord
!
! Revision 2.204  2014/09/04 23:47:44  vsnyder
! More complete and accurate allocate/deallocate size tracking.
! Add some tracing.
!
! Revision 2.203  2014/07/23 21:59:53  pwagner
! Fixed bugs in ExtractL2GPRecord
!
! Revision 2.202  2014/07/21 21:59:27  pwagner
! Respect options when dumping AscDescMode
!
! Revision 2.201  2014/04/07 17:24:33  pwagner
! Added new l2gp data field AscDescMode
!
! Revision 2.200  2014/04/02 23:03:11  pwagner
! Removed redundant open_ and close_MLSFile
!
! Revision 2.199  2014/01/09 00:24:29  pwagner
! Some procedures formerly in output_m now got from highOutput
!
! Revision 2.198  2013/10/15 23:53:40  pwagner
! May copy quantity values to a file global attribute
!
! Revision 2.197  2013/09/25 00:45:44  pwagner
! Convert must set quantity geolocations, too
!
! Revision 2.196  2013/09/24 00:53:53  pwagner
! Can now convert an l2gpData type to a vector quantity
!
! Revision 2.195  2013/09/17 22:35:40  pwagner
! Changed api of Embed, Extract arrays to match hyperslab
!
! Revision 2.194  2013/08/31 01:24:53  vsnyder
! Replace MLSMessageCalls with trace_begin and trace_end
!
! Revision 2.193  2013/08/12 23:47:25  pwagner
! FindSomethings moved to MLSFinds module
!
! Revision 2.192  2013/04/05 23:17:52  pwagner
! Uses function from dates_module
!
! Revision 2.191  2013/02/26 00:11:49  pwagner
! Dumps times both as tai and hoursInDay
!
! Revision 2.190  2013/02/21 22:22:53  pwagner
! New optional args ti ContractL2GPRecord
!
! Revision 2.189  2013/01/02 21:00:37  pwagner
! Warn instead dying if no swaths to copy
!
! Revision 2.188  2012/10/09 00:34:52  pwagner
! Fixed bugs in OutputL2GP_attributes_MF, AppendL2GPData
!
! Revision 2.187  2012/09/18 18:50:25  pwagner
! Dont write over already-there global attributes
!
! Revision 2.186  2012/07/10 15:18:05  pwagner
! Adapted to new api for GetUniqueList
!
! Revision 2.185  2012/07/05 23:49:04  pwagner
! Sets swath Status field to "crashed" if cpl2gpdata options contains "f" and geolocations contain Fills
!
! Revision 2.184  2012/05/10 00:45:24  pwagner
! Better way for AppendL2GPData to know it must create swath
!
! Revision 2.183  2012/05/01 23:12:34  pwagner
! Maneuver around an unnecessary and sometimes unfortunate if
!
! Revision 2.182  2012/01/11 17:35:20  pwagner
! Turn off extra debugging; snip commented-out lines
!
! Revision 2.181  2011/12/07 01:18:14  pwagner
! Added new 'd' and 'g' options to diff
!
! Revision 2.180  2011/07/07 00:31:03  pwagner
! Treats diffs of geolocation fields with periods
!
! Revision 2.179  2011/02/05 01:33:44  pwagner
! Automatically sheds rank when nFreqs=1; passes options to dump routines
!
! Revision 2.178  2010/11/17 01:19:13  pwagner
! Fixed bug when diffing files with different swathnames; units_name for iwc now 'g/m^3'
!
! Revision 2.177  2010/11/10 02:03:06  pwagner
! Copies correct TAI93At0zOfGranule global attribute
!
! Revision 2.176  2010/06/22 16:53:15  pwagner
! Consistent with new behavior of switchDetail: must include 'f' if override default options
!
! Revision 2.175  2010/06/19 00:11:05  pwagner
! Diffs were checking wrong default field names; fixed
!
! Revision 2.174  2010/02/12 00:25:46  pwagner
! Fixed bug when number of dims is 4 (but what if 5 or more?)
!
! Revision 2.173  2010/02/04 23:08:00  vsnyder
! Remove USE or declaration for unused names
!
! Revision 2.172  2009/11/04 23:15:53  pwagner
! Restored MAXCHUNKTIMES to 120; at 1 DGG file became enormous (who could have guessed)
!
! Revision 2.171  2009/10/22 00:52:39  pwagner
! Second try at fixing hdfeos5 bug; still no luck
!
! Revision 2.170  2009/10/05 23:39:24  pwagner
! Moved use hdf5 statements from module scope to speedup Lahey; this is the last time we do that
!
! Revision 2.169  2009/09/29 23:35:43  pwagner
! Changes needed by 64-bit build
!
! Revision 2.168  2009/08/26 16:45:55  pwagner
! Public global flag WRITEMASTERSFILEATTRIBUTES determines whether to overwrite slave file attributes with masters
!
! Revision 2.167  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.166  2009/06/16 17:21:08  pwagner
! Changed api for dump, diff routines; now rely on options for most optional behavior
!
! Revision 2.165  2009/06/02 17:49:17  cvuu
! Write NRT Lat and Lon to L2Metadata
!
! Revision 2.164  2009/05/14 22:01:21  pwagner
! dumps now take width optional arg
!
! Revision 2.163  2009/05/08 00:42:58  pwagner
! Shows StatusBitNames when dumping Status bits
!
! Revision 2.162  2008/12/02 23:11:41  pwagner
! mls_io_gen_[openF,closeF] functions now private; use MLSFile_T interfaces instead
!
! Revision 2.161  2008/09/10 00:44:58  pwagner
! diffing files can now subset according to geolocation boxes
!
! Revision 2.160  2008/09/09 00:24:49  pwagner
! Fix bug in dumpRange
!
! Revision 2.159  2008/09/03 20:43:09  pwagner
! Added ContractL2GPRecord, diffRange, dumpRange
!
! Revision 2.158  2008/07/09 16:36:31  pwagner
! Fixed "sleepy chunk" syndrome
!
! Revision 2.157  2008/04/22 18:57:33  pwagner
! Fixed bug in ExtractL2GPData that left quality undefined
!
! Revision 2.156  2008/04/18 21:06:59  pwagner
! Fixed error introduced last time
!
! Revision 2.155  2008/02/28 01:27:06  pwagner
! Worked around an apparent bug in HDFEOS5 when adding to swath fields
!
! Revision 2.154  2008/01/24 23:33:04  pwagner
! when diffing, matchtimes will try to match either times or geod. angle
!
! Revision 2.153  2008/01/07 21:37:23  pwagner
! Replace DEFAULTUNDEFINEDVALUE with user-settable undefinedValue
!
! Revision 2.152  2007/12/19 01:29:37  pwagner
! Removed unused args
!
! Revision 2.151  2007/10/10 00:00:51  pwagner
! DumpL2GPData ought not to override optional parameter details if supplied
!
! Revision 2.150  2007/08/13 17:36:38  pwagner
! Push some procedures onto new MLSCallStack
!
! Revision 2.149  2007/06/26 00:20:14  pwagner
! May specify pressure levels on which to show diffs
!
! Revision 2.148  2007/05/30 22:03:41  pwagner
! HuntRange does not work correctly when list has missing values; workaround employed
!
! Revision 2.147  2007/02/26 23:58:48  pwagner
! New optional arg diffs only at matching profile times
!
! Revision 2.146  2007/02/13 21:57:08  pwagner
! Added convergence to defaultfields
!
! Revision 2.145  2006/10/11 00:16:33  pwagner
! Added new data field to hold convergence ratio
!
! Revision 2.144  2006/09/29 23:56:49  pwagner
! Got rid of unused GetGeolocUnits function
!
! Revision 2.143  2006/04/12 20:49:10  pwagner
! nTimesTotal component now dumped, too
!
! Revision 2.142  2006/04/06 23:04:52  pwagner
! Optionally cp only ranges of freq, level, profile
!
! Revision 2.141  2006/03/13 23:40:28  pwagner
! verbose now an optional arg to diff
!
! Revision 2.140  2006/02/28 21:44:44  pwagner
! Diff works better, added IsL2GPSetUp
!
! Revision 2.139  2006/02/21 19:09:25  pwagner
! GetHashElement is now a generic
!
! Revision 2.138  2006/02/16 00:09:23  pwagner
! Show how l2gp, Hgrid shapes differ
!
! Revision 2.137  2006/02/03 21:25:00  pwagner
! Finally ended unnecessary debug printing
!
! Revision 2.136  2006/01/27 01:02:19  pwagner
! Does better at picking units name for column abundances
!
! Revision 2.135  2006/01/26 00:33:09  pwagner
! demoted more use statements from module level to speed Lahey compiles
!
! Revision 2.134  2006/01/19 00:29:10  pwagner
! Units attribute for column abundances corrected
!
! Revision 2.133  2006/01/17 17:49:26  pwagner
! Fixed bug in call to HuntRange
!
! Revision 2.132  2006/01/14 00:54:08  pwagner
! May diff, dump specified chunks only
!
! Revision 2.131  2006/01/04 20:31:19  pwagner
! Diff procedures may keep silent, returning num of diffs only
!
! Revision 2.130  2005/12/16 00:05:57  pwagner
! Changes to reflect new MLSFillValues module; diff from dump0
!
! Revision 2.129  2005/11/04 18:51:19  pwagner
! Non-substantial stylistic tweaks of GEO, DATA_FIELDS
!
! Revision 2.128  2005/10/28 23:12:02  pwagner
! Removed unnecessary HGridField from repairL2GPData
!
! Revision 2.127  2005/10/18 23:06:49  pwagner
! Belated yet crude hack to insure Tropopause units_name is hPa
!
! Revision 2.126  2005/10/11 17:36:55  pwagner
! Added MLSFile interface to cpL2GPData, diff procedures
!
! Revision 2.125  2005/09/23 23:38:30  pwagner
! rename only effective if non-blank
!
! Revision 2.124  2005/09/21 23:16:40  pwagner
! Improvements to RepairL2GP_HGrid; they may even be correct
!
! Revision 2.123  2005/09/14 00:08:15  pwagner
! Dispense with some debug prints; observe HGRid%noProfsLowerOverlap
!
! Revision 2.122  2005/08/19 23:37:22  pwagner
! May use HGrid to repair l2gp geolocations
!
! Revision 2.121  2005/08/15 10:40:38  hcp
! in AppendL2GPData_fileID, mls_swdetach was called with an undefined
! variable as hdfVersion, causing failures. Now fixed.
!
! Revision 2.120  2005/08/08 23:55:04  pwagner
! Never leave status undefined if L2GPFile%errorCode assigned to it
!
! Revision 2.119  2005/08/05 20:32:56  pwagner
! Added RepairL2GP
!
! Revision 2.118  2005/07/12 17:15:47  pwagner
! Dropped global attribute InputVersion
!
! Revision 2.117  2005/07/06 00:31:06  pwagner
! More robust determination of column swaths; quality data field has no units
!
! Revision 2.116  2005/06/29 00:41:50  pwagner
! Passes MLSFiles to mls_swattach,create
!
! Revision 2.115  2005/06/22 17:25:49  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.114  2005/06/14 22:33:45  pwagner
! Moved rdonly error check from read procedure to write
!
! Revision 2.113  2005/06/14 20:38:27  pwagner
! Interfaces changed to accept MLSFile_T args
!
! Revision 2.112  2004/12/14 21:38:59  pwagner
! New verticalCoordinate component, attribute
!
! Revision 2.111  2004/10/27 00:32:35  pwagner
! Missing/Fill value for l2gp%status changed to 513
!
! Revision 2.110  2004/08/04 23:19:01  pwagner
! Much moved from MLSStrings to MLSStringLists
!
! Revision 2.109  2004/08/03 17:59:35  pwagner
! Gets DEFAULTUNDEFINEDVALUE from MLSCommon
!
! Revision 2.108  2004/07/22 17:07:15  pwagner
! Fixed set fill values
!
! Revision 2.107  2004/06/11 19:07:02  pwagner
! Fixed small bug in diffl2gpfiles
!
! Revision 2.106  2004/06/09 00:04:51  pwagner
! New optional args to diff
!
! Revision 2.105  2004/06/02 19:37:52  pwagner
! Added stats option to diff
!
! Revision 2.104  2004/05/26 23:57:29  pwagner
! Added new diff routines
!
! Revision 2.103  2004/05/26 20:31:29  pwagner
! Raised MAXNUMSWATHPERFILE to 250
!
! Revision 2.102  2004/05/05 21:29:53  pwagner
! May change copied swath names with optional arg swathlist2
!
! Revision 2.101  2004/04/23 01:43:34  livesey
! Added Paul's fix for reading status
!
! Revision 2.100  2004/04/06 01:17:43  livesey
! Changed some allocatables to pointers.
!
! Revision 2.99  2004/03/29 18:51:18  pwagner
! Remembered to tell hdf5 status was int, not char
!
! Revision 2.98  2004/03/24 23:53:02  pwagner
! Switched from HE5T_NATIVE_SCHAR to MLS_CHARTYPE
!
! Revision 2.97  2004/03/12 00:38:19  pwagner
! cpL2GPData can convert hdfVersions; hdf version default increased to 5
!
! Revision 2.96  2004/02/26 22:02:36  pwagner
! Acts more gracefully if l2gp file lacks global attributes
!
! Revision 2.95  2004/02/13 00:18:58  pwagner
! Added DumpL2GP_attributes_hdf5 which doesnt work yet
!
! Revision 2.94  2004/02/11 23:05:05  pwagner
! Undid effect of 2nd-to-last fix; seems to have been wrong
!
! Revision 2.93  2004/02/11 23:04:38  pwagner
! Undid effect of 2nd-to-last fix; seems to have been wrong
!
! Revision 2.92  2004/02/11 22:59:59  pwagner
! Undid effect of 2nd-to-last fix; seems to have been wrong
!
! Revision 2.91  2004/02/11 17:23:25  pwagner
! l2gp status an integer, not a char
!
! Revision 2.90  2004/02/10 18:45:43  pwagner
! lastProfile among optional args to AppendL2GPData
!
! Revision 2.89  2004/02/05 23:33:21  pwagner
! Some bug fixes in cpL2GPData and its relatives
!
! Revision 2.88  2004/01/23 01:14:04  pwagner
! Can cp l2gps from one file to another
!
! Revision 2.87  2004/01/09 00:20:56  pwagner
! Added avoidUnlimitedDims to allow bypassing bug directWriting range of chunks
!
! Revision 2.86  2003/12/03 17:51:14  pwagner
! L2GP tracks both nTimes (for this slave) and nTimesTotal (done by all)
!
! Revision 2.85  2003/11/19 22:14:08  livesey
! Added option (not invoked at the moment) for 'limited' dimensions in
! AppendL2GP
!
! Revision 2.84  2003/11/05 21:44:38  pwagner
! Skips trying to read 4 mls-specific fields from non-mls files
!
! Revision 2.83  2003/10/31 12:35:57  hcp
! HCP made some small fixes to AppendL2GPData_fileID so that it doesn't
! crash if optional arg TotNumProfs is not present
!
! Revision 2.82  2003/10/28 21:41:06  pwagner
! Removed swath-level MissingValue attribute; renamed -Precision softlink
!
! Revision 2.81  2003/10/28 00:39:00  pwagner
! Fixed bug where character-vlaued attributes were only 1 char long
!
! Revision 2.80  2003/09/09 23:03:46  livesey
! More rigorous 'empty l2gp' check in append.
!
! Revision 2.79  2003/09/08 22:51:54  pwagner
! Simplified and unified read/write of hdfVersions
!
! Revision 2.78  2003/09/04 03:02:53  livesey
! Removed some print statements
!
! Revision 2.77  2003/09/03 22:39:44  pwagner
! Uses actual num of profiles for nTimes even if Unlim
!
! Revision 2.76  2003/08/28 23:51:03  livesey
! Various bug fixes to the AppendL2GP stuff
!
! Revision 2.75  2003/08/27 20:06:26  livesey
! Removed some print statements
!
! Revision 2.74  2003/07/21 23:32:09  pwagner
! Will write attributes, create alias when appending hdfeos5 if mustcreate
!
! Revision 2.73  2003/07/15 23:35:11  pwagner
! Disabled most printing; uses mls_SWdetach
!
! Revision 2.72  2003/07/10 22:18:49  livesey
! More minor bug fixes, but lots of print statements.
!
! Revision 2.71  2003/07/09 21:49:06  pwagner
! Wont try swattaching just to see if swath already there
!
! Revision 2.70  2003/07/08 00:43:17  livesey
! Bug fix in zero length chunk test
!
! Revision 2.69  2003/07/07 21:04:55  pwagner
! Tries to deal sensibly with profile-less chunks
!
! Revision 2.68  2003/07/02 00:55:27  pwagner
! Some improvements in DirectWrites of l2aux, l2gp
!
! Revision 2.67  2003/06/26 00:04:46  pwagner
! Added optional DONTFAIL arg to MLS_SWATTACH
!
! Revision 2.66  2003/06/20 19:31:39  pwagner
! Changes to allow direct writing of products
!
! Revision 2.65  2003/06/09 22:44:09  pwagner
! Uses mls_swcreate
!
! Revision 2.64  2003/04/21 19:34:59  pwagner
! Restored reading/writing char-valued Status datafield
!
! Revision 2.63  2003/04/17 23:07:02  pwagner
! Now uses MLSHDFEOS more thoroughly; can read HIRDLS L2GP files
!
! Revision 2.62  2003/04/15 23:14:50  pwagner
! New chunking so hdfeos5 files dont balloon; some swsetfill tweaking
!
! Revision 2.61  2003/04/11 23:35:10  pwagner
! Added new UniqueFieldDefinition attribute; sets fill and MissingValue 
! attributes for all fields
!
! Revision 2.60  2003/04/03 22:58:40  pwagner
! Alias now set in lib/L2GPData instead of l2/write_meta
!
! Revision 2.59  2003/04/02 23:53:56  pwagner
! Checks for FILENOTFOUND
!
! Revision 2.58  2003/03/07 00:40:23  pwagner
! Call HE5_SWsetfill; removed some spaces from attribute names
!
! Revision 2.57  2003/02/26 17:36:11  pwagner
! Repaired ReadL2GPData to close swath when given WILDCARDHDFVERSION
!
! Revision 2.56  2003/02/21 23:41:53  pwagner
! Also writes Fill Value attribute
!
! Revision 2.55  2003/02/12 21:50:33  pwagner
! New api for ReadL2GPData lets you supply fileName and hdfversion, which can be wildcard
!
! Revision 2.54  2003/02/10 22:04:49  pwagner
! ChunkNumber correctly HE5_SWdefgfld-ed as an int
!
! Revision 2.53  2003/02/08 00:33:32  pwagner
! Writes he5 attributes w/o bombing
!
! Revision 2.52  2003/02/07 01:05:12  pwagner
! rgp now public
!
! Revision 2.512003/02/06 00:24:48  pwagner
! Squashed a (last?) bug in hdfeos2 stuff
!
! Revision 2.50  2003/02/04 22:43:29  pwagner
! Fixed another serious bug
!
! Revision 2.49  2003/02/03 21:33:15  pwagner
! Having fixed (most)bugs, brought back from he5lib
!
! Revision 1.28  2003/01/30 00:58:53  pwagner
! Writing first attributes for hdfeos5
!
! Revision 1.27  2003/01/15 23:23:34  pwagner
! data types for L2GPData_T now adjustable
!
! Revision 1.26  2003/01/15 19:13:30  pwagner
! Smane no monkeying fix, but for hdf5
!
! Revision 1.25  2003/01/15 00:05:34  pwagner
! Now can write L2GPData w/o monkeying around
!
! Revision 1.24  2003/01/09 01:02:27  pwagner
! Moved some use statements in attempt to work around Lahey long compile time bug
!
! Revision 1.23  2002/10/25 15:21:05  hcp
! Nasty hack for unlimited swaths removed. Local requirement for this worked
! around in a different way. HDF-EOS5 team admit that problem is really
! caused by a bug in HDF-EOS5 so it should get fixed sometime anyway.
!
! Revision 1.22  2002/10/23 16:59:14  hcp
! Added _VERY_ cheesy work-around for reading files where the time
! dimension is unlimited.  This should not have any effect on files where
! the time dimension is not an unlimited dimension.  This is a
! work-around for some pathetically poor behaviour in HDF-EOS5.
!
! Revision 1.21  2002/10/09 14:05:30  hcp
! replaced HE5S_UNLIMITED with   HE5S_UNLIMITED_F
!
! Revision 1.20  2002/08/15 22:26:22  pwagner
! Added he5_swdefchunk to .._createFile_hdf5
!
! Revision 1.19  2002/06/10 12:00:43  hcp
! Removed use hdf5_params as hdf-eos5 now has an include file with the
! parameters you need ready-defined. References to various types changed
! to be consistent with that.
!
! Revision 1.18  2002/05/01 13:05:53  hcp
! Changed a warning to a debug so I didn't have to see it
!
! Revision 1.17  2002/05/01 09:28:10  hcp
! Some print statements commented
!
! Revision 1.16  2002/03/15 23:02:29  pwagner
! Gets HDFVERSION_4 and 5 from MLSFiles; checks for illegal hdfversions
!
! Revision 1.15  2002/02/01 21:32:34  pwagner
! offset treated properly for appendl2gp for hdf4; untested
!
! Revision 1.14  2002/01/29 23:47:21  pwagner
! Repaired bugs relating to hdf4 compatibility
!
! Revision 1.13  2002/01/29 00:48:43  pwagner
! Now should handle both hdfVersions; not tested yet
!
! Revision 1.12  2002/01/26 00:18:05  pwagner
! (Read)(Write)L2GPData accepts optional hdfVersion arg
!
! Revision 1.11  2001/11/29 10:07:13  pumphrey
! L2GPData (HDF-EOS5 version) brought up to date with HDF-EOS4 version
!
! Revision 1.10  2001/11/28 17:46:28  pumphrey
! In the middle of syncing up l2gpdata with HDF4 version in lib.
! Compiles, but not tested. Hack to write char swaths cut-n-pasted but not
! examined for sanity
!
! Revision 1.9  2001/11/28 16:17:16  pumphrey
! Syncing prior to major re-sync with HDF4 version.
!
! Revision 1.8  2001/07/31 11:26:19  archie
! Corrected case for ChunkNumber
! .
!
! Revision 1.7  2001/07/11 19:01:16  pumphrey
! quality->Quality, status->Status.
!
! Revision 1.6  2001/04/27 07:48:54  pumphrey
! Many nested loops in l3ascii replaced with array ops. Small fixes
! (e.g. spelling mistakes) in other modules.
!
! Revision 1.5  2001/04/06 20:16:38  pumphrey
! Not much, just keeping in sync
!
! Revision 1.4  2001/03/29 17:33:27  pumphrey
! Huge changes to L2GPData to sync with the HDF4 version and add unlimited
! dimension along the track
!
! Revision 1.3  2001/03/20 14:00:30  pumphrey
! fixing inconsistencies -- nothing important
!
! Revision 1.2  2001/02/23 13:33:14  pumphrey
! Fixed type definition for L2GPData_T so all the pointers are => NULL()
! This error was detected by nagf95 on Solaris but not on Linux. Odd.
!
! Revision 1.1  2000/12/22 15:55:53  pumphrey
! Initial commit of HDF-EOS5 versions of L2GP interface.
!
! Revision 2.8  2000/12/04 23:43:59  vsnyder
! Move more of addItemToDatabase into the include
!
! Revision 2.7  2000/09/22 14:29:42  pumphrey
! OutputL2GP_createFile was setting LocalSolarTime as byte (should be
! REAL)  and Time as REAL (should be Double prec.). Fixed.
!
! Revision 2.6  2000/09/21 13:48:18  pumphrey
! fixed a bug in the write routine.
!
! Revision 2.5  2000/09/19 12:42:11  pumphrey
! added chunkNumber to SetupNewL2GPRecord and other bug fixes
!
! Revision 2.4  2000/09/18 10:19:49  pumphrey
! Removed some debugging statements.
!
! Revision 2.1  2000/09/15 21:50:18  livesey
! New version of L2GP data, moved some stuff from l2 to lib
!
! Revision 2.0  2000/09/05 18:57:03  ahanzel
! Changing file revision to 2.0.
!
! Revision 1.1  2000/09/02 02:05:04  vsnyder
! Initial entry
!
